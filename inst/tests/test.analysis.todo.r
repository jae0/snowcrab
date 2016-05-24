

# 1 .generate Z estimates and possibly relate them to changes in prey and pred abundance.
  
	loadfunctions( "snowcrab", functionname="initialise.local.environment.r") 
  q = get.population.projection.parameters( assessment.year = p$current.assessment.year )
  
  scenario = "mature.only"
  
  sizes = mean.weights.by.category( p )
  tmatrix = make.transition.matrices(p) 
 
#  fp = forward.project.main ( p, q, scenario, exploit.rate, start.projection.year, nyears.projection, good.years, source="file" )  # load the saved version

  eps = 1e-4  # assume this is min resolution


    nnodes = length(nodes)
    fall.surveys = c(2003)
    spring.surveys = c(1998:2001)
  
    
    sex = "male"    
    region=  q$regions[1]


    
    sur = get.annual.timeseries( p, sex=sex, outtype="yearclass", region=region)

    fm = landings.decomposed ( p, q, sex=sex, region )
   
    # correct timing of landings to account for change in survey season (by offsetting landings)
    fm[,,"2002"] = fm[,,"2002"] + fm[,,"2001"]
    for (y in 2001:1999)  fm[,,as.character(y)] = fm[,,as.character(y-1)]
    fm[,,"1998"] = fm[,,"1998"] * NA


    N.pre = sur$TS + fm   # ... need to fix this 
    N.pre[ N.pre < eps ] = 0
    N.post = sur$TS 
  
##################

    cl = make.classes (sex)
    yclass = cl$yclass
    cats = cl$cats
    varmap = cl$varmap
    yrs =  dimnames(N.pre)[[3]]

    k = get.cohort.data (yclass, yrs, cats, varmap, N.pre, outvar="total")
    l = get.cohort.data (yclass, yrs, cats, varmap, N.post, outvar="total")
    

    plot(rowSums(k[,,"1998"]))
    lines(rowSums(l[,,"1998"]))
 
#################

    for (i in q$nodeyears) {
        k2 = apply(k,c(1,3), sum, na.rm=T)
        k2[ which(k2==0) ] = NA
        
        k3 = apply(k2, 1, mean, na.rm=T) 

        plot ( as.numeric(as.character(rownames(k2))), k2[,"1998"])

        Z =  extract.mortality( k, i, q$nodes, fishing, sex="male")

      }
    }



# 2. try a vpa 


###################

# 3. opengl interactive 3d plots:

  require(rgl)

  data(volcano)

     z <- 2 * volcano        # Exaggerate the relief

     x <- 10 * (1:nrow(z))   # 10 meter spacing (S to N)
     y <- 10 * (1:ncol(z))   # 10 meter spacing (E to W)

     zlim <- range(y)
     zlen <- zlim[2] - zlim[1] + 1

     colorlut <- terrain.colors(zlen) # height color lookup table

     col <- colorlut[ z-zlim[1]+1 ] # assign colors to heights for each point

     open3d()
     surface3d(x, y, z, color=col, back="lines")




  # -------------------------------------------------------------------------------------
  # source required functions

	loadfunctions( "snowcrab", functionname="initialise.local.environment.r") 


  # -------------------------------------------------------------------------------------
  # not yet complete ....  +++++++++++++++++++++++++++++++++++++
    if (compute.cummulative.landings) {
       x = get.fisheries.data ( source="file")
       xc = cummulative.landings(x)
    }


# ---------------------------------
# exploratory data plots ... generic timeseries of landings ..
  landings =  landings.db()

#vars = c("tac.tons","landings.tons","cpue.kg.trap","effort.100traps")
#for (v in vars) plot.fisheries.ts(landings, v, "ts.fisheries")

#ts.landings.and.tacs(out$yr, out$landings.kt, out$tac.tons, outdir="ts.fisheries", outfile="landings", title="Landings")

#ts.landings.and.tacs(out$yr, out$cpue, out$tac.tons, outdir="ts.fisheries", outfile="cpue", title="CPUE")

#ts.landings.and.tacs(out$yr, out$effort, out$tac.tons, outdir="ts.fisheries", outfile="effort", title="Fishing effort")


# -----------------------------------
# 6. Biomass estimates from simple interpolations

  #  key flags
  p$interp.method ="tps"
  # p$interp.method ="inv.dist.gstat"

  p$ofname = paste("B", p$interp.method, p$tension, p$maskres, p$interpres, "rdata", sep=".")
  K = interpolated.estimates ( list (
        outfile = file.path( project.datadirectory("snowcrab"), "R", p$ofname),
       make.biomass.grids = T,
       biomass.regrid = T,
       biomass.extract =T,
       years = p$years.to.model,
       vars = p$vars.to.model,
       regions = p$regions.to.model
  ) )

  # analysis
  load (B, file=fname, compress=T)

  i = which( B$var=="totmass.male.com" & B$region=="cfaall" )
  print( B[i ,] )
  plot(B$yr[i] , B$B[i])


# -----------------------------------
# 7. to test influence of block size
  load("temp.fix.pred.surface.rdata")

  P.regridded = block.xyz2xyz( xyz=PS, params=p, engine="R2")

  P.regridded = P.regridded[ is.finite(P.regridded[,"z.m"]) ,]
  bk = P.regridded[, c("plon", "plat")]

gc()
  i = NULL
  z = NULL
  i = filter.region.polygon(x=bk[, c("plon", "plat")], region=r, planar=T)
  o = bk[i, c("plon", "plat")]
  surfacearea = length(i) * (p$pres*p$pres) / 1000 # divide by 1000 to make kt or kn
  z = predict(object=g, newdata=data.frame(plon=0, plat=0), block=o)
  #   to do .. cokriging here ... too slow to do yet
  meanval = z[, k.pname]
  stdev = sqrt(z[, k.vname])
  total =  decode.variable(x=meanval, v, offset) * surfacearea
  lbound = decode.variable(x=(meanval- 1.96*stdev), v, offset) * surfacearea
  ubound = decode.variable(x=(meanval+ 1.96*stdev), v, offset) * surfacearea
  sse = ifelse(is.null(attr(final.vgm, "SSErr")), NA, attr(final.vgm, "SSErr") )
  K0 = data.frame( yr=y, vars=v, region=r, meanval, stdev, total, lbound, ubound,
    surfacearea, vario.model=final.vgm$model[2], vario.sse=sse,
    psill=final.vgm$psill[2], nugget=final.vgm$psill[1], range=final.vgm$range[2]  )
  print( K0 )

  a = read.table( file.path( project.datadirectory("snowcrab"), "R", "block.size.dat"), header=T) # size influence of block dimensions
  errbar(log10(a$res), a$totno, a$ub, a$lb, xlab="Log10(block size; km)", ylab="No. x 1000", lty=1, lwd=2, ylim=c(0,10000) )
  lines(log10(a$res), a$totno, lty=2, col="orange")
  dev.print(pdf, "block.size.pdf")



