
  # -------------------------
  # forward projection
  
  # ++++++++++ must check timing of 4X to make sure that the season is not an issue ++++++++++
  p=list()
  p$init.files = loadfunctions( "snowcrab", functionname="initialise.local.environment.r") 
 
  p$clusters = "localhost"
  # p$clusters = rep("localhost",2)
  # p$clusters = rep("tethys",8)
  # p$clusters = rep( c("kaos", "nyx", "tartarus") ,24) 

  p = get.population.projection.parameters( p )
  
# start.projection.year = 2004
  p$start.projection.year = 2009
  p$nyears.projection = 6
  p$pyears.projection = c( p$start.projection.year : (p$start.projection.year + p$nyears.projection))

  p$good.years = as.character( c(2003:p$current.assessment.year) )  # use only the autumn surveys ...
# p$scenario = "historic"
  p$scenario = "mature.only"
  p$exploit.rate = c(0, 0.1, 0.2, 0.4, 0.6, 1)



  # mean sizes from trawl surveys ... slow worth parallelizing in future
  update.files=F # redo transition matrices, and size estimates?
  if (update.files) {
    sizes = mean.weights.by.category( p, redo=T ) 
    tmatrix = transition.matrices( p, redo=T) 
  }

  p = make.list( list( er = 1:length(p$exploit.rate))  , Y=p )
  parallel.run( sim.markov, p=p, DS="redo" ) 
  sim.markov( p=p, DS="redo" ) 
 
  fp = sim.markov( p=p, DS="collect" )
  fp = sim.markov( p=p, DS="saved" )



  td = create.projection.summary( p, fp$FB, rescale.results=F )
 
  for (ir in 1:length(p$regions)) {
    x11()
    regs= p$regions[ir]
    reg.label = p$regions.labels[ir]
    pj = td[ which( td$region==regs ) ,]
    plot.projections( pj, reg=regs, outdir=file.path(p$annual.results,"ts"), reg.label )
  }


#tables for output

# the mean size by instar
  sizes = mean.weights.by.category( p )  # use the previously saved file
  
  latex( round( apply(sizes[,"cw.mean",,],c(1,3),mean, na.rm=T),1), file="", title="", label="", rowlabel="Stage", cgroup="Mean carapace width (cm)", na.blank=T, caption="Caption cw")
  latex( round( apply(sizes[,"mass.mean",,],c(1,3),mean, na.rm=T),1), file="", title="", label="", rowlabel="Stage", cgroup="cg", na.blank=T, caption="Caption cwsd")


  tmatrix = make.transition.matrices(p) 
  

# fishing mortality  
  eps = 0.001
  FM = tmatrix$FM
  FM[ which(!is.finite(FM)) ] = NA
  FM[ which(FM==1) ] = NA
  FM[ which(FM<eps) ] = NA

  three.yr.retro = as.character( p$current.assessment.year + c(0, -1, -2))
  FM = FM[, three.yr.retro, ]
  fm = apply( FM, c(1,3), mean, na.rm=T )
  fm = apply( FM, c(1,3), mean, na.rm=T )
  fm.sd = apply( FM, c(1,3), sd, na.rm=T )


fm = round( fm,3); fm[ !is.finite(fm)] = "-"; fm[ fm < eps] = "-"
fm.sd = round( fm.sd,3); fm.sd[ !is.finite(fm.sd)] = "-" ; fm.sd[ fm.sd < eps] = "-"


latex(fm , file="", title="", label="", rowlabel="Stage", collabel="Stage", na.blank=T, caption="Exploitation rates; three year mean")
latex(fm.sd , file="", title="", label="", rowlabel="Stage", collabel="Stage", na.blank=T, caption="Exploitation rates; sd")


# transition matrices

  TM = tmatrix$TM[,, p$good.years ,] # good TM estimates
  tm.n = choose.transition.matrix( TM, region="cfanorth", threshold.to.delete=3 )$XM # threshold.to.delete is the upper bound of a TM to believe ..  this is what is used for projection 
  tm.s = choose.transition.matrix( TM, region="cfasouth", threshold.to.delete=3 )$XM # threshold.to.delete is the upper bound of a TM to believe ..  this is what is used for projection 
  tm.x = choose.transition.matrix( TM, region="cfa4x", threshold.to.delete=3 )$XM # threshold.to.delete is the upper bound of a TM to believe ..  this is what is used for projection 

eps =0.001
tm.n = round( tm.n,2); tm.n[ tm.n < eps] = "-"
tm.s = round( tm.s,2); tm.s[ tm.s < eps] = "-"
tm.x = round( tm.x,2); tm.x[ tm.x < eps] = "-"


latex(tm.n , file="", title="", label="", rowlabel="Stage", collabel="Stage", na.blank=T, caption="Pseudo-Markov transition matrix used for projections in N-ENS")
latex(tm.s , file="", title="", label="", rowlabel="Stage", collabel="Stage", na.blank=T, caption="Pseudo-Markov transition matrix used for projections in S-ENS")
latex(tm.x , file="", title="", label="", rowlabel="Stage", collabel="Stage", na.blank=T, caption="Pseudo-Markov transition matrix used for projections in 4X")
 

# scenarios 
  i = c(1:10)
  n = array(data=NA, dim=length(i))
  n[1] = 100
  ex = 0.2
  nm = 0.25
  for(j in i[-1]) {
    n[j] = n[j-1] * (1-(ex+nm))
  }
  plot(n)


