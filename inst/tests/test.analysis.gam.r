


Test:     Biomass estimation via GAM / LME 



  require(mgcv)
  require(nlme4)
  
	loadfunctions( "snowcrab", functionname="initialise.local.environment.r") 
 

      
    p$model.type = "gam.full"
      # p$model.type = "gamm.full"

    # model habitat using GAM
    
    # quantile at which to consider zero-valued abundance
    p$habitat.threshold.quantile = 0.05 
      
    set = snowcrab.db( DS="set.logbook" )
    set$Y = set$R0.mass

    set$plon= jitter(set$plon)
    set$plat= jitter(set$plat)
 
    Q = quantile( set$Y[ which(set$Y>0)],  p$habitat.threshold.quantile )
    
    g =set
    g = g[ which(set$Y> p$habitat.threshold.quantile ) ,]
    g = g[ , c("Y", "yr", "cfa", "plon", "plat", "t", "tamp", "wmin", "z", "substrate.mean", "dZ", "ddZ"  ) ]
    
    Vrange = NULL
    Vpsill = NULL
    Zannual = NULL
 
    for ( iy in sort(unique(g$yr)) ) {
      ii = which(g$yr==iy)

      Z = var( g$Y[ii], na.rm=T )
      V = variogram( Y ~ 1, locations=~plon + plat, data=g[ii ,] , cutoff=100, boundaries=c(10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80, 100 ) , cressie=T)
      V$gamma = V$gamma / Z  # normalise it
      vf = fit.variogram(V, vgm(psill=0.5, model="Sph", range=50, nugget=0.5 ))
      plot( V, model=vf )

      Vrange = c( Vrange, vf$range[2] )
      Vpsill = c( Vpsill, vf$psill[2] / (vf$psill[1] + vf$psill[2] ) )
      Zannual = c( Zannual, Z )

    }
    
    Vrange = mean(Vrange[which(Vrange< 100)] , na.rm=T) 
    Vpsill = mean(Vpsill[which(Vrange< 100)], na.rm=T)

     
    # g = g[ which( g$yr %in% c(2008:2009)) ,]
        #  g = g[ sample(1:nrow(g), 1000 ) ,]


    # g$yr =  as.factor( g$yr )
    # g$cfa =  as.factor( g$cfa )
   
    # rm(set); gc()

    Q = gamm( 
          formula = Y ~  yr + s(plon,plat) +  s(t)  + s( tamp) + s( wmin) +  s( z) +  s( dZ) +s(ddZ), 
          correlation=corSpher(value=c( Vrange, Vpsill ), form=~plon+plat|yr, nugget=T ), 
          data=g, family=Gamma()
    )

    Ql = gamm( 
          formula = Y ~  s(yr) + s(plon,plat)  + s( tamp) + s( wmin) +  s(t)+ s( z) + s( dZ), 
          correlation=corSpher(value=c( Vrange, Vpsill ), form=~plon+plat|yr, nugget=T ), 
          data=g, family=Gamma("log")
    )


    save( Ql, file=file.path( project.datadirectory("snowcrab"), "test.Ql.gam.rdata", compress=T )


    cor( predict( Q$gam, g, type="response"), g$Y, use="pairwise.complete.obs")^2 # = 0.038

R = gls( Y ~  yr * cfa + ( bs( t )  + bs( tamp) + bs( wmin) ) + bs( z) + bs( substrate.mean) + bs( ddZ) +bs( dZ), 
correlation = corGaus(c( Vrange, Vpsill ), form=~plon+plat | yr, nugget=T ), data=g )
cor( predict( R, g, type="response"), g$Y, use="pairwise.complete.obs")^2 # = 0.378962
o = allEffects(R)



C = snowcrab.db( DS="cat.initial" )
S = sort( unique( C$spec ) )
rm (C)

ms.mass = paste( "ms.mass", S, sep="." )
ms.no = paste( "ms.no", S, sep="." )

other = c( "total.effort" ,"total.cpue" , "total.landings", "total.visits", "cpue", "landings", "visits" ) 
hab = c( "substrate.mean" , "dZ", "ddZ" , "tmean" , "tamp", "wmin" , "thp" , "tsd.H", "t" , "tamp.cl" , "wmin.cl", "tsd" , "t" , "tsd", "z"  )
crab = c( "cw.male.mat.mean","cw.fem.mat.mean","totmass.female.berried","totmass.female.primiparous","totmass.female.multiparous","mi123.no","mi4.no","mi5.no","mi6.no","mi7.no","mi8.no","mi9.no","mi10.no","mi11.no","mi12.no","fi1234.no", "fi5.no", "fi6.no" ,"fi7.no", "fi8.no", "fi9.no", "fi10.no", "R0.no", "R1.no", "R2.no", "R3.no", "R4.no", "R5p.no", "dwarf.no", "ma9.no", "ma10.no" ,"ma11.no", "ma12.no", "ma13.no", "fa7.no", "fa8.no" ,"fa9.no","fa10.no"  )

log10Y = log10( set$Y + min(set$Y[set$Y>0])/2 )


m = lm( log10Y ~ ., data=set[,other ] )
summary(m)
n = stepAIC(m)

i=which( is.finite( rowSums( set[,hab] ) ))
m = lm( log10Y[i] ~ ., data=set[i,hab ] )
summary(m)
n = stepAIC(m)

m = lm( log10Y ~ ., data=set[,crab ] )
summary(m)
n = stepAIC(m)

m = lm( log10Y ~ ., data=set[,ms.no ] )
summary(m)
n = stepAIC(m)


