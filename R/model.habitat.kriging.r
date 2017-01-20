

  model.habitat.kriging = function( model.type="gam.full", p=NULL, hvar="R0.mass", plotdata=T,
    dat= snowcrab.db( DS="set.complete"),
    outdir = file.path( project.datadirectory("bio.snowcrab"), "R", "habitat", "models" )
    ) {

    dir.create(path=outdir, recursive=T, showWarnings=F)

    qq =  grep( "redo", model.type)
    if ( length(qq) != 0  ) {
      set = dat
      q = quantile( set[which(set[,hvar]>0), hvar], p$habitat.threshold.quantile )
      set$Y = 0
      set$Y[ which( set[,hvar] > q ) ] = 1
    }

    if ( model.type %in% c("gam.full", "gam.full.redo") ) {

			fn = file.path( outdir, paste( "habitat.gam.full", hvar, "rdata", set="." ) )

      if ( model.type == "gam.full" ) {
        load(fn)
        return( Q)
      }
      require(mgcv)
      require(arm)

#      .M.gam = formula( Y ~  s( yr ) + s( tmean ) + s( I(t-tmean) ) + s( tamp) + s( plon, plat) + s( z)
#        + s( log.substrate.grainsize) + s( ddZ) +s( dZ)  )
      .M.gam = formula( Y ~  s( yr ) + s( t )  + s( tamp) + s( plon, plat) + s( z)
        + s( log.substrate.grainsize) + s( ddZ) +s( dZ)  )
      Q = gam( .M.gam, data=set, na.action="na.pass", family=binomial() )
      AIC (Q) # = 3367.65
      summary(Q)
      # P = predict( Q, set, type="response")
      # cor(P,set$Y)^2 # = 0.41
      save( Q, file=fn, compress=T )
      if (plotdata) plot(Q, all.terms=T, rug=T, jit=T, seWithMean=T, trans=invlogit, pers=T , scale=0)
      return ( "Complete" )
    }



    if ( model.type %in% c("glm.splines.redo", "glm.splines") ) {
      fn = file.path( outdir, "habitat.gam.full.rdata" )
      if ( model.type == "glm.splines" ) {
        load(fn)
        return( Q)
      }

      # interaction terms were significant but were heaviliy influenced by extremes in data .. they were dropped
      .M.splines = formula( Y ~  bs(yr, df=4)+  tmean + I(t-tmean) + bs(tamp, df=3) + bs( plon, df=3)* bs(plat,df=4)+  bs(z, df=4) + bs(log.substrate.grainsize, df=3) + bs(dZ, df=3)+  bs(ddZ, df=3)
        + bs(dZ,df=3)
      )
      Q = glm( .M.splines, data=set, na.action="na.omit", family=binomial() )
      AIC (Q) # 3573.2
      summary(Q)
      # p = predict( Q, set, type="response")
      cor(p,set$Y)^2 #
      if (plotdata) plot(Q, all.terms=T, rug=T, jit=T, seWithMean=T, trans=invlogit, pers=T )
      save( Q, file=fn, compress=T )
      return ( "Complete" )
    }

    if (model.type == "gamm.redo" ) {
          # GAMM
          ii = sample(1:nrow(set), 200 )
          set$id = factor( paste( set$plon, set$plat, sep="~") )
          set = set[, c("id", "plon", "plat", "Y", "yr.factor", "t", "z", "log.substrate.grainsize", "ddZ", "dZ") ]
          .M.gamm = formula( Y ~ yr.factor + s(t) + s(z, log.substrate.grainsize) + s(ddZ, dZ) - 1 )
          # Q = gamm( .M.gam, random=list(yr.factor=~1), data=set, na.action="na.omit", family=binomial() )
          # deviance = 3014.19       2721  df; ~ 36% of deviance explained
          Q = gamm( .M.gamm, random=list(yr.factor=~1+id), data=set, na.action="na.omit", family=binomial(), subset=ii)
          save( Q, file=file.path( outdir, "q.gamm.rdata" ), compress=T )
    }

    if (model.type =="gamm.correl" ) {
          # GAMM  with correlation
          set$plon = jitter(set$plon, amount=1)
          set$plat = jitter(set$plat, amount=1)
          # ~ 24 hr
          Q = gamm(  Y ~ s(yr) + s(z) + s(log.substrate.grainsize) + s(dZ) + s(ddZ) + s( tamp) + s( log.substrate.grainsize) , correlation=corGaus(form=~plon+plat),
          data=set, na.action="na.omit", family=binomial())
          # r = predict.gam( q.gamm.cor$gam, set, type="response")
          # cor(r,set$Y)^2 =  0.182554
          save( Q, file=file.path( outdir, "q.gamm.cor.rdata" ), compress=T )
    }

    if (model.type =="gamm.correl" ) {
          # ~ 68 hrs
          Q = gamm(  Y ~ yr.factor + s(t) + s(plon, plat) + s(z, log.substrate.grainsize) + s(dZ, ddZ),
            correlation=corGaus(form=~plon+plat),
            data=set, na.action="na.omit", family=binomial())
            r2 = predict.gam( q.gamm.cor2$gam, set, type="response")
            cor(r2,set$Y)^2
          # = 0.377711
            save( Q, file=file.path( outdir, "q.gamm.cor2.rdata" ), compress=T )
    }

  }


