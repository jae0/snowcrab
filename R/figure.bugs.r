
  figure.bugs = function( vname="", type="density", sb=NULL, y=NULL, fn=NULL, labs=c("N-ENS","S-ENS","4X") ,save.plot=T) {
 
    ntacs = sb$nProj
    yrs0 = as.numeric( as.character( rownames(sb$IOA) ) )
    yrs = c( yrs0, (max(yrs0)+c(1:sb$M) ) ) 
    yrs.last = max(yrs0) + 0.5
    ndata = length(yrs0)
    hdat = 1:ndata


    if (vname =="r.ts") {
      # catch this first as the layout is different
        br = 75

        x11()
        layout( matrix(c(1:(sb$N*3)), ncol=3, nrow=sb$N ))
        par(mar = c(1., 1., 0.65, 0.75))

        for (i in 1:3) {
        for (yr in 1:sb$N) { 
          dta = y$r[ yr,i,,]
          qs = apply( dta, 2, quantile, probs=c(0.025, 0.5, 0.975) )
          qs = signif( qs, 3 )
          pdat = as.vector( dta)
          xrange = range( pdat, na.rm=T )
          postdat = hist( pdat, breaks=br, plot=FALSE )
          yrange = range( 0, postdat$density, na.rm=T ) * 1.02
          hist( pdat, freq=FALSE, breaks=br, xlim=xrange, ylim=yrange, main="", xlab="", ylab="Density", col="lightgray", border="gray")  
          YR = rownames(sb$IOA) [yr]
          legend( "topright", bty="n", legend=paste( labs[i], YR, "r", " = ", qs[2,i], " {", qs[1,i], ", ",  qs[3,i], "}  ", sep="" ), cex=0.9 )   
        }}
        
     if(save.plot) savePlot( filename=fn, type="png" )
      return( fn)

    }


    x11()
    layout( matrix(c(1,2,3), 3, 1 ))
    par(mar = c(4.4, 4.4, 0.65, 0.75))
   
    prr = NULL
    prr$class ="none"  # no priors by default


    if ( type=="density" ) {  # default
    
      if ( vname=="K" ) {
        qs = apply( y$K[,,], 1, quantile, probs=c(0.025, 0.5, 0.975) )
        qs = signif( qs, 3 )
        for (i in 1:3) {
          pdat = as.vector(y$K[i,,])
           prr=NULL
           prr$class="lognormal"

          # E(X) = exp(mu + 1/2 sigma^2)
          # med(X) = exp(mu)
          # Var(X) = exp(2*mu + sigma^2)*(exp(sigma^2) - 1)
          # CV = sqrt( Var(X) ) / E(X) = sqrt(exp(sigma^2) - 1) ~ sigma; sigma < 1/2
          # SD(X) = sqrt( exp(2*mu + sigma^2)*(exp(sigma^2) - 1) )
          #  or   = CV * E(X) 

           prr$meanlog= sb$K.mu[i]
           prr$sdlog = sqrt(sb$K.sd[i])
          plot.freq.distribution.prior.posterior( prior=prr, posterior=pdat )
          
          legend( "topright", bty="n", legend=paste( labs[i], "\n", vname, " = ", qs[2,i], " {", qs[1,i], ", ",  qs[3,i], "}  ", sep="" ))
      }}

      if ( vname=="r" ) {
        qs = apply( y$r[,,], 1, quantile, probs=c(0.025, 0.5, 0.975) )
        qs = signif( qs, 3 )
        for (i in 1:3) {    prr=NULL
          prr=NULL
          prr$class='normal'
          prr$mean=sb$r.mu[i]
          prr$sd=sqrt(sb$r.sd[i])
          pdat = as.vector(y$r[i,,])
          plot.freq.distribution.prior.posterior( prior=prr, posterior=pdat )
          legend( "topright", bty="n", legend=paste( labs[i], "\n", vname, " = ", qs[2,i], " {", qs[1,i], ", ",  qs[3,i], "}  ", sep="" ) )   
      }}



      if ( vname=="q" ) {
        qs = apply( y$q[,,], 1, quantile, probs=c(0.025, 0.5, 0.975) )
        qs = signif( qs, 3 )
        for (i in 1:3) {
          pdat = as.vector(y$q[i,,])
        prr=NULL
        prr$class='uniform'
        prr$max=1
        prr$min=0.001
          plot.freq.distribution.prior.posterior( prior=prr, posterior=pdat )
          legend( "topright", bty="n", legend=paste( labs[i], "\n", vname, " = ", qs[2,i], " {", qs[1,i], ", ",  qs[3,i], "}  ", sep="" )   
      )}}

      if ( vname=="qs" ) {
        QQ = apply( y$q[,,], 1, quantile, probs=c(0.025, 0.5, 0.975) )
        QQ = signif( QQ, 3 )
        for (i in 1:3) {
          pdat = as.vector(y$qs[i,,])
          # prr=NULL
          # prr$class="normal"
          # prr$mean=sb$q0x[i]
          # prr$sd=sb$q0x[i]*sb$cv
          plot.freq.distribution.prior.posterior( prior=prr, posterior=pdat )
          legend( "topright", bty="n", legend=paste( labs[i], "\n", vname, " = ", QQ[2,i], " {", QQ[1,i], ", ",  QQ[3,i], "}  ", sep="" )   
      )}}


      if ( vname=="BMSY" ) {
        qs = apply( y$BMSY[,,], 1, quantile, probs=c(0.025, 0.5, 0.975) )
        qs = signif( qs, 3 )
        for (i in 1:3) {
          pdat = as.vector(y$BMSY[i,,])
          prr=NULL
          prr$class="none"
          plot.freq.distribution.prior.posterior( prior=prr, posterior=pdat )
          legend( "topright", bty="n",
            legend=paste( labs[i], " ", vname, " = ", qs[2,i], " {", qs[1,i], ", ",  qs[3,i], "}", sep="" )   
      )}}
  
      if ( vname=="FMSY" ) {
        qs = apply( y$FMSY[,,], 1, quantile, probs=c(0.025, 0.5, 0.975) )
        qs = signif( qs, 3 )
        for (i in 1:3) {
          pdat = as.vector(y$FMSY[i,,])
          prr=NULL
          prr$class="none"
          plot.freq.distribution.prior.posterior( prior=prr, posterior=pdat )
          legend( "topright", bty="n",
            legend=paste( labs[i], " ", vname, " = ", qs[2,i], " {", qs[1,i], ", ",  qs[3,i], "}", sep="" )   
      )}}
	if ( vname=="bo.sd" ) {
        qs = apply( y$bo.sd[,,], 1, quantile, probs=c(0.025, 0.5, 0.975) )
        qs = signif( qs, 3 )
        for (i in 1:3) {
          pdat = as.vector(y$bo.sd[i,,])
          prr=NULL
          prr$class="lognormal"
          prr$meanlog=sb$bo.mup
          prr$sdlog=sqrt(sb$bo.sdp)
          plot.freq.distribution.prior.posterior( prior=prr, posterior=pdat )
          legend( "topright", bty="n",
            legend=paste( labs[i], " ", vname, " = ", qs[2,i], " {", qs[1,i], ", ",  qs[3,i], "}", sep="" )   
      )}}
          if ( vname=="bp.sd" ) {
        qs = apply( y$bp.sd[,,], 1, quantile, probs=c(0.025, 0.5, 0.975) )
        qs = signif( qs, 3 )
        for (i in 1:3) {
          pdat = as.vector(y$bp.sd[i,,])
          prr=NULL
          prr$class="lognormal"
          prr$meanlog=sb$bp.mup
          prr$sdlog=sqrt(sb$bp.sdp)
          plot.freq.distribution.prior.posterior( prior=prr, posterior=pdat )
          legend( "topright", bty="n",
            legend=paste( labs[i], " ", vname, " = ", qs[2,i], " {", qs[1,i], ", ",  qs[3,i], "}", sep="" )   
      )}}

     
    }

    # --------------

    if ( type=="timeseries" ) {
      
      if (vname=="biomass") { 

        SI =  apply( y$q, 1, mean, na.rm=T  )

        for (i in 1:3) {
          qIOA = sb$IOA[,i] / SI[i]
          IOA = sb$IOA[,i] 
          meanval = apply( y$B[,i,,], 1, mean, na.rm=T  )

          prs = seq( from=0.025, to=0.975, length.out=600)
          Bq =  apply( y$B[,i,,], 1, quantile, probs=prs, na.rm=T  )

          yran = range(c(0, Bq, sb$IOA[,i] ), na.rm=T )*1.01
          plot( yrs, Bq[1,], type="n", ylim=yran, xlab="", ylab=""  )
          cols = gray.colors( floor(length( prs)/2) )
          cols2 = c(cols[length(cols):1], cols )
          for ( j in 1:length(prs) ) {
            lines ( yrs, Bq[j,], lwd=4, col=cols2[j] )
          }
          # lines( yrs, B, lwd=3, col="darkgreen" )
          abline (v=yrs.last , lwd=2, lty="dashed" )
          if (i==2) title( ylab="Fishable biomass (kt)" ) 
          if (i==3) title( xlab="Year" ) 
          points( yrs0, qIOA, pch=20, col="darkgreen" )
          lines ( yrs0, qIOA, lwd=3, col="darkgreen", lty="dashed" )
          lines ( yrs, meanval, lwd=2, col="blue", lty="dotted" )
          points( yrs0, IOA, pch=20, col="darkred" )
          lines( yrs0, IOA, lwd=3, lty="dashed", col="red" )
          legend( "topright", bty="n", legend=labs[i])
      }}

      if (vname=="fishingmortality") { 
        Fmsy = apply( y$FMSY, 1, mean, na.rm=T ) 
        for (i in 1:3) {
          prs = seq( from=0.025, to=0.975, length.out=600)
          Fi = apply( y$F[1:sb$N,i,,], 1, quantile, probs=prs, na.rm=T )
          yran = range(c(0, max(c(Fi,Fmsy))), na.rm=T )*1.01
          yran = pmin( yran, 1.2 )
          plot( yrs0, Fi[1,], type="n", ylim=yran, xlab="", ylab="" )
          cols = gray.colors( floor(length( prs)/2) )
          cols2 = c(cols[length(cols):1], cols )
          for ( j in 1:length(prs) ) {
            lines ( yrs0, Fi[j,], lwd=4, col=cols2[j] )
          }
          if (i==2) title( ylab="Fishing mortality" ) 
          if (i==3) title( xlab="Year" ) 
          legend( "topright", bty="n", legend=labs[i])
          abline (h=-log(1-0.2), lwd=2, lty="dashed" )
          abline (h=Fmsy[i], lwd=2, lty="solid", col="red" )
      }}
    }

    if (type=="hcr") {
      if (vname=="default") {
          B =  apply( y$B, c(1,2), mean, na.rm=T  )
          F =  apply( y$F, c(1,2), mean, na.rm=T  )
          K =  apply( y$K, c(1), mean, na.rm=T  )
          FMSY = apply( y$FMSY, c(1), mean, na.rm=T  )
          BMSY = apply( y$BMSY, c(1), mean, na.rm=T  )

        for (i in 1:3 ) {
          ylims = c(0, min( 1, max( FMSY[i] * 1.25, F[hdat,i] ) ) )
          plot( B[hdat,i], F[hdat,i],  type="b", xlim=c(0, K[i] * 1.1 ),  
            ylim=ylims, col="darkorange", cex=0.8, lwd=2, xlab="", ylab="", pch=20 )
 

          # nn = as.matrix( cbind( Bx=as.vector( y$B[ndata,i,,] ), Fx = as.vector( y$F[ndata,i,,] ) ))
          # ellipse.2d(nn[,1], nn[,2], pv=0.05, sc=30)
          
          if (i==3) title( xlab="Fishable biomass (kt)" ) 
          if (i==2) title( ylab="Fishing mortality" ) 

          F30 = -log(1-0.3)
          F10 = -log(1-0.1)

          Fref =  0.22
          Bmsy = K[i] * 0.5
          Bref = K[i] * 0.2 
          BK = K[i] 
          BK25 = K[i] * .25 
          Fhistorical = mean( F[hdat,i], na.rm=T )
          Bhistorical = mean( B[hdat,i], na.rm=T ) 
          yl = 0.05
         
            polygon(x=c(Bmsy,Bmsy*2,Bmsy*2, Bmsy),y=c(-0.1,-0.1,FMSY[i],FMSY[i]),col='lightgreen',border=NA)
          polygon(x=c(Bmsy/2,Bmsy,Bmsy, Bmsy/2),y=c(-0.1,-0.1,FMSY[i],FMSY[i]),col='lightgoldenrod',border=NA)
          polygon(x=c(0,Bmsy/2,Bmsy/2, 0),y=c(-0.1,-0.1,FMSY[i],FMSY[i]),col='darksalmon',border=NA)

        lines( B[hdat,i], F[hdat,i],  type="b", xlim=c(0, K[i] * 1.1 ),  
            ylim=ylims, col="darkblue", cex=0.8, lwd=2, xlab="", ylab="", pch=20 )
        


          abline (h=Fref, lty="solid", col="gray", lwd=2 )
    
          abline (h=F10, lty="dotted", col="gray")
          # text( 0.05*K[i], F10, "10% HR", pos=1 )
          
          abline (h=F30, lty="dotted", col="gray")
          # text( 0.05*K[i], F30, "30% HR", pos=1 )
 

          abline (h=FMSY[i], lty="dashed", col="red" )

          # abline (h=Fhistorical, lty="dashed")
          # text( 0.05*K[i], Fhistorical, "Mean", pos=1, lwd=2 )
          
          # abline (v=Bref, lty="dotted")
          # text( Bref-0.2, 0.25, "Lower biomass reference point\n (LBRP = 0.2 * BMSY)" , srt=90, pos=3)

          abline (v=Bmsy, lty="dotted")
        
          abline (v=BK, lty="dotted")
    
          abline (v=BK25, lty="dotted")
   
          text( Bmsy-0.01*K[i], yl, "K/2" , srt=90, pos=3)
          text( BK-0.01*K[i], yl, "K" , srt=90, pos=3)
          text( BK25-0.01*K[i], yl, "K/4" , srt=90, pos=3)
          text( 0.05*K[i], Fref, "20% HR", pos=1 )
          text( 0.05*K[i], FMSY[i], "FMSY", pos=3, lwd=2, col="red" )
          text( B[hdat,i], F[hdat,i],  labels=yrs0, pos=3, cex= 0.8 )
  
          text( 0, ylims[2]*0.9,  labels=labs[i], pos=3, cex= 0.85 )

          # abline (v=Bhistorical, lty="dashed")
          # text( Bhistorical-0.01*K[i], yl, "Mean" , srt=90, pos=3,  lwd=2)

        }
      }
      
      if (vname=="default.unmodelled") {
      
        B =  sb$IOA
          F =  apply( y$F, c(1,2), mean, na.rm=T  )
            
          areas = c("cfa4x", "cfasouth", "cfanorth" )
          regions = c("4X", "S-ENS", "N-ENS")
    
          td = exploitationrates(p=p, areas=areas, labels=regions, CFA4X.exclude.current.year=FALSE )
   
          K =  apply( y$K, c(1), mean, na.rm=T  )
          FMSY = apply( y$FMSY, c(1), mean, na.rm=T  )
          BMSY = apply( y$BMSY, c(1), mean, na.rm=T  )

        for (i in 1:3 ) {
          ylims = c(0, FMSY[i] * 1.25)
          plot( B[hdat,i], F[hdat,i],  type="b", xlim=c(0, K[i] * 1.1 ),  
            ylim=ylims, col="darkorange", cex=0.8, lwd=2, xlab="", ylab="", pch=20 )

          # nn = as.matrix( cbind( Bx=as.vector( y$B[ndata,i,,] ), Fx = as.vector( y$F[ndata,i,,] ) ))
          # ellipse.2d(nn[,1], nn[,2], pv=0.05, sc=30)
          
          if (i==3) title( xlab="Fishable biomass (kt)" ) 
          if (i==2) title( ylab="Fishing mortality" ) 

          F30 = -log(1-0.3)
          F10 = -log(1-0.1)

          Fref =  0.22
          Bmsy = BMSY[i]
          Bref = K[i] * 0.2 
          BK = K[i] 
          BK25 = K[i] * .25 
          Fhistorical = mean( F[hdat,i], na.rm=T )
          Bhistorical = mean( B[hdat,i], na.rm=T ) 
          yl = 0.05
          
       
          abline (h=Fref, lty="solid", col="gray", lwd=2 )
    
          abline (h=F10, lty="dotted", col="gray")
          # text( 0.05*K[i], F10, "10% HR", pos=1 )
          
          abline (h=F30, lty="dotted", col="gray")
          # text( 0.05*K[i], F30, "30% HR", pos=1 )
 

          abline (h=FMSY[i], lty="dashed", col="red" )

          # abline (h=Fhistorical, lty="dashed")
          # text( 0.05*K[i], Fhistorical, "Mean", pos=1, lwd=2 )
          
          # abline (v=Bref, lty="dotted")
          # text( Bref-0.2, 0.25, "Lower biomass reference point\n (LBRP = 0.2 * BMSY)" , srt=90, pos=3)

          abline (v=Bmsy, lty="dotted")

          abline (v=BK, lty="dotted")
    
          abline (v=BK25, lty="dotted")
        
          text( 0.05*K[i], Fref, "20% HR", pos=1 )
          text( 0.05*K[i], FMSY[i], "FMSY", pos=3, lwd=2, col="red" )
          text( BK-0.01*K[i], yl, "K" , srt=90, pos=3)
          text( Bmsy-0.01*K[i], yl, "K/2" , srt=90, pos=3)
          text( BK25-0.01*K[i], yl, "K/4" , srt=90, pos=3)
          text( B[hdat,i], F[hdat,i],  labels=yrs0, pos=3, cex= 0.8 )
  
          text( 0, ylims[2]*0.9,  labels=labs[i], pos=3, cex= 0.85 )



          # abline (v=Bhistorical, lty="dashed")
          # text( Bhistorical-0.01*K[i], yl, "Mean" , srt=90, pos=3,  lwd=2)

        }
      }
   
      if (vname=="simple") {
        require(car)
     
          B =  apply( y$B, c(1,2), mean, na.rm=T  )
          F =  apply( y$F, c(1,2), mean, na.rm=T  )
          C =  apply( y$C, c(1,2), mean, na.rm=T  )
          K =  apply( y$K, c(1), mean, na.rm=T  )
#          for (i in 1:3) C[,i] = C[,i]  * K[i]
          FMSY = apply( y$FMSY, c(1), mean, na.rm=T  )
          BMSY = apply( y$BMSY, c(1), mean, na.rm=T  )
        
          labs = c("N-ENS", "S-ENS", "4X")

        for (i in 1:3 ) {
          ylims = max(C[,i] )* c(0, 1.1)
          plot( B[hdat,i], C[hdat,i],  type="l", xlim=c(0, K[i]*1.05  ),  
            ylim=ylims, xlab="", ylab="",  lwd=2, col="darkorange" )

          abline(0,0.1, lty="dotted", lwd=2, col="gray" )
          abline(0,0.2, lwd=3, col="gray" )
          abline(0,0.3, lty="dotted", lwd=2, col="gray" )

          points( B[hdat,i], C[hdat,i], col="orange", cex=0.8, pch=20 )
          text( B[hdat,i], C[hdat,i],  labels=yrs0, pos=3, cex= 0.85 )
          
          text( 0, ylims[2]*0.9,  labels=labs[i], pos=3, cex= 0.85 )

          # nn = as.matrix( cbind( Bx=as.vector( y$B[ndata,i,,] ), Fx = as.vector( y$F[ndata,i,,] ) ))
          # ellipse.2d(nn[,1], nn[,2], pv=0.05, sc=30)

          if (i==3) title( xlab="Fishable biomass (kt)" ) 
          if (i==2) title( ylab="Catch (kt)" ) 

          Cmsy = ( exp( FMSY[i] ) - 1) 
          Cref = ( exp( FMSY[i] * 0.2 ) - 1) * K[i]
          Bmsy = BMSY[i]
          Bref = K[i] * 0.2 
          BK = K[i] 
          BK25 = K[i] * .25 
          Chistorical = mean( C[hdat,i], na.rm=T )
          Bhistorical = mean( B[hdat,i], na.rm=T ) 
          yl = 0.1 * max(C[hdat,i])
          
          # abline (h=Fref, lty="dotted")
          # text( 0.25, Fref, "Target\n (0.2 * FMSY) ", pos=1 )

          abline (0, Cmsy, lty="dotted", col="red")
          # text( 0.1*K[i], Cmsy, "FMSY", pos=1 )

          # abline (h=Chistorical, lty="dashed")
          # text( 0.1*K[i], Chistorical, "Mean", pos=3, lwd=2 )
          
          # abline (v=Bref, lty="dotted")
          # text( Bref-0.2, 0.25, "Lower biomass reference point\n (LBRP = 0.2 * BMSY)" , srt=90, pos=3)

        
          abline (v=BK, lty="dotted")
          text( BK-0.01*K[i], yl, "K" , srt=90, pos=3)
  
          abline (v=BK/2, lty="dotted")
          text( BK/2-0.01*K[i], yl, "K/2" , srt=90, pos=3)
    
          abline (v=BK25, lty="dotted")
          text( BK25-0.01*K[i], yl, "K/4" , srt=90, pos=3)

        }


      }
    }

    if (type=="diagnostic.catch") {
      
    }

    
    if (type=="diagnostic.phase") {
      
      B =  apply( y$B, c(1,2), mean, na.rm=T  )
      K =  apply( y$K, c(1), mean, na.rm=T  )
         
      for (i in 1:3 ) {
        plot( B[1:ndata-1,i], B[2:ndata,i],  type="b", xlab="t", ylab="t+1",
          xlim=c(0, K[i] * 1.25 ), ylim=max(K[i] )* c(0, 1.25), lwd=2, col="darkorange" )

         # abline(0,0.1, lty="dotted", lwd=2, col="gray" )
         abline( coef=c(0,1) )
       #text( B[1:ndata-1,], B[2:ndata,],  labels=yrs4 , pos=4, cex=0.8 )
      }
    }


    if (type=="diagnostic.errors") {

      # observation vs process error
      graphics.off()
      layout( matrix(c(1,2,3), 3, 1 ))
      par(mar = c(5, 4, 0, 2))
      require(car)
      
      eP = y$bp.sd
      eO = y$bo.sd
      for (i in 1:3 ) {
          plot( eP[,i,,], eO[,i,,],  type="p", pch=22 ) 
          if (i==2) title( ylab="Process error (SD)" ) 
          if (i==3) title( xlab="Observation error (SD)" ) 

      }


    }


    if (type=="diagnostic.production") {
       
      B =  apply( y$B, c(1,2), mean, na.rm=T  )
      P =  apply( y$P, c(1,2), mean, na.rm=T  )

      MSY = apply( y$MSY, c(1), mean, na.rm=T  )
      FMSY = apply( y$FMSY, c(1), mean, na.rm=T  )
      BMSY = apply( y$BMSY, c(1), mean, na.rm=T  )
      C =  apply( y$C, c(1,2), mean, na.rm=T  )
      K =  apply( y$K, c(1), mean, na.rm=T  )
         
      # production vs biomass
      x11()
      layout( matrix(c(1,2,3), 3, 1 ))
      par(mar = c(5, 4, 0, 2))
      for (i in 1:3) {
        plot( B[,i], P[,i], type="n", pch=20, ylim=c(0, max( c(P[,i], MSY[i]))*1.1), xlim=c(0,K[i]*1.05), xlab="Biomass; kt", ylab="Yield; kt"  )
        a = MSY[i] / (BMSY[i])^2 
        curve( -a*(x-BMSY[i])^2 + MSY[i], from=0, to=K[i], add=TRUE, lwd=3, col="gray" )
        abline(v=BMSY[i], lty="dotted", lwd=3, col="gray")
        points( B[,i], P[,i], type="p", pch=20  )
        text( B[,i], P[,i], yrs0, pos=3, cex=0.8 )
        # abline(h=0)
      }
   
    }


    if (is.null(fn)) fn = paste(vname, "tmp", ".png", sep="" )
  if(save.plot) savePlot( filename=fn, type="png" )
    return( fn)

  }


 


