
  net.configuration = function( N, t0=NULL, t1=NULL, set_timestamp=NULL, yr=NULL, plotdata=TRUE ) {

    # N is netmind data
    # t0 is current best estimate of start and end time
    # set_timestamp is timestamp from set-log .. alternate if there is no other useful time marker
    # note:: NETMIND data stream calls the wingspread as "DOORSPREAD".
    # It is not. It is wingspread.

    plotdir = project.datadirectory("bio.snowcrab", "data", "netmind", "figures" )

    # create default output should the following fail
    out = data.frame( slon=NA, slat=NA, distance=NA, spread=NA, spread_sd=NA,
      surfacearea=NA, vel=NA, vel_sd=NA, netmind_n=NA, t0=NA, t1=NA, dt=NA, yr=NA )

    n.req = 10

    if (!is.null(t0)) {
      if (length(t0)>1) {t0 = NULL}
      if (is.na(t0)) t0 = NULL
    }

    if (!is.null(t1)) {
      if (is.na(t1)) t1 = NULL
    }

    #changed this switch from depgth filed to lat as if there is not depth info can still run script
	  if ( length( which( is.finite( N$lat))) < n.req ) print(N[1,])

    problem = F

    # time checks
    if ( is.null(t0) ){
      if ( !is.null(set_timestamp) ) {
        t0 = set_timestamp  # no data in t0 ,, use set_timestamp as alternate
        set_timestamp =NULL # remove to cause no more effects
      }
    }

    #bad.list = c('netmind.S26092014.9.541.17.48.304',
    #             'netmind.S20092007.8.333.17.27.241' )
    bad.list = NULL
    bad.list = unique( c(bad.list, p$netmensuration.problems) )
    if ( N$netmind_uid[1] %in% bad.list) {
      return(out)
    }

    if ( any( is.null( t1 ) || is.null(t0) ) )  {
      # try to determine from netmind data if minilog/seadbird data methods have failed. .. not effective due to noise/and small data stream

      M = N[, c("timestamp", "depth") ]

      oo = which( is.finite( M$depth ))
      if ( length(oo) < n.req ) return(out)

      if(!is.null(set_timestamp)) settimestamp=set_timestamp
      if(is.null(set_timestamp)) settimestamp=t0
      time.gate =  list( t0=settimestamp - dminutes(6), t1=settimestamp + dminutes(12) )

      bcp = list(id=N$netmind_uid[1], nr=nrow(M), YR=yr, tdif.min=3, tdif.max=11, time.gate=time.gate,
                   depth.min=20, depth.range=c(-25,15), eps.depth=3,
                   smooth.windowsize=5, modal.windowsize=5,
                   noisefilter.trim=0.05, noisefilter.target.r2=0.8, noisefilter.quants=c(0.05, 0.95) )

      if(yr<2007)bcp$from.manual.archive=FALSE # manual touchdown only done since 2007


      bcp = bottom.contact.parameters( bcp ) # add other default parameters
      bcp$user.interaction = FALSE
      bc = NULL
      bc = bottom.contact( x=M, bcp=bcp )
      if ( is.null(bc) || (!is.null( bc$res) && ( ( !is.finite(bc$res$t0 ) || !is.finite(bc$res$t1 ) ) ) )) {
         bcp$noisefilter.target.r2=0.8
         if(bcp$id=="netmind.S20092011.14.527.0.47.334") browser()
         bc = bottom.contact( x=M, bcp=bcp )
      }
      if ( is.null(bc) || (!is.null( bc$res) && ( ( !is.finite(bc$res$t0 ) || !is.finite(bc$res$t1 ) ) ) )) {
        # try once more with random settings
        bcp$noisefilter.inla.h =  0.05
        bcp$noisefilter.target.r2= 0.75
        bc = bottom.contact( x=M, bcp=bcp )
      }
      if ( is.null(bc) || (!is.null( bc$res) && ( ( !is.finite(bc$res$t0 ) || !is.finite(bc$res$t1 ) ) ) )) {
        # try once more with more random settings .. noise is high in netminds data
        bcp$eps.depth = 5
        M$depth = jitter( M$depth, amount = bcp$eps.depth/10 )
        bcp$noisefilter.inla.h =  0.1
        bc = bottom.contact( x=M, bcp=bcp )
      }
      if ( !is.null(bc) )  {
        if (plotdata) {
          bottom.contact.plot( bc )
          plotfn = file.path( plotdir, paste(bcp$id, "pdf", sep="." ) )
          print (plotfn)
          dev.flush()
          dev.copy2pdf( file=plotfn )
          graphics.off()
        }
      }
      if (is.null(t0) & !is.null(bc$bottom0) ) t0 = bc$bottom0
      if (is.null(t1) & !is.null(bc$bottom1) ) t1 = bc$bottom1
      N = N[ bc$bottom.contact , ]
    }

    if (all(is.na(t0))) t0=NA
    if (all(is.na(t1))) t1=NA

    if ( is.null(t1) || !is.finite(t1) ) {
      t1_tmp = NA  # rather than guessing, flag and then fill later
    } else {
      t1_tmp = t1
    }

    # if we are here, it is because a reasonable amount of data is present ..
    # do some more checks and get a first estimate of some parameters in case other errors/limits are found
    out$slon=N$lon[1]
    out$slat=N$lat[1]
    out$spread=mean( N$doorspread, na.rm=T ) / 1000  # though called doorspread it is actually wingspread
    out$spread_sd=sd( N$doorspread, na.rm=T ) /1000

    # first set the impossible door spreads to NA and then interpolate based upon a 5 and 95 % quantiles
    # Netmind has very noisy data
    hl = c( 3, 16 ) # m
    ihl = which( N$doorspread < hl[1] | N$doorspread > hl[2] )
    if (length (ihl)>0 ) N$doorspread[ihl] = NA
    quantiles.to.trim = c(0.05, 0.95)
    qnt =  quantile( N$doorspread[N$doorspread>0], quantiles.to.trim, na.rm=T)
    iqnt = which( N$doorspread < qnt[1] | N$doorspread > qnt[2] )
    N$doorspread[ iqnt ] = NA
    gooddoor =  which(   is.finite(N$doorspread) )
    baddoor =   which( ! is.finite(N$doorspread) )

    if (length(baddoor)>0) N$doorspread[ baddoor ] = NA
    if ( length( gooddoor) < n.req ) {
      #problem =T turned this off December 23, 2013 as if we have decent gps data now we can calc distance
      if ( length( gooddoor) >0 ) {
      	out$netmind_n <- length(gooddoor)
        out$spread = mean( N$doorspread[ gooddoor ], na.rm=T ) / 1000
        out$spread_sd = sd( N$doorspread[ gooddoor ], na.rm=T ) / 1000
      }
    }

    if ( length(t0)>1) t0 = NA
    out$t0 = t0
    out$t1 = t1

    out$dt = difftime( as.POSIXct(t1), as.POSIXct(t0) ) # minutes

    if(!is.na(out$t0))    out$yr = lubridate::year( out$t0)
    if(is.na(out$t0))    out$yr = lubridate::year(N$timestamp[1])

    itime =  which( N$timestamp >= t0  &  N$timestamp <= t1 )
    if ( length( itime) < n.req ) problem = T

    if (problem) {
      out$t0 = t0
      out$t1 = NA
      out$dt = NA
      return(out)
    }

    # eOR checks

    
      N = N[ itime, ]
      
      # Higher GPS data accuracy in 2019 resulted in high approximation of distance travelled. This was due to frequent occasions of backwards
      # motion (assuming from boat roll and pitch) which was added to the total distance. When looking at previous years this also seem to have
      # but not to the same extent however removing the backwards motion would result in better data accuracy from year using e-sonar 2014+.
      # The solution is to trim the data by 3/4, or, with the current recording rate of every 4 seconds instead of every 1 second.
      if(yr>=2014){
        N = N[unique(c(seq(1, nrow(N), by = 4), nrow(N))),]
      }
      
      
      if(nrow(N)>1) {

        # t1 gives the approximate time of net lift-off from bottom
        # this is not good enough as there is a potential backdrift period before net lift off
        # this means distance trawled calculations must use the geo-positioning of the boat
        # at maximum tension and not the position at net movemnent up
        # (otherwise, this would introduce a net bias towards smaller swept areas).
        pos = c( "lon", "lat" )
        distance.from.start = geosphere::distCosine(N[1, pos], N[, pos])/1000  #  in km^2
        end = which.max( distance.from.start)
        if( !is.finite(end) ) end=nrow(N)
        n = N[ 1:end , ]

    
	      out$vel = mean(n$speed, na.rm=T, trim=0.1)
	      out$vel_sd = sd(n$speed, na.rm=T)
	      out$slon=n$lon[1]
	      out$slat=n$lat[1]
	      out$netmind_n=end

      delta.distance = NULL
	    n$distance = NA

	  
    if (nrow(n) > 10) {
      # integrate area:: piece-wise integration is used as there is curvature of the fishing track (just in case)
      delta.distance = geosphere::distGeo ( n[ 1:(end-1), pos ], n[ 2:end, pos ] ) / 1000 ## in meters convert to km
      n$distances = c( 0, cumsum( delta.distance  ) ) # cumsum used to do piecewise integration of distance

      
      # model/smooth/interpolate the spreads
      n$doorspread.predicted = NA
      ii = which( is.finite( n$doorspread ) )
      if ( length(ii) > n.req ) {
          # recall that doorspread is actually wingspread ...
         n$doorspread.predicted = approx( x=n$distances, y=n$doorspread, xout=n$distances, method="linear", rule=2 )$y
       #turned off gam model in December 20, 2013 giving unrealistic values for spread as the
       # new esnoar files have 0 and NA whereas older netmind are filled with previous value
			      	#gam.model = try( gam( doorspread ~ s(distances, k=5, bs="ts"), data=n[ii,], optimizer=c("outer", "nlm")), silent = T )
			        #if ( ! "try-error" %in% class( gam.model )) {
			        #  n$doorspread.predicted = predict( gam.model, newdata=n, newdata.guaranteed=T )
			        #} else {
			        #  n$doorspread.predicted = approx( x=n$distances, y=n$doorspread, xout=n$distances, method="linear", rule=2 )$y
			        #}
      }
      if ( length( which( is.finite( n$doorspread.predicted ) ) ) < 10 ) {
        n$doorspread.predicted = mean( n$doorspread , na.rm=T, trim=0.1 )
      }
      mean.doorspreads = ( n$doorspread.predicted[1:(end-1)] + n$doorspread.predicted[2:(end)] ) / 2 / 1000  # mean between two ends
      partial.area =  delta.distance * mean.doorspreads
      out$surfacearea = sum( partial.area )  # km^2
      out$surfacearea = abs(  out$surfacearea )

      out$spread = mean(n$doorspread.predicted, na.rm=T, trim=0.1)/1000  # in km
     spread_sd = sd(n$doorspread.predicted, na.rm=T )/1000
     if(!is.na(spread_sd) & spread_sd!=0) out$spread_sd = spread_sd #if just using the mean from above do not over write spread_sd
     out$distance=n$distances[end]
   }
	}
    return (out)
  }


