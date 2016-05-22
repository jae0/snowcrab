

  logbook.db = function( DS, prorate=T, p=NULL, yrs=NULL ) {
    		
		if (DS %in% c("odbc.logbook", "odbc.logbook.redo")) {
			fn.root =  file.path( project.datadirectory("snowcrab"), "data", "logbook", "datadump" )
			dir.create( fn.root, recursive = TRUE, showWarnings = FALSE )
				
			if (DS=="odbc.logbook") {
				out = NULL
				for ( YR in yrs ) {
					fny = file.path( fn.root, paste( YR, "rdata", sep="."))
					if (file.exists(fny)) {
						load (fny)
						out = rbind( out, logbook )
					}
				}
				return (out)
			}
         
			require(RODBC)
      con=odbcConnect(oracle.snowcrab.server , uid=oracle.snowcrab.user, pwd=oracle.snowcrab.password, believeNRows=F)

			for ( YR in yrs ) {
				fny = file.path( fn.root, paste( YR,"rdata", sep="."))
				query = paste(
					"SELECT * from marfissci.marfis_crab ", 
					"where target_spc=705", 
					"AND EXTRACT(YEAR from DATE_LANDED) = ", YR )
				logbook = NULL
				logbook = sqlQuery(con, query )
				save( logbook, file=fny, compress=T)
				gc()  # garbage collection
				print(YR)
			}
			odbcClose(con)
			return (yrs)
			
		}


    if (DS %in% c("odbc.licence.redo", "odbc.licence" ) ) {
               
      filename.licence = file.path( project.datadirectory("snowcrab"), "data", "logbook", "lic.datadump.rdata" )

      if (DS=="odbc.licence") {
        load(filename.licence)
        return (lic)
      }

      require(RODBC)
      con=odbcConnect(oracle.snowcrab.server , uid=oracle.snowcrab.user, pwd=oracle.snowcrab.password, believeNRows=F)
      
      lic = sqlQuery(con, "select * from marfissci.licence_areas")
      save(lic, file=filename.licence, compress=T) 
     
		}


    if (DS %in% c("odbc.areas.redo", "odbc.areas" ) ) {
               
      filename.areas = file.path( project.datadirectory("snowcrab"), "data", "logbook", "areas.datadump.rdata" )

      if (DS=="odbc.areas") {
        load(filename.areas)
        return (areas)
      }

      require(RODBC)
      con=odbcConnect(oracle.snowcrab.server , uid=oracle.snowcrab.user, pwd=oracle.snowcrab.password, believeNRows=F)
      areas = sqlQuery(con, "select * from marfissci.areas")
      save(areas, file=filename.areas, compress=T) 
      
      return ("Complete")

    }
    

  # -------

    if (DS %in% c("logbook.filtered.positions", "logbook.filtered.positions.redo")) {
    
      # exclude data that have positions that are incorrect

      filename = file.path( project.datadirectory("snowcrab"), "data", "logbook", "logbook.filtered.positions.rdata" )

      if (DS=="logbook.filtered.positions") {
        load( filename )
        return(lgbk)
      }

      lgbk = logbook.db( DS="logbook" )

      h = which( !is.na( lgbk$lon + lgbk$lat ) ) 
      lgbk = lgbk[h,]
      
      i = filter.region.polygon( lgbk[, c("lon", "lat")], region="cfaall")
      lgbk = lgbk[i,]

      j = filter.region.polygon( lgbk[, c("lon", "lat")], region="isobath1000m")
      lgbk = lgbk[j,]
      
      # additional constraint ..  
      # only accept "correctly" positioned data within each subarea ... in case miscoded data have a large effect
     
        icfa4x = filter.region.polygon( lgbk[, c("lon", "lat")], "cfa4x")
        icfanorth = filter.region.polygon( lgbk[, c("lon", "lat")], "cfanorth")
        icfa23 = filter.region.polygon( lgbk[, c("lon", "lat")], "cfa23")
        icfa24 = filter.region.polygon( lgbk[, c("lon", "lat")], "cfa24")
      
      gooddata = sort( unique( c(icfa4x, icfanorth, icfa23, icfa24 ) ) )
      lgbk = lgbk[gooddata, ]
      
      save( lgbk, file=filename, compress=T )

      return(filename)

    }


  # -------


    if (DS %in% c("logbook", "logbook.redo")) {

      filename = file.path( project.datadirectory("snowcrab"), "data", "logbook", "logbook.rdata" )

      if (DS=="logbook") {
        load( filename )
        return(logbook)
      }
    
			
			# stop(" MUST ADD EXTERNAL DATA -- landings in GULF ")
      
			
			# logbooks from historical data tables (1996-2003; but 2002 and 2003 seem to be partial records) 
      lb.historical = logbook.db( DS="fisheries.historical" )
      lb.historical$cfa = NA
      lb.historical$date.fished = paste( lb.historical$date.fished, "01:00:00" )
      lb.historical$date.landed = as.POSIXct( lb.historical$date.landed  )

      # logbooks from marfissci tables 
      x = logbook.db( DS="odbc.logbook", yrs=1996:p$current.assessment.year )
 
      names( x ) = tolower( names( x ) )
      names( x ) = rename.snowcrab.variables(names( x ))
      to.char = c( "cfa", "licence_id_old", "licence", "date.fished", "date.landed",
                   "captain", "vessel", "doc_id", "doc_type")
      for (i in to.char) x[,i] = as.character(x[,i])

      iy = which(!is.finite(x$year))
      if (length(iy) > 0) {
        x$year[iy] = years(x$date.landed) [iy]
        iy = which(!is.finite(x$year))
        if (length(iy) > 0)  x = x[ -iy, ] 
      }

      prorate=FALSE

      # if (prorate) x = logbook.prorate( x )  # assume pro-rating not required for historical data as it was maintained manually by Moncton
      
      x$soak.time = x$soak_days * 24  # make into hours
      x$trap.type = "" # dummy value until this can be added to the views .. check with Alan
      x$status = ""  # "" "" ""

      datelanded =  matrix(unlist(strsplit(x$date.landed, " ", fixed=T)), ncol=2, byrow=T)[,1]
      x$date.landed = lubridate::ymd( datelanded )
      x$landings = x$pro_rated_slip_wt_lbs * 0.454  # convert to kg
      x$cpue = x$landings / x$effort
      x$depth = x$depth_fm*1.83

      x$lat =   round( as.numeric(substring(x$lat, 1,2)) + as.numeric(substring(x$lat, 3,6))/6000 ,6)
      x$lon = - round((as.numeric(substring(x$lon, 1,2)) + as.numeric(substring(x$lon, 3,6))/6000), 6)
    
      to.extract = c( "year","lat","lon","depth","landings","effort","soak.time",
                      "cpue","trap.type","cfv","status","licence",
                      "date.landed", "date.fished", "cfa" )
      
      lb.marfis = x[, to.extract ]
      
      lb.historical$date.landed = as.POSIXct(lb.historical$date.landed)
      a = lb.historical$date.fished
      b = paste(substr(a, 1, 4), "/", substr(a, 5, 6), "/", substr(a, 7, 8), sep="")
      lb.historical$date.fished = as.POSIXct(b)

      
      x = NULL
      x = rbind( lb.historical, lb.marfis )
      
      dups = which(duplicated(x))
      toremove = sort(unique(c(iy, dups)))
      if (length(toremove) > 0) x = x[-toremove,]


      # known errors:  manual fixes
      ix = which( round(x$lat)==46 & round(x$lon)==-5930 )

      if ( length(ix > 0 ))  x$lon[ ix ]  = x$lon[ ix ] / 100
   
      x = logbook.determine.region(x)  # using licence info and geographics

      i.cfa4x = which( x$cfa == "cfa4x" )
      i.offset = which( lubridate::month(x$date.landed) >= 1 & lubridate::month(x$date.landed) <= 7 )
      to.offset = intersect( i.cfa4x, i.offset)

      x$yr = x$year
      x$yr[to.offset] = x$yr[to.offset] - 1
      # x$yr[i.cfa4x] = x$yr[i.cfa4x] + 1  # ie:: fishery from 1999-2000 in 4X is now coded as 2000
      
      a= x[which(x$cfa0=='cfa4x'),]
      head(a)
     
      # enforce bounds in effort and cpue
      oo = which( x$cpue > 650 * 0.454 )  # 600 - 650 lbs / trap is a real/reasonable upper limit
      if ( length(oo) > 0 ) x$cpue[oo] = NA  

      pp = which ( x$effort > 240 ) # small traps were used at times with large 240 trap compliments 
      if ( length(pp) > 0 ) x$effort[pp] = NA

      logbook = x
      
      save (logbook, file=filename, compress=T )  # this is for plotting maps, etc

      return( "Complete" )
    
    }

  # ---------------------

    if (DS %in% c( "fishing.grounds.annual", "fishing.grounds.global", "fishing.grounds.redo")) {
      
      fn1 = file.path(  project.datadirectory("snowcrab"), "data", "logbook", "fishing.grounds.global.rdata")
      fn2 = file.path(  project.datadirectory("snowcrab"), "data", "logbook", "fishing.grounds.annual.rdata")
      
      if (DS=="fishing.grounds.global") {
        load( fn1 )
        return (fg) 
      }

      if (DS=="fishing.grounds.annual") {
        load( fn2 )
        return (fg) 
      }

      out = NULL
      fg.res = p$fisheries.grid.resolution
      
      x = logbook.db( DS="logbook.filtered.positions" )
    
      x = lonlat2planar( x,  proj.type=p$internal.projection ) 
      x$plon = grid.internal( x$plon, p$plons )
      x$plat = grid.internal( x$plat, p$plats )
      yrs = c(T,F)
      for ( Y in yrs ) {
        if (Y) {
          fn = fn2
          x$gridid = paste(x$plon%/%fg.res*fg.res, x$plat%/%fg.res*fg.res, x$year, sep="." )
          ncols=3
        } else {
          fn = fn1
          x$gridid = paste(x$plon%/%fg.res*fg.res, x$plat%/%fg.res*fg.res, sep="." )
          ncols=2
        }
        v = "visits"
        x$visits=1
        w = x[is.finite(x[,v]),]
        tmp = as.data.frame( xtabs( as.integer(x[,v]) ~ as.factor(x[,"gridid"]), exclude="" ) )
        names(tmp) = c("gridid", "total.visits")
        out = tmp

        v = "landings"
        w = x[is.finite(x[,v]),]
        tmp = as.data.frame( xtabs( as.integer(x[,v]) ~ as.factor(x[,"gridid"]), exclude="" ) )
        names(tmp) = c("gridid", "total.landings")
        out = merge( out, tmp, by="gridid", all=T, sort=F)

        v = "effort"
        w = x[is.finite(x[,v]),]
        tmp = as.data.frame( xtabs( as.integer(x[,v]) ~ as.factor(x[,"gridid"]), exclude="" ) )
        names(tmp) = c("gridid", "total.effort")
        out = merge( out, tmp, by="gridid", all=T, sort=F)

        out$total.cpue = out$total.landings / out$total.effort

        tmp = matrix(unlist(strsplit(as.character(out$gridid), ".", fixed=T)), ncol=ncols, byrow=T)
        out$plat = as.numeric(tmp[,2])
        out$plon = as.numeric(tmp[,1])

        out$gridid = as.character( out$gridid )

        if (Y) out$yr = as.numeric(tmp[,3])
        
        fg = out
        save( fg, file=fn, compress=T )
      }

      return( "Complete")
    }


  # -----------------------------


    if (DS %in% c("fisheries.historical", "fisheries.historical.redo" )) {

      fn = file.path( project.datadirectory("snowcrab"), "data", "logbook", "logbook.historical.rdata" ) 

      if (DS=="fisheries.historical") {
        logs = NULL
        if (file.exists(fn)) load( fn)
        return( logs )
      }
    
      historicaldataloc = file.path( project.datadirectory("snowcrab"), "data", "logbook", "archive") 
      files = c("logbooks1996.csv", "logbooks1997.csv", "logbooks1998.csv", "logbooks1999.csv", "logbooks2000.csv", "logbooks2001.csv")
      out = NULL
      for (g in files) {
        f = file.path(historicaldataloc, g )
        a = read.table(file=f, skip=1, as.is=T, strip.white=T, sep=";")
        a = a[,c(1:14)]
        names(a) = c("cfv", "areas", "status", "licence", "date.landed", "date.fished", "lat", "lon", 
            "landings.kg", "landings.lbs", "effort", "soak.time", "cpue", "trap.type" )


        a$file = f
        a$yr = as.numeric(gsub("[[:alpha:]]", "", f))
        out = rbind(out, a)
      }

       f4x = read.table(file=file.path(historicaldataloc, "logbooks4x.csv"), skip=1, as.is=T, strip.white=T, sep=";")
       names(f4x) = c("cfv", "areas",  "date.landed", "date.fished", "lat", "lon", 
              "landings.kg", "landings.lbs", "effort", "soak.time", "cpue", "trap.type" )
       f4x$areas="4X"
       f4x$status = NA
       f4x$licence = NA
       f4x$file = "logbooks4x.csv"
       f4x$yr = trunc(f4x$date.landed/10000)
       f4x = f4x[, names(out) ]
       
       logs = rbind (out, f4x)
       logs$lon = -logs$lon
        
       logs$depth = NA
       logs$year = logs$yr
       logs$landings = logs$landings.kg
       logs$cpue = logs$landings / logs$effort

       i = which( nchar(logs$date.landed) <10 & !is.na(logs$date.landed) )

       dt = logs[i, "date.landed"]
       yr = substring(dt,1,4)
       mon = substring(dt,5,6)
       da = substring(dt,7,8)

       logs[i, "date.landed"] = paste( mon, da, yr, sep="/")
       logs$date.landed = lubridate::mdy( logs$date.landed)

       to.extract = c( "year","lat","lon","depth","landings","effort","soak.time",
                        "cpue","trap.type","cfv","status","licence",
                        "date.landed", "date.fished")
       logs = logs[, to.extract]
       logs = logs[ -which( !is.finite(logs$year) ), ]

       save(logs, file=fn, compress=T)

       return( "completed")
    
    }   # end if historical


    if (DS %in% c("fisheries.complete", "fisheries.complete.redo" )) {

      fn = file.path( project.datadirectory("snowcrab"), "data", "logbook", "logbook.complete.rdata" ) 

      if (DS=="fisheries.complete") {
        load( fn)
        return( logbook )
      }

			logbook = logbook.db(DS="logbook.filtered.positions")

      nl0 = nrow( logbook ) 
      logbook = lonlat2planar( logbook,  proj.type=p$internal.projection )
      
      logbook$plon = grid.internal( logbook$plon, p$plons )
      logbook$plat = grid.internal( logbook$plat, p$plats )

			logbook$timestamp = as.POSIXct( logbook$date.landed, tz="America/Halifax", origin=lubridate::origin  )  # required for temperature lookups
      logbook$dyear = lubridate::decimal_date( logbook$timestamp ) - lubridate::year(logbook$timestamp )

      logbook = logbook[which(logbook$yr<=p$current.assessment.year),]
			# bring in time invariant features:: depth
      logbook$z = logbook$depth
			logbook$depth = NULL
      oo =  which( logbook$z < 10 | logbook$z > 500 ) # screen out large z's
      if (length(oo) > 0 )  logbook$z[ oo ] = NA  
			logbook = habitat.lookup( logbook, p=p, DS="depth" )
      logbook$z = log( logbook$z )
			
      ii = which( ! is.finite( logbook$z) )  
      if (length(ii)>0) logbook = logbook[ -ii, ]

		  # bring in time varing features:: temperature
			logbook$t = NA
      logbook$t = habitat.lookup( logbook, p=p, DS="temperature" )

			# bring in habitat variables
			logbook = habitat.lookup( logbook, p=p, DS="all.data" )
      logbook$z = log(logbook$z) 

			save( logbook, file=fn, compress=T )

      return  ("Complete")
    } 
  
  
    if (DS %in% c("logbook.gridded",  "logbook.gridded.redo" ) ) {
      
      loc = file.path( project.datadirectory("snowcrab"), "data", "logbook", "gridded.fishery.data" )
      dir.create( path=loc, recursive=T, showWarnings=F)

      if (DS == "logbook.gridded") {
        fn = file.path(loc, paste( "gridded.fishery", yrs, "rdata", sep=".") )
        load(fn)
        return(gridded.fishery.data)
      }
     
      yy = logbook.db(DS="logbook")
      yrs = sort( unique( yy$year))

      for ( y in yrs ) {

        fn = file.path(loc, paste( "gridded.fishery", y, "rdata", sep=".") )
        print (fn)

        # load logbook info: global summary
        fg = logbook.db( DS="fishing.grounds.global" )  # in dataframe fg
        fg0 = regrid.lonlat(old=fg, res=p$fisheries.grid.resolution, vr.to.sum=c("total.landings", "total.effort", "total.visits"))
        fg0$core.visits0 = ifelse( fg0$total.visits >= 5, 1, 0 )
        fg0$total.cpue = log( fg0$total.landings / fg0$total.effort + 1)
        fg0$total.landings = log( fg0$total.landings+1 )
        fg0$total.effort = log( fg0$total.effort+1 )
        fg0$total.visits = log( fg0$total.visits+1 )
        fg0 = fg0[ , c( "gridid", "core.visits0", "total.effort", "total.cpue", "total.landings", "total.visits" ) ]
        rm(fg)

        # load logbook info: annual
        fg = logbook.db( DS="fishing.grounds.annual"  )  # in dataframe fg
        fg = fg[ which(fg$yr==y),]
        fg = regrid.lonlat(old=fg, res=p$fisheries.grid.resolution, vr.to.sum=c("total.landings", "total.effort", "total.visits") )
        fg$core.visits = ifelse(fg$total.visits >= 3, 1, 0)
        fg$core.visits = ifelse(fg$total.visits >= 3, 1, 0)
        fg$core.landings = ifelse(fg$total.landings >= 5, 1, 0)
        fg$core.effort = ifelse(fg$total.effort >= 100, 1, 0)
        fg$cpue = log( fg$total.landings / fg$total.effort + 1)
        fg$landings = log( fg$total.landings+1 )
        fg$effort = log( fg$total.effort+1 )
        fg$visits = log( fg$total.visits +1 )

        fg = fg[ , c( "gridid", "core.visits", "core.landings", "core.effort", "cpue", "landings", "effort", "visits" ) ]
        fg = fg[ is.finite(fg$cpue * fg$landings * fg$effort),]
        fg = fg[ which(fg$core.visits==1) , ]

        gridded.fishery.data = merge(fg0, fg, by="gridid", all.x=T, all.y=T, sort=F)
        save( gridded.fishery.data, file=fn, compress=T)
      }
    } # end gridded fishieries data
 
  
  }


