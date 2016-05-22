  observer.db = function( DS, p=NULL, yrs=NULL ) {
		

    if (DS %in% c("odbc.redo", "odbc") ) {
			
			if (  Sys.info()["sysname"] == "Windows" ) {
				.Library.site <- "D://R//library-local"
				.libPaths("D://R//library-local")
			}

			fn.root =  file.path( project.datadirectory("snowcrab"), "data", "observer", "datadump" )
			dir.create( fn.root, recursive = TRUE, showWarnings = FALSE )
		
			if (DS=="odbc") {
				out = NULL
				for ( YR in yrs ) {
					fny = file.path( fn.root, paste( YR, "rdata", sep="."))
					if (file.exists(fny)) {
						load (fny)
						out = rbind( out, odb )
					}
				}
				return (out)
			}

      # for the full list of tables:
      # tbls = sqlTables(connect)
      # gs.tables = tbls[ which(tbls[,2] == "GROUNDFISH"),]
      # print(gs.tables)

      print ("Run on windows sessions as the linux odbc session is not happy")
      require(RODBC)
      con=odbcConnect(oracle.snowcrab.server , uid=oracle.snowcrab.user, pwd=oracle.snowcrab.password, believeNRows=F)
      
			for ( YR in yrs ) {
				fny = file.path( fn.root, paste( YR,"rdata", sep="."))
				odbq = paste(
        "SELECT s.LATITUDE, s.LONGITUDE, s.LANDING_DATE, s.SET_NO, s.PRODCD_ID, s.EST_CATCH," ,
               "s.NUM_HOOK_HAUL, d.BOARD_DATE, d.FISH_NO, d.SEXCD_ID, d.FISH_LENGTH, " ,
               "d.FEMALE_ABDOMEN, d.CHELA_HEIGHT, d.SHELLCOND_CD, d.DUROMETRE, d.TRIP_ID, d.TRIP  " ,
         "FROM SNOWCRAB.SNCRABDETAILS_OBS d, SNOWCRAB.SNCRABSETS_OBS s " ,
         "WHERE d.TRIP_ID = s.TRIP_ID  " ,
         "AND d.SET_NO = s.SET_NO  " ,
         "AND d.FISH_NO Is Not Null" ,
				 "AND EXTRACT(YEAR from d.BOARD_DATE) = ", YR )
				odb = NULL
				odb = sqlQuery(con, odbq )
				save( odb, file=fny, compress=T)
				gc()  # garbage collection
				print(YR)
			}
			odbcClose(con)
			return (yrs)
			
    }

    # ---------------------

    if (DS %in% c("odb", "odb.redo")) {
      fn = file.path( project.datadirectory("snowcrab"), "data", "observer", "odb.rdata" )
      if (DS=="odb") {
        load( fn )
        return(odb)
      }

      mod1 = allometry.snowcrab ( "cw.mass", "male")
      mod2 = allometry.snowcrab ( "chela.mass", "male")
      mod3 = allometry.snowcrab ( "cw.chela.mat", "male")
      
      odb = observer.db( DS="odbc", yrs=1996:p$current.assessment.year )
      names(odb) = rename.snowcrab.variables( names(odb) )

      i.m = which( odb$sex==1)
      i.f = which( odb$sex==2)
      i.o = setdiff( 1:nrow( odb ), c(i.m, i.f) )
      odb$sex[ i.m ] = male
      odb$sex[ i.f ] = female
      odb$sex[ i.o ] = sex.unknown

      odb$cw[ which(odb$cw < 40 | odb$cw > 185)  ] = NA
      odb$log.cw = log( odb$cw )
      odb$mass =NA
      odb$log.mass = predict( mod1, odb )
      kk = intersect( which( !is.finite( odb$mass ) ), i.m )
      odb$mass[kk] = exp( odb$log.mass[kk] )

      kk = intersect( which( !is.finite( odb$mass ) ), i.m )
      odb$log.chela = log( odb$chela )
      log.mass2 = predict( mod2, odb )
      odb$mass[kk] = exp( log.mass2[kk] )
      rm (mod1, mod2)


      odb$mat_tmp = predict( mod3, odb, type="response" )
      odb$mat = NA
      odb$mat [ which(odb$mat_tmp <0.5 & odb$sex==male) ] = immature
      odb$mat [ which(odb$mat_tmp >=0.5 & odb$sex==male) ] = mature

      odb$log.cw = NULL
      odb$log.chela = NULL
      odb$log.mass = NULL
      odb$mat_tmp = NULL

      odb$lon = -odb$lon
      odb$timestamp =  odb$sdate
      odb$chron = as.chron(odb$sdate)
      odb$yr = as.numeric(as.character(years(odb$chron)))

      # cfa 4X has a fishing season that spans two years recode "yr" to "fishyr" to accomodate this
      cfa4x = filter.region.polygon(odb, recode.areas("cfa4x"))
      to.offset = which( as.numeric(lubridate::month(odb$chron)) >= 1 & as.numeric( lubridate::month(odb$chron)) <= 7 )
      to.offset = sort(intersect(cfa4x, to.offset))
      odb$fishyr = odb$yr
      odb$fishyr[to.offset]  = odb$fishyr[to.offset] - 1

      odb = odb[odb$fishyr >= 1996 ,] # years for which observer database are good
      #  odb$cw[odb$cw>175] = NA
      odb$tripset = paste( odb$trip, odb$set_no, sep="~")
      odb$cpue.kg.trap = ( odb$totmass*1000)/odb$num_hook_haul
      save(odb, file=fn, compress=T)

      return( "complete" )
    }


  }


