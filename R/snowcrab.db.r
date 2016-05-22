
	snowcrab.db = function( DS, p=NULL, yrs=NULL) {
		# long!
		# handles all basic data tables, etc. ... can be broken into smaller pieces to make easier to maintain

		if (DS %in% c("set.odbc.redo", "set.odbc") ) {

	    fn.root =  file.path( project.datadirectory("snowcrab"), "data", "trawl", "SNCRABSETS" )
			dir.create( fn.root, recursive = TRUE, showWarnings = FALSE )

	    if (DS=="set.odbc") {
				out = NULL
				fl = list.files( path=fn.root, pattern="*.rdata", full.names=T ) 
        for ( fny in fl ) {
					load (fny)
					out = rbind( out, SNCRABSETS )
				}
				return (out)
			}

			require(RODBC)
			con = odbcConnect(oracle.snowcrab.server , uid=oracle.snowcrab.user, pwd=oracle.snowcrab.password, believeNRows=F)
			# believeNRows=F required for oracle db's

			for ( YR in yrs ) {
				fny = file.path( fn.root, paste( YR,"rdata", sep="."))
				SNCRABSETS = NULL
				SNCRABSETS = sqlQuery(con, 
									paste("select * from SNCRABSETS where EXTRACT(YEAR from BOARD_DATE) = ", YR) )
				save( SNCRABSETS, file=fny, compress=T)
				gc()  # garbage collection
				print(YR)
			}
			odbcClose(con)
			return (yrs)
		}



		if (DS %in% c("det.odbc.redo", "det.odbc") ) {
      fn.root =  file.path( project.datadirectory("snowcrab"), "data", "trawl", "SNCRABDETAILS" )
			dir.create( fn.root, recursive = TRUE, showWarnings = FALSE  )
  
	    if (DS=="det.odbc") {
				out = NULL
        fl = list.files( path=fn.root, pattern="*.rdata", full.names=T ) 
				for ( fny in fl ) {
					load (fny)
					out = rbind( out, SNCRABDETAILS )
				}
				return (out)
			}

			require(RODBC)
			con = odbcConnect(oracle.snowcrab.server , uid=oracle.snowcrab.user, pwd=oracle.snowcrab.password, believeNRows=F)
			# believeNRows=F required for oracle db's

			for ( YR in yrs ) {
				fny = file.path( fn.root, paste( YR,"rdata", sep="."))
				SNCRABDETAILS = NULL
				SNCRABDETAILS = sqlQuery(con, 
                  paste("select * from SNCRABDETAILS where EXTRACT(YEAR from BOARD_DATE) = ", YR) )
				save( SNCRABDETAILS, file=fny, compress=T)
				gc()  # garbage collection
				print(YR)
			}
			odbcClose(con)
              
      return (yrs)

		}




		if (DS %in% c("cat.odbc.redo", "cat.odbc") ) {

      fn.root =  file.path( project.datadirectory("snowcrab"), "data", "trawl", "SNTRAWLBYCATCH" )
			dir.create( fn.root, recursive = TRUE, showWarnings = FALSE  )
  
	    if (DS=="cat.odbc") {
				out = NULL
        fl = list.files( path=fn.root, pattern="*.rdata", full.names=T ) 
				for ( fny in fl ) {
					load (fny)
					out = rbind( out, SNTRAWLBYCATCH )
				}
				return (out)
			}

			require(RODBC)
			con = odbcConnect(oracle.snowcrab.server , uid=oracle.snowcrab.user, pwd=oracle.snowcrab.password, believeNRows=F)
			# believeNRows=F required for oracle db's

			for ( YR in yrs ) {
				fny = file.path( fn.root, paste( YR,"rdata", sep="."))
				SNTRAWLBYCATCH = NULL
				SNTRAWLBYCATCH = sqlQuery(con, 
                  paste("select * from SNTRAWLBYCATCH where EXTRACT(YEAR from BOARD_DATE) = ", YR) )
				save( SNTRAWLBYCATCH, file=fny, compress=T)
				gc()  # garbage collection
				print(YR)
			}
			odbcClose(con)
              
      return (yrs)
		}


		# ---------------------


    if (DS %in% c("setInitial.redo", "setInitial")) {
      fn = file.path( project.datadirectory( "snowcrab", "data" ), "set.initial.rdata" )

      if (DS =="setInitial") {
        set = NULL
				if (file.exists( fn ) ) load( fn) 
        return(set)
      }
        
      # field protocol:
      # Durometer measure if CW >= 60 mm  .. as in past
      # Chela height measure if CW >= 35 mm .. changed in 2007: previously it was 30 mm CW
      # Female abdomen width measure if CW >= 30 mm .. changed in 2007 .. previously all were measured.

      # data dump from the observer system
      # August 2015 added in setcd_id from observer system to address the MPA survey sets (type=11) and regular fix station sets (type=4)
      set = snowcrab.db( DS="set.odbc")
      names( set ) = rename.snowcrab.variables(names( set))
      setvars = c("trip", "set", "station", "stime", "observer", "cfa", "lon", "lat", "lon1", "lat1", "towquality", "Zx", "Tx", "gear", "sa", "dist", "dist0" )
      print('need to addin the mpa station index')

	    set$trip = as.character(set$trip)
      set$set  = as.integer(set$set)

      set = set[ order( set$sdate ), ]

      # first estimate of distance of trawl track based upon logged start and end locations
      set$dist0 = geosphere::distGeo( set[,c("lon", "lat" )], set[,c( "lon1", "lat1") ] )  # m
      ii = which( set$dist0 > 1000 | set$dist0 < 100 )
      if ( length(ii) > 0 ) {
        print( "Positional information incorrect or suspect in the following:" )
        oo = set[ ii, c("trip", "set", "dist", "dist0", "lon", "lat", "lon1", "lat1", "towquality") ]
        print( oo ) 
      }
      set$dist0[ii] = NA

      # vs set$dist which are hand computed by Moncton for the early time series
      #  plot( dist0 ~ dist, set) shows good correspondence

      ############
      # We had discussed that one station got missed in the assessment last year as there was a 
      # typo in the Haul Code ID which indicated that the tow was bad when it was actually a good tow. 
      # You can easily build a double check into your script to catch that. When the tow is "bad", there
      # is an estimated catch of "0". In all cases where the tow is good, the estimated catch is >0. 
      # You can add a line that says if the estimated catch is >0 but the HAUL_CD_ID =4 (bad tow),
      # flag that entry for us to check before proceeding.
      oo = which( set$towquality ==4 & set$totmass >0 )
      if (length(oo)>1) {
        print( "The following needs some checking as there is a mass estimate but towquality set to be bad" )
        print( set[ oo, ] )
      }
  
      dbug.2011 = T
      if ( dbug.2011 ) {
        # one-off error corrections until database gets refreshed
        i = which(set$trip=="S15092006" & set$set==3)
        if (length(i)==1) set$totmass[i] = 0

        i = which(set$trip=="S21092007" & set$set==12)
        if (length(i)==1) set$towquality[i] = 1
      }
      dbug.2009 = T 
      if(dbug.2009) {
        i = which(set$trip =='S05092009' & set$set == 9 & set$station==189)
        set$stime[i] <- 1532
      }

      dbug.2014 = T
      if(dbug.2014) {
        i = which(set$trip == 'S26112014' &set$set==7)
        set$station[i] = 212
        i = which(set$trip == 'S07102014' &set$set==13)
        set$stime[i] = 2146
        i = which(set$trip == 'S06112014' &set$set==8)
        set$lon[i] = 60.08317
      }

      debug.2015= T
      if(debug.2015) {
        i= which(set$trip_id == '100045549' & set$set==1 & set$station ==1)
        set$stime[i] = 0449
        i= which(set$trip_id == '100045553' & set$set==4 & set$station ==615)
        set$stime[i] = 0944
        i= which(set$trip_id == '100045757' & set$set==4 & set$station ==210)
        set$stime[i] = 0845
        i= which(set$trip_id == '100045764' & set$set ==5 & set$station == 344)
        set$stime[i] = 0947
	
	      i= which(set$trip == 'S21112015' & set$set ==6 & set$station == 378)
        set$stime[i] = 1102
	      i= which(set$trip == 'S26102015' & set$set ==5 & set$station == 724)
        set$stime[i] = 1948

	      i=which(set$trip_id == '100045765' & set$set==15 & set$station==911)
        set$sdate[i]= '2015-12-12 01:00:00'

        i=which(set$trip_id == '100045765' & set$set ==15 & set$station==809)
        set$sdate[i] = '2015-12-12 01:00:00'

        i=which(set$trip== 'S12092015' & set$set==7)
        set$towquality[i] = 4

        i=which(set$trip== 'S05092015' & set$set==9 & set$station==187)
        set$lon[i]= 60.2575
        print("2015 debugging complete")
      }

       
      set = set[,setvars]
      set$sa[ which(set$sa==0) ] = NA
      set$sa = set$sa / 10^6    # convert to km2 ... sa was stored as m2

      set = set[ which(set$towquality==1) , ]
      nset0 = nrow(set)
      
      set$station = as.numeric(set$station)

      set$stime = ifelse( nchar(set$stime)==4,
        paste(substring(set$stime, 1,2),    ":", substring(set$stime, 3,4), ":", "00", sep=""), set$stime )
      set$stime = ifelse( nchar(set$stime)==3,
        paste("0", substring(set$stime, 1,1), ":", substring(set$stime, 2,4), ":", "00", sep=""), set$stime)

      set$stime = ifelse( nchar(set$stime)==2, NA,  set$stime )  # these are indeterminate

      # ---------------------
      # local lookuptable for netmind/minilog data (historical data only)
      sntows = read.table(file.path( project.datadirectory("snowcrab"), "data", "lookuptables", "sntow.csv"),  sep=";", as.is=T, colClasses="character", header=T)
      sntow.vars = c("trip",  "set", "patchtime", "netmind", "minilog")

      set = merge(set, sntows[,sntow.vars], by=c("trip", "set"), all.x=T, all.y=F, sort=F)
      if ( nrow( set) != nset0 ) { print("Merge error"); stop() }
      
      set$lon = -set$lon # the data are stored as degress West convert to standard degrees E notation
      set$lon1 = -set$lon1 # the data are stored as degress West convert to standard degrees E notation
     
      overwrite = which( is.na(set$stime))
      set$stime[overwrite] = set$patchtime[overwrite]
      set$patchtime = NULL
      
      i.na = which(is.na(set$stime))
      if (length( i.na) > 0 ) {
        set$stime[i.na] = "12:00:00"  # set to noon by default if there are any more na's due to the merge
      }
      
      # "chron" is the best estimate of sampling time 
      # sdate (POSIXct, America/Halifax AST/ADT) does not seem to be reliable 
      # and so we use minilog data where possible in the recent period
      # and then records from the stime and trip id where minilog data are absent 
      set$chron = tripcode.to.chron( set$trip, set$stime )  # using chron .. deprecated 
      set$timestamp = tripcode.to.timestamp( set$trip, set$stime, tzone="America/Halifax" )  # using lubridate/POSIXct

      set$stime = NULL ### --need to check timezone!!! TODO .... seems to be "America/Halifax" .. compare with seabird/minilog
      
      i = which(is.na(set$chron))
      if (length(i)>0) set$chron[i] = tripcode.to.chron( set$trip[i], "12:00:00" )
      
      set$julian =  convert.datecodes(set$chron, "julian")
      set$yr = convert.datecodes(set$chron, "year")
     
      save( set, file=fn, compress=TRUE )

      ## michelle:: please do not place hard-links into the code as this will force a fail for others ..
      ## this is probably better created as a functions and you can send to the data for a save into OGR format
      ## to a location of your choice? ..
      michelle = FALSE
      if (michelle) {
        #MG Save the trawl file to a shapefile to bring into ArcGIS
        shape.set <- set
        shape.set$lon <- -shape.set$lon
        shape.set$chron <- as.character(shape.set$chron)

        set.cords <- shape.set[, c("lon", "lat")]
        sdf.set <- SpatialPointsDataFrame(set.cords, data=shape.set)
        proj4string(sdf.set) <- CRS(p$geog.proj)
        shpdir = file.path(project.datadirectory("snowcrab"), "maps", "shapefiles", "survey")
        setwd(shpdir)
        
        writeOGR(sdf.set, ".", "SurveyDataUpdate", driver="ESRI Shapefile", overwrite=T)
        setwd("/home/michelle/tmp")  
        shp.path <- paste("SurveyDataUpdate shapefile created at", shpdir, sep=" ")
        print(shp.path)
      }

      return ( fn )
    }
     

    # --------------------

    if (DS %in% c("det.initial", "det.initial.redo") ) {
      fn = file.path(project.datadirectory("snowcrab"), "R", "det.initial.rdata")
      if (DS =="det.initial") {
        load(fn)
        return(det)
      }
      
      X = snowcrab.db( DS="set.clean" )
      det = snowcrab.db( DS="det.odbc"  )

      names( det ) = rename.snowcrab.variables(names(det) )
      detvars = c( "trip", "set", "crabno", "sex", "cw", "mass", "abdomen", "chela", "mat",
                   "shell", "gonad", "eggcol", "eggPr", "durometer", "legs")
      
      # merge in the sa which is used as a weighting factor of most analyses
      det = merge(x=det[,detvars], y=X[,c("trip","set","sa")], by=c("trip", "set"), all.x=T, all.y=F)
     
      # Trips and sets with no matching SA ...  
      ii = which( !is.finite(det$sa) )
      #print ("DET without matching SA ... due to bad tows? ... check these ") 
      #print (det[ii, c("trip","set")])
      #print (det[ii,])
      det$sa[ii] = median(det$sa, na.rm=T )

      i.mat = which(det$mat==1)
      i.imm = which(det$mat==2) 
      i.other = setdiff( 1:nrow( det), c(i.mat, i.imm) ) 
      
      det$mat[ i.mat ] = mature
      det$mat[ i.imm ] = immature
      det$mat[ i.other ] = mat.unknown
      
      i.male = which( det$sex == 1 )
      i.female = which( det$sex == 2 )
      i.other =  setdiff( 1:nrow(det), c(i.male, i.female) )

      det$sex [ i.male ] = male  # male defined as a gloabl parameter
      det$sex [ i.female ] = female  # female defined as a gloabl parameter
      det$sex [ i.other ] = sex.unknown  # sex codes defined as a gloabl parameter
      

      #Identify morphology errors and print, save to CSV
      yr.e <- p$current.assessment.year
      fn.e = file.path(project.datadirectory("snowcrab"), "data", "trawl", "morphology.errors")
      dir.create(fn.e, recursive=T, showWarnings=F)
      outfile.e =  file.path( fn.e, paste("morphologyerrors", yr.e, ".csv", sep=""))
      
      #Sex.e: Unknown Sex
      sex.e <- det[which(det$sex==0),]
      sex.e.2015 <- sex.e[grep("2015", sex.e$trip),]
      sex.e$error <- 'sex.e'
      #Cw.e: Carapace Width below 5 or greater than 185
      cw.e <- det[ which(det$cw<5 | det$cw>185 ),]
      cw.e$error <- 'cw.e'
      #Chela.e: Chela less than 1 or greater than 50
      chela.e <- det[which(det$chela < 1 | det$chela > 50  ),]
      chela.e$error <- 'chela.e'
      #Abdomen.e:Abdomen less than 1 and greater than 66
      abdomen.e <- det[which(det$abdomen < 1 | det$abdomen > 66 ),]
      abdomen.e$error <- 'abdomen.e'
      #Mass.e: Mass less than 1 or greater than 1500
      mass.e <- det[which( det$mass < 1 | det$mass > 1500  ),]
      mass.e$error <- 'mass.e'
      #Sex.a: Indeterminate sex based on measurements taken (abdomen values where sex=male)
      sex.a <- det[which(is.finite( det$abdomen ) & det$sex==1),]
      sex.a$error <- 'sex.a'      
      #Sex.c: Indeterminate sex based on measurements taken (chela values where sex=female
      sex.c <- det[which(is.finite( det$chela ) & det$sex==2),]
      sex.c$error <- 'sex.c'
      
     
      det$cw [ which(det$cw<5 | det$cw>185 ) ] = NA  # a few zero-values
      det$chela [ which(det$chela < 1 | det$chela > 50  )] = NA # remove 0's and unreliably small values
      det$abdomen [ which(det$abdomen < 1 | det$abdomen > 66 ) ] = NA # remove 0's and unreliably small values
      det$mass  [ which( det$mass < 1 | det$mass > 1500  )]= NA # remove 0's and unreliably small /large values
   
      # indeterminate sex based upon measurements taken
      iii = which( is.finite( det$abdomen ) & det$sex==male )
      det$sex[iii] = sex.unknown

      iii = which( is.finite( det$chela ) & det$sex==female )
      det$sex[iii] = sex.unknown

      # assume a reading error of +/- 0.25 mm and +/- 0.25 g
      # changes in reading resolution occurs over time
      # .. this helps to smooth the data
      
      # det$cw = jitter(det$cw, amount=0.2)
      # det$chela = jitter(det$chela, amount=0.2)
      # det$abdomen = jitter(det$abdomen, amount=0.2)
      # det$mass =  jitter(det$mass, amount=0.2)  # mass in grams


      det = predictweights (det )
     
      unreliable = which( det$mass < 0.25 | det$mass > 1800  )
      det$mass  [ unreliable ]= NA # remove 0's and unreliably small /large values
      det$cw  [ unreliable ]= NA # remove as these cw were used to create the above unreliable masses
      
      det = predictmaturity (det, method="logistic.regression")
      
      #Mat.e: Unknown Maturity
      mat.e <- det[which(det$mat ==0),]
      mat.e$error <- 'mat.e'
      
      primiparous = filter.class( det, "primiparous")
      multiparous = filter.class( det, "multiparous")

      det$fecundity = NA
      det$fecundity[primiparous] = fecundity.allometric( cw=det$cw[primiparous], method="moncton.primiparous" )
      det$fecundity[multiparous] = fecundity.allometric( cw=det$cw[multiparous], method="moncton.multiparous" )
      det$fecundity[ which(det$fecundity> 250000) ] = NA
      
      save(det, file=fn, compress=T)

      # do only after the above save 
      
      allometry.snowcrab( "cw.mass", "male", redo=T )
      allometry.snowcrab( "chela.mass", "male", redo=T  )
      allometry.snowcrab( "cw.chela.mat", "male", redo=T  )
      allometry.snowcrab( "cw.mass", "female", redo=T  )
      allometry.snowcrab( "chela.mass", "female", redo=T  )
      allometry.snowcrab( "cw.chela.mat", "female", redo=T  )
      
      names.e <- list(mat.e, sex.e, cw.e, chela.e, abdomen.e, mass.e, sex.a, sex.c)
      errors = NULL
      for (e in names.e){
        if (nrow(e) > 0)
          errors <- rbind(errors, e)
      }
      
      errors.yearly <- errors[grep(yr.e, errors$trip),]
      
      write.csv(errors.yearly, file=outfile.e)
      print("Morphology Errors saved to file")
      print(outfile.e)
      cat("ERROR CODES\
      Mat.e: Unknown Maturity\
      Sex.e: Unknown Sex\
      Cw.e: Carapace Width below 5 or greater than 185\
      Chela.e: Chela less than 1 or greater than 50\
      Abdomen.e:Abdomen less than 1 and greater than 66\
      Mass.e: Mass less than 1 or greater than 1500\
      Sex.a: Indeterminate sex based on measurements taken (abdomen values where sex=male)\
      Sex.c: Indeterminate sex based on measurements taken (chela values where sex=female\n")
      

      return ( "Complete" )
    }

    # ------------------------------


    if (DS %in% c("cat.initial", "cat.initial.redo") ) {
      fn = file.path(project.datadirectory("snowcrab"), "R", "cat.initial.rdata")
      if(DS =="cat.initial" ) {
        load(fn)
        return(cat)
      }     

			# two steps:: bycatch from the cat tables (missing snow crab) 
      # and determine totals from the snow crab det tables
      
      det = snowcrab.db( DS="det.initial" )
      
      cat = snowcrab.db( DS="cat.odbc" )
      names( cat ) = rename.snowcrab.variables(names( cat ) )
      
      # two different time periods (pre and post Moncton)
      # the earlier was saved as totmass.kept and the latter as discarded 
      cat$totmass = cat$totmass.discarded  # bycatch weights are stored here
			cat$totmass.discarded = NULL # clean up 

			## note: historical data prior to 2005 did not capture total weights only numbers, occasionally
			##
			gc()
			
      ii = which( cat$totmass <= 0.00001 )
      cat$totmass[ ii] = 0  # observer database does not allow storage of zero values
      
      catvars =  c("trip", "set", "spec", "totno", "totmass")

			# clean species codes ... this causes multiple entries for some species that need to be summed up
      # cat$spec = taxonomy.parsimonious( spec=cat$spec )
      # --- no longer here ... only when integrated into bio.db

      # remove data where species codes are ambiguous, or missing or non-living items
      xx = which( !is.finite( cat$spec) ) 
      if (length(xx)>0) cat = cat[ -xx, ] 
      cat = cat[ taxonomy.filter.taxa( cat$spec, taxafilter="living.only", outtype="groundfishcodes" ) , ]


      # update catch biomass/numbers due to altering of species id's
				cat = cat[,catvars]
				catn = as.data.frame( xtabs( cat$totno ~ as.factor(trip) + as.factor(set) + as.factor(spec), data=cat ) )
				catb = as.data.frame( xtabs( cat$totmass ~ as.factor(trip) + as.factor(set) + as.factor(spec), data=cat ) )
			
				names(catn) = c("trip", "set", "spec", "totno")
				names(catb) = c("trip", "set", "spec", "totmass")

				chars = "trip"
				numbers = c("set", "spec")
      
				for (j in chars) catn[,j] = as.character(catn[,j] )
				for (j in numbers) catn[,j] = as.numeric(as.character(catn[,j] ))
      
				for (j in chars) catb[,j] = as.character(catb[,j] )
				for (j in numbers) catb[,j] = as.numeric(as.character(catb[,j] ))
    

				catn = catn[ which(catn$totno > 0) , ] 
				catb = catb[ which(catb$totmass > 0) , ] 

        # this contains all the data from the by-catch tables
        x = merge( catn, catb, by=c("trip", "set", "spec"), all.x=T, all.y=T, sort=F )
          oo = which( !is.finite(x$totno) & !is.finite(x$totmass) )
          if (length(oo)>0) x = x[-oo,]

      # compute snow crab abundance from det tables
      numbers = as.data.frame( xtabs( ~ as.factor(trip) + as.factor(set), data=det ) )
      names(numbers) = c("trip", "set", "totno")
			numbers$trip = as.character( numbers$trip )
			numbers$set = as.numeric( as.character(numbers$set ))

      good = which(is.finite(det$mass))
      biomass = as.data.frame(xtabs( mass ~ as.factor(trip) + as.factor(set), data=det[good,], exclude="" ) )
      names(biomass) = c("trip", "set", "totmass")
      biomass$trip = as.character( biomass$trip )
			biomass$set = as.numeric( as.character(biomass$set ))
			biomass$totmass = biomass$totmass / 1000  # convert from grams to kg
      # !!!!!!!! must verify units of other species from observer system

      snowcrab = merge(x=numbers, y=biomass, by=c("trip", "set"), all=T)
      snowcrab = snowcrab[ which( as.character(snowcrab$trip) != "-1") , ]
			
       snowcrab$spec = taxonomy.recode(to='parsimonious',from='spec', tolookup = 2526 )  # 2526 is the code used in the groundfish/snow crab surveys .. convert to internally consistent state
      # longer here -- in bio db only

			final = snowcrab[,names(x)]  # get the right sequence of variables
			
				# strip zeros when both  no and mass are 0
				z = which( final$totno==0 & final$totmass==0)
				if (length(z) > 0 ) final = final[-z, ]

				# data for which mass estimates were not recorded -- mostly crab not assigned a fishno --> fragments of crab that were not measureable .. asign a small weight
				o = which(final$totmass==0 & final$totno == 1 ) 
				if (length(o) >0 ) final$totmass[o] = median( det$mass / 1000, na.rm=T )  # kg 
		
				# catch any strange data
				o = which(final$totmass==0 & final$totno > 0 ) 
				if (length(o) >0 ) final$totmass[o] = min( final$totmass[ which(final$totmass>0) ] )  # kg 
		
      # -----------------------
      # merge in the snowcrab weights
      cat = rbind(x, final)



			# estimate number from size and weight
	

			# fix missing numbers and mass estimates:
			# zeros for one while nonzeros for correpsonding records
			meanwgt = cat$totmass / cat$totno
			good = which( is.finite( meanwgt) & is.finite( 1/meanwgt ) )
			mw = as.data.frame( unlist( (tapply( log( meanwgt[good]), cat$spec[good], mean )) ))
			names(mw) = "meanweight"
			mw$spec= as.numeric( as.character( rownames(mw) ) )
			mw$meanweight = exp( mw$meanweight )
			mw = mw[which(is.finite(mw$meanweight)) ,]
		
			# add groundfish data if it does not exist already
			loadfunctions( "groundfish" )
			gs = groundfish.db( "cat" )  
			meanwgt = gs$totwgt / gs$totno
			good = which( is.finite( meanwgt) & is.finite( 1/meanwgt ) )
			mw2 = as.data.frame( unlist( (tapply( log( meanwgt[good]), gs$spec[good], mean )) ))
			names(mw2) = "meanweight"
			mw2$spec= as.numeric( as.character( rownames(mw2) ) )
			mw2$meanweight = exp( mw2$meanweight )
			mw2 = mw2[which(is.finite(mw2$meanweight)) ,]
			
			mw = merge(mw, mw2, by="spec", all=T, suffixes=c("","gs") )
			rm(gs, mw2, meanwgt, good); gc()

			i = which( !is.finite( mw$meanweight) & is.finite(mw$meanweightgs) )
			mw$meanweight[i] = mw$meanweightgs[i]



      ii = which( is.na(cat$totno) & cat$totmass >  0 ) 
      if (length(ii)>0) {
        # replace each number estimate with a best guess based upon average body weight in the historical record
        uu = unique( cat$spec[ii] )
        for (u in uu ) {
          os =  which( mw$spec==u ) 
          if (length( os)==0 ) next()
          toreplace = intersect( ii, which( cat$spec==u) )
          cat$totno[toreplace] = cat$totmass[toreplace] / mw$meanweight[os]
        }
      }

      jj = which( cat$totno >  0 & is.na(cat$totmass) ) 
      if (length(jj)>0) {
        # replace each number estimate with a best guess based upon average body weight in the historical record
        uu = unique( cat$spec[jj] )
        for (u in uu ) {
          os =  which( mw$spec==u ) 
          if (length( os)==0 ) next()
          toreplace = intersect( jj, which( cat$spec==u) )
          cat$totmass[toreplace] = cat$totno[toreplace] * mw$meanweight[os]
        }
      }

      save( cat, file=fn, compress=T)
      return("Complete")
    }


    # -------------

    if ( DS %in% c("set.merge.det","set.merge.det.redo") ) {
    # browser()

      fn = file.path( project.datadirectory("snowcrab"), "R", "set.biologicals.rdata")

      if (DS=="set.merge.det" ) {
        load(fn) 
        return( set)
      }

      factors = c("trip", "set")

      X = snowcrab.db( DS="set.clean" )
      #X2015 = X[which(X$yr==2015),]
      #head(X2015)
      Y = snowcrab.db( DS="det.initial" )
 
      # add various variables to set-level data

      # add fecunity estimates
        fecund = as.data.frame.table( 
          tapply(Y$fecundity, INDEX=Y[,factors], FUN=sum, na.rm=T, simplify=T ) 
        )
        names(fecund) = c(factors, "fecundity")
        fecund = factor2character(fecund, factors)
        X = merge(x=X, y=fecund, by=factors, all.x=T, sort=F )
        X$fecundity = X$fecundity / X$sa / 10^6   # rescale due to large numbers

      # add sex ratios of all crabs
        y=sex.ratios(Y[,c(factors, "sex")], factors)
        names(y) = c(factors, "no.male.all", "no.female.all", "sexratio.all")
        X = merge(x=X, y=y, by=factors, all.x=T )

      # add sex ratios of all mature crabs
        y = sex.ratios(Y[filter.class(Y, "mat"), c(factors, "sex")], factors)
        names(y) = c(factors, "no.male.mat", "no.female.mat", "sexratio.mat")
        X = merge(x=X, y=y, by=factors, all.x=T )

      # add sex ratios of all immature crabs
        y = sex.ratios(Y[filter.class(Y, "imm"), c(factors, "sex")], factors)
        names(y) = c(factors, "no.male.imm", "no.female.imm", "sexratio.imm")
        X = merge(x=X, y=y, by=factors, all.x=T )

      # ------------------------------------------------------------------------------------------------
      # add (mean,var,count) of cw
      # all snowcrabs
      # y = bodysize(Y[,c(factors, "mass")], factors, "mass", logtransform=T)  # <<< example for extracting mean mass
        y = bodysize(Y[,c(factors, "cw")], factors, "cw", logtransform=T)
        names(y) = c(factors, "cw.mean", "cw.var", "cw.n")
        X = merge(x=X, y=y, by=factors, all.x=T )
        X$cw.n = X$cw.n / X$sa

      # commercially sized male snowcrabs
        y = bodysize(Y[filter.class(Y, "m.com"), c(factors, "cw")], factors, "cw", logtransform=T)
        names(y) = c(factors, "cw.comm.mean", "cw.comm.var", "cw.comm.n")
        X = merge(x=X, y=y,  by=factors, all.x=T )
        X$cw.comm.n = X$cw.comm.n / X$sa

      # noncommercial male snowcrabs
        y = bodysize(Y[filter.class(Y, "m.ncom"), c(factors, "cw")], factors, "cw", logtransform=T)
        names(y) = c(factors, "cw.notcomm.mean", "cw.notcomm.var", "cw.notcomm.n")
        X = merge(x=X, y=y,  by=factors, all.x=T )
        X$cw.notcomm.n = X$cw.notcomm.n / X$sa

      # mature female snowcrabs
        y = bodysize(Y[filter.class(Y, "f.mat"), c(factors, "cw")], factors, "cw", logtransform=T)
        names(y) = c(factors, "cw.fem.mat.mean", "cw.fem.mat.var", "cw.fem.mat.n")
        X = merge(x=X, y=y,  by=factors, all.x=T )
        X$cw.fem.mat.n = X$cw.fem.mat.n / X$sa

      # immature female snowcrabs
        y = bodysize(Y[filter.class(Y, "f.imm"), c(factors, "cw")], factors, "cw", logtransform=T)
        names(y) = c(factors, "cw.fem.imm.mean", "cw.fem.imm.var", "cw.fem.imm.n")
        X = merge(x=X, y=y,  by=factors, all.x=T )
        X$cw.fem.imm.n = X$cw.fem.imm.n / X$sa


      # mature male snowcrabs
        y = bodysize(Y[filter.class(Y, "m.mat"), c(factors, "cw")], factors, "cw", logtransform=T)
        names(y) = c(factors, "cw.male.mat.mean", "cw.male.mat.var", "cw.male.mat.n")
        X = merge(x=X, y=y,  by=factors, all.x=T )
        X$cw.male.mat.n = X$cw.male.mat.n / X$sa

      # immature male snowcrabs
        y = bodysize(Y[filter.class(Y, "m.imm"), c(factors, "cw")], factors, "cw", logtransform=T)
        names(y) = c(factors, "cw.male.imm.mean", "cw.male.imm.var", "cw.male.imm.n")
        X = merge(x=X, y=y,  by=factors, all.x=T )
        X$cw.male.imm.n = X$cw.male.imm.n / X$sa

    # ------------------------------------------------------------------------------------------------
    # add biomass of various components of the snowcrab population
    #      X = setmerge(X, det, varname="totmass.all", filter="all", variable="mass")
    #      ... better to use the total catch tables as subsampling may be used in the future

      print( "Biomass density estimates complete" )
      
      vars = lookup.biomass.vars()
      for (i in 1:nrow(vars)) {
        print(vars[i,])
        X=setmerge(X, Y, varname=vars[i,1], filter=vars[i,2], variable="mass")
        X[, vars[i,1] ] = X[, vars[i,1] ] / 10^6 # grams .. convert to metric tons
      }

    # ------------------------------------------------------------------------------------------------
    # add numbers of various components of the snowcrab population
    #      X = setmerge(X, Y, varname="totno.all", filter="all", variable="number")
    #       ... better to use the total catch tables as subsampling may be used in the future
      
      vars = lookup.numbers.vars()

      for (i in 1:nrow(vars)) {
        print(vars[i,])
        X=setmerge(X, Y, varname=vars[i,1], filter=vars[i,2], variable="number")
      }

      print( "Numerical density estimates complete" )

    # ------------------------------------------------------------------------------------------------
    # add biomass and numbers directly from the catch (cat) tables (e.g. for multi-species analysis)
    # using a separate system of analysis and access is probably better

      rm(Y); gc()
      set = X

      if ( nrow( snowcrab.db( DS="set.clean" )) != nrow( set) ) {   print( "Merge failure ... " );  stop()    }

      save( set, file=fn, compress=T )

    return ( "Complete" )
  }


    # ---------------------------

    if ( DS %in% c("set.merge.cat","set.merge.cat.redo") ) {

      fn = file.path( project.datadirectory("snowcrab"), "R", "set.cat.rdata")

      if (DS %in% c("set.merge.cat" )) {
        load(fn) 
        return( set )
      }

      factors = c("trip", "set")

      X = snowcrab.db( DS="set.merge.det" )
      X2015 = X[which(X$yr == 2015),]
      print(head(X2015))
      
      cat = snowcrab.db( DS="cat.initial" )
      cat2015 = cat[grep("2015", cat$trip),]
      print(head(cat2015))
      
      cat0 = cat[ taxonomy.filter.taxa( cat$spec, taxafilter="snowcrab", outtype="groundfishcodes"), c(factors, "totno")]
      names(cat0) = c(factors, "totno.all")
      X = merge(x=X, y=cat0, by=factors, all.x=T )
      X$totno.all   = X$totno.all   / X$sa
      X$totno.all[!is.finite(X$totno.all)] = 0  # convert na's to zero

      cat0 = cat[ taxonomy.filter.taxa( cat$spec, taxafilter="snowcrab", outtype="groundfishcodes" ), c(factors, "totmass")]
      names(cat0) = c(factors, "totmass.all")
      X = merge(x=X, y=cat0, by=factors, all.x=T )
      X$totmass.all = X$totmass.all / X$sa
      X$totmass.all[!is.finite(X$totmass.all)] = 0  # convert na's to zero
      # the above masses are in kilograms .. convert to metric tons
      var="totmass.all";  X[,var] = X[,var] / 10^3

      set = X
      
      if ( nrow( snowcrab.db( DS="setInitial" )) != nrow( set) ) {   print( "Merge failure ... " );  stop()    }
      X2015 = X[which(X$yr == 2015),]
      print(head(X2015))
      
      save( set, file=fn, compress=T )
      return ( "Complete" )
    }



    # --------------------------------
    
  
    if (DS %in% c("set.minilog.seabird.retired", "set.minilog.seabird.redo.retired")) {
      # merge setInitial with minilog stats, seabird stats to generate sensible start and end times 
      # used for creating netmind stats /metrics
      fn = file.path( project.datadirectory( "snowcrab", "data"), "set.minilog.seabird.rdata" )
    
      if (DS=="set.minilog.seabird.redo.retired") {
        set = NULL
        if ( file.exists( fn) ) load (fn)
        return (set)
      }

      if (DS=="set.minilog.seabird.retired") {
        tzone = "America/Halifax"
        set = snowcrab.db( DS="setInitial") 
        set.names= names(set)
        #set2015 = set[which(set$yr==2015),]
        
        #These have 2015 data in them
        sb = seabird.db( DS="set.seabird.lookuptable" )
        ml = minilog.db( DS="set.minilog.lookuptable" )
        # ml2015 = ml[grep("2015", ml$trip),]
        nm = netmind.db( DS="set.netmind.lookuptable" )

        set = merge( set, sb, by=c("trip","set"), all.x=T, all.y=F, sort=F, suffixes=c("", ".seabird") )
        set = merge( set, ml, by=c("trip","set"), all.x=T, all.y=F, sort=F, suffixes=c("", ".minilog") )
        set = merge( set, nm, by=c("trip","set"), all.x=T, all.y=F, sort=F, suffixes=c("", ".netmind") )
        
        set = set[ , c( set.names, "seabird_uid", "minilog_uid", "netmind_uid" ) ]

        #These do not have 2015 data in them 
        sbStats =  seabird.db( DS="stats" )
   #     sbStats = sbStats[ , c("seabird_uid", "z", "zsd", "t", "tsd", "n", "t0", "t1", "dt" ) ]
        sbStats = sbStats[ , c("seabird_uid",'trip','set', "z", "zsd", "t", "tsd", "n", "t0", "t1", "dt" ) ]
        #sbStats2015= sbStats[grep("2015", sbStats$trip),]
        mlStats =  minilog.db( DS="stats" )
   #     mlStats = mlStats[ , c("minilog_uid", "z", "zsd", "t", "tsd", "n", "t0", "t1", "dt" ) ]
        mlStats = mlStats[ , c("minilog_uid",'trip','set', "z", "zsd", "t", "tsd", "n", "t0", "t1", "dt" ) ]
         #mlStats2015= mlStats[grep("2015", mlStats$trip),]



        names( mlStats ) = c("minilog_uid",'trip','set', "z.ml", "zsd.ml", "t.ml", "tsd.ml", "n.ml", "t0.ml", "t1.ml", "dt.ml" )

#      set = merge( set, sbStats, by="seabird_uid", all.x=TRUE, all.y=FALSE, sort=FALSE )
#      set = merge( set, mlStats, by="minilog_uid", all.x=TRUE, all.y=FALSE, sort=FALSE )
        set = merge( set, sbStats, by=c("trip","set"), all.x=TRUE, all.y=FALSE, sort=FALSE )
        set = merge( set, mlStats, by=c("trip","set"), all.x=TRUE, all.y=FALSE, sort=FALSE )
        set$t0 = as.POSIXct(set$t0,format="%Y-%m-%d %H:%M:%S", tz=tzone,origin=lubridate::origin )
        set$t1 = as.POSIXct(set$t1,format="%Y-%m-%d %H:%M:%S", tz=tzone ,origin=lubridate::origin)
        set = toNums(set,c('dt','t0.ml','t1.ml', 'dt.ml'))

        # use seabird data as the standard, replace with minilog data where missing
        ii = which(!is.finite( set$t0) )
        if (length(ii) > 0 )  set$t0[ ii] = as.POSIXct(set$t0.ml[ii],origin=lubridate::origin, tz=tzone)
   
        ii = which(!is.finite( set$t1) )
        if (length(ii) > 0 )  set$t1[ ii] = as.POSIXct(set$t1.ml[ii],origin=lubridate::origin, tz=tzone)
        
        ii = which(!is.finite( set$z) )
        if (length(ii) > 0 )  set$z[ ii] = set$z.ml[ii]
      
        ii = which(!is.finite( set$zsd) )
        if (length(ii) > 0 )  set$zsd[ ii] = set$zsd.ml[ii]
         
        ii = which(!is.finite( set$t) )
        if (length(ii) > 0 )  set$t[ ii] = set$t.ml[ii]
      
        ii = which(!is.finite( set$tsd) )
        if (length(ii) > 0 )  set$tsd[ ii] = set$tsd.ml[ii]
        


        ii = which(!is.finite( set$dt) )
        if (length(ii) > 0 )  set$dt[ ii] = set$dt.ml[ii]
      
        set = set[ ,c(set.names, "netmind_uid", "z", "zsd", "t", "tsd", "t0", "t1", "dt" ) ]
     #   set2015 = set[which(set$yr==2015),]  ## why are we doing this over here? and if you want to do this make it more generic rather than for a specific year? (Jae)
     #   print(head(set2015))
        save( set, file=fn, compress=TRUE)
        return (fn)
      }
    }

    # --------------------------------
 
  

    if ( DS %in% c("set.clean", "set.clean.redo") ) {
    
      # merge seabird, minilog and netmind data and do some checks and cleaning
      fn = file.path( project.datadirectory( "snowcrab" ), "data", "set.clean.rdata" )

      if ( DS=="set.clean" ) {
        set= NULL
        if (file.exists( fn) ) load( fn )
        return (set) 
      }
      
      # the beginning here is identical to the netmind.db( "stat.redo" ) .. simpler to keep it this way (jae)
      tzone = "America/Halifax"
      set = snowcrab.db( DS="setInitial") 

      sbStats =  seabird.db( DS="stats" )
      sbv = c('trip','set', "z", "zsd", "t", "tsd", "n", "t0", "t1", "dt" )
      set_sb = merge( set[, c("trip", "set") ], sbStats[,sbv], by=c("trip","set"), all.x=TRUE, all.y=FALSE, sort=FALSE )
      # tapply( as.numeric(set_sb$dt), year(set_sb$t1), mean, na.rm=T )
      # tapply( as.numeric(set_sb$dt), year(set_sb$t1), function(x) length(which(is.finite(x))) )

      mlStats =  minilog.db( DS="stats" )
       # mlStats$dt = as.numeric(mlStats$dt )
      mlv =  c('trip','set', "z",    "zsd",    "t",    "tsd",    "n",    "t0",    "t1",    "dt" ) 
      set_ml = merge( set[, c("trip", "set") ], mlStats[,mlv], by=c("trip","set"), all.x=TRUE, all.y=FALSE, sort=FALSE )
      # tapply( as.numeric(set_ml$dt), lubridate::year(set_ml$t1), mean, na.rm=T )
      # tapply( as.numeric(set_ml$dt), year(set_ml$t1), function(x) length(which(is.finite(x))) )

      set = merge( set, set_sb, by=c("trip", "set" ), all.x=TRUE, all.y=FALSE, sort=FALSE )
      set = merge( set, set_ml, by=c("trip", "set" ), all.x=TRUE, all.y=FALSE, sort=FALSE, suffixes=c("", ".ml" ))

      # use seabird data as the standard, replace with minilog data where missing
      ii = which(!is.finite( set$t0) )
      if (length(ii) > 0 )  set$t0[ ii] = set$t0.ml[ii]

      ii = which(!is.finite( set$t1) )
      if (length(ii) > 0 )  set$t1[ ii] = set$t1.ml[ii]
      
      ii = which(!is.finite( set$z) )
      if (length(ii) > 0 )  set$z[ ii] = set$z.ml[ii]
    
      ii = which(!is.finite( set$zsd) )
      if (length(ii) > 0 )  set$zsd[ ii] = set$zsd.ml[ii]
       
      ii = which(!is.finite( set$t) )
      if (length(ii) > 0 )  set$t[ ii] = set$t.ml[ii]
    
      ii = which(!is.finite( set$tsd) )
      if (length(ii) > 0 )  set$tsd[ ii] = set$tsd.ml[ii]

      ii = which(!is.finite( set$dt) )
      if (length(ii) > 0 )  set$dt[ ii] = set$dt.ml[ii]
     
      tokeep = grep( "\\.ml$", colnames(set), invert=TRUE )
      set = set[, tokeep]
      set$n = NULL
      # tapply( as.numeric(set$dt), year(set$t1), mean, na.rm=T )
      # tapply( as.numeric(set$dt), year(set$t1), function(x) length(which(is.finite(x))) )

      # this is repeated to return to the same state as just prior to the netmind operations 
      # merging there would have been easier but it is possible to merge here to make things more modular

      nm = netmind.db( DS="stats" )
      set = merge( set, nm, by =c("trip","set"), all.x=TRUE, all.y=FALSE, suffixes=c("", ".nm") )
  
      # last resort: use netmind data to fill
      ii = which(!is.finite( set$t0) )
      if (length(ii) > 0 )  set$t0[ ii] = as.POSIXct( set$t0.nm[ii], origin=lubridate::origin, tz=tzone)
      set$t0.nm = NULL

      ii = which(!is.finite( set$t1) )
      if (length(ii) > 0 )  set$t1[ ii] = as.POSIXct( set$t1.nm[ii], origin=lubridate::origin, tz=tzone)
      set$t1.nm = NULL
   
      ii = which( !is.finite( set$dt) )
      if (length(ii) > 0 )  set$dt[ ii] =  set$dt.nm[ii]
      set$dt.nm = NULL

      # historical data do not have these fields filled .. fill
      ii = which( is.na( set$t0 ) )   
      if ( length (ii) > 0 ) {
        set$t0[ii] = set$timestamp[ii]
      }
      # fix t1
      ii = which( is.na( set$t1 ) )  # historical data do not have these fields filled .. fill 
      if ( length (ii) > 0 ) {
        set$t1[ii] = set$t0[ii] + median(set$dt, na.rm=TRUE )
      }

      # positional data obtained directly from Netmind GPS and Minilog T0
      # overwrite all, where available
      ilon = which( is.finite( set$slon)  )
      set$lon[ilon] = set$slon[ilon]
      
      ilat = which( is.finite( set$slat) )
      set$lat[ilat] = set$slat[ilat]

      set = lonlat2planar(set, proj.type=p$internal.projection) # get planar projections of lon/lat in km
      set$plon = grid.internal( set$plon, p$plons )
      set$plat = grid.internal( set$plat, p$plats )

      # merge surfacearea from net mesnuration into the database 
      set = clean.surface.area( set )
         
      set$slon = NULL
      set$slat = NULL
      set$Tx = NULL
      set$Zx = NULL
      set$observer = NULL
      set$cfa = NULL
      set$gear = NULL

    #  set2015= set[which(set$yr==2015),]
    #  print(head(set2015))
   
      save( set, file=fn, compress=TRUE )
      
      return(fn)
    }


    # --------------------------------
    

    if (DS %in% c("set", "set.complete", "set.complete.redo") ) {
      
      fn = file.path( project.datadirectory("snowcrab"), "R", "set.complete.rdata") 

      if (DS %in% c("set", "set.complete") ){
        load( fn )
        return ( set )
      }
      
      set = snowcrab.db( DS="set.merge.cat" )
      # set2015 = set[which(set$yr == 2015),]
      # print(head(set2015))

			# bring in time invariant features:: depth
			print ("Bring in depth")
      set = habitat.lookup( set,  p=p, DS="depth" )
      set$z = log( set$z )
			
		  # bring in time varing features:: temperature
			print ("Bring in temperature")
      set = habitat.lookup( set, p=p, DS="temperature" )

			# bring in all other habitat variables, use "z" as a proxy of data availability
			# and then rename a few vars to prevent name conflicts
		  set = habitat.lookup( set,  p=p, DS="all.data" )
		
      # return planar coords to correct resolution
      set = lonlat2planar( set, proj.type=p$internal.projection )
      
      # complete area designations  
      set = fishing.area.designations(set, type="lonlat")

      # ----add other species
      print( "Adding other species to 'set' ")
      cat = snowcrab.db( DS="cat.initial" )
      cat2015 = cat[grep("2015", cat$trip),]
      print(head(cat2015))
			
			cat$uid = paste(cat$trip, cat$set, sep="~")
			set$uid = paste(set$trip, set$set, sep="~")
			suid = unique(sort( set$uid)) # necessary as some sp are found in sets that are dropped (bad tows)

			ns = nrow(set)

      for ( i in sort( unique( cat$spec ) ) ) {
        print(i)
        tmp = NULL
        tmp = cat[ which(cat$spec==i & cat$uid %in% suid ) , c("uid","totno","totmass")  ]
        tmp$meansize = tmp$totmass / tmp$totno
        names(tmp) = c("uid", paste( c("ms.no", "ms.mass", "ms.size"), i, sep="." ) )
        o = merge( set, tmp, by=c("uid"), all.x=T, all.y=F, sort=F )
				if ( nrow(o) == nrow(set) ) {
					set = o
				} else {
					print (nrow(o))
					stop()
				}
			}
      
      j = unique( c(grep("ms.mass", names(set)), grep("ms.no.", names(set)) ))
      for ( k in j ) {
        l = which( !is.finite( set[,k] ) )
        set[l,k] = 0
        set[,k] = set[,k] / set$sa
      }
      
      set2015 = set[which(set$yr ==2015), ]
      print(head(set2015))
      
      save(set, file=fn, compress=T)
    
    }

    if (DS %in% c("det.georeferenced", "det.georeferenced.redo" ) ) {
      fn = file.path( project.datadirectory("snowcrab"), "R", "det.georef.rdata" )  
      if (DS=="det.georeferenced") {
        load(fn)
        return(det)
      }
      set = snowcrab.db( "set.clean")
      set  = set[, c("trip", "set", "lon", "lat", "plon", "plat", "yr")]
      det = snowcrab.db("det.initial")
      det = merge( det, set, by=c("trip", "set"), all.x=T, all.y=F, sort=F, suffixes=c("",".set") )
      det$sa.set = NULL
      save(det, file=fn,compress=T)
    }
   
    if (DS %in% c("cat.georeferenced", "cat.georeferenced.redo" ) ) {
      fn = file.path( project.datadirectory("snowcrab"), "R", "cat.georef.rdata" )  
      if (DS=="cat.georeferenced") {
        load(fn)
        return(cat)
      }

      set = snowcrab.db( "set.clean")  #require SA estimates
      set  = set[, c("trip", "set", "lon", "lat", "plon", "plat", "yr", "sa")]
      cat =  snowcrab.db("cat.initial")
      cat = merge( cat, set, by=c("trip", "set"), all.x=T, all.y=F, sort=F, suffixes=c("",".set") )
      cat$totmass = cat$totmass / cat$sa
      cat$totno = cat$totno / cat$sa

      cat$sa.set = NULL
      save(cat, file=fn,compress=T)
    }

    
    if (DS %in% c("set.logbook", "set.logbook.redo" )) {
      
      outdir = file.path( project.datadirectory("snowcrab"), "data" )
      dir.create(path=outdir, recursive=T, showWarnings=F)
      fn = file.path( outdir, "set.logbook.rdata" ) 
      
      if (DS =="set.logbook") {
        if ( file.exists( fn) ) load(fn)
        return(set)
      }
   
      set0 = snowcrab.db( DS="set.complete" )
      set = logbook.fisheries.stats.merge( set0 )
      save ( set, file=fn, compress=T )
      return (fn)
    }
  

  }  ## end snowcrab.db


