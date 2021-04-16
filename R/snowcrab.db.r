
snowcrab.db = function( DS, p=NULL, yrs=NULL, fn.root=NULL, redo=FALSE, extrapolation_limit=NA, extrapolation_replacement="extrapolation_limit", ... ) {

	# handles all basic data tables, etc. ...

  # sex codes
  male = 0
  female = 1
  sex.unknown = 2

  # maturity codes
  immature = 0
  mature = 1
  mat.unknown = 2


	if (DS %in% c("set.rawdata.redo", "set.rawdata") ) {
    if (is.null(fn.root)) fn.root =  file.path( project.datadirectory("bio.snowcrab"), "data", "trawl", "SNCRABSETS" )
		dir.create( fn.root, recursive = TRUE, showWarnings = FALSE )

    if (DS=="set.rawdata") {
			out = NULL
			fl = list.files( path=fn.root, pattern="*.rdata", full.names=T )
      for ( fny in fl ) {
				load (fny)
				out = rbind( out, SNCRABSETS )
			}
			return (out)
		}

		con=ROracle::dbConnect(DBI::dbDriver("Oracle"),dbname=oracle.snowcrab.server , username=oracle.snowcrab.user, password=oracle.snowcrab.password, believeNRows=F)
					# believeNRows=F required for oracle db's

		for ( YR in yrs ) {
			fny = file.path( fn.root, paste( YR,"rdata", sep="."))
			SNCRABSETS = NULL
			SNCRABSETS = ROracle::dbGetQuery(con, paste("select * from SNCRABSETS
			                                            where EXTRACT(YEAR from BOARD_DATE) = ", YR , "
			                                            OR (EXTRACT(YEAR from BOARD_DATE) = ", YR+1 , " AND EXTRACT(MONTH FROM Board_DATE)=1)") )
			save( SNCRABSETS, file=fny, compress=TRUE)
			gc()  # garbage collection
			print(YR)
		}

		ROracle::dbDisconnect(con)
		return (yrs)
	}


  # -------------------------------


	if (DS %in% c("det.rawdata.redo", "det.rawdata") ) {
    if (is.null(fn.root)) fn.root =  file.path( project.datadirectory("bio.snowcrab"), "data", "trawl", "SNCRABDETAILS" )
		dir.create( fn.root, recursive = TRUE, showWarnings = FALSE  )

    if (DS=="det.rawdata") {
			out = NULL
      fl = list.files( path=fn.root, pattern="*.rdata", full.names=TRUE )
			for ( fny in fl ) {
				load (fny)
				out = rbind( out, SNCRABDETAILS )
			}
			return (out)
		}

		con=ROracle::dbConnect(DBI::dbDriver("Oracle"),dbname=oracle.snowcrab.server , username=oracle.snowcrab.user, password=oracle.snowcrab.password, believeNRows=F)

		for ( YR in yrs ) {
			fny = file.path( fn.root, paste( YR,"rdata", sep="."))
			SNCRABDETAILS = NULL
			#in following line replaced sqlQuery (Rrawdata) with  ROracle::dbGetQuery (ROracle)
			SNCRABDETAILS = ROracle::dbGetQuery(con,
                paste("select * from SNCRABDETAILS
                      where EXTRACT(YEAR from BOARD_DATE) = ", YR , "
			                                            OR (EXTRACT(YEAR from BOARD_DATE) = ", YR+1 , " AND EXTRACT(MONTH FROM Board_DATE)=1)") )
			save( SNCRABDETAILS, file=fny, compress=TRUE)
			gc()  # garbage collection
			print(YR)
		}

		ROracle::dbDisconnect(con)
    return (yrs)

	}

  # -------------------------------

	if (DS %in% c("cat.rawdata.redo", "cat.rawdata") ) {

    if (is.null(fn.root)) fn.root =  file.path( project.datadirectory("bio.snowcrab"), "data", "trawl", "SNTRAWLBYCATCH" )
		dir.create( fn.root, recursive = TRUE, showWarnings = FALSE  )

    if (DS=="cat.rawdata") {
			out = NULL
      fl = list.files( path=fn.root, pattern="*.rdata", full.names=TRUE )
			for ( fny in fl ) {
				load (fny)
				out = rbind( out, SNTRAWLBYCATCH )
			}
			return (out)
		}

		con=ROracle::dbConnect(DBI::dbDriver("Oracle"),dbname=oracle.snowcrab.server , username=oracle.snowcrab.user, password=oracle.snowcrab.password, believeNRows=F)

		for ( YR in yrs ) {
			fny = file.path( fn.root, paste( YR,"rdata", sep="."))
			SNTRAWLBYCATCH = NULL
			#in following line replaced sqlQuery (Rrawdata) with  ROracle::dbGetQuery (ROracle)
			SNTRAWLBYCATCH = ROracle::dbGetQuery(con,
                paste("select * from SNTRAWLBYCATCH
                      where EXTRACT(YEAR from BOARD_DATE) = ", YR , "
			                                            OR (EXTRACT(YEAR from BOARD_DATE) = ", YR+1 , " AND EXTRACT(MONTH FROM Board_DATE)=1)") )
			save( SNTRAWLBYCATCH, file=fny, compress=TRUE)
			gc()  # garbage collection
			print(YR)
		}

		ROracle::dbDisconnect(con)
    return (yrs)
	}


	# ---------------------


  if (DS %in% c("setInitial.redo", "setInitial")) {
    fn = file.path( project.datadirectory( "bio.snowcrab", "data" ), "set.initial.rdata" )

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
    # August 2015 added in setcd_id from observer system to address the MPA survey sets (type=11) and regular fix station sets (type=4) .. renamed to set_type
    set = snowcrab.db( DS="set.rawdata")
    names( set ) = rename.bio.snowcrab.variables(names( set))
    setvars = c("trip", "set", "set_type", "station", "stime", "observer", "cfa", "lon", "lat", "lon1", "lat1", "towquality", "Zx", "Tx", "gear", "sa", "dist", "dist0" )
    print('need to addin the mpa station index')

    set$trip = as.character(set$trip)
    set$set  = as.integer(set$set)

    set = set[ order( set$sdate ), ]

    # first estimate of distance of trawl track based upon logged start and end locations
    set$dist0 = geosphere::distGeo( set[,c("lon", "lat" )], set[,c( "lon1", "lat1") ] )  # m
    ii = which( set$dist0 > 1000 | set$dist0 < 100 )
    if ( length(ii) > 0 ) {
      print( "Positional information from 'set' incorrect or suspect in the following as the naive tow distance, 'dist0' seems out of range (100 to 1000m):" )
      oo = set[ ii, c("trip", "set", "dist", "dist0", "lon", "lat", "lon1", "lat1", "towquality") ]
      print( oo )
    }
    set$dist0[ii] = NA

    # vs set$dist which are hand computed by Moncton for the early time series
    #  plot( dist0 ~ dist, set) shows good correspondence

    ############
    # As a simple double check, totmass cannot be 0 on a good tow (quality 1) or >0 on a bad tow (quality 4)
    oo = which( set$towquality ==4 & set$totmass >0 )
    if (length(oo)>1) {
      print( "The following needs some checking as there is a mass estimate but towquality set to be bad" )
      print( set[ oo, ] )
    }

  # The following section were "on the fly" fixes to minor problems before data base corrections can be made
  # 2018- BZ. These have all been corrected in the database. One "debug" is left as an example in case needed in the future

    #dbug.2011 = TRUE
    #if ( dbug.2011 ) {
      # one-off error corrections until database gets refreshed
     # i = which(set$trip=="S15092006" & set$set==3)
      #if (length(i)==1) set$totmass[i] = 0

      #i = which(set$trip=="S21092007" & set$set==12)
      #if (length(i)==1) set$towquality[i] = 1
    #}

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
    sntows = read.table(file.path( project.datadirectory("bio.snowcrab"), "data", "lookuptables", "sntow.csv"),  sep=";", as.is=TRUE, colClasses="character", header=TRUE)
    sntow.vars = c("trip",  "set", "patchtime", "netmind", "minilog")

    set = merge(set, sntows[,sntow.vars], by=c("trip", "set"), all.x=TRUE, all.y=FALSE, sort=FALSE)
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

    # "timestamp" is the best estimate of sampling time
    # sdate (POSIXct, America/Halifax AST/ADT) does not seem to be reliable
    # and so we use minilog data where possible in the recent period
    # and then records from the stime and trip id where minilog data are absent
    set$timestamp = tripcode.to.timestamp( set$trip, set$stime )  # using lubridate/POSIXct
    set$stime = NULL ### --need to check timezone!!! TODO .... seems to be "America/Halifax" .. compare with seabird/minilog

    i = which(is.na(set$timestamp))
    if (length(i)>0) set$timestamp[i] = tripcode.to.timestamp( set$trip[i], "12:00:00" )

    # from this point forwards all timestamps are in UTC for set
    set$timestamp = with_tz( set$timestamp, "UTC")

    set$julian = lubridate::yday( set$timestamp )
    set$yr = lubridate::year( set$timestamp )  # "survey year"

    # some surveys extend into January (e.g., 2020) force them to be part of the correct "survey year", i.e., "yr"
    i = which(lubridate::month(set$timestamp)==1)
    if (length(i) > 0) set$yr[i] = set$yr[i] - 1

    # set$timestamp[i] = set$timestamp[i] - 1382400  # should not touch timestmp as this is a key index

    save( set, file=fn, compress=TRUE )

    return ( fn )
  }


  # --------------------

  if (DS %in% c("det.initial", "det.initial.redo") ) {
    fn = file.path(project.datadirectory("bio.snowcrab"), "data", "det.initial.rdata")
    if (DS =="det.initial") {
      load(fn)
      return(det)
    }

    X = snowcrab.db( DS="set.clean" )
    det = snowcrab.db( DS="det.rawdata"  )


    names( det ) = rename.bio.snowcrab.variables(names(det) )
    detvars = c( "trip", "set", "crabno", "sex", "cw", "mass", "abdomen", "chela", "mat",
                 "shell", "gonad", "eggcol", "eggPr", "durometer", "legs")
    det=det[is.finite(det$crabno),]
    # merge in the sa which is used as a weighting factor of most analyses
    det = merge(x=det[,detvars], y=X[,c("trip","set","sa")], by=c("trip", "set"), all.x=TRUE, all.y=FALSE)

    # Trips and sets with no matching SA ...
    ii = which( !is.finite(det$sa) )
    #print ("DET without matching SA ... due to bad tows? ... check these ")
    #print (det[ii, c("trip","set")])
    #print (det[ii,])
    det$sa[ii] = median(det$sa, na.rm=TRUE )

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
    yr.e <- p$year.assessment
    fn.e = file.path(project.datadirectory("bio.snowcrab"), "output", "morphology.errors")
    dir.create(fn.e, recursive=TRUE, showWarnings=F)
    outfile.e =  file.path( fn.e, paste("morphologyerrors", yr.e, ".csv", sep=""))
    outfile.e2 =  file.path( fn.e, paste("morphologyerrors.allyears", yr.e, ".csv", sep=""))

    #Sex.e: Unknown Sex
    sex.e <- det[which(det$sex==sex.unknown),]
    if ( !is.na(sex.e$trip[1])) sex.e$error <- 'sex.e'
    #Cw.e: Carapace Width below 5 or greater than 185
    cw.e <- det[ which(det$cw<5 | det$cw>185 ),]
    if ( !is.na(cw.e$trip[1])) cw.e$error <- 'cw.e'
    #Chela.e: Chela less than 1 or greater than 50
    chela.e <- det[which(det$chela < 1 | det$chela > 50  ),]
    if ( !is.na(chela.e$trip[1])) chela.e$error <- 'chela.e'
    #Abdomen.e:Abdomen less than 1 and greater than 66
    abdomen.e <- det[which(det$abdomen < 1 | det$abdomen > 66 ),]
    #abdomen.e$error <- 'abdomen.e' #BZ 2018 no abdomen lengths met "error" condition, broke script #
    if ( !is.na(abdomen.e$trip[1]))  abdomen.e$test='abdomen.e' #replaced above statement

    #Mass.e: Mass less than 1 or greater than 1500
    mass.e <- det[which( det$mass < 1 | det$mass > 1500  ),]
    if ( !is.na(mass.e$trip[1]))  mass.e$error <- 'mass.e'

    #Sex.a: Indeterminate sex based on measurements taken (abdomen values where sex=male)
    sex.a <- det[which(is.finite( det$abdomen ) & det$sex==male),]
    if ( !is.na(sex.a$trip[1])) sex.a$error <- 'sex.a'
    #Sex.c: Indeterminate sex based on measurements taken (chela values where sex=female
    sex.c <- det[which(is.finite( det$chela ) & det$sex==female),]
    if ( !is.na(sex.c$trip[1])) sex.c$error <- 'sex.c'


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
    mat.e <- det[which(det$mat ==2 & (is.finite(det$chela) | is.finite(det$abdomen))),]
    if ( !is.na(mat.e$trip[1])) mat.e$error <- 'mat.e'


    primiparous = filter.class( det, "primiparous")
    multiparous = filter.class( det, "multiparous")

    det$fecundity = NA
    det$fecundity[primiparous] = fecundity.allometric( cw=det$cw[primiparous], method="moncton.primiparous" )
    det$fecundity[multiparous] = fecundity.allometric( cw=det$cw[multiparous], method="moncton.multiparous" )
    det$fecundity[ which(det$fecundity> 250000) ] = NA

    save(det, file=fn, compress=TRUE)

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
    errors <<- errors
    message("check dataframe 'errors' for the errors")
    if ( !is.na(errors.yearly$trip[1]))  {
    print(errors.yearly)
    write.csv(errors.yearly, file=outfile.e)
    print("Current Year Morphology Errors saved to file")
    print(outfile.e)
}

    write.csv(errors, file=outfile.e2)
    print("All Years Morphology Errors saved to file")
    print(outfile.e2)

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
    fn = file.path(project.datadirectory("bio.snowcrab"), "data", "cat.initial.rdata")
    if(DS =="cat.initial" ) {
      load(fn)
      return(cat)
    }

		# two steps:: bycatch from the cat tables (missing snow crab)
    # and determine totals from the snow crab det tables

    det = snowcrab.db( DS="det.initial" )

    cat = snowcrab.db( DS="cat.rawdata" )
    names( cat ) = rename.bio.snowcrab.variables(names( cat ) )

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
    # --- no longer here ... only when integrated into survey_db

    # remove data where species codes are ambiguous, or missing or non-living items
    xx = which( !is.finite( cat$spec) )
    if (length(xx)>0) cat = cat[ -xx, ]
    cat = cat[ taxonomy.filter.taxa( cat$spec, taxafilter="living.only", outtype="rvsurveycodes" ) , ]

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
    # longer here -- in survey_db only

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
    gp = aegis.survey::groundfish_parameters( year.assessment=p$year.assessment )
    gs = aegis.survey::groundfish_survey_db( p=gp, DS="gscat" )  # raw data .. does not need vessel corrections
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

    if ( DS=="areal_units_input" ) {

      fn = file.path( p$datadir,  "areal_units_input.rdata" )
      if ( !file.exists(p$datadir)) dir.create( p$datadir, recursive=TRUE, showWarnings=FALSE )

      xydata = NULL
      if (!redo)  {
        if (file.exists(fn)) {
          load( fn)
          return( xydata )
        }
      }

      xydata = snowcrab.db( p=p, DS="set.clean"  )  #
      xydata = xydata[ , c("lon", "lat", "yr" )]
      save(xydata, file=fn, compress=TRUE )
      return( xydata )
    }


  # --------------------------------


  if ( DS %in% c("set.clean", "set.clean.redo") ) {

    # merge seabird, minilog and netmind data and do some checks and cleaning
    fn = file.path( project.datadirectory( "bio.snowcrab" ), "data", "set.clean.rdata" )

    if ( DS=="set.clean" ) {
      set= NULL
      if (file.exists( fn) ) load( fn )
      return (set)
    }

    # the beginning here is identical to the netmind.db( "stat.redo" ) .. simpler to keep it this way (jae)
    set = snowcrab.db( DS="setInitial")  # timestamp in UTC
    nI = nrow(set)

    sbStats =  seabird.db( DS="stats" )  # timestamp in UTC

    sbv = c('trip','set', "z", "zsd", "t", "tsd", "n", "t0", "t1", "dt", "seabird_uid" )
    set_sb = merge( set[, c("trip", "set") ], sbStats[,sbv], by=c("trip","set"), all.x=TRUE, all.y=FALSE, sort=FALSE )
    # tapply( as.numeric(set_sb$dt), year(set_sb$t1), mean, na.rm=T )
    # tapply( as.numeric(set_sb$dt), year(set_sb$t1), function(x) length(which(is.finite(x))) )
    if (nrow(set) !=nI ) stop( "merge error with seabird0" )

    mlStats =  minilog.db( DS="stats" )
    ids = paste(mlStats$trip, mlStats$set, sep=".")
    uu = which( duplicated( ids ) )
    if (length(uu)>0 ) {
      message( "Duplicated minilog data (mlStats) trip/set:" )
      toshow = which( ids %in% ids[uu] )
      print( mlStats[ toshow,] )
      message("Dropping for now ... ")
      mlStats = mlStats[-uu,]
    }

     # mlStats$dt = as.numeric(mlStats$dt )
    mlv =  c('trip', 'set', "z",    "zsd",    "t",    "tsd",    "n",    "t0",    "t1",    "dt", "minilog_uid" )
    set_ml = merge( set[, c("trip", "set") ], mlStats[,mlv], by=c("trip","set"), all.x=TRUE, all.y=FALSE, sort=FALSE )
    # tapply( as.numeric(set_ml$dt), lubridate::year(set_ml$t1), mean, na.rm=T )
    # tapply( as.numeric(set_ml$dt), year(set_ml$t1), function(x) length(which(is.finite(x))) )

    if (nrow(set) !=nI ) stop( "merge error with minilogs0" )

    set = merge( set, set_sb, by=c("trip", "set" ), all.x=TRUE, all.y=FALSE, sort=FALSE )
    if (nrow(set) !=nI ) stop( "merge error with seabird" )

    set = merge( set, set_ml, by=c("trip", "set" ), all.x=TRUE, all.y=FALSE, sort=FALSE, suffixes=c("", ".ml" ))
    if (nrow(set) !=nI ) stop( "merge error with minilogs" )

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
    if (nrow(set) !=nI ) stop( "merge error with netmind" )

    # last resort: use netmind data to fill
    ii = which(!is.finite( set$t0) )
    # if (length(ii) > 0 )  set$t0[ ii] = as.POSIXct( set$t0.nm[ii], origin=lubridate::origin, tz="UTC" )

    if (length(ii) > 0 )  set$t0[ ii] = set$t0.nm[ii]
    set$t0.nm = NULL

    ii = which(!is.finite( set$t1) )
    # if (length(ii) > 0 )  set$t1[ ii] = as.POSIXct( set$t1.nm[ii], origin=lubridate::origin, tz="UTC")
    if (length(ii) > 0 )  set$t1[ ii] =  set$t1.nm[ii]
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

    set = lonlat2planar(set, proj.type=p$aegis_proj4string_planar_km) # get planar projections of lon/lat in km

    grid = spatial_grid(p=p, DS="planar.coords")

    message("probably do not need to grid any longer")

    set$plon = grid_internal( set$plon, grid$plon )
    set$plat = grid_internal( set$plat, grid$plat )

    # merge surfacearea from net mesnuration into the database
    set = clean.surface.area( set, qreject = c( 0, 1 ))

    zmod = glm( Zx ~ z - 1, data=set)
    zres = residuals( zmod)
    # hist(abs(zres), "fd")
    not.reliable = which( abs(zres) > 25 )
    set$z[not.reliable] = NA  # force these to a default lookup from from bathymetry_db

    set$slon = NULL
    set$slat = NULL
    set$Tx = NULL
    set$Zx = NULL
    set$observer = NULL
    set$cfa = NULL
    set$gear = NULL

    save( set, file=fn, compress=TRUE )
    return(fn)
  }


  # --------------------------------


  if (DS %in% c("det.georeferenced", "det.georeferenced.redo" ) ) {
    fn = file.path( project.datadirectory("bio.snowcrab"), "data", "det.georef.rdata" )
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


  # -------------------------------


  if (DS %in% c("cat.georeferenced", "cat.georeferenced.redo" ) ) {
    fn = file.path( project.datadirectory("bio.snowcrab"), "data", "cat.georef.rdata" )
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



  # -------------

  if ( DS %in% c("set.biologicals", "set.biologicals.redo") ) {

    fn = file.path( project.datadirectory("bio.snowcrab"), "data", "set.biologicals.rdata")

    if (DS=="set.biologicals" ) {
      load(fn)
      return( set)
    }

    factors = c("trip", "set")

    X = snowcrab.db( DS="set.clean" )
    nsInit = nrow(X)

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

    if ( nsInit != nrow( set) ) {   print( "Merge failure 1... " );  stop()    }

    # X2015 = X[which(X$yr == 2015),]
    # print(head(X2015))

    cat = snowcrab.db( DS="cat.initial" )
    # cat2015 = cat[grep("2015", cat$trip),]
    # print(head(cat2015))

    cat0 = cat[ taxonomy.filter.taxa( cat$spec, taxafilter="snowcrab", outtype="rvsurveycodes"), c(factors, "totno")]
    names(cat0) = c(factors, "totno.all")
    X = merge(x=X, y=cat0, by=factors, all.x=T )
    X$totno.all   = X$totno.all   / X$sa
    X$totno.all[!is.finite(X$totno.all)] = 0  # convert na's to zero

    cat0 = cat[ taxonomy.filter.taxa( cat$spec, taxafilter="snowcrab", outtype="rvsurveycodes" ), c(factors, "totmass")]
    names(cat0) = c(factors, "totmass.all")
    X = merge(x=X, y=cat0, by=factors, all.x=T )
    X$totmass.all = X$totmass.all / X$sa
    X$totmass.all[!is.finite(X$totmass.all)] = 0  # convert na's to zero
    # the above masses are in kilograms .. convert to metric tons
    var="totmass.all";  X[,var] = X[,var] / 10^3

    set = X
    if ( nsInit != nrow( set) ) {   print( "Merge failure 2... " );  stop()    }

    # complete area designations
    set = fishing.area.designations(set, type="lonlat")

    # ----add other species
    print( "Adding other species to 'set' ")
    cat = snowcrab.db( DS="cat.initial" )
    # cat2015 = cat[grep("2015", cat$trip),]
    # print(head(cat2015))

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
    if ( nsInit != nrow( set) ) {   print( "Merge failure 3... " );  stop()    }

    j = unique( c(grep("ms.mass", names(set)), grep("ms.no.", names(set)) ))
    for ( k in j ) {
      l = which( !is.finite( set[,k] ) )
      set[l,k] = 0
      set[,k] = set[,k] / set$sa
    }

    if ( nsInit != nrow( set) ) {   print( "Merge failure ... " );  stop()    }
    # X2015 = X[which(X$yr == 2015),]
    # print(head(X2015))

    save( set, file=fn, compress=T )

    return ( "Complete" )
  }

  # -------------------------------


  if (DS %in% c( "set.complete", "set.complete.redo") ) {

    fn = file.path( project.datadirectory("bio.snowcrab"), "data", "set.complete.rdata")

    if (DS %in% c("set", "set.complete") ){
      load( fn )
      return ( set )
    }

    set = snowcrab.db( DS="set.biologicals" )
    # set2015 = set[which(set$yr == 2015),]
    # print(head(set2015))
    # return planar coords to correct resolution
    set = lonlat2planar( set, proj.type=p$aegis_proj4string_planar_km )

    # bring in time invariant features:: depth
    ii = which(!is.finite(set$z))
    if (length(ii)>0){
      set$z[ii] =  bathymetry_lookup( LOCS=set[ ii, c("lon", "lat")],  lookup_from="core", lookup_to="points" , lookup_from_class="aggregated_data" ) # core=="rawdata"
 
    }

    set$z = log( set$z )
    # as of 2016, there are 11 locations where there are missing depths, because they are outside the area defined for snow crab ... they are all bad sets too (set_type=4) in NENS ... ignoring for now

    # bring in time varing features:: temperature
    ii = which(!is.finite(set$t))
    if (length(ii)>0){
      set$t[ii] = temperature_lookup( LOCS=set[ ii, c("lon", "lat", "timestamp")],lookup_from="core", lookup_to="points", lookup_from_class="aggregated_data", tz="America/Halifax"  )

    }


    # set2015 = set[which(set$yr ==2015), ]
    # print(head(set2015))

    set = logbook.fisheries.stats.merge( set ) # add gridded logbook data

    if ( nrow( snowcrab.db( DS="setInitial" )) != nrow( set) ) {   print( "Merge failure ... " );  stop()    }

    save(set, file=fn, compress=T)

  }

  # -------------------------------


  if (DS %in% c("data.transforms", "data.transforms.redo") ) {
    REPOS = NULL
    if (is.null(p)) p = bio.snowcrab::snowcrab_parameters()

    if (DS=="data.transforms") {
      if (file.exists( p$transform_lookup ) ) load (p$transform_lookup)
      return(REPOS)
    }

    log.transform = bio.snowcrab::snowcrab.variablelist("log.transform")
    sn = bio.snowcrab::snowcrab.variablelist("all.data")
    set = bio.snowcrab::snowcrab.db(DS="set.complete")
    logs = bio.snowcrab::logbook.db(DS='logbook')
    scaled.centered = bio.snowcrab::snowcrab.variablelist("scaled.centered")

    dataset.names = unique( c(names(set), names(logs)) )

    for (si in 1:length(sn)) {
      transform = offset = scaling =NA
      varname = sn[si]
      if (! varname %in% dataset.names ) next()
      if(varname %in% names(set))  x = set[, varname]
      if(varname %in% names(logs)) x = logs[, varname]
      if (varname %in% log.transform) {
        transform="log10"
        offset = 0
        y = x[ is.finite(x) & x>0 ]
        if( length(y) > 0) offset = min(y)
      } else if (varname %in% scaled.centered) {
        transform = "scaled+centered"
        y = scale( x )
        offset = attr(y,"scaled:center") # mean
        scaling = attr(y,"scaled:scale") # RMS error  .. i.e. a Z-transform
      } else {
        transform = "none"
        y = x
        offset = 0
        scaling = 1
      }
      # add more as needed
      REPOS = rbind( REPOS,  cbind( varname, transform, offset, scaling  )
      )
    }
    REPOS = data.frame( REPOS, stringsAsFactors=F )
    REPOS$offset = as.numeric(REPOS$offset)
    REPOS$scaling = as.numeric(REPOS$scaling)

    save( REPOS, file=p$transform_lookup, compress=TRUE )
    return( REPOS )
  }


  # ----------------------



  if ( DS=="biological_data") {

    set = aegis.survey::survey_db( p=p, DS="filter" ) 
    
    if ( p$selection$type=="number") {
      # should be snowcrab survey data only taken care of p$selection$survey = "snowcrab"
      # robustify input data: .. upper bound trim
      if (exists("quantile_bounds", p)) {
        highestpossible = quantile( set$totno_adjusted, probs=p$quantile_bounds[2], na.rm=TRUE )
        set$totno_adjusted[ set$totno_adjusted > highestpossible ] = highestpossible
        # keep "zero's" to inform spatial processes but only as "lowestpossible" value
        jj = which( set$totno_adjusted > 0 )
        lowestpossible =  quantile( set$totno_adjusted[jj], probs=p$quantile_bounds[1], na.rm=TRUE )
        ii = which( set$totno_adjusted < lowestpossible )
        set$totno_adjusted[ii] = 0
      }
      set$data_offset  = 1 / set[, "cf_set_no"]
    }

    if ( p$selection$type=="biomass") {
      # should be snowcrab survey data only taken care of p$selection$survey = "snowcrab"
      # robustify input data: .. upper bound trim
      if (exists("quantile_bounds", p)) {
        highestpossible = quantile( set$totwgt_adjusted, probs=p$quantile_bounds[2], na.rm=TRUE )
        set$totwgt_adjusted[ set$totwgt_adjusted > highestpossible ] = highestpossible

        # keep "zero's" to inform spatial processes but only as "lowestpossible" value
        jj = which( set$totwgt_adjusted > 0 )
        lowestpossible =  quantile( set$totwgt_adjusted[jj], probs=p$quantile_bounds[1], na.rm=TRUE )
        ii = which( set$totwgt_adjusted < lowestpossible )
        set$totwgt_adjusted[ii] = 0 ## arbitrary but close to detection limit
      }
      set$data_offset  = 1 / set[, "cf_set_mass"]
    }

    if ( p$selection$type=="presence_absence") {
      # must run here as we need the wgt from this for both PA and abundance
      if ( "logbook" %in% p$selection$survey$data.source ) {

        if (p$selection$biologicals$sex == 0 &
            p$selection$biologicals$mat == 1 &
            min(p$selection$biologicals$len) >= 95/10 &
            max(p$selection$biologicals$len) <= 200/10
        ) {
          # add commerical fishery data --
          # depth data is problematic ... drop for now
          lgbk = logbook.db( DS="fisheries.complete", p=p )
          lgbk = lgbk[ which( is.finite( lgbk$landings)), ]
          lgbk = lgbk[ which( lgbk$year > 2005), ]  # previous to this all sorts of traps were used
          lgbk = lgbk[ which( as.numeric(lgbk$soak.time) >= 12 & as.numeric(lgbk$soak.time) <= 48), ]   # avoid nonlinearity in catch with time
          lgbk$cpue_time = lgbk$cpue / as.numeric(lgbk$soak.time)  # approx with realtive catch rate in time

          lgbk$qm = NA   # default when no data
          oo = which( lgbk$cpue_time == 0 )  # retain as zero values
          if (length(oo)>0 ) lgbk$qm[oo] = 0
          ii = which( lgbk$cpue_time != 0 )
          lgbk$qm[ii] = quantile_estimate( lgbk$cpue_time[ii]  )  # convert to quantiles
          lgbk$zm = quantile_to_normal( lgbk$qm )
          lgbk$data.source = "logbook"
          lgbk$z = exp( lgbk$z )

          # transparently create NA filled vars to pass all variables through
          missingvars = setdiff( names(set) , names( lgbk) )
          for (nn in missingvars) lgbk[,nn] = NA
          nms = intersect( names(set) , names( lgbk) )
          set = rbind( set[, nms], lgbk[,nms] )
        }
      }

      pa = presence.absence( X=set$zm, px=p$habitat.threshold.quantile )  # determine presence absence and weighting
      set[, p$variabletomodel] = pa$pa
      set$data_offset  = pa$probs  # just a dummy value to make sure offsets are filled (with another dummy value)
      pa = NULL
    }

    set = set[ which(is.finite(set$plon + set$plat)),]

    crs_lonlat = st_crs(projection_proj4string("lonlat_wgs84"))
    inside = st_points_in_polygons(
      pts = st_as_sf( set[, c("lon", "lat")], coords=c("lon","lat"), crs=crs_lonlat ),
      polys = st_transform( coastline_db( p=p ), crs_lonlat )
    )

    onland =which (is.finite(inside))
    if (length(onland)>0) set = set[-onland, ]

    set$tiyr = lubridate::decimal_date( set$timestamp )

    ## lookup data

    # set = aegis_db_lookup(
    #   X=set,
    #   lookupvars=p$stmv_variables$COV,
    #   xy_vars=c("lon", "lat"),
    #   time_var="timestamp"
    # )

    # if (!alldata) {
    #  set = set[, which(names(set) %in% c( p$stmv_variables$LOCS, p$stmv_variables$COV, p$stmv_variables$Y, p$stmv_variables$TIME, "dyear", "yr",  "wt") ) ]  # a data frame
    #   oo = setdiff( c( p$stmv_variables$LOCS, p$stmv_variables$COV ), names(set))
    #   if (length(oo) > 0 ) {
    #     print(oo )
    #     warning("Some variables are missing in the input data")
    #   }
    #   set = na.omit(set)
    # }

    # cap quantiles of dependent vars
    # if (exists("quantile_bounds", p)) {
    #   dr = list()
    #   for (pvn in p$stmv_variables$COV) {
    #     dr[[pvn]] = quantile( set[,pvn], probs=p$quantile_bounds, na.rm=TRUE ) # use 95%CI
    #     il = which( set[,pvn] < dr[[pvn]][1] )
    #     if ( length(il) > 0 ) set[il,pvn] = dr[[pvn]][1]
    #     iu = which( set[,pvn] > dr[[pvn]][2] )
    #     if ( length(iu) > 0 ) set[iu,pvn] = dr[[pvn]][2]
    #   }
    # }

    #set$Y = set$totno  # unadjusted value is used as we are usinmg offsets ...

    set$data_offset[which(!is.finite(set$data_offset))] = median(set$data_offset, na.rm=TRUE )  # just in case missing data
    set$wt = set$data_offset

    return( set )
  }


  # -------------------------


  if ( DS=="carstm_inputs") {

    # prediction surface
    crs_lonlat = st_crs(projection_proj4string("lonlat_wgs84"))
    sppoly = areal_units( p=p )  # will redo if not found
    sppoly = st_transform(sppoly, crs=crs_lonlat )
    areal_units_fn = attributes(sppoly)[["areal_units_fn"]]

    fn = carstm_filenames( p=p, returntype="carstm_inputs", areal_units_fn=areal_units_fn )

    # inputs are shared across various secneario using the same polys
    #.. store at the modeldir level as default
    outputdir = dirname( fn )
    if ( !file.exists(outputdir)) dir.create( outputdir, recursive=TRUE, showWarnings=FALSE )

    if (!redo)  {
      if (file.exists(fn)) {
        message( "Loading previously saved carstm_inputs ... ", fn)
        load( fn )
        return( M )
      }
    }
    message( "Generating carstm_inputs ... ", fn)


    # do this immediately to reduce storage for sppoly (before adding other variables)
    M = snowcrab.db( p=p, DS="biological_data" )  # will redo if not found .. not used here but used for data matching/lookup in other aegis projects that use bathymetry
  
    # some survey timestamps extend into January (e.g., 2020) force them to be part of the correct "survey year", i.e., "yr"
    i = which(lubridate::month(M$timestamp)==1)
    if (length(i) > 0) M$timestamp[i] = M$timestamp[i] - lubridate::duration(month=1)

    M$tiyr = lubridate::decimal_date(M$timestamp)

    # reduce size
    M = M[ which( M$lon > p$corners$lon[1] & M$lon < p$corners$lon[2]  & M$lat > p$corners$lat[1] & M$lat < p$corners$lat[2] ), ]
    # levelplot(z.mean~plon+plat, data=M, aspect="iso")

    M$AUID = st_points_in_polygons(
      pts = st_as_sf( M, coords=c("lon","lat"), crs=crs_lonlat ),
      polys = sppoly[, "AUID"],
      varname = "AUID"
    )
    M = M[!is.na(M$AUID),]
    M$AUID = as.character( M$AUID )  # match each datum to an area


    names(M)[which(names(M)=="yr") ] = "year"
    # M = M[ which(M$year %in% p$yrs), ]
    # M$tiyr = lubridate::decimal_date ( M$timestamp )
    # M$dyear = M$tiyr - M$year



    # --------------------------
     # bathymetry lookup
    pB = bathymetry_parameters( p=parameters_reset(p), project_class="carstm"  )
    vnB = pB$variabletomodel
    if ( !(exists(vnB, M ))) {
      vnB2 = paste(vnB, "mean", sep=".")
      if ((exists(vnB2, M ))) {
        names(M)[which(names(M) == vnB2 )] = vnB
      } else {
        M[,vnB] = NA
      }
    }
    iM = which(!is.finite( M[, vnB] ))
    if (length(iM > 0)) {
      M[iM, vnB] = bathymetry_lookup( LOCS=M[ iM, c("lon", "lat")],  lookup_from="core", lookup_to="points" , lookup_from_class="aggregated_data" ) # core=="rawdata"
 
    }

    M = M[ is.finite(M[ , vnB]  ) , ]

    if ( exists("spatial_domain", p)) {
        M = M[ geo_subset( spatial_domain=p$spatial_domain, Z=M ) , ] # need to be careful with extrapolation ...  filter depths
    }

    if ( p$carstm_inputdata_model_source$bathymetry %in% c("stmv", "hybrid") ) {
      pBD = bathymetry_parameters(  spatial_domain=p$spatial_domain, project_class=p$carstm_inputdata_model_source$bathymetry )  # full default
      LU = bathymetry_db( p=pBD, DS="baseline", varnames="all" )
      LU_map = array_map( "xy->1", LU[,c("plon","plat")], gridparams=p$gridparams )
      M_map  = array_map( "xy->1", M[, c("plon","plat")], gridparams=p$gridparams )
      iML = match( M_map, LU_map )
      vns = intersect(  c( "z", "dZ", "ddZ", "b.sdSpatial", "b.sdObs", "b.phi", "b.nu", "b.localrange" ), names(LU) )
      for (vn in setdiff( vns, "z") ) {
        M[, vn] = LU[ iML, vn ]
    }
      M = M[ is.finite( rowSums( M[ , vns])  ) , ]
    }


  
    # --------------------------
   # substrate lookup
    pS = substrate_parameters( p=parameters_reset(p), project_class="carstm"  )
    if (!(exists(pS$variabletomodel, M ))) M[,pS$variabletomodel] = NA
    iM = which(!is.finite( M[, pS$variabletomodel] ))
    if (length(iM > 0)) {
      M[iM, pS$variabletomodel] = substrate_lookup( LOCS=M[iM, c("lon", "lat")], lookup_from="core", lookup_to="points" , lookup_from_class="aggregated_data" ) # core=="rawdata"
    }

    M = M[ is.finite(M[ , pS$variabletomodel]  ) , ]



    # --------------------------
    # temperature observations lookup
    pT = temperature_parameters( p=parameters_reset(p), project_class="carstm", year.assessment=p$year.assessment  )
    if (!(exists(pT$variabletomodel, M ))) M[,pT$variabletomodel] = NA
    iM = which(!is.finite( M[, pT$variabletomodel] ))
    if (length(iM > 0)) {
      M[iM, pT$variabletomodel] = temperature_lookup(  LOCS=M[ iM, c("lon", "lat", "timestamp")],lookup_from="core", lookup_to="points", lookup_from_class="aggregated_data", tz="America/Halifax",
          year.assessment=p$year.assessment
        )
    }


    M = M[ is.finite(M[ , pT$variabletomodel]  ) , ]
    M = M[ which( M[, pT$variabletomodel]  < 14 ) , ]  #


    pPC1 = speciescomposition_parameters( p=parameters_reset(p), project_class="carstm", variabletomodel="pca1" , year.assessment=p$year.assessment)
    if (!(exists(pPC1$variabletomodel, M ))) M[,pPC1$variabletomodel] = NA
    iM = which(!is.finite( M[, pPC1$variabletomodel] ))
    if (length(iM > 0)) {
      M[iM, pPC1$variabletomodel] = speciescomposition_lookup(  LOCS=M[ iM, c("lon", "lat", "timestamp")],lookup_from="core", lookup_to="points", lookup_from_class="aggregated_data", tz="America/Halifax" ,
          year.assessment=p$year.assessment ,
          vnames=pPC1$variabletomodel
        )

    }

    M = M[ which(is.finite(M[, pPC1$variabletomodel] )),]


    pPC2 = speciescomposition_parameters( p=parameters_reset(p), project_class="carstm", variabletomodel="pca2", year.assessment=p$year.assessment )
    if (!(exists(pPC2$variabletomodel, M ))) M[,pPC2$variabletomodel] = NA
    iM = which(!is.finite( M[, pPC2$variabletomodel] ))
    if (length(iM > 0)) {
      M[iM, pPC2$variabletomodel] = speciescomposition_lookup( LOCS=M[ iM, c("lon", "lat", "timestamp")], lookup_from="core", lookup_to="points", lookup_from_class="aggregated_data", tz="America/Halifax" ,
          year.assessment=p$year.assessment,
          vnames=pPC2$variabletomodel
        )

    }
    M = M[ which(is.finite(M[, pPC2$variabletomodel] )),]


    M$plon = NULL
    M$plat = NULL
    M$lon = NULL
    M$lat = NULL
 
    # AUID is character; auid is factor -> numeric 

    M = M[ which(!is.na(M$AUID)),]
    M$AUID = as.character( M$AUID )  # match each datum to an area

    M$tag = "observations"



      # end observations
      # ----------

      # predicted locations (APS)


    region.id = slot( slot(sppoly, "nb"), "region.id" )
    APS = st_drop_geometry(sppoly)

    APS$AUID = as.character( APS$AUID )
    APS$tag ="predictions"
    APS[,p$variabletomodel] = NA
    
    APS[, pB$variabletomodel] = bathymetry_lookup(  LOCS=sppoly, 
      lookup_from = p$carstm_inputdata_model_source$bathymetry,
      lookup_to = "areal_units", 
      vnames="z" 
    ) 
    
    APS[, pS$variabletomodel] = substrate_lookup(  LOCS=sppoly, 
      lookup_from = p$carstm_inputdata_model_source$substrate,
      lookup_to = "areal_units", 
      vnames="substrate.grainsize" 
    ) 

    # to this point APS is static, now add time dynamics (teperature)
    # ---------------------

    vn = c( p$variabletomodel, pB$variabletomodel,  pS$variabletomodel, "tag", "AUID" )
    APS = APS[, vn]

    # expand APS to all time slices
    n_aps = nrow(APS)
    APS = cbind( APS[ rep.int(1:n_aps, p$nt), ], rep.int( p$prediction_ts, rep(n_aps, p$nt )) )
    names(APS) = c(vn, "tiyr")
    APS$year = aegis_floor( APS$tiyr)
    APS$dyear = APS$tiyr - APS$year
    APS$timestamp = lubridate::date_decimal( APS$tiyr, tz=p$timezone )


    APS[, pT$variabletomodel] = temperature_lookup(  LOCS=APS[ , c("AUID", "timestamp")], AU_target=sppoly, 
      lookup_from = p$carstm_inputdata_model_source$temperature,
      lookup_to = "areal_units", 
      vnames_from= paste( pT$variabletomodel, "predicted", sep="."), 
      vnames=pT$variabletomodel ,
      year.assessment=p$year.assessment
    )


    APS[, pPC1$variabletomodel] = speciescomposition_lookup(  LOCS=APS[ , c("AUID", "timestamp")], AU_target=sppoly, 
      lookup_from = p$carstm_inputdata_model_source$speciescomposition,
      lookup_to = "areal_units", 
      vnames_from= paste( pPC1$variabletomodel, "predicted", sep="."), 
      vnames=pPC1$variabletomodel ,
      year.assessment=p$year.assessment
    )


    APS[, pPC2$variabletomodel] = speciescomposition_lookup(  LOCS=APS[ , c("AUID", "timestamp")], AU_target=sppoly, 
      lookup_from = p$carstm_inputdata_model_source$speciescomposition,
      lookup_to = "areal_units", 
      vnames_from= paste( pPC2$variabletomodel, "predicted", sep="."), 
      vnames=pPC2$variabletomodel,
      year.assessment=p$year.assessment
    )




    # useful vars to have for analyses outside of carstm_summary
    varstoadd = c( "totwgt", "totno", "sa", "data_offset",  "zn", "qn" )

    for (vn in varstoadd) if (!exists( vn, APS)) APS[,vn] = NA
    APS$data_offset = 1  # force to solve for unit area


    M = rbind( M[, names(APS)], APS )

    APS = NULL


    M$AUID  = as.character(M$AUID)  # revert to factors
    M$auid  = as.numeric( factor(M$AUID) )
  
    M$auid_main = M$auid

    M$zi  = discretize_data( M[, pB$variabletomodel], p$discretization[[pB$variabletomodel]] )
    M$ti  = discretize_data( M[, pT$variabletomodel], p$discretization[[pT$variabletomodel]] )
    M$gsi = discretize_data( M[, pS$variabletomodel], p$discretization[[pS$variabletomodel]] )

    M$pca1i = discretize_data( M[, pPC1$variabletomodel], p$discretization[[pPC1$variabletomodel]] )
    M$pca2i = discretize_data( M[, pPC2$variabletomodel], p$discretization[[pPC2$variabletomodel]] )

    M$tiyr  = aegis_floor( M$tiyr / p$tres )*p$tres    # discretize for inla .. midpoints

    M$yr = aegis_floor( M$tiyr)

    M$dyear =  M$tiyr - M$yr   # revert dyear to non-discretized form

    M$dyri = discretize_data( M[, "dyear"], p$discretization[["dyear"]] )

    # M$seasonal = (as.numeric(M$year_factor) - 1) * length(p$dyears)  + as.numeric(M$dyear)

    save( M, file=fn, compress=TRUE )

    if (redo) M=fn  # when redo'ing .. return file name aftert the save

    return( M )
  }



  # -------------------------


  if ( DS=="carstm_inputs_hybrid") {

    # various global data sources

    # prediction surface
    crs_lonlat = st_crs(projection_proj4string("lonlat_wgs84"))
    sppoly = areal_units( p=p )  # will redo if not found
    sppoly = st_transform(sppoly, crs=crs_lonlat )
    areal_units_fn = attributes(sppoly)[["areal_units_fn"]]

    # shared accross various secneario using the same polys
    #.. store at the modeldir level as default
    # outputdir = file.path(p$modeldir, p$carstm_model_label)
    outputdir = file.path(p$modeldir )
    if ( !file.exists(outputdir)) dir.create( outputdir, recursive=TRUE, showWarnings=FALSE )

    fn = file.path( outputdir,
      paste( "snowcrab", "carstm_inputs", areal_units_fn,
        p$variabletomodel, paste0(p$selection$survey$data.source, collapse=""),
        p$inputdata_spatial_discretization_planar_km,
        round(p$inputdata_temporal_discretization_yr, 6),
        "rdata",
        sep="."
      )
    )

    if (!redo)  {
      if (file.exists(fn)) {
        load( fn )
        return( M )
      }
    }
    message( "Generating carstm_inputs ... ")


    # do this immediately to reduce storage for sppoly (before adding other variables)
    M = snowcrab.db( p=p, DS="biological_data" )  # will redo if not found .. not used here but used for data matching/lookup in other aegis projects that use bathymetry

    # some survey timestamps extend into January (e.g., 2020) force them to be part of the correct "survey year", i.e., "yr"
    i = which(lubridate::month(M$timestamp)==1)
    if (length(i) > 0) M$timestamp[i] = M$timestamp[i] - lubridate::duration(month=1)

    M$tiyr = lubridate::decimal_date(M$timestamp)

    # reduce size
    M = M[ which( M$lon > p$corners$lon[1] & M$lon < p$corners$lon[2]  & M$lat > p$corners$lat[1] & M$lat < p$corners$lat[2] ), ]
    # levelplot(z.mean~plon+plat, data=M, aspect="iso")

    M$AUID = st_points_in_polygons(
      pts = st_as_sf( M, coords=c("lon","lat"), crs=crs_lonlat ),
      polys = sppoly[, "AUID"],
      varname = "AUID"
    )
    M = M[!is.na(M$AUID),]

    # names(M)[which(names(M)=="yr") ] = "year"
    # M = M[ which(M$year %in% p$yrs), ]
    # M$tiyr = lubridate::decimal_date ( M$timestamp )
    # M$dyear = M$tiyr - M$year


  # --------------------------
    # bathymetry observations  lookup
    pB = bathymetry_parameters( project_class="core"  )
    vnB = pB$variabletomodel
    if ( !(exists(vnB, M ))) {
      vnB2 = paste(vnB, "mean", sep=".")
      if ((exists(vnB2, M ))) {
        names(M)[which(names(M) == vnB2 )] = vnB
      } else {
        M[,vnB] = NA
      }
    }
    iM = which(!is.finite( M[, vnB] ))
    if (length(iM > 0)) {
      M[iM, vnB] = bathymetry_lookup( LOCS=M[ iM, c("lon", "lat")], lookup_from="core", lookup_to="points" , lookup_from_class="aggregated_data" ) # core=="rawdata"
 
    }

    M = M[ is.finite(M[ , vnB]  ) , ]


    if (p$carstm_inputs_aggregated) {
      if ( exists("spatial_domain", p)) {
        M = M[ geo_subset( spatial_domain=p$spatial_domain, Z=M ) , ] # need to be careful with extrapolation ...  filter depths
      }
    }

    if ( exists("spatial_domain", p)) {
        M = M[ geo_subset( spatial_domain=p$spatial_domain, Z=M ) , ] # need to be careful with extrapolation ...  filter depths
    }


    # --------------------------
    # substrate observations  lookup
    pS = substrate_parameters( project_class="core"  )
    if (!(exists(pS$variabletomodel, M ))) M[,pS$variabletomodel] = NA
    iM = which(!is.finite( M[, pS$variabletomodel] ))
    if (length(iM > 0)) {
      M[iM, pS$variabletomodel] = substrate_lookup( LOCS=M[iM, c("lon", "lat")], lookup_from="core", lookup_to="points" , lookup_from_class="aggregated_data" ) # core=="rawdata"
    }

    # M = M[ is.finite(M[ , pS$variabletomodel]  ) , ]



    # --------------------------
    # temperature observations lookup
    pT = temperature_parameters( project_class="core"  )
    if (!(exists(pT$variabletomodel, M ))) M[,pT$variabletomodel] = NA
    iM = which(!is.finite( M[, pT$variabletomodel] ))
    if (length(iM > 0)) {
      M[iM, pT$variabletomodel] = temperature_lookup( LOCS=M[ iM, c("lon", "lat", "timestamp")],lookup_from="core", lookup_to="points", lookup_from_class="aggregated_data", tz="America/Halifax" ,
          year.assessment=p$year.assessment )
    }

    # M = M[ is.finite(M[ , pT$variabletomodel]  ) , ]
    # M = M[ which( M[, pT$variabletomodel]  < 14 ) , ]  #


    # --------------------------

    pPC1 = speciescomposition_parameters(  project_class="core", variabletomodel="pca1" )
    if (!(exists(pPC1$variabletomodel, M ))) M[,pPC1$variabletomodel] = NA
    iM = which(!is.finite(M[, pPC1$variabletomodel]))
    if (length(iM) > 0 ) {
      M[iM, pPC1$variabletomodel] = speciescomposition_lookup( M=M[iM, c("lon", "lat", "timestamp")], sppoly=sppoly, vnames=pPC1$variabletomodel , lookup_from="core", lookup_to="points", lookup_from_class="aggregated_data", tz="America/Halifax" ,
          year.assessment=p$year.assessment ,
          vnames=pPC1$variabletomodel )
    }
    # M = M[ which(is.finite(M[, pPC1$variabletomodel] )),]

    # --------------------------

    pPC2 = speciescomposition_parameters(  project_class="core", variabletomodel="pca2" )
    if (!(exists(pPC2$variabletomodel, M ))) M[,pPC2$variabletomodel] = NA
    iM = which(!is.finite(M[, pPC2$variabletomodel]))
    if (length(iM) > 0 ) {
      M[iM, pPC2$variabletomodel] = speciescomposition_lookup( M=M[iM, c("lon", "lat", "timestamp")], sppoly=sppoly, vnames=pPC2$variabletomodel, lookup_from="core", lookup_to="points", lookup_from_class="aggregated_data", tz="America/Halifax" ,
          year.assessment=p$year.assessment,
          vnames=pPC2$variabletomodel  )
    }
    # M = M[ which(is.finite(M[, pPC2$variabletomodel] )),]

    M$plon = NULL
    M$plat = NULL
    M$lon = NULL
    M$lat = NULL


    M = M[ which(!is.na(M$AUID)),]
    M$AUID = as.character( M$AUID )  # match each datum to an area

    M$tag = "observations"


    # end observations
    # ----------

    # predicted locations (APS)


    region.id = slot( slot(sppoly, "nb"), "region.id" )
    APS = st_drop_geometry(sppoly)

    APS$AUID = as.character( APS$AUID )
    APS$tag ="predictions"
    APS[,p$variabletomodel] = NA


    APS[, pB$variabletomodel] = bathymetry_lookup(  LOCS=sppoly, 
      lookup_from = p$carstm_inputdata_model_source$bathymetry,
      lookup_to = "areal_units", 
      vnames="z" 
    ) 
    
    APS[, pS$variabletomodel] = substrate_lookup(  LOCS=sppoly, 
      lookup_from = p$carstm_inputdata_model_source$substrate,
      lookup_to = "areal_units", 
      vnames="substrate.grainsize" 
    ) 


    # to this point APS is static, now add time dynamics (teperature)
    # ---------------------

    vn = c( p$variabletomodel, pB$variabletomodel,  pS$variabletomodel, "tag", "AUID" )
    APS = APS[, vn]

    # expand APS to all time slices
    n_aps = nrow(APS)
    APS = cbind( APS[ rep.int(1:n_aps, p$nt), ], rep.int( p$prediction_ts, rep(n_aps, p$nt )) )
    names(APS) = c(vn, "tiyr")
    APS$yr = aegis_floor( APS$tiyr)
    APS$dyear = APS$tiyr - APS$yr


    APS[, pT$variabletomodel] = temperature_lookup(  LOCS=APS[ , c("AUID", "timestamp")], AU_target=sppoly, 
      lookup_from = p$carstm_inputdata_model_source$temperature,
      lookup_to = "areal_units", 
      vnames_from= paste(pT$variabletomodel, "predicted", sep="."),
      vnames=pT$variabletomodel,
      year.assessment=p$year.assessment 
    ) 


    APS[, pPC1$variabletomodel] = speiescomposition_lookup(  LOCS=APS[ , c("AUID", "timestamp")], AU_target=sppoly, 
      lookup_from = p$carstm_inputdata_model_source$speciescomposition,
      lookup_to = "areal_units", 
      vnames_from= paste(pPC1$variabletomodel, "predicted", sep="."),
      vnames=pPC1$variabletomodel,
      year.assessment=p$year.assessment 
    ) 


    APS[, pPC2$variabletomodel] = speiescomposition_lookup(  LOCS=APS[ , c("AUID", "timestamp")], AU_target=sppoly, 
      lookup_from = p$carstm_inputdata_model_source$speciescomposition,
      lookup_to = "areal_units", 
      vnames_from= paste(pPC2$variabletomodel, "predicted", sep="."),
      vnames=pPC2$variabletomodel,
      year.assessment=p$year.assessment 
    ) 

    # useful vars to have for analyses outside of carstm_summary
    varstoadd = c( "totwgt", "totno", "sa", "data_offset",  "zn", "qn" )

    for (vn in varstoadd) if (!exists( vn, APS)) APS[,vn] = NA
    APS$data_offset = 1  # force to solve for unit area

    M = rbind( M[, names(APS)], APS )
    APS = NULL


    M$AUID  = as.character(M$AUID)  # revert to factors
    M$auid = match( M$AUID, region.id )

    M$zi  = discretize_data( M[, pB$variabletomodel], p$discretization[[pB$variabletomodel]] )
    M$ti  = discretize_data( M[, pT$variabletomodel], p$discretization[[pT$variabletomodel]] )
    M$gsi = discretize_data( M[, pS$variabletomodel], p$discretization[[pS$variabletomodel]] )

    M$pca1i = discretize_data( M[, pPC1$variabletomodel], p$discretization[[pPC1$variabletomodel]] )
    M$pca2i = discretize_data( M[, pPC2$variabletomodel], p$discretization[[pPC2$variabletomodel]] )

    M$tiyr  = aegis_floor( M$tiyr / p$tres )*p$tres    # discretize for inla .. midpoints

    M$yr = aegis_floor( M$tiyr)
    M$year_factor = as.numeric( factor( M$yr, levels=p$yrs))

    M$dyear =  M$tiyr - M$yr   # revert dyear to non-discretized form

    M$dyri = discretize_data( M[, "dyear"], p$discretization[["dyear"]] )

    # M$seasonal = (as.numeric(M$year_factor) - 1) * length(p$dyears)  + as.numeric(M$dyear)

    save( M, file=fn, compress=TRUE )

    if (redo) M=fn  # when redo'ing .. return file name aftert the save

    return( M )
  }


  # --------------------------


  if ( any( grepl("carstm_output", DS) ) ) {

    # ie. usually run by "carstm_output_compute" which will bring you here

    sppoly = areal_units( p=p )
    areal_units_fn = attributes(sppoly)[["areal_units_fn"]]

    aufns = carstm_filenames( p=p, returntype="carstm_modelled_fit", areal_units_fn=areal_units_fn )
    # same file naming as in carstm ..
    outputdir = dirname( aufns )

    if ( !file.exists(outputdir)) dir.create( outputdir, recursive=TRUE, showWarnings=FALSE )

    fn     =  paste( gsub(".rdata", "", aufns), "aggregated_timeseries", "rdata", sep="." )  
    fn_no  =  paste( gsub(".rdata", "", aufns), "space_timeseries_number", "rdata", sep="." ) 
    fn_bio =  paste( gsub(".rdata", "", aufns), "space_timeseries_biomass", "rdata", sep="." ) 
    fn_pa  =  paste( gsub(".rdata", "", aufns), "space_timeseries_pa", "rdata", sep="." )  

    if ( DS=="carstm_output_timeseries" ) {
      RES = NA
      if (file.exists(fn)) load( fn)
      return( RES )
    }

    if ( DS=="carstm_output_spacetime_number" ) {
      nums = NA
      if (file.exists(fn_no)) load( fn_no )
      return( nums )
    }

    if ( DS=="carstm_output_spacetime_biomass" ) {
      biom = NA
      if (file.exists(fn_bio)) load( fn_bio )
      return( biom )
    }

    if ( DS=="carstm_output_spacetime_pa")  {
      pa = NA
      if (file.exists(fn_pa)) load( fn_pa )
      return( pa )
    }

    # construct meanweights matrix used to convert number to weight

    fit = carstm_model( p=p, DS="carstm_modelled_fit" ) # to load currently saved res
    res = carstm_model( p=p, DS="carstm_modelled_summary"  )

    if (p$carstm_modelengine == "inla") {
      ps = inla.posterior.sample(n=p$nsims, fit, selection=list(Predictor=0))  # only predictions, 0== all predictions
    }

    if (p$carstm_modelengine %in% c( "glm", "gam") ) {
      # sample from marginal distributions as iid assumed
      mu = res[[ paste( p$variabletomodel, "predicted", sep=".")]]
      sigma = res[[ paste( p$variabletomodel, "predicted_se", sep=".")]]
      n = length( c(mu) )
      ncolres = ncol( mu )
      nrowres = nrow( mu )
      ps = tapply( 1:p$nsims, INDEX=1:p$nsims, FUN = function(x) { rnorm( n, mean=c(mu), sd=c(sigma) ) } )
   }


    if (p$selection$type %in% c("presence_absence") ) {
      pa =  res[[ paste( p$variabletomodel, "predicted", sep=".")]]
      pa[!is.finite(pa)] = NA
      pa = inverse.logit(pa)
      pa[!is.finite(pa)] = NA
      # if (is.na(extrapolation_limit)) extrapolation_limit = c(0,1)
      save( pa, file=fn_pa, compress=TRUE )

      sims = sapply( ps,
        function(x) {

          if (p$carstm_modelengine %in% c("glm", "gam") ) {
            pa = matrix(x, nrow=nrowres, ncol=ncolres)
          }

          if (p$carstm_modelengine == "inla") {
            input = x$latent[res$i_preds]
            input[!is.finite(input)] = NA
            input = inverse.logit( input )
            pa = reformat_to_array( input=input , matchfrom=res$matchfrom, matchto=res$matchto )
          }

          pa[!is.finite(pa)] = NA

          o = list()
          o$cfaall    = colSums( pa * sppoly$au_sa_km2/ sum(sppoly$au_sa_km2), na.rm=TRUE )
          o$cfanorth  = colSums( pa * sppoly$cfanorth_surfacearea/ sum(sppoly$cfanorth_surfacearea), na.rm=TRUE )
          o$cfasouth  = colSums( pa * sppoly$cfasouth_surfacearea/ sum(sppoly$cfasouth_surfacearea), na.rm=TRUE )
          o$cfa23     = colSums( pa * sppoly$cfa23_surfacearea/ sum(sppoly$cfa23_surfacearea), na.rm=TRUE )
          o$cfa24     = colSums( pa * sppoly$cfa24_surfacearea/ sum(sppoly$cfa24_surfacearea), na.rm=TRUE )
          o$cfa4x     = colSums( pa * sppoly$cfa4x_surfacearea/ sum(sppoly$cfa4x_surfacearea), na.rm=TRUE )
          return(o)
        }, simplify=TRUE
      )

    }


    if (p$selection$type %in% c("biomass", "number") ) {

      # M = snowcrab.db( p=p, DS="carstm_inputs" )

      M = res$M
      M$yr = M$year  # req for meanweights

      wgts = meanweights_by_arealunit(
        set=M[M$tag=="observations",],
        AUID=as.character( sppoly$AUID ),
        yrs=p$yrs,
        fillall=TRUE,
        annual_breakdown=TRUE,
        robustify_quantiles=p$quantile_bounds  # high upper bounds are more dangerous
      )

      if (p$selection$type == "biomass") {
        biom = res[[ paste( p$variabletomodel, "predicted", sep=".")]]
        biom[!is.finite(biom)] = NA

        NA_mask = NULL
        nnn = which( !is.finite(biom ))
        if (length(nnn)>0 ) NA_mask = nnn

        if (is.na(extrapolation_limit)) extrapolation_limit = quantile( M$totwgt/M$data_offset, probs=p$quantile_bounds[2], na.rm=T) # 28921.8426

        uu = which( biom > extrapolation_limit )
        if (length(uu) > 0 ) {
          if (is.character(extrapolation_replacement)) if (extrapolation_replacement=="extrapolation_limit" ) extrapolation_replacement = extrapolation_limit
          biom[ uu] = extrapolation_replacement
          warning("\n Extreme-valued predictions were found, capping them to max observed rates .. \n you might want to have more informed priors, or otherwise set extrapolation_replacement=NA to replacement value \n")
        }
        biom[biom > extrapolation_limit ] = extrapolation_limit
        
        biom = biom / 10^6  # kg / km^2 -> kt / km^2

        nums = biom / wgts   # (n * 10^6) / km^2
        save( biom, file=fn_bio, compress=TRUE )
        save( nums, file=fn_no, compress=TRUE )


        sims = sapply( ps,
          function(x) {
            if (p$carstm_modelengine %in% c("glm", "gam") ) {
              biom = matrix(x, nrow=nrowres, ncol=ncolres)
            }

            if (p$carstm_modelengine == "inla") {
              input = exp( x$latent[res$i_preds])
              biom = reformat_to_array( input=input , matchfrom=res$matchfrom, matchto=res$matchto )
              if (!is.null(NA_mask)) biom[NA_mask] = NA
            }

            biom[!is.finite(biom)] = NA
            biom = biom / 10^6  # kg / km^2 -> kt / km^2
            
            o = list()
            o$cfaall    = colSums( biom * sppoly$au_sa_km2, na.rm=TRUE )
            o$cfanorth  = colSums( biom * sppoly$cfanorth_surfacearea, na.rm=TRUE )
            o$cfasouth  = colSums( biom * sppoly$cfasouth_surfacearea, na.rm=TRUE )
            o$cfa23     = colSums( biom * sppoly$cfa23_surfacearea, na.rm=TRUE )
            o$cfa24     = colSums( biom * sppoly$cfa24_surfacearea, na.rm=TRUE )
            o$cfa4x     = colSums( biom * sppoly$cfa4x_surfacearea, na.rm=TRUE )
            return(o)
          }, simplify=TRUE
        )
      }


      if (p$selection$type == "number") {

        nums = res[[ paste( p$variabletomodel, "predicted", sep=".")]]
        nums[!is.finite(nums)] = NA
        NA_mask = NULL
        nnn = which( !is.finite(nums ))
        if (length(nnn)>0 ) NA_mask = nnn

        if (is.na(extrapolation_limit)) extrapolation_limit = quantile( M$totno/M$data_offset, probs=p$quantile_bounds[2], na.rm=T) # 28921.8426
        uu = which( nums > extrapolation_limit )
        if (length(uu) > 0 ) {
          if (is.character(extrapolation_replacement)) if (extrapolation_replacement=="extrapolation_limit" ) extrapolation_replacement = extrapolation_limit
          nums[ uu] = extrapolation_replacement
          warning("\n Extreme-valued predictions were found, capping them to max observed rates .. \n you might want to have more informed priors, or otherwise set extrapolation=NA to replacement value \n")
        }
        nums[nums > extrapolation_limit] = extrapolation_limit
        
        biom = nums * wgts / 10^6  # kg / km^2 -> kt / km^2
        nums = nums / 10^6  # n * 10^6 / km^2

        save( biom, file=fn_bio, compress=TRUE )
        save( nums, file=fn_no, compress=TRUE )

        sims = sapply( ps,
          function(x) {

            if (p$carstm_modelengine %in% c("glm", "gam") ) {
              nums = matrix(x, nrow=nrowres, ncol=ncolres  )
            }

            if (p$carstm_modelengine == "inla") {
              input = exp( x$latent[res$i_preds])
              nums = reformat_to_array( input=input , matchfrom=res$matchfrom, matchto=res$matchto )
              if (!is.null(NA_mask)) nums[NA_mask] = NA
            }

            nums[!is.finite(nums)] = NA
            biom = nums * wgts / 10^6   # kg / km^2 -> kt / km^2

            o = list()
            o$cfaall    = colSums( biom * sppoly$au_sa_km2, na.rm=TRUE )
            o$cfanorth  = colSums( biom * sppoly$cfanorth_surfacearea, na.rm=TRUE )
            o$cfasouth  = colSums( biom * sppoly$cfasouth_surfacearea, na.rm=TRUE )
            o$cfa23     = colSums( biom * sppoly$cfa23_surfacearea, na.rm=TRUE )
            o$cfa24     = colSums( biom * sppoly$cfa24_surfacearea, na.rm=TRUE )
            o$cfa4x     = colSums( biom * sppoly$cfa4x_surfacearea, na.rm=TRUE )
            return(o)
          }, simplify=TRUE
        )

      }
    }

    RES = data.frame( yrs = p$yrs )
    RES$cfaall = apply( simplify2array(sims["cfaall",]), 1, mean )
    RES$cfaall_sd = apply( simplify2array(sims["cfaall",]), 1, sd )
    RES$cfaall_median = apply( simplify2array(sims["cfaall",]), 1, median )
    RES$cfaall_lb = apply( simplify2array(sims["cfaall",]), 1, quantile, probs=0.025 )
    RES$cfaall_ub = apply( simplify2array(sims["cfaall",]), 1, quantile, probs=0.975 )

    RES$cfanorth = apply( simplify2array(sims["cfanorth",]), 1, mean )
    RES$cfanorth_sd = apply( simplify2array(sims["cfanorth",]), 1, sd )
    RES$cfanorth_median = apply( simplify2array(sims["cfanorth",]), 1, median )
    RES$cfanorth_lb = apply( simplify2array(sims["cfanorth",]), 1, quantile, probs=0.025 )
    RES$cfanorth_ub = apply( simplify2array(sims["cfanorth",]), 1, quantile, probs=0.975 )

    RES$cfasouth = apply( simplify2array(sims["cfasouth",]), 1, mean )
    RES$cfasouth_sd = apply( simplify2array(sims["cfasouth",]), 1, sd )
    RES$cfasouth_median = apply( simplify2array(sims["cfasouth",]), 1, median )
    RES$cfasouth_lb = apply( simplify2array(sims["cfasouth",]), 1, quantile, probs=0.025 )
    RES$cfasouth_ub = apply( simplify2array(sims["cfasouth",]), 1, quantile, probs=0.975 )

    RES$cfa23 = apply( simplify2array(sims["cfa23",]), 1, mean )
    RES$cfa23_sd = apply( simplify2array(sims["cfa23",]), 1, sd )
    RES$cfa23_median = apply( simplify2array(sims["cfa23",]), 1, median )
    RES$cfa23_lb = apply( simplify2array(sims["cfa23",]), 1, quantile, probs=0.025 )
    RES$cfa23_ub = apply( simplify2array(sims["cfa23",]), 1, quantile, probs=0.975 )

    RES$cfa24 = apply( simplify2array(sims["cfa24",]), 1, mean )
    RES$cfa24_sd = apply( simplify2array(sims["cfa24",]), 1, sd )
    RES$cfa24_median = apply( simplify2array(sims["cfa24",]), 1, median )
    RES$cfa24_lb = apply( simplify2array(sims["cfa24",]), 1, quantile, probs=0.025 )
    RES$cfa24_ub = apply( simplify2array(sims["cfa24",]), 1, quantile, probs=0.975 )

    RES$cfa4x = apply( simplify2array(sims["cfa4x",]), 1, mean )
    RES$cfa4x_sd = apply( simplify2array(sims["cfa4x",]), 1, sd )
    RES$cfa4x_median = apply( simplify2array(sims["cfa4x",]), 1, median )
    RES$cfa4x_lb = apply( simplify2array(sims["cfa4x",]), 1, quantile, probs=0.025 )
    RES$cfa4x_ub = apply( simplify2array(sims["cfa4x",]), 1, quantile, probs=0.975 )

    save( RES, file=fn, compress=TRUE )

    # repeat file extraction here  in case computation was required
    if ( DS=="carstm_output_timeseries" ) {
      RES = NA
      if (file.exists(fn)) load( fn)
      return( RES )
    }

    if ( DS=="carstm_output_spacetime_number" ) {
      nums = NA
      if (file.exists(fn_no)) load( fn_no )
      return( nums )
    }

    if ( DS=="carstm_output_spacetime_biomass" ) {
      biom = NA
      if (file.exists(fn_bio)) load( fn_bio )
      return( biom )
    }

    if ( DS=="carstm_output_spacetime_pa")  {
      pa = NA
      if (file.exists(fn_pa)) load( fn_pa )
      return( pa )
    }

    return(fn)
  }


}  ## end snowcrab.db
