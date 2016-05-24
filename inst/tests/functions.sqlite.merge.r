# 
#
# This is an initial attempt at interfacing with an SQLITE backend.
# Stored for posterity.
#
#
#-------------------------------------------------------------------------------------------------------------------------

# Brent Cameron 
# R code Developed as part of a service contract for Fisheries and Oceans Canada  
# Created and tested with R version 2.10.1 in Microsoft Windows Vista Business SP2

# REFERENCES
# * R code provided by Jae Choi found in ecomod/common and ecomod/snowcrab/src directories  
# * Tutorials and ref card found at http://cran.r-project.org/
# * Tutorials found at http://www.sqlite.org/
# * Program used for quick database views http://sqliteman.com/

# CONTACT
# Brent Cameron
# 107 Gibbons Rd
# Huntington NS
# B1K 1V1
# Telephone: 1 902 727 2157
# Email: brentcameron1@gmail.com

# ---------------------------------------------------------------------------------------------------------------------------


  # Function that combines the netmind and minilog bottom and stats data with 
  # the trawl data. Running this will create or overwrite a table entry in the 
  # set database called merged_yr where yr is a desired year. This function first 
  # merges the minilog and netmind data on station, unique_ids and timestamps within 
  # 10 min. The result of this is merged with the trawl data on trip, set num and 
  # occurance time within 10 minuits of each other. Failure of any merge attempt will 
  # be printed. This function also formats the merged table to have the proper field 
  # names in the proper places. NOTE: Problem in trawl data stationid, exa. stationid
  # of 1012 instead of 101r2 code found in ncleandbs.r lines 236 - 243 can resolve this   

mergeTrawl = function(path, yearDesired, proj.type){
	
	print("            MERGE NETMIND AND MINILOG DATA WITH TRAWL/SET DATA                   ") 
	
	# Open connections to netmind and minilog databases 
	minicon <- dbConnect(dbDriver("SQLite"), "minilog.db")
	netcon  <- dbConnect(dbDriver("SQLite"), "netmind.db")

    # Get table names from these databases
	res <- dbSendQuery(minicon, "SELECT name FROM sqlite_master WHERE type='table' ORDER BY name")
	minitables <- fetch(res)
	dbClearResult(res)

	res <- dbSendQuery(netcon, "SELECT name FROM sqlite_master WHERE type='table' ORDER BY name")
	nettables <- fetch(res)
	dbClearResult(res)

	minitables = minitables[,1]
	nettables = nettables[,1]

	ntables = length(minitables)
	mtables = length(nettables)

	# Loop that loops through each table and looks at the desired tables 
	for(i in 1:ntables){ 
		yr = substr(minitables[i], 16, 19)
		new = NULL
		if(yr %in% yearDesired || length(yearDesired) == 0)	
		if(substr(minitables[i], 1, 15) == "minilog_bottom_"){
			if(!paste("netmind_stats_", yr, sep = "") %in% nettables){
				print(paste("No netmind stats table for", yr,  sep = " "))
				next()
			}
            
			#Get data from desired tables
			minilog.bottom = dbReadTable(minicon, minitables[i])
			netmind.stats = dbReadTable(netcon, paste("netmind_stats_", yr, sep = ""))
			
			# OPTIONAL Remove the stations that have been redone and should not have trawl data for anyway, 
				netmind.stats = removeredos(netmind.stats)
				minilog.bottom = removeredos(minilog.bottom)
			
			# Remove or rename field entries 
			netmind.stats$t0_n = netmind.stats$t0
			netmind.stats$t0 = NULL
			netmind.stats$t1_n = netmind.stats$t1
			netmind.stats$t1 = NULL
			netmind.stats$dt_n = netmind.stats$dt
			netmind.stats$dt = NULL
			netmind.stats$stationid_netconfig = netmind.stats$stationid
			netmind.stats$stationid = NULL
			netmind.stats$filename_n = netmind.stats$filename
			netmind.stats$filename = NULL
			netmind.stats$unique_idn = netmind.stats$unique_id
			netmind.stats$unique_id = NULL
			netmind.stats$trip_netconfig = netmind.stats$trip
			netmind.stats$trip = NULL
			netmind.stats$yr_n = netmind.stats$yr
			netmind.stats$yr = NULL
			
			print(paste("Merging Netmind and Minilog bottom data for", yr,  sep = " "))
			
			# Stop if unequal numer of stations to merge. Indicates a problem with the previous 
			# data entry
			if(nrow(minilog.bottom) != nrow(netmind.stats)){
				print(paste("Aborting, please run cleandbs. Unequal number of good stations for ", yr))
				stop() 
			}

			# Loop through each station and match on station number, unique_ids and occurance time
			for(i in 1:nrow(minilog.bottom)){
				same = which(minilog.bottom$stationid[i] == netmind.stats$stationid_n)
				uninum = gsub(".minilog", "", minilog.bottom$unique_id[i], ignore.case=T)
				
				same2 = which(grepl(uninum, netmind.stats$unique_idn))
				same3 = which(abs(as.numeric(        string2chron(minilog.bottom$timestamp[i]) - string2chron(netmind.stats$netmind_timestamp)    ))  < 20/60/24) #10 min
				same = append(same, same2)
				same = append(same, same3)
				
				same = unique(same[duplicated(same)])
				
				# No match, problem exists.
				if(length(same) == 0){
					print(paste("No matching stats file for" , minilog.bottom$stationid[i], ", try running cleandbs function.", sep = " "))
					next()
				}
				#More than one match, problem exists
				if(length(same) > 1) {
					print(paste("Cannot determine correct stats file match for" , minilog.bottom$stationid[i], "between files:", netmind.stats$stationid_netconfig[same],  sep = " "))
					next()
				}
				
				#add the entry and continue to the next station
				new = rbind(new , cbind(minilog.bottom[i ,] , netmind.stats[same ,]))

				
			}

            #Update touchdown, liftoff and down time fields with netmind data 
			new$t0 = new$t0_n
			new$t0_n = NULL 
			new$t1 = new$t1_n
			new$t1_n = NULL 
			new$dt = new$dt_n
			new$dt_n = NULL 
			
			# Duplicate yr feilds, remove one
			new$yr = NULL
			
			#Check to see if set database has been created if not then create
			setcon <- dbConnect(dbDriver("SQLite"), "set.db")
			if(!dbExistsTable(setcon, "set_date_time")){
				RTrawldata2db(path)
				setdatetime()
				
			}
		    
			# Get trawl set information
			set = dbReadTable(setcon, "set_date_time")
			dbDisconnect(setcon)
			
           
			set$chron = as.chron(as.numeric(set$chron))
			out = NULL
			
			# Get data from the trawl data where the years are the same, The fullupdate variable 
			# will be nessesary in order to show which stations from the trawl data have not yet been merged 
			fullupdate = which(yr == set$yr)
			fullupdate = cbind(set$trip[fullupdate], set$set[fullupdate], set$station[fullupdate])
			fullupdate = as.data.frame(fullupdate)	
			names(fullupdate) = c("trip", "set", "station")

			print(paste("Merging bottom data with trawl/set data for", yr,  sep = " "))

            # For loop that goes through each netmind/minilog entry and merges with the correct trawl entry.
			# The merge is made on trip and set numbers along with the time of occurance 
			for(i in 1:nrow(new)){
				
				same = which(new$trip[i] == set$trip  )
				same2 = which(new$setno[i] == set$set )
				same3 = which( abs(as.numeric(  as.chron(as.numeric(new$t0[i])) - set$chron )) < 20/60/24 )
				same = append(same, same2)
				same = append(same, same3)
				
				same = unique(same[duplicated(same)])
				
				# Move without merging if no match found
				if(length(same) == 0){
					print("---------------------------------------------------------------------------")
					print(paste("No matching Trawl file for" , new$filename[i], ", there may not have had a good tow for this station or the station has not yet been entered into the trawl database.", sep = " "))
					print("---------------------------------------------------------------------------")	
					next()
				}
				# Move on without merging if more than one match found 
				if(length(same) > 1) {
					print("---------------------------------------------------------------------------")
					print(paste("Cannot determine correct trawl file match for" , new$filename[i], "between files:", set$station[same],  sep = " "))
					print("---------------------------------------------------------------------------")
					next()
				}
				
				
				# add the new entry
				out = rbind(out , cbind(new[i ,] , set[same ,]))
				
				
			}
			
			# Call to fix the fields of minilog data similar to how they were fixed in the past they have
			out = fixminifields(out)
			# Call to fix the netmind feilds similar to how they were in the past
			out = fixnetfields(out, proj.type)
			
			# Call to check of additional errors that the user should be made aware of and also 
			# clean up the feilds so that it is of the same format as in the past
			out = checkandclean(out)
			
			#for loop that removes matched trawl data from fullupdate list, the remaining stations 
            #if any will be printed			
			for(j in 1:nrow(out)){
				pos = which(out$trip[j] == fullupdate$trip)
				pos2 = which(out$set[j] == fullupdate$set)
				pos = append(pos, pos2)
				pos = unique(pos[duplicated(pos)])
				fullupdate = fullupdate[- pos,]
				
			}
			if(nrow(fullupdate) == 0){
				print("---------------------------------------------------------------------------")
				print("All trawl data matched for this year")
				print("---------------------------------------------------------------------------")
			}
			else{
				print("---------------------------------------------------------------------------")
				print("Trawldata for which no match was found:")
				print(fullupdate)
				print("---------------------------------------------------------------------------")
			}
			#Write merged table to the set database
			if(!is.null(out)){
				setcon  <- dbConnect(dbDriver("SQLite"), "set.db")
				dbWriteTable(setcon, paste("merged_", yr, sep="") , out, overwrite = T, row.names=F)
				dbDisconnect(setcon)
				print("---------------------------------------------------------------------------")
				print(paste( nrow(out), " entries were written to Table merged_", yr, " in the set database.", sep = ""))  
				print("---------------------------------------------------------------------------")
			}
			
		}

	}

}


# Function that is called to remove the stations that have ben redone and should have 
# no trawl data match. If it is found that some trawl stations are not matching with the  
# minilog/netmind data then try commenting these calls out to see if a match is found from
# a deemed bad station 
removeredos = function(x){
	stations = NULL
	if(is.null(x$stationid)){
		x$stationid = x$station
		x$station = NULL
	}

	stations = x$stationid
	stations = gsub("[[:alpha:]]*", "", stations, ignore.case=T)
	stations = gsub("[[:punct:]]*", "", stations, ignore.case=T)
	stations = unique(stations)
    
	for(i in 1:length(stations)){

		stat = paste("ep", stations[i],"r3", sep = "")
		if( stat %in% x$stationid){
			x = x[- which(paste("ep", stations[i], "r2", sep = "") == x$stationid ) ,]        
			x = x[- which(paste("ep", stations[i], "r", sep = "") == x$stationid ) ,]
			x = x[- which(paste("ep", stations[i], sep = "") == x$stationid ) ,]  
			
		}
		stat = paste("ep", stations[i],"r2", sep = "")
		if( stat %in% x$stationid){
			x = x[- which(paste("ep", stations[i], "r", sep = "") == x$stationid ) ,]
			x = x[- which(paste("ep", stations[i], sep = "") == x$stationid ) ,]  
			
		}
		stat = paste("ep", stations[i],"r", sep = "")
		if( stat %in% x$stationid){
			x = x[- which(paste("ep", stations[i], sep = "") == x$stationid ) ,]  
			
		}



	}

	return(x)

}

# Function that is called to fix minilog fields according to methods used in the past
fixminifields = function(set){

	
	set$minilog_uid = as.character( set$unique_id )
	set$unique_id = NULL
	set$surveytype = NULL
	set$SN = NULL
	set$studyid = NULL
	set$minilog_n = set$n
	set$n = NULL
	set$stationid = gsub("[[:alpha:]]*", "", set$stationid, ignore.case=T)
	set$stationid = gsub("[[:punct:]]*", "", set$stationid, ignore.case=T)
	set$stationid = as.numeric( set$stationid )
	
	set$timestamp = as.character(set$timestamp)
	set$timestamp = string2chron(set$timestamp)
    set$t1 = as.chron(as.numeric(set$t1))
	set$t0 = as.chron(as.numeric(set$t0))
	set$startdate = convert.datecodes(set$t0, "tripcode")
	
	
	set$sn = NULL
	set$study_id = NULL

	# merge in historical temp and depth records when no data are obtained from minilogs
	q = which( !is.finite(set$t) )
	if ( length (q) > 0 ) {
		set$t[q] = set$Tx[q]
	}

	q = which(!is.finite( set$z) )
	if ( length (q) > 0 ) {
		set$z[q] = set$Zx[q]
	}

	set$Zx = NULL
	set$Tx = NULL
	
	#Can add additional checks here if building databases from historical data
	
	return(set)
}

# Function that is called to fix netmind fields and provide error checks according to methods used in the past
fixnetfields = function(set, proj.type){

	set$netmind_uid = as.character(set$unique_idn)
	set$unique_idn = NULL
	
	set$trip = as.character( set$trip ) 
	set$netmind_timestamp = string2chron(set$netmind_timestamp)
	

	iii = which( is.finite(set$surfacearea) )
	set$sa[iii] = set$surfacearea[iii]

	
	time.diff = set$netmind_timestamp - set$chron
	time.thresh = 30/60/24  # in days
	i = which( abs( as.numeric( time.diff)) > time.thresh )  
	if(length(i)>0){
		print ("Potential date/time mismatches::")
		print( set[i, ] )
	}
	## check positional information
	
	k = which( abs(set$slon - set$lon) > 0.2 )
	if(length(k)>0){ 
		
		print("---------------------------------------------------------------------------")
		print( "The following have discrepancies in longitude (SNCARSETS vs NETMIND): ")
		print( "They will be overridden by NETMIND info as they are directly entered from GPS entries: ")
		print( set[k , c("trip", "set", "station", "lon", "slon")])
		print("---------------------------------------------------------------------------")	 
		set$lon[k] = set$slon[k]
		set$lat[k] = set$slat[k]
	}
	
	l = which( abs(set$slat - set$lat) > 0.2 )
	if(length(l)>0){
		print("---------------------------------------------------------------------------")	 
		print( "The following have discrepancies in longitude (SNCARSETS vs NETMIND): ")
		print( "They will be overridden by NETMIND info as they are directly entered from GPS entries: ")
		print( set[l , c("trip", "set", "station", "lat", "slat")])
		print("---------------------------------------------------------------------------")	      
	}
	
	j = which( is.finite( set$slon) & is.finite(set$slat) ) # positional data obtained directly from Netmind GPS and Minilog T0
	set$lon[j] = set$slon[j]
	set$lat[j] = set$slat[j]

	set = lonlat2planar(set, proj.type=proj.type) # get planar projections of lon/lat in km
	
	set$slon = NULL
	set$slat = NULL

	return(set)
}

# Function that does additional error checks and further formats the merged table to a desired format
checkandclean = function(set){


	set$date.check = as.numeric( set$timestamp - set$t0 ) * 24 *60
	e0 = which( (abs( set$date.check) > 7)  & 
	( as.numeric(as.character(years(set$t0)))>=2004 ) )
	set$dt = as.chron(as.numeric(set$dt))
	if  (length(e0)>0 ) {
		print("---------------------------------------------------------------------------")	 
		print( "The following have rather long tow times (dt):" )
		print( set[e0, c("stationid", "yr", "t0", "t1", "dt", "date.check")] )
		print("---------------------------------------------------------------------------")	 
	}
	set$date.check = NULL

	pos = which(set$yr != set$yr_n)
	if(length(pos)>0){
		print("---------------------------------------------------------------------------")	 
		print("Merge failure! Year mismatches at:")
		print(set$filename[pos])
		print("---------------------------------------------------------------------------")	  
		stop()
	}
	set$yr_n = NULL


	pos = which(set$filename != set$filename_n)
	if(length(pos)>0){
		print("---------------------------------------------------------------------------")	  
		print("Filename mismatches at:")
		print(paste(set$filename[pos], " filenames did not match with ", set$filename_n[pos], sep = ""))
		print("---------------------------------------------------------------------------")	  
	}

	pos = which(set$set != set$setno)
	if(length(pos)>0){
		print("---------------------------------------------------------------------------")	 
		print("set number mismatches at:")
		print(paste(set$filename[pos], " set numbers did not match with ", set$filename_n[pos], sep = ""))
		print("---------------------------------------------------------------------------")	 
	}
    nset = NULL

	ord = match(c('netmind_uid', 'minilog_uid', "trip", "set", "station", "observer", "cfa", "lon", "lat", "towquality", "gear", "sa", "netmind", "minilog", "chron", "julian", "yr", "z", "t", "zsd", "tsd", "minilog_n", "t0", "t1", "timestamp", "stationid", "startdate", "dt", "trip_netconfig", "setno", "stationid_netconfig", "filename", "netmind_timestamp", "distance", "spread", "spread_sd", "surfacearea", "vel", "vel_sd", "netmind_n", "plon", "plat"), names(set))
    ord = na.omit(ord)
	set <- set[,unlist(ord)] 

	return(set)
}




# This function is called in order to maintain consistency between netmind and minilog databases.
# Calling this function renames the station names so that they follow the same structure between years. The 
# redone stations are renamed based on their timestamp. This ensures that the latest station towed gets 
# the higest r# value. This function also maintains consistency by doing two things if a station is in one 
# database and not in the other. First the user will be told what station is needed in which database. 
# Second it will remove the station so that both contain the same stations. If a station
# is removed then the unique ids are updated in order to keep further additions unique.
# Updating unique ids is timely so running this function should only be done periodically
cleandbs = function(yearDesired){
	print("             CLEANING NETMIND AND MINILOG DATABASES             ")
	#Open connections to netmind aand minilog databases
	minicon = dbConnect(dbDriver("SQLite"), "minilog.db")
	netcon  = dbConnect(dbDriver("SQLite"), "netmind.db")

    #Get all table names from these databases
	res <- dbSendQuery(minicon, "SELECT name FROM sqlite_master WHERE type='table' ORDER BY name")
	minitables <- fetch(res)
	dbClearResult(res)

	res <- dbSendQuery(netcon, "SELECT name FROM sqlite_master WHERE type='table' ORDER BY name")
	nettables <- fetch(res)
	dbClearResult(res)

	minitables = minitables[,1]
	nettables = nettables[,1]

	ntables = length(minitables)
	mtables = length(nettables)
    #Loop through each table and pick out desired tables
	for(i in 1:ntables){ 
		yr = substr(minitables[i], 18, 21)
		if(yr %in% yearDesired || length(yearDesired) == 0)	
		if(substr(minitables[i], 1, 17) == "minilog_metadata_"){
			#Only continue if netmind table also exists to match the minilog table
			if(!paste("netmind_metadata_", yr, sep = "") %in% nettables){
				print(paste("netmind database contains no data for", yr, sep = " "))
				next()
			}
			print(yr)
			# Read the tables to dataframes, convert the timestamp fields to a chron objects
			# and order the meta frames by when the tow occured. This allows the redone stations
			# to be named appropriately
			minilog.metadata = dbReadTable(minicon, minitables[i])
			minilog.metadata$timestamp = string2chron(minilog.metadata$timestamp)
			minilog.metadata = minilog.metadata[order(minilog.metadata$timestamp),] 	
			
			minilog.basedata = dbReadTable(minicon, paste("minilog_basedata_", yr, sep =""))
			
			netmind.metadata = dbReadTable(netcon, paste("netmind_metadata_", yr, sep = ""))
			netmind.metadata$netmind_timestamp = string2chron(netmind.metadata$netmind_timestamp)
			netmind.metadata = netmind.metadata[order(netmind.metadata$netmind_timestamp),]
		
			netmind.basedata = dbReadTable(netcon, paste("netmind_basedata_", yr, sep =""))
			
			netmind.metadata$station = fileform(netmind.metadata$station)
            minilog.metadata$stationid = fileform(minilog.metadata$stationid)

            #Variables to keep track of which stations are in which databases			
			needNet = c()
			haveMini = c()
			needMini = c()
			netmis = ""
			#Loop through each station of minilog metadata checking its existance in the netmind metadata
			#and preforming the nessesery tasks
			j=1
			while(!is.na(minilog.metadata$stationid[j])){ 
				if(!minilog.metadata$stationid[j] %in%  netmind.metadata$station){
					#Add missing file to a list of missing netmind files
					netmis = minilog.metadata$filename[j]
					needNet = append(needNet, netmis)
					#Remove basedata for this station as it does not exist in the netmind database
					minilog.basedata = subset(minilog.basedata , minilog.metadata$unique_id[j]!= minilog.basedata$unique_id)
					#Remove metadata for this station as it does not exist in the netmind database
					minilog.metadata=minilog.metadata[-j,]
					
					print(paste("Please add file", netmis , "for", yr, "as it was not in the netmind database. This file has been removed from the minilog database.", sep = " ")) 
					j=j-1
				}
				#Add station to haveMini list. 
				else{ haveMini = append(haveMini, minilog.metadata$stationid[j]) }
				j=j+1
			}
			# if(length(needNet) > 0){
				# print(paste("Files needed in netmind.db for ", yr, " :", sep = ""))
				# print(needNet)
			#}
			
			if(nrow(minilog.basedata) < 1){
				print("ABORTING. no similar files found, check data")
				next()
			}
			
			#Loop through each netmind station to check if it is also in minilog database
			j=1
			while(!is.na(netmind.metadata$station[j])){ 
				if(!netmind.metadata$station[j] %in%  haveMini){
					minimis = netmind.metadata$filename[j]
					#needMini = append(needMini, minimis)
					# Remove basedata for this station as this station does not exist in the minilog database
					netmind.basedata = subset(netmind.basedata , netmind.metadata$unique_id[j]!= netmind.basedata$unique_id)
					# Remove metadata for this station as this station does not exist in the minilog database
					netmind.metadata=netmind.metadata[-j,]
					
					print(paste("Please add file", minimis , "for", yr, "as it was not in the minilog database. This file has been removed from the netmind database.", sep = " ")) 
					j=j-1
				}
				j=j+1
			}
			# if(length(needMini) > 0){
				# print(paste("Files needed in minilog.db for ", yr, " :", sep = ""))
				# print(needMini)
			# }
			
			
			newbase = NULL
			newbase2 = NULL
			# Incase update did not properly work however this should not happen
			if(nrow(minilog.metadata) != nrow(netmind.metadata)){
				print("FILES MISSING, cannot proceed to update unique ids for this year")
				next()
			}
			
			# The following tests if the unique ids need to be remade. If files have been removed
			# above then renaming will be nesessary else no rename done as it is time and resource
			# consuming
			hole = F
			#Get numeric part of id 
			sequ = gsub("[[:alpha:]]*", "", minilog.metadata$unique_id, ignore.case=T)
			sequ = gsub("[[:punct:]]*", "", sequ, ignore.case=T)
			sequ = substr(sequ, 5, length(sequ))
			sequ = as.numeric(sequ)
			#Test if numeric part of id is in sequence
			for(m in 1:length(sequ))
			if(sequ[m] != m) hole = T
			#Update the ids if they are not in sequence
			if(hole == T){
			
            uni = c()
			uni2 = c()
			newuni = c()
			newuni2 = c()
            # Loop through each row and assign a new unique id for basedata for both netmind and 
			# minilog entries
			for(h in 1:nrow(minilog.metadata)){
				#Make a unique id list for easy update of metadata ids
				uni = append(uni, paste(yr, h , "minilog", sep = "."))
				uni2 = append(uni2, paste(yr, h , "netmind", sep = "."))
				
				#Pick out minilog basedata for the station and change the ids
				wanted =  subset(minilog.basedata , minilog.metadata$unique_id[h]== minilog.basedata$unique_id)
				wanted$unique_id = paste(yr, h , "minilog", sep = ".")
				#Bind the new data to a new frame
				newbase = rbind(newbase, wanted)

				#Pick out minilog basedata for the station and change the ids
			    wanted =  subset(netmind.basedata , netmind.metadata$unique_id[h]== netmind.basedata$unique_id)
				wanted$unique_id = paste(yr, h , "netmind", sep = ".")
				#Bind the new data to a frame
				newbase2 = rbind(newbase2, wanted)
				
			}
			#Update the metadata uniquje ids
			minilog.metadata$unique_id = uni
			netmind.metadata$unique_id = uni2
			minilog.basedata = newbase
			netmind.basedata = newbase2
			}
			
			minilog.metadata$timestamp = as.character(minilog.metadata$timestamp)
			netmind.metadata$netmind_timestamp = as.character(netmind.metadata$netmind_timestamp)
			
			#Write updated tables to the proper databases
			
			dbWriteTable(minicon, paste("minilog_basedata", yr, sep = "_"), minilog.basedata, overwrite = T, row.names=F)
			dbWriteTable(minicon, paste("minilog_metadata", yr, sep = "_"), minilog.metadata, overwrite = T, row.names=F)
			rm(minilog.metadata, minilog.basedata, newbase)

			dbWriteTable(netcon, paste("netmind_basedata", yr, sep = "_"), netmind.basedata, overwrite = T, row.names=F)
			dbWriteTable(netcon, paste("netmind_metadata", yr, sep = "_"), netmind.metadata, overwrite = T, row.names=F)
			rm(netmind.metadata, netmind.basedata, newbase2)

		}
		
	}

	dbDisconnect(minicon)		
	dbDisconnect(netcon)		

	#makebottomTables(yearDesired)
	#makestatsTables(yearDesired)	
}





# Function returns a formated list station names

fileform = function(filenames) {
	filenames = unlist(filenames)
	h = 1
	while(!is.na(filenames[h])){
		station = filenames[h]
		#Removes ep from station name
		station = gsub("ep", "", station)
		station = unlist(strsplit(station, ""))
		
		i = 1
		stationnum = ""
		#Loop that gets the station number out of the station name works better
		#than just removing everything but numbers since some stations are named 
		#such things as 433bad2 or 433r2 resulting in station to be 4332 instead
		#of a more accurate 433
		while(!is.na(as.numeric(station[i]))){
			stationnum = paste(stationnum, station[i], sep = "")
			i = i+1
		}
		#Some station have leading zero this handels this case by removing the leading zero
		if(nchar(stationnum)>3){
			
			filename2 = strsplit(stationnum, "")
			
			if(unlist(filename2)[1] == 0){
				old = filenames[h]
				filenames[h] = paste( substr(filenames[h], 1, 2), substr(filenames[h], 4, nchar(filenames[h])))
				filenames[h] = gsub(" ", "", filenames[h])
				
			}		  
		}  

		old = filenames[h]
		#Pick out places where station numbers are the same
		ind = which(grepl(stationnum, filenames))
		string = paste("ep", stationnum, "r", sep = "") 
		x = "r"
		
		i = 1
		#loops through same station numbers and reassigns r where nessesary
		while(i <= length(ind)){
			if(i == 1) filenames[ind[i]] = gsub("r", "", string)
			else if(i==2) filenames[ind[i]] =  string
			else filenames[ind[i]] = gsub("r", paste("r", i-1, sep = ""), string)
			i = i+1
			
		}
		h = h+1
	}
	return (filenames)
}


    #  This function takes all .rdata files found in a directory and attempts to add them to the databases,
    #  It will only add if all the required information is also found and if specific stations have not 
	#  already been added. This code only loads the minilog and netmind raw data and bottom tables should
    #  be reproduced after any new data is added. 
	


	Rf2db2 = function(filepath, yearDesired){
		#Optional years list to only years specified
		# = list(...)
		
		#Open temporary database to store datafiles
		
		con <- dbConnect(dbDriver("SQLite"), "temp.db")
		
		print("             RFILES TO NETMIND & MINILOG DATABASE             ")
		
		#Get .rdata files out of a directory
		dirlist = c(list.files(filepath , pattern = NULL, all.files = FALSE, full.names = TRUE, recursive = TRUE))
		dirlist = grep(".rdata", dirlist, value = T)

		filesAdded = c()
		nfiles = length(dirlist)
		#if files exist loop through and add each to the temp database and names to a list
		if(nfiles>0)
		for(f in 1:nfiles) {
			fpath = dirlist[f]
			filename  = basename(fpath)
			fpath = get(load(fpath))
			filename = gsub("\\.", "_", filename)
			filesAdded = append(filesAdded, filename)
			dbWriteTable(con, filename, fpath, overwrite=T, row.names = T)
		}
		
		nfiles = length(filesAdded)
		# if .rdata files exist them check each for usefull data. 
		if(nfiles > 0)
		for (i in 1:nfiles) {
			if(substr(filesAdded[i], 1, 16) == "minilog_metadata"){
				yr = substr(filesAdded[i], 18, 21)
				#Don't continue if year undisered  
				if(!(yr %in% yearDesired || length(yearDesired) == 0)) next()
				added = 0
				matching = paste("minilog_basedata", yr, "rdata", sep = "_")
				#Only add data if the base data has also been found
				if(matching %in% filesAdded){
					con2 <- dbConnect(dbDriver("SQLite"), "minilog.db")
					#Load previously stored metadata to try to add to minilog database
					x = dbReadTable(con, filesAdded[i])
					x$filename = tolower(x$filename)
					x$timestamp = as.chron(as.numeric(x$timestamp))
					#Load data from minilog data base if data exists
					if(dbExistsTable(con2, substr(filesAdded[i], 1, 21))){ 
						y = dbReadTable(con2, substr(filesAdded[i], 1, 21))
						y$timestamp = string2chron(as.character(y$timestamp))
					}
					#if data does not yet exist then set up to write info found 
					else{ 
						y = x
						added = nrow(y)
					}
					#Load previously stored basedata to add to minilog database
					a = dbReadTable(con, matching)
					a$chron = string2chron(as.character(a$chron))
					#Load data from minilog data base if data exists
					if(dbExistsTable(con2,substr(matching, 1, 21))){ 
						b = dbReadTable(con2, substr(matching, 1, 21))
						names(b) = c( "mdate", "mtime",  "temperature", "depth", "chron", "unique_id" )
						b$chron = string2chron(as.character(b$chron))
					}
					#if data does not yet exist then set up to write info found 
					else{ 
						b = a
					}
					# For each station found add to database if it doesn't already exist
				for(j in 1:nrow(x)){
					if(!(x$filename[j] %in% y$filename)){
						unipos = nrow(y)+1 #Variable to keep track of unique ids
						newrow = x[j,]
						prevuni = newrow$unique_id
						#Update unique id
						newrow$unique_id = paste(yr, unipos, "minilog", sep = ".")
						y = rbind(y, newrow)
						added = added + 1
						#Add to metadata frame
						d = data.frame()
						#get relavent base data
						wanted =  subset(a , a$unique_id == prevuni)
						d = rbind(d, wanted)
						#update basedata unique id
						d$unique_id = paste(yr, unipos, "minilog", sep = ".")
						#Add to basedata frame
						b = rbind(b, d)
						
					}
				}				
				
				
				#if new stations found add the updated data frames to the databases
				if(added > 0){
					y$filename = tolower(y$filename)
					y$timestamp = as.character(y$timestamp)
					b$chron = as.character(b$chron )
					if(! is.null(y) || is.null(b)){
						dbWriteTable(con2, substr(filesAdded[i], 1, 21), y, overwrite = T, row.names = F )
						dbWriteTable(con2, substr(matching, 1, 21), b, overwrite = T, row.names = F )
						print(paste(added, " stations were added to Tables minilog_metadata_", yr, " and minilog_basedata_", yr, " in database minilog.db"  , sep = "")) 
					}
					else print("No previous Netmind data entered and no stations found in RFiles")
				}
				
				else print(paste("0 new stations were found to add to Tables minilog_metadata_", yr, " and minilog_basedata_", yr, " in database minilog.db", sep = ""))
				
				dbDisconnect(con2)	

			}
			else print(paste("Found file minilog_metadata_", yr, " however no matching basedata file found, metadata will not be added to the database", sep =""))  
			
		}
		#The following segment of code does similar actions as the above except for netmind data instead of minilog data
		if(substr(filesAdded[i], 1, 16) == "netmind_metadata"){
			yr = substr(filesAdded[i], 18, 21)
			if(!(yr %in% yearDesired || length(yearDesired) == 0)) next()          
			added = 0
			matching = paste("netmind_basedata", yr, "rdata", sep = "_")
			if(matching %in% filesAdded){
				con2 <- dbConnect(dbDriver("SQLite"), "netmind.db")
				x = dbReadTable(con, filesAdded[i])
				x$filename = tolower(x$filename)
				x$netmind_timestamp = as.chron(as.numeric(x$netmind_timestamp))
				if(dbExistsTable(con2,substr(filesAdded[i], 1, 21))){ 
					y = dbReadTable(con2, substr(filesAdded[i], 1, 21))
					names( y ) = c("filename", "unique_id", "yr", "netmind_timestamp", "trip", "setno", "station", "comments")
					y$netmind_timestamp = string2chron(as.character(y$netmind_timestamp))
				}
				else {
					y = x
					added = nrow(y)   
				}
				
				a = dbReadTable(con, matching)
				a$chrono = string2chron(as.character(a$chrono))
				
				if(dbExistsTable(con2,substr(matching, 1, 21))){ 
					b = dbReadTable(con2, substr(matching, 1, 21))
					names( b) = c( "ndate", "ntime",  "lat", "lon", "speed", "primary", "secondary", "doorspread", "depth", "chrono", "unique_id" )
					b$chrono = string2chron(as.character(b$chrono))
				}
				else b = a
				
				for(j in 1:nrow(x)){
					if(!(x$filename[j] %in% y$filename)){
						station = unlist(strsplit(gsub("ep", "", x$filename[j]), ""))
						z = 1
						stationnum = ""
						
						#Historical station num storage somtimes incorrect, updating with this method
						while(!is.na(as.numeric(station[z]))){
							stationnum = paste(stationnum, station[z], sep = "")
							z = z+1
						}
						x$station[j] = paste("ep", stationnum, substr(x$filename[j], z+2, nchar(x$filename[j])-4), sep = "") 
						unipos = nrow(y)+1
						newrow = x[j,]
						prevuni = newrow$unique_id
						newrow$unique_id = paste(yr, unipos, "netmind", sep = ".")
						y = rbind(y, newrow)
						added = added+1
						d = data.frame()
						wanted =  subset(a , a$unique_id == prevuni)
						d = rbind(d, wanted)
						d$unique_id = paste(yr, unipos, "netmind", sep = ".")
						b = rbind(b, d)
					}
				}				
				if(added > 0){
					y$filename = tolower(y$filename)
					y$netmind_timestamp = as.character(y$netmind_timestamp)
					b$chrono = as.character(b$chrono )
					if(! is.null(y) || is.null(b)){
						dbWriteTable(con2, substr(filesAdded[i], 1, 21), y, overwrite = T, row.names = F)
						dbWriteTable(con2, substr(matching, 1, 21), b, overwrite = T, row.names = F)
						print(paste(added, " stations were added to Tables netmind_metadata_", yr, " and netmind_basedata_", yr, " in database netmind.db"  , sep = "")) 
					}
					else print("No previous Netmind data entered and no stations found in RFiles")
				}
				
				
				else print(paste("0 new stations were found to add to Tables netmind_metadata_", yr, " and netmind_basedata_", yr, " in database netmind.db", sep = ""))
				
				dbDisconnect(con2)
			}
			
			else print(paste("Found file netmind_metadata_", yr, " however no matching basedata file found, metadata will not be added to the database", sep =""))  
			
		}
	}		 
	if(nfiles<1)print("There were no .rdata files in this directory")
	
	dbDisconnect(con)
	file.remove("temp.db")	
}


