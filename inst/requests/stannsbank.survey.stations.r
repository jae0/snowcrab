#Set up working environment
require(aegis.env)

#if (!exists("year.assessment")) {
  #year.assessment=lubridate::year(Sys.Date())      # year.assessment
  year.assessment=lubridate::year(Sys.Date()) -1   # or year previous to current
#}

p = bio.snowcrab::load.environment( year.assessment=year.assessment )

#Import data set from most recent assessment
dat = snowcrab.db('set.biologicals')
d = dat[,c('trip', 'set', 'lon','lat','set_type','distance', 'sa','timestamp', names(dat)[grep('ms.no',names(dat))])] #Required variables

#Recode species codes to common name
nn = as.numeric(substr(names(dat)[grep('ms.no',names(dat))],7,10))
nn = taxonomy.recode(from='spec',tolookup=nn)$vern
names(d)[9:ncol(d)] = nn

d[,c(9:ncol(d))]=round(d[,c(9:ncol(d))]*d$sa, digits=0)#Convert density estimates to raw numbers
d$yr = lubridate::year( d$timestamp )#add year

#Define SAB polygon
pp = read.csv(aegis::polygon_file('StAnnsMPA.csv'))
require(PBSmapping)
d$X = d$lon
d$Y = d$lat
d$EID = 1:nrow(d)
p1 = findPolys(d,pp)

dal = d[which(d$EID %in% p1$EID),] #only sets within SAB
dal=dal[dal$yr==max(dal$yr),] #Extract most recent year
dal[,9:ncol(dal)][dal[, 9:ncol(dal)]== 0]= NA #subsitute NA for 0 values
dal = Filter(function(x)!all(is.na(x)),dal) #retain only columns who have at least one value
dal[,9:ncol(dal)][is.na(dal[, 9:ncol(dal)])]= 0  #subsitute 0 for NA values
dal$tripset=paste(dal$trip, dal$set)

#Create output directory
outloc= file.path( p$project.outputdir, "requests", "SAB", year.assessment)
dir.create(outloc, recursive=T, warnings=F)

#create summary table for the year
sum.table= colSums(dal[,9:(ncol(dal)-6)])
fn=paste(year.assessment, "SAB.all.species.csv", sep=".")
write.csv(sum.table,file=paste(outloc, fn, sep="/"), row.names=TRUE)


#Determine number of Stomach Sampled within MPA
require(ROracle)
con=ROracle::dbConnect(DBI::dbDriver("Oracle"),dbname=oracle.snowcrab.server , username=oracle.snowcrab.user, password=oracle.snowcrab.password, believeNRows=F)

stomach=dbGetQuery(con, "SELECT SNSTOMACHDETAILS.TRIP, SNSTOMACHDETAILS.BOARD_DATE, SNSTOMACHDETAILS.SET_NO,
SNCRABSETS.SETCD_ID SET_TYPE, SNCRABSETS.START_LAT LATITUDE, SNCRABSETS.START_LONG LONGITUDE, SNSTOMACHDETAILS.EST_NUM_CAUGHT,
SNSTOMACHDETAILS.EST_DISCARD_WT, SNSTOMACHDETAILS.FISH_NO, SNSTOMACHDETAILS.SEXCD_ID,
SNSTOMACHDETAILS.FISH_LENGTH, SNSTOMACHDETAILS.SPECCD_ID, SNSTOMACHDETAILS.MEASURED_WGT,
SNSTOMACHDETAILS.QUANT_VALUE
FROM SNOWCRAB.SNCRABSETS SNCRABSETS, SNOWCRAB.SNSTOMACHDETAILS SNSTOMACHDETAILS
WHERE SNCRABSETS.TRIP = SNSTOMACHDETAILS.TRIP AND SNCRABSETS.SET_NO = SNSTOMACHDETAILS.SET_NO")
names(stomach)=tolower(names(stomach))
stomach$yr = lubridate::year( stomach$board_date )
stomach$tripset=paste(stomach$trip, stomach$set_no)
stom=stomach[stomach$yr==year.assessment,]

mpastom=stom[stom$tripset %in% dal$tripset,] #stomach samples within SAB
not.sab=stom[!stom$tripset %in% dal$tripset,]
i=which(not.sab$set_type=="11" & not.sab$latitude>45)
xtra.sampl=not.sab[i,] #stomach samples adjacent to SAB


print(paste("Number of total stations in SAB ", year.assessment,"= ", nrow(dal)))
print(paste("Number of extended sampling stations in SAB ", year.assessment,"= ", nrow(dal[which(dal$set_type=="11"),])))
print(paste("Number of extended sampling stations adjacent to SAB ", year.assessment,"= ", length(unique(xtra.sampl$tripset))))
print(paste("Number of species in all stations within SAB ", year.assessment,"= ", length(sum.table)))
print(paste("Number of stomach samples within SAB ", year.assessment,"= ", length(mpastom$trip)))
print(paste("Number of stomach samples adjacent to SAB ", year.assessment,"= ", length(xtra.sampl$trip)))
print(paste("Start date of SAB sampling in ", year.assessment,"= ", min(dal$timestamp)))
print(paste("End date of SAB sampling in ", year.assessment,"= ", max(dal$timestamp)))
print(paste("Number of actual sampling days in SAB in ", year.assessment,"= ", length(unique(lubridate::day(dal$timestamp)))))





