
loadfunctions( "snowcrab", functionname="initialise.local.environment.r") 

dat = snowcrab.db('set.complete')
pp = read.csv(find.ecomod.gis('StAnnsMPA.csv'))
require(PBSmapping)

dat$X = dat$lon
dat$Y = dat$lat
dat$EID = 1:nrow(dat)

p1 = findPolys(dat,pp)

da = dat[which(dat$EID %in% p1$EID),]

dal = da[,c('lon','lat','distance','chron', names(da)[grep('ms.size',names(da))])]
nn = as.numeric(substr(names(da)[grep('ms.size',names(da))],9,15))
nn = taxonomy.recode(from='spec',tolookup=nn)$vern
names(dal)[5:ncol(dal)] <- nn
 dal$yr = year(dal$chron)
  dal = dal[which(dal$yr>2003),]

dal = Filter(function(x)!all(is.na(x)),dal)
dal$yr = NULL

#species data
write.csv(dal,file='/home/ecomod_data/snowcrab/R/Requests/StAnnsSCSurveySpeciesData.csv',row.names=F)

dal$NSp = apply(dal[,5:ncol(dal)],1,function(x) length(x[!is.na(x)]))

x = dal[,c('lon','lat','distance','chron','Nsp')]
xx = x[which(years(x$chron)>2003),]
write.csv(xx,'/home/adam/ecomod/snowcrab/R/Requests/StAnnsSnowCrabSurvey.csv',row.names=F)