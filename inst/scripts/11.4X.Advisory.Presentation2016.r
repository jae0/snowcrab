#logbook 4X

if (0) { 
  # to delete
  require(bio.base)
  require(bio.snowcrab)
  require(bio.utilities)
  require(bio.spacetime)
  require(bio.polygons)
}

  
  if (!exists("current.year")) current.year=lubridate::year(Sys.Date())
  p = bio.snowcrab::load.environment( year.assessment=current.year)


outdir = file.path(project.datadirectory('bio.snowcrab'),'assessments',  p$year.assessment, "presentations", '4X')
dir.create(outdir,showWarnings=F)

#Map the Area
    require(PBSmapping)
    require(SpatialHub)
    bioLibrary ( 'spacetime','bio.utilities','bio.polygons' )
    logs = logbook.db('logbook')
    #logs$yr = logs$yr -1 # to make fishing year start of season ie march 2015 is fishing year 2014
    logs = makePBS(logs,polygon=F)
    logs = logs[which(logs$cfa=='cfa4x'),]
    lp = logs[,c('X','Y','EID','effort','landings','cpue','yr')]
    lp = na.omit(lp)
  	yy=sort(unique(lp$yr))
    effortgrids=list()
    catchgrids=list()
    cpuegrids=list()
    grid.polyData=list()
    grid.polyData[[1]]<-list()
    grid.polyData[[2]]<-list()
    grid.polyData[[3]]<-list()

	effort.levels=c(1,20,50,100,200,500,1000)
	catch.levels=c(1,200,500,1000,2000,5000,10000)
	cpue.levels=c(1,2,5,10,20,50,100)
	Min_lon=-66
	Max_lon=-63.33
	Min_lat=42
	Max_lat=44.5

   	for( y in 1:length(yy)) {
        effortgrids[[y]]<-gridData(subset(lp,yr==yy[y],c('EID','X','Y','effort')),lvls=effort.levels,bcol="YlOrRd",FUN=sum,border=NA,grid.size=2,sx=Min_lon,sy=Min_lat,ex=Max_lon,ey=Max_lat)
        grid.polyData[[1]][[y]] = effortgrids[[y]][[2]]

        catchgrids[[y]]<-gridData(subset(lp,yr==yy[y],c('EID','X','Y','landings')),lvls=catch.levels,bcol="YlOrRd",FUN=sum,border=NA,grid.size=2,sx=Min_lon,sy=Min_lat,ex=Max_lon,ey=Max_lat)
        grid.polyData[[2]][[y]] = catchgrids[[y]][[2]]

        cpuegrids[[y]]<-gridData(subset(lp,yr==yy[y],c('EID','X','Y','cpue')),lvls=cpue.levels,bcol="YlOrRd",FUN=mean,border=NA,grid.size=2,sx=Min_lon,sy=Min_lat,ex=Max_lon,ey=Max_lat)
        grid.polyData[[2]][[y]] = cpuegrids[[y]][[2]]

   		x11()
   		bioMap(xlim=c(-66.5,-63.1),ylim=c(42.7,44.8),boundaries='snowcrab',poly.lst=effortgrids[[y]][1:2],title=yy[y])
		ContLegend('bottomleft',lvls=effortgrids[[y]]$lvls,Cont.data=effortgrids[[y]],title='Effort (trap hauls)',inset=0.02,cex=0.8,bg='white')

   			savePlot(file.path(outdir,paste('logbook.effort',yy[y],'png',sep=".")),type='png')
		dev.off()

   		x11()
   		bioMap(xlim=c(-66.5,-63.1),ylim=c(42.7,44.8),boundaries='snowcrab',poly.lst=catchgrids[[y]][1:2],title=yy[y])
		ContLegend('bottomleft',lvls=catchgrids[[y]]$lvls,Cont.data=catchgrids[[y]],title='Catch (kg)',inset=0.02,cex=0.8,bg='white')

   			savePlot(file.path(outdir,paste('logbook.catch',yy[y],'png',sep=".")),type='png')
		dev.off()

   		x11()
   		bioMap(xlim=c(-66.5,-63.1),ylim=c(42.7,44.8),boundaries='snowcrab',poly.lst=cpuegrids[[y]][1:2],title=yy[y])
		ContLegend('bottomleft',lvls=cpuegrids[[y]]$lvls,Cont.data=cpuegrids[[y]],title='CPUE (kg/th)',inset=0.02,cex=0.8,bg='white')

   			savePlot(file.path(outdir,paste('logbook.cpue',yy[y],'png',sep=".")),type='png')
		dev.off()
	}





#Plot annual landings

	landings = with(logs,tapply(landings,yr,sum,na.rm=T))/1000
	x11()
	plot(landings, type="n", ylim=c(0,500), ylab="Landings (mt)", main="4X Landings by Year", xaxt="n", xlab="Year" )
	axis(1, at=1:length(landings), labels=names(landings))
	points(landings, col="red", pch=20)
	lines(landings, col="red")
	savePlot(file.path(outdir,'annual.landings.png'),type='png')


#Plot annual EFFORT
	traps = with(logs,tapply(effort,yr,sum,na.rm=T))
	x11()
	plot(traps, type="n", ylim=c(0,max(traps)), ylab="Trap Hauls", main="4X Trap Hauls by Year", xaxt="n", xlab="Year" )
	axis(1, at=1:length(traps), labels=names(traps))
	points(traps, col="red", pch=20)
	lines(traps, col="red")
	savePlot(file.path(outdir,'annual.effort.png'),type='png')



#Plot annual CPUE
	cpue.catch = with(subset(logs,!is.na(effort)&!is.na(landings)),tapply(landings,yr,sum,na.rm=T))
	cpue.traps = with(subset(logs,!is.na(effort)&!is.na(landings)),tapply(effort,yr,sum,na.rm=T))

	cpue = cpue.catch/cpue.traps 
	x11()
	plot(cpue, type="n", ylim=c(0,max(cpue)), ylab="Kg / Trap Haul", main="4X CPUE by Year", xaxt="n", xlab="Year" )
	axis(1, at=1:length(cpue), labels=names(cpue))
	points(cpue, col="red", pch=20)
	lines(cpue, col="red")
	savePlot(file.path(outdir,'annual.cpue.png'),type='png')





#--------------------------------------
# Clean for Catch Rate Data


   logs = logbook.db('rawdata.logbook',yrs=2003:2016)
		names(logs) = tolower(names(logs))
    	logs = logs[which(logs$cfa=='24W'),]
 		logs$lat =   round( as.numeric(substring(logs$latitude, 1,2)) + as.numeric(substring(logs$latitude, 3,6))/6000 ,6)
    	logs$lon = - round((as.numeric(substring(logs$longitude, 1,2)) + as.numeric(substring(logs$longitude, 3,6))/6000), 6)
   		logs$date_landed=logs$date_landed +8035200 #add 93 days to all dates (8,035,200 seconds)
		logs$yr=as.character( lubridate::year( logs$date_landed ))        #determine year
		logs$date_landed=logs$date_landed-8035200 #subtract 61 days
		logs$yr=as.character(as.numeric(logs$yr)-1) # subtract one year to give starting year of season rather than end year

#### *NB- Year refers to the starting year of season


		logs=logs[logs$yr!="2001",] #remove 2001 / 2002 season data as it is incomplete
		yrs=as.character(sort(as.numeric(unique(logs$yr))))
		logs$kg=logs$pro_rated_slip_wt_lbs/2.204626


		mts= c("October","November", "December", "January", "February", "March", "April","May")
		logs$month=as.character(months(logs$date_landed))  #populate month field
		logs=logs[logs$month %in% mts,]  #remove any data not within expected months

			landings=logs
			cleanlogsna=landings[is.finite(landings$num_of_traps),]
			cleanlogsna=landings[is.finite(landings$pro_rated_slip_wt_lbs),]
			cleanlogsna$lbspertrap=(cleanlogsna$pro_rated_slip_wt_lbs)/(cleanlogsna$num_of_traps)
			cleanlogsna=cleanlogsna[is.finite(cleanlogsna$lbspertrap),]
			cleanlogscr=cleanlogsna[cleanlogsna$lbspertrap < 800,]
			logs.fixed = cleanlogscr
			logs.fixed$catch = logs.fixed$pro_rated_slip_wt_lbs
			logs.fixed$effort = logs.fixed$num_of_traps


	logs.fixed$moy <- format(logs.fixed$date_fished,"%m")
		logs.fixed$area  = NA
	    logs.fixed[which(logs.fixed$lon> -64) , 'area'] <- 'East'
		logs.fixed[which(logs.fixed$lon<= -64 & logs.fixed$lon> -64.7) , 'area'] <- 'Central'
		logs.fixed[which(logs.fixed$lon<= -64.7) , 'area'] <- 'West'

##standardized catch rates adam catch rates by Month, by vessel, by year
standardized.catch.rates=T
if(standardized.catch.rates) {
						cp = logs.fixed
						cp$cp <- cp$pro_rated_slip_wt_lbs/cp$num_of_traps
						cp$logcp <- log(cp$cp)
						ssp = cp[,c('cp','logcp','yr','moy','area','vr_number')]
						ssp = ssp[complete.cases(ssp),]

						pred.grid <- expand.grid(yr=as.factor(unique(ssp$yr)),moy=as.factor(unique(ssp$moy)),
							vr_number=as.factor(unique(ssp$vr_number)),area=as.factor(unique(ssp$area)))

				linear.model=F
					if(linear.model) {
							glm.cpue <- glm(log(cp)~yr+as.factor(moy)+as.factor(vr_number)+as.factor(area),data=ssp)
							preds <-predict(glm.cpue,pred.grid,se=T,type='response')
						#estimate fitted values from the lognormal model using the delta method as per Faraway 2006, Extending the  linear model with R page 155
							preds$fit <- exp(preds$fit)
							preds$se.fit <- preds$fit * preds$se.fit
							devianceGLM(glm.cpue,ssp)
						}

				mixed.effects=T
					if(mixed.effects) {
							lme.cpue <- lme(log(cp)~yr+as.factor(moy)+as.factor(area),random=~1|vr_number, data=ssp)
							preds = exp(predict(lme.cpue,newdata=pred.grid,level=0,type='response'))
							#SE conditional on the random effects, ie only incorporates the variability in fixed effects
							DM <- model.matrix(eval(eval(lme.cpue$call$fixed)[-2]), pred.grid)
							predvar <- diag(DM %*% lme.cpue$varFix %*% t(DM))
							SE <- sqrt(preds+exp(predvar))
							preds = list(fit=preds,se.fit=SE)
						}

				#estimates for plot with 95%CI
					preds.m <-	aggregate(preds$fit,by=list(pred.grid$yr),FUN=mean)
					preds.s <- aggregate(preds$se.fit,by=list(pred.grid$yr),FUN=mean)

				#	plot(2002:2015,preds.m$x,type='b',lwd=2,xlab='Year',ylab=paste(expression(CPUE, 'lbs/trap')),ylim=c(0,140))
				#	arrows(x0=2002:2015,x1=2002:2015,y0=preds.m$x,y1=preds.m$x+preds.s$x,angle=90,length=0.03)
				#	arrows(x0=2002:2015,x1=2002:2015,y0=preds.m$x,y1=preds.m$x-preds.s$x,angle=90,length=0.03)
				#		savePlot(file.path(outdir,paste('standardized.catch.rate.png',sep=".")),type='png')
	}

x11()
	plot(cpue, type="n", ylim=c(0,max(cpue)), ylab="Kg / Trap Haul", main="4X CPUE by Year", xaxt="n", xlab="Year" )
	axis(1, at=1:length(cpue), labels=names(cpue))
	points(cpue, col="red", pch=20)
	lines(2:15,preds.m$x*0.453592, type='b',lwd=2)
	arrows(x0=2:15,x1=2:15,y0=preds.m$x*0.453592,y1=(preds.m$x+preds.s$x)*0.453592,angle=90,length=0.03)
	arrows(x0=2:15,x1=2:15,y0=preds.m$x*0.453592,y1=(preds.m$x-preds.s$x)*0.453592,angle=90,length=0.03)
	savePlot(file.path(outdir,paste('standardized.catch.rate.png',sep=".")),type='png')






# Import Observer Data for histograms
#------------------------------------------------------------

l = observer.db('rawdata',yrs=2014:2016)

#need to remove potential non-4x sets, interim, fix by using longitude.
l=l[l$LONGITUDE>63,]

# Remove CW's outside norms and remove production (pre-sorted) samples
a = l[l$FISH_LENGTH>50 & l$FISH_LENGTH<170,]

# convert lat's and long's to recognizable format for recode.areas
h=names(a)
h[h=="LATITUDE"] = "lat"
h[h=="LONGITUDE"] = "lon"
names(a) = h
a$lon=-a$lon
a$year=as.character(lubridate::year(a$BOARD_DATE))

#----------------------------------------
# create columns for area and CFA

a = a[filter.region.polygon(a,'cfa4x'),]
a$mn = months(a$BOARD_DATE)
a = a[which((a$mn %in% c('November','December') & a$year %in% 2015) | (a$mn %in% c('January','February','March') & a$year %in% 2016)),]
x=a
# --------------------------------------
# divide into 5 CC's and create histograms of CW and combine into one table

xCC1=x[x$SHELLCOND_CD==1,]
xCC2=x[x$SHELLCOND_CD==2,]
xCC3=x[x$SHELLCOND_CD==3,]
xCC4=x[x$SHELLCOND_CD==4,]
xCC5=x[x$SHELLCOND_CD==5,]
xhistCC1= hist(xCC1$FISH_LENGTH, breaks=seq(50, 170, by=3), freq=FALSE, plot=FALSE)
xhistCC2= hist(xCC2$FISH_LENGTH, breaks=seq(50, 170, by=3), freq=FALSE, plot=FALSE)
xhistCC3= hist(xCC3$FISH_LENGTH, breaks=seq(50, 170, by=3), freq=FALSE, plot=FALSE)
xhistCC4= hist(xCC4$FISH_LENGTH, breaks=seq(50, 170, by=3), freq=FALSE, plot=FALSE)
xhistCC5= hist(xCC5$FISH_LENGTH, breaks=seq(50, 170, by=3), freq=FALSE, plot=FALSE)
xplot= rbind(xhistCC1$counts, xhistCC2$counts, xhistCC3$counts, xhistCC4$counts, xhistCC5$counts)

xcc1perc= round((sum(xhistCC1$counts)/(sum(xhistCC1$counts)+ sum(xhistCC2$counts)+ sum(xhistCC3$counts)+ sum(xhistCC4$counts)+ sum(xhistCC5$counts)))*100, 1)
xcc2perc= round((sum(xhistCC2$counts)/(sum(xhistCC1$counts)+ sum(xhistCC2$counts)+ sum(xhistCC3$counts)+ sum(xhistCC4$counts)+ sum(xhistCC5$counts)))*100, 1)
xcc3perc= round((sum(xhistCC3$counts)/(sum(xhistCC1$counts)+ sum(xhistCC2$counts)+ sum(xhistCC3$counts)+ sum(xhistCC4$counts)+ sum(xhistCC5$counts)))*100, 1)
xcc4perc= round((sum(xhistCC4$counts)/(sum(xhistCC1$counts)+ sum(xhistCC2$counts)+ sum(xhistCC3$counts)+ sum(xhistCC4$counts)+ sum(xhistCC5$counts)))*100, 1)
xcc5perc= round((sum(xhistCC5$counts)/(sum(xhistCC1$counts)+ sum(xhistCC2$counts)+ sum(xhistCC3$counts)+ sum(xhistCC4$counts)+ sum(xhistCC5$counts)))*100, 1)

# --------------------------------------
# create stacked barplots with legends
# change years in main title

CFA4xplot=barplot(xplot [c(5:1),], space=0,names.arg=seq(50, 170, by=3)[-1],
 main="2015 / 2016", legend.text=c(paste("CC5 (",xcc5perc,"%)"),
 paste("CC4 (",xcc4perc,"%)"), paste("CC3 (",xcc3perc,"%)"), paste("CC2 (",xcc2perc,"%)"),
  paste("CC1 (",xcc1perc,"%)")), xlab="Carapace Width in mm", ylab="Number of Crab")
savePlot(file.path(outdir,paste('cc.histogram.2015','png',sep=".")),type='png')
#---------------------------------------------------------------------------------------------
# Mapping
#---------------------------------------------------------------------------------------------

# Format to allow conversion to EventData for pbsMapping
#--------------------------------------

# Observer Location Data

obs=l
    obs=obs[obs$LONGITUDE>63,]  #remove non-4X data
    obs=obs[is.finite(obs$LONGITUDE),]
    names(obs) = tolower(names(obs))
  		obs$landing_date=obs$landing_date +8035200 #add 93 days to all dates (8,035,200 seconds)
		obs$yr=as.character(lubridate::year(obs$landing_date) )        #determine year
		obs$landing_date=obs$landing_date-8035200 #subtract 61 days
		obs$yr=as.character(as.numeric(obs$yr)-1) # subtract one year to give starting year of season rather than end year

# Format to allow conversion to EventData for pbsMapping
#--------------------------------------

obs$X=-obs$longitude
obs$Y=obs$latitude
obs$EID=1:nrow(obs)
syr=min(as.numeric(as.character(obs$yr)))
eyr=max(as.numeric(as.character(obs$yr)))

x11()
makeMap(area='4X',addSummerStrata=F,main=2015)
oo = aggregate(cbind(X,Y)~trip_id,data=obs[which(obs$yr==2015),],FUN=mean)
oo$EID = 1:nrow(oo)
addPoints(lp[which(lp$yr==2015),],col='red',pch=16)
addPoints(oo,bg='green',pch=21)
cover()
savePlot(file.path(outdir,paste('observer.map.2015','png',sep=".")),type='png')

x11()
makeMap(area='4X',addSummerStrata=F,main=2014)
oo = aggregate(cbind(X,Y)~trip_id,data=obs[which(obs$yr==2014),],FUN=mean)
oo$EID = 1:nrow(oo)
addPoints(lp[which(lp$yr==2014),],col='red',pch=16)
addPoints(oo,bg='green',pch=21)
cover()
savePlot(file.path(outdir,paste('observer.map.2014','png',sep=".")),type='png')




#survey Index
set = snowcrab.db('set.biologicals')

d = set[filter.region.polygon(set,'cfa4x'),c('yr','totmass.male.com',"totmass.male.ncom","totmass.female.mat","totmass.female.imm")]
d = subset(d,yr>2001)

# zero inflated "hurdle" model
#d$nonzero = ifelse(d$totmass.male.com > 0, 1, 0)
#m1 <- glm(nonzero ~ as.factor(yr)-1, data = d, family = binomial(link = logit))
#m2 <- glm(totmass.male.com ~ as.factor(yr)-1, data = subset(d, nonzero == 1), family = Gamma(link = log))

 
set16 = read.csv(file.path(outdir,"4X_report_2016.csv"))
set16[set16=='NULL'] = 0
set16$yr=2016
set16$lon = -set16$longitude
set16$lat = set16$latitude
set16$totmass.male.com=as.numeric(set16$mature.male..kg.)/set16$surace_area.km.2/1000
set16$totmass.male.ncom=as.numeric(set16$immature.male..kg.)/set16$surace_area.km.2/1000
set16$totmass.female.mat=as.numeric(set16$mature.female..kg.)/set16$surace_area.km.2/1000
set16$totmass.female.imm=as.numeric(set16$immature.female..kg.)/set16$surace_area.km.2/1000

d = rbind(d,set16[filter.region.polygon(set16,'cfa4x'),c('yr','totmass.male.com',"totmass.male.ncom","totmass.female.mat","totmass.female.imm")])

Year = 2002:2016
offset=0.1
gmmm = sapply(Year,function(y){with(subset(d,yr==y),exp(mean(log(totmass.male.com+offset)))-offset)})
gmim = sapply(Year,function(y){with(subset(d,yr==y),exp(mean(log(totmass.male.ncom+offset)))-offset)})
gmmf = sapply(Year,function(y){with(subset(d,yr==y),exp(mean(log(totmass.female.mat+offset)))-offset)})
gmif = sapply(Year,function(y){with(subset(d,yr==y),exp(mean(log(totmass.female.imm+offset)))-offset)})

x11()
plot(Year,gmmm,type='b',col='blue',xlab='Year',ylab='Geometric mean t / km^2',main='Male biomass',ylim=c(0,0.26),xlim=c(2002,2016),pch=16)
lines(Year,gmim,type='b',col='blue',pch=17,lty=2)
legend('topleft',c('mature','immature'),col=c('blue'),lty=1:2,pch=16:17)
savePlot(file.path(outdir,paste('survey.trend.males','png',sep=".")),type='png')


x11()
plot(Year,gmmf,type='b',col='red',xlab='Year',ylab='Geometric mean t / km^2',main='Female biomass',ylim=c(0,0.26),xlim=c(2002,2016),pch=16)
lines(Year,gmif,type='b',col='red',pch=17,lty=2)
legend('topleft',c('mature','immature'),col=c('red'),lty=1:2,pch=16:17)
savePlot(file.path(outdir,paste('survey.trend.females','png',sep=".")),type='png')



bub.ex=1
bub.min=0.3
leg=c(1,2,5,10)

bioMap(xlim=c(-66.5,-63.1),ylim=c(42.7,44.8),boundaries='snowcrab',main='2016')
with(subset(set16,totmass.male.mat>0),points(X,Y,cex=bub.min+sqrt(totmass.male.mat)*bub.ex,pch=21,bg=rgb(0,1,0,0.2)))
with(subset(set16,totmass.male.mat==0),points(X,Y,pch=4,cex=bub.min))
legend('bottomleft',legend=c(0,leg),title= expression(t/km^2), pt.cex=c(bub.min,bub.min+sqrt(leg)*bub.ex),pch=c(4,rep(21,4)),pt.bg=rgb(0,1,0,0.2),bty='o',bg='white',box.col='white',inset=0.03,cex=1.2)
savePlot(file.path(outdir,paste('surveyMMbubbles.2016','png',sep=".")),type='png')

set = snowcrab.db('set.biologicals')
for (i in 2015:2002){
	x11()
	bioMap(xlim=c(-66.5,-63.1),ylim=c(42.7,44.8),boundaries='snowcrab',main=i)
	with(subset(set,lon<(-63)&yr==i&totmass.male.mat>0),points(lon,lat,cex=bub.min+sqrt(totmass.male.mat)*bub.ex,pch=21,bg=rgb(0,1,0,0.2)))
	with(subset(set,lon<(-63)&yr==i&totmass.male.mat==0),points(lon,lat,pch=4,cex=bub.min))
	legend('bottomleft',legend=c(0,leg),title= expression(t/km^2), pt.cex=c(bub.min,bub.min+sqrt(leg)*bub.ex),pch=c(4,rep(21,4)),pt.bg=rgb(0,1,0,0.2),bty='o',bg='white',box.col='white',inset=0.03,cex=1.2)
	savePlot(file.path(outdir,paste('surveyMMbubbles',i,'png',sep=".")),type='png')

}


set = snowcrab.db('set.biologicals')
for (i in 2015:2002){
	x11()
	bioMap('not4X',boundaries='snowcrab',main=i)
	with(subset(set,yr==i&totmass.male.mat>0),points(lon,lat,cex=bub.min+sqrt(totmass.male.mat)*bub.ex,pch=21,bg=rgb(0,1,0,0.2)))
	with(subset(set,yr==i&totmass.male.mat==0),points(lon,lat,pch=4,cex=bub.min))
	legend('bottomleft',legend=c(0,leg),title= expression(t/km^2), pt.cex=c(bub.min,bub.min+sqrt(leg)*bub.ex),pch=c(4,rep(21,4)),pt.bg=rgb(0,1,0,0.2),bty='o',bg='white',box.col='white',inset=0.03,cex=1.2)
	savePlot(file.path(outdir,paste('surveyMMbubblesNot4X',i,'png',sep=".")),type='png')

}


