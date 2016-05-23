#logbook 4X
loadfunctions( "snowcrab", functionname="initialise.local.environment.r") 
require(chron)
outdir = file.path(project.datadirectory('snowcrab'),'4X')
dir.create(outdir,showWarnings=F)
#Map the Area  	
if(map.logs) {
    require(PBSmapping)
    loadfunctions(c('spacetime','utility','polygons'))
    logs = logbook.db('logbook')
    logs$yr = logs$yr -1 # to make fishing year start of season ie march 2015 is fishing year 2014
    logs = makePBS(logs,polygon=F)
    logs = logs[which(logs$cfa=='cfa4x'),]
    lp = logs[,c('X','Y','EID','cpue','yr')]
     lp = na.omit(lp)
  	yy=unique(lp$yr)
   	for( y in yy) {
   				x11()
   				makeMap(area='4X',addSummerStrata=F)	
		   		b=lp[which(lp$yr==y),]
addPoints(b,pch=16, col='red')
		   		title(y)
		   		cover()
	   			savePlot(file.path(outdir,paste('logbook.locations',y,'png',sep=".")),type='png')
	dev.off()
		   			}
			 	}




if(map.logs.cpue) {
    require(PBSmapping)
    loadfunctions(c('spacetime','utility','polygons'))
    logs = logbook.db('logbook')

    res.lat = 0.045855 # = 2km on SS
    res.lon = 0.017938 # = 2km on SS
    xr = range(logs$lon,na.rm=T)
    yr = range(logs$lat,na.rm=T)

    grdx = seq(floor(xr[1]*10000)/10000,ceiling(xr[2]*10000)/10000,by=res.lon)
    grdy = seq(floor(yr[1]*10000)/10000,ceiling(yr[2]*10000)/10000,by=res.lat)
    grr = makeGrid(grdx,grdy,projection='LL')

    logs = makePBS(logs,polygon=F)
    logs = logs[which(logs$cfa=='cfa4x'),]
    lp = logs[,c('X','Y','EID','cpue','yr')]
     lp = na.omit(lp)
    lp$Z = log(lp$cpue)
    quants = c(0,seq(.05,.95,by=0.1))
   	qq = quantile(lp$Z,quants)
   	cols = color.code(n=length(qq))
  	yy=unique(lp$yr)
   	for( y in yy) {
   				x11()
   				makeMap(area='4X',addSummerStrata=F)	
		   		b=lp[which(lp$yr==y),]
			    fp = findCells(b,grr)
		        pdata = combineEvents(b,fp,FUN=mean)
		   		pp = makeProps(pdata,breaks=qq,propName='col',propVal=cols)
		   		addPolys(grr,polyProps = na.omit(pp), border=NULL)
		   		title(y)
	   			savePlot(file.path(outdir,paste('logbook.locations',y,'.png',sep=".")),type='png')
	
		   			}
		 	}



   logs = logbook.db('odbc.logbook',yrs=2003:2015)
		names(logs) = tolower(names(logs))
    	logs = logs[which(logs$cfa=='24W'),]
 		logs$lat =   round( as.numeric(substring(logs$latitude, 1,2)) + as.numeric(substring(logs$latitude, 3,6))/6000 ,6)
    	logs$lon = - round((as.numeric(substring(logs$longitude, 1,2)) + as.numeric(substring(logs$longitude, 3,6))/6000), 6)
   		logs$date_landed=logs$date_landed +8035200 #add 93 days to all dates (8,035,200 seconds)
		logs$yr=as.character(years(as.chron(logs$date_landed)))        #determine year
		logs$date_landed=logs$date_landed-8035200 #subtract 61 days
		logs$yr=as.character(as.numeric(logs$yr)-1) # subtract one year to give starting year of season rather than end year

#### *NB- Year refers to the starting year of season


		logs=logs[logs$yr!="2001",] #remove 2001 / 2002 season data as it is incomplete
		yrs=as.character(sort(as.numeric(unique(logs$yr))))
		logs$kg=logs$pro_rated_slip_wt_lbs/2.204626

 
		mts= c("October","November", "December", "January", "February", "March", "April","May")
		logs$month=as.character(months(logs$date_landed))  #populate month field
		logs=logs[logs$month %in% mts,]  #remove any data not within expected months


#determine landings by year

			yearly = aggregate(cbind(pro_rated_slip_wt_lbs,num_of_traps)~yr,data=logs,FUN=sum)
			names(yearly)= c("Year", "lbs",'effort')
			yearly$kg= yearly$lbs / 2.204626
			yearly$mt=(yearly$kg)/1000

#Plot annual landings
			x11()
			plot(yearly$mt, type="n", ylim=c(0,500), ylab="Landings (mt)", main="4X Landings by Year", xaxt="n", xlab="Year" )
			axis(1, at=1:length(yearly$kg), labels=yearly$Year)
			points(yearly$mt, col="red", pch=20)
			lines(yearly$mt, col="red")
	savePlot(file.path(outdir,paste('annual.landings.png',sep=".")),type='png')
	
#determine monthly landings by year
			monthly = aggregate(kg~month+yr,data=logs,FUN=sum)
			monthly$fm = recode(monthly$month,"'November'=1;'December'=2;'January'=3;'February'=4;'March'=5;'April'=6")
			monthly$kg[!is.finite(monthly$kg)]=0
			mts= c("November", "December", "January", "February", "March", "April")
			monthly=monthly[monthly$month %in% mts,]
			monthly$mt=monthly$kg/1000
			monthly = monthly[order(monthly$yr,monthly$fm),]

#plot monthly landings
			x11()
			plot(monthly$mt, type="n", ylim=c(0,max(monthly$mt)), ylab="Landings (mt)", xlim=c(1,length(mts)),
			 main="4X Landings by Month", xaxt="n", xlab="Month" )
			axis(1, at=1:length(mts), labels=mts)
			yu = as.numeric(max(as.numeric(monthly$yr)))
			yrs=(yu-2):yu
			cols = c("red", "green4", "black")
			for(i in 1:length(yrs)) {
			q = which(monthly$yr == yrs[i] & monthly$mt>0)
			lines(monthly$mt[q], col=cols[i],lty=1 )
			points(monthly$mt[q], col=cols[i], pch=20)
			}
		legend( "bottomright",legend=yrs,bty="n",lty=c(1,2),lwd=2, col=cols)
	savePlot(file.path(outdir,paste('monthly.landings.png',sep=".")),type='png')
	
#Plot annual EFFORT
	x11()
	traps = aggregate(num_of_traps~yr,data=logs,FUN=sum)
	plot(traps$yr,traps$num_of_traps, type="n", ylim=c(0,max(traps$num_of_traps)), ylab="Trap Hauls", main="4X Trap Hauls by Year", xlab="Year" )
		points(traps$yr,traps$num_of_traps, col="red", pch=20,type='b')
	savePlot(file.path(outdir,paste('trap.hauls.per.year','png',sep=".")),type='png')



#--------------------------------------
# Clean for Catch Rate Data
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
		
	
cpue = jackknifeCPUE(logs.fixed[,c('yr','catch','effort')],grouping=c('yr','area'))
cpue = cpue[which(cpue$area=='East'),]
# calculate mean catch rates for the season by CFA

		x11()
		ylims = c(0,max(cpue$cpue+cpue$cpue.var)*1.2)
		plot(cpue$yr, cpue$cpue,type="n", ylim=ylims, ylab="Lbs / Trap", main="Sambro", xlab="Year" )
		with(cpue,points(yr,cpue, col="red", pch=20,type='b'))
		with(cpue,arrows(x0=as.numeric(yr),y0=cpue-cpue.var,y1=(cpue+cpue.var), col="red", length=0))
		
savePlot(file.path(outdir,paste('yearly.cpue.sambro.png',sep=".")),type='png') 

##--
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
					
					plot(2002:2014,preds.m$x,type='b',lwd=2,xlab='Year',ylab=paste(expression(CPUE, 'lbs/trap')),ylim=c(0,140))
					arrows(x0=2002:2014,x1=2002:2014,y0=preds.m$x,y1=preds.m$x+preds.s$x,angle=90,length=0.03)
					arrows(x0=2002:2014,x1=2002:2014,y0=preds.m$x,y1=preds.m$x-preds.s$x,angle=90,length=0.03)
						savePlot(file.path(outdir,paste('standardized.catch.rate.png',sep=".")),type='png')
	}			


#Delury 
	if(Delury) {
				logs.fixed$weekno = as.numeric(round((logs.fixed$date_fished - as.POSIXct(paste(logs.fixed$yr,11,1,sep="-")))/604800)+1)
				logs.fixed$landings  = logs.fixed$catch
				logs.fixed$effort  = logs.fixed$num_of_traps
				yr = unique(logs.fixed[,c('yr','area')])
				outd = list()
				outl = list()
				for(i in 1:nrow(yr)){
				outd[[i]] <- delury.leslie(x= logs.fixed[which(logs.fixed$yr==yr[i,1] & logs.fixed$area==yr[i,2]),],estimate=c('delury'),day.or.week='week')
				if(outd[[i]]$coef[1] != 1234) {
				ti = paste(yr[i,1],yr[i,2],sep="-")
				title(ti)
				savePlot(file.path(outdir,paste('delury',ti,'png',sep=".")),type='png')
				dev.off()
				outl[[i]] <- delury.leslie(x= logs.fixed[which(logs.fixed$yr==yr[i,1] & logs.fixed$area==yr[i,2]),],estimate=c('leslie'),day.or.week='week')
				if(outl[[i]]$coef[1] != 1234) {
				ti = paste(yr[i,1],yr[i,2],sep="-")
				title(ti)
				savePlot(file.path(outdir,paste('leslie',ti,'png',sep=".")),type='png')
				dev.off()
				}
			}
			}
}
		
il = length(outl)

	
 ##########- Following script calculates monthly catch rates by area (in lbs)

#sep log data into east and west
		logsW = logs.fixed
		monthcr = jackknifeCPUE( logsW[,c('catch','effort','month','yr')], grouping=c("month", "yr")  )
		mts= c("November", "December", "January", "February", "March", "April")
		monthcr=monthcr[monthcr$month %in% mts,]
		monthcr = na.omit(monthcr)
		monthcr$fm = recode(monthcr$month,"'November'=1;'December'=2;'January'=3;'February'=4;'March'=5;'April'=6")
		monthcr = monthcr[order(monthcr$yr,monthcr$fm),]
		x11()
		ylims = c(0,200)
		cols = c("red", "green4", "black",'blue')
		plot(1:6, 1:6,type="n", ylim=ylims, ylab="Lbs / Trap", main="Monthly 4X Catch Rates by Year", xlab="Month" ,
			xaxt='n')
		axis(1, at=1:length(mts), labels=mts)
		ny = 2011:2014
		for(y in 1:length(ny)) {
			with(monthcr[monthcr$yr==ny[y],],points(fm,cpue,type='b',col=cols[y]))
			with(monthcr[monthcr$yr==ny[y],],arrows(x0=fm,y0=cpue-sqrt(cpue.var),y1=cpue+sqrt(cpue.var),length=0,col=cols[y]))
		}
		legend('topleft',legend=ny,col=cols,lty=rep(1,4),pch=rep(1,4),bty='n')
	savePlot(file.path(outdir,'monthly.landings.by.year.png'),type='png')
			

if(map.logs.month) {
    require(PBSmapping)
    loadfunctions(c('spacetime','utility','polygons'))
   logs = logs.fixed
    logs = makePBS(logs,polygon=F)
    logs$fm = recode(logs$month,"'November'=1;'December'=2;'January'=3;'February'=4;'March'=5;'April'=6")
	
     lp = logs[,c('X','Y','EID','yr','fm','month')]
     lp = na.omit(lp)
     lp = lp[which(lp$yr==2014),]
  	yy=unique(lp$month)
   	for( y in yy) {
   				x11()
   				makeMap(area='4X',addSummerStrata=F)	
		   		b=lp[which(lp$month==y),]
addPoints(b,pch=16, col='red')
		   		title(y)
		   		cover()
	   			savePlot(file.path(outdir,paste('monthly.logbook.locations',y,'png',sep=".")),type='png')
	dev.off()
		   			}
			 	}

#effort by month
x11()
	  logs.fixed$fm = recode(logs.fixed$month,"'November'=1;'December'=2;'January'=3;'February'=4;'March'=5;'April'=6")

	traps = aggregate(num_of_traps~yr+fm+month,data=logs.fixed,FUN=sum)
		cols = c("red", "green4", "black",'blue')
		ylims = c(0,9000)
traps = traps[order(traps$fm),]
		plot(1:6, 1:6,type="n", ylim=ylims, ylab="Trap Hauls", main="Monthly 4X Trap Hauls", xlab="Month" ,
			xaxt='n')
		axis(1, at=1:length(mts), labels=mts)
	ny = 2011:2014
		for(y in 1:length(ny)) {
			with(traps[traps$yr==ny[y],],points(fm,num_of_traps,type='b',col=cols[y]))
		}
				legend('topleft',legend=ny,col=cols,lty=rep(1,4),pch=rep(1,4),bty='n')

	savePlot(file.path(outdir,paste('trap.hauls.per.month','png',sep=".")),type='png')


### Determine number of active vessel by month
x11()
	boats = aggregate(vr_number~month+yr,data=unique(logs[,c('vr_number','month','yr')]),FUN=length)
	mts= c("November", "December", "January", "February", "March", "April")
	boats=boats[boats$month %in% mts,]
	boats$fm = recode(boats$month,"'November'=1;'December'=2;'January'=3;'February'=4;'March'=5;'April'=6")
	boats = boats[order(boats$yr,boats$fm),]
	ylims = c(0,12)
	cols = c("red", "green4", "black",'blue')

	plot(1:6, 1:6,type="n", ylim=ylims, ylab="# of Vessels", main="4X Vessels Active By Month", xlab="Month" ,
	xaxt='n')
	axis(1, at=1:length(mts), labels=mts)
	ny = 2011:2014
	for(y in 1:length(ny)) {
		with(boats[boats$yr==ny[y],],points(fm,pch=16, vr_number,type='b',col=cols[y]))
	}
legend('topleft',legend=ny,col=cols,lty=rep(1,4),pch=rep(1,4),bty='n')
savePlot(file.path(outdir,paste('active.vessels.by.month','png',sep=".")),type='png')





# Import Observer Data for histograms
#------------------------------------------------------------

l = observer.db('odbc',yrs=2014:2015)

#need to remove potential non-4x sets, interim, fix by using longitude.
	l=l[l$LONGITUDE>63,]

	require (chron)
	require(lattice)




# --------------------------------------
# Remove CW's outside norms and remove production (pre-sorted) samples

		a = l[l$FISH_LENGTH>50 & l$FISH_LENGTH<170,]
		# --------------------------------------
# convert lat's and long's to recognizable format for recode.areas
	h=names(a)
	h[h=="LATITUDE"] = "lat"
	h[h=="LONGITUDE"] = "lon"
	names(a) = h
	a$lon=-a$lon
	a$year=as.character(years(as.chron(a$BOARD_DATE)))

#----------------------------------------
# create columns for area and CFA

a = a[filter.region.polygon(a,'cfa4x'),]
a$mn = months(a$BOARD_DATE)
a = a[which((a$mn %in% c('November','December') & a$year %in% 2013) | (a$mn %in% c('January','February','March') & a$year %in% 2014)),]
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
 main="2013 / 2014", legend.text=c(paste("CC5 (",xcc5perc,"%)"),
 paste("CC4 (",xcc4perc,"%)"), paste("CC3 (",xcc3perc,"%)"), paste("CC2 (",xcc2perc,"%)"),
  paste("CC1 (",xcc1perc,"%)")), xlab="Carapace Width in mm", ylab="Number of Crab")
savePlot(file.path(outdir,paste('cc.histogram.2013','png',sep=".")),type='png')
#---------------------------------------------------------------------------------------------
# Mapping
#---------------------------------------------------------------------------------------------

# Format to allow conversion to EventData for pbsMapping
#--------------------------------------
makeMap(area='4X',addSummerStrata=F)	
lp = makePBS(logs.fixed,polygon=F)
lp = na.omit(lp[,c('X','Y','EID','yr')])
addPoints(lp[which(lp$yr==2014),],col='red',pch=16)
cover()
savePlot(file.path(outdir,paste('logbook.map.2014','png',sep=".")),type='png')

# Observer Location Data

obs=l
    obs=obs[obs$LONGITUDE>63,]  #remove non-4X data
    obs=obs[is.finite(obs$LONGITUDE),]
    names(obs) = tolower(names(obs))
  		obs$landing_date=obs$landing_date +8035200 #add 93 days to all dates (8,035,200 seconds)
		obs$yr=as.character(years(as.chron(obs$landing_date)))        #determine year
		obs$landing_date=obs$landing_date-8035200 #subtract 61 days
		obs$yr=as.character(as.numeric(obs$yr)-1) # subtract one year to give starting year of season rather than end year
  
# Format to allow conversion to EventData for pbsMapping
#--------------------------------------

obs$X=-obs$longitude
obs$Y=obs$latitude
obs$EID=1:nrow(obs)
syr=min(as.numeric(as.character(obs$yr)))
eyr=max(as.numeric(as.character(obs$yr)))

makeMap(area='4X',addSummerStrata=F)
oo = aggregate(cbind(X,Y)~trip_id,data=obs[which(obs$yr==2014),],FUN=mean)
oo$EID = 1:nrow(oo)
addPoints(oo,col='red',pch=16)
cover()
savePlot(file.path(outdir,paste('observer.map.2014','png',sep=".")),type='png')




#survey Index

set = snowcrab.db('set.complete')   

st = unique(set[,c('yr','lon','lat','station')])
st2013 = st[which(st$yr==2013),]
st2014 = st[which(st$yr==2014),]

s3 = makePBS(st2013,polygon=F)
s4 = makePBS(st2014,polygon=F)

makeMap(area='4X',addSummerStrata=F)
addPoints(s3,pch=16,col='green')
addPoints(s4,pch=16,col='red')
savePlot(file.path(outdir,paste('survey.stations.2014','png',sep=".")),type='png')


   td =  get.time.series ( from.file=T,reduced.stations=T )
  td = td[which(td$region=='cfa4x'),]
td = td[order(td$year),]
  plot(td$year,td$mean,type='b',col='red',xlab='Year',ylab='Geometric mean t / km^2',ylim=c(0,0.8),pch=16)
arrows(x0=td$year,y0=td$mean-td$lb,y1=td$mean+td$ub,col='red',length=0)
savePlot(file.path(outdir,paste('survey.R0.trend.reduced.stations','png',sep=".")),type='png')



ser = set[which(set$station %in% st),]
aggregate(R0.mass~yr,data=ser,FUN=geomean)