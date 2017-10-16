#logbook 4X

warning( "The 11.4X.Advisory.Presentation* files need to be merged and cleaned up .. Ben")

if (!exists("year.assessment")) year.assessment=lubridate::year(Sys.Date())
p = bio.snowcrab::load.environment( year.assessment=year.assessment )

outdir = file.path(project.datadirectory('bio.snowcrab'),'assessments',  p$year.assessment, "presentations", '4X')
dir.create(outdir,showWarnings=T)

#Map the Area
    require(PBSmapping)
    require(SpatialHub)
    bioLibrary ( 'spacetime','bio.utilities','bio.polygons' )
    
#do new db pull, as needed, landings data from past winter, often not there from assessment time
    logbook.db(DS='rawdata.logbook.redo',yrs=1996:p$year.assessment)
    logbook.db(  DS="rawdata.licence.redo" ) 
    logbook.db(  DS="rawdata.areas.redo" )
    logbook.db('logbook.redo', p=p)
    
#Import logbook record for use
#*NB: yr variable refers to starting year of season so winter 2017 is yr=2016
    logs = logbook.db('logbook')
    logs = logs[which(logs$cfa=='cfa4x'),] #only 4X logs
    
    ys=c((max(logs$yr)-11):max(logs$yr)) #last 12 years, can as desired  
    logs=logs[logs$yr %in% ys,]
    
   # mts= c("October","November", "December", "January", "February", "March", "April","May")
    logs$month=as.character(lubridate::month(logs$date.landed))  #populate month field
    #logs=logs[logs$month %in% mts,]  #remove any data not within expected months
    
    
#Prep for mapping
    logs = makePBS(logs,polygon=F)
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

   		plot.new()
   		bioMap(xlim=c(-66.5,-63.1),ylim=c(42.7,44.8),boundaries='snowcrab',poly.lst=effortgrids[[y]][1:2],title=yy[y])
		ContLegend('bottomleft',lvls=effortgrids[[y]]$lvls,Cont.data=effortgrids[[y]],title='Effort (trap hauls)',inset=0.02,cex=0.8,bg='white')

   			savePlot(file.path(outdir,paste('logbook.effort',yy[y],'png',sep=".")),type='png')
		dev.off()

   		plot.new()
   		bioMap(xlim=c(-66.5,-63.1),ylim=c(42.7,44.8),boundaries='snowcrab',poly.lst=catchgrids[[y]][1:2],title=yy[y])
		ContLegend('bottomleft',lvls=catchgrids[[y]]$lvls,Cont.data=catchgrids[[y]],title='Catch (kg)',inset=0.02,cex=0.8,bg='white')

   			savePlot(file.path(outdir,paste('logbook.catch',yy[y],'png',sep=".")),type='png')
		dev.off()

   		plot.new()
   		bioMap(xlim=c(-66.5,-63.1),ylim=c(42.7,44.8),boundaries='snowcrab',poly.lst=cpuegrids[[y]][1:2],title=yy[y])
		ContLegend('bottomleft',lvls=cpuegrids[[y]]$lvls,Cont.data=cpuegrids[[y]],title='CPUE (kg/th)',inset=0.02,cex=0.8,bg='white')

   			savePlot(file.path(outdir,paste('logbook.cpue',yy[y],'png',sep=".")),type='png')
		dev.off()
	}

#Plot annual landings

	landings = aggregate(landings~yr,data=logs,FUN=sum)
	landings$landings=landings$landings/1000
	plot.new()
	plot(landings$landings, type="n", ylim=c(0,500), ylab="Landings (mt)", main="4X Landings by Year", xaxt="n", xlab="Year" )
	axis(1, at=1:length(landings$landings), labels=landings$yr)
	points(landings$landings, col="red", pch=20)
	lines(landings$landings, col="red")
	savePlot(file.path(outdir,'annual.landings.png'),type='png')

	#determine monthly landings by year
	monthly = aggregate(landings~month+yr,data=logs,FUN=sum)
	monthly$fm = recode(monthly$month,"'11'='November';'12'='December';'1'='January';'2'='February';'3'='March'; '4'='April'")
	mts= c("November", "December", "January", "February", "March")
	monthly$fm=factor(monthly$fm, levels=mts)
	monthly$landings[!is.finite(monthly$landings)]=0
	
	monthly$mt=monthly$landings/1000
	monthly = monthly[order(monthly$yr,monthly$fm),]
	
	#plot monthly landings
	plot.new()
	yu = as.numeric(max(as.numeric(monthly$yr)))
	yrs=(yu-2):yu
	tbz=monthly[monthly$yr %in% c(yrs),]
	plot(monthly$mt, type="n", ylim=c(0,max(tbz$mt)), ylab="Landings (mt)", xlim=c(1,length(mts)),
	     main="4X Landings by Month", xaxt="n", xlab="Month" )
	axis(1, at=1:length(mts), labels=mts)
	cols = c("red", "green4", "black")
	for(i in 1:length(yrs)) {
	  q = which(monthly$yr == yrs[i] & monthly$mt>0)
	  lines(monthly$mt[q], col=cols[i],lty=1 )
	  points(monthly$mt[q], col=cols[i], pch=20)
	}
	legend( "topright",legend=yrs,bty="n",lty=c(1,2),lwd=2, col=cols, cex=1.3)
	
	savePlot(file.path(outdir,'monthly.landings.png'),type='png')
	
	### Determine number of active vessel by month
	plot.new()
	boats = aggregate(cfv~month+yr,data=unique(logs[,c('cfv','month','yr')]),FUN=length)
	mts= c("November", "December", "January", "February", "March")
	boats$fm = recode(boats$month,"'11'='November';'12'='December';'1'='January';'2'='February';'3'='March'")
	boats$fm=factor(boats$fm, levels=mts)
	boats$cfv[!is.finite(boats$cfv)]=0
	boats = boats[order(boats$yr,boats$fm),]
	ylims = c(0,7)
	cols = c("red", "green4", "black",'blue')
	
	plot(1:5, 1:5,type="n", ylim=ylims, ylab="# of Vessels", main="4X Vessels Active By Month", xlab="Month" ,
	     xaxt='n')
	axis(1, at=1:length(mts), labels=mts)
	yu = as.numeric(max(as.numeric(boats$yr)))
	ny=(yu-2):yu
	for(y in 1:length(ny)) {
	  with(boats[boats$yr==ny[y],],points(fm,pch=16, cfv,type='b',col=cols[y]))
	}
	legend('topright',legend=ny,col=cols,lty=rep(1,4),pch=rep(1,4),bty='n')
	savePlot(file.path(outdir,paste('active.vessels.by.month','png',sep=".")),type='png')
	
	

#Plot annual EFFORT
	traps = with(logs,tapply(effort,yr,sum,na.rm=T))
	plot.new()
	plot(traps, type="n", ylim=c(0,max(traps)), ylab="Trap Hauls", main="4X Trap Hauls by Year", xaxt="n", xlab="Year" )
	axis(1, at=1:length(traps), labels=names(traps))
	points(traps, col="red", pch=20)
	lines(traps, col="red")
	savePlot(file.path(outdir,'annual.effort.png'),type='png')



#Plot annual CPUE
#maybe jackknife this?
	cpue.catch = with(subset(logs,!is.na(effort)&!is.na(landings)),tapply(landings,yr,sum,na.rm=T))
	cpue.traps = with(subset(logs,!is.na(effort)&!is.na(landings)),tapply(effort,yr,sum,na.rm=T))

	cpue = cpue.catch/cpue.traps 
	plot.new()
	plot(cpue, type="n", ylim=c(0,max(cpue)), ylab="Kg / Trap Haul", main="4X CPUE by Year", xaxt="n", xlab="Year" )
	axis(1, at=1:length(cpue), labels=names(cpue))
	points(cpue, col="red", pch=20)
	lines(cpue, col="red")
	abline(h=mean(cpue), col="blue", lty=3, lwd=1)
	savePlot(file.path(outdir,'annual.cpue.png'),type='png')

#add jackknife estimates of error
	jack=logs[logs$cfa0=="cfa4x",] #create new dataframe (direct copy of logs)
	names(jack)[names(jack) == 'landings'] <- 'catch'
	names(jack)[names(jack) == 'cfa'] <- 'area'
	jack$yr=as.character(jack$yr)
	jack=jack[,c('yr','catch','effort','area')]
	jack=na.omit(jack)
	
	cpue = jackknifeCPUE(jack,grouping=c('yr'))
	
	plot.new()
	ylims = c(0,max(cpue$cpue+cpue$cpue.var)*1.2)
	plot(cpue$yr, cpue$cpue,type="n", ylim=ylims, ylab="Kgs / Trap", main="4X Catch Rates by Year", xlab="Year" )
	with(cpue,points(yr,cpue, col="red", pch=20,type='b'))
	with(cpue,arrows(x0=as.numeric(yr),y0=cpue-cpue.var,y1=(cpue+cpue.var), col="red", length=0))
	
	savePlot(file.path(outdir,paste('annual.cpue.jack.png',sep=".")),type='png')
	
	rm(jack)
	
# Import Observer Data for histograms
#------------------------------------------------------------
#Need to do new data pull if first time being run since assessment (~15 minutes)
observer.db( DS="rawdata.redo", yrs=2003:p$year.assessment )

yu = p$year.assessment
yrs=(yu-4):yu

l = observer.db('rawdata',yrs=yrs)
l$tripset=paste(l$TRIP, l$SET_NO, sep=":")
#need to remove potential non-4x sets, interim, fix by using longitude.
#l=l[l$LONGITUDE>63,]
# Remove CW's outside norms and remove production (pre-sorted) samples

a = l[l$FISH_LENGTH>50 & l$FISH_LENGTH<170,]
# --------------------------------------
# convert lat's and long's to recognizable format for recode.areas
	h=names(a)
	h[h=="LATITUDE"] = "lat"
	h[h=="LONGITUDE"] = "lon"
	names(a) = h
	a$lon=-a$lon
	a$LANDING_DATE=a$LANDING_DATE +8035200 #add 93 days to all dates (8,035,200 seconds)
  a$yr=as.character( lubridate::year(a$LANDING_DATE ))        #determine year
  a$LANDING_DATE=a$LANDING_DATE-8035200 #subtract 93 days
  a$yr=as.character(as.numeric(a$yr)-1) # subtract one year to give starting year of season rather than end year

#----------------------------------------
# create columns for area and CFA
a = a[filter.region.polygon(a,'cfa4x'),]
  yu = p$year.assessment
  yrs=(yu-4):(yu-1)

for (y in yrs)  {
plot.new()
x=a[a$yr==y,]
# --------------------------------------
# divide into 5 CC's and create histograms of CW and combine into one table

xCC1=x[x$SHELLCOND_CD==1,]
xCC2=x[x$SHELLCOND_CD==2,]
xCC3=x[x$SHELLCOND_CD==3,]
xCC4=x[x$SHELLCOND_CD==4,]
xCC5=x[x$SHELLCOND_CD==5,]
xhistCC1= hist(xCC1$FISH_LENGTH, breaks=seq(50, 170, by=3),  plot=FALSE)
xhistCC2= hist(xCC2$FISH_LENGTH, breaks=seq(50, 170, by=3),  plot=FALSE)
xhistCC3= hist(xCC3$FISH_LENGTH, breaks=seq(50, 170, by=3),  plot=FALSE)
xhistCC4= hist(xCC4$FISH_LENGTH, breaks=seq(50, 170, by=3),  plot=FALSE)
xhistCC5= hist(xCC5$FISH_LENGTH, breaks=seq(50, 170, by=3),  plot=FALSE)
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
 main=paste(y, as.character (as.numeric(y)+1), sep="/"), legend.text=c(paste("CC5 (",xcc5perc,"%)"),
 paste("CC4 (",xcc4perc,"%)"), paste("CC3 (",xcc3perc,"%)"), paste("CC2 (",xcc2perc,"%)"),
  paste("CC1 (",xcc1perc,"%)")), xlab="Carapace Width in mm", ylab="Number of Crab")
savePlot(file.path(outdir,paste('cc.histogram',y,'png',sep=".")),type='png')
}
#Compile some observer stats for presentation
  

  out=NULL
  
  for (y in yrs){
      j=a[a$yr==y,]
      j=j[!duplicated(j$tripset),]
      # determine total mt observed by area
      observedmt= (sum(j$EST_KEPT_WT, na.rm=T)/1000)
      #observedmt
      
      # determine # of traps sampled
      sampledtraps= length(j$SET_NO)
      #sampledtraps
      
      
      # --------------------------------------
      # determine # of traps observed (bycatch, landings, etc)
      
      observedtraps= sum(j$NUM_HOOK_HAUL, na.rm=T)
      observedtraps
      
      # --------------------------------------
      # determine # of traps observed (bycatch, landings, etc)
      
      observedtraps= sum(j$NUM_HOOK_HAUL, na.rm=T)
      observedtraps
      
      #print(paste(y,c,"Observed mt=", observedmt, "Traps Sampled=", sampledtraps, "Observed Traps=", observedtraps, sep=" "))
      out=rbind(out,c(y, observedmt, sampledtraps, observedtraps))
    }
  
  out=as.data.frame(out)
  b=out
  out=b
  
  names(out)=c("Year", "Observed", "Traps Sampled", "Traps Observed")
  
  ys=yrs
  
  out$Landings=NA
  cfa=unique(out$Area)
  for (y in ys){
      out$Landings[out$Year==y]=round(landings$landings[landings$yr==y],0)
  }

  
  
  
  out$Observed=round(as.numeric(as.character(out$Observed)))
  out$Percent=round(out$Observed/out$Landings*100,1)
  names(out)=c("Year", "Observed (mt)", "Traps Sampled", "Traps Observed", "Landings (mt)", "% Observed")
  
 output=out[,c(-2, -5)]
 
  plot.new()
  gridExtra::grid.table(out, theme=ttheme_default(), rows=NULL)
  savePlot(file.path(outdir,paste('observersummary.png',sep=".")),type='png')
#BZ TODO- Add number of observed trips?
  
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
		obs$landing_date=obs$landing_date-8035200 #subtract 93 days
		obs$yr=as.character(as.numeric(obs$yr)-1) # subtract one year to give starting year of season rather than end year

# Format to allow conversion to EventData for pbsMapping
#--------------------------------------
obs$tripset=paste(obs$trip, obs$set_no, sep=":")
obs$X=-obs$longitude
obs$Y=obs$latitude
obs$EID=1:nrow(obs)
syr=min(as.numeric(as.character(obs$yr)))
eyr=max(as.numeric(as.character(obs$yr)))

#the follow (by year), plots observed sets over commercial log sets
#comparison of overlap of fishing and observed locations

plot.new()

for (y in yrs){
    makeMap(area='4X',addSummerStrata=F,main=y)
    oo = aggregate(cbind(X,Y)~tripset,data=obs[which(obs$yr==y),],FUN=mean)
    oo$EID = 1:nrow(oo)
    addPoints(lp[which(lp$yr==y),],col='red',pch=16)
    addPoints(oo,bg='green',pch=21, cex=0.6)
    
    cover()
    savePlot(file.path(outdir,paste('observer.map',y,'png',sep=".")),type='png')
}

#---------------------------------------
#Plotting of raw survey catches
#---------------------------------------

#Set some plotting variables
bub.ex=1
bub.min=0.3
leg=c(1,2,5,10)


set = snowcrab.db('set.biologicals')
for (i in yrs){
      plot.new()
      bioMap(xlim=c(-66.5,-63.1),ylim=c(42.7,44.8),boundaries='snowcrab',main=i)
      with(subset(set,lon<(-63)&yr==i&totmass.male.mat>0),points(lon,lat,cex=bub.min+sqrt(totmass.male.mat)*bub.ex,pch=21,bg=rgb(0,1,0,0.2)))
      with(subset(set,lon<(-63)&yr==i&totmass.male.mat==0),points(lon,lat,pch=4,cex=bub.min))
      legend('bottomleft',legend=c(0,leg),title= expression(t/km^2), pt.cex=c(bub.min,bub.min+sqrt(leg)*bub.ex),pch=c(4,rep(21,4)),pt.bg=rgb(0,1,0,0.2),bty='o',bg='white',box.col='white',inset=0.03,cex=1.2)
  savePlot(file.path(outdir,paste('survey.mat.male.bub',i,'png',sep=".")),type='png')
}



#---------------------------------------
# Following Script allows for "on-the fly" comparison
# from current survey data (raw data from survey vessel)
#---------------------------------------
#The following 
# #survey Index
# set = snowcrab.db('set.biologicals')
# 
# d = set[filter.region.polygon(set,'cfa4x'),c('yr','totmass.male.com')]
# d = subset(d,yr>2001)
# 
# # zero inflated "hurdle" model
# #d$nonzero = ifelse(d$totmass.male.com > 0, 1, 0)
# #m1 <- glm(nonzero ~ as.factor(yr)-1, data = d, family = binomial(link = logit))
# #m2 <- glm(totmass.male.com ~ as.factor(yr)-1, data = subset(d, nonzero == 1), family = Gamma(link = log))
# 

# set16 = read.csv(file.path(outdir,"4X_report_2016.csv"))
# set16[set16=='NULL'] = 0
# set16$yr=2016
# set16$lon = -set16$longitude
# set16$lat = set16$latitude
# set16$totmass.male.com=as.numeric(set16$mature.male..kg.)/set16$surace_area.km.2/1000
# set16$totmass.male.ncom=as.numeric(set16$immature.male..kg.)/set16$surace_area.km.2/1000
# set16$totmass.female.mat=as.numeric(set16$mature.female..kg.)/set16$surace_area.km.2/1000
# set16$totmass.female.imm=as.numeric(set16$immature.female..kg.)/set16$surace_area.km.2/1000
# 
# d = rbind(d,set16[filter.region.polygon(set16,'cfa4x'),c('yr','totmass.male.com',"totmass.male.ncom","totmass.female.mat","totmass.female.imm")])
# 
# Year = 2002:(p$year.assessment-1)
# offset=0.1
# gmmm = sapply(Year,function(y){with(subset(d,yr==y),exp(mean(log(totmass.male.com+offset)))-offset)})
# gmim = sapply(Year,function(y){with(subset(d,yr==y),exp(mean(log(totmass.male.ncom+offset)))-offset)})
# gmmf = sapply(Year,function(y){with(subset(d,yr==y),exp(mean(log(totmass.female.mat+offset)))-offset)})
# gmif = sapply(Year,function(y){with(subset(d,yr==y),exp(mean(log(totmass.female.imm+offset)))-offset)})
# 
# plot.new()
# plot(Year,gmmm,type='b',col='blue',xlab='Year',ylab='Geometric mean t / km^2',main='Male biomass',ylim=c(0,0.26),xlim=c(2002,2016),pch=16)
# lines(Year,gmim,type='b',col='blue',pch=17,lty=2)
# legend('topleft',c('mature','immature'),col=c('blue'),lty=1:2,pch=16:17)
# savePlot(file.path(outdir,paste('survey.trend.males','png',sep=".")),type='png')
# 
# 
# plot.new()
# plot(Year,gmmf,type='b',col='red',xlab='Year',ylab='Geometric mean t / km^2',main='Female biomass',ylim=c(0,0.26),xlim=c(2002,2016),pch=16)
# lines(Year,gmif,type='b',col='red',pch=17,lty=2)
# legend('topleft',c('mature','immature'),col=c('red'),lty=1:2,pch=16:17)
# savePlot(file.path(outdir,paste('survey.trend.females','png',sep=".")),type='png')
# 
# 
# 
# bub.ex=1
# bub.min=0.3
# leg=c(1,2,5,10)
# 
# bioMap(xlim=c(-66.5,-63.1),ylim=c(42.7,44.8),boundaries='snowcrab',main='2016')
# with(subset(set16,totmass.male.mat>0),points(X,Y,cex=bub.min+sqrt(totmass.male.mat)*bub.ex,pch=21,bg=rgb(0,1,0,0.2)))
# with(subset(set16,totmass.male.mat==0),points(X,Y,pch=4,cex=bub.min))
# legend('bottomleft',legend=c(0,leg),title= expression(t/km^2), pt.cex=c(bub.min,bub.min+sqrt(leg)*bub.ex),pch=c(4,rep(21,4)),pt.bg=rgb(0,1,0,0.2),bty='o',bg='white',box.col='white',inset=0.03,cex=1.2)
# savePlot(file.path(outdir,paste('surveyMMbubbles.2016','png',sep=".")),type='png')
# 
# set = snowcrab.db('set.biologicals')
# for (i in 2015:2002){
# 	plot.new()
# 	bioMap(xlim=c(-66.5,-63.1),ylim=c(42.7,44.8),boundaries='snowcrab',main=i)
# 	with(subset(set,lon<(-63)&yr==i&totmass.male.mat>0),points(lon,lat,cex=bub.min+sqrt(totmass.male.mat)*bub.ex,pch=21,bg=rgb(0,1,0,0.2)))
# 	with(subset(set,lon<(-63)&yr==i&totmass.male.mat==0),points(lon,lat,pch=4,cex=bub.min))
# 	legend('bottomleft',legend=c(0,leg),title= expression(t/km^2), pt.cex=c(bub.min,bub.min+sqrt(leg)*bub.ex),pch=c(4,rep(21,4)),pt.bg=rgb(0,1,0,0.2),bty='o',bg='white',box.col='white',inset=0.03,cex=1.2)
# 	savePlot(file.path(outdir,paste('surveyMMbubbles',i,'png',sep=".")),type='png')
# 
# }
# 
# 
# set = snowcrab.db('set.biologicals')
# for (i in 2015:2002){
# 	plot.new()
# 	bioMap('not4X',boundaries='snowcrab',main=i)
# 	with(subset(set,yr==i&totmass.male.mat>0),points(lon,lat,cex=bub.min+sqrt(totmass.male.mat)*bub.ex,pch=21,bg=rgb(0,1,0,0.2)))
# 	with(subset(set,yr==i&totmass.male.mat==0),points(lon,lat,pch=4,cex=bub.min))
# 	legend('bottomleft',legend=c(0,leg),title= expression(t/km^2), pt.cex=c(bub.min,bub.min+sqrt(leg)*bub.ex),pch=c(4,rep(21,4)),pt.bg=rgb(0,1,0,0.2),bty='o',bg='white',box.col='white',inset=0.03,cex=1.2)
# 	savePlot(file.path(outdir,paste('surveyMMbubblesNot4X',i,'png',sep=".")),type='png')
# 
# }
# 
# 
