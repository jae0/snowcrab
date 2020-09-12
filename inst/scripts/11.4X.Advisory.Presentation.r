#logbook 4X

year.assessment=2019

p = bio.snowcrab::load.environment( year.assessment=year.assessment )

outdir = file.path(project.datadirectory('bio.snowcrab'),'assessments',  p$year.assessment, "presentations", '4X')
dir.create(outdir,showWarnings=T, recursive = T)

warning( "This maping section does not like RStudio, run directly in R")

#Map the Area
    require(PBSmapping)
    require(SpatialHub)
    project.library ( 'stmv', 'aegis' )

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
	landings$tac=NA
	landings$tac[landings$yr=="2006"]=338
	landings$tac[landings$yr=="2007"]=230
	landings$tac[landings$yr=="2008"]=230
	landings$tac[landings$yr=="2009"]=230
	landings$tac[landings$yr=="2010"]=346
	landings$tac[landings$yr=="2011"]=346
	landings$tac[landings$yr=="2012"]=263
	landings$tac[landings$yr=="2013"]=80
	landings$tac[landings$yr=="2014"]=80
	landings$tac[landings$yr=="2015"]=150
	landings$tac[landings$yr=="2016"]=80
	landings$tac[landings$yr=="2017"]=110 #Need to Update list each year

		landings$landings=landings$landings/1000
	landplot=function(){
    plot(landings$landings, type="n", ylim=c(0,500), ylab="Landings (mt)", main="4X Landings by Year", xaxt="n", xlab="Year" )
	  axis(1, at=1:length(landings$landings), labels=landings$yr)
	  lines(landings$tac, col="blue", lty=2)
	  points(landings$landings, col="red", pch=20)
	  lines(landings$landings, col="red")
	  legend(x=1, y=80, c("TAC", "Landings"), lty= c(2, 1), col=c("blue", "red"))
	  }


	pdf(file=file.path(outdir,"annual.landings.pdf"))
	landplot()
	dev.off()

	require(car)
	#determine monthly landings by year
	monthly = aggregate(landings~month+yr,data=logs,FUN=sum)
	monthly$fm = recode(monthly$month,"'11'='November';'12'='December';'1'='January';'2'='February';'3'='March'; '4'='April'")
	mts= c("November", "December", "January", "February", "March")
	monthly$fm=factor(monthly$fm, levels=mts)
	monthly$landings[!is.finite(monthly$landings)]=0

	monthly$mt=monthly$landings/1000
	monthly = monthly[order(monthly$yr,monthly$fm),]

	#plot monthly landings
	#plot.new()
monthplot=function(){
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
}

pdf(file=file.path(outdir,"monthly.landings.pdf"))
monthplot()
dev.off()

	### Determine number of active vessel by month
	boats = aggregate(cfv~month+yr,data=unique(logs[,c('cfv','month','yr')]),FUN=length)
	mts= c("November", "December", "January", "February", "March")
	boats$fm = recode(boats$month,"'11'='November';'12'='December';'1'='January';'2'='February';'3'='March'")
	boats$fm=factor(boats$fm, levels=mts)
	boats$cfv[!is.finite(boats$cfv)]=0
	boats = boats[order(boats$yr,boats$fm),]
	boat.plot=function(){
	  ylims = c(0,7)
	  cols = c("red", "green4", "black",'blue')
		plot(1:5, 1:5,type="n", ylim=ylims, ylab="# of Vessels", main="4X Vessels Active By Month", xlab="Month" ,
	     xaxt='n')
	  axis(1, at=1:length(mts), labels=mts)
	  yu = as.numeric(max(as.numeric(boats$yr)))
	  ny=(yu-2):yu
	  for(y in 1:length(ny)) {
	  with(boats[boats$yr==ny[y],],points(fm,pch=17,cex=1.2, cfv,col=cols[y]))
	    with(boats[boats$yr==ny[y],],lines(fm,pch=16, cfv,col=cols[y]))
	  }
	legend('topright',legend=ny,col=cols,lty=rep(1,4),pch=rep(1,4),bty='n')
	}

	pdf(file=file.path(outdir,"annual.vessels.pdf"))
	boat.plot()
	dev.off()



#Plot annual EFFORT
	traps = with(logs,tapply(effort,yr,sum,na.rm=T))

	trap.plot=function(){
	  plot(traps, type="n", ylim=c(0,max(traps)), ylab="Trap Hauls", main="4X Trap Hauls by Year", xaxt="n", xlab="Year" )
	  axis(1, at=1:length(traps), labels=names(traps))
	  points(traps, col="red", pch=20)
	  lines(traps, col="red")
	}

	pdf(file=file.path(outdir,"annual.effort.pdf"))
	trap.plot()
	dev.off()


#Plot annual CPUE
#maybe jackknife this?
	cpue.catch = with(subset(logs,!is.na(effort)&!is.na(landings)),tapply(landings,yr,sum,na.rm=T))
	cpue.traps = with(subset(logs,!is.na(effort)&!is.na(landings)),tapply(effort,yr,sum,na.rm=T))

	cpue = cpue.catch/cpue.traps

	cpue.plot=function(){
	  plot(cpue, type="n", ylim=c(0,max(cpue)), ylab="Kg / Trap Haul", main="4X CPUE by Year", xaxt="n", xlab="Year" )
	  axis(1, at=1:length(cpue), labels=names(cpue))
	  points(cpue, col="red", pch=20)
  	lines(cpue, col="red")
	  abline(h=mean(cpue), col="blue", lty=3, lwd=1)
	}

	pdf(file=file.path(outdir,"cpue.pdf"))
	cpue.plot()
	dev.off()

#add jackknife estimates of error
	jack=logs[logs$cfa0=="cfa4x",] #create new dataframe (direct copy of logs)
	names(jack)[names(jack) == 'landings'] <- 'catch'
	names(jack)[names(jack) == 'cfa'] <- 'area'
	jack$yr=as.character(jack$yr)
	jack=jack[,c('yr','catch','effort','area')]
	jack=na.omit(jack)

	cpue = jackknifeCPUE(jack,grouping=c('yr'))

	jack.cpue.plot=function(){
	ylims = c(0,max(cpue$cpue+cpue$cpue.var)*1.2)
	plot(cpue$yr, cpue$cpue,type="n", ylim=ylims, ylab="Kgs / Trap", main="4X Catch Rates by Year", xlab="Year" )
	with(cpue,arrows(x0=as.numeric(yr),y0=cpue-cpue.var,y1=(cpue+cpue.var), col="red", length=0))
	with(cpue,points(yr,cpue, col="blue", pch=23,type='b'))
	}

	pdf(file=file.path(outdir,"annual.cpue.jack.pdf"))
	jack.cpue.plot()
	dev.off()

	rm(jack)



# Import Observer Data
#------------------------------------------------------------
#Need to do new data pull if first time being run since assessment (~7 minutes)
observer.db( DS="rawdata.redo", yrs=2003:p$year.assessment )

yu = p$year.assessment
yrs=(yu-4):yu

l = observer.db('rawdata',yrs=yrs)
l$tripset=paste(l$TRIP, l$SET_NO, sep=":")

# Remove CW's outside norms and remove production (pre-sorted) samples
a = l[l$FISH_LENGTH>50 & l$FISH_LENGTH<170,]

# convert lat's and long's to recognizable format for aegis.polygons::polygon_internal_code
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
a = a[polygon_inside(a,'cfa4x'),]
  yu = p$year.assessment
  yrs=(yu-8):(yu-1)



cc.hist.plot=function(y=y){
    x=a[a$yr==y,]
    avg=mean(x$FISH_LENGTH)

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

    # create stacked barplots with legends

    barplot(xplot [c(5:1),], space=0,names.arg=seq(50, 170, by=3)[-1],
     main=paste(y, as.character (as.numeric(y)+1), sep="/"), legend.text=c(paste("CC5 (",xcc5perc,"%)"),
     paste("CC4 (",xcc4perc,"%)"), paste("CC3 (",xcc3perc,"%)"), paste("CC2 (",xcc2perc,"%)"),
      paste("CC1 (",xcc1perc,"%)")), xlab="Carapace Width in mm", ylab="Number of Crab")
    abline(v=(95-50)/3, lty=2, lwd=2) #minimum legal size
    #abline(v=(avg-50)/3, lty=3, lwd=2, col="red") #minimum legal size
  }

#Need to increment years below. NB- starting year of season

pdf(file=file.path(outdir,"cc.obs.present.pdf"))
cc.hist.plot(y=2017)
dev.off()

pdf(file=file.path(outdir,"cc.obs.past.pdf"))
cc.hist.plot(y=2016)
dev.off()

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

      #determine number of observed trips
      trips= length(unique(j$TRIP))
      trips


     out=rbind(out,c(y, observedmt, sampledtraps, observedtraps, trips))
    }

  out=as.data.frame(out)
  b=out
  out=b

  names(out)=c("Year", "Observed", "Traps Sampled", "Traps Observed", "Trips")

  ys=yrs

  out$Landings=NA
  cfa=unique(out$Area)
  for (y in ys){
      out$Landings[out$Year==y]=round(landings$landings[landings$yr==y],0)
  }

  ### TODO BZ- Make sure trip get carried from last step to here


  out$Observed=round(as.numeric(as.character(out$Observed)))
  out$Percent=round(out$Observed/out$Landings*100,1)
  names(out)=c("Year", "Observed (mt)", "Traps Sampled", "Traps Observed","Trips", "Landings (mt)", "% Observed")

 output=out[,c(-2, -6)]

 require(gridExtra)
  obs.stats=function(){
   gridExtra::grid.table(output, theme=ttheme_default(), rows=NULL)
  }

  pdf(file=file.path(outdir,"observersummary.pdf"))
  obs.stats()
  dev.off()

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
set.plot=function(){
for (i in yrs){
      plot.new()
      bioMap(xlim=c(-66.5,-63.1),ylim=c(42.7,44.8),boundaries='snowcrab',main=i)
      with(subset(set,lon<(-63)&yr==i&totmass.male.mat>0),points(lon,lat,cex=bub.min+sqrt(totmass.male.mat)*bub.ex,pch=21,bg=rgb(0,1,0,0.2)))
      with(subset(set,lon<(-63)&yr==i&totmass.male.mat==0),points(lon,lat,pch=4,cex=bub.min))
      legend('bottomleft',legend=c(0,leg),title= expression(t/km^2), pt.cex=c(bub.min,bub.min+sqrt(leg)*bub.ex),pch=c(4,rep(21,4)),pt.bg=rgb(0,1,0,0.2),bty='o',bg='white',box.col='white',inset=0.03,cex=1.2)

  pdf(file=file.path(outdir,paste("survey.mat.male.bub",i,"pdf", sep=".")))
  set.plot()
  dev.off()
}

}
set.plot()

pdf(file=file.path(outdir,paste("survey.mat.male.bub",i,"pdf", sep=".")))
set.plot()
dev.off()


