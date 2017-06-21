
#Welcome to BCFC. Ben’s Clunky Fucking Code! It’s ugly but it works.
  
  require(bio.base)
  
  if (!exists("current.year")) current.year=year(Sys.Date())
  p = bio.snowcrab::load.environment( year.assessment=current.year)

  fp = file.path( p$annual.results, "BZfigures")


#Figure 37
 
#This is a panel plot of maps where Bitter Crab has been observed. You can easily create this yourself using level plot. There is a field in the SNCRABDETAILS view that is called BCD. If a crab is bitter it has a value of 1. Otherwise it is an NA. So all you need to do is plot all survey stations in black and then plot and stations which have a bitter crab over the top in red. Let me know if this is confusing as hell.
 
#Figure 39
 
# Calculate the percentage of spring landings for each area by year
# You will likely have to change the dataset from logs to whatever your is called
# need to ensure that the data frame has year, cfa (aka cfa0) and Quarter of the year (q)

	logs = logbook.db(DS="logbook" )
	logs$q = quarter(logs$date.landed)
 
    compute.sums = function (x, var, index) {
    for (i in index) x[,i] = as.factor(x[,i])
    res = as.data.frame.table( tapply( X=x[,var], INDEX=x[,index],
                FUN=function(q) { sum(q, na.rm=T)}, simplify=T))
    names(res) = c(index, var)
    for (i in index) { res[,i] = as.character( res[,i] ) }
    return(res) }
 
    lbyq = compute.sums( x=logs, var="landings", index=c("q", "year", "cfa0")  )
   
    
     names(lbyq)=c("q", "year", "cfa", "kg")
     lbyq=lbyq[(as.numeric(lbyq$year)>2006),]
     lbyq=lbyq[lbyq$cfa %in% c("cfanorth", "cfa23", "cfa24"),]
     lbyq$kg[which(!is.finite(lbyq$kg))]=0
    
     lan=compute.sums(x=lbyq, var="kg", index=c("cfa", "year"))
    
     lbyq=merge(lbyq, lan, by=c("year", "cfa"))
     names(lbyq)=c("year", "cfa", "q", "kg", "total")
     lbyq$perc=NA
     lbyq$perc=(lbyq$kg/lbyq$total)*100
    
    
     spring=lbyq[lbyq$q==2,]
     spring$sorter=NA
     spring$sorter[spring$cfa=="cfanorth"]=1
     spring$sorter[spring$cfa=="cfa23"]=2
     spring$sorter[spring$cfa=="cfa24"]=3
         
      require(doBy)
      spring=orderBy(~cfa, spring)
     
     spring=spring[order(spring$cfa),]
      
    
 
# Save Figure
graphic='pdf'
 
filename=file.path(fp,paste("percent_spring_landings", graphic, sep="."))
require(devEMF)
  if(graphic=='emf')emf(file=filename, bg="white")
  if(graphic=='pdf')pdf(file=filename, bg="white")
        cf=unique(spring$cfa)
        cols=c("red", "green", "black")
        point=c(1,2,4)
        
        plot(spring$year, spring$perc, type="n", main= "Spring Landings", cex.main=1.8,
        ylab="Percent of Total", xlab="Year", ylim=c(0,100), lty=1, col="red", pch=19, bg="white")
        
        for (y in 1:length(cf)) {
         c=spring[spring$cfa==cf[1],]
         c=spring[spring$cfa==cf[y],]
         lines(x=c$year, y=c$perc, col=cols[y], lwd=2 )
         points(x=c$year, y=c$perc, col=cols[y], pch=point[y])
       
        }
        legend("bottomright",paste(cf), bty="n", col=cols, lty=1, pch=point)
 
  dev.off()
print(paste("Find file here: ", filename))
 
  
#Figure 41: Vessels Active
 
#Modify the following to fit
for (y in yrs){
 
vessels=logs[logs$year==y,]
 
compute.vessels = function (x, var, index) {
for (i in index) x[,i] = as.factor(x[,i])
res = as.data.frame.table( by( data=x, INDICES=x[,index],
                FUN=function(q) {
                                length(unique(q$vr_number))
                    } ) )
names(res) = c(index, var)
for (i in index) { res[,i] = as.character( res[,i] ) }
return(res) }
 
 
boats = compute.vessels( x=logs, var="vessels", index=c("year", "cfa0")  )
boats=boats[is.finite(boats$vessels) & boats$year>2004,]
 
areas=unique(boats$cfa0)
cols = c("red", "green4", "black", "cornflowerblue")
point=c(16, 17, 8, 15)
 
lim=c(1:4)
iye=max(boats$year[boats$cfa0=="cfa23"])
 
 
 
filename=paste(iye,"_vessels_per_year", ".emf", sep="")
        win.metafile(file=filename, width=10)
        plot(boats$year, boats$vessels, type="n", ylim=c(0,(max(boats$vessels)+1)), ylab="Vessels",
        main="Vessels Active by Year", xlab="Year" )
        for (l in lim){
        q = which(boats$cfa0 == areas[l] )
        lines(boats$year[q],boats$vessels[q], col=cols[l],lty=1, pch=point[l], type="o" )
        }
        legend("topright",paste(areas), bty="n", col=cols, pch=point, lty=1, ncol=2)
dev.off()
print(paste("Find file here: ", wd, filename,sep="/"))
 
filename=paste("vessels_per_year", ".pdf", sep="")
        pdf(file=filename)
        plot(boats$year, boats$vessels, type="n", ylim=c(0,(max(boats$vessels)+1)), ylab="Vessels",
        main="Vessels Active by Year", xlab="Year" )
        for (l in lim){
        q = which(boats$cfa0 == areas[l] )
        lines(boats$year[q],boats$vessels[q], col=cols[l],lty=1, pch=point[l], type="o" )
        }
        legend("topright",paste(areas), bty="n", col=cols, pch=point, lty=1, ncol=2)
dev.off()
print(paste("Find file here: ", wd, filename,sep="/"))
 
Figure 43- Weekly Smoothed Catch Rates
 
#Following script calculates weekly catch rates by area (in lbs)
#uses a three week running average of lbs and traps to smooth
#########################################
 
cleanlogsna=logs[! is.na(logs$num_of_traps),]
cleanlogscr=cleanlogsna[cleanlogsna$lbspertrap < 800,]
logs.fixed = cleanlogscr
rm(cleanlogscr)
rm(cleanlogsna)
 
logs.fixed$julian = as.integer(julian(logs.fixed$date_landed))
logs.fixed$dayofseason= ((logs.fixed$julian)-(min(logs.fixed$julian)))+1
logs.fixed$weekofseason=ceiling(logs.fixed$dayofseason/7)
 
nweek <- function(x, format="%Y-%m-%d", origin){
        if(missing(origin)){
                as.integer(format(strptime(x, format=format), "%W"))
        }else{
                x <- as.Date(x, format=format)
                o <- as.Date(origin, format=format)
                w <- as.integer(format(strptime(x, format=format), "%w"))
                2 + as.integer(x - o - w) %/% 7
        }
}
 
logs.fixed$weekofyear=NA
logs.fixed$weekofyear=nweek(logs.fixed$date_fished)
logs.fixed$weekofyear=factor(logs.fixed$weekofyear,levels=c(40:52, 1:39),ordered=TRUE)
logs.fixed$cfa0 = as.character(logs.fixed$cfa0)
areas = unique(logs.fixed$cfa0)
 
 
weekly=as.data.frame(xtabs(logs.fixed$pro_rated_slip_wt_lbs~logs.fixed$weekofyear+logs.fixed$cfa0+logs.fixed$year), stringsAsFactors=F)
effort=as.data.frame(xtabs(logs.fixed$num_of_traps~logs.fixed$weekofyear+logs.fixed$cfa0+logs.fixed$year),stringsAsFactors=F)
 
names(weekly)=c("week", "area", "year", "tot_lbs")
names(effort)=c("week", "area", "year", "tot_traps")
 
wk=merge(weekly, effort)
 
wk$run_lbs=NA
wk$run_traps=NA
wk$week=as.numeric(wk$week)
 
lags =c(-1,0,1)
areas=unique(wk$area)
 
for (j in 1:nrow(wk)){
   w=wk$week[j]
   y=wk$year[j]
   a=wk$area[j]
   wks=w+lags
   i=which(wk$week %in% wks & wk$year==y & wk$area==a)
   wk$run_lbs[j]=sum(wk$tot_lbs[i])
   wk$run_traps[j]=sum(wk$tot_traps[i])
    }
 
wk$cpue=wk$run_lbs/wk$run_traps
wk$cpuekg=wk$cpue/2.204626
wk$time=as.numeric((as.numeric(wk$year)+ (wk$week/52-.0001)))
wk$area[which(wk$area=="north")]="N-ENS"
wk$area[which(wk$area=="cfa23")]="CFA 23"
wk$area[which(wk$area=="cfa24")]="CFA 24"
wk$area[which(wk$area=="cfa4X")]="4X"
 
#Separate out last 3 years
back=c(-2,-1, 0)
past=as.character(max(as.numeric(logs.fixed$year)[which(logs.fixed$cfa0=="cfa23")])+ back)
i=which(wk$year %in% past)
recent=wk[i,]
recent=recent[is.finite(recent$cpue),]
cf=unique(recent$area)
non4x=recent[recent$area %in% c("CFA 23", "CFA 24", "N-ENS"),]
fc=unique(non4x$area)
 
cols=c("red","green","black", "blue")
point=c(16, 17, 8, 15)
 
smoothcatch=function(x){
     plot(recent$time,recent$cpuekg , type="n", main= "Weekly CPUE", cex.main=1.8, ylab="Catch Rate (kg/ trap)",
     xlab="Week", ylim=c(0,(max(recent$cpuekg))), xaxp=c(min(as.numeric(past)), max(as.numeric(past)), 2))
     legend("topright",paste(cf), bty="n", col=cols, pch=point)
   for (y in 1:length(cf)) {
         c=recent[recent$area==cf[y],]
         c=c[order(c$time),]
         c=c[c$run_lbs>quantile(c$run_lbs,.05),]
         points(x=c$time, y=c$cpuekg, col=cols[y], pch=point[y], cex=0.6)
   for (y in 1:length(fc) )
     {for (p in past){     
       c=recent[recent$area==fc[y],]
       c=c[order(c$time),]
       c=c[c$run_lbs>quantile(c$run_lbs,.05),]   
       this=c[c$year==p,]
       lines(x=this$time, y=this$cpuekg, col=cols[y+1], cex=0.6)
       }}}}
   
#Produce emf file
require(devEMF)
filename=paste("weekly_cpue_smoothed.emf", sep="")
emf(file=filename, bg="white")
smoothcatch(recent)
dev.off()
print(paste("Find file here: ", wd,"/", filename,sep=""))
#Produce pdf file
filename=paste("weekly_cpue_smoothed.pdf", sep="")
pdf(file=filename, width=10)
smoothcatch(recent)
dev.off()
print(paste("Find file here: ", wd,"/", filename,sep=""))
 
#Figure 47: Soft by Month
 
#----Full disclosure- All the tapplies and merges in the middle section are ugly…always mean to clean this up…haven’t..so you can clean or run as is….
 
# Use below to plot percent soft by month by year by cfa
 
require (chron)
require(lattice)
 
  connect=odbcConnect("BANK.CANSO3", uid="snowcrab", pwd="opilio99")
setsobs=sqlQuery(connect, ("SELECT * FROM SNCRABSETS_OBS, SNCRABDETAILS_OBS
WHERE SNCRABDETAILS_OBS.TRIP_ID = SNCRABSETS_OBS.TRIP_ID AND
    SNCRABDETAILS_OBS.SET_NO = SNCRABSETS_OBS.SET_NO"))
a=setsobs
 
# setting "hardlines" below lets you switch between multiple durometer levels to be considered hard
hl=68
 
 
# --------------------------------------
# convert lat's and long's to recognizable format for recode.areas
a=setsobs
h=names(a)
h[h=="LATITUDE"] = "lat"
h[h=="LONGITUDE"] = "lon"
names(a) = h
a$lon=-a$lon
 
#if no landing date, populate with board date
a$LANDING_DATE[is.na(a$LANDING_DATE)]=a$BOARD_DATE[is.na(a$LANDING_DATE)]
 
 
# Add a column for year
#--------------------------------------------
a$year=NA
a$year=as.character(years(as.chron(a$LANDING_DATE)))
 
a$date=a$LANDING_DATE
 
# Add month
#--------------------------------------------
a$month=NA
a$month=as.character(months(as.chron(a$LANDING_DATE)))
 
# Add season (Q1 =winter, Q2=spring........)
#--------------------------------------------
a$season=NA
a$season=as.character(quarters(a$date))
a$season[which(a$season=="Q1")]="winter"
a$season[which(a$season=="Q2")]="spring"
a$season[which(a$season=="Q3")]="summer"
a$season[which(a$season=="Q4")]="fall"
 
 
# --------------------------------------
# create a column unique to each set called tripset
 
a$tripset= paste(a$TRIP,a$SET_NO,sep="~")
a$kgpertrap= (a$EST_CATCH*1000)/(a$NUM_HOOK_HAUL)
a=a[!is.na(a$kgpertrap),]
 
 
 
# --------------------------------------
# create field for soft/hard
 
a= a[is.finite(a$DUROMETRE),]
a$dummy= 1
a$hardness= ifelse(a$DUROMETRE>=hl, 1, 0)
 
# --------------------------------------
# calculate counts of hard and total crab for each trip set
 
b=xtabs(as.integer(a$hardness)~as.factor(a$tripset))
b=as.data.frame(b)
names(b)= c("tripset", "no.hard")
c=xtabs(as.integer(a$dummy)~as.factor(a$tripset))
c=as.data.frame(c)
names(c)= c("tripset", "total")
 
# --------------------------------------
# merge two new data frames
 
d=merge(x=b, y=c, by="tripset", all.x=F, all.y=T, sort=F)
d$percenthard=(d$no.hard/d$total)*100
d$percentsoft=100-d$percenthard
 
 
# --------------------------------------
# import catch rate from a
 
e=tapply(X=a$kgpertrap, INDEX= a$tripset, FUN= function(q){unique(q)[1]})
e=as.data.frame(e)
f=data.frame(kgpertrap=as.vector(e$e), tripset=dimnames(e)[[1]] )
f$tripset =  as.character(f$tripset)
 
g=merge(x=d, y=f, by="tripset", all.x=F, all.y=T, sort=F)
 
# --------------------------------------
# year from a
 
 
ee=tapply(X=a$year, INDEX= a$tripset, FUN= function(q){unique(q)[1]})
ee=as.data.frame(ee)
ff=data.frame(year=as.vector(ee$ee), tripset=dimnames(ee)[[1]] )
ff$tripset =  as.character(ff$tripset)
ff$year = as.character(ff$year)
 
 
gg=merge(x=g, y=ff, by="tripset", all.x=F, all.y=T, sort=F)
 
 
# --------------------------------------
# import latitude from a
 
h=tapply(X=a$lat, INDEX= a$tripset, FUN= function(q){unique(q)[1]})
h=as.data.frame(h)
i=data.frame(lat=as.vector(h$h), tripset=dimnames(h)[[1]] )
i$tripset =  as.character(i$tripset)
 
j=merge(x=gg, y=i, by="tripset", all.x=F, all.y=T, sort=F)
 
# --------------------------------------
# import longitude from a
 
k=tapply(X=a$lon, INDEX= a$tripset, FUN= function(q){unique(q)[1]})
k=as.data.frame(k)
l=data.frame(lon=as.vector(k$k), tripset=dimnames(k)[[1]] )
l$tripset =  as.character(l$tripset)
 
mm=merge(x=j, y=l, by="tripset", all.x=F, all.y=T, sort=F)
 
# --------------------------------------
# import season from a
 
kk=tapply(X=a$season, INDEX= a$tripset, FUN= function(q){unique(q)[1]})
kk=as.data.frame(kk)
ll=data.frame(season=as.vector(kk$kk), tripset=dimnames(kk)[[1]] )
ll$tripset =  as.character(ll$tripset)
 
mmm=merge(x=mm, y=ll, by="tripset", all.x=F, all.y=T, sort=F)
 
# --------------------------------------
# import month from a
 
kkk=tapply(X=a$month, INDEX= a$tripset, FUN= function(q){unique(q)[1]})
kkk=as.data.frame(kkk)
lll=data.frame(month=as.vector(kkk$kkk), tripset=dimnames(kkk)[[1]] )
lll$tripset =  as.character(lll$tripset)
 
m=merge(x=mmm, y=lll, by="tripset", all.x=F, all.y=T, sort=F)
# --------------------------------------
# import estimated catch from a
 
n=tapply(X=a$EST_CATCH, INDEX= a$tripset, FUN= function(q){unique(q)[1]})
n=as.data.frame(n)
p=data.frame(EST_CATCH=as.vector(n$n), tripset=dimnames(n)[[1]] )
p$tripset =  as.character(p$tripset)
 
x=merge(x=m, y=p, by="tripset", all.x=F, all.y=T, sort=F)
 
x$discardkg=(x$EST_CATCH*1000)*(x$percentsoft/100)
 
 
# --------------------------------------
# divide into north, south, and 4X by positions
 
source("C:/Scripts/geo.filter.r")
 
cfa=c("cfanorth", "cfa23", "cfa24", "cfa4x")
 
x$cfa=NA
for  (a in cfa){
     rowindex= filter.region.polygon(x,recode.areas(a))
     x$cfa[rowindex]=a
}
 
area=c("cfanorth", "cfasouth", "cfa4x")
 
x$area=NA
for (a in area){
rowindex= filter.region.polygon(x,recode.areas(a))
x$area[rowindex]=a
}
 
#4X fishing activity on 4X line doesn't get recoded properly by recode.areas
#force thes into 4X
x$long=NA
x$long=-(as.numeric(x$lon))
i4x=which(x$month %in% c("Nov", "Dec", "Jan", "Feb") & x$long>63.3)
x$cfa[i4x]="cfa4x"
 
#set up month and year as ordered factors
x$area=factor(x$area, levels=c("cfa4x","cfasouth","cfanorth"), labels=c("4X","South","North"))
x=x[!is.na(x$area),]
 
x$cfa=factor(x$cfa, levels=c("cfa4x","cfa24","cfa23","cfanorth"), labels=c("4X","CFA 24","CFA 23","North"))
 
####### Need to resolve month script to become an ordered factor
x$month=factor(x$month, levels=c("Jan","Feb","Mar","Apr","May", "Jun", "Jul", "Aug", "Sep","Oct","Nov", "Dec"),
               labels=c("1","2","3","4","5","6", "7","8", "9", "10", "11", "12"))
 
yrs=as.character(sort(as.numeric(unique(x$year))))
x$year=factor(x$year,levels=yrs, labels=yrs)
 
x$areaseason=NA
x$areaseason=paste(x$area ,":", x$season)
 
xnorth=x[x$area=="North", ]
xsouth=x[x$area=="South",]
x23=x[x$cfa=="CFA 23",]
x24=x[x$cfa=="CFA 24",]
 
jj=x
 
 
softbymonth=xtabs(percentsoft~month+year+cfa, data=x)
#softbymonth
 
#-----------------------------------------------
# Use below to plot percent soft by month by year by cfa
 
compute.sums = function (x, var, index) {
for (i in index) x[,i] = as.factor(x[,i])
res = as.data.frame.table( tapply( X=x[,var], INDEX=x[,index],
FUN=function(q) { sum(q, na.rm=T)}, simplify=T))
names(res) = c(index, var)
for (i in index) { res[,i] = as.character( res[,i] ) }
return(res) }
 
z=compute.sums(x=x, var=c("discardkg"), index=c("month","year","cfa"))
z=na.omit(z)
 
v=compute.sums(x=x, var=c("EST_CATCH"), index=c("month","year","cfa"))
v=na.omit(v)
 
w=merge(x=z, y=v, by=c("month","year","cfa"), all.x=T, all.y=F, sort=F)
 
w$totalkg=NA
w$totalkg=w$EST_CATCH*1000
w$percsoft=NA
w$percsoft=(w$discardkg/w$totalkg)*100
w$month=as.numeric(as.character(w$month))
w$year=as.numeric(as.character(w$year))
 
baddata=which(w$cfa=="4X" & w$month %in% c(5,6,7,8,9,10))
w$bad=1
w$bad[baddata]=0
w=w[w$bad>0,]
 
 
#-----------------------------------------------
# Use below to plot soft by month for years for all areas by month
 
w=w[w$year>(iy-5),]
w$year=factor(w$year)
 
 
#Save as metafile
filename=paste("soft_crab_by_month.emf", sep="")
win.metafile(file=filename, width=10)
setup.lattice.options = function () {
      trellis.par.set(col.whitebg())
      s.bkg = trellis.par.get("strip.background")
      s.bkg$col="gray95"
      trellis.par.set("strip.background", s.bkg)
    }
        default.options = options()
  setup.lattice.options()
 
xyplot(
    percsoft~month|year+cfa, data=w,  main=paste("Percent Soft Shell by Month", sep=""),
    ylim=c(0, 100), xlab="Month", ylab="Percent Soft",
        panel = function(x, y, ...) {
       panel.xyplot(x, y, type="p",  pch=20,
col="red", ...)
panel.abline(h=20, col="red", lty=1, lwd=1, ...) }
)
options(default.options)
 
dev.off()
 
#save as pdf
filename=paste("soft_crab_by_month.pdf", sep="")
pdf(file=filename)
setup.lattice.options = function () {
      trellis.par.set(col.whitebg())
      s.bkg = trellis.par.get("strip.background")
      s.bkg$col="gray95"
      trellis.par.set("strip.background", s.bkg)
    }
        default.options = options()
  setup.lattice.options()
 
xyplot(
    percsoft~month|year+cfa, data=w,  main=paste("Percent Soft Shell by Month", sep="")
    ,ylim=c(0, 100), xlab="Month", ylab="Percent Soft",
        panel = function(x, y, ...) {
       panel.xyplot(x, y, type="p",  pch=20,
col="red", ...)
panel.abline(h=20, col="red", lty=1, lwd=1, ...) }
)
options(default.options)
 
dev.off()
 
 
