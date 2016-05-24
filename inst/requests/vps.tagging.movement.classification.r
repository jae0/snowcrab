#Emera Tagging
        
loadfunctions( "snowcrab", functionname="initialise.local.environment.r") 

dn = file.path(project.datadirectory('snowcrab'),'data','tagging','Emera')

if(redo.base.data) {} 
        a = dir(dn)
        a = a[grep('crab',a)] 

        out = NULL
                for(i in a) {
                         h = read.csv(file.path(dn,i),header=T)
                         out = rbind(h,out)
                        }
         
        names(out) = tolower(names(out))
        out$dt = strptime(out$datetime, format= "%m/%d/%Y %H:%M")
         
        #N Animals detects
                        b= unique(out$transmitter) #58
                         lonlat2planar(oo)
        #N Detects per anmimal
                        a = aggregate(detectedid~transmitter,data=out,FUN=length)
         

        #each animal Filter to reasonable speeds 
        filt.data = NULL
        #hpe is the error of the data points <8 is good #need to implement
        for(i in 1:length(b)) {
                                oo = out[which(out$transmitter==b[i]),]
                        oo = oo[order(oo$dt),]
                        dt=1
                        thres = 15
                        oo1 = oo
                        while(dt==1) {
                                oo1$dist = c(NA,distVincentyEllipsoid(p1=oo1[1:(nrow(oo1)-1),c('lon','lat')],p2=oo1[2:(nrow(oo1)),c('lon','lat')]))
                                oo1$time = c(NA,oo1[2:(nrow(oo1)),'dt']-oo1[1:(nrow(oo1)-1),'dt'])
                                oo1$speed  = oo1$dist/oo1$time # m / min
                                if(any(oo1$speed>thres,na.rm=T)) {
                                        io = which.max(oo1$speed)
                                        oo1 = oo1[-io,]
                                } else {
                                 dt=0       
                                }
                            }
                hist(oo1$speed,'fd',main=b[i])
        savePlot(file.path(dn,paste('hist.of.speed',b[i],"png",sep=".")),type='png')
                                        
        filt.data = rbind(filt.data,oo1)                
        loadfunctions(c('utility','polygons'))
        xx = range(oo$lon)
        yy = range(oo$lat)
        xx[1] = xx[1]-0.02; xx[2] = xx[2] +0.02
        yy[1] = yy[1]-0.02; yy[2] = yy[2] +0.02

        makeMap(xlim=xx,ylim=yy,addEmera=T,addSummerStrata=F)
                        with(oo,lines(lon,lat,type='l',lwd=0.6))
                        with(oo1,lines(lon,lat,type='l',lwd=0.6,col='green'))

        title(paste(b[i],'N days =', round(as.numeric(max(oo$dt)-min(oo$dt)[1]),1)        ))
        savePlot(file.path(dn,paste('tag',b[i],"png",sep=".")),type='png')
            }
     save(filt.data,file=file.path(dn,'filtered.movement.data.2015.rdata'))
        }
load(file=file.path(dn,'filtered.movement.data.2015.rdata'))



#Cumulative weighted distribution of movements
plot(1,1,xlim=c(0,15),ylim=c(0,1),type='n')
filt.data$round.speed = round(filt.data$speed*100)/100
ff = aggregate(time~round.speed+transmitter,data=filt.data,FUN=length)
ui = unique(ff$transmitter)
for(j in ui){
        w = ff[which(ff$transmitter == j),]
        w$cspeed = cumsum(w$time) / sum(w$time)
        lines(w$round.speed,w$cspeed,col=1)
}
hist()

#build a polygon covering the data
 ch = filt.data[chull(filt.data$lon,filt.data$lat),c('lon','lat')]
 fd = makePBS(filt.data[,c(-4,-5)],polygon=F)
 bdpoly = data.frame(PID=1,POS=1:4,X=c(min(ch$lon),min(ch$lon),max(ch$lon),max(ch$lon)),Y=c(min(ch$lat),max(ch$lat),max(ch$lat),min(ch$lat)))
 pp = data.frame(PID=rep(1:100,each=4),POS=rep(1:4,times=100),Y=rep(c(bdpoly[1,'Y'],bdpoly[2,'Y'],bdpoly[2,'Y'],bdpoly[1,'Y']),times=100),X=c(min(ch$lon),min(ch$lon)),Xinc = c(0,0, rep(0.00023392,198)),XX = c(0,0,rep(1:99,each=4),100,100))
pp$X = pp$X+pp$Xinc*pp$XX

a = findPolys(fd,pp,maxRows=1300000)
 ff=  merge(fd,a,by='EID')
 a = unique(cbind(ff$transmitter,ff$PID))
 a = aggregate(transmitter~PID,data=a,FUN=length)
 barplot(a$transmitter)


