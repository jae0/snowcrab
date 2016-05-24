loadfunctions( "snowcrab", functionname="initialise.local.environment.r") 


makeMap(area='NENS',addSurvey=T,addEmera=T)
makeMap(addSurvey=T)
 
  
 areas =  "cfanorth" 
 regions =  "N-ENS"
    
    n.regions = length(regions)
    n.areas = length(areas)
    varnames = names(set)
#    vars = varnames[ grep ( vars, varnames) ]
fdir <-file.path( project.directory("snowcrab"), "R", "gam","habitat" )
        fs <- dir(fdir)
    fs <- fs[setdiff(grep('K.R0.mass',fs) , grep('environmentals.only',fs))]
        lo <- c()
        for(i in 1:length(fs)) {
          load(file.path(fdir,fs[i]))
        lo <- rbind(lo,K)
        rm(K)
        }
        
        areas=c("cfanorth")
        
        td = lo[ which( lo$region %in% areas & lo$yr>1995) , c('yr','region','sa.region') ]
        td1 =data.frame(yr=rep(2013,3),region=c('cfa4x','cfanorth','cfasouth'),sa.region=aggregate(sa.region~region,data=td[which(td$yr>2008),],FUN=mean)[2])
        
        td = rbind(td,td1)
        td$year =td$yr
        td$sa.region = td$sa.region/1000
        load(file.path(project.directory('snowcrab'),"R","ts.rdata"))
    ts1 <- ts[which(ts$variable=='totno.male.com' & ts$region %in% areas & ts$year > 2003),c('year','region','mean','ub','lb')]
      td1 <- merge(ts1,td, by=c('year','region'),all.x=T)
    
     ylim=c(0,max(c(td1$mean+td1$ub,td1$mean-td1$lb), na.rm=T))
 
pdf('NumberCommercial.pdf')    
with(td1,plot(year,mean,ylim=ylim,ylab='Number Commerical per km2',pch=16))
with(td1,arrows(x0=year,y0=mean-lb,y1=mean+ub,length=0))
dev.off()



set = snowcrab.db( DS="set.merge.det")
set = set[which(set$station %in% c(173,184,188,612,3,703)),]

a=merge(aggregate(totno.male+totno.female~yr,data=set,FUN= mean),aggregate(totno.male+totno.female~yr,data=set,FUN= sd),by='yr')

yl = c(0,max(a[,2]+a[,3]))
pdf('TotalNsubset.pdf')
plot(a$yr,a[,2],pch=16,ylab='Number per km2',ylim=yl,main='Subset Stations')
arrows(x0=a$yr,y0=a[,2]-a[,3],y1=a[,2]+a[,3],length=0)
dev.off()