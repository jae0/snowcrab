load("/home/ben/bio.data/bio.snowcrab/data/set.complete.rdata")
set$P = 1
require(PBSmapping)

a = aggregate(P~cfa+yr,subset(set,totmass.male.com>0),FUN=sum)
b = aggregate(P~cfa+yr,subset(set),FUN=sum)
g=merge(a,b,by=c('cfa','yr'))
g$prop = g$P.x/g$P.y
nN = list()
m=0
for(i in 2004:2018){
      m = m+1 
       nor=set[which(set$cfa=='cfanorth'&set$yr==i),]
      #AREA OF NORTH BASED ON CONVEX HULL OF SURVEY
           nor$x = nor$lon
            nor$y = nor$lat
            norcord = nor[chull(nor$x,nor$y),c('x','y')]
            plot(norcord)
            norcord$X = norcord$x;
            norcord$Y = norcord$y
            norcord$PID = 1
            norcord$POS = 1:nrow(norcord)
            plotPolys(norcord)
            attr(norcord,'projection') <- "LL"
            nN[[m]] = c(i, calcArea(norcord))
}
  nN = as.data.frame(do.call(rbind,nN))
  mean(unlist(nN$area) )
      #AREA OF NORTH 8349
gN = subset(g,cfa=='cfanorth')
gN$Area = gN$prop *8349 #PROXY FOR HABITAT

require(PBSmapping)

###SOUTH
sE = list()
sW = list()
m=0
for(i in 2004:2018) {
  m=m+1
        nor=set[which(set$cfa=='cfasouth' & set$lon>-62 & set$yr==i),]
        #AREA OF NORTH BASED ON CONVEX HULL OF SURVEY
        nor$x = nor$lon
        nor$y = nor$lat
        norcord = nor[chull(nor$x,nor$y),c('x','y')]
        plot(norcord)
        norcord$X = norcord$x;
        norcord$Y = norcord$y
        norcord$PID = 1
        norcord$POS = 1:nrow(norcord)
        plotPolys(norcord)
        attr(norcord,'projection') <- "LL"
        sE[[m]] = c(i,calcArea(norcord))
        
        nor=set[which(set$cfa=='cfasouth' & set$lon< -62& set$yr==i),]
        nor$x = nor$lon
        nor$y = nor$lat
        norcord = nor[chull(nor$x,nor$y),c('x','y')]
        plot(norcord)
        norcord$X = norcord$x;
        norcord$Y = norcord$y
        norcord$PID = 1
        norcord$POS = 1:nrow(norcord)
        plotPolys(norcord)
        attr(norcord,'projection') <- "LL"
        sW[[m]] = c(i,calcArea(norcord))
}
  sE = as.data.frame(do.call(rbind,sE))
  sW = as.data.frame(do.call(rbind,sW))
  sS = merge(sE,sW,by=c('V1','PID'))
  sS$area = unlist(sS$area.x)+unlist(sS$area.y)
  mean(sS$area)
#AREA OF South 71190
gS = subset(g,cfa=='cfasouth')
gS$Area = gS$prop *71190 #PROXY FOR HABITAT

#4X
sX = list()
    m=0
for(i in c(2004:2013,2015:2018)) {
  m=m+1
        nor=set[which(set$cfa=='cfa4x' & set$yr==i),]
        #AREA OF NORTH BASED ON CONVEX HULL OF SURVEY
        require(PBSmapping)
        nor$x = nor$lon
        nor$y = nor$lat
        norcord = nor[chull(nor$x,nor$y),c('x','y')]
        plot(norcord)
        norcord$X = norcord$x;
        norcord$Y = norcord$y
        norcord$PID = 1
        norcord$POS = 1:nrow(norcord)
        plotPolys(norcord)
        attr(norcord,'projection') <- "LL"
        sX[[m]] = c(i,calcArea(norcord))
        }
        sX = as.data.frame(do.call(rbind,sX))
        mean(unlist(sX$area))
        #AREA OF 4X 5854
        gX = subset(g,cfa=='cfa4x')
        gX$Area = gX$prop *5854 #PROXY FOR HABITAT
        td = reshape(td,idvar='yr',timevar='region',direction='wide')
        
        Area = rbind(rbind(gX[,c('cfa','yr','Area')],gS[,c('cfa','yr','Area')]),gN[,c('cfa','yr','Area')])
        Area = reshape(Area,idvar='yr',timevar='cfa',direction='wide')
        Area = Area[order(Area$yr),]
        Area = subset(Area,yr>1998)
        Area$Area.cfa4x[which(Area$yr<2004)] <- 0
        
