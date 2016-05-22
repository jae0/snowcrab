
  figure.timeseries.gam.potential.habitat.survey.biomass= function( outdir=file.path(project.datadirectory('snowcrab'), "assessments","2014"), all.areas=T ) {
 
    set = snowcrab.db( DS="set.merge.det")
  
    if (all.areas) {
      areas = c("cfa4x", "cfasouth", "cfanorth" )
      regions = c("4X", "S-ENS", "N-ENS")
    } else {
      areas = c("cfasouth", "cfanorth" )
      regions = c("S-ENS", "N-ENS")
    }

    n.regions = length(regions)
    n.areas = length(areas)
    varnames = names(set)
#    vars = varnames[ grep ( vars, varnames) ]
fdir <-file.path( project.datadirectory("snowcrab"), "R", "gam","habitat" )
        fs <- dir(fdir)
    fs <- fs[setdiff(grep('K.R0.mass',fs) , grep('environmentals.only',fs))]
        lo <- c()
        for(i in 1:length(fs)) {
          load(file.path(fdir,fs[i]))
        lo <- rbind(lo,K)
        rm(K)
        }
        
        areas=c( "cfa4x","cfasouth", "cfanorth")
        
        td = lo[ which( lo$region %in% areas & lo$yr>1995) , c('yr','region','sa.region') ]
        td1 =data.frame(yr=rep(2013,3),region=c('cfa4x','cfanorth','cfasouth'),sa.region=aggregate(sa.region~region,data=td[which(td$yr>2008),],FUN=mean)[2])
        
        td = rbind(td,td1)
         td1 =data.frame(yr=rep(2014,3),region=c('cfa4x','cfanorth','cfasouth'),sa.region=aggregate(sa.region~region,data=td[which(td$yr>2009),],FUN=mean)[2])
        td = rbind(td,td1)
      
        td$year =td$yr
        td$sa.region = td$sa.region/1000
        load(file.path(project.datadirectory('snowcrab'),"R","ts.rdata"))
    ts <- ts[which(ts$variable=='R0.mass' & ts$region %in% areas),c('year','region','mean','ub','lb')]
      td <- merge(ts,td, by=c('year','region'),all.x=T)
    
         td = td[ order(td$region, td$year) , ]
      td$region = factor(td$region, levels=areas, labels =regions)
      td[which(td$mean<0),c('mean','ub','lb')] <- NA
      td$mean <- td$mean*td$sa.region
td$lb <- td$lb*td$sa.region
td$ub <- td$ub*td$sa.region


      ylim=range(c(td$mean+td$ub), na.rm=T); ylim[1]=ylim[1]-0.1*ylim[2]; ylim[2] = ylim[2]+ylim[2]*0.2
      xlim=range(td$year); xlim[1]=xlim[1]-0.5; xlim[2]=xlim[2]+0.5

      dir.create( outdir, recursive=T, showWarnings=F )
      fn = file.path( outdir, paste( 'fishablebiomasssurveygam', sep="." ) )
  
        cex.main = 1.4
    cex.lab = 1.3
    cex.axis = 1.3

     Cairo( file=fn, type="pdf", bg="white",  units="in", width=6, height=8 )
    setup.lattice.options()
    pl = xyplot( mean~year|region, data=td, ub=td$ub, lb=td$lb,
          layout=c(1,n.regions), xlim=xlim, 
              main="Fishable Biomass", xlab="Year", ylab=list("Fishable Biomass kt"),
              cex.lab=cex.lab, cex.axis=cex.axis, cex.main = cex.main,scales = list(y = "free"),
              panel = function(x, y, subscripts, ub, lb, ...) {
             larrows(x, lb[subscripts],
                     x, ub[subscripts],
                     angle = 90, code = 3, length=0.05)
             panel.xyplot(x, y, type="b", lty=1, lwd=2, pch=20, col="black", ...)
             panel.abline(h=median(y), col="gray75", ...)
         }
      )
   
       
      print(pl)
      dev.off()
    cmd( "convert   -trim -quality 9  -geometry 200% -frame 2% -mattecolor white -antialias ", paste(fn, "pdf", sep="."),  paste(fn, "png", sep=".") )
   
    }



