figure.crab.groundfish.survey <- function(species , outdir=file.path(project.datadirectory('snowcrab'), 
	"assessments",p$current.assessment.year)) {
		require(latticeExtra)
    require(boot)
      areas.crab = c("cfa4x", "cfasouth", "cfanorth" )
      areas.groundfish = list(cfa4x = c(470,483),cfasouth = c(443,467),cfanorth=c(440,442))
      regions = c("4X", "S-ENS", "N-ENS")
      n.regions = length(regions)
      n.areas = length(areas.crab)
	      ft2m = 0.3048
	      m2km = 1/1000
	      nmi2mi = 1.1507794
	      mi2ft = 5280
	      sakm2 = (41 * ft2m * m2km ) * ( 1.75 * nmi2mi * mi2ft * ft2m * m2km )  # surface area sampled in km^2
		
#	      species = c( "cod", "skates", "thornyskate", "northernshrimp",
#	                   "brittlestar", "basketstar",
#	                   "snowcrab",  "hermitcrab", "jonahcrab",
#	                   "lessertoadcrab", "northernstonecrab", "porcupinestonecrab", "portlyspidercrab",
#	                   "redcrab", "atlanticrockcrab", "toadcrab", "snowcrab",'atlanticwolffish','halibut','smoothskate','winterskate')
	      #cod, halibut, thornyskate, wolfish, lessertoadcrab, jonahcrab, smoothskate, winterskate, northernshrimp, 
	      species = c(10, 30, 201, 50, 2521, 2511, 202, 204, 2211) 
 set = snowcrab.db( DS="set.complete")
 set = set[which(set$yr>2003),]
 ii  = names(set)[grep(paste("^ms.mass.",species,"$",sep=""),names(set))]
  
    if(length(ii)==1){
		 set = set[,c('yr','cfa',ii)]
		 names(set)[3] = 'mass'
		 ou = unique(set[,c('yr','cfa')])
		 ou = ou[order(ou$cfa,ou$yr),]
		 ou$mean = ou$upper = ou$lower = 0
		 boot.mean <- function(x,i){mean(x[i])}
     #get the bs estimates for snowcrab
     
		 for(j in 1:nrow(ou)) {
  	  lp = set[which(set$yr==ou[j,'yr'] & set$cfa==ou[j,'cfa']),'mass']
      if(all(lp==0)) next()
 	    ir = boot(lp,statistic=boot.mean,R=1000)
 	    irr = boot.ci(ir,type='bca')
 	    ou[j,c('lower','upper','mean')] = c(irr$bca[c(4,5)],ir$t0)
      }
      } else {
      ou = NULL
      }

 gr.dir = file.path(project.datadirectory('groundfish'),'analysis')
 gr.list = list()

 for(i in 1:n.regions) {
 		io = paste('stratified',species,'summer.strata',areas.groundfish[[c(i,1)]],areas.groundfish[[c(i,2)]],sep=".")
 		ioi = grep(io,dir(gr.dir))
 		if(length(ioi)==0) stop(paste(paste(gr.dir,io,sep="/"),'not found'))
 		if(length(ioi)>1) stop(paste(paste(gr.dir,io,sep="/"),'has multiple entries'))
 		load(paste(gr.dir,dir(gr.dir)[ioi],sep="/"))
 		out = out[,c('yr','w.yst','w.ci.yst.l','w.ci.yst.u')]
 		names(out) = c('yr','mean','lower','upper')
 		out$cfa = areas.crab[i]
 		out$mean = out$mean / sakm2
 		out$lower = out$lower / sakm2
 		out$upper = out$upper / sakm2
 		gr.list[[i]] = out
 		}
gr.list = do.call(rbind,gr.list) 
gr.list$grp = 'groundfish'
ou$grp = 'snowcrab'

xlim=c(1968,p$current.assessment.year+2)
oi = gr.list
 if(!is.null(ou)) oi =rbind(ou,gr.list)

###TO DO## this figure is not working properly but can run through one at a time.  AMC Feb2015

for(l in areas.crab) {

ooi = oi[oi$cfa==l,]
     fn = file.path( outdir, paste( "groundfish.snowcrab.survey.species.one.axis",species,l, "png", sep="." ) )

    x11()
    browser()
    panel.ci <- function(x, y, ly, uy, subscripts, pch = 20, col=c('black','red'), ...)   { 
        x <- as.numeric(x) 
        y <- as.numeric(y) 
        ly <- as.numeric(ly[subscripts]) 
        uy <- as.numeric(uy[subscripts]) 
        panel.arrows(x, ly, x, uy, col = col, 
                 length = 0.25, unit = "native", 
                 angle=90,code=3)
        panel.xyplot(x, y, pch = pch, col= col,, ...) 
        } 
xyplot(mean ~ yr, 
        groups=grp, 
        data=ooi, 
        ly = ooi$lower, 
        uy = ooi$upper, 
        prepanel = function(x,y,uy)
         list(ylim=c(0,max(uy)*0.75))
        ,ylab=expression("kg per km"^2),
        xlab='Year',
        main=unique(ooi$cfa),
        panel = panel.superpose, 
        panel.groups = panel.ci, 
        type="b",col=c('black','red'),pch=20,lwd=3
       ) 
  # savePlot(fn,type='png')
   print(fn)
 }
 }  

