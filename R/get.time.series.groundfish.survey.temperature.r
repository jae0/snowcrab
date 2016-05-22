
  get.time.series.groundfish.survey.temperature = function( from.file=F, outfile=file.path(p$annual.results, "timeseries", "survey","groundfish.t.rdata") ) {
    if (is.null(outfile)) outfile = file.path( project.datadirectory("snowcrab"), "R", "ts.rdata" )

    
    if( from.file) {
      if (file.exists(outfile) ) load(outfile)
      return(ts)
    }
   
     out = data.frame(r=NA,yrs=NA,meanval=NA,se=NA,ub=NA,lb=NA,n=NA)
	h = groundfish.db('gshyd')
	g = groundfish.db('gsinf')
	g = g[,c('id','sdate','lon','lat','bottom_temperature')]
	h = h[,c('id','temp')]
	names(h)[2] <- 'bottom_temperature'
	f <-merge(g,h,by='id',all.x=T) 	
	 i <- which(is.na(f$bottom_temperature.x) & !is.na(f$bottom_temperature.y))  
	f[i,'bottom_temperature.x'] <- f[i,'bottom_temperature.y']
	f$yr <- as.numeric(format(f$sdate,'%Y'))
	f <- fishing.area.designations(f)
	ar <- unique(f$cfa)
	yy <- unique(f$yr)
	yy <- yy[order(yy)]
      for (r in ar) {
              for (yrs in yy) {
          y = f[which(f$yr == yrs & f$cfa ==r & !is.na(f$bottom_temperature.x)),]
if(nrow(y)>3) {
	ym <- min(y$bottom_temperature.x[y$bottom_temperature.x>0])
             q = log(y$bottom_temperature.x+ym)

              m =  mean (q, na.rm=T)
              n = length(q)
              se = sd(q, na.rm=T)/ sqrt(n-1)
              meanval = exp(m)-ym
              ub = exp(m+se*1.96)-ym
              lb = exp(m-se*1.96)-ym
              j = as.data.frame(cbind(r, yrs, meanval, se, ub, lb, n))
            out <- rbind(out,j)
		}
          }
}
          colnames(out) = c("region", "year", "mean", "se", "ub", "lb", "n")
        numbers = c("year", "mean", "se", "n", "ub", "lb")
        out = factor2number(out, numbers)
       out <- out[!is.na(out$year),]

        save(out, file=outfile, compress=T)
print(outfile)
return(out)
      }
  

 

