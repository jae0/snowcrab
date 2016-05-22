
  get.time.series = function(x=NULL, regions=NULL, vars=NULL, trim=0, from.file=F, outfile=NULL,reduced.stations=F) {

 #browser()
    if (is.null(outfile) & !reduced.stations) outfile = file.path( project.datadirectory("snowcrab"), "R", "ts.rdata" )
    if (is.null(outfile) & reduced.stations) outfile = file.path( project.datadirectory("snowcrab"), "R", "ts.reduced.rdata" )
    
   ts= NULL

    if( from.file) {
      if (file.exists(outfile) ) load(outfile)
      return(ts)
    }
    # ra.jit(2)  # JIT optimisation
 if(reduced.stations) {
    b = x[which(x$yr==2014),'station']
    x = x[which(x$station %in% b),]
    }
     # x2015 = x[which(x$yr == 2015),] # check to see if x has data in it
    #  print (head(x2015))
      
      out = NULL
      tmpfile = file.path( tempdir(), make.random.string(".ts.csv"))

      for (r in regions) {
        print (r)
        qr = filter.region.polygon(x, r)
        #print(head(qr))
        for (yrs in unique(x$yr)) {
          print(yrs)
          qy = which( x$yr == yrs )
         # print(head(qy))
          filter = sort(unique(intersect(qr, qy)))
          
          if (length(filter)>0) {
            for (v in vars) {
              print(v)
              z = length( which( is.finite( x[ filter, v] ) ) )
              #print(head(z))
              if (z < 1) next()
              y = x[ filter, ]
              #print(head(y))
              #xi = x[which(x$trip == 'S01122015' & x$set==1),]
              #print(xi$R0.mass)
              #if(v=='t' & yrs==2003) 
              a = variable.recode( y[,v], v, direction="forward", db="snowcrab" ) # transform variables where necessary
              m =  mean (a, trim=trim, na.rm=T)
              n = length(a)
              se = sd(a, na.rm=T)/ sqrt(n-1)
              meanval = variable.recode (m , v,  direction="backward", db="snowcrab" )
              ub = variable.recode ( m+se*1.96, v,  direction="backward", db="snowcrab" )
              lb = variable.recode ( m-se*1.96, v,  direction="backward", db="snowcrab" )
              j = as.data.frame(cbind(r, yrs, v, meanval, se, ub, lb, n))
              #print(j)
              write.table(j, file=tmpfile, append=T, row.names=F, col.names=F, quote=F, sep=";")
            }
          }
      }}

      ts = NULL
      if (file.exists( tmpfile)) {
        ts = read.table(tmpfile, sep=";", as.is=T, colClasses="character", header=F)
        colnames(ts) = c("region", "year", "variable", "mean", "se", "ub", "lb", "n")
        numbers = c("year", "mean", "se", "n", "ub", "lb")
        ts = factor2number(ts, numbers)
        #print ('ts whole')
        #print(ts)
        #ts2015 = ts[which(ts$year == 2015),]
        #print('2015')
        #print(head(ts2015))
        save(ts, file=outfile, compress=T)
      }
      remove.files(tmpfile)
     return(ts)

  }


