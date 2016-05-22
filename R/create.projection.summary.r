
create.projection.summary = function( p, pFB, rescale.results=F ) {
  
  outdata = NULL

  for (rr in c(1:p$nregions)) {
    td = pFB[ which(pFB$region==p$regions[rr]) , ]

    td0 = td[ which(as.character(td$vars)=="fb.prefishery"),] ; td0$fbpre = td0$val; td0$val=NULL; td0$vars=NULL
    td1 = td[ which(as.character(td$vars)=="error"),] ; td1$error = td1$val; td1$val=NULL; td1$vars=NULL
    td2 = td[ which(as.character(td$vars)=="landings"),] ; td2$landings = td2$val; td2$val=NULL; td2$vars=NULL

    td = merge( td0, td1, by=c("yr","region","ER") )
    td = merge( td,  td2, by=c("yr","region","ER") )

    td$ER = as.factor(td$ER)


    scale.factor = 1000 # convert to X 1000 t

    if ( rescale.results ) scale.factor = td[ which( td$yr==p$start.projection.year ), "fbpre"] [1]

    td$fbpre = td$fbpre / scale.factor
    td$error = td$error / scale.factor
    td$landings = td$landings / scale.factor

    td = td[ order(td$region, td$ER, td$yr) , ]
    td$ub = td$fbpre + td$error  # already 2SD
    td$lb = td$fbpre - td$error
     
    i.bigerror = which( td$error >= td$fbpre) 
    td$lb[ i.bigerror ] = 0
    td$ub[ i.bigerror ] = 0 
    
    td$region = p$regions[rr]
    
    outdata = rbind(outdata, td)
  }

  return(outdata)
} 


