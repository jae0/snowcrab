  make.timeseries.data = function( areas=NULL, p=NULL, outfile=NULL, reduced.stations=F,vars=NULL ) {
    set = snowcrab.db( DS ="set.complete", p=p )
    #set2015 = set[which(set$yr == 2015),]
    #print (head(set2015))
   
    
   # if (is.null(areas)) areas = c( "cfa20", "cfa21", "cfa22", "cfa23", "cfa24", "cfa4x", "cfa23slope", "cfa24slope", "cfaslope", "cfanorth", "cfasouth", "cfaall" )
   if (is.null(areas)) areas = c( "cfa4x", "cfanorth", "cfasouth")
   #if (is.null(areas)) areas = "cfasouth"

   if(is.null(vars)) variables =  variable.list.expand("all.data")
   if(!is.null(vars)) variables=vars
    print("Variables Created")
    print(variables)
    tsdata =  get.time.series (set, areas, variables, from.file=F, outfile= outfile, reduced.stations=reduced.stations)  # this returns 1.96SE as "se"
    return( paste("Saved: ts.rdata in",  file.path(project.datadirectory("snowcrab"), "R") ) )
  }
