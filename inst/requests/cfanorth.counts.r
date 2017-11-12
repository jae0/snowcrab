

p = bio.snowcrab::load.environment()

  w = snowcrab.db( DS ="set.biologicals", p=p )
  w = w[ emaf::polygon_inside( w, "cfanorth"),]

  out = data.frame( cbind(
    no.male = tapply( w$totno.male, w$yr, sum, na.rm=T) ,
    no.female = tapply( w$totno.female, w$yr, sum, na.rm=T)
  ) )

  out$year = rownames( out )

  write.csv( out, file="cfanorth.counts.csv")



