
# to obtain the kriging weights

require (gstat)
data (meuse)
data(meuse.grid)
  
  coordinates(meuse) = ~x+y
  gridded(meuse.grid) = ~x+y
  m <- vgm(.59, "Sph", 874, .04)
  
  # ordinary kriging:
  ok <- krige(log(zinc)~1, meuse, meuse.grid, model = m)
   

kriging.weights = function(x, formula, newdata, model) {
  weighti = function(x, i, formula,...) {
     ret =rep(0,nrow(x))
     ret[i]=1
     x[[1]]=ret
     krige(formula=formula, locations=x, model=model, newdata=newdata )
 }
 ret = sapply(1:nrow(x), weighti, x=x, newdata=newdata[1,], model=model,formula=formula)
 ret = t(sapply(ret, as.data.frame))
 unlist(ret[,3])
}


kw = kriging.weights(meuse["zinc"], zinc~1, as(meuse.grid, "SpatialPoints")[1], model=vgm(1,"Exp",300))

# kriging.weights(x=meuse["zinc"], formula=log(zinc)~1, newdata=SpatialPoints(meuse.grid[,1:2]), model=vgm(1,"Exp",300))


 CC[,"kv"] = 0
 CC[1,"kv"] = 1

 ll = krige(formula=p$kformula, locations=~plat+plon, model=vgm.m, newdata=PS, data=CC)
 ll[3]



