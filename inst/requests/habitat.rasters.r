p = list()
p$init.files = loadfunctions( "snowcrab", functionname="initialise.local.environment.r") 
  p$libs = RLibrary( "mgcv", "chron", "lattice", "lattice", "grid", "fields", "parallel"  ) 


loadfunctions(c('snowcrab','spacetime'))
   basedir = file.path( project.directory("snowcrab"), "R", "gam" )

      loc.map = file.path( basedir, "maps" )
      loc.sol = file.path( basedir, "predictions" )
      loc.res = file.path( basedir, "results" )
  


makeRaster  <- function(yr) {
require(raster)
require(gstat)
  fn.res = file.path( loc.sol, paste( "PS.simulation.means", 'R0.mass', yr, "rdata", sep="." ) )
	load(fn.res)
	d = PS[,c('plon','plat','habitat.mean')]
	d = planar2lonlat(x=d,proj.type=p$internal.projection )
	names(d) = c('plon','plat','z','x','y')
	e = extent(d[,4:5])
	r = raster(e,resolution=0.005,crs="+proj=longlat +datum=WGS84")

mod = gstat(id = 'means', formula = z~1,locations= ~x+y,data=d, nmax = 10, set=list(idp= 0.5),maxdist=0.01)
x = interpolate(r, mod)

writeRaster(x, filename=paste("/home/adam/ecomod/snowcrab/R/Requests/CommercialMaleHabitat",yr,".tif",sep=""), format="GTiff", overwrite=TRUE)
return('Done')
}
