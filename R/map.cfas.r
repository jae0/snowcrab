map.cfas = function( p, conversions=c("ps2png") ) {
  #Import coastline
  #MG: Switch this to the smaller coastline with no islands
  library('rgdal')
  internal.crs <- "+proj=utm +zone=20 ellps=WGS84"
  geog.proj <- CRS("+proj=longlat +ellps=WGS84")
  #seis <- colorRampPalette(c("darkblue","blue3", "green", "yellow", "orange","red3", "darkred"), space = "Lab")
  #loadfunctions("bathymetry")
  
  polydir = file.path(project.datadirectory("polygons"), "data", "Basemaps", "Marine", "Coastline")
  
  setwd(polydir)
  
  coast<-readOGR(".", "NY_to_Nova_UTM20")
  coast<-spTransform(coast, geog.proj)
  coast <- gSimplify(coast, tol=0.01, topologyPreserve=TRUE)
  
  
  cfadir =  file.path(project.datadirectory("polygons"), "data", "Management_Areas", "Fisheries", "Snowcrab")
  setwd(cfadir)
  
  cfa20 = read.table("cfa20.dat")
  cfa21 = read.table("cfa21.dat")
  cfa22 = read.table("cfa22.dat")
  cfa23 = read.table("cfa23.dat")
  cfa24 = read.table("cfa24.dat")
  cfa4x = read.table("cfa4x.dat")
  cfaall = read.table('cfaall.dat')
  cfasouth = read.table("cfasouth.dat")
  cfanorth = read.table("cfanorth.dat")
  
  cfas = list(cfa20, cfa21, cfa22, cfa23, cfa24, cfa4x, cfaall, cfasouth, cfanorth)
  names(cfas) = c("cfa20", "cfa21", "cfa22", "cfa23", "cfa24", "cfa4x", "cfaall", "cfasouth", "cfanorth")
  
  
  for (i in 1:length(cfas)){
    name = names[i]
    print(name)
    print(cfas[i])
    p = Polygon(cfas[i])
    p2 = Polygons(list(p),1)
    sp = SpatialPolygons(list(p2))
    plot(sp)
    plot(coast, col='lightgrey', add=TRUE)
    df = data.frame(id=getSpPPolygonsIDSlots(sp))
    sdf = SpatialPolygonsDataFrame(sp, data=df)
    writeOGR(sdf, ".", name,"ESRI Shapefile", overwrite=TRUE)
  }
  
  cfa23.p = Polygon(cfa23)
  cfa23.p2 = Polygons(list(cfa23.p),1)
  cfa23.sp = SpatialPolygons(list(cfa23.p2))
  df<- data.frame(id=getSpPPolygonsIDSlots(cfa23.sp))
  cfa23.df = SpatialPolygonsDataFrame(cfa23.sp, data=df)
  writeOGR(cfa23.df, ".", "CFA23", "ESRI Shapefile")

  plot(spcfa23)
  plot(coast, col='lightgrey', add=TRUE)
  
}
  