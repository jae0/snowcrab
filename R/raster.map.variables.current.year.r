
#Function to summarize and map raster variables
#Cell size can also be adjusted to change the default cell size that is mapped, default = 6km

	raster.map.variables.current.year = function(p=p, grid.fun, variables, cell=NULL, years=years) {
	  
	  browser()
	  
	  K = logbook.db( DS="logbook" )
	  #x$landings = x$landings/1000  # convert kg to ton #MG code this out so that it calculates quantiles properly
	  K$sa = 1  # this a dummy variable required by the mapping routine
	  K = K [filter.region.polygon( K, region="isobath1000m"),]
	  K = K[ which(K$effort <= 300) ,]
	  K = K[ which(K$cpue < 500),]
	  
	  polydir = file.path(project.datadirectory("polygons"), "data", "Basemaps", "Marine", "Coastline")
	  rasdir = file.path(project.datadirectory("snowcrab"), "maps", "rasters", "logbook")
	  mapdir = file.path(project.datadirectory("snowcrab"), "maps", "images", "logbook", "currentyear")
	  
	  plottimes = p$plottimes
		extent=p$extUTM
		p.ext2 = extent(matrix(c(-66.4, 42.2, -57.2, 47.4), nrow=2, ncol=2))
		cell=NULL
		internal.crs <- "+proj=utm +zone=20 ellps=WGS84"
		geog.proj <- CRS("+proj=longlat +ellps=WGS84")
		seis <- colorRampPalette(c("darkblue","blue3", "green", "yellow", "orange","red3", "darkred"), space = "Lab")
		loadfunctions("bathymetry")
		# in meters
		ib = isobath.db( depths=c(100, 200) 
		                 #, return.lonlat=TRUE 
		                 )
		
		#Import coastline
		#MG: Switch this to the smaller coastline with no islands
		setwd(polydir)
		coast<-readOGR(".", "NY_to_Nova_UTM20")
		coast<-spTransform(coast, geog.proj)
		coast <- gSimplify(coast, tol=0.01, topologyPreserve=TRUE)
		
		#Designate Cell size
		if(is.null(cell)) {
			# 1 minute grid = 0.166, 2 minute grid = 0.033, 3minute grid = 0.05, 4 minute grid = 0.0666, 5 minute grid = 0.083
			cell<- 0.083
			cell.big <- 0.083
		}

		#Set rows and columns for blank grid
		ncols <- length(p.ext2[2]:p.ext2[1])/cell
		nrows <- (length(p.ext2[4]:p.ext2[3])-1)/cell
		
		#Set rows and columns for the raster stack. Larger cell size because of smaller plots
		ncols.big <- length(p.ext2[2]:p.ext2[1])/cell.big
		nrows.big <- (length(p.ext2[4]:p.ext2[3])-1)/cell.big

		#Create the color palette for each variable
		for (v in variables) {
			print(v)
			#Extract data for the raster creation
			M = K[, c("yr", "lon", "lat", v)]
			M = M[is.finite(M[,v] *M[,"lon"] *M[,"lat"]),]

			#Create a blank grid for plotting
			grid <- raster(nrows=nrows, ncols=ncols, ext=p.ext2, crs=geog.proj)
			grid[]<- 1:ncell(grid)
			grid.big <- raster(nrows=nrows.big, ncols=ncols.big, ext=p.ext2, crs=geog.proj)
			grid.big[]<- 1:ncell(grid.big)
			rstack <- stack()
	
			#Loop through and create a raster-based plot for each year
			for (i in 1:length(sort(unique(years)))) {
			  y= years
			  #y= sort(unique(M$yr))
				name<- y[i]
				print(name)
				s= subset(M, M$yr==name)

				#Extract coordinates and create spatial points data frame to be gridded
				cor = s[, c("lon", "lat")]
				sPDF <- SpatialPointsDataFrame (cor, data=s)
				proj4string(sPDF)<-(geog.proj)

				ras.name <- paste(v, name, ".tif", sep="")
				grid.sum <- rasterize(x=sPDF, y=grid, field=variables, fun=grid.fun)
				grid.sum.big <- rasterize(x=sPDF, y=grid.big, field=variables, fun=grid.fun)
				setwd(rasdir)
				writeRaster(grid.sum, filename=ras.name, datatype= 'GTiff', overwrite=T)
				
				#Create raster stack to calculate the max and min values for the colour scale
				rstack <- stack(rstack, grid.sum.big)
				names(rstack[[i]]) <- name
				
				#Define colour scale
				z <- getValues(rstack)
				z <- z[is.finite(z)]
				z <- round(z, digits=0)
				
				if (v == "landings") {
					z <- z/100
					z <- round(z, digits=0)
					z <- z*100
				} else {
					z <- z/10
					z <- round(z, digits=0)
					z <- z*10
				}

				#MG might want to convert landings /1000 to calculate tons for landings, right now it's using KG
				quant <- unique(quantile(z, seq(0,1, length.out=75)))
			  #MG check to see if quantiles are being calculated properly for effort, there doesn't seem to be a lot of red on the maps
				quant.small <- unique(quantile(z, seq(0,1, length.out=5)))
				ckey <- list(at=quant, labels=list(at=quant.small))
				#Plot the variable	
				setwd(mapdir)
			    fig.name <- paste(v, name, ".pdf", sep="")
				#png(filename=fig.name, width=6, height=5, units="in", res=450)
				pdf(file=fig.name, width=6.5, height=4.8, bg='white')
				#print(levelplot(grid.sum, at=quant, colorkey=ckey, col.regions=seis, 
				#margin=F, xlab="", ylab="", main=name, scales = list(x=list(cex=0.7), y=list(cex=0.7))) + layer(sp.polygons(coast, fill='lightgrey'))+ layer(sp.polygons(coast, fill='lightgrey')) + layer(sp.lines(ib, col='gray10', alpha= 0.6, lwd= 0.6)))
				
				print(levelplot(grid.sum, at=quant, colorkey=ckey, col.regions=seis, alpha.regions=1,
				                margin=F, xlab="", ylab="", main=name, scales = list(x=list(cex=0.7), y=list(cex=0.7))) 
				      + layer(sp.polygons(coast, fill='lightgrey')) 
				      + layer(sp.polygons(coast, fill='lightgrey')) 
				      + layer(sp.lines(ib, col='gray10', 
				          alpha= 0.3, 
				          lwd= 0.8)))
				dev.off()
		
		
			}
			}
		return("mapping completed")
		}

	