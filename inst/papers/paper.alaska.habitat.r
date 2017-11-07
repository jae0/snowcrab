
	p = bio.snowcrab::load.environment()

  setwd( project.datadirectory("bio.snowcrab") )

  det = snowcrab.db( DS ="det.georeferenced" )
  set = snowcrab.db( DS ="set.complete" )

  xdet = merge( det, set[,c("trip","set", "t", "z")], by=c("trip", "set"), all.x=T, all.y=F )
  xdet = xdet[ which( xdet$sex %in% c(female,male) ), ]

  xdet$sex = as.factor( xdet$sex)
  xdet$size = "small"
  xdet$size[which( xdet$cw > 95) ] =  "big"
  xdet$size[which( xdet$cw > 60 & xdet$cw < 95) ] =  "med"

  xdet$size = as.factor(  xdet$size )
  levels(xdet$mat) = c("Immature", "Mature")
  levels(xdet$sex) = c("Male", "Female")


  setup.lattice.options()

  histogram( ~ t | sex + mat , data = xdet,
    xlab = "Bottom temperature (C)", type = "density",
    panel = function(x, ...) {
      panel.histogram(x, col="white", ...)
      panel.mathdensity( dmath=dnorm, col = "black", args=list(mean=mean(x),sd=sd(x)))
   } )

   Pr( "png", dname=project.datadirectory("bio.snowcrab", "output"), fname="temp.hist", trim=F , res=144)


  histogram( ~ z | sex + mat , data = xdet,
    xlab = "Bottom temperature (C)", type = "density",
    panel = function(x, ...) {
      panel.histogram(x, col="white", ...)
      panel.mathdensity( dmath=dnorm, col = "black", args=list(mean=mean(x),sd=sd(x)))
   } )

   Pr( "png", dname=project.datadirectory("bio.snowcrab", "output"), fname="temp.hist", trim=F , res=144)



  histogram( ~ t | sex + size , data = xdet,
    xlab = "depth", type = "density",
    panel = function(x, ...) {
      panel.histogram(x, col="white", ...)
      panel.mathdensity( dmath=dnorm, col = "black", args=list(mean=mean(x),sd=sd(x)))
   } )
`
   Pr( "png", dname=project.datadirectory("bio.snowcrab", "output"), fname="temp.hist", trim=F , res=144)




