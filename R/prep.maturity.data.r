

  prep.maturity.data = function(x, sex="male", quantiles=c(0,0.25,0.5,0.75,1) ) {

    # sex codes
    male = 0
    female = 1
    sex.unknown = 2

    # maturity codes
    immature = 0
    mature = 1
    mat.unknown = 2


      x=x[is.finite(x$z),]
      x=x[is.finite(x$t),]
      x=x[is.finite(x$cw),]

      x=x[is.finite(x$yr),]
      x$yr=as.factor(x$yr)

      x$reg = NA # initialise the column / field
      x$reg[ emgis::polygon_inside(x,emgis::polygon_internal_code("cfanorth")) ] = "cfanorth"
      x$reg[ emgis::polygon_inside(x,emgis::polygon_internal_code("cfasouth")) ] = "cfasouth"
      x$reg[ emgis::polygon_inside(x,emgis::polygon_internal_code("cfa4x")) ] = "cfa4x"
      x = x[ which( x$reg %in% c("cfanorth", "cfasouth" ) ) ,]
      x$reg=as.factor(x$reg)

      if (sex=="male") {
        x = x[ which( x$sex==1) ,]
        x = x[,c( "cw", "chela", "mat", "cfa", "z", "t", "yr", "reg", "julian")]
        x = x[ which(is.finite(x$chela)),]
      }
      if (sex=="female") {
        x = x[ which( x$sex==2) ,]
        x = x[,c( "cw", "abdomen", "mat", "cfa", "z", "t", "yr", "reg", "julian")]
        x = x[ which(is.finite(x$abdomen)),]
      }

    ### Switch Maturity Codes such that 0=immature   1=mature

      x$mat = as.numeric( as.character( x$mat) )
      x$mat[which(x$mat==2)]= 0
      x$mat = as.factor(x$mat )
      x = x[ which( x$mat %in% c(0,1)) ,]

      # temp fix as the cut function seems to drop the lowest data points.
      t.q = quantile(x$t, quantiles)
      t.q[1] = t.q[1] - 0.1

      z.q = quantile(x$z, quantiles)
      z.q[1] = z.q[1] - 10

      # discrete codings for temp and depth
      x$tq=cut(x$t, t.q, include.lowest=T)
      x$zq=cut(x$z, z.q, include.lowest=T)

    return ( x )

  }

