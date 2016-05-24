
require(geoR)
require(geoRglm)



      vg = variog( coords=Q[,c("plon", "plat")], data=Q$t, max.dist=100 )
        vgm = variofit(vg); rm(vg) ;gc()

        # requires 32 GB
        pred.grid = expand.grid( plon=p$plons, plat=p$plats)
        inside = filter.region.polygon( pred.grid[,c(1:2)], region="cfaall", planar=T,  proj.type=p$internal.projection )
        pred.grid = pred.grid [ inside ,]
        rm(inside); gc()
       

        ee = krige.conv( coords=Q[,c("plon", "plat")], data=Q$t, locations=pred.grid, krige = krige.control( obj.model = vgm) )





data(elevation)

x11();points(elevation, cex.min=1, cex.max=4, pt.div="quint")
x11();plot(elevation, lowess=T)

str(elevation)
names(elevation)


data(ca20)
ca20.3 = subset(ca20, area=3)
ca20.x = subarea(ca20, xlim=c(5200,5600))

str(ca20)
names(ca20)







# ----------------------
  
	loadfunctions( "snowcrab", functionname="initialise.local.environment.r") 
 
kvar = "R0.mass"
kyear = 2007

formatted.data = get.PS.S.gridded (p, kyear, kvar )


CC = formatted.data$CC
cnames = names( CC )
ydata = which( cnames=="kv" )
coords = c( which( cnames=="plon"), which( cnames=="plat") )
# covars = c( which(cnames=="z"), which(cnames=="t"), which(cnames=="dZ"), which(cnames=="ddZ") )
covars = c( which(cnames=="z"), which(cnames=="t") )
CC = CC[ is.finite(rowSums( CC[, c(coords,ydata) ])) ,]  # do not rely upon na removal within routines
CC = as.geodata( CC, coords.col=coords, data.col=ydata, covar.col=covars, covar.names=c("z","t"), na.action="ifany", rep.data.action=mean )
rm(coords, covars)


PS = formatted.data$PS
totalsurfacearea = dim(PS)[1]
psnames = names( PS )
pscoords = c( which(psnames=="plon"), which( psnames=="plat") )
pscovars = c( which(psnames=="z"), which(psnames=="t")) 
PS = as.geodata( PS, coords.col=pscoords, covar.col=pscovars, covar.names=c("z","t"), na.action="ifany", rep.data.action=mean )
rm(pscoords, pscovars)
rm(formatted.data); gc()


uvec = c( 10, seq( 20, 160, by=20 ) )

# cte= mean constant
vgm0 = variog( CC, trend="cte", uvec=uvec, option="bin", estimator="classical", max.dist=150, pairs.min=3 )
plot(vgm0, type="b", lty=1, ylim=c(0,1))

# 1st = mean is a first order polynomial of the coords
# mean ~ x1 + x2
vgm1 = variog( CC, trend="1st", uvec=uvec, option="bin", estimator="classical", max.dist=150, pairs.min=3 )
lines(vgm1, lty=1, col="orange")

# 2nd = mean is a second order polynomial of the coords
# mean ~ x1 + x2 + I(x1^2) + I(x2^2) + I(x1*x2) 
vgm2 = variog( CC, trend="2nd", uvec=uvec, option="bin", estimator="classical", max.dist=150, pairs.min=3 )
lines(vgm2, lty=2, col="blue")


# with covariates
vgm.z = variog( CC, trend=~CC$coords+CC$covar$z, uvec=uvec, option="bin", estimator="classical", max.dist=150, pairs.min=3 )
lines(vgm.z, lty=3, col="red")

vgm.t = variog( CC, trend=~CC$coords+CC$covar$t, uvec=uvec, option="bin", estimator="classical", max.dist=150, pairs.min=3 )
lines(vgm.t, lty=3, col="brown")

vgm.dZ = variog( CC, trend=~CC$coords+CC$covar$dZ, uvec=uvec, option="bin", estimator="classical", max.dist=150, pairs.min=3 )
lines(vgm.dZ, lty=3, col="pink", pty=3)

vgm.ddZ = variog( CC, trend=~CC$coords+CC$covar$ddZ, uvec=uvec, option="bin", estimator="classical", max.dist=150, pairs.min=3 )
lines(vgm.ddZ, lty=3, col="gray", pty=4)

vgm5 = variog( geodata=CC, trend=~CC$coords+CC$covar$t+CC$covar$z, uvec=uvec, option="bin", estimator="classical", max.dist=150, pairs.min=3 )
lines( vgm5, lty=3, col="blue")

vgm6 = variog( CC, trend=~CC$coords+CC$covar$t+CC$covar$z+CC$covar$dZ, uvec=uvec, option="bin", estimator="classical", max.dist=150, pairs.min=3 )
lines( vgm6, lty=3, col="blue")

vgm7 = variog( CC, trend=~CC$coords+CC$covar$t+CC$covar$z+CC$covar$dZ+CC$covar$ddZ, uvec=uvec, option="bin", estimator="classical", max.dist=150, pairs.min=3 )
lines( vgm7, lty=6, col="green", pty=6)


# modelled variograms
vm0 = variofit( vgm5, cov.model = "matern" )  # wls regression estaimtes

vm1 = likfit(  CC, trend=~CC$coords+CC$covar$t+CC$covar$z, ini.cov.pars=vm0$cov.pars, lik.method = "ML", cov.model = "matern", lambda=1 )

vm1.prof = proflik(vm1, CC, uni.only = T)


ocontrol = output.control(n.posterior=10, moments=T, mean.var=T, quantile=c(0.025, 0.5, 0.975) )

kcontrol = krige.control(type.krige = "sk", 
  trend.d =~CC$coords+CC$covar$t+CC$covar$z, 
  trend.l = ~PS$coords+PS$covar$t+PS$covar$z,
  obj.model = vm0,
  lambda=1
)

gc()
k = krige.conv( geodata=CC, locations=PS$coords, borders=NULL, krige=kcontrol, output=ocontrol )

k = ksline( geodata=CC, locations=PS$coords, borders=NULL, cov.pars=vm0$cov.pars, krige=kcontrol, output=ocontrol, m0="kte", nwin=30, ktedata=PS$covar, ktelocations=PS$coords )



set.seed(1) 


krige.glm.control( type.krige = "sk", 
  trend.d =~CC$coords+CC$covar$t, 
  trend.l= ~PS$coords+PS$t+..., 
  obj.model = NULL, beta, cov.model, cov.pars, kappa,
  nugget, micro.scale, dist.epsilon = 1e-10, 
  aniso.pars, lambda 
)


krweights

