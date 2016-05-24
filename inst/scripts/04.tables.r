

  loadfunctions( "snowcrab", functionname="initialise.local.environment.r") 



# Tables for the SSR/RESDOC
# ---------------------------------------------
  # Tables of CC from trawl surveys > 95mmCW

  outvars = c("yr", "region", "vars", "total", "lbound", "ubound")
  yy = c(1999:p$current.assessment.year)
  rr = c("cfasouth", "cfanorth", "cfa4x")
  vv = c( "totno.male.com.CC1", "totno.male.com.CC2", "totno.male.com.CC3",
          "totno.male.com.CC4", "totno.male.com.CC5" )

  p = make.list( list(y=yy, v=vv ), Y=p )
  L = interpolation.db( DS="interpolation.simulation", p=p )

  L$vars = as.character(L$vars)

  for (y in yy) {
  for (r in rr) {
    i = which( L$yr==y & L$region==r )
#    isum = sum(L$total[i], na.rm=T)
    L[i, c("total", "lbound", "ubound")] =round(  L[i, c("total", "lbound", "ubound")], 4)
  }}

  scale.factor = 1e2
  L[, c("total", "lbound", "ubound")] = L[, c("total", "lbound", "ubound")] * scale.factor
  L$vars = factor(x=L$vars, levels=vv, ordered=T)

  Mn = xtabs(formula=total~yr+vars, data=L[which(L$region=="cfanorth"),], drop.unused.levels = T ) / scale.factor
  dimnames(Mn)$vars = c(1:5)

  Ms = xtabs(formula=total~yr+vars, data=L[which(L$region=="cfasouth"),], drop.unused.levels = T ) / scale.factor
  dimnames(Ms)$vars = c(1:5)

  Mx = xtabs(formula=total~yr+vars, data=L[which(L$region=="cfa4x"),], drop.unused.levels = T ) / scale.factor
  dimnames(Mx)$vars = c(1:5)

 (Mnp = round(Mn/rowSums(Mn)*100,2))
 (Msp = round(Ms/rowSums(Ms)*100,2))
 (Mxp = round(Mx/rowSums(Mx)*100,2))

  latex( Mnp, file="", title="", label="table.CC.north.trawl", rowlabel="Year", cgroup="Carapace condition", na.blank=T, caption="Carapace condition of crab larger than 95 mm CW (percent by number) over time for N-ENS from trawl surveys. The transition from a spring to a fall survey occurred in 2002/2003.")

  latex(Msp, file="", title="", label="table.CC.south.trawl", rowlabel="Year", cgroup="Carapace condition", na.blank=T, caption="Carapace condition of crab larger than 95 mm CW (percent by number) over time for S-ENS from trawl surveys. The transition from a spring to a fall survey occurred in 2002/2003.")

  latex(Mxp, file="", title="", label="table.CC.4x.trawl", rowlabel="Year", cgroup="Carapace condition", na.blank=T, caption="Carapace condition (percent by number) over time for CFA 4X from trawl surveys.")




# ---------------------------------------------
  # Tables of fishable biomass from trawl surveys

 
  outvars = c("yr", "region", "vars", "total", "lbound", "ubound")
  yy = c(1999:p$current.assessment.year)
  rr = c("cfasouth", "cfanorth", "cfa4x")
  vv = c( "totmass.male.com", "R0.mass", "R0.no", "R1.no", "R2.no", "R3.no", "totno.male.imm", "totno.male.com.CC1to2", "totno.male.com.CC3to4",
          "totno.male.com.CC5" )
 
  p = make.list( list(y=yy, v=vv ), Y=p )
  L = interpolation.db( DS="interpolation.simulation", p=p )
  # L = K[K$yr %in% yy & K$vars %in% vv & K$region %in% rr, outvars]
  L$vars = as.character(L$vars)

  for (y in yy) {
  for (r in rr) {
    i = which( L$yr==y & L$region==r )
#    isum = sum(L$total[i], na.rm=T)
    L[i, c("total", "lbound", "ubound")] =round(  L[i, c("total", "lbound", "ubound")], 4)
  }}


   scale.factor = 1
  L[, c("total", "lbound", "ubound")] = L[, c("total", "lbound", "ubound")] * scale.factor
  L$vars = factor(x=L$vars, levels=vv, ordered=T)

  Mn = xtabs(formula=total~yr+vars, data=L[which(L$region=="cfanorth"),], drop.unused.levels = T ) / scale.factor
  dimnames(Mn)$vars = vv

  Ms = xtabs(formula=total~yr+vars, data=L[which(L$region=="cfasouth"),], drop.unused.levels = T ) / scale.factor
  dimnames(Ms)$vars = vv

  Mx = xtabs(formula=total~yr+vars, data=L[which(L$region=="cfa4x"),], drop.unused.levels = T ) / scale.factor
  dimnames(Mx)$vars = vv

  options(width=240)
  options(digits=1)
  Mn
  Ms
  Mx



  latex(Mn, file="", title="", label="table.stats.north.trawl", rowlabel="Year", cgroup="Various statistics", na.blank=T, caption="Various statistics of crab for N-ENS from trawl surveys. The transition from a spring to a fall survey occurred in 2002/2003.")

  latex(Ms, file="", title="", label="table.stats.south.trawl", rowlabel="Year", cgroup="Various statistics", na.blank=T, caption="Various statistics for S-ENS from trawl surveys. The transition from a spring to a fall survey occurred in 2002/2003.")

  latex(Mx, file="", title="", label="table.stats.4x.trawl", rowlabel="Year", cgroup="Various statistics", na.blank=T, caption="Various statistics for CFA 4X from trawl surveys.")


