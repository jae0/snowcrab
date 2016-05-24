 
    	
		loadfunctions( "snowcrab", functionname="initialise.local.environment.r") 
  

  # maturity temporal transition

  det = snowcrab.db("det.georef")
  to.keep = c("trip", "set", "sex", "cw", "mass", "abdomen", "chela", "mat",
    "shell", "gonad", "eggcol", "eggPr", "durometer", "legs",
    "station", "t0", "lon", "lat",
    "chron", "julian", "yr", "z", "t",
    "zsd", "tsd", "sa"
  )
  x = det[,to.keep]
  rm (det)
  gc()

  # female:
  female = 2 
  x = x[ x$sex == female , ]
  x$ mat = as.numeric(as.character(x$mat))
  x$jul = cut(x$julian, breaks=seq(0,370, 10) )
  x$cw2 = cut(x$cw, breaks=seq(0,100, 10) )


#   if (type=="f.berried")  i = which(x$sex==2 & x$eggPr %in% c(1:4) )
#    if (type=="primiparous") i = which(x$sex==2 & x$mat==1 & x$shell<=2 )
#    if (type=="multiparous") i = which(x$sex==2 & x$mat==1 & x$shell>=3 )
#    if (type=="senile") i = which(x$sex==2 & x$mat==1 & x$shell >=4 & x$eggPr <=1 )
#    if (type=="adolescent") i = which(x$sex==2 & x$mat==2 & x$gonad<=2)
#    if (type=="preprimiparous") i = which(x$sex==2 & x$mat==2 & x$gonad==3)

  pre.primi = which(x$sex==2 & x$mat==2 & x$gonad>=3)
  primi = which( x$sex==2 & x$mat==1 & x$shell<=3 )  
  population =  sort( unique( c( primi, pre.primi)  ))

  # recode maturity to: 0 (imm) 1 (mat)
  tmp = x$mat
  x$mat = NA
  x$mat[ which(tmp==1)] = 0
  x$mat[ which(tmp==2)] = 1


  x = x[ which(x$mat %in% c(0,1)) ,]
  x = x[ which(x$cw >= 10) ,]
  

  mm = glm( mat ~ jul + cw2 , data=x, weight=sa, family=binomial(link = "logit") ) 
  summary(mm)

  car::Anova(mm)
  mme <- effects::all.effects(mm, se=F )
  mmep = effects:plot(mme,  type=c("response"))
  epicalc::logistic.display(mm)
 
