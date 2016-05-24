
# species composition analysis from snow crab data
  	
	loadfunctions( "sorted.ordination") 
	loadfunctions( "snowcrab", functionname="initialise.local.environment.r") 


  C = snowcrab.db( DS="cat.georeferenced" )
  C$totno = log10( C$totno )
  C$totmass = log10( C$totmass )

  sp = sort( unique( C$spec ) )

  C$id = paste(C$trip, C$set, sep=".")
 
  m = matrix.form( C, type="number" )

  gstaxa = taxonomy.db( "gstaxa" )
  
  ss = data.frame( spec= as.numeric(colnames( m))  )
  ss = merge(ss, gstaxa, by="spec", all.x=T, all.y=F, sort=F )
  colnames( m) = ss$namecom
  ll =  rownames(m)


    corel = cor( m, use="pairwise.complete.obs" ) # set up a correlation matrix ignoring NAs
    corel[ is.na(corel) ] = 0
    
    s = svd(corel)  # eigenanalysis via singular value decomposition
    
    scores = matrix.multiply (m, s$v)  # i.e., b %*% s$v  .. force a multiplication ignoring NAs
    evec = s$v
    eval = s$d
  
    x = cbind( scores[,1] / sqrt(eval[1] ), scores[,2] / sqrt( eval[2]) )
    y = cbind( evec[,1] * sqrt(eval[1] ) , evec[,2] * sqrt( eval[2]) )

  
    fn = file.path( tempdir(), "temp.pdf" )
    Cairo( file=fn, type="pdf", dpi=150 )
      plot(y, pch=".", col="red" )
      text(y,labels= colnames(m), cex=0.6 )
    dev.off()
    sys( "inkscape &", fn )

    plot(x, pch=".", col="blue" )
    text(x,labels= rownames(m), cex=0.6 )
 

  sm = merge(sm, scores, by="id", all.x=T, all.y=F, sort=F)
  save (sm, file=fn, compress=T)


  
  ######## density modelling GAM/GLM
  
  loadfunctions( "snowcrab", functionname="initialise.local.environment.r" )

  S = snowcrab.db( DS ="set.logbook" ) 
  Sn = colnames( S )

  iv0 = Sn[ grep( "^ms.no." , Sn ) ]
  isp = as.numeric( as.character( gsub( "ms.no.", "", iv0 ) ))
  tx = data.frame( spec=isp )
  tx = merge(tx, taxonomy.db("gstaxa"), by="spec", all.x=T, all.y=F, sort=F )
  #  tx = tx [ which(tx$spec<3000) ,]  # an initial lm suggests that those with spec# >3000  have minimal linear influence .. drop as there are too many vars
  txsp= taxonomy.db( from="spec", to="parsimonious", tolookup=2526 )
  tx = tx [ - which(tx$spec==txsp), ] #get rid of snow crab
    
  iv = gsub("[[:space:]]+", "", tx$namecom )
  iv = gsub("[[:punct:]]+", "", iv )
  iv = paste(iv, tx$spec, sep="_"  )
  iv = make.names( iv ) 

  # rename variables 
  for ( i in 1:length(iv) ) {
    names(S)[which(names(S)==iv0[i])] = iv[i]
  }
  
  corevars = c( "totno.all", "yr", "cfa", "plon", "plat", "t", "tamp", "wmin", "z", "substrate.mean", "dZ", "ddZ" )
  tolog =  c("R0.mass", "totno.all", iv )
  for( i in tolog) S[,i] = log(S[,i] + 1 )
  allvars = c( iv, corevars )
  S = S[,allvars]
 
  iv = c("AtlanticCod_10", "Haddock_11", "WhiteHake_12", "RedHake_13", "SilverHake_14", "Cusk_15","Pollock_16",
    "AtlanticHalibut_30", "AmericanPlaice_40","WitchFlounder_41", "YellowtailFlounder_42", "WinterFlounder_43",
    "StripedAtlanticWolffish_50", "SpottedWolffish_51", "NorthernWolffish_52","Capelin_64", "BarndoorSkate_200",
    "ThornySkate_201","SmoothSkate_202", "WinterSkate_204", "SpinyDogfish_220","NorthernHagfish_241", "LonghornSculpin_300", 
    "AmericanSandLance_599", "PandalusBorealis_2211", "ToadCrab_2520","AtlanticRockCrab_2513", "NorthernStone_2523")
  depvars = c(  "t", "tamp", "wmin", "z", "substrate.mean", "dZ", "ddZ",  iv )
  
  # depvars = depvars[ - which(depvars=="SnowCrabQueen_2526") ]
  
  yvar = "totno.all"
 
  # linear model
   
  fl = "totno.all~1 +  t  +  tamp+  wmin+  z  +  substrate.mean  +  dZ  +  ddZ  +  AtlanticCod_10  +  Haddock_11  +  WhiteHake_12  +  RedHake_13  +  SilverHake_14  +   Pollock_16  +  AtlanticTomcod_17  +  Hakesp_18  +   Redfishsp_23  +  AtlanticHalibut_30  +  TurbotGreenlandHalibut_31  +  AmericanPlaice_40  +    YellowtailFlounder_42  +  WinterFlounder_43  +  GulfStreamFlounder_44  +  NorthernSennet_46  +  WhiteMullet_47  +  StripedAtlanticWolffish_50  +  SpottedWolffish_51  +  NorthernWolffish_52  +  Wolffish_59  +  HerringAtlantic_60  +    Capelin_64  +  AtlanticMackerel_70  +  CrevalleJack_86  +    LongfinHake_112  +  FourbeardRockling_114  +  BlueWhiting_117  +  GreenlandCod_118  +  RosefishBlackBelly_123  +  BrillWindowpane_143  +    BarndoorSkate_200  +  ThornySkate_201  +  LittleSkate_203  +  WinterSkate_204  +   NorthernHagfish_241  +  Dogfish_274  +    Sculpinfamily_311  +   MonkfishGoosefishAngler_400  +  MarlinSpikeGrenadier_410  +  RoughheadGrenadier_411  +  RockGrenadierRoundnose_414  +  Seasnails_500  +  SeasnailGelatinous_505   +  SandLances_590  +  AmericanSandLance_599  +  AmericanEel_600  +  EelpoutNewfoundland_619  +  RadiatedShanny_625   +  CuskEels_660  +   AtlanticSauryNeedlefish_720  +  WhiteBarracudina_727  +    SeaPotato_1823  +  PandalusBorealis_2211  +  PandalusMontagui_2212 +  Argissp_2410  +  Crangonsp_2416  +  JonahCrab_2511 "


  fl = paste( yvar, "1", sep="~" )
  for (i in 1:length(depvars) ) fl = paste( fl, "+", depvars[i]  )
  Qlm = lm( formula=as.formula(fl), data=S )
  summary(Qlm)
  require(MASS)
  y = stepAIC( Qlm)
  fy = labels(terms(y))  # accepted terms 

  save( Qlm, file=file.path( project.datadirectory("snowcrab"), "R", "density.model.lm.rdata"), compress=T )
  save( y, file=file.path( project.datadirectory("snowcrab"), "R", "density.model.lm.stepwise.rdata"), compress=T )


  # gam model
  require(mgcv)
  fl = paste( yvar, "1", sep="~" )
     for ( i in 1:length(fy) ) fl = paste( fl, "+ s(", fy[i], ")" )
  
  
  fl = "totno.all~1 + s( t ) + s( tamp) + s( z ) + s( substrate.mean )  + s( AtlanticCod_10 ) + s( Haddock_11 ) + s( WhiteHake_12 ) + s( RedHake_13 ) + s( SilverHake_14 )  + s( Pollock_16 ) + s( AtlanticTomcod_17 ) + s( Hakesp_18 )  + s( Redfishsp_23 ) + s( AtlanticHalibut_30 ) + s( TurbotGreenlandHalibut_31 ) + s( AmericanPlaice_40 ) + s( YellowtailFlounder_42 ) + s( WinterFlounder_43 ) + s( GulfStreamFlounder_44 ) + s( NorthernSennet_46 ) + s( WhiteMullet_47 ) + s( StripedAtlanticWolffish_50 ) + s( SpottedWolffish_51 ) + s( NorthernWolffish_52 ) + s( Wolffish_59 ) + s( HerringAtlantic_60 ) +  s( Capelin_64 ) + s( AtlanticMackerel_70 ) + s( CrevalleJack_86 ) + s( LongfinHake_112 ) + s( FourbeardRockling_114 ) + s( BlueWhiting_117 ) + s( GreenlandCod_118 ) + s( RosefishBlackBelly_123 ) + s( BrillWindowpane_143 ) + s( BarndoorSkate_200 ) + s( ThornySkate_201 ) + s( LittleSkate_203 ) + s( WinterSkate_204 ) +  s( NorthernHagfish_241 ) + s( Dogfish_274 ) + s( Sculpinfamily_311 ) + s( MonkfishGoosefishAngler_400 ) + s( MarlinSpikeGrenadier_410 ) + s( RoughheadGrenadier_411 ) + s( RockGrenadierRoundnose_414 ) + s( Seasnails_500 ) + s( SeasnailGelatinous_505 ) + s( SandLances_590 )  + s( AmericanSandLance_599 ) + s( AmericanEel_600 )  + s( EelpoutNewfoundland_619 ) + s( RadiatedShanny_625 )   + s( CuskEels_660 )  + s( AtlanticSauryNeedlefish_720 ) + s( WhiteBarracudina_727 )  + s( SeaPotato_1823 ) + s( PandalusBorealis_2211 ) + s( PandalusMontagui_2212 )+ s( Argissp_2410 ) + s( Crangonsp_2416 ) + s( JonahCrab_2511 )"
  Qgam = gam( formula=as.formula(fl), data=S )
  summary(Qgam)
  plot( Qgam)
  

  save( Qgam, file=file.path( project.datadirectory("snowcrab"), "R", "density.model.gam.rdata", compress=T )



