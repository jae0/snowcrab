
  allometry.snowcrab = function( id, sex, redo=F ) {

    M = NULL
    fn = file.path( project.datadirectory("snowcrab"), "R", paste( sex, id, "rdata", sep=".") ) 
    
    if (!redo)  {
      if (file.exists(fn)) load(fn)
      return(M)
    }

    det = snowcrab.db( DS="det.initial" )
    
    if (sex == "male" ) {
      det = det[ which( det$sex==male) , ]
      det$log.cw = log( det$cw )
      det$log.mass = log( det$mass )
      det$log.chela = log( det$chela )
      
      if (id == "cw.mass") M = glm( log.mass ~ log.cw, data=det, na.action='na.exclude')
      if (id == "chela.mass" ) M = glm( log.mass ~ log.chela, data=det, na.action='na.exclude')
      if (id =="cw.chela.mat" ) M = glm( mat ~ log.cw + log.chela, data=det, na.action='na.exclude', family=binomial(link="logit") )
    }

    if (sex == "female" ) {
      det = det[ which( det$sex==female) , ]
      det$log.cw = log( det$cw )
      det$log.mass = log( det$mass )
      det$log.abdomen = log( det$abdomen )

      if (id == "cw.mass") M = glm( log.mass ~ log.cw, data=det, na.action='na.exclude')
      if (id == "abdomen.mass" ) M = glm( log.mass ~ log.abdomen, data=det, na.action='na.exclude')
      if (id =="cw.chela.mat" ) M = glm( mat ~ log.cw + abdomen, data=det, na.action='na.exclude', family=binomial(link="logit") )
    }
   
    save( M, file=fn, compress=T )
    return (M)
  }


