
  fishery.model.jags = function( DS="", yr=NULL ) {

    warning( "This function is deprecated")

    out = NULL
    bugsdir= system.file( "bugs", package="bio.snowcrab" )
    if (DS=="biomass.dynamic" ) {
      out = file.path( bugsdir, "biomassdynamic.bugs" )
    }
    if (DS=="biomass.dynamic.illegal") {
      out = file.path( bugsdir, "biomassdynamic.illegal.bugs" )
    }
    if (DS=="delay.difference") {
      out = file.path( bugsdir, "delaydifference.bugs" )
    }
    if (DS=="delay.difference.illegal") {
      out = file.path( bugsdir, "delaydifference.illegal.bugs" )
    }
    if (DS=="biomass.dynamic.candidate" ) {
      out = file.path( bugsdir, paste("biomassdynamic_",yr,"_candidate.bugs",sep="") )
    }

    if (is.null(out)) out = file.path( bugsdir, DS )
    return(out)
  }



