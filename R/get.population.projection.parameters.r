 
  # functions used for the forward projection model of snow crab
  # -----------------------
  
  get.population.projection.parameters = function( p ) {
    # parameters for fishable biomass calculations and projections
    
    assessment.year = p$current.assessment.year
    
    p$regions = c("cfanorth", "cfasouth", "cfa4x")
    p$regions.labels = c( "N-ENS", "S-ENS", "4X")
    p$nregions = length(p$regions)
    p$fisheryyears = c(1998:assessment.year)
    p$nodeyears = c( 1999:assessment.year )  # previous data seem strange
    p$nyears = length(p$nodeyears)
    p$nfisheryyears = length(p$fisheryyears)
    p$nodes = c(
      paste("imm", c(5:12), sep="."),
      paste("imm.sm", c(9:12), sep="."),
      paste("CC1to2", c(9:13), sep="."),
      paste("CC3to4", c(9:13), sep="."),
      paste("CC5",    c(9:13), sep=".")
    )
    p$nnodes = length(p$nodes)
    # map between variable name and category
    p$nodelabels = c( 
      paste("mi", c(5:12), sep=""),
      paste("mi", c(9:12), ".skip.moulter", sep=""),
      paste("ma", c(9:13), ".CC1to2", sep=""),
      paste("ma", c(9:13), ".CC3to4", sep=""),
      paste("ma", c(9:13), ".CC5", sep="")
    )
    return(p)
  }

