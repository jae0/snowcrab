
  observer.data.formats = function( year=2003) {
    ## I had thought they were uniform in the data series, obviously they have changed
    ## ... this may be a good mechanism to switch, depending upon year
    ##.. the year would have to be added as another argument to this function
   if ( year %in% c(1996, 1999:2003)) {
     w = c(1, 8, 1, 2, 1, 9, 3, 1, 1, 1, 3, 1,
          4, 2, 4, 4, 1,
          1, 2, 1, 6, 1, 14, 1, 2, 1, 6, 3, 6,
          3, 3, 2, 6, 1, 2, 1, 1, 1, 10, 5, 23, 8, 18 )
     v = c("B1", "sdate", "d0", "trap", "dx1", "tripno", "crabno", "d4", "sex", "d0", "carapacewidth", "d5",
          "ndays", "trapno", "clawheightRight", "clawheightLeft", "shell_condition",
          "shell_condition_3M", "durometer", "dx2", "cfv", "dx3", "comment", "dx4", "ll", "dx5", "latitude", "dx6", "longitude",
          "d8", "depth", "clawlength", "weight", "d9", "zone", "d10", "data", "d11", "digits", "d12", "observer", "d13",
          "vessel" )
     keep = c( "sdate", "trap", "tripno", "crabno", "sex", "carapacewidth",
          "ndays", "trapno", "clawheightRight", "clawheightLeft", "shell_condition",
          "shell_condition_3M", "durometer", "cfv", "comment", "ll", "latitude", "longitude",
          "depth", "clawlength", "weight", "zone", "data", "digits", "observer",
          "vessel" )
     }

     if ( year %in% c(1997, 1998)) {
     w = c(9,1, 2, 1, 9, 3, 1, 1, 1, 3, 1,
          4, 2, 4, 4, 1,
          1, 2, 1, 6, 1, 14, 1, 2, 1, 6, 3, 6,
          3, 3, 2, 6, 1, 2, 1, 1, 1, 10, 5, 23, 8, 18 )
     v = c("sdate", "d0", "trap", "dx1", "tripno", "crabno", "d4", "sex", "d0", "carapacewidth", "d5",
          "ndays", "trapno", "clawheightRight", "clawheightLeft", "shell_condition",
          "shell_condition_3M", "durometer", "dx2", "cfv", "dx3", "comment", "dx4", "ll", "dx5", "latitude", "dx6", "longitude",
          "d8", "depth", "clawlength", "weight", "d9", "zone", "d10", "data", "d11", "digits", "d12", "observer", "d13",
          "vessel" )
     keep = c( "sdate", "trap", "tripno", "crabno", "sex", "carapacewidth",
          "ndays", "trapno", "clawheightRight", "clawheightLeft", "shell_condition",
          "shell_condition_3M", "durometer", "cfv", "comment", "ll", "latitude", "longitude",
          "depth", "clawlength", "weight", "zone", "data", "digits", "observer",
          "vessel" )
     }


    out = list( w=w, v=v, keep=keep )

    return (out)
   }


