
  rename.snowcrab.variables = function(oldnames) {
    
    lookupnames = snowcrab.lookupnames()
    oldnames = tolower( oldnames )
    newnames = oldnames
    nvars = length(oldnames)
    for (i in 1:nvars) {
      j = which(lookupnames[,1] == oldnames[i])
      if (length(j)>0) newnames[i] = lookupnames[j,2]
    }

   return (newnames)
  }
	

