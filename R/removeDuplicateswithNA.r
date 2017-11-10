
removeDuplicateswithNA = function(x,cols = c('trip','set'),idvar='dt'){
	 f = do.call(rbind,lapply(split(x,x[,cols]), function(rms) rms[which(!is.na(rms[,idvar])),]))
	 return(f)
}
