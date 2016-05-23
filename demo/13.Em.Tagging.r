#Emera Tagging

dn = file.path(project.datadirectory("snowcrab"),'data','tagging','Emera')

a = dir(dn)
b = a[grep('tags',a)]
a = a[grep('meta',a)]
a = read.csv(file.path(dn,a),header=T)
out = NULL
for(i in b) {
	 h = read.csv(file.path(dn,i),header=T)
	 out = rbind(h,out)
}

names(out) = tolower(names(out))
out$dt = strptime(out$datetime, format= "%d/%m/%Y %H:%M")

#N Animals detects
	b= unique(out$transmitter) #58

#N Detects per anmimal
	a = aggregate(detectedid~transmitter,data=out,FUN=length)

#each animal

	with(out[which(out$transmitter==b[1]),])
for(i in 1:length(b)) {
	oo = out[which(out$transmitter==b[i]),]
	oo = oo[order(oo$dt),]
	ll = nrow(oo)
	oo$ms = c(0,oo$dt[2:(ll)]-oo$dt[1:(ll-1)])
 loadfunctions(c('utility','polygons'))
 makeMap(xlim=c(-60.3,-60.1),ylim=c(46.5,46.6))
	with(oo,plot(lon,lat,type='l',lwd=0.4))
title(paste(b[i],'N days =', round(as.numeric(max(oo$dt)-min(oo$dt)[1]),1)	))
savePlot(file.path(dn,paste('tag',b[i],".pdf",sep=".")),type='pdf')
}
