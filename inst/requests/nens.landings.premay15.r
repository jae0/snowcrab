#nens landings pre may 15 

a = logbook.db('logbook')

a = a[which(a$cfa=='cfanorth'),]

b = aggregate(landings~yr,data=a,FUN=sum)

y=2005:2015

x=data.frame(yr=NA,tot.l=NA,ear.l=NA)
for(i in 1:length(y)){
	r = a[which(a$yr==y[i]),]
	dd = paste(y[i],'05','16',sep="-")
	r = r[which(r$date.fished<dd),]
	g=0
	if(nrow(r)>0) {g = aggregate(landings~yr,data=r,FUN=sum)[2]}
	x[i,] = c(y[i],b[which(b$yr==y[i]),'landings'],g)
}

table.view(x)