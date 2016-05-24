loadfunctions( "snowcrab", functionname="initialise.local.environment.r")
a = logbook.db('logbook.filtered.positions')
 
ii = which(a$cfa =='cfa4x')
a = a[ii,]
a$datf = strptime(a$date.fished,format='%Y-%m-%d')
ij = which(is.na(a$datf))
a$datf[ij] = strptime(a$date.fished[ij],format='%Y%m%d')
b= aggregate(date.fished~yr,data = a, FUN=min)
b$datf = strptime(b$date.fished,format='%Y-%m-%d')
 
out = data.frame(yr=NA,nblocks=NA,nsamples=NA,meancpue=NA,sdcpue=NA)
for(i in 1:nrow(b)) {
                r = a[which(a$yr == b[i,'yr']),]
                r = r[which(r$datf>= b$[i,'datf'] & r$datf <= 60*60*24*12),] #12 days 
                out[i,'yr'] = b[i,'yr']
                out[i,'nblocks'] = nrow(unique(cbind(r[which(!is.na(r$landings)),c('lat','lon')])))
                out[i,'nsamples'] = nrow(r[which(!is.na(r)),])
                out[i,'meancpue'] = mean(r$cpue,na.rm=T)
                out[i,'sdcpue'] = sd(r$cpue,na.rm=T)
}
 
#plot the cpue
ylim=c(min(out$meancpue-out$sdcpue),max(out$meancpue+out$sdcpue))
with(out,plot(yr,meancpue,type='b',ylim=ylim,xlab='Year',ylab='CPUE'))
with(out,arrows(x0=yr,y0=meancpue,y1=meancpue+sdcpue,length=0.05,angle=90))
with(out,arrows(x0=yr,y0=meancpue,y1=meancpue-sdcpue,length=0.05,angle=90))
 
#plot the area of the fishery
 
with(out,plot(yr,nblocks,xlab='Year',ylab='Number of grid cells fished',type='b'))
 