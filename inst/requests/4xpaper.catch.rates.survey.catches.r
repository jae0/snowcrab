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
                
                out[i,'yr'] = b[i,'yr']
                out[i,'nblocks'] = nrow(unique(cbind(r[which(!is.na(r$landings)),c('lat','lon')])))
                out[i,'nsamples'] = nrow(r[which(!is.na(r)),])
                out[i,'meancpue'] = mean(r$cpue,na.rm=T)
                out[i,'sdcpue'] = sd(r$cpue,na.rm=T)/sqrt(length(r$cpue))
}
 
#plot the cpue
par(mar=c(5,4,4,5))
ylim=c(min(out$meancpue-out$sdcpue),max(out$meancpue+out$sdcpue))
with(out,plot(yr-1,meancpue,type='b',ylim=ylim,xlab='Year',ylab=expression(Fishery (kg / trap)))) #to put on same year as survey
with(out,arrows(x0=yr-1,y0=meancpue,y1=meancpue+sdcpue,length=0.05,angle=90))
with(out,arrows(x0=yr-1,y0=meancpue,y1=meancpue-sdcpue,length=0.05,angle=90))
 

# survey data
 k = snowcrab.db(DS = 'set.complete', p = p)

 l = aggregate(totmass.male.com~cfa+yr,data = k, FUN = mean)
 m = aggregate(totmass.male.com~cfa+yr,data = k, FUN = sd)
 n = aggregate(totmass.male.com~cfa+yr,data = k, FUN = length)
 names(m)[3]  = 'totmass.male.com.sd' 
 names(n)[3]  = 'totmass.male.com.n' 

a = merge(merge(l,m),n)
a$se = a$totmass.male.com.sd / sqrt(a$totmass.male.com.n)
b = a[which(a$cfa=='cfa4x'),]

b$upper = b$totmass.male.com +b$se
b$lower = b$totmass.male.com -b$se
b[which(b$yr<2004),c(3:8)]<-NA
par(mar=c(5,4,4,5),new=T)
with(b,plot(yr,totmass.male.com*1000,type='b',xlab='Year',ylim=c(0,max(b$upper,na.rm=T)*1000),yaxt='n',col='red',ylab=""))
with(b,arrows(x0=yr,y0=totmass.male.com*1000,y1=upper*1000,length=0.05,angle=90,col='red'))
with(b,arrows(x0=yr,y0=totmass.male.com*1000,y1=lower*1000,length=0.05,angle=90,col='red'))
axis(side=4,at=seq(round(min(b$lower,na.rm=T)*1000),round(max(b$upper,na.rm=T)*1000),100))
mtext(side=4,line=3,expression('Survey' (kg/km^2)))
legend('topleft',lty=c(1,1),pch=c(1,1),cex=0.8, bty='n',c('Fishery','Survey'),col=c('black','red'))
savePlot('4xfisherysurvey.png',type='png')
