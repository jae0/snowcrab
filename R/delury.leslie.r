#delury estimate
delury.leslie <-	function(x,estimate=c('delury','leslie'),day.or.week='week') {
dl = list(coef=1234)
if(nrow(x)>10) {
	#x$lncpue <- log(x$cpue)
	
		#if sum across week before average cpue
		
			if(day.or.week=='day')  y <- aggregate(cbind(landings,effort)~date.fished,data=x,FUN=sum)
			if(day.or.week=='week') y <- aggregate(cbind(landings,effort)~weekno,data=x,FUN=sum)
			y$cpue <- y$landings/y$effort
			y$lncpue <- log(y$cpue)
			y$sumef <- cumsum(y$effort)
			y$sumlan <- cumsum(y$landings)
			y$yr = unique(x$yr)
			y$area = unique(x$area)

		#yi <- max(x[which(x$date_fished<as.Date(paste(unique(x$yr),06,30,sep="-"))),'weekno'])
		
		if(estimate=='delury') {
		y$sumef <- cumsum(y$effort)
		delu <- summary(lm(lncpue~sumef,data=y))$coefficients
		with(y,plot(sumef,lncpue,xlab='Effort',ylab='ln(CPUE)'))#,main=paste('Delury',unique(x$yr),sep="-")
		abline(a=delu[1,1],b=delu[2,1],col='red')
		#abline(v=y[which(y$weekno==yi),'sumef'],col='red')
		delu <- round(delu, 5) # extract coefficients
	
	mtext(bquote(lncpue == .(delu[2,1])*sumU + .(delu[1,1])-- pslope*.(delu[2,4])), side=1,line=4,cex=0.7) # display equation
		
		return(list(coef = delu, data = y))
		}
		if(estimate=='leslie') {
		y$sumlan <- cumsum(y$landings)
		lesl <- summary(lm(lncpue~sumlan,data=y))$coefficients
		with(y,plot(sumlan,lncpue,xlab='Landings',ylab='ln(CPUE)'))#main=paste('Leslie',unique(x$yr),sep="-")
		abline(a=lesl[1,1],b=lesl[2,1],col='red')
		#abline(v=y[which(y$weekno==yi),'sumlan'],col='red')
		return(list(coef = lesl, data = y))
		}
	}
	return(dl)
	}
