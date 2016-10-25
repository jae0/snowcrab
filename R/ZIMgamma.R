ZIMgamma = function(d,response='totmass.male.com'){

	d$nonzero = ifelse(d[response] > 0, 1, 0)

	formula1 = paste("nonzero ~ as.factor(yr)-1")
	formula2 = paste(response,"~ as.factor(yr)-1")


	m1 <- glm(formula1, data = d, family = binomial(link = logit))
	m2 <- glm(formula2, data = subset(d, nonzero == 1), family = Gamma(link = log))


	bin_coef <- plogis(coef(m1))
	gamma_coef <- exp(coef(m2))

	pred <- exp(log(bin_coef) + log(gamma_coef))

	return(pred)

}
