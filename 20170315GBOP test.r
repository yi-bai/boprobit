# This program estimates a bivariate binary-ordered probit model
# using R's optim function
# Final@2017-3-16

library(MASS)
library(foreign)
library(mvtnorm)
library(optimx)
#library(ggplot2)

set.seed(12315)

# If TEST_DATA = TRUE, use sample data for model testing
# If TEST_DATA = FALSE, use your own data for model testing
TEST_DATA = FALSE

if(TEST_DATA==TRUE){
	n  <- 5000							
	x1 <- rnorm(n)
	x2 <- rnorm(n)
	c1 <- rnorm(n) # variable for threshold
	sig <- diag(2)
	sig[1,2] <- sig[2,1] <- .5
	mu <- c(0,0)
	
	# set parameter
	e <- rmvnorm(n, mu, sig)
	a <- .5	        #cons for threshold
	b <- .2         #coef for y1
	c <- -.2        #cons for y2
	d <- -.9        #coef for y2
	tao <- -0.2     #coef for treshold

	# main function
	y1 <- b*x1 + e[,1]
	y2 <- c + d*x2 + e[,2]

	# add threshold function
	cut1 <- 0
	cut2 <- cut1 + exp(a + tao*c1)

	y1star <- (y1 <= cut1)*1 + (y1 > cut1 & y1 <= cut2)*2  + (y1 > cut2)*3
	y2star <- (y2 <= 0)*0 + (y2 > 0)*1

	mydata <- cbind(y1star, y2star, x1, x2, c1)
	mydata <- as.data.frame(mydata)

	y1.start <- as.factor(y1star)
	y2.start <- as.factor(y2star)

}else{
	# set your own dataset for estimation
	mydata <- read.csv("C:/Users/Bai/Documents/boprobit/20170301data.csv")
	#mydata <- read.csv("C:/Users/liyanyan/Desktop/matlab/fatigue/20170301data.csv")
	mydata <- mydata[1:2000,]  # small sample for testing real data
	
	column_offset = 1
	# Data preperation
	y1star <- mydata[,column_offset+1]    # inj
    y2star <- mydata[,column_offset+2]    # fatig
	x1 <- mydata[,c((column_offset+3):(column_offset+21))]   # variable for inj
	x2 <- mydata[,c((column_offset+22):(column_offset+31))]  # variable for fatig
	c1 <- mydata[,c((column_offset+3):(column_offset+6))]  #variable for threshold function

}

# calculate iteration time
# start time
iterate_time <<- 0


# The bivariate ordered probit function
log_lik <- function(par , X1 , X2 , Y1 , Y2, C1) {
	
	# iterate_time start
	start_time <- Sys.time()
	
	# dependent variables
	X1 <- as.matrix( X1 )		# dependent variables for y1
	X2 <- as.matrix(cbind(1,X2)) 	# depedent variables for y2
	C1 <- as.matrix(cbind(1,C1)) 	# dependent variables for thresholds

	# parameters for estimation
	# the order of parameters: c(b, c, d, gamma, a, tao)
	Beta1 <- par[1 : ncol(X1)]		#parameters for the y1
	Beta2 <- par[(ncol(X1) + 1) : (ncol(X1) + ncol(X2))]   #parameters for y2
	gamma <- par[ncol(X1) + ncol(X2) + 1]  # for calculate rho
	theta <- par[(ncol(X1) + ncol(X2) + 2) : (ncol(X1) + ncol(X2) + 6)]  # for thresholds
	
	#	multiply gamma by a column of 1's, and then transform: (exp(rho)-1)/(1+exp(rho))
	rho <- (exp(gamma) - 1) / (1 + exp( gamma)) 

	mu1 <- as.matrix(X1) %*% Beta1  # xb for y1
	mu2 <- as.matrix(X2) %*% Beta2	# xb for y2
	
	if(TEST_DATA==TRUE){
		cu1 <- as.matrix(C1) %*% theta  # xb for threshold funtion
		cut1 <- 0  
		cut2 <- cut1 + exp(cu1) 
	
		cut_y1 <- c(-10000, as.double(cut1), as.double(cut2), 10000)
		cut_y2 <- c(-10000, 0, 10000)
		
	}else{
	
		cu1 <- as.matrix(C1) %*% theta  # xb for threshold funtion
		cut1 <- 0
		cut2 <- cut1 + exp(cu1) 
		cut3 <- cut2 + exp(cu1)

		matrix_dim <- dim(cu1)[1]
		
		cut_y1 <- matrix(
			c(rep(-10000, matrix_dim), rep(cut1, matrix_dim), cut2, cut3, rep(10000, matrix_dim)),
			nrow = dim(cu1)[1],
			ncol = 5
			)
		cut_y2 <- c(-10000, 0, 10000)	
	}

	
	
	llik <- 0
	
	for (i in 1:nrow(mu1)){
		Sigma <- matrix(c(1, rho, rho, 1), 2, 2) 
		cut_left <- cut_y1[i, Y1[i]]
		cut_right <- cut_y1[i, Y1[i]+1]
		cut_down <- cut_y2[Y2[i]+1]
		cut_up <- cut_y2[Y2[i]+2]
		#print(pmvnorm(lower = c(cut_left, cut_down), upper = c(cut_right, cut_up), mean = c(mu1[i,], mu2[i,]), corr = Sigma))
		llik <- llik + log(pmvnorm(lower = c(cut_left, cut_down), upper = c(cut_right, cut_up), mean = c(mu1[i,], mu2[i,]), corr = Sigma))
	}
	
	#if(is.finite(llik)==TRUE){
	#	llik <- llik
	#}else{
	#	llik <- 1e+20
	#}
	
	time_taken <- Sys.time() - start_time
	
	if(iterate_time %% 1000 == 0){
    	print(llik)
	}
	iterate_time <<- iterate_time + 1
	print(iterate_time)
	print(time_taken)
	print(llik)
	print(cut_y1[1:10, ])
	
	return(llik)
}
 
# Generate starting values with an ordinary probit
if(TEST_DATA==TRUE){

	op.result1 <- polr(y1.start ~ x1, method = c("probit"))
	op.result2 <- glm(y2.start ~ x2, family = binomial(link = "probit"), data = mydata)
	# the order of parameters: c(b, c, d, gamma, a, tao)
	# try not to set parameter of threshold in 0
	start.val  <- c(op.result1$coef[1], op.result2$coef, 0.4, 0.01, 0.01) 
	
}else{
	
	# initial value for y1
	op.result1 <- polr(as.factor(inj1) ~ iso_y + ins + lic_no + terrain_head + ope_y + light_ngt_y + light_ngt_n + age2635 + age4655 + age5665 + age6600 + gender_m + form_tag + form_head + form_side + driexp0002 + roadtype_exp + roadtype_u + fatig, data=mydata, method = "probit")
	
	# initial value for y2
	op.result2 <- glm(fatig ~ ins.1 + terrain_head.1 + ope_y.1 + time0006 + time0708 + time1719 + light_ngt_y.1 + light_ngt_n.1 + roadtype_exp.1 + roadtype_u.1, family = binomial(link = "probit"), data = mydata)
	
	# initial value for threshold
	thres <- c(rep(0.01, (ncol(c1)+1))) 
	
	# the order of parameters: c(y1, y2, gamma, thres)
	start.val <- c(op.result1$coef, op.result2$coef, 0.01, thres)
}

if(TEST_DATA==TRUE){

	lower_1 = c(rep(-Inf,6))
	upper_1 = c(rep(10,6))

}else{

	lower_1 = c(rep(-Inf, 28), -2, rep(-Inf, (ncol(c1)+1)))
	upper_1 = c(rep(Inf, 28), 3, rep(Inf, (ncol(c1)+1)))

}


res <- optim(start.val, log_lik, method = "L-BFGS-B", hessian = TRUE, control = list(fnscale = -1, ndeps= rep(1e-1, 36)), X1 = x1, X2 = x2, Y1 = y1star, Y2 = y2star, C1 = c1, lower = lower_1, upper = upper_1) 

#print(res)
#gamma <- res$par[5]
#rho <- (exp(res$par[6]) - 1) / (1 + exp(res$par[6]))
#print(rho)
#res$par  # Return the parameters
#-solve(res$hessian) # Return the covariance matrix
#rho.convert(res$par[7])	 	# Return the correlation coefficint
	

