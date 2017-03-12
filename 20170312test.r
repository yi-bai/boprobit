#This program estimates a bivariate ordered probit model with each dependent variable having 
# 3 outcomes using R's optim function

library(MASS)
library(foreign)
library(mvtnorm)
library(ggplot2)

#Generating data to test the model
n <- 1000								
x1 <- rnorm(n)
x2 <- rnorm(n)
sig <- diag(2)
sig[1,2] <- sig[2,1] <- .5
mu <- c(0,0)

e <- rmvnorm(n,mu,sig)
a <- .5
b <- .2
c <- -.2
d <- -.9

y1 <- b*x1+e[,1]
y2 <- c+d*x2+e[,2]

y1 <- (y1 <= 0)*1 + (y1 > 0 & y1 <= 1)*2  + (y1 > 1)*3
y2 <- (y2 <= 0)*0 + (y2 > 0)*1

mydata <- cbind(y1,y2,x1,x2)
mydata <- as.data.frame(mydata)
#write.dta(mydata, file="/Users/tiernay/Desktop/probit/bioprobit_r.dta")

y1.start <- as.factor(y1)
y2.start <- as.factor(y2)

ggplot(mydata, aes(x = x2, y = y2)) + geom_point()



# The bivariate ordered probit function
log.lik <- function(par , X1 , X2 , Y1 , Y2) {
#X1 <- cbind(1,X1)
X1 <- as.matrix(X1)
X2 <- as.matrix(cbind(1,X2))

Beta1 <- par[ncol(X1)]		#parameters for the first equation
Beta2 <- par[(ncol(X1)+1):(ncol(X1)+2)]   #parameters for the second equation

cut1 <- par[ncol(X1)+3]	#Second cut point for the first equation
cut2 <- par[ncol(X1)+4]
#cut22 <- par[(ncol(X1)+ncol(X2)+2)]	#Second cut point for the second equation

gamma <- par[ncol(X1)+5]	  #parameter for the correlation coefficient

#multiply gamma by a column of 1's, and then transform: (exp(rho)-1)/(1+exp(rho))
rho <- (exp(gamma) - 1) / (1 + exp( gamma)) 

mu1 <- as.matrix(X1) %*% Beta1 
mu2 <- as.matrix(X2) %*% Beta2



llik <- 0 

for (i in 1:nrow(mu1)){
	cut_y1 <- c(-10000, cut1, cut2, 10000)
	cut_y2 <- c(-10000, 0, 10000)
	
	for
	
}

for (i in 1:nrow(mu1)){ 
	Sigma <- matrix(c(1, rho, rho, 1), 2, 2) 
  
	if (Y1[i]==0){	#All the calculations if y1 = 0				
		if (Y2[i]==0){    #Y1=0 and Y2=0
			llik <- llik + log(pmvnorm(lower = c(-Inf, -Inf), upper = c(0, 0), mean = c((cut1 - mu1[i,]), 0 - mu2[i,]), corr = Sigma))
		}else{	#Y1=0 and Y2=1
			llik <- llik + log(pmvnorm(lower = c(-Inf, 0), upper = c(Inf, Inf), mean = c((cut1 - mu1[i,]),0), corr = Sigma) - pmvnorm(lower = c(-Inf, -Inf), upper = c(Inf, 0), mean = c((cut1 - mu1[i,]), -mu2[i,]), corr = Sigma))
		}
	}else if (Y1[i]==1){	#All the calculations if y1 = 1	
		if (Y2[i]==0){  #y1=1 and y2=0
			llik <- llik + log(pmvnorm(lower = c(-Inf, -Inf), upper = c(Inf, 0), mean = c((cut2 - mu1[i,]), -mu2[i,]), corr = Sigma) - pmvnorm(lower = c(-Inf, -Inf), upper = c(Inf, 0), mean = c((cut1 - mu1[i,]), -mu2[i,]), corr = Sigma))	
		}else{	#y1=1 and y2=1
			llik <- llik + log(pmvnorm(lower = c(-Inf, 0), upper = c(Inf, Inf), mean=c((cut2 - mu1[i,]), 0), corr = Sigma) - pmvnorm(lower = c(-Inf, 0), upper = c(Inf, Inf), mean = c((cut2 - mu1[i,]), -mu2[i,]), corr = Sigma) - pmvnorm(lower = c(-Inf, -Inf), upper = c(Inf, Inf), mean = c((cut1 - mu1[i,]), 0), corr = Sigma) + pmvnorm(lower = c(-Inf, -Inf), upper = c(Inf, 0), mean = c((cut1 - mu1[i,]), -mu2[i,]), corr = Sigma))
		}					
	}else{		#All the calculations if y1 = 2
		if (Y2[i]==0){  #y1=2 and y2=0		
			llik <- llik +  log(pmvnorm(lower = c(-Inf, -Inf), upper = c(Inf,0),mean=c(0,-mu2[i,]), corr = Sigma) - pmvnorm(lower = c(-Inf, -Inf), upper = c(Inf, 0), mean = c((cut2 - mu1[i,]), -mu2[i,]), corr = Sigma)) 	
		}else{		#y1=2 and y2=1
			llik <- llik + log(1 - pmvnorm(lower = c(-Inf, -Inf),upper = c(Inf, Inf), mean = c((cut2 - mu1[i,]), 0),corr = Sigma) + pmvnorm(lower = c(-Inf, -Inf), upper = c(Inf, 0), mean = c((cut2 - mu1[i,]), -mu2[i,]), corr = Sigma))
		}
    }
	return(llik)}

}
 





# Generate starting values with an ordinary probit
op.result1 <- polr(y1.start ~ x1, method=c("probit"))
#op.result2 <- polr(y2.start ~ x2, method=c("probit"))
op.result2 <- glm(y2.start ~ x2, family=binomial(link="probit"), data=mydata)
start.val <- c(op.result1$coef[1],op.result2$coef,op.result1$zeta,0) 

res <- optim(start.val, log.lik, method = "L-BFGS-B", hessian = TRUE, control = list(fnscale = -1), X1 = x1, X2 = x2, Y1 = y1, Y2 = y2, lower = c(-Inf, -Inf, -Inf, -2, -2, -2), upper = c(Inf, Inf, Inf, 3, 3, 2)) 

print(res)
#res$par  # Return the parameters
#-solve(res$hessian) # Return the covariance matrix
#rho.convert(res$par[7])	 	# Return the correlation coefficint
	

