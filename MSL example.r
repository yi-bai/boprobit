set.seed(10234)
nobs <- 1000

x1 <- rnorm(nobs)*.15^.5
x2 <- rnorm(nobs)*.35^.5
z <- rnorm(nobs)*.25^.5
y <- round(runif(nobs, 1,5), 0)
x <- cbind(x1, x2)

#### Generate Halton Sequences
library("randtoolbox")
R <- 200
a <- halton(n=nobs, dim=R, normal=T, init=T)

# Likelihood for 5 category ordered probit
llk.oprobit5 <- function(param, x, y) {
    # preliminaries
    x <- as.matrix(x)
    os <- rep(1, nrow(x))
    x <- cbind(os, x)
    b <- param[1:ncol(x)]
    t2 <- param[(ncol(x)+1)]
    t3 <- param[(ncol(x)+2)]
    t4 <- param[(ncol(x)+3)]
    sigma_a <- param[ncol(x)+4]

    # probabilities and penalty function
    xb <- x %*% b %*% rep(1, R)
    asigma <- a*sigma_a

    p1 <- pnorm(-xb-asigma)
    p1 <- log(apply(p1, MARGIN=1, FUN=sum)/R)

    if (t2 <= 0) {
        p2 <- -(abs(t2) * 10000)    # penalty function to keep t2>0
    } else {
        p2 <- pnorm(t2-xb-asigma)-pnorm(-xb-asigma)
        p2 <- log(apply(p2, MARGIN=1, FUN=sum)/R)
    }
    if (t3 <= t2) {
        p3 <- -((t2-t3)*10000)    # penalty to keep t3>t2
    } else {
        p3 <- pnorm(t3-xb-asigma)-pnorm(t2-xb-asigma)   
        p3 <- log(apply(p3, MARGIN=1, FUN=sum)/R)
    }
    if (t4 <= t3) {
        p4 <- -((t3-t4)*10000)
    } else {
        p4 <- pnorm(t4-xb-asigma)-pnorm(t3-xb-asigma) 
        p4 <- log(apply(p4, MARGIN=1, FUN=sum)/R)
    }
    p5 <- 1 - pnorm(t4-xb-asigma)
    p5 <- log(apply(p5, MARGIN=1, FUN=sum)/R)

    # -1 * log likelihood (optim is a minimizer)
    -sum(cbind(y==1,y==2,y==3,y==4, y==5) * cbind(p1,p2,p3,p4,p5))
}

# Use optim directly
ls.result <- lm(y~x)                    # use ls estimates as starting   values
stval <- c(ls.result$coefficients,1,2,3,2)  # initial guesses

oprobit.result <- optim(stval, llk.oprobit5, method="BFGS", x=x, y=y, hessian=T, control = list(trace = 10, REPORT = 1))
