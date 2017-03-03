# Testing BOP with fatigue data
# Final@2017-3-2

library(MASS)
library(maxLik)

#-----------------------
# import data
#-----------------------
mydata <- read.csv("D:/road_safety/fatigue/20170301data.csv")

#set.seed(33567)
sed.seed(56789)



#-----------------------
# Likelihood function
#-----------------------

BOP <- function(param, dat){
  
  # assign which part of data should be used in this model
  inj <- dat[,1]
  fatig <- dat[,2]
  x <- dat[,c(3:22)]
  z <- dat[,23:32]
  
  z <- cbind(1, z) # add an constant vector to Eq2
  
  # storage of parameters
  alpha <- as.matrix(param[1:ncol(x)]) # coefficients of Eq1
  bet <- as.matrix(param[ (ncol(x)+1) : (ncol(x)+11) ])  # coefficients of Eq2
  ramda <- param[ncol(x)+12] 
                        
  cut1 <- param[ncol(x)+13]
  cut2 <- param[ncol(x)+14]
  cut3 <- param[ncol(x)+15]                      
                          
  # generate common error
  miu <- rnorm(nrow(x),mean = 0, sd = 1) # generate common error term
  #print(x)
  xa  <- as.matrix(x)%*%alpha
                       
  inj_star <- xa + ramda*miu
                          
  zb  <- as.matrix(z)%*%bet
  fatig <- zb + miu
                          
  cut1_inj_star = cut1-inj_star
  cut2_inj_star = cut2-inj_star
  cut3_inj_star = cut3-inj_star
                          
  # probability of Eq1
  p11 = pnorm(cut1_inj_star)	# Injury==1
  p12 = pnorm(cut2_inj_star) - pnorm(cut1_inj_star)   # Injury == 2
  p13 = pnorm(cut3_inj_star) - pnorm(cut2_inj_star)   # Injury == 3
  p14 = 1 - pnorm(cut3_inj_star)    # Injury == 4
                            
  if(cut2 <= cut1){ p12 <- pnorm(-10000) }
  if(cut3 <= cut2){ p13 <- pnorm(-10000) }
                          
  ll_1 <- (inj==1)*log(p11) + (inj==2)*log(p12) + (inj==3)*log(p13) + (inj==4)*log(p14)  # loglikelihood function
  print(sum((inj==1)*log(p11)))                       
  #probability of Eq2
  p21 <- pnorm(fatig)
  ll_2 <- log(fatig*(p21) + (1-fatig)*(1-p21))
                          
  ll <- ll_1 + ll_2
  print(p13)
  
  return(-sum(ll))
}

#-------------------------
# Set starting value
#-------------------------

#xb_star <- iso_y + ins + lic_no + terrain_head + ope_y + light_ngt_y + light_ngt_n + age2635 + age3645 + age4655','age5665 + age6600 + gender_m + form_tag + form_head + form_side + driexp0002 + roadtype_exp + roadtype_u + x_m_fatig
#ins + terrain_head + ope_y + time0006 + time0708 + time1719 + light_ngt_y + light_ngt_n + roadtype_exp + roadtype_u

st_inj <- polr(as.factor(inj1) ~ iso_y + ins + lic_no + terrain_head + ope_y + light_ngt_y + light_ngt_n + age2635 + age3645 + age4655 + age5665 + age6600 + gender_m + form_tag + form_head + form_side + driexp0002 + roadtype_exp + roadtype_u + fatig, data=mydata, method = "probit")
stval_alpha <- as.matrix(st_inj$coefficients)
stval_cut <- as.matrix(st_inj$zeta)
#stval_cut <- c(0.1, 1, 1.3)

st_fatig <- glm(fatig ~ ins.1 + terrain_head.1 + ope_y.1 + time0006 + time0708 + time1719 + light_ngt_y.1 + light_ngt_n.1 + roadtype_exp.1 + roadtype_u.1, family=binomial(link="probit"), data=mydata)
stval_beta <- as.matrix(st_fatig$coefficients)
stval <- c(stval_alpha, stval_beta, -0.5, stval_cut)

#stval <- c(0.04, -0.27, 0.67, -0.1, -0.6, 0.1, 0.2, -0.1, -0.1, 0.1, 0.4, 0.6, -0.1, 0.2, 0.1, -0.1, 0.1, 0.8, -0.2, 0.9, -4, 0.5, -0.2, 0.1, 0.5, 0.2, -0.1, -0.2, 0.1, 0.4, -0.3, -0.2, 0.01,  1,  1.3)


#---------------------
# use optim to do MLE
#---------------------
result1 <- optim(stval, BOP, dat = mydata, method="BFGS", hessian = TRUE) # Estimate the model
+result1 <- optim(stval, BOP, dat=mydata, method='Nelder-Mead', hessian = TRUE)
#result2 <- maxLik(BOP, start = stval, dat=mydata)
print(result1)



#------------------
# display results
#------------------
df <- length(mydata$inj1) - length(result1$par) # degree of freedom
se <- sqrt(diag(abs(solve(result1$hessian))))  # standard error
t <- result1$par/se # t value
p <- (1-pt(abs(t),df))*1.96 # p value
rho <- (result1$par[32]/sqrt(2*(1+result1$par[32])^2))

display1 <- cbind(result1$par, se, t, p)
#display1(colnames) <- c('b', 'se', 't-value', 'p-value')
#display1(rownames) <- c('inj_x1', 'inj_x2', 'fat_con', 'z1', 'z2', 'ramda', 'cu1', 'cut2', 'cut3')

print(display1, digits = 3) # Displays the coefficients and standard errors of the parameters
print(rho)
