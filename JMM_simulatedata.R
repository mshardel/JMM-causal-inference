####################################
#12/16
#simulating data for joint longitudinal modeling of 
#outcome and exposure (shared parameter model),
#which will be used to estimate causal effects.
#Includes a baseline 'visit' and two follow-up visits
#############################################

#install.packages('mvtnorm')
library(mvtnorm)

#for saving simulated data
path=[insert path to save simulated data]

#number of participants
n <- 1000

#number of simulations
n.sim <- 1000 


#variance of random slope
var.b1.vec <- c(0, 0.3, 1, 3)

#variance of random intercept
var.b0 <- 1

#correlation of random slope and intercept
cor.b0.b1 <-0.5

#coefficient of random intercept or slope in exposure model
alpha.b.vec <- c(0,1)


###############outcome model:
#####outcome model (equation 3):
#####Y.t = beta0 + beta1*X + beta2*Y.(t-1)*(t>0) + beta30*A0*(t>0) + beta31*A1*(t>1)  + beta40*A0*X*(t>0)  + beta41*A1*X*(t>1) +
#####      beta5*t + b0 + b1*X*(1 + A0*(t>0) + A1*(t>1)) + error.t

####coefficients:
beta0 <- 1
beta1 <- -1
beta2 <- 1
beta30 <- 1
beta31 <- 2
beta40 <- 0
beta41 <- 0
beta5 <- -1


##############exposure model:
#####exposure model (equation 4):
#####logit(E[A.(t-1) | X, Y.(t-1), \bar{A.(t-2)}, b0, b1] = alpha0 + alpha1*X + alpha2*Y.(t-1)*(t>0) + alpha3*A0*(t>1) +   
#####                                                         alpha4*X*Y.(t-1)*(t>0) + alpha5*(t-1) + alpha6*b0 + alpha7*b1

alpha0 <- 0
alpha1 <- 0.5
alpha2 <- -0.75
alpha3 <- 1.0
alpha4 <- 0.25
alpha5 <- 0
######alpha6 and alpha7 depend on value from var.b1.vec and alpha.b.vec (see below)
#alpha6 <- alpha.b*(var.b1==0) #if no random slope, shared parameter is random intercept 
#alpha7 <- alpha.b*(var.b1>0)  #if random slope, shared parameter is random slope



for(var.b1 in var.b1.vec){     ##variance of random slope

for(alpha.b in alpha.b.vec){   ##coefficient of shared parameter in exposure model

#coefficient for random intercept and slope in exposure model
alpha6 <- alpha.b*(var.b1==0) #if no random slope, shared parameter is random intercept 
alpha7 <- alpha.b*(var.b1>0)  #if random slope, shared parameter is random slope




#initializing dataset
totaldata <- NULL


for(i in 1:n.sim){  ##start simulating

#baseline covariate: X ~ normal(0, 1)
X <- rnorm(n)
X <- X-mean(X) ##centering X


if(var.b1>0){
#random intercept and slope
b <- rmvnorm(n, mean=rep(0,2),sigma=matrix(c(var.b0,rep(cor.b0.b1*sqrt(var.b0*var.b1),2),var.b1),nrow=2,ncol=2)) 
b0 <- b[,1]
b1 <- b[,2]
}else{
#random intercept only  
b0 <- rnorm(n, mean=0,sd=sqrt(var.b0)) 
b1 <- rep(0,n)
}



####baseline outcome value: 
##from equation 3: Y0 = beta0 + beta1*X + b0 + b1*X + error.0
Y0 <- beta0 + beta1*X + b0 + b1*X + rnorm(n)


####baseline exposure value:
##from equation 4: logit(E[A0 | X, Y0, b0, b1]) =  alpha0 + alpha1*X + alpha2*Y0 + alpha4*X*Y0 + alpha6*b0 + alpha7*b1  
logit.A0 <- alpha0 + alpha1*X + alpha2*Y0 + alpha4*X*Y0 + alpha6*b0 + alpha7*b1  
A0 <- sapply(logit.A0,function(x){rbinom(1,1,prob=(exp(x)/(1 + exp(x))))})


####time 1 outcome value: 
##from equation 3: Y1 = beta0 + beta1*X + beta2*Y0 + beta30*A0 + beta40*A0*X +
##                      beta5*t + b0 + b1*X*(1 + A0) + error.1
Y1 <- beta0 + beta1*X + beta2*Y0 + beta30*A0 + beta40*A0*X + beta5 + b0 + b1*X*(1 + A0) + rnorm(n)         


#####time 1 exposure value: 
##from equation 4: logit(E[A1 | X, Y0, A0, b0, b1] = alpha0 + alpha1*X + alpha2*Y1 + alpha3*A0 +   
##                                                     alpha4*X*Y1 + alpha5 + alpha6*b0 + alpha7*b1
logit.A1 <- alpha0 + alpha1*X + alpha2*Y1 + alpha3*A0 + alpha4*X*Y1 + alpha5 + alpha6*b0 + alpha7*b1
A1 <- sapply(logit.A1,function(x){rbinom(1,1,prob=(exp(x)/(1 + exp(x))))})


####time 2 outcome value: 
##from equation 3: Y2 = beta0 + beta1*X + beta2*Y1 + beta30*A0 + beta31*A1 + beta40*A0*X + beta41*A1*X +
##                      beta5*2 + b0 + b1*X*(1 + A0 + A1) + error.2
Y2 <- beta0 + beta1*X + beta2*Y1 + beta30*A0 + beta31*A1 + beta40*A0*X + beta41*A1*X + beta5*2 + b0 + b1*X*(1 + A0 + A1) + rnorm(n)         





######output dataset into .csv file

sim <- rep(i,n) #simulation number
id <- c(1:n) #id for simulated participant

the.data <- cbind(sim,id,X, Y0, A0, Y1, A1, Y2) #dataset for current simulation
the.data <- data.frame(the.data)
totaldata <- rbind(totaldata,the.data)

} #end dataset




##filename=JMMsim_coef[coefficient of random effect in exposuremodel]_var[variance of random slope]
filename <- paste0("JMMsim_coef", round(alpha.b),"_var",var.b1,".csv")
write.csv(totaldata, paste0(path,filename),row.names=FALSE)

} #end coefficient of random effect in exposure model
} #end random slope

