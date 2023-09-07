################################################################################
## Water Maze Paper Codes
## Bianca Luedeker, Jacob Turner, Monnie McGee
## Created: 12/21/2022
## Date Updated: 8/21/2023
################################################################################

#### Packages ####
library(compositions)  ##analysis of compositional data from book
library(progress)      ## Creates progress bars for looped functions.
library(xtable)        ## exports tables and matrices in latex format.
library(sirt)          ##Estimates Dirichlet alpha parameters
library(gtools)  ##Generates all possible permutations, rdirichlet here
library(DirichletReg)  ## Compare two Dirichlet Models
library(DescTools)     ## Hotelling's T^2 test

#### Table 1 Functions and Simulations ####


##Necessary User Created Functions

## Log-likelihood function for four distinct alphas seen in alternative Hypothesis.
log.like.alt <- function(alphas, dat) {
  log.p.bar <- apply(log(dat), 2, mean)
  sum.alpha <- sum(alphas)
  result <- 7*log(gamma(sum.alpha)) - 7*sum(log(gamma(alphas))) + 7*sum((alphas-1)*log.p.bar)
  return(result)
}

## Log-likelihood function for the null of equal alphas.
log.like.null <- function(alpha, dat){
  sum.p.bar <- sum(apply(log(dat), 2, mean))
  result <- 7*log(gamma(4*alpha)) - 28*log(gamma(alpha)) + 7*(alpha-1)*sum.p.bar
  return(result)
}

## Derivative of null log-likelihood
deriv.null.log.like <- function(alpha, data){
  sum.p.bar<- sum(apply(log(data), 2, mean))
  result <- 28*digamma(4*alpha) - 28*digamma(alpha) + 7*sum.p.bar
  return(result)
}


## Water Maze Simulation Study Function.
## Input a vector of Dirichlet Parameters alpha and number of simulations n.
water.maze.sim <- function(n, alphas){
  sim.data <- array(NA, dim=c(14, 4)) ## Matrix to hold the raw simulated data.
  reject.both <- 0
  accept.both <- 0
  not.the.same <- 0
  pb <- progress_bar$new(total = n)
  for(i in 1:n){
    ##Create data and divide into two groups.
    sim.data <- rdirichlet(14, alphas)
    tg.data <- sim.data[1:7,]
    wild.data <- sim.data[8:14,]
    
    ## Looking for MLES under the NUll
    deriv.null.tg <- function(alpha) deriv.null.log.like(alpha, tg.data)
    deriv.null.wild <- function(alpha) deriv.null.log.like(alpha, wild.data)
    #where the maximums occur under the NULL
    ifelse (sign(deriv.null.tg(0.1))==sign(deriv.null.tg(100)), max.alpha.tg.null <- NA,
            max.alpha.tg.null <- uniroot(deriv.null.tg, c(0.1,100))$root )
    ifelse(sign(deriv.null.wild(0.1))==sign(deriv.null.wild(100)), max.alpha.wild.null <- NA,
           max.alpha.wild.null <- uniroot(deriv.null.wild, c(0.1, 100))$root )
    ## Maximums under the null.
    max.tg.null <- log.like.null(max.alpha.tg.null, tg.data)
    max.wild.null <- log.like.null(max.alpha.wild.null, wild.data)
    
    # Parameter MLES  under the alternative.
    ## Where the maxs occur with no restrictions (alternative hypothesis)
    max.alpha.tg.alt <- sirt::dirichlet.mle(tg.data)$alpha
    max.alpha.wild.alt <- sirt::dirichlet.mle(wild.data)$alpha
    ## Maximums under the alternative.
    max.tg.alt <- log.like.alt(max.alpha.tg.alt, tg.data)
    max.wild.alt <- log.like.alt(max.alpha.wild.alt, wild.data)
    
    
    ss.correction <- 3/(3+5.4*7^-1.4)  ##small sample size correction.
    LR.stat.tg <- 2*(max.tg.alt - max.tg.null)*ss.correction
    p.tg <- pchisq(LR.stat.tg, 3, lower.tail=FALSE)
    LR.stat.wild <- 2*(max.wild.alt - max.wild.null)*ss.correction
    p.wild <- pchisq(LR.stat.wild, 3, lower.tail=FALSE)
    
    if(p.tg < 0.05 & p.wild < 0.05){ reject.both <- reject.both + 1}
    else if(p.tg >= 0.05 & p.wild >= 0.05){ accept.both <- accept.both + 1}
    else {not.the.same <- not.the.same + 1 }
    
    results <- c(reject.both, accept.both, not.the.same)
    
    pb$tick()
    Sys.sleep(1 / n)
  }
  return(results)
}

## Table 1 Simulations
water.maze.data <- read.csv("water_maze_data.csv")
(data.alphas<- sirt::dirichlet.mle(water.maze.data[,2:5])$alpha) ##Estimate appropriate alpha values using sample data.

sim.alphas.1 <- c(9.8, 6.1, 5.4, 5.9)  ## Rounded values close to those seen in the data.  Precision 27.2
sim.alphas.2 <- c(6.8, 6.8, 6.8, 6.8 ) ## Uniform, precision 27.2
sim.alphas.3 <- c(9, 6, 6, 6)          ## Rounded similar to data
sim.alphas.4 <- c(8, 7, 7, 5)
sim.alphas.5 <- c(12, 5, 5, 5)  ## Most extreme in one value

set.seed(1983)
(water.maze.sim(10000, sim.alphas.1))  ## (5501, 663, 3836) (reject both, accept both, reject one/accept one)
set.seed(1983)
(water.maze.sim(10000, sim.alphas.2))  ## (26, 9033, 941) 
#20,000 single tests, reject = 26*2+941 = 993, type one = 993/20000 = 0.04965)
set.seed(1983)
(water.maze.sim(10000, sim.alphas.3))  ##(2307, 2675, 5018)
set.seed(1983)
(water.maze.sim(10000, sim.alphas.4))  ##(1812, 3269, 4919)
set.seed(1983)
(water.maze.sim(10000, sim.alphas.5))  ##(9918, 0, 82)

all.sim.results <- data.frame( 
  "alpha" = c(9.8, 6.8, 9, 8, 12), 
  "reject both" = c(5501, 26, 2307, 1812, 9918),
  "fail to reject both" = c(663, 9033, 2675, 3269, 0),
  "reject only one" = c(3836, 941, 5018, 4919, 82), 
  "P(type I error" = c(0.3836, 0.0941, 0.5018, 0.4919, 0.0082))
xtable(all.sim.results)





#############################################################3
# Dirichlet Testing for the Water Maze Data set from Maugard
#
#Reading in the Water Maze Data set from Maugard.
#The data can be found at the following links and is also reference in the Maugard paper:
#https://raw.githubusercontent.com/xuod/dirichlet/master/example/3Tg.csv
#https://raw.githubusercontent.com/xuod/dirichlet/master/example/wt.csv
## If downloading data, make sure your setwd() is set to the directory where these files reside
## Please keep these object names - they are used later in the code.
# Tg3 <- read.csv("Tg3.csv") 
# Wt <- read.csv("Wt.csv")

## Alternative method (grab data directly from GitHub without downloading to local machine)
site1 <- "https://raw.githubusercontent.com/xuod/dirichlet/master/example/3Tg.csv"
site2 <- "https://raw.githubusercontent.com/xuod/dirichlet/master/example/wt.csv"
Tg3 <- read.csv(site1)
Wt <- read.csv(site2)

#Scaling to proportions and combining
Tg3<-Tg3/100
Wt<-Wt/100
Full.data<-rbind(Tg3,Wt)
#Group Variable
Group<-factor(rep(c("Tg3","Wt"),each=7))

#Ensuring data is a true composition
Full.data<-DR_data(Full.data)

##########################################################################################
#Creating Figure 3
plot(acomp(Full.data), pch=20, col = c("red", "blue")[Group], main="3 Part Subcompositions for 3tg and Wild Groups")


##########################################################################################
#Performing two sample test for means of two Dirichlet samples with unequal precision parameters
#
#Estimating mean and precision by group (unrestricted case)
mod1 <- DirichletReg::DirichReg(Full.data ~ Group | Group, model="alternative")
#Estimating common mean, precisions still estimated by group (null case)
mod2 <- DirichletReg::DirichReg(Full.data ~ 1 | Group, model = "alternative")
#Computing LRT test statistic and p-value
anova(mod1,mod2)  #Test: 7.1491 p-value: .06729



##########################################################
#Type-I error Simulation code shown in Table 2
#Note: This code will produce type-I error estimates of 
#the test of equal means using the Dirichlet model and proposed LRT.
#We've provided the code for the first situation to keep
#the code as concise as possible for reviewers.
#
#Note: to replicate our exact simulations, the user should reset
#the seed for each scenario. The simulation time is quite long
#for each scenario due to running two Dirichlet regression calls
#which typically require 15 to 30 iterations to converge.

#Scenario 1
alpha=c(9.8,6.1,5.4,5.9)
ss=7
n.sims=10000

#Type-I Simulation

set.seed(1234)
type1<-c()  #Storing type-I error rates by sample size scenario
cnt<-0 #initializing count for # of rejections under H0
XX<-factor(rep(c("A","B"),each=ss))
for(i in 1:n.sims){
  YY<-DR_data(rdirichlet(2*ss,alpha))
  mod1 <- DirichletReg::DirichReg(YY ~ XX | XX, model="alternative")
  mod2 <- DirichletReg::DirichReg(YY ~ 1 | XX, model = "alternative")
  LRT<-2*(mod1$logLik-mod2$logLik)
  cnt<-cnt+ifelse(pchisq(LRT,df=3,lower.tail=F)<.05,1,0)  #Pulling LRT p-value and comparing to sig level.
}
type1<-cnt/n.sims
type1

#Run the above code for the remaining 4 scenarios listed below.  Change ss to different sample sizes.
alpha=c(6.8,6.8,6.8,6.8)
alpha=c(9,6,6,6)
alpha=c(8,7,7,5)
alpha=c(12,5,5,5)



#################################################################
# Power simulation as a function of sample size using Maugard
#estimates as reference parameters.
#
Tg3.mles<-dirichlet.mle(Full.data[1:7,]) 
Wt.mles<-dirichlet.mle(Full.data[8:14,])
#Note: these are the same estimates when using Dirichlet Reg
Tg3.mles  # means .301,.255,.216,.228, precision 41.678
Wt.mles   # means .423,.194,.181,.202, precision 27.025

#Parameter values for power sim
alpha1<-Tg3.mles$alpha
alpha2<-Wt.mles$alpha

#Simulation for empirical power
n<-c(5,7,10,15,20)
pow<-c()
set.seed(1234)

for(j in 1:length(n)){
cnt<-0 #initializing count for # of rejections under HA
XX<-factor(rep(c("A","B"),each=n[j]))
for(i in 1:n.sims){
  x1<-rdirichlet(n[j],alpha1)
  x2<-rdirichlet(n[j],alpha2)
  YY<-DR_data(rbind(x1,x2))
  mod1sim <- DirichletReg::DirichReg(YY ~ XX | XX, model="alternative")
  mod2sim <- DirichletReg::DirichReg(YY ~ 1 | XX, model = "alternative")
  LRT<-2*(mod1sim$logLik-mod2sim$logLik)
  cnt<-cnt+ifelse(pchisq(LRT,df=3,lower.tail=F)<.05,1,0)  #Pulling LRT p-value and comparing to sig level.
}
sim.pow<-cnt/n.sims
pow[j]<-sim.pow
}

pow.results<-data.frame(N=n,Power=pow)
pow.results
plot(pow.results)

###################################################################
#Test for difference in mean compositions using CLR log ratio transformations
#and performing Hotelling's T^2 test
#
#
wmclr <- clr(Full.data)
wmclr.df <- as.data.frame(wmclr)
wmclr.df$Group <- Group

HotellingsT2Test(cbind(TQ,AQ1,OQ) ~ Group, data=wmclr.df)






####################################################################
#LRT Test for the NDD model applied to Water Maze data
#
#


#Transforming to 3 sub-tree compositions based on the tree given 
#in Figure 6
#

#Subtree under Root
b_root<-cbind(apply(Full.data[,c(2,3)],1,sum),apply(Full.data[,c(1,4)],1,sum))
b_root<-DR_data(b_root)
#Subtree under N1
b_N1<-Full.data[,c(2,3)]
b_N1<-DR_data(b_N1)
#Subtree under N2
b_N2<-Full.data[,c(4,1)]
b_N2<-DR_data(b_N2)


#Computing LRT for each subtree
#XX<-factor(rep(c("WT","AD"),each=7))
mod_root <- DirichletReg::DirichReg(b_root ~ Group | Group, model="alternative")
mod_rootnull <- DirichletReg::DirichReg(b_root ~ 1 | Group, model = "alternative")
LRT_root<-2*(mod_root$logLik-mod_rootnull$logLik)
LRT_root

mod_N1 <- DirichletReg::DirichReg(b_N1 ~ Group | Group, model="alternative")
mod_N1null <- DirichletReg::DirichReg(b_N1 ~ 1 | Group, model = "alternative")
LRT_N1<-2*(mod_N1$logLik-mod_N1null$logLik)
LRT_N1

mod_N2 <- DirichletReg::DirichReg(b_N2 ~ Group | Group, model="alternative")
mod_N2null <- DirichletReg::DirichReg(b_N2 ~ 1 | Group, model = "alternative")
LRT_N2<-2*(mod_N2$logLik-mod_N2null$logLik)
LRT_N2

#Global LRT and p-value for NDD test
LRT<-LRT_root+LRT_N1+LRT_N2
LRT
pchisq(LRT,df=3,lower.tail=F)



#############################################
#Results derived from Appendix (Confidence Intervals)
#Table 5 and Table 6

#Function to compute intervals for the difference in 
#Dirichlet means.  x is a result object from using the 
#predict function on a Dirichlet Reg object that contains the mean and 
#precision estimates for each group.  
#x is a list with two objects mu and phi. mu is a 2xk
#matrix of mean estimates (each row corresponds to a group)
#phi corrsponds to a (2x1 matrix) containing the two precision estimates
#
#Function returns data frame with estimates, se and CI
#for the differences for each component.
confint.DirMeans<-function(x,n1,n2,level=0.95){
  #n1<-nrow(data1)
  #n2<-nrow(data2)
  #logbar1<-apply(log(data1),2,sum)
  #logbar2<-apply(log(data2),2,sum)
  k=ncol(x$mu)
  A1<-c(x$phi[1,])
  A2<-c(x$phi[2,])
  pi.hat1<-c(x$mu[1,])
  pi.hat2<-c(x$mu[2,])
  
  
  I11.g1<- diag(A1^2*trigamma(A1*pi.hat1[-k]),k-1,k-1)+A1^2*trigamma(A1*pi.hat1[k])
  I12.g1<- A1*pi.hat1[-k]*trigamma(A1*pi.hat1[-k])-A1*pi.hat1[k]*trigamma(A1*pi.hat1[k])
  #I12.g1<- digamma(A1*pi.hat1[-k])+A1*pi.hat1[-k]-digamma(A1*pi.hat1[k])+trigamma(A1*pi.hat1[k])-logbar1[-k]+logbar1[k]
  I22.g1<- sum(pi.hat1^2*trigamma(A1*pi.hat1))-trigamma(A1)
  Fisher1<- n1*rbind(cbind(I11.g1,I12.g1),c(I12.g1,I22.g1))
  
  I11.g2<- diag(A2^2*trigamma(A2*pi.hat2[-k]),k-1,k-1)+A2^2*trigamma(A2*pi.hat2[k])
  I12.g2<- A2*pi.hat2[-k]*trigamma(A2*pi.hat2[-k])-A2*pi.hat2[k]*trigamma(A2*pi.hat2[k])
  #I12.g2<- digamma(A2*pi.hat2[-k])+A2*pi.hat2[-k]-digamma(A2*pi.hat1[k])+trigamma(A2*pi.hat1[k])-logbar2[-k]+logbar2[k]
  I22.g2<- sum(pi.hat2^2*trigamma(A2*pi.hat2))-trigamma(A2)
  Fisher2<- n2*rbind(cbind(I11.g2,I12.g2),c(I12.g2,I22.g2))
  
  Sigma1<-solve(Fisher1)
  Sigma2<-solve(Fisher2)
  if(k>2){
  sum.cont<-matrix(c(rep(1,k-1),0),k,1)
  sigma1k<-t(sum.cont) %*% Sigma1 %*% sum.cont
  sigma2k<-t(sum.cont) %*% Sigma2 %*% sum.cont
  
  Var1<-c(diag(Sigma1)[-k],sigma1k)
  Var2<-c(diag(Sigma2)[-k],sigma2k)
  }
  
  if(k==2){
  Var1<-c(diag(Sigma1)[-k],diag(Sigma1)[-k])
  Var2<-c(diag(Sigma2)[-k],diag(Sigma2)[-k])
  }
  
  Var.dif=Var1+Var2
  
  int.results<-data.frame(Mean1=pi.hat1,SE1=sqrt(Var1),Mean2=pi.hat2,SE2=sqrt(Var2),Dif=pi.hat1-pi.hat2,SE.dif=sqrt(Var.dif))
  z.crit<- abs(qnorm((1-level)/2))
  int.results$L<-int.results$Dif-z.crit*int.results$SE.dif
  int.results$U<-int.results$Dif+z.crit*int.results$SE.dif
  
  return(int.results)
}

#Dirichlet Confidence interval
results<-predict(mod1,newdata=data.frame(Group=factor((c("Tg3","Wt")))),mu=T,phi=T)

DDfinal<-confint.DirMeans(results,n1=7,n2=7,level=1-.05/4)
rownames(DDfinal)<-colnames(Full.data)
DDfinal


#Function which takes a data frame with 4 columns to compute a Nested Dirichlet CI for the difference in means
#Columns 1-2: Population 1 Mean estimates and SE's corresponding to the branches of the tree that correspond
#             to the component of interested.
#Columns 3-4: Population 2 Mean estimates and SE's corresponding to the branches of the tree that correspond
#             to the component of interested.

confint.NDirMeans<-function(Estimates,level){
  n_nodes<-nrow(Estimates)
  mean1<-prod(Estimates[,1])
  mean2<-prod(Estimates[,3])
  Var1<-diag(Estimates[,2]^2,nrow=n_nodes,ncol=n_nodes)
  Var2<-diag(Estimates[,4]^2,nrow=n_nodes,ncol=n_nodes)
  grad.h1<-c()
  grad.h2<-c()
  for(i in 1:n_nodes){
    grad.h1[i]<-prod(Estimates[-i,1])
    grad.h2[i]<-prod(Estimates[-i,3])
  }
  grad.h1<-matrix(grad.h1,nrow=1)
  grad.h2<-matrix(grad.h2,nrow=1)
  SE1<- sqrt(grad.h1 %*%Var1 %*% t(grad.h1))
  SE2<- sqrt(grad.h2 %*%Var2 %*% t(grad.h2))
  
  Var.dif<-SE1^2+SE2^2
  int.results<-data.frame(Mean1=mean1,SE1,Mean2=mean2,SE2,Dif=mean1-mean2,SE.dif=sqrt(Var.dif))
  z.crit<- abs(qnorm((1-level)/2))
  int.results$L<-int.results$Dif-z.crit*int.results$SE.dif
  int.results$U<-int.results$Dif+z.crit*int.results$SE.dif
  
  return(int.results)
  
}


#Extracting means and SEs for each node of the nested Dirichlet tree for
#the Maurgard data. Note mod_root, mod_N1, mod_N1 are the Dirichlet regression objects 
#previously used to compute the LRT statistic.
mleN0<-predict(mod_root,newdata=data.frame(Group=factor((c("Tg3","Wt")))),mu=T,phi=T)
mleN1<-predict(mod_N1,newdata=data.frame(Group=factor((c("Tg3","Wt")))),mu=T,phi=T)
mleN2<-predict(mod_N2,newdata=data.frame(Group=factor((c("Tg3","Wt")))),mu=T,phi=T)

resultN0<-confint.DirMeans(mleN0,n1=7,n2=7,level=.95) #note level here is irrelevant
resultN1<-confint.DirMeans(mleN1,n1=7,n2=7,level=.95)
resultN2<-confint.DirMeans(mleN2,n1=7,n2=7,level=.95)

#Creating data frames "estimates" to compute final confidence intervals for each
#component.

branch.info<-list()
branch.info$TQ<-data.frame(rbind(resultN0[2,1:4],resultN2[2,1:4]))
branch.info$AQ1<-data.frame(rbind(resultN0[1,1:4],resultN1[1,1:4]))
branch.info$OQ<-data.frame(rbind(resultN0[1,1:4],resultN1[2,1:4]))
branch.info$AQ2<-data.frame(rbind(resultN0[2,1:4],resultN2[1,1:4]))




NDDfinal<-do.call(rbind,lapply(branch.info,confint.NDirMeans,level=1-.05/4))
NDDfinal


