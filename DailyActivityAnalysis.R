library(compositions)
library(coda.base)
library(DirichletReg)
View(bmi_activity)


## Cleaning up Activity Data Set 
## Healthy have zBMI<1.03, Obese defined to have zBMI>1.96
bmi_activity2<-bmi_activity[-which(bmi_activity$zBMI>1.03 & bmi_activity$zBMI<1.96),]
bmi_activity2<-bmi_activity2[which(bmi_activity2$gender=="girl"),]
dim(bmi_activity2)

## Final Cleaned Data 
bmi<-data.frame(bmi_activity2[,3:7],BMI=ifelse(bmi_activity2$zBMI>1.03,"Obese","Healthy"))
names(bmi)<-c("Sleep","Sedent","Lpa","Mpa","Vpa","BMI")

## Visualizing
plot(acomp(bmi[,-6]), pch=20, col = c("red", "blue")[bmi$BMI], main="")
plot(acomp(bmi[,1:5]),pch=ifelse(bmi$BMI=="Obese",3,20),col=ifelse(bmi$BMI=="Healthy","red","blue"))
plot(acomp(bmi[,1:5]),cex=.5,pch=ifelse(bmi_activity$gender=="girl",20,3),col=ifelse(bmi_activity$gender=="girl","red","blue"))

## Summary Statistics
apply(bmi[,1:5],2,mean)
apply(bmi[,1:5],2,sd)

apply(bmi[which(bmi$BMI=="Obese"),-6],2,mean)
apply(bmi[which(bmi$BMI=="Obese"),-6],2,sd)

apply(bmi[which(bmi$BMI=="Healthy"),-6],2,mean)
apply(bmi[which(bmi$BMI=="Healthy"),-6],2,sd)


cor(bmi[,-6])
cor(bmi[which(bmi$BMI=="Obese"),-6])
cor(bmi[which(bmi$BMI=="Healthy"),-6])

library(GGally)
ggpairs(bmi,aes(color=BMI))
ggpairs(bmi_activity[,3:8],aes(color=bmi_activity$gender))


## Performing a fit using a standard Dirichlet
Y<-DR_data(bmi[,-6])
BMI<-bmi$BMI
mod1 <- DirichletReg::DirichReg(Y ~ BMI | BMI, model="alternative")
#Estimating common mean, precisions still estimated by group (null case)
mod2 <- DirichletReg::DirichReg(Y~ 1 | BMI, model = "alternative")
#Computing LRT test statistic and p-value
anova(mod1,mod2)  #Test: 34.408 p-value: 6.146e-07

## Hotellings T-squared using CLR transform
library(DescTools)
bmiclr <- clr(DR_data(bmi[,1:5]))
bmiclr.df <- as.data.frame(bmiclr)
bmiclr.df$Group <- BMI

HotellingsT2Test(cbind(Sleep,Sedent,Lpa,Mpa) ~ Group, data=bmiclr.df)
## T2 stat:  9.61 , p-value 5.97e-07


## Dirichlet Confidence interval

## Definining the function
## x is a list containing the mle estimates for the mean and precision parameters
## each element in the list is a matrix with two rows (Row 1 group A, row 2 Group B)
## n1 and n2 are sample sizes for Group A and Group B

confint.DirMeans<-function(x,n1,n2,level=0.95){
  k=ncol(x$mu)
  A1<-c(x$phi[1,])
  A2<-c(x$phi[2,])
  pi.hat1<-c(x$mu[1,])
  pi.hat2<-c(x$mu[2,])
  
  
  I11.g1<- diag(A1^2*trigamma(A1*pi.hat1[-k]),k-1,k-1)+A1^2*trigamma(A1*pi.hat1[k])
  I12.g1<- A1*pi.hat1[-k]*trigamma(A1*pi.hat1[-k])-A1*pi.hat1[k]*trigamma(A1*pi.hat1[k])
  I22.g1<- sum(pi.hat1^2*trigamma(A1*pi.hat1))-trigamma(A1)
  Fisher1<- n1*rbind(cbind(I11.g1,I12.g1),c(I12.g1,I22.g1))
  
  I11.g2<- diag(A2^2*trigamma(A2*pi.hat2[-k]),k-1,k-1)+A2^2*trigamma(A2*pi.hat2[k])
  I12.g2<- A2*pi.hat2[-k]*trigamma(A2*pi.hat2[-k])-A2*pi.hat2[k]*trigamma(A2*pi.hat2[k])
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




## Computing the DD CI's
results<-predict(mod1,newdata=data.frame(BMI=factor((c("Healthy","Obese")))),mu=T,phi=T)

DDfinal<-confint.DirMeans(results,n1=135,n2=21,level=1-.05/5)
rownames(DDfinal)<-colnames(bmi[,1:5])
DDfinal

## Computing NDD fit using tree described in manuscript
###################################################
## Subtree under Root
b_root<-cbind(apply(Y[,1:2],1,sum),apply(Y[,3:5],1,sum))
b_root<-DR_data(b_root)
## Subtree under N1
b_N1<-Y[,c(1,2)]
b_N1<-DR_data(b_N1)
## Subtree under N2
b_N2<-cbind(Y[,3],apply(Y[,4:5],1,sum))
b_N2<-DR_data(b_N2)
## Subtree under N3
b_N3<-Y[,4:5]
b_N3<-DR_data(b_N3)
Group2<-BMI

## Computing LRT for each subtree
mod_root <- DirichletReg::DirichReg(b_root ~ Group2 | Group2, model="alternative")
mod_rootnull <- DirichletReg::DirichReg(b_root ~ 1 | Group2, model = "alternative")
LRT_root<-2*(mod_root$logLik-mod_rootnull$logLik)
LRT_root

mod_N1 <- DirichletReg::DirichReg(b_N1 ~ Group2 | Group2, model="alternative")
mod_N1null <- DirichletReg::DirichReg(b_N1 ~ 1 | Group2, model = "alternative")
LRT_N1<-2*(mod_N1$logLik-mod_N1null$logLik)
LRT_N1

mod_N2 <- DirichletReg::DirichReg(b_N2 ~ Group2 | Group2, model="alternative")
mod_N2null <- DirichletReg::DirichReg(b_N2 ~ 1 | Group2, model = "alternative")
LRT_N2<-2*(mod_N2$logLik-mod_N2null$logLik)
LRT_N2

mod_N3 <- DirichletReg::DirichReg(b_N3 ~ Group2 | Group2, model="alternative")
mod_N3null <- DirichletReg::DirichReg(b_N3 ~ 1 | Group2, model = "alternative")
LRT_N3<-2*(mod_N3$logLik-mod_N3null$logLik)
LRT_N3




## Global LRT and p-value for NDD test
LRT<-LRT_root+LRT_N1+LRT_N2+LRT_N3
LRT #46.32877
pchisq(LRT,df=4,lower.tail=F) #p-value: 2.103833e-09

## Confidence interval for the NDD model
## Defining the function.
## "estimates" input is a data frame with 4 colums
##  First 2 are branch means and standard errors for a single NDD variable for group A 
##  Second 2 are branch means and standard errors for a single NDD variable for group B
##  Computes confidence interval for the difference in means using the delta method
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







## Computing NDD CI's
mleN0<-predict(mod_root,newdata=data.frame(Group2=factor((c("Healthy","Obese")))),mu=T,phi=T)
mleN1<-predict(mod_N1,newdata=data.frame(Group2=factor((c("Healthy","Obese")))),mu=T,phi=T)
mleN2<-predict(mod_N2,newdata=data.frame(Group2=factor((c("Healthy","Obese")))),mu=T,phi=T)
mleN3<-predict(mod_N3,newdata=data.frame(Group2=factor((c("Healthy","Obese")))),mu=T,phi=T)


resultN0<-confint.DirMeans(mleN0,n1=135,n2=21,level=.95) #note level here is irrelevant
resultN1<-confint.DirMeans(mleN1,n1=135,n2=21,level=.95) #note level here is irrelevant
resultN2<-confint.DirMeans(mleN2,n1=135,n2=21,level=.95) #note level here is irrelevant
resultN3<-confint.DirMeans(mleN3,n1=135,n2=21,level=.95) #note level here is irrelevant


resultN0
resultN1
resultN2
resultN3

## Creating data frames "estimates" to compute final confidence intervals for each
## component.

branch.info<-list()
branch.info$Sleep<-data.frame(rbind(resultN0[1,1:4],resultN1[1,1:4]))
branch.info$Sedent<-data.frame(rbind(resultN0[1,1:4],resultN1[2,1:4]))
branch.info$Lpa<-data.frame(rbind(resultN0[2,1:4],resultN2[1,1:4]))
branch.info$Mpa<-data.frame(rbind(resultN0[2,1:4],resultN2[2,1:4],resultN3[1,1:4]))
branch.info$Vpa<-data.frame(rbind(resultN0[2,1:4],resultN2[2,1:4],resultN3[2,1:4]))

## Printing final results
NDDfinal<-do.call(rbind,lapply(branch.info,confint.NDirMeans,level=1-.05/5))
NDDfinal

