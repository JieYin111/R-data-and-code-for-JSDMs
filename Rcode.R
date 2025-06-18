################################################################################################################################
## Improving multi-species habitat identification through species weighting assignment using Joint Species Distribution Model ##
################################################################################################################################
## authors: "Yunlei Zhang,Jie Yin,Yupeng Ji, Chongliang Zhang, Binduo Xu, Yiping Ren, Ying Xue"
## Corresponding author: Ying Xue
## Email: xueying@ouc.edu.cn 
## date: "May, 2025"
## require R-4.0.3

### More details on modelling and cross-validation can be found in the paper "Impact of life history stages on fish species interactions and spatio-temporal distribution"(Zhang et al., 2023).


rm(list = ls())
setwd("E:/ecological informatics")##Before running the code below, please set the working directory  
load("study_data.RData")  ##data

####Required packages & functions
library(boral)
library(coda)
library(foreach) 
library(doSNOW)  # for parallel computing
library(pacman)

Tjur= function(obs, pred){
  categories <- unique(obs)
  m1 <- mean(pred[which(obs == categories[1])], na.rm = T)
  m2 <- mean(pred[which(obs == categories[2])], na.rm = T)
  cod = abs(m2 - m1)
  names(cod) <- "D"
  return(cod)
} 

ptest <- function(obs,pred){
  if(is.na(pred[1]) ) return(NA) else{
    D_total = Tjur(c(obs), c(pred))
    A_total = pROC::auc(c(obs), c(pred))
    D_spec= sapply(1:ncol(obs),function(i) Tjur(obs[,i],pred[,i]) )
    freq = colSums(obs >0 )
    del= which(freq==0 | freq== nrow(obs))
    forcomp= 1:ncol(obs)
    if (length(del)>0) forcomp=forcomp[-del]
    A_spec= sapply(forcomp,function(i) pROC::auc(obs[,i],pred[,i]) )
    return(list(TjurD= D_total, AUC= A_total, Dspecies= D_spec,AUCspecies= A_spec, miss= del))
  }
}


################################################
#    model fit,   Stacked species distribution     
###############################################
##error data have been removed
envdata=envdata[ ,c("SBT","SBS","Depth","Distance")]
env= scale(envdata)
#env=envdata
X=cbind(env,env^2)
colnames(X)[5:8]= paste(colnames(X)[1:4],2,sep="")
summary(X)
  
set.seed(123)
mcmc.control <- list(n.burnin = 100, n.iteration =1000,  n.thin = 10) 
n=nrow(Y)
fakedistmat <- as.matrix(dist(1:n))
m2b <- boral(y = Y, X=X,family = "binomial", lv.control = list(num.lv = 2, type = "independent", distmat = fakedistmat),
             mcmc.control = mcmc.control,save.model = T)  #,calc.ics = TRUE  tweedie normal binomial
summary(m2b)
b2 <- fitted(m2b)$out
probability=list(b2)
summary(b2)

###weight
summary(weight)
weight_AUC=weight$weight_AUC
weight_prevalence=weight$weight_prevalence
weight_trophic=weight$weight_Trophic


Y1a=Y[,c(2,4,6,8,10,12,14,16,18,20)]
b5_a=b2[,c(2,4,6,8,10,12,14,16,18,20)]

###S1
EI_ob1=rowSums(Y1a);EI_p1=rowSums(b5_a)
EI_ob1=(EI_ob1-min(EI_ob1))/(max(EI_ob1)-min(EI_ob1))
EI_p1=(EI_p1-min(EI_p1))/(max(EI_p1)-min(EI_p1))
newww=data.frame(data[,1],EI_ob1,EI_p1)
colnames(newww)=c("year","Observed","Predicated")

groupcu_mean=newww %>% group_by(year) %>% summarise_all(mean, na.rm=TRUE)
groupcu_sd=newww %>% group_by(year) %>% summarise_all(sd, na.rm=TRUE)
qushi1=data.frame(groupcu_mean,groupcu_sd[,2:3],Scenario=rep("Scenario1",8))
colnames(qushi1)=c("year","Observed","Predicated","Observed_sd","Predicated_sd","Scenario")

n1=nrow(Y)
###S2,S3,S4
EI_ob22=matrix(NA,n1,10);EI_p22=matrix(NA,n1,10)
EI_ob33=matrix(NA,n1,10);EI_p33=matrix(NA,n1,10)
EI_ob44=matrix(NA,n1,10);EI_p44=matrix(NA,n1,10)
for(k in 1:10){
  EI_ob22[,k]=Y1a[,k]*weight_AUC[k]; EI_p22[,k]=b5_a[,k]*weight_AUC[k]
  EI_ob33[,k]=Y1a[,k]*weight_AUC[k]*weight_prevalence[k]; EI_p33[,k]=b5_a[,k]*weight_AUC[k]*weight_prevalence[k]
  EI_ob44[,k]=Y1a[,k]*weight_AUC[k]*weight_trophic[k]; EI_p44[,k]=b5_a[,k]*weight_AUC[k]*weight_trophic[k]
}

EI_ob2=rowSums(EI_ob22);EI_p2=rowSums(EI_p22)
EI_ob2=(EI_ob2-min(EI_ob2))/(max(EI_ob2)-min(EI_ob2))
EI_p2=(EI_p2-min(EI_p2))/(max(EI_p2)-min(EI_p2))
newww=data.frame(data[,1],EI_ob2,EI_p2)
colnames(newww)=c("year","Observed","Predicated")
groupcu_mean=newww %>% group_by(year) %>% summarise_all(mean, na.rm=TRUE)
groupcu_sd=newww %>% group_by(year) %>% summarise_all(sd, na.rm=TRUE)
qushi2=data.frame(groupcu_mean,groupcu_sd[,2:3],Scenario=rep("Scenario2",8))
colnames(qushi2)=c("year","Observed","Predicated","Observed_sd","Predicated_sd","Scenario")


EI_ob3=rowSums(EI_ob33);EI_p3=rowSums(EI_p33)
EI_ob3=(EI_ob3-min(EI_ob3))/(max(EI_ob3)-min(EI_ob3))
EI_p3=(EI_p3-min(EI_p3))/(max(EI_p3)-min(EI_p3))
newww=data.frame(data[,1],EI_ob3,EI_p3)
colnames(newww)=c("year","Observed","Predicated")
groupcu_mean=newww %>% group_by(year) %>% summarise_all(mean, na.rm=TRUE)
groupcu_sd=newww %>% group_by(year) %>% summarise_all(sd, na.rm=TRUE)
qushi3=data.frame(groupcu_mean,groupcu_sd[,2:3],Scenario=rep("Scenario3",8))
colnames(qushi3)=c("year","Observed","Predicated","Observed_sd","Predicated_sd","Scenario")

EI_ob4=rowSums(EI_ob44);EI_p4=rowSums(EI_p44)
EI_ob4=(EI_ob4-min(EI_ob4))/(max(EI_ob4)-min(EI_ob4))
EI_p4=(EI_p4-min(EI_p4))/(max(EI_p4)-min(EI_p4))
newww=data.frame(data[,1],EI_ob4,EI_p4)
colnames(newww)=c("year","Observed","Predicated")
groupcu_mean=newww %>% group_by(year) %>% summarise_all(mean, na.rm=TRUE)
groupcu_sd=newww %>% group_by(year) %>% summarise_all(sd, na.rm=TRUE)
qushi4=data.frame(groupcu_mean,groupcu_sd[,2:3],Scenario=rep("Scenario4",8))
colnames(qushi4)=c("year","Observed","Predicated","Observed_sd","Predicated_sd","Scenario")


###S5
EI_ob55=matrix(NA,n1,10);EI_p55=matrix(NA,n1,10)
for(k in 1:10){
  EI_ob55[,k]=Y1a[,k]*weight_AUC[k]*weight_prevalence[k]*weight_trophic[k];EI_p55[,k]=b5_a[,k]*weight_AUC[k]*weight_prevalence[k]*weight_trophic[k]
}


EI_ob5=rowSums(EI_ob55);EI_p5=rowSums(EI_p55)
EI_ob5=(EI_ob5-min(EI_ob5))/(max(EI_ob5)-min(EI_ob5))
EI_p5=(EI_p5-min(EI_p5))/(max(EI_p5)-min(EI_p5))
newww=data.frame(data[,1],EI_ob5,EI_p5)
colnames(newww)=c("year","Observed","Predicated")
groupcu_mean=newww %>% group_by(year) %>% summarise_all(mean, na.rm=TRUE)
groupcu_sd=newww %>% group_by(year) %>% summarise_all(sd, na.rm=TRUE)
qushi5=data.frame(groupcu_mean,groupcu_sd[,2:3],Scenario=rep("Scenario5",8))
colnames(qushi5)=c("year","Observed","Predicated","Observed_sd","Predicated_sd","Scenario")


my_datafr_qushi=rbind(qushi1,qushi2,qushi3,qushi4,qushi5)
OB=data.frame(my_datafr_qushi[,c(1,2,4,6)],variable=rep("Observed",length(my_datafr_qushi[,1])))
colnames(OB)=c("year","Mean_EI","SD_EI","Scenario","variable")
PRE=data.frame(my_datafr_qushi[,c(1,3,5,6)],variable=rep("Predicated",length(my_datafr_qushi[,1])))
colnames(PRE)=c("year","Mean_EI","SD_EI","Scenario","variable")
my_datafr_qushi1=rbind(OB,PRE)


####EI_pre_ob
newww1=data.frame(data[,1],EI_ob1,EI_p1,Scenario=rep("Scenario1",length(data[,1])))
colnames(newww1)=c("year","Observed","Predicated","Scenario")
newww2=data.frame(data[,1],EI_ob2,EI_p2,Scenario=rep("Scenario2",length(data[,1])))
colnames(newww2)=c("year","Observed","Predicated","Scenario")
newww3=data.frame(data[,1],EI_ob3,EI_p3,Scenario=rep("Scenario3",length(data[,1])))
colnames(newww3)=c("year","Observed","Predicated","Scenario")
newww4=data.frame(data[,1],EI_ob4,EI_p4,Scenario=rep("Scenario4",length(data[,1])))
colnames(newww4)=c("year","Observed","Predicated","Scenario")
newww5=data.frame(data[,1],EI_ob5,EI_p5,Scenario=rep("Scenario5",length(data[,1])))
colnames(newww5)=c("year","Observed","Predicated","Scenario")


EI_pre_ob=rbind(newww1,newww2,newww3,newww4,newww5)
OB=data.frame(EI_pre_ob[,c(1,2,4)],variable=rep("Observed",length(EI_pre_ob[,1])))
colnames(OB)=c("year","Mean_EI","Scenario","variable")
PRE=data.frame(EI_pre_ob[,c(1,3,4)],variable=rep("Predicated",length(EI_pre_ob[,1])))
colnames(PRE)=c("year","Mean_EI","Scenario","variable")
my_datafr_EI_pre_ob=rbind(OB,PRE)
