#script to construct a logistic regression model to predict reaction essentiality based on centrality metrics

#load libraries
rm(list=ls())
library(caret)
library(cowplot)
library(devtools)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(readr)
library(readxl)
library(rms)
library(stats)
library(tibble)
library(tidyr)

#load data with reaction essentiality and centrality metric information
df<-read_csv("LogisticRegr_Data.csv")
df$is_essential<-as.factor(df$is_essential) #convert outcome variable to factor

####Variable Selection####
#Use a 4/5, 1/5 split of the data for training and test set
#later can extend this by doing k-fold cross-validation
#First, randomly shuffle the data
shuffled<-df[sample(nrow(df)),]
#Create 5 equally size folds
folds <- cut(seq(1,nrow(shuffled)),breaks=5,labels=FALSE) #breaks data into folds of 250 when breaks=5
#Use the first fold as test data and the other four folds as training data
shuffled<-shuffled %>% add_column(folds)
test<-shuffled %>% filter(folds==1)
train<-shuffled %>% filter(folds>1)




####Code chunk to perform stepwise model selection####
{
#retrain the model on training set
mod0 = glm(is_essential~log10_bridging_centrality, data=train, family="binomial")
mod2 = glm(is_essential~log10_bridging_centrality + 
             log10_betweeness_centrality +
             clustering_coefficient +
             log10_degree +
             cascadeNum, data=train, family = "binomial")
stepResults<-step(mod0,scope=list(lower=mod0,upper=mod2),direction=c("both"))
#looks like the final model again contains all centrality metrics
modFinal<-glm(formula = is_essential ~ log10_bridging_centrality + cascadeNum + 
                log10_betweeness_centrality + log10_degree, data=train, family="binomial")
summary(modFinal) #look at coefficients and p values

#now run the test data through this final model, and get the probability that each reaction is essential
probabilities<-predict(modFinal, newdata=test, type="response")

#0.5 cutoff
{#cutoff is set at 0.5 for whether or not a reaction is essential (1) or nonessential (0)
predicted.classes <- ifelse(probabilities > 0.5, "1", "0")
tab<-data.frame("reaction" = test$reaction, "observed_essentiality" = test$is_essential, "predicted_essentiality" = predicted.classes)

####Performance evaluation####
tp<-length((tab %>% filter(observed_essentiality==1) %>% filter(predicted_essentiality==1))$reaction) #true positives
tn<- length((tab %>% filter(observed_essentiality==0) %>% filter(predicted_essentiality==0))$reaction) #true negatives
fp<- length((tab %>% filter(observed_essentiality==0) %>% filter(predicted_essentiality==1))$reaction)#false positives
fn<- length((tab %>% filter(observed_essentiality==1) %>% filter(predicted_essentiality==0))$reaction)#false negatives
t<-c(tp,fp,fn,tn)
#performance metric functions
recall=function(t){t[1]/(t[1]+t[3])}
specificity=function(t){t[4]/(t[2]+t[4])}
fpr=function(t){1-specificity(t)}
precision=function(t){t[1]/(t[1]+t[2])}
fdr=function(t){1-precision(t)}
#performance metric calculations
val_recall<-recall(t)
val_specificity<-specificity(t)
val_fpr<-fpr(t)
val_precision<-precision(t)
val_fdr<-fdr(t)
performance<-data.frame("Recall" = val_recall, "Specificity" = val_specificity,
                        "FPR" = val_fpr, "Precision" = val_precision, "FDR" = val_fdr)
#construct plots that show performance metrics
PRC<-ggplot()+geom_point(data=performance, aes(x=Recall, y=Precision), col="red", fill="red", size=4) + 
  theme_classic()+xlim(c(0,1))+ylim(c(0,1))
ROC<-ggplot()+geom_point(data=performance, aes(x=FPR, y=Recall), col="red", fill="red", size=4) + 
  theme_classic()+xlim(c(0,1))+ylim(c(0,1))+ylab("TPR")
}

#0.4 cutoff
{predicted.classes <- ifelse(probabilities > 0.4, "1", "0")
  tab<-data.frame("reaction" = test$reaction, "observed_essentiality" = test$is_essential, "predicted_essentiality" = predicted.classes)
  
  ####Performance evaluation####
  tp<-length((tab %>% filter(observed_essentiality==1) %>% filter(predicted_essentiality==1))$reaction) #true positives
  tn<- length((tab %>% filter(observed_essentiality==0) %>% filter(predicted_essentiality==0))$reaction) #true negatives
  fp<- length((tab %>% filter(observed_essentiality==0) %>% filter(predicted_essentiality==1))$reaction)#false positives
  fn<- length((tab %>% filter(observed_essentiality==1) %>% filter(predicted_essentiality==0))$reaction)#false negatives
  t<-c(tp,fp,fn,tn)
  #performance metric functions
  recall=function(t){t[1]/(t[1]+t[3])}
  specificity=function(t){t[4]/(t[2]+t[4])}
  fpr=function(t){1-specificity(t)}
  precision=function(t){t[1]/(t[1]+t[2])}
  fdr=function(t){1-precision(t)}
  #performance metric calculations
  val_recall<-recall(t)
  val_specificity<-specificity(t)
  val_fpr<-fpr(t)
  val_precision<-precision(t)
  val_fdr<-fdr(t)
  performance<-data.frame("Recall" = val_recall, "Specificity" = val_specificity,
                          "FPR" = val_fpr, "Precision" = val_precision, "FDR" = val_fdr)
  #construct plots that show performance metrics
  PRC<-ggplot()+geom_point(data=performance, aes(x=Recall, y=Precision), col="red", fill="red", size=4) + 
    theme_classic()+xlim(c(0,1))+ylim(c(0,1))
  ROC<-ggplot()+geom_point(data=performance, aes(x=FPR, y=Recall), col="red", fill="red", size=4) + 
    theme_classic()+xlim(c(0,1))+ylim(c(0,1))+ylab("TPR")
}

#0.3 cutoff
{predicted.classes <- ifelse(probabilities > 0.3, "1", "0")
  tab<-data.frame("reaction" = test$reaction, "observed_essentiality" = test$is_essential, "predicted_essentiality" = predicted.classes)
  
  ####Performance evaluation####
  tp<-length((tab %>% filter(observed_essentiality==1) %>% filter(predicted_essentiality==1))$reaction) #true positives
  tn<- length((tab %>% filter(observed_essentiality==0) %>% filter(predicted_essentiality==0))$reaction) #true negatives
  fp<- length((tab %>% filter(observed_essentiality==0) %>% filter(predicted_essentiality==1))$reaction)#false positives
  fn<- length((tab %>% filter(observed_essentiality==1) %>% filter(predicted_essentiality==0))$reaction)#false negatives
  t<-c(tp,fp,fn,tn)
  #performance metric functions
  recall=function(t){t[1]/(t[1]+t[3])}
  specificity=function(t){t[4]/(t[2]+t[4])}
  fpr=function(t){1-specificity(t)}
  precision=function(t){t[1]/(t[1]+t[2])}
  fdr=function(t){1-precision(t)}
  #performance metric calculations
  val_recall<-recall(t)
  val_specificity<-specificity(t)
  val_fpr<-fpr(t)
  val_precision<-precision(t)
  val_fdr<-fdr(t)
  performance<-data.frame("Recall" = val_recall, "Specificity" = val_specificity,
                          "FPR" = val_fpr, "Precision" = val_precision, "FDR" = val_fdr)
  #construct plots that show performance metrics
  PRC<-ggplot()+geom_point(data=performance, aes(x=Recall, y=Precision), col="red", fill="red", size=4) + 
    theme_classic()+xlim(c(0,1))+ylim(c(0,1))
  ROC<-ggplot()+geom_point(data=performance, aes(x=FPR, y=Recall), col="red", fill="red", size=4) + 
    theme_classic()+xlim(c(0,1))+ylim(c(0,1))+ylab("TPR")
  }
}




####Validate the stepwise results manually####
##Manual model selection, forward##
{
glm(is_essential~log10_bridging_centrality, data=train, family="binomial") #AIC=955.9
glm(is_essential~log10_bridging_centrality+log10_betweeness_centrality, data=train, family="binomial") #AIC=953.5
glm(is_essential~log10_bridging_centrality+log10_betweeness_centrality+clustering_coefficient, data=train, family="binomial") #AIC=954.8
glm(is_essential~log10_bridging_centrality+log10_betweeness_centrality+log10_degree, data=train, family="binomial") #AIC=929.4
glm(is_essential~log10_bridging_centrality+log10_betweeness_centrality+log10_degree+cascadeNum, data=train, family="binomial") #AIC=920.7
    #based on the above AIC values, it looks like all centrality metrics except clustering coefficient are useful for determining reaction essentiality
}
##Now test manual model selection again, but backward##
{
glm(is_essential~log10_bridging_centrality+log10_betweeness_centrality+log10_degree+cascadeNum+clustering_coefficient, data=train, family="binomial") #AIC=921
glm(is_essential~log10_bridging_centrality+log10_betweeness_centrality+cascadeNum+clustering_coefficient, data=train, family="binomial") #AIC=934
glm(is_essential~log10_bridging_centrality+log10_betweeness_centrality+log10_degree+cascadeNum, data=train, family="binomial") #AIC=920.7
glm(is_essential~log10_bridging_centrality+log10_degree+cascadeNum, data=train, family="binomial") #AIC=934.8
glm(is_essential~log10_betweeness_centrality+log10_degree+cascadeNum, data=train, family="binomial") #AIC=921.9
glm(is_essential~log10_degree+cascadeNum, data=train, family="binomial") #AIC=963.7
glm(is_essential~log10_betweeness_centrality+cascadeNum, data=train, family="binomial") #AIC=935.5
glm(is_essential~log10_betweeness_centrality+log10_degree, data=train, family="binomial") #AIC=934.9
    #these "backward" stepwise results agree with the forward stepwise results
}

