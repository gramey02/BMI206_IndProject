#script to compare correlations (between reaction essentiality and centrality measures) calculated with either Pearson's or Spearman's method

rm(list=ls())
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
library(tidyr)

#load correlations for each centrality metric--Pearson's correlation is in [[2]] of each list, and Spearman's is in [[3]]
load('CorResults.R')

#for each centrality metric, compute the absolute value of the difference between pearson's and spearman's correlation coefficients
dif_betwCentr<-(BetwCentrResult[[2]] - BetwCentrResult[[3]])
dif_brdgCentr<-(BrdgCentrResult[[2]] - BrdgCentrResult[[3]])
dif_cc<-(CCResult[[2]] - CCResult[[3]])
dif_deg<-(DegreeResult[[2]] - DegreeResult[[3]])
dif_cn<-(CascNumResult[[2]] - CascNumResult[[3]])

#parse out the dataframes for each centrality metric
df_brC<-BrdgCentrResult[[4]]
df_betC<-BetwCentrResult[[4]]
df_cn<-CCResult[[4]]
df_cn<-DegreeResult[[4]]
df_cn<-CascNumResult[[4]]

#construct a null distribution for each centrality metric by permuting the values and calculating correlation differences

#first for bridging centrality
{
  P=10000 #number of times we'll run the permutation
  perm=NULL
  for(p in 1:P){
    #randomly rearrange values in each column of centrality metrics' dataframes
    p_cor<-cor(sample(df_brC$x), sample(df_brC$yPercent), method="pearson")
    s_cor<-cor(sample(df_brC$x), sample(df_brC$yPercent), method="spearman")
    
    perm[p]<-(p_cor - s_cor)
  }
  
  brCplot<-ggplot()+geom_histogram(data=as.data.frame(perm), aes(perm), bins=100, color="black", fill="black")+
    geom_vline(xintercept=mean(perm), color="red")+
    geom_vline(xintercept=dif_brdgCentr, color="blue")+
    geom_vline(xintercept=mean(perm)+(1.96*sd(perm)), linetype="dashed")+
    geom_vline(xintercept=mean(perm)-(1.96*sd(perm)), linetype="dashed")+
    theme_classic()+
    scale_y_continuous(expand=c(0,0), limits=c(0, 350))+
    xlab("|Difference in correlations|")+
    ylab("Frequency")+
    labs(title="Bridging Centrality")+
    theme(plot.title = element_text(hjust = 0.5))+
    annotate(geom="text", x=-0.17, y=325, label="-0.03", color="blue", fontface="bold")
  
  
  #calculate a p value associated with this
  print(((sum(perm<dif_brdgCentr) + sum(perm>abs(dif_brdgCentr)))/10001))
}

#next for betweenness centrality
{
  P=10000 #number of times we'll run the permutation
  perm=NULL
  for(p in 1:P){
    #randomly rearrange values in each column of centrality metrics' dataframes
    p_cor<-cor(sample(df_betC$x), sample(df_betC$yPercent), method="pearson")
    s_cor<-cor(sample(df_betC$x), sample(df_betC$yPercent), method="spearman")
    
    perm[p]<-(p_cor - s_cor)
  }
  
  betCplot<-ggplot()+geom_histogram(data=as.data.frame(perm), aes(perm), bins=100, color="black", fill="black")+
    geom_vline(xintercept=mean(perm), color="red")+
    geom_vline(xintercept=dif_betwCentr, color="blue")+
    geom_vline(xintercept=mean(perm)+(1.96*sd(perm)), linetype="dashed")+
    geom_vline(xintercept=mean(perm)-(1.96*sd(perm)), linetype="dashed")+
    theme_classic()+
    scale_y_continuous(expand=c(0,0), limits=c(0, 350))+
    xlab("|Difference in correlations|")+
    ylab("Frequency")+
    labs(title="Betweenness Centrality")+
    theme(plot.title = element_text(hjust = 0.5))+
    annotate(geom="text", x=0.17, y=325, label="0.02", color="blue", fontface="bold")
  
  #calculate a p value associated with this
  print(((sum(perm>dif_betwCentr) + sum(perm<(-dif_betwCentr)))/10001))
  }

#next for correlation coefficient
{
  P=10000 #number of times we'll run the permutation
  perm=NULL
  for(p in 1:P){
    #randomly rearrange values in each column of centrality metrics' dataframes
    p_cor<-cor(sample(df_cn$x), sample(df_cn$yPercent), method="pearson")
    s_cor<-cor(sample(df_cn$x), sample(df_cn$yPercent), method="spearman")
    
    perm[p]<-(p_cor - s_cor)
  }
  
  CCplot<-ggplot()+geom_histogram(data=as.data.frame(perm), aes(perm), bins=100, color="black", fill="black")+
    geom_vline(xintercept=mean(perm), color="red")+
    geom_vline(xintercept=dif_cc, color="blue")+
    geom_vline(xintercept=mean(perm)+(1.96*sd(perm)), linetype="dashed")+
    geom_vline(xintercept=mean(perm)-(1.96*sd(perm)), linetype="dashed")+
    theme_classic()+
    scale_y_continuous(expand=c(0,0), limits=c(0, 350))+
    xlab("|Difference in correlations|")+
    ylab("Frequency")+
    labs(title="Clustering Coefficient")+
    theme(plot.title = element_text(hjust = 0.5))+
    annotate(geom="text", x=0.17, y=325, label="0.01", color="blue", fontface="bold")
  
  
  #calculate a p value associated with this
  print(((sum(perm>dif_cc) + sum(perm<(-dif_cc)))/10001))
}

#next for degree
{
  P=10000 #number of times we'll run the permutation
  perm=NULL
  for(p in 1:P){
    #randomly rearrange values in each column of centrality metrics' dataframes
    p_cor<-cor(sample(df_cn$x), sample(df_cn$yPercent), method="pearson")
    s_cor<-cor(sample(df_cn$x), sample(df_cn$yPercent), method="spearman")
    
    perm[p]<-(p_cor - s_cor)
  }
  
  DEGplot<-ggplot()+geom_histogram(data=as.data.frame(perm), aes(perm), bins=100, color="black", fill="black")+
    geom_vline(xintercept=mean(perm), color="red")+
    geom_vline(xintercept=dif_deg, color="blue")+
    geom_vline(xintercept=mean(perm)+(1.96*sd(perm)), linetype="dashed")+
    geom_vline(xintercept=mean(perm)-(1.96*sd(perm)), linetype="dashed")+
    theme_classic()+
    scale_y_continuous(expand=c(0,0), limits=c(0, 350))+
    xlab("|Difference in correlations|")+
    ylab("Frequency")+
    labs(title="Degree")+
    theme(plot.title = element_text(hjust = 0.5))+
    annotate(geom="text", x=0.25, y=325, label="0.10", color="blue", fontface="bold")
  
  
  #calculate a p value associated with this
  print(((sum(perm>dif_deg) + sum(perm<(-dif_deg)))/10001))
}

#last for cascade number
{
  P=10000 #number of times we'll run the permutation
  perm=NULL
  for(p in 1:P){
    #randomly rearrange values in each column of centrality metrics' dataframes
    p_cor<-cor(sample(df_cn$x), sample(df_cn$yPercent), method="pearson")
    s_cor<-cor(sample(df_cn$x), sample(df_cn$yPercent), method="spearman")
    
    perm[p]<-(p_cor - s_cor)
  }
  
  CNplot<-ggplot()+geom_histogram(data=as.data.frame(perm), aes(perm), bins=100, color="black", fill="black")+
    geom_vline(xintercept=mean(perm), color="red")+
    geom_vline(xintercept=dif_deg, color="blue")+
    geom_vline(xintercept=mean(perm)+(1.96*sd(perm)), linetype="dashed")+
    geom_vline(xintercept=mean(perm)-(1.96*sd(perm)), linetype="dashed")+
    theme_classic()+
    scale_y_continuous(expand=c(0,0), limits=c(0, 350))+
    xlab("|Difference in correlations|")+
    ylab("Frequency")+
    labs(title="Cascade Number")+
    theme(plot.title = element_text(hjust = 0.5))+
    annotate(geom="text", x=0.28, y=325, label="0.03", color="blue", fontface="bold")
  
  
  #calculate a p value associated with this
  print(((sum(perm>dif_cn) + sum(perm<(-dif_cn)))/10001))
}

finalPlot<-grid.arrange(arrangeGrob(brCplot, betCplot, CCplot, ncol = 3),
                        arrangeGrob(DEGplot, CNplot, ncol=2),
                        nrow=2)
