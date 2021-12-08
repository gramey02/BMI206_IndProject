#set working directory load libraries
rm(list=ls())
setwd("~/Desktop/UCSF/Biostats/Project") #change this to current working directory
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

#load data
#first file contains 2,583 reaction names in the first column and whether or not the reaction is essential in the second column
essentialRxns<-read.csv("EssentialRxns.csv", header=FALSE)
#second file contains 1,251 reaction names in the first column and cascade number for each reaction in the second column
cascadeNums<-read.csv("CascadeNums.csv", header=FALSE)
#third file contains 1,251 reactions and their centrality metrics (excluding cascade number)
centrMetrics<-read.csv("e.coli(iJ01366)-metrics-95.csv")


#######Clean and merge the data#######
#remove last column of centrMetrics (since this info is contained in EssentialRxns.csv)
drop<-"is_essential"
centrMetrics<-centrMetrics[,!(names(centrMetrics) %in% drop)]
#change(fix) column names of centrality metrics df (since these were logged already)
colnames(centrMetrics) <- c("reaction","log10_bridging_centrality","log10_betweeness_centrality", "clustering_coefficient", "log10_degree")
#change column names of cascadeNums and essentialRxns
colnames(cascadeNums)<- c("reaction", "cascadeNum")
colnames(essentialRxns)<- c("reaction", "is_essential")
#reactions c("PPK2r", "SULRi", "THRA2i", "THRAi") in cascadeNums or centrMetrics are equivalent to c("PPK2", "SULR", "THRA2", "THRA") in essentialRxns
essentialRxns$reaction[2188]<-"PPK2r"
essentialRxns$reaction[2345]<-"THRA2i"
essentialRxns$reaction[2346]<-"THRAi"
essentialRxns$reaction[2391]<-"SULRi"
#merge data by reaction name
df<-merge(centrMetrics, cascadeNums, by="reaction")
df<-merge(df, essentialRxns, by="reaction")
#take care of NaNs, non-finite values
for(i in 1:length(df$reaction)){
  if(is.finite(df$log10_bridging_centrality[i])==FALSE ||
     is.nan(df$log10_bridging_centrality[i])==TRUE){
    df$log10_bridging_centrality[i]<-0
  }
  
  if(is.finite(df$log10_betweeness_centrality[i])==FALSE ||
     is.nan(df$log10_betweeness_centrality[i])==TRUE){
    df$log10_betweeness_centrality[i]<-0
  }
}
#take care of negative values
rm(i)
for(i in 1:length(df$log10_bridging_centrality)){
  if(df$log10_bridging_centrality[i]<0){df$log10_bridging_centrality[i]<-0}
  if(df$log10_betweeness_centrality[i]<0){df$log10_betweeness_centrality[i]<-0}
  if(df$log10_degree[i]<0){df$log10_degree[i]<0}
  if(df$clustering_coefficient[i]<0){df$clustering_coefficient[i]<0}
}



####Source functions to make new plot with 5 panels, including cascade number####
#source the functions to return centrality metric plots
source("plot_BrdgCentr.R") #bridging cxentrality plot
source("plot_BetwCentr.R") #betweenness centrality plot
source('plot_CC.R') #correlation coefficient plot
source('plot_Degree.R') #degree plot
source('plot_CascNum.R') #cascade number plot



####Create 5 individual figures with linear model fits####
BrdgCentrResult<-plot_BrdgCentr(df, method="linear") #Bridging centrality results
BetwCentrResult<-plot_BetwCentr(df, method="linear") #Betweenness centrality results
CCResult<-plot_CC(df, method="linear") #Clustering coefficient results
DegreeResult<-plot_Degree(df, method="linear") #Degree results
CascNumResult<-plot_CascNum(df, method="linear") #Cascade number results



####Final linear plot arrangement####
#ggplot object is the result first element of each results list
#first row has three plots, second row has 2
plotFinal<-grid.arrange(arrangeGrob(BrdgCentrResult[[1]], BetwCentrResult[[1]], CCResult[[1]],
                                    ncol = 3),
                        arrangeGrob(DegreeResult[[1]], CascNumResult[[1]], ncol=2),
                        nrow=2)


####See if cascade number correlation is significantly different from zero####
r<-NULL
P=100000
for(p in 1:P){
  r[p]<-cor(sample((CascNumResult[[4]])$x), sample((CascNumResult[[4]])$y), method="pearson")
}
#calculate a p value associated with this
print(((sum(r>CascNumResult[[2]]) + sum(r<(-(CascNumResult[[2]]))))/100001)) #this is significant


####Save final R objects to be used in CorCompare.R and ####
save(BetwCentrResult, BrdgCentrResult, CascNumResult, CCResult, DegreeResult, file="CorResults.R")
write_csv(df, "LogisticRegr_Data.csv")
