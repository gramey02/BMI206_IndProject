#clustering coefficient plot script
#Input: column vector of reation names, column vector of log10_bridging_centrality, and column vector of whether each reaction is essential
#Output: 
  #1) plot with briding centrality on x axis, %essential reactions on left y-axis, and total number of reactions in each bin on right y-axis
  #2) Pearson's correlation between number of essential reactions and binned centrality metric
  #3) same as #2 but Spearman's method

plot_CC<-function(df, method="linear"){
  #create a dataframe of only essential rxns
  df_essential<- df %>% filter(is_essential==1)
  #plot essential and nonessential rxns' centrality metrics
  cc<-ggplot()+geom_histogram(data=df, aes(x=clustering_coefficient), binwidth = 0.1,
                                col="black", fill="light gray") +
    geom_histogram(data=df_essential, aes(x=clustering_coefficient), binwidth=0.1, col="black",
                   fill="black") +
    scale_x_continuous(name="Clustering coefficient", breaks=c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7,
                                                               0.8, 0.9, 1),
                       labels=c("0.0", "0.1","0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9",
                                "1.0")) +
    scale_y_continuous(name="Number of reactions", limits=c(0,600), 
                       breaks=c(0, 100, 200, 300, 400, 500, 600), 
                       labels=c("0", "100", "200", "300", "400", "500", "600"),
                       sec.axis = sec_axis(trans= ~.*30/600,name = "% Essential reactions", 
                                           breaks=seq(from=0, to=30, by=5)), expand = c(0,0)) +
    theme(panel.background = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border = element_rect(colour = "black", fill=NA, size=1))
  
  #get counts of the essential reactions
  #initialize vectors
  x_essent<-c()
  y_essent<-c()
  x<-c()
  y<-c()
  #anything indexed by [[1]] is all reactions, anything by [[2]] is essential rxns
  for(i in 1:length(ggplot_build(cc)$data[[1]]$count)){
    y<-c(y,ggplot_build(cc)$data[[1]]$count[i]) #index the brdg95 object to get the count and add to vector end
    x<-c(x, ggplot_build(cc)$data[[1]]$x[i]) #index the brdg95 object to get the x values and add to vector end
    y_essent<-c(y_essent, ggplot_build(cc)$data[[2]]$count[i])
    x_essent<-c(x_essent, ggplot_build(cc)$data[[2]]$x[i])
  }
  #convert y to a percentage by dividing by total # of reactions within each bin
  essentialRxns<-data.frame(x,y, y_essent) #combine x and y into a single dataframe; x and x_essent should be same
  essentialRxns$yPercent<-100*essentialRxns$y_essent/essentialRxns$y
  #now scale this column so that is looks correct on the 0-400 range of the right-hand y-axis
  essentialRxns$yPercentScale<-essentialRxns$yPercent * (600/30) #600 and 30 are the limits of the two y axes
  
  #take care of NaNs in correlation data
  for(i in 1:length(essentialRxns$yPercentScale)){
    if(is.nan(essentialRxns$yPercentScale[i])){
      essentialRxns$yPercentScale[i]<-0
      essentialRxns$yPercent[i]<-0
    }
  }
  
  #compute correlation coefficient - looks like they ended up using pearson method in the paper (even though they don't mention their exact method) because the result is very similar
  r_ccPearson<-cor(essentialRxns$x, essentialRxns$yPercent, method="pearson")
  r_ccSpearman<-cor(essentialRxns$x, essentialRxns$yPercent, method="spearman")
  
  if(method=="linear"){
    #run a linear regression on the bins of the essential reactions and the essential reaction percentages
    #have to run it on the scaled results in order to plot this line
    mod<-lm(essentialRxns$yPercentScale~essentialRxns$x)
    #now add the new information to the plot
    cc<-cc + geom_point(data=essentialRxns, aes(x=x, y=yPercentScale), col="black", fill="black")+
      annotate(geom="text", x=0.3, y=450, label=paste("r =", format(round(r_ccPearson, 2), nsmall=2), sep=" "),
             color="black", fontface="bold") +
      geom_abline(intercept=mod$coefficients[[1]], slope=mod$coefficients[[2]])
  }
  else if(method=="logistic"){
    #fit a logistic model
    mod<-glm(essentialRxns$yPercentScale~essentialRxns$x, family=binomial)
    #plot the glm on the bridging centrality plot
    cc<-cc + geom_point(data=essentialRxns, aes(x=x, y=yPercentScale), col="black", fill="black") +
      stat_smooth(method="glm", se=FALSE, method.args = list(family=binomial))
    }
  
  my_return<-list(cc, r_ccPearson, r_ccSpearman)
  return(my_return)
}