#degree plot script
#Input: column vector of reation names, column vector of log10_degree, and column vector of whether each reaction is essential
#Output: 
  #1) plot with briding centrality on x axis, %essential reactions on left y-axis, and total number of reactions in each bin on right y-axis
  #2) Pearson's correlation between number of essential reactions and binned centrality metric
  #3) same as #2 but Spearman's method

plot_Degree<-function(df, method="linear"){
  #create a dataframe of only essential rxns
  df_essential<- df %>% filter(is_essential==1)
  #plot essential and nonessential rxns' centrality metrics
  deg<-ggplot()+geom_histogram(data=df, aes(x=log10_degree), binwidth = 0.25,
                                 col="black", fill="light gray") +
    geom_histogram(data=df_essential, aes(x=log10_degree), binwidth=0.25, col="black",
                   fill="black") +
    scale_x_continuous(name=expression('Log'[10]*'(Degree)'), breaks=c(0, 0.25, 0.50, 0.75, 1.00, 1.25,
                                                                       1.50, 1.75, 2.00, 2.25),
                       labels=c("0", "0.25","0.50", "0.75", "1.00", "1.25", "1.50", "1.75", "2.00",
                                "2.25")) +
    scale_y_continuous(name="Number of reactions", limits=c(0,300), 
                       breaks=c(0, 50, 100, 150, 200, 250, 300), 
                       labels=c("0", "50", "100", "150", "200", "250", "300"),
                       sec.axis = sec_axis(trans= ~.*35/300,name = "% Essential reactions", 
                                           breaks=seq(from=0, to=35, by=5)), expand = c(0,0)) +
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
  for(i in 1:length(ggplot_build(deg)$data[[1]]$count)){
    y<-c(y,ggplot_build(deg)$data[[1]]$count[i]) #index the brdg95 object to get the count and add to vector end
    x<-c(x, ggplot_build(deg)$data[[1]]$x[i]) #index the brdg95 object to get the x values and add to vector end
    y_essent<-c(y_essent, ggplot_build(deg)$data[[2]]$count[i])
    x_essent<-c(x_essent, ggplot_build(deg)$data[[2]]$x[i])
  }
  #convert y to a percentage by dividing by total # of reactions within each bin
  essentialRxns<-data.frame(x,y, y_essent) #combine x and y into a single dataframe; x and x_essent should be same
  essentialRxns$yPercent<-100*essentialRxns$y_essent/essentialRxns$y
  #now scale this column so that is looks correct on the 0-400 range of the right-hand y-axis
  essentialRxns$yPercentScale<-essentialRxns$yPercent * (300/35) #300 and 35 are the limits of the two y axes
  
  #take care of NaNs
  for(i in 1:length(essentialRxns$yPercentScale)){
    if(is.nan(essentialRxns$yPercentScale[i])){
      essentialRxns$yPercentScale[i]<-0
      essentialRxns$yPercent[i]<-0
    }
  }
  
  #compute correlation coefficient - looks like they ended up using pearson method in the paper (even though they don't mention their exact method) because the result is very similar
  r_degPearson<-cor(essentialRxns$x, essentialRxns$yPercent, method="pearson")
  r_degSpearman<-cor(essentialRxns$x, essentialRxns$yPercent, method="spearman")
  
  if(method=="linear"){
    #run a linear regression on the bins of the essential reactions and the essential reaction percentages
    #have to run it on the scaled results in order to plot this line
    mod<-lm(essentialRxns$yPercentScale~essentialRxns$x)
    #now add the new information to the plot
    deg<-deg + geom_point(data=essentialRxns, aes(x=x, y=yPercentScale), col="black", fill="black")+
      annotate(geom="text", x=2.25, y=230, label=paste("r =", format(round(r_degPearson, 2), nsmall=2), 
                                                     sep=" "),
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
  
  my_return<-list(deg, r_degPearson, r_degSpearman, essentialRxns)
  return(my_return)
}