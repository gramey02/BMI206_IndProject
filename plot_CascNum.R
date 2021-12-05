#cascade number plot script
#Input: column vector of reation names, column vector of cascade number, and column vector of whether each reaction is essential
#Output: 
  #1) plot with cascade number on x axis, %essential reactions on left y-axis, and total number of reactions in each bin on right y-axis
  #2) Pearson's correlation between number of essential reactions and binned centrality metric
  #3) same as #2 but Spearman's method

plot_CascNum<-function(df, method="linear"){
  #create a dataframe of only essential rxns
  df_essential<- df %>% filter(is_essential==1)
  #plot essential and nonessential rxns' centrality metrics
  cn<-ggplot()+geom_histogram(data=df, aes(x=cascadeNum), binwidth = 1,
                                 col="black", fill="light gray") +
    geom_histogram(data=df_essential, aes(x=cascadeNum), binwidth=1, col="black",
                   fill="black") +
    scale_x_continuous(name=expression('Cascade number'), breaks=c(0, 1, 2, 3, 4, 5,
                                                                       6, 7),
                       labels=c("0", "1","2", "3", "4", "5", "6", "7")) +
    scale_y_continuous(name="Number of reactions", limits=c(0,1000), 
                       breaks=c(0, 200, 400, 600, 800, 1000), 
                       labels=c("0", "200", "400", "600", "800", "1000"),
                       sec.axis = sec_axis(trans= ~.*100/1000,name = "% Essential reactions", 
                                           breaks=seq(from=0, to=100, by=20)), expand = c(0,0)) +
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
  for(i in 1:length(ggplot_build(cn)$data[[1]]$count)){
    y<-c(y,ggplot_build(cn)$data[[1]]$count[i]) #index the brdg95 object to get the count and add to vector end
    x<-c(x, ggplot_build(cn)$data[[1]]$x[i]) #index the brdg95 object to get the x values and add to vector end
    y_essent<-c(y_essent, ggplot_build(cn)$data[[2]]$count[i])
    x_essent<-c(x_essent, ggplot_build(cn)$data[[2]]$x[i])
  }
  #convert y to a percentage by dividing by total # of reactions within each bin
  essentialRxns<-data.frame(x,y, y_essent) #combine x and y into a single dataframe; x and x_essent should be same
  essentialRxns$yPercent<-100*essentialRxns$y_essent/essentialRxns$y
  #now scale this column so that is looks correct on the 0-400 range of the right-hand y-axis
  essentialRxns$yPercentScale<-essentialRxns$yPercent * (1000/100) #1000 and 35 are the limits of the two y axes
  
  #take care of NaNs
  for(i in 1:length(essentialRxns$yPercentScale)){
    if(is.nan(essentialRxns$yPercentScale[i])){
      essentialRxns$yPercentScale[i]<-0
      essentialRxns$yPercent[i]<-0
    }
  }
  
  #compute correlation coefficient - looks like they ended up using pearson method in the paper (even though they don't mention their exact method) because the result is very similar
  r_cnPearson<-cor(essentialRxns$x, essentialRxns$yPercent, method="pearson")
  r_cnSpearman<-cor(essentialRxns$x, essentialRxns$yPercent, method="spearman")
  
  if(method=="linear"){
    #run a linear regression on the bins of the essential reactions and the essential reaction percentages
    #have to run it on the scaled results in order to plot this line
    mod<-lm(essentialRxns$yPercentScale~essentialRxns$x)
    #now add the new information to the plot
    cn<-cn + geom_point(data=essentialRxns, aes(x=x, y=yPercentScale), col="black", fill="black")+
      annotate(geom="text", x=2.25, y=230, label=paste("r =", format(round(r_cnPearson, 2), nsmall=2), 
                                                     sep=" "),
             color="black", fontface="bold") +
      geom_abline(intercept=mod$coefficients[[1]], slope=mod$coefficients[[2]])
  }
  else if(method=="logistic"){
    #fit a logistic model
    mod<-glm(essentialRxns$yPercentScale~essentialRxns$x, family=binomial)
    #plot the glm on the bridging centrality plot
    cn<-cn + geom_point(data=essentialRxns, aes(x=x, y=yPercentScale), col="black", fill="black") +
      stat_smooth(method="glm", se=FALSE, method.args = list(family=binomial))
    }
  
  my_return<-list(cn, r_cnPearson, r_cnSpearman)
  return(my_return)
}