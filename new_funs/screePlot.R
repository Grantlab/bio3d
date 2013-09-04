## screePlot.R :: create principal component scree plots
## last edited sep-02-2013
##
## sep-02-2013 (exe) created

screePlot <- function(compVar) {

  ##----- calculate the cumulative sum and variance index
  compVar <- round(compVar[1:20], digits=1)

  ##----- set the component indices and updated component variance vector
  compIdx <- 1:length(compVar)
  compVarSum <- cumsum(compVar)
  ##--- create the label vector
  labelValues <- format(compVarSum, digits=3)
  labelIdx <- compIdx[c(1:5,10,13,16,19)]

  ##----- determine the range of the y-axis
  yMax <- max(compVar) + 1
  yMinMax <- c(0, yMax)
  
  ##----- create the dataframe of data to be plotted
  screeDF <- data.frame(compIdx, compVar, label=labelValues)

  ##----- create the titles
  MainTitle <- paste("Component Variance")
  xTitle <- paste("Component")
  yTitle <- paste("Percent Variance (%)")

  ##----- make the plot
  screePlot <- ggplot(screeDF, aes(x=compIdx, y=compVar, label=label) )
  screePlot <- screePlot + geom_line(colour="navy")
  screePlot <- screePlot + geom_point(size=3, colour="navy")
  screePlot <- screePlot + ylim(yMinMax)
  screePlot <- screePlot + geom_text(data=screeDF[labelIdx, ],
                                     hjust=-0.20, vjust=-0.4, angle=15, colour="navy", size=4.5)
  screePlot <- screePlot + labs(title=MainTitle, x=xTitle, y=yTitle)
  screePlot <- screePlot + theme(plot.title = element_text(face="bold", size=14),
                                 axis.text.x=element_text(size=10, colour="grey20"),
                                 axis.text.y=element_text(size=10, colour="grey20"),
                                 axis.title.x=element_text(size=12),
                                 axis.title.y=element_text(size=12, angle=90))  

  ##----- return plot information to user
  screePlot
}

