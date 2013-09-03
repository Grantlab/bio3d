## scorePlot.R :: create principal component score plots
## last edited sep-02-2013
##
## sep-02-2013 (exe) created

scorePlot <- function (scores, comp1=1, comp2=2, MainTitleYes=FALSE, classification=NULL) {

  ##----- create the scores dataframe
  if ( is.null(classification) ) {
    scoresDF <- data.frame(compX=scores[, comp1], compY=scores[, comp2])
  } else {
    scoresDF <- data.frame(Classification=as.integer(classification), compX=scores[, comp1], compY=scores[, comp2])
  }

  ##----- determine the color scheme
  dotColors <- c("red", "purple", "yellow", "blue", "orange", "green", "black", "brown")

  ##----- create the titles
  if ( MainTitleYes == TRUE ) {
    MainTitle <- paste("Component", comp1, "versus Component", comp2, sep=" ")      
  } else {
    MainTitle <- NULL
  }
  xTitle <- paste("Component", comp1, sep=" ")
  yTitle <- paste("Component", comp2, sep=" ")

  ##----- make the plot
  if ( is.null(classification) ) {
    scorePlot <- ggplot(scoresDF, aes(x=compX, y=compY))
  } else {
    scorePlot <- ggplot(scoresDF, aes(x=compX, y=compY, colour=Classification))
  }
  scorePlot <- scorePlot + geom_point(size=3)
  scorePlot <- scorePlot + scale_colour_identity()
  scorePlot <- scorePlot + labs(title=MainTitle,
                                x=xTitle, y=yTitle)
  scorePlot <- scorePlot + theme(plot.title = element_text(face="bold", size=14),
                                 axis.text.x=element_text(size=10, colour="grey20"),
                                 axis.text.y=element_text(size=10, colour="grey20"),
                                 axis.title.x=element_text(size=12),
                                 axis.title.y=element_text(size=12, angle=90))
  scorePlot <- scorePlot + theme(legend.title=element_text(face="bold", size=11),
                                 legend.title.align=1,            ## center the legend title
                                 legend.text=element_text(size=10),
                                 legend.key.size=unit(0.6, "cm"))
  scorePlot <- scorePlot + theme(legend.position="none")

  ##----- return plot information to user
  scorePlot  
}
