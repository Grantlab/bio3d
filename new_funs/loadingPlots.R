## loadingPlots.R :: create loading plots
## last edited sep-03-2013
##
## sep-03-2013 (exe) created

loadingPlots <- function(bio3d.pca, secondaryStructure=NULL, resIndices, maxLoadings=3, MainTitlesYes=FALSE, compTitleLoc="top") {

  ##----- make sure the user is using the right function
  if ( maxLoadings <= 1 ) stop(
           "This function requires the loading data for at least two (2) components. Please use the dataVSresidue function to create an individual loadings plot.")

  ##-----
  compTitleLoc <- tolower(compTitleLoc)
  stopifnot ( compTitleLoc == "top" | compTitleLoc == "side" )
  #stop("The component information can either be placed on the top or to the left of each component.")

  ##----- get the component variance
  compVar.raw <- bio3d.pca$L
  compVar <- (compVar.raw/(sum(compVar.raw))) * 100
  compVar <- format( round(compVar, digits=1), digits=3)
  compVarCum <- format( cumsum(compVar), digits=3)

  ##----- get the loadings
  loadings <- abs(bio3d.pca$au[, 1:maxLoadings])

  ##----- extract needed information from the loadings
  colnames(loadings) <- paste("Component ", 1:maxLoadings,
                              " (Variance ", compVar[1:maxLoadings],
                              "%; Total Variance ", compVarCum[1:maxLoadings], "%)", sep="")
  rownames(loadings) <- resIndices
  
  ##----- transform the dataframe into the long format
  loadDF <- melt(loadings)
  ##--- rename the columns
  colnames(loadDF) <- c("index", "component", "value")

  ##----- create the plot
  loadPlot <- ggplot(loadDF, aes(x=index, y=value, group=component))
  ##--- if secondary structure information provided, use it
  if ( !is.null(secondaryStructure) ) {
    for (pos in 1:length(secondaryStructure$helix$start)) {
      loadPlot <- loadPlot + annotation_raster("pink",
                                               xmin=secondaryStructure$helix$start[pos],
                                               xmax=secondaryStructure$helix$end[pos], ymin=-Inf, ymax=Inf)
    }
    for (pos in 1:length(secondaryStructure$sheet$start)) {
      loadPlot <- loadPlot + annotation_raster(raster="#ffff99",
                                               xmin=secondaryStructure$sheet$start[pos],
                                               xmax=secondaryStructure$sheet$end[pos], ymin=-Inf, ymax=Inf)
    }
    for (pos in 1:length(secondaryStructure$turn$start)) {
      loadPlot <- loadPlot + annotation_raster("#a6cee3",
                                               xmin=secondaryStructure$turn$start[pos],
                                               xmax=secondaryStructure$turn$end[pos], ymin=-Inf, ymax=Inf)
    }
  }
  ##--- start making the plot
  loadPlot <- loadPlot + geom_bar(stat="identity")
  ##--- determine where to place the component variance information
  if ( compTitleLoc == "top" ) {
    loadPlot <- loadPlot + facet_wrap(~component, ncol=1, scale="free_y")
  } else {
    loadPlot <- loadPlot + facet_grid(component ~ ., scale="free_y")
  }
  ##--- format the axes 
  loadPlot <- loadPlot + labs(title="", x="Residue Index", y="")
  loadPlot <- loadPlot + scale_x_continuous(breaks=seq(from=0, to=length(resIndices), by=100))
  loadPlot <- loadPlot + theme(legend.position="none")
  loadPlot <- loadPlot + theme(plot.title=element_text(face="bold", size=14),
                               axis.text.x=element_text(angle=0, hjust=0.5, size=8, colour="grey20"),
                               axis.text.y=element_text(size=8, colour="grey20"),
                               axis.title.y=element_text(size=10, angle=90))

  ##----- return plot information to user
  loadPlot
}
