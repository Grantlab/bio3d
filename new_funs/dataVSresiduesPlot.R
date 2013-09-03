## dataVSresiduesPlot.R :: create data versus residues plot
## last edited sep-03-2013
##
## sep-03-2013 (exe) created

dataVSresiduesPlot <- function(data, resIndices, secondaryStructure=NULL, MainTitle="RMSF per Residue", yTitle=expression(RMSF(ring(A))) ) {

  plotData <- data[resIndices]
  plotData <- data.frame(resIndices=resIndices, data=plotData)
  
  ##----- create the plot
  dataVSresPlot <- ggplot(plotData, aes(x=resIndices, y=data))
  ##--- if secondary structure information provided, use it
  if ( !is.null(secondaryStructure) ) {
    for (pos in 1:length(secondaryStructure$helix$start)) {
      dataVSresPlot <- dataVSresPlot + annotation_raster("pink",
                                               xmin=secondaryStructure$helix$start[pos],
                                               xmax=secondaryStructure$helix$end[pos], ymin=-Inf, ymax=Inf)
    }
    for (pos in 1:length(secondaryStructure$sheet$start)) {
      dataVSresPlot <- dataVSresPlot + annotation_raster(raster="#ffff99",
                                               xmin=secondaryStructure$sheet$start[pos],
                                               xmax=secondaryStructure$sheet$end[pos], ymin=-Inf, ymax=Inf)
    }
    for (pos in 1:length(secondaryStructure$turn$start)) {
      dataVSresPlot <- dataVSresPlot + annotation_raster("#a6cee3",
                                               xmin=secondaryStructure$turn$start[pos],
                                               xmax=secondaryStructure$turn$end[pos], ymin=-Inf, ymax=Inf)
    }
  }
  ##--- start making the plot
  dataVSresPlot <- dataVSresPlot + geom_bar(stat="identity")
  ##--- format the axes 
  dataVSresPlot <- dataVSresPlot + labs(title=MainTitle, x="Residue Index", y=yTitle)
  dataVSresPlot <- dataVSresPlot + scale_x_continuous(breaks=seq(from=0, to=length(resIndices), by=100))
  dataVSresPlot <- dataVSresPlot + theme(legend.position="none")
  dataVSresPlot <- dataVSresPlot + theme(plot.title=element_text(face="bold", size=14),
                               axis.text.x=element_text(angle=0, hjust=0.5, size=8, colour="grey20"),
                               axis.text.y=element_text(size=8, colour="grey20"),
                               axis.title.y=element_text(size=10, angle=90))

  ##----- return plot information to user
  dataVSresPlot
}
