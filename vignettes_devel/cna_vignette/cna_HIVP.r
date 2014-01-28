#' # Protein Structure Networks with Bio3D 

#+ setup, include=FALSE
library(knitr)

#+ date, echo=FALSE
date()

#+ preamble, include=TRUE, eval=FALSE, echo=FALSE
#library(knitr)
#spin('cna_HIVP.r')
#system("pandoc -o cna_HIVP.pdf cna_HIVP.md")
#'
#' ## Overview
#' The Bio3D package (version 2.1 and above) contains functionality for the creation and analysis of protein structure networks. Here we will introduce this functionality by building and analyzing networks obtained from atomic correlation data derived from molecular dynamics, normal mode analysis, and large ensembles of crystal structures. We will also explore network visualization options and discuss some limitations of the current approach. 
#'
#' Key functions we will cover in this vignette include:
#'
#' - **cna()**. **C**orrelation **N**etwork **A**nalysis (CNA), create a protein structure network from an input correlation matrix and perform a community clustering analysis.
#'
#' - **summary.cna()** and **print.cna()**. Summarize and print information about a given protein structure network.
#'
#' - **plot.cna()**, **plot3d.cna()**, **layout.cna()**, **identify.cna()**, and **view.cna()**. Utilities for protein structure network visualization in 2D or 3D with structure derived community coordinates, as well as rendering network community structure into a VMD molecular graphics session. 
#'
#' - **prune.cna()** Remove small communities from a network object.
#'
#' - **dccm()** and **lmi()**, Determine the correlations of atomic displacements. 
#'
#' - **dccm.mean()**, **cmap.filter()** Calculate the consensus of dynamical correlation or contact map matrices from multiple inputs.
#'
#'
#' ### Network Generation From Normal Mode Analysis 
#' Bio3D contains extensive functionality for normal mode analysis (NMA). For an introduction to these features please refer to the [NMA vignette](http://thegrantlab.org/bio3d/tutorials/normal-mode-analysis) available online. Here we will focus on correlation network analysis derived from NMA. The code below first reads an example PDB structure, performs NMA, then dynamic cross-correlation analysis, and finally network construction.

#+ NMAbasic, results="hide"
library(bio3d)
pdb <- read.pdb("1rv7")
modes <- nma(pdb)
cij <- dccm(modes)
net <- cna(cij)

#+ NMAsmmary
# Summary information
net

#' Next we will plot the network
#+ NMAplot
par(mfcol=c(1,2), mar=c(0,0,0,0))
plot(net, pdb, full=TRUE)
plot(net, pdb)

#+ viewNMA, eval=FALSE
# View the network in VMD
view.cna(net, pdb, launch=TRUE)

#+NMAcij
plot.dccm(cij, margin.segments=net$communities$membership, sse=pdb)

#'
#' ### Network Generation From Molecular Dynamics Trajectory Data
#' For this example section we apply correlation network analysis to a short molecular dynamics trajectory of Human Immunodeficiency Virus aspartic protease (HIVpr). This trajectory is included with the Bio3D package and stored in CHARMM/NAMD DCD format. Note that all solvent and non C-alpha protein atoms have been excluded to reduce overall file size.

#' The code snippet below loads the Bio3D package, sets the file paths for the example HIVpr starting structure (pdbfile) and trajectory data (dcdfile), and finally reads these files (producing the objects dcd and pdb).

#+ readtrj, results="hide"
library(bio3d)
dcdfile <- system.file("examples/hivp.dcd", package="bio3d")
pdbfile <- system.file("examples/hivp.pdb", package="bio3d")
dcd <- read.dcd(dcdfile)
pdb <- read.pdb(pdbfile)

#' For further details of trajectory IO and subsequent analysis please refer to the Trajectory Analysis vignette available [online](http://thegrantlab.org/bio3d/html/index.html).
#'
#' The steps for correlation network analysis typically entail 
#'
#' ### Cross-Correlation Analysis
#' The extent to which the atomic fluctuations/displacements of a system are correlated with one another can be assessed by examining the magnitude of all pairwise cross-correlation coefficients. The Bio3D `dccm()` function returns a matrix of all atom-wise cross-correlations whose elements may be displayed in a graphical representation frequently termed a dynamical cross-correlation map, or DCCM. Prior to this calculation we will superpose the trajectory frames based on the the C-alpha atoms of residues 24 to 27 and 85 to 90 in both chains. 

#+ fitting, results="hide"
inds <- atom.select(pdb, resno=c(24:27,85:90), elety="CA")
xyz <- fit.xyz(fixed=pdb$xyz, mobile=dcd,
               fixed.inds=inds$xyz,
               mobile.inds=inds$xyz)

#'
ca.inds <- atom.select(pdb, "calpha")
cij <- dccm(xyz[,ca.inds$xyz])

#+ plotDCCM, echo=TRUE, fig.cap="Correlation analysis of HIV protease."
plot(cij)

#'
#' ### Note:
#' The **dccm()** function can be slow to run on larger trajectories. In such cases you may want to use multiple CPU cores for the calculation by setting the 'ncore' option, see help(dccm).
#'
#' Figure 1 indicates that the two chains of HIVpr have somewhat similar correlation patterns (residues 1 to 99 for chain A, and 100 to 198 for chain B). However, further dissecting the pattern of off-diagonal correlations may require considerable (often manual) inspection. One can annotate the plot with structural features such as secondary structure element locations etc. (see help(dccm) for details) and view the a rendering of the pair-wise correlations on a representative structure (Figure 2.) can be informative. However, we believe that network analysis can be particularly helpful in this regard.  

#' A 3D visualization of these correlations can be provided through the function `view.dccm()`

#+ viewdccm, eval=FALSE
# View the correlations in pymol
view.dccm(cij, pdb, launch=TRUE)


#' \includegraphics{figure/dccm_hivpr.png}

#'
#' ### Note:
#' You can also use linear mutual information (see help(lmi)) and or take a consensus of cij matrices obtained from replica simulations with the function **dccm.mean()**. A mean of the couplings whose absolute values are higher than a selected cutoff in a specified number of simulations is calculated. See help(dccm.mean) for a complete explanation.
#'
#' ## Correlation network construction and analysis
#'
#' Now we will use the cross-correlation data to partition our protein structure into dynamically consistent communities of residues (dynamic subdomains) by clustering an undirected network of residue nodes that are edge weighted by their correlations.
#'
#' ### Cna function input
#'
#' Here we will use the **cna()** function. This function takes as input the cij matrix and creates a network where each residue is a node, and two nodes are connected by the -log(abs(cij)) values. Filters on the cij values can be applied. The option cij.cutoff lets you exclude all the cij whose absolute value is lower than the specified cutoff (by default abs(cij)<0.4 are excluded).
#' Another filter is the multiplication of the cij matrix by a contact matrix. This option is useful to build a network using only the correlations of residues that are in contact during the simulation.
#' To have more information on the filters, please see help(cna).
#'
#' The **cna()** function will also cluster the network into communities and produce a clustered network for further analysis.
#' Different methods are available to identify communities, the default one uses the betweenness measure with the girvan-newman betweenness clustering algorithm to determine the communities (cluster.method="btwn").
#' Other methods available include random walk partitioning (cluster.method="walk") and fast-greedy modularity optimization (cluster.method="greed"). For more information, please see help(cna).
#'

#+ Network1
# Build and betweeness cluster a correlation network graph 
net <- cna(cij)
net

#' ### cna function output
#'
#' The output of the **cna()** function consists of an object containing the following attributes:
#'
#' - **$network**: The **full** protein structure network with a node per residue and connecting edges weighted by their corresponding filtered correlation values (i.e. input cij values filtered by 'cutoff.cij'). This filtered cij matrix is also returned as **$cij**.
#'
#' - **$communities**: A community clustering object detailing the results of clustering the $network. Notable attributes include **$communities$membership** and **$communities$algorithm** that detail **$network** nodes (i.e. residues) belonging to each community and the clustering algorithm used respectively.  
#'
#' - **$cij**: The filtered cij matrix used to build the network.
#'
#' - **$community.network**: A coarse-grained community network object with the number of nodes equal to the number of communities present in **$communities** (i.e. the clustering of the full network) and edges based on the inter-community coupling as determined by the input 'collapse.method' (by default this is the maximum coupling between the original nodes of the respective communities).
#'
#' - **$community.cij**: A cij matrix obtained by applying a "collapse.method" on $cij. The rows and columns match the number of communities in $communities. The values are based on the cij couplings of inter-communities residue. (collapse.method available are: max, mean, trimmed and median).
#'
#' Data access
head( V(net$network)$color )
V(net$community.network)$size
E(net$community.network)$weight
V(net$community.network)$name

#' 
#' ## Plot the network
#'
#' The function **plot.cna()** plots the obtained network.

# Plot a simple network graph highlighting community clusters
#+ coordinates
x <- plot(net)

#' To pick out the the residue corespondance of each community you can use the **identify.cna function()**
# identify.cna(x, net) ## to add labels to plot


#' Different plot layout schemes are available including those based on PDB coordinates
par(mfcol=c(1,2), mar=c(0,0,0,0))
x <- layout.cna(net, pdb)
plot(net, layout=x)

y <- layout.cna(net, pdb, full=TRUE)
plot(net, layout=y, full=TRUE, vertex.size=10, weights=1)

#' You could multiply y by -1 to enforce a similar layout in the two plots.
#plot3d.cna(net, pdb, layout=layout.cna(net, pdb, k=3))
# plot3d.cna(net, pdb)

#+ viewcna, eval=FALSE
# View the correlations in pymol
view.cna(net, pdb, launch=TRUE)

#' In addition, a cij matrix can be plotted with the community colors calling the function **plot.dccm()** and declaring the margin.segments option equal to the community memberships obtained from the network analysis.
#'
# This plot will be anotated with the same colors as the network plots and VMD output.
#+ fig.cap="cij couplings. Communities are shown on the bottom and left axes, secondary structure elements on the top and right ones."
plot.dccm(cij, margin.segments=net$communities$membership)
#'
#' ### Sidenote:
#'
#' To make the two plots with the same layout, we used the option 'layout=layout.pdb()', see below for further explanations on this.
#' Anyway, you can also make a plot without that and using a random layout with 'plot.cna(net)'.
#'
#' There are different levels to represent a network. One possibility is to plot each residue as a node, as shown in the left panel of Figure 2.
#' The network looks messy to interpret, but still with a careful looking it is possible to identify some features. In our example you can identify the two protein lobes colored in purple and dark green, plus the region where the two domain closely interact (which is clustered into many communities). The residues are colored by community membership (i.e. residues with the same color belong to the same community).
#' The plot on the right is the community network, which has a node for each community. Each of these nodes is scaled by the number of residues within that community and shows a clearer representation of the community clustering.
#'
#'
#' \includegraphics{figure/HIVP_comms_colors.png}
#'
#'
#' ### Sidenote: On network graph visualization
#' Now is probably a good time to take an interlude to talk about network graph visualization.
#' It's useful to think of there being two kinds of network visualizations: _pretty plots_ that encode 'soft' information and _ugly plots_ that encode 'hard' information. Most of what you see in journals etc falls under the first category.
#'
#' These 'soft' plots try to display a network so you can see broad structure and patterns in a totally non-rigourous way. Most of them are based on a sort of a physical metaphor: imagine the nodes of the graphs are little balls that are contantly repelling on another like two mis-aligned magnets, but somebody went and attached a bunch of them together with springs. You throw them on the table and see how they settle. This often gives really nice results, with vertices well spread out, and with highly connected vertices tending to be closer together. However, it is based on somewhat heuristic layout methods and it can be easy to read too much into a particular layout if you do not appreciate these heuristics. Plus each time you run it it may have a different orientation/layout.
#'
#' However, as briefly mentioned previously, we have a much better way to layout our protein network plots by simply taking advantage of the Cartesian coordinates of the atoms representing a particular community node. For this purpose we can use the **layout.pdb()** function (see below).
#'
#'
#' ## Schematic Network Table
#'
#' Besides plotting the network, it is also useful to summarize its structure more schematically. **summary.cna()** and **print.cna** functions allow us to examine the community structure. They print on the screen a schematic table with the network elements.
#' **summary.cna()** function lets you to save the table in a R object.
#' In the output, the first column represents the community IDs, the second the community sizes (in terms of number of residues belonging to that community) and the last columns the residue members.

# summary of HIVP protein network
x <- summary.cna(net)

attributes(x)
x$tbl

#' In our example, there are 18 communities. The two biggest ones (ids #6 and #12) correspond to the two HIVP lobes forming the active site of the protein. There are also other smaller communities in the two domains contact region. Overall, the simmetry of the protein structure is well represented by the community clustering.
#'
#' ### Note on attributes
#'
#' You can store whatever annotation you want along with our network graph for example:
net$community.network$title="Dynamic Network of HIVP"
list.graph.attributes(net$community.network)
net$community.network$title
#'
#' Now let's move on other functions useful to plot the network data generated.
#'
#' ## Other plot functions
#'
#' ### Color settings
#' The function **vmd.colors()** can be used to change the default color of the nodes. This function takes in input the number of different colors you want to generate, and gives an output a vector of corresponding VMD colors (ordered as they are in VMD molecular graphics program).
#'
colbar.nodes <- vmd.colors()[net$communities$membership]
colbar.comms <- vmd.colors(max(net$communities$membership), alpha=0.5)

# Creating a list object where the residues are grouped by their community membership
grp.col <- list()
for(i in 1:max(net$communities$membership)){
  grp.tmp <- which(net$communities$membership == i)
  grp.col[[i]] <- grp.tmp
}

#'
#' The flags 'mark.groups' and 'mark.col' let you to draw and color areas around groups of residues corresponding to the community clustering
#+ fig.cap="$network colored according to communities grouping"
plot(net$network, mark.groups=grp.col, mark.col=colbar.comms, vertex.size=10, layout=coords.network)
#'
#'
#' ### Identifying nodes/communities and their composite residues
#' The **identify.cna()** function allows one to click on the nodes of a plot.cna() network plot to identify the residues within a give set of nodes. This can also be a useful first step to determine the placment of text labes to a plot. You should try it in an active R session, see help(identify.cna) for details.
#'
#'
#+ IdentifyPlot, eval=FALSE
# Click with the right mouse to lable, left button to exit.
#+ fig.cap="btwn network plot"
xy<-plot.cna(net)
identify.cna(xy, cna=net)

#' ### Plot layout
#' We can use **layout.pdb()** function to detemine the _geometric center_ of protein residues and communities and also dedermine 3D and 2D network and community graph layouts that most closely match the protein systems 3D coords. See figure 2 for **layout.pdb()** example.
#'

#+ layout.coords, echo=FALSE
layout3D <- layout.pdb(pdb, net) ## EDIT these!!
layout2D <- layout.pdb(pdb, net, k=2)

#' ### Network as 3D plots
#' You can make 3D plots using the \underline{rgl} environment with the function **plot3d.cna()**. An example is shown in Figure 6.
#'
#+ plot3D, eval=TRUE, custom.plot=TRUE
knit_hooks$set(custom.plot=hook_plot_custom)
plot3d.cna(net, pdb)
rgl.viewpoint(260,-90, zoom=0.9)
rgl.snapshot( paste(fig_path(),"png", sep=".") )
rgl.quit()


#' ### Note:
#'
#' You can use the function **view.cna()** to draw spheres in VMD corresponding to the identified communities and color the protein "by chain" in order to match community colors. We are not running that command in this example, but it is reported below with a picture of its output.
#+ eval=FALSE
## view.cna(net, layout=layout.pdb(pdb,net), launch=TRUE)
#'
#+ fig.cap="example of **view.cna()**"
#' \includegraphics{figure/output_view_cna.png}
#'
#'
#' ## Collapse and/or prune nodes.
#' The **prune.cna()** function allows us to remove communities under a certain size or with less less than a certain number of edges to other communities pior to visualization of further analysis. For example, communities composed of only 2 residues.
#'
#+ PruneNetwork
#+ fig.cap="Original (left) and pruned (right) networks"
net.pruned <- prune.cna(net, size.min=3)
par(mfcol=c(1,2), mar=c(0,0,0,0))
plot.cna(net)
plot.cna(net.pruned)

#' Applying this command to our network example, you can see how the communities composed by less than 3 residues are deleted from the graph.
#'
#'
#' ## Contact map filtering
#' We can also filter the dccm by a given contact map similar to what was done in the work from the Luthy-Shulten group [1]. This filtering removes a lot of signal from the correlation matrix, as you can see from the plot of the contact map matrix.

#+ cmap, cache=TRUE
# Contact map based on calpha
cm <- cmap(trj[,ca.inds$xyz], pcut=0.75, scut=0, dcut=5, mask.lower=FALSE, ncore=8)
#+ fig.cap="contact map matrix"
plot.dccm(cm, main="cm")

#' ### Note:
#'
#' As for the cij matrix, it is possible to use a consensus of contact map matrices from replica simulations with the function **cmap.filter()**. See **?cmap.filter** for further details.
#'
net.cmap <- cna(cij, cluster.method="btwn", cutoff.cij=0.4, cm)
layout2D.cmap <- layout.pdb(pdb, net.cmap, k=2)

#+ fig.cap="network without (left) and with (right) cm filter"
par(mfcol=c(1,2), mar=c(0,0,0,0))
plot.cna(net, layout=layout2D)
plot.cna(net.cmap, layout=layout2D.cmap)

#' The cm filter removes a lot of cij couplings, as consequence it might exclude useful information. In our example, the application of the cm filter affects the resulting network breaking the protein into small communities.
#'
#' ### Network Generation From Normal Mode Analysis
#' For this example section we apply normal mode analysis to the Human Immunodeficiency Virus aspartic protease (HIVpr). This trajectory is included with the Bio3D package and stored in CHARMM/NAMD DCD format. Note that all solvent and non C-alpha protein atoms have been excluded to reduce overall file size.
modes <- nma(pdb)

## Calculate correlation matrix
cm <- dccm.nma(modes)

#+ nmaDCCMplot
## Plot correlation map
plot(cm)
#plot(cm, contour = FALSE, col.regions = bwr.colors(20), at = seq(-1, 1, 0.1))

## Network analysis
net2 <- cna(cm)

#+ nmaNETplot
par(mfcol=c(1,2), mar=c(0,0,0,0))
plot(net2, layout=layout.cna(net2, pdb))
plot(net2, layout=layout.cna(net2, pdb, full=T), full=T, vertex.size=12)


#+ viewCNA, eval=FALSE
# View the community network in VMD
view.cna(net2, pdb, launch=TRUE)

#+ fig.cap="example of NMA derived community structure network produced by **view.cna()**"
#' \includegraphics{figure/nma_network.png}

#'
#' ## Conclusions
#'
#' Now you have played with the main functions of feature-cna branch and should have an idea of the different output that can help you to analyse the information in a correlation matrix.
#' But we are aware of the limitations intrisically present in the cna analysis, such as the parameter chosen to build the network, the clustering of the nodes, etc. For these aspects, we trust the user to discern what in the output can be biologically meaningful and what it is not...
#'
#'
#' ## References
#'
#' [1] Sethi A, Eargle J, Black AA, Luthey-Schulten Z. Proc Natl Acad Sci U S A.2009, 106(16):6620-5.
