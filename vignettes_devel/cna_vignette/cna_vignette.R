#' # Protein Structure Networks with Bio3D
#' **Xin-Qiu Yao, Guido Scarabelli, Lars Skj\ae rven & Barry J. Grant**

#+ setup, include=FALSE
opts_chunk$set(dev='pdf')

#' ### Overview
#' The Bio3D package (version 2.1 and above) contains functionality for the creation and analysis of protein structure networks. Here we will introduce this functionality by building and analyzing networks obtained from atomic correlation data derived from normal mode analysis (NMA) and molecular dynamics (MD) simulation data (**Figure 1**). 
#' 
#' ![**Overview of a typical Bio3D network analysis protocol with key functions detailed in blue**. Normal mode and/or molecular dynamics input first undergoes **correlation analysis** (with the functions **dccm()** or **lmi()**). The output from these functions is typically a matrix of residue-by-residue cross-correlations. **Correlation network analysis** can then be performed (with the **cna()** function) to generate a correlation network with residues as nodes that are linked by weighted edges that are proportional to their degree of correlated motion. The **cna()** function also performs a **community clustering analysis** to identify _communities_ of highly correlated residues. The full residue network and coarse-grained community network can be visualized and analyzed further. This methodology can provide valuable insight not readily available from the measures of correlation alone.](figure/cna_fig1-v2.png)
#' 
#' 
#' 
#' ### Basic network generation and visualization
#' Bio3D contains extensive functionality for NMA and MD analysis. An introduction to each of these areas is detailed in separate [NMA](http://thegrantlab.org/bio3d/tutorials/normal-mode-analysis) and [MD](http://thegrantlab.org/bio3d/tutorials/trajectory-analysis) vignettes available online. In this section we will focus on correlation network analysis derived from NMA, and in subsequent sections MD.
#' 
#' The code snippet below first loads the Bio3D package, then reads an example PDB structure from the RCSB online database. We then follow the flow in Figure 1 by performing NMA, then dynamic cross-correlation analysis, and finally network analysis with the **cna()** function.
#' 
#' 
#+ example1, results="hide", warn="hide"
library(bio3d)
pdb <- read.pdb("1rv7")
modes <- nma(pdb)
cij <- dccm(modes)
net <- cna(cij, cutoff.cij=0.3)
 
#' To enable quick examination of the network analysis results we provide **cna** specific variants of the **print()**, **summary()** and **plot()** functions, as demonstrated below:
 
#+ print1
print(net)

# Summary information as a table
summary(net)$tbl


#+ figure1, results="hide", fig.cap="Full all-residue network and simplified community network"
# Plot both the ‘full’ all-residue network and simplified community network
par(mfcol = c(1, 2), mar = c(0, 0, 0, 0))
plot(net, pdb, full = TRUE)
plot(net, pdb)


 
#' 
#' Note that the above code plots both a *full* all-residue network as well as a more coarse grained community network - see **Figure 2**. In this example, and by default, the Girvan-Newman clustering method was used to map the full network into communities of highly intra-connected but loosely inter-connected nodes that were in turn used to derive the simplified network displayed in **Figure 2**. We will discus community detection in further detail below and simply highlight here that a comparison of the distribution and connectivity of communities between similar biomolecular systems (e.g. in the presence and absence of a binding partner or mutation) has the potential to highlight differences in coupled motions that may be indicative of potential allosteric pathways. This approach has previously been applied to investigate allostery in tRNA–protein complexes, kinesin motor domains, G-proteins, thrombin, and other systems (refs).
#' 
#' Another useful function for the inspection of community structure is the **view.cna()** function, which permits interactive network visualization in the VMD molecular graphics program (see **Figure 3**).
#'   


#+ viewCNA, eval=FALSE
## # View the network in VMD
## view.cna(net, pdb, launch = TRUE)


#' 
#' ![Example of **view.cna()** output](figure/output_view_cna.png)
 
#' 
#' One of the main advantages of correlation network analysis is that it often facilitates the interoperation of correlated motions that can be challenging to rationalize from the output of a dynamic correlation map alone (e.g. the output of the dccm() function). For example, the identified communities of highly intra-connected residues can be tied into the visualization of the raw correlation data (**Figure 4**) and cross referenced to the network visualizations presented in **Figures 2** and **3**:
#' 
#' 
#+ figNMAcij, results="hide", fig.cap="DCCM plot. Note that unique colors are assigned to communities and that these correspond to those annotated in all previous Figures and those used by VMD."
plot.dccm(cij, margin.segments = net$communities$membership)


#' 
#' To more fully appreciate the results of network analysis one should become somewhat familiar with the structure of the object returned from the **cna()** function. This will allow you to generate customized visualizations beyond those in **Figures 2** to **4** and also enable you to more fully examine the effect of different parameters on network structure and topology. 
#' 
#' 
#' #### The structure of cna() network objects
#' The output of the **cna()** function is a list object of class cna with both network and community structure attributes. These are summarized at the bottom of the **print(net)** results above and can be explicitly listed with the **attributes(cna)** command:
#' 
#' 

attributes(net)
 
#' 
#' Here we introduce each of these components in turn:
#' 
#' * **\$network**: The *full* protein structure network typically with a node per residue and connecting edges weighted by their corresponding correlation values (i.e. input cij values optionally filtered by the user defined 'cutoff.cij'). This filtered cij matrix is also returned in the output as **\$cij**, as listed below.
#' 
#' * **\$communities**: A community clustering object detailing the results of clustering the full $network. This is itself a list object with a number of attributes (see *attributes(net\$communities)* for details). Notable attributes include **\$communities$membership**, **$communities$modularity** and **\$communities$algorithm**. These components detail the nodes (typically residues) belonging to each community; the modularity (which represents the difference in probability of intra- and intercommunity connections for a given network division); and the clustering algorithm used respectively. 
#' * **\$cij**: The filtered correlation matrix used to build the full $network.
#' * **\$community.network**: A coarse-grained community network object with the number of nodes equal to the number of communities present in **\$communities** (i.e. the community clustering of the full network) and edges based on the inter-community coupling as determined by the input 'collapse.method' (by default this is the maximum coupling between the original nodes of the respective communities).
#' * **\$community.cij**: A cij matrix obtained by applying a "collapse.method" on \$cij. The rows and columns match the number of communities in $communities. The values are based on the cij couplings of inter-communities residue. (collapse.method available are: max, mean, trimmed and median).
#' 
#' #### Additional network annotation
#' Both the **\$network** and **$community.network** attributes are [igraph package network objects](http://igraph.sourceforge.net) and contain additional node (or vertex) and edge annotations that can be accessed with the **V()** and **E()** functions in the following way:
#'  
#' 
#+ dataLook
#' Data access
# list.graph.attributes(net$community.network)
head(V(net$network)$color)
V(net$community.network)$size
E(net$community.network)$weight
V(net$community.network)$name

#' 
#' In addition to various Bio3D functions, the full set of [igraph package](http://igraph.sourceforge.net) features can be used to further analyze these network objects.
#' 
#' #### Additional community annotation:
#' The **\$communities** attribute returned by the **cna()** function provides details of the community clustering procedure. This procedure aims to split the full network into highly correlated local substructures (referred to as communities) within which the network connections are dense but between which they are sparser.  A number of community detection algorithms can be used to solve the graph-partitioning problem. The default Girvan-Newman edge-betweenness approach is a divisive algorithm that is based on the use of the edge betweenness as a partitioning criterion. "Betweenness" is calculated by finding the shortest path(s) between a pair of vertexes and scoring each of the edges on this/these path(s) with the inverse value of the number of shortest paths. (So if there was only one path of the shortest length, each edge on it would score 1 and if there were 10 paths of that length, each edge would score 1/10.) This is done for every pair of vertexes. In this way each edge accumulates a "betweenness" score for the whole network. The network is separated into clusters by removing the edge with the highest "betweenness", then recalculating betweenness and repeating. The method is fully described in (Girvan and Newman, 2002 PNAS).
#' 
#' A useful quantity used to measure the quality of a given community partition is the **modularity** (detailed in **\$communities\$modularity** component of **cna()** output). The modularity represents the difference in probability of intra- and inter-community connections for a given network division. Modularity values fall in the range from 0 to 1, with higher values indicating higher quality of the community structure. The optimum community structures obtained for typical correlation networks fall in the 0.4 to 0.7 range.
#' 
#' 

head(net$communities$modularity)
max(net$communities$modularity)

 
#' 
#' #### To do: Demonstrate the clustering progress and choice of modularity cut-point as well as showing the as.hclust() and dendPlot() functions.
 
#' 
#' ### Network Generation From Molecular Dynamics Data
#' For this example section we apply correlation network analysis to a short molecular dynamics trajectory of Human Immunodeficiency Virus aspartic protease (HIVpr). This trajectory is included with the Bio3D package and stored in CHARMM/NAMD DCD format. To reduce file size and speed up example execution, non C-alpha atoms have been previously removed and the frames superposed. Note that typically an input trajectory would require superposition with the **fit.xyz()** function prior to correlation analysis with the **dccm()** function (see the [Trajectory Analysis vignette](http://thegrantlab.org/bio3d/html/index.html) for full details).
#' 
#' The code snippet below first sets the file paths for the example HIVpr starting structure (pdbfile) and trajectory data (dcdfile), then reads these files (producing the objects trj and pdb).
#' 
#' 
#+ infiles, results="hide"
dcdfile <- system.file("examples/hivp.dcd", package = "bio3d")
pdbfile <- system.file("examples/hivp.pdb", package = "bio3d")

 
#' 
#+ readMD, results="hide"
# Read MD data
dcd <- read.dcd(dcdfile)
pdb <- read.pdb(pdbfile)


#' 
#+ fitTraj
inds <- atom.select(pdb, resno = c(24:27, 85:90), elety = "CA")
trj <- fit.xyz(fixed = pdb$xyz, mobile = dcd, fixed.inds = inds$xyz, mobile.inds = inds$xyz)

 
#' 
#' Once we have the superposed trajectory frames we can asses the extent to which the atomic fluctuations of individual residues are correlated with one another and build a network from this data:
#' 
#' 
#+ example2MD, fig.cap="Network analysis of MD data"
cij <- dccm(trj)
net <- cna(cij)
plot(net, pdb)

 
#' 
#' Note that the **dccm()** function can be slow to run on very large trajectories. In such cases you may want to use multiple CPU cores for the calculation by setting the 'ncore' option to an appropriate value for your computer - see *help(dccm)* for details.
#' 
#' 
#+ viewDCCM, eval=FALSE
## # View the correlations in pymol
## view.dccm(cij, pdb, launch = TRUE)

#' 
#' 
#' 
#' ![alt 1](figure/HIVP_comms_colors.png)
#' ![alt 2](figure/dccm_hivpr.png)
#' 
#' 
#' #### Using a contact map filter:
#' In the original Luthy-Shulten and co-workers approach to correlation network analysis [1] edges were only included between “in contact” nodes. Contacting nodes were defined as those having any heavy atom within 5.0 Å for greater than 75% of the simulation. The weight of these contact edges was then taken from the cross-correlation data. This approach effectively ignores long range correlations between residues that are not in physical contact for the majority of the trajectory. To replicate this approach one can simply define a contact map with the **cmap()** function and provide this to the **cna()** function along with the required cross-correlation matrix, e.g.: 
#' 
#' 
#+ cmfliter
cm <- cmap(trj, dcut = 5, scut = 0, pcut = 0.75, mask.lower = FALSE)
net.cut <- cna(cij, cm = cm)

#' 
#' 
#' 
#' One could also simply multiple the correlation matrix by the contact map to zero out these long-range correlations. However, we currently prefer not to use this filtering step as we have found that this approach can remove potentially interesting correlations that might prove to be important for the types of long-range coupling we are most interested in exploring. For example, compare the plot of correlation values below to those in **Figure X**.
#' 
#' 
#+ plotDCCM, results="hide", fig.cap="DCCM plot with annotated communities"
par(mfcol = c(1, 2), mar = c(0, 0, 0, 0))
plot.dccm((cij * cm), margin.segments = net.cut$communities$membership)
plot.dccm(cij, margin.segments = net$communities$membership)

#' 
### ![DCCM plot with annotated communities](figure/unnamed-chunk-101.png) 
### ![DCCM plot with annotated communities](figure/unnamed-chunk-102.png) 
#' 
#' 
#' #### NOTE: Need to fix potential error message if mask.lower=TRUE is used in the cmap() call.
#' 
#' It is also possible to use a consensus of contact map and dynamical correlation matrices from replica simulations with the functions **cmap.filter()** and **dccm.mean()** - see their respective help pages for details.
#' 
#' 
#' ### Network visualization
#' The **layout.cna()** function can be used to determine the _geometric center_ of protein residues and communities and thus provide 3D and 2D network and community graph layouts that most closely match the protein systems 3D coords. 
#' 

layout3D <- layout.cna(net, pdb, k=3) 
layout2D <- layout.cna(net, pdb, k=2)

#' 
#' Other network layout options are available via the igraph package, see *?igraph::layout* in R for details.
#' 
#' #### Identifying nodes/communities and their composite residues
#' The **identify.cna()** function allows one to click on the nodes of a **plot.cna()** network plot to identify the residues within a give set of nodes. This can also be a useful first step to determine the placement of text labels for a plot -see *help(identify.cna)* for details.
#' 
#' 
## plotCNA,  eval=FALSE
## # Click with the right mouse to label, left button to exit.
## xy <- plot.cna(net)
## identify.cna(xy, cna = net)


#' 
#' You can also make 3D plots using the \underline{rgl} environment with the function **plot3d.cna()**. An example is shown in Figure 8.
#' 
#' 
#+ plotFalse, eval=FALSE
## plot3d.cna(net, pdb)

#' 
#' 
#' ![Interactive 3D community network plot](figure/plot3D.png) 
#' 
#' 
#' #### Color settings and node/edge labels
#' The **vmd.colors()** function is used by default to define node color. This can be changed to any input color vector as demonstrated below. Note that the **vmd.colors()** function outputs a vector of colors that is ordered as they are in the VMD molecular graphics program and changing the colors will likely break this correspondence.
#' 
#' 
bwr.colors(20)[net$communities$membership]

#' 
#' 
#' 
#' Useful options in **plot.cna()** to be aware of include the optional flags 'mark.groups' and 'mark.col', which let you to draw and color areas around groups of residues corresponding to the community clustering or any other residue partitioning of interest:
#' 

#+ plotNET, fig.cap="$network colored according to communities grouping" 
colbar.full <- vmd.colors()[net$communities$membership]
colbar.comms <- vmd.colors(max(net$communities$membership), alpha = 0.5)
plot(net, pdb, full = TRUE, mark.groups = grp.col, mark.col = colbar.comms)

#' 
#' 
#### ![$network colored according to communities grouping](figure/unnamed-chunk-12.png) 
#' 
#' 
#' #### Network manipulation
#' The **prune.cna()** function allows one to remove communities under a certain size or with less less than a certain number of edges to other communities prior to visualization or further analysis. For example, to remove communities composed of only 2 residues:
#' 
#' 
net.pruned <- prune.cna(net, size.min = 3)

#' 
#' 
#+ pruneNET, fig.cap="Original (left) and pruned (right) networks" 
par(mfcol = c(1, 2), mar = c(0, 0, 0, 0))
plot.cna(net)
plot.cna(net.pruned)

#' 
#' Applying this command to our network example, you can see how the communities composed by less than 3 residues are deleted from the graph.
#' 
#' 
#' ### Summary
#' In this vignette we have confined ourselves to introducing relatively simple correlation network generation protocols. However, the individual Bio3D functions for network analysis permit extensive customization of the network generation procedure. For example, one can chose different node representations (including all heavy atoms, residue center of mass, Calpha atoms for each residue, or collections of user defined atoms), different correlation methods (including Pearson and mutual information), different community detection methods (including Girvan-Newman betweenness, Random walk, and Greedy modularity optimization approaches), as well as customized network layouts and result visualizations. We encourage users to explore these options and provide feedback as appropriate.
#' 
#' Although the protocols detailed herein can provide compeling results in a number of cases, defining an optimal procedure for network generation is an area of active research. For example rearrangement of some residues in the communities and changes in the total number of communities may be observed among different community structures based on different independent MD simulations of the same system. These variations of community structures are partially due to the small changes in the calculated cij values that reflect statistical fluctuations among independent simulations. However, changes in mutual information or correlation values may still be relatively small and differences among independent community structures may be largely due to limitations of the modularity metric, a methodology that has been developed for binary adjacency matrices, where edges are counted without taking into account their lengths. We are currently developing a modularity methodology based on natural familial networks, where the strength of a community is defined by the number of close connections between members.
#' 
#' #### To analyze the interactions that generate the changes in the community network, we must focus on specific hydrophobic, ionic, and hydrogen-bonding interactions associated with the amino acid residues identified by the community network analysis.
#' 
#' 
#' ### References
#' 
#' [1] Sethi A, Eargle J, Black AA, Luthey-Schulten Z. Proc Natl Acad Sci U S A.2009, 106(16):6620-5.
#' 
#' #### TO BE FINISHED
