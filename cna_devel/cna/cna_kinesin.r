#' # Correlation Network Analysis with Bio3D 

#+ setup, include=FALSE
library(knitr)
#opts_chunk$set(fig.path='figure/cna-', cache.path='cache/cna-', dev='pdf', fig.align='center', fig.width=5, fig.height=5, cache=TRUE, par=TRUE)


#+ date, echo=FALSE
date()

#' ## Key Function List
#' - cna()  _SORT OUT RAW AND CLUSTERED NOTATION_
#' - summary.cna(), print.cna()
#' - plot.cna(), plot3d.cna(), prune.cna()  _COMBINE THESE_
#' - view.cna(), layout.pdb()
#' - dccm.mean(), cmap.filter() _NEED TO REVIST THESE_
#' - Others: plot.dccm2(), vmd.colors(), identify.cna()   _ADD annotation to dccm_
#'
#' ## To Do List:
#' - Incoporate kinesin multi trajectory data into vignette.
#' - Incoporate nma() results into vignette.
#' - Incoporate contract.vertices() into vignette.
#'          y=simplify(contract.vertices(x$raw.network, mapping = x$membership))
#'
#' - Write/improve a prune.cna() function and once working
#'      incoporate into plot.cna() and plot3d.cna() functions
#' - Write/improve view.cna() for both vmd and pymol
#' - Improve identify.cna() so 2nd argumet could be cna object.
#'
#' ## Check List:
#' - Test with larger trajectory and structure files and check for bugs and useful results
#'
#' - Consider changing the name of cna() "correlation network analysis" to
#'         dccn() "dynamic cross-correlation analysis" or cnet() "correlation network"
#'
#' - Compare different clustering methods and look at the results of cluster.sim()
#' -- For other community detection methods see "modularity()" and
#'     fastgreedy.community(), spinglass.community(),
#'     leading.eigenvector.community(), edge.betweenness.community() for
#'     other community detection methods.
#' - Compare these results to those from VMD
#'
#' - Explore the compare.communities() and similarity.jaccard() functions!
#' - Explore the graph.intersection() and graph.difference() functions.
#'


#'
#' ## Preamble
#' This input file is a basic R script, called "cna_kinesin.r" that contains roxygen style comment lines. Lines that begin with #' will appear as document text (text chunks). Lines that begin with #+ or #- set attributes of code chunks.
#'   See:  http://yihui.name/knitr/demo/stitch/
#'
#'  Steps to make a PDF report include:

#+ preamble, include=TRUE, eval=FALSE
library(knitr)
spin('cna_kinesin.r')
system("pandoc -o cna_kinesin.pdf cna_kinesin.md")


#'
#' ## Example HIVp trajectory data
#' We begin by loading the cna package (eventually this will be part of the main bio3d package).


# Load the cna package and check its version number
library(cna)
packageVersion("cna")

#' Next we will load and supperpose the example HIVp trajectory data that ships with the bio3d package. 

#+ hivp_load, cache=TRUE
# Read example trajectory file
trtfile <- system.file("examples/hivp.dcd", package="bio3d")
trj <- read.dcd(trtfile, verbose=FALSE)

# Read the starting PDB file to determine atom correspondence
pdbfile <- system.file("examples/hivp.pdb", package="bio3d")
pdb <- read.pdb(pdbfile)

# Select residues 24 to 27 and 85 to 90 in both chains
inds <- atom.select(pdb,"///24:27,85:90///CA/", verbose=FALSE)

# Supperpose trajectory frames onto PDB
xyz <- fit.xyz(pdb$xyz, trj, fixed.inds=inds$xyz, mobile.inds=inds$xyz)

#' The object 'xyz' now contains our fitted trajectory data
dim(xyz) == dim(trj)

#+ rmsdTs, echo=FALSE, eval=FALSE
# RMSD of trj frames from PDB
r1 <- rmsd(a=pdb, b=xyz)
plot(r1, typ="o", ylab="RMSD", xlab="Frame Number")




#' ### Determine the cross-correlations of atomic displacements
#'
#' Now we have superposed data we can calculate the dynamical cross-correlation matrix for this C-alpha only trajectory with the dccm() function in bio3d. The dccm() function can be slow to run on larger trajectories. In such cases you may want to use multiple CPU cores for the calculation by setting the 'ncore' option (see ?dccm for details). 

#+ plotDCCM, cache=TRUE
## Calculate dynamical cross-correlation matrix  
cij <- dccm(xyz)
plot(cij)

#+ pymolCij, eval=FALSE
# View a representation of this dccm data in pymol
view.dccm(cij, pdb, launch=TRUE)

## < Inert Figure From File Here >


#' ### Correlation network construction and analysis
#'
#' Now we will use the cross-correlation data to partition our protein into dynamically consistent communities of residues (dynamic subdomains) by clustering an undirected network of residue nodes that are edge weighted by their correlations. Here we will use the 'cna()' function with its default -log(abs(cij[cij>0.4])) as the input adjacency matrix weights for our initial network. The function will also cluster this network into communities and produce a clustered network for further analysis.

#+ Network1, cache=TRUE
# Build and betweeness cluster a correlation network graph 
net <- cna(cij)
net

# Plot a simple network graph highlighting communities/clusters pattern
par(mfcol=c(1,2), mar=c(0,0,0,0))
plot(net$raw.network, vertex.size=5) ## Slow!!
plot.cna(net)

#' Note that the first plot is our (raw) residue wise network colored by community and second the (clustered) coarse grained community network, which has a node for each community that is scaled by the number of residues within that community. The coloring should be consistent but the orientation may vary.
#'

# If we use all residues layout from PDB coorsds we get an awesome view!!
col <- vmd.colors()[net$raw.communities$membership]
plot(net$raw.network, layout=layout.pdb(pdb, 1:sum(pdb$calpha),k=2),
     vertex.color=col, vertex.size=8 ) ## addedge.width and community color etc...

plot.communities(net$raw.communities, net$raw.network,
                 layout=layout.pdb(pdb, 1:sum(pdb$calpha),k=2),
                 vertex.color=col, vertex.size=8 )

plot.cna(net, layout=layout.pdb(pdb, net,k=2))
## Why are the colors different??


#' The residue-wise network (net$raw.network) is rather slow to plot as it nearly 2000 edges for this example, which obviously can obscure interoperation. The clustered network (net$clustered.network) aims to simplify this representation by showing a node per community with a scaled radius representing the number of residues (or nodes) within each community, as well as simplified weighted edges.

#'
#' ### NOTE: Get sugetions on output names from cna() and whether we can fix layouts for raw and clustered network plots - (one option is to take from 3D PDB coords?)
# Has the same colors as network plots and VMD output
plot.dccm2(cij, margin.segments=net$raw.communities$membership)


#' Now is probably a good time to take an interlude to talk about network graph visualization.
#' It's useful to think of there being two kinds of network visualizations: _pretty plots_ that encode 'soft' information and _ugly plots_ that encode 'hard' information. Most of what you see in journals etc falls under the first category.
#'
#' These 'soft' plots try to display a network so you can see broad structure and patterns in a totally non-rigourous way. Most of them are based on a sort of a physical metaphor: imagine the nodes of the graphs are little balls that are contantly repelling on another like two mis-aligned magnets, but somebody went and attached a bunch of them together with springs. You throw them on the table and see how they settle. This often gives really nice results, with vertices well spread out, and with highly connected vertices tending to be closer together. But it's also totally heuristic and very easy to read too much into. Plus each time you run it it may have a different layout.


#+ DCCM_with_grid
# Remove some of the small communities from the plot and add a grid
#source("http://bitbucket.org/Grantlab/bio3d/raw/ce6b7519f455b1094864b3ebe059984c27484d75/new_funs/add.dccm.grid.R")

##plot.dccm2(cij, margin.segments=net$raw.communities$membership, segment.min=25)
#plot.dccm2(cij, margin.segments=net$raw.communities$membership, segment.col=vmd.colors(), segment.min=25 )
#add.dccm.grid( net$raw.communities$membership, segment.min=25, fill.col="gray")

#'
#' ### NOTE: Need to improve the plot.dccm() function in general to help with the above task of annotating these matrix plots
#'

#' Output a PDB file for viewing community structure (i.e. clustering) that can be colored in VMD by chain to have the same colors as the previous plots

#+ eval=FALSE
# Use the chain field to store cluster membership data for color in VMD
ch <- vec2resno(vec=net$raw.communities$membership, resno=pdb$atom[,"resno"])
write.pdb(pdb, chain=LETTERS[ch], file="tmp.pdb")

## < Inert Figure From File Here >

#'
#' ### NOTE: We need a better summary.cna() function to simplify network interoperation
#' Lets examine the community structure
# Print and store a summary of whats in each community.
x <- summary(net)

attributes(x)
x$size
x$members$"3"
x$members[[3]]





#'
#' ### Collapse and/or prune nodes. 
#'
#' ### NOTE: We should have an option to not plot communities under a certain size  or with less less than a certain number of edges to other communities. For example, a community composed of only 2 residues that is not linked to any other community).
#'
#'### NOTE: Need to fix and formalize this kind of prune.cna() function
#'
# Define a function to prune network
prune.cna <- function(x, edges.min=1) {
  ##-- Prune nodes with less than 'edge.min' edges
  ##     prune.cna(net)

  network=x$clustered.network
  ## Identify nodes with less than 'edges.min' to other nodes.
  zero.inds <- degree(network) < edges.min

  ## Bug Input 'size.min=2' size are scaled V(net$network)$size
  ## Identify nodes with size less than 'size.min'
  ##zero.inds <- c(zero.inds, (sizes(x) > size.min))

  zero.vs <- V(network)[zero.inds]
  cat( paste("Removing Nodes:", paste(zero.vs, collapse=", ")),"\n")
  ## Should be able to print residue membership of these nodes that will be removed
  d <- delete.vertices(network, zero.vs)

  comms <- edge.betweenness.community(d, directed = FALSE)
  output <- list(prune.network=d, prune.communities=comms)
  class(output) <- c("community", "cna")
  return(output)
}

#+ PruneNetwork
dnet <- prune.cna(net)
par(mfcol=c(1,2), mar=c(0,0,0,0))
plot(net)
plot(dnet$prune.network)

#' ### NOTE: Need to fix this mess of names for networks $raw... $clustered... and now $prune.... Especialy neeed to do this to help plot and summary functions etc...

#' You can store whatever annotation you want along with our network graph for example:
net$clustered.network$title="Dynamic Network of HIVp"
list.graph.attributes(net$clustered.network)
net$clustered.network$title


#' ### NOTE: We actually want to zoom through at different clustering levels to show a hierarchy of collective coupled motions and dynamic domains.
#'
#' There is also the igraph function 'contract.vertices()'
#'
#' ### NOTE: Need to double check that cna.collapse is working coreclty as using 'mean' and 'median' as collapse methods give an unconected network!
## cnet.ave <- cna.collapse(net, "mean")
## xy <- plot.cna(cnet.ave)

#+ Cuttoff0_48
# Or change the cutoff and examine effects the conections
##xy <- plot.cna( cna.collapse(net, "max", cutoff.cij=0.48) )

#'
#' ### Identifying nodes/communities and their composite residues
#'
#' ### NOTE: How do we get the coords of the nodes so we can use _identify()_ function or add extra labels to the plot in some useful way? 

## we can get the coordinates but how to use them is still not clear
##xy<-plot.cna(cnet)
xy
layout(xy) #??

#' See ?igraph.plotting


#'
#' ## Networks as 3D plots
#' We could take the C.O.M. for each community and use as 3D coords for community graph
#'
#' To do this we need to"
#' - get residues within each community
#' - get their coords
#' - get their mean values
#' - scale these values to plotting frame
#' - try passing to rglplot along with node widths 


#+ plot3D, eval=TRUE, custom.plot=TRUE
##rglplot(net$clustered.network)
knit_hooks$set(custom.plot=hook_plot_custom)
rglplot(dnet$prune.network)
rgl.viewpoint(10,10/4)
## paste(fig_path(1),"png", sep=".") # Need a knit_hook for this
rgl.snapshot( paste(fig_path(),"png", sep=".") )
#ll <- layout.fruchterman.reingold(dnet$network,dim=3)
#rglplot(dnet$network,layout=ll)
#rgl.quit()


#' Note: par3d can obtain the projection matrix ‘P = par3d("projMatrix")’
#' Can use par3d etc. to orient plot and even make an animation.


# A 90 degree rotation about the x axis:
##rotationMatrix(pi/2, 1, 0, 0)

#'
#' ## Try the contract.vertices() function to replace cna.collapse 
#'

# Use igraphs contract.vertices() function
tmp1 <- simplify(contract.vertices(net$raw.network, mapping=net$raw.communities$membership))
par(mfcol=c(1,3), mar=c(0,0,0,0))
plot(tmp1)
## V(tmp1) ## names of nodes now in communities
V(tmp1)$name = LETTERS[1:vcount(tmp1)]
V(tmp1)$size = sizes(net$raw.communities)
## E(tmp1)$weight
plot(tmp1)
plot(tmp1, edge.width=E(tmp1)$weight)

#' And lets add communities
eb<-edge.betweenness.community(tmp1)
plot(eb, tmp1, edge.width=E(tmp1)$weight)

#' Lets paly with 'attribute.combination' 
# Second try collapsing names
tmp2 <- simplify(contract.vertices(net$raw.network, mapping=net$raw.communities$membership,
                        list(weight="mean", name=toString, "ignore") ))
par(mfcol=c(1,3), mar=c(0,0,0,0))
V(tmp2)$name = LETTERS[1:vcount(tmp1)]
plot(tmp2)
## V(tmp2) ## names of nodes now in communities
V(tmp2)$size = sizes(net$raw.communities)
## E(tmp2)$weight
plot(tmp2)
plot(tmp2, edge.width=E(tmp2)$weight)

#' Need to see the help page for ?attribute.combination to figure all this out

#sizes(tmp1)
#neighbors(tmp1, 12)



#'
#' ## Contact map filtering
#' We can also filter the dccm by a given contact map similar to what was done in the PNAS paper from the Luthy-Shulten group. However, this is probably not a good idea for HIVp as the cmap is so sparse...

#+ cmapHIVp
## Contact map
cm <- cmap( xyz, pcut=0.75, scut=0, dcut=5, mask.lower=FALSE)
plot.dccm(cm, main="cm")


#' The igraph.comms() function can take as input an OPTIONAL contact map filter (called 'cm').


# layout.mds looks cool but does it take weight into account??
xy <- layout.mds(net$clustered.network)
plot.cna(net, layout=xy)


#' Some cluster dendogram
dend <- as.dendrogram(cnet, use.modularity=TRUE)
plot(dend, nodePar=list(pch=c(NA, 20)))



# ?layout.fruchterman.reingold
l <- layout.fruchterman.reingold(net$clustered.network,niter=500,
                                 area=vcount(net$clustered.network)^2.3,
                                 repulserad=vcount(net$clustered.network)^2.8)
                                        #weights=E(net$clustered.network)$weight)
par(mfcol=c(1,2), mar=c(0,0,0,0))
plot(net$clustered.network, layout=l)
plot(net$clustered.network, layout=layout.mds(net$clustered.network), edge.width=E(net$clustered.network)$weight*15)


#h <- comms.collaps(g)
#plot(h)


#'
#' # Still to Explore
#' - graph.intersection:    Creates the intersection of two or more graphs (only edges present#'                          in all graphs will be included. The corresponding operator is %s%.
#'
#' - graph.difference:      Creates the difference of two graphs. Only edges present in the
#'                          first graph but not in the second will be be included in the new
#'                          graph. The corresponding operator is %m%.
#'
#' - graph.union:          Creates the union of two or more graphs.
#'
#' - compare.communities     Compares community structures using various metrics
#'                            compare(sg, le, method="nmi")
#'
#' - dendPlot                Plot dendrograms of community clustering
#'                            karate <- graph.famous("Zachary")
#'                            eb=edge.betweenness.community(karate); dendPlot(eb)
#'                            ?dendPlot.communities
#'
#' - communities           Functions to deal with the result of network community detection
#'                         as.hclust(x, hang = -1,use.modularity = FALSE, ...)
#'                         cutat(communities, no, steps)
#'                         length(communities) number of communities
#'                         size(communities)

#' ## Need to explore the layout functions more and how to use them
#' With the example code in layout we can actually use identify()!!
#'
#' - layout                  Generate coordinates for plotting graphs
#' - layout.mds              Graph layout by multidimensional scaling
#' - layout.svd(graph, d=shortest.paths(graph), ...) d=> distance matrix of graph!!
#'     ‘layout.svd’ is a currently _experimental_ layout function based
#'     on singular value decomposition. It does not have the usual
#'     ‘params’ argument, but take a single argument, the distance matrix
#'     of the graph. This function generates the layout separately for
#'     each graph component and then merges them via ‘layout.merge’.

#'
#' - clusters
#'
#' - articulation.points:  vertex whose removal disconnects the graph
#'
#' - graph.maxflow           Maximum flow in a network


#' It is possibly to call the ‘plot’ function on ‘communities’ objects. This will plot the graph (and uses ‘plot.igraph’ internally), with the communities shown. By default it colores the vertices according to their communities, and also marks the vertex groups corresponding to the communities. It passes additional arguments to ‘plot.igraph’, please see that and also ‘igraph.plotting’ on how to change the plot.
#'
#' See ‘dendPlot’ for plotting community structure dendrograms.
#'
#' See ‘compare.communities’ for comparing two community structures on the same graph.
#'
#' ## Network Flow 
#' - assortativity:         The assortativity coefficient is the Pearson correlation coefficient of degree between pairs of linked nodes
#' - attributes
#' - authority.score:         Kleinberg's centrality scores
#' - cliques                 The functions find cliques, ie. complete subgraphs in a graph
#' - cohesive.blocks         Calculate Cohesive Blocks

#' - contract.vertices       Contract several vertices into a single one
#' - community.to.membership Common functions supporting community detection algorithms
#' - dendPlot                Plot dendrograms of community clustering
#' - dendPlot.communities    Community structure dendrogram plots
#' - graph.de.bruijn         De Bruijn graphs.

#' - graph.motifs            Graph motifs
#' - igraph.plotting         Drawing graphs

#' - leading.eigenvector.community  Community structure detecting based on the leading eigenvector of the community matrix
#' - minimum.spanning.tree   Minimum spanning tree
#' - modularity              Modularity of a community structure of a graph
#' - neighborhood            Neighborhood of graph vertices
#' - optimal.community       Optimal community structure
#' - print.igraph            Print graphs to the terminal
#' - read.graph              Reading foreign file formats
#' - revolver                Measuring the driving force in evolving networks
#' - rewire                  Graph rewiring
#' - shortest.paths          Shortest (directed or undirected) paths between vertices
#' - similarity.jaccard      Similarity measures of two vertices


#'
#' ## Working with example kinesin data
#' We will begin with example MD trajectory output, calculate correlation matrices, contact residencey maps and then build simple networks from this data and investigate their properties.
#'   
#' Looks like the sse options on this data are messed up - can you find out why please?
#' 

## Load kinesin.pdb, kin-trj1.ncdf, kin-trj2.ncdf, kin-trj3.ncdf,
##pth <- system.file("data/", package="cna")
##pdb <- read.pdb( paste0(pth, "kinesin.pdb") )
##trj.1 <- read.ncdf( paste0(pth, "kin-trj1.ncdf") )

#+ atom.selectKin
## Setup selection indices SSE definition
#noh.inds <- atom.select(pdb, "noh")
#ca.inds <- atom.select(pdb, "calpha")
#sse <- dssp(pdb)


#'------------------
#' ## Other output of use??
#' Need output for viewing in VMD - network and communities
#' Need plot of 2d community graph - Please include these functions ASAP!!

#'
#' ## Information About the Current R Session
print(sessionInfo(), FALSE)



cnet <- NULL
cnet$network = net$clustered.network

# Problem is node size is scaled from cna.collapse
list.vertex.attributes(cnet$network)
V(cnet$network)$size
sizes(net)
degree(cnet$network)

# examples of soft plots:
ug <- cnet$network
plot(ug,layout=layout.fruchterman.reingold)
plot(ug,layout=layout.fruchterman.reingold.grid)
plot(ug,layout=layout.kamada.kawai)

#' 'Hard' plots try to represent some underlying feature of the data through the layout. These are often based on the idea of trying to find a way to embed an n-dimensional object in two dimensions, preserving as much information as possible.

# examples of hard plots:
plot(ug,layout=layout.svd)
plot(ug,layout=layout.mds)


#' igraph is fantastic with plotting graphs. Everything is customizable.
#' help('igraph.plotting') is your friend
#' a major feature of this is the huge library of layout algorithms it has
#' use help('layout.random') to see a list of these.

# layout can be specified a few ways
myLayout <- layout.fruchterman.reingold(ug)
myLayout
plot(ug,layout=myLayout)

plot(ug,layout=layout.fruchterman.reingold)

ug$layout=layout.fruchterman.reingold
plot(ug)
plot(ug)

#' you can make your own layout algoritms or build a layout matrix by hand

#' the easiest way to change the attributes of vertices/edges is to set vertex/edge attributes
V(ug)$color <- V(ug)$favColor
V(ug)$size <- V(ug)$favNumber+3
V(ug)$label <- c('alice','bob','carol','derek','erica','anon','anon')
plot(ug)
V(ug)$frame.color <- "#00000000"

list.edge.attributes(ug)
E(ug)$width=E(ug)$weight*4
E(ug)[type==0]$color <- 'dark orange'
E(ug)[type==1]$color <- 'purple'
E(ug)[type==0]$label <- 'love'
E(ug)[type==1]$label <- 'hate'
E(ug)$label.cex <- .7
E(ug)$curved <- TRUE
plot(ug)

#' and now you have the world's most garish network plot. (except sadly that's not true)
