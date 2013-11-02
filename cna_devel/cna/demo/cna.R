## Package Demo
options(digits=2)
require(cna)


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

## Calculate dynamical cross-correlation matrix  
cij <- dccm(xyz)
plot(cij)

#+ pymolCij, eval=FALSE
# View a representation of this dccm data in pymol
view.dccm(cij, pdb, launch=TRUE)

net <- cna(cij)
net

# Plot a simple network graph highlighting communities/clusters pattern
##par(mfcol=c(1,2), mar=c(0,0,0,0))
##plot(net$raw.network, vertex.size=5) ## Slow!!
plot.cna(net)

#col <- vmd.colors()[net$raw.communities$membership]
#layout <- layout.pdb(pdb, 1:sum(pdb$calpha),k=2)
#plot(net$raw.network, layout=layout,
#     vertex.color=col, vertex.size=10 )

plot.cna(net, layout=layout.pdb(pdb, net,k=2))

plot.dccm2(cij, margin.segments=net$raw.communities$membership)

view.cna(net, alpha=0.7, launch=TRUE)

plot3d.cna(net, pdb)

