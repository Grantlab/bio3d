###################################################
### chunk number 1: 
###################################################
#line 33 "Bio3D_trajectory.Rnw"
library(bio3d)
lbio3d()


###################################################
### chunk number 2: 
###################################################
#line 44 "Bio3D_trajectory.Rnw"
dcdfile <- system.file("examples/hivp.dcd", package="bio3d")
pdbfile <- system.file("examples/hivp.pdb", package="bio3d")


###################################################
### chunk number 3: 
###################################################
#line 50 "Bio3D_trajectory.Rnw"
mydcdfile <- "/path/to/my/data/myfile.dcd"


###################################################
### chunk number 4: 
###################################################
#line 54 "Bio3D_trajectory.Rnw"
dcd <- read.dcd(dcdfile)
pdb <-  read.pdb(pdbfile)


###################################################
### chunk number 5: 
###################################################
#line 60 "Bio3D_trajectory.Rnw"
pdb.summary(pdb)
length(pdb$xyz)
dim(dcd)


###################################################
### chunk number 6: 
###################################################
#line 69 "Bio3D_trajectory.Rnw"
ca.inds <- atom.select(pdb, elety="CA")


###################################################
### chunk number 7: 
###################################################
#line 73 "Bio3D_trajectory.Rnw"
xyz <- fit.xyz(fixed=pdb$xyz, mobile=dcd,
               fixed.inds=ca.inds$xyz,
               mobile.inds=ca.inds$xyz)


###################################################
### chunk number 8: 
###################################################
#line 83 "Bio3D_trajectory.Rnw"
rd <- rmsd(xyz[1,ca.inds$xyz], xyz[,ca.inds$xyz])

plot(rd, typ="l", ylab="RMSD", xlab="Frame No.")
points(lowess(rd), typ="l", col="red", lty=2, lwd=2)


###################################################
### chunk number 9: 
###################################################
#line 91 "Bio3D_trajectory.Rnw"
hist(rd, breaks=40, freq=FALSE, main="RMSD Histogram")
lines(density(rd), col="gray", lwd=3)


###################################################
### chunk number 10: 
###################################################
#line 96 "Bio3D_trajectory.Rnw"
summary(rd)


###################################################
### chunk number 11: 
###################################################
#line 102 "Bio3D_trajectory.Rnw"
rf <- rmsf(xyz[,ca.inds$xyz])
plot(rf, ylab="RMSF", xlab="Residue Position", typ="l")


###################################################
### chunk number 12: 
###################################################
#line 113 "Bio3D_trajectory.Rnw"
pc <- pca.xyz(xyz[,ca.inds$xyz])
plot(pc, col=bwr.colors(nrow(xyz)) )


###################################################
### chunk number 13: 
###################################################
#line 119 "Bio3D_trajectory.Rnw"
hc <- hclust(dist(pc$z[,1:2]))
grps <- cutree(hc, k=2)
plot(pc, col=grps)


###################################################
### chunk number 14: 
###################################################
#line 127 "Bio3D_trajectory.Rnw"
plot.bio3d(pc$au[,1], ylab="PC1 (A)", xlab="Residue Position", typ="l")
points(pc$au[,2], typ="l", col="blue")


###################################################
### chunk number 15: 
###################################################
#line 133 "Bio3D_trajectory.Rnw"
p1 <- mktrj.pca(pc, pc=1, b=pc$au[,1], file="pc1.pdb")
p2 <- mktrj.pca(pc, pc=2,b=pc$au[,2], file="pc2.pdb")


###################################################
### chunk number 16:  eval=FALSE
###################################################
## #line 139 "Bio3D_trajectory.Rnw"
## write.ncdf(p1, "trj_pc1.nc")


###################################################
### chunk number 17: 
###################################################
#line 150 "Bio3D_trajectory.Rnw"
cij<-dccm(xyz[,ca.inds$xyz])
plot.dccm(cij)


###################################################
### chunk number 18: 
###################################################
#line 156 "Bio3D_trajectory.Rnw"
toLatex(sessionInfo())


