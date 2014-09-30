library(bio3d)
library(testthat)

# If update=FALSE, start from the very beginning...
update = TRUE

if(!update) {
   aa <- get.seq("1tag")
   blast <- blast.pdb(aa)
   ids <- plot.blast(blast, cutoff=332)
} else {
   ids <- basename(transducin$pdbs$id)
   annotation <- transducin$annotation
} 

files <- get.pdb(ids, split = TRUE, ncore=8)
pdbs <- pdbaln(files, fit=TRUE, ncore=8)

core <- core.find(pdbs)
xyz <- fit.xyz(pdbs$xyz[1,], pdbs, core$c1A.xyz, core$c1A.xyz)
pdbs$xyz <- xyz
pdbs$id <- ids
rownames(pdbs$ali) <- ids

expect_identical(annotation, transducin$annotation)
expect_equivalent(pdbs$ali, transducin$pdbs$ali)
expect_equivalent(core$c1A.xyz, transducin$core$c1A.xyz)
expect_equivalent(core$c0.5A.xyz, transducin$core$c0.5A.xyz)
expect_equal(as.vector(pdbs$xyz), as.vector(transducin$pdbs$xyz), tolerance=1.e-5)

transducin = list(pdbs=pdbs, core=core, annotation=annotation)

save(transducin, file="transducin.RData")
