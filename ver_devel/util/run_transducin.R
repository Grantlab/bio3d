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
rownames(pdbs$xyz) <- ids
rownames(pdbs$resno) <- ids
rownames(pdbs$resid) <- ids
rownames(pdbs$b) <- ids
rownames(pdbs$chain) <- ids
rownames(pdbs$sse) <- ids

try(expect_identical(annotation, transducin$annotation))
try(expect_equal(pdbs, transducin$pdbs))
try(expect_equal(core, transducin$core))
#expect_equivalent(pdbs$ali, transducin$pdbs$ali)
#expect_equivalent(pdbs$xyz, transducin$pdbs$xyz)
#expect_equivalent(pdbs$resno, transducin$pdbs$resno)
#expect_equivalent(pdbs$resid, transducin$pdbs$resid)
#expect_equivalent(core$c1A.xyz, transducin$core$c1A.xyz)
#expect_equivalent(core$c0.5A.xyz, transducin$core$c0.5A.xyz)
#expect_equal(as.vector(pdbs$xyz), as.vector(transducin$pdbs$xyz), tolerance=1.e-5)

transducin = list(pdbs=pdbs, core=core, annotation=annotation)

save(transducin, file="transducin.RData")
