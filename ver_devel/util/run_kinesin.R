library(bio3d)
library(testthat)
attach(kinesin)

pdbs = read.fasta.pdb(pdbs, prefix="scop_pdbs/", pdbext=".ent", ncore=8)

core <- core.find(pdbs)

xyz = fit.xyz(pdbs$xyz[1,], pdbs, core$c1A.xyz, core$c1A.xyz)
pdbs$xyz <- xyz

try(expect_identical(annotation, kinesin$annotation))
try(expect_equal(pdbs, kinesin$pdbs))
try(expect_equal(core, kinesin$core))
#expect_equivalent(pdbs$ali, kinesin$pdbs$ali)
#expect_equivalent(core$c1A.xyz, kinesin$core$c1A.xyz)
#expect_equivalent(core$c0.5A.xyz, kinesin$core$c0.5A.xyz)
#expect_equal(as.vector(pdbs$xyz), as.vector(kinesin$pdbs$xyz), tolerance=1.e-5)

kinesin = list(pdbs=pdbs, core=core, annotation=annotation)

save(kinesin, file="kinesin.RData")
