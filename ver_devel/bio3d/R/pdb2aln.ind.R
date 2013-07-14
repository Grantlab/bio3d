#  Find the best alignment between a PDB structure and an 
#  existing alignment. Then, given a set of residue indices
#  defined for the original alignment, return the equivalent 
#  CA atom indices in the PDB coordinates.
"pdb2aln.ind" <-
function(aln, pdb, inds, ...) {
   # get the new alignment
   naln <- pdb2aln(aln=aln, pdb=pdb, ...)

   ninds <- which(naln$ref["ali.pos",] %in% inds)
   ca.inds <- naln$ref["ca.inds", ninds]

   if(any(is.na(ca.inds))) {
      warning("Gaps are found in equivalent positions in PDB")
   }

   return(ca.inds)
}
