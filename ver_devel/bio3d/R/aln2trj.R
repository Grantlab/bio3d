# Find equivalent positions between alignment and trajectory
"aln2trj" <-
function(aln, pdb, inds, id="seq", file="aln2trj.fa") {
   # build reference resno
   resno.ref <- rep(NA, ncol(aln$ali))
   resno.ref[aln$ali[1, ]!="-"] <- 1:length(which(aln$ali[1, ]!="-"))
   aa1 <- seq.pdb(pdb)
   findid <- grep(id, aln$id)
   if(length(findid) == 0) {
      naln <- seq2aln(seq2add=aa1, aln=aln, id=id, file=file)
      idhit <- nrow(naln$ali)
   } else {
      if(length(findid) > 1) {
         warning(paste("Multiple entities found in alignment for id=", id, 
               ". Use the first one...", sep=""))
      }
      naln <- aln 
      idhit <- findid[1]
   }
   # build reference resno after realignment
   nresno.ref <- rep(NA, ncol(naln$ali))
   nresno.ref[naln$ali[1, ]!="-"] <- resno.ref[!is.na(resno.ref)]
   # resno for PDB after realignment
   nresno <- rep(NA, ncol(naln$ali))
   ca.inds <- atom.select(pdb, elety="CA") 
   nresno[naln$ali[idhit, ]!="-"] <- pdb$atom[ca.inds$atom, "resno"]
   
   if(!is.list(inds)) inds <- list(inds)
   new.inds.all <- lapply(inds, function(i){
      new.inds.ref <- which(nresno.ref %in% resno.ref[i])
      nresno.i <- nresno[new.inds.ref]
      if(any(is.na(nresno.i))) {
         warning("Gaps are found in finding equivalent positions in PDB. Ignore gaps...")
      }
      new.inds <- atom.select(pdb, resno=nresno.i[!is.na(nresno.i)], elety="CA")
   } )
   if(length(new.inds.all) == 1) new.inds.all <- new.inds.all[[1]]
   return (new.inds.all)
}
