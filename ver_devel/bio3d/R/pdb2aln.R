#  Find the best alignment between a PDB structure and an 
#  existing alignment. Then, given a set of residue indices
#  defined for the original alignment, return the equivalent 
#  CA atom indices in the PDB coordinates.
"pdb2aln" <-
function(aln, pdb, inds, aln.id=NULL, id="seq", file="pdb2aln.fa") {
   # build reference resno
   resno.ref <- rep(NA, ncol(aln$ali))
   ii <- which(!is.gap(aln$ali[1, ]))
   resno.ref[ii] <- 1:length(ii)
   aa1 <- seq.pdb(pdb)
   
   if(!is.null(aln.id)) findid <- grep(aln.id, aln$id)
   if(is.null(aln.id) || length(findid)==0) {
      # do sequence-profile alignment
      naln <- seq2aln(seq2add=aa1, aln=aln, id=id, file=file)
   } else {
      # do pairwise sequence alignment     
      if(length(findid) > 1) {
         warning(paste("Multiple entities found in alignment for id=", 
            aln.id, ". Use the first one...", sep=""))
      }
      idhit <- findid[1]

      ##- Align seq to masked template from alignment
      tmp.msk <- aln$ali[idhit, ]
      tmp.msk[is.gap(tmp.msk)] <- "X"
      seq2tmp <- seqaln.pair(seqbind(aa1, tmp.msk), id=c("seq","tmp"), file=file)
      
      ##- check sequence identity
      ii <- which(seq2tmp$ali[2,]=='X' & is.gap(seq2tmp$ali[1,]))
      idt <- identity(seq2tmp$ali[, -ii])[1,2]
      if(idt < 0.4) {
         warning(paste("Sequence identity is too low (<40%).",
             "You may want profile alignment (aln.id=NULL)", sep=" "))
      }
 
      ##- Insert gaps to adjust alignment
      ins <- which(is.gap( seq2tmp$ali[2,] ))
      if( length(ins)==0 ) {
        ntmp <- aln$ali
      } else {
        ntmp <- matrix("-", nrow=nrow(aln$ali), ncol=(ncol(aln$ali)+length(ins)))
        ntmp[,-ins] <- aln$ali
      }
    
      ## Add seq to bottom of adjusted alignment
      naln.ali <- seqbind(ntmp, seq2tmp$ali[1,])
      rownames(naln.ali) <- c(rownames(aln$ali), id)
      naln <- list(id=c(aln$id, id), ali=naln.ali) 
   }

   # build reference resno after realignment
   nresno.ref <- rep(NA, ncol(naln$ali))
   nresno.ref[!is.gap(naln$ali[1, ])] <- resno.ref[!is.na(resno.ref)]
   # resno for PDB after realignment
   nresno <- rep(NA, ncol(naln$ali))
   ca.inds <- atom.select(pdb, elety="CA") 
   nresno[!is.gap(naln$ali[nrow(naln$ali), ])] <- pdb$atom[ca.inds$atom, "resno"]
  
   # use list to get indices in batch 
   if(!is.list(inds)) inds <- list(inds)
   new.inds.all <- lapply(inds, function(i){
      new.inds.ref <- which(nresno.ref %in% resno.ref[i])
      nresno.i <- nresno[new.inds.ref]
      if(any(is.na(nresno.i))) {
         warning(paste("Gaps are found in finding equivalent positions in PDB.",
                  "Ignore gaps...", sep=" "))
      }
      new.inds <- atom.select(pdb, resno=nresno.i[!is.na(nresno.i)], elety="CA")
   } )
   if(length(new.inds.all) == 1) new.inds.all <- new.inds.all[[1]]
   return (new.inds.all)
}
