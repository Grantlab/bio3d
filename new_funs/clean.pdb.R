"clean.pdb" <-
function(pdb, consecutive=TRUE, force.renumber = FALSE, fix.chain = FALSE, 
    fix.aa = FALSE, rm.wat = FALSE, rm.lig = FALSE, rm.h = FALSE, verbose=TRUE) {

  if(!is.pdb(pdb)) 
     stop("Input should be a 'pdb' object")

  cl <- match.call()

  ## processing message
  ## stored as an N-by-3 matrix with columns: 
  ##     FACT, OPERATION, IMPORTANT NOTE
  msg <- NULL

  ## Recognized amino acid names
  prot.aa <- c("ALA", "CYS", "ASP", "GLU",
               "PHE", "GLY", "HIS", "ILE", "LYS", "LEU", "MET", "ASN",
               "PRO", "GLN", "ARG", "SER", "THR", "VAL", "TRP", "TYR",
               "SEP", "TPO", "MLY", "MSE", "IAS", "ABA", "CSO", "CSD",
               "CYM", "CME", "CSX", "CMT", "CYX", "HIE", "HIP", "HID",
               "HSD", "HSE", "HSP", "DDE", "MHO", "ASX", "CIR", "PFF")
  
  ## for residues and atoms renumbering
  first.eleno = 1
  first.resno = 1

  ## remove water
  if(rm.wat) {
     wat <- atom.select(pdb, "water", verbose = FALSE)
     if(length(wat$atom) > 0) {
        pdb$atom <- pdb$atom[-wat$atom, ]
        pdb$xyz <- pdb$xyz[, -wat$xyz]
        msg <- rbind(msg, c(paste("Found", length(wat$atom), "water atoms"), "REMOVED", ""))
     }
  }

  ## remove ligands 
  if(rm.lig) {
     lig <- atom.select(pdb, "ligand", verbose = FALSE)
     if(length(lig$atom) > 0) {
        pdb$atom <- pdb$atom[-lig$atom, ]
        pdb$xyz <- pdb$xyz[, -lig$xyz]
        msg <- rbind(msg, c(paste("Found", length(lig$atom), "ligand atoms"), "REMOVED", ""))
     }
  }

  ## remove hydrogens
  if(rm.h) {
     h.inds <- atom.select(pdb, "h", verbose = FALSE)
     if(length(h.inds$atom) > 0) {
        pdb$atom <- pdb$atom[-h.inds$atom, ]
        pdb$xyz <- pdb$xyz[, -h.inds$xyz]
        msg <- rbind(msg, c(paste("Found", length(h.inds$atom), "hydrogens"), "REMOVED", ""))
     }
  }
 
  ## check if 'alt' coords exist
  if(any(rm.p <- !is.na(pdb$atom$alt))) {
     pdb$atom <- pdb$atom[!rm.p, ]
     pdb$xyz <- pdb$xyz[, -atom2xyz(which(rm.p))]
     msg <- rbind(msg, c(paste("Found", sum(rm.p), "ALT records"), "REMOVED", ""))
  }

  ## following operations are on an independent object
  npdb <- pdb

  ca.inds <- atom.select(npdb, "calpha", verbose = FALSE)

  ## check chain breaks and missing chain ids
  has.fixed.chain <- FALSE
  capture.output( new.chain <- chain.pdb(npdb) )
  chain <- npdb$atom[, "chain"]  
  if(any(is.na(npdb$atom[, "chain"]))) {
     msg <- rbind(msg, c("Found empty chain IDs", 
               ifelse(fix.chain, "FIXED", "NO CHANGE"), 
               ifelse(fix.chain, "ALL CHAINS ARE RELABELED", "")) )
     if(fix.chain) {
        npdb$atom[, "chain"] <- new.chain
        has.fixed.chain <- TRUE 
     } 
  }
  else {
     ## check if new chain id assignment is consistent to original one
     chn.brk <- bounds(chain[ca.inds$atom], dup.inds=TRUE, pre.sort=FALSE)
     new.chn.brk <- bounds(new.chain[ca.inds$atom], dup.inds=TRUE, pre.sort=FALSE)
     if(!isTRUE(all.equal(chn.brk, new.chn.brk))) {
        msg <- rbind(msg, c("Found inconsistent chain breaks",
                 ifelse(fix.chain, "FIXED", "NO CHANGE"),
                 ifelse(fix.chain, "ALL CHAINS ARE RELABELED", "")))
     
        msg <- rbind(msg, c("Original chain breaks:", "", ""))
        pre.ca <- ca.inds$atom[chn.brk[, "end"]]
        pre.msg <- capture.output(npdb$atom[pre.ca, c("resid", "resno", "chain")])
        msg <- rbind(msg, cbind(pre.msg, "", ""))
        msg <- rbind(msg, c("", "", ""))
 
        msg <- rbind(msg, c("New chain breaks:", "", "") )
        new.ca <- ca.inds$atom[new.chn.brk[, "end"]]
        new.msg <- capture.output(npdb$atom[new.ca, c("resid", "resno", "chain")])
        msg <- rbind(msg, cbind(new.msg, "", ""))

        if(fix.chain) {
           npdb$atom[, "chain"] <- new.chain
           has.fixed.chain <- TRUE 
        }
     } 
  }
   
  ## Renumber residues and atoms
  if( any(!is.na(npdb$atom[, "insert"])) ) {
     renumber <- TRUE
     msg <- rbind(msg, c("Found INSERT records", "RENUMBERED", ""))
     npdb$atom[, "insert"] <- as.character(NA)
  }
  else if(force.renumber) {
     renumber <- TRUE 
     msg <- rbind(msg, c("force.renumber = TRUE", "RENUMBERED", ""))
  }
  else {
     renumber <- FALSE
  }
  if(renumber) {

    ## Assign consecutive atom numbers 
    npdb$atom[,"eleno"] <- seq(first.eleno, length=nrow(npdb$atom))

    ## Determine what chain ID we have
    chain <- unique(npdb$atom[, "chain"])
    
    ##- Assign new (consecutive) residue numbers for each chain
    prev.chain.res = 0  ## Number of residues in previous chain
    for(i in 1:length(chain)) {
       inds <- which(npdb$atom[, "chain"] == chain[i])
       
       ## Combination of chain id, resno and insert code uniquely defines a residue (wwpdb.org)
       ## Here we use original pdb because we assume it should at least 
       ## distinguish different residues by above combination.
       ## We don't use the modified pdb (npdb) because all non-protein residues 
       ## are assigned a chain ID as "X" after calling chain.pdb(); 
       ## These residues could have the same resno (which are still in original 
       ## form) as they may be assigned different chain IDs in the original pdb. 
       res <- paste(pdb$atom[inds, "chain"], pdb$atom[inds, "resno"], 
                    pdb$atom[inds, "insert"], sep="_")

       n.chain.res <- length(unique(res))

       new.nums <- (first.resno+prev.chain.res):(first.resno+n.chain.res-1+prev.chain.res)
       npdb$atom[inds, "resno"] <- vec2resno(new.nums, res)

       if(consecutive) {
         ## Update prev.chain.res for next iteration 
         prev.chain.res = prev.chain.res + n.chain.res
       }
    }
  } 
  else {
     if(has.fixed.chain) {
        warning(paste("Chain IDs are relabeled but residues are not renumbered.",
                    "Overlappling residue numbers may exist in the same chain!", sep="\n"))
     }
  }

  ## update SSE
  if(length(pdb$helix$start) > 0 || length(pdb$sheet$start) > 0) {
     if(has.fixed.chain || renumber) {
        sse <- .pdb2sse(pdb)
        names(sse) <- sub(".*_.*_(.*)", "\\1", names(sse))
        names(sse) <- paste(npdb$atom[ca.inds$atom, "resno"], 
               npdb$atom[ca.inds$atom, "chain"], names(sse), sep="_")
        sse.ind <- .bounds.sse(sse)
#        if(length(sse.ind$helix$start) == length(pdb$helix$start))
#           sse.ind$helix$type <- pdb$helix$type
#        if(length(sse.ind$sheet$start) == length(pdb$sheet$start))
#           sse.ind$sheet$sense <- pdb$sheet$sense
           
        npdb$helix <- sse.ind$helix
        npdb$sheet <- sse.ind$sheet
 
        msg <- rbind(msg, c("SSE annotation", "UPDATED", ""))
     }
  }

  ## update amino acid name
  naa.atom <- which(npdb$atom[, "resid"] %in% prot.aa[-c(1:20)])
  naa.res <- intersect(ca.inds$atom, naa.atom)
  if(length(naa.res) > 0) {
     msg <- rbind(msg, c(paste("Found", length(naa.res), "non-standard amino acids"), 
              ifelse(fix.aa, "FIXED", "NO CHANGE"), 
              ifelse(fix.aa, "AMINO ACID NAMES ARE CHANGED", "")) )
     if(fix.aa) {
         npdb$atom[naa.atom, "resid"] <- aa123(aa321(npdb$atom[naa.atom, "resid"])) 
     }
  }

  ## update pdb$calpha
  npdb$calpha <- seq(1, nrow(npdb$atom)) %in% ca.inds$atom
  if(!identical(pdb$calpha, npdb$calpha)) 
     msg <- rbind(msg, c("The component calpha", "UPDATED", ""))

  ## update pdb$call
  npdb$call <- cl

  ## update class
  if(!inherits(npdb, "pdb") || !inherits(npdb, "sse")) {
     class(npdb) <- c("pdb", "sse")
     msg <- rbind(msg, c("PDB object class", "UPDATED", ""))
  }

  ## log
  if(!is.null(msg)) {
     chk.log <- apply(msg, 1, function(x) {
         if(nchar(x[2]) == 0)
            sprintf("%-40s", x[1])
         else if(nchar(x[3]) == 0)
            paste(sprintf("%-40s", x[1]), "->", sprintf("%-15s", x[2]))
         else
            paste(paste(sprintf("%-40s", x[1]), "->", sprintf("%-15s", x[2]), "!!", x[3], "!!") )
     } )
#     chk.log <- paste(paste(chk.log, collapse = "\n"), "\n")
     if(verbose) print(chk.log)
     npdb$chk.log <- chk.log
  }

  return(npdb)
} 
     
.pdb2sse <- function(pdb) {
  ##- Function to obtain an SSE sequence vector from a PDB object
  ##   Result similar to that returned by stride(pdb)$sse and dssp(pdb)$sse
  ##   This could be incorporated into read.pdb() if found to be more generally useful
  ## Modified to include helix "type" and sheet "sense"
  ## - Jan 15, 2015

  if(is.null(pdb$helix) & is.null(pdb$sheet)) {
    warning("No helix and sheet defined in input 'sse' PDB object: try using dssp()")
    ##ss <- try(dssp(pdb)$sse)
    ## Probably best to get user to do this separately due to possible 'exefile' problems etc..
    return(NULL)
  }
  rn <- pdb$atom[pdb$calpha, c("resno", "chain")]
  ss <- rep(" ", nrow(rn))
  names(ss) = paste(rn$resno,rn$chain, "", sep="_")

  for(i in 1:length(pdb$helix$start)) {
    ind <- (rn$chain==pdb$helix$chain[i] &
         rn$resno >= pdb$helix$start[i] &
         rn$resno <= pdb$helix$end[i])
    ss[ind] = "H"
    names(ss)[ind] = sub("_$", paste("_", pdb$helix$type[i], sep=""), 
                     names(ss)[ind])
  }
  for(i in 1:length(pdb$sheet$start)) {
    ind <- (rn$chain==pdb$sheet$chain[i] &
         rn$resno >= pdb$sheet$start[i] &
         rn$resno <= pdb$sheet$end[i])
    ss[ind] = "E"
    names(ss)[ind] = sub("_$", paste("_", pdb$sheet$sense[i], sep=""), 
                     names(ss)[ind])
  }
  return(ss)
}

.bounds.sse <- function(sse) {
  ## - reverse operation of .pdb2sse()
  string0 <- strsplit(names(sse), split="_")
  chain <- sapply(string0, "[", 2)
  resno <- as.numeric(sapply(string0, "[", 1))
  type <- sapply(string0, function(x) ifelse(length(x)>2, x[3], NA))

  string <- paste(sse, chain, sep="_")

  ends <- c(which(string[-length(string)] != string[-1]), length(string))
  starts <- c(1, ends[-length(ends)] + 1)
  inds <- cbind(starts, ends)
  colnames(inds) <- c("start", "end")

  h <- which(sse[inds[, "start"]] == "H")
  h.inds <- list(start = resno[inds[h, "start"]], end = resno[inds[h, "end"]], 
     chain = chain[inds[h, "start"]], type = type[inds[h, "start"]])
  e <- which(sse[inds[, "start"]] == "E")
  e.inds <- list(start = resno[inds[e, "start"]], end = resno[inds[e, "end"]], 
     chain = chain[inds[e, "start"]], sense = type[inds[e, "start"]])

  return(list(helix = h.inds, sheet = e.inds) )
}
