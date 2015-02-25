# A function to automate the generation
# of requested multimeric alignment from a single chain based 'pdbs'
play.biounit <- function(pdbs, nmer = 1, dir.rawpdb = "raw_pdbs", 
      dir.output = "biounit", cofactor.ext = "_cofactor", file = "aln.fa", ncore=NULL) {
   require(bio3d)
   source("pdb.biounit.R")

   # raw pdb ids to process
   ids = file.path(dir.rawpdb, paste(unique(substr(basename(pdbs$id), 1, 4)), ".pdb", sep=""))
   cat("Found ", length(ids), " raw pdb files\n")

   # loop to process each id
   files <- lapply(ids, function(x) {
      cat("\n     Generating biological unit for ", x, "...\n")

      biounit <- pdb.biounit(x)
      if(is.pdb(biounit) ) {
         cat("          No valid biological unit was found\n")
         return(NA)
      }
      cat("          Found ", length(biounit), " biological units\n")

      invisible( capture.output( pdb <- read.pdb(x) ))
      chain <- unique(pdb$atom[, "chain"])

      # check if the biounit is what we want 
      new.biounit <- lapply(1:length(biounit), function(i) {
          cat("               checking biounit No. ", i, "...\n")
          y <- biounit[[i]]

          # Is it correct multimeric state?
          new.chain <- unique(y$atom[, "chain"])
          cat("               The biounit contains ", length(new.chain), " chains\n")
          # temporarily split chains
          tpdbs <- lapply(new.chain, function(x) trim.pdb(y, chain=x, verbose=FALSE))
          # which chain is our target?
          chk.id <- paste(substr(basename(x), 1, 4), "_", new.chain, sep="")
          ref.inds <- match(chk.id, basename.pdb(pdbs$id))
          # with the reference chain id for target obtained above, 
          # find more chains that are identical copies of the target
          match.ind <- sapply(tpdbs, function(x) {
             aa <- as.character(pdbseq(x))
             for(j in ref.inds) {
                aa1 <- pdbs$ali[j, !is.gap(pdbs$ali[j, ])]
                if(identical(aa, aa1)) return (TRUE)
             }
             FALSE
          } )
          # multimeric state of current biounit
          cur.nmer <- sum(match.ind)
          cat("               ", cur.nmer, " chains are the target (request is ", nmer, " chains)\n")
          target.ind <- which(match.ind)
          if(cur.nmer>0) cat("                  ", paste(chk.id[target.ind], collapse=" "), "\n")

          if(cur.nmer == nmer) {
            # place target chains the first
            if(any(!match.ind[1:nmer])) {
               cat("               Rearrange chain order...\n")
               tpdbs <- c(tpdbs[target.ind], tpdbs[-target.ind])
               cl <- y$call
               y <- do.call(cat.pdb, tpdbs)
               y$call <- cl
            }
          }
          else
            y <- NA
          y 
      } )
      new.biounit <- new.biounit[!is.na(new.biounit)]

      if(length(new.biounit) == 0) {
         cat("          No valid biological unit was found\n")
         return(NA)
      }
     
      # check if duplicate hits exist (duplication may occur in the case that 
      # biounits are distinguished by cofactor binding state)
      if(length(new.biounit) > 1) {
         keep.ind <- 1:length(new.biounit)
         np <- pairwise(length(new.biounit))
         for(i in 1:nrow(np)) {
            bio1 <- new.biounit[[np[i, 1]]]
            bio2 <- new.biounit[[np[i, 2]]]
            chs1 <- unique(bio1$atom[, "chain"])
            chs2 <- unique(bio2$atom[, "chain"])
            dup <- all( sapply(chs1, function(x) 
              any( sapply(chs2, function(y) 
                     identical(bio1$atom[bio1$atom[, "chain"] == x, ], 
                               bio2$atom[bio2$atom[, "chain"] == y, ])
                   ) )
            ) )
           if(dup) {
#              cat("In ", x, ", biounit ", np[i, 1], " (", length(chs1), 
#                 " chains) has the same coords to biounit ", np[i, 2], " (", length(chs2), " chains)\n")
#              cat("Keep the one has the largest number of chains\n")
              cat("          Removing duplicated biological units (keep the one containing most chains)...\n")
              keep.ind <- keep.ind[-np[i, which.min(length(chs1), length(chs2))]]
           }
         }
         new.biounit <- new.biounit[keep.ind]
      }

      cat("          Writing ", length(new.biounit), " biological units to ", dir.output, "\n")
      
      if(!file.exists(dir.output)) dir.create(dir.output)

      files <- sapply( new.biounit, function(y) 
         file.path(dir.output, paste(substr(basename(x), 1, 4), "_", y$atom[1, "chain"], ".pdb", sep="")) )

      # if file names are conflicting, simply add sequencial numbers to distinguish them
      nn <- bounds(files, dup.inds=TRUE, pre.sort=FALSE)
      for(i in 1:nrow(nn)) {
         if(nn[i, "length"] > 1)
            files[nn[i, "start"]:nn[i, "end"]] <- 
               paste(files[nn[i, "start"]:nn[i, "end"]], 1:nn[i, "length"], ".pdb", sep="")
      }

      # write separated biounit and cofactor pdb files 
      sapply(1:length(new.biounit), function(i) {
         chs <- unique(new.biounit[[i]]$atom[, "chain"])
         opdb <- trim.pdb(new.biounit[[i]], chain=chs[1:nmer], verbose=FALSE)
         write.pdb(opdb, file=files[i])
         if(length(chs) > nmer) {
            co.pdb <- trim.pdb(new.biounit[[i]], chain=chs[(nmer+1):length(chs)], verbose=FALSE)
            write.pdb(co.pdb, file=file.path(dir.output, paste(basename.pdb(files[i]), cofactor.ext, ".pdb", sep="")))
         }
      } )

      files
   } )
   files <- files[!is.na(files)]
 
   # Continue to do alignment
   cat("\nGenerating alignment for biological units...\n") 
   new.pdbs <- pdbaln(unlist(files), ncore = ncore)

   if(!is.null(file)) write.fasta(new.pdbs, file = file)
   
   return(new.pdbs)
}
