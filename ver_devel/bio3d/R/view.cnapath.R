view.cnapath <- function(x, pdb, out.prefix = "view.cnapath", launch = FALSE) {

   if(!inherits(x, "cnapath")) 
      stop("Input x is not a 'cnapath' object")

   file = paste(out.prefix, ".vmd", sep="")
   pdbfile = paste(out.prefix, ".pdb", sep="")

   res <- unique(unlist(x$path))
   ind.source <- match(x$path[[1]][1], res)
   ind.sink <- match(x$path[[1]][length(x$path[[1]])], res)
   
   rmin <- min(x$dist)
   rmax <- max(x$dist)
   rad <- function(r, rmin, rmax, radmin = 0.01, radmax = 0.5) {
      (rmax - r) / (rmax - rmin) * (radmax - radmin) + radmin
   }

   cols <- colorRampPalette(c("blue", "red"))(256)
   col.mat <- matrix(NA, length(res), length(res))
   conn <- matrix(0, length(res), length(res))
   rr <- conn
   for(j in 1:length(x$path)) {
      y = x$path[[j]]
      for(i in 1:(length(y)-1)) {
         i1 = match(y[i], res)
         i2 = match(y[i+1], res)
         if(conn[i1, i2] == 0) conn[i1, i2] = conn[i2, i1] = 1
         r = rad(x$dist[j], rmin, rmax)
         ic = floor((rmax - x$dist[j]) / (rmax - rmin) * 255) + 1
         col = cols[ic]
         if(r > rr[i1, i2]) {
            rr[i1, i2] = rr[i2, i1] = r
            col.mat[i1, i2] = col.mat[i2, i1] = col
         }
      }
   }
 
   ca.inds <- atom.select(pdb, elety="CA", verbose = FALSE)
   res.pdb <- pdb$atom[ca.inds$atom[res], "resno"] 
   
   rownames(conn) <- res
   colnames(conn) <- res
   rownames(rr) <- res
   colnames(rr) <- res

   # Draw molecular structures
   cat("mol new ", pdbfile, " type pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all
mol delrep 0 top
mol representation NewCartoon 0.300000 10.000000 4.100000 0
mol color colorID 8
mol selection {all}
mol material Opaque
mol addrep top
mol representation Licorice 0.300000 10.000000 10.000000
mol color name
mol selection {(resid ", res.pdb[ind.source], " ", res.pdb[ind.sink], ")}
mol material Opaque
mol addrep top 
mol representation VDW 0.4 10
mol color colorID 2
mol selection {(resid ", paste(res.pdb, collapse=' '), ") and name CA}
mol material Opaque
mol addrep top
", file=file)

   # Draw paths 
   k = 500
   for(i in 1:(nrow(conn)-1)) {
      for(j in (i+1):ncol(conn)) {
         if(conn[i, j] == 1) {
            col = as.numeric(col2rgb(col.mat[i, j]))/255
            cat("color change rgb ", k, " ", paste(col, collapse=" "), "\n", sep="", file=file, append=TRUE)
            cat("graphics top color ", k, "\n", sep="", file=file, append=TRUE)
            cat("draw cylinder {", pdb$xyz[atom2xyz(ca.inds$atom[res[i]])], 
               "} {", pdb$xyz[atom2xyz(ca.inds$atom[res[j]])], "} radius", rr[i, j], 
               " resolution 6 filled 0\n", sep=" ", file=file, append=TRUE)
            k = k + 1
         }
      }
   } 

   write.pdb(pdb, file=pdbfile)
   
   if(launch) {
      
      cmd <- paste("vmd -e", file)

      os1 <- .Platform$OS.type
      if (os1 == "windows") {
        shell(shQuote(cmd))
      } else{
        if(Sys.info()["sysname"]=="Darwin") {
          system(paste("/Applications/VMD\\ 1.9.*app/Contents/MacOS/startup.command -e", file))
        } else {
          system(cmd)
        }
      }
   }
}
