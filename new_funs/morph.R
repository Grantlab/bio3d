morph <- function(x, y, step=0.125, file='morph.pdb', rock=TRUE) {
   if(!is.pdb(x) || !is.pdb(y)) {
      stop('"x" and "y" must be a pdb object')
   }

   if(sum(x$calpha) != sum(y$calpha)) {
      stop('"x" and "y" must have the same length')
   }

   x <- trim(x, 'calpha'); y <- trim(y, 'calpha')

   ## - make a dummy 'pca' object
   L <- sum((y$xyz - x$xyz)^2)
   mean <- as.vector((y$xyz + x$xyz) / 2)
   U <- matrix((y$xyz - x$xyz) / sqrt(L), ncol=1)
   out <- list(L=L, U=U, mean=mean)
   class(out) <- "pca"

   mktrj(out, pc=1, mag=0.5, step=step, pdb=x, file=file, rock=rock)
}
