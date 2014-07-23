"print.pca" <-
  function(x, nmodes=6, ...) {

    cn <- class(x)
    
    cat("\nCall:\n  ", paste(deparse(x$call), sep = "\n", collapse = "\n"),
        "\n\n", sep = "")
    
    cat("Class:\n  ", cn, "\n\n", sep = "")
    
    cat("Number of eigenvalues:\n  ", length(x$L), 
        "\n\n", sep="")

    inds <- 1:nmodes
    eign <- round(x$L[inds], 3)
    cat("Eigenvalues:\n", sep="")
    
    i <- 0
    for ( f in eign ) {
      i <- i+1
      cat("  PC ", i, ": \t", f,  "\n", sep="")
    }
    
    cat("\n")

  i <- paste( attributes(x)$names, collapse=", ")
  cat(strwrap(paste(" + attr:",i,"\n"),width=60, exdent=8), sep="\n")

    invisible(x)
  }
