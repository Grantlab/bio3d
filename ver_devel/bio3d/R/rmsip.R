"rmsip" <-
  function(modes.a, modes.b, subset = 10, row.name="a", col.name="b") {
    
    if(missing(modes.a))
      stop("rmsip: 'modes.a' must be prodivded")
    if(missing(modes.b))
      stop("rmsip: 'modes.b' must be prodivded")
    
    if ("pca" %in% class(modes.a)) {
      ev.a <- modes.a$U
      first.mode <- 1
      inds.a <- seq(first.mode, subset)
    }
    else if("nma" %in% class(modes.a)) {
      ev.a <- modes.a$U ## Use orthogonal raw vectors
      #ev.a <- modes.a$modes ## Otherwise modes needs to be mass-w re-normed
      mass.a <- modes.a$mass
      first.mode <- modes.a$triv.modes+1
      nmodes <- modes.a$triv.modes + subset
      inds.a <- seq(first.mode, nmodes)
    }
    else {
      if(class(modes.a)!="matrix" && class(modes.a)!="pca.loadings")
        stop("overlap: 'modes.a' must be an object of type 'pca', 'nma', or 'matrix'")
      ev.a <- modes.a
      first.mode <- 1
      inds.a <- seq(first.mode, subset)
    }

    if ("pca" %in% class(modes.b)) {
      ev.b <- modes.b$U
      first.mode <- 1
      inds.b <- seq(first.mode, subset)
    }
    else if("nma" %in% class(modes.b)) {
      ev.b <- modes.b$U
      #ev.b <- modes.b$modes
      mass.b <- modes.b$mass
      first.mode <- modes.b$triv.modes+1
      nmodes <- modes.b$triv.modes + subset
      inds.b <- seq(first.mode, nmodes)
    }
    else {
      if(class(modes.b)!="matrix" && class(modes.b)!="pca.loadings")
        stop("overlap: 'modes.b' must be an object of type 'pca', 'nma', or 'matrix'")
      ev.b <- modes.b
      first.mode <- 1
      inds.b <- seq(first.mode, subset)
    }

    mass.a <- NULL; mass.b <- NULL;
    x <- normalize.vector(ev.a[,inds.a], mass.a)
    y <- normalize.vector(ev.b[,inds.b], mass.b)

    if(is.null(mass.a))
      o <- ( t(x) %*% y ) **2
    else
      o <- t(apply(x, 2, inner.prod, y, mass.a) **2)

    if (!is.null(row.name)) {
      rownames(o) <- paste(row.name, c(1:subset), sep="")
    }
    
    if (!is.null(col.name)) {
      colnames(o) <- paste(col.name, c(1:subset), sep="")
    }
    
    rmsip <- sqrt(sum(o)/subset)
    out <- list(overlap=round(o,3), rmsip=rmsip)

    class(out) <- "rmsip"
    return( out )
  }
