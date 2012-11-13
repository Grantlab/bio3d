`dm.xyz` <-
function(xyz, grpby=NULL, scut=NULL, mask.lower=TRUE) {
  ##-- New distance matrix function with 'grpby' option
  ##  ndm(pdb$xyz, grpby=pdb$atom[,"resno"], scut=3)

  if( !is.vector(xyz) )
    stop(" Input 'xyz' should be a numeric vector" )
    
  ##- Full Distance matrix (could use 'dm' or 'dist.xyz')
  d <- as.matrix(dist(matrix(xyz, ncol = 3, byrow = TRUE)))
  
  ##- Mask lower.tri  
  if( mask.lower )
    d[lower.tri(d)] = NA

  ##- Mask concetive atoms
  if( is.null(grpby) ) {
    if (!is.null(scut))
      d[diag.ind(d, n = scut)] = NA

    return(d)
    
  } else {

    ##- Group by concetive numbers in 'grpby'
    if( length(xyz) != (length(grpby)*3) )
      stop("dimension miss-match in 'xyz' and 'grpby', check lengths")

    ##- Bounds of 'grpby' numbers
    inds <- bounds(grpby, dup=TRUE)
    nres <- nrow(inds)
    
    ##- Per-residue matrix
    m <- matrix(, ncol=nres, nrow=nres)
    ij <- pairwise(nres)
    
    ##  Ignore concetive groups (residues)
    if (!is.null(scut))
      ij <- ij[ij[,2]-ij[,1] > (scut-1),]
    
    ##- Min per residue
    for(k in 1 : nrow(ij) ) {
      m[ij[k,1],ij[k,2]] <-
        min( d[ (inds[ij[k,1],"start"]:inds[ij[k,1],"end"]),
                  (inds[ij[k,2],"start"]:inds[ij[k,2],"end"])],
            na.rm=TRUE )
    }
    if( !mask.lower )
      m[lower.tri(m)] = t(m)[lower.tri(m)]
    
    return(m)
  
  }
}

