sse.vector <- function(pdb, sse=pdb, coil="-", helix="H", sheet="E", calpha=TRUE) {

	##-- Output a secondary structure vector for a given pdb input.
	##   Can be used to output sse sequence or color vector for visualisation.
	##    E.G.
	##     paste( sse.vector(pdb), collapse="")  ## "-----EEEEEEEE----HHHHHH-
	##     sse.vector(pdb, coil="gray", helix="purple", sheet="yellow")
	##

	if(!is.pdb(pdb))
    stop("Input 'pdb' must be a 'pdb' class object as returned by the read.pdb function")

	rn.ch <- function(x=sse$helix) {
		## Return SSE resno_chain reference vector
		if(length(x$start) > 0) {
			nres <- (x$end-x$start+1)
			ch <- rep(x$chain, nres)
			rn <- unbound(x$start, x$end)
			rn.ch <- paste(rn, ch, sep="_")
			attributes(rn.ch) <- list(nrep=nres)
		} else {
			rn.ch <- NULL
		}
		return(rn.ch)
	}

	ref.h <- rn.ch(sse$helix)
	ref.e <- rn.ch(sse$sheet)

	## Note. in the future attributes(ref.h)$nrep can be used to 
	## replicate $type, etc. E.G. 
	##   "rep(x$type, attributes(ref.h)$nrep)""

	if(is.null(ref.h)) {
		warning("No alpha helices defined in $helix 'sse' input")
	}
	if(is.null(ref.e)) {
		warning("No beta strands defined in $sheet 'sse' input")
	}

	## All residue resno_chain reference vector
	if(calpha) {
		ca.inds <- atom.select(pdb, "calpha", verbose=FALSE)
  	ref <- paste(pdb$atom[ca.inds$atom, "resno"],
    	           pdb$atom[ca.inds$atom, "chain"], sep="_")
  } else {
  	ref <- paste(pdb$atom$resno, pdb$atom$chain, sep="_")
  }
	ss <- rep(coil, length(ref))
	names(ss) = ref
	ss[(ref %in% ref.h)] = helix
	ss[(ref %in% ref.e)] = sheet
	return(ss)
}

