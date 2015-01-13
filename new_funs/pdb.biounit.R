#' Construct Biological Units Based on A PDB File
#'
#' Construct biological assemblies/units based on a local 
#' PDB file, a 4-letter PDB code, or a url to a PDB file online.
#'
#' @details
#' A valid structural/simulation study should be performed on the biological 
#' unit of a protein system. For example, the alpha2-beta2 tetramer form of
#' hemoglobin. However, canonical PDB files usually contain the asymmetric unit of 
#' the crystal cell, which can be:
#' \enumerate{
#'      \item One biological unit 
#'      \item A portion of a biological unit 
#'      \item Multiple biological units
#' }
#' The function performs symmetry operations to the coordinates based on the 
#' transformation matrices stored in the PDB file, and returns the 
#' \code{pdb} object of biological unit for later analysis.
#'
#' @param file a single element character vector containing the name of the
#'          PDB file to be read, or the four letter PDB identifier for
#'          online file access. 
#' @param ncore number of CPU cores used to do the calculation. By default
#'          (\code{ncore=NULL}), use all available CPU cores.
#' @param ... additional arguments passed to \code{\link{read.pdb}}.
#'
#' @return 
#'    a \code{pdb} object if one biological unit is constructed, or
#'    a list of \code{pdb} objects with each representing an individual biological unit.
#'
#'    If the number of generated copies is small (<=10), each copy is represented as
#'    independent chains with distinct chain IDs.
#'    If the number of copies is large (>10), copies are represented as multiple 'models'.
#'
#' @seealso \code{\link{read.pdb}}
#'
#' @author Xin-Qiu Yao
#'
#' @examples
#' pdb <- read.pdb("2dn1")
#' biounit <- pdb.biounit("2dn1")
#' pdb
#' biounit
#'
#' \dontrun{
#' pdb <- pdb.biounit("2bfu")
#' write.pdb(pdb, file="biounit.pdb")
#' # open the pdb file in VMD to have a look on the biological unit
#' } 
pdb.biounit <- function(file, ncore = NULL, ...) {

    require(bio3d)
    require(parallel)   
    ncore = setup.ncore(ncore) 

    # This also does initial check on 'file' 
    pdb <- read.pdb(file, ...)

    dots <- list(...)
    if("maxlines" %in% names(dots)) 
       maxlines = dots$maxlines
    else 
       maxlines = eval(formals(read.pdb)$maxlines)

    ## same check as in read.pdb()
    toread <- file.exists(file)
    if (substr(file, 1, 4) == "http") {
        toread <- TRUE
    }
    if (!toread) {
        if (nchar(file) == 4) {
            file <- get.pdb(file, URLonly = TRUE)
        }
    }
    #####

    raw.lines <- readLines(file, n = maxlines)
    remarks <- .parse.pdb.remarks(raw.lines)
 
    if(!is.null(remarks)) { 
       # number of copies in each biological unit
       ncopy <- sapply(remarks$mat, length)
   
       if(max(ncopy) <= 10) {    
          ## save copies as individual chains 

          # available chain ID repository
          chains0 <- setdiff(c(LETTERS, letters, 0:9), unique(pdb$atom[, "chain"]))
       
          # number of chains to apply transformation in each biological unit
          nchain <- sapply(remarks$chain, length)

          # start index in the chain ID repository for each biological unit
          ch.ind <- cumsum(c(1, (ncopy - 1) * nchain))

          biounits <- lapply(1:remarks$num, function(i) {

             # The original copy
             biounit <- trim.pdb(pdb, chain=remarks$chain[[i]])

             if(ncopy[i] >= 2) {    
                xyz = rbind(matrix(biounit$xyz, nrow=3), 1)
                chain <- biounit$atom[, "chain"] 
                rl <- rle(chain)

                ll <- mclapply(2:ncopy[i], function(j) {

                   mt <- remarks$mat[[i]][[j]]
                   xyz = matrix(mt %*% xyz, nrow=1)

                   ch <- chains0[(ch.ind[i] + nchain[i] * (j-2)) :
                              (ch.ind[i] + nchain[i] * (j-1) - 1)]
                   names(ch) <- remarks$chain[[i]]
                   rl$values <- ch[rl$values]
                   chain <- inverse.rle(rl)
                   return(list(xyz=xyz, chain=chain))
                }, mc.cores = ncore )

                biounit1 <- biounit
                for(j in 2:ncopy[i]) {
                   biounit1$atom <- rbind(biounit1$atom, biounit$atom)
                }
                xyz <- lapply(ll, "[[", "xyz")
                xyz <- cbind(biounit$xyz, do.call(cbind, xyz))
                chain <- c(chain, sapply(ll, "[[", "chain"))
            
                # temporarily write the pdb of biounit, change ATOM (chain id, xyz, etc), and read it again
                tmpf <- tempfile()
                write.pdb(biounit1, xyz = xyz, chain = chain, file=tmpf)
                biounit = read.pdb(tmpf, verbose=FALSE)
             }
             return(biounit)
          } )
       }
       else {
          ## save copies as multi-models
          biounits <- lapply(1:remarks$num, function(i) {

             # The original copy
             biounit <- trim.pdb(pdb, chain=remarks$chain[[i]])

             if(ncopy[i] >= 2) {    
                xyz = rbind(matrix(biounit$xyz, nrow=3), 1)
                ll <- mclapply(2:ncopy[i], function(j) {

                   mt <- remarks$mat[[i]][[j]]
                   xyz = matrix(mt %*% xyz, nrow=1)
                   xyz
                }, mc.cores = ncore )
                biounit$xyz <- rbind(biounit$xyz, do.call(rbind, ll))
                class(biounit$xyz) <- "xyz"
             }
             return(biounit)
          } )
       }
       if(length(biounits) == 1) biounits = biounits[[1]]  
    }
    else {
      warning("Problems occurred in building biological unit. Return 'pdb' as is in the file")
      biounits <- pdb
    }

    return(biounits)
}

.parse.pdb.remarks <- function(x) {

    raw.lines <- x

    # How many copies of biological unit?
    biolines <- grep("^REMARK\\s+350\\s+BIOMOLECULE", raw.lines)
    nbios <- length(biolines)

    if(nbios == 0) {
       warning("REMARK 350 is incomplete.")
       return(NULL)
    }
    
    # End line number of each biological unit
    biolines2 <- c(biolines[-1], length(raw.lines))

    chain <- lapply(1:nbios, function(i) {
       # Get chain IDs to apply the transformation
       chlines <- intersect(grep("^REMARK\\s+350\\s+APPLY THE FOLLOWING TO CHAINS", raw.lines), 
                            biolines[i]:biolines2[i])
       if(length(chlines) == 1) {
          chs <- gsub("\\s*", "", sub("^.*:", "", raw.lines[chlines]))
          chs <- strsplit(chs, split=",")[[1]]
       } 
       else {
          warning(paste("Can't determine chain IDs from REMARK 350 for biological unit", 
               i, sep=""))
          chs = NA
       }
       return(chs)
    } )
    if(any(is.na(chain))) return(NULL)

    mat <- lapply(1:nbios, function(i) {
       # Get transformation matrices
       mtlines <- intersect(grep("^REMARK\\s+350\\s+BIOMT", raw.lines), 
                            biolines[i]:biolines2[i])
        
       if(length(mtlines) == 0 || length(mtlines) %% 3 != 0) {
          warning("Incomplete transformation matrix")
          mat <- NA
       }
       else {
          mat <- lapply(seq(1, length(mtlines), 3), function(j) {
             mt <- matrix(NA, 3, 4)
             for(k in 1:3) {
                vals <- sub("^REMARK\\s+350\\s+BIOMT[123]\\s*", "", raw.lines[mtlines[j+k-1]])
                vals <- strsplit(vals, split="\\s+")[[1]]
                mt[k, ] <- as.numeric(vals[-1])
             }
             mt
          } )
       }
       return(mat)
    } )
    if(any(is.na(mat))) return(NULL)
  
    out <- list(num=nbios, chain=chain, mat=mat) 
    return(out)
}
