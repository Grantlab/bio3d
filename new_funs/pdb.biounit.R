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
#' \donttest{
#'    pdb <- read.pdb("2dn1")
#'    biounit <- pdb.biounit("2dn1")
#'    pdb
#'    biounit
#' }
#' \dontrun{
#'    biounit <- pdb.biounit("2bfu")
#'    write.pdb(biounit[[1]], file="biounit.pdb")
#'    # open the pdb file in VMD to have a look on the biological unit
#' } 
pdb.biounit <- function(file, ncore = NULL, ...) {

    require(bio3d)
    require(parallel)   
    ncore = setup.ncore(ncore) 

    cl <- match.call()

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

       biounits <- lapply(1:remarks$num, function(i) {
          # the transformation matrices
          mats <- remarks$mat[[i]]

          # applied to the chains
          chain <- remarks$chain[[i]]
 
          # number of copies
          ncopy <- length(mats)
  
          # Are chains treated differently?
          nn <- length(unique(names(mats)))

          if(ncopy <= 10 || nn > 1) {    
             ## save copies as individual chains 

             # The original copy stored as spearated chains
             biounit0 <- lapply(chain, function(x) trim.pdb(pdb, chain=x, verbose=FALSE) )
             # available chain ID repository
             chains0 <- setdiff(c(LETTERS, letters, 0:9), chain)
             
             jch <- 1
             used.chain <- NULL
             biounit <- NULL
             for(j in 1:ncopy) {
                mt <- mats[[j]]
                chs <- strsplit(names(mats)[j], split=" ")[[1]]
                for(k in chs) {
                   bio <- biounit0[[match(k, chain)]]
                   xyz <- rbind(matrix(bio$xyz, nrow=3), 1)
                   xyz <- matrix(mt %*% xyz, nrow = 1)
                   if(! k %in% used.chain) {
                      ch <- k
                      used.chain <- c(used.chain, k)
                   } else {
                      ch <- chains0[jch]
                      jch = jch + 1
                   }
                   bio$xyz <- xyz
                   bio$atom[, "chain"] <- ch
                   bio$atom[, c("x", "y", "z")] <- round(matrix(xyz, ncol=3, byrow=TRUE), digits=3)
                   
                   biounit <- c(biounit, list(bio))
                }
             } 
             biounit <- do.call(cat.pdb, biounit)

#             # temporarily write the pdb of biounit and re-read it
#             tmpf <- tempfile()
#             write.pdb(biounit, file=tmpf)
#             biounit = read.pdb(tmpf, verbose=FALSE)
          } 
          else {
             ## save copies as multi-models

             # The original copy
             biounit <- trim.pdb(pdb, chain=chain, verbose=FALSE)

             xyz = rbind(matrix(biounit$xyz, nrow=3), 1)
             ll <- mclapply(2:ncopy, function(j) {

                mt <- mats[[j]]
                xyz = matrix(mt %*% xyz, nrow=1)
                xyz
             }, mc.cores = ncore )
             biounit$xyz <- rbind(biounit$xyz, do.call(rbind, ll))
             class(biounit$xyz) <- "xyz"
          }
          
          biounit$call <- cl
          return(biounit)
       } ) # end of lapply(1:remarks$num)

       ## multimeric state
       nchs <- sapply(biounits, function(x) length(unique(x$atom[, "chain"])) * nrow(x$xyz))
       mer <- c("monomer", "dimer", "trimer", "tetramer", "multimer")
       names(biounits) <- paste(remarks$method, ".determined.", 
            mer[ifelse(nchs>5, 5, nchs)], " (",  nchs, " chains)", sep="") 
#       if(length(biounits) == 1) biounits = biounits[[1]]  
    }
    else {
      warning("Problems occurred in building biological unit. Return 'pdb' as is in the file")
      biounits <- pdb
    }

    return(biounits)
}

.parse.pdb.remarks <- function(x) {

    raw.lines <- x

    # How many lines of REMARK 350?
    remark350 <- grep("^REMARK\\s+350", raw.lines)
    nlines <- length(remark350)

    # How many distinct biological unit?
    biolines <- grep("^REMARK\\s+350\\s+BIOMOLECULE", raw.lines)
    nbios <- length(biolines)
    
    if(nbios == 0) {
       warning("REMARK 350 is incomplete.")
       return(NULL)
    }
    
    # End line number of each biological unit
    biolines2 <- c(biolines[-1], remark350[nlines])

    # How the biological unit was determined?
    method <- sapply(1:nbios, function(i) {
       author <- intersect(grep("^REMARK\\s+350\\s+AUTHOR DETERMINED BIOLOGICAL UNIT", raw.lines), 
                            biolines[i]:biolines2[i])
       if(length(author) >= 1) return("AUTHOR")
       else return("SOFTWARE")
    } )
    
    # Get chain IDs to apply the transformation
    chain <- lapply(1:nbios, function(i) {
       chlines <- intersect(grep("^REMARK\\s+350\\s+APPLY THE FOLLOWING TO CHAINS", raw.lines), 
                            biolines[i]:biolines2[i])
       if(length(chlines) >= 1) {
          chs <- gsub("\\s*", "", sub("^.*:", "", raw.lines[chlines]))
          chs <- unlist(strsplit(chs, split=","))
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
       # Get chain ID again: different trans matrices may be applied to different chains
       chlines <- intersect(grep("^REMARK\\s+350\\s+APPLY THE FOLLOWING TO CHAINS", raw.lines), 
                            biolines[i]:biolines2[i])
       chs <- gsub("\\s*", "", sub("^.*:", "", raw.lines[chlines]))
       chs <- strsplit(chs, split=",")
        
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
          chs.pos <- findInterval(mtlines[seq(1, length(mtlines), 3)], chlines)
          names(mat) <- sapply(chs[chs.pos], paste, collapse=" ") ## apply each mat to specific chains
       }
       return(mat)
    } )
    if(any(is.na(mat))) return(NULL)
  
    out <- list(num=nbios, chain=chain, mat=mat, method=method) 
    return(out)
}
