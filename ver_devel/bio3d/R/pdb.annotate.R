"pdb.annotate" <- function(ids, anno.terms=NULL, unique=FALSE, verbose=FALSE, 
                           extra.terms=NULL) {
  
  oops <- !requireNamespace("httr", quietly = TRUE)
  if(oops) {
    stop("Please install the httr package from CRAN")
  }
  
  if(!is.null(extra.terms)) {
    message("Currently 'extra.terms' is not supported")
    extra.terms <- NULL
  }
  
  if(inherits(ids, "blast")) ids = ids$pdb.id
  
  if(!is.vector(ids)) {
    stop("Input argument 'ids' should be a vector of PDB identifiers/accession codes")
  }
  
  ## Basic annotation terms (note 'citation' is a meta term)
  anno.basicterms <- c("structureId", "chainId", "macromoleculeType",
                       "chainLength", "experimentalTechnique",  "resolution", 
                       "scopDomain", "pfam", "ligandId", "ligandName",
                       "source", "structureTitle", "citation", "rObserved",
                       "rFree", "rWork", "spaceGroup")  
  if(is.null(anno.terms)) {
    anno.terms <- anno.basicterms
  } 
  else {
    anno.terms <- match.arg(anno.terms, anno.basicterms, several.ok=TRUE)
    anno.terms <- unique(anno.terms)
  }
  anno.terms.input <- anno.terms
  
  ## Check if we have any valid terms remaining
  if( length(anno.terms) == 0 ) {
    stop( paste("No valid anno.terms specified. Please select from:\n\t ",
                paste(anno.basicterms, collapse=", ")) )
  }
  
  ## force the structureId and chainId terms to be present
  req.terms <- c("structureId", "chainId")
  
  if(any(c("ligandId", "ligandName") %in% anno.terms)) {
    ## force ligandChainId
    req.terms <- c(req.terms, "ligandChainId")
  }
  
  inds <- req.terms %in% anno.terms

  if(!all(inds)) {
    anno.terms <- c(req.terms[!inds], anno.terms)
  }
  
  if (missing(ids)) {
    stop("please specify PDB ids for annotating")
  }
  
  if (any(nchar(ids) != 4)) {
#    warning("ids should be standard 4 character PDB-IDs: trying first 4 characters...")

#    if(unique) {
#      ids <- unique(substr(basename(ids), 1, 4))
#    }
    
    ## first 4 chars should be upper
    ## any chainId should remain untouched - see e.g. PDB ID 3R1C
    mysplit <- function(x) {
      str <- unlist(strsplit(x, "_"))
      if(length(str)>1) {
        paste(toupper(str[1]), "_", str[2], sep="")
      }
      else {
        toupper(str[1])
      }
    }
  
    ids <- unlist(lapply(ids, mysplit))
  } else {
    ids <- toupper(ids)
  }
  ids.short <- unique( substr(basename(ids), 1, 4) )

  ## prepare query
  baseurl <- "https://data.rcsb.org/graphql"

  anno.terms.new <- unlist(sapply(anno.terms, .map_terms))
  anno.terms.new <- c(anno.terms.new, extra.terms)
  
  query <- "query($id: [String!]!){
     entries(entry_ids: $id){
  "
  query <- paste(query, 
    paste( sapply(anno.terms.new, .string2json), collapse="\n"), 
    "\n}}", sep="")

  resp <- httr::POST(baseurl,
                     httr::accept_json(),
                     body = list(query=query,
                                 variables=list(id=ids.short)), 
                     encode="json")
  
  if(httr::http_error(resp)) {
    stop('Access to PDB server failed')
  }
  else {
    ret <- httr::content(resp)
  }
  
  if("error" %in% names(ret)) {
    stop('Retrieving data from PDB failed')
  }
  
  ## Generate a formatted table
  ## Also taking care of merging data for unique structureId, 
  ##   excluding non-requested chain IDs, formatting citation, etc.
  data <- .format_tbl(ret, ids, anno.terms, unique=unique)
  
  if(unique) {
    rownames(data) <- data$structureId
  } 
  else {
    rownames(data) <- paste(data$structureId, data$chainId, sep="_")
  }
  
  ## include only requested terms 
  ## (NOTE: need to modify for future support of 'extra.terms')
  col.inds <- which(colnames(data) %in% anno.terms.input)
  data <- data[, col.inds, drop=FALSE]

  return(data)
}

## map a string to JSON-like input parameters
.string2json <- function(x) {
  x <- strsplit(x, split="\\.")[[1]]
  if(length(x)>1) {
    paste( paste(x, collapse="{"), paste(rep("}", length(x)-1), collapse=""), sep="")
  }
  else if(length(x)==0) {
    ""
  }
  else {
    x
  }
}

## map from old terms to the new ones
.map_terms <- function(x) {
  switch(x,
         "structureId" = "entry.id",
         "chainId" = "polymer_entities.polymer_entity_instances.rcsb_polymer_entity_instance_container_identifiers.auth_asym_id",
         "macromoleculeType" = "polymer_entities.polymer_entity_instances.polymer_entity.entity_poly.rcsb_entity_polymer_type",
         "chainLength" = "polymer_entities.polymer_entity_instances.polymer_entity.entity_poly.rcsb_sample_sequence_length",
         "experimentalTechnique" = "rcsb_entry_info.experimental_method",
         "resolution" = "rcsb_entry_info.resolution_combined",
         "scopDomain" = paste(
           "polymer_entities.polymer_entity_instances.rcsb_polymer_instance_feature", 
           c("name", "type"), 
           sep="."),
         "pfam" = paste(
           "polymer_entities.polymer_entity_instances.polymer_entity.rcsb_polymer_entity_annotation", 
           c("annotation_id", "name", "type"), 
           sep="."),
         "ligandChainId" = "nonpolymer_entities.nonpolymer_entity_instances.rcsb_nonpolymer_entity_instance_container_identifiers.auth_asym_id",
         "ligandId" = "nonpolymer_entities.nonpolymer_entity_instances.nonpolymer_entity.nonpolymer_comp.chem_comp.id",
         "ligandName" = "nonpolymer_entities.nonpolymer_entity_instances.nonpolymer_entity.nonpolymer_comp.chem_comp.name",
         "source" = "polymer_entities.polymer_entity_instances.polymer_entity.rcsb_entity_source_organism.ncbi_scientific_name", 
         "structureTitle" = "struct.title",
         "citation" = paste("rcsb_primary_citation", 
                            c( "rcsb_authors",
                               "rcsb_journal_abbrev",
                               "year"), 
                            sep="."),
         "rObserved" = "refine.ls_R_factor_obs", 
         "rFree" = "refine.ls_R_factor_R_free", 
         "rWork" = "refine.ls_R_factor_R_work",
         "spaceGroup" = "symmetry.space_group_name_H_M"
  )
}

## Return a formatted table/data.frame
.format_tbl <- function(x, query.ids, anno.terms, unique=FALSE) {
  if(!"data" %in% names(x)) {
    stop("No data retrieved")
  }
  x <- x$data$entries
  pdb.ids <- sapply(x, function(x) x$entry$id)
  chain.ids <- lapply(x, function(x) {
    lapply(x$polymer_entities, function(y) {
      sapply(y$polymer_entity_instances, function(z) {
        id <- z$rcsb_polymer_entity_instance_container_identifiers$auth_asym_id
        if(is.null(id)) {
          id<- as.character(NA)
        }
        id
      })
    })
  })
  nchains <- sapply(chain.ids, function(x) length(unlist(x)))
  nchains.entity <- sapply(chain.ids, sapply, length)
#  ids <- paste(rep(pdb.ids, nchains), unlist(chain.ids), sep="_")
  
  ids <- rep(pdb.ids, nchains)
  chainId <- unlist(chain.ids)
  out <- data.frame(structureId=ids, chainId=chainId, stringsAsFactors=FALSE)


  if("chainLength" %in% anno.terms) {
    cl <- lapply(x, function(x) {
      lapply(x$polymer_entities, function(y) {
        sapply(y$polymer_entity_instances, function(z) {
           cl <- z$polymer_entity$entity_poly$rcsb_sample_sequence_length
           if(is.null(cl)) {
              cl <- as.integer(NA)
           }
           cl
        })
      })
    })

    #    cl <- rep(unlist(cl), unlist(nchains.entity))
    cl <- unlist(cl)
    out$chainLength <- cl
  }
  
  if("experimentalTechnique" %in% anno.terms) {
    em <- sapply(x, function(x) {
      em <- x$rcsb_entry_info$experimental_method
      if(is.null(em)) {
        em <- as.character(NA)
      }
      em
    })
    em <- rep(em, nchains)
    out$experimentalTechnique <- em
  }
  
  if("resolution" %in% anno.terms) {
    reso <- sapply(x, function(x) {
      reso <- x$rcsb_entry_info$resolution_combined[[1]]
      if(is.null(reso)) {
        reso <- as.numeric(NA)
      }
      reso
    })
    reso <- rep(reso, nchains)
    out$resolution <- reso
  }
  
  if("macromoleculeType" %in% anno.terms) {
    moltype <- lapply(x, function(x) {
      lapply(x$polymer_entities, function(y) {
        sapply(y$polymer_entity_instances, function(z) {
           typ <- z$polymer_entity$entity_poly$rcsb_entity_polymer_type
           if(is.null(typ)) {
              typ <- as.character(NA)
           }
           typ
        })
      })
    })
    #    moltype <- rep(unlist(moltype), unlist(nchains.entity))
    moltype <- unlist(moltype)
    out$macromoleculeType <- moltype
  }
  
  if("scopDomain" %in% anno.terms) {
    scop <- lapply(x, function(x) {
      lapply(x$polymer_entities, function(y) {
        sapply(y$polymer_entity_instances, function(z) {
          types <- sapply(z$rcsb_polymer_instance_feature, "[[", "type")
          s <- z$rcsb_polymer_instance_feature[[which(types=="SCOP")[1]]]$name
          if(is.null(s)) {
            s <- as.character(NA)
          }
          s
        })
      })
    })
    scop <- unlist(scop)
    out$scopDomain <- scop
  }
  
  if("pfam" %in% anno.terms) {
    pfam <- lapply(x, function(x) {
      lapply(x$polymer_entities, function(y) {
        sapply(y$polymer_entity_instances, function(z) {
          types <- sapply(z$polymer_entity$rcsb_polymer_entity_annotation, "[[", "type")
          p <- z$polymer_entity$rcsb_polymer_entity_annotation[[which(types=="Pfam")[1]]]$name
          if(is.null(p)) {
            p <- as.character(NA)
          }
          p
        })
      })
    })
    pfam <- unlist(pfam)
    out$pfam <- pfam
  }
  
  if("ligandChainId" %in% anno.terms) {
    lch <- lapply(x, function(x) {
      unlist( lapply(x$nonpolymer_entities, function(y) {
        sapply(y$nonpolymer_entity_instances, function(z) {
          id <- z$rcsb_nonpolymer_entity_instance_container_identifiers$auth_asym_id
          if(is.null(id)) {
             id <- as.character(NA)
          }
          id
        })
      }) )
    })
  }
  
  if("ligandId" %in% anno.terms) {
    lid <- lapply(1:length(x), function(i) {
      x <- x[[i]]
      lch <- lch[[i]]
      chain.ids <- unlist(chain.ids[[i]])
      if(!is.null(x$nonpolymer_entities)) {
        id <- unlist( lapply(x$nonpolymer_entities, function(y) {
          sapply(y$nonpolymer_entity_instances, function(z) {
            id <- z$nonpolymer_entity$nonpolymer_comp$chem_comp$id
            if(is.null(id)) {
              id <- as.character(NA)
            }
            id
          })
        }) )
        id <- tapply(id, lch, function(x) {
          count <- table(x)
          x <- unique(x)
          count <- count[x]
          x[count>1] <- paste(x[count>1], " (", count[count>1], ")", sep="")
          paste(x, collapse=",")
        })
        id <- id[chain.ids]
      }
      else {
        id <- rep(as.character(NA), length(chain.ids))
      }
    })
    lid <- unlist(lid)
    out$ligandId <- lid
  }
  
  if("ligandName" %in% anno.terms) {
    lname <- lapply(1:length(x), function(i) {
      x <- x[[i]]
      lch <- lch[[i]]
      chain.ids <- unlist(chain.ids[[i]])
      if(!is.null(x$nonpolymer_entities)) {
        nam <- unlist( lapply(x$nonpolymer_entities, function(y) {
          sapply(y$nonpolymer_entity_instances, function(z) {
            nam <- z$nonpolymer_entity$nonpolymer_comp$chem_comp$name
            if(is.null(nam)) {
              nam <- as.character(NA)
            }
            nam
          })
        }) )
        nam <- tapply(nam, lch, function(x) {
          count <- table(x)
          x <- unique(x)
          count <- count[x]
          x[count>1] <- paste(x[count>1], " (", count[count>1], ")", sep="")
          paste(x, collapse=",")
        })
        nam <- nam[chain.ids]
      }
      else {
        nam <- rep(as.character(NA), length(chain.ids))
      }
    })
    lname <- unlist(lname)
    out$ligandName <- lname  
  }
  
  if("source" %in% anno.terms) {
    src <- lapply(x, function(x) {
      lapply(x$polymer_entities, function(y) {
        sapply(y$polymer_entity_instances, function(z) {
           src <- sapply(z$polymer_entity$rcsb_entity_source_organism, function(z2) {
              src <- z2$ncbi_scientific_name
              if(is.null(src)) {
                 src <- as.character(NA)
              }
              src
           })
           if(length(src)>1) {
             paste(src, collapse="/")
           }
           else {
             src
           }
        })
      })
    })
    #    src <- rep(unlist(src), unlist(nchains.entity))
    src <- unlist(src)
    out$source <- src
  } 
  
  if("structureTitle" %in% anno.terms) {
    title <- sapply(x, function(x) {
      title <- x$struct$title
      if(is.null(title)) {
        title <- as.character(NA)
      }
      title
    })
    title <- rep(title, nchains)
    out$structureTitle <- title
  }
  
  if("citation" %in% anno.terms) {
    citation <- sapply(x, function(x) {
      aut <- x$rcsb_primary_citation$rcsb_authors
      if(!is.null(aut) && length(aut)>1) {
        aut <- paste(aut[[1]], ", et al.", sep="")
      }
      jrnl <- x$rcsb_primary_citation$rcsb_journal_abbrev
      year <- x$rcsb_primary_citation$year
      if(!is.null(year)) {
        year <- paste("(", year, ")", sep="")
      }
      citation <- paste(aut, jrnl, year)
      if(length(citation)==0) {
        citation <- as.character(NA)
      }
      citation
    })
    citation <- rep(citation, nchains)
    out$citation <- citation
  }
  
  if("rObserved" %in% anno.terms) {
    robs <- sapply(x, function(x) {
      robs <- x$refine[[1]]$ls_R_factor_obs
      if(is.null(robs)) {
        robs <- as.numeric(NA)
      }
      robs
    })
    robs <- rep(robs, nchains)
    out$rObserved <- robs    
  }
  
  if("rFree" %in% anno.terms) {
    rfree <- sapply(x, function(x) {
      rfree <- x$refine[[1]]$ls_R_factor_R_free
      if(is.null(rfree)) {
        rfree <- as.numeric(NA)
      }
      rfree
    })
    rfree <- rep(rfree, nchains)
    out$rFree <- rfree       
  }
  
  if("rWork" %in% anno.terms) {
    rwork <- sapply(x, function(x) {
      rwork <- x$refine[[1]]$ls_R_factor_R_work
      if(is.null(rwork)) {
        rwork <- as.numeric(NA)
      }
      rwork
    })
    rwork <- rep(rwork, nchains)
    out$rWork <- rwork       
  }
  
  if("spaceGroup" %in% anno.terms) {
    sg <- sapply(x, function(x) {
      sg <- x$symmetry$space_group_name_H_M
      if(is.null(sg)) {
        sg <- as.character(NA)
      }
      sg
    })
    sg <- rep(sg, nchains)
    out$spaceGroup <- sg     
  }
  
  # Filter the table
  query.ids <- sub("\\.", "_", query.ids) ## allow PDBId.chainId
  inds <- lapply(query.ids, function(x) {
    if(nchar(x)==4) {
      which(out$structureId %in% x) 
    }
    else {
      which(paste(out$structureId, out$chainId, sep="_") %in% x)
    }
  })
  chk <- sapply(inds, length)
  if(any(chk==0)) {
    warning(paste("Annotation data could not be found for PDB ids:\n  ",
                  paste(unique(query.ids[chk==0]), collapse=", ")))
  }
  query.ids <- query.ids[chk>0]
  out <- out[unique(unlist(inds)), , drop=FALSE]
  
  if(unique) {
    # Fold the table based on PDB IDs
    out <- tapply(1:nrow(out), factor(out$structureId, levels=unique(out$structureId)),
     function(i){
       out <- out[i, , drop=FALSE]
       labs <- colnames(out); names(labs) <- labs
       cols <- lapply(labs, function(j) {
          if(j %in% c("chainId", "macromoleculeType", "chainLength", 
                      "scopDomain", "pfam", "ligandId", "ligandName", "source")) {
            paste(out[, j], collapse=";")
          }
          else {
            out[1, j]
          }
       })
       as.data.frame(cols, stringsAsFactors=FALSE)
    }, simplify=FALSE)
    out <- do.call(rbind, out)
    rownames(out) <- NULL
  }
  
  req.terms <- c("structureId", "chainId")
  out <- out[, c(req.terms, setdiff(anno.terms, c(req.terms, "ligandChainId"))), 
             drop=FALSE]
  out
}
