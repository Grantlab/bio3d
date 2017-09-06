uniprot <- function(accid) {
  oops <- requireNamespace("XML", quietly = TRUE)
  if(!oops)
    stop("Please install the XML package from CRAN")
  
  url <- paste('http://www.uniprot.org/uniprot/', accid, '.xml', sep="")

  tmpfile <- tempfile()
  download.file(url, tmpfile) ##, method="wget")
  xml <- XML::xmlRoot(XML::xmlParse(tmpfile))
  
  node.names <- XML::xmlSApply(xml[[1]], XML::xmlName)
  
  ## acession
  inds <- which(node.names=="accession")
  accession <- NULL
  for(i in 1:length(inds))
    accession <- c(accession, XML::xmlValue(xml[[1]][[inds[i]]]))
    
  ## and name
  inds <- which(node.names=="name")
  name <- NULL
  for(i in 1:length(inds))
    name <- c(name, XML::xmlValue(xml[[1]][[inds[i]]]))
  
  ## sequence
  inds <- which(node.names=="sequence")
  sequence <-  gsub("\n", "", XML::xmlValue(xml[[1]][[inds]]))
  
  ## organism
  inds <- which(node.names=="organism")
  node <- xml[[1]][[inds]]
  organism <- NULL
  tmpl <- unlist(XML::xmlApply(node, XML::xmlAttrs))
  if("scientific" %in% tmpl)
      organism <- XML::xmlValue(node[[which(tmpl %in% "scientific")]])
  if("common" %in% tmpl)
      organism <- c(organism, XML::xmlValue(node[[which(tmpl %in% "common")]]))
    
  ## taxon
  inds <- which(node.names=="organism")
  node <- xml[[1]][[inds]]
  taxon <- NULL
  for ( i in 1:XML::xmlSize(node[['lineage']]) ) {
    taxon <- c(taxon, XML::xmlValue(node[['lineage']][[i]]))
  }
  
  ## protein
  node <- xml[[1]][['protein']]
  fullName <- XML::xmlValue(node[['recommendedName']][['fullName']])
  shortName <- XML::xmlValue(node[['recommendedName']][['shortName']])
  
  ## gene
  node <- xml[[1]][['gene']]
  gene <- XML::xmlValue(node[[1]])

  ## dbReference
  inds <- which(node.names=="dbReference")
  dbref <- list()
  for(i in 1:length(inds)) {
      node <- xml[[1]][[inds[i]]]
      dbref[[i]] <- XML::xmlAttrs(node)
  }
  
  dbref <- unlist(dbref)
  type.inds <- names((dbref)) == "type"
  id.inds <- names((dbref)) == "id"
  dbref <- data.frame(type=dbref[type.inds], id=dbref[id.inds],
                      stringsAsFactors=FALSE)
  
  out <- list(accession = accession, name = name,
              fullName = fullName, shortName = shortName,
              sequence = sequence, gene = gene,
              organism = organism, taxon = taxon,
              dbref=dbref)
  
  return(out)
}
