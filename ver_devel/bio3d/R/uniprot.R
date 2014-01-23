uniprot <- function(accid) {
  
  url <- paste('http://www.uniprot.org/uniprot/', accid, '.xml', sep="")

  tmpfile <- tempfile()
  download.file(url, tmpfile) ##, method="wget")
  xml <- xmlRoot(xmlParse(tmpfile))
  
  node.names <- xmlSApply(xml[[1]], xmlName)
  
  ## acession
  inds <- which(node.names=="accession")
  accession <- NULL
  for(i in 1:length(inds))
    accession <- c(accession, xmlValue(xml[[1]][[inds[i]]]))
    
  ## and name
  inds <- which(node.names=="name")
  name <- NULL
  for(i in 1:length(inds))
    name <- c(name, xmlValue(xml[[1]][[inds[i]]]))
  
  ## sequence
  inds <- which(node.names=="sequence")
  sequence <-  gsub("\n", "", xmlValue(xml[[1]][[inds]]))
  
  ## organism
  inds <- which(node.names=="organism")
  node <- xml[[1]][[inds]]
  organism <- xmlValue(node[[1]])
  
  ## taxon
  inds <- which(node.names=="organism")
  node <- xml[[1]][[inds]]
  taxon <- NULL
  for ( i in 1:xmlSize(node[['lineage']]) ) {
    taxon <- c(taxon, xmlValue(node[['lineage']][[i]]))
  }
  
  ## protein
  node <- xml[[1]][['protein']]
  fullName <- xmlValue(node[['recommendedName']][['fullName']])
  shortName <- xmlValue(node[['recommendedName']][['shortName']])
  
  ## gene
  node <- xml[[1]][['gene']]
  gene <- xmlValue(node[[1]])
  
  out <- list(accession = accession, name = name,
              fullName = fullName, shortName = shortName,
              sequence = sequence, gene = gene,
              organism = organism, taxon = taxon)
  
  return(out)
}
