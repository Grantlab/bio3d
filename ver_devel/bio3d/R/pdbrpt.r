# function: pdbrpt (pdb reports)
# date: 08/06/2013
# library(bio3d)
# library(XML)

## pdbrpt can return pdb reports in csv form.

pdbrpt <- function(ids, file="pdbrpt.csv", terms) 
{
    if (missing(ids)) 
    {
        stop("please specify PDB ids for reporting")
    }
    if (any(nchar(ids) != 4)) 
    {
        warning("ids should be standard 4 character PDB formart: trying first 4 char...")
        ids <- substr(basename(ids), 1, 4)
    }
    
    ids <- unique(ids)
    ids <- paste(ids, collapse=",")
    c = c("structureId", "experimentalTechnique", "resolution", "chainId", "ligandId", "ligandName", "source", "scopDomain", "classification", "compound", "title", "citationAuthor", "journalName", "publicationYear")
    
    c1 <- paste(c, collapse=",")
    cmd <- paste("http://www.rcsb.org/pdb/rest/customReport?pdbids=", ids, "&customReportColumns=", c1, "&ssa=n&primaryOnly=1", sep = "")
    tmpfile <- tempfile()
    test <- download.file(cmd, tmpfile)
    if(test!=0) stop("Files could not be downloaded from PDB, check returned value")
    doc = xmlRoot(xmlTreeParse(tmpfile))
    tmp1 = xmlApply(doc, function(x) xmlApply(x, xmlValue))
    tmp2 <- tmp1[[1]]
    for(i in 2:length(tmp1))
    {
        tmp2 <- mapply(c, tmp2, tmp1[[i]], SIMPLIFY=FALSE)
    }

# list tmp2 to matrix tmp3
    tmp3 <- sapply(tmp2, unlist)

# information of unique ids
    pdb.ids <- as.vector(tmp2[[1]])
    dup.len <- rle(pdb.ids)
    unq.ids <- dup.len$values[ (dup.len$lengths==1) ]
    new.tbl <- as.matrix(tmp3[pdb.ids %in% unq.ids,])

# if there was only one unique id, we need the transformation of the matrix new.tbl
    if(length(unq.ids)==1)
    {
        new.tbl <- t(new.tbl)
    }
    dup.rows <- bounds(which(duplicated(pdb.ids)))

# information of duplicated ids
    for( i in 1:nrow(dup.rows))
    {
        l <- tmp3[ (dup.rows[i,"start"]-1):dup.rows[i,"end"], ]
        row.store <- NULL
        for(j in 1:ncol(l))
        {
         row.store <- c(row.store, paste(unique(as.vector(l[,j])), collapse=", "))
        }
    new.tbl <- rbind(new.tbl, row.store)
    }

# format the colnames
    colnames(new.tbl) <- colnames(tmp3)

# change colnames (dimEntity.structureId -> structureId)
    for( i in 1:ncol(l))
    {
    a <- unlist(strsplit(colnames(new.tbl)[i], "\\."))
    colnames(new.tbl)[i] <- a[2]
    }

    citation <- NULL
    lig.auth <- new.tbl[,"citationAuthor"]
    lig.year <- new.tbl[,"publicationYear"]
    lig.jnal <- new.tbl[,"journalName"]

    for(i in 1:length(lig.auth))
    {
        citation <- c(citation, paste( unlist(strsplit(lig.auth[i], ","))[1], " et al. ", lig.jnal[i], " (", lig.year[i],")",sep=""))
    }
    new.tbl <- cbind(new.tbl, citation)
    rownames(new.tbl) <- new.tbl[, "structureId"]
    allterms <- c("structureId", "experimentalTechnique", "resolution", "chainId", "ligandId", "ligandName", "source", "scopDomain", "classification", "compound", "title", "citation", "citationAuthor", "journalName", "publicationYear")
    terms <- match.arg(terms, allterms, several.ok=TRUE)
    out.tbl <- as.matrix(new.tbl[, terms])
    colnames(out.tbl) <- terms
# output csv file
    write.csv(out.tbl, row.names=TRUE, file=file)
}


