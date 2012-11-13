"wiki.tbl" <-
function(mat) {
  cat(paste("| **",
            paste(colnames(mat), collapse = "** | **"),
            "** |",sep=""), sep="\n")
  for(i in 1:nrow(mat)) {
    cat(paste("|",paste(mat[i,], collapse = " | "),"|"),sep="\n")
  }
}

