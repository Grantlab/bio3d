check.utility <- function(x = c("muscle", "clustalo", "dssp", "stride", "mustang", "makeup"), 
   quiet = TRUE) {
  
  utilities <- match.arg(x, several.ok = TRUE)

  ##- Check on missing utility programs
  if('dssp' %in% utilities) {
    utilities <- c(utilities, 'mkdssp')
    missing.util <- nchar(Sys.which(utilities)) == 0
    missing.util['dssp'] <- all(missing.util[c('dssp', 'mkdssp')])
    missing.util <- missing.util[-length(missing.util)]
  } else {
    missing.util <- nchar(Sys.which(utilities)) == 0
  }
  if( any(missing.util) ) {
    if(!quiet) {
       warning(paste0("  Checking for external utility programs failed\n",
         "    Please make sure '", paste(names(missing.util[missing.util]), collapse="', '"),
         "' is in your search path, see:\n",
         "    http://thegrantlab.org/bio3d/tutorials/installing-bio3d#utilities"))
    }
    pass = FALSE
  } else {
    if(!quiet) cat("External utility programs found\n")
    pass = TRUE
  }
  invisible(pass)
}  
