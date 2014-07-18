"read.crd" <- function(file, ...) {
  require(tools)
  ext <- file_ext(file)

  if(ext %in% c("crd")) {
    class(file)=c("character", "charmm")
  }

  if(ext %in% c("inpcrd", "rst")) {
    class(file)=c("character", "amber")
  }

  UseMethod("read.crd", file)
}

