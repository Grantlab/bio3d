formula2mass <- function(form, sum.mass=TRUE) {
  errmsg <- "Please specify a formula on the form: 'C3 H5 N O1'"
  
  if(class(form)!="character" || missing(form))
    stop(errmsg)
    
  eles <- unlist(strsplit(form, "\\ "))
  eles <- eles[sapply(eles, nzchar)]
  no.num <- !grepl("[0-9]",eles)
  eles[no.num] <- paste(eles[no.num], "1", sep="")

  elemass <- function(ele) {
    num <- gsub("[A-z]","",ele)
    cha <- gsub("[0-9]","",ele)
    M <- atom2mass(cha, na.rm = FALSE)*as.numeric(num)
    return(M)
  }
  M <- sapply(eles, elemass)
  if(any(is.na(M))) stop(errmsg)
  if(sum.mass) M <- sum(M)
  return(M)
}
