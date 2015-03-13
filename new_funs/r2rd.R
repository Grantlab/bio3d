r2rd <- function(x, test=FALSE) {
   devtools::create("tmp.package")
   system(paste("cp", x, "tmp.package/R"))
   roxygen2::roxygenize("tmp.package")
   system(paste("cp tmp.package/man/", sub(".R$",".Rd",x), " .", sep="")) 
   if(test) {
      devtools::load_all("tmp.package")
      help(sub(".R$","",x), package="tmp.package")
   }
   unlink("tmp.package", recursive=TRUE)
}

#r2rd(x="dist.dccm.R", test=TRUE)
