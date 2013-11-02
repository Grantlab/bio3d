\name{pdb}
\alias{pdb}
\docType{data}
\title{
Pdb Example File
}
\description{
This is a pdb file used in the examples
}
\usage{data(pdb)}
\format{
  The format is:
List of 8
 $ atom      : chr [1:136, 1:13] "1" "2" "3" "4" ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : NULL
  .. ..$ : Named chr [1:13] "eleno" "elety" "alt" "resid" ...
  .. .. ..- attr(*, "names")= chr [1:13] "2" "4" "5" "6" ...
 $ het       : chr[0 , 1:13] 
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : NULL
  .. ..$ : Named chr [1:13] "eleno" "elety" "alt" "resid" ...
  .. .. ..- attr(*, "names")= chr [1:13] "2" "4" "5" "6" ...
 $ helix     :List of 4
  ..$ start: num(0) 
  ..$ end  : num(0) 
  ..$ chain: chr(0) 
  ..$ type : chr(0) 
 $ sheet     :List of 4
  ..$ start: num(0) 
  ..$ end  : num(0) 
  ..$ chain: chr(0) 
  ..$ sense: chr(0) 
 $ seqres    : NULL
 $ xyz       : num [1:408] 51.9 52.8 42.8 51.3 51.5 ...
 $ xyz.models: NULL
 $ calpha    : logi [1:136] FALSE TRUE FALSE FALSE FALSE TRUE ...
 - attr(*, "class")= chr "pdb"
}
\keyword{datasets}
