
## Wed Apr 11 13:41:53 EDT 2012
## Barry
##
## Utility functions for working with strings
## see also ls("package:base", pattern="str")
##     and packages "stringr", "gsubfn" and "brew"
##

##
##> strhead("egghead", 3)
##[1] "egg"
##> strhead("beagle", -1) # negative index
##[1] "beagl"
##> strtail(c("bowl", "snowboard"), 3) # vector-able in the first argument
##[1] "owl" "ard"
##

## "sprintf()", "format()" and "pretty()" are powerful functions for formatting strings.
## However, sometimes I miss the named template syntax as in Python or in Makefiles.
##
##> strsubst(
##+   "$(WHAT) is $(HEIGHT) meters high.", 
##+   list(
##+     WHAT="Berlin's teletower",
##+     HEIGHT=348
##+   )
##+ )
##[1] "Berlin's teletower is 348 meters high."
##> d <- strptime("2012-03-18", "%Y-%m-%d")
##> strsubst(c(
##+   "Be careful with dates.",
##+   "$(NO_CONV) shows a list.",
##+   "$(CONV) is more helpful."),
##+   list(
##+     NO_CONV=d,
##+     CONV= as.character(d)
##+   )
##+ )
##[1] "Be careful with dates."                                                                                        
##[2] "list(sec = 0, min = 0, hour = 0, mday = 18, mon = 2, year = 112, wday = 0, yday = 77, isdst = 0) shows a list."
##[3] "2012-03-18 is more helpful."                                                                               

## PARSE RAW TEXT
## > lines <- c(
##+     'VARIABLE LABELS weight "weight".',
##+     'VARIABLE LABELS altq "Year of birth".',
##+     'VARIABLE LABELS hhg "Household size".',
##+     'missing values all (-1).',
##+     'EXECUTE.'
##+ )
##> pat <- 'VARIABLE LABELS (?<name>[^\\s]+) \\"(?<lbl>.*)\\".$'
##> matches <- grepl(pat, lines, perl=TRUE)
##> strparse(pat, lines[matches])
##name     lbl             
##[1,] "weight" "weight"        
##[2,] "altq"   "Year of birth" 
##[3,] "hhg"    "Household size"





strhead <- function(s,n=1) {
  if(n<0) 
    substr(s,1,nchar(s)+n) 
  else 
    substr(s,1,n)
}

strtail <- function(s,n=1) {
  if(n<0) 
    substring(s,1-n) 
  else 
    substring(s,nchar(s)-n+1)
}

strsubst <- function(template, map, verbose=getOption("verbose")) {
  pat <- "\\$\\([^\\)]+\\)"
  res <- template
  map[["$"]] <- "$"
  m <- gregexpr(pat, template)
  idx <- which(sapply(m, function(x) x[[1]]!=-1)) # faster than 1:length(template)?
  for (i in idx) {
    line <- template[[i]]
    if(verbose) cat("input: |", template[[i]], "|\n")
    starts <- m[[i]]
    ml <- attr(m[[i]], "match.length")
    sym <- substring(line, starts+2, starts+ml-2)
    repl <- map[sym]
    idx1 <- is.null(repl)
    repl[idx1] <- sym[idx1]
    norepl <- substring(line, c(1, starts+ml), c(starts-1, nchar(line)))
    res[[i]] <- paste(norepl, c(repl, ""), sep="", collapse="") # more elegant?
    if (verbose) cat("output: |", res[[i]], "|\n")
  }
  return(res)
}


strparse <- function(pat, x) {
  parsed <- regexpr(pat, x, perl=TRUE)
  if (length(x)==1) {
    if(parsed[1]==-1) return(NULL)
    st <- attr(parsed, "capture.start")[1,]
    m <- substring(x, st, st + attr(parsed, "capture.length")[1,]-1)
    names(m) <- attr(parsed, "capture.names")
  } else {
    m <- do.call(rbind, lapply(seq_along(parsed), function(i) {
      if(parsed[i] == -1) return("")
      st <- attr(parsed, "capture.start")[i, ]
      substring(x[i], st, st + attr(parsed, "capture.length")[i, ] - 1) }))
    colnames(m) <- attr(parsed, "capture.names")
  }
  return(m)
}


strrecode <- function(pats, repls, x, ...) {
  res <- rep(NA, length(x))
  hits <- rep(FALSE, length(x))
  for (i in seq_along(pats)) {
    ##browser()
    new_hits <- grepl(pats[[i]],x[!hits],...)
    res[!hits][new_hits] <- repls[[i]]
    hits[!hits][new_hits] <- TRUE
    if(all(hits)) break
  }
  return(res)
}
