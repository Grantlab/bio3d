.arg.filter <- function(new.args, FUN=NULL, ...) {
  ##-- Simple list filtering for duplicate 
  ##    function input argument removal and validation.
  ##
  ## new.args = The new default args that can be overwritten 
  ##              by those in 'dots' (i.e. user supplied "...") 
  ##               E.G. "new.args=list(col=mydefualtcol, lwd=3)"
  ##
  ## FUN      = Function name from which allowed arguments are checked
  ##              and used to limit output of this function.
  ##              This is typically only required if there are multiple 
  ##              (sub)functions to be called each with other specific  
  ##              things in /dots.
  ##               E.G. allowed=names(formals( mysubfunction2call ))
  ##
  ## dots     = Full user supplied updated values typically 
  ##               this is the extra /dots values, i.e. list(...) 
  ##    
  ##   sse.default <- list(coil="gray", helix="purple", sheet="yellow")
  ##   sse.args <- arg.filter( sse.default, FUN=sse.vector )
  ##   col <- do.call('sse.vector', c(list(pdb=pdb), sse.args) )
  ##
  ## Returns entries of 'dots' updated with those in 'new.args'
  ##   that intersect with allowed FUN input args.
  
  dots <- list(...)
  ans <- c(dots, new.args[!names(new.args) %in% names(dots)])
  if(!is.null(FUN)) { ans <- ans[names(ans) %in% names(formals(FUN))] }
  return(ans)
}
