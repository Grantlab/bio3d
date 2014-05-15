combine.sel <- 
 function(sel1=NULL, sel2=NULL, op="AND", verbose=TRUE) {
   if(is.null(c(sel1, sel2))) {
      stop("Invalid selection")
   }
   if(is.null(sel1))
      return(sel2)
   if(is.null(sel2))
      return(sel1)
 
   op.tbl <- c(rep("AND",3), rep("OR",4), rep("NOT",4))
   op <- op.tbl[match(op, c("AND","and","&","OR","or","|","+","NOT","not","!","-"))]
   sel <- switch(op, 
      "AND"= {
         if(verbose) { cat(" Intersection of sel1 and sel2\n") }
         list(atom = intersect(sel1$atom, sel2$atom), 
              xyz = intersect(sel1$xyz, sel2$xyz))
      },
      "OR" = {
        if(verbose) { cat("Union of sel1 and sel2\n") }
         list(atom = sort(union(sel1$atom, sel2$atom)), 
              xyz = sort(union(sel1$xyz, sel2$xyz)))
      },
      "NOT" = {
        if(verbose) { cat("Sel2 is subtracted from sel1\n") }
         list(atom = setdiff(sel1$atom, sel2$atom), 
               xyz = setdiff(sel1$xyz, sel2$xyz))
      },
      stop("Unknown operation") )
   if(verbose) { cat(paste(" *  Selected a total of:", length(sel$atom), "atoms  *\n")) }
   class(sel)="select"
   return(sel)
}
