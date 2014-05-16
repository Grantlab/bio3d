are.symb <- function(x) {
  to.return <- (x %in% elements$symb)
  to.return[is.na(x)] <- NA
  return(to.return)
}

connectivity <- function(...)
  UseMethod("connectivity")

connectivity.default <- function(eleno.1, eleno.2, ...){
  if(missing(eleno.1) | missing(eleno.2))
    stop("Please specify 'eleno.1' and 'eleno.2'")
  if(length(eleno.1) != length(eleno.2))
    stop("'eleno.1' and 'eleno.2' must have the same length")
  con <- data.frame(eleno.1 = eleno.1, eleno.2 = eleno.2, ...)
  con <- con[order(eleno.1, eleno.2),]
  attr(con, "class") <- c("connectivity", "data.frame")
  return(con)
}

connectivity.connectivity <- function(x, ...)
  connectivity.default(x$eleno.1, x$eleno.2, ...)

connectivity.xyz <- function(x, elesy, safety = 1.2, by.block = FALSE, ...){
  if(missing(elesy))
    stop("Please specify 'elesy'")
  if(!inherits(x, "xyz"))
    stop("'x' must be an object of class 'xyz'")
  if(nrow(x) != 1)
    stop("'x' must be a single row 'xyz' matrix")
  if(length(elesy) != length(x)/3)
    stop("'x' and 'elesy' must have matching lengths")
  if(!all(are.symb(elesy) & !is.na(elesy)))
    elesy[!are.symb(elesy) & is.na(elesy)] <- "Xx"
  
  radii <- elements[match(elesy, elements$symb), "rcov"]*safety
  x <- as.data.frame(matrix(x, ncol=3, byrow=TRUE,
                            dimnames = list(NULL, c("x1","x2","x3"))))
  data <- cbind(x, radii)

  findCon <- function(data, order = TRUE) {
    nat <- nrow(data)
    if(nat==0) return(NULL)
    r <- sqrt(
        outer(data$x1, data$x1, "-")^2 +
        outer(data$x2, data$x2, "-")^2 +
        outer(data$x3, data$x3, "-")^2
    )
    bond.dist <- outer(data$radii, data$radii, "+")
    M <- lower.tri(r) & (r < bond.dist)
    if(all(!M)) return(NULL)
    eleno <- matrix(rownames(data), nrow = nat, ncol = nat)   
    eleno.1 <- as.integer(t(eleno)[M])
    eleno.2 <- as.integer(  eleno [M])
    con <- data.frame(eleno.1=eleno.1, eleno.2=eleno.2)
    if(order) con <- con[order(con$eleno.1, con$eleno.2),]
    return(con)
  }
  
  if(!by.block) {
    con <- findCon(data)
  }
  else {
    get.con <- function(shift = c(0,0,0), x, radii, step) {
      coords.range <- apply(x, 2, range)
      dimnames(coords.range) <- list(c("min", "max"), c("x", "y", "z"))
      coords.range <- t(t(coords.range) - shift)
      
      x.cuts <- seq(coords.range["min","x"], coords.range["max","x"] + step, step)
      y.cuts <- seq(coords.range["min","y"], coords.range["max","y"] + step, step)
      z.cuts <- seq(coords.range["min","z"], coords.range["max","z"] + step, step)
      
      x.index <- cut(x$x1, x.cuts, include.lowest = TRUE)
      y.index <- cut(x$x2, y.cuts, include.lowest = TRUE)
      z.index <- cut(x$x3, z.cuts, include.lowest = TRUE)
      
      f <- interaction(x.index, y.index, z.index)
      
      data <- cbind(x, radii)
      
      con <- split(data, f)
      con <- lapply(con, findCon, order = FALSE)
      con <- do.call(rbind, con)
      
      return(con)
    }
    
    width <-10
    shift <- expand.grid(0:1,0:1,0:1)*width/2
    
    con <- apply(shift, 1, get.con, x, radii, width)
    if(!is.null(con)){
      con <- unique(do.call(rbind, con))
      con <- con[order(con$eleno.1, con$eleno.2),]
    }
    
    rownames(con) <- NULL
  }
  return(con)
}

connectivity.pdb <- function(x, safety = 1.2, by.block = TRUE, ...) {
  if(!is.pdb(x))
    stop("'x' must be an object of class 'pdb'")   
  x <- connectivity.xyz(x = x$xyz, elesy = x$atom$elesy,
                        safety = safety, by.block = by.block)
  x$con$eleno.1 <- x$atom$eleno[x$con$eleno.1]
  x$con$eleno.2 <- x$atom$eleno[x$con$eleno.2]
  return(x)
}

## Setter
"connectivity<-" <- function(object, value)
  UseMethod("connectivity<-", object)

# Setter for object of class pdb
# Set the connectivity component of an object of class pdb
# Not really need as the "$" syntaxe can be used.
# However, the syntaxe "x$con <- value" will completly remove the cell component
# if value is NULL. This is not really a probleme as x$con will still return NULL 
# if this component is removed. To avoid this, we define a setter function
"connectivity<-.pdb" <- function(object, value, ...) {
  object["con"] <- list(value)
  return(object)
}  

