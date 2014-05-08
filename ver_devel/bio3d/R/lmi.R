lmi <- function (trj, grpby = NULL, ncore=1) {
    ncore <- setup.ncore(ncore)

# rm:r-value matrix
    cm <- var(trj)
    l <- dim(cm)[1]/3
    rm <- matrix(nrow=l, ncol=l)
    d <- 3
    ij <- pairwise(l)

# mclapply or lapply
    if (ncore > 1) {
        lmiapply = mclapply
    } else {
        lmiapply = lapply
    }  

# list1: marginal-covariance 
    list1 <- lmiapply(1:l, function(i) det(cm[(3*i-2):(3*i), (3*i-2):(3*i)]) )
    dm <- unlist(list1)

# list2: pair-covariance
    list2 <- lmiapply(1:nrow(ij), function(i) {
        x <- det(cm[c((3*ij[i,1]-2):(3*ij[i,1]),(3*ij[i,2]-2):(3*ij[i,2])), c((3*ij[i,1]-2):(3*ij[i,1]),(3*ij[i,2]-2):(3*ij[i,2]))])
        y <- 1/2 * (log(dm[ij[i,1]]) + log(dm[ij[i,2]]) - log(x))
        (1 - exp(-2 * y / d))^(1/2)
        }
    )
    list2 <- unlist(list2)

    for (k in 1:nrow(ij)) {
        rm[ij[k, 1], ij[k, 2]] <- list2[k]
    }
    rm[lower.tri(rm)] = t(rm)[lower.tri(rm)]
    diag(rm) <- 1

# group by or not
    if (!is.null(grpby)) {
        if (ncol(trj) != (length(grpby) * 3)) 
            stop("dimension miss-match in 'trj' and 'grpby', check lengths")
        inds <- bounds(grpby, dup.inds = TRUE)
        l <- dim(inds)[1]
        m <- matrix(, ncol = l, nrow = l)
        ij <- pairwise(l)
# list3: lmi 
        list3 <- lmiapply(1:nrow(ij), function(k) max(rm[(inds[ij[k, 1], "start"]:inds[ij[k, 1], "end"]), (inds[ij[k, 2], "start"]:inds[ij[k, 2], "end"])], na.rm = TRUE))
        list3 <- unlist(list3)

        for (k in 1:nrow(ij)) {
            m[ij[k, 1], ij[k, 2]] <- list3[k]
        }
        m[lower.tri(m)] = t(m)[lower.tri(m)]
        diag(m) <- 1; rm=m
    }
    class(rm) = c("dccm", "matrix")
    return(rm)
}


