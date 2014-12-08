`seqaln.pair` <-
function(aln, extra.args = "", ...) {
  cl <- match.call()
  l <- seqaln(aln,
              extra.args= paste("-matrix",
                system.file("matrices/custom.mat", package="bio3d"),
                "-gapopen -3.0 ",
                "-gapextend -0.5",
                "-center 0.0", extra.args), ... )

  if(!all((seqidentity(l))==1)) {
    warning("Sequences are not identical, use seqaln()")
  }
  l$call <- cl
  return(l)
}

