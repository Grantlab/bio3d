
####################################
####     webGL functions        ####
####################################

output$nmaWebGL  <- renderWebGL({
  pdb <- get_pdb()
  modes <- calcModes()

  trj <- mktrj(modes, mode=as.numeric(input$mode_choice), rock=FALSE)
  n <- nrow(trj)
  
  amalcol <- function(x) {
    col <- rep("grey50", length(x))
    col[1] <- "blue"
    col[length(col)] <- "red"
    return(col)
  }
  
  magcol <- function() {
    rf <- rmsf(trj)
    return(t(replicate(n, vec2color(rf, c('blue', 'red')),simplify=TRUE)))
  }

  class(trj)  <- 'xyz'
  col <- switch(input$viewColor2,
                'mag' = magcol(), # vec2coxlor(rmsf(m)), #!! col=col, type=2
                'amalgam' = amalcol(1:n),
                'default' = colorRampPalette(c('blue', 'gray', 'red'))(n)
                )
  
  typ <- switch(input$viewColor2,
                'mag' = 2,
                'amalgam' = 1,
                'default' = 1
                )
  
    view.xyz(trj, bg.col=input$viewBGcolor2, col=col, add=TRUE, type=typ)
})

