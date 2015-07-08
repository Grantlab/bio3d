
####################################
####     webGL functions        ####
####################################

output$nmaWebGL  <- renderWebGL({
  path <- data_path()
  pdb <- rv$final_pdb
  pdb <- trim(pdb, "calpha")
  modes <- rv$modes
  mag <- as.numeric(input$mag)
  step <- mag/8
  
  i <- as.numeric(input$mode_choice)
  fname  <- paste0(path, '/', 'mode', input$mode_choice, '.pdb')
  trj <- mktrj(modes, mode=i, file=fname,
             mag=mag, step=step,
             b=modes$fluctuations,
             resno=pdb$atom$resno,
             resid=pdb$atom$resid,
             chain=pdb$atom$chain,
             rock=FALSE)
  
  
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
  
    view.xyz(trj, bg.col=input$viewBGcolor2, col=col, d.cut=8)
})

