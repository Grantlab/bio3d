
############################
## Geostas domain analysis
###########################
geostas2 <- reactive({
  pdb <- get_pdb()
  modes <- calcModes()
  print(modes)
  
  progress <- shiny::Progress$new(session, min=1, max=5)
  on.exit(progress$close())
  
  progress$set(message = 'Calculation in progress',
               detail = 'This may take a moment...')
  
  progress$set(value = 2)
  gs <- geostas(modes, k=input$ndomains, m.inds=seq(7, input$nmodes+6), ncore=1)
  print(gs)
  return(gs)
  
  ##file <- "gstrj7.pdb"
  #  mktrj(modes, mode=7,
  #        chain=gs$grps,
  #        resno=pdb.ca$atom$resno,
  #        resid=pdb.ca$atom$resid,
  #        rock=FALSE,
  #        file=file)
  
  #zip("tmp_gs.zip", file, flags="-FS")
  #return("tmp_gs.zip")
})

  output$geostas2zip = downloadHandler(
    filename = "domains.zip",
    content = function(file) {
      file.copy(geostas2(), file)
    })
  

output$geostasWebGL  <- renderWebGL({
  pdb <- get_pdb()
  modes <- calcModes()
  gs <- geostas2()
  col <- gs$grps
  view.xyz(modes$xyz, col=col, add=TRUE, type=2)
})
