
############################
## Geostas domain analysis
###########################
geostas2 <- reactive({
  pdb <- final_pdb()
  modes <- calcModes()
  
  progress <- shiny::Progress$new(session, min=1, max=5)
  on.exit(progress$close())
  
  progress$set(message = 'Finding domains',
               detail = 'Please wait')
  
  progress$set(value = 2)
  gs <- geostas(modes, k=input$ndomains, m.inds=seq(7, input$nmodes+6), ncore=1)
  print(gs)
  return(gs)
})

gstrj <- reactive({
  path <- data_path()
  modes <- calcModes()
  pdb <- trim(final_pdb(), "calpha")
  gs <- geostas2()
  
  i <- 7
  file <- paste0(path, "/gs-mode_", i, ".pdb")
  
  mktrj(modes, mode=i,
        chain=gs$grps,
        resno=pdb$atom$resno,
        resid=pdb$atom$resid,
        rock=TRUE,
        file=file)
  return(file)
})

output$geostas2zip = downloadHandler(
  filename = "geostas.zip",
  content = function(file) {
    zip(file, files=gstrj(), flags = "-9Xj")
  })


output$geostasWebGL  <- renderWebGL({
  pdb <- final_pdb()
  modes <- calcModes()
  gs <- geostas2()
  col <- gs$grps
  view.xyz(modes$xyz, col=col, add=TRUE, type=2)
})
