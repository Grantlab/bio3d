tabPanel("HELP", icon=icon("wheelchair"),
         tags$style(type="text/css", "body {padding-top: 80px;}"),

         
         fluidRow(
           column(4,
                  wellPanel(
                    h3("Multiple structure analysis with Bio3D"),
                    hr(),
                    p("Bio3D@web provides a rapid and rigorous tool for comparative structure analysis of protein families."),
                    p("Start by entering a PDB code of interest to perform structure similarity search. Proceed to sequence/structre analysis and structure analysis by navigating in the above tabs."),
                    img(src="http://thegrantlab.org/bio3d/images/geostas_nma-v6-small.png",
                        width=250, style="display: block; margin-left: auto; margin-right: auto;")
                    )
                  ),
           
           column(4,
                  wellPanel(
                    h3("Contact... "),
                    hr()
                    
                    )
                  ),
           
           column(4,
                  wellPanel(
                    h3("to be filled in..."),
                    hr()
                    
                    )
                  )
           )


         
         
         )
        
