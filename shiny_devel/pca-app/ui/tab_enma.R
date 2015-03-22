
tabPanel("4. NMA",
         fluidRow(
           plotOutput("nma_plot")
           ),

         fluidRow(
           column(3,
                  checkboxInput('show_options', 'More options', value=FALSE)
                  ),
           column(3),
           column(3),
           column(3,
                  downloadButton('nmaplot2pdf', "Download Plot PDF")
                  )
           ),

         fluidRow(
           plotOutput("rmsip_plot")
           ),
         
         fluidRow(
           column(3,
                  checkboxInput('show_options', 'More options', value=FALSE)
                  ),
           column(3),
           column(3),
           column(3,
                  downloadButton('rmsipplot2pdf', "Download Plot PDF")
                  )
           )
           
         )
