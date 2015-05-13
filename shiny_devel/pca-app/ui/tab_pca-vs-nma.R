tabPanel("6. PCAvsNMA", icon=icon("arrow-right"),
         tags$style(type="text/css", "body {padding-top: 80px;}"),

         
    fluidRow(
      column(4,
             wellPanel(
               h4('PCA vs NMA comparison'),
               uiOutput('struct_dropdown3')
               
               )
             ),

      column(4,
             plotOutput("rmsip_plot3")
             ),

      column(4,
             verbatimTextOutput("rmsip_print3")
             )
      

      )
         
 )

