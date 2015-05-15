tabPanel("5. eNMA", icon=icon("arrow-right"),
         tags$style(type="text/css", "body {padding-top: 80px;}"),
         

  fluidRow(
    column(4,
           wellPanel(
             h4('NM Trajectory Viewing Options'),
             checkboxInput('show_trj2', 'Show NM Trajectory', value=FALSE),
             
             selectInput('viewMode', 'Choose Mode:', choices=c(1:10)),
             uiOutput('struct_dropdown2'),
             
             radioButtons('viewColor2', label='Structure color',
                          choices=list(
                            'Amalgam' = 'amalgam',
                            'Magnitude'='mag',
                            'By Frame (blue->gray->red)'='default'
                            ),
                          selected='amalgam'),
             radioButtons('viewBGcolor2', label='Background color',
                          choices=list('Black'='black', 'White'='white'),
                          selected='white'),
             br(),
             #actionButton('viewUpdate2', label='Refresh', icon=icon('undo')),
             downloadButton('nmtraj', label='Download PDB Trajectory')
             )
           ),

    column(8,
           conditionalPanel(
             condition='input.show_trj2 == true',
             webGLOutput('nmaWebGL')
             )
           )
    ),

         
   fluidRow(
     column(4,
            wellPanel(
              h4('Residue fluctuations'),
              ##selectInput('', 'Choose Mode #:', choices=c("all", 1:10)),
              ##checkboxInput("toggle_rmsf2", "Show RMSF", FALSE)
              
              checkboxInput('spread', 'Spread lines', value=FALSE),
              checkboxInput('seqide', 'Sequence identity', value=FALSE),
              ##checkboxInput('signif', 'Show fluct signif', value=FALSE),
              checkboxInput('rm.gaps', 'Omit gaps', value=TRUE),
              
              checkboxInput('cluster', 'Color by clustering', value=FALSE),
              radioButtons("group_by", "group by",
                           c("RMSD" = "rmsd",
                             "RMSIP" = "rmsip"),
                           ##"bhat" = "bhat"),
                           inline=TRUE),
              sliderInput("nclusts", "N clusters:",
                          min = 1, max = 10, value = 3),
              checkboxInput('show_options1', 'More options', value=FALSE)
              )
            ),
     column(8,
            plotOutput("nma_fluctplot")
            )
     ),
         
     fluidRow(
       column(3),
       column(3),
       column(3),
       column(3,
              downloadButton('nmaplot2pdf', "Download Plot PDF")
              )
       ),
         
      conditionalPanel(
        condition = "input.show_options1 == true",
        fluidRow(
          column(3,
                 wellPanel(
                   sliderInput("width1", "width",
                               min = 4, max = 12, value = 7, step=0.5),
                   sliderInput("height1", "height",
                               min = 4, max = 12, value = 7, step=0.5)
                   )
                 )
          ) 
        ),

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
      ),
           
         
            
         br(),br(),
         h3("Cluster dendrogram"),
         hr(),
         
         fluidRow(
           column(12,
                  plotOutput("dendrogram2")
                  )
           ),

         fluidRow(
           column(3,
                  checkboxInput('show_options2', 'More options', value=FALSE)
                  ),
           column(3),
           column(3),
           column(3,
                  downloadButton('nmadendrogram2pdf', "Download Plot PDF")
                  )
           ),

         conditionalPanel(
           condition = "input.show_options2 == true",
           fluidRow(
             column(3,
                    wellPanel(
                      sliderInput("cex2", "Label size",
                                  min = 0.1, max = 3, value = 1, step=0.1),
                      sliderInput("margins2", "Margins",
                                  min = 3, max = 10, value = 5, step=1)
                      
                      )
                    ),
             
             column(3,
                    wellPanel(
                      sliderInput("width2", "Width",
                                  min = 4, max = 12, value = 7, step=0.5),
                      sliderInput("height2", "Height",
                                  min = 4, max = 12, value = 7, step=0.5)
                      )
                    )
             )
           ),

      
  

         br(),br(),
         h3("Heatmaps"),
         hr(),
         
         fluidRow(
           column(6,
                  checkboxInput('show_rmsd_heatmap', 'RMSD Heatmap', value=TRUE),
                  conditionalPanel(
                    condition = "input.show_rmsd_heatmap == true",
                    plotOutput("rmsd_heatmap2")
                    )
                  ),

           column(6,
                  checkboxInput('show_rmsip_heatmap', 'RMSIP Heatmap', value=TRUE),
                  conditionalPanel(
                    condition = "input.show_rmsip_heatmap == true",
                    plotOutput("rmsip_heatmap2")
                    )
                  )
           ),
         
         fluidRow(
           column(3,
                  checkboxInput('show_options3', 'More options', value=FALSE)
                  ),
           column(3,
                  downloadButton('nma_rmsd_heatmap2pdf', "Download Plot PDF")
                  ),
           column(3),
           column(3,
                  downloadButton('nma_rmsip_heatmap2pdf', "Download Plot PDF")
                  )
           ),

         
         conditionalPanel(
           condition = "input.show_options3 == true",
           fluidRow(
             column(3,
                    wellPanel(
                      sliderInput("cex3", "Label size",
                                  min = 0.1, max = 3, value = 1, step=0.1),
                      sliderInput("margins3", "Margins",
                                  min = 3, max = 10, value = 5, step=1)
                      
                      )
                    ),
             
             column(3,
                    wellPanel(
                      sliderInput("width3", "Width and height",
                                  min = 4, max = 12, value = 7, step=0.5)
                      #sliderInput("height3", "height",
                      #            min = 4, max = 12, value = 7, step=0.5)
                      )
                    )
             )
           )
         
         )
