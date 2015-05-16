tabPanel("5. eNMA", icon=icon("arrow-right"),
         tags$style(type="text/css", "body {padding-top: 80px;}"),
         

         ### WebGL visualization
         fluidRow(
           column(4,
                  wellPanel(
                    h4('Normal Modes Visualization'),
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
                    ##actionButton('viewUpdate2', label='Refresh', icon=icon('undo')),
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
         
         ### Fluctuation plot
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
                    
                    checkboxInput('cluster', 'Color by clustering', value=TRUE),
                    radioButtons("group_by", "Cluster by",
                                 c("RMSD" = "rmsd",
                                   "RMSIP" = "rmsip",
                                   "PC distance" = "pc_space"
                                   ),
                                 ##"bhat" = "bhat"),
                                 inline=TRUE),
                    sliderInput("nclusts", "N clusters:",
                                min = 1, max = 10, value = 3),
                    checkboxInput('show_options1', 'More options', value=FALSE),
                    downloadButton('nmaplot2pdf', "Download Plot PDF")
                    )
                  ),
           column(8,
                  plotOutput("nma_fluctplot")
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
                    ),

             column(3,
                    wellPanel(
                      checkboxInput("toggle_all2", "Toggle all", TRUE),
                      uiOutput("checkboxgroup_label_ids2")
                      )
                    )
                 )
           ),
         
         
         ### Heatmaps
         fluidRow(
           column(6,
                  h4("RMSIP heatmap"),
                  plotOutput("rmsip_heatmap2")
                  ),
           column(6,
                  h4("RMSD heatmap"),
                  plotOutput("rmsd_heatmap2")
                  )
           ),
         
         
         fluidRow(
           column(3,
                  checkboxInput('show_options3', 'More options', value=FALSE)
                  ),
           column(3,
                  downloadButton('nma_rmsip_heatmap2pdf', "Download Plot PDF")
                  ),
           column(3),
           column(3,
                  downloadButton('nma_rmsd_heatmap2pdf', "Download Plot PDF")
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
                      )
                    )
             )
           ),
         
         
         ## Cluster dendrogram
         fluidRow(
           column(4,
                  wellPanel(
                    h4('Cluster dendrogram'),
                    checkboxInput('show_confplot2', 'Show PC conformer plot', value=FALSE),
                    checkboxInput('show_options2', 'More options', value=FALSE),

                    radioButtons("group_by", "Cluster by",
                                 c("RMSD" = "rmsd",
                                   "RMSIP" = "rmsip",
                                   "PC distance" = "pc_space"
                                   ),
                                 ##"bhat" = "bhat"),
                                 inline=TRUE),
                    
                    downloadButton('nmadendrogram2pdf', "Download Plot PDF")
                    )
                  ),
           column(8,
                  plotOutput("dendrogram2")
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
                      sliderInput("width-pcload", "Width",
                                  min = 4, max = 12, value = 7, step=0.5),
                      sliderInput("height-pcload", "Height",
                                  min = 4, max = 12, value = 7, step=0.5)
                      )
                    )
             )
           ),


         ### PCA vs NMA (RMSIP)
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
