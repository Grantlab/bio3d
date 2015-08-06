tabPanel(
  "5. eNMA", icon=icon("arrow-right"),
  tags$style(type="text/css", "body {padding-top: 80px;}"),
  
  
  actionButton3("about_nmatab", "About this tab", icon=icon("comment"), cl="btn btn-warn btn-input action-button", style = "position: fixed; top: 14px; right: 16px; z-index: 2000;"),
  
  bsModal("modal_nma", "Ensemble Normal Mode Analysis", "about_nmatab", size = "large", 
          content=tags$div(
            
            p(HTML("In the NMA tab, normal mode analysis (NMA) for each individual structure in the ensemble is performed. The motions described by each normal mode can be visualized either in the browser, or by downloading a PDB trajectory file or a PyMOL session (.pse) containing a vector field representation. ")),

            p(HTML("The analyses of the resulting modes allows for direct comparison of predicted fluctuations,")),
            
            img(src="./images/enma-fluctuations.png", style="display: block; margin-left: auto; margin-right: auto;"),
            
            p(HTML("as well as clustering based on the similarity of the described motions.")),

          img(src="./images/enma-clustering.png", style="display: block; margin-left: auto; margin-right: auto;")
            )
          
          
          ),

  
  
### Filter structures
  fluidRow(
    column(4,
           wellPanel(
             h4("Filter structures"),
             helpText("Filter similar structures (by RMSD) to reduce the computational load for NMA. PDB IDs colored red in the dendrogram will be omitted from the calculation."), 
             
             ##numericInput("filter_rmsd","RMSD Cutoff", value = 1.0, step = 0.25),
             uiOutput("filter_rmsd"),

             verbatimTextOutput("filter_summary"),
             
             actionButton("run_nma", "Run Ensemble NMA", icon=icon("gears"),
                          style = "display: block; margin-left: auto; margin-right: auto;", 
                          class = "btn btn-success")

             )
           ),

    column(8,
           plotOutput("rmsd_dendrogram2")
           )
    ),
                     


    
### WebGL visualization
  hr(),
  fluidRow(
           column(4,
                  wellPanel(
                    bsPopover("popnma1",
                              "Normal mode visualization",
                              "Visualize the normal modes by toggeling the the <b>Show NM Trajectory</b> checkbox. <br><br>Download buttons enable visualization of the motions described by the principal component in external viewers such as PyMOL or VMD.",
                              placement = "right", trigger = "hover",
                              options = list(container = "body")),
                    
                    tags$div(id = "popnma1", icon("question-circle"),
                             style = "position: absolute; right: 25px; top: 5px;"
                             ),
                    
                    
                    h4('Normal Modes Visualization'),
                    checkboxInput('show_trj2', 'Show NM Trajectory', value=FALSE),
                    
                    selectInput('viewMode_nma', 'Choose Mode:', choices=c(1:10)),
                    uiOutput('struct_dropdown2'), ## viewStruct_nma

                    sliderInput("mag2", "Magnification factor:",
                                min = 1, max = 20, value = 5),

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
                    condition='input.show_trj2 == true && output.nmaIsDone',
                    webGLOutput('nmaWebGL')
                    )
                  )
    ),

         ### Fluctuation plot
         fluidRow(
           column(4,
                  wellPanel(
                    bsPopover("popnma2",
                              "Residue fluctuations",
                              "The fluctuations are calculated based on all 3N-6 normal modes (where N is the number of atoms). Magnitudes should be used by care as normal mode vectors are by definition without magnitude. <br><br>Use the option <b>Spread lines</b> provides fluctuation profiles not in top of eachother.",
                              placement = "right", trigger = "hover",
                              options = list(container = "body")),
                    
                    tags$div(id = "popnma2", icon("question-circle"),
                             style = "position: absolute; right: 25px; top: 5px;"
                             ),
                    
                    
                    
                    h4('Residue fluctuations'),
                    
                    checkboxInput('spread', 'Spread lines', value=FALSE),
                    checkboxInput('seqide', 'Sequence identity', value=FALSE),
                    ##checkboxInput('signif', 'Show fluct signif', value=FALSE),
                    checkboxInput('rm_gaps', 'Omit gap regions', value=TRUE),

                    checkboxInput('cluster', 'Color by clustering', value=TRUE),
                    radioButtons("cluster_by2", "Cluster by",
                                 c("RMSIP" = "rmsip",
                                   "RMSD" = "rmsd",
                                   "PC subspace" = "pc_space",
                                   "Sequence" = "sequence"
                                   ),
                                 ##"bhat" = "bhat"),
                                 inline=TRUE),

                    conditionalPanel(
                      condition = "input.cluster_by2 == 'rmsip'",
                      
                      ## K-selecter
                      uiOutput("kslider_rmsip"),
                      actionButton("setk_rmsip", "Auto set number of K groups",
                                   icon=icon("cogs"))
                      ),
                    
                    checkboxInput('show_options3', 'More options', value=FALSE),
                    downloadButton('nmaplot2pdf', "Download Plot PDF")
                    )
                  ),
           column(8,
                  conditionalPanel(
                    condition='output.nmaIsDone',
                    plotOutput("nma_fluctplot")
                    )
                  )
           ),

         conditionalPanel(
           condition = "input.show_options3 == true",
           fluidRow(
             column(3,
                    wellPanel(
                      selectInput("hclustMethod_rmsip", label="Clustering method", 
                                  choices=list(
                                    "single"="single","complete"="complete","average"="average",
                                    "mcquitty"="mcquitty","median"="median","centroid"="centroid",
                                    "ward.D"="ward.D","ward.D2"="ward.D2"
                                    ),selected="ward.D2"), 
                      
                      numericInput("minDistance_rmsip","Minimum branching gap", value = 0.1, step = 0.2)
                      
                      )
                    ),
             
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
           hr(),
           column(6,
                  conditionalPanel(
                    condition = "input.rm_gaps == true && output.nmaIsDone",
                    h4("RMSIP heatmap"),
                    plotOutput("rmsip_heatmap2")
                    )
                  ),
           column(6,
                  conditionalPanel(
                    condition = "output.nmaIsDone",
                    h4("RMSD heatmap"),
                    plotOutput("rmsd_heatmap2")
                    )
                  )
           ),

          fluidRow(
            conditionalPanel(
              condition = "output.nmaIsDone",
              column(3,
                     checkboxInput('show_options4', 'More options', value=FALSE)
                     ),
              column(3,
                     downloadButton('nma_rmsip_heatmap2pdf', "Download Plot PDF")
                     ),
              column(3),
              column(3,
                     downloadButton('nma_rmsd_heatmap2pdf', "Download Plot PDF")
                     )
              )
           ),

         conditionalPanel(
           condition = "input.show_options4== true",
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
           hr(),
           column(4,
                  wellPanel(
                    h4('Cluster dendrogram'),
                    checkboxInput('show_confplot2', 'Show PC conformer plot', value=FALSE),
                    checkboxInput('show_options5', 'More options', value=FALSE),

                    #radioButtons("group_by", "Cluster by",
                    #             c("RMSD" = "rmsd",
                    #               "RMSIP" = "rmsip",
                    #               "PC distance" = "pc_space"
                    #               ),
                    #             ##"bhat" = "bhat"),
                    #             inline=TRUE),

                    downloadButton('nmadendrogram2pdf', "Download Plot PDF")
                    )
                  ),
           column(8,
                  conditionalPanel(
                    condition = "output.nmaIsDone",
                    plotOutput("dendrogram2")
                    )
                  )
           ),

         conditionalPanel(
           condition = "input.show_options5 == true",
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
           )
  
  
         ### PCA vs NMA (RMSIP)
  #fluidRow(
  #         hr(),
  #         column(4,
  #                wellPanel(
  #                  bsPopover("popnma4",
  #                            "Comparison with PCA", 
  #                            "The normal mode vectors from one PDB structure can be compared with the eigenvectors obtained from principal component analysis (from the PCA tab). The values provided in the table to the right correspond to the overlap (i.e. the dot product) between the mode vectors. The RMSIP value provides an overall score for the similarity between the principal components and the normal mode vectors. Orthogonal vectors give the score of 0, while identical vectors give a score of 1. ", 
  #                            placement = "right", trigger = "hover",
  #                            options = list(container = "body")),
                    
  #                  tags$div(id = "popnma4", icon("question-circle"),
  #                           style = "position: absolute; right: 25px; top: 5px;"
  #                           ),
                    
  #                  h4('PCA vs NMA comparison'),
  #                  uiOutput('struct_dropdown3')
  #                  )
  #                )

           #column(4,
                  #conditionalPanel(
                  #  condition = "output.nmaIsDone",
                  #  plotOutput("rmsip_plot3")
                  #  )
           #       ),

           #column(4,
                  #conditionalPanel(
                  #  condition = "output.nmaIsDone",
                  #  DT::dataTableOutput("rmsip_table")
                  #  )
           #       )
   #        )


         )
