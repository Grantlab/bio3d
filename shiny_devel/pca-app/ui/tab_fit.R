tabPanel("3. FIT", icon=icon("arrow-right"),
         fluidRow(
           column(4,
                  wellPanel(
             h4("Initial structure analysis"),
             
             radioButtons("fit_type", "Superimpose to",
                          c("invariant core" = "core",
                            "all c-alpha atoms" = "full"),
                          inline=TRUE),

             
             radioButtons("str_plot", "Plot options",
                          c("heatmap" = "heatmap",
                            "dendrogram" = "dendrogram",
                            "rmsf" = "rmsf",
                            "hist" = "hist"),
                          inline=TRUE),
             
             sliderInput("cex", "cex",
                         min = 0.1, max = 3, value = 1, step=0.1),
         
             sliderInput("clusters", "Cluster by pairwise RMSD",
                         min = 1, max = 10, value = 3, step=1),

             hr(),
             downloadButton('pdbsZIP', "Download Aligned PDBs"),
             downloadButton('rmsdZIP', "Download RMSD matrix")

                    )
                  ),
           
           column(8,
                  
             conditionalPanel(
               condition = "input.str_plot == 'heatmap'",
               plotOutput("rmsd_heatmap"),
               downloadButton('rmsd_heatmap2pdf', "Download PDF")
               ),
             
             conditionalPanel(
               condition = "input.str_plot == 'dendrogram'",
               plotOutput("rmsd_dendrogram"),
               downloadButton('rmsd_dendrogram2pdf', "Download PDF")
               ),
             
             conditionalPanel(
               condition = "input.str_plot == 'rmsf'",
               plotOutput("rmsf_plot"),
               downloadButton('rmsf2pdf', "Download PDF")
               ),
             
             conditionalPanel(
               condition = "input.str_plot == 'hist'",
               plotOutput("rmsd_hist"),
               downloadButton('rmsd_hist2pdf', "Download PDF")
               )
             )
           ),

         fluidRow(
           column(4,
                  h4("Summary of invariant core"),
                  dataTableOutput("print_core")
                  ),

           column(4,
                  h4("RMSD summary"),
                  uiOutput("reference_selector"),
                  dataTableOutput("rmsd_table")
                  )
           )

         
         )
