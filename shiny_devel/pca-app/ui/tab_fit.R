tabPanel("3. FIT", icon=icon("arrow-right"),
         tags$style(type="text/css", "body {padding-top: 80px;}"),

         fluidRow(
           column(4,
                  wellPanel(
                    h4("Initial structure analysis"),
                    hr(),

                    radioButtons("fit_type", "Superimpose to",
                                 c("Invariant core" = "core",
                                   "All c-alpha atoms" = "full"),
                                 inline=TRUE),

                    radioButtons("str_plot", "Plot options",
                                 c("Heatmap" = "heatmap",
                                   "Dendrogram" = "dendrogram",
                                   "RMSF" = "rmsf",
                                   "RMSD Histogram" = "hist"),
                                 inline=TRUE),

                    sliderInput("clusters", "Cluster by pairwise RMSD",
                                min = 1, max = 10, value = 3, step=1),

                    checkboxInput('show_options', 'More options', value=FALSE)
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

         conditionalPanel(
           condition = "input.show_options == true",
           fluidRow(
             column(4,
                    wellPanel(
                    sliderInput("cex", "Label size",
                                min = 0.1, max = 3, value = 1, step=0.1),
                    sliderInput("margins", "Plot margins",
                                min = 3, max = 10, value = 5, step=1)
                    )
                    ),
             column(4,
                    wellPanel(
                    sliderInput("width", "width",
                                min = 4, max = 12, value = 7, step=0.5),
                    sliderInput("height", "height",
                                min = 4, max = 12, value = 7, step=0.5)
                    )
                    )
             )
           ),


         br(),
         fluidRow(
           column(4,
                  wellPanel(
                    h4('PDBs Viewing Options'),
                    checkboxInput('show_pdbs', 'Show PDBs', value=FALSE),

                    radioButtons('viewColor1', label='Structure color',
                                 choices=list(
                                   'By cluster ID'='cluster',
                                   'By structure ID'='struct',
                                   'Invariant core'='core',
                                   'Gap regions'='gaps'
                                   ),
                                 selected='cluster'),

                    radioButtons('viewBGcolor1', label='Background color',
                                 choices=list('Black'='black', 'White'='white'),
                                 selected='white'),
                    br(),
                    actionButton('viewUpdate1', label='Refresh', icon=icon('undo')),
                    downloadButton('pdbsZIP', "Download Aligned PDBs")
                    )
                  ),

           column(6,
                  conditionalPanel(
                    condition='input.show_pdbs == true',
                    webGLOutput('pdbsWebGL')
                    )
                  )
           ),


         hr(),

         fluidRow(
           column(4,
                  wellPanel(
                    h4("Summary of invariant core"),
                    dataTableOutput("print_core")
                    )
                  ),

           column(4,
                  wellPanel(
                    h4("RMSD summary"),
                    uiOutput("reference_selector"),
                    dataTableOutput("rmsd_table")
                    )
                  ),

           column(4,
                  wellPanel(
                    h4("Download"),
                    #downloadButton('pdbsZIP', "Download Aligned PDBs"),
                    hr(),
                    downloadButton('rmsdZIP', "Download RMSD matrix"),
                    hr()
                    ),

                  wellPanel(
                    h4("Cluster representatives"),
                    verbatimTextOutput("representatives")
                    )
                  )
           )
         )
