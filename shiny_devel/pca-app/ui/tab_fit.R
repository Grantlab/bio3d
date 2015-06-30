tabPanel("3. FIT", icon=icon("arrow-right"),
         tags$style(type="text/css", "body {padding-top: 80px;}"),

         actionButton3("about_fittab", "About this tab", icon=icon("comment"), cl="btn btn-warn btn-input action-button", style = "position: fixed; top: 14px; right: 16px; z-index: 2000;"),
         
         bsModal("modal_fit", "Structure superposition", "about_fittab", size = "large", 
                 content=tags$div(
                   p(HTML("In the FIT tab all structures are superimposed based on the sequence alignment in the ALIGN tab. The algorithm performs iterated rounds of structural superposition to identify the most invariant region in an aligned set of protein structures. Basic structural analysis in this tab entails pair-wise structural deviations (RMSD), fluctuation analysis (RMSF), and structure visualization.")),

                   img(src="./images/4q21-core.png", width=400, style="display: block; margin-left: auto; margin-right: auto;"), 
                   
                   p(HTML("The fitting algorithm attempts to iteratively refine the initial structural superposition determined from the multiple alignment. This involves iterated rounds of superposition, where at each round the position(s) displaying the largest differences is(are) excluded from the dataset. The identified core can be visualized either in the browser, or in PyMOL by downloading the PyMOL state (.pse) file. A summary of the core can be found in the bottom of the page."))
                   
                                    
                   ) 
                 ),
         
         
         fluidRow(
           column(4,
                  wellPanel(
                    
                    h4("Initial structure analysis"),

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

                    checkboxInput('show_options', 'More options', value=FALSE),

                    conditionalPanel(
                      condition = "input.str_plot == 'heatmap'",
                      downloadButton('rmsd_heatmap2pdf', "Download PDF")
                      ),

                    conditionalPanel(
                      condition = "input.str_plot == 'dendrogram'",
                      downloadButton('rmsd_dendrogram2pdf', "Download PDF")
                      ),

                    conditionalPanel(
                      condition = "input.str_plot == 'rmsf'",
                      downloadButton('rmsf2pdf', "Download PDF")
                      ),

                    conditionalPanel(
                      condition = "input.str_plot == 'hist'",
                      downloadButton('rmsd_hist2pdf', "Download PDF")
                      ),
                    downloadButton('rmsdZIP', "Download RMSD matrix")
                    )
                  ),

           column(8,

             conditionalPanel(
               condition = "input.str_plot == 'heatmap'",
               plotOutput("rmsd_heatmap")
               ##iplotCorr_output("qtl_rmsd_heatmap")
               ),

             conditionalPanel(
               condition = "input.str_plot == 'dendrogram'",
               plotOutput("rmsd_dendrogram")
               ),

             conditionalPanel(
               condition = "input.str_plot == 'rmsf'",
               plotOutput("rmsf_plot")
               ),

             conditionalPanel(
               condition = "input.str_plot == 'hist'",
               plotOutput("rmsd_hist")
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
                    br(),
                    ##downloadButton('pdbsRData', "Download PDBs RData"),
                    downloadButton('pdbsZIP', "Download Aligned PDBs"),
                    downloadButton('pdbs2pymol', "Download PyMOL session file")
                    )
                  ),

           column(8,
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
                    style="background: #FFFFFF;",
                    uiOutput("reference_selector")
                    )
                  )
           ),
         
         fluidRow(
           column(4,
                  wellPanel(
                    style="background: #FFFFFF;",

                    h4("Summary of invariant core"),
                    dataTableOutput("print_core")
                    )
                  ),

           column(4,
                  wellPanel(
                    style="background: #FFFFFF;",

                    h4("RMSD summary"),
                    ##uiOutput("reference_selector"),
                    dataTableOutput("rmsd_table")
                    )
                  ),

           column(4,
                  wellPanel(
                    style="background: #FFFFFF;",
                    
                    h4("Cluster representatives"),
                    verbatimTextOutput("representatives")
                    )
                  )
           )
         )
