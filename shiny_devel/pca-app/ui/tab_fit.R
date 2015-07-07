tabPanel("3. FIT", icon=icon("arrow-right"),
         tags$style(type="text/css", "body {padding-top: 80px;}"),

         actionButton3("about_fittab", "About this tab", icon=icon("comment"), cl="btn btn-warn btn-input action-button", style = "position: fixed; top: 14px; right: 16px; z-index: 2000;"),
         
         bsModal("modal_fit", "Structure superposition", "about_fittab", size = "large", 
                 content=tags$div(
                   p(HTML("In the FIT tab all structures are superimposed based on the sequence alignment obtained in the ALIGN tab. By default, the structures are superimposed on the identified <i>invariant core</i>. Alternatively, a conventional superposition based on all protein residues can be performed.")),

                   p(HTML("Basic structural analysis in this tab entails the calculation of all pair-wise structural deviations (RMSD), fluctuation analysis (RMSF), and structure visualization. Based on the pair-wise RMSD values the structures can be clustered into groups of similar conformations. The results of the RMSD and clustering analyses are available as a heatmap, clustering dendrogram, and RMSD histogram.")),
                   
                   p(HTML("The algorithm for identifying the invariant core attempts to iteratively refine the initial structural superposition determined from the multiple alignment. This involves iterated rounds of superposition, where at each round the position(s) displaying the largest differences is(are) excluded from the dataset. The identified core can be visualized either in the browser, or in PyMOL by downloading the PyMOL state (.pse) file. A summary of the core can be found in the bottom of the page.")),
                   
                   img(src="./images/4q21-core.png", width=400, style="display: block; margin-left: auto; margin-right: auto;")

                                    
                   ) 
                 ),
         
         
         fluidRow(
           column(4,
                  wellPanel(
                    bsPopover("popfit1",
                              "Structure analysis",
                              "All structures are (by default) superimposed on the <b>invariant core</b> - a region with low structural variability within the ensemble. </br></br>A hierarchical cluster analysis is performed on all pair-wise root mean square deviation (RMSD) values. Results can be visualized as a clustering <b>dendrogram</b>, alternatively in combination with a <b>heatmap</b>. </br></br><b>RMSF</b>: plot the residue-wise root mean square fluctuations (RMSF). </br></br> <b>RMSD histogram</b>: display the histogram of the pair-wise RMSD values.",
                              placement = "right", trigger = "hover",
                              options = list(container = "body")),
                    
                    tags$div(id = "popfit1", icon("question-circle"),
                             style = "position: absolute; right: 25px; top: 5px;"
                             ),
                    

                    
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
                    bsPopover("popfit2",
                              "Viewing options",
                              "Visualize the superimposed PDB structures by toggeling the <b>Show PDBs</b> checkbox. The <i>Download</i> buttons allow visualizing the structures in an external program such as PyMOL or VMD. <br><br>Coloring options include <b>Cluster ID</b>: colors the structures according to their cluster membership; <b>Invariant core</b>: by the region with low structural variation within the ensemble; and <b>Gap regions</b>: residues belonging to columns in the alignment (see ALIGN tab) which contain one or more gaps. ",
                              placement = "right", trigger = "hover",
                              options = list(container = "body")),
                    
                    tags$div(id = "popfit2", icon("question-circle"),
                             style = "position: absolute; right: 25px; top: 5px;"
                             ),
                    
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

                    bsPopover("popfit3",
                              "Invariant core",
                              "This panel provides a summary of which amino acid residues belongs to the <i>invariant core</i>. The residue numbering correspond to the PDB ID chosen in the <b>Reference PDB id</b> drop down menu.",
                              placement = "right", trigger = "hover",
                              options = list(container = "body")),
                    
                    tags$div(id = "popfit3", icon("question-circle"),
                             style = "position: absolute; right: 25px; top: 5px;"
                             ),
                    

                    h4("Summary of invariant core"),
                    dataTableOutput("print_core")
                    )
                  ),

           column(4,
                  wellPanel(
                    style="background: #FFFFFF;",

                    bsPopover("popfit4",
                              "RMSD",
                              "Structural deviation (root mean square deviation; RMSD) between a chosen PDB structure (see drop down meny <b>Reference PDB id</b>), and all other PDBs in the ensemble. The RMSD values are calculated based on all C-alpha atoms in the two proteins.",
                              placement = "right", trigger = "hover",
                              options = list(container = "body")),
                    
                    tags$div(id = "popfit4", icon("question-circle"),
                             style = "position: absolute; right: 25px; top: 5px;"
                             ),
                    
                    
                    h4("RMSD summary"),
                    ##uiOutput("reference_selector"),
                    dataTableOutput("rmsd_table")
                    )
                  ),

           column(4,
                  wellPanel(
                    style="background: #FFFFFF;",

                    bsPopover("popfit5",
                              "Cluster representatives",
                              "One structure from each cluster with the minimal distance to all the other cluster members is defined to be the <i>cluster representative</i>.",
                              placement = "left", trigger = "hover",
                              options = list(container = "body")),
                    
                    tags$div(id = "popfit5", icon("question-circle"),
                             style = "position: absolute; right: 25px; top: 5px;"
                             ),
                    
                    h4("Cluster representatives"),
                    dataTableOutput("representatives")
                    )
                  )
           )
         )
