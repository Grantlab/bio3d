tabPanel(
  "4. PCA", icon=icon("arrow-right"),
  tags$style(type="text/css", "body {padding-top: 80px;}"),

  
  actionButton3("about_pcatab", "About this tab", icon=icon("comment"), cl="btn btn-warn btn-input action-button", style = "position: fixed; top: 14px; right: 16px; z-index: 2000;"),
  
  bsModal("modal_pca", "Principal Component Analysis", "about_pcatab", size = "large", 
          content=tags$div(
            
            p(HTML("In the PCA tab, a principal component analysis (PCA) is performed based on the coordinates of the superimposed structures in the FIT tab. PCA is an effective approach to capture and characterize inter-conformer relationships. The motions described by the principal components (PCs) can be visualized either in the browser, or by downloading a PDB trajectory file or a PyMOL session (.pse) containing a vector field representation.  ")),
            
            img(src="./images/aaah-pc1.gif", width=500, style="display: block; margin-left: auto; margin-right: auto;")                      
            ),
          
          p(HTML("The figure above shows the first principal component obtained from the 32 PDB structures of the enzyme family <i>aromatic amino acid hydroxylases</i>. The analysis reveals the presence of a sub-domain within the catalytic domain which obtains a more compact conformation upon substrate binding.")),

          br(),
          img(src="./images/pca-conformerplot.png", style="display: block; margin-left: auto; margin-right: auto;"),
          
          p(HTML("The figure above shows the conformer plot of 90 <i>E. coli</i> DHFR structures. The PCA reveals that the ensemble can be divided into three major groups along their first two principal components. These conformers display either a closed, occluded,oranopenconformationof two active site loops. The conformer plot displays a low-dimensional representation of conformational variability described by two selected PCs (typically PC-1 and PC-2). "))
          ),
  

  ##- PC Visualization 
  fluidRow(
    column(4,
           wellPanel(
             bsPopover("poppca1",
                       "Visualization",
                       "Visualize the motions described by the individual principal components by toggeling the the <b>Show PC Trajectory</b> checkbox. <br><br>Download buttons enable visualization of the motions described by the principal component in external viewers such as PyMOL or VMD.",
                       placement = "right", trigger = "hover",
                       options = list(container = "body")),
             
             tags$div(id = "poppca1", icon("question-circle"),
                      style = "position: absolute; right: 25px; top: 5px;"
                      ),
             
             
             h4('Principal Component Visualization'),
             checkboxInput('show_trj', 'Show PC Trajectory', value=FALSE),

             selectInput('viewPC', 'Choose Principal Component:', choices=c(1:10)),

             sliderInput("mag1", "Magnification factor:",
                         min = 1, max = 10, value = 1),
             
             radioButtons('viewColor', label='Structure color',
                             choices=list(
                               'Amalgam' = 'amalgam',
                               'Magnitude'='mag',
                               'By Frame (blue->gray->red)'='default'
                               ),
                             selected='amalgam'
                          ),
             
             radioButtons('viewBGcolor', label='Background color',
                          choices=list('Black'='black', 'White'='white'),
                          selected='white'),
             br(),
             actionButton('viewUpdate', label='Refresh', icon=icon('undo')),
             downloadButton('pctraj', label='Download PDB Trajectory'),
             downloadButton('pca2pymol', "Download PyMOL session file")

             )
           ),

    column(8,
           conditionalPanel(
             condition='input.show_trj == true',
             webGLOutput('pcaWebGL')
             )
           )
  ),


  ##- Conformer plots
  fluidRow(
    column(4,
           wellPanel(
             bsPopover("poppca2",
                       "Conformer plot",
                       "A conformer plot is a low-dimensional representation of the conformational variability within the ensemble of PDB structures. The plot is obtained by projecting individual structures onto two or three selected principal components. ",
                       placement = "right", trigger = "hover",
                       options = list(container = "body")),
             
             tags$div(id = "poppca2", icon("question-circle"),
                      style = "position: absolute; right: 25px; top: 5px;"
                      ),
                    
                     
        h4("Conformer plot"),
             popSelectInput('pcx', 'PC on X-axis:', choices=c(1:10), selected=1),
        selectInput('pcy', 'PC on Y-axis:', choices=c(1:10), selected=2),

        conditionalPanel(
          condition = "input.plot_type == '3dscatter1' || input.plot_type == '3dscatter2'",
          selectInput('pcz', 'PC on Z-axis:', choices=c(1:10), selected=3)
          ),

        radioButtons('cluster_by', label='Cluster by',
                          choices=list('RMSD'='rmsd', 'PC subspace'='pc_space'),
                     selected='rmsd', inline=TRUE),

        conditionalPanel(
          condition = "input.cluster_by == 'pc_space'",
          sliderInput("clust_npcs", "PCs in subspace",
                      min = 1, max = 10, value = 2)
          ),

        sliderInput("nclust", "Clusters",
          min = 1, max = 10, value = 3),

        radioButtons("plot_type", "",
                     c("2D scatter" = "normal",
                       "3D scatter" = "3dscatter1",
                       "Interactive" = "fancy"
                       ##"3D scatter (three-js)" = "3dscatter2",
                       ),
                     inline=TRUE),

        conditionalPanel(
          condition = "input.plot_type == 'normal'",
          checkboxInput("show_options1", "More options", value=FALSE)
          ),

        conditionalPanel(
          condition = "input.plot_type == '3dscatter2'",
          selectInput("renderer", label="Rendering method",
                      choices = list("Auto"="auto", "Canvas"="canvas", "WebGL"="webgl"),
                      selected = 1),
          checkboxInput("grid", label = "Grid", value = TRUE)
          )
        )
           ),

    conditionalPanel(
      condition = "input.plot_type == '3dscatter1'",
      column(width=8,
             webGLOutput("scatterplot3d_webgl")
             )
      ),

    conditionalPanel(
      condition = "input.plot_type == '3dscatter2'",
      column(width=8,
             scatterplotThreeOutput("scatterplot3d_rthreejs")
             )
      ),

    conditionalPanel(
      condition = "input.plot_type != '3dscatter'",

      column(width=4,
             conditionalPanel(
               condition = "input.plot_type == 'normal'",
               plotOutput("pca_plot1_conf")
               ),

             conditionalPanel(
               condition = "input.plot_type == 'fancy'",
               showOutput("pca_plot2_conf","dimple")
               )
             ),

      column(4,
             conditionalPanel(
               condition = "input.plot_type == 'normal'",
               plotOutput("pca_plot1_scree")
               ),
             conditionalPanel(
               condition = "input.plot_type == 'fancy'",
               showOutput("pca_plot2_scree","nvd3")
               )
             )
      )

    ),

 
    



  ##- Additional plotting options for conformer plot
  fluidRow(
    conditionalPanel(
    condition = "input.show_options1 == true && input.plot_type == 'normal'",
      column(3,
             wellPanel(
               sliderInput("cex_points", "Point size",
                           min = 0.1, max = 3, value = 1, step=0.1),
               sliderInput("inner_margin", "Scale axes",
                           min = 1, max = 2, value = 1.2, step=0.1)
               )
             ),

      column(3,
             wellPanel(
               checkboxInput("labelplot", "Label plot", value=FALSE),

               conditionalPanel(
                 condition = "input.labelplot == true",

                 checkboxInput("distribute_labels", "Distribute labels", value=FALSE),
                 sliderInput("cex_labels", "Label size",
                             min = 0.1, max = 3, value = 1, step=0.1),

                 sliderInput("offset", "label offset",
                             min = 0, max = 2, value = 0.5, step=0.1)
                 )
               )
             ),

      column(6,
             wellPanel(
               conditionalPanel(
                 condition = "input.labelplot == true",

                 checkboxInput("toggle_all", "Toggle all", TRUE),
                 uiOutput("checkboxgroup_label_ids")
                 )
               )
             )
    )
    ),

  fluidRow(
    column(12,
           wellPanel(
             h4("Hits annotation"),
             helpText("Highlight structures in conformer plot by clicking their entries in the below table (only for plot type '2D Scatter')."),
             DT::dataTableOutput("pdbs_table")
             )
           )
    ),
           

  
  ##- PC Loadings - fluctuations
  fluidRow(
    column(4,
           wellPanel(
             bsPopover("poppca3",
                       "Residue contributions",
                       "Residue-wise contributions to the chosen principal component.",
                       placement = "right", trigger = "hover",
                       options = list(container = "body")),
             
             tags$div(id = "poppca3", icon("question-circle"),
                      style = "position: absolute; right: 25px; top: 5px;"
                      ),

             
             h4('Residue contributions'),
             selectInput('loadings_pc', 'Choose Principal Component:',
                         choices=c(1:10), selected=1, multiple=TRUE),
             ##checkboxInput("toggle_rmsf1", "Show RMSF", TRUE),
             checkboxInput("toggle_rmsf1", "Show RMSF", TRUE),
             checkboxInput("multiplot", "Multiline plot", FALSE),
             checkboxInput("spread_pcload", "Spread lines", FALSE),
             downloadButton('pcloadings2pdf', label='Download PDF'),
             checkboxInput("show_options2", "More options", value=FALSE)
             )
           ),
    
    column(8,
           plotOutput("loadings_plot")
           )
    ),
  
  conditionalPanel(
    condition = "input.show_options2 == true",
    
    fluidRow(
      column(4,
             wellPanel(
               sliderInput("width_pcload", "Width",
                           min = 4, max = 12, value = 8, step=0.5),
               sliderInput("height_pcload", "Height",
                           min = 4, max = 12, value = 5, step=0.5)
               )
             )
      )
    )

  ##- Annotation datatable
  #fluidRow(
  #  column(12,
  #         wellPanel(
  #           DT::dataTableOutput("pdbs_table")
  #           )
  #         )
  #  )
  )
