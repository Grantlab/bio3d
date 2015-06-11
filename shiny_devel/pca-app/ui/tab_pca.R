tabPanel("4. PCA", icon=icon("arrow-right"),
  tags$style(type="text/css", "body {padding-top: 80px;}"),

  fluidRow(
    column(4,
           wellPanel(
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
                          selected='amalgam'),
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

  fluidRow(
    column(4,
      wellPanel(
        h4("Conformer plot"),
        helpText("Two dimensional representation of conformational variability described by the two principal components ..."),

        ##textInput("pcx", "PC on X-axis", value=1),
        ##textInput("pcy", "PC on Y-axis", value=2),
        selectInput('pcx', 'PC on X-axis:', choices=c(1:10), selected=1),
        selectInput('pcy', 'PC on Y-axis:', choices=c(1:10), selected=2),

        conditionalPanel(
          condition = "input.plot_type == '3dscatter1' || input.plot_type == '3dscatter2'",
          selectInput('pcz', 'PC on Z-axis:', choices=c(1:10), selected=3)
          ),

        radioButtons('cluster_by', label='Cluster by',
                          choices=list('RMSD'='rmsd', 'PC subspace'='pc_space'),
                     selected='rmsd'),

        conditionalPanel(
          condition = "input.cluster_by == 'pc_space'",
          sliderInput("clust_npcs", "PCs in subspace",
                      min = 1, max = 10, value = 2)
          ),

        sliderInput("nclust", "Clusters",
          min = 1, max = 10, value = 3),

        radioButtons("plot_type", "",
                   c("2D scatter" = "normal",
                     "3D scatter (rgl)" = "3dscatter1",
                     "3D scatter (three-js)" = "3dscatter2",
                     "Interactive" = "fancy"
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

  fluidRow(
#    column(3,
#           wellPanel(
#             h4("Trajectory download"),
#             helpText("Trajectories of the first 5 principal components can be downloaded and visualized in your favorite visualization program (e.g. PyMOL or VMD)."),
#             downloadButton('pctrajZIP', "Download")
#             )
#           ),

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
           column(4,
                  wellPanel(
                    h4('Residue contributions'),
                    selectInput('loadings_pc', 'Choose Principal Component:',
                                choices=c(1:10), selected=1, multiple=TRUE),
                    checkboxInput("toggle_rmsf1", "Show RMSF", TRUE),
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
                      sliderInput("width-pcload", "Width",
                                  min = 4, max = 12, value = 7, step=0.5),
                      sliderInput("height-pcload", "Height",
                                  min = 4, max = 12, value = 7, step=0.5)
                      )
                    )
             )
           ),

         fluidRow(
           column(12,
                  wellPanel(
                    DT::dataTableOutput("pdbs_table")
                    )
                  )
           )
         )
