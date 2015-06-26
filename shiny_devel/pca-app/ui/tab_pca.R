tabPanel(
  "4. PCA", icon=icon("arrow-right"),
  tags$style(type="text/css", "body {padding-top: 80px;}"),

  
  modalBox(id="3", button_label = "Help ", icon = "question",
           heading="Principal component analysis",
           content = tags$div(
             HTML("<p>In this tab the collected structures are superimposed on each other either based on the <strong>identified invariant core</strong>, or on all C-alpha atoms. The invariant core is the region ...</p>"),
             
             p("In this panel you can perform simple structure analysis such as calculating all pair-wise RMSD values ... ")
             )
           ),
  

  ##- PC Visualization 
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


  ##- Conformer plots
  fluidRow(
    column(4,
      wellPanel(
        ##-- !! This needs to come just before the table that allows structure selection !! --##
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
             helpText("Highlight structures in conformer plot by clicking their entries in the below table (only for plot type '2D Scatter'."),
             DT::dataTableOutput("pdbs_table")
             )
           )
    ),
           

  
  ##- PC Loadings - fluctuations
  fluidRow(
    column(4,
           wellPanel(
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
