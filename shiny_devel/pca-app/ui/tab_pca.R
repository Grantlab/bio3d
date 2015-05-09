tabPanel("4. PCA", icon=icon("arrow-right"),
  tags$style(type="text/css", "body {padding-top: 80px;}"),
  fluidRow(
    column(4,
      wellPanel(
        h4("Conformer plot"),
        helpText("Two dimensional representation of conformational variability described by the two principal components ..."),
        textInput("pcx", "PC on X-axis", value=1),
        textInput("pcy", "PC on Y-axis", value=2),
        sliderInput("nclust", "Cluster by pairwise RMSD",
          min = 1, max = 10, value = 3),

        radioButtons("plot_type", "",
                   c("Normal" = "normal",
                     "Interactive" = "fancy"),
                     inline=TRUE),
        checkboxInput('show_trj', 'Show PC Trajectory', value=FALSE),
        conditionalPanel(
             condition = "input.plot_type == 'normal'",
             checkboxInput("show_options", "More options", value=FALSE)
        )
        )

    ),

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
    condition = "input.show_options == true && input.plot_type == 'normal'",
      column(3,
             wellPanel(
               sliderInput("cex_points", "Point size",
                           min = 0.1, max = 3, value = 1, step=0.1)
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

      column(3,
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


  conditionalPanel(
    condition='input.show_trj == true',
    column(3,
           wellPanel(
                     h4('PC Trajectory Viewing Options'),
                     selectInput('viewPC', 'Choose Principal Component:', choices=c(1:10)),
                     radioButtons('viewColor', label='Structure color',
                                  choices=list('Magnitude'='mag', 'By Frame (blue->gray->red)'='default'),
                                  selected='mag'),
                     radioButtons('viewBGcolor', label='Background color',
                                  choices=list('Black'='black', 'White'='white'),
                                  selected='white'),
                     br(),
                     actionButton('viewUpdate', label='Refresh', icon=icon('undo')),
                     downloadButton('pctraj', label='Download PDB Trajectory')
                     )
           ),
    column(4,
            webGLOutput('pcaWebGL')
          )

  ),
  fluidRow(
    column(12,
      wellPanel(
        dataTableOutput("pdbs_table")
      )
    )
  )
)
