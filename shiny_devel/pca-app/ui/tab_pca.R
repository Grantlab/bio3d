tabPanel("4. PCA", icon=icon("arrow-right"), 
  tags$style(type="text/css", "body {padding-top: 80px;}"), 
  fluidRow(
    column(4, 
      wellPanel(
        helpText("TODO: trajectory download"),
        textInput("pcx", "PC on X-axis", value=1),
        textInput("pcy", "PC on Y-axis", value=2),
        sliderInput("nclust", "Cluster by pairwise RMSD", 
          min = 1, max = 10, value = 1),
        radioButtons("plot_type", "Plot type:",
                   c("Normal" = "normal",
                     "Fancy" = "fancy"),
                 inline=TRUE), 
        checkboxInput("labelplot", "Label plot", value=FALSE), 

        conditionalPanel(
          condition = "input.labelplot == true", 
          sliderInput("offset", "label offset", 
            min = 0, max = 2, value = 0.5, step=0.1), 
          uiOutput("checkboxgroup_label_ids")
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
    column(12, 
      wellPanel(
        dataTableOutput("pdbs_table")
      )
    )
  )
)
