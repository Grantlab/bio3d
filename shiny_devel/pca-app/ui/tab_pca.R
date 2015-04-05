 
tabPanel("4. PCA",
         
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
             
             checkboxInput("screeplot", "Scree plot", value=TRUE),
             
             checkboxInput("labelplot", "Label plot", value=FALSE),

             conditionalPanel(
               condition = "input.labelplot == true",
               sliderInput("offset", "label offset",
                           min = 0, max = 2, value = 0.5, step=0.1),
               uiOutput("checkboxgroup_label_ids")
               )
            )
                  ),
           
           column(8, 
           
             conditionalPanel(
               condition = "input.plot_type == 'normal'",
               plotOutput("pca_plot1")
               ),

             conditionalPanel(
               condition = "input.plot_type == 'fancy'",
               showOutput("pca_plot2","dimple")
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
