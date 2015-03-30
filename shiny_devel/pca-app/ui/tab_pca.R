 
tabPanel("3. PCA",

         sidebarLayout(
           sidebarPanel(
             checkboxInput("screeplot", "Scree plot:", value=TRUE),
             textInput("pcx", "PC on X-axis", value=1),
             textInput("pcy", "PC on Y-axis", value=2),
                         
                         
             sliderInput("nclust", "Cluster by pairwise RMSD",
                         min = 1, max = 10, value = 1),

             checkboxInput("labelplot", "Label plot", value=FALSE),

             conditionalPanel(
               condition = "input.labelplot == true",
               sliderInput("offset", "label offset",
                           min = 0, max = 2, value = 0.5, step=0.1),
               uiOutput("checkboxgroup_label_ids")
               )
             ),
           mainPanel(
             showOutput("pca_plot","dimple"),
             dataTableOutput("myTable")
             )
           )
         )
