library(shiny)
library(bio3d)
library(shinyRGL)


## Could open some saved results as an RData file here for first time display 
## Calculations would then only run once single SUBMIT button was pressed.

shinyUI(fluidPage(

  title = "Dummy Example (3 Column Input)",

  h3("Input"),  
  hr(),

  fluidRow(
    column(4,
      wellPanel(
        h4("1. PDB Input Selection"),
        tags$hr(),

        ##-PDB input (moved to server.R)
        uiOutput('resetable_pdb_input'),
        #textInput("pdbid", label="Enter RCSB PDB code/ID:", value = "4Q21"),

        ##- Chain selection
        h5("Detected chain IDs:"),
        verbatimTextOutput("chains1"),
        
        checkboxInput("limit", "Limit calculation to a subset of chains?"),
        helpText("Note: Use this option to exclude particular chains form further consideration."),
        
        conditionalPanel(
          condition = "input.limit == true",
          uiOutput("chains2"),
          helpText("Note: Only selected chains will be analyzed.")
        )

        ##-TODO Have this panel auto-update only!
        #,actionButton("pdbaction", "Fetch PDB")
        #,submitButton("Submit")
        ,actionButton("reset_pdb_input", "Reset PDB input")
      )
    ),

    column(4,
      wellPanel(
        h4("2. NMA Calculation"),
        tags$hr(),

        uiOutput('resetable_nma_input'),
        
        helpText("Note: Some help text here. "),
 
        actionButton("reset_nma_input", "Reset NMA inputs")
        ##submitButton("Submit") - will effect first PDB input panel also
      )
    ),
    
    column(4,
      wellPanel(
        h4("3. Result Visualization"),
        tags$hr(),
        checkboxInput('fluxs1', 'Fluctuations', value=TRUE),
        #checkboxInput('fluxs2', 'Fluctuations with B-factors', value=FALSE),
        checkboxInput('cijs', 'Correlations', value=FALSE),
        checkboxInput('mktrj', 'Mode Displacements Trajectory', value=FALSE),
        checkboxInput('log', 'Calculation Summary Log', value=TRUE)
      )
    )
  ),



  #########################
  ##-- Results Section --##
  #########################
  hr(),
  h3("Results"),

  ##-A. Fluctuations Panel
  conditionalPanel(
    condition = "input.fluxs1 == true",
    plotOutput("fluctPlot"),

    ## TODO: Add more plot customization options for typ, lwd, col.
    checkboxInput('fluxs2', 'Show B-factors', value=FALSE),
 
    ## TODO: Download button (NOT YET WORKING)
    actionButton("plot1pdf", label = "Download Plot PDF")
  ),
  br(),br(),



  ##-B. 3D Molecular Viewer Panel
  conditionalPanel(
    condition = "input.mktrj == true",
    column(4,
      webGLOutput("myWebGL") #,
    ),
    column(3,
      wellPanel(
        h4("Mode Trajectory Viewing Options"),
        selectInput("viewMode", "Choose a mode number:", choices = c(7:20)),
        #radioButtons("viewMode", label = "Structure representation",
        #  choices = list("default" = "default", 
        #  "All"="all", 
        #  "Calpha" = "calpha", 
        #  "Backbone" = "back"), 
        #  selected = "default"),

        ## TODO: Some bug fixing in view.xyz() required before these wok properly
        radioButtons("viewColor", label = "Structure color",
          choices = list("By Frame" = "default", 
                         "Magnitude"="all", 
                         "Gray" = "gray"), 
          selected = "default"),

        radioButtons("viewBGcolor", label = "Background Color",
          choices = list("Black" = "black", "White"="white"), 
          selected = "black"),

        br(),
        actionButton("viewUpdate", label = "Refresh"),
        actionButton("trj2pdb", label = "Download PDB Trajectory")
      )
    )
  ), 



  ##-C. DCCM Panel
  conditionalPanel(
    condition = "input.cijs == true",
    column(5,
      plotOutput("dccmPlot"),
      actionButton("plot2pdf", label = "Download Plot PDF")
    )
  ),
  br(),br(),


  ##-D. Log Report Panel
  conditionalPanel(
    condition = "input.log == true",
    column(12,
    h4("Calculation Summary Log"),
    verbatimTextOutput("pdbSummary"),
    verbatimTextOutput("modeSummary")# ),
    )
  )
))
