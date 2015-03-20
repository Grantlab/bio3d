library(shiny)
library(bio3d)
library(shinyRGL)


## Could open some saved results as an RData file here for first time display 
## Calculations would then only run once single SUBMIT button was pressed.

shinyUI(fluidPage(

  title = "Bio3D NMA",

  h1("Bio3D NMA"),  
  h2("Input"),  
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
  h2("Results"),

  ##-A. Fluctuations Panel
  conditionalPanel(
    condition = "input.fluxs1 == true",

    h3("Fluctuations"),
    fluidRow(
      plotOutput("fluct_plot"),
      column(3,
             checkboxInput('show_options', 'More options', value=FALSE)
             ),
      column(3,
             checkboxInput('fluxs2', 'Show B-factors', value=FALSE)
             ),
      column(3),
      column(3,
             downloadButton('plot1pdf', "Download Plot PDF")
             )
      ),

    conditionalPanel(
      condition = "input.show_options == true",
      fluidRow(
        column(3,
               h4("Plot options"),
               checkboxInput('fluxs3', 'Show', value=TRUE),
               selectInput("typ1", "Type:",
                           c("hist" = "h",
                             "lines" = "l",
                             "points" = "p",
                             "both" = "b")),
               sliderInput("lty1", "Line type:",
                           min = 1, max = 6, value = 1),
               sliderInput("lwd1", "Line width:",
                           min = 0.1, max = 2, value = 1, step=0.1)
               ),

        column(3,
               h4("Plot options"),
               sliderInput("pch1", "Point type:",
                           min = 1, max = 25, value = 1),
               sliderInput("cex1", "Point size:",
                           min = 0.1, max = 2, value = 1, step=0.1),
               sliderInput("col1", "Color:",
                           min = 1, max = 8, value = 1)
               ),
        
        column(3,
               h4("B-factors"),
               checkboxInput('fluxs2', 'Show B-factors', value=FALSE),
               selectInput("typ2", "Type:",
                           c("lines" = "l",
                             "points" = "p",
                             "both" = "b",
                             "hist" = "h")),
               
               sliderInput("lty2", "Line type:",
                           min = 1, max = 6, value = 1),
               sliderInput("lwd2", "Line width:",
                           min = 0.1, max = 2, value = 1, step=0.1)
               ),

        column(3,
               h4("B-factors"),
               sliderInput("pch2", "Point type:",
                           min = 1, max = 25, value = 1),
               sliderInput("cex2", "Point size:",
                           min = 0.1, max = 2, value = 1, step=0.1),
               sliderInput("col2", "Color:",
                           min = 1, max = 8, value = 4)
               )
        ) ## end fluid row
      ) ## end conditional panel
    
   
    ),
  
  br(),br(),





  ##-C. DCCM Panel
  conditionalPanel(
    condition = "input.cijs == true",
    h3("Correlations"),
    fluidRow(
      plotOutput("dccm_plot")
      ),
    fluidRow(
      column(3,
             checkboxInput('contourplot', 'Contourplot', value=TRUE)
             ),
      column(3,
             checkboxInput('sse', 'Show SSE', value=TRUE)
             ),
      column(3,
             checkboxInput('colorkey', 'Colorkey', value=TRUE)
             ),
      column(3,
             downloadButton('plot2pdf', "Download Plot PDF")
             )
      )
    ),
  br(),br(),

  

  hr(),
  ##-D. Log Report Panel
  conditionalPanel(
    condition = "input.log == true",
    column(6,
           h4("PDB Summary"),
           verbatimTextOutput("pdbSummary")
           ),
    column(6,
           h4("Calculation Summary Log"),
           verbatimTextOutput("modeSummary")
           )
    
    )
))
