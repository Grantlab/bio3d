library(shiny)
library(bio3d)

## Could open some saved results here for first display 

shinyUI(fluidPage(

  title = "Dummy Example (3 Column Input)",

  h3("Input"),  
  hr(),

  fluidRow(
    column(4,
      wellPanel(
        h4("1. PDB Input Selection"),
        tags$hr(),

        ##-PDB input
        textInput("pdbid", label="Enter RCSB PDB code/ID:", value = "4Q21"),

        ##- Chain selection
        h5("Detected chain IDs:"),
        verbatimTextOutput("chains"),
      
        checkboxInput("limit", "Limit calculation to a subset of chains?"),
       helpText("Note: Use this option to exclude particular chains form further consideration."),

        conditionalPanel(
          condition = "input.limit == true",
          checkboxGroupInput("chain", label = "Limit to chain IDs:", 
          choices = list("A" = 1, "B" = 2, "C" = 3),
          selected = 1:3),
          ## N.B. should use choices list from "chainIDs" defined in server.R

        helpText("Note: Only selected chains will be analyzed. _NOT IMPLIMENTED_")

        )
      )
    ),

    column(4, #offset = 1,
      wellPanel(
      h4("2. NMA Calculation"),
      tags$hr(),

      ##-- N.B. Not currently wired up to server.R
       selectInput("forcefield", "Choose a forcefield:",
         choices = c("calpha", "anm", "pfanm", "sdenm")),


       sliderInput("cutoff", "Cutoff value:",
         min = 6, max = 50, value = 10),

       
      # radioButtons("mass", "Mass-weighting:",
      #   c("Yes", "Nope")),

       sliderInput("temp", "Temperature scaling:",
         min = 0, max = 350, value = 300),

       checkboxInput("mass", "Mass-weighting", value=TRUE),

      helpText("Note: Some help text here. _NOT IMPLIMENTED_")
       ##submitButton("Submit") - will effect first PDB input panel also
      )
    ),
    
    column(4,
      wellPanel(
      h4("3. Result Visualization"),
      tags$hr(),

      checkboxInput('fluxs', 'Fluctuations', value=TRUE),
      checkboxInput('fluxs', 'Fluctuations with B-factors', value=FALSE),
      checkboxInput('cijs', 'Correlations', value=FALSE),
      checkboxInput('mktrj', 'Mode Displacements Trajectory', value=FALSE),
      checkboxInput('log', 'Calculation Summary Log', value=TRUE)

      )
    )
  ),

  hr(),
  h3("Results"),

  plotOutput("fluctPlot"),

  ## Download button
  actionButton("plot1pdf", label = "Download Plot PDF"),

  br(),br(),

#  column(6,
#    h4("DCCM"),
#    verbatimTextOutput("pdbSummary") ),
#  column(6,
#    h4("Mode Trajectory"),
#    verbatimTextOutput("modeSummary"),# ),


  ##-- The below should be an optional display that is turned off by default...
  h5("Calculation Summary Log"),
  #verbatimTextOutput("pdbSummary"),
  verbatimTextOutput("modeSummary")# ),


))
