tabPanel("NMA", icon=icon("home"),
         tags$style(type="text/css", "body {padding-top: 80px;}"),
         
  fluidRow(
    column(4,
           wellPanel(
             checkboxInput('show_pdb', 'View Input PDB', value=FALSE),

             conditionalPanel(
               condition = "input.show_pdb == true",
               h3("Input PDB"),
               webGLOutput('pdbWebGL')
               ),
             
             conditionalPanel(
               condition = "input.show_pdb == false",
               h3("Normal modes analysis with Bio3D"),
               p("Bio3D@web provides a ... "),
               p("Start by entering a PDB code of interest .. ."),
               img(src="geostas_250x182.png",
                   width=250, style="display: block; margin-left: auto; margin-right: auto;")
               )
             
             
             )
           ),
    
    
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
        h4("2. NMA parameters"),
        tags$hr(),

        uiOutput('resetable_nma_input'),
        
        helpText("Note: Cutoff applies only to 'ANM' and 'pfANM'. Recommended values are 15 and 50 Å, respectively."),
 
        actionButton("reset_nma_input", "Reset NMA inputs")
        )
    )
    
    #column(4,
    #  wellPanel(
    #    h4("3. Result Visualization"),
    #    tags$hr(),
    #    checkboxInput('fluxs1', 'Fluctuations', value=TRUE),
    #    #checkboxInput('fluxs2', 'Fluctuations with B-factors', value=FALSE),
    #    checkboxInput('cijs', 'Correlations', value=FALSE),
    #    checkboxInput('mktrj', 'Mode Displacements Trajectory', value=FALSE),
    #    checkboxInput('domains', 'Domain analysis', value=FALSE),
    #    checkboxInput('log', 'Calculation Summary Log', value=TRUE),
    #
    #    helpText("Scroll down to see results.")
    #  )
    #)
  ),



  #########################
  ##-- Results Section --##
  #########################
  h2("Results"),
  hr(),

  ##-A. Fluctuations Panel
  fluidRow(
    column(4,
           wellPanel(
             h4('Residue fluctuations'),
             checkboxInput('fluxs3', 'Show mode fluctuations', value=TRUE),
             checkboxInput('fluxs2', 'Show B-factors', value=FALSE),
             checkboxInput('show_options1', 'More options', value=FALSE),
             downloadButton('fluctplot2pdf', "Download Plot PDF")
             )
           ),
    column(8,
           plotOutput("fluct_plot")
           )
    ),
         

    conditionalPanel(
      condition = "input.show_options1 == true",
      fluidRow(
      column(3,
             h4("Plot options"),
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
      ), ## end condition
  
  br(),br(),


  ### WebGL visualization
  fluidRow(
    column(4,
           wellPanel(
             h4('Normal Modes Visualization'),
             checkboxInput('show_trj2', 'Show NM Trajectory', value=FALSE),
             
             selectInput('mode_choice', 'Choose Mode:', choices=c(7:50)),
             
             radioButtons('viewColor2', label='Structure color',
                          choices=list(
                            'Amalgam' = 'amalgam',
                            'Magnitude'='mag',
                            'By Frame (blue->gray->red)'='default'
                                   ),
                          selected='amalgam'),
             
             radioButtons('viewBGcolor2', label='Background color',
                          choices=list('Black'='black', 'White'='white'),
                          selected='white'),
             br(),
             downloadButton('nmtraj', label='Download PDB Trajectory')
             )
           ),
    
    column(8,
           conditionalPanel(
             condition='input.show_trj2 == true',
             webGLOutput('nmaWebGL')
             )
           )
    ),
         

  ##-C. DCCM Panel
  fluidRow(
    column(4,
           wellPanel(
             h4('DCCM'),
             checkboxInput('calc_dccm', 'Calculate DCCM', value=FALSE),
             ##checkboxInput('show_options2', 'More options', value=FALSE),
             checkboxInput('contourplot', 'Contourplot', value=TRUE),
             checkboxInput('sse', 'Show SSE', value=TRUE),
             checkboxInput('colorkey', 'Colorkey', value=TRUE),
             downloadButton('dccm2zip', "Download PyMOL Visualization script"),
             downloadButton('dccmplot2pdf', "Download Plot PDF")
             )
           ),
    column(8,
           plotOutput("dccm_plot")
           )
    ),     
   br(),br(),

         
  ##-C. Overlap analysis
  fluidRow(
    column(4,
           wellPanel(
             h4('Overlap analysis '),
             textInput("pdbid2", "Enter PDB ID", value=""),

             ##- Chain selection
             h5("Detected chain IDs:"),
             verbatimTextOutput("chains3"),
             
             checkboxInput("limit2", "Limit calculation to a subset of chains?"),
             conditionalPanel(
               condition = "input.limit2 == true",
               uiOutput("chains4")
               )
             
             )
           ),
    column(8,
           plotOutput("overlap_plot")
           )
    ),     
   br(),br(),

         
  ##-F. Domain analysis
  fluidRow(
    column(4,
           wellPanel(
             h4('Domain analysis '),
             checkboxInput('domains', 'Do domain analysis', value=FALSE),
             sliderInput("ndomains", "Number of domains:",
                         min = 2, max = 10, value = 3),
             sliderInput("nmodes", "Number of modes:",
                         min = 1, max = 5, value = 5),
             downloadButton('geostas2zip', "Download PDB Trajectory")
             )
           ),

    column(8, 
           conditionalPanel(
             condition = "input.domains == true",
             webGLOutput('geostasWebGL')
             )
           )
    ),
  br(),br(),

  

  hr(),
  ##-G. Log Report Panel
  h3("Summary Log"),
    column(6,
           h4("PDB Summary"),
           verbatimTextOutput("pdbSummary")
           ),
    column(6,
           h4("Calculation Summary Log"),
           verbatimTextOutput("modeSummary")
           )
         
)
