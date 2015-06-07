tabPanel("NMA", icon=icon("home"),
         tags$style(type="text/css", "body {padding-top: 80px;}"),
         
  fluidRow(
    column(4,
           wellPanel(
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

        ##-PDB input
        textInput("pdbid", label="Enter RCSB PDB code/ID:", value = "4Q21"),
        
        ##- Chain selection
        uiOutput("chain_input"),

        ## reset PDB input
        actionButton("reset_pdbid", "Reset PDB input", icon=icon("undo")),
        checkboxInput('show_pdb', 'View Input PDB', value=FALSE)
        )
   ),

    column(4,
      wellPanel(
        h4("2. NMA parameters"),
        tags$hr(),

        selectInput("forcefield", "Choose a forcefield:",
                    choices = c("calpha", "sdenm", "reach", "anm", "pfanm")),
        
        sliderInput("cutoff", "Cutoff value:",
                    min = 7, max = 50, value = 15),
        
        helpText("Note: Cutoff applies only to 'ANM' and 'pfANM'. Recommended values are 15 and 50 Å, respectively."),
 
        actionButton("reset_nma_input", "Reset NMA inputs", icon=icon("undo"))
        )
    )
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
             selectInput('mode_inds', 'Choose Mode indices:',
                         choices=c("all", 7:50), selected="all", multiple=TRUE),
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
                         min = 0.1, max = 2, value = 1, step=0.1),
             
             sliderInput("height1", "PDF height:",
                         min = 4, max = 12, value = 5, step=1)
             ),
      
      column(3,
             h4("Plot options"),
             sliderInput("pch1", "Point type:",
                         min = 1, max = 25, value = 1),
             sliderInput("cex1", "Point size:",
                         min = 0.1, max = 2, value = 1, step=0.1),
             sliderInput("col1", "Color:",
                         min = 1, max = 8, value = 1),
             sliderInput("width1", "PDF width:",
                         min = 4, max = 12, value = 7, step=1)
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
             sliderInput("mag", "Magnification factor:",
                         min = 1, max = 15, value = 5),
             
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
             downloadButton('trj2zip', label='Download PDB Trajectory'),
             downloadButton('nma2pymol', label='Download PyMOL vector field')
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
             actionButton("run_dccm", "Run DCCM", icon=icon("cog")),
             ##checkboxInput('calc_dccm', 'Calculate DCCM', value=FALSE),
             ##checkboxInput('show_options2', 'More options', value=FALSE),
             checkboxInput('contourplot', 'Contourplot', value=TRUE),
             checkboxInput('sse', 'Show SSE', value=TRUE),
             checkboxInput('colorkey', 'Colorkey', value=TRUE),

                
             sliderInput("height2", "PDF height:",
                         min = 4, max = 12, value = 5, step=1),
             sliderInput("width2", "PDF width:",
                         min = 4, max = 12, value = 7, step=1),
             
             ##downloadButton('dccm2py', "Download PyMOL Visualization script"),
             downloadButton('dccm2pymol', "Download PyMOL session"),
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
             helpText("Run BLAST and then select PDB IDs for overlap analysis"),
             actionButton("run_blast", "Run BLAST", icon=icon("cog")),
             actionButton("run_overlap", "Run Overlap", icon=icon("cog")),
             
             DT::dataTableOutput('blast_table'),
             downloadButton('overlapplot2pdf', "Download Plot PDF"),
             checkboxInput('show_options2', 'More options', value=FALSE)
             
             )
           ),
    column(8,
           plotOutput("overlap_plot")
           )
    ),
         
    conditionalPanel(
      condition = "input.show_options2 == true",
        fluidRow(
          column(3,
                 h4("Plot options"),
                 checkboxInput('show_legend3', 'Show legend', value=TRUE),
           
                 
                 sliderInput("height3", "PDF height:",
                             min = 4, max = 12, value = 5, step=1),
                 sliderInput("width3", "PDF width:",
                             min = 4, max = 12, value = 7, step=1)
                 ),
          column(3,
                 h4("Overlap values"),
                 checkboxInput('show_overlap', 'Plot overlap values', value=TRUE),
                 
                 selectInput("typ3", "Type:",
                             choices=c(
                               "hist" = "h",
                               "lines" = "l",
                               "points" = "p"
                               ), selected=c("p", "h"), multiple=TRUE),
                 sliderInput("cex3", "Point size:",
                             min = 0.1, max = 2, value = 1, step=0.1),
                 sliderInput("lty3", "Line type:",
                             min = 1, max = 6, value = 1),
                 sliderInput("lwd3", "Line width:",
                             min = 0.1, max = 2, value = 1, step=0.1)
                 ),
          column(3,
                 h4("Cumulative overlap values"),
                 checkboxInput('show_overlap_cum', 'Plot cumulative values', value=TRUE),
                 selectInput("typ4", "Type:",
                             choices=c(
                               "hist" = "h",
                               "lines" = "l",
                               "points" = "p"
                               ), selected=c("p", "l"), multiple=TRUE),
                 sliderInput("cex4", "Point size:",
                             min = 0.1, max = 2, value = 1, step=0.1),
                 sliderInput("lty4", "Line type:",
                             min = 1, max = 6, value = 1),
                 sliderInput("lwd4", "Line width:",
                             min = 0.1, max = 2, value = 1, step=0.1)
                 )
          )
      ),

         
   br(),br(),

         
  ##-F. Domain analysis
  fluidRow(
    column(4,
           wellPanel(
             h4('Domain analysis '),
             actionButton("run_geostas", "Run Geostas", icon=icon("cog")),
             ##checkboxInput('domains', 'Show domain analysis', value=FALSE),
             sliderInput("ndomains", "Number of domains:",
                         min = 2, max = 10, value = 3),
             sliderInput("nmodes", "Number of modes:",
                         min = 1, max = 5, value = 3),
             downloadButton('geostas2zip', "Download PDB Trajectory")
             )
           ),

    column(8, 
           #conditionalPanel(
           #  condition = "input.domains == true",
             webGLOutput('geostasWebGL')
           #  )
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
           h4("NMA Summary"),
           verbatimTextOutput("modeSummary")
           )
         
)
