tabPanel("1. BLAST",
         
         fluidRow(
           column(4,
                  wellPanel(
                    h4("1. PDB Input Selection"),
                    tags$hr(),
                    
                    ##-PDB input (moved to server.R)
                    uiOutput('resetable_pdb_input'),
                    
                    ##- Chain selection
                    h5("Detected chain IDs:"),
                    verbatimTextOutput("chains1"),
                    
                    checkboxInput("limit", "Limit calculation to a subset of chains?"),
                    helpText("Note: Use this option to exclude particular chains form further consideration."),
                    
                    conditionalPanel(
                      condition = "input.limit == true",
                      uiOutput("chains2"),
                      helpText("Note: Only selected chains will be analyzed.")
                      ),
                    
                    actionButton("reset_pdb_input", "Reset PDB input")
                    )
                  ),
           
           column(4,
                  wellPanel(
                    h4("2. BLAST"),
                    radioButtons("blast_prg", "",
                                 c("HMMER" = "hmmer",
                                   "BLAST" = "blast")),
                    helpText("currently only hmmer implemented"),
                    
                    tags$hr(),
                    sliderInput("cutoff", "Cutoff:",
                                min = 0, max = 700, value = 250),
                    
                    helpText("TODO: default value should be obtained from plot.blast / plot.hmmer")
                    )
                  ),
           
           column(4,
                  wellPanel(
                    h4("3. PDB select"),
                    tags$hr(),

                    uiOutput("hits_slider"),
                    uiOutput("pdbids_checkboxgroup")
                    

                    )
                  )
           ),
         
         #########################
         ##-- Results Section --##
         #########################
         h2("Blast plot"),
         hr(),
         
         ##-A. blast panel
         fluidRow(
           plotOutput("blast_plot")
           )
         
         )
