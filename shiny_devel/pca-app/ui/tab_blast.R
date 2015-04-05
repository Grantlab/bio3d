tabPanel("1. BLAST",
         
         fluidRow(
           column(6,
                  wellPanel(
                    h4("1. Sequence Input"),
                    ##tags$hr(),

                    radioButtons("input_type", "",
                                 c("Enter PDB code" = "pdb",
                                   "Paste a sequence" = "sequence"),
                                 inline=TRUE),
                    
                    conditionalPanel(
                      condition = "input.input_type == 'sequence'",
                      tags$textarea(id="sequence", rows=4, cols=40, ""),
                      actionButton("action_input", "Go")
                      ),
                    
                    conditionalPanel(
                      condition = "input.input_type == 'pdb'",
                      
                      ##-PDB input 
                      uiOutput('resetable_pdb_input'),
                    
                      ##- Chain selection
                      uiOutput("chains2"),
                      actionButton("reset_pdb_input", "Reset PDB input")
                      )
                    )
                  ),
           
           column(6,
                  wellPanel(
                    h4("2. BLAST"),
                    
                    tags$hr(),
                    uiOutput("cutoff_slider"),
                    uiOutput("hits_slider")
                    )
                  )
           
           #column(4,
           #       wellPanel(
           #         h4("3. PDB select"),
           #         tags$hr(),
           #         uiOutput("hits_slider"),
           #         uiOutput("pdbids_checkboxgroup")
           #         )
           #       )
           ),
         
         #########################
         ##-- Results Section --##
         #########################
         h2("Blast results"),
         hr(),
         
         ##-A. blast panel
         fluidRow(
           column(12, 
                  plotOutput("blast_plot")
                  )
           ),

         fluidRow(
           column(12,
                  wellPanel(
                    dataTableOutput("blast_table")
                    )
                  )
           ),
         

         fluidRow(
           column(12,
                  wellPanel(
                    h4("Uncheck PDB ID to exclude from analysis"),
                    uiOutput("pdbids_checkboxgroup")
                    )
                  )
           )
         
         
         )
