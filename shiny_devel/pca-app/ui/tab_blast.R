tabPanel("1. BLAST", icon=icon("home"), 
         tags$style(type="text/css", "body {padding-top: 80px;}"),
         
         fluidRow(
           column(6,
                  wellPanel(
                    h4("A) Input query structure or sequence"),
                    hr(),

                    radioButtons("input_type", "",
                                 c("Enter PDB code" = "pdb",
                                   "Paste a sequence" = "sequence"),
                                 inline=TRUE),
                    
                    conditionalPanel(
                      condition = "input.input_type == 'sequence'",
                      tags$textarea(id="sequence", rows=4, cols=40, "")
                      #actionButton("action_input", "Go")
                      ),
                    
                    conditionalPanel(
                      condition = "input.input_type == 'pdb'",
                      
                      ##-PDB input 
                      uiOutput('resetable_pdb_input'),
                    
                      ##- Chain selection
                      uiOutput("pdb_chains"),
                      actionButton("reset_pdb_input", "Reset PDB input", icon=icon("undo"))
                      )
                    )
                  ),
           
           column(6,
                  wellPanel(
                    h4("B) Hit selection for further analysis"),
                    
                    tags$hr(),
                    uiOutput("resetable_cutoff_slider"),
                    uiOutput("hits_slider"),
                    actionButton("reset_cutoff", "Reset cutoff", icon=icon("undo"))
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
#                   plotOutput("blast_plot")
                  showOutput("blast_plot","nvd3"),
                    tags$script(HTML(
                      'var css = document.createElement("style");
                      css.type = "text/css";
                      css.innerHTML = ".nv-x .nv-axislabel { font-size: 20px; }";
                      document.body.appendChild(css);
                      css = document.createElement("style");
                      css.type = "text/css";
                      css.innerHTML = ".nv-y .nv-axislabel { font-size: 20px; }";
                      document.body.appendChild(css);')
                    )
                  )
           ),

#         fluidRow(
#           column(12, 
#                  showOutput("blast_plot2", "dimple")
#                  )
#           ),

         fluidRow(
           column(12,
                  wellPanel(
                    dataTableOutput("blast_table")
                    )
                  )
           )        
         
         )
