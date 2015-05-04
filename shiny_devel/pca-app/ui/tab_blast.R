tabPanel("1. BLAST", icon=icon("home"),
         tags$style(type="text/css", "body {padding-top: 80px;}"),

         ## experimental load and save results (from previous runs)
         #verbatimTextOutput("saveText"),
         #actionButton("saveButton", "Save", icon=icon("undo")),

         #uiOutput('datafileslist'),
         #actionButton("loadButton", "Load", icon=icon("undo")),
         #verbatimTextOutput("loadText"),


         fluidRow(
           column(6,
                  wellPanel(
                    h4("A) Input query structure or sequence"),
                    hr(),

                    radioButtons("input_type", "",
                                 c("Enter PDB code" = "pdb",
                                   "Paste a sequence" = "sequence",
                                   "Enter mutliple PDB codes" = "multipdb"),
                                 inline=TRUE),

                    conditionalPanel(
                      condition = "input.input_type == 'multipdb'",
                      tags$textarea(id="pdb_codes", rows=4, cols=40, "")
                      ),

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
                      verbatimTextOutput("input_pdb_summary"),

                      actionButton("reset_pdb_input", "Reset PDB input", icon=icon("undo"))
                      )
                    )
                  ),


           column(6,
                  wellPanel(

                    conditionalPanel(
                      condition = "input.input_type != 'multipdb'",

                      h4("B) Hit selection for further analysis"),

                      tags$hr(),
                      uiOutput("resetable_cutoff_slider"),
                      uiOutput("hits_slider"),
                      actionButton("reset_cutoff", "Reset cutoff", icon=icon("undo"))
                      ),

                    conditionalPanel(
                      condition = "input.input_type == 'multipdb'",
                      h4("BLAST options not available for multiple PDB input")
                      )
                    )
                  )
           ),


         #########################
         ##-- Results Section --##
         #########################
         conditionalPanel(
           condition = "input.input_type != 'multipdb'",
           h2("Blast results")
           ),
         hr(),

         ##-A. blast panel
         conditionalPanel(
           condition = "input.input_type != 'multipdb'",

           fluidRow(
             column(12,
                    uiOutput("blast_plot"),

                    tags$script(HTML(
                      'var css = document.createElement("style");
                      css.type = "text/css";
                      css.innerHTML = ".nv-x .nv-axislabel { font-size: 20px; }";
                      document.body.appendChild(css);
                      css = document.createElement("style");
                      css.type = "text/css";
                      css.innerHTML = ".nv-y .nv-axislabel { font-size: 20px; }";
                      document.body.appendChild(css);
                      document.getElementById("blast_plot2").focus();'
                      ))
                    )
             )
           ),


         fluidRow(
           column(12,
                  wellPanel(
                    dataTableOutput("blast_table")
                    )
                  )
           )

         )
