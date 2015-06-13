tabPanel("1. SEARCH", icon=icon("home"),
         tags$style(type="text/css", "body {padding-top: 80px;}"),
         ##shinyjs::useShinyjs(),

         
         
         fluidRow(
           column(4,
                  wellPanel(

                    h4("A)  Input Structure(s) or Sequence"),
                    helpText("Please enter either a single PDB code of interest, multiple related PDB codes or a single protein sequence (see the help page for more details)."),
                    
                    radioButtons("input_type", "",
                                 c("Enter a single PDB code" = "pdb",
                                   "Paste a single sequence" = "sequence",
                                   "Enter mutliple PDB codes" = "multipdb"),
                                 inline=FALSE),
                    
                    conditionalPanel(
                      condition = "input.input_type == 'multipdb'",
                      tags$textarea(id="pdb_codes", rows=4, cols=40, ""),
                      helpText("Separate PDB ids (4 character codes) with a comma ','")
                      ),
                    
                    conditionalPanel(
                      condition = "input.input_type == 'sequence'",
                      tags$textarea(id="sequence", rows=4, cols=40,
                                    "MQYKLVINGKTLKGETTTKAVDAETAEKAFKQYANDNGVDGVWTYDDATKTFTVTE")
                      ),
                    
                    conditionalPanel(
                      condition = "input.input_type == 'pdb'",
                      
                      ##-PDB input
                      textInput("pdbid", label="Enter a 4 character RCSB PDB code/ID:", value = "2LUM")
                      ),

                    actionButton("reset_pdbid", "Reset PDB input", icon=icon("undo"))
                    )
                  ),
           
           
           column(4,
                  wellPanel(
                    style="background: #FFFFFF;",

                    conditionalPanel(
                      condition = "input.input_type == 'multipdb'",
                      h4("Multiple PDB Summary"),
                      helpText("Annotation table of requested PDBs.")
                      #tags$textarea(id="pdb_codes", rows=4, cols=40, ""),
                      #helpText("Seperate PDB ids with ','")
                      ),
                    
                    conditionalPanel(
                      condition = "input.input_type == 'sequence'",
                      h4("Protein Sequence Summary"),
                      helpText("Report on whether sequence was acceptable (protein, no strange characters, sufficient length, top PDB hit id and link for visualization option?).")
                      #tags$textarea(id="sequence", rows=4, cols=40, "")
                      ),
                    
                    conditionalPanel(
                      condition = "input.input_type == 'pdb'",
                      h4("Structure Summary and Visualization"),

                      ##-PDB input
                      #textInput("pdbid", label="Enter RCSB PDB code/ID:", value = "2LUM"),
                                            
                      ##- PDB summary
                      tags$label("PDB Summary:"),
                      verbatimTextOutput("input_pdb_summary"),
                      
                      ##- Chain selection
                      uiOutput("pdb_chains"),
                      
                      ##submitButton("Update"),
                      checkboxInput('show_pdb', 'View Input PDB', value=FALSE),

                      checkboxInput('show_pdblog', 'View PDB Read Log', value=FALSE)
                      )
                    
                    )
                  ),
           
           
           column(4,
                  wellPanel(
                    style="background: #FFFFFF;",

                    conditionalPanel(
                      condition = "input.show_pdb == true",
                      h4("Input Structure Visualization"),
                      webGLOutput('pdbWebGL'),
                      radioButtons("view_inpdb_as", "View as",
                                   c("Overview" = "overview",
                                     "Calpha trace" = "calpha"),
                                   inline=TRUE),
                      radioButtons("view_inpdb_col", "Color by",
                                   c("SSE" = "sse",
                                     "Index" = "index"),
                                   inline=TRUE)
                      
                      ),

                    conditionalPanel(
                      condition = "input.show_pdb == false",
                      h4("Bio3D PCA/eNMA WebApp"),
                      p("This Bio3D WebApp provides a rapid and rigorous tool for comparative structure analysis of protein families. Methods include inter-conformer characterization with principal component analysis (PCA) and ensemble normal mode analysis (eNMA)."),
                      p("Start by entering a PDB code of interest then proceed by navigating through the above tabs."),
                      img(src="geostas_250x182.png",
                          width=250, style="display: block; margin-left: auto; margin-right: auto;")
                      )


                    )
                  )
           ),


         
         #########################
         ##-- Results Section --##
         #########################
         fluidRow(
           #conditionalPanel(
           #  condition = "input.input_type != 'multipdb'",
           #  ##h2("Blast results")
           #  ),
           hr(),
           
           ##-A. blast panel
           conditionalPanel(
             condition = "input.input_type != 'multipdb'",
             
             column(4,
                    wellPanel(
                      
                      conditionalPanel(
                        condition = "input.input_type != 'multipdb'",
                        
                        h4("B) Hit selection for further analysis"),
                        helpText("Optional refinement (via exclusion and selection) of related structures. The final selected set of structures will from the ensemble for analysis in subsequent steps."),

                        ##tags$hr(),
                        uiOutput("cutoff_slider"),
                        
                        uiOutput("hits_slider"),
                        actionButton("reset_cutoff", "Reset cutoff", icon=icon("undo"))
                        ),
                      
                      conditionalPanel(
                        condition = "input.input_type == 'multipdb'",
                        h4("BLAST options not available for multiple PDB input")
                        )
                      )
                    ),
             
             column(8,
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
                      document.getElementById("blast_plot").focus();
                      '
                      ))
                    )
             )
           ),

         hr(),
         fluidRow(
           column(12,
                  wellPanel(
                    h4("Manual hit filtering"),
                    uiOutput("include_hits"),
                    radioButtons("filter_sorting", "Order by",
                                 c("BLAST score" = "blast",
                                   "PDB ID" = "pdbid"),
                                 inline=TRUE)
                    ))
           ),
         
         hr(),
         fluidRow(
           column(12,
                  wellPanel(
                    h4("C) Optional filtering of related structures for further analysis"),
                    helpText("Select structures by clicking their entries in the below table."),
                     
                    DT::dataTableOutput('blast_table')
                    )
             
             )
           )

         )
