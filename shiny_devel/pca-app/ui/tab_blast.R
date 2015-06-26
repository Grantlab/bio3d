
tabPanel("1. SEARCH", icon=icon("home"),
         tags$style(type="text/css", "body {padding-top: 80px;}"),
         
         fluidRow(
           column(4,
                  wellPanel(
                    popoverQuestion(id="popQues1", content="test", trigger="focus",
                                    data_toggle = "pop_blast_input"),
                    
                    h4("A)  Input Structure(s) or Sequence"),
                    
                    helpText("Please enter either a single PDB code of interest, multiple related PDB codes or a single protein sequence (see the help page for more details)."),
                    
                    popRadioButtons(inputId = "input_type", label = "",
                                 choices = c("Enter a single PDB code" = "pdb",
                                   "Paste a single sequence" = "sequence",
                                   "Enter mutliple PDB codes" = "multipdb"),
                                 selected = NULL,inline=FALSE,
                                 placement = "right", data_toggle = "pop_blast_input",
                                 title = "Structure or sequence input"),

                    conditionalPanel(
                      condition = "input.input_type == 'multipdb'",
                      tags$textarea(id="pdb_codes", rows=4, cols=40, "1TND, 1KJY"),
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
                      popTextInput("pdbid", label="Enter a 4 character RCSB PDB code/ID:", value = "2LUM", data_toggle = "pop_blast_input", title = "RCSB PDB ID", content = "Please enter a four letter PDB code. We parse the first four characters only.")
                      ),

                    actionButton3("page1_hits", "Next (Hit selection)", icon=icon("arrow-down"), cl="btn btn-primary btn-input action-button"),
                   tags$script(HTML(
                       '$(".btn-input").click(function(){',
                       'document.getElementById("blast_plot").scrollIntoView({block: "start", behavior: "smooth"});',
                       'window.scrollBy(0,-100);',
                       '$("#blast_plot").parent().siblings().find(".well").addClass("show-border");',
                       'window.setTimeout(function(){',
                       '$("#blast_plot").parent().siblings().find(".well").removeClass("show-border");',
                       '}, 2500);',
                       '});'
                                    )),
                    actionButton("reset_pdbid", "Reset PDB", icon=icon("undo"))
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
                                            ##- pfam_table (To be beautified!!)
                      ,tags$label("PFAM Annotation:")
                      #dataTableOutput("pfam_table")
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

                      ##- pfam_table (To be beautified!!)
                      tags$label("PFAM chain annotation:"),
                      dataTableOutput("pfam_table"),

                      ##submitButton("Update"),
                      radioButtons("logviewer", "PDB Log:",
                                   c("About Bio3D" = "bio3d",
                                     "View structure" = "pdb",
                                     "PDB Read Log" = "pdblog"), inline=TRUE)
                      )

                    )
                  ),

           column(4,
                  wellPanel(
                    style="background: #FFFFFF;",

                    conditionalPanel(
                      condition = "input.logviewer == 'pdb' && input.input_type == 'pdb'",
                      h4("Input Structure Visualization"),
                      webGLOutput('pdbWebGL'),

                      selectInput("view_inpdb_as", "View mode",
                                  c("C-alpha Trace" = "calpha",
                                    "Overview" = "overview"),
                                  selected = "calpha",
                                  multiple = FALSE), 
                      
                      selectInput("view_inpdb_col", "Color C-alpha trace by",
                                   c("Secondary structure elemets" = "sse",
                                     "Residue Index" = "index"),
                                  multiple=FALSE)

                      ),

                    conditionalPanel(
                      condition = "input.logviewer == 'pdblog' && input.input_type == 'pdb'",
                      h4("Input PDB Read Log"),
                      verbatimTextOutput('pdb_log')
                      ),

                    conditionalPanel(
                      condition = "input.logviewer == 'bio3d' || input.input_type != 'pdb'",
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
                        actionButton3("page1_table", "Next (Further selection)", icon=icon("arrow-down"), cl="btn btn-primary btn-hits action-button"),
                        tags$script(HTML(
                       '$(".btn-hits").click(function(){',
                       '$("#blast_table").parent()[0].scrollIntoView(false);',
                       #'window.scrollBy(0,-50);',
                       '$("#blast_table").parent().addClass("show-border");',
                       'window.setTimeout(function(){',
                       '$("#blast_table").parent().removeClass("show-border");',
                       '}, 2500);',
                       '});'
                                    )),
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

         fluidRow(
           hr(),
           column(12,
                  wellPanel(
                    h4("C) Optional filtering of related structures for further analysis"),
                    helpText("Optionaly select or de-select structures by clicking their entries in the below table."),

                    DT::dataTableOutput('blast_table'),
                    hr(),
                    actionButton3("page2", "Next (Alignment)", icon=icon("arrow-right"), cl="btn btn-primary btn-next-blast action-button"),
                    tags$script(HTML(
                       '$(".btn-next-blast").click(function(){',
                       'tabs = $(".nav.navbar-nav li");',
                       'tabs[1].childNodes[1].click()',
                       '});'
                    ))
                    )

             )
           )

         )
