
tabPanel("1. SEARCH", icon=icon("home"),
         tags$style(type="text/css", "body {padding-top: 80px;}"),
         
         fluidRow(
           column(4,
                  wellPanel(
                    popoverQuestion(id="popQues1", 
                                    content="For both <b>single structure</b> and <b>single sequence</b> options a search will be performed to find related PDB structures upon which subsequent analysis will be based. The results of this search will be presented below along with options to restrict subsequent analysis to certain chains. </br></br>With <b>multiple structure</b> input, analysis will be confined to the specified structures and only their annotation will be displayed below. </br></br>To continue analysis proceed by navigating through the <b>NEXT</b> buttons. </br>Please refer to the main <a href='http://thegrantlab.org'>Help</a> page for further details.", 
                                    trigger="focus",
                                    data_toggle = "pop_blast_input"),
                    
                    h4("A)  Input Structure(s) or Sequence"),

                    helpText("Please enter either a single PDB code of interest, a single protein sequence, or multiple related PDB codes (see the ",
                    a(href="http://thegrantlab.org", target="_blank", "Help"), " page for more details)."),

                    popRadioButtons(inputId = "input_type", label = "",
                                 choices = c("Enter a single PDB structure code" = "pdb",
                                   "Paste a single protein sequence" = "sequence",
                                   "Enter multiple PDB structure codes" = "multipdb"),
                                 selected = NULL, inline=FALSE,
                                 placement = "right", data_toggle = "pop_blast_input",
                                 title = "Select your input data type"),##-- SEE popQues1 above for content... 

                    br(),

                    conditionalPanel(
                      condition = "input.input_type == 'multipdb'",
                      tags$textarea(id="pdb_codes", rows=4, cols=40, "1TND, 1KJY_A"),
                      helpText("Enter multiple comma ',' separated PDB IDs (4 character RCSB PDB codes with optional underscore chain, e.g. '1KJY_A')")
                      ),

                    conditionalPanel(
                      condition = "input.input_type == 'sequence'",
                      tags$textarea(id="sequence", rows=4, cols=40,
                                    "MQYKLVINGKTLKGETTTKAVDAETAEKAFKQYANDNGVDGVWTYDDATKTFTVTE"),
                      helpText("Paste a protein sequence with no identifiers or FASTA headers.")
          
                      ),

                    conditionalPanel(
                      condition = "input.input_type == 'pdb'",

                      ##- PDB input
                      textInput("pdbid", label="Enter a 4 character RCSB PDB code/ID:", value = "2LUM")

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

                    popoverQuestion(id="popQues2", 
                                    content="This panel presents summary information about your input specified in panel <b>A</b> to the left. </br></br> Options here allow you to limit further analysis to a certain chain as well as <b>View</b> the <b>3D structure</b> of each chain and inspect their <b>PFAM chain annotations</b>. </br></br>To continue analysis proceed by navigating through the <b>NEXT</b> buttons.", trigger="focus",
                                    data_toggle = "pop_summary_input"),

                    conditionalPanel(
                      condition = "input.input_type == 'multipdb'",
                      h4("Multiple PDB Summary"),
                      helpText("Annotation table of requested PDBs.")
                      #tags$textarea(id="pdb_codes", rows=4, cols=40, ""),
                      #helpText("Seperate PDB ids with ','")
                                            ##- pfam_table (To be beautified!!)
                      ,tags$label("PFAM Annotation:")
                      #dataTableOutput("pfam_table")

                      ##-- ToDo:
                      ## Need a table here PFAM annotations

                      ),

                    conditionalPanel(
                      condition = "input.input_type == 'sequence'",
                      h4("Protein Sequence Summary"),
                      helpText("Report on whether sequence was acceptable (protein, no strange characters, sufficient length, top PDB hit id and link for visualization option?).")
                      #tags$textarea(id="sequence", rows=4, cols=40, "")

                      ##-- ToDo:
                      ## Need a table here of sequence properties and PFAM annotation?
                      ## See function mk.pfam.tbl() in server/util.R for a starting point
                      ),

                    conditionalPanel(
                      condition = "input.input_type == 'pdb'",
                      h4("Structure Summary and Visualization"),


                      ##- PDB summary
                      tags$label("PDB Summary:"),
                      verbatimTextOutput("input_pdb_summary"),

                      ##- Chain selection
                      uiOutput("pdb_chains"),

                      ##- pfam_table (To be beautified!!)
                      tags$label("PFAM chain annotation:"),
                      dataTableOutput("pfam_table"),

                      br(),

                      # radioButtons("logviewer", "PDB Log:",
                      #              c("App Info" = "bio3d",
                      #                "3D structure" = "pdb",
                      #                "PDB processing log" = "pdblog"), inline=TRUE),
                      #
                      # bsPopover("logviewer",
                      #            title = "Chain selection, annotation and visualization options",
                      #            content = "Report on whether sequence was acceptable (protein, no strange characters, sufficient length, top PDB hit id and link for visualization option?).",
                      #            trigger = "hover",
                      #            placement = "right")


                      popRadioButtons("logviewer", "View:",
                                    c("App Info" = "bio3d",
                                      "3D structure" = "pdb",
                                      "PDB processing log" = "pdblog"), 
                                     inline=TRUE,
                                     placement = "right", 
                                     data_toggle = "pop_summary_input",
                                     title = "Chain selection, annotation and visualization options")
                      
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
                                   c("Secondary structure elements" = "sse",
                                     "Residue Index" = "index"),
                                  multiple=FALSE)

                      ,bsTooltip("view_inpdb_as", 
                                  title="Change the structure view and coloring options. </br></br>Click and drag on the display to rotate, middle mouse button to zoom.",
                                  placement = "left", options = list(container = "body"))

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
                      p("Start by entering a PDB code of interest then proceed by navigating through the above tabs or following the <b>NEXT</b> buttons."),
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
                        helpText("Optional filtering and refinement (via similarity threshold specification) of related structures for further analysis."),

##                      ##-- Probably not needed...?...
##                      popoverQuestion(id="popQues3", content="This panel and the subsequent table allow for refinement of the structure set used for further analysis. </br></br> Move the sliders to increase or decrease the selected structures. Selected structures will be reflected in both the plot to the right and the table below.  </br></br>To continue analysis proceed by navigating through the <b>NEXT</b> buttons.", trigger="focus",
##                                    data_toggle = "pop_blast_plot"),

                        ##- Sliders
                        uiOutput("cutoff_slider"),
                        uiOutput("hits_slider"),

                        ##- Slider ToolTips
                        bsTooltip("cutoff_slider", 
                                  title="Sets the bitscore similarity threshold. Structures above this similarity threshold may be subject to further analysis and are indicated in red in the plot to the right and with a red icon in the table below.",
                                  placement = "right", options = list(container = "body")),

                        bsTooltip("hits_slider", 
                                  title="Sets the maximum number of structures to be used for further analysis. </br></br>The maximum permitted value is dependent on the bitscore cutoff chosen above. </br></br>Note that selected structures are indicated in blue in the plot to the right and also highlighted in the table below.",
                                  placement = "right", options = list(container = "body")),


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
                ## PopOver probably not needed...?...
                #popoverQuestion(id="popQues4", content="This table provides annotation of identified structures and allows for selection (and de-selection) via <b>clicking to highlight</b> individual rows. This allows finer grained selection than the sliders in panel B above.", 
                #                data_toggle = "blast_table"),

                  wellPanel(
                    h4("C) Optional filtering of related structures for further analysis"),
                    helpText("Optionally select (or de-select) structures by <b>clicking to highlight</b> their entries in the below table. This allows for finer grained selection than the sliders in panel B above."),

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
