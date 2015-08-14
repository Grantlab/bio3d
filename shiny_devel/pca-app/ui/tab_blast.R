
tabPanel("1. SEARCH", icon=icon("home"),
         tags$style(type="text/css", "body {padding-top: 80px;}"),
         div(
           h3("Structure Search", style="border: 1px solid #e3e3e3; border-radius: 4px; padding-top: 10px; padding-bottom: 10px; padding-left: 5px; margin-top: -10px; background-color: white;")
           ),


         actionButton3("about_blasttab", "About this tab", icon=icon("comment"), cl="btn btn-warn btn-input action-button", style = "position: fixed; top: 14px; right: 16px; z-index: 2000;"),

         bsModal("modal_blast", "About the Bio3D PCA/eNMA WebApp", "about_blasttab", size = "large",
                 content=tags$div(
                   ##h3("Bio3D PCA/eNMA WebApp"),
                   p(HTML("This <a href=\"http://thegrantlab.org/bio3d/index.php\">Bio3D</a> WebApp provides a rapid and rigorous tool for comparative structure analysis of protein families. Methods include inter-conformer characterization with <a href=\"http://thegrantlab.org/bio3d/tutorials/principal-component-analysis\">principal component analysis</a> (PCA) and <a href=\"http://thegrantlab.org/bio3d/tutorials/ensemble-nma-part-1\">ensemble normal mode analysis</a> (eNMA).")),
                   p(HTML("Start by entering a PDB code of interest then proceed by navigating through the above tabs or following the <b>NEXT</b> buttons.")),
                   img(src="./images/geostas_250x182.png", width=250, style="display: block; margin-left: auto; margin-right: auto;")
                   )
                 ),

         fluidRow(
           column(4,
                  wellPanel(
                    bsPopover("popblast1",
                              "Select your input data type",
                              "For both <b>single structure</b> and <b>single sequence</b> options a search will be performed to find related PDB structures upon which subsequent analysis will be based. The results of this search will be presented below along with options to restrict subsequent analysis to certain chains. </br></br>With <b>multiple structure</b> input, analysis will be confined to the specified structures and only their annotation will be displayed below. </br></br>To continue analysis proceed by navigating through the <b>NEXT</b> buttons.",
                              placement = "right", trigger = "hover",
                              options = list(container = "body")),

                    tags$div(id = "popblast1", icon("question-circle"),
                             style = "position: absolute; right: 25px; top: 5px;"
                             ),

                    h4("A)  Input Structure(s) or Sequence"),
                    helpText("Please enter either a single PDB code of interest, a single protein sequence, or multiple related PDB codes (see the ",
                             a(href="http://thegrantlab.org/bio3d", target="_blank", "Help"), " page for more details)."),

                    radioButtons(inputId = "input_type", label = "",
                                 choices = c("Enter a single PDB structure code" = "pdb",
                                   "Paste a single protein sequence" = "sequence",
                                   "Enter multiple PDB structure codes" = "multipdb"),
                                 selected = NULL, inline=FALSE),
                    br(),

                    conditionalPanel(
                      condition = "input.input_type == 'multipdb'",
                      tags$textarea(id="pdb_codes", rows=4, cols=40, "1TND, 1KJY_A"),
                      ##textInput("pdb_codes", label = "hei", value = "1TND, 1KJY_A"),
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
                       '$("#page1_hits").click(function(){',
                       '$("html, body").animate({scrollTop:$("#blast_row_plot").position().top - (0.1 * $(window).height())}, "smooth");',
                       '$("#blast_plot").parent().siblings().find(".well").addClass("show-border");',
                       'window.setTimeout(function(){',
                       '$("#blast_plot").parent().siblings().find(".well").removeClass("show-border");',
                       '}, 2500);',
                       '});'
                                    )),
                    #actionButton("reset_pdbid", "Reset", icon=icon("undo")),
                    actionButton3("reset_pdbid", HTML("<b>Reset</b>"), icon=icon("undo"), cl="btn btn-link btn-input action-button")

                    )
                  ),

           column(4,
                  wellPanel(
                    style="background: #FFFFFF;",

                    bsPopover("popblast2",
                              "Chain selection, annotation and visualization options",
                              "This panel presents summary information about your input specified in panel <b>A</b> to the left. <br><br> Options here allow you to limit further analysis to a certain chain as well as <b>View</b> the <b>3D structure</b> of each chain and inspect their <b>PFAM chain annotations</b>. <br><br>To continue analysis proceed by navigating through the <b>NEXT</b> buttons.",
                              placement = "right", trigger = "hover",
                              options = list(container = "body")),

                    tags$div(id = "popblast2", icon("question-circle"),
                             style = "position: absolute; right: 25px; top: 5px;"
                             ),

                    conditionalPanel(
                      condition = "input.input_type == 'multipdb'",
                      h4("Multiple PDB Summary"),
                      helpText("Annotation table of requested PDBs.")
                      #tags$textarea(id="pdb_codes", rows=4, cols=40, ""),
                      #helpText("Seperate PDB ids with ','")
                                            ##- pfam_table (To be beautified!!)
                      ,tags$label("PFAM Annotation:")
                      ,dataTableOutput("pfam_table_multi")

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
                      h4("Structure Summary"),


                      ##- PDB summary
                      tags$label("PDB Summary:"),
                      verbatimTextOutput("input_pdb_summary"),

                      br(),

                      ##- pfam_table (To be beautified!!)
                      ##tags$label("PFAM chain annotation:"),
                      #dataTableOutput("pfam_table"),
                      DT::dataTableOutput('pfam_table'),

                      br(),br(),

                      ##- Chain selection
                      uiOutput("pdb_chains")

                      #,br()

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


                      #popRadioButtons("logviewer2", "View:",
                      #              c("3D structure" = "pdb",
                      #                "PDB processing log" = "pdblog",
                      #                "App Info" = "bio3d"),
                      #                selected="pdb",   ##<--- Change to "pdb"/"bio3d" ???
                      #               inline=TRUE,
                      #               placement = "right",
                      #               data_toggle = "pop_summary_input",
                      #               title = "Chain selection, annotation and visualization options")

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

                      selectInput("view_inpdb_as", "Display options:",
                                  c("Overview" = "overview",
                                    "C-alpha Trace" = "calpha",
                                    "All atoms" = "allatoms"),
                                  selected = "overview",
                                  multiple = FALSE),

                      selectInput("view_inpdb_col", "Color options:",
                                   c("Secondary structure" = "sse",
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
                      p(HTML("This <a href=\"http://thegrantlab.org/bio3d/index.php\">Bio3D</a> WebApp provides a rapid and rigorous tool for comparative structure analysis of protein families. Methods include inter-conformer characterization with <a href=\"http://thegrantlab.org/bio3d/tutorials/principal-component-analysis\">principal component analysis</a> (PCA) and <a href=\"http://thegrantlab.org/bio3d/tutorials/ensemble-nma-part-1\">ensemble normal mode analysis</a> (eNMA).")),
                      p(HTML("Start by entering a PDB code of interest then proceed by navigating through the above tabs or following the <b>NEXT</b> buttons.")),##<font color=\"red\">NEXT</font> buttons.")),
                      img(src="./images/geostas_250x182.png",
                          width=250, style="display: block; margin-left: auto; margin-right: auto;")
                      )
###
                    ,conditionalPanel(
                                      condition = "input.input_type == 'pdb'",
                                      radioButtons("logviewer", "View options:",
                                    c("3D structure" = "pdb",
                                      "PDB processing log" = "pdblog",
                                      "App Info" = "bio3d"),
                                      selected="pdb",   ##<--- Change to "pdb"/"bio3d" ???
                                     inline=TRUE)
                    )


###

                    )
                  )
           ),

         #########################
         ##-- Results Section --##
         #########################
         fluidRow(
           id = 'blast_row_plot',
           #conditionalPanel(
           #  condition = "input.input_type != 'multipdb'",
           #  ##h2("Blast results")
           #  ),
           hr(),

           ##-- B. blast panel
           conditionalPanel(
             condition = "input.input_type != 'multipdb'",

             column(4,
                    wellPanel(

                      conditionalPanel(
                        condition = "input.input_type != 'multipdb'",

                        bsPopover("popblast3",
                                  "Hit selection",
                                  "The plot on the right shows a representation of the hits in the PDB identified in the search. Each point represent a hit - or a PDB structure. All points above the red cutoff line can be chosen for further analysis. Note that only those marked with a blue circle is currently chosen. ",
                                  placement = "right", trigger = "hover",
                                  options = list(container = "body")),

                        tags$div(id = "popblast3", icon("question-circle"),
                                 style = "position: absolute; right: 25px; top: 5px;"
                                 ),


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
                       '$("html, body").animate({scrollTop:$("#blast_row_table").position().top }, "smooth");',
                       #'$("#blast_table").parent()[0].scrollIntoView({block: "end", behavior: "smooth"});',
                       #'window.scrollBy(0,-50);',
                       '$("#blast_table").parent().addClass("show-border");',
                       'window.setTimeout(function(){',
                       '$("#blast_table").parent().removeClass("show-border");',
                       '}, 2500);',
                       '});'
                                    )),
##                        actionButton("reset_cutoff", "Reset cutoff", icon=icon("undo"))
                        actionButton3("reset_cutoff", HTML("<b>Reset</b>"), icon=icon("undo"), cl="btn btn-link btn-input action-button")


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
           id = 'blast_row_table',
           hr(),
           column(12,
                ## PopOver probably not needed...?...
                #popoverQuestion(id="popQues4", content="This table provides annotation of identified structures and allows for selection (and de-selection) via <b>clicking to highlight</b> individual rows. This allows finer grained selection than the sliders in panel B above.",
                #                data_toggle = "blast_table"),

                  wellPanel(
                    h4("C) Optional filtering of related structures for further analysis"),
                    helpText(HTML("Optionally select (or de-select) structures by <b>clicking to highlight</b> their entries in the below table. This allows for finer grained selection than the sliders in panel B above.")),

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
