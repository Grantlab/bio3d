tabPanel("NMA", icon=icon("home"),
         tags$style(type="text/css", "body {padding-top: 80px;}"),
         
         tags$button(id = "about_app", type = "button",
                class = "btn btn-warn btn-input action-button",
                style = "position: fixed; top: 14px; right: 16px; z-index: 2000;",
                list(icon = icon("comment"), label = "About this app")),
         
         bsModal("modal_blast", "About the Bio3D NMA WebApp", "about_app", size = "large",
                 content=tags$div(
                   
                   p(HTML("This <a href=\"http://thegrantlab.org/bio3d/index.php\">Bio3D</a> WebApp provides a rapid and rigorous tool for comparative structure analysis of protein families. Methods include inter-conformer characterization with <a href=\"http://thegrantlab.org/bio3d/tutorials/principal-component-analysis\">principal component analysis</a> (PCA) and <a href=\"http://thegrantlab.org/bio3d/tutorials/ensemble-nma-part-1\">ensemble normal mode analysis</a> (eNMA).")),
                   p(HTML("Start by entering a PDB code of interest then proceed by navigating through the above tabs or following the <b>NEXT</b> buttons.")),
                   img(src="./images/geostas_250x182.png", width=250, style="display: block; margin-left: auto; margin-right: auto;")
                   )
                 ),
         
         
  fluidRow(
    
    column(4,
      wellPanel(
        bsPopover("pop1",
                  "Input PDB",
                  "Enter a 4-letter PDB code of interest. The normal modes will be calculated on the selected chain ID(s). </br></br>To continue analysis proceed by navigating through the <b>NEXT</b> buttons.", 
                  placement = "right", trigger = "hover",
                  options = list(container = "body")),
        
        tags$div(id = "pop1", icon("question-circle"),
                 style = "position: absolute; right: 25px; top: 5px;"
                 ),
        
         
        h4("A) PDB Input Selection"),
        tags$hr(),

        helpText("Please enter either a single PDB code of interest, a single protein sequence, or multiple related PDB codes (see the ",
                 a(href="http://thegrantlab.org/bio3d", target="_blank", "Help"), " page for more details)."),
        
        ##-PDB input
        textInput("pdbid", label="Enter RCSB PDB code/ID:", value = "4Q21"),
        
        ##- Chain selection
        uiOutput("chain_input"),

        conditionalPanel(
          condition = "output.pdb_isok == false",

          div(
            p(style = "color: red;",
              strong("Error:"),
              "PDB with < 10 or > 600 c-alpha atoms not allowed"
              )
            )
          ),
          
        actionButton3("next-btn-1", "Next (Results)", icon=icon("arrow-down"), cl="btn btn-primary btn-input action-button"),
        
        tags$script(HTML(
          '$(".btn-input").click(function(){',
          '$("html, body").animate({scrollTop:$("#resultsdiv").position().top - (0.1 * $(window).height())}, "smooth");',

          
          '$("#fluct_plot").parent().siblings().find(".well").addClass("show-border");',
          'window.setTimeout(function(){',
          '$("#fluct_plot").parent().siblings().find(".well").removeClass("show-border");',
          '}, 2500);',
          '});'
          
          )),

        

        ## reset PDB input
        actionButton("reset_pdbid", "Reset PDB input", icon=icon("undo"))
        )
   ),

    column(4,
      wellPanel(
        style="background: #FFFFFF;",
        
        bsPopover("pop2",
                  "Force field",
                  "Five popular elastic network model (ENM) force fields are available. <br><br>The <b>calpha</b> force field - originally developed by Konrad Hinsen - is the recommended one for most applications. It employs a spring force constant differentiating between nearest-neighbour pairs along the backbone and all other pairs. The force constant function was parameterized by fitting to a local minimum of a crambin model using the AMBER94 force field.", 
                  placement = "right", trigger = "hover",
                  options = list(container = "body")),
        
        tags$div(id = "pop2", icon("question-circle"),
                 style = "position: absolute; right: 25px; top: 5px;"
                 ),
        
        h4("B) NMA parameters"),
             tags$hr(),
             
             helpText("Use the drop down menu to change force field. The C-alpha force field is recommended for most applications. "),
             

        selectInput("forcefield", "Choose a forcefield:",
                    choices = c("calpha", "sdenm", "reach", "anm", "pfanm")),

        conditionalPanel(
          condition = "input.forcefield == 'anm' || input.forcefield == 'pfanm' ",
          sliderInput("cutoff", "Cutoff value:",
                      min = 7, max = 50, value = 15)
        
          ##helpText("Note: Cutoff applies only to 'ANM' and 'pfANM'. Recommended values are 15 and 50 Å, respectively."),
          ),
 
        actionButton("reset_nma_input", "Reset NMA inputs", icon=icon("undo"))
        )
           ),

    column(4,
           wellPanel(
             style="background: #FFFFFF;",
             
             conditionalPanel(
               condition = "input.logviewer == 'pdb'",

               h4("Input Structure Visualization"),
               webGLOutput('pdbWebGL'),
               
               selectInput("view_inpdb_as", "Display options:",
                           c("C-alpha Trace" = "calpha",
                             "Overview" = "overview"),
                           selected = "calpha",
                           multiple = FALSE),
               
               selectInput("view_inpdb_col", "Color options:",
                           c("Secondary structure elements" = "sse",
                             "Residue Index" = "index"),
                           multiple=FALSE)
               
               ,bsTooltip("view_inpdb_as",
                          title="Change the structure view and coloring options. </br></br>Click and drag on the display to rotate, middle mouse button to zoom.",
                          placement = "left", options = list(container = "body"))              
               
               ),

             conditionalPanel(
               condition = "input.logviewer == 'pdblog'",
               h4("Input PDB Read Log"),
               verbatimTextOutput('pdb_log')
               ),
             
             conditionalPanel(
               condition = "input.logviewer == 'bio3d'",

               h3("Normal mode analysis with Bio3D"),
               p(HTML("This <a href=\"http://thegrantlab.org/bio3d/index.php\">Bio3D</a> WebApp provides a rapid and rigorous tool for normal mode analysis of protein structures. Options include multiple popular elastic network models (ENMs), as well as enhanced analyses including residue fluctuations, mode visualization, dynamic cross-correlations, and overlap analysis.")),
               p(HTML("Start by entering a PDB code of interest then proceed by navigating through following the <b>NEXT</b> buttons.")),
               
               img(src="./images/geostas_250x182.png",
                   width=250, style="display: block; margin-left: auto; margin-right: auto;")
               ),
             
             radioButtons("logviewer", "View options:",
                          c("3D structure" = "pdb",
                            "PDB processing log" = "pdblog",
                            "App Info" = "bio3d"),
                          selected="bio3d",
                          inline=TRUE)
             
             )
           )
    ),



  #########################
  ##-- Results Section --##
  #########################
  tags$div(id = "resultsdiv"),
  h2("Results"),
  hr(),

  ##-A. Fluctuations Panel
  fluidRow(
    id = "fluct_row",
    column(4,
           wellPanel(
             bsPopover("pop3",
                       "Residue fluctuations",
                       "The fluctuations are calculated based on all 3N-6 normal modes (where N is the number of atoms). Magnitudes should be used by care as normal mode vectors are by definition without magnitude. ",
                       placement = "right", trigger = "hover",
                       options = list(container = "body")),
             
             tags$div(id = "pop3", icon("question-circle"),
                      style = "position: absolute; right: 25px; top: 5px;"
                      ),
             

             
             h4('Residue fluctuations'),
             checkboxInput('fluxs3', 'Show mode fluctuations', value=TRUE),
             checkboxInput('fluxs2', 'Show B-factors', value=FALSE),
             selectInput('mode_inds', 'Choose Mode indices:',
                         choices=c("all", 7:50), selected="all", multiple=TRUE),
             checkboxInput('show_options1', 'More options', value=FALSE),

             downloadButton('fluctplot2pdf', "Download Plot PDF"),

             br(),
             
             actionButton3("next-btn-2", "Next (Visualize)", icon=icon("arrow-down"),
                           cl="btn btn-primary btn-fluct action-button",
                           style = "margin-top: 5px; width: 100%;"),

             tags$script(HTML(
               '$(".btn-fluct").click(function(){',
               '$("html, body").animate({scrollTop:$("#visrow").offset().top - (0.1 * $(window).height()) }, "smooth");',
               '});'
               ))
             
             )
           ),
    column(8,
           plotOutput("fluct_plot"),

           tags$script(HTML(
             'var css = document.createElement("style");
              css.type = "text/css";
              css.innerHTML = ".nv-x .nv-axislabel { font-size: 20px; }";
              document.body.appendChild(css);
              css = document.createElement("style");
              css.type = "text/css";
              css.innerHTML = ".nv-y .nv-axislabel { font-size: 20px; }";
              document.body.appendChild(css);
              document.getElementById("fluct_plot").focus();
              '
             ))
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
             id = "visrow",

             bsPopover("pop4",
                       "Normal mode visualization",
                       "Visualize the normal modes by toggeling the the <b>Show NM Trajectory</b> checkbox. <br><br>Download buttons enable visualization of the motions described by the principal component in external viewers such as PyMOL or VMD.",
                       placement = "right", trigger = "hover",
                       options = list(container = "body")),
             
             tags$div(id = "pop4", icon("question-circle"),
                      style = "position: absolute; right: 25px; top: 5px;"
                      ),
             
             h4('Normal Modes Visualization'),
             ##checkboxInput('show_trj2', 'Show NM Trajectory', value=FALSE),
             
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
             downloadButton('nma2pymol', label='Download PyMOL vector field'),

             br(), 
             actionButton3("next-btn-3", "Next (DCCM)", icon=icon("arrow-down"),
                           cl="btn btn-primary btn-visu action-button",
                           style = "margin-top: 5px; width: 100%;"
                           ),
             

             tags$script(HTML(
               '$(".btn-visu").click(function(){',
               '$("html, body").animate({scrollTop:$("#dccmrow").offset().top - (0.1 * $(window).height()) }, "smooth");',
               '});'
               ))
             )
           ),
    
    column(8,
           conditionalPanel(
             condition = "output.pdb_isok",
             webGLOutput('nmaWebGL'),
             
             tags$script(HTML(
               'var css = document.createElement("style");
                css.type = "text/css";
                css.innerHTML = ".nv-x .nv-axislabel { font-size: 20px; }";
                document.body.appendChild(css);
                css = document.createElement("style");
                css.type = "text/css";
                css.innerHTML = ".nv-y .nv-axislabel { font-size: 20px; }";
                document.body.appendChild(css);
                document.getElementById("nmaWebGL").focus();
                '
               )
              )
             
             )
           )
    ),
         

  ##-C. DCCM Panel
  fluidRow(
    column(4,
           wellPanel(
             id="dccmrow",

             bsPopover("pop5",
                       "Cross correlations",
                       "The extent to which the atomic fluctuations/displacements of a system are correlated with one another can be assessed by examining the magnitude of all pairwise cross-correlation coefficients. <br><br> This panel calculates a matrix of all atom-wise cross-correlations whose elements are displayed in a graphic representation frequently termed a dynamical cross-correlation map, or DCCM.",
                       placement = "right", trigger = "hover",
                       options = list(container = "body")),
             
             tags$div(id = "pop5", icon("question-circle"),
                      style = "position: absolute; right: 25px; top: 5px;"
                      ),
             
             h4('Cross correlation analysis'),

             #conditionalPanel(
             #  condition = "output.pdb_isok",
             hr(),
             actionButton("run_dccm", "Calculate correlations", icon=icon("gears"),
                          style = "display: block; margin-left: auto; margin-right: auto;", 
                          class = "btn btn-success"),
             hr(),
             
             
             conditionalPanel(
               condition = "output.dccm_isdone",
               checkboxInput('contourplot', 'Contourplot', value=TRUE),
               checkboxInput('sse', 'Show SSE', value=TRUE),
               checkboxInput('colorkey', 'Colorkey', value=TRUE),
               
               sliderInput("height2", "PDF height:",
                           min = 4, max = 12, value = 5, step=1),
               sliderInput("width2", "PDF width:",
                           min = 4, max = 12, value = 7, step=1),
             
               ##downloadButton('dccm2py', "Download PyMOL Visualization script"),
               downloadButton('dccm2pymol', "Download PyMOL session"),
               downloadButton('dccmplot2pdf', "Download Plot PDF"),

               br(),

               actionButton3("next-btn-4", "Next (Overlap)", icon=icon("arrow-down"),
                             cl="btn btn-primary btn-dccm action-button",
                             style = "margin-top: 5px; width: 100%;"),

               tags$script(HTML(
                 '$(".btn-dccm").click(function(){',
                 '$("html, body").animate({scrollTop:$("#overlaprow").offset().top - (0.1 * $(window).height()) }, "smooth");',
                 '});'
                 ))
               )
             )
           ),
    column(8,
           conditionalPanel(
             condition = "output.dccm_isdone",
             plotOutput("dccm_plot")
             )
           )
    ),     
   br(),br(),

         
  ##-C. Overlap analysis
  fluidRow(
    column(4,
           wellPanel(
             id="overlaprow",

             
             bsPopover("pop6",
                       "Overlap analysis",
                       "This panel will search the PDB for structures with a similar sequence as the input PDB (hit the <b>Launch PDB SEARCH</b> button). When identified, select a structure, and hit the <b>Calculate overlap</b> button. <br><br>The overlap, or dot product, measures the similarity between a normal mode vector and a vector describing the conformational differenec between two structures. An overlap of 1 corresponds to identical vectors, while an overlap of 0 correspond to orthogonal vectors.",
                       placement = "right", trigger = "hover",
                       options = list(container = "body")),
             
             tags$div(id = "pop6", icon("question-circle"),
                      style = "position: absolute; right: 25px; top: 5px;"
                      ),
             
             h4('Overlap analysis '),

             conditionalPanel(
               condition = "output.blast_isdone == false",
               helpText("Press the", strong("Launch PDB SEARCH"),
                        "button to enable overlap analysis")  
               ),
             
             conditionalPanel(
               condition = "output.blast_isdone == true",
               helpText("Great! Now select a PDB ID from the table below and press the ",
                        strong("Calculate overlap"), " button")
               ),
             

             hr(),
             
             div(
               style = "margin: 0px auto 0px auto; text-align: center;",
               
               actionButton("run_blast", "1) Launch PDB SEARCH", icon=icon("search"),
                            disabled = FALSE,
                            class = "btn btn-success"),
               
               actionButton("run_overlap", "2) Calculate overlap", icon=icon("gears"),
                            disabled = TRUE,
                            class = "btn btn-success")
               ),

             hr(),

             conditionalPanel(
               condition = "output.blast_isdone",
               DT::dataTableOutput('blast_table')
               ),
             
             conditionalPanel(
               condition = "output.overlap_isdone",
               downloadButton('overlapplot2pdf', "Download Plot PDF"),
               checkboxInput('show_options2', 'More options', value=FALSE)
               )
             )
           ),
    column(8,
           conditionalPanel(
             condition = "output.overlap_isdone",
             plotOutput("overlap_plot")
             )
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
      )

         
   ##br(),br(),

         
  ##-F. Domain analysis
  #fluidRow(
  #  column(4,
  #         wellPanel(
  #           h4('Domain analysis '),
  #           actionButton("run_geostas", "Run Geostas", icon=icon("cog")),
  #           ##checkboxInput('domains', 'Show domain analysis', value=FALSE),
  #           sliderInput("ndomains", "Number of domains:",
  #                       min = 2, max = 10, value = 3),
  #           sliderInput("nmodes", "Number of modes:",
  #                       min = 1, max = 5, value = 3),
  #           downloadButton('geostas2zip', "Download PDB Trajectory")
  #           )
  #         ),
  #
  #  column(8, 
  #         #conditionalPanel(
  #         #  condition = "input.domains == true",
  #           webGLOutput('geostasWebGL')
  #         #  )
  #         )
  #  )
)
