tabPanel(
  "4. PCA", icon=icon("arrow-right"),
  tags$style(type="text/css", "body {padding-top: 80px;}"),
  div(
    h3("Principal Component Analysis", style="border: 1px solid #e3e3e3; border-radius: 4px; padding-top: 10px; padding-bottom: 10px; padding-left: 5px; margin-top: -10px; background-color: white;")
    ),
  
  actionButton3("about_pcatab", "About this tab", icon=icon("comment"), cl="btn btn-warn btn-input action-button", style = "position: fixed; top: 14px; right: 16px; z-index: 2000;"),
  
  bsModal("modal_pca", "Principal Component Analysis", "about_pcatab", size = "large", 
          content=tags$div(
            
            p(HTML("In the PCA tab, a principal component analysis (PCA) is performed based on the coordinates of the superimposed structures in the FIT tab. PCA is an effective approach to capture and characterize inter-conformer relationships. The motions described by the principal components (PCs) can be visualized either in the browser, or by downloading a PDB trajectory file or a PyMOL session (.pse) containing a vector field representation.  ")),
            
            img(src="./images/aaah-pc1.gif", width=500, style="display: block; margin-left: auto; margin-right: auto;")                      
            ),
          
          p(HTML("The figure above shows the first principal component obtained from the 32 PDB structures of the enzyme family <i>aromatic amino acid hydroxylases</i>. The analysis reveals the presence of a sub-domain within the catalytic domain which obtains a more compact conformation upon substrate binding.")),

          br(),
          img(src="./images/pca-conformerplot.png", style="display: block; margin-left: auto; margin-right: auto;"),
          
          p(HTML("The figure above shows the conformer plot of 90 <i>E. coli</i> DHFR structures. The PCA reveals that the ensemble can be divided into three major groups along their first two principal components. These conformers display either a closed, occluded,oranopenconformationof two active site loops. The conformer plot displays a low-dimensional representation of conformational variability described by two selected PCs (typically PC-1 and PC-2). "))
          ),
  

  ##- PC Visualization 
  fluidRow(
    column(4,
           wellPanel(
             bsPopover("poppca1",
                       "Visualization",
                       "Visualize the motions described by the individual principal components by toggeling the the <b>Show PC Trajectory</b> checkbox. <br><br>Download buttons enable visualization of the motions described by the principal component in external viewers such as PyMOL or VMD.",
                       placement = "right", trigger = "hover",
                       options = list(container = "body")),
             
             tags$div(id = "poppca1", icon("question-circle"),
                      style = "position: absolute; right: 25px; top: 5px;"
                      ),
             
             
             h4('A) Principal component visualization'),
             checkboxInput('show_trj', 'Show PC Trajectory', value=FALSE),

             selectInput('viewPC', 'Choose Principal Component:', choices=c(1:10)),

             sliderInput("mag1", "Magnification factor:",
                         min = 1, max = 10, value = 1),
             
             selectInput('viewColor', 'Color options',
                             choices=list(
                               'Amalgam' = 'amalgam',
                               'Magnitude'='mag',
                               'By Frame (blue->gray->red)'='default'
                               ),
                             selected='amalgam'
                          ),
             
             #radioButtons('viewBGcolor', label='Background color',
             #             choices=list('Black'='black', 'White'='white'),
             #             selected='white'),
             
             selectInput("viewBGcolor", "Background color:",
                         c('White'='white', 'Black'='black'),
                         selected = 'white'),
             
             br(),
             actionButton('viewUpdate', label='Refresh', icon=icon('undo')),
             downloadButton('pctraj', label='Download PDB Trajectory'),
             downloadButton('pca2pymol', "Download PyMOL session file"),
             br(),
             actionButton3("next-btn-pca1", "Next (Plot)", icon=icon("arrow-down"), cl="btn btn-primary btn-input action-button", style = "margin-top: 0.9%;"),
             tags$script(HTML(
               '$("#next-btn-pca1").click(function(){',
               '$("html, body").animate({scrollTop:$("#pca_conformer_row").position().top - (0.1 * $(window).height())}, "smooth");',
               'well = $("#pca_conformer_row").children().find(".well");',
               'well.addClass("show-border");',
               'window.setTimeout(function(){',
               'well.removeClass("show-border");',
               '}, 2500);',
               '});'
             ))
             )
           ),

    column(8,
           conditionalPanel(
             condition='input.show_trj == true',
             webGLOutput('pcaWebGL')
             )
           )
  ),


  ##- Conformer plots
  fluidRow(
    id = 'pca_conformer_row',
    column(4,
           wellPanel(
             bsPopover("poppca2",
                       "Conformer plot",
                       "A conformer plot is a low-dimensional representation of the conformational variability within the ensemble of PDB structures. The plot is obtained by projecting individual structures onto two or three selected principal components. ",
                       placement = "right", trigger = "hover",
                       options = list(container = "body")),

             tags$div(id = "poppca2", icon("question-circle"),
                      style = "position: absolute; right: 25px; top: 5px;"
                      ),

        h4("B) Conformer plot"),
             popSelectInput('pcx', 'PC on X-axis:', choices=c(1:10), selected=1),
        selectInput('pcy', 'PC on Y-axis:', choices=c(1:10), selected=2),

        conditionalPanel(
          condition = "input.plot_type == '3dscatter1' || input.plot_type == '3dscatter2'",
          selectInput('pcz', 'PC on Z-axis:', choices=c(1:10), selected=3)
          ),

        radioButtons('cluster_by', label='Cluster by',
                     choices=list(
                       'RMSD' = 'rmsd',
                       'PC subspace' = 'pc_space',
                       'Sequence' = 'sequence'
                       ),
                     selected='rmsd', inline=TRUE),

        conditionalPanel(
          condition = "input.cluster_by == 'pc_space'",

          ## K-selecter
          uiOutput("kslider_pca"),
          actionButton("setk_pca", "Auto set number of K groups",
                       icon=icon("cogs")),

          br(),br(),
          sliderInput("clust_npcs", "PCs in subspace",
                      min = 1, max = 10, value = 2)
          ),

        radioButtons("plot_type", "Plot type",
                     c("2D scatter" = "normal",
                       "3D scatter" = "3dscatter1",
                       "Interactive" = "fancy"
                       ##"3D scatter (three-js)" = "3dscatter2",
                       ),
                     inline=TRUE),

        conditionalPanel(
          condition = "input.plot_type == 'normal'",
          checkboxInput("show_options1", "More clustering and output options", value=FALSE)
          ),

        conditionalPanel(
          condition = "input.plot_type == '3dscatter2'",
          selectInput("renderer", label="Rendering method",
                      choices = list("Auto"="auto", "Canvas"="canvas", "WebGL"="webgl"),
                      selected = 1),
          checkboxInput("grid", label = "Grid", value = TRUE)
          ),
          actionButton3("next-btn-pca2", "Next (Residue contributions)", icon=icon("arrow-down"), cl="btn btn-primary btn-input action-button"),
          tags$script(HTML(
            '$("#next-btn-pca2").click(function(){',
            '$("html, body").animate({scrollTop:$("#pca_contrib_row").position().top - (0.1 * $(window).height())}, "smooth");',
            'var well = $("#pca_contrib_row").children().find(".well");',
            'well.addClass("show-border");',
            'window.setTimeout(function(){',
            'well.removeClass("show-border");',
            '}, 2500);',
            '});'
            ))
        )
           ),

    conditionalPanel(
      condition = "input.plot_type == '3dscatter1'",
      column(width=8,
             webGLOutput("scatterplot3d_webgl")
             )
      ),

    conditionalPanel(
      condition = "input.plot_type == '3dscatter2'",
      column(width=8,
             scatterplotThreeOutput("scatterplot3d_rthreejs")
             )
      ),

    conditionalPanel(
      condition = "input.plot_type != '3dscatter'",

      column(width=4,
             conditionalPanel(
               condition = "input.plot_type == 'normal'",
               plotOutput("pca_plot1_conf")
               ),

             conditionalPanel(
               condition = "input.plot_type == 'fancy'",
               showOutput("pca_plot2_conf","dimple")
               )
             ),

      column(4,
             conditionalPanel(
               condition = "input.plot_type == 'normal'",
               plotOutput("pca_plot1_scree")
               ),
             conditionalPanel(
               condition = "input.plot_type == 'fancy'",
               showOutput("pca_plot2_scree","nvd3")
               )
             )
      )

    ),


  ##- Additional plotting options for conformer plot
  fluidRow(
    conditionalPanel(
      condition = "input.show_options1 == true",
      column(3,
             wellPanel(
               selectInput("hclustMethod_pca", label="Clustering method", 
                           choices=list(
                             "single"="single","complete"="complete","average"="average",
                             "mcquitty"="mcquitty","median"="median","centroid"="centroid",
                             "ward.D"="ward.D","ward.D2"="ward.D2"
                             ),selected="ward.D2"), 
               
               numericInput("minDistance_pca","Minimum branching gap", value = 0.5, step = 0.5)
               )
             )
      ),
    
    conditionalPanel(
    condition = "input.show_options1 == true && input.plot_type == 'normal'",
      column(3,
             wellPanel(
               sliderInput("cex_points", "Point size",
                           min = 0.1, max = 3, value = 1, step=0.1),
               sliderInput("inner_margin", "Scale axes",
                           min = 1, max = 2, value = 1.2, step=0.1)
               )
             ),

      column(3,
             wellPanel(
               checkboxInput("labelplot", "Label plot", value=FALSE),

               conditionalPanel(
                 condition = "input.labelplot == true",

                 checkboxInput("distribute_labels", "Distribute labels", value=FALSE),
                 sliderInput("cex_labels", "Label size",
                             min = 0.1, max = 3, value = 1, step=0.1),

                 sliderInput("offset", "label offset",
                             min = 0, max = 2, value = 0.5, step=0.1)
                 )
               )
             ),

      column(3,
             wellPanel(
               conditionalPanel(
                 condition = "input.labelplot == true",

                 checkboxInput("toggle_all", "Toggle all", TRUE),
                 uiOutput("checkboxgroup_label_ids")
                 )
               )
             )
    )
    ),

  fluidRow(
    conditionalPanel(
      condition = "input.plot_type == 'normal'",
      column(12,
             wellPanel(
               h4("PCA conformer plot annotation"),
               helpText("Highlight structures in conformer plot by clicking their entries in the below table (only for plot type '2D Scatter')."),
               DT::dataTableOutput("pdbs_table")
               )
             )
      )
    ),
           

  
  ##- PC Loadings - fluctuations
  fluidRow(
    id = 'pca_contrib_row',
    column(4,
           wellPanel(
             bsPopover("poppca3",
                       "Residue contributions",
                       "Residue-wise contributions to the chosen principal component.",
                       placement = "right", trigger = "hover",
                       options = list(container = "body")),

             tags$div(id = "poppca3", icon("question-circle"),
                      style = "position: absolute; right: 25px; top: 5px;"
                      ),

             h4('C) Residue contributions'),
             selectInput('loadings_pc', 'Choose Principal Component:',
                         choices=c(1:10), selected=1, multiple=TRUE),
             ##checkboxInput("toggle_rmsf1", "Show RMSF", TRUE),
             checkboxInput("toggle_rmsf1", "Show RMSF", TRUE),
             checkboxInput("multiplot", "Multiline plot", FALSE),
             checkboxInput("spread_pcload", "Spread lines", FALSE),
             downloadButton('pcloadings2pdf', label='Download PDF'),
             checkboxInput("show_options2", "More options", value=FALSE),
             actionButton3("page5", "Next (eNMA)", icon=icon("arrow-right"), cl="btn btn-primary btn-next-blast action-button"),
             tags$script(HTML(
               '$("#page5").click(function(){',
               'tabs = $(".nav.navbar-nav li");',
               'tabs.children()[4].click();',
               '});'
               ))
             )
           ),
    
    column(8,
           plotOutput("loadings_plot")
           )
    ),
  
  conditionalPanel(
    condition = "input.show_options2 == true",
    
    fluidRow(
      column(4,
             wellPanel(
               sliderInput("width_pcload", "Figure width (PDF only)",
                           min = 4, max = 12, value = 8, step=0.5),
               sliderInput("height_pcload", "Figure height (PDF only)",
                           min = 4, max = 12, value = 5, step=0.5)
               )
             )
      )
    )

  ##- Annotation datatable
  #fluidRow(
  #  column(12,
  #         wellPanel(
  #           DT::dataTableOutput("pdbs_table")
  #           )
  #         )
  #  )
  )
