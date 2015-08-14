tabPanel("2. ALIGN",
         icon=icon("arrow-right"),
         tags$style(type="text/css", "body {padding-top: 80px;}"),

         div(
           h3("Sequence Alignment", style="border: 1px solid #e3e3e3; border-radius: 4px; padding-top: 10px; padding-bottom: 10px; padding-left: 5px; margin-top: -10px; background-color: white;")
           ),

         
         actionButton3("about_aligntab", "About this tab",
                       icon=icon("comment"),
                       cl="btn btn-warn btn-input action-button",
                       style = "position: fixed; top: 14px; right: 16px; z-index: 2000;"),
         
         bsModal("modal_align", "Sequence Alignment and Analysis", "about_aligntab", size = "large", 
                 content=tags$div(
                   p(HTML("In this tab all PDB structures selected in the previous (SEARCH) tab are parsed and their sequences aligned. The tab displays the final alignment of the selected PDB structures, as well as basic analysis of sequence identity.")),

                   p(HTML("You can optionally exclude structures, either manually (see 'Exlcude / include hits'), or automatically omit structures with missing in-structure residues. Functionality for uploading a corrected / revised sequence alignment (FASTA format) is also provided."))
                   ),
                 
                 p(HTML("The sequence alignment is performed with MUSCLE and is displayed in blocks of up to 80 columns with a row per PDB structure. Conserved columns are annotated with an asterisk (*) while columns containing similar amino acids are depicted with a hat (^). Detailed information (including PDB residue number, and position in alignment) on each residue in the alignment can be obtained by hovering over the the residue. ")), 


                 img(src="./images/alignment.png", width=700, style="display: block; margin-left: auto; margin-right: auto;"),
                   
                 p(HTML("Basic analyses of sequence identity entails clustering analysis, as well as entropy and sequence conservation. These analyses are available as plot dendrograms, heatmap and histogram plot. ")), 
                 img(src="./images/seqide_heatmap.png", width=600, style="display: block; margin-left: auto; margin-right: auto;")
                 
                 ),
         
         fluidRow(
           column(4,
                  wellPanel(
                    bsPopover("popalign3",
                              "Alignment summary", 
                              "Grey bars in the schematic alignment represents non-gap positions. Red bar above the alignment denote sequence conservation (white to red scale). The sequences are grouped/clustered (dendrogram on the left) based on their pair-wise sequence identities. ", 
                              placement = "right", trigger = "hover",
                              options = list(container = "body")),
                    
                    tags$div(id = "popalign3", icon("question-circle"),
                             style = "position: absolute; right: 25px; top: 5px;"
                             ),
                    
                    
                    h4("Alignment summary"),

                    htmlOutput("alignment_summary"),
                    checkboxInput('cluster_alignment', "Alignment overview clustering", value=TRUE),
                    htmlOutput("omitted_pdbs_summary"),
                    htmlOutput("missres_summary"),

                    conditionalPanel(
                      condition = "output.npdbs_with_missres > 0 || input.omit_missing == true",
                      checkboxInput("omit_missing", "Omit PDBs with missing in-structure residues",
                                    value=FALSE)
                      ),

                    hr(),
                    checkboxInput('show_alignment_edit', 'Filter/edit alignment', value=FALSE),
                    
                    actionButton3("next-btn-align1", "Next (Analysis)", icon=icon("arrow-down"), cl="btn btn-primary btn-input action-button"),
        
                    tags$script(HTML(
                      '$("#next-btn-align1").click(function(){',
                      '$("html, body").animate({scrollTop:$("#seqanalysis_row").position().top - (0.1 * $(window).height())}, "smooth");',
                      '$("#seqanalysis_row").children().find(".well").addClass("show-border");',
                      'window.setTimeout(function(){',
                      '$("#seqanalysis_row").children().find(".well").removeClass("show-border");',
                      '}, 2500);',
                      '});'
                      
                      ))


                    )
                  ),

           column(8,
                  plotOutput("schematic_alignment")
                  )
           ),


         conditionalPanel(
           condition = "input.show_alignment_edit == true",
           
           fluidRow(
             column(12,
                    wellPanel(
                      style="overflow: auto; background: #FFFFFF;",
                      
                      
                      h4("Edit alignment"),
                    
                      radioButtons("toggle_editing", "Action",
                                   c("Filter structures" = "filter",
                                     "Upload alignment" = "upload"),
                                   inline = TRUE),
                      
                      conditionalPanel(
                        condition = "input.toggle_editing == 'filter'",
                        uiOutput("include_hits")
                        
                        #checkboxInput("omit_missing",
                        #            "Omit PDBs with missing in-structure residues",
                        #              value=FALSE)
                        ),
                      
                      conditionalPanel(
                        condition = "input.toggle_editing == 'upload'",
                        
                        div(
                          div(style="float: left; width: 45%;
                                         padding: 10px;",
                              
                              strong("Download FASTA file"),
                              helpText("To edit the alignment, download the FASTA alignment file,
                                            and upload in the box to the right. "),
                              
                              br(),
                              downloadButton('fastafile', "Download FASTA alignment file")
                              ),
                          
                          div(style="float: right; width: 45%;
                                         padding: 10px; margin-left: 5px;",
                              
                              strong("Upload FASTA file"),
                              helpText(strong("Important:"), "when editing the FASTA file, do not edit the
                                            sequence identifiers, and do preserve the amino acid sequences."),
                              
                              fileInput('fastafile_upload', '',
                                        accept=c('text/fasta', 'text/plain', '.fasta'))
                              
                            )
                          
                          )
                        ),
                      
                      actionButton("reset_fasta", "Reset alignment", icon=icon("undo"))
                      )
                    )
             )
           ),

         
         fluidRow(
           id = "seqanalysis_row",
           column(4,
                  wellPanel(
                    style="overflow: auto;",
                    
                    bsPopover("popalign2",
                              "Sequence alignment analysis",
                              "Basic sequence analysis entails clustering of the selected PDB structures based on their pair-wise sequence similarity, which can be visualized as a <b>dendrogram</b>, alternatively in combination with a <b>heatmap</b>. The sequence <b>entropy</b> uses the 10-letter alphabet and conserved (low entropy) columns score 1 and diverse (high entropy) columns score 0.",
                              placement = "right", trigger = "hover",
                              options = list(container = "body")),
                    
                    tags$div(id = "popalign2", icon("question-circle"),
                             style = "position: absolute; right: 25px; top: 5px;"
                             ),
                    
                    
                    
                    h4("Sequence alignment analysis"),

                    radioButtons("seq_plot", "Alignment overview plot options",
                                 c("Dendrogram" = "dendrogram",
                                   "Heatmap" = "heatmap",
                                   "Conservation" = "conservation"),
                                 inline=TRUE),

                    ## K-selecter
                    uiOutput("kslider"),
                    actionButton("setk", "Auto set number of K groups",
                                 icon=icon("cogs")),


                    conditionalPanel(
                      condition = "input.seq_plot == 'conservation'",
                      radioButtons("conserv_method", "Method",
                                 c("Similarity" = "similarity",
                                   "Identity" = "identity", 
                                   "Entropy 22" = "entropy22",
                                   "Entropy 10" = "entropy10"),
                                 inline=TRUE)
                      ),
                    

                   
                    checkboxInput('show_options0', 'More clustering and output options', value=FALSE),
                    
                    #conditionalPanel(
                    #  condition = "input.seq_plot == 'heatmap'",
                    #  downloadButton('seqide_heatmap2pdf', "Download Figure (PDF)")
                     # ),

                    #conditionalPanel(
                    #  condition = "input.seq_plot == 'dendrogram'",
                    #  downloadButton('seqide_dendrogram2pdf', "Download Figure (PDF)")
                    #  ),

                    #conditionalPanel(
                    #  condition = "input.seq_plot == 'conservation'",
                    #  downloadButton('conservation2pdf', "Download Figure (PDF)")
                    #  ),
                    
                    #downloadButton('seqideZIP', "Sequence identity matrix"),

                    actionButton3("next-btn-align2", "Next (Alignment)", icon=icon("arrow-down"), cl="btn btn-primary btn-input action-button"),
        
                    tags$script(HTML(
                      '$("#next-btn-align2").click(function(){',
                      '$("html, body").animate({scrollTop:$("#alignment_row").position().top - (0.1 * $(window).height())}, "smooth");',

                      
                      '$("#alignment").addClass("show-border");',
                      'window.setTimeout(function(){',
                      '$("#alignment").removeClass("show-border");',
                      '}, 2500);',
                      '});'
                      
                      ))
                    )
                  ),
            
            column(8,
                   
                   conditionalPanel(
                     condition = "input.seq_plot == 'heatmap'",
                     plotOutput("seqide_heatmap")
                     ),
                   
                   conditionalPanel(
                    condition = "input.seq_plot == 'dendrogram'",
                     plotOutput("seqide_dendrogram")
                     ),
                   
                   conditionalPanel(
                     condition = "input.seq_plot == 'conservation'",
                     plotOutput("conservation")
                     )
                   
                   )
            ),
         
         conditionalPanel(
           condition = "input.show_options0 == true",
           fluidRow(

             column(3,
                    wellPanel(
                      selectInput("hclustMethod", label="Clustering method", 
                                  choices=list(
                                    "single"="single","complete"="complete","average"="average",
                                    "mcquitty"="mcquitty","median"="median","centroid"="centroid",
                                    "ward.D"="ward.D","ward.D2"="ward.D2"
                                    ),selected="ward.D2"), 
                      
                      numericInput("minDistance","Minimum branching gap", value = 0.1, step = 0.05)
                      
                      #numericInput("splitTreeAt","Or split tree at height",value="",min=0,max=100,step=1),
                      #numericInput("splitTreeK","Or split tree into K groups (not implimented!)", value="", min=1, max=100, step=1)
                      )
                    ),
             
             column(3,
                    wellPanel(
                    sliderInput("cex0", "Label size",
                                min = 0.1, max = 3, value = 1, step=0.1),
                    sliderInput("margins0", "Plot margins",
                                min = 3, max = 10, value = 5, step=1)
                    )
                    ),
             
             column(3,
                    wellPanel(
                    sliderInput("width0", "Figure width (PDF output only)",
                                min = 4, max = 12, value = 7, step=0.5),
                    sliderInput("height0", "Figure height (PDF output only)",
                                min = 4, max = 12, value = 7, step=0.5)
                    )
                    ),
            
                      
                                        # ###
             column(3,
                    wellPanel(
                      
                      h4("Download PDF Figures"),
                      downloadButton('seqide_dendrogram2pdf', "Dendrogram (PDF)"),
                      downloadButton('seqide_heatmap2pdf', "Heatmap (PDF)"),
                      
                      downloadButton('conservation2pdf', "Conservation (PDF)"),
                      downloadButton('seqideZIP', "Sequence identity matrix")

                      )
                    )

             ###
             )
           ),
         
         
         
         
         hr(),

         fluidRow(
           id = "alignment_row",
           column(12,
                  h3("Sequence alignment"),
                  helpText("Rendering the alignment might be time consuming for large data sets (e.g. > 50 PDB IDs depending on the sequence lengths)."),
                  
                  radioButtons("show_alignment", "Show alignment",
                               c("Show" = "yes",
                                  "Hide" = "no"),
                               selected = "no",
                               inline = TRUE),
                  
                  conditionalPanel(
                    condition = "input.show_alignment == 'yes'",
                    
                    uiOutput("alignment"),
                    tags$head(tags$script(src='tooltip.js')),
                    tags$head(tags$script(src='popover.js')),
                    tags$head(tags$script(HTML('$(document).ready(function(){
                      $("body")
                      .popover({html: true,
                      selector: "[data-toggle=\'popover\']",
                      container: "body",
                      title: "Residue info",
                      trigger: "hover",
                      delay: { "show": 200, "hide": 40 },
                      placement: "auto bottom"
                      });
                      });
                      ')))
                    )
                    )
           )

         
         )

