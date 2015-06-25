tabPanel("2. ALIGN",
         icon=icon("arrow-right"),
         tags$style(type="text/css", "body {padding-top: 80px;}"),

         
         fluidRow(
           column(5,
                  wellPanel(
                    h4("A) Alignment summary"),
                    
                    verbatimTextOutput("alignment_summary"),
                    verbatimTextOutput("missres_summary")
                  
                    )
                  ),

           column(7,
                  wellPanel(
                    style="overflow: auto; background: #FFFFFF;",

                            
                    modalBox(id="1", button_label = "Help ", icon = "question",
                             heading="Sequence and structure alignment",
                             content = tags$div(
                               HTML("<p>In this tab the collected structures are superimposed on each other either based on the <strong>identified invariant core</strong>, or on all C-alpha atoms. The invariant core is the region ...</p>"),
                               
                               p("In this panel you can perform simple structure analysis such as calculating all pair-wise RMSD values ... ")
                               )
                             ),
                    
                    
                    h4("B) Edit alignment (optional)"),
                    
                    radioButtons("toggle_editing", "Action",
                                 c("Filter structures" = "filter",
                                   "Upload alignment" = "upload"),
                                 inline = TRUE),
                    
                    conditionalPanel(
                      condition = "input.toggle_editing == 'filter'",
                      uiOutput("include_hits"),
                      
                      checkboxInput("omit_missing",
                                    "Omit PDBs with missing in-structure residues",
                                    value=FALSE)
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
           ),
         
         
         
         
         hr(),
         fluidRow(
           column(12,
                  
                  h2("Final alignment"),
                  p("(slow for large aligments)"),
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
           ),
         
          fluidRow(
           column(4,
                  wellPanel(
                    style="overflow: auto;",
                    
                    h4("Sequence analysis"),

                    radioButtons("seq_plot", "Plot options",
                                 c("Heatmap" = "heatmap",
                                   "Dendrogram" = "dendrogram"),
                                 inline=TRUE),

                    sliderInput("clusters_seq", "Cluster by pairwise sequence identity",
                                min = 1, max = 10, value = 3, step=1),
                    
                    checkboxInput('show_options0', 'More options', value=FALSE),
                    
                    conditionalPanel(
                      condition = "input.seq_plot == 'heatmap'",
                      downloadButton('seqide_heatmap2pdf', "Download PDF")
                      ),

                    conditionalPanel(
                      condition = "input.seq_plot == 'dendrogram'",
                      downloadButton('seqide_dendrogram2pdf', "Download PDF")
                      ),

                    downloadButton('seqideZIP', "Download Seq ide matrix")
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
               )

             )
           ),

         conditionalPanel(
           condition = "input.show_options0 == true",
           fluidRow(
             column(4,
                    wellPanel(
                    sliderInput("cex0", "Label size",
                                min = 0.1, max = 3, value = 1, step=0.1),
                    sliderInput("margins0", "Plot margins",
                                min = 3, max = 10, value = 5, step=1)
                    )
                    ),
             column(4,
                    wellPanel(
                    sliderInput("width0", "width",
                                min = 4, max = 12, value = 7, step=0.5),
                    sliderInput("height0", "height",
                                min = 4, max = 12, value = 7, step=0.5)
                    )
                    )
             )
           )


         
         )

