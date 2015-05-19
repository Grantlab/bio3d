tabPanel("2. ALIGN",
         icon=icon("arrow-right"),
         tags$style(type="text/css", "body {padding-top: 80px;}"),

         fluidRow(
           column(5,
                  wellPanel(
                    h4("A) Alignment summary"),
                    tags$hr(),

                    verbatimTextOutput("alignment_summary"),

                    checkboxInput("omit_missing",
                                  "Omit PDBs with missing in-structure residues",
                                  value=FALSE)

                    )
                  ),

           column(7,
                  wellPanel(style="overflow: auto;",

                            h4("B) Edit alignment (optional)"),
                            hr(),

                            div(
                              div(style="float: left; width: 45%;
                                         padding: 10px;",

                                  strong("Download FASTA file"),
                                  helpText("To edit the alignment, download the FASTA alignment file,
                                            and upload in the box to the right. "),

                                  br(),
                                  downloadButton('fastafile', "Download FASTA alignment file"),
                                  actionButton("reset_fasta", "Reset alignment", icon=icon("undo"))

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
                            )
                  )
           ),

         hr(),
         fluidRow(
           column(12,
                  h2("Final alignment"),
                  p("(slow for large aligments)"),
                  uiOutput("alignment")
                  )
           )
         )

