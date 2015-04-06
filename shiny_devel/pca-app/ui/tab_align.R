tabPanel("2. ALIGN", icon=icon("arrow-right"),

         
         fluidRow(
           column(4,
                  wellPanel(
                    h4("Alignment summary"),
                    tags$hr(),

                    verbatimTextOutput("alignment_summary"),

                    checkboxInput("omit_missing",
                                  "Omit PDBs with missing in-structure residues",
                                  value=FALSE)

                    )
                  ),
                    
                    
           column(4,
                  wellPanel(
                    h4("Download Alignment"),
                    hr(),
                    
                    ##downloadButton('pdbsRData', "Download pdbs object"),
                    downloadButton('fastafile', "Download FASTA file"),
                    helpText("Download the FASTA alignment file ... ")
                    )
                  ),
                                        
           column(4,
                  wellPanel(
                    h4("Upload Alignment"),
                    tags$hr(),

                    helpText("Upload manually edited alignment file"),
                    fileInput('fastafile_upload', 'Upload FASTA File',
                              accept=c('text/fasta', 'text/plain', '.fasta')),
                    
                    actionButton("reset_fasta", "Reset alignment")
                    )
                  )
           ),

         hr(),
         fluidRow(
           column(12,
                  h2("Final alignment"),
                  uiOutput("alignment")
                  )
           )
         )
         
