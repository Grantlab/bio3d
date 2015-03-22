tabPanel("2. ALIGN",

         
         fluidRow(
           column(12,
                  wellPanel(
                    h4("4. PDBs Alignment"),
                    helpText("TODO: download and upload alignment file for use with read.fasta.pdbs"),
                    helpText("TODO: download pdbs"),
                    helpText("TODO: nicer format of alignment?")
                    )
                  ),
           column(12,
                  verbatimTextOutput("alignment")
                  )
           )
         
         
         )
         
