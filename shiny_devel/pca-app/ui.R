library(shiny)
library(rCharts)

## Could open some saved results as an RData file here for first time display 
## Calculations would then only run once single SUBMIT button was pressed.

shinyUI(navbarPage("Bio3D",
                   source("ui/tab_blast.R", local=TRUE)$value,
                   source("ui/tab_align.R", local=TRUE)$value,
                   source("ui/tab_pca.R", local=TRUE)$value,
                   source("ui/tab_enma.R", local=TRUE)$value
                   )
        )
         
