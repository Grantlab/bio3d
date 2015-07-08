options(rgl.useNULL=TRUE)

library(shiny)
library(shinyRGL)
library(DT)
library(shinyBS)

source('ui/ui_utils.R')


## Could open some saved results as an RData file here for first time display 
## Calculations would then only run once single SUBMIT button was pressed.

##shinyUI(fluidPage(
shinyUI(
  navbarPage(windowTitle="Bio3D NMA",
             title=img(src="./images/bio3d_logo.png", height=35, alt="Bio3D"),
             theme = "bootstrap.css",
             position = "fixed-top",
             
             source("ui/tab_nma.R", local=TRUE)$value,
             source("ui/tab_help.R", local=TRUE)$value,
             includeCSS("www/styles.css")
             
             )
  )
 
