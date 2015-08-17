options(rgl.useNULL=TRUE)

library(shiny)
library(DT)
library(shinyRGL)
library(rCharts)
library(threejs)
library(shinyBS)
library(shinyjs)

source('ui/ui_utils.R')

shinyUI(
   
  navbarPage(windowTitle="Bio3D PCA App",
             title=img(src="./images/bio3d_logo.png", height=35, alt="Bio3D"),
             theme = "bootstrap.css",
             position = "fixed-top",

             ##useShinyjs(),
             source("ui/tab_blast.R", local=TRUE)$value,
             source("ui/tab_align.R", local=TRUE)$value,
             source("ui/tab_fit.R", local=TRUE)$value,
             source("ui/tab_pca.R", local=TRUE)$value,
             source("ui/tab_enma.R", local=TRUE)$value,             
             ##source("ui/tab_help.R", local=TRUE)$value,

             includeCSS("www/styles.css")
             
             )
  )
         
