
library(shiny)

# Define UI for application that draws a histogram
shinyUI(fluidPage(

  # Application title
  titlePanel("NMA on your favorite protein!"),

  # Sidebar with a slider input for the number of bins
  sidebarLayout(
    sidebarPanel(
      textInput("pdbid", "Ender a PDB code:", "1hel"),

      selectInput("forcefield", "Choose a forcefield:", 
                  choices = c("calpha", "anm", "pfanm", "sdenm")),

      sliderInput("cutoff",
                  "Cutoff value:",
                  min = 6,
                  max = 50,
                  value = 10),

      radioButtons("mass", "Mass-weighting:",
                  c("Yes", "Nope")),

      sliderInput("temp",
                  "Temperature scaling:",
                  min = 0,
                  max = 500,
                  value = 300),

      radioButtons("plot", "What to plot:",
                  c("Fluctuations", "Cross correlations"))
    ),

    # Show a plot of the generated distribution
    mainPanel(
      plotOutput("distPlot")
    )
  )
))

