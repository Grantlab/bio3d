
library(shiny)

# Define UI for application that draws a histogram
shinyUI(navbarPage("Bio3D NMA",

                   tabPanel("PDB",
                            sidebarLayout(

                                sidebarPanel(
                                    textInput("pdbid", "Ender a PDB code:", "1hel"),

                                    textInput("chain", "Chain ID:", "A"),

                                    submitButton("Submit")
                                    ),


                                mainPanel(
                                    h4("Summary"),
                                    verbatimTextOutput("pdbSummary")
                                    )
                                )
                            ),

                   tabPanel("NMA",
                            sidebarLayout(
                                sidebarPanel(
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
                                                max = 350,
                                                value = 300),

                                    submitButton("Submit")
                                    ),

                                mainPanel(
                                    h4("Summary"),
                                    verbatimTextOutput("modeSummary")
                                    )
                                )
                            ),

                   tabPanel("Fluctuations",
                            plotOutput("fluctPlot")
                            ),

                     tabPanel("DCCM",
                              sidebarLayout(
                                sidebarPanel(
                                    radioButtons("style", "Style:",
                                                 c("Cyan", "Red-blue")),

                                    checkboxInput('colorkey', 'Color key', TRUE),

                                    checkboxInput('sse', 'Show SSE', TRUE),

                                    submitButton("Submit")
                                    ),

                                  mainPanel(
                                      plotOutput("dccmPlot")
                                      )
                                  )
                              )
                   )
        )


