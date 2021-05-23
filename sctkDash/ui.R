#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinydashboard)
library(shinydashboardPlus)
library(shinyBS)
library(shinyWidgets)
library(shiny)
library(shinyjs)
library(shinyFiles)
library(ComplexHeatmap)
library(limma)
library(ggplot2)
library(plotly)
library(data.table)
library(colourpicker)
library(gridExtra)
library(cluster)
library(ggtree)
library(ape)
library(GSVA)
library(GSVAdata)
library(shinyalert)
library(enrichR)
library(matrixStats)
library(Biobase)
library(base)
library(SingleCellExperiment)
library(singleCellTK)
library(celda)
library(shinycssloaders)
library(shinythemes)
library(shinyBS);
library(shinyjqui);
library(Seurat);
library(ggplotify);
library(ggplot2);
library(cowplot);
library(tidyverse)
library(dplyr)
library(readxl)
library(broom)
library(RColorBrewer)
library(grDevices)
library(shinyWidgets)
library(stringr)
library(Hmisc)

# Define UI for application that draws a histogram
shinyUI(
    dashboardPage(
        title = "Box API",
        dashboardHeader(),
        sidebar = dashboardSidebar(
            sidebarMenu(
                menuItem("Dashboard", tabName = "dashboard", icon = icon("dashboard")),
                menuItem("Widgets", icon = icon("th"), tabName = "widgets",
                         badgeLabel = "new", badgeColor = "green")
            )
        ),
        body = dashboardBody(
            tabItems(
                tabItem(tabName = "dashboard",
                        h2("Dashboard tab content"),
                        fluidRow(
                            tabBox(
                                title = "",
                                # The id lets us use input$tabset1 on the server to find the current tab
                                id = "tabset1", height = "250px", width = NULL,
                                tabPanel("Feature Selection", "First tab content",
                                         fluidRow(
                                             box(title = "Plot", plotOutput("abc"), background = "gray", sidebar = boxSidebar(
                                                 id = "mycardsidebar",
                                                 width = 25,
                                                 sliderInput(
                                                     "obs", 
                                                     "Number of observations:",
                                                     min = 0, 
                                                     max = 1000, 
                                                     value = 500
                                                 )
                                             )),
                                             box(title = "Options", textInput("asd", "Label:"), actionButton("abcs", "Run"), background = "gray", collapsible = TRUE)
                                         )),
                                tabPanel("Dimensionality Reduction", "Tab content 2")
                            )
                        )
                ),
                
                tabItem(tabName = "widgets",
                        h2("Widgets tab content"),
                        bsCollapse(id = "SeuratUI", open = "Data Input",
                                   bsCollapsePanel("Normalize Data"),
                                   bsCollapsePanel("Scale Data"),
                                   bsCollapsePanel("Highly Variable Genes"),
                                   bsCollapsePanel("Dimensionality Reduction",
                                                   tabsetPanel(type = "tabs",
                                                               tabPanel("PCA",
                                                                        br(),
                                                                        fluidRow(
                                                                            column(4,
                                                                                   fluidRow(
                                                                                       column(12,
                                                                                              panel(heading = "PCA",
                                                                                                    textInput(inputId = "pca_no_components", label = "Select number of components to compute: ", value = "50"),
                                                                                                    materialSwitch(inputId = "pca_compute_elbow", label = "Compute ElbowPlot?", value = TRUE),
                                                                                                    materialSwitch(inputId = "pca_compute_jackstraw", label = "Compute JackStrawPlot?", value = FALSE),
                                                                                                    materialSwitch(inputId = "pca_compute_heatmap", label = "Compute Heatmap?", value = TRUE),
                                                                                                    conditionalPanel(
                                                                                                        condition = 'input.pca_compute_heatmap == true',
                                                                                                        numericInput(inputId = "pca_compute_heatmap_nfeatures",
                                                                                                                     label = "Set number of features for heatmap:", value = 30, step = 1),
                                                                                                    ),
                                                                                                    actionButton(inputId = "run_pca_button", "Run PCA")
                                                                                              ),
                                                                                              panel(heading = "Select No. of Components",
                                                                                                    htmlOutput(outputId = "pca_significant_pc_output", inline = FALSE),
                                                                                                    numericInput(inputId = "pca_significant_pc_counter", label = "Select number of components for downstream analysis: ", min = 1, max = 20, value = 10)
                                                                                              )
                                                                                       )
                                                                                   )
                                                                            ),
                                                                            column(8,
                                                                                   fluidRow(
                                                                                       column(12,
                                                                                              hidden(
                                                                                                  tags$div(class = "seurat_pca_plots", tabsetPanel(id = "seuratPCAPlotTabset", type = "tabs"
                                                                                                  )
                                                                                                  ))
                                                                                       )
                                                                                       
                                                                                   )
                                                                            )
                                                                        )
                                                                        
                                                               ),
                                                               tabPanel("ICA",
                                                                        br(),
                                                                        fluidRow(
                                                                            column(4,
                                                                                   fluidRow(
                                                                                       column(12,
                                                                                              panel(heading = "ICA",
                                                                                                    textInput(inputId = "ica_no_components", label = "Select number of components to compute: ", value = "20"),
                                                                                                    materialSwitch(inputId = "ica_compute_heatmap", label = "Compute Heatmap?", value = TRUE),
                                                                                                    conditionalPanel(
                                                                                                        condition = 'input.ica_compute_heatmap == true',
                                                                                                        numericInput(inputId = "ica_compute_heatmap_nfeatures",
                                                                                                                     label = "Set number of features for heatmap:", value = 30, step = 1),
                                                                                                    ),
                                                                                                    actionButton(inputId = "run_ica_button", "Run ICA")
                                                                                              ),
                                                                                              panel(heading = "Select No. of Components",
                                                                                                    #h5("Number of components suggested by ElbowPlot: "),
                                                                                                    #verbatimTextOutput(outputId = "ica_significant_pc_output", placeholder = TRUE),
                                                                                                    numericInput(inputId = "ica_significant_ic_counter", label = "Select number of components for downstream analysis: ", min = 1, max = 20, value = 10)
                                                                                              )
                                                                                       )
                                                                                   )
                                                                            ),
                                                                            column(8,
                                                                                   fluidRow(
                                                                                       column(12,
                                                                                              hidden(
                                                                                                  tags$div(class = "seurat_ica_plots", tabsetPanel(id="seuratICAPlotTabset", type = "tabs"
                                                                                                  ))
                                                                                              )
                                                                                       )
                                                                                   )
                                                                            )
                                                                        )
                                                                        
                                                               )
                                                   ),
                                                   style = "primary")
                        )
                )
            )
        )
        
    )
)


