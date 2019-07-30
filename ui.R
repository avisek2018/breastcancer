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
library(tidyverse)
library(DT)
library(plotly)
library(shinycssloaders)

source("source.R")


ui <- dashboardPage( 
    
    dashboardHeader(title = "ST 558 - Project 3 - Avisek Choudhury", titleWidth=500),
    ## Sidebar content
    dashboardSidebar(
        sidebarMenu(
            menuItem("About", tabName = "about", icon = icon("question")),
            menuItem("Data", tabName = "data", icon = icon("database")),
            menuItem("Unsupervised Learning", tabName = "unsuper", icon = icon("cogs")),
            menuItem("Supervised Learning", tabName = "super", icon = icon("hornbill"), 
                     menuSubItem('k-Nearest Neighbor', tabName = 'knn'), 
                     menuSubItem('Logistic Regression', tabName = 'logistic'),
                     startExpanded = TRUE)
        )
    ),
    
    ## Body content
    dashboardBody(
        
        tabItems(
            
            # About the Application Tab
            # First tab content
            tabItem(tabName = "about",
                    fluidRow(
                        #add in latex functionality if needed
                        withMathJax(),
                        
                        #two columns for each of the two items
                        column(6,
                               #box to contain description
                               box(style = "background-color:#3c8dbc; color:white",width=12,
                                     h3(shiny::tags$u("About the Dataset")),
                                    includeHTML("dataset.html"),
                                    h3(shiny::tags$u("About the Application")),
                                    includeHTML("application.html")
                        )
                        ),
                        column(6,
                               #How to use the app
                               #h3("About the Application"),
                               #box to contain description
                               box(style = "background-color:#3c8dbc; color:white",width=12,
                                   includeHTML("theory.html")
                               )
                        )
                    )
            ),
            
            # Second tab content - Data Tab
            tabItem(tabName = "data",
                    fluidRow(
                        column(3,
                               box(width=12,title="Select predictors for Visualization",
                                   
                                   selectizeInput("predictor", "Predictor", selected = "Clump_Thickness", 
                                                  choices = levels(as.factor(names(select(data, -Class))))),
                                   
                                   selectizeInput("diagram", "Diagram", selected = "Bar Chart", 
                                                  choices = c("Bar Chart", "Histogram", "BoxPlot")),
                                   conditionalPanel(condition = "input.diagram == 'Histogram' || input.diagram == 'Bar Chart'",
                                                    selectizeInput("position", "Position", selected = "stack",
                                                                   choices = c("dodge","stack"))),
                                   conditionalPanel(condition = "input.diagram == 'Histogram'",
                                                    checkboxInput("density", h6("Overlay Density Plot", style = "color:blue;")))
                               ),
                               box(width = 12, title = "Five Point Number Summary",
                                   verbatimTextOutput("summary")),
                               box(width = 12, title = "Classification",
                                   h4("1 = Malignant"),
                                   h4("0 = Benign"))
                        ),
                           
                        column(9,
                               tabsetPanel(type = "tabs",
                                           #Tab for plotting of data.
                                           tabPanel("Data Visualization", icon = icon("chart-line"),
                                                    fluidRow(
                                                        withSpinner(plotlyOutput("plot1"), type = 1), 
                                                        br(),
                                                        box(width = 7, title = "Boxplot of All Predictors",
                                                            withSpinner(plotlyOutput("allboxplot"), type = 7)), 
                                                        box(width = 5, title = "Correlation Plot",
                                                            withSpinner(plotlyOutput("corrplot"), type = 7))
                                                    )
                                           ), 
                                           #Tab for display of data.
                                           tabPanel("Data Display", icon = icon("chalkboard"),
                                                    fluidRow(
                                                        br(),
                                                        downloadButton("downloadData", "Download", 
                                                                       style = "background-color:#3c8dbc; color:white" ),
                                                        br(),
                                                        br(),
                                                        dataTableOutput("cancerdata")
                                                    )
                                           )
                               )
                        )
                    )
            ),
            
            #Application Tab - Unsupervised Learning
            tabItem(tabName = "unsuper",
                    h3("Principle Component Analysis"),
                    br(),
                    withSpinner(dataTableOutput("pcaresult"), type = 1),
                    fluidRow(
                        column(2,
                               box(width=12,title="Select PCs for Biplot",
                                   selectizeInput("pcs1", "First PC",
                                                  choices =  c(paste0("PC", seq(1, dim(PCs$rotation)[2]))),
                                                  selected = "PC1"),
                                   selectizeInput("pcs2", "Second PC",
                                                  choices =  c(paste0("PC", seq(1, dim(PCs$rotation)[2]))),
                                                  selected = "PC2")
                                   # uiOutput("firstPC"),
                                   # uiOutput("secondPC")
                               )
                        ),
                        column(10,
                               box(width = 5, title = "Biplot for Selected PCs",
                                   withSpinner(plotlyOutput("biplot"), type = 7)
                               ),
                               box(width = 7, title = "Variability with Fewer Uncorrelated Variables",
                                   withSpinner(plotOutput("varpcplot"), type = 7)
                               )
                        )
                    )
            ),
            #kNN Clasifier Tab
            tabItem(tabName = "knn",
                    fluidRow(
                        column(3,
                               box(width=12, title=h4("Breast Cancer Classification w/kNN"),
                                   h5(strong("Change k and Dataset to Reflect the Change Below"), style = "color:#3c8dbc"),
                                   br(),
                                   sliderInput("k",
                                               "Number of Neighbors: ",
                                               min = 1,
                                               max = 30,
                                               value = 5),
                                   checkboxGroupInput("checkGroup", label = h4("Dataset Features: "), 
                                                      choices = colnames(select(data, -Class)) , inline = F,
                                                      selected = colnames(select(data, -Class)))
                               ),
                               box(width = 12, title = h4("Confusion Matrix"),
                                   tableOutput('confusionMatrix'),
                                   verbatimTextOutput("classError"))
                        ),
                        column(9,
                               tabsetPanel(type = "tabs",
                                           tabPanel("Classifier", icon = icon("chevron-circle-right"),
                                                    fluidRow(
                                                        box(width = 5, title = h4(strong("kNN Model Fit from 10-fold 
                                                        Repeated CV with 3 Repetition")),
                                                            withSpinner(verbatimTextOutput("knnfrom10Fold"), type = 5)),
                                                        box(width = 7, 
                                                            title = strong("kNN Accuracy Plot - From 10 fold Repeated CV"),
                                                            withSpinner(plotlyOutput("knnAccuracy"), type = 7),
                                                            hr(),
                                                            h4("Confusion Matrix from 10-fold Repeated CV"),
                                                            tableOutput('confusionMatrix10'),
                                                            verbatimTextOutput("classError10"))
                                                    )
                                           ),
                                           tabPanel("Compare Accuracy", icon = icon("gitter"),
                                                    # fluidRow(                           
                                                    #     column(4, selectInput("featureDisplay_x", 
                                                    #                           label = h4("X-Axis"), 
                                                    #                           choices = colnames(select(data, -Class)),
                                                    #                           selected = colnames(select(data, -Class))[1])),
                                                    #     column(4, selectInput("featureDisplay_y", 
                                                    #                           label = h4("Y-Axis"), 
                                                    #                           choices = colnames(select(data, -Class)),
                                                    #                           selected = colnames(select(data, -Class))[2]))
                                                    #     
                                                    # ),
                                                    fluidRow(
                                                        column(12, 
                                                               h4("Left graph shows the Accuracy plot using all predictors 
                                                                  using 10-fold repeated cross validation with 3 repetition 
                                                                  and right graph shows the same but with selected parameters.
                                                                  Please click the button below to re-draw the 
                                                                  right plot. Everytime the selected parameter changes we
                                                                  need to Click to re-draw the plot.", style = "color:#3c8dbc"),
                                                               br(),
                                                               h4(strong("It's a 10-fold Repeated CV so it might take some time 
                                                                         to draw after the Click!"), style = "color:#3c8dbc"),
                                                               hr()
                                                        )
                                                    ),
                                                    fluidRow(
                                                        column(6, 
                                                               br(),
                                                               br(),
                                                               h4(strong("Accuracy Plot with 10-fold Repeated CV 
                                                                         Using ALL Predictors")),
                                                               plotlyOutput("knnAccPlot1")
                                                        ),
                                                        column(6,
                                                               actionButton("plotAcc", icon = icon("arrow-alt-circle-right"),
                                                                            " Plot Accuracy Using Selected Predictors",
                                                                            style = "background-color:#3c8dbc; color: white"),
                                                               h4(strong("Accuracy Plot with 
                                                                         10-fold Repeated CV Using SELECTED Predictors")),
                                                               withSpinner(plotlyOutput("knnAccPlot2"), type = 7)
                                                        )
                                                    )
                                                    # fluidRow(
                                                    #     column(12,
                                                    #            plotlyOutput("scatterPlotAB")
                                                    #     )
                                                    # )
                                           )
                               )
                        )
                    )
            ),
            #Logistic Regression Tab
            tabItem(tabName = "logistic",
                    fluidRow(
                        column(3,
                               box(width=12, title=h4(strong("Breast Cancer Classification w/Logistic Regression")),
                                   h5(strong("Change Dataset to Reflect the Change"), style = "color:#3c8dbc"),
                                   br(),
                                   checkboxGroupInput("checkGroupLog", label = h4("Dataset Features: "), 
                                                      choices = colnames(select(data, -Class)) , inline = F,
                                                      selected = colnames(select(data, -Class)))
                               ),
                               box(width = 12, title = h4(strong("Select Threshold")),
                                   selectizeInput("threshold", "Threshold", selected = 0.5, 
                                                  choices = c(0.3,0.4,0.5,0.6,0.7,0.8))
                                )
                        ),
                        column(9,
                               box(width=12, title = h4(strong("Regression Model Output")),
                                   verbatimTextOutput("logRegOut"),
                                   hr(),
                                   h4(strong("Confusion Matrix from Logistic Regression")),
                                   br(),
                                   tableOutput('confusionMatrixLog'),
                                   verbatimTextOutput("classErrorLog"))
                        )
                    )
                    
            )
        )
    )
)