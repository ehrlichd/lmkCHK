
library(shiny)

# Define UI for application 
shinyUI(fluidPage(
    
    titlePanel("lmkCHK Data Screening"),
    
    sidebarLayout(

    #### Sidebar 
        sidebarPanel(
            
            tabsetPanel(
                tabPanel("Data",
                         
                         ## Select DTA from file
                         fluidRow(fileInput("dta", "Select Data",accept = c(".nts", ".NTS", ".ntsys", ".NTSYS", ".dta", ".DTA", ".bsc", ".BSC"))),
                         
                         fluidRow(fileInput("lmkpairs","Define paired lmks", accept = c(".csv", ".CSV", ".txt", ".TXT"))),
                         
                         fluidRow(fileInput("lmknames", "Specify lmk names",accept = c(".csv", ".CSV", ".txt", ".TXT")))
                         ),
                
                
                tabPanel("Align",
                         fluidRow(radioButtons("gpa", "Define Alignment Criteria", choices = c("Align (no scale)" = "nogpa", "GPA (no sliding)" = "gpa", "GPA (with sliding)" = "sgpa"))),
                         fluidRow(helpText("Define semi landmarks to be treated as contours or surfaces")),
                         fluidRow(fileInput("curves", "Define Curves")),
                         fluidRow(fileInput("surf", "Define Surfaces")),
                         
                ),
                
                tabPanel("Toggle LMKs",
                        fluidRow(uiOutput("lmkcontrol"))
                ),
                
                tabPanel("Run",
                         fluidRow(actionButton("goButton", "Run!"))
                ),
                
                type = "tabs")
            ),
            
        #### Main panel
        
        mainPanel(
            fluidPage(
                titlePanel(textOutput("deets")),
                tabsetPanel(
                    tabPanel("active LMK",
                             #h2("PLACEHOLDER") 
                             #plotOutput("plot2") ## call 3D plot
                             rglwidgetOutput("plot2")),
                    tabPanel("PCA",
                             plotOutput("plot1")),
                    tabPanel("Var",
                             plotOutput("variance")),
                    tabPanel("Outliers",
                             plotOutput("outliers"))
                    )
                ),
            

            
            ) ## close mainPanel
        ) ## close sidebarLayout
    
    ) ## close whole app fluidPage
    
    ) ## Close ui


    
    



