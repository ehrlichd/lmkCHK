#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

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
                
                
                tabPanel("Analyze",
                        fluidRow(uiOutput("lmkcontrol"))
                ),
                
                tabPanel("Run",
                         fluidRow(actionButton("goButton", "Run!"))
                ),
                
                type = "tabs"
                
                
            )
            ),
            
        #### Main panel
        
        mainPanel(  
            # Output: Plot of the requested variable against mpg
            h2(textOutput("deets")),
            plotOutput("plot1"),
            plotOutput("plot2")
        )
            
    
    )
    )
    )


