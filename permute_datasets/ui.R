#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

# Define UI for application that draws a histogram
shinyUI(fluidPage(
    
    titlePanel("Basic widgets"),
    
    sidebarLayout(
        
    #### Sidebar 
        sidebarPanel(
            
            ## First row
            fluidRow(column(3, checkboxGroupInput("checkGroup",h3("Checkbox group"), choices = list("Choice 1" = 1, "Choice 2" = 2, "Choice 3" = 3), selected = 1))),
            
            ## Select dta
            fluidRow(column(1, fileInput("file.dta", h3("File input")))),
            ## select lmk names
            fluidRow(column(2, fileInput("file.lmkname", h3("File input")))),
            
            
            ## GPA Y/N 
            fluidRow(column(1,radioButtons("gpa", h3("Radio buttons"),choices = list("Choice 1" = 1, "Choice 2" = 2,"Choice 3" = 3),selected = 1))),
            
            ## lmks
            fluidRow(column(2,selectInput("select", h3("Select box"), choices = list("Choice 1" = 1, "Choice 2" = 2, "Choice 3" = 3), selected = 1))
                     ),
            
            fluidRow(actionButton("goButton", "Go!"))
        ),
        
        #### Main panel
        
        mainPanel(  
            # Output: Plot of the requested variable against mpg
            plotOutput("plot1")
            
            )
    )
    
    ))

