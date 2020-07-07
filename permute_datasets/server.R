#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

# Define server logic required to draw a histogram
shinyServer(function(input, output) {
    
    library(lmkCHK)
    
    

    output$plot1 <- renderPlot({
        
        if (input$goButton == 0){
            return()
        } else {
            isolate({
                dat <- input$file.dta
                
                names <- input$file.lmkname
                
                gpa <- gpagen(dat)
                pca <- prcomp(gpa$coords)
            })
            plot(pca$x)
        }
        
        
    })

})
