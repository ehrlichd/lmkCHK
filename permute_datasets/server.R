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
    
    options(shiny.maxRequestSize=30*1024^2) ## change max file size to 30 MBs
    
    library(lmkCHK)
    
    v <- reactiveValues(dat = NULL)
    observeEvent(input$goButton, {
        v$dat <- readland.nts(input$file$datapath)
        v$gpa <- gpagen(v$dat)
        v$pca <- prcomp(two.d.array(v$gpa$coords))
    })

    
    output$plot1 <- renderPlot({
        if (is.null(v$pca)) return()
        plot(v$pca$x)
    })
    
    output$deets <- renderText( {
        if (is.null(v$dat)) return()
        dd <- dim(v$dat)
        
        paste("This dataset contains", dd[1], "landmarks for", dd[3], "individuals in", dd[2], "dimensions")
        })
    
    output$plot2 <- renderPlot({
        if (is.null(v$dat)) return()
        clear3d()
        points3d(v$gpa$consensus)
    })
}
)
