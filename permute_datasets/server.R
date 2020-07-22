#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

# Define server logic 
shinyServer(function(input, output) {
    
    options(shiny.maxRequestSize=30*1024^2) ## change max file size to 30 MBs
    
    library(lmkCHK)
    
    v <- reactiveValues(dat = NULL)
    param <- reactiveValues(crv = NULL, 
                            srf = NULL)
    lmklist <- reactiveValues(all = NULL)
    
    ### Real issue seems to be initializing empty
    
    ### PCA should be calcuated OUTSIDE of observe button; possible as a reactive value
    ### GPA (and loading) need to be controlled by the button. Could this also be saved to a reactive() object???
    
    observeEvent(input$goButton, {
        v$dat <- readland.nts(input$dta$datapath)
        
        
        
        if (input$gpa=="nogpa"){
            v$coords <- ProcGPA(v$dat, scale = F, CSinit = F)$rotated ## Morpho nomenclature
            v$mshape <- mshape(v$coords)
        } 
        
        else {
            if (!is.null(input$curves)) param$crv <- as.matrix(read.csv(input$curves$datapath))
            if (!is.null(input$surf)) param$srf <- as.matrix(read.csv(input$surf$datapath))
            
            tt <- gpagen(v$dat, curves = param$crv, surfaces = param$srf)
            
            v$coords <- tt$coords
            v$mshape <- tt$consensus
        }
        
        })
    


    ##### PCA #####
    
    output$plot1 <- renderPlot({
        if (is.null(v$dat)) return()
        v$pca <- prcomp(two.d.array(v$coords[lmklist$all,,]))
        plot(v$pca$x)
    })
    
    
    ##### Text #####
    output$deets <- renderText( {
        if (is.null(v$dat)) return()
        v$dims <- dim(v$dat)
        
        paste("This dataset contains", sum(!is.na(lmklist$all)), "landmarks for", v$dims[3], "individuals in", v$dims[2], "dimensions")
        
        })
    
    ##### 3D Plot #####
    output$plot2 <- renderPlot({
        if (is.null(v$dat)) return()
        clear3d()
        lmklist$all <- as.numeric(input$lmkall)
        points3d(v$mshape)
        spheres3d(v$mshape[lmklist$all,], col = "red", radius = .2)
            
        
    })
    
    
    #### Variance Plot ####
    
    output$variance <- renderPlot({
        if(is.null(v$dat)) return()
        
        
        
        plotVar <- function(a, num = TRUE,tabonly=FALSE,...){
            p <- dim(a)[[1]]
            k <- dim(a)[[2]]
            n <- dim(a)[[3]]
            
            lmkVar <- matrix(nrow = p, ncol = k)
            for(i in 1:p){
                for(j in 1:k){ 
                    lmkVar[i,j] <- var(a[i,j,])
                }
            }
            means <- apply(lmkVar,1,mean)
            if(tabonly==TRUE){
                return(means)
            }else if (num ==TRUE) {
                plot(means, ylab = "Var", type = "n",...)
                text(means, labels = 1:190, cex = .6, col = "red",...)
                abline(h = mean(means), col = "green")
                abline(h = mean(means)*1.67, col = "dark green")
            }else {
                plot(means, ylab = "Var", type = "n",...)
                text(means, labels = lmknames$lmkName, cex = .4, col = "red",...)
                abline(h = mean(means), col = "green")
                abline(h = mean(means)*1.67, col = "dark green")
            }
            return(means)
        }
        
        plotVar(v$coords[lmklist$all,,])
    })
    
    ##### Error Plot #####
    
    output$outliers <- renderPlot({
        if(is.null(v$dat)) return()

        LMK_plotoutliers(v$coords[lmklist$all,,], gpa = FALSE,plotALL = F)
    })
    
    
    
    
    
    #### Toggle Landmarks #####
    output$lmkcontrol <- renderUI({
        if (is.null(v$dat)) return()
        
        if (is.null(input$lmknames)) {
            lmkchoice <- paste("lmk",1:v$dims[1],sep=".")
        } else {
            lmkchoice <- input$lmknames
        }
        
        checkboxGroupInput("lmkall", "Include (red) / Exclude (black) lmks",
                           #selected = 1:v$dims[[1]],
                           choiceNames = lmkchoice,
                           choiceValues = 1:v$dims[[1]])
    })
}
)
