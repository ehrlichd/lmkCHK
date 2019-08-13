#'Plot Outliers
#'
#'Function to identify outliers visually
#'
#'
#'@param A Raw 3D coordinate data, procrustes aligned coordinate data, or an object of class == gpagen
#'@param gpa Logical, should a GPA be calculated. If TRUE, a is assumed to be raw coordinate data and GPA is implemented via geomorph::gpagen. If FALSE, no imposition will be done. If TRUE, gpa-aligned coordinates and mshape will also be exported
#'@param plotALL Logical, plot all specimens (and mean) in 3D space
#'
#'@param ... Additional parameters to pass to gpagen
#'
#'@export
#' 
#'@return Returns a list containing summary information of the dataset, individual procrustes distances, GPA aligned Procrustes Coordinates and mean shape configuration

LMK_plotoutliers <- function(A, gpa = TRUE, plotALL = TRUE, ...){
  name <- deparse(substitute(A))
  p <- dim(A)[[1]]  #p landmarks
  k <- dim(A)[[2]]  #k dimensions
  n <- dim(A)[[3]]  #n observations
  
  if (is.null(dimnames(A)[[3]])){
    lbl <- paste("Ind",1:n, sep = ".")
  } else { lbl <- dimnames(A)[[3]] }#Get IDs 
  
  
  if (class(A)=="gpagen"){A <- A$coords}
  if (gpa == TRUE) { A <- gpagen(A, ...)$coords}
  
  ##Create tables and variables for later
  dist <- matrix(NA, nrow = n, ncol = 3)
  colnames(dist) <- c("ind", "proc.D", "index")
  dist[,1] <- lbl
  
  
  tab <- NULL
  
  ##Grand mean
  grandM <- mshape(A)
  
  ##Range of data
  r <- range(A)
  
  
  ####Plot procrustes aligned GrandMean shape and all observations
  
  if (plotALL == TRUE){
    open3d()
    plot3d(grandM, xlim = r, ylim = r, zlim = r, axes = F, type = "s", col = "red", size = 1)
    text3d(colMeans(grandM), texts = "Outliers of", adj = c(.5,0))
    text3d(colMeans(grandM), texts = name, adj = c(.5,1))
    
    for (i in 1:n){
      points3d(A[,,i], col = hsv(0,0,.8, alpha = .8))
      }
  
    }
  
  
  ####Calculate average procrustes distance ("error") from GrandMean
  
  for(i in 1:n){
    for(j in 1:p){
      tab[j] <- sqrt(sum(((grandM[j,1] - A[j,1,i])^2), ((grandM[j,2] - A[j,2,i])^2), ((grandM[j,3] - A[j,3,i])^2))) ##Should use LMK_euD
      if (j == p){
        dist[i,2] <- sum(tab)
        dist[i,3] <- i}
    }
    
  }

#  ###Alternative calc using euD  
#  for (i in 1:n){
#    dist[i,2] <- sum(LMK_euD(grandM, A[,,i]))
#    dist[i,3] <- i
#  }
  
  ##Re-order table
  dist <- dist[order(dist[,2], decreasing = T),]

  y <- range(dist[,2])
  x <- c(1,n)
  
  ####Plot 
  plot(1:n, dist[,2], type = "n", main = paste("Outliers plot of",name), xlab = "Rank", ylab = "Proc Distastance")
  
  for (i in 1:n){
    
  if (dist[i,2] > mean(as.numeric(dist[,2])) + 3*sd(as.numeric(dist[,2]))){
    points(i, dist[i,2], pch = 16, col = "red")
    #text(i, dist[i,2],labels = dist[i,1], col = "red", adj = c(1,1)) ##No idea why this doesn't work
    
    } else {
      points(i, dist[i,2], pch = 16, col = "grey")
    }
    
  }
  sum <- matrix(
    c("p", "k", "n", "mean.ProcD", "sd.ProcD", "min.ProcD", "max.ProcD", 
      p, k , n, 
      mean(as.numeric(dist[,2]), na.rm = T), 
      sd(as.numeric(dist[,2]), na.rm = T), 
      range(as.numeric(dist[,2]), na.rm = T)), nrow = 7, ncol = 2)
  
  out <- list("summary.info" = sum, "ind.info" = dist, "proc.coords" = A, "mshape" = grandM)
  
  return(out)
}


###Outliers identified from this plot in the SOF data set do not seem to be outliers?