#'Plot Outliers
#'
#'Function to identify outliers visually
#'
#'
#'@param a Raw 3D coordinate data, procrustes aligned coordinate data, or an object of class == gpagen
#'@param gpa Logical, should a GPA be calculated. If TRUE, a is assumed to be raw coordinate data and GPA is implemented via geomorph::gpagen. If FALSE, no imposition will be done.
#'@param plotALL Logical, plot all specimens (and mean) in 3D space
#'
#'@param ... Additional parameters to pass to gpagen
#'
#'@export
#' 
#'@return Returns a list of specimens and their squared procrustes distance to the mean shape

LMK_plotoutliers <- function(a, gpa = TRUE, plotALL = TRUE, ...){
  
  p <- dim(a)[[1]]  #p landmarks
  k <- dim(a)[[2]]  #k dimensions
  n <- dim(a)[[3]]  #n observations
  
  if (is.null(dimnames(a)[[3]])){
    lbl <- paste("Ind",1:n, sep = ".")
  } else { lbl <- dimnames(a)[[3]] }#Get IDs 
  
  
  if (class(a)=="gpagen"){a <- a$coords}
  if (gpa == TRUE) { a <- gpagen(a, ...)$coords}
  
  ##Create tables and variables for later
  dist <- matrix(NA, nrow = n, ncol = 2, dimnames = list(NULL, c("ind", "proc.D")))
  dist[,1] <- lbl
  
  
  tab <- NULL
  
  ##Grand mean
  grandM <- mshape(a)
  
  ##Range of data
  r <- range(a)
  
  
  ####Plot procrustes aligned GrandMean shape and all observations
  
  if (plotALL == TRUE){
    plot3d(grandM, xlim = r, ylim = r, zlim = r, xlab = "", axes = F, type = "s", col = "red", size = 1)
    
    for (i in 1:n){
      points3d(a[,,i], col = hsv(0,0,.8, alpha = .8))
      }
  
    }
  
  
  ####Calculate average procrustes distance ("error") from GrandMean
  
  for(i in 1:n){
    for(j in 1:p){
      tab[j] = sqrt(sum(((grandM[j,1] - a[j,1,i])^2), ((grandM[j,1] - a[j,1,i])^2), ((grandM[j,1] - a[j,1,i])^2)))
      if (j == p){dist[i,2] = sum(tab)}
    }
    
  }
  
  ##Re-order table
  dist <- dist[order(dist[,2], decreasing = T),]

  y <- range(dist[,2])
  x <- c(1,n)
  
  ####Plot 
  plot(1:n, dist[,2], type = "n")
  
  for (i in 1:n){
    
  if (dist[i,2] > mean(as.numeric(dist[,2])) + 3*sd(as.numeric(dist[,2]))){
    points(i, dist[i,2], pch = 16, col = "red")
    text(i*10, dist[i,2],labels = dist[i,1], col = "red", adj = c(0,0))
    
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
  
  out <- list(sum, dist)
  names(out) <- c("summary.info")
  return(out)
}
