

#' Compare Two landmark arrays
#' 
#' Plot 2 landmark configurations in 3D space.
#' 
#' Presumably, this will be procrustes aligned coordinates, but any 2 p x k x n arrays will work
#' 
#' @param A1 A 3D array of landmark data; A reference 
#' @param A2 A second 3D array of landmark data to compare to the reference
#' @param new Should a new graphing window be opened (default is to overwrite the previous)
#' @param main A character string to display as the title. If main= NULL (default) function will attempt to deparse object names.
#' 
#' 
#' @export
#' 
LMK_compare_two <- function(A1, A2, new = FALSE, main = NULL){
  #if (all(dim(A1)) != dim(A2)){stop("Landmark configurations must be the same")} ##No idea why this won't run "In addition: Warning message: In if (all(dim(A1)) != dim(A2)) { : the condition has length > 1 and only the first element will be used"
  if (new == TRUE){open3d()}
  
  
  ###Instead of clearing 3D, find a better way to control plots/outputs. Simple way to plot all shape changes?
  
  r <- range(range(A1),range(A2))
  p <- dim(A1)[[1]]
  k <- dim(A1)[[2]]
  if (length(dim(A1))==3){ ## to compare groups
    n <- dim(A1)[[3]]
    
    for (i in 1:n){
      plot3d(A1, col = "red", axes = F, xlim = r, ylim = r, zlim = r, size = 4, xlab = "X", ylab = "Y", zlab = "Z")
      points3d(A2[,,i], col = "grey", size = 5)
      
      if (is.null(main)){
      text3d(colMeans(A1), pos = 3, texts = paste(names(A2), "compared to"))
      text3d(colMeans(A1), pos = 1, texts = names(A1), col = "red")
      } else {
        text3d(colMeans(A1), texts = as.character(main))
      }
      for (j in 1:p){
        lines3d(rbind(A1[j,], A2[j,,i]))
      } 
    }
    
    } else {
      plot3d(A1, col = "red", axes = F, xlim = r, ylim = r, zlim = r, size = 4, xlab = "X", ylab = "Y", zlab = "Z")
      points3d(A2, col = "grey", size = 5)
      text3d(colMeans(A1), pos = 3, texts = paste(deparse(substitute(A2)), "compared to"))
      text3d(colMeans(A1), pos = 1, texts = deparse(substitute(A1)), col = "red")
      for (j in 1:p){
        lines3d(rbind(A1[j,], A2[j,]))
      }
    }
  }



#' Compare Wireframes
#'
#'Function to plot two landmark configurations as wireframes, and compare directly
#' 
#' 
#'Plots two LMK configurations on top of each other as some combination of fixed points and wireframe. Presumably this will be PC min/max, two group means, or mshape and a target individual, but any two matching configurations should work
#' 
#' @param A1 a reference
#' @param A2 a target
#' @param links a wireframe listing landmarks to be linked
#' @param cols a list of two colors to use for plots (default Red, Blue)
#' @param pts Should landmarks be plotted as points, spheres, or not at all.
#' @param rad radius to plot for spheres, default is .001. Unclear how project specific
#' @param vectors Should vectors be drawn linking homologous landmarks. This may help show the direction of change, or could look clunky
#' @param ... Additional arguments to pass to plotting functions.
#' 
#' 
#' @export
#' 
#' 

LMK_wireframe <- function(A1, A2, cols = c("red","blue"), links = NULL, pts = c("point", "sphere", "none"), rad = c(.001, .001), vectors = FALSE, ...){
  
  r <- range(A1,A2)
  p <- dim(A1)[[1]]
  k <- dim(A1)[[2]]
  l <- dim(links)[[1]]
  
  
  plot3d(rbind(A1,A2), type = "n", xlim = r, ylim = r, zlim = r, axes = F)
  if(pts == "point"){
    points3d(A1, col = cols[1])
    points3d(A2, col = cols[2])
  } else if (pts =="sphere"){
    spheres3d(A1, col = cols[1], radius = rad)
    spheres3d(A2, col = cols[2], radius = rad)
  }
  if(!is.null(links)){
    
    for (i in 1:l){
      segments3d(A1[links[i,],], col = cols[1],...)
      segments3d(A2[links[i,],], col = cols[2],...)

    }
  }
  if(vectors == TRUE){
    
  }
}
