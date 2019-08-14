

#' Compare Two landmark arrays
#' 
#' Plot 2 landmark configurations in 3D space.
#' 
#' Presumably, this will be procrustes aligned coordinates, but any 2 p x k x n arrays will work
#' 
#' @param A1 A 3D array of landmark data; A reference 
#' @param A2 A second 3D array of landmark data to compare to the reference
#' @param new Should a new graphing window be opened (Default is to overwrite the previous)
#' 
#' 
#' @export
#' 
LMK_compare_two <- function(A1, A2, new = FALSE){
  #if (all(dim(A1)) != dim(A2)){stop("Landmark configurations must be the same")} ##No idea why this won't run "In addition: Warning message: In if (all(dim(A1)) != dim(A2)) { : the condition has length > 1 and only the first element will be used"
  if (new == TRUE){open3d()}
  clear3d()
  
  r <- range(range(A1),range(A2))
  p <- dim(A1)[[1]]
  k <- dim(A1)[[2]]
  if (length(dim(A1))==3){
    n <- dim(A1)[[3]]
    
    for (i in 1:n){
      plot3d(A1, col = "red", axes = F, xlim = r, ylim = r, zlim = r, size = 4, xlab = "X", ylab = "Y", zlab = "Z")
      points3d(A2[,,i], col = "grey", size = 5)
      text3d(colMeans(A1), pos = 3, texts = paste(names(A2), "compared to"))
      text3d(colMeans(A1), pos = 1, texts = names(A1), col = "red")
      for (j in 1:p){
        lines3d(rbind(A1[j,], A2[j,,i]))
      } 
    }
    }
  else {
    plot3d(A1, col = "red", axes = F, xlim = r, ylim = r, zlim = r, size = 4, xlab = "X", ylab = "Y", zlab = "Z")
    points3d(A2, col = "grey", size = 5)
    text3d(colMeans(A1), pos = 3, texts = paste(deparse(substitute(A2)), "compared to"))
    text3d(colMeans(A1), pos = 1, texts = deparse(substitute(A1)), col = "red")
    for (j in 1:p){lines3d(rbind(A1[j,], A2[j,]))
    }
    }
    
  }



#'Screen Data
#'
#'Interactive function to compare individuals in an array to grandmean (or other refrence) to spot misplaced landmarks/other outliers.
#'
#'Visually inspect a 3D landmark array for outliers/other data inconsistencies. Create a list of potential individuals to inspect further. 
#'
#'@param A a 3D array of landmark data in the form p x k x n
#'@param ref a known standard to compare to (if null (default), function will calculate group mean)
#'@param gpa logical, should a GPA be implemented first. If TRUE, data are assumed to be raw landmark coordinates. If GPA = FALSE, data are assumed to be already aligned procrustes coordinates
#'
#'@export
#'@return Returns a vector of individuals
#'
#'

###Make seperate function, plot.comp, to compare 2 individuals. Use LMK_screen as a wrapper for LMK_compare -- this might be useful for ii?

LMK_screen <- function(A, ref=NULL, gpa = TRUE ){
  
  
  p <- dim(A)[[1]]
  k <- dim(A)[[2]]
  n <- dim(A)[[3]]
  
  
  alist <- dimnames(A)[[3]]
  
  if (is.null(alist)){alist <- paste("Ind",1:n, sep=".")}
  
  flist <- NULL
  
  if (is.null(ref)){
    ref <- mshape(A)
    names(ref) <- "GrandMean"
  }
  
  r <- range(A, na.rm = T)
  
  if (gpa ==TRUE){
    ##Set up gpagen 
  }
  
  ##otherwise, assume procrustes coordinates
  c <- 1
  t <- 1
 
  
  for (i in 1:n){
    LMK_compare_two(ref, A[,,i])
    
    
    ans <- readline(menu(c("Flag", "Don't Flag", "Previous (No Action)", "Next (No Action)", "Esc"), title = "Flag for additional inspection? (pres Esc to exit)"))
    if (ans == 1){
      flist <- rbind(flist, c(alist[i], i)) ###This doesn't seem to be working
      c <- c+1
      t <- t+1
      clear3d()
      
    } else if (ans == 3){
      t <- t-1
      clear3d()
      
    } else if (ans == 5){
     
      return(flist)
    }else{
      t<- t+1
      clear3d()
    }
    
    
  }
 
  return(flist)
}