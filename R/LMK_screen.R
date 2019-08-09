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
LMK_screen <- function(A, ref=NULL, gpa = TRUE ){
  open3d()
  
  p <- dim(A)[[1]]
  k <- dim(A)[[2]]
  n <- dim(A)[[3]]
  
  
  alist <- dimnames(A)[[3]]
  
  if (is.null(alist)){alist <- paste("Ind",1:n, sep=".")}
  
  flist <- matrix(NULL)
  
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
  t <- t
  for (i in 1:n){
    i <- t
    plot3d(ref, col = "red", axes = F, xlim = r, ylim = r, zlim = r, size = 4, xlab = "X", ylab = "Y", zlab = "Z")
    points3d(A[,,i], col = "grey", size = 5)
    text3d(colMeans(ref), pos = 3, texts = paste(alist[i], "compared to"))
    text3d(colMeans(ref), pos = 1, texts = names(ref), col = "red")
    for (j in 1:p){
      lines3d(rbind(ref[j,], A[j,,i]))
    }
    ans <- readline(menu(c("Flag", "Don't Flag", "Previous (No Action)", "Next (No Action)"), title = "Flag for additional inspection? (pres Esc to exit)"))
    if (ans == 1){
      flist[c,1] <- alist[i]
      flist[c,2] <- i
      c <- c+1
      t <- t+1
      clear3d()
    } else if (ans == 3){
      t <- t-1
      clear3d()
    } else {
      t<- t+1
      clear3d()
    }
    
    
  }
  return(flist)
}