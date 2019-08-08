#'@name lmkCHK
#'
#'
#'@title Landmark Check
#'@author Daniel E. Ehrlich
#'
#'@description Helper functions for visualizing and screening 2D/3D landmark data.
#'
#'@import rgl
#'@import geomorph
#'@import Morpho
#'@import stats
#'@import utils
#'@import graphics
#'@import psych
#'@import lme4
#'@import grDevices
#'
#'
NULL

#####Need to find sample LMK data to include; MWU skulls?#####

#'Write DTA/NTS
#'
#'Write 3D array as a .dta/.nts file
#'
#'@param A a 3D array in the form of [p x k x n]
#'@param filepath Filename/Filepath, including extension: .dta, or .nts
#'
#'
#'@export
#'
#'@author Daniel Ehrlich
#'
#'
#'
LMK_writeland_nts <- function(A, filepath){
  p <- dim(x)[1]
  k <- dim(x)[2]
  n <- dim(x)[3]
  lbls <- dimnames(x)[[3]]
  
  #check/format missing data
  if (anyNA(x)){
    x[is.na(x)] = 9999
    h = paste(1, paste(n,"L", sep=""), (p*k), 1, 9999, "Dim=3", sep = " ")} else {
      h = paste(1, paste(n,"L", sep=""), (p*k), 0, "Dim=3", sep = " ")
    }
  #write header line
  write(c(h, "\n"), file=filepath, ncolumns = 1, append= FALSE)
  
  #write specimen labels
  write(c(lbls, "\n"), file=filepath, ncolumns = 1, append= TRUE , sep = "\t")
  
  #write coords
  for(i in 1:n){
    write(c(t(x[,,i]), "\n"), file=filepath, append=T, sep = "\t", ncolumns = 3)
    
  }
}

#####


#' Swap Landmarks
#' 
#' 
#' Function to swap landmark configurations, or set as missing
#' 
#' @param a a p x k x n matrix of landmark coordinates
#' @param l1 a numeric vector of landmark index (sequence number) to be changed
#' @param l2 a numeric vector of landmark index (Sequence number) to change to
#' 
#' @author Daniel Ehrlich
#' 
LMK_swap <- function(a, l1, l2){
  #function to swap landmark configurations either by list or single landmarks
  #change the order, or set LMKs to NA
  #code written by DEE
  #WIP
  t.a <- a
  if (length(l1) != length(l2)){stop("Vectors must be the same length")}
  
  for (i in 1:length(l1)){
    
    t.a[l1[i],,] <- a[l2[i],,]
    t.a[l2[i],,] <- a[l1[i],,]
    
  }
  }

#####
#' Get sliders
#' 
#' Define a sliders table (geomorph) for a contour with fixed enpoints
#' 
#' Assumes all points slide sequentially between endpoints
#'  
#'@param c1 A vector of landmark numbers
#' 
#'@return Returns a table of sliders suitable for use with geomorph::gpagen(). This function assumes all points slide sequentially between the max/min landmarks in the sequence.
#' 
#' @export
#' 
#' @author Daniel Ehrlich
#'    
LMK_get_sliders <- function(c1){
  
  t <- as.numeric(as.matrix(c1)) # convert list to vector
  l <- length(t)
  l2 <- l-1
  
  sl <- matrix(data = cbind(t[-c(l2,l)], 
                            t[-c(1,l)], 
                            t[-c(1,2)]), 
              nrow = l-2, 
              ncol = 3, 
              dimnames = list("row" = paste("lmk", seq(2, l-1), sep = "." ), "col" = c("before", "slide", "after")))
  
  return(sl)
  
}

#####

#' Euclidean Distance
#' 
#' Calculate Euclidean distance between two sets of 3D landmark data. Data should be in the typical p x k x n array.
#' 
#' 
#' @param A1 A 3D landmark array in the form p x k x n
#' @param A2 A second 3D landmark array in the form p x k x n
#' 
#'  
#' @return Returns a list (3L) containing the direct Euclidean distance for all individuals/landmarks, as well as sumarized by individual and landmark
#' 
#' @export
#' 
LMK_euD <- function(A1, A2){
  if ( length(dim(A1))==3 & length(dim(A2))==3){
    if (dim(A1) != dim(A2)){stop("Arrays must have the same extent")}
  
  
  l <- dim(A1)[[1]]
  n <- dim(A1)[[3]]
  
  o1 <- matrix(data = NA,nrow = l, ncol = n, dimnames = list(dimnames(A1)[[1]], dimnames(A1)[[2]]))
  
  
  for (c in 1:n){
    #for each individual
    for (r in 1:l){
      #for each lmk
      
      #calculate distance
      o1[r,c] <- sqrt(sum(((A1[r,1,c]-A2[r,1,c])^2), ((A1[r,2,c]-A2[r,2,c])^2), ((A1[r,3,c]-A2[r,3,c])^2)))
    }
    
  }
  
  all.dif <- o1
  lmk.dif <- data.frame("avg.dif" = rowMeans(o1))
  ind.dif <- data.frame("avg.dif" = colMeans(o1))
  
  out <- list("all.dif"= all.dif, "by.lmk" = lmk.dif, "by.ind" = ind.dif)
  return(out)
  } ###figure out how to generalize input and internal calcs
}


#####
#'Color Ramp
#'
#'Create a color ramp based on a grouping factor
#'
#'Creates color ramp for a given factor and returns a vector of collor assignments for each individual using grDevices::rainbow(). Useful for plotting by group.
#'
#'@param grp A grouping variable, should be a factor
#'@param mute Logical, should color be muted. If true (default), color values will be generated with saturation value of .4, if false, values will have full saturation (1.0)
#'
#'@return Returns a vector of color assignments for each individual based on a grouping factor. 
#'@export
#'
#'@author Daniel Ehrlich
#'


LMK_colramp <- function(grp, mute = TRUE){
  
  f.grp <- as.factor(grp)
  #ensure vector is a factor to match index
  
  c.grp <- as.character(grp)
  #generate character list to recieve colors
  
  l <- length(levels(grp))
  if( mute == TRUE){
    cr <- grDevices::rainbow(length(levels(f.grp)), s = 0.4, v = 1)
    #generate color ramp based on levels of data
    
  }else{
    cr <- grDevices::rainbow(length(levels(f.grp)))
    #generate color ramp based on levels of data
    
  }
  
  for (i in 1:l){
    c.grp[as.integer(f.grp) == i] <- cr[i]
  }
  return(c.grp)
}

#####

#'Backscale an array
#'
#'Scale each individual of a landmark array by an arbitrary value (e.g. centroid size)
#'
#'@param A A 3D landmark array
#'@param cs a vector of centroid sizes
#'
#'@return Returns a scaled 3D landmark array
#'
#'@export

LMK_bscale <- function(A, cs){
  
  al <- dim(A)[[3]]
  cl <- length(cs)
  if (length(al) != length(cl)){
    stop("data mismatch")
  } else {
    for (i in 1:al){
      A[,,i] <- A[,,i]*cs[i]
    }
  }  
  return(A)
}

#####

#'Group Apply
#'
#'Apply a function (f) to a set of data (d) by a grouping variable (grp).
#'
#'@param d a dataset
#'@param f a function
#'@param grp a grouping varaible
LMK_grp_apply <- function(d, f, grp){
  
  out <- as.data.frame(matrix(nrow = (length(levels(grp))+1), ncol = length(d)*2))
  c.names <- c(colnames(d), paste(colnames(d),"n=", sep = "."))
  id <- c(1:length(d), (1:length(d)+0.5))
  colnames(out) <- c.names[order(id)]
  row.names(out) <- c(levels(grp), paste("G", deparse(substitute(f)), sep = "."))
  
  for (i in 1:length(levels(grp))){
    #for each group
    
    for (j in 1:length(d)){
      #for each column
      
      out[i,((j*2)-1)] = f(na.omit(d[as.integer(grp)==i,j]))
      #perform the operation f, dropping NAs
      
      out[i,(j*2)] = length(na.omit(d[as.integer(grp)==i,j]))
      #count the individuals used in the operation
    }
    
  }
  #calculate f for whole data (eg Grand Mean, Grand Median, etc)
  for (k in seq(from = 1, to = length(out), by = 2)){
    out[(length(levels(grp))+1),k] = f(out[1:length(levels(grp)),k])
    out[(length(levels(grp))+1),(k+1)] = sum(out[1:length(levels(grp)),(k+1)])
  }
  
  return(out)
}

#' Intra class correlations and Cronbach's Alpha
#' 
#' Compute ICC and CA for each dimension (X,Y,Z) of arrays
#' 
#' @param obs1 a 3D landmark array
#' @param obs2 a 3D landmark array
#' 
#' @return Returns a table of ICC coefficients, associated p-values, and Chronbach's Alpha
#' 
#' @export
#' 
#' 
LMK_iccA <- function(obs1, obs2){
  #inter- or intra- oberserer error
  #where obs1 and obs2 are two identical landmark data sets
  #function to perform inter class correlation and reports F- and p- values, chronbach's alpha, as well as euclidean distance for a set of landmarks
  #currently assumes: .dta format, 3 dimensions, 34 lmks, object symmetry
  #just get it from psych package
 
  l <- dim(obs1)[[1]]
  d <- dim(obs1)[[2]]
  out <- data.frame("ICC.X" = numeric(l), "f.X" = numeric(l),"sig.X" = numeric(l), "ICC.Y" = numeric(l),"f.Y" = numeric(l), "sig.Y" = numeric(l), "ICC.Z" = numeric(l), "f.Z"=numeric(l),"sig.Z" = numeric(l),"CrA.X" = numeric(l), "CrA.Y" = numeric(l), "CrA.Z" = numeric(l),row.names = dimnames(obs1)[[1]])
  
  for (r in 1:l){
    tx <- psych::ICC(as.matrix(cbind(obs1[r,1,], obs2[r,1,])), lmer=FALSE)
    ty <- psych::ICC(as.matrix(cbind(obs1[r,2,], obs2[r,2,])), lmer=FALSE)
    tz <- psych::ICC(as.matrix(cbind(obs1[r,3,], obs2[r,3,])), lmer=FALSE)
    ax <- psych::alpha(as.matrix(cbind(obs1[r,1,], obs2[r,1,])))
    ay <- psych::alpha(as.matrix(cbind(obs1[r,2,], obs2[r,2,])))
    az <- psych::alpha(as.matrix(cbind(obs1[r,3,], obs2[r,3,])))
    
    out[r,1] <- tx$results$ICC[3]
    out[r,2] <- tx$results$F[3]
    out[r,3] <- tx$results$p[3]
    
    
    out[r,4] <- ty$results$ICC[3]
    out[r,5] <- ty$results$F[3]
    out[r,6] <- ty$results$p[3]
    
    out[r,7] <- tz$results$ICC[3]
    out[r,8] <- tz$results$F[3]
    out[r,9] <- tz$results$p[3]
    
    out[r,10] <- ax$total[1]
    out[r,11] <- ay$total[1]
    out[r,12] <- az$total[1]
    
  }
  
  return(out)
}


#####

#'Shapespace PCA
#'
#'Plot 2D or 3D scatter plot of tangent shape-space for a GPA object (of gpagen())
#'
#'@param A A 3d landmark array
#'@param xPC an integer indicating which PC to plot, default = 1
#'@param yPC an integer indicating which PC to plot, default = 2
#'@param zPC an integer indiciating which PC to plot, default = NULL (2D)
#'@param grp a grouping factor
#'
#'@export
#'
#'
LMK_PCA_plot <- function(A, xPC = 1, yPC = 2, zPC = NULL, grp){
  f.PCA <- geomorph::plotTangentSpace(A, warpgrids = F)
  grp <- as.factor(grp)
  l <- length(levels(grp))
  
  ##Calculate means
  Mean.tab <- matrix(nrow = l+1, ncol = 10)
  for (i in 1:l){
    Mean.tab[i,] <- colMeans(f.PCA$pc.scores[as.integer(grp) == i, 1:10 ])
  }
  Mean.tab[l+1,] <- colMeans(f.PCA$pc.scores[,1:10]) 
  rownames(Mean.tab) <- c(levels(grp),"GrandM")
  colnames(Mean.tab) <- paste("PC",1:10, sep=".")
  
  if (is.null(zPC)){
    plot(f.PCA$pc.scores[,c(xPC, yPC)], col = LMK_colramp(grp), pch = 16)
    points(Mean.tab[,c(xPC, yPC)], pch = 21, bg = c(rainbow(l, s = .4), "grey"), cex = 1.5)
    legend("topleft", legend = c(levels(grp), "G.mean"), pch = 16, col = c(rainbow(l, s=.4), "grey"))
  } #implement 3D plot
  
}

