#'@name lmkCHK
#'
#'
#'@title Landmark Check
#'@author Daniel E. Ehrlich
#'
#'@description Helper functions for visualizing and screening 2D/3D landmark data.
#'
#'@import Morpho
#'@import geomorph
#'@import rgl
#'@import stats
#'@import utils
#'@import graphics
#'@import grDevices
#'@import shiny
#'
#'
#'
NULL



#'Symmetrize
#'
#'Function to provide easy specification and handling of multistep symmetrization process
#'
#'Under full sequence landmarks are first estimated using bilateral reflection. If missing data are still present, landmarks are estimated using TPS interpolation of n=3 nearest neighbors. Now complete LMK arrays are reflected a final time and averaged creating perfectly symmetrical configurations
#'
#'@param A A landmark configuration in the form p x k x n
#'@param LMpair A matrix listing paried landmarks, those not listed are assumed to be midline
#'
#'
#'@export
#'
#'@return Returns a symmetrized landmark configuration with no missing data
#'
LMK_sym <- function (A, LMpair){
  ##Add on toggles for each stage to customize how symmetrization happens
  a.mir <- fixLMmirror(A, LMpair)
  
  if (anyNA(a.mir)){a.mir <- fixLMtps(a.mir)$out}
  
  sym <- symmetrize(a.mir, LMpair)
  
  return(sym)
}


##### Plotting




#' Graphical Color Picker
#' 
#' Interactive function to choose precise color from a pallete of tints or shades
#' 
#' Select any number of colors from palette
#' 
#' @param type One of c("shades, "tints") to plot
#' @param transparency Set color transparency (alpha) from 0 (invisible) to 1 (opaque). Default is opaque
#' @param  n Number of colors to pick
#'
#'  
#' @return Returns a matrix of n colors listed in 3 formats: Rcolor, RGB, and HSV.
#' 
#'
#' @export
#'
#'

LMK_colpick <- function(n=1, type = c("shades", "tints"), transparency = 1){
  
  tr <- transparency 
  out <- matrix(nrow = n, ncol = 8)
  colnames(out) <- c("R.color", "R", "G","B","H","S","V","alpha")
  
  ## plot shades or tints
  
  plot(1:100, type = "n")
  if (type == "tints"){
    for( i in 1:100){
      for(j in 1:100){
        points(i,j,col = hsv(i/100, s=j/100, v = 100/100), pch = 16, cex = 2.5)
      }
    }
  } else {
    for( i in 1:100){
      for(j in 1:100){
        points(i,j,col = hsv(i/100, s=100/100, v = j/100), pch = 16, cex = 2.5)
      } 
    }  
  }
  
  ### Pick the colors
  print(paste("Choose", n, "colors"))
  print("Press [Esc] to exit")
  
  cols <- locator(n=n)
  
  if (type =="tints"){
    h <- cols$x/100
    s <- cols$y/100
    v <- rep(1,n) 
  } else {
    h <- cols$x/100
    s <- rep(1,n)
    v <- cols$y/100
  } 
  
  out[,1] <- hsv(h, s, v, tr) ## Rcol
  out[,2:4] <- col2rgb(out[,1]) ## RBG
  out[,5] <- h
  out[,6] <- s
  out[,7] <- v
  out[,8] <- tr
  
  return(out)
  
}






#' limset
#'
#'
#' Function to easily set plot buffer
#'
#' Takes a vector of data, calculates the range of data and scales it by a factor. Returns new range of length=2
#'
#' @param x A vector to be plotted
#' @param factor to expand the data (default=1.2)
#'
#' 
#' @examples 
#' 
#' vec <- 1:10
#' extend_vec <- LMK_limset(vec)
#' extend_vec
#' 
#' @export




LMK_limset <- function(x, factor = 1.2){
  r <- range(x)
  int <- abs(r[1] - r[2])
  int2 <- int * factor
  
  d <- (int2 - int)/2
  r1 <- r
  low <- r[1]-d
  hi <- r[2]+d
  out <- c(low,hi)
  return(out)
}




#'Confidence Ellipses
#'
#'Draw confidence ellipses around data
#'
#'@param dat data to be plotted in matrix form [X,Y]
#'@param ci Confidence interval to be plotted. Must be one of c(67.5, 90, 95,99)
#'@param linesCol Color for the line. Currently takes hsv() format 
#'@param fillCol color for the fill. Currently only takes hsv() format, set NULL for no fill
#'@param smoothness Smoothness for ellipses. Default should be sufficient but is customizable.
#'
#' @examples 
#' set.seed(1)
#' 
#' datx <- rnorm(100,0,1)
#' daty <- rnorm(100,0,10)
#' 
#' dat <- cbind(datx, daty)
#' 
#' plot(dat)
#' LMK_ellipse(dat, ci = 90)
#' 
#'@export


LMK_ellipse <- function(dat, ci=c(67.5,90,95,99), linesCol = "black", fillCol = "grey", smoothness = 20){

##set smoothness of the circle (manually), the greater the number, the smoother it will be
sm <- smoothness
##90% inside curve = 4.605
##95% inside curve = 5.991
##99% inside curve = 9.210

if ( ci == 90){
  chi.v <- 4.605
} else if (ci ==95) {
  chi.v <- 5.991
} else if (ci == 99){
  chi.v <- 9.210
} else if (ci == 67.5){
  chi.v <- 2.25
}  else {
  stop("Invalid CI, please choose either 90,95, or 99")
}


cov.dat <- cov(dat)
cent <- t(colMeans(dat))

##trace and determinate of cov matrix
tr <- sum(cov.dat[1,1], cov.dat[2,2])

det <- ((cov.dat[1,1]*cov.dat[2,2]) - (cov.dat[1,2]*cov.dat[2,1]))

##determine eigen values
ei1 <- (tr + ((tr^2) - 4*det)^.5)/2
ei2 <- tr-ei1

##determine length of radii
ei.a <- (ei1^.5) * (chi.v^.5)
ei.b <- (ei2^.5) * (chi.v^.5) 


##calculate rotation of ellipse
###NOTE: these terms may need to be flipped
th <- atan2((ei1-cov.dat[1,1]), cov.dat[1,2])

##rotation matrix
q <- matrix(nrow = 2, ncol = 2)
q[1,1] <- cos(th)
q[2,2] <- cos(th)
q[2,1] <- sin(th)
q[1,2] <- sin(th)*-1

##generate points along ellipse
#first, generate equidistant points in radians along unit circle

circ <- matrix(nrow = 2*sm+1, ncol = 3)

for(i in 1:dim(circ)[1]){
  circ[i,1] <- (i-1) * (pi/sm)
  circ[i,2] <- (q[1,1] * ei.a * cos(circ[i,1]) + q[1,2] * ei.b * sin(circ[i,1])) + cent[1,1]
  circ[i,3] <- (q[2,1] * ei.a * cos(circ[i,1]) + q[2,2] * ei.b * sin(circ[i,1])) + cent[1,2]
}


if (is.null(fillCol)){
  lines(circ[,c(2,3)], col = linesCol)
  } else {
    polygon(circ[,c(2,3)], col = fillCol) ## Make sure to set alfa 
    lines(circ[,c(2,3)], col = linesCol, lwd = 1.5)
  }


}




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
  p <- dim(A)[1]
  k <- dim(A)[2]
  n <- dim(A)[3]
  lbls <- dimnames(A)[[3]]
  
  #check/format missing data
  if (anyNA(A)){
    A[is.na(A)] = 9999
    h = paste(1, paste(n,"L", sep=""), (p*k), 1, 9999, "Dim=3", sep = " ")} else {
      h = paste(1, paste(n,"L", sep=""), (p*k), 0, "Dim=3", sep = " ")
    }
  #write header line
  write(c(h, "\n"), file=filepath, ncolumns = 1, append= FALSE)
  
  #write specimen labels
  write(c(lbls, "\n"), file=filepath, ncolumns = 1, append= TRUE , sep = "\t")
  
  #write coords
  for(i in 1:n){
    write(c(t(A[,,i]), "\n"), file=filepath, append=T, sep = "\t", ncolumns = 3)
    
  }
  
}

#####


#' Swap Landmarks
#' 
#' 
#' Function to swap landmark sequences within an individual or across an entire array 
#' 
#' @param a a p x k x n matrix of landmark coordinates
#' @param l1 a numeric vector of landmark index (sequence number) to be changed
#' @param l2 a numeric vector of landmark index (Sequence number) to change to
#' 
#' 
#' @examples 
#' 
#' lmks <- cbind(1:10, 1:10, 1:10)
#' plot(lmks, col = rainbow(10), pch = 16) ## plot sequence (in 2 dimensions)
#' lmks2 <- LMK_swap(lmks, l1 = c(2,8), l2 = c(8,2)) ## flip landmarks 2,8
#' plot(lmks2, col = rainbow(10), pch = 16) ## show flipped landmarks
#' 
#' 
#' @export
#' @author Daniel Ehrlich
#' 
LMK_swap <- function(a, l1, l2){
  
  t.a <- a
  
  if (length(l1) != length(l2)){stop("Vectors must be the same length")}
  
  if(length(dim(t.a))==3){
    for (i in 1:length(l1)){
      
      t.a[l1[i],,] <- a[l2[i],,]
      t.a[l2[i],,] <- a[l1[i],,]
      
    }
  } else {
    for (i in 1:length(l1)){
      
      t.a[l1[i],] <- a[l2[i],]
      t.a[l2[i],] <- a[l1[i],]
    }
    }
  return(t.a)
  }

#####
#' Get sliders
#' 
#' Define a sliders table (geomorph) for a contour with fixed endpoints
#' 
#' Assumes all points slide sequentially between endpoints
#'  
#'@param c1 A vector of landmark numbers
#' 
#'@return Returns a table of sliders suitable for use with geomorph::gpagen(). This function assumes all points slide sequentially between the max/min landmarks in the sequence.
#' 
#' 
#' @examples 
#' 
#' lmknums <- 1:100
#' 
#' sli <- LMK_get_sliders(lmknums[80:100]) ## to treat landmarks 80 through 100 as sliders
#' 
#' sli
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
#' @return Returns a list (3L) containing the direct Euclidean distance for all individuals/landmarks, as well as summarized by individual and landmark
#' 
#' 
#' @examples 
#' 
#' lmks1 <- cbind(1:10, 1:10, 1:10)
#' lmks2 <- lmks1*10
#' 
#' LMK_euD(lmks1,lmks2) ## distance between each landmark (1 through 10)
#' 
#' LMK_euD(lmks1[1,], lmks2[1,]) ## distance between a single landmark
#' 
#' 
#' @export
#' 
LMK_euD <- function(A1, A2){
  
  #####For Arrays####
  
  ###For inter landmark distances (ie 2 points)
  if (length(A1)==3 & length(A2)==3){
    out <- ((A1[1] - A2[1])^2 + (A1[2] - A2[2])^2 + (A1[3]-A2[3])^2)^.5
    return(out)
      }
  
  if (length(dim(A1))==3 & length(dim(A2))==3){
    if (all(dim(A1) != dim(A2))){stop("Arrays must have the same extent")}
  
  
  l <- dim(A1)[[1]]
  k <- dim(A1)[[2]]
  n <- dim(A1)[[3]]
  
  if (is.null(dimnames(A1)[[1]])){
    lmklbl <- paste("lmk", 1:l, sep = ".")
  }
  lmklbl <- dimnames(A1)[[1]]
  
  if (is.null(dimnames(A1)[[3]])){
    indlbl <- paste("ind", 1:n, sep = ".")
  }
  indlbl <- dimnames(A1)[[3]]
  
  o1 <- matrix(data = NA,nrow = l, ncol = n, dimnames = list("lmk" = lmklbl, "ind" = indlbl))
  
  
  for (i in 1:n){
    #for each individual
    for (j in 1:l){
      #for each lmk
      
      #calculate distance
      if (k == 3){
        ## for 3 dimensional data
      o1[j,i] <- sqrt(sum(((A1[j,1,i]-A2[j,1,i])^2), ((A1[j,2,i]-A2[j,2,i])^2), ((A1[j,3,i]-A2[j,3,i])^2)))
      } else {
        ## for 2 dimensional data
        o1[j,i] <- sqrt(sum(((A1[j,1,i]-A2[j,1,i])^2), ((A1[j,2,i]-A2[j,2,i])^2)))
        
      }
    }
    
  }
  
  all.dif <- o1
  lmk.dif <- data.frame("avg.dif" = rowMeans(o1))
  ind.dif <- data.frame("avg.dif" = colMeans(o1))
  
  out <- list("all.dif"= all.dif, "by.lmk" = lmk.dif, "by.ind" = ind.dif)
  return(out)
  
  } 
  
  
  #####For Matrices#####
  ### note: this refers to the structure of the data (array vs matrix) rather than the dimensionality of the data (2D vs 3D landmarks)
  
  else if (length(dim(A1))==2 & length(dim(A2))==2){
    if (all(dim(A1) != dim(A2))){stop("Matrices must have the same extent")}
    
    l <- dim(A1)[[1]]
    
    o1 <- numeric(l)
    for (i in 1:l){
    #for each lmk
      
    #calculate distance
    o1[i] <- sqrt(sum(((A1[i,1]-A2[i,1])^2), ((A1[i,2]-A2[i,2])^2), ((A1[i,3]-A2[i,3])^2)))
        }
      return (o1)  
      }
      
  }



#####
#'Color Ramp
#'
#'Create a color ramp based on a grouping factor
#'
#'Creates color ramp for a given factor and returns a vector of color assignments for each individual using grDevices::rainbow(). Useful for plotting by group.
#'
#'@param grp A grouping variable, should be a factor
#'@param mute Logical, should color be muted. If true (default), color values will be generated with saturation value of .4, if false, values will have full saturation (1.0)
#'
#'@return Returns a vector of color assignments for each individual based on a grouping factor. 
#'
#'
#'
#'@examples
#'
#' data <- cbind(1:20, 20:1)
#' grps <- as.factor(rep(c("A","B","C","D","E"),4)) 
#' plot(data, col = LMK_colramp(grps), pch = 16) ##convert factor to rainbow pallete
#'@export
#'
#'@author Daniel Ehrlich
#'


LMK_colramp <- function(grp, mute = TRUE){
  
  f.grp <- as.factor(grp)
  #ensure vector is a factor to match index
  
  c.grp <- as.character(grp)
  #generate character list to recieve colors
  
  l <- length(levels(f.grp))
  if( mute == TRUE){
    cr <- grDevices::rainbow(l, s = 0.4, v = 1)
    #generate color ramp based on levels of data
    
  }else{
    cr <- grDevices::rainbow(l)
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
#'@examples
#'
#'require(geomorph)
#'data(plethodon)
#'
#'plotAllSpecimens(plethodon$land) ## plot all unaligned specimens
#'gpa <- gpagen(plethodon$land) ## align AND SCALE data
#'plotAllSpecimens(gpa$coords)
#'
#'backscaled <- LMK_bscale(gpa$coords, gpa$Csize)
#'
#'plotAllSpecimens(backscaled) ## data are aligned, but in their original size
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
#'Apply a function (f) to a set of data (d) by a grouping variable (grp). Calculates values within group (eg mean), as well as for the entire ungrouped sample (eg, grand mean). 
#'
#'@param d a dataset
#'@param f a function
#'@param grp a grouping variable
#'
#'@examples
#'data(iris)
#'LMK_grp_apply(d = iris[,1:4], f = mean, grp = iris$species) 
#'

#'@export
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

#' Intra class correlations and Chronbach's Alpha
#' 
#' Compute ICC and CA for each dimension (X,Y,Z) of arrays
#' 
#' @param obs1 a 3D landmark array
#' @param obs2 a 3D landmark array
#' 
#' @return Returns a table of ICC coefficients, associated p-values, and Chronbach's Alpha
#' 
#' 
#' 
#' @export
#' 
#' 
LMK_iccA <- function(obs1, obs2){
  
 
  
  #inter- or intra- oberserer error
  #where obs1 and obs2 are two identical landmark data sets
  #function to perform inter class correlation and reports F- and p- values, Chronbach's alpha, as well as euclidean distance for a set of landmarks
  #currently assumes: .dta format, 3 dimensions, 34 lmks, object symmetry
  #just get it from psych package
  
  
 
  l <- dim(obs1)[[1]]
  k <- dim(obs1)[[2]]
  n <- dim(obs1)[[3]]
  
  if (k==3){
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
   
  } else {
    
    ##### need to streamline 2D data steps ####
    
    out <- data.frame("ICC.X" = numeric(l), "f.X" = numeric(l),"sig.X" = numeric(l), "ICC.Y" = numeric(l),"f.Y" = numeric(l), "sig.Y" = numeric(l), "ICC.Z" = numeric(l), "f.Z"=numeric(l),"sig.Z" = numeric(l),"CrA.X" = numeric(l), "CrA.Y" = numeric(l), "CrA.Z" = numeric(l),row.names = dimnames(obs1)[[1]])
    for (r in 1:l){
      tx <- psych::ICC(as.matrix(cbind(obs1[r,1,], obs2[r,1,])), lmer=FALSE)
      ty <- psych::ICC(as.matrix(cbind(obs1[r,2,], obs2[r,2,])), lmer=FALSE)
      tz <- NA
      ax <- psych::alpha(as.matrix(cbind(obs1[r,1,], obs2[r,1,])))
      ay <- psych::alpha(as.matrix(cbind(obs1[r,2,], obs2[r,2,])))
      az <- NA
      
      out[r,1] <- tx$results$ICC[3]
      out[r,2] <- tx$results$F[3]
      out[r,3] <- tx$results$p[3]
      
      out[r,4] <- ty$results$ICC[3]
      out[r,5] <- ty$results$F[3]
      out[r,6] <- ty$results$p[3]
      
      out[r,7] <- NA
      out[r,8] <- NA
      out[r,9] <- NA
      
      out[r,10] <- ax$total[1]
      out[r,11] <- ay$total[1]
      out[r,12] <- NA
    }
    
  }
  return(out)
}

#' Launch interactive screening plots
#' 
#' GUI interface to choose DTA file, and various classifiers to conduct preliminary analysis.
#' 
#' @export
#' 

LMK_Screen <- function() {
  appDir <- system.file("apps", "LMK_Screen", package = "lmkCHK")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `lmkCHK`.", call. = FALSE)
  }
  
  shiny::runApp(appDir, display.mode = "normal")
}


