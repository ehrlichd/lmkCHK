#' Inter-/Intra-Oberver Error
#' 
#' Calculate the inter- or intra-observer error between two sets of LMK data
#' 
#' Computes three measures of error between two sets of LMK data. Direct Euclidean distance, intraclass correlation (ICC), and Chronbachs alpha. All values are summarized by landmark, and by individual to investigate/identify systemic error. Output is a summary table, or a list (3L) containing the summary table, as well as additional raw and transformed values.
#' 
#' @param obs1 A set of observations of landmark data
#' @param obs2 A second set of observations of landmark data 
#' @param lmk.lbl A charactor vector containing informative names for landmarks. If Null (default) LMKs will be labeld as 'lmk.1, lmk.2, lmk.3, ...'
#' @param full Logical, should a full output be returned. By default (FALSE), only a summary table is returned. If TRUE, a list (3L) is returned with additional infomation.
#'  
#' @return Test returns a table summarizing euclidean distance, ICC and chronbach's alpha. If full = TRUE, then the function returns a list (3L) of additional raw/transformed values used in calculations
#' 
#' @export
#' 
#' @author Daniel Ehrlich
#' 
#' 
obs_error <- function(obs1, obs2, lmk.lbl = NULL, full = FALSE){
 
  #1 euclidian difference between individual landmarks and summarized by individua and landmark
  
  #2inter class correlation  accross all individuals by each dimension of each landmark
  ##ie lmk1.x1 , lmk1.y1, lmk.z1, lmk2.x1, etc
  #3chronbach's alpha
  
  

  #Flatten arrays
  raw <- rbind(geomorph::two.d.array(obs1),geomorph::two.d.array(obs2))
  
  #l landmarks
  l <- dim(obs1)[1]
  #n individuals 
  n <- dim(obs1)[3]
  
  if (!is.null(lmk.lbl)){
    lbl <- list("lmk" = lmk.lbl, "xyz" = c("x","y","z"), "ind" = dimnames(raw)[[1]]) 
  }else{
    lbl <- list("lmk" = paste("lmk",1:l), "xyz" = c("x","y","z"), "ind" = dimnames(raw)[[1]])
  }
  
  
  
  ###1: procrustes superimposition WITHOUT scaling
  araw <- geomorph::arrayspecs(raw, l, 3)
  #requires shapes package
  
  aln <- geomorph::gpagen(araw,ProcD = F)
  aln <- bscale(aln$coords, aln$Csize)
  
  dimnames(aln) <- lbl
  
  ###2:Euclidean Distances
  ar1 <- aln[,,1:n]
  ar2 <- aln[,,(1:n)+n]
  
  
  d.out <- euD(ar1, ar2)
  
  
  ###3: Intraclass Correlation
  
  
  
  icc.out <- iccA(ar1, ar2)
  
  
  ###5. Formatting
  pretty <- as.data.frame(cbind(d.out$lmk.dif, icc.out[,c(1,3,4,6,7,9:12)]))
  pretty$avg.ICC <- rowMeans(pretty[,c(2,4,6)])
  pretty$avg.sig <- rowMeans(pretty[,c(3,5,7)])
  pretty$avg.ChrA <- rowMeans(pretty[,c(8:10)])
  if (full == TRUE){
    out <- list("eu.D" = d.out, "icc" = icc.out, "quick" = pretty)
  } else {out <- pretty}
  
  return(out)
}

