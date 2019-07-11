#Organizing notes for LMK.chk, an R package to assist in 2D and 3D geometric morphometrics, and other landmark based applications

#draws on:
#Morpho
#geomorph
#shapes   -  does it???
#psych

####description####



####functions####

#####
writeland.nts = function(x, dta){
  #write 3D array landmark file to .dta/NTSys file
  #where x is data and dta is the file to be written, including file exstension
  #code written by DEE
  p = dim(x)[1]
  k = dim(x)[2]
  n = dim(x)[3]
  lbls = dimnames(x)[[3]]
  
  #check/format missing data
  if (anyNA(x)){
    x[is.na(x)] = 9999
    h = paste(1, paste(n,"L", sep=""), (p*k), 1, 9999, "Dim=3", sep = " ")} else {
      h = paste(1, paste(n,"L", sep=""), (p*k), 0, "Dim=3", sep = " ")
    }
  write(c(h, "\n"), file=dta, ncol=1, append=F)
  
  #write specimen labels
  write(c(lbls, "\n"), file=dta, ncol=1, append=T, sep = "\t")
  
  #write coords
  for(i in 1:n){
    write(c(t(x[,,i]), "\n"), file=dta, append=T, sep = "\t", ncol=3)
    
  }
}

#####

lmk.swap = function(l1, l2){
  #function to swap landmark configurations either by list or single landmarks
  #change the order, or set LMKs to NA
  #code written by DEE
  
  }
#####

get.sliders = function(c1){
  #function to define a sliders table (for geomorph) based on a contour with fixed enpoints
  #code written by DEE
  
  t = as.numeric(as.matrix(c1)) # convert list to vector
  l = length(t)
  l2 = l-1
  
  sl = matrix(data = cbind(t[-c(l2,l)], 
                           t[-c(1,l)], 
                           t[-c(1,2)]), 
              nrow = l-2, 
              ncol = 3, 
              dimnames = list("row" = paste("lmk", seq(2, l-1), sep = "." ), "col" = c("before", "slide", "after")))
  
  return(sl)
  
}

#####
rm(multi.lmk)

obs.error = function(d1, d2, p = NULL){
    #inter-/intra- observer for any 2 sets of landmark data
    #code written by DEE
    
    #function to compile all parts to LMK reliability. Expects2 3d arrays of any number of landmarks
    #1 euclidian difference between individual landmarks and summarized by individual/landmark
    #2inter class correlation (paired anova(???)) accross all individuals by each dimension of each landmark
    ##ie lmk1.x1 , lmk1.y1, lmk.z1, lmk2.x1, etc
    #3chronbach's alpha
    
    
    library(geomorph)
    library(shapes)
    library(psych)
    library(lme4)
    
    
    #Fallten arrays
    raw = rbind(two.d.array(d1),two.d.array(d2))
    
    ##Miller Cranial Dataset lmks: c("g", "n","prn","sn","ls","sto","li","sl","pg","gn", "en_R", "ps_R","ex_R","pi_R","en_L", "ps_L","ex_L","pi_L", "al_R","ac_R","sbal_R","c_R", "al_L","ac_L","sbal_L","c_L","ch_R","cph_R", "ch_L", "cph_L","obi_R","obi_L", "x7","x8")
    #l landmarks
    l = dim(d1)[1]
    #n individuals 
    n = dim(d1)[3]
    
    if (!is.null(p)){
     lbl = list("lmk" = p, "xyz" = c("x","y","z"), "ind" = dimnames(raw)[[1]]) 
    }else{
      lbl = list("lmk" = paste("lmk",1:l), "xyz" = c("x","y","z"), "ind" = dimnames(raw)[[1]])
    }
    

    
    ###1: procrustes superimposition WITHOUT scaling
    a.raw = arrayspecs(raw, l, 3)
    #requires shapes package
    aln = procGPA(a.raw,scale = FALSE)$rotated
    dimnames(aln) = lbl
    
    ###2:Euclidean Distances
    ar1 = aln[,,1:n]
    ar2 = aln[,,(1:n)+n]
    eu.d = function(d1, d2){
      #function to calculate the euclidean distances between 2 arrays of LMK data
      #assumes arrays are the same dimensions and in the format [p,k,n] where p = landmarks, k = dimensions, and n = individuals
      #assumes dimnames are set
      
      
      #create output matrix of l rows, n columns
      o1 = matrix(data = NA,l,n)
      
      lbl = dimnames(d1)[[3]]
      lmk = dimnames(d1)[[1]]
      
      for (c in 1:n){
        #for each individual
        for (r in 1:l){
          #for each lmk
          o1[r,c] = sqrt(sum(((d1[r,1,c]-d2[r,1,c])^2), ((d1[r,2,c]-d2[r,2,c])^2), ((d1[r,3,c]-d2[r,3,c])^2)))
        }
        
      }
      colnames(o1) = lbl
      rownames(o1) = lmk
      all.dif = o1
      lmk.dif = data.frame("avg.dif" = rowMeans(o1))
      ind.dif = data.frame("avg.dif" = colMeans(o1))
      out = list("all.dif"= all.dif, "lmk.dif" = lmk.dif, "ind.dif" = ind.dif)
      return(out)
    }
    
    d.out = eu.d(ar1, ar2)
    
    
    ###3: Intraclass Correlation
    icc.tab = function(d1, d2){
      #inter- or intra- oberserer error
      #where d1 and d2 are two identical landmark data sets
      #function to perform inter class correlation and reports F- and p- values, chronbach's alpha, as well as euclidean distance for a set of landmarks
      #currently assumes: .dta format, 3 dimensions, 34 lmks, object symmetry
      #just get it from psych package
      library(psych, lme4)
      
      l = dim(d1)[1]
      d = 3
      out = data.frame("ICC.X" = numeric(l), "f.X" = numeric(l),"sig.X" = numeric(l), "ICC.Y" = numeric(l),"f.Y" = numeric(l), "sig.Y" = numeric(l), "ICC.Z" = numeric(l), "f.Z"=numeric(l),"sig.Z" = numeric(l),"CrA.X" = numeric(l), "CrA.Y" = numeric(l), "CrA.Z" = numeric(l))
      row.names(out) = dimnames(d1)[[1]]
      
      
      for (r in 1:l){
        t.x = ICC(as.matrix(cbind(d1[r,1,], d2[r,1,])), lmer=FALSE)
        t.y = ICC(as.matrix(cbind(d1[r,2,], d2[r,2,])), lmer=FALSE)
        t.z = ICC(as.matrix(cbind(d1[r,3,], d2[r,3,])), lmer=FALSE)
        a.x = alpha(as.matrix(cbind(d1[r,1,], d2[r,1,])))
        a.y = alpha(as.matrix(cbind(d1[r,2,], d2[r,2,])))
        a.z = alpha(as.matrix(cbind(d1[r,3,], d2[r,3,])))
        
        out[r,1] = t.x$results$ICC[3]
        out[r,2] = t.x$results$F[3]
        out[r,3] = t.x$results$p[3]
        
        
        out[r,4] = t.y$results$ICC[3]
        out[r,5] = t.y$results$F[3]
        out[r,6] = t.y$results$p[3]
        
        out[r,7] = t.z$results$ICC[3]
        out[r,8] = t.z$results$F[3]
        out[r,9] = t.z$results$p[3]
        
        out[r,10] = a.x$total[1]
        out[r,11] = a.y$total[1]
        out[r,12] = a.z$total[1]
        
      }
      
      return(out)
    }
    
    icc.out = icc.tab(ar1, ar2)
    
    
    ###5. Formatting
    pretty = as.data.frame(cbind(d.out$lmk.dif, icc.out[,c(1,3,4,6,7,9:12)]))
    pretty$avg.ICC = rowMeans(pretty[,c(2,4,6)])
    pretty$avg.sig = rowMeans(pretty[,c(3,5,7)])
    pretty$avg.ChrA = rowMeans(pretty[,c(8:10)])
    out = list("eu.d" = d.out, "icc" = icc.out, "quick" = pretty)
    
    return(out)
  }
  


#####
col.ramp = function(grp, mute = TRUE){
  #generate a rainbow color ramp for any number of levels and output a vector of colors for the data 
  #add in conditional statement to lighten, ie col.ramp(grp, opaq=T) #if true, then default rainbow, if F then s,v ~ 0.8
  
  
  f.grp = as.factor(grp)
  #ensure vector is a factor to match index
  
  c.grp = as.character(grp)
  #generate character list to recieve colors
  l = length(levels(grp))
  if( mute == TRUE){
    cr = rainbow(length(levels(f.grp)), s = 0.4, v = 1)
    #generate color ramp based on levels of data
    
  }else{
    cr = rainbow(length(levels(f.grp)))
    #generate color ramp based on levels of data
    
  }
  
  for (i in 1:l){
    c.grp[as.integer(f.grp) == i] = cr[i]
  }
  return(c.grp)
}

#####

b.scale = function(A, cs){
  #where A is a pxkxn array, and cs is a vector of centroid sizes
  al = dim(A)[[3]]
  cl = length(cs)
  if (length(al) != length(cl)){
    print("data mismatch")
  } else {
    for (i in 1:al){
      A[,,i] = A[,,i]*cs[i]
    }
  }  
  return(A)
}

#####
grp.apply = function(d, f, grp){
  #apply a function (f) to a dataframe (d) using a grouping variable (grp)
  #future implementation: n=T/F to indicate if you want sample size calculated
  out = as.data.frame(matrix(nrow = (length(levels(grp))+1), ncol = length(d)*2))
  c.names = c(colnames(d), paste(colnames(d),"n=", sep = "."))
  id = c(1:length(d), (1:length(d)+0.5))
  colnames(out) = c.names[order(id)]
  row.names(out) = c(levels(grp), paste("G", deparse(substitute(f)), sep = "."))
  
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

CV.2D.shape = function(CVA, variate = 1, scale = 5){
  #Function to visualize CVA shape changes of landmarks 
  #where CVA is an object of class 'CVA' from Morpho
  
  ref = CVA$Grandm
  target = scale*matrix(CVA$CVvis[,variate], nrow(CVA$Grandm), ncol(CVA$Grandm)) + CVA$Grandm
  
  plot(ref, asp=1, xlab = "", ylab = "", axes = F)
  for (i in 1:nrow(ref)){
    lines(rbind(ref[i,], target[i,]), col = "red")
  }
  
}


####Note: scale makes spheres too big, add an inverse scaling factor to those plots.
####Also, open3d means multi plots are impossible
CV.3D.shape = function(CVA, variate = 1, scale = 5, links = NULL){
  ref = CVA$Grandm
  target = scale*matrix(CVA$CVvis[,variate], nrow(CVA$Grandm), ncol(CVA$Grandm)) + CVA$Grandm
  r = range(ref[,1])*1.1
  
  open3d()
  rgl.clear()
  plot3d(ref, xlim = r, ylim = r, zlim = r, axes = F, type = "p", size = abs(scale)*1.5)
  
  for (i in 1:nrow(ref)){
    lines3d(rbind(ref[i,],target[i,]), col = "red", lwd = 1.5)
  }
  
}


form.PCA.plot = function(A, xPC = 1, yPC = 2, zPC = NULL, grp){
  ###Function to plot formspace PCA for any GPA object (of gpagen())
  ##Where A is a 3D array of GPA aligned coordinates
  f.PCA = plotTangentSpace(A, warpgrids = F)
  grp = as.factor(grp)
  l = length(levels(grp))
  
  ##Calculate means
  Mean.tab = matrix(nrow = l+1, ncol = 10)
  for (i in 1:l){
    Mean.tab[i,] = colMeans(f.PCA$pc.scores[as.integer(grp) == i, 1:10 ])
  }
  Mean.tab[l+1,] = colMeans(f.PCA$pc.scores[,1:10]) 
  rownames(Mean.tab) = c(levels(grp),"GrandM")
  colnames(Mean.tab) = paste("PC",1:10, sep=".")
  
  if (is.null(zPC)){
    plot(f.PCA$pc.scores[,c(xPC, yPC)], col = col.ramp(grp), pch = 16)
    points(Mean.tab[,c(xPC, yPC)], pch = 21, bg = c(rainbow(l, s = .4), "grey"), cex = 1.5)
    legend("topleft", legend = c(levels(grp), "G.mean"), pch = 16, col = c(rainbow(l, s=.4), "grey"))
  }
  
}

