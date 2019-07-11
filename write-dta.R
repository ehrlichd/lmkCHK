writeland.nts = function(x, dta){
  #write 3D array landmark file to .dta/NTSys file
  #where x is data and dta is the file to be written, including file exstension
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
