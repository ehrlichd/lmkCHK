#Function to plot out all HSV values to assist in choosing colors

########################
#####HSV Color Tool#####
########################
#First, load the helper functions by running the two lines below

hsv.plot = function(h, q){
  #where h is the desired Hue (0:1) and q is the plotting increment (<1)
  ## x = .01 provides a seemingly continous color spectrum
  #q = .1
  x = c(rep(seq(0, 1, by = q), (1/q)+1), seq(0,1, by = q))
  y = c(rep(seq(0, 1, by = q), each = (1/q+1)), seq(1,0, by = -q))
  col.tab = data.frame(cbind(x,y))
  plot(0:1, 0:1, xlab = "Saturation", ylab = "Value", type = "n", main = c("Hue =", as.character(h)))
  for (i in 1:length(col.tab[,1])){
    #points(as.list(col.tab[i,1]), col = hsv(h, as.list(col.tab[i,])), pch = 16)
    points(col.tab[i,1], col.tab[i,2], col = hsv(h, col.tab[i,1], col.tab[i,2]), pch = 16)
  }
}
hsv.all = function(x){
  #where x is the increment
  library(rgl)
  
  x = .1
  h = rep(seq(0, 1, by = x), each = ((1/x)+1)^2)
  s = rep(rep(seq(0, 1, by = x), each = 1+(1/x), (1/x)+1))
  v = rep(seq(0, 1, by = x), ((1/x)+1)^2)
  
  #h = c(rep(seq(0, 1, by = x), 1/x), rep(seq(0, 1, by = x), each = 1/x), rep(seq(0, 1, by = x), each = 1/x), seq(1,0, by = -x))
  #s = c(rep(seq(0, 1, by = x), each = 1/x), rep(seq(0, 1, by = x), 1/x), rep(seq(0, 1, by = x), each = 1/x), seq(1,0, by = -x))
  #v = c(rep(seq(0, 1, by = x), each = 1/x), rep(seq(0, 1, by = x), each = 1/x), rep(seq(0, 1, by = x), 1/x), seq(1,0, by = -x))
  col.tab = data.frame(cbind(h,s,v))
  plot3d(0:1,0:1,0:1, xlab = "Hue", ylab = "Saturation", zlab = "Value", type = "n")
  
  for (i in 1:length(col.tab[,1])){
    #points3d(col.tab[i,1], col.tab[i,2], col.tab[i,3], col = hsv(col.tab[i,1], col.tab[i,2], col.tab[i,3]))
    spheres3d(col.tab[i,1], col.tab[i,2], col.tab[i,3], col = hsv(col.tab[i,1], col.tab[i,2], col.tab[i,3]), radius = x)
  }
}

####Intro####

#This script contains a handful of tools that allow users to explore color space to aid in creating quality visuals in R

#R allows several ways of specifying colors including by name, RGB, and Hex-code. This script uses another, HSV.
#
#HSV stands for Hue, Saturation, Value and can be specified using the hsv() function
#Press Ctrl+enter to execute a line of code
?hsv

#As you can see, hsv() takes 3 arguments, Hue(h), Saturation(s), and Value(v)
#Each of these values ranges from 0 to 1 and represent different colors, tints, and shades, respectively. S and V values of 1 correspond to a "pure" color.
#EG hsv(h=0, s=1, v=1) is pure red, hsv(h=0, s=.5, v=1) is light red, and hsv(h=0, s=1, v=.5) is dark red.

#Run each line below and watch what happens in the plot
plot(0:1,0:1, type = "n", xlab = "Saturation", ylab = "Value", main = "Cyan, h=0.5")
points(1,1, col = hsv(h=.5,s=1,v=1), pch=16)
points(.5,1, col = hsv(h=.5,s=.5,v=1), pch=16)
points(1,.5, col = hsv(h=.5,s=1,v=.5), pch=16)

#You can explore the HSV color space in either 2D or 3D, check out the following sections below

####2D Colorspace#####

#Choosing a Hue
hsv.2d.all = function(x){
  #where x is the color increment
  x=.1
  h = rep(seq(0, 1, by = x), each = (2*(1/x)+1))
  y.s = rep(seq(-1,1,by=x), ((1/x)+1))
 
  
  col.tab = cbind(h, y.s, s, v)
  
  plot(x= 0:10, y= 0:10, type = "n", xlab = "Hue", ylab = "Value:Saturation" , main = "HSV Colors")
  
points(seq(0,10, by = .1),seq(0,10, by=.1), col = terrain.colors(101), pch =16)  
  
}









hsv.plot(.5, .01)
#####Plot3d#####

hsv.all(.1)

