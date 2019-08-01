#'Plot Outliers
#'
#'Function to identify outliers visually
#'
library(geomorph)
dat = readland.nts("cran.dta")

dim(dat)
a.dat = gpagen(dat, ProcD = F)

p = dim(dat)[[1]]
n = dim(dat)[[3]]

a.dat$consensus
r = range(a.dat$coords)
plot3d(NULL, xlim = r, ylim = r, zlim = r, xlab = "", axes = F)
points3d(a.dat$consensus, col = "red")

for (i in 1:n){
  points3d(a.dat$coords[,,i], col = "grey")
}
spheres3d(a.dat$consensus, col = "red",radius = 1)
dim(a.dat$coords)

dist = matrix(NA, nrow = 610, ncol = 2)
dist[,1] = dimnames(dat)[[3]]
rm(i,j)
tab = NULL
for(i in 1:n){
  for(j in 1:p){
    tab[j] = sqrt(sum(((a.dat$consensus[j,1] - a.dat$coords[j,1,i])^2), ((a.dat$consensus[j,1] - a.dat$coords[j,1,i])^2), ((a.dat$consensus[j,1] - a.dat$coords[j,1,i])^2)))
    if (j == p){dist[i,2] = sum(tab)}
  }
  
  
}
View(dist[order(dist[,2], decreasing = T),])
y = range(dist[,2])
x = c(1,n)
rm(i,j)
plot(1:610, dist[order(dist[,2], decreasing = T),2], type = "n")
for (i in 1:610){
  d = as.numeric(dist[order(dist[,2], decreasing = T),2][i])
  if (d > mean(as.numeric(dist[,2])) + 3*sd(as.numeric(dist[,2]))){
    points(i, d, pch = 16, col = "red")
    text(i*10,d,labels = dist[i,1], col = "red", adj = c(0,0))
    } else {
      
      points(i, d, pch = 16, col = "grey")
      }
  }

 abline(h=mean(as.numeric(dist[,2])) + 3*sd(as.numeric(dist[,2])), col = "red")

View(dist[,2])
dist
t.dat = dist[order(dist[,2],decreasing = T),]
View(t.dat)
View(dist[sort])
length(a.dat$consensus[1,]) == length(a.dat$coords[1,,1])
