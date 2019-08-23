###Prelim screening and analysis for Cranial SOF dataset. Fossil + Papionin
library(lmkCHK) #loads geomorph and morpho

#SOF = list("LMpairs" = SOF_pair, "curves" = sli, "surf" = surf)
#save(SOF, file= "SOF.RData")
load("SOF.RData")
#save(clas, file = "clas.RData")
load("clas.RData")
###NOTE: FOR SOF DATA COLLECTION facial foramen should be filled prior to landmark colelction#####

###LMK screening sequence
#1. read in dtas
#2. symmetrize
#3. plot/calc outliers
#4. for the top n outliers, LMK_screen
####symbolize vectors based on length? PC loadings???
###highlight landmarks with high loadings in top T PCs

#5. flag ones for further inspection in LME
#6. PCA plot for outleirs as well

normscale = function(x){
  t = (x - min(x))/(max(x)-min(x))
  return(t)
}


#####Run through screening criteria individually first to spot outliers, then combine for actual analysis


#####
#####Preliminary Screening#####
#####

#####1: Read in each DTA indvidually#####


fossil <- readland.nts("Fossil ZM Patch.dta")
fossil <- fossil[,,-c(1,2,35)] ##Remove extant atlas; fragmentary inds

PHK_F <- readland.nts("SOF_PHK-F.dta")
PHK_J <- readland.nts("SOF_PHK-J.dta")
PHK_J<- PHK_J[,,-15] ##Remove stage 0

PHK_M <- readland.nts("SOF_PHK-M.dta")

PHU_F <- readland.nts("SOF_PHU-F.dta")
PHU_J <- readland.nts("SOF_PHU-J.dta")
PHU_M <- readland.nts("SOF_PHU-M.dta")



#####1.1: Restructure data#####
#restructure raw data if needed, create/read in classifiers

all.ls = list("fossil" = fossil, "PHK_F" = PHK_F, "PHK_J" = PHK_J, "PHK_M" = PHK_M, "PHU_F" = PHU_F, "PHU_J" = PHU_J, "PHU_M" = PHU_M)

a = c(
  dimnames(fossil)[[3]],
  dimnames(PHK_F)[[3]],
  dimnames(PHK_J)[[3]],
  dimnames(PHK_M)[[3]],
  
  dimnames(PHU_F)[[3]],
  dimnames(PHU_J)[[3]],
  dimnames(PHU_M)[[3]]
)





###set parsable names
#######NOTE: This step should occur now, but will be implemented later to get on with SOP#####

###pars names can be set once everything is combined in 1 dataset, until then, screen with museum names

###Additional Classifiers

class = c(rep("fos", 25), rep("Kf", 20), rep("Kj", 20), rep("Km", 21), rep("Uf", 6), rep("Uj",31), rep("Um",21))

#####1.2: Define Landmarks####

##Paired LMs
#SOF_pair = matrix(c(10:15, 241:465, 4:9, 16:240), nrow = 231, ncol = 2)
#now contained in SOF.RData

##Define Countour semi-landmarks
#sli = psych::read.clipboard() ##slider table made by hand in excel
#now contained in SOF.RData
##Surface semi-landmarks

#surf = 1:225 #all landmarks
#surf = surf[!(1:225 %in% sli[,2])] #drop countours
#surf = surf[-c(1,15,106,211, 225)] #drop fixed LMKs
#now contianed in SOF.Rdata

#####





###2: Cleaning and evaluating seperately#####
sym.ls = NULL

for (i in 1:7){
  sym.ls[[i]] <- lmkCHK::LMK_sym(all.ls[[i]], SOF$LMpairs)
}

names(sym.ls) <- names(all.ls)

###Manualy tweak some fosil landmarks to aid symmestrization

t = fixLMmirror(fossil, SOF$LMpairs)
t[1:3,,15] = NA ##TP11
t[c(1,3),,31] = NA ##56614
t[1,,32] = NA ##56615
t[3,,21] = NA ##KA1562326

t = fixLMtps(t)
sym.ls$fossil = symmetrize(t$out, SOF$LMpairs)
rm(t)


cmb = NULL
for (i in 1:7){
  cmb = rbind(cmb, two.d.array(sym.ls[[i]]))
}
cmb = arrayspecs(cmb, p = 465, k = 3)

View(dimnames(cmb)[[3]])
#####Running outlier screening by hand for some reason
#####Manually itterate through sym.ls
#screen = function(dta, sym = TRUE, gpa = TRUE){

#t = readland.nts("fossil.bsc")

#t = readland.nts("PHK_F.bsc")
#t = readland.nts("PHK_J.bsc")
#t = readland.nts("PHK_M.bsc")

#t = readland.nts("PHU_F.bsc")
#t = readland.nts("PHU_J.bsc")
#t = readland.nts("PHU_M.bsc")
t = readland.nts("cmb.bsc")
#####Parsname to specname#####
for (i in 1:dim(t)[3]){
  dimnames(t)[[3]][i] = as.character(clas[clas$parsname==dimnames(t)[[3]][i],2])
}
dimnames(t)[[3]]

#####specnum to parsname#####
for (i in 1:dim(t)[3]){
  dimnames(t)[[3]][i] = as.character(clas[clas$specnum==dimnames(t)[[3]][i],3])
}

dimnames(t)[[3]]

 
#####Screening#####
  a = gpagen(t, curves = SOF$curves, surfaces = SOF$surf) 
  
  pca = prcomp(two.d.array(a$coords))
  
  
  r = range(a$coords)
  p = dim(a$coords)[1]
  n = dim(a$coords)[3]
  
  
  out1 = LMK_plotoutliers(a$coords, gpa = F, plotALL = F)
  legend("topright", legend = c(round(as.numeric(out1$ind.info[1:5,2]),2), out1$ind.info[1:5,1]), cex = .6, ncol = 2, bty = "n")
  dev.copy(png, "out1_shape_.png")
  dev.off()
  
  #View(out1$ind.info)
  d = as.numeric(out1$ind.info[,2])
  col = normscale(d)
  rm(i)
  for (i in 1:n){
    col[i] = rainbow(1, .25+(normscale(d)[i]/2), .25+(normscale(d)[i]/2),start = 0, end = 0.1)
  }
  
  
  plot(1:10, summary(pca)$importance[2,1:10], main = "PC Loadings", ylab = "Variance", xlab = "PC", las = 1)
  legend("topright", legend = c(round(summary(pca)$importance[3,1:5],2), paste("PC",1:5)), cex = .8, ncol = 2)
  dev.copy(png, "out2_loadings_.png")
  dev.off()
  
  for (i in 1:2){
    plot(pca$x[,c(i,i+1)], col = col, pch = 16, main = paste("Outliers on PC", i,"and" ,i+1, sep = " "),type = "n")
    text(pca$x[,c(i, i+1)], labels = dimnames(a$coords)[[3]], cex = .6, col = col)
    dev.copy(png, paste("out", i+2, "_text_",".png", sep = ""))
    dev.off()
    
    plot(pca$x[,c(i,i+1)], col = col, pch = 16, main = paste("Outliers on PC", i,"and" ,i+1, sep = " "))
    dev.copy(png, paste("out", i+2,"_dots_" ,".png", sep = ""))
    dev.off()
  }
  
  #return(out1)
  LMK_writeland_nts(t, "pars.bsc")
  BSC.cmb = t
  
  
  #####
  #####
  #####End Basic#####
  
  as.factor(substr(dimnames(cmb)[[3]],start = 1,stop = 3))
  
  collist = LMK_colramp(as.factor(substr(dimnames(cmb)[[3]],start = 1,stop = 3)), mute = F)
  
  collist = c(
    rep("grey", 35), #fossil
    rep(rainbow(3,start = .1, end = 0.15)[1],20), #phk f
    rep(rainbow(3,start = .1, end = 0.15)[2],19), #phk j
    rep(rainbow(3,start = .1, end = 0.15)[3],19), #phk m
    
    rep(rainbow(3,start = .5, end = 0.6)[1],6),
    rep(rainbow(3,start = .5, end = 0.6)[2],31),
    rep(rainbow(3,start = .5, end = 0.6)[3],21))
  
  dim(BSC.cmb)
  
  
  ####Plot individuals as text
  plot(pca$x[,1:2], pch = 16, col = collist, type = "n")
  text(pca$x[,1:2], labels = substr(dimnames(BSC.cmb)[[3]],start = 1,stop = 3), cex = .8, col = collist)
  
  
  ###
  plot(pca$x[,c(2,3)], type = "n")
  points(pca$x[36:151,c(2,3)], pch = 16, col = collist[36:151], cex = 2)
  
  text(pca$x[1:35,c(2,3)], labels = dimnames(t)[[3]][1:35], cex = .4)
  legend("topleft", legend = levels(clas$dataset), pch = 16, col = unique(collist), cex = .8, bty = "n")
  
  pairs(pca$x[,1:5], col = collist, pch = 16)
  LMK_compare_two(a$consensus, a$coords[,,1])
  
  
  View(dimnames(BSC.cmb)[[3]])
  View(pca$rotation)
  pc.shape = NULL
  
  pc.shape = BSC.cmb
  for (i in 1:151){
    pc.shape[,,i] = matrix(pca$rotation[,i], ncol = 3, byrow = T)
  }
  View(pc.shape[,,1])
  LMK_compare_two(a$consensus, pc.shape[,,1])
  
  t = matrix(pca$rotation[,1], ncol = 3, byrow = T)
  View(t)
  #####SAVE GOOD DATA#####
  names[132:152,2] == dimnames(BSC.PHU_M)[[3]]
  
  
  
  BSC.cmb = LMK_bscale(a$coords, a$Csize)
  View(dimnames(BSC.cmb)[[3]])
  dimnames(BSC.cmb)[[3]] %in% names$specnum
  
  names[dimnames(BSC.cmb)[[3]][2] %in% names$SOFname,]

  names$specnum %in% dimnames(BSC.cmb)[[3]]
  dimnames(BSC.cmb)[[3]] = names$parsname
  
  #####One final combined analysis#####
  
  cmb = arrayspecs(rbind(
    two.d.array(BSC.fossil),
    two.d.array(BSC.PHK_F),
    two.d.array(BSC.PHK_J),
    two.d.array(BSC.PHK_M),
    two.d.array(BSC.PHU_F),
    two.d.array(BSC.PHU_J),
    two.d.array(BSC.PHU_M)), p = 465,k = 3
    
  )
  
  
  anyNA(cmb)
  
  
  
  LMK_writeland_nts(BSC.cmb[,,1:35], "fossil.bsc")
  LMK_writeland_nts(BSC.cmb[,,36:55], "PHK_F.bsc")
  LMK_writeland_nts(BSC.cmb[,,56:74], "PHK_J.bsc")
  LMK_writeland_nts(BSC.cmb[,,75:93], "PHK_M.bsc")
  
  LMK_writeland_nts(BSC.cmb[,,94:99], "PHU_F.bsc")
  LMK_writeland_nts(BSC.cmb[,,100:130], "PHU_J.bsc")
  LMK_writeland_nts(BSC.cmb[,,131:151], "PHU_M.bsc")
  
  dimnames(BSC.cmb)[[3]]
  #####Plot ALL####
  
  
  
  plot3d(NULL, xlim = r, ylim = r, zlim = r, axes = F)
  
  for (i in 1:n){
    points3d(a$coords[,,i], col = rainbow(n, s = .8)[i])
    }
  
  #####Plot A specific indvidual####
  
  d = 11
  clear3d()
  plot3d(sym.ls$fossil[,,d], axes = F, xlim = range(sym.ls$fossil[,,d]), ylim = range(sym.ls$fossil[,,d]), zlim = range(sym.ls$fossil[,,d]))
  
  LMK_writeland_nts(sym.ls$fossil, "fossil.sym")
  
  #Write out BSC data
  #Combine all BSC into BSC all for a final gpagen/screening
  
 # }

#####END####
  
  
  
  
#####Analysis####
  #1 regress on CS
  

#####
#####Outliers list#####
#####
#####Running through outlier list doesn't seem to work
outliers <- NULL
for ( i in 1:7){
  outliers <- append(outliers, LMK_plotoutliers(sym.ls[[i]], plotALL = F, curves = SOF$curves, surfaces = SOF$surf))
  
}

#save plots
dev.copy(png, 'outliers-fossils.png')
dev.off()

dev.copy(png, 'outliers-phk-F.png')
dev.off()

dev.copy(png, 'outliers-phk-J.png')
dev.off()

dev.copy(png, 'outliers-phk-M.png')
dev.off()

dev.copy(png, 'outliers-phu-F.png')
dev.off()


dev.copy(png, 'outliers-phu-J.png')
dev.off()


dev.copy(png, 'outliers-phu-M.png')
dev.off()

#####When implementing through a for loop, deparse(substitute(x)) fails####
###Not sure if there's anything to be done about this, it'd work in any other cases

###Any better way to export these graphs?

##Reset names
names(outliers) <- paste(rep(names(all.ls), each = 4), rep(c("summary.info", "ind.info", "proc.coords", "mshape"), 7), sep = ".")
names(outliers)


View(outliers$PHK_M.ind.info)
t = (outliers$fossil.ind.info)
for (i in 1:10){
  LMK_compare_two(outliers$fossil.mshape, outliers$fossil.proc.coords[,,as.numeric(outliers$fossil.ind.info[i,3])])
  a <-readline(menu(c("Yes","No"), title = "Continue"))
  if (a == 2){stop()}
}

LMK_compare_two(outliers$fossil.mshape, outliers$fossil.proc.coords[,,8])

###writeland doesn't work within a for loop
for (i in 1:7){
  LMK_writeland_nts(sym.ls[[i]], paste(names(sym.ls)[i],".sym", sep=""))
}

LMK_writeland_nts(sym.ls[[1]], "fos.sym")

LMK_writeland_nts(sym.ls[[2]], "phkF.sym")
LMK_writeland_nts(sym.ls[[3]], "phkJ.sym")
LMK_writeland_nts(sym.ls[[4]], "phkM.sym")

LMK_writeland_nts(sym.ls[[5]], "phuF.sym")
LMK_writeland_nts(sym.ls[[6]], "phuJ.sym")
LMK_writeland_nts(sym.ls[[7]], "phuM.sym")


class(sym.ls$fossil)
class(sym.ls$PHK_F)




plot(1:5,1:5, pch = 16, col = rainbow(5, 1, 1, start = .2, end = .21))
normscale(n)

normscale(n)
n
plot(1:20, 1:20, col = col, pch = 16)
View(col)
plot(p$x[,1:2], pch = 16, col = col)
range(n)
scale(1:20)
n1 = as.numeric(as.character(out1$ind.info[,2]))
n1 = out1$ind.info[,2]
n2 = as.numeric(as.character(out1$ind.info[,3]))

n = cbind(n1,n2)
n = as.numeric(as.character(out1$ind.info[order(as.numeric(as.character(out1$ind.info[,3]))),2]))


out1$ind.info[out1$ind.info[,3] == 1,2] > 1 
dim(p$x)
pairs(p$x[,1:3], pch = 16, col = LMK_colramp(as.factor(n)))

View(out1$ind.info)
summary(p)$importance[2,10]
p$sdev
View(p$rotation)
sum(abs(p$rotation[,1]))
plot(p)
plotTangentSpace(a$coords)
plot
plot(prcomp(two.d.array(a$coords))$x[,1:3])
plot(t$x)
plotTangentSpace(gpagen(sym.ls$PHK_F))
dim(sym.ls$PHK_F)
#####
#####Manually Inspect Outliers####
#####
###How to automate this?


###comparetwo does what it should. it doesn't need additional parameters. lmk_screen is the wrapper that will take specific parameters such as (gpa align = T, bscale = F, comp.all = T)

##need to figure out way to jump to specific index
##readline for index number and just make sure to display current slide


View(outliers$PHK_M.ind.info)  ##outliers are already gpa aligned

t.gpa = gpagen(sym.ls$PHK_J_,curves = sli, surfaces = surf)

fix.gpa = gpagen(sym.ls$PHK_J_)

LMK_compare_two(sli.gpa$consensus, sli.gpa$coords[,,1])
LMK_compare_two(sli.gpa$consensus, fix.gpa$consensus)


writeClipboard(as.matrix(outliers$PHU_M.ind.info))

###might be helpful for screening####
cat("Press [1] to continue")
line <- readline()

a = readline(prompt = invisible(menu(c("A", "B", "C", "D"), title = "Choose one")))
a = readline(prompt = "press [enter] to continue")
######

#####
#####Combined Analysis####
#####

#####1A: Combine seperate arrays into one

all.dat = arrayspecs(
  rbind(
    two.d.array(fossil),
    two.d.array(PHK_F),
    two.d.array(PHK_J),
    two.d.array(PHK_M),
    
    two.d.array(PHU_F),
    two.d.array(PHU_J),
    two.d.array(PHU_M)
  ),p = 465, k = 3
)


#####1B: Read then all into one array#####

{all.dat = NULL
for (i in 1:7){
  all.dat = rbind(all.dat, two.d.array(readland.nts(list.files(pattern = ".dta")[i])))
}
all.dat = arrayspecs(all.dat, p = 465, k = 3)
rm(i)
}

##Check if first specimen is labeled atlas and, if so, remove

dimnames(all.dat)[[3]][[1]]
all.dat = all.dat[,,-1]



sym.dat = symmetrize(fixLMtps(fixLMmirror(all.dat, SOF_pair))$out, SOF_pair)

t = plotOutliers(gpa$coords)



gpa = gpagen(sym.dat, curves = sli, surfaces = surf)

###Check for outliers
t1 = LMK_plotoutliers(gpa$coords, gpa = F)
t2 = LMK_plotoutliers(sym.dat)



##Get PC Scores ###Maybe there's a better way for this
ts = plotTangentSpace(gpa$coords, warpgrids = F)
plot(1:20, ts$pc.summary$importance[1,1:20])
sum(ts$pc.summary$importance[1,1:10])

###Ontogenetic sample affecting outlier inspection


###Screen by PCA
View(ts$pc.summary$importance)
LMK_PCA_plot(gpa$coords, grp = class)
t <- LMK_screen(gpa$coords)
t <-LMK_plotoutliers(gpa$coords)
View(t[[2]])

points(ts$pc.scores[c(60,131,132,58),1:2], col = "black", pch = 16)
text(ts$pc.scores[c(60,131,132,58),1:2], labels = dimnames(gpa$coords)[[3]][c(60,131,132,58)], adj = c(.5,1))
View(ts$pc.scores)
dimnames(gpa$coords[,,60])

plot3d()
rm(screen)

library(lmkCHK)  
