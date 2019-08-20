###Prelim screening and analysis for Cranial SOF dataset. Fossil + Papionin
library(lmkCHK) #loads geomorph and morpho


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


#####Run through screening criteria individually first to spot outliers, then combine for actual analysis


#####
#####Preliminary Screening#####
#####

#####1: Read in each DTA indvidually#####


fossil <- readland.nts("Fossil_ZMPatch.dta")
fossil <- fossil[,,-1] ##Remove extant atlas

PHK_F <- readland.nts("PHK-F_zMPatch.dta")
PHK_J <- readland.nts("PHK-J_zMPatch.dta")
PHK_M <- readland.nts("PHK-M_zMPatch.dta")

PHU_F <- readland.nts("PHU-F_zMPatch.dta")
PHU_J <- readland.nts("PHU-J_zMPatch.dta")
PHU_M <- readland.nts("PHU-M_zMPatch.dta")



#####1.1: Restructure data#####
#restructure raw data if needed, create/read in classifiers

all.ls = list("fossil" = fossil, "PHK_F" = PHK_F, "PHK_J_" = PHK_J, "PHK_M" = PHK_M, "PHU_F" = PHU_F, "PHU_J" = PHU_J, "PHU_M" = PHU_M)



###set parsable names
#######NOTE: This step should occur now, but will be implemented later to get on with SOP#####

###Additional Classifiers

class = c(rep("fos", 25), rep("Kf", 20), rep("Kj", 20), rep("Km", 21), rep("Uf", 6), rep("Uj",31), rep("Um",21))

#####1.2: Define Landmarks####

##Paired LMs
SOF_pair = matrix(c(10:15, 241:465, 4:9, 16:240), nrow = 231, ncol = 2)

##Define Countour semi-landmarks
sli = psych::read.clipboard() ##slider table made by hand in excel

##Surface semi-landmarks

surf = 1:225 #all landmarks
surf = surf[!(1:225 %in% sli[,2])] #drop countours
surf = surf[-c(1,15,106,211, 225)] #drop fixed LMKs


#####





###2: Cleaning and evaluating seperately#####
sym.ls = NULL


for (i in 1:7){
  sym.ls[[i]] <- lmkCHK::LMK_sym(all.ls[[i]], SOF_pair)
}

outliers <- NULL

for ( i in 1:7){
  outliers <- append(outliers, LMK_plotoutliers(sym.ls[[i]], plotALL = F, curves = sli, surfaces = surf))
  
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


View(outliers$fossil.ind.info)

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


#####
#####Manually Inspect Outliers####
#####
###How to automate this?


###comparetwo does what it should. it doesn't need additional parameters. lmk_screen is the wrapper that will take specific parameters such as (gpa align = T, bscale = F, comp.all = T)

##need to figure out way to jump to specific index
##readline for index number and just make sure to display current slide


View(outliers$PHU_M.ind.info)  ##outliers are already gpa aligned

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
