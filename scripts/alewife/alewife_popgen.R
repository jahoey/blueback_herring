library(ade4)
library(adegenet)
library(hierfstat)
library(poppr)
library(ape)
library(ggtree)
library(ggplot2)
library(viridis)

#### Reading in data ####
# Reading in genepop file for 92 loci across 2408 landlocked and anadromous alewife
ale <- read.genepop("~/Documents/UCSC_postdoc/river_herring/data/alewife/ALE_land_origins_2408_92_Genepop.gen")

#### PCA of all anadromous and landlocked populations ####
sum(is.na(ale$tab)) #8676
X <- scaleGen(ale, NA.method = "mean")
dim(X)
class (X)

# make PCA
pca1 <- dudi.pca(X,cent=FALSE,scale=FALSE,scannf=FALSE,nf=3)
barplot(pca1$eig[1:50],main="PCA eigenvalues", col=heat.colors(50))
eig_percent <- round((pca1$eig/(sum(pca1$eig)))*100,2)
eig_percent [1:3]

# Plotting PC1 and PC2
s.label(pca1$li)
title("PCA of all alewife\naxes 1-2")
add.scatter.eig(pca1$eig[1:20], 3,1,2)

s.class(pca1$li, pop(ale)) 
title("PCA of all alewife\naxes 1-2")
add.scatter.eig(pca1$eig[1:20], 3,1,2)

png(file="~/Documents/UCSC_postdoc/river_herring/results/alewife/pca_all_pops.png", width=8, height=7, res=300, units="in")

par(
  mar=c(4, 3, 3, 1), # panel margin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=14, # point size, which is the font size
  bg=NA,
  bty = 'n'
)

s.class(pca1$li, pop(ale), xax=1,yax=2, col = viridis(43, alpha = 1), axesell=TRUE, cellipse=1.5, cstar=1,cpoint=1.75, grid=FALSE, addaxes = FALSE, xlim = c(-11,15), ylim = c(-15,14), clabel = 0)
axis(1, at=seq(-10,15, by=5), labels=seq(-10,15, by= 5), line = -1)
axis(2, at=seq(-10,10, by = 5), labels=seq(-10,10, by= 5), line = -3, las = 2)
mtext("PC1 (9.97%)", side = 1, line = 1.5)
mtext("PC2 (3.77%)", side = 2, line = -0.5)

legend(15, 12,
       legend=unique(pop(ale)),
       pch=19,
       col = viridis(43),
       bty = "n",
       y.intersp = 0.7,
       cex = 0.7)

dev.off()

1907:1954 #EGL
1955:2002 #LCH
2003:2013 # OSL
2014:2058 #CYL
2059:2081 #SNL
2082:2117 #CAL
2118:2156 #LON
2157:2211 #LMI
2212:2215 #LSU
2216:2265 #PAL
2266:2315 #ROG
2316:2362 #QUO
2363:2408 #KER

#### PCA for all landlocked and SE Atlantic anadromous populations ####
# Subset the original genind object to landlocked ALE only
ale_sub <- as.genind(ale@tab[1907:2408,]) #502 individuals, 92 loci
pop(ale_sub) <- factor(pop(ale)[1907:2408])

# PCA
sum(is.na(ale_sub$tab)) #2120
Y <- scaleGen(ale_sub, NA.method = "mean") #landlocked ALE
dim(Y)
class (Y)

pca2 <- dudi.pca(Y,cent=FALSE,scale=FALSE,scannf=FALSE,nf=3) #landlocked ALE
barplot(pca2$eig[1:50],main="PCA eigenvalues", col=heat.colors(50))
eig_percent <- round((pca2$eig/(sum(pca2$eig)))*100,2)
eig_percent [1:3]

png(file="~/Documents/UCSC_postdoc/river_herring/results/alewife/pca_land.png", width=8, height=7, res=300, units="in")

par(
  mar=c(4, 3, 3, 1), # panel magin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=14, # point size, which is the font size
  bg='white',
  bty = 'n'
)

s.class(pca2$li, pop(ale_sub), xax=1,yax=2, col = viridis(13, alpha = 0.5), axesell=TRUE, cellipse=1.5, cstar=1,cpoint=1.75, grid=FALSE, addaxes = FALSE, xlim = c(-10,6), ylim = c(-20,14), clabel = 0)
axis(1, at=seq(-10,10, by=5), labels=seq(-10,10, by= 5), line = -1.5)
axis(2, at=seq(-15,10, by = 5), labels=seq(-15,10, by= 5), line = -7, las = 2)
mtext("PC1 (20.2%)", side = 1, line = 1.5)
mtext("PC2 (12.3%)", side = 2, line = -4)
legend(10, 11,
       legend=unique(pop(ale_sub)),
       pch=19,
       col = viridis(13, alpha = 0.5),
       bty = "n",
       y.intersp = 0.9,
       cex = 1)

dev.off()

