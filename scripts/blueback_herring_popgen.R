library(ade4)
library(adegenet)
library(hierfstat)
library(poppr)
library(ape)
library(ggtree)
library(viridis)

#### Reading in data ####
# Reading in genepop file for 95 loci across 405 landlocked blueback herring
bbh <- read.genepop("~/Documents/UCSC_postdoc/blueback_herring/data/BBH_genepop_LLMS_modifiedlandonly.gen")

# Reading in genepop file for 95 loci/190 alleles across 1228 anadromous and landlocked blueback herring (28 populations)
# bbh_all <- read.genepop("~/Documents/UCSC_postdoc/blueback_herring/data/BBH_genepop_LLMS_modified.gen")
bbh_all <- read.genepop("~/Documents/UCSC_postdoc/blueback_herring/data/BBH_land_origins_1970_95_genepop.gen") # additional anadromous populations

# PCA
sum(is.na(bbh$tab)) #1159
X <- scaleGen(bbh, NA.method = "mean") #landlocked only
dim(X)
class (X)

sum(is.na(bbh_all$tab)) #3094; 4570
Y <- scaleGen(bbh_all, NA.method = "mean") #landlocked and anadromous
dim(Y)
class (Y)

# make PCA
pca1 <- dudi.pca(X,cent=FALSE,scale=FALSE,scannf=FALSE,nf=3) #landlocked
barplot(pca1$eig[1:50],main="PCA eigenvalues", col=heat.colors(50))
eig_percent <- round((pca1$eig/(sum(pca1$eig)))*100,2)
eig_percent [1:3]

pca2 <- dudi.pca(Y,cent=FALSE,scale=FALSE,scannf=FALSE,nf=3) #landlocked and anadromous
barplot(pca2$eig[1:50],main="PCA eigenvalues", col=heat.colors(50))
eig_percent <- round((pca2$eig/(sum(pca2$eig)))*100,2)
eig_percent [1:3]

# Plotting PC1 and PC2
s.label(pca1$li) #landlocked
title("PCA of landlocked blueback herring\naxes 1-2")
add.scatter.eig(pca1$eig[1:20], 3,1,2)

s.class(pca1$li, pop(bbh)) 
title("PCA of landlocked blueback herring\naxes 1-2")
add.scatter.eig(pca1$eig[1:20], 3,1,2)

s.label(pca2$li) #landlocked and anadromous
title("PCA of all blueback herring\naxes 1-2")
add.scatter.eig(pca2$eig[1:20], 3,1,2)

s.class(pca2$li, pop(bbh_all))
title("PCA of all blueback herring\naxes 1-2")
add.scatter.eig(pca2$eig[1:20], 3,1,2)

### To make a nice plot of the PCA broken down by population ###
# Landlocked only by population
png(file="~/Documents/UCSC_postdoc/blueback_herring/results/pca_landlocked_pops.png", width=8, height=7, res=300, units="in")

par(
  mar=c(4, 3, 3, 1), # panel margin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=14, # point size, which is the font size
  bg=NA,
  bty = 'n'
)

s.class(pca1$li, pop(bbh), xax=1,yax=2, col = viridis(9, alpha = 1), axesell=TRUE, cellipse=1.5, cstar=1,cpoint=1.75, grid=FALSE, addaxes = FALSE, xlim = c(-11,11), ylim = c(-15,14), clabel = 0)
axis(1, at=seq(-10,10, by=5), labels=seq(-10,10, by= 5), line = -1)
axis(2, at=seq(-10,10, by = 5), labels=seq(-10,10, by= 5), line = -3, las = 2)
mtext("PC1 (8.76%)", side = 1, line = 1.5)
mtext("PC2 (3.66%)", side = 2, line = -0.5)

legend(7, 10,
       legend=c("LLA", "LHA", "LCT", "LNO", "LSE", "LYO", "LTU", "LBU", "LRA"),
       pch=19,
       col = viridis(9),
       bty = "n",
       y.intersp = 0.8,
       cex = 1)

# 1-46,47-92,93-139, 140-174,175-220,221-265,266-311,312-358,359-405

dev.off()

# All landlocked and anadromous by population
png(file="~/Documents/UCSC_postdoc/blueback_herring/results/pca_all_pops.png", width=8, height=7, res=300, units="in")

par(
  mar=c(4, 3, 3, 1), # panel magin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=14, # point size, which is the font size
  bg=NA,
  bty = 'n'
)

s.class(pca2$li, pop(bbh_all), xax=1,yax=2, col = plasma(35, alpha = 1), axesell=TRUE, cellipse=1.5, cstar=1,cpoint=1.75, grid=FALSE, addaxes = FALSE, xlim = c(-11,11), ylim = c(-15,14), clabel = 0)
axis(1, at=seq(-10,10, by=5), labels=seq(-10,10, by= 5), line = -1)
axis(2, at=seq(-10,10, by = 5), labels=seq(-10,10, by= 5), line = -3.5, las = 2)
mtext("PC1 (7.33%)", side = 1, line = 1.5)
mtext("PC2 (2.94%)", side = 2, line = -1)

legend(10, 12,
       legend=c("MAR", "PET", "EMA", "LOC","OYS", "EXE", "PAR", "MYS","HER", "MON", "GIL", "MET","DEL", "SUS", "POT", "RAP","YOR", "JAM", "CHO", "ROA", "NEU", "CF","SAN", "SAV", "ALT", "STR", "LLA", "LHA", "LCT", "LNO", "LSE", "LYO", "LTU", "LBU", "LRA"),
       pch=19,
       col = plasma(35),
       bty = "n",
       y.intersp = 0.7,
       cex = 0.85)

dev.off()

# All landlocked and anadromous by life history
png(file="~/Documents/UCSC_postdoc/blueback_herring/results/pca_all_life_hist.png", width=8, height=7, res=300, units="in")

par(
  mar=c(4, 3, 3, 1), # panel magin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=14, # point size, which is the font size
  bg=NA,
  bty = 'n'
)

# Create vector with life history code: 1 = anadromous, 2 = landlocked
# life_hist <- rep(1,1228) # 1 = anadromous (823) for smaller bbh_all dataset
# life_hist[776:1180] <- 2 # 2 = landlocked (405) for smaller bbh_all dataset

life_hist <- rep(1,1970) # 1 = anadromous (1565)
life_hist[1566:1970] <- 2 # 2 = landlocked (405)

pop(bbh_all) <- life_hist

s.class(pca2$li, pop(bbh_all), xax=1,yax=2, col = viridis(2, alpha = 0.5), axesell=TRUE, cellipse=1.5, cstar=1,cpoint=1.75, grid=FALSE, addaxes = FALSE, xlim = c(-11,11), ylim = c(-15,14), clabel = 0)
axis(1, at=seq(-10,10, by=5), labels=seq(-10,10, by= 5), line = -1)
axis(2, at=seq(-10,10, by = 5), labels=seq(-10,10, by= 5), line = -3.5, las = 2)
mtext("PC1 (7.33%)", side = 1, line = 1.5)
mtext("PC2 (2.94%)", side = 2, line = -1)
legend(7, 10,
       legend=c("anadromous", "landlocked"),
       pch=19,
       col = viridis(2),
       bty = "n",
       y.intersp = 1,
       cex = 1)

dev.off()

# All landlocked and anadromous by Reid et al 2018 regions 
png(file="~/Documents/UCSC_postdoc/blueback_herring/results/pca_all_region.png", width=8, height=7, res=300, units="in")

par(
  mar=c(4, 3, 3, 1), # panel magin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=14, # point size, which is the font size
  bg=NA,
  bty = 'n'
)

# Create vector with region code: 1 = NE, 2 = Mid-Atlantic, 3 = SE, 4 = landlocked
region <- rep(1,1970) # 1 = canada-nne (165)
region[166:312] <- 2 # 2 = MNE (147)
region[313:491] <- 3 # 3 = SNE (179)
region[492:1060] <- 4 # 4 = MAT (569)
region[1061:1565] <- 5 # 5 = SAT (505)
region[1566:1970] <- 6 # 6 = landlocked (405)

pop(bbh_all) <- region

s.class(pca2$li, pop(bbh_all), xax=1,yax=2, col = viridis(6, alpha = 0.5), axesell=TRUE, cellipse=1.5, cstar=1,cpoint=1.75, grid=FALSE, addaxes = FALSE, xlim = c(-11,11), ylim = c(-15,14), clabel = 0)
axis(1, at=seq(-10,10, by=5), labels=seq(-10,10, by= 5), line = -1)
axis(2, at=seq(-10,10, by = 5), labels=seq(-10,10, by= 5), line = -3.5, las = 2)
mtext("PC1 (7.33%)", side = 1, line = 1.5)
mtext("PC2 (2.94%)", side = 2, line = -1)
legend(5, -5,
       legend=c("anadromous - canada/nne","anadromous - mne","anadromous - sne", "anadromous - mat", "anadromous - sat", "landlocked"),
       pch=19,
       col = viridis(6),
       bty = "n",
       y.intersp = 0.8,
       cex = 0.9)
dev.off()

# PCA for all landlocked and SE Atlantic anadromous populations
# Subset the original genind object to landlocked and SE
bbh_all_sub <- as.genind(bbh_all@tab[1061:1970,])
pop(bbh_all_sub) <- pop(bbh_all)[1061:1970]

# PCA
sum(is.na(bbh_all_sub$tab)) #2122
Z <- scaleGen(bbh_all_sub, NA.method = "mean") #landlocked and SE anadromous
dim(Z)
class (Z)

pca3 <- dudi.pca(Z,cent=FALSE,scale=FALSE,scannf=FALSE,nf=3) #landlocked and SE anadromous
barplot(pca3$eig[1:50],main="PCA eigenvalues", col=heat.colors(50))
eig_percent <- round((pca3$eig/(sum(pca3$eig)))*100,2)
eig_percent [1:3]

png(file="~/Documents/UCSC_postdoc/blueback_herring/results/pca_land_and_se.png", width=8, height=7, res=300, units="in")

par(
  mar=c(4, 3, 3, 1), # panel magin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=14, # point size, which is the font size
  bg=NA,
  bty = 'n'
)

s.class(pca3$li, pop(bbh_all_sub), xax=1,yax=2, col = viridis(14, alpha = 0.5), axesell=TRUE, cellipse=1.5, cstar=1,cpoint=1.75, grid=FALSE, addaxes = FALSE, xlim = c(-11,11), ylim = c(-15,14), clabel = 0)
axis(1, at=seq(-10,10, by=5), labels=seq(-10,10, by= 5), line = -1)
axis(2, at=seq(-10,10, by = 5), labels=seq(-10,10, by= 5), line = -3.5, las = 2)
mtext("PC1 (6.20%)", side = 1, line = 1.5)
mtext("PC2 (3.25%)", side = 2, line = -1)
legend(10, 11,
       legend=c("CF","SAN", "SAV", "ALT", "STR", "LLA", "LHA", "LCT", "LNO", "LSE", "LYO", "LTU", "LBU", "LRA"),
       pch=19,
       col = viridis(14),
       bty = "n",
       y.intersp = 1,
       cex = 1)

dev.off()


#### FST calculations with each population representing a unique population ####
bbh_all_hier <- genind2hierfstat(bbh_all) # converts genind object into heirfstat object (each row is an individual, each column is a locus with cells denoting the alleles an individual has at each locus)
pairwise.WCfst(bbh_all_hier,diploid=TRUE) # calculates W-C Fst
fst <- genet.dist(bbh_all_hier, method = 'WC84')

write.table(as.matrix(round(fst, 4)), '~/Documents/UCSC_postdoc/blueback_herring/results/fst_table.txt')

#### Basic popgen stats ####
# Find columns that have zero variance (all the same values = monomorphic in this sense)
bbh_all_df <- cbind(data.frame(bbh_all@pop), bbh_all@tab) # make data frame with population & allele counts

var.table <- aggregate(bbh_all_df[, -1], list(bbh_all_df$bbh_all.pop), function(x) var(x,na.rm=T)==0)

table(as.character(var.table[1,-1])) # gives counts of number of polymorphic (FALSE), monomorphic (TRUE) or NA (check if all individuals are NA or not); each row is a population
var.table[1,-1]

# Heterozygosity
# Option 1
stats <- basic.stats(bbh_all_hier) # hierfstat package

bbh_all_split <- split(bbh_all_hier, f = bbh_all_hier$pop) #split out dataframe based on population
basic.stats(bbh_all_split[[1]][,-1]) #list number refers to population; doesn't work for all populations

# Option 2: Divide up genind object by population using poppr package
MAR <- popsub(bbh_all, sublist = 'MAR041')
MAR.he <- summary(MAR)
mean(MAR.he$Hobs)

PET <- popsub(bbh_all, sublist = 'PET032')
PET.he <- summary(PET)
mean(PET.he$Hobs)

EMA <- popsub(bbh_all, sublist = 'EMA046')
EMA.he <- summary(EMA)
mean(EMA.he$Hobs)

LOC <- popsub(bbh_all, sublist = 'LOC046')
LOC.he <- summary(LOC)
mean(LOC.he$Hobs)

OYS <- popsub(bbh_all, sublist = 'OYS034')
OYS.he <- summary(OYS)
mean(OYS.he$Hobs)

EXE <- popsub(bbh_all, sublist = 'EXE047')
EXE.he <- summary(EXE)
mean(EXE.he$Hobs)

PAR <- popsub(bbh_all, sublist = 'PAR046')
PAR.he <- summary(PAR)
mean(PAR.he$Hobs)

MYS <- popsub(bbh_all, sublist = 'MYS039')
MYS.he <- summary(MYS)
mean(MYS.he$Hobs)

HER <- popsub(bbh_all, sublist = 'HER047')
HER.he <- summary(HER)
mean(HER.he$Hobs)

MON <- popsub(bbh_all, sublist = 'MON047')
MON.he <- summary(MON)
mean(MON.he$Hobs)

GIL <- popsub(bbh_all, sublist = 'GIL046')
GIL.he <- summary(GIL)
mean(GIL.he$Hobs)

MET <- popsub(bbh_all, sublist = 'MET034')
MET.he <- summary(MET)
mean(MET.he$Hobs)

DEL <- popsub(bbh_all, sublist = 'DEL047')
DEL.he <- summary(DEL)
mean(DEL.he$Hobs)

SUS <- popsub(bbh_all, sublist = 'SUS038')
SUS.he <- summary(SUS)
mean(SUS.he$Hobs)

POT <- popsub(bbh_all, sublist = 'POT051')
POT.he <- summary(POT)
mean(POT.he$Hobs)

CPF <- popsub(bbh_all, sublist = 'CPF045')
CPF.he <- summary(CPF)
mean(CPF.he$Hobs)

SAN <- popsub(bbh_all, sublist = 'SAN042')
SAN.he <- summary(SAN)
mean(SAN.he$Hobs)

SAV <- popsub(bbh_all, sublist = 'SAV048')
SAV.he <- summary(SAV)
mean(SAV.he$Hobs)

LLA <- popsub(bbh_all, sublist = 'LLA046')
LLA.he <- summary(LLA)
mean(LLA.he$Hobs)

LHA <- popsub(bbh_all, sublist = 'LHA046')
LHA.he <- summary(LHA)
mean(LHA.he$Hobs)

LCT <- popsub(bbh_all, sublist = 'LCT047')
LCT.he <- summary(LCT)
mean(LCT.he$Hobs)

LNO <- popsub(bbh_all, sublist = 'LNO035')
LNO.he <- summary(LNO)
mean(LNO.he$Hobs)

LSE <- popsub(bbh_all, sublist = 'LSE046')
LSE.he <- summary(LSE)
mean(LSE.he$Hobs)

LYO <- popsub(bbh_all, sublist = 'LYO045')
LYO.he <- summary(LYO)
mean(LYO.he$Hobs)

LTU <- popsub(bbh_all, sublist = 'LTU046')
LTU.he <- summary(LTU)
mean(LTU.he$Hobs)

LBU <- popsub(bbh_all, sublist = 'LBU047')
LBU.he <- summary(LBU)
mean(LBU.he$Hobs)

LRA <- popsub(bbh_all, sublist = 'LRA047')
LRA.he <- summary(LRA)
mean(LRA.he$Hobs)

STJ <- popsub(bbh_all, sublist = 'STJ048')
STJ.he <- summary(STJ)
mean(STJ.he$Hobs)

#### Generating and plotting a neighbor joining tree ####
bbh_dist <- genet.dist(bbh_all_hier, method = 'Ds')

temp <- as.data.frame(as.matrix(bbh_dist))
table.paint(temp, cleg=0, clabel.row = 0.5, clabel.col = 0.5) #darker shades mean greater distances

tree <- nj(bbh_dist) #unrooted
tree$edge.length[27] <- abs(tree$edge.length[27]) # some branches are negative because of nj method, so options are to leave them, set them to zero, or take absolute value
tree$edge.length[31] <- abs(tree$edge.length[31])
tree$tip.label <- c('MAR', 'PET', 'EMA', 'LOC', 'OYS', 'EXE', 'PAR', 'MYS', 'HER', 'MON', 'GIL', 'MET', 'DEL', 'SUS', 'POT', 'CPF', 'SAN', 'SAV', 'LLA', 'LHA','LCT', 'LNO', 'LSE', 'LYO', 'LTU', 'LBU', 'LRA', 'STJ') # modify tip labels for easier plotting

plot(tree, show.tip = FALSE, cex = 0.8)
tiplabels(unique(bbh_all@pop), cex= 0.6, bg = 'white', frame = 'none', adj = -0.1)
plot(tree,type='unrooted',show.tip=FALSE)
title('Unrooted NJ tree')
tiplabels(unique(bbh_all@pop), cex= 0.6, bg = 'white', frame = 'none')

# tree.root <- root(tree, out = 1) #rooted, but what is the root?
# plot(tree.root, show.tip = FALSE, cex = 0.8)
# tiplabels(unique(bbh_all@pop), cex= 0.6, bg = 'white', frame = 'none', adj = -0.1)

# Two trees
# Circular cladogram tree
png(file="~/Documents/UCSC_postdoc/blueback_herring/results/bbh_all_tree_circle.png", width=8, height=7, res=300, units="in")

par(
  mar=c(4, 3, 3, 1), # panel magin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=14, # point size, which is the font size
  bg=NA,
  bty = 'n'
)

ggtree(tree, layout = 'circular', branch.length = 'none') + geom_treescale(x = 7, y =0) + geom_tiplab2(size = 4, hjust = -0.1, aes(angle=angle))

dev.off()

# Unrooted fan cladogram tree
png(file="~/Documents/UCSC_postdoc/blueback_herring/results/bbh_all_tree_fan.png", width=8, height=7, res=300, units="in")

par(
  mar=c(4, 3, 3, 1), # panel magin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=14, # point size, which is the font size
  bg=NA,
  bty = 'n'
)

ggtree(tree, layout = 'daylight', branch.length = 'none') + geom_treescale(x=2,y=-25) + geom_tiplab2(size = 3, hjust = -0.1, aes(angle=angle)) + ggplot2::xlim(2,26) + ggplot2::ylim(-25,0)

dev.off()

# Unrooted daylight
png(file="~/Documents/UCSC_postdoc/blueback_herring/results/bbh_all_tree_unrooted_phylo_fan.png", width=8, height=7, res=300, units="in")

par(
  mar=c(4, 3, 3, 1), # panel magin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=14, # point size, which is the font size
  bg=NA,
  bty = 'n'
)

ggtree(tree, layout = 'daylight') + geom_tiplab2(size = 2.5, hjust = -0.1, aes(angle=angle)) + ggplot2::xlim(-0.06,0.05) + ggplot2::ylim(-0.05,0.03) + geom_treescale(x=0.0,y=0.01, width = 0.01) + annotate("text", label = "0.01", x = 0.005, y = 0.013) 

dev.off()

#### Exploring how genetic diversity scales with dam age and size ####
data <- read.table("~/Documents/UCSC_postdoc/blueback_herring/results/bbh_all_%poly.txt", header = TRUE)

plot(data$Dam_year,data$X._polymorphic)
plot(data$Dam_year,data$H_exp)

plot(data$Dam_size,data$X._polymorphic)
plot(data$Dam_size,data$H_exp)
