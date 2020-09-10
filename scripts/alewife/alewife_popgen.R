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

png(file="~/Documents/UCSC_postdoc/river_herring/results/alewife/pca_regionalAND_uniqLAND.png", width=8, height=7, res=300, units="in")

par(
  mar=c(4, 3, 3, 1), # panel margin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=14, # point size, which is the font size
  bg=NA,
  bty = 'n'
)

# Create vector with region code: 1 = CAN, 2 = NNE, 3 = SNE, 4 = MAT, followed by landlocked pops
groups <- as.vector(pop(ale))

groups[c(1:150, 226:272)] <- 'CAN' # 1 = CAN (197)
groups[c(151:225, 273:872)] <- 'NNE' # 2 = NNE (675)
groups[873:1314] <- 'SNE' # 3 = SNE (442)
groups[1315:1906] <- 'MAT' # 4 = MAT (592)

pop(ale) <- groups

s.class(pca1$li[1:1906,], pop(ale)[1:1906], xax=1,yax=2, col = viridis(17, alpha = 0.6), axesell=TRUE, cellipse=1.5, cstar=1,cpoint=1.2, grid=FALSE, addaxes = FALSE, xlim = c(-11,19), ylim = c(-15,14), clabel = 0, pch = 17)
s.class(pca1$li[1907:2408,], pop(ale)[1907:2408], xax=1,yax=2, col = viridis(17, alpha = 0.6), add.plot = TRUE, clabel = 0, cpoint=1.75)
axis(1, at=seq(-5,15, by=5), labels=seq(-5,15, by= 5), line = -1)
axis(2, at=seq(-10,10, by = 5), labels=seq(-10,10, by= 5), line = -3, las = 2)
mtext("PC1 (9.97%)", side = 1, line = 1.5)
mtext("PC2 (3.77%)", side = 2, line = -0.5)

legend(15, 10,
       legend=unique(pop(ale))[1:4],
       title = 'Anadromous',
       pch=17,
       col = viridis(17, alpha = 0.8)[1:4],
       bty = "n",
       y.intersp = 0.8,
       cex = 0.9)

legend(15.3, 5.5,
       legend=unique(pop(ale))[5:17],
       title = 'Landlocked',
       pch=19,
       col = viridis(17, alpha = 0.8)[5:17],
       bty = "n",
       y.intersp = 0.8,
       cex = 0.9)

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

#### PCA for all landlocked populations ####
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

#### FST calculations with each population representing a unique population ####
ale_hier <- genind2hierfstat(ale) # converts genind object into heirfstat object (each row is an individual, each column is a locus with cells denoting the alleles an individual has at each locus)
pairwise.WCfst(ale_hier,diploid=TRUE) # calculates W-C Fst
fst <- genet.dist(ale_hier, method = 'WC84')

write.table(as.matrix(round(fst, 4)), '~/Documents/UCSC_postdoc/river_herring/results/alewife/fst_table.txt')

# Mean FST among all pops, landlocked only and anadromous only
mean_fst <- read.table('~/Documents/UCSC_postdoc/river_herring/results/alewife/fst_table.txt')
mean_fst[upper.tri(mean_fst)] <- NA
mean(data.matrix(mean_fst), na.rm = TRUE)

land <- mean_fst[31:43,31:43]
mean(data.matrix(land), na.rm = TRUE)

anad <- mean_fst[1:30,1:30]
mean(data.matrix(anad), na.rm = TRUE)

#### Basic popgen stats ####
# Find columns that have zero variance (all the same values = monomorphic in this sense)
ale_df <- cbind(data.frame(ale@pop), ale@tab) # make data frame with population & allele counts

var.table <- aggregate(ale_df[, -1], list(ale_df$ale.pop), function(x) var(x,na.rm=T)==0)

table(as.character(var.table[1,-1])) # gives counts of number of polymorphic (FALSE), monomorphic (TRUE) or NA (check if all individuals are NA or not); each row is a population
var.table[1,-1]

# Observed heterozygosity is calculated as part of mstoolkit, but can also do something like below
MIR <- popsub(ale, sublist = 'MIR')
MIR.he <- summary(MIR)
mean(MIR.he$Hobs)

PET <- popsub(ale, sublist = 'PET')
PET.he <- summary(PET)
mean(PET.he$Hobs)

#### Generating and plotting a neighbor joining tree ####
ale_dist <- genet.dist(ale_hier, method = 'Ds')

temp <- as.data.frame(as.matrix(ale_dist))
table.paint(temp, cleg=0, clabel.row = 0.5, clabel.col = 0.5) #darker shades mean greater distances

tree <- nj(ale_dist) #unrooted
tree$edge.length[12] <- abs(tree$edge.length[12]) # some branches are negative because of nj method, so options are to leave them, set them to zero, or take absolute value; easier to read when plotted as regular looking tree
tree$edge.length[25] <- abs(tree$edge.length[25])
tree$edge.length[49] <- abs(tree$edge.length[49])
tree$edge.length[51] <- abs(tree$edge.length[51])
tree$edge.length[53] <- abs(tree$edge.length[53])
tree$edge.length[54] <- abs(tree$edge.length[54])
tree$edge.length[55] <- abs(tree$edge.length[55])
tree$tip.label <- as.character(unique(ale_hier$pop))

plot(tree, show.tip = FALSE, cex = 0.8)
tiplabels(unique(ale@pop), cex= 0.6, bg = 'white', frame = 'none', adj = -0.1)
plot(tree,type='unrooted',show.tip=FALSE)
title('Unrooted NJ tree')
tiplabels(unique(ale@pop), cex= 0.6, bg = 'white', frame = 'none')

x <- as.vector(ale_dist)
y <- as.vector(as.dist(cophenetic(tree)))
plot(x, y, xlab="original pairwise distances", ylab="pairwise distances on the tree", main="Is NJ appropriate?", pch=20, col=transp("black",.1), cex=3)
abline(lm(y~x), col="red")
cor(x,y)^2

# Unrooted tree
png(file="~/Documents/UCSC_postdoc/river_herring/results/alewife/ale_tree_unrooted_phylo.png", width=8, height=7, res=300, units="in")

par(
  mar=c(4, 3, 3, 1), # panel magin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=14, # point size, which is the font size
  bg=NA,
  bty = 'n'
)

ggtree(tree) + geom_tiplab(size = 3, hjust = -0.15) + geom_treescale(x=0.05,y=0, width = 0.01)

dev.off()

# Unrooted daylight
png(file="~/Documents/UCSC_postdoc/river_herring/results/alewife/ale_tree_unrooted_phylo_fan.png", width=8, height=7, res=300, units="in")

par(
  mar=c(4, 3, 3, 1), # panel magin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=14, # point size, which is the font size
  bg=NA,
  bty = 'n'
)

ggtree(tree, layout = 'daylight') + geom_tiplab2(size = 2.3, hjust = -0.15, aes(angle=angle)) + ggplot2::xlim(-0.1,0.05) + ggplot2::ylim(-0.15,0.02) + geom_treescale(x=-0.025,y=-0.05, width = 0.01) + annotate("text", label = "0.01", x = -0.02, y = -0.06) 

dev.off()

# Circular cladogram tree
png(file="~/Documents/UCSC_postdoc/river_herring/results/alewife/ale_tree_circle.png", width=8, height=7, res=300, units="in")

par(
  mar=c(4, 3, 3, 1), # panel magin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=14, # point size, which is the font size
  bg=NA,
  bty = 'n'
)

ggtree(tree, layout = 'circular', branch.length = 'none') + geom_treescale(x = -2.5, y =0) + geom_tiplab2(size = 4, hjust = -0.1, aes(angle=angle))

dev.off()
