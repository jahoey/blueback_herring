#### Generate text file in two column format for selected anadromous populations and landlocked populations ####
library(tidyr)

# Read in genetic data
anadromous <- read.csv('~/Documents/UCSC_postdoc/river_herring/data/blueback_herring/BBH_ANAD_baseline_95_2201.csv', header = TRUE) #2201 fish, 95 biallelic snps
land <- read.table('~/Documents/UCSC_postdoc/river_herring/data/blueback_herring/BBH_landlocked_twocol.txt', header = TRUE) #405 fish, 95 biallelic snps

# Create column loci names for later
col_names <- read.table('~/Documents/UCSC_postdoc/river_herring/data/blueback_herring/BBH_landlocked_twocol.txt', header = FALSE)
col_names <- data.frame(col_names[1,])
colnames(col_names) <- colnames(land)

# Check to make sure columns are in same order
colnames(anadromous) == colnames(land) # yes

#### Subset the anadromous dataset to those populations that I've selected ####
# First read in which populations I want
anadromous_latlon <- read.delim("~/Documents/UCSC_postdoc/river_herring/maps/anadromous_latlon.txt", header = TRUE) # Many anadromous populations from Reid et al 2018, but not all
pops <- anadromous_latlon$Code

# Next, separate sample name by letters and numbers so that I can subset the anadromous dataframe
anadromous.samplenames <- data.frame(Sample = anadromous$Sample)

names <- anadromous.samplenames %>%
  separate(Sample, 
           into = c("pop", "individual"), 
           sep = "(?<=[A-Za-z])(?=[0-9])"
  )

# Finally, subset the anadromous dataframe based on population identifier by subsetting the names to only those in the populations I'm interested in. Then subset anadromous dataframe using the names
names.subset <- which(names$pop %in% pops) # indices of the populations I want

anadromous.subset <- anadromous[names.subset,] #1565 x 191

#### Append landlocked population to anadromous data ####
colnames(anadromous.subset) == colnames(land) # double check column names are the same

BBH_land_origins <- rbind(anadromous.subset, land) #1970 x 191
BBH_land_origins2 <- rbind(col_names, BBH_land_origins) #add correct column names

write.table(BBH_land_origins2, "~/Documents/UCSC_postdoc/river_herring/data/blueback_herring/BBH_land_origins_1970_95.txt", col.names = FALSE, row.names = FALSE, sep = '\t')
