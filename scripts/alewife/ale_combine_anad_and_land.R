#### Generate text file in two column format for selected anadromous populations and landlocked populations ####
library(tidyr)

# Read in genetic data
start <- read.csv('~/Documents/UCSC_postdoc/river_herring/data/alewife/ALE_origins_firstpick.csv', header = TRUE) #starting dataset that I'll add additional samples to
add <- read.csv('~/Documents/UCSC_postdoc/river_herring/data/alewife/ALE_ANAD_baseline_92_5517.csv', header = TRUE) #baseline alewife dataset that I will select other populations from

# Create column loci names for later
col_names <- read.csv('~/Documents/UCSC_postdoc/river_herring/data/alewife/ALE_ANAD_baseline_92_5517.csv', header = FALSE)
col_names <- data.frame(col_names[1,])
colnames(col_names) <- colnames(add)

# Check to make sure columns are in same order
colnames(start) == colnames(add) # yes

#### Subset the anadromous dataset to those populations that I've selected ####
# First read in which populations I want
anadromous_latlon <- read.delim("~/Documents/UCSC_postdoc/river_herring/data/alewife/ALE_anadromous_latlon.txt", header = TRUE) # Many anadromous populations from Reid et al 2018, but not all
landlocked_latlon <- read.delim("~/Documents/UCSC_postdoc/river_herring/data/alewife/ALE_landlocked_latlon.txt", header = TRUE) # Many anadromous populations from Reid et al 2018, but not all

pops <- unique(c(anadromous_latlon$Location, landlocked_latlon$Code))

# Next, separate sample name by letters and numbers so that I can subset the anadromous dataframe
add.samplenames <- data.frame(Sample = add$Code_name)

names <- add.samplenames %>%
  separate(Sample, 
           into = c("pop", "individual"), 
           sep = "_"
           # sep = "(?<=[A-Za-z])(?=[0-9])"
  )

# Finally, subset the anadromous dataframe based on population identifier by subsetting the names to only those in the populations I'm interested in. Then subset anadromous dataframe using the names
names.subset <- which(names$pop %in% pops) # indices of the populations I want

add.subset <- add[names.subset,] #1906 x 185, but this has some duplicated population
# Remove underscore from names in this data frame
add.subset.names <- add$Code_name[names.subset]
add.subset.names.split <- do.call(rbind, strsplit(as.character(add$Code_name[names.subset]), '_'))
add.subset.new.names <- paste(add.subset.names.split[,1], add.subset.names.split[,2], sep = '') # fish names without _
add.subset$Code_name <- add.subset.new.names # replace names with underscore with names that don't have underscore

#### Append first dataset to additional data ####
colnames(add.subset) == colnames(start) # double check column names are the same

ALE_land_origins <- rbind(add.subset, start) #3192 x 185, use merge instead of rbind to take care of duplicated populations or rbind then remove duplicated rows

#### Remove duplicate rows based on fish identifier
ALE_land_origins2 <- ALE_land_origins[!duplicated(ALE_land_origins$Code_name),] #2408 x 185
ALE_land_origins3 <- rbind(col_names, ALE_land_origins2) #add correct column names

write.table(ALE_land_origins3, "~/Documents/UCSC_postdoc/river_herring/data/alewife/ALE_land_origins_2408_92.txt", col.names = FALSE, row.names = FALSE, sep = '\t') # then remove all the ""'s and order the populations as I see fit
