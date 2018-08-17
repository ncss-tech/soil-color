library(soilDB)
library(aqp)


# get all CA630 pedons
x <- fetchNASIS(rmHzErrors = FALSE)

## Andrew already did the work to associate pedon records with intersecting map units
# load some suggested peiids from Andrew's files
p.1 <- readLines('S:/NRCS/Archive_Andrew_Brown/Scripts/GroupedProfilePlotByMU/pedonlists/6202.csv')
p.2 <- readLines('S:/NRCS/Archive_Andrew_Brown/Scripts/GroupedProfilePlotByMU/pedonlists/6070.csv')

# split text string into vector of characters
p.1 <- strsplit(p.1, ',')[[1]]
p.2 <- strsplit(p.2, ',')[[1]]

# combine into single character vector
p <- unique(c(p.1, p.2))

# make an index to those profiles
idx <- which(profile_id(x) %in% p)

# check: looks OK
plot(x[idx, ], label='pedon_id')

# subset
x <- x[idx, ]

# save SPC for later use
save(x, file='pedon-data.rda')

# combine site + horizon data
d <- as(x, 'data.frame')

# subset columns
labsheet <- d[, c('peiid', 'phiid', 'site_id', 'pedon_id', 'hzdept', 'hzdepb', 'hzname')]

# generate a serial number: this is specific to a "lab sheet"
labsheet$serial <- 1:nrow(labsheet)

# add columns for LAB coordinates
labsheet$L <- NA
labsheet$A <- NA
labsheet$B <- NA

# each batch should have its own lab sheet and associated "folder" in the Nix Sensor app
# save
write.csv(labsheet, file='soil-color-labsheet-2017-07-19.csv', row.names = FALSE, na = '')

