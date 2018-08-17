# load required libraries
library(colorspace)
library(aqp)
library(plyr)

# load CSV files  
x <- read.csv('Tile_Calibration_Test.csv', stringsAsFactors = FALSE)

# check the results
str(x)
head(x)

# create LAB object
x.lab <- LAB(x$L, x$A, x$B, names=x$Color.Name)

# check graphically
plot(x.lab, cex=5)

# convert to sRGB
x.srgb <- as(x.lab, 'sRGB')
plot(x.srgb, cex=5)

# convert to closest Munsell chip
# the @coords extracts just the color coordinates as a matrix
x.munsell <- rgb2munsell(x.srgb@coords)

# combine with source data
# row-order and length must be the same
y <- cbind(x, x.munsell)

y

