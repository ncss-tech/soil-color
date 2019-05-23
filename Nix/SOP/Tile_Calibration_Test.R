# load required libraries
library(colorspace)
library(aqp)
library(farver)

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

# graphical comparison
x.hex <- hex(x.srgb)
par(bg='black')
swatchplot(x.hex)

# CIE delta-E00
d <- compare_colour(coords(x.lab)[, 1:3], from_space='LAB', method='CIE2000', white_from='D65')
d[d == 0] <- NA
par(bg='white')
hist(d, breaks = 10)

# ... these are all very small differences, below the level of perception

# convert to closest Munsell chip
# the @coords extracts just the color coordinates as a matrix
x.munsell <- rgb2munsell(x.srgb@coords)

# combine with source data
# row-order and length must be the same
y <- cbind(x, x.munsell)

y

