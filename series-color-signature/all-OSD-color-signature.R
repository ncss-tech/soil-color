library(aqp)
library(soilDB)
library(sharpshootR)
library(plyr)
library(colorspace)

## why not use profile compare with pigments at some regular interval?
##
##

## !!! more efficient distance calculation
## -----> https://cran.r-project.org/web/packages/gower/vignettes/intro.html


# OSD colors
s <- read.csv('osd_colors.csv.gz', stringsAsFactors = FALSE)
# remove crap
s <- subset(s, subset = series != '' & !is.na(series))
# re-order just in case
s <- s[order(s$series, s$top), ]

# MLRA: multiple records / series
mlra <- read.csv('series-mlra-overlap.csv', stringsAsFactors = FALSE)

## TODO: moist / dry?
# manually convert Munsell -> RGB
rgb.data <- with(s, munsell2rgb(matrix_wet_color_hue, matrix_wet_color_value, matrix_wet_color_chroma, return_triplets = TRUE))
s$r <- rgb.data$r
s$g <- rgb.data$g
s$b <- rgb.data$b

# convert to SPC
depths(s) <- series ~ top + bottom

# ## TODO: appropriate rescaling of L?
# # this takes a while...
# # get signatures
# pig <- soilColorSignature(s, RescaleLightnessBy = 5)
# row.names(pig) <- pig$series
# pig$series <- NULL

## possibly best compromise
## ~ 15 minutes
## depth-slicing: NAs possible where there are missing data or hz errors
pig <- soilColorSignature(s, RescaleLightnessBy = 5, method='depthSlices')
row.names(pig) <- pig$series
pig$series <- NULL

# ## seems to be the best (?)
# ## much slower: ~ 3 hours
# ## PAM
# pig <- soilColorSignature(s, RescaleLightnessBy = 5, method='pam', pam.k = 3)
# row.names(pig) <- pig$series
# pig$series <- NULL
# 

# iterate over series
l.res <- list()
n <- nrow(pig)

## very fast: 5 minutes
## vectorized distance calc via Euclidean distance
##   -> only those series that occur in the same MLRA (could be several MLRA) are included in the distance calc
##   -> consider using only those MLRA with score > x

## TODO: optimal distance metric?
## https://cran.r-project.org/web/packages/gower/vignettes/intro.html

pb <- txtProgressBar(min=1, max=n, style=3)
for(i in 1:n) {
  p.i <- pig[i, ]
  name.i <- row.names(p.i)
  
  # filter out only those series that occur in the same MLRA as the current series
  which.mlra <- mlra$mlra[which(mlra$series == name.i)]
  
  ## note: there may be some series with no MLRA overlap information: skip these
  if(length(which.mlra) < 1)
    next
  
  mlra.idx <- which(mlra$mlra %in% which.mlra)
  mlra.filter <- unique(mlra$series[mlra.idx])
  idx <- which(row.names(pig) %in% mlra.filter)
  # apply filter
  pig.mlra <- pig[idx, ]
  
  # remove the current row
  current.idx <- which(row.names(pig.mlra) == name.i)
  pig.mlra <- pig.mlra[-current.idx, ]
  
  # compute pair-wise distances between current row and all others
  d <- sweep(pig.mlra, MARGIN = 2, STATS = as.matrix(p.i), FUN = '-')
  # NA allowed in data matrix... not sure how this will affect things
  d <- sqrt(rowSums(d^2, na.rm = TRUE))
  
  # assign top N to current series
  ## keep top 25
  top.n <- sort(d)[1:25]
  l.res[[name.i]] <- top.n
  setTxtProgressBar(pb, i)
}
close(pb)

## should scaling be done over all series or within logical strata?
# rescale distances by max distance within collection
max.dist <- max(sapply(l.res, max, na.rm=TRUE), na.rm=TRUE)
l.res <- lapply(l.res, function(i) i / max.dist)


## write out to file
save(l.res, pig, file='cached-pigment-data-depthSlices.Rda')

# distribution?
quantile(unlist(l.res), na.rm=TRUE)

## check closest 10
s.name <- toupper('amador')
s.list <- c(s.name, names(l.res[[s.name]])[1:10])
s.ranks <- c(-1, l.res[[s.name]][1:10])
rank.data <- data.frame(id=s.list, rank=s.ranks, stringsAsFactors = FALSE)

s.osd <- fetchOSD(s.list)
site(s.osd) <- rank.data

SoilTaxonomyDendrogram(s.osd)

par(mar=c(3,0,0,1))
new.order <- order(s.osd$rank)
plot(s.osd, plot.order=new.order)
axis(1, at=1:length(s.osd), labels = c(NA, round( s.osd$rank[new.order][-1], 3)) , line=-1, las=2, cex.axis=0.85)


## how can we better describe distance in the profile plot?
## make more room, and re-scale distances to integers from 1:n + fudge factor
##



## this won't work, no coordinates...
# plotTransect(s.osd, 'rank', plot.order=new.order)


## TODO: split by MLRA
# plot all OSD colors
lab.colors <- as(RGB(rgb.data[['r']], rgb.data[['g']], rgb.data[['b']]), 'LAB')@coords
cols <- cbind(rgb.data, lab.colors)
cols <- na.omit(cols)
cols <- as.data.frame(cols)
png(file='all-osd-colors-LAB.png', width=900, height=900, type='cairo', antialias = 'subpixel')
pairs(~ L + A + B, data=cols, pch=21, cex=2, col='white', bg=rgb(cols$r, cols$g, cols$b))
dev.off()


