library(aqp)
library(plyr)
library(cluster)
library(ape)
library(colorspace)
library(soilDB)
library(sharpshootR)

## TODO: how can we better describe change in color with depth?
##  --> slicing and relative position?


## TODO: integrate ideas from colordistance package
# https://cran.r-project.org/web/packages/colordistance/vignettes/colordistance-introduction.html



### trivial example, not very interesting
data(sp1)
sp1$soil_color <- with(sp1, munsell2rgb(hue, value, chroma))
depths(sp1) <- id ~ top + bottom

# manually convert Munsell -> RGB
rgb.data <- munsell2rgb(sp1$hue, sp1$value, sp1$chroma, return_triplets = TRUE)
sp1$r <- rgb.data$r
sp1$g <- rgb.data$g
sp1$b <- rgb.data$b

pig <- soilColorSignature(sp1)
row.names(pig) <- pig$id
d <- daisy(pig[, 2:6])
dd <- diana(d)

plotProfileDendrogram(sp1, dd, dend.y.scale = 0.25, scaling.factor = 0.001, y.offset = 0.02, width=0.15)
 

### use OSDs
s.list <- c('amador', 'redding', 'pentz', 'willows', 'pardee', 'yolo', 'hanford', 'cecil', 'sycamore', 'KLAMATH', 'MOGLIA', 'boomer', 'vleck', 'drummer', 'CANEYHEAD', 'musick', 'sierra', 'HAYNER', 'zook', 'argonaut', 'PALAU')
# s.list <- c('amador', 'redding', 'pentz', 'willows', 'pardee', 'yolo', 'hanford', 'cecil', 'sycamore')

s <- fetchOSD(s.list)

plot(s)

# manually convert Munsell -> sRGB
rgb.data <- munsell2rgb(s$hue, s$value, s$chroma, return_triplets = TRUE)
s$r <- rgb.data$r
s$g <- rgb.data$g
s$b <- rgb.data$b

## colorBucket method
# what is appropriate rescaling of L?
pig <- soilColorSignature(s, RescaleLightnessBy = 5, method = 'colorBucket')

# move row names over for distance matrix
row.names(pig) <- pig$id

d <- daisy(pig[, -1])
dd <- diana(d)

png(file='soils-by-pigments.png', width=1400, height=600, res=120, type='cairo', antialias = 'subpixel', pointsize = 14)
par(mar=c(1,1,1,1))
plotProfileDendrogram(s, dd, dend.y.scale = max(d) * 2, scaling.factor = max(d) / max(s), y.offset = max(d) / 10, width=0.15, cex.names=0.45)
dev.off()


## depthSlices method
# what is appropriate rescaling of L?
pig <- soilColorSignature(s, RescaleLightnessBy = 5, method='depthSlices')


# move row names over for distance matrix
row.names(pig) <- pig$id

d <- daisy(pig[, -1])
dd <- diana(d)

png(file='soils-by-depth-slice-pigments.png', width=1400, height=600, res=120, type='cairo', antialias = 'subpixel', pointsize = 14)
par(mar=c(1,1,1,1))
plotProfileDendrogram(s, dd, dend.y.scale = max(d) * 2, scaling.factor = 0.25, y.offset = 6, width=0.15, cex.names=0.45)
dev.off()



## PAM method
# what is appropriate rescaling of L?
pig <- soilColorSignature(s, RescaleLightnessBy = 5, method = 'pam', pam.k = 3)

# move row names over for distance matrix
row.names(pig) <- pig$id

d <- daisy(pig[, -1])
dd <- diana(d)

png(file='soils-by-PAM-pigments.png', width=1400, height=600, res=120, type='cairo', antialias = 'subpixel', pointsize = 14)
par(mar=c(1,1,1,1))
plotProfileDendrogram(s, dd, dend.y.scale = max(d) * 2, scaling.factor = 0.25, y.offset = 6, width=0.15, cex.names=0.45)
dev.off()







## tinkering
rgb.colors <- munsell2rgb(s$hue, s$value, s$chroma, return_triplets = TRUE)
lab.colors <- as(RGB(rgb.colors[['r']], rgb.colors[['g']], rgb.colors[['b']]), 'LAB')@coords
cols <- cbind(rgb.colors, lab.colors)
cols <- na.omit(cols)
cols <- as.data.frame(cols)
pairs(~ L + A + B, data=cols, pch=16, cex=2, col=rgb(cols$r, cols$g, cols$b))


## transition probabilities of colors
library(sharpshootR)
library(igraph)
library(plyr)
library(markovchain)


s$color <- paste0(s$hue, ' ', s$value, '/', s$chroma)
sort(table(s$color))

# remove missing colors
h <- horizons(s)
h <- subset(h, subset=h$hue != '')
horizons(s) <- h

tp <- hzTransitionProbabilities(s, "color", loopTerminalStates = FALSE)
par(mar = c(1, 1, 1, 1))
g <- plotSoilRelationGraph(tp, graph.mode = "directed", edge.arrow.size = 0.5, edge.scaling.factor = 2, vertex.label.cex = 0.75, vertex.label.family = "sans")


d <- data.frame(color=V(g)$name, stringsAsFactors = FALSE)
d <- join(d, horizons(s), by = 'color', type='left', match='first')
V(g)$color <- d$soil_color

# convert colors to HSV
hsv.cols <- t(rgb2hsv(col2rgb(V(g)$color)))
hsv.cols[, 1] <- 0
hsv.cols[, 2] <- 0
hsv.cols[, 3] <- ifelse(hsv.cols[, 3] > 0.5, 0, 1)
V(g)$label.color <- hsv(hsv.cols[, 1], hsv.cols[, 2], hsv.cols[, 3])

# remove loops only
g <- simplify(g, remove.loops = TRUE, remove.multiple=FALSE)

png(file='soil-color-transition-probabilities-1.png', width=1200, height=1200, res=90, type='cairo', antialias = 'subpixel', pointsize = 14)
par(mar=c(0,0,0,0), bg=grey(0.85))
set.seed(1010101)
plot(g, edge.arrow.size = 0.25, vertex.label.cex = 0.55, vertex.label.family = "sans", vertex.label.font=2, edge.color='black')
dev.off()

tp <- hzTransitionProbabilities(s, "color", loopTerminalStates = TRUE)
mc <- new("markovchain", states = dimnames(tp)[[1]], transitionMatrix = tp)

rmarkovchain(10, mc, include.t0 = TRUE, t0 = "10YR 3/1")
mostLikelyHzSequence(mc, "10YR 3/1")

# simulate colors: must use custom transition matrix
i <- 5
x <- s[i, ]
x$new.color <- rmarkovchain(nrow(x)-1, mc, include.t0 = TRUE, t0 = x$color[1])
x$new.color <- parseMunsell(x$new.color)

plotMultipleSPC(list(x, s[i, ]), group.labels=c('MC', 'Original'), color='new.color', name='hzname', max.d=150)

## old stuff


# 
# data("loafercreek")
# x <- loafercreek
# 
# x.moist.LAB <- as(with(horizons(x), RGB(m_r, m_g, m_b)), 'LAB')
# plot(x.moist.LAB)
# 
# x$L <- x.moist.LAB@coords[, 1]
# x$A <- x.moist.LAB@coords[, 2]
# x$B <- x.moist.LAB@coords[, 3]
# # 
# # get data fraction
# x$data.fraction <- evalMissingData(x, vars=c('L', 'A', 'B'))
# 
# 
# # another data set
# data(sp2)
# x <- sp2
# depths(x) <- id ~ top + bottom
# 
# x.moist.LAB <- as(with(horizons(x), RGB(r, g, b)), 'LAB')
# plot(x.moist.LAB)
# 
# x$L <- x.moist.LAB@coords[, 1]
# x$A <- x.moist.LAB@coords[, 2]
# x$B <- x.moist.LAB@coords[, 3]
# 
# # get data fraction
# x$data.fraction <- evalMissingData(x, vars=c('L', 'A', 'B'), name = 'name')
# 
# 
# 
# 
# # check: OK
# par(mar=c(0,0,3,0))
# plot(x, color='A', plot.order=order(x$data.fraction))
# 
# # subset profiles with nearly complete data
# x <- subsetProfiles(x, s = "data.fraction > 0.95")
# 
# # distance matrix and divisive clustering
# # rescale D to {0,1}
# d <- profile_compare(x, vars=c('L', 'A', 'B'), max_d=100, k=0, rescale.result=TRUE)
# 
# plotProfileDendrogram(x, d, scaling.factor = 0.005, max.depth=100)
# 



