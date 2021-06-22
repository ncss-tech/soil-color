library(aqp)
library(soilDB)
library(sharpshootR)
library(cluster)
library(latticeExtra)
library(reshape2)
library(plyr)
library(grDevices)
library(cluster)
library(MASS)

# all MLRA
x <- read.csv('mlra_drainage_colors.csv.gz', stringsAsFactors = FALSE)

# check data
xtabs(~ mlra + drainagecl, data = x)

# remove records missing drainage class
x <- subset(x, subset=drainagecl != '')

# pick an MLRA
# mlra.set <- c('108A', '110', '108B', '111D') # nice example, based on drummer
mlra.set <- c('136', '133A', '137', '130B', '148') # great example based on cecil
x <- subset(x, mlra %in% mlra.set)

# set drainage class levels
x$drainagecl <- factor(x$drainagecl, levels = c('very poorly', 'poorly', 'somewhat poorly', 'moderately well', 'well', 'somewhat excessively', 'excessively'))

table(x$drainagecl, useNA = 'always')

# convert moist colors to HEX notation
x$soil_color <- munsell2rgb(x$matrix_wet_color_hue, x$matrix_wet_color_value, x$matrix_wet_color_chroma)

# convert moist colors to sRGB
x.rgb <- munsell2rgb(x$matrix_wet_color_hue, x$matrix_wet_color_value, x$matrix_wet_color_chroma, return_triplets = TRUE)

# convert sRGB -> LAB for mixing
x.lab <- convertColor(x.rgb, from='sRGB', to='Lab', from.ref.white = 'D65', to.ref.white = 'D65', clip = FALSE)

# save back to original DF
x$m_L <- x.lab[, 1]
x$m_A <- x.lab[, 2]
x$m_B <- x.lab[, 3]

# init SPC
depths(x) <- series ~ top + bottom
site(x) <- ~ drainagecl

# # compute depth class and join to site data
# dc <- getSoilDepthClass(x, name='hzname', top='top', bottom='bottom')
# site(x) <- dc

# aggregate data by drainage class, via slice-wise median
a.colors <- slab(x, drainagecl ~ m_L + m_A + m_B)

# throw out aggregate data that are deeper than 150cm
a.colors <- subset(a.colors, subset = bottom < 150)

# convert long -> wide format
x.colors <- dcast(a.colors, drainagecl + top + bottom ~ variable, value.var = 'p.q50')

# # check: OK
# head(a.colors)
# head(x.colors)

# composite sRGB triplets into an R-compatible color
# note that missing colors must be padded with NA
x.colors$soil_color <- NA
not.na <- which(complete.cases(x.colors[, c('m_L', 'm_A', 'm_B')]))
cols.srgb <- data.frame(convertColor(cbind(x.colors$m_L, x.colors$m_A, x.colors$m_B), from='Lab', to='sRGB', from.ref.white = 'D65', to.ref.white = 'D65', clip = FALSE))
names(cols.srgb) <- c('R', 'G', 'B')

x.colors$soil_color[not.na] <- with(cols.srgb[not.na, ], rgb(R, G, B, maxColorValue = 1))

# init a new SoilProfileCollection from aggregate data
depths(x.colors) <- drainagecl ~ top + bottom

# generate index for new ordering
new.order <- match(c('very poorly', 'poorly', 'somewhat poorly', 'moderately well', 'well', 'somewhat excessively', 'excessively'), profile_id(x.colors))


mlra.set.text <- paste(mlra.set, collapse = ', ')
mlra.set.text.2 <- paste(mlra.set, collapse = '_')
fname <- paste0('example-', mlra.set.text.2, '.png')

png(file=fname, res=100, width=800, height=400, type='windows', antialias = 'cleartype')

par(mar=c(0,0,1,2))

plotSPC(x.colors, divide.hz=FALSE, name = NA, plot.order=new.order, col.label='Soil Color', lwd=1.25, axis.line.offset=-2, cex.depth.axis=1, cex.id=0.85, id.style='side')
text(1:length(new.order), 0, table(x$drainagecl)[new.order], pos=3, cex=0.85)
title(paste0("Soil Series from MLRAs: ", mlra.set.text), line=0)
title(sub='median moist color, CIELAB color space', line=0)

dev.off()


## aggregate representation of colors by drainage class, at a couple of depth slices
top.seq <- c(5, 10, 25, 50, 75, 100)
s <- slice(x, c(5, 10, 25, 50, 75, 100) ~ ., strict=FALSE)


fname <- paste0('slices-', mlra.set.text.2, '.png')
png(file=fname, res=100, width=800, height=900, type='windows', antialias = 'cleartype')

par(mfrow=c(6,1), mar=c(0.1,10,0.1,0.5))
for(i in seq_along(top.seq)) {
  depth.i <- top.seq[i]
  a <- aggregateColor(s[, i], groups='drainagecl', k=8)
  aggregateColorPlot(a, print.label = FALSE, rect.border = NA, horizontal.borders = TRUE, horizontal.border.lwd=1, label.cex=0.5, x.axis = FALSE)
  mtext(paste0(depth.i, 'cm'), side = 4, line=-1)
}

dev.off()


## same thing, at genhz
## TODO: aggregateColorPlot() will fail when there are too many NA -- e.g. O horizons
## ... skipping those for now

gehz.rules <- c('O', 'A', 'E', 'B', 'C')
x$genhz <- generalize.hz(x$hzname, new=gehz.rules, pat=gehz.rules)

fname <- paste0('genhz-', mlra.set.text.2, '.png')
png(file=fname, res=100, width=800, height=800, type='windows', antialias = 'cleartype')

par(mfrow=c(4,1), mar=c(0.1,10,0.1,0.5))

## note: skipping O horizons
for(i in seq_along(gehz.rules[-1])) {
  hz.i <- gehz.rules[-1][i]
  # make a temporary working copy of full SPC
  xx <- x
  # extract / subset horizons from original
  h <- horizons(x)
  h <- h[which(h$genhz == hz.i), ]
  # replace in temp copy
  horizons(xx) <- h
  # aggregate and plot
  a <- aggregateColor(xx, groups='drainagecl', k=8)
  aggregateColorPlot(a, print.label = FALSE, rect.border = NA, horizontal.borders = TRUE, horizontal.border.lwd=1, label.cex=0.5, x.axis = FALSE)
  mtext(paste0(hz.i), side = 4, line=-1, font=2)
}

dev.off()



##
## adapted from "2016-NSSC-DB-meeting" presentation / OSD reboot ideas
##
h <- as(s, 'data.frame')

# don't need these 
h$m_munsell <- paste0(h$matrix_wet_color_hue, ' ', h$matrix_wet_color_value, '/', h$matrix_wet_color_chroma)
h$matrix_wet_color_hue <- factor(h$matrix_wet_color_hue, levels = c('10Y', '2.5Y', '10YR', '7.5YR', '5YR', '2.5YR'))


## using LAB coordinates

pp <- xyplot(m_L ~ m_B | drainagecl + factor(top), as.table=TRUE, data=h, subscripts = TRUE,  scales=list(alternating=3, draw=FALSE), xlab='CIE A-coordinate', ylab='CIE L-coordinate', panel=function(x, y, subscripts=subscripts, ...) {
  
  p.data <- data.frame(x=x, y=y, col=h$soil_color[subscripts], m=h$m_munsell[subscripts], stringsAsFactors = FALSE)
  tab <- prop.table(table(p.data$m, useNA = 'always'))
  tab <- as.data.frame(tab)
  names(tab) <- c('m', 'freq')
  p.data <- join(p.data, tab, by='m', type='left')
  p.data <- na.omit(p.data)
  p.data <- subset(p.data, subset=freq > 0.05)
  panel.grid(-1, -1)
  panel.xyplot(p.data$x, p.data$y, pch=15, col=scales::alpha(p.data$col, 0.25), cex=2)
})

fname <- paste0('panels-LAB-', mlra.set.text.2, '.png')
png(file=fname, res=100, width=900, height=700, type='windows', antialias = 'cleartype')

useOuterStrips(pp, strip=strip.custom(bg=grey(0.85), par.strip.text=list(cex=0.65)), strip.left = strip.custom(bg=grey(0.85), par.strip.text=list(rot=0)))

dev.off()


##
## even better, use ordination of LAB
##

idx <- which(complete.cases(h[, c('m_L', 'm_A', 'm_B')]))
h <- h[idx, ]
d <- daisy(h[, c('m_L', 'm_A', 'm_B')], stand=TRUE)

## slow!
mds <- cmdscale(d)
h$mds1 <- mds[, 1]
h$mds2 <- mds[, 2]

pp <- xyplot(mds1 ~ mds2 | drainagecl + factor(top), as.table=TRUE, data=h, subscripts = TRUE,  scales=list(alternating=3, draw=FALSE), xlab='', ylab='', panel=function(x, y, subscripts=subscripts, ...) {
  
  p.data <- data.frame(x=x, y=y, col=h$soil_color[subscripts], m=h$m_munsell[subscripts], stringsAsFactors = FALSE)
  tab <- prop.table(table(p.data$m, useNA = 'always'))
  tab <- as.data.frame(tab)
  names(tab) <- c('m', 'freq')
  p.data <- join(p.data, tab, by='m', type='left')
  p.data <- na.omit(p.data)
  p.data <- subset(p.data, subset=freq > 0.05)
  panel.grid(-1, -1)
  panel.xyplot(p.data$x, p.data$y, pch=15, col=scales::alpha(p.data$col, 0.25), cex=2)
})


fname <- paste0('panels-MDS-', mlra.set.text.2, '.png')
png(file=fname, res=100, width=900, height=700, type='windows', antialias = 'cleartype')

useOuterStrips(pp, strip=strip.custom(bg=grey(0.85), par.strip.text=list(cex=0.65)), strip.left = strip.custom(bg=grey(0.85), par.strip.text=list(rot=0)))

dev.off()



##
## even better, very simple hz generalization
##

pp <- xyplot(mds1 ~ mds2 | drainagecl + factor(genhz, levels=gehz.rules), as.table=TRUE, data=h, subscripts = TRUE,  scales=list(alternating=3, draw=FALSE), xlab='', ylab='', panel=function(x, y, subscripts=subscripts, ...) {
  
  p.data <- data.frame(x=x, y=y, col=h$soil_color[subscripts], m=h$m_munsell[subscripts], stringsAsFactors = FALSE)
  tab <- prop.table(table(p.data$m, useNA = 'always'))
  tab <- as.data.frame(tab)
  names(tab) <- c('m', 'freq')
  p.data <- join(p.data, tab, by='m', type='left')
  p.data <- na.omit(p.data)
  p.data <- subset(p.data, subset=freq > 0.05)
  panel.grid(-1, -1)
  panel.xyplot(p.data$x, p.data$y, pch=15, col=scales::alpha(p.data$col, 0.25), cex=2)
})

fname <- paste0('panels-MDS-genhz-', mlra.set.text.2, '.png')
png(file=fname, res=100, width=900, height=700, type='windows', antialias = 'cleartype')

useOuterStrips(pp, strip=strip.custom(bg=grey(0.85), par.strip.text=list(cex=0.65)), strip.left = strip.custom(bg=grey(0.85), par.strip.text=list(rot=0)))

dev.off()



