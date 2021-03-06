library(soilDB)
library(plyr)
library(cluster)
library(sharpshootR)
library(reshape2)
library(maps)
library(colorspace)


## TODO:
## using simplifyColors=TRUE makes the code much simpler

# state soils
ss <- read.csv('state_soils.csv', stringsAsFactors = FALSE)

getPedons <- function(x) {
  suppressMessages(fetchKSSL(series=x, returnMorphologicData = TRUE))
}

# extract a morphologic table from each record set
# NULL data are filtered out
# converted to data.frame
extractMorphTable <- function(x, table='phcolor') {
  m <- lapply(x, function(i) i[['morph']][[table]])
  # index pedons with data
  idx <- which(sapply(m, length) > 0) 
  if(length(idx) < 1)
    return(NULL)
  # filter non-NULL and convert to DF
  m <- ldply(m[idx])
  return(m)
}


# get KSSL + morph for all state soils
x <- lapply(ss$series, getPedons)

# compile color data
phcolor <- extractMorphTable(x, table='phcolor')

# extract pedons from list
pedons <- lapply(x, function(i) i$SPC)

# combine into a single SPC object
pedons <- aqp::union(pedons)

## TODO: there may be some issues here with labsampnum
# simplify combined color data
colors <- simplifyColorData(phcolor, id.var = 'labsampnum', wt='colorpct')

# join combined color data to SPC
h <- horizons(pedons)
h <- join(h, colors, by='labsampnum', type='left', match='first')
horizons(pedons) <- h

# normalize taxonname / series
pedons$taxonname <- toupper(pedons$taxonname)
table(pedons$taxonname)

# re-name for join
pedons$series <- pedons$taxonname
pedons$taxonname <- NULL

# rename state derived from GIS intersection
pedons$gis_state <- pedons$state
pedons$state <- NULL

# join state from original LUT back in
ss$series <- toupper(ss$series)
site(pedons) <- ss


previewColors(pedons$moist_soil_color, method='MDS')



# check: looks good
png(file='state-soil-kssl-data-eval-colors.png', width = 800, height=400, type = 'cairo', antialias = 'subpixel', res = 70)

par(mar=c(0.5,0.5,1,0.5))
groupedProfilePlot(pedons[1:40, ], groups = 'state', color = "moist_soil_color", print.id=FALSE, max.depth=150)

dev.off()

## aggregate moist color by state
# convert soil colors to LAB colorspace
# moist colors
cols.lab <- convertColor(cbind(pedons$m_r, pedons$m_g, pedons$m_b), from='sRGB', to = 'Lab', from.ref.white = 'D65', to.ref.white = 'D65', clip = FALSE)
pedons$m_L <- cols.lab[, 1]
pedons$m_A <- cols.lab[, 2]
pedons$m_B <- cols.lab[, 3]


# aggregate data by normalized taxonname, via slice-wise median
a.colors <- slab(pedons, state ~ m_r + m_g + m_b + m_L + m_A + m_B, slab.fun = median, na.rm=TRUE)

# throw out aggregate data that are deeper than 150cm
a.colors <- subset(a.colors, subset=bottom < 150)

# convert long -> wide format
x.colors <- dcast(a.colors, state + top + bottom ~ variable, value.var = 'value')


# composite sRGB triplets into an R-compatible color
# note that missing colors must be padded with NA
x.colors$soil_color <- NA
not.na <- which(complete.cases(x.colors[, c('m_L', 'm_A', 'm_B')]))
cols.srgb <- data.frame(convertColor(cbind(x.colors$m_L, x.colors$m_A, x.colors$m_B), from='Lab', to = 'sRGB', from.ref.white = 'D65', to.ref.white = 'D65', clip = FALSE))
names(cols.srgb) <- c('R', 'G', 'B')


x.colors$soil_color[not.na] <- with(cols.srgb[not.na, ], rgb(R, G, B, maxColorValue = 1))

# init a new SoilProfileCollection from aggregate data
depths(x.colors) <- state ~ top + bottom

# not bad
par(mar=c(1,0,3,4))
plot(x.colors, divide.hz=FALSE, name='', col.label='Soil Color', lwd=1.25, axis.line.offset=0, cex.depth.axis=1, cex.id=1)




## simple color signature
pig <- soilColorSignature(x.colors, r = 'm_r', g = 'm_g', b='m_b')


# account for missing data
idx <- which(complete.cases(pig[, -1]))
pig <- pig[idx, ]

# the first column is the ID
row.names(pig) <- pig[, 1]
d <- daisy(pig[, 2:6])
dd <- diana(d)

# index to those profiles present in `d`
idx <- which(profile_id(x.colors) %in% pig$state)
sp <- x.colors[idx, ]

plotProfileDendrogram(sp, dd, dend.y.scale = 0.5, divide.hz=FALSE, scaling.factor = 0.001, y.offset = 0.01, width=0.15)

dd$order.lab
dd$order

# 
# ## all colors from all pedons
# a <- aggregateColor(pedons, groups = 'state', col = 'moist_soil_color')
# 
# 
# 
# par(mar=c(2, 6, 1, 1))
# aggregateColorPlot(a, print.label = FALSE)
# 



# # subset / slice for aggregation
# pedons$genhz <- generalize.hz()
# 
# 
# 
# 
# 
# 
# 
# s <- slice(pedons, c(5, 10, 15, 25, 50, 75) ~ moist_soil_color, strict = FALSE)
# # make horizon labels based on slice depth
# s$slice <- paste0(s$hzn_top, " cm")
# s$slice <- factor(s$slice, levels = guessGenHzLevels(s, "slice")$levels)
# 
# ## Error: SPC object corruption. This shouldn't happen and will be fixed in aqp 2.0
# sp <- subsetProfiles(s, h="slice == '10 cm'")
# 
# a <- aggregateColor(sp, groups = 'state', col = 'moist_soil_color')
# par(mar=c(2, 6, 1, 1))
# aggregateColorPlot(a, print.label = FALSE)
# 





### work with OSDs

# get these soil series
s <- fetchOSD(unique(ss$series))

# join in state
s$series <- profile_id(s)
site(s) <- ss

# manually convert Munsell -> sRGB
rgb.data <- munsell2rgb(s$hue, s$value, s$chroma, return_triplets = TRUE)
s$r <- rgb.data$r
s$g <- rgb.data$g
s$b <- rgb.data$b

# check
par(mar=c(1,1,1,1))
plot(s)



rgb.colors <- munsell2rgb(s$hue, s$value, s$chroma, return_triplets = TRUE)
lab.colors <- as(RGB(rgb.colors[['r']], rgb.colors[['g']], rgb.colors[['b']]), 'LAB')@coords
cols <- cbind(rgb.colors, lab.colors)
cols <- na.omit(cols)
cols <- as.data.frame(cols)

png(file='state-soils-osd-color-LAB-palette.png', width=800, height=800, res=90, type='cairo', antialias = 'subpixel')

pairs(~ L + A + B, data=cols, pch=16, cex=3, col=rgb(cols$r, cols$g, cols$b))

dev.off()

# generate color signatures using OSDs and PAM method
pig <- soilColorSignature(s, RescaleLightnessBy = 5, method='pam')

# move row names over for distance matrix
row.names(pig) <- pig[, 1]
d <- daisy(pig[, -1])
dd <- diana(d)


par(mar=c(1,1,1,1))

plotProfileDendrogram(s, dd, dend.y.scale = max(d) * 2, scaling.factor = 0.25, y.offset = 6, width=0.15, cex.names=0.45, label='state', name='')


par(mar=c(1,1,1,1))
plot(s, plot.order=dd$order, label='state', name='')


# set order of states based on clustering order
ll <- ss$state[match(dd$order.lab, ss$series)]
s$state <- factor(s$state, levels=ll)


# aggregate soil color based on sate
a.osd <- aggregateColor(s, groups='state', col='soil_color')

# plot, vertical axis will be in order of dendrogram leaves
png(file='state-soils-osd-signatures.png', width = 1000, height=800, type = 'cairo', antialias = 'subpixel', res = 90)

par(mar=c(0.5, 6, 1, 0.5))
aggregateColorPlot(a.osd, print.label = FALSE, x.axis = FALSE, rect.border = NA, horizontal.borders = TRUE, horizontal.border.lwd = 1)
title(main='State Soil Color Signatures: OSD', line=-1, cex.main=2)

dev.off()


png(file='state-soils-osd-signatures-inverse.png', width = 1000, height=900, type = 'quartz', antialias = 'subpixel', res = 90)

par(mar=c(0.5, 6, 1, 0.5), bg='black', fg='white')
aggregateColorPlot(a.osd, print.label = FALSE, x.axis = FALSE, rect.border = NA, horizontal.borders = TRUE, horizontal.border.lwd = 1)

dev.off()





## all colors from KSSL pedons
pedons$state <- factor(pedons$state, levels=ll)
a.kssl <- aggregateColor(pedons, groups = 'state', col = 'moist_soil_color')
a.kssl.15 <- aggregateColor(pedons, groups = 'state', col = 'moist_soil_color', k = 15)

png(file='state-soils-kssl-signatures.png', width = 1000, height=800, type = 'cairo', antialias = 'subpixel', res = 90)

par(mar=c(0.5, 6, 1, 0.5))
aggregateColorPlot(a.kssl, print.label = FALSE, x.axis = FALSE, rect.border = NA, horizontal.borders = TRUE, horizontal.border.lwd = 1)
title(main='State Soil Color Signatures: KSSL', line=-1, cex.main=2)

dev.off()


png(file='state-soils-kssl-signatures-15.png', width = 1000, height=800, type = 'cairo', antialias = 'subpixel', res = 90)

par(mar=c(0.5, 6, 1, 0.5), bg='black', fg='white')
aggregateColorPlot(a.kssl.15, print.label = FALSE, x.axis = FALSE, rect.border = NA, horizontal.borders = TRUE, horizontal.border.lwd = 1)

dev.off()







## single color / state based on KSSL morph

a.aggregate <- a.kssl$aggregate.data

a.aggregate$munsell <- paste0(a.aggregate$munsell.hue, ' ', a.aggregate$munsell.value, '/', a.aggregate$munsell.chroma)

# make a grid for plotting
n <- ceiling(sqrt(nrow(a.aggregate)))
# read from top-left to bottom-right
g <- expand.grid(x=1:n, y=n:1)[1:nrow(a.aggregate),]




png(file='state-soils-single-color.png', width = 1000, height=900, type = 'cairo', antialias = 'subpixel', res = 90)

par(mar=c(1,0,1,1))
plot(g$x, g$y, pch=15, cex=12, axes=FALSE, xlab='', ylab='', col=a.aggregate$col, xlim=c(0.5,8.5), ylim=c(1.5,8.5))
text(g$x, g$y, a.aggregate$state, adj=c(0.45,5), cex=1, font=2)
text(g$x, g$y, a.aggregate$munsell, col='white', pos=1, cex=0.85, font=2)
title(main='State Soil Colors', line=-1, cex.main=2)

dev.off()



## simple map with maps package

# get state names in order of plotting
state.map <- map('state', plot=FALSE)$names
# clean cruft from names
state.map <- strsplit(state.map, ':')
state.map <- sapply(state.map, function(i) i[[1]])

# index mapping states to colors
col.idx <- match(state.map, tolower(a.aggregate$state))


png(file='state-soils-single-color-map.png', width = 1000, height=900, type = 'cairo', antialias = 'subpixel', res = 90)

par(mar=c(1,0,1,1))
map('state', col=a.aggregate$col[col.idx], fill=TRUE, mar=c(1,1,2,1))
title(main='State Soil Colors', line=1, cex.main=2)

dev.off()


##
## geofacet with state soil colours + Munsell
##
library(ggplot2)
library(geofacet)

state_cols <- a.aggregate$col
names(state_cols) <- a.aggregate$state

geoplot <- ggplot(a.aggregate) + 
  geom_rect(aes(xmin = 0, xmax = 1, ymin = 0, ymax = 1, fill = state)) + 
  geom_text(aes(x = 0.5, y = 0.5, label = munsell), colour = "#ffffff") +
  facet_geo(~state, strip.position = "bottom") + 
  scale_fill_manual(
    guide = FALSE, 
    values = state_cols
  ) + 
  coord_equal() +
  theme_bw() + 
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    panel.border = element_rect(fill = NA, colour = NA),
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank()
  )

png("./geofacet-soils.png", res = 200, width = 2500, height = 2000, type='cairo', antialias = 'subpixel')
print(geoplot)
dev.off()


## geofacet with treemaps of main soil colours per state

library(treemapify)

soil_data_states <- do.call(rbind, a.kssl$scaled.data)
soil_data_states$state <- stringr::str_remove(stringr::str_extract(row.names(soil_data_states), '.*[.]'), "\\.")

treemap_cols <- soil_data_states$moist_soil_color
names(treemap_cols) <- munsell

geoplot_treemap <- ggplot(data = soil_data_states) + 
  geom_treemap(aes(area = weight, fill = munsell)) +
  facet_geo(~state, strip.position = "bottom") + 
  scale_fill_manual(
    guide = FALSE, 
    values = treemap_cols
  ) + 
  coord_equal() +
  theme_bw() + 
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    panel.border = element_rect(fill = NA, colour = NA),
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank()
  )

png("./geofacet-treemap-soils.png", res = 300, width = 5000, height = 4000)
print(geoplot_treemap)
dev.off()
