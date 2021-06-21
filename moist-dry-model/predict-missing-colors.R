
##
## 2020-03-30: no longer using this, switching to OLS modeling
##

stop('no longer using RF, see OLS model')

##
## 2015-06-23
## fill-in missing / incorrectly parsed OSD colors using brute-force supervised classification
##

## TODO: save model object for times when an update would suffice and we don't have all of the samples
## TODO: better model accuracy reporting based on CV or side-by-side eval of common colors

# toggle model re-fitting, takes ~ 10 minutes
reFit <- TRUE


library(randomForest)
library(hexbin)
library(viridis)
library(aqp)
library(sharpshootR)

# from OSDs
# keep factors so that models function as expected
d <- read.csv('parsed-data.csv.gz', stringsAsFactors=TRUE)

# initial conditions
summary(d)

## brute force prediction
## takes ~ 10 minutes for model fitting
## consider use of LAB colorspace for models

if(reFit) {
  # moist from dry
  mv.rf <- randomForest(moist_value ~ dry_value + dry_chroma + dry_hue, data=d, na.action=na.omit)
  mc.rf <- randomForest(moist_chroma ~ dry_value + dry_chroma + dry_hue, data=d, na.action=na.omit)
  
  # dry from moist
  dv.rf <- randomForest(dry_value ~ moist_value + moist_chroma + moist_hue, data=d, na.action=na.omit)
  dc.rf <- randomForest(dry_chroma ~ moist_value + moist_chroma + moist_hue, data=d, na.action=na.omit)
  
  # save for next time
  save(mv.rf, mc.rf, dv.rf, dc.rf, file='models/missing-color-RF-models.rda')
} else {
  # load the models from last time
  load('models/missing-color-models.rda')
}


## save a record of model accuracy
sink(file = 'QC/rf-model-accuracy.txt')
# moist colors
print(mv.rf)
print(mc.rf)
# dry colors
print(dv.rf)
print(dc.rf)
sink()

## graphical eval: seems reasonable

## RF model
png(filename = 'figures/mv-model.png', width=800, height=800, res=90)
hexbinplot(predict(mv.rf, newdata=d) ~ d$moist_value, trans = log, inv=exp, colramp=viridis, asp=1, xbins=15, xlab='Observed Moist Value', ylab='Predicted Moist Value')
dev.off()

png(filename = 'figures/dv-model.png', width=800, height=800, res=90)
hexbinplot(predict(dv.rf, newdata=d) ~ d$dry_value, trans = log, inv=exp, colramp=viridis, asp=1, xbins=15, xlab='Observed Dry Value', ylab='Predicted Dry Value')
dev.off()

png(filename = 'figures/mc-model.png', width=800, height=800, res=90)
hexbinplot(predict(mc.rf, newdata=d) ~ d$moist_chroma, trans = log, inv=exp, colramp=viridis, asp=1, xbins=15, xlab='Observed Moist Chroma', ylab='Predicted Moist Chroma')
dev.off()

png(filename = 'figures/dc-model.png', width=800, height=800, res=90)
hexbinplot(predict(dc.rf, newdata=d) ~ d$dry_chroma, trans = log, inv=exp, colramp=viridis, asp=1, xbins=15, xlab='Observed Dry Chroma', ylab='Predicted Dry Chroma')
dev.off()

## make a copy of some of the data,  
x.original <- subset(d, subset = seriesname %in% c('AMADOR', 'DRUMMER', 'CECIL', 'REDDING', 'AVA', 'MIAMI', 'FRISCO'))

# promote to SPC and convert colors
depths(x.original) <- seriesname ~ top + bottom
x.original$dry_soil_color <- munsell2rgb(x.original$dry_hue, x.original$dry_value, x.original$dry_chroma)
x.original$moist_soil_color <- munsell2rgb(x.original$moist_hue, x.original$moist_value, x.original$moist_chroma)

# label
x.original$group <- rep('Original', times=length(x.original))


## fill missing color components via models

# value
d$moist_value[which(is.na(d$moist_value))] <- round(predict(mv.rf, d[which(is.na(d$moist_value)), ]))
d$dry_value[which(is.na(d$dry_value))] <- round(predict(dv.rf, d[which(is.na(d$dry_value)), ]))

# chroma
d$moist_chroma[which(is.na(d$moist_chroma))] <- round(predict(mc.rf, d[which(is.na(d$moist_chroma)), ]))
d$dry_chroma[which(is.na(d$dry_chroma))] <- round(predict(dc.rf, d[which(is.na(d$dry_chroma)), ]))

## convert factors -> character
d$moist_hue <- as.character(d$moist_hue)
d$dry_hue <- as.character(d$dry_hue)
d$name <- as.character(d$name)
d$seriesname <- as.character(d$seriesname)

# copy vs. prediction of hue, use moist / dry hue
d$moist_hue[which(is.na(d$moist_hue))] <- d$dry_hue[which(is.na(d$moist_hue))]
d$dry_hue[which(is.na(d$dry_hue))] <- d$moist_hue[which(is.na(d$dry_hue))]


## filling missing O horizon colors requires fixing 0->O OCR errors
idx <- grep('^0', d$name)
sort(table(d$name[idx]), decreasing = TRUE)

# repalce 0 with O
d$name[idx] <- gsub('0', 'O', d$name[idx])


## O horizon colors: moist and dry colors missing

# find some to eval
x.o <- d[grep('^O', d$name), ]
nrow(x.o)
head(x.o)
sort(table(x.o$name), decreasing = TRUE)

# generalize into a 3 classes + everything else
x.o$genhz <- generalize.hz(x.o$name, new=c('Oi', 'Oe', 'Oa'), pat = c('Oi', 'Oe', 'Oa'))

# convert colors
x.o$dry_color <- munsell2rgb(x.o$dry_hue, x.o$dry_value, x.o$dry_chroma)
x.o$moist_color <- munsell2rgb(x.o$moist_hue, x.o$moist_value, x.o$moist_chroma)

# split and upgrade to SPC
x.o.d <- subset(x.o, subset=! is.na(dry_color) & !is.na(top) & !is.na(bottom))
x.o.m <- subset(x.o, subset=! is.na(moist_color) & !is.na(top) & !is.na(bottom))

depths(x.o.d) <- seriesname ~ top + bottom
depths(x.o.m) <- seriesname ~ top + bottom

# aggregate colors
a.d <- aggregateColor(x.o.d, groups='genhz', col='dry_color', k=10)
a.m <- aggregateColor(x.o.d, groups='genhz', col='moist_color', k=10)

png(file='figures/O-hz-colors-dry.png', width = 900, height=550, res=90)
aggregateColorPlot(a.d, main='Dry Colors')
dev.off()

png(file='figures/O-hz-colors-moist.png', width = 900, height=550, res=90)
aggregateColorPlot(a.m, main='Moist Colors')
dev.off()

knitr::kable(a.d$aggregate.data)
knitr::kable(a.m$aggregate.data)

## find O horizons that are missing colors, and use these ones

# Oi / dry
idx <- which(grepl('Oi', d$name) & is.na(d$dry_hue))
d$dry_hue[idx] <- '10YR'
d$dry_value[idx] <- 4
d$dry_chroma[idx] <- 2

# Oi / moist
idx <- which(grepl('Oi', d$name) & is.na(d$moist_hue))
d$moist_hue[idx] <- '7.5YR'
d$moist_value[idx] <- 2
d$moist_chroma[idx] <- 2

# Oe / dry
idx <- which(grepl('Oe', d$name) & is.na(d$dry_hue))
d$dry_hue[idx] <- '7.5YR'
d$dry_value[idx] <- 4
d$dry_chroma[idx] <- 2

# Oe / moist
idx <- which(grepl('Oe', d$name) & is.na(d$moist_hue))
d$moist_hue[idx] <- '7.5YR'
d$moist_value[idx] <- 2
d$moist_chroma[idx] <- 2

# Oa / dry
idx <- which(grepl('Oa', d$name) & is.na(d$dry_hue))
d$dry_hue[idx] <- '7.5YR'
d$dry_value[idx] <- 4
d$dry_chroma[idx] <- 2

# Oa / moist
idx <- which(grepl('Oa', d$name) & is.na(d$moist_hue))
d$moist_hue[idx] <- '10YR'
d$moist_value[idx] <- 2
d$moist_chroma[idx] <- 1

# everything else, dry
idx <- which(grepl('O', d$name) & is.na(d$dry_hue))
d$dry_hue[idx] <- '7.5YR'
d$dry_value[idx] <- 4
d$dry_chroma[idx] <- 2

# everything else, moist
idx <- which(grepl('O', d$name) & is.na(d$moist_hue))
d$moist_hue[idx] <- '7.5YR'
d$moist_value[idx] <- 2
d$moist_chroma[idx] <- 2


##
## extract same series and compare original vs. filled colors
##
x <- subset(d, subset = seriesname %in% c('AMADOR', 'DRUMMER', 'CECIL', 'REDDING', 'AVA', 'MIAMI', 'FRISCO'))
x$seriesname <- paste0(x$seriesname, '-filled')
depths(x) <- seriesname ~ top + bottom
x$dry_soil_color <- munsell2rgb(x$dry_hue, x$dry_value, x$dry_chroma)
x$moist_soil_color <- munsell2rgb(x$moist_hue, x$moist_value, x$moist_chroma)

# label
x$group <- rep('Filled', times=length(x))

# stack
g <- union(list(x.original, x))

## graphical comparison... still needs some work

png(file='figures/dry-original-vs-filled-example.png', width = 900, height=800, res=90)

par(mar=c(1,1,3,1), mfrow=c(2,1))
groupedProfilePlot(g, groups='group', color='dry_soil_color', id.style='side') ; title('Dry Colors')
groupedProfilePlot(g, groups='group', color='moist_soil_color') ; title('Moist Colors')

dev.off()


png(file='figures/original-dry-vs-moist.png', width = 900, height=800, res=90)

par(mar=c(2,1,3,1), mfrow=c(2,1))
plot(x.original, color='dry_soil_color', max.depth=165)
title('Original Dry Colors')
plot(x.original, color='moist_soil_color', max.depth=165)
title('Original Moist Colors')

dev.off()


png(file='figures/filled-dry-vs-moist.png', width = 900, height=800, res=90)

par(mar=c(2,1,3,1), mfrow=c(2,1))
plot(x, color='dry_soil_color', max.depth=165)
title('Filled Dry Colors')
plot(x, color='moist_soil_color', max.depth=165)
title('Filled Moist Colors')

dev.off()


## TODO: illustrate missing colors / filled colors / predictions

par(mar=c(1,1,1,1))
plot(expand.grid(x=1:36, y=1:2), xlim=c(0.5,36.5), ylim=c(0.5, 5), axes=FALSE, type='n')
points(expand.grid(x=1:36, y=1), pch=22, cex=3, bg=x.original$dry_soil_color)
points(expand.grid(x=1:36, y=2), pch=22, cex=3, bg=x$dry_soil_color)

par(mar=c(1,1,1,1))
plot(expand.grid(x=1:36, y=1:2), xlim=c(0.5,36.5), ylim=c(0.5, 5), axes=FALSE, type='n')
points(expand.grid(x=1:36, y=1), pch=22, cex=3, bg=x.original$moist_soil_color)
points(expand.grid(x=1:36, y=2), pch=22, cex=3, bg=x$moist_soil_color)




## save results
write.csv(d, file=gzfile('parsed-data-est-colors.csv.gz'), row.names=FALSE)




