library(aqp)
library(farver)
library(Hmisc)
library(rms)

# load model objects
load('models/missing-color-models.rda')

# from OSDs
d <- read.csv('parsed-data.csv.gz', stringsAsFactors=FALSE)

# full data table
x <- na.omit(d[, c('dry_hue', 'dry_value', 'dry_chroma', 'moist_hue', 'moist_value', 'moist_chroma')])


## predictions from full set

# dry
pp.value.dry <- predict(m.value.dry)
pp.chroma.dry <- predict(m.chroma.dry)

# moist
pp.value.moist <- predict(m.value.moist)
pp.chroma.moist <- predict(m.chroma.moist)

# combine
pp.dry <- data.frame(pp.value.dry, pp.chroma.dry)
pp.moist <- data.frame(pp.value.moist, pp.chroma.moist)
z <- cbind(x, pp.dry, pp.moist)
head(z)



## dE00 comparisons

# put all colors in Munsell notation for dE00 metrics
dry.cols <- sprintf("%s %s/%s", z$dry_hue, z$dry_value, z$dry_chroma)
moist.cols <- sprintf("%s %s/%s", z$moist_hue, z$moist_value, z$moist_chroma)

pred.dry.cols <- sprintf("%s %s/%s", z$dry_hue, round(z$pp.value.dry), round(z$pp.chroma.dry))
pred.moist.cols <- sprintf("%s %s/%s", z$moist_hue, round(z$pp.value.moist), round(z$pp.chroma.moist))

## sample for viz
set.seed(35)
s <- sample(1:nrow(z), size = 8)

png(filename = 'figures/dE00-eval-example-001.png', width=1000, height=800, res=90)
par(mfrow=c(3,1))
colorContrastPlot(dry.cols[s], moist.cols[s], labels = c('dry ', 'moist'))

colorContrastPlot(dry.cols[s], pred.dry.cols[s], labels = c('dry ', 'pred '))

colorContrastPlot(moist.cols[s], pred.moist.cols[s], labels = c('moist', 'pred '))
dev.off()


## just dry -> moist contrast, all colors
cc.dry_moist <- colorContrast(dry.cols, moist.cols)

hist(cc.dry_moist$dE00, breaks = 50, las=1, xlab='dE00')
table(cc.dry_moist$cc)
quantile(cc.dry_moist$dE00, na.rm = TRUE, probs = c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1))


## dry vs. predicted dry
cc.dry <- colorContrast(dry.cols, pred.dry.cols)

hist(cc.dry$dE00, breaks = 50, las=1, xlab='dE00')
table(cc.dry$cc)
quantile(cc.dry$dE00, na.rm = TRUE, probs = c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1))

# how does it compare with latent variance?
var(cc.dry$dE00, na.rm = TRUE) / var(cc.dry_moist$dE00, na.rm = TRUE)


## moist vs. predicted moist
