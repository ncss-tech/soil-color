library(hexbin)
library(viridis)
library(aqp)
library(sharpshootR)
library(rms)
library(farver)

# dry colors
dry.lab <- with(d, munsell2rgb(the_hue = dry_hue, the_value = dry_value, the_chroma = dry_chroma, returnLAB = TRUE))
names(dry.lab) <- c('dry_L', 'dry_A', 'dry_B')

# moist colors
moist.lab <- with(d, munsell2rgb(the_hue = moist_hue, the_value = moist_value, the_chroma = moist_chroma, returnLAB = TRUE))
names(moist.lab) <- c('moist_L', 'moist_A', 'moist_B')

# check: OK
summary(dry.lab)
summary(moist.lab)

# combine: ok
lab <- cbind(dry.lab, moist.lab)
nrow(lab) == nrow(d)

# remove NA: ok
lab <- na.omit(lab)
nrow(lab)

# model
dd <- datadist(lab)
options(datadist="dd")

# dry from moist
(m.L.dry <- ols(dry_L ~ rcs(moist_L) + moist_A + moist_B, data=lab))
(m.A.dry <- ols(dry_A ~ moist_L + rcs(moist_A) + moist_B, data=lab))
(m.B.dry <- ols(dry_B ~ moist_L + moist_A + rcs(moist_B), data=lab))

plot(Predict(m.L.dry))
plot(Predict(m.A.dry))
plot(Predict(m.B.dry))

anova(m.L.dry)
anova(m.A.dry)
anova(m.B.dry)

# moist from dry
(m.L.moist <- ols(moist_L ~ rcs(dry_L) + dry_A + dry_B, data=lab))
(m.A.moist <- ols(moist_A ~ dry_L + rcs(dry_A) + dry_B, data=lab))
(m.B.moist <- ols(moist_B ~ dry_L + dry_A + rcs(dry_B), data=lab))

plot(Predict(m.L.moist))
plot(Predict(m.A.moist))
plot(Predict(m.B.moist))

# errors, why?
# anova(m.L.moist)
anova(m.A.moist)
anova(m.B.moist)

## predictions from full set
pp.L.dry <- predict(m.L.dry)
pp.A.dry <- predict(m.A.dry)
pp.B.dry <- predict(m.B.dry)

pp.L.moist <- predict(m.L.moist)
pp.A.moist <- predict(m.A.moist)
pp.B.moist <- predict(m.B.moist)

pp.dry <- data.frame(pp.L.dry, pp.A.dry, pp.B.dry)
pp.moist <- data.frame(pp.L.moist, pp.A.moist, pp.B.moist)


# combine
z <- cbind(lab, pp.dry, pp.moist)

# operate by row, no need for pair-wise comparisons
dE00.dry <- vector(mode = 'numeric', length = nrow(z))
dE00.moist <- vector(mode = 'numeric', length = nrow(z))

for(i in 1:nrow(z)) {
  dE00.dry[i] <- compare_colour(z[i, c('dry_L', 'dry_A', 'dry_B')], to = z[i, c('pp.L.dry', 'pp.A.dry', 'pp.B.dry')], from_space = 'lab', method = 'CIE2000')
  dE00.moist[i] <- compare_colour(z[i, c('moist_L', 'moist_A', 'moist_B')], to = z[i, c('pp.L.moist', 'pp.A.moist', 'pp.B.moist')], from_space = 'lab', method = 'CIE2000')
}


hist(dE00.dry, breaks = 100, las=1)
hist(dE00.moist, breaks = 100, las=1)

quantile(dE00.dry, na.rm = TRUE, probs = c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1))

# eval back-transformed colors
dry.cols <- rgb2munsell(convert_colour(z[, c('dry_L', 'dry_A', 'dry_B')], from = 'lab', to = 'rgb') / 255)
pred.dry.cols <- rgb2munsell(convert_colour(z[, c('pp.L.dry', 'pp.A.dry', 'pp.B.dry')], from = 'lab', to = 'rgb') / 255)

dry.cols <- sprintf("%s %s/%s", dry.cols$hue, dry.cols$value, dry.cols$chroma)
pred.dry.cols <- sprintf("%s %s/%s", pred.dry.cols$hue, pred.dry.cols$value, pred.dry.cols$chroma)

cc <- colorContrast(dry.cols, pred.dry.cols)

hist(cc$dE00, breaks = 100)
table(cc$cc)
length(which((cc$m1 == cc$m2))) / nrow(cc)
quantile(cc$dE00, na.rm = TRUE, probs = c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1))



## sample for viz
s <- z[sample(1:nrow(z), size = 8), ]

## TODO: back-transformation is still limited to issues with rgb2munsell

dry.cols <- rgb2munsell(convert_colour(s[, c('dry_L', 'dry_A', 'dry_B')], from = 'lab', to = 'rgb') / 255)
pred.dry.cols <- rgb2munsell(convert_colour(s[, c('pp.L.dry', 'pp.A.dry', 'pp.B.dry')], from = 'lab', to = 'rgb') / 255)

dry.cols <- sprintf("%s %s/%s", dry.cols$hue, dry.cols$value, dry.cols$chroma)
pred.dry.cols <- sprintf("%s %s/%s", pred.dry.cols$hue, pred.dry.cols$value, pred.dry.cols$chroma)

colorContrastPlot(dry.cols, pred.dry.cols, labels = c('dry colors', 'predicted dry colors'))


moist.cols <- rgb2munsell(convert_colour(s[, c('moist_L', 'moist_A', 'moist_B')], from = 'lab', to = 'rgb') / 255)
pred.moist.cols <- rgb2munsell(convert_colour(s[, c('pp.L.moist', 'pp.A.moist', 'pp.B.moist')], from = 'lab', to = 'rgb') / 255)

moist.cols <- sprintf("%s %s/%s", moist.cols$hue, moist.cols$value, moist.cols$chroma)
pred.moist.cols <- sprintf("%s %s/%s", pred.moist.cols$hue, pred.moist.cols$value, pred.moist.cols$chroma)

colorContrastPlot(moist.cols, pred.moist.cols, labels = c('moist\ncolors', 'predicted\nmoist colors'))





hexbinplot(pp.L.dry ~ lab$dry_L, trans = log, inv=exp, colramp=viridis, asp=1, xbins=10, xlab='Observed Dry L', ylab='Predicted Dry L')
hexbinplot(pp.A.dry ~ lab$dry_A, trans = log, inv=exp, colramp=viridis, asp=1, xbins=10, xlab='Observed Dry A', ylab='Predicted Dry A')
hexbinplot(pp.B.dry ~ lab$dry_B, trans = log, inv=exp, colramp=viridis, asp=1, xbins=10, xlab='Observed Dry B', ylab='Predicted Dry B')

