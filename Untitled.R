library(isoWater)
library(assignR)
library(terra)
library(gstat)

# Read and prep data ----
## Get tap water data for USA
tw = wiDB_data(countries = "US", types = "Tap")
twd = tw$data

## Remove a few proprietary data
max(twd$d18O, na.rm = TRUE)
min(twd$d18O, na.rm = TRUE)
twd$d18O[twd$d18O == 9999] = NA

max(twd$d2H, na.rm = TRUE)
min(twd$d2H, na.rm = TRUE)
twd$d2H[twd$d2H == 9999] = NA

## How many records?
sum(!is.na(twd$d2H))
sum(!is.na(twd$d18O))

## Check w/ plot, looks good
plot(twd$d18O, twd$d2H)

## Remove any sites w/ no data
twd = twd[!(is.na(twd$d2H) & is.na(twd$d18O)),]

## Find sites w/ same coords and average them
tw.dups = duplicated(twd[, 3:4])
twd.keep = twd[!tw.dups,]
for(i in 1:length(twd.keep)){
  kr = twd.keep[i,]
  rr = twd[twd[, 3] == twd.keep[i, 3] & twd[, 4] == twd.keep[i, 4], ]
  twd.keep$d2H[i] = mean(rr$d2H, na.rm = TRUE)
  twd.keep$d18O[i] = mean(rr$d18O, na.rm = TRUE)
}

## Make spatial
tw.sp = vect(twd.keep[,c(1,3,4,15,16)], geom = c("Longitude", "Latitude"), crs = "WGS84")

## Check again
points(tw.sp$d18O, tw.sp$d2H, col = "red")
## Map
plot(states)
points(tw.sp)

# Modeling ----
## Get precipitation grids
iso = getIsoscapes("USPrecipMA")
states.proj = project(states, iso)

##Map of NorthAmerica
worldvect <- world(path=tempdir())
usmap <- subset(worldvect, worldvect$NAME_0 == "United States")
plot(usmap)
usmap1 = crop(namap, c(-125, -24, 66, 50))
usmap1 = project (namap1, "ESRI:102008")
bb=ext(namap1)
xmin(bb)=-5e6
namap1=crop(namap1, bb)
tw.sp =project(tw.sp, namap1)
plot(usmap1)
points(tw.sp)
library(sf)

# Load the world shapefile
worldvect <- st_read(system.file("shape/world", package = "maptools"))

# Subset to get the shapefile for the United States
usmap <- subset(worldvect, ADMIN == "United States")

# Plot the shapefile for the United States
plot(usmap)

# Crop to the contiguous United States
usmap1 <- st_crop(usmap, xmin = -125, ymin = 24, xmax = -66, ymax = 50)

# Project to an appropriate coordinate system (e.g., USA Contiguous Albers Equal Area Conic)
usmap1 <- st_transform(usmap1, crs = st_crs(usmap1)$USA_Contiguous_Albers_Equal_Area_Conic)

# Plot the cropped and projected shapefile
plot(usmap1)


## Project and spatial subset data
tw.sp = project(tw.sp, iso)
tw.sp = tw.sp[states.proj,]
plot(states.proj)
points(tw.sp)

## Extract values to calculate residuals
tw.sp$d2Hpcp = extract(iso$d2h_MA, tw.sp, ID = FALSE)
tw.sp$d18Opcp = extract(iso$d18o_MA, tw.sp, ID = FALSE)

## Look at results
plot(tw.sp$d2Hpcp, tw.sp$d2H)
abline(0, 1)
plot(tw.sp$d18Opcp, tw.sp$d18O)
abline(0, 1)

## Residuals
tw.sp$d2Hres = tw.sp$d2H - tw.sp$d2Hpcp
tw.sp$d18Ores = tw.sp$d18O - tw.sp$d18Opcp

## Variograms
hsub = data.frame(geom(tw.sp[!is.na(tw.sp$d2Hres),])[, 3:4])
hsub = cbind(hsub, values(tw.sp[!is.na(tw.sp$d2Hres),]))
hvar = variogram(d2Hres ~ 1, ~x+y, hsub, cutoff = 2e6,
                 width = 2e6 / 20)
plot(hvar)

osub = data.frame(geom(tw.sp[!is.na(tw.sp$d18Ores),])[, 3:4])
osub = cbind(osub, values(tw.sp[!is.na(tw.sp$d18Ores),]))
ovar = variogram(d18Ores ~ 1, ~x+y, osub, cutoff = 2e6,
                 width = 2e6 / 20)
plot(ovar)

## Starter variogram models
hvgm = vgm(250, "Sph", 1.7e6, 80)
plot(hvar, hvgm)
ovgm = vgm(4, "Sph", 1.7e6, 2)
plot(ovar, ovgm)

## Fit them
hvgm.f = fit.variogram(hvar, hvgm, fit.method = 6)
plot(hvar, hvgm.f)

ovgm.f = fit.variogram(ovar, ovgm, fit.method = 6)
plot(ovar, ovgm.f)

## Prediction grid 25 km x 25 km
pgg = aggregate(iso$d2h_MA, 25)

## Krige
kh = gstat(formula = d2Hres ~ 1, locations = ~x + y, data = hsub, model = hvgm.f)
ph = interpolate(pgg, kh)
plot(ph)

ko = gstat(formula = d18Ores ~ 1, locations = ~x + y, data = osub, model = ovgm.f)
po = interpolate(pgg, ko)
plot(po)

## Downscale
ph = resample(ph, iso)
po = resample(po, iso)

## Mask
ph = mask(ph, iso$d2h_MA)
po = mask(po, iso$d18o_MA)

## Sum isoscape and interpolated residuals
htap = iso$d2h_MA + ph$var1.pred
otap = iso$d18o_MA + po$var1.pred

# Uncertainty ----
## Tap isoscape residuals 
hres.tap = tw.sp$d2H - extract(htap, tw.sp, ID = FALSE)
ores.tap = tw.sp$d18O - extract(otap, tw.sp, ID = FALSE)
hres.tap = hres.tap[!is.na(hres.tap)]
ores.tap = ores.tap[!is.na(ores.tap)]

## Compare to normal...heavy tailed
dev.off()
plot(density(hres.tap))
lines(density(rnorm(1e6, 0, sd(hres.tap))), col = "red")

## Additive variance
htap.sd = sqrt(ph$var1.var + var(hres.tap))
otap.sd = sqrt(po$var1.var + var(ores.tap))

# Combine and write ----
UStap = c(htap, sqrt(ph$var1.var), htap.sd,
          otap, sqrt(po$var1.var), otap.sd)
names(UStap) = c("d2h.m", "d2h.se", "d2h.sd", "d18o.m", "d18o.se", "d18o.sd")

#Write out
writeRaster(UStap, "out/UStap.tif")
