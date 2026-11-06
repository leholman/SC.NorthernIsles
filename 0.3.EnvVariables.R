###############################################################################
## ENV data output
###############################################################################

## Libraries 

library(ncdf4)
library(ggOceanMaps)
library(sf)
library(ggspatial)


sites <- read.csv("metadata_site.csv")


###############################################################################
## LGMR palaeoclimate (Osman et al. 2021) - one CSV per site, per field
## Nearest grid cell WITH DATA, so coastal sites never come back NaN.
###############################################################################

dir.create("cleandata/Climate", recursive = TRUE, showWarnings = FALSE)

## great-circle distance (km) from one point to every cell of a grid
gc_dist <- function(lon0, lat0, lon, lat) {
  p <- pi / 180
  a <- sin((lat - lat0) * p / 2)^2 +
    cos(lat0 * p) * cos(lat * p) * sin((lon - lon0) * p / 2)^2
  2 * 6371 * asin(pmin(sqrt(a), 1))
}

## force a coordinate variable into a matrix matching x[,,1].
## covers 2D (either orientation) and 1D (vector or 1D array, either axis).
grid_coord <- function(v, x, what) {
  n1 <- dim(x)[1]
  n2 <- dim(x)[2]
  d  <- as.integer(dim(v))
  if (identical(d, c(n1, n2))) return(v)
  if (identical(d, c(n2, n1))) return(t(v))
  v <- as.vector(v)                                  # drop any dim attribute
  if (length(v) == n1) return(matrix(v, n1, n2))              # varies by row
  if (length(v) == n2) return(matrix(v, n1, n2, byrow = TRUE)) # varies by col
  stop(what, ": coords are ", paste(dim(v), collapse = "x"), " / length ",
       length(v), ", data grid is ", n1, "x", n2)
}

## load one LGMR field: values, SD, ages, coordinate matrices, land mask
read_lgmr <- function(file, var) {
  nc  <- nc_open(file)
  x   <- ncvar_get(nc, var)
  sd  <- ncvar_get(nc, paste0(var, "_std"))
  age <- ncvar_get(nc, "age")
  lon <- grid_coord(ncvar_get(nc, "lon"), x, paste(var, "lon"))
  lat <- grid_coord(ncvar_get(nc, "lat"), x, paste(var, "lat"))
  nc_close(nc)
  
  ## cells with a complete time series; everything else is land or masked
  wet <- as.vector(rowSums(is.finite(x) & x < 1e30, dims = 2) == dim(x)[3])
  
  cat(sprintf("%s: %d x %d grid, %d cells with data\n",
              var, dim(x)[1], dim(x)[2], sum(wet)))
  list(x = x, sd = sd, age = age, lon = lon, lat = lat, wet = wet)
}

## time series at the nearest cell that has data
extract <- function(f, lon0, lat0, nm) {
  d <- as.vector(gc_dist(lon0, lat0, f$lon, f$lat))
  d[!f$wet] <- NA
  i  <- which.min(d)                              # which.min skips NA
  ix <- ((i - 1L) %%  nrow(f$lon)) + 1L
  iy <- ((i - 1L) %/% nrow(f$lon)) + 1L
  
  m <- f$x[ix, iy, ]
  s <- f$sd[ix, iy, ]
  stopifnot(all(is.finite(m)))                    # never write a column of NAs
  out <- data.frame(age = f$age, m = m, sd = s, upr = m + 2 * s, lwr = m - 2 * s,
                    cell_lon = f$lon[i], cell_lat = f$lat[i], dist_km = d[i])
  names(out)[2:5] <- paste0(nm, c("", "_sd", "_upr", "_lwr"))
  out
}

SST <- read_lgmr("rawdata/Climate/LGMR_SST_climo.nc", "sst")
SAT <- read_lgmr("rawdata/Climate/LGMR_SAT_climo.nc", "sat")

for (i in seq_len(nrow(sites))) {
  
  id   <- sites$CoreID[i]
  lon0 <- sites$lon[i]
  lat0 <- sites$lat[i]
  
  sst <- extract(SST, lon0, lat0, "sst")
  sat <- extract(SAT, lon0, lat0, "sat")
  
  write.csv(sst, paste0("cleandata/Climate/", id, "_SST.csv"), row.names = FALSE)
  write.csv(sat, paste0("cleandata/Climate/", id, "_SAT.csv"), row.names = FALSE)
  
  cat(sprintf("%-8s SST cell %.0f km away, SAT cell %.0f km away\n",
              id, sst$dist_km[1], sat$dist_km[1]))
}

## plot
cores  <- sites$CoreID
cutoff <- 11000   # set to 11 if `age` is in ka rather than years

dir.create("figures/Climate/", recursive = TRUE, showWarnings = FALSE)

plot_core <- function(sst, sat, id, xmax = NULL) {
  
  if (!is.null(xmax)) {
    sst <- sst[sst$age <= xmax, ]
    sat <- sat[sat$age <= xmax, ]
  }
  
  yr <- range(sst$sst_lwr, sst$sst_upr, sat$sat_lwr, sat$sat_upr, na.rm = TRUE)
  xr <- rev(range(c(sst$age, sat$age), na.rm = TRUE))
  
  ttl <- if (is.null(xmax)) id else paste0(id, " (", xmax, " BP to present)")
  
  plot(sst$age, sst$sst, type = "l", lwd = 2, col = "dodgerblue",
       xlim = xr, ylim = yr,
       xlab = "Age (years BP)", ylab = "Temperature (\u00b0C)", main = ttl)
  lines(sst$age, sst$sst_lwr, lty = 3, col = "dodgerblue")
  lines(sst$age, sst$sst_upr, lty = 3, col = "dodgerblue")
  
  lines(sat$age, sat$sat, lwd = 2, col = "red4")
  lines(sat$age, sat$sat_lwr, lty = 3, col = "red4")
  lines(sat$age, sat$sat_upr, lty = 3, col = "red4")
  
  legend("bottomright", c("SST", "SAT"), lwd = 2, bty = "n",
         col = c("dodgerblue", "red4"))
}

for (id in cores) {
  
  sst <- read.csv(paste0("cleandata/Climate/", id, "_SST.csv"))
  sat <- read.csv(paste0("cleandata/Climate/", id, "_SAT.csv"))
  
  pdf(paste0("figures/Climate/", id, "_climate.pdf"), width = 7, height = 8)
  par(mfrow = c(2, 1), mar = c(4, 4, 3, 1))
  plot_core(sst, sat, id)
  plot_core(sst, sat, id, xmax = cutoff)
  dev.off()
}

###############################################################################
## Which grid cells were actually used?
## Table + map of the requested sites against the LGMR cells they mapped to.
## Runs after the extraction section (needs SST, SAT and gc_dist).
###############################################################################

wrap180 <- function(x) ((x + 180) %% 360) - 180   # 0-360 -> -180-180 for plotting

## the nearest cell with data, coordinates only
nearest <- function(f, lon0, lat0) {
  d <- as.vector(gc_dist(lon0, lat0, f$lon, f$lat))
  d[!f$wet] <- NA
  i <- which.min(d)
  list(lon = wrap180(f$lon[i]), lat = f$lat[i], dist = d[i])
}

used <- do.call(rbind, lapply(seq_len(nrow(sites)), function(i) {
  s <- nearest(SST, sites$lon[i], sites$lat[i])
  a <- nearest(SAT, sites$lon[i], sites$lat[i])
  data.frame(CoreID   = sites$CoreID[i],
             site_lon = sites$lon[i], site_lat = sites$lat[i],
             sst_lon  = s$lon, sst_lat = s$lat, sst_km = round(s$dist, 1),
             sat_lon  = a$lon, sat_lat = a$lat, sat_km = round(a$dist, 1),
             stringsAsFactors = FALSE)
}))

write.csv(used, "rawdata/Climate/site_grid_cells.csv", row.names = FALSE)
print(used, row.names = FALSE)



## every ocean cell of the SST grid in view - gives the model's coastline
glon <- wrap180(as.vector(SST$lon))
glat <- as.vector(SST$lat)
k <- which(SST$wet & glon > xlim[1] & glon < xlim[2] &
             glat > ylim[1] & glat < ylim[2])

gused <- data.frame("lat"= glat[k],"lon"= glon[k])

plot(glon[k], glat[k], pch = 15, cex = 0.6, col = "grey85",
     xlim = xlim, ylim = ylim, asp = 1 / cos(mean(ylim) * pi / 180),
     xlab = "Longitude", ylab = "Latitude",
     main = "Requested sites vs LGMR cells used")

segments(used$site_lon, used$site_lat, used$sst_lon, used$sst_lat, col = "dodgerblue")
segments(used$site_lon, used$site_lat, used$sat_lon, used$sat_lat, col = "red4")

points(used$sst_lon,  used$sst_lat,  pch = 17, col = "dodgerblue", cex = 1.2)
points(used$sat_lon,  used$sat_lat,  pch = 15, col = "red4",       cex = 1.2)
points(used$site_lon, used$site_lat, pch = 19, col = "black",      cex = 1.1)
text(used$site_lon, used$site_lat, used$CoreID, pos = 3, cex = 0.7)

legend("topleft", bty = "n", pt.cex = 1.1,
       pch = c(19, 17, 15, 15), cex = 0.8,
       col = c("black", "dodgerblue", "red4", "grey85"),
       legend = c("requested site", "SST cell used", "SAT cell used",
                  "ocean cells (SST grid)"))


## ---- map --------------------------------------------------------------------
xlim <- range(c(used$site_lon, used$sst_lon, used$sat_lon)) + c(-2, 2)
ylim <- range(c(used$site_lat, used$sst_lat, used$sat_lat)) + c(-1, 1)

basemap(limits = c(xlim , ylim ), bathy.style = "rcb", rotate = TRUE, glaciers = TRUE) +
  ggspatial::geom_spatial_segment(
    aes(x = site_lon, y = site_lat, xend = sst_lon, yend = sst_lat),
    data = used, crs = 4326, great_circle = FALSE, detail = 100,
    color = "grey30", lwd = 0.5)+
  ggspatial::geom_spatial_segment(
    aes(x = site_lon, y = site_lat, xend = sat_lon, yend = sat_lat),
    data = used, crs = 4326, great_circle = FALSE, detail = 100,
    color = "grey30", lwd = 0.5)+
  ggspatial::geom_spatial_point(aes(x = lon, y = lat), data = sites, color = "black", cex = 3.5, crs = 4326) +
  ggspatial::geom_spatial_point(aes(x = lon, y = lat), data = sites,
                                color = c("#1F5C8B", "#C1663A", "#4E8C57"),
                                cex = 2.5, crs = 4326) +
  ggspatial::geom_spatial_point(aes(x = sst_lon, y = sst_lat), data = used, color = "darkblue", cex = 3.5, crs = 4326) +
  ggspatial::geom_spatial_point(aes(x = sat_lon, y = sat_lat), data = used, color = "darkred", cex = 3.5, crs = 4326) +
  
    xlab("") + ylab("") +
  annotation_scale(location = "br") +
  annotation_north_arrow(location = "tr", which_north = "true",
                         height = unit(0.9, "cm"), width = unit(0.9, "cm")) 





###############################################################################
##  Monte Carlo trend through the Rockall Trough cold-water coral eNd record
##
##  Draws n_iter simulated versions of the record. In each one, every coral is
##  jittered in age and in eNd according to its own uncertainty, and a smooth
##  curve is fitted through the jittered points. The spread of those curves is
##  the uncertainty on the trend.
##
##  Base R only - no packages required.
##  Input : coral_eNd_plot_data.csv
##  Runtime: ~1 min for 10000 iterations.
###############################################################################

## ---- settings ---------------------------------------------------------
csv_in      <- "rawdata/Climate/Coral_eNd_all.csv"
png_out     <- "figures/Climate/coral_eNd_montecarlo.png"
curve_out   <- "cleandata/Climate/coral_eNd_curve.csv"

n_iter      <- 10000   # number of simulated records
method      <- "spline"# "kernel" (Gaussian, base R ksmooth) or "spline" (smooth.spline)
bandwidth   <- 0.40    # kernel width, ka. Larger = smoother
spline_df   <- 14      # only used when method = "spline"
sigma_level <- 2       # published errors are 2 sigma; divided by this before sampling
bootstrap   <- FALSE    # also resample corals with replacement (see note at the bottom)
max_gap     <- 0.50    # blank the curve where the nearest coral is further away, ka
grid_n      <- 400     # resolution of the output curve
n_show      <- 300     # how many individual routes to draw on the figure
drop_group  <- "Porcupine Seabight"   # site_group values to exclude
seed        <- 42

## ---- load and filter --------------------------------------------------
istrue <- function(v) !is.na(v) & as.character(v) %in% c("TRUE", "true", "T", "1")

d <- read.csv(csv_in, stringsAsFactors = FALSE)
d <- d[istrue(d$use_for_plot) & !(d$site_group %in% drop_group), ]
d <- d[!is.na(d$age_ka_BP) & !is.na(d$eNd) &
         !is.na(d$age_err_ka) & !is.na(d$eNd_err), ]

x  <- d$age_ka_BP
y  <- d$eNd
sx <- d$age_err_ka / sigma_level   # 1 sigma
sy <- d$eNd_err    / sigma_level   # 1 sigma
n  <- length(x)

cat(sprintf("%d corals, %.2f to %.2f ka BP\n", n, min(x), max(x)))
cat(sprintf("1 sigma: age median %.3f ka (max %.3f), eNd median %.2f (max %.2f)\n",
            median(sx), max(sx), median(sy), max(sy)))

## ---- Monte Carlo ------------------------------------------------------
grid <- seq(min(x), max(x), length.out = grid_n)
fits <- matrix(NA_real_, nrow = n_iter, ncol = grid_n)

set.seed(seed)
pb <- txtProgressBar(min = 0, max = n_iter, style = 3)
for (i in seq_len(n_iter)) {
  
  k  <- if (bootstrap) sample.int(n, n, replace = TRUE) else seq_len(n)
  xi <- rnorm(n, x[k], sx[k])          # jitter the ages
  yi <- rnorm(n, y[k], sy[k])          # jitter the eNd values
  o  <- order(xi); xi <- xi[o]; yi <- yi[o]
  
  s <- if (method == "kernel") {
    ksmooth(xi, yi, kernel = "normal", bandwidth = bandwidth, x.points = grid)$y
  } else {
    f <- try(smooth.spline(xi, yi, df = spline_df), silent = TRUE)
    if (inherits(f, "try-error")) rep(NA_real_, grid_n) else predict(f, grid)$y
  }
  
  ## blank the curve wherever it would be interpolating across a real gap
  j    <- findInterval(grid, xi)
  dlo  <- ifelse(j == 0, Inf, abs(grid - xi[pmax(j, 1)]))
  dhi  <- ifelse(j == n, Inf, abs(xi[pmin(j + 1, n)] - grid))
  s[pmin(dlo, dhi) > max_gap] <- NA
  
  fits[i, ] <- s
  if (i %% 200 == 0) setTxtProgressBar(pb, i)
}
close(pb)

## ---- summarise --------------------------------------------------------
qs   <- function(p) apply(fits, 2, quantile, probs = p, na.rm = TRUE, names = FALSE)
med  <- apply(fits, 2, median, na.rm = TRUE)
mn   <- apply(fits, 2, mean,   na.rm = TRUE)
lo95 <- qs(0.025); hi95 <- qs(0.975)
lo68 <- qs(0.160); hi68 <- qs(0.840)
nval <- colSums(!is.na(fits))

out <- data.frame(age_ka_BP = grid, median = med, mean = mn,
                  p2.5 = lo95, p16 = lo68, p84 = hi68, p97.5 = hi95,
                  n_valid = nval)
write.csv(out, curve_out, row.names = FALSE)

cat(sprintf("\nmedian curve spans %.2f to %.2f eNd\n", min(med, na.rm=TRUE), max(med, na.rm=TRUE)))
cat(sprintf("mean 95%% band width: %.2f epsilon units\n", mean(hi95 - lo95, na.rm = TRUE)))
res <- y - approx(grid, med, x, rule = 2)$y
cat(sprintf("scatter of the corals about the median curve: RMS %.2f epsilon units\n", sd(res)))

## ---- figure -----------------------------------------------------------
runs <- function(ok) {                       # contiguous non-NA stretches
  r <- rle(ok); e <- cumsum(r$lengths)
  cbind(e - r$lengths + 1, e)[r$values, , drop = FALSE]
}

src  <- factor(d$source)
pal  <- c("#1F5C8B", "#C1663A", "#4E8C57")[as.integer(src)]

png(png_out, width = 2400, height = 1400, res = 220)
par(mar = c(4.4, 4.6, 3.0, 1.2), mgp = c(2.7, 0.7, 0), las = 1, bty = "n")

plot(NA, xlim = range(grid), ylim = range(c(y - d$eNd_err, y + d$eNd_err)),
     xlab = "Age (cal ka BP)", ylab = expression(epsilon[Nd]),
     main = "SW Rockall Trough cold-water corals: Monte Carlo trend in eNd")
#grid(col = "grey92", lty = 1)

## a subset of the individual routes
for (i in sample.int(n_iter, min(n_show, n_iter)))
  lines(grid, fits[i, ], col = adjustcolor("grey35", alpha.f = 0.035), lwd = 1)

## credible bands
for (r in seq_len(nrow(runs(!is.na(med))))) {
  ix <- runs(!is.na(med))[r, 1]:runs(!is.na(med))[r, 2]
  polygon(c(grid[ix], rev(grid[ix])), c(lo95[ix], rev(hi95[ix])),
          col = adjustcolor("#1F5C8B", 0.16), border = NA)
#  polygon(c(grid[ix], rev(grid[ix])), c(lo68[ix], rev(hi68[ix])),
 #         col = adjustcolor("#1F5C8B", 0.26), border = NA)
}

## data with 2 sigma error bars in both axes
segments(x - d$age_err_ka, y, x + d$age_err_ka, y, col = adjustcolor(pal, 0.55))
segments(x, y - d$eNd_err, x, y + d$eNd_err, col = adjustcolor(pal, 0.55))
points(x, y, pch = 21, bg = adjustcolor(pal, 0.85), col = "white", cex = 0.95, lwd = 0.6)

lines(grid, med, col = "#0B2E4A", lwd = 3)
#lines(grid, mn,  col = "#0B2E4A", lwd = 1.4, lty = 3)

legend("bottomright", bty = "n", cex = 0.8,
       legend = c(levels(src), "median of 10000 routes", "95% of routes"),
       pch = c(rep(21, nlevels(src)), NA, 22),
       pt.bg = c(c("#1F5C8B", "#C1663A", "#4E8C57")[seq_len(nlevels(src))], NA,
                 adjustcolor("#1F5C8B", 0.26)),
       lty = c(rep(NA, nlevels(src)), 1, NA),
       lwd = c(rep(NA, nlevels(src)), 3, NA),
       col = c(rep("white", nlevels(src)), "#0B2E4A", NA))
dev.off()
cat(sprintf("\nwrote %s and %s\n", png_out, curve_out))

##PALTIDE

pal_SF <- read.csv("rawdata/Climate/Paltide_RSL_SF.csv")
pal_FB <- read.csv("rawdata/Climate/Paltide_RSL_FB.csv")
pal_MB <- read.csv("rawdata/Climate/Paltide_RSL_MB.csv")

#where are we 
str(pal_MB)

pal <- list(PC012 = pal_FB, PC014 = pal_MB, GC017 = pal_SF)

## unique grid location(s) per core
pal_loc <- do.call(rbind, lapply(names(pal), function(id) {
  d <- pal[[id]]
  u <- unique(d[, c("longitude", "latitude")])
  i <- match(id, sites$CoreID)
  data.frame(CoreID   = id,
             site_lon = sites$lon[i], site_lat = sites$lat[i],
             pal_lon  = u$longitude,  pal_lat  = u$latitude,
             stringsAsFactors = FALSE)
}))

pal_loc$dist_km <- round(gc_dist(pal_loc$site_lon, pal_loc$site_lat,
                                 pal_loc$pal_lon,  pal_loc$pal_lat), 1)
pal_loc$loc <- ave(seq_len(nrow(pal_loc)), pal_loc$CoreID, FUN = seq_along)

write.csv(pal_loc, "rawdata/Climate/site_paltide_cells.csv", row.names = FALSE)
pxlim <- range(c(pal_loc$site_lon, pal_loc$pal_lon)) + c(-1, 1)
pylim <- range(c(pal_loc$site_lat, pal_loc$pal_lat)) + c(-0.5, 0.5)


pdf("figures/Climate/paltide_sitesused.pdf",height=5,width = 4)
basemap(limits = c(pxlim, pylim), bathy.style = "rcb",
        rotate = TRUE, glaciers = TRUE) +
  ggspatial::geom_spatial_segment(
    aes(x = site_lon, y = site_lat, xend = pal_lon, yend = pal_lat),
    data = pal_loc, crs = 4326, great_circle = FALSE,
    color = "grey30", lwd = 0.5) +
  ggspatial::geom_spatial_point(aes(x = pal_lon, y = pal_lat), data = pal_loc,
                                color = "darkorange3", cex = 3, crs = 4326) +
  ggspatial::geom_spatial_point(aes(x = lon, y = lat), data = sites,
                                color = "black", cex = 3.5, crs = 4326) +
  ggspatial::geom_spatial_point(aes(x = lon, y = lat), data = sites,
                                color = c("#1F5C8B", "#C1663A", "#4E8C57"),
                                cex = 2.5, crs = 4326) +
  xlab("") + ylab("") +
  annotation_scale(location = "br") +
  annotation_north_arrow(location = "tr", which_north = "true",
                         height = unit(0.9, "cm"), width = unit(0.9, "cm"))
dev.off()


