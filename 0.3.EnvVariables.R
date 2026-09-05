###############################################################################
## ENV data output
###############################################################################

## Libraries 

library(ncdf4)


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
  
  write.csv(sst, paste0("cleandata//Climate/", id, "_SST.csv"), row.names = FALSE)
  write.csv(sat, paste0("cleandata//Climate/", id, "_SAT.csv"), row.names = FALSE)
  
  cat(sprintf("%-8s SST cell %.0f km away, SAT cell %.0f km away\n",
              id, sst$dist_km[1], sat$dist_km[1]))
}

## optional check: plot the last site extracted
plot(sst$age, sst$sst, type = "l", lwd = 2, col = "dodgerblue",
     xlim = rev(range(sst$age)), ylim = range(sst$sst_lwr, sat$sat_upr),
     xlab = "Age (years BP)", ylab = "Temperature (\u00b0C)", main = id)
lines(sat$age, sat$sat, lwd = 2, col = "red4")
legend("bottomright", c("SST", "SAT"), lwd = 2, bty = "n",
       col = c("dodgerblue", "red4"))

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

write.csv(used, "output/Climate/site_grid_cells.csv", row.names = FALSE)
print(used, row.names = FALSE)


## ---- map --------------------------------------------------------------------
xlim <- range(c(used$site_lon, used$sst_lon, used$sat_lon)) + c(-2, 2)
ylim <- range(c(used$site_lat, used$sst_lat, used$sat_lat)) + c(-1, 1)

## every ocean cell of the SST grid in view - gives the model's coastline
glon <- wrap180(as.vector(SST$lon))
glat <- as.vector(SST$lat)
k <- which(SST$wet & glon > xlim[1] & glon < xlim[2] &
             glat > ylim[1] & glat < ylim[2])

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


