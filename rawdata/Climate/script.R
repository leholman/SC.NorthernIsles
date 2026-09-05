###### SST recovery
library(ncdf4)
dat <- nc_open("LGMR_SST_climo.nc")

print(dat)
names(dat$var)      # variable names
names(dat$dim)      # dimension names
sst <- ncvar_get(dat, "sst")          # 4D array: lon × lat × age 
sst_std <- ncvar_get(dat, "sst_std")

dim(sst)
dim(sst_std)

lon <- ncvar_get(dat, "lon")
lat <- ncvar_get(dat, "lat")
age <- ncvar_get(dat, "age")
ncatt_get(dat, "age")

## find the data - 320 384 120 lon lat time

# target - Kalø Vig
lon0 <- 10.397342
lat0 <- 56.249136



lon0n <- ifelse(lon0 < 0, lon0 + 360, lon0)

# lon
ilon <- which.min(abs(((lon[,1] - (lon0n)) )))
#nearest lon 
lon[ilon,1]

# lat 
ilat <- which.min(abs(lat[1,] - lat0))
#nearest lat
lat[1,ilat]


# extract SST timeseries
sst_ts <- sst[ilon,ilat,]
sst_std <- sst_std[ilon,ilat,]


data <- data.frame("age"=age,"sst"=sst_ts ,"sstSD" =sst_std)
data$CIupp <- data$sst + (2*data$sstSD)
data$CIlwr <- data$sst - (2*data$sstSD)

##### NOW SAT
dat <- nc_open("LGMR_SAT_climo.nc")

print(dat)
names(dat$var)      # variable names
names(dat$dim)      # dimension names
sat <- ncvar_get(dat, "sat")          # 4D array: lon × lat × age 
sat_std <- ncvar_get(dat, "sat_std")

dim(sat)
dim(sat_std)

lon <- ncvar_get(dat, "lon")
lat <- ncvar_get(dat, "lat")
age <- ncvar_get(dat, "age")
ncatt_get(dat, "age")

## find the data

# target - uses same lon0/lat0 as above

lon0n <- ifelse(lon0 < 0, lon0 + 360, lon0)

# lon
ilon <- which.min(abs(((lon - (lon0n)) )))
#nearest lon 
lon[ilon]

# lat 
ilat <- which.min(abs(lat - lat0))
#nearest lat
lat[ilat]

# extract SAT timeseries
sat_ts <- sat[ilon,ilat,]
sat_std <- sat_std[ilon,ilat,]

dat <- nc_open("LGMR_SAT_climo.nc")

print(dat)
names(dat$var)      # variable names
names(dat$dim)      # dimension names
sat <- ncvar_get(dat, "sat")          # 4D array: lon × lat × age × nEns
sat_std <- ncvar_get(dat, "sat_std")
dim(sat)


lon <- ncvar_get(dat, "lon")
lat <- ncvar_get(dat, "lat")
age <- ncvar_get(dat, "age")

## find the data

# target
lon0 <- 10.397342
lat0 <- 56.249136



ilon <- which.min(abs(((lon - (lon0)) )))
lon[ilon]
ilat <- which.min(abs(lat - lat0))
lat[ilat]




# extract SAT timeseries
sat_ts <- sat[ilon,ilat,]
sat_std <- sat_std[ilon,ilat,]


data$sat <- sat_ts 
data$sat_upr <- sat_ts + (2*sat_std) 
data$sat_lwr <- sat_ts - (2*sat_std) 



# plot (index vs age labels)
plot(data$age,data$sst, pch=16, xlab="n", ylab="SAT (°C)",ylim=c(0,26),cex=0)
lines(data$age,data$sst, lty=1,col="dodgerblue")
lines(data$age,data$CIupp, lty=2,col="dodgerblue")
lines(data$age,data$CIlwr, lty=2,col="dodgerblue")

lines(data$age,data$sat, lty=1,col="red4")
lines(data$age,data$sat_upr, lty=2,col="red4")
lines(data$age,data$sat_lwr, lty=2,col="red4")










