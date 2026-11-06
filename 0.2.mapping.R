### Mapping
library(ggOceanMaps)
library(sf)
library(ggspatial)


sites <- read.csv("metadata_site.csv")
sites$side <- c("right", "left", "left")
sites$hj    <- ifelse(sites$side == "left", 1, 0)
sites$lon_l <- sites$lon + ifelse(sites$side == "left", -0.12, 0.12)
sites$lat_l <- sites$lat + 0.035
sites$lat_n <- sites$lat - 0.035

pdf("mapping/mapv1.pdf",width = 4,height = 3.5)
basemap(limits = c(xlim , ylim ), bathy.style = "rcb", rotate = TRUE, glaciers = TRUE) +
  ggspatial::geom_spatial_point(aes(x = lon, y = lat), data = sites, color = "black", cex = 3.5, crs = 4326) +
  ggspatial::geom_spatial_point(aes(x = lon, y = lat), data = sites,
                                color = c("#1F5C8B", "#C1663A", "#4E8C57"),
                                cex = 2.5, crs = 4326)+ 
  ggspatial::geom_spatial_text(aes(x = lon_l, y = lat_l, label = CoreID, hjust = hj),
                               data = sites, crs = 4326,
                               color = "black", fontface = "bold", size = 3.2) +
  #ggspatial::geom_spatial_text(aes(x = lon_l, y = lat_n, label = name, hjust = hj),
   #                            data = sites, crs = 4326,
    #                          color = "grey35", fontface = "italic", size = 2.8) +
  xlab("") + ylab("") +
  annotation_scale(location = "br") +
  annotation_north_arrow(location = "tr", which_north = "true",
                         height = unit(0.9, "cm"), width = unit(0.9, "cm")) 

dev.off()
