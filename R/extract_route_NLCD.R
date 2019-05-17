# 2016 NLCD raster is too large to include in this directory
# It can be downloaded here: https://www.mrlc.gov/data/nlcd-2016-land-cover-conus 
# Modify the path to fit your needs...
nlcd <- raster("../Geodata/NLCD_2016/NLCD_2016_Land_Cover_L48_20190424.img")

# Compile separate shapefiles of MABM routes, convert to NCLD projection, and save
# source("./R/compile_routes_sf.R")
# sites <- st_transform(sites, proj4string(nlcd))
# saveRDS(sites, file = "./Output/compiled_routes.rds")
sites <- readRDS("./Output/compiled_routes.rds")

# Function to generate Gaussian smoothing kernel weights at different sigmas
# Based on Chandler and Hepinstall-Cymermann (doi:10.1007/s10980-016-0380-z)
# Note we do not ignore habitat that falls right on the transect
gauss_kern_w <- function(D, sigma, normalize = TRUE) {
    w <- exp(-D^2 / (2*sigma^2))
    if (normalize) w <- w/sum(w, na.rm = TRUE)
    # ## Ignore habitat in focal site?
    # w[which(D == 0)] <- 0
    return(w)
}

n_sites <- nrow(sites)
vars <- paste(rep(c("wood", "upl_wood", "upl_dec", "upl_con", "upl_mix", 
                    "wet_wood", "water", "urban"), each = 3),
              c(1500, 250, 100), sep = "_")

# Create dataframe to receive zonal stats from loop
MABM_site_env_data <- data.frame(matrix(NA, nrow = n_sites, ncol = length(vars) + 1))
names(MABM_site_env_data) <- c("site", vars)

pb <- txtProgressBar(min = 0, max = n_sites, style = 3)
for (i in seq_len(n_sites)) {
    rte <- sites[i, ]
    rte_name <- rte$site
    # message(rte_name)
    rte_b <- st_buffer(rte, dist = 4500)
    
    # Distance to route grid
    rte_nlcd <- crop(nlcd, as(rte_b, "Spatial"))
    rte_g <- rasterize(as(rte, "Spatial"), rte_nlcd) 
    rte_dist <- distance(rte_g)
    
    # Set up Gaussian kernel-weighted grids at three SD (100, 250, & 1500 m)
    # Truncate to 3 SD and trim to reduce computations later

    rte_1500 <- rte_dist
    rte_1500[rte_dist[] > 4500] <- NA
    rte_1500[] <- gauss_kern_w(rte_1500[], 1500)

    rte_250 <- rte_dist
    rte_250[rte_dist[] > 750] <- NA
    rte_250[] <- gauss_kern_w(rte_250[], 250)

    rte_100 <- rte_dist
    rte_100[rte_dist[] > 300] <- NA
    rte_100[] <- gauss_kern_w(rte_100[], 100)
    
    all <- brick(rte_1500, rte_250, rte_100)
    
    # Calculate various derived products from NLCD to be weighted-averaged
    # based on the kernel-weighted smoothing grids
    
    ## % upland deciduous forest
    upl_dec <- rte_nlcd
    upl_dec[] <- ifelse(upl_dec[] == 41, 1, 0)
    p_upl_dec <- round(cellStats(upl_dec * all, sum), 3)
    names(p_upl_dec) <- paste("upl_dec", c(1500, 250, 100), sep = "_")
    
    ## % upland coniferous forest
    upl_con <- rte_nlcd
    upl_con[] <- ifelse(upl_con[] == 42, 1, 0)
    p_upl_con <- round(cellStats(upl_con * all, sum), 3)
    names(p_upl_con) <- paste("upl_con", c(1500, 250, 100), sep = "_")
    
    ## % upland mixed forest
    upl_mix <- rte_nlcd
    upl_mix[] <- ifelse(upl_mix[] == 43, 1, 0)
    p_upl_mix <- round(cellStats(upl_mix * all, sum), 3)
    names(p_upl_mix) <- paste("upl_mix", c(1500, 250, 100), sep = "_")
    
    ## % upland forest
    p_upl_wood <- round(p_upl_dec + p_upl_con + p_upl_mix, 3)
    names(p_upl_wood) <- paste("upl_wood", c(1500, 250, 100), sep = "_")
    
    ## % woody wetland
    wet_wood <- rte_nlcd
    wet_wood[] <- ifelse(wet_wood[] == 90, 1, 0)
    p_wet_wood <- round(cellStats(wet_wood * all, sum), 3)
    names(p_wet_wood) <- paste("wet_wood", c(1500, 250, 100), sep = "_")
    
    ## % forested
    p_wood <- round(p_upl_wood + p_wet_wood, 3)
    names(p_wood) <- paste("wood", c(1500, 250, 100), sep = "_")
    
    ## % "open" water, includes emergent wetland
    water <- rte_nlcd
    water[] <- ifelse(water[] %in% c(11, 95), 1, 0)
    p_water <- round(cellStats(water * all, sum), 3)
    names(p_water) <- paste("water", c(1500, 250, 100), sep = "_")
    
    ## % developed
    urban <- rte_nlcd
    urban[] <- ifelse(urban[] %in% 21:24, 1, 0)
    p_urban <- round(cellStats(urban * all, sum), 3)
    names(p_urban) <- paste("urban", c(1500, 250, 100), sep = "_")
    
    MABM_site_env_data[i, ] <- c(rte_name, p_wood, p_upl_wood, p_upl_dec, p_upl_con, p_upl_mix,
                            p_wet_wood, p_water, p_urban)
    setTxtProgressBar(pb, i)
}
close(pb)

MABM_site_env_data <- MABM_site_env_data %>%
    mutate_at(vars(wood_1500:urban_100), as.numeric)
saveRDS(MABM_site_env_data, file = "./Output/MABM_site_env_data.rds")
