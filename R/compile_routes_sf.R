shps <- list.files("../Geodata/_MABM_canonical_routes",
                   pattern = ".shp$", full.names = TRUE)

for (i in seq_along(shps)) {
    rte_nm <- gsub("_.*$","",basename(shps[i]))
    tmp <- st_read(shps[i], quiet = TRUE) %>%
        select(geometry) %>%
        dplyr::mutate(site = rte_nm) %>%
        group_by(site) %>%
        summarise()
    if (i == 1)
        sites <- tmp
    else
        sites <- rbind(sites, tmp)
}
