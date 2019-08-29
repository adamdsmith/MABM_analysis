pacman::p_load(raster, dplyr, ggplot2, ggmap, sf, readxl, USAboundaries)
ggmap::register_google(Sys.getenv("goog_key"))

vec_proxy.sfc <- function(x, ...) {
  unclass(x)
}

usa <- us_states() %>%
  filter(!state_abbr %in% c("AK", "HI", "PR"))

routes <- readxl::read_excel("./Output/Raw_DB_output/MABM_routes.xlsx") %>%
  select(station = ORGNAME,
         site = Site_Name,
         site_notes = Loc_Notes,
         len_mi = Rte_length,
         lat = Centroid_Y_Coord,
         lon = Centroid_X_Coord,
         state = State) %>%
  filter(!is.na(lat)|!is.na(lon)) %>%
  arrange(site)
routes_sf <- st_as_sf(routes, coords = c("lon", "lat"), crs = 4326) %>%
  st_transform(3857)

sa <- st_bbox(c(xmin = -100, xmax = -75, ymax = 41, ymin = 24), crs = st_crs(4326))
sa3857 <- st_as_sfc(sa) %>% st_transform(3857)

bm <- get_map(location = unname(sa), source = "stamen",
              zoom = 6, maptype = "terrain-background")

# overwrite the bbox of the ggmap object with that from transformed study area
attr(bm, "bb")$ll.lat <- st_bbox(sa3857)["ymin"]
attr(bm, "bb")$ll.lon <- st_bbox(sa3857)["xmin"]
attr(bm, "bb")$ur.lat <- st_bbox(sa3857)["ymax"]
attr(bm, "bb")$ur.lon <- st_bbox(sa3857)["xmax"]

p <- ggmap(bm) +
  geom_sf(data = usa, lwd = 1, fill = NA, color = "black", inherit.aes = FALSE) +
  geom_sf(data = routes_sf, size = 1, 
          pch = 21, fill = "red", inherit.aes = FALSE) +
  coord_sf(crs = st_crs(3857),
           xlim = st_bbox(sa3857)[c(1,3)], ylim = st_bbox(sa3857)[c(2,4)]) +
  labs(x = NULL, y = NULL) +
  theme_bw()

inset <- ggplot() + 
  geom_sf(data = sa3857) +
  geom_sf(data = usa, fill = NA) +
  coord_sf() + theme_bw() +
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = NA, color = NA))

p +
  annotation_custom(
    grob = ggplotGrob(inset),
    xmin = st_bbox(sa3857)["xmin"],
    xmax = st_bbox(sa3857)["xmin"] + abs(diff(c(st_bbox(sa3857)["xmax"], st_bbox(sa3857)["xmin"])))/3,
    ymin = 2550000,
    # ymin = st_bbox(sa3857)["ymin"],
    ymax = st_bbox(sa3857)["ymin"] + abs(diff(c(st_bbox(sa3857)["ymin"], st_bbox(sa3857)["ymax"])))/3)

ggsave("Output/MABM_route_centroids.png", dpi = 600, width = 6.5, height = 4.5)

ggplot(routes) +
  ggspatial::annotation_map_tile(zoom = 5) +
  geom_sf(inherit.aes = FALSE) +
  coord_sf(datum = NA) +
  theme_minimal()
