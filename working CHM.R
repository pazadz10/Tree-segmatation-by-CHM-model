
library(lidR)
library(raster)

# Read the LAS file
las <- readLAS("C:\\Users\\pazad\\Documents\\lac_230907_l1_modified.las")
# Define CSF parameters
mycsf <- csf(
  sloop_smooth = TRUE,         # Enable slope smoothing
  class_threshold = 0.5,       # Distance threshold for classification
  cloth_resolution = 0.5,      # Distance between particles in the cloth
  rigidness = 1L,              # Rigidness of the cloth (1 = very soft)
  iterations = 500L,           # Maximum iterations for cloth simulation
  time_step = 0.65             # Time step for cloth simulation under gravity
)

# Run CSF classification
las <- classify_ground(las, mycsf, last_returns = TRUE)

# Filter ground points
las_ground <- filter_ground(las)

# Normalize height
las_normalized <- normalize_height(las_ground, knnidw())
# Downsample the LAS file
las_downsampled <- decimate_points(las, homogenize(3))  # Adjust the factor as needed



# Rasterize the canopy using the downsampled data
chm <- rasterize_canopy(las_downsampled, 0.25, pitfree(subcircle = 0.1))

col <- height.colors(50)


# Visualize the result
plot(chm, col = col)

# Generate kernel and smooth chm
kernel <- matrix(1, 3, 3)
schm <- terra::focal(x = chm, w = kernel, fun = median, na.rm = TRUE)
plot(schm, col = col)

ttops <- locate_trees(las = schm, algorithm = lmf(ws = 2.5))
ttops
plot(chm, col = col)
plot(ttops, col = "black", add = TRUE, cex = 0.5)


algo <- dalponte2016(chm, ttops)
las_segmented <- segment_trees(las_downsampled, algo) # segment point cloud
plot(las_segmented, bg = "white", size = 4, color = "treeID") # visualize trees

# Calculate crown metrics
crowns <- crown_metrics(las_segmented, func = .stdtreemetrics, geom = "convex")
plot(crowns["convhull_area"], main = "Crown area (convex hull)")

# Calculate additional crown metrics
crowns_dalponte <- crown_metrics(las_segmented, func = .stdtreemetrics, geom = "concave")
crowns_li <- crown_metrics(las_segmented, func = .stdtreemetrics, geom = "concave")

# Visualize crown metrics
par(mfrow=c(1,2),mar=rep(0,4))
plot(sf::st_geometry(crowns_dalponte), reset = FALSE)
plot(sf::st_geometry(crowns_li), reset = FALSE)

# Export segmented tree crowns to GeoPackage
output_file <- "C:\\Users\\pazad\\Documents\\tree_crowns2.gpkg"
sf::st_write(crowns, output_file, layer = "tree_crowns2", delete_layer = TRUE)