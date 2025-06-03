library(sf)
library(ggplot2)
library(raster)

pffw_sites <- read.csv("C:/Users/hl51981/OneDrive - University of Georgia/Leavitt_Herbert/PFFW/Manuscripts/Habscale/Repository/landscape_analysis/input/drop_field2308.csv")

# Replace this with the path to your shapefile (without file extension)
shapefile_path <- "C:/Users/hl51981/OneDrive - University of Georgia/Leavitt_Herbert/PFFW/Pt Fourchon Food Webs R Directory/raw_input/shapefiles/"

# Read the shapefile (it will look for .shp, .shx, .dbf automatically)
polygon <- st_read(paste0(shapefile_path, "flight9_waterpoly.shp"))
full_polygon <- st_read("C:/Users/hl51981/OneDrive - University of Georgia/Leavitt_Herbert/PFFW/Manuscripts/Walkers/input/Flight_9.shp")



pffw_sites <- st_as_sf(pffw_sites, coords = c("longitude", "latitude"), crs = 4326)
pffw_sites <- st_transform(pffw_sites, st_crs(polygon))



# Optionally dissolve to one polygon if you don’t care which feature
polygon_union <- st_union(polygon)

# Calculate distance from each point to the polygon
distances <- st_distance(pffw_sites, polygon_union)

# Logical vector: within 1 meter
within_1m <- as.numeric(distances) <= 1

# Subset points
points_near_poly <- pffw_sites[within_1m, ]

#Get nearest point on the polygon boundary for each site
nearest_lines <- st_nearest_points(points_near_poly, polygon_union)

# Extract the endpoint on the polygon for each nearest line
nearest_coords <- lapply(nearest_lines, function(line) {
  coords <- st_coordinates(line)
  coords[nrow(coords), ]  # the point on the polygon
})

# Convert to sf points
nearest_sf <- do.call(rbind, nearest_coords) %>%
  as.data.frame() %>%
  st_as_sf(coords = c("X", "Y"), crs = st_crs(points))

# Convert to coordinates
start_site <- st_geometry(nearest_sf)[[1]]  # gives you a POINT object
start_site_coords <- c(start_site[[1]], start_site[[2]])  # x and y as numeric

library(sf)
library(fasterize)

# Assume 'polygon' is your water body sf object
polygon_union_sf <- st_sf(geometry = polygon_union)
polygon_union_sf$value <- 1  # add dummy field to rasterize

# Create a numeric ID column for class names
full_polygon$Class_id <- as.numeric(factor(full_polygon$Class_name))



# Now use this new numeric column in fasterize
full_raster <- fasterize(full_polygon, template_raster, field = "Class_id")
# Create a raster template
template_raster <- raster(extent(full_polygon), res = 0.1, crs = st_crs(full_polygon)$proj4string)

# Rasterize using the 'value' column
water_mask <- fasterize(polygon_union_sf, template_raster, field = "value")
full_raster <- fasterize(full_polygon, template_raster, field = "Class_name") 

is_within_water <- function(xy, raster_mask) {
  # Ensure input is a 1-row matrix (2 columns)
  if (length(xy) != 2 || any(!is.finite(xy))) return(FALSE)
  cell <- cellFromXY(raster_mask, matrix(xy, nrow = 1))
  !is.na(cell) && !is.na(raster_mask[cell])
}

plot_random_walk_step <- function(current_pos, step_num, targets, water_mask, window_size = 200) {
  # Define zoomed-in extent around the current position
  xlim <- c(current_pos[1] - window_size / 2, current_pos[1] + window_size / 2)
  ylim <- c(current_pos[2] - window_size / 2, current_pos[2] + window_size / 2)
  
  # Crop water mask for faster plotting (only show visible area)
  cropped_mask <- crop(water_mask, extent(xlim[1], xlim[2], ylim[1], ylim[2]))
  
  # Plot water raster first
  plot(cropped_mask, col = "lightblue", legend = FALSE, xlim = xlim, ylim = ylim,
       main = paste("Step", step_num), xlab = "X", ylab = "Y", asp = 1)
  
  # Plot current position
  points(current_pos[1], current_pos[2], pch = 19, col = "blue", cex = 1.5)
  
  # Plot targets
  points(targets[,1], targets[,2], pch = 4, col = "red", cex = 1.2)
} 

compute_target_distances <- function(pos, targets) {
  if (inherits(targets, "sf")) {
    targets <- st_coordinates(targets)
  }
  if (is.vector(targets)) {
    targets <- matrix(targets, ncol = 2, byrow = TRUE)
  }
  matrix_pos <- matrix(pos, nrow = nrow(targets), ncol = 2, byrow = TRUE)
  sqrt(rowSums((targets - matrix_pos)^2))
}

simulate_walk_correlated <- function(start, targets, water_mask, step_size = 3,
                                      max_steps = 600, n_sim = 100,
                                      gamma = 0.8, beta = 0, alpha = 15,
                                      sleep = 0.05, full_raster = NULL, buffer_meters = 50) {
    par(mar = c(3, 3, 2, 1))  # bottom, left, top, right
    if (inherits(targets, "sf")) {
    targets <- st_coordinates(targets)
  } else if (is.list(targets)) {
    targets <- do.call(rbind, targets)
  }
  
  track_list <- vector("list", n_sim)
  mean_compositions <- vector("list", n_sim)
  
  R <- gamma * matrix(c(cos(beta), sin(beta),
                        -sin(beta), cos(beta)), 2, 2)
  D <- diag(alpha, 2)
  
  reached_counts <- numeric(nrow(targets))
  pb <- txtProgressBar(min = 0, max = n_sim * max_steps, style = 3)
  
  for (sim in 1:n_sim) {
    pos <- start
    steps_taken <- 0
    v_prev <- runif(2, -1, 1)
    v_prev <- v_prev / sqrt(sum(v_prev^2))
    track_coords <- matrix(NA, nrow = max_steps + 1, ncol = 2)
    track_coords[1, ] <- start
    already_hit <- FALSE
    stuck_attempts <- 0
    max_attempts <- 10000
    
    # Initialize habitat tracking
    class_counts <- integer()
    total_pixels <- 0
    
    while (steps_taken < max_steps && stuck_attempts < max_attempts) {
      v_rotated <- R %*% v_prev
      z <- as.vector(mvtnorm::rmvnorm(1, mean = c(0, 0), sigma = D))
      v <- v_rotated + z
      v <- step_size * (v / sqrt(sum(v^2)))
      new_pos <- pos + v
      
      if (is_within_water(new_pos, water_mask)) {
        steps_taken <- steps_taken + 1
        stuck_attempts <- 0
        pos <- new_pos
        v_prev <- v
        track_coords[steps_taken + 1, ] <- new_pos
        
        # ✅ Extract and accumulate habitat composition
        if (!is.null(full_raster)) {
          vals <- raster::extract(full_raster,
                                  SpatialPoints(matrix(pos, nrow = 1), proj4string = CRS(projection(full_raster))),
                                  buffer = 5)[[1]]
          vals <- vals[!is.na(vals)]
          
          if (length(vals) > 0) {
            new_counts <- table(vals)
            for (cls in names(new_counts)) {
              if (!cls %in% names(class_counts)) {
                class_counts[cls] <- 0
              }
              class_counts[cls] <- class_counts[cls] + new_counts[[cls]]
            }
            total_pixels <- total_pixels + length(vals)
          }
        }
        
        if (!already_hit) {
          dists <- compute_target_distances(pos, targets)
          hit_targets <- which(dists <= step_size*3)
          if (length(hit_targets) > 0) {
            reached_counts[hit_targets[1]] <- reached_counts[hit_targets[1]] + 1
            already_hit <- TRUE
            cat("  🎯 Reached target", hit_targets[1], "\n")
          }
        }
      } else {
        stuck_attempts <- stuck_attempts + 1
      }
      
      setTxtProgressBar(pb, ((sim - 1) * max_steps) + steps_taken)
    }
    
    track_list[[sim]] <- track_coords[1:(steps_taken + 1), ]
    
    # ✅ Store mean composition
    if (total_pixels > 0) {
      mean_compositions[[sim]] <- class_counts / total_pixels
    } else {
      mean_compositions[[sim]] <- NA  # handle if no pixels were collected
    }
  }
  
  close(pb)
  return(list(
    probs = reached_counts / n_sim,
    tracks = track_list,
    mean_compositions = mean_compositions
  ))
}

# 3. Run and plot
set.seed(42)

targets_clean <- nearest_sf[!(round(st_coordinates(nearest_sf)[,1], 6) == round(start_site_coords[1], 6) &
                                round(st_coordinates(nearest_sf)[,2], 6) == round(start_site_coords[2], 6)), ]

result <- simulate_walk_correlated(
  start = start_site_coords,
  targets = targets_clean,
  water_mask = water_mask,
  full_raster = full_raster,       # <- your habitat raster
  buffer_meters = 5,
  step_size = 6,
  n_sim = 10,
  gamma = .9,
  beta = 0,
  alpha = 2, 
  max_steps = 500
)



tracks <- result$tracks
probs <- result$probs
probs
# Base R plot of all tracks
plot(water_mask, col = "lightblue", main = "All walker tracks", legend = FALSE)
cols <- rainbow(length(tracks))

for (i in seq_along(tracks)) {
  lines(tracks[[i]], col = adjustcolor(cols[i], alpha.f = 0.25), lwd = 2)}

points(st_coordinates(targets_clean), pch = 19, col = "black", cex = 1)
text(target_coords[,1], target_coords[,2], labels = seq_len(nrow(target_coords)),
     pos = 3, cex = 0.8, col = "black")
points(st_coordinates(start_site), pch = 17, col = "black", cex = 1)


# Convert mean compositions to data frame
comp_df <- do.call(rbind, lapply(result$mean_compositions, function(x) {
  if (is.null(x)) return(rep(NA, length(unique_vals)))
  x_vals <- as.numeric(x)
  names(x_vals) <- names(x)
  return(x_vals)
}))
comp_df <- as.data.frame(comp_df)
class_table <- data.frame(
  code = c(1, 2, 3, 4),
  name = c("Mangrove", "Marsh", "Artificial", "Water")
)
# Rename columns in comp_df
colnames(comp_df) <- class_table$name[match(colnames(comp_df), class_table$code)]

means <- colMeans(comp_df, na.rm = TRUE)
vars  <- apply(comp_df, 2, var, na.rm = TRUE)

summary_stats <- data.frame(
  Class = names(means),
  Mean  = means,
  Variance = vars
)
print(summary_stats)

# Combine all steps into one matrix
all_steps <- do.call(rbind, tracks)

# Calculate Euclidean distance from origin
distances <- sqrt((all_steps[,1] - start_site_coords[1])^2 +
                    (all_steps[,2] - start_site_coords[2])^2)

# Plot histogram
hist(distances, breaks = 50, col = "skyblue", main = "Distance from Start",
     xlab = "Distance (meters)", ylab = "Step count")


