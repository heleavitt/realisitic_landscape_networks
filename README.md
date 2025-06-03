realistic_landscape_networks
Author: Herbert Leavitt
Last Updated: [4/15/25]

Overview:
This repository simulates aquatic movement paths (random walks) of nekton across coastal landscapes while remaining constrained to water bodies. It calculates the probabilities of reaching specific habitat targets and estimates the composition of habitat types encountered along each simulated path.

Core Features:
- Reads spatial data including points, polygons, and classified rasters.
- Identifies closest water boundaries to sample sites.
- Rasterizes habitat shapefiles to create masks for constrained movement.
- Implements a correlated random walk simulation bounded by a water mask.
- Computes the likelihood of reaching targets and the composition of traversed habitats.
- Visualizes all random walk paths and distance distributions.

Usage:
The main R script:
- Loads site data and habitat polygons.
- Converts polygons to raster and masks non-water areas.
- Runs a parameterized random walk with options for step size, autocorrelation (gamma), and noise (alpha).
- Produces spatial plots and summary statistics of habitat use and distances.

Requirements (R packages):
- sf
- ggplot2
- raster
- fasterize
- mvtnorm

Inputs:
- Shapefiles: water polygons (e.g., `flight9_waterpoly.shp`) and classified habitat polygons (`Flight_9.shp`)
- CSV: site coordinates (e.g., `drop_field2308.csv`)

Outputs:
- Reached target probabilities
- Composition of habitat types along paths
- Maps of walker tracks and target sites
- Histogram of distances from start
