# Import necessary libraries
import numpy as np
import richdem as rd

# Generate a random landscape based on perlin noise
dem = rd.generate_perlin_terrain(20, 20)

# Get a simple labels array indicating all the edge cells belong to the ocean
# and all the interior cells are not yet assigned to a depression
labels = rd.get_new_depression_hierarchy_labels(dem.shape)

# Generate the Depression Hierarchy
dephier, flowdirs = rd.get_depression_hierarchy(dem, labels)

# Array to hold the water depth
water_depth = rd.rdarray(
  np.zeros(dem.shape), no_data=-9999, geotransform=rd.STANDARD_GEOTRANSFORM
)

# A simple model of progressively adding water and redistributing it
for t in range(10):
  # Make it rain uniformly across the landscape
  water_depth += 0.01
  # Use Fill-Spill-Merge to redistribute water into depressions
  rd.fill_spill_merge(dem, labels, flowdirs, dephier, water_depth)
  # Display the output for inspection
  rd.rdShow(water_depth, vmin=0, vmax=0.175)