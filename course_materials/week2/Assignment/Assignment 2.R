
# 1.  Try making a map of infection prevalence at admin 2 level (shapefile found in folder)
# Hint be careful - not all admin 2 areas have observations.
 

# 2. Use this file of population density to calculate a table of 
# positives per admin 2 level
BF_pop <- raster("BF_pop.tif") # originally from worlpop.org

# 3. use the getData function to get hold of more bioclimatic variables.
#    Add 2-3 more variables to each survey point and explore relationships with plots
bioclim_global <- raster::getData('worldclim', var='bio', res=10) # low res but easy to work with for the moment
bioclim_global # is a raster stack. Bascially a cube of rasters
plot(bioclim_global)
# See http://www.worldclim.org/bioclim for more details on these layers
# to get the first raster 'bio1'
bio1 <- bioclim_global[["bio1"]] # bioclim_global[[1]] will also work
plot(bioclim_global[["bio1"]])

# 4. reclassify the BF_land_use layer according to 'LCCS Entry' column at the bottom of the
# doc http://due.esrin.esa.int/files/GLOBCOVER2009_Validation_Report_2.2.pdf
# And extract values. plot prevalence by category


# 5. Calculate distance to nearest health facility in km for each survey point
# Health facility shapefile available here
# https://data.humdata.org/dataset/burkina-faso-healthsites
# Hint1: the distm function generates a distance matrix between all pairs of points
# Hint2: the apply() function can find the minimum value across rows or columns







