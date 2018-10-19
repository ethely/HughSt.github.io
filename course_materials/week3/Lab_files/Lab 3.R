# Load relevant libraries
library(spatstat)
library(raster)
library(sp)
library(lgcp)
library(geoR)
library(gtools)
library(lme4)
library(leaflet)
library(oro.nifti)


# Open Namibia malaria case data
setwd("...")
CaseControl<-read.csv("CaseControl.csv")
Cases<-CaseControl[CaseControl$case==1,]
Controls<-CaseControl[CaseControl$case==0,]

# And boundary file
NAM_Adm0<-raster::getData('GADM',country='NAM',level=0)

# Convert to a SPDF
CaseControl_SPDF <- SpatialPointsDataFrame(coords = CaseControl[,c("long", "lat")],
                                           data = CaseControl[,c("household_id", "case")])

# Let's plot and see what we have
case_color_scheme <- colorNumeric(c("blue", "red"), CaseControl_SPDF$case)
leaflet() %>% addTiles() %>% addCircleMarkers(data=CaseControl_SPDF, 
                                              color = case_color_scheme(CaseControl_SPDF$case),
                                              radius = 2)


# To generate a kernel density estimate, we first 
# need to generate point pattern object of points
# First need to define a window
Nam_Owin <- owin(xrange=range(CaseControl$long),yrange=range(CaseControl$lat))

# Now define a ppp of the case data
Cases_ppp <- ppp(Cases$long, Cases$lat, window = Nam_Owin)
plot(Cases_ppp)

# Plot kernel density estimate
plot(density(Cases_ppp)) # Units are intensity of points per unit square

# Try with different bandwidths
plot(density(Cases_ppp,0.02))
plot(density(Cases_ppp,0.1))
plot(density(Cases_ppp,bw.ppl)) # automatic bandwidth selection based on cross-validation

# Map on leaflet (needs to be a rasterLayer object)
density_raster <- raster(density(Cases_ppp, bw.ppl), crs = crs(NAM_Adm0))
leaflet() %>% addTiles() %>% addRasterImage(density_raster, opacity=0.6)

# But this is just a density of cases, e.g. 
# it doesn't account for the denominator
# To do this, we can use the kelsall & diggle method, 
# which calculates the ratio of the
# density estimate of cases:controls

# First we have to add 'marks' to the points
# marks are just values associated with each point
# such as case or control (1/0), or magnitude of 
# point pattern of earthquakes. 
CaseControl_ppp <- ppp(CaseControl$long, CaseControl$lat, 
                       window = Nam_Owin, 
                       marks=as.factor(CaseControl$case))


risk_est <-  relrisk(CaseControl_ppp)
plot(risk_est)

# to plot on a web map, first specify the projection
risk_raster <- raster(risk_est, crs = crs(NAM_Adm0))

# Then plot using the leaflet package
pal = colorNumeric(palette=tim.colors(64), domain=exp(values(risk_raster)), na.color = NA)
leaflet() %>% addTiles("http://{s}.basemaps.cartocdn.com/light_all/{z}/{x}/{y}.png") %>% 
  addRasterImage(exp(risk_raster), opacity=0.6, col = pal)


# Interpolation of point (prevalence etc.) data
# Open BF malaria data
BF_malaria_data <- read.csv("BF_malaria_data.csv",
                            header=T)

BF_Adm_1 <- raster::getData("GADM", country="BFA", level=1)

# Calc prevalence
BF_malaria_data$prevalence <- BF_malaria_data$positives / BF_malaria_data$examined

# Inverse distance weighting
BF_malaria_window<-owin(xrange=range(BF_malaria_data$longitude),yrange=range(BF_malaria_data$latitude))
BF_malaria_data_ppp<-ppp(BF_malaria_data$longitude,BF_malaria_data$latitude,
                         marks=BF_malaria_data$prevalence,window=BF_malaria_window)

par(mfrow=c(2,2))
plot(idw(BF_malaria_data_ppp, power=0.2, at="pixels"),col=heat.colors(20), main="power = 0.2")
plot(idw(BF_malaria_data_ppp, power=0.5, at="pixels"),col=heat.colors(20), main="power = 0.5")
plot(idw(BF_malaria_data_ppp, power=1, at="pixels"),col=heat.colors(20), main="power = 0.1")
plot(idw(BF_malaria_data_ppp, power=2, at="pixels"),col=heat.colors(20), main="power = 2") # Larger power puts more weight on nearer values

# Plot using leaflet
BF_malaria_prev_idw_raster <- raster(idw(BF_malaria_data_ppp, power=0.2, at="pixels"),
                                     crs= crs(BF_Adm_1))

colPal <- colorNumeric(tim.colors(), BF_malaria_prev_idw_raster[], na.color = NA)
leaflet() %>% addTiles() %>% addRasterImage(BF_malaria_prev_idw_raster, col = colPal, opacity=0.7) %>%
  addLegend(pal = colPal, values = BF_malaria_prev_idw_raster[])


# To calculate the 'best' power to use, you can use cross-validation. 
CV_idw_1<-idw(BF_malaria_data_ppp, power=1, at="points")
plot(BF_malaria_data_ppp$marks, CV_idw_1, asp=1) # Not very good!

# Calc MSE
library(Metrics)
mse(BF_malaria_data_ppp$marks,CV_idw_1) # Mean squared error


# Kriging
# Before Kriging, good to transform prevalence 
# data to something vaguely normal
# Here we use the logistic transformation (log odds)
BF_malaria_data$log_odds <- logit(BF_malaria_data$prevalence)
hist(BF_malaria_data$log_odds)


# First have to create a geodata object with the package GeoR
# this wants dataframe of x,y and data 
BF_malaria_data_geo<-as.geodata(BF_malaria_data[,c("longitude","latitude","log_odds")])

# We can plot a summary plot
plot(BF_malaria_data_geo, lowes=T) # the lowes option gives us lowes curves for relationship with x and y
# potentially a trend on x and y

plot(BF_malaria_data_geo, lowes=T,trend="2nd") # Trend option regresses on x and y

# Now generate and plot a variogram
MaxDist=max(dist(BF_malaria_data[,c("longitude","latitude")]))  /2 # the max distance you should estimate is half max interpoint distance 
VarioCloud<-variog(BF_malaria_data_geo,option="cloud",max.dist=3)
plot(VarioCloud) # all pairwise comparisons

# To produce binned variogram
Vario<-variog(BF_malaria_data_geo, max.dist = MaxDist)
plot(Vario)
Vario<-variog(BF_malaria_data_geo,max.dist=MaxDist,uvec=seq(0.01,MaxDist,0.12))
Vario$n # Shows you the number in each bin
min(Vario$n)# should be at least 30 pairs in each bin
plot(Vario,pch=16)

# Fit variogram model by minimized least sqaures
VarioMod_lin<-variofit(VarioCloud, cov.model = "linear")
VarioMod_sph<-variofit(VarioCloud, cov.model = "sph")
VarioMod_exp<-variofit(VarioCloud, cov.model = "exp")

# plot results
lines(VarioMod_lin,col="green",lwd=2)
lines(VarioMod_sph,col="blue",lwd=2)
lines(VarioMod_exp,col="red",lwd=2) # In this example, all models converge on essentially the same line

VarioMod_lin
VarioMod_sph
VarioMod_exp
# lines model has lower sum of squares so 'better'

# Use variogram to Krig values at prediction locations 
# First get grid of points from the IDW example for comparison
# could use the expand.grid function 
IDW<-idw(BF_malaria_data_ppp, power=0.2, at="pixels")
pred_grid_x<-rep(IDW$xcol,length(IDW$yrow))
pred_grid_y<-sort(rep(IDW$yrow,length(IDW$xcol)))
pred_grid<-cbind(pred_grid_x,pred_grid_y)


KrigPred <- krige.conv(BF_malaria_data_geo, loc=pred_grid,
                       krige=krige.control(obj.model=VarioMod_lin))

# Visualize predictions
image(KrigPred,col=heat.colors(50))

# Back transform to prevalence
KrigPred_prev<-inv.logit(KrigPred$predict)
KrigPred_raster <- rasterFromXYZ(data.frame(x=pred_grid_x,
                                 y=pred_grid_y,
                                 z=KrigPred_prev))

plot(KrigPred_raster)
points(BF_malaria_data[,c("longitude","latitude")],
       cex = BF_malaria_data$prevalence)


# Its straightforward to CV kriged predictions in geoR
xvalid_result <- xvalid(BF_malaria_data_geo, model = VarioMod_lin) # By default it xvalidates point by point
plot(xvalid_result$data,xvalid_result$predicted, asp=1) # log odds scale
abline(0,1)







