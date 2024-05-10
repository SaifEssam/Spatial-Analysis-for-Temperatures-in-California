read.csv("ENV.csv")
data =read.csv("ENV.csv")
#calculating AUTOCORRELATION USING MOARNS I 
library(ape)
data.dists<-as.matrix(dist(cbind(data$LAT,data$LONG)))
View(data.dists)
data.dist.inv<-1/data.dists
diag(data.dist.inv)<-0
Moran.I(data$AVG,data.dist.inv)
l.moransI =Moran.I(data$AVG,data.dist.inv)

# localized morans and morans plot 
library(lctools)
l.moran2<-l.moransI(data.dists,100,data$AVG, WType='Bi-square', scatter.plot = TRUE, family = "adaptive")

library(spdep)
library(spatialEco)
library(sp)

morans.plot(data$AVG, coords =coordinates(data))
crs    <- CRS("+init=epsg:28992") 

coords <- data[ , c("LONG", "LAT")]

spdf <- SpatialPointsDataFrame(coords      = coords,
                               data        = data, 
                               proj4string = crs)

knn<-knearneigh(spdf, k=100, longlat = NULL)
knn2nb<-knn2nb(knn)
mp <- moran.plot(data$AVG, nb2listw(knn2nb))



# lisa plot 
library(ncf)
lisa <- lisa(data$LONG, data$LAT, data$AVG, neigh = 100, latlon=TRUE)
lisa
plot(lisa)

local <- localmoran(x = data$AVG, listw = nb2listw(knn2nb))
moran.map <- cbind(data, local)



boxplot(data$AVG)

# diffrent plots(bubble plot)
library(spatial)
library(plot3D)
library(plot3Drgl)
symbols(data$LONG, data$LAT, circles=data$AVG)


radius <- sqrt(data$AVG/ pi )
symbols(data$LONG, data$LAT, circles=radius)
symbols(data$LONG, data$LAT, circles=radius, inches=0.15, fg="white", bg="red")

# diffrent plots(3D plot)

scatter3D(data$LONG,data$LAT,data$AVG, zcol=data$AVG)
scatter3D(data$LONG,data$LAT,data$AVG, zcol=data$AVG,ticktype="detailed")
scatter3D(data$LONG,data$LAT,data$AVG, zcol=data$AVG,pty="g",ticktype="detailed")

library("plot3Drgl")
plotrgl()


##fitting a trend surface model by least squares
library(spatial)
x<-data$LONG
y<-data$LAT
z<-data$AVG

fit.sfc3 <- surf.ls(3,x,y,z)
summary(fit.sfc3)
fit.sfc3$beta
trsurf3 <- trmat(fit.sfc3, min(x), max(x), min(y), max(y), 50)
scatter3D(x,y,z,surf = list(x = trsurf3$x, y = trsurf3$y, z = trsurf3$z,
                            NAcol = "grey", shade = 0.1))
plotrgl()


fit.sfc2 <- surf.ls(2,x,y,z)
summary(fit.sfc2)
fit.sfc2$beta
trsurf2 <- trmat(fit.sfc2, min(x), max(x), min(y), max(y), 50)
scatter3D(x,y,z,surf = list(x = trsurf2$x, y = trsurf2$y, z = trsurf2$z,
                            NAcol = "grey", shade = 0.1))
plotrgl()


fit.sfc1 <- surf.ls(1,x,y,z)
summary(fit.sfc1)
fit.sfc2$beta
trsurf1 <- trmat(fit.sfc1, min(x), max(x), min(y), max(y), 50)
scatter3D(x,y,z,surf = list(x = trsurf1$x, y = trsurf1$y, z = trsurf1$z,
                            NAcol = "grey", shade = 0.1))
plotrgl()

fit.sfc4 <- surf.ls(4,x,y,z)
summary(fit.sfc4)
trsurf4 <- trmat(fit.sfc4, min(x), max(x), min(y), max(y), 50)
scatter3D(x,y,z,surf = list(x = trsurf4$x, y = trsurf4$y, z = trsurf4$z,
                            NAcol = "grey", shade = 0.1))
plotrgl()


library(plot3Drgl)
contour(trsurf4)

#fitting the variogram 


library(gstat)
#exponential variogram
evgm <- variogram( data[, c("AVG")]~1,spdf)
plot(evgm)

fvgm3 <- fit.variogram(evgm,vgm("Exp"))
fvgm3

summary(fvgm3)
plot(evgm,model=fvgm3)
attr(fvgm3, "SSErr")

#spherical variogram 
fvgm2 <- fit.variogram(evgm,vgm("Sph"))
fvgm2
plot(evgm,model=fvgm2)
attr(fvgm2, "SSErr")

#gaussian variogram 
fvgm1 <- fit.variogram(evgm,vgm("Gau"))
fvgm1
plot(evgm,model=fvgm1)
attr(fvgm1, "SSErr")


#fitting the grid 
s.grid2 <- spsample(spdf, type = "regular", n = 6000)
krig.est <- krige(data[, c("AVG")] ~ 1, spdf, newdata = s.grid2, model = fvgm3)
spplot(krig.est["var1.pred"])
spplot(krig.est["var1.var"])

krig.grid <- SpatialPixelsDataFrame(krig.est, krig.est@data)
levs <- c(10, 50, 60, 70, 80, Inf)
var.levs <- c(30, 70, 90, 100, 120, Inf)
#using t map but proplem in variance map
library(tmap)
krig.map.est <- tm_shape(krig.grid) +
  tm_raster(col = 'var1.pred', breaks = levs, title = 'ozone', palette = 'Reds') +
  tm_layout(legend.bg.color = 'white', legend.frame = TRUE)

krig.map.var <- tm_shape(krig.grid) +
  tm_raster(col = 'var1.var', breaks = var.levs, title = 'Estimate Variance', palette = 'Reds') +
  tm_layout(legend.bg.color = 'white', legend.frame = TRUE)

tmap_arrange(krig.map.est, krig.map.var)





library(phylin)
Long <- seq(from=-130,to=-112,by=0.03)
Lat <- seq(from=30,to=43,by=0.03)
grid<-expand.grid(Long,Lat)
gridsp<- SpatialPoints(grid,proj4string = crs)


krig.est <- krige(data[, c("AVG")] ~ 1, spdf, newdata = gridsp, model = fvgm3)
kriging<-as.matrix(krig.est$var1.pred)
grid.image(kriging, grid, main='ordinary kriging', xlab='Longitude', 
           ylab='Latitude', sclab="Genetic distance to sample s2", colFUN = terrain.colors)


variance<-as.matrix(krig.est$var1.var)
grid.image(variance, grid, main='ordinary kriging', xlab='Longitude', 
           ylab='Latitude', sclab="Genetic distance to sample s2",  colFUN = terrain.colors)


#geographical weighted regression model
library(spgwr)
library(ggplot2)

model1 <- lm(data$AVG ~ data$LONG + data$LAT)
summary(model1)
plot(model1)

par(mfrow=c(2,2))
plot(model1)

resids<-residuals(model1)
colours <- c("dark blue", "blue", "red", "dark red")
map.resids <- SpatialPointsDataFrame(data=data.frame(resids), coords=cbind(x,y)) 
spplot(map.resids, cuts=quantile(resids), col.regions=colours, cex=1) 


#calculate kernel bandwidth
GWRbandwidth <- gwr.sel(data$AVG ~ data$LAT+data$LONG, data=data, coords=cbind(x,y),adapt=T) 

gwr.model = gwr(data$AVG ~ data$LAT+data$LONG, data=data, coords=cbind(x,y), adapt=GWRbandwidth, hatmatrix=TRUE, se.fit=TRUE) 

gwr.model

results<-as.data.frame(gwr.model$SDF)
head(results)

data$coefLONG<-results$LONG
data$coefLAT<-results$LAT
gwr.point1<-ggplot(data, aes(x=x,y=y))+geom_point(aes(colour=data$coefLONG))+scale_colour_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0, space = "rgb", na.value = "grey50", guide = "colourbar", guide_legend(title="Coefs"))
gwr.point3+geom_path(data,aes(Long, Lat, group=id), colour="grey")+coord_equal()


gwr.point3<-ggplot(data, aes(x=x,y=y))+geom_point(aes(colour=data$coefLAT))+scale_colour_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0, space = "rgb", na.value = "grey50", guide = "colourbar", guide_legend(title="Coefs"))





