library(rgdal)
library(gdalUtils)
library(doParallel)
library(sf)
library(rgdal)

#'data' is a dataframe of latitude and longitude coordinates.

properties<-c("cec","clay","silt","sand","soc","nitrogen","phh2o")

cl <- makeCluster(detectCores()-1,outfile=paste("outfile",Sys.Date(),".txt",sep=""))
registerDoParallel(cl)
start<-Sys.time()

all_soil_propoerties<-foreach(i=properties,.combine="cbind",.packages=c("rgdal","gdalUtils","sf"),.export="data") %dopar% {

  print(paste("Starting",i,Sys.time()))
  voi = i # variable of interest
depth = "5-15cm"
layer = "mean"
voi_layer = paste(voi,depth,layer, sep="_") # layer of interest 
webdav_path="/vsicurl/https://files.isric.org/soilgrids/latest/data/"

spdata=st_as_sf(data,coords = c("Longitude", "Latitude"), crs = 4326)

igh='+proj=igh +lat_0=0 +lon_0=0 +datum=WGS84 +units=m +no_defs'
spdata_igh=st_transform(spdata, igh)

data_igh=data.frame(st_coordinates(spdata_igh),id=spdata_igh$Env)

fun_pixel_values=function(rowPX,data,VOI,VOI_LYR){
  as.numeric(
    gdallocationinfo(
      srcfile=paste0(webdav_path,"/",VOI,"/", VOI_LYR,".vrt"),
      x=data[rowPX,"X"],
      y=data[rowPX,"Y"],
      geoloc=TRUE,
      valonly=TRUE))
}

print("Unlisting values....")
unlist(lapply(1:nrow(data_igh),function(x){fun_pixel_values(x,data_igh,voi,voi_layer)}))
}
stopImplicitCluster()
stopCluster(cl)

colnames(all_soil_propoerties)<-properties
rownames(all_soil_propoerties)<-data$Env
all_soil_propoerties[all_soil_propoerties==-32768]<-NA


print("Finished!")
print(Sys.time()-start)


