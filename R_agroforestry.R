
#= = = = = = = = = = = = = = = =
#- - - - - - - - - - - - - - - -
#       PROJET MITIGATE +
#- - - - - - - - - - - - - - - -
#= = = = = = = = = = = = = = = =

#Code created by Paul Bostyn the 26/10/2023 under the behalf of the Mitigate+
#Project with the CGIAR 

#-_-_-_-_-_-_-_-_-_-_ NDC AGROFORESTRY -_-_-_-_-_-_-_-__-_--_-__--_-_

## 1 - Packages to install  ----
library(tidyverse)
library(readxl)
library(raster)
library(leaflet)
library(stars)
library(sf)
library(ggplot2)
library(dplyr)
library(rlist)
library(terra)
library(igraph)
library(writexl)
library(sp)
library(ncdf4)
library(quantreg)
library(caret)
library(gdata)
#S'il y a un problème avec les packages rgdal et rgeos, il faut utilisr les package terra
# et sf. Par exemple pour importer un shapefil, il est possible d'utiliser :
#st_read("Data/Agroforestry2/PAL-polygon.shp")


#For this code the use of R_code_1 is necessary to have the df2All data base.

#2 - Import and clean data ----

### a) For the NasaDem Elevation model Data ----

"Define the emplacement of the folder where there are the raster"
emplacementfolder2 <- "Data/Agroforestry"
"Define the source of the folder"
origin2 <- "Data/Agroforestry/"

"Define the territory to which the rasters will be adjusted"
Rastdecoupe <- read_sf("Data/Data_Kenyaboundary/kenboundary.shp")
crsref <-crs(Rastdecoupe)

"Stock the rasters in a list"
files <- list.files(emplacementfolder2)#Give the list of all name of the rasters files in the selected folder
list_raster <- list() #Create the list to stock the rasters
for (i in 1:length(files)){
  name_file <- files[[i]]#Assigne the name of the ième element of the list "files"
  real_name <- paste(origin2,name_file, sep="")#Attribute a name of each element
  list_raster <- append(list_raster,raster(real_name))#Add each ieme element in list_raster
  names(list_raster)[i] <- paste("R", i, sep = "")#Rename each element added to list_raster
}
remove(files, i, name_file, origin2, real_name,emplacementfolder2)

"Assign the same crs to each raster"
for (i in 1:length(list_raster)) {
  crs(list_raster[[i]]) <- crs(Rastdecoupe)
}

res(list_raster[[2]])#Check the difference of resolution between raster
res(list_raster[[69]])

"Changing the resolution of raster files "
numdebutgrouprast <- 1
numfingrouprast <- 68
for(i in numdebutgrouprast:numfingrouprast){
  list_raster[[i]] <- aggregate(list_raster[[i]],fact = 299, fun = mean)
}
res(list_raster[[22]]) #Check the resolution that should be approximately the same to HYDE raster
for(i in numdebutgrouprast:numfingrouprast){
  list_raster[[i]] <- resample(list_raster[[i]], list_raster[[69]]) #list_rast[69] is the resolution of HYDE
}#The resolution of the HYDE raster is used as the reference.
res(list_raster[[22]])
res(list_raster[[69]])
"Warning resample function can alter the pixel layout it is necessary to check after her use."

list_dfofraster <- lapply(list_raster, function(x) {as.data.frame(rasterToPoints(x,fun=function(x){x>=0}))})
for (i in 1:length(list_dfofraster)) {
  names(list_dfofraster[[i]]) <- c("x", "y", "NasaDEM")#Rename the column
}
list_dfofraster <- list_dfofraster[- c(69)]#remove the HYDE raster, used as a reference for the resolution
NasaDEM <- list_dfofraster[[1]] #This raster is defined as the reference 
for (i in 2:length(list_dfofraster)) {
  NasaDEM <- bind_rows(list_dfofraster[[i]],NasaDEM) #Assembling raster
}
"Cleaning if necessary"
remove(list_raster,list_dfofraster)

"Checking the results in a GIS software of the NasaDEM transformed into vector data"
dfnasadem <- NasaDEM[,c("x","y","NasaDEM")]
coordinates(dfnasadem) <- ~ x + y
proj4string(dfnasadem) <- proj4string(Rastdecoupe)
dfnasadem  <- st_as_sf(dfnasadem)
#Extract shapefile to check if it corresponds to the initial raster layers
st_write(dfnasadem, "NasaDEMmapcheck.shp")
remove(dfnasadem)


### b) For the NasaDem Slope Data ----
"Define the emplacement of the folder where there are the raster"
emplacementfolder3 <- "Data/Agroforestry2/nasademslope"
"Define the source of the folder"
origin3 <- "Data/Agroforestry2/nasademslope/"
"Stock the rasters in a list"
files2 <- list.files(emplacementfolder3)#Give the list of all name of the rasters files in the selected folder
list_raster_slope <- list() #Create the list to stock the rasters
for (i in 1:length(files2)){
  name_file2 <- files2[[i]]#Assigne the name of the ième element of the list "files"
  real_name2 <- paste(origin3,name_file2, sep="")#Attribute a name of each element
  list_raster_slope <- append(list_raster_slope,raster(real_name2))#Add each ieme element in list_raster
  names(list_raster_slope)[i] <- paste("R", i, sep = "")#Rename each element added to list_raster
}
remove(files2, i, name_file2, origin3, real_name2,emplacementfolder3)
"Assign the same crs to each raster"
for (i in 1:length(list_raster_slope)) {
  crs(list_raster_slope[[i]]) <- crs(Rastdecoupe)
}
res(list_raster_slope[[2]])#Check the difference of resolution between raster
res(list_raster_slope[[69]])#Hyde raster as reference
"Changing the resolution of raster files "
numdebutgrouprast <- 1
numfingrouprast <- 68
for(i in numdebutgrouprast:numfingrouprast){
  list_raster_slope[[i]] <- aggregate(list_raster_slope[[i]],fact = 299, fun = mean)
}
res(list_raster_slope[[22]]) #Check the resolution that should be approximately the same to HYDE raster
for(i in numdebutgrouprast:numfingrouprast){
  list_raster_slope[[i]] <- resample(list_raster_slope[[i]], list_raster_slope[[69]]) #list_rast[69] is the resolution of HYDE
}#The resolution of the HYDE raster is used as the reference.
res(list_raster_slope[[22]])
res(list_raster_slope[[69]])
"Warning resample function can alter the pixel layout it is necessary to check after her use."
list_dfofraster2 <- lapply(list_raster_slope, function(x) {as.data.frame(rasterToPoints(x,fun=function(x){x>=0}))})
for (i in 1:length(list_dfofraster2)) {
  names(list_dfofraster2[[i]]) <- c("x", "y", "NasaSLOPE")#Rename the column
}
list_dfofraster2 <- list_dfofraster2[- c(69)]#remove the HYDE raster, used as a reference for the resolution
NasaSLOPE <- list_dfofraster2[[1]] #This raster is defined as the reference 
for (i in 2:length(list_dfofraster2)) {
  NasaSLOPE <- bind_rows(list_dfofraster2[[i]],NasaSLOPE) #Assembling raster
}
"Cleaning if necessary"
remove(list_raster_slope,list_dfofraster2)
"Checking the results in a GIS software of the NasaDEM transformed into vector data"
dfnasaslope <- NasaSLOPE[,c("x","y","NasaSLOPE")]
coordinates(dfnasaslope) <- ~ x + y
proj4string(dfnasaslope) <- proj4string(Rastdecoupe)
dfnasaslope  <- st_as_sf(dfnasaslope)
#Extract shapefile to check if it corresponds to the initial raster layers
st_write(dfnasaslope, "NasaSLOPEmapcheck.shp")
remove(dfnasaslope)

### c) For the World Clim - Minim Temperature data ----
"Define the emplacement of the folder where there are the rasters files"
emplacementfolder4 <- "Data/Agroforestry2/worldclim"
"Define the source of the folder"
origin4 <- "Data/Agroforestry2/worldclim/"
"Stock the rasters in a list"
files3 <- list.files(emplacementfolder4)#Give the list of all name of the rasters files in the selected folder
list_raster_temp <- list() #Create the list to stock the rasters
for (i in 1:length(files3)){
  name_file3 <- files3[[i]]#Assigne the name of the ième element of the list "files"
  real_name3 <- paste(origin4,name_file3, sep="")#Attribute a name of each element
  list_raster_temp <- append(list_raster_temp,raster(real_name3))#Add each ieme element in list_raster
  names(list_raster_temp)[i] <- paste("R", i, sep = "")#Rename each element added to list_raster
}
remove(files3, i, name_file3, origin4, real_name3,emplacementfolder4)
"Assign the same crs to each raster"
for (i in 1:length(list_raster_temp)) {
  crs(list_raster_temp[[i]]) <- crs(Rastdecoupe)
}
res(list_raster_temp[[2]])#Check the difference of resolution between raster - World Clim Raster vs HYDE
res(list_raster_temp[[13]])#HYDE raster
"Cut the raster with the Kenya boundaries"
list_raster_temp2 <- list()
for (i in 1:length(list_raster_temp)){
  list_raster_temp2 <- append(list_raster_temp2, crop(list_raster_temp[[i]],Rastdecoupe))
}#crop() allow to cut the raster according to the Kenya territory shp
list_raster_temp <- list_raster_temp2
remove(list_raster_temp2)
"Changing the resolution of raster files "
numdebutgrouprast <- 1
numfingrouprast <- 12
for(i in numdebutgrouprast:numfingrouprast){
  list_raster_temp[[i]] <- resample(list_raster_temp[[i]], list_raster_temp[[13]]) #list_rast[13] is the resolution of HYDE
}#The resolution of the HYDE raster is used as the reference.
remove(numdebutgrouprast,numfingrouprast)
res(list_raster_temp[[11]])
res(list_raster_temp[[13]])
"Warning resample function can alter the pixel layout it is necessary to check after her use."
list_dfofraster3 <- lapply(list_raster_temp, function(x) {as.data.frame(rasterToPoints(x,fun=function(x){x>=0}))})
"Group the data of World Clim"
Dfreference <- list_dfofraster3[[1]]
coordinates(Dfreference)=~x+y
proj4string(Dfreference) <- proj4string(Rastdecoupe)
Dfreference <- Dfreference[Rastdecoupe,]
dftempmin <- as.data.frame(Dfreference)
dftempmin$x <- round(x = dftempmin$x, 10)
dftempmin$y <- round(x = dftempmin$y, 10)
names(dftempmin)[3] <- "Rast1"
remove(Dfreference)
for (i in 2:length(list_dfofraster3)){
  Dfreference <- list_dfofraster3[[i]]
  coordinates(Dfreference)=~x+y
  proj4string(Dfreference) <- proj4string(Rastdecoupe)
  Dfreference <- Dfreference[Rastdecoupe,]
  Dfreference <- as.data.frame(Dfreference)
  nam <- paste("Rast", i , sep = "")#donne un nom Rast1 au premier elemnt de la listdf_rast
  names(Dfreference)[3] <- nam
  Dfreference$x <- round(x = Dfreference$x, 10)
  Dfreference$y <- round(x = Dfreference$y, 10)
  dftempmin <- full_join(dftempmin, Dfreference, by = c("x","y"))
  remove(Dfreference,nam)
}
"Rename the columns"
dftempmin <- dftempmin %>% rename(WCtmin012021 = Rast1,WCtmin022021 = Rast2,
                                  WCtmin032021 = Rast3,WCtmin042021 = Rast4,
                                  WCtmin052021 = Rast5,WCtmin062021 = Rast6,
                                  WCtmin072021 = Rast7,WCtmin082021 = Rast8,
                                  WCtmin092021 = Rast9,WCtmin102021 = Rast10,
                                  WCtmin112021 = Rast11,WCtmin122021 = Rast12)#01 = January...
dftempmin <- dftempmin[ , !(names(dftempmin) %in% c("Rast13"))]#remove the col of HYDE
"Cross the data"
df2All <- full_join(df2All,dftempmin , by = c("x","y"), multiple = "all")
remove(dftempmin)

"Crossing the data of NasaDEM to have only one dataframe"
### d) For the NasaDem Elevation model Data (crossing data) ----
NasaDEM2 <- NasaDEM
NasaDEM2$x <- round(x = NasaDEM2$x, 10)
NasaDEM2$y <- round(x = NasaDEM2$y, 10)
df2All <- full_join(df2All, NasaDEM2, by = c("x","y"), multiple = "all")
df2All <- df2All[is.na(df2All$GFW2010) == "FALSE", ]#Permet de supprimer les pixels en plus
#ne correspondant pas au raster de base de GFW2010
#The for function below take 2 hours 
for(i in 1:length(df2All$x)){
  id1<-df2All$IDpixel[i]
  for(n in 1:length(df2All$IDpixel)){
    df2All$checkID[n]<-ifelse(df2All$IDpixel[n]==id1,1,0)
  }
  dfsub<-subset(df2All, checkID ==1)
  meansub<-mean(dfsub$NasaDEM)
  df2All$NasaDEM[i]<-meansub
}#allow to give the mean value of the duplicate before removing them
#Indeed, The NAsadem data can have two or more values for the same pixel of GFW, so
#we create this for function above to have the mean value for each pixel where
#there are many values.
remove(dfsub,meansub)
df2All = df2All[!duplicated(df2All$IDpixel),]#Allow to remove the duplicate value
"Checking the results in a GIS software of the NasaDEM assimilation to the df2All"
dfnasademcheck2 <- df2All[,c("x","y","NasaDEM")]
coordinates(dfnasademcheck2) <- ~ x + y
proj4string(dfnasademcheck2) <- proj4string(Rastdecoupe)
dfnasademcheck2  <- st_as_sf(dfnasademcheck2)
#Extract shapefile to check if it corresponds to the initial raster layers
st_write(dfnasademcheck2, "NasaDEMmapcheck2.shp")
remove(dfnasademcheck2)

### e) For the NasaDem Slope Data (crossing data) ----
NasaSLOPE2 <- NasaSLOPE
NasaSLOPE2$x <- round(x = NasaSLOPE2$x, 10)
NasaSLOPE2$y <- round(x = NasaSLOPE2$y, 10)
df2All <- full_join(df2All, NasaSLOPE2, by = c("x","y"), multiple = "all")
df2All <- df2All[is.na(df2All$GFW2010) == "FALSE", ]#Permet de supprimer les pixels en plus
#ne correspondant pas au raster de base de GFW2010
#The for function below take 2 hours 
for(i in 1:length(df2All$x)){
  id1<-df2All$IDpixel[i]
  for(n in 1:length(df2All$IDpixel)){
    df2All$checkID[n]<-ifelse(df2All$IDpixel[n]==id1,1,0)
  }
  dfsub<-subset(df2All, checkID ==1)
  meansub<-mean(dfsub$NasaSLOPE)
  df2All$NasaSLOPE[i]<-meansub
}#allow to give the mean value of the duplicate before removing them (because some data can have more pixel than the df
#so there are duplicate and we need to remove them(for example in 1 pixel it can have 2 informations because it's coming from 
#a raster layer that was having more pixels informations, so we calculate the mean value of them and after remove the duplicate.))
df2All = df2All[!duplicated(df2All$IDpixel),]#Allow to remove the duplicate value
"Checking the results in a GIS software of the NasaSLOPE assimilation to the df2All"
dfnasaslopecheck<- df2All[,c("x","y","NasaSLOPE")]
coordinates(dfnasaslopecheck) <- ~ x + y
proj4string(dfnasaslopecheck) <- proj4string(Rastdecoupe)
dfnasaslopecheck  <- st_as_sf(dfnasaslopecheck)
#Extract shapefile to check if it corresponds to the initial raster layers
st_write(dfnasaslopecheck, "NasaSLOPEcheck2.shp")
remove(dfnasaslopecheck)

### f) For the SPAM data physical area ----
"Define the emplacement of the folder where there are the raster"
emplacementfolder3 <- "Data/Agroforestry2/spam_physicalarea"
"Define the source of the folder"
origin5 <- "Data/Agroforestry2/spam_physicalarea/"
"Stock the rasters in a list"
files3 <- list.files(emplacementfolder3)#Give the list of all name of the rasters files in the selected folder
list_raster_spamphysic <- list() #Create the list to stock the rasters
for (i in 1:length(files3)){
  name_file3 <- files3[[i]]#Assigne the name of the ième element of the list "files"
  real_name3 <- paste(origin5,name_file3, sep="")#Attribute a name of each element
  list_raster_spamphysic <- append(list_raster_spamphysic,raster(real_name3))#Add each ieme element in list_raster
  names(list_raster_spamphysic)[i] <- paste("R", i, sep = "")#Rename each element added to list_raster
}
remove(files3, i, name_file3, origin5, real_name3,emplacementfolder3)
"Assign the same crs to each raster"
for (i in 1:length(list_raster_spamphysic)) {
  crs(list_raster_spamphysic[[i]]) <- crs(Rastdecoupe)
}
res(list_raster_spamphysic[[2]])#Check the difference of resolution between raster
res(list_raster_spamphysic[[43]])#Here it's the raster of grazing from HYDE

"Changing the resolution of raster files "
#There is 42 raster of SPAM
for(i in 1:42){
  list_raster_spamphysic[[i]] <- resample(list_raster_spamphysic[[i]], list_raster_spamphysic[[43]]) #list_rast[43] is the resolution of HYDE
}#The resolution of the HYDE raster is used as the reference.
res(list_raster_spamphysic[[2]])
res(list_raster_spamphysic[[43]])
"Warning resample function can alter the pixel layout it is necessary to check after her use."
list_raster_spamphysic <- list_raster_spamphysic[- c(43)]#remove the HYDE raster, used as a reference for the resolution
#Transform into df
list_dfofraster5 <- lapply(list_raster_spamphysic, function(x) {as.data.frame(rasterToPoints(x,fun=function(x){x>=0}))})


#Join all the df of SPAM coming from the rasters
#this data set below is used as a first reference to then join the other to it.
dfspamref <- list_dfofraster5[[1]]
coordinates(dfspamref)=~x+y
proj4string(dfspamref) <- proj4string(Rastdecoupe)
dfspamref <- dfspamref[Rastdecoupe,]
dfspamref1 <- as.data.frame(dfspamref)
dfspamref1$x <- round(x = dfspamref1$x, 10)
dfspamref1$y <- round(x = dfspamref1$y, 10)
names(dfspamref1)[3] <- "Rast1"
remove(dfspamref)
for (i in 2:length(list_dfofraster5)){
  dfspamref <- list_dfofraster5[[i]]
  coordinates(dfspamref)=~x+y
  proj4string(dfspamref) <- proj4string(Rastdecoupe)
  dfspamref <- dfspamref[Rastdecoupe,]
  dfspamref <- as.data.frame(dfspamref)
  nam <- paste("Rast", i , sep = "")#donne un nom Rast2 au premier elemnt de la listdf_rast
  names(dfspamref)[3] <- nam
  dfspamref$x <- round(x = dfspamref$x, 10)
  dfspamref$y <- round(x = dfspamref$y, 10)
  dfspamref1 <- full_join(dfspamref1, dfspamref, by = c("x","y"))
  remove(dfspamref,nam)
}
dfspamref1 <- dfspamref1 %>% rename(
  SPAMACOFPhysarea = Rast1, SPAMBANAPhysarea = Rast2, 
  SPAMBARLPhysarea = Rast3, SPAMBEANPhysarea = Rast4, SPAMCASSPhysarea = Rast5,
  SPAMCHICPhysarea = Rast6, SPAMCNUTPhysarea = Rast7, SPAMCOCOPhysarea = Rast8,
  SPAMCOTTPhysarea = Rast9, SPAMCOWPPhysarea = Rast10, SPAMGROUPhysarea = Rast11,
  SPAMLENTPhysarea = Rast12, SPAMMAIZPhysarea = Rast13, SPAMOCERPhysarea = Rast14,
  SPAMOFIBPhysarea = Rast15, SPAMOILPPhysarea = Rast16, SPAMOOILPhysarea = Rast17,
  SPAMOPULPhysarea = Rast18, SPAMORTSPhysarea = Rast19, SPAMPIGEPhysarea = Rast20,
  SPAMPLNTPhysarea = Rast21, SPAMPMILPhysarea = Rast22, SPAMPOTAPhysarea = Rast23,
  SPAMRAPEPhysarea = Rast24, SPAMRCOFPhysarea = Rast25,
  SPAMRESTPhysarea = Rast26, SPAMRICEPhysarea = Rast27,
  SPAMSESAPhysarea = Rast28, SPAMSMILPhysarea = Rast29,
  SPAMSORGPhysarea = Rast30, SPAMSOYBPhysarea = Rast31,
  SPAMSUGBPhysarea = Rast32, SPAMSUGCPhysarea = Rast33,
  SPAMSUNFPhysarea = Rast34, SPAMSWPOPhysarea = Rast35,
  SPAMTEASPhysarea = Rast36, SPAMTEMFPhysarea = Rast37,
  SPAMTOBAPhysarea = Rast38, SPAMTROFPhysarea = Rast39, SPAMVEGEPhysarea = Rast40,
  SPAMWHEAPhysarea = Rast41, SPAMYAMSPhysarea = Rast42)
df2All <- full_join(df2All, dfspamref1, by = c("x","y"))
remove(dfspamref1,list_dfofraster5,list_raster_spamphysic)



for (i in 1:length(list_dfofraster5)) {
  names(list_dfofraster2[[i]]) <- c("x", "y", "SPAMphys")#Rename the column
}
NasaSLOPE <- list_dfofraster2[[1]] #This raster is defined as the reference 
for (i in 2:length(list_dfofraster2)) {
  NasaSLOPE <- bind_rows(list_dfofraster2[[i]],NasaSLOPE) #Assembling raster
}
"Cleaning if necessary"
remove(list_raster_slope,list_dfofraster2)
"Checking the results in a GIS software of the NasaDEM transformed into vector data"
dfnasaslope <- NasaSLOPE[,c("x","y","NasaSLOPE")]
coordinates(dfnasaslope) <- ~ x + y
proj4string(dfnasaslope) <- proj4string(Rastdecoupe)
dfnasaslope  <- st_as_sf(dfnasaslope)
#Extract shapefile to check if it corresponds to the initial raster layers
st_write(dfnasaslope, "NasaSLOPEmapcheck.shp")
remove(dfnasaslope)

### g) For the soil classification data ----
"Import and cut the shapefile of soil classification"
shape2 <- st_read("Data/Agroforestry2/PAL-polygon.shp")%>%
  st_transform(st_crs(Rastdecoupe))#Same crs is applied
Shpcutting  <- st_as_sf(Rastdecoupe)#Transforming into sf
st_crs(shape2)==st_crs(Rastdecoupe)#Checking the same crs 
shape2cropped <- shape2[Shpcutting,]#Cutting the shapefile of soil
remove(shape2,Shpcutting)
"Remove the useless column"
shape2cropped<-shape2cropped[,-(1:8)]
shape2cropped<-shape2cropped[,-(2:3)]
soil_categories<-unique(shape2cropped$DESCRIP)#get the soil categories object
soil_categories[2]<-"All"


#3 - Gain of Agroforestry calcul and location of planting ----

#To launch the loop to get the gain from agroforestry (in $) and the
#location where the agroforestry specie can grow, launch the code from the 
# "Launching starting point" to the "Launching ending point"

#L'utilisateur doit lancer un algorithme avec le R_code1 et ensuite lancer 
#l'algorithme ci-dessous pour l'agroforesterie

#L'utilisateur doit remplir l'excel DF1plant,les règles à suivre sont les 
#suivantes : 

#-- For the maximum degree of slope indicated by levelslopemax  : 
#The maximum degree of slope needed (if there are no informations
#on the slope needed for the plant and according to The Barcelona Field Studies Centre,
#above 8.5 degrees of slope it's hard slope and so it's more complicated to grow smthing)
#So --> 8,5 need to be indicated in the first excel.

#-- For the soil type indicated by the column name "soiltype" : 
#If only one type of soil is needed for the plant or if all the type of soil are usable
#for this plant  we should put 0 on the 'Conditionformanysoiltype and 1 if there are many
#type of soil usable for this plant but not all type of soil.
#IF you don't know and to have results we use all the type of soil so we should add 'All' 
#to the value of the column 'soiltype'. 
#Then for the column soiltype2 to soiltype6, it's just in the case that there is different
#soil type needed. So if there is only one soil type needed or if we don't know or we
#considere that all the soil type are ok for the plant we should put "No" in value for
#these columns of from soiltype2 to soiltype6. 

"-- !! Be careful to not give the same soil type name to different soiltype column !!
(for example, soiltype 2 should not have the same name as soiltype3 or 4 or 5 or 6.)
Also, the columns soiltype2, soiltype3,soiltype4,soiltype5,soiltype6, can't take the value
'All' because there aim is only to have name of soil needed by the plant or 'No' 
if they are useless."


#-- In the column "Price" each species concerning by the comumn Speciescrop should
#have the value in $ of one tone. The source of this price should be inclue in the column
# 'Source2' 

#-- In the column "Increaseyields" each species should have a RR factor at which the use of
#Agroforestery increase yields and in the column "Source1" it should have the source of
#this increase in yield. For example, for Sesbania Sesban
#The RR factor is = (yield in agroforestry)/(yield in monoculture) = 1,6604 on maiz crop
#meaning that in average, an area of maiz in agroforestry could produce 1,6604 more 
#than in monoculture with the specie of Sesbania sesban used.

#-- In column "Speciescrop" it should have the specie or species in which, the use of the 
#specific agroforestry specie increase the yield. For example, for Sesbania Sesban
#the specific agroforestry specie is Sesbania Sesban and the "Speciescrop" = Maiz

#-- It's not possible to have two examples of the same specific agroforestry specie with
#the same Speciescrop, for examples, not possible to have 2 lines in the DF1plant excel
#with Sesbania Sesban and Maiz

#-- Don't put some space area in the name of the plant, please
#fullfill this space area with "_" 

#-- The name of the crop species should be included in the SPAM data, so, please
#give the same name of crop species as in SPAM, data for example it's in english in SPAM
#so put MAIZ in the "Speciescrop" column.  
#In the notepad "Species_name_for_R_agroforestry" there is the name of the species
#according to the SPAM data, please use it.

#---- Other specifications  : 


###  Launching starting point -- -- -- -- ---- 

"Insert which alogirhtm was used for the reforestatin (1 for algo minimizing loss,
2 for algo 2 ecological continuum part a, 3 for algo 2 ecological continuum part b) :"
Algouseagrofo <- 1
"If you want the production of the map containing all the pixel where the 
specific agroforestry specie can grow : 1 for yes, 2 for no 
(the name of the file extract should be change before doing it again)"
Havethemap <- 2

newvalagrofo <- 0#let it like that
AddedvalSPAMprod<-0#let it like that

#Please, you have to launch an alorithm with the R_code1 with a specific tree cover
#to get data on the reforestation area for each pixels
#and be sure that diffalgo1 or diffalgo4 or diffalgo42 are still in df2All.


#/!\
#IT IS ALSO VERY IMPORTANT TO LAUNCH THE CLEANING CODE after the use of the
#for loop below
#/!\




"Import the data used after to store data from each species and DF1plant need to be
completed."
DF1plant <- read_excel("Data/Data_Agrofo/DF1plant.xlsx")#Ici cet excel
#contient l'ensemble des espèces pour lesquels on souhaite connaître leur aire de répartition au Kenya
DF2plant <- read_excel("Data/Data_Agrofo/DF2plant.xlsx")#Ici ce df
#permet de stocker les données à l'intérieur pour chaque boucle donc pour chaque espèce.

Checking_namecol <- read_excel("Data/Data_Agrofo/Checking_namecol.xlsx")#Ici cet excel
#POur chaque espèce de DF1plant, le nom de la colonne SPAM sélectionnée affiliée à cette espèce
#est intégré dans ce dataset et permet en regardant DF1plant et Checking_namecol une vérification.

#Import dataset to store the benefits of each species :
#(Dans ce dataset, chaque ligne correspond à la ligne de DF1plant donc par exemple
#si la première ligne dans DF1plant est l'espèce Sesbania sesban, alors le benefice
#total généré par la mise en place de système agroforestier partout où l'espèce
#Sesbania sesban peut pousser sera la première ligne de Dfbenefcroptotal)
Dfbenefcroptotal <- read_excel("Data/Data_Agrofo/Results_gain_agrofo_Total.xlsx")#Ici cet excel
"Les trois premières lignes de Dfbenefcroptotal sont inutiles d'un point de vue résultats
(utile uniquement pour le traitement des données) et pourront être supprimé après l'utilisation
de la boucle"
#(Il en est de même pour le DfbenefNDC mais pour les benefices générés par
#l'atteinte des objectifs de NDC) :
DfbenefcropNDC <- read_excel("Data/Data_Agrofo/Results_gain_agrofo_NDC.xlsx")#Ici cet excel
"Il en est de même pour les 3 premières lignes de DfbenefcropNDC"

#Si besoin, pour vérifier la sélection des pixels aléatoires et leurs données:
Checkrandompixselect1 <- data.frame(IDpix = 1)
Checkrandompixselect2 <- data.frame(SPAMPhysarea = 1)
Checkrandompixselect3 <- data.frame(SPAMProd = 1)

"Le code ci-dessous fonctionne seulement pour les systèmes agroforestiers avec un RRlevel1>1"

#Remove the columns specific of each algorithm : 
if(Algouseagrofo == 1){
  df2All <- df2All[ , !(names(df2All) %in% c("loosevalalgo1","loosequantialgo1",
                                             "GFW2010_new_1","km2plantalgo1",
                                             "diffalgo4","diffalgo42"))]
}

if(Algouseagrofo == 2){
  df2All <- df2All[ , !(names(df2All) %in% c("kmplantalgo4","GFW2010_new_41","diffreforesta41",
                                             "diffalgo1","diffalgo42"))]
}

if(Algouseagrofo == 3){
  df2All <- df2All[ , !(names(df2All) %in% c("GFW2010_new_42","kmplantalgo42",
                                             "Plantationcondalgo4","loosequantialgo41",
                                             "diffalgo41","km2plantalgo41",
                                             "Plantcondialgo4","Plantalgo42",
                                             "celladjalgo4","Plantindicealgo42",
                                             "diffalgo1","diffalgo4"))]
}


for(i in 1:length(DF1plant$nameofplant)){
  Nameplantcrop <- DF1plant$Speciescrop[i]
  Nameofplant <- DF1plant$nameofplant[i]
  DF2plant$nameofplant[1]<-DF1plant$nameofplant[i]
  DF2plant$tempmin[1]<-DF1plant$tempmin[i]#The minimum (monthly) temperature needed (based on 2021)
  DF2plant$tempmax[1]<-DF1plant$tempmax[i]#The maximum (monthly) temperature needed (based on 2021)
  DF2plant$precipmin[1]<-DF1plant$precipmin[i]#The minimum annual average precipitation level needed (based on 2021)
  DF2plant$precipmax[1]<-DF1plant$precipmax[i]#The maximum annual average precipitation level needed (based on 2021)
  DF2plant$altitudmin[1]<-DF1plant$altitudmin[i]#The minimum altitude level needed
  DF2plant$altitudmax[1]<-DF1plant$altitudmax[i]#The maximum altitude level needed
  DF2plant$levelslopemax[1]<-DF1plant$levelslopemax[i]#The maximum degree of slope needed (if there are no informations
  #on the slope needed for the plant and according to The Barcelona Field Studies Centre,
  #above 8.5 degrees of slope it's hard slope and so it's more complicated to grow smthing)
  #So --> 8,5 need to be indicated in the first excel. 
  DF2plant$soiltype[1]<-DF1plant$soiltype[i]
  DF2plant$soiltype2[1]<-DF1plant$soiltype2[i]
  DF2plant$soiltype3[1]<-DF1plant$soiltype3[i]
  DF2plant$soiltype4[1]<-DF1plant$soiltype4[i]
  DF2plant$soiltype5[1]<-DF1plant$soiltype5[i]
  DF2plant$soiltype6[1]<-DF1plant$soiltype6[i]
  
  print("--------------------------------")
  print("L'espèce étudiée ici est :")
  print(Nameofplant)
  print("--------------------------------")
  
  "Insert a level of Ratio at which the use of the specie can increased the agricultural
  production (RRlevel1 = production with agroforestry/production without agroforestry"
  RRlevel1<-DF1plant$Increaseyields[i]
  "Insert a level of price for the specie of crop for 1 ton"
  priceofcropspecie<-DF1plant$Price[i]
  #dataset to check the good selection of the SPAM used :
  
  Conditionformanysoiltype <- DF1plant$Conditionformanysoiltype[i] 
  
  #see the "soil_categories" object which include the soil categories present in Kenya
  #So there is 7 type of soil given by the Rastdecoupe in Kenya.
  
  #To have more information on the soil type needed for the plant see :
  #this code below allows to check if the information written in the column
  #of "soiltype" is good and can correspond to soil_categories.
  if(DF2plant$soiltype[1]==soil_categories[1] || 
     DF2plant$soiltype[1]==soil_categories[2] || 
     DF2plant$soiltype[1]==soil_categories[3] || 
     DF2plant$soiltype[1]==soil_categories[4] || 
     DF2plant$soiltype[1]==soil_categories[5] || 
     DF2plant$soiltype[1]==soil_categories[6] || 
     DF2plant$soiltype[1]==soil_categories[7]){
    print("Great")
  } 
  else { print("This soil is not present in Kenya or you didn't write it the 
               exact way as 'soil_categories' object.")
    stop()}
  if(DF2plant$soiltype2[1]== "All" || DF2plant$soiltype3[1]== "All" ||
     DF2plant$soiltype4[1]== "All" || DF2plant$soiltype5[1]== "All" ||
     DF2plant$soiltype6[1]== "All"){
    print("You can't give to soiltype2,3,4,5,6 the value All. soiltype2,3,4,5,6 are
        only there if there are differents soil types needed for the plant. (Be careful to
        not give the same soil type). ")
  } 
  else { print("Ok (Be careful that you didn't give two times the same soil type name.)")
    }
  
  "Get the pixel answering the conditions"
  df2All$Tempmaxagrofocond<-df2All$Tempmaxcond#Create the columns to avoid problem
  df2All$Tempminagrofocond<-df2All$Tempmaxcond
  df2All$Precipagrofocond<-df2All$Tempmaxcond
  df2All$Altiagrofocond<-df2All$Tempmaxcond
  df2All$Slopeagrofocond<-df2All$Tempmaxcond
  df2All$Allcondagrofo<-df2All$Tempmaxcond
  for(n in 1:length(df2All$x)){
    df2All$Tempmaxagrofocond[n] <- ifelse(df2All$WCtmax012021[n]<DF2plant$tempmax[1], 
                                          ifelse(df2All$WCtmax022021[n]<DF2plant$tempmax[1],
                                                 ifelse(df2All$WCtmax032021[n]<DF2plant$tempmax[1],
                                                        ifelse(df2All$WCtmax042021[n]<DF2plant$tempmax[1],
                                                               ifelse(df2All$WCtmax052021[n]<DF2plant$tempmax[1],
                                                                      ifelse(df2All$WCtmax062021[n]<DF2plant$tempmax[1],
                                                                             ifelse(df2All$WCtmax072021[n]<DF2plant$tempmax[1],
                                                                                    ifelse(df2All$WCtmax082021[n]<DF2plant$tempmax[1],
                                                                                           ifelse(df2All$WCtmax092021[n]<DF2plant$tempmax[1],
                                                                                                  ifelse(df2All$WCtmax102021[n]<DF2plant$tempmax[1],
                                                                                                         ifelse(df2All$WCtmax112021[n]<DF2plant$tempmax[1],
                                                                                                                ifelse(df2All$WCtmax122021[n]<DF2plant$tempmax[1],1,0),0),0),0)
                                                                                           ,0),0),0),0),0),0),0),0)}
  
  for(w in 1:length(df2All$x)){
    df2All$Tempminagrofocond[w] <- ifelse(df2All$WCtmin012021[w]>=DF2plant$tempmin[1], 
                                          ifelse(df2All$WCtmin022021[w]>=DF2plant$tempmin[1],
                                                 ifelse(df2All$WCtmin032021[w]>=DF2plant$tempmin[1],
                                                        ifelse(df2All$WCtmin042021[w]>=DF2plant$tempmin[1],
                                                               ifelse(df2All$WCtmin052021[w]>=DF2plant$tempmin[1],
                                                                      ifelse(df2All$WCtmin062021[w]>=DF2plant$tempmin[1],
                                                                             ifelse(df2All$WCtmin072021[w]>=DF2plant$tempmin[1],
                                                                                    ifelse(df2All$WCtmin082021[w]>=DF2plant$tempmin[1],
                                                                                           ifelse(df2All$WCtmin092021[w]>=DF2plant$tempmin[1],
                                                                                                  ifelse(df2All$WCtmin102021[w]>=DF2plant$tempmin[1],
                                                                                                         ifelse(df2All$WCtmin112021[w]>=DF2plant$tempmin[1],
                                                                                                                ifelse(df2All$WCtmin122021[w]>=DF2plant$tempmin[1],1,0),0),0),0)
                                                                                           ,0),0),0),0),0),0),0),0)}
  
  for(d in 1:length(df2All$x)){
    df2All$Precipagrofocond[d] <- ifelse(df2All$WCprecimeanall[d]>=DF2plant$precipmin[1],
                                         ifelse(df2All$WCprecimeanall[d]<=DF2plant$precipmax[1],1,0),0)}
  
  for(n in 1:length(df2All$x)){
    df2All$Altiagrofocond[n] <- ifelse(df2All$NasaDEM[n]>=DF2plant$altitudmin[1],
                                       ifelse(df2All$NasaDEM[n]<=DF2plant$altitudmax[1],1,0),0)}
  
  for(z in 1:length(df2All$x)){
    df2All$Slopeagrofocond[z] <- ifelse(df2All$NasaSLOPE[z]<=DF2plant$levelslopemax[1],1,0)}
  
  for(w in 1:length(df2All$x)){
    df2All$Allcondagrofo[w] <- ifelse(df2All$Tempmaxagrofocond[w] == 1,
                                      ifelse(df2All$Tempminagrofocond[w] == 1,
                                             ifelse(df2All$Precipagrofocond[w] == 1,
                                                    ifelse(df2All$Altiagrofocond[w] == 1,
                                                           ifelse(df2All$Slopeagrofocond[w] == 1,
                                                                  1,0),0),0),0),0)}
  #Column Allcondagrofo allows to know where the specific species for agroforestry is
  #able to grow but in this column the soil type is not taking into account. 
  print("ok from here")
  
  "-- Some stats on the conditions :"
  
  "Number of pixels respecting the condition of temperature max"
  nrow(subset(df2All[c("Tempmaxagrofocond")], Tempmaxagrofocond == 1))
  "Number of pixels respecting the condition of temperature min"
  nrow(subset(df2All[c("Tempminagrofocond")], Tempminagrofocond == 1))
  "Number of pixels respecting the condition of Precipitation"
  nrow(subset(df2All[c("Precipagrofocond")], Precipagrofocond == 1))
  "Number of pixels respecting the condition of Altitude"
  nrow(subset(df2All[c("Altiagrofocond")], Altiagrofocond == 1))
  "Number of pixels respecting the condition of Slope"
  nrow(subset(df2All[c("Slopeagrofocond")], Slopeagrofocond == 1))
  
  
  "Get the df with the conditions and transform it into shp"
  dfcrossingmap <- subset(df2All[c("x","y","Allcondagrofo","IDpixel")])
  coordinates(dfcrossingmap) <- ~ x + y
  dfcrossingmap  <- st_as_sf(dfcrossingmap)
  crs(Rastdecoupe)
  st_crs(dfcrossingmap) <- st_crs(Rastdecoupe)
  
  #Extract le shapefile if necessary (line below)
  #st_write(dfcrossingmap, "mapcondtioncroisement.shp")
  
  "Join the shapefil to get for each type of soil if there are pixel answering the conditions"
  Shpsoiljoin <- st_join(dfcrossingmap, shape2cropped, left = FALSE, largest = TRUE)
  "Transform in df and make some arrangements"
  dfcheckingsoil <- as.data.frame(Shpsoiljoin)
  dfcheckingsoil <- dfcheckingsoil[ , !(names(dfcheckingsoil) %in% c("geometry"))]
  dfcheckingsoil$removagrofo <- dfcheckingsoil$Allcondagrofo
  for(z in 1:length(dfcheckingsoil$Allcondagrofo)){
    dfcheckingsoil$removagrofo[z] <- ifelse(dfcheckingsoil$DESCRIP[z] == "Waterbodies",1,0)
  }
  dfcheckingsoil <- subset(dfcheckingsoil,removagrofo == 0)
  dfcheckingsoil <- dfcheckingsoil[ , !(names(dfcheckingsoil) %in% c("removagrofo"))]
  remove(dfcrossingmap, Shpsoiljoin)
  #During the process of join and transform into df it looses 120 pixels where some waterbodies
  #which are not interesting for the analysis
  
  "To find out which pixels meet the conditions regarding the soil :"
  #Dans cette boucle ci-dessous c'est s'il y a un seul type de sol ou si
  #tous les types de sol sont possibles pour l'espèce rentrée dans DF1plant.
  dfcheckingsoil$Allcondagrofo2 <- dfcheckingsoil$Allcondagrofo
  if(Conditionformanysoiltype == 0){
    for(n in 1:length(dfcheckingsoil$Allcondagrofo)){
      dfcheckingsoil$Allcondagrofo2[n] <- ifelse(DF2plant$soiltype[1] == "All",
                                                     ifelse(dfcheckingsoil$Allcondagrofo[n] == 1,1,0)
                                                     ,ifelse(dfcheckingsoil$DESCRIP[n]==DF2plant$soiltype[1],
                                                             ifelse(dfcheckingsoil$Allcondagrofo[n] == 1,1,0),0))
    }
  }
  
  "If there are different types of soil"
  for(q in 1:length(dfcheckingsoil$Allcondagrofo)){
    if(Conditionformanysoiltype == 1){
      if(is.na(dfcheckingsoil$Allcondagrofo[q])=="FALSE"){
        if(dfcheckingsoil$Allcondagrofo[q] == 1){
          if(dfcheckingsoil$DESCRIP[q]==DF2plant$soiltype2[1] || 
           dfcheckingsoil$DESCRIP[q]==DF2plant$soiltype3[1] || 
           dfcheckingsoil$DESCRIP[q]==DF2plant$soiltype4[1] || 
           dfcheckingsoil$DESCRIP[q]==DF2plant$soiltype5[1] || 
           dfcheckingsoil$DESCRIP[q]==DF2plant$soiltype[1] || 
           dfcheckingsoil$DESCRIP[q]==DF2plant$soiltype6[1]){
            dfcheckingsoil$Allcondagrofo2[q]<-1
            }else{dfcheckingsoil$Allcondagrofo2[q]<-0}
         }else{dfcheckingsoil$Allcondagrofo2[q]<-0}
      }else{dfcheckingsoil$Allcondagrofo2[q]<-0}
    }
  }
  
  
  #The colum Allcondagrofo2 gives the pixels where the specie selected can grow.
  
  #Add the column of Allcondagrofo2 to df2All :
  dfchecksoljoin <- dfcheckingsoil[,c("IDpixel","Allcondagrofo2")]
  df2All <- full_join(df2All, dfchecksoljoin, by ="IDpixel")
  remove(dfchecksoljoin)
  
  # 4 - Results from the area where the plant can grow --
  "Number of pixels respecting all the conditions for the plant to grow"
  nrow(subset(dfcheckingsoil[c("Allcondagrofo2")], Allcondagrofo2 == 1))
  
  "Get the map of the pixel respecting all the conditions for the plant to grow"
  if(Havethemap==1){
    dfcheck1a <- df2All[,c("x","y","IDpixel")]
    dfchecksol <- dfcheckingsoil[,c("IDpixel","Allcondagrofo2")]
    dfcheck1b <- full_join(dfcheck1a, dfchecksol, by ="IDpixel")
    df2All <- full_join(df2All, dfchecksol, by ="IDpixel")
    coordinates(dfcheck1b) <- ~ x + y
    dfcheck1b  <- st_as_sf(dfcheck1b)
    st_crs(dfcheck1b) <- st_crs(Rastdecoupe)
    #Extract le shapefile
    st_write(dfcheck1b, "Soilchecking.shp")#Here it's the map of the pixel respecting the
    #requirements for the plant to grow
    remove(dfcheck1a,dfcheck1b,dfchecksol,dfcrossingmap)
  }
  
  
  print("How many km² of cropland are concerned with the potential implantation of agroforestry ? (Reforestation is not taking into account)")
  sum1agro<-sum(subset(df2All[c("croplandH2016","Allcondagrofo2")], Allcondagrofo2 == 1)$croplandH2016, na.rm = T)
  print(sum1agro)#Ici on ne prends que en compte Allcondagrofo2 et non reforestagrofo
  #En gros, ici on sait quels sont les pixels où l'espèce nommée ici espèce A
  #indiquée dans DF1plant (avec toutes les caractéristique de l'espèce A) peut pousser.
  #Et donc en sélectionnant ces pixels où l'espèce A peut pousser on peut également sélectionner leur
  #surface de culture. Du coup ici sum1agro représente la somme des surfaces de culture agricole
  #au Kenya donné par HYDE où l'on pourrait mettre l'espèce A et donc faire de l'agroforesterie
  #avec l'espèce A. ATTENTION : l'espèce A ne correspond pas forcément à la culture agricole
  #Mise en place dans le pixel, c'est dans une situation idéale.
  print("How many % it represent on the total cropland area in Kenya ?")
  print(sum1agro/sum(df2All$croplandH2016,na.rm = T)*100)
  #Ici ça représente le % de terre agricole du Kenya où on peut insérer de l'espèce A
  #étant donné que le pixel possède toutes les caractéristiques agricoles pour cette espèce.
  remove(sum1agro)
  
  
  #Selection de la bonne colonne par rapport à l'espèce précisée dans la colonne "Speciescrop"
  dfselect1 <- df2All[,grep("^SPAM",names(df2All))]
  dfselect2 <- dfselect1[,grepl(Nameplantcrop, names(dfselect1))]
  dfselect2 <- dfselect2[,-2]#remove the useless column 
  nameofSPAMprdctn <- colnames(dfselect2)
  colnames(dfselect2)[1] = "Column2"
  remove(dfselect1)
  df2All$ColumnAgrof1 <- dfselect2$Column2
  #Ici ColumnAgrof1 correspond à la colonne de SPAM donnant la quantité de prodctn de l'éspèce 
  #précisé dans Nameplantcrop en tonne
  
  print("What is the quantity of the crop specy where the agroforestry specy can grow ? (in ton)")
  somcropspecyinagrofquanti <- sum(subset(df2All[c("ColumnAgrof1","Allcondagrofo2")],Allcondagrofo2 == 1)$ColumnAgrof1, na.rm = T)
  print(somcropspecyinagrofquanti)
  #La somcropspecyinagrofquanti donne la quantité totale de crop specy (indiquée dans DF1plant et
  #provenant des données SPAM) produite au Kenya sur les pixels pouvant accueillir la "nameofplant"
  #de DF1plant en tonne
  remove(somcropspecyinagrofquanti)
  
  #Creating a column to select only pixel where there is no reforestation
  if(Algouseagrofo == 1){
    dfAgroforeforestation <- df2All$diffalgo1
  }
  if(Algouseagrofo == 2){
    dfAgroforeforestation <- df2All$diffalgo4
  }
  if(Algouseagrofo == 3){
    dfAgroforeforestation <- df2All$diffalgo42
  }
  for(i in 1:length(df2All$x)){
    df2All$Reforestagrofo[i] <- ifelse(dfAgroforeforestation[i]>0,0,1)
  }
  #La colonne Reforestagrofo permet de connaître si le pixel est reforesté (=0) ou pas
  #reforesté (= 1)
  
  #On cherche à retirer les zones reforestés car leur zone de culture peut être atteinte et donc
  #les résultats pourraient être biaisés.
  
  df2All$SPAMAgrofospecies <- df2All$ColumnAgrof1
  for(z in 1:length(df2All$x)){
    ifelse(is.na(df2All$ColumnAgrof1[z])=="FALSE",
           ifelse(is.na(df2All$Allcondagrofo2[z])=="FALSE",
                  ifelse(is.na(df2All$Reforestagrofo[z])=="FALSE",
                         ifelse(df2All$Reforestagrofo[z]==1,
                                ifelse(df2All$Allcondagrofo2[z] == 1,
                                       df2All$SPAMAgrofospecies[z] <- df2All$ColumnAgrof1[z]*RRlevel1,
                                       df2All$SPAMAgrofospecies[z] <- 0),
                                df2All$SPAMAgrofospecies[z] <- 0),
                        df2All$SPAMAgrofospecies[z] <- NA),
                  df2All$SPAMAgrofospecies[z] <- NA),
           df2All$SPAMAgrofospecies[z] <- NA)
  }
  
  #Sortir les colonnes à chaque boucle (comprise dans la boucle for ci-dessous)
  #pour vérifier que le code fonctionne bien.
  
  #Les colonnes "NameESPECE" permettent de vérifier que la sélection des colonnes a été fait correctement.

  #La colonne SPAMAgrofospecies donne la quantité de maïs(ou autre specie crop) potentielle après l'augmentation
  #du niveau de production due à l'instauration d'un système agroforestier par pixel (en tonne) dans les 
  #pixels où il est possible d'instaurer un système agroforestier selon les conditions définies avant(precip,temp,altitude..;)
  #Pour les pixels ne respectant pas les conditions pour la plante la valeur est de 0.
  
  
  print("The quantity of the crop specy produce with agroforestry is (in tonne) :")
  sumcropspecywithagroquanti <- sum(subset(df2All[c("SPAMAgrofospecies","Allcondagrofo2")],Allcondagrofo2 == 1)$SPAMAgrofospecies, na.rm = T)
  print(sumcropspecywithagroquanti)
  remove(sumcropspecywithagroquanti)
  #La sumcropspecywithagroquanti donne la quantité totale de crop specy en tonne(indiquée dans DF1plant et
  #provenant des données SPAM) produite au Kenya sur les pixels pouvant accueillir la "nameofplant"
  #de DF1plant et ayant été augmenté par l'instauration de l'agroforesterie en utilisant la "nameofplant"
  
  print("How many pixels have a new value of SPAM production (increased by agroforestry) ?
        (It needs data for SPAM Production)")
  print(length(subset(df2All[c("IDpixel","SPAMAgrofospecies")],SPAMAgrofospecies > 0)$SPAMAgrofospecies))
  #Ici on a le nombre de pixel  où l'on a augmenté la production agricole de la
  #crop specy ("crop specy" étant la même chose que la "Speciescrop" indiqué dans DF1plant)
  #avec l'instauration de la nameofplant (de DF1plant)
  
  #Dans le bloc de code ci-dessous on obtient le nom de la colonne
  #comportant les surfaces de cultures pour la crop specy indiquée dans DF1plant
  #La colonne en question provient des données de SPAM
  dfcropspecyselection <- df2All[,grep("^SPAM",names(df2All))]
  dfcropspecyselection2 <- dfcropspecyselection[,grepl(Nameplantcrop, names(dfcropspecyselection))]
  dfcropspecyselection2 <- dfcropspecyselection2[,-1]#remove the useless column 
  nameofcolspecycrop <- colnames(dfcropspecyselection2)
  remove(dfcropspecyselection,dfcropspecyselection2)
  #Le seul problème étant qu'il y a des pixels où il y a une production agricole d'indiquée 
  #mais il y a NA pour la surface agricole 
  
  print("How many km² of cropland  of the crop specy are concerned within this number of pixel (Reforestated areas removed) ?")
  dfsomme33agrof <- subset(df2All[c("SPAMAgrofospecies",nameofcolspecycrop)],SPAMAgrofospecies > 0)
  dfsomme33agrof <- dfsomme33agrof[,-1] #Ici on supprime la colonne SPAMAgrofospecies qui a 
  #permis de faire le filtre pour avoir seulement les pixels où on a mis de l'agroforesterie.
  somme33agrof <- sum(dfsomme33agrof, na.rm = T)/100
  print(somme33agrof)
  #Ici cette somme33agrof représente la surface totale de culture agricole au Kenya de la crop specy où les rendements
  #ont été augmenté par l'instauration de la nameofplant (indiqué dans DF1plant)
  remove(somme33agrof,dfsomme33agrof)
  
  
  df2All$SPAMAgrofospeciesval <- df2All$croplandH2016#création de la colonne pour accueillir
  #les données
  for(z in 1:length(df2All$x)){
    ifelse(is.na(df2All$ColumnAgrof1[z]) == "FALSE",
           ifelse(is.na(df2All$SPAMAgrofospecies[z]) == "FALSE",
                  ifelse(df2All$SPAMAgrofospecies[z] > 0,
                         df2All$SPAMAgrofospeciesval[z] <- (df2All$SPAMAgrofospecies[z]-df2All$ColumnAgrof1[z])*priceofcropspecie,
                         df2All$SPAMAgrofospeciesval[z] <- 0),
                  df2All$SPAMAgrofospeciesval[z] <- NA),
           df2All$SPAMAgrofospeciesval[z] <- NA)}
  #La colonne SPAMAgrofospeciesval donne la valeur du surplus de production engendré par la mise
  #en place d'un système agroforestier pour les pixels où c'est possible et sur la crop specy concernée
  
  
  print("Le gain engendré par la mise en place de système agroforestier 
  sur l'ensemble du territoire Kényan (où cette mise en place est possible) est de ($): ")
  valgainagrofototal <- sum(df2All$SPAMAgrofospeciesval,na.rm=T)
  print(valgainagrofototal)
  
  "Calcul du bénéfice pour l'objectif de 809km² inscrit dans les NDCs du Kenya"
  #In order to have the benefit coming from an agroforestry system in line
  #with the objective of the NDC of Kenya which is 200 000 acres (=~809km²) of 
  #agroforestry system. We have to compute the benefit for 809 km² of Agroforestry.
  
  #To do so, we need to first select randomly pixels where there is no reforestation and the
  #species can grow on it : 
  
  
  #Get the SPAM physical area depending on the specie added into Nameplantcrop in DF1plant
  dfselectother <- df2All[,grep("^SPAM",names(df2All))]
  dfselectother2 <- dfselectother[,grepl(Nameplantcrop, names(dfselectother))]
  dfselectother2 <- dfselectother2[,-1]#remove the useless column 
  nameSPAMPhysarea <- colnames(dfselectother2)
  colnames(dfselectother2)[1] = "Physarea"
  remove(dfselectother)
  #POur ce bloc de code haut dessus on sélectionne la colonne comportant pour chaque pixel
  #l'aire de culture de la specy crop précisée dans DF1plant
  
  #Obtain df where all the pixels in it are available to put agroforestry on  : 
  DfAgrofosub <- subset(df2All[c("IDpixel",nameSPAMPhysarea,nameofSPAMprdctn,"Allcondagrofo2","Reforestagrofo")],
                        Allcondagrofo2 == 1 & Reforestagrofo == 1)
  DfAgrofosub <- na.omit(DfAgrofosub)#Suppression des NA pour ne pas tirer aux
  #hasard ci-dessous, des pixels avec des NA pour SPAMProduction et Aire physique.
  
  if(is.na(DfAgrofosub$IDpixel[1])=="FALSE"){
    if((sum(DfAgrofosub[2],na.rm=T)/100)>809){
      
      while(newvalagrofo<839){
        
        #selection random of 1 pixel to put agroforestry on
        valagroforandom <- DfAgrofosub %>% sample_n(1, replace = FALSE)
        
        #Get and Keep the value of the physical area,SPAM production and ID pixel of the selected pixel
        valIDpixselected <- valagroforandom$IDpixel[1]
        Valspamphysareaselected <- dfselectother2$Physarea[valIDpixselected]
        Valspamprodselected <- dfselect2$Column2[valIDpixselected]
        
        #Pour vérification si besoin (sur les prdctns et aires sélectionnées au hasard) 
        #Cette vérification s'opère seulement sur la dernière espèce rentrée dans DF1plant
        Checkrandompixselect1 <- rbind(Checkrandompixselect1, valIDpixselected)
        Checkrandompixselect2 <- rbind(Checkrandompixselect2, Valspamphysareaselected)
        Checkrandompixselect3 <- rbind(Checkrandompixselect3, Valspamprodselected)
        
        #Remove the pixel to avoid double selection
        DfAgrofosub<-DfAgrofosub[!(DfAgrofosub$IDpixel== valIDpixselected),]
        
        #Check if the area of the pixel is not bigger than the objectiv :
        
        #La fonction ci-dessous permet d'ajouter à chaque fois l'aire du pixel sélectionné
        #de façon aléatoire ainsi que sa production à une valeur total qui ne doit pas
        #dépasser 809km² puisque c'est l'objectif des NDC au Kenya pour l'agrofo.
        #De plus, dans le code ci-dessous, lorsqu'on augmente l'aire où on met de l'agroforesterie
        #au bout d'un moment on peut dépasser ces 809km², il faut alors seulement la valeur
        #de production qui permet d'atteindre 809km² tout pile et non la valeur de production
        #qui permet d'atteindre haut delà de 809km². Par exemple, si je lance plusieurs fois
        #la sélection aléatoire de pixel, j'arrive à une newvalagrofo (réprésentant la surface
        #où l'on va mettre de l'agroforesterie) de 800km² et le prochain pixel sélectionné a une
        #surface agricole pour mettre de l'agroforesterie de 20 km². 
        #Il faut seulement la valeur de production pour ces 9 km² c'est pour cette raison qu'il y a
        #val1agrofspam.
        
        #newvalagrofo représente l'aire allouée à l'agroforesterie au total
        #AddedvalSPAMprod représente la production par l'agroforesterie
        
        if(newvalagrofo<=809){
          if((newvalagrofo+(Valspamphysareaselected/100))>809){
            val1agrofspam<-(Valspamphysareaselected/100)-((newvalagrofo+(Valspamphysareaselected/100))-809)
            newvalagrofo<-809.0001#Ici j'ajoute 809,001 dans le cas où la selectn aléatoire donne pile 809.
            AddedvalSPAMprod<-AddedvalSPAMprod+((val1agrofspam/(Valspamphysareaselected/100))*Valspamprodselected)
          }
          else{
            newvalagrofo<-newvalagrofo+(Valspamphysareaselected/100)
            AddedvalSPAMprod<-AddedvalSPAMprod+Valspamprodselected
          }
        }
        else{
          print("The total value win by agroforestry only in the objectiv of 809km² is ($) :")
          print(ifelse(RRlevel1>1,(RRlevel1-1)*AddedvalSPAMprod*priceofcropspecie,"Erreur RRlevel not >1"))
          ifelse(RRlevel1>1,
                 valgainagrofoNDC<-((RRlevel1-1)*AddedvalSPAMprod*priceofcropspecie),
                 "Erreur RRlevel not >1")
          newvalagrofo<-840}
        remove(val1agrofspam)
      } #Ici on ferme la boucle qui va tourner pour avoir une somme qui est inférieure à 809km²
      
    }#Fin de la condition permettant de vérifier s'il y a suffisamment d'espace disponible
    #de l'espece de crop où on peut mettre l'espèce d'arbre pour faire de l'agroforesterie.
    else {
      valgainagrofoNDC <- 0
      print("Il n'y a pas  assez de surface disponible pour atteindre l'objectif de 809km²
            de systèmes agroforestiers mis en place.")
    }
  }#Fin de la condition permettant de ne pas effectuer la selection aléatoire 
  #s'il n'y a pas de valeur dans DfAgrofosub, en effet, il se peut que certaines
  #espèces n'aient pas d'espace pour l'agroforesterie étant donné que l'on reforeste
  #sur plusieurs de leurs pixels (où cette espèce peut être plantée)
  else {
    valgainagrofoNDC <- 0
    print("Il n'y a pas de surface disponible pour mettre en place de l'agroforesterie
          car occupée par la reforestation.")
  }
  
  Dfbenefcroptotal <- rbind(Dfbenefcroptotal, valgainagrofototal)
  #Dfbenefcroptotal permet d'avoir le gain  ($) du surplus de production
  #engendré par l'instauration de système agroforestier à l'échelle du Kenya
  
  DfbenefcropNDC <- rbind(DfbenefcropNDC, valgainagrofoNDC)
  #DfbenefcropNDC permet d'avoir le gain ($) du surplus de production
  #engendré par l'instauration de système agroforestier dans la limite des 809km²
  
  #In order to have the column for each species included in DF1plant: 
  
  Changingname1 <- paste("Allcondagrofo",Nameplantcrop,sep = "_")
  colnames(df2All)[colnames(df2All) == "Allcondagrofo"] = Changingname1 
  df2All <- df2All[ , !(names(df2All) %in% c("Allcondagrofo"))]
  #This change in the name allow to have for each species the column Allcondagrofo
  #which correspond for a value = 1 of all the pixel where the specific agroforestry specie
  #Can grow. If the value is = 0 the specific agroforestry specie can't grow. Here
  #Allcondagrofo don't take into account the soil conditions.
  
  Changingname2 <- paste("Allcondagrofo2",Nameplantcrop,sep = "_")
  colnames(df2All)[colnames(df2All) == "Allcondagrofo2"] = Changingname2
  df2All <- df2All[ , !(names(df2All) %in% c("Allcondagrofo2"))]
  #Allcondagrofo give (if 1 in value) where the specific agroforestry specie can grow.
  #If the value is = 0 the specific agroforestry specie can't grow. Here if there are infor
  #-mation on the soil available, this column take into account the soil requirements.
  
  Changingname3 <- paste("SPAMAgrofospecies",Nameplantcrop,sep = "_")
  colnames(df2All)[colnames(df2All) == "SPAMAgrofospecies"] = Changingname3
  #La colonne SPAMAgrofo permet de connaître la production totale par pixel en
  #ayant appliqué de l'agroforesterie dans ce pixel et donc en ayant augmenté
  #les rendements.
  
  Changingname4 <- paste("SPAMAgrofospeciesval",Nameplantcrop,sep = "_")
  colnames(df2All)[colnames(df2All) == "SPAMAgrofospeciesval"] = Changingname4
  #La colonne SPAMAgrofoval permet de connaître pour chaque pixel, le bénéfice
  #créé en plus par l'agroforesterie en dollars(donc vraiment juste le surplus de production
  #produit par l'agroforesterie et valorisé).
  
  #remove column for the loop : 
  df2All <- df2All[ , !(names(df2All) %in% c("SPAMAgrofospecies","SPAMAgrofospeciesval"))]
  df2All <- df2All[ , !(names(df2All) %in% c("ColumnAgrof1"))]
  
  #Pour vérification concernant les colonnes SPAM sélectionnées:
  Checking_namecol <- rbind(Checking_namecol, nameofSPAMprdctn)
  
  
  remove(nameofSPAMprdctn,Changingname4,Changingname3,Changingname2,Changingname1,
         nameSPAMPhysarea,AddedvalSPAMprod,valIDpixselected,Valspamphysareaselected,
         Valspamprodselected,valagroforandom,dfselect2,dfselectother2,valgainagrofototal,
         valgainagrofototal)
  
  print("-----------------------------------------")
  print("----------NEXT--------------------------")
  
}
Checking_namecol <- Checking_namecol[-c(1,2),]#remove the useless lines
DfbenefcropNDC <- DfbenefcropNDC[-c(1,2,3),]#remove the useless lines
Dfbenefcroptotal <- Dfbenefcroptotal[-c(1,2,3),]#remove the useless lines
Checkrandompixselect1 <- Checkrandompixselect1[-1,]
Checkrandompixselect2 <- Checkrandompixselect2[-1,]
Checkrandompixselect3 <- Checkrandompixselect3[-1,]


#Vérification de la sélection des pixels aléatoires pour la dernière espèce données dans DF1plant: 
Checkrandompixselect <- data.frame(Checkrandompixselect1,Checkrandompixselect2,
                                   Checkrandompixselect3)
colnames(Checkrandompixselect) = c("IDpix","SPAMPhysarea","SPAMProd")
sum(Checkrandompixselect$SPAMPhysarea)/100#Si haut-dessus de 899, il y a un problème
"Pour la somme ci-dessous il faut préciser le prix de l'espèce de crop indiqué dans DF1plant"
sum(Checkrandompixselect$SPAMProd)*165.5*(RRlevel1-1) #doit être à peu près équivalent à 
#"valgainagrofoNDC" , sauf que dans la somme juste haut-dessus on prends en compte l'entièreté
#de la production du dernier pixel sélectionné dans Checkrandompixselect de plus il y en a
#un pixel en + qui est ajouté dans cette somme par rapport au calcul de gain en respectant
#l'objectif de la NDC, car pour le tout dernier pixel sélectionné au hasard, qu'il vienne ajouter
#ou pas le niveau d'atteinte des 809km² défini par le NDC, il est sélectioné et ajouté 
#dans Checkrandompixelselect.

### Launching ending point -- -- -- -- ----



#Information supplémentaire : 
#Pour l'obtention du gain Dfbenefcroptotal, c'est seulement si l'on transforme toutes
# les aires de production de namespeciecrop en agroforesterie et en posant l'hypothèse
#qu'il n'y avait pas de système agroforestier similaire avant. 

#(((( Parfois la production (dans SPAM) de différentes espèces crop est la même mais ça représente
#qu'un tout petit pourcentage sur 6931 pixels. Par exemple entre MAIZ et COTT, il y a 316
#pixels où ces deux cultures ont la même production. )))))




#4 - Cleaning/Nettoyage ----
remove(DF1plant,DF2plant,Checking_namecol,dfselect1,dfselect2,dfselectother,
       dfselectother2,Test23,DfAgrofosub)
df2All <- df2All[ , !(names(df2All) %in% c("Allcondagrofo","Allcondagrofo2",
                                           "ColumnAgrof1","Reforestagrofo"))]
#It is needed, before launch again the loop to remove the column of the species as
#SPAMAgrofospecies_nameofthespecie and SPAMAgrofospeciesval_nameofthespecie and
# Allcondagrofo2_nameofthespecie and Allcondagrofo_nameofthespecie, for example for MAIZ and COTT:
df2All <- df2All[ , !(names(df2All) %in% c("SPAMAgrofospecies_COTT","SPAMAgrofospecies_MAIZ"))]
df2All <- df2All[ , !(names(df2All) %in% c("SPAMAgrofospeciesval_COTT","SPAMAgrofospeciesval_MAIZ"))]
df2All <- df2All[ , !(names(df2All) %in% c("Allcondagrofo_COTT","Allcondagrofo_MAIZ",
                                           "Allcondagrofo2_COTT","Allcondagrofo2_MAIZ"))]
df2All <- df2All[ , !(names(df2All) %in% c("Allcondagrofo"))]
df2All <- df2All[ , !(names(df2All) %in% c("SPAMAgrofospecies","SPAMAgrofospeciesval"))]
remove(Checking_namecol,Checkrandompixselect1,Checkrandompixselect2,Checkrandompixselect3)
remove(DF1plant,DF2plant,DfbenefcropNDC,Dfbenefcroptotal,dfcheckingsoil)
remove(dfcheckingsoil,dfcrossingmap,Checkrandompixselect)
remove(Shpsoiljoin)
remove(valagroforandom, Dfbenefcroptotal, DfbenefcropNDC)
remove(Changingname1,Changingname2,Changingname3,Changingname4)
remove(dfAgroforeforestation,compteurpromp1,compteurpromp2,Havethemap,
       Algouseagrofo,nameofSPAMprdctn,Nameplant1,
       Nameplantcrop,nameSPAM,newvalagrofo,priceofcropspecie,RRlevel1,AddedvalSPAMprod)
remove(list_raster,list_dfofraster,NasaDEM)
remove(rastbatprecitemp,AddedvalSPAMprod)
remove(NameCOTT,NameMAIZ)#Remove the NameSpecie df
#SI besoin de relancer tout :
df2All <- df2All[ , !(names(df2All) %in% c("diffalgo1","diffalgo4","diffalgo42"))]
df2All <- df2All[ , !(names(df2All) %in% c("Reforestagrofo"))]
remove(DfAgrofosub)




#Pour vérifier les valeurs il est possible de sortir df2All après lancement des boucles.
#Attention à bien supprimer s'il y a déjà une base existante.
write_xlsx(df2All,"Data/df2All1.xlsx")

