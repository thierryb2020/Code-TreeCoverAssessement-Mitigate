
#= = = = = = = = = = = = = = = =
#- - - - - - - - - - - - - - - -
#       PROJET MITIGATE +
#- - - - - - - - - - - - - - - -
#= = = = = = = = = = = = = = = =

#Code created by Paul Bostyn the 26/10/2023 under the behalf of the Mitigate+
#Project with the CGIAR 


"Pour une utilisation optimale du code, l'utilisateur doit lancer partie par partie
en lisant les instructions et explications propres à chaque partie."

#  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
# I - Importation et mise en forme de la base de données ----
#  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -


## 1 - Liste des packages à installer et chargement  ----
library(tidyverse)
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
library(car)
library(readxl)
library(VIF)
library(stargazer)
library(lwgeom)

library(rgeos)
library(rgdal)

## 2 - Définition de l'origine de l'importation    ----

setwd("D:/Mes Donnees/Projets/CGIAR - Mitigate+/Modèle CostKenya")
"Définir l'emplacement du fichier dans lequel se trouve l'ensemble des rasters"
emplacementfolder <- "Data/All"
"Définir la source du fichier dans lequel se trouve l'ensemble des rasters"
origin <- "Data/All/"
df2All <- read_excel("Data/df2All.xlsx")

"Possibilité de définir le CRS commun et la résolution pour chaque raster"
CRS <- 1
RESOL <- 1
"Dans le code ci-dessous le CRS et la résolution des rasters sont basés
sur celui du premier raster."

"Définition des valeurs pour les rasters à assembler (pour que tous les raster aient la même zone d'étude) :"
#Premier groupe de raster
numrastassembl1min <- 1 # numéro du premier raster à assembler en Rast 1 dans list_rast
numrastassembl1max <- 4 # numéro du dernier raster à assembler en Rast 1 dans list_rast
name1 <- "GFWKenya"
#Plus d'informations voir Annexe 1

"Possibilité de définir un territoire de référence auquel les rasters seront découpés."
Rastdecoupe <- read_sf("Data/Data_Kenyaboundary/kenboundary.shp")
crsref <-crs(Rastdecoupe)

"Pour les autres rasters (ici ceux pour obtenir la moyenne des précipitations sur la période
2017-2021 inclus) :"
"Définir l'emplacement du fichier"
emplacementfolder31 <- "Data/Worldclimimport"
"Définir la source du fichier"
origin31 <- "Data/Worldclimimport/"


## 3 - Importation des rasters  ----
files <- list.files(emplacementfolder)#donne un vecteur sous forme de liste avec tous les noms des fichier dans le document choisi
list_rast <- list() #création de la liste pour stocker les rasters
for (i in 1:length(files)){
  name_file <- files[[i]]#ici ça assigne le nom au i ème élément de la list "files"
  real_name <- paste(origin,name_file, sep="") #ici la fonction paste permet de combiner le nom de l element i de la ligne haut dessus
  #et la source du fichier pour avoir le vrai nom de la source du fichier et le sep c'est pour qu'il n'y ait pas d'espace
  list_rast <- append(list_rast,raster(real_name))#ici ça ajoute chaque nouveau ieme element dans list_rast par la fonction append
  names(list_rast)[i] <- paste("R", i, sep = "")#ici le code renomme chaque element ajouté dans list_rast
}
remove(files, i, name_file, origin, real_name)


## 4 - Modification du CRS et de la résolution des rasters ----
### a) - CRS ----
for (i in 1:length(list_rast)) {
  crs(list_rast[[i]]) <- crs(Rastdecoupe)
}
remove(i)
"list_rast représente l'ensemble des rasters du fichiers avec un CRS commun"


### b) - Résolution ----
"Définition du groupe de raster à modifier de la même façon pour la résolution"
"Modification de la résolution des rasters de World CLim"
numdebutgrouprast <- 20
numfingrouprast <- 43
for(i in numdebutgrouprast:numfingrouprast){
  list_rast[[i]] <- aggregate(list_rast[[i]],fact = 2, fun = mean)
}
res(list_rast[[22]])


## 5 - Assemblage de certains rasters ----
#(Information supplémentaire cf Annexe 1)

### a) Global Forest Watch ----
listGFW <- list()
listGFW <- list_rast[c(numrastassembl1min:numrastassembl1max)] #extraire les rasters du premier groupe
#de raster de list_rastreso

#Modification de la résolution des rasters de GFW 
"La résolution de référence est celle de SPAM et HYDE, mais cette résolution de 0.0833333 pose problème
donc il est nécessaire de passer par le formatage ci-dessous pour que les rasters de GFW ayant une
résolution initiale de 0.00025 arrivent à une résolution ~ de 0.083333 pour qu'on puisse ensuite les
fusionner ensemble pour obtenir une df."
"Information : Code prenant du temps :"
for(i in 1:length(listGFW)) {
  listGFW[[i]] <- aggregate(listGFW[[i]],fact = 333, fun = mean)
}
"ATTENTION AU RESAMPLE POUVANT ALTERER LA DISPOSITION DES PIXELS
donc vérifier en regardant sur un logiciel SIG"
for(i in 1:length(listGFW)){
  listGFW[[i]] <- resample(listGFW[[i]], list_rast[[6]])
}
#La référence de list_rast 6 ici correspond à la résolution de HYDE

#Pour assembler les rasters en un seul df : 
listGFWrast <- listGFW
listGFW <-lapply(listGFW, function(x) {as.data.frame(rasterToPoints(x,fun=function(x){x>=0}))})
#listGFW devient ici la liste des dataframes obtenus à partir des rasters du premier groupe de raster
for (i in 1:length(listGFW)) {
  names(listGFW[[i]]) <- c("x", "y", name1)#renomme les colonnes des df du premier groupe de raster pour les fusionner
}

GFWdf <- listGFW[[1]] #Définition d'un des raster pour assembler les autres à celui-ci
for (i in 2:length(listGFW)) {
  GFWdf <- bind_rows(listGFW[[i]],GFWdf) #assemblage des rasters
}


GFWrast <- rasterFromXYZ(GFWdf)#modification du dataframe de GFW en raster 

#Supprimer les anciens rasters du premier groupe de raster de la liste :
list_rast <- list_rast[- c(numrastassembl1min:numrastassembl1max)] #ici

#Ajouter le raster assemblé de GFW : 
list_rast <- append(list_rast, GFWrast)
#il se positionne en dernier dans list_rastreso


## 6 - Découpe des rasters au territoire de référence ----

list_rastreso2<-list()
for (i in 1:length(list_rast)){
  list_rastreso2 <- append(list_rastreso2, crop(list_rast[[i]],Rastdecoupe))
}#crop() permet de découper les rasters de la liste par rapport à Rastdecoupe

#permet de resample tous les raster à la résolution et à la forme de HYDE
#sauf les globalforestwatch
#!! /!\ ATTENTION AU resample qui peut altérer la disposition des pixels --> 
#besoin d'une vérification en observant le raster initial et après traitement
for(i in 1:(length(list_rastreso2))){
  list_rastreso2[[i]] <- resample(list_rastreso2[[i]], list_rast[[6]])
}

#Sortir le raster de GFW
writeRaster(list_rastreso2[[82]], "rastGFW2", format = "GTiff")#le raster sera exporté
#dans le document de R_studio

"list_rastreso2 regroupe tous les rasters avec la même résolution, même CRS et découpage à la 
forme du Kenya."


## 7 - Vérification des rasters (conformité,résolution..) ----
res(list_rastreso2[[1]])#HYDE
res(list_rastreso2[[15]])#SPAM_valeur
res(list_rastreso2[[20]])#World CLIM
res(list_rastreso2[[78]])#SPAM_prod
res(list_rastreso2[[82]])#GFW
extent(list_rastreso2[[1]])#HYDE
extent(list_rastreso2[[15]])#SPAM_valeur
extent(list_rastreso2[[20]])#World CLIM
extent(list_rastreso2[[78]])#SPAM_prod
extent(list_rastreso2[[82]])#GFW

remove(GFWdf,listGFW,listGFWrast,rastsurfacehyde,list_rast)


## 8 - Transformation des rasters en données vecteur et fusion en dataframe ----

### a) Conversion des rasters en vecteur ----
#Pour avoir une liste des rasters en format dataframe :
listdf_rast <-lapply(list_rastreso2, function(x) {as.data.frame(rasterToPoints(x, fun = function(x){x>=0}))})#le as data frame les transforme dans la fonction })


### b) Jointure des données ----

#Création d'un premier df découper de référence
dfshpfile <- listdf_rast[[1]]
coordinates(dfshpfile)=~x+y
proj4string(dfshpfile) <- proj4string(Rastdecoupe)
dfshpfile <- dfshpfile[Rastdecoupe,]
df2All <- as.data.frame(dfshpfile)
df2All$x <- round(x = df2All$x, 10)
df2All$y <- round(x = df2All$y, 10)
names(df2All)[3] <- "Rast1"
remove(dfshpfile)
#Les x et y des rasters sont légèrement arrondis pour une bonne concordance.

for (i in 2:length(listdf_rast)){
  dfshpfile <- listdf_rast[[i]]
  coordinates(dfshpfile)=~x+y
  proj4string(dfshpfile) <- proj4string(Rastdecoupe)
  dfshpfile <- dfshpfile[Rastdecoupe,]
  dfshpfile <- as.data.frame(dfshpfile)
  nam <- paste("Rast", i , sep = "")#donne un nom Rast1 au premier elemnt de la listdf_rast
  names(dfshpfile)[3] <- nam
  dfshpfile$x <- round(x = dfshpfile$x, 10)
  dfshpfile$y <- round(x = dfshpfile$y, 10)
  df2All <- full_join(df2All, dfshpfile, by = c("x","y"))
  remove(dfshpfile,nam)
}

### c) Exportation des rasters pour vérifier ----
#Ici l'objectif est de re-transformer les rasters (que l'on a déjà
#transformés en dataframe) en raster et vérifier dans un logiciel SIG
#qu'il n'y a pas eu de modification.

#Export d'un raster de HYDE
dfHYDE <- subset(df2All[c("x","y","Rast1")], )
coordinates(dfHYDE)=~x+y
proj4string(dfHYDE) <- proj4string(Rastdecoupe)
writeOGR(dfHYDE, dsn = '.', layer = 'hyde', driver = "ESRI Shapefile")
#Export de SPAM 
dfSPAM <- subset(df2All[c("x","y","Rast15")], )
coordinates(dfSPAM)=~x+y
proj4string(dfSPAM) <- proj4string(Rastdecoupe)
writeOGR(dfSPAM, dsn = '.', layer = 'SPAM', driver = "ESRI Shapefile")
#Export de World Clim
dfWC <- subset(df2All[c("x","y","Rast20")], )
coordinates(dfWC)=~x+y
proj4string(dfWC) <- proj4string(Rastdecoupe)
writeOGR(dfWC, dsn = '.', layer = 'WlrdCLim', driver = "ESRI Shapefile")
#Export de GFW
dfGFW <- subset(df2All[c("x","y","Rast82")], )
coordinates(dfGFW)=~x+y
proj4string(dfGFW) <- proj4string(Rastdecoupe)
writeOGR(dfGFW, dsn = '.', layer = 'GFW', driver = "ESRI Shapefile")



## 9 - Importation des autres types de rasters et fusion ----

### a) Accessibilité ----
New_file_rastaccess <- nc_open("mkt_access_5m.nc")#Si l'importation ne fonctionne
#pas, il est nécessaire de préciser l'endroit. Ce fichier se trouve dans le dossier de R_studio.
mkt_access_tmp1 <- ncvar_get(New_file_rastaccess, varid = "Band1")
RastAccess <- t(raster(mkt_access_tmp1))
res(RastAccess)
remove(New_file_rastaccess,mkt_access_tmp1)
extent(RastAccess) <- extent(list_rastreso2[[1]])
"modifier l'extent en premier pour pouvoir verifier les df, la modification
de l'extent permet de mettre meme resolution"
crs(RastAccess) <- crs(Rastdecoupe)
RastAccess <- crop(RastAccess,Rastdecoupe)
dftest2 <- as.data.frame(rasterToPoints(RastAccess))
coordinates(dftest2)=~x+y
proj4string(dftest2) <- proj4string(Rastdecoupe)
dfRastAccess <- dftest2[Rastdecoupe,]
dfRastAccess <- as.data.frame(dfRastAccess)
dfRastAccess$x <- round(x = dfRastAccess$x, 10)
dfRastAccess$y <- round(x = dfRastAccess$y, 10)
remove(dftest2,RastAccess)

df2All <- full_join(df2All, dfRastAccess, by = c("x","y"))
remove(dfRastAccess)


### b) World Clim Précipitation ----

#Permet de créer les colonnes WCprecimean2017, WCprecimean2018, WCprecimean2019
# WCprecimean2020, WCprecimean2021 et WCprecimeanall dans df2All

#1 - Importation des rasters  
files31 <- list.files(emplacementfolder31)
list_rast31 <- list() 
for (i in 1:length(files31)){
  name_file31 <- files31[[i]]
  real_name31 <- paste(origin31,name_file31, sep="") 
  list_rast31 <- append(list_rast31,raster(real_name31))
  names(list_rast31)[i] <- paste("R", i, sep = "")
}
remove(files31, i, name_file31, origin31, real_name31)
#2 - Modification du CRS et de la résolution des rasters
for (i in 1:length(list_rast31)) {
  crs(list_rast31[[i]]) <- crs(Rastdecoupe)
}
remove(i)
#Modification de la résolution
for(i in 2:49){
  list_rast31[[i]] <- aggregate(list_rast31[[i]],fact = 2, fun = mean)
}
res(list_rast31[[22]])#Vérification si résolution même que HYDE
res(list_rast31[[1]])
extent(list_rast31[[22]])#Vérification si extent même que HYDE
extent(list_rast31[[1]])
"découper à la forme du kenya"
list_rastreso31<-list()
for (i in 1:length(list_rast31)){
  list_rastreso31 <- append(list_rastreso31, crop(list_rast31[[i]],Rastdecoupe))
}
"Resample pr avoir même resolutn ect"
for(i in 2:length(list_rastreso31)){
  list_rastreso31[[i]] <- resample(list_rastreso31[[i]], list_rastreso31[[1]])
}
res(list_rastreso31[[1]])
res(list_rastreso31[[22]])
extent(list_rastreso31[[1]])
extent(list_rastreso31[[22]])
"supprimer HYDE qui servait avant de référence"
list_rastreso31 <- list_rastreso31[-1]
"Transformation en df"
listdf_rast31 <-lapply(list_rastreso31, function(x) {as.data.frame(rasterToPoints(x, fun = function(x){x>=0}))})
"1 df sert de référence"
dfWClim1 <- listdf_rast31[[1]]
coordinates(dfWClim1)=~x+y
proj4string(dfWClim1) <- proj4string(Rastdecoupe)
dfWClim1 <- dfWClim1[Rastdecoupe,]
dfWClim2 <- as.data.frame(dfWClim1)
dfWClim2$x <- round(x = dfWClim2$x, 10)
dfWClim2$y <- round(x = dfWClim2$y, 10)
names(dfWClim2)[3] <- "Rast1"
remove(dfWClim1)
"Appliquer ces transformations à tout les rasters"
for (i in 2:length(listdf_rast31)){
  dfWClim1 <- listdf_rast31[[i]]
  coordinates(dfWClim1)=~x+y
  proj4string(dfWClim1) <- proj4string(Rastdecoupe)
  dfWClim1 <- dfWClim1[Rastdecoupe,]
  dfWClim1 <- as.data.frame(dfWClim1)
  nam <- paste("Rast", i , sep = "")#donne un nom Rast1 au premier elemnt de la listdf_rast
  names(dfWClim1)[3] <- nam
  dfWClim1$x <- round(x = dfWClim1$x, 10)
  dfWClim1$y <- round(x = dfWClim1$y, 10)
  dfWClim2 <- full_join(dfWClim2, dfWClim1, by = c("x","y"))
  remove(dfWClim1,nam)
}
remove(list_rast31,list_rastreso31,listdf_rast31)

"--renommer les colonnes--"
dfWClim2 <- dfWClim2 %>% rename(
  WCpreci012017 = Rast1, WCpreci022017 = Rast2, 
  WCpreci032017 = Rast3, WCpreci042017 = Rast4, WCpreci052017 = Rast5,
  WCpreci062017 = Rast6, WCpreci072017 = Rast7, WCpreci082017 = Rast8,
  WCpreci092017 = Rast9, WCpreci102017 = Rast10, WCpreci112017 = Rast11,
  WCpreci122017 = Rast12, WCpreci012018 = Rast13, WCpreci022018 = Rast14,
  WCpreci032018 = Rast15, 
  WCpreci042018 = Rast16, WCpreci052018 = Rast17,
  WCpreci062018 = Rast18, WCpreci072018 = Rast19,
  WCpreci082018 = Rast20, WCpreci092018 = Rast21,
  WCpreci102018 = Rast22, WCpreci112018 = Rast23,
  WCpreci122018 = Rast24, WCpreci012019 = Rast25,
  WCpreci022019 = Rast26, WCpreci032019 = Rast27,
  WCpreci042019 = Rast28, WCpreci052019 = Rast29,
  WCpreci062019 = Rast30, WCpreci072019 = Rast31,
  WCpreci082019 = Rast32, WCpreci092019 = Rast33,
  WCpreci102019 = Rast34, WCpreci112019 = Rast35,
  WCpreci122019 = Rast36, WCpreci012020 = Rast37,
  WCpreci022020 = Rast38, WCpreci032020 = Rast39, WCpreci042020 = Rast40,
  WCpreci052020 = Rast41, WCpreci062020 = Rast42, 
  WCpreci072020 = Rast43, WCpreci082020 = Rast44, WCpreci092020 = Rast45,
  WCpreci102020 = Rast46, WCpreci112020 = Rast47, WCpreci122020 = Rast48,
)
dfWClim2$WCprecimean2017 <- dfWClim2$WCpreci012017
dfWClim2$WCprecimean2018 <- dfWClim2$WCprecimean2017
dfWClim2$WCprecimean2019 <- dfWClim2$WCprecimean2017
dfWClim2$WCprecimean2020 <- dfWClim2$WCprecimean2017
for(i in 1:length(dfWClim2$x)){
  dfWClim2$WCprecimean2017[i]<- (dfWClim2$WCpreci012017[i]+dfWClim2$WCpreci022017[i]+
                                   dfWClim2$WCpreci032017[i]+dfWClim2$WCpreci042017[i]+
                                   dfWClim2$WCpreci052017[i]+dfWClim2$WCpreci062017[i]+
                                   dfWClim2$WCpreci072017[i]+dfWClim2$WCpreci082017[i]+
                                   dfWClim2$WCpreci092017[i]+dfWClim2$WCpreci102017[i]+dfWClim2$WCpreci112017[i]+
                                   dfWClim2$WCpreci122017[i])/12
}
for(i in 1:length(dfWClim2$x)){
  dfWClim2$WCprecimean2018[i]<- (dfWClim2$WCpreci012018[i]+dfWClim2$WCpreci022018[i]+
                                   dfWClim2$WCpreci032018[i]+dfWClim2$WCpreci042018[i]+
                                   dfWClim2$WCpreci052018[i]+dfWClim2$WCpreci062018[i]+
                                   dfWClim2$WCpreci072018[i]+dfWClim2$WCpreci082018[i]+
                                   dfWClim2$WCpreci092018[i]+dfWClim2$WCpreci102018[i]+dfWClim2$WCpreci112018[i]+
                                   dfWClim2$WCpreci122018[i])/12
}
for(i in 1:length(dfWClim2$x)){
  dfWClim2$WCprecimean2019[i]<- (dfWClim2$WCpreci012019[i]+dfWClim2$WCpreci022019[i]+
                                   dfWClim2$WCpreci032019[i]+dfWClim2$WCpreci042019[i]+
                                   dfWClim2$WCpreci052019[i]+dfWClim2$WCpreci062019[i]+
                                   dfWClim2$WCpreci072019[i]+dfWClim2$WCpreci082019[i]+
                                   dfWClim2$WCpreci092019[i]+dfWClim2$WCpreci102019[i]+dfWClim2$WCpreci112019[i]+
                                   dfWClim2$WCpreci122019[i])/12
}
for(i in 1:length(dfWClim2$x)){
  dfWClim2$WCprecimean2020[i]<- (dfWClim2$WCpreci012020[i]+dfWClim2$WCpreci022020[i]+
                                   dfWClim2$WCpreci032020[i]+dfWClim2$WCpreci042020[i]+
                                   dfWClim2$WCpreci052020[i]+dfWClim2$WCpreci062020[i]+
                                   dfWClim2$WCpreci072020[i]+dfWClim2$WCpreci082020[i]+
                                   dfWClim2$WCpreci092020[i]+dfWClim2$WCpreci102020[i]+dfWClim2$WCpreci112020[i]+
                                   dfWClim2$WCpreci122020[i])/12
}
dfWClimmean <- subset(dfWClim2[c("x","y","WCprecimean2020","WCprecimean2019",
                                 "WCprecimean2017","WCprecimean2018")])
df2All <- full_join(df2All, dfWClimmean, by = c("x","y"))
remove(dfWClim2,dfWClimmean)


### c) Données Agro-climatiques ----
#L'objectif avec cette import est pour la régression précipitation ~ couverture arborée.
#Pour contrôler par type de biome. 

#Source des données : https://geoportal.icpac.net/layers/geonode%3Aken_aczones
"1) -- Importer et avoir le shapefile de référence "
#Import d'un shpfile contenant les différents biomes en se centrant sur les agroclimatic zone: 
shape <- readOGR(dsn = "Data/Data_agroclim/ken_aczones.shp", layer = "ken_aczones")
#ajout colonne df2All
df2All$IDbiome <- c(1:length(df2All$GFW2010))
"2) -- avoir le format shapefile de df2All"
dftest1 <- df2All[,c("x","y","GFW2010","IDbiome")]
coordinates(dftest1) <- ~ x + y
proj4string(dftest1) <- proj4string(shape)
dftest1 <- st_as_sf(dftest1)
shape <- st_as_sf(shape)
"3) -- pour voir les intersection entre les deux :"
out <- st_intersection(dftest1,shape)
datafrme1 <- as.data.frame(out)
"4) -- Insérer le type de biome dans la base "
for(i in 1:length(datafrme1$IDbiome)){
  valueindix <- datafrme1$IDbiome[i]
  valeurKACZ <- datafrme1$KACZ_[i]
  df2All$Biometype[valueindix]<-valeurKACZ
}
"5) -- Insérer des NA lorsqu'il n'y a pas de valeur "
df2All$Biometype[1]#vérification de la composition des cellules vides
for(i in 1:length(df2All$GFW2010)){
  df2All$Biometype[i] <- ifelse(df2All$Biometype[i]=="",
                                NA,df2All$Biometype[i])}
"6) -- Vérification de la conversion "
dftest <- subset(df2All[c("x","y","Biometype")])
coordinates(dftest) <- ~ x + y
proj4string(dftest) <- proj4string(Rastdecoupe)
dftest <- st_as_sf(dftest)
#Extract le shapefile
st_write(dftest, "Classmaps02.shp")#Mettre le shapefile dans QGIS et comparer
#Avec le qgis initial ken_aczones si besoin d'une vérification


### d) Séquestration carbone ----
# -- Données de séquestration carbone naturelle par reforestation naturelle (Cook-Patton et al. 2020 -
#Mapping carbon accumulation potential from global natural forest regrowth)

"Importation du raster"
rastGFWRegrow <- raster("Data/Data_cookpatton/regrowthcookpatton.tif")
rastref <- raster("Data/Data_rastref/grazing2016AD.asc")
"Modification de sa résolution et de son crs"
newproj1 <- crs(Rastdecoupe)
rastGFWRegrow<-projectRaster(rastGFWRegrow,crs = newproj1)
rastGFWRegrow <- resample(rastGFWRegrow, rastref) #ATTENTION AU RESAMPLE POUVANT
#ALTERER LA DISPOSITION DES PIXELS ET LEUR VALEUR (ça modifie la résolution donc
#la taille d'un pixel longueurxlargeur)
extent(rastref)
extent(rastGFWRegrow)
crs(rastGFWRegrow)
crs(Rastdecoupe)
res(rastGFWRegrow)
res(rastref)
"Transformation du raster en shp et en df et découpe taille Kenya"
rastGFWRegrow <- rasterToPoints(rastGFWRegrow,fun=function(x){x>=0})
dfRastGFWregrow <- as.data.frame(rastGFWRegrow)
coordinates(dfRastGFWregrow)=~x+y
proj4string(dfRastGFWregrow) <- proj4string(Rastdecoupe)
dfRastGFWregrow <- dfRastGFWregrow[Rastdecoupe,]
shpRastgfw <- dfRastGFWregrow
dfRastGFWregrow <- as.data.frame(dfRastGFWregrow)
"Vérification de la bonne transformation du raster sur logiciel SIG"
verifrast <- rasterFromXYZ(dfRastGFWregrow)
writeRaster(verifrast, "verifrastregrowthcookpat", format = "GTiff")#Ici un raster sort et peut être mis dans
#un logiciel SIG pour vérifier si le resample fait haut-dessus est correct.
"Somme de la métrique :"
sum(dfRastGFWregrow$regrowthcookpatton, na.rm = T)
remove(rastref,rastGFWRegrow,dfRastGFWregrow,shpRastgfw,verifrast,newproj1)
#Croisement des données entre df2All et dfRastGFWregrow
dfRastGFWregrow$x <- round(x = dfRastGFWregrow$x, 10)
dfRastGFWregrow$y <- round(x = dfRastGFWregrow$y, 10)
df2All <- full_join(df2All, dfRastGFWregrow, by = c("x","y"))
remove(CookPatton,dfRastGFWregrow,dftest,newproj1,
       rastGFWRegrow,rastGFWtest,rastref,shpRastgfw)
#Ici la colonne regrowthcookpatton représente le niveau de Mg carbon/ha/year
#pouvant être séquestré en laissant les forêts du pixels repousser naturellement
#pour les années 2020-2050 ici  1 Mg = 1 tonne car c'est dit dans l'article que
#ce sont des mégagrammes.

### e) Counties Kenya ----
shapecunties <- st_read("Data/Data_Kenyacounties/ke_county.shp")%>%
  st_transform(st_crs(Rastdecoupe))#Same crs is applied


## 10 - Mise en forme des données ----

#Modification des noms de colonne
df2All <- df2All %>% rename(Accessi = layer)
df2All <- df2All %>% rename(
    croplandH2016 = Rast1, grazingH2016 = Rast2, 
    irnoriceH2016 = Rast3, irriceH2016 = Rast4, landlakeH2016 = Rast5,
    maxlnH2016 = Rast6, pastureH2016 = Rast7, rangelandH2016 = Rast8,
    rfnoriceH2016 = Rast9, rfriceH2016 = Rast10, totirriH2016 = Rast11,
    totrainfedH2016 = Rast12, totriceH2016 = Rast13, uoppH2016 = Rast14,
    SPAMProdvalha2010 = Rast15, 
    WCpreci012021 = Rast16, WCpreci022021 = Rast17,
    WCpreci032021 = Rast18, WCpreci042021 = Rast19,
    WCpreci052021 = Rast20, WCpreci062021 = Rast21,
    WCpreci072021 = Rast22, WCpreci082021 = Rast23,
    WCpreci092021 = Rast24, WCpreci102021 = Rast25,
    WCpreci112021 = Rast26, WCpreci122021 = Rast27,
    WCtmax012021 = Rast28, WCtmax022021 = Rast29,
    WCtmax032021 = Rast30, WCtmax042021 = Rast31,
    WCtmax052021 = Rast32, WCtmax062021 = Rast33,
    WCtmax072021 = Rast34, WCtmax082021 = Rast35,
    WCtmax092021 = Rast36, WCtmax102021 = Rast37,
    WCtmax112021 = Rast38, WCtmax122021 = Rast39,
    SPAMACOFP = Rast40, SPAMBANAP = Rast41, 
    SPAMBARLP = Rast42, SPAMBEANP = Rast43, SPAMCASSP = Rast44,
    SPAMCHICP = Rast45, SPAMCNUTP = Rast46, SPAMCOCOP = Rast47,
    SPAMCOTTP = Rast48, SPAMCOWPP = Rast49, SPAMGROUP = Rast50,
    SPAMLENTP = Rast51, SPAMMAIZP = Rast52, SPAMOCERP = Rast53,
    SPAMOFIBP = Rast54, SPAMOILPP = Rast55, SPAMOOILP = Rast56,
    SPAMOPULP = Rast57, SPAMORTSP = Rast58, SPAMPIGEP = Rast59,
    SPAMPLNTP = Rast60, SPAMPMILP = Rast61, SPAMPOTAP = Rast62,
    SPAMRAPEP = Rast63, SPAMRCOFP = Rast64,
    SPAMRESTP = Rast65, SPAMRICEP = Rast66,
    SPAMSESAP = Rast67, SPAMSMILP = Rast68,
    SPAMSORGP = Rast69, SPAMSOYBP = Rast70,
    SPAMSUGBP = Rast71, SPAMSUGCP = Rast72,
    SPAMSUNFP = Rast73, SPAMSWPOP = Rast74,
    SPAMTEASP = Rast75, SPAMTEMFP = Rast76,
    SPAMTOBAP = Rast77, SPAMTROFP = Rast78, SPAMVEGEP = Rast79,
    SPAMWHEAP = Rast80, SPAMYAMSP = Rast81, GFW2010 = Rast82)
stop(print("Il est très important que les colonnes aient bien le bon nom du fichier dans All"))
#En effet, par exemple, dans le fichier All, cropland est le 5eme raster dans l'ordre, mais
#en enlevant les 4 premiers rasters de Global Forest Watch - GFW, il devient le premier, donc
#pour vérifier il faut regarder dans le fichier All à quel niveau se situe le raster et faire - 4
#pour connaître son numéro dans "Rastnuméro".

#Fonction pour avoir la moyenne de la precipitation pour chaque pixel pour l'année 2021
for (i in 1:length(df2All$x)) {
  df2All$WCprecimean2021[i] <- (df2All$WCpreci012021[i]+df2All$WCpreci022021[i]+
    df2All$WCpreci032021[i]+df2All$WCpreci042021[i]+df2All$WCpreci052021[i]+df2All$WCpreci062021[i]+
    df2All$WCpreci072021[i]+df2All$WCpreci082021[i]+df2All$WCpreci092021[i]+df2All$WCpreci102021[i]+
    df2All$WCpreci112021[i]+df2All$WCpreci122021[i])/12
}
remove(i)


df2All$WCprecimeanall <- df2All$WCpreci012021
for(i in 1:length(df2All$x)){
  df2All$WCprecimeanall[i] <- (df2All$WCprecimean2017[i]+df2All$WCprecimean2018[i]+
                                 df2All$WCprecimean2019[i]+df2All$WCprecimean2020[i]+
                                 df2All$WCprecimean2021[i])/5
}


#Ajout de la colonne de km2 de forêt par pixel 
for(i in 1:length(df2All$x)){
  df2All$Forestkm2[i]<-(df2All$GFW2010[i]/100)*df2All$maxlnH2016[i]
}

df2All <- transform(df2All, Forestkm2 = as.numeric(Forestkm2))

"Adding ID to each pixel"
for(i in 1:length(df2All$x)){
  df2All$IDpixel[i] <- i
}


#Ajout de la colonne de production de l'élevage :
Valeurfaolostelvage <- 3902886*1000 #Rentrer la valeur donnée par FAostat pour l'élevage en 2017
#cette valeur est en $ interntnal pour le Kenya pour l'item aggregated "Livestock + (Total)
#et pour l'élement "gross production value (constant 2014-2016 thousand I$) sur le site 
#https://www.fao.org/faostat/en/#data/QV
somelevgeKenya <- sum(df2All$grazingH2016,na.rm = T)#Surface total accordée à l'élevage au kenya
#en km2
valelvge <- Valeurfaolostelvage/somelevgeKenya #Valeur moyenne de l'élevage ($) au Kenya par km2
for(i in 1:length(df2All$x)){
  df2All$valelvge2017[i] <- df2All$grazingH2016[i]*valelvge
}
remove(Valeurfaolostelvage,somelevgeKenya)
#La colonne valelvge2017 permet de connaître la valeur totale de l'élevage pour chaque
#pixel, en sachant que la valeur total de l'élevage du kenya est répartie
#de façon homogène. 
for(i in 1:length(df2All$x)){
  df2All$valeurelvgha2017[i] <- df2All$valelvge2017[i]/((df2All$grazingH2016[i])*100)#passage des km² en ha
  ifelse(df2All$valeurelvgha2017[i]=="NaN", df2All$valeurelvgha2017[i]<- 0,"")
}
#La colonne valeurelvgha2017 donne la valeur de production à l'ha pour chaque pixel


## 11 - Vérification de la concordance des données ----
"Vérifier si les espaces occupés par la forêt des données provenant de GFW sont en accords
avec l'utilisation des terres de HYDE"
dfverif <- df2All
dfverif$verif <- ifelse((dfverif$maxlnH2016-dfverif$croplandH2016-dfverif$grazingH2016
                        -dfverif$uoppH2016)<(dfverif$maxlnH2016*(dfverif$GFW2010/100)),1,0)
"Le nombre de pixels où les espaces considérés comme des forêts par GFW,
sont ou en partie, également des espaces de cropland,grazing ou surface batie d'après HYDE :"
print(nrow(subset(dfverif[c("verif")], verif == 1)))

dfverif$verif2 <- ifelse((dfverif$maxlnH2016-dfverif$uoppH2016)<(dfverif$maxlnH2016*
                                                                   (dfverif$GFW2010/100)),1,0)
"Nombre de pixel où les forêts d'après GFW, sont sur des surfaces bâties"
nrow(subset(dfverif[c("verif2")], verif2 == 1))

#Analyse et calcul des données de séquestration carbone 
mean(df2All$regrowthcookpatton, na.rm = T)
valmeanseqcarbcookpatton <- mean(df2All$regrowthcookpatton, na.rm = T)
max(df2All$regrowthcookpatton, na.rm = T)
min(df2All$regrowthcookpatton, na.rm = T)
median(df2All$regrowthcookpatton, na.rm = T)
quantile(subset(df2All[c("regrowthcookpatton")]), probs=seq(0, 1, 0.1), na.rm = T)


## 12 - Nettoyage ----
remove(dfverif)
remove(CRS,emplacementfolder,i,name1,name2,
       numrastassembl1max,numrastassembl1min,numrastassembl2max,numrastassembl2min,
       RESOL, listdf_rast,numcellrast,nam,name_file,files,origin,real_name,
       dfGFW,dfHYDE,dfHydesurface,dfSPAM,dfWC,GFWrast)


#Base de sauvegarde 
df2sauve <- df2All


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# II - Définition des seuils et variables ----
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#La définition des seuils ci-dessous permet de faire varier les résultats. Ces seuils
#sont à définir par l'utilisateur. 
#Attention aux formats des seuils, les seuils de Global Forest Watch donnent des valeurs en %.
#Si les seuils de Global Forest Watch sont modifiés il est nécessaire des les mettre au format % (80 et non 0.8 par exemple).

#Seuils initiaux :
# -- Zone urbaine
seuilbati <-quantile(subset(df2All[c("uoppH2016")],uoppH2016>0), probs=seq(0, 1, 0.05), na.rm = T)[[19]]#Sans les valeurs nulles

# -- Température
seuiltemprturemax <- 35

#Seuils à définir par algorithme : 
# -- Algorithme n°1 
seuilforestalgo1 <- quantile(subset(df2All[c("GFW2010")],GFW2010>0), probs=seq(0, 1, 0.01), na.rm = T)[[99]] #Sans les valeurs nulles
seuil2forestalgo1 <- quantile(subset(df2All[c("GFW2010")],GFW2010>0), probs=seq(0, 1, 0.01), na.rm = T)[[99]] #Sans les valeurs nulles
# -- Algorithme n°2 - partie a)
seuilforestalgo4 <- quantile(subset(df2All[c("GFW2010")],GFW2010>0), probs=seq(0, 1, 0.0001), na.rm = T)[[9958]]
seuil2forestalgo4 <- quantile(subset(df2All[c("GFW2010")],GFW2010>0), probs=seq(0, 1, 0.01), na.rm = T)[[99]]
seuil3forestalgo4 <- quantile(subset(df2All[c("GFW2010")],GFW2010>0), probs=seq(0, 1, 0.01), na.rm = T)[[99]]
seuil4forestalgo4 <- quantile(subset(df2All[c("GFW2010")],GFW2010>0), probs=seq(0, 1, 0.01), na.rm = T)[[99]]
# -- Algorithme n°2 - partie b)
seuilforestalgo42 <- quantile(subset(df2All[c("GFW2010")],GFW2010>0), probs=seq(0, 1, 0.01), na.rm = T)[[99]]
seuilforest41 <- quantile(subset(df2All[c("GFW2010")],GFW2010>0), probs=seq(0, 1, 0.01), na.rm = T)[[99]]
seuilforest42 <- quantile(subset(df2All[c("GFW2010")],GFW2010>0), probs=seq(0, 1, 0.01), na.rm = T)[[99]]

# -- Précipitation 
dfprecipthreshold <- subset(df2All[c("WCprecimeanall","GFW2010")], GFW2010 >= seuilforestalgo1)
dfprecipthreshold <- dfprecipthreshold[order(dfprecipthreshold$WCprecimeanall, decreasing = F),]
seuilpreci <- dfprecipthreshold$WCprecimeanall[1]
remove(dfprecipthreshold)

#Seuils spécifiques externes :
# -- Pour le calcul des rendements (en valeurs $) voisin d'une forêt
seuilforestgain1 <- quantile(subset(df2All[c("GFW2010")],GFW2010>0), probs=seq(0, 1, 0.01), na.rm = T)[[99]]


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# III - Ajustement des données relatives aux seuils ----
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

"Ajout de la colonne donnant une indication pour la température max mensuelle en 2021 en °C "
df2All$Tempmaxcond <- ifelse(df2All$WCtmax012021<seuiltemprturemax, 
                             ifelse(df2All$WCtmax022021<seuiltemprturemax,
                                    ifelse(df2All$WCtmax032021<seuiltemprturemax,
                                           ifelse(df2All$WCtmax042021<seuiltemprturemax,
                                                  ifelse(df2All$WCtmax052021<seuiltemprturemax,
                                                         ifelse(df2All$WCtmax062021<seuiltemprturemax,
                                                                ifelse(df2All$WCtmax072021<seuiltemprturemax,
                                                                       ifelse(df2All$WCtmax082021<seuiltemprturemax,
                                                                              ifelse(df2All$WCtmax092021<seuiltemprturemax,
                                                                                     ifelse(df2All$WCtmax102021<seuiltemprturemax,
                                                                                            ifelse(df2All$WCtmax112021<seuiltemprturemax,
                                                                                                   ifelse(df2All$WCtmax122021<seuiltemprturemax,1,0),0),0),0)
                                                                              ,0),0),0),0),0),0),0),0)
"Ajout de la colonne donnant une indication sur la température mensuelle moyenne maximale par pixel en 2021 en °C"
df2All$WCtmaxmean2021 <- (df2All$WCtmax012021+df2All$WCtmax022021+df2All$WCtmax032021+
                         df2All$WCtmax042021+df2All$WCtmax052021+df2All$WCtmax062021+
                         df2All$WCtmax072021+df2All$WCtmax082021+df2All$WCtmax092021+
                         df2All$WCtmax102021+df2All$WCtmax112021+df2All$WCtmax122021)/12


"Ajout de la colonne donnant une valeur 1 si le pixel respecte les conditions définies en termes
de surfaces urbanisées dans le pixel(km²), température mensuelle maximale en 2021(°C), précipitation 
mensuelle moyenne pour la période 2017-2021 (mm)"
df2All$Condtnbatprecitemp <- ifelse(df2All$uoppH2016<seuilbati,
                                    ifelse(df2All$Tempmaxcond == 1,
                                           ifelse(df2All$WCprecimeanall>=seuilpreci,1,0),0),0)

"Ajout de la colonne de production totale SPAM en quantité (en tonne metric) :"
df2All$SPAMProdquant2010_mtunit <- df2All$SPAMACOFP+df2All$SPAMBANAP+df2All$SPAMBARLP+
  +df2All$SPAMBEANP+df2All$SPAMCASSP+df2All$SPAMCHICP+df2All$SPAMCNUTP+
  df2All$SPAMCOCOP+df2All$SPAMCOTTP+df2All$SPAMCOWPP+df2All$SPAMGROUP+
  df2All$SPAMLENTP+df2All$SPAMMAIZP+df2All$SPAMOCERP+df2All$SPAMOFIBP+
  df2All$SPAMOILPP+df2All$SPAMOOILP+df2All$SPAMOPULP+df2All$SPAMORTSP+
  df2All$SPAMPIGEP+df2All$SPAMPLNTP+df2All$SPAMPMILP+df2All$SPAMPOTAP+
  df2All$SPAMRAPEP+df2All$SPAMRCOFP+df2All$SPAMRESTP+df2All$SPAMRICEP+
  df2All$SPAMSESAP+df2All$SPAMSMILP+df2All$SPAMSORGP+df2All$SPAMSOYBP+
  df2All$SPAMSUGBP+df2All$SPAMSUGCP+df2All$SPAMSUNFP+df2All$SPAMSWPOP+
  df2All$SPAMTEASP+df2All$SPAMTEMFP+df2All$SPAMTOBAP+df2All$SPAMTROFP+
  df2All$SPAMVEGEP+df2All$SPAMWHEAP+df2All$SPAMYAMSP



"Création d'une nouvelle colonne où les NA de SPAM sont remplacés par la valeur de
production de valeurelvgha2017"
for(i in 1:length(df2All$x)){
  df2All$SPAM2[i] <- ifelse(is.na(df2All$SPAMProdvalha2010[i]) == "FALSE",
                            df2All$SPAMProdvalha2010[i], df2All$valeurelvgha2017[i])
}

for(i in 1:length(df2All$x)){
  df2All$SPAM2[i] <- ifelse(is.na(df2All$valeurelvgha2017[i]) == "FALSE",
                            ifelse(df2All$valeurelvgha2017[i] > df2All$SPAM2[i],
                            df2All$valeurelvgha2017[i],df2All$SPAM2[i]),df2All$SPAM2[i])
}

#La colonne SPAM2 possède les valeurs de production de l'élevage par ha (basées sur 
#valeurelvgha2017) lorsqu'il n'y a pas de valeur dans SPAM. 
#De plus, on modifie SPAM2 pour que la colonne SPAM2, prenne aussi les valeurs maximales de production
#de l'agriculture. Donc s'il y a une valeur pour SPAM est qu'elle est inférieure à une valeur de production
#de l'élevage, on prend la valeur la plus élevée donc, ici, la valeur d'élevage. Car la colonne SPAM2 dans les algorithmes
#vise à mesurer la perte potentielle de production agricole maximale, en d'autres termes, Quelle pourrait
#être la meilleure utilisation du pixel si on gardait l'espace et quelle valeur on pourrait en tirer au max ?

#(((De ce fait, on considère que les NA SPAM s'apparentent à des valeurs 0 en production agricole de culture))).
#Mais on n'ajoute pas les valeurs de valeurelvgha2017 à chaque fois dans SPAM car SPAM possède des coûts
#d'opportunité intégrés dans ces valeurs.




# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# IV Analyse des données ----
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# -- Analyse générale des données :
"Nombre de pixel correspondant aux critères de températue :"
nrow(subset(df2All[c("Tempmaxcond")], Tempmaxcond == 1))
"Nombre de pixel repondant critère précipitation : "
df2All$preci <- ifelse(df2All$WCprecimeanall>seuilpreci,1,0)
nrow(subset(df2All[c("preci")], preci == 1))
df2All <- df2All[ , !(names(df2All) %in% c("preci"))]
"Nombre de pixel répondant au critère de préci, temp et surface bâtie :"
nrow(subset(df2All[c("Condtnbatprecitemp")], Condtnbatprecitemp == 1))
"Surface totale de forêt au kenya actuelle en km2:"
sum(df2All$Forestkm2, na.rm = T)

"Number of pixel with NA for SPAM :"
sum(is.na(df2All$SPAMProdvalha2010))
"Number of pixel with NA for Global Forest Watch :"
sum(is.na(df2All$GFW2010))
"Number of pixel with NA for SPAM2 :"
sum(is.na(df2All$SPAM2))
"Number of pixel with NA for HYDE - cropland :"
sum(is.na(df2All$croplandH2016))
"Number of pixel with NA for HYDE - land area of the pixel :"
sum(is.na(df2All$maxlnH2016))
"Number of pixel with NA for HYDE - grazing :"
sum(is.na(df2All$grazingH2016))
"Number of pixel with 0 in value for HYDE - cropland :"
nrow(subset(df2All[c("croplandH2016")], croplandH2016 == 0))
"Number of pixel with 0 in value for SPAM2 :"
nrow(subset(df2All[c("SPAM2")], SPAM2 == 0))
"Production moyenne par ha (en tonne/ha) :"
sum(df2All$SPAMProdquant2010_mtunit, na.rm = T)/(sum(df2All$maxlnH2016, na.rm = T)*100)



## I - Algorithme n°1 ----
"Algorithme de minimisation des pertes de production agricole potentielle"

#Algorithme fonctionnant en minimisant la perte de production agricole potentielle (étant en $/ha).
#Un seuil de production agricole ($/ha) va être défini au début de l'algorithme par seuilSPAMvalalgo1.
#Les pixels étant inférieurs à ce seuil et respectant les conditions vont être reforestés. 


#Rappel : les données agricoles utilisées pour le calcul des pertes de production viennent
#de la colonne SPAM2 qui prend en compte les valeurs d'élevage et les données de SPAM. Ainsi,
#si un pixel a un NA SPAM et 0 espace consacré à l'élevage on considère qu'il a 0 en valeur de production
#agricole.


#Explication de l'algorithme n°1 : 
#L'utilisateur rentre un objectif de % tree cover au Kenya à atteindre dans "objectifalgo1".
#Le % de tree cover représente le % de couverture arborée sur l'ensemble du territoire Kenyan.
#Il rentre également un niveau de millentile (percentile X 10) permettant de définir le 
#niveau de production agricole en valeurs (en $ dans SPAM) en dessous duquel on sélectionne les pixels
#pour les reforester. L'utilisateur rentre ce niveau dans "number". En sachant qu'à chaque
#boucle de l'algorithme ce niveau va être augmenté jusqu'à atteindre le niveau souhaité de
#tree cover. En l'occurence, ce niveau va être augmenté par "echelajoutalgo1". Donc si par
#example, l'utilisateur définit "number" = 500 et "echelajoutalgo1" = 20, à chaque boucle "number"
#va augmenter de 20, permettant de sélectionner des pixels avec des productions agricoles SPAM
#en valeurs ($) plus élevées et donc de sélectionner plus de pixels afin d'atteindre une tree cover
#plus élevée.

#/!\
#Il est nécessaire de lancer la partie "Nettoyage" avant de relancer la boucle
#de l'algorithme avec d'autres conditions.




### A) Calcul de la Tree cover ----

objectifalgo1 <- 7
compteuralgo1 <- 1
number <- 800
echelajoutalgo1 <-100
kad <- 0 #si l'on souhaite atteindre le maximum de couverture arborée possible il faut que kad = 1 
#si kad = 1 alors on prend même le pixel étant le plus forestés.

#La boucle ci-dessous permet d'obtenir l'objectifalgo1 en élevant à chaque tour le seuilSPAMvalalgo1
while(compteuralgo1 < objectifalgo1){
  number <- number+echelajoutalgo1
  seuilSPAMvalalgo1 <- quantile(df2All[c("SPAM2")], probs=seq(0, 1, 0.001), na.rm = T)[[number]]
  
  print("Les résultats ci-dessous sont calculés en fonction des conditions définies haut-dessus par l'utilisateur.")
  for(i in 1:length(df2All$x)){
    df2All$Plantalgo1[i] <- ifelse(df2All$Condtnbatprecitemp[i] == 1,
                                   ifelse(is.na(df2All$SPAM2[i]) == "FALSE",
                                          ifelse(df2All$SPAM2[i]<(seuilSPAMvalalgo1+kad),
                                                 ifelse(df2All$GFW2010[i]<(seuilforestalgo1+kad),1,2),3),4),0)
  }#Obtention des pixels respectant les conditions de précipitation, température, urbanisation,
  #production agricole potentielle et ayant une couverture arborée inférieure à seuilforestalgo1.
  
  print("Nombre de pixel répondant correctement aux conditions (aire urbanisée,précipitation
  température,valeur de production agricole potentielle et surface arborée) :")
  print(nrow(subset(df2All[c("Plantalgo1")], Plantalgo1 == 1)))

  #Dans la colonne Plantalgo1, pixels = 1 = pixel à reforester car répondant à l'ensemble des 
  #conditions, pixel = 2 = pixel ayant déjà des surfaces de forêt élevés, pixel = 3 = pixel
  #répondant condtions bat preci temp et ayant valeur pour SPAM mais valeur SPAM trop élevée,
  #pixel = 4 = pixel ok conditions bat preci temp mais NA pour SPAM (très peu avec SPAM2) 
  #et pixel = 0 = pixel répondant pas condition bat preci temp
  
  #### 1) Calcul des surfaces de forêt ----
  for(i in 1:length(df2All$x)){
    df2All$km2plantalgo1[i] <- ifelse(df2All$Plantalgo1[i] == 1,
                                      ifelse(df2All$GFW2010[i]<seuil2forestalgo1,
                                             (df2All$maxlnH2016[i]-df2All$uoppH2016[i])*(seuil2forestalgo1/100),
                                             df2All$Forestkm2[i]),
                                      ifelse(df2All$Plantalgo1[i] == 2,
                                             df2All$Forestkm2[i],
                                             ifelse(df2All$Plantalgo1[i] == 3,
                                                    df2All$Forestkm2[i],
                                                    ifelse(df2All$Plantalgo1[i] == 4,
                                                           df2All$Forestkm2[i],df2All$Forestkm2[i]))))
  }#Augmentation de la couverture arborée des pixels ayant plantalgo1 = 1 au seuil2forestalgo1.
  #Les autres pixels gardent leur surface initiale de forêt. 
  
  for(i in 1:length(df2All$x)){
    df2All$km2plantalgo1[i] <- ifelse(df2All$Plantalgo1[i] == 0,
                                      df2All$Forestkm2[i],df2All$km2plantalgo1[i])
  }#Ajout dans la colonne km2plantalgo1 des surfaces de forêts initiales pour les pixels ne répondant pas
  #correctement aux conditions fixés.
  for(i in 1:length(df2All$x)){
    df2All$km2plantalgo1[i] <- ifelse(is.na(df2All$km2plantalgo1[i]) == "FALSE",
                                     ifelse(df2All$km2plantalgo1[i]>df2All$Forestkm2[i],
                                            df2All$km2plantalgo1[i],df2All$Forestkm2[i]),
                                     df2All$Forestkm2[i])
  }
  
  #La colonne km2plantalgo1 donne pour chaque pixel, la surface de forêt après reforestation
  #au seuil2forestalgo1, en km²
  
  
  #Obtention du nouveau pourcentage de forêt du pixel 
  for(i in 1:length(df2All$x)){
    df2All$GFW2010_new_1[i] <- (df2All$km2plantalgo1[i]/df2All$maxlnH2016[i])*100}
  print("% de couverture arborée atteint :")
  print(sum(df2All$km2plantalgo1, na.rm = T)/sum(df2All$maxlnH2016, na.rm = T)*100)
  
  #Surface totale de forêt en km² qu'on peut planter (selon les critères haut dessus)
  for (i in 1:length(df2All$x)){
    df2All$diffalgo1[i] <- ifelse(df2All$km2plantalgo1[i]>df2All$Forestkm2[i],
                                  (df2All$km2plantalgo1[i]-df2All$Forestkm2[i]),0)
  }
  #diffalgo1 donne pour chaque pixel répondant aux conditions, la surface de
  #forêt (km²) que l'on peut planter avec la reforestation (pour les pixels ayant valeur SPAM2)
  print("Nombre de km² de forêt plantés :")
  print(sum(df2All$diffalgo1, na.rm = T))
  
  ### 2) Calcul de la perte de production ----
  #Perte de production agricole SPAM2
  for (i in 1:length(df2All$x)){
    df2All$loosevalalgo1[i] <- ifelse(df2All$diffalgo1[i] > 0,
                                      df2All$SPAM2[i]*(df2All$diffalgo1[i]*100),0)
  } 
  #La colonne loosevalalgo1 donne, pour chaque pixel, la perte totale de production agricole
  #(incluant production des champs et élevage) potentielle en $ sur le pixel engendrée par la
  #reforestation. 
  
  print("La perte de production agricole potentielle totale pour l'ensemble du Kenya pour une année ($):")
  print(sum(df2All$loosevalalgo1, na.rm = T))
  print("Perte de production agricole potentielle moyenne par ha ($/an):")#SPAM2 avec NA non pris en compte.
  print(sum(df2All$loosevalalgo1, na.rm = T)/(sum(df2All$diffalgo1, na.rm = T)*100))
  
  #Calcul de la perte de Quantitée produite (en tonne) et en se basant 
  #sur les données SPAM : 
  for (i in 1:length(df2All$x)){
    df2All$loosequantialgo1[i]<- ifelse(df2All$diffalgo1[i]>0,
                                        df2All$SPAMProdquant2010_mtunit[i]*(df2All$diffalgo1[i]/df2All$maxlnH2016[i]),0)
  }
  #loosequantialgo1 donne pour chque pixel la perte de production agricole 
  #en termes de qttée (en tonne métrique) engendré par l'espace pris par la reforestation
  print("La perte de production agricole potentielle totale pour l'ensemble du Kenya pour une année
  (en tonne et seulement basée sur les productions données par SPAM) :")
  print(sum(df2All$loosequantialgo1, na.rm = T))
  print("La perte de production agricole potentielle moyenne par hectare pour
        une année est (en tonne et basée juste sur données SPAM) :")
  print(sum(df2All$loosequantialgo1, na.rm = T)/(sum(df2All$diffalgo1, na.rm = T)*100))
  
  print("Indication sur le percentile :")
  print((number-1)/10)
  print("Indication sur la valeur SPAM2 affiliée à ce percentile (en $):")
  print(quantile(df2All[c("SPAM2")], probs=seq(0, 1, 0.001), na.rm = T)[[number]])
  
  #Compteur pour la boucle : 
  compteuralgo1 <- sum(df2All$km2plantalgo1, na.rm = T)/sum(df2All$maxlnH2016, na.rm = T)*100
  
  
  print("-_-_-_-_-_ NEXT-_-_-_-_-_-")
  print("-_-_-_-_-_-_-_-_-_-_-_-_-_")
  if(number>=1001){stop()}
}


#Nettoyage
df2All <- df2All[ , !(names(df2All) %in% c("Plantalgo1",
                                           "km2plantalgo12","loosevalalgo12",
                                           "loosequantialgo12","km2forplantalg1"))]


### B) Calcul d'une Tree cover spécifique ----

#La boucle ci-dessous permet d'obtenir un niveau de tree cover spécifique. Par exemple,
#si en lançant l'algorithme n°1 haut-dessus la tree cover obtenue est 18,46% et que l'utilisateur
#souhaite obtenir 18% il peut lancer la boucle ci-dessous qui va rétablir les surfaces de forêts
#initiales sur certains pixels.

#Explication du fonctionnement du code : 
#(La boucle ci-dessous, trie les pixels reforestés en fonction de leur perte agricole.
#Les pixels avec les plus grosses pertes agricoles sont enlevés des pixels reforestés jusqu'à
#l'atteinte du niveau spécifique de tree cover souhaité.)
#L'utilisateur doit rentrer l'objectif à atteindre (% tree cover au Kenya) dans compteurpromp1
#L'utilisateur doit rentrer combien de pixels sont sélectionnés pour être retirés
#des pixels reforestés dans echelleretireralgo1.
#Après la boucle while, il est nécessaire de lancer le reste des codes en dessous.
df2Allalg1 <- df2All
compteurpromp1 <- 20.06 #doit être un tout petit peu haut dessus de l'objectif
compteurpromp2 <- 25 #doit être > à compteurpromp1
while(compteurpromp1<compteurpromp2){
  df2Allalg1 <- df2Allalg1 %>% arrange(desc(loosevalalgo1))
  echelleretireralgo1 <- 10#Nombre de pixels qu'on remet aux niveaux de forêt initial
  #à chaque tour de la boucle.
  for(i in 1:echelleretireralgo1){
    df2Allalg1$km2plantalgo1[i]<-df2Allalg1$Forestkm2[i]
    df2Allalg1$diffalgo1[i]<-0
    df2Allalg1$loosequantialgo1[i]<-0
    df2Allalg1$loosevalalgo1[i]<-0
  }
  compteurpromp2<-sum(df2Allalg1$km2plantalgo1,na.rm = T)/sum(df2Allalg1$maxlnH2016,na.rm=T)*100
  print("Le niveau de tree cover atteint est")
  print(compteurpromp2)
}

#De plus, dans certains cas où on vise une couverture arborée faible, il est possible
#qu'on ne puisse pas la diminuer. Cet algorithme se base sur la perte de production
#agricole pour diminuer les surfaces forestées, sauf que certains pixels peuvent
#être reforestés sans avoir de perte de production agricole, car certains pixels ont
#une SPAM2 = 0 soit provenant d'un NA SPAM soit provenant d'une surface de 0 en grazing.
#Dans ce cas, par exemple pour 7% de tree cover, on a une perte totale de production agricole
#potentielle valant 0 pour l'ensemble du Kenya. 

#/!\

#Pour ajuster les données de diffalgo (étant donné qu'on retire des pixels reforestés
#pour atteindre un certain seuil il faut relancer diffalgo pour avoir les surfaces de forêts reforestées mises à jour) : 
for(i in 1:length(df2Allalg1$x)){
  df2Allalg1$diffalgo1[i] <- ifelse(df2Allalg1$km2plantalgo1[i]>df2Allalg1$Forestkm2[i],
                                   (df2Allalg1$km2plantalgo1[i]-df2Allalg1$Forestkm2[i]),0)
}
#Obtention du nouveau pourcentage de forêt du pixel 
for(i in 1:length(df2Allalg1$x)){
  df2Allalg1$GFW2010_new_1[i] <- (df2Allalg1$km2plantalgo1[i]/df2Allalg1$maxlnH2016[i])*100}

"Avec ce niveau il faut planter combien de km2 de forêt :"
print(sum(df2Allalg1$diffalgo1, na.rm = T))
print("Pourcentage de tree cover :")
((sum(df2Allalg1$Forestkm2,na.rm=T)+sum(df2Allalg1$diffalgo1,na.rm=T))/sum(df2Allalg1$maxlnH2016,na.rm=T))*100
print("Avec ce niveau la perte de production totale agricole en valeur ($) en plantant
des forêts selon les conditions définies s'élève pour une année à :")
print(sum(df2Allalg1$loosevalalgo1, na.rm = T))
print("Cette perte de production agricole potentielle (en $) s'élève en moyenne pour une année pour 1 hectare à :")
print(sum(df2Allalg1$loosevalalgo1, na.rm = T)/(sum(df2Allalg1$diffalgo1, na.rm = T)*100))
print("La perte de production agricole potentielle totale pour l'ensemble du Kenya pour une année
  (en tonne et seulement basée sur les productions données par SPAM) :")
print(sum(df2Allalg1$loosequantialgo1, na.rm = T))
print("La perte de production agricole potentielle moyenne par hectare pour
        une année est (en tonne et basée juste sur données SPAM) :")
print(sum(df2Allalg1$loosequantialgo1, na.rm = T)/(sum(df2Allalg1$diffalgo1, na.rm = T)*100))
#Pour avoir la différence de pourcentage de forêt d'un pixel entre avant reforestation
#et après :
for(i in 1:length(df2Allalg1$x)){
  df2Allalg1$diffreforesta1[i] <- ifelse(df2Allalg1$GFW2010_new_1[i]>df2Allalg1$GFW2010[i],
                                         (df2Allalg1$GFW2010_new_1[i]-df2Allalg1$GFW2010[i]),0)
}
df2All$diffalgo1 <- df2Allalg1$diffalgo1
df2All$loosevalalgo1 <- df2Allalg1$loosevalalgo1


### C) Réalisation des cartes & graphs ----
"Attention à supprimer les fichiers shp si l'on souhaite relancer les codes ci-dessous.
De plus, attention à bien les nommer dans la partie 'layer = ....'  "

# -- Carte représentant les espaces reforestés (km²) :
df1algo1 <- subset(df2Allalg1[c("x","y","diffalgo1")])
dshp <- df1algo1
coordinates(dshp)=~x+y
crs(dshp) <- crsref
dshp <- st_as_sf(dshp)
#Sortie carte
st_write(dshp, dsn = "diffalgo1.shp", layer = "diffalgo1", driver = "ESRI Shapefile")
remove(df1algo1,dshp)
#ou en raster :
df1algo1dif  <- subset(df2All[c("x","y","diffalgo1")])
rastdiffalgo1 <- rasterFromXYZ(df1algo1dif)
writeRaster(rastdiffalgo1, "rastalgo110%", format = "GTiff")
remove(rastdiffalgo1,df1algo1dif)

# -- Carte représentant les pertes de production agricole potentielle en $  :
df1algo1 <- subset(df2Allalg1[c("x","y","loosevalalgo1")])
dshp <- df1algo1
coordinates(dshp)=~x+y
proj4string(dshp) <- proj4string(Rastdecoupe)
#Sortie carte
writeOGR(dshp, dsn = 'coucheshpR', layer = "loosevalalgo1", driver = "ESRI Shapefile")
remove(df1algo1,dshp)

# -- Graph représentant les pertes de production agricole potentielle par counties
"Get the df with the conditions and transform it into shp"
dfcrossdatacounty1 <- subset(df2All[c("x","y","IDpixel","loosevalalgo1")])
coordinates(dfcrossdatacounty1) <- ~ x + y
shapecunties <- sf:::as_Spatial(shapecunties)
proj4string(dfcrossdatacounty1) <- proj4string(shapecunties)
dfcrossdatacounty1  <- st_as_sf(dfcrossdatacounty1)
shapecunties <- st_as_sf(shapecunties)
"Join the shapefil to get for each type of soil if there are pixel answering the conditions"
dfcheckcounties1 <- st_join(dfcrossdatacounty1, shapecunties,left = TRUE, largest = TRUE)
"Create the graphique"
#Faire la somme des loosevalalgo1 par county : 
uniquecountykenya <- unique(dfcheckcounties1$county)
dfsumcounty <- data.frame(c(1,1,1))
colnames(dfsumcounty)[colnames(dfsumcounty) == 'c.1..1..1.'] <- 'Sumagriloos'
dfcountyname <- data.frame(c(1,1,1))
colnames(dfcountyname)[colnames(dfcountyname) == 'c.1..1..1.'] <- 'Countyname'
for(i in 1:length(uniquecountykenya)){
  name<-uniquecountykenya[i]
  dfcountyname <- rbind(dfcountyname,name)
  sum1<-sum(subset(dfcheckcounties1[c("loosevalalgo1","county")], county == name)$loosevalalgo1)
  dfsumcounty <- rbind(dfsumcounty,sum1)
}
dfsumcounties<-data.frame(dfcountyname,dfsumcounty)
dfsumcounties <- dfsumcounties[-c(1,2,3),]
#erreur :
sum(df2All$loosevalalgo1,na.rm=T)-sum(dfsumcounties$Sumagriloos,na.rm=T)
sum(dfcheckcounties1$loosevalalgo1,na.rm=T)
#l'erreur vient du fait que certains pixels durant la jointure avec le shapefile
#des counties du Kenya ne sont pas référencés dans un county, et donc ont une NA
#dans la colonne county et donc les pertes de production agricole potentielle de ces
#pixels ne sont pas attribués à leur county.

#dfsumcounties donne, pour un niveau de tree cover (%) précisé
#par l'utilisateur et après que l'utilisateur ait lancé l'algorithme
# le niveau de perte totale de production agricole en $ par county du Kénya.
remove(dfsumcounty,dfcountyname,dfsumcounties,
       uniquecountykenya,dfcheckcounties1,shapecunties,
       dfcrossdatacounty1)



### D) Nettoyage ----
remove(df2Allalg1)
df2All <- df2All[ , !(names(df2All) %in% c("loosevalalgo1",
                                           "loosequantialgo1","GFW2010_new_1",
                                           "diffreforesta1","km2plantalgo1"))]
df2All <- df2All[ , !(names(df2All) %in% c("Plantalgo1","km2plantalgo12",
                                           "km2forplantalg1","loosevalalgo12",
                                           "loosequantialgo12","pertemoyparha"))]
#Si l'utilisateur n'utilise pas le code R_agroforestry après, il peut lancer
#le nettoyage de la ligne juste ci-dessous :
df2All <- df2All[ , !(names(df2All) %in% c("diffalgo1"))]




## II - Algorithme n°2 - partie A ----
"Algorithme de continuité écologique"

#Algorithme fonctionnant en visant la continuité écologique. Cet algorithme va chercher à
#reforester les pixels étant adjacents aux forêts existantes. Il va dans un premier temps sélectionner
#des pixels répondant correctement aux conditions (précipitation, surface urbanisée, température) et ayant
#un certain niveau de couverture arborée. Ces pixels vont donc être considérés comme des forêts.
#Ensuite il va regarder les pixels adjacents (méthode = "Investigating Commuting Time in a Metropolitan
#Statistical Area Using Spatial Autocorrelation Analysis" -> Queen's case - Hessam Miri) à ces pixels-forêt et
#regarder si ces premiers pixels adjacents respectent les conditions (précip, urbain, température).
#Ensuite l'algorithme va regarder les pixels adjacents des premiers pixels adjacents respectant les conditions.
#Dans ces "deuxième pixels adjacents", l'algorithme va regarder les pixels qui respectent les conditions (précip,urbain,temp).
#Enfin, l'algorithme augmente la couverture arborée de ces 3 types de pixels : 1) pixels étant considérés comme des forêts
#2) pixels étant les premiers adjacents et respectant les conditions   3) Deuxième pixels adjacents et respectant conditions.

#Explication de l'algorithme n°2 - partie a) : 
#L'utilisateur rentre un objectif de % tree cover au Kenya à atteindre dans "objalgo4part1".
#Le % de tree cover représente le % de couverture arborée sur l'ensemble du territoire Kenyan.
#L'utilisateur rentre également un niveau de dix millentile (percentile X 100) de couverture 
#arborée (% Tree cover) à partir duquel les premiers pixels
#étant considérés comme des forêts vont être sélectionnés. L'utilisateur rentre ce niveau dans "numberalgo4part1".
#En sachant qu'à chaque boucle de l'algorithme, ce niveau va être diminué jusqu'à l'atteinte
#de la tree cover souhaité à l'échelle du Kenya. En l'occurence, ce niveau va être diminué par
# "stepalgo4part1". Donc, par exemple, si l'utilisateur définit "numberalgo4part1" = 500 et
# "stepalgo4part1" = 20, pour atteindre l'objalgo4part1 le numberalgo4part1 diminuera de 20 à 
#chaque tour de l'algorithme. Cette diminution de numberalgo4part1 permettra de sélectionner
# plus de pixels 'car le niveau requis de couverture arborée sera plus faible) et donc
#d'avoir plus de pixels adjacents et donc d'atteindre une tree cover pour le Kenya plus élevée.


#/!\
#Il est nécessaire de lancer la partie "Nettoyage" avant de relancer la boucle
#de l'algorithme avec d'autres conditions.

### A) Calcul de la Tree cover ----

objalgo4part1 <- 1
compteuralgo4part1 <- 1 #Doit être < à objalgo4part1
stepalgo4part1 <- 40
numberalgo4part1 <- 5163
while(compteuralgo4part1 < objalgo4part1 ){
  print("Les résultats ci-dessous sont calculés en fonction des conditions définies haut-dessus par l'utilisateur.")
  numberalgo4part1<- numberalgo4part1 - stepalgo4part1
  seuilforestalgo4 <- quantile(subset(df2All[c("GFW2010")],GFW2010>0), probs=seq(0, 1, 0.0001), na.rm = T)[[numberalgo4part1]]
  quantile(subset(df2All[c("GFW2010")],GFW2010>0), probs=seq(0, 1, 0.0001), na.rm = T)[[numberalgo4part1]]
  print("Percentile de GFW utilisé (distribution sans valeur nulle) :")
  print((numberalgo4part1-1)/100)
  print("Valeur de tree cover (%) affilié à ce percentile :")
  print(quantile(subset(df2All[c("GFW2010")],GFW2010>0), probs=seq(0, 1, 0.0001), na.rm = T)[[numberalgo4part1]])
  
  #### 1) Calcul des surfaces de forêt ----
  "Avoir les pixels ayant une couverture arborée supérieure au seuilforestalgo4 et respectant conditions:"
  df2All$Plantindicealgo4 <- ifelse(df2All$Condtnbatprecitemp == 1,
                                    ifelse(df2All$GFW2010 > seuilforestalgo4,1,0),0)
  
  
  #Obtention des premiers pixels adjacents par rapport aux cellules ayant 1 dans Plantindicealgo4
  # -- a) Convertir en raster
  dfalgo4 <- df2All[,c("x","y","Plantindicealgo4")]
  Rastercomparalgo4 <-rasterFromXYZ(dfalgo4)#Raster de comparaison
  Rasterdfalgo4 <- rasterFromXYZ(dfalgo4)
  numcellrast <- nrow(dfalgo4)
  # -- b) Création du df
  from <- c("1")
  to <- c("1")
  Dfadjacent <- data.frame(from,to)
  # -- c) Chaque cellule du raster est vérifiée, si elle correspond a Plantindicealgo4 ==1
  #alors on sélectionne avec adjacent() les cellules étant adjacentes à celle-ci.
  #Ensuite on intègre cela dans un df
  for(i in 1:numcellrast){
    ifelse(is.na(extract(Rasterdfalgo4,c(1,i))[2]) == "FALSE",
           ifelse(extract(Rasterdfalgo4,c(1,i))[2] == 1,
                  Dfadjacent<-rbind(Dfadjacent,adjacent(Rasterdfalgo4, cells=c(i),directions = 8, pairs = T)),0),0)
  }
  Dfadjacent <- Dfadjacent[-(1),] #suppression de la première ligne du df créée avec from = c(1)
  remove(from,to,i)
  
  # -- d) Supprimer les doublons et mise en forme
  Dfadjacent<-Dfadjacent[!duplicated(Dfadjacent[,2]),]
  Dfadjacent <- transform(Dfadjacent, to = as.numeric(to))
  
  # -- e) Modification des valeurs du raster par rapport au numéro de cellule
  "Code pour compter le nombre de NA dans les cellules adjacentes :"
  for(i in 1:length(Dfadjacent$to)){
    n <- Dfadjacent$to[i]
    Dfadjacent$count[i] <- ifelse(is.na(Rasterdfalgo4[n]) == "TRUE", 1, 0)
    remove(n)
  }
  "Nombre de NA :"
  nrow(subset(Dfadjacent[c("count")], count == 1))
  
  "Cette fonction permet de convertir les numéros des cellules adjacentes récupérées
  (dans le raster) et de les intégrer dans un raster pour avoir les coordonnées x 
  et y de ces cellules comme c'est le cas dans df2All (en sachant que les
  cellules adjacentes peuvent seulement prendre la valeur 2 si pas de NA
  sur le pixel initialement et 3 s'il y avait une NA sur le pixel initialement)."
  for(i in 1:length(Dfadjacent$to)){
    n <- Dfadjacent$to[i]
    ifelse(is.na(Rasterdfalgo4[n]) == "FALSE",
           ifelse(Rasterdfalgo4[n] == 1,"",
                  Rasterdfalgo4[n] <- 2),
           Rasterdfalgo4[n] <- 3)
    remove(n)
  }
  remove(i)
  
  # -- f) Transformation du raster en dataframe
  dfalgo4rast <- as.data.frame(rasterToPoints(Rasterdfalgo4,fun = NULL))
  
  # -- g) Jointure
  dfalgo4rast <- dfalgo4rast %>% rename(celladjalgo4 = Plantindicealgo4)
  for (i in 1:length(dfalgo4rast$x)){
    dfalgo4rast$x[i]<- round(x = dfalgo4rast$x[i],10)
    dfalgo4rast$y[i]<- round(x = dfalgo4rast$y[i],10)
  }#Permet d'arrondir les valeurs pour permettre la jointure
  df2All <- full_join(df2All, dfalgo4rast, by = c("x","y"))
  #On enlève les pixels qui se sont rajoutés :
  df2All <- df2All[is.na(df2All$GFW2010) == "FALSE", ] 
  
  "celladjalgo4 permet de donner les cellules adjacentes avec 2 et 3 en valeurs
  dans la colonne."
  "Nombre de première cellule adjacente :"
  nrow(subset(df2All[c("celladjalgo4")], celladjalgo4 > 1))
  
  "Vérification des conditions pour ces cellules adjacentes :"
  df2All$Plantalgo4 <- ifelse(df2All$celladjalgo4 >1,
                              ifelse(df2All$Condtnbatprecitemp == 1,
                                     1,0),0)
  "Nombre de pixel répondant à cette condition : "
  nrow(subset(df2All[c("Plantalgo4")], Plantalgo4 == 1))
  ntry <- nrow(subset(df2All[c("Plantalgo4")], Plantalgo4 == 1))
  #Création d'une nouvelle colonne pour que les cellules adjacentes des cellules adjacentes
  #ne puissent pas contenir les pixels de base que l'on a utilisé pour déterminer les
  #premières cellules adjacentes.
  df2All$Plantcondialgo4 <- ifelse(df2All$Plantalgo4 == 1, 1, 0)
  df2All$Plantcondialgo4 <- ifelse(df2All$Plantcondialgo4 == 0,
                                   ifelse(df2All$Plantindicealgo4 == 1, 1,
                                          df2All$Plantcondialgo4),
                                   df2All$Plantcondialgo4)
  "Plantalgo4 donne la liste  (si = 1) des pixels étant des cellules adjacentes et répondant
  aux conditions de précipitation, température et surface bâtie."
  "Plantcondialgo4 permet d'éviter le dédoublement des cellules adjacentes pour la suite."
  remove(Dfadjacent, Rastercomparalgo4, Rasterdfalgo4,dfalgo4,dfalgo4rast)
  
  #Vérification, problème si ce n'est pas = 0:
  print("Vérification devant être égale à 0 : ")
  print(nrow(subset(df2All[c("Plantindicealgo4")], Plantindicealgo4 == 1)) - (nrow(subset(df2All[c("Plantcondialgo4")], Plantcondialgo4 == 1))-ntry))
  remove(ntry)
  
  #Obtention des deuxièmes pixels adjacents : 
  dfalgo4bis <- df2All[,c("x","y","Plantcondialgo4")]
  Rastercomparalgo4bis <-rasterFromXYZ(dfalgo4bis)#Raster de comparaison
  Rasterdfalgo4bis <- rasterFromXYZ(dfalgo4bis)
  numcellrast <- nrow(dfalgo4bis)
  from <- c("1")
  to <- c("1")
  Dfadjacent <- data.frame(from,to)
  for(i in 1:numcellrast){
    ifelse(is.na(extract(Rasterdfalgo4bis,c(1,i))[2]) == "FALSE",
           ifelse(extract(Rasterdfalgo4bis,c(1,i))[2] == 1,
                  Dfadjacent<-rbind(Dfadjacent,adjacent(Rasterdfalgo4bis, cells=c(i),directions = 8, pairs = T)),0),0)
  }
  Dfadjacent <- Dfadjacent[-(1),] #suppression de la première ligne du df créée avec from = c(1)
  remove(from,to,i)
  Dfadjacent<-Dfadjacent[!duplicated(Dfadjacent[,2]),]
  Dfadjacent <- transform(Dfadjacent, to = as.numeric(to))
  for(i in 1:length(Dfadjacent$to)){
    n <- Dfadjacent$to[i]
    Dfadjacent$count[i] <- ifelse(is.na(Rasterdfalgo4bis[n]) == "TRUE", 1, 0)
    remove(n)
  }
  "Nombre de NA dans les nouvelles cellules adjacentes :"
  nrow(subset(Dfadjacent[c("count")], count == 1))
  for(i in 1:length(Dfadjacent$to)){
    n <- Dfadjacent$to[i]
    ifelse(is.na(Rasterdfalgo4bis[n]) == "FALSE",
           ifelse(Rasterdfalgo4bis[n] == 1,"",
                  Rasterdfalgo4bis[n] <- 4),
           Rasterdfalgo4bis[n] <- 5)
    remove(n)
  }
  remove(i)
  dfalgo4rastbis <- as.data.frame(rasterToPoints(Rasterdfalgo4bis,fun = NULL))
  dfalgo4rastbis <- dfalgo4rastbis %>% rename(celladjalgo4bis = Plantcondialgo4)
  for (i in 1:length(dfalgo4rastbis$x)){
    dfalgo4rastbis$x[i]<- round(x = dfalgo4rastbis$x[i],10)
    dfalgo4rastbis$y[i]<- round(x = dfalgo4rastbis$y[i],10)
  }#Permet d'arrondir les valeurs pour permettre la jointure
  df2All <- full_join(df2All, dfalgo4rastbis, by = c("x","y"))
  #On enlève les pixels qui se sont rajoutés :
  df2All <- df2All[is.na(df2All$GFW2010) == "FALSE", ]
  "Nombre de pixel étant adja aux pixels adja :"
  nrow(subset(df2All[c("celladjalgo4bis")], celladjalgo4bis >1))
  #La colonne celladjalgo4bis donne : pour les pixels = 4 ce sont les 
  #cellules adjacentes aux première cellules adjacentes et qui avaient des valeurs.
  #pour les pixels = 5 ce sont les cellules adja aux premières cellules adja
  #qui avaient des valeurs NA.
  #pour les pixels = 1 ce sont les premières cellules adjacentes ainsi que
  #les premières cellules d'après lesquels on a déterminé les premières cellules
  #adjacentes.
  "Vérification des conditions pour ces cellules adjacentes :"
  df2All$Condi2algo4 <- ifelse(df2All$celladjalgo4bis >1,
                               ifelse(df2All$Condtnbatprecitemp == 1,
                                      1,0),0)
  #En mettant celladjalgo4bis on prends en compte tous les types de pixels (pixels-forêts initaux,
  #premiers pixels adjacents, deuxième pixels adjacents)
  "Nombre de pixel répondant à cette condition : "
  nrow(subset(df2All[c("Condi2algo4")], Condi2algo4 == 1))
  "Colonne Condi2algo4 donne si =1 les cellules qui vont être reforestées."
  
  
  #1) compter pour les premiers pixels répondant conditions et étant considérés comme des forêts
  for(i in 1:length(df2All$x)){
    df2All$kmplantalgo4[i]<- ifelse(df2All$Plantindicealgo4[i] == 1,
                                    ifelse(df2All$GFW2010[i]<seuil2forestalgo4,
                                           (df2All$maxlnH2016[i]-df2All$uoppH2016[i])*(seuil2forestalgo4/100),
                                           df2All$Forestkm2[i]),0)
  }
  #2) compter pour les premiers pixels adja répondant conditions
  for(i in 1:length(df2All$x)){
    df2All$kmplantalgo4[i]<- ifelse(df2All$kmplantalgo4[i] == 0,
                                    ifelse(df2All$Plantalgo4[i] == 1,
                                           ifelse((df2All$GFW2010[i]/100)<(seuil3forestalgo4/100),
                                                  (df2All$maxlnH2016[i]-df2All$uoppH2016[i])*(seuil3forestalgo4/100),
                                                  df2All$Forestkm2[i]),df2All$kmplantalgo4[i]),df2All$kmplantalgo4[i])
  }
  #3) compter pour les deuxiemes pixels adja répondant conditions
  for(i in 1:length(df2All$x)){
    df2All$kmplantalgo4[i]<- ifelse(df2All$kmplantalgo4[i] == 0,
                                    ifelse(df2All$Condi2algo4[i] == 1,
                                           ifelse((df2All$GFW2010[i]/100)<(seuil4forestalgo4/100),
                                                  (df2All$maxlnH2016[i]-df2All$uoppH2016[i])*(seuil4forestalgo4/100),
                                                  df2All$Forestkm2[i]),df2All$Forestkm2[i]),df2All$kmplantalgo4[i])
  }
  #4) Ajouter les cellules NA ayant des valeurs pour Forestkm2
  for(i in 1:length(df2All$x)){
    df2All$kmplantalgo4[i] <- ifelse(is.na(df2All$kmplantalgo4[i]) == "FALSE",
                                     ifelse(df2All$kmplantalgo4[i]>df2All$Forestkm2[i],
                                            df2All$kmplantalgo4[i],df2All$Forestkm2[i]),
                                     df2All$Forestkm2[i])
  }
  
  print("% de couverture arborée atteint :")
  print((sum(df2All$kmplantalgo4, na.rm = T)/sum(df2All$maxlnH2016, na.rm = T))*100)
  
  #Obtention du nouveau pourcentage de forêt de chaque pixel
  for(i in 1:length(df2All$x)){
    df2All$GFW2010_new_41[i] <- (df2All$kmplantalgo4[i]/df2All$maxlnH2016[i])*100
  }
  
  #Pour avoir la différence de pourcentage de forêt d'un pixel entre avant reforestation
  #et après :
  for(i in 1:length(df2All$x)){
    df2All$diffreforesta41[i] <- ifelse(df2All$GFW2010_new_41[i]>df2All$GFW2010[i],
                                        (df2All$GFW2010_new_41[i]-df2All$GFW2010[i]),0)
  }
  
  #calcul des surfaces à planter
  for (i in 1:length(df2All$x)){
    df2All$diffalgo4[i] <- ifelse(df2All$kmplantalgo4[i] > df2All$Forestkm2[i],
                                  df2All$kmplantalgo4[i] - df2All$Forestkm2[i],0)
  }
  
  print("Nombre de km² de forêt plantés :")
  print(sum(df2All$diffalgo4, na.rm = T))
  
  ### 2) Calcul de la perte de production  ----
  for (i in 1:length(df2All$x)){
    df2All$loosevalalgo4[i] <- ifelse(df2All$diffalgo4[i]>0,
                                      df2All$SPAM2[i]*(df2All$diffalgo4[i]*100),0)
  }
  
  print("La perte de production agricole potentielle totale pour l'ensemble du Kenya pour une année ($):")
  print(sum(df2All$loosevalalgo4, na.rm = T))
  
  print("Perte de production agricole potentielle moyenne par ha ($/an):")
  print(sum(df2All$loosevalalgo4, na.rm = T)/(sum(df2All$diffalgo4, na.rm = T)*100))
  
  #Calcul de la perte de Quantitée produite (en tonne) et en se basant sur données SPAM uniquement: 
  for (i in 1:length(df2All$x)){
    df2All$loosequantialgo4[i]<- ifelse(df2All$diffalgo4[i]>0,
                                        df2All$SPAMProdquant2010_mtunit[i]*(df2All$diffalgo4[i]/df2All$maxlnH2016[i]),0)
  }
  
  #loosequantialgo4 donne pour chque pixel la perte de production agricole 
  #en termes de qttée (en tonne métrique) pour l'espace reforesté en se basant seulement
  #sur les données SPAM.
  
  print("La perte de production agricole potentielle totale pour l'ensemble du Kenya pour une année
  (en tonne et seulement basée sur les productions données par SPAM) :")
  print(sum(df2All$loosequantialgo4, na.rm = T))
  print("La perte de production agricole potentielle moyenne par hectare pour
        une année est (en tonne et basée juste sur données SPAM) :")
  print(sum(df2All$loosequantialgo4, na.rm = T)/(sum(df2All$diffalgo4, na.rm = T)*100))
  #ATTENTION : loosequantialgo4 ne contient pas la perte de production de l'élevage
  #il y a uniquement la perte de production agricole et provenant des données SPAM
  
  compteuralgo4part1 <- sum(df2All$kmplantalgo4, na.rm = T)/sum(df2All$maxlnH2016, na.rm = T)*100
  
  print("-_-_-_-__-_-_-_-_-_-_-__-_-_")
  print("-_-_-_-_-_ N E X T -_-_-__-_-")
  df2All <- df2All[ , !(names(df2All) %in% c("celladjalgo4bis","Condi2algo4","loosequantialgo4"))]
  df2All <- df2All[ , !(names(df2All) %in% c("celladjalgo4"))]
  df2All <- df2All[ , !(names(df2All) %in% c("Plantcondialgo4","Plantalgo4","Plantindicealgo4"))]
  remove(Dfadjacent, numcellrast,i)
  remove(Rastercomparalgo4bis,Rasterdfalgo4bis,dfalgo4bis,dfalgo4rastbis)
  if(numberalgo4part1<1){stop}
}


### B) Réalisation des cartes & graphs  ----
"Attention à supprimer les fichiers shp si l'on souhaite relancer les codes ci-dessous.
De plus, attention à bien les nommer dans la partie 'layer = ....'  "

# -- Carte représentant les espaces reforestés (km²) :
df4algo41 <- subset(df2All[c("x","y","diffalgo4")])
dshp <- df4algo41
coordinates(dshp)=~x+y
crs(dshp) <- crsref
dshp <- st_as_sf(dshp)
#Sortie carte
st_write(dshp, dsn = "diffalgo2a.shp", layer = "diffalgo4", driver = "ESRI Shapefile")
remove(df4algo41,dshp)

#ou en raster
dfalgo2dif  <- subset(df2All[c("x","y","diffalgo4")])
rastdiffalgo4 <- rasterFromXYZ(dfalgo2dif)
writeRaster(rastdiffalgo4, "rastalgo410%", format = "GTiff")
remove(rastdiffalgo4,dfalgo2dif)


# -- Carte représentant les pertes de production agricole potentielle en $  :
df4algo41loose <- subset(df2All[c("x","y","loosevalalgo4")])
dshp <- df4algo41loose
coordinates(dshp)=~x+y
proj4string(dshp) <- proj4string(Rastdecoupe)
#Sortie carte
writeOGR(dshp, dsn = 'coucheshpR', layer = "loosevalalgo4", driver = "ESRI Shapefile")
remove(df4algo41loose,dshp)

# -- Graph représentant les pertes de production agricole potentielle par counties
"Get the df with the conditions and transform it into shp"
dfcrossdatacounty2 <- subset(df2All[c("x","y","IDpixel","loosevalalgo4")])
coordinates(dfcrossdatacounty2) <- ~ x + y
shapecunties <- sf:::as_Spatial(shapecunties)
proj4string(dfcrossdatacounty2) <- proj4string(shapecunties)
dfcrossdatacounty2  <- st_as_sf(dfcrossdatacounty2)
shapecunties <- st_as_sf(shapecunties)
"Join the shapefil to get for each type of soil if there are pixel answering the conditions"
dfcheckcounties2 <- st_join(dfcrossdatacounty2, shapecunties,left = TRUE, largest = TRUE)
"Create the graphique"
#Faire la somme des loosevalalgo1 par county : 
uniquecountykenya <- unique(dfcheckcounties2$county)
dfsumcounty <- data.frame(c(1,1,1))
colnames(dfsumcounty)[colnames(dfsumcounty) == 'c.1..1..1.'] <- 'Sumagriloos'
dfcountyname <- data.frame(c(1,1,1))
colnames(dfcountyname)[colnames(dfcountyname) == 'c.1..1..1.'] <- 'Countyname'
for(i in 1:length(uniquecountykenya)){
  name<-uniquecountykenya[i]
  dfcountyname <- rbind(dfcountyname,name)
  sum1<-sum(subset(dfcheckcounties2[c("loosevalalgo4","county")], county == name)$loosevalalgo4)
  dfsumcounty <- rbind(dfsumcounty,sum1)
}
dfsumcounties<-data.frame(dfcountyname,dfsumcounty)
dfsumcounties <- dfsumcounties[-c(1,2,3),]
#erreur :
sum(df2All$loosevalalgo4,na.rm=T)-sum(dfsumcounties$Sumagriloos,na.rm=T)
sum(dfcheckcounties2$loosevalalgo4,na.rm=T)
#l'erreur vient du fait que certains pixels durant la jointure avec le shapefile
#des counties du Kenya ne sont pas référencés dans un county, et donc ont une NA
#dans la colonne county et donc les pertes de production agricole potentielle de ces
#pixels ne sont pas attribués à leur county.

#dfsumcounties donne, pour un niveau de tree cover (%) précisé
#par l'utilisateur et après que l'utilisateur ait lancé l'algorithme
# le niveau de perte totale de production agricole en $ par county du Kénya.
remove(dfsumcounty,dfcountyname,dfsumcounties,
       uniquecountykenya,dfcheckcounties2,shapecunties,
       dfcrossdatacounty2)


### C) Nettoyage ----
remove(dfalgo4, dfalgo4rast,Rastercomparalgo4,Rasterdfalgo4,dfalgo4bis,dfalgo4rastbis,
       Dfadjacent,Rastercomparalgo4bis,Rasterdfalgo4bis)
df2All <- df2All[ , !(names(df2All) %in% c("Plantationcondalgo4","celladjalgo4","diffalgo42"))]
df2All <- df2All[ , !(names(df2All) %in% c("diffalgo41","loosequantialgo4","km2plantalgo41"))]
df2All <- df2All[ , !(names(df2All) %in% c("GFW2010_new_41","diffreforesta41"))]
df2All <- df2All[ , !(names(df2All) %in% c("kmplantalgo4","Plantindicealgo4"))]
df2All <- df2All[ , !(names(df2All) %in% c("loosevalalgo4"))]
#Cleaning supplémentaire seulement si on n'utilise pas R_agroforestry ou R_code2 :
df2All <- df2All[ , !(names(df2All) %in% c("diffalgo4"))]




## III - Algorithme n°2 - partie B ----
"Algorithme de continuité écologique"


#Algorithme fonctionnant en visant la continuité écologique. Cet algorithme va chercher à
#reforester les pixels étant adjacents aux forêts existantes. Il va dans un premier temps sélectionner
#des pixels répondant correctement aux conditions (précipitation, surface urbanisée, température) et ayant
#un certain niveau de couverture arborée (seuilforestalgo42). Ces pixels vont donc être considérés comme des forêts.
#Ensuite il va regarder les pixels adjacents (méthode = "Investigating Commuting Time in a Metropolitan
#Statistical Area Using Spatial Autocorrelation Analysis" -> Queen's case - Hessam Miri) à ces pixels-forêt et
#regarder si ces premiers pixels adjacents respectent les conditions (précip, urbain, température). 
#Ensuite l'algorithme va regarder les pixels adjacents des premiers pixels adjacents respectant les conditions.
#Dans ces "deuxième pixels adjacents", l'algorithme va regarder les pixels qui respectent les conditions (précip,urbain,temp).
#Ce mécanisme de sélection des pixels adjacents et de vérification (relatif aux conditions) va être répété, toujours en se basant
#sur les pixels adjacents déterminés juste avant, jusqu'à l'atteinte de l'objectif de tree cover du Kenya souhaitée.
#Pour déterminer la nouvelle tree cover, cet algorithme augmente la tree cover (jusqu'à seuilforest42) de tous les pixels sélectionnés par l'adjacence et
#qui respectent les conditions (surface urbanisée,précipitation,température).
#Le seuil de forêt "seuilforest41" sert uniquement dans la première étape de cet algorithme et sert juste à donner
#une information sur les premiers pixels adjacents. Il n'est pas vraiment utile dans l'obtention d'une tree cover pour le Kenya.


#Explication de l'algorithme n°2 - partie b) : 
#L'utiliser lance les lignes ci-dessous de la première étape jusqu'à l'étape 2. 
#Ensuite, l'utilisateur rentre un objectif de % tree cover au Kenya à atteindre dans "obj".
#Le % de tree cover représente le % de couverture arborée sur l'ensemble du territoire Kenyan.
#Ensuite l'utilisateur lance la boucle. A chaque boucle l'algorithme sélectionne plus de pixel
# et donc, est capable d'atteindre une couverture arborée plus importante. 

#/!\
#Il est nécessaire de lancer la partie "Nettoyage" avant de relancer la boucle
#de l'algorithme avec d'autres conditions.


### A) Calcul de la Tree cover ----

#### 1) Première étape ----

"Avoir les pixels ayant une couverture arborée supérieure au seuilforestalgo42 :"
df2All$Plantindicealgo42 <- ifelse(df2All$Condtnbatprecitemp == 1,
                                 ifelse(df2All$GFW2010 > seuilforestalgo42,1,0),0)

"Nombre de pixel correspondant aux premiers critères de l'algorithme 4"
nrow(subset(df2All[c("Plantindicealgo42")], Plantindicealgo42 == 1))
nalgo41 <- nrow(subset(df2All[c("Plantindicealgo42")], Plantindicealgo42 == 1))
"Nombre de NA dans la colonne Plantindicealgo42 :"
sum(is.na(df2All$Plantindicealgo42))

#Obtention des premiers pixels adjacents par rapport aux cellules ayant 1 dans Plantindicealgo42
# -- a) Convertir en raster
dfalgo42 <- df2All[,c("x","y","Plantindicealgo42")]
Rastercomparalgo42 <-rasterFromXYZ(dfalgo42)#Raster de comparaison
Rasterdfalgo42 <- rasterFromXYZ(dfalgo42)
numcellrast <- nrow(dfalgo42)
# -- b) Création du df
from <- c("1")
to <- c("1")
Dfadjacent <- data.frame(from,to)
# -- c) Chaque cellule du raster est vérifiée, si elle correspond a Plantindicealgo42 ==1
#alors on sélectionne avec adjacent() les cellules étant adjacentes à celle-ci.
#Ensuite on intègre cela dans un df
for(i in 1:numcellrast){
  ifelse(is.na(extract(Rasterdfalgo42,c(1,i))[2]) == "FALSE",
         ifelse(extract(Rasterdfalgo42,c(1,i))[2] == 1,
                Dfadjacent<-rbind(Dfadjacent,adjacent(Rasterdfalgo42, cells=c(i),directions = 8, pairs = T)),0),0)
}
Dfadjacent <- Dfadjacent[-(1),] #suppression de la première ligne du df créée avec from = c(1)
remove(from,to,i)

# -- d) Supprimer les doublons et mise en forme
Dfadjacent<-Dfadjacent[!duplicated(Dfadjacent[,2]),]
Dfadjacent <- transform(Dfadjacent, to = as.numeric(to))

# -- e) Modification des valeurs du raster par rapport au numéro de cellule
"Code pour compter le nombre de NA dans les cellules adjacentes :"
for(i in 1:length(Dfadjacent$to)){
  n <- Dfadjacent$to[i]
  Dfadjacent$count[i] <- ifelse(is.na(Rasterdfalgo42[n]) == "TRUE", 1, 0)
  remove(n)
}

"Nombre de NA :"
nrow(subset(Dfadjacent[c("count")], count == 1))

"Cette fonction permet de convertir les numéros des cellules adjacentes récupérés (dans le raster)
et de les intégrer dans un raster pour avoir les coordonnées x et y de
ces cellules comme c'est le cas dans df2All (en sachant que les cellules adjacentes 
peuvent seulement prendre la valeur 2 si pas de NA sur le pixel initialement et
3 s'il y avait une NA sur le pixel initialement)."
for(i in 1:length(Dfadjacent$to)){
  n <- Dfadjacent$to[i]
  ifelse(is.na(Rasterdfalgo42[n]) == "FALSE",
         ifelse(Rasterdfalgo42[n] == 1,"",
                Rasterdfalgo42[n] <- 2),
         Rasterdfalgo42[n] <- 3)
  remove(n)
}
remove(i)

# -- f) Transformation du raster en dataframe
dfalgo42rast <- as.data.frame(rasterToPoints(Rasterdfalgo42,fun = NULL))

# -- g) Jointure
names(dfalgo42rast)[names(dfalgo42rast) == "Plantindicealgo42"] <- "celladjalgo4"
for (i in 1:length(dfalgo42rast$x)){
  dfalgo42rast$x[i]<- round(x = dfalgo42rast$x[i],10)
  dfalgo42rast$y[i]<- round(x = dfalgo42rast$y[i],10)
}#Permet d'arrondir les valeurs pour permettre la jointure
df2All <- full_join(df2All, dfalgo42rast, by = c("x","y"))
#On enlève les pixels qui se sont rajoutés :
df2All <- df2All[is.na(df2All$GFW2010) == "FALSE", ] 

"celladjalgo4 permet de donner les cellules adjacentes avec 2 et 3 en valeurs
dans la colonne."
"Nombre de première cellule adjacente :"
nrow(subset(df2All[c("celladjalgo4")], celladjalgo4 > 1))

"Vérification des conditions pour ces cellules adjacentes :"
df2All$Plantalgo42 <- ifelse(df2All$celladjalgo4 >1,
                            ifelse(df2All$Condtnbatprecitemp == 1,
                                   1,0),0)
"Nombre de pixel répondant à cette condition : "
nrow(subset(df2All[c("Plantalgo42")], Plantalgo42 == 1))
ntry <- nrow(subset(df2All[c("Plantalgo42")], Plantalgo42 == 1))
#Création d'une nouvelle colonne pour que les cellules adjacentes des cellules adjacentes
#ne puissent pas contenir les pixels de base que l'on a utilisé pour déterminer les
#premières cellules adjacentes.
df2All$Plantcondialgo4 <- ifelse(df2All$Plantalgo42 == 1,1,0)
df2All$Plantcondialgo4 <- ifelse(df2All$Plantcondialgo4 == 0,
                                 ifelse(df2All$Plantindicealgo42 == 1, 1,
                                        df2All$Plantcondialgo4),
                                 df2All$Plantcondialgo4)
"Plantalgo42 donne la liste  (si = 1) des pixels étant des cellules adjacentes et répondant
aux conditions de précipitation, température et surface bâtie."
"Plantcondialgo4 est pareille que Plantalgo42, mais dans Plantcondialgo4, les pixels = 1 comprennent
également les tout premiers pixels sélectionnés étant considérés comme des forêts."
"Plantcondialgo4 permet d'éviter le dédoublement des cellules adjacentes pour la suite."
remove(Dfadjacent, Rastercomparalgo42, Rasterdfalgo42,dfalgo42,dfalgo42rast)

#Calcul de la perte de production (donnant info sur premiers pixels+premiers pixel adja)
for(i in 1:length(df2All$x)){
  df2All$km2plantalgo41[i] <- ifelse(df2All$Plantcondialgo4[i] == 1,
                                     ifelse(df2All$GFW2010[i]<seuilforest41,
                                            (df2All$maxlnH2016[i]- df2All$uoppH2016[i])*(seuilforest41/100),
                                            df2All$Forestkm2[i]),df2All$Forestkm2[i])
}
"% de forêt au kenya avec ce seuil :"
sum(df2All$km2plantalgo41, na.rm = T)/sum(df2All$maxlnH2016, na.rm =T)*100

df2All$diffalgo41 <- ifelse(df2All$km2plantalgo41>df2All$Forestkm2,
                            df2All$km2plantalgo41-df2All$Forestkm2,
                            0)

#Vérification, problème si ce n'est pas = 0:
print("Vérification devant être égale à 0 : ")
print(nrow(subset(df2All[c("Plantindicealgo42")], Plantindicealgo42 == 1)) - (nrow(subset(df2All[c("Plantcondialgo4")], Plantcondialgo4 == 1))-ntry))
remove(ntry)




#### 2) Deuxième étape ----
#La boucle while ci-dessous va se répéter jusqu'a l'atteinte de l'objectif de tree cover 
#souhaité pour le Kenya. 
mo <- 1
obj <- 10
while(mo < obj){
  dfalgo4bis <- df2All[,c("x","y","Plantcondialgo4")]
  Rasterdfalgo4bis <- rasterFromXYZ(dfalgo4bis)
  numcellrast <- nrow(dfalgo4bis)
  from <- c("1")
  to <- c("1")
  Dfadjacent <- data.frame(from,to)
  for(i in 1:numcellrast){
    ifelse(is.na(extract(Rasterdfalgo4bis,c(1,i))[2]) == "FALSE",
           ifelse(extract(Rasterdfalgo4bis,c(1,i))[2] == 1,
                  Dfadjacent<-rbind(Dfadjacent,adjacent(Rasterdfalgo4bis, cells=c(i),directions = 8, pairs = T)),0),0)
  }
  Dfadjacent <- Dfadjacent[-(1),]
  remove(from,to,i)
  Dfadjacent<-Dfadjacent[!duplicated(Dfadjacent[,2]),]
  Dfadjacent <- transform(Dfadjacent, to = as.numeric(to))
  for(i in 1:length(Dfadjacent$to)){
    n <- Dfadjacent$to[i]
    ifelse(is.na(Rasterdfalgo4bis[n]) == "FALSE",
           ifelse(Rasterdfalgo4bis[n] == 1,"",
                  Rasterdfalgo4bis[n] <- 4),
           Rasterdfalgo4bis[n] <- 5)
    remove(n)
  }
  remove(i)
  dfalgo4rastbis <- as.data.frame(rasterToPoints(Rasterdfalgo4bis,fun = NULL))
  dfalgo4rastbis <- dfalgo4rastbis %>% rename("celladjalgo4bis" = "Plantcondialgo4")
  for (i in 1:length(dfalgo4rastbis$x)){
    dfalgo4rastbis$x[i]<- round(x = dfalgo4rastbis$x[i],10)
    dfalgo4rastbis$y[i]<- round(x = dfalgo4rastbis$y[i],10)
  }#Permet d'arrondir les valeurs pour permettre la jointure
  df2All <- full_join(df2All, dfalgo4rastbis, by = c("x","y"))
  df2All <- df2All[is.na(df2All$GFW2010) == "FALSE", ]#Supr pixel en trop
  #Vérification des condtitions cell adja 
  for(i in 1:length(df2All$x)){
  df2All$Plantationcondalgo4[i] <- ifelse(df2All$celladjalgo4bis[i] >1,
                                       ifelse(df2All$Condtnbatprecitemp[i] == 1,
                                              1,0),0)
  }
  #Plantationcondalgo4 ne se remplit qu'avec les dernieres cellules adjacentes, trouvées
  #par la dernière boucle qui a tourné, et respectant les conditions. Et la boucle va mettre
  #les nouvelles cellules adjacentes & respectant condtns dans la colonne Plantcondialgo4
  
  #il faut que dans Plantalgo42cond il y ait toutes les cellules adjacentes
  #à savoir les premiers pixels + les anciennes cellules adjacentes + les nouvelles cellules adjacentes
  for(i in 1:length(df2All$x)){
    df2All$Plantalgo42cond[i] <- ifelse(df2All$Plantcondialgo4[i] == 1,
                                        1,ifelse(df2All$Plantationcondalgo4[i] == 1,
                                                 1,0))
  }
  
  ##### a) Calcul des surfaces de forêt ----
  for(i in 1:length(df2All$x)){
    df2All$kmplantalgo42[i] <- ifelse(df2All$Plantalgo42cond[i] == 1,
                                        ifelse(df2All$GFW2010[i]<seuilforest42,
                                               (df2All$maxlnH2016[i]-df2All$uoppH2016[i])*(seuilforest42/100),
                                               df2All$Forestkm2[i]),df2All$Forestkm2[i])
  }
  #La colonne kmplantalgo42 permet de connaître pour tous les pixels (les premiers pixels considérées comme des forêts
  # + 1er pixels adjacents + tous les autres pixels adjacents) leur niveau de forêt en km² après
  #qu'ils aient été reforestés.
  print("% de couverture arborée atteint :")
  mo <- sum(df2All$kmplantalgo42, na.rm = T)/sum(df2All$maxlnH2016, na.rm = T)*100
  print(mo)
  
  #Obtention du nouveau pourcentage de couverture arborée pour chaque pixel
  for(i in 1:length(df2All$x)){
    df2All$GFW2010_new_42[i] <- (df2All$kmplantalgo42[i]/df2All$maxlnH2016[i])*100
  }
  
  #Obtention des nouvelles surfaces de forêts plantées par la reforestation
  for(i in 1:length(df2All$x)){
    df2All$diffalgo42[i] <- ifelse(df2All$kmplantalgo42[i]>df2All$Forestkm2[i],
                                   df2All$kmplantalgo42[i]-df2All$Forestkm2[i],0)
  }#La colonne diffalgo42 permet de connaître le nombre de km² de forêt planté par la reforestation
  #sur chaque pixel
  print("Nombre de km² de forêt plantés :")
  difdifalgo4 <- sum(df2All$diffalgo42,na.rm = T)
  print(difdifalgo4)
  #### b) Calcul de la perte de production agricole ----

  #Calcul de la perte de production en valeur :
  for (i in 1:length(df2All$x)){
    df2All$loosevalalgo42[i] <- ifelse(df2All$diffalgo42[i]>0,
                                      df2All$SPAM2[i]*(df2All$diffalgo42[i]*100),0)
  }
  perte1 <- sum(df2All$loosevalalgo42, na.rm = T)
  print("La perte de production agricole potentielle totale pour l'ensemble du Kenya pour une année ($):")
  print(perte1)
  print("Perte de production agricole potentielle moyenne par ha ($/an):")
  perte2 <- sum(df2All$loosevalalgo42, na.rm = T)/(sum(df2All$diffalgo42, na.rm = T)*100)
  print(perte2)
  
  #Calcul de la perte de production agricole en quantité (en tonne) :
  for (i in 1:length(df2All$x)){
    df2All$loosequantialgo42[i]<- ifelse(df2All$diffalgo42[i]>0,
                                        df2All$SPAMProdquant2010_mtunit[i]*(df2All$diffalgo42[i]/df2All$maxlnH2016[i]),0)
  }
  perte3 <- sum(df2All$loosequantialgo42,na.rm = T)
  print("La perte de production agricole potentielle totale pour l'ensemble du Kenya pour une année
  (en tonne et seulement basée sur les productions données par SPAM) :")
  print(perte3)
  print("La perte de production agricole potentielle moyenne par hectare pour
        une année est (en tonne et basée juste sur données SPAM) :")
  perte4 <- sum(df2All$loosequantialgo42, na.rm = T)/(sum(df2All$diffalgo42, na.rm = T)*100)
  print(perte4)
  print(" ------------------------------------")
  print("--------------__NEXT__---------------")
  df2All$Plantcondialgo4 <- df2All$Plantalgo42cond
  df2All <- df2All[ , !(names(df2All) %in% c("celladjalgo4bis",
                                             "loosequantialgo42",
                                            "Plantalgo42cond"))]
  remove(Rasterdfalgo4bis,numcellrast,dfalgo4bis,dfalgo4rastbis,Dfadjacent,i,perte1,perte2,
         perte3,perte4,difdifalgo4)
}


### B) Calcul d'une Tree cover spécifique ----


#La boucle ci-dessous permet d'obtenir un niveau de tree cover spécifique. Par exemple,
#si en lançant l'algorithme n°2-b) haut-dessus la tree cover obtenue est 18,46% et que l'utilisateur
#souhaite obtenir 18% il peut lancer la boucle ci-dessous qui va rétablir les surfaces de forêts
#initiales sur certains pixels.

#Explication du fonctionnement du code : 
#(La boucle ci-dessous, trie les pixels reforestés. Ensuite, la boucle va retirer le nombre
#de pixel définit pas l'utilisateur avec scalepoint, jusqu'à l'atteinte de l'objectif de tree cover (définit
#également par l'utilisateur dans "obj2").
#Après la boucle while, il est nécessaire de lancer le reste des codes en dessous.

df2All2 <- df2All
df2All2 <- df2All2 %>% arrange(desc(Plantationcondalgo4))
obj2 <- 9
scalepoint <- 6
while(mo > obj2){
  n <- 1
  while(n < scalepoint){
    if(df2All2$Plantationcondalgo4[n] == 1){
      df2All2$Plantcondialgo4[n] <- 0
      df2All2$Plantationcondalgo4[n] <- 0
    }
    n <- n+1
    df2All2 <- df2All2 %>% arrange(desc(Plantationcondalgo4))
  }
  for(i in 1:length(df2All2$x)){
    df2All2$kmplantalgo42[i] <- ifelse(df2All2$Plantcondialgo4[i] == 1,
                                      ifelse(df2All2$GFW2010[i]<seuilforest42,
                                             (df2All2$maxlnH2016[i]-df2All2$uoppH2016[i])*(seuilforest42/100),
                                             df2All2$Forestkm2[i]),df2All2$Forestkm2[i])
  }
  mo <- sum(df2All2$kmplantalgo42, na.rm = T)/sum(df2All2$maxlnH2016, na.rm = T)*100
  print("% de couverture arborée atteint :")
  print(mo)
  #Obtention du nouveau %age de forêt sur le pixel 
  for(i in 1:length(df2All2$x)){
    df2All2$GFW2010_new_42[i] <- (df2All2$kmplantalgo42[i]/df2All2$maxlnH2016[i])*100
  }
  #Obtention de la différence de %age de forêt avant reforestation et après par pixel
  for(i in 1:length(df2All2$x)){
    df2All2$diffreforesta42[i] <- ifelse(df2All2$GFW2010_new_42[i]>df2All2$GFW2010[i],
                                         (df2All2$GFW2010_new_42[i]-df2All2$GFW2010[i]),0)}
  
  #Obtention des nouvelles surfaces de forêts plantées 
  for(i in 1:length(df2All2$x)){
    df2All2$diffalgo42[i] <- ifelse(df2All2$kmplantalgo42[i]>df2All2$Forestkm2[i],
                                   df2All2$kmplantalgo42[i]-df2All2$Forestkm2[i],0)
  }
  print("Nombre de km² de forêt plantés :")
  difdifalgo4 <- sum(df2All2$diffalgo42,na.rm = T)
  print(difdifalgo4)
  #Calcul de la perte de production en valeur :
  for (i in 1:length(df2All2$x)){
    df2All2$loosevalalgo42[i] <- ifelse(df2All2$diffalgo42[i]>0,
                                        df2All2$SPAM2[i]*(df2All2$diffalgo42[i]*100),0)
  }
  perte1 <- sum(df2All2$loosevalalgo42, na.rm = T)
  print("La perte de production agricole potentielle totale pour l'ensemble du Kenya pour une année ($):")
  print(perte1)
  print("Perte de production agricole potentielle moyenne par ha ($/an):")
  perte2 <- sum(df2All2$loosevalalgo42, na.rm = T)/(sum(df2All2$diffalgo42, na.rm = T)*100)
  print(perte2)
  #Calcul de la perte de production agricole en quantité (en tonne) :
  for (i in 1:length(df2All2$x)){
    df2All2$loosequantialgo42[i]<- ifelse(df2All2$diffalgo42[i]>0,
                                          df2All2$SPAMProdquant2010_mtunit[i]*(df2All2$diffalgo42[i]/df2All2$maxlnH2016[i]),0)
  }
  perte3 <- sum(df2All2$loosequantialgo42,na.rm = T)
  print("La perte de production agricole potentielle totale pour l'ensemble du Kenya pour une année
  (en tonne et seulement basée sur les productions données par SPAM) :")
  print(perte3)
  print("La perte de production agricole potentielle moyenne par hectare pour
        une année est (en tonne et basée juste sur données SPAM) :")
  perte4 <- sum(df2All2$loosequantialgo42, na.rm = T)/(sum(df2All2$diffalgo42, na.rm = T)*100)
  print(perte4)
  print(" ------------------------------------")
  print("---__NEXT__---")
  df2All2 <- df2All2[ , !(names(df2All2) %in% c("loosequantialgo42"))]
  remove(i,perte1,perte2,perte3,perte4,difdifalgo4)
}
df2All$diffalgo42 <- df2All2$diffalgo42


### C) Réalisation des cartes & graphs ----
"Attention à supprimer les fichiers shp si l'on souhaite relancer les codes ci-dessous.
De plus, attention à bien les nommer dans la partie 'layer = ....'  "

# -- Carte représentant les pertes de production agricole potentielle en $  :
dfalgo4part2 <- subset(df2All2[c("x","y","loosevalalgo42")])
nam <- names(dfalgo4part2[3])
coordinates(dfalgo4part2)=~x+y
proj4string(dfalgo4part2) <- proj4string(Rastdecoupe)
#Sortie carte
writeOGR(dfalgo4part2, dsn = 'coucheshpR', layer = nam, driver = "ESRI Shapefile")
# -- Carte représentant les surfaces arborées après la reforestation en km²  :
dfalgo4part2a <- subset(df2All2[c("x","y","kmplantalgo42")])
nam <- names(dfalgo4part2a[3])
coordinates(dfalgo4part2a)=~x+y
proj4string(dfalgo4part2a) <- proj4string(Rastdecoupe)
#Sortie carte
writeOGR(dfalgo4part2a, dsn = 'coucheshpR', layer = nam, driver = "ESRI Shapefile")
# -- Carte représentant les espaces reforestés (km²) :
dfalgo4part3 <- subset(df2All[c("x","y","diffalgo42")])
dshp <- dfalgo4part3
coordinates(dshp)=~x+y
crs(dshp) <- crsref
dshp <- st_as_sf(dshp)
#Sortie carte
st_write(dshp, dsn = "diffalgo2b.shp", layer = "diffalgo42", driver = "ESRI Shapefile")
remove(dfalgo4part3,dshp)

#ou en raster :
dfalgo42dif  <- subset(df2All[c("x","y","diffalgo42")])
rastdiffalgo42 <- rasterFromXYZ(dfalgo42dif)
writeRaster(rastdiffalgo42, "rastalgo2partb10%", format = "GTiff")
remove(rastdiffalgo42,dfalgo42dif)

# -- Graph représentant les pertes de production agricole potentielle par counties
"Get the df with the conditions and transform it into shp"
dfcrossdatacounty3 <- subset(df2All[c("x","y","IDpixel","loosevalalgo42")])
coordinates(dfcrossdatacounty3) <- ~ x + y
shapecunties <- sf:::as_Spatial(shapecunties)
proj4string(dfcrossdatacounty3) <- proj4string(shapecunties)
dfcrossdatacounty3  <- st_as_sf(dfcrossdatacounty3)
shapecunties <- st_as_sf(shapecunties)
"Join the shapefil to get for each type of soil if there are pixel answering the conditions"
dfcheckcounties3 <- st_join(dfcrossdatacounty3, shapecunties,left = TRUE, largest = TRUE)
"Create the graphique"
#Faire la somme des loosevalalgo1 par county : 
uniquecountykenya <- unique(dfcheckcounties3$county)
dfsumcounty <- data.frame(c(1,1,1))
colnames(dfsumcounty)[colnames(dfsumcounty) == 'c.1..1..1.'] <- 'Sumagriloos'
dfcountyname <- data.frame(c(1,1,1))
colnames(dfcountyname)[colnames(dfcountyname) == 'c.1..1..1.'] <- 'Countyname'
for(i in 1:length(uniquecountykenya)){
  name<-uniquecountykenya[i]
  dfcountyname <- rbind(dfcountyname,name)
  sum1<-sum(subset(dfcheckcounties3[c("loosevalalgo42","county")], county == name)$loosevalalgo42)
  dfsumcounty <- rbind(dfsumcounty,sum1)
}
dfsumcounties<-data.frame(dfcountyname,dfsumcounty)
dfsumcounties <- dfsumcounties[-c(1,2,3),]
#erreur :
sum(df2All$loosevalalgo42,na.rm=T)-sum(dfsumcounties$Sumagriloos,na.rm=T)
sum(dfcheckcounties3$loosevalalgo42,na.rm=T)
#l'erreur vient du fait que certains pixels durant la jointure avec le shapefile
#des counties du Kenya ne sont pas référencés dans un county, et donc ont une NA
#dans la colonne county et donc les pertes de production agricole potentielle de ces
#pixels ne sont pas attribués à leur county.

#dfsumcounties donne, pour un niveau de tree cover (%) précisé
#par l'utilisateur et après que l'utilisateur ait lancé l'algorithme
# le niveau de perte totale de production agricole en $ par county du Kénya.
remove(dfsumcounty,dfcountyname,dfsumcounties,
       uniquecountykenya,dfcheckcounties3,shapecunties,
       dfcrossdatacounty3)


### D) Nettoyage ----
remove(Dfadjacent,dfalgo4bis,dfalgo4rastbis,Rasterdfalgo4bis,df2All2,n,nalgo41,
       obj,obj2,mo,ma,n)
df2All <- df2All[ , !(names(df2All) %in% c("Plantindicealgo42","celladjalgo4","Plantalgo4",
                                           "Plantcondialgo4","Plantationcondalgo4","celladjalgo4bis"))]
df2All <- df2All[ , !(names(df2All) %in% c("celladjalgo4bis","kmplantalgo42",
                                           "loosevalalgo42",
                                           "loosevalalgo42","loosequantialgo42",
                                           "km2plantalgo41","loosequantialgo41",
                                           "diffalgo41","Plantalgo42cond"))]
df2All <- df2All[ , !(names(df2All) %in% c("celladjalgo4.x","celladjalgo4.y"))]
df2All <- df2All[ , !(names(df2All) %in% c("Plantalgo42","GFW2010_new_42"))]
#Si l'utilisateur ne souhaite pas utiliser le code R_Agroforestry après, il peut
#lancer le nettoyage ci-dessous :
df2All <- df2All[ , !(names(df2All) %in% c("diffalgo42"))]



## IV - Catégorisation des pixels (non utilisé) ----
"L'objectif de cette partie 5 est de vérifier que l'augmentation
de forêt ne vient pas dépasser la surface disponible de chaque
pixel."


### a) - Première catégorisation des pixels  ----
#Cette catégorisation permet de voir s'il y a des problèmes de concordance 
#entre les données de HYDE et de Global Forest Watch étant donné qu'elles ne
#sont pas de la même année et que Global Forest Watch donne un % de forêt par pixel


#-- Catégorisation des éléments selon docRtreecover
for(i in 1:length(df2All$x)){
  df2All$forestselonHYDE[i] <- df2All$maxlnH2016[i] - df2All$croplandH2016[i] -
    df2All$grazingH2016[i] - df2All$uoppH2016[i]
}
sum(subset(df2All[c("forestselonHYDE")], forestselonHYDE < 0), na.rm = T)
nrow(subset(df2All[c("forestselonHYDE")], forestselonHYDE < 0))
#La colonne forestselonHYDE donne la surface de forêt/désert/Savanne qu'il
#reste dans le pixel en km2 après avoir retiré les terres agricoles et
#d'élevage et les surfaces bâties (il y a des erreurs de HYDE mais 
#ça ne représente que 5km2 sur les 583 311 du Kenya, ces erreurs ont
#forestselonHYDE <0) 
#Sur les 26 pixels concernés par cette erreur, tous ont des surfaces 
#de forêt > 1%  : 
dfforesthyde <- subset(df2All[c("x","forestselonHYDE","GFW2010",
                                "Forestkm2","croplandH2016",
                                "grazingH2016","uoppH2016","maxlnH2016")], forestselonHYDE < 0)
#Pour connaître nombre pixel ayant couv foret dépassant grazing dans les erreurs :
for(i in 1:length(dfforesthyde$x)){
  dfforesthyde$forestsurf1[i] <- ifelse(dfforesthyde$Forestkm2[i] > dfforesthyde$grazingH2016[i],
                                        1,0)}
nrow(subset(dfforesthyde[c("forestsurf1")], forestsurf1 == 1))
#Pour connaître le nombre de pixel ayant une couverture arborée dépassant la 
#surface de cropland et grazing au sein du pixel et donc allant sur la surface
#bâtie 
for(i in 1:length(dfforesthyde$x)){
  dfforesthyde$forestsurf[i] <- ifelse(dfforesthyde$Forestkm2[i] > (dfforesthyde$croplandH2016[i]+
                                                                      dfforesthyde$grazingH2016[i]),
                                       1,0)}
nrow(subset(dfforesthyde[c("forestsurf")], forestsurf == 1))
#Il n'y en a aucun donc pour ces pixels, on pose le fait que la surface de forêt
#se répartit sur la surface utilisée pour l'élevage et sur la surface agricole
#avec en priorité l'usage de la surface de l'élevage.
remove(dfforesthyde)

# Catégorisation - Pixel A - Case A, Case B, Case Z (Cas n°1)
for(i in 1:length(df2All$x)){
  df2All$Pixeltype[i] <- ifelse(df2All$forestselonHYDE[i]<0,"A",0)}
df2All$Casetype <- df2All$Pixeltype
for(i in 1:length(df2All$x)){df2All$Casetype[i]<-0}
for(i in 1:length(df2All$x)){
  ifelse(df2All$Pixeltype[i] == "A",
         ifelse(df2All$Forestkm2[i]>df2All$grazingH2016[i],
                ifelse(df2All$Forestkm2[i]>(df2All$grazingH2016[i]+
                                              df2All$croplandH2016[i]),
                       df2All$Casetype[i] <- "Z",
                       df2All$Casetype[i]<- "A"),
                df2All$Casetype[i]<-"B"),"")}

# Catégorisation - Pixel B - Case C, Case D, Case Z (Cas n°2)
for(i in 1:length(df2All$x)){
  df2All$Pixeltype[i] <- ifelse(df2All$forestselonHYDE[i] == 0,"B",df2All$Pixeltype[i])}
for(i in 1:length(df2All$x)){
  ifelse(df2All$Pixeltype[i] == "B",
         ifelse(df2All$Forestkm2[i]>df2All$grazingH2016[i],
                ifelse(df2All$Forestkm2[i]>(df2All$grazingH2016[i]+
                                              df2All$croplandH2016[i]),
                       df2All$Casetype[i] <- "Z",
                       df2All$Casetype[i]<- "D"),
                df2All$Casetype[i]<-"C"),"")}

# Catégorisation - Pixel C - Case E, Case F, Case G, Case Z (Cas n°3)
for(i in 1:length(df2All$x)){
  df2All$Pixeltype[i] <- ifelse(df2All$forestselonHYDE[i] > 0,"C",df2All$Pixeltype[i])}
for(i in 1:length(df2All$x)){
  ifelse(df2All$Pixeltype[i] == "C",
         ifelse(df2All$Forestkm2[i]>df2All$forestselonHYDE[i],
                ifelse(df2All$Forestkm2[i]>(df2All$forestselonHYDE[i]+df2All$grazingH2016[i]),
                       ifelse(df2All$Forestkm2[i]>(df2All$grazingH2016[i]+
                                                     df2All$croplandH2016[i]+df2All$forestselonHYDE[i]),
                              df2All$Casetype[i] <- "Z",
                              df2All$Casetype[i]<- "G"),
                       df2All$Casetype[i]<-"F"),
                df2All$Casetype[i]<-"E"),"")}

# Obtention des zones mixtes (Cropland & Forêt, Grazing & Forêt) :
#  --  Définition des colonnes de zones mixtes
df2All$surfacemixtpat <- df2All$grazingH2016 #surfacemixtpat = surface en km2 où il y a des forêts
#sur des patures
for(i in 1:length(df2All$x)){df2All$surfacemixtpat[i] <- 0}
df2All$surfacemixtcrop <- df2All$surfacemixtpat    #surfacemixtcrop = surface en km2 où il y a des
#forêts sur des terres agricoles
df2All$surfaceplus <- df2All$surfacemixtpat#pour les case Z permet de connaître la surface en
#km2 qui empiète sur les surfaces bâties (en sachant que surface crop + grazing occupée)
df2All$Celluletype <- df2All$surfacemixtpat # la colonne Celluletype permet de catégoriser
#les différents types de pixel une fois que l'on a calculé l'augmentation du couvert
#arboré
df2All$surfacemixtcropplus <- df2All$surfacemixtpat #surfacemixtcropplus permet d'avoir la surface
#en km2 où il y a des forêts sur des terres agricoles APRES AVOIR AUGMENTER LE COUVERT arboré
df2All$surfacemixtpatplus <- df2All$surfacemixtpat #surfacemixtpatplus permet d'avoir la surface
#en km2 où il y a des forêts sur des terres agricoles APRES AVOIR AUGMENTER LE COUVERT arboré
df2All$surfaceplus2 <- df2All$surfacemixtpat #pour les cases Z permet de connaître la surface
#en km2 qui empiète sur les surfaces bâties (en sachant que surface crop+grazing occupée)


#  --  Calcul des zones mixtes (Cas n°1)
for(i in 1:length(df2All$x)){
  ifelse(df2All$Casetype[i] == "B",df2All$surfacemixtpat[i] <- df2All$Forestkm2[i],"")}
for(i in 1:length(df2All$x)){
  ifelse(df2All$Casetype[i] == "A", df2All$surfacemixtpat[i] <- df2All$grazingH2016[i],"")}
for(i in 1:length(df2All$x)){
  ifelse(df2All$Casetype[i] == "A", df2All$surfacemixtcrop[i]<-(df2All$Forestkm2[i]-df2All$grazingH2016[i]),"")}

#  --  Calcul des zones mixtes (Cas n°2)
for(i in 1:length(df2All$x)){
  ifelse(df2All$Casetype[i] == "C",df2All$surfacemixtpat[i] <- df2All$Forestkm2[i],"")}
for(i in 1:length(df2All$x)){
  ifelse(df2All$Casetype[i] == "D",df2All$surfacemixtpat[i] <- df2All$grazingH2016[i],"")}
for(i in 1:length(df2All$x)){
  ifelse(df2All$Casetype[i] == "D",df2All$surfacemixtcrop[i] <- (df2All$Forestkm2[i]-df2All$grazingH2016[i]),"")}

#  --  Calcul des zones mixtes (Cas n°3)
for(i in 1:length(df2All$x)){
  ifelse(df2All$Casetype[i] == "E",df2All$surfacemixtpat[i] <- 0,"")}
for(i in 1:length(df2All$x)){
  ifelse(df2All$Casetype[i] == "F",df2All$surfacemixtpat[i] <- (df2All$Forestkm2[i]-df2All$forestselonHYDE[i]),"")}
for(i in 1:length(df2All$x)){
  ifelse(df2All$Casetype[i] == "G",df2All$surfacemixtpat[i] <- df2All$grazingH2016[i],"")}
for(i in 1:length(df2All$x)){
  ifelse(df2All$Casetype[i] == "G",df2All$surfacemixtcrop[i] <- (df2All$Forestkm2[i]-df2All$forestselonHYDE[i]-
                                                                   df2All$grazingH2016[i]),"")}

# Pour case Z où Forest km2 > cropland + grazing :
#Pour info , pour les cases Z la forêt prend l'ensemble de l'espace de cropland et de grazing et plus
for(i in 1:length(df2All$x)){
  ifelse(df2All$Casetype[i] == "Z",df2All$surfaceplus[i]<- (df2All$Forestkm2[i]-df2All$grazingH2016[i]-
                                                              df2All$croplandH2016[i]),"")}
for(i in 1:length(df2All$x)){
  ifelse(df2All$Pixeltype[i] == "C",
         ifelse(df2All$Casetype[i] == "Z",df2All$surfaceplus[i] <- (df2All$Forestkm2[i]-df2All$forestselonHYDE[i]-
                                                                      df2All$grazingH2016[i]-df2All$croplandH2016[i]),""),"")}
for(i in 1:length(df2All$x)){
  if(df2All$Casetype[i] == "Z"){
    df2All$surfacemixtpat[i] <- df2All$grazingH2016[i]
    df2All$surfacemixtcrop[i] <- df2All$croplandH2016[i]}}

#Ajout des NA pour les valeurs vides:
for(i in 1:length(df2All$x)){ifelse(df2All$Casetype[i] == 0, 
                                    df2All$Casetype[i] <- NA,"")}


#Vérification pour la répartition de Forestkm2 :
#Pixel A, si résultat pas = 0, il y a un problème
sum(subset(df2All[c("Pixeltype","surfacemixtpat",
                    "surfacemixtcrop","surfaceplus")], Pixeltype == "A")$surfacemixtcrop, na.rm = T)+
  sum(subset(df2All[c("Pixeltype","surfacemixtpat",
                      "surfacemixtcrop","surfaceplus")], Pixeltype == "A")$surfacemixtpat, na.rm = T)+
  sum(subset(df2All[c("Pixeltype","surfacemixtpat",
                      "surfacemixtcrop","surfaceplus")], Pixeltype == "A")$surfaceplus, na.rm = T)-
  sum(subset(df2All[c("Pixeltype","Forestkm2")], Pixeltype == "A")$Forestkm2, na.rm = T)
#Pixel B, si résultat pas = 0, il y a un problème
sum(subset(df2All[c("Pixeltype","surfacemixtpat",
                    "surfacemixtcrop","surfaceplus")], Pixeltype == "B")$surfacemixtcrop, na.rm = T)+
  sum(subset(df2All[c("Pixeltype","surfacemixtpat",
                      "surfacemixtcrop","surfaceplus")], Pixeltype == "B")$surfacemixtpat, na.rm = T)+
  sum(subset(df2All[c("Pixeltype","surfacemixtpat",
                      "surfacemixtcrop","surfaceplus")], Pixeltype == "B")$surfaceplus, na.rm = T)-
  sum(subset(df2All[c("Pixeltype","Forestkm2")], Pixeltype == "B")$Forestkm2, na.rm = T)
#Pixel C, si résultat pas = 0, il y a un problème
sum(subset(df2All[c("Pixeltype","surfacemixtpat",
                    "surfacemixtcrop","surfaceplus")], Pixeltype == "B")$surfacemixtcrop, na.rm = T)+
  sum(subset(df2All[c("Pixeltype","surfacemixtpat",
                      "surfacemixtcrop","surfaceplus")], Pixeltype == "B")$surfacemixtpat, na.rm = T)+
  sum(subset(df2All[c("Pixeltype","surfacemixtpat",
                      "surfacemixtcrop","surfaceplus")], Pixeltype == "B")$surfaceplus, na.rm = T)+
  sum(subset(df2All[c("Pixeltype","surfacemixtpat","forestselonHYDE",
                      "surfacemixtcrop","surfaceplus")], Pixeltype == "B")$forestselonHYDE, na.rm = T)-
  sum(subset(df2All[c("Pixeltype","Forestkm2")], Pixeltype == "B")$Forestkm2, na.rm = T)

"Combien de pixel catégorie A :"
nrow(subset(df2All[c("Pixeltype")], Pixeltype == "A"))
"Combien de pixel catégorie B :"
nrow(subset(df2All[c("Pixeltype")], Pixeltype == "B"))
"Combien de pixel catégorie C :"
nrow(subset(df2All[c("Pixeltype")], Pixeltype == "C"))
"Combien de pixel de la catégorie A, empiète sur les terres agricoles ?"
nrow(subset(df2All[c("Casetype")], Casetype == "A"))
"Combien de km² ça représente sur ces terres agricoles ?"
sum(subset(df2All[c("Casetype","surfacemixtcrop")], Casetype == "A")$surfacemixtcrop, na.rm = T)
#Par contre l'ensemble de la surface du paturage est occupé par les forêts
"Combien de pixel de la catégorie B, empiète sur les terres agricoles ?"
nrow(subset(df2All[c("Casetype")], Casetype == "D"))
"Combien de km² ça représente sur ces terres agricoles ?"
sum(subset(df2All[c("Casetype","surfacemixtcrop")], Casetype == "D")$surfacemixtcrop, na.rm = T)
"Combien de pixel de la catégorie C, empiète sur les terres agricoles ?"
nrow(subset(df2All[c("Casetype")], Casetype == "G"))
"Combien de km² ça représente sur ces terres agricoles ?"
sum(subset(df2All[c("Casetype","surfacemixtcrop")], Casetype == "G")$surfacemixtcrop, na.rm = T)
"Combien de pixel empiète sur les surfaces bâties ?"
nrow(subset(df2All[c("Casetype")], Casetype == "Z"))
"Combien de km² ça représente sur ces surfaces bâties ?"
sum(subset(df2All[c("Casetype","surfaceplus")], Casetype == "Z")$surfaceplus, na.rm = T)
"Nombre de km2 total de forêt empiétant sur les terres agricoles ?"
sum(subset(df2All[c("Casetype","surfacemixtcrop")], Casetype == "G")$surfacemixtcrop, na.rm = T)+
  sum(subset(df2All[c("Casetype","surfacemixtcrop")], Casetype == "A")$surfacemixtcrop, na.rm = T)+
  sum(subset(df2All[c("Casetype","surfacemixtcrop")], Casetype == "D")$surfacemixtcrop, na.rm = T)
"Combien de km2 total sont sur les surfaces de pâturage ?"
sommepaturkm2 <- sum(subset(df2All[c("Casetype","surfacemixtpat")], Casetype == "A")$surfacemixtpat, na.rm = T)+
  sum(subset(df2All[c("Casetype","surfacemixtpat")], Casetype == "B")$surfacemixtpat, na.rm = T)+
  sum(subset(df2All[c("Casetype","surfacemixtpat")], Casetype == "C")$surfacemixtpat, na.rm = T)+
  sum(subset(df2All[c("Casetype","surfacemixtpat")], Casetype == "D")$surfacemixtpat, na.rm = T)+
  sum(subset(df2All[c("Casetype","surfacemixtpat")], Casetype == "F")$surfacemixtpat, na.rm = T)+
  sum(subset(df2All[c("Casetype","surfacemixtpat")], Casetype == "G")$surfacemixtpat, na.rm = T)
sommepaturkm2
"Quel pourcentage ça représente sur l'ensemble de la surface en pâturage ?"
sommepaturkm2/(sum(df2All$grazingH2016, na.rm = T))*100





### b) - Deuxième catégorisation des pixels ----

#Le code suivant est à faire tourner après l'utilisation d'un algorithme
#ayant permis de déterminer les surfaces de forêts en plus à reforester
#______________________________________________________________________
#Dans ce code, avant chaque utilisation il est nécessaire d'ajuster
#la colonne "diffalgo" en fonction de l'algorithme par exemple
#diffalgo1 concerne  l'algorithme 1 ainsi que le nom de la base
#par exemple df2Allalg1 pour l'algorithme 1.
#diffalgo4 pour algo 2 part a, diffalgo42 pour algo 2 part b 


nombase <- df2All #permet de donner le nom de la base pour chaque algorithme

# -- Catégorisation des pixels avec diffalgo - Pixel A, Case B, Cellule Y,H,I :
for(i in 1:length(nombase$x)){
  ifelse(nombase$Pixeltype[i] == "A",
         ifelse(nombase$Casetype[i] == "B",
                ifelse((nombase$diffalgo42[i]+nombase$Forestkm2[i])>nombase$grazingH2016[i],
                       ifelse((nombase$diffalgo42[i]+nombase$Forestkm2[i])>(nombase$grazingH2016[i]+
                                                                             nombase$croplandH2016[i]),
                              nombase$Celluletype[i] <- "Y",
                              nombase$Celluletype[i] <- "H"),
                       nombase$Celluletype[i] <- "I"),""),"")
}
# -- Catégorisation des pixels avec diffalgo - Pixel A, Case A, Cellule Y,J : 
for(i in 1:length(nombase$x)){
  ifelse(nombase$Pixeltype[i] == "A",
         ifelse(nombase$Casetype[i] == "A",
                ifelse((nombase$diffalgo42[i]+nombase$Forestkm2[i])>(nombase$grazingH2016[i]+
                                                                      nombase$croplandH2016[i]),
                       nombase$Celluletype[i] <- "Y",
                       nombase$Celluletype[i] <- "J"),""),"")
}
# -- Catégorisation des pixels avec diffalgo - Pixel B, Case C, Cellule K,L
for(i in 1:length(nombase$x)){
  ifelse(nombase$Pixeltype[i] == "B",
         ifelse(nombase$Casetype[i] == "C",
                ifelse((nombase$diffalgo42[i]+nombase$Forestkm2[i])>nombase$grazingH2016[i],
                       ifelse((nombase$diffalgo42[i]+nombase$Forestkm2[i])>(nombase$grazingH2016[i]+
                                                                             nombase$croplandH2016[i]),
                              nombase$Celluletype[i] <- "Y",
                              nombase$Celluletype[i] <- "L"),
                       nombase$Celluletype[i] <- "K"),""),"")
}
# -- Catégorisation des pixels avec diffalgo - Pixel B, Case D, Cellule Y,M
for(i in 1:length(nombase$x)){
  ifelse(nombase$Pixeltype[i] == "B",
         ifelse(nombase$Casetype[i] == "D",
                ifelse((nombase$diffalgo42[i]+nombase$Forestkm2[i])>(nombase$grazingH2016[i]+
                                                                      nombase$croplandH2016[i]),
                       nombase$Celluletype[i] <- "Y",
                       nombase$Celluletype[i] <- "M"),""),"")
}
# -- Catégorisation des pixels avec diffalgo - Pixel C, Case E, Cellule Y,O,P,N
for(i in 1:length(nombase$x)){
  ifelse(nombase$Pixeltype[i] == "C",
         ifelse(nombase$Casetype[i] == "E",
                ifelse((nombase$diffalgo42[i]+nombase$Forestkm2[i])>nombase$forestselonHYDE[i],
                       ifelse((nombase$diffalgo42[i]+nombase$Forestkm2[i])>(nombase$grazingH2016[i]+
                                                                             nombase$forestselonHYDE[i]),
                              ifelse((nombase$diffalgo42[i]+nombase$Forestkm2[i])>(nombase$grazingH2016[i]+
                                                                                    nombase$forestselonHYDE[i]+
                                                                                    nombase$croplandH2016[i]),
                                     nombase$Celluletype[i] <- "Y",
                                     nombase$Celluletype[i] <- "P"),
                              nombase$Celluletype[i] <- "O"),
                       nombase$Celluletype[i] <- "N"),""),"")
}
# -- Catégorisation des pixels avec diffalgo - Pixel C, Case F, Cellule Y,R,Q
for(i in 1:length(nombase$x)){
  ifelse(nombase$Pixeltype[i] == "C",
         ifelse(nombase$Casetype[i] == "F",
                ifelse((nombase$diffalgo42[i]+nombase$Forestkm2[i])>(nombase$grazingH2016[i]+nombase$forestselonHYDE[i]),
                       ifelse((nombase$diffalgo42[i]+nombase$Forestkm2[i])>(nombase$grazingH2016[i]+
                                                                             nombase$croplandH2016[i]+
                                                                             nombase$forestselonHYDE[i]),
                              nombase$Celluletype[i] <- "Y",
                              nombase$Celluletype[i] <- "R"),
                       nombase$Celluletype[i] <- "Q"),""),"")
}
# -- Catégorisation des pixels avec diffalgo - Pixel C, Case G, Cellule Y,S
for(i in 1:length(nombase$x)){
  ifelse(nombase$Pixeltype[i] == "C",
         ifelse(nombase$Casetype[i] == "G",
                ifelse((nombase$diffalgo42[i]+nombase$Forestkm2[i])>(nombase$grazingH2016[i]+
                                                                      nombase$croplandH2016[i]+
                                                                      nombase$forestselonHYDE[i]),
                       nombase$Celluletype[i] <- "Y",
                       nombase$Celluletype[i] <- "S"),""),"")
}

# -- Catégorisation pour les pixels ayant Case Z : 
for(i in 1:length(nombase$x)){ ifelse(nombase$Casetype[i] == "Z",
                                      nombase$Celluletype[i] <- "Y","")}

# -- Catégorisation des pixels NA : 
#Ajout des NA pour les valeurs vides:
for(i in 1:length(nombase$x)){ifelse(nombase$Celluletype[i] == 0, 
                                    nombase$Celluletype[i] <- NA,"")}

#  --  Calcul des zones mixtes (Situation n°1)
for(i in 1:length(nombase$x)){
  ifelse(nombase$Celluletype[i] == "I",nombase$surfacemixtpatplus[i] <- (nombase$Forestkm2[i]+
                                                                         nombase$diffalgo42[i]),"")}
for(i in 1:length(nombase$x)){ifelse(nombase$Celluletype[i] == "H",nombase$surfacemixtpatplus[i] <- nombase$grazingH2016[i],"")}
for(i in 1:length(nombase$x)){
  ifelse(nombase$Celluletype[i] == "H",nombase$surfacemixtcropplus[i] <- (nombase$Forestkm2[i]+
                                                                          nombase$diffalgo42[i]-
                                                                          nombase$grazingH2016[i]),"")}
for(i in 1:length(nombase$x)){
  ifelse(nombase$Celluletype[i] == "J",nombase$surfacemixtpatplus[i] <- nombase$grazingH2016[i],"")}
for(i in 1:length(nombase$x)){
  ifelse(nombase$Celluletype[i] == "J",nombase$surfacemixtcropplus[i] <- (nombase$Forestkm2[i]+
                                                                          nombase$diffalgo42[i]-
                                                                          nombase$grazingH2016[i]),"")}
#  --  Calcul des zones mixtes (Situation n°2)
for(i in 1:length(nombase$x)){
  ifelse(nombase$Celluletype[i] == "K",nombase$surfacemixtpatplus[i] <- (nombase$Forestkm2[i]+
                                                                         nombase$diffalgo42[i]),"")}
for(i in 1:length(nombase$x)){
  ifelse(nombase$Celluletype[i] == "L",nombase$surfacemixtpatplus[i] <- nombase$grazingH2016[i],"")}
for(i in 1:length(nombase$x)){
  ifelse(nombase$Celluletype[i] == "L",nombase$surfacemixtcropplus[i] <- (nombase$Forestkm2[i]+
                                                                          nombase$diffalgo42[i]-
                                                                          nombase$grazingH2016[i]),"")}
for(i in 1:length(nombase$x)){
  ifelse(nombase$Celluletype[i] == "M",nombase$surfacemixtpatplus[i] <- nombase$grazingH2016[i],"")}
for(i in 1:length(nombase$x)){
  ifelse(nombase$Celluletype[i] == "M",nombase$surfacemixtcropplus[i] <- (nombase$Forestkm2[i]+
                                                                          nombase$diffalgo42[i]-
                                                                          nombase$grazingH2016[i]),"")}
#  --  Calcul des zones mixtes (Situation n°3)
for(i in 1:length(nombase$x)){
  ifelse(nombase$Celluletype[i] == "Q",nombase$surfacemixtpatplus[i] <- (nombase$diffalgo42[i]+
                                                                         nombase$Forestkm2[i]-
                                                                         nombase$forestselonHYDE[i]),"")}
for(i in 1:length(nombase$x)){
  ifelse(nombase$Celluletype[i] == "R",nombase$surfacemixtpatplus[i] <- nombase$grazingH2016[i],"")}
for(i in 1:length(nombase$x)){
  ifelse(nombase$Celluletype[i] == "R",nombase$surfacemixtcropplus[i] <- (nombase$Forestkm2[i]+
                                                                          nombase$diffalgo42[i]
                                                                        -nombase$grazingH2016[i]
                                                                        -nombase$forestselonHYDE[i]),"")}
for(i in 1:length(nombase$x)){
  ifelse(nombase$Celluletype[i] == "S",nombase$surfacemixtpatplus[i] <- nombase$grazingH2016[i],"")}
for(i in 1:length(nombase$x)){
  ifelse(nombase$Celluletype[i] == "S",nombase$surfacemixtcropplus[i] <- (nombase$Forestkm2[i]+
                                                                          nombase$diffalgo42[i]
                                                                        -nombase$grazingH2016[i]
                                                                        -nombase$forestselonHYDE[i]),"")}
for(i in 1:length(nombase$x)){
  ifelse(nombase$Celluletype[i] == "N",nombase$surfacemixtpatplus[i] <- 0,"")}                                                                      
for(i in 1:length(nombase$x)){
  ifelse(nombase$Celluletype[i] == "O",nombase$surfacemixtpatplus[i] <- (nombase$Forestkm2[i]+
                                                                         nombase$diffalgo42[i]-
                                                                         nombase$forestselonHYDE[i]),"")} 
for(i in 1:length(nombase$x)){
  ifelse(nombase$Celluletype[i] == "P",nombase$surfacemixtpatplus[i] <- nombase$grazingH2016[i],"")}
for(i in 1:length(nombase$x)){
  ifelse(nombase$Celluletype[i] == "P",nombase$surfacemixtcropplus[i] <- (nombase$Forestkm2[i]+
                                                                          nombase$diffalgo42[i]-
                                                                          nombase$forestselonHYDE[i]
                                                                        -nombase$grazingH2016[i]),"")}
# -- Calcul des zones mixtes pour les cellules Y (erreur)

#Pour l'ensemble des pixels : 

for(i in 1:length(nombase$x)){
  ifelse(nombase$Celluletype[i] == "Y",nombase$surfacemixtpatplus[i] <- nombase$grazingH2016[i],"")}
for(i in 1:length(nombase$x)){
  ifelse(nombase$Celluletype[i] == "Y",nombase$surfacemixtcropplus[i] <- nombase$croplandH2016[i],"")}
#Ici on sait que Y est une erreur et donc que les surfaces en crop ou paturage sont remplies

#Pour les pixels C :
for(i in 1:length(nombase$x)){
  ifelse(nombase$Pixeltype[i] == "C",
         ifelse(nombase$Celluletype[i] == "Y",nombase$surfaceplus2[i] <- (nombase$diffalgo42[i]+nombase$Forestkm2[i]-
                                                                          nombase$grazingH2016[i]-
                                                                          nombase$croplandH2016[i]-
                                                                          nombase$forestselonHYDE[i]),""),"")}

#Pour les pixels B :
for(i in 1:length(nombase$x)){
  ifelse(nombase$Pixeltype[i] == "B",
         ifelse(nombase$Celluletype[i] == "Y",nombase$surfaceplus2[i] <- (nombase$diffalgo42[i]+nombase$Forestkm2[i]-
                                                                          nombase$grazingH2016[i]-
                                                                          nombase$croplandH2016[i]),""),"")}
#Pour les pixels A :
for(i in 1:length(nombase$x)){
  ifelse(nombase$Pixeltype[i] == "A",
         ifelse(nombase$Celluletype[i] == "Y",nombase$surfaceplus2[i] <- (nombase$diffalgo42[i]+nombase$Forestkm2[i]-
                                                                          nombase$grazingH2016[i]-
                                                                          nombase$croplandH2016[i]),""),"")}
#La colonne "surfaceplus2" permet de connaître de combien de km2 les aires rentrés
#par les données dépassent sur les surfaces bâties.

#La colonne "surfacemixtcropplus" permet de connaître la surface en km² où il y a des
#culture et de la forêt
#La colonne "surfacemixtpatplus" permet de connaître la surface en km² où il y a des
#pâturage et de la forêt

#Information supplémentaire sur la catégorisation des pixels après la reforestation :
"Nombre de cellules Y sur l'ensemble des pixels du Kenya :"
nrow(subset(nombase[c("Celluletype")], Celluletype == "Y"))
"La surface totale provenant des cellules Y et empiétant sur les surfaces bâties :"
sum(subset(nombase[c("Celluletype","surfaceplus2")], Celluletype == "Y")$surfaceplus2, na.rm = T)

"Pour l'algorithme 1 l'erreur semble très faible : 1 pixels ayant Y, pareille pour l'algo 2 part a,
pareille pour algo 2 part b"


#Nettoyage :
remove(nombase)



## V - Autres gains externes ----
#Dans cette partie il y a les gains qu'on a cherché a déterminé pour obtenir de l'information et
#pour évalué les gains dues à la hausse de la précipitation. 

### a) - Gain de production voisin forêt ----
"Objectif : déterminer les gains de production des cellules étant adjacentes
au pixel ayant une forte couverture arborée (> au seuilforestgain1."
seuilforestgain1 <- quantile(subset(df2All[c("GFW2010")],GFW2010>0), probs=seq(0, 1, 0.01), na.rm = T)[[76]]
# a) Calcul des conditions -
df2All$celladjagain1 <- ifelse(df2All$GFW2010 > seuilforestgain1,1,0)
# b) Cellules adjacentes -
dfgain1 <- df2All[,c("x","y","celladjagain1")]
Rastergain1 <- rasterFromXYZ(dfgain1)
numcellrast <- nrow(dfgain1)
from <- c("1")
to <- c("1")
Dfadjacent <- data.frame(from,to)
for(i in 1:numcellrast){
  ifelse(is.na(extract(Rastergain1,c(1,i))[2]) == "FALSE",
         ifelse(extract(Rastergain1,c(1,i))[2] == 1,
                Dfadjacent<-rbind(Dfadjacent,adjacent(Rastergain1, cells=c(i),directions = 8, pairs = T)),0),0)
}
Dfadjacent <- Dfadjacent[-(1),]
remove(from,to,i)
Dfadjacent<-Dfadjacent[!duplicated(Dfadjacent[,2]),]
Dfadjacent <- transform(Dfadjacent, to = as.numeric(to))
for(i in 1:length(Dfadjacent$to)){
  n <- Dfadjacent$to[i]
  ifelse(is.na(Rastergain1[n]) == "FALSE",
         ifelse(Rastergain1[n] == 1,"",
                Rastergain1[n] <- 2),
         Rastergain1[n] <- 3)
  remove(n)
}
remove(i)
dfgain1rast <- as.data.frame(rasterToPoints(Rastergain1,fun = NULL))
dfgain1rast <- dfgain1rast %>% rename(cellgain1 = celladjagain1)
#dans la colonne cellgain1 il y a les premiers pixels valant 1 (c-à-d ceux
#qui avaient GFW > à seuilforestgain1), les pixels valant 2 (c_à_d les pixels
#étant adjacent aux premiers pixels et ayant pas de NA) et les pixels valant 3
#Qui eux ont une valeur NA mais sont adjacents à un premier pixel
for (i in 1:length(dfgain1rast$x)){
  dfgain1rast$x[i]<- round(x = dfgain1rast$x[i],10)
  dfgain1rast$y[i]<- round(x = dfgain1rast$y[i],10)
}
df2All <- full_join(df2All, dfgain1rast, by = c("x","y"))
df2All <- df2All[is.na(df2All$GFW2010) == "FALSE", ] 
remove(Dfadjacent,i,numcellrast,Rastergain1,dfgain1,dfgain1rast)

for(i in 1:length(df2All$x)){
  df2All$celladjagain2[i] <- ifelse(df2All$cellgain1[i]>1,1,0)
}
"celladjagain2 donne les pixels étant adjacents aux pixels qui répondent
au critère de seuilforestgain1."
# c) Analyse des cellules -
"La moyenne SPAM2 des cellules adjacentes aux cellules où il y a des forêts est de : "
mean(subset(df2All[c("SPAM2","celladjagain2")], celladjagain2==1)$SPAM2, na.rm = T)
df2All %>% count(celladjagain2)
"La moyenne SPAM2 de l'ensemble des cellules du kenya est de : "
mean(df2All$SPAM2, na.rm = T) 

#Différence de production agricole potentielle (en $) entre les pixels proches
#des forêts et l'ensemble des pixels : 
mean(subset(df2All[c("SPAM2","celladjagain2")], celladjagain2==1)$SPAM2, na.rm = T)-
  mean(df2All$SPAM2, na.rm = T) 


### b) - Gain potentiel lié à l'accessibilité ----

dfagriaccess<-subset(df2All[c("SPAM2","celladjagain2",
                              "Accessi")],celladjagain2 == 1)
#dfagriaccess contient les pixels étant adjacents aux premiers pixels où la couv
#arborée est supérieure à un seuil défini 
Niveaudefaccessi <- 0.6 #A partir de la courbe de distribution sans
#valeurs nulles on peut définir 0,6 comme seuil à partir du quel on
#considère que les pixels sont accessibles.
mean(df2All$Accessi)


#-- Relation between accessibility and agricultural yield (grazing and croping)
#Pour les pixel étant adjacents aux pixels de forêt et accessibles:
mean(subset(dfagriaccess[c("SPAM2","Accessi")], Accessi > Niveaudefaccessi)$SPAM2, na.rm = T)
nrow(subset(dfagriaccess[c("SPAM2","Accessi")], Accessi > Niveaudefaccessi))

#Pour les pixels accessibles :
mean(subset(df2All[c("SPAM2","Accessi")], Accessi > Niveaudefaccessi)$SPAM2, na.rm = T)
nrow(subset(df2All[c("SPAM2","Accessi")], Accessi > Niveaudefaccessi))
#Pour tous les pixels du kenya :
mean(df2All$SPAM2, na.rm = T)

#Régression pour observer l'impact d'être une cellule adjacente sur les 
#rendements
summary(lm(df2All$SPAM2 ~ df2All$celladjagain2))

#Nettoyage si nécessaire : 
remove(dfagriaccess)
remove(dfgain1rast, Dfadjacent,numcellrast,Rastergain1,
       seuilforestgain1,dfgain1)
df2All <- df2All[ , !(names(df2All) %in% c("celladjagain2","cellgain1",
                                           "celladjagain1"))]

#Information sur la distribution de l'accessibilité des pixels :
mean(df2All$Accessi, na.rm = T)
ggplot(df2All, aes(x=Accessi)) + geom_density(fill="lightblue") +
  geom_vline(aes(xintercept=mean(Accessi)),
              color="blue", linetype="dashed", size=1)+labs(y = "Distribution",
                                                            x = "Level of Accessibility of land")+
  ggtitle(label ="Distribution of the level of Accessibility")
#distribution de l'accessibilité en enlevant les valeurs nulles : 
dfaccessi1 <- subset(df2All[c("SPAM2","Accessi")], Accessi > 0.00001)
ggplot(dfaccessi1, aes(x=Accessi)) + geom_density(fill="lightblue") +
  geom_vline(aes(xintercept=mean(Accessi)),
             color="blue", linetype="dashed", size=1)+labs(y = "Distribution",
                                                           x = "Level of Accessibility of land (1 = Accessible)")+
  ggtitle(label ="Distribution of the level of Accessibility")

#Carte de la donnée d'accessibilité des territoires au kenya : 
dftest <- subset(df2All[c("x","y","Accessi")])
nam <- names(dftest[3])
coordinates(dftest)=~x+y
proj4string(dftest) <- proj4string(Rastdecoupe)
writeOGR(dftest, dsn = 'coucheshpR', layer = nam, driver = "ESRI Shapefile")

#Nettoyage si nécessaire :
remove(dfaccess1, dfaccessi1)



### c) - Gain de précipitation forêt ----
"Objectif : Montrer que les forêts augmentent le niveau de 
précipitation de la zone est dans quelle mesure."
#((Certaines régressions doivent utiliser la partie V - 2 - A))


#Regression pour voir l'impact des forêts sur les précipitation
#Avec la valeur des biomes en dummies
"1) - Avoir seulement les valeurs uniques"
datafrunique <- unique(df2All$Biometype)
datafrunique <- as.data.frame(datafrunique)
"2) - Avoir en colonne toutes les valeurs"
dfregpreci <- subset(df2All[c("WCprecimeanall","GFW2010","Biometype")])
for(i in 2:length(datafrunique$datafrunique)){
  namtouy <- paste("Biome",datafrunique$datafrunique[i], sep = "")
  dfregpreci$new <- ifelse(df2All$Biometype==datafrunique$datafrunique[i],
                           datafrunique$datafrunique[i],NA)
  names(dfregpreci)[names(dfregpreci) == "new"] <- namtouy
}#ici le début de la boucle est 2 car le 1 = NA
#La fonction ci-dessus permet d'assigner à chaque pixel ayant une couverture arborée
#et un niveau de précipitation, un type de biome spécifique et d'avoir dans 
#La colonne correspondant à ce numéro de biome spécifique les valeurs 1 lorsqu'il
#Y a ce type de biome dans le pixel et 0 sinon pour pouvoir ensuite utiliser cette colonne
#pour une régression avec des dummies. 
"3) - Ajouter une valeur 1 ou 0 "
fun1 <- function(x){ifelse(is.na(x)==FALSE,x <- 1,x <- 0)}
dfregpreci[4:262] <- lapply(dfregpreci[4:262], fun1)
#Ici [4:262] permet de sélectionner les colonnes où on veut qu'il y ait 0 ou 1 
#donc les colonnes des biometypes pour les utiliser comme dummies.

#### 2) Régressions ----
"4) - Régression "
y <- dfregpreci$WCprecimeanall
dfregpreci <- dfregpreci[ , !(names(dfregpreci) %in% c("WCprecimeanall","Biometype","Biome2"))]#suppression du dummy de ref
#et des colonnes inutiles 
Regpreci1<-lm(y ~.,dfregpreci)
summary(Regpreci1)
#Ici le dummy de comparaison est Biome2, donc la comparaison avec Biome2 doit
#être prise en compte pour l'interprétation des résultats.
coefregpreci1 <- summary(Regpreci1)$coefficients[2,1]#Donne le coefficient de GFW
#dans la régression
stargazer(Regpreci1, type = "text", out = "Regpreci1_R_code1")


#Vérification de la multicolinéarité
vif(Regpreci1)
remove(dfregpreci)

#Régression quantile pour voir la relation entre la précipitation et le niveau de forêt : 
# Mettre tau = 0.5 pour une régression à la médiane, tu peux aussi rentrer une matrice si tu veux plusieurs quantiles
PrecipMod <-  rq(WCprecimeanall ~ GFW2010,tau=0.5 ,data=df2All) 
summary(PrecipMod)
"p-value = 0 donc c'est significatif"
t.test(PrecipMod)
remove(PrecipMod)


#Régression plus simple entre WCprecimean et GFW2010 en considérant qu'il n'y a pas
#beaucoup d'auto-corrélation entre WCprecimean et GFW
Regpreci2 <- lm(df2All$WCprecimeanall ~ df2All$GFW2010)
summary(Regpreci2)
coefregpreci2 <- (summary(Regpreci2)$coefficients[2,1])


#Régression pour connaître l'impact des précipitations sur le PIB du secteur agricole
dfPIBagriandPrecip <- read_excel("Data/Data_PIBprecip/PIBetprecip.xlsx")
Regpreci3 <- lm(dfPIBagriandPrecip$PIBsecteuragricole_Wbank ~ dfPIBagriandPrecip$Precipmean_Wbank)
summary(Regpreci3)
coefregpreci3 <- (summary(Regpreci3)$coefficients[2,1])


#Regression pour voir la relation entre la production agricole en valeur
# et ma précipitation 
Regagri <- lm(df2All$SPAMProdvalha2010 ~ df2All$WCprecimeanall)
summary(Regagri)
coefregagri1 <- summary(Regagri)$coefficient[2,1]
remove(Regagri)
#Pour ajouter plus de puissance statistique : 
y <- df2All$SPAM2
"Here we use dfregpreci because it has all Biome variables"
dfregpreci$WCprecimeanall <- df2All$WCprecimeanall
dfregpreci <- dfregpreci[ , !(names(dfregpreci) %in% c("GFW2010"))]
dfregpreci <- dfregpreci %>% relocate(WCprecimeanall)
Regagri2<-lm(y ~.,dfregpreci)
summary(Regagri2)
#Ici le dummy de comparaison est Biome2, donc la comparaison avec Biome2 doit
#être prise en compte pour l'interprétation des résultats.
coefregagri2 <- summary(Regagri2)$coefficients[2,1]
#Cleaning 
remove(namtouy, datafrunique,my_data,dfregpreci,fun1,y,Regpreci1)

#Nettoyage
remove(dfdatahere,dfchoose,Regagri,Regagri2,Regpreci1,Regpreci2)









#- - - - - - - - - - - - - - - - - - - - - 
# V Informations complémentaires ----
#- - - - - - - - - - - - - - - - - - - - -

## 1 - Graphiques & informations ----

# _a) -- Distribution
#Cropland
mean(df2All$croplandH2016, na.rm = T)
Graph1cropland1 <- ggplot(df2All, aes(x=croplandH2016)) +
  geom_histogram(binwidth=.5, colour="black", fill="white") +
  geom_vline(aes(xintercept=mean(croplandH2016, na.rm=T)),  #le geom_vline permet de mettre ligne de moyenne
             color="red", linetype="dashed", linewidth=1) + 
  ggtitle (label = "Distribution de la surface cultivée par pixel en km²")
Graph1cropland2 <- ggplot(subset(df2All[c(1,2,4)],croplandH2016 > 0), aes(x=croplandH2016)) +
  geom_histogram(binwidth=.5, colour="black", fill="white") +
  geom_vline(aes(xintercept=mean(croplandH2016, na.rm=T)),  #le geom_vline permet de mettre ligne de moyenne
             color="red", linetype="dashed", linewidth=1) + 
  ggtitle (label = "Distribution de la surface cultivée par pixel en km² en enlevant les valeur 0")
#Grazing
mean(df2All$grazingH2016, na.rm = T)
Graph2Grazing1 <- ggplot(df2All, aes(x=grazingH2016)) +
  geom_histogram(binwidth=.5, colour="black", fill="white") +
  geom_vline(aes(xintercept=mean(grazingH2016, na.rm=T)),  #le geom_vline permet de mettre ligne de moyenne
             color="red", linetype="dashed", linewidth=1) + 
  ggtitle (label = "Distribution de la surface utilisée pour le bétail par pixel en km²")
Graph2Grazing2 <- ggplot(subset(df2All[c(1,2,5)],grazingH2016 > 0), aes(x=grazingH2016)) +
  geom_histogram(binwidth=.5, colour="black", fill="white") +
  geom_vline(aes(xintercept=mean(grazingH2016, na.rm=T)),  #le geom_vline permet de mettre ligne de moyenne
             color="red", linetype="dashed", linewidth=1) + 
  ggtitle (label = "Distribution de la surface utilisée pour le bétail par pixel en km² en enlevant les 0")
#Pasture
mean(df2All$pastureH2016, na.rm = T)
Dfcropland <- subset(df2All[c(1,2,10)], pastureH2016 > 0)
Graph2Pasture <- ggplot(Dfcropland, aes(x=pastureH2016)) +
  geom_histogram(binwidth=.5, colour="black", fill="white") +
  geom_vline(aes(xintercept=mean(pastureH2016, na.rm=T)),  #le geom_vline permet de mettre ligne de moyenne
             color="red", linetype="dashed", linewidth=1) + 
  ggtitle (label = "Distribution de la surface utilisée pour les pâturages du bétail \npar pixel en km² en enlevant les 0")
#Rangeland/parcours
mean(df2All$rangelandH2016, na.rm = T)
Dfcropland <- subset(df2All[c(1,2,11)], rangelandH2016 > 0)
Graph2rangeland <- ggplot(Dfcropland, aes(x=rangelandH2016)) +
  geom_histogram(binwidth=.5, colour="black", fill="white") +
  geom_vline(aes(xintercept=mean(rangelandH2016, na.rm=T)),  #le geom_vline permet de mettre ligne de moyenne
             color="red", linetype="dashed", linewidth=1) + 
  ggtitle (label = "Distribution de la surface utilisée pour les parcours du bétail \npar pixel en km² en enlevant les 0")
#Surface batie
mean(df2All$uoppH2016, na.rm = T)
median(df2All$uoppH2016, na.rm = T)
quantile(subset(df2All[c("uoppH2016")]), probs=seq(0, 1, 0.1), na.rm = T)
quantile(subset(df2All[c("uoppH2016")],uoppH2016>0), probs=seq(0, 1, 0.1), na.rm = T)

Graph3Bati1 <- ggplot(df2All, aes(x=uoppH2016)) +
  geom_histogram(binwidth=.5, colour="black", fill="white") +
  geom_vline(aes(xintercept=mean(uoppH2016, na.rm=T)),  #le geom_vline permet de mettre ligne de moyenne
             color="red", linetype="dashed", linewidth=1) + xlab("Urbanized area (km²)") +
  ylab("Frequency") +
  ggtitle (label = "Distribution of the urbanized area in km² per pixel (with 0 values)")
plot(Graph3Bati1)

Graph3Bati2 <- ggplot(subset(df2All[c("uoppH2016")],uoppH2016 > 0), aes(x=uoppH2016)) +
  geom_histogram(binwidth=.5, colour="black", fill="white") +
  geom_vline(aes(xintercept=mean(uoppH2016, na.rm=T)),  #le geom_vline permet de mettre ligne de moyenne
             color="red", linetype="dashed", linewidth=1) + xlab("Urbanized area (km²)") +
  ylab("Frequency") +
  ggtitle (label = "Distribution of the urbanized area in km² per pixel (without 0 values)")
plot(Graph3Bati2)
remove(Graph3Bati1,Graph3Bati2)

#Valeurdeprod SPAM
mean(df2All$SPAMProdvalha2010, na.rm = T)
median(df2All$SPAMProdvalha2010, na.rm = T)
Graph4Valprod1 <- ggplot(df2All, aes(x=SPAMProdvalha2010, na.rm = T)) +
  geom_histogram(binwidth=.5, colour="black", fill="white") +
  geom_vline(aes(xintercept=mean(SPAMProdvalha2010, na.rm=T)),  
             color="red", linetype="dashed", linewidth=1) + 
  ggtitle (label = "Distribution de la valeur de production par hectare par pixel")
plot(Graph4Valprod1)
#SPAM2 - Agricultural and grazing/livestock yields 
mean(df2All$SPAM2, na.rm = T)
median(df2All$SPAM2, na.rm = T)
Graph4Valprod3 <- ggplot(df2All, aes(x=SPAM2, na.rm = T)) +
  geom_histogram(binwidth=.5, colour="black", fill="white") +
  geom_vline(aes(xintercept=mean(SPAM2, na.rm=T)),  
             color="red", linetype="dashed", linewidth=1) + 
  ggtitle (label = "Distribution de la valeur de production (culture et élevage) par hectare par pixel")
plot(Graph4Valprod1)

#Graphique de la relation entre précipitation moyenne et pourcentage de couverture arborée
plot(df2All$WCprecimeanall, df2All$GFW2010, xlab = "Précipitations moyennes mensuelles de 2021 (en mm)"
     ,ylab = "Pourcentage de couverture arborée", main = "Relation entre les précipitations moyennes mensuelles \n et la couverture arborée au Kenya")

#Précipitation
mean(df2All$WCprecimeanall, na.rm = T)
median(df2All$WCprecimeanall, na.rm = T)
mediantreshforest1 <- median(df2All$WCprecimeanall, na.rm = T)
quantile(subset(df2All[c("WCprecimeanall")]), probs=seq(0, 1, 0.1), na.rm = T)
length(subset(df2All[c("WCprecimeanall")],WCprecimeanall<mediantreshforest1)$WCprecimeanall)
sum(is.na(df2All$WCprecimeanall))
Graphpreci1 <- ggplot(df2All, aes(x=WCprecimeanall)) +
  geom_histogram(binwidth=.5, colour="black", fill="white") +
  geom_vline(aes(xintercept=mean(WCprecimeanall, na.rm=T)),  #le geom_vline permet de mettre ligne de moyenne
             color="red", linetype="dashed", linewidth=1) + xlab("Average monthly precipitation (mm)") +
  ylab("Frequency") +
  ggtitle (label = "Distribution of average monthly precipitation in mm per pixel")
plot(Graphpreci1)
remove(Graphpreci1,mediantreshforest1)
#Temperature
mean(df2All$WCtmaxmean2021, na.rm = T)
median(df2All$WCtmaxmean2021, na.rm = T)
mediantreshforest2 <- median(df2All$WCtmaxmean2021, na.rm = T)
quantile(subset(df2All[c("WCtmaxmean2021")]), probs=seq(0, 1, 0.1), na.rm = T)
length(subset(df2All[c("WCtmaxmean2021")],WCtmaxmean2021<mediantreshforest2)$WCtmaxmean2021)
sum(is.na(df2All$WCtmaxmean2021))
Graphpreci1 <- ggplot(df2All, aes(x=WCtmaxmean2021)) +
  geom_histogram(binwidth=.5, colour="black", fill="white") +
  geom_vline(aes(xintercept=mean(WCtmaxmean2021, na.rm=T)),  #le geom_vline permet de mettre ligne de moyenne
             color="red", linetype="dashed", linewidth=1) + xlab("Average monthly maximum temperatures (°C)") +
  ylab("Frequency") +
  ggtitle (label = "Distribution of average monthly maximum temperatures in degrees per pixel")
plot(Graphpreci1)
remove(Graphpreci1,mediantreshforest2)
#Distribution de Forestkm²
mean(df2All$Forestkm2, na.rm = T)
median(df2All$Forestkm2, na.rm = T)
mediantreshforest3 <- median(df2All$Forestkm2, na.rm = T)
quantile(subset(df2All[c("Forestkm2")]), probs=seq(0, 1, 0.1), na.rm = T)
length(subset(df2All[c("Forestkm2")],Forestkm2<mediantreshforest3)$Forestkm2)
sum(is.na(df2All$Forestkm2))
Graphpreci1 <- ggplot(df2All, aes(x=Forestkm2)) +
  geom_histogram(binwidth=.5, colour="black", fill="white") +
  geom_vline(aes(xintercept=mean(Forestkm2, na.rm=T)),  #le geom_vline permet de mettre ligne de moyenne
             color="red", linetype="dashed", linewidth=1) + xlab("Average monthly maximum temperatures (°C)") +
  ylab("Frequency") +
  ggtitle (label = "Distribution of tree cover in km² based on Global Forest Watch and HYDE data per pixel")
plot(Graphpreci1)
remove(Graphpreci1,mediantreshforest3)

#Dsitribution de la couverture arborée des études
"Pour tout les pixels : "
mean(df2All$GFW2010)
median(df2All$GFW2010)
quantile(df2All[c("GFW2010")], probs=seq(0, 1, 0.1), na.rm = T)
GraphGFW1 <- ggplot(df2All, aes(x=GFW2010)) +
  geom_histogram(binwidth=.5, colour="black", fill="white") +
  geom_vline(aes(xintercept=mean(GFW2010, na.rm=T)),  #le geom_vline permet de mettre ligne de moyenne
             color="red", linetype="dashed", linewidth=1) + xlab("Tree cover (%)") +
  ylab("Frequency") +
  ggtitle (label = "Distribution of the tree cover per pixel in percent")
plot(GraphGFW1)

"Tous les pixels ayant GFW > 0,9999999"
mean(subset(df2All[c("GFW2010")],GFW2010>0.9999)$GFW2010)
median(subset(df2All[c("GFW2010")],GFW2010>0.9999)$GFW2010)
quantile(subset(df2All[c("GFW2010")],GFW2010>0.9999), probs=seq(0, 1, 0.1), na.rm = T)
GraphGFW2 <- ggplot(subset(df2All[c("GFW2010")],GFW2010>0.9999), aes(x=GFW2010)) +
  geom_histogram(binwidth=.5, colour="black", fill="white") +
  geom_vline(aes(xintercept=mean(GFW2010, na.rm=T)),  #le geom_vline permet de mettre ligne de moyenne
             color="red", linetype="dashed", linewidth=1) + xlab("Tree cover (%)") +
  ylab("Frequency") +
  ggtitle (label = "Distribution of the tree cover per pixel in percent (without [0;1[ values)")
plot(GraphGFW2)


# -- Boxplot de distribution
boxplot1 <- boxplot(df2All[c(4,5,17)],na.rm = T, 
                    main = "Boxplot des surfaces utilisées pour la culture,\n l'élevage et les constructions en termes de km²\n par pixel en gardant les 0",names = c("crop","graz","baties"),
                    las = 2,col = c("orange","red"),border = "brown",
                    notch = TRUE)
boxplot2crop <- boxplot(subset(df2All[c("croplandH2016")],croplandH2016>0),na.rm = T, 
                        main = "Boxplot des surfaces cultivées en km²\n par pixel (sans les 0)",
                        las = 2,col = c("orange","red"),border = "brown",
                        notch = TRUE)
boxplot3graz <- boxplot(subset(df2All[c("grazingH2016")],grazingH2016>0),na.rm = T, 
                        main = "Boxplot des surfaces utilisées pour l'élevage,\n en km² par pixel (sans les 0)",
                        las = 2,col = c("orange","red"),border = "brown",
                        notch = TRUE)
boxplot4bati <- boxplot(subset(df2All[c("uoppH2016")],uoppH2016>0),na.rm = T, 
                        main = "Boxplot des surfaces utilisées pour les constructions,\n en km² par pixel (sans les 0)",
                        las = 2,col = c("orange","red"),border = "brown",
                        notch = TRUE)
remove(boxplot1,boxplot2crop,boxplot3graz,boxplot4bati,dfboxplot)


# -- Distribution de la séquestration carbone - Cook Patton et al. 2020
ggplot(df2All, aes(x=regrowthcookpatton)) +
  geom_histogram(binwidth=.5, colour="black", fill="white") +
  geom_vline(aes(xintercept=mean(regrowthcookpatton, na.rm=T)), 
             color="red", linetype="dashed", linewidth=1) + 
  ggtitle (label = "Distribution des niveaux de séquestration carbone 
           en laissant la forêt naturelle repousser (en Mg carbone/ha/year)") +
  labs(y= "Nombre d'observations", x = "Niveau Mg carbone/ha/year")

#Analyse de la data de Cook Patton et al - Seq Carbone
#La somme totale de carbon pouvant être séquestrée par le Kenya par an : 
df2All$sequestrbyyr <- df2All$maxlnH2016*100*df2All$regrowthcookpatton 
#La colonne sequestrbyyr donne le montant total de carbone pouvant être séquestré
#par an pour chaque pixel en Mg (ou tonne) de carbone.
sum(df2All$sequestrbyyr, na.rm = T)#donne la quantité totale de carbone
#pouvant être séquestrée par le Kenya en laissant pousser la forêt naturelle
df2All$sequestrbyyrinco2 <- df2All$sequestrbyyr*3.67 #according to https://www.epa.gov/energy/greenhouse-gases-equivalencies-calculator-calculations-and-references
#Colonne sequestrbyyrinco2 donne la quantité de CO2 pouvant être absorbé par pixel par an

#Pour déterminer quels sont les espaces où il peut y avoir des forêts naturelles
#au Kenya si on veut savoir avec précision quels sont les pixels où il y a
#un certain niveau de séquestration carbone.
"Niveau de séquestration d'une forêt naturelle :"
seuilfornatcookpatt <- 1.955 #1.955 représente la moyenne ici pour le Kenya de
#la seq carbone donnée par la colonne regrowthcookpatton, ce seuil peut être modifié
#pour obtenir d'autres informtions
for(i in 1:length(df2All$GFW2010)){
  df2All$forcookpattype[i] <- ifelse(is.na(df2All$regrowthcookpatton[i])==FALSE,
                                     ifelse(df2All$regrowthcookpatton[i]>=seuilfornatcookpatt,
                                            1,0),NA)
}
#forcookpattype donne si le niveau de séquestration du carbone est haut dessus
#du seuil défini et donc peut être considéré comme une forêt, donc si forcookpattype = 1 le pixel
#peut-être considéré comme un endroit où il peut y avoir une forêt naturelle importante
df2All <- df2All[ , !(names(df2All) %in% c("forcookpattype"))]




## 2 - Informations précipitations ----

#Information générale :
min(df2All$WCprecimeanall)
max(df2All$WCprecimeanall)
length(subset(df2All[c("x","y","WCprecimeanall")], WCprecimeanall < 90))
length(subset(df2All[c("x","y","WCprecimeanall")], WCprecimeanall < 100))
length(subset(df2All[c("x","y","WCprecimeanall")], WCprecimeanall > 131))

#Information décile
df2All$Classprecipi <- ifelse(df2All$WCprecimeanall < quantile(df2All[c("WCprecimeanall")], probs=seq(0, 1, 0.1), na.rm = T)[[2]],1,
                              ifelse(df2All$WCprecimeanall < quantile(df2All[c("WCprecimeanall")], probs=seq(0, 1, 0.1), na.rm = T)[[3]],2,
                                     ifelse(df2All$WCprecimeanall < quantile(df2All[c("WCprecimeanall")], probs=seq(0, 1, 0.1), na.rm = T)[[4]],3,
                                            ifelse(df2All$WCprecimeanall < quantile(df2All[c("WCprecimeanall")], probs=seq(0, 1, 0.1), na.rm = T)[[5]],4,
                                                   ifelse(df2All$WCprecimeanall < quantile(df2All[c("WCprecimeanall")], probs=seq(0, 1, 0.1), na.rm = T)[[6]],5,
                                                          ifelse(df2All$WCprecimeanall < quantile(df2All[c("WCprecimeanall")], probs=seq(0, 1, 0.1), na.rm = T)[[7]],6,
                                                                 ifelse(df2All$WCprecimeanall < quantile(df2All[c("WCprecimeanall")], probs=seq(0, 1, 0.1), na.rm = T)[[8]],7,
                                                                        ifelse(df2All$WCprecimeanall < quantile(df2All[c("WCprecimeanall")], probs=seq(0, 1, 0.1), na.rm = T)[[9]],8,
                                                                               ifelse(df2All$WCprecimeanall < quantile(df2All[c("WCprecimeanall")], probs=seq(0, 1, 0.1), na.rm = T)[[10]],9,
                                                                                      ifelse(df2All$WCprecimeanall <= quantile(df2All[c("WCprecimeanall")], probs=seq(0, 1, 0.1), na.rm = T)[[11]],10,0))))))))))

#Information sur les déciles de la précipitation moyenne
mean(df2All$WCprecimeanall)
quantile(df2All[c("WCprecimeanall")], probs=seq(0, 1, 0.1), na.rm = T)
quantile(subset(df2All, WCprecimeanall > 0.99999)[c("WCprecimeanall")], probs=seq(0, 1, 0.1), na.rm = T)
#Nettoyage : 
df2All <- df2All[ , !(names(df2All) %in% c("Classprecipi"))]


#Indication sur le territoire du Kenya en fonction de la précipitation moyenne de ce territoire
for(i in 1:10){
  print("Nombre de km² de surface par décile et indiquant le nombre de pixel:")
  print("Précipitation moyenne mensuelle < à :")
  print(quantile(df2All[c("WCprecimeanall")], probs=seq(0, 1, 0.1), na.rm = T)[[i+1]])
  print("Nombre de km² de surface :")
  print(sum(subset(df2All, Classprecipi == i)$maxlnH2016,na.rm = T))
  print("Nombre de pixel relatifs :")
  print(nrow(subset(df2All[c("Classprecipi")], Classprecipi == i)))
  print("")
  print("--NEXT--")
}

# -- précipitations des forêts
DFforest1<-subset(df2All, GFW2010 >= quantile(df2All[c("GFW2010")],probs=seq(0, 1, 0.1), na.rm = T)[[9]])
for(i in 1:10){
  print("-__Pour les territoires avec une couverture arborée > 8eme décile de la couverture arborée")
  print("Nombre de km² de surface par décile et indiquant le nombre de pixel:")
  print("Précipitation moyenne mensuelle < à :")
  print(quantile(DFforest1[c("WCprecimeanall")], probs=seq(0, 1, 0.1), na.rm = T)[[i+1]])
  print("Nombre de km² de surface :")
  print(sum(subset(DFforest1, Classprecipi == i)$maxlnH2016,na.rm = T))
  print("Nombre de pixel relatifs :")
  print(nrow(subset(DFforest1[c("Classprecipi")], Classprecipi == i)))
  print("")
  print("--NEXT--")
}
remove(DFforest1)


#-_-_-_-_-_-_-_-__-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-__-_-_-_-_-_-
#Relation entre précipitation et GFW :
#__--__--___--_-_-_-__-_-_-_-_-__-_-_-__-_-___-_-_-_-___-_-_-_-_

mean(df2All$WCprecimeanall, na.rm = T)
median(df2All$WCprecimeanall,na.rm = T)
#graphique de relation de base entre précipitation et forêt
plot1precip<-plot(df2All$WCprecimeanall, df2All$GFW2010,
     main="Relation between tree cover (%) and precipitations (mm)",
     xlab="X axis title",
     ylab="Y axis title")
plot(df2All$WCprecimeanall, df2All$GFW2010,
     main="Relation between tree cover (%) and precipitations (mm)",
     xlab="X axis title",
     ylab="Y axis title")
#graphique de relation entre précipitation et forêt en enlevant 
#les pixels avec couverture arborée = 0 sont retirés
dfforest2<-subset(df2All, GFW2010 > 0)
plot(dfforest2$WCprecimeanall, dfforest2$GFW2010)
remove(dfforest2)
#graphique de relation entre précipitation et forêt en enlevant 
#les pixels avec couverture arborée <3 sont retirés
dfforest3<-subset(df2All, GFW2010 > 3)
plot(dfforest3$WCprecimeanall, dfforest3$GFW2010)
remove(dfforest3)

#Calcul représentatif
#Pour des pixels avec forêt supérieur à 47%
dfforest4<-subset(df2All, GFW2010 > 10)
dfforest4<-subset(dfforest4, WCprecimeanall > 57)
#la moyenne des précipitations de ces pixels est de 
mean(dfforest4$WCprecimeanall, na.rm = T)
plot(dfforest4$WCprecimeanall, dfforest4$GFW2010,
     main="Relation between tree cover (%) and precipitations (mm)
     (without values < 10 to better see the interaction)",
     xlab="Precipitations (mm)",
     ylab="Tree cover (%)")
plotprecip2<-plot(dfforest4$WCprecimeanall, dfforest4$GFW2010,
     main="Relation between tree cover (%) and precipitations (mm)
     (without values < 10 to better see the interaction)",
     xlab="Precipitations (mm)",
     ylab="Tree cover (%)")
remove(dfforest4)
#Pour connaître le nombre de pixel ayant couverture arborée >= 47% 
# et ayant précipitations > 62 mm
df2All$precipnumber2 <- df2All$Condtnbatprecitemp
for(i in 1:length(df2All$x)){
  if(df2All$GFW2010[i] >= 47.86 &
     df2All$WCprecimeanall[i] >= 62){
  df2All$precipnumber2[i]<-1}
  else{df2All$precipnumber2[i]<-0}
}
dftestprecipita3 <- subset(df2All[c("GFW2010","WCprecimeanall","precipnumber2")],
                           precipnumber2 == 1)
df2All <- df2All[ , !(names(df2All) %in% c("precipnumber2"))]
remove(dftestprecipita3)

#Pour connaître le nombre de pixel relatif à une précipitation comprise entre 2 niveaux 
df2All$precipnumber4 <- df2All$Condtnbatprecitemp
for(i in 1:length(df2All$x)){
  if(df2All$WCprecimeanall[i] < 68 &
     df2All$WCprecimeanall[i] >= 62.13){
    df2All$precipnumber4[i]<-1}
  else{df2All$precipnumber4[i]<-0}
}
sum(df2All$precipnumber4, na.rm = T)

df2All <- df2All[ , !(names(df2All) %in% c("precipnumber4"))]

#Pour connaître les pixels compris entre deux seuils de précip et avec une couverture
#arborée > à un certain seuil :
df2All$precipnumber5 <- df2All$Condtnbatprecitemp
for(i in 1:length(df2All$x)){
  if(df2All$GFW2010[i] >= 51 &
     df2All$WCprecimeanall[i] >= 75 &
     df2All$WCprecimeanall[i] < 100){
    df2All$precipnumber5[i]<-1}
  else{df2All$precipnumber5[i]<-0}
}
dftestprecipita12 <- subset(df2All[c("GFW2010","WCprecimeanall","precipnumber5")],
                           precipnumber5 == 1)
df2All <- df2All[ , !(names(df2All) %in% c("precipnumber5"))]
remove(dftestprecipita12)


#Pour des précipitation entre 51 et 52mm (représentant la mediane de la précipitation
#moyenne mensuelle pour l'année 2021 au kenya)
dfforest5<-subset(df2All, WCprecimeanall>51)
dfforest5<-subset(df2All, WCprecimeanall<52)
#la moyenne GFW de ces précipitations est :
mean(dfforest5$GFW2010,na.rm = T)
#le max GFW de ces précipitations est 
max(dfforest5$GFW2010, na.rm = T)
remove(dfforest5)
#Pour les pixels ayant une couverture de forêt > à 5%
dfforest6<-subset(df2All, GFW2010 > 5)
quantile(dfforest6[c("WCprecimeanall")], probs=seq(0, 1, 0.1), na.rm = T)
remove(dfforest6)
#Courbe de distribution des précipitation pour les pixels
#où la couverture arborée est supérieure au 1er décile de couverture (en enlevant
# les valeurs nulles de GFW) (Condition 1k2)
dfforest7 <- subset(df2All, GFW2010 >  0)
quantile(dfforest7[c("GFW2010")], probs=seq(0, 1, 0.1), na.rm = T)
premdec <- quantile(dfforest7[c("GFW2010")], probs=seq(0, 1, 0.1), na.rm = T)[2]
dfforest8 <- subset(df2All, GFW2010 >  premdec)
ggplot(dfforest8, aes(x=WCprecimeanall)) +
  geom_histogram(binwidth=.5, colour="black", fill="white") +
  geom_vline(aes(xintercept=mean(WCprecimeanall, na.rm=T)),  #le geom_vline permet de mettre ligne de moyenne
             color="red", linetype="dashed", linewidth=1) + 
  ggtitle (label = "Distribution de la précipitation selon condition 1k2 en mm")


#Boxplot de la répartition de la précipitation en fonction de GFW
dfpreciind <- subset(df2All[c("GFW2010","WCprecimeanall")])
for(i in 1:length(dfpreciind$GFW2010)){
  dfpreciind$indice[i] <- ifelse(dfpreciind$GFW2010[i]<=25,
                             "[0;25]",ifelse(dfpreciind$GFW2010[i]>75,"]75;max]",0))
  dfpreciind$indice[i] <- ifelse(dfpreciind$indice[i] == 0,
                             ifelse(dfpreciind$GFW2010[i]<=50,"]25;50]","]50;75]"),dfpreciind$indice[i])
}

#Faire un boxplot :
boxplot(WCprecimeanall~indice,data=dfpreciind,
        xlab = "Level of forest (in % of space)*",
        col =(c("azure","cadetblue1","cornflowerblue","blue1")),
        main="Precipitation regarding the level of forest",
        ylab = "Level of precipitation (in mm)")
#Moyenne, médiane par groupe : 
# - Pour le premier quantile
mean(subset(dfpreciind[c("WCprecimeanall","GFW2010")], dfpreciind$indice == "[0;25]")$WCprecimeanall, na.rm = T)
median(subset(dfpreciind[c("WCprecimeanall","GFW2010")], dfpreciind$indice == "[0;25]")$WCprecimeanall, na.rm = T)
# - Pour le deuxième quantile
mean(subset(dfpreciind[c("WCprecimeanall","GFW2010")], dfpreciind$indice == "]25;50]")$WCprecimeanall, na.rm = T)
median(subset(dfpreciind[c("WCprecimeanall","GFW2010")], dfpreciind$indice == "]25;50]")$WCprecimeanall, na.rm = T)
# - Pour le troisième quantile
mean(subset(dfpreciind[c("WCprecimeanall","GFW2010")], dfpreciind$indice == "]50;75]")$WCprecimeanall, na.rm = T)
median(subset(dfpreciind[c("WCprecimeanall","GFW2010")], dfpreciind$indice == "]50;75]")$WCprecimeanall, na.rm = T)
# - Pour le quatrième quantile
mean(subset(dfpreciind[c("WCprecimeanall","GFW2010")], dfpreciind$indice == "]75;max]")$WCprecimeanall, na.rm = T)
median(subset(dfpreciind[c("WCprecimeanall","GFW2010")], dfpreciind$indice == "]75;max]")$WCprecimeanall, na.rm = T)

remove(DF4_precifor,DF3_precifor,DF2_precifor,DF1_precifor)
remove(dfpreciind,combined)


## 3 - Informations accessibilité ----

#Number of NA in accessibility data (Market Access)
sum(is.na(df2All$Accessi))
#For the entire data of accessibility:
#Mean and median :
mean(df2All$Accessi,na.rm = T)
median(df2All$Accessi,na.rm = T)
quantile(subset(df2All[c("Accessi")],Accessi>0.0001), probs=seq(0, 1, 0.1), na.rm = T)

#For data of accessibility without value = 0: 
mean(subset(df2All[c("Accessi")],Accessi>0.0001)$Accessi)
median(subset(df2All[c("Accessi")],Accessi>0.0001)$Accessi)
quantile(subset(df2All[c("Accessi")],Accessi>0.0001), probs=seq(0, 1, 0.1), na.rm = T)

#Boxplot de la répartition de l'accessibilité en fonction de SPAM
dfaccess1 <- subset(df2All[c("Accessi","SPAM2")])
seuil50access <- quantile(subset(df2All[c("Accessi")]), probs=seq(0, 1, 0.01), na.rm = T)[[51]]
seuil75access <- quantile(subset(df2All[c("Accessi")]), probs=seq(0, 1, 0.01), na.rm = T)[[76]]
for(i in 1:length(dfaccess1$Accessi)){
  dfaccess1$indice[i] <- ifelse(dfaccess1$Accessi[i]<=seuil50access,
                                 "[0 - Médiane]",ifelse(dfaccess1$Accessi[i]<=seuil75access,"]Médiane - 3e Quartile]","]3e Quartile - 4e Quartile]"))
}
#Faire un boxplot :
boxplot(SPAM2~indice,data=dfaccess1,
        xlab = "Quartile repartition*",
        col =(c("cornsilk","bisque3","cornsilk4")),
        main="Level of agriculture production regarding the level of accessibility of land",
        ylab = "Value of agriculture production ($)")
#Moyenne de SPAM pour les valeur inférieur à la médiane
mean(subset(dfaccess1[c("indice","SPAM2")], indice = "[0 - Médiane]")$SPAM2, na.rm = T)
#Moyenne de SPAM pour les valeur entre la médiane et 3eme quartile
mean(subset(dfaccess1[c("indice","SPAM2")], indice = "]Médiane - 3e Quartile]")$SPAM2, na.rm = T)
#Moyenne de SPAM pour les valeur entre le 3eme quartile et 4eme quartile
mean(subset(dfaccess1[c("indice","SPAM2")], indice = "]3e Quartile - 4e Quartile]")$SPAM2, na.rm = T)

remove(dfaccess1,seuil50access,seuil75access)






## 4 - Réalisation de carte ----

#-- Précipitation moyenne :
dfprecimean <- df2All[,c("x","y","WCprecimeanall")]
raspreci <- rasterFromXYZ(dfprecimean)
writeRaster(raspreci, "rastprecimean", format = "GTiff")#le raster sera exporté
remove(raspreci,dfprecimean)

#-- Température maximale moyenne : 
dftempmean <- df2All[,c("x","y","WCtmaxmean2021")]
rastemp <- rasterFromXYZ(dftempmean)
writeRaster(rastemp, "rastempmean", format = "GTiff")#le raster sera exporté
remove(rastemp,dftempmean)

#-- Carte des pixels répondant aux conditions (Surface bâtie,précipitation,
#Température)
dfbatprecitemp <- df2All[,c("x","y","Condtnbatprecitemp")]
rastbatprecitemp <- rasterFromXYZ(dfbatprecitemp)
writeRaster(rastbatprecitemp,"Rastconditions", format = "GTiff")
remove(dfbatprecitemp,rastbatprecitemp)

#-- Carte des pixels répondant aux conditions (Surface bâtie,précipitation,
#Température) mais ayant des valeurs NA pour SPAM
dfSPAM <- df2All[,c("x","y","SPAMProdvalha2010")]
rastSPAM <- rasterFromXYZ(dfSPAM)
writeRaster(rastSPAM,"rastSPAM2", format = "GTiff")
remove(dfSPAM,rastSPAM)
plot(rastSPAM)

#-- Carte des pixels ayant 0 km² de surface agricole-cropland
for(i in 1:length(df2All$GFW2010)){
  df2All$hydetest[i]<-ifelse(df2All$croplandH2016[i] == 0,1,0)
}
dfHYDEcheck <- df2All[,c("x","y","hydetest")]
rastHYDE1 <- rasterFromXYZ(dfHYDEcheck)
writeRaster(rastHYDE1 ,"rastHYDE1", format = "GTiff")
remove(dfHYDEcheck,rastHYDE1)
df2All <- df2All[, !(names(df2All) %in% c("hydetest"))]

#-- Carte des pixels ayant 0 km² de surface agricole-grazing
df2All$hydetest2<-df2All$croplandH2016
for(i in 1:length(df2All$GFW2010)){
  df2All$hydetest2[i]<-ifelse(df2All$grazingH2016[i] == 0,1,0)
}
dfHYDEcheck2 <- df2All[,c("x","y","hydetest2")]
rastHYDEgraz1 <- rasterFromXYZ(dfHYDEcheck2)
writeRaster(rastHYDEgraz1 ,"rastHYDEverifsigrazing", format = "GTiff")
remove(dfHYDEcheck2,rastHYDEgraz1)
df2All <- df2All[, !(names(df2All) %in% c("hydetest2"))]

#-- Carte des pixels ayant NA pour SPAM-valeur de prod
for(i in 1:length(df2All$GFW2010)){
  df2All$SPAMisna[i]<-ifelse(is.na(df2All$SPAMProdvalha2010[i])==TRUE,1,0)
}
dfSPAM2 <- df2All[,c("x","y","SPAMisna")]
rastSPAM2 <- rasterFromXYZ(dfSPAM2)
writeRaster(rastSPAM2 ,"rastSPAMverif", format = "GTiff")
remove(dfSPAM2,rastSPAM2)
df2All <- df2All[, !(names(df2All) %in% c("SPAMisna"))]

#-- Carte des pixels ayant NA pour SPAM-quantité produite (en tonne)
df2All$SPAMisnaquantity<-df2All$Condtnbatprecitemp
for(i in 1:length(df2All$GFW2010)){
  df2All$SPAMisnaquantity[i]<-ifelse(is.na(df2All$SPAMProdquant2010_mtunit[i])==TRUE,1,0)
}
dfSPAMquantityNA <- df2All[,c("x","y","SPAMisnaquantity")]
rastSPAMquantityNA <- rasterFromXYZ(dfSPAMquantityNA)
writeRaster(rastSPAMquantityNA ,"rastSPAMquantiverif", format = "GTiff")
remove(dfSPAMquantityNA,rastSPAMquantityNA)
df2All <- df2All[, !(names(df2All) %in% c("SPAMisnaquantity"))]

#-- Carte des pixels des surfaces bâties :
dfhydebat2 <- df2All[,c("x","y","uoppH2016")]
rasthydebat2 <- rasterFromXYZ(dfhydebat2)
writeRaster(rasthydebat2 ,"rasthydebat2", format = "GTiff")
remove(dfhydebat2,rasthydebat2)

#-- Carte des rendements agricole de SPAM : 
dfspammap <- df2All[,c("x","y","SPAMProdvalha2010")]
rastSPAMmap2 <- rasterFromXYZ(dfspammap)
writeRaster(rastSPAMmap2 ,"rastSPAMrendement$", format = "GTiff")
remove(dfspammap,rastSPAMmap2)

#-- Carte des rendements agricoles SPAM en termes de quantité produite (en tonne) :
dfmapSPAMqtty <- df2All[,c("x","y","SPAMProdquant2010_mtunit")]
coordinates(dfmapSPAMqtty ) <- ~ x + y
proj4string(dfmapSPAMqtty ) <- proj4string(Rastdecoupe)
dfmapSPAMqtty  <- st_as_sf(dfmapSPAMqtty)
#Extract le shapefile
st_write(dfmapSPAMqtty , "SpamProdQuanti.shp")
remove(dfmapSPAMqtty )

#-- Carte de la distribution du km² de forêt par pixel :
dfForestmap <- df2All[,c("x","y","Forestkm2")]
coordinates(dfForestmap) <- ~ x + y
proj4string(dfForestmap) <- proj4string(Rastdecoupe)
dfForestmap  <- st_as_sf(dfForestmap)
#Extract le shapefile
st_write(dfForestmap, "Forestkm2parpixelmap.shp")
remove(dfForestmap)

#-- Carte des rendements agricole de SPAM2 : 
dfspammap2 <- df2All[,c("x","y","SPAM2")]
rastSPAM2map1 <- rasterFromXYZ(dfspammap2)
writeRaster(rastSPAM2map1 ,"rastSPAM2$", format = "GTiff")
remove(dfspammap2,rastSPAM2map1)

#-- Map of land accessibility :
dfaccessibilty <- df2All[,c("x","y","Accessi")]
coordinates(dfaccessibilty) <- ~ x + y
proj4string(dfaccessibilty) <- proj4string(Rastdecoupe)
dfaccessibilty <- st_as_sf(dfaccessibilty)
#Extract le shapefile
st_write(dfaccessibilty, "Accessimaps02.shp")
remove(dfaccessibilty)

#-- Carte de la distribution des SPAM2 = 0 (donc les pixels où SPAM = NA et
#où la surface consacré à l'élevage = 0)
#POur voir la carte SPAM2 et voir où sont situés les pixels avec SPAM2 = 0 et donc
#malheureusement et logiquement les pixels où SPAM2=0 sont présent sur les zones de
#reforestation de l'algorithme 1 étant donné que cleui-ci vise à diminuer la perte
#de production agricole lors de la reforestation
nrow(subset(df2All[c("SPAM2")], SPAM2 == 0))
df1algo1 <- subset(df2All[c("x","y","SPAM2")])
dshp <- df1algo1
coordinates(dshp)=~x+y
proj4string(dshp) <- proj4string(Rastdecoupe)
writeOGR(dshp, dsn = 'coucheshpR', layer = "SPAM2", driver = "ESRI Shapefile")
remove(df1algo1,dshp)

#-- Carte de la séquestration carbone - Cook Patton et al. 2020 :
dfcookpattonverif <- df2All[,c("x","y","regrowthcookpatton")]
rastcookpatton <- rasterFromXYZ(dfcookpattonverif)
writeRaster(rastcookpatton,"rastcookpatton", format = "GTiff")
remove(dfcookpattonverif,rastcookpatton)

#-- Carte des % de couverture arborée par pixel d'après Global Forest Watch:
dfglobalforestw <- df2All[,c("x","y","GFW2010")]
rastGlobalForestW <- rasterFromXYZ(dfglobalforestw)
writeRaster(rastGlobalForestW,"rastglobalforestwatch", format = "GTiff")
remove(dfglobalforestw,rastGlobalForestW)



#- - - - - - - - - - - - - - - - - - - - - 
# VI - Annexes ----
#- - - - - - - - - - - - - - - - - - - - -

#Annexe 1 - Assemblage des rasters : 
#Certains rasters ne représentent qu'une partie d'un territoire,
#dans ce cas ces rasters vont être assemblés pour ensuite suivre les mêmes
#traitement que les autres rasters. Par exemple pour Global Forest Watch il 
#y a 4 rasters et au lieu d'avoir 4 colonnes dans une base on obtient à la fin 1 seule
#colonne (représentant l'entièreté du Kenya) permettant de faire des traitements et analyse avec les autres rasters 
#qui représentent, eux, directement l'entièreté du Kenya.


#Annexe 2 - Obtention de la valeur de production totale de crop au Kenya d'après SPAM  (en $)
rastsomme <- raster("Data/Data_rastSPAMProduction/spam2010V2r0_global_V_agg_VP_CROP_A.tif")
crs(rastsomme) <- crs(Rastdecoupe)
rastsomme <- resample(rastsomme,list_rastreso2[[1]])#Attention au resample
#pouvant altérer la disposition des pixels donc nécessaire de faire une
#vérification en observant le raster initiale et le raster provenant du df
#de R qui a été traité et découpé à la taille du Kenya
rastsomme <- crop(rastsomme,Rastdecoupe)
dfRastsomme <- as.data.frame(rasterToPoints(rastsomme))
coordinates(dfRastsomme)=~x+y
proj4string(dfRastsomme) <- proj4string(Rastdecoupe)
dfRastsomme <- as.data.frame(dfRastsomme[Rastdecoupe,])
names(dfRastsomme)[3] <- "prod"
"La valeur totale de production agricole du Kenya selon SPAM est de ($) :"
sum(dfRastsomme$prod, na.rm = T)
remove(dfRastsomme,rastsomme)

#Sortir df2All :
write_xlsx(df2All,"Data/df2All1.xlsx")
