#Analysis examples:
# https://onlinelibrary.wiley.com/doi/full/10.1111/mec.13426
# https://onlinelibrary.wiley.com/doi/full/10.1111/jeb.12915#support-information-section
# https://onlinelibrary.wiley.com/doi/full/10.1111/j.1558-5646.2008.00425.x 

library(vegan)
#mantel and partial mantel tests, http://cc.oulu.fi/~jarioksa/softhelp/vegan/html/mantel.html
#make dissimilarity matrix: vegdist
#http://rfunctions.blogspot.com/2016/10/simple-and-partial-mantel-tests.html

#READ DATA
setwd("C:\\Users\\Buckley\\Google Drive\\Buckley\\Work\\GrasshopperGenetics\\data\\matrices")

elev= read.csv("elev.csv", row.names="Site")
envi= read.csv("envi.csv", row.names="Site")
geo= read.csv("geo.csv", row.names="Site")

#phenotypes
#phen= read.csv("phen.csv", row.names="Site")
phen_boul= read.csv("phen_boul.csv", row.names="Site")
phen_clav= read.csv("phen_clav.csv", row.names="Site")
phen_pell= read.csv("phen_pell.csv", row.names="Site")
phen_sang= read.csv("phen_sang.csv", row.names="Site")
#drop pell RF
phen_pell= phen_pell[c("a1","b1","c1"),]

#genetics
gen_boul= read.csv("genetics_boul.csv", row.names="X")
gen_clav= read.csv("genetics_clav.csv", row.names="X")
gen_pell= read.csv("genetics_pell.csv", row.names="X")
gen_sang= read.csv("genetics_sang.csv", row.names="X")

#======================================================================
#SET UP SITE MATRICES

#geo
geo_boul= geo[c("a1","b1","c1","d1"),]
geo_clav= geo[c("RF","a1","b1","d1"),]
geo_pell= geo[c("a1","b1","c1"),]
geo_sang= geo[c("RF","a1","b1","c1","d1"),]

#elev
elev_boul= elev[c("a1","b1","c1","d1"),]
elev_clav= elev[c("RF","a1","b1","d1"),]
elev_pell= elev[c("a1","b1","c1"),]
elev_sang= elev[c("RF","a1","b1","c1","d1"),]

#envi
envi_boul= envi[c("a1","b1","c1","d1"),]
envi_clav= envi[c("RF","a1","b1","d1"),]
envi_pell= envi[c("a1","b1","c1"),]
envi_sang= envi[c("RF","a1","b1","c1","d1"),]

#phen dist
phen_boul= phen_boul[,c("PBT.m","Ctmin.m","Ctmax.m","Mean.Femur.Length","Clutch.Size","Mean.Egg.Mass","Number.Ovarioles","clutch.weight.g")]
phen_clav= phen_clav[,c("PBT.m","Ctmin.m","Ctmax.m","Mean.Femur.Length","Clutch.Size","Mean.Egg.Mass","Number.Ovarioles","clutch.weight.g")]
phen_pell= phen_pell[,c("PBT.m","Ctmin.m","Ctmax.m","Mean.Femur.Length","Clutch.Size","Mean.Egg.Mass","Number.Ovarioles","clutch.weight.g")]
phen_sang= phen_sang[,c("PBT.m","Ctmin.m","Ctmax.m","Mean.Femur.Length","Clutch.Size","Mean.Egg.Mass","Number.Ovarioles","clutch.weight.g")]

#======================================================================
#make distance matrices
#geodist<- vegdist(geo, method="euclidean") 
#elevdist<- vegdist(elev, method="euclidean") 
#envidist<- vegdist(envi, method="euclidean") 

geodist_boul<- vegdist(geo_boul, method="euclidean", na.rm=TRUE) 
geodist_clav<- vegdist(geo_clav, method="euclidean", na.rm=TRUE)
geodist_pell<- vegdist(geo_pell, method="euclidean", na.rm=TRUE)
geodist_sang<- vegdist(geo_sang, method="euclidean", na.rm=TRUE)

elevdist_boul<- vegdist(elev_boul, method="euclidean", na.rm=TRUE) 
elevdist_clav<- vegdist(elev_clav, method="euclidean", na.rm=TRUE)
elevdist_pell<- vegdist(elev_pell, method="euclidean", na.rm=TRUE)
elevdist_sang<- vegdist(elev_sang, method="euclidean", na.rm=TRUE)

envidist_boul<- vegdist(envi_boul, method="euclidean", na.rm=TRUE) 
envidist_clav<- vegdist(envi_clav, method="euclidean", na.rm=TRUE)
envidist_pell<- vegdist(envi_pell, method="euclidean", na.rm=TRUE)
envidist_sang<- vegdist(envi_sang, method="euclidean", na.rm=TRUE)

phendist_boul<- vegdist(phen_boul, method="euclidean", na.rm=TRUE) 
phendist_clav<- vegdist(phen_clav, method="euclidean", na.rm=TRUE)
phendist_pell<- vegdist(phen_pell, method="euclidean", na.rm=TRUE)
phendist_sang<- vegdist(phen_sang, method="euclidean", na.rm=TRUE)

#======================================================================
#MANTEL TESTS

#boul
mantel(phendist_boul, geodist_boul)
mantel(phendist_boul, elevdist_boul)
mantel(phendist_boul, envidist_boul)
mantel(phendist_boul, gen_boul)

#partial mantel
mantel.partial(phendist_boul, gen_boul, elevdist_boul)
mantel.partial(phendist_boul, elevdist_boul, gen_boul)
#-----

#clav
mantel(phendist_clav, geodist_clav)
mantel(phendist_clav, elevdist_clav)
mantel(phendist_clav, envidist_clav)
mantel(phendist_clav, gen_clav)

#partial mantel
mantel.partial(phendist_clav, gen_clav, elevdist_clav)
mantel.partial(phendist_clav, elevdist_clav, gen_clav)
#-----

#pell
mantel(phendist_pell, geodist_pell)
mantel(phendist_pell, elevdist_pell)
mantel(phendist_pell, envidist_pell)
mantel(phendist_pell, gen_pell)

#partial mantel
mantel.partial(phendist_pell, gen_pell, elevdist_pell)
mantel.partial(phendist_pell, elevdist_pell, gen_pell)
#-----

#sang
mantel(phendist_sang, geodist_sang)
mantel(phendist_sang, elevdist_sang)
mantel(phendist_sang, envidist_sang)
mantel(phendist_sang, gen_sang)

#partial mantel
mantel.partial(phendist_sang, gen_sang, elevdist_sang)
mantel.partial(phendist_sang, elevdist_sang, gen_sang)
#-----