# LOAD DATA BY RUNNING PhenFigs_20May2018.R

require(nlme)
require(lme4)

#--------------
#STATS

#Fig 1. Elevation patterns

#GDDS~elevation
clim2= clim1[,c("Site","Year","elevation","period","Cdd_seas") ]
clim2= na.omit(clim2)
mod1= lme(Cdd_seas~elevation*period,random=~1|Year , data=clim2)

#---
#Phenology~elevation

hop.ave= hop4[which(hop4$quantile==50),]
hop.ave$elevation= as.numeric(as.character(hop.ave$elevation))
hop.ave= as.data.frame( hop.ave[,c("elevation", "cdd_sum","site","year","ordinal","cdd_seas","species","period")] )

#ordinal
mod1= lme(ordinal~elevation*species*period,random=~1|year , data=hop.ave)
#gdd
mod2= lme(cdd_sum~elevation*species*period,random=~1|year , data=hop.ave)

#----------------
#Fig 2. Phenology

hop4q= subset(hop4, hop4$quantile==50)
hop4q.2= hop4q[,c("site","year","elevation","period","cdd_seas","cdd_sum","ordinal","species") ]
hop4q.2= na.omit(hop4q.2)

#ordinal
mod1= lm(ordinal~cdd_seas*species*elevation*period , data=hop4q.2)
#gdd
mod2= lm(cdd_sum~cdd_seas*species*elevation*period, data=hop4q.2)

#----------------
#Fig 3. DI

di.dat= dat[,c("DI","site","year","elev","period","cdd_seas","cdd_sumfall","ordinal","species","Cdd_siteave") ]
di.dat= na.omit(di.dat)

#make random var
di.dat$spsiyr= paste(di.dat$species, di.dat$site, di.dat$year, sep="")
  
#ordinal
mod1= lme(DI~ordinal*species*period*cdd_seas*elev, random=~1|spsiyr, data=di.dat)
mod1= lme(DI~ordinal*species*period*Cdd_siteave, random=~1|spsiyr, data=di.dat)
#gdd
mod2= lme(DI~cdd_sumfall*species*period*cdd_seas*elev, random=~1|spsiyr, data=di.dat)
mod2= lme(DI~cdd_sumfall*species*period*Cdd_siteave, random=~1|spsiyr, data=di.dat)
#FIX ******

#divide by elev
di.dat2= subset(di.dat, di.dat$site=="C1")

#ordinal
mod1= lme(DI~ordinal*species*period*cdd_seas, random=~1|spsiyr, data=di.dat2)
#gdd
mod2= lme(DI~cdd_sumfall*species*period*cdd_seas, random=~1|spsiyr, data=di.dat2)

#divide by species
di.dat2= subset(di.dat, di.dat$species=="Melanoplus boulderensis")

#ordinal
mod1= lme(DI~ordinal*elev*period*cdd_seas, random=~1|spsiyr, data=di.dat2)
#gdd
mod2= lme(DI~cdd_sumfall*elev*period*cdd_seas, random=~1|spsiyr, data=di.dat2)

#----------------
#Fig 4. Composition plot
#No stats?

#============================
##ordinal date
mod1= lm(ordinal~cdd_yr*species+elevation+period, data=hop4)

#by species
i=1
hop.sp= subset(hop4, hop4$species=specs[i])
mod1= lm(ordinal~cdd_yr*elevation+period, data=hop.sp)
summary(mod1)

##GDD
mod1= lm(cdd_sum~cdd_yr*species+elevation+period, data=hop4)

#by species
i=1
hop.sp= subset(hop4, hop4$species=specs[i])
mod1= lm(cdd_sum~cdd_yr*elevation+period, data=hop.sp)
summary(mod1)

#-------------
#STATS
hopq= subset(hop4, hop4$quantile==50)
mod1= lm(ordinal~cdd_yr*elevation+species+period,  data=hopq)

#split by species
i=1
hopq.sp= subset(hopq, hopq$species== specs[i])
mod1= lm(ordinal~cdd_yr*elevation+period,  data=hopq.sp)
mod1= lm(cdd_sum~cdd_yr*elevation+period,  data=hopq.sp)
summary(mod1)

#------------
