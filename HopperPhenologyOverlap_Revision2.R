library(reshape)
library(reshape2)
library(tidyr)
library(cowplot)
require(nlme)
require(lme4)
require(plyr)
library(MuMIn)
library(RColorBrewer)
library(lemon)
library(ggplot2)
library(tidyverse)
library(r2glmm)
library(sjPlot)
library(lmerTest)

## DATA WITH ALL NYMPHS
setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/GrasshopperPhenSynch/data/")
dat.all= read.csv("HopperData_Sept2019_forPhenOverlap.csv")

#Prepare data-----------------

#LOAD SPECIES DATA
setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/GrasshopperPhenSynch/data/")
spec.diapause= read.csv("SpeciesDiapause.csv")
#dat.all1= read.csv("HopperClimateData.csv")
#data currently from phenology file PhenFigs_20May2018.R
timing.mat= read.csv("timingmat_2019.csv" )
colnames(timing.mat)[1]<-"species"
rownames(timing.mat)<- timing.mat$species
#alphabetize by group
timing.mat= timing.mat[order(timing.mat$timing, timing.mat$species),]

#climate data
clim1= read.csv("Clim1Data.csv")

setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/GrasshopperPhenSynch/figures/")

#---
#compare egg diapausers vs non
dat.all$diapause="egg"
dat.all$species= as.character(dat.all$species)
dat.all$diapause[which(dat.all$species %in% c("Arphia conspersa","Eritettix simplex","Pardalaphoa apiculata","Xanthippus corallipes"))]="nymph"

#fix names
dat.all$species[which(dat.all$species=="Melanoplus bivatattus")]="Melanoplus bivittatus"

#focal species
specs= c("Aeropedellus clavatus","Camnula pellucida","Melanoplus dawsoni","Melanoplus boulderensis","Melanoplus sanguinipes","Arphia conspersa","Eritettix simplex","Pardalaphoa apiculata","Xanthippus corallipes","Melanoplus fasciatus") #"Chloealtis abdominalis",
specs= unique(dat.all$species)

#subset to focal species
##dat= subset(dat.all, dat.all$species %in% specs)
#keep all species
dat=dat.all

years= unique(na.omit(dat$year))
sites= unique(dat$site)

#elevs and sites
sites=c("CHA", "A1", "B1", "C1")
elevs= c(1752, 2195, 2591, 3048)

#Calculate development index
dat$DIp=0
inds=which(dat$total>0)  
dat$DIp[inds]= (dat$in1[inds] +dat$in2[inds]*2 +dat$in3[inds]*3 +dat$in4[inds]*4 +dat$in5[inds]*5 +dat$in6[inds]*6)/6

dat$spsiteyear= paste(dat$siteyear, dat$species, sep="")

#calculate proportion
dat.sum= ddply(dat, c("site", "year","species","spsiteyear"), summarise,
               DIp = sum(DIp, na.rm=TRUE), TotMean=mean(total, na.rm=TRUE), TotSum=sum(total, na.rm=TRUE), AdultSum=sum(in6, na.rm=TRUE), Ndates=length(total[DIp>0]) )
dat.sum= dat.sum[order(dat.sum$site,dat.sum$species, dat.sum$year),]

dat$DIptotal=NA
dat$Ndates= NA
dat$TotSum= NA
dat$AdultSum= NA
match1= match(dat$spsiteyear, dat.sum$spsiteyear)
matched= which(!is.na(match1))
dat$DIptotal[matched]<- dat.sum$DIp[match1[matched]]  
dat$Ndates[matched]<- dat.sum$Ndates[match1[matched]]
dat$TotSum[matched]<- dat.sum$TotSum[match1[matched]]
dat$AdultSum[matched]<- dat.sum$AdultSum[match1[matched]]
#normalize
dat$DIpNorm= dat$DIp/dat$DIptotal
dat$AdultNorm= dat$in6/dat$AdultSum

#drop cases with TotSum<30 or Ndates<5
dat= dat[which(dat$Ndates>5),]
dat= dat[which(dat$TotSum>30),]

#all species combinations
sp.combs= combn(specs, 2)
#switch focal species
sp.combs= cbind(sp.combs, sp.combs[c(2,1),])

#set up data storage
po= array(data = 0, dim = c(length(sp.combs[1,]), length(years), length(sites),10) ) #add dimension for multiple metrics

#combine species combinations names
sp.combs.n= paste(sp.combs[1,],sp.combs[2,],sep="_")

#add po names
dimnames(po)[[1]] <- sp.combs.n
dimnames(po)[[2]] <- years
dimnames(po)[[3]] <- sites

#dataset stats
sum(dat$total)

siteyro= paste(dat$siteyear, dat$ordinal, sep="_")
length(unique(siteyro))
#---
# plot out densities
#made season warmth a factor
#make cdd_july index
#dat$swarm= as.factor(round(dat$Cdd_siteave,2))
dat$swarm= as.factor(round(dat$Cdd_july_siteave,2))

#adults y=In6
#DIp y= DIp

#plot as proportion of total
dat$elev.lab= paste(dat$elev,"m",sep="")
dat$elev.lab= factor(dat$elev.lab, levels=c("1752m","2195m","2591m","3048m"))
dat$diapause= as.factor(dat$diapause)

#add seasonal timing (calculated below)
dat$timing= timing.mat$timing[match(dat$species, timing.mat$species)]
##drop species without timing
dat= dat[which(!is.na(dat$timing)),]

#order species by timing
dat$species= factor(dat$species, levels=timing.mat$species )

setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/GrasshopperPhenSynch/figures/")

# FIGURE 1-------------------------
#No normalization
pdf("Fig1_distyear_byspec.pdf",height = 12, width = 12)
ggplot(data=dat, aes(x=ordinal, y = DIp, group=spsiteyear, color=Cdd_siteave))+ 
  geom_smooth(method="loess", se=FALSE)+geom_point()+facet_grid(species~elev.lab, scales="free")+
  theme_bw()+theme(strip.text.y = element_text(angle = 0))+ 
  theme(legend.position="bottom", legend.key.width=unit(3,"cm"), axis.title=element_text(size=16))+
  scale_color_viridis_c()+
  #scale_color_gradientn(colours = c('blue', 'cadetblue', 'orange')) +
  xlab("ordinal date") +ylab("abundance")+ 
   labs(color = "seasonal GDDs")+ theme(strip.text = element_text(face = "italic")) 
dev.off()

# FIGURE 2---------------------------
#ABUNDANCE DISTRIBUTION ANALYSIS

#for each species, site, year
#spline interpolation?

#15th percentile

#normalize dates to wed each week
date= paste(dat$year,format(strptime(dat$ordinal, format="%j"), format="%m-%d"),sep="-")
wd= as.POSIXlt(date)$wday
dat$ord.wed= dat$ordinal +(3-wd)

#cummulative DIp
df.max<- dat %>% 
  group_by(species, site, year) %>%
  arrange(ordinal, .by_group = TRUE) %>%
  summarise(ordinal=ordinal,ord.wed=ord.wed, DIp=DIp, elev=elev, elev.lab=elev.lab, cdd_seas=cdd_seas, timing=timing,
            DIptot=DIptotal, maxDIp = max(DIp), cumDIp= cumsum(DIp),
            ord.p15=ordinal[which.max( cumDIp > DIptotal*0.15 )],
            ord.p85=ordinal[which.min( cumDIp < DIptotal*0.85 )],
            ord.p15.w=ord.wed[which.max( cumDIp > DIptotal*0.15 )],
            ord.p85.w=ord.wed[which.min( cumDIp < DIptotal*0.85 )]
            )

df.max= as.data.frame(df.max)
#aggregate
df.c= aggregate(df.max, list(df.max$species, df.max$site, df.max$year), FUN="head", 1)
names(df.c)[1:3]=c("species","site","year")
df.c=df.c[,-(4:6)]

#breadth 85th-15th percentile
df.c$breadth= df.c$ord.p85 - df.c$ord.p15
df.c$breadth.w= df.c$ord.p85.w - df.c$ord.p15.w

#check
df.2= df.max[c(df.max$species=="Melanoplus boulderensis" & df.max$site=="C1" & df.max$year==2008),]
df.2= df.2[order(df.2$ordinal),]
which.max( df.2$cumDIp > df.2$DIptot*0.15 )
which.min( df.2$cumDIp < df.2$DIptot*0.85 )
plot(df.2$ordinal, df.2$DIp)

#PLOT
#ABUND FIGURE 
figaa <- ggplot(df.c, aes(x=cdd_seas, y=log(maxDIp), colour=species, group=species)) + geom_point()+geom_smooth(method="lm", se=FALSE)+
  theme_classic() +  #geom_line() +
  facet_grid(~elev.lab, drop=TRUE, scales="free") +
  xlab("")+ylab("peak abundance")+
  scale_color_viridis_d()+ theme(legend.text = element_text(face = "italic")) #+ theme(legend.position="none")

figab <- ggplot(df.c, aes(x=cdd_seas, y=ord.p15, colour=species, group=species)) + geom_point()+geom_smooth(method="lm", se=FALSE)+
  theme_classic() +  #geom_line() +
  facet_grid(~elev.lab, drop=TRUE, scales="free") +
  xlab("")+ylab("15th percentile of abundance (day)")+
  scale_color_viridis_d()+ theme(legend.text = element_text(face = "italic")) #+ theme(legend.position="none")

figac <- ggplot(df.c, aes(x=cdd_seas, y=breadth, colour=species, group=species)) + geom_point()+geom_smooth(method="lm", se=FALSE)+
  theme_classic() +  #geom_line() +
  facet_grid(~elev.lab, drop=TRUE, scales="free") +
  xlab("")+ylab("breadth of abundance distribution (day)")+
  scale_color_viridis_d()+ theme(legend.text = element_text(face = "italic")) #+ theme(legend.position="none")

figad <- ggplot(df.c, aes(x=cdd_seas, y=log(DIptot), colour=species, group=species)) + geom_point()+geom_smooth(method="lm", se=FALSE)+
  theme_classic() +  #geom_line() +
  facet_grid(~elev.lab, drop=TRUE, scales="free") +
  xlab("seasonal GDDs (°C) ")+ylab("log total seasonal abundance")+
  scale_color_viridis_d()+ theme(legend.text = element_text(face = "italic")) #+ theme(legend.position="none")

#---
pdf("Fig2_AbundMet.pdf", height = 10, width = 8)
grid_arrange_shared_legend(figab, figac, figad, ncol = 1, nrow = 3, position='right')
dev.off()

#FIGURE 2 STATS---------------------------------------

#lmer model 
df.c$elev.ord= factor(df.c$elev.lab, ordered=FALSE, levels=c("1752m", "2195m", "2591m", "3048m") )
df.c$timing= factor(df.c$timing, ordered=TRUE, levels=c(1, 2, 3) )
df.c$sp_site= paste(df.c$species, df.c$site, sep="_")
df.c$logDIptot= log(df.c$DIptot) #to improve distribution of residuals

mo=lmer(ord.p15 ~ cdd_seas +timing +elev.ord +cdd_seas:timing +cdd_seas:elev.ord+(1|species:year), REML=TRUE, data=df.c)
mb=lmer(breadth ~ cdd_seas +timing +elev.ord +cdd_seas:timing +cdd_seas:elev.ord+(1|species:year), REML=TRUE, data=df.c)
mt=lmer(logDIptot ~ cdd_seas +timing +elev.ord +cdd_seas:timing +cdd_seas:elev.ord+(1|species:year), REML=TRUE, data=df.c)

mod1=mo
anova(mod1)
summary(mod1)$coefficients
r.squaredGLMM(mod1)
coef(mod1)

#save coefficients
mo.p=summary(mo)$coefficients
mb.p=summary(mb)$coefficients
mt.p=summary(mt)$coefficients
mo.c= rbind(mo.p,mb.p,mt.p)  
mo.c[,1:2]= signif(mo.c[,1:2],2)
mo.c[,3]= round(mo.c[,3],0)
mo.c[,4]= round(mo.c[,4],1)
mo.c[,5]= signif(mo.c[,5],2)
write.csv(mo.c, "table1coef.csv")

#save anova
mo.a= rbind(anova(mo), anova(mb),anova(mt))
mo.a[,c(1:2,5)]=round( mo.a[,c(1:2,5)],1)
mo.a[,4]=round( mo.a[,4],0)
mo.a[,6]= round(mo.a[,6],2)
write.csv(mo.a, "table1anova.csv")

#plots
mo.p=plot_model(mo, type="pred",terms=c("cdd_seas","timing","elev.ord"), show.data=TRUE, title="",axis.title=c("seasonal growing degree days (C)","Ordinal date of 15th quantile"))
mo.b=plot_model(mb, type="pred",terms=c("cdd_seas","timing","elev.ord"), show.data=TRUE, title="",axis.title=c("seasonal growing degree days (C)","Breadth of abundance distribution (days)"))
mo.t=plot_model(mt, type="pred",terms=c("cdd_seas","timing","elev.ord"), show.data=TRUE, title="",axis.title=c("seasonal growing degree days (C)","Total number individuals"))

pdf("Fig2_mods.pdf",height = 6, width = 8)
par(mfrow=c(3,1))
mo.p
mo.b
mo.t
#plot_grid(c(mo.p,mo.b,mo.t))
dev.off()

#Calculate overlap metrics---------------------

for(year in 1:length(years)){
  for(site in 1:length(sites)){
    
    #sum all individuals by species and year 
    dat.sub= dat[dat$year==years[year],]
    dat.sub= dat.sub[ which(dat.sub$site==sites[site]),]
    #restrict to adults
    #dat.sub= subset(dat.sub, dat.sub$in6>0 )
    
    #change to GDD
    #dat.sub$ordinal= dat.sub$cdd_sumfall
    
    if( nrow(dat.sub)>0 ){ #check data exists
      
      ords= unique(dat.sub$ordinal)
      
      #totals by ordinal date
      dat.agg= dat.sub[,c("species", "in6")] #extract adults
      dat.agg$Species= as.factor(dat.sub$species)
      dat.agg= aggregate(dat.agg$in6, by=list(dat.agg$species), FUN=sum) #sum number adults
      names(dat.agg)= c("Species", "YrTot")
      
      #turn into proportions
      match1= match( dat.sub$species, dat.agg$Species)
      dat.sub$YrTot= dat.agg$YrTot[match1] #Add total number adults per species per year
      dat.sub= dat.sub[which(dat.sub$YrTot>0),]
      dat.sub$AdultProp= dat.sub$in6/dat.sub$YrTot #Add proportion of annual individuals for each ordinal date, normalize species with different abundances 
      #subset to dates with adults
      dat.sub= dat.sub[which(dat.sub$AdultProp>0),]
      
      #find species combinations potentially present
      specs1= unique(dat.sub$species)
      #find species combinations present
      inds=which(sp.combs[1,]%in%specs1 & sp.combs[2,]%in%specs1)
      
      #loop through species combinations
      for(ind.k in 1:length(inds) ){
        sp1= sp.combs[1,inds[ind.k] ]
        sp2= sp.combs[2,inds[ind.k] ]
        #subset to combinations
        dat.sub2= dat.sub[which(dat.sub$species==sp1 | dat.sub$species==sp2),]
        
        #   #Pianka (1974) metric, http://www.jstor.org/stable/pdf/4217327.pdf, 
        
        #new metric?: sum across dates: (max proportion-overlap)/total, On the Analysis of Phenological Overlap
        #   #po[inds[c], year, site]= po[inds[c], year, site]+ (max(prop1,prop2)-abs(prop1-prop2))/(prop1+prop2)
       
        #DI overlap
       
        #check each species has at least 4 dates
        if( length(which(dat.sub2$species==sp1))>3 & length(which(dat.sub2$species==sp2))>3 ){
        
        #loess curve for each species
        d.sp1= dat.sub2[which(dat.sub2$species==sp1),]
        d.sp2= dat.sub2[which(dat.sub2$species==sp2),]
        
        # create loess object and prediction function
        l1 <- loess(DIp~ordinal, d.sp1)
        l2 <- loess(DIp~ordinal, d.sp2)
        
        #set negative values and NA to zero
        f1 <- function(x) {f=predict(l1,newdata=x); f[f<0]=0;f[is.na(f)]=0; return(f)}
        f2 <- function(x) {f=predict(l2,newdata=x); f[f<0]=0;f[is.na(f)]=0; return(f)}
        
        #find range
        low= min(dat.sub2$ordinal)
        up= max(dat.sub2$ordinal)
       
        # perform integration
        a1=integrate(f1,low,up)$value
        a2=integrate(f2,low,up)$value
     
        #area overlap
        #normalize to 1
        sp1.n= f1(low:up)/a1
        sp2.n= f2(low:up)/a2
        #pianka's index
        po[inds[ind.k], year, site,1]= sum(sp1.n*sp2.n)/sqrt(sum(sp1.n^2)*sum(sp2.n^2))
        
        #overlapping area normalized to total area
        #normalize to 1
        sp1.n= f1(low:up)
        sp2.n= f2(low:up)
        dinds= which(sp1.n>0 & sp2.n>0)
        if(length(dinds>0)) po[inds[ind.k], year, site,2]= sum(sapply(dinds, FUN=function(x) min(sp1.n[x],sp2.n[x])) )/(a1+a2)
        #overlapping area normalized to total area of focal species
        if(length(dinds>0)) po[inds[ind.k], year, site,3]= sum(sapply(dinds, FUN=function(x) min(sp1.n[x],sp2.n[x])) )/(a1)
        #not normalized
        if(length(dinds>0)) po[inds[ind.k], year, site,4]= sum(sapply(dinds, FUN=function(x) min(sp1.n[x],sp2.n[x])) )
        
        #area of focal species
        po[inds[ind.k], year, site,5]= a1
        
        #time difference metrics
        #days focal
        po[inds[ind.k], year, site,6]= length(which(sp1.n>0))
        #days overlap
        po[inds[ind.k], year, site,7]= length(dinds)
        #days overlap normalized
        po[inds[ind.k], year, site,8]= length(dinds)/length(which(sp1.n>0))
        #days difference in peak
        po[inds[ind.k], year, site,9]= abs(which.max(sp1.n)-which.max(sp2.n))
        #days difference in 20% of peak
        po[inds[ind.k], year, site,10]= abs(which.max(sp1.n>max(sp1.n)*0.2) - which.max(sp2.n>max(sp2.n)*0.2) )
        
        } #end check each species has at least 4 dates
      } #end loop species combinations
      
    } #end check data
  } #end loop year
} #end loop sites

#replace 0 values with NA
po[po==0]=NA

#---
#PLOTS
# PO: sp comb x years x sites

po1= melt(po, varnames=c("sp","year","site","metric"))

#drop NAs
po1= po1[which(!is.na(po1$value)),]

#add time period
po1$period="initial"
po1[which(po1$year>1960),"period"]<-"resurvey"

#ADD cdd and temp
clim1$siteyear= paste(clim1$Site, clim1$Year, sep="")
po1$siteyear= paste(po1$site, po1$year, sep="")

#select cdd metric
clim1$Cdd= clim1$Cdd_seas

po1$Tmean=NA
po1$cdd=NA
po1$cdd_july=NA
match1= match(po1$siteyear, clim1$siteyear)
matched= which(!is.na(match1))
po1$Tmean[matched]<- clim1$Mean[match1[matched]]  
po1$cdd[matched]<- clim1$Cdd[match1[matched]]
po1$cdd_july[matched]<- clim1$Cdd_july[match1[matched]]

#---
#separate species names for plotting
po1$sp= as.character(po1$sp)

sps= colsplit(po1$sp, "_", names = c('sp1', 'sp2'))
po1$sp1<- sps$sp1
po1$sp2<- sps$sp2

#make species factor for plotting
po1$sp1= factor(po1$sp1, levels=specs)
po1$sp2= factor(po1$sp2, levels=specs)

rownames(timing.mat)<- timing.mat$species
timing.mat$fits1= as.numeric(as.character(timing.mat$fits1))

#add seasonal timing (calculated below)
po1$sp1_timing= timing.mat$timing[match(po1$sp1, rownames(timing.mat))]
po1$sp2_timing= timing.mat$timing[match(po1$sp2, rownames(timing.mat))]
po1$sp1_timing_doy= round(timing.mat$fits1[match(po1$sp1, rownames(timing.mat))])
po1$sp2_timing_doy= round(timing.mat$fits1[match(po1$sp2, rownames(timing.mat))])

##drop species without timing
po1= po1[which(!is.na(po1$sp1_timing)),]
po1= po1[which(!is.na(po1$sp2_timing)),]

#PLOT
#add elevation
po1$elevation= as.factor(elevs[match(po1$site, sites)])

#FIGURE 2 OVERLAP--------------------------
#Plot overlap

#assess significance
#calculate significant regressions
#apply through combinations of species and elevations
po1$elevspec= paste(po1$elevation, po1$sp, sep="")
elevspec= matrix(unique(po1$elevspec))

#add average year CDD
fit.df$siteyear= paste(fit.df$site, fit.df$year, sep="")
po1$Cdd_siteave = fit.df[match(po1$siteyear, fit.df$siteyear),"Cdd_siteave"]

#---
#EXTRACT COEFFS
#2: overlap area normalized to 1
#7: days overlap
#9: days difference in peak
po2= po1#[which(po1$metric==metrics[metric.ind]),]

#PLOTS

#drop if less than 5 years 
po2.counts= ddply(po2, c("site", "sp","elevation","elevspec"), summarise, count=length(unique(year)) )
po2.keep= po2.counts[which(po2.counts$count>3),]
po2= po2[which(po2$elevspec %in% unique(po2.keep$elevspec)),]

#order species by timing
po2$sp1= factor(po2$sp1, levels=row.names(timing.mat) )

#labels
#add elevation label
po2$elev.lab= paste(po2$elevation, "m", sep="")
metric.lab= c("overlapping area (proportion)", "days of overlap", "days difference in peak abundance" )

metric=2 #2,7,9

#plot by elevation
if(metric==2) pdf("Fig3S_PhenOverlap_byGDD_overlaparea.pdf", height = 10, width = 10)
if(metric==7) pdf("Fig3S_PhenOverlap_byGDD_daysoverlap.pdf", height = 10, width = 10)
if(metric==9) pdf("Fig3S_PhenOverlap_byGDD_daysdiffpeak.pdf", height = 10, width = 10)
ggplot(data=po2[po2$metric==metric,], aes_string(x=x1, y = "value", group="sp2", color="sp2_timing_doy"))+geom_point(aes(shape=period, fill=period), size=2)+
  facet_grid(sp1~elev.lab, drop=TRUE, scales="free_x")+theme_bw()+geom_smooth(method="lm", se=FALSE)+
  theme(strip.text.y = element_text(angle = 0)) + 
  scale_color_gradientn(colours = c('blue', 'cadetblue', 'orange')  ) +
  xlab("seasonal growing degree days") +ylab(metric.lab[metric.ind])+ 
  theme(strip.text = element_text(face = "italic")) +
  labs(color = "seasonal timing \n (doy)")
#, aes(linetype=sig.gdd)
dev.off()

#+ylab("Phenological overlap (days)")+xlab("Growing degree days")

#FIGURE 2 OVERLAP STATS--------------------------------------
#STATS ACROSS SPECIES
#drop symmetric cases
spsort= t(apply(po2[,c("sp1","sp2")], MARGIN=1, FUN=sort))
po2$sps= paste(spsort[,1],spsort[,2], sep="_")
spbind= paste(spsort[,1],spsort[,2], po2$site, po2$year, po2$metric, sep="_")
reps= duplicated(spbind)
po3= po2[which(reps==FALSE),]

po3$sp1_sp2= paste(po3$sp1_timing,po3$sp2_timing, sep="_")
#relabel
po3[which(po3$sp1_sp2=="2_1"),"sp1_sp2"]="1_2"
po3[which(po3$sp1_sp2=="2_3"),"sp1_sp2"]="3_2"
po3[which(po3$sp1_sp2=="3_1"),"sp1_sp2"]="1_3"
po3[which(po3$sp1_sp2=="3_2"),"sp1_sp2"]="2_3"
#rename
po3[which(po3$sp1_sp2=="1_1"),"sp1_sp2"]="nymphal_nymphal"
po3[which(po3$sp1_sp2=="1_2"),"sp1_sp2"]="nymphal_early"
po3[which(po3$sp1_sp2=="1_3"),"sp1_sp2"]="nymphal_late"
po3[which(po3$sp1_sp2=="2_2"),"sp1_sp2"]="early_early"
po3[which(po3$sp1_sp2=="2_3"),"sp1_sp2"]="early_late"
po3[which(po3$sp1_sp2=="3_3"),"sp1_sp2"]="late_late"
#arrange
po3$sp1_sp2= factor(po3$sp1_sp2, levels=c("nymphal_nymphal", "nymphal_early", "nymphal_late", "early_early",
                                          "early_late","late_late") )

#make factors
po3$elev.ord= factor(po3$elevation, ordered=FALSE, levels=c(1752, 2195, 2591, 3048) )
po3$sp1_timing= factor(po3$sp1_timing, ordered=TRUE, levels=c(1,2,3) )
po3$sp2_timing= factor(po3$sp1_timing, ordered=TRUE, levels=c(1,2,3) )

options(contrasts=c("contr.sum","contr.poly"))

mover.2=lmer(value ~ cdd +sp1_sp2+ elev.ord +
            cdd:sp1_sp2 + cdd:elev.ord + elev.ord:sp1_sp2 +
            (1|sps), REML=TRUE, data=po3[po3$metric==2,])
mover.7=lmer(value ~ cdd +sp1_sp2+ elev.ord +
               cdd:sp1_sp2 + cdd:elev.ord + elev.ord:sp1_sp2 +
               (1|sps), REML=TRUE, data=po3[po3$metric==7,])
mover.9=lmer(value ~ cdd +sp1_sp2+ elev.ord +
               cdd:sp1_sp2 + cdd:elev.ord + elev.ord:sp1_sp2 +
               (1|sps), REML=TRUE, data=po3[po3$metric==9,])

mod1= mover.7
anova(mod1)
summary(mod1)$coefficients
r.squaredGLMM(mod1)
coef(mod1)
lmerTest(mod1)

#save coefficients
mo.2=summary(mover.2)$coefficients
mo.7=summary(mover.7)$coefficients
mo.9=summary(mover.9)$coefficients
mo.c= rbind(mo.2,mo.7, mo.9)  
mo.c[,1:2]= signif(mo.c[,1:2],2)
mo.c[,3]= round(mo.c[,3],0)
mo.c[,4]= round(mo.c[,4],1)
mo.c[,5]= signif(mo.c[,5],2)
write.csv(mo.c, "table2coef.csv")

#save anova
mo.a= rbind(anova(mover.2), anova(mover.7),anova(mover.9))
mo.a[,c(1:2,5)]=round( mo.a[,c(1:2,5)],1)
mo.a[,4]=round( mo.a[,4],0)
mo.a[,6]= round(mo.a[,6],2)
write.csv(mo.a, "table2anova.csv")

#plot
labs= c("overlapping area (proportion)", "days of overlap", "days difference in peak abundance" )
overlap.plot=plot_model(mover.2, type="pred",terms=c("cdd","elev.ord","sp1_sp2"), show.data=FALSE, 
                        title="", legend.title = "elevation (m)",
                        axis.title=c("seasonal growing degree days (C)",labs[1]))


pdf("SpOverlapMod2.pdf", height = 10, width = 10)
overlap.plot
dev.off()

#---
#STATS

#MAKE TIMING MAT
#C1 timing
# fits1= fits[,,3,1]
# fits1= sort(rowMeans(fits1, na.rm=TRUE))
# 
# #Calculate average timing
# #early vs late species
# fits1= fits[,,2,1]
# fits1= sort(rowMeans(fits1, na.rm=TRUE))
# timing.mat= as.data.frame( cbind(colnames(fits1), fits1) )
# timing.mat$timing= c(rep(1,3),rep(2,2),rep(3,8) )
# #Drop Trimerotropis cincta due to limited data
# #timing.mat= timing.mat[1:13,]
# #ADD M. fasciatus 
# timing.mat= rbind(timing.mat, c(NA, 3))
# rownames(timing.mat)[nrow(timing.mat)]<- "Melanoplus fasciatus"
# 
# #write out
# setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/GrasshopperPhenSynch/data/")
# #write.csv(timing.mat, "timingmat_2019.csv" )
setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/GrasshopperPhenSynch/figures/")
#---

#Average overlap across community-----------------------

#AVERAGE
#by timing groups
po.tim= ddply(po3, c("site", "year","metric","period","siteyear","sp1_sp2","sp1_timing", "sp2_timing"), summarise,
               value = mean(value, na.rm=TRUE), Tmean= Tmean[1], cdd=cdd[1],cdd_july=cdd_july[1],elevation=elevation[1], Cdd_siteave= Cdd_siteave[1])
po.tim$elev.lab= paste(po.tim$elevation,"m",sep="")

#Average across community
po.comm= ddply(po3, c("site", "year","metric","period","siteyear"), summarise,
               value = mean(value, na.rm=TRUE), Tmean= Tmean[1], cdd=cdd[1],cdd_july=cdd_july[1],elevation=elevation[1], Cdd_siteave= Cdd_siteave[1])
po.comm$elev.lab= paste(po.comm$elevation,"m",sep="")

#FIG 3 COMMUNITY OVERLAP---------------------------
pdf("Fig3_CommOverlap_bytiming2.pdf", height = 10, width = 10)
ggplot(data=po.tim[which(po.tim$metric==2),], aes(x=cdd, y = value, color=elev.lab))+geom_point(aes(), size=2)+
  theme_bw()+geom_smooth(method="lm", se=FALSE)+facet_wrap(~sp1_sp2, drop=TRUE, scales="free_x")+
  scale_color_manual(breaks = c("1752m", "2195m", "2591m","3048m"),
                     values=c("darkorange3", "darkorange", "cornflowerblue","blue3")) +
  theme(legend.position="bottom",axis.title=element_text(size=16))+
  ylab(metric.lab[1])+xlab("seasonal growing degree days (C)")+ 
  guides(color=guide_legend(title="elevation (m)" ))
dev.off()

#FIGURE 4: plot across community--------------------------------------
#by year
#ggplot(data=po.comm[which(po.comm$metric==3),], aes(x=year, y = value, color=elevation))+geom_point(aes(shape=period, fill=period), size=2)+theme_bw()+geom_smooth(method="lm", se=FALSE)

#focus on 2,7,9?
plot.m2= ggplot(data=po.comm[which(po.comm$metric==2),], aes(x=cdd, y = value, color=elev.lab))+geom_point(aes(), size=2)+
  theme_bw()+geom_smooth(method="lm", se=FALSE) +theme(legend.position="none") +ylab(metric.lab[1])+xlab("")+
 scale_color_manual(breaks = c("1752m", "2195m", "2591m","3048m"),
                     values=c("darkorange3", "darkorange", "cornflowerblue","blue3"))+ 
  guides(color=guide_legend(title="elevation (m)" ))

plot.m7= ggplot(data=po.comm[which(po.comm$metric==7),], aes(x=cdd, y = value, color=elev.lab))+geom_point(aes(), size=2)+
  theme_bw()+geom_smooth(method="lm", se=FALSE) +theme(legend.position="none") +ylab(metric.lab[2])+xlab("")+
  scale_color_manual(breaks = c("1752m", "2195m", "2591m","3048m"),
                     values=c("darkorange3", "darkorange", "cornflowerblue","blue3"))+ 
  guides(color=guide_legend(title="elevation (m)" ))

plot.m9= ggplot(data=po.comm[which(po.comm$metric==9),], aes(x=cdd, y = value, color=elev.lab))+geom_point(aes(), size=2)+
  theme_bw()+geom_smooth(method="lm", se=FALSE) +theme(legend.position="bottom") +ylab(metric.lab[3])+xlab("seasonal growing degree days (C)")+
  scale_color_manual(breaks = c("1752m", "2195m", "2591m","3048m"),
                     values=c("darkorange3", "darkorange", "cornflowerblue","blue3"))+ 
  guides(color=guide_legend(title="elevation (m)" ))

#FIGURE 4
pdf("Fig4_CommOverlap.pdf", height = 10, width = 6)
cowplot::plot_grid(plot.m2, plot.m7, plot.m9, nrow=3, labels=c("a.","b.","c."), label_x=0.1)
dev.off()

#FIGURE 4 STATS--------------------------------------
library(car)
po.comm$elev.ord= factor(po.comm$elevation, ordered=FALSE, levels=c(1752, 2195, 2591, 3048) )

mco.2= lm(value~ cdd +elevation +cdd:elevation, data=po.comm[which(po.comm$metric==2),])


mco.2= lm(value~ cdd +elev.ord +cdd:elev.ord, data=po.comm[which(po.comm$metric==2),])
mco.7= lm(value~ cdd +elev.ord +cdd:elev.ord, data=po.comm[which(po.comm$metric==7),])
mco.9= lm(value~ cdd +elev.ord +cdd:elev.ord, data=po.comm[which(po.comm$metric==9),])
anova(mco.7)
summary(mco.7)
Anova(mco.9, type="II")   # Type II tests

#save coefficients
mo.2=summary(mco.2)$coefficients
mo.7=summary(mco.7)$coefficients
mo.9=summary(mco.9)$coefficients
mo.c= rbind(mo.2,mo.7, mo.9)  
mo.c[,1:2]= signif(mo.c[,1:2],2)
mo.c[,3]= round(mo.c[,3],1)
mo.c[,4]= signif(mo.c[,4],2)
write.csv(mo.c, "table3coef.csv")

#plot
community.plot=plot_model(mco.9, type="pred",terms=c("cdd","elev.ord"), show.data=TRUE, 
                        title="", legend.title = "elevation (m)",
                        axis.title=c("seasonal growing degree days (C)","overlap"))

#single elevation
mod1= lm(value~ cdd, data=po.comm[which(po.comm$metric==9 & po.comm$elevation==3048),])

#ABUNDANCE CONSEQUENCES OF OVERLAP-------------------------------
#Check abundance consequences
#is it possible to see whether species abundances decline at higher levels of overlap?

#add parameter fits to po1
po2$spsiteyear= paste(po2$sp1, po2$year, po2$site,  sep="_")
po.sp.agg= aggregate(po2, list(po2$spsiteyear, po2$metric, po2$sp1, po2$elevation, po2$sp1_timing),FUN=mean) 
colnames(po.sp.agg)[1:5]<-c("spsiteyear","metric","sp1","elevation","sp1_timing")

#CHECK MATCH
df.c$spsiteyear= paste(df.c$species, df.c$year, df.c$site,  sep="_")
match1= match(po.sp.agg$spsiteyear, df.c$spsiteyear)
is.matched= which(!is.na(match1))
po.sp.agg$maxDIp=NA
po.sp.agg$ord.p15=NA
po.sp.agg$breadth=NA
po.sp.agg$DIptot=NA

po.sp.agg[is.matched,"maxDIp"]= df.c[match1[is.matched],"maxDIp"]
po.sp.agg[is.matched,"ord.p15"]= df.c[match1[is.matched],"ord.p15"]
po.sp.agg[is.matched,"breadth"]= df.c[match1[is.matched],"breadth"]
po.sp.agg[is.matched,"DIptot"]= df.c[match1[is.matched],"DIptot"]

po.sp.agg1=po.sp.agg[,c("spsiteyear","value","metric","sp1","sp1_timing","cdd","maxDIp","DIptot","elevation")]
po.sp.agg1=na.omit(po.sp.agg1)
po.sp.agg1$sp_site= paste(po.sp.agg1$sp1, po.sp.agg1$elevation  ,sep="_")
po.sp.agg1$sp_site= paste(po.sp.agg1$sp1, po.sp.agg1$elevation  ,sep="_")

#order species by timing
po.sp.agg1$sp1_timing= factor(po.sp.agg1$sp1_timing, ordered=TRUE, levels=c(1,2,3)  )
#log metrics
po.sp.agg1$logmaxDIp= log(po.sp.agg1$maxDIp)
po.sp.agg1$logDIptot= log(po.sp.agg1$DIptot)

#plot
plot5a<- ggplot(data=po.sp.agg1[po.sp.agg1$metric==2,], aes(x=value, y = logDIptot, color=cdd))+geom_point()+
  facet_wrap(~sp1)+theme_bw()+scale_color_viridis_c(name="seasonal growing degree days (C)")+
  xlab(metric.lab[1])+ylab("log(total seasonal abundance) (individuals)")+
  theme(legend.position="bottom")+ theme(strip.text = element_text(face = "italic")) + theme(legend.key.width=unit(2,"cm"))

#models
#metrics 2,7,9
dat.scale= po.sp.agg1[po.sp.agg1$metric==2,]

#mod.scale=lme(maxDIp~ value*sp1_timing*elevation, random=~1|sp1, data=dat.scale)
mod.scale=lmer(logDIptot~ value*sp1_timing +(1|sp1), data=dat.scale)

summary(mod.scale)
#marginal R squared values are those associated with your fixed effects, the conditional ones are those of your fixed effects plus the random effects
r.squaredGLMM(mod.scale)
anova(mod.scale)

#plot
set_theme(
  base = theme_bw(), 
  legend.inside = TRUE,         # legend inside plot
  legend.pos = "bottom")  # legend position inside plot

plot5b=plot_model(mod.scale, type="pred",terms=c("value","sp1_timing"), show.data=TRUE, 
                        title="", legend.title = "species' timing",
                        axis.title=c(metric.lab[1],"log(total seasonal abundance)") )

#FIGURE 5
pdf("Fig5_CommOverlap.pdf", height = 8, width = 10)
cowplot::plot_grid(plot5a, plot5b, nrow=1, labels=c("a.","b."), rel_widths = c(1,0.6), rel_heights = c(1,0.6))
dev.off()

#Interaction surface plot-----------------------------------------------------

library(akima)

po.t<- po.tim[which(po.tim$metric==2),]
po.t$elevation= as.numeric(as.character(po.t$elevation))

#combined
fld <- with(po.t, interp(x = cdd, y = elevation, z = value, duplicate=TRUE))

gdat <- interp2xyz(fld, data.frame=TRUE)

p3d= ggplot(gdat) + 
  aes(x = x, y = y, z = z, fill = z) + 
  geom_tile() + 
  coord_equal() +
  geom_contour(color = "white", alpha = 0.5) + 
  scale_fill_distiller(palette="Spectral", na.value="white", name="overlapping area") +
  theme_bw(base_size=18)+xlab("seasonal GDDs (°C)")+ylab("elevation (m)") +theme(legend.position="bottom")

#+ coord_fixed(ratio = 0.01) 
#,trans = "log", breaks=c(1,2, 5,10,20)
#+ylim(c(0,1500))+xlim(c(-2,23))

#---
focs= c("focal: nymphal diapausers","focal: early season","focal: late season")
comps= c("compared: nymphal diapausers","compared: early season","compared: late season")

foc.lab= c("nymphal diapause:","early season:","late season:")
comp.lab= c("nymphal diapause","early season","late season")

#restrict to focal timing
for(k in 1:6){

  if(k==1) {po.sub<- po.t[which(po.t$fst.lab==focs[1] & po.t$st.lab==comps[1]),];lab<-paste(foc.lab[1],comp.lab[1],sep=" ")} 
  if(k==2) {po.sub<- po.t[which(po.t$fst.lab==focs[1] & po.t$st.lab==comps[2]),];lab<-paste(foc.lab[1],comp.lab[2],sep=" ")} 
  if(k==3) {po.sub<- po.t[which(po.t$fst.lab==focs[1] & po.t$st.lab==comps[3]),];lab<-paste(foc.lab[1],comp.lab[3],sep=" ")} 
  if(k==4) {po.sub<- po.t[which(po.t$fst.lab==focs[2] & po.t$st.lab==comps[2]),];lab<-paste(foc.lab[2],comp.lab[2],sep=" ")} 
  if(k==5) {po.sub<- po.t[which(po.t$fst.lab==focs[2] & po.t$st.lab==comps[3]),];lab<-paste(foc.lab[2],comp.lab[3],sep=" ")} 
  if(k==6) {po.sub<- po.t[which(po.t$fst.lab==focs[3] & po.t$st.lab==comps[3]),];lab<-paste(foc.lab[3],comp.lab[3],sep=" ")} 

  #combined
  fld <- with(po.sub, interp(x = cdd, y = elevation, z = value, duplicate=TRUE))
  
  gdat <- interp2xyz(fld, data.frame=TRUE)
  
  sp= ggplot(gdat) + 
    aes(x = x, y = y, z = z, fill = z) + 
    geom_tile() + 
    coord_equal() +
    geom_contour(color = "white", alpha = 0.5) + 
    scale_fill_distiller(palette="Spectral", na.value="white", name="overlap") +
    theme_bw(base_size=12)+xlab("seasonal GDDs (°C)")+ylab("elevation (m)") +theme(legend.position="bottom", legend.key.width=unit(1,"cm"))+
    ggtitle(lab)
  
  if(k==1) po.sub1<- sp
  if(k==2) po.sub2<- sp
  if(k==3) po.sub3<- sp
  if(k==4) po.sub4<- sp
  if(k==5) po.sub5<- sp
  if(k==6) po.sub6<- sp
  
  }

library(cowplot)

pdf("FigSx_CommOverlapContour.pdf", height = 10, width = 12)
plot_grid(po.sub1,po.sub2,po.sub3,po.sub4,po.sub5,po.sub6,ncol=3)
dev.off()
