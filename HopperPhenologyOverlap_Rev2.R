ggplot(data=dat, aes(x=ordinal, y = DIp, group=spsiteyear, color=Cdd_siteave))+ 
  geom_smooth(method="loess", se=FALSE)+geom_point()+facet_grid(species~elev.lab, scales="free")+
  theme_bw()+theme(strip.text.y = element_text(angle = 0))+ 
  theme(legend.position="bottom", legend.key.width=unit(3,"cm"), axis.title=element_text(size=16))+
  scale_color_viridis_c()+
  #scale_color_gradientn(colours = c('blue', 'cadetblue', 'orange')) +
  xlab("ordinal date") +ylab("abundance")+ 
  labs(color = "seasonal GDDs")+ theme(strip.text = element_text(face = "italic"))

#--------------------------

#for each species, site, year
#spline interpolation?

#15th percentile
library(tidyverse)

#cummulative DIp
df.max<- dat %>% 
  group_by(species, site, year) %>%
  arrange(ordinal, .by_group = TRUE) %>%
  summarise(ordinal=ordinal, DIp=DIp, elev=elev, elev.lab=elev.lab, cdd_seas=cdd_seas, timing=timing,
            DIptot=DIptotal, maxDIp = max(DIp), cumDIp= cumsum(DIp),
            ord.p15=ordinal[which.max( cumDIp > DIptotal*0.15 )],
            ord.p85=ordinal[which.min( cumDIp < DIptotal*0.85 )] )

df.max= as.data.frame(df.max)
#aggregate
df.c= aggregate(df.max, list(df.max$species, df.max$site, df.max$year), FUN="head", 1)
names(df.c)[1:3]=c("species","site","year")
df.c=df.c[,-(4:6)]

#breadth 85th-15th percentile
df.c$breadth= df.c$ord.p85 - df.c$ord.p15

#PLOT
figaa <- ggplot(df.c, aes(x=cdd_seas, y=maxDIp, colour=species, group=species)) + geom_point()+geom_smooth(method="lm", se=FALSE)+
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

figad <- ggplot(df.c, aes(x=cdd_seas, y=DIptot, colour=species, group=species)) + geom_point()+geom_smooth(method="lm", se=FALSE)+
  theme_classic() +  #geom_line() +
  facet_grid(~elev.lab, drop=TRUE, scales="free") +
  xlab("")+ylab("breadth of abundance distribution (day)")+
  scale_color_viridis_d()+ theme(legend.text = element_text(face = "italic")) #+ theme(legend.position="none")


#----
pdf("Fig2_AbundMet.pdf", height = 10, width = 8)
grid_arrange_shared_legend(figaa, figab, figac, ncol = 1, nrow = 3, position='right')
dev.off()

#-----
#STATS

#lme model 
df.c$elev.ord= factor(df.c$elev, ordered=TRUE, levels=c(1752, 2195, 2591, 3048) )

mod1=lme(maxDIp~ cdd_seas +timing +elev.ord + cdd_seas:timing + cdd_seas:elev.ord  +timing:elev.ord+ cdd_seas:timing:elev.ord, random=~1|species, data=df.c )
mod1=lme(ord.p15~ cdd_seas +timing +elev.ord + cdd_seas:timing + cdd_seas:elev.ord  +timing:elev.ord+ cdd_seas:timing:elev.ord, random=~1|species, data=df.c )
mod1=lme(breadth~ cdd_seas +timing +elev.ord + cdd_seas:timing + cdd_seas:elev.ord  +timing:elev.ord+ cdd_seas:timing:elev.ord, random=~1|species, data=df.c )

summary(mod1)
anova(mod1)

#======================
#abund vs overlap


