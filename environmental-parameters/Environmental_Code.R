#Compile environmental data:
setwd("/Users/sarahhu/Desktop/Projects/metaT_ALOHA_SPOT/Metadata/")
#
base<-read.csv("All_metadata.csv")
head(base[1:2,])
base$X<-NULL
colnames(base)<-c("Site", "Depth", "Temp", "Sal", "Oxy", "Flour")
head(base)
#
f1<-read.delim("flor_env_ALOHAJuly_SPOT.txt", header=TRUE)
head(f1[1:2,])
f2<-read.delim("flour_ALOHAMarch.txt",header=TRUE)
head(f2[1:2,])
library(dplyr)
f12<-inner_join(f1,f2,by="Depth")
dim(f1);dim(f2);dim(f12)
head(f12[1:2,])
range(f12$Depth)
#
head(base[1:4,])
#
range(base$Depth)
head(f12[1:2,])
library(reshape2)
long_f<-melt(f12,id.vars = "Depth")
head(long_f[1:3,])
colnames(long_f)<-c("Depth", "Site", "Fluor")
head(base[1:4,])
all<-left_join(base,long_f, by=c("Depth","Site"))
head(all)
# write.csv(all, file="tmp.csv")
#
#
# Fluoresence values were not present for all depths. So many depths are filled in an "NA". Will fix this later.
# Surface Catalina dn Port of LA, need to be filled in manually:
# Catalina: 0.40ug/L
# Port of LA: 6.60ug/L
head(all[1:3,])
str(all)
all$Fluor[all$Site=="Catalina"]=0.40
all$Fluor[all$Site=="PortofLA"]=6.60
head(all[1:4,]) #check
str(all) #remained numeric? Yes.
#
# Get average of July and March at ALOHA
aloha_avg<-subset(all, grepl("ALOHA", all$Site))
unique(aloha_avg$Site)
head(aloha_avg)
tmp<-melt(aloha_avg, id.vars = c("Depth", "Site"))
tmp_noNA<-dplyr::filter(tmp,  !is.na(value))
str(tmp_noNA)
#
aloha_calavg<-tmp_noNA %>%
  group_by(Depth, variable) %>%
  summarise(ALOHA=mean(value))%>%
  as.data.frame
head(aloha_calavg)
dim(aloha_calavg)
#
head(all)
aloha_w<-dcast(aloha_calavg, Depth~variable, fill=NA)
head(aloha_w)
aloha_w$Site<-"ALOHA"
head(aloha_w)
head(all)
# Generate environmental document with both ALOHA averages and as individual months
all_wavg<-rbind(all,aloha_w)
unique(all_wavg$Site)
#
# Make second layer with only points where samples are from:
spot<-c(1,150,885)
spotsite<-c("SPOT", "Catalina", "PortofLA")
aloha<-c(1,120,150,1000)
pt_depth<-subset(all_wavg, (Site %in% spotsite & Depth %in% spot) | (Site %in% "ALOHA" & Depth %in% aloha))
pt_depth
#
# write.csv(all_wavg, file="All_metadata.csv")
# write.csv(pt_depth, file="Discrete_metadata.csv")
#

## Commented out sections below show where in the plot I tried to place a data point for the discrete sample depths
#fill=Site, shape=Site, linetype=Site
unique(all_wavg$Site)
head(all_wavg)
# all_melt<-melt(all_wavg, id.vars = c("Site", "Depth"))
# head(all_melt)
#
plot_profile<-function(DATA, X, Y, ylabel, val_line){
  ggplot(DATA, aes(x=X,y=Y, linetype=Site))+
    geom_smooth(size=1, color="black")+
    scale_linetype_manual(values = val_line)+
    labs(title="", x="Depth (m)",y=ylabel)+
    coord_flip()+scale_x_reverse()+scale_y_continuous(position="top")+
    # geom_vline(xintercept=0, linetype=2, color="#525252")+
    # geom_vline(xintercept=120, linetype=2, color="#525252")+
    # geom_vline(xintercept=150, linetype=2, color="#525252")+
    # geom_vline(xintercept=885, linetype=2, color="#525252")+
    # geom_vline(xintercept=1000, linetype=2, color="#525252")+
    theme_bw()+
    theme(legend.position = "none")+
    NULL
}
plot_points<-function(DATA, X, Y, ylabel, val_shape){
  ggplot(DATA, aes(x=X,y=Y, shape=Site))+
    labs(title="", x="Depth (m)",y=ylabel)+
    coord_flip()+scale_x_reverse()+scale_y_continuous(position="top")+
    # geom_vline(xintercept=0, linetype=2, color="#525252")+
    # geom_vline(xintercept=120, linetype=2, color="#525252")+
    # geom_vline(xintercept=150, linetype=2, color="#525252")+
    # geom_vline(xintercept=885, linetype=2, color="#525252")+
    # geom_vline(xintercept=1000, linetype=2, color="#525252")+
    geom_point(size=3, color="black", fill="white")+
    scale_shape_manual(values = val_shape)+
    theme_bw()+
    theme(legend.position = "none")+
    NULL
}
#
# Fluoresence:
fluor_noNA<-subset(all_wavg, !(is.na(all_wavg$Fluor)))
# Vertical profile only:
verts<-c("ALOHA", "SPOT")
stations<-subset(all_wavg, Site %in% verts)
stations_fl<-subset(fluor_noNA, Site %in% verts)
#
# df with only the points, all stations and sites
month<-c("ALOHAJuly", "ALOHAMarch")
pt_station<-subset(pt_depth, !(Site %in% month))
#
library(cowplot)
point<-get_legend(plot_points(pt_station, pt_station$Depth, pt_station$Sal, "Salinity ppt", c(21,22,23,24))+theme(legend.position = "right"))
prof<-get_legend(plot_profile(stations, stations$Depth, stations$Oxy, expression("Oxygen"~mu~"mol/L"), c(1,2))+theme(legend.position = "right"))
#
# Environmental profiles for each station
plot_grid(
  plot_profile(stations, stations$Depth, stations$Temp, expression("Temperature"^o~"C"), c(1,2))+theme(axis.line.x = element_blank()),
  plot_profile(stations, stations$Depth, stations$Sal, "Salinity ppt", c(1,2))+theme(axis.line.x = element_blank()),
  plot_profile(stations, stations$Depth, stations$Oxy, expression("Oxygen"~mu~"mol/L"), c(1,2))+theme(axis.line.x = element_blank()),
  plot_profile(stations_fl, stations_fl$Depth, stations_fl$Fluor, expression("Chlorophyll"~mu~"g/L"), c(1,2))+theme(axis.line.x = element_blank())+scale_y_continuous(position = "top", limits = c(0,1))+geom_line(),
          prof,
          ncol=5, align="hv")
# Save: W:1300, H: 430
# +theme(axis.line.x = element_blank())+scale_y_continuous(position = "top", limits = c(0,1))
#
# Discrete points:
plot_grid(plot_points(pt_station, pt_station$Depth, pt_station$Temp, expression("Temperature"^o~"C"), c(21,22,23,24)),
          plot_points(pt_station, pt_station$Depth, pt_station$Sal, "Salinity ppt", c(21,22,23,24))+theme(axis.line.x = element_blank()),
          plot_points(pt_station, pt_station$Depth, pt_station$Oxy, expression("Oxygen"~mu~"mol/L"), c(21,22,23,24))+theme(axis.line.x = element_blank()),
          plot_points(pt_station, pt_station$Depth, pt_station$Fluor, expression("Chlorophyll"~mu~"g/L"), c(21,22,23,24))+theme(axis.line.x = element_blank()),
          point,
          ncol=5, align="hv")
# Save: W:1300, H: 430
#
# Chlorophyll with modified axis
plot_grid(plot_points(pt_station, pt_station$Depth, pt_station$Temp, expression("Temperature"^o~"C"), c(21,22,23,24)),
          plot_points(pt_station, pt_station$Depth, pt_station$Sal, "Salinity ppt", c(21,22,23,24))+theme(axis.line.x = element_blank()),
          plot_points(pt_station, pt_station$Depth, pt_station$Oxy, expression("Oxygen"~mu~"mol/L"), c(21,22,23,24))+theme(axis.line.x = element_blank()),
          plot_points(pt_station, pt_station$Depth, pt_station$Fluor, expression("Chlorophyll"~mu~"g/L"), c(21,22,23,24))+theme(axis.line.x = element_blank())+scale_y_continuous(position = "top", limits = c(0,1)),
          point,
          ncol=5, align="hv")
#
# last updated 11-24-2019


## Modify supplementary table S3
raw <- read.delim("environ-raw.txt")
head(raw[1:2,])
library(tidyverse)
str(raw)

unique(raw$Type)
unique(raw$Unit)

tmp <- raw %>% 
  select(-Source, -Note) %>% 
  filter(!(Type == "Mixed Layer Depth"), !(Unit == "uG/L"), !(Type == "Nitrogen")) %>% 
  unite(SAMPLE, Location:Depth, sep = " ") %>%
  pivot_wider(id_cols = c(Type, Unit), names_from = SAMPLE, values_from = Value, values_fill = list(Value = NA)) %>% 
  data.frame
head(tmp)
# View(tmp)
write.csv(tmp, file = "raw-TableS3.csv")
