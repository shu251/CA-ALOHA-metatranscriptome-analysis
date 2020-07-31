library(reshape2); library(dplyr)
library(ggplot2); library(pheatmap)
library(vegan);library(cowplot)
library(tidyverse)

load("Normed_avg_annotated_08022018.RData", verbose=TRUE) # From Zenodo

# Format to plot at Supergroup level
# unique(df_wtax$Supergroup)
df_wtax$sample<-paste(df_wtax$Location, df_wtax$Depth, sep=" ")
#
# Remove unassigned taxonomy.
df_wtax_noNA<-subset(df_wtax, !(Taxonomy %in% "Not assigned"))

# Supergroup level
super<-aggregate(df_wtax_noNA$mean_CPM, by=list(Location=df_wtax_noNA$Location, Depth=df_wtax_noNA$Depth, sample=df_wtax_noNA$sample, Supergroup=df_wtax_noNA$Supergroup), sum)

## Re-compile taxonomy breakdown for visualization purposes
compile_tax<-function(df){
    df$Nextlevel<-df$Phylum
    df$Nextlevel[df$Phylum==""]="Other"
    df$Nextlevel[df$Class=="Foraminifera"]="Foraminifera"
    df$Nextlevel[df$Class=="Acantharia"]="Acantharia"
    df$Nextlevel[df$Class=="Polycystinea"]="Polycystinea" 
    df$Nextlevel[df$Class=="Syndinians"]="Syndiniales" 
    df$Nextlevel[df$Class=="Bacillariophyceae"]="Bacillariophyceae" 
    df$Nextlevel[df$Class=="Pelagophyceae"]="Pelagophyceae" 
  df$tax_compiled<-df$Supergroup
  df$tax_compiled[df$Supergroup == "Alveolate"]<-"Other Alveolate"
  df$tax_compiled[df$Nextlevel=="Ciliate"]="Ciliate"
  df$tax_compiled[df$Nextlevel=="Dinoflagellate"]="Dinoflagellate"
  df$tax_compiled[df$Nextlevel=="Syndiniales"]="Syndiniales"
  #
  df$tax_compiled[df$Supergroup == "Archaeplastida"]<-"Other Archaeplastida"
  df$tax_compiled[df$Nextlevel=="Chlorophyta"]="Chlorophyta"
  #
  #df$tax_compiled[df$Supergroup=="Rhizaria"]<-"Other Rhizaria"
  #df$tax_compiled[df$Nextlevel=="Cercozoa"]="Cercozoa"
  #df$tax_compiled[df$Nextlevel=="Retaria"]="Retaria"
  #df$tax_compiled[df$Nextlevel=="Acantharia"]="Acantharia"
  #df$tax_compiled[df$Nextlevel=="Foraminifera"]="Foraminifera"
  #df$tax_compiled[df$Nextlevel=="Polycystinea"]="Polycystinea"
  #
  df$tax_compiled[df$Supergroup=="Stramenopile"]<-"Other Stramenopile"
  df$tax_compiled[df$Nextlevel=="MAST"]="MAST"
  df$tax_compiled[df$Nextlevel=="Bacillariophyceae"]="Diatom"
  df$tax_compiled[df$Nextlevel=="Pelagophyceae"]="Pelagophytes"
  #
  other<-c("Amoebozoa", "Cryptista", "Discoba")
  df$tax_compiled[df$Supergroup %in% other]="Other"
  return(df)
}
df_tax<-compile_tax(df_wtax_noNA)
unique(df_tax$tax_compiled)

# Generate taxonomy reference
tax_ref <- df_tax %>% 
    select(Taxonomy, Nextlevel, tax_compiled) %>% 
    distinct() %>% 
    data.frame
# write_delim(tax_ref, path = "taxonomic-assign-reference.txt", delim = "\t")
# ^ Use for downstream curation of taxonomic assignments when necessary

# Remove zero counts and sum by curated taxa levels
df_tax_no0 <- subset(df_tax, mean_CPM > 0)
df_tax_sum<-df_tax_no0 %>%
  group_by(Location, Depth, sample, tax_compiled) %>%
  summarise(SUM = sum(mean_CPM), COUNT = n()) %>%
  as.data.frame

# Factor:
tax_order_compiled<-c("Dinoflagellate","Ciliate","Syndiniales","Other Alveolate","MAST","Diatom","Pelagophytes","Other Stramenopile","Chlorophyta","Other Archaeplastida","Haptophytes","Rhizaria","Opisthokont", "Other")
tax_order_label<-c("Dinoflagellates","Ciliates","Syndiniales","Other Alveolates","MAST","Diatoms","Pelagophytes","Other Stramenopiles","Chlorophytes","Other Archaeplastid","Haptophytes","Rhizaria","Opisthokonts", "Other")
tax_order_color<-c("#612741","#b74a70","#b7757c","#eecfbf","#92462f","#bb603c","#dfa837","#ccc050","#33431e","#93b778","#61ac86","#657abb","#1c1949","#8a8d84")
df_tax_sum$TAX_ORDER<-factor(df_tax_sum$tax_compiled, levels = rev(tax_order_compiled), labels = rev(tax_order_label))
names(tax_order_color)<-(tax_order_label)
head(df_tax_sum[1:2,])
#
sample_list<-c("Catalina surface", "PortofLA surface", "SPOT surface", "ALOHA surface","ALOHA DCM", "ALOHA 150m", "SPOT 150m","SPOT 890m","ALOHA 1000m")
df_tax_sum$SAMPLE_ORDER<-factor(df_tax_sum$sample, levels=rev(sample_list))
#
plot_tax_compiled<-ggplot(df_tax_sum,aes(y = SUM,x = SAMPLE_ORDER,fill=(TAX_ORDER)))+
  geom_bar(stat="identity", position="fill", color="white")+labs(title="", x="",y="Relative abundance CPM")+
  scale_x_discrete(limits=c(), expand = c(0, 0))+
  scale_fill_manual(values=rev(tax_order_color))+
  coord_flip()+
  scale_y_continuous(expand=c(0,0))+
  theme(legend.title=element_blank(),legend.position="bottom",legend.text.align=0, axis.text = element_text(color="black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),panel.border = element_blank(), axis.line = element_line())+
  guides(fill=guide_legend(reverse=TRUE))

options(repr.plot.width = 6, repr.plot.height = 5)
# svg("taxa-barplot.svg", w = 6, h = 5)
plot_tax_compiled #w: 645  h: 490
# dev.off()

# Re-aggregate
super_phy<-aggregate(df_tax$mean_CPM, by=list(Location=df_tax$Location, 
                                                    Depth=df_tax$Depth, 
                                                    sample=df_tax$sample, 
                                                    Supergroup=df_tax$Supergroup, 
                                                    Nextlevel=df_tax$Nextlevel), sum)
# head(super_phy)
#
# Factor/order Next level
unique(super_phy$Nextlevel)
next_order<-c("Dinoflagellate","Ciliate","Syndiniales","Apicomplexa","Ochrophyta","Bacillariophyceae","Pelagophyceae","MAST","Bicosoecida","Labyrinthulomycota","Chlorophyta","Rhodophyta","Glaucophyta","Acantharia","Retaria","Cercozoa","Polycystinea","Foraminifera","Choanoflagellatea","Cryptophyta","Discosea","Euglenida","Fungi","Haptophyta","Heterolobosea","Kinetoplastea","Metazoa","Tubulinea","Variosea","Other")
next_color<-c("#67001f","#e7298a","#c994c7","#e7e1ef","#bd0026","#fc4e2a","#feb24c","#d94801","#ffeda0","#fdae6b","#238443","#78c679","#d9f0a3","#08306b","#2171b5","#6baed6","#c6dbef","#9e9ac8","#54278f","#54278f","#54278f","#54278f","#54278f","#54278f","#54278f","#54278f","#54278f","#54278f","#54278f","#525252")
super_phy$Nextlevel_order<-factor(super_phy$Nextlevel, levels = next_order)
names(next_color)<-next_order
#
library(ggplot2)
library(RColorBrewer)
#
head(super_phy)
plot_pies<-function(data, supergroup, loc, col){
  ggplot(subset(data, Supergroup %in% supergroup & sample %in% loc),aes(y=x,x=sample,fill=Nextlevel_order))+
    geom_bar(stat="identity", position="fill", color="white")+labs(title="", x=loc,y="")+
    scale_x_discrete(limits=c(), expand = c(0, 0))+
    scale_fill_manual(values = next_color)+
    coord_flip()+
    scale_y_continuous(position = "top", expand=c(0,0))+
    theme(legend.title=element_blank(),legend.position="none",legend.text.align=0, axis.text = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),panel.border = element_blank(), axis.line = element_blank(), axis.ticks = element_blank())+
    coord_polar(theta='y')+
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+
    NULL
}

#
# Generate pie charts for individual taxonomic groups
#

# Alveolates pie plots
leg<-get_legend(plot_pies(super_phy, "Alveolate", "SPOT 150m", next_color)+theme(legend.title=element_blank(),legend.position="bottom", legend.direction = "vertical"))
# svg("alveolate.svg",h=15)
plot_grid(plot_pies(super_phy, "Alveolate", "Catalina surface", next_color),
          plot_pies(super_phy, "Alveolate", "PortofLA surface", next_color),
          plot_pies(super_phy, "Alveolate", "SPOT surface", next_color),
          plot_pies(super_phy, "Alveolate", "ALOHA DCM", next_color),
          plot_pies(super_phy, "Alveolate", "ALOHA 150m", next_color),
          plot_pies(super_phy, "Alveolate", "ALOHA surface", next_color),
          plot_pies(super_phy, "Alveolate", "SPOT 890m", next_color),
          plot_pies(super_phy, "Alveolate", "SPOT 150m", next_color),
          plot_pies(super_phy, "Alveolate", "ALOHA 1000m", next_color),
          leg,
          ncol=1, nrow=10, hjust=-2,vjust=5)

# leg<-get_legend(plot_pies(super_phy, "Stramenopile", "SPOT 150m", next_color)+theme(legend.title=element_blank(),legend.position="bottom", legend.direction = "vertical"))
# # svg("Stramenopile_pies.svg",h=15)
# plot_grid(plot_pies(super_phy, "Stramenopile", "Catalina surface", next_color),
#           plot_pies(super_phy, "Stramenopile", "PortofLA surface", next_color),
#           plot_pies(super_phy, "Stramenopile", "SPOT surface", next_color),
#           plot_pies(super_phy, "Stramenopile", "ALOHA DCM", next_color),
#           plot_pies(super_phy, "Stramenopile", "ALOHA 150m", next_color),
#           plot_pies(super_phy, "Stramenopile", "ALOHA surface", next_color),
#           plot_pies(super_phy, "Stramenopile", "SPOT 890m", next_color),
#           plot_pies(super_phy, "Stramenopile", "SPOT 150m", next_color),
#           plot_pies(super_phy, "Stramenopile", "ALOHA 1000m", next_color),
#           leg,
#           ncol=1, nrow=10, rel_widths = 0.1, rel_heights = 0.01,hjust=-2,vjust=5)
# # dev.off()
# #
# leg<-get_legend(plot_pies(super_phy, "Archaeplastida", "SPOT 150m", next_color)+theme(legend.title=element_blank(),legend.position="bottom", legend.direction = "vertical"))
# # svg("Archaeplastida_pies.svg",h=15)
# plot_grid(plot_pies(super_phy, "Archaeplastida", "Catalina surface", next_color),
#           plot_pies(super_phy, "Archaeplastida", "PortofLA surface", next_color),
#           plot_pies(super_phy, "Archaeplastida", "SPOT surface", next_color),
#           plot_pies(super_phy, "Archaeplastida", "ALOHA DCM", next_color),
#           plot_pies(super_phy, "Archaeplastida", "ALOHA 150m", next_color),
#           plot_pies(super_phy, "Archaeplastida", "ALOHA surface", next_color),
#           plot_pies(super_phy, "Archaeplastida", "SPOT 890m", next_color),
#           plot_pies(super_phy, "Archaeplastida", "SPOT 150m", next_color),
#           plot_pies(super_phy, "Archaeplastida", "ALOHA 1000m", next_color),
#           leg,
#           ncol=1, nrow=10, rel_widths = 0.1, rel_heights = 0.01,hjust=-2,vjust=5)
# # dev.off()
# #
# leg<-get_legend(plot_pies(super_phy, "Rhizaria", "SPOT 150m", next_color)+theme(legend.title=element_blank(),legend.position="bottom", legend.direction = "vertical"))
# # svg("Rhizaria_pies.svg",h=15)
# plot_grid(plot_pies(super_phy, "Rhizaria", "Catalina surface", next_color),
#           plot_pies(super_phy, "Rhizaria", "PortofLA surface", next_color),
#           plot_pies(super_phy, "Rhizaria", "SPOT surface", next_color),
#           plot_pies(super_phy, "Rhizaria", "ALOHA DCM", next_color),
#           plot_pies(super_phy, "Rhizaria", "ALOHA 150m", next_color),
#           plot_pies(super_phy, "Rhizaria", "ALOHA surface", next_color),
#           plot_pies(super_phy, "Rhizaria", "SPOT 890m", next_color),
#           plot_pies(super_phy, "Rhizaria", "SPOT 150m", next_color),
#           plot_pies(super_phy, "Rhizaria", "ALOHA 1000m", next_color),
#           leg,
#           ncol=1, nrow=10, rel_widths = 0.1, rel_heights = 0.01,hjust=-2,vjust=5)

## Supplementary barplots with selected taxa at higher resolution:
head(df_wtax_noNA[1:2,])
#
class<-aggregate(df_tax$mean_CPM, by=list(Location=df_tax$Location, Depth=df_tax$Depth, 
                                          sample=df_tax$sample, Nextlevel=df_tax$Nextlevel, 
                                          Class=df_tax$Class), sum)

class_lev<-c("Dinoflagellate","Ciliate","Bacillariophyceae","Chlorophyta")
class_sub<-subset(class, Nextlevel %in% class_lev)

# Ciliate plot
cil<-subset(class_sub, Nextlevel %in% "Ciliate")
sample_list<-c("Catalina surface", "PortofLA surface", "SPOT surface", "ALOHA DCM", "ALOHA 150m", "ALOHA surface", "SPOT 150m", "ALOHA 1000m", "SPOT 890m")
cil$sample_order<-factor(cil$sample, levels=rev(sample_list))
#
plot_ciliate<-ggplot(cil,aes(y=x,x=sample_order,fill=Class))+
  geom_bar(stat="identity", position="fill", color="white")+labs(title="", x="",y="Relative abundance CPM")+
  scale_x_discrete(limits=c(), expand = c(0, 0))+
  # scale_fill_brewer(palette = "PuRd")+
  coord_flip()+
  scale_y_continuous(expand=c(0,0))+
  theme(legend.title=element_blank(),legend.position="bottom",legend.text.align=0, axis.text = element_text(color="black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),panel.border = element_blank(), axis.line = element_line())+
  guides(fill=guide_legend(reverse=TRUE))
# plot_ciliate+scale_fill_brewer(palette = "Set3")

# non-ciliate groups go to the Order level
order<-aggregate(df_tax$mean_CPM, by=list(Location=df_tax$Location, Depth=df_tax$Depth, 
                                          sample=df_tax$sample, Nextlevel=df_tax$Nextlevel,
                                          Order=df_tax$Order), sum)
order_sub<-subset(order, Nextlevel %in% class_lev)
unique(order_sub$Order)
noncil<-c("Dinoflagellate","Bacillariophyceae","Chlorophyta")
# head(order_sub)
# Plot others by order:
order_sub$sample_order<-factor(order_sub$sample, levels=rev(sample_list))
#
plot_order<-ggplot(order_sub,aes(y=x,x=sample_order,fill=Order))+
  geom_bar(stat="identity", position="fill", color="white")+labs(title="", x="",y="Relative abundance CPM")+
  scale_x_discrete(limits=c(), expand = c(0, 0))+
  guides(fill=guide_legend(reverse=TRUE))+
  coord_flip()+
  scale_y_continuous(expand=c(0,0))+
  theme(legend.title=element_blank(),legend.position="bottom",legend.text.align=0, 
        axis.text = element_text(color="black"),panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),panel.background = element_blank(),
        panel.border = element_blank(), axis.line = element_line())

# plot_order %+% subset(order_sub, Nextlevel %in% "Dinoflagellate")+scale_fill_brewer(palette = "Set3")
# #
# plot_order %+% subset(order_sub, Nextlevel %in% "Bacillariophyceae")
# #
# plot_order %+% subset(order_sub, Nextlevel %in% "Chlorophyta")+scale_fill_brewer(palette = "Set3")

options(repr.plot.width = 10, repr.plot.height = 7)

# svg("Suppl_tax_res.svg", h=11,w=17)
plot_grid(plot_order %+% subset(order_sub, Nextlevel %in% "Dinoflagellate")+scale_fill_brewer(palette = "Set3"),
          plot_ciliate+scale_fill_brewer(palette = "Set3"),
          plot_order %+% subset(order_sub, Nextlevel %in% "Chlorophyta")+scale_fill_brewer(palette = "Set3"),
          plot_order %+% subset(order_sub, Nextlevel %in% "Bacillariophyceae"),
          labels=c('a. Dinoflagellates', 'b. Ciliates', 'c. Chlorophytes', 'd. Diatoms'),
         align = "hv")
# dev.off()

load("Normed_avg_annotated_08022018.RData", verbose=T) #available from Zenodo
load("ReNorm_bytax_08022018.RData", verbose=T) #available from Zenodo

# import Kegg list information:

# Derived from available Kegg information
keggs <- read.delim("compile-raw-data/K0_geneIDs_fulllist_08192019.txt",sep="\t",fill=T)

# Curated kegg identities to hone in on specific protistan nutritional modes
custom<-read.delim("Custom_KO_list_30-04-2020.txt",header=TRUE, sep="\t")

tax_ref<-read.csv("taxonomic-assign-reference.txt",row.names=NULL)


# Join whole community data:
wholecomm_wK0<-left_join(avg_CPM, keggs,by="KO")
whole_wcat <- left_join(wholecomm_wK0, custom, by="KO")
# head(whole_wcat)

# Option to save and use R object
# save(whole_wcat, file="NormedbyWhole_annotated.RData")

# Join with data that is normalized by individual taxonomic groups
combo_tax_wK0<-left_join(comboTax, keggs,by="KO")

# tax_wcat - KO IDs with taxonomic ID
tax_wcat<-left_join(combo_tax_wK0, custom, by="KO")

# Option to save
# save(tax_wcat, file="Normedbytax_annotated.RData")

library(UpSetR)

melt_wKO<-melt(df_wKO_wdups)
melt_wKO$variable<-NULL
# head(melt_wKO)
melt_wKO_noNA<-subset(melt_wKO, !(Taxonomy %in% "Not assigned" | KO %in% "Not assigned" | value == 0))
# head(melt_wKO_noNA) 
# removed zeros and unassigned, KO IDs are duplicated for all applicable modules.
tmp<-melt_wKO_noNA[c(1:5,10)]
melt_wKO_noNA_noDup<-tmp[!duplicated(tmp),]
# head(melt_wKO_noNA_noDup)
toUpset<-melt_wKO_noNA_noDup[c(1:2,3:4,6)]
head(toUpset) # Use below to add annotation

# Format for UpSetR plots:
toUpset$Uniq<-paste(toUpset$Taxonomy, toUpset$KO, sep="_")
toUpset$sample<-paste(toUpset$Location, toUpset$Depth, sep="_")
# Change to binary
summed<-aggregate(toUpset$value, by=list(sample=toUpset$sample, Uniq=toUpset$Uniq),sum)
summed$bin<-ifelse(summed$x > 0, 1, 0)
binary_wide<-dcast(summed[c(1,2,4)], Uniq~sample, fill=0)
row.names(binary_wide)<-binary_wide$Uniq; binary_wide$Uniq<-NULL
head(binary_wide)

options(repr.plot.width = 16, repr.plot.height = 7)
upset(binary_wide, keep.order = T, group.by=c("both"), nsets=13, number.angles = 45, text.scale=2, intersections=list(list("ALOHA_surface","Catalina_surface","PortofLA_surface","SPOT_150m","SPOT_890m","SPOT_surface","ALOHA_150m","ALOHA_DCM","ALOHA_1000m"),list("ALOHA_surface","ALOHA_150m","ALOHA_DCM","ALOHA_1000m"),list("ALOHA_surface","ALOHA_150m","ALOHA_DCM"),list("ALOHA_150m","ALOHA_DCM"),list("ALOHA_surface","ALOHA_DCM"),list("Catalina_surface","PortofLA_surface","SPOT_150m","SPOT_890m","SPOT_surface"),list("Catalina_surface","PortofLA_surface","SPOT_surface"),list("Catalina_surface","PortofLA_surface"),list("SPOT_150m","SPOT_890m","SPOT_surface"),list("SPOT_150m","SPOT_890m"),list("ALOHA_surface","Catalina_surface","PortofLA_surface","SPOT_surface"),list("ALOHA_surface","Catalina_surface","PortofLA_surface","SPOT_surface","ALOHA_150m","ALOHA_DCM"),list("SPOT_150m","SPOT_890m","ALOHA_1000m"),list("SPOT_150m","SPOT_890m","ALOHA_150m","ALOHA_DCM","ALOHA_1000m"),list("ALOHA_surface"),list("Catalina_surface"),list("PortofLA_surface"),list("SPOT_150m"),list("SPOT_890m"),list("SPOT_surface"),list("ALOHA_150m"),list("ALOHA_DCM"),list("ALOHA_1000m")))

# Dissect binary wide further - color by specific categories:
#
tax <- read.delim("taxonomic-assign-reference.txt")
# head(tax)
# head(binary_wide[1:2,])
binary_wide$annot<-row.names(binary_wide) # pull out row names as a category to annotate
tmp<-colsplit(binary_wide$annot, "_", c("Taxonomy", "KO"));head(tmp)
tmp_wtax<-left_join(tmp, tax, by="Taxonomy")

df<-binary_wide[1:9]
# head(df)
df$Intersect <- apply(df > 0, 1, function(x){toString(names(df)[x])})
# head(df[1:3,])

binary_tax <- df %>%
    rownames_to_column(var = "uniq") %>% 
    separate(uniq, c("Taxonomy", "KEGGID"), sep = "_", remove = FALSE) %>% 
    inner_join(tax, by = "Taxonomy") %>% 
    column_to_rownames(var = "uniq") %>% 
    data.frame
head(binary_tax[1:2,])
unique(binary_tax$tax_compiled)

# Factor plot
intersect_order<-c("ALOHA_1000m, ALOHA_150m, ALOHA_DCM, ALOHA_surface, Catalina_surface, PortofLA_surface, SPOT_150m, SPOT_890m, SPOT_surface","ALOHA_1000m, ALOHA_150m, ALOHA_DCM, ALOHA_surface","ALOHA_150m, ALOHA_DCM, ALOHA_surface","ALOHA_150m, ALOHA_DCM","ALOHA_DCM, ALOHA_surface","Catalina_surface, PortofLA_surface, SPOT_150m, SPOT_890m, SPOT_surface","Catalina_surface, PortofLA_surface, SPOT_surface","Catalina_surface, PortofLA_surface","SPOT_150m, SPOT_890m, SPOT_surface","SPOT_150m, SPOT_890m","ALOHA_surface, Catalina_surface, PortofLA_surface, SPOT_surface","ALOHA_150m, ALOHA_DCM, ALOHA_surface, Catalina_surface, PortofLA_surface, SPOT_surface","ALOHA_1000m, SPOT_150m, SPOT_890m","ALOHA_1000m, ALOHA_150m, ALOHA_DCM, SPOT_150m, SPOT_890m","ALOHA_surface","Catalina_surface","PortofLA_surface","SPOT_150m","SPOT_890m","SPOT_surface","ALOHA_150m","ALOHA_DCM","ALOHA_1000m")
binary_tax$order_x<-factor(binary_tax$Intersect, levels = intersect_order) 
tax_order<-c("Dinoflagellate","Ciliate","Syndiniales","Other Alveolate","MAST","Diatom","Pelagophytes","Other Stramenopile","Chlorophyta","Other Archaeplastida","Haptophytes","Rhizaria","Opisthokont", "Other")
tax_order_label<-c("Dinoflagellates","Ciliates","Syndiniales","Other Alveolates","MAST","Diatoms","Pelagophytes","Other Stramenopile","Chlorophytes","Other Archaeplastid","Haptophytes","Rhizaria","Opisthokonts", "Other")
tax_order_color<-c("#612741","#b74a70","#b7757c","#eecfbf","#92462f","#bb603c","#dfa837","#ccc050","#33431e","#93b778","#61ac86","#657abb","#1c1949","#8a8d84")
binary_tax$TAX_ORDER<-factor(binary_tax$tax_compiled, levels = (tax_order), labels = (tax_order_label))
names(tax_order_color)<-(tax_order_label)

plot_tax_bin<-ggplot(binary_tax, aes(x=order_x, fill=TAX_ORDER))+
  geom_bar(stat = "count", position="stack")+
  theme(axis.text.x = element_text(angle=45, size=0.8))+
  scale_x_discrete(limits=c(), expand = c(0, 0))+
  scale_fill_manual(values=(tax_order_color))+
  scale_y_continuous(position = "top", expand=c(0,0))
#
options(repr.plot.width = 16, repr.plot.height = 7)
plot_tax_bin %+% subset(binary_tax, Intersect %in% intersect_order) #w:1080 h:650

custom_paths <- read.delim("Custom_KO_list_30-04-2020.txt", header = TRUE)
load("Normedbytax_annotated.RData", verbose = T)
tax_wcat$Category<- NULL
tax_wcat_subset <- inner_join(tax_wcat, custom_paths, by="KO")
# head(tax_wcat_subset)

# Write function to subset specific taxa
# library(tidyverse)
generate_pcoa <- function(df, tax) {
  tax_wide <- df %>%
    filter(taxa == tax) %>%
    group_by(sample, KO) %>%
    summarise(MEAN = mean(mean_count)) %>%
    pivot_wider(names_from = KO, values_from = MEAN) %>%
    data.frame
  row.names(tax_wide) <- tax_wide$sample
  tax_wide$sample <- NULL
  jac_bytax <- vegdist(tax_wide, method = "jaccard")
  pcoa_bytax <- princomp(jac_bytax)
  plot(pcoa_bytax)
  return(pcoa_bytax)
}
##
##
# Function to plot output PCOA
plot_pcoa <- function(pcoa_in, tax) {
  tmp <- data.frame(pcoa_in$scores)
  tmp$SAMPLE <- row.names(tmp)
  df <- separate(tmp, col = SAMPLE, c("Location", "Depth"), sep = "_", remove = FALSE)
  df$TAXA <- tax
  eigenvalues<-round(pcoa_in$sdev, 4)*100
  # Factoring
  depth_order <- c("surface", "DCM", "150m", "890m", "1000m")
  df$DEPTH_ORDER <-factor(df$Depth, levels = depth_order)
  DEPTH_COLOR <- c("#f768a1","#fe9929","#41ab5d", "#084594", "#084594")
  names(DEPTH_COLOR)<-depth_order
  loc_order <-c("ALOHA", "Catalina", "PortofLA", "SPOT")
  df$LOC_ORDER <- factor(df$Location, levels = loc_order)
  LOC_SHAPE <- c(5, 18, 15, 19)
  names(LOC_SHAPE)<-loc_order
  # Factoring
#   sample_list<-c("Catalina_surface", "PortofLA_surface", "SPOT_surface", "ALOHA_surface","ALOHA_DCM", "ALOHA_150m", "SPOT_150m","SPOT_890m","ALOHA_1000m")
#   df$SAMPLE_ORDER<-factor(df$SAMPLE, levels=sample_list)
#   sample_alpha <- c(1,1,1,0.5,0.5,0.5,1,1,0.5)
#   depth_order <- c("surface", "DCM", "150m", "890m", "1000m")
#   df$DEPTH_ORDER <-factor(df$Depth, levels = depth_order)
#   DEPTH_SHAPE <- c(24,25,21,22,23)
#   loc_order <-c("ALOHA", "Catalina", "PortofLA", "SPOT")
#   df$LOC_ORDER <- factor(df$Location, levels = loc_order)
#   LOC_COL <-c("#c51b7d","#fee08b","#1a9850","#4575b4")
#   names(LOC_COL)<-loc_order
#   tax_order<-c("Dinoflagellate","Ciliate","MAST","Bacillariophyceae","Pelagophyceae","Chlorophyta","Haptophyta","Rhizaria")
#   df$TAX_ORDER <- factor(df$TAXA, levels= tax_order)
#   TAX_COL<-c("#612741","#b74a70","#92462f","#bb603c","#dfa837","#33431e","#61ac86","#657abb")
#   names(TAX_COL)<-tax_order
  ggplot(df, aes(x=Comp.1, y=Comp.2, fill=DEPTH_ORDER, shape = LOC_ORDER)) +
    geom_hline(yintercept=0, color="black")+
    geom_vline(xintercept=0, color="black")+
    geom_point(stat="identity", size=4, aes(fill=DEPTH_ORDER, color = DEPTH_ORDER, shape=LOC_ORDER))+
    scale_fill_manual(values = DEPTH_COLOR)+
    scale_color_manual(values = DEPTH_COLOR)+
    scale_shape_manual(values = LOC_SHAPE)+
    labs(x=eigenvalues[1], y=eigenvalues[2])+
    theme_minimal()+
    theme(axis.text = element_text(color="black", face="bold"),
          panel.grid.major = element_line(color="grey"),
         legend.position = "none")+
    guides(fill = guide_legend(override.aes = list(shape = 21)),
           shape = guide_legend(override.aes = list(fill = "black")))+
    scale_x_continuous(limits = c(-1.5, 1.5))+
    scale_y_continuous(limits = c(-1, 1))
}
unique(tax_wcat_subset$taxa)

options(repr.plot.width = 7, repr.plot.height = 6)
pc_dinos <- generate_pcoa(tax_wcat_subset, "Dinoflagellate")
pc_ciliate <- generate_pcoa(tax_wcat_subset, "Ciliate")
pc_mast <- generate_pcoa(tax_wcat_subset, "MAST")
pc_rhiz <- generate_pcoa(tax_wcat_subset, "Rhizaria")
pc_hapto <- generate_pcoa(tax_wcat_subset, "Haptophyta")
pc_diat <- generate_pcoa(tax_wcat_subset, "Bacillariophyceae")
pc_chloro <- generate_pcoa(tax_wcat_subset, "Chlorophyta")
pc_pela <- generate_pcoa(tax_wcat_subset, "Pelagophyceae")

options(repr.plot.width = 14, repr.plot.height = 6)

# svg("pca-bytax-supp.svg", w = 14, h = 6)
plot_grid(plot_pcoa(pc_dinos, "Dinoflagellate"),
          plot_pcoa(pc_ciliate, "Ciliate"),
          plot_pcoa(pc_mast, "MAST"),
          plot_pcoa(pc_rhiz, "Rhizaria"),
          plot_pcoa(pc_hapto, "Haptophyta"),
          plot_pcoa(pc_diat, "Bacillariophyceae"),
          plot_pcoa(pc_chloro, "Chlorophyta"),
          plot_pcoa(pc_pela, "Pelagophyceae"),
          ncol = 4, align = "hv", 
          labels = c("a. Dinoflagellates", "b. Ciliates", "c. MAST", 
                     "d. Rhizaria", "e. Haptophytes", "f. Diatoms", "g. Chlorophytes", "h. Pelagophytes"))
# dev.off()

load("NormedbyWhole_annotated.RData", verbose = T)
tax <- read.delim("taxonomic-assign-reference.txt")

# head(whole_wcat)
# head(tax)

unigene_prof <- whole_wcat %>% 
    unite(SAMPLE, Location, Depth, sep = "_", remove = FALSE) %>% 
    unite(KO_TAX, KO, Taxonomy, sep = "_", remove = FALSE) %>% 
    # If grouping at curated taxonomic level, import tax ref and use to group.
#     left_join(tax, by = "Taxonomy") %>% 
#     filter(!is.na(tax_compiled)) %>% 
    select(SAMPLE, KO_TAX, mean_CPM) %>% 
    distinct() %>% 
    filter(mean_CPM > 0) %>% 
    pivot_wider(id_cols = SAMPLE, names_from = KO_TAX, values_from = mean_CPM, 
                values_fill = list(mean_CPM = 0),
                values_fn = sum) %>% 
    column_to_rownames(var = "SAMPLE") %>% 
    data.frame
head(unigene_prof[1:5])

# Jaccard dist metric:
library(vegan)
jac<-vegdist(unigene_prof, method="jaccard")
pcoa<-princomp(jac)
plot(pcoa)

tmp<-data.frame(pcoa$scores)
tmp$SAMPLE<-row.names(tmp)
out<-colsplit(tmp$SAMPLE, "_", c("Loc", "Depth"))
pcoa_all_df<-data.frame(tmp, out)
eigenvalues<-round(pcoa$sdev, 4)*100
unique(pcoa_all_df$Depth)
# Factoring
depth_order <- c("surface", "DCM", "150m", "890m", "1000m")
pcoa_all_df$DEPTH_ORDER <-factor(pcoa_all_df$Depth, levels = depth_order)
DEPTH_COLOR <- c("#f768a1","#fe9929","#41ab5d", "#084594", "#084594")
names(DEPTH_COLOR)<-depth_order
# DEPTH_SHAPE <- c(24,25,21,22,23)
#
loc_order <-c("ALOHA", "Catalina", "PortofLA", "SPOT")
pcoa_all_df$LOC_ORDER <- factor(pcoa_all_df$Loc, levels = loc_order)
# LOC_COL <-c("#c51b7d","#fee08b","#1a9850","#4575b4")
LOC_SHAPE <- c(5, 18, 15, 19)
names(LOC_SHAPE)<-loc_order
# loc_alpha <- c(0.5, 1, 1, 1)

byloc_PCA <- ggplot(pcoa_all_df, aes(x=Comp.1, y=Comp.2, fill=DEPTH_ORDER, shape=LOC_ORDER)) +
  geom_hline(yintercept=0, color="#525252")+
  geom_vline(xintercept=0, color="#525252")+
  geom_point(stat="identity", size=6, aes(fill=DEPTH_ORDER, color = DEPTH_ORDER, shape=LOC_ORDER))+
  scale_fill_manual(values = DEPTH_COLOR)+
  scale_color_manual(values = DEPTH_COLOR)+
  scale_shape_manual(values = LOC_SHAPE)+
  labs(x=paste(eigenvalues[1], "%"), y=paste(eigenvalues[2], "%"))+
  theme_minimal()+
  theme(axis.text = element_text(color="black", face="bold", size = 12),
        axis.title=element_text(size=14,face="bold"),
        panel.grid.major = element_line(color="grey"),
        legend.title = element_blank())+
#   scale_x_continuous(limits = c(-0.5, 0.75))+
#   scale_y_continuous(limits = c(-0.5, 0.25))+
  guides(fill = guide_legend(override.aes=list(shape=22)))

options(repr.plot.width = 6, repr.plot.height = 5)
# svg("pca-metaT-all.svg", w = 6, h = 5)
byloc_PCA
# dev.off()

library(plotly)

# save(pcoa_all_df, eigenvalues, pcoa, file = "pcoa-input-2020.RData")

pca_3d <- plot_ly(pcoa_all_df,x=~Comp.1,y=~Comp.2,z=~Comp.3,symbol=~Loc, color=~Depth,
      colors=c("#084594","#41ab5d","#084594","#fe9929","#f768a1"), symbols = c(5,18,15,19))%>%
  layout(title='PCA unigene profiles - ALOHA & CA',
         scene=list(xaxis=list(title=paste0('PC1 ',eigenvalues[1],'%'),
                               scale=pcoa$sdev[1]),
                    yaxis=list(title=paste0('PC2 ',eigenvalues[2],'%'),
                               scale=pcoa$sdev[2]),
                    zaxis=list(title=paste0('PC3 ',eigenvalues[3],'%'),
                               scale=pcoa$sdev[3])))

options(repr.plot.width = 6, repr.plot.height = 20)
pca_3d

# load("NormedbyWhole_annotated.RData", verbose=T)
head(whole_wcat[1:2,])

unique(whole_wcat$Target) #Upper level categories
unique(whole_wcat$Target_2) 
unique(whole_wcat$Category) # Use for categories to plot

by_category_WHOLE <- whole_wcat %>% 
    unite(SAMPLE, Location, Depth, sep = "_", remove = FALSE) %>% 
    select(SAMPLE, Taxonomy, Target, Category, KO, Gene_identification, mean_CPM) %>% 
    filter(!is.na(Category)) %>% 
    filter(mean_CPM > 0) %>% 
    distinct() %>% 
    group_by(SAMPLE, Taxonomy, Target, Category, KO, Gene_identification) %>%
    summarise(unigene_sum=sum(mean_CPM)) %>%
    na.omit() %>%
    data.frame
head(by_category_WHOLE[1:3,])
unique(by_category_WHOLE$Target)

unique(by_category_WHOLE$Target)
unique(by_category_WHOLE$Category)

# Total KOs
length(unique(by_category_WHOLE$KO))
#
# Subset to relevant categories
unique(by_category_WHOLE$Target)
rm <- c("butyryl", "other", "else", "additional breakdown")
whole_targetsonly <- by_category_WHOLE %>% 
    filter(!(Target == "Phagotrophy-other")) %>% 
    filter(!(Category %in% rm)) %>% 
    data.frame
    
# Sum to show whole community profiles for targetted pathways
whole_targetSUM <- whole_targetsonly %>%
  group_by(SAMPLE, Target, Category) %>%
  summarise(SUM = sum(unigene_sum)) %>%
  data.frame
head(whole_targetSUM)

unique(whole_targetSUM$Target)
unique(whole_targetSUM$Category)
length(unique(whole_targetSUM$Category))

options(repr.plot.width = 6, repr.plot.height = 6)
hist(log(whole_targetSUM$SUM)) #Normal distribution

# Format for pheatmap plots
# Make heat map:
whole_targetSUM$Target <- NULL
#
df_w<-dcast(whole_targetSUM,Category~SAMPLE,fill=0, fun.aggregate = sum)
head(df_w[1:2,])
row.names(df_w)<-df_w$Category

#Categories of interest:
all <- c("Catalina_surface", "PortofLA_surface", "SPOT_surface", "ALOHA_surface","ALOHA_DCM", "ALOHA_150m", "SPOT_150m", "ALOHA_1000m", "SPOT_890m")
#
ALOHA <- c("ALOHA_surface","ALOHA_DCM", "ALOHA_150m", "ALOHA_1000m")
#
CA <- c("Catalina_surface", "PortofLA_surface", "SPOT_surface","SPOT_150m", "SPOT_890m")
#
surface_only <- c("Catalina_surface", "PortofLA_surface", "SPOT_surface", "ALOHA_surface")
surface_wDCM <- c("Catalina_surface", "PortofLA_surface", "SPOT_surface", "ALOHA_surface","ALOHA_DCM")
all_euphotic <- c("Catalina_surface", "PortofLA_surface", "SPOT_surface", "ALOHA_surface","ALOHA_DCM", "ALOHA_150m")
subsurface <- c("ALOHA_DCM", "ALOHA_150m", "SPOT_150m", "ALOHA_1000m", "SPOT_890m")
subeuphotic <- c("SPOT_150m", "ALOHA_1000m", "SPOT_890m")

# Function to plot pheatmaps
pheat_bySAMPLE<-function(df, subset_samples, title){
  df_subset <- subset(df, SAMPLE %in% subset_samples)
  df_w<-dcast(df_subset,Category~SAMPLE,fill=0, fun.aggregate = sum)
  row.names(df_w)<-df_w$Category
  df_w<-df_w[subset_samples]
  df_m<-as.matrix(df_w)
  out<-pheatmap(df_m, scale = "row", cluster_cols = FALSE,cluster_rows = TRUE, cellwidth=9, cellheight = 9, main= title)
  return(out)
}

options(repr.plot.width = 6, repr.plot.height = 10)
a<-pheat_bySAMPLE(whole_targetSUM, all, "All")
# b<-pheat_bySAMPLE(whole_targetSUM, ALOHA, "ALOHA")
# c<-pheat_bySAMPLE(whole_targetSUM, CA, "CA")
# d<-pheat_bySAMPLE(whole_targetSUM, surface_only, "Surface only")
# e<-pheat_bySAMPLE(whole_targetSUM, surface_wDCM, "Surface w/ DCM")
# f<-pheat_bySAMPLE(whole_targetSUM, all_euphotic, "All euphotic")
# g<-pheat_bySAMPLE(whole_targetSUM, subsurface, "All subsurface")
# h<-pheat_bySAMPLE(whole_targetSUM, subeuphotic, "Subeuphotic")

# ## Save all pheatmap outputs
# save_pheatmap_png <- function(x, filename, width=1200, height=1000, res = 150) {
#   png(filename, width = width, height = height, res = res)
#   grid::grid.newpage()
#   grid::grid.draw(x$gtable)
#   dev.off()
# }
# save_pheatmap_png(a, "wholecommunity.png")
# save_pheatmap_png(b, "NPSG.png")
# save_pheatmap_png(c, "CA.png")
# save_pheatmap_png(d, "surface.png")
# save_pheatmap_png(e, "surface_wDCM.png")
# save_pheatmap_png(f, "euphotic.png")
# save_pheatmap_png(g, "subsurface.png")
# save_pheatmap_png(h, "subeuphotic.png")

load("Normedbytax_annotated.RData",verbose=T)
# data includes normalized BY individual groups
head(tax_wcat[1:2,])

unique(tax_wcat$taxa); unique(tax_wcat$Category)
tax_wcat$Category <- as.character(tax_wcat$Category)

# Subset to relevant categories
rm <- c("butyryl", "other", "else", "additional breakdown", "por")
tax_plot_paths <- tax_wcat %>% 
    filter(!(Target == "Phagotrophy-other")) %>% 
    filter(!(Category %in% rm)) %>% 
    # Generate combined gluconeogenesus-glycolysis
    mutate(Category = case_when(
        grepl("Gluconeogenesis", Category) ~ "Gluconeogenesis-Glycolysis",
        grepl("Glycolysis", Category) ~ "Gluconeogenesis-Glycolysis",
        TRUE ~ Category)) %>% 
    select(Taxonomy, taxa, KO, mean_count, sample, Target, Category) %>% 
    filter(mean_count > 0) %>% 
    distinct() %>% 
    group_by(taxa, sample, Target, Category) %>%
    summarise(SUM_CPM = sum(mean_count)) %>%
    data.frame
head(tax_plot_paths[1:2,])

# Factoring
sample_list<-c("Catalina_surface", "PortofLA_surface", "SPOT_surface", "ALOHA_surface","ALOHA_DCM", "ALOHA_150m", "SPOT_150m","SPOT_890m","ALOHA_1000m")
sample_label<-c("Catalina surface", "Port of LA surface", "SPOT surface", "ALOHA surface","ALOHA DCM", "ALOHA 150m", "SPOT 150m","SPOT 890m","ALOHA 1000m")
tax_plot_paths$SAMPLE_ORDER<-factor(tax_plot_paths$sample, levels=rev(sample_list), labels = rev(sample_label))

category_list <- c("Photosynthesis","Calvin cycle","Gluconeogenesis-Glycolysis","Glyoxylate cycle","TCA cycle","GS/GOGAT","AMT","NRT","Inorganic N uptake and assimilation","Organic N uptake and assimilation","Urea cycle","Nitrate reduction (assimilatory)","P metabolism","PDH","Actin polymerization","Lysosome binding and processing","Chitinase","Phagosome maturation","Endocytosis","Fatty acid biosynthesis","Fatty acid breakdown","Motility and prey recognition","SNARE complex","V-type ATPase")
category_color <- c("#ffffb2","#fed976","#feb24c","#fd8d3c","#f03b20","#bd0026","#f7fcb9","#addd8e","#41ab5d","#238443","#006837","#004529","#02818a","#a6bddb","#41b6c4","#1d91c0","#225ea8","#253494","#081d58","#f2f0f7","#cbc9e2","#9e9ac8","#756bb1","#54278f")
tax_plot_paths$CAT_ORDER <-factor(tax_plot_paths$Category, levels = rev(category_list))
names(category_color) <- (category_list)

# taxonomic order and label
tax_order<-c("Dinoflagellate","Ciliate","MAST","Rhizaria","Haptophyta","Bacillariophyceae","Pelagophyceae","Chlorophyta")
tax_label<-c("a. Dinoflagellates","b. Ciliates","c. MAST","d. Rhizaria","e. Haptophytes","f. Diatoms","g. Pelagophytes","h. Chlorophytes")
tax_plot_paths$TAX_ORDER <- factor(tax_plot_paths$taxa, levels= tax_order, labels = tax_label)

bar_cat <- ggplot(tax_plot_paths, aes(y=(SUM_CPM), x=SAMPLE_ORDER))+
  geom_bar(stat = "identity", position="fill", aes(fill=(CAT_ORDER)), color="#525252", width = 0.70, size = 0.2)+
  coord_flip() +
  scale_fill_manual(values = (category_color))+
  facet_wrap(.~TAX_ORDER, ncol = 4) +
  theme_minimal()+
  theme(legend.position = "bottom", 
        axis.text.y = element_text(color="black", hjust=1, face = "bold"),
        axis.title = element_text(color = "black", face = "bold"),
        strip.text = element_text(h=0, face="bold"), 
        legend.title = element_blank())+
  guides(fill = guide_legend(reverse = TRUE)) +
  labs(x="", y="Relative abundance CPM")

options(repr.plot.width = 12, repr.plot.height = 6)
# svg("barplot-bytax-allfxn.svg", w = 12, h = 6)
bar_cat
# dev.off()

# load library ggpubr
library(ggpubr)

load("Normedbytax_annotated.RData", verbose = T)
head(tax_wcat[1:2,]); dim(tax_wcat)

rm <- c("butyryl", "other", "else", "additional breakdown", "por")

# input_wilcox_df <- sample_n(tax_wcat, 500000) %>% #Test dataset
input_wilcox_df <- tax_wcat %>%  #FULL
    filter(!(Target == "Phagotrophy-other")) %>% 
    filter(!(Category %in% rm)) %>% 
    unite(tax_cat, Category, taxa, sep = "_", remove = FALSE) %>% 
    data.frame

# hist(input_wilcox_df$mean_count)

wilcox_output <- compare_means(mean_count ~ sample, 
                     data = input_wilcox_df, 
                     method = "wilcox.test", 
                     group.by = "tax_cat",
                     p.adjust.method = "fdr")
# tmp <- head(wilcox_output)
table(wilcox_output$p.signif)

wilcox_output_wide <- wilcox_output %>% 
    separate(tax_cat, c("Category", "Taxa"), sep = "_") %>% 
    select(-.y., -p, -p.format, -method, -p.signif) %>% 
    mutate(p.sig = case_when(
        p.adj > 0.05 ~ "ns",
        (p.adj <= 0.05 & p.adj > 0.001) ~ "*",
        p.adj <= 0.001 ~ "***")) %>% 
    unite(COMPARE, Taxa, group1, group2, sep = "-") %>% 
    pivot_wider(id_cols = COMPARE, names_from = Category, 
                values_from = c(p.adj, p.sig),
                names_glue = "{Category}_{.value}") %>% 
    column_to_rownames("COMPARE") %>% 
    select(sort(colnames(.))) %>% 
    rownames_to_column("COMPARE") %>% 
    separate(COMPARE, c("Taxa", "SampleA", "SampleB"), sep = "-") %>% 
    data.frame
head(wilcox_output_wide)

dim(wilcox_output_wide)

write_delim(wilcox_output_wide, path = "Wilcox-paired-output.txt", delim = "\t")

# load("Normedbytax_annotated.RData", verbose = T)
head(tax_wcat[1:2,]); dim(tax_wcat)

top10_bytax <- tax_wcat %>% 
    filter(!(Target == "Phagotrophy-other")) %>% 
    filter(!(Category %in% rm)) %>% 
    # Sum CPM across location and depth to find top transcripts per taxa
    group_by(Taxonomy, taxa, Gene_identification, KO, mean_count) %>% 
    distinct() %>% # Unduplicate kegg entries that fall into more than 1 category
    summarise(SUM = sum(mean_count)) %>% 
    group_by(taxa) %>% 
    arrange(taxa, desc(SUM)) %>% 
    top_n(10, SUM) %>% 
    data.frame
head(top10_bytax)

write_delim(top10_bytax, path = "top10-bytax.txt", delim = "\t")
