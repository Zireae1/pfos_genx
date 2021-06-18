# Clean code to generate figures in the paper
rm(list=ls(all=TRUE))

library(data.table) # for fread function, which is the fastest way to load big tables
library(RColorBrewer) # for nice colour palettes
library(scales) # allows to adjust transparency of colors on plots
library(ggplot2)
library(ggsignif)
library(gplots)
library(vegan)
library(reshape2)
library(ggpubr)
library(ecodist)
library(ape)

# set default working directory. I also prefer to have the following folders within:
# src - all scripts
# data - all data tables
# graphs - folder to store saved plots
# output - folder to store all generated output such as tables with processed counts etc.:
setwd("/Users/Zireael/Desktop/16S_analysis/Clean_code/") # replace with your working directory

# load sample metadata
meta<-fread("data/Metadata.csv", stringsAsFactors = F, header = T, data.table = F)
meta$Chemical<-NA
meta$Concentration<-NA
for(i in 1:nrow(meta)){
  meta$Chemical[i]<-strsplit(meta$Treatment[i], "_")[[1]][1]
  meta$Concentration[i]<-strsplit(meta$Treatment[i], "_")[[1]][2]
}
meta<-meta[order(meta$`Sample ID`),]
meta<-meta[-which(meta$Barcode==""),]

### set colors      PFOS        GenX
main_palette<-c("#0073C2FF", "#EFC000FF")

# load microbial abundances
abundances<-fread("data/OTU_abundances.txt", stringsAsFactors = F, header = T, data.table = F)
which(is.na(abundances))

################################################################################################   
###                Fig 1C          calculate diversity indices and plot    
################################################################################################  

meta$Sh_diversity<-diversity(t(abundances[,3:ncol(abundances)]))

cairo_pdf(paste("graphs/Diversity_plots.pdf", sep=""), 
          width = 6, height = 4) # save plot
par(mfrow=c(1,2))

sub<-meta[25:48,]
sub<-sub[order(sub$Treatment),]

div.cm<-c(mean(sub$Sh_diversity[which(sub$Concentration=="C"&sub$Chemical=="GenX")]), 
          mean(sub$Sh_diversity[which(sub$Concentration=="10mg"&sub$Chemical=="GenX")]),
          mean(sub$Sh_diversity[which(sub$Concentration=="20mg"&sub$Chemical=="GenX")]),
          mean(sub$Sh_diversity[which(sub$Concentration=="100mg"&sub$Chemical=="GenX")]))
div.csd<-c(sd(sub$Sh_diversity[which(sub$Concentration=="C")]), 
           sd(sub$Sh_diversity[which(sub$Concentration=="10mg"&sub$Chemical=="GenX")]),
           sd(sub$Sh_diversity[which(sub$Concentration=="20mg"&sub$Chemical=="GenX")]),
           sd(sub$Sh_diversity[which(sub$Concentration=="100mg"&sub$Chemical=="GenX")]))


plot(c(1,3,4,5), div.cm,
     ylim=range(c(div.cm-div.csd, div.cm+div.csd)),
     xlim=c(1,5),
     xlab="Concentration, mg/L",
     pch=19, 
     ylab="Shannon diversity",
     #main="Colon",
     main="Small intestine",
     col=alpha(main_palette[2], 1.0), xaxt="n"
)
axis(1, at=c(1,2,3,4,5), labels=c("Ctrl", "5", "10", "20","100") )
# hack: we draw arrows but with very special "arrowheads"
arrows(c(1,3,4,5), div.cm-div.csd, c(1,3,4,5), div.cm+div.csd, length=0.05, angle=90, code=3, 
       alpha(main_palette[2], 1.0), lwd=2)
lines(c(1,3,4,5), div.cm, col=alpha(main_palette[2], 1.0), lwd=2)

div.cm<-c(mean(sub$Sh_diversity[which(sub$Concentration=="C"&sub$Chemical=="PFOS")]), 
          mean(sub$Sh_diversity[which(sub$Concentration=="5mg"&sub$Chemical=="PFOS")]),
          mean(sub$Sh_diversity[which(sub$Concentration=="10mg"&sub$Chemical=="PFOS")]),
          mean(sub$Sh_diversity[which(sub$Concentration=="20mg"&sub$Chemical=="PFOS")]))
div.csd<-c(sd(sub$Sh_diversity[which(sub$Concentration=="C"&sub$Chemical=="PFOS")]), 
           sd(sub$Sh_diversity[which(sub$Concentration=="5mg"&sub$Chemical=="PFOS")]),
           sd(sub$Sh_diversity[which(sub$Concentration=="10mg"&sub$Chemical=="PFOS")]),
           sd(sub$Sh_diversity[which(sub$Concentration=="20mg"&sub$Chemical=="PFOS")]))

points(c(1,2,3,4), div.cm, col=alpha(main_palette[1], 1.0), pch=16)

# hack: we draw arrows but with very special "arrowheads"
arrows(c(1,2,3,4), div.cm-div.csd, c(1,2,3,4), div.cm+div.csd, 
       length=0.05, angle=90, code=3, lwd=2,
       alpha(main_palette[1], 1.0))
lines(c(1,2,3,4), div.cm, col=alpha(main_palette[1], 1.0), lwd=2)

legend("bottomleft", #inset=c(-0.18,0), 
       bty="n", legend = c("GenX", "PFOS"), 
       col = c(main_palette[2], main_palette[1]), lty= 0, pch = 16)


sub<-meta[1:24,]
sub<-sub[order(sub$Treatment),]

div.cm<-c(mean(sub$Sh_diversity[which(sub$Concentration=="C"&sub$Chemical=="GenX")]), 
          mean(sub$Sh_diversity[which(sub$Concentration=="10mg"&sub$Chemical=="GenX")]),
          mean(sub$Sh_diversity[which(sub$Concentration=="20mg"&sub$Chemical=="GenX")]),
          mean(sub$Sh_diversity[which(sub$Concentration=="100mg"&sub$Chemical=="GenX")]))
div.csd<-c(sd(sub$Sh_diversity[which(sub$Concentration=="C")]), 
           sd(sub$Sh_diversity[which(sub$Concentration=="10mg"&sub$Chemical=="GenX")]),
           sd(sub$Sh_diversity[which(sub$Concentration=="20mg"&sub$Chemical=="GenX")]),
           sd(sub$Sh_diversity[which(sub$Concentration=="100mg"&sub$Chemical=="GenX")]))


plot(c(1,3,4,5), div.cm,
     ylim=range(c(div.cm-div.csd, div.cm+div.csd)),
     xlim=c(1,5),
     xlab="Concentration, mg/L",
     pch=19, 
     ylab="Shannon diversity",
     main="Colon",
     #main="Small intestine",
     col=alpha(main_palette[2], 1.0), xaxt="n"
)
axis(1, at=c(1,2,3,4,5), labels=c("Ctrl", "5", "10", "20","100") )
# hack: we draw arrows but with very special "arrowheads"
arrows(c(1,3,4,5), div.cm-div.csd, c(1,3,4,5), div.cm+div.csd, length=0.05, angle=90, code=3, 
       alpha(main_palette[2], 1.0), lwd=2)
lines(c(1,3,4,5), div.cm, col=alpha(main_palette[2], 1.0), lwd=2)

div.cm<-c(mean(sub$Sh_diversity[which(sub$Concentration=="C"&sub$Chemical=="PFOS")]), 
          mean(sub$Sh_diversity[which(sub$Concentration=="5mg"&sub$Chemical=="PFOS")]),
          mean(sub$Sh_diversity[which(sub$Concentration=="10mg"&sub$Chemical=="PFOS")]),
          mean(sub$Sh_diversity[which(sub$Concentration=="20mg"&sub$Chemical=="PFOS")]))
div.csd<-c(sd(sub$Sh_diversity[which(sub$Concentration=="C"&sub$Chemical=="PFOS")]), 
           sd(sub$Sh_diversity[which(sub$Concentration=="5mg"&sub$Chemical=="PFOS")]),
           sd(sub$Sh_diversity[which(sub$Concentration=="10mg"&sub$Chemical=="PFOS")]),
           sd(sub$Sh_diversity[which(sub$Concentration=="20mg"&sub$Chemical=="PFOS")]))

points(c(1,2,3,4), div.cm, col=alpha(main_palette[1], 1.0), pch=16)

# hack: we draw arrows but with very special "arrowheads"
arrows(c(1,2,3,4), div.cm-div.csd, c(1,2,3,4), div.cm+div.csd, 
       length=0.05, angle=90, code=3, lwd=2,
       alpha(main_palette[1], 1.0))
lines(c(1,2,3,4), div.cm, col=alpha(main_palette[1], 1.0), lwd=2)

legend("bottomleft", #inset=c(-0.18,0), 
       bty="n", legend = c("GenX", "PFOS"), 
       col = c(main_palette[2], main_palette[1]), lty= 0, pch = 16)

dev.off()

length(which(rowSums(abundances[,3:ncol(abundances)])<0.001))
abundances<-abundances[-which(rowSums(abundances[,3:ncol(abundances)])<0.001),]

################################################################################################   
###                Fig 1B          PCoA for Bray-Curtis distance   
################################################################################################  

#### create two separate pcoa plots for SI and Colon
######   SI
subset_abund<-abundances[,c(1,which(colnames(abundances)%in%meta$`Sample ID`[which(meta$`Body site`=="SI")]))]
subset_abund<-subset_abund[which(rowSums(subset_abund[3:ncol(subset_abund)])>0.005),]
dm<-as.matrix(bcdist(t(subset_abund[,2:ncol(subset_abund)])))
pc<-pcoa(dm)

# Axis 1&2
#coord<-as.data.frame(pc$vectors[,1:2])
coord<-as.data.frame(cbind(-pc$vectors[,1],pc$vectors[,2]))
colnames(coord)<-c("Coord 1", "Coord 2")

perc<-pc$values$Cum_corr_eig[1:2]
perc[2]<-perc[2]-perc[1]
perc<-round(perc*100, digits = 2)


colnames(coord)<-c(paste("Axis 1 ", "(", perc[1], "%)", sep=""), 
                   paste("Axis 2 ", "(", perc[2], "%)", sep=""))

coord$Chemical<-NA
coord$Concentration<-NA
coord$Part<-NA
for(i in 1:nrow(coord)){
  #i<-1
  coord$Chemical[i]<-meta$Chemical[which(meta$`Sample ID`==rownames(coord)[i])]
  coord$Concentration[i]<-meta$Concentration[which(meta$`Sample ID`==rownames(coord)[i])]
  coord$Part[i]<-strsplit(rownames(coord)[i], "_")[[1]][1]
}

coord$color<-NA
coord$pch<-NA
coord$size<-NA

#rownames(coord)<-gsub("Colon_", "", rownames(coord))
rownames(coord)<-gsub("SI_", "", rownames(coord))

coord$color[which(coord$Chemical=="GenX")]<-main_palette[2]
coord$color[which(coord$Chemical=="PFOS")]<-main_palette[1]
coord$pch[which(coord$Part=="Colon")]<-21
coord$pch[which(coord$Part=="SI")]<-24
coord$color[which(coord$Concentration=="C")]<-alpha(coord$color[which(coord$Concentration=="C")],0.25)
coord$color[which(coord$Concentration=="5mg")]<-alpha(coord$color[which(coord$Concentration=="5mg")],0.5)
coord$color[which(coord$Concentration=="10mg"&coord$Chemical=="PFOS")]<-alpha(coord$color[which(coord$Concentration=="10mg"&coord$Chemical=="PFOS")],0.75)
coord$color[which(coord$Concentration=="20mg"&coord$Chemical=="PFOS")]<-alpha(coord$color[which(coord$Concentration=="20mg"&coord$Chemical=="PFOS")],1)
coord$color[which(coord$Concentration=="10mg"&coord$Chemical=="GenX")]<-alpha(coord$color[which(coord$Concentration=="10mg"&coord$Chemical=="GenX")],0.5)
coord$color[which(coord$Concentration=="20mg"&coord$Chemical=="GenX")]<-alpha(coord$color[which(coord$Concentration=="20mg"&coord$Chemical=="GenX")],0.75)
coord$color[which(coord$Concentration=="100mg")]<-alpha(coord$color[which(coord$Concentration=="100mg")],1)
coord$size[which(coord$Concentration=="C")]<-1
coord$size[which(coord$Concentration=="5mg")]<-1.25
coord$size[which(coord$Concentration=="10mg")]<-1.5
coord$size[which(coord$Concentration=="20mg")]<-1.75
coord$size[which(coord$Concentration=="100mg")]<-2


cairo_pdf(paste("graphs/pcoa_12_SI.pdf", sep=""), 
          width = 6.5, height = 6.5) # save plot

plot(coord[,1:2], col="black", bg=coord$color, pch=coord$pch, cex=coord$size, ylim=c(-.5,0.4)
)
text(coord[,1:2], rownames(coord), pos = 3)
legend("topright", inset=c(0,0), title = c("GenX"),
       bty="n", legend = c("Control","10mg/L", "20mg/L", "100mg/L"), col="black",
       pt.bg = c(alpha(main_palette[2],0.25), alpha(main_palette[2],0.5),alpha(main_palette[2],0.75), alpha(main_palette[2],1)), 
       lty= 0, pch = unique(coord$pch), pt.cex=c(1,1.5,1.75,2))
legend("bottomright", inset=c(0,0), title = c("PFOS"),
       bty="n", legend = c("Control", "5mg/L", "10mg/L", "20mg/L"), col="black",
       pt.bg = c(alpha(main_palette[1],0.25),alpha(main_palette[1],0.5),
                 alpha(main_palette[1],0.75), alpha(main_palette[1],1)), 
       lty= 0, pch = unique(coord$pch), pt.cex=c(1,1.25,1.5,1.75))



dev.off()

######   Colon
subset_abund<-abundances[,c(1,which(colnames(abundances)%in%meta$`Sample ID`[which(meta$`Body site`=="Colon")]))]
subset_abund<-subset_abund[which(rowSums(subset_abund[3:ncol(subset_abund)])>0.005),]
dm<-as.matrix(bcdist(t(subset_abund[,2:ncol(subset_abund)])))
pc<-pcoa(dm)

# Axis 1&2
coord<-as.data.frame(pc$vectors[,1:2])
#coord<-as.data.frame(cbind(-pc$vectors[,1],pc$vectors[,2]))
colnames(coord)<-c("Coord 1", "Coord 2")

perc<-pc$values$Cum_corr_eig[1:2]
perc[2]<-perc[2]-perc[1]
perc<-round(perc*100, digits = 2)


colnames(coord)<-c(paste("Axis 1 ", "(", perc[1], "%)", sep=""), 
                   paste("Axis 2 ", "(", perc[2], "%)", sep=""))

coord$Chemical<-NA
coord$Concentration<-NA
coord$Part<-NA
for(i in 1:nrow(coord)){
  #i<-1
  coord$Chemical[i]<-meta$Chemical[which(meta$`Sample ID`==rownames(coord)[i])]
  coord$Concentration[i]<-meta$Concentration[which(meta$`Sample ID`==rownames(coord)[i])]
  coord$Part[i]<-strsplit(rownames(coord)[i], "_")[[1]][1]
}

coord$color<-NA
coord$pch<-NA
coord$size<-NA

rownames(coord)<-gsub("Colon_", "", rownames(coord))
#rownames(coord)<-gsub("SI_", "", rownames(coord))

coord$color[which(coord$Chemical=="GenX")]<-main_palette[2]
coord$color[which(coord$Chemical=="PFOS")]<-main_palette[1]
coord$pch[which(coord$Part=="Colon")]<-21
coord$pch[which(coord$Part=="SI")]<-24
coord$color[which(coord$Concentration=="C")]<-alpha(coord$color[which(coord$Concentration=="C")],0.25)
coord$color[which(coord$Concentration=="5mg")]<-alpha(coord$color[which(coord$Concentration=="5mg")],0.5)
coord$color[which(coord$Concentration=="10mg"&coord$Chemical=="PFOS")]<-alpha(coord$color[which(coord$Concentration=="10mg"&coord$Chemical=="PFOS")],0.75)
coord$color[which(coord$Concentration=="20mg"&coord$Chemical=="PFOS")]<-alpha(coord$color[which(coord$Concentration=="20mg"&coord$Chemical=="PFOS")],1)
coord$color[which(coord$Concentration=="10mg"&coord$Chemical=="GenX")]<-alpha(coord$color[which(coord$Concentration=="10mg"&coord$Chemical=="GenX")],0.5)
coord$color[which(coord$Concentration=="20mg"&coord$Chemical=="GenX")]<-alpha(coord$color[which(coord$Concentration=="20mg"&coord$Chemical=="GenX")],0.75)
coord$color[which(coord$Concentration=="100mg")]<-alpha(coord$color[which(coord$Concentration=="100mg")],1)
coord$size[which(coord$Concentration=="C")]<-1
coord$size[which(coord$Concentration=="5mg")]<-1.25
coord$size[which(coord$Concentration=="10mg")]<-1.5
coord$size[which(coord$Concentration=="20mg")]<-1.75
coord$size[which(coord$Concentration=="100mg")]<-2


cairo_pdf(paste("graphs/pcoa_12_Colon.pdf", sep=""), 
          width = 6.5, height = 6.5) # save plot

plot(coord[,1:2], col="black", bg=coord$color, pch=coord$pch, cex=coord$size, ylim=c(-.5,0.4)
)
text(coord[,1:2], rownames(coord), pos = 3)
legend("topright", inset=c(0,0), title = c("GenX"),
       bty="n", legend = c("Control","10mg/L", "20mg/L", "100mg/L"), col="black",
       pt.bg = c(alpha(main_palette[2],0.25), alpha(main_palette[2],0.5),alpha(main_palette[2],0.75), alpha(main_palette[2],1)), 
       lty= 0, pch = unique(coord$pch), pt.cex=c(1,1.5,1.75,2))
legend("bottomright", inset=c(0,0), title = c("PFOS"),
       bty="n", legend = c("Control", "5mg/L", "10mg/L", "20mg/L"), col="black",
       pt.bg = c(alpha(main_palette[1],0.25),alpha(main_palette[1],0.5),
                 alpha(main_palette[1],0.75), alpha(main_palette[1],1)), 
       lty= 0, pch = unique(coord$pch), pt.cex=c(1,1.25,1.5,1.75))



dev.off()

##### check that controls are closer to each other than to others:
ctrl_ids<-which(colnames(dm)%in%c("SI_C1", "SI_C2", "SI_C3","SI_C4", "SI_C5", "SI_C6"))
ctrl_ids<-which(colnames(dm)%in%c("Colon_C1", "Colon_C2", "Colon_C3","Colon_C4", "Colon_C5", "Colon_C6"))

dist_cvsc<-as.numeric(dm[ctrl_ids, ctrl_ids])
dist_cvsc<-dist_cvsc[-which(dist_cvsc==0)]
dist_cvsa<-as.numeric(dm[-ctrl_ids, ctrl_ids])

t.test(dist_cvsc,dist_cvsa)
ks.test(dist_cvsc,dist_cvsa)


################################################################################################   
###                   Heatmap top genera separate colon/SI       (Fig 2A)
################################################################################################
library(ComplexHeatmap)
ss_abund<-fread("data/RDP_80_genus_abundance.txt", data.table = F, stringsAsFactors = F, header = T)
ss_abund<-ss_abund[which(rowSums(ss_abund[,2:ncol(ss_abund)])>0.005),]

meta$Concentration_num<-gsub("mg","",meta$Concentration)
meta$Concentration_num[which(meta$Concentration_num=="C")]<-0
meta$Concentration_num<-as.numeric(meta$Concentration_num)


my_palette<-c(colorRampPalette(rev(brewer.pal(11,"RdBu")))(300))
df<-as.matrix(ss_abund[,2:ncol(ss_abund)])

rownames(df)<-ss_abund$rdp

group_col<-character(ncol(df))
group_col[which(colnames(df)%in%meta$`Sample ID`[which(meta$Chemical=="PFOS")])]<-main_palette[1]
group_col[which(colnames(df)%in%meta$`Sample ID`[which(meta$Chemical=="GenX")])]<-main_palette[2]

split <- factor(meta$`Body site`, levels = c("Colon", "SI"))


ids<-which(colnames(df)%in%meta$`Sample ID`[which(meta$Concentration_num==0)])
group_col[ids]<-alpha(group_col[ids],0.25)
ids<-which(colnames(df)%in%meta$`Sample ID`[which(meta$Concentration_num==5)])
group_col[ids]<-alpha(group_col[ids],0.5)
ids<-which(colnames(df)%in%meta$`Sample ID`[which(meta$Concentration_num==10&meta$Chemical=="GenX")])
group_col[ids]<-alpha(group_col[ids],0.5)
ids<-which(colnames(df)%in%meta$`Sample ID`[which(meta$Concentration_num==10&meta$Chemical=="PFOS")])
group_col[ids]<-alpha(group_col[ids],0.75)
ids<-which(colnames(df)%in%meta$`Sample ID`[which(meta$Concentration_num==20&meta$Chemical=="GenX")])
group_col[ids]<-alpha(group_col[ids],0.75)
ids<-which(colnames(df)%in%meta$`Sample ID`[which(meta$Concentration_num==20&meta$Chemical=="PFOS")])
group_col[ids]<-alpha(group_col[ids],1)

colnames(df)<-gsub("Colon_", "", colnames(df))
colnames(df)<-gsub("SI_", "", colnames(df))

#group_col<-as.data.frame(group_col)
names(group_col)<-meta$`Sample ID`

cairo_pdf(paste("graphs/heatmap_top_genera_split.pdf", sep=""), 
          width = 13, height = 9.5) # save plot

ncol(df)

column_ha = HeatmapAnnotation(treatment=as.vector(names(group_col)), col = list(treatment=group_col))
Heatmap(log10(df+10^-5), column_split = split, col = my_palette, top_annotation = column_ha)

dev.off()

################################################################################################   
###                   Heatmap top family separate colon/SI       (Fig S1)
################################################################################################

ss_abund<-fread("data/RDP_80_family_abundance.txt", data.table = F, stringsAsFactors = F, header = T)
ss_abund<-ss_abund[which(rowSums(ss_abund[,2:ncol(ss_abund)])>0.005),]

meta$Concentration_num<-gsub("mg","",meta$Concentration)
meta$Concentration_num[which(meta$Concentration_num=="C")]<-0
meta$Concentration_num<-as.numeric(meta$Concentration_num)


my_palette<-c(colorRampPalette(rev(brewer.pal(11,"RdBu")))(300))
df<-as.matrix(ss_abund[,2:ncol(ss_abund)])

rownames(df)<-ss_abund$rdp

group_col<-character(ncol(df))
group_col[which(colnames(df)%in%meta$`Sample ID`[which(meta$Chemical=="PFOS")])]<-main_palette[1]
group_col[which(colnames(df)%in%meta$`Sample ID`[which(meta$Chemical=="GenX")])]<-main_palette[2]

split <- factor(meta$`Body site`, levels = c("Colon", "SI"))


ids<-which(colnames(df)%in%meta$`Sample ID`[which(meta$Concentration_num==0)])
group_col[ids]<-alpha(group_col[ids],0.25)
ids<-which(colnames(df)%in%meta$`Sample ID`[which(meta$Concentration_num==5)])
group_col[ids]<-alpha(group_col[ids],0.5)
ids<-which(colnames(df)%in%meta$`Sample ID`[which(meta$Concentration_num==10&meta$Chemical=="GenX")])
group_col[ids]<-alpha(group_col[ids],0.5)
ids<-which(colnames(df)%in%meta$`Sample ID`[which(meta$Concentration_num==10&meta$Chemical=="PFOS")])
group_col[ids]<-alpha(group_col[ids],0.75)
ids<-which(colnames(df)%in%meta$`Sample ID`[which(meta$Concentration_num==20&meta$Chemical=="GenX")])
group_col[ids]<-alpha(group_col[ids],0.75)
ids<-which(colnames(df)%in%meta$`Sample ID`[which(meta$Concentration_num==20&meta$Chemical=="PFOS")])
group_col[ids]<-alpha(group_col[ids],1)

colnames(df)<-gsub("Colon_", "", colnames(df))
colnames(df)<-gsub("SI_", "", colnames(df))

#group_col<-as.data.frame(group_col)
names(group_col)<-meta$`Sample ID`

cairo_pdf(paste("graphs/heatmap_top_family_split.pdf", sep=""), 
          width = 15, height = 10) # save plot

ncol(df)

column_ha = HeatmapAnnotation(sample=as.vector(names(group_col)), col = list(sample=group_col))
Heatmap(log10(df+10^-5), column_split = split, col = my_palette, top_annotation = column_ha)

dev.off()


################################################################################################   
###                   Plots for individual genera       (Fig S2-5)
################################################################################################
# load microbial abundances
abundance<-fread("data/RDP_80_genus_abundance.txt", stringsAsFactors = F, header = T, data.table = F)
colSums(abundance[,2:ncol(abundance)])
abundance<-abundance[which(rowSums(abundance[,2:ncol(abundance)])>0.005),]
rownames(abundance)<-abundance$rdp
abundance<-abundance[,-1]
rownames(abundance)[which(rownames(abundance)=="Escherichia/Shigella")]<-"Escherichia"

# A function to add arrows on the chart
error.bar <- function(x, y, upper, lower=0, length=0.1, color){
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, col=color)
}

df<-t(abundance)
for(i in 1:ncol(df)){
  #i<-1
  bac<-df[,i]
  for (bs in unique(meta$`Body site`)){
    #bs<-"Colon"
    it<-2
    for(ch in unique(meta$Chemical)){
      
      #ch<-"PFOS"
      means<-vector()
      sds<-vector()
      k<-1
      for(tr in unique(meta$Treatment[which(meta$Chemical==ch)])){
        #tr<-"GenX_10mg"
        means[k]<-mean(bac[which(names(bac)%in%meta$`Sample ID`[which(meta$Treatment==tr&meta$`Body site`==bs&meta$Chemical==ch)])])
        sds[k]<-sd(bac[which(names(bac)%in%meta$`Sample ID`[which(meta$Treatment==tr&meta$`Body site`==bs&meta$Chemical==ch)])])
        k<-k+1
      }
      names(means)<-unique(meta$Treatment[which(meta$Chemical==ch)])
      names(sds)<-unique(meta$Treatment[which(meta$Chemical==ch)])
      print(colnames(df)[i])
      print(sds)
      means<-means[c(3,1,2,4)]
      sds<-sds[c(3,1,2,4)]
      
      if(ch=="PFOS"){
        names(means)<-c("Ctrl", "5mg", "10mg", "20mg")
      }else{
        names(means)<-c("Ctrl", "10mg", "20mg", "100mg")
      }
      
      cairo_pdf(paste("graphs/",ch,"_",bs,"_",colnames(df)[i],".pdf", sep=""), 
                width = 3.8, height = 4) # save plot
      
      ze_barplot<-barplot(means*100, ylim=c(0,max((means+sds)*100)), 
                          xlab = c("Concentration, mg/L"), col=alpha(main_palette[it], 0.8), ylab="Relative abundance, %",
                          main=colnames(df)[i])
      
      error.bar(ze_barplot, means*100, sds*100, color="black")
      dev.off()
      it<-it-1
    }
  }
  print(i)
}

################################################################################################   
###                   Plots for phylum and B/F ratio       (Fig S6)
################################################################################################
abundance<-fread("data/RDP_50_phylum_abundance.txt", stringsAsFactors = F, header = T, data.table = F)
rownames(abundance)<-abundance$rdp
abundance<-abundance[,-1]
rownames(abundance)[which(rownames(abundance)=="Cyanobacteria/Chloroplast")]<-"Cyanobacteria"


for (bs in unique(meta$`Body site`)){
  #bs<-"Colon"
  ### PFOS
  sub<-abundance[which(colnames(abundance)%in%meta$`Sample ID`[which(meta$Chemical=="PFOS"&meta$`Body site`==bs)])]
  colnames(sub)
  p.div.cm<-cbind(rowMeans(sub[,7:9]), 
                rowMeans(sub[,1:3]),
                rowMeans(sub[,4:6]),
                rowMeans(sub[,10:12]))
  p.div.csd<-cbind(apply(sub[,7:9],1, sd), 
                 apply(sub[,1:3],1, sd),
                 apply(sub[,4:6],1, sd),
                 apply(sub[,10:12],1, sd))
  ### GenX
  sub<-abundance[which(colnames(abundance)%in%meta$`Sample ID`[which(meta$Chemical=="GenX"&meta$`Body site`==bs)])]
  colnames(sub)
  g.div.cm<-cbind(rowMeans(sub[,7:9]), 
                  rowMeans(sub[,1:3]),
                  rowMeans(sub[,4:6]),
                  rowMeans(sub[,10:12]))
  g.div.csd<-cbind(apply(sub[,7:9],1, sd), 
                   apply(sub[,1:3],1, sd),
                   apply(sub[,4:6],1, sd),
                   apply(sub[,10:12],1, sd))
  for(i in 1:nrow(abundance)){
    #i<-1
    cairo_pdf(paste("graphs/Phylum_",bs,"_",rownames(abundance)[i],".pdf", sep=""), 
              width = 3, height = 4) # save plot
    plot(c(1,2,3,4), p.div.cm[i,]*100,
         #ylim=range(c(div.cm[i,]-div.csd[i,], div.cm[i,]+div.csd[i,])),
         xlim=c(1,5),
         xlab="Concentration, mg/L",
         pch=19, 
         ylab="Relative abundance, %",
         #main="Colon",
         ylim=c(0,max(c(p.div.cm[i,]+p.div.csd[i,],g.div.cm[i,]+g.div.csd[i,]))*100),
         main=rownames(sub)[i],
         col=alpha(main_palette[1], 1), xaxt="n"
    )
    axis(1, at=c(1,2,3,4,5), labels=c("Ctrl", "5", "10", "20","100") )
    # hack: we draw arrows but with very special "arrowheads"
    arrows(c(1,2,3,4), (p.div.cm[i,]-p.div.csd[i,])*100, c(1,2,3,4), (p.div.cm[i,]+p.div.csd[i,])*100, 
           length=0.05, angle=90, code=3, 
           alpha(main_palette[1], 1), lwd=2)
    lines(c(1,2,3,4), p.div.cm[i,]*100, col=alpha(main_palette[1], 1), lwd=2)
    
    
    points(c(1,3,4,5), g.div.cm[i,]*100, col=alpha(main_palette[2], 1), pch=16)
    
    # hack: we draw arrows but with very special "arrowheads"
    arrows(c(1,3,4,5), (g.div.cm[i,]-g.div.csd[i,])*100, c(1,3,4,5), (g.div.cm[i,]+g.div.csd[i,])*100, 
           length=0.05, angle=90, code=3, lwd=2,
           alpha(main_palette[2], 1))
    lines(c(1,3,4,5), g.div.cm[i,]*100, col=alpha(main_palette[2], 1), lwd=2)
    
    legend("topleft", #inset=c(-0.18,0), 
           bty="n", legend = c("GenX", "PFOS"), 
           col = c(main_palette[2], main_palette[1]), lty= 0, pch = 16)
    
    dev.off()
  }
  
}

#### look at B/F ratio
sub.abund<-abundance[,1:24] # colon
sub.abund<-abundance[,25:48] # SI

ratio<-sub.abund[2,]/sub.abund[6,]

### create a plot:
pfos<-ratio[which(colnames(ratio)%in%meta$`Sample ID`[which(meta$Chemical=="GenX")])]
#GenX
colnames(pfos)
div.cm<-cbind(rowMeans(pfos[,7:9]), 
              rowMeans(pfos[,1:3]),
              rowMeans(pfos[,4:6]),
              rowMeans(pfos[,10:12]))
div.csd<-cbind(apply(pfos[,7:9],1, sd), 
               apply(pfos[,1:3],1, sd),
               apply(pfos[,4:6],1, sd),
               apply(pfos[,10:12],1, sd))
cairo_pdf(paste("graphs/SI_B.F_ratio.pdf", sep=""), 
          width = 3, height = 4) # save plot


plot(c(1,3,4,5), div.cm,
     #ylim=range(c(div.cm[i,]-div.csd[i,], div.cm[i,]+div.csd[i,])),
     xlim=c(1,5),
     xlab="Concentration, mg/L",
     pch=19, 
     ylab="Bacteroides/Firmicutes ratio",
     #main="Colon",
     ylim=c(0,2),
     main="Colon",
     col=alpha(main_palette[2], 1), xaxt="n"
)
axis(1, at=c(1,2,3,4,5), labels=c("Ctrl", "5", "10", "20","100") )
# hack: we draw arrows but with very special "arrowheads"
arrows(c(1,3,4,5), div.cm-div.csd, c(1,3,4,5), div.cm+div.csd, 
       length=0.05, angle=90, code=3, 
       alpha(main_palette[2], 1), lwd=2)
lines(c(1,3,4,5), div.cm, col=alpha(main_palette[2], 1), lwd=2)
#PFOS
pfos<-ratio[which(colnames(ratio)%in%meta$`Sample ID`[which(meta$Chemical=="PFOS")])]
colnames(pfos)
div.cm<-cbind(rowMeans(pfos[,7:9]), 
              rowMeans(pfos[,1:3]),
              rowMeans(pfos[,4:6]),
              rowMeans(pfos[,10:12]))
div.csd<-cbind(apply(pfos[,7:9],1, sd), 
               apply(pfos[,1:3],1, sd),
               apply(pfos[,4:6],1, sd),
               apply(pfos[,10:12],1, sd))

points(c(1,2,3,4), div.cm, col=alpha(main_palette[1], 1), pch=16)

# hack: we draw arrows but with very special "arrowheads"
arrows(c(1,2,3,4), div.cm-div.csd, c(1,2,3,4), div.cm+div.csd, 
       length=0.05, angle=90, code=3, lwd=2,
       alpha(main_palette[1], 1))
lines(c(1,2,3,4), div.cm, col=alpha(main_palette[1], 1), lwd=2)

legend("topleft", #inset=c(-0.18,0), 
       bty="n", legend = c("GenX", "PFOS"), 
       col = c(main_palette[2], main_palette[1]), lty= 0, pch = 16)

abline(h=1, lty=2, col="grey")

dev.off()


################################################################################################   
###                   Plots for regression analysis      (Fig 2B-C)
################################################################################################
abundance<-fread("data/RDP_80_genus_abundance.txt", stringsAsFactors = F, header = T, data.table = F)
colSums(abundance[,2:ncol(abundance)])
abundance<-abundance[which(rowSums(abundance[,2:ncol(abundance)])>0.005),]
rownames(abundance)<-abundance$rdp
abundance<-abundance[,-1]
rownames(abundance)[which(rownames(abundance)=="Escherichia/Shigella")]<-"Escherichia"

### create regression plots for each metabolite, pick best:
# meta$OtherID<-gsub("SI_","",meta$`Sample ID`)
# meta$OtherID<-gsub("Colon_","",meta$`Sample ID`)

sub.meta<-meta[which(meta$`Body site`=="Colon"&meta$Chemical=="GenX"),]
sub.abund<-abundance[,which(colnames(abundance)%in%sub.meta$`Sample ID`)]
sub.abund<-sub.abund[which(rowSums(sub.abund)>0),]
sub.meta<-sub.meta[which(sub.meta$`Sample ID`%in%colnames(sub.abund)),]
sub.meta<-sub.meta[order(sub.meta$`Sample ID`),]

lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}


score<-vector()
pval<-vector()
res<-data.frame()
sub.abund<-sub.abund+0.00005
for (i in 1:nrow(sub.abund)){
  #i<-1
  ab<-log2(as.numeric(sub.abund[i,])*100)
  conc<-sub.meta$Concentration_num/max(sub.meta$Concentration_num)
  reg<-lm(ab~conc)
  # grab p-value
  pval[i]<-lmp(reg)
  # grab multiplier
  score[i]<-reg$coefficients[2]
  # save plots on the way
  # plot(sub.meta$Concentration_num,sub.abund[i,], pch=21, col=main_palette[2], 
  #      xlab="Concentration, mg/L", xaxt="n")
  # axis(1, at=c(1,2,3,4,5), labels=c("Ctrl", "5", "10", "20","100") )
  
}
##### create log ratio barplots

diff_si_genx<-cbind(rownames(sub.abund)[which(pval<0.05)],score[which(pval<0.05)],main_palette[2])
diff_si_pfos<-cbind(rownames(sub.abund)[which(pval<0.05)],score[which(pval<0.05)],main_palette[1])

diff_colon_pfos<-cbind(rownames(sub.abund)[which(pval<0.05)],score[which(pval<0.05)],main_palette[1])

diff_colon_genx<-cbind(rownames(sub.abund)[which(pval<0.05)],score[which(pval<0.05)],main_palette[2])

diff_si_pfos<-as.data.frame(diff_si_pfos)
#diff_si_pfos$V3<-"royal blue"
diff_si_genx<-as.data.frame(diff_si_genx)
all_si<-rbind(diff_si_pfos,diff_si_genx)
#all_si$V2<-1/as.numeric(all_si$V2)  ### reorder to represent size of effect
all_si<-all_si[order(as.numeric(all_si$V2)),]

cairo_pdf('graphs/Barplot_genus_SI_0.05.pdf', width = 6, height = 4) # save plot
par(mai=c(1,3,1,1))
barplot(as.numeric(all_si[,2]), col=all_si[,3], horiz=TRUE, 
        names.arg = all_si[,1], las=1, xlim = c(-10,10), 
        xlab = "Log2 Fold Change")

dev.off()

all_colon<-rbind(diff_colon_pfos,diff_colon_genx)
all_colon<-as.data.frame(all_colon)
#all_si$V2<-1/as.numeric(all_si$V2)  ### reorder to represent size of effect
all_colon<-all_colon[order(as.numeric(all_colon$V2)),]

cairo_pdf('graphs/Barplot_genus_Colon_0.05.pdf', width = 6, height = 4) # save plot
par(mai=c(1,3,1,1))
barplot(as.numeric(all_colon[,2]), col=all_colon[,3], horiz=TRUE, 
        names.arg = all_colon[,1], las=1, xlim = c(-10,10), 
        xlab = "Log2 Fold Change")

dev.off()


################################################################################################   
###                   Plots for top PICRUST pathways     (Fig 2D-E) finish
################################################################################################

# load pathway abundances and descriptions
pathw<-fread("data/picrust_path_abun_unstrat.tsv", stringsAsFactors = F, header = T, data.table = F)
desc<-fread("data/picrust_path_abun_unstrat_descrip.tsv", stringsAsFactors = F, header = T, data.table = F)[,1:2]

# hist(rowSums(log10(pathw[,2:ncol(pathw)]+10^-5)))
# barplot(colSums(pathw[,2:ncol(pathw)]))

rownames(pathw)<-pathw$pathway
pathw<-pathw[,-1]
path_reg<-function(pathw, meta, bs, chem){
  sub.meta<-meta[which(meta$`Body site`==bs&meta$Chemical==chem),]
  #sub.meta<-meta[which(meta$`Body site`=="Colon"&meta$Chemical=="GenX"),]
  
  #sub.meta<-meta[which(meta$`Body site`=="SI"&meta$Chemical=="PFOS"),]
  #sub.meta<-meta[which(meta$`Body site`=="SI"&meta$Chemical=="GenX"),]
  
  
  sub.pathw<-pathw[,which(colnames(pathw)%in%sub.meta$`Sample ID`)]
  sub.pathw<-sub.pathw[which(rowSums(sub.pathw)>0),]
  sub.meta<-sub.meta[which(sub.meta$`Sample ID`%in%colnames(sub.pathw)),]
  sub.meta<-sub.meta[order(sub.meta$`Sample ID`),]
  
  ### let's first subtract control??? - too many negative intersepts
  #sub.pathw<-t(apply(sub.pathw,1,function(x) x- mean(x[7:9])))
  
  score<-vector()
  pval<-vector()
  fc_effect<-vector()
  fc_means<-vector()
  sub.pathw<-sub.pathw+0.00005
  for (i in 1:nrow(sub.pathw)){
    #i<-3
    #reg<-lm(sub.meta$Concentration_num/max(sub.meta$Concentration_num)~log10(as.numeric(sub.pathw[i,])))
    #reg<-lm(sub.meta$Concentration_num/max(sub.meta$Concentration_num)~as.numeric(sub.pathw[i,]))
    renorm<-sub.meta$Concentration_num/max(sub.meta$Concentration_num)
    
    reg<-lm(as.numeric(sub.pathw[i,])~renorm)
    #reg<-lm(log2(as.numeric(sub.pathw[i,]))~renorm)
    # grab p-value
    pval[i]<-lmp(reg)
    # grab multiplier
    score[i]<-reg$coefficients[2]
    print(reg$coefficients[1])
    fc_effect[i]<-(reg$coefficients[2]+reg$coefficients[1])/abs(reg$coefficients[1])
    fc_means[i]<-mean(as.numeric(sub.pathw[i,10:12]))/mean(as.numeric(sub.pathw[i,7:9])) # max by ctrl
  }
}


diff_colon_pfos<-cbind(rownames(sub.pathw)[which(pval<0.01)],
                       score[which(pval<0.01)],fc_effect[which(pval<0.01)],log2(fc_means[which(pval<0.01)]),"royal blue")
#diff_colon_pfos[,3]<-"royal blue"
diff_colon_genx<-cbind(rownames(sub.pathw)[which(pval<0.01)],
                       score[which(pval<0.01)],fc_effect[which(pval<0.01)],log2(fc_means[which(pval<0.01)]),"red")

all_colon<-rbind(as.data.frame(diff_colon_pfos), as.data.frame(diff_colon_genx))
all_colon<-all_colon[order(as.numeric(all_colon$V2)),]
all_colon$V1[which(duplicated(all_colon$V1))]

#all_colon<-all_colon[order(as.numeric(all_colon$V3)),]
all_colon<-all_colon[order(as.numeric(all_colon$V2)),]

all_colon$V3<-log2(as.numeric(all_colon$V3))

### add descriptions:
for(i in 1:nrow(all_colon)){
  all_colon$V1[i]<-desc$description[which(desc$pathway==all_colon$V1[i])]
}
#which(all_colon$V1%in%desc$pathway)
#all_colon<-all_colon[c(which(log2(as.numeric(all_colon[,3]))<(-1)),which(log2(as.numeric(all_colon[,3]))>1)),]
all_colon<-all_colon[c(which(as.numeric(all_colon[,2])<(-1)),which(as.numeric(all_colon[,2])>1)),]


#which(log2(as.numeric(all_colon[,3]))<(-1))
all_colon<-all_colon[order(as.numeric(all_colon$V3)),]

cairo_pdf('graphs/Barplot_Metacyc_pathway_SI_0.01_fc2_reg_all_new.pdf', width = 9, height = 8) # save plot

par(mai=c(1,6,1,1))
# barplot(log2(as.numeric(all_colon[,3])), col=all_colon[,5], horiz=TRUE, 
#         names.arg = all_colon[,1], las=1, xlim = c(-2,6),
#         xlab = "Log2 Fold change")
barplot(as.numeric(all_colon[,3]), col=all_colon[,5], horiz=TRUE, 
        names.arg = all_colon[,1], las=1, xlim = c(-5,15),
        xlab = "Log2 Fold Change")


dev.off()

# df<-t(abundance)
# for(i in 1:ncol(df)){
#   #i<-1
#   chem<-df[,i]
#   #str(chem)
#   #names(chem)<-colnames(df)
#   means<-vector()
#   sds<-vector()
#   k<-1
#   for(tr in unique(meta$Treatment)){
#     #tr<-"GenX_10mg"
#     means[k]<-mean(chem[which(names(chem)%in%meta$`Sample ID`[which(meta$Treatment==tr)])])
#     sds[k]<-sd(chem[which(names(chem)%in%meta$`Sample ID`[which(meta$Treatment==tr)])])
#     k<-k+1
#   }
#   unique(meta$Treatment)
#   
#   cairo_pdf(paste("graphs/Micr_",i,".pdf", sep=""), 
#             width = 3, height = 4) # save plot
#   
#   plot(c(1,3,4,5), means[c(5,1,3,6)],
#        ylim=range(c(means-sds, means+sds)),
#        xlim=c(1,5),
#        xlab="Concentration, mg/L",
#        pch=19, 
#        ylab="Metabolite z-score",
#        main=colnames(df)[i],
#        col=alpha("red", 0.7), xaxt="n"
#   )
#   axis(1, at=c(1,2,3,4,5), labels=c("Ctrl", "5", "10", "20","100") )
#   # hack: we draw arrows but with very special "arrowheads"
#   arrows(c(1,3,4,5), means[c(5,1,3,6)]-sds[c(5,1,3,6)], c(1,3,4,5), means[c(5,1,3,6)]+sds[c(5,1,3,6)], length=0.05, angle=90, code=3, 
#          alpha("red", 0.7), lwd=2)
#   lines(c(1,3,4,5), means[c(5,1,3,6)], col=alpha("red", 0.7), lwd=2)
#   
#   points(c(1,2,3,4), means[c(5,2,4,7)], col=alpha("royal blue", 0.7), pch=16)
#   
#   # hack: we draw arrows but with very special "arrowheads"
#   arrows(c(1,2,3,4), means[c(5,2,4,7)]-sds[c(5,2,4,7)], c(1,2,3,4), means[c(5,2,4,7)]+sds[c(5,2,4,7)], 
#          length=0.05, angle=90, code=3, lwd=2,
#          alpha("royal blue", 0.7))
#   lines(c(1,2,3,4), means[c(5,2,4,7)], col=alpha("royal blue", 0.7), lwd=2)
#   
#   dev.off()
#   
# }
# 





################################################################################################   
###                   Plots for regression analysis      (Fig 2B-C)
################################################################################################
