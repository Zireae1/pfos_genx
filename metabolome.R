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
setwd("/Users/veronika/Desktop/pfos_genx-main/") # replace with your working directory

# load sample metadata
meta<-fread("data/Metadata.csv", stringsAsFactors = F, header = T, data.table = F)
meta$Chemical<-NA
meta$Concentration<-NA
for(i in 1:nrow(meta)){
  meta$Chemical[i]<-strsplit(meta$Treatment[i], "_")[[1]][1]
  meta$Concentration[i]<-strsplit(meta$Treatment[i], "_")[[1]][2]
}
meta<-meta[order(meta$`Sample ID`),]

### set colors      PFOS        GenX
main_palette<-c("#0073C2FF", "#EFC000FF")

# load microbial abundances
liver<-fread("data/Liver_PFOS_GenX_raw_all_processed.txt", 
                  stringsAsFactors = F, header = T, data.table = F)
si<-fread("data/SI_PFOS_GenX_raw_all_processed.txt", 
             stringsAsFactors = F, header = T, data.table = F)
colon<-fread("data/Colon_PFOS_GenX_raw_all_processed.txt", 
          stringsAsFactors = F, header = T, data.table = F)


which(is.na(liver))

################################################################################################   
###                Fig 1E          PCA for metabolome samples      
################################################################################################  

metab<-liver[,-c(1:7)]
colnames(metab)<-gsub("Area: ", "", colnames(metab))
colnames(metab)<-gsub("\\..*", "", colnames(metab))

metab<-si[,-c(1:7)]
colnames(metab)<-gsub("Area: ", "", colnames(metab))
colnames(metab)<-gsub("\\-.*", "", colnames(metab))

meta$`Sample ID`<-gsub("SI_", "", meta$`Sample ID`)

metab<-colon[,-c(1:7)]
colnames(metab)<-gsub("Area: ", "", colnames(metab))
colnames(metab)<-gsub("\\-.*", "", colnames(metab))

meta$`Sample ID`<-gsub("Colon_", "", meta$`Sample ID`)


### log and z-score normalization
metab<-t(apply(metab, 1, function(x) log10(x)))
metab<-t(apply(metab, 1, function(x) (x-mean(x))/sd(x)))


pc<-prcomp(t(metab))
summary(pc)

coord<-as.data.frame(cbind(pc$x[,1],pc$x[,2]))
colnames(coord)<-c("Coord 1", "Coord 2")

perc<-summary(pc)$importance
perc<-round(perc[2, 1:2]*100, digits = 2)

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
#rownames(coord)<-gsub("SI_", "", rownames(coord))

coord$color[which(coord$Chemical=="GenX")]<-main_palette[2]
coord$color[which(coord$Chemical=="PFOS")]<-main_palette[1]
c#oord$pch[which(coord$Part=="Colon")]<-21
#coord$pch[which(coord$Part=="SI")]<-24
coord$pch<-22
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

cairo_pdf(paste("graphs/pcoa_12_Metab_Colon.pdf", sep=""), 
          width = 6.5, height = 6.5) # save plot

plot(coord[,1:2], col="black", bg=coord$color, pch=coord$pch, cex=coord$size, ylim=c(-40,60)
)
text(coord[,1:2], rownames(coord), pos = 3)
legend("topleft", inset=c(0,0), title = c("GenX"),
       bty="n", legend = c("Control","10mg/L", "20mg/L", "100mg/L"), col="black",
       pt.bg = c(alpha(main_palette[2],0.25), alpha(main_palette[2],0.5),alpha(main_palette[2],0.75), alpha(main_palette[2],1)), 
       lty= 0, pch = unique(coord$pch), pt.cex=c(1,1.5,1.75,2))
legend("topleft", inset=c(0.25,0), title = c("PFOS"),
       bty="n", legend = c("Control", "5mg/L", "10mg/L", "20mg/L"), col="black",
       pt.bg = c(alpha(main_palette[1],0.25),alpha(main_palette[1],0.5),
                 alpha(main_palette[1],0.75), alpha(main_palette[1],1)), 
       lty= 0, pch = unique(coord$pch), pt.cex=c(1,1.25,1.5,1.75))

dev.off()

### try euclidean distance:
dm<-as.matrix(dist(t(metab), method = "euclidean"))

### try correlation distance:
dm<-matrix(0, ncol=ncol(metab), nrow=ncol(metab))
for(i in 1:(ncol(dm)-1)){
  for(j in (i+1):ncol(dm)){
    dm[i,j]<-1-cor(metab[,i], metab[,j], method = "spearman")
    dm[j,i]<-dm[i,j]
  }
}
colnames(dm)<-colnames(metab)
rownames(dm)<-colnames(metab)
hist(dm)
max(dm)
dm<-dm/max(dm)

#Run PERMANOVA on distances.
lab<-meta$Treatment[which(meta$`Sample ID`%in%colnames(dm))]
adonis(as.dist(dm) ~ lab,  permutations = 1000)

### Group dissimilarity comparison:
### compare ctrl vs ctrl, ctrl vs PFOS, ctrl vs GenX, PFOS vs GenX (by conc???)
trt<-c("PFOS_C","GenX_C")
rep.ids<-which(colnames(dm)%in%meta$`Sample ID`[which(meta$Treatment%in%trt)])
dist_cvsc<-as.numeric(dm[rep.ids, rep.ids])
dist_cvsc<-dist_cvsc[-which(dist_cvsc==0)]
#dist_cvsa<-as.numeric(dm[-rep.ids, rep.ids])

df<-cbind(unique(dist_cvsc), "ctrl vs ctrl")

trt2<-c("PFOS_5mg","PFOS_10mg","PFOS_20mg")
rep.ids2<-which(colnames(dm)%in%meta$`Sample ID`[which(meta$Treatment%in%trt2)])
dist_cvsa<-as.numeric(dm[rep.ids2, rep.ids])
df<-rbind(df, cbind(dist_cvsa, "ctrl vs PFOS"))

trt3<-c("GenX_10mg","GenX_20mg","GenX_100mg")
rep.ids3<-which(colnames(dm)%in%meta$`Sample ID`[which(meta$Treatment%in%trt3)])
dist_cvsa<-as.numeric(dm[rep.ids3, rep.ids])
df<-rbind(df, cbind(dist_cvsa, "ctrl vs GenX"))

dist_cvsa<-as.numeric(dm[rep.ids3, rep.ids2])
df<-rbind(df, cbind(dist_cvsa, "PFOS vs GenX"))

### add replicates with replicates:
# trts<-c("PFOS_C","PFOS_5mg","PFOS_10mg","PFOS_20mg", "GenX_10mg","GenX_20mg","GenX_100mg")
# for(trt in trts){
#   rep.ids<-which(colnames(dm)%in%meta$`Sample ID`[which(meta$Treatment%in%trt)])
#   dist_cvsc<-as.numeric(dm[rep.ids, rep.ids])
#   dist_cvsc<-dist_cvsc[-which(dist_cvsc==0)]
#   df<-rbind(df, cbind(unique(dist_cvsc), "cRep vs cRep"))
#   
# }



df<-as.data.frame(df)
df$dist_cvsa<-as.numeric(df$dist_cvsa)

cairo_pdf(paste("graphs/Dist_Colon_metab_cord.pdf", sep=""), 
          width = 2.5, height = 5) # save plot

theme_set(
  theme_classic() +
    theme(legend.position = "top")
)

unique(df$V2)
df <- df %>%
  mutate( V2=factor(V2,levels=c("ctrl vs ctrl", "ctrl vs PFOS", "ctrl vs GenX", "PFOS vs GenX")) )


ggplot(df, aes(x = V2, y = dist_cvsa, fill = V2)) + xlab("Compared sample groups") + 
  ylab("Normalized Correlation Distance") +
  ylim(0.2,1.1) +
  geom_violin(alpha = 0.5) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.8, 
               #dotsize=1, 
               binwidth=0.01)+
  #geom_point(position = position_jitter(seed = 1, width = 0.2)) +
  theme(legend.position = "none") +
  scale_fill_manual(values = c(#"violet",
                               "#FC4E07", "#0073C2FF", "#EFC000FF", "olivedrab3")) +
  geom_signif(comparisons = list(c("ctrl vs ctrl", "ctrl vs PFOS"), 
                                 c("ctrl vs ctrl", "ctrl vs GenX"), 
                                 c("ctrl vs GenX", "ctrl vs PFOS")), test = "wilcox.test",
              #y_position=c(0.8, 0.9, 0.85),   
              y_position=c(0.96,1.04, 1.00), #vjust = 0.1,
              map_signif_level=TRUE)

dev.off()


