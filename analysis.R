#source("http://bioconductor.org/biocLite.R")
#biocLite("phyloseq")
#install.packages("plyr")
#install.packages("ggplot2")
#install.packages("data.table")
library(phyloseq)
library(plyr)
library(ggplot2)
library(data.table)
setwd("/Users/adina/scratch/allen")

==============
#LOAD FUNCTIONAL DATA
abund <- read.delim(sep='\t', file="./summary-count-cleaned-norm-bp.tsv",header=TRUE, strip.white=TRUE, row.names=1)
abund <- read.delim(sep='\t', file="./summary-pres-abs.tsv",header=TRUE, strip.white=TRUE, row.names=1)

#Use the following to normalize the unnormalized abundance file
#d2 <- abund/rep(colSums(abund),each=length(abund$X77C_4.2_IND))
#abund <- d2
#write.table(x=d2, file="./abundance_table.norm.txt", sep="\t", col.names=TRUE, row.names=TRUE, quote=FALSE)
meta <- read.delim(sep=',', file="./metadata.txt", header=TRUE, strip.white=TRUE, row.names=1)
ann <- read.delim(sep='\t', file="./primer_annotations_cleaned2.txt",header=FALSE, strip.white=TRUE, row.names=1)
ann_data_matrix <- as.matrix(ann)
annotation <- tax_table(ann_data_matrix)
abundance_data_norm_matrix <- as.matrix(abund)
abundance <- otu_table(abundance_data_norm_matrix, taxa_are_rows=TRUE)
metadata <- sample_data(meta)

#phyloseq object with annotations
all <- phyloseq(metadata, abundance, annotation)
plot_bar(all, "Location",facet_grid=~Medicated, fill="V6")
plot_bar(all, "V6")
plot_heatmap(all, sample.label="Location")
heatmap(otuTable(all))


mdf <- psmelt(all)
#summarizing average over experiment
f <- ddply(mdf, .(Medicated, Location, OTU, V6), summarise, SUM=sum(Abundance))
f2 <- ddply(f, .(Medicated, Location, V6), summarise, MEAN=mean(SUM), SE=sd(SUM)/sqrt(length(SUM)))
f <- ddply(mdf, .(Medicated, Location, OTU), summarise, SUM=sum(Abundance))
f2 <- ddply(f, .(Medicated, Location), summarise, MEAN=mean(SUM), SE=sd(SUM)/sqrt(length(SUM)))
limits<-aes(ymin=MEAN-SE, ymax=MEAN+SE)
p = ggplot(f2, aes_string(x="V6", y="MEAN", color="Medicated"))
p + geom_point()+geom_errorbar(limits, width=0)+facet_grid(Medicated~Location)+theme_bw()+theme(strip.text.x=element_text(angle=90), axis.text.x=element_text(angle=90))
p = ggplot(f2, aes_string(x="Location", y="MEAN", color="Medicated"))
p + geom_point()+geom_errorbar(limits, width=0)+theme_bw()+theme(strip.text.x=element_text(angle=90), axis.text.x=element_text(angle=90))

#NMDS All fractions
all <- phyloseq(metadata, abundance)
#all <- phyloseq(metadata, abundance, annotation)
GPdist = phyloseq::distance(all, "bray")
GPNMDS = ordinate(all, "NMDS", GPdist)
#3/6 3/12 , 3/27, 4/2
#all_t <- subset_samples(all, date == "3/27/13" | date == "4/2/13")
#p2 = plot_ordination(all_t, GPNMDS, color="ext_frac", shape="cage_str")
p2 = plot_ordination(all, GPNMDS, color="ext_frac", shape="cage_str")+theme_bw()+theme(text=element_text(size=20))+geom_point(size=4)+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
p2 + guides(colour=guide_legend(title="DNA Fraction"), shape=guide_legend(title="Baseline diet"))



#Stats on NMDS
x3 <- subset(x2, cage_str == "77")
x2 <- x3
vlp = subset(x2, ext_frac == "ind")
vlp$day_c <- as.character(vlp$day)
summary(aov(dist~day_c, data=vlp))
aov.result<-aov(dist~day_c, data=vlp)
TukeyHSD(aov.result)
ind = subset(x2, ext_frac == "vlp")
ind$day_c <- as.character(ind$day)
summary(aov(dist~day_c, data=ind))
aov.result<-aov(dist~day_c, data=ind)
TukeyHSD(aov.result)



#Looking at clustering
veganotu = function(physeq) {
    require("vegan")
    OTU = otu_table(physeq)
    if (taxa_are_rows(OTU)) {
        OTU = t(OTU)
    }
    return(as(OTU, "matrix"))
}

#Examples on clustering
all <- phyloseq(metadata, abundance, annotation)
all <- phyloseq(metadata, abundance)
all_vlp <- all
all_vlp <- subset_samples(all, diet != 'LF-MF')
all_vlp <- subset_samples(all_vlp, diet != 'MF-LF')
all_vlp <- subset_samples(all, Age == "0" | Age == "43")
#all_vlp <- subset_samples(all_vlp, cage_str == "77")
#all_vlp <- subset_samples(all_vlp, diet2 != 'LF*')
#all_vlp <- subset_samples(all_vlp, diet2 != 'MF*')
meta <- sample_data(all_vlp)
v_all <- veganotu(all_vlp)
bd <- vegdist(v_all, method="bray")
adonis(bd ~ Age, perm=9999, as(sample_data(all_vlp),"data.frame"))
adonis(bd ~ ext_frac, perm=9999, as(sample_data(all_vlp),"data.frame"))
groups <- meta$Age
mod <- betadisper(bd, groups)
anova(mod)
permutest(mod, pairwise = TRUE, permutations=9999)
(mod.HSD <- TukeyHSD(mod))
clust <- hclust(bd, method="average")
setEPS()
postscript("suppfig4.eps")
plot(clust)
dev.off()

#Functions that are statistically significant
all <- phyloseq(metadata, abundance, annotation)
mdf <- psmelt(all)
f <- ddply(mdf, .(Medicated, V6, OTU), summarise, SUM=sum(Abundance))
f2 <- subset(f, V6 == "Bacteroidetes/Chlorobi group")
summary(aov.result <- aov(SUM ~ Medicated, data = f2))

f3 <- subset(f, ext_frac == "ind") #change to vlp also
l = length(levels(f3$l3))
stat_sig = rep(0, l)
for (i in 1:l){
	f_dat <- subset(f3, l3 == levels(f3$l3)[i])
	#print(f_dat)
	summary(aov.result <- aov(SUM ~ diet2, data = f_dat))
	if(summary(aov.result)[[1]][["Pr(>F)"]][1] <= 0.05 & summary(aov.result)[[1]][["Pr(>F)"]][1] != "NaN"){
			print(f_dat)
			print(summary(aov.result <- aov(SUM ~ diet2, data = f_dat)))
			stat_sig[i] = as.character(f_dat$l3[1])
			print(TukeyHSD(aov.result))}
				}
				
#Contigs that are significant				
mdf_phage = subset(mdf, l3 == "Phages, Prophages, Transposable elements, Plasmids")
#mdf_phage = subset(mdf, l3 == "Nitrogen Metabolism")
#mdf_phage = subset(mdf, l3 == "Motility and Chemotaxis")
f <- ddply(mdf_phage, .(OTU, diet2, Age, ext_frac, Sample), summarise, SUM=sum(Abundance))
f3 <- subset(f, ext_frac == "ind")
l = length(unique(f3$OTU))
stat_sig = rep(0, l)
for (i in 1:l){
	f_dat <- subset(f3, OTU == unique(f3$OTU)[i])
	#print(f_dat)
	summary(aov.result <- aov(SUM ~ diet2, data = f_dat))
	if(summary(aov.result)[[1]][["Pr(>F)"]][1] <= 0.05 & summary(aov.result)[[1]][["Pr(>F)"]][1] != "NaN"){
			print(f_dat)
			print(summary(aov.result <- aov(SUM ~ diet2, data = f_dat)))
			pvalue = summary(aov.result)[[1]][["Pr(>F)"]][1]
			stat_sig[i] = as.character(f_dat$OTU[1])
			print(TukeyHSD(aov.result))
			x<-c(unique(f_dat[1]),as.character(f_dat$SUM[1]), as.character(f_dat$SUM[2]),as.character(f_dat$SUM[3]),as.character(f_dat$SUM[4]), as.character(f_dat$SUM[5]),as.character(f_dat$SUM[6]), as.character(f_dat$SUM[7]),as.character(f_dat$SUM[8]), as.character(f_dat$SUM[9]), pvalue)
			write.table(x, file="ind.tsv", append=TRUE, quote=FALSE, sep="\t ", eol="\n", row.names = FALSE, col.names=FALSE)}
	}

#mdf_phage = subset(mdf, l3 == "Phages, Prophages, Transposable elements, Plasmids")
#mdf_phage = subset(mdf, l3 == "Nitrogen Metabolism")
mdf_phage = subset(mdf, l3 == "Motility and Chemotaxis")
f <- ddply(mdf_phage, .(OTU, diet2, Age, ext_frac, Sample), summarise, SUM=sum(Abundance))
f3 <- subset(f, ext_frac == "vlp")
l = length(unique(f3$OTU))
stat_sig = rep(0, l)
for (i in 1:l){
	f_dat <- subset(f3, OTU == unique(f3$OTU)[i])
	#print(f_dat)
	summary(aov.result <- aov(SUM ~ diet2, data = f_dat))
	if(summary(aov.result)[[1]][["Pr(>F)"]][1] <= 0.05 & summary(aov.result)[[1]][["Pr(>F)"]][1] != "NaN"){
			print(f_dat)
			print(summary(aov.result <- aov(SUM ~ diet2, data = f_dat)))
			pvalue = summary(aov.result)[[1]][["Pr(>F)"]][1]
			stat_sig[i] = as.character(f_dat$OTU[1])
			print(TukeyHSD(aov.result))
			x<-c(unique(f_dat[1]),as.character(f_dat$SUM[1]), as.character(f_dat$SUM[2]),as.character(f_dat$SUM[3]),as.character(f_dat$SUM[4]), as.character(f_dat$SUM[5]),as.character(f_dat$SUM[6]), as.character(f_dat$SUM[7]),as.character(f_dat$SUM[8]), as.character(f_dat$SUM[9]), pvalue)
			write.table(x, file="vlp.tsv", append=TRUE, quote=FALSE, sep="\t ", eol="\n", row.names = FALSE, col.names=FALSE)}
	}

	
# Exploring core proteins (Figure 8)
all <- phyloseq(metadata, abundance)
all_virus <- subset_samples(all, ext_frac != "total")
all_virus <- subset_samples(all, ext_frac != "total" & cage_str == "50")
abund <- otu_table(all_virus)
virus_core <- abund[apply(abund, MARGIN=1, function(x) all(x > 0)),]
#virus_core_phy <- phyloseq(virus_core, metadata, annotation)
virus_core_phy <- phyloseq(virus_core, metadata)
mdf <- psmelt(virus_core_phy)
f <- ddply(mdf, .(OTU, cage_str, ext_frac, Age), summarise, MEAN=mean(Abundance), SE=sd(Abundance)/sqrt(length(Abundance)))
limits<-aes(ymin=MEAN-SE, ymax=MEAN+SE)
f2 = subset(f, ext_frac == "ind")
#labels annotated
subsetx <- f2[which(f2$OTU == "10397_10396_IND" | f2$OTU == "10398_10397_IND" | f2$OTU == "11291_11290_IND" | f2$OTU == "11296_11295_IND" | f2$OTU == "11327_11326_IND" | f2$OTU == "11343_11342_IND" | f2$OTU == "1392833_97998_IND" | f2$OTU == "1541070_99692_IND" | f2$OTU == "1541404_99769_IND" | f2$OTU == "1858_1857_VLP" | f2$OTU == "2033_2032_VLP" | f2$OTU == "2250_2249_VLP" | f2$OTU == "22767_6492_TOT" | f2$OTU == "2366_2365_VLP" | f2$OTU == "2368_2367_VLP" | f2$OTU == "26_6575_TOT" | f2$OTU == "2650_1950_TOT" | f2$OTU == "2701_2700_VLP" | f2$OTU == "276638_14923_VLP" | f2$OTU == "28_3490_VLP" | f2$OTU == "3327_3326_VLP" | f2$OTU == "3349_3348_VLP" | f2$OTU == "4585_5883_TOT" | f2$OTU == "6916_6915_IND" | f2$OTU == "74376_10776_VLP" | f2$OTU == "9907_9906_IND"),]
temp <- subset(f2, cage_str == "50" & Age == "0")
temp_foo <- subset(subsetx, cage_str == "50" & Age == "0")
f2$OTU <- factor(f2$OTU, levels=temp$OTU[order(-temp$MEAN)])
subsetx$OTU <- factor(subsetx$OTU, levels=temp_foo$OTU[order(-temp_foo$MEAN)])
p = ggplot(f2, aes_string(x="OTU", y="MEAN"))
p+theme_bw()+geom_point(stackt="identity", aes(color=ext_frac))+geom_errorbar(limits, width=0)+facet_grid(Age~cage_str)

#clustering significant contigs
abund_sig <- abund[which(rownames(abund) ==  "11070_5415_VLP" |  rownames(abund) ==  "11392_11391_IND" |  rownames(abund) ==  "11495_11494_IND" |  rownames(abund) ==  "11500_11499_IND" |  rownames(abund) ==  "11573_11572_IND" |  rownames(abund) ==  "11577_11576_IND" |  rownames(abund) ==  "1168485_95035_IND" |  rownames(abund) ==  "11895_11894_IND" |  rownames(abund) ==  "1301_1300_VLP" |  rownames(abund) ==  "1395_1394_VLP" |  rownames(abund) ==  "1460_1459_VLP" |  rownames(abund) ==  "1541036_99683_IND" |  rownames(abund) ==  "1545310_100359_IND" |  rownames(abund) ==  "1546282_100483_IND" |  rownames(abund) ==  "1568425_102182_IND" |  rownames(abund) ==  "189_188_IND" |  rownames(abund) ==  "2052_2051_IND" |  rownames(abund) ==  "2281_2280_VLP" |  rownames(abund) ==  "268318_14671_VLP" |  rownames(abund) ==  "278888_15308_VLP" |  rownames(abund) ==  "279518_15392_VLP" |  rownames(abund) ==  "282627_15745_VLP" |  rownames(abund) ==  "3314_3313_VLP" |  rownames(abund) ==  "3398_3397_VLP" |  rownames(abund) ==  "345_344_IND" |  rownames(abund) ==  "47_46_IND" |  rownames(abund) ==  "485381_86166_IND" |  rownames(abund) ==  "6233_6232_IND" |  rownames(abund) ==  "81561_11054_VLP" |  rownames(abund) ==  "8182_8181_IND" |  rownames(abund) ==  "8233_8232_IND" |  rownames(abund) ==  "9169_9168_IND" |  rownames(abund) ==  "9619_9618_IND"),]
bcd <- vegdist(abund_sig, method="bray")
#bcd<-vegdist(foo, method="bray")
clust <- hclust(bcd, method="average")
plot(clust)
#saved cluster names to cluster.txt file
cluster <-read.delim(sep='\t', file="./cluster.txt",header=FALSE, strip.white=TRUE, row.names=1)
colnames(cluster) <- c("g")
cluster <- data.frame(cluster)
cluster$g <- as.character(cluster$g)
bcd <- vegdist(abund_sig, method="bray")
adonis(decostand(abund_sig,"total") ~ g,data=cluster,perm=9999)
groups <- cluster$g
mod <- betadisper(bcd, groups)
anova(mod)
permutest(mod, pairwise = TRUE, permutations=9999)
(mod.HSD <- TukeyHSD(mod))

#looking into 33 contigs
head(annotation)
ann_sig <- ann[which(rownames(ann) ==  "11070_5415_VLP" |  rownames(ann) ==  "11392_11391_IND" |  rownames(ann) ==  "11495_11494_IND" |  rownames(ann) ==  "11500_11499_IND" |  rownames(ann) ==  "11573_11572_IND" |  rownames(ann) ==  "11577_11576_IND" |  rownames(ann) ==  "1168485_95035_IND" |  rownames(ann) ==  "11895_11894_IND" |  rownames(ann) ==  "1301_1300_VLP" |  rownames(ann) ==  "1395_1394_VLP" |  rownames(ann) ==  "1460_1459_VLP" |  rownames(ann) ==  "1541036_99683_IND" |  rownames(ann) ==  "1545310_100359_IND" |  rownames(ann) ==  "1546282_100483_IND" |  rownames(ann) ==  "1568425_102182_IND" |  rownames(ann) ==  "189_188_IND" |  rownames(ann) ==  "2052_2051_IND" |  rownames(ann) ==  "2281_2280_VLP" |  rownames(ann) ==  "268318_14671_VLP" |  rownames(ann) ==  "278888_15308_VLP" |  rownames(ann) ==  "279518_15392_VLP" |  rownames(ann) ==  "282627_15745_VLP" |  rownames(ann) ==  "3314_3313_VLP" |  rownames(ann) ==  "3398_3397_VLP" |  rownames(ann) ==  "345_344_IND" |  rownames(ann) ==  "47_46_IND" |  rownames(ann) ==  "485381_86166_IND" |  rownames(ann) ==  "6233_6232_IND" |  rownames(ann) ==  "81561_11054_VLP" |  rownames(ann) ==  "8182_8181_IND" |  rownames(ann) ==  "8233_8232_IND" |  rownames(ann) ==  "9169_9168_IND" |  rownames(ann) ==  "9619_9618_IND"),]
summary(ann_sig)
write.table(ann_sig, file="annotation_of_sig.txt", quote=FALSE,sep="\t")
ann_sig_g1 <- ann[which(rownames(ann) ==  "11070_5415_VLP" |  rownames(ann) ==  "11495_11494_IND" |  rownames(ann) ==  "1168485_95035_IND" |  rownames(ann) ==  "11895_11894_IND" |  rownames(ann) ==  "1301_1300_VLP" |  rownames(ann) ==  "1541036_99683_IND" |  rownames(ann) ==  "1546282_100483_IND" |  rownames(ann) ==  "2052_2051_IND" |  rownames(ann) ==  "3314_3313_VLP" |  rownames(ann) ==  "345_344_IND" |  rownames(ann) ==  "47_46_IND" |  rownames(ann) ==  "485381_86166_IND" |  rownames(ann) ==  "8233_8232_IND" |  rownames(ann) ==  "9619_9618_IND"),]
summary(ann_sig_g1)
ann_sig_g2 <- ann[which(rownames(ann) ==  "11392_11391_IND" |  rownames(ann) ==  "11500_11499_IND" |  rownames(ann) ==  "11573_11572_IND" |  rownames(ann) ==  "11577_11576_IND" |  rownames(ann) ==  "1395_1394_VLP" |  rownames(ann) ==  "1460_1459_VLP" |  rownames(ann) ==  "1545310_100359_IND" |  rownames(ann) ==  "1568425_102182_IND" |  rownames(ann) ==  "189_188_IND" |  rownames(ann) ==  "2281_2280_VLP" |  rownames(ann) ==  "268318_14671_VLP" |  rownames(ann) ==  "278888_15308_VLP" |  rownames(ann) ==  "279518_15392_VLP" |  rownames(ann) ==  "282627_15745_VLP" |  rownames(ann) ==  "3398_3397_VLP" |  rownames(ann) ==  "6233_6232_IND" |  rownames(ann) ==  "81561_11054_VLP" |  rownames(ann) ==  "8182_8181_IND" |  rownames(ann) ==  "9169_9168_IND"),]
summary(ann_sig_g2)

ann_data_matrix_g1 <- as.matrix(ann_sig_g1)
annotation_g1 <- tax_table(ann_data_matrix_g1)
g1 <- phyloseq(annotation_g1, abundance, metadata)
mdf <- psmelt(g1)
f <- ddply(mdf, .(l3, diet2, Age, ext_frac, cage_str, Sample), summarise, SUM=sum(Abundance))
f2 <- ddply(f, .(l3, diet2, Age, ext_frac, cage_str), summarise, MEAN=mean(SUM), SE=sd(SUM)/sqrt(length(SUM)))
limits<-aes(ymin=MEAN-SE, ymax=MEAN+SE)
ggsave("group1.eps")
p = ggplot(f2, aes_string(x="Age", y="MEAN", shape="ext_frac", color="cage_str"))+facet_grid(ext_frac~cage_str)
p+theme_bw() + theme(strip.background = element_rect(colour='white',fill = 'white'))+geom_point(stat="identity", sdize=7)+geom_errorbar(limits, width=0)+ theme(panel.grid.major=element_blank(), panel.grid.minor = element_blank())+ylab("Relative Abundance")+theme(text=element_text(size=25, family="Helvetica"), legend.position="none")+theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5,size=15))+scale_colour_manual(values=c("blue","orange"))
dev.off()





ann_data_matrix_g2 <- as.matrix(ann_sig_g2)
annotation_g2 <- tax_table(ann_data_matrix_g2)
g2 <- phyloseq(annotation_g2, abundance, metadata)
mdf <- psmelt(g2)
f <- ddply(mdf, .(l3, diet2, Age, ext_frac, cage_str, Sample), summarise, SUM=sum(Abundance))
f2 <- ddply(f, .(l3, diet2, Age, ext_frac, cage_str), summarise, MEAN=mean(SUM), SE=sd(SUM)/sqrt(length(SUM)))
limits<-aes(ymin=MEAN-SE, ymax=MEAN+SE)
ggsave("group2.eps")
p = ggplot(f2, aes_string(x="Age", y="MEAN", shape="ext_frac", color="cage_str"))+facet_grid(ext_frac~cage_str)
p+theme_bw() + theme(strip.background = element_rect(colour='white',fill = 'white'))+geom_point(stat="identity", sdize=7)+geom_errorbar(limits, width=0)+ theme(panel.grid.major=element_blank(), panel.grid.minor = element_blank())+ylab("Relative Abundance")+theme(text=element_text(size=25, family="Helvetica"), legend.position="none")+theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5,size=15))+scale_colour_manual(values=c("blue","orange"))
dev.off()

abundance_data_norm_matrix <- as.matrix(abund_sig)
abundance <- otu_table(abundance_data_norm_matrix, taxa_are_rows=TRUE)
metadata <- sample_data(meta)

#phyloseq object with annotations
all <- phyloseq(metadata, abundance)
mdf<-psmelt(all)
mdf <- subset(mdf, cage_str == "77" & ext_frac == "vlp" & Age == "0")
f <- ddply(mdf, .(OTU, Age, ext_frac, Sample), summarise, SUM=sum(Abundance))
f2 <- ddply(f, .(OTU, Age, ext_frac), summarise, MEAN=mean(SUM), SE=sd(SUM)/sqrt(length(SUM)))

