# Haemonchus genome - genome variation analyses

## Table of contents

1.
2. [World map of sampling sites](#global_sampling_map)
3. [Genome wide nucleotide diversity analysis](#nuc_div)



# PCA of mtDNA genotypes
#--- filter
```shell
bcftools-1.9 view -e 'FORMAT/DP[0]<10 | MQ[*]<30' 7.hcontortus_chr_mtDNA_arrow_pilon.cohort.vcf.gz | bcftools-1.9 view -i 'TYPE="snp" & AF>0.01' -O z -o allsamples.mtDNA.filtered.vcf.gz

awk -F '[_]' '{print $0,$1,$2}' OFS="\t" samples.list > samples.pops.list
```
```R
R-3.5.0
library(gdsfmt)
library(SNPRelate)
library(ggplot2)

vcf.in <- "allsamples.mtDNA.filtered.vcf.gz"
gds<-snpgdsVCF2GDS(vcf.in, "mtDNA.gds", method="biallelic.only")

genofile <- snpgdsOpen(gds)

pca	<-	snpgdsPCA(genofile, num.thread=2,autosome.only = F)

pops<-	read.table("samples.pops.list",header=F)

tab <- data.frame(sample.id = pca$sample.id,
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],
                  EV3 = pca$eigenvect[,3],
                  EV4 = pca$eigenvect[,4],
                  EV5 = pca$eigenvect[,5],
                  EV6 = pca$eigenvect[,6],    # the second eigenvector
                  COUNTRY = pops$V2,
                  POP = pops$V3,
                  stringsAsFactors = FALSE)

plot(tab$EV2, tab$EV1, xlab="eigenvector 2", ylab="eigenvector 1",pch=20,cex=2,col=pops$V2)
```



```R
R-3.5.0
library(vcfR)
library(poppr)
library(ape)
library(RColorBrewer)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(patchwork)

metadata<-read.table("sample_metadata_colours.list",header=T,comment.char="")

rubi.VCF <- read.vcfR("allsamples.mtDNA.filtered.vcf.gz")
pop.data <- read.table("samples.pops.list", sep = "\t", header = F)
gl.rubi <- vcfR2genlight(rubi.VCF)
ploidy(gl.rubi) <- 1

pop(gl.rubi) <- metadata$country



# distance matrix from genlight object
x.dist <- poppr::bitwise.dist(gl.rubi)



# make a tree
tree <- aboot(gl.rubi, tree = "upgma", distance = bitwise.dist, sample = 100, showtree = F, cutoff = 50, quiet = T)
write.tree(tree, file="MyNewickTreefile.nwk")


cols <- brewer.pal(n = nPop(gl.rubi), name = "Dark2")
plot.phylo(tree, cex = 0.3, font = 2, adj = 0)
nodelabels(tree$node.label, adj = c(1.3, -0.5), frame = "n", cex = 0.3,font = 3, xpd = TRUE)
#legend(35,10,c("CA","OR","WA"),cols, border = FALSE, bty = "n")
legend('topleft', legend = c("CA","OR","WA"),fill = cols, border = FALSE, bty = "n", cex = 2)
axis(side = 1)
title(xlab = "Genetic distance (proportion of loci that are different)")





# pca


rubi.pca <- glPca(gl.rubi, nf = 10)
var_frac <- rubi.pca$eig/sum(rubi.pca$eig)*100
rubi.pca.scores <- as.data.frame(rubi.pca$scores)
rubi.pca.scores$pop <- pop(gl.rubi)
rubi.pca.scores$strain <- metadata$strain
set.seed(9)


#--- plot eigenvectors
barplot(100*rubi.pca$eig/sum(rubi.pca$eig), col = heat.colors(50), main="PCA Eigenvalues")
title(ylab="Percent of variance\nexplained", line = 2)
title(xlab="Eigenvalues", line = 1)


#--- plot PCA


p12 <- ggplot(rubi.pca.scores, aes(x=PC1, y=PC2, colour=pop, label=pop)) + geom_point(size=2)+ theme_bw() + geom_text_repel(data = subset(rubi.pca.scores, pop == "ZAI" ))
p34 <- ggplot(rubi.pca.scores, aes(x=PC3, y=PC4, colour=pop, label=pop)) + geom_point(size=2)+ theme_bw() + geom_text_repel(data = subset(rubi.pca.scores, pop == "ZAI" ))
p56 <- ggplot(rubi.pca.scores, aes(x=PC5, y=PC6, colour=pop, label=pop)) + geom_point(size=2)+ theme_bw() + geom_text_repel(data = subset(rubi.pca.scores, pop == "ZAI" ))
p12 + p34 + p56


ggplot(rubi.pca.scores, aes(x=PC1, y=PC2, colour=pop, label=pop)) + geom_point(size=3)+ theme_bw()


# netview - this all works, but is not very informative
# https://github.com/esteinig/netview/blob/master/tutorials/PearlOysterTutorial.md
#
# library("netview")
# library(DT)
# library(networkD3)
#
# rubi.dist <- bitwise.dist(gl.rubi)
# rubi.dist.matrix <- as.matrix(rubi.dist)
#
# metadata <- read.table("sample_metadata_colours.list",header=T,comment.char="")
#
# oysterOptions <- netviewOptions(selectionTitle="k-Selection", nodeID="sample_id", nodeGroup="country", nodeColour="country_colour", communityAlgorithms=c("Walktrap", "Infomap", "Fast-Greedy"))
#
# graphs <- netview(rubi.dist.matrix, metadata, k=1:60, cluster = TRUE, options=oysterOptions)
# kPlot <- plotSelection(graphs, options=oysterOptions)
#
# k20 <- graphs$k20
# plot(k20, vertex.size=7, vertex.label=NA)
# legend('topleft',legend=levels(as.factor(metadata$country)),col=levels(as.factor(metadata$country_colour)),pch=20)




# DAPC

pnw.dapc <- dapc(gl.rubi, n.pca = 3, n.da = 2)

scatter(pnw.dapc, cex = 2, legend = TRUE, clabel = F, posi.leg = "bottomleft", scree.pca = TRUE, posi.pca = "topleft", cleg = 0.75)



dapc.results <- as.data.frame(pnw.dapc$posterior)
dapc.results$pop <- pop(gl.rubi)
dapc.results$indNames <- rownames(dapc.results)
library(reshape2)
dapc.results <- melt(dapc.results)
colnames(dapc.results) <- c("Original_Pop","Sample","Assigned_Pop","Posterior_membership_probability")
p <- ggplot(dapc.results, aes(x=Sample, y=Posterior_membership_probability, fill=Assigned_Pop))
p <- p + geom_bar(stat='identity')
p <- p + facet_grid(~Original_Pop, scales = "free")
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))
p


### to do - look at variant contribution from dapc to ID SNPs with greatest discrimination power - could rank SNPs this way
#eg.

contrib <- loadingplot(pnw.dapc$var.contr, axis = 2, thres = 0.07, lab.jitter = 1)
# setting axis between 1 or 2 will ID SNPs with greatest impact -  dapc was made with only two axes originally, however this could be increased.
heatmap.2(pnw.dapc$var.contr)



# PCA plot comparison between SNPrealte and DAPC
dapc_PCA <- ggplot(rubi.pca.scores, aes(x=PC1, y=PC2, colour=pop, label=pop)) + geom_point(size=3)+ theme_bw()+stat_ellipse(level = 0.95, size = 1)
snprelate_PCA <- ggplot()+geom_point(aes(tab$EV1*-1, tab$EV2,group=metadata$country,col=metadata$country),size=3)+scale_fill_manual(values=metadata$country_colour)+ theme_bw()
snprelate_PCA + dapc_PCA


# check subpopulaitons within each country - simply change the pop code in the geom_text_repel section
dapc_PCA <- ggplot(rubi.pca.scores, aes(x=PC1, y=PC2, colour=pop, label=strain))+
				geom_point(size=3)+ theme_bw()+
				geom_text_repel(data = subset(rubi.pca.scores, pop == "GB"  ))+
				scale_fill_manual(values=metadata$country_colour)


# check subpopulaitons within each country - simply change the pop code in the geom_text_repel section

ch_data	<-	rubi.pca.scores[(rubi.pca.scores$pop=="CH"),]
gb_data	<-	rubi.pca.scores[(rubi.pca.scores$pop=="GB"),]
pk_data	<-	rubi.pca.scores[(rubi.pca.scores$pop=="PK"),]
us_data	<-	rubi.pca.scores[(rubi.pca.scores$pop=="US"),]
new_data <- dplyr::bind_rows(ch_data,gb_data,pk_data,us_data)

final_PCA <- ggplot()+
			geom_point(aes(rubi.pca.scores$PC1, rubi.pca.scores$PC2, colour=rubi.pca.scores$pop),alpha=1,size=2,stroke = NA)+
			geom_point(aes(new_data$PC1, new_data$PC2, colour=new_data$pop),size=2,stroke = NA)+
			theme_bw()+
			scale_fill_manual(values=metadata$country_colour)+
			xlab(paste("PC1: variance = ",var_frac[1]))+ylab(paste("PC2: variance = ",var_frac[2]))

```			












### 02 - World map of sampling sites <a name="global_sampling_map"></a>
Make a map of H.contortus sampling sites from global population set

Working environment
```shell
cd /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/POPULATION_DIVERSITY
```



```R
R-3.5.0

library(ggplot2)
library(ggmap)
library(maps)
library(mapdata)
library(dplyr)
library(ggrepel)

metadata = read.delim("global_sampling_coords.txt",header=TRUE,sep="\t")

palette(c("#31A197","#E15956","#EF724B","#D35E5C","#606EB8","#6570B0","#34AFE7","#6973A8","#3C9C93","#6E75A0","#C56462","#B76968","#A64EB4","#727898","#A96F6E","#9B7474","#3FA8D8","#8D7A7A"))

pdf("global_sampling_map1.pdf",useDingbats=FALSE)
par(fg = "black")
map("world",col="grey85",fill=TRUE, border=FALSE)
map.axes()
points(metadata$lon, metadata$lat, cex=1, pch=c(16,17)[as.numeric(metadata$dataset)],col=metadata$country_code)
legend( x="bottomright", legend=c("New data; n = 74","Salle et al (2018); n = 264"),col=c("black"), lwd="1", lty=c(0,0), pch=c(17,16),box.lwd = 0,cex = 0.9)

dev.off()

pdf("global_sampling_map_inset.pdf",useDingbats=FALSE)
par(fg = "white")
map("world", col="grey85",fill=TRUE, border=TRUE, xlim=c(-25,25), ylim=c(35,65))
#map.axes()
points(metadata$lon, metadata$lat, cex=1.5, pch=c(16,17)[as.numeric(metadata$dataset)],col=metadata$country_code)
dev.off()


```

```shell
scp sd21@pcs5.internal.sanger.ac.uk:/nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/POPULATION_DIVERSITY/global_sampling_map*  ~/Documents/workbook/hcontortus_genome/04_analysis

global_sampling_map*

```





### 03 - Genome wide nucleotide diversity analysis <a name="nuc_div"></a>

Working environment
```shell
# working dir
cd /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/POPULATION_DIVERSITY/VARIANTS/VCFTOOLS

# get the refrence files
ln -s /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/TRANSCRIPTOME/TRANSCRIPTOME_CURATION/HCON_V4_WBP11plus_190118.ips.gff3 ANNOTATION.gff

ln -s /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/REF/HAEM_V4_final.chr.fa REF.fa
```


```shell
# calculate fst for all populaitons in 100 kbp windows
vcftools-0.1.14 --gzvcf ../1.hcontortus_chr1_Celeg_TT_arrow_pilon.cohort.vcf.gz --fst-window-size 100000 --weir-fst-pop ZAI.population.list --weir-fst-pop AUS.population.list --weir-fst-pop GB.population.list --weir-fst-pop CAP.population.list --weir-fst-pop FRG.population.list --weir-fst-pop ITA.population.list --weir-fst-pop PK.population.list --weir-fst-pop US.population.list --weir-fst-pop BEN.population.list --weir-fst-pop CH.population.list --weir-fst-pop MOR.population.list --weir-fst-pop POR.population.list --weir-fst-pop STO.population.list --weir-fst-pop STA.population.list --weir-fst-pop BRA.population.list --weir-fst-pop FRA.population.list --weir-fst-pop IND.population.list --weir-fst-pop NAM.population.list --out chr1_fst_100k_allpop

vcftools-0.1.14 --gzvcf ../2.hcontortus_chr2_Celeg_TT_arrow_pilon.cohort.vcf.gz --fst-window-size 100000 --weir-fst-pop ZAI.population.list --weir-fst-pop AUS.population.list --weir-fst-pop GB.population.list --weir-fst-pop CAP.population.list --weir-fst-pop FRG.population.list --weir-fst-pop ITA.population.list --weir-fst-pop PK.population.list --weir-fst-pop US.population.list --weir-fst-pop BEN.population.list --weir-fst-pop CH.population.list --weir-fst-pop MOR.population.list --weir-fst-pop POR.population.list --weir-fst-pop STO.population.list --weir-fst-pop STA.population.list --weir-fst-pop BRA.population.list --weir-fst-pop FRA.population.list --weir-fst-pop IND.population.list --weir-fst-pop NAM.population.list --out chr2_fst_100k_allpop

vcftools-0.1.14 --gzvcf ../3.hcontortus_chr3_Celeg_TT_arrow_pilon.cohort.vcf.gz --fst-window-size 100000 --weir-fst-pop ZAI.population.list --weir-fst-pop AUS.population.list --weir-fst-pop GB.population.list --weir-fst-pop CAP.population.list --weir-fst-pop FRG.population.list --weir-fst-pop ITA.population.list --weir-fst-pop PK.population.list --weir-fst-pop US.population.list --weir-fst-pop BEN.population.list --weir-fst-pop CH.population.list --weir-fst-pop MOR.population.list --weir-fst-pop POR.population.list --weir-fst-pop STO.population.list --weir-fst-pop STA.population.list --weir-fst-pop BRA.population.list --weir-fst-pop FRA.population.list --weir-fst-pop IND.population.list --weir-fst-pop NAM.population.list --out chr3_fst_100k_allpop

vcftools-0.1.14 --gzvcf ../4.hcontortus_chr4_Celeg_TT_arrow_pilon.cohort.vcf.gz --fst-window-size 100000 --weir-fst-pop ZAI.population.list --weir-fst-pop AUS.population.list --weir-fst-pop GB.population.list --weir-fst-pop CAP.population.list --weir-fst-pop FRG.population.list --weir-fst-pop ITA.population.list --weir-fst-pop PK.population.list --weir-fst-pop US.population.list --weir-fst-pop BEN.population.list --weir-fst-pop CH.population.list --weir-fst-pop MOR.population.list --weir-fst-pop POR.population.list --weir-fst-pop STO.population.list --weir-fst-pop STA.population.list --weir-fst-pop BRA.population.list --weir-fst-pop FRA.population.list --weir-fst-pop IND.population.list --weir-fst-pop NAM.population.list --out chr4_fst_100k_allpop

vcftools-0.1.14 --gzvcf ../5.hcontortus_chr5_Celeg_TT_arrow_pilon.cohort.vcf.gz --fst-window-size 100000 --weir-fst-pop ZAI.population.list --weir-fst-pop AUS.population.list --weir-fst-pop GB.population.list --weir-fst-pop CAP.population.list --weir-fst-pop FRG.population.list --weir-fst-pop ITA.population.list --weir-fst-pop PK.population.list --weir-fst-pop US.population.list --weir-fst-pop BEN.population.list --weir-fst-pop CH.population.list --weir-fst-pop MOR.population.list --weir-fst-pop POR.population.list --weir-fst-pop STO.population.list --weir-fst-pop STA.population.list --weir-fst-pop BRA.population.list --weir-fst-pop FRA.population.list --weir-fst-pop IND.population.list --weir-fst-pop NAM.population.list --out chr5_fst_100k_allpop

vcftools-0.1.14 --gzvcf ../6.hcontortus_chrX_Celeg_TT_arrow_pilon.cohort.vcf.gz --fst-window-size 100000 --weir-fst-pop ZAI.population.list --weir-fst-pop AUS.population.list --weir-fst-pop GB.population.list --weir-fst-pop CAP.population.list --weir-fst-pop FRG.population.list --weir-fst-pop ITA.population.list --weir-fst-pop PK.population.list --weir-fst-pop US.population.list --weir-fst-pop BEN.population.list --weir-fst-pop CH.population.list --weir-fst-pop MOR.population.list --weir-fst-pop POR.population.list --weir-fst-pop STO.population.list --weir-fst-pop STA.population.list --weir-fst-pop BRA.population.list --weir-fst-pop FRA.population.list --weir-fst-pop IND.population.list --weir-fst-pop NAM.population.list --out chrX_fst_100k_allpop

vcftools-0.1.14 --gzvcf ../7.hcontortus_chr_mtDNA_arrow_pilon.cohort.vcf.gz --fst-window-size 1000 --weir-fst-pop ZAI.population.list --weir-fst-pop AUS.population.list --weir-fst-pop GB.population.list --weir-fst-pop CAP.population.list --weir-fst-pop FRG.population.list --weir-fst-pop ITA.population.list --weir-fst-pop PK.population.list --weir-fst-pop US.population.list --weir-fst-pop BEN.population.list --weir-fst-pop CH.population.list --weir-fst-pop MOR.population.list --weir-fst-pop POR.population.list --weir-fst-pop STO.population.list --weir-fst-pop STA.population.list --weir-fst-pop BRA.population.list --weir-fst-pop FRA.population.list --weir-fst-pop IND.population.list --weir-fst-pop NAM.population.list --out chrMT_fst_1k_allpop


# determine the position of singleton SNPs
vcftools-0.1.14 --gzvcf ../1.hcontortus_chr1_Celeg_TT_arrow_pilon.cohort.vcf.gz --singletons --out singletons_chr1
vcftools-0.1.14 --gzvcf ../2.hcontortus_chr2_Celeg_TT_arrow_pilon.cohort.vcf.gz --singletons --out singletons_chr2
vcftools-0.1.14 --gzvcf ../3.hcontortus_chr3_Celeg_TT_arrow_pilon.cohort.vcf.gz --singletons --out singletons_chr3
vcftools-0.1.14 --gzvcf ../4.hcontortus_chr4_Celeg_TT_arrow_pilon.cohort.vcf.gz --singletons --out singletons_chr4
vcftools-0.1.14 --gzvcf ../5.hcontortus_chr5_Celeg_TT_arrow_pilon.cohort.vcf.gz --singletons --out singletons_chr5
vcftools-0.1.14 --gzvcf ../6.hcontortus_chrX_Celeg_TT_arrow_pilon.cohort.vcf.gz --singletons --out singletons_chrX

# make a position file out of the singletons
for i in singletons_*; do cut -f1,2 ${i} > ${i}.2; done


# calculate window Pi 100kb of singletons only
echo -e "
vcftools-0.1.14 --gzvcf ../1.hcontortus_chr1_Celeg_TT_arrow_pilon.cohort.vcf.gz --positions singletons_chr1.singletons.2 --window-pi 100000 --out chr1.singletons.pi

vcftools-0.1.14 --gzvcf ../2.hcontortus_chr2_Celeg_TT_arrow_pilon.cohort.vcf.gz --positions singletons_chr2.singletons.2 --window-pi 100000 --out chr2.singletons.pi

vcftools-0.1.14 --gzvcf ../3.hcontortus_chr3_Celeg_TT_arrow_pilon.cohort.vcf.gz --positions singletons_chr3.singletons.2 --window-pi 100000 --out chr3.singletons.pi

vcftools-0.1.14 --gzvcf ../4.hcontortus_chr4_Celeg_TT_arrow_pilon.cohort.vcf.gz --positions singletons_chr4.singletons.2 --window-pi 100000 --out chr4.singletons.pi

vcftools-0.1.14 --gzvcf ../5.hcontortus_chr5_Celeg_TT_arrow_pilon.cohort.vcf.gz --positions singletons_chr5.singletons.2 --window-pi 100000 --out chr5.singletons.pi

vcftools-0.1.14 --gzvcf ../6.hcontortus_chrX_Celeg_TT_arrow_pilon.cohort.vcf.gz --positions singletons_chrX.singletons.2 --window-pi 100000 --out chrX.singletons.pi
" > run_nucdiv_singletons
chmod a+x run_nucdiv_singletons
bsub.py --queue yesterday 10 get_singletons ./run_nucdiv_singletons


# calculate window Pi 100kb of non-singletons only
echo -e "
vcftools-0.1.14 --gzvcf ../1.hcontortus_chr1_Celeg_TT_arrow_pilon.cohort.vcf.gz --exclude-positions singletons_chr1.singletons.2 --window-pi 100000 --out chr1.shared.pi

vcftools-0.1.14 --gzvcf ../2.hcontortus_chr2_Celeg_TT_arrow_pilon.cohort.vcf.gz --exclude-positions singletons_chr2.singletons.2 --window-pi 100000 --out chr2.shared.pi

vcftools-0.1.14 --gzvcf ../3.hcontortus_chr3_Celeg_TT_arrow_pilon.cohort.vcf.gz --exclude-positions singletons_chr3.singletons.2 --window-pi 100000 --out chr3.shared.pi

vcftools-0.1.14 --gzvcf ../4.hcontortus_chr4_Celeg_TT_arrow_pilon.cohort.vcf.gz --exclude-positions singletons_chr4.singletons.2 --window-pi 100000 --out chr4.shared.pi

vcftools-0.1.14 --gzvcf ../5.hcontortus_chr5_Celeg_TT_arrow_pilon.cohort.vcf.gz --exclude-positions singletons_chr5.singletons.2 --window-pi 100000 --out chr5.shared.pi

vcftools-0.1.14 --gzvcf ../6.hcontortus_chrX_Celeg_TT_arrow_pilon.cohort.vcf.gz --exclude-positions singletons_chrX.singletons.2 --window-pi 100000 --out chrX.shared.pi
" > run_nucdiv_shared
chmod a+x run_nucdiv_shared
bsub.py --queue yesterday 10 get_shared ./run_nucdiv_shared


# make a gene coords file for chromosome
cut -f1 ANNOTATION.gff  | grep -v "#" |  sort | uniq | while read -r CHR; do awk -v CHR=$CHR '{if($1==CHR && $3=="gene") print $1,($5+$4)/2,$7}' OFS="\t" ANNOTATION.gff  > ${CHR}.genepos.list;  done

python3 ~alt/python/bin/fasta_gaps_to_bed.py REF.fa > HAEM_V4_final.chr.Ns.bed
```


```R
R-3.5.0
library(ggplot2)
library(patchwork)

chr1<-read.table("hcontortus_chr1_Celeg_TT_arrow_pilon.genepos.list",header=F)
chr2<-read.table("hcontortus_chr2_Celeg_TT_arrow_pilon.genepos.list",header=F)
chr3<-read.table("hcontortus_chr3_Celeg_TT_arrow_pilon.genepos.list",header=F)
chr4<-read.table("hcontortus_chr4_Celeg_TT_arrow_pilon.genepos.list",header=F)
chr5<-read.table("hcontortus_chr5_Celeg_TT_arrow_pilon.genepos.list",header=F)
chrX<-read.table("hcontortus_chrX_Celeg_TT_arrow_pilon.genepos.list",header=F)
NNN <- read.table("HAEM_V4_final.chr.Ns.bed",header=F)

chr1pos<-chr1[chr1$V3=="+",]
chr1neg<-chr1[chr1$V3=="-",]
chr2pos<-chr2[chr2$V3=="+",]
chr2neg<-chr2[chr2$V3=="-",]
chr3pos<-chr3[chr3$V3=="+",]
chr3neg<-chr3[chr3$V3=="-",]
chr4pos<-chr4[chr4$V3=="+",]
chr4neg<-chr4[chr4$V3=="-",]
chr5pos<-chr5[chr5$V3=="+",]
chr5neg<-chr5[chr5$V3=="-",]
chrXpos<-chrX[chrX$V3=="+",]
chrXneg<-chrX[chrX$V3=="-",]
chr1N<-NNN[NNN$V1=="hcontortus_chr1_Celeg_TT_arrow_pilon",]
chr2N<-NNN[NNN$V1=="hcontortus_chr2_Celeg_TT_arrow_pilon",]
chr3N<-NNN[NNN$V1=="hcontortus_chr3_Celeg_TT_arrow_pilon",]
chr4N<-NNN[NNN$V1=="hcontortus_chr4_Celeg_TT_arrow_pilon",]
chr5N<-NNN[NNN$V1=="hcontortus_chr5_Celeg_TT_arrow_pilon",]
chrXN<-NNN[NNN$V1=="hcontortus_chrX_Celeg_TT_arrow_pilon",]

chr1_single_pi<-read.table("chr1.singletons.pi.windowed.pi",header=T)			
chr1_shared_pi<-read.table("chr1.shared.pi.windowed.pi",header=T)			
chr1_fst <- read.table("chr1_fst_100k_allpop.windowed.weir.fst",header=T)

chr2_single_pi<-read.table("chr2.singletons.pi.windowed.pi",header=T)			
chr2_shared_pi<-read.table("chr2.shared.pi.windowed.pi",header=T)			
chr2_fst <- read.table("chr2_fst_100k_allpop.windowed.weir.fst",header=T)


chr3_single_pi<-read.table("chr3.singletons.pi.windowed.pi",header=T)			
chr3_shared_pi<-read.table("chr3.shared.pi.windowed.pi",header=T)			
chr3_fst <- read.table("chr3_fst_100k_allpop.windowed.weir.fst",header=T)


chr4_single_pi<-read.table("chr4.singletons.pi.windowed.pi",header=T)			
chr4_shared_pi<-read.table("chr4.shared.pi.windowed.pi",header=T)			
chr4_fst <- read.table("chr4_fst_100k_allpop.windowed.weir.fst",header=T)


chr5_single_pi<-read.table("chr5.singletons.pi.windowed.pi",header=T)			
chr5_shared_pi<-read.table("chr5.shared.pi.windowed.pi",header=T)			
chr5_fst <- read.table("chr5_fst_100k_allpop.windowed.weir.fst",header=T)


chrX_single_pi<-read.table("chrX.singletons.pi.windowed.pi",header=T)			
chrX_shared_pi<-read.table("chrX.shared.pi.windowed.pi",header=T)			
chrX_fst <- read.table("chrX_fst_100k_allpop.windowed.weir.fst",header=T)

# calculate genome wide average Fst, and confidence intervals
fst_all <- rbind(chr1_fst,chr2_fst,chr3_fst,chr4_fst,chr5_fst,chrX_fst)
error <- qt(0.975,df=length(fst_all$WEIGHTED_FST)-1)*sd(fst_all$WEIGHTED_FST)/sqrt(length(fst_all$WEIGHTED_FST))
fst_upper_ci <- mean(fst_all$WEIGHTED_FST)+error
fst_lower_ci <- mean(fst_all$WEIGHTED_FST)-error
quantile<-quantile(chr1_fst$WEIGHTED_FST, c(.05,.95))


chr1_div_plot <- ggplot()+
               geom_point(aes(chr1_single_pi$BIN_START,chr1_single_pi$PI),alpha=0.3,colour="blue",size=0.5)+
               geom_smooth(aes(chr1_single_pi$BIN_START,chr1_single_pi$PI),span=0.1,colour="blue")+
               geom_point(aes(chr1_shared_pi$BIN_START,chr1_shared_pi$PI),alpha=0.3,colour="red",size=0.5)+
               geom_smooth(aes(chr1_shared_pi$BIN_START,chr1_shared_pi$PI),span = 0.1,colour="red")+
               xlim(0,5.2e7)+ylim(0,0.008)+
               xlab("")+ylab(expression(pi))+
               theme_classic()+theme(axis.title.x=element_blank(),
               axis.text.x=element_blank(),
               axis.ticks.x=element_blank())

chr1_fst_plot <- ggplot()+
               geom_rect(aes(xmin=0,ymin=quantile[1],xmax=5.2e7,ymax=quantile[2]),fill="grey90")+
               geom_point(aes(chr1_fst$BIN_START,chr1_fst$WEIGHTED_FST),colour="red",alpha=0.3,size=0.5)+
               theme_classic()+
               xlim(0,5.2e7)+ylim(0,0.4)+
               xlab("")+ylab("Fst")+
               theme(axis.title.x=element_blank(),
               axis.text.x=element_blank(),
               axis.ticks.x=element_blank())


chr1_gene_plot	<-	ggplot()+
               geom_rect(aes(xmin=1,ymin=-1,xmax=max(chr1$V2),ymax=1),fill="#b2182b",alpha=0.5)+
               geom_rect(aes(xmin=chr1N$V2,ymin=-1,xmax=chr1N$V3,ymax=1),fill="white",linetype=0)+
               geom_linerange(aes(chr1pos$V2,ymin=0,ymax=1),size=0.1,alpha=0.2)+
               geom_linerange(aes(chr1neg$V2,ymin=-1,ymax=0),size=0.1,alpha=0.2)+
               theme_classic()+ scale_y_continuous(breaks=NULL)+
               xlim(0,5.2e7)+ylab("I")+xlab("")

chr1_blank <- ggplot()+geom_blank()+xlim(0,5.2e7)+theme_classic()+
               theme(axis.title.x=element_blank(),
               axis.text.x=element_blank(),
               axis.ticks.x=element_blank(),
               axis.line = element_blank())



chr2_div_plot <- ggplot()+
               geom_point(aes(chr2_single_pi$BIN_START,chr2_single_pi$PI),alpha=0.3,colour="blue",size=0.5)+
               geom_smooth(aes(chr2_single_pi$BIN_START,chr2_single_pi$PI),span=0.1,colour="blue")+
               geom_point(aes(chr2_shared_pi$BIN_START,chr2_shared_pi$PI),alpha=0.3,colour="red",size=0.5)+
               geom_smooth(aes(chr2_shared_pi$BIN_START,chr2_shared_pi$PI),span = 0.1,colour="red")+
               xlim(0,5.2e7)+ylim(0,0.008)+
               xlab("")+ylab(expression(pi))+
               theme_classic()+theme(axis.title.x=element_blank(),
               axis.text.x=element_blank(),
               axis.ticks.x=element_blank())

chr2_fst_plot <- ggplot()+
               geom_rect(aes(xmin=0,ymin=quantile[1],xmax=5.2e7,ymax=quantile[2]),fill="grey90")+
               geom_point(aes(chr2_fst$BIN_START,chr2_fst$WEIGHTED_FST),colour="red",alpha=0.3,size=0.5)+
               theme_classic()+
               xlim(0,5.2e7)+ylim(0,0.4)+
               xlab("")+ylab("Fst")+
               theme(axis.title.x=element_blank(),
               axis.text.x=element_blank(),
               axis.ticks.x=element_blank())

chr2_gene_plot <- ggplot()+
               geom_rect(aes(xmin=1,ymin=-1,xmax=max(chr2$V2),ymax=1),fill="#fc8d59",alpha=0.5)+
               geom_rect(aes(xmin=chr2N$V2,ymin=-1,xmax=chr2N$V3,ymax=1),fill="white",linetype=0)+
               geom_linerange(aes(chr2pos$V2,ymin=0,ymax=1),size=0.1,alpha=0.3)+
               geom_linerange(aes(chr2neg$V2,ymin=-1,ymax=0),size=0.1,alpha=0.3)+
               theme_classic()+
               scale_y_continuous(breaks=NULL)+
               xlim(0,5.2e7)+ylab("II")+xlab("")

chr2_blank <- ggplot()+geom_blank()+xlim(0,5.2e7)+theme_classic()+
               theme(axis.title.x=element_blank(),
               axis.text.x=element_blank(),
               axis.ticks.x=element_blank(),
               axis.line = element_blank())

chr3_div_plot <- ggplot()+
               geom_point(aes(chr3_single_pi$BIN_START,chr3_single_pi$PI),alpha=0.3,colour="blue",size=0.5)+
               geom_smooth(aes(chr3_single_pi$BIN_START,chr3_single_pi$PI),span=0.1,colour="blue")+
               geom_point(aes(chr3_shared_pi$BIN_START,chr3_shared_pi$PI),alpha=0.3,colour="red",size=0.5)+
               geom_smooth(aes(chr3_shared_pi$BIN_START,chr3_shared_pi$PI),span = 0.1,colour="red")+
               xlim(0,5.2e7)+ylim(0,0.008)+
               xlab("")+ylab(expression(pi))+
               theme_classic()+theme(axis.title.x=element_blank(),
               axis.text.x=element_blank(),
               axis.ticks.x=element_blank())

chr3_fst_plot <- ggplot()+
               geom_rect(aes(xmin=0,ymin=quantile[1],xmax=5.2e7,ymax=quantile[2]),fill="grey90")+
               geom_point(aes(chr3_fst$BIN_START,chr3_fst$WEIGHTED_FST),colour="red",alpha=0.3,size=0.5)+
               theme_classic()+
               xlim(0,5.2e7)+ylim(0,0.4)+
               xlab("")+ylab("Fst")+
               theme(axis.title.x=element_blank(),
               axis.text.x=element_blank(),
               axis.ticks.x=element_blank())

chr3_gene_plot <- ggplot()+
               geom_rect(aes(xmin=1,ymin=-1,xmax=max(chr3$V2),ymax=1),fill="#fee090",alpha=0.5)+
               geom_rect(aes(xmin=chr3N$V2,ymin=-1,xmax=chr3N$V3,ymax=1),fill="white",linetype=0)+
               geom_linerange(aes(chr3pos$V2,ymin=0,ymax=1),size=0.1,alpha=0.3)+
               geom_linerange(aes(chr3neg$V2,ymin=-1,ymax=0),size=0.1,alpha=0.3)+
               theme_classic()+ scale_y_continuous(breaks=NULL)+
               xlim(0,5.2e7)+ylab("III")+xlab("")

chr3_blank <- ggplot()+geom_blank()+xlim(0,5.2e7)+theme_classic()+
               theme(axis.title.x=element_blank(),
               axis.text.x=element_blank(),
               axis.ticks.x=element_blank(),
               axis.line = element_blank())


chr4_div_plot <- ggplot()+
               geom_point(aes(chr4_single_pi$BIN_START,chr4_single_pi$PI),alpha=0.3,colour="blue",size=0.5)+
               geom_smooth(aes(chr4_single_pi$BIN_START,chr4_single_pi$PI),span=0.1,colour="blue")+
               geom_point(aes(chr4_shared_pi$BIN_START,chr4_shared_pi$PI),alpha=0.3,colour="red",size=0.5)+
               geom_smooth(aes(chr4_shared_pi$BIN_START,chr4_shared_pi$PI),span = 0.1,colour="red")+
               xlim(0,5.2e7)+ylim(0,0.008)+
               xlab("")+ylab(expression(pi))+
               theme_classic()+theme(axis.title.x=element_blank(),
               axis.text.x=element_blank(),
               axis.ticks.x=element_blank())

chr4_fst_plot <- ggplot()+
               geom_rect(aes(xmin=0,ymin=quantile[1],xmax=5.2e7,ymax=quantile[2]),fill="grey90")+
               geom_point(aes(chr4_fst$BIN_START,chr4_fst$WEIGHTED_FST),colour="red",alpha=0.3,size=0.5)+
               theme_classic()+
               xlim(0,5.2e7)+ylim(0,0.4)+
               xlab("")+ylab("Fst")+
               theme(axis.title.x=element_blank(),
               axis.text.x=element_blank(),
               axis.ticks.x=element_blank())
chr4_gene_plot <- ggplot()+
               geom_rect(aes(xmin=1,ymin=-1,xmax=max(chr4$V2),ymax=1),fill="#d1e5f0",alpha=0.5)+
               geom_rect(aes(xmin=chr4N$V2,ymin=-1,xmax=chr4N$V3,ymax=1),fill="white",linetype=0)+
               geom_linerange(aes(chr4pos$V2,ymin=0,ymax=1),size=0.1,alpha=0.3)+
               geom_linerange(aes(chr4neg$V2,ymin=-1,ymax=0),size=0.1,alpha=0.3)+
               theme_classic()+ scale_y_continuous(breaks=NULL)+
               xlim(0,5.2e7)+ylab("IV")+xlab("")

chr4_blank <- ggplot()+geom_blank()+xlim(0,5.2e7)+theme_classic()+
               theme(axis.title.x=element_blank(),
               axis.text.x=element_blank(),
               axis.ticks.x=element_blank(),
               axis.line = element_blank())


chr5_div_plot <- ggplot()+
               geom_point(aes(chr5_single_pi$BIN_START,chr5_single_pi$PI),alpha=0.3,colour="blue",size=0.5)+
               geom_smooth(aes(chr5_single_pi$BIN_START,chr5_single_pi$PI),span=0.1,colour="blue")+
               geom_point(aes(chr5_shared_pi$BIN_START,chr5_shared_pi$PI),alpha=0.3,colour="red",size=0.5)+
               geom_smooth(aes(chr5_shared_pi$BIN_START,chr5_shared_pi$PI),span = 0.1,colour="red")+
               xlim(0,5.2e7)+ylim(0,0.008)+
               xlab("")+ylab(expression(pi))+
               theme_classic()+theme(axis.title.x=element_blank(),
               axis.text.x=element_blank(),
               axis.ticks.x=element_blank())

chr5_fst_plot <- ggplot()+
               geom_rect(aes(xmin=0,ymin=quantile[1],xmax=5.2e7,ymax=quantile[2]),fill="grey90")+
               geom_point(aes(chr5_fst$BIN_START,chr5_fst$WEIGHTED_FST),colour="red",alpha=0.3,size=0.5)+
               theme_classic()+
               xlim(0,5.2e7)+ylim(0,0.4)+
               xlab("")+ylab("Fst")+
               theme(axis.title.x=element_blank(),
               axis.text.x=element_blank(),
               axis.ticks.x=element_blank())

chr5_gene_plot <- ggplot()+
               geom_rect(aes(xmin=1,ymin=-1,xmax=max(chr5$V2),ymax=1),fill="#67a9cf",alpha=0.5)+
               geom_rect(aes(xmin=chr5N$V2,ymin=-1,xmax=chr5N$V3,ymax=1),fill="white",linetype=0)+
               geom_linerange(aes(chr5pos$V2,ymin=0,ymax=1),size=0.1,alpha=0.3)+
               geom_linerange(aes(chr5neg$V2,ymin=-1,ymax=0),size=0.1,alpha=0.3)+
               theme_classic()+ scale_y_continuous(breaks=NULL)+
               xlim(0,5.2e7)+ylab("V")+xlab("")

chr5_blank <- ggplot()+geom_blank()+xlim(0,5.2e7)+theme_classic()+
               theme(axis.title.x=element_blank(),
               axis.text.x=element_blank(),
               axis.ticks.x=element_blank(),
               axis.line = element_blank())

chrX_div_plot <- ggplot()+
               geom_point(aes(chrX_single_pi$BIN_START,chrX_single_pi$PI),alpha=0.3,colour="blue",size=0.5)+
               geom_smooth(aes(chrX_single_pi$BIN_START,chrX_single_pi$PI),span=0.1,colour="blue")+
               geom_point(aes(chrX_shared_pi$BIN_START,chrX_shared_pi$PI),alpha=0.3,colour="red",size=0.5)+
               geom_smooth(aes(chrX_shared_pi$BIN_START,chrX_shared_pi$PI),span = 0.1,colour="red")+
               xlim(0,5.2e7)+ylim(0,0.008)+
               xlab("")+ylab(expression(pi))+
               theme_classic()+theme(axis.title.x=element_blank(),
               axis.text.x=element_blank(),
               axis.ticks.x=element_blank())

chrX_fst_plot <- ggplot()+
               geom_rect(aes(xmin=0,ymin=quantile[1],xmax=5.2e7,ymax=quantile[2]),fill="grey90")+
               geom_point(aes(chrX_fst$BIN_START,chrX_fst$WEIGHTED_FST),colour="red",alpha=0.3,size=0.5)+
               theme_classic()+
               xlim(0,5.2e7)+ylim(0,0.4)+
               xlab("")+ylab("Fst")+
               theme(axis.title.x=element_blank(),
               axis.text.x=element_blank(),
               axis.ticks.x=element_blank())

chrX_gene_plot <- ggplot()+
               geom_rect(aes(xmin=1,ymin=-1,xmax=max(chrX$V2),ymax=1),fill="#4575b4",alpha=0.5)+
               geom_rect(aes(xmin=chrXN$V2,ymin=-1,xmax=chrXN$V3,ymax=1),fill="white",linetype=0)+
               geom_linerange(aes(chrXpos$V2,ymin=0,ymax=1),size=0.1,alpha=0.3)+
               geom_linerange(aes(chrXneg$V2,ymin=-1,ymax=0),size=0.1,alpha=0.3)+
               theme_classic()+ scale_y_continuous(breaks=NULL)+
               xlim(0,5.2e7)+ylab("X")+xlab("")


#patchwork
chr1_div_plot +
chr1_fst_plot +
chr1_gene_plot +
#chr1_blank +
chr2_div_plot +
chr2_fst_plot +
chr2_gene_plot +
#chr2_blank +
chr3_div_plot +
chr3_fst_plot +
chr3_gene_plot +
#chr3_blank +
chr4_div_plot +
chr4_fst_plot +
chr4_gene_plot +
#chr4_blank +
chr5_div_plot +
chr5_fst_plot +
chr5_gene_plot +
#chr5_blank +
chrX_div_plot +
chrX_fst_plot +
chrX_gene_plot +
plot_layout(ncol=1)

ggsave("genome_diversity.pdf")
```






# R-3.4.0
# library(ggplot2)
# library(patchwork)		
# 			
# a<-read.table("chr1.singletons.pi.windowed.pi",header=T)			
# b<-read.table("chr1.shared.pi.windowed.pi",header=T)			
# c <- read.table("chr1_fst_100k_allpop.windowed.weir.fst",header=T)			
# 			
# bz <- read.table("../../../../XQTL/04_VARIANTS/XQTL_BZ/XQTL_BZ.merged.fst",header=F)
# bz_1 <- bz[bz$V1=="hcontortus_chr1_Celeg_TT_arrow_pilon",]
#
# 			
# div_plot <- ggplot()+
# 				geom_point(aes(a$BIN_START,a$PI),alpha=0.3,colour="blue")+
# 				geom_smooth(aes(a$BIN_START,a$PI),span=0.1,colour="blue")+
# 				geom_point(aes(b$BIN_START,b$PI),alpha=0.3,colour="red")+
# 				geom_smooth(aes(b$BIN_START,b$PI),span = 0.1,colour="red")+
# 				theme_bw()
# fst_plot <- ggplot()+geom_vline(xintercept = 7029311,linetype = "longdash")+geom_point(aes(c$BIN_START,c$WEIGHTED_FST),colour="red",alpha=0.3)+theme_bw()
# bz_plot <- ggplot()+geom_vline(xintercept = 7029311,linetype = "longdash")+geom_point(aes(bz_1$V2,bz_1$V49),colour="red",alpha=0.3)+theme_bw()		
# 			
#
#
# div_plot + fst_plot + bz_plot + plot_layout(ncol = 1)


```shell
cd /nfs/users/nfs_s/sd21/lustre118_link/play
cut -f1 caenorhabditis_elegans.PRJNA13758.WBPS11.annotations.gff3 | grep -v "#" |  sort | uniq | while read -r CHR; do awk -v CHR=$CHR '{if($1==CHR && $3=="gene" && $2=="WormBase") print $1,($5+$4)/2,$7,$9}' OFS="\t" caenorhabditis_elegans.PRJNA13758.WBPS11.annotations.gff3 | grep "protein_coding" > ${CHR}.genepos.list;  done
```

```R
R-3.5.0
library(ggplot2)
library(patchwork)

chr1<-read.table("I.genepos.list",header=F)
chr2<-read.table("II.genepos.list",header=F)
chr3<-read.table("III.genepos.list",header=F)
chr4<-read.table("IV.genepos.list",header=F)
chr5<-read.table("V.genepos.list",header=F)
chrX<-read.table("X.genepos.list",header=F)

chr1pos<-chr1[chr1$V3=="+",]
chr1neg<-chr1[chr1$V3=="-",]
chr2pos<-chr2[chr2$V3=="+",]
chr2neg<-chr2[chr2$V3=="-",]
chr3pos<-chr3[chr3$V3=="+",]
chr3neg<-chr3[chr3$V3=="-",]
chr4pos<-chr4[chr4$V3=="+",]
chr4neg<-chr4[chr4$V3=="-",]
chr5pos<-chr5[chr5$V3=="+",]
chr5neg<-chr5[chr5$V3=="-",]
chrXpos<-chrX[chrX$V3=="+",]
chrXneg<-chrX[chrX$V3=="-",]



chr1_gene_plot	<-	ggplot()+
               geom_rect(aes(xmin=1,ymin=-1,xmax=max(chr1$V2),ymax=1),fill="white")+
               geom_linerange(aes(chr1pos$V2,ymin=0,ymax=1),size=0.1,alpha=0.2)+
               geom_linerange(aes(chr1neg$V2,ymin=-1,ymax=0),size=0.1,alpha=0.2)+
               theme_classic()+ scale_y_continuous(breaks=NULL)+
               xlim(0,2.2e7)+ylab("I")+xlab("")


chr2_gene_plot	<-	ggplot()+
               geom_rect(aes(xmin=1,ymin=-1,xmax=max(chr2$V2),ymax=1),fill="white")+
               geom_linerange(aes(chr2pos$V2,ymin=0,ymax=1),size=0.1,alpha=0.2)+
               geom_linerange(aes(chr2neg$V2,ymin=-1,ymax=0),size=0.1,alpha=0.2)+
               theme_classic()+ scale_y_continuous(breaks=NULL)+
               xlim(0,2.2e7)+ylab("II")+xlab("")


chr3_gene_plot	<-	ggplot()+
               geom_rect(aes(xmin=1,ymin=-1,xmax=max(chr3$V2),ymax=1),fill="white")+
               geom_linerange(aes(chr3pos$V2,ymin=0,ymax=1),size=0.1,alpha=0.2)+
               geom_linerange(aes(chr3neg$V2,ymin=-1,ymax=0),size=0.1,alpha=0.2)+
               theme_classic()+ scale_y_continuous(breaks=NULL)+
               xlim(0,2.2e7)+ylab("III")+xlab("")


chr4_gene_plot	<-	ggplot()+
               geom_rect(aes(xmin=1,ymin=-1,xmax=max(chr4$V2),ymax=1),fill="white")+
               geom_linerange(aes(chr4pos$V2,ymin=0,ymax=1),size=0.1,alpha=0.2)+
               geom_linerange(aes(chr4neg$V2,ymin=-1,ymax=0),size=0.1,alpha=0.2)+
               theme_classic()+ scale_y_continuous(breaks=NULL)+
               xlim(0,2.2e7)+ylab("IV")+xlab("")


chr5_gene_plot	<-	ggplot()+
               geom_rect(aes(xmin=1,ymin=-1,xmax=max(chr5$V2),ymax=1),fill="white")+
               geom_linerange(aes(chr5pos$V2,ymin=0,ymax=1),size=0.1,alpha=0.2)+
               geom_linerange(aes(chr5neg$V2,ymin=-1,ymax=0),size=0.1,alpha=0.2)+
               theme_classic()+ scale_y_continuous(breaks=NULL)+
               xlim(0,2.2e7)+ylab("V")+xlab("")

chrX_gene_plot	<-	ggplot()+
               geom_rect(aes(xmin=1,ymin=-1,xmax=max(chrX$V2),ymax=1),fill="white")+
               geom_linerange(aes(chrXpos$V2,ymin=0,ymax=1),size=0.1,alpha=0.2)+
               geom_linerange(aes(chrXneg$V2,ymin=-1,ymax=0),size=0.1,alpha=0.2)+
               theme_classic()+ scale_y_continuous(breaks=NULL)+
               xlim(0,2.2e7)+ylab("X")+xlab("")

chr1_gene_plot +
chr2_gene_plot +
chr3_gene_plot +
chr4_gene_plot +
chr5_gene_plot +
chrX_gene_plot +
plot_layout(ncol=1)

```
