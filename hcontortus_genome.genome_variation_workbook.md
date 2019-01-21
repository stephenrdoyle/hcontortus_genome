# Haemonchus genome - genome variation analyses

## Table of contents

1.
2.
3. [Genome wide nucleotide diversity analysis](#nuc_div)













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
               geom_rect(aes(xmin=1,ymin=-1,xmax=max(chr1$V2),ymax=1),fill="#b2182b")+
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
               geom_rect(aes(xmin=1,ymin=-1,xmax=max(chr2$V2),ymax=1),fill="#fc8d59")+
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
               geom_rect(aes(xmin=1,ymin=-1,xmax=max(chr3$V2),ymax=1),fill="#fee090")+
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
               geom_rect(aes(xmin=1,ymin=-1,xmax=max(chr4$V2),ymax=1),fill="#d1e5f0")+
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
               geom_rect(aes(xmin=1,ymin=-1,xmax=max(chr5$V2),ymax=1),fill="#67a9cf")+
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
               geom_rect(aes(xmin=1,ymin=-1,xmax=max(chrX$V2),ymax=1),fill="#4575b4")+
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
