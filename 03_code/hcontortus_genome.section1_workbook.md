# Haemonchus genome paper
## Section 1: Chromosome structure of Haemonchus contortus

1. [Genome assembly](#genome)
2. [Completion of X chromosome](#xchromosome)
3. [Genome polishing](#polishing)
4. [Chromosomal](#chromosomesynteny)
     * [Figure 1a](#figure1a)
     * [Supplementary Figure 3b](#figureS2b)
5. [Microsynteny](#microsynteny)
     * [Figure 1b](#figure1b)
     * [Supplementary Figure 3b](#figureS3b)
     * [Supplementary Figure 3c](#figureS3c)
     * [Supplementary Figure 3d](#figureS3d)
     * [Supplementary Figure 3e](#figureS3e)
8. [Genome completeness - CEGMA & BUSCO](#cegmabusco)
     * [Supplementary Figure 4](#figureS4)
7. [Repeats](#repeats)
     * [Supplementary Figure 5](#figureS5)



<a name="figure1b"></a>
******
## Genome assembly <a name="genome"></a>
From the paper:

Initial manual improvement on the V1 genome focused on iterative scaffolding with SSPACE (Boetzer et al., 2011) and gap-filling with IMAGE (Tsai et al., 2010) using Illumina 500 bp and 3 kbp libraries, with additional low coverage data from 3, 8 and 20 kbp libraries generated using Roche 454 sequencing. These improvements were performed alongside breaking of discordant joins using Reapr (Hunt et al., 2013), and visual inspection using Gap5 (Bonfield and Whitwham, 2010). Substantial genetic variation was present in the sequence data due to the sequencing of DNA derived from a pool of individuals, resulting in a high frequency of haplotypes that assembled separately and therefore present as multiple copies of unique regions in the assembly. We surmised that much of the assembly fragmentation was due to the scaffolding tools not being able to deal with the changing rates of haplotypic variation so we attempted to solve this manually in gap5. We were aware that we did not have sufficient information to correctly phase these haplotypes, so instead, we chose the longest scaffold paths available, accepting that our scaffolds would not represent single haplotypes but would rather be an amalgamation of haplotypes representing single chromosomal regions. This approach was initially difficult and time-consuming, and was further confounded by a large number of repetitive sequences present in each haplotype.

Significant improvements in scaffold length were gained by the integration of OpGen (http://www.opgen.com/) optical mapping data. Optical mapping was performed following methods described previously (Tsai et al., 2013) with the following exceptions: high molecular weight DNA was prepared by proteinase K lysis of pools of ~500 H. contortus L3 embedded in agarose plugs, after which one of three restriction enzymes (KpnI, AflII and KpnI) were used to generate three separate restriction map datasets. Initial attempts to generate a de novo assembly using optical molecules alone was unsuccessful, and therefore, optical contigs were generated using DNA sequence seeds from the genome assembly using GenomeBuilder (OpGen) and visualised and manually curated using AssemblyViewer (OpGen). Although this approach was successful, it was limited by the quality and integrity of the gap-dense scaffolds and arbitrary nature of the haplotype scaffolding.

Subsequent integration of PacBio long-read data alongside the optical mapping data resulted in major increases in contiguity. PacBio sequencing libraries were prepared and sequenced using the PacBio RSII system. A total of 32.3 Gbp raw subreads (n = 4,085,541, N50 = 10,299 bp) were generated from 33 flow cells, representing approximately 133.8× coverage (based on the estimated genome size of 283 Mbp). A de novo PacBio assembly was generated using Sprai (v0.9.9.18; http://zombie.cb.k.u-tokyo.ac.jp/sprai/index.html), which was mapped to the assembly and used to manually resolve many gaps and resolve some of the phasing between haplotypes. Using the Sprai PacBio de novo assembly, we were also able to incorporate contigs that were previously missing from the working assembly. The increase in contiguity of the PacBio assemblies, further improved using canu v1.3 (Koren et al., 2017), revealed two major but diverse haplotype groups segregating at approximately 65% and 30% frequency in the pooled individuals sequenced; the presence of such diverse haplotypes resulted in a significantly expanded assembly over 500 Mbp. The major haplotype was more contiguous and therefore was chosen as the primary sequence to incorporate and improve the assembly. This approach was supported by competitive mapping of a single worm short read sequencing library (ENA: ERS086777), which was found to preferentially map to different scaffolds and/or different contigs within scaffolds that we inferred were different haplotypes in this single worm. Regions in the chromosome with no ERS086777 coverage, but for which an alternate haplotype was identified with coverage from the PacBio assembly, were manually removed from the main chromosomes and replaced. Further, this selective mapping correlated well with the optical contigs, and once these sequences were removed, much better optical alignments were produced further improving the assembly. Alternate haplotypes from the PacBio assembly, and those removed from the main assembly, were analysed as a separate haplotype assembly.

The increase in contiguity and resolution of shorter repetitive regions using PacBio began to reveal chromosome-specific repetitive units. Although these repeats were highly collapsed in the assembly and were typically not spanned by optical molecules, we were able to iteratively identify and join contigs/scaffolds flanking large tandem repeats that had clear read-pair evidence that they only occurred in a single location in the genome (i.e. read pairs were all mapped locally once the join was made). These were further supported by local assemblies of subsets of PacBio reads that contained copies of the repeat regions followed by de novo assembly using canu (Koren et al., 2017) to reconstruct the flanking unique regions surrounding a repeat. This iterative process resulted in the production of chromosome-scale scaffolds, each terminating with a 6 bp repeat consistent with being a telomeric sequence (sequence motif: TTAGGC).


[↥ **Back to top**](#top)
******
## Completion of the X chromosome <a name="xchromosome"></a>
The V3 version of the genome consisted of assembled autosomes, however, the X chromosome was still in pieces that to date were not resolved, althoguh we had a good idea of what pieces were left in the assembly that belonged in the X chromosome based on male vs female coverage.

>hcontortus_6_chrX_Celeg_Ta  
>hcontortus_7_chrX_Celeg_Tb  

>CONTAM_scf7180000021165  
>CONTAM_scf7180000022481  
>CONTAM_scf7180000025901#1  
>U_alt_001073  
>U_Celeg_1_Poss_join_XTb  
>U_Contig4229_quiver  
>U_Contig986_quiver  
>U_contig_no_celeg_promer_poss_join_XTa  



```shell
#
working dir: /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/FIX_XCHR

# get barcodes 10X reads
barcoded.fastq.gz -> /lustre/scratch118/infgen/team133/sd21/hc/10X_genomics/chromium/hc_inbred/hc_inbred/outs/barcoded.fastq.gz

# note - to get the barcoded.fastq.gz, needed to run 10X longranger , eg.
longranger-2.1.2/longranger basic --id=hc_inbred --fastqs=/nfs/users/nfs_s/sd21/lustre118_link/hc/10X_genomics/chromium/hc_inbred --readgroup=hc_inbred --localcores=4

# map reads
bwa mem \
     -t 30 xchr_TaTb_plus_bits.sized.renamed.fa \
     -p CHROMIUM_interleaved.fastq |\
     samtools-1.3 view -Sb - |\
     samtools-1.3 sort -n -o ./CHROMIUM-sorted.bam -


# make a bam list from the mapped reads
ls -1 CHROMIUM-sorted.bam > bam.list

# run ARCS
arcs -f xchr_TaTb_plus_bits.sized.renamed.fa -a bam.list -s 95 -c 3 -l 0 -z 500 -m 5-10000 -d 0 -e 50000 -r 0.05 -i 16 -v 1


# run LINKS
LINKS -f xchr_TaTb_plus_bits.sized.renamed.fa -s empty.fof -k 15 -b xchr_TaTb_plus_bits.sized.renamed.fa.scaff_s95_c3_l0_d0_e50000_r0.05_original -l 1 -t 2 -a 0.9 -v 1
```



[↥ **Back to top**](#top)
******
## Genome Polishing <a name="polishing"></a>
### Arrow
Used Shane McCarthy's run-arrow

```shell
# run Arrow
prepend_path PERL5LIB /software/vertres/runner/modules
prepend_path PATH /software/vertres/runner/scripts

run-arrow +loop 60 +mail sd21 +maxjobs 1000 +retries 4 -a HAEM_V4.fa -f /nfs/users/nfs_s/sd21/lustre118_link/hc/raw/pacbio/bams.fofn -o ./POLISH_ARROW_HcV4
```

### Pilon
Used the output of arrow as input for pilon
```shell
# map reads from single worm to reference

~sd21/bash_scripts/run_bwamem_splitter Pilon_prep_singlefemale2HaemV4 $PWD/HAEM_V4_arrow.fa $PWD/single_adult_female_19220_merge_R1.fq $PWD/single_adult_female_19220_merge_R2.fq

# run Pilon
working dir: /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/POLISH/POLISH_ARROW_HcV4_haplo_contam


java -Xmx200G -jar pilon-1.22/pilon-1.22.jar \   
     --genome HAEM_V4_arrow.fa \
     --bam Pilon_prep_singlefemale2HaemV4.merged.sorted.marked.bam \
     --fix all \
     --verbose \
     --diploid \
     --changes \
     --chunksize 1000000 \
     --nostrays \
     --threads 30
```



 [↥ **Back to top**](#top)
 ******
## Chromosomal synteny <a name="chromosomesynteny"></a>
### CIRCOS <a name="figure1a"></a>
Using circos to highlight broad chromosomal similarities between Haemonchus and Celegans
```shell   
working dir: /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/CIRCOS

# get data
ln -s ../REF/HAEM_V4_final.chr.fa
ln -s ../../../REFERENCE_SEQUENCES/caenorhabditis_elegans/caenorhabditis_elegans.PRJNA13758.WBPS9.genomic_softmasked.fa

# run promer
promer --mum --prefix Ce_V_HcV4 caenorhabditis_elegans.PRJNA13758.WBPS9.genomic_softmasked.fa HAEM_V4_final.chr.fa

# run James' nucmer to circos params script
perl /nfs/users/nfs_j/jc17/bin/Nucmer.2.circos.pl --promer --ref_order=relw --debug --ribbons --coord-file Ce_V_HcV4.coords --flipquery --no_ref_labels --no_query_labels --colour_links_by_query --min_chr_len=20000 --min_hit_len=10000

# run circos
perl /nfs/users/nfs_j/jc17/software/circos-0.67-pre5/bin/circos

# this generated the circos plot used in Figure 1a. Note that I modified in conf file manually to include the correct chromosome colours as shown in the manuscript
```

### male and female genome coverage_plot - used in Supplementary Figure 2b <a name="figureS2b"></a>
cd ~/lustre118_link/hc/GENOME/POPULATION_DIVERSITY/GENOME_COVERAGE

```R
# load libraries
library(ggplot2)
library(patchwork)

# load data
male<-read.table("GB_ISEN1_001.merged.100000_window.cov",header=F)
male<-male[male$V1!="hcontortus_chr_mtDNA_arrow_pilon",]
female<-read.table("GB_ISEN1_006.merged.100000_window.cov",header=F)
female<-female[female$V1!="hcontortus_chr_mtDNA_arrow_pilon",]

# set colours
chr_colours<-c("#b2182b","#fc8d59","#fee090","#d1e5f0","#67a9cf","#4575b4")

# plot
plot_female<-ggplot(female,aes(V2/10^6,log10(V5),col=V1))+
	geom_point(alpha=0.8)+
	scale_colour_manual(values=chr_colours, guide = FALSE)+
	facet_grid(.~V1)+
	labs(title="Sample: MHco3(ISE).N1_006 - female XX",x="Genomic coordinate (Mbp)",y="Coverage per 100 kbp window (log10)")+
	theme_bw()+
	ylim(0.5,2.25)

plot_male<-ggplot(male,aes(V2/10^6,log10(V5),col=V1))+
	geom_point(alpha=0.8)+
	scale_colour_manual(values=chr_colours, guide = FALSE)+
	facet_grid(.~V1)+
	labs(title="Sample: MHco3(ISE).N1_001 - male XO",x="Genomic coordinate (Mbp)",y="Coverage per 100 kbp window (log10)")+
	theme_bw()+
	ylim(0.5,3)

# put it together
plot_female + plot_male + plot_layout(ncol=1)

ggsave("male_v_female_genomecoverage.pdf", useDingbats = FALSE)
```



[↥ **Back to top**](#top)


******
## Microsynteny <a name="microsynteny"></a>

### Compare conservation of orthologs between Ce and Hc chromosomes
```bash
cd /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/ORTHOLOGY/SYNTENY


ln -s /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/SELECTION/PROTEIN_FASTAs/Results_Jan25/Orthologues_Jan25/Orthologues/Orthologues_hc_V4.proteins.unique/hc_V4.proteins.unique__v__ce.proteins.unique.csv

ln -fs /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/TRANSCRIPTOME/TRANSCRIPTOME_CURATION/HCON_V4_WBP11plus_190125.ips.gff3 hc.gff3
ln -fs /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/TRANSCRIPTOME/TRANSCRIPTOME_SUMMARY_STATS/CELEGANS/wormbase.20240.complete.gff3 ce.gff3


# make lists of 1 to 1 orthologs
awk 'NF==3 {print $1,$2,$3}' OFS="\t" hc_V4.proteins.unique__v__hp.proteins.unique.csv > hc_hp_1to1.list
#> 9970 hc_hp_1to1.list

awk 'NF==3 {print $2,$3}' OFS="\t" hc_V4.proteins.unique__v__ce.proteins.unique.csv > hc_ce_1to1.list
#> 7361 hc_ce_1to1.list

# get gene coordinates
awk -F '[\t;]' '$3=="gene" {print $1,$4,$5,$7,$10}' OFS="\t" ce.gff3 | sed 's/Name=//g' > ce.gff.coords
awk -F '[\t;]' '$3=="mRNA" {print $1,$4,$5,$7,$9}' OFS="\t" hc.gff3 | sed 's/ID=//g' > hc.gff.coords

# extract positional data for both Hc and Ce
while read HC_ID CE_ID; do COORDS1=$( grep ${HC_ID} hc.gff.coords); COORDS2=$( grep ${CE_ID} ce.gff.coords); echo -e "${COORDS1}\t${COORDS2}"; done <hc_ce_1to1.list > hc_ce_1to1.hccoords
sort -k1,1 -k2,2n hc_ce_1to1.hccoords | awk '{if($2+0==$2 && $7+0==$7) print $0}' > hc_ce_1to1.coords.chrsorted

# count the number of concordant vs non concordant chromosomal orthologs
cut -f1,6 hc_ce_1to1.coords.chrsorted | sort | uniq -c | awk '{print $2,$3,$1}' OFS="\t" > hc_ce_orthologcount_perchr.txt
```



### Make the plot - Figure 1 b <a name="figure1b"></a>
```R
# load libraries
library(ggplot2)
library(patchwork)
library(dplyr)

# load data
data <- read.table("hc_ce_orthologcount_perchr.txt",header=F)
data.1 <- data %>% group_by(V1) %>%  mutate(per=paste0(round(V3/sum(V3)*100, 2), "%")) %>% ungroup

# plot number of genes per chromosome
plot1<- ggplot(data,aes(data$V2,data$V3,group=data$V1))+
	geom_bar(aes(fill=data$V2), stat = "identity", position = "dodge")+
	theme_bw()+
	labs(x="C .elegans chromosome name",y="1:1 ortholog count")+
	facet_grid(.~data$V1)

	ggplot(data.1,aes(data.1$V2,data.1$V3,group=data.1$V1))+
	geom_bar(aes(fill=data.1$V2), stat = "identity", position = "dodge")+
	geom_text(aes(label=data.1$per), position=position_dodge(width=0.9), vjust=-0.25,size=1)+
	theme_bw()+
	labs(x="C .elegans chromosome name",y="1:1 ortholog count")+
	facet_grid(.~data.1$V1)

# plot scatter between hc and ce orthologs by position, by hc chromosome, coloured by ce chromosome
data2<-read.table("hc_ce_1to1.coords.chrsorted",header=F)
plot2<-ggplot(data2,aes(data2$V2,data2$V7,group=data2$V1))+
	geom_point(aes(col=data2$V6),alpha=0.5,size=0.5)+
	facet_grid(.~data2$V1)+
	labs(x="H. contortus chromosome position (bp)",y="C. elegans chromosome position (bp)")+
	theme_bw()

plot1 + plot2 + plot_layout(ncol=1)
ggsave("hc_ce_orthology.pdf",useDingbats=F)
```




### compare distance between pairs of orthologs between Ce and Hc
```bash
awk '{if($1=="hcontortus_chr1_Celeg_TT_arrow_pilon" && $6=="I") print $0}' hc_ce_1to1.coords.chrsorted > hc_ce_1to1.coords.chr1.sorted
awk '{if($1=="hcontortus_chr2_Celeg_TT_arrow_pilon" && $6=="II") print $0}' hc_ce_1to1.coords.chrsorted > hc_ce_1to1.coords.chr2.sorted
awk '{if($1=="hcontortus_chr3_Celeg_TT_arrow_pilon" && $6=="III") print $0}' hc_ce_1to1.coords.chrsorted > hc_ce_1to1.coords.chr3.sorted
awk '{if($1=="hcontortus_chr4_Celeg_TT_arrow_pilon" && $6=="IV") print $0}' hc_ce_1to1.coords.chrsorted > hc_ce_1to1.coords.chr4.sorted
awk '{if($1=="hcontortus_chr5_Celeg_TT_arrow_pilon" && $6=="V") print $0}' hc_ce_1to1.coords.chrsorted > hc_ce_1to1.coords.chr5.sorted
awk '{if($1=="hcontortus_chrX_Celeg_TT_arrow_pilon" && $6=="X") print $0}' hc_ce_1to1.coords.chrsorted > hc_ce_1to1.coords.chrX.sorted

#
paste <(head -n -1 hc_ce_1to1.coords.chr1.sorted) <(tail -n +2 hc_ce_1to1.coords.chr1.sorted) > hc_ce_1to1.coords.chr1.sorted.offset
paste <(head -n -1 hc_ce_1to1.coords.chr2.sorted) <(tail -n +2 hc_ce_1to1.coords.chr2.sorted) > hc_ce_1to1.coords.chr2.sorted.offset
paste <(head -n -1 hc_ce_1to1.coords.chr3.sorted) <(tail -n +2 hc_ce_1to1.coords.chr3.sorted) > hc_ce_1to1.coords.chr3.sorted.offset
paste <(head -n -1 hc_ce_1to1.coords.chr4.sorted) <(tail -n +2 hc_ce_1to1.coords.chr4.sorted) > hc_ce_1to1.coords.chr4.sorted.offset
paste <(head -n -1 hc_ce_1to1.coords.chr5.sorted) <(tail -n +2 hc_ce_1to1.coords.chr5.sorted) > hc_ce_1to1.coords.chr5.sorted.offset
paste <(head -n -1 hc_ce_1to1.coords.chrX.sorted) <(tail -n +2 hc_ce_1to1.coords.chrX.sorted) > hc_ce_1to1.coords.chrX.sorted.offset

cat *.offset > all.offset



#Dist between Hc1-Hc2 =  $12+($13-$12)  - $2+($3-$2)
#Direction of Hc1-Hc2 =  $4"_"$14
#
#Dist between Ce1-Ce2 =  $17+($18-$17)  - $7+($8-$7)
#Direction of Ce1-Ce2 =  $9"_"$19

for i in *.offset; do \
distSP1=$( awk '{print ($12+($13-$12))-($2+($3-$2))}' ${i} )
dirSP1=$(awk '{print $4"_"$14}' ${i})
distSP2=$( awk '{print ($17+($18-$17))-($7+($8-$7))}' ${i} )
dirSP2=$( awk '{print $9"_"$19 }' ${i} )
dir=$( paste <(echo "$dirSP1") <(echo "$dirSP2") | while read -r dirSP1 dirSP2 ; do if [ $dirSP1 == $dirSP2 ] ; then echo YES ; else echo NO ; fi ;  done )
paste <(echo "$distSP1") <(echo "$dirSP1") <(echo "$distSP2") <(echo "$dirSP2")  <(echo "$dir") > ${i%.offset}.plotdata;  done




# EXAMPLE - where orthologs are intentionally on difference chromosomes
awk '{if($1=="hcontortus_chr1_Celeg_TT_arrow_pilon" && $6=="II") print $0}' hc_ce_1to1.coords.chrsorted > mismatch_hc1_ce2
paste <(head -n -1 mismatch_hc1_ce2) <(tail -n +2 mismatch_hc1_ce2) > mismatch_hc1_ce2.offset

for i in mismatch_hc1_ce2.offset; do \
distHc=$( awk '{print ($12+($13-$12))-($2+($3-$2))}' ${i} )
dirHc=$(awk '{print $4"_"$14}' ${i})
distCe=$( awk '{print ($17+($18-$17))-($7+($8-$7))}' ${i} )
dirCe=$( awk '{print $9"_"$19 }' ${i} )
dir=$( paste <(echo "$dirHc") <(echo "$dirCe") | while read -r hcdir cedir ; do if [ $hcdir == $cedir ] ; then echo YES ; else echo NO ; fi ;  done )
paste <(echo "$distHc") <(echo "$dirHc") <(echo "$distCe") <(echo "$dirCe")  <(echo "$dir") > mismatch_hc1_ce2.plotdata;  done


# EXAMPLE - generate inputdata for all chromosomes combined
for i in all.offset; do \
distHc=$( awk '{print ($12+($13-$12))-($2+($3-$2))}' ${i} )
dirHc=$(awk '{print $4"_"$14}' ${i})
distCe=$( awk '{print ($17+($18-$17))-($7+($8-$7))}' ${i} )
dirCe=$( awk '{print $9"_"$19 }' ${i} )
dir=$( paste <(echo "$dirHc") <(echo "$dirCe") | while read -r hcdir cedir ; do if [ $hcdir == $cedir ] ; then echo YES ; else echo NO ; fi ;  done )
paste <(echo "$distHc") <(echo "$dirHc") <(echo "$distCe") <(echo "$dirCe")  <(echo "$dir") > all.offset.plotdata;  done





#
for i in *.offset; do \
distSP1=$( awk '{print ($12+($13-$12))-($2+($3-$2))}' ${i} )
dirSP1=$(awk '{print $4"_"$14}' ${i})
distSP2=$( awk '{print ($17+($18-$17))-($7+($8-$7))}' ${i} )
dirSP2=$( awk '{print $9"_"$19 }' ${i} )
dir=$( paste <(echo "$dirSP1") <(echo "$dirSP2") | while read -r dirSP1 dirSP2 ; do \
if [ $dirSP1 == $dirSP2 ] ; then echo YES_same ; \
elif [ "$dirSP1" == '-_+' ] && [ "$dirSP2" == '+_-' ]; then echo NO_large_rearrangement ; \
elif [ "$dirSP1" == '+_-' ] && [ "$dirSP2" == '-_+' ]; then echo NO_large_rearrangement ; \
elif [ "$dirSP1" == '-_-' ] && [ "$dirSP2" == '+_+' ]; then echo NO_large_rearrangement ; \
elif [ "$dirSP1" == '+_+' ] && [ "$dirSP2" == '-_-' ]; then echo NO_large_rearrangement ; \
else echo NO_recombined ; fi ;  done )
paste <(echo "$distSP1") <(echo "$dirSP1") <(echo "$distSP2") <(echo "$dirSP2")  <(echo "$dir") > "${i%.offset}".plotdata_V2;  done
```


### make plot - Supplementary Figure 3b <a name="figureS3b"></a>
```R
# plot pairwise distance between orthologs comparing position on Ce with position on Hc chromosomes

# # load libraries
library(ggplot2)
library(patchwork)
library(ggExtra)

# included this to as example of loading individual chromosomes. Did not use this in the paper, but explored the data
#data_1<-read.table("hc_ce_1to1.coords.chr1.sorted.plotdata_V2", header=F)
#data_1<-read.table("hc_ce_1to1.coords.chr2.sorted.plotdata_V2", header=F)
#data_1<-read.table("hc_ce_1to1.coords.chr3.sorted.plotdata_V2", header=F)
#data_1<-read.table("hc_ce_1to1.coords.chr4.sorted.plotdata_V2", header=F)
#data_1<-read.table("hc_ce_1to1.coords.chr5.sorted.plotdata_V2", header=F)

# load data
data <-read.table("all.plotdata_V2", header=F)

# correlations
small <- data[abs(data$V3) < 100000 & abs(data$V1) < 100000,]
large <- data[abs(data$V3) > 100000 & abs(data$V1) > 100000,]

cor.test((abs(small$V1)),(abs(small$V3)),method=c("spearman"))
> rho = 0.45999
> p-value < 2.2e-16

cor.test((abs(large$V1)),(abs(large$V3)),method=c("spearman"))
> rho = 0.0168747
> p-value = 0.6595


small_lm<-lm(abs(small$V3)~abs(small$V1))

#plot
plot <- ggplot(data,aes((abs(V1)),(abs(V3)),colour = factor(V5)))+
			geom_point(alpha=0.5)+
			scale_x_log10(limits=c(1E2,5E7))+scale_y_log10(limits=c(1E2,5E7))+
			theme_bw()+theme(legend.position = "bottom")+
			labs(x="Distance between H. contortus orthologs (Log10[bp])", y="Distance between C. elegans orthologs (Log10[bp])")

plot2 <- ggMarginal(plot,groupFill = TRUE, groupColour = TRUE,type = "histogram",binwidth = 0.05)
plot2

ggsave(plot=plot2,"Hc_Ce_pairwise_ortholog_distance.pdf",useDingbats=FALSE)

```




### make plot - Supplementary Figure 3c <a name="figureS3c"></a>
```R
# plot proportion of gene arrangements for distances less than 100kbp vs all

# load libraries
library(reshape2)
library(ggplot2)

# proportions determined - just reformat for plotting
hc_100 <- c(0.434196891,0.428497409,0.137305699)
hc_all <- c(0.306791955,0.311407847,0.381800198)
data<-rbind(hc_100,hc_all)
colnames(data)<-c("Conserved","Reorientated","Recombined")


data<- melt(data)

col<-c("#2600FF","#FF0000","#00FF00")

#plot
ggplot(data,aes(x=Var1,y=value,fill=Var2))+
	geom_bar(stat="identity",position=position_dodge())+
	labs(x="Pairwise ortholog order between H. contortus and C. elegans",y="Proportion of total")+
	scale_fill_manual(values=col)+
	ylim(0,0.5)+
	theme(legend.position = "")+
	theme_bw()
ggsave("hc_ce_ortholog_summary_proportions.pdf")
```



### Check of runs of colinear genes

```bash
cat all.offset | cut -f1,2,6,7 > synteny_count.data

for i in {2..10}; do
perl /nfs/users/nfs_j/jc17/bin/SyntenyChecker.pl --A_chr_col=0 --A_pos_col=1 --B_chr_col=2 --B_pos_col=3 --window_size=${i} --window_skip=1 --synteny_table=${i}_genes_synteny_table --detailed_table=${i}_genes_detailed_table $PWD/synteny_count.data;
done
```

### make plot - Supplementary Figure 3d <a name="figureS3d"></a>
```R
# load libraries
library(reshape2)
library(ggplot2)
library(ggrepel)

# manually made the table, using the data from the synteny counter output.
colinear <- c(3040,1405,652,303,131,53,22,9,1)
mixed <- c(3020,4649,5396,5739,5905,5977,6002,6009,6011)
gene_window <-c(2,3,4,5,6,7,8,9,10)
data<-data.frame(cbind(gene_window,colinear,mixed))

# plot
ggplot(data,aes(gene_window,((colinear/(colinear+mixed))*100),label=colinear))+
	geom_line()+
	geom_point()+
	ylim(0,100)+
	theme_bw()+
	labs(x="Number of colinear genes", y="Conserved gene order (percent)")+
	geom_text_repel()

ggsave("hc_ce_ortholog_colinear_counts.pdf")

```

### determine distribution of pairwise distances across Hc and Ce genomes
```bash
# calculate difference between Hc only orthologs
for i in hc_ce_1to1.coords.chr*.sorted; do paste <( cut -f1,2,3 ${i} | sort -k1,1 -k2,2n | head -n -1 ) <( cut -f1,2,3 ${i}  | sort -k1,1 -k2,2n | tail -n +2 ) >${i}_hconly; done
cat *hconly > hconly_increment_pairs


# calculate difference between Ce only orthologs
for i in hc_ce_1to1.coords.chr*.sorted; do paste <( cut -f6,7,8 ${i}  | sort -k1,1 -k2,2n | head -n -1) <( cut -f6,7,8 ${i}  | sort -k1,1 -k2,2n | tail -n +2 ) >${i}_ceonly; done
cat *ceonly > ceonly_increment_pairs
```

### make plot
```R
# load libraries
library(ggplot2)
library(patchwork)

# load data
ce <- read.table("ceonly_increment_pairs",header=F)
ce$species <- "ce"

hc <- read.table("hconly_increment_pairs",header=F)
hc$species <- "hc"

data <- rbind(hc,ce)

# plot distributions on chromosomes by species
a<-ggplot(data,aes((V5-V2),col=species))+geom_density()
b<-ggplot(data,aes(log10(V5-V2),col=species))+geom_density()
a+b+plot_layout(ncol=2)


# distance between orthologs along genome
ggplot(data,aes(V2,log10(V5-V2),col=species))+geom_point()+facet_grid(V1~.)

# didnt use this in the paper in the end
```


### Comparing the length of runs of syntenic orthologs between C. elegans and H. contortus  
```shell
working dir: /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/ORTHOLOGY/SYNTENY

# run James' synteny checker script
./run_jc_syntenychecker.sh rerun

for i in `ls rerun* | sort -V `; do grep --with-filename "SC" ${i}; done | sed -e 's/:/\t/g' -e 's/rerun_//g' -e 's/_genes_detailed_table//g' > colinear_genes.data
```

### Make plot - Supplementary Figure 3e <a name="figureS3e"></a>
```R
# load libraries
library(ggplot2)

# import data
data <- read.table("colinear_genes.data",header=F)

# extract colinear genesets with 5 or more genes
data2 <- data[data$V1>=5,]

# fix labels for facet grid
chr.labels <- c("1","2","3","4","5", "X")
names(chr.labels) <- c("hcontortus_chr1_Celeg_TT_arrow_pilon","hcontortus_chr2_Celeg_TT_arrow_pilon","hcontortus_chr3_Celeg_TT_arrow_pilon","hcontortus_chr4_Celeg_TT_arrow_pilon","hcontortus_chr5_Celeg_TT_arrow_pilon","hcontortus_chrX_Celeg_TT_arrow_pilon")

# make plot
ggplot(data2)+
     geom_rect(aes(xmin=V3/10E5,ymin=0,xmax=V4/10E5,ymax=1,fill=factor(V1)))+
     facet_grid(V2~.,switch="y",labeller = labeller(V2 = chr.labels))+
     scale_fill_brewer(palette="Reds")+
     labs(x="Genomic position (Mb)",fill="Colinear genes (n)") +
     theme_classic()+
     theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

ggsave("colinear_genes_per_chromosome.pdf")
```


[↥ **Back to top**](#top)

****
## BUSCO and CEGMA <a name="cegmabusco"></a>

working dir: /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/GENOME_QC

I have some scripts I use to run BUSCO and CEGMA shown below. Both tools were run on
- V4 genome
- V4 haplotype genome
- V1 genome
- McMaster genome
- NZ genome


run_busco_nematode.sh
```bash

#!/usr/local/bin/bash

# run nematode busco

export AUGUSTUS_CONFIG_PATH=/nfs/users/nfs_s/sd21/software/augustus-3.2.1/config

PREFIX=$1
REF=$2


/nfs/users/nfs_s/sd21/lustre118_link/software/ASSEMBLY_QC/busco_v3/scripts/run_BUSCO.py \
--in ${REF} \
--out ${PREFIX} \
--mode genome \
--lineage_path /nfs/users/nfs_s/sd21/databases/busco/nematoda_odb9 \
--species caenorhabditis \
--cpu 8 --force --restart --long --blast_single_core \
--tmp_path ${REF}.tmp
```

run_cegma.sh
```bash
#!/usr/local/bin/bash

# run cegma

export CEGMA=/nfs/users/nfs_s/sd21/lustre118_link/software/ASSEMBLY_QC/CEGMA_v2.5
export CEGMATMP=$PWD/tmp
export PERL5LIB=$PERL5LIB:$CEGMA/lib


PREFIX=$1
REF=$2


# Step 1: fix fasta

#fastaq enumerate_names --suffix ${PREFIX} ${REF} REF.fa.tmp
#fastaq acgtn_only REF.fa.tmp REF.fa
#rm REF.fa.tmp
#cp $REF REF.fa
awk '/^>/{print ">sequence" ++i; next}{print}' < ${REF} > REF.fa

# Step 2: run cegma
/nfs/users/nfs_s/sd21/lustre118_link/software/ASSEMBLY_QC/CEGMA_v2.5/bin/cegma --genome REF.fa --output ${PREFIX} --verbose --threads 1

```




## make a heatmap of BUSCO data to compare missingness across clade 5 nematodes
- this is used in Supplementary Figure 4 of the Supplementary INformation
- the data comes from running BUSCO on all of the nematodes in WBP

```bash
working dir: /nfs/users/nfs_s/sd21/lustre118_link/WBP/WBP_GENOMES/BUSCO_ANALYSIS
mkdir CLADE_V

#get the full output from busco
while read NAME; do cp ../BUSCO_COMPLETE/run_${NAME}*/full* CLADE_V/ ; done < clade5.list

while read NAME; do cut -f1,2 CLADE_V/full_table_${NAME}* | grep -v "#" > CLADE_V/${NAME}.busco.data; done < clade5.list

for i in *data; do mv ${i} ${i%.data}.txt; done

for i in *.txt; do sed -i -e 's/Complete/1/g' -e 's/Duplicated/2/g' -e 's/Fragmented/0.5/g' -e 's/Missing/0/g' ${i}; done
```

### Make some plots <a name="figureS4"></a>
```R
# read libraries
library(ggplot2)
library(dplyr)
library(gplots)
library(tidyverse)


filenames <- gsub("\\.txt$","", list.files(pattern="\\.txt$"))

for(i in filenames){
  assign(paste("busco_",i,sep=""), read.delim(paste(i, ".txt", sep=""), header=F))
}


inner_join(c(busco_ancylostoma_caninum.PRJNA72585.busco, busco_ancylostoma_ceylanicum.PRJNA231479.busco,busco_ancylostoma_ceylanicum.PRJNA72583.busco), by="V1") %>% distinct()

join_all(c(busco_ancylostoma_caninum.PRJNA72585.busco, busco_ancylostoma_ceylanicum.PRJNA231479.busco,busco_ancylostoma_ceylanicum.PRJNA72583.busco), by='V1', type='left')


busco_data <- left_join(busco_ancylostoma_caninum.PRJNA72585.busco, busco_ancylostoma_ceylanicum.PRJNA231479.busco, by='V1') %>%
                left_join(.,  busco_ancylostoma_ceylanicum.PRJNA72583.busco, by='V1') %>%
                left_join(.,  busco_ancylostoma_duodenale.PRJNA72581.busco, by='V1') %>%
                left_join(.,  busco_angiostrongylus_cantonensis.PRJEB493.busco, by='V1') %>%
                left_join(.,  busco_angiostrongylus_costaricensis.PRJEB494.busco, by='V1') %>%
                left_join(.,  busco_caenorhabditis_angaria.PRJNA51225.busco, by='V1') %>%
                left_join(.,  busco_caenorhabditis_brenneri.PRJNA20035.busco, by='V1') %>%
                left_join(.,  busco_caenorhabditis_briggsae.PRJNA10731.busco, by='V1') %>%
                left_join(.,  busco_caenorhabditis_elegans.PRJNA13758.busco, by='V1') %>%
                left_join(.,  busco_caenorhabditis_japonica.PRJNA12591.busco, by='V1') %>%
                left_join(.,  busco_caenorhabditis_latens.PRJNA248912.busco, by='V1') %>%
                left_join(.,  busco_caenorhabditis_nigoni.PRJNA384657.busco, by='V1') %>%
                left_join(.,  busco_caenorhabditis_remanei.PRJNA248909.busco, by='V1') %>%
                left_join(.,  busco_caenorhabditis_remanei.PRJNA248911.busco, by='V1') %>%
                left_join(.,  busco_caenorhabditis_remanei.PRJNA53967.busco, by='V1') %>%
                left_join(.,  busco_caenorhabditis_sinica.PRJNA194557.busco, by='V1') %>%
                left_join(.,  busco_caenorhabditis_sp34.PRJDB5687.busco, by='V1') %>%
                left_join(.,  busco_caenorhabditis_tropicalis.PRJNA53597.busco, by='V1') %>%
                left_join(.,  busco_cylicostephanus_goldi.PRJEB498.busco, by='V1') %>%
                left_join(.,  busco_dictyocaulus_viviparus.PRJEB5116.busco, by='V1') %>%
                left_join(.,  busco_dictyocaulus_viviparus.PRJNA72587.busco, by='V1') %>%
                left_join(.,  busco_diploscapter_coronatus.PRJDB3143.busco, by='V1') %>%
                left_join(.,  busco_diploscapter_pachys.PRJNA280107.busco, by='V1') %>%
                left_join(.,  busco_haemonchus_contortus.PRJEB506.busco, by='V1') %>%
                left_join(.,  busco_haemonchus_contortus.PRJNA205202.busco, by='V1') %>%
                left_join(.,  busco_haemonchus_placei.PRJEB509.busco, by='V1') %>%
                left_join(.,  busco_heligmosomoides_polygyrus.PRJEB1203.busco, by='V1') %>%
                left_join(.,  busco_heligmosomoides_polygyrus.PRJEB15396.busco, by='V1') %>%
                left_join(.,  busco_heterorhabditis_bacteriophora.PRJNA13977.busco, by='V1') %>%
                left_join(.,  busco_necator_americanus.PRJNA72135.busco, by='V1') %>%
                left_join(.,  busco_nippostrongylus_brasiliensis.PRJEB511.busco, by='V1') %>%
                left_join(.,  busco_oesophagostomum_dentatum.PRJNA72579.busco, by='V1') %>%
                left_join(.,  busco_oschieus_tipulae.PRJEB15512.busco, by='V1') %>%
                left_join(.,  busco_pristionchus_exspectatus.PRJEB6009.busco, by='V1') %>%
                left_join(.,  busco_pristionchus_pacificus.PRJNA12644.busco, by='V1') %>%
                left_join(.,  busco_strongylus_vulgaris.PRJEB531.busco, by='V1') %>%
                left_join(.,  busco_teladorsagia_circumcincta.PRJNA72569.busco, by='V1') %>% distinct()


colnames(busco_data) <- c("BUSCO",
"ancylostoma_caninum.PRJNA72585",
"ancylostoma_ceylanicum.PRJNA231479",
"ancylostoma_ceylanicum.PRJNA72583",
"ancylostoma_duodenale.PRJNA72581",
"angiostrongylus_cantonensis.PRJEB493",
"angiostrongylus_costaricensis.PRJEB494",
"caenorhabditis_angaria.PRJNA51225",
"caenorhabditis_brenneri.PRJNA20035",
"caenorhabditis_briggsae.PRJNA10731",
"caenorhabditis_elegans.PRJNA13758",
"caenorhabditis_japonica.PRJNA12591",
"caenorhabditis_latens.PRJNA248912",
"caenorhabditis_nigoni.PRJNA384657",
"caenorhabditis_remanei.PRJNA248909",
"caenorhabditis_remanei.PRJNA248911",
"caenorhabditis_remanei.PRJNA53967",
"caenorhabditis_sinica.PRJNA194557",
"caenorhabditis_sp34.PRJDB5687",
"caenorhabditis_tropicalis.PRJNA53597",
"cylicostephanus_goldi.PRJEB498",
"dictyocaulus_viviparus.PRJEB5116",
"dictyocaulus_viviparus.PRJNA72587",
"diploscapter_coronatus.PRJDB3143",
"diploscapter_pachys.PRJNA280107",
"haemonchus_contortus.PRJEB506",
"haemonchus_contortus.PRJNA205202",
"haemonchus_placei.PRJEB509",
"heligmosomoides_polygyrus.PRJEB1203",
"heligmosomoides_polygyrus.PRJEB15396",
"heterorhabditis_bacteriophora.PRJNA13977",
"necator_americanus.PRJNA72135",
"nippostrongylus_brasiliensis.PRJEB511",
"oesophagostomum_dentatum.PRJNA72579",
"oschieus_tipulae.PRJEB15512",
"pristionchus_exspectatus.PRJEB6009",
"pristionchus_pacificus.PRJNA12644",
"strongylus_vulgaris.PRJEB531",
"teladorsagia_circumcincta.PRJNA72569")

busco_data <- busco_data %>% remove_rownames %>% column_to_rownames(var="BUSCO")




heatmap.2(as.matrix(busco_data), scale = "none",
          trace = "none", density.info = "none")
```

[↥ **Back to top**](#top)


******
## Repeats <a name="repeats"></a>
- repeat masker was run on:
     - Haem V4
     - Haem V4 HAPLOTYPES
     - Haem V1
- below is example command used.
```bash
# RepeatMasker
/nfs/users/nfs_s/sd21/lustre118_link/software/REPEATMASKER/RepeatModeler-open-1.0.11/BuildDatabase -name HAEM_V4 -engine ncbi HAEM_V4_final.chr.fa

/nfs/users/nfs_s/sd21/lustre118_link/software/REPEATMASKER/RepeatModeler-open-1.0.11/RepeatModeler -pa 20 -engine ncbi -database HAEM_V4

/nfs/users/nfs_s/sd21/lustre118_link/software/REPEATMASKER/RepeatMasker/RepeatMasker -e ncbi -pa 7 -s -dir ./RM_V4_OUT -small -gff -lib RM_43484.TueMar271017292018/consensi.fa.classified HAEM_V4_final.chr.fa


# LTRfinder

# Pipeline to run LTRharvest and LTRdigest
# largely based on Avrils notes on Avrilomics: http://avrilomics.blogspot.com/2015/09/ltrharvest.html


# Step 1 : create a suffix array index for the genome assembly

bsub.py 5 01_ltr_suffixerator gt suffixerator -db HAEM_V4_final.chr.fa -indexname HAEM_V4_final.chr -tis -suf -lcp -des -ssp -sds -dna



# Step 2: run LTRharvest and sort the output

bsub.py 5 02_ltr_harvest gt ltrharvest -index HAEM_V4_final.chr -gff3 HAEM_V4_final.chr.ltr.gff -out HAEM_V4_final.chr.ltr.fa -seqids yes

bsub.py 1 03_sort_ltrs "gt gff3 -sort  HAEM_V4_final.chr.ltr.gff \> HAEM_V4_final.chr.ltr.sorted.gff"


# Step 3: run LTRdigest and filter output

bsub.py 5 04_ltrdigest gt ltrdigest -hmms hmm/*hmm -aaout -outfileprefix HAEM_V4_final.chr -seqfile HAEM_V4_final.chr.fa -matchdescstart \< HAEM_V4_final.chr.ltr.sorted.gff \> HAEM_V4_final
.chr.ltrdigest.gff

#--- remove all LTR retrotransposon candidates that don't have any domain hit at all (to help get rid of things that might not be LTR retrotransposon insertions)
gt select -rule_files ~sd21/bin/filter_protein_match.lua -- <  HAEM_V4_final.chr.ltrdigest.gff > HAEM_V4_final.chr.ltrdigest.gff.2

#--- Extra the DNA sequences for the remaining full-length elements
gt extractfeat -type LTR_retrotransposon -matchdescstart -retainids -encseq HAEM_V4_final.chr.fa HAEM_V4_final.chr.ltrdigest.gff.2 > HAEM_V4_final.chr.ltrdigest.gff.2.fa
```






### Make some plots of genome wide repeat distribution - Supplementary Figure 5 <a name="figureS5"></a>
```
cd /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/REPEATS/CHR/RM_V4_OUT/FIGURE

# get data
ln -s ../HAEM_V4_final.chr.fa.out

ln -s ../../../../REF/HAEM_V4_final.chr.fa

# get repeat types
awk '{print $11}'  HAEM_V4_final.chr.fa.out | sort | uniq -c

    138 DNA
   3055 DNA/CMC-EnSpm
    412 DNA/hAT-Ac
    112 DNA/hAT-Tag1
     59 DNA/Kolobok-T2
    169 DNA/Merlin
    195 DNA/MuLE-MuDR
    640 DNA/MULE-MuDR
    169 DNA/PIF-Harbinger
   1629 DNA/PiggyBac
    218 DNA/TcMar-ISRm11
  16968 DNA/TcMar-Mariner
  11527 DNA/TcMar-Tc1
   2800 DNA/TcMar-Tc4
   5354 LINE/CR1
    343 LINE/CR1-Zenon
     18 LINE/Dong-R4
   2537 LINE/L1
     49 LINE/L1-Tx1
   2264 LINE/L2
  10652 LINE/Penelope
      1 LINE/R1-LOA
     27 LINE/R2
   7270 LINE/RTE-BovB
  16844 LINE/RTE-RTE
   4575 Low_complexity
    217 LTR
    330 LTR/Copia
    109 LTR/DIRS
   1231 LTR/ERV1
    125 LTR/ERV4
   1915 LTR/ERVK
   5313 LTR/Gypsy
     64 LTR/Ngaro
   5178 LTR/Pao
      1 position
   3411 RC/Helitron
    118 rRNA
    359 Satellite
  30629 Simple_repeat
   1523 SINE?
  10474 SINE/tRNA
 205593 Unknown



# get coordinates of repeat types, output as a bed file

awk '{print $11}'  HAEM_V4_final.chr.fa.out | sort | uniq | while read RPT; do OUT="$( echo $RPT | sed 's/\//_/g' )" ; grep  "$RPT" HAEM_V4_final.chr.fa.out | awk '{print $5,$6,$7}' OFS="\t" > $OUT.repeats.bed ; done


# make some windows in the genome - 100kb
samtools faidx  HAEM_V4_final.chr.fa
cut -f1,2 HAEM_V4_final.chr.fa.fai > HAEM_V4.genome
bedtools-2 makewindows -g HAEM_V4.genome -w 200000 > HAEM_V4.200kb_windows.bed


# calculate coverage of repeats
for i in *repeats.bed; do
	NAME="$( echo $i | sed 's/.repeats.bed//g' )" ;
	COUNT="$( cat ${i} | wc -l )";
	bedtools-2 coverage -a HAEM_V4.200kb_windows.bed -b ${i} | awk -v NAME=$NAME -v COUNT=$COUNT '{print $0,NAME" n="COUNT}' OFS="\t";
	done > genome_repeat_coverage.200k.data
```

### make some plots
```R
library(ggplot2)

data <- read.table("genome_repeat_coverage.200k.data",header=F,sep="\t")
#--- remove mtDNA
data <- data[data$V1!="hcontortus_chr_mtDNA_arrow_pilon",]

chromosome.labels <- c("I","II","III","IV","V", "X" )
names(chromosome.labels) <- c("hcontortus_chr1_Celeg_TT_arrow_pilon",
	"hcontortus_chr2_Celeg_TT_arrow_pilon",
	"hcontortus_chr3_Celeg_TT_arrow_pilon",
	"hcontortus_chr4_Celeg_TT_arrow_pilon",
	"hcontortus_chr5_Celeg_TT_arrow_pilon",
	"hcontortus_chrX_Celeg_TT_arrow_pilon")




ggplot(data,aes(V2,V7,fill=V8))+geom_area()+facet_grid(V1~.)

#--- remove the "Unknown" class, as they are abundant, but dont add much

data2<-data[data$V8!="Unknown n=205593",]
ggplot(data2,aes(V2,V7,fill=V8))+
	geom_area()+
	facet_grid(V1~.,labeller = labeller(V1 = chromosome.labels))+
	theme_bw()+
	labs(x="Genomic position (bp)",
       	y="Fraction of 200 kb window containing repeat type",
       	fill="Repeat type")+
	guides(fill=guide_legend(ncol=1))+
	ylim(0,1)
```
[↥ **Back to top**](#top)
******
