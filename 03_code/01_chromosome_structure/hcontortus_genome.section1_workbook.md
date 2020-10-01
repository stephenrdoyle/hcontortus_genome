# Haemonchus genome paper
## Section 1: Chromosome structure of Haemonchus contortus

1. [Genome assembly](#genome)
2. [Completion of X chromosome](#xchromosome)
3. [Genome polishing](#polishing)
4. [Circos plot - Figure 1A](#circos)
5. [Orthology conservation and order per chromosome - Figure 1B top and bottom](#orthology)
6. [Microsynteny](#microsynteny)
7. [Genome stats](#genomestats)
8. [Genome completeness - CEGMA & BUSCO](#cegmabusco)
9. [Comparative analysis of the NZ Haemonchus genome](#nzgenome)



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
/nfs/users/nfs_s/sd21/lustre118_link/software/10X_GENOMICS/longranger-2.1.2/longranger basic --id=hc_inbred --fastqs=/nfs/users/nfs_s/sd21/lustre118_link/hc/10X_genomics/chromium/hc_inbred --readgroup=hc_inbred --localcores=4

# map reads
bwa mem \
-t 30 xchr_TaTb_plus_bits.sized.renamed.fa \
-p CHROMIUM_interleaved.fastq |\
samtools-1.3 view -Sb - |\
samtools-1.3 sort -n -o ./CHROMIUM-sorted.bam -


# make a bam list from the mapped reads
ls -1 CHROMIUM-sorted.bam > bam.list

# run ARCS
~sd21/lustre118_link/software/10X_GENOMICS/arcs/Arcs/arcs -f xchr_TaTb_plus_bits.sized.renamed.fa -a bam.list -s 95 -c 3 -l 0 -z 500 -m 5-10000 -d 0 -e 50000 -r 0.05 -i 16 -v 1


# run LINKS
LINKS -f xchr_TaTb_plus_bits.sized.renamed.fa -s empty.fof -k 15 -b xchr_TaTb_plus_bits.sized.renamed.fa.scaff_s95_c3_l0_d0_e50000_r0.05_original -l 1 -t 2 -a 0.9 -v 1
```



[↥ **Back to top**](#top)
******
## Polishing <a name="polishing"></a>
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

java -Xmx200G -jar /nfs/users/nfs_s/sd21/lustre118_link/software/GENOME_IMPROVEMENT/pilon-1.22/pilon-1.22.jar --genome HAEM_V4_arrow.fa --bam Pilon_prep_singlefemale2HaemV4.merged.sorted.marked.bam --fix all --verbose --diploid --changes --chunksize 1000000 --nostrays --threads 30
```



 [↥ **Back to top**](#top)
 ******
## Circos plot - Figure 1A <a name="circos"></a>
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
```



[↥ **Back to top**](#top)
******
## Microsynteny <a name="microsynteny"></a>

### Comparing runs of syntenic orthologs between C. elegans and H. contortus  
```shell
working dir: /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/ORTHOLOGY/SYNTENY

# run James' synteny checker script
./run_jc_syntenychecker.sh rerun

for i in `ls rerun* | sort -V `; do grep --with-filename "SC" ${i}; done | sed -e 's/:/\t/g' -e 's/rerun_//g' -e 's/_genes_detailed_table//g' > colinear_genes.data
```

Load R
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



#

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



# Make the plot - Figure 1 b
```R
# load
R-3.5.0
library(ggplot2)
library(patchwork)
library(dplyr)

# plot number of genes per chromosome
data<-read.table("hc_ce_orthologcount_perchr.txt",header=F)
data.1<- data %>% group_by(V1) %>%  mutate(per=paste0(round(V3/sum(V3)*100, 2), "%")) %>% ungroup


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





#––––––– PLAY
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

# plot pairwise distance between orthologs comparing position on Ce with position on Hc chromosomes
```R
library(ggplot2)
library(patchwork)
library(ggExtra)


#data_1<-read.table("hc_ce_1to1.coords.chr1.sorted.plotdata_V2", header=F)
#data_1<-read.table("hc_ce_1to1.coords.chr2.sorted.plotdata_V2", header=F)
#data_1<-read.table("hc_ce_1to1.coords.chr3.sorted.plotdata_V2", header=F)
#data_1<-read.table("hc_ce_1to1.coords.chr4.sorted.plotdata_V2", header=F)
#data_1<-read.table("hc_ce_1to1.coords.chr5.sorted.plotdata_V2", header=F)
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


plot<-ggplot(data,aes((abs(V1)),(abs(V3)),colour = factor(V5)))+
			geom_point(alpha=0.5)+
			scale_x_log10(limits=c(1E2,5E7))+scale_y_log10(limits=c(1E2,5E7))+
			theme_bw()+theme(legend.position = "bottom")+
			labs(x="Distance between H. contortus orthologs (Log10[bp])", y="Distance between C. elegans orthologs (Log10[bp])")

plot2<-ggMarginal(plot,groupFill = TRUE, groupColour = TRUE,type = "histogram",binwidth = 0.05)
plot2
ggsave(plot=plot2,"Hc_Ce_pairwise_ortholog_distance.pdf",useDingbats=FALSE)

# used as Supplementary Figure 3b
```





### plot proportion of gene arrangements for distances less than 100kbp vs all
```R
hc_100 <- c(0.434196891,0.428497409,0.137305699)
hc_all<-c(0.306791955,0.311407847,0.381800198)
data<-rbind(hc_100,hc_all)
colnames(data)<-c("Conserved","Reorientated","Recombined")

library(reshape2)
data<- melt(data)

col<-c("#2600FF","#FF0000","#00FF00")

ggplot(data,aes(x=Var1,y=value,fill=Var2))+
	geom_bar(stat="identity",position=position_dodge())+
	labs(x="Pairwise ortholog order between H. contortus and C. elegans",y="Proportion of total")+
	scale_fill_manual(values=col)+
	ylim(0,0.5)+
	theme(legend.position = "")+
	theme_bw()
ggsave("hc_ce_ortholog_summary_proportions.pdf")

# used as Supplementary Figure 3c
```



### Check of runs of colinear genes

```bash
cat all.offset | cut -f1,2,6,7 > synteny_count.data

for i in {2..10}; do
perl /nfs/users/nfs_j/jc17/bin/SyntenyChecker.pl --A_chr_col=0 --A_pos_col=1 --B_chr_col=2 --B_pos_col=3 --window_size=${i} --window_skip=1 --synteny_table=${i}_genes_synteny_table --detailed_table=${i}_genes_detailed_table $PWD/synteny_count.data;
done
```

### in R
```R
library(reshape2)
library(ggplot2)
library(ggrepel)

colinear <- c(3040,1405,652,303,131,53,22,9,1)
mixed <- c(3020,4649,5396,5739,5905,5977,6002,6009,6011)
gene_window <-c(2,3,4,5,6,7,8,9,10)
data<-data.frame(cbind(gene_window,colinear,mixed))

ggplot(data,aes(gene_window,((colinear/(colinear+mixed))*100),label=colinear))+
	geom_line()+
	geom_point()+
	ylim(0,100)+
	theme_bw()+
	labs(x="Number of colinear genes", y="Conserved gene order (percent)")+
	geom_text_repel()

ggsave("hc_ce_ortholog_colinear_counts.pdf")

# used as Supplementary Figure 3d
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

#in R
```R
library(ggplot2)
library(patchwork)

ce <- read.table("ceonly_increment_pairs",header=F)
ce$species <- "ce"

hc <- read.table("hconly_increment_pairs",header=F)
hc$species <- "hc"

data <- rbind(hc,ce)

# distribution of ortholog distances
a<-ggplot(data,aes((V5-V2),col=species))+geom_density()
b<-ggplot(data,aes(log10(V5-V2),col=species))+geom_density()
a+b+plot_layout(ncol=2)


# distance between orthologs along genome
ggplot(data,aes(V2,log10(V5-V2),col=species))+geom_point()+facet_grid(V1~.)

# didnt use this in the paper in the end
```
