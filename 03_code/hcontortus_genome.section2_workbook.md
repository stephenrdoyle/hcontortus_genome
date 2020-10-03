# Haemonchus genome paper

## Section 2: Resolving haplotypic diversity and repeat distribution within the chromosomes

1.  [Genome graph of chromosomes and haplotypes](#genomegraph)
    -   [Figure 2a](#figure2a)
    -   [Figure 2b](#figure2b)
2.  [Haplotype density and distribution](#haplotypes)
    -   [Figure 2c](#figure2c)
3.  [Haplotype switching between individuals](#haploswitching)
    -   [Figure 2d](#figure2d)
4.  [Repeat analyses](#repeats)
    -   [Supplementary Figure 5](#figureS5)

## Genome graph of chromosomes and haplotypes <a name="genomegraph"></a>

-   Note: some non commandline tools were used to create the figures
    -   Figure 2a was visualised using bandage <a name="figureS2a"></a>
    -   Figure 2b was visualised using bandage and GenomeRibbon <a name="figureS2b"></a>

```shell
working dir: /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/GRAPH/MINIGRAPH

# get genomes
ln -s ../../REF/HAEM_V4_final.chr.fa
ln -s ../../REF/HAEM_V4_final.haplocomtam_only.fa


# run minigraph
/nfs/users/nfs_s/sd21/lustre118_link/software/GENOME_GRAPH/minigraph/minigraph -xggs -t16 HAEM_V4_final.chr.fa HAEM_V4_final.haplocomtam_only.fa > out.gfa

# The GFA can be visualised using bandage - this was used to make FIgure 2a and the top part of Figure 2b which is zoomed in

# want to show an example of haplotypic diversity. Extracted a random haplotype on chromosome X and saved it as "test_X_other.fa"

# run Nucmer
nucmer --maxmatch --coords -p other HAEM_V4_final.chr.fa test_X_other.fa
show-coords -lTH other.delta | grep "chrX" | awk '{if($1>=7300000 && $2<=7500000) print }' > nucmer_other_GR.out

# get genome annotation for visualisation in genome ribbons
awk '{if($1=="hcontortus_chrX_Celeg_TT_arrow_pilon" && $3=="exon" && $4>=7300000 && $5<7400000) print $1, $4, $5, $9, ".", $7, $3}' OFS="\t" haemonchus_contortus.PRJEB506.WBPS13.annotations.gff3 > chrX.exons.bed

# used "nucmer_other_GR.out" and "chrX.exons.bed" as input to genome ribbon - this was used to make the bottom part of Figure 2b
```

## 2. Haplotype density and distribution <a name="haplotypes"></a>

### make plot <a name="figure2c"></a>

```R
load(file = "hcontortus_genome.workbook.Rdata")

setwd("02_data/")

# load libraries
library(circlize)
library(cnv)

# load data
chr  <-  read.table("HAEM_V4_final.chr.fa.genome", header=F)
gc_100k  <-  read.table("HAEM_V4_final.chr.fa.100k.nuc", skip = 1, header=F)
gc_100k  <-  gc_100k[(gc_100k$V5>0.41 & gc_100k$V5<0.5), ]
line_1M <- read.table("LINE_repeats.1M.counts", header=F)
sine_1M <- read.table("SINE_repeats.1M.counts", header=F)
dna_1M <- read.table("DNA_repeats.1M.counts", header=F)
ltr_1M <- read.table("LTR_repeats.1M.counts", header=F)
haplotypes <- read.table("haplotype_hits_5k.coords2", header=T)
haplotypes <- haplotypes[haplotypes$similarity>0.2, ]
haplotype_cov <- read.table("haplotype_hits_5k.1M.coverage2", header=F)

# set colours
chr_colours <- c("#b2182b", "#fc8d59", "#fee090", "#d1e5f0", "#67a9cf", "#4575b4")

#pdf("figure1_circos.pdf", useDingbats = FALSE, height = 7, width = 7)
png("../04_analysis/figure1_circos.png")
circos.par("track.height" = 0.15, start.degree = 90 , gap.degree=c(1, 1, 1, 1, 1, 12))
circos.initialize(factors = chr$V1,  x = chr$V2)


# Track 1 - GC_CONTENT
gc_mean <- mean(gc_100k$V5)
circos.track(ylim = c(0.41,  0.46), factors = gc_100k$V1,  x=gc_100k$V2, y=gc_100k$V5,
             panel.fun = function(x,  y) {
             circos.axis(labels.cex = 0.5, major.at = c(1, 1e7, 2e7, 3e7, 4e7, 5e7))
             cell.xlim = get.cell.meta.data("cell.xlim")
             circos.lines(cell.xlim,  c(gc_mean, gc_mean),  col="black", lwd=1, lty="dashed")
               })
circos.trackPoints(factors = gc_100k$V1,  x=gc_100k$V3,  y=gc_100k$V5, pch = 20, cex = 0.5, col=chr_colours)
circos.yaxis(labels.cex=0.4, side = 'left', tick = T, sector.index = 'hcontortus_chr1_Celeg_TT_arrow_pilon')


# Track 2 - Haplotypes - length vs identity
rbPal  <-  colorRampPalette(c('lightgrey', 'black'))
haplotypes$Col  <-  rbPal(10)[as.numeric(cut(haplotypes$similarity, breaks = 10))]

circos.track(ylim = c(3.5, 5.5), factors=haplotypes$chromosome,  haplotypes$start,  haplotypes$length)
#circos.rect(haplotypes$start, haplotypes$length, haplotypes$end, haplotypes$length, lty="solid", lwd=1)
circos.trackPoints(factors =haplotypes$chromosome,  x=haplotypes$start,  y=log10(haplotypes$length), pch = 20, cex = 0.5, col = haplotypes$Col)
circos.yaxis(labels.cex=0.4, side = 'left', tick = T, sector.index = 'hcontortus_chr1_Celeg_TT_arrow_pilon')

# Track 3 - haplotype - coverage per 1 Mb window
circos.track(ylim = c(0, 1), factors=haplotype_cov$V1,  haplotype_cov$V2,  haplotype_cov$V7)
#circos.rect(haplotypes$start, haplotypes$length, haplotypes$end, haplotypes$length, lty="solid", lwd=1)
circos.trackLines(factors = haplotype_cov$V1,  haplotype_cov$V3-500000,  haplotype_cov$V7, col=chr_colours, lwd=2)
circos.yaxis(labels.cex=0.4, side = 'left', tick = T, sector.index = 'hcontortus_chr1_Celeg_TT_arrow_pilon')



# Track 4 - Repeats - LINEs,  SINEs,  DNA,  LTRs
circos.track(factors=line_1M$V1,  line_1M$V2,  line_1M$V4)
circos.trackLines(factors = line_1M$V1, line_1M$V2,  line_1M$V4, col="green", lwd=2)
circos.trackLines(factors = sine_1M$V1, sine_1M$V2,  sine_1M$V4, col="red", lwd=2)
circos.trackLines(factors = dna_1M$V1, dna_1M$V2,  dna_1M$V4, col="blue", lwd=2)
circos.trackLines(factors = ltr_1M$V1, ltr_1M$V2,  ltr_1M$V4, col="orange", lwd=2)
circos.yaxis(labels.cex=0.4, side = 'left', tick = T, sector.index = 'hcontortus_chr1_Celeg_TT_arrow_pilon')


# finish
circos.clear()
dev.off()
```

## Haplotype switching <a name="haploswitching"></a>

working dir: /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/HAPLOTYPES/HAPLO_SWITCHING

```bash
# get data -
pathfind -t lane -i 21766_7 -l ./ --rename --filetype fastq


# make samples_lanes.list containing:
MHco3(ISE).N1_1	21766_7_1
MHco3(ISE).N1_10	21766_7_10
MHco3(ISE).N1_11	21766_7_11
MHco3(ISE).N1_2	21766_7_2
MHco3(ISE).N1_3	21766_7_3
MHco3(ISE).N1_4	21766_7_4
MHco3(ISE).N1_5	21766_7_5
MHco3(ISE).N1_6	21766_7_6
MHco3(ISE).N1_7	21766_7_7
MHco3(ISE).N1_8	21766_7_8
MHco3(ISE).N1_9	21766_7_9



# get reference, which contains chromosomes AND haplotypes
cp ../../REF/HAEM_V4_final.all.fa .


# set off mapping
echo "while read name lane; do \
~sd21/bash_scripts/run_bwamem_splitter $"{name}"_$"{lane}" $PWD/HAEM_V4_final.all.fa \
$PWD/$"{lane}"_1.fastq.gz \
$PWD/$"{lane}"_2.fastq.gz; \
done < $PWD/samples_lanes.list" > run_mapping
chmod a+x run_mapping
screen
./run_mapping &


# once mapping has finished,  generate coverage stats in 100 kbp windows

/nfs/users/nfs_s/sd21/bash_scripts/run_cov_stats 100000


##########################################################################################
# run_cov_stats
##########################################################################################

# Usage: ~sd21/bash_scripts/run_cov_stats < window size >


WINDOW=$1

for i in *.bam
do

bamtools header -in ${i} | grep "^@SQ" | awk -F'[:\t]' '{printf $3"\t"1"\t"$5"\n"}' OFS="\t" > ${i%.bam}.chr.bed
bamtools header -in ${i} | grep "^@SQ" | awk -F'[:\t]' '{printf $3"\t"$5"\n"}' OFS="\t" > ${i%.bam}.chr.genome

bedtools makewindows -g ${i%.bam}.chr.genome -w ${WINDOW} > ${i%.bam}.${WINDOW}_window.bed

samtools-1.6 bedcov -Q 20 ${i%.bam}.chr.bed ${i} | awk -F'\t' '{printf $1"\t"$2"\t"$3"\t"$4"\t"$4/($3-$2)"\n"}' OFS="\t" > ${i%.bam}.chr.cov
samtools-1.6 bedcov -Q 20 ${i%.bam}.${WINDOW}_window.bed ${i} | awk -F'\t' '{printf $1"\t"$2"\t"$3"\t"$4"\t"$4/($3-$2)"\n"}' OFS="\t" > ${i%.bam}.${WINDOW}_window.cov

rm ${i%.bam}.chr.bed ${i%.bam}.${WINDOW}_window.bed ${i%.bam}.chr.genome

done

for i in *.chr.cov; do printf "${i}\n" > ${i}.tmp | awk '{print $5}' OFS="\t" ${i} >> ${i}.tmp; done
paste *.tmp > coverage_stats.summary
rm *.tmp
```

### make the plot <a name="figure2d"></a>
```R
# load libraries
require(gtools)
library(ggplot2)
library(reshape2)
library(dplyr)

# load data
filenames  <-  gsub("\\.cov.2$", "",  list.files(pattern="\\.cov.2$"))

#--- note to self - mixed sorts will sort them numerically
filenames  <-  mixedsort(filenames)

for(i in filenames){
      assign(paste(i, sep=""),  read.delim(paste(i,  ".cov.2",  sep=""),  header=F))
}

# select the columns needed from the raw data - chromosome, position, and coverage in window
MHco3_ISE.N1_1_21766_7_1.merged.sorted.marked.100000_window.chr  <-  select(MHco3_ISE.N1_1_21766_7_1.merged.sorted.marked.100000_window.chr,  c(V1, V2, V5))
MHco3_ISE.N1_2_21766_7_2.merged.sorted.marked.100000_window.chr  <-  select(MHco3_ISE.N1_2_21766_7_2.merged.sorted.marked.100000_window.chr,  c(V1, V2, V5))
MHco3_ISE.N1_3_21766_7_3.merged.sorted.marked.100000_window.chr  <-  select(MHco3_ISE.N1_3_21766_7_3.merged.sorted.marked.100000_window.chr,  c(V1, V2, V5))
MHco3_ISE.N1_4_21766_7_4.merged.sorted.marked.100000_window.chr  <-  select(MHco3_ISE.N1_4_21766_7_4.merged.sorted.marked.100000_window.chr,  c(V1, V2, V5))
MHco3_ISE.N1_5_21766_7_5.merged.sorted.marked.100000_window.chr  <-  select(MHco3_ISE.N1_5_21766_7_5.merged.sorted.marked.100000_window.chr,  c(V1, V2, V5))
MHco3_ISE.N1_6_21766_7_6.merged.sorted.marked.100000_window.chr  <-  select(MHco3_ISE.N1_6_21766_7_6.merged.sorted.marked.100000_window.chr,  c(V1, V2, V5))
MHco3_ISE.N1_7_21766_7_7.merged.sorted.marked.100000_window.chr  <-  select(MHco3_ISE.N1_7_21766_7_7.merged.sorted.marked.100000_window.chr,  c(V1, V2, V5))
MHco3_ISE.N1_8_21766_7_8.merged.sorted.marked.100000_window.chr  <-  select(MHco3_ISE.N1_8_21766_7_8.merged.sorted.marked.100000_window.chr,  c(V1, V2, V5))
MHco3_ISE.N1_9_21766_7_9.merged.sorted.marked.100000_window.chr  <-  select(MHco3_ISE.N1_9_21766_7_9.merged.sorted.marked.100000_window.chr,  c(V1, V2, V5))
MHco3_ISE.N1_10_21766_7_10.merged.sorted.marked.100000_window.chr  <-  select(MHco3_ISE.N1_10_21766_7_10.merged.sorted.marked.100000_window.chr,  c(V1, V2, V5))
MHco3_ISE.N1_11_21766_7_11.merged.sorted.marked.100000_window.chr  <-  select(MHco3_ISE.N1_11_21766_7_11.merged.sorted.marked.100000_window.chr,  c(V1, V2, V5))

#normalise per sample to 95% quantile
MHco3_ISE.N1_1_21766_7_1.merged.sorted.marked.100000_window.chr$V5  <-  MHco3_ISE.N1_1_21766_7_1.merged.sorted.marked.100000_window.chr$V5/quantile(MHco3_ISE.N1_1_21766_7_1.merged.sorted.marked.100000_window.chr$V5, 0.95) MHco3_ISE.N1_2_21766_7_2.merged.sorted.marked.100000_window.chr$V5  <-  MHco3_ISE.N1_2_21766_7_2.merged.sorted.marked.100000_window.chr$V5/quantile(MHco3_ISE.N1_2_21766_7_2.merged.sorted.marked.100000_window.chr$V5, 0.95)
MHco3_ISE.N1_3_21766_7_3.merged.sorted.marked.100000_window.chr$V5  <-  MHco3_ISE.N1_3_21766_7_3.merged.sorted.marked.100000_window.chr$V5/quantile(MHco3_ISE.N1_3_21766_7_3.merged.sorted.marked.100000_window.chr$V5, 0.95)
MHco3_ISE.N1_4_21766_7_4.merged.sorted.marked.100000_window.chr$V5  <-  MHco3_ISE.N1_4_21766_7_4.merged.sorted.marked.100000_window.chr$V5/quantile(MHco3_ISE.N1_4_21766_7_4.merged.sorted.marked.100000_window.chr$V5, 0.95)
MHco3_ISE.N1_5_21766_7_5.merged.sorted.marked.100000_window.chr$V5  <-  MHco3_ISE.N1_5_21766_7_5.merged.sorted.marked.100000_window.chr$V5/quantile(MHco3_ISE.N1_5_21766_7_5.merged.sorted.marked.100000_window.chr$V5, 0.95)
MHco3_ISE.N1_6_21766_7_6.merged.sorted.marked.100000_window.chr$V5  <-  MHco3_ISE.N1_6_21766_7_6.merged.sorted.marked.100000_window.chr$V5/quantile(MHco3_ISE.N1_6_21766_7_6.merged.sorted.marked.100000_window.chr$V5, 0.95)
MHco3_ISE.N1_7_21766_7_7.merged.sorted.marked.100000_window.chr$V5  <-  MHco3_ISE.N1_7_21766_7_7.merged.sorted.marked.100000_window.chr$V5/quantile(MHco3_ISE.N1_7_21766_7_7.merged.sorted.marked.100000_window.chr$V5, 0.95)
MHco3_ISE.N1_8_21766_7_8.merged.sorted.marked.100000_window.chr$V5  <-  MHco3_ISE.N1_8_21766_7_8.merged.sorted.marked.100000_window.chr$V5/quantile(MHco3_ISE.N1_8_21766_7_8.merged.sorted.marked.100000_window.chr$V5, 0.95)
MHco3_ISE.N1_9_21766_7_9.merged.sorted.marked.100000_window.chr$V5  <-  MHco3_ISE.N1_9_21766_7_9.merged.sorted.marked.100000_window.chr$V5/quantile(MHco3_ISE.N1_9_21766_7_9.merged.sorted.marked.100000_window.chr$V5, 0.95)
MHco3_ISE.N1_10_21766_7_10.merged.sorted.marked.100000_window.chr$V5  <-  MHco3_ISE.N1_10_21766_7_10.merged.sorted.marked.100000_window.chr$V5/quantile(MHco3_ISE.N1_10_21766_7_10.merged.sorted.marked.100000_window.chr$V5, 0.95)
MHco3_ISE.N1_11_21766_7_11.merged.sorted.marked.100000_window.chr$V5  <-  MHco3_ISE.N1_11_21766_7_11.merged.sorted.marked.100000_window.chr$V5/quantile(MHco3_ISE.N1_11_21766_7_11.merged.sorted.marked.100000_window.chr$V5, 0.95)

# bring all of the individual datasets together, and some sensible column names
data  <-  MHco3_ISE.N1_1_21766_7_1.merged.sorted.marked.100000_window.chr %>%
      left_join(.,  MHco3_ISE.N1_2_21766_7_2.merged.sorted.marked.100000_window.chr,  by=c('V1', 'V2')) %>%
      left_join(.,  MHco3_ISE.N1_3_21766_7_3.merged.sorted.marked.100000_window.chr,  by=c('V1', 'V2')) %>%
      left_join(.,  MHco3_ISE.N1_4_21766_7_4.merged.sorted.marked.100000_window.chr,  by=c('V1', 'V2')) %>%
      left_join(.,  MHco3_ISE.N1_5_21766_7_5.merged.sorted.marked.100000_window.chr,  by=c('V1', 'V2')) %>%
      left_join(.,  MHco3_ISE.N1_6_21766_7_6.merged.sorted.marked.100000_window.chr,  by=c('V1', 'V2')) %>%
      left_join(.,  MHco3_ISE.N1_7_21766_7_7.merged.sorted.marked.100000_window.chr,  by=c('V1', 'V2')) %>%
      left_join(.,  MHco3_ISE.N1_8_21766_7_8.merged.sorted.marked.100000_window.chr,  by=c('V1', 'V2')) %>%
      left_join(.,  MHco3_ISE.N1_9_21766_7_9.merged.sorted.marked.100000_window.chr,  by=c('V1', 'V2')) %>%
      left_join(.,  MHco3_ISE.N1_10_21766_7_10.merged.sorted.marked.100000_window.chr,  by=c('V1', 'V2')) %>%
      left_join(.,  MHco3_ISE.N1_11_21766_7_11.merged.sorted.marked.100000_window.chr,  by=c('V1', 'V2')) %>%
      distinct()

colnames(data)  <-  c("CHR", "START", "MHco3_ISE.N1_1", "MHco3_ISE.N1_2", "MHco3_ISE.N1_3", "MHco3_ISE.N1_4", "MHco3_ISE.N1_5", "MHco3_ISE.N1_6", "MHco3_ISE.N1_7", "MHco3_ISE.N1_8", "MHco3_ISE.N1_9", "MHco3_ISE.N1_10", "MHco3_ISE.N1_11")

# reformat the data
data2  <-  melt(data, id.vars = c("CHR", "START"))
data2  <-  data2[data2$CHR!="hcontortus_7_chr_mtDNA_arrow_pilon", ]

# remove high coverage regions by settign any window higher than 1 to 1
data2$value  <-  ifelse(data2$value > 1,  1,  data2$value)

ggplot(data2, aes(START/10^6, variable, fill=value)) +
     geom_tile() +
     labs(x="Genomic position in chromosome (Mbp)",  y="",  fill="Relative\ncoverage") +
     scale_y_discrete(limits=rev) +
     facet_grid(.~CHR) +
     scale_fill_gradient2(low="#b2182b", mid="#fee090", high="#4575b4", midpoint=0.5) +
     theme_bw() + theme(legend.position="bottom")

ggsave("MHco3_ISE.N1_haplotype_coverage_plot_by_chromosome.pdf")
```


* * *

## Repeats <a name="repeats"></a>

-   repeat masker was run on:
    -   Haem V4
    -   Haem V4 HAPLOTYPES
    -   Haem V1
-   below is example command used.

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



    # get coordinates of repeat types,  output as a bed file

    awk '{print $11}'  HAEM_V4_final.chr.fa.out | sort | uniq | while read RPT; do OUT="$( echo $RPT | sed 's/\//_/g' )" ; grep  "$RPT" HAEM_V4_final.chr.fa.out | awk '{print $5, $6, $7}' OFS="\t" > $OUT.repeats.bed ; done


    # make some windows in the genome - 100kb
samtools faidx  HAEM_V4_final.chr.fa
cut -f1, 2 HAEM_V4_final.chr.fa.fai > HAEM_V4.genome
bedtools-2 makewindows -g HAEM_V4.genome -w 200000 > HAEM_V4.200kb_windows.bed


    # calculate coverage of repeats
for i in *repeats.bed; do
	NAME="$( echo $i | sed 's/.repeats.bed//g' )" ;
	COUNT="$( cat ${i} | wc -l )";
	bedtools-2 coverage -a HAEM_V4.200kb_windows.bed -b ${i} | awk -v NAME=$NAME -v COUNT=$COUNT '{print $0, NAME" n="COUNT}' OFS="\t";
	done > genome_repeat_coverage.200k.data

### make the figure in R
```R
# load libraries
library(ggplot2)

# read data and remove mtDNA
data  <-  read.table("genome_repeat_coverage.200k.data", header=F, sep="\t")
data  <-  data[data$V1!="hcontortus_chr_mtDNA_arrow_pilon", ]

# fix the chromosome labels
chromosome.labels  <-  c("I", "II", "III", "IV", "V",  "X" )
names(chromosome.labels)  <-  c("hcontortus_chr1_Celeg_TT_arrow_pilon",
	"hcontortus_chr2_Celeg_TT_arrow_pilon",
	"hcontortus_chr3_Celeg_TT_arrow_pilon",
	"hcontortus_chr4_Celeg_TT_arrow_pilon",
	"hcontortus_chr5_Celeg_TT_arrow_pilon",
	"hcontortus_chrX_Celeg_TT_arrow_pilon")


ggplot(data, aes(V2, V7, fill=V8)) + geom_area() + facet_grid(V1~.)

#--- remove the "Unknown" class,  as they are abundant,  but dont add much
data2 <- data[data$V8!="Unknown n=205593", ]

# make the plot
ggplot(data2, aes(V2, V7, fill=V8)) +
	geom_area() +
	facet_grid(V1~., labeller = labeller(V1 = chromosome.labels)) +
	theme_bw() +
	labs(x="Genomic position (bp)",
       	y="Fraction of 200 kb window containing repeat type",
       	fill="Repeat type") +
	guides(fill=guide_legend(ncol=1)) +
	ylim(0, 1)
```
