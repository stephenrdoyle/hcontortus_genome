# Haemonchus contortus genome paper
## Section 3: Generation of a high-quality transcriptome annotation incorporating short and long reads

1. [Short-read RNAseq](#srRNAseq)
2. [Long-read RNAseq](#lrRNAseq)
3. [PASA](#pasa)
4. [Exonerate](#exonerate)
5. [EvidenceModeller](#evidencemodeller)
6. [Annotation QC - Sensitivity and specificity](#qc_ss)
7. [Transcriptome Summary Stats](#summarystats)
8. [Gene model plotter](#gene_model_plotter)
9. [Orthology](#orthology)
10. [Other](#other)



******
## Short-read RNAseq <a name="srRNAseq"></a>
```bash
# --- RNAseq data - info
cat /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/TRANSCRIPTOME/RAW/lanes.list

eggs1_236476_3881079	7059_6#1
eggs2_236476_3881080	7059_6#2
eggs3_236476_3881081	7059_6#3
L1_1_236476_3881082	7059_6#4
L1_2_236476_3881083	7059_6#5
L1_3_236476_3881084	7059_6#6
L4_1_236476_3881085	7059_6#7code
L4_2_236476_3881086	7059_6#8
L4_3_236476_3881087	7059_6#9
AdultF1_236476_3881088	7059_6#10
AdultF2_236476_3881089	7059_6#11
AdultF3_236476_3881090	7059_6#12
ISE_BXWF1_236476_3881091	7059_6#13
ISE_BXWF2_236476_3881092	7059_6#14
ISE_BXWF3_236476_3881093	7059_6#15
AdultM1_236476_2305_3914303	7062_6#1
AdultM2_236476_2305_3914304	7062_6#2
AdultM3_236476_2305_3914305	7062_6#3
BXC1_236476_2292_3914306	7062_6#4
BXC2_236476_2292_3914307	7062_6#5
BXC3_236476_2292_3914308	7062_6#6
EXL3_1_236476_2305_3914309	7062_6#7
EXL3_2_236476_2305_3914310	7062_6#8
EXL3_3_236476_2305_3914311	7062_6#9
SHL3_1_236476_2305_3914312	7062_6#10
SHL3_2_236476_2305_3914313	7062_6#11
SHL3_3_236476_2305_3914314	7062_6#12
gut1_236476_1517_3914315	7062_6#13
gut2_236476_1589_3914316	7062_6#14
gut3_236476_635J_3914317	7062_6#15



#-----------------------------------------------------------------------------------------

mkdir STAR_MAP_ALL

cd STAR_MAP_ALL

ln -s /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/REF/HAEM_V4_final.all.fa

mkdir /nfs/users/nfs_s/sd21/lustre118_link/REFERENCE_SEQUENCES/haemonchus_contortus/hc_v4_rnaseq_star_index

#--- make reference index
bsub.py --threads 8 20 01_star_index \
~sd21/lustre118_link/software/bin/star_2.5.2 \
--runMode genomeGenerate \
--runThreadN 8 \
--genomeDir /nfs/users/nfs_s/sd21/lustre118_link/REFERENCE_SEQUENCES/haemonchus_contortus/hc_v4_rnaseq_star_index \
--genomeFastaFiles /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/REF/HAEM_V4_final.all.fa


#--- map reads

for i in ` cd ../RAW/ && ls -1d */ | sed -e 's/\///g' `; do \
bsub.py 10 --threads 8 starmap_${i} /nfs/users/nfs_s/sd21/lustre118_link/software/TRANSCRIPTOME/STAR/bin/Linux_x86_64/STAR \
--runThreadN 8 \
--genomeDir /nfs/users/nfs_s/sd21/lustre118_link/REFERENCE_SEQUENCES/haemonchus_contortus/hc_v4_rnaseq_star_index \
--readFilesIn ../RAW/${i}/${i}_1.fastq.gz ../RAW/${i}/${i}_2.fastq.gz \
--readFilesCommand zcat \
--alignIntronMin 10 \
--outTmpDir starmap_${i}_tmp \
--outFileNamePrefix starmap_${i}_other_out \
--outSAMtype BAM SortedByCoordinate \
; done



# merge bams into a single bam for braker

ls -1 *bam > bams.list
bsub.py 5 03_bammerge "samtools-1.3 merge -b bams.list starmap_merge.bam"


cd ../


#--- repeat for chromosome only

mkdir STAR_MAP_CHR

cd STAR_MAP_CHR

ln -s /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/REF/HAEM_V4_final.chr.fa

 mkdir /nfs/users/nfs_s/sd21/lustre118_link/REFERENCE_SEQUENCES/haemonchus_contortus/hc_v4chr_rnaseq_star_index

#--- make reference index
bsub.py --threads 8 20 01_star_index_chr \
~sd21/lustre118_link/software/bin/star_2.5.2 \
--runMode genomeGenerate \
--runThreadN 8 \
--genomeDir /nfs/users/nfs_s/sd21/lustre118_link/REFERENCE_SEQUENCES/haemonchus_contortus/hc_v4chr_rnaseq_star_index \
--genomeFastaFiles /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/REF/HAEM_V4_final.chr.fa


#--- map reads

for i in ` cd ../RAW/ && ls -1d */ | sed -e 's/\///g' `; do \
bsub.py 10 --threads 8 starmap_${i} /nfs/users/nfs_s/sd21/lustre118_link/software/TRANSCRIPTOME/STAR/bin/Linux_x86_64/STAR \
--runThreadN 8 \
--genomeDir /nfs/users/nfs_s/sd21/lustre118_link/REFERENCE_SEQUENCES/haemonchus_contortus/hc_v4chr_rnaseq_star_index \
--readFilesIn ../RAW/${i}/${i}_1.fastq.gz ../RAW/${i}/${i}_2.fastq.gz \
--readFilesCommand zcat \
--alignIntronMin 10 \
--outTmpDir starmap_${i}_tmp \
--outFileNamePrefix starmap_${i}_other_out \
--outSAMtype BAM SortedByCoordinate \
; done



# merge bams into a single bam for braker

ls -1 *bam > bams.list
bsub.py 5 03_bammerge "samtools-1.3 merge -b bams.list starmap_merge.bam"


cd ../


# make intro hints from RNAseq

#- make hints files
bsub.py 10 01_bam2hints_RNAseq \
"/lustre/scratch118/infgen/archive/ss34/SCHISTO/augustus/augustus-3.2.2//bin/bam2hints --intronsonly --minintronlen=20 --source=E --in=starmap_merge.chr.bam --out=starmap_merge.chr.intron-hints.gff"


samtools-1.3 faidx HAEM_V4_final.chr.fa
awk '{print $1,$2}' OFS="\t" HAEM_V4_final.chr.fa.fai >chromosome_size.list

samtools-1.3 index -b starmap_merge.chr.bam
/software/pathogen/external/apps/usr/local/Python-2.7.13/bin//bam2wig.py -i starmap_merge.chr.bam -s chromosome_size.list -o starmap_merge.chr_bam2wig_out; \
cat starmap_merge.chr_bam2wig_out | \
/software/pathogen/external/apps/usr/bin/wig2hints.pl --width=10 --margin=10 --minthresh=2 --minscore=4 --prune=0.1 --src=W --type=ep --radius=4.5 --pri=4 --strand="." > RNAseq.merge.exon-parts.gff
```


```bash
#-----------------------------------------------------------------------------------------
# Braker
#-----------------------------------------------------------------------------------------

mkdir BRAKER_ALL
mkdir BRAKER_CHR

cd BRAKER_ALL
ln -s /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/TRANSCRIPTOME/STAR_MAP/starmap_merge.bam
ln -s /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/REF/HAEM_V4_final.all.fa

# paper: http://www.ncbi.nlm.nih.gov/pubmed/26559507



export AUGUSTUS_CONFIG_PATH=/nfs/users/nfs_s/sd21/software/augustus-3.2.1/config
export AUGUSTUS_SCRIPTS_PATH=/nfs/users/nfs_s/sd21/software/augustus-3.2.1/scripts
export BAMTOOLS_PATH=/nfs/users/nfs_s/sd21/lustre118_link/software/bamtools/bin
export GENEMARK_PATH=/nfs/users/nfs_s/sd21/lustre118_link/software/gm_et_linux_64/gmes_petap
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/software/gcc-4.9.2/lib64/libstdc++.so.6
export BAMTOOLS_PATH=/nfs/users/nfs_s/sd21/lustre118_link/software/bamtools/bin

### need to copy the gm_key to home directory for it to work. Only needs to be done once. Has a 400 day expiry
#cp ~sd21/lustre118_link/software/gm_et_linux_64/gmes_petap/gm_key ~/.gm_key


bsub.py --queue hugemem --threads 30 200 01_braker_all \
/lustre/scratch118/infgen/team133/sd21/software/TRANSCRIPTOME/BRAKER_v2.0/braker.pl \
--genome=HAEM_V4_final.all.fa \
--bam=starmap_merge.all.bam \
--cores 30 \
--gff3 \
--species=Hc_V4_all \
--UTR=off \
--overwrite \
--workingdir=$PWD \
--useexisting


cd ../BRAKER_CHR/

ln -s /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/REF/HAEM_V4_final.chr.fa
ln -s /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/TRANSCRIPTOME/STAR_MAP_CHR/starmap_merge.chr.bam

bsub.py --queue long --threads 30 30 01_braker_chr \
/lustre/scratch118/infgen/team133/sd21/software/TRANSCRIPTOME/BRAKER_v2.0/braker.pl \
--genome=HAEM_V4_final.chr.fa \
--bam=starmap_merge.chr.bam \
--cores 30 \
--gff3 \
--species=Hc_V4_chr \
--UTR=off \
--overwrite \
--workingdir=$PWD \
--useexisting
```











[↥ **Back to top**](#top)



******
## Long-read RNAseq <a name="lrRNAseq"></a>
```bash
#-----------------------------------------------------------------------------------------
# Running IsoSeq3 pipeline
#-----------------------------------------------------------------------------------------

# Step 1: run ccs
# --- rawdata
cd /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/TRANSCRIPTOME/ISOSEQ_ISOFORMS/RAW

ls -1 *.subreads.bam > bams.list

cd /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/TRANSCRIPTOME/ISOSEQ_ISOFORMS/POOLED

samtools merge -b bams.list merged_subreads.bam

bsub.py --queue long --threads 10 3 01_ccs_pool "ccs merged_subreads.bam merged_subreads.ccs.bam  --noPolish --minPasses 1 --numThreads 10"


# Step 2. run lima
bsub.py --threads 7 1 02_lima_pool "lima merged_subreads.ccs.bam primers.fa merged_subreads.demux.bam --isoseq --no-pbi --num-threads 7";


# Step 3. run isoseq3 cluster
for i in *.demux.primer_5p--primer_3p.bam; do
bsub.py --threads 7 2 03_isoseq3_cluster "isoseq3 cluster ${i} ${i%.demux.primer_5p--primer_3p.bam}.isoseq3_unpolished.bam --verbose --num-threads 7";
done


# Step4. run isoseq3 polish
for i in *isoseq3_unpolished.bam; do
bsub.py --threads 7 10 04_isoseq3_polish "isoseq3 polish ${i} ${i%.isoseq3_unpolished.bam}*subreads.bam ${i%.isoseq3_unpolished.bam}.isoseq3_polished.bam --verbose --num-threads 7";
done


# Step 5 . summarise data
for i in *.isoseq3_polished.bam; do
isoseq3 summarize $i $(i}.summary.csv; done
```




[↥ **Back to top**](#top)





******
## PASA <a name="pasa"></a>
```bash

# connect to mysql server
/usr/bin/mysql -usd21 -pIgebie5Ahbah --port=3307 -hutlt-db
# handy sql cheatsheet : https://gist.github.com/hofmannsven/9164408

cd /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/TRANSCRIPTOME/RAW/ISOSEQ
cat RSII_polished_high_qv_consensus_isoforms.fasta sequel_polished_high_qv_consensus_isoforms.fasta > all_HQ_isoseq.fasta
fastaq enumerate_names --suffix _HQ_isoform all_HQ_isoseq.fasta all_HQ_isoseq.renamed.fasta

cd /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/TRANSCRIPTOME
mkdir PASA_CHR
cd PASA_CHR

ln -s ../BRAKER_CHR/POST_BRAKER/augustus.filtered.gff3
ln -s ../../REF/HAEM_V4_final.chr.fa
ln -s /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/TRANSCRIPTOME/RAW/ISOSEQ/all_isoseq.fasta
fastaq enumerate_names --suffix _hc_isoseq all_isoseq.fasta all_isoseq.renamed.fasta
ln -s /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/TRANSCRIPTOME/RAW/ISOSEQ/all_HQ_isoseq.renamed.fasta

grep ">" all_HQ_isoseq.renamed.fasta | sed 's/>//g' > all_HQ_isoseq.renamed.names

# need pasa config files. modify the MySQL database: sd21_pasa_HcV4
cp ../../V3/TRANSCRIPTOME/PASA/FINAL_V3/alignAssembly.config .
cp ../../V3/TRANSCRIPTOME/PASA/FINAL_V3/annotCompare.config .

bsub.py --queue yesterday --threads 7 20 01_pasa_isoseq_HcV4_chr \
"/nfs/users/nfs_s/sd21/lustre118_link/software/TRANSCRIPTOME/PASApipeline-pasa-v2.2.0/scripts/Launch_PASA_pipeline.pl \
-c alignAssembly.config \
-C -R \
-g HAEM_V4_final.chr.fa \
-t all_HQ_isoseq.renamed.fasta \
-f all_HQ_isoseq.renamed.names \
--ALIGNERS blat,gmap --CPU 7"


bsub.py --queue yesterday 1 02_pasa_add_evm_annotations \
"/nfs/users/nfs_s/sd21/lustre118_link/software/TRANSCRIPTOME/PASApipeline-pasa-v2.2.0/scripts/Load_Current_Gene_Annotations.dbi \
-c alignAssembly.config \
-g HAEM_V4_final.chr.fa \
-P augustus.filtered.gff3"



bsub.py --queue yesterday 1 03_pasa_evm_compare_update \
"/nfs/users/nfs_s/sd21/lustre118_link/software/TRANSCRIPTOME/PASApipeline-pasa-v2.2.0/scripts/Launch_PASA_pipeline.pl \
-c annotCompare.config \
-A -L \
--annots_gff3 augustus.filtered.gff3 \
--ALT_SPLICE \
-g HAEM_V4_final.chr.fa \
-t all_HQ_isoseq.renamed.fasta \
-f all_HQ_isoseq.renamed.names"



bsub.py --queue yesterday 1 02_pasa_add_evm_annotations_round2 \
"/nfs/users/nfs_s/sd21/lustre118_link/software/TRANSCRIPTOME/PASApipeline-pasa-v2.2.0/scripts/Load_Current_Gene_Annotations.dbi \
-c alignAssembly.config \
-g HAEM_V4_final.chr.fa \
-P sd21_pasa_HcV4.gene_structures_post_PASA_updates.39527.gff3"


bsub.py --queue yesterday 1 03_pasa_evm_compare_update_round2 \
"/nfs/users/nfs_s/sd21/lustre118_link/software/TRANSCRIPTOME/PASApipeline-pasa-v2.2.0/scripts/Launch_PASA_pipeline.pl \
-c annotCompare.config \
-A -L \
--annots_gff3 sd21_pasa_HcV4.gene_structures_post_PASA_updates.39527.gff3 \
--ALT_SPLICE \
-g HAEM_V4_final.chr.fa \
-t all_HQ_isoseq.renamed.fasta \
-f all_HQ_isoseq.renamed.names"
```





[↥ **Back to top**](#top)



******
## Exonerate <a name="exonerate"></a>

```bash
cd /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/TRANSCRIPTOME
mkdir EXONERATE_CHR
cd EXONERATE_CHR
ln -s /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/REF/HAEM_V4_final.chr.fa
mkdir V1_2_V4
cd V1_2_V4

wget ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS9/species/haemonchus_contortus/PRJEB506/haemonchus_contortus.PRJEB506.WBPS9.protein.fa.gz
gunzip haemonchus_contortus.PRJEB506.WBPS9.protein.fa.gz


bsub.py --queue yesterday 10 01_exonserate \
/nfs/users/nfs_s/sd21/bash_scripts/run_exonerate_splitter \
../HAEM_V4_final.chr.fa \
haemonchus_contortus.PRJEB506.WBPS9.protein.fa
```



[↥ **Back to top**](#top)





******
## EvidenceModeller <a name="evidencemodeller"></a>
```bash
#-----------------------------------------------------------------------------------------
# Evidence Modeller
#-----------------------------------------------------------------------------------------
# want to merge
#--- pasa output, braker(filter), transposon filter

mkdir /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/TRANSCRIPTOME/EVM_CHR
cd /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/TRANSCRIPTOME/EVM_CHR

ln -s ../../REF/HAEM_V4_final.chr.fa




cat exonerate.V1_2_V4.gff | sed '/^$/d' | awk '{$2 = "hcv1_exonerate"; print}' OFS="\t" > exonerate.V1_2_V4.renamed.gff
cat ce_2_v4.exonerate.gff3 | sed '/^$/d' | awk '{$2 = "ce_exonerate"; print}' OFS="\t" > ce_2_v4.exonerate.renamed.gff3
cat HC_V4_augustus_merge.AAfiltered.gff | sed '/^$/d' | awk '{$2 = "braker_augustus"; print}' OFS="\t" > braker.renamed.gff3
cat sd21_pasa_HcV4_2.gene_structures_post_PASA_updates.18752.gff3 | sed '/^$/d' | grep -v "#" | awk '{$2 = "pasa"; print}' OFS="\t" > pasa.renamed.gff3


# weights.txt - second column in weights file must match second column in data files
PROTEIN	hcv1_exonerate	2
PROTEIN	ce_exonerate	1
PROTEIN curated_exonerate 2
TRANSCRIPT	pasa	5
ABINITIO_PREDICTION	braker_augustus	2

bsub.py 5 01_evm \
"/nfs/users/nfs_s/sd21//lustre118_link/software/TRANSCRIPTOME/EVidenceModeler-1.1.1/EvmUtils/partition_EVM_inputs.pl \
--genome HAEM_V4_final.chr.fa \
--gene_predictions braker.renamed.gff3 \
--transcript_alignments pasa.renamed.gff3 \
--protein_alignments exonerate.V1_2_V4.renamed.gff \
--protein_alignments ce_2_v4.exonerate.renamed.gff3 \
--protein_alignments AUGUSTUS_JN_CURATED.exonerate.renamed.gff \
--segmentSize 100000 \
--overlapSize 10000 \
--partition_listing partitions_list.out"




bsub.py 5 02_evm_write_commands \
"/nfs/users/nfs_s/sd21//lustre118_link/software/TRANSCRIPTOME/EVidenceModeler-1.1.1/EvmUtils/write_EVM_commands.pl \
--weights /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/TRANSCRIPTOME/EVM_CHR/weights.txt \
--genome HAEM_V4_final.chr.fa \
--gene_predictions braker.renamed.gff3 \
--transcript_alignments pasa.renamed.gff3 \
--protein_alignments exonerate.V1_2_V4.renamed.gff \
--protein_alignments ce_2_v4.exonerate.renamed.gff3 \
--protein_alignments AUGUSTUS_JN_CURATED.exonerate.renamed.gff \
--output_file_name evm.out \
--partitions partitions_list.out \> commands.list"

chmod a+x commands.list
bsub.py 5 03_evm_run "./commands.list"

bsub.py 5 04_evm_recombine_partitions \
"/nfs/users/nfs_s/sd21//lustre118_link/software/TRANSCRIPTOME/EVidenceModeler-1.1.1/EvmUtils/recombine_EVM_partial_outputs.pl \
--partitions partitions_list.out \
--output_file_name evm.out"

bsub.py 5 05_evm_make_gff3 \
"/nfs/users/nfs_s/sd21//lustre118_link/software/TRANSCRIPTOME/EVidenceModeler-1.1.1/EvmUtils/convert_EVM_outputs_to_GFF3.pl \
--partitions partitions_list.out \
--output evm.out \
--genome HAEM_V4_final.chr.fa"


find . -name "evm.out.gff3" | sort | while read -r line; do cat $line >> HAEM_V4.chr.evm_merge.gff; done

```




[↥ **Back to top**](#top)





******
## Annotation QC - Sensitivity and specificity <a name="qc_ss"></a>

```bash
Want to compare the final annotation to steps along the way. These include
- V1 genome vs manually curated V1 genes
_ V4 final vs AUGUSTUS
- V4 final vs BRAKER
- V4 final vs PASA w Isoseq R1
- V4 final vs EVM w Isoseq R2

Need to do this with the 110 manually curated genes, and the whole V4 final.

### Working environment
```shell
cd /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/TRANSCRIPTOME
mkdir TRANSCRIPTOME_QC
cd TRANSCRIPTOME_QC
```


### Get some data
```shell
# V1 curated genes
cp /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/V3/TRANSCRIPTOME/V1_MANUAL_CURATION/ABC_LGic.2.gff ABC_LGic.2.gff

# V1 GENOME
cp ~sd21/lustre118_link/REFERENCE_SEQUENCES/haemonchus_contortus/V1/haemonchus_contortus.PRJEB506.WBPS8.genomic.fa HCON_V1.fa

# V1 Genome annotation
cp /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/V3/TRANSCRIPTOME/V1_MANUAL_CURATION/Hc_rztk_1+2+8+9.augustus.gff3 Hc_rztk_1+2+8+9.augustus.gff3

cp ~sd21/lustre118_link/REFERENCE_SEQUENCES/haemonchus_contortus/V1/haemonchus_contortus.PRJEB506.WBPS8.annotations.gff3 HCON_V1.annotation.gff3

# V4 final annotaiton
cp ~sd21/lustre118_link/hc/GENOME/TRANSCRIPTOME/TRANSCRIPTOME_CURATION/HCON_V4_WBP11plus_190125.ips.gff3 HCON_V4_FINAL.gff3

# V4 AUGUSTUS

# V4 BRAKER
cp /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/TRANSCRIPTOME/BRAKER_CHR/braker/Hc_V4_chr/augustus.gff3 HCON_V4_BRAKER.gff3

# V4 PASA
cp /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/TRANSCRIPTOME/PASA_CHR/sd21_pasa_HcV4_2.gene_structures_post_PASA_updates.18752.gff3 HCON_V4_PASA.gff3

# V4 EMV
cp /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/TRANSCRIPTOME/PASA_CHR_R2/HCON_V4.renamed.gff3 HCON_V4_EVM.gff3
```



### Transcriptome comparisons
```shell
# V1 vs curated V1 genes
gffcompare -R -r ABC_LGic.2.gff -o V1cutated_vs_V1annotation Hc_rztk_1+2+8+9.augustus.gff3

# Result
#-----------------| Sensitivity | Precision  |
        Base level:    92.9     |     2.7    |
        Exon level:    86.3     |     2.6    |
      Intron level:    89.5     |     2.7    |
Intron chain level:    33.3     |     0.7    |
  Transcript level:    34.3     |     0.7    |

       Locus level:    34.8     |     0.8    |


# V4 Final vs BRAKER
gffcompare -R -r HCON_V4_FINAL.gff3 -o V4_FINAL_vs_BRAKER HCON_V4_BRAKER.gff3
# Result

       #-----------------| Sensitivity | Precision  |
               Base level:    85.4     |    79.7    |
               Exon level:    89.6     |    86.7    |
             Intron level:    94.8     |    93.9    |
       Intron chain level:    67.5     |    61.4    |
         Transcript level:    72.7     |    56.3    |
              Locus level:    78.0     |    60.6    |

# V4 Final vs PASA_CHR
gffcompare -R -r HCON_V4_FINAL.gff3 -o V4_FINAL_vs_PASA HCON_V4_PASA.gff3

       #-----------------| Sensitivity | Precision  |
               Base level:    85.0     |    97.7    |
               Exon level:    89.6     |    95.4    |
             Intron level:    90.0     |    96.5    |
       Intron chain level:    33.3     |    32.9    |
         Transcript level:    33.3     |    30.8    |
              Locus level:    30.3     |    29.8    |

# V4 Final vs EVM
gffcompare -R -r HCON_V4_FINAL.gff3 -o V4_FINAL_vs_EVM HCON_V4_EVM.gff3

       #-----------------| Sensitivity | Precision  |
               Base level:    92.9     |    98.7    |
               Exon level:    94.7     |    96.1    |
             Intron level:    96.3     |    98.0    |
       Intron chain level:    84.1     |    84.1    |
         Transcript level:    86.7     |    86.6    |
              Locus level:    87.5     |    86.4    |
```

Don't think it is worth including the V1 comparison, as these curated genes would have been incorporated into the final annotation. Not really a good comparison. V1 precision is low due to only a subset of genes being used.

Results suggest:
- no increase in sensitivity, but big increase in precision from BRAKER to PASA
- increase in Sensitivity and Precision from PASA to EVM
```

[↥ **Back to top**](#top)





******
## Transcriptome summary stats <a name="summarystats"></a>

```bash

## Annotation quantitative quantitative data
cd ~/lustre118_link/hc/GENOME/TRANSCRIPTOME/TRANSCRIPTOME_CURATION/
mkdir HCON_V4_WBP11plus_190125_ANALYSIS
cd HCON_V4_WBP11plus_190125_ANALYSIS

ln -sf ~sd21/lustre118_link/hc/GENOME/TRANSCRIPTOME/TRANSCRIPTOME_CURATION/HCON_V4_WBP11plus_190125.ips.gff3
gag.py -f ../HAEM_V4_final.chr.fa -g HCON_V4_WBP11plus_190125.ips.gff3
```




#--------------------------------
# Transcriptome Summary Stats summary stats

cd /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/TRANSCRIPTOME/TRANSCRIPTOME_SUMMARY_STATS

for i in V1 MCMASTER V4 CELEGANS V4_190114; do awk -v name="$i" '$3=="gene" {print name,"gene",$5-$4}' OFS="\t" ${i}/*.gff* > ${i}/${i}.genelength.txt; awk -v name="$i" '$3=="mRNA" {print name,"mrna",$5-$4}' OFS="\t" ${i}/*.gff* > ${i}/${i}.mrnalength.txt; awk -v name="$i" '$3=="exon" {print name,"exon",$5-$4}' OFS="\t" ${i}/*.gff* > ${i}/${i}.exonlength.txt; awk -v name="$i" '$3=="CDS" {print name,"cds",$5-$4}' OFS="\t" ${i}/*.gff* > ${i}/${i}.cdslength.txt; done

for i in V4_190114; do awk -v name="$i" '$3=="gene" {print name,"gene",$5-$4}' OFS="\t" ${i}/*.gff* > ${i}/${i}.genelength.txt; awk -v name="$i" '$3=="mRNA" {print name,"mrna",$5-$4}' OFS="\t" ${i}/*.gff* > ${i}/${i}.mrnalength.txt; awk -v name="$i" '$3=="exon" {print name,"exon",$5-$4}' OFS="\t" ${i}/*.gff* > ${i}/${i}.exonlength.txt; awk -v name="$i" '$3=="CDS" {print name,"cds",$5-$4}' OFS="\t" ${i}/*.gff* > ${i}/${i}.cdslength.txt; done



# get intron lengths
gt gff3 -addintrons haemonchus_contortus.PRJEB506.WBPS8.annotations.gff3 | awk '$3=="intron" {print "V1","intron",$5-$4}' OFS="\t" > V1.intronlength.txt
#gt gff3 -tidy -addintrons HCON_V4.renamed.gff3 | awk '$3=="intron" {print "V4","intron",$5-$4}' OFS="\t" > V4.intronlength.txt
gt gff3 -tidy -addintrons wormbase.20240.complete.gff3 | awk '$3=="intron" {print "CELEGANS","intron",$5-$4}' OFS="\t" > celegans.intronlength.txt
gt gff3 -addintrons haemonchus_contortus.PRJNA205202.WBPS9.annotations.gff3 | awk '$3=="intron" {print "MCMASTER","intron",$5-$4}' OFS="\t" > MCMASTER.intronlength.txt
gt gff3 -tidy -addintrons HCON_V4_WBP11plus_190114.gff3 | awk '$3=="intron" {print "V4","intron",$5-$4}' OFS="\t" > V4_190114.intronlength.txt

# collate data
cat CELEGANS/CELEGANS.genelength.txt MCMASTER/MCMASTER.genelength.txt V1/V1.genelength.txt V4_190114/V4_190114.genelength.txt > genelength.txt
cat CELEGANS/CELEGANS.mrnalength.txt MCMASTER/MCMASTER.mrnalength.txt V1/V1.mrnalength.txt V4_190114/V4_190114.mrnalength.txt > mrnalength.txt
cat CELEGANS/CELEGANS.exonlength.txt MCMASTER/MCMASTER.exonlength.txt V1/V1.exonlength.txt V4_190114/V4_190114.exonlength.txt > exonlength.txt
cat CELEGANS/CELEGANS.cdslength.txt MCMASTER/MCMASTER.cdslength.txt V1/V1.cdslength.txt V4_190114/V4_190114.cdslength.txt > cdslength.txt
cat CELEGANS/celegans.intronlength.txt MCMASTER/MCMASTER.intronlength.txt V1/V1.intronlength.txt V4_190114/V4_190114.intronlength.txt > intronlength.txt
```

```R
# Summary stats

library(ggplot2)
library(patchwork)

gene<-read.table("genelength.txt",header=F)
mRNA<-read.table("mrnalength.txt",header=F)
exon<-read.table("exonlength.txt",header=F)
intron<-read.table("intronlength.txt",header=F)


gene_plot <- ggplot()+geom_density(aes(log10(gene$V3),col=gene$V1,fill=gene$V1),alpha = 0.2)+theme_bw()+theme(legend.position="bottom")+labs(title ="Gene length", x = "Length (log10[bp])", y = "Density")
mRNA_plot <- ggplot()+geom_density(aes(log10(mRNA$V3),col=mRNA$V1,fill=mRNA$V1),alpha = 0.2)+theme_bw()+ theme(legend.position="none")+labs(title ="mRNA length", x = "Length (log10[bp])", y = "Density")
exon_plot <- ggplot()+geom_density(aes(log10(exon$V3),col=exon$V1,fill=exon$V1),alpha = 0.2)+theme_bw()+ theme(legend.position="none")+labs(title ="Exon length", x = "Length (log10[bp])", y = "Density")
cds_plot <- ggplot()+geom_density(aes(log10(cds$V3),col=cds$V1,fill=cds$V1),alpha = 0.2)+theme_bw()+ theme(legend.position="none")+labs(title ="CDS length", x = "Length (log10[bp])", y = "Density")
intron_plot <- ggplot()+geom_density(aes(log10(intron$V3),col=intron$V1,fill=intron$V1),alpha = 0.2)+theme_bw()+ theme(legend.position="none")+labs(title ="Intron length", x = "Length (log10[bp])", y = "Density")



gene_plot + mRNA_plot + exon_plot + intron_plot + plot_layout(ncol = 4)
ggsave(filename="transcriptome_stats.pdf", width=35,height=10,units="cm")
```


```
library(ggplot2)
data<-read.table("summary_stats.txt",sep="\t",header=T)


ggplot()+
     geom_point(aes(x=log10(data$count),y=log10(data$mean),col=data$class,shape=data$Species),size=3)+
     theme_bw()+
     labs(y="Mean length (log10[bp])",x="Feature count (log10[total])")

ggsave("annotation_comparison_4species_scatter.pdf",useDingbats=F)
ggsave("annotation_comparison_4species_scatter.png")
```



## 03 - Gene model plotter <a name="gene_model_plotter"></a>

Working environment and data

```shell
cd /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/TRANSCRIPTOME/GENE_MODEL_PLOTS

#gff
ln -fs /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/TRANSCRIPTOME/TRANSCRIPTOME_CURATION/HCON_V4_WBP11plus_190125.ips.gff3 ANNOTATION.gff
#ccs subread data

# generate bed files of isoseq reads in bam file : /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/TRANSCRIPTOME/ISOSEQ_ISOFORMS/POOLED
#bedtools-2 bamtobed -cigar -split -i pooled_ccs_subreads.sorted.bam > pooled_ccs_subreads.sorted.bed
#bedtools-2 bamtobed -cigar -i pooled_ccs_subreads.sorted.bam > pooled_ccs_subreads.sorted.whole.bed

ln -s /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/TRANSCRIPTOME/ISOSEQ_ISOFORMS/POOLED/pooled_ccs_subreads.sorted.bed
ln -s /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/TRANSCRIPTOME/ISOSEQ_ISOFORMS/POOLED/pooled_ccs_subreads.sorted.whole.bed

# old V1 annotations
ln -s /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/TRANSCRIPTOME/EXONERATE_CHR/V1_2_V4/exonerate.V1_2_V4.gff OLD_ANNOTATIONS.gff

```



```R
#install.packages("data.table")
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("ggbio", version = "3.8")

library(data.table)
library(ggplot2)
library(dplyr)

# load data
# --- genome annotation
gff<-fread(cmd="grep ^hcontortus ANNOTATION.gff")
colnames(gff) <- c("chr","source","feature","start","end","point1","strand","frame","info")

# --- isoseq data - gaps
iso_gap<-fread("pooled_ccs_subreads.sorted.bed")
colnames(iso_gap) <- c("chr","start","end","readname","score","frame")

# --- isoseq data - full length
iso_whole<-fread("pooled_ccs_subreads.sorted.whole.bed")
colnames(iso_whole) <- c("chr","start","end","readname","score","frame","cigar")

# --- V1 genome annotation GFF
old_gff <- fread(cmd="grep exonerate OLD_ANNOTATIONS.gff",sep="\t", sep2=";")
colnames(old_gff) <- c("chr","source","feature","start","end","point1","strand","frame","info")



# select gene ID
gene='HCON_00107450'

#btub1 HCON_00005260


gene_model_plot <- function(gene){

# filter data to select chromosome and mRNA
mrna_data <- gff[grep(gene, gff$info), ]
mrna_data <- mrna_data[mrna_data$feature=='mRNA',]
mrna_data <- cbind(mrna_data, read.table(text = as.character(mrna_data$info), sep = ";"))
mrna_data <- cbind(mrna_data, read.table(text = as.character(mrna_data$V1), sep = "=",col.names=c("ID","unique_ID")))
mrna_id <- head(data.frame(mrna_data$unique_ID),1)  # gives 1st isoform if multiple
chromosome <- mrna_data[1,1]
colnames(chromosome)<-c("chromosome_ID")

data <- gff[grep(mrna_id$mrna_data.unique_ID, gff$info), ]

# filter by feature type
cds <- data[data$feature=="CDS",]
mrna <- data[data$feature=="mRNA",]
#gene <- data[data$feature=="gene",]


intron<-data.frame(head(cds$end,-1),tail(cds$start,-1),(tail(cds$start,-1)-head(cds$end,-1))/2)
colnames(intron)<-c("start","end","midpoint")

utr5<-data.frame(head(mrna$start,1),head(sort(cds$start),1))
colnames(utr5)<-c("start","end")
utr3<-data.frame(head(mrna$end,1),tail(sort(cds$end),1))
colnames(utr3)<-c("start","end")
utr<-rbind(utr5,utr3)



#longest 20 full length ccs reads
iso_gap2 <- iso_gap[(iso_gap$chr==chromosome$chromosome_ID) & (iso_gap$start > (mrna$start-(0.1*(mrna$end-mrna$start)))) & (iso_gap$end < (mrna$end+(0.1*(mrna$end-mrna$start)))),]
iso_whole2 <- iso_whole[(iso_whole$chr==chromosome$chromosome_ID) & (iso_whole$start > (mrna$start-(0.1*(mrna$end-mrna$start)))) & (iso_whole$end < (mrna$end+(0.1*(mrna$end-mrna$start)))),]
iso_whole2$length <- (iso_whole2$end-iso_whole2$start)
iso_whole2_10 <- tail(iso_whole2[order(iso_whole2$length)],20)
iso_whole2_10$rank <- 1:nrow(iso_whole2_10)*1/20*5


test_join <- dplyr:::inner_join(iso_gap2, iso_whole2_10, by = "readname")

# filter V1 gff file
old_gff2 <- old_gff[old_gff$chr==chromosome$chromosome_ID & old_gff$start > (mrna$start-(0.1*(mrna$end-mrna$start))) & old_gff$end < (mrna$end+(0.1*(mrna$end-mrna$start))),]
old_gff2 <- cbind(old_gff2, read.table(text = as.character(old_gff2$info), sep = ";",col.names=c("exonerate_ID","V1_ID")))

if(mrna$strand=="+"){
  arrow <- data.frame(mrna$end,mrna$end+(0.02*(mrna$end-mrna$start)))
  intron<-data.frame(head(cds$end,-1),tail(cds$start,-1),(tail(cds$start,-1)-head(cds$end,-1))/2)

  } else {
  arrow <- data.frame(mrna$start,mrna$start-(0.02*(mrna$end-mrna$start)))
  intron<-data.frame(tail(cds$end,-1),head(cds$start,-1),(head(cds$start,-1)-tail(cds$end,-1))/2)
  }

colnames(arrow) <-  c("start","end")
colnames(intron)  <-  c("start","end","midpoint")

# make plot
ggplot()+
  #geom_rect(data=mrna,aes(xmin=mrna$V4,ymin=0,xmax=mrna$V5,ymax=1),fill="grey90")+
  # new gene model
  geom_rect(data=utr,aes(xmin=utr$start,ymin=0.5,xmax=utr$end,ymax=1.5),fill=NA,col="black",size=0.4)+
  geom_segment(data=intron,aes(x=intron$start,xend=intron$start+intron$midpoint,y=1,yend=0.5),size=0.5)+
  geom_segment(data=intron,aes(x=intron$start+intron$midpoint,xend=intron$end,y=0.5,yend=1),size=0.5)+
  geom_rect(data=cds,aes(xmin=cds$start,ymin=0.5,xmax=cds$end,ymax=1.5),fill="black",col=NA)+
  geom_segment(data=arrow,aes(x=arrow$start,xend=arrow$end,y=1.5,yend=1),size=0.4)+
  geom_segment(data=arrow,aes(x=arrow$start,xend=arrow$end,y=0.5,yend=1),size=0.4)+
  geom_text(aes(x=mrna$end+(0.15*(mrna$end-mrna$start)), y=1, label = "Haem V4 Gene model"))+
  # old gene model
  geom_rect(data=old_gff2,aes(xmin=old_gff2$start,ymin=2.5,xmax=old_gff2$end,ymax=3.5,fill=old_gff2$V1_ID))+
  geom_text(aes(x=mrna$end+(0.15*(mrna$end-mrna$start)), y=3, label = "Haem V1 CDS \n (Exonerate: protein2genome)"))+
  # isoseq
  geom_rect(data=iso_whole2_10,aes(xmin=iso_whole2_10$start,ymin=4+as.numeric(iso_whole2_10$rank)-0.02,xmax=iso_whole2_10$end,ymax=4+as.numeric(iso_whole2_10$rank)+0.02),fill="grey90")+
  geom_rect(data=test_join,aes(xmin=test_join$start.x,ymin=4+test_join$rank-0.08,xmax=test_join$end.x,ymax=4+test_join$rank+0.08),fill="cornflowerblue")+
  geom_text(aes(x=mrna$end+(0.15*(mrna$end-mrna$start)), y=5, label = "IsoSeq cDNA reads \n (minimap2 splice)"))+
  # plot layout
  theme_bw()+
  #xlab("Genome position (bp)")+
  labs(title= paste("Gene ID: ",gene), x =paste("Chromosome: ",chromosome," position (bp)"))+
  xlim(mrna$start-(0.1*(mrna$end-mrna$start)),mrna$end+(0.25*(mrna$end-mrna$start)))+
  scale_y_reverse(lim=c(10,0))+ scale_fill_discrete(guide=FALSE)+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

}

gene_model_plot('')

```




---
## 04 - Orthology <a name="orthology"></a>
---
### Working envoronment
```shell

cd /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/SELECTION
```

#### get data
```shell
wget ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS11/species/haemonchus_placei/PRJEB509/haemonchus_placei.PRJEB509.WBPS11.protein.fa.gz
gunzip haemonchus_placei.PRJEB509.WBPS11.protein.fa.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS11/species/caenorhabditis_elegans/PRJNA13758/caenorhabditis_elegans.PRJNA13758.WBPS11.protein.fa.gz
gunzip caenorhabditis_elegans.PRJNA13758.WBPS11.protein.fa.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS11/species/haemonchus_contortus/PRJNA205202/haemonchus_contortus.PRJNA205202.WBPS11.protein.fa.gz
gunzip haemonchus_contortus.PRJNA205202.WBPS11.protein.fa.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS10/species/haemonchus_contortus/PRJEB506/haemonchus_contortus.PRJEB506.WBPS10.protein.fa.gz
gunzip haemonchus_contortus.PRJEB506.WBPS10.protein.fa.gz

ln -fs ../TRANSCRIPTOME/TRANSCRIPTOME_CURATION/HCON_V4_WBP11plus_190125.ips.gff3
gffread HCON_V4_WBP11plus_190125.ips.gff3 -g HAEM_V4_final.chr.fa -y HCON_V4_WBP11plus_190125.ips.proteins.fa

```

```shell
# curate data for input into orthofinder
# get one coding sequence per gene - certainly Ce and HcV4 has multiple  isoforms, and therefore multiple coding sequneces per gene
fastaq to_fasta -l0 caenorhabditis_elegans.PRJNA13758.WBPS11.protein.fa caenorhabditis_elegans.PRJNA13758.WBPS11.protein.fa2
fastaq to_fasta -l0 haemonchus_placei.PRJEB509.WBPS11.protein.fa haemonchus_placei.PRJEB509.WBPS11.protein.fa2
fastaq to_fasta -l0 haemonchus_contortus.PRJEB506.WBPS10.protein.fa haemonchus_contortus.PRJEB506.WBPS10.protein.fa2
fastaq to_fasta -l0 haemonchus_contortus.PRJNA205202.WBPS11.protein.fa haemonchus_contortus.PRJNA205202.WBPS11.protein.fa2
fastaq to_fasta -l0 HCON_V4_WBP11plus_190125.ips.proteins.fa  HCON_V4_WBP11plus_190125.ips.proteins.fa2

# fix fasta header - likely only c elegans that had excessive informaiton in the header, but fixed to make them all consistent
awk '{if($1 ~ /^>/) print ">"$3,$2; else print $0}' caenorhabditis_elegans.PRJNA13758.WBPS11.protein.fa2 > ce.proteins.fa
awk '{if($1 ~ /^>/) print ">"$3,$2; else print $0}' haemonchus_placei.PRJEB509.WBPS11.protein.fa2 > hp.proteins.fa
awk '{if($1 ~ /^>/) print ">"$3,$2; else print $0}' haemonchus_contortus.PRJEB506.WBPS10.protein.fa2 > hc_V1.proteins.fa
awk '{if($1 ~ /^>/) print ">"$3,$2; else print $0}' haemonchus_contortus.PRJNA205202.WBPS11.protein.fa2 > hc_McM.proteins.fa
cat HCON_V4_WBP11plus_190125.ips.proteins.fa2 > hc_V4.proteins.fa

#remove stop codons from hc_V4
sed -i 's/[A-Z]*\.\.*[A-Z]//g' hc_V4.proteins.fa
sed -i 's/\.$//g' hc_V4.proteins.fa

# make unique gene lists
grep ">" ce.proteins.fa | cut -f 1 -d " " | sort | uniq > unique.Ce_genes.list
#> 20208 unique.Ce_genes.list

grep ">" hp.proteins.fa | cut -f 1 -d " " | sort | uniq > unique.Hp_genes.list
#> 21928 unique.Hp_genes.list

grep ">" hc_McM.proteins.fa | cut -f 1 -d " " | sort | uniq > unique.Hc_McM_genes.list
#> 23610 unique.Hc_McM_genes.list

grep ">" hc_V1.proteins.fa | cut -f 1 -d " " | sort | uniq > unique.Hc_V1_genes.list
#> 21869 unique.Hc_V1_genes.list


grep ">" hc_V4.proteins.fa | cut -f 2 -d " " | sort | uniq > unique.Hc_V4_genes.list
#19438 unique.Hc_V4-1901140_genes.list


while read -r gene; do grep -m1 -A1 ${gene} ce.proteins.fa; done < unique.Ce_genes.list > ce.proteins.unique.fa &
while read -r gene; do grep -m1 -A1 ${gene} hp.proteins.fa; done < unique.Hp_genes.list > hp.proteins.unique.fa &
while read -r gene; do grep -m1 -A1 ${gene} hc_McM.proteins.fa; done < unique.Hc_McM_genes.list > hc_McM.proteins.unique.fa &
while read -r gene; do grep -m1 -A1 ${gene} hc_V1.proteins.fa; done < unique.Hc_V1_genes.list > hc_V1.proteins.unique.fa &
while read -r gene; do grep -m1 -A1 ${gene} hc_V4.proteins.fa; done < unique.Hc_V4_genes.list > hc_V4.proteins.unique.fa &


# run OrthoFinder
#--- setup data
mkdir PROTEIN_FASTAs
mv *.unique.fa PROTEIN_FASTAs
cd PROTEIN_FASTAs
for i in *.fa; do cut -f1 -d " " $i | sed '/^$/d' > tmp; mv tmp $i; done
cd ../

# run OF
bsub.py --queue long --threads 20 20 orthofinder_2.2.7 "python2.7 /nfs/users/nfs_s/sd21/lustre118_link/software/POPGEN/OrthoFinder-2.2.7_source/orthofinder/orthofinder.py -t 20 -a 20 -S diamond -M msa -A mafft -f PROTEIN_FASTAs/"



1:1 orthologs


ce.proteins.unique      hc_McM.proteins.unique  hc_V1.proteins.unique   hc_V4.proteins.unique   hp.proteins.unique

ce.proteins.unique      0.0     4424.0  4529.0  7361.0  6371.0
hc_McM.proteins.unique  4424.0  0.0     6559.0  7581.0  7861.0
hc_V1.proteins.unique   4529.0  6559.0  0.0     9595.0  7991.0
hc_V4.proteins.unique   7361.0  7581.0  9595.0  0.0     9970.0
hp.proteins.unique      6371.0  7861.0  7991.0  9970.0  0.0

#--- KINFIN

cd /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/SELECTION/PROTEIN_FASTAs/Results_Jan25
mkdir KINFIN
cd KINFIN
ln -s ../WorkingDirectory/SpeciesIDs.txt
ln -s ../WorkingDirectory/SequenceIDs.txt
ln -s ../Orthogroups.txt
ln -s ../Orthologues_*/SpeciesTree_rooted.txt


#---- kinfin
#-- made config.txt file containing
echo -e "#IDX,TAXON,OUT
0,ce.proteins.unique,1
1,hc_mcm.proteins.unique,0
2,hc_V1.proteins.unique,0
3,hc_V4.proteins.unique,0
4,hp.proteins.unique,0" > config.txt



# run KINFIN
/nfs/users/nfs_s/sd21/lustre118_link/software/POPGEN/kinfin/kinfin \
--cluster_file Orthogroups.txt \
--config_file config.txt \
--sequence_ids_file SequenceIDs.txt \
--fasta_dir /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/SELECTION/PROTEIN_FASTAs \
--species_ids_file SpeciesIDs.txt \
--tree_file SpeciesTree_rooted.txt

#-- useful - kinfin generated a 1-to-1 ortholog list, with fuzzy matching, which basically means it allows some missingness, and not required to be in all samples

/nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/SELECTION/PROTEIN_FASTAs/Results_Jan23_1/KINFIN/kinfin_results/all/all.all.cluster_1to1s.txt
#> 5746 all.all.cluster_1to1s.txt
#>





# plotting orthogroups - all and 1to1 - using UpSetR
cd /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/SELECTION/PROTEIN_FASTAs/Results_*/KINFIN/kinfin_results/TAXON

cat TAXON.cluster_summary.txt | awk '{print $1,$9,$10,$11,$12,$13}' OFS="\t" | awk 'NR>1 {for(i=2;i<=NF;i++)if($i>0)$i=1}1' OFS="\t" > all_orthogroups.upsetr.data

cat TAXON.cluster_summary.txt | awk '{print $1,$9,$10,$11,$12,$13}' OFS="\t" | head -n1 > 1to1_orthogroups.upsetr.data
cat TAXON.cluster_summary.txt | awk '{print $1,$9,$10,$11,$12,$13}' OFS="\t" | awk 'NR>1{if ($2<=1 && $3 <=1 && $4 <= 1 && $5 <=1 && $6 <=1) print}' OFS="\t" >> 1to1_orthogroups.upsetr.data
```

```R
R-3.5.0
library(UpSetR)
all_orthogroups<-read.table("all_orthogroups.upsetr.data",header=T,comment.char="")
pdf("all_orthogroups_plot.upsetr.pdf",height=5,width=10,useDingbats=FALSE)
upset(all_orthogroups)
dev.off()

png("all_orthogroups_plot.upsetr.png",height=5,width=10)
upset(all_orthogroups)
dev.off()

one2one_orthogroups <- read.table("1to1_orthogroups.upsetr.data",header=T,comment.char="")
pdf("one2one_orthogroups_plot.upsetr.pdf",height=5,width=10,useDingbats=FALSE)
upset(one2one_orthogroups)
dev.off()

png("one2one_orthogroups_plot.upsetr.png",height=5,width=10)
upset(one2one_orthogroups)
dev.off()

```

```shell
scp sd21@pcs5.internal.sanger.ac.uk:/nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/SELECTION/PROTEIN_FASTAs/Results_Jan25/KINFIN/kinfin_results/TAXON/*pdf ~/Documents/workbook/hcontortus_genome/04_analysis
```

















## Other  <a name="other"></a>
This section contains some analysis that I tried by didnt end up using for one reason or another. Likely becasue I focun a better way of doing it.
- running Augustus post braker to incorportate the isoseq
     - didnt work as gene models were too long with mis-formatted 5' and 3' ends
- some filtering based on amino acid composition and distribution via PCA
     - found that there were some models with long runs of biased amino acids and so was trying to correct for these. This did work, but did also remove some valid gene models that were naturally an edge case on distribution  
- post annotation removeal of repeats
- mapping of Iso-Seq data
     - this was used for various things, including in Apollo
- making a GO term database
     - this also worked, but didnt use it in the paper in the end
- manual curation in apollo
     - also worked, but done mostly post annotation freeze

```bash
#-----------------------------------------------------------------------------------------
# AUGUSTUS post Braker
#-----------------------------------------------------------------------------------------
#--- generate intron hints for augustus

# merge bams
ls -1 *qv.sorted.bam > bams.hq-lq.list
ls -1 *fasta.sorted.bam > bams.fl-nfl.list

bsub.py 5 bammerge1 "samtools-1.3 merge -b bams.hq-lq.list hq-lq.merged.bam"
bsub.py 5 bammerge2 "samtools-1.3 merge -b bams.fl-nfl.list fl-nfl.merged.bam"


#- make hints files
bsub.py 10 01_bam2hints_hq-lq \
"/lustre/scratch118/infgen/archive/ss34/SCHISTO/augustus/augustus-3.2.2//bin/bam2hints --intronsonly --minintronlen=20 --source=PB --in=hq-lq.merged.bam --out=hq-lq.merged.intron-hints.gff"


bsub.py 10 01_bam2hints_fl-nfl \
"/lustre/scratch118/infgen/archive/ss34/SCHISTO/augustus/augustus-3.2.2//bin/bam2hints --intronsonly --minintronlen=20 --source=PB --in=fl-nfl.merged.bam --out=fl-nfl.merged.intron-hints.gff"

samtools-1.3 faidx HAEM_V4_final.chr.fa
awk '{print $1,$2}' OFS="\t" HAEM_V4_final.chr.fa.fai >chromosome_size.list

samtools-1.3 index -b fl-nfl.merged.bam
/software/pathogen/external/apps/usr/local/Python-2.7.13/bin//bam2wig.py -i fl-nfl.merged.bam -s chromosome_size.list -o fl-nfl_bam2wig_out




cd /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/TRANSCRIPTOME/
mkdir AUGUSTUS_POST_BRAKER
cd AUGUSTUS_POST_BRAKER

ln -s ../STAR_MAP_CHR/RNAseq.merge.exon-parts.gff
ln -s ../STAR_MAP_CHR/starmap_merge.chr.intron-hints.gff
ln -s ../ISOSEQ_MAP_CHR/hq-lq.merged.intron-hints.gff
ln -s ../ISOSEQ_MAP_CHR/fl-nfl.merged.intron-hints.gff

cat *gff > hints.gff


# prepare reference
ln -s ../../REF/HAEM_V4_final.chr.fa

/software/pathogen/external/apps/usr/bin/splitMfasta.pl HAEM_V4_final.chr.fa --outputpath=./
for f in *.split.*; do NAME=`grep ">" $f`; mv $f ${NAME#>}.fa; done

/software/pathogen/external/apps/usr/bin/summarizeACGTcontent.pl HAEM_V4_final.chr.fa > basesummary.out

# prepare augustus run files

grep "bases" basesummary.out | awk -v PWD=$PWD -v HINTS=hints.gff '{print PWD"/"$3".fa",PWD"/"HINTS,"1",$1}' OFS="\t" > sequences.list


/software/pathogen/external/apps/usr/bin/createAugustusJoblist.pl --sequences=sequences.list --wrap="#" --overlap=50000 --chunksize=3000000 --outputdir=augustus_split_out --joblist=jobs.lst --jobprefix=augsplit --command \
"/nfs/users/nfs_s/sd21/lustre118_link/software/TRANSCRIPTOME/augustus_v3.3/bin/augustus --species=Hc_V4_chr --strand=both --genemodel=partial --protein=on --introns=on --start=on --stop=on --cds=on --codingseq=on --UTR=off --nc=off --gff3=on --alternatives-from-evidence=on --extrinsicCfgFile=/nfs/users/nfs_s/sd21/lustre118_link/software/TRANSCRIPTOME/augustus_v3.3/config/species/Hc_V4_chr/extrinsic.Hc_V4_chr.modified.cfg --AUGUSTUS_CONFIG_PATH=/nfs/users/nfs_s/sd21/lustre118_link/software/TRANSCRIPTOME/augustus_v3.3/config/"


mkdir augustus_split_out

# run augustus
for i in augsplit*; do echo -e "bsub.py 4 augsplit_log ./${i}" >> run_augsplit; done        #Max memory used was 1.2Gb
chmod a+x run_augsplit
bsub.py --queue yesterday 1 run_splitter ./run_augsplit

mkdir augsplit_run_dir
mv augsplit* augsplit_run_dir

cat augustus_split_out/*gff | /software/pathogen/external/apps/usr/bin/join_aug_pred.pl > HC_V4_augustus_merge.gff


echo -e "##gff-version 3" > HC_V4_augustus_merge.filtered.gff; grep "AUGUSTUS" HC_V4_augustus_merge.gff >> HC_V4_augustus_merge.filtered.gff


gffread HC_V4_augustus_merge.filtered.gff -g HAEM_V4_final.chr.fa -y HC_V4_augustus_merge.filtered.aa.fa
```


```bash
#-----------------------------------------------------------------------------------------
# Filtering post braker  AUGUSTUS output - MANUAL
#-----------------------------------------------------------------------------------------

cd /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/TRANSCRIPTOME/BRAKER_CHR/POST_BRAKER

ln -s ../braker/Hc_V4_chr/augustus.gff3
ln -s ../braker/Hc_V4_chr/augustus.aa
ln -s ../HAEM_V4_final.chr.fa



cut -f1 -d ' ' augustus.aa > augustus.aa2
~sd21/bash_scripts/AAResidueFreqCalculator.pl


grep -v "#" augustus.aa2.freq | head -n 31 | cut -f2- > augustus.aa2.freq2
```

```R
a<-read.table("augustus.aa2.freq2",header=T,sep="\t")
c<-prcomp(na.omit(t(a[1:ncol(a)])))
d<-as.data.frame(c$x[,1:30])


library(ggplot2)
library(gplots)

# Plot pairwise comparisons of PCs - do this iteratively to look for outlier clusters
ggplot()+geom_point(aes(d$PC2,d$PC3,alpha=0.2),size=0.5)+theme_bw()

# Plot rotations to look for drivers of PC variation - heatmap will short regression of factors with PCs, showing which part of the amino acid composition is driving partiular PCs
heatmap.2(c$rotation,Rowv=FALSE,Colv=FALSE,trace="none",dendrogram="none")



# apply filter based on PC distibutions - these need to be manually set based on looking at pairwise distributions of PCs and some blasting
filter<-d[(d$PC2>=-15 & d$PC2<=20 &
	d$PC3>=-25 & d$PC3<=20 &  
	d$PC4>=-25 & d$PC4<=25 &
	d$PC5>=-15 & d$PC5<=15 &
	d$PC6>=-20 & d$PC6<=20 &
	d$PC7>=-15 & d$PC7<=15 &
	d$PC8>=-15 & d$PC8<=15 &
	d$PC9>=-10 & d$PC9<=10 &
	d$PC10>=-10 & d$PC10<=12.5),]


filter<-d[(d$PC2>=-20 & d$PC2<=20 &
	d$PC3>=-30 & d$PC3<=30 &  
	d$PC4>=-25 & d$PC4<=35 &
	d$PC5>=-25 & d$PC5<=20 &
	d$PC6>=-25 & d$PC6<=30 &
	d$PC7>=-20 & d$PC7<=15 &
	d$PC8>=-20 & d$PC8<=15 &
	d$PC9>=-20 & d$PC9<=20 &
	d$PC10>=-15 & d$PC10<=12.5),]


# check filtered
ggplot()+geom_point(aes(filter$PC2,filter$PC3,alpha=0.2),size=0.5)+theme_bw()

write.table(row.names(filter),file="freq_filtered_transcripts.list",quote=F,row.names=F,col.names=F)
```


```bash
#--- filter transposon-like sequences

#- diamond filtered repeatmasker output - from V3, but will be fine
ln -s /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/V3/REPEATS/REPEATMASKER/consensi.fa.classified_SDfiltered.fa

# make fasta of CDS
gffread -g HAEM_V4_final.chr.fa -x HC_V4_augustus_merge.AAfiltered.CDS.fa HC_V4_augustus_merge.AAfiltered.gff

# make blast database
makeblastdb -in consensi.fa.classified_SDfiltered.fa -parse_seqids -dbtype nucl

# run blast
bsub.py --queue yesterday 10 01_blastn_repeats "/software/pubseq/bin/ncbi_blast+/blastn -outfmt 6 -db consensi.fa.classified_SDfiltered.fa  -query HC_V4_augustus_merge.AAfiltered.CDS.fa \> HC_V4_augustus_merge.AAfiltered.CDS.repeat_pos.out"

python3 ~alt/python/bin/blast_parser_for_filtering_retrotrans_by_length_of_cds.py HC_V4_augustus_merge.AAfiltered.CDS.repeat_pos.out HC_V4_augustus_merge.AAfiltered.CDS.fa > repeats_vs_cds.out

cat repeats_vs_cds.out |  awk '{if($4 >= 95){print$1}}' | sort > genes_to_remove.txt




cat freq_filtered_transcripts.list genes_to_remove.txt | sort -V | uniq -c | awk '{if($1=="1") print $2}' > filtered_transcripts_to_keep.list


while read id; do grep "${id}[;$]" augustus.gff3 >> filtered.augustus.gff.tmp; done < filtered_transcripts_to_keep.list
grep "mRNA" filtered.augustus.gff.tmp | cut -f2 -d ';' | sed 's/Parent=//g' | sort -V | uniq > gene_list.tmp

while read gene_id; do grep "ID=${gene_id};$" augustus.gff3; done < gene_list.tmp > genes.filtered.gff.tmp


echo "##gff-version 3" > HC_V4_augustus_merge.AAfiltered.gff; cat filtered.augustus.gff.tmp genes.filtered.gff.tmp | sort -k1,1 -k4,4n >> HC_V4_augustus_merge.AAfiltered.gff
while read name; do samtools-1.3 faidx HC_V4_augustus_merge.filtered.aa.fa2 ${name} >> HC_V4_augustus_merge.AAfiltered.aa; done < freq_filtered_transcripts.list



#--- post evm
while read id; do grep "${id}[;.]" HAEM_V4.chr.evm_merge.gff >> filtered.augustus.gff.tmp; done < filtered_transcripts_to_keep.list
grep "mRNA" filtered.augustus.gff.tmp | cut -f2 -d ';' | sed 's/Parent=//g' | sort -V | uniq > gene_list.tmp

while read gene_id; do grep "ID=${gene_id};" HAEM_V4.chr.evm_merge.gff >> genes.filtered.gff.tmp; done < gene_list.tmp
echo "##gff-version 3" > HC_V4_augustus_merge.AAfiltered.gff; cat filtered.augustus.gff.tmp genes.filtered.gff.tmp | sort -k1,1 -k4,4n >> HC_V4_augustus_merge.AAfiltered.gff

#-----------------------------------------------------------------------------------------
# Filtering post braker  AUGUSTUS output - MANUAL
#-----------------------------------------------------------------------------------------


mkdir POST_AUGUSTUS_FILTER
cd POST_AUGUSTUS_FILTER

ln -s ../HC_V4_augustus_merge.filtered.aa.fa
ln -s ../HC_V4_augustus_merge.filtered.gff

cut -f1 -d ' ' HC_V4_augustus_merge.filtered.aa.fa > HC_V4_augustus_merge.filtered.aa.fa2

 ~sd21/bash_scripts/AAResidueFreqCalculator.pl
 # input file = HC_V4_augustus_merge.filtered.aa.fa
# output file = HC_V4_augustus_merge.filtered.aa.freq

# Calculations
# 1	Length of amino acid
# 2	Percentage of non-polar aliphatic (GAVIL) residues
# 3	Percentage of non-polar aromatic (FW) residues
# 4	Percentage of non-polar cyclic (P) residues
# 5	Percentage of polar sulphur containing (CM) residues
# 6	Percentage of polar hydroxyl (ST) residues
# 7	Percentage of polar aromatic (Y) residues
# 8	Percentage of polar acidic-amide (NQ) residues
# 9	Percentage of acidic (DE) residues
# 10	Percentage of basic (RHK) residues
# 11	Percentage of alanine (A) residues
# 12	Percentage of cysteine (C) residues
# 13	Percentage of aspartate (D) residues
# 14	Percentage of glutamate (E) residues
# 15	Percentage of phenylalanine (F) residues
# 16	Percentage of glycine (G) residues
# 17	Percentage of histidine (H) residues
# 18	Percentage of isoleucine (I) residues
# 19	Percentage of lysine (K) residues
# 20	Percentage of leucine (L) residues
# 21	Percentage of methionine (M) residues
# 22	Percentage of asparagine (N) residues
# 23	Percentage of proline (P) residues
# 24	Percentage of glutamine (Q) residues
# 25	Percentage of arginine (R) residues
# 26	Percentage of serine (S) residues
# 27	Percentage of threonine (T) residues
# 28	Percentage of valine (V) residues
# 29	Percentage of tryptophan (W) residues
# 30	Percentage of tyrosine (Y) residues



grep -v "#" HC_V4_augustus_merge.filtered.aa.freq | head -n 31 | cut -f2- > HC_V4_augustus_merge.filtered.aa.freq2
```

```R
R
a<-read.table("HC_V4_augustus_merge.filtered.aa.freq2",header=T,sep="\t")
c<-prcomp(na.omit(t(a[1:ncol(a)])))
d<-as.data.frame(c$x[,1:30])


library(ggplot2)
library(gplots)


# Plot pairwise comparisons of PCs - do this iteratively to look for outlier clusters
ggplot()+geom_point(aes(d$PC2,d$PC3,alpha=0.2),size=0.5)+theme_bw()

# Plot rotations to look for drivers of PC variation - heatmap will short regression of factors with PCs, showing which part of the amino acid composition is driving partiular PCs
heatmap.2(c$rotation,Rowv=FALSE,Colv=FALSE,trace="none",dendrogram="none")



# apply filter based on PC distibutions - these need to be manually set based on looking at pairwise distributions of PCs and some blasting
filter<-d[(d$PC2>=-20 & d$PC2<=15 &
	d$PC3>=-25 & d$PC3<=30 &  
	d$PC4>=-20 & d$PC4<=20 &
	d$PC5>=-20 & d$PC5<=20 &
	d$PC6>=-15 & d$PC6<=20 &
	d$PC7>=-15 & d$PC7<=15 &
	d$PC8>=-12.5 & d$PC8<=15 &
	d$PC9>=-15 & d$PC9<=15 &
	d$PC10>=-15 & d$PC10<=15),]

write.table(row.names(filter),file="freq_filtered_transcripts.list",quote=F,row.names=F,col.names=F)
# exit R

while read id; do grep ${id} HC_V4_augustus_merge.filtered.gff >> filtered.augustus.gff.tmp; done < freq_filtered_transcripts.list
grep "transcript" filtered.augustus.gff.tmp | cut -f2 -d ';' | sed 's/Parent=//g' | sort -V | uniq > gene_list.tmp

while read gene_id; do grep "ID=${gene_id}$" HC_V4_augustus_merge.filtered.gff; done < gene_list.tmp > genes.filtered.gff.tmp


echo "##gff-version 3" > HC_V4_augustus_merge.AAfiltered.gff; cat filtered.augustus.gff.tmp genes.filtered.gff.tmp | sort -k1,1 -k4,4n >> HC_V4_augustus_merge.AAfiltered.gff
while read name; do samtools-1.3 faidx HC_V4_augustus_merge.filtered.aa.fa2 ${name} >> HC_V4_augustus_merge.AAfiltered.aa; done < freq_filtered_transcripts.list

```


## Post filtering of repeats

```bash

#--- ALans approach to removing repeats
# Alans instructions of Hmic

#My working is in /lustre/scratch118/infgen/team133/alt/HMIC/braker2/braker/HMN_v3
#
#
# Here are the steps:
#
# /software/pubseq/bin/ncbi_blast+/blastn -outfmt 6 -db HMN_V3-families_mobilome_only.fa -query HMN_v3_raw_cds.fa > HMN_V3-families_mobilome_only.fa_out
#
# python3 blast_parser_for_filtering_retrotrans_by_length_of_cds.py HMN_V3-families_mobilome_only.fa_out HMN_v3_raw_cds.fa > HMN_V3-families_mobilome_only.fa_vs_HMN_v3_raw_cds.fa_out
#
# #Some sort of filtering - histogram suggested cutoff of 97% blast hit:
#
# cat HMN_V3-families_mobilome_only.fa_vs_HMN_v3_raw_cds.fa_out | awk '{if($4 >= 97){print$1}}' | sort > genes_to_remove_hmnv3_rptlib.txt
#
# python3 gff3_parser.py augustus.gff3 genes_to_remove_combined_sorted.txt   #I used a combination of repeat libraries generated by 50HGI and by me using RepeatModeller to get the final list of transcripts to remove
#


#- diamond filtered repeatmasker output - from V3, but will be fine
ln -s /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/V3/REPEATS/REPEATMASKER/consensi.fa.classified_SDfiltered.fa

# make fasta of CDS
gffread -g ../HAEM_V4_final.chr.fa -x HC_V4_augustus_merge.AAfiltered.CDS.fa HC_V4_augustu_merge.AAfiltered.gff

# make balst database
makeblastdb -in consensi.fa.classified_SDfiltered.fa -parse_seqids -dbtype nucl

# run blast
bsub.py --queue yesterday 10 01_blastn_repeats "/software/pubseq/bin/ncbi_blast+/blastn -outfmt 6 -db consensi.fa.classified_SDfiltered.fa  -query HC_V4_augustus_merge.AAfiltered.CDS.fa \> HC_V4_augustus_merge.AAfiltered.CDS.repeat_pos.out"

python3 ~alt/python/bin/blast_parser_for_filtering_retrotrans_by_length_of_cds.py HC_V4_augustus_merge.AAfiltered.CDS.repeat_pos.out HC_V4_augustus_merge.AAfiltered.CDS.fa > repeats_vs_cds.out

cat repeats_vs_cds.out |  awk '{if($4 >= 50){print$1}}' | sort > genes_to_remove.txt

python3  /lustre/scratch118/infgen/team133/alt/HMIC/braker2/braker/HMN_v3/gff3_parser.py HC_V4_augustus_merge.AAfiltered.gff genes_to_remove.txt
```


```bash
#-----------------------------------------------------------------------------------------
# mapping isoseq data
#-----------------------------------------------------------------------------------------

# get data
cd /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/TRANSCRIPTOME/RAW/ISOSEQ

wget -O sequel_polished_high_qv_consensus_isoforms.fastq https://sf2-farm-srv1.internal.sanger.ac.uk:8243/SMRTLink/1.0.0/secondary-analysis/datastore-files/3e563bad-4a17-47f9-b9e0-f61fe09c5640/download
fastaq to_fasta sequel_polished_high_qv_consensus_isoforms.fastq sequel_polished_high_qv_consensus_isoforms.fasta

wget -O sequel_polished_low_qv_consensus_isoforms.fastq https://sf2-farm-srv1.internal.sanger.ac.uk:8243/SMRTLink/1.0.0/secondary-analysis/datastore-files/3be12200-50d7-4167-9ba3-def01292c10c/download
fastaq to_fasta sequel_polished_low_qv_consensus_isoforms.fastq sequel_polished_low_qv_consensus_isoforms.fasta

wget -O RSII_polished_high_qv_consensus_isoforms.fasta http://sf2-farm-srv2.internal.sanger.ac.uk:8080/smrtportal/api/jobs/20636/contents/data/polished_high_qv_consensus_isoforms.fasta
wget -O RSII_polished_low_qv_consensus_isoforms.fasta http://sf2-farm-srv2.internal.sanger.ac.uk:8080/smrtportal/api/jobs/20636/contents/data/polished_low_qv_consensus_isoforms.fasta

cat *consensus_isoforms.fasta > all_isoseq.fasta

#--- mapping isoseq data
cd /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/TRANSCRIPTOME/ISOSEQ_MAP
mkdir ISOSEQ_MAP_ALL
mkdir ISOSEQ_MAP_CHR

cd ISOSEQ_MAP_ALL
ln -s /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/REF/HAEM_V4_final.all.fa
for i in `(cd ../RAW/ISOSEQ/ ; ls -1 *fasta)`; do \
bsub.py --threads 3 10 01_minimap2_isoseq "~sd21/bash_scripts/run_minimap2_splice ../RAW/ISOSEQ/${i} HAEM_V4_final.all.fa ${i}"; \
done

cd ../
cd ISOSEQ_MAP_CHR
ln -s /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/REF/HAEM_V4_final.chr.fa
for i in `(cd ../RAW/ISOSEQ/ ; ls -1 *fasta)`; do \
bsub.py --threads 3 10 01_minimap2_isoseq "~sd21/bash_scripts/run_minimap2_splice ../RAW/ISOSEQ/${i} HAEM_V4_final.chr.fa ${i}"; \
done



# TO CHECK AND FIX
#https://github.com/Magdoll/cDNA_Cupcake/wiki/Cupcake-ToFU:-supporting-scripts-for-Iso-Seq-after-clustering-step

#gmap -D /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/TRANSCRIPTOME/ISOSEQ -d haem_v3.3.chr -f samse -n 0 -t 12 -z sense_force all.HQ.isoforms.fastq > hq_isoforms.fastq.sam
#sort -k 3,3 -k 4,4n hq_isoforms.fastq.sam > hq_isoforms.fastq.sorted.sam
#source activate anaCogent
#bsub.py 10 03_collapse_isoforms collapse_isoforms_by_sam.py --input all.HQ.isoforms.fastq --fq -s hq_isoforms.fastq.sorted.sam --dun-merge-5-shorter -o test
#bsub.py 10 04_fusion_isoforms fusion_finder.py --input all.HQ.isoforms.fastq --fq -s hq_isoforms.fastq.sorted.sam -o lq_isoforms.fasta.fusion --cluster_report_csv cluster_report.csv

export PATH="/nfs/users/nfs_s/sd21/lustre118_link/software/anaconda2/bin:$PATH"
source activate anaCogent

cat sequel_polished_high_qv.sam RSII_polished_high_qv.sam | sort -k 3,3 -k 4,4n > polished_high_qv.sorted.sam

cat ../RAW/ISOSEQ/sequel_polished_high_qv_consensus_isoforms.fasta ../RAW/ISOSEQ/RSII_polished_high_qv_consensus_isoforms.fasta > high_qv_consensus_isoforms.fasta
bsub.py 1 03_collapse_isoforms collapse_isoforms_by_sam.py --input high_qv_consensus_isoforms.fasta -s polished_high_qv.sorted.sam --dun-merge-5-shorter -o high_qv_consensus

bsub.py 10 04_fusion_isoforms fusion_finder.py --input high_qv_consensus_isoforms.fasta -s polished_high_qv.sorted.sam -o polished_high_qv.fusion --cluster_report_csv polished_high_qv_cluster_report.csv
```




### Making GO term databases
- lifting C elegans GO terms to 1:1 orthologs
- lifting over many:1 or 1:many, providing all GO terms for the many are all the same
- extracting interpro GO terms from IPS annotation in gff
- sorting / filtering to get a unique set.




Get 1:1s
```shell

cd /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/GO_ANALYSIS

ln -sf /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/SELECTION/PROTEIN_FASTAs/Results_Jan25/Orthologues_Jan25/Orthologues/Orthologues_ce.proteins.unique/ce.proteins.unique__v__hc_V4.proteins.unique.csv
ln -sf /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/SELECTION/PROTEIN_FASTAs/Results_Jan25/Orthologues_Jan25/Or
thologues/Orthologues_hc_V4.proteins.unique/hc_V4.proteins.unique__v__ce.proteins.unique.csv

ln -sf /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/
ln -sf /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/TRANSCRIPTOME/TRANSCRIPTOME_CURATION/HCON_V4_WBP11plus_1901
25.ips.gff3

# from local computer - transfer downloaded Ce GO terms from WBP
scp 02_data/WBP_Ce_GOterms_download.txt sd21@pcs5.internal.sanger.ac.uk:/nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/GO_ANALYSIS/

# extract 1:1s
awk 'NF==3 {print $2,$3}' OFS="\t" ce.proteins.unique__v__hc_V4.proteins.unique.csv | sed -e 's/gene=//g' -e 's/\-.*//g' -e 's/\..*//g' | grep -v "unique" > ce_hcV4_1to1.list
#--- 7361 ce_hc4_1to1.list


# extract GO terms from mRNAs in GFF
awk '$3=="mRNA" {print $0}' OFS="\t" HCON_V4_WBP11plus_190125.ips.gff3 | grep "Ontology_term" | cut -f 9 | cut -f3,4 -d ";" | sed -e 's/;/\t/g' -e 's/Name=//g' -e 's/"//g' -e 's/Ontology_term=//g' -e 's/-.*\t/\t/g' -e 's/,/\t/g' -e 's/\..*\t/\t/g' > annotation_GO_per_gene.txt

# work through columns to split multiple go terms per gene into one term per gene, repeating the gene name if multiple GO terms present
awk '{for(i=2; i<=NF; i++) {print $1,$i}}' OFS="\t" annotation_GO_per_gene.txt > annotation_GO_per_gene_split.txt

# convert Ce terms into Hc terms.
# ce genes with GO terms: 14638
# total ce GO terms: 195609

cut -f2,3 WBP_Ce_GOterms_download.txt | grep "GO:" > Ce_genes_GOterms.txt

# run the real substitution
export PATH="/nfs/users/nfs_s/sd21/lustre118_link/software/anaconda2/bin:$PATH"
# fsed - https://github.com/wroberts/fsed

fsed --pattern-format=tsv --output Ce_genes_GOterms_Hc_sub.txt  ce_hcV4_1to1.list Ce_genes_GOterms.txt &

grep "HCON" Ce_genes_GOterms_Hc_sub.txt > Hc_genes_Ce_GOterms.txt

# bring it all together and sorted
cat Hc_genes_Ce_GOterms.txt annotation_GO_per_gene_split.txt | sort | uniq -c

cat Hc_genes_Ce_GOterms.txt annotation_GO_per_gene_split.txt | sort | uniq > HCON_V4_GOterm.db

# Before liftover - original IPS annotation
#--- total GO terms: 20418
#--- total genes with GO term: 7824

# After liftover
#--- total GO terms: 60001
#--- total genes with GO term: 9627
```



# find GO terms for multi Ce: single Hc genes, for which all GO terms in the multi Ce are conserved.
cat hc_V4.proteins.unique__v__ce.proteins.unique.csv | cut -f2 | awk '{print NF}' OFS="\t"  > hc.fields
cat hc_V4.proteins.unique__v__ce.proteins.unique.csv | cut -f3 | awk '{print NF}' OFS="\t"   > ce.fields

cat ce.proteins.unique__v__hc_V4.proteins.unique.csv | awk -F '[\t]'  '{print $2}' | sed -e 's/gene=//g' -e 's/,//g' > ce.genes
paste ce.fields hc.fields ce.genes  | awk '$1>1 && $2==1 {print $0}' > hc_1toMulti.txt
# number of Hc genes: 616

cat hc_V4.proteins.unique__v__ce.proteins.unique.csv | cut -f2 > hc.genes
paste hc.fields ce.fields hc.genes ce.genes > hc.ce.data

awk '$1==1 && $2>1 {print $0}' hc.ce.data | sed 's/-.*\t/\t/g' | cut -f3,4 | sort -k3 | uniq > hc.ce.data.filtered

>hcgenes.cemultiGOs
while read hgene cgene; do
     echo -e "$hgene $cgene" | awk '{for(i=2; i<=NF; i++) {print $i}}' > ${hgene}.tmp
     count=$(wc -l ${hgene}.tmp | cut -f1 -d " ")
     grep -f ${hgene}.tmp WBP_Ce_GOterms_download.txt | cut -f 3 | sort | uniq -c | awk -v count="${count}" -v hgene="${hgene}" '{if($1>=count && $2!="") print hgene,$2}' OFS="\t" >> hcgenes.cemultiGOs;
     done < hc.ce.data.filtered

rm *tmp*

cat HCON_V4_GOterm.db hcgenes.cemultiGOs | sort | uniq > tmp; mv tmp HCON_V4_GOterm.db

# multiCe genes : single Hc genes
#--- additional GO terms: 3694
#--- additional genes with GO term: 517



#Total
#--- genes with GO terms: 9739
#--- GO terms: 62733
```




## 03 - Manual Curation in Apollo <a name="manual_curation_apollo"></a>

The genome annotaiton has been maually curated in apollo.

Tracks used
- Genome annotaiton
- RNAseq per lifestage
- splice leaders (SL1 & SL2)
- Pacbio IsoSeq (CCS subreads and HQ isoforms)



Once out of Apollo, some curation needs to be done to clean things up a little. The main problem is that Apollo uses a unique code ID per feature to keep track of informaiton and to make sure there are no clashes in IDs. While this is important in Apollo, it makes it confusing in downstream analyses that use the annotaiton. Decided to replace these so they are consistent throughout the whole annotaiton.

*NOTE* this approach below will have to be modified for subsequent apollo updates to ensure consistent naming of features. Eg, if a new isoform is added.

### Working environment
```shell
cd /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/TRANSCRIPTOME/TRANSCRIPTOME_CURATION

```
Get dumped GFF form local computer
```shell
scp Desktop/hc_v4.gff3.xz sd21@pcs5.internal.sanger.ac.uk:/nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/TRANSCRIPTOME/TRANSCRIPTOME_CURATION/
```

```shell
#unzip it - it is in xz format
unxz hc_v4.gff3.xz

# make a copy to work on
mv hc_v4.gff3 HCON_V4_WBP11plus_190125.gff3


# get mRNAs IDs and NAMES
awk '$3=="mRNA" {print $0}' OFS="\t"  HCON_V4_WBP11plus_190125.gff3 | sed -e 's/Note=Manually dissociate transcript from gene;//g' | cut -f3,5 -d ";" | sed -e 's/ID=//g' -e 's/;Name=/\t/g' > mRNA_IDs_NAMEs.txt

# remove transcript extensions
cat mRNA_IDs_NAMEs.txt | awk  '{ $2 = substr($2, 1,13); print }' OFS="\t" | sort -k2 > mRNA_IDs_NAMEs_trimmed-sorted.txt

# make a unique list
cut -f2  mRNA_IDs_NAMEs_trimmed-sorted.txt | uniq > mRNA_IDs_NAMEs_unique.txt

# add unique ID to each transcript
>mRNA_IDs_NAMEs_transcriptIDs.txt
while read NAME; do grep -w ${NAME} mRNA_IDs_NAMEs_trimmed-sorted.txt | cat -n | awk '{print $2,$3,$3"-0000"$1}' OFS="\t" >> mRNA_IDs_NAMEs_transcriptIDs.txt; done < mRNA_IDs_NAMEs_unique.txt

# run the real substitution
export PATH="/nfs/users/nfs_s/sd21/lustre118_link/software/anaconda2/bin:$PATH"
# fsed - https://github.com/wroberts/fsed
awk '{print $1,$3}' OFS="\t" mRNA_IDs_NAMEs_transcriptIDs.txt > mRNA_IDs_NAMEs_transcriptIDs.2.txt

fsed --pattern-format=tsv --output HCON_V4_WBP11plus_190125.renamed.gff3  mRNA_IDs_NAMEs_transcriptIDs.2.txt HCON_V4_WBP11plus_190125.gff3 &

```


Fixing GFF to prepare for interproscan. Stripping out info from existing interproscan, as is is incorrectly formatted.

```shell
sed -e 's/;info.*$//g' -e 's/method.*//g' -e '/^$/d' HCON_V4_WBP11plus_190125.renamed.gff3 > tmp.gff

echo '##FASTA' >> tmp.gff
ln -sf ../../REF/HAEM_V4_final.chr.fa
cat tmp.gff HAEM_V4_final.chr.fa > tmp.gff2; mv tmp.gff2 tmp.gff
```

interproscan
```shell
# generate a protein fasta from annotation and reference genome
gffread -y PROTEINS.fa -g HAEM_V4_final.chr.fa HCON_V4_WBP11plus_190125.renamed.gff3
sed -e 's/\.//g' PROTEINS.fa > tmp; mv tmp PROTEINS.fa

# run interproscan
farm_interproscan -a PROTEINS.fa -o IPS.output.gff

# lift over GO terms from interproscan to GFF
extract_interproscan_go_terms -i IPS.output.gff -e HCON_V4_WBP11plus_190125.renamed.gff3

# filter and rename
grep ^'\#\#\|hc' tmp.gff.go.gff | grep -v "mtDNA" | grep -v "FASTA" | grep -v ">" > HCON_V4_WBP11plus_190125.ips.gff3


```
