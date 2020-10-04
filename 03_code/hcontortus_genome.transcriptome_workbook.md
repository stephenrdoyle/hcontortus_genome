
---
## 04 - Kallisto <a name="kallisto"></a>
---


Date 190116

Using Kallisto to quantify transcripts, and use in the gene expression clustering analyses.


### Working environment
```shell
cd /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/TRANSCRIPTOME
mkdir KALLISTO
cd KALLISTO
```


### Get the GFF to work on
```shell
ln -fs /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/REF/HAEM_V4_final.chr.fa REF.fa
ln -fs /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/TRANSCRIPTOME/TRANSCRIPTOME_CURATION/HCON_V4_WBP11plus_190125.ips.gff3 ANNOTATION.gff3
```

### get some raw data
```shell
cd ~/lustre118_link/hc/GENOME/TRANSCRIPTOME/RAW

kinit # log into iRODs
icd /seq/7059
ils | grep "7059_6" | grep -v "phi" | while read -r name; do iget /seq/7059/${name} . ; done
ils | grep "7062_6" | grep -v "phi" | while read -r name; do iget /seq/7062/${name} . ; done
rm *168.ba* 7059_6#0.bam 7062_6#0.bam

# convert bam to fastqs
for i in *.bam; do samtools fastq -1 ${i%.bam}_1.fastq.gz -2 ${i%.bam}_2.fastq.gz ${i}; done
rm *.bam
rename "s/#/_/g" *
```




### Run Kallisto
```shell
cd /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/TRANSCRIPTOME/KALLISTO

# make a transcripts fasta
gffread -x TRANSCRIPTS.fa -g REF.fa ANNOTATION.gff3


# index the transcripts
kallisto index --index HCON_V4.TRANSCRIPTS.ixd TRANSCRIPTS.fa

# run kallisto
for i in ` cd ../RAW/ ; ls -1 *_1.fastq.gz | sed -e "s/_1.fastq.gz//g" `; do \
kallisto quant \
--bias \
--index HCON_V4.TRANSCRIPTS.ixd \
--output-dir kallisto_${i}_out \
--bootstrap-samples 100 \
--threads 7 \
--fusion \
../RAW/${i}_1.fastq.gz ../RAW/${i}_2.fastq.gz; done

mkdir KALLISTO_MAPPED_SAMPLES
mv kallisto_* KALLISTO_MAPPED_SAMPLES/

```

To run sleuth, a metadata file is needed with all samples IDs, conditions, and paths.





### Load R and environment

```R

R-3.5.0
#load(file = "hcontortus_genome.workbook.Rdata")
library("sleuth")
library(ggplot2)
library(patchwork)
```

### Run Sleuth
```R
hc_metadata <- read.table("sample_name_path.list", header = TRUE, stringsAsFactors=FALSE)
hc_so <- sleuth_prep(hc_metadata, extra_bootstrap_summary = TRUE)
hc_so <- sleuth_fit(hc_so, ~name, 'full')
hc_so <- sleuth_fit(hc_so, ~1, 'reduced')
hc_so <- sleuth_lrt(hc_so, 'reduced', 'full')

sleuth_table <- sleuth_results(hc_so, 'reduced:full', 'lrt', show_all = FALSE)
sleuth_significant <- dplyr::filter(sleuth_table, qval <= 0.05)
head(sleuth_significant, 20)
```

Generate some plots for QC
```R
pcaplot_allsamples	<-	plot_pca(hc_so, color_by = 'name')
heatmap_allsamples   <-    plot_sample_heatmap(hc_so)


# PCA of all samples shows gut a both variable and outlier - will remove just to have a look
hc_so_noGUT_meta   <-    hc_metadata[(hc_metadata$name!="GUT"),]
hc_so_noGUT   <-    sleuth_prep(hc_so_noGUT_meta, extra_bootstrap_summary = TRUE)
hc_so_noGUT   <-    sleuth_fit(hc_so_noGUT, ~name, 'full')
hc_so_noGUT   <-    sleuth_fit(hc_so_noGUT, ~1, 'reduced')
hc_so_noGUT   <-    sleuth_lrt(hc_so_noGUT, 'reduced', 'full')

pcaplot_allsamples_minusgut	<-	plot_pca(hc_so_noGUT, color_by = 'name')



# PCA of L3 samples - SHL3 and EXL3
hc_so_L3_meta    <-    hc_metadata[(hc_metadata$name=="EXL3" | hc_metadata$name=="SHL3"),]
hc_so_L3    <-    sleuth_prep(hc_so_L3_meta, extra_bootstrap_summary = TRUE)
hc_so_L3    <-    sleuth_fit(hc_so_L3, ~name, 'full')
hc_so_L3    <-    sleuth_fit(hc_so_L3, ~1, 'reduced')
hc_so_L3    <-    sleuth_lrt(hc_so_L3, 'reduced', 'full')


pcaplot_L3	<-	plot_pca(hc_so_L3, color_by = 'name')
pc_varianceplot_L3    <-    plot_pc_variance(hc_so_L3)
heatmap_L3   <-    plot_sample_heatmap(hc_so_L3)

# patchwork
kallistoQC_L3_plots <- (pcaplot_L3 | pc_varianceplot_L3) / heatmap_L3 + plot_layout(ncol = 1)
ggsave("kallistoQC_L3_plots.pdf",width = 28, height = 28, units = "cm")
ggsave("kallistoQC_L3_plots.png",width = 28, height = 28, units = "cm")

```
- Copy to local dir - run this from local machine
```shell
scp sd21@pcs5.internal.sanger.ac.uk:/nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/TRANSCRIPTOME/KALLISTO/kallistoQC_L3_plots.* ~/Documents/workbook/hcontortus_genome/04_analysis
```

![Kallisto QC - L3 plots](04_analysis/kallistoQC_L3_plots.png)
Fig - Kalliso QC - L3 QC plots (i) PCA (ii) Loading plot of each PC, (iii) Heatmap


Because the the L3 samples have been mixed, need to regenerate the metadata file to reflect the switch of L3 IDs. The new IDs are as follows:
- SHL3
    - 7062_6_10
    - 7062_6_8
    - 7062_6_11
    - 7062_6_7 (probable - may drop)
- EXL3
    - 7062_6_12
    - 7062_6_9

Made a new file called "sample_name_path_L3fixed.list" with correct IDs, paths


```R
R-3.5.0
#load(file = "hcontortus_genome.workbook.Rdata")
library("sleuth")
library(ggplot2)
library(patchwork)

hc_metadata_L3fixed   <-    read.table("sample_name_path_L3fixed.list", header = TRUE, stringsAsFactors=FALSE)
hc_so   <-    sleuth_prep(hc_metadata_L3fixed, extra_bootstrap_summary = TRUE)
hc_so   <-    sleuth_fit(hc_so, ~name, 'full')
hc_so   <-    sleuth_fit(hc_so, ~1, 'reduced')
hc_so   <-    sleuth_lrt(hc_so, 'reduced', 'full')

sleuth_table    <-    sleuth_results(hc_so, 'reduced:full', 'lrt', show_all = FALSE)
sleuth_significant    <-    dplyr::filter(sleuth_table, qval <= 0.05)
head(sleuth_significant, 20)

pcaplot_allsamples2	<-	plot_pca(hc_so, color_by = 'name')
heatmap_allsamples2   <-    plot_sample_heatmap(hc_so)
kallistoQC_allsamples2_plots <- pcaplot_allsamples2 + heatmap_allsamples2 + plot_layout(ncol = 2)
ggsave("kallistoQC_allsamples2_plots.pdf",width = 28, height = 10, units = "cm")
ggsave("kallistoQC_allsamples2_plots.png",width = 28, height = 10, units = "cm")
```
- Copy to local dir - run this from local machine
```shell
scp sd21@pcs5.internal.sanger.ac.uk:/nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/TRANSCRIPTOME/KALLISTO/kallistoQC_allsamples2_plots.* ~/Documents/workbook/hcontortus_genome/04_analysis
```

![Kallisto QC - All samples with L3 fixed](04_analysis/kallistoQC_allsamples2_plots.png)
Fig - Kalliso QC - All samples with L3 IDs fixed


### Run pairwise comparisons
Want to generate tables of most significantly DE genes per pair, comparing sensible transitions throughout the life cycle

```R
# EGG vs L1 only
hc_so_EGGvL1	<-	hc_metadata_L3fixed[(hc_metadata_L3fixed$name=="EGG" | hc_metadata_L3fixed$name=="L1"),]
hc_so_EGGvL1 <- sleuth_prep(hc_so_EGGvL1, extra_bootstrap_summary = TRUE,num_cores=1)
hc_so_EGGvL1 <- sleuth_fit(hc_so_EGGvL1, ~name, 'full')
hc_so_EGGvL1 <- sleuth_fit(hc_so_EGGvL1, ~1, 'reduced')
hc_so_EGGvL1 <- sleuth_lrt(hc_so_EGGvL1, 'reduced', 'full')

#sleuth_table_EGGvL1 <- sleuth_results(hc_so_EGGvL1, 'reduced:full', 'lrt', show_all = FALSE)
#sleuth_significant_EGGvL1 <- dplyr::filter(sleuth_table_EGGvL1, qval <= 0.05)
#head(sleuth_significant_EGGvL1, 20)
#write.table(sleuth_table_EGGvL1,file="sleuth_table_EGGvL1.txt",sep="\t",quote=FALSE, row.names=FALSE)

#sleuth_live(hc_so_EGGvL1)

hc_so_EGGvL1_wt<-sleuth_wt(hc_so_EGGvL1,'nameL1',which_model = "full")
sleuth_table_EGGvL1_wt <- sleuth_results(hc_so_EGGvL1_wt,test="nameL1",which_model = "full",test_type = 'wt')

# up
sleuth_significant_EGGvL1_wt_up <- dplyr::filter(sleuth_table_EGGvL1_wt, qval <= 0.01,b>=2)
write.table(sleuth_significant_EGGvL1_wt_up,file="sleuth_table_EGGvL1_up.txt",sep="\t",quote=FALSE, row.names=FALSE)
# down
sleuth_significant_EGGvL1_wt_down <- dplyr::filter(sleuth_table_EGGvL1_wt, qval <= 0.01,b<=-2)
write.table(sleuth_significant_EGGvL1_wt_down,file="sleuth_table_EGGvL1_down.txt",sep="\t",quote=FALSE, row.names=FALSE)



# L1 vs SHL3
hc_so_L1vSHL3_meta	<-	hc_metadata_L3fixed[(hc_metadata_L3fixed$name=="L1" | hc_metadata_L3fixed$name=="SHL3"),]
hc_so_L1vSHL3 <- sleuth_prep(hc_so_L1vSHL3_meta, extra_bootstrap_summary = TRUE,num_cores=1)
hc_so_L1vSHL3 <- sleuth_fit(hc_so_L1vSHL3, ~name, 'full')
hc_so_L1vSHL3 <- sleuth_fit(hc_so_L1vSHL3, ~1, 'reduced')
hc_so_L1vSHL3 <- sleuth_lrt(hc_so_L1vSHL3, 'reduced', 'full')

sleuth_table_L1vSHL3 <- sleuth_results(hc_so_L1vSHL3, 'reduced:full', 'lrt', show_all = FALSE)
#sleuth_significant_L1vSHL3 <- dplyr::filter(sleuth_table_L1vSHL3, qval <= 0.05)
#head(sleuth_significant_L1vSHL3, 20)
#write.table(sleuth_table_L1vSHL3,file="sleuth_table_L1vSHL3.txt",sep="\t",quote=FALSE, row.names=FALSE)


hc_so_L1vSHL3_wt<-sleuth_wt(hc_so_L1vSHL3,'nameSHL3',which_model = "full")
sleuth_table_L1vSHL3_wt <- sleuth_results(hc_so_L1vSHL3_wt,test="nameSHL3",which_model = "full",test_type = 'wt')

# up
sleuth_significant_L1vSHL3_wt_up <- dplyr::filter(sleuth_table_L1vSHL3_wt, qval <= 0.01,b>=2)
write.table(sleuth_significant_L1vSHL3_wt_up,file="sleuth_table_L1vSHL3_up.txt",sep="\t",quote=FALSE, row.names=FALSE)
# down
sleuth_significant_L1vSHL3_wt_down <- dplyr::filter(sleuth_table_L1vSHL3_wt, qval <= 0.01,b<=-2)
write.table(sleuth_significant_L1vSHL3_wt_down,file="sleuth_table_L1vSHL3_down.txt",sep="\t",quote=FALSE, row.names=FALSE)







# SHL3 vs EXL3
hc_so_SHL3vEXL3	<-	hc_metadata_L3fixed[(hc_metadata_L3fixed$name=="SHL3" | hc_metadata_L3fixed$name=="EXL3"),]
hc_so_SHL3vEXL3 <- sleuth_prep(hc_so_SHL3vEXL3, extra_bootstrap_summary = TRUE,num_cores=1)
hc_so_SHL3vEXL3 <- sleuth_fit(hc_so_SHL3vEXL3, ~name, 'full')
hc_so_SHL3vEXL3 <- sleuth_fit(hc_so_SHL3vEXL3, ~1, 'reduced')
hc_so_SHL3vEXL3 <- sleuth_lrt(hc_so_SHL3vEXL3, 'reduced', 'full')

sleuth_table_SHL3vEXL3 <- sleuth_results(hc_so_SHL3vEXL3, 'reduced:full', 'lrt', show_all = FALSE)
#sleuth_significant_SHL3vEXL3 <- dplyr::filter(sleuth_table_SHL3vEXL3, qval <= 0.05)
#head(sleuth_significant_SHL3vEXL3, 20)


hc_so_SHL3vEXL3_wt<-sleuth_wt(hc_so_SHL3vEXL3,'nameSHL3',which_model = "full")
sleuth_table_SHL3vEXL3_wt <- sleuth_results(hc_so_SHL3vEXL3_wt,test="nameSHL3",which_model = "full",test_type = 'wt')

# up
sleuth_significant_SHL3vEXL3_wt_up <- dplyr::filter(sleuth_table_SHL3vEXL3_wt, qval <= 0.01,b>=2)
write.table(sleuth_significant_SHL3vEXL3_wt_up,file="sleuth_table_SHL3vEXL3_up.txt",sep="\t",quote=FALSE, row.names=FALSE)
# down
sleuth_significant_SHL3vEXL3_wt_down <- dplyr::filter(sleuth_table_SHL3vEXL3_wt, qval <= 0.01,b<=-2)
write.table(sleuth_significant_SHL3vEXL3_wt_down,file="sleuth_table_SHL3vEXL3_down.txt",sep="\t",quote=FALSE, row.names=FALSE)





# EXL3 vs L4
hc_so_EXL3vL4	<-	hc_metadata_L3fixed[(hc_metadata_L3fixed$name=="EXL3" | hc_metadata_L3fixed$name=="L4"),]
hc_so_EXL3vL4 <- sleuth_prep(hc_so_EXL3vL4, extra_bootstrap_summary = TRUE,num_cores=1)
hc_so_EXL3vL4 <- sleuth_fit(hc_so_EXL3vL4, ~name, 'full')
hc_so_EXL3vL4 <- sleuth_fit(hc_so_EXL3vL4, ~1, 'reduced')
hc_so_EXL3vL4 <- sleuth_lrt(hc_so_EXL3vL4, 'reduced', 'full')

#sleuth_table_EXL3vL4 <- sleuth_results(hc_so_EXL3vL4, 'reduced:full', 'lrt', show_all = FALSE)
#sleuth_significant_EXL3vL4 <- dplyr::filter(sleuth_table_EXL3vL4, qval <= 0.05)
#head(sleuth_significant_EXL3vL4, 20)
#write.table(sleuth_table_EXL3vL4,file="sleuth_table_EXL3vL4.txt",sep="\t",quote=FALSE, row.names=FALSE)


hc_so_EXL3vL4_wt<-sleuth_wt(hc_so_EXL3vL4,'nameL4',which_model = "full")
sleuth_table_EXL3vL4_wt <- sleuth_results(hc_so_EXL3vL4_wt,test="nameL4",which_model = "full",test_type = 'wt')

# up
sleuth_significant_EXL3vL4_wt_up <- dplyr::filter(sleuth_table_EXL3vL4_wt, qval <= 0.01,b>=2)
write.table(sleuth_significant_EXL3vL4_wt_up,file="sleuth_table_EXL3vL4_up.txt",sep="\t",quote=FALSE, row.names=FALSE)
# down
sleuth_significant_EXL3vL4_wt_down <- dplyr::filter(sleuth_table_EXL3vL4_wt, qval <= 0.01,b<=-2)
write.table(sleuth_significant_EXL3vL4_wt_down,file="sleuth_table_EXL3vL4_down.txt",sep="\t",quote=FALSE, row.names=FALSE)






# L4 vs Adult Male
hc_so_L4vADULTM	<-	hc_metadata_L3fixed[(hc_metadata_L3fixed$name=="L4" | hc_metadata_L3fixed$name=="ADULT_M"),]
hc_so_L4vADULTM <- sleuth_prep(hc_so_L4vADULTM, extra_bootstrap_summary = TRUE,num_cores=1)
hc_so_L4vADULTM <- sleuth_fit(hc_so_L4vADULTM, ~name, 'full')
hc_so_L4vADULTM <- sleuth_fit(hc_so_L4vADULTM, ~1, 'reduced')
hc_so_L4vADULTM <- sleuth_lrt(hc_so_L4vADULTM, 'reduced', 'full')

#sleuth_table_L4vADULTM <- sleuth_results(hc_so_L4vADULTM, 'reduced:full', 'lrt', show_all = FALSE)
#sleuth_significant_L4vADULTM <- dplyr::filter(sleuth_table_L4vADULTM, qval <= 0.05)
#head(sleuth_significant_L4vADULTM, 20)
#write.table(sleuth_table_L4vADULTM,file="sleuth_table_L4vADULTM.txt",sep="\t",quote=FALSE, row.names=FALSE)


hc_so_L4vADULTM_wt<-sleuth_wt(hc_so_L4vADULTM,'nameL4',which_model = "full")
sleuth_table_L4vADULTM_wt <- sleuth_results(hc_so_L4vADULTM_wt,test="nameL4",which_model = "full",test_type = 'wt')


# up
sleuth_significant_L4vADULTM_wt_up <- dplyr::filter(sleuth_table_L4vADULTM_wt, qval <= 0.01,b>=2)
write.table(sleuth_significant_L4vADULTM_wt_up,file="sleuth_table_L4vADULTM_up.txt",sep="\t",quote=FALSE, row.names=FALSE)
# down
sleuth_significant_L4vADULTM_wt_down <- dplyr::filter(sleuth_table_L4vADULTM_wt, qval <= 0.01,b<=-2)
write.table(sleuth_significant_L4vADULTM_wt_down,file="sleuth_table_L4vADULTM_down.txt",sep="\t",quote=FALSE, row.names=FALSE)





# L4 vs Adult Female
hc_so_L4vADULTF	<-	hc_metadata_L3fixed[(hc_metadata_L3fixed$name=="L4" | hc_metadata_L3fixed$name=="ADULT_F"),]
hc_so_L4vADULTF <- sleuth_prep(hc_so_L4vADULTF, extra_bootstrap_summary = TRUE,num_cores=1)
hc_so_L4vADULTF <- sleuth_fit(hc_so_L4vADULTF, ~name, 'full')
hc_so_L4vADULTF <- sleuth_fit(hc_so_L4vADULTF, ~1, 'reduced')
hc_so_L4vADULTF <- sleuth_lrt(hc_so_L4vADULTF, 'reduced', 'full')

#sleuth_table_L4vADULTF <- sleuth_results(hc_so_L4vADULTF, 'reduced:full', 'lrt', show_all = FALSE)
#sleuth_significant_L4vADULTF <- dplyr::filter(sleuth_table_L4vADULTF, qval <= 0.05)
#head(sleuth_significant_L4vADULTF, 20)
#write.table(sleuth_table_L4vADULTF,file="sleuth_table_L4vADULTF.txt",sep="\t",quote=FALSE, row.names=FALSE)


hc_so_L4vADULTF_wt<-sleuth_wt(hc_so_L4vADULTF,'nameL4',which_model = "full")
sleuth_table_L4vADULTF_wt <- sleuth_results(hc_so_L4vADULTF_wt,test="nameL4",which_model = "full",test_type = 'wt')

# up
sleuth_significant_L4vADULTF_wt_up <- dplyr::filter(sleuth_table_L4vADULTF_wt, qval <= 0.01,b>=2)
write.table(sleuth_significant_L4vADULTF_wt_up,file="sleuth_table_L4vADULTF_up.txt",sep="\t",quote=FALSE, row.names=FALSE)
# down
sleuth_significant_L4vADULTF_wt_down <- dplyr::filter(sleuth_table_L4vADULTF_wt, qval <= 0.01,b<=-2)
write.table(sleuth_significant_L4vADULTF_wt_down,file="sleuth_table_L4vADULTF_down.txt",sep="\t",quote=FALSE, row.names=FALSE)








# Adult Male vs Adult Female
hc_so_ADULTMvADULTF	<-	hc_metadata_L3fixed[(hc_metadata_L3fixed$name=="ADULT_M" | hc_metadata_L3fixed$name=="ADULT_F"),]
hc_so_ADULTMvADULTF <- sleuth_prep(hc_so_ADULTMvADULTF, extra_bootstrap_summary = TRUE,num_cores=1)
hc_so_ADULTMvADULTF <- sleuth_fit(hc_so_ADULTMvADULTF, ~name, 'full')
hc_so_ADULTMvADULTF <- sleuth_fit(hc_so_ADULTMvADULTF, ~1, 'reduced')
hc_so_ADULTMvADULTF <- sleuth_lrt(hc_so_ADULTMvADULTF, 'reduced', 'full')

#sleuth_table_ADULTMvADULTF <- sleuth_results(hc_so_ADULTMvADULTF, 'reduced:full', 'lrt', show_all = FALSE)
#sleuth_significant_ADULTMvADULTF <- dplyr::filter(sleuth_table_ADULTMvADULTF, qval <= 0.05)
#head(sleuth_significant_ADULTMvADULTF, 20)
#write.table(sleuth_table_ADULTMvADULTF,file="sleuth_table_ADULTMvADULTF.txt",sep="\t",quote=FALSE, row.names=FALSE)


hc_so_ADULTMvADULTF_wt<-sleuth_wt(hc_so_ADULTMvADULTF,'nameADULT_M',which_model = "full")
sleuth_table_ADULTMvADULTF_wt <- sleuth_results(hc_so_ADULTMvADULTF_wt,test="nameADULT_M",which_model = "full",test_type = 'wt')

# up
sleuth_significant_ADULTMvADULTF_wt_up <- dplyr::filter(sleuth_table_ADULTMvADULTF_wt, qval <= 0.01,b>=2)
write.table(sleuth_significant_ADULTMvADULTF_wt_up,file="sleuth_table_ADULTMvADULTF_up.txt",sep="\t",quote=FALSE, row.names=FALSE)
# down
sleuth_significant_ADULTMvADULTF_wt_down <- dplyr::filter(sleuth_table_ADULTMvADULTF_wt, qval <= 0.01,b<=-2)
write.table(sleuth_significant_ADULTMvADULTF_wt_down,file="sleuth_table_ADULTMvADULTF_down.txt",sep="\t",quote=FALSE, row.names=FALSE)


#####–-------------- TESTING REVERSE, with female first and in model rather than male

hc_so_ADULTMvADULTF	<-	hc_metadata_L3fixed[(hc_metadata_L3fixed$name=="ADULT_M" | hc_metadata_L3fixed$name=="ADULT_F"),]
hc_so_ADULTFvADULTM <- arrange(hc_so_ADULTMvADULTF, rev(rownames(hc_so_ADULTMvADULTF)))
hc_so_ADULTFvADULTM <- sleuth_prep(hc_so_ADULTFvADULTM, extra_bootstrap_summary = TRUE,num_cores=1)
hc_so_ADULTFvADULTM <- sleuth_fit(hc_so_ADULTFvADULTM, ~name, 'full')
hc_so_ADULTFvADULTM <- sleuth_fit(hc_so_ADULTFvADULTM, ~1, 'reduced')
hc_so_ADULTFvADULTM <- sleuth_lrt(hc_so_ADULTFvADULTM, 'reduced', 'full')


hc_so_ADULTFvADULTM_wt<-sleuth_wt(hc_so_ADULTFvADULTM,'nameADULT_F',which_model = "full")
sleuth_table_ADULTFvADULTM_wt <- sleuth_results(hc_so_ADULTFvADULTM_wt,test="nameADULT_F",which_model = "full",test_type = 'wt')

# reverse
hc_so_ADULTFvADULTM_wt<-sleuth_wt(hc_so_ADULTMvADULTF,'nameADULT_F',which_model = "full")
sleuth_table_ADULTFvADULTM_wt <- sleuth_results(hc_so_ADULTFvADULTM_wt,test="nameADULT_F",which_model = "full",test_type = 'wt')

# up
sleuth_significant_ADULTFvADULTM_wt_up <- dplyr::filter(sleuth_table_ADULTFvADULTM_wt, qval <= 0.01,b>=2)
write.table(sleuth_significant_ADULTFvADULTM_wt_up,file="sleuth_table_ADULTFvADULTM_up.txt",sep="\t",quote=FALSE, row.names=FALSE)
# down
sleuth_significant_ADULTFvADULTM_wt_down <- dplyr::filter(sleuth_table_ADULTFvADULTM_wt, qval <= 0.01,b<=-2)
write.table(sleuth_significant_ADULTFvADULTM_wt_down,file="sleuth_table_ADULTFvADULTM_down.txt",sep="\t",quote=FALSE, row.names=FALSE)

#####–-------------- TESTING REVERSE






# Adult Female vs Gut
hc_so_ADULTFvGUT	<-	hc_metadata_L3fixed[(hc_metadata_L3fixed$name=="ADULT_F" | hc_metadata_L3fixed$name=="GUT"),]
hc_so_ADULTFvGUT <- sleuth_prep(hc_so_ADULTFvGUT, extra_bootstrap_summary = TRUE,num_cores=1)
hc_so_ADULTFvGUT <- sleuth_fit(hc_so_ADULTFvGUT, ~name, 'full')
hc_so_ADULTFvGUT <- sleuth_fit(hc_so_ADULTFvGUT, ~1, 'reduced')
hc_so_ADULTFvGUT <- sleuth_lrt(hc_so_ADULTFvGUT, 'reduced', 'full')

#sleuth_table_ADULTFvGUT <- sleuth_results(hc_so_ADULTFvGUT, 'reduced:full', 'lrt', show_all = FALSE)
#sleuth_significant_ADULTFvGUT <- dplyr::filter(sleuth_table_ADULTFvGUT, qval <= 0.05)
#head(sleuth_significant_ADULTFvGUT, 20)
#write.table(sleuth_table_ADULTFvGUT,file="sleuth_table_ADULTFvGUT.txt",sep="\t",quote=FALSE, row.names=FALSE)


hc_so_GUTvADULTF_wt<-sleuth_wt(hc_so_ADULTFvGUT,'nameGUT',which_model = "full")
sleuth_table_GUTvADULTF_wt <- sleuth_results(hc_so_ADULTFvGUT_wt,test="nameGUT",which_model = "full",test_type = 'wt')
#sleuth_significant_ADULTFvGUT_wt <- dplyr::filter(sleuth_table_ADULTFvGUT_wt, qval <= 0.05)
#head(sleuth_significant_ADULTFvGUT_wt, 100)

# up
sleuth_significant_GUTvADULTF_wt_up <- dplyr::filter(sleuth_table_GUTvADULTF_wt, qval <= 0.01,b>=2)
write.table(sleuth_significant_GUTvADULTF_wt_up,file="sleuth_table_GUTvADULTF_up.txt",sep="\t",quote=FALSE, row.names=FALSE)
# down
sleuth_significant_GUTvADULTF_wt_down <- dplyr::filter(sleuth_table_GUTvADULTF_wt, qval <= 0.01,b<=-2)
write.table(sleuth_significant_GUTvADULTF_wt_down,file="sleuth_table_GUTvADULTF_down.txt",sep="\t",quote=FALSE, row.names=FALSE)





# Adult Female vs Egg
hc_so_ADULTFvEGG	<-	hc_metadata_L3fixed[(hc_metadata_L3fixed$name=="EGG" | hc_metadata_L3fixed$name=="ADULT_F"),]
hc_so_ADULTFvEGG <- sleuth_prep(hc_so_ADULTFvEGG, extra_bootstrap_summary = TRUE,num_cores=1)
hc_so_ADULTFvEGG <- sleuth_fit(hc_so_ADULTFvEGG, ~name, 'full')
hc_so_ADULTFvEGG <- sleuth_fit(hc_so_ADULTFvEGG, ~1, 'reduced')
hc_so_ADULTFvEGG <- sleuth_lrt(hc_so_ADULTFvEGG, 'reduced', 'full')

#sleuth_table_ADULTFvEGG <- sleuth_results(hc_so_ADULTFvEGG, 'reduced:full', 'lrt', show_all = FALSE)
#sleuth_significant_ADULTFvEGG <- dplyr::filter(sleuth_table_ADULTFvEGG, qval <= 0.05)
#head(sleuth_significant_ADULTFvEGG, 20)

#write.table(sleuth_table_ADULTFvEGG,file="sleuth_table_ADULTFvEGG.txt",sep="\t",quote=FALSE, row.names=FALSE)

hc_so_EGGvADULTF_wt	<-	sleuth_wt(hc_so_ADULTFvEGG,'nameEGG',which_model = "full")
sleuth_table_EGGvADULTF_wt	<-	sleuth_results(hc_so_EGGvADULTF_wt,test="nameEGG",which_model = "full",test_type = 'wt')
#sleuth_significant_EGGvADULTF_wt <- dplyr::filter(sleuth_table_EGGvADULTF_wt, qval <= 0.05)
#head(sleuth_significant_EGGvADULTF_wt, 100)
# up
sleuth_significant_EGGvADULTF_wt_up <- dplyr::filter(sleuth_table_EGGvADULTF_wt, qval <= 0.01,b>=2)
write.table(sleuth_significant_EGGvADULTF_wt_up,file="sleuth_table_EGGvADULTF_up.txt",sep="\t",quote=FALSE, row.names=FALSE)
# down
sleuth_significant_EGGvADULTF_wt_down <- dplyr::filter(sleuth_table_EGGvADULTF_wt, qval <= 0.01,b<=-2)
write.table(sleuth_significant_EGGvADULTF_wt_down,file="sleuth_table_EGGvADULTF_down.txt",sep="\t",quote=FALSE, row.names=FALSE)



save.image(file = "hc_genome_kallisto.RData")
```



# make a heatmap of top 1000 variable genes across all life stages

```shell
cd /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/TRANSCRIPTOME/KALLISTO/KALLISTO_MAPPED_SAMPLES
```

### generate the tpm data
```shell
# extract TPMs per sample
for i in ` ls -1d *out `; do echo $i > ${i}.tpm ; cat ${i}/abundance.tsv | cut -f5 | sed '1d' >> ${i}.tpm; done

# generate a "transcripts" list, taken from the TRANSCRIPTS.fa file
#echo "ID" > transcripts.list; grep ">" ../TRANSCRIPTS.fa | cut -f1 -d  " " | sed 's/>//g' >> transcripts.list
# due to Apollo giving long unique codes, the transcript IDs are obscure. Here is the fix
#awk '$3=="mRNA" {print $9}' ../ANNOTATION.gff3 | cut -f3,5 -d";" | sed -e 's/ID=//g' -e 's/;Name=/\t/g' > mRNA_IDtoNAME_conversion.txt

#while read ID NAME; do sed -i "s/${ID}/${NAME}/g" transcripts.list; done < mRNA_IDtoNAME_conversion.txt &

# ALTERNATE WAY, direct from the annotaiton
echo "ID" > transcripts.list; grep ">" ../TRANSCRIPTS.fa | cut -f1 -d" " | sed -e 's/>//g' >> transcripts.list



# make a data frame containing all TMP values from all samples
paste transcripts.list \
kallisto_7059_6_1_out.tpm \
kallisto_7059_6_2_out.tpm \
kallisto_7059_6_3_out.tpm \
kallisto_7059_6_4_out.tpm \
kallisto_7059_6_5_out.tpm \
kallisto_7059_6_6_out.tpm \
kallisto_7062_6_8_out.tpm \
kallisto_7062_6_10_out.tpm \
kallisto_7062_6_11_out.tpm \
kallisto_7062_6_9_out.tpm \
kallisto_7062_6_12_out.tpm \
kallisto_7059_6_7_out.tpm \
kallisto_7059_6_8_out.tpm \
kallisto_7059_6_9_out.tpm \
kallisto_7062_6_1_out.tpm \
kallisto_7062_6_2_out.tpm \
kallisto_7062_6_3_out.tpm \
kallisto_7059_6_10_out.tpm \
kallisto_7059_6_11_out.tpm \
kallisto_7059_6_12_out.tpm \
kallisto_7062_6_13_out.tpm \
kallisto_7062_6_14_out.tpm \
kallisto_7062_6_15_out.tpm \
> kallisto_allsamples.tpm.table

#kallisto_7062_6_7_out.tpm  has been removed
```

Curate the data, including:
  - setting minimum TMP at 1
  - transforming to log10 scale
  - calculating variance per row, of which the top 1000 most variable rows are selected
  - plot heatmap, of most variable transcripts



```R
R-3.5.0
library(gplots)
library(tibble)
library(RColorBrewer)

data<-read.table("kallisto_allsamples.tpm.table",header=T,row.names=1)

# set a TPM cutoff,
data<-(data > 1) * (data - 1) + 1
data<-log10(data)
data<-as.matrix(data)

# fix infinite values to NA
is.na(data) <- sapply(data, is.infinite)

# calculate variance per row
#https://stackoverflow.com/questions/25099825/row-wise-variance-of-a-matrix-in-r
RowVar <- function(x, ...) {
  rowSums(na.rm=TRUE,(x - rowMeans(x, ...))^2, ...)/(dim(x)[2] - 1)
}
var<-as.matrix(RowVar(data))
var<-as.data.frame(var)
var <- var[order(var), ,drop = FALSE]
var_filter <- tail(var,1000)   # set number of genes here
var_filter <- rownames_to_column(var_filter)

# bring all data together
data<-as.data.frame(data)
data<-rownames_to_column(data)

data_filtered <- dplyr::semi_join(data, var_filter, by = "rowname")
data_filtered <- column_to_rownames(data_filtered,'rowname')
data_filtered<-as.matrix(data_filtered)
#var<-as.matrix(RowVar(data))
#data<-cbind(data, variance = var )

# make heatmap
pdf("top1000variablegenes_allstages_minTPM1.pdf")
heatmap.2(data_filtered,trace="none",na.color="grey",labRow=F,dendrogram='row',Colv=FALSE,col= colorRampPalette(brewer.pal(8, "Blues"))(25))
dev.off()

png("top1000variablegenes_allstages_minTPM1.png")
heatmap.2(data_filtered,trace="none",na.color="grey",labRow=F,dendrogram='row',Colv=FALSE,col= colorRampPalette(brewer.pal(8, "Blues"))(25))
dev.off()

```

- Copy to local dir - run this from local machine
```shell
scp sd21@pcs5.internal.sanger.ac.uk:/nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/TRANSCRIPTOME/KALLISTO/KALLISTO_MAPPED_SAMPLES/top1000variablegenes_allstages_minTPM1.* ~/Documents/workbook/hcontortus_genome/04_analysis
```

![Kallisto - top 1000 most variable genes across lifestages](04_analysis/top1000variablegenes_allstages_minTPM1.png)
Fig - Kalliso  top 1000 most variable genes across lifestages




### Run clustering analysis of gene expression across the life stages

Need to generate a replicates file that tells CLUST what samples to group.

### collate data for clust
```shell
mkdir clust_data
cp kallisto_allsamples.tpm.table clust_data/
```

Run clust
- requires a replicates dataset - see "replicates.txt"

kallisto_allsamples.tpm.table   EGG     kallisto_7059_6_1_out   kallisto_7059_6_2_out   kallisto_7059_6_3_out
kallisto_allsamples.tpm.table   L1      kallisto_7059_6_4_out   kallisto_7059_6_5_out   kallisto_7059_6_6_out
kallisto_allsamples.tpm.table   SHL3    kallisto_7062_6_10_out  kallisto_7062_6_11_out  kallisto_7062_6_8_out
kallisto_allsamples.tpm.table   EXL3    kallisto_7062_6_12_out  kallisto_7062_6_9_out
kallisto_allsamples.tpm.table   L4      kallisto_7059_6_7_out   kallisto_7059_6_8_out   kallisto_7059_6_9_out
kallisto_allsamples.tpm.table   ADULT_F kallisto_7059_6_10_out  kallisto_7059_6_11_out  kallisto_7059_6_12_out
kallisto_allsamples.tpm.table   ADULT_M kallisto_7062_6_1_out   kallisto_7062_6_2_out   kallisto_7062_6_3_out
kallisto_allsamples.tpm.table   GUT     kallisto_7062_6_13_out  kallisto_7062_6_14_out  kallisto_7062_6_15_out


```shell
# run clust
clust $PWD/clust_data -r replicates.txt
```
Clust produces a cluster profile PDF containing all Clusters. However, it is not that nice, and so will make better ones using ggplot. Need to convert PDF to PNG to post however.

```shell
# pdf to png conversion
cd Results
convert Clusters_profiles.pdf -quality 200  Clusters_profiles.png
```

- Copy to local dir - run this from local machine
```shell
scp sd21@pcs5.internal.sanger.ac.uk:/nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/TRANSCRIPTOME/KALLISTO/KALLISTO_MAPPED_SAMPLES/Results_*/Clusters_profiles*.png ~/Documents/workbook/hcontortus_genome/04_analysis/
```

![Clust - cluster profiles ](04_analysis/Clusters_profiles-0.png)
![Clust - cluster profiles ](04_analysis/Clusters_profiles-1.png)
Fig - Clust - correlated gene expression across lifestages

Results
- The largest cluster expression profiles split Egg/L1/L3 and L4/Adults
-






Make nice clust plots using ggplot2
```shell
cd /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/TRANSCRIPTOME/KALLISTO/KALLISTO_MAPPED_SAMPLES/Results_26_Jan_19

# make gene lists for each cluster set
for i in {1..19}; do cut -f "${i}" Clusters_Objects.tsv | sed '2d' | cut -f1 -d " " | sed '/^$/d' > cluster_${i}.list; done
```
```R
R-3.5.0
library(dplyr)
library(reshape2)
library(ggplot2)
library(gtools)

myFiles <- list.files(pattern="*.list")
myFiles <- mixedsort(sort(myFiles))   # get files in order
data<-read.table("Processed_Data/kallisto_allsamples.tpm.table_processed.tsv",header=T)


c0<-read.table("c1.list",header=T)

for(i in 1:length(myFiles)){
		cluster_data <- read.table(myFiles[i],header=TRUE)
		data<-read.table("Processed_Data/kallisto_allsamples.tpm.table_processed.tsv",header=T)
		data_cluster<-dplyr::semi_join(data, cluster_data, by = c("Genes" = colnames(cluster_data[1])))
		df_melted <- melt(data_cluster, id.vars = 'Genes')
  		ggplot()+geom_line(aes(df_melted$variable,df_melted$value,group = df_melted$Genes),alpha=0.1)+theme_bw()+ylab("Normalised TMP")+ylim(-2.5,2.5)
		ggsave(paste("cluster",i,".pdf",sep=""),width = 6, height = 6, units = c("cm"))
}

# example plot - c0
data_c0<-dplyr::semi_join(data, c0, by = c("Genes" = "C0"))
df_melted <- melt(data_c0, id.vars = 'Genes')
ggplot()+geom_line(aes(df_melted$variable,df_melted$value,group = df_melted$Genes),alpha=0.1)+theme_bw()

```

Import into Illustrator to make a nice figure.

- Copy to local dir - run this from local machine
```shell
scp sd21@pcs5.internal.sanger.ac.uk:/nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/TRANSCRIPTOME/KALLISTO/KALLISTO_MAPPED_SAMPLES/Results_*/cluster*.pdf ~/Documents/workbook/hcontortus_genome/04_analysis
```

### Distribution of clustered gene expression profiles across the genome

Want to explore whether genome locate has any impact on the coexpression of genes.

# get genome coordinates of genes in clusters
```shell
ln -s /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/TRANSCRIPTOME/TRANSCRIPTOME_CURATION/HCON_V4_WBP11plus_190125.ips.gff3 ANNOTATION.gff

while read NAME; do grep "ID=${NAME};" ANNOTATION.gff | awk '{if($3=="mRNA") print $1,$4,$5,$7,$9}' OFS="\t" ; done < cluster_1.list > cluster_1.coords
while read NAME; do grep "ID=${NAME};" ANNOTATION.gff | awk '{if($3=="mRNA") print $1,$4,$5,$7,$9}' OFS="\t" ; done < cluster_2.list > cluster_2.coords
while read NAME; do grep "ID=${NAME};" ANNOTATION.gff | awk '{if($3=="mRNA") print $1,$4,$5,$7,$9}' OFS="\t" ; done < cluster_3.list > cluster_3.coords
while read NAME; do grep "ID=${NAME};" ANNOTATION.gff | awk '{if($3=="mRNA") print $1,$4,$5,$7,$9}' OFS="\t" ; done < cluster_4.list > cluster_4.coords
while read NAME; do grep "ID=${NAME};" ANNOTATION.gff | awk '{if($3=="mRNA") print $1,$4,$5,$7,$9}' OFS="\t" ; done < cluster_5.list > cluster_5.coords
while read NAME; do grep "ID=${NAME};" ANNOTATION.gff | awk '{if($3=="mRNA") print $1,$4,$5,$7,$9}' OFS="\t" ; done < cluster_6.list > cluster_6.coords
while read NAME; do grep "ID=${NAME};" ANNOTATION.gff | awk '{if($3=="mRNA") print $1,$4,$5,$7,$9}' OFS="\t" ; done < cluster_7.list > cluster_7.coords
while read NAME; do grep "ID=${NAME};" ANNOTATION.gff | awk '{if($3=="mRNA") print $1,$4,$5,$7,$9}' OFS="\t" ; done < cluster_8.list > cluster_8.coords
while read NAME; do grep "ID=${NAME};" ANNOTATION.gff | awk '{if($3=="mRNA") print $1,$4,$5,$7,$9}' OFS="\t" ; done < cluster_9.list > cluster_9.coords
while read NAME; do grep "ID=${NAME};" ANNOTATION.gff | awk '{if($3=="mRNA") print $1,$4,$5,$7,$9}' OFS="\t" ; done < cluster_10.list > cluster_10.coords
while read NAME; do grep "ID=${NAME};" ANNOTATION.gff | awk '{if($3=="mRNA") print $1,$4,$5,$7,$9}' OFS="\t" ; done < cluster_11.list > cluster_11.coords
while read NAME; do grep "ID=${NAME};" ANNOTATION.gff | awk '{if($3=="mRNA") print $1,$4,$5,$7,$9}' OFS="\t" ; done < cluster_12.list > cluster_12.coords
while read NAME; do grep "ID=${NAME};" ANNOTATION.gff | awk '{if($3=="mRNA") print $1,$4,$5,$7,$9}' OFS="\t" ; done < cluster_13.list > cluster_13.coords
while read NAME; do grep "ID=${NAME};" ANNOTATION.gff | awk '{if($3=="mRNA") print $1,$4,$5,$7,$9}' OFS="\t" ; done < cluster_14.list > cluster_14.coords
while read NAME; do grep "ID=${NAME};" ANNOTATION.gff | awk '{if($3=="mRNA") print $1,$4,$5,$7,$9}' OFS="\t" ; done < cluster_15.list > cluster_15.coords
while read NAME; do grep "ID=${NAME};" ANNOTATION.gff | awk '{if($3=="mRNA") print $1,$4,$5,$7,$9}' OFS="\t" ; done < cluster_16.list > cluster_16.coords
while read NAME; do grep "ID=${NAME};" ANNOTATION.gff | awk '{if($3=="mRNA") print $1,$4,$5,$7,$9}' OFS="\t" ; done < cluster_17.list > cluster_17.coords
while read NAME; do grep "ID=${NAME};" ANNOTATION.gff | awk '{if($3=="mRNA") print $1,$4,$5,$7,$9}' OFS="\t" ; done < cluster_18.list > cluster_18.coords
while read NAME; do grep "ID=${NAME};" ANNOTATION.gff | awk '{if($3=="mRNA") print $1,$4,$5,$7,$9}' OFS="\t" ; done < cluster_19.list > cluster_19.coords
while read NAME; do grep "ID=${NAME};" ANNOTATION.gff | awk '{if($3=="mRNA") print $1,$4,$5,$7,$9}' OFS="\t" ; done < cluster_20.list > cluster_20.coords
```


# stat testing the clusters
```shell
ln -s ../../../../REF/HAEM_V4_final.chr.fa
samtools faidx HAEM_V4_final.chr.fa
cut -f1,2 HAEM_V4_final.chr.fa.fai > HAEM_V4_final.chr.genome

# make bed files - note per cluster bed files are already made, ie. cluster_20.coords
bedtools-2 makewindows -g  HAEM_V4_final.chr.genome -w 500000 > HAEM_V4_final.chr.500k.bed
awk '$3=="mRNA" {print $1,$4,$5,$9}' OFS="\t" ANNOTATION.gff >HCON_V4.mRNA.bed


bedtools-2 coverage -a HAEM_V4_final.chr.500k.bed -b HCON_V4.mRNA.bed -counts > HCON_V4.mRNA.counts

for i in *.coords; do bedtools-2 coverage -a HAEM_V4_final.chr.500k.bed -b <( sort -k1,1 -k2,2n ${i}) -counts > ${i%.coords}.counts; done
```




Make the plots
```R
R-3.5.0
library(ggplot2)
library(patchwork)
library(viridis)

cluster_1 <- read.table("cluster_1.coords",header=F)
cluster_2 <- read.table("cluster_2.coords",header=F)
cluster_3 <- read.table("cluster_3.coords",header=F)
cluster_4 <- read.table("cluster_4.coords",header=F)
cluster_5 <- read.table("cluster_5.coords",header=F)
cluster_6 <- read.table("cluster_6.coords",header=F)
cluster_7 <- read.table("cluster_7.coords",header=F)
cluster_8 <- read.table("cluster_8.coords",header=F)
cluster_9 <- read.table("cluster_9.coords",header=F)
cluster_10 <- read.table("cluster_10.coords",header=F)
cluster_11 <- read.table("cluster_11.coords",header=F)
cluster_12 <- read.table("cluster_12.coords",header=F)
cluster_13 <- read.table("cluster_13.coords",header=F)
cluster_14 <- read.table("cluster_14.coords",header=F)
cluster_15 <- read.table("cluster_15.coords",header=F)
cluster_16 <- read.table("cluster_16.coords",header=F)
cluster_17 <- read.table("cluster_17.coords",header=F)
cluster_18 <- read.table("cluster_18.coords",header=F)
cluster_19 <- read.table("cluster_19.coords",header=F)
cluster_20 <- read.table("cluster_20.coords",header=F)

cluster_1$cluster <- 1
cluster_2$cluster <- 2
cluster_3$cluster <- 3
cluster_4$cluster <- 4
cluster_5$cluster <- 5
cluster_6$cluster <- 6
cluster_7$cluster <- 7
cluster_8$cluster <- 8
cluster_9$cluster <- 9
cluster_10$cluster <- 10
cluster_11$cluster <- 11
cluster_12$cluster <- 12
cluster_13$cluster <- 13
cluster_14$cluster <- 14
cluster_15$cluster <- 15
cluster_16$cluster <- 16
cluster_17$cluster <- 17
cluster_18$cluster <- 18
cluster_19$cluster <- 19
cluster_20$cluster <- 20


clusters <- rbind(cluster_1,cluster_2,cluster_3,cluster_4,cluster_5,cluster_6,cluster_7,cluster_8,cluster_9,cluster_10,cluster_11,cluster_12,cluster_13,cluster_14,cluster_15,cluster_16,cluster_17,cluster_18,cluster_19)

#ggplot(clusters,aes(clusters$V2,clusters$cluster))+ geom_jitter(alpha=0.1,size=0.5)+facet_grid(.~clusters$V1)+theme_bw()+ scale_y_continuous(breaks=seq(1,20,1),trans = 'reverse')


all<-read.table("HCON_V4.mRNA.counts",header=F)
c1<-read.table("cluster_1.counts",header=F)
c2<-read.table("cluster_2.counts",header=F)
c3<-read.table("cluster_3.counts",header=F)
c4<-read.table("cluster_4.counts",header=F)
c5<-read.table("cluster_5.counts",header=F)
c6<-read.table("cluster_6.counts",header=F)
c7<-read.table("cluster_7.counts",header=F)
c8<-read.table("cluster_8.counts",header=F)
c9<-read.table("cluster_9.counts",header=F)
c10<-read.table("cluster_10.counts",header=F)
c11<-read.table("cluster_11.counts",header=F)
c12<-read.table("cluster_12.counts",header=F)
c13<-read.table("cluster_13.counts",header=F)
c14<-read.table("cluster_14.counts",header=F)
c15<-read.table("cluster_15.counts",header=F)
c16<-read.table("cluster_16.counts",header=F)
c17<-read.table("cluster_17.counts",header=F)
c18<-read.table("cluster_18.counts",header=F)
c19<-read.table("cluster_19.counts",header=F)
#c20<-read.table("cluster_20.counts",header=F)


# calculate the expected number of clustered transcripts per window with the null hypothesis that they are equally distributed throughout the genome.
c1_stat<-all
c1_stat$observed<-c1$V4
c1_stat$expected<-all$V4*(sum(c1$V4)/21007)
c1_stat$chqsq<-(c1$V4 - all$V4*(sum(c1$V4)/21007))^2/(all$V4*sum(c1$V4)/21007)
c1_stat$pvalue <- pchisq(c1_stat$chqsq,df=1,lower.tail=FALSE)
c1_stat$cluster <- 1

c2_stat<-all
c2_stat$observed<-c2$V4
c2_stat$expected<-all$V4*(sum(c2$V4)/21007)
c2_stat$chqsq<-(c2$V4 - all$V4*(sum(c2$V4)/21007))^2/(all$V4*sum(c2$V4)/21007)
c2_stat$pvalue <- pchisq(c2_stat$chqsq,df=1,lower.tail=FALSE)
c2_stat$cluster <- 2

c3_stat<-all
c3_stat$observed<-c3$V4
c3_stat$expected<-all$V4*(sum(c3$V4)/21007)
c3_stat$chqsq<-(c3$V4 - all$V4*(sum(c3$V4)/21007))^2/(all$V4*sum(c3$V4)/21007)
c3_stat$pvalue <- pchisq(c3_stat$chqsq,df=1,lower.tail=FALSE)
c3_stat$cluster <- 3

c4_stat<-all
c4_stat$observed<-c4$V4
c4_stat$expected<-all$V4*(sum(c4$V4)/21007)
c4_stat$chqsq<-(c4$V4 - all$V4*(sum(c4$V4)/21007))^2/(all$V4*sum(c4$V4)/21007)
c4_stat$pvalue <- pchisq(c4_stat$chqsq,df=1,lower.tail=FALSE)
c4_stat$cluster <- 4

c5_stat<-all
c5_stat$observed<-c5$V4
c5_stat$expected<-all$V4*(sum(c5$V4)/21007)
c5_stat$chqsq<-(c5$V4 - all$V4*(sum(c5$V4)/21007))^2/(all$V4*sum(c5$V4)/21007)
c5_stat$pvalue <- pchisq(c5_stat$chqsq,df=1,lower.tail=FALSE)
c5_stat$cluster <- 5

c6_stat<-all
c6_stat$observed<-c6$V4
c6_stat$expected<-all$V4*(sum(c6$V4)/21007)
c6_stat$chqsq<-(c6$V4 - all$V4*(sum(c6$V4)/21007))^2/(all$V4*sum(c6$V4)/21007)
c6_stat$pvalue <- pchisq(c6_stat$chqsq,df=1,lower.tail=FALSE)
c6_stat$cluster <- 6

c7_stat<-all
c7_stat$observed<-c7$V4
c7_stat$expected<-all$V4*(sum(c7$V4)/21007)
c7_stat$chqsq<-(c7$V4 - all$V4*(sum(c7$V4)/21007))^2/(all$V4*sum(c7$V4)/21007)
c7_stat$pvalue <- pchisq(c7_stat$chqsq,df=1,lower.tail=FALSE)
c7_stat$cluster <- 7

c8_stat<-all
c8_stat$observed<-c8$V4
c8_stat$expected<-all$V4*(sum(c8$V4)/21007)
c8_stat$chqsq<-(c8$V4 - all$V4*(sum(c8$V4)/21007))^2/(all$V4*sum(c8$V4)/21007)
c8_stat$pvalue <- pchisq(c8_stat$chqsq,df=1,lower.tail=FALSE)
c8_stat$cluster <- 8

c9_stat<-all
c9_stat$observed<-c9$V4
c9_stat$expected<-all$V4*(sum(c9$V4)/21007)
c9_stat$chqsq<-(c9$V4 - all$V4*(sum(c9$V4)/21007))^2/(all$V4*sum(c9$V4)/21007)
c9_stat$pvalue <- pchisq(c9_stat$chqsq,df=1,lower.tail=FALSE)
c9_stat$cluster <- 9

c10_stat<-all
c10_stat$observed<-c10$V4
c10_stat$expected<-all$V4*(sum(c10$V4)/21007)
c10_stat$chqsq<-(c10$V4 - all$V4*(sum(c10$V4)/21007))^2/(all$V4*sum(c10$V4)/21007)
c10_stat$pvalue <- pchisq(c10_stat$chqsq,df=1,lower.tail=FALSE)
c10_stat$cluster <- 10

c11_stat<-all
c11_stat$observed<-c11$V4
c11_stat$expected<-all$V4*(sum(c11$V4)/21007)
c11_stat$chqsq<-(c11$V4 - all$V4*(sum(c11$V4)/21007))^2/(all$V4*sum(c11$V4)/21007)
c11_stat$pvalue <- pchisq(c11_stat$chqsq,df=1,lower.tail=FALSE)
c11_stat$cluster <- 11

c12_stat<-all
c12_stat$observed<-c12$V4
c12_stat$expected<-all$V4*(sum(c12$V4)/21007)
c12_stat$chqsq<-(c12$V4 - all$V4*(sum(c12$V4)/21007))^2/(all$V4*sum(c12$V4)/21007)
c12_stat$pvalue <- pchisq(c12_stat$chqsq,df=1,lower.tail=FALSE)
c12_stat$cluster <- 12

c13_stat<-all
c13_stat$observed<-c13$V4
c13_stat$expected<-all$V4*(sum(c13$V4)/21007)
c13_stat$chqsq<-(c13$V4 - all$V4*(sum(c13$V4)/21007))^2/(all$V4*sum(c13$V4)/21007)
c13_stat$pvalue <- pchisq(c13_stat$chqsq,df=1,lower.tail=FALSE)
c13_stat$cluster <- 13

c14_stat<-all
c14_stat$observed<-c14$V4
c14_stat$expected<-all$V4*(sum(c14$V4)/21007)
c14_stat$chqsq<-(c14$V4 - all$V4*(sum(c14$V4)/21007))^2/(all$V4*sum(c14$V4)/21007)
c14_stat$pvalue <- pchisq(c14_stat$chqsq,df=1,lower.tail=FALSE)
c14_stat$cluster <- 14

c15_stat<-all
c15_stat$observed<-c15$V4
c15_stat$expected<-all$V4*(sum(c15$V4)/21007)
c15_stat$chqsq<-(c15$V4 - all$V4*(sum(c15$V4)/21007))^2/(all$V4*sum(c15$V4)/21007)
c15_stat$pvalue <- pchisq(c15_stat$chqsq,df=1,lower.tail=FALSE)
c15_stat$cluster <- 15

c16_stat<-all
c16_stat$observed<-c16$V4
c16_stat$expected<-all$V4*(sum(c16$V4)/21007)
c16_stat$chqsq<-(c16$V4 - all$V4*(sum(c16$V4)/21007))^2/(all$V4*sum(c16$V4)/21007)
c16_stat$pvalue <- pchisq(c16_stat$chqsq,df=1,lower.tail=FALSE)
c16_stat$cluster <- 16


c17_stat<-all
c17_stat$observed<-c17$V4
c17_stat$expected<-all$V4*(sum(c17$V4)/21007)
c17_stat$chqsq<-(c17$V4 - all$V4*(sum(c17$V4)/21007))^2/(all$V4*sum(c17$V4)/21007)
c17_stat$pvalue <- pchisq(c17_stat$chqsq,df=1,lower.tail=FALSE)
c17_stat$cluster <- 17

c18_stat<-all
c18_stat$observed<-c18$V4
c18_stat$expected<-all$V4*(sum(c18$V4)/21007)
c18_stat$chqsq<-(c18$V4 - all$V4*(sum(c18$V4)/21007))^2/(all$V4*sum(c18$V4)/21007)
c18_stat$pvalue <- pchisq(c18_stat$chqsq,df=1,lower.tail=FALSE)
c18_stat$cluster <- 18

c19_stat<-all
c19_stat$observed<-c19$V4
c19_stat$expected<-all$V4*(sum(c19$V4)/21007)
c19_stat$chqsq<-(c19$V4 - all$V4*(sum(c19$V4)/21007))^2/(all$V4*sum(c19$V4)/21007)
c19_stat$pvalue <- pchisq(c19_stat$chqsq,df=1,lower.tail=FALSE)
c19_stat$cluster <- 19

#c20_stat<-all
#c20_stat$observed<-c20$V4
#c20_stat$expected<-all$V4*(sum(c20$V4)/21007)
#c20_stat$chqsq<-(c20$V4 - all$V4*(sum(c20$V4)/21007))^2/(all$V4*sum(c20$V4)/21007)
#c20_stat$pvalue <- pchisq(c20_stat$chqsq,df=1,lower.tail=FALSE)
#c20_stat$cluster <- 20

cluster_stats <- rbind(c1_stat,c2_stat,c3_stat,c4_stat,c5_stat,c6_stat,c7_stat,c8_stat,c9_stat,c10_stat,c11_stat,c12_stat,c13_stat,c14_stat,c15_stat,c16_stat,c17_stat,c18_stat,c19_stat)
cluster_stats <- cluster_stats[cluster_stats$V1!="hcontortus_chr_mtDNA_arrow_pilon",]
cluster_stats[is.na(cluster_stats)] <- 1
cluster_stats$neglogP <- -log10(cluster_stats$pvalue)


# make the final plots
plot_clusts<-ggplot(clusters,aes(clusters$V2,clusters$cluster))+
		geom_jitter(alpha=0.1,size=0.5)+facet_grid(.~clusters$V1)+
		theme_bw()+
		scale_y_continuous(breaks=seq(1,20,1),trans = 'reverse')

plot_stats<-ggplot(cluster_stats)+
		geom_rect(aes(xmin=cluster_stats$V2,ymin=cluster_stats$cluster-0.5,xmax=cluster_stats$V3,ymax=cluster_stats$cluster+0.5,fill=-log10(cluster_stats$pvalue)))+
		facet_grid(.~cluster_stats$V1)+
		theme_bw()+
		scale_fill_viridis(direction=-1)+
		scale_y_continuous(breaks=seq(1,20,1),trans = 'reverse')

plot_clusts + plot_stats + plot_layout(ncol=2)

ggsave("clust_profiles_across_genome.pdf",width = 28, height = 10, units = "cm")
ggsave("clust_profiles_across_genome.png",width = 28, height = 10, units = "cm")
```
- Copy to local dir - run this from local machine
```shell
scp sd21@pcs5.internal.sanger.ac.uk:/nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/TRANSCRIPTOME/KALLISTO/KALLISTO_MAPPED_SAMPLES/Results_*/clust_profiles_across_genome.* ~/Documents/workbook/hcontortus_genome/04_analysis
```

![CLust - expression profiles across genome](04_analysis/clust_profiles_across_genome.png)
Fig - CLust - expression profiles across genome



##########################################################################################
# normalised gene expression per chromosome per lifestage
##########################################################################################


cd /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/TRANSCRIPTOME/KALLISTO/KALLISTO_MAPPED_SAMPLES/Results_28_Oct_18_1/Processed_Data

#--- get mRNA sequences from gff
ln -s ../../../../HCON_V4.renamed.gff3
awk '$3=="mRNA" {print $1,$4,$5,$7,$9}' HCON_V4.renamed.gff3 > mRNA.coords

#--- extract mRNA ids from normalised expresison dataset, and  pull coordinates
cut -f1 all.tpm_processed.tsv | while read -r NAME; do grep "ID=${NAME};" mRNA.coords; done > coords
echo -e "CHR\tSTART\tEND\tSTRAND\tID" > tmp; cat tmp coords > coords.tmp

#--- merge coordinates and normalised expression data
paste coords.tmp all.tpm_processed.tsv | sort -k1,1 -k2,2n > coords.tmp.txt

#--- extract data per chromosome
grep "hcontortus_chr1_Celeg_TT_arrow_pilon" coords.tmp.txt > chr1.coords.tmp.txt
grep "hcontortus_chr2_Celeg_TT_arrow_pilon" coords.tmp.txt > chr2.coords.tmp.txt
grep "hcontortus_chr3_Celeg_TT_arrow_pilon" coords.tmp.txt > chr3.coords.tmp.txt
grep "hcontortus_chr4_Celeg_TT_arrow_pilon" coords.tmp.txt > chr4.coords.tmp.txt
grep "hcontortus_chr5_Celeg_TT_arrow_pilon" coords.tmp.txt > chr5.coords.tmp.txt
grep "hcontortus_chrX_Celeg_TT_arrow_pilon" coords.tmp.txt > chrX.coords.tmp.txt





R-3.5.0
library(ggplot2)
library(patchwork)

c1 <- read.table("chr1.coords.tmp.txt",header=F)
c2 <- read.table("chr2.coords.tmp.txt",header=F)
c3 <- read.table("chr3.coords.tmp.txt",header=F)
c4 <- read.table("chr4.coords.tmp.txt",header=F)
c5 <- read.table("chr5.coords.tmp.txt",header=F)
cX <- read.table("chrX.coords.tmp.txt",header=F)


egg_plot <- ggplot()+
		geom_jitter(aes(c1$V2,1,col=c1$V7),size=0.5)+
		geom_jitter(aes(c2$V2,2,col=c2$V7),size=0.5)+
		geom_jitter(aes(c3$V2,3,col=c3$V7),size=0.5)+
		geom_jitter(aes(c4$V2,4,col=c4$V7),size=0.5)+
		geom_jitter(aes(c5$V2,5,col=c5$V7),size=0.5)+
		geom_jitter(aes(cX$V2,6,col=cX$V7),size=0.5)+
		scale_colour_gradient2()+
		scale_y_continuous(breaks=seq(1,6,1),trans = 'reverse')+
		ggtitle("Egg")+xlab("Genome position (bp)")+ylab("Normalised expression")+
		theme_bw()
l1_plot <- ggplot()+
		geom_jitter(aes(c1$V2,1,col=c1$V8),size=0.5)+
		geom_jitter(aes(c2$V2,2,col=c2$V8),size=0.5)+
		geom_jitter(aes(c3$V2,3,col=c3$V8),size=0.5)+
		geom_jitter(aes(c4$V2,4,col=c4$V8),size=0.5)+
		geom_jitter(aes(c5$V2,5,col=c5$V8),size=0.5)+
		geom_jitter(aes(cX$V2,6,col=cX$V8),size=0.5)+
		scale_colour_gradient2()+
		scale_y_continuous(breaks=seq(1,6,1),trans = 'reverse')+
		ggtitle("L1")+xlab("Genome position (bp)")+ylab("Normalised expression")+
		theme_bw()
shl3_plot <- ggplot()+
		geom_jitter(aes(c1$V2,1,col=c1$V9),size=0.5)+
		geom_jitter(aes(c2$V2,2,col=c2$V9),size=0.5)+
		geom_jitter(aes(c3$V2,3,col=c3$V9),size=0.5)+
		geom_jitter(aes(c4$V2,4,col=c4$V9),size=0.5)+
		geom_jitter(aes(c5$V2,5,col=c5$V9),size=0.5)+
		geom_jitter(aes(cX$V2,6,col=cX$V9),size=0.5)+
		scale_colour_gradient2()+
		scale_y_continuous(breaks=seq(1,6,1),trans = 'reverse')+
		ggtitle("SHL3")+xlab("Genome position (bp)")+ylab("Normalised expression")+
		theme_bw()
exl3_plot <- ggplot()+
		geom_jitter(aes(c1$V2,1,col=c1$V10),size=0.5)+
		geom_jitter(aes(c2$V2,2,col=c2$V10),size=0.5)+
		geom_jitter(aes(c3$V2,3,col=c3$V10),size=0.5)+
		geom_jitter(aes(c4$V2,4,col=c4$V10),size=0.5)+
		geom_jitter(aes(c5$V2,5,col=c5$V10),size=0.5)+
		geom_jitter(aes(cX$V2,6,col=cX$V10),size=0.5)+
		scale_colour_gradient2()+
		scale_y_continuous(breaks=seq(1,6,1),trans = 'reverse')+
		ggtitle("EXL3")+xlab("Genome position (bp)")+ylab("Normalised expression")+
		theme_bw()
l4_plot <- ggplot()+
		geom_jitter(aes(c1$V2,1,col=c1$V11),size=0.5)+
		geom_jitter(aes(c2$V2,2,col=c2$V11),size=0.5)+
		geom_jitter(aes(c3$V2,3,col=c3$V11),size=0.5)+
		geom_jitter(aes(c4$V2,4,col=c4$V11),size=0.5)+
		geom_jitter(aes(c5$V2,5,col=c5$V11),size=0.5)+
		geom_jitter(aes(cX$V2,6,col=cX$V11),size=0.5)+
		scale_colour_gradient2()+
		scale_y_continuous(breaks=seq(1,6,1),trans = 'reverse')+
		ggtitle("L4")+xlab("Genome position (bp)")+ylab("Normalised expression")+
		theme_bw()
adultm_plot <- ggplot()+
		geom_jitter(aes(c1$V2,1,col=c1$V12),size=0.5)+
		geom_jitter(aes(c2$V2,2,col=c2$V12),size=0.5)+
		geom_jitter(aes(c3$V2,3,col=c3$V12),size=0.5)+
		geom_jitter(aes(c4$V2,4,col=c4$V12),size=0.5)+
		geom_jitter(aes(c5$V2,5,col=c5$V12),size=0.5)+
		geom_jitter(aes(cX$V2,6,col=cX$V12),size=0.5)+
		scale_colour_gradient2()+
		scale_y_continuous(breaks=seq(1,6,1),trans = 'reverse')+
		ggtitle("Adult Male")+xlab("Genome position (bp)")+ylab("Normalised expression")+
		theme_bw()
adultf_plot <- ggplot()+
		geom_jitter(aes(c1$V2,1,col=c1$V13),size=0.5)+
		geom_jitter(aes(c2$V2,2,col=c2$V13),size=0.5)+
		geom_jitter(aes(c3$V2,3,col=c3$V13),size=0.5)+
		geom_jitter(aes(c4$V2,4,col=c4$V13),size=0.5)+
		geom_jitter(aes(c5$V2,5,col=c5$V13),size=0.5)+
		geom_jitter(aes(cX$V2,6,col=cX$V13),size=0.5)+
		scale_colour_gradient2()+
		scale_y_continuous(breaks=seq(1,6,1),trans = 'reverse')+
		ggtitle("Adult Female")+xlab("Genome position (bp)")+ylab("Normalised expression")+
		theme_bw()
gut_plot <- ggplot()+
		geom_jitter(aes(c1$V2,1,col=c1$V14),size=0.5)+
		geom_jitter(aes(c2$V2,2,col=c2$V14),size=0.5)+
		geom_jitter(aes(c3$V2,3,col=c3$V14),size=0.5)+
		geom_jitter(aes(c4$V2,4,col=c4$V14),size=0.5)+
		geom_jitter(aes(c5$V2,5,col=c5$V14),size=0.5)+
		geom_jitter(aes(cX$V2,6,col=cX$V14),size=0.5)+
		scale_colour_gradient2()+
		scale_y_continuous(breaks=seq(1,6,1),trans = 'reverse')+
		ggtitle("Adult Female Gut")+xlab("Genome position (bp)")+ylab("Normalised expression")+
		theme_bw()


#--- make plot
egg_plot +
l1_plot +
shl3_plot +
exl3_plot +
l4_plot +
adultm_plot +
adultf_plot +
gut_plot +
plot_layout(ncol=2)






# GO term analysis of cluster_stats

```R
R-3.5.0
library(topReviGO)
genes<-c(t(read.table("~/lustre118_link/hc/GENOME/GO_ANALYSIS/all_genes",header=F)))
mapfile<-read.table("HCON_V4_GOterm.db")


cluster1<-c(t(read.table("clust_cluster1.genelist",header=F)))
c1_allGenes <-  factor(as.integer((genes) %in% cluster1))
names(c1_allGenes) <- (genes)

topReviGO(c1_allGenes, "cluster1_MF", "/nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/GO_ANALYSIS/HCON_V4_GOterm.db", mapOrDb = "map", p = 0.05,ontology="MF")
topReviGO(c1_allGenes, "cluster1_BP", "/nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/GO_ANALYSIS/HCON_V4_GOterm.db", mapOrDb = "map", p = 0.05,ontology="BP")


cluster8<-c(t(read.table("clust_cluster8.genelist",header=F)))
c8_allGenes <-  factor(as.integer((genes) %in% cluster8))
names(c8_allGenes) <- (genes)

topReviGO(c8_allGenes, "cluster8_MF", "/nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/GO_ANALYSIS/HCON_V4_GOterm.db", mapOrDb = "map", p = 0.05,ontology="MF")
topReviGO(c8_allGenes, "cluster8_BP", "/nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/GO_ANALYSIS/HCON_V4_GOterm.db", mapOrDb = "map", p = 0.05,ontology="BP")




# plot cluster 8
revigo.names <- c("term_ID","description","frequency_%","plot_X","plot_Y","plot_size","log10_p_value","uniqueness","dispensability");
revigo.data <- rbind(c("GO:0004175","endopeptidase activity", 1.918,-0.793, 6.971, 5.431,-8.2757,0.713,0.000),
c("GO:0005328","neurotransmitter:sodium symporter activity", 0.057, 5.352, 4.778, 3.907,-2.5528,0.795,0.000),
c("GO:0019211","phosphatase activator activity", 0.013,-5.142,-3.744, 3.247,-1.3143,0.892,0.000),
c("GO:0030246","carbohydrate binding", 0.723, 5.887,-4.605, 5.007,-2.0605,0.885,0.000),
c("GO:0004329","formate-tetrahydrofolate ligase activity", 0.023,-6.798, 0.072, 3.510,-1.3143,0.866,0.023),
c("GO:0004609","phosphatidylserine decarboxylase activity", 0.041,-3.990,-0.303, 3.757,-1.3143,0.866,0.024),
c("GO:0008410","CoA-transferase activity", 0.076, 2.166,-6.048, 4.030,-2.1675,0.839,0.026),
c("GO:0000150","recombinase activity", 0.105,-5.806, 3.755, 4.171,-1.3143,0.864,0.026),
c("GO:0016646","oxidoreductase activity, acting on the CH-NH group of donors, NAD or NADP as acceptor", 0.239, 7.170, 0.594, 4.527,-1.3778,0.815,0.029),
c("GO:0016874","ligase activity", 3.540, 1.894,-2.047, 5.698,-1.6757,0.859,0.039),
c("GO:0035639","purine ribonucleoside triphosphate binding",15.815,-0.700, 1.487, 6.348,-2.5528,0.901,0.077),
c("GO:0003831","beta-N-acetylglucosaminylglycopeptide beta-1,4-galactosyltransferase activity", 0.001,-2.123,-5.397, 2.053,-1.3143,0.848,0.130),
c("GO:0004001","adenosine kinase activity", 0.008, 0.275,-6.320, 3.060,-1.3143,0.844,0.148),
c("GO:0004067","asparaginase activity", 0.022,-2.704, 6.251, 3.483,-1.3143,0.820,0.188),
c("GO:0004365","glyceraldehyde-3-phosphate dehydrogenase (NAD+) (phosphorylating) activity", 0.014, 6.875,-0.518, 3.282,-1.3143,0.824,0.240),
c("GO:0005337","nucleoside transmembrane transporter activity", 0.064, 4.925, 5.319, 3.953,-1.3778,0.811,0.398),
c("GO:0008121","ubiquinol-cytochrome-c reductase activity", 0.043, 5.942, 3.133, 3.782,-1.3143,0.728,0.554),
c("GO:0004185","serine-type carboxypeptidase activity", 0.143,-0.249, 7.172, 4.303,-1.6478,0.736,0.611));

one.data <- data.frame(revigo.data);
names(one.data) <- revigo.names;
one.data <- one.data [(one.data$plot_X != "null" & one.data$plot_Y != "null"), ];
one.data$plot_X <- as.numeric( as.character(one.data$plot_X) );
one.data$plot_Y <- as.numeric( as.character(one.data$plot_Y) );
one.data$plot_size <- as.numeric( as.character(one.data$plot_size) );
one.data$log10_p_value <- as.numeric( as.character(one.data$log10_p_value) );
one.data$frequency <- as.numeric( as.character(one.data$frequency) );
one.data$uniqueness <- as.numeric( as.character(one.data$uniqueness) );
one.data$dispensability <- as.numeric( as.character(one.data$dispensability) );
#head(one.data);

# --------------------------------------------------------------------------
# Names of the axes, sizes of the numbers and letters, names of the columns,
# etc. can be changed below

p1 <- ggplot( data = one.data );
p1 <- p1 + geom_point( aes( plot_X, plot_Y, colour = log10_p_value, size = plot_size), alpha = I(0.6) ) + scale_size_area();
p1 <- p1 + scale_colour_gradientn( colours = c("blue", "green", "yellow", "red"), limits = c( min(one.data$log10_p_value), 0) );
p1 <- p1 + geom_point( aes(plot_X, plot_Y, size = plot_size), shape = 21, fill = "transparent", colour = I (alpha ("black", 0.6) )) + scale_size_area();
p1 <- p1 + scale_size( range=c(5, 25)) + theme_bw(); # + scale_fill_gradientn(colours = heat_hcl(7), limits = c(-300, 0) );
ex <- one.data [ one.data$dispensability < 0.15, ];
p1 <- p1 + geom_label_repel( data = ex, aes(plot_X, plot_Y, label = description), colour = I(alpha("black", 0.85)), size = 2.5 ,force=2);
p1 <- p1 + labs (y = "semantic space x", x = "semantic space y");
p1 <- p1 + theme(legend.key = element_blank()) ;
one.x_range = max(one.data$plot_X) - min(one.data$plot_X);
one.y_range = max(one.data$plot_Y) - min(one.data$plot_Y);
p1 <- p1 + xlim(min(one.data$plot_X)-one.x_range/10,max(one.data$plot_X)+one.x_range/10);
p1 <- p1 + ylim(min(one.data$plot_Y)-one.y_range/10,max(one.data$plot_Y)+one.y_range/10);

# --------------------------------------------------------------------------
# Output the plot to screen
p1;





# cluster 1
revigo.names <- c("term_ID","description","frequency_%","plot_X","plot_Y","plot_size","log10_p_value","uniqueness","dispensability");
revigo.data <- rbind(c("GO:0000977","RNA polymerase II regulatory region sequence-specific DNA binding", 0.160,-5.206,-4.503, 4.352,-2.4089,0.783,0.000),
c("GO:0000981","RNA polymerase II transcription factor activity, sequence-specific DNA binding", 0.640, 6.209, 0.321, 4.954,-2.1487,0.791,0.000),
c("GO:0004888","transmembrane signaling receptor activity", 0.962,-0.366,-6.201, 5.132,-10.0969,0.792,0.000),
c("GO:0015075","ion transmembrane transporter activity", 3.726,-5.688, 1.839, 5.720,-2.6778,0.473,0.000),
c("GO:0016413","O-acetyltransferase activity", 0.066, 3.105, 4.825, 3.966,-2.0969,0.761,0.000),
c("GO:0019825","oxygen binding", 0.057, 4.144,-4.291, 3.904,-1.5498,0.783,0.037),
c("GO:0015018","galactosylgalactosylxylosylprotein 3-beta-glucuronosyltransferase activity", 0.009, 0.218, 6.442, 3.086,-1.6459,0.761,0.148),
c("GO:0005344","oxygen transporter activity", 0.033,-5.524, 2.740, 3.662,-1.6271,0.548,0.501));

one.data <- data.frame(revigo.data);
names(one.data) <- revigo.names;
one.data <- one.data [(one.data$plot_X != "null" & one.data$plot_Y != "null"), ];
one.data$plot_X <- as.numeric( as.character(one.data$plot_X) );
one.data$plot_Y <- as.numeric( as.character(one.data$plot_Y) );
one.data$plot_size <- as.numeric( as.character(one.data$plot_size) );
one.data$log10_p_value <- as.numeric( as.character(one.data$log10_p_value) );
one.data$frequency <- as.numeric( as.character(one.data$frequency) );
one.data$uniqueness <- as.numeric( as.character(one.data$uniqueness) );
one.data$dispensability <- as.numeric( as.character(one.data$dispensability) );
#head(one.data);

# --------------------------------------------------------------------------
# Names of the axes, sizes of the numbers and letters, names of the columns,
# etc. can be changed below

p1 <- ggplot( data = one.data );
p1 <- p1 + geom_point( aes( plot_X, plot_Y, colour = log10_p_value, size = plot_size), alpha = I(0.6) ) + scale_size_area();
p1 <- p1 + scale_colour_gradientn( colours = c("blue", "green", "yellow", "red"), limits = c( min(one.data$log10_p_value), 0) );
p1 <- p1 + geom_point( aes(plot_X, plot_Y, size = plot_size), shape = 21, fill = "transparent", colour = I (alpha ("black", 0.6) )) + scale_size_area();
p1 <- p1 + scale_size( range=c(5, 25)) + theme_bw(); # + scale_fill_gradientn(colours = heat_hcl(7), limits = c(-300, 0) );
ex <- one.data [ one.data$dispensability < 0.15, ];
p1 <- p1 + geom_label_repel( data = ex, aes(plot_X, plot_Y, label = description), colour = I(alpha("black", 0.85)), size = 2.5 ,force=2);
p1 <- p1 + labs (y = "semantic space x", x = "semantic space y");
p1 <- p1 + theme(legend.key = element_blank()) ;
one.x_range = max(one.data$plot_X) - min(one.data$plot_X);
one.y_range = max(one.data$plot_Y) - min(one.data$plot_Y);
p1 <- p1 + xlim(min(one.data$plot_X)-one.x_range/10,max(one.data$plot_X)+one.x_range/10);
p1 <- p1 + ylim(min(one.data$plot_Y)-one.y_range/10,max(one.data$plot_Y)+one.y_range/10);

# --------------------------------------------------------------------------
# Output the plot to screen
p1;
```

```shell
scp sd21@pcs5.internal.sanger.ac.uk:/nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/GO_ANALYSIS/revigo_cluster.* ~/Documents/workbook/hcontortus_genome/04_analysis
```





























---
## 05 - Splice Leader analyses


Setup Working directory and reference sequences
```shell

cd /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME
mkdir SPLICE_LEADERS
cd SPLICE_LEADERS

ln -s ../../REF/HAEM_V4_final.chr.fa REF.fa
ln -s ../TRANSCRIPTOME_CURATION/HCON_V4_WBP11plus_190125.ips.gff3 ANNOTATION.gff3

```

Slice leader sequences
SL1 - GGTTTAATTACCCAAGTTTGAG
SL2 - GGTTTTAACCCAGTATCTCAAG





Run the splice leader analysis tools

```shell
# run splice leader finder scripts for SL1 and SL2 sequences
bsub.py --queue yesterday 1 merged_AdultF_SL1 ~sd21/bash_scripts/run_spliceleader_finder.sh merged_AdultF_SL1 $PWD/REF.fa $PWD/ANNOTATION.gff3 GGTTTAATTACCCAAGTTTGAG 10 /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/TRANSCRIPTOME/RAW/merged_AdultF_R1.fastq.gz /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/TRANSCRIPTOME/RAW/merged_AdultF_R2.fastq.gz
bsub.py --queue yesterday 1 merged_AdultM_SL1 ~sd21/bash_scripts/run_spliceleader_finder.sh merged_AdultM_SL1 $PWD/REF.fa $PWD/ANNOTATION.gff3 GGTTTAATTACCCAAGTTTGAG 10 /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/TRANSCRIPTOME/RAW/merged_AdultM_R1.fastq.gz /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/TRANSCRIPTOME/RAW/merged_AdultM_R2.fastq.gz
bsub.py --queue yesterday 1 merged_BXC_SL1 ~sd21/bash_scripts/run_spliceleader_finder.sh merged_BXC_SL1 $PWD/REF.fa $PWD/ANNOTATION.gff3 GGTTTAATTACCCAAGTTTGAG 10 /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/TRANSCRIPTOME/RAW/merged_BXC_R1.fastq.gz /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/TRANSCRIPTOME/RAW/merged_BXC_R2.fastq.gz
bsub.py --queue yesterday 1 merged_EXL3_SL1 ~sd21/bash_scripts/run_spliceleader_finder.sh merged_EXL3_SL1 $PWD/REF.fa $PWD/ANNOTATION.gff3 GGTTTAATTACCCAAGTTTGAG 10 /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/TRANSCRIPTOME/RAW/merged_EXL3_R1.fastq.gz /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/TRANSCRIPTOME/RAW/merged_EXL3_R2.fastq.gz
bsub.py --queue yesterday 1 merged_gut_SL1 ~sd21/bash_scripts/run_spliceleader_finder.sh merged_gut_SL1 $PWD/REF.fa $PWD/ANNOTATION.gff3 GGTTTAATTACCCAAGTTTGAG 10 /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/TRANSCRIPTOME/RAW/merged_gut_R1.fastq.gz /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/TRANSCRIPTOME/RAW/merged_gut_R2.fastq.gz
bsub.py --queue yesterday 1 merged_ISE_BXWF_SL1 ~sd21/bash_scripts/run_spliceleader_finder.sh merged_ISE_BXWF_SL1 $PWD/REF.fa $PWD/ANNOTATION.gff3 GGTTTAATTACCCAAGTTTGAG 10 /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/TRANSCRIPTOME/RAW/merged_ISE_BXWF_R1.fastq.gz /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/TRANSCRIPTOME/RAW/merged_ISE_BXWF_R2.fastq.gz
bsub.py --queue yesterday 1 merged_L1_SL1 ~sd21/bash_scripts/run_spliceleader_finder.sh merged_L1_SL1 $PWD/REF.fa $PWD/ANNOTATION.gff3 GGTTTAATTACCCAAGTTTGAG 10 /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/TRANSCRIPTOME/RAW/merged_L1_R1.fastq.gz /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/TRANSCRIPTOME/RAW/merged_L1_R2.fastq.gz
bsub.py --queue yesterday 1 merged_L4_SL1 ~sd21/bash_scripts/run_spliceleader_finder.sh merged_L4_SL1 $PWD/REF.fa $PWD/ANNOTATION.gff3 GGTTTAATTACCCAAGTTTGAG 10 /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/TRANSCRIPTOME/RAW/merged_L4_R1.fastq.gz /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/TRANSCRIPTOME/RAW/merged_L4_R2.fastq.gz
bsub.py --queue yesterday 1 merged_SHL3_SL1 ~sd21/bash_scripts/run_spliceleader_finder.sh merged_SHL3_SL1 $PWD/REF.fa $PWD/ANNOTATION.gff3 GGTTTAATTACCCAAGTTTGAG 10 /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/TRANSCRIPTOME/RAW/merged_SHL3_R1.fastq.gz /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/TRANSCRIPTOME/RAW/merged_SHL3_R2.fastq.gz


bsub.py --queue yesterday 1 merged_AdultF_SL2 ~sd21/bash_scripts/run_spliceleader_finder.sh merged_AdultF_SL2 $PWD/REF.fa $PWD/ANNOTATION.gff3 GGTTTTAACCCAGTATCTCAAG 10 /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/TRANSCRIPTOME/RAW/merged_AdultF_R1.fastq.gz /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/TRANSCRIPTOME/RAW/merged_AdultF_R2.fastq.gz
bsub.py --queue yesterday 1 merged_AdultM_SL2 ~sd21/bash_scripts/run_spliceleader_finder.sh merged_AdultM_SL2 $PWD/REF.fa $PWD/ANNOTATION.gff3 GGTTTTAACCCAGTATCTCAAG 10 /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/TRANSCRIPTOME/RAW/merged_AdultM_R1.fastq.gz /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/TRANSCRIPTOME/RAW/merged_AdultM_R2.fastq.gz
bsub.py --queue yesterday 1 merged_BXC_SL2 ~sd21/bash_scripts/run_spliceleader_finder.sh merged_BXC_SL2 $PWD/REF.fa $PWD/ANNOTATION.gff3 GGTTTTAACCCAGTATCTCAAG 10 /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/TRANSCRIPTOME/RAW/merged_BXC_R1.fastq.gz /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/TRANSCRIPTOME/RAW/merged_BXC_R2.fastq.gz
bsub.py --queue yesterday 1 merged_EXL3_SL2 ~sd21/bash_scripts/run_spliceleader_finder.sh merged_EXL3_SL2 $PWD/REF.fa $PWD/ANNOTATION.gff3 GGTTTTAACCCAGTATCTCAAG 10 /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/TRANSCRIPTOME/RAW/merged_EXL3_R1.fastq.gz /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/TRANSCRIPTOME/RAW/merged_EXL3_R2.fastq.gz
bsub.py --queue yesterday 1 merged_gut_SL2 ~sd21/bash_scripts/run_spliceleader_finder.sh merged_gut_SL2 $PWD/REF.fa $PWD/ANNOTATION.gff3 GGTTTTAACCCAGTATCTCAAG 10 /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/TRANSCRIPTOME/RAW/merged_gut_R1.fastq.gz /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/TRANSCRIPTOME/RAW/merged_gut_R2.fastq.gz
bsub.py --queue yesterday 1 merged_ISE_BXWF_SL2 ~sd21/bash_scripts/run_spliceleader_finder.sh merged_ISE_BXWF_SL2 $PWD/REF.fa $PWD/ANNOTATION.gff3 GGTTTTAACCCAGTATCTCAAG 10 /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/TRANSCRIPTOME/RAW/merged_ISE_BXWF_R1.fastq.gz /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/TRANSCRIPTOME/RAW/merged_ISE_BXWF_R2.fastq.gz
bsub.py --queue yesterday 1 merged_L1_SL2 ~sd21/bash_scripts/run_spliceleader_finder.sh merged_L1_SL2 $PWD/REF.fa $PWD/ANNOTATION.gff3 GGTTTTAACCCAGTATCTCAAG 10 /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/TRANSCRIPTOME/RAW/merged_L1_R1.fastq.gz /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/TRANSCRIPTOME/RAW/merged_L1_R2.fastq.gz
bsub.py --queue yesterday 1 merged_L4_SL2 ~sd21/bash_scripts/run_spliceleader_finder.sh merged_L4_SL2 $PWD/REF.fa $PWD/ANNOTATION.gff3 GGTTTTAACCCAGTATCTCAAG 10 /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/TRANSCRIPTOME/RAW/merged_L4_R1.fastq.gz /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/TRANSCRIPTOME/RAW/merged_L4_R2.fastq.gz
bsub.py --queue yesterday 1 merged_SHL3_SL2 ~sd21/bash_scripts/run_spliceleader_finder.sh merged_SHL3_SL2 $PWD/REF.fa $PWD/ANNOTATION.gff3 GGTTTTAACCCAGTATCTCAAG 10 /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/TRANSCRIPTOME/RAW/merged_SHL3_R1.fastq.gz /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/TRANSCRIPTOME/RAW/merged_SHL3_R2.fastq.gz
```

```shell
awk '{if($6=="+" || $6=="-") print $0}' ${PREFIX}_transcript_start_windows.tmp.bed | sed 's/-[0-9][0-9]*/1/' transcript_start_windows.bed
```


Run subsequence analyses to collate and summarise data
```shell
# collate all splice window coordinates
cat *SL1_SL_ANALYSIS_out/*SL-only.bed > SL1.SL-only.bed
cat *SL2_SL_ANALYSIS_out/*SL-only.bed > SL2.SL-only.bed

# get the transcript start and internal window coordinates - these are identical in each of the sub-analyses, so just choose one set
cp merged_EXL3_SL1_SL_ANALYSIS_out/merged_EXL3_SL1_transcript_start_windows.bed transcript_start_windows.bed
cp merged_EXL3_SL1_SL_ANALYSIS_out/merged_EXL3_SL1_internal_cds_windows.bed internal_cds_windows.bed

# get SL coverage per transcript start for both SL1 and SL2 sequences
bedtools coverage -s -b transcript_start_windows.bed -a SL1.SL-only.bed | sort -k1,1 -k2,2n > SL1_transcript_start_windows.SL.coverage
bedtools coverage -s -b transcript_start_windows.bed -a SL2.SL-only.bed | sort -k1,1 -k2,2n > SL2_transcript_start_windows.SL.coverage


# get SL coverage per internal CDSs (skipping terminal CDSs) for both SL1 and SL2 sequences
bedtools coverage -s -b internal_cds_windows.bed -a SL1.SL-only.bed | sort -k1,1 -k2,2n > SL1_internal_cds_windows.SL.coverage
bedtools coverage -s -b internal_cds_windows.bed -a SL2.SL-only.bed | sort -k1,1 -k2,2n > SL2_internal_cds_windows.SL.coverage


awk '{if($7>=1) print $4}' SL1_transcript_start_windows.SL.coverage | sed -e 's/;/\t/g' -e 's/ID=//g' | cut -f1 > transcript_start_windows.SL1.genelist
awk '{if($7>=1) print $4}' SL1_transcript_start_windows.SL.coverage | wc -l
#>> 7895 transcripts
cut -f1 -d "-" transcript_start_windows.SL1.genelist | sort | uniq | wc -l
#>> 7293 genes



# awk '{if($7>=1) print $4}' SL1_internal_cds_windows.SL.coverage | sed -e 's/;/\t/g' -e 's/ID=//g' | cut -f1 > internal_cds_windows.SL1.genelist
# awk '{if($7>=1) print $4}' SL1_internal_cds_windows.SL.coverage | wc -l
# #>> 5933 transcripts
# cut -f1 -d "-" internal_cds_windows.SL1.genelist | sort | uniq | wc -l
# #>> 3824 genes



awk '{if($7>=1) print $4}' SL2_transcript_start_windows.SL.coverage | sed -e 's/;/\t/g' -e 's/ID=//g' | cut -f1 > transcript_start_windows.SL2.genelist
awk '{if($7>=1) print $4}' SL2_transcript_start_windows.SL.coverage | wc -l
#>> 1222 transcripts
cut -f1 -d "-" transcript_start_windows.SL2.genelist | sort | uniq | wc -l
#>> 1139 genes

# awk '{if($7>=1) print $4}' SL2_internal_cds_windows.SL.coverage | sed -e 's/;/\t/g' -e 's/ID=//g' | cut -f1 > internal_cds_windows.SL2.genelist
# awk '{if($7>=1) print $4}' SL2_internal_cds_windows.SL.coverage | wc -l
# #>> 950 transcripts
# cut -f1 -d "-" internal_cds_windows.SL2.genelist | sort | uniq | wc -l
# #>> 766 genes


# #>> 3824 genes
```


# comparion SL counts per gene, to check coverage of genes that have both SL1 and SL2 sequences in reads
paste SL1_transcript_start_windows.SL.coverage SL2_transcript_start_windows.SL.coverage > test
awk '{if($7>0 || $17>0) print $7,$17}' OFS="\t" test > SL1v2_comparison.counts

R.3.5.0
a<-read.table("SL1v2_comparison.counts",header=F,sep="\t")
ggplot(a,aes(log10(V1),log10(V2)))+geom_hex(bins = 30)+labs(x="Number of SL1 reads [log10(count)]", y="Number of SL2 reads [log10(count)]")
ggplot(a,aes((V1+0.001),(V2+0.001)))+geom_hex(bins = 30)+labs(x="Number of SL1 reads [log10(count)]", y="Number of SL2 reads [log10(count)]")+scale_y_log10()+scale_x_log10()
##########################################################################################
# distance and density of genes, and their relationship between SL1 and SL2
```shell


awk -F '[\t=;]' '{if($3=="mRNA" && $7=="+") print $1,$4,$5,$10,$6,$7}' OFS="\t" ANNOTATION.gff3 | sort -k1,1 -k2,2n > transcripts.pos_strand.bed
awk -F '[\t=;]' '{if($3=="mRNA" && $7=="-") print $1,$4,$5,$10,$6,$7}' OFS="\t" ANNOTATION.gff3 | sort -k1,1 -k2,2n > transcripts.neg_strand.bed


bedtools-2 spacing -i transcripts.pos_strand.bed | sed 's/\.$/NA/g' > transcripts.pos_strand.spacing.bed
bedtools-2 spacing -i transcripts.neg_strand.bed | sed 's/\.$/NA/g' > transcripts.neg_strand.spacing.bed

awk '{print $4,$7}' OFS="\t" transcripts.pos_strand.spacing.bed > transcripts.pos_strand.spacing
awk '{print $4,$7}' OFS="\t" transcripts.neg_strand.spacing.bed > transcripts.neg_strand.spacing

# fix neg strand coords to correct for orientation
cut -f 1 transcripts.neg_strand.spacing | head -n -1 > neg.1.tmp
cut -f 2 transcripts.neg_strand.spacing | tail -n +2 > neg.2.tmp
paste neg.1.tmp neg.2.tmp > transcripts.neg_strand.spacing


awk '{print $1,"SL1"}' OFS="\t" transcript_start_windows.SL1.genelist > mRNA_SL1.list
awk '{print $1,"SL2"}' OFS="\t" transcript_start_windows.SL2.genelist > mRNA_SL2.list

cat mRNA_SL1.list mRNA_SL2.list > mRNA_SLall.list
```


# hybrid SL1/SL2 shared sites
bedtools-2 intersect -s -a SL1.SL-only.bed -b SL2.SL-only.bed > SL1SL2_sharedsites.bed

bedtools coverage -s -b transcript_start_windows.bed -a SL1SL2_sharedsites.bed | sort -k1,1 -k2,2n > SL1SL2_sharedsites_transcript_start_windows.SL.coverage

awk '{if($7>=1) print $4}' SL1SL2_sharedsites_transcript_start_windows.SL.coverage | sed -e 's/;/\t/g' -e 's/ID=//g' | cut -f1 > SL1SL2_sharedsites.genelist
awk '{if($7>=1) print $4}' SL1SL2_sharedsites_transcript_start_windows.SL.coverage | wc -l
# #>> 1021 transcripts
cut -f1 -d "-" SL1SL2_sharedsites.genelist | sort | uniq | wc -l
# #>> 1019 genes



```shell
cd /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/TRANSCRIPTOME/SPLICE_LEADERS
```


```R
R-3.5.0
library(ggplot2)
library(patchwork)
library("ggsci")

pos<-read.table("transcripts.pos_strand.spacing",header=F)
neg<-read.table("transcripts.neg_strand.spacing",header=F)
sl<-read.table("mRNA_SLall.list",header=F)

data_pos<-dplyr::left_join(pos,sl,by="V1")
data_neg<-dplyr::left_join(neg,sl,by="V1")

sl1_pos<-data_pos[(data_pos$V2.y=="SL1"),]
sl2_pos<-data_pos[(data_pos$V2.y=="SL2"),]

sl1_neg<-data_neg[(data_neg$V2.y=="SL1"),]
sl2_neg<-data_neg[(data_neg$V2.y=="SL2"),]


par(mfrow=c(2,2))
plot(density(log10(sl1_pos$V2.x),na.rm=T),xlim=c(0,6))
plot(density(log10(sl1_neg$V2.x),na.rm=T),xlim=c(0,6))
plot(density(log10(sl2_pos$V2.x),na.rm=T),xlim=c(0,6))
plot(density(log10(sl2_neg$V2.x),na.rm=T),xlim=c(0,6))


plot_pos_all <- ggplot()+geom_histogram(aes(log10(data_pos$V2.x),fill=data_pos$V2.y),bins = 100, alpha = 0.5,position="identity")+theme_bw()+xlab("Distance to upstream gene (log10(bp))")+ylab("Count")+scale_fill_npg(na.value="grey90")+xlim(0,6)+ylim(0,500)
plot_neg_all <- ggplot()+geom_histogram(aes(log10(data_neg$V2.x),fill=data_neg$V2.y),bins = 100, alpha = 0.5,position="identity")+theme_bw()+xlab("Distance to upstream gene (log10(bp))")+ylab("Count")+scale_fill_npg(na.value="grey90")+xlim(0,6)+ylim(0,500)
plot_pos <- ggplot()+geom_histogram(aes(log10(sl1_pos$V2.x),stat(ndensity)),bins = 100,fill = "#E64B35B2", alpha = 0.5)+geom_histogram(aes(log10(sl2_pos$V2.x),stat(ndensity)),bins = 100,fill = "#4DBBD5B2", alpha = 0.5)+theme_bw()+xlab("Distance to upstream gene (log10(bp))")+ylab("Density")+xlim(0,6)
plot_neg <- ggplot()+geom_histogram(aes(log10(sl1_neg$V2.x),stat(ndensity)),bins = 100,fill = "#E64B35B2", alpha = 0.5)+geom_histogram(aes(log10(sl2_neg$V2.x),stat(ndensity)),bins = 100,fill = "#4DBBD5B2", alpha = 0.5)+theme_bw()+xlab("Distance to upstream gene (log10(bp))")+ylab("Density")+xlim(0,6)


# combined pos and neg
data_all <- dplyr::bind_rows(data_pos,data_neg)
sl1_all<-data_all[(data_all$V2.y=="SL1"),]
sl2_all<-data_all[(data_all$V2.y=="SL2"),]


p1 <- ggplot()+
     geom_histogram(aes(log10(data_all$V2.x),fill=(data_all$V2.y)),bins=100, alpha = 0.7,position="identity")+
     theme_bw()+
     xlab("Distance to upstream gene (log10(bp))")+ylab("Count")+scale_fill_npg(na.value="grey80")+
     xlim(0,6)+ylim(0,500)

p2 <- ggplot()+
     geom_histogram(aes(log10(sl1_all$V2.x),stat(ndensity)),bins = 100,fill = "#E64B35B2", alpha = 0.5)+
     geom_histogram(aes(log10(sl2_all$V2.x),stat(ndensity)),bins = 100,fill = "#4DBBD5B2", alpha=0.5)+
     theme_bw()+xlab("Distance to upstream gene (log10[bp])")+
     guides(fill = guide_legend(reverse = TRUE))+
     ylab("Density")+xlim(0,6)

p1 + p2 + plot_layout(ncol=2)
ggsave("SL_gene_distance.pdf")

```

```shell
scp sd21@pcs5.internal.sanger.ac.uk:/nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/TRANSCRIPTOME/SPLICE_LEADERS/SL_gene_distance.pdf /Users/sd21/Documents/workbook/hcontortus_genome/04_analysis




### weblogo plots

```shell
*_adapters.fixedlength.fa
*_SL1.twentymer.fasta

cat *SL1_SL_ANALYSIS_out/*_adapters.fixedlength.fa > SL1_adapters_fixedlength.fa
cat *SL2_SL_ANALYSIS_out/*_adapters.fixedlength.fa > SL2_adapters_fixedlength.fa

~sd21/lustre118_link/software/bin/weblogo --format pdf --datatype fasta --ignore-lower-case <  SL1_adapters_fixedlength.fa >  SL1_adapters_fixedlength.weblogo.pdf
~sd21/lustre118_link/software/bin/weblogo --format pdf --datatype fasta --ignore-lower-case <  SL2_adapters_fixedlength.fa >  SL2_adapters_fixedlength.weblogo.pdf

cat *SL1_SL_ANALYSIS_out/*_SL1.twentymer.fasta > SL1_splicesite_20mer.fa
cat *SL2_SL_ANALYSIS_out/*_SL2.twentymer.fasta > SL2_splicesite_20mer.fa

~sd21/lustre118_link/software/bin/weblogo --format pdf --datatype fasta --ignore-lower-case <  SL1_splicesite_20mer.fa >  SL1_splicesite_20mer.weblogo.pdf
~sd21/lustre118_link/software/bin/weblogo --format pdf --datatype fasta --ignore-lower-case <  SL2_splicesite_20mer.fa >  SL2_splicesite_20mer.weblogo.pdf
```



---
## 05 - Differential splicing w Leafcutter <a name="ds_leafcutter"></a>
---

LEAFCUTTER
following http://davidaknowles.github.io/leafcutter/articles/Usage.html


##### NOTE - made modification to leafcutter scripts
# modified the following scripts
# --  ~sd21/lustre118_link/software/TRANSCRIPTOME/leafcutter/scripts/leafcutter_ds.R
# --  ~sd21/lustre118_link/software/TRANSCRIPTOME/leafcutter/leafviz/prepare_results.R
# -- ~sd21/lustre118_link/software/TRANSCRIPTOME/leafcutter/scripts/ds_plots.R
# -- ~sd21/lustre118_link/software/TRANSCRIPTOME/leafcutter/leafviz/run_leafviz.R

# changed the header from:
#  #!/usr/bin/env Rscript
#
# to
#
# #!/usr/bin/env /software/R-3.5.0/bin/Rscript
#
# to ensure the correct Rscript version was being used.


```shell
export PATH=$PATH:~sd21/lustre118_link/software/TRANSCRIPTOME/leafcutter/scripts/

# raw data: /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/TRANSCRIPTOME/RAW

cd /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/TRANSCRIPTOME/LEAFCUTTER

ln -sf /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/TRANSCRIPTOME/TRANSCRIPTOME_CURATION/HCON_V4_WBP11plus_190125.ips.gff3 ANNOTATION.gff3
ln -s /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/REF/HAEM_V4_final.chr.fa REF.fa
```

Need to make a gtf of the annotation and curate it for leafcutter
```shell
# make gtf from gff
gffread -T -g REF.fa -o ANNOTATION.gtf ANNOTATION.gff3

# use leafcutter script to curate
perl ~sd21/lustre118_link/software/TRANSCRIPTOME/leafcutter/leafviz/gtf2leafcutter.pl -o HCON_V4 ANNOTATION.gtf
```





```shell
# map RNAseq reads to reference using STAR
 while read name lane; do \
 		/nfs/users/nfs_s/sd21/lustre118_link/software/TRANSCRIPTOME/STAR/bin/Linux_x86_64/STAR \
 		--twopassMode Basic \
 		--runThreadN 8 \
 		--genomeDir /nfs/users/nfs_s/sd21/lustre118_link/REFERENCE_SEQUENCES/haemonchus_contortus/hc_v4chr_rnaseq_star_index \
 		--readFilesIn ../RAW/${lane}_1.fastq.gz ../RAW/${lane}_2.fastq.gz \
 		--readFilesCommand zcat \
 		--outSAMstrandField intronMotif \
 		--outSAMtype BAM Unsorted \
 		--outFileNamePrefix ${name}; \
 done < lanes.list



# make junc files
mkdir LEAFCUTTER_JUNCS

for bamfile in `ls STAR_MAPPING/*bam | awk -F "/" '{print $NF}' `; do \
		echo Converting $PWD/STAR_MAPPING/$bamfile to $PWD/LEAFCUTTER_JUNCS/$bamfile.junc; \
		sh bam2junc.sh $PWD/STAR_MAPPING/$bamfile $PWD/LEAFCUTTER_JUNCS/$bamfile.junc ; \
		echo $PWD/STAR_MAPPING/$bamfile.junc >> $PWD/LEAFCUTTER_JUNCS/test_juncfiles.txt; \
done



# make exon list

echo "chr start end strand gene_name" > HCON_V4.renamed.exons.leafcutter.list; \
     awk -F'[\t=;]' '{if($3=="exon") print $1,$4,$5,$7,$12}'  ANNOTATION.gff3 | sort -k1,1 -k 2,2n >> HCON_V4.renamed.exons.leafcutter.list
gzip HCON_V4.renamed.exons.leafcutter.list
```


Run full analysis
```shell
# make comparison directories
mkdir LC_EGG_L1
mkdir LC_L1_SHL3
mkdir LC_EXL3_L4
mkdir LC_L4_AdultF
mkdir LC_EXL3_L4
mkdir LC_L4_AdultF
mkdir LC_L4_AdultM
mkdir LC_GUT_AdultF
mkdir LC_FEMALE_EGG
mkdir LC_SHL3_EXL3
mkdir LC_ADULTM_ADULTF


# copy junc files to each comparison dirctory
cp LEAFCUTTER_JUNCS/AdultF*junc LC_L4_AdultF/
cp LEAFCUTTER_JUNCS/AdultF*junc LC_GUT_AdultF
cp LEAFCUTTER_JUNCS/AdultF*junc LC_FEMALE_EGG
cp LEAFCUTTER_JUNCS/AdultM*junc LC_L4_AdultM
cp LEAFCUTTER_JUNCS/Adult*junc LC_ADULTM_ADULTF
cp LEAFCUTTER_JUNCS/egg*junc LC_EGG_L1
cp LEAFCUTTER_JUNCS/egg*junc LC_FEMALE_EGG
cp LEAFCUTTER_JUNCS/L1*junc LC_EGG_L1
cp LEAFCUTTER_JUNCS/L1*junc LC_L1_SHL3
cp LEAFCUTTER_JUNCS/L4*junc LC_EXL3_L4
cp LEAFCUTTER_JUNCS/L4*junc LC_L4_AdultF
cp LEAFCUTTER_JUNCS/L4*junc LC_L4_AdultM/
cp LEAFCUTTER_JUNCS/gut*junc LC_GUT_AdultF

cp LEAFCUTTER_JUNCS/SHL3*junc LC_SHL3_EXL3
cp LEAFCUTTER_JUNCS/EXL3*junc LC_SHL3_EXL3

# fix for renaming of SHL3 and EXL3 as per kallisto analysis
cd LC_SHL3_EXL3
mv SHL3_3_236476_2305_3914314Aligned.out.bam.junc EXL3_kallisto_7062_6_12_out.junc
mv EXL3_2_236476_2305_3914310Aligned.out.bam.junc SHL3_kallisto_7062_6_8_out.junc
rm EXL3_1_236476_2305_3914309Aligned.out.bam.junc

# move correct samples to other directories
cp EXL3*junc ../LC_EXL3_L4
cp SHL3*junc ../LC_L1_SHL3

cd ../


# make list of junc files in each comparison directory
for i in LC_* ; do \
		ls -1  $i/*junc | \
		xargs -n1 basename > ${i}/juncfiles.list; \
		done

#


# cluster splice junctions
#--- min depth 5
for i in LC_* ; do \
     cd ${i} && \
     python2.7 ~sd21/lustre118_link/software/TRANSCRIPTOME/leafcutter/clustering/leafcutter_cluster.py -j juncfiles.list -m 5 -o ${i}_cov5 -l 500000 && \
     cd ../; \
     done

#--- min depth 10
for i in LC_* ; do \
     cd ${i} && \
     python2.7 ~sd21/lustre118_link/software/TRANSCRIPTOME/leafcutter/clustering/leafcutter_cluster.py -j juncfiles.list -m 10 -o ${i}_cov10 -l 500000 && \
     cd ../; \
     done

#--- min depth 30 (default)
for i in LC_* ; do \
          cd ${i} && \
          python2.7 ~sd21/lustre118_link/software/TRANSCRIPTOME/leafcutter/clustering/leafcutter_cluster.py -j juncfiles.list -m 30 -o ${i}_cov30 -l 500000 && \
          cd ../; \
          done



# all samples
mkdir LEAFCUTTER_JUNCS_V4
cp LC_*/*junc LEAFCUTTER_JUNCS_V4
cd LEAFCUTTER_JUNCS_V4

ls -1 *junc > juncfiles.list

python2.7 ~sd21/lustre118_link/software/TRANSCRIPTOME/leafcutter/clustering/leafcutter_cluster.py -j juncfiles.list -m 5 -o all_samples_cov5 -l 500000

python2.7 ~sd21/lustre118_link/software/TRANSCRIPTOME/leafcutter/clustering/leafcutter_cluster.py -j juncfiles.list -m 10 -o all_samples_cov10 -l 500000

python2.7 ~sd21/lustre118_link/software/TRANSCRIPTOME/leafcutter/clustering/leafcutter_cluster.py -j juncfiles.list -m 30 -o all_samples_cov30 -l 500000


cd ..





# perform differential analysis
#--- made "samples_groups.list" files for each comparison manually, as the sample names were not really parsable.
#--- min depth 5
for i in LC_FEMALE_EGG LC_GUT_AdultF LC_L1_SHL3 LC_L4_AdultF LC_L4_AdultM LC_ADULTM_ADULTF LC_EGG_L1; do \
     		cd ${i} && \
     		~sd21/lustre118_link/software/TRANSCRIPTOME/leafcutter/scripts/leafcutter_ds.R \
     		--min_samples_per_intron=3 --min_samples_per_group=3 --min_coverage=10 --exon_file=../HCON_V4_all_exons.txt.gz ${i}_cov5_perind_numers.counts.gz samples_groups.list --output_prefix=${i}_cov5 && \
     		cd ../ ; \
     		done

for i in LC_SHL3_EXL3 LC_EXL3_L4 ; do \
               cd ${i} && \
               ~sd21/lustre118_link/software/TRANSCRIPTOME/leafcutter/scripts/leafcutter_ds.R \
               --min_samples_per_intron=2 --min_samples_per_group=2 --min_coverage=10 --exon_file=../HCON_V4_all_exons.txt.gz ${i}_cov5_perind_numers.counts.gz samples_groups.list --output_prefix=${i}_cov5 && \
               cd ../ ; \
               done



#--- min depth 10
for i in LC_FEMALE_EGG LC_GUT_AdultF LC_L1_SHL3 LC_L4_AdultF LC_L4_AdultM LC_ADULTM_ADULTF LC_EGG_L1; do \
          cd ${i} && \
          ~sd21/lustre118_link/software/TRANSCRIPTOME/leafcutter/scripts/leafcutter_ds.R \
          --min_samples_per_intron=3 --min_samples_per_group=3 --min_coverage=10 --exon_file=../HCON_V4_all_exons.txt.gz ${i}_cov10_perind_numers.counts.gz samples_groups.list --output_prefix=${i}_cov10 && \
          cd ../ ; \
          done

for i in LC_SHL3_EXL3 LC_EXL3_L4 ; do \
          cd ${i} && \
          ~sd21/lustre118_link/software/TRANSCRIPTOME/leafcutter/scripts/leafcutter_ds.R \
          --min_samples_per_intron=2 --min_samples_per_group=2 -exon_file=../HCON_V4_all_exons.txt.gz ${i}_cov10_perind_numers.counts.gz samples_groups.list --output_prefix=${i}_cov10 && \
          cd ../ ; \
          done


#--- min depth 30
for i in LC_FEMALE_EGG LC_GUT_AdultF LC_L1_SHL3 LC_L4_AdultF LC_L4_AdultM LC_ADULTM_ADULTF LC_EGG_L1; do \
                    cd ${i} && \
                    ~sd21/lustre118_link/software/TRANSCRIPTOME/leafcutter/scripts/leafcutter_ds.R \
                    --min_samples_per_intron=3 --min_samples_per_group=3 --exon_file=../HCON_V4_all_exons.txt.gz ${i}_cov30_perind_numers.counts.gz samples_groups.list --output_prefix=${i}_cov30 && \
                    cd ../ ; \
                    done

for i in LC_SHL3_EXL3 LC_EXL3_L4 ; do \
                    cd ${i} && \
                    ~sd21/lustre118_link/software/TRANSCRIPTOME/leafcutter/scripts/leafcutter_ds.R \
                    --min_samples_per_intron=2 --min_samples_per_group=2 --exon_file=../HCON_V4_all_exons.txt.gz ${i}_cov30_perind_numers.counts.gz samples_groups.list --output_prefix=${i}_cov30 && \
                    cd ../ ; \
                    done


#--- made "samples_groups.list" files for each comparison manually, as the sample names were not really parsable
for i in LC*; do cd $i; ls -1 *junc | awk -F'[_]' '{print $0,$1}' OFS="\t" > samples_groups.list; cd ../ ;  done


# prepare results for visualisation
#--- min depth 5
for i in LC_* ; do \
		cd ${i} && \
		~sd21/lustre118_link/software/TRANSCRIPTOME/leafcutter/leafviz/prepare_results.R ${i}_cov5_perind_numers.counts.gz ${i}_cov5_cluster_significance.txt ${i}_cov5_effect_sizes.txt ../HCON_V4 --meta_data_file samples_groups.list --output=${i}_leafviz_min5.RData && \
		cd ../ ; \
		done

#--- min depth 10
for i in LC_* ; do \
		cd ${i} && \
		~sd21/lustre118_link/software/TRANSCRIPTOME/leafcutter/leafviz/prepare_results.R ${i}_cov10_perind_numers.counts.gz ${i}_cov5_cluster_significance.txt ${i}_cov5_effect_sizes.txt ../HCON_V4 --meta_data_file samples_groups.list --output=${i}_leafviz_min10.RData && \
		cd ../ ; \
		done



#--- min depth 30
for i in LC_* ; do \
          cd ${i} && \
          ~sd21/lustre118_link/software/TRANSCRIPTOME/leafcutter/leafviz/prepare_results.R ${i}_cov30_perind_numers.counts.gz ${i}_cov30_cluster_significance.txt ${i}_cov30_effect_sizes.txt ../HCON_V4 --meta_data_file samples_groups.list --output=${i}_leafviz_min30.RData && \
          cd ../ ; \
          done


# copy data to laptop for visualisation

for i in LC_* ; do cp ${i}/${i}_leafviz_min30.RData ~sd21/desktop_link/${i}_leafviz_min30.RData; done


#- on laptop, move data here: /Users/sd21/Documents/Work/aaa_projects/haemonchus_contortus/project_hc_genome/data/Section_4/leafcutter_data
cp /Volumes/sd21/desktop_link/*.RData .

# obviosuly needs leafcutter installed locally, which is ehre: http://davidaknowles.github.io/leafcutter/articles/Installation.html

# also copied

# to visualise on laptop
cd /Users/sd21/Documents/Work/aaa_projects/haemonchus_contortus/project_hc_genome/data/Section_4/leafcutter_data
cp ~/Documents/bioinformatics_bin/leafcutter/leafviz/server.R .
cp ~/Documents/bioinformatics_bin/leafcutter/leafviz/ui.R .

~/Documents/bioinformatics_bin/leafcutter/leafviz/run_leafviz.R LF_L4_AdultM_leafviz.RData



########################################################################################
# Meta analysis - compare differential intro usage across all replicates for all lifestages
########################################################################################
cd /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/TRANSCRIPTOME/LEAFCUTTER
mkdir METAANALYSIS_V4

# get sig clusters from  cluster_significance.txt  - top 50 based of FDR pvalue
for i in LC_* ; do \
          cd ${i} && \
          awk -F'[:\t]' '{if($7<0.05) print $2,$7}' OFS="\t" ${i}_cov30_cluster_significance.txt | sort -k2 -g | head -n 100 > ../METAANALYSIS_V4/${i}_sigclusters.list && \
          cd ../ ; \
          done


# extract introns from clusters that are signficant from effect_sizes.txt

rm METAANALYSIS_V4/sigclusters_introncoords.list

for i in LC_* ; do \
          cd ${i} && \
          while read cluster sig; do \
          grep $cluster ${i}_cov5_effect_sizes.txt; done < ../METAANALYSIS_V4/${i}_sigclusters.list  | awk -F'[:\t]' '{if($4>=5 || $4<=-5) print $1,$2,$3}' OFS=":" >> ../METAANALYSIS_V4/sigclusters_introncoords.list && \
          cd ../ ; \
          done

cd METAANALYSIS_V4
cat sigclusters_introncoords.list | sort -k1n | uniq > sigclusters_introncoords.sorted.list


# extract counts using cluster, intro coords

while read coords; do grep $coords <(zcat ../LEAFCUTTER_JUNCS_V4/all_samples_cov30_perind.counts.gz) ; done < sigclusters_introncoords.sorted.list > sigclusters_introncoords.junccounts


# ---- Sample columns in all_samples_cov30_perind.counts.gz

eggs1_236476_3881079Aligned.out.bam 20
eggs2_236476_3881080Aligned.out.bam 19
eggs3_236476_3881081Aligned.out.bam 14

L1_1_236476_3881082Aligned.out.bam 24
L1_2_236476_3881083Aligned.out.bam 10
L1_3_236476_3881084Aligned.out.bam 4

SHL3_1_236476_2305_3914312Aligned.out.bam 6
SHL3_2_236476_2305_3914313Aligned.out.bam 12
SHL3_kallisto_7062_6_8_out 22

EXL3_3_236476_2305_3914311Aligned.out.bam 5
EXL3_kallisto_7062_6_12_out 18

L4_1_236476_3881085Aligned.out.bam 17
L4_2_236476_3881086Aligned.out.bam 21
L4_3_236476_3881087Aligned.out.bam 11


AdultF1_236476_3881088Aligned.out.bam 2
AdultF2_236476_3881089Aligned.out.bam 23
AdultF3_236476_3881090Aligned.out.bam 8

AdultM1_236476_2305_3914303Aligned.out.bam 7
AdultM2_236476_2305_3914304Aligned.out.bam 9
AdultM3_236476_2305_3914305Aligned.out.bam 3

gut3_236476_635J_3914317Aligned.out.bam 13
gut2_236476_1589_3914316Aligned.out.bam 15
gut1_236476_1517_3914315Aligned.out.bam 16




# reorder the columns to fit the haem lifecycle, and then calc the splice frequency
awk '{print $20,$19,$14,$24,$10,$4,$6,$12,$22,$5,$18,$17,$21,$11,$2,$23,$8,$7,$9,$3,$13,$15,$16}' sigclusters_introncoords.junccounts |\
	sed 's/\//\t/g' |\
	awk '{print $1/($2+=0.0001),$3/($4+=0.0001),$5/($6+=0.0001),$7/($8+=0.0001),$9/($10+=0.0001),$11/($12+=0.0001),$13/($14+=0.0001),$15/($16+=0.0001),$17/($18+=0.0001),$19/($20+=0.0001),$21/($22+=0.0001),$23/($24+=0.0001),$25/($26+=0.0001),$27/($28+=0.0001),$29/($30+=0.0001),$31/($32+=0.0001),$33/($34+=0.0001),$35/($36+=0.0001),$37/($38+=0.0001),$39/($40+=0.0001),$41/($42+=0.0001),$43/($44+=0.0001),$45/($46+=0.0001)}' OFS="\t" > sigclusters_introncoords.juncfreq
	awk '{if(($1+$2+$3+$4+$5+$6+$7+$8+$9+$10+$11+$12+$13+$14+$15+$16+$17+$18+$19+$20+$21+$22+$23)>=0.5) print $0}' OFS="\t" sigclusters_introncoords.juncfreq > sigclusters_introncoords.juncfreq2


# R commands to make heatmap
R-3.5.0
pdf("allsamples.leafcutter_ds.top100.pdf")
library(gplots)
library(RColorBrewer)
a<-read.table("sigclusters_introncoords.juncfreq2",header=F)
colnames(a)<- c("Egg_1","Egg_2","Egg_3","L1_1","L1_2","L1_3","SHL3_1","SHL3_2","SHL3_3","EXL3_1","EXL3_2","L4_1","L4_2","L4_3","Adult_Female_1","Adult_Female_2","Adult_Female_3","Adult_Male_1","Adult_Male_2","Adult_Male_3","Gut_1","Gut_2","Gut_3")
heatmap.2(as.matrix(a),Colv=NULL,trace="none",col= colorRampPalette(brewer.pal(8, "Blues"))(25),labRow=FALSE,margins = c(10, 10))
heatmap.2(as.matrix(a),Colv=NULL,trace="none",col= colorRampPalette(brewer.pal(8, "Blues"))(25),labRow=FALSE,labCol=FALSE,margins = c(10, 10),dendrogram='none',density.info='none')
dev.off()
quit(save = "no")

cd ../


########################################################################################
# repeat extract significant clusters for each pairwise comparison


cd /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/TRANSCRIPTOME/LEAFCUTTER


# get sig clusters from  cluster_significance.txt  - top 100 based of FDR pvalue
# --  done in the meta analysis step above

# extract intron coords per sample, sort
for i in LC_* ; do \
		cd ${i} && \
		awk -F'[:\t]' '{if($6<0.05) print $2,$7,$8}' OFS="\t" ${i}_cov30_cluster_significance.txt | sort -k2 -g | head -n 100 > ${i}_sigclusters.list
		while read cluster sig mrna_id; do grep $cluster ${i}_cov30_effect_sizes.txt; done < ${i}_sigclusters.list  | awk -F'[:\t]' '{if($4>=5 || $4<=-5) print $1,$2,$3}' OFS=":" >  ${i}_sigclusters_introncoords.list
		cat ${i}_sigclusters_introncoords.list | sort -k1n | uniq > ${i}_sigclusters_introncoords.sorted.list
		while read coords; do grep $coords <(zcat ${i}_cov30_perind.counts.gz) ; done < ${i}_sigclusters_introncoords.sorted.list > ${i}_sigclusters_introncoords.junccounts
		cd ../ ; \
		done

# eggs to L1
cd LC_EGG_L1
awk '{print $6,$5,$4,$3,$2,$7,$1}' LC_EGG_L1_sigclusters_introncoords.junccounts |\
sed -e 's/\//\t/g' -e 's/:/\t/g' |\
awk '{print $16,$1/($2+=0.0001),$3/($4+=0.0001),$5/($6+=0.0001),$7/($8+=0.0001),$9/($10+=0.0001),$11/($12+=0.0001)}' OFS="\t" > LC_EGG_L1_sigclusters_introncoords.juncfreq
while read cluster sig mrna; do sed -i "s/${cluster}/${mrna}/g" LC_EGG_L1_sigclusters_introncoords.juncfreq; done < LC_EGG_L1_sigclusters.list
awk '{if(($2+$3+$4+$5+$6+$7)>=0.5) print $2,$3,$4,$5,$6,$7}' OFS="\t" LC_EGG_L1_sigclusters_introncoords.juncfreq > LC_EGG_L1_sigclusters_introncoords.juncfreq2


R-3.5.0
pdf("LC_EGG_L1.leafcutter_ds.top100.pdf")
library(gplots)
library(RColorBrewer)
a<-read.table("LC_EGG_L1_sigclusters_introncoords.juncfreq2",header=F)
colnames(a)<- c("Egg_1","Egg_2","Egg_3","L1_1","L1_2","L1_3")
#heatmap.2(as.matrix(a),Colv=NULL,trace="none",col= colorRampPalette(brewer.pal(8, "Blues"))(25),labRow=FALSE,margins = c(10, 10))
heatmap.2(as.matrix(a),Colv=NULL,trace="none",col= colorRampPalette(brewer.pal(8, "Blues"))(25),labRow=FALSE,labCol=FALSE,margins = c(10, 10),dendrogram='none',density.info='none')
dev.off()
quit(save = "no")

cd ..



# L1 to SHL3
cd LC_L1_SHL3
awk '{print $5,$2,$3,$6,$4,$7,$1}' LC_L1_SHL3_sigclusters_introncoords.junccounts |\
sed -e 's/\//\t/g' -e 's/:/\t/g' |\
awk '{print $16,$1/($2+=0.0001),$3/($4+=0.0001),$5/($6+=0.0001),$7/($8+=0.0001),$9/($10+=0.0001),$11/($12+=0.0001)}' OFS="\t" > LC_L1_SHL3_sigclusters_introncoords.juncfreq
while read cluster sig mrna; do sed -i "s/${cluster}/${mrna}/g" LC_L1_SHL3_sigclusters_introncoords.juncfreq ; done < LC_L1_SHL3_sigclusters.list
awk '{if(($2+$3+$4+$5+$6+$7)>=0.5) print $2,$3,$4,$5,$6,$7}' OFS="\t" LC_L1_SHL3_sigclusters_introncoords.juncfreq > LC_L1_SHL3_sigclusters_introncoords.juncfreq2


R-3.5.0
pdf("LC_L1_SHL3.leafcutter_ds.top100.pdf")
library(gplots)
library(RColorBrewer)
a<-read.table("LC_L1_SHL3_sigclusters_introncoords.juncfreq2",header=F)
colnames(a)<- c("L1_1","L1_2","L1_3","SHL3_1","SHL3_2","SHL3_3")
#heatmap.2(as.matrix(a),Colv=NULL,trace="none",col= colorRampPalette(brewer.pal(8, "Blues"))(25),labRow=FALSE,margins = c(10, 10))
heatmap.2(as.matrix(a),Colv=NULL,trace="none",col= colorRampPalette(brewer.pal(8, "Blues"))(25),labRow=FALSE,labCol=FALSE,margins = c(10, 10),dendrogram='none',density.info='none')
dev.off()
quit(save = "no")

cd ..



# SHL3 to EXL3
cd LC_SHL3_EXL3
awk '{print $2,$3,$6,$4,$5,$1}' LC_SHL3_EXL3_sigclusters_introncoords.junccounts |\
sed -e 's/\//\t/g' -e 's/:/\t/g' |\
awk '{print $14,$1/($2+=0.0001),$3/($4+=0.0001),$5/($6+=0.0001),$7/($8+=0.0001),$9/($10+=0.0001)}' OFS="\t" > LC_SHL3_EXL3_sigclusters_introncoords.juncfreq
while read cluster sig mrna; do sed -i "s/${cluster}/${mrna}/g" LC_SHL3_EXL3_sigclusters_introncoords.juncfreq ; done < LC_SHL3_EXL3_sigclusters.list
awk '{if(($2+$3+$4+$5+$6)>=0.5) print $2,$3,$4,$5,$6}' OFS="\t" LC_SHL3_EXL3_sigclusters_introncoords.juncfreq > LC_SHL3_EXL3_sigclusters_introncoords.juncfreq2

R-3.5.0
pdf("LC_SHL3_EXL3.leafcutter_ds.top100.pdf")
library(gplots)
library(RColorBrewer)
a<-read.table("LC_SHL3_EXL3_sigclusters_introncoords.juncfreq2",header=F)
colnames(a)<- c("SHL3_1","SHL3_2","SHL3_3","EXL3_1","EXL3_2")
#heatmap.2(as.matrix(a),Colv=NULL,trace="none",col= colorRampPalette(brewer.pal(8, "Blues"))(25),labRow=FALSE,margins = c(10, 10))
heatmap.2(as.matrix(a),Colv=NULL,trace="none",col= colorRampPalette(brewer.pal(8, "Blues"))(25),labRow=FALSE,labCol=FALSE,margins = c(10, 10),dendrogram='none',density.info='none')
dev.off()
quit(save = "no")

cd ../


# EXL3 to L4
cd LC_EXL3_L4
awk '{print $4,$5,$2,$6,$3,$1}' LC_EXL3_L4_sigclusters_introncoords.junccounts |\
sed -e 's/\//\t/g' -e 's/:/\t/g' |\
awk '{print $14,$1/($2+=0.0001),$3/($4+=0.0001),$5/($6+=0.0001),$7/($8+=0.0001),$9/($10+=0.0001)}' OFS="\t" > LC_EXL3_L4_sigclusters_introncoords.juncfreq
while read cluster sig mrna; do sed -i "s/${cluster}/${mrna}/g" LC_EXL3_L4_sigclusters_introncoords.juncfreq ; done < LC_EXL3_L4_sigclusters.list
awk '{if(($2+$3+$4+$5+$6)>=0.5) print $2,$3,$4,$5,$6}' OFS="\t" LC_EXL3_L4_sigclusters_introncoords.juncfreq > LC_EXL3_L4_sigclusters_introncoords.juncfreq2


R-3.5.0
pdf("LC_EXL3_L4.leafcutter_ds.top100.pdf")
library(gplots)
library(RColorBrewer)
a<-read.table("LC_EXL3_L4_sigclusters_introncoords.juncfreq2",header=F)
colnames(a)<- c("EXL3_1","EXL3_2","L4_1","L4_2","L4_3")
#heatmap.2(as.matrix(a),Colv=NULL,trace="none",col= colorRampPalette(brewer.pal(8, "Blues"))(25),labRow=FALSE,margins = c(10, 10))
heatmap.2(as.matrix(a),Colv=NULL,trace="none",col= colorRampPalette(brewer.pal(8, "Blues"))(25),labRow=FALSE,labCol=FALSE,margins = c(10, 10),dendrogram='none',density.info='none')
dev.off()
quit(save = "no")
cd ../



# L4 to AdultMale
cd LC_L4_AdultM
awk '{print $2,$5,$3,$6,$7,$4,$1}' LC_L4_AdultM_sigclusters_introncoords.junccounts |\
sed -e 's/\//\t/g' -e 's/:/\t/g' |\
awk '{print $16,$1/($2+=0.0001),$3/($4+=0.0001),$5/($6+=0.0001),$7/($8+=0.0001),$9/($10+=0.0001),$11/($12+=0.0001)}' OFS="\t" > LC_L4_AdultM_sigclusters_introncoords.juncfreq
while read cluster sig mrna; do sed -i "s/${cluster}/${mrna}/g" LC_L4_AdultM_sigclusters_introncoords.juncfreq ; done < LC_L4_AdultM_sigclusters.list
awk '{if(($2+$3+$4+$5+$6+$7)>=0.5) print $2,$3,$4,$5,$6,$7}' OFS="\t" LC_L4_AdultM_sigclusters_introncoords.juncfreq > LC_L4_AdultM_sigclusters_introncoords.juncfreq2


R-3.5.0
pdf("LC_L4_AdultM.leafcutter_ds.top100.pdf")
library(gplots)
library(RColorBrewer)
a<-read.table("LC_L4_AdultM_sigclusters_introncoords.juncfreq2",header=F)
colnames(a)<- c("L4_1","L4_2","L4_3","Adult_Male_1","Adult_Male_2","Adult_Male_3")
#heatmap.2(as.matrix(a),Colv=NULL,trace="none",col= colorRampPalette(brewer.pal(8, "Blues"))(25),labRow=FALSE,margins = c(10, 10))
heatmap.2(as.matrix(a),Colv=NULL,trace="none",col= colorRampPalette(brewer.pal(8, "Blues"))(25),labRow=FALSE,labCol=FALSE,margins = c(10, 10),dendrogram='none',density.info='none')
dev.off()
quit(save = "no")

cd ../




# L4 to AdultFemale
cd LC_L4_AdultF
awk '{print $2,$6,$4,$3,$5,$7,$1}' LC_L4_AdultF_sigclusters_introncoords.junccounts |\
sed -e 's/\//\t/g' -e 's/:/\t/g' |\
awk '{print $16,$1/($2+=0.0001),$3/($4+=0.0001),$5/($6+=0.0001),$7/($8+=0.0001),$9/($10+=0.0001),$11/($12+=0.0001)}' OFS="\t" > LC_L4_AdultF_sigclusters_introncoords.juncfreq
while read cluster sig mrna; do sed -i "s/${cluster}/${mrna}/g" LC_L4_AdultF_sigclusters_introncoords.juncfreq ; done < LC_L4_AdultF_sigclusters.list
awk '{if(($2+$3+$4+$5+$6+$7)>=0.5) print $2,$3,$4,$5,$6,$7}' OFS="\t" LC_L4_AdultF_sigclusters_introncoords.juncfreq > LC_L4_AdultF_sigclusters_introncoords.juncfreq2


R-3.5.0
pdf("LC_L4_AdultF.leafcutter_ds.top100.pdf")
library(gplots)
library(RColorBrewer)
a<-read.table("LC_L4_AdultF_sigclusters_introncoords.juncfreq2",header=F)
colnames(a)<- c("L4_1","L4_2","L4_3","Adult_Female_1","Adult_Female_2","Adult_Female_3")
#heatmap.2(as.matrix(a),Colv=NULL,trace="none",col= colorRampPalette(brewer.pal(8, "Blues"))(25),labRow=FALSE,margins = c(10, 10))		
heatmap.2(as.matrix(a),Colv=NULL,trace="none",col= colorRampPalette(brewer.pal(8, "Blues"))(25),labRow=FALSE,labCol=FALSE,margins = c(10, 10),dendrogram='none',density.info='none')
dev.off()
quit(save = "no")

cd ../



# AdultMale to AdultFemale
cd LC_ADULTM_ADULTF
awk '{print $5,$7,$3,$2,$4,$6,$1}' LC_ADULTM_ADULTF_sigclusters_introncoords.junccounts |\
     sed -e 's/\//\t/g' -e 's/:/\t/g' |\
     awk '{print $16,$1/($2+=0.0001),$3/($4+=0.0001),$5/($6+=0.0001),$7/($8+=0.0001),$9/($10+=0.0001),$11/($12+=0.0001)}' OFS="\t" > LC_ADULTM_ADULTF_sigclusters_introncoords.juncfreq
while read cluster sig mrna; do sed -i "s/${cluster}/${mrna}/g" LC_ADULTM_ADULTF_sigclusters_introncoords.juncfreq ; done < LC_ADULTM_ADULTF_sigclusters.list
awk '{if(($2+$3+$4+$5+$6+$7)>=0.5) print $2,$3,$4,$5,$6,$7}' OFS="\t" LC_ADULTM_ADULTF_sigclusters_introncoords.juncfreq > LC_ADULTM_ADULTF_sigclusters_introncoords.juncfreq2

while read cluster sig mrna; do sed -i "s/${cluster}/${mrna}_${cluster}/g" LC_ADULTM_ADULTF_sigclusters_introncoords.juncfreq ; done < LC_ADULTM_ADULTF_sigclusters.list
awk '{if(($2+$3+$4+$5+$6+$7)>=0.5) print NR"_"$1,$2,$3,$4,$5,$6,$7}' OFS="\t" LC_ADULTM_ADULTF_sigclusters_introncoords.juncfreq > LC_ADULTM_ADULTF_sigclusters_introncoords.juncfreq2




R-3.5.0
pdf("LC_ADULTM_ADULTF.leafcutter_ds.top100.pdf")
library(gplots)
library(RColorBrewer)
a<-read.table("LC_ADULTM_ADULTF_sigclusters_introncoords.juncfreq2",header=F,row.names=1)
colnames(a)<- c("Adult_Male_1","Adult_Male_2","Adult_Male_3","Adult_Female_1","Adult_Female_2","Adult_Female_3")
#heatmap.2(as.matrix(a),Colv=NULL,trace="none",col= colorRampPalette(brewer.pal(8, "Blues"))(25),labRow=FALSE,margins = c(10, 10))
heatmap.2(as.matrix(a),Colv=NULL,trace="none",col= colorRampPalette(brewer.pal(8, "Blues"))(25),labRow=FALSE,labCol=FALSE,margins = c(10, 10),dendrogram='none',density.info='none')
dev.off()
quit(save = "no")
cd ../

# AdultFemale to eggs
cd LC_FEMALE_EGG
awk '{print $2,$3,$7,$6,$5,$4,$1}' LC_FEMALE_EGG_sigclusters_introncoords.junccounts |\
sed -e 's/\//\t/g' -e 's/:/\t/g' |\
awk '{print $16,$1/($2+=0.0001),$3/($4+=0.0001),$5/($6+=0.0001),$7/($8+=0.0001),$9/($10+=0.0001),$11/($12+=0.0001)}' OFS="\t" > LC_FEMALE_EGG_sigclusters_introncoords.juncfreq
while read cluster sig mrna; do sed -i "s/${cluster}/${mrna}/g" LC_FEMALE_EGG_sigclusters_introncoords.juncfreq ; done < LC_FEMALE_EGG_sigclusters.list
awk '{if(($2+$3+$4+$5+$6+$7)>=0.5) print $2,$3,$4,$5,$6,$7}' OFS="\t" LC_FEMALE_EGG_sigclusters_introncoords.juncfreq > LC_FEMALE_EGG_sigclusters_introncoords.juncfreq2





R-3.5.0
pdf("LC_FEMALE_EGG.leafcutter_ds.top100.pdf")
library(gplots)
library(RColorBrewer)
a<-read.table("LC_FEMALE_EGG_sigclusters_introncoords.juncfreq2",header=F)
colnames(a)<- c("Adult_Female_1","Adult_Female_2","Adult_Female_3","Egg_1","Egg_2","Egg_3")
#heatmap.2(as.matrix(a),Colv=NULL,trace="none",col= colorRampPalette(brewer.pal(8, "Blues"))(25),labRow=FALSE,margins = c(10, 10))
heatmap.2(as.matrix(a),Colv=NULL,trace="none",col= colorRampPalette(brewer.pal(8, "Blues"))(25),labRow=FALSE,labCol=FALSE,margins = c(10, 10),dendrogram='none',density.info='none')
dev.off()
quit(save = "no")
cd ../

# AdultFemale to gut
cd LC_GUT_AdultF
awk '{print $2,$3,$7,$6,$5,$4,$1}' LC_GUT_AdultF_sigclusters_introncoords.junccounts |\
sed -e 's/\//\t/g' -e 's/:/\t/g' |\
awk '{print $16,$1/($2+=0.0001),$3/($4+=0.0001),$5/($6+=0.0001),$7/($8+=0.0001),$9/($10+=0.0001),$11/($12+=0.0001)}' OFS="\t" > LC_GUT_AdultF_sigclusters_introncoords.juncfreq
while read cluster sig mrna; do sed -i "s/${cluster}/${mrna}/g" LC_GUT_AdultF_sigclusters_introncoords.juncfreq ; done < LC_GUT_AdultF_sigclusters.list
awk '{if(($2+$3+$4+$5+$6+$7)>=0.5) print $2,$3,$4,$5,$6,$7}' OFS="\t" LC_GUT_AdultF_sigclusters_introncoords.juncfreq > LC_GUT_AdultF_sigclusters_introncoords.juncfreq2


R-3.5.0
pdf("LC_GUT_AdultF.leafcutter_ds.top100.pdf",)
library(gplots)
library(RColorBrewer)
a<-read.table("LC_GUT_AdultF_sigclusters_introncoords.juncfreq2",header=F)
colnames(a)<- c("Adult_Female_1","Adult_Female_2","Adult_Female_3","Gut_1","Gut_2","Gut_3")
#heatmap.2(as.matrix(a),Colv=NULL,trace="none",col= colorRampPalette(brewer.pal(8, "Blues"))(25),labRow=FALSE,margins = c(10, 10))
heatmap.2(as.matrix(a),Colv=NULL,trace="none",col= colorRampPalette(brewer.pal(8, "Blues"))(25),labRow=FALSE,labCol=FALSE,margins = c(10, 10),dendrogram='none',density.info='none')
dev.off()
quit(save = "no")
cd ../

# calculated summary stats
for i in LC_*; do echo ${i}; \
# total genes
echo "total genes"; \
cut -f7 ${i}/*cov30_cluster_significance.txt | sort | uniq | wc -l ; \
echo "total significant genes"; \
awk '{if($6<0.05) print $7}' ${i}/*cov30_cluster_significance.txt | sort | uniq | wc -l ; \
# total introns
echo "total introns"; \
cat ${i}/*cov30_cluster_significance.txt | wc -l ; \
# total significant introns
echo "total significant introns"; \
awk '{if($6<0.05) print $0}' ${i}/*cov30_cluster_significance.txt | wc -l; \
echo " "
done


# calculate the number of genes in total impacted by differential splicing

cat LC_*/*cov30_cluster_significance.txt | awk '{if($6<0.05) print }' | cut -f7 | sed 's/,/\n/g' | sort | uniq | wc -l
1055






# Motif enrichment in 5' utr


dreme -o cluster_10_DREMEOUT -p cluster_10.single_transcript.5UTR.100bp.fa -n mrna.5UTR.100bp.fasta -eps && fimo -o cluster_10_DREMEOUT/FIMO cluster_10_DREMEOUT/dreme.txt cluster_10.single_transcript.5UTR.100bp.fa
dreme -o cluster_11_DREMEOUT -p cluster_11.single_transcript.5UTR.100bp.fa -n mrna.5UTR.100bp.fasta -eps && fimo -o cluster_11_DREMEOUT/FIMO cluster_11_DREMEOUT/dreme.txt cluster_11.single_transcript.5UTR.100bp.fa
dreme -o cluster_12_DREMEOUT -p cluster_12.single_transcript.5UTR.100bp.fa -n mrna.5UTR.100bp.fasta -eps && fimo -o cluster_12_DREMEOUT/FIMO cluster_12_DREMEOUT/dreme.txt cluster_12.single_transcript.5UTR.100bp.fa
dreme -o cluster_13_DREMEOUT -p cluster_13.single_transcript.5UTR.100bp.fa -n mrna.5UTR.100bp.fasta -eps && fimo -o cluster_13_DREMEOUT/FIMO cluster_13_DREMEOUT/dreme.txt cluster_13.single_transcript.5UTR.100bp.fa
dreme -o cluster_14_DREMEOUT -p cluster_14.single_transcript.5UTR.100bp.fa -n mrna.5UTR.100bp.fasta -eps && fimo -o cluster_14_DREMEOUT/FIMO cluster_14_DREMEOUT/dreme.txt cluster_14.single_transcript.5UTR.100bp.fa
dreme -o cluster_15_DREMEOUT -p cluster_15.single_transcript.5UTR.100bp.fa -n mrna.5UTR.100bp.fasta -eps && fimo -o cluster_15_DREMEOUT/FIMO cluster_15_DREMEOUT/dreme.txt cluster_15.single_transcript.5UTR.100bp.fa
dreme -o cluster_16_DREMEOUT -p cluster_16.single_transcript.5UTR.100bp.fa -n mrna.5UTR.100bp.fasta -eps && fimo -o cluster_16_DREMEOUT/FIMO cluster_16_DREMEOUT/dreme.txt cluster_16.single_transcript.5UTR.100bp.fa
dreme -o cluster_17_DREMEOUT -p cluster_17.single_transcript.5UTR.100bp.fa -n mrna.5UTR.100bp.fasta -eps && fimo -o cluster_17_DREMEOUT/FIMO cluster_17_DREMEOUT/dreme.txt cluster_17.single_transcript.5UTR.100bp.fa
dreme -o cluster_18_DREMEOUT -p cluster_18.single_transcript.5UTR.100bp.fa -n mrna.5UTR.100bp.fasta -eps && fimo -o cluster_18_DREMEOUT/FIMO cluster_18_DREMEOUT/dreme.txt cluster_18.single_transcript.5UTR.100bp.fa
dreme -o cluster_19_DREMEOUT -p cluster_19.single_transcript.5UTR.100bp.fa -n mrna.5UTR.100bp.fasta -eps && fimo -o cluster_19_DREMEOUT/FIMO cluster_19_DREMEOUT/dreme.txt cluster_19.single_transcript.5UTR.100bp.fa
dreme -o cluster_1_DREMEOUT -p cluster_1.single_transcript.5UTR.100bp.fa -n mrna.5UTR.100bp.fasta -eps && fimo -o cluster_1_DREMEOUT/FIMO cluster_1_DREMEOUT/dreme.txt cluster_1.single_transcript.5UTR.100bp.fa
dreme -o cluster_2_DREMEOUT -p cluster_2.single_transcript.5UTR.100bp.fa -n mrna.5UTR.100bp.fasta -eps && fimo -o cluster_2_DREMEOUT/FIMO cluster_2_DREMEOUT/dreme.txt cluster_2.single_transcript.5UTR.100bp.fa
dreme -o cluster_3_DREMEOUT -p cluster_3.single_transcript.5UTR.100bp.fa -n mrna.5UTR.100bp.fasta -eps && fimo -o cluster_3_DREMEOUT/FIMO cluster_3_DREMEOUT/dreme.txt cluster_3.single_transcript.5UTR.100bp.fa
dreme -o cluster_4_DREMEOUT -p cluster_4.single_transcript.5UTR.100bp.fa -n mrna.5UTR.100bp.fasta -eps && fimo -o cluster_4_DREMEOUT/FIMO cluster_4_DREMEOUT/dreme.txt cluster_4.single_transcript.5UTR.100bp.fa
dreme -o cluster_5_DREMEOUT -p cluster_5.single_transcript.5UTR.100bp.fa -n mrna.5UTR.100bp.fasta -eps && fimo -o cluster_5_DREMEOUT/FIMO cluster_5_DREMEOUT/dreme.txt cluster_5.single_transcript.5UTR.100bp.fa
dreme -o cluster_6_DREMEOUT -p cluster_6.single_transcript.5UTR.100bp.fa -n mrna.5UTR.100bp.fasta -eps && fimo -o cluster_6_DREMEOUT/FIMO cluster_6_DREMEOUT/dreme.txt cluster_6.single_transcript.5UTR.100bp.fa
dreme -o cluster_7_DREMEOUT -p cluster_7.single_transcript.5UTR.100bp.fa -n mrna.5UTR.100bp.fasta -eps && fimo -o cluster_7_DREMEOUT/FIMO cluster_7_DREMEOUT/dreme.txt cluster_7.single_transcript.5UTR.100bp.fa
dreme -o cluster_8_DREMEOUT -p cluster_8.single_transcript.5UTR.100bp.fa -n mrna.5UTR.100bp.fasta -eps && fimo -o cluster_8_DREMEOUT/FIMO cluster_8_DREMEOUT/dreme.txt cluster_8.single_transcript.5UTR.100bp.fa
dreme -o cluster_9_DREMEOUT -p cluster_9.single_transcript.5UTR.100bp.fa -n mrna.5UTR.100bp.fasta -eps && fimo -o cluster_9_DREMEOUT/FIMO cluster_9_DREMEOUT/dreme.txt cluster_9.single_transcript.5UTR.100bp.fa






for i in *.single_transcript.list ; do for j in $( cut -c-13 ${i} | sort | uniq ); do MOTIF=$( grep "${j}" ${i%.single_transcript.list}_DREMEOUT/FIMO/fimo.tsv | cut -f1 | sort | uniq | awk 'BEGIN { ORS = "," } { print }'); echo -e "${j}\t${MOTIF}" ; done >> ${i%.single_transcript.list}.genes.motifs ; done
