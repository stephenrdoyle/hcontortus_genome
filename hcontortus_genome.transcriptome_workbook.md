# Haemonchus genome - transcriptome analyses

## Table of contents

1.
2.
3. [Transcriptome QC](#transcriptome_qc)
4. [Kallisto](#kallisto)
5. [Differential splicing w Leafcutter](#ds_leafcutter)

















---

## 03 - Transcriptome QC <a name="transcriptome_qc"></a>

Want to compare the final annotation to steps along the way. These include
- V1 genome vs manually curated V1 genes
_ V4 final vs AUGUSTUS
- V4 final vs BRAKER
- V4 final vs PASA w Isoseq R1
- V4 final vs EVM w Isoseq R2

Need to do this with the 11 manually curated genes, and the whole V4 final.

### Working environment
```shell
cd /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/TRANSCRIPTOME
mkdir TRANSCRIPTOME_QC
cd TRANSCRIPTOME_QC
```


### Get some data
```shell
ln -s






gffcompare -R -r AUGUSTUSTraining.gtf -o jn_training_VS_breaker sd.augustus.gff3
```













---
## 04 - Kallisto <a name="kallisto"></a>
---








---
## 05 - Differential splicing w Leafcutter <a name="ds_leafcutter"></a>
---
