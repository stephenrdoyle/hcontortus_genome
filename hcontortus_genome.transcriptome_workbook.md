# Haemonchus genome - transcriptome analyses

## Table of contents

1.
2.
3. [Transcriptome QC](#transcriptome_qc)
4. [Kallisto](#kallisto)
5. [Differential splicing w Leafcutter](#ds_leafcutter)

















---

## 03 - Transcriptome QC <a name="transcriptome_qc"></a>

Date 190116

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
cp ~sd21/lustre118_link/hc/GENOME/TRANSCRIPTOME/TRANSCRIPTOME_CURATION/HCON_V4_WBP11plus_190114.gff3 HCON_V4_FINAL.gff3

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
        Base level:    85.5     |    79.7    |
        Exon level:    89.7     |    86.7    |
      Intron level:    94.8     |    93.9    |
Intron chain level:    67.5     |    61.5    |
  Transcript level:    72.7     |    56.3    |

       Locus level:    78.0     |    60.6    |


# V4 Final vs PASA_CHR
gffcompare -R -r HCON_V4_FINAL.gff3 -o V4_FINAL_vs_PASA HCON_V4_PASA.gff3

#-----------------| Sensitivity | Precision  |
        Base level:    85.0     |    97.6    |
        Exon level:    89.7     |    95.4    |
      Intron level:    90.1     |    96.5    |
Intron chain level:    33.3     |    32.9    |

  Transcript level:    33.2     |    30.8    |
       Locus level:    30.3     |    29.8    |


# V4 Final vs EVM

#-----------------| Sensitivity | Precision  |
        Base level:    92.9     |    98.7    |
        Exon level:    94.8     |    96.2    |
      Intron level:    96.4     |    98.1    |
Intron chain level:    84.2     |    84.3    |
  Transcript level:    86.8     |    86.8    |

       Locus level:    87.6     |    86.5    |

```

Don't think it is worth including the V1 comparison, as these curated genes would have been incorporated into the final annotation. Not really a good comparison. V1 precision is low due to only a subset of genes being used.

Results suggest:
- no increase in sensitivity, but big increase in precision from BRAKER to PASA
- increase in Sensitivity and Precision from PASA to EVM




##### HACK - didn't use this in the end, but keeping for reference
- There are some additionaly lines contaminating the final gff, and so while that is being fixed, need a workaround to extract sd21_modified genes from final gff3
- plan: grep relevant mRNA lines, and extract mRNA name (HCON etc) and ID (apollo unique ID), and use this as a list to grep out relevant lines back out from the original GFF

```shell  
grep "^hcon" HCON_V4_FINAL.gff3 | awk '$3=="mRNA" {print $0}' OFS="\t"  | cut -f 9 | sed -e 's/=/\t/g' -e 's/;/\t/g' | awk '{print $6,$10}' OFS="\n" > sd21_genes.list

while read NAME; do grep "=$NAME;" HCON_V4_FINAL.gff3; done  < sd21_genes.list >  sd21_genes.gff
sort sd21_genes.gff | uniq > sd21_genes.uniq.gff
```









---
## 04 - Kallisto <a name="kallisto"></a>
---








---
## 05 - Differential splicing w Leafcutter <a name="ds_leafcutter"></a>
---
