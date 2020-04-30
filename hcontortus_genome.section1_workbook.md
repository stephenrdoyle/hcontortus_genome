# Haemonchus genome paper
## Section 1: Chromosome structure of Haemonchus contortus

1. [Genome assembly](#genome)
2. [Completion of X chromosome](#xchromosome)
3. [Genome polishing](#polishing)
3. [Circos plot - Figure 1A](#circos)
4. [Orthology conservation and order per chromosome - Figure 1B top and bottom](#orthology)
5. [Microsynteny](#microsynteny)
6. [Genome stats](#genomestats)
7. [Genome completeness - CEGMA & BUSCO](#cegmabusco)
8. [Comparative analysis of the NZ Haemonchus genome](#nzgenome)



## Genome assembly
From the paper:

Initial manual improvement on the V1 genome focused on iterative scaffolding with SSPACE (Boetzer et al., 2011) and gap-filling with IMAGE (Tsai et al., 2010) using Illumina 500 bp and 3 kbp libraries, with additional low coverage data from 3, 8 and 20 kbp libraries generated using Roche 454 sequencing. These improvements were performed alongside breaking of discordant joins using Reapr (Hunt et al., 2013), and visual inspection using Gap5 (Bonfield and Whitwham, 2010). Substantial genetic variation was present in the sequence data due to the sequencing of DNA derived from a pool of individuals, resulting in a high frequency of haplotypes that assembled separately and therefore present as multiple copies of unique regions in the assembly. We surmised that much of the assembly fragmentation was due to the scaffolding tools not being able to deal with the changing rates of haplotypic variation so we attempted to solve this manually in gap5. We were aware that we did not have sufficient information to correctly phase these haplotypes, so instead, we chose the longest scaffold paths available, accepting that our scaffolds would not represent single haplotypes but would rather be an amalgamation of haplotypes representing single chromosomal regions. This approach was initially difficult and time-consuming, and was further confounded by a large number of repetitive sequences present in each haplotype.

Significant improvements in scaffold length were gained by the integration of OpGen (http://www.opgen.com/) optical mapping data. Optical mapping was performed following methods described previously (Tsai et al., 2013) with the following exceptions: high molecular weight DNA was prepared by proteinase K lysis of pools of ~500 H. contortus L3 embedded in agarose plugs, after which one of three restriction enzymes (KpnI, AflII and KpnI) were used to generate three separate restriction map datasets. Initial attempts to generate a de novo assembly using optical molecules alone was unsuccessful, and therefore, optical contigs were generated using DNA sequence seeds from the genome assembly using GenomeBuilder (OpGen) and visualised and manually curated using AssemblyViewer (OpGen). Although this approach was successful, it was limited by the quality and integrity of the gap-dense scaffolds and arbitrary nature of the haplotype scaffolding.

Subsequent integration of PacBio long-read data alongside the optical mapping data resulted in major increases in contiguity. PacBio sequencing libraries were prepared and sequenced using the PacBio RSII system. A total of 32.3 Gbp raw subreads (n = 4,085,541, N50 = 10,299 bp) were generated from 33 flow cells, representing approximately 133.8× coverage (based on the estimated genome size of 283 Mbp). A de novo PacBio assembly was generated using Sprai (v0.9.9.18; http://zombie.cb.k.u-tokyo.ac.jp/sprai/index.html), which was mapped to the assembly and used to manually resolve many gaps and resolve some of the phasing between haplotypes. Using the Sprai PacBio de novo assembly, we were also able to incorporate contigs that were previously missing from the working assembly. The increase in contiguity of the PacBio assemblies, further improved using canu v1.3 (Koren et al., 2017), revealed two major but diverse haplotype groups segregating at approximately 65% and 30% frequency in the pooled individuals sequenced; the presence of such diverse haplotypes resulted in a significantly expanded assembly over 500 Mbp. The major haplotype was more contiguous and therefore was chosen as the primary sequence to incorporate and improve the assembly. This approach was supported by competitive mapping of a single worm short read sequencing library (ENA: ERS086777), which was found to preferentially map to different scaffolds and/or different contigs within scaffolds that we inferred were different haplotypes in this single worm. Regions in the chromosome with no ERS086777 coverage, but for which an alternate haplotype was identified with coverage from the PacBio assembly, were manually removed from the main chromosomes and replaced. Further, this selective mapping correlated well with the optical contigs, and once these sequences were removed, much better optical alignments were produced further improving the assembly. Alternate haplotypes from the PacBio assembly, and those removed from the main assembly, were analysed as a separate haplotype assembly.

The increase in contiguity and resolution of shorter repetitive regions using PacBio began to reveal chromosome-specific repetitive units. Although these repeats were highly collapsed in the assembly and were typically not spanned by optical molecules, we were able to iteratively identify and join contigs/scaffolds flanking large tandem repeats that had clear read-pair evidence that they only occurred in a single location in the genome (i.e. read pairs were all mapped locally once the join was made). These were further supported by local assemblies of subsets of PacBio reads that contained copies of the repeat regions followed by de novo assembly using canu (Koren et al., 2017) to reconstruct the flanking unique regions surrounding a repeat. This iterative process resulted in the production of chromosome-scale scaffolds, each terminating with a 6 bp repeat consistent with being a telomeric sequence (sequence motif: TTAGGC).



## Completion of the X chromosome
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


#

working dir: /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/FIX_XCHR



# get barcodes 10X reads

barcoded.fastq.gz -> /lustre/scratch118/infgen/team133/sd21/hc/10X_genomics/chromium/hc_inbred/hc_inbred/outs/barcoded.fastq.gz

# note - to get the barcoded.fastq.gz, needed to run 10X longranger , eg.
bsub.py --queue yesterday --threads 4 10 hc_fq /nfs/users/nfs_s/sd21/lustre118_link/software/10X_GENOMICS/longranger-2.1.2/longranger basic --id=hc_inbred --fastqs=/nfs/users/nfs_s/sd21/lustre118_link/hc/10X_genomics/chromium/hc_inbred --readgroup=hc_inbred --localcores=4

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
