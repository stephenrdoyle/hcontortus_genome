# Haemonchus contortus genome paper
## Section 6: Distribution of global genetic diversity throughout the chromosomes

1. [mtDNA analysis](#mtdna_analysis)
2. [World map of sampling sites](#global_sampling_map)
3. [Genome wide nucleotide diversity analysis](#nuc_div)

# Get Raw data
```shell
cd /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/POPULATION_DIVERSITY/RAW

for i in 16693_1#1 16693_2#1 16720_1#1 16720_2#1; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 16693_1#2 16693_2#2 16720_1#2 16720_2#2; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 16693_1#3 16693_2#3 16720_1#3 16720_2#3; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 16693_1#4 16693_2#4 16720_1#4 16720_2#4; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 16693_1#5 16693_2#5 16720_1#5 16720_2#5; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 16693_1#6 16693_2#6 16720_1#6 16720_2#6; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 16693_1#7  16693_2#7  16720_1#7  16720_2#7; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 16693_1#8  16693_2#8  16720_1#8  16720_2#8; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 16693_1#9  16693_2#9  16720_1#9  16720_2#9; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 16693_1#10  16693_2#10  16720_1#10  16720_2#10; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 16693_1#11  16693_2#11  16720_1#11  16720_2#11; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 16693_1#12  16693_2#12  16720_1#12  16720_2#12; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 16693_1#13  16693_2#13  16720_1#13  16720_2#13; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 16693_1#14  16693_2#14  16720_1#14  16720_2#14; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 16693_1#15  16693_2#15  16720_1#15  16720_2#15; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 16693_1#16  16693_2#16  16720_1#16  16720_2#16; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 16693_1#17  16693_2#17  16720_1#17  16720_2#17; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 16693_1#18  16693_2#18  16720_1#18  16720_2#18; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 16693_1#19  16693_2#19  16720_1#19  16720_2#19; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 16693_1#20  16693_2#20  16720_1#20  16720_2#20; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 16693_1#21  16693_2#21  16720_1#21  16720_2#21; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 16693_1#22  16693_2#22  16720_1#22  16720_2#22; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 16693_1#23  16693_2#23  16720_1#23  16720_2#23; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 16693_1#24  16693_2#24  16720_1#24  16720_2#24; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 21766_7#1; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 21766_7#10; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 21766_7#11; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 21766_7#2; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 21766_7#3; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 21766_7#4; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 21766_7#5; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 21766_7#6; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 21766_7#7; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 21766_7#8; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 21766_7#9; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 20601_8#4; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 20601_8#5; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 20601_8#6; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 20601_8#7; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 20601_8#8; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 15150_1#34  15150_2#34  16553_1#34  16553_2#34; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 15150_1#35  15150_2#35  16553_1#35  16553_2#35; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 15150_1#36  15150_2#36  16553_1#36  16553_2#36; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 15150_1#19  15150_2#19  16553_1#19  16553_2#19; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 15150_1#20  15150_2#20  16553_1#20  16553_2#20; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 15150_1#21  15150_2#21  16553_1#21  16553_2#21; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 15150_1#22  15150_2#22  16553_1#22  16553_2#22; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 15150_1#23  15150_2#23  16553_1#23  16553_2#23; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 15150_1#24  15150_2#24  16553_1#24  16553_2#24; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 15150_1#25  15150_2#25  16553_1#25  16553_2#25; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 15150_1#26  15150_2#26  16553_1#26  16553_2#26; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 15150_1#27  15150_2#27  16553_1#27  16553_2#27; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 15150_1#28  15150_2#28  16553_1#28  16553_2#28; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 15150_1#29  15150_2#29  16553_1#29  16553_2#29; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 15150_1#30  15150_2#30  16553_1#30  16553_2#30; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 15150_1#31  15150_2#31  16553_1#31  16553_2#31; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 15150_1#32  15150_2#32  16553_1#32  16553_2#32; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 15150_1#33  15150_2#33  16553_1#33  16553_2#33; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 15045_1#1  15045_2#1  16552_1#1  16552_2#1; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 15045_1#2  15045_2#2  16552_1#2  16552_2#2; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 15045_1#3  15045_2#3  16552_1#3  16552_2#3; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 15045_1#4  15045_2#4  16552_1#4  16552_2#4; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 15045_1#5  15045_2#5  16552_1#5  16552_2#5; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 15045_1#6  15045_2#6  16552_1#6  16552_2#6; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 15045_1#7  15045_2#7  16552_1#7  16552_2#7; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 15045_1#8  15045_2#8  16552_1#8  16552_2#8; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 15045_1#9  15045_2#9  16552_1#9  16552_2#9; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 15045_1#10  15045_2#10  16552_1#10  16552_2#10; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 15045_1#11  15045_2#11  16552_1#11  16552_2#11; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 15045_1#12  15045_2#12  16552_1#12  16552_2#12; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 15045_1#13  15045_2#13  16552_1#13  16552_2#13; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 15045_1#14  15045_2#14  16552_1#14  16552_2#14; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 15045_1#15  15045_2#15  16552_1#15  16552_2#15; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
```


# mapping
```shell
#--- run_bwamem_splitter automatically submits jobs to LSF, so no need to bsub it. Will slowly work though the list when running screen
#--- samples_lanes.list was curated in excel (See: "Genome_paper_population_diversity_samples.xls"), but is a tab delimited text file containing sample name
 and sequencing lane ID for all samples to be mapped
# CH_SWI_001	16693_1#13
# CH_SWI_001	16693_2#13
# CH_SWI_001	16720_1#13
# CH_SWI_001	16720_2#13
# CH_SWI_002	16693_1#14

cd /lustre/scratch118/infgen/team133/sd21/hc/GENOME/POPULATION_DIVERSITY/MAPPING

echo "while read name lane; do \
~sd21/bash_scripts/run_bwamem_splitter $"{name}"_$"{lane}" /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/POPULATION_DIVERSITY/HAEM_V4_final.chr.fa \
/nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/POPULATION_DIVERSITY/RAW/$"{lane}"_1.fastq.gz \
/nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/POPULATION_DIVERSITY/RAW/$"{lane}"_2.fastq.gz; \
done < /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/POPULATION_DIVERSITY/samples_lanes.list" > run_mapping
chmod a+x run_mapping
screen
./run_mapping &



GB_IRE_002_20601_8_5_bwasplitter_out
PK_PAK_011_15045_1_11_bwasplitter_out
PK_PAK_011_15045_2_11_bwasplitter_out
PK_PAK_011_16552_1_11_bwasplitter_out





# merge bams
cd /lustre/scratch118/infgen/team133/sd21/hc/GENOME/POPULATION_DIVERSITY/MAPPING

echo "while read name lane; do ls -1 $"{name}"_*_out/$"{name}"*merged.sorted.marked.bam > $"{name}".bamlist ; \
samtools-1.3 merge -c -b $"{name}".bamlist $"{name}".merge.bam; done < ../samples_lanes.list" > run_merge

chmod a+x run_merge

bsub.py --queue yesterday 5 merge_bams ./run_merge


# cleanup
mkdir BAM_STATS
mv *out/*stat* BAM_STATS/
rm -r *out



# getting reads from guillumes global data
cd /nfs/users/nfs_g/gs12/link_lustre/ANGSD_223

cat /nfs/users/nfs_g/gs12/link_lustre/ANGSD_223/*.list > /lustre/scratch118/infgen/team133/sd21/hc/GENOME/POPULATION_DIVERSITY/GS_ORIGINAL/samples.bamlist

cd /lustre/scratch118/infgen/team133/sd21/hc/GENOME/POPULATION_DIVERSITY/GS_ORIGINAL/

while read list; do name=$( echo ${list} | cut -f 2 -d "." ); samtools view -ub --threads 4 ${list} | samtools sort -n - | samtools fastq -1 ${name}.R1.fast
q.gz -2 ${name}.R2.fastq.gz - ; done < samples.bamlist  &



# extracting name and lane IDs from Guillaumes BAM files
while read file; do
name=$( echo ${file} | cut -f 2 -d "." )
lane=$( samtools view -H $file | grep map_splitter | cut -f3 | awk -F " " '{print $NF}' | awk -F "/" '{print $NF}' | sed -e 's/_2.fastq.gz//g' )
echo -e "$name\t$lane" >> sample_lanes.list;
done < samples.bamlist


cd /lustre/scratch118/infgen/team133/sd21/hc/GENOME/POPULATION_DIVERSITY/GS_ORIGINAL/RAW
while read sample lane; do pathfind -t lane -i $lane --symlink ./ --rename --filetype fastq ; done < ../sample_lanes.list


cd ..

sort sample_lanes.list | uniq -c | awk '{print $2,$3}' OFS="\t" | sed -e "s/\#/_/g" > sample_lanes.list2

mkdir MAPPING
cd MAPPING

echo "while read name lane; do \
~sd21/bash_scripts/run_bwamem_splitter $"{name}"_$"{lane}" /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/POPULATION_DIVERSITY/HAEM_V4_final.chr.fa \
/lustre/scratch118/infgen/team133/sd21/hc/GENOME/POPULATION_DIVERSITY/GS_ORIGINAL/RAW/$"{lane}"_1.fastq.gz \
/lustre/scratch118/infgen/team133/sd21/hc/GENOME/POPULATION_DIVERSITY/GS_ORIGINAL/RAW/$"{lane}"_2.fastq.gz; \
done < /lustre/scratch118/infgen/team133/sd21/hc/GENOME/POPULATION_DIVERSITY/GS_ORIGINAL/sample_lanes.list2" > run_mapping

chmod a+x run_mapping
./run_mapping


# looks like some samples failed to map properly, which I think is due to hitting a memory limit. To check and print the suspect dirs, will check to see if
the flagstat file is made
# in this case, 870 samples were ok, 170 failed

for i in *out; do if [ ! -f ${i}/*flag* ]; then echo ${i} ; fi; done | wc -l

#--- once happy those samples are no good, remove them with a modification to the above script, and then restart mapping

for i in *out; do if [ ! -f ${i}/*flag* ]; then rm -r ${i} ; fi; done

# merge bams
cd /lustre/scratch118/infgen/team133/sd21/hc/GENOME/POPULATION_DIVERSITY/GS_ORIGINAL/MAPPING

echo "while read name lane; do ls -1 $"{name}"_*_out/$"{name}"*merged.sorted.marked.bam > $"{name}".bamlist ; \
samtools merge --threads 7 -c -b $"{name}".bamlist $"{name}".merge.bam; done < ../sample_lanes.list2" > run_merge

chmod a+x run_merge

bsub.py --queue yesterday 1 --threads 7 merge_bams ./run_merge

# cleanup
mkdir BAM_STATS
mv *out/*stat* BAM_STATS/
rm -r *out
```




```shell
cd /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/POPULATION_DIVERSITY/RAW

for i in 16693_1#1 16693_2#1 16720_1#1 16720_2#1; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 16693_1#2 16693_2#2 16720_1#2 16720_2#2; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 16693_1#3 16693_2#3 16720_1#3 16720_2#3; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 16693_1#4 16693_2#4 16720_1#4 16720_2#4; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 16693_1#5 16693_2#5 16720_1#5 16720_2#5; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 16693_1#6 16693_2#6 16720_1#6 16720_2#6; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 16693_1#7  16693_2#7  16720_1#7  16720_2#7; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 16693_1#8  16693_2#8  16720_1#8  16720_2#8; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 16693_1#9  16693_2#9  16720_1#9  16720_2#9; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 16693_1#10  16693_2#10  16720_1#10  16720_2#10; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 16693_1#11  16693_2#11  16720_1#11  16720_2#11; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 16693_1#12  16693_2#12  16720_1#12  16720_2#12; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 16693_1#13  16693_2#13  16720_1#13  16720_2#13; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 16693_1#14  16693_2#14  16720_1#14  16720_2#14; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 16693_1#15  16693_2#15  16720_1#15  16720_2#15; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 16693_1#16  16693_2#16  16720_1#16  16720_2#16; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 16693_1#17  16693_2#17  16720_1#17  16720_2#17; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 16693_1#18  16693_2#18  16720_1#18  16720_2#18; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 16693_1#19  16693_2#19  16720_1#19  16720_2#19; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 16693_1#20  16693_2#20  16720_1#20  16720_2#20; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 16693_1#21  16693_2#21  16720_1#21  16720_2#21; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 16693_1#22  16693_2#22  16720_1#22  16720_2#22; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 16693_1#23  16693_2#23  16720_1#23  16720_2#23; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 16693_1#24  16693_2#24  16720_1#24  16720_2#24; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 21766_7#1; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 21766_7#10; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 21766_7#11; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 21766_7#2; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 21766_7#3; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 21766_7#4; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 21766_7#5; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 21766_7#6; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 21766_7#7; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 21766_7#8; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 21766_7#9; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 20601_8#4; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 20601_8#5; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 20601_8#6; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 20601_8#7; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 20601_8#8; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 15150_1#34  15150_2#34  16553_1#34  16553_2#34; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 15150_1#35  15150_2#35  16553_1#35  16553_2#35; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 15150_1#36  15150_2#36  16553_1#36  16553_2#36; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 15150_1#19  15150_2#19  16553_1#19  16553_2#19; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 15150_1#20  15150_2#20  16553_1#20  16553_2#20; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 15150_1#21  15150_2#21  16553_1#21  16553_2#21; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 15150_1#22  15150_2#22  16553_1#22  16553_2#22; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 15150_1#23  15150_2#23  16553_1#23  16553_2#23; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 15150_1#24  15150_2#24  16553_1#24  16553_2#24; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 15150_1#25  15150_2#25  16553_1#25  16553_2#25; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 15150_1#26  15150_2#26  16553_1#26  16553_2#26; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 15150_1#27  15150_2#27  16553_1#27  16553_2#27; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 15150_1#28  15150_2#28  16553_1#28  16553_2#28; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 15150_1#29  15150_2#29  16553_1#29  16553_2#29; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 15150_1#30  15150_2#30  16553_1#30  16553_2#30; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 15150_1#31  15150_2#31  16553_1#31  16553_2#31; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 15150_1#32  15150_2#32  16553_1#32  16553_2#32; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 15150_1#33  15150_2#33  16553_1#33  16553_2#33; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 15045_1#1  15045_2#1  16552_1#1  16552_2#1; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 15045_1#2  15045_2#2  16552_1#2  16552_2#2; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 15045_1#3  15045_2#3  16552_1#3  16552_2#3; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 15045_1#4  15045_2#4  16552_1#4  16552_2#4; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 15045_1#5  15045_2#5  16552_1#5  16552_2#5; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 15045_1#6  15045_2#6  16552_1#6  16552_2#6; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 15045_1#7  15045_2#7  16552_1#7  16552_2#7; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 15045_1#8  15045_2#8  16552_1#8  16552_2#8; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 15045_1#9  15045_2#9  16552_1#9  16552_2#9; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 15045_1#10  15045_2#10  16552_1#10  16552_2#10; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 15045_1#11  15045_2#11  16552_1#11  16552_2#11; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 15045_1#12  15045_2#12  16552_1#12  16552_2#12; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 15045_1#13  15045_2#13  16552_1#13  16552_2#13; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 15045_1#14  15045_2#14  16552_1#14  16552_2#14; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done
for i in 15045_1#15  15045_2#15  16552_1#15  16552_2#15; do pathfind -t lane -i $i --symlink ./ --rename --filetype fastq ; done



# mapping
#--- run_bwamem_splitter automatically submits jobs to LSF, so no need to bsub it. Will slowly work though the list when running screen
#--- samples_lanes.list was curated in excel (See: "Genome_paper_population_diversity_samples.xls"), but is a tab delimited text file containing sample name and sequencing lane ID for all samples to be mapped
# CH_SWI_001	16693_1#13
# CH_SWI_001	16693_2#13
# CH_SWI_001	16720_1#13
# CH_SWI_001	16720_2#13
# CH_SWI_002	16693_1#14

cd /lustre/scratch118/infgen/team133/sd21/hc/GENOME/POPULATION_DIVERSITY/MAPPING

echo "while read name lane; do \
~sd21/bash_scripts/run_bwamem_splitter $"{name}"_$"{lane}" /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/POPULATION_DIVERSITY/HAEM_V4_final.chr.fa \
/nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/POPULATION_DIVERSITY/RAW/$"{lane}"_1.fastq.gz \
/nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/POPULATION_DIVERSITY/RAW/$"{lane}"_2.fastq.gz; \
done < /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/POPULATION_DIVERSITY/samples_lanes.list" > run_mapping
chmod a+x run_mapping
screen
./run_mapping &



GB_IRE_002_20601_8_5_bwasplitter_out
PK_PAK_011_15045_1_11_bwasplitter_out
PK_PAK_011_15045_2_11_bwasplitter_out
PK_PAK_011_16552_1_11_bwasplitter_out





# merge bams
cd /lustre/scratch118/infgen/team133/sd21/hc/GENOME/POPULATION_DIVERSITY/MAPPING

echo "while read name lane; do ls -1 $"{name}"_*_out/$"{name}"*merged.sorted.marked.bam > $"{name}".bamlist ; \
samtools-1.3 merge -c -b $"{name}".bamlist $"{name}".merge.bam; done < ../samples_lanes.list" > run_merge

chmod a+x run_merge

bsub.py --queue yesterday 5 merge_bams ./run_merge

# index bams
echo -e "for i in *bam; do samtools index -@ 7 -b $"{i}"; done" > run_index_bams
chmod a+x run_index_bams
bsub.py --queue yesterday 1 --threads 7 index_bams ./run_index_bams


# cleanup
mkdir BAM_STATS
mv *out/*stat* BAM_STATS/
rm -r *out



# getting reads from guillumes global data
cd /nfs/users/nfs_g/gs12/link_lustre/ANGSD_223

cat /nfs/users/nfs_g/gs12/link_lustre/ANGSD_223/*.list > /lustre/scratch118/infgen/team133/sd21/hc/GENOME/POPULATION_DIVERSITY/GS_ORIGINAL/samples.bamlist

cd /lustre/scratch118/infgen/team133/sd21/hc/GENOME/POPULATION_DIVERSITY/GS_ORIGINAL/

while read list; do name=$( echo ${list} | cut -f 2 -d "." ); samtools view -ub --threads 4 ${list} | samtools sort -n - | samtools fastq -1 ${name}.R1.fastq.gz -2 ${name}.R2.fastq.gz - ; done < samples.bamlist  &



# extracting name and lane IDs from Guillaumes BAM files
while read file; do
name=$( echo ${file} | cut -f 2 -d "." )
lane=$( samtools view -H $file | grep map_splitter | cut -f3 | awk -F " " '{print $NF}' | awk -F "/" '{print $NF}' | sed -e 's/_2.fastq.gz//g' )
echo -e "$name\t$lane" >> sample_lanes.list;
done < samples.bamlist


cd /lustre/scratch118/infgen/team133/sd21/hc/GENOME/POPULATION_DIVERSITY/GS_ORIGINAL/RAW
while read sample lane; do pathfind -t lane -i $lane --symlink ./ --rename --filetype fastq ; done < ../sample_lanes.list


cd ..

sort sample_lanes.list | uniq -c | awk '{print $2,$3}' OFS="\t" | sed -e "s/\#/_/g" > sample_lanes.list2

mkdir MAPPING
cd MAPPING

echo "while read name lane; do \
~sd21/bash_scripts/run_bwamem_splitter $"{name}"_$"{lane}" /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/POPULATION_DIVERSITY/HAEM_V4_final.chr.fa \
/lustre/scratch118/infgen/team133/sd21/hc/GENOME/POPULATION_DIVERSITY/GS_ORIGINAL/RAW/$"{lane}"_1.fastq.gz \
/lustre/scratch118/infgen/team133/sd21/hc/GENOME/POPULATION_DIVERSITY/GS_ORIGINAL/RAW/$"{lane}"_2.fastq.gz; \
done < /lustre/scratch118/infgen/team133/sd21/hc/GENOME/POPULATION_DIVERSITY/GS_ORIGINAL/sample_lanes.list2" > run_mapping

chmod a+x run_mapping
./run_mapping


# looks like some samples failed to map properly, which I think is due to hitting a memory limit. To check and print the suspect dirs, will check to see if the flagstat file is made
# in this case, 870 samples were ok, 170 failed

for i in *out; do if [ ! -f ${i}/*flag* ]; then echo ${i} ; fi; done | wc -l

#--- once happy those samples are no good, remove them with a modification to the above script, and then restart mapping

for i in *out; do if [ ! -f ${i}/*flag* ]; then rm -r ${i} ; fi; done

# merge bams
cd /lustre/scratch118/infgen/team133/sd21/hc/GENOME/POPULATION_DIVERSITY/GS_ORIGINAL/MAPPING

echo "while read name lane; do ls -1 $"{name}"_*_out/$"{name}"*merged.sorted.marked.bam > $"{name}".bamlist ; \
samtools merge --threads 7 -c -b $"{name}".bamlist $"{name}".merge.bam; done < ../sample_lanes.list2" > run_merge

chmod a+x run_merge

bsub.py --queue yesterday 1 --threads 7 merge_bams ./run_merge

# index bams
echo -e "for i in *bam; do samtools index -@ 7 -b $"{i}"; done" > run_index_bams
chmod a+x run_index_bams
bsub.py --queue yesterday 1 --threads 7 index_bams ./run_index_bams

# cleanup
mkdir BAM_STATS
mv *out/*stat* BAM_STATS/
rm -r *out
rm *bamlist




echo -e "for i in *.bam; do gatk-4.0.3.0 AddOrReplaceReadGroups --INPUT=$"{i}" --OUTPUT=$"{i%.bam}d.bam" --RGID=$"{i%.merge.bam}" --RGLB=$"{i%.merge.bam}" --RGPL=HS --RGPU=WSI --RGSM=$"{i%.merge.bam}"; done" > run_fix_bam_header
chmod a+x run_fix_bam_header
bsub.py --queue yesterday 1 --threads 7 fix_bam_header ./run_fix_bam_header

gatk-4.0.3.0 HaplotypeCaller --input ZAI_ZAI_OA_015.merge2.bam --output ZAI_ZAI_OA_015.merge.vcf --reference ../../HAEM_V4_final.chr.fa


# create bam list using full path to bams - this allos bams to be anywhere
#ls $PWD/*bam > bam.list   
BAM_LIST=/lustre/scratch118/infgen/team133/sd21/hc/GENOME/POPULATION_DIVERSITY/MAPPING/MERGED_BAMS/bam.list   # new samples

# specify data using full paths to bam list and fasta refernce sequence
BAM_LIST=/lustre/scratch118/infgen/team133/sd21/hc/GENOME/POPULATION_DIVERSITY/GS_ORIGINAL/MAPPING/bam.list   # guillaumes samples
REFERENCE=/lustre/scratch118/infgen/team133/sd21/hc/GENOME/POPULATION_DIVERSITY/HAEM_V4_final.chr.fa

mkdir GATK_HC_GVCF
cd GATK_HC_GVCF
mkdir REF_sequences

# make reference seqeunces
cp ${REFERENCE} REF_sequences/REF.fa
samtools faidx REF_sequences/REF.fa
samtools dict REF_sequences/REF.fa > REF_sequences/REF.dict
cp ${BAM_LIST} bam.list
grep ">" REF_sequences/REF.fa | sed -e 's/>//g' > sequences.list

ulimit -c unlimited
#/software/pathogen/external/apps/usr/local/gatk-4.0.3.0/gatk-package-4.0.3.0-local.jar

# make jobs
while read BAM; do \
	n=1
	SAMPLE=$( echo ${BAM} | awk -F '/' '{print $NF}' | sed -e 's/.bam//g' )  
	mkdir ${SAMPLE}_GATK_HC_GVCF
	mkdir ${SAMPLE}_GATK_HC_GVCF/LOGFILES
	echo "gatk-4.0.3.0 GatherVcfsCloud \\" > ${SAMPLE}_GATK_HC_GVCF/run_gather_${SAMPLE}_gvcf
	while read SEQUENCE; do
	echo -e "java -Xmx15G -jar /software/pathogen/external/apps/usr/local/gatk-4.0.3.0/gatk-package-4.0.3.0-local.jar HaplotypeCaller --input ${BAM} --output ${SAMPLE}_GATK_HC_GVCF/${n}.${SAMPLE}.${SEQUENCE}.tmp.gvcf.gz --reference $PWD/REF_sequences/REF.fa --intervals ${SEQUENCE} --emit-ref-confidence GVCF " > ${SAMPLE}_GATK_HC_GVCF/run_hc_${SAMPLE}.${SEQUENCE}.tmp.job_${n};
	echo -e "--input ${PWD}/${SAMPLE}_GATK_HC_GVCF/${n}.${SAMPLE}.${SEQUENCE}.tmp.gvcf.gz \\" >> ${SAMPLE}_GATK_HC_GVCF/run_gather_${SAMPLE}_gvcf;
	let "n+=1"; done < sequences.list;

	echo -e "--output ${PWD}/${SAMPLE}_GATK_HC_GVCF/${SAMPLE}.gvcf.gz; tabix -p vcf ${PWD}/${SAMPLE}_GATK_HC_GVCF/${SAMPLE}.gvcf.gz" >> ${SAMPLE}_GATK_HC_GVCF/run_gather_${SAMPLE}_gvcf;

	echo -e "rm ${PWD}/${SAMPLE}_GATK_HC_GVCF/*.tmp.* && mv ${PWD}/${SAMPLE}_GATK_HC_GVCF/*.[oe] ${SAMPLE}_GATK_HC_GVCF/LOGFILES && cd ${PWD} && mv ${PWD}/${SAMPLE}_GATK_HC_GVCF ${PWD}/${SAMPLE}_GATK_HC_GVCF_complete" > ${SAMPLE}_GATK_HC_GVCF/run_clean_${SAMPLE};

	chmod a+x ${SAMPLE}_GATK_HC_GVCF/run_*

	# setup job conditions
	JOBS=$( ls -1 ${SAMPLE}_GATK_HC_GVCF/run_hc_* | wc -l )
	ID="U$(date +%s)"

	#submit job array to call variants put scaffold / contig
	bsub -q long -R'span[hosts=1] select[mem>15000] rusage[mem=15000]' -n 6 -M15000 -J GATK_HC_${ID}_[1-$JOBS]%100 -e ${SAMPLE}_GATK_HC_GVCF/GATK_HC_${ID}_[1-$JOBS].e -o ${SAMPLE}_GATK_HC_GVCF/GATK_HC_${ID}_[1-$JOBS].o "./${SAMPLE}_GATK_HC_GVCF/run_hc_${SAMPLE}.*job_\$LSB_JOBINDEX"

	#submit job to gather gvcfs into a single, per sample gvcf
	bsub -q normal -w "done(GATK_HC_${ID}_[1-$JOBS])" -R'span[hosts=1] select[mem>500] rusage[mem=500]' -n 1 -M500 -J GATK_HC_${ID}_gather_gvcfs -e ${SAMPLE}_GATK_HC_GVCF/GATK_HC_${ID}_gather_gvcfs.e -o ${SAMPLE}_GATK_HC_GVCF/GATK_HC_${ID}_gather_gvcfs.o "./${SAMPLE}_GATK_HC_GVCF/run_gather_${SAMPLE}_gvcf"

	# clean up
	bsub -q normal -w "done(GATK_HC_${ID}_gather_gvcfs)" -R'span[hosts=1] select[mem>500] rusage[mem=500]' -n 1 -M500 -J GATK_HC_${ID}_clean -e ${SAMPLE}_GATK_HC_GVCF/GATK_HC_${ID}_clean.e -o ${SAMPLE}_GATK_HC_GVCF/GATK_HC_${ID}_clean.o "./${SAMPLE}_GATK_HC_GVCF/run_clean_${SAMPLE}"

	sleep 1
done < ${BAM_LIST}







#echo "gatk-4.0.3.0 GatherVcfsCloud \\" > run_gather_${SAMPLE}_gvcf
#echo -e "--input ${SAMPLE}_GATK_HC_GVCF/${n}.${SAMPLE}.${SEQUENCE}.g.vcf.tmp \\" >> run_gather_${SAMPLE}_gvcf
#
# echo -e "--output ${SAMPLE}.g.vcf" >> run_gather_${SAMPLE}_gvcf
# chmod a+x run_gather_${SAMPLE}_gvcf

mkdir GATK_HC_MERGED
cd GATK_HC_MERGED

ls -1 ../GATK_HC/*complete/*gz > gvcf.list
ls -1 ../../GS_ORIGINAL/GATK_HC/*complete/*gz >> gvcf.list

GVCF_LIST=$PWD/gvcf.list
REFERENCE=/lustre/scratch118/infgen/team133/sd21/hc/GENOME/POPULATION_DIVERSITY/HAEM_V4_final.chr.fa

# echo -e "java -Xmx15G -jar /software/pathogen/external/apps/usr/local/gatk-4.0.3.0/gatk-package-4.0.3.0-local.jar CombineGVCFs -R ${REFERENCE} \\" > run_merge_gvcfs
# while read SAMPLE; do
# echo -e "--variant ${SAMPLE} \\" >> run_merge_gvcfs;
#    done < ${GVCF_LIST}
#    echo -e "--output cohort.g.vcf.gz" >> run_merge_gvcfs

# chmod a+x run_merge_gvcfs
# bsub.py --queue hugemem --threads 30 200 merge_vcfs "./run_merge_gvcfs"
# threads make a big difference, even thoguh they are not a parameter in the tool


grep ">" ${REFERENCE} | sed -e 's/>//g' > sequences.list

n=1
while read SEQUENCE; do
echo -e "java -Xmx15G -jar /software/pathogen/external/apps/usr/local/gatk-4.0.3.0/gatk-package-4.0.3.0-local.jar CombineGVCFs -R ${REFERENCE} --intervals ${SEQUENCE} \\" > ${n}.run_merge_gvcfs_${SEQUENCE}
while read SAMPLE; do
echo -e "--variant ${SAMPLE} \\" >> ${n}.run_merge_gvcfs_${SEQUENCE};
   done < ${GVCF_LIST}
   echo -e "--output ${SEQUENCE}.cohort.g.vcf.gz" >> ${n}.run_merge_gvcfs_${SEQUENCE};
   let "n+=1"; done < sequences.list

chmod a+x *run_merge_gvcfs*

for i in *run_merge_gvcfs*; do
bsub.py --queue hugemem --threads 30 200 merge_vcfs "./${i}"; done


# single genotyping call,
# bsub.py --queue hugemem --threads 30 200 genotype_vcfs "java -Xmx15G -jar /software/pathogen/external/apps/usr/local/gatk-4.0.3.0/gatk-package-4.0.3.0-local.jar GenotypeGVCFs \
#    -R /lustre/scratch118/infgen/team133/sd21/hc/GENOME/POPULATION_DIVERSITY/HAEM_V4_final.chr.fa \
#    -V cohort.g.vcf.gz \
#    -O cohort.vcf.gz"


# split each chromosome up into separate jobs, and run genotyping on each individually.   
n=1
while read SEQUENCE; do
echo -e "java -Xmx15G -jar /software/pathogen/external/apps/usr/local/gatk-4.0.3.0/gatk-package-4.0.3.0-local.jar GenotypeGVCFs \
-R /lustre/scratch118/infgen/team133/sd21/hc/GENOME/POPULATION_DIVERSITY/HAEM_V4_final.chr.fa \
-V ${SEQUENCE}.cohort.g.vcf.gz \
--intervals ${SEQUENCE} \
-O ${n}.${SEQUENCE}.cohort.vcf.gz" > run_hc_genotype.${SEQUENCE}.tmp.job_${n};
let "n+=1"; done < sequences.list

chmod a+x run_hc_genotype*

mkdir LOGFILES

	# setup job conditions
	JOBS=$( ls -1 run_hc_* | wc -l )
	ID="U$(date +%s)"

# ln -s /lustre/scratch118/infgen/team133/sd21/hc/GENOME/POPULATION_DIVERSITY/MAPPING/cohort.g.vcf.gz
# ln -s /lustre/scratch118/infgen/team133/sd21/hc/GENOME/POPULATION_DIVERSITY/MAPPING/cohort.g.vcf.gz.tbi

bsub -q yesterday -R'span[hosts=1] select[mem>20000] rusage[mem=20000]' -n 6 -M20000 -J GATK_HC_GENOTYPE_${ID}_[1-$JOBS] -e LOGFILES/GATK_HC_GENOTYPE_${ID}_[1-$JOBS].e -o LOGFILES/GATK_HC_GENOTYPE_${ID}_[1-$JOBS].o "./run_hc_*\$LSB_JOBINDEX"

```

## MtDNA analysis
```shell

###########################################################################################
# MTDNA ANALYSIS


cd /lustre/scratch118/infgen/team133/sd21/hc/GENOME/POPULATION_DIVERSITY
cat samples_lanes.list GS_ORIGINAL/sample_lanes.list2 | cut -f1 | sort | uniq > MAPPING/GATK_HC_MERGED/MTDNA_ANALYSIS/samples.list

cd /lustre/scratch118/infgen/team133/sd21/hc/GENOME/POPULATION_DIVERSITY/MAPPING/GATK_HC_MERGED/MTDNA_ANALYSIS

samtools faidx HAEM_V4_final.chr.fa hcontortus_chr_mtDNA_arrow_pilon > hcontortus_chr_mtDNA_arrow_pilon.fa
samtools faidx hcontortus_chr_mtDNA_arrow_pilon.fa
samtools dict hcontortus_chr_mtDNA_arrow_pilon.fa > hcontortus_chr_mtDNA_arrow_pilon.dict

mkdir INDV_VCFS


# rm run_select_variants
# while read NAME; do \
# echo -e "java -Xmx15G -jar /software/pathogen/external/apps/usr/local/gatk-4.0.3.0/gatk-package-4.0.3.0-local.jar SelectVariants \
# --variant 7.hcontortus_chr_mtDNA_arrow_pilon.cohort.vcf.gz \
# --reference hcontortus_chr_mtDNA_arrow_pilon.fa \
# --sample-name ${NAME} \
# -select 'vc.getGenotype(\"${NAME}\").isHomVar()' \
# -select-type SNP \
# -select-type INDEL \
# -xl-select-type MNP \
# --output INDV_VCFS/${NAME}.raw.vcf.gz" >> run_select_variants ; done < samples.list

rm run_select_variants
while read NAME; do \
echo -e "java -Xmx15G -jar /software/pathogen/external/apps/usr/local/gatk-4.0.3.0/gatk-package-4.0.3.0-local.jar SelectVariants \
--variant 7.hcontortus_chr_mtDNA_arrow_pilon.cohort.vcf.gz \
--reference hcontortus_chr_mtDNA_arrow_pilon.fa \
--sample-name ${NAME} \
-select-type SNP \
-select-type INDEL \
-xl-select-type MNP \
-select 'vc.getGenotype(\"${NAME}\").isHomVar()' \
--output INDV_VCFS/${NAME}.raw.vcf.gz" >> run_select_variants ; done < samples.list


chmod a+x run_select_variants

bsub.py --queue yesterday 1 selectvars2 ./run_select_variants

mkdir INDV_FASTAS

for i in ` ls -1 *.raw.vcf.gz | sed -e 's/.raw.vcf.gz//g' ` ; do bsub -I java -jar ~sd21/lustre118_link/software/SNP_CALLING/GATK-3.6/GenomeAnalysisTK.jar -T FastaAlternateReferenceMaker --variant ${i}.raw.vcf.gz -o ${i}.mtDNA_consensus.fa -R ../hcontortus_chr_mtDNA_arrow_pilon.fa; sed -i 's/^>.*.$/>'"${i}"'/g' ${i}.mtDNA_consensus.fa; done &


java -jar GenomeAnalysisTK.jar \
-T VariantFiltration \
-R hs37d5.fa \
-V unfiltered_snps.vcf \
--filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
--filterName "my_snp_filter" \
-o filtered_snps.vc






export PATH="/nfs/users/nfs_s/sd21/lustre118_link/software/anaconda2/bin:$PATH"





# SNP filtering using BCFTOOLS,
# generates a mask file of positions with low depth (<10) and MQ (<30) score, whcih it will convert to Ns when the consensus is made
# only considers SNPs that are homozygous variant, MQ > 30, and population AF > 0.01

bcftools-1.9 query --list-samples 7.hcontortus_chr_mtDNA_arrow_pilon.cohort.vcf.gz > samples.list

while read SAMPLE; do
bcftools-1.9 view --samples $SAMPLE 7.hcontortus_chr_mtDNA_arrow_pilon.cohort.vcf.gz | bcftools-1.9 view -i 'FORMAT/DP[0]<10 | MQ[0]<30' > $SAMPLE.mask
bcftools-1.9 view --samples $SAMPLE 7.hcontortus_chr_mtDNA_arrow_pilon.cohort.vcf.gz | bcftools-1.9 view -i 'TYPE="snp" & FORMAT/GT[0]="1/1" & MQ[0]>30 & AF[0]>0.01' > $SAMPLE.keep
bcftools-1.9 consensus  --sample $SAMPLE -i 'TYPE="snp" & FORMAT/GT[0]="1/1" & MQ[0]>30 & AF[0]>0.01' -f hcontortus_chr_mtDNA_arrow_pilon.fa --mask $SAMPLE.mask 7.hcontortus_chr_mtDNA_arrow_pilon.cohort.vcf.gz > $SAMPLE.pseudomtDNA.fa;
sed -e "s/hcontortus_chr_mtDNA_arrow_pilon/$SAMPLE/g" $SAMPLE.pseudomtDNA.fa | sed -e 's/_OA_/_/g' -e 's/_CH_/_/g' -e 's/_//g' > $SAMPLE.pseudomtDNA.fa.tmp; mv $SAMPLE.pseudomtDNA.fa.tmp $SAMPLE.pseudomtDNA.fa;
done < samples.list





cat *pseudomtDNA.fa hplacei_mtDNA_genome_AP017687.1.fa tcircumcincta_mtDNA_genome_AP017699.1.fa > all_pseudoref.outgroups.fa




```




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

metadata<-read.table("sample_metadata_colours.list",header=T,comment.char="",sep="\t")

rubi.VCF <- read.vcfR("allsamples.mtDNA.filtered.vcf.gz")
pop.data <- read.table("samples.pops.list", sep = "\t", header = F)
gl.rubi <- vcfR2genlight(rubi.VCF)
ploidy(gl.rubi) <- 1

pop(gl.rubi) <- metadata$country



# distance matrix from genlight object
#x.dist <- poppr::bitwise.dist(gl.rubi)



# make a tree
#tree <- aboot(gl.rubi, tree = "upgma", distance = bitwise.dist, sample = 100, showtree = F, cutoff = 50, quiet = T)
#write.tree(tree, file="MyNewickTreefile.nwk")


#cols <- brewer.pal(n = nPop(gl.rubi), name = "Dark2")
#plot.phylo(tree, cex = 0.3, font = 2, adj = 0)
#nodelabels(tree$node.label, adj = c(1.3, -0.5), frame = "n", cex = 0.3,font = 3, xpd = TRUE)
#legend(35,10,c("CA","OR","WA"),cols, border = FALSE, bty = "n")
#legend('topleft', legend = c("CA","OR","WA"),fill = cols, border = FALSE, bty = "n", cex = 2)
#axis(side = 1)
#title(xlab = "Genetic distance (proportion of loci that are different)")





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


#p12 <- ggplot(rubi.pca.scores, aes(x=PC1, y=PC2, colour=pop, label=pop)) + geom_point(size=2)+ theme_bw() + geom_text_repel(data = subset(rubi.pca.scores, pop == "ZAI" ))
#p34 <- ggplot(rubi.pca.scores, aes(x=PC3, y=PC4, colour=pop, label=pop)) + geom_point(size=2)+ theme_bw() + geom_text_repel(data = subset(rubi.pca.scores, pop == "ZAI" ))
#p56 <- ggplot(rubi.pca.scores, aes(x=PC5, y=PC6, colour=pop, label=pop)) + geom_point(size=2)+ theme_bw() + geom_text_repel(data = subset(rubi.pca.scores, pop == "ZAI" ))
#p12 + p34 + p56


country_colours_shapes1 <- c("#31A197","#E15956","#EF724B","#D35E5C","#606EB8","#6570B0","#34AFE7","#6973A8","#3C9C93","#6E75A0","#C56462","#B76968","#A64EB4","#727898","#A96F6E","#9B7474","#3FA8D8","#8D7A7A")
country_colours_shapes2 <-c("16","16","16","16","17","16","16","17","16","16","16","16","17","16","16","16","17","16")
country_colours_shapes3 <-c("0.5","0.5","0.5","0.5","1","0.5","0.5","1","0.5","0.5","0.5","0.5","1","0.5","0.5","0.5","1","0.5")

country_colours_shapes <- rbind(country_colours_shapes1,country_colours_shapes2,country_colours_shapes3)

#add names to the data
colnames(country_colours_shapes) <- c("Australia","Benin","Brazil","Cape_Verde","Switzerland","France","Guadeloupe","United_Kingdom","Indonesia","Italy","Morocco","Namibia","Pakistan","Portugal","South_Africa","São_Tomé","USA","DRC")

# sort the columns by name - this is really important.
country_colours_shapes <- country_colours_shapes[ , order(names(as.data.frame(country_colours_shapes)))]

PC1_variance <- formatC(head(rubi.pca$eig)[1]/sum(rubi.pca$eig)*100)
PC2_variance <- formatC(head(rubi.pca$eig)[2]/sum(rubi.pca$eig)*100)
#--- note: formatC() limits output to two decimal places

ggplot(rubi.pca.scores, aes(x=PC1, y=PC2, colour=pop, label=pop, shape=pop)) +
          geom_point(data = subset(rubi.pca.scores, pop == "Australia" | pop == "Benin" | pop == "Brazil" | pop == "Cape_Verde" | pop == "France" | pop == "Guadeloupe" | pop == "Indonesia" | pop == "Italy" | pop == "Morocco" | pop == "Namibia" | pop == "Portugal" | pop == "South_Africa" | pop == "São_Tomé"| pop == "DRC"), size=3,alpha=0.8)+
          geom_point(data = subset(rubi.pca.scores, pop == "Switzerland" | pop == "USA" | pop== "Pakistan" | pop == "United_Kingdom"), size=3,alpha=0.8)+
          theme_bw()+
          scale_colour_manual(values = country_colours_shapes[1,])+
          scale_shape_manual(values = as.numeric(country_colours_shapes[2,]))+
          labs(x=paste("PC1: ",PC1_variance,"% variance"),y=paste("PC2: ",PC2_variance,"% variance"))+
          geom_text_repel(data = subset(rubi.pca.scores[grep("^GB_ISEN1_001", row.names(rubi.pca.scores)),]),label="ISE.N1",show.legend = FALSE,point.padding=1)

ggsave("global_diversity_mtDNA_SNPs.pdf",height=6,width=7.5,useDingbats = FALSE)
ggsave("global_diversity_mtDNA_SNPs.png",height=6,width=7.5)

```

```shell
scp sd21@pcs5.internal.sanger.ac.uk:/nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/POPULATION_DIVERSITY/VARIANTS/MTNDA/global_diversity_mtDNA_SNPs.* ~/Documents/workbook/hcontortus_genome/04_analysis
```

![Global diversity - mtDNA](04_analysis/global_diversity_mtDNA_SNPs.png)
Fig - Global genetic diversity based on mtDNA variation




```R
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
vcftools-0.1.14 --gzvcf ../1.hcontortus_chrX_Celeg_TT_arrow_pilon.cohort.vcf.gz --fst-window-size 100000 --weir-fst-pop ZAI.population.list --weir-fst-pop AUS.population.list --weir-fst-pop GB.population.list --weir-fst-pop CAP.population.list --weir-fst-pop FRG.population.list --weir-fst-pop ITA.population.list --weir-fst-pop PK.population.list --weir-fst-pop US.population.list --weir-fst-pop BEN.population.list --weir-fst-pop CH.population.list --weir-fst-pop MOR.population.list --weir-fst-pop POR.population.list --weir-fst-pop STO.population.list --weir-fst-pop STA.population.list --weir-fst-pop BRA.population.list --weir-fst-pop FRA.population.list --weir-fst-pop IND.population.list --weir-fst-pop NAM.population.list --out chrX_fst_100k_allpop

vcftools-0.1.14 --gzvcf ../2.hcontortus_chr2_Celeg_TT_arrow_pilon.cohort.vcf.gz --fst-window-size 100000 --weir-fst-pop ZAI.population.list --weir-fst-pop AUS.population.list --weir-fst-pop GB.population.list --weir-fst-pop CAP.population.list --weir-fst-pop FRG.population.list --weir-fst-pop ITA.population.list --weir-fst-pop PK.population.list --weir-fst-pop US.population.list --weir-fst-pop BEN.population.list --weir-fst-pop CH.population.list --weir-fst-pop MOR.population.list --weir-fst-pop POR.population.list --weir-fst-pop STO.population.list --weir-fst-pop STA.population.list --weir-fst-pop BRA.population.list --weir-fst-pop FRA.population.list --weir-fst-pop IND.population.list --weir-fst-pop NAM.population.list --out chr2_fst_100k_allpop

vcftools-0.1.14 --gzvcf ../3.hcontortus_chr3_Celeg_TT_arrow_pilon.cohort.vcf.gz --fst-window-size 100000 --weir-fst-pop ZAI.population.list --weir-fst-pop AUS.population.list --weir-fst-pop GB.population.list --weir-fst-pop CAP.population.list --weir-fst-pop FRG.population.list --weir-fst-pop ITA.population.list --weir-fst-pop PK.population.list --weir-fst-pop US.population.list --weir-fst-pop BEN.population.list --weir-fst-pop CH.population.list --weir-fst-pop MOR.population.list --weir-fst-pop POR.population.list --weir-fst-pop STO.population.list --weir-fst-pop STA.population.list --weir-fst-pop BRA.population.list --weir-fst-pop FRA.population.list --weir-fst-pop IND.population.list --weir-fst-pop NAM.population.list --out chr3_fst_100k_allpop

vcftools-0.1.14 --gzvcf ../4.hcontortus_chr4_Celeg_TT_arrow_pilon.cohort.vcf.gz --fst-window-size 100000 --weir-fst-pop ZAI.population.list --weir-fst-pop AUS.population.list --weir-fst-pop GB.population.list --weir-fst-pop CAP.population.list --weir-fst-pop FRG.population.list --weir-fst-pop ITA.population.list --weir-fst-pop PK.population.list --weir-fst-pop US.population.list --weir-fst-pop BEN.population.list --weir-fst-pop CH.population.list --weir-fst-pop MOR.population.list --weir-fst-pop POR.population.list --weir-fst-pop STO.population.list --weir-fst-pop STA.population.list --weir-fst-pop BRA.population.list --weir-fst-pop FRA.population.list --weir-fst-pop IND.population.list --weir-fst-pop NAM.population.list --out chr4_fst_100k_allpop

vcftools-0.1.14 --gzvcf ../5.hcontortus_chr5_Celeg_TT_arrow_pilon.cohort.vcf.gz --fst-window-size 100000 --weir-fst-pop ZAI.population.list --weir-fst-pop AUS.population.list --weir-fst-pop GB.population.list --weir-fst-pop CAP.population.list --weir-fst-pop FRG.population.list --weir-fst-pop ITA.population.list --weir-fst-pop PK.population.list --weir-fst-pop US.population.list --weir-fst-pop BEN.population.list --weir-fst-pop CH.population.list --weir-fst-pop MOR.population.list --weir-fst-pop POR.population.list --weir-fst-pop STO.population.list --weir-fst-pop STA.population.list --weir-fst-pop BRA.population.list --weir-fst-pop FRA.population.list --weir-fst-pop IND.population.list --weir-fst-pop NAM.population.list --out chr5_fst_100k_allpop

vcftools-0.1.14 --gzvcf ../6.hcontortus_chrX_Celeg_TT_arrow_pilon.cohort.vcf.gz --fst-window-size 100000 --weir-fst-pop ZAI.population.list --weir-fst-pop AUS.population.list --weir-fst-pop GB.population.list --weir-fst-pop CAP.population.list --weir-fst-pop FRG.population.list --weir-fst-pop ITA.population.list --weir-fst-pop PK.population.list --weir-fst-pop US.population.list --weir-fst-pop BEN.population.list --weir-fst-pop CH.population.list --weir-fst-pop MOR.population.list --weir-fst-pop POR.population.list --weir-fst-pop STO.population.list --weir-fst-pop STA.population.list --weir-fst-pop BRA.population.list --weir-fst-pop FRA.population.list --weir-fst-pop IND.population.list --weir-fst-pop NAM.population.list --out chrX_fst_100k_allpop

vcftools-0.1.14 --gzvcf ../7.hcontortus_chr_mtDNA_arrow_pilon.cohort.vcf.gz --fst-window-size 1000 --weir-fst-pop ZAI.population.list --weir-fst-pop AUS.population.list --weir-fst-pop GB.population.list --weir-fst-pop CAP.population.list --weir-fst-pop FRG.population.list --weir-fst-pop ITA.population.list --weir-fst-pop PK.population.list --weir-fst-pop US.population.list --weir-fst-pop BEN.population.list --weir-fst-pop CH.population.list --weir-fst-pop MOR.population.list --weir-fst-pop POR.population.list --weir-fst-pop STO.population.list --weir-fst-pop STA.population.list --weir-fst-pop BRA.population.list --weir-fst-pop FRA.population.list --weir-fst-pop IND.population.list --weir-fst-pop NAM.population.list --out chrMT_fst_1k_allpop


# determine the position of singleton SNPs
vcftools-0.1.14 --gzvcf ../1.hcontortus_chrX_Celeg_TT_arrow_pilon.cohort.vcf.gz --singletons --out singletons_chrX
vcftools-0.1.14 --gzvcf ../2.hcontortus_chr2_Celeg_TT_arrow_pilon.cohort.vcf.gz --singletons --out singletons_chr2
vcftools-0.1.14 --gzvcf ../3.hcontortus_chr3_Celeg_TT_arrow_pilon.cohort.vcf.gz --singletons --out singletons_chr3
vcftools-0.1.14 --gzvcf ../4.hcontortus_chr4_Celeg_TT_arrow_pilon.cohort.vcf.gz --singletons --out singletons_chr4
vcftools-0.1.14 --gzvcf ../5.hcontortus_chr5_Celeg_TT_arrow_pilon.cohort.vcf.gz --singletons --out singletons_chr5
vcftools-0.1.14 --gzvcf ../6.hcontortus_chrX_Celeg_TT_arrow_pilon.cohort.vcf.gz --singletons --out singletons_chrX

# make a position file out of the singletons
for i in singletons_*; do cut -f1,2 ${i} > ${i}.2; done


# calculate window Pi 100kb of singletons only
echo -e "
vcftools-0.1.14 --gzvcf ../1.hcontortus_chrX_Celeg_TT_arrow_pilon.cohort.vcf.gz --positions singletons_chrX.singletons.2 --window-pi 100000 --out chrX.singletons.pi

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
vcftools-0.1.14 --gzvcf ../1.hcontortus_chrX_Celeg_TT_arrow_pilon.cohort.vcf.gz --exclude-positions singletons_chrX.singletons.2 --window-pi 100000 --out chrX.shared.pi

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


# check gene density per 100 k window
ln -s ../../../REF/HAEM_V4_final.chr.fa
samtools faidx HAEM_V4_final.chr.fa
cut -f1,2 HAEM_V4_final.chr.fa.fai > HAEM_V4_final.chr.genome

bedtools-2 makewindows -g  HAEM_V4_final.chr.genome -w 100000 > HAEM_V4_final.chr.500k.bed
awk '$3=="gene" {print $1,$4,$5,$9}' OFS="\t" ANNOTATION.gff >HCON_V4.gene.bed

bedtools-2 coverage -a HAEM_V4_final.chr.500k.bed -b HCON_V4.gene.bed -counts > HCON_V4.gene.counts

R-3.5.0
library(ggplot2)
a<-read.table("HCON_V4.gene.counts",header=F)
ggplot(a,aes(a$V2,a$V4))+geom_point(aes(col=a$V1))+facet_grid(a$V1~.)+geom_smooth(span=0.1)



# Alan wrote a script to make a bed file of gaps in the genome
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
chr1_single_pi$group <- "singleton SNPs"
chr1_shared_pi$group <- "shared SNPs"
chr1_pi <- rbind(chr1_single_pi,chr1_shared_pi)
chr1_fst <- read.table("chr1_fst_100k_allpop.windowed.weir.fst",header=T)

chr2_single_pi<-read.table("chr2.singletons.pi.windowed.pi",header=T)			
chr2_shared_pi<-read.table("chr2.shared.pi.windowed.pi",header=T)			
chr2_single_pi$group <- "singleton SNPs"
chr2_shared_pi$group <- "shared SNPs"
chr2_pi <- rbind(chr2_single_pi,chr2_shared_pi)
chr2_fst <- read.table("chr2_fst_100k_allpop.windowed.weir.fst",header=T)


chr3_single_pi<-read.table("chr3.singletons.pi.windowed.pi",header=T)			
chr3_shared_pi<-read.table("chr3.shared.pi.windowed.pi",header=T)
chr3_single_pi$group <- "singleton SNPs"
chr3_shared_pi$group <- "shared SNPs"
chr3_pi <- rbind(chr3_single_pi,chr3_shared_pi)			
chr3_fst <- read.table("chr3_fst_100k_allpop.windowed.weir.fst",header=T)


chr4_single_pi<-read.table("chr4.singletons.pi.windowed.pi",header=T)			
chr4_shared_pi<-read.table("chr4.shared.pi.windowed.pi",header=T)			
chr4_single_pi$group <- "singleton SNPs"
chr4_shared_pi$group <- "shared SNPs"
chr4_pi <- rbind(chr4_single_pi,chr4_shared_pi)
chr4_fst <- read.table("chr4_fst_100k_allpop.windowed.weir.fst",header=T)


chr5_single_pi<-read.table("chr5.singletons.pi.windowed.pi",header=T)			
chr5_shared_pi<-read.table("chr5.shared.pi.windowed.pi",header=T)
chr5_single_pi$group <- "singleton SNPs"
chr5_shared_pi$group <- "shared SNPs"
chr5_pi <- rbind(chr5_single_pi,chr5_shared_pi)			
chr5_fst <- read.table("chr5_fst_100k_allpop.windowed.weir.fst",header=T)


chrX_single_pi<-read.table("chrX.singletons.pi.windowed.pi",header=T)			
chrX_shared_pi<-read.table("chrX.shared.pi.windowed.pi",header=T)
chrX_single_pi$group <- "singleton SNPs"
chrX_shared_pi$group <- "shared SNPs"
chrX_pi <- rbind(chrX_single_pi,chrX_shared_pi)			
chrX_fst <- read.table("chrX_fst_100k_allpop.windowed.weir.fst",header=T)

# calculate genome wide average Fst, and confidence intervals
fst_all <- rbind(chrX_fst,chr2_fst,chr3_fst,chr4_fst,chr5_fst,chrX_fst)
error <- qt(0.975,df=length(fst_all$WEIGHTED_FST)-1)*sd(fst_all$WEIGHTED_FST)/sqrt(length(fst_all$WEIGHTED_FST))
fst_upper_ci <- mean(fst_all$WEIGHTED_FST)+error
fst_lower_ci <- mean(fst_all$WEIGHTED_FST)-error
quantile<-quantile(chrX_fst$WEIGHTED_FST, c(.05,.95))


chr1_div_plot <- ggplot(chr1_pi,aes(BIN_START,PI,fill=group))+geom_area()+scale_fill_manual(values=c("#b2182b","#851220"))+
               coord_cartesian(xlim = c(0,5.2e7), ylim = c(0,0.01))+
               xlab("")+ylab(expression(pi))+
               theme_classic()+theme(axis.title.x=element_blank(),
               axis.text.x=element_blank(),
               axis.ticks.x=element_blank(),
               legend.position = "none")

chr1_fst_plot <-ggplot(chr1_fst,aes(BIN_START,WEIGHTED_FST))+geom_area(fill="#9C1526")+theme_bw()+
               theme_classic()+
               coord_cartesian(xlim = c(0,5.2e7), ylim = c(0,0.35))+
               xlab("")+ylab("Fst")+
               theme(axis.title.x=element_blank(),
               axis.text.x=element_blank(),
               axis.ticks.x=element_blank())


chr1_gene_plot	<-	ggplot()+
               geom_rect(aes(xmin=1,ymin=-1,xmax=max(chr1$V2),ymax=1),fill="#b2182b",alpha=0.5)+
               geom_rect(aes(xmin=chr1N$V2,ymin=-1,xmax=chr1N$V3,ymax=1),fill="white",linetype=0)+
               geom_linerange(aes(chr1pos$V2,ymin=0,ymax=1),size=0.05,alpha=0.1)+
               geom_linerange(aes(chr1neg$V2,ymin=-1,ymax=0),size=0.05,alpha=0.1)+
               theme_classic()+ scale_y_continuous(breaks=NULL)+
               xlim(0,5.2e7)+ylab("I")+xlab("")

chr2_div_plot <- ggplot(chr2_pi,aes(BIN_START,PI,fill=group))+geom_area()+theme_bw()+scale_fill_manual(values=c("#FC8D59","#BD6A43"))+
                              coord_cartesian(xlim = c(0,5.2e7), ylim = c(0,0.01))+
                              xlab("")+ylab(expression(pi))+
                              theme_classic()+theme(axis.title.x=element_blank(),
                              axis.text.x=element_blank(),
                              axis.ticks.x=element_blank(),
                              legend.position = "none")

chr2_fst_plot <-ggplot(chr2_fst,aes(BIN_START,WEIGHTED_FST))+geom_area(fill="#DD7B4E")+theme_bw()+
                              theme_classic()+
                              coord_cartesian(xlim = c(0,5.2e7), ylim = c(0,0.35))+
                              xlab("")+ylab("Fst")+
                              theme(axis.title.x=element_blank(),
                              axis.text.x=element_blank(),
                              axis.ticks.x=element_blank())

chr2_gene_plot <- ggplot()+
               geom_rect(aes(xmin=1,ymin=-1,xmax=max(chr2$V2),ymax=1),fill="#fc8d59",alpha=0.5)+
               geom_rect(aes(xmin=chr2N$V2,ymin=-1,xmax=chr2N$V3,ymax=1),fill="white",linetype=0)+
               geom_linerange(aes(chr2pos$V2,ymin=0,ymax=1),size=0.05,alpha=0.1)+
               geom_linerange(aes(chr2neg$V2,ymin=-1,ymax=0),size=0.05,alpha=0.1)+
               theme_classic()+
               scale_y_continuous(breaks=NULL)+
               xlim(0,5.2e7)+ylab("II")+xlab("")

chr3_div_plot <- ggplot(chr3_pi,aes(BIN_START,PI,fill=group))+geom_area()+theme_bw()+scale_fill_manual(values=c("#FEE090","#BFA86C"))+
                              coord_cartesian(xlim = c(0,5.2e7), ylim = c(0,0.01))+
                              xlab("")+ylab(expression(pi))+
                              theme_classic()+theme(axis.title.x=element_blank(),
                              axis.text.x=element_blank(),
                              axis.ticks.x=element_blank(),
                              legend.position = "none")

chr3_fst_plot <- ggplot(chr3_fst,aes(BIN_START,WEIGHTED_FST))+geom_area(fill="#DEC47E")+theme_bw()+
                              theme_classic()+
                              coord_cartesian(xlim = c(0,5.2e7), ylim = c(0,0.35))+
                              xlab("")+ylab("Fst")+
                              theme(axis.title.x=element_blank(),
                              axis.text.x=element_blank(),
                              axis.ticks.x=element_blank())

chr3_gene_plot <- ggplot()+
               geom_rect(aes(xmin=1,ymin=-1,xmax=max(chr3$V2),ymax=1),fill="#fee090",alpha=0.5)+
               geom_rect(aes(xmin=chr3N$V2,ymin=-1,xmax=chr3N$V3,ymax=1),fill="white",linetype=0)+
               geom_linerange(aes(chr3pos$V2,ymin=0,ymax=1),size=0.05,alpha=0.1)+
               geom_linerange(aes(chr3neg$V2,ymin=-1,ymax=0),size=0.05,alpha=0.1)+
               theme_classic()+ scale_y_continuous(breaks=NULL)+
               xlim(0,5.2e7)+ylab("III")+xlab("")


chr4_div_plot <- ggplot(chr4_pi,aes(BIN_START,PI,fill=group))+geom_area()+theme_bw()+scale_fill_manual(values=c("#D1E5F0","#9DACB4"))+
                              coord_cartesian(xlim = c(0,5.2e7), ylim = c(0,0.01))+
                              xlab("")+ylab(expression(pi))+
                              theme_classic()+theme(axis.title.x=element_blank(),
                              axis.text.x=element_blank(),
                              axis.ticks.x=element_blank(),
                              legend.position = "none")

chr4_fst_plot <-ggplot(chr4_fst,aes(BIN_START,WEIGHTED_FST))+geom_area(fill="#B7C8D2")+theme_bw()+
                              theme_classic()+
                              coord_cartesian(xlim = c(0,5.2e7), ylim = c(0,0.35))+
                              xlab("")+ylab("Fst")+
                              theme(axis.title.x=element_blank(),
                              axis.text.x=element_blank(),
                              axis.ticks.x=element_blank())
chr4_gene_plot <- ggplot()+
               geom_rect(aes(xmin=1,ymin=-1,xmax=max(chr4$V2),ymax=1),fill="#d1e5f0",alpha=0.5)+
               geom_rect(aes(xmin=chr4N$V2,ymin=-1,xmax=chr4N$V3,ymax=1),fill="white",linetype=0)+
               geom_linerange(aes(chr4pos$V2,ymin=0,ymax=1),size=0.05,alpha=0.1)+
               geom_linerange(aes(chr4neg$V2,ymin=-1,ymax=0),size=0.05,alpha=0.1)+
               theme_classic()+ scale_y_continuous(breaks=NULL)+
               xlim(0,5.2e7)+ylab("IV")+xlab("")



chr5_div_plot <- ggplot(chr5_pi,aes(BIN_START,PI,fill=group))+geom_area()+theme_bw()+scale_fill_manual(values=c("#67A9CF","#4D7F9B"))+
                              coord_cartesian(xlim = c(0,5.2e7), ylim = c(0,0.01))+
                              xlab("")+ylab(expression(pi))+
                              theme_classic()+theme(axis.title.x=element_blank(),
                              axis.text.x=element_blank(),
                              axis.ticks.x=element_blank(),
                              legend.position = "none")

chr5_fst_plot <-ggplot(chr5_fst,aes(BIN_START,WEIGHTED_FST))+geom_area(fill="#5A94B5")+theme_bw()+
                              theme_classic()+
                              coord_cartesian(xlim = c(0,5.2e7), ylim = c(0,0.35))+
                              xlab("")+ylab("Fst")+
                              theme(axis.title.x=element_blank(),
                              axis.text.x=element_blank(),
                              axis.ticks.x=element_blank())

chr5_gene_plot <- ggplot()+
               geom_rect(aes(xmin=1,ymin=-1,xmax=max(chr5$V2),ymax=1),fill="#67a9cf",alpha=0.5)+
               geom_rect(aes(xmin=chr5N$V2,ymin=-1,xmax=chr5N$V3,ymax=1),fill="white",linetype=0)+
               geom_linerange(aes(chr5pos$V2,ymin=0,ymax=1),size=0.05,alpha=0.1)+
               geom_linerange(aes(chr5neg$V2,ymin=-1,ymax=0),size=0.05,alpha=0.1)+
               theme_classic()+ scale_y_continuous(breaks=NULL)+
               xlim(0,5.2e7)+ylab("V")+xlab("")


chrX_div_plot <- ggplot(chrX_pi,aes(BIN_START,PI,fill=group))+geom_area()+theme_bw()+scale_fill_manual(values=c("#4575B4","#345887"))+
               coord_cartesian(xlim = c(0,5.2e7), ylim = c(0,0.01))+
               xlab("")+ylab(expression(pi))+
               theme_classic()+theme(axis.title.x=element_blank(),
               axis.text.x=element_blank(),
               axis.ticks.x=element_blank(),
               legend.position = "none")

chrX_fst_plot <-ggplot(chrX_fst,aes(BIN_START,WEIGHTED_FST))+geom_area(fill="#3C669E")+theme_bw()+
               theme_classic()+
               coord_cartesian(xlim = c(0,5.2e7), ylim = c(0,0.35))+
               xlab("")+ylab("Fst")+
               theme(axis.title.x=element_blank(),
               axis.text.x=element_blank(),
               axis.ticks.x=element_blank())

chrX_gene_plot <- ggplot()+
               geom_rect(aes(xmin=1,ymin=-1,xmax=max(chrX$V2),ymax=1),fill="#4575b4",alpha=0.5)+
               geom_rect(aes(xmin=chrXN$V2,ymin=-1,xmax=chrXN$V3,ymax=1),fill="white",linetype=0)+
               geom_linerange(aes(chrXpos$V2,ymin=0,ymax=1),size=0.05,alpha=0.1)+
               geom_linerange(aes(chrXneg$V2,ymin=-1,ymax=0),size=0.05,alpha=0.1)+
               theme_classic()+ scale_y_continuous(breaks=NULL)+
               xlim(0,5.2e7)+ylab("X")+xlab("")


#patchwork
chr1_div_plot + chr1_fst_plot + (chr1_gene_plot / plot_spacer()) +
chr2_div_plot + chr2_fst_plot + (chr2_gene_plot / plot_spacer()) +
chr3_div_plot + chr3_fst_plot + (chr3_gene_plot / plot_spacer()) +
chr4_div_plot + chr4_fst_plot + (chr4_gene_plot / plot_spacer()) +
chr5_div_plot + chr5_fst_plot + (chr5_gene_plot / plot_spacer()) +
chrX_div_plot + chrX_fst_plot + (chrX_gene_plot / plot_spacer()) + plot_layout(ncol=2,byrow=F)


ggsave("global_diversity_genomewide_chromosome_plots.pdf", height=15, width=15, useDingbats = FALSE)
ggsave("global_diversity_genomewide_chromosome_plots.png", height=15, width=15)
```

```shell
scp sd21@pcs5.internal.sanger.ac.uk:/nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/POPULATION_DIVERSITY/VARIANTS/VCFTOOLS/global_diversity_genomewide_chromosome_plots.*  ~/Documents/workbook/hcontortus_genome/04_analysis

```

# Fst outlier analysis






#- plot top 1% of windows for each chromosome
library(ggplot2)
library(dplyr)

a <- read.table("10k_allpop.windowed.weir.fst",header=F)


pcent <- 0.01

ggplot(a,aes(V2,V5))+
	geom_point(alpha=0.5,size=0.5)+
	geom_point(data=subset(a %>% filter(V1=="hcontortus_chr1_Celeg_TT_arrow_pilon") %>% filter(V5> quantile(V5,prob=1-pcent))),aes(V2,V5),col="red")+
	geom_point(data=subset(a %>% filter(V1=="hcontortus_chr2_Celeg_TT_arrow_pilon") %>% filter(V5> quantile(V5,prob=1-pcent))),aes(V2,V5),col="red")+
	geom_point(data=subset(a %>% filter(V1=="hcontortus_chr3_Celeg_TT_arrow_pilon") %>% filter(V5> quantile(V5,prob=1-pcent))),aes(V2,V5),col="red")+
	geom_point(data=subset(a %>% filter(V1=="hcontortus_chr4_Celeg_TT_arrow_pilon") %>% filter(V5> quantile(V5,prob=1-pcent))),aes(V2,V5),col="red")+
	geom_point(data=subset(a %>% filter(V1=="hcontortus_chr5_Celeg_TT_arrow_pilon") %>% filter(V5> quantile(V5,prob=1-pcent))),aes(V2,V5),col="red")+
	geom_point(data=subset(a %>% filter(V1=="hcontortus_chrX_Celeg_TT_arrow_pilon") %>% filter(V5> quantile(V5,prob=1-pcent))),aes(V2,V5),col="red")+
	facet_grid(.~V1)+
	labs(x="Genomic position", y="Genetic differentiation between populations (Fst)")


chr1_fstoutlier <- a %>% filter(V1=="hcontortus_chr1_Celeg_TT_arrow_pilon") %>% filter(V5> quantile(V5,prob=1-pcent))
chr2_fstoutlier <- a %>% filter(V1=="hcontortus_chr2_Celeg_TT_arrow_pilon") %>% filter(V5> quantile(V5,prob=1-pcent))
chr3_fstoutlier <- a %>% filter(V1=="hcontortus_chr3_Celeg_TT_arrow_pilon") %>% filter(V5> quantile(V5,prob=1-pcent))
chr4_fstoutlier <- a %>% filter(V1=="hcontortus_chr4_Celeg_TT_arrow_pilon") %>% filter(V5> quantile(V5,prob=1-pcent))
chr5_fstoutlier <- a %>% filter(V1=="hcontortus_chr5_Celeg_TT_arrow_pilon") %>% filter(V5> quantile(V5,prob=1-pcent))
chrX_fstoutlier <- a %>% filter(V1=="hcontortus_chrX_Celeg_TT_arrow_pilon") %>% filter(V5> quantile(V5,prob=1-pcent))

fstoutlier <- dplyr::bind_rows(chr1_fstoutlier,chr2_fstoutlier,chr3_fstoutlier,chr4_fstoutlier,chr5_fstoutlier,chrX_fstoutlier)

write.table(fstoutlier,"fstoutlier_10k_top1pc.txt",quote=FALSE, row.names = FALSE, col.names = FALSE, sep="\t")


cut -f1,2,3 fstoutlier_10k_top1pc.txt > fstoutlier_10k_top1pc.bed

bedtools intersect -b fstoutlier_10k_top1pc.bed -a ANNOTATION.gff | grep "mRNA" | cut -f9 | cut -c 4-16 | sort | uniq









# nucleotide diversity summary stats
```
cat *.shared.pi.windowed.pi | grep "hcontortus" | sed 's/chr[12345]/chr/g' | datamash -g1 median 5 q1 5 q3 5
hcontortus_chr_Celeg_TT_arrow_pilon	0.00419839	0.00354622	0.00486901
hcontortus_chrX_Celeg_TT_arrow_pilon	0.00152453	0.00114663	0.00198185

cat *.sin*.pi.windowed.pi | grep "hcontortus" | sed 's/chr[12345]/chr/g' | datamash -g1 median 5 q1 5 q3 5
hcontortus_chr_Celeg_TT_arrow_pilon	0.0012937	0.00106099	0.00156967
hcontortus_chrX_Celeg_TT_arrow_pilon	0.000517191	0.000393264	0.000722199

```



```
# R-3.4.0
# library(ggplot2)
# library(patchwork)		
# 			
# a<-read.table("chrX.singletons.pi.windowed.pi",header=T)			
# b<-read.table("chrX.shared.pi.windowed.pi",header=T)			
# c <- read.table("chrX_fst_100k_allpop.windowed.weir.fst",header=T)			
# 			
# bz <- read.table("../../../../XQTL/04_VARIANTS/XQTL_BZ/XQTL_BZ.merged.fst",header=F)
# bz_1 <- bz[bz$V1=="hcontortus_chrX_Celeg_TT_arrow_pilon",]
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
```

```shell
cd /nfs/users/nfs_s/sd21/lustre118_link/play
cut -f1 caenorhabditis_elegans.PRJNA13758.WBPS11.annotations.gff3 | grep -v "#" |  sort | uniq | while read -r CHR; do awk -v CHR=$CHR '{if($1==CHR && $3=="gene" && $2=="WormBase") print $1,($5+$4)/2,$7,$9}' OFS="\t" caenorhabditis_elegans.PRJNA13758.WBPS11.annotations.gff3 | grep "protein_coding" > ${CHR}.genepos.list;  done
```

```R
R-3.5.0
library(ggplot2)
library(patchwork)

chrX<-read.table("I.genepos.list",header=F)
chr2<-read.table("II.genepos.list",header=F)
chr3<-read.table("III.genepos.list",header=F)
chr4<-read.table("IV.genepos.list",header=F)
chr5<-read.table("V.genepos.list",header=F)
chrX<-read.table("X.genepos.list",header=F)

chrXpos<-chrX[chrX$V3=="+",]
chrXneg<-chrX[chrX$V3=="-",]
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



chrX_gene_plot	<-	ggplot()+
               geom_rect(aes(xmin=1,ymin=-1,xmax=max(chrX$V2),ymax=1),fill="white")+
               geom_linerange(aes(chrXpos$V2,ymin=0,ymax=1),size=0.1,alpha=0.2)+
               geom_linerange(aes(chrXneg$V2,ymin=-1,ymax=0),size=0.1,alpha=0.2)+
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

chrX_gene_plot +
chr2_gene_plot +
chr3_gene_plot +
chr4_gene_plot +
chr5_gene_plot +
chrX_gene_plot +
plot_layout(ncol=1)

```





## Genetic Map V4
```shell
cd /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/GENETICMAP_V4/

mkdir 00_RAW 01_REFERENCE 02_MAPPING 03_VARIANTS 04_ANALYSIS


# curated sample lanes and sample names
/nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/GENETICMAP_V4/00_RAW/samples_lanes_IDs.list

# get raw data

while read NAME LANE; do pf data -t lane -i ${LANE} -f fastq -r -l ./ ; done < samples_lanes_IDs.list

# mapping

cd ../02_MAPPING

screen

# repeat this a few times and it'll kick of multiple mapping runs
while read NAME LANE; do ~sd21/bash_scripts/run_bwamem_splitter ${NAME}_${LANE} /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/GENETICMAP_V4/01_REFERENCE/HAEM_V4_final.chr.fa /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/GENETICMAP_V4/00_RAW/${LANE}_1.fastq.gz /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/GENETICMAP_V4/00_RAW/${LANE}_2.fastq.gz; done < ../00_RAW/samples_lanes_IDs.list &


# once mapping is complete, need to merge
#--- note this was placed in a script and bsubbed
while read name lane; do ls -1 ${name}*merged.sorted.marked.bam > ${name}.bamlist ; samtools merge --threads 7 -c -b ${name}.bamlist ${name}.merge.bam; done < ../00_RAW
/samples_lanes_IDs.list

# clean up
rm *bamlist
rm *marked.bam.bai
rm *marked.bam



# fix read headers
for i in *.bam; do bsub.py --threads 3 10 fix_bam_header "gatk-4.0.3.0 AddOrReplaceReadGroups \
     --INPUT=${i} \
     --OUTPUT=${i%.bam}d.bam \
     --RGID=${i%.merge.bam} \
     --RGLB=${i%.merge.bam} \
     --RGPL=HS --RGPU=WSI --RGSM=${i%.merge.bam}";\
done

# remake indices
for i in *.merged.bam; do \
     bsub.py 1 index_bams "samtools index ${i}";
done

rm *.merge.bam
rm *.merge.bam.bai


# variant celling

ls -1 $PWD/*merged.bam > bams.list

cd ../03_VARIANTS

# specify data using full paths to bam list and fasta refernce sequence
BAM_LIST=/nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/GENETICMAP_V4/02_MAPPING/bams.list
REFERENCE=/nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/GENETICMAP_V4/01_REFERENCE/HAEM_V4_final.chr.fa

mkdir GATK_HC_GVCF
cd GATK_HC_GVCF
mkdir REF_sequences

# make reference seqeunces
cp ${REFERENCE} REF_sequences/REF.fa
samtools faidx REF_sequences/REF.fa
samtools dict REF_sequences/REF.fa > REF_sequences/REF.dict
cp ${BAM_LIST} bam.list
grep ">" REF_sequences/REF.fa | sed -e 's/>//g' > sequences.list

ulimit -c unlimited
#/software/pathogen/external/apps/usr/local/gatk-4.0.3.0/gatk-package-4.0.3.0-local.jar

# make jobs
while read BAM; do \
	n=1
	SAMPLE=$( echo ${BAM} | awk -F '/' '{print $NF}' | sed -e 's/.bam//g' )
	if [ -d "${SAMPLE}_GATK_HC_GVCF_complete" ]; then
	continue
	fi
	mkdir ${SAMPLE}_GATK_HC_GVCF
	mkdir ${SAMPLE}_GATK_HC_GVCF/LOGFILES
	echo "gatk-4.0.3.0 GatherVcfsCloud \\" > ${SAMPLE}_GATK_HC_GVCF/run_gather_${SAMPLE}_gvcf
	while read SEQUENCE; do
	echo -e "java -Xmx15G -jar /software/pathogen/external/apps/usr/local/gatk-4.0.3.0/gatk-package-4.0.3.0-local.jar HaplotypeCaller --input ${BAM} --output ${SAMPLE}_GATK_HC_GVCF/${n}.${SAMPLE}.${SEQUENCE}.tmp.gvcf.gz --reference $PWD/REF_sequences/REF.fa --intervals ${SEQUENCE} --emit-ref-confidence GVCF " > ${SAMPLE}_GATK_HC_GVCF/run_hc_${SAMPLE}.${SEQUENCE}.tmp.job_${n};
	echo -e "--input ${PWD}/${SAMPLE}_GATK_HC_GVCF/${n}.${SAMPLE}.${SEQUENCE}.tmp.gvcf.gz \\" >> ${SAMPLE}_GATK_HC_GVCF/run_gather_${SAMPLE}_gvcf;
	let "n+=1"; done < sequences.list;

	echo -e "--output ${PWD}/${SAMPLE}_GATK_HC_GVCF/${SAMPLE}.gvcf.gz; tabix -p vcf ${PWD}/${SAMPLE}_GATK_HC_GVCF/${SAMPLE}.gvcf.gz" >> ${SAMPLE}_GATK_HC_GVCF/run_gather_${SAMPLE}_gvcf;

	echo -e "rm ${PWD}/${SAMPLE}_GATK_HC_GVCF/*.tmp.* && mv ${PWD}/${SAMPLE}_GATK_HC_GVCF/*.[oe] ${SAMPLE}_GATK_HC_GVCF/LOGFILES && cd ${PWD} && mv ${PWD}/${SAMPLE}_GATK_HC_GVCF ${PWD}/${SAMPLE}_GATK_HC_GVCF_complete" > ${SAMPLE}_GATK_HC_GVCF/run_clean_${SAMPLE};

	chmod a+x ${SAMPLE}_GATK_HC_GVCF/run_*

	# setup job conditions
	JOBS=$( ls -1 ${SAMPLE}_GATK_HC_GVCF/run_hc_* | wc -l )
	ID="U$(date +%s)"

	#submit job array to call variants put scaffold / contig
	bsub -q long -R'span[hosts=1] select[mem>15000] rusage[mem=15000]' -n 6 -M15000 -J GATK_HC_${ID}_[1-$JOBS]%100 -e ${SAMPLE}_GATK_HC_GVCF/GATK_HC_${ID}_[1-$JOBS].e -o ${SAMPLE}_GATK_HC_GVCF/GATK_HC_${ID}_[1-$JOBS].o "./${SAMPLE}_GATK_HC_GVCF/run_hc_${SAMPLE}.*job_\$LSB_JOBINDEX"

	#submit job to gather gvcfs into a single, per sample gvcf
	bsub -q normal -w "done(GATK_HC_${ID}_[1-$JOBS])" -R'span[hosts=1] select[mem>500] rusage[mem=500]' -n 1 -M500 -J GATK_HC_${ID}_gather_gvcfs -e ${SAMPLE}_GATK_HC_GVCF/GATK_HC_${ID}_gather_gvcfs.e -o ${SAMPLE}_GATK_HC_GVCF/GATK_HC_${ID}_gather_gvcfs.o "./${SAMPLE}_GATK_HC_GVCF/run_gather_${SAMPLE}_gvcf"

	# clean up
	bsub -q normal -w "done(GATK_HC_${ID}_gather_gvcfs)" -R'span[hosts=1] select[mem>500] rusage[mem=500]' -n 1 -M500 -J GATK_HC_${ID}_clean -e ${SAMPLE}_GATK_HC_GVCF/GATK_HC_${ID}_clean.e -o ${SAMPLE}_GATK_HC_GVCF/GATK_HC_${ID}_clean.o "./${SAMPLE}_GATK_HC_GVCF/run_clean_${SAMPLE}"

	sleep 1
done < ${BAM_LIST}







#echo "gatk-4.0.3.0 GatherVcfsCloud \\" > run_gather_${SAMPLE}_gvcf
#echo -e "--input ${SAMPLE}_GATK_HC_GVCF/${n}.${SAMPLE}.${SEQUENCE}.g.vcf.tmp \\" >> run_gather_${SAMPLE}_gvcf
#
# echo -e "--output ${SAMPLE}.g.vcf" >> run_gather_${SAMPLE}_gvcf
# chmod a+x run_gather_${SAMPLE}_gvcf

mkdir GATK_HC_MERGED
cd GATK_HC_MERGED


ls -1 *complete/*gz >> gvcf.list

GVCF_LIST=$PWD/gvcf.list
REFERENCE=/lustre/scratch118/infgen/team133/sd21/hc/GENOME/POPULATION_DIVERSITY/HAEM_V4_final.chr.fa

# echo -e "java -Xmx15G -jar /software/pathogen/external/apps/usr/local/gatk-4.0.3.0/gatk-package-4.0.3.0-local.jar CombineGVCFs -R ${REFERENCE} \\" > run_merge_gvcfs
# while read SAMPLE; do
# echo -e "--variant ${SAMPLE} \\" >> run_merge_gvcfs;
#    done < ${GVCF_LIST}
#    echo -e "--output cohort.g.vcf.gz" >> run_merge_gvcfs

# chmod a+x run_merge_gvcfs
# bsub.py --queue hugemem --threads 30 200 merge_vcfs "./run_merge_gvcfs"
# threads make a big difference, even thoguh they are not a parameter in the tool


grep ">" ${REFERENCE} | sed -e 's/>//g' > sequences.list

n=1
while read SEQUENCE; do
echo -e "java -Xmx15G -jar /software/pathogen/external/apps/usr/local/gatk-4.0.3.0/gatk-package-4.0.3.0-local.jar CombineGVCFs -R ${REFERENCE} --intervals ${SEQUENCE} \\" > ${n}.run_merge_gvcfs_${SEQUENCE}
while read SAMPLE; do
echo -e "--variant ${SAMPLE} \\" >> ${n}.run_merge_gvcfs_${SEQUENCE};
   done < ${GVCF_LIST}
   echo -e "--output ${SEQUENCE}.cohort.g.vcf.gz" >> ${n}.run_merge_gvcfs_${SEQUENCE};
   let "n+=1"; done < sequences.list

chmod a+x *run_merge_gvcfs*

for i in *run_merge_gvcfs*; do
bsub.py --queue hugemem --threads 30 200 merge_vcfs "./${i}"; done


# single genotyping call,
# bsub.py --queue hugemem --threads 30 200 genotype_vcfs "java -Xmx15G -jar /software/pathogen/external/apps/usr/local/gatk-4.0.3.0/gatk-package-4.0.3.0-local.jar GenotypeGVCFs \
#    -R /lustre/scratch118/infgen/team133/sd21/hc/GENOME/POPULATION_DIVERSITY/HAEM_V4_final.chr.fa \
#    -V cohort.g.vcf.gz \
#    -O cohort.vcf.gz"


# split each chromosome up into separate jobs, and run genotyping on each individually.   
n=1
while read SEQUENCE; do
echo -e "java -Xmx15G -jar /software/pathogen/external/apps/usr/local/gatk-4.0.3.0/gatk-package-4.0.3.0-local.jar GenotypeGVCFs \
-R /lustre/scratch118/infgen/team133/sd21/hc/GENOME/POPULATION_DIVERSITY/HAEM_V4_final.chr.fa \
-V ${SEQUENCE}.cohort.g.vcf.gz \
--intervals ${SEQUENCE} \
-O ${n}.${SEQUENCE}.cohort.vcf.gz" > run_hc_genotype.${SEQUENCE}.tmp.job_${n};
let "n+=1"; done < sequences.list

chmod a+x run_hc_genotype*

mkdir LOGFILES

	# setup job conditions
	JOBS=$( ls -1 run_hc_* | wc -l )
	ID="U$(date +%s)"

# ln -s /lustre/scratch118/infgen/team133/sd21/hc/GENOME/POPULATION_DIVERSITY/MAPPING/cohort.g.vcf.gz
# ln -s /lustre/scratch118/infgen/team133/sd21/hc/GENOME/POPULATION_DIVERSITY/MAPPING/cohort.g.vcf.gz.tbi

bsub -q yesterday -R'span[hosts=1] select[mem>20000] rusage[mem=20000]' -n 6 -M20000 -J GATK_HC_GENOTYPE_${ID}_[1-$JOBS] -e LOGFILES/GATK_HC_GENOTYPE_${ID}_[1-$JOBS].e -o LOGFILES/GATK_HC_GENOTYPE_${ID}_[1-$JOBS].o "./run_hc_*\$LSB_JOBINDEX"



# bring it all together
ls -1v *.cohort.vcf.gz > vcffiles.fofn

vcf-concat -f vcffiles.fofn | gzip -c >  HCV4_GENETICMAP_GATKHC.raw.vcf.gz


# clean up
rm *tmp*
rm *.cohort.vcf.gz
rm *.cohort.vcf.gz.tbi
rm *.cohort.g.vcf.gz
rm *.cohort.g.vcf.gz.tbi
rm *run_merge_gvcfs*
mv *.[oe] LOGFILES











#--------------------------------------------------------------------------------
# coverage variance berween samples, with particular focus on X chromosome coverage

cd /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/POPULATION_DIVERSITY/GENOME_COVERAGE

# extract coverage data for chr1 and chrX
for i in *100000_window.cov; do data=$( cat $i | grep "chr[1|X]" | datamash median 5 -g 1 q1 5 q3 5 | awk 'BEGIN { ORS = " " } { print }' OFS="\t" ) ; country=$( echo $i | cut -c -3 ); echo ${i%.merged.100000_window.cov} ${country} ${data} ; done > genome_coverage.data




R-3.5.0
library(ggplot2)
library(patchwork)

data<-read.table("genome_coverage.data",header=F)

# coverage plot comparing chr1 and chrX
coverage_plot <- ggplot(data,aes(V4,V8))+
     geom_pointrange(aes(ymin=V9, ymax=V10),alpha=0.2)+
     geom_errorbarh(aes(xmax = V6, xmin = V5, height = 0),alpha=.2)+
     scale_y_log10(limits = c(0.001,100))+
     scale_x_log10(limits = c(0.001,100))+
     geom_point(data=subset(data,(V10-V9)<(V6-V5)),aes(V4,V8),colour="orange",alpha=0.2)+
     geom_point(data=subset(data,(V10-V9)>(V6-V5)),aes(V4,V8),colour="blue",alpha=0.2)+
     labs(x="Log10(median) coverage of 100,000 bp window of Chr 1",y="Log10(median) coverage of 100,000 bp window of Chr X")+
     theme_bw()

variance_plot <- ggplot(data,aes(x=reorder(V2,(V10-V9)/(V6-V5),FUN = median),y=(V10-V9)/(V6-V5)))+
     geom_boxplot()+
     geom_jitter(data=subset(data,(V10-V9)<(V6-V5)),colour="orange",alpha=0.2)+
     geom_jitter(data=subset(data,(V10-V9)>(V6-V5)),colour="blue",alpha=0.2)+
     labs(x="Country",y="Coverage variance between chr X and chr1")+
     theme_bw()


coverage_plot + variance_plot + plot_layout(ncol=2)
ggsave("X-to_autosome_coverage_ratio.pdf", useDingbats=F)






a_plot <-ggplot(data,aes(x=reorder(V2,(V6-V5),FUN = median),y=(V6-V5)))+
     geom_boxplot()+
     geom_jitter(data=subset(data,(V10-V9)<(V6-V5)),colour="orange",alpha=0.2)+
     geom_jitter(data=subset(data,(V10-V9)>(V6-V5)),colour="blue",alpha=0.2)+
     labs(x="Country",y="Coverage variance between chr X and chr1 \n Coverage ratio = (chrX / chr1)")+
     theme_bw()


b_plot <-ggplot(data,aes(x=reorder(V2,(V10-V9),FUN = median),y=(V10-V9)))+
	     geom_boxplot()+
	     geom_jitter(data=subset(data,(V10-V9)<(V6-V5)),colour="orange",alpha=0.2)+
	     geom_jitter(data=subset(data,(V10-V9)>(V6-V5)),colour="blue",alpha=0.2)+
	     labs(x="Country",y="Coverage variance between chr X and chr1 \n Coverage ratio = (chrX[Q3 - Q1] / chr1[Q3 - Q1])")+
	     theme_bw()
a_plot + b_plot + plot_layout(ncol=2)






#----- SNPeff

# setup for HCON_V4

cd /nfs/users/nfs_s/sd21/lustre118_link/software/COMPARATIVE_GENOMICS/snpEff

mkdir data/HCON_V4

cd  data/HCON_V4
ln -s /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/REF/HAEM_V4_final.chr.fa HCON_V4.fa
ln -s /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/TRANSCRIPTOME/TRANSCRIPTOME_CURATION/HCON_V4_WBP11plus_190125.ips.gff3 genes.gff
gffread genes.gff -g HCON_V4.fa -y protein.fa


# modify config file
cd /nfs/users/nfs_s/sd21/lustre118_link/software/snpEff

echo "

# Haemonchus contortus chromosomes V4
HCON_V4_20200130.genome : HCON_V4_20200130

" > new.genome

cat snpEff.config new.genome > tmp; mv tmp snpEff.config


# build database
java -jar snpEff.jar build -v HCON_V4




# run SNPeff on population diversity VCFs to calcukated KnKs



bsub.py 20 vcf_maf_chr1 "vcftools-0.1.14 --gzvcf 1.hcontortus_chr1_Celeg_TT_arrow_pilon.cohort.vcf.gz --maf 0.01 --recode --out chr1_maf0.01"
bsub.py 20 vcf_maf_chr2 "vcftools-0.1.14 --gzvcf 2.hcontortus_chr2_Celeg_TT_arrow_pilon.cohort.vcf.gz --maf 0.01 --recode --out chr2_maf0.01"
bsub.py 20 vcf_maf_chr3 "vcftools-0.1.14 --gzvcf 3.hcontortus_chr3_Celeg_TT_arrow_pilon.cohort.vcf.gz --maf 0.01 --recode --out chr3_maf0.01"
bsub.py 20 vcf_maf_chr4 "vcftools-0.1.14 --gzvcf 4.hcontortus_chr4_Celeg_TT_arrow_pilon.cohort.vcf.gz --maf 0.01 --recode --out chr4_maf0.01"
bsub.py 20 vcf_maf_chr5 "vcftools-0.1.14 --gzvcf 5.hcontortus_chr5_Celeg_TT_arrow_pilon.cohort.vcf.gz --maf 0.01 --recode --out chr5_maf0.01"
bsub.py 20 vcf_maf_chr6 "vcftools-0.1.14 --gzvcf 6.hcontortus_chrX_Celeg_TT_arrow_pilon.cohort.vcf.gz --maf 0.01 --recode --out chrX_maf0.01"

bsub.py --done vcf_maf_chr1 20 snpeff_chr1 "java -Xmx4g -jar /nfs/users/nfs_s/sd21/lustre118_link/software/COMPARATIVE_GENOMICS/snpEff/snpEff.jar -no-downstream -no-intergenic -no-intron -no-upstream -no-utr HCON_V4 chr1_maf0.01.recode.vcf \> chr1.maf0.01.snpeff.vcf.gz"
bsub.py --done vcf_maf_chr2 20 snpeff_chr2 "java -Xmx4g -jar /nfs/users/nfs_s/sd21/lustre118_link/software/COMPARATIVE_GENOMICS/snpEff/snpEff.jar -no-downstream -no-intergenic -no-intron -no-upstream -no-utr HCON_V4 chr2_maf0.01.recode.vcf \> chr2.maf0.01.snpeff.vcf.gz"
bsub.py --done vcf_maf_chr3 20 snpeff_chr3 "java -Xmx4g -jar /nfs/users/nfs_s/sd21/lustre118_link/software/COMPARATIVE_GENOMICS/snpEff/snpEff.jar -no-downstream -no-intergenic -no-intron -no-upstream -no-utr HCON_V4 chr3_maf0.01.recode.vcf \> chr3.maf0.01.snpeff.vcf.gz"
bsub.py --done vcf_maf_chr4 20 snpeff_chr4 "java -Xmx4g -jar /nfs/users/nfs_s/sd21/lustre118_link/software/COMPARATIVE_GENOMICS/snpEff/snpEff.jar -no-downstream -no-intergenic -no-intron -no-upstream -no-utr HCON_V4 chr4_maf0.01.recode.vcf \> chr4.maf0.01.snpeff.vcf.gz"
bsub.py --done vcf_maf_chr5 20 snpeff_chr5 "java -Xmx4g -jar /nfs/users/nfs_s/sd21/lustre118_link/software/COMPARATIVE_GENOMICS/snpEff/snpEff.jar -no-downstream -no-intergenic -no-intron -no-upstream -no-utr HCON_V4 chr5_maf0.01.recode.vcf \> chr5.maf0.01.snpeff.vcf.gz"
bsub.py --done vcf_maf_chr6 20 snpeff_chrX "java -Xmx4g -jar /nfs/users/nfs_s/sd21/lustre118_link/software/COMPARATIVE_GENOMICS/snpEff/snpEff.jar -no-downstream -no-intergenic -no-intron -no-upstream -no-utr HCON_V4 chrX_maf0.01.recode.vcf \> chrX.maf0.01.snpeff.vcf.gz"





grep '^#\|missense_variant\|synonymous_variant' chr1.maf0.01.snpeff.vcf.gz chr2.maf0.01.snpeff.vcf.gz chr3.maf0.01.snpeff.vcf.gz chr4.maf0.01.snpeff.vcf.gz chr5.maf0.01.snpeff.vcf.gz chrX.maf0.01.snpeff.vcf.gz hcon.snpeff.vcf.gz > mis_syn.txt


cat HCON_V4_WBP11plus_190125.ips.gff3 | grep "gene" | cut -f9 | cut -f2 -d ";" | cut -c6-20 > gene_list
echo -e "NAME\tka_count\tka_freq\tks_count\tks_freq\tka_ks" > kaks_out

while read GENE; do bsub.py --queue small 0.1 kaks_snps "./run_kaks_from_vcf_pergene snpeff.vcf.kaks mis_syn.txt ${GENE}" ; done < gene_list
