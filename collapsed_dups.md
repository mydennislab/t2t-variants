# GRCh38 collapsed duplications identification

This repository describes the evaluation of putative problematic regions in GRCh38 and T2T-CHM13 associated with missing or collapsed gene copies.

## 1. Variant calling

```bash
cd /share/dennislab/projects/t2t/variants
```

### 1.1 PacBio HiFi WGS

#### 1.1.1 CHM13 HiFi versus Hg38

Reads were linked from `/share/dennislab/users/dcsoto/ms_pilot/5_pipelines/chm13/chm13_pacbio/reads/`.

Calling variants in GRCh38:
```bash
cd /share/dennislab/projects/t2t/variants/pacbio/chm13_pacbio_vs_hg38

module load gatk/4.1.8.1
conda activate varcalling

/share/dennislab/users/dcsoto/Miniconda3/bin/snakemake \
--snakefile varcalling_pacbio.smk \
--config \
filename="chm13_vs_hg38" \
reference="../../references/hg38.noalt.fa" \
refbasename="../../references/hg38.noalt" \
-p -j 20
```

#### 1.1.2 CHM13 HiFi versus T2T-CHM13

Reads were linked from `/share/dennislab/users/dcsoto/ms_pilot/5_pipelines/chm13/chm13_pacbio/reads/`.

Calling variants in T2T:
```bash
cd /share/dennislab/projects/t2t/variants/pacbio/chm13_pacbio_vs_t2t

module load gatk/4.1.8.1
conda activate varcalling

/share/dennislab/users/dcsoto/Miniconda3/bin/snakemake \
--snakefile varcalling_pacbio.smk \
--config \
filename="chm13_vs_t2t" \
reference="../../references/chm13.draft_v1.0.plusY.fasta" \
refbasename="../../references/chm13.draft_v1.0.plusY" \
-p -j 40
```

### 1.2. CHM13 Illumina WGS

#### 1.2.1 CHM13 Illumina vs GRCh38

Reads were linked from `/share/dennislab/users/dcsoto/ms_pilot/5_pipelines/chm13/chm13_illumina/reads/`.

Reads are 100 bp long.

Calling variants in GRCh38:
```bash
cd /share/dennislab/projects/t2t/variants/illumina/chm13_illumina_vs_hg38

module load gatk/4.1.8.1
conda activate varcalling

/share/dennislab/users/dcsoto/Miniconda3/bin/snakemake \
--snakefile varcalling_illumina.smk \
--config \
filename="chm13_illumina_vs_hg38" \
reference="../../references/hg38.noalt.fa" \
refbasename="../../references/hg38.noalt" \
-p -j 10
```

#### 1.2.2 CHM13 Illumina vs T2T-CHM13

Reads were linked from `/share/dennislab/users/dcsoto/ms_pilot/5_pipelines/chm13/chm13_illumina/reads/`.

Calling variants in T2T:
```bash
cd /share/dennislab/projects/t2t/variants/illumina/chm13_illumina_vs_t2t

module load gatk/4.1.8.1
conda activate varcalling

/share/dennislab/users/dcsoto/Miniconda3/bin/snakemake \
--snakefile varcalling_illumina.smk \
--config \
filename="chm13_illumina_vs_t2t" \
reference="../../references/chm13.draft_v1.0.fasta" \
refbasename="../../references/chm13.draft_v1.0" \
-p -j 10
```

### 1.3. SimReads from T2T-CHM13

```bash
cd /share/dennislab/projects/t2t/variants/simulated
```

Reads length is 150 bp.

Downloading alignments to Hg38:
```bash
globus ls 9db1f0a6-a05a-11ea-8f06-0a21f750d19b:team-variants/read_aligns/bams/chm13_v1.0_simReads_to_hg38

globus transfer 9db1f0a6-a05a-11ea-8f06-0a21f750d19b:team-variants/read_aligns/bams/chm13_v1.0_simReads_to_hg38/chm13_v1.0_simReadsTohg38.bam $(globus endpoint local-id):/share/dennislab/users/dcsoto/other/globus_connect/chm13_v1.0_simReadsTohg38.bam
globus transfer 9db1f0a6-a05a-11ea-8f06-0a21f750d19b:team-variants/read_aligns/bams/chm13_v1.0_simReads_to_hg38/chm13_v1.0_simReadsTohg38.bam.bai $(globus endpoint local-id):/share/dennislab/users/dcsoto/other/globus_connect/chm13_v1.0_simReadsTohg38.bam.bai
```

#### 1.3.1. SimReads from T2T-CHM13 vs GRCh38

Filtering alignments:
```bash
picard AddOrReplaceReadGroups I=chm13_v1.0_simReadsTohg38.bam O=chm13_v1.0_simReadsTohg38.rg.bam RGID=1 RGLB=lib RGPL=SIM RGPU=unit1 RGSM=Sample1
samtools index chm13_v1.0_simReadsTohg38.rg.bam
samtools view -b chm13_v1.0_simReadsTohg38.rg.bam chr{1..22} chrX chrY > chm13_v1.0_simReadsTohg38.rg.no_alt.bam
samtools index chm13_v1.0_simReadsTohg38.rg.no_alt.bam
```

Calling variants in Hg38:
```bash
module load gatk/4.1.8.1

gatk HaplotypeCaller --minimum-mapping-quality 30 -ploidy 2 -R ../references/hg38.noalt.fa -I chm13_v1.0_simReadsTohg38.rg.no_alt.bam -O chm13_v1.0_simReadsTohg38.rg.no_alt.gatk_30.vcf
```

#### 1.3.2 SimReads from T2T-CHM13 vs T2T-CHM13

Mapping and variant calling to T2T:
```bash
samtools sort -@ 30 -n chm13_v1.0_simReadsTohg38.rg.bam > chm13_v1.0_simReadsTohg38.rg.sortn.bam
bedtools bamtofastq -i chm13_v1.0_simReadsTohg38.rg.sortn.bam -fq chm13_v1.0_simReads_1.fastq -fq2 chm13_v1.0_simReads_2.fastq

bwa mem -t 40 ../references/chm13.draft_v1.0.plusY.fasta chm13_v1.0_simReads_1.fastq chm13_v1.0_simReads_2.fastq > chm13_v1.0_simReadsTot2tplusY.sam
samtools view -Sb -@ 30 chm13_v1.0_simReadsTot2tplusY.sam | samtools sort -@ 30 > chm13_v1.0_simReadsTot2tplusY.bam

picard AddOrReplaceReadGroups I=chm13_v1.0_simReadsTot2tplusY.bam O=chm13_v1.0_simReadsTot2tplusY.rg.bam RGID=1 RGLB=lib RGPL=ILLUMINA RGPU=unit1 RGSM=Sample1
samtools index chm13_v1.0_simReadsTot2tplusY.rg.bam

samtools faidx ../references/chm13.draft_v1.0.plusY.fasta
picard CreateSequenceDictionary R=../references/chm13.draft_v1.0.plusY.fasta O=../references/chm13.draft_v1.0.plusY.dict

module load gatk/4.1.8.1
gatk HaplotypeCaller --minimum-mapping-quality 30 -ploidy 2 -R ../references/chm13.draft_v1.0.plusY.fasta -I chm13_v1.0_simReadsTot2tplusY.rg.bam -O chm13_v1.0_simReadsTot2tplusY.rg.vcf
```

## 2. FP hets analyses

### 2.1 Hets distribution

```bash
cd /share/dennislab/projects/t2t/variants/analysis/hets
```

Variant files:
```bash
chm13_pacbio_vs_hg38="/share/dennislab/projects/t2t/variants/pacbio/chm13_pacbio_vs_hg38/results/variants_gatk/chm13_vs_hg38.mm2.primary.rg.gatk_30.vcf"
simulated_vs_hg38="/share/dennislab/projects/t2t/variants/simulated/chm13_v1.0_simReadsTohg38.rg.no_alt.gatk_30.vcf"
chm13_illumina_vs_hg38="/share/dennislab/projects/t2t/variants/illumina/chm13_illumina_vs_hg38/results/variants_gatk/chm13_illumina_vs_hg38.bwa.primary.rg.srt.markdup.gatk_30.vcf"

chm13_pacbio_vs_t2t="/share/dennislab/projects/t2t/variants/pacbio/chm13_pacbio_vs_t2t/results/variants_gatk/chm13_vs_t2t.mm2.primary.rg.gatk_30.vcf"
simulated_vs_t2t="/share/dennislab/projects/t2t/variants/simulated/chm13_v1.0_simReadsTot2tplusY.rg.vcf"
chm13_illumina_vs_t2t="/share/dennislab/projects/t2t/variants/illumina/chm13_illumina_vs_t2t/results/variants_gatk/chm13_illumina_vs_t2t.bwa.primary.rg.srt.markdup.gatk_30.vcf"
```

Obtaining heterozygous SNPs:
```bash
# GRCh38
bcftools view --exclude-types indels ${chm13_pacbio_vs_hg38} | grep "^#\|0/1" > chm13_pacbio_vs_hg38.het_snps.vcf
bcftools view --exclude-types indels ${simulated_vs_hg38} | grep "^#\|0/1" > chm13_v1.0_simReadsTohg38.het_snps.vcf
bcftools view --exclude-types indels ${chm13_illumina_vs_hg38} | grep "^#\|0/1" > chm13_illumina_vs_hg38.het_snps.vcf

# T2T
bcftools view --exclude-types indels ${chm13_pacbio_vs_t2t} | grep "^#\|0/1" > chm13_pacbio_vs_t2t.het_snps.vcf
bcftools view --exclude-types indels ${simulated_vs_t2t} | grep "^#\|0/1" > chm13_v1.0_simReadsTot2t.het_snps.vcf
bcftools view --exclude-types indels ${chm13_illumina_vs_t2t} | grep "^#\|0/1" > chm13_illumina_vs_t2t.het_snps.vcf
```

Merging files (sample names were manually changed):
```bash
# GRCh38
bcftools sort chm13_pacbio_vs_hg38.het_snps.vcf | bgzip -c > chm13_pacbio_vs_hg38.het_snps.vcf.gz
tabix chm13_pacbio_vs_hg38.het_snps.vcf.gz
bcftools sort chm13_v1.0_simReadsTohg38.het_snps.vcf | bgzip -c > chm13_v1.0_simReadsTohg38.het_snps.vcf.gz
tabix chm13_v1.0_simReadsTohg38.het_snps.vcf.gz
bcftools merge chm13_pacbio_vs_hg38.het_snps.vcf.gz chm13_v1.0_simReadsTohg38.het_snps.vcf.gz > merged.het_snps.vcf
awk '{if(($10~"^0/1" && $11~"^0/1")){print}}' merged.het_snps.vcf | wc -l # 57246 shared variants
awk '{if($0~"^#" || ($10~"^0/1" && $11~"^0/1")){print}}' merged.het_snps.vcf > merged.het_snps.both.vcf # identified by both technologies

# T2T
bcftools sort chm13_pacbio_vs_t2t.het_snps.vcf | bgzip -c > chm13_pacbio_vs_t2t.het_snps.vcf.gz
tabix chm13_pacbio_vs_t2t.het_snps.vcf.gz
bcftools sort chm13_v1.0_simReadsTot2t.het_snps.vcf | bgzip -c > chm13_v1.0_simReadsTot2t.het_snps.vcf.gz
tabix chm13_v1.0_simReadsTot2t.het_snps.vcf.gz
bcftools merge chm13_pacbio_vs_t2t.het_snps.vcf.gz chm13_v1.0_simReadsTot2t.het_snps.vcf.gz > merged_toT2T.het_snps.vcf
```

Generating non-overlapping windows:
```bash
# GRCh38
~/scripts/fasta_len.awk /share/dennislab/projects/t2t/variants/references/hg38.noalt.fa > hg38.noalt.len
awk '{print $1"\t"0"\t"$2}' hg38.noalt.len > hg38.noalt.bed
bedtools makewindows -w 1000 -s 1000 -b hg38.noalt.bed > hg38.noalt.windows_1000.bed

# T2T
~/scripts/fasta_len.awk /share/dennislab/projects/t2t/variants/references/chm13.draft_v1.0.fasta > chm13.draft_v1.0.len
awk '{print $1"\t"0"\t"$2}' chm13.draft_v1.0.len > chm13.draft_v1.0.bed
bedtools makewindows -w 1000 -s 1000 -b chm13.draft_v1.0.bed > chm13.draft_v1.0.windows_1000.bed
```

Counting number of variants per window:
```bash
bedtools coverage -a hg38.noalt.windows_1000.bed -b merged.het_snps.vcf > merged.het_snps.windows_1000.txt
bedtools coverage -a chm13.draft_v1.0.windows_1000.bed -b merged_toT2T.het_snps.vcf > merged_toT2T.het_snps.windows_1000.txt
```

> Number of het clusters reported in the manuscript excluded chrY.

### 2.2 Problematic regions in Hg38

```bash
cd /share/dennislab/projects/t2t/variants/analysis/problematic
```

Considering the above distribution, we obtained problematic regions in Hg38 defined as >= 2 kbp with >= 2 het calls.

We calculated problematic regions from variants merged from PB and SimReads:
```bash
awk '{if($4>=2){print}}' ../hets/merged.het_snps.windows_1000.txt | bedtools merge | awk '{ if(($3-$2)>=2000){print}}' > merged.flagged_regions.bed
```

```bash
wc -l merged.flagged_regions.bed # 4178
awk '{sum+=$3-$2}END{print sum}' merged.flagged_regions.bed # 22736000
```

We obtained the distance between regions:
```bash
bedtools closest -io -d -a merged.flagged_regions.bed -b merged.flagged_regions.bed > merged.flagged_regions.distance.txt
```

Based on the distance between problematic regions, we connected those closer than 5 kbp:
```bash
bedtools merge -d 5000 -i merged.flagged_regions.bed > merged.flagged_regions.connected.bed
```

```bash
wc -l merged.flagged_regions.connected.bed # 2848
awk '{sum+=$3-$2}END{print sum}' merged.flagged_regions.connected.bed # 25640000
```

Based on the size distribution we decided to procede with regions larger than a 5 kbp:
```bash
awk '{print $0"\t"$3-$2}' merged.flagged_regions.connected.bed | awk '{if($4>=5000){print}}' | cut -f1-3 > merged.flagged_regions.connected.over_5k.bed
```

```bash
wc -l merged.flagged_regions.connected.over_5k.bed # 921
awk '{sum+=$3-$2}END{print sum}' merged.flagged_regions.connected.over_5k.bed # 21232000
```

We annotated platforms supporting each problematic region:
```bash
bedtools coverage -a merged.flagged_regions.connected.over_5k.bed -b ../hets/chm13_v1.0_simReadsTohg38.het_snps.vcf | cut -f4 > tmp1
bedtools coverage -a merged.flagged_regions.connected.over_5k.bed -b ../hets/chm13_pacbio_vs_hg38.het_snps.vcf | cut -f4 > tmp2
paste merged.flagged_regions.connected.over_5k.bed tmp1 tmp2 | awk '{if($4>=1 && $5>=1){print $1"\t"$2"\t"$3"\tSimReads,PacBio"}else if($4>=1 && $5<1){print $1"\t"$2"\t"$3"\tSimReads"}else if($4<1 && $5>=1){print $1"\t"$2"\t"$3"\tPacBio"}}' > merged.flagged_regions.connected.over_5k.platform.bed
```

We then obtained a VCF of the het calls within the 921 problematic regions:
```bash
bedtools intersect -header -a ../hets/merged.het_snps.vcf -b merged.flagged_regions.connected.over_5k.bed > merged.het_snps.flagged_regions.connected.over_5k.vcf
grep -v "^#" merged.het_snps.flagged_regions.connected.over_5k.vcf | wc -l # 299923
```

> Again, the file shared in the manuscript had chrY removed.

### 2.3 Excess heterozygosity

```bash
cd /share/dennislab/projects/t2t/variants/analysis/gnomad
```

First, we selected only SNPs in InbreedingCoeff file:
```bash
bcftools view --exclude-types indels gnomad.genomes.r3.0.sites.InbreedingCoeff.vcf > gnomad.genomes.r3.0.sites.InbreedingCoeff.snps.vcf
bgzip -c gnomad.genomes.r3.0.sites.InbreedingCoeff.snps.vcf > gnomad.genomes.r3.0.sites.InbreedingCoeff.snps.vcf.gz
tabix gnomad.genomes.r3.0.sites.InbreedingCoeff.snps.vcf.gz
```

**Direct comparison of variants**

We directly compared InbreedingCoeff variants with all het calls:
```bash
bgzip ../hets/merged.het_snps.vcf
tabix ../hets/merged.het_snps.vcf.gz

bcftools isec -p isec_output -Ov ../hets/merged.het_snps.vcf.gz gnomad.genomes.r3.0.sites.InbreedingCoeff.snps.vcf.gz
for file in isec_output/*vcf; do grep -v "^#" $file | wc -l ; done
```

- Private to merged.het_snps.vcf.gz: 288891 (without Y: 281569)
- Private to gnomad.genomes.r3.0.sites.InbreedingCoeff.snps.vcf.gz: 205460 (without Y: 201487)
- Shared records: 87877 (without Y: 87005)
- 23.32% of variants are flagged with excess of heterozygosity

We repeated the analysis using variants identified by both platforms only:
```bash
bgzip ../hets/merged.het_snps.both.vcf
tabix ../hets/merged.het_snps.both.vcf.gz

bcftools isec -p isec_output_both -Ov ../hets/merged.het_snps.both.vcf.gz gnomad.genomes.r3.0.sites.InbreedingCoeff.snps.vcf.gz
for file in isec_output_both/*vcf; do grep -v "^#" $file | wc -l ; done
```

- Private to merged.het_snps.both.vcf.gz: 20587
- Private to gnomad.genomes.r3.0.sites.InbreedingCoeff.snps.vcf.gz: 256678
- Shared records: 36659
- 64.04% of variants identified by both technologies are flagged with excess of heterozygosity

Finally, we repeated the analysis with het calls within delineated problematic regions only:
```bash
bgzip ../problematic/merged.het_snps.flagged_regions.connected.over_5k.vcf
tabix ../problematic/merged.het_snps.flagged_regions.connected.over_5k.vcf.gz

bcftools isec -p isec_output_problematic -Ov ../problematic/merged.het_snps.flagged_regions.connected.over_5k.vcf.gz gnomad.genomes.r3.0.sites.InbreedingCoeff.snps.vcf.gz
for file in isec_output_problematic/*vcf; do grep -v "^#" $file | grep -v "chrY" | wc -l ; done
```

- Private to merged.het_snps.flagged_regions.connected.over_5k.vcf.gz: 214480
- Private to gnomad.genomes.r3.0.sites.InbreedingCoeff.snps.vcf.gz: 210458
- Shared records: 78034

**Enrichment of InbreedingCoeff**

Then, we answered if there were more InbreedingCoeff variants than expected by random chance inside problematic regions.

Then, we calculated the number of observed InbreedingCoeff variants:
```bash
bedtools intersect -a gnomad.genomes.r3.0.sites.InbreedingCoeff.snps.vcf -b ../problematic/merged.flagged_regions.connected.over_5k.platform.bed | cut -f1,2 | sort | uniq | wc -l # 134842
```

We then calculated the number of these variants within random regions of similar size:
```bash
for num in `seq 1 10000`; do bedtools shuffle -i ../problematic/merged.flagged_regions.connected.over_5k.bed -g hg38.genome -excl gaps.hg38.noalt.bed  > random/random_${num}.bed; done

seq 1 10000 | xargs -n1 -P10 bash -c 'bedtools intersect -a gnomad.genomes.r3.0.sites.InbreedingCoeff.snps.vcf -b random/random_$0.bed | cut -f1,2 | sort | uniq | wc -l' > random_InbreedingCoeff.txt
```

### 2.4 Known GRC issues

```bash
cd /share/dennislab/projects/t2t/variants/analysis/issues
```

Downloading last GRCh38 issues file:
```bash
globus transfer 9db1f0a6-a05a-11ea-8f06-0a21f750d19b:/team-variants/grch38_issues/hg38.parsedissues.bed $(globus endpoint local-id):hg38.parsedissues.bed
cut -f1-3 hg38.parsedissues.bed | bedtools merge > hg38.parsedissues.ideogram.bed
```

There were 418 known issues of Hg38.

```bash
module load R/4.0.1

bedtools intersect -wao -a ../problematic/merged.flagged_regions.connected.over_5k.bed -b hg38.parsedissues.bed | \
awk 'BEGIN{OFS="\t"}{if($12>0){print $1,$2,$3,$9","$10}else{print $1,$2,$3"\tNone"}}'| \
/software/R/4.0.1/lssc0-linux/bin/Rscript ~/scripts/collapse_bed.R /dev/stdin merged.flagged_regions.connected.over_5k.hg38_issue.bed
```

We found that 344 out of 921 problematic regions have been previously reported as problematic.

### 2.5 Lifting to T2T

```bash
cd /share/dennislab/projects/t2t/variants/analysis/lifting
```

**LiftOver**

> Old assembly (Target). New assembly (Query).

Lifting Hg38 connected problematic regions to T2T coordinates using liftover:
```bash
wget http://t2t.gi.ucsc.edu/chm13/hub/t2t-chm13-v1.0/hg38Lastz/hg38.t2t-chm13-v1.0.over.chain.gz
wget http://t2t.gi.ucsc.edu/chm13/hub/t2t-chm13-v1.0/hg38Lastz/hg38.t2t-chm13-v1.0.all.chain.gz
gunzip hg38.t2t-chm13-v1.0.over.chain.gz
gunzip hg38.t2t-chm13-v1.0.all.chain.gz
```

Lifting with LiftOver:
```bash
# we added the original coordinates as an extra column
liftOver <(awk '{print $0"\t"$1":"$2"-"$3}' ../problematic/merged.flagged_regions.connected.over_5k.bed) hg38.t2t-chm13-v1.0.over.chain merged.flagged_regions.connected.over_5k.liftover_lifted.bed merged.flagged_regions.connected.over_5k.lifover_unmapped.bed

module load R/4.0.1
bedtools intersect -wao -a ../problematic/merged.flagged_regions.connected.over_5k.bed -b <(awk 'BEGIN{OFS="\t"}{print $4,$1":"$2"-"$3}' merged.flagged_regions.connected.over_5k.liftover_lifted.bed | awk 'BEGIN{OFS="\t"}{gsub(":","\t",$1); gsub("-","\t",$1); print}') | \
awk 'BEGIN{OFS="\t"}{if($8>0){print $1,$2,$3,$7}else{print $1,$2,$3"\tUnmapped"}}' | \
/software/R/4.0.1/lssc0-linux/bin/Rscript ~/scripts/collapse_bed.R /dev/stdin merged.flagged_regions.connected.over_5k.liftover.bed
```

We were unable to lift 623/921 regions using this approach.

**Minimap2**

Lifting with Minimap2:
```bash
bedtools getfasta -fi ../../references/hg38.noalt.fa -bed ../problematic/merged.flagged_regions.connected.over_5k.bed > merged.flagged_regions.connected.over_5k.fasta
minimap2 -x asm5 -t 25 ../../references/chm13.draft_v1.0.fasta merged.flagged_regions.connected.over_5k.fasta > merged.flagged_regions.connected.over_5k.mm2.paf

cat merged.flagged_regions.connected.over_5k.mm2.paf | cut -f1,6,8,9,12 > minimap2_hits.txt
```

All problematic regions have at least one hit with minimap2.

**CHM13 coordinates**

We processed hits following a "decision tree":
```r
library(tidyverse)

# Loading mm2 hits
mm2 <- read.table("minimap2_hits.txt", sep="\t", col.names=c("coords","chr_hit","start_hit","end_hit","mapq"))
mm2 <- mm2 %>% mutate(coords_hit=paste(chr_hit,":",start_hit,"-",end_hit, sep=""))

# Annotating liftover coordinates
liftover <- read.table("merged.flagged_regions.connected.over_5k.liftover.bed") %>%
  unite("coords",V1,V2,sep=":") %>% unite("coords",coords,V3,sep="-") %>% rename(liftover = V4)
data <- list(mm2, liftover) %>% reduce(left_join, by = "coords")

# Adding sizes to each region
data_sizes <- data %>%
  mutate(coords2=coords) %>%
  separate(coords2, c("chr_hg38","temp"), sep=":") %>%
  separate(temp, c("start_hg38","end_hg38"), sep="-") %>%
  mutate(start_hg38=as.numeric(start_hg38), end_hg38=as.numeric(end_hg38)) %>%
  mutate(size_hg38=end_hg38-start_hg38) %>%
  mutate(size_hit=end_hit-start_hit) %>%
  mutate(size_hit_percent=(size_hit/size_hg38)*100) %>%
  mutate(liftover2=liftover) %>%
  mutate(liftover2=ifelse(liftover2=="Unmapped",NA,liftover2)) %>%
  separate(liftover2, c("chr_liftover","temp"), sep=":") %>%
  separate(temp, c("start_liftover","end_liftover"), sep="-") %>%
  mutate(start_liftover=as.numeric(start_liftover), end_liftover=as.numeric(end_liftover)) %>%
  mutate(size_liftover=end_liftover-start_liftover) %>%
  mutate(size_liftover_percent=(size_liftover/size_hg38)*100)

# Annotating hit distance to original coordinates
# If coordinates are in different chromosomes then, "NA"
data_dist <- data_sizes %>%
 mutate(hit_dist=ifelse(chr_hg38==chr_hit, abs((as.numeric(start_hg38)+as.numeric(end_hg38))/2 - (start_hit+end_hit)/2), NA))

data <- data_dist %>% select(coords, liftover, coords_hit, size_hg38, size_liftover, size_liftover_percent, size_hit, size_hit_percent, hit_dist)

# Selecting LiftOVer hits within 80 to 120 percent of the original size
set1 <- data %>%
    filter(liftover!="Unmapped" & size_liftover_percent>=80 & size_liftover_percent<=120) %>%
    select(coords, liftover) %>% distinct() %>% rename("CHM13_coords"=liftover)

# Selecting entries with mm2 hit within 80-120 percent of the original size.
# If more than one good hit was found by Minimap2, then we chose the "closest one".
set2 <- data %>%
    filter(liftover=="Unmapped" | (size_liftover_percent<=80 | size_liftover_percent>=120)) %>%
    filter(size_hit_percent>=80 & size_hit_percent<=120) %>%
    group_by(coords) %>%
    mutate(min_hit_dist=min(hit_dist)) %>%
    filter(hit_dist==min_hit_dist) %>%
    select(coords, coords_hit) %>%
    distinct() %>% rename("CHM13_coords"=coords_hit)


lifted_coords <- rbind(set1,set2) %>% pull(coords)
hits_tocheck <- data %>% filter(!coords %in% lifted_coords)

write.table(rbind(set1,set2), "new_lifted_coords.tsv", quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
write.table(hits_tocheck, "minimap2_hits_toCheck.tsv", quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
```

New lifted coordinates (Hg38->CHM13) were saved in `merged.flagged_regions.connected.over_5k.lifted_toT2T_v2.bed`

The new lifted coordinates (CHM13 only) were saved in `CHM13_coords.v2.bed`

### 2.6 Copy number

```bash
cd /share/dennislab/projects/t2t/variants/analysis/cn
```

**WSSD subset in CHM13 reference**

Copy-number estimates from the Eichler lab were downloaded from the browser to `/share/dennislab/projects/t2t/wssd/hgdp_sub/` and renamed it with:
```bash
for f in *; do mv "$f" `echo $f | tr ' ' '_'`; done
```

Converting bigwig to bed:
```bash
for file in /share/dennislab/projects/t2t/wssd/hgdp_sub/*; do
  filename=$(basename $file | cut -d"_" -f1-3 | cut -d"-" -f1)
  tail -n +2 $file | cut -f1-3,10 > input/$filename.CN.bed
done
```

**WSSD complete dataset in CHM13 reference**

We also downloaded the complete dataset from the ftp to `/share/dennislab/projects/t2t/wssd/hgdp_complete/` and converted them to bed:
```bash
ls /share/dennislab/projects/t2t/wssd/hgdp_complete/*.bed | cut -d"/" -f8 | cut -d"." -f1 > samples.txt
cat samples.txt | xargs -n1 -P20 bash -c 'cut -f1-3,10 /share/dennislab/projects/t2t/wssd/hgdp_complete/$0.bed > input/$0.CN.bed'
```

From this dataset, we found that sample `LP6005442-DNA_A08_wssd.CN.bed` had strange copy numbers compared to the rest of the samples, so we removed it from input 
ory.

**Genotyping**

Genotyping Hg38 problematic regions in CHM13:
```bash
/share/dennislab/users/dcsoto/Miniconda3/envs/igc/bin/python genotype_cn_parallel.py -t 15 --path input \
--genes ../lifting/CHM13_coords.v1.bed --output CHM13_coords.v1.CN.bed
```

```bash
/share/dennislab/users/dcsoto/Miniconda3/envs/igc/bin/python genotype_cn_parallel.py -t 40 --path input \
--genes ../lifting/CHM13_coords.v2.bed --output CHM13_coords.v2.CN.bed
```

Adding stats (min/max/median):
```r
# module load R/4.0.1
library(tidyverse)

df <- read.table("CHM13_coords.v2.CN.bed", header = TRUE)

df_stats_1 <- df %>% select(!starts_with("hg38")) %>% select(!starts_with("chm13")) %>%
  gather(key = "sample", value = "cn", -chrom, -chromStart, -chromEnd) %>%
  group_by(chrom, chromStart, chromEnd) %>% summarise(Mean_Pop_CN = mean(cn), Median_Pop_CN=median(cn), Min_Pop_CN=min(cn), Max_Pop_CN=max(cn)) %>%
  unite("CHM13_coords", chrom, chromStart, chromEnd)

df_stats_2 <- df %>% select(chrom, chromStart, chromEnd, chm13_wssd) %>%
  unite("CHM13_coords", chrom, chromStart, chromEnd)

df_stats_3 <- df %>% select(!starts_with("chm13")) %>%
  gather(key = "sample", value = "cn", -chrom, -chromStart, -chromEnd, -hg38_wssd) %>%
  group_by(chrom, chromStart, chromEnd, hg38_wssd) %>% summarise(count=sum(cn<hg38_wssd)) %>%
  unite("CHM13_coords", chrom, chromStart, chromEnd)

joined <- list(df_stats_1, df_stats_2, df_stats_3) %>%
  reduce(left_join, by = "CHM13_coords") %>%
  rename(CHM13_WSSD=chm13_wssd, Hg38_WSSD=hg38_wssd, Below_Hg38=count)

write.table(joined, "CHM13_coords.v2.CN.stats.bed", quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
```

### 2.7 Gene annotations

```bash
cd /share/dennislab/projects/t2t/variants/analysis/genes
```

**Gencode v35**

GTF to bed:
```bash
wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_35/gencode.v35.annotation.gff3.gz
gunzip gencode.v35.annotation.gff3.gz
gff2bed < gencode.v35.annotation.gff3 | awk '{if($8=="gene"){print}}' |  cut -f1-3,10 | sed 's/;/\t/g' | cut -f1-3,6,7 | sed 's/=/\t/g' | awk '{print $1"\t"$2"\t"$3"\t"$7";"$5}' > gencode.v35.annotation.genesOnly.bed
```

Genes within problematic regions in Hg38:
```bash
module load R/4.0.1

# All gene features
bedtools intersect -wao -a ../problematic/merged.flagged_regions.connected.over_5k.bed -b <(sed 's/;/,/g' gencode.v35.annotation.genesOnly.bed) | \
awk 'BEGIN{OFS="\t"}{if($8>0){print $1,$2,$3,$7}else{print $1,$2,$3"\tNone"}}'| \
/software/R/4.0.1/lssc0-linux/bin/Rscript ~/scripts/collapse_bed.R /dev/stdin merged.flagged_regions.connected.over_5k.gencode_v35.bed
```

**Medically relevant**

Medically relevant genes in roblematic regions in Hg38:
```bash
module load R/4.0.1

bedtools intersect -wao -a ../problematic/merged.flagged_regions.connected.over_5k.bed -b GRCh38_ENSEMBL_genes_biomart_version_100_Medical_Gene_Coordinates.bed | \
awk 'BEGIN{OFS="\t"}{if($8>0){print $1,$2,$3,$7}else{print $1,$2,$3"\tNone"}}'| \
/software/R/4.0.1/lssc0-linux/bin/Rscript ~/scripts/collapse_bed.R /dev/stdin merged.flagged_regions.connected.over_5k.medical.bed
```

**Genes intolerant to loss of function**

To further refine the impact of these problematic regions, we obtained how many of the impacted genes are loss of function intolerant and, therefore, likely functional.

```bash
wget https://storage.googleapis.com/gnomad-public/legacy/exac_browser/forweb_cleaned_exac_r03_march16_z_data_pLI_CNV-final.txt.gz
gunzip -c forweb_cleaned_exac_r03_march16_z_data_pLI_CNV-final.txt.gz | cut -f2,20 > gene2pli.txt
awk '{if($2>=0.9){print}}' gene2pli.txt | cut -f1 | tail -n +2 > pli_scores_09.txt
```

Obtaining pLI≥0.9 genes impacted by problematic regions:
```bash
bedtools intersect -wao -a ../problematic/merged.flagged_regions.connected.over_5k.bed -b gencode.v36.annotation.genesOnly.bed | \
grep -Fwf pli_scores_09.txt > pli_scores_09.impacted.txt
```

**CHM13 gene annotations**

```bash
gff2bed < /share/dennislab/projects/t2t/annotation/cat_v4/CHM13.combined.v4.gff3 | awk '{if($8=="gene"){print}}' | cut -f1-3,10 | sed 's/;/\t/g' | cut -f1-3,4,6 | sed 's/=/\t/g' | awk '{print $1"\t"$2"\t"$3"\t"$5";"$7}' > CHM13.combined.v4.genesOnly.bed
```

### 2.8 Evaluating problematic regions in T2T-CHM13

```bash
cd /share/dennislab/projects/t2t/variants/analysis/fixed
```

**Identifying problematic regions**

First, we obtained problematic regions in T2T:
```bash
awk '{if($4>=2){print}}' ../hets/merged_toT2T.het_snps.windows_1000.txt | bedtools merge | awk '{ if(($3-$2)>=2000){print}}' > problematic_T2T.bed
wc -l problematic_T2T.bed # 2724
awk '{sum+=$3-$2}END{print sum}' problematic_T2T.bed # 15231000
```

We connected those closer than a threshold:
```bash
bedtools merge -d 5000 -i problematic_T2T.bed > problematic_T2T.connected.bed
```

We proceded with regions larger than a size threshold:
```bash
awk '{print $0"\t"$3-$2}' problematic_T2T.connected.bed | awk '{if($4>=5000){print}}' | cut -f1-3 > problematic_T2T.connected.over_5k.bed
```

We annotated platforms supporting each problematic region:
```bash
bedtools coverage -a problematic_T2T.connected.over_5k.bed -b ../hets/chm13_v1.0_simReadsTot2t.het_snps.vcf | cut -f4 > tmp1
bedtools coverage -a problematic_T2T.connected.over_5k.bed -b ../hets/chm13_pacbio_vs_t2t.het_snps.vcf | cut -f4 > tmp2
paste problematic_T2T.connected.over_5k.bed tmp1 tmp2 | \
awk '{if($4>=2 && $5>=2){print $1"\t"$2"\t"$3"\tSimReads,PacBio"}else if($4>=2 && $5<2){print $1"\t"$2"\t"$3"\tSimReads"}else if($4<2 && $5>=2){print $1"\t"$2"\t"$3"\tPacBio"}}' > problematic_T2T.connected.over_5k.platform.bed
```

**Intersecting with Hg38 problematic regions**

Hg38 problematic regions and CHM13 problematic regions:
```bash
bedtools intersect -wao -a ../lifting/../lifting/CHM13_coords.v2.bed -b problematic_T2T.connected.over_5k.platform.bed | \
awk 'BEGIN{OFS="\t"}{if($8>0){print $1,$2,$3,"Yes"}else{print $1,$2,$3"\tNo"}}' | \
/software/R/4.0.1/lssc0-linux/bin/Rscript ~/scripts/collapse_bed.R /dev/stdin CHM13_coords.v2.T2T_problematic.bed
```

Hg38 problematic regions and CHM13 known issues:
```bash
bedtools intersect -wao -a ../lifting/merged.flagged_regions.connected.over_5k.lifted_toT2T.CHM13_coords.bed -b chm13_issues.bed | \
awk 'BEGIN{OFS="\t"}{if($13>0){print $1,$2,$3,$4":"$5"-"$6","$7}else{print $1,$2,$3"\tNo"}}' | \
/software/R/4.0.1/lssc0-linux/bin/Rscript ~/scripts/collapse_bed.R /dev/stdin merged.flagged_regions.connected.over_5k.lifted_toT2T.CHM13_coords.CHM13_issues.bed
```

### 2.9 Other annotations

#### 2.9.1  Alt haplotypes (GRCh38)

```bash
cd /share/dennislab/projects/t2t/variants/analysis/alt
```

Obtaining alt coordinates:
```bash
wget https://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/altLocations.txt.gz
gunzip altLocations.txt.gz
cut -f2- altLocations.txt > altLocations.bed
```

```bash
module load R/4.0.1

bedtools intersect -wao -a ../problematic/merged.flagged_regions.connected.over_5k.bed -b altLocations.bed | \
awk 'BEGIN{OFS="\t"}{if($8>0){print $1,$2,$3,$7}else{print $1,$2,$3"\tNone"}}'| \
/software/R/4.0.1/lssc0-linux/bin/Rscript ~/scripts/collapse_bed.R /dev/stdin merged.flagged_regions.connected.over_5k.alt.bed
```

#### 2.9.2 Repeat masker (GRCh38)

```bash
cd /share/dennislab/projects/t2t/variants/analysis/repeat
```

We obtained RepeatMasker annotations from UCSC table browser:
```bash
awk 'BEGIN{OFS="\t"}{print $6,$7,$8,$11","$12}' repeatmasker_hg38.txt | tail -n +2 > repeatmasker_hg38.bed
```

```bash
bedtools intersect -f 0.5 -wao -a ../problematic/merged.flagged_regions.connected.over_5k.bed -b repeatmasker_hg38.bed | \
awk 'BEGIN{OFS="\t"}{if($8>0){print $1,$2,$3,$7}else{print $1,$2,$3"\tNone"}}'| \
/software/R/4.0.1/lssc0-linux/bin/Rscript ~/scripts/collapse_bed.R /dev/stdin merged.flagged_regions.connected.over_5k.repeatmasker.bed
```

#### 2.9.3 CenSat (T2T-CHM13)

```bash
cd /share/dennislab/projects/t2t/variants/analysis/censat
```

We intersected lifted coordinates with centromeres:
```bash
module load R/4.0.1

bedtools intersect -wao -a ../lifting/CHM13_coords.v2.bed -b /share/dennislab/projects/t2t/coordinates/cenSats.bed | \
awk 'BEGIN{OFS="\t"}{if($13>0){print $1,$2,$3"\tcenSat"}else{print $1,$2,$3"\tNo"}}' | \
/software/R/4.0.1/lssc0-linux/bin/Rscript ~/scripts/collapse_bed.R /dev/stdin CHM13_coords.v2.CenSat.bed
```

#### 2.9.4 SegDups overlap

```bash
cd /share/dennislab/projects/t2t/variants/analysis/segdups
```

**GRCh38 reference**

We evaluated overlap with known SegDups:
```bash
# SegDup file
tail -n +2 genomicSuperDups | cut -f 2-4 | bedtools sort | bedtools merge > genomicSuperDups.bed

# Overlap
bedtools intersect -wo -a ../problematic/merged.flagged_regions.connected.over_5k.bed -b genomicSuperDups.bed | sort | uniq | wc -l # 769
bedtools intersect -a ../problematic/merged.flagged_regions.connected.over_5k.bed -b genomicSuperDups.bed | awk '{sum+=$3-$2}END{print sum}' # 16,166,390
```

- 769/921*100 = 83.5%

```bash
module load R/4.0.1

bedtools intersect -wao -a ../problematic/merged.flagged_regions.connected.over_5k.bed -b genomicSuperDups.bed | \
awk 'BEGIN{OFS="\t"}{if($7>0){print $1,$2,$3"\tSegDup("$4":"$5"-"$6")"}else{print $1,$2,$3"\tNo"}}' | \
/software/R/4.0.1/lssc0-linux/bin/Rscript ~/scripts/collapse_bed.R /dev/stdin merged.flagged_regions.connected.over_5k.segdup.bed
```

**T2T CHM13 reference**

We, then, overlapped the lifted CHM13 coordinates with SegDup annotations in CHM13:
```bash
# Downloading SegDups
wget http://t2t.gi.ucsc.edu/chm13/hub/t2t-chm13-v1.0/sedefSegDups/chm13.draft_v1.0_plus38Y.SDs.bed.bb
bigBedToBed chm13.draft_v1.0_plus38Y.SDs.bed.bb chm13.draft_v1.0_plus38Y.SDs.bed
cut -f1-3 chm13.draft_v1.0_plus38Y.SDs.bed | bedtools sort | bedtools merge > chm13.draft_v1.0_plus38Y.SDs.merged.bed

# Overlap
module load R/4.0.1

bedtools intersect -wao -a ../lifting/CHM13_coords.v2.bed -b chm13.draft_v1.0_plus38Y.SDs.merged.bed | \
awk 'BEGIN{OFS="\t"}{if($7>0){print $1,$2,$3"\tSegDup"}else{print $1,$2,$3"\tNo"}}' | \
/software/R/4.0.1/lssc0-linux/bin/Rscript ~/scripts/collapse_bed.R /dev/stdin CHM13_coords.v2.SDs.bed
```

**Structurally variable SegDup**

```bash
bedtools intersect -a chm13.draft_v1.0_plus38Y.SDs.merged.bed -b ../synteny/chm13.draft_v1.0_plus38Y.no_snyteny_1Mbp.bed > structurally_variable_SDs.t2t.bed
```

```bash
module load R/4.0.1

bedtools intersect -wao -a ../lifting/CHM13_coords.v2.bed -b structurally_variable_SDs.t2t.bed | \
awk 'BEGIN{OFS="\t"}{if($7>0){print $1,$2,$3"\tSD_diff"}else{print $1,$2,$3"\tNo"}}' | \
/software/R/4.0.1/lssc0-linux/bin/Rscript ~/scripts/collapse_bed.R /dev/stdin CHM13_coords.v2.SD_diff.bed
```

#### 2.9.5 Non-syntenic regions

```bash
cd /share/dennislab/projects/t2t/variants/analysis/synteny
```

There are multiple definitions of non-syntenic regions.

- 1Mbp Mitchell's non-syntenic regions
- Sergey's wm non-syntenic regions
- ns_wm (`bedtools intersect -a <NON-SYNTENIC BED> -b <WINNOWMAP BED>`)
- nw_wm_100mer (`bedtools intersect -a ns_wm.bed -b <MAPPABLE_REGIONS BED>`)

```bash
module load R/4.0.1

bedtools intersect -wao -a ../lifting/CHM13_coords.v2.bed -b chm13.draft_v1.0_plus38Y.no_snyteny_1Mbp.bed | \
awk 'BEGIN{OFS="\t"}{if($7>0){print $1,$2,$3"\tNon-syntenic"}else{print $1,$2,$3"\tNo"}}' | \
/software/R/4.0.1/lssc0-linux/bin/Rscript ~/scripts/collapse_bed.R /dev/stdin CHM13_coords.v2.ns_1mb.bed

bedtools intersect -wao -a ../lifting/CHM13_coords.v2.bed -b chm13_v1.0_uncoveredByGRCh38WinnowmapAlignments.bed | \
awk 'BEGIN{OFS="\t"}{if($7>0){print $1,$2,$3"\tNon-syntenic"}else{print $1,$2,$3"\tNo"}}' | \
/software/R/4.0.1/lssc0-linux/bin/Rscript ~/scripts/collapse_bed.R /dev/stdin CHM13_coords.v2.ns_wm.bed
```

#### 2.9.6 
heterozygosity (T2T)

```bash
cd /share/dennislab/projects/t2t/variants/analysis/inbreeding
```

```bash
bcftools view --max-alleles 2 --exclude-types indels -i 'INFO/InbreedingCoeff<-0.3' /share/dennislab/databases/data/1KG_highcov_t2t/1kgp.recalibrated.snp_indel.pass.nogenos.vcf.gz > 1kgp.recalibrated.snp_indel.pass.nogenos.snps.InbreedingCoeff.vcf
bgzip 1kgp.recalibrated.snp_indel.pass.nogenos.snps.InbreedingCoeff.vcf
tabix 1kgp.recalibrated.snp_indel.pass.nogenos.snps.InbreedingCoeff.vcf.gz
```

### 2.10 Ideogram

```bash
cd /share/dennislab/projects/t2t/variants/analysis/ideogram
```

**Cluster hets**

Bases impacted:
```bash
awk '{sum+=$3-$2}END{print sum}' hg38_cluster_hets.bed # 20821000
awk '{sum+=$3-$2}END{print sum}' hg38_cluster_hets.bed # 20821000
```

Genes impacted:
```bash
bedtools intersect -wo -a <(cut -d";" -f1 /share/dennislab/projects/t2t/variants/analysis/genes/gencode.v35.annotation.genesOnly.bed) -b hg38_cluster_hets.bed | cut -f4 | sort | uniq | wc -l # 987

bedtools intersect -wo -a <(cut -d";" -f1 /share/dennislab/projects/t2t/variants/analysis/genes/CHM13.combined.v4.genesOnly.bed) -b chm13_cluster_hets.bed | cut -f4 | sort | uniq | wc -l # 430
```

GWAS hits:
```bash
bedtools intersect -a /share/dennislab/users/dcsoto/ms_hsd/1_variability/GWAS/GWAScatalogue.2.bed -b hg38_cluster_hets.bed | wc -l # 324
```

**Missing copies**

Bases impacted:
```bash
awk '{sum+=$3-$2}END{print sum}' hg38_missing_dups.bed # 8041000
awk '{sum+=$3-$2}END{print sum}' chm13_missing_dups.bed # 121000
```

```bash
bedtools intersect -wo -a <(cut -d";" -f1 /share/dennislab/projects/t2t/variants/analysis/genes/gencode.v35.annotation.genesOnly.bed) -b hg38_missing_dups.bed | cut -f4 | sort | uniq | wc -l # 307

bedtools intersect -wo -a <(cut -d";" -f1 /share/dennislab/projects/t2t/variants/analysis/genes/CHM13.combined.v4.genesOnly.bed) -b chm13_missing_dups.bed | cut -f4 | sort | uniq | wc -l # 7
```

GWAS hits in GRCh38:
```bash
bedtools intersect -a /share/dennislab/users/dcsoto/ms_hsd/1_variability/GWAS/GWAScatalogue.2.bed -b hg38_missing_dups.bed | wc -l # 30
```

**Falsely duplicated**

Bases impacted:
```bash
awk '{sum+=$3-$2}END{print sum}' hg38_false_dups.bed # 1930833
```

```bash
bedtools intersect -wo -a <(cut -d";" -f1 /share/dennislab/projects/t2t/variants/analysis/genes/gencode.v35.annotation.genesOnly.bed) -b hg38_false_dups.bed | cut -f4 | sort | uniq | wc -l # 91
```

**LD-discordant**

Lifting over CHM13 ld-discordant coordinates to CHM13:
```bash
liftOver ld_discordant_haplotypes_CHM13.txt ../lifting/hg38.t2t-chm13-v1.0.over.chain ld_discordant_haplotypes_CHM13.to-CHM13.txt ld_discordant_haplotypes_CHM13.to-CHM13_unmapped.txt

liftOver ld_discordant_snp_pairs_CHM13.txt ../lifting/hg38.t2t-chm13-v1.0.over.chain ld_discordant_snp_pairs_CHM13.to-CHM13.txt ld_discordant_snp_pairs_CHM13.to-CHM13_unmapped.txt
```

```bash
wc -l ld_discordant_snp_pairs_GRCh38.bed # 18813
wc -l ld_discordant_snp_pairs_CHM13.txt # 209
```

```bash
bedtools intersect -F 1 -wo -a <(cut -d";" -f1 /share/dennislab/projects/t2t/variants/analysis/genes/gencode.v35.annotation.genesOnly.bed) -b ld_discordant_snp_pairs_GRCh38.bed | cut -f4 | sort | uniq | wc -l # 796

bedtools intersect -F 1 -wo -a <(cut -d";" -f1 /share/dennislab/projects/t2t/variants/analysis/genes/CHM13.combined.v4.genesOnly.bed) -b ld_discordant_snp_pairs_CHM13.to-CHM13.txt | cut -f4 | sort | uniq | wc -l # 34

bedtools intersect -F 1 -wo -a <(cut -d";" -f1 /share/dennislab/projects/t2t/variants/analysis/genes/gencode.v35.annotation.genesOnly.bed) -b ld_discordant_snp_pairs_CHM13.txt | cut -f4 | sort | uniq | wc -l # 34
```

**liftOver issues**

```bash
cut -f4 liftover_failures_dbsnp.all.bed | sort | uniq -c
cut -f4 liftover_failures_clinvar.all.bed | sort | uniq -c
cut -f4 liftover_failures_gwas.all.bed | sort | uniq -c

grep -v "MismatchedRefAllele" liftover_failures_dbsnp.all.bed | wc -l # 13061295
grep -v "MismatchedRefAllele" liftover_failures_clinvar.all.bed | wc -l # 1732
grep -v "MismatchedRefAllele" liftover_failures_gwas.all.bed | wc -l # 914

bedtools intersect -F 1 -wo -a <(cut -d";" -f1 /share/dennislab/projects/t2t/variants/analysis/genes/gencode.v35.annotation.genesOnly.bed) -b <(grep -v "MismatchedRefAllele" liftover_failures_dbsnp.all.bed) | cut -f4 | sort | uniq | wc -l # 27755

bedtools intersect -F 1 -wo -a <(cut -d";" -f1 /share/dennislab/projects/t2t/variants/analysis/genes/gencode.v35.annotation.genesOnly.bed) -b <(grep -v "MismatchedRefAllele" liftover_failures_clinvar.all.bed) | cut -f4 | sort | uniq | wc -l # 673

bedtools intersect -F 1 -wo -a <(cut -d";" -f1 /share/dennislab/projects/t2t/variants/analysis/genes/gencode.v35.annotation.genesOnly.bed) -b <(grep -v "MismatchedRefAllele" liftover_failures_gwas.all.bed) | cut -f4 | sort | uniq | wc -l # 606
```

- dbSNP: 736,178,420 (no Y chromosome or any patches/fixes for y chromosome)
- GWAS Catalog: 151,876 (subset from dbSNP, only primary chromosomes, no Y)
- Clinvar: 802,674 (no Y chromosome)

- 13061295/736178420*100 = 1.78%
- 1732/151876*100 = 1.14%
- 914/802674*100 = 0.11%

**Non-syntenic regions**

Bases impacted:
```bash
awk '{sum+=$3-$2}END{print sum}' GRCh38.no_snyteny_1Mbp.bed # 282173048
awk '{sum+=$3-$2}END{print sum}' chm13.draft_v1.0_plus38Y.no_snyteny_1Mbp.bed # 240044315
```

```bash
bedtools intersect -wo -a <(cut -d";" -f1 /share/dennislab/projects/t2t/variants/analysis/genes/gencode.v35.annotation.genesOnly.bed) -b GRCh38.no_snyteny_1Mbp.bed | cut -f4 | sort | uniq | wc -l # 2138

bedtools intersect -wo -a <(cut -d";" -f1 /share/dennislab/projects/t2t/variants/analysis/genes/CHM13.combined.v4.genesOnly.bed) -b chm13.draft_v1.0_plus38Y.no_snyteny_1Mbp.bed| cut -f4 | sort | uniq | wc -l # 1441
```
