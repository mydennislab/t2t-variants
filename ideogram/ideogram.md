# Fig 1. Ideogram

```bash
cd /share/dennislab/projects/t2t/variants/analysis/ideogram
```

## 1. Cluster hets

Bases impacted:
```bash
awk '{sum+=$3-$2}END{print sum}' hg38_cluster_hets.bed # 20821000
awk '{sum+=$3-$2}END{print sum}' chm13_cluster_hets.bed # 67000
```

Genes impacted:
```bash
bedtools intersect -wo -a <(cut -d";" -f1 /share/dennislab/projects/t2t/variants/analysis/genes/gencode.v35.annotation.genesOnly.bed) -b hg38_cluster_hets.bed | cut -f4 | sort | uniq | wc -l # 994
bedtools intersect -wo -a <(cut -d";" -f1 /share/dennislab/projects/t2t/variants/analysis/genes/CHM13.combined.v4.genesOnly.bed) -b chm13_cluster_hets.bed | cut -f4 | sort | uniq | wc -l # 0
```

GWAS hits:
```bash
bedtools intersect -a /share/dennislab/users/dcsoto/ms_hsd/1_variability/GWAS/GWAScatalogue.2.bed -b hg38_cluster_hets.bed | wc -l # 324
```

## 2. Missing copies

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

## 3. Falsely duplicated

Bases impacted:
```bash
awk '{sum+=$3-$2}END{print sum}' hg38_false_dups.bed # 1930833
```

```bash
bedtools intersect -wo -a <(cut -d";" -f1 /share/dennislab/projects/t2t/variants/analysis/genes/gencode.v35.annotation.genesOnly.bed) -b hg38_false_dups.bed | cut -f4 | sort | uniq | wc -l # 91
```

## 4. LD-discordant

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

## 5. liftOver issues

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

## 6. Non-syntenic regions

Bases impacted:
```bash
awk '{sum+=$3-$2}END{print sum}' GRCh38.no_snyteny_1Mbp.bed # 282173048
awk '{sum+=$3-$2}END{print sum}' chm13.draft_v1.0_plus38Y.no_snyteny_1Mbp.bed # 240044315
```

```bash
bedtools intersect -wo -a <(cut -d";" -f1 /share/dennislab/projects/t2t/variants/analysis/genes/gencode.v35.annotation.genesOnly.bed) -b GRCh38.no_snyteny_1Mbp.bed | cut -f4 | sort | uniq | wc -l # 2138
bedtools intersect -wo -a <(cut -d";" -f1 /share/dennislab/projects/t2t/variants/analysis/genes/CHM13.combined.v4.genesOnly.bed) -b chm13.draft_v1.0_plus38Y.no_snyteny_1Mbp.bed| cut -f4 | sort | uniq | wc -l # 1441
```
