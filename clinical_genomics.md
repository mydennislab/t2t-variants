# Impact of T2T-CHM13 in clinical genomics

### 3.1 Additional analyses

#### 3.1.1 CHM13 DipCall

```bash
cd /share/dennislab/projects/t2t/variants/dipcall
```

> Manually added `##INFO=<ID=CSQ,Number=1,Type=String,Descriotion="VEP annotation">` to `chm13.202000921_with38Y-align2-GRCh38.dip.VEP.dipbed.vcf`

CHM13 DipCall to Hg38 variants:
```bash
bedtools intersect -header -a chm13.202000921_with38Y-align2-GRCh38.dip.VEP.vcf.gz -b chm13.202000921_with38Y-align2-GRCh38.dip.bed > chm13.202000921_with38Y-align2-GRCh38.dip.VEP.dipbed.vcf

bgzip chm13.202000921_with38Y-align2-GRCh38.dip.VEP.dipbed.vcf
tabix chm13.202000921_with38Y-align2-GRCh38.dip.VEP.dipbed.vcf.gz

bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\t%INFO/CSQ\n" chm13.202000921_with38Y-align2-GRCh38.dip.VEP.dipbed.vcf.gz | \
awk '{print $1"_"$2"_"$3"_"$4"\t"$5}' > chm13_dipcall_vep_toHg38.tsv
```

1kgp Hg38-called variants:
```bash
bcftools norm -m - /share/dennislab/projects/t2t/variants/analysis/med_genes/variants_1kgp/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_all.recalibrated_variants.pass.nogenos.sort.vcf.gz > 1kgp_variants_toHg38.norm.vcf

bgzip 1kgp_variants_toHg38.norm.vcf
tabix 1kgp_variants_toHg38.norm.vcf.gz

bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\t%INFO/AC\t%INFO/AF\n" 1kgp_variants_toHg38.norm.vcf.gz | \
awk '{print $1"_"$2"_"$3"_"$4"\t"$5"\t"$6}' > 1kgp_variants_toHg38.tsv
```

Annotating CHM13 dipcall to Hg38 variants with 1kgp AF:
```bash
join -1 1 -2 1 -t $'\t' <(sort -k1,1 chm13_dipcall_vep_toHg38.tsv) <(sort -k1,1 1kgp_variants_toHg38.tsv) > chm13_dipcall_vep_toHg38.1kgp_af.tsv
```

#### 3.1.2 1KGP AF=1

```bash
cd /share/dennislab/projects/t2t/variants/af_1/grch38_af_1

ls *vcf.gz | tr '\n' ' '

vcf-concat 20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr10.recalibrated_variants.pass.singleton.af_1.vcf.gz 20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr11.recalibrated_variants.pass.singleton.af_1.vcf.gz 20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr12.recalibrated_variants.pass.singleton.af_1.vcf.gz 20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr13.recalibrated_variants.pass.singleton.af_1.vcf.gz 20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr14.recalibrated_variants.pass.singleton.af_1.vcf.gz 20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr15.recalibrated_variants.pass.singleton.af_1.vcf.gz 20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr16.recalibrated_variants.pass.singleton.af_1.vcf.gz 20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr17.recalibrated_variants.pass.singleton.af_1.vcf.gz 20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr18.recalibrated_variants.pass.singleton.af_1.vcf.gz 20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr19.recalibrated_variants.pass.singleton.af_1.vcf.gz 20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr1.recalibrated_variants.pass.singleton.af_1.vcf.gz 20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr20.recalibrated_variants.pass.singleton.af_1.vcf.gz 20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr21.recalibrated_variants.pass.singleton.af_1.vcf.gz 20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr22.recalibrated_variants.pass.singleton.af_1.vcf.gz 20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr2.recalibrated_variants.pass.singleton.af_1.vcf.gz 20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr3.recalibrated_variants.pass.singleton.af_1.vcf.gz 20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr4.recalibrated_variants.pass.singleton.af_1.vcf.gz 20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr5.recalibrated_variants.pass.singleton.af_1.vcf.gz 20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr6.recalibrated_variants.pass.singleton.af_1.vcf.gz 20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr7.recalibrated_variants.pass.singleton.af_1.vcf.gz 20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr8.recalibrated_variants.pass.singleton.af_1.vcf.gz 20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr9.recalibrated_variants.pass.singleton.af_1.vcf.gz 20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chrX.recalibrated_variants.pass.singleton.af_1.vcf.gz > 20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_all.recalibrated_variants.pass.singleton.af_1.vcf

bgzip 20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_all.recalibrated_variants.pass.singleton.af_1.vcf
tabix 20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_all.recalibrated_variants.pass.singleton.af_1.vcf.gz
```

```bash
cd /share/dennislab/projects/t2t/variants/af_1/chm13_af_1

ls *vcf.gz | tr '\n' ' '
# Excluding chromosome Y
vcf-concat 1kgp.chr10.recalibrated.snp_indel.pass.singleton.af_1.vcf.gz 1kgp.chr11.recalibrated.snp_indel.pass.singleton.af_1.vcf.gz 1kgp.chr12.recalibrated.snp_indel.pass.singleton.af_1.vcf.gz 1kgp.chr13.recalibrated.snp_indel.pass.singleton.af_1.vcf.gz 1kgp.chr14.recalibrated.snp_indel.pass.singleton.af_1.vcf.gz 1kgp.chr15.recalibrated.snp_indel.pass.singleton.af_1.vcf.gz 1kgp.chr16.recalibrated.snp_indel.pass.singleton.af_1.vcf.gz 1kgp.chr17.recalibrated.snp_indel.pass.singleton.af_1.vcf.gz 1kgp.chr18.recalibrated.snp_indel.pass.singleton.af_1.vcf.gz 1kgp.chr19.recalibrated.snp_indel.pass.singleton.af_1.vcf.gz 1kgp.chr1.recalibrated.snp_indel.pass.singleton.af_1.vcf.gz 1kgp.chr20.recalibrated.snp_indel.pass.singleton.af_1.vcf.gz 1kgp.chr21.recalibrated.snp_indel.pass.singleton.af_1.vcf.gz 1kgp.chr22.recalibrated.snp_indel.pass.singleton.af_1.vcf.gz 1kgp.chr2.recalibrated.snp_indel.pass.singleton.af_1.vcf.gz 1kgp.chr3.recalibrated.snp_indel.pass.singleton.af_1.vcf.gz 1kgp.chr4.recalibrated.snp_indel.pass.singleton.af_1.vcf.gz 1kgp.chr5.recalibrated.snp_indel.pass.singleton.af_1.vcf.gz 1kgp.chr6.recalibrated.snp_indel.pass.singleton.af_1.vcf.gz 1kgp.chr7.recalibrated.snp_indel.pass.singleton.af_1.vcf.gz 1kgp.chr8.recalibrated.snp_indel.pass.singleton.af_1.vcf.gz 1kgp.chr9.recalibrated.snp_indel.pass.singleton.af_1.vcf.gz 1kgp.chrX.recalibrated.snp_indel.pass.singleton.af_1.vcf.gz > 1kgp.all.recalibrated.snp_indel.pass.singleton.af_1.vcf

bgzip 1kgp.all.recalibrated.snp_indel.pass.singleton.af_1.vcf
tabix 1kgp.all.recalibrated.snp_indel.pass.singleton.af_1.vcf.gz
```

Total variants:
```bash
cd /share/dennislab/projects/t2t/variants/af_1

gunzip -c grch38_af_1/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_all.recalibrated_variants.pass.singleton.af_1.vcf.gz | grep -v "^#" | wc -l # 20529 all
bcftools view --exclude-types indels grch38_af_1/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_all.recalibrated_variants.pass.singleton.af_1.vcf.gz | grep -v "^#" | wc -l # 13290 SNVs

gunzip -c chm13_af_1/1kgp.all.recalibrated.snp_indel.pass.singleton.af_1.vcf.gz | grep -v "^#" | wc -l # 10317 all
bcftools view --exclude-types indels chm13_af_1/1kgp.all.recalibrated.snp_indel.pass.singleton.af_1.vcf.gz | grep -v "^#" | wc -l # 9880 SNVs
```

### 3.2 Med genes impacted

```bash
cd /share/dennislab/projects/t2t/med_genes
```

1KGP SNV AF = 1:
```bash
bedtools intersect -wo -a GRCh38_mrg_full_gene.bed -b ../variants/af_1/grch38_af_1/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_all.recalibrated_variants.pass.singleton.af_1.vcf.gz > snvs_af_1/GRCh38_mrg_full_gene.annotated_af_1.bed

cut -f4 snvs_af_1/GRCh38_mrg_full_gene.annotated_af_1.bed | sort | uniq -c | tr -s " " | cut -d" " -f2- | awk '{print $2"\t"$1}' > snvs_af_1/snvs_af_1.hg38.txt

bedtools intersect -wo -a GRCh38_mrg_full_gene.t2t-chm13-v1.0.picard.sort.bed -b ../variants/af_1/chm13_af_1/1kgp.all.recalibrated.snp_indel.pass.singleton.af_1.vcf.gz > snvs_af_1/GRCh38_mrg_full_gene.t2t-chm13-v1.0.picard.sort.annotated_af_1.bed

cut -f4 snvs_af_1/GRCh38_mrg_full_gene.t2t-chm13-v1.0.picard.sort.annotated_af_1.bed | sort | uniq -c | tr -s " " | cut -d" " -f2- | awk '{print $2"\t"$1}' > snvs_af_1/snvs_af_1.chm13.txt
```

Cluster hets:
```bash
bedtools intersect -wo -a GRCh38_mrg_full_gene.bed -b /share/dennislab/projects/t2t/variants/analysis/ideogram/hg38_cluster_hets.bed > cluster_hets/GRCh38_mrg_full_gene.annotated_cluster_hets.bed

cut -f4 cluster_hets/GRCh38_mrg_full_gene.annotated_cluster_hets.bed | sort | uniq | awk '{print $1"\tYes"}' > cluster_hets/cluster_hets.hg38.txt

bedtools intersect -wo -a GRCh38_mrg_full_gene.t2t-chm13-v1.0.picard.sort.bed -b /share/dennislab/projects/t2t/variants/analysis/ideogram/chm13_cluster_hets.bed > cluster_hets/GRCh38_mrg_full_gene.t2t-chm13-v1.0.picard.sort.annotated_cluster_hets.bed

cut -f4 cluster_hets/GRCh38_mrg_full_gene.t2t-chm13-v1.0.picard.sort.annotated_cluster_hets.bed | sort | uniq | awk '{print $1"\tYes"}' > cluster_hets/cluster_hets.chm13.txt
```

Missing dups:
```bash
bedtools intersect -wo -a GRCh38_mrg_full_gene.bed -b /share/dennislab/projects/t2t/variants/analysis/ideogram/hg38_missing_dups.bed > missing_dups/GRCh38_mrg_full_gene.annotated_missing_dups.bed

cut -f4 missing_dups/GRCh38_mrg_full_gene.annotated_missing_dups.bed | sort | uniq | awk '{print $1"\tYes"}' > missing_dups/missing_dups.hg38.txt

bedtools intersect -wo -a GRCh38_mrg_full_gene.t2t-chm13-v1.0.picard.sort.bed -b /share/dennislab/projects/t2t/variants/analysis/ideogram/chm13_missing_dups.bed > missing_dups/GRCh38_mrg_full_gene.t2t-chm13-v1.0.picard.sort.annotated_missing_dups.bed

cut -f4 missing_dups/GRCh38_mrg_full_gene.t2t-chm13-v1.0.picard.sort.annotated_missing_dups.bed | sort | uniq | awk '{print $1"\tYes"}' > missing_dups/missing_dups.chm13.txt
```

False dups:
```bash
bedtools intersect -wo -a GRCh38_mrg_full_gene.bed -b /share/dennislab/projects/t2t/variants/analysis/ideogram/hg38_false_dups.bed > false_dups/GRCh38_mrg_full_gene.annotated_false_dups.bed

cut -f4 false_dups/GRCh38_mrg_full_gene.annotated_false_dups.bed | sort | uniq | awk '{print $1"\tYes"}' > false_dups/false_dups.hg38.txt
```

LD-discordant SNPs pairs:
```bash
bedtools intersect -F 1 -wo -a GRCh38_mrg_full_gene.bed -b /share/dennislab/projects/t2t/variants/analysis/ideogram/ld_discordant_snp_pairs_GRCh38.bed > ld_discordant/GRCh38_mrg_full_gene.annotated_ld_discordant.bed

cut -f4 ld_discordant/GRCh38_mrg_full_gene.annotated_ld_discordant.bed | sort | uniq | awk '{print $1"\tYes"}' > ld_discordant/ld_discordant.hg38.txt

bedtools intersect -F 1 -wo -a GRCh38_mrg_full_gene.t2t-chm13-v1.0.picard.sort.bed -b /share/dennislab/projects/t2t/variants/analysis/ideogram/ld_discordant_snp_pairs_CHM13.to-CHM13.txt > ld_discordant/GRCh38_mrg_full_gene.t2t-chm13-v1.0.picard.sort.annotated_ld_discordant.bed

cut -f4 ld_discordant/GRCh38_mrg_full_gene.t2t-chm13-v1.0.picard.sort.annotated_ld_discordant.bed | sort | uniq | awk '{print $1"\tYes"}' > ld_discordant/ld_discordant.chm13.txt
```

Non-syntenic:
```bash
bedtools intersect -wo -a GRCh38_mrg_full_gene.bed -b /share/dennislab/projects/t2t/variants/analysis/synteny/GRCh38.no_snyteny_1Mbp.bed > nonsyntenic/GRCh38_mrg_full_gene.annotated_nonsyntenic.bed

cut -f4 nonsyntenic/GRCh38_mrg_full_gene.annotated_nonsyntenic.bed | sort | uniq | awk '{print $1"\tYes"}' > nonsyntenic/nonsyntenic.hg38.txt

bedtools intersect -wo -a GRCh38_mrg_full_gene.t2t-chm13-v1.0.picard.sort.bed -b /share/dennislab/projects/t2t/variants/analysis/synteny/chm13.draft_v1.0_plus38Y.no_snyteny_1Mbp.bed > nonsyntenic/GRCh38_mrg_full_gene.t2t-chm13-v1.0.picard.sort.annotated_nonsyntenic.bed

cut -f4 nonsyntenic/GRCh38_mrg_full_gene.t2t-chm13-v1.0.picard.sort.annotated_nonsyntenic.bed | sort | uniq | awk '{print $1"\tYes"}' > nonsyntenic/nonsyntenic.chm13.txt
```

LiftOver failures:
```bash
# dbSNP
bedtools intersect -wo -a GRCh38_mrg_full_gene.bed -b /share/dennislab/projects/t2t/variants/analysis/ideogram/liftover_failures_dbsnp.all.bed > liftover_failures/GRCh38_mrg_full_gene.annotated_liftover_dbsnp.bed

cut -f4 liftover_failures/GRCh38_mrg_full_gene.annotated_liftover_dbsnp.bed | grep -v "MismatchedRefAllele" | sort | uniq -c | tr -s " " | cut -d" " -f2- | awk '{print $2"\t"$1}' > liftover_failures/liftover_failures_dbsnp.hg38.txt

# ClinVar
bedtools intersect -wo -a GRCh38_mrg_full_gene.bed -b /share/dennislab/projects/t2t/variants/analysis/ideogram/liftover_failures_clinvar.all.bed > liftover_failures/GRCh38_mrg_full_gene.annotated_liftover_clinvar.bed

cut -f4 liftover_failures/GRCh38_mrg_full_gene.annotated_liftover_clinvar.bed | grep -v "MismatchedRefAllele" | sort | uniq -c | tr -s " " | cut -d" " -f2- | awk '{print $2"\t"$1}' > liftover_failures/liftover_failures_clinvar.hg38.txt

# GWAS
bedtools intersect -wo -a GRCh38_mrg_full_gene.bed -b /share/dennislab/projects/t2t/variants/analysis/ideogram/liftover_failures_gwas.all.bed > liftover_failures/GRCh38_mrg_full_gene.annotated_liftover_gwas.bed

cut -f4 liftover_failures/GRCh38_mrg_full_gene.annotated_liftover_gwas.bed | grep -v "MismatchedRefAllele" | sort | uniq -c | tr -s " " | cut -d" " -f2- | awk '{print $2"\t"$1}' > liftover_failures/liftover_failures_gwas.hg38.txt
```

GRC issues:
```bash
bedtools intersect -wo -a GRCh38_mrg_full_gene.bed -b /share/dennislab/projects/t2t/variants/analysis/ideogram/hg38.parsedissues.bed > grc_issues/GRCh38_mrg_full_gene.annotated_grc_issues.bed

cut -f4,10,11 grc_issues/GRCh38_mrg_full_gene.annotated_grc_issues.bed | awk '{print $1"\t"$2";"$3}' | sort | uniq > grc_issues/grc_issues.hg38.txt
```

### 3.3 KCNJ18 and KCNE1 close-up

```bash
cd /share/dennislab/projects/t2t/variants/analysis/kcnj18/
```

#### 3.4.1 1KGP variant files

```bash
cd /share/dennislab/projects/t2t/variants/analysis/kcnj18/variants_1kgp
```

Downloading CHM13-called 1KGP original files (labeled) and concatenating:
```bash
for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX; do
  globus transfer 9db1f0a6-a05a-11ea-8f06-0a21f750d19b:/team-variants/grch38_t2t_liftover/20200921/variant_set_round_trip_tests/1kg_links/1kgp.$chr.recalibrated.snp_indel.pass.nogenos.labeled.vcf.gz $(globus endpoint local-id):/share/dennislab/users/dcsoto/other/globus_connect/1kgp.$chr.recalibrated.snp_indel.pass.nogenos.labeled.vcf.gz
done

vcf-concat 1kgp.chr10.recalibrated.snp_indel.pass.nogenos.labeled.vcf.gz 1kgp.chr11.recalibrated.snp_indel.pass.nogenos.labeled.vcf.gz 1kgp.chr12.recalibrated.snp_indel.pass.nogenos.labeled.vcf.gz 1kgp.chr13.recalibrated.snp_indel.pass.nogenos.labeled.vcf.gz 1kgp.chr14.recalibrated.snp_indel.pass.nogenos.labeled.vcf.gz 1kgp.chr15.recalibrated.snp_indel.pass.nogenos.labeled.vcf.gz 1kgp.chr16.recalibrated.snp_indel.pass.nogenos.labeled.vcf.gz 1kgp.chr17.recalibrated.snp_indel.pass.nogenos.labeled.vcf.gz 1kgp.chr18.recalibrated.snp_indel.pass.nogenos.labeled.vcf.gz 1kgp.chr19.recalibrated.snp_indel.pass.nogenos.labeled.vcf.gz 1kgp.chr1.recalibrated.snp_indel.pass.nogenos.labeled.vcf.gz 1kgp.chr20.recalibrated.snp_indel.pass.nogenos.labeled.vcf.gz 1kgp.chr21.recalibrated.snp_indel.pass.nogenos.labeled.vcf.gz 1kgp.chr22.recalibrated.snp_indel.pass.nogenos.labeled.vcf.gz 1kgp.chr2.recalibrated.snp_indel.pass.nogenos.labeled.vcf.gz 1kgp.chr3.recalibrated.snp_indel.pass.nogenos.labeled.vcf.gz 1kgp.chr4.recalibrated.snp_indel.pass.nogenos.labeled.vcf.gz 1kgp.chr5.recalibrated.snp_indel.pass.nogenos.labeled.vcf.gz 1kgp.chr6.recalibrated.snp_indel.pass.nogenos.labeled.vcf.gz 1kgp.chr7.recalibrated.snp_indel.pass.nogenos.labeled.vcf.gz 1kgp.chr8.recalibrated.snp_indel.pass.nogenos.labeled.vcf.gz 1kgp.chr9.recalibrated.snp_indel.pass.nogenos.labeled.vcf.gz 1kgp.chrX.recalibrated.snp_indel.pass.nogenos.labeled.vcf.gz > 1kgp.all.recalibrated.snp_indel.pass.nogenos.labeled.vcf
```

- GRCh38-called variants: `/share/dennislab/projects/t2t/variants/analysis/kcnj18/variants_1kgp/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_all.recalibrated_variants.pass.nogenos.sort.vcf.gz`
- CHM13-called variants: `/share/dennislab/projects/t2t/variants/analysis/kcnj18/variants_1kgp/1kgp.all.recalibrated.snp_indel.pass.nogenos.labeled.sort.vcf.gz`
- Lifted variants: `/share/dennislab/projects/t2t/variants/analysis/kcnj18/variants_1kgp/1kgp.allvars.recalibrated.snp_indel.pass.nogenos.labeled.grch38.sort.vcf.gz`

#### 3.4.2 KCNJ18/17

```bash
cd /share/dennislab/projects/t2t/variants/analysis/kcnj18/kcnj18_17
```

Coding region KCNJ18:
- CHM13	chr17:21651942-21653243
- Hg38	chr17:21702787-21704088

Coordinates of new copy in CHM13 (KCNJ17):
- KNCJ18-new: chr17:22634421-22683415

First, we extracted only KCNJ18 variants in CHM13 (chr17:21651942-21653243):
```bash
bcftools view --exclude-types indels --regions chr17:21651942-21653243 ../variants_1kgp/1kgp.all.recalibrated.snp_indel.pass.nogenos.labeled.sort.vcf.gz | bcftools norm -m - > 1kgp.chm13_all_snps.kcnj18.vcf
```

- Lines total/split/realigned/skipped: 112/5/0/0

Extracting information about variants and allele frequencies:
```bash
bcftools query -f "%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/AF\t%INFO/AC\n" 1kgp.chm13_all_snps.kcnj18.vcf | awk '{print $1"_"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' > kcnj18_af.chm13.tsv
```

Then we searched for KCNJ18 CHM13 variants that were successfully lifted to GRCh38. We did it by looking at the original CHM13 variants IDs in the lifted to GRCh38 variants VCF file:
```bash
# KCNJ18 CHM13-called variants
cat 1kgp.chm13_all_snps.kcnj18.vcf | grep -v "^#" | cut -f3 > kcnj18_chm13_ids.list

bcftools view --exclude-types indels --include ID==@kcnj18_chm13_ids.list ../variants_1kgp/1kgp.allvars.recalibrated.snp_indel.pass.nogenos.labeled.grch38.sort.vcf.gz | bcftools norm -m - > 1kgp.chm13_all_snps_toHg38.kcnj18.vcf
```

- Lines total/split/realigned/skipped: 112/5/0/0

Extracting information of CHM13-called variants lifted to Hg38:
```bash
bcftools query -f "%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/AF\t%INFO/AC\n" 1kgp.chm13_all_snps_toHg38.kcnj18.vcf | awk '{print $1"_"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' > kcnj18_af.chm13_toHg38.tsv
```

Finally, we extracted only KCNJ18 variants in Hg38-called VCF (chr17:21702787-21704088):
```bash
bcftools view --exclude-types indels --regions chr17:21702787-21704088 ../variants_1kgp/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_all.recalibrated_variants.pass.nogenos.sort.vcf.gz | bcftools norm -m - > 1kgp.hg38_all_snps.kcnj18.vcf
```

- Lines total/split/realigned/skipped: 129/8/0/0

Extracting information of Hg38-called variants:
```bash
bcftools query -f "%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/AF\t%INFO/AC\n" 1kgp.hg38_all_snps.kcnj18.vcf | awk '{print $1"_"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' > kcnj18_af.hg38.tsv
```

Additionally, we extracted only KCNJ17 CHM13-called variants (chr17:22634421-22683415):
```bash
bcftools view --exclude-types indels --regions chr17:22634421-22683415 ../variants_1kgp/1kgp.all.recalibrated.snp_indel.pass.nogenos.labeled.sort.vcf.gz | bcftools norm -m - > 1kgp.chm13_all_snps.kcnj17.vcf
```

- Lines total/split/realigned/skipped: 2704/112/0/0

Extracting information about variants and allele frequencies:
```bash
bcftools query -f "%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/AF\t%INFO/AC\n" 1kgp.chm13_all_snps.kcnj17.vcf | awk '{print $1"_"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' > kcnj17_af.chm13.tsv
```

Then we look for their lifted positions in GRCh38:
```bash
# KCNJ18 CHM13-called variants
cat 1kgp.chm13_all_snps.kcnj17.vcf | grep -v "^#" | cut -f3 > kcnj17_chm13_ids.list

bcftools view --exclude-types indels --include ID==@kcnj17_chm13_ids.list ../variants_1kgp/1kgp.allvars.recalibrated.snp_indel.pass.nogenos.labeled.grch38.sort.vcf.gz | bcftools norm -m - > 1kgp.chm13_all_snps_toHg38.kcnj17.vcf
```

- Lines total/split/realigned/skipped: 2629/81/0/0

Extracting information of Hg38-called variants:
```bash
bcftools query -f "%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/AF\t%INFO/AC\n" 1kgp.chm13_all_snps_toHg38.kcnj17.vcf | awk '{print $1"_"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' > kcnj17_af.chm13_toHg38.tsv
```

```r
# module load R/4.0.1
# R
library(tidyverse)

kcnj18_chm13  <- read.table("kcnj18_af.chm13.tsv",
  col.names=c("pos_chm13_kcnj18","id_chm13_kcnj18","ref_chm13_kcnj18","alt_chm13_kcnj18","af_chm13_kcnj18","ac_chm13_kcnj18"))

kcnj18_lifted <- read.table("kcnj18_af.chm13_toHg38.tsv",
  col.names=c("pos_hg38_kcnj18","id_chm13_kcnj18","V3","V4","V5","V6")) %>% select(pos_hg38_kcnj18, id_chm13_kcnj18)

kcnj18_hg38 <- read.table("kcnj18_af.hg38.tsv",
  col.names=c("pos_hg38_kcnj18","id_hg38_kcnj18","ref_hg38_kcnj18","alt_hg38_kcnj18","af_hg38_kcnj18","ac_hg38_kcnj18"))

df1 <- list(kcnj18_chm13, kcnj18_lifted) %>% reduce(full_join, by = "id_chm13_kcnj18") %>% distinct()
df2 <- list(df1, kcnj18_hg38) %>% reduce(full_join, by = "pos_hg38_kcnj18") %>% distinct()

kcnj17_chm13 <- read.table("kcnj17_af.chm13.tsv",
  col.names=c("pos_chm13_kcnj17","id_chm13_kcnj17","ref_chm13_kcnj17","alt_chm13_kcnj17","af_chm13_kcnj17","ac_chm13_kcnj17"))

kcnj17_lifted <- read.table("kcnj17_af.chm13_toHg38.tsv",
  col.names=c("pos_hg38_kcnj17","id_chm13_kcnj17","V3","V4","V5","V6")) %>% select(pos_hg38_kcnj17, id_chm13_kcnj17)

df3 <- list(kcnj17_chm13, kcnj17_lifted) %>% reduce(full_join, by = "id_chm13_kcnj17") %>% distinct()
df4 <- full_join(df2, df3, by = c("pos_hg38_kcnj18" = "pos_hg38_kcnj17"))

write.table(df4, "kcnj18_17.tsv", quote=FALSE, row.names=FALSE, col.names=TRUE)
```

#### 3.4.3 KCNJ12/17

KCNJ12:
- KCNJ12_Hg38: chr17:21415343-21416644
- KCNJ12_CHM13: chr17:21364558-21365859

KCNJ17 (KCNJ18-new):
- KCNJ17_CHM13: chr17:22666582-22667880

KCNJ12 Hg38-called VCF (chr17:21415343-21416644):
```bash
bcftools view --exclude-types indels --regions chr17:21415343-21416644 ../variants_1kgp/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_all.recalibrated_variants.pass.nogenos.sort.vcf.gz | bcftools norm -m - > 1kgp.hg38_all_snps.kcnj12.vcf

bcftools query -f "%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/AF\t%INFO/AC\n" 1kgp.hg38_all_snps.kcnj12.vcf | awk '{print $1"_"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' > kcnj12_af.hg38.tsv
```

> Lines total/split/realigned/skipped: 61/2/0/0

KCNJ12 CHM13-called VCF (chr17:21364558-21365859):
```bash
bcftools view --exclude-types indels --regions chr17:21364558-21365859 ../variants_1kgp/1kgp.all.recalibrated.snp_indel.pass.nogenos.labeled.sort.vcf.gz | bcftools norm -m - > 1kgp.chm13_all_snps.kcnj12.vcf

bcftools query -f "%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/AF\t%INFO/AC\n" 1kgp.chm13_all_snps.kcnj12.vcf | awk '{print $1"_"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' > kcnj12_af.chm13.tsv
```

> Lines total/split/realigned/skipped: 45/0/0/0

Lifted KCNJ12 CHM13-called variants:
```bash
cat 1kgp.chm13_all_snps.kcnj12.vcf | grep -v "^#" | cut -f3 > kcnj12_chm13_ids.list

bcftools view --exclude-types indels --include ID==@kcnj12_chm13_ids.list ../variants_1kgp/1kgp.allvars.recalibrated.snp_indel.pass.nogenos.labeled.grch38.sort.vcf.gz | bcftools norm -m - > 1kgp.chm13_all_snps_toHg38.kcnj12.vcf

bcftools query -f "%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/AF\t%INFO/AC\n" 1kgp.chm13_all_snps_toHg38.kcnj12.vcf | awk '{print $1"_"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' > kcnj12_af.chm13_toHg38.tsv
```

> Lines total/split/realigned/skipped: 45/0/0/0

KCNJ17 CDS CHM13-called VCF (chr17:22666582-22667880):
```bash
bcftools view --exclude-types indels --regions chr17:22666582-22667880 ../variants_1kgp/1kgp.all.recalibrated.snp_indel.pass.nogenos.labeled.sort.vcf.gz | bcftools norm -m - > 1kgp.chm13_all_snps.kcnj17.vcf

bcftools query -f "%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/AF\t%INFO/AC\n" 1kgp.chm13_all_snps.kcnj17.vcf | awk '{print $1"_"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' > kcnj17_af.chm13.tsv
```

> Lines total/split/realigned/skipped: 150/4/0/0

Lifted KCNJ17 CHM13-called variants:
```bash
cat 1kgp.chm13_all_snps.kcnj17.vcf | grep -v "^#" | cut -f3 > kcnj17_chm13_ids.list

bcftools view --exclude-types indels --include ID==@kcnj17_chm13_ids.list ../variants_1kgp/1kgp.allvars.recalibrated.snp_indel.pass.nogenos.labeled.grch38.sort.vcf.gz | bcftools norm -m - > 1kgp.chm13_all_snps_toHg38.kcnj17.vcf

bcftools query -f "%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/AF\t%INFO/AC\n" 1kgp.chm13_biallelic_snps_toHg38.kcnj17.vcf | awk '{print $1"_"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' > kcnj17_af.chm13_toHg38.tsv
```

```r
# module load R/4.0.1
# R
library(tidyverse)

kcnj12_hg38 <- read.table("kcnj12_af.hg38.tsv",
  col.names=c("pos_hg38_kcnj12","id_hg38_kcnj12","ref_hg38_kcnj12","alt_hg38_kcnj12","af_hg38_kcnj12","ac_hg38_kcnj12"))
kcnj12_chm13  <- read.table("kcnj12_af.chm13.tsv",
  col.names=c("pos_chm13_kcnj12","id_chm13_kcnj12","ref_chm13_kcnj12","alt_chm13_kcnj12","af_chm13_kcnj12","ac_chm13_kcnj12"))
kcnj12_lifted <- read.table("kcnj12_af.chm13_toHg38.tsv",
  col.names=c("pos_hg38_kcnj12","id_chm13_kcnj12","V3","V4","V5","V6")) %>% select(pos_hg38_kcnj12, id_chm13_kcnj12)

df1 <- list(kcnj12_hg38, kcnj12_lifted) %>% reduce(full_join, by = "pos_hg38_kcnj12")
df2 <- list(df1, kcnj12_chm13) %>% reduce(full_join, by = "id_chm13_kcnj12")

kcnj17_chm13 <- read.table("kcnj17_af.chm13.tsv",
  col.names=c("pos_chm13_kcnj17","id_chm13_kcnj17","ref_chm13_kcnj17","alt_chm13_kcnj17","af_chm13_kcnj17","ac_chm13_kcnj17"))
kcnj17_lifted <- read.table("kcnj17_af.chm13_toHg38.tsv",
  col.names=c("pos_hg38_kcnj17","id_chm13_kcnj17","V3","V4","V5","V6")) %>% select(pos_hg38_kcnj17, id_chm13_kcnj17)

df3 <- list(kcnj17_chm13, kcnj17_lifted) %>% reduce(full_join, by = "id_chm13_kcnj17")
df4 <- full_join(df2, df3, by = c("pos_hg38_kcnj12" = "pos_hg38_kcnj17"))

write.table(df4, "kcnj12_17.tsv", quote=FALSE, row.names=FALSE, col.names=TRUE)
```

#### 3.4.4 KCNJ18/17/12 SFS

```bash
cd /share/dennislab/projects/t2t/variants/analysis/kcnj18/sfs
```

Coding region KCNJ18:
- Hg38	chr17:21702787-21704088
- CHM13	chr17:21651942-21653243

KCNJ12:
- KCNJ12_Hg38: chr17:21415343-21416644
- KCNJ12_CHM13: chr17:21364558-21365859

KCNJ17 (KCNJ18-new):
- KCNJ17_CHM13: chr17:22666582-22667880

```bash
bcftools view --min-ac=1 --max-alleles 2 --exclude-types indels --regions chr17:21702787-21704088 ../variants_1kgp/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_all.recalibrated_variants.pass.nogenos.sort.vcf.gz | bcftools query -f "%INFO/AF\n" > KCNJ18_Hg38_AF.txt

bcftools view --min-ac=1 --max-alleles 2 --exclude-types indels --regions chr17:21415343-21416644 ../variants_1kgp/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_all.recalibrated_variants.pass.nogenos.sort.vcf.gz | bcftools query -f "%INFO/AF\n" > KCNJ12_Hg38_AF.txt
```

```bash
bcftools view --min-ac=1 --max-alleles 2 --exclude-types indels --regions chr17:21651942-21653243 ../variants_1kgp/1kgp.all.recalibrated.snp_indel.pass.nogenos.labeled.sort.vcf.gz | bcftools query -f "%INFO/AF\n" > KCNJ18_CHM13_AF.txt

bcftools view --min-ac=1 --max-alleles 2 --exclude-types indels --regions chr17:21364558-21365859 ../variants_1kgp/1kgp.all.recalibrated.snp_indel.pass.nogenos.labeled.sort.vcf.gz | bcftools query -f "%INFO/AF\n" > KCNJ12_CHM13_AF.txt

bcftools view --min-ac=1 --max-alleles 2 --exclude-types indels --regions chr17:22666582-22667880 ../variants_1kgp/1kgp.all.recalibrated.snp_indel.pass.nogenos.labeled.sort.vcf.gz | bcftools query -f "%INFO/AF\n" > KCNJ17_CHM13_AF.txt
```

```r
# module load R/4.0.1
# R
library(tidyverse)

# Read data
kcnj18_hg38 <- read.table("KCNJ18_Hg38_AF.txt", col.names=c("AF")) %>% mutate(gene="KCNJ18", ref="Hg38")
kcnj12_hg38 <- read.table("KCNJ12_Hg38_AF.txt", col.names=c("AF")) %>% mutate(gene="KCNJ12", ref="Hg38")

kcnj18_chm13 <- read.table("KCNJ18_CHM13_AF.txt", col.names=c("AF")) %>% mutate(gene="KCNJ18", ref="CHM13")
kcnj12_chm13 <- read.table("KCNJ12_CHM13_AF.txt", col.names=c("AF")) %>% mutate(gene="KCNJ12", ref="CHM13")
kcnj17_chm13 <- read.table("KCNJ17_CHM13_AF.txt", col.names=c("AF")) %>% mutate(gene="KCNJ17", ref="CHM13")

data <- rbind(kcnj18_hg38, kcnj12_hg38, kcnj18_chm13, kcnj12_chm13, kcnj17_chm13)

# SFS AF
ggplot(data, aes(AF, color=ref, fill=ref)) +
  geom_histogram(position=position_dodge(), binwidth=0.05) +
  xlab("SFS") +
  theme_bw() + facet_wrap(~gene)
ggsave("KCNJ18_12_17.SFS.pdf", width=9, height=3)
```

#### 3.4.4 HG002 variants

**Illumina variants**

```bash
cd /share/dennislab/projects/t2t/variants/analysis/kcnj18/variants_illm
```

```bash
bcftools view -f .,PASS --exclude-types indels -s HG002 --min-ac 1 GIABTrios.grch38.recalibrated.snp_indel.vcf.gz > HG002.grch38.recalibrated.pass.all_snps.vcf
bgzip HG002.grch38.recalibrated.pass.all_snps.vcf
tabix HG002.grch38.recalibrated.pass.all_snps.vcf.gz

bcftools view -f .,PASS --exclude-types indels -s HG002 --min-ac 1 GIABTrios.chm13.recalibrated.snp_indel.vcf.gz > HG002.chm13.recalibrated.pass.all_snps.vcf
bgzip HG002.chm13.recalibrated.pass.all_snps.vcf
tabix HG002.chm13.recalibrated.pass.all_snps.vcf.gz
```

**PBCCS variants**

```bash
cd /share/dennislab/projects/t2t/variants/analysis/kcnj18/variants_lrs
```

```bash
bcftools view -f .,PASS --exclude-types indels HG002vGRCh38_DeepVariantv1_PBCCS.vcf.gz > HG002vGRCh38_DeepVariantv1_PBCCS.pass.all_snps.vcf
bgzip HG002vGRCh38_DeepVariantv1_PBCCS.pass.all_snps.vcf
tabix HG002vGRCh38_DeepVariantv1_PBCCS.pass.all_snps.vcf.gz

bcftools view -f .,PASS --exclude-types indels HG002vCHM13_T2Tv1_DeepVariantv1_PBCCS.vcf.gz > HG002vCHM13_T2Tv1_DeepVariantv1_PBCCS.pass.all_snps.vcf
bgzip HG002vCHM13_T2Tv1_DeepVariantv1_PBCCS.pass.all_snps.vcf
tabix HG002vCHM13_T2Tv1_DeepVariantv1_PBCCS.pass.all_snps.vcf.gz
```

**Benchmark variants**

```bash
cd /share/dennislab/projects/t2t/variants/analysis/benchmark
```

```bash
bcftools view --exclude-types indels HG002_CHM13v1.0_CMRG_smallvar_v1.00_draft.vcf.gz > HG002_CHM13v1.0_CMRG_smallvar_v1.00_draft.snp.vcf
bgzip HG002_CHM13v1.0_CMRG_smallvcd ar_v1.00_draft.snp.vcf
tabix HG002_CHM13v1.0_CMRG_smallvar_v1.00_draft.snp.vcf.gz

bcftools view --exclude-types indels HG002_GRCh38_CMRG_smallvar_v1.00.vcf.gz > HG002_GRCh38_CMRG_smallvar_v1.00.snp.vcf
bgzip HG002_GRCh38_CMRG_smallvar_v1.00.snp.vcf
tabix HG002_GRCh38_CMRG_smallvar_v1.00.snp.vcf.gz

bcftools view --exclude-types indels HG002v11-align2-GRCh38.dip.vcf.gz > HG002v11-align2-GRCh38.dip.snp.vcf
bgzip HG002v11-align2-GRCh38.dip.snp.vcf
tabix HG002v11-align2-GRCh38.dip.snp.vcf.gz
```
