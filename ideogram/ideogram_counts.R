library(tidyverse)
library(scales)
library(data.table)
library(glue)

setwd("/Users/dcsoto/Documents/PhD/Manuscripts/submitted/ms_T2T_variants/Plots/summaries")

mytheme <- theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + theme(legend.position = "none")

## FP het clusters

data.frame(bases=c(20821000, 67000),
           ref=factor(c("hg38","chm13"), level=c("hg38","chm13")),
           type=factor(c("cluster_hets","cluster_hets"),
                       levels=c("cluster_hets"))) %>%
  ggplot() +  geom_bar(aes(x=ref, y=bases, fill=type), stat="identity") + facet_wrap(~type) +
  scale_fill_manual(values=c("#a6cee3")) + scale_y_continuous(label=comma) + mytheme

ggsave("cluster_hets_bases.pdf", width=3, height=5)

data.frame(genes=c(994, 0),
           ref=factor(c("hg38","chm13"), level=c("hg38","chm13")),
           type=factor(c("cluster_hets","cluster_hets"),
                       levels=c("cluster_hets"))) %>%
  ggplot() +  geom_bar(aes(x=ref, y=genes, fill=type), stat="identity") + facet_wrap(~type) +
  scale_fill_manual(values=c("#a6cee3","#2F7DDE","#48D1CC")) + scale_y_continuous(label=comma) + mytheme

ggsave("cluster_hets_genes.pdf", width=3, height=5)

## Duplication errors

data.frame(bases=c(8041000, 7000, 1930833),
           ref=factor(c("hg38","chm13","hg38"), level=c("hg38","chm13")),
           type=factor(c("missing_dups", "missing_dups","false_dups"),
                       levels=c("missing_dups","false_dups"))) %>%
  ggplot() +  geom_bar(aes(x=ref, y=bases, fill=type), stat="identity") + facet_wrap(~type) +
  scale_fill_manual(values=c("#2F7DDE","#48D1CC")) + scale_y_continuous(label=comma, limits=c(0,10000000)) + mytheme

ggsave("copies_bases.pdf", width=3, height=5)

data.frame(genes=c(309, 0, 91),
           ref=factor(c("hg38","chm13","hg38"), level=c("hg38","chm13")),
           type=factor(c("missing_dups", "missing_dups","false_dups"),
                       levels=c("missing_dups","false_dups"))) %>%
  ggplot() +  geom_bar(aes(x=ref, y=genes, fill=type), stat="identity") + facet_wrap(~type) +
  scale_fill_manual(values=c("#2F7DDE","#48D1CC")) + scale_y_continuous(label=comma, limits=c(0,400)) + mytheme

ggsave("copies_genes.pdf", width=3, height=5)

# LD-discordant haps

data.frame(counts=c(18813, 209),ref=factor(c("hg38","chm13"), level=c("hg38","chm13")), type=c("LD-discordant SNPs")) %>%
  ggplot() +  geom_bar(aes(x=ref, y=counts, fill=type), stat="identity") + facet_wrap(~type) +
  scale_fill_manual(values=c("#7a52a1")) + scale_y_continuous(label=comma) + mytheme

ggsave("ld_discordant_snps.pdf", width=3, height=5)

data.frame(genes=c(797, 34),ref=factor(c("hg38","chm13"), level=c("hg38","chm13")), type=c("LD-discordant SNPs")) %>%
  ggplot() +  geom_bar(aes(x=ref, y=genes, fill=type), stat="identity") + facet_wrap(~type) +
  scale_fill_manual(values=c("#7a52a1")) + scale_y_continuous(label=comma) + mytheme

ggsave("ld_discordant_genes.pdf", width=3, height=5)

# Lift over issues plot

data.frame(counts=c(1.77, 0.21, 0.6),
           ref=factor(c("hg38","hg38","hg38"), level=c("hg38","chm13")), 
           type=factor(c("dbSNP","ClinVar","GWAS"), level=c("dbSNP","ClinVar","GWAS"))) %>%
  ggplot() +  geom_bar(aes(x=ref, y=counts, fill=type), stat="identity") + facet_wrap(~type) +
  scale_fill_manual(values=c("#f69999","#169099","#ffc961")) + scale_y_continuous(label=comma) + mytheme

ggsave("lifover_issues_counts.pdf", width=2, height=5)

data.frame(genes=c(27811, 674, 606),
           ref=factor(c("hg38","hg38","hg38"), level=c("hg38","chm13")), 
           type=factor(c("dbSNP","ClinVar","GWAS"), level=c("dbSNP","ClinVar","GWAS"))) %>%
  ggplot() +  geom_bar(aes(x=ref, y=genes, fill=type), stat="identity") + facet_wrap(~type) +
  scale_fill_manual(values=c("#f69999","#169099","#ffc961")) + scale_y_continuous(label=comma, trans='log10') + mytheme

ggsave("lifover_issues_genes.pdf", width=2, height=5)

## Non-syntenic

data.frame(bases=c(282173048, 240044315, 189036735),
           ref=factor(c("hg38","chm13","chm13"), level=c("hg38","chm13")), 
           type=c("non-syntenic","non-syntenic","novel")) %>%
  ggplot() +  geom_bar(aes(x=ref, y=bases, fill=type), stat="identity") + facet_wrap(~type) +
  scale_fill_manual(values=c("#aa773f","#aa773f")) + scale_y_continuous(label=comma, breaks = seq(0, 300000000, by = 50000000))  + mytheme

ggsave("non-syntenic_bases.pdf", width=3, height=5)

data.frame(genes=c(2170, 2905, 1791),
           ref=factor(c("hg38","chm13", "chm13"), level=c("hg38","chm13")), 
           type=c("non-syntenic", "non-syntenic", "novel")) %>%
  ggplot() +  geom_bar(aes(x=ref, y=genes, fill=type), stat="identity") + facet_wrap(~type) +
  scale_fill_manual(values=c("#aa773f","#aa773f")) + scale_y_continuous(label=comma)  + mytheme

ggsave("non-syntenic_genes.pdf", width=3, height=5)
