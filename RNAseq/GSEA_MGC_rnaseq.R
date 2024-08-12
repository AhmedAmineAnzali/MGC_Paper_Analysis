library(tidyverse)
#library(tximportData)
#library(tximport)
library(readr)
library(rhdf5)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(biomaRt)
library(DESeq2)
library(data.table)
library(clusterProfiler)
BiocManager::install("TCGAbiolinks")
library(TCGAbiolinks)
BiocManager::install("zellkonverter")
BiocManager::install("AnnotationDbi")
BiocManager::install("EDASeq")
library(EDASeq)
#BiocManager::install("ensembldb")
library(ensembldb)
#BiocManager::install("EnsDb.Hsapiens.v86")
library(EnsDb.Hsapiens.v86)
#BiocManager::install("DESeq2")
#BiocManager::install("clusterProfiler")
#BiocManager::install("enrichplot", force = T)
library(enrichplot)

Sig_MGC = read.csv("/Users/kmulder/data/nanostring/MGC_genes/Visium_MGC_Signature.csv", header = 1) %>% 
  dplyr::arrange(p_val_adj) %>% 
  head(100) %>% pull(gene)

# Query platform for HNSC bulk RNA sequencing data
query <- TCGAbiolinks::GDCquery(
  project = "TCGA-HNSC", 
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
)


# Download the data for HNSC patients
TCGAbiolinks::GDCdownload(query,files.per.chunk = 10)

# Prepare expression matrix with geneID in the rows and samples (barcode) in the columns
# rsem.genes.results as values
HNSC.Rnaseq.SE <- TCGAbiolinks::GDCprepare(query)

# For gene expression if you need to see a boxplot correlation and AAIC plot to define outliers you can run
HNSC.RNAseq_CorOutliers <- TCGAbiolinks::TCGAanalyze_Preprocessing(HNSC.Rnaseq.SE)

# normalization of genes
dataNorm <- TCGAbiolinks::TCGAanalyze_Normalization(
  tabDF = HNSC.RNAseq_CorOutliers, 
  geneInfo =  TCGAbiolinks::geneInfoHT
)

# quantile filter of genes
dataFilt <- TCGAbiolinks::TCGAanalyze_Filtering(
  tabDF = dataNorm,
  method = "quantile", 
  qnt.cut =  0.25
)

# selection of normal samples "NT"
samplesNT <- TCGAbiolinks::TCGAquery_SampleTypes(
  barcode = colnames(dataFilt),
  typesample = c("NT")
)

# selection of tumor samples "TP"
samplesTP <- TCGAbiolinks::TCGAquery_SampleTypes(
  barcode = colnames(dataFilt), 
  typesample = c("TP")
)

datasampleNT <- dataFilt[,samplesNT]
datasampleTP <- dataFilt[,samplesTP]

# write into csv
write.csv(datasampleTP, 'Tumor_HNSC_bulk.csv',row.names = T)
write.csv(datasampleNT, 'RNAseq/Normal_HNSC_bulk.csv',row.names = T)



#Reading clinical data
clinical <- read.csv('Clinical/TCGA_HNSC_Surival.csv', row.names = 1, sep = ';')

#Adding clinical data descriptions
#Keratin state
clinical <- clinical %>% mutate(K_state = case_when(
  as.numeric(as.character(clinical$K)) >= 10 ~ "Keratin High",
  as.numeric(as.character(clinical$K)) < 10 ~ "Keratin Low")
)
#MGC Content
clinical <- clinical %>% mutate(MGC_Three = case_when(
  as.numeric(as.character(clinical$MCG_total_mm2)) >= 1 ~ "MGC High",
  as.numeric(as.character(clinical$MCG_total_mm2)) < 0.2 ~ "MGC Low",
  as.numeric(as.character(clinical$MCG_total_mm2)) < 1 & as.numeric(as.character(clinical$MCG_total_mm2)) >= 0.2 ~ "MGC Medium")
)
# Transform data
rna <- read.csv('Tumor_HNSC_bulk.csv', row.names = 1)

rownames(rna) <- gsub("\\.\\d+$", "", rownames(rna))
ensembl.genes <- rownames(rna)

geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86, keys= ensembl.genes, keytype = "GENEID", columns = c("SYMBOL","GENEID"))

# Keep Only Unique ID Names
geneIDs1 <- geneIDs1 %>% distinct(SYMBOL, .keep_all= TRUE)

#Dropping geneIDs with no available conversion
rna <- rna[rownames(rna) %in% geneIDs1$GENEID,]
#Changing the rownames to symbols
rownames(rna) <- geneIDs1$SYMBOL

#Keeping only short sample IDs
colnames(rna) <- gsub('\\.','-',substr(colnames(rna),1,12))

# write into csv
write.csv(rna, 'Tumor_HNSC_bulk_PreProcessed.csv',row.names = T)

#Reading RNAseq data
rna <- read.csv('Tumor_HNSC_bulk_PreProcessed.csv', row.names = 1)
rna <- type.convert(rna, as.is = TRUE)
colnames(rna) <- gsub('\\.','-',substr(colnames(rna),1,12))

#Keeping only relevant bulk data
rna <- rna[,colnames(rna) %in% rownames(clinical)]

#Keep  only the patients with available bulk RNA seq
clinical_rna <- clinical[rownames(clinical) %in% colnames(rna),]

#Need to reorder the colnames and rownames to run DESeq2
rna = rna[,order(colnames(rna))]
clinical_rna = clinical_rna[order(rownames(clinical_rna)),]


ddsTxi <- DESeqDataSetFromMatrix(rna,
                                   colData = clinical_rna,
                                   design = ~MGC_Three)

keep <- rowSums(counts(ddsTxi)) >= 10
ddsTxi <- ddsTxi[keep,]

# run all DGE analysis

dds <- DESeq(ddsTxi)
MGC_Hi_Lo <- lfcShrink(dds, contrast=c('MGC_Three',"MGC High","MGC Low"), type = 'ashr')
MGC_Hi_Med <- lfcShrink(dds, contrast=c('MGC_Three',"MGC High","MGC Medium"), type = 'ashr')
MGC_Med_Lo <- lfcShrink(dds, contrast=c('MGC_Three',"MGC Medium","MGC Low"), type = 'ashr')

# bind all comparison analysis results

statistics <- as.data.frame(MGC_Hi_Lo) %>% rownames_to_column('symbol') %>% 
  dplyr::select(1,3,6) %>% dplyr::rename(log2FC_MGC_Hi_Lo = 2, fdr_MGC_Hi_Lo = 3) %>%
  left_join(as.data.frame(MGC_Hi_Med) %>% rownames_to_column('symbol') %>% 
          dplyr::select(1,3,6) %>% dplyr::rename(log2FC_MGC_Hi_Med = 2, fdr_MGC_Hi_Med = 3)) %>%
  left_join(as.data.frame(MGC_Med_Lo) %>% rownames_to_column('symbol') %>% 
          dplyr::select(1,3,6) %>% dplyr::rename(log2FC_MGC_Med_Lo = 2, fdr_MGC_Med_Lo = 3))


# MGC HIGH vs Low
hi_lo_rank <- as.data.frame(MGC_Hi_Lo) %>%
  dplyr::mutate(rank = -log10(pvalue)*sign(log2FoldChange)) %>% na.omit() %>% rownames_to_column("Gene") %>%
  dplyr::select(Gene, rank) %>%
  arrange(desc(rank)) %>%
  deframe()

# MGC Int vs Low
hi_med_rank <- as.data.frame(MGC_Hi_Med) %>%
  dplyr::mutate(rank = -log10(pvalue)*sign(log2FoldChange)) %>% na.omit() %>% rownames_to_column("Gene") %>%
  dplyr::select(Gene, rank) %>%
  arrange(desc(rank)) %>%
  deframe()

# MGC High vs Int
med_lo_rank <- as.data.frame(MGC_Med_Lo) %>%
  dplyr::mutate(rank = -log10(pvalue)*sign(log2FoldChange)) %>% na.omit() %>% rownames_to_column("Gene") %>%
  dplyr::select(Gene, rank) %>%
  arrange(desc(rank)) %>%
  deframe()

# calculate GSEA and generate plots 

MGC_visium <- data.frame(list = 'Visum_Genes', symbol =Sig_MGC)

set.seed(123)
gsea_hi_lo <- GSEA(geneList = hi_lo_rank, TERM2GENE = MGC_visium, seed = 1, nPerm = 10000)
anno <- gsea_hi_lo[1, c("NES", "p.adjust")]
lab <- paste0(names(anno), "=",  round(anno, 3), collapse="\n")

hi_lo_plot <- gseaplot2(gsea_hi_lo, geneSetID = 1, title = "High vs Low MGC (Visium Signature)") #+ annotate("text", 0, gsea_hi_lo[1, "enrichmentScore"], label = lab, hjust=-1.2, vjust=1.2)

gsea_hi_med <- GSEA(geneList = hi_med_rank, TERM2GENE = MGC_visium, seed = 1, nPerm = 10000)
anno <- gsea_hi_med[1, c("NES", "p.adjust")]
lab <- paste0(names(anno), "=",  round(anno, 3), collapse="\n")

hi_med_plot <- gseaplot2(gsea_hi_med, geneSetID = 1, title = "High vs Med MGC (Visium Signature)") +
  annotate("text", 0, gsea_hi_med[1, "enrichmentScore"], label = lab, hjust=-1.2, vjust=1)

gsea_med_lo <- GSEA(geneList = med_lo_rank, TERM2GENE = MGC_visium, seed = 1, nPerm = 10000, pvalueCutoff = 1)
anno <- gsea_med_lo[1, c("NES", "p.adjust")]
lab <- paste0(names(anno), "=",  round(anno, 3), collapse="\n")

med_lo_plot <- gseaplot2(gsea_med_lo, geneSetID = 1, title = "Med vs Low MGC (Visium Signature)") +
  annotate("text", 0, gsea_med_lo[1, "enrichmentScore"], label = lab, hjust=-1.2, vjust=0)

data.frame(list = c("Hi_Lo", 'Hi_Med'),
           NES = c(gsea_hi_lo@result$NES, gsea_hi_med@result$NES),
           FDR = c(gsea_hi_lo@result$p.adjust, gsea_hi_med@result$p.adjust))

cowplot::plot_grid(hi_lo_plot)
cowplot::plot_grid(hi_lo_plot, hi_med_plot, ncol=1, labels=LETTERS[1:2])
dev.off()
# should work !


fwrite(Print,"Blank_rank.rnk", sep = '\t')



fwrite(statistics, file = "stats_deseq_groups.txt", sep = '\t')
saveRDS(dds2, file = "dds.RDS")
