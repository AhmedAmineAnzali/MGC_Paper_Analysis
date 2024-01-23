# MGC_Paper_Analysis

## Code to Reproduce Bulk RNAseq Data Analysis

To reproduce the bulk RNAseq data analysis described in the paper titled "Trem2 Giant Macrophages: Biomarker of Good Prognosis in SCC," follow the steps below:

1. Access the rendered HTML for the analysis and figures at: [MGC_Paper_Analysis](https://ahmedamineanzali.github.io/MGC_Paper_Analysis/RNAseq).

2. **Data Retrieval:**
   - Head and neck SCC patients with tongue SCC data were retrieved using TCGAbiolinks.
   - Clinical data for the TCGA cohort was obtained from the supplementary material of the paper by Liu et al., 2018[^1^].
     Progression-Free Interval (PFI) and Overall Survival (OS) were used as time points for plotting, as recommended in the same paper.
   - Survival plots were constructed using the survival and survminer R packages.

3. **Data Preprocessing:**
   - Bulk RNAseq data were normalized using the EDASeq package based on gene length.
   - Genes with zero expression in more than 25% of the samples were removed.

4. **Differential Expression Analysis:**
   - The DESeq2 package was used for differential expression testing between MGC High and Low patients.
   - Volcano plots were constructed using the EnhancedVolcano package.
   - The R package clusterprofiler was employed to analyze the Differentially Expressed Genes (DEGs) between different groups.

For detailed code and implementation, refer to the analysis [here](https://ahmedamineanzali.github.io/MGC_Paper_Analysis/RNAseq).

## Code to Reproduce Visium Data Analysis

Nine patients were selected for spatial transcriptomics, including 6 MGC high and 3 MGC low cases. Additionally, giant cell-only spots were annotated by a pathologist. The different slides were preprocessed using the Seurat V4 pipeline, Sctransformed, and integrated by RPCA.

For detailed code and implementation, refer to the analysis [here](https://ahmedamineanzali.github.io/MGC_Paper_Analysis/Visium).

[^1^]: Liu J, Lichtenberg T, Hoadley KA, Poisson LM, Lazar AJ, Cherniack AD, Kovatich AJ, Benz CC, Levine DA, Lee AV, Omberg L, Wolf DM, Shriver CD, Thorsson V; Cancer Genome Atlas Research Network; Hu H. An Integrated TCGA Pan-Cancer Clinical Data Resource to Drive High-Quality Survival Outcome Analytics. Cell. 2018 Apr 5;173(2):400-416.e11. doi: 10.1016/j.cell.2018.02.052. PMID: 29625055; PMCID: PMC6066282.


