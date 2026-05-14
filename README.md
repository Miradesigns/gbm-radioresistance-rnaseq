# GBM Radioresistance — Transcriptomic Analysis

Differential gene expression analysis of glioblastoma multiforme (GBM) RNA-seq data, focused on identifying novel genes associated with radioresistance through the integrin αvβ3 and thyroid hormone signalling axis.

This project was conducted independently as self-directed research. All code was written from scratch without institutional support or supervision.

---

## Background

Radioresistance is one of the principal reasons GBM treatment fails. The integrin αvβ3 receptor and thyroid hormone (L-T4) signalling axis have been implicated in promoting tumour cell survival after radiation, with downstream effectors including PARP15 in AML and GBM contexts (Glinsky et al., *Endocrine Research*, 2024).

This analysis asks: which genes are differentially expressed in GBM tumours versus normal controls, and which of those genes have not previously been linked to GBM in any major disease database?

---

## Repository Structure

```
gbm-radioresistance-rnaseq/
│
├── README.md
│
├── 01_preprocessing_and_DESeq2/
│   └── RNA-Seq_Data_Analysis.ipynb       # Full R pipeline: load, QC, DESeq2, volcano, heatmap
│
├── 02_gene_id_annotation/
│   ├── ensembl_gene_id_to_name_converter.ipynb          # ENSG → gene symbol via MyGene.info
│   ├── ensembl_gene_id_to_name_converter_FIXED.ipynb    # Improved: Ensembl REST + MyGene fallback
│   └── unmapped_ensembl_gene_name_recovery.ipynb        # 4-source recovery for unannotated IDs
│
├── 03_go_pathway_analysis/
│   ├── GO_Pathway_Analysis_DEGs.ipynb    # GO/KEGG/Reactome enrichment (original)
│   ├── GO_Pathway_Analysis_DEGs_v2.ipynb # Full 3-tier analysis including RP11 loci
│   └── Reduction_GO.ipynb               # GO semantic similarity reduction + psoriasis filtering
│
├── data/
│   ├── DESeq2_filtered.txt
│   ├── GO_Biological_Process_2026_table.txt
│   ├── GO_Molecular_Function_2026_table.txt
│   └── GO_Cellular_Component_2026_table.txt
│
├── results/
│   ├── figures/
│   │   ├── raw_counts_barplot.png
│   │   ├── pca_plot.png
│   │   ├── pca_plot_final.png
│   │   └── Dendrogram_Heatmap.png
│   └── tables/
│       └── Psoriasis_relevance_filtered_enrichment_results.xlsx
│
└── docs/
    └── venn_diagram.png
```

---

## Analysis Pipeline

### Step 1 — Preprocessing & Quality Control
**Notebook:** `01_preprocessing_and_DESeq2/RNA-Seq_Data_Analysis.ipynb`

- Loaded raw RNA-seq count matrix (GBM tumour samples vs normal controls)
- Removed low-expression genes (rowSums ≤ 5) and duplicate gene entries
- Saved raw read count distribution across all samples
- Performed PCA for dimensionality reduction and sample quality assessment
- Identified and removed outlier samples before differential expression analysis

### Step 2 — Differential Expression Analysis (DESeq2)
**Notebook:** `01_preprocessing_and_DESeq2/RNA-Seq_Data_Analysis.ipynb`

- Ran DESeq2 with GBM vs Control comparison
- Generated log2 fold-change, p-values, and Benjamini-Hochberg adjusted p-values
- Initial significance filter: padj < 0.05, |log2FC| > 1 → **401 significant DEGs**
- Stringent filter: log2FC > 4 or < −2, padj < 0.01 → **367 upregulated, 34 downregulated**
- Visualized results via volcano plot, MA plot, and dendrogram heatmap of top 10 up/downregulated genes

### Step 3 — Database Comparison & Novel DEG Identification
- Exported DEG list and cross-referenced against:
  - [CTD](http://ctdbase.org/) — Comparative Toxicogenomics Database
  - [Open Targets](https://www.opentargets.org/) — GBM-associated genes
  - [GeneCards](https://www.genecards.org/) — curated gene function
- Generated four-way Venn diagram using [Venny](https://bioinfogp.cnb.csic.es/tools/venny/)
- **12 DEGs were exclusive** to our dataset — absent from all three databases

### Step 4 — Gene ID Annotation
**Notebooks:** `02_gene_id_annotation/`

Of the 12 exclusive DEGs:
- 5 carried valid HGNC-approved gene symbols
- 1 had an Ensembl ID only (snoU13, ENSG00000238683) — no Entrez record
- 6 carried retired RP11 clone-based identifiers, absent from all current databases

Multiple annotation strategies were applied in sequence:
1. Ensembl REST API (current GRCh38)
2. Ensembl GRCh37 archive
3. MyGene.info
4. NCBI Entrez Gene via Biopython

### Step 5 — GO & Pathway Enrichment Analysis
**Notebooks:** `03_go_pathway_analysis/`

A three-tier strategy was applied based on gene annotation status:

| Tier | Genes | Method |
|------|-------|--------|
| 1 | PARP15, TMIGD3 | Direct GO enrichment (protein-coding) |
| 1 | LUADT1, LINC02728, LOC100128906 | Guilt-by-association (literature surrogates + co-expression) |
| 2 | snoU13 | Surrogate analysis using snoRNA/ribosome biogenesis partners |
| 3 | RP11-22P4.1, RP11-473M20.9, RP11-706O15.1, RP11-34P13.13, RP11-152L20.3, RP11-521D12.5 | Genomic neighborhood analysis via biomaRt (±500kb) |

Tools used: `clusterProfiler`, `ReactomePA`, `enrichplot`, `biomaRt`, `org.Hs.eg.db`

### Step 6 — GO Reduction & Psoriasis Relevance Filtering
**Notebook:** `03_go_pathway_analysis/Reduction_GO.ipynb`

- Applied semantic similarity reduction using `rrvgo` to collapse redundant GO terms
- Filtered enriched terms against a curated psoriasis-relevant keyword list
- TMIGD3 mapped to inflammatory response and cell proliferation regulation — both core pathways in psoriasis pathophysiology
- PARP15 mapped to RNA Pol II transcriptional regulation, flagged for further review

---

## Key Findings

| Gene | Type | Top GO Term | Biological Relevance |
|------|------|-------------|----------------------|
| PARP15 | Protein-coding | Protein poly-ADP-ribosylation (padj = 0.019) | DNA damage response; links to AML/GBM radioresistance |
| TMIGD3 | Protein-coding | Activation of Adenylate Cyclase Activity (padj = 0.019) | NF-κB/Akt tumour suppressor; inflammatory response |
| LUADT1 | lncRNA | PRC2/SUZ12 axis | Epigenetic p27 suppression; cell cycle regulation |
| snoU13 | snoRNA | Ribosome biogenesis (surrogate) | Pre-rRNA processing; DDX5/DDX17 regulated |
| RP11-521D12.5 | Unannotated | Genomic neighborhood pending | Not found in any current database — novel locus |

---

## Requirements

### R packages
```r
BiocManager::install(c(
  "DESeq2", "clusterProfiler", "ReactomePA", "enrichplot",
  "org.Hs.eg.db", "AnnotationDbi", "biomaRt",
  "EnhancedVolcano", "pheatmap", "rrvgo", "GOSemSim"
))
install.packages(c("ggplot2", "dplyr", "stringr"))
```

### Python packages
```bash
pip install pandas openpyxl requests mygene biopython
```

---

## Data Availability

The raw count matrix used in this analysis is GBM RNA-seq data. If you need access to the source dataset, please open an issue or contact via the email on my profile.

The filtered DESeq2 results (`data/DESeq2_filtered.txt`) and GO input tables are included in this repository.

---

## Reference

Glinsky et al. (2024). *Additional considerations in cancer cell radioresistance, integrin αvβ3 and thyroid hormones.* Endocrine Research. DOI: 10.1080/07435800.2024.2361152

---

## Author

**Miracle Okonkwo**
Enugu, Nigeria
B.Sc. Computer Science and Statistics, University of Nigeria Nsukka

*Independent research project*
