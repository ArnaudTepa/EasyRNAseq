<img width="700" height="700" alt="Designer (14)" src="https://github.com/user-attachments/assets/21f655eb-ac41-4ad5-9b33-dc67d1bb6143" />

# EasyRNAseq: A Complete Pipeline for RNA-seq and Polymorphism Analysis

Welcome to **EasyRNAseq**, an open-source pipeline designed to simplify and standardize RNA-seq analysis for insecticide resistance studies.  

This repository contains all scripts, workflows, and documentation used in the study:

**Insights from transcriptomic profiling identify the CYP6Z gene family as key drivers of pyrethroid resistance escalation in *Anopheles gambiae* from Cameroon**  
Published in *BMC Genomics* (2025)

---

## 🚀 What EasyRNAseq Offers

- **End-to-End Workflow**: From raw FASTQ files to publication-ready figures.
- **Modules Included**:
  - ✅ Quality Control (FastQC, MultiQC)
  - ✅ Read Mapping (HISAT2)
  - ✅ Gene Expression Analysis (DESeq2)
  - ✅ Functional Enrichment (GO term analysis)
  - ✅ Polymorphism Detection (Samtools, VarScan2, SnpEff)
  - ✅ Population Genetics Metrics (FST, Tajima’s D, π)
  - ✅ Visualization Tools (Volcano plots, heatmaps, selection scans)
- **Reproducibility**: All scripts and parameter files are included.
- **Flexibility**: Modular design for RNA-seq, GO analysis, and SNP-based population genetics.

---

## 📂 Data & Code Availability

- Raw RNA-seq data: ENA accession **[PRJEB97187]**
- Pipeline scripts & configs: Available in this repository
- Environment setup: Provided via `environment.yml` for Conda
- Full workflow documentation: See **`Wiki`**

---

## ⚡ Quick Start

### 1. Clone the repository
```bash
git clone https://github.com/ArnaudTepa/EasyRNAseq.git
cd EasyRNAseq
```

### 2. Create the environment
```bash
conda env create -f environment.yml
conda activate easyrnaseq_env
```

### 3. Run the full pipeline
```bash
bash run_all.sh
```

---

## ❓ Why EasyRNAseq?

- **Transparent**: Complete workflow for RNA-seq and polymorphism analysis.  
- **Open Science**: Code and data are publicly available for reuse and adaptation.  
- **Optimized for Resistance Studies**: Includes candidate gene analysis and population genetics metrics.  
```
