## 🚀 Pipeline for RNA‑seq & Polymorphism Analysis

Below is a step‑by‑step command‑line workflow with inline comments (`#`).  
**Prerequisites**: Activate your conda env (e.g., `conda activate agam_resist_env`) and verify all tools are on PATH.

---

### 0) Recommended folder layout and variables

# Set once per session to make commands shorter and reproducible
PROJECT_ROOT=/ABSOLUTE/PATH/TO/PROJECT
RAW=${PROJECT_ROOT}/data/fastq
REF=${PROJECT_ROOT}/config/reference
OUT=${PROJECT_ROOT}/outputs
BAM=${OUT}/bam
COUNTS=${OUT}/counts
VAR=${OUT}/variants
PLOTS=${PROJECT_ROOT}/figures
META=${PROJECT_ROOT}/data/metadata

# Create output dirs if they don't exist
mkdir -p "$OUT" "$BAM" "$COUNTS" "$VAR" "$PLOTS" "$META"

---

### 1) Quality checking

# NOTE: It is advised to (re)name raw files consistently at the beginning,
# e.g., SAMPLEID_R1.fastq.gz / SAMPLEID_R2.fastq.gz

# Run FastQC on all FASTQ files and summarize with MultiQC
# USAGE: fastqc.sh INPUT_DIR OUTPUT_DIR
bash path/to/fastqc.sh "$RAW" "$OUT/qc"

# Optional: run MultiQC if not included inside fastqc.sh
multiqc "$OUT/qc" -o "$OUT/qc"

---

### 2) Mapping with HISAT2

# Proceed only if QC is acceptable (adapter/quality issues resolved)
# USAGE: mapping_hisat2.sh HISAT2_INDEX_DIR READS_DIR OUTPUT_BAM_DIR
bash path/to/mapping_hisat2.sh "$REF/hisat2_index" "$RAW" "$BAM"

# Ensure coordinate-sorted, indexed BAMs exist (samtools sort/index if needed)
# Example:
# samtools sort -@ 8 -o ${BAM}/SAMPLE.sorted.bam ${BAM}/SAMPLE.bam
# samtools index ${BAM}/SAMPLE.sorted.bam

---

### 3) Generate counts with featureCounts and curate

# Thread count (-T), paired-end (-p), multi-mapping (-M), multi-overlap (-O)
# Count CDS features (-t CDS), grouping by gene_id attribute (-g gene_id)
# Adjust -a to GTF/GFF file as appropriate
featureCounts -T 40 -p -M -O \
  -t CDS -g gene_id \
  -a ${REF}/annotation.gff \
  -o ${COUNTS}/COUNTS_raw.txt \
  --countReadPairs ${BAM}/*.sorted.bam

# Curate table:
# - remove first header line, keep feature column + sample columns
# - clean path prefixes so column headers are just sample names
cat ${COUNTS}/COUNTS_raw.txt \
  | tail -n +2 \
  | cut -f1,6- \
  | perl -ne '$_=~s/\.\.\/alignment\///g; print $_' \
  > ${COUNTS}/COUNTS_parsed.txt

---

### 4) Differential expression with DESeq2 and figure prep

# NOTE: Check parameters (contrasts, thresholds) and ensure all required
# inputs (sample sheet, group labels) follow the provided templates.

# Run DESeq2 contrasts (R–S, R–C, C–S; and dose: 10×–1×, 5×–1×, 10×–5×)
# This script should read counts + metadata and produce tables for each contrast.
Rscript path/to/Script_DESEQ2.R \
  --counts ${COUNTS}/COUNTS_parsed.txt \
  --meta ${META}/sample_sheet.csv \
  --out ${OUT}/deseq2

# Inspect output tables; shortlist candidate genes as needed.

# Generate volcano plots and Venn diagrams for the selected contrasts
Rscript path/to/Script_volcano_venn_plot.R \
  --deseq_dir ${OUT}/deseq2 \
  --out ${PLOTS}/fig1_fig2

---

### 5) (Optional) Merge BAM files

# If you need merged BAMs per group/condition (e.g., for mpileup from pooled data)
# USAGE: Merged_bam.sh INPUT_BAM_DIR
bash path/to/Merged_bam.sh "$BAM"

---

### 6) Variant calling: mpileup → sync → VCF

# This wrapper should:
#  - create mpileup using samtools (with appropriate min MAPQ/baseQ)
#  - create a sync file (pooled format) if using PoPoolation2 tools
#  - call variants with VarScan2 to produce a VCF
# USAGE:
# Variant_call.sh POP2_MP2SYNC SCRIPT  REF_FASTA REF_FASTA_FAI BAM_DIR WORK_DIR VARSCAN_JAR
bash path/to/Variant_call.sh \
  path/to/popoolation2_1201/mpileup2sync.pl \
  ${REF}/AgamP4.fa \
  ${REF}/AgamP4.fa.fai \
  "$BAM" \
  "$VAR" \
  path/to/VarScan.jar

---

### 7) Fisher’s exact test (PoPoolation2)

# Compare allele frequencies between groups (e.g., resistant vs control)
# USAGE: fisher.sh SYNC_FILE POP2_FISHER_TEST OUTPUT_DIR SAMPLE_DEFINITION
# SAMPLE_DEFINITION lists groups/contrasts (see your template)
bash path/to/fisher.sh \
  ${VAR}/all_samples.sync \
  path/to/popoolation2_1201/fisher-test.pl \
  ${VAR}/fisher \
  ${META}/sample_contrasts.txt

---

### 8) Annotation & filtering with SnpEff / SnpSift

# 8.1 Rename samples in VCF for clarity
# USAGE: rename_vcf.sh INPUT_VCF SAMPLE_MAP OUTPUT_VCF
# SAMPLE_MAP: two columns (old_id \t new_id)
bash path/to/rename_vcf.sh \
  ${VAR}/merged.vcf \
  ${META}/sample_map.txt \
  ${VAR}/renamed.vcf

# 8.2 Annotate with SnpEff
# -c: config file; -stats/-csvStats: summary outputs
# -chr: restrict to chromosomes of interest if desired (e.g., 2L, 2R, 3L, 3R, X)
# Ref_Name must match a valid SnpEff database entry configured in snpEff.config
cat ${VAR}/renamed.vcf \
  | java -Xmx16g -jar path/to/snpEff/snpEff.jar \
      -c path/to/snpEff/snpEff.config \
      -csvStats ${VAR}/renamed.snpeff.csv \
      -stats    ${VAR}/renamed.snpeff.html \
      -chr 2L,2R,3L,3R,X \
      -v AgamP4 \
  > ${VAR}/renamed.ann.vcf

# 8.3 Extract TSV with annotation + per-sample allele frequencies
cat ${VAR}/renamed.ann.vcf \
  | path/to/snpEff/scripts/vcfEffOnePerLine.pl \
  | java -jar path/to/snpEff/SnpSift.jar extractFields -e "NA" - \
      CHROM POS REF ALT \
      "ANN[*].EFFECT" "ANN[*].IMPACT" "ANN[*].HGVS_C" "ANN[*].HGVS_P" "ANN[*].GENE" "ANN[*].GENEID" \
      "GEN[*].FREQ" \
  > ${VAR}/renamed.ann.tsv

---

### 9) PCA and clustering (candidate transcripts subset + PLINK)

# Filter VCF by candidate transcript regions using a BED (chrom start end)
bgzip -c ${VAR}/renamed.vcf > ${VAR}/renamed.vcf.gz
bcftools index -t ${VAR}/renamed.vcf.gz

bcftools view -R ${META}/selected_candidates.bed \
  ${VAR}/renamed.vcf.gz -O v -o ${VAR}/renamed_candidates.vcf

# PLINK pipeline
VCF=${VAR}/renamed_candidates.vcf

# 9.1 Convert to PLINK bed/bim/fam (retain extra chromosomes if needed)
plink \
  --vcf "$VCF" \
  --double-id \
  --allow-extra-chr \
  --set-missing-var-ids @:# \
  --make-bed \
  --out ${VAR}/plink_nofilter

# 9.2 LD pruning (window=50 SNPs, step=10, r^2 threshold=0.1; tune as needed)
plink \
  --bfile ${VAR}/plink_nofilter \
  --indep-pairwise 50 10 0.1 \
  --out ${VAR}/plink_prune \
  --allow-extra-chr

# 9.3 PCA on pruned variants
plink \
  --bfile ${VAR}/plink_nofilter \
  --extract ${VAR}/plink_prune.prune.in \
  --pca 20 \
  --out ${VAR}/plink_PCA \
  --allow-extra-chr

# 9.4 Distance matrix (optional)
plink \
  --bfile ${VAR}/plink_nofilter \
  --extract ${VAR}/plink_prune.prune.in \
  --distance square \
  --const-fid \
  --out ${VAR}/plink_distance \
  --allow-extra-chr

# 9.5 Plot structure / PCA in R (provide your plotting script)
Rscript path/to/Script_population_structure.R \
  --pca ${VAR}/plink_PCA.eigenvec \
  --out ${PLOTS}/structure

---

### 10) Diversity metrics for candidate genes (Grenedalf)

# Ensure a clean sync file (sorted, unique POS per CHR)
sort -k1,1 -k2,2n ${VAR}/all_samples.sync > ${VAR}/all_samples.sorted.sync
awk '!seen[$1,$2]++' ${VAR}/all_samples.sorted.sync > ${VAR}/all_samples.sorted.clean.sync

# FST across candidate intervals
# USAGE: fst_grenedalf_sync.sh SYNC OUT_DIR POOL_SIZE THREADS PREFIX GRENEDALF_PATH BED FILE_SUFFIX SAMPLE_FILE
bash path/to/fst_grenedalf_sync.sh \
  ${VAR}/all_samples.sorted.clean.sync \
  ${VAR}/grenedalf \
  POOL_SIZE_FILE_OR_VALUE \
  8 \
  fst_candidates \
  path/to/grenedalf \
  ${META}/candidates.bed \
  _cand \
  ${META}/sample_groups.txt

# Nucleotide diversity (π) and Tajima’s D across candidates
bash path/to/diversity_grenedalf_sync.sh \
  ${VAR}/all_samples.sorted.clean.sync \
  ${VAR}/grenedalf \
  POOL_SIZE_FILE_OR_VALUE \
  8 \
  div_candidates \
  path/to/grenedalf \
  ${META}/candidates.bed \
  _cand \
  ${META}/sample_groups.txt

---

### 11) Visualization: FST, π, Tajima’s D, SNP heatmaps

# Produce integrated plots (selection scans, diversity summaries, SNP heatmap with p-values)
Rscript path/to/visual_TWAS.R \
  --fst ${VAR}/grenedalf/fst_candidates_results.tsv \
  --pi  ${VAR}/grenedalf/div_candidates_pi.tsv \
  --tajimasd ${VAR}/grenedalf/div_candidates_tajd.tsv \
  --snp ${VAR}/renamed.ann.tsv \
  --out ${PLOTS}/fig3_fig5

---

## Notes & Best Practices

# DEG thresholds used in the paper: FDR < 0.05; abs(FC) ≥ 2 (field vs lab) and ≥ 1.5 (within-field dose response).
# RNA‑seq SNPs are limited to expressed regions and may be affected by allele‑specific expression (ASE); interpret selection cautiously.
# For VGSC/kdr, RNA‑seq may miss some variants; confirm with targeted genotyping where possible.
# Keep a config/params.yaml with paths, thresholds, and sample grouping for reproducibility.
# Version‑control your outputs (figures/tables) by date-stamped folders or Git tags.

