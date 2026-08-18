# Watershed-SV

Watershed-SV extends [Watershed](https://github.com/BennyStrobes/Watershed) to model the
impact of rare structural variants (DUP, DEL, DUP-CNV, DEL-CNV, INV, INS) on nearby
molecular outliers. This repository packages the annotation pipeline, the multi-omic
integration steps, and the model-training/visualization wrapper used to run Watershed-SV
across expression and protein traits.

This repository is a fork of, and builds on, the original Watershed-SV implementation by
[**jasonbhn/Watershed-SV**](https://github.com/jasonbhn/Watershed-SV). The base annotation
and integration engine (`scripts/executable_scripts/`, packaged as the `watershed-sv` conda
tools, and the `.github/` build recipes) is their work — please cite the original authors
for that core. This fork adds a cell-type–specific (ABC/E2G) annotation layer, a Polars-based
multi-omic integration step, and an end-to-end orchestration + visualization wrapper used in
the ADRC long-read WGS flagship study (see [Citation](#citation)).

> **Downstream analysis & figures.** The result-integration and plotting code that consumes
> the posterior probabilities produced here lives in
> [tjense25/ADRC-Flagship-lrGS-Code/WatershedSV_plotting](https://github.com/tjense25/ADRC-Flagship-lrGS-Code/tree/main/WatershedSV_plotting).

---

## What this fork contributes

The files below are the contribution of this fork; everything else is inherited from the
upstream repository and credited above.

| Area | File | Description |
|------|------|-------------|
| Cell-type–specific annotations | [`modified_scripts/generate_annotations_ABC`](modified_scripts/generate_annotations_ABC) | Extends the base annotation pipeline with a `-z/--cell-type` option that intersects rare SVs with cell-type ABC/E2G enhancer–gene maps ([e2g.stanford.edu](https://e2g.stanford.edu/)). |
| Multi-omic integration | [`modified_scripts/combine_all_annotations_ABC_polars`](modified_scripts/combine_all_annotations_ABC_polars) | Polars re-implementation of annotation collapsing to gene-level and SV-level tables, including the ABC/E2G features, for fast integration on large cohorts. |
| Supporting steps | [`modified_scripts/extract_gene_exec.sh`](modified_scripts/extract_gene_exec.sh), [`modified_scripts/run_extract_sv_vep_annotations.sh`](modified_scripts/run_extract_sv_vep_annotations.sh) | Gene-model extraction (pyranges/Polars) and the VEP annotation wrapper used by the pipeline. |
| Orchestration | [`WS_wrapper/driver_Watershed.sh`](WS_wrapper/driver_Watershed.sh) | End-to-end, checkpointed SLURM driver: annotation generation → multi-omic integration → Watershed training/evaluation/prediction → visualization. |
| Visualization | [`WS_wrapper/watershed_model_weights_rna.R`](WS_wrapper/watershed_model_weights_rna.R), [`WS_wrapper/watershed_model_weights_GAM_WS_performance_comparison.R`](WS_wrapper/watershed_model_weights_GAM_WS_performance_comparison.R) | Plots learned Watershed annotation weights and the Watershed-vs-GAM performance comparison. |

Although the annotations were designed for expression outliers, the integration and modeling
steps are trait-agnostic: any per-sample molecular z-score table can be supplied via
`--expressions` / `--expression-field`. In the flagship study this pipeline was run
independently for three omes — **RNA expression, plasma proteins, and CSF proteins**.

---

## Repository layout

```
Watershed-SV/
├── WS_wrapper/                 # This fork: end-to-end driver + visualization (R)
├── modified_scripts/           # This fork: ABC/E2G annotation + Polars integration
├── scripts/executable_scripts/ # Base pipeline (upstream jasonbhn/Watershed-SV)
├── .github/                    # Conda/container build recipes (upstream)
├── WatershedSV.yml             # Full conda environment used to collect annotations
└── env.yml                     # Minimal environment (watershed-sv + ensembl-vep)
```

---

## Installation

The base annotation pipeline is distributed as a conda package. We recommend `pixi`:

```bash
pixi global install -c dnachun -c conda-forge -c bioconda watershed-sv=0.1.22
pixi global install -c bioconda ensembl-vep
```

Or with conda:

```bash
conda install -c dnachun -c conda-forge -c bioconda watershed-sv=0.1.22
conda install -c bioconda ensembl-vep
```

To reproduce the full environment used to collect annotations, see
[`WatershedSV.yml`](WatershedSV.yml); a minimal environment is in [`env.yml`](env.yml).

**Reference data.** Before running, download:

- A VEP offline cache (GRCh38):
  ```bash
  wget https://ftp.ensembl.org/pub/current_variation/indexed_vep_cache/homo_sapiens_vep_113_GRCh38.tar.gz
  ```
- A GENCODE transcript model GTF (e.g. `gencode.v32.annotation.gtf.gz`).
- A genome-bound file (contig name, start, end).
- For the cell-type layer: ABC/E2G enhancer–gene maps for each cell type of interest.

Paths to these references are set as command-line arguments (see below) and, for the ABC
layer, in a small paths file described under [Cell-type–specific model](#cell-typespecific-model).

---

## Running the pipeline

### 1. Generate annotations

```bash
generate_annotations_ABC \
  -p population \        # or "smallset" for n < ~100
  -v <input.vcf> \       # VCF with ≥1 sample column; SVTYPEs DUP/DEL/DUP_CNV/DEL_CNV/CNV/INS/INV
  -f PASS \              # variant FILTERs to keep
  -k 100000 \            # flank up/downstream of genes (bp)
  -r 0.01 \              # rareness (MAF) threshold when --filter-rare
  -o <outdir> \
  -b <genome_bound_file> \
  -g <gencode.gtf.gz> \
  -c <vep_cache_dir> \
  -e False \             # filter by ethnicity (GTEx relic)
  -i True                # filter to rare variants (population model)
```

Key parameters:

| Flag | Meaning |
|------|---------|
| `-p, --pipeline` | `population` (n > ~100; enables `--filter-ethnicity`, `--filter-rare`) or `smallset`. |
| `-v, --input-vcf` | Input VCF, ≥1 sample column. Considered SVTYPEs: DUP, DEL, DUP_CNV, DEL_CNV, CNV, INS, INV. |
| `-f, --filters` | VCF FILTER values to retain. |
| `-k, --flank` | bp up/downstream of genes to consider (e.g. 100000 or 10000). |
| `-r, --rareness` | MAF threshold for rare variants when `--filter-rare` (≤ 0.01 recommended). |
| `-l, --liftover-bed` | Optional crossmap/liftover BED to convert VCF coordinates to GRCh38. |
| `-o, --outdir` | Output directory for annotations. |
| `-b, --genome-bound-file` | Contig name, start, end. |
| `-g, --gencode-genes` | GENCODE transcript model GTF. |
| `-c, --vep-cache-dir` | VEP offline cache directory. |
| `-a, --metadata` | Metadata file for ethnicity filtering (e.g. GTEx dbGaP metadata). |
| `-e, --filter-ethnicity` | Filter by ethnicity (GTEx relic). |
| `-i, --filter-rare` | Filter to rare variants (`population` model). |
| `-z, --cell-type` | *(this fork)* Comma-separated cell types for ABC/E2G annotation. |

### 2. Integrate annotations with molecular outliers

`combine_all_annotations_ABC_polars` collapses the annotations to **gene-level** and
**SV-level** tables and merges them with a per-sample molecular z-score table
(`--expressions`, `--expression-field`, `--expression-id-field`). `eval_watershed_prep`
then formats the merged data and assigns N2 pairs for Watershed. See
[`WS_wrapper/driver_Watershed.sh`](WS_wrapper/driver_Watershed.sh) for a fully worked
invocation of both steps.

### 3. Train, evaluate, and predict

The driver calls `evaluate_watershed.R` and `predict_watershed.R` from the upstream
[Watershed](https://github.com/BennyStrobes/Watershed) package (models: `Watershed_exact`,
`Watershed_approximate`, `RIVER`). Posterior probabilities from prediction are the primary
output consumed by the [downstream plotting repo](https://github.com/tjense25/ADRC-Flagship-lrGS-Code/tree/main/WatershedSV_plotting).

### 4. Visualize

- [`watershed_model_weights_rna.R`](WS_wrapper/watershed_model_weights_rna.R) — learned annotation weights from the prediction object.
- [`watershed_model_weights_GAM_WS_performance_comparison.R`](WS_wrapper/watershed_model_weights_GAM_WS_performance_comparison.R) — Watershed vs. GAM (WGS-only) performance from the evaluation object.

---

## Cell-type–specific model

To incorporate cell-type–specific regulatory information (ABC/E2G annotations from
[e2g.stanford.edu](https://e2g.stanford.edu/)):

1. Download the relevant ABC/E2G enhancer–gene map files.
2. Collect them into a single directory and point the pipeline at them via the paths file
   sourced in [`modified_scripts/generate_annotations_ABC`](modified_scripts/generate_annotations_ABC)
   (the `e2g_paths.sh` reference near line 359) — update this path for your environment.
3. Run `generate_annotations_ABC` with `-z <cell_type[,cell_type...]>`.

**Notes**
- Annotations must be generated separately for each cell type of interest.
- All reference/input paths must be set correctly before running.

---

## End-to-end automation

[`WS_wrapper/driver_Watershed.sh`](WS_wrapper/driver_Watershed.sh) runs the complete
workflow with per-step checkpointing (a `.done` marker lets re-runs resume):

1. Annotation generation
2. Multi-omic integration
3. Watershed training, evaluation, and prediction
4. Performance assessment and visualization

The driver's input/output paths and dataset parameters are environment-specific — edit the
variables at the top (VCF path, dataset name, reference paths, output directories, and the
molecular-trait/z-score inputs) before submitting.

---

## Citation

If you use this pipeline, please cite the original Watershed and Watershed-SV work:

- Watershed: [BennyStrobes/Watershed](https://github.com/BennyStrobes/Watershed).
- Watershed-SV (base implementation): [jasonbhn/Watershed-SV](https://github.com/jasonbhn/Watershed-SV).

The cell-type–specific extension and multi-omic application in this fork were developed for
the ADRC long-read WGS flagship study; the downstream integration and figures are in
[tjense25/ADRC-Flagship-lrGS-Code](https://github.com/tjense25/ADRC-Flagship-lrGS-Code).
