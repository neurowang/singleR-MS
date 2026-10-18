# single cell analysis

Analysis of four ALS/control microglia samples, covering integration, cell-state annotation, differential expression, BSG comparisons, enrichment, and glycolysis expression summaries. Each stage runs separately and exports numerical results.

## Files

| File | Purpose |
| --- | --- |
| `single cell analysis.R` | Analysis parameters, helpers, and stage functions |
| `scripts/run_hub.R` | Command-line entry point |
| `scripts/check_syntax.R` | R syntax and definition-loading check |
| [DATA_REQUIREMENTS.md](DATA_REQUIREMENTS.md) | Input objects, sample mapping, and data-access details |
| [PARAMETERS.md](PARAMETERS.md) | Analysis settings and statistical definitions |
| [OUTPUTS.md](OUTPUTS.md) | Output files and their interpretation |
| [DEPENDENCIES.md](DEPENDENCIES.md) | Software requirements and installation |
| [VALIDATION.md](VALIDATION.md) | Validation results and scope |
| [CHANGELOG.md](CHANGELOG.md) | Release notes |

## Analysis stages

| Stage | Computation | Default input |
| --- | --- | --- |
| `integration` | RNA normalization, SCTransform, PCA, Harmony, clustering, and UMAP coordinates | `work_dir/ALS_Microglia_QC_Annotated.rds` |
| `annotation` | Fixed cluster-to-state mapping, marker summaries, and cell composition | `work_dir/sc_hd_wgcna_completed.rds` |
| `de` | All-cell ALS/control and ALS subgroup contrasts | `work_dir/sc_hd_wgcna_completed.rds` |
| `bsg` | State-specific BSG expression summaries and Wilcoxon tests | `output_dir/annotation/ALS_Microglia_State_Annotated.rds` |
| `enrichment` | Late/Core disease KEGG enrichment and BSG Hallmark GSEA | Same annotated RDS |
| `metabolism` | Late/Core BSG grouping and glycolysis expression summaries | Same annotated RDS |

The default data directory is `~/rawData`, or `STUDY_WORK_DIR` set before loading the script. The default output directory is `work_dir/publication_outputs/hub_v3`. Use explicit absolute paths for each run.

Annotation uses a fixed study-specific map. Confirm the input's clustering and active assay against [DATA_REQUIREMENTS.md](DATA_REQUIREMENTS.md), especially when using a newly integrated object.

## Command-line execution

Install the dependencies in [DEPENDENCIES.md](DEPENDENCIES.md). From the package directory:

```sh
Rscript --vanilla scripts/run_hub.R --help
Rscript --vanilla scripts/check_syntax.R
```

Replace the example paths with your input paths and a fresh result directory. `--vanilla` bypasses local R startup files. For an existing clustered study object:

```sh
Rscript --vanilla scripts/run_hub.R --stage annotation \
  --work-dir /absolute/path/to/study_data \
  --output-dir /absolute/path/to/results/hub_v3 \
  --input-rds /absolute/path/to/study_data/sc_hd_wgcna_completed.rds

Rscript --vanilla scripts/run_hub.R --stage de \
  --work-dir /absolute/path/to/study_data \
  --output-dir /absolute/path/to/results/hub_v3 \
  --input-rds /absolute/path/to/study_data/sc_hd_wgcna_completed.rds

Rscript --vanilla scripts/run_hub.R --stage bsg \
  --work-dir /absolute/path/to/study_data \
  --output-dir /absolute/path/to/results/hub_v3

Rscript --vanilla scripts/run_hub.R --stage enrichment \
  --work-dir /absolute/path/to/study_data \
  --output-dir /absolute/path/to/results/hub_v3

Rscript --vanilla scripts/run_hub.R --stage metabolism \
  --work-dir /absolute/path/to/study_data \
  --output-dir /absolute/path/to/results/hub_v3
```

The last three commands use the annotation output. Use `--input-rds` to select an alternative annotated object.

To compute integration from a merged Seurat object with RNA counts:

```sh
Rscript --vanilla scripts/run_hub.R --stage integration \
  --work-dir /absolute/path/to/study_data \
  --output-dir /absolute/path/to/results/hub_v3 \
  --input-rds /absolute/path/to/study_data/ALS_Microglia_QC_Annotated.rds
```

Integration writes `integration/data/processed/ALS_Microglia_Harmony_Clustered.rds`. After verifying its cluster-to-state mapping and assay context, select that file explicitly with `--input-rds` for both `annotation` and `de`.

## Interactive R execution

Load the definitions and select one stage:

```r
source("single cell analysis.R")

study_data_dir <- "/absolute/path/to/study_data"
study_results_dir <- "/absolute/path/to/results/hub_v3"
study_article_rds <- file.path(study_data_dir, "sc_hd_wgcna_completed.rds")

run_hub("annotation", work_dir = study_data_dir,
        output_dir = study_results_dir, input_rds = study_article_rds)
```

Use `run_hub(stage, work_dir, output_dir, input_rds)` for the other stages. Scientific settings are defined in the main script and listed in [PARAMETERS.md](PARAMETERS.md).

## Results and validation

Results are written under `output_dir/<stage>/`, together with input, parameter, session, random-number-state, and completion records. Inspect `stage_status.csv` and the per-analysis status tables described in [OUTPUTS.md](OUTPUTS.md). Use a fresh output directory for each run because reruns overwrite generated files while older conditional files can remain.

Static checks are documented in [VALIDATION.md](VALIDATION.md). No real study RDS was used for these checks; numerical agreement with study results remains unverified.
