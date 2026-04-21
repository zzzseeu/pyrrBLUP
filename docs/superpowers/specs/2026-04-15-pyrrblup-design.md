# Python rrBLUP Alignment Design

## Goal

Build a Python implementation of the core `rrBLUP` workflow that can run end-to-end on the provided genotype and phenotype CSV files, using R `rrBLUP` outputs as the baseline for incremental numerical alignment.

## Scope

The first delivery targets a minimal but usable pipeline:

- `A.mat`
- `mixed.solve`
- `kin.blup`
- Real-data end-to-end execution on one phenotype location at a time
- Automated comparison against baseline outputs produced by R `rrBLUP`

The first delivery does not try to cover the full public surface area of the R package. It prioritizes:

1. Running the workflow successfully in the Python environment.
2. Establishing tests and baseline comparison infrastructure.
3. Tightening numerical agreement with R in later iterations.

## Inputs And Constraints

### Environments

- R baseline environment: `rrblup`
- Python development environment: `evo2`
- Environment activation must use `conda activate`
- Do not use `conda run`

### Data Files

- Genotype: `../../project/genomic_selection_project/data/processed/model_dataset/agront_ds_202603/genotype_012.csv`
- Phenotype: `../../project/genomic_selection_project/data/processed/phenotype/phenotype_ETN_processed_20260227.csv`

### Observed Data Shape

Genotype file:

- First column: `sample`
- Remaining columns: marker columns encoded as `0/1/2`
- One row per sample

Phenotype file:

- First column: `sample`
- Second column: `value`
- Remaining columns: one-hot location columns such as `loc_BeiJ`
- One sample can appear multiple times, one row per location-specific phenotype record

### First-Pass Modeling Choice

Development and testing will use one location at a time. The pipeline will filter phenotype rows by a selected `loc_*` column, then align the filtered phenotype samples with the genotype samples by sample name intersection.

## Recommended Approach

Use an R-baseline-driven, layered implementation strategy.

### Why This Approach

- It matches the requested priority: run the workflow first, then align numerically.
- It avoids guessing `rrBLUP` internals from memory.
- It creates a reusable regression harness for future compatibility work.

### Alternative Approaches Considered

#### 1. Python-first direct implementation

Implement the Python formulas immediately and only compare with R later.

Trade-off:

- Faster first code
- Higher risk of hidden formula or ordering mismatches

#### 2. Full API compatibility first

Attempt to mirror the R package interface broadly from the beginning.

Trade-off:

- Better long-term parity
- Slower delivery
- Too much surface area for the first milestone

#### 3. R-baseline-driven layered implementation

Generate baseline artifacts from R on the real data, then build Python to match those artifacts.

Trade-off:

- Slightly more upfront setup
- Best balance for fast end-to-end delivery plus later numerical alignment

This design chooses option 3.

## Architecture

The implementation is split into four layers:

### 1. Data Preparation Layer

Responsibilities:

- Read genotype and phenotype CSV files
- Validate required columns
- Select one phenotype location by location column
- Drop rows not matching the selected location
- Align phenotype and genotype by sample name
- Preserve deterministic sample ordering for all downstream outputs

Reasoning:

Most alignment failures against R are caused by mismatched filtering or ordering rather than solver logic. This layer needs to be explicit and testable.

### 2. Baseline Generation Layer

Responsibilities:

- Run R `rrBLUP` in the `rrblup` environment
- Produce baseline artifacts for the same filtered dataset
- Serialize matrices, vectors, metadata, and sample ordering for reuse by Python tests

Reasoning:

The R package is the source of truth for iterative alignment. The baseline artifacts need to be stable, inspectable, and reusable.

### 3. Core Algorithm Layer

Responsibilities:

- Compute additive relationship matrix via Python `A.mat`
- Solve mixed model equations via Python `mixed.solve`
- Build the single-location prediction flow via Python `kin.blup`

Reasoning:

These functions are the minimal core needed for a working rrBLUP-style pipeline and align with the requested initial scope.

### 4. Verification Layer

Responsibilities:

- Small deterministic unit tests
- Real-data smoke tests
- R baseline regression tests with staged tolerance tightening

Reasoning:

The first milestone must prove both operational correctness and future readiness for numerical alignment.

## Repository Structure

```text
pyrrBLUP/
├── baselines/
├── docs/
│   └── superpowers/
│       └── specs/
├── scripts/
├── src/
│   └── pyrrblup/
└── tests/
```

Planned responsibilities:

- `src/pyrrblup/amat.py`: additive relationship matrix implementation
- `src/pyrrblup/mixed_solve.py`: mixed model solver implementation
- `src/pyrrblup/kin_blup.py`: single-location prediction pipeline
- `src/pyrrblup/io.py`: CSV loading, location filtering, sample alignment
- `src/pyrrblup/models.py`: structured result containers
- `scripts/build_r_baseline.R`: generate baseline outputs from R `rrBLUP`
- `scripts/run_python_pipeline.py`: execute the Python pipeline on the same inputs
- `tests/test_io.py`: data filtering and alignment tests
- `tests/test_amat.py`: matrix property and baseline tests
- `tests/test_mixed_solve.py`: solver behavior and baseline tests
- `tests/test_kin_blup.py`: pipeline behavior tests
- `tests/test_against_r_baseline.py`: real-data regression tests

## Baseline Artifact Design

The baseline outputs produced by R should include:

- Selected location name
- Filtered sample list in exact order
- Filtered phenotype vector
- Marker sample list in exact order after alignment
- `A.mat` output matrix
- `mixed.solve` key results:
  - fixed effects estimates
  - random effects estimates
  - variance component estimates
- `kin.blup` predictions
- Metadata describing:
  - selected location column
  - filtering rules
  - sample counts before and after filtering
  - any missing-value handling choices

The main purpose of these artifacts is to distinguish solver differences from data-preparation differences quickly.

## Functional Design

### `A.mat`

First-pass expectations:

- Input: aligned marker matrix with samples on rows and markers on columns
- Marker values: current dataset uses `0/1/2`
- Output: additive relationship matrix with sample labels preserved
- First milestone only covers the standard path needed for this dataset
- Advanced options can be added later behind explicit parameters

First-pass success criteria:

- Executes on the real aligned dataset
- Produces a square sample-by-sample matrix
- Matrix ordering matches the aligned sample order
- Matrix can be compared elementwise against the R baseline

### `mixed.solve`

First-pass expectations:

- Support the minimal common path needed by the project:
  - response vector `y`
  - fixed-effect design matrix `X`
  - random-effect design matrix `Z`
  - covariance matrix `K`
- Output at least:
  - fixed effects estimates
  - random effects estimates
  - variance component estimates
  - enough information to derive fitted values

First-pass success criteria:

- Can fit the model configuration required by the single-location pipeline
- Produces structured output suitable for direct comparison with the R baseline

### `kin.blup`

First-pass expectations:

- Accept phenotype input with sample column, value column, and selected location column
- Use genotype-derived `K` or marker-derived path needed for the first milestone
- Filter one location at a time
- Align samples explicitly
- Fit the model and return per-sample predictions

First-pass success criteria:

- Runs end-to-end on the provided real dataset for one location
- Produces predictions indexed by sample
- Produces outputs that can be compared directly with R baseline predictions

## Error Handling

The first milestone should prefer explicit failures over implicit guesses.

Required checks:

- Missing required columns in genotype or phenotype file
- Unknown location column requested by the user or script
- Empty sample intersection after filtering
- Illegal marker values outside the supported first-pass encoding
- Too few samples after filtering to fit the requested workflow safely

Error messages must say what failed and which columns or counts caused the failure.

## Testing Strategy

### 1. Unit Tests

Use small constructed inputs to verify:

- phenotype location filtering
- sample alignment behavior
- `A.mat` shape and symmetry
- basic mixed model output structure
- `kin.blup` output indexing

### 2. Real-Data Smoke Test

Use the provided genotype and phenotype files to verify:

- one location can be selected
- aligned data can be constructed
- the full Python pipeline runs without error
- result files can be written

### 3. R Baseline Regression Tests

Read baseline artifacts generated by the R script and compare:

- sample ordering
- phenotype vector shape
- `A.mat`
- variance components
- predictions

Tolerance should be staged:

1. Structure and shape match
2. Summary statistics and correlations are reasonable
3. Elementwise numeric agreement tightens over time

## Development Sequence

1. Initialize the repository and ignore worktree directories.
2. Write and commit this design spec.
3. Review and approve the spec.
4. Create a project-local worktree for implementation.
5. Scaffold the Python package, test layout, and script entry points.
6. Implement the R baseline script and generate baseline artifacts on one selected location.
7. Implement genotype and phenotype loading, location filtering, and sample alignment.
8. Implement Python `A.mat`.
9. Implement Python `mixed.solve`.
10. Implement Python `kin.blup` by combining the earlier pieces.
11. Add real-data smoke tests.
12. Add R baseline regression tests.
13. Iterate on numeric alignment with R.

## Worktree Plan

Implementation should happen in a project-local `.worktrees/` directory so the main repository root remains clean. The directory must be ignored by git before worktree creation.

## Open Decisions Resolved In This Spec

- Initial scope is limited to `A.mat`, `mixed.solve`, and `kin.blup`.
- Development uses one phenotype location at a time.
- Priority order is:
  1. make the pipeline run,
  2. establish tests and baselines,
  3. tighten numerical agreement with R.
- R `rrBLUP` outputs are the numerical reference for Python development.

## Out Of Scope For The First Milestone

- Broad reproduction of the full R package API
- Multi-location joint modeling
- Advanced options not required by the provided dataset
- Silent data correction heuristics

## Acceptance Criteria For The First Milestone

- A Python package exists in the `evo2` environment and can run on the provided CSV files.
- The pipeline can filter one location, align samples, and execute `A.mat`, `mixed.solve`, and `kin.blup`.
- An R baseline script exists and runs in the `rrblup` environment.
- Baseline artifacts are written in a stable format and reused by tests.
- Python tests cover both workflow execution and baseline comparison.
