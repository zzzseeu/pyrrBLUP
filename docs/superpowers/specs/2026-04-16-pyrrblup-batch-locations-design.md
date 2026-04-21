# pyrrBLUP Multi-Location Batch Design

## Goal

Extend the current single-location `rrBLUP` alignment workflow so that the project can run R baselines and Python comparisons across multiple phenotype locations in one batch, while preserving per-location isolation and producing a summary of success, failure, and numeric agreement.

## Scope

This design adds a batch orchestration layer on top of the existing single-location workflow.

The batch layer must:

- Discover all phenotype location columns matching `loc_*`
- Optionally filter to a user-provided subset of locations
- Run the existing R baseline export for each location
- Run the existing Python pipeline for each location
- Compare per-location R and Python outputs
- Continue processing after per-location failures
- Produce a per-location directory structure and a cross-location summary

The batch layer does not change the core mathematical goals of `A.mat`, `mixed.solve`, or `kin.blup`. It reuses those pieces.

## Confirmed Product Decisions

### Batch Selection

Default behavior will automatically discover every phenotype column with the `loc_*` prefix.

An optional location filter will allow running only a subset of locations.

### Output Structure

Results will be stored per location, with separate subdirectories for R and Python outputs:

- `runs/<location>/r/`
- `runs/<location>/python/`

### Failure Strategy

Batch execution will continue after a per-location failure.

The final batch result will include a success or failure summary across all attempted locations.

## Approaches Considered

### 1. Script Enhancement Only

Extend the existing single-location scripts directly to support looping over all locations.

Trade-offs:

- Smallest immediate code change
- Higher risk of script sprawl
- Harder to test and maintain once batch logic grows

### 2. Lightweight Orchestration Layer

Keep single-location scripts as reusable building blocks and add a separate batch orchestration layer for discovery, iteration, status tracking, comparison, and summary generation.

Trade-offs:

- Slightly more upfront structure
- Best separation of concerns
- Easiest to test and extend

### 3. Full CLI Framework

Refactor everything into a larger multi-command CLI system.

Trade-offs:

- Clean long-term interface
- Too much overhead for the current milestone

This design chooses option 2.

## Architecture

The batch workflow is split into five layers.

### 1. Location Discovery Layer

Responsibilities:

- Read phenotype columns
- Identify all `loc_*` columns
- Apply optional user-specified filtering
- Fail fast if requested locations do not exist or if filtering removes every location

Reasoning:

Location selection is independent from model execution and should be deterministic and easy to test.

### 2. Single-Location Execution Layer

Responsibilities:

- Reuse the current single-location R baseline script
- Reuse the current single-location Python pipeline
- Keep single-location behavior stable while batch orchestration is added around it

Reasoning:

The project already has a validated single-location core. The batch layer should call into that instead of reimplementing it.

### 3. Per-Location Artifact Layer

Responsibilities:

- Create a dedicated directory per location
- Store R outputs under `runs/<location>/r/`
- Store Python outputs under `runs/<location>/python/`
- Store per-tool status metadata in each output directory

Reasoning:

Per-location isolation makes reruns, debugging, and cross-location comparisons much simpler.

### 4. Comparison Layer

Responsibilities:

- Compare R and Python outputs per location
- Compute numeric metrics such as:
  - `amat_rmse`
  - `prediction_rmse`
  - variance component differences
- Record when comparison is not possible because one side failed

Reasoning:

Batch mode is only useful if it not only runs locations, but also tells us where Python still diverges from R.

### 5. Summary Layer

Responsibilities:

- Aggregate all per-location status and comparison metrics
- Write a machine-readable summary file
- Preserve partial success information when some locations fail

Reasoning:

The user needs a single place to inspect overall batch health and prioritize further alignment work.

## Repository Changes

The current algorithm modules stay focused on single-location math. New files will be added for orchestration.

### New Python Modules

- `src/pyrrblup/locations.py`
  - discover `loc_*` columns
  - apply optional location filtering
- `src/pyrrblup/batch.py`
  - main batch orchestration logic
  - iterate over locations
  - invoke per-location execution steps
  - collect comparison and status results
- `src/pyrrblup/compare.py`
  - compare R and Python artifacts
  - compute per-location numeric metrics
- `src/pyrrblup/status.py`
  - write and read `status.json`
  - define status serialization helpers

### New Script Entrypoints

- `scripts/run_batch_pipeline.py`
  - Python batch entrypoint
  - auto-discover locations by default
  - support optional location subset filtering
- `scripts/build_r_baseline_batch.R`
  - R batch entrypoint
  - reuse current single-location baseline generation logic

### Existing Files With Limited Changes

- `src/pyrrblup/models.py`
  - add dataclasses for batch run status and comparison summaries if needed
- `scripts/build_r_baseline.R`
  - keep single-location behavior
  - make it callable or reusable from the batch R script
- `scripts/run_python_pipeline.py`
  - keep single-location behavior
  - allow reuse from the batch Python entrypoint

## Output Layout

For each location, batch mode will create a directory shaped like:

```text
runs/
└── loc_BeiJ/
    ├── r/
    └── python/
```

### R Output Directory

`runs/<location>/r/` should contain:

- `amat.csv`
- `kin_blup_predictions.csv`
- `mixed_solve_beta.csv`
- `mixed_solve_u.csv`
- `mixed_solve_varcomp.csv`
- `samples.csv`
- `phenotype.csv`
- `location.txt`
- `status.json`

### Python Output Directory

`runs/<location>/python/` should contain:

- `amat.csv`
- `kin_blup_predictions.csv`
- `mixed_solve_varcomp.json`
- `status.json`

Python output can later be expanded to include additional files such as `mixed_solve_u.csv`, but this is not required for the first batch milestone.

## Status Metadata Design

Each tool-specific output directory should contain a `status.json` file with at least:

- `location`
- `tool`
- `success`
- `error_message`
- `started_at`
- `finished_at`

This allows failures to be inspected after batch completion without relying on terminal logs.

## Batch Summary Design

Batch mode should write:

- `runs/summary.csv`

The summary should include at least:

- `location`
- `r_success`
- `python_success`
- `sample_count`
- `amat_rmse`
- `prediction_rmse`
- `vu_r`
- `vu_py`
- `ve_r`
- `ve_py`
- `error_stage`
- `error_message`

This summary provides a single table for triaging failures and ranking locations by alignment quality.

## Batch Execution Flow

### 1. Discover Locations

- Read phenotype columns
- Identify all columns beginning with `loc_`
- Apply optional location filter
- Raise an explicit error if no locations remain

### 2. Execute Each Location

For each chosen location:

- Create `runs/<location>/r/`
- Create `runs/<location>/python/`
- Run the R baseline export into the R directory
- Run the Python pipeline into the Python directory
- Run a comparison step if both sides succeeded
- Write `status.json` for each side regardless of success or failure

### 3. Aggregate Results

- Collect per-location statuses and numeric comparison metrics
- Write `runs/summary.csv`

### 4. Exit Behavior

The batch process should continue after per-location failures.

If any location fails, the batch command should still finish the full run and return a non-zero exit code so callers know the batch is not fully clean.

## Error Handling

The batch system should prefer explicit status reporting over silent skips.

Required checks include:

- no `loc_*` columns found
- requested filtered locations not present in the phenotype file
- no locations remain after filtering
- per-location R execution failure
- per-location Python execution failure
- comparison failure because required artifacts are missing

Per-location errors should be written into status metadata and summary output.

## Testing Strategy

### 1. Unit Tests

Add focused tests for:

- location discovery
- location filtering
- summary row construction
- comparison metric calculation
- continue-on-failure behavior in orchestrator logic

### 2. Integration Tests With Small Synthetic Data

Create tests that simulate 2-3 locations and verify:

- output directories are created correctly
- `status.json` is written for each location/tool
- summary generation works
- one failing location does not stop the batch

### 3. Real-Data Smoke Tests

Start with a small subset of real locations, then expand to all discovered locations.

Verify:

- directory layout is correct
- summary output is created
- per-location comparisons are computed

## Success Criteria

The first batch milestone is complete when:

- the system can discover all `loc_*` columns automatically
- the user can optionally filter to a subset of locations
- a batch run creates `runs/<location>/r/` and `runs/<location>/python/`
- failures in one location do not stop other locations from running
- `runs/summary.csv` is produced
- tests cover location discovery, continue-on-failure behavior, and batch comparison output

## Out Of Scope For This Milestone

- redesigning the single-location core API
- large CLI framework refactors
- multi-location joint statistical modeling
- parallel batch execution optimization
- replacing the current per-location artifact formats
