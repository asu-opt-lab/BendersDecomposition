# BendersX Paper Experiments

This is the standalone paper experiment harness. It is separate from the
existing `experiments/` test-style scripts so the paper runs have a stable
matrix, raw result files, and analysis scripts.

## Layout

- `Project.toml`: dedicated Julia environment.
- `scripts/common.jl`: shared solver factories, instance lists, model builders,
  SCFLP generator, baseline MIP models, and CSV logging.
- `scripts/01_correctness.jl`: UFLP/CFLP `p1:p71` correctness against freshly
  solved extensive-form MIP baselines.
- `scripts/02_oracle_ablation.jl`: fixed `BendersSeq`, varying oracle.
- `scripts/03_environment_ablation.jl`: fixed oracle, varying environment.
- `scripts/04_parallel_scflp.jl`: SCFLP `100x200` separable-oracle scaling.
- `results/raw/`: append-only run summaries and iteration traces.
- `results/processed/`: table-ready summaries from `analysis/`.
- `analysis/*.jl`: scripts that summarize raw results.

## Fixed Instance Plan

- Correctness:
  - UFLP: `p1:p71`
  - CFLP: `p1:p71`
- Runtime / ablation:
  - UFLP: `ga250a-1:5`, `ga250b-1:5`
  - CFLP: `T200x200_5_1:5`, `T200x200_10_1:5`
  - SCFLP: generated `100 facilities x 200 customers`, 15 instances:
    `f100-c200-s256-r5-1:5`, `f100-c200-s512-r5-1:5`, and
    `f100-c200-s1024-r5-1:5`.

The packaged artifacts do not currently contain `f100-c200-*` SCFLP files, so
`scripts/common.jl` generates deterministic `SCFLPData` instances using the
original stochastic CFLP generator logic. The current `SCFLPData` type does not
store scenario probabilities, so the uniform `prob_scenarios` vector is not
written into the data object.

## Commands

Instantiate:

```bash
julia --project=paper_experiments -e 'using Pkg; Pkg.instantiate()'
```

Run experiments:

```bash
julia --project=paper_experiments paper_experiments/scripts/01_correctness.jl
julia --project=paper_experiments paper_experiments/scripts/02_oracle_ablation.jl
julia --project=paper_experiments paper_experiments/scripts/03_environment_ablation.jl
julia --threads=1 --project=paper_experiments paper_experiments/scripts/04_parallel_scflp.jl
julia --threads=2 --project=paper_experiments paper_experiments/scripts/04_parallel_scflp.jl
julia --threads=4 --project=paper_experiments paper_experiments/scripts/04_parallel_scflp.jl
julia --threads=8 --project=paper_experiments paper_experiments/scripts/04_parallel_scflp.jl
```

Or use the shell wrappers:

```bash
paper_experiments/scripts/run_01_correctness.sh
paper_experiments/scripts/run_02_oracle_ablation.sh
paper_experiments/scripts/run_03_environment_ablation.sh
paper_experiments/scripts/run_04_parallel_scflp.sh
```

The wrappers accept environment-variable overrides, for example:

```bash
TIME_LIMIT=60 REPEATS=1 paper_experiments/scripts/run_02_oracle_ablation.sh
THREADS_LIST="1 4" REPEATS=1 paper_experiments/scripts/run_04_parallel_scflp.sh
```

## Slurm / Supercomputer Submission

The `hpc/` folder follows the existing repository pattern: each submit script
sets experiment metadata, copies the Julia script into the output directory,
generates one temporary job script per case, submits it with `sbatch`, and then
removes the temporary job script.

```bash
paper_experiments/hpc/submit_01_correctness.sh
paper_experiments/hpc/submit_02_oracle_ablation.sh
paper_experiments/hpc/submit_03_environment_ablation.sh
paper_experiments/hpc/submit_04_parallel_scflp.sh
```

Raw CSVs are written under nested directories in `results/raw/<round>/<version>/`.
The analysis scripts recursively merge those files.

Generate processed tables:

```bash
julia --project=paper_experiments paper_experiments/analysis/summarize_correctness.jl
julia --project=paper_experiments paper_experiments/analysis/summarize_oracle_ablation.jl
julia --project=paper_experiments paper_experiments/analysis/summarize_environment_ablation.jl
julia --project=paper_experiments paper_experiments/analysis/summarize_parallel_scflp.jl
```
