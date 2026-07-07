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
- `scripts/02_oracle_ablation.jl`: UFLP/CFLP `p1:p71`, fixed `BendersSeq`,
  varying oracle.
- `scripts/03_environment_ablation.jl`: UFLP only, fixed oracle, varying
  environment across `seq`, `seq_inout`, and `callback`.
- `scripts/04_parallel_scflp.jl`: callback-based SCFLP `50x100`
  separable-oracle scaling.
- `scripts/05_parallel_scflp_lp_relax.jl`: LP-relaxed SCFLP `50x50x512`
  separable-oracle thread scaling with `BendersSeq`.
- `results/raw/`: append-only run summaries and iteration traces.
- `results/processed/`: table-ready summaries from `analysis/`.
- `analysis/*.jl`: scripts that summarize raw results.

## Fixed Instance Plan

- Correctness:
  - UFLP: `p1:p71`
  - CFLP: `p1:p71`
- Oracle ablation:
  - UFLP: `p1:p71`
  - CFLP: `p1:p71`
- Environment ablation:
  - UFLP: `ga250c-1:5`, `gs250c-1:5`
- Parallel scaling:
  - SCFLP: generated `50 facilities x 100 customers`, 15 instances:
    `f50-c100-s128-r5-1:5`, `f50-c100-s256-r5-1:5`, and
    `f50-c100-s512-r5-1:5`.
- LP relaxation parallel scaling:
  - SCFLP: generated `50 facilities x 50 customers x 512 scenarios`, 5
    instances: `f50-c50-s512-r5-1:5`.

The packaged artifacts do not currently contain `f50-c100-*` SCFLP files, so
`scripts/common.jl` generates deterministic `SCFLPData` instances using the
original stochastic CFLP generator logic. The current `SCFLPData` type does not
store scenario probabilities, so the uniform `prob_scenarios` vector is not
written into the data object.

## Commands

Instantiate:

```bash
julia --project=paper_experiments -e 'using Pkg; Pkg.instantiate()'
```

The paper harness defaults to `--solver gurobi`. Use `--solver cplex` on the
Julia scripts, or `SOLVER=cplex` with the shell wrappers, to switch back to
CPLEX.

Run experiments:

```bash
julia --project=paper_experiments paper_experiments/scripts/01_correctness.jl
julia --project=paper_experiments paper_experiments/scripts/02_oracle_ablation.jl
julia --project=paper_experiments paper_experiments/scripts/03_environment_ablation.jl
julia --threads=1 --project=paper_experiments paper_experiments/scripts/04_parallel_scflp.jl
julia --threads=2 --project=paper_experiments paper_experiments/scripts/04_parallel_scflp.jl
julia --threads=4 --project=paper_experiments paper_experiments/scripts/04_parallel_scflp.jl
julia --threads=8 --project=paper_experiments paper_experiments/scripts/04_parallel_scflp.jl
julia --threads=16 --project=paper_experiments paper_experiments/scripts/04_parallel_scflp.jl
julia --threads=8 --project=paper_experiments paper_experiments/scripts/05_parallel_scflp_lp_relax.jl
```

`04_parallel_scflp.jl` defaults to `--env callback`, which uses `BendersBnB`.
Use `--env seq` to reproduce the old sequential environment run.
`05_parallel_scflp_lp_relax.jl` defaults to `--env seq`, generated
`f50-c50-s512-r5-1:5` instances, `--oracle cfl_knapsack`, and relaxes the
master integrality before solving.

Or use the shell wrappers:

```bash
paper_experiments/scripts/run_01_correctness.sh
paper_experiments/scripts/run_02_oracle_ablation.sh
paper_experiments/scripts/run_03_environment_ablation.sh
paper_experiments/scripts/run_04_parallel_scflp.sh
paper_experiments/scripts/run_05_parallel_scflp_lp_relax.sh
```

The wrappers accept environment-variable overrides, for example:

```bash
TIME_LIMIT=60 REPEATS=1 paper_experiments/scripts/run_02_oracle_ablation.sh
SOLVER=cplex TIME_LIMIT=60 REPEATS=1 paper_experiments/scripts/run_02_oracle_ablation.sh
THREADS_LIST="1 4" REPEATS=1 paper_experiments/scripts/run_04_parallel_scflp.sh
ENV_NAME=seq THREADS_LIST="1 4" REPEATS=1 paper_experiments/scripts/run_04_parallel_scflp.sh
THREADS_LIST="1 2 4 8 16" REPEATS=1 paper_experiments/scripts/run_05_parallel_scflp_lp_relax.sh
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
paper_experiments/hpc/submit_05_parallel_scflp_lp_relax.sh
```

Raw CSVs are written under nested directories in `results/raw/<round>/<version>/`.
The analysis scripts recursively merge those files.

Generate processed tables:

```bash
julia --project=paper_experiments paper_experiments/analysis/summarize_correctness.jl
julia --project=paper_experiments paper_experiments/analysis/summarize_oracle_ablation.jl
julia --project=paper_experiments paper_experiments/analysis/summarize_environment_ablation.jl
julia --project=paper_experiments paper_experiments/analysis/summarize_parallel_scflp.jl
julia --project=paper_experiments paper_experiments/analysis/summarize_parallel_scflp_lp_relax.jl
```
