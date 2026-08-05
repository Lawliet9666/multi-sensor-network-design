# STGPKF Network Design

Julia code for the paper “Kalman-Bucy Filtering with Randomized Sensing: Fundamental Limits and Sensor Network Design for Field Estimation.” 

## Setup

Install Julia 1.11, then run:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

## Run

Generate all paper results:

```bash
julia --project=. scripts/reproduce_paper.jl --profile paper
```

This runs all seven experiments. Data, figures, resolved configurations, and
run metadata are saved under `results/paper-<timestamp>/`. Every saved file is printed with
its absolute path.

Run one result only:

```bash
julia --project=. experiments/03_discrete_continuous_bounds.jl --profile paper
```

Single-experiment outputs are saved under
`results/<experiment>/paper-<timestamp>/`.


| Script                             | Result                                      |
| ---------------------------------- | ------------------------------------------- |
| `01_field_reconstruction.jl`       | Field reconstruction and clarity map        |
| `02_clarity_vs_time.jl`            | Mean clarity over time                      |
| `03_discrete_continuous_bounds.jl` | Discrete and continuous covariance bounds   |
| `04_grid_convergence.jl`           | Grid-convergence study                      |
| `05_sensor_number_table.jl`        | Minimum sensor-number table                 |
| `06_clarity_vs_sensors.jl`         | Clarity versus sensor count                 |
| `07_noise_rate_tradeoff.jl`        | Measurement-noise and sensing-rate tradeoff |


Shared parameters are in `config/common.toml`; each script's `paper` parameters
are in the matching file under `config/experiments/`. Output details are in
`results/README.md`.

## Repository

```text
config/       Shared and experiment-specific parameters
experiments/  Seven result scripts
scripts/      All-results runner
src/          STGPKF and covariance-bound implementation
results/      Generated data, figures, and metadata
```

Citation metadata is in `CITATION.cff`. The code is released under the MIT
license; see `LICENSE` and `NOTICE.md`.