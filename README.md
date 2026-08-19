# STGPKF Network Design

Paper "Kalman-Bucy Filtering with Randomized Sensing: Fundamental Limits and Sensor Network Design for Field Estimation."

## Prerequisites

- Julia 1.11

## Getting Started

Clone the repository, enter it, and instantiate the Julia environment:

```bash
git clone https://github.com/Lawliet9666/multi-sensor-network-design.git
cd multi-sensor-network-design
julia --startup-file=no --project=. -e 'using Pkg; Pkg.instantiate()'
```

## Quick Start

```bash
julia --startup-file=no --project=. \
  experiments/03_discrete_continuous_bounds.jl \
  --output results/03_discrete_continuous_bounds
```

### Experiment Outputs


| Experiment                                                                         | Result                                           | Reference output                                                                                            |
| ---------------------------------------------------------------------------------- | ------------------------------------------------ | ----------------------------------------------------------------------------------------------------------- |
| `[01_field_reconstruction.jl](experiments/01_field_reconstruction.jl)`             | Field reconstruction and clarity map             | `[estimate_1.pdf](results/01_field_reconstruction/figures/estimate_1.pdf)`                                  |
| `[02_clarity_vs_time.jl](experiments/02_clarity_vs_time.jl)`                       | Mean clarity over time for 1, 6, and 20 sensors  | `[estimate_2.pdf](results/02_clarity_vs_time/figures/estimate_2.pdf)`                                       |
| `[03_discrete_continuous_bounds.jl](experiments/03_discrete_continuous_bounds.jl)` | Finite-horizon covariance and bound comparison   | `[continuous_vs_discrete2.pdf](results/03_discrete_continuous_bounds/figures/continuous_vs_discrete2.pdf)`  |
| `[04_grid_convergence.jl](experiments/04_grid_convergence.jl)`                     | Spatial-grid convergence study                   | `[ng_converge_clarity.pdf](results/04_grid_convergence/figures/ng_converge_clarity.pdf)`                    |
| `[05_sensor_number_table.jl](experiments/05_sensor_number_table.jl)`               | Minimum sensor count for each target clarity     | `[sensor_number_table.csv](results/05_sensor_number_table/data/sensor_number_table.csv)`                    |
| `[06_clarity_vs_sensors.jl](experiments/06_clarity_vs_sensors.jl)`                 | Expected clarity versus sensor count             | `[clarity_vs_nr_compact.pdf](results/06_clarity_vs_sensors/figures/clarity_vs_nr_compact.pdf)`              |
| `[07_noise_rate_tradeoff.jl](experiments/07_noise_rate_tradeoff.jl)`               | Measurement-noise and sampling-interval tradeoff | `[clarity_heatmap_Nr1_compact.pdf](results/07_noise_rate_tradeoff/figures/clarity_heatmap_Nr1_compact.pdf)` |


## Citation

If you use this code, please cite the corresponding paper:

```bibtex
@misc{wang2025kalmanbucy,
  title         = {Kalman-Bucy Filtering with Randomized Sensing: Fundamental Limits and Sensor Network Design for Field Estimation},
  author        = {Xinyi Wang and Devansh R. Agrawal and Dimitra Panagou},
  year          = {2025},
  eprint        = {2511.03740},
  archivePrefix = {arXiv},
  primaryClass  = {eess.SY},
  doi           = {10.48550/arXiv.2511.03740},
  url           = {https://arxiv.org/abs/2511.03740}
}
```

## Project Structure

```text
.
├── config/                 # Shared and experiment-specific TOML parameters
├── experiments/            # Seven supported paper experiment entry points
├── results/                # Committed reference data, PDFs, and configurations
├── scripts/
│   └── main.jl             # Runs and validates all seven experiments
├── src/                    # STGPKF, simulation, plotting, and bound code
├── test/
│   └── runtests.jl         # Core and paper-regression tests
├── Manifest-v1.11.toml     # Julia 1.11 dependency manifest
└── Project.toml            # Package metadata and compatibility constraints
```
