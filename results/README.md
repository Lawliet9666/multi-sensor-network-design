# Paper results

This directory contains the latest result snapshot used by the paper. The
experiment scripts update one fixed directory per experiment and print every
saved file with its absolute path.

```text
results/<experiment>/
├── data/       # CSV or JLD2 numerical result
├── figures/    # Paper PDF, when the experiment produces a figure
└── configs/    # Fully resolved paper configuration
```

Running the same experiment again updates the same files. Historical runs,
SVG copies, and figures not used by the paper are intentionally excluded.

| Experiment | Saved result |
|---|---|
| `01_field_reconstruction.jl` | Field-reconstruction JLD2 and `estimate_1.pdf` |
| `02_clarity_vs_time.jl` | One clarity JLD2 per agent count and `estimate_2.pdf` |
| `03_discrete_continuous_bounds.jl` | Finite-horizon error/trajectory CSVs and `continuous_vs_discrete2.pdf` |
| `04_grid_convergence.jl` | Grid-convergence CSV and `ng_converge_clarity.pdf` |
| `05_sensor_number_table.jl` | Minimum-sensor table CSV |
| `06_clarity_vs_sensors.jl` | Expected-clarity CSV and `clarity_vs_nr_compact.pdf` |
| `07_noise_rate_tradeoff.jl` | Tradeoff-grid JLD2 and `clarity_heatmap_Nr1_compact.pdf` |
