# Generated results

The experiment scripts write all generated data and figures under this
directory and print the absolute path immediately after each file is saved.

All runs use one fixed directory per experiment:

```text
results/<experiment>/
├── data/          # CSV and JLD2 numerical outputs
├── figures/       # SVG and PDF figures
├── configs/       # Resolved configuration for each experiment
└── metadata.toml  # Added by the all-results runner
```

Running the same experiment again updates its generated files at the same
paths. Existing historical timestamped directories are not used by new runs.

| Experiment | Saved result |
|---|---|
| `01_field_reconstruction.jl` | Field-reconstruction JLD2 and `estimate_1` SVG/PDF |
| `02_clarity_vs_time.jl` | One clarity JLD2 per sensor count and `estimate_2` SVG/PDF |
| `03_discrete_continuous_bounds.jl` | Bound-comparison CSV and SVG/PDF |
| `04_grid_convergence.jl` | Grid-convergence CSV and SVG/PDF |
| `05_sensor_number_table.jl` | Minimum-sensor table CSV |
| `06_clarity_vs_sensors.jl` | Clarity-versus-sensors CSV and SVG/PDF |
| `07_noise_rate_tradeoff.jl` | Clarity-grid JLD2 and heatmap SVG/PDF |
