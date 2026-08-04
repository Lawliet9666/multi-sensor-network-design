function save_paper_plot(plot_object, output_directory, stem)
    mkpath(output_directory)
    svg_path = joinpath(output_directory, "$stem.svg")
    pdf_path = joinpath(output_directory, "$stem.pdf")
    savefig(plot_object, svg_path)
    println("Saved figure (SVG): $(abspath(svg_path))")
    savefig(plot_object, pdf_path)
    println("Saved figure (PDF): $(abspath(pdf_path))")
    return (svg=abspath(svg_path), pdf=abspath(pdf_path))
end

function plot_field_reconstruction(data::SpatiotemporalData2D, problem::STGPKFProblem,
                                   state::KFState, time_index::Int, output_directory)
    1 <= time_index <= length(data.ts) || throw(BoundsError(data.ts, time_index))
    ground_truth = data.data[:, :, time_index]'
    estimate = reshape(get_estimate(problem, state), length(data.xs), length(data.ys))'
    clarity = reshape(get_estimate_clarity(problem, state), length(data.xs), length(data.ys))'
    common = (aspect_ratio=:equal, framestyle=:box, colorbar=false)
    p1 = heatmap(data.xs, data.ys, ground_truth; common..., color=:bluesreds, title="Ground truth")
    p2 = heatmap(data.xs, data.ys, estimate; common..., color=:bluesreds, title="Estimate")
    p3 = heatmap(data.xs, data.ys, clarity; common..., clims=(0, 1), title="Clarity")
    combined = plot(p1, p2, p3; layout=(1, 3), size=(1050, 320), fontfamily="Times")
    return save_paper_plot(combined, output_directory, "estimate_1")
end

function plot_mean_clarity(times, curves::AbstractDict, output_directory)
    figure = plot(; xlabel="Time [min]", ylabel="Clarity", ylims=(0, 1),
                  legend=:bottomright, size=(620, 410), fontfamily="Times",
                  framestyle=:box)
    colors = palette(:plasma, length(curves))
    for ((label, values), color) in zip(sort(collect(curves); by=first), colors)
        plot!(figure, times, values; label=label, linewidth=2, color=color)
    end
    hline!(figure, [0.7]; label="Target clarity", linestyle=:dash, color=:gray)
    return save_paper_plot(figure, output_directory, "estimate_2")
end

function plot_bound_comparison(steps, continuous_values, discrete_values, output_directory)
    left = plot(steps, continuous_values; label="Continuous bound", marker=:circle,
                xlabel="Sampling interval", ylabel="Mean covariance", linewidth=2,
                fontfamily="Times", framestyle=:box)
    plot!(left, steps, discrete_values; label="Discrete bound", marker=:diamond, linewidth=2)
    relative_error = abs.((continuous_values .- discrete_values) ./ discrete_values)
    right = plot(steps, relative_error; label=false, marker=:circle,
                 xlabel="Sampling interval", ylabel="Relative error", linewidth=2,
                 fontfamily="Times", framestyle=:box)
    combined = plot(left, right; layout=(1, 2), size=(900, 360))
    return save_paper_plot(combined, output_directory, "continuous_vs_discrete2")
end

function plot_grid_convergence(rows, output_directory)
    figure = plot(; xlabel="Number of grid points", ylabel="Clarity lower bound", ylims=(0, 1),
                  fontfamily="Times", framestyle=:box)
    sensor_counts = sort(unique(row.N_robots for row in rows))
    for count in sensor_counts
        selected = sort(filter(row -> row.N_robots == count, rows); by=row -> row.Ng)
        plot!(figure, [row.Ng for row in selected], [row.clarity for row in selected];
              label="Sensors = $count", marker=:circle, linewidth=2)
    end
    return save_paper_plot(figure, output_directory, "ng_converge_clarity")
end

function plot_sensor_table_curve(rows, output_directory)
    ordered = sort(rows; by=row -> row.N_robots)
    counts = [row.N_robots for row in ordered]
    bound = [row.bound_clarity for row in ordered]
    figure = plot(counts, bound; label="Analytic lower bound", marker=:circle,
                  xlabel="Number of sensors", ylabel="Clarity", ylims=(0, 1), linewidth=2,
                  fontfamily="Times", framestyle=:box)
    if hasproperty(first(ordered), :empirical_clarity)
        empirical = [row.empirical_clarity for row in ordered]
        plot!(figure, counts, empirical; label="Empirical clarity", marker=:diamond, linewidth=2)
    end
    return save_paper_plot(figure, output_directory, "clarity_vs_nr")
end

function plot_noise_rate_tradeoff(noise_variances, steps, clarity, sensor_count, output_directory)
    figure = contourf(noise_variances, steps, clarity;
        xlabel="Measurement noise variance", ylabel="Sampling interval / sensors",
        color=:plasma, colorbar_title="Clarity lower bound",
        fontfamily="Times", framestyle=:box, size=(720, 520), linewidth=0)
    for intensity in (2.0, 5.0, 10.0)
        visible_x = [x for x in noise_variances if minimum(steps) <= sensor_count / (intensity * x) <= maximum(steps)]
        isempty(visible_x) && continue
        plot!(figure, visible_x, sensor_count ./ (intensity .* visible_x);
              color=:white, linestyle=:dash, linewidth=1.5, label=false)
    end
    return save_paper_plot(figure, output_directory, "clarity_heatmap_Nr1")
end
