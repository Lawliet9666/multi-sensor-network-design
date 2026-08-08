const PAPER_COLORS = (
    continuous=:orange,
    discrete=:steelblue,
    empirical=:steelblue,
    analytic=:orange,
    absolute_error=:crimson,
    relative_error=:purple,
    target=:gray,
)

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
    common = (
        aspect_ratio=:equal,
        framestyle=:box,
        colorbar=true,
        xlabel="x [km]",
        ylabel="y [km]",
        titlefontsize=12,
        tickfontsize=9,
        guidefontsize=11,
        colorbar_tickfontsize=9,
        left_margin=4Plots.mm,
    )
    p1 = heatmap(data.xs, data.ys, ground_truth;
                 common..., color=:bluesreds, clims=(-5, 5), title="Ground Truth")
    p2 = heatmap(data.xs, data.ys, estimate;
                 common..., color=:bluesreds, clims=(-5, 5), title="Estimated")
    p3 = heatmap(data.xs, data.ys, clarity;
                 common..., color=:inferno, clims=(0, 1), title="Clarity")
    combined = plot(p1, p2, p3; layout=(1, 3), size=(1050, 340), fontfamily="Times")
    return save_paper_plot(combined, output_directory, "estimate_1")
end

function plot_mean_clarity(times, curves::AbstractDict, output_directory)
    figure = plot(; xlabel="Time [min]", ylabel="Clarity", ylims=(0, 1),
                  legend=:bottomright, size=(600, 400), fontfamily="Times",
                  framestyle=:box, grid=true, gridalpha=0.3)
    ordered = sort(collect(curves); by=first)
    # Reuse the Plasma positions from the original paper figure.
    original_color_positions = Dict(
        1 => 0.4869846698753592,
        6 => 0.639367550096278,
        20 => 0.7423411802387613,
    )
    fallback_positions = range(0.2, 0.8; length=length(ordered))
    color_map = cgrad(:plasma)
    for (index, (sensor_count, values)) in enumerate(ordered)
        color_position = get(original_color_positions, sensor_count, fallback_positions[index])
        color = get(color_map, color_position)
        label = sensor_count == 1 ? "1 sensor" : "$sensor_count sensors"
        plot!(figure, times, values; label=label, linewidth=2, color=color)
    end
    hline!(figure, [0.7]; label="Target clarity", linestyle=:dash,
           linewidth=2, color=PAPER_COLORS.target)
    return save_paper_plot(figure, output_directory, "estimate_2")
end

function plot_bound_comparison(steps, continuous_values, discrete_values, output_directory)
    maximum_bound = maximum(vcat(continuous_values, discrete_values))
    bound_limit = maximum_bound == 0 ? 1.0 : 1.1 * maximum_bound
    bounds = plot(steps, continuous_values;
        label="Continuous bound", marker=:circle, markersize=7,
        xlabel="Sampling interval", ylabel="Mean covariance", linewidth=3,
        color=PAPER_COLORS.continuous, legend=:topleft, gridalpha=0.3,
        fontfamily="Times", framestyle=:box, ylims=(0, bound_limit),
        left_margin=4Plots.mm)
    plot!(bounds, steps, discrete_values;
        label="Discrete bound", marker=:diamond, markersize=7, linewidth=3,
        color=PAPER_COLORS.discrete)

    absolute_error = abs.(continuous_values .- discrete_values)
    relative_error = 100 .* absolute_error ./ abs.(discrete_values)
    maximum_absolute_error = maximum(absolute_error)
    absolute_error_limit = maximum_absolute_error == 0 ? 1.0 : 1.1 * maximum_absolute_error
    errors = plot(steps, absolute_error;
        label="Absolute error", marker=:square, markersize=7,
        xlabel="Sampling interval", ylabel="Absolute error", linewidth=3,
        color=PAPER_COLORS.absolute_error, legend=:topleft, gridalpha=0.3,
        fontfamily="Times", framestyle=:box, ylims=(0, absolute_error_limit),
        left_margin=4Plots.mm)
    plot!(errors, [NaN], [NaN];
        label="Relative error (%)", marker=:star, markersize=10, linewidth=3,
        color=PAPER_COLORS.relative_error)
    relative_axis = twinx(errors)
    plot!(relative_axis, steps, relative_error;
        label=false, marker=:star, markersize=10,
        ylabel="Relative error (%)", linewidth=3,
        color=PAPER_COLORS.relative_error, grid=false, ylims=(0, 100),
        right_margin=4Plots.mm)

    combined = plot(bounds, errors; layout=(2, 1), size=(650, 760))
    return save_paper_plot(combined, output_directory, "continuous_vs_discrete2")
end

function plot_grid_convergence(rows, output_directory)
    sensor_counts = sort(unique(row.N_robots for row in rows))
    colors = palette(:plasma, max(length(sensor_counts), 7))
    figure = plot(; xlabel="Number of grid points", ylabel="Clarity lower bound",
                  legend=:bottomright, gridalpha=0.3, size=(620, 460),
                  fontfamily="Times", framestyle=:box)
    for (index, count) in enumerate(sensor_counts)
        selected = sort(filter(row -> row.N_robots == count, rows); by=row -> row.Ng)
        plot!(figure, [row.Ng for row in selected], [row.clarity for row in selected];
              label="Sensors = $count", marker=:circle, markersize=5, linewidth=2,
              color=colors[index])
    end
    return save_paper_plot(figure, output_directory, "ng_converge_clarity")
end

function plot_sensor_table_curve(rows, output_directory)
    ordered = sort(rows; by=row -> row.N_robots)
    counts = [row.N_robots for row in ordered]
    bound = [row.bound_clarity for row in ordered]
    figure = plot(; xlabel="Number of sensors", ylabel="Clarity", ylims=(0, 0.99),
                  legend=:bottomright, gridalpha=0.3,
                  fontfamily="Times", framestyle=:box)
    if hasproperty(first(ordered), :empirical_clarity)
        empirical = [row.empirical_clarity for row in ordered]
        plot!(figure, counts, empirical;
              label="Empirical clarity", marker=:diamond, markersize=8, linewidth=2,
              color=PAPER_COLORS.empirical)
    end
    plot!(figure, counts, bound;
          label="Analytic lower bound", marker=:square, markersize=5, linewidth=2,
          color=PAPER_COLORS.analytic)
    return save_paper_plot(figure, output_directory, "clarity_vs_nr")
end

function plot_noise_rate_tradeoff(
    noise_variances, steps, clarity, sensor_count, output_directory;
    fix_sigma_c=false,
)
    noise_label = fix_sigma_c ?
        "Continuous-time noise intensity" : "Measurement noise variance"
    figure = contourf(noise_variances, steps, clarity;
        xlabel=noise_label, ylabel="Sampling interval",
        color=:plasma, colorbar_title="Clarity lower bound",
        fontfamily="Times", framestyle=:box, size=(800, 600), linewidth=0)
    for intensity in (2.0, 5.0, 10.0)
        if fix_sigma_c
            noise_intensity = sensor_count / intensity
            minimum(noise_variances) <= noise_intensity <= maximum(noise_variances) || continue
            vline!(figure, [noise_intensity];
                   color=:white, linestyle=:dash, linewidth=2, label=false)
        else
            visible_x = [
                x for x in noise_variances
                if minimum(steps) <= sensor_count / (intensity * x) <= maximum(steps)
            ]
            isempty(visible_x) && continue
            plot!(figure, visible_x, sensor_count ./ (intensity .* visible_x);
                  color=:white, linestyle=:dash, linewidth=2, label=false)
        end
    end
    return save_paper_plot(figure, output_directory, "clarity_heatmap_Nr1")
end
