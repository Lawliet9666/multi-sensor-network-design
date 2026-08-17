const PAPER_COLORS = (
    continuous=:orange,
    discrete=:steelblue,
    empirical=:darkgreen,
    analytic=:orange,
    absolute_error=:crimson,
    relative_error=:purple,
    target=:gray,
)

function save_paper_plot(plot_object, output_directory, stem)
    mkpath(output_directory)
    pdf_path = joinpath(output_directory, "$stem.pdf")
    savefig(plot_object, pdf_path)
    println("Saved figure (PDF): $(abspath(pdf_path))")
    return abspath(pdf_path)
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
        top_margin=2Plots.mm,
    )
    title_x = (first(data.xs) + last(data.xs)) / 2
    title_y = last(data.ys) + 0.105 * (last(data.ys) - first(data.ys))
    p1 = heatmap(data.xs, data.ys, ground_truth;
                 common..., color=:bluesreds, clims=(-5, 5))
    annotate!(p1, title_x, title_y, text("Ground Truth", 12, :black, :center, "Times"))
    p2 = heatmap(data.xs, data.ys, estimate;
                 common..., color=:bluesreds, clims=(-5, 5))
    annotate!(p2, title_x, title_y, text("Estimated", 12, :black, :center, "Times"))
    p3 = heatmap(data.xs, data.ys, clarity;
                 common..., color=:inferno, clims=(0, 1))
    annotate!(p3, title_x, title_y, text("Clarity", 12, :black, :center, "Times"))
    combined = plot(p1, p2, p3; layout=(1, 3), size=(1050, 340), fontfamily="Times")
    return save_paper_plot(combined, output_directory, "estimate_1")
end

function plot_mean_clarity(times, curves::AbstractDict, output_directory)
    figure = plot(; xlabel="Time [min]", ylabel="Clarity", ylims=(0, 1),
                  legend=:bottomright, size=(450, 300), fontfamily="Times",
                  framestyle=:box, grid=false,
                  guidefontsize=12, tickfontsize=10, legendfontsize=10,
                  left_margin=6Plots.mm, bottom_margin=6Plots.mm)
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
        label = sensor_count == 1 ? "1 Agent" : "$sensor_count Agents"
        plot!(figure, times, values; label=label, linewidth=2, color=color)
    end
    hline!(figure, [0.7]; label="Target clarity", linestyle=:dash,
           linewidth=2, color=PAPER_COLORS.target)
    return save_paper_plot(figure, output_directory, "estimate_2")
end

function plot_bound_comparison(trajectory, error_rows, output_directory)
    maximum_points = 500
    display_stride = max(1, cld(length(trajectory.times), maximum_points))
    display_indices = collect(1:display_stride:length(trajectory.times))
    last(display_indices) == length(trajectory.times) || push!(
        display_indices, length(trajectory.times),
    )
    display_times = trajectory.times[display_indices]
    empirical_values = trajectory.empirical_mean_covariance[display_indices]
    discrete_values = trajectory.discrete_bound_mean_covariance[display_indices]
    continuous_values = trajectory.continuous_bound_mean_covariance[display_indices]

    maximum_bound = maximum(vcat(empirical_values, discrete_values, continuous_values))
    bound_limit = maximum_bound == 0 ? 1.0 : 1.1 * maximum_bound
    bounds = plot(display_times, empirical_values;
        label=L"\mathcal{L}(\widehat{\mathrm{E}}[\Sigma_k])", xlabel="Time [min]",
        ylabel=L"\mathcal{L}(\cdot)", linewidth=3, color=PAPER_COLORS.empirical,
        legend=:topright, grid=false,
        fontfamily="Times", framestyle=:box, ylims=(0, bound_limit),
        guidefontsize=22, tickfontsize=18, legendfontsize=20,
        left_margin=9Plots.mm, bottom_margin=10Plots.mm)
    plot!(bounds, display_times, discrete_values;
        label=L"\mathcal{L}(\Delta_k)", linewidth=3, linestyle=:dash,
        color=PAPER_COLORS.discrete)
    plot!(bounds, display_times, continuous_values;
        label=L"\mathcal{L}(\Delta(t_k))", linewidth=3, linestyle=:dot,
        color=PAPER_COLORS.continuous)

    ordered_rows = sort(collect(error_rows); by=row -> row.dt)
    steps = getproperty.(ordered_rows, :dt)
    empirical_errors = getproperty.(
        ordered_rows, :empirical_max_relative_linear_operator_error_percent,
    )
    discrete_errors = getproperty.(
        ordered_rows, :discrete_max_relative_linear_operator_error_percent,
    )
    errors = plot(steps, empirical_errors;
        label=L"e_{\mathrm{exp}}^T(\Delta t)",
        marker=:circle, markersize=7,
        xlabel="Sampling interval [min]",
        ylabel=L"e^T(\Delta t)\ [\%]",
        linewidth=3, color=PAPER_COLORS.empirical, legend=:topleft, grid=false,
        fontfamily="Times", framestyle=:box, ylims=(0, 100), yticks=0:20:100,
        guidefontsize=22, tickfontsize=18, legendfontsize=20,
        left_margin=9Plots.mm, bottom_margin=10Plots.mm)
    plot!(errors, steps, discrete_errors;
        label=L"e_{\mathrm{bound}}^T(\Delta t)",
        marker=:diamond, markersize=7,
        linewidth=3, color=PAPER_COLORS.discrete)

    combined = plot(bounds, errors; layout=(1, 2), size=(1100, 460))
    return save_paper_plot(combined, output_directory, "continuous_vs_discrete2")
end

function plot_grid_convergence(rows, output_directory)
    sensor_counts = sort(unique(row.N_robots for row in rows))
    colors = palette(:plasma, max(length(sensor_counts), 7))
    figure = plot(; xlabel="Number of grid points", ylabel="Clarity lower bound",
                  legend=:bottomright, grid=false, size=(620, 350),
                  fontfamily="Times", framestyle=:box,
                  guidefontsize=15, tickfontsize=12, legendfontsize=11,
                  left_margin=7Plots.mm, bottom_margin=7Plots.mm)
    for (index, count) in enumerate(sensor_counts)
        selected = sort(filter(row -> row.N_robots == count, rows); by=row -> row.Ng)
        plot!(figure, [row.Ng for row in selected], [row.clarity for row in selected];
              label="Agents = $count", linewidth=2,
              color=colors[index])
    end
    return save_paper_plot(figure, output_directory, "ng_converge_clarity")
end

function plot_sensor_table_curve(
    rows,
    output_directory;
    target_clarity=nothing,
)
    ordered = sort(rows; by=row -> row.N_robots)
    counts = [row.N_robots for row in ordered]
    bound = [row.bound_clarity for row in ordered]
    all(row -> hasproperty(row, :expected_spatial_mean_clarity_estimate), ordered) ||
        error("Sensor-clarity rows must contain expected spatial-mean clarity estimates.")
    monte_carlo_trials = unique(row.monte_carlo_trials for row in ordered)
    length(monte_carlo_trials) == 1 || error(
        "Sensor-clarity rows must use a common Monte Carlo trial count.",
    )
    trial_count = only(monte_carlo_trials)
    isinteger(trial_count) || error("Monte Carlo trial count must be an integer.")
    expected_clarity = [row.expected_spatial_mean_clarity_estimate for row in ordered]
    label_expected = L"\widehat{\mathrm{E}}[\bar q]"
    common = (
        xlabel=L"N_r",
        ylabel="Clarity",
        xlims=(0, 52),
        ylims=(0, 0.95),
        xticks=0:10:50,
        yticks=0:0.2:0.8,
        legend=:bottomright,
        grid=false,
        size=(760, 560),
        fontfamily="Times",
        framestyle=:box,
        guidefontsize=36,
        tickfontsize=28,
        legendfontsize=23,
        left_margin=10Plots.mm,
        bottom_margin=10Plots.mm,
    )
    figure = plot(; common...)
    plot!(figure, counts, expected_clarity;
          label=label_expected, marker=:diamond,
          markersize=10, linewidth=4,
          color=:steelblue)
    plot!(figure, counts, bound;
          label=L"\bar q_{\Delta^\Pi_\infty}", marker=:square,
          markersize=9, linewidth=4,
          color=PAPER_COLORS.analytic)
    if !isnothing(target_clarity)
        target_index = findfirst(>=(target_clarity), bound)
        isnothing(target_index) && error(
            "No configured sensor count reaches target clarity $target_clarity.",
        )
        target_count = counts[target_index]
        plot!(figure, [0, target_count], [target_clarity, target_clarity];
              label=false, color=:red, linestyle=:dash,
              linewidth=3)
        plot!(figure, [target_count, target_count], [0, target_clarity];
              label=false, color=:red, linestyle=:dash,
              linewidth=3)
    end
    return save_paper_plot(figure, output_directory, "clarity_vs_nr_compact")
end

function plot_noise_rate_tradeoff(
    measurement_variances,
    normalized_intervals,
    clarity,
    output_directory;
    contour_targets=(0.4, 0.5, 0.6),
)
    size(clarity) == (length(normalized_intervals), length(measurement_variances)) ||
        throw(DimensionMismatch(
            "Clarity grid must have one row per normalized interval and one column per measurement variance.",
        ))
    all(>(0), measurement_variances) || throw(ArgumentError(
        "Measurement variances must be positive.",
    ))
    all(>(0), normalized_intervals) || throw(ArgumentError(
        "Normalized sampling intervals must be positive.",
    ))

    levels = range(minimum(clarity), maximum(clarity); length=13)
    figure = contourf(measurement_variances, normalized_intervals, clarity;
        levels=levels,
        xlabel=L"\sigma_m^2",
        ylabel=L"\Delta t/N_r",
        title=L"\bar q_{\Delta^\Pi_\infty}",
        color=:plasma,
        colorbar_ticks=[0.3, 0.5, 0.7, 0.9],
        fontfamily="Times",
        framestyle=:box,
        size=(760, 560),
        linewidth=0,
        xlims=extrema(measurement_variances),
        ylims=extrema(normalized_intervals),
        xticks=[0.5, 2.0, 3.5],
        yticks=[0.5, 1.0, 1.5, 2.0],
        guidefontsize=36,
        tickfontsize=28,
        titlefontsize=32,
        colorbar_tickfontsize=28,
        left_margin=10Plots.mm,
        bottom_margin=10Plots.mm,
        top_margin=4Plots.mm,
        right_margin=12Plots.mm,
    )

    x_min, x_max = extrema(measurement_variances)
    y_min, y_max = extrema(normalized_intervals)
    for target in contour_targets
        index = argmin(abs.(clarity .- target))
        row, column = Tuple(index)
        theta = inv(measurement_variances[column] * normalized_intervals[row])
        visible_min = max(x_min, inv(theta * y_max))
        visible_max = min(x_max, inv(theta * y_min))
        visible_min < visible_max || continue
        curve_x = range(visible_min, visible_max; length=300)
        curve_y = inv.(theta .* curve_x)
        plot!(figure, curve_x, curve_y;
              color=:white, linestyle=:dash, linewidth=3, label=false)

        label_x = visible_min + 0.7 * (visible_max - visible_min)
        label_y = inv(theta * label_x) + 0.08 * (y_max - y_min)
        annotate!(figure, label_x, label_y, text(
            latexstring("\\bar{q}=", round(target; digits=1)),
            24,
            :white,
        ))
    end
    return save_paper_plot(figure, output_directory, "clarity_heatmap_Nr1_compact")
end
