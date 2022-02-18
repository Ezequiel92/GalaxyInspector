####################################################################################################
# Post-processing functions
####################################################################################################
# Every function is expected to have the following signature:
#
#   foo!(figure, args...; kwargs...) -> ([marker_object, ...], [label, ...]) or nothing
#
# Where:
#
#   - figure::Makie.Figure
#   - marker_object::LegendElement
#   - label::String
#
# So, it must take a figure, add something to it, and return how to label the additions to the plot. 
# It should return nothing if there are no new labels to add.
####################################################################################################

"""
    verticalFlags!(
        figure::Makie.Figure,
        flags::Tuple{Vector{<:Real},Vector{<:AbstractString}}, 
    )::Union{Tuple{Vector{<:LegendElement},Vector{String}},Nothing}

Draw vertical lines at specified positions.

# Arguments
- `figure::Makie.Figure`: Figure to be modified.
- `flags::Tuple{Vector{<:Real},Vector{<:AbstractString}}`: The first vector has the positions of
  the vertical lines. The second has the corresponding labels.

# Returns
- Tuple with two vectors:
  - `LineElement` to be used as markers in the legend.
  - Labels.

"""
function verticalFlags!(
    figure::Makie.Figure,
    flags::Tuple{Vector{<:Real},Vector{<:AbstractString}},
)::Union{Tuple{Vector{<:LegendElement},Vector{String}},Nothing}

    # Copy input data
    v_lines = copy(flags[1])
    line_labels = copy(flags[2])

    # Filter out values outside of the original plot's range
    rangeCut!(v_lines, line_labels, xlimits(figure); keep_edges = false)

    # Don't draw anything if all the flags are outside the original plot
    !isempty(v_lines) || return nothing

    # Compute as many distinguishable colors as vertical lines are left
    colors = distinguishable_colors(length(v_lines), [RGB(1, 1, 1), RGB(0, 0, 0)], dropseed = true)

    # Draw the vertical lines
    for (v_line, color, linestyle) in zip(v_lines, colors, LINE_STYLES)
        vlines!(figure.current_axis.x, v_line; color, linestyle, linewidth = 3)
    end

    # Construct markers for the legend
    line_elements = [
        LineElement(; color, linestyle, linewidth = 4) for
        (color, linestyle) in zip(colors, LINE_STYLES)
    ]

    return line_elements, line_labels

end

@doc raw"""
    ksLine!(
        figure::Makie.Figure; 
        <keyword arguments>
    )::Tuple{Vector{<:LegendElement},Vector{String}}

Draw two line plots and three lines of text. 

One line plot is the fit of the data in `figure`, while the other is the result of Kennicutt (1998). 
The lines of text will be drawn in the upper left corner of the plot, the first with the equation 

```math
\Sigma_\mathrm{SFR} = \mathrm{A}\,\Sigma_\mathrm{gas}^{\,\mathrm{N}} \, ,
```
and the second and third giving the values of ``\mathrm{N}`` and ``\mathrm{log}_{10}(\mathrm{A})``
respectively.

!!! note
    `x_unit`, `y_unit`, `x_log` and `y_log` have to be set with the same values as 
    they have been in [`kennicuttSchmidt`](@ref), otherwise the results will not be correct.

# Arguments
- `figure::Makie.Figure`: Figure to be modified.
- `x_unit::Unitful.Units = UnitfulAstro.Msun * UnitfulAstro.pc^-2`: Unit of the surface density of 
  gas used in `figure`.
- `y_unit::Unitful.Units = UnitfulAstro.Msun * UnitfulAstro.yr^-1 * UnitfulAstro.kpc^-2`: Unit of 
  the surface density of star formation rate used in `figure`. 
- `x_log::Bool = true`: If the x axis will be ``\log(\Sigma_\mathrm{gas})`` (default) or just ``\Sigma_\mathrm{gas}``.
- `y_log::Bool = true`: If the y axis will be ``\log(\Sigma_\mathrm{sfr})`` (default) or just ``\Sigma_\mathrm{sfr}``.
- `error_formating::String = "std_error"`: What to print as the error of the fit. The options are:
  - `"std_error"` ⟶ mean ± standard_error.
  - `"conf_interval"` ⟶ mean ± max(upper``_{95\%}`` - mean, mean - lower``_{95\%}``). 
- `colors::NTuple{2,ColorInput} = (:red, :blue)`: The color of the two line plots. The first is for 
  the result of Kennicutt (1998) and the second is for the linear fit.
- `linestyles::NTuple{2,LineStyleInput} = (nothing, :dash)`: The line styles of the two line plots. 
  The first is for the result of Kennicutt (1998) and the second is for the linear fit.

# Returns
- Tuple with two elements:
  - `LineElement` to be used as markers in the legend.
  - Labels.

"""
function ksLine!(
    figure::Makie.Figure;
    x_unit::Unitful.Units = UnitfulAstro.Msun * UnitfulAstro.pc^-2,
    y_unit::Unitful.Units = UnitfulAstro.Msun * UnitfulAstro.yr^-1 * UnitfulAstro.kpc^-2,
    x_log::Bool = true,
    y_log::Bool = true,
    error_formating::String = "std_error",
    colors::NTuple{2,ColorInput} = (:red, :blue),
    linestyles::NTuple{2,LineStyleInput} = (nothing, :dash),
)::Tuple{Vector{<:LegendElement},Vector{String}}

    # Get the axis to be modified
    axis = figure.current_axis.x
    # Get the data points in the plot
    data = pointData(figure)

    # Get the x coordinate of the points
    x_data = [Float64(point[1]) for point in data]
    # Get the y coordinate of the points
    y_data = [Float64(point[2]) for point in data]

    # Construct markers for the legend
    line_elements = [
        LineElement(; color, linestyle, linewidth = 4) for
        (color, linestyle) in zip(colors, linestyles)
    ]

    ################################################################################################
    # Line using the fit by Kennicutt (1998)
    ################################################################################################

    x_values = x_log ? 10 .^ x_data : x_data
    x_KS_98 = ustrip.(KS98_RHO_UNIT, x_values .* x_unit)

    y_values = ustrip.(y_unit, KS98_INTERCEPT .* x_KS_98 .^ KS98_SLOPE)
    y_KS_98 = y_log ? log10.(y_values) : y_values

    lines!(axis, x_data, y_KS_98, color = colors[1], linestyle = linestyles[1])

    ################################################################################################
    # Linear fit of the data in `figure`
    ################################################################################################

    x_fit_data = x_log ? x_data : log10.(x_data)
    y_fit_data = y_log ? y_data : log10.(y_data)

    # Compute linear fit
    X = [ones(length(x_fit_data)) x_fit_data]
    linear_model = lm(X, y_fit_data)

    # Mean values
    coeff = coef(linear_model)
    mean_intercept = coeff[1]
    mean_slope = coeff[2]

    # Get predicted values
    y_fit = y_log ? predict(linear_model, X) : exp10.(predict(linear_model, X))

    # Plot linear fit
    lines!(axis, x_data, y_fit, color = colors[2], linestyle = linestyles[2])

    ################################################################################################
    # Annotations
    ################################################################################################

    # Errors
    if error_formating == "conf_interval"
        conf_int = confint(linear_model)
        intercept_err = max(conf_int[1, 2] - mean_intercept, mean_intercept - conf_int[1, 1])
        slope_err = max(conf_int[2, 2] - mean_slope, mean_slope - conf_int[2, 1])
    elseif error_formating == "std_error"
        intercept_err, slope_err = stderror(linear_model)
    else
        throw(ArgumentError("`error_formating` has to be 'conf_interval' or 'std_error'"))
    end

    intercept, intercept_error = formatError(mean_intercept, intercept_err)
    slope, slope_error = formatError(mean_slope, slope_err)

    # Draw annotations with the fitted parameters and their errors
    annotations!(
        axis,
        [
            "Σsfr = A * Σρ^N",
            "N = $slope ± $slope_error",
            "log10(A) = $intercept ± $intercept_error",
        ],
        [
            Point(absCoor(figure, 0.05, 0.93)...),
            Point(absCoor(figure, 0.05, 0.88)...),
            Point(absCoor(figure, 0.05, 0.83)...),
        ],
    )

    return line_elements, ["Kennicutt et al. (1998)", "Fit"]

end

"""
    cross!(
        figure::Makie.Figure; 
        <keyword arguments>
    )::Nothing

Draw two lines, one horizontal and one vertical, crossing at `cross`.

# Arguments
- `figure::Makie.Figure`: Figure to be modified.
- `cross::NTuple{2,Float64} = (0.0, 0.0)`: Crossing point of the lines.
- `color::Symbol = :black`: Color of the lines.

"""
function cross!(
    figure::Makie.Figure;
    cross::NTuple{2,Float64} = (0.0, 0.0),
    color::Symbol = :black,
)::Nothing

    # Plot limits
    xlims = xlimits(figure)
    ylims = ylimits(figure)

    if xlims[1] < cross[1] < xlims[2]
        # Draw the vertical line
        vlines!(figure.current_axis.x, cross[1]; color, linewidth = 2)
    end
    if ylims[1] < cross[2] < ylims[2]
        # Draw the horizontal line
        hlines!(figure.current_axis.x, cross[2]; color, linewidth = 2)
    end

    return nothing

end

"""
    compareExperiment!(
        figure::Makie.Figure,
        quantity::String,
        x_unit::Unitful.Units,
        y_unit::Unitful.Units; 
        <keyword arguments>
    )::Tuple{Vector{<:LegendElement},Vector{String}}

Draw the experimental data compiled by Mollá et al. (2015).

# Arguments
- `figure::Makie.Figure`: Figure to be modified.
- `quantity:String`: The physical quantity for the y axis. The options are "ΣHI", "ΣHII", "logΣ*", 
  "logΣsfr", "C/H", "N/H", "O/H". See the documentation on the function `getMollá2015` in 
  `src/data_acquisition.jl` for details. The x axis will be the radial distance.
- `x_unit::Unitful.Units`: Unit of lenght for the x axis of `figure`.
- `y_unit::Unitful.Units`: Unit of the quantity in the y axis of `figure`. If the quantity is 
  inside a logarithm, you have to put the unit used to normalize it.
- `source_path::String = "./"`: Path to the directory containing the `data.csv` file.
- `color::ColorInput = :red`: The color of the markers and error whiskers.
- `markersize::Int64 = 6`: Size of the markers.
- `whiskerwidth::Int64 = 10`: Width of the error whiskers.   

# Returns
- Tuple with two elements:
    - `MarkerElement` to be used as markers in the legend.
    - Labels.

"""
function compareExperiment!(
    figure::Makie.Figure,
    quantity::String,
    x_unit::Unitful.Units,
    y_unit::Unitful.Units;
    source_path::String = "./",
    color::ColorInput = :red,
    markersize::Int64 = 6,
    whiskerwidth::Int64 = 10,
)::Tuple{Vector{<:LegendElement},Vector{String}}

    quantity ∈ ["ΣHI", "ΣHII", "logΣ*", "logΣsfr", "C/H", "N/H", "O/H"] || throw(
        ArgumentError("The quantity $quantity is not supported, see the documentation for a \
        list of supported quantities.")
    )

    # Get the axis to be modified
    axis = figure.current_axis.x

    # Load the experimental data
    data, units = getMollá2015(source_path)

    # Convert to target units
    x_data = ustrip.(x_unit, data[!, "R"] .* units["R"][1])
    if units[quantity][2]
        y_data = log10.(ustrip.(y_unit, (10 .^ data[!, quantity]) .* units[quantity][1]))
        error = log10.(ustrip.(y_unit, (10 .^ data[!, quantity*"_error"]) .* units[quantity][1]))
    else
        y_data = ustrip.(y_unit, data[!, quantity] .* units[quantity][1])
        error = ustrip.(y_unit, data[!, quantity*"_error"] .* units[quantity][1])
    end

    # Plot the data
    errorbars!(axis, x_data, y_data, error; color, whiskerwidth)
    scatter!(axis, x_data, y_data; markersize, color)

    return [MarkerElement(; color, marker = :circle, markersize)], ["Mollá et al. (2015)"]

end
