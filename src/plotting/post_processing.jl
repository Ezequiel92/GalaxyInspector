####################################################################################################
# Post-processing functions
####################################################################################################
#
# A post-processing function must take a Makie figure, add something to it, and return how to label
# the additions (or `nothing` when no new labels should be drawn).
#
# Expected signature:
#
#   post_processing(figure, args...; kwargs...) -> ([marker, ...], [label, ...]) or `nothing`
#
# where:
#
#   - figure::Makie.Figure
#   - marker::LegendElement
#   - label::AbstractString
#
####################################################################################################

"""
    ppVerticalFlags!(
        figure::Makie.Figure,
        positions::Vector{<:Real};
        <keyword arguments>
    )::Union{Tuple{Vector{<:LegendElement},Vector{AbstractString}},Nothing}

Draw vertical lines.

# Arguments

  - `figure::Makie.Figure`: Makie figure.
  - `positions::Vector{<:Real}`: The x coordinates of the lines.
  - `colors::Vector{<:ColorType}=[WONG_RED]`: Colors of the lines.
  - `line_styles::Vector{<:LineStyleType}=[:solid]`: Styles of the lines.
  - `labels::Union{Vector{AbstractString},Nothing}=nothing`: Labels for the lines. If set to `nothing` no label is printed.

# Returns

  - A tuple with the elements for the legend:

      + `LineElement`s to be used as the marker.
      + The labels.
"""
function ppVerticalFlags!(
    figure::Makie.Figure,
    positions::Vector{<:Real};
    colors::Vector{<:ColorType}=[WONG_RED],
    line_styles::Vector{<:LineStyleType}=[:solid],
    labels::Union{Vector{AbstractString},Nothing}=nothing,
)::Union{Tuple{Vector{<:LegendElement},Vector{AbstractString}},Nothing}

    # Filter out values outside the range of the original plot
    rangeCut!(positions, xlimits!(figure); keep_edges=false)

    if isempty(positions)

        # Don't draw anything if all the flags are outside the original plot window
        logging[] && @warn("ppVerticalFlags!: All vertical lines lie outside the plot range")

        return nothing

    else

        # Draw the vertical lines
        for (i, position) in pairs(positions)

            color     = ring(colors, i)
            linestyle = ring(line_styles, i)

            vl = vlines!(figure.current_axis.x, position; color, linestyle)

            translate!(Accum, vl, 0, 0, -10)

        end

    end

    isnothing(label) && return nothing

    return (
        [
            LineElement(; color=ring(colors, i), linestyle=ring(line_styles, i)) for
            i in eachindex(labels)
        ],
        labels,
    )

end

"""
    ppHorizontalFlags!(
        figure::Makie.Figure,
        positions::Vector{<:Real};
        <keyword arguments>
    )::Union{Tuple{Vector{<:LegendElement},Vector{AbstractString}},Nothing}

Draw horizontal lines.

# Arguments

  - `figure::Makie.Figure`: Makie figure.
  - `positions::Vector{<:Real}`: The y coordinates of the lines.
  - `colors::Vector{<:ColorType}=[WONG_RED]`: Colors of the lines.
  - `line_styles::Vector{<:LineStyleType}=[:solid]`: Styles of the lines.
  - `labels::Union{Vector{AbstractString},Nothing}=nothing`: Labels for the lines. If set to `nothing` no label is printed.

# Returns

  - A tuple with the elements for the legend:

      + `LineElement`s to be used as the marker.
      + The labels.
"""
function ppHorizontalFlags!(
    figure::Makie.Figure,
    positions::Vector{<:Real};
    colors::Vector{<:ColorType}=[WONG_RED],
    line_styles::Vector{<:LineStyleType}=[:solid],
    labels::Union{Vector{AbstractString},Nothing}=nothing,
)::Union{Tuple{Vector{<:LegendElement},Vector{AbstractString}},Nothing}

    # Filter out values outside the range of the original plot
    rangeCut!(positions, ylimits!(figure); keep_edges=false)

    if isempty(positions)

        # Don't draw anything if all the flags are outside the original plot window
        logging[] && @warn("ppHorizontalFlags!: All horizontal lines lie outside the plot range")

        return nothing

    else

        # Draw the horizontal lines
        for (i, position) in pairs(positions)

            color     = ring(colors, i)
            linestyle = ring(line_styles, i)

            hl = hlines!(figure.current_axis.x, position; color, linestyle)

            translate!(Accum, hl, 0, 0, -10)

        end

    end

    isnothing(labels) && return nothing

    return (
        [
            LineElement(; color=ring(colors, i), linestyle=ring(line_styles, i)) for
            i in eachindex(labels)
        ],
        labels,
    )

end

"""
    ppCross!(
        figure::Makie.Figure,
        cross_point::Tuple{<:Real,<:Real};
        <keyword arguments>
    )::Union{Tuple{Vector{<:LegendElement},Vector{AbstractString}},Nothing}

Draw two lines, one horizontal and one vertical.

# Arguments

  - `figure::Makie.Figure`: Makie figure.
  - `cross_point::Tuple{<:Real,<:Real}`: Crossing point of the lines.
  - `color::ColorType=WONG_RED`: Color of the lines.
  - `linestyle::LineStyleType=:solid`: Style of the lines.
  - `label::Union{AbstractString,Nothing}=nothing`: Label for the lines. If set to `nothing` no label is printed.

# Returns

  - A tuple with the elements for the legend:

      + A `LineElement` to be used as the marker.
      + The label.
"""
function ppCross!(
    figure::Makie.Figure,
    cross_point::Tuple{<:Real,<:Real};
    color::ColorType=WONG_RED,
    linestyle::LineStyleType=:solid,
    label::Union{AbstractString,Nothing}=nothing,
)::Union{Tuple{Vector{<:LegendElement},Vector{AbstractString}},Nothing}

    x_limits = xlimits!(figure)
    y_limits = ylimits!(figure)

    # Draw the vertical line
    if x_limits[1] < cross_point[1] < x_limits[2]

        vl = vlines!(figure.current_axis.x, cross_point[1]; color, linestyle)

        translate!(Accum, vl, 0, 0, -10)

    else

        logging[] && @warn("ppCross!: The vertical line lies outside the plot range")

    end

    # Draw the horizontal line
    if y_limits[1] < cross_point[2] < y_limits[2]

        hl = hlines!(figure.current_axis.x, cross_point[2]; color, linestyle)

        translate!(Accum, hl, 0, 0, -10)

    else

        logging[] && @warn("ppCross!: The horizontal line lies outside the plot range")

    end

    isnothing(label) && return nothing

    return [LineElement(; color, linestyle)], [label]

end

"""
    ppABFlags!(
        figure::Makie.Figure,
        slopes::Vector{<:Real}
        intercepts::Vector{<:Real};
        <keyword arguments>
    )::Nothing

Draw lines defined by f(x) = slope * x + intercept.

# Arguments

  - `figure::Makie.Figure`: Makie figure to be drawn over.
  - `slopes::Vector{<:Real}`: List of slopes.
  - `intercepts::Vector{<:Real}`: List of intercepts.
  - `colors::Vector{<:ColorType}=[:red]`: Colors of the lines.
  - `line_styles::Vector{<:LineStyleType}=[:solid]`: Styles of the lines.
"""
function ppABFlags!(
    figure::Makie.Figure,
    slopes::Vector{<:Real},
    intercepts::Vector{<:Real};
    colors::Vector{<:ColorType}=[:red],
    line_styles::Vector{<:LineStyleType}=[:solid],
)::Nothing

    x_limits = xlimits!(figure)
    y_limits = ylimits!(figure)

    for (i, (slope, intercept)) in enumerate(zip(slopes, intercepts))

        y_low = slope * x_limits[1] + intercept
        y_high = slope * x_limits[2] + intercept

        if min(y_low, y_high) > y_limits[2] || max(y_low, y_high) < y_limits[1]
            (
                logging[] &&
                @warn("ppABFlags!: The line y = $(slope) * x + $(intercept) lies outside the plot \
                range")
            )
            continue
        end

        color     = ring(colors, i)
        linestyle = ring(line_styles, i)

        abl = ablines!(figure.current_axis.x, intercept, slope; color, linestyle)

        translate!(Accum, abl, 0, 0, -10)

    end

    return nothing

end

"""
    ppFillBelowLine!(
        figure::Makie.Figure;
        <keyword arguments>
    )::Union{Tuple{Vector{<:LegendElement},Vector{AbstractString}},Nothing}

Fill the space below a line plot with a solid `color` down to a `lower_limit`.

# Arguments

  - `figure::Makie.Figure`: Makie figure.
  - `lower_limit::Number=0.0`: Lower bound.
  - `color::ColorType=WONG_RED`: Color.
  - `alpha::Float64=0.3`: Level of transparency. 1.0 is completely opaque and 0.0 is completely transparent.
  - `label::Union{AbstractString,Nothing}=nothing`: Label for the shaded region. If set to `nothing` no label is printed.

# Returns

  - A tuple with the elements for the legend:

      + A `MarkerElement` to be used as the marker.
      + The label.
"""
function ppFillBelowLine!(
    figure::Makie.Figure;
    lower_limit::Number=0.0,
    color::ColorType=WONG_RED,
    alpha::Float64=0.3,
    label::Union{AbstractString,Nothing}=nothing,
)::Union{Tuple{Vector{<:LegendElement},Vector{AbstractString}},Nothing}

    # Read the data points in the plot
    points = pointData(figure)

    if isempty(points)
        logging[] && @warn("ppFillBelowLine!: There are no points in the figure")
        return nothing
    end

    xs = Vector{Float64}(undef, length(points))
    ys = Vector{Float64}(undef, length(points))

    for (i, point) in pairs(points)
        xs[i] = point[1]
        ys[i] = point[2]
    end

    ys_low = fill(lower_limit, length(xs))

    bp = band!(figure.current_axis.x, xs, ys_low, ys; color, alpha)

    translate!(Accum, bp, 0, 0, -10)

    isnothing(label) && return nothing

    return [MarkerElement(; color, alpha, marker=:rect)], [label]

end

"""
    ppArrows!(
        figure::Makie.Figure,
        positions::Vector{NTuple{4,Float64}};
        <keyword arguments>
    )::Nothing

Draw arrows.

# Arguments

  - `figure::Makie.Figure`: Makie figure.
  - `positions::Vector{NTuple{4,Float64}}`: The x, y, u, and v coordinates of the arrows. x and y indicate the position of the tails, and u and v the position of the point of the arrows.
  - `colors::Vector{<:ColorType}=[WONG_RED]`: Colors of the arrows.
"""
function ppArrows!(
    figure::Makie.Figure,
    positions::Vector{NTuple{4,Float64}};
    colors::Vector{<:ColorType}=[WONG_RED],
)::Nothing

    for (i, position) in pairs(positions)

        x, y, u, v = position
        color = ring(colors, i)

        # Draw the arrows
        arrows2d!(figure.current_axis.x, [x], [y], [u], [v]; color)

    end

    return nothing

end

"""
    ppAnnotation!(figure::Makie.Figure, text::AbstractString; <keyword arguments>)::Nothing

Add an annotation to the plot.

# Arguments

  - `figure::Makie.Figure`: Makie figure.
  - `text::AbstractString`: Text to be written.
  - `position::Tuple{<:Real,<:Real}=(0.04, 0.98)`: Relative position of the top left corner of the text box.
  - `color=:black`: Text color.
  - `fontsize::Int=35`: Font size.
"""
function ppAnnotation!(
    figure::Makie.Figure,
    text::AbstractString;
    position::Tuple{<:Real,<:Real}=(0.04, 0.98),
    color=:black,
    fontsize::Int=35,
)::Nothing

    text!(
        figure.current_axis.x,
        position[1],
        position[2];
        text,
        align=(:left, :top),
        color,
        space=:relative,
        fontsize,
    )

    return nothing

end

"""
    ppBarPlotLabels(
        ::Makie.Figure,
        components::Vector{Symbol};
        <keyword arguments>
    )::Tuple{Vector{<:LegendElement},Vector{AbstractString}}

Return the legend elements for the plot made by [`gasBarPlot`](@ref).

# Arguments

  - `::Makie.Figure`: Makie figure.
  - `components::Vector{Symbol}`: List of gas components to be considered. See [`COMPONENTS`](@ref) for options.
  - `colors=Makie.wong_colors()`: Colors of the bars.

# Returns

  - A tuple with the elements for the legend:

      + `PolyElement`s to be used in the legend.
      + The label strings.
"""
function ppBarPlotLabels(
    ::Makie.Figure,
    components::Vector{Symbol};
    colors=Makie.wong_colors(),
)::Tuple{Vector{<:LegendElement},Vector{AbstractString}}

    labels = [plotParams(Symbol(component, :_fraction)).var_name for component in components]

    return [PolyElement(polycolor=colors[i]) for i in 1:length(labels)], labels

end

@doc raw"""
    ppFitLine!(
        figure::Makie.Figure;
        <keyword arguments>
    )::Union{Tuple{Vector{<:LegendElement},Vector{AbstractString}},Nothing}

Draw a linear fit for the data in `figure`.

An annotation with the equation $y = a \, x + b$, and the fitted values for $a$ and $b$, will be positioned in the upper right corner of the plot.

# Arguments

  - `figure::Makie.Figure`: Makie figure.
  - `text_position::Tuple{<:Real,<:Real}=(0.04, 0.98)`: Relative position of the top left corner of the text box.
  - `text_align::NTuple{2,Symbol}=(:left, :top)`: Corner of the legend box whose positions is set with `text_position`.
  - `wts::Union{Vector{Float64},Nothing}=nothing`: Weights for the fits. Set to `nothing` for a non-weighted fit.
  - `error_formating::Symbol=:std_error`: Error format for the annotation. The options are:

      + `:std_error`     -> mean ± standard_error.
      + `:conf_interval` -> mean ± max(upper$_{95\%}$ - mean, mean - lower$_{95\%}$).

  - `color::ColorType=WONG_RED`: Color of the line.
  - `linestyle::LineStyleType=:solid`: Style of the line.
  - `linewidth::Int=3`: Line width.

# Returns

  - A tuple with the elements for the legend:

      + A `LineElement` to be used as the marker.
      + The label.
"""
function ppFitLine!(
    figure::Makie.Figure;
    text_position::Tuple{<:Real,<:Real}=(0.04, 0.98),
    text_align::NTuple{2,Symbol}=(:left, :top),
    wts::Union{Vector{Float64},Nothing}=nothing,
    error_formating::Symbol=:std_error,
    color::ColorType=WONG_RED,
    linestyle::LineStyleType=:solid,
    linewidth::Int=3,
)::Union{Tuple{Vector{<:LegendElement},Vector{AbstractString}},Nothing}

    # Read the data points in the plot
    points = pointData(figure)

    if isempty(points)

        logging[] && @warn("ppFitLine!: There are no points in the figure")

        return nothing

    end

    sort!(points; by=x -> x[1])

    # Get the scaling of each axis
    x_scaling = xscale(figure)
    y_scaling = yscale(figure)

    ################################################################################################
    # Linear fit
    ################################################################################################

    # Get the x coordinates of the points
    x_points = x_scaling.(Float64[point[1] for point in points])

    # Get the y coordinates of the points
    y_points = y_scaling.(Float64[point[2] for point in points])

    # Find NaN points
    xnan_idxs = map(isnan, x_points)
    ynan_idxs = map(isnan, y_points)

    nan_idxs = xnan_idxs ∪ ynan_idxs

    # Delete NaN points
    deleteat!(x_points, nan_idxs)
    deleteat!(y_points, nan_idxs)

    # Compute the linear fit
    X = [ones(length(x_points)) x_points]

    if isnothing(wts)
        linear_model = lm(X, y_points)
    else
        linear_model = lm(X, y_points; wts)
    end

    # Read the fitted coeficients
    coeff = coef(linear_model)
    intercept_mean = coeff[1]
    slope_mean = coeff[2]

    y_limits = [extrema(y_points)...]
    x_limits = @. (y_limits - intercept_mean) / slope_mean

    # Plot the linear fit
    lp = lines!(
        figure.current_axis.x,
        Makie.inverse_transform(x_scaling).(x_limits),
        Makie.inverse_transform(y_scaling).(y_limits);
        color,
        linestyle,
        linewidth,
    )

    ################################################################################################
    # Annotation
    ################################################################################################

    # Compute the fitting errors
    if error_formating == :conf_interval

        conf_int = confint(linear_model)
        intercept_error = max(conf_int[1, 2] - intercept_mean, intercept_mean - conf_int[1, 1])
        slope_error = max(conf_int[2, 2] - slope_mean, slope_mean - conf_int[2, 1])

    elseif error_formating == :std_error

        intercept_error, slope_error = stderror(linear_model)

    else

        throw(ArgumentError("ppFitLine!: `error_formating` can only be :conf_interval or \
        :std_error, but I got :$(error_formating)"))

    end

    # Format the mean and error values
    intercept, δintercept = formatError(intercept_mean, intercept_error)
    slope, δslope = formatError(slope_mean, slope_error)

    # Set the position of the selected top corner of the legend box
    x_tp = text_position[1]
    y_tp = text_position[2]

    # Draw the annotation
    text!(
        figure.current_axis.x,
        [(x_tp, y_tp), (x_tp, y_tp - 0.05), (x_tp, y_tp - 0.1)];
        text=[L"y = a \, x + b", L"a = %$slope \pm %$δslope", L"b = %$intercept \pm %$δintercept"],
        align=text_align,
        color,
        space=:relative,
    )

    ################################################################################################
    # Band
    ################################################################################################

    # Compute the intercept and slope as numbers with uncertainties
    band_intercept = intercept_mean ± intercept_error
    band_slope     = slope_mean ± slope_error

    # Compute the values of the y axis as numbers with uncertainties
    x_band = collect(range(x_limits[1], x_limits[2], 100))
    y_band = @. Makie.inverse_transform(y_scaling)(x_band * band_slope + band_intercept)

    values        = Measurements.value.(y_band)
    uncertainties = Measurements.uncertainty.(y_band)

    # Draw the uncertainty band
    bp = band!(
        figure.current_axis.x,
        x_band,
        values .- uncertainties,
        values .+ uncertainties;
        color,
    )

    translate!(Accum, bp, 0, 0, -10)
    translate!(Accum, lp, 0, 0, -10)

    return nothing

end

###############
# Observations
###############

"""
    ppKennicutt1998!(
        figure::Makie.Figure;
        <keyword arguments>
    )::Union{Tuple{Vector{<:LegendElement},Vector{AbstractString}},Nothing}

Draw a line plot with the fit for the Kennicutt-Schmidt relation from Kennicutt (1998).

# Arguments

  - `figure::Makie.Figure`: Makie figure.
  - `x_unit::Unitful.Units=u"Msun * pc^-2"`: Unit for the area density of gas used in `figure`.
  - `y_unit::Unitful.Units=u"Msun * yr^-1 * kpc^-2"`: Unit for the area density of star formation rate used in `figure`.
  - `x_log::Bool=true`: If the x axis is ``\\log_{10}(\\Sigma_\\mathrm{HI + H_2})`` (`x_log` = true) or just ``\\Sigma_\\mathrm{HI + H_2}`` (`x_log` = false).
  - `y_log::Bool=true`: If the y axis is ``\\log_{10}(\\Sigma_\\mathrm{SFR})`` (`y_log` = true) or just ``\\Sigma_\\mathrm{SFR}`` (`y_log` = false).
  - `extend::Float64=0.0`: By default the y axis limits of the line will be the vertical range of points in the plot. This can be extended by the fraction `extend` of the vertical range.
  - `color::ColorType=WONG_RED`: Color for the line.
  - `linestyles::Vector{<:LineStyleType}=[:solid, :dash],`: Styles for the lines. The first style will indicate the range for which there are experimental data, and the second one will be for the extrapolation.
  - `linewidth::Int=3`: Line width.

# Returns

  - A tuple with the elements for the legend:

      + A `LineElement` to be used as the marker.
      + The label.

# References

R. C. Kennicutt (1998). *The Global Schmidt Law in Star-forming Galaxies*. The Astrophysical Journal, **498(2)**, 541-552. [doi:10.1086/305588](https://doi.org/10.1086/305588)
"""
function ppKennicutt1998!(
    figure::Makie.Figure;
    x_unit::Unitful.Units=u"Msun * pc^-2",
    y_unit::Unitful.Units=u"Msun * yr^-1 * kpc^-2",
    x_log::Bool=true,
    y_log::Bool=true,
    extend::Float64=0.0,
    color::ColorType=WONG_RED,
    linestyles::Vector{<:LineStyleType}=[:solid, :dash],
    linewidth::Int=3,
)::Union{Tuple{Vector{<:LegendElement},Vector{AbstractString}},Nothing}

    # Read the data points in the plot
    points = pointData(figure)

    if isempty(points)
        logging[] && @warn("ppKennicutt1998!: There are no points in the figure")
        return nothing
    end

    ################################################################################################
    # Range of Σsfr in the plot
    ################################################################################################

    # Get the extrema of the y coordinates
    y_limits = [extrema(Float64[point[2] for point in points])...]

    # Compute the vertical range, and if it is too small extend it
    y_range = y_limits[2] - y_limits[1]

    if iszero(extend) && (abs(y_range / y_limits[1]) < 0.1)
        extend = 0.4 / abs(y_range / y_limits[1])
    end

    # Extend the y limits the required amount
    if extend > 0.0

        extension = y_range * extend

        y_limits[1] -= extension
        y_limits[2] += extension

    end

    ################################################################################################
    # Range of Σsfr in Kennicutt (1998)
    ################################################################################################

    # Set the correct unit and scale for the range of experimental values
    if y_log
        kennicutt_range = log10.(ustrip.(y_unit, KS98_SFR_RANGE))
    else
        kennicutt_range = ustrip.(y_unit, KS98_SFR_RANGE)
    end

    # Set the y ranges to be plotted
    y_ranges = [kennicutt_range]

    # If there is extrapolation add new ranges
    if y_limits[2] > kennicutt_range[2]
        push!(y_ranges, [kennicutt_range[2], y_limits[2]])
    end
    if y_limits[1] < kennicutt_range[1]
        push!(y_ranges, [y_limits[1], kennicutt_range[1]])
    end

    ################################################################################################
    # Plot the fit from Kennicutt (1998)
    ################################################################################################

    for (y_zone, linestyle) in zip(y_ranges, [linestyles..., linestyles[2]])

        # Compute the extrema of the star formation area density
        if y_log
            Σsfr = @. exp10(y_zone) * y_unit
        else
            Σsfr = @. y_zone * y_unit
        end

        # Compute the extrema of the gas mass area density
        Σgas     = invKSLawKennicutt1998(Σsfr; log_output=false)
        x_limits = Measurements.value.(Σgas)

        # Compute the values for the x axis
        x_points = collect(range(x_limits[1], x_limits[2], 100))

        if x_log
            x_axis = @. log10(ustrip(x_unit, x_points))
        else
            x_axis = @. ustrip(x_unit, x_points)
        end

        # Compute the values for the y axis
        y_points = KSLawKennicutt1998(x_points; log_output=false)

        if y_log
            y_axis = @. log10(ustrip(y_unit, y_points))
        else
            y_axis = @. ustrip(y_unit, y_points)
        end

        values        = Measurements.value.(y_axis)
        uncertainties = Measurements.uncertainty.(y_axis)

        bp = band!(
            figure.current_axis.x,
            x_axis,
            values .- uncertainties,
            values .+ uncertainties;
            color,
        )

        lp = lines!(figure.current_axis.x, x_axis, values; color, linestyle, linewidth)

        # Put the post processing elements at the back of the plot
        translate!(Accum, bp, 0, 0, -10)
        translate!(Accum, lp, 0, 0, -10)

    end

    return [LineElement(; color, linestyle=linestyles[1], linewidth)], ["Kennicutt (1998)"]

end

"""
    ppBigiel2008!(
        figure::Makie.Figure,
        molecular::Bool;
        <keyword arguments>
    )::Union{Tuple{Vector{<:LegendElement},Vector{AbstractString}},Nothing}

Draw a line plot with the fit for the Kennicutt-Schmidt relation from Bigiel et al. (2008).

!!! note

    The resolution used in Bigiel et al. (2008) is 750 pc (see Section 1).

# Arguments

  - `figure::Makie.Figure`: Makie figure.
  - `molecular::Bool`: If the x axis will be the area mass density of molecular hydrogen, or, if set to `false`, the area mass density of neutral hydrogen.
  - `x_unit::Unitful.Units=u"Msun * pc^-2"`: Unit for the area density of gas used in `figure`.
  - `y_unit::Unitful.Units=u"Msun * yr^-1 * kpc^-2"`: Unit for the area density of star formation rate used in `figure`.
  - `x_log::Bool=true`: If the x axis is ``\\log_{10}(\\Sigma_\\mathrm{H})`` (`x_log` = true) or just ``\\Sigma_\\mathrm{H}`` (`x_log` = false).
  - `y_log::Bool=true`: If the y axis is ``\\log_{10}(\\Sigma_\\mathrm{SFR})`` (`y_log` = true) or just ``\\Sigma_\\mathrm{SFR}`` (`y_log` = false).
  - `extend::Float64=0.0`: By default the y axis limits will be the vertical range of the points in the plot. This can be extended by the multiplicative factor `extend` of the vertical range.
  - `color::ColorType=WONG_RED`: Color for the line.
  - `linestyles::Vector{<:LineStyleType}=[:solid, :dash],`: Styles for the lines. The first style will indicate the range for which there are experimental data, and the second one will be for the extrapolation.
  - `linewidth::Int=3`: Line width.

# Returns

  - A tuple with the elements for the legend:

      + A `LineElement` to be used as the marker.
      + The label.

# References

F. Bigiel et al. (2008). *THE STAR FORMATION LAW IN NEARBY GALAXIES ON SUB-KPC SCALES*. The Astrophysical Journal, **136(6)**, 2846. [doi:10.1088/0004-6256/136/6/2846](https://doi.org/10.1088/0004-6256/136/6/2846)
"""
function ppBigiel2008!(
    figure::Makie.Figure,
    molecular::Bool;
    x_unit::Unitful.Units=u"Msun * pc^-2",
    y_unit::Unitful.Units=u"Msun * yr^-1 * kpc^-2",
    x_log::Bool=true,
    y_log::Bool=true,
    extend::Float64=0.0,
    color::ColorType=WONG_RED,
    linestyles::Vector{<:LineStyleType}=[:solid, :dash],
    linewidth::Int=3,
)::Union{Tuple{Vector{<:LegendElement},Vector{AbstractString}},Nothing}

    # Read the data points in the plot
    points = pointData(figure)

    if isempty(points)
        logging[] && @warn("ppBigiel2008!: There are no points in the figure")
        return nothing
    end

    ################################################################################################
    # Range of Σsfr in the plot
    ################################################################################################

    # Get the extrema of the y coordinates
    y_limits = [extrema(Float64[point[2] for point in points])...]

    # Compute the vertical range, and if it is too small extend it
    y_range = y_limits[2] - y_limits[1]

    if iszero(extend) && (abs(y_range / y_limits[1]) < 0.1)
        extend = 0.4 / abs(y_range / y_limits[1])
    end

    # Extend the y limits the required amount
    if extend > 0.0

        extension = y_range * extend

        y_limits[1] -= extension
        y_limits[2] += extension

    end

    ################################################################################################
    # Range of Σsfr in Bigiel et al. (2008)
    ################################################################################################

    # Set the correct unit and scale for the range of experimental values
    if y_log
        bigiel_range = log10.(ustrip.(y_unit, BIGIEL2008_SFR_RANGE))
    else
        bigiel_range = ustrip.(y_unit, BIGIEL2008_SFR_RANGE)
    end

    # Set the y ranges to be plotted
    y_ranges = [bigiel_range]

    # If there is extrapolation add new ranges
    if y_limits[2] > bigiel_range[2]
        push!(y_ranges, [bigiel_range[2], y_limits[2]])
    end
    if y_limits[1] < bigiel_range[1]
        push!(y_ranges, [y_limits[1], bigiel_range[1]])
    end

    ################################################################################################
    # Plot the fit from Bigiel et al. (2008)
    ################################################################################################

    for (y_zone, linestyle) in zip(y_ranges, [linestyles..., linestyles[2]])

        # Compute the extrema of the star formation area density
        if y_log
            Σsfr = @. exp10(y_zone) * y_unit
        else
            Σsfr = @. y_zone * y_unit
        end

        # Compute the extrema of the gas mass area density
        ΣH       = invKSLawBigiel2008(Σsfr; molecular, log_output=false)
        x_limits = Measurements.value.(ΣH)

        # Compute the values for the x axis
        x_points = collect(range(x_limits[1], x_limits[2], 100))

        if x_log
            x_axis = @. log10(ustrip(x_unit, x_points))
        else
            x_axis = @. ustrip(x_unit, x_points)
        end

        # Compute the values for the y axis
        y_points = KSLawBigiel2008(x_points; molecular, log_output=false)

        if y_log
            y_axis = @. log10(ustrip(y_unit, y_points))
        else
            y_axis = @. ustrip(y_unit, y_points)
        end

        values        = Measurements.value.(y_axis)
        uncertainties = Measurements.uncertainty.(y_axis)

        bp = band!(
            figure.current_axis.x,
            x_axis,
            values .- uncertainties,
            values .+ uncertainties;
            color,
        )

        lp = lines!(figure.current_axis.x, x_axis, values; color, linestyle, linewidth)

        # Put the post processing elements at the back of the plot
        translate!(Accum, bp, 0, 0, -10)
        translate!(Accum, lp, 0, 0, -10)

    end

    return [LineElement(; color, linestyle=linestyles[1], linewidth)], ["Bigiel et al. (2008)"]

end

"""
    ppBigiel2010!(
        figure::Makie.Figure;
        <keyword arguments>
    )::Tuple{Vector{<:LegendElement},Vector{AbstractString}}

Draw a scatter plot of the SFR surface density vs gas surface density (Kennicutt-Schmidt law) for a given galaxy, using the data from Bigiel et al. (2010).

!!! note

    The resolution used in Bigiel et al. (2010) is 600 pc (see Section 2).

# Arguments

  - `figure::Makie.Figure`: Makie figure.
  - `galaxy::Union{String,Symbol}="NGC 628"`: Target galaxy. The options are:

      + With molecular and and atomic data (Table 2): "NGC 628", "NGC 3184", "NGC 3521", "NGC 4736", "NGC 5055", "NGC 5194", "NGC 6946".
      + With only atomic data (Table 3): "NGC 925", "NGC 2403", "NGC 2841", "NGC 2903", "NGC 3198", "NGC 3351", "NGC 3621", "NGC 3627", "NGC 5236", "NGC 5457", "NGC 7331", "NGC 7793".
      + :all: Every galaxy that is available for the given `quantity`. For more information on each galaxy see Bigiel et al. (2010).
  - `quantity::Symbol=:molecular`: Gas quantity for the x axis. The options are:

      + `:molecular` -> Surface density of molecular gas.
      + `:neutral`   -> Surface density of neutral gas.
      + `:atomic`    -> Surface density of atomic gas.
  - `x_unit::Unitful.Units=u"Msun * pc^-2"`: Unit for the area density of gas used in `figure`.
  - `y_unit::Unitful.Units=u"Msun * yr^-1 * kpc^-2"`: Unit for the area density of star formation rate used in `figure`.
  - `x_log::Bool=true`: If the x axis will be plotted as the ``\\log_{10}`` of the gas surface density.
  - `y_log::Bool=true`: If the y axis will be plotted as the ``\\log_{10}`` of the SFR surface density.
  - `color::ColorType=WONG_RED`: Color of the markers.

# Returns

  - A tuple with the elements for the legend:

      + A `MarkerElement` to be used as the marker.
      + The label.

# References

F. Bigiel et al. (2010). *EXTREMELY INEFFICIENT STAR FORMATION IN THE OUTER DISKS OF NEARBY GALAXIES*. The Astrophysical Journal, **140(5)**, 1194. [doi:10.1088/0004-6256/140/5/1194](https://doi.org/10.1088/0004-6256/140/5/1194)
"""
function ppBigiel2010!(
    figure::Makie.Figure;
    galaxy::Union{String,Symbol}="NGC 628",
    quantity::Symbol=:molecular,
    x_unit::Unitful.Units=u"Msun * kpc^-2",
    y_unit::Unitful.Units=u"Msun * yr^-1 *  kpc^-2",
    x_log::Bool=true,
    y_log::Bool=true,
    color::ColorType=WONG_RED,
)::Tuple{Vector{<:LegendElement},Vector{AbstractString}}

    ################################################################################################
    # Load table 2 from Bigiel et al. 2010
    ################################################################################################

    raw_data_2 = CSV.read(
        BIGIEL2010_TABLE_2,
        DataFrame;
        delim='\t',
        skipto=50,
        header=false,
        ignorerepeated=true,
        silencewarnings=true,
    )

    table_2 = DataFrame(
        gtype=String[],
        name=String[],
        HI=Union{Quantity,Missing}[],
        H2=Union{Quantity,Missing}[],
        SFR=Union{Quantity,Missing}[],
    )

    for row in eachrow(raw_data_2)

        data = row[1]

        gtype = strip(data[1:13])
        name  = strip(data[14:24])
        HI    = exp10(parserWS(data[26:29])) * u"Msun * pc^-2"
        H2    = exp10(parserWS(data[36:39])) * u"Msun * pc^-2"
        SFR   = exp10(parserWS(data[46:50])) * u"Msun * yr^-1 * kpc^-2"

        push!(table_2, [gtype name HI H2 SFR])

    end

    spirals_2 = filter(:gtype => isequal("Spirals"), table_2)

    ################################################################################################
    # Load table 3 from Bigiel et al. 2010
    ################################################################################################

    raw_data_3 = CSV.read(
        BIGIEL2010_TABLE_3,
        DataFrame;
        delim='\t',
        skipto=48,
        header=false,
        ignorerepeated=true,
        silencewarnings=true,
    )

    table_3 = DataFrame(
        gtype=String[],
        name=String[],
        HI=Union{Quantity,Missing}[],
        SFR=Union{Quantity,Missing}[],
    )

    for row in eachrow(raw_data_3)

        data = row[1]

        gtype = strip(data[1:8])
        name  = strip(data[9:19])
        HI    = exp10(parserWS(data[21:25])) * u"Msun * pc^-2"
        SFR   = parserWS(data[32:37]) * exp10(-5.0) * u"Msun * yr^-1 * kpc^-2"

        push!(table_3, [gtype name HI SFR])

    end

    spirals_3 = filter(:gtype => isequal("Spirals"), table_3)

    ################################################################################################
    # Find the target galaxy and read its values
    ################################################################################################

    if quantity == :molecular

        if galaxy == :all
            target_galaxy = spirals_2
        else
            target_galaxy = filter(:name => isequal(galaxy), spirals_2)
        end

        (
            isempty(target_galaxy) &&
            throw(ArgumentError("ppBigiel2010!: `galaxy` = $(galaxy) is not a spiral galaxy in \
            Table 2 of Bigiel et al. 2010 (notice that you have `quantity` = $(quantity))."))
        )

        Σg   = target_galaxy[!, :H2]
        Σsfr = target_galaxy[!, :SFR]

    elseif quantity == :neutral

        if galaxy == :all
            target_galaxy = spirals_2
        else
            target_galaxy = filter(:name => isequal(galaxy), spirals_2)
        end

        (
            isempty(target_galaxy) &&
            throw(ArgumentError("ppBigiel2010!: `galaxy` = $(galaxy) is not a spiral galaxy in \
            Table 2 of Bigiel et al. 2010 (notice that you have `quantity` = $(quantity))."))
        )

        ΣH2 = target_galaxy[!, :H2]
        ΣHI = target_galaxy[!, :HI]

        Σg   = ΣH2 .+ ΣHI
        Σsfr = target_galaxy[!, :SFR]

    elseif quantity == :atomic

        if galaxy == :all

            # Find the galaxies that are exclusive to table 3
            galaxies_t3 = setdiff(unique(spirals_3[!, :name]), unique(spirals_2[!, :name]))

            # Join the data in table 2 and table 3, ignoring repeated galaxies
            target_galaxy = vcat(
                spirals_2[!, [:HI, :SFR]],
                filter(:name => in(galaxies_t3), spirals_3)[!, [:HI, :SFR]],
            )

        else

            target_galaxy = filter(:name => isequal(galaxy), spirals_2)

        end

        if isempty(target_galaxy)

            target_galaxy = filter(:name => isequal(galaxy), spirals_3)

            (
                isempty(target_galaxy) &&
                throw(ArgumentError("ppBigiel2010!: `galaxy` $(galaxy) is not a spiral galaxy in \
                Table 2 or 3 of Bigiel et al. 2010."))
            )

        end

        Σg   = target_galaxy[!, :HI]
        Σsfr = target_galaxy[!, :SFR]

    else

        throw(ArgumentError("ppBigiel2010!: `quantity` can only be :molecular, :neutral or \
        :atomic, but I got :$(quantity)."))

    end

    # Delete missing data
    idxs = map(ismissing, Σg) ∪ map(ismissing, Σsfr)

    deleteat!(Σg, idxs)
    deleteat!(Σsfr, idxs)

    # Set the correct scale and units
    if x_log
        Σg = log10.(ustrip.(x_unit, Σg))
    else
        Σg = ustrip.(x_unit, Σg)
    end

    if y_log
        Σsfr = log10.(ustrip.(y_unit, Σsfr))
    else
        Σsfr = ustrip.(y_unit, Σsfr)
    end

    ################################################################################################
    # Plot the galactic data
    ################################################################################################

    sp = scatter!(figure.current_axis.x, Σg, Σsfr; color=(color, 0.5), marker=:star4, markersize=10)

    translate!(Accum, sp, 0, 0, -10)

    if galaxy == :all
        label = "Bigiel et al. 2010"
    else
        label = "$(galaxy) - Bigiel et al. 2010"
    end

    return ([MarkerElement(; color, marker=:star4)], [label])

end

@doc raw"""
    ppSun2023!(
	    figure::Makie.Figure;
        <keyword arguments>
    )::Tuple{Vector{<:LegendElement},Vector{AbstractString}}

Draw a scatter plot of the SFR surface density vs molecular surface density (molecular Kennicutt-Schmidt law) for a given galaxy, using the data from Sun et al. (2023).

!!! note
    The resolution used in Sun et al. (2023) is 1.5 kpc (see Section 2).

# Arguments

  - `figure::Makie.Figure`: Makie figure.
  - `galaxy::Union{String,Symbol}=:main`: Target galaxy. The options are:

      + One of the 80 galaxies in the dataset, e.g. "ESO097-013", "IC1954", "IC5273", "NGC1546", "NGC1559", "NGC1566", etc. For a full list see the reference below.
      + :main: Every galaxy with $-2.0 < \log_{10}(t_\mathrm{dep} \, / \, \mathrm{Gyr}) < 2.0$, where $t_\mathrm{dep} = \Sigma_\mathrm{H_2} / \Sigma_\mathrm{SFR}$ is the depletion time.
      + :all: All 80 galaxies in the dataset. For more information on each galaxy see Sun et al. (2023).
  - `sfr_calibration::Symbol=:Halpha`: SFR calibration for combining UV, optical, and/or IR data. The options are: `:Halpha`, `:FUV`, and `:AV_corrected_Halpha`. For an explanation of each one see section 2 of Sun et al. (2023).
  - `h2_prescription::Symbol=:S20`: Prescription for the CO-to-H₂ conversion factor. The options are: `:S20`, `:Mw`, `:B13`, and `:G20`. For an explanation of each one see section 2 of Sun et al. (2023).
  - `x_unit::Unitful.Units=u"Msun * pc^-2"`: Unit for the area density of gas used in `figure`.
  - `y_unit::Unitful.Units=u"Msun * yr^-1 * kpc^-2"`: Unit for the area density of star formation rate used in `figure`.
  - `x_log::Bool=true`: If the x axis will be plotted as the ``\\log_{10}`` of the gas surface density.
  - `y_log::Bool=true`: If the y axis will be plotted as the ``\\log_{10}`` of the SFR surface density.
  - `color::ColorType=WONG_RED`: Color of the markers.

# Returns

  - A tuple with the elements for the legend:

      + A `MarkerElement` to be used as the marker.
      + The label.

# References

J. Sun et al. (2023). *Star Formation Laws and Efficiencies across 80 Nearby Galaxies*. The Astrophysical Journal Letters, **945(2)**, L19. [doi:10.3847/2041-8213/acbd9c](https://doi.org/10.3847/2041-8213/acbd9c)
"""
function ppSun2023!(
    figure::Makie.Figure;
    galaxy::Union{String,Symbol}=:main,
    sfr_calibration::Symbol=:Halpha,
    h2_prescription::Symbol=:S20,
    x_unit::Unitful.Units=u"Msun * pc^-2",
    y_unit::Unitful.Units=u"Msun * yr^-1 *  kpc^-2",
    x_log::Bool=true,
    y_log::Bool=true,
    color::ColorType=WONG_RED,
)::Tuple{Vector{<:LegendElement},Vector{AbstractString}}

    ################################################################################################
    # Load table A1 from Sun et al. 2023
    ################################################################################################

    raw_data = CSV.read(
        SUN2023_TABLE,
        DataFrame;
        delim=' ',
        skipto=58,
        header=false,
        ignorerepeated=true,
        silencewarnings=true,
    )

    clean_data = ifelse.(raw_data .=== missing, Inf, raw_data)

    # Shift in the column indices correspondig to each SFR calibration
    sfr_calibrations = Dict(:Halpha => 0, :FUV => 2, :AV_corrected_Halpha => 4)

    # Shift in the column indices correspondig to each CO-to-H₂ prescription
    h2_prescriptions = Dict(:S20 => 0, :Mw => 2, :B13 => 4, :G20 => 6)

    # List of available galaxies
    galaxies = unique(raw_data[:, 1])

    ################################################################################################
    # Find the target galaxy and read its values
    ################################################################################################

    if isa(galaxy, String)

        (
            galaxy ∈ galaxies ||
            throw(ArgumentError("ppSun2023!: `galaxy` = $(galaxy) is not a valid galaxy"))
        )

        data = filter(:Column1 => isequal(galaxy), clean_data)

    else

        (
            galaxy ∈ [:main, :all] ||
            throw(ArgumentError("ppSun2023!: `galaxy` can only be :main or :all but I got \
            :$(galaxy)"))
        )

        data = clean_data

    end

    Σsfr = data[:, 4 + sfr_calibrations[sfr_calibration]] .* u"Msun * yr^-1 * kpc^-2"
    Σh2  = data[:, 10 + h2_prescriptions[h2_prescription]] .* u"Msun * pc^-2"

    # Set the correct scale and units
    if x_log
        x_data = @. log10(ustrip(x_unit, Σh2))
    else
        x_data = @. ustrip(x_unit, Σh2)
    end

    if y_log
        y_data = @. log10(ustrip(y_unit, Σsfr))
    else
        y_data = @. ustrip(y_unit, Σsfr)
    end

    # Delete missing data
    filter = x -> isnan(x) || isinf(x)
    idxs   = map(filter, x_data) ∪ map(filter, y_data)

    if galaxy == :main
        # Compute the depletion time
        tdep = @. log10(ustrip(u"Gyr", Σh2 ./ Σsfr))

        # Filter galaxies with a depletion time outside the range [-2.0, 2.0]
        # to reproduce Fig. 1 of Sun et al. (2023)
        tdep_filter = x -> isnan(x) || isinf(x) || x < -2.0 || x > 2.0

        idxs = idxs ∪ map(tdep_filter, tdep)
    end

    deleteat!(x_data, idxs)
    deleteat!(y_data, idxs)

    ################################################################################################
    # Plot the galactic data
    ################################################################################################

    sp = scatter!(
        figure.current_axis.x,
        x_data,
        y_data;
        color=(color, 0.5),
        marker=:star4,
        markersize=10,
    )

    # Put the post processing elements at the back of the plot
    translate!(Accum, sp, 0, 0, -10)

    if isa(galaxy, String)
        label = "$(galaxy) - Sun et al. (2023)"
    else
        label = "Sun et al. (2023)"
    end

    return ([MarkerElement(; color, marker=:star4)], [label])

end

"""
    ppLeroy2008!(
        figure::Makie.Figure;
        <keyword arguments>
    )::Tuple{Vector{<:LegendElement},Vector{AbstractString}}

Draw a scatter plot of the SFR surface density vs gas surface density (Kennicutt-Schmidt law) for a given galaxy, using the data from Leroy et al. (2008).

!!! note

    The resolution used in Leroy et al. (2008) is 800 pc for spirals and 400 pc for dwarf galaxies (Section 3.1).

# Arguments

  - `figure::Makie.Figure`: Makie figure.
  - `quantity::Symbol=:molecular`: Gas quantity for the x axis. The options are:

      + `:molecular` -> Surface density of molecular gas.
      + `:neutral`   -> Surface density of neutral gas.
      + `:atomic`    -> Surface density of atomic gas.
  - `galaxies::Vector=["MW"]`: Target galaxies. The options are:

      + One of the 23 galaxies in Leroy et al. (2008), e.g. "DDO154", "HOI", "HOII", "IC2574", "NGC0628", "NGC0925", etc. For a full list see Leroy et al. (2008).
      + :all: All 23 galaxies in Leroy et al. (2008).
  - `x_unit::Unitful.Units=u"Msun * pc^-2"`: Unit for the area density of gas used in `figure`.
  - `y_unit::Unitful.Units=u"Msun * yr^-1 * kpc^-2"`: Unit for the area density of star formation rate used in `figure`.
  - `x_log::Bool=true`: If the x axis will be plotted as the ``\\log_{10}`` of the gas surface density.
  - `y_log::Bool=true`: If the y axis will be plotted as the ``\\log_{10}`` of the SFR surface density.
  - `error_bars::Bool=false`: If error bars will be plotted.
  - `color::ColorType=WONG_RED,`: Color for the markers.

# Returns

  - A tuple with the elements for the legend:

      + A `MarkerElement` to be used as the marker.
      + The label.

# References

A. K. Leroy et al. (2008). *THE STAR FORMATION EFFICIENCY IN NEARBY GALAXIES: MEASURING WHERE GAS FORMS STARS EFFECTIVELY*. The Astronomical Journal **136(6)**, 2782–2845. [doi:10.1088/0004-6256/136/6/2782](https://doi.org/10.1088/0004-6256/136/6/2782)
"""
function ppLeroy2008!(
    figure::Makie.Figure;
    quantity::Symbol=:molecular,
    galaxies::Vector{<:Union{Symbol,String}}=[:all],
    x_unit::Unitful.Units=u"Msun * pc^-2",
    y_unit::Unitful.Units=u"Msun * kpc^-2 * yr^-1",
    x_log::Bool=true,
    y_log::Bool=true,
    error_bars::Bool=false,
    color::ColorType=WONG_RED,
)::Tuple{Vector{<:LegendElement},Vector{AbstractString}}

    ################################################################################################
    # Load the data from Leroy et al. (2008)
    ################################################################################################

    leroy2008 = load(LEROY2008_DATA_PATH)["dataframe"]

    # List of available galaxies in Leroy et al. (2008)
    leroy_galaxies = unique(leroy2008[!, "Name"])

    ################################################################################################
    # Find the target galaxies and read their values
    ################################################################################################

    legend_elements = Vector{LegendElement}(undef, length(galaxies))
    labels = Vector{String}(undef, length(galaxies))

    for galaxy in galaxies

        if galaxy ∈ leroy_galaxies

            leroy_data = filter(:Name => isequal(galaxy), leroy2008)

        elseif galaxy == :all

            leroy_data = leroy2008

        else

            throw(ArgumentError("ppLeroy2008!: $(galaxy) is not a valid argument for a galaxy"))

        end

        Σsfr = leroy_data[!, "FUV+24"] .± leroy_data[!, "e_FUV+24"]
        ΣH2  = leroy_data[!, "SigmaH2"] .± leroy_data[!, "e_SigmaH2"]
        ΣHI  = leroy_data[!, "SigmaHI"] .± leroy_data[!, "e_SigmaHI"]

        if quantity == :molecular

            Σgas = ΣH2

        elseif quantity == :atomic

            Σgas = ΣHI

        elseif quantity == :neutral

            Σgas = ΣH2 .+ ΣHI

        else

            throw(ArgumentError("ppLeroy2008!: `quantity` can only be :molecular, :atomic or \
            :neutral, but I got :$(quantity)"))

        end

        # Set the correct scale and units for the x axis
        if x_log
            x_data  = @. Measurements.value(log10(ustrip(x_unit, Σgas)))
            x_error = @. Measurements.uncertainty(log10(ustrip(x_unit, Σgas)))
        else
            x_data  = @. Measurements.value(ustrip(x_unit, Σgas))
            x_error = @. Measurements.uncertainty(ustrip(x_unit, Σgas))
        end

        # Set the correct scale and units for the y axis
        if y_log
            y_data  = @. Measurements.value(log10(ustrip(y_unit, Σsfr)))
            y_error = @. Measurements.uncertainty(log10(ustrip(y_unit, Σsfr)))
        else
            y_data  = @. Measurements.value(ustrip(y_unit, Σsfr))
            y_error = @. Measurements.uncertainty(ustrip(y_unit, Σsfr))
        end

        ############################################################################################
        # Plot the galactic data
        ############################################################################################

        sp = scatter!(
            figure.current_axis.x,
            x_data,
            y_data;
            color=(color, 0.7),
            marker=:star4,
            markersize=10,
        )

        translate!(Accum, sp, 0, 0, -10)

        if error_bars

            xerr = errorbars!(
                figure.current_axis.x,
                x_data,
                y_data,
                x_error,
                direction=:x;
                color=(color, 0.7),
            )

            yerr = errorbars!(figure.current_axis.x, x_data, y_data, y_error; color=(color, 0.7))

            translate!(Accum, xerr, 0, 0, -10)
            translate!(Accum, yerr, 0, 0, -10)

        end

    end

    return [MarkerElement(; color, marker=:star4)], ["Leroy et al. (2008)"]

end

"""
    ppMolla2015!(
        figure::Makie.Figure,
        quantity::Symbol,
        <keyword arguments>
    )::Tuple{Vector{<:LegendElement},Vector{AbstractString}}

Draw a profile for the Milky Way using the data compiled by Mollá et al. (2015).

# Arguments

  - `figure::Makie.Figure`: Makie figure.
  - `quantity::Symbol`: Quantity for the y axis. The options are:

      + `:stellar_area_density`   -> Stellar mass surface density, as ``\\log10(\\Sigma_\\star)``.
      + `:molecular_area_density` -> Molecular mass surface density, as ``\\log10(\\Sigma_\\text{H2})``.
      + `:atomic_area_density`    -> Atomic mass surface density, as ``\\log10(\\Sigma_\\text{HI})``.
      + `:sfr_area_density`       -> Star formation rate surface density, as ``\\log10(\\Sigma_\\text{SFR})``.
      + `:X_stellar_abundance`    -> Stellar abundance of element ``\\mathrm{X}``, as [`ABUNDANCE_SHIFT`](@ref) + ``\\log_{10}(\\mathrm{X \\, / \\, H})``. ``\\mathrm{X}`` can be O (oxygen), N (nitrogen), or C (carbon).
  - `y_unit::Unitful.Units=Unitful.NoUnits`: Target unit for `quantity`.
  - `color::ColorType=WONG_RED`: Color of the line.
  - `linestyle::LineStyleType=:solid`: Style of the line.
  - `error_bars::Bool=true`: If the error bars will be plotted.

# Returns

  - A tuple with the elements for the legend:

      + A `MarkerElement` to be used as the marker.
      + The label.

# References

M. Mollá et al. (2015). *Galactic chemical evolution: stellar yields and the initial mass function*. Monthly Notices of the Royal Astronomical Society **451(4)**, 3693–3708. [doi:10.1093/mnras/stv1102](https://doi.org/10.1093/mnras/stv1102)
"""
function ppMolla2015!(
    figure::Makie.Figure,
    quantity::Symbol;
    y_unit::Unitful.Units=Unitful.NoUnits,
    color::ColorType=WONG_RED,
    linestyle::LineStyleType=:solid,
    error_bars::Bool=true,
)::Tuple{Vector{<:LegendElement},Vector{AbstractString}}

    ################################################################################################
    # Load the data from Mollá et al. (2015)
    ################################################################################################

    raw = CSV.read(
        MOLLA2015_DATA_PATH,
        DataFrame;
        header=[
            "R",
            "ΣHI",
            "ΣHI error",
            "ΣH2",
            "ΣH2 error",
            "logΣ*",
            "logΣ* error",
            "logΣsfr",
            "logΣsfr error",
            "C/H",
            "ΔC/H",
            "N/H",
            "ΔN/H",
            "O/H",
            "ΔO/H",
        ],
    )

    x_values = raw[!, "R"]

    ################################################################################################
    # Select the quantity for the y axis
    ################################################################################################

    if quantity == :stellar_area_density
        # M⊙ pc^-2
        factor = log10(ustrip(y_unit, 1.0u"Msun * pc^-2"))
        y_data = (raw[!, "logΣ*"] .± raw[!, "logΣ* error"]) .+ factor
    elseif quantity ∈ [:ode_molecular_area_density, :br_molecular_area_density]
        # M⊙ pc^-2
        factor = ustrip(y_unit, 1.0u"Msun * pc^-2")
        y_data = log10.((raw[!, "ΣH2"] .± raw[!, "ΣH2 error"]) .* factor)
    elseif quantity ∈ [:ode_atomic_area_density, :br_atomic_area_density]
        # M⊙ pc^-2
        factor = ustrip(y_unit, 1.0u"Msun * pc^-2")
        y_data = log10.((raw[!, "ΣHI"] .± raw[!, "ΣHI error"]) .* factor)
    elseif quantity == :sfr_area_density
        # M⊙ pc^-2 Gyr^-1
        factor = log10(ustrip(y_unit, 1.0u"Msun * pc^-2 * Gyr^-1"))
        y_data = (raw[!, "logΣsfr"] .± raw[!, "logΣsfr error"]) .+ factor
    elseif quantity == :O_stellar_abundance
        # dimensionless
        y_data = @. (raw[!, "O/H"] - 12.0 + ABUNDANCE_SHIFT[:O]) ± raw[!, "ΔO/H"]
    elseif quantity == :N_stellar_abundance
        # dimensionless
        y_data = @. (raw[!, "N/H"] - 12.0 + ABUNDANCE_SHIFT[:N]) ± raw[!, "ΔN/H"]
    elseif quantity == :C_stellar_abundance
        # dimensionless
        y_data = @. (raw[!, "C/H"] - 12.0 + ABUNDANCE_SHIFT[:C]) ± raw[!, "ΔC/H"]
    else
        throw(ArgumentError("ppMolla2015: `x_quantity` can only be  :stellar_area_density, \
        :ode_molecular_area_density, :br_molecular_area_density, :ode_atomic_area_density, \
        :br_atomic_area_density, :sfr_area_density, :O_stellar_abundance, :N_stellar_abundance or :C_stellar_abundance, but I got :$(quantity)"))
    end

    y_values = Measurements.value.(y_data)
    y_uncertainties = Measurements.uncertainty.(y_data)

    ################################################################################################
    # Plot the galactic profiles
    ################################################################################################

    # Plot the mean values
    slp = scatterlines!(figure.current_axis.x, x_values, y_values; color, linestyle, marker=:utriangle)

    translate!(Accum, slp, 0, 0, -10)

    if error_bars
        # Plot the error bars
        ep = errorbars!(figure.current_axis.x, x_values, y_values, y_uncertainties; color)
        translate!(Accum, ep, 0, 0, -10)
    end

    return [MarkerElement(; color, marker=:utriangle)], ["Mollá et al. (2015)"]

end

"""
    ppAgertz2021!(
        figure::Makie.Figure;
        <keyword arguments>
    )::Tuple{Vector{<:LegendElement},Vector{AbstractString}}

Draw the stellar density profiles from McMillan (2011) and Leroy et al. (2008), show in Agertz et al. (2021).

# Arguments

  - `figure::Makie.Figure`: Makie figure.
  - `galaxies::Vector=["MW"]`: Target galaxies. The options are:

      + One of the 23 galaxies in Leroy et al. (2008), e.g. "DDO154", "HOI", "HOII", "IC2574", "NGC0628", "NGC0925", etc. For a full list see Agertz et al. (2021).
      + :all: All 23 galaxies in the dataset of Leroy et al. (2008) (this options will plot the data as an scatter plot with transparency instead of a line plot).
      + "MW": The Milky Way fits from McMillan (2011).
  - `x_unit::Unitful.Units=u"kpc"`: Unit for the x axis in `figure`.
  - `y_unit::Unitful.Units=u"Msun * kpc^-2"`: Unit for the y axis in `figure`.
  - `x_log::Bool=false`: If the x axis will be plotted as the ``\\log_{10}`` of the galactocentric radius.
  - `y_log::Bool=true`: If the y axis will be plotted as the ``\\log_{10}`` of the stellar surface density.
  - `error_band::Bool=true`: If the error band will be plotted.
  - `colors::Vector{<:ColorType}=[WONG_RED]`: Colors for stellar profiles.
  - `linestyle::LineStyleType=:solid`: Style of the line.
  - `linewidth::Int=3`: Width of the lines.

# Returns

  - A tuple with the elements for the legend:

      + A `LineElement` to be used as the marker.
      + The label.

# References

O. Agertz et al. (2021). *VINTERGATAN – I. The origins of chemically, kinematically, and structurally distinct discs in a simulated Milky Way-mass galaxy*. Monthly Notices of the Royal Astronomical Society. **503(4), 5826–5845. [doi:10.1093/mnras/stab322](https://doi.org/10.1093/mnras/stab322)

A. K. Leroy et al. (2008). *THE STAR FORMATION EFFICIENCY IN NEARBY GALAXIES: MEASURING WHERE GAS FORMS STARS EFFECTIVELY*. The Astronomical Journal **136(6)**, 2782–2845. [doi:10.1088/0004-6256/136/6/2782](https://doi.org/10.1088/0004-6256/136/6/2782)

P. J. McMillan (2011). *Mass models of the Milky Way*. Monthly Notices of the Royal Astronomical Society **414(3)**, 2446–2457. [doi:10.1111/j.1365-2966.2011.18564.x](https://doi.org/10.1111/j.1365-2966.2011.18564.x)
"""
function ppAgertz2021!(
    figure::Makie.Figure;
    galaxies::Vector=["MW"],
    x_unit::Unitful.Units=u"kpc",
    y_unit::Unitful.Units=u"Msun * kpc^-2",
    x_log::Bool=false,
    y_log::Bool=true,
    error_band::Bool=true,
    colors::Vector{<:ColorType}=[WONG_RED],
    linestyle::LineStyleType=:solid,
    linewidth::Int=3,
)::Tuple{Vector{<:LegendElement},Vector{AbstractString}}

    ################################################################################################
    # Load the data from McMillan (2011) and Leroy et al. (2008)
    ################################################################################################

    mcmillan2011 = load(MCMILLAN2011_DATA_PATH)["dataframe"]
    leroy2008 = load(LEROY2008_DATA_PATH)["dataframe"]

    # List of available galaxies in Leroy et al. (2008)
    leroy_galaxies = unique(leroy2008[!, "Name"])

    ################################################################################################
    # Find the target galaxies and read its values
    ################################################################################################

    legend_elements = Vector{LegendElement}(undef, length(galaxies))
    labels = Vector{String}(undef, length(galaxies))

    for (i, galaxy) in pairs(galaxies)

        if galaxy ∈ leroy_galaxies

            leroy_data = filter(:Name => isequal(galaxy), leroy2008)
            r  = leroy_data[!, "Rad"]
            Σs = leroy_data[!, "Sigma*"] .± leroy_data[!, "e_Sigma*"]

        elseif galaxy == :all

            r  = leroy2008[!, "Rad"]
            Σs = leroy2008[!, "Sigma*"] .± leroy2008[!, "e_Sigma*"]

        elseif galaxy == "MW"

            r  = mcmillan2011[!, :r]
            Σs = mcmillan2011[!, :Sigma]

        else

            throw(ArgumentError("ppAgertz2021!: $(galaxy) is not a valid argument for a galaxy"))

        end

        # Set the correct scale and units for the x axis
        if x_log
            x_data = @. log10(ustrip(x_unit, r))
        else
            x_data = @. ustrip(x_unit, r)
        end

        # Set the correct scale and units for the y axis
        if y_log
            y_data  = @. Measurements.value(log10(ustrip(y_unit, Σs)))
            y_error = @. Measurements.uncertainty(log10(ustrip(y_unit, Σs)))
        else
            y_data  = @. Measurements.value(ustrip(y_unit, Σs))
            y_error = @. Measurements.uncertainty(ustrip(y_unit, Σs))
        end

        color = ring(colors, i)

        ############################################################################################
        # Plot the stellar profiles
        ############################################################################################

        if galaxy == :all

            sp = scatter!(figure.current_axis.x, x_data, y_data; color=(color, 0.4))

            # Put the post processing elements at the back of the plot
            translate!(Accum, sp, 0, 0, -10)

            legend_elements[i] = MarkerElement(; color, marker=:circle)
            labels[i] = "Agertz et al. (2021)"

        else

            if error_band
                bp = band!(
                    figure.current_axis.x,
                    x_data,
                    y_data .- y_error,
                    y_data .+ y_error;
                    color,
                )
            end

            lp = lines!(figure.current_axis.x, x_data, y_data; color, linestyle, linewidth)

            # Put the post processing elements at the back of the plot
            if error_band
                translate!(Accum, bp, 0, 0, -10)
            end
            translate!(Accum, lp, 0, 0, -10)

            legend_elements[i] = LineElement(; color, linestyle, linewidth)
            labels[i] = "$(galaxy) - Agertz et al. (2021)"

        end

    end

    return legend_elements, labels

end

"""
    ppFeldmann2020!(
        figure::Makie.Figure,
        x_quantity::Symbol,
        y_quantity::Symbol;
        <keyword arguments>
    )::Tuple{Vector{<:LegendElement},Vector{AbstractString}}

Draw the experimental data from xGASS and xCOLD GASS, compiled by Feldmann (2020), as a line or scatter plot.

# Arguments

  - `figure::Makie.Figure`: Makie figure.
  - `x_quantity::Symbol`: Quantity for the x axis. The options are:

      + `:stellar`   -> Stellar mass.
      + `:molecular` -> Molecular mass.
      + `:atomic`    -> Atomic mass.
      + `:sfr`       -> Star formation rate of the last `AGE_RESOLUTION`.
  - `y_quantity::Symbol`: Quantity for the y axis. The options are:

      + `:stellar`   -> Stellar mass.
      + `:molecular` -> Molecular mass.
      + `:atomic`    -> Atomic mass.
      + `:sfr`       -> Star formation rate of the last `AGE_RESOLUTION`.
  - `scatter::Bool=false`: If the data will be presented as a line plot with error bands (default), or alternatively, a scatter plot.
  - `xlog::Bool=true`: If the x axis will be plotted as the ``\\log_{10}`` of the x quantity.
  - `ylog::Bool=true`: If the y axis will be plotted as the ``\\log_{10}`` of the y quantity.

# Returns

  - A tuple with the elements for the legend:

      + A `MarkerElement` to be used as the marker.
      + The label.

# References

R. Feldmann (2020). *The link between star formation and gas in nearby galaxies*. Communications Physics **3(226)**. [doi:10.1038/s42005-020-00493-0](https://doi.org/10.1038/s42005-020-00493-0)
"""
function ppFeldmann2020!(
    figure::Makie.Figure,
    x_quantity::Symbol,
    y_quantity::Symbol;
    scatter::Bool=false,
    xlog::Bool=true,
    ylog::Bool=true,
)::Tuple{Vector{<:LegendElement},Vector{AbstractString}}

    ################################################################################################
    # Load the data from McMillan (2011) and Feldmann (2020)
    ################################################################################################

    raw = CSV.read(FELDMANN2020_DATA_PATH, DataFrame, comment="#")

    ################################################################################################
    # Select the quantity for the x axis
    ################################################################################################

    if x_quantity == :stellar

        x_data = exp10.(raw[!, "lgMstar"])

    elseif x_quantity == :molecular

        x_data = raw[!, "MH2"]

    elseif x_quantity == :atomic

        x_data = raw[!, "MHI"]

    elseif x_quantity == :sfr

        x_data = raw[!, "SFR"]

    else

        throw(ArgumentError("ppFeldmann2020!: `x_quantity` can only be :stellar, \
        :molecular, :atomic or :sfr, but I got :$(x_quantity)"))

    end

    ################################################################################################
    # Select the quantity for the y axis, with its uncertainty
    ################################################################################################

    if y_quantity == :stellar

        # Compute the mean "error" for the stellar mass
        err_low  = raw[!, "lgMstar"] .- raw[!, "lgMstar_p16"]
        err_high = raw[!, "lgMstar_p84"] .- raw[!, "lgMstar"]
        err_mean = @. (err_low + err_high) / 2.0

        y_data = exp10.(raw[!, "lgMstar"] .± err_mean)

    elseif y_quantity == :molecular

        y_data = raw[!, "MH2"] .± raw[!, "e_MH2"]

    elseif y_quantity == :atomic

        y_data = raw[!, "MHI"] .± raw[!, "e_MHI"]

    elseif y_quantity == :sfr

        y_data = raw[!, "SFR"] .± raw[!, "e_SFR"]

    else

        throw(ArgumentError("ppFeldmann2020!: `y_quantity` can only be :stellar_mass, \
        :molecular_mass, :br_molecular_mass, :atomic_mass, or :observational_sfr, \
        but I got :$(y_quantity)"))

    end

    ################################################################################################
    # Delete missing data
    ################################################################################################

    missing_idxs = map(ismissing, x_data) .| map(ismissing, y_data)
    deleteat!(x_data, missing_idxs)
    deleteat!(y_data, missing_idxs)

    ################################################################################################
    # Get the scaling of each axis
    ################################################################################################

    if xlog
        x_scale = log10
    else
        x_scale = identity
    end

    if ylog
        y_scale = log10
    else
        y_scale = identity
    end

    ################################################################################################
    # Plot the galactic data as a scatter plot
    ################################################################################################

    if scatter

        y_mean_values = Measurements.value.(y_data)

        sp = scatter!(
            figure.current_axis.x,
            x_scale.(x_data),
            y_scale.(y_mean_values);
            color=WONG_RED,
            marker=:star4,
            markersize=12,
        )

        translate!(Accum, sp, 0, 0, -10)

        return [MarkerElement(; color=WONG_RED, marker=:star4)], ["Feldmann (2020)"]

    end

    ################################################################################################
    # Plot the galactic data as a line plot
    ################################################################################################

    # Scale the raw data
    scaled_xs = x_scale.(x_data)
    scaled_ys = y_scale.(y_data)

    # Set up the grid
    x_min, x_max = extrema(scaled_xs)
    grid         = CircularGrid(x_max, 50; shift=x_min)
    n_bins       = length(grid.grid)
    bin_width    = (x_max - x_min) / n_bins

    histogram = [Int[] for _ in 1:n_bins]

    # Compute the indices that fall within each bin,
    # ignoring missings and values outside the grid range
    for (x_idx, x_value) in pairs(scaled_xs)

        if x_value < x_min || x_max < x_value
            continue
        elseif x_value == x_min
            hist_idx = 1
        elseif x_value == x_max
            hist_idx = n_bins
        else
            hist_idx = ceil(Int, (x_value - x_min) / bin_width)
        end

        push!(histogram[hist_idx], x_idx)

    end

    ################################################################################################
    # Values for the y axis
    ################################################################################################

    y_axis = Vector{Union{Measurement{Float64},Missing}}(undef, n_bins)

    for (i, hist_idxs) in pairs(histogram)

        # For each bin select the corresponding y measurements
        binned_y_measurements = scaled_ys[hist_idxs]

        if isempty(binned_y_measurements)

            y_axis[i] = missing

        else

            # Compute the weighted mean of all the y measurements in this bin
            weighted_mean = weightedmean(binned_y_measurements)

            # Compute the mean
            y_mean = Measurements.value(weighted_mean)

            # Compute the uncertainty
            spread = std(Measurements.value.(binned_y_measurements); corrected=false)
            stderr = Measurements.uncertainty(weighted_mean)
            y_err  = sqrt(spread^2 + stderr^2)

            y_axis[i] = y_mean ± y_err

        end

    end

    # Find the empty bins
    missing_idxs = map(ismissing, y_axis)

    # Delete the empty bins
    deleteat!(y_axis, missing_idxs)

    ################################################################################################
    # Values for the x axis
    ################################################################################################

    x_axis = grid.grid

    # Delete the empty bins
    deleteat!(x_axis, missing_idxs)

    ################################################################################################
    # Line and band plots
    ################################################################################################

    y_axis_value       = Measurements.value.(y_axis)
    y_axis_uncertainty = Measurements.uncertainty.(y_axis)

    # Plot the central line
    lp = lines!(
        figure.current_axis.x,
        x_axis,
        y_axis_value;
        color=(WONG_RED, 0.6),
    )

    # Construct the 1σ band
    yupper_1σ = y_axis_value .+ y_axis_uncertainty
    ylower_1σ = y_axis_value .- y_axis_uncertainty

    # Plot the 1σ band
    bp1 = band!(figure.current_axis.x, x_axis, ylower_1σ, yupper_1σ; color=WONG_RED)

    # Construct the 2σ upper band
    yupper_2σ = y_axis_value .+ 2 * y_axis_uncertainty

    # Plot the 2σ upper band
    bp2u = band!(figure.current_axis.x, x_axis, yupper_1σ, yupper_2σ; color=WONG_ORANGE)

    # Construct the 2σ lower band
    ylower_2σ = y_axis_value .- 2 * y_axis_uncertainty

    # Plot the 2σ lower band
    bp2l = band!(figure.current_axis.x, x_axis, ylower_2σ, ylower_1σ; color=WONG_ORANGE)

    # Put the post processing elements at the back of the plot
    translate!(Accum, bp1, 0, 0, -10)
    translate!(Accum, bp2u, 0, 0, -10)
    translate!(Accum, bp2l, 0, 0, -10)
    translate!(Accum, lp, 0, 0, -10)

    return (
        [PolyElement(; color=(WONG_RED, 0.5)), PolyElement(; color=(WONG_ORANGE, 0.5))],
        [
            L"\mathrm{Feldmann \,\, (2020)} \,\, 1\sigma",
            L"\mathrm{Feldmann \,\, (2020)} \,\, 2\sigma",
        ],
    )

end
