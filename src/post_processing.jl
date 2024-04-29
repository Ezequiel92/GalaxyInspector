####################################################################################################
# Post-processing functions.
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
    )::Nothing

Draw vertical lines.

# Arguments

  - `figure::Makie.Figure`: Makie figure to be drawn over.
  - `positions::Vector{<:Real}`: The x coordinates of the lines.
  - `colors::Vector{<:ColorType}=[:red]`: Colors of the lines.
  - `line_styles::Vector{<:LineStyleType}=[nothing]`: Styles of the lines. `nothing` will produce a solid line.
  - `warnings::Bool=true`: If a warning will be raised when all the lines are outside the plot range.
"""
function ppVerticalFlags!(
    figure::Makie.Figure,
    positions::Vector{<:Real};
    colors::Vector{<:ColorType}=[:red],
    line_styles::Vector{<:LineStyleType}=[nothing],
    warnings::Bool=true,
)::Nothing

    # Filter out values outside the range of the original plot
    rangeCut!(positions, xlimits!(figure); keep_edges=false)

    if isempty(positions)

        # Don't draw anything if all the flags are outside the original plot window
        !warnings || @warn("ppVerticalFlags!: All vertical lines lie outside the plot range")

    else

        # Draw the vertical lines
        @inbounds for (i, position) in pairs(positions)
            color = ring(colors, i)
            linestyle = ring(line_styles, i)
            vlines!(figure.current_axis.x, position; color, linestyle)
        end

    end

    return nothing

end

"""
    ppHorizontalFlags!(
        figure::Makie.Figure,
        positions::Vector{<:Real};
        <keyword arguments>
    )::Nothing

Draw horizontal lines.

# Arguments

  - `figure::Makie.Figure`: Makie figure to be drawn over.
  - `positions::Vector{<:Real}`: The y coordinates of the lines.
  - `colors::Vector{<:ColorType}=[:red]`: Colors of the lines.
  - `line_styles::Vector{<:LineStyleType}=[nothing]`: Styles of the lines. `nothing` will produce a solid line.
  - `warnings::Bool=true`: If a warning will be raised when all the lines are outside the plot range.
"""
function ppHorizontalFlags!(
    figure::Makie.Figure,
    positions::Vector{<:Real};
    colors::Vector{<:ColorType}=[:red],
    line_styles::Vector{<:LineStyleType}=[nothing],
    warnings::Bool=true,
)::Nothing

    # Filter out values outside the range of the original plot
    rangeCut!(positions, ylimits!(figure); keep_edges=false)

    if isempty(positions)

        # Don't draw anything if all the flags are outside the original plot window
        !warnings || @warn("ppHorizontalFlags!: All horizontal lines lie outside the plot range")

    else

        # Draw the horizontal lines
        @inbounds for (i, position) in pairs(positions)
            color = ring(colors, i)
            linestyle = ring(line_styles, i)
            hlines!(figure.current_axis.x, position; color, linestyle)
        end

    end

    return nothing

end

"""
    ppCross!(
        figure::Makie.Figure,
        cross_point::Tuple{<:Real,<:Real};
        <keyword arguments>
    )::Nothing

Draw two lines, one horizontal and one vertical.

# Arguments

  - `figure::Makie.Figure`: Makie figure to be drawn over.
  - `cross_point::Tuple{<:Real,<:Real}`: Crossing point of the lines.
  - `color::ColorType=:red`: Color of the lines.
  - `linestyle::LineStyleType=nothing`: Style of the lines. `nothing` will produce a solid line.
  - `warnings::Bool=true`: If a warning will be raised when at least one of the lines is outside the plot range.
"""
function ppCross!(
    figure::Makie.Figure,
    cross_point::Tuple{<:Real,<:Real};
    color::ColorType=:red,
    linestyle::LineStyleType=nothing,
    warnings::Bool=true,
)::Nothing

    x_limits = xlimits!(figure)
    y_limits = ylimits!(figure)

    if x_limits[1] < cross_point[1] < x_limits[2]
        # Draw the vertical line
        vlines!(figure.current_axis.x, cross_point[1]; color, linestyle)
    else
        !warnings || @warn("ppCross!: The vertical line lies outside the plot range")
    end

    if y_limits[1] < cross_point[2] < y_limits[2]
        # Draw the horizontal line
        hlines!(figure.current_axis.x, cross_point[2]; color, linestyle)
    else
        !warnings || @warn("ppCross!: The horizontal line lies outside the plot range")
    end

    return nothing

end

"""
    ppAnnotation!(
        figure::Makie.Figure,
        text::String;
        <keyword arguments>
    )::Nothing

Add an annotation to the plot.

# Arguments

  - `figure::Makie.Figure`: Makie figure to be drawn over.
  - `text::String`: Text to be written.
  - `color::ColorType=:black`: Color of the text.
  - `rel_pos::Union{NTuple{2,Real},Nothing}=nothing`: Relative position for the top left corner of the text box, within the plot. If set to `nothing`, the text will be printed at the top left corner of the plot.
"""
function ppAnnotation!(
    figure::Makie.Figure,
    text::String;
    color::ColorType=:black,
    rel_pos::Union{NTuple{2,Real},Nothing}=nothing,
)::Nothing

    if isnothing(rel_pos)
        pos = absCoor(figure, 0.03, 0.15)
    else
        if rel_pos[1] < 0.0 || rel_pos[2] < 0.0 || 1.0 < rel_pos[1] || 1.0 < rel_pos[2]
            throw(ArgumentError("ppAnnotation!: The values in `rel_pos` should be between 0 and 1"))
        else
            pos = absCoor(figure, rel_pos...)
        end
    end

    text!(
        figure.current_axis.x,
        pos[1],
        pos[2];
        text,
        color,
        align=(:left, :top),
        fontsize=40,
        font=:bold,
    )

    return nothing

end

@doc raw"""
    ppFitLine!(
        figure::Makie.Figure;
        <keyword arguments>
    )::Union{Tuple{Vector{<:LegendElement},Vector{AbstractString}},Nothing}

Draw a linear fit for the data in `figure`.

An annotation with the equation $y = a \, x + b$, and the fitted values for $a$ and $b$, will be positioned in the upper right corner of the plot.

# Arguments

  - `figure::Makie.Figure`: Makie figure to be drawn over.
  - `error_formating::Symbol=:std_error`: Error format for the annotation. The options are:

      + `:std_error`     -> mean ± standard_error.
      + `:conf_interval` -> mean ± max(upper$_{95\%}$ - mean, mean - lower$_{95\%}$).
  - `color::ColorType=:red`: Color of the line.
  - `linestyle::LineStyleType=nothing`: Style of the line. `nothing` will produce a solid line.
  - `warnings::Bool=true`: If a warning will be raised when there are no points to fit.

# Returns

  - A tuple with the elements for the legend:

      + A `LineElement` to be used as the marker.
      + The label string.
"""
function ppFitLine!(
    figure::Makie.Figure;
    error_formating::Symbol=:std_error,
    color::ColorType=:red,
    linestyle::LineStyleType=nothing,
    warnings::Bool=true,
)::Union{Tuple{Vector{<:LegendElement},Vector{AbstractString}},Nothing}

    # Read the data points in the plot
    points = pointData(figure)

    if isempty(points)
        !warnings || @warn("ppFitLine!: There are no points in the figure")
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

    # Compute the linear fit in scaled (e.g. log10) space
    X = [ones(length(x_points)) x_points]
    linear_model = lm(X, y_points)

    # Read the fitted coeficients
    coeff = coef(linear_model)
    intercept_mean = coeff[1]
    slope_mean = coeff[2]

    # Revert to unscaled values
    x_line = Makie.inverse_transform(x_scaling).(x_points)
    y_line = Makie.inverse_transform(y_scaling).(x_points .* slope_mean .+ intercept_mean)

    # Plot the linear fit
    lines!(figure.current_axis.x, x_line, y_line; color, linestyle)

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

    # Draw the annotation
    text!(
        figure.current_axis.x,
        [absCoor(figure, 0.81, 0.11), absCoor(figure, 0.81, 0.06), absCoor(figure, 0.81, 0.01)];
        text=[L"y = a \, x + b", L"a = %$slope \pm %$δslope", L"b = %$intercept \pm %$δintercept"],
        align=(:left, :bottom),
        fontsize=30,
        color,
    )

    return nothing

end

"""
    ppKennicutt1998!(
        figure::Makie.Figure;
        <keyword arguments>
    )::Union{Tuple{Vector{<:LegendElement},Vector{AbstractString}},Nothing}

Draw a line plot with the fit for the KS relation in Kennicutt (1998).

# Arguments

  - `figure::Makie.Figure`: Makie figure to be drawn over.
  - `x_unit::Unitful.Units=u"Msun * pc^-2"`: Unit for the area density of gas used in `figure`.
  - `y_unit::Unitful.Units=u"Msun * yr^-1 * kpc^-2"`: Unit for the area density of star formation rate used in `figure`.
  - `x_log::Bool=true`: If the x axis is ``\\log_{10}(\\Sigma_\\mathrm{gas})`` (`x_log` = true) or just ``\\Sigma_\\mathrm{gas}`` (`x_log` = false).
  - `y_log::Bool=true`: If the y axis is ``\\log_{10}(\\Sigma_\\mathrm{SFR})`` (`y_log` = true) or just ``\\Sigma_\\mathrm{SFR}``  (`y_log` = false).
  - `color::ColorType=:red`: Color of the line.
  - `linestyle::LineStyleType=nothing`: Style of the line. `nothing` will produce a solid line.

# Returns

  - A tuple with the elements for the legend:

      + A `LineElement` to be used as the marker.
      + The label string.

# References

R. C. Kennicutt (1998). *The Global Schmidt Law in Star-forming Galaxies*. The Astrophysical Journal, **498(2)**, 541-552. [doi:10.1086/305588](https://doi.org/10.1086/305588)
"""
function ppKennicutt1998!(
    figure::Makie.Figure;
    x_unit::Unitful.Units=u"Msun * pc^-2",
    y_unit::Unitful.Units=u"Msun * yr^-1 * kpc^-2",
    x_log::Bool=true,
    y_log::Bool=true,
    color::ColorType=:red,
    linestyle::LineStyleType=nothing,
)::Union{Tuple{Vector{<:LegendElement},Vector{AbstractString}},Nothing}

    # Read the data points in the plot
    points = pointData(figure)

    if isempty(points)
        !warnings || @warn("ppKennicutt1998!: There are no points in the figure")
        return nothing
    end

    # Get the x coordinates of the points
    x_points = sort!(Float64[point[1] for point in points])

    # Compute the area density of gas
    if x_log
        Σg = @. ustrip(KS98_RHO_UNIT, exp10(x_points) * x_unit)
    else
        Σg = @. ustrip(KS98_RHO_UNIT, x_points * x_unit)
    end

    # Compute the values for the y axis
    if y_log
        y_points = @. log10(ustrip(y_unit, KS98_INTERCEPT * Σg^KS98_SLOPE))
    else
        y_points = @. ustrip(y_unit, KS98_INTERCEPT * Σg^KS98_SLOPE)
    end

    # Plot the fit for the KS relation
    lines!(figure.current_axis.x, x_points, y_points; color, linestyle)

    return [LineElement(; color, linestyle)], ["Kennicutt (1998)"]

end

"""
    ppBigiel2008!(
        figure::Makie.Figure,
        quantity::Symbol;
        <keyword arguments>
    )::Tuple{Vector{<:LegendElement},Vector{AbstractString}}

Draw a line plot with the fit for the KS relation in Bigiel et al. (2008).

# Arguments

  - `figure::Makie.Figure`: Makie figure to be drawn over.
  - `quantity::Symbol`: Quantity for the x axis. The options are:

      + `:molecular_area_density` -> Molecular hydrogen area mass density, for a radius of `FILTER_R`.
      + `:neutral_area_density`   -> Neutral hydrogen area mass density, for a radius of `FILTER_R`.
  - `x_unit::Unitful.Units=u"Msun * pc^-2"`: Unit for the area density of gas used in `figure`.
  - `y_unit::Unitful.Units=u"Msun * yr^-1 * kpc^-2"`: Unit for the area density of star formation rate used in `figure`.
  - `x_log::Bool=true`: If the x axis is ``\\log_{10}(\\Sigma_\\mathrm{H})`` (`x_log` = true) or just ``\\Sigma_\\mathrm{H}`` (`x_log` = false).
  - `y_log::Bool=true`: If the y axis is ``\\log_{10}(\\Sigma_\\mathrm{SFR})`` (`y_log` = true) or just ``\\Sigma_\\mathrm{SFR}``  (`y_log` = false).
  - `color::ColorType=:red`: Color of the line.
  - `linestyle::LineStyleType=nothing`: Style of the line. `nothing` will produce a solid line.

# Returns

  - A tuple with the elements for the legend:

      + A `LineElement` to be used as the marker.
      + The label string.

# References

F. Bigiel et al. (2008). *THE STAR FORMATION LAW IN NEARBY GALAXIES ON SUB-KPC SCALES*. The Astrophysical Journal, **136(6)**, 2846. [doi:10.1088/0004-6256/136/6/2846](https://doi.org/10.1088/0004-6256/136/6/2846)
"""
function ppBigiel2008!(
    figure::Makie.Figure,
    quantity::Symbol;
    x_unit::Unitful.Units=u"Msun * pc^-2",
    y_unit::Unitful.Units=u"Msun * yr^-1 * kpc^-2",
    x_log::Bool=true,
    y_log::Bool=true,
    color::ColorType=:red,
    linestyle::LineStyleType=nothing,
)::Tuple{Vector{<:LegendElement},Vector{AbstractString}}

    # Read the data points in the plot
    points = pointData(figure)

    if isempty(points)
        !warnings || @warn("ppBigiel2008!: There are no points in the figure")
        return nothing
    end

    # Get the x coordinates of the points
    x_points = sort!(Float64[point[1] for point in points])

    # Read the file with the fit data
    raw = readdlm(BIGIEL2008_DATA_PATH, skipstart=5, header=true)[1]

    # Parse the fits
    if quantity == :molecular_area_density

        A = parse(Float64, split(raw[end, 2], "+or-")[1])

        b08_intercept = exp10(A) * u"Msun * yr^-1 * kpc^-2"
        b08_slope = parse(Float64, split(raw[end, 3], "+or-")[1])

    elseif quantity == :neutral_area_density

        A = parse(Float64, split(raw[end, 5], "+or-")[1])

        b08_intercept = exp10(A) * u"Msun * yr^-1 * kpc^-2"
        b08_slope = parse(Float64, split(raw[end, 6], "+or-")[1])

    else

        throw(ArgumentError("ppBigiel2008!: `quantity` can only be :molecular_area_density or \
        :neutral_area_density, but I got :$(quantity)"))

    end

    # Compute the area density of gas
    if x_log
        Σg = @. ustrip(u"Msun * pc^-2", exp10(x_points) * 0.1 * x_unit)
    else
        Σg = @. ustrip(u"Msun * pc^-2", x_points * 0.1 * x_unit)
    end

    # Compute the values for the y axis
    if y_log
        y_points = @. log10(ustrip(y_unit, b08_intercept * Σg^b08_slope))
    else
        y_points = @. ustrip(y_unit, b08_intercept * Σg^b08_slope)
    end

    # Plot the fit for the KS relation
    lines!(figure.current_axis.x, x_points, y_points; color, linestyle)

    return [LineElement(; color, linestyle)], ["Bigiel et al. (2008)"]

end

"""
    ppMolla2015!(
        figure::Makie.Figure,
        quantity::Symbol,
        <keyword arguments>
    )::Tuple{Vector{<:LegendElement},Vector{AbstractString}}

Draw a profile for the Milky Way using the data compiled by Mollá et al. (2015).

# Arguments

  - `figure::Makie.Figure`: Makie figure to be drawn over.
  - `quantity::Symbol`: Quantity for the y axis. The options are:

      + `:stellar_area_density`   -> Stellar area mass density.
      + `:molecular_area_density` -> Molecular hydrogen area mass density.
      + `:atomic_area_density`    -> Atomic hydrogen area mass density.
      + `:sfr_area_density`       -> Star formation rate area density.
      + `:X_stellar_abundance`    -> Stellar abundance of element ``\\mathrm{X}``, as ``12 + \\log_{10}(\\mathrm{X \\, / \\, H})``. ``\\mathrm{X}`` can be O (oxygen), N (nitrogen), or C (carbon).
  - `color::ColorType=:red`: Color of the line.
  - `linestyle::LineStyleType=nothing`: Style of the line. `nothing` will produce a solid line.
  - `error_bars::Bool=true`: If the error bars will be plotted.

# Returns

  - A tuple with the elements for the legend:

      + A `MarkerElement` to be used as the marker.
      + The label string.

# References

M. Mollá et al. (2015). *Galactic chemical evolution: stellar yields and the initial mass function*. Monthly Notices of the Royal Astronomical Society **451(4)**, 3693–3708. [doi:10.1093/mnras/stv1102](https://doi.org/10.1093/mnras/stv1102)
"""
function ppMolla2015!(
    figure::Makie.Figure,
    quantity::Symbol;
    color::ColorType=:red,
    linestyle::LineStyleType=nothing,
    error_bars::Bool=true,
)::Tuple{Vector{<:LegendElement},Vector{AbstractString}}

    # Read the file with the compiled data
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

    # Select the quantity for the y axis
    if quantity == :stellar_area_density
        y_data = 10 .^ (raw[!, "logΣ*"] .± raw[!, "logΣ* error"])
    elseif quantity == :molecular_area_density
        y_data = raw[!, "ΣH2"] .± raw[!, "ΣH2 error"]
    elseif quantity == :atomic_area_density
        y_data = raw[!, "ΣHI"] .± raw[!, "ΣHI error"]
    elseif quantity == :sfr_area_density
        y_data = 10 .^ (raw[!, "logΣsfr"] .± raw[!, "logΣsfr error"])
    elseif quantity == :O_stellar_abundance
        y_data = raw[!, "O/H"] .± raw[!, "ΔO/H"]
    elseif quantity == :N_stellar_abundance
        y_data = raw[!, "N/H"] .± raw[!, "ΔN/H"]
    elseif quantity == :C_stellar_abundance
        y_data = raw[!, "C/H"] .± raw[!, "ΔC/H"]
    else
        throw(ArgumentError("ppMolla2015: `x_quantity` can only be  :stellar_area_density, \
        :molecular_area_density, :atomic_area_density, :sfr_area_density, :O_stellar_abundance, \
        :N_stellar_abundance, or :C_stellar_abundance, but I got :$(quantity)"))
    end

    y_values = Measurements.value.(y_data)
    y_uncertainties = Measurements.uncertainty.(y_data)

    # Plot the mean values
    scatterlines!(figure.current_axis.x, x_values, y_values; color, linestyle, marker=:utriangle)

    if error_bars
        # Plot the error bars
        errorbars!(figure.current_axis.x, x_values, y_values, y_uncertainties; color)
    end

    return ([MarkerElement(; color, marker=:utriangle)], ["Mollá et al. (2015)"])

end

"""
    ppFeldmann2020!(
        figure::Makie.Figure,
        x_quantity::Symbol,
        y_quantity::Symbol;
        <keyword arguments>
    )::Tuple{Vector{<:LegendElement},Vector{AbstractString}}

Draw a line, or scatter, plot using the experimental data from the xGASS and xCOLD GASS collaborations, processed by Feldmann (2020).

# Arguments

  - `figure::Makie.Figure`: Makie figure to be drawn over.
  - `x_quantity::Symbol`: Quantity for the x axis. The options are:

      + `:stellar_mass`      -> Stellar mass.
      + `:molecular_mass`    -> Molecular hydrogen (``\\mathrm{H_2}``) mass.
      + `:atomic_mass`       -> Atomic hydrogen (``\\mathrm{HI}``) mass.
      + `:observational_sfr` -> The star formation rate of the last `AGE_RESOLUTION`.
  - `y_quantity::Symbol`: Quantity for the y axis. The options are:

      + `:stellar_mass`      -> Stellar mass.
      + `:molecular_mass`    -> Molecular hydrogen (``\\mathrm{H_2}``) mass.
      + `:atomic_mass`       -> Atomic hydrogen (``\\mathrm{HI}``) mass.
      + `:observational_sfr` -> The star formation rate of the last `AGE_RESOLUTION`.
  - `scatter::Bool=false`: If the data will be presented as a line plot with error bands (default), or alternatively, a scatter plot.

# Returns

  - A tuple with the elements for the legend:

      + A `MarkerElement` to be used as the marker.
      + The label string.

# References

R. Feldmann (2020). *The link between star formation and gas in nearby galaxies*. Communications Physics **3(226)**. [doi:10.1038/s42005-020-00493-0](https://doi.org/10.1038/s42005-020-00493-0)
"""
function ppFeldmann2020!(
    figure::Makie.Figure,
    x_quantity::Symbol,
    y_quantity::Symbol;
    scatter::Bool=false,
)::Tuple{Vector{<:LegendElement},Vector{AbstractString}}

    # Read the CSV file with the raw data
    raw = CSV.read(FELDMANN2020_DATA_PATH, DataFrame, comment="#")

    # Select the quantity for the x axis
    if x_quantity == :stellar_mass

        x_data = exp10.(raw[!, "lgMstar"])

    elseif x_quantity == :molecular_mass

        x_data = raw[!, "MH2"]

    elseif x_quantity == :atomic_mass

        x_data = raw[!, "MHI"]

    elseif x_quantity == :observational_sfr

        x_data = raw[!, "SFR"]

    else

        throw(ArgumentError("ppFeldmann2020!: `x_quantity` can only be :stellar_mass, \
        :molecular_mass, :atomic_mass, or :observational_sfr, but I got :$(x_quantity)"))

    end

    # Select the quantity for the y axix; with its uncertainty
    if y_quantity == :stellar_mass

        # Compute the mean "error" for the stellar mass
        err_low = raw[!, "lgMstar"] .- raw[!, "lgMstar_p16"]
        err_high = raw[!, "lgMstar_p84"] .- raw[!, "lgMstar"]
        err_mean = @. (err_low + err_high) / 2.0

        y_data = exp10.(raw[!, "lgMstar"] .+ err_mean)

    elseif y_quantity == :molecular_mass

        y_data = raw[!, "MH2"] .± raw[!, "e_MH2"]

    elseif y_quantity == :atomic_mass

        y_data = raw[!, "MHI"] .± raw[!, "e_MHI"]

    elseif y_quantity == :observational_sfr

        y_data = raw[!, "SFR"] .± raw[!, "e_SFR"]

    else

        throw(ArgumentError("ppFeldmann2020!: `y_quantity` can only be :stellar_mass, \
        :molecular_mass, :atomic_mass, or :observational_sfr, but I got :$(y_quantity)"))

    end

    if scatter
        y_mean_values = Measurements.value.(y_data)

        # Plot the selected values as a scatter plot
        scatter!(
            figure.current_axis.x,
            x_data,
            y_mean_values;
            color=:red,
            marker=:star4,
            markersize=12,
        )

        return [MarkerElement(; color=:red, marker=:star4)], ["Feldmann (2020)"]
    end

    # Get the scaling of each axis
    x_scaling = xscale(figure)
    y_scaling = yscale(figure)

    # Scale the raw data
    scaled_xs = x_scaling.(x_data)
    scaled_ys = y_scaling.(y_data)

    # Set up the grid
    x_min, x_max = extrema(skipmissing(scaled_xs))
    grid         = CircularGrid(x_max, 30; shift=x_min)
    n_bins       = length(grid.grid)
    bin_width    = (x_max - x_min) / n_bins

    # Allocate memory for the histogram of indices
    histogram = Vector{Int}[Vector{Int}[] for _ in 1:n_bins]

    # Compute the histogram; ignoring missings and values outside the grid range
    @inbounds for (x_idx, x_value) in pairs(scaled_xs)

        if ismissing(x_value)
            continue
        elseif x_value < x_min || x_max < x_value
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

    # Allocate memory
    y_axis = Vector{Union{Measurement{Float64},Missing}}(undef, n_bins)

    @inbounds for (i, hist_idxs) in pairs(histogram)

        # For each bin select the corresponding y measurements
        binned_y_measurements = skipmissing(scaled_ys[hist_idxs])

        if isempty(binned_y_measurements)

            y_axis[i] = missing

        else

            # Compute the weighted mean of all the y measurements in this bin
            weighted_mean = weightedmean(binned_y_measurements)

            # Compute the y mean
            y_mean = Measurements.value(weighted_mean)

            # Compute y uncertainty
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

    x_axis = Makie.inverse_transform(x_scaling).(grid.grid)

    # Delete the empty bins
    deleteat!(x_axis, missing_idxs)

    ################################################################################################
    # Plots
    ################################################################################################

    y_axis_value       = Measurements.value.(y_axis)
    y_axis_uncertainty = Measurements.uncertainty.(y_axis)

    # Plot the mean line
    lines!(
        figure.current_axis.x,
        x_axis,
        Makie.inverse_transform(y_scaling).(y_axis_value);
        color=(:red, 0.4),
    )

    # Construct the 1σ band
    yupper_1σ = Makie.inverse_transform(y_scaling).(y_axis_value .+ y_axis_uncertainty)
    ylower_1σ = Makie.inverse_transform(y_scaling).(y_axis_value .- y_axis_uncertainty)
    # Plot the 1σ band
    band!(figure.current_axis.x, x_axis, ylower_1σ, yupper_1σ; color=(:red, 0.15))

    # Construct the 2σ upper band
    high_band_2σ = Makie.inverse_transform(y_scaling).(y_axis_value .+ 2 * y_axis_uncertainty)
    # Plot the 2σ upper band
    band!(figure.current_axis.x, x_axis, yupper_1σ, high_band_2σ; color=(:orange, 0.15))

    # Construct the 2σ lower band
    low_band_2σ = Makie.inverse_transform(y_scaling).(y_axis_value .- 2 * y_axis_uncertainty)
    # Plot the 2σ lower band
    band!(figure.current_axis.x, x_axis, low_band_2σ, ylower_1σ; color=(:orange, 0.15))

    return (
        [PolyElement(; color=(:red, 0.3)), PolyElement(; color=(:orange, 0.3))],
        [L"\mathrm{Feldmann \,\, (2020)} \,\, 1σ", L"\mathrm{Feldmann \,\, (2020)} \,\, 2σ"],
    )

end
