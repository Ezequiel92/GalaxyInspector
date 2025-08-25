####################################################################################################
# Plotting utilities
####################################################################################################

"""
Extract the limits of the x axis, from a [Makie](https://docs.makie.org/stable/) plot, axis, or figure. In the case of a figure, it will take the limits from the current axis object.
"""
function xlimits!(axis::Makie.Axis)::NTuple{2,Float32}
    reset_limits!(axis)
    return axis.xaxis.attributes.limits[]
end
xlimits!(plot::Makie.FigureAxisPlot)::NTuple{2,Float32} = xlimits!(plot.axis)
xlimits!(fig::Makie.Figure)::NTuple{2,Float32} = xlimits!(fig.current_axis.x)

"""
Extract the limits of the y axis, from a [Makie](https://docs.makie.org/stable/) plot, axis, or figure. In the case of a figure, it will take the limits from the current axis object.
"""
function ylimits!(axis::Makie.Axis)::NTuple{2,Float32}
    reset_limits!(axis)
    return axis.yaxis.attributes.limits[]
end
ylimits!(plot::Makie.FigureAxisPlot)::NTuple{2,Float32} = ylimits!(plot.axis)
ylimits!(fig::Makie.Figure)::NTuple{2,Float32} = ylimits!(fig.current_axis.x)

"""
Extract the scale function of the x axis, from a [Makie](https://docs.makie.org/stable/) plot, axis, or figure. In the case of a figure, it will take the scale from the current axis object.
"""
xscale(axis::Makie.Axis)::Function = axis.xaxis.attributes.scale[]
xscale(plot::Makie.FigureAxisPlot)::Function = xscale(plot.axis)
xscale(fig::Makie.Figure)::Function = xscale(fig.current_axis.x)

"""
Extract the scale function of the y axis, from a [Makie](https://docs.makie.org/stable/) plot, axis, or figure. In the case of a figure, it will take the scale from the current axis object.
"""
yscale(axis::Makie.Axis)::Function = axis.yaxis.attributes.scale[]
yscale(plot::Makie.FigureAxisPlot)::Function = yscale(plot.axis)
yscale(fig::Makie.Figure)::Function = yscale(fig.current_axis.x)

"""
Extract all the data points in a [Makie](https://docs.makie.org/stable/) plot, axis, or figure. In the case of a figure, it will only take the data from the current axis object. It only works for scatter, line, and scatterline plots.
"""
function pointData(axis::Makie.Axis)::Vector{Point{2,Float32}}

    series = copy(axis.scene.plots)

    filter!(x -> isa(x, Union{Scatter,Lines,ScatterLines}), series)

    !isempty(series) || return Point{2,Float32}[]

    return collect(Iterators.flatten(serie[1][] for serie in series))

end
pointData(plot::Makie.FigureAxisPlot)::Vector{Point{2,Float32}} = pointData(plot.axis)
pointData(fig::Makie.Figure)::Vector{Point{2,Float32}} = pointData(fig.current_axis.x)

"""
    absCoor(
        plot::Union{Makie.FigureAxisPlot,Makie.Axis,Makie.Figure},
        r_x::Real,
        r_y::Real,
    )::NTuple{2,Float64}

Compute the absolute x and y coordinates of a plot, from the relative ones.

# Arguments

  - `plot::Union{Makie.FigureAxisPlot,Makie.Axis,Makie.Figure}`: Plot, axis, or figure for which the absolute coordinates will be calculated. In the case of a figure, it will use the limits from the current axis object.
  - `r_x::Real`: Relative x coordinate.
  - `r_y::Real`: Relative y coordinate.

# Returns

  - A tuple with the absolute coordinates, (x, y).

# Examples

```julia-repl
julia> absCoor(lines(rand(100)), 0.5, 0.5)
(50.50000071525574, 0.48792968317866325)
```
"""
function absCoor(
    plot::Union{Makie.FigureAxisPlot,Makie.Axis,Makie.Figure},
    r_x::Real,
    r_y::Real,
)::NTuple{2,Float64}

    # Get the scaling functions
    x_scale = xscale(plot)
    y_scale = yscale(plot)

    # Get the limits of the axes
    x_limits = x_scale.(xlimits!(plot))
    y_limits = y_scale.(ylimits!(plot))

    # Compute the absolute coordinates
    a_x = Makie.inverse_transform(x_scale).(x_limits[1] + abs(r_x) * (x_limits[2] - x_limits[1]))
    a_y = Makie.inverse_transform(y_scale).(y_limits[1] + abs(r_y) * (y_limits[2] - y_limits[1]))

    return sign(r_x) * a_x, sign(r_y) * a_y

end

"""
    cleanPlot!(figure::Makie.Figure)::Nothing

Delete all the legends of a figure and empty all its axes.

# Arguments

  - `figure::Makie.Figure`: Figure to be cleaned.
"""
function cleanPlot!(figure::Makie.Figure)::Nothing

    # Compute the number of elements (axes and legends) in the figure
    n_elements = length(figure.content)
    i = 1

    for _ in 1:n_elements
        i += cleanPlot!(figure.content[i])
    end

    return nothing

end

"""
    cleanPlot!(ax::Makie.Axis)::Bool

Empty an axis.

# Arguments

  - `ax::Makie.Axis`: Axis to be emptied.

# Returns

  - Flag to indicate that an axis has been emptied.
"""
function cleanPlot!(ax::Makie.Axis)::Bool

    empty!(ax)

    return true

end

"""
    cleanPlot!(legend::Union{Makie.Legend,Makie.Colorbar})::Bool

Delete a legend or colorbar.

# Arguments

  - `legend::Union{Makie.Legend,Makie.Colorbar}`: Legend or colorbar to be deleted.

# Returns

  - Flag to indicate that a legend or colorbar has been deleted.
"""
function cleanPlot!(legend::Union{Makie.Legend,Makie.Colorbar})::Bool

    delete!(legend)

    return false

end

"""
Default function to end `cleanPlot!` recursion if an unknown type is encountered.
"""
cleanPlot!(default) = error("cleanPlot!: I cannot clean elements of type $(typeof(default))")

"""
    getUnitLabel(factor::Int, unit::Unitful.Units; <keyword arguments>)::AbstractString

Construct the unit part of an axis label.

# Arguments

  - `factor::Int`: Exponential factor to scale down the units. If different from 0, a term of the form 10^`factor` will be added to the label.
  - `unit::Unitful.Units`: Unit of the axis.
  - `latex::Bool=true`: If the output will be a `LaTeXString` or a plain `String`.

# Returns

  - The `LaTeXString` or `String`: "10^`factor` `unit`". The `factor` term only appears if `factor` != 0, the unit term only appears if `unit` != `Unitful.NoUnits`.
"""
function getUnitLabel(factor::Int, unit::Unitful.Units; latex::Bool=true)::AbstractString

    if latex

        # Replace special characters for their LaTeX counterpart
        str_unit = replace(
            string(unit),
            "M⊙" => raw"M_{\odot}",
            r"\^(-?[0-9]+)" => s"^{\1}",
            " " => raw"\,",
        )

        if iszero(factor)
            if isempty(str_unit)
                out_str = ""
            else
                out_str = L"\mathrm{%$str_unit}"
            end
        else
            if isempty(str_unit)
                out_str = L"10^{%$factor}"
            else
                out_str = L"10^{%$factor} \, \mathrm{%$str_unit}"
            end
        end

    else

        str_unit = string(unit)

        if iszero(factor)
            if isempty(str_unit)
                out_str = ""
            else
                out_str = str_unit
            end
        else
            if isempty(str_unit)
                out_str = "10^$factor"
            else
                out_str = "10^$factor $str_unit"
            end
        end

    end

    return out_str

end

"""
    getLabel(
        label::AbstractString,
        factor::Int,
        unit::Unitful.Units;
        <keyword arguments>
    )::AbstractString

Construct an axis label.

# Arguments

  - `label::AbstractString`: Variable name.
  - `factor::Int`: Exponential factor to scale down the units. If different from 0, a term of the form 10^`factor` will be added to the label.
  - `unit::Unitful.Units`: Unit of the axis.
  - `latex::Bool=true`: If the output will be a `LaTeXString` or a plain `String`.

# Returns

  - The `LaTeXString` or `String`: "`label` [10^`factor` `unit`]". If `label` is "", an empty string is returned. The `factor` term only appears if `factor` != 0, the unit term only appears if `unit` != `Unitful.NoUnits`, and the brackets only appears if there are a factor and/or a unit term.
"""
function getLabel(
    label::AbstractString,
    factor::Int,
    unit::Unitful.Units;
    latex::Bool=true,
)::AbstractString

    !isempty(label) || return ""

    unit_label = getUnitLabel(factor, unit; latex)

    if isempty(unit_label)
        return latex ? L"%$label" : label
    end

    return latex ? L"%$label [%$unit_label]" : "$label [$unit_label]"

end

"""
    barPlotLabelFormater(x::Number)::LaTeXString

Format a number to be a barplot label.

For values between 0 and 0.01 the label will be "< 0.01", otherwise it will be the value itself with 2 digits.

# Arguments

  - `x::Number`: Value to be formatted.

# Returns

  - The bar label.
"""
function barPlotLabelFormater(x::Number)::LaTeXString

    if 0 < x < 0.01
        return L"< \, 0.01"
    end

    return latexstring(round(x; digits=2))

end

"""
    barPlotLabelFormater(x::LaTeXString)::LaTeXString

Format a number to be a barplot label.

Method for compatibility with the barplot! function of [Makie](https://docs.makie.org/stable/).

# Arguments

  - `x::Number`: Value to be formatted.

# Returns

  - The bar label.
"""
barPlotLabelFormater(x::LaTeXString)::LaTeXString = x

"""
    reduceTicks(hr_ticks::Vector{<:Number}, factor::Int)::Vector{<:Number}

Reduce the length of a given list of axis ticks, while keeping the total length of the axis the same.

# Arguments

  - `hr_ticks::Vector{<:Number}`: Original "high resolution" list of ticks. It has to be regularly spaced
  - `factor::Int`: Factor by which the number of values will be reduced. It has to divide the size of `hr_ticks` exactly.

# Returns

  - The new shorter list of ticks.
"""
function reduceTicks(hr_ticks::Vector{<:Number}, factor::Int)::Vector{<:Number}

    !isone(factor) || return hr_ticks

    l = length(hr_ticks)

    (
        l % factor == 0 ||
        throw(ArgumentError("reduceTicks: `factor` must divide the size of `hr_ticks` \
        exactly, but I got length of `hr_ticks` / `factor` = $(l / factor)"))
    )

    (
        factor >= 1 ||
        throw(ArgumentError("reduceTicks: `factor` must be >= 1, but I got `factor` = $(factor)"))
    )

    # Compute the size of the new vector
    new_size = l ÷ factor

    # Allocate memory
    lr_ticks = similar(hr_ticks, new_size)

    if iseven(factor)

        shift = factor ÷ 2

        for i in eachindex(lr_ticks)

            idx = (i - 1) * factor + shift

            lr_ticks[i] = (hr_ticks[idx] + hr_ticks[idx + 1]) / 2.0

        end

    else

        shift = ceil(Int, factor / 2)

        for i in eachindex(lr_ticks)

            idx = (i - 1) * factor + shift

            lr_ticks[i] = hr_ticks[idx]

        end

    end

    return lr_ticks

end

"""
    formatError(q_mean::Number, q_error::Number)::NTuple{2,<:Number}

Nicely format a magnitude with uncertainty.

It follows the traditional rules for error presentation: the error has only one significant digit, unless such digit is a one, in which case two significant digits are used. The mean will have as many digits as to match the last significant position of the error. An error equal to 0 will leave the mean unchanged.

# Arguments

  - `q_mean::Number`: Mean value.
  - `q_error::Number`: Error value. It must be positive.

# Returns

  - A tuple with the formatted mean and error values.

# Examples

```julia-repl
julia> formatError(69.42069, 0.038796)
(69.42, 0.04)

julia> formatError(69.42069, 0.018796)
(69.421, 0.019)

julia> formatError(15.42, 0.00004)
(15.42, 4.0e-5)

julia> formatError(69.42069, 0.0)
(69.42069, 0.0)

julia> formatError(69.42069, 93.4)
(70.0, 90.0)

julia> formatError(69.42069, 123.4)
(70.0, 120.0)

julia> formatError(15.42069, 16.4)
(15.0, 16.0)
```
"""
function formatError(q_mean::Number, q_error::Number)::NTuple{2,<:Number}

    # Positive error check
    (
        q_error >= zero(q_error) ||
        throw(ArgumentError("formatError: `q_error` must be positive, but I got \
        `q_error` = $(q_error)"))
    )

    # Use the values without units
    mean = ustrip(q_mean)
    error = ustrip(q_error)

    if iszero(error)

        round_mean = mean
        round_error = error

    else

        sigdigit_pos = abs(log10(error))

        if error < 1.0
            first_digit = trunc(error * 10.0^(floor(sigdigit_pos) + 1.0))
            extra = first_digit == 1.0 ? 1 : 0
            digits = ceil(Int, sigdigit_pos) + extra
            round_mean = round(mean; digits)
            round_error = round(error, sigdigits=1 + extra)
        else
            first_digit = trunc(error * 10.0^(-floor(sigdigit_pos)))
            extra = first_digit == 1.0 ? 2 : 1
            sigdigits = ceil(Int, log10(abs(mean))) - ceil(Int, sigdigit_pos) + extra
            round_mean = round(mean; sigdigits)
            round_error = round(error, sigdigits=extra)
        end

    end

    return round_mean * unit(q_mean), round_error * unit(q_error)

end

















#TODO
function derivedQuantitiesPlotParams(magnitude::Symbol, component::Symbol)::PlotParams

    ######################
    # Particle components
    ######################
    if component == :stellar

        c_label = "\\star"

    elseif component == :dark_matter

        c_label = "\\text{DM}"

    elseif component == :black_hole

        c_label = "\\text{BH}"

    #################
    # Gas components
    #################
    elseif component == :gas

        c_label = "\\text{gas}"

    elseif component == :hydrogen

        c_label = "\\text{H}"

    elseif component == :helium

        c_label = "\\text{He}"

    elseif component == :metals

        c_label = "\\text{Z}"

    elseif component == :ionized

        c_label = "\\text{HII}"

    elseif component == :neutral

        c_label = "\\mathrm{HI + H2}"

    elseif component == :br_atomic

        c_label = "\\text{BR HI}"

    elseif component == :br_molecular

        c_label = "\\text{BR H2}"

    ###########################################
    # Components from our star formation model
    ###########################################
    elseif component == :ode_ionized

        c_label = "\\text{i}"

    elseif component == :ode_atomic

        c_label = "\\text{a}"

    elseif component == :ode_molecular

        c_label = "\\text{m}"

    elseif component == :ode_stellar

        c_label = "\\text{s}"

    elseif component == :ode_metals

        c_label = "\\text{Z}"

    elseif component == :ode_dust

        c_label = "\\text{d}"

    elseif component == :ode_neutral

        c_label = "\\mathrm{a + m}"

    elseif component == :generic

        c_label = ""

    else

        throw(ArgumentError("derivedQuantitiesPlotParams: I don't recognize the component \
        :$(component)"))

    end

    #######################
    # Mass-like magnitudes
    #######################
    if magnitude == :mass

        request = Dict{Symbol,Vector{String}}()

        if isempty(c_label)
            var_name = L"M"
        else
            var_name = L"M_%$(c_label)"
        end

        exp_factor = 10

        unit = u"Msun"

    elseif magnitude == :mass_density

        request = Dict{Symbol,Vector{String}}()

        if isempty(c_label)
            var_name = L"\rho"
        else
            var_name = L"\rho_%$(c_label)"
        end

        exp_factor = 0

        unit = u"Msun * kpc^-3"

    elseif magnitude == :number_density

        request = Dict{Symbol,Vector{String}}()

        if isempty(c_label)
            var_name = L"n"
        else
            var_name = L"n_{\,%$(c_label)}"
        end

        exp_factor = 0

        unit = u"cm^-3"

    elseif magnitude == :area_density

        request = Dict{Symbol,Vector{String}}()

        if isempty(c_label)
            var_name = L"\Sigma"
        else
            var_name = L"\Sigma_%$(c_label)"
        end

        exp_factor = 0

        unit = u"Msun * pc^-2"

    elseif magnitude == :number

        request = Dict{Symbol,Vector{String}}()

        if isempty(c_label)
            var_name = L"N"
        else
            var_name = L"N_%$(c_label)"
        end

        exp_factor = 0

        unit = Unitful.NoUnits

    elseif magnitude == :fraction

        request = Dict{Symbol,Vector{String}}()

        if isempty(c_label)
            var_name = L"f"
        else
            var_name = L"f_{\,%$(c_label)}"
        end

        exp_factor = 0

        unit = Unitful.NoUnits

    #######################
    # Cinematic magnitudes
    #######################
    elseif magnitude == :specific_z_angular_momentum

        # See computeSpecificAngularMomentum() in ./src/analysis/compute_quantities/velocities.jl
        blocks = ["VEL ", "POS "]
        types  = [:gas, :dark_matter, :stellar, :black_hole]

        if c_label ∈ types
            request = Dict(c_label => blocks)
        elseif component == :generic
            request = Dict(type => blocks for type in types)
        else
            request = Dict(:gas => blocks)
        end

        if isempty(c_label)
            var_name = L"j_z"
        else
            var_name = L"j_{z, \,\,%$(c_label)}"
        end

        exp_factor = 0

        unit = u"kpc^2 * s^-1"

    elseif magnitude == :z_angular_momentum

        # See computeAngularMomentum() in ./src/analysis/compute_quantities/velocities.jl
        blocks = ["VEL ", "POS ", "MASS"]
        types  = [:gas, :dark_matter, :stellar, :black_hole]

        if c_label ∈ types
            request = Dict(c_label => blocks)
        elseif component == :generic
            request = Dict(type => blocks for type in types)
        else
            request = Dict(:gas => blocks)
        end

        if isempty(c_label)
            var_name = L"l_z"
        else
            var_name = L"l_{z, \,\,%$(c_label)}"
        end

        exp_factor = 0

        unit = u"Msun * kpc^2 * s^-1"

    elseif magnitude == :eff

        request = Dict{Symbol,Vector{String}}()

        if isempty(c_label)
            var_name = L"\epsilon_\text{ff}"
        else
            var_name = L"\epsilon_{\text{ff}, \,\,%$(c_label)}"
        end

        exp_factor = 0

        unit = Unitful.NoUnits

    elseif magnitude == :spin_parameter

        # See computeSpinParameter() in ./src/analysis/compute_quantities/velocities.jl
        blocks = ["VEL ", "POS ", "MASS"]
        types  = [:gas, :dark_matter, :stellar, :black_hole]

        if c_label ∈ types
            request = Dict(c_label => blocks)
        elseif component == :generic
            request = Dict(type => blocks for type in types)
        else
            request = Dict(:gas => blocks)
        end

        if isempty(c_label)
            var_name = L"\lambda"
        else
            var_name = L"\lambda_%$(c_label)"
        end

        exp_factor = 0

        unit = Unitful.NoUnits

    elseif magnitude == :circularity

        # See computeCircularity() in ./src/analysis/compute_quantities/velocities.jl
        request = Dict(type => ["VEL", "POS ", "MASS"] for type in keys(PARTICLE_INDEX))

        if isempty(c_label)
            var_name = L"\epsilon"
        else
            var_name = L"\epsilon_{\,%$(c_label)}"
        end

        exp_factor = 0

        unit = Unitful.NoUnits

    elseif magnitude == :circular_velocity

        # See computeVcirc() in ./src/analysis/compute_quantities/velocities.jl
        request = Dict(type => ["MASS", "POS "] for type in keys(PARTICLE_INDEX))

        if isempty(c_label)
            var_name = L"v_\text{circ}"
        else
            var_name = L"v_{\text{circ}, \,\,%$(c_label)}"
        end

        exp_factor = 0

        unit = u"km * s^-1"

    elseif magnitude == :radial_velocity

        # See computeVpolar() in ./src/analysis/compute_quantities/velocities.jl
        blocks = ["VEL ", "POS "]
        types  = [:gas, :dark_matter, :stellar, :black_hole]

        if c_label ∈ types
            request = Dict(c_label => blocks)
        elseif component == :generic
            request = Dict(type => blocks for type in types)
        else
            request = Dict(:gas => blocks)
        end

        if isempty(c_label)
            var_name = L"v_r"
        else
            var_name = L"v_{r, \,\,%$(c_label)}"
        end

        exp_factor = 0

        unit = u"km * s^-1"

    elseif magnitude == :tangential_velocity

        # See computeVpolar() in ./src/analysis/compute_quantities/velocities.jl
        blocks = ["VEL ", "POS "]
        types  = [:gas, :dark_matter, :stellar, :black_hole]

        if c_label ∈ types
            request = Dict(c_label => blocks)
        elseif component == :generic
            request = Dict(type => blocks for type in types)
        else
            request = Dict(:gas => blocks)
        end

        if isempty(c_label)
            var_name = L"v_\theta"
        else
            var_name = L"v_{\theta, \,\,%$(c_label)}"
        end

        exp_factor = 0

        unit = u"km * s^-1"

    elseif magnitude == :zstar_velocity

        # See computeVpolar() in ./src/analysis/compute_quantities/velocities.jl
        blocks = ["VEL ", "POS "]
        types  = [:gas, :dark_matter, :stellar, :black_hole]

        if c_label ∈ types
            request = Dict(c_label => blocks)
        elseif component == :generic
            request = Dict(type => blocks for type in types)
        else
            request = Dict(:gas => blocks)
        end

        if isempty(c_label)
            var_name = L"v_z \,\, \mathrm{sign}(z)"
        else
            var_name = L"v_{z, \,\,%$(c_label)} \,\, \mathrm{sign}(z)"
        end

        exp_factor = 0

        unit = u"km * s^-1"

    ########
    # Other
    ########
    elseif magnitude == :depletion_time

        request = Dict{Symbol,Vector{String}}()

        if isempty(c_label)
            var_name = L"\tau_\text{dep}"
        else
            var_name = L"\tau_{\text{dep}, \,\,%$(c_label)}"
        end

        exp_factor = 0

        unit = u"Gyr"

    elseif magnitude == :xy_distance

        request = Dict{Symbol,Vector{String}}()

        if isempty(c_label)
            var_name = L"d_{xy}"
        else
            var_name = L"d_{xy, \,\,%$(c_label)}"
        end

        exp_factor = 0

        unit = u"kpc"

    elseif magnitude == :radial_distance

        request = Dict{Symbol,Vector{String}}()

        if isempty(c_label)
            var_name = L"r"
        else
            var_name = L"r_{\,%$(c_label)}"
        end

        exp_factor = 0

        unit = u"kpc"

    else

        throw(ArgumentError("derivedQuantitiesPlotParams: I don't recognize the magnitude \
        :$(magnitude)"))

    end

    return PlotParams(; request, var_name, exp_factor, unit)

end






#TODO
function plotParams2(quantity::Symbol)::PlotParams


    if quantity == :ode_metallicity

        plot_params = PlotParams(;
            request  = Dict(:gas => ["MASS", "POS ", "FRAC", "NH  ", "NHP ", "RHO ", "GZ  "]),
            var_name = L"Z_\mathrm{gas}^\star \, [\mathrm{Z_\odot}]",
        )

    elseif quantity == :molecular_neutral_fraction

        plot_params = PlotParams(;
            request  = Dict(
                :gas => ["RHO ", "MASS", "POS ", "FRAC", "NH  ", "NHP ", "PRES", "RHO ", "GZ  "],
            ),
            var_name = L"f_\mathrm{H_2} \, / f_\mathrm{n}",
        )

    elseif quantity == :ionized_neutral_fraction

        plot_params = PlotParams(;
            request  = Dict(
                :gas => ["RHO ", "MASS", "POS ", "FRAC", "NH  ", "NHP ", "PRES", "RHO ", "GZ  "],
            ),
            var_name = L"f_\mathrm{HII} \, / f_\mathrm{n}",
        )

    elseif quantity == :sfr_area_density

        plot_params = PlotParams(;
            request  = Dict(:stellar => ["MASS", "POS ", "GAGE"]),
            var_name = L"\Sigma_\mathrm{SFR}",
            unit     = u"Msun * yr^-1 * kpc^-2",
        )

    elseif quantity == :gas_metallicity

        plot_params = PlotParams(;
            request  = Dict(:gas => ["MASS", "POS ", "GMET", "GZ  "]),
            var_name = L"Z_\mathrm{gas} \, [\mathrm{Z_\odot}]",
        )

    elseif quantity == :stellar_metallicity

        plot_params = PlotParams(;
            request  = Dict(:stellar => ["MASS", "POS ", "GME2", "GZ2 "]),
            var_name = L"Z_\star \, [\mathrm{Z_\odot}]",
        )

    elseif quantity == :gas_sfr

        plot_params = PlotParams(;
            request  = Dict(:gas => ["SFR "]),
            var_name = L"\mathrm{SFR_{gas}}",
            unit     = u"Msun * yr^-1",
        )

    elseif quantity == :mass_accretion

        plot_params = PlotParams(;
            request  = Dict(
                :gas        => ["ID  ", "MASS"],
                :stellar    => ["ID  ", "MASS"],
                :black_hole => ["ID  ", "MASS"],
                :group      => ["G_R_Crit200", "G_M_Crit200"],
                :tracer     => ["PAID", "TRID"],
            ),
            var_name = "Mass accretion",
            unit     = u"Msun * yr^-1",
        )

    elseif quantity == :stellar_age

        plot_params = PlotParams(;
            request  = Dict(:stellar => ["GAGE"]),
            var_name = L"\mathrm{Stellar \,\, age}",
            unit     = u"Gyr",
        )

    elseif quantity == :sfr

        plot_params = PlotParams(;
            request  = Dict(:stellar => ["MASS", "POS ", "GAGE"]),
            var_name = L"\mathrm{SFR}",
            unit     = u"Msun * yr^-1",
        )

    elseif quantity == :ssfr

        plot_params = PlotParams(;
            request  = Dict(:stellar => ["MASS", "POS ", "GAGE"]),
            var_name = L"\mathrm{sSFR}",
            unit     = u"yr^-1",
        )

    elseif quantity == :observational_sfr

        plot_params = PlotParams(;
            request  = Dict(:stellar => ["MASS", "POS ", "GAGE"]),
            var_name = L"\mathrm{SFR}",
            unit     = u"Msun * yr^-1",
        )

    elseif quantity == :observational_ssfr

        plot_params = PlotParams(;
            request  = Dict(:stellar => ["MASS", "POS ", "GAGE"]),
            var_name = L"\mathrm{sSFR}",
            unit     = u"yr^-1",
        )

    elseif quantity == :temperature

        plot_params = PlotParams(;
            request    = Dict(:gas => ["MASS", "POS ", "TEMP"]),
            axis_label = L"\log_{10} \, T \, [\mathrm{K}]",
        )

    elseif quantity == :pressure

        plot_params = PlotParams(; request = Dict(:gas => ["PRES"]), var_name = L"P", unit = u"Pa")

    elseif quantity == :scale_factor

        plot_params = PlotParams(; var_name=L"a")

    elseif quantity == :redshift

        plot_params = PlotParams(; var_name=L"z")

    elseif quantity == :physical_time

        plot_params = PlotParams(; var_name=L"t", unit=u"Gyr")

    elseif quantity == :lookback_time

        plot_params = PlotParams(; var_name = L"\mathrm{Lookback \,\, time}", unit = u"Gyr")

    elseif quantity == :time_step

        plot_params = PlotParams(; var_name=L"\mathrm{Number \,\, of \,\, time \,\, steps}")

    elseif quantity == :clock_time_s

        plot_params = PlotParams(; var_name = L"\mathrm{Wallclock \,\, time}", unit = u"s")

    elseif quantity == :clock_time_percent

        plot_params = PlotParams(; axis_label=L"\mathrm{Wallclock \,\, time \, [\%]}")

    elseif quantity == :tot_clock_time_s

        plot_params =
            PlotParams(; var_name = L"\mathrm{Cumulative \,\, wallclock \,\, time}", unit = u"s")

    elseif quantity == :tot_clock_time_percent

        plot_params = PlotParams(; axis_label=L"\mathrm{Cumulative \,\, wallclock \,\, time \, [\%]}")

    elseif quantity ∈ GAS_ABUNDANCE

        element = get(GAS_ABUNDANCE_SPLITS, quantity)

        plot_params = PlotParams(;
            request    = Dict(:gas => ["MASS", "POS ", "GMET"]),
            axis_label = L"12 + \log_{10}(\mathrm{%$element} \, / \, \mathrm{H})",
        )

    elseif quantity ∈ STELLAR_ABUNDANCE

        element = get(STELLAR_ABUNDANCE_SPLITS, quantity)

        plot_params = PlotParams(;
            request    = Dict(:stellar => ["MASS", "POS ", "GME2"]),
            axis_label = L"12 + \log_{10}(\mathrm{%$element} \, / \, \mathrm{H})",
        )

     elseif quantity ∈ DERIVED_QUANTITIES

        magnitude, component = splitQuantity(quantity)

        plot_params = derivedQuantitiesPlotParams(magnitude, component)

    elseif quantity ∈ MAGNITUDES

        plot_params = derivedQuantitiesPlotParams(quantity, :generic)

    else

        throw(ArgumentError("plotParams: I don't recognize the quantity :$(quantity)"))

    end

    return plot_params

end





















"""
    plotParams(quantity::Symbol)::PlotParams

Return the plotting parameters for a given `quantity`.

# Arguments

  - `quantity::Symbol`: The options are:

      + `:stellar_mass`                -> Stellar mass.
      + `:gas_mass`                    -> Gas mass.
      + `:hydrogen_mass`               -> Hydrogen mass.
      + `:dm_mass`                     -> Dark matter mass.
      + `:bh_mass`                     -> Black hole mass.
      + `:molecular_mass`              -> Molecular hydrogen (``\\mathrm{H_2}``) mass.
      + `:br_molecular_mass`           -> Molecular hydrogen (``\\mathrm{H_2}``) mass, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_mass`                 -> Atomic hydrogen (``\\mathrm{HI}``) mass.
      + `:ionized_mass`                -> Ionized hydrogen (``\\mathrm{HII}``) mass.
      + `:neutral_mass`                -> Neutral hydrogen (``\\mathrm{HI + H_2}``) mass.
      + `:stellar_gas_mass`            -> Stellar gas mass (according to our SF model).
      + `:ode_metal_mass`              -> Metal mass (according to our SF model).
      + `:ode_metallicity`             -> Metallicity (according to our SF model).
      + `:dust_mass`                   -> Dust mass.
      + `:generic_mass`                -> Parameters for plots with several different masses.
      + `:stellar_number`              -> Number of stellar particles.
      + `:gas_number`                  -> Number of gas cells.
      + `:dm_number`                   -> Number of dark matter particles.
      + `:bh_number`                   -> Number of black hole particles.
      + `:molecular_fraction`          -> Gas mass fraction of molecular hydrogen.
      + `:br_molecular_fraction`       -> Gas mass fraction of molecular hydrogen, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_fraction`             -> Gas mass fraction of atomic hydrogen.
      + `:ionized_fraction`            -> Gas mass fraction of ionized hydrogen.
      + `:neutral_fraction`            -> Gas mass fraction of neutral hydrogen.
      + `:molecular_neutral_fraction`  -> Fraction of molecular hydrogen in the neutral gas.
      + `:ionized_neutral_fraction`    -> Fraction of ionized gas to neutral gas.
      + `:stellar_gas_fraction`        -> Stellar gas fraction (according to our SF model).
      + `:metal_gas_fraction`          -> Metallicity (according to our SF model).
      + `:dust_fraction`               -> Dust mass fraction.
      + `:generic_fraction`            -> Parameters for plots with several different fraction.
      + `:gas_mass_density`            -> Gas mass density.
      + `:hydrogen_mass_density`       -> Hydrogen mass density.
      + `:gas_number_density`          -> Gas number density.
      + `:molecular_number_density`    -> Molecular hydrogen number density.
      + `:br_molecular_number_density` -> Molecular hydrogen number density, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_number_density`       -> Atomic hydrogen number density.
      + `:ionized_number_density`      -> Ionized hydrogen number density.
      + `:neutral_number_density`      -> Neutral hydrogen number density.
      + `:stellar_area_density`        -> Stellar area mass density, for a radius of `DISK_R`.
      + `:gas_area_density`            -> Gas mass surface density, for a radius of `DISK_R`.
      + `:molecular_area_density`      -> Molecular mass surface density, for a radius of `DISK_R`.
      + `:br_molecular_area_density`   -> Molecular mass surface density, for a radius of `DISK_R`, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_area_density`         -> Atomic hydrogen area mass density, for a radius of `DISK_R`.
      + `:ionized_area_density`        -> Ionized hydrogen area mass density, for a radius of `DISK_R`.
      + `:neutral_area_density`        -> Neutral mass surface density, for a radius of `DISK_R`.
      + `:sfr_area_density`            -> Star formation rate area density, for the last `AGE_RESOLUTION` and a radius of `DISK_R`.
      + `:generic_area_density`        -> Parameters for plots with several different area densities.
      + `:gas_td`                      -> Total gas depletion time.
      + `:molecular_td`                -> Molecular hydrogen (``\\mathrm{H_2}``) depletion time.
      + `:br_molecular_td`             -> Molecular hydrogen (``\\mathrm{H_2}``) depletion time, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_td`                   -> Atomic hydrogen (``\\mathrm{HI}``) depletion time.
      + `:ionized_td`                  -> Ionized hydrogen (``\\mathrm{HII}``) depletion time.
      + `:neutral_td`                  -> Neutral hydrogen (``\\mathrm{HI + H_2}``) depletion time.
      + `:gas_metallicity`             -> Mass fraction of all elements above He in the gas (solar units).
      + `:stellar_metallicity`         -> Mass fraction of all elements above He in the stars (solar units).
      + `:X_gas_abundance`             -> Gas abundance of element ``\\mathrm{X}``, as ``12 + \\log_{10}(\\mathrm{X \\, / \\, H})``. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
      + `:X_stellar_abundance`         -> Stellar abundance of element ``\\mathrm{X}``, as ``12 + \\log_{10}(\\mathrm{X \\, / \\, H})``. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
      + `:stellar_radial_distance`     -> Distance of every stellar particle to the origin.
      + `:gas_radial_distance`         -> Distance of every gas cell to the origin.
      + `:dm_radial_distance`          -> Distance of every dark matter particle to the origin.
      + `:stellar_xy_distance`         -> Projected distance of every stellar particle to the origin.
      + `:gas_xy_distance`             -> Projected distance of every gas cell to the origin.
      + `:dm_xy_distance`              -> Projected distance of every dark matter particle to the origin.
      + `:gas_sfr`                     -> SFR associated to each gas particle/cell within the code.
      + `:mass_accretion`              -> Gas accretion rate. Positive values mean gas infall into the virial radius ``R_{200}``, and negative values mean outflow.
      + `:stellar_specific_am`         -> Norm of the stellar specific angular momentum.
      + `:gas_specific_am`             -> Norm of the gas specific angular momentum.
      + `:dm_specific_am`              -> Norm of the dark matter specific angular momentum.
      + `:stellar_circularity`         -> Stellar circularity.
      + `:stellar_vcirc`               -> Stellar circular velocity.
      + `:stellar_vradial`             -> Stellar radial speed.
      + `:stellar_vtangential`         -> Stellar tangential speed.
      + `:stellar_vzstar`              -> Stellar speed in the z direction, computed as ``v_z \\, \\mathrm{sign}(z)``.
      + `:stellar_age`                 -> Stellar age.
      + `:sfr`                         -> Star formation rate.
      + `:ssfr`                        -> The specific star formation rate.
      + `:observational_sfr`           -> Star formation rate of the last `AGE_RESOLUTION`.
      + `:observational_ssfr`          -> The specific star formation rate of the last `AGE_RESOLUTION`.
      + `:stellar_eff`                 -> Star formation efficiency per free-fall time for the gas that has turn into stars.
      + `:gas_eff`                     -> Star formation efficiency per free-fall time for the gas.
      + `:molecular_eff`               -> Star formation efficiency per free-fall time for the molecular hydrogen (``\\mathrm{H_2}``) gas.
      + `:br_molecular_eff`            -> Star formation efficiency per free-fall time for the molecular hydrogen (``\\mathrm{H_2}``) gas, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_eff`                  -> Star formation efficiency per free-fall time for the atomic hydrogen (``\\mathrm{HI}``) gas.
      + `:ionized_eff`                 -> Star formation efficiency per free-fall time for the ionized hydrogen (``\\mathrm{HII}``) gas.
      + `:neutral_eff`                 -> Star formation efficiency per free-fall time for the neutral hydrogen (``\\mathrm{HI + H_2}``) gas.
      + `:temperature`                 -> Gas temperature, as ``\\log_{10}(T \\, / \\, \\mathrm{K})``.
      + `:pressure`                    -> Gas pressure.
      + `:scale_factor`                -> Scale factor.
      + `:redshift`                    -> Redshift.
      + `:physical_time`               -> Physical time since the Big Bang.
      + `:lookback_time`               -> Physical time left to reach the last snapshot.
      + `:ode_gas_it`                  -> Integration time.
      + `:ode_gas_tau_s`               -> Star formation time scale, ``\\tau_\\mathrm{S}``.
      + `:ode_gas_eta_d`               -> Photodissociation efficiency, ``\\eta_\\mathrm{diss}``.
      + `:ode_gas_eta_i`               -> Photoionization efficiency, ``\\eta_\\mathrm{ion}``.
      + `:ode_gas_r`                   -> Mass recycling parameter, ``R``.
      + `:ode_gas_cold_mf`             -> Cold gas mass fraction.
      + `:ode_stellar_it`              -> Integration time, for the gas that form the stars.
      + `:ode_stellar_tau_s`           -> Star formation time scale, ``\\tau_\\mathrm{S}``, for the gas that form the stars.
      + `:ode_stellar_eta_d`           -> Photodissociation efficiency, ``\\eta_\\mathrm{diss}``, for the gas that form the stars.
      + `:ode_stellar_eta_i`           -> Photoionization efficiency, ``\\eta_\\mathrm{ion}``, for the gas that form the stars.
      + `:ode_stellar_r`               -> Mass recycling parameter, ``R``, for the gas that form the stars.
      + `:ode_stellar_cold_mf`         -> Cold gas mass fraction, for the gas that form the stars.
      + `:ode_stellar_gas_rho`         -> Gas mass density, for the gas that form the stars.
      + `:ode_stellar_gas_Z`           -> Gas metallicity, for the gas that form the stars (solar units).
      + `:ode_stellar_gas_mass`        -> Cell mass, for the gas that form the stars.
      + `:ode_stellar_gas_sfr`         -> SFR associated to the gas particles/cells within the code, for the gas that form the stars.
      + `:ode_stellar_gas_P`           -> Gas pressure, for the gas that form the stars.

# Returns

  - A [`PlotParams`](@ref) object, with entries:

      + `request::Dict{Symbol,Vector{String}}` -> Data request for [`readSnapshot`](@ref).
      + `var_name::AbstractString`             -> Name of the quantity for the plot axis.
      + `exp_factor::Int`                      -> Numerical exponent to scale down the axis.
      + `unit::Unitful.Units`                  -> Target unit for the axis.
      + `axis_label::AbstractString`           -> Label for the axis.

# References

L. Blitz et al. (2006). *The Role of Pressure in GMC Formation II: The H2-Pressure Relation*. The Astrophysical Journal, **650(2)**, 933. [doi:10.1086/505417](https://doi.org/10.1086/505417)
"""
function plotParams(quantity::Symbol)::PlotParams

    if quantity == :stellar_mass

        plot_params = PlotParams(;
            request    = Dict(:stellar => ["MASS", "POS "]),
            var_name   = L"M_\star",
            exp_factor = 10,
            unit       = u"Msun",
        )

    elseif quantity == :gas_mass

        plot_params = PlotParams(;
            request    = Dict(:gas => ["MASS", "POS ", "RHO "]),
            var_name   = L"M_\mathrm{gas}",
            exp_factor = 10,
            unit       = u"Msun",
        )

    elseif quantity == :hydrogen_mass

        plot_params = PlotParams(;
            request    = Dict(:gas => ["MASS", "POS ", "RHO "]),
            var_name   = L"M_\mathrm{H}",
            exp_factor = 10,
            unit       = u"Msun",
        )

    elseif quantity == :dm_mass

        plot_params = PlotParams(;
            request    = Dict(:dark_matter => ["MASS", "POS "]),
            var_name   = L"M_\mathrm{DM}",
            exp_factor = 10,
            unit       = u"Msun",
        )

    elseif quantity == :bh_mass

        plot_params = PlotParams(;
            request  = Dict(:black_hole => ["MASS", "POS "]),
            var_name = L"M_\mathrm{BH}",
            unit     = u"Msun",
        )

    elseif quantity == :molecular_mass

        plot_params = PlotParams(;
            request    = Dict(:gas => ["MASS", "POS ", "FRAC", "RHO ", "PRES", "NH  ", "NHP "]),
            var_name   = L"M_\mathrm{H_2}",
            exp_factor = 10,
            unit       = u"Msun",
        )

    elseif quantity == :br_molecular_mass

        plot_params = PlotParams(;
            request    = Dict(:gas => ["MASS", "POS ", "PRES"]),
            var_name   = L"M_\mathrm{H_2^{BR}}",
            exp_factor = 10,
            unit       = u"Msun",
        )

    elseif quantity == :atomic_mass

        plot_params = PlotParams(;
            request    = Dict(
                :gas => ["MASS", "POS ", "FRAC", "NH  ", "NHP ", "RHO ", "GZ  ", "PRES"],
            ),
            var_name   = L"M_\mathrm{HI}",
            exp_factor = 10,
            unit       = u"Msun",
        )

    elseif quantity == :ionized_mass

        plot_params = PlotParams(;
            request    = Dict(:gas => ["MASS", "POS ", "FRAC", "NH  ", "NHP ", "RHO ", "GZ  "]),
            var_name   = L"M_\mathrm{HII}",
            exp_factor = 10,
            unit       = u"Msun",
        )

    elseif quantity == :neutral_mass

        plot_params = PlotParams(;
            request    = Dict(:gas => ["MASS", "POS ", "FRAC", "NH  ", "NHP ", "RHO ", "GZ  "]),
            var_name   = L"M_\mathrm{HI + H_2}",
            exp_factor = 10,
            unit       = u"Msun",
        )

    elseif quantity == :stellar_gas_mass

        plot_params = PlotParams(;
            request    = Dict(:gas => ["MASS", "POS ", "FRAC", "RHO "]),
            var_name   = L"M_\star^\mathrm{gas}",
            exp_factor = 10,
            unit       = u"Msun",
        )

    elseif quantity == :ode_metal_mass

        plot_params = PlotParams(;
            request  = Dict(:gas => ["MASS", "POS ", "FRAC", "NH  ", "NHP ", "RHO ", "GZ  "]),
            var_name = L"M_Z^\mathrm{gas}",
            unit     = u"Msun",
        )

    elseif quantity == :ode_metallicity

        plot_params = PlotParams(;
            request  = Dict(:gas => ["MASS", "POS ", "FRAC", "NH  ", "NHP ", "RHO ", "GZ  "]),
            var_name = L"Z_\mathrm{gas}^\star \, [\mathrm{Z_\odot}]",
        )

    elseif quantity == :dust_mass

        plot_params = PlotParams(;
            request  = Dict(:gas => ["MASS", "POS ", "FRAC", "NH  ", "NHP ", "RHO ", "GZ  "]),
            var_name = L"M_\mathrm{d}",
            unit     = u"Msun",
        )

    elseif quantity == :generic_mass

        plot_params = PlotParams(;
            request    = Dict(
                :stellar    => ["MASS", "POS "],
                :gas        => ["MASS", "POS ", "FRAC", "NH  ", "NHP ", "PRES", "RHO ", "GZ  "],
                :dark_matter       => ["MASS", "POS "],
                :black_hole => ["MASS", "POS "],
            ),
            var_name   = L"M",
            exp_factor = 10,
            unit       = u"Msun",
        )

    elseif quantity == :stellar_number

        plot_params = PlotParams(;
            request  = Dict(:stellar => ["MASS", "POS "]),
            var_name = L"\mathrm{Number \,\, of \,\, stellar \,\, particles}",
        )

    elseif quantity == :gas_number

        plot_params = PlotParams(;
            request  = Dict(:gas => ["MASS", "POS "]),
            var_name = L"\mathrm{Number \,\, of \,\, gas \,\, cells}",
        )

    elseif quantity == :dm_number

        plot_params = PlotParams(;
            request  = Dict(:dark_matter => ["MASS", "POS "]),
            var_name = L"\mathrm{Number \,\, of \,\, DM \,\, particles}",
        )

    elseif quantity == :bh_number

        plot_params = PlotParams(;
            request  = Dict(:black_hole => ["MASS", "POS "]),
            var_name = L"\mathrm{Number \,\, of \,\, BH \,\, particles}",
        )

    elseif quantity == :molecular_fraction

        plot_params = PlotParams(;
            request  = Dict(:gas => ["RHO ", "MASS", "POS ", "FRAC", "PRES", "NH  ", "NHP "]),
            var_name = L"f_\mathrm{H_2}",
        )

    elseif quantity == :br_molecular_fraction

        plot_params = PlotParams(;
            request  = Dict(:gas => ["MASS", "POS ", "PRES"]),
            var_name = L"f_\mathrm{H_2}^\mathrm{BR}",
        )

    elseif quantity == :atomic_fraction

        plot_params = PlotParams(;
            request  = Dict(
                :gas => ["RHO ", "MASS", "POS ", "FRAC", "NH  ", "NHP ", "GZ  ", "PRES"],
            ),
            var_name = L"f_\mathrm{HI}",
        )

    elseif quantity == :ionized_fraction

        plot_params = PlotParams(;
            request  = Dict(:gas => ["RHO ", "MASS", "POS ", "FRAC", "NH  ", "NHP ", "GZ  "]),
            var_name = L"f_\mathrm{HII}",
        )

    elseif quantity == :neutral_fraction

        plot_params = PlotParams(;
            request  = Dict(:gas => ["RHO ", "MASS", "POS ", "FRAC", "NH  ", "NHP ", "GZ  "]),
            var_name = L"f_\mathrm{H_I + H_2}",
        )

    elseif quantity == :molecular_neutral_fraction

        plot_params = PlotParams(;
            request  = Dict(
                :gas => ["RHO ", "MASS", "POS ", "FRAC", "NH  ", "NHP ", "PRES", "RHO ", "GZ  "],
            ),
            var_name = L"f_\mathrm{H_2} \, / f_\mathrm{n}",
        )

    elseif quantity == :ionized_neutral_fraction

        plot_params = PlotParams(;
            request  = Dict(
                :gas => ["RHO ", "MASS", "POS ", "FRAC", "NH  ", "NHP ", "PRES", "RHO ", "GZ  "],
            ),
            var_name = L"f_\mathrm{HII} \, / f_\mathrm{n}",
        )

    elseif quantity == :stellar_gas_fraction

        plot_params = PlotParams(; request  = Dict(:gas => ["RHO ", "FRAC"]), var_name = L"f_\star")

    elseif quantity == :metal_gas_fraction

        plot_params = PlotParams(;
            request  = Dict(:gas => ["RHO ", "MASS", "POS ", "FRAC", "NH  ", "NHP ", "GZ  "]),
            var_name = L"f_Z",
        )

    elseif quantity == :dust_fraction

        plot_params = PlotParams(;
            request  = Dict(:gas => ["RHO ", "MASS", "POS ", "FRAC", "NH  ", "NHP ", "GZ  "]),
            var_name = L"f_\mathrm{d}",
        )

    elseif quantity == :generic_fraction

        plot_params = PlotParams(;
            request  = Dict(
                :gas => ["RHO ", "MASS", "POS ", "FRAC", "NH  ", "NHP ", "PRES", "RHO ", "GZ  "],
            ),
            var_name = L"f",
        )

    elseif quantity == :gas_mass_density

        plot_params = PlotParams(;
            request  = Dict(:gas => ["MASS", "POS ", "RHO "]),
            var_name = L"\rho_\mathrm{gas}",
            unit     = u"Msun * kpc^-3",
        )

    elseif quantity == :hydrogen_mass_density

        plot_params = PlotParams(;
            request  = Dict(:gas => ["MASS", "POS ", "RHO "]),
            var_name = L"\rho_\mathrm{H}",
            unit     = u"Msun * kpc^-3",
        )

    elseif quantity == :gas_number_density

        plot_params = PlotParams(;
            request  = Dict(:gas => ["MASS", "POS ", "RHO "]),
            var_name = L"n_\mathrm{gas}",
            unit     = u"cm^-3",
        )

    elseif quantity == :molecular_number_density

        plot_params = PlotParams(;
            request  = Dict(:gas => ["MASS", "POS ", "FRAC", "RHO ", "PRES", "NH  ", "NHP "]),
            var_name = L"n_\mathrm{H_2}",
            unit     = u"cm^-3",
        )

    elseif quantity == :br_molecular_number_density

        plot_params = PlotParams(;
            request  = Dict(:gas => ["MASS", "POS ", "PRES", "RHO "]),
            var_name = L"n_\mathrm{H_2}^{BR}",
            unit     = u"cm^-3",
        )

    elseif quantity == :atomic_number_density

        plot_params = PlotParams(;
            request  = Dict(
                :gas => ["MASS", "POS ", "FRAC", "NH  ", "NHP ", "RHO ", "GZ  ", "PRES"],
            ),
            var_name = L"n_\mathrm{HI}",
            unit     = u"cm^-3",
        )

    elseif quantity == :ionized_number_density

        plot_params = PlotParams(;
            request  = Dict(:gas => ["MASS", "POS ", "FRAC", "NH  ", "NHP ", "RHO ", "GZ  "]),
            var_name = L"n_\mathrm{HII}",
            unit     = u"cm^-3",
        )

    elseif quantity == :neutral_number_density

        plot_params = PlotParams(;
            request  = Dict(:gas => ["MASS", "POS ", "FRAC", "NH  ", "NHP ", "RHO ", "GZ  "]),
            var_name = L"n_\mathrm{HI + H_2}",
            unit     = u"cm^-3",
        )

    elseif quantity == :stellar_area_density

        plot_params = PlotParams(;
            request  = Dict(:stellar => ["MASS", "POS "]),
            var_name = L"\Sigma_\star",
            unit     = u"Msun * pc^-2",
        )

    elseif quantity == :gas_area_density

        plot_params = PlotParams(;
            request  = Dict(:gas => ["MASS", "POS ", "RHO "]),
            var_name = L"\Sigma_\mathrm{gas}",
            unit     = u"Msun * pc^-2",
        )

    elseif quantity == :molecular_area_density

        plot_params = PlotParams(;
            request  = Dict(:gas => ["MASS", "POS ", "FRAC", "RHO ", "PRES", "NH  ", "NHP "]),
            var_name = L"\Sigma_\mathrm{H_2}",
            unit     = u"Msun * pc^-2",
        )

    elseif quantity == :br_molecular_area_density

        plot_params = PlotParams(;
            request  = Dict(:gas => ["MASS", "POS ", "PRES", "RHO "]),
            var_name = L"\Sigma_\mathrm{H_2}^\mathrm{BR}",
            unit     = u"Msun * pc^-2",
        )

    elseif quantity == :atomic_area_density

        plot_params = PlotParams(;
            request  = Dict(
                :gas => ["MASS", "POS ", "FRAC", "NH  ", "NHP ", "RHO ", "GZ  ", "PRES"],
            ),
            var_name = L"\Sigma_\mathrm{HI}",
            unit     = u"Msun * pc^-2",
        )

    elseif quantity == :ionized_area_density

        plot_params = PlotParams(;
            request  = Dict(:gas => ["MASS", "POS ", "FRAC", "NH  ", "NHP ", "RHO ", "GZ  "]),
            var_name = L"\Sigma_\mathrm{HII}",
            unit     = u"Msun * pc^-2",
        )

    elseif quantity == :neutral_area_density

        plot_params = PlotParams(;
            request  = Dict(:gas => ["MASS", "POS ", "FRAC", "NH  ", "NHP ", "RHO ", "GZ  "]),
            var_name = L"\Sigma_\mathrm{HI + H_2}",
            unit     = u"Msun * pc^-2",
        )

    elseif quantity == :sfr_area_density

        plot_params = PlotParams(;
            request  = Dict(:stellar => ["MASS", "POS ", "GAGE"]),
            var_name = L"\Sigma_\mathrm{SFR}",
            unit     = u"Msun * yr^-1 * kpc^-2",
        )

    elseif quantity == :generic_area_density

        plot_params = PlotParams(;
            request  = Dict(
                :stellar => ["MASS", "POS ", "GAGE"],
                :gas     => ["MASS", "POS ", "FRAC", "NH  ", "NHP ", "RHO ", "PRES", "GZ  "],
            ),
            var_name = L"\Sigma",
            unit     = u"Msun * pc^-2",
        )

    elseif quantity == :gas_td

        plot_params = PlotParams(;
            request  = Dict(:gas => ["MASS", "SFR "]),
            var_name = L"t_d^\mathrm{gas}",
            unit     = u"Gyr",
        )

    elseif quantity == :molecular_td

        plot_params = PlotParams(;
            request  = Dict(
                :gas => ["MASS", "POS ", "FRAC", "RHO ", "PRES", "NH  ", "NHP ", "SFR "],
            ),
            var_name = L"t_d^\mathrm{H_2}",
            unit     = u"Gyr",
        )

    elseif quantity == :br_molecular_td

        plot_params = PlotParams(;
            request  = Dict(:gas => ["MASS", "PRES", "RHO ", "SFR "]),
            var_name = L"t_d^\mathrm{H_2^\mathrm{BR}}",
            unit     = u"Gyr",
        )

    elseif quantity == :atomic_td

        plot_params = PlotParams(;
            request  = Dict(
                :gas => ["MASS", "POS ", "FRAC", "NH  ", "NHP ", "RHO ", "GZ  ", "PRES", "SFR "],
            ),
            var_name = L"t_d^\mathrm{HI}",
            unit     = u"Gyr",
        )

    elseif quantity == :ionized_td

        plot_params = PlotParams(;
            request  = Dict(
                :gas => ["MASS", "POS ", "FRAC", "NH  ", "NHP ", "RHO ", "GZ  ", "SFR "],
            ),
            var_name = L"t_d^\mathrm{HII}",
            unit     = u"Gyr",
        )

    elseif quantity == :neutral_td

        plot_params = PlotParams(;
            request  = Dict(
                :gas => ["MASS", "POS ", "FRAC", "NH  ", "NHP ", "RHO ", "GZ  ", "SFR "],
            ),
            var_name = L"t_d^\mathrm{HI + H_2}",
            unit     = u"Gyr",
        )

    elseif quantity == :gas_metallicity

        plot_params = PlotParams(;
            request  = Dict(:gas => ["MASS", "POS ", "GMET", "GZ  "]),
            var_name = L"Z_\mathrm{gas} \, [\mathrm{Z_\odot}]",
        )

    elseif quantity == :stellar_metallicity

        plot_params = PlotParams(;
            request  = Dict(:stellar => ["MASS", "POS ", "GME2", "GZ2 "]),
            var_name = L"Z_\star \, [\mathrm{Z_\odot}]",
        )

    elseif quantity ∈ GAS_ABUNDANCE

        element_string = first(split(string(quantity), "_"))

        plot_params = PlotParams(;
            request    = Dict(:gas => ["MASS", "POS ", "GMET"]),
            axis_label = L"12 + \log_{10}(\mathrm{%$element_string} \, / \, \mathrm{H})",
        )

    elseif quantity ∈ STELLAR_ABUNDANCE

        element_string = first(split(string(quantity), "_"))

        plot_params = PlotParams(;
            request    = Dict(:stellar => ["MASS", "POS ", "GME2"]),
            axis_label = L"12 + \log_{10}(\mathrm{%$element_string} \, / \, \mathrm{H})",
        )

    elseif quantity == :stellar_radial_distance

        plot_params = PlotParams(;
            request  = Dict(:stellar => ["POS "]),
            var_name = L"r",
            unit     = u"kpc",
        )

    elseif quantity == :gas_radial_distance

        plot_params = PlotParams(;
            request  = Dict(:gas => ["POS "]),
            var_name = L"r",
            unit     = u"kpc",
        )

    elseif quantity == :dm_radial_distance

        plot_params = PlotParams(;
            request  = Dict(:dark_matter => ["POS "]),
            var_name = L"r",
            unit     = u"kpc",
        )

    elseif quantity == :stellar_xy_distance

        plot_params = PlotParams(;
            request  = Dict(:stellar => ["POS "]),
            var_name = L"r_{xy}",
            unit     = u"kpc",
        )

    elseif quantity == :gas_xy_distance

        plot_params = PlotParams(;
            request  = Dict(:gas => ["POS "]),
            var_name = L"r_{xy}",
            unit     = u"kpc",
        )

    elseif quantity == :dm_xy_distance

        plot_params = PlotParams(;
            request  = Dict(:dark_matter => ["POS "]),
            var_name = L"r_{xy}",
            unit     = u"kpc",
        )

    elseif quantity == :gas_sfr

        plot_params = PlotParams(;
            request  = Dict(:gas => ["SFR "]),
            var_name = L"\mathrm{SFR_{gas}}",
            unit     = u"Msun * yr^-1",
        )

    elseif quantity == :mass_accretion

        plot_params = PlotParams(;
            request  = Dict(
                :gas        => ["ID  ", "MASS"],
                :stellar    => ["ID  ", "MASS"],
                :black_hole => ["ID  ", "MASS"],
                :group      => ["G_R_Crit200", "G_M_Crit200"],
                :tracer     => ["PAID", "TRID"],
            ),
            var_name = "Mass accretion",
            unit     = u"Msun * yr^-1",
        )

    elseif quantity == :stellar_specific_am

        plot_params = PlotParams(;
            request  = Dict(:stellar => ["MASS", "POS ", "VEL"]),
            var_name = L"j_\star",
            unit     = u"kpc^2 * s^-1",
        )

    elseif quantity == :gas_specific_am

        plot_params = PlotParams(;
            request  = Dict(:gas => ["MASS", "POS ", "VEL"]),
            var_name = L"j_\mathrm{gas}",
            unit     = u"kpc^2 * s^-1",
        )

    elseif quantity == :dm_specific_am

        plot_params = PlotParams(;
            request  = Dict(:dark_matter => ["MASS", "POS ", "VEL"]),
            var_name = L"j_\mathrm{DM}",
            unit     = u"kpc^2 * s^-1",
        )

    elseif quantity == :stellar_circularity

        # `daBandProfile` expects that the first element in the request dictionary is for the stars
        plot_params = PlotParams(;
            request  = Dict(
                :stellar    => ["MASS", "POS ", "VEL "],
                :gas        => ["MASS", "POS "],
                :dark_matter       => ["MASS", "POS "],
                :black_hole => ["MASS", "POS "],
            ),
            var_name = L"\epsilon",
        )

    elseif quantity == :stellar_vcirc

        # `daBandProfile` expects that the first element in the request dictionary is for the stars
        plot_params = PlotParams(;
            request  = Dict(
                :stellar    => ["MASS", "POS "],
                :gas        => ["MASS", "POS "],
                :dark_matter       => ["MASS", "POS "],
                :black_hole => ["MASS", "POS "],
            ),
            var_name = L"v_\mathrm{circ}",
            unit     = u"km * s^-1",
        )

    elseif quantity == :stellar_vradial

        plot_params = PlotParams(;
            request  = Dict(:stellar => ["POS ", "VEL "]),
            var_name = L"v_r",
            unit     = u"km * s^-1",
        )

    elseif quantity == :stellar_vtangential

        plot_params = PlotParams(;
            request  = Dict(:stellar => ["POS ", "VEL "]),
            var_name = L"v_\theta",
            unit     = u"km * s^-1",
        )

    elseif quantity == :stellar_vzstar

        plot_params = PlotParams(;
            request  = Dict(:stellar => ["POS ", "VEL "]),
            var_name = L"v_z \,\, \mathrm{sign}(z)",
            unit     = u"km * s^-1",
        )

    elseif quantity == :stellar_age

        plot_params = PlotParams(;
            request  = Dict(:stellar => ["GAGE"]),
            var_name = L"\mathrm{Stellar \,\, age}",
            unit     = u"Gyr",
        )

    elseif quantity == :sfr

        plot_params = PlotParams(;
            request  = Dict(:stellar => ["MASS", "POS ", "GAGE"]),
            var_name = L"\mathrm{SFR}",
            unit     = u"Msun * yr^-1",
        )

    elseif quantity == :ssfr

        plot_params = PlotParams(;
            request  = Dict(:stellar => ["MASS", "POS ", "GAGE"]),
            var_name = L"\mathrm{sSFR}",
            unit     = u"yr^-1",
        )

    elseif quantity == :observational_sfr

        plot_params = PlotParams(;
            request  = Dict(:stellar => ["MASS", "POS ", "GAGE"]),
            var_name = L"\mathrm{SFR}",
            unit     = u"Msun * yr^-1",
        )

    elseif quantity == :observational_ssfr

        plot_params = PlotParams(;
            request  = Dict(:stellar => ["MASS", "POS ", "GAGE"]),
            var_name = L"\mathrm{sSFR}",
            unit     = u"yr^-1",
        )

    elseif quantity == :stellar_eff

        plot_params = PlotParams(;
            request  = Dict(:stellar => ["GMAS", "GSFR", "RHOC"]),
            var_name = L"\epsilon_\mathrm{ff}^\star",
        )

    elseif quantity == :gas_eff

        plot_params = PlotParams(;
            request  = Dict(:gas => ["MASS", "SFR ", "RHO "]),
            var_name = L"\epsilon_\mathrm{ff}^\mathrm{gas}",
        )

    elseif quantity == :molecular_eff

        plot_params = PlotParams(;
            request  = Dict(
                :gas => ["MASS", "POS ", "FRAC", "RHO ", "PRES", "NH  ", "NHP ", "SFR "],
            ),
            var_name = L"\epsilon_\mathrm{ff}^\mathrm{H_2}",
        )

    elseif quantity == :br_molecular_eff

        plot_params = PlotParams(;
            request  = Dict(:gas => ["MASS", "PRES", "RHO ", "SFR "]),
            var_name = L"\epsilon_\mathrm{ff}^\mathrm{H_2^\mathrm{BR}}",
        )

    elseif quantity == :atomic_eff

        plot_params = PlotParams(;
            request  = Dict(
                :gas => ["MASS", "POS ", "FRAC", "NH  ", "NHP ", "RHO ", "GZ  ", "PRES", "SFR "],
            ),
            var_name = L"\epsilon_\mathrm{ff}^\mathrm{HI}",
        )

    elseif quantity == :ionized_eff

        plot_params = PlotParams(;
            request  = Dict(
                :gas => ["MASS", "POS ", "FRAC", "NH  ", "NHP ", "RHO ", "GZ  ", "SFR "],
            ),
            var_name = L"\epsilon_\mathrm{ff}^\mathrm{HII}",
        )

    elseif quantity == :neutral_eff

        plot_params = PlotParams(;
            request  = Dict(
                :gas => ["MASS", "POS ", "FRAC", "NH  ", "NHP ", "RHO ", "GZ  ", "SFR "],
            ),
            var_name = L"\epsilon_\mathrm{ff}^\mathrm{HI + H_2}",
        )

    elseif quantity == :temperature

        plot_params = PlotParams(;
            request    = Dict(:gas => ["MASS", "POS ", "TEMP"]),
            axis_label = L"\log_{10} \, T \, [\mathrm{K}]",
        )

    elseif quantity == :pressure

        plot_params = PlotParams(; request = Dict(:gas => ["PRES"]), var_name = L"P", unit = u"Pa")

    elseif quantity == :scale_factor

        plot_params = PlotParams(; var_name=L"a")

    elseif quantity == :redshift

        plot_params = PlotParams(; var_name=L"z")

    elseif quantity == :physical_time

        plot_params = PlotParams(; var_name=L"t", unit=u"Gyr")

    elseif quantity == :lookback_time

        plot_params = PlotParams(; var_name = L"\mathrm{Lookback \,\, time}", unit = u"Gyr")

    elseif quantity == :time_step

        plot_params = PlotParams(; var_name=L"\mathrm{Number \,\, of \,\, time \,\, steps}")

    elseif quantity == :clock_time_s

        plot_params = PlotParams(; var_name = L"\mathrm{Wallclock \,\, time}", unit = u"s")

    elseif quantity == :clock_time_percent

        plot_params = PlotParams(; axis_label=L"\mathrm{Wallclock \,\, time \, [\%]}")

    elseif quantity == :tot_clock_time_s

        plot_params =
            PlotParams(; var_name = L"\mathrm{Cumulative \,\, wallclock \,\, time}", unit = u"s")

    elseif quantity == :tot_clock_time_percent

        plot_params =
            PlotParams(; axis_label=L"\mathrm{Cumulative \,\, wallclock \,\, time \, [\%]}")

    elseif quantity == :ode_gas_it

        plot_params = PlotParams(;
            request  = Dict(:gas => ["ODIT"]),
            var_name = L"\mathrm{Integration\,\, time}",
            unit     = u"Myr",
        )

    elseif quantity == :ode_gas_tau_s

        plot_params = PlotParams(;
            request  = Dict(:gas => ["TAUS"]),
            var_name = L"\tau_\mathrm{S}",
            unit     = u"Myr",
        )

    elseif quantity == :ode_gas_eta_d

        plot_params = PlotParams(;
            request  = Dict(:gas => ["ETAD"]),
            var_name = L"\eta_\mathrm{diss}",
        )

    elseif quantity == :ode_gas_eta_i

        plot_params = PlotParams(;
            request  = Dict(:gas => ["ETAI"]),
            var_name = L"\eta_\mathrm{ion}",
        )

    elseif quantity == :ode_gas_r

        plot_params = PlotParams(; request  = Dict(:gas => ["PARR"]), var_name = L"R")

    elseif quantity == :ode_gas_cold_mf

        plot_params = PlotParams(; request  = Dict(:gas => ["COLF"]), var_name = L"c_f")

    elseif quantity == :ode_stellar_it

        plot_params = PlotParams(;
            request  = Dict(:stellar => ["ODIT"]),
            var_name = L"\mathrm{it}^\star",
            unit     = u"Myr",
        )

    elseif quantity == :ode_stellar_tau_s

        plot_params = PlotParams(;
            request  = Dict(:stellar => ["TAUS"]),
            var_name = L"\tau_\mathrm{S}^\star",
            unit     = u"Myr",
        )

    elseif quantity == :ode_stellar_eta_d

        plot_params = PlotParams(;
            request  = Dict(:stellar => ["ETAD"]),
            var_name = L"\eta_\mathrm{diss}^\star",
        )

    elseif quantity == :ode_stellar_eta_i

        plot_params = PlotParams(;
            request  = Dict(:stellar => ["ETAI"]),
            var_name = L"\eta_\mathrm{ion}^\star",
        )

    elseif quantity == :ode_stellar_r

        plot_params = PlotParams(; request  = Dict(:stellar => ["PARR"]), var_name = L"R^\star")

    elseif quantity == :ode_stellar_cold_mf

        plot_params = PlotParams(; request  = Dict(:stellar => ["COLF"]), var_name = L"c_f^\star")

    elseif quantity == :ode_stellar_gas_rho

        plot_params = PlotParams(;
            request  = Dict(:stellar => ["RHOC"]),
            var_name = L"\rho_\mathrm{gas}^\star",
            unit     = u"Msun * kpc^-3",
        )

    elseif quantity == :ode_stellar_gas_Z

        plot_params = PlotParams(;
            request  = Dict(:stellar => ["PARZ", "GMAS"]),
            var_name = L"Z_\mathrm{gas}^\star \, [\mathrm{Z_\odot}]",
        )

    elseif quantity == :ode_stellar_gas_mass

        plot_params = PlotParams(;
            request    = Dict(:stellar => ["GMAS"]),
            var_name   = L"M_\mathrm{gas}^\star",
            exp_factor = 10,
            unit       = u"Msun",
        )

    elseif quantity == :ode_stellar_gas_sfr

        plot_params = PlotParams(;
            request  = Dict(:stellar => ["GSFR"]),
            var_name = L"\mathrm{SFR}_\mathrm{gas}^\star",
            unit     = u"Msun * yr^-1",
        )

    elseif quantity == :ode_stellar_gas_P

        plot_params = PlotParams(;
            request  = Dict(:stellar => ["GPRE"]),
            var_name = L"P^\star",
            unit     = u"Pa",
        )

    else

        throw(ArgumentError("plotParams: I don't recognize the quantity :$(quantity)"))

    end

    return plot_params

end
