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

@doc raw"""
    KSLawBigiel2008(
        ΣH::Vector{<:SurfaceDensity};
        <keyword arguments>
    )::Vector{<:Number}

Evaluate the Kennicutt-Schmidt law for the molecular or neutral gas, taken from Bigiel et al. (2008).

From Bigiel et al. (2008) (Section 3.1), we have

```math
\Sigma_\mathrm{SFR} = a \left( \frac{\Sigma_\mathrm{HI, H_2, gas}}{10 \, \mathrm{M_\odot \, pc^{-2}}} \right)^{\!N} \, ,
```
where N is the power-law index, and $A = \log_{10}(a)$, where $a$ is $\Sigma_\mathrm{SFR}$ at the fiducial gas surface density of $10 \, \mathrm{M_\odot \, pc^{-2}}$.

# Arguments

  - `ΣH::Vector{<:SurfaceDensity}`: Values of the molecular or neutral gas surface density.
  - `molecular::Bool=true`: If `ΣH` is the surface density of molecular hydrogen, or of neutral hydrogen.
  - `log_output::Bool=true`: If the output will the $\log_{10} \, \Sigma_\text{SFR}$, or $\Sigma_\text{SFR}$. If `log_output` = true, the implied units are $\mathrm{M_\odot \, yr^{-1} \, kpc^{-2}}$

# Returns

  - The SFR surface density.

# References

F. Bigiel et al. (2008). *THE STAR FORMATION LAW IN NEARBY GALAXIES ON SUB-KPC SCALES*. The Astrophysical Journal, **136(6)**, 2846. [doi:10.1088/0004-6256/136/6/2846](https://doi.org/10.1088/0004-6256/136/6/2846)
"""
function KSLawBigiel2008(
    ΣH::Vector{<:SurfaceDensity};
    molecular::Bool=true,
    log_output::Bool=true,
)::Vector{<:Number}

    log10ΣH = @. log10(uconvert(Unitful.NoUnits, ΣH / 10.0u"Msun * pc^-2"))

    if molecular
        log10Σsfr = @. A_BIGIEL2008_BF_MOLECULAR + log10ΣH * N_BIGIEL2008_BF_MOLECULAR
    else
        log10Σsfr = @. A_BIGIEL2008_NEUTRAL + log10ΣH * N_BIGIEL2008_NEUTRAL
    end

    if log_output
        return log10Σsfr
    else
        return @. exp10(log10Σsfr) * u"Msun * yr^-1 * kpc^-2"
    end

end

@doc raw"""
    invKSLawBigiel2008(
        Σsfr ::Vector{<:MassFlowDensity};
        <keyword arguments>
    )::Vector{<:Number}

Evaluate the inverse Kennicutt-Schmidt law for the molecular or neutral gas, taken from Bigiel et al. (2008).

From Bigiel et al. (2008) (Section 3.1, Eq. 2), we have

```math
\Sigma_\mathrm{SFR} = a \left( \frac{\Sigma_\mathrm{HI, H_2, gas}}{10 \, \mathrm{M_\odot \, pc^{-2}}} \right)^{\!N} \, ,
```
where N is the power-law index, and $A = \log_{10}(a)$, where $a$ is $\Sigma_\mathrm{SFR}$ at the fiducial gas surface density of $10 \, \mathrm{M_\odot \, pc^{-2}}$.

# Arguments

  - `Σsfr ::Vector{<:MassFlowDensity}`: Values of the SFR surface density.
  - `molecular::Bool=true`: If the output will be the surface density of molecular hydrogen, or of neutral hydrogen.
  - `log_output::Bool=true`: If the output will the $\log_{10} \, \Sigma_\text{H}$, or $\Sigma_\text{H}$. If `log_output` = true, the implied units are $10 \, \mathrm{M_\odot \, pc^{-2}}$

# Returns

  - The molecular or neutral gas surface density.

# References

F. Bigiel et al. (2008). *THE STAR FORMATION LAW IN NEARBY GALAXIES ON SUB-KPC SCALES*. The Astrophysical Journal, **136(6)**, 2846. [doi:10.1088/0004-6256/136/6/2846](https://doi.org/10.1088/0004-6256/136/6/2846)
"""
function invKSLawBigiel2008(
    Σsfr::Vector{<:MassFlowDensity};
    molecular::Bool=true,
    log_output::Bool=true,
)::Vector{<:Number}

    log10Σsfr = @. log10(ustrip(u"Msun * yr^-1 * kpc^-2", Σsfr))

    if molecular
        log10ΣH = @. (log10Σsfr - A_BIGIEL2008_BF_MOLECULAR) / N_BIGIEL2008_BF_MOLECULAR
    else
        log10ΣH = @. (log10Σsfr - A_BIGIEL2008_NEUTRAL) / N_BIGIEL2008_NEUTRAL
    end

    if log_output
        return log10ΣH
    else
        return @. exp10(log10ΣH) * 10.0u"Msun * pc^-2"
    end

end

@doc raw"""
    KSLawKennicutt1998(Σgas::Vector{<:SurfaceDensity}; <keyword arguments>)::Vector{<:Number}

Evaluate the Kennicutt-Schmidt law, taken from Kennicutt (1998).

From Kennicutt (1998) (Section 4, Eq. 4), we have

```math
\Sigma_\mathrm{SFR} = a \left( \frac{\Sigma_\mathrm{gas}}{1 \, \mathrm{M_\odot \, pc^{-2}}} \right)^{\!N} \mathrm{M_\odot \, yr^{-1] \, kpc^{-2}} \, ,
```
where N is the power-law index and $a$ is $\Sigma_\mathrm{SFR}$ at the fiducial gas surface density of $1 \, \mathrm{M_\odot \, pc^{-2}}$.

# Arguments

  - `Σgas::Vector{<:SurfaceDensity}`: Values of the gas mass surface density.
  - `log_output::Bool=true`: If the output will the $\log_{10} \, \Sigma_\text{SFR}$, or $\Sigma_\text{SFR}$. If `log_output` = true, the implied units are $\mathrm{M_\odot \, yr^{-1} \, kpc^{-2}}$

# Returns

  - The SFR surface density.

# References

R. C. Kennicutt (1998). *The Global Schmidt Law in Star-forming Galaxies*. The Astrophysical Journal, **498(2)**, 541-552. [doi:10.1086/305588](https://doi.org/10.1086/305588)
"""
function KSLawKennicutt1998(Σgas::Vector{<:SurfaceDensity}; log_output::Bool=true)::Vector{<:Number}

    log10Σgas = @. log10(ustrip(u"Msun * pc^-2", Σgas))
    log10Σsfr = @. log10(a_KS98) + log10Σgas * N_KS98

    if log_output
        return log10Σsfr
    else
        return @. exp10(log10Σsfr) * u"Msun * yr^-1 * kpc^-2"
    end

end

@doc raw"""
    invKSLawKennicutt1998(
        Σsfr::Vector{<:MassFlowDensity};
        <keyword arguments>
    )::Vector{<:Number}

Evaluate the inverse Kennicutt-Schmidt law, taken from Kennicutt (1998).

From Kennicutt (1998) (Section 4, Eq. 4), we have

```math
\Sigma_\mathrm{SFR} = a \left( \frac{\Sigma_\mathrm{gas}}{1 \, \mathrm{M_\odot \, pc^{-2}}} \right)^{\!N} \mathrm{M_\odot \, yr^{-1] \, kpc^{-2}} \, ,
```
where N is the power-law index and $a$ is $\Sigma_\mathrm{SFR}$ at the fiducial gas surface density of $1 \, \mathrm{M_\odot \, pc^{-2}}$.

# Arguments

  - `Σsfr::Vector{<:MassFlowDensity}`: Values of the SFR surface density..
  - `log_output::Bool=true`: If the output will the $\log_{10} \, \Sigma_\text{gas}$, or $\Sigma_\text{gas}$. If `log_output` = true, the implied units are $\mathrm{M_\odot \, pc^{-2}}$

# Returns

  - The gas mass surface density.

# References

R. C. Kennicutt (1998). *The Global Schmidt Law in Star-forming Galaxies*. The Astrophysical Journal, **498(2)**, 541-552. [doi:10.1086/305588](https://doi.org/10.1086/305588)
"""
function invKSLawKennicutt1998(
    Σsfr::Vector{<:MassFlowDensity};
    log_output::Bool=true,
)::Vector{<:Number}

    log10Σsfr = @. log10(ustrip(u"Msun * yr^-1 * kpc^-2", Σsfr))
    log10Σgas = @. (log10Σsfr - log10(a_KS98)) / N_KS98

    if log_output
        return log10Σgas
    else
        return @. exp10(log10Σgas) * u"Msun * pc^-2"
    end

end

"""
    derivedQtyPlotParams(magnitude::Symbol, component::Symbol)::PlotParams

Return the plotting parameters for a given `magnitude` of `component`.

# Arguments

  - `magnitude::Symbol`: One of the physical magnitudes in [`MAGNITUDES`](@ref).
  - `component::Symbol`: One of the physical components in [`COMPONENTS`](@ref).

# Returns

  - A [`PlotParams`](@ref) object, with entries:

      + `request::Dict{Symbol,Vector{String}}` -> Data request for [`readSnapshot`](@ref).
      + `var_name::AbstractString`             -> Name of the quantity for the plot axis.
      + `exp_factor::Int`                      -> Numerical exponent to scale down the axis.
      + `unit::Unitful.Units`                  -> Target unit for the axis.
      + `axis_label::AbstractString`           -> Label for the axis.
"""
function derivedQtyPlotParams(magnitude::Symbol, component::Symbol)::PlotParams

    ################################################################################################
    # Physical components in a simulation
    ################################################################################################

    ############################
    # Particle-based components
    ############################

    if component == :stellar

        c_label = "\\star"
        cp_type = component

    elseif component == :dark_matter

        c_label = "\\text{DM}"
        cp_type = component

    elseif component == :black_hole

        c_label = "\\text{BH}"
        cp_type = component

    elseif component == :Z_stellar

        c_label = "\\text{Z\\!\\star}"
        cp_type = :stellar

    #######################
    # Gas-based components
    #######################

    elseif component == :gas

        c_label = "\\text{gas}"
        cp_type = component

    elseif component == :hydrogen

        c_label = "\\text{H}"
        cp_type = :gas

    elseif component == :helium

        c_label = "\\text{He}"
        cp_type = :gas

    elseif component == :Z_gas

        c_label = "\\text{Z\\,gas}"
        cp_type = :gas

    elseif component == :ionized

        c_label = "\\text{HII}"
        cp_type = :gas

    elseif component == :neutral

        c_label = "\\mathrm{HI + H2}"
        cp_type = :gas

    elseif component == :br_atomic

        c_label = "\\text{BR HI}"
        cp_type = :gas

    elseif component == :br_molecular

        c_label = "\\text{BR H2}"
        cp_type = :gas

    ###############################
    # Components from our SF model
    ###############################

    elseif component == :ode_ionized

        c_label = "\\text{ODE i}"
        cp_type = :gas

    elseif component == :ode_atomic

        c_label = "\\text{ODE a}"
        cp_type = :gas

    elseif component == :ode_molecular

        c_label = "\\text{ODE m}"
        cp_type = :gas

    elseif component == :ode_stellar

        c_label = "\\text{ODE s}"
        cp_type = :gas

    elseif component == :ode_metals

        c_label = "\\text{ODE Z}"
        cp_type = :gas

    elseif component == :ode_dust

        c_label = "\\text{ODE d}"
        cp_type = :gas

    elseif component == :ode_neutral

        c_label = "\\mathrm{ODE a + m}"
        cp_type = :gas

    ######################
    # Wild card component
    ######################

    elseif component == :generic

        c_label = ""
        cp_type = nothing

    else

        throw(ArgumentError("derivedQtyPlotParams: I don't recognize the component :$(component)"))

    end

    ################################################################################################
    # Physical magnitudes
    ################################################################################################

    #######################
    # Mass-like magnitudes
    #######################

    if magnitude == :mass

        # See computeMass() in ./src/analysis/compute_quantities/masses.jl
        if component ∈ [:stellar, :dark_matter, :black_hole, :gas]
            request = Dict(component => ["MASS"])
        elseif component == :Z_stellar
            request = Dict(:stellar => ["MASS", "GZ2 "])
        elseif component ∈ [:hydrogen, :helium]
            request = Dict(:gas => ["MASS"])
        elseif component == :Z_gas
            request = Dict(:gas => ["MASS", "GZ  "])
        elseif component ∈ [:ionized, :neutral]
            request = Dict(:gas => ["MASS", "NH  ", "NHP "])
        elseif component ∈ [:br_atomic, :br_molecular]
            request = Dict(:gas => ["MASS", "NH  ", "NHP ", "PRES"])
        elseif component ∈ [:ode_ionized, :ode_atomic, :ode_metals, :ode_dust, :ode_neutral]
            request = Dict(:gas => ["MASS", "NH  ", "NHP ", "FRAC", "RHO ", "GZ  "])
        elseif component ∈ [:ode_molecular, :ode_stellar]
            request = Dict(:gas => ["MASS", "FRAC", "RHO "])
        else
            # Generic component
            request = Dict(
                :gas         => ["MASS", "NH  ", "NHP ", "FRAC", "RHO ", "GZ  ", "PRES"],
                :stellar     => ["MASS", "GZ2 "],
                :dark_matter => ["MASS"],
                :black_hole  => ["MASS"],
            )
        end

        if isempty(c_label)
            var_name = L"M"
        else
            var_name = L"M_%$(c_label)"
        end

        exp_factor = 10

        unit = u"Msun"

    elseif magnitude == :mass_density

        # See computeMassDensity() in ./src/analysis/compute_quantities/masses.jl
        if component ∈ [:gas, :hydrogen, :helium]
            request = Dict(:gas => ["MASS", "RHO "])
        elseif component == :Z_gas
            request = Dict(:gas => ["MASS", "RHO ", "GZ  "])
        elseif component ∈ [:ionized, :neutral]
            request = Dict(:gas => ["MASS", "RHO ", "NH  ", "NHP "])
        elseif component ∈ [:br_atomic, :br_molecular]
            request = Dict(:gas => ["MASS", "RHO ", "NH  ", "NHP ", "PRES"])
        elseif component ∈ [:ode_ionized, :ode_atomic, :ode_metals, :ode_dust, :ode_neutral]
            request = Dict(:gas => ["MASS", "NH  ", "NHP ", "FRAC", "RHO ", "GZ  "])
        elseif component ∈ [:ode_molecular, :ode_stellar]
            request = Dict(:gas => ["MASS", "FRAC", "RHO "])
        elseif component == :generic
            request = Dict(:gas => ["MASS", "NH  ", "NHP ", "FRAC", "RHO ", "GZ  ", "PRES"])
        else
            throw(ArgumentError("derivedQtyPlotParams: `component` for magnitude \
            :$(magnitude) can only be one of the gas elements of `COMPONENTS` \
            (see `./src/constants/globals.jl`), but I got :$(component)"))
        end

        if isempty(c_label)
            var_name = L"\rho"
        else
            var_name = L"\rho_%$(c_label)"
        end

        exp_factor = 0

        unit = u"Msun * kpc^-3"

    elseif magnitude == :number_density

        # See computeNumberDensity() in ./src/analysis/compute_quantities/masses.jl
        if component ∈ [:gas, :hydrogen, :helium]
            request = Dict(:gas => ["MASS", "RHO "])
        elseif component == :Z_gas
            request = Dict(:gas => ["MASS", "RHO ", "GZ  "])
        elseif component ∈ [:ionized, :neutral]
            request = Dict(:gas => ["MASS", "RHO ", "NH  ", "NHP "])
        elseif component ∈ [:br_atomic, :br_molecular]
            request = Dict(:gas => ["MASS", "RHO ", "NH  ", "NHP ", "PRES"])
        elseif component ∈ [:ode_ionized, :ode_atomic, :ode_metals, :ode_dust, :ode_neutral]
            request = Dict(:gas => ["MASS", "NH  ", "NHP ", "FRAC", "RHO ", "GZ  "])
        elseif component ∈ [:ode_molecular, :ode_stellar]
            request = Dict(:gas => ["MASS", "FRAC", "RHO "])
        elseif component == :generic
            request = Dict(:gas => ["MASS", "NH  ", "NHP ", "FRAC", "RHO ", "GZ  ", "PRES"])
        else
            throw(ArgumentError("derivedQtyPlotParams: `component` for magnitude \
            :$(magnitude) can only be one of the gas elements of `COMPONENTS` \
            (see `./src/constants/globals.jl`), but I got :$(component)"))
        end

        if isempty(c_label)
            var_name = L"n"
        else
            var_name = L"n_{\,%$(c_label)}"
        end

        exp_factor = 0

        unit = u"cm^-3"

    elseif magnitude == :area_density

        # See density2DProjection() in ./src/analysis/compute_quantities/masses.jl
        if component ∈ [:stellar, :dark_matter, :black_hole, :gas]
            request = Dict(component => ["MASS", "POS ", "RHO "])
        elseif component == :Z_stellar
            request = Dict(:stellar => ["MASS", "GZ2 ", "POS ", "RHO "])
        elseif component ∈ [:hydrogen, :helium]
            request = Dict(:gas => ["MASS", "POS ", "RHO "])
        elseif component == :Z_gas
            request = Dict(:gas => ["MASS", "GZ  ", "POS ", "RHO "])
        elseif component ∈ [:ionized, :neutral]
            request = Dict(:gas => ["MASS", "NH  ", "NHP ", "POS ", "RHO "])
        elseif component ∈ [:br_atomic, :br_molecular]
            request = Dict(:gas => ["MASS", "NH  ", "NHP ", "PRES", "POS ", "RHO "])
        elseif component ∈ [:ode_ionized, :ode_atomic, :ode_metals, :ode_dust, :ode_neutral]
            request = Dict(:gas => ["MASS", "NH  ", "NHP ", "FRAC", "RHO ", "GZ  ", "POS "])
        elseif component ∈ [:ode_molecular, :ode_stellar]
            request = Dict(:gas => ["MASS", "FRAC", "RHO ", "POS "])
        else
            # Generic component
            request = Dict(
                :gas         => ["MASS", "NH  ", "NHP ", "FRAC", "RHO ", "GZ  ", "PRES", "POS "],
                :stellar     => ["MASS", "GZ2 ", "POS ", "RHO "],
                :dark_matter => ["MASS", "POS ", "RHO "],
                :black_hole  => ["MASS", "POS ", "RHO "],
            )
        end

        if isempty(c_label)
            var_name = L"\Sigma"
        else
            var_name = L"\Sigma_%$(c_label)"
        end

        exp_factor = 0

        unit = u"Msun * pc^-2"

    elseif magnitude == :number

        # See computeNumber() in ./src/analysis/compute_quantities/masses.jl
        if component ∈ [:stellar, :dark_matter, :black_hole, :gas]
            request = Dict(component => ["MASS"])
        elseif component == :Z_stellar
            request = Dict(:stellar => ["MASS"])
        elseif component ∈ [:hydrogen, :helium]
            request = Dict(:gas => ["MASS"])
        elseif component == :Z_gas
            request = Dict(:gas => ["MASS", "GZ  "])
        elseif component ∈ [:ionized, :neutral]
            request = Dict(:gas => ["MASS", "NH  ", "NHP "])
        elseif component ∈ [:br_atomic, :br_molecular]
            request = Dict(:gas => ["MASS", "NH  ", "NHP ", "PRES"])
        elseif component ∈ [:ode_ionized, :ode_atomic, :ode_metals, :ode_dust, :ode_neutral]
            request = Dict(:gas => ["MASS", "NH  ", "NHP ", "FRAC", "RHO ", "GZ  "])
        elseif component ∈ [:ode_molecular, :ode_stellar]
            request = Dict(:gas => ["MASS", "FRAC", "RHO "])
        else
            request = Dict(
                :gas         => ["MASS", "NH  ", "NHP ", "FRAC", "RHO ", "GZ  ", "PRES"],
                :stellar     => ["MASS"],
                :dark_matter => ["MASS"],
                :black_hole  => ["MASS"],
            )
        end

        if isempty(c_label)
            var_name = L"N"
        else
            var_name = L"N_%$(c_label)"
        end

        exp_factor = 0

        unit = Unitful.NoUnits

    elseif magnitude == :fraction

        # See computeFraction() in ./src/analysis/compute_quantities/masses.jl
        if component == :Z_stellar
            request = Dict(:stellar => ["GZ2 "])
        elseif component ∈ [:gas, :hydrogen, :helium]
            request = Dict(:gas => ["MASS"])
        elseif component == :Z_gas
            request = Dict(:gas => ["MASS", "GZ  "])
        elseif component ∈ [:ionized, :neutral]
            request = Dict(:gas => ["MASS", "NH  ", "NHP "])
        elseif component ∈ [:br_atomic, :br_molecular]
            request = Dict(:gas => ["MASS", "NH  ", "NHP ", "PRES"])
        elseif component ∈ [:ode_ionized, :ode_atomic, :ode_metals, :ode_dust, :ode_neutral]
            request = Dict(:gas => ["MASS", "NH  ", "NHP ", "FRAC", "RHO ", "GZ  "])
        elseif component ∈ [:ode_molecular, :ode_stellar]
            request = Dict(:gas => ["MASS", "FRAC", "RHO "])
        elseif component == :generic
            request = Dict(
                :gas     => ["MASS", "NH  ", "NHP ", "FRAC", "RHO ", "GZ  ", "PRES"],
                :stellar => ["GZ2 "],
            )
        else
            throw(ArgumentError("derivedQtyPlotParams: `component` for the magnitude \
            :$(magnitude) can only be one of the gas elements of `COMPONENTS` \
            (see `./src/constants/globals.jl`) or :Z_stellar, but I got :$(component)"))
        end

        if isempty(c_label)
            var_name = L"f"
        else
            var_name = L"f_{\,%$(c_label)}"
        end

        exp_factor = 0

        unit = Unitful.NoUnits

    elseif magnitude == :eff

        # See computeEfficiencyFF() in ./src/analysis/compute_quantities/masses.jl
        if component == :stellar
            request = Dict(:stellar => ["RHOC", "GMAS", "GSFR"])
        elseif component ∈ [:gas, :hydrogen, :helium]
            request = Dict(:gas => ["SFR ", "MASS", "RHO "])
        elseif component == :Z_gas
            request = Dict(:gas => ["SFR ", "MASS", "RHO ", "GZ  "])
        elseif component ∈ [:ionized, :neutral]
            request = Dict(:gas => ["SFR ", "MASS", "RHO ", "NH  ", "NHP "])
        elseif component ∈ [:br_atomic, :br_molecular]
            request = Dict(:gas => ["SFR ", "MASS", "RHO ", "NH  ", "NHP ", "PRES"])
        elseif component ∈ [:ode_ionized, :ode_atomic, :ode_metals, :ode_dust, :ode_neutral]
            request = Dict(:gas => ["SFR ", "MASS", "NH  ", "NHP ", "FRAC", "RHO ", "GZ  "])
        elseif component ∈ [:ode_molecular, :ode_stellar]
            request = Dict(:gas => ["SFR ", "MASS", "FRAC", "RHO "])
        elseif component == :generic
            request = Dict(
                :gas => ["SFR ", "MASS", "NH  ", "NHP ", "FRAC", "RHO ", "GZ  ", "PRES"],
                :stellar => ["RHOC", "GMAS", "GSFR"],
            )
        else
            throw(ArgumentError("derivedQtyPlotParams: `component` for magnitude \
            :$(magnitude) can only be one of the gas elements of `COMPONENTS` \
            (see `./src/constants/globals.jl`), but I got :$(component)"))
        end

        if isempty(c_label)
            var_name = L"\epsilon_\text{ff}"
        else
            var_name = L"\epsilon_{\text{ff}, %$(c_label)}"
        end

        exp_factor = 0

        unit = Unitful.NoUnits

    elseif magnitude == :clumping_factor

        # See computeClumpingFactor() in ./src/analysis/compute_quantities/masses.jl
        if component ∈ [:gas, :hydrogen, :helium]
            request = Dict(:gas => ["MASS", "RHO "])
        elseif component == :Z_gas
            request = Dict(:gas => ["MASS", "RHO ", "GZ  "])
        elseif component ∈ [:ionized, :neutral]
            request = Dict(:gas => ["MASS", "RHO ", "NH  ", "NHP "])
        elseif component ∈ [:br_atomic, :br_molecular]
            request = Dict(:gas => ["MASS", "RHO ", "NH  ", "NHP ", "PRES"])
        elseif component ∈ [:ode_ionized, :ode_atomic, :ode_metals, :ode_dust, :ode_neutral]
            request = Dict(:gas => ["MASS", "NH  ", "NHP ", "FRAC", "RHO ", "GZ  "])
        elseif component ∈ [:ode_molecular, :ode_stellar]
            request = Dict(:gas => ["MASS", "FRAC", "RHO "])
        elseif component == :generic
            request = Dict(:gas => ["MASS", "NH  ", "NHP ", "FRAC", "RHO ", "GZ  ", "PRES"])
        else
            throw(ArgumentError("derivedQtyPlotParams: `component` for magnitude \
            :$(magnitude) can only be one of the gas elements of `COMPONENTS` \
            (see `./src/constants/globals.jl`), but I got :$(component)"))
        end

        if isempty(c_label)
            var_name = L"C_\rho"
        else
            var_name = L"C_{\rho, \,%$(c_label)}"
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

        if component ∈ types
            request = Dict(component => blocks)
        elseif component == :generic
            request = Dict(type => blocks for type in types)
        elseif component == :Z_stellar
            request = Dict(:stellar => blocks)
        else
            # Gas-based components
            request = Dict(:gas => blocks)
        end

        if isempty(c_label)
            var_name = L"j_z"
        else
            var_name = L"j_{z, \,%$(c_label)}"
        end

        exp_factor = 0

        unit = u"kpc^2 * s^-1"

    elseif magnitude == :z_angular_momentum

        # See computeAngularMomentum() in ./src/analysis/compute_quantities/velocities.jl
        blocks = ["VEL ", "POS ", "MASS"]
        types  = [:gas, :dark_matter, :stellar, :black_hole]

        if component ∈ types
            request = Dict(component => blocks)
        elseif component == :generic
            request = Dict(type => blocks for type in types)
        elseif component == :Z_stellar
            request = Dict(:stellar => blocks)
        else
            # Gas-based components
            request = Dict(:gas => blocks)
        end

        if isempty(c_label)
            var_name = L"l_z"
        else
            var_name = L"l_{z, \,%$(c_label)}"
        end

        exp_factor = 0

        unit = u"Msun * kpc^2 * s^-1"

    elseif magnitude == :spin_parameter

        # See computeSpinParameter() in ./src/analysis/compute_quantities/velocities.jl
        blocks = ["VEL ", "POS ", "MASS"]
        types  = [:gas, :dark_matter, :stellar, :black_hole]

        if component ∈ types
            request = Dict(component => blocks)
        elseif component == :generic
            request = Dict(type => blocks for type in types)
        elseif component == :Z_stellar
            request = Dict(:stellar => blocks)
        else
            # Gas-based components
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
        request = Dict(type => ["VEL ", "POS ", "MASS"] for type in keys(PARTICLE_INDEX))

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
            var_name = L"v_{\text{circ}, \,%$(c_label)}"
        end

        exp_factor = 0

        unit = u"km * s^-1"

    elseif magnitude == :radial_velocity

        # See computeVpolar() in ./src/analysis/compute_quantities/velocities.jl
        blocks = ["VEL ", "POS "]
        types  = [:gas, :dark_matter, :stellar, :black_hole]

        if component ∈ types
            request = Dict(component => blocks)
        elseif component == :generic
            request = Dict(type => blocks for type in types)
        elseif component == :Z_stellar
            request = Dict(:stellar => blocks)
        else
            # Gas-based components
            request = Dict(:gas => blocks)
        end

        if isempty(c_label)
            var_name = L"v_r"
        else
            var_name = L"v_{r, \,%$(c_label)}"
        end

        exp_factor = 0

        unit = u"km * s^-1"

    elseif magnitude == :tangential_velocity

        # See computeVpolar() in ./src/analysis/compute_quantities/velocities.jl
        blocks = ["VEL ", "POS "]
        types  = [:gas, :dark_matter, :stellar, :black_hole]

        if component ∈ types
            request = Dict(component => blocks)
        elseif component == :generic
            request = Dict(type => blocks for type in types)
        elseif component == :Z_stellar
            request = Dict(:stellar => blocks)
        else
            # Gas-based components
            request = Dict(:gas => blocks)
        end

        if isempty(c_label)
            var_name = L"v_\theta"
        else
            var_name = L"v_{\theta, \,%$(c_label)}"
        end

        exp_factor = 0

        unit = u"km * s^-1"

    elseif magnitude == :zstar_velocity

        # See computeVpolar() in ./src/analysis/compute_quantities/velocities.jl
        blocks = ["VEL ", "POS "]
        types  = [:gas, :dark_matter, :stellar, :black_hole]

        if component ∈ types
            request = Dict(component => blocks)
        elseif component == :generic
            request = Dict(type => blocks for type in types)
        else
            # Gas-based components
            request = Dict(:gas => blocks)
        end

        if isempty(c_label)
            var_name = L"v_z \, \mathrm{sign}(z)"
        else
            var_name = L"v_{z, \,%$(c_label)} \, \mathrm{sign}(z)"
        end

        exp_factor = 0

        unit = u"km * s^-1"

    elseif magnitude == :kinetic_energy

        # See computeKineticEnergy() in ./src/analysis/compute_quantities/energies.jl
        if component ∈ [:stellar, :dark_matter, :black_hole, :gas]
            request = Dict(component => ["VEL ", "MASS"])
        elseif component == :Z_stellar
            request = Dict(:stellar => ["VEL ", "MASS", "GZ2 "])
        elseif component ∈ [:hydrogen, :helium]
            request = Dict(:gas => ["VEL ", "MASS"])
        elseif component == :Z_gas
            request = Dict(:gas => ["VEL ", "MASS", "GZ  "])
        elseif component ∈ [:ionized, :neutral]
            request = Dict(:gas => ["VEL ", "MASS", "NH  ", "NHP "])
        elseif component ∈ [:br_atomic, :br_molecular]
            request = Dict(:gas => ["VEL ", "MASS", "NH  ", "NHP ", "PRES"])
        elseif component ∈ [:ode_ionized, :ode_atomic, :ode_metals, :ode_dust, :ode_neutral]
            request = Dict(:gas => ["VEL ", "MASS", "NH  ", "NHP ", "FRAC", "RHO ", "GZ  "])
        elseif component ∈ [:ode_molecular, :ode_stellar]
            request = Dict(:gas => ["VEL ", "MASS", "FRAC", "RHO "])
        else
            # Generic component
            request = Dict(
                :gas         => ["VEL ", "MASS", "NH  ", "NHP ", "FRAC", "RHO ", "GZ  ", "PRES"],
                :stellar     => ["VEL ", "MASS", "GZ2 "],
                :dark_matter => ["VEL ", "MASS"],
                :black_hole  => ["VEL ", "MASS"],
            )
        end

        if isempty(c_label)
            var_name = L"E_k"
        else
            var_name = L"E_{k, \,%$(c_label)}"
        end

        exp_factor = 51

        unit = u"erg"

    elseif magnitude == :potential_energy

        # See computePotentialEnergy() in ./src/analysis/compute_quantities/energies.jl
        if component ∈ [:stellar, :dark_matter, :black_hole, :gas]
            request = Dict(component => ["POT ", "MASS"])
        elseif component == :Z_stellar
            request = Dict(:stellar => ["POT ", "MASS", "GZ2 "])
        elseif component ∈ [:hydrogen, :helium]
            request = Dict(:gas => ["POT ", "MASS"])
        elseif component == :Z_gas
            request = Dict(:gas => ["POT ", "MASS", "GZ  "])
        elseif component ∈ [:ionized, :neutral]
            request = Dict(:gas => ["POT ", "MASS", "NH  ", "NHP "])
        elseif component ∈ [:br_atomic, :br_molecular]
            request = Dict(:gas => ["POT ", "MASS", "NH  ", "NHP ", "PRES"])
        elseif component ∈ [:ode_ionized, :ode_atomic, :ode_metals, :ode_dust, :ode_neutral]
            request = Dict(:gas => ["POT ", "MASS", "NH  ", "NHP ", "FRAC", "RHO ", "GZ  "])
        elseif component ∈ [:ode_molecular, :ode_stellar]
            request = Dict(:gas => ["POT ", "MASS", "FRAC", "RHO "])
        else
            # Generic component
            request = Dict(
                :gas         => ["POT ", "MASS", "NH  ", "NHP ", "FRAC", "RHO ", "GZ  ", "PRES"],
                :stellar     => ["POT ", "MASS", "GZ2 "],
                :dark_matter => ["POT ", "MASS"],
                :black_hole  => ["POT ", "MASS"],
            )
        end

        if isempty(c_label)
            var_name = L"E_p"
        else
            var_name = L"E_{p, \,%$(c_label)}"
        end

        exp_factor = 51

        unit = u"erg"

    elseif magnitude == :total_energy

        # See computeTotalEnergy() in ./src/analysis/compute_quantities/energies.jl
        if component ∈ [:stellar, :dark_matter, :black_hole, :gas]
            request = Dict(component => ["VEL ", "POT ", "MASS"])
        elseif component == :Z_stellar
            request = Dict(:stellar => ["VEL ", "POT ", "MASS", "GZ2 "])
        elseif component ∈ [:hydrogen, :helium]
            request = Dict(:gas => ["VEL ", "POT ", "MASS"])
        elseif component == :Z_gas
            request = Dict(:gas => ["VEL ", "POT ", "MASS", "GZ  "])
        elseif component ∈ [:ionized, :neutral]
            request = Dict(:gas => ["VEL ", "POT ", "MASS", "NH  ", "NHP "])
        elseif component ∈ [:br_atomic, :br_molecular]
            request = Dict(:gas => ["VEL ", "POT ", "MASS", "NH  ", "NHP ", "PRES"])
        elseif component ∈ [:ode_ionized, :ode_atomic, :ode_metals, :ode_dust, :ode_neutral]
            request = Dict(:gas => ["VEL ", "POT ", "MASS", "NH  ", "NHP ", "FRAC", "RHO ", "GZ  "])
        elseif component ∈ [:ode_molecular, :ode_stellar]
            request = Dict(:gas => ["VEL ", "POT ", "MASS", "FRAC", "RHO "])
        else
            # Generic component
            request = Dict(
                :gas         => ["VEL ", "POT ", "MASS", "NH  ", "NHP ", "FRAC", "RHO ", "GZ  ", "PRES"],
                :stellar     => ["VEL ", "POT ", "MASS", "GZ2 "],
                :dark_matter => ["VEL ", "POT ", "MASS"],
                :black_hole  => ["VEL ", "POT ", "MASS"],
            )
        end

        if isempty(c_label)
            var_name = L"E"
        else
            var_name = L"E_{%$(c_label)}"
        end

        exp_factor = 51

        unit = u"erg"

    ########
    # Other
    ########

    elseif magnitude == :depletion_time

        # See computeDepletionTime() in ./src/analysis/compute_quantities/times.jl
        if component ∈ [:gas, :hydrogen, :helium]
            request = Dict(:gas => ["SFR ", "MASS"])
        elseif component == :Z_gas
            request = Dict(:gas => ["SFR ", "MASS", "GZ  "])
        elseif component ∈ [:ionized, :neutral]
            request = Dict(:gas => ["SFR ", "MASS", "NH  ", "NHP "])
        elseif component ∈ [:br_atomic, :br_molecular]
            request = Dict(:gas => ["SFR ", "MASS", "NH  ", "NHP ", "PRES"])
        elseif component ∈ [:ode_ionized, :ode_atomic, :ode_metals, :ode_dust, :ode_neutral]
            request = Dict(:gas => ["SFR ", "MASS", "NH  ", "NHP ", "FRAC", "RHO ", "GZ  "])
        elseif component ∈ [:ode_molecular, :ode_stellar]
            request = Dict(:gas => ["SFR ", "MASS", "FRAC", "RHO "])
        else
            # Generic component
            request = Dict(:gas => ["SFR ", "MASS", "NH  ", "NHP ", "FRAC", "RHO ", "GZ  ", "PRES"])
        end

        if isempty(c_label)
            var_name = L"\tau_\text{dep}"
        else
            var_name = L"\tau_{\text{dep}, \,%$(c_label)}"
        end

        exp_factor = 0

        unit = u"Gyr"

    elseif magnitude == :xy_distance

        # See computeXYDistance() in ./src/analysis/compute_quantities/positions.jl
        blocks = ["POS "]
        types  = [:gas, :dark_matter, :stellar, :black_hole]

        if component ∈ types
            request = Dict(component => blocks)
        elseif component == :generic
            request = Dict(type => blocks for type in types)
        elseif component == :Z_stellar
            request = Dict(:stellar => blocks)
        else
            # Gas-based components
            request = Dict(:gas => blocks)
        end

        if isempty(c_label)
            var_name = L"d_{xy}"
        else
            var_name = L"d_{xy, \,%$(c_label)}"
        end

        exp_factor = 0

        unit = u"kpc"

    elseif magnitude == :radial_distance

        # See computeRadialDistance() in ./src/analysis/compute_quantities/positions.jl
        blocks = ["POS "]
        types  = [:gas, :dark_matter, :stellar, :black_hole]

        if component ∈ types
            request = Dict(component => blocks)
        elseif component == :generic
            request = Dict(type => blocks for type in types)
        elseif component == :Z_stellar
            request = Dict(:stellar => blocks)
        else
            # Gas-based components
            request = Dict(:gas => blocks)
        end

        if isempty(c_label)
            var_name = L"r"
        else
            var_name = L"r_{\,%$(c_label)}"
        end

        exp_factor = 0

        unit = u"kpc"

    else

        throw(ArgumentError("derivedQtyPlotParams: I don't recognize the magnitude :$(magnitude)"))

    end

    return PlotParams(; request, var_name, exp_factor, unit, cp_type)

end

"""
    sfmQtyPlotParams(quantity::Symbol)::PlotParams

Return the plotting parameters for a given base code quantity of our SF model.

# Arguments

  - `quantity::Symbol`: One of the quantities in [`SFM_QTY`](@ref).

# Returns

  - A [`PlotParams`](@ref) object, with entries:

      + `request::Dict{Symbol,Vector{String}}` -> Data request for [`readSnapshot`](@ref).
      + `var_name::AbstractString`             -> Name of the quantity for the plot axis.
      + `exp_factor::Int`                      -> Numerical exponent to scale down the axis.
      + `unit::Unitful.Units`                  -> Target unit for the axis.
      + `axis_label::AbstractString`           -> Label for the axis.
"""
function sfmQtyPlotParams(quantity::Symbol)::PlotParams

    (
        quantity ∈ SFM_QTY ||
        throw(ArgumentError("sfmQtyPlotParams: I don't recognize the quantity :$(quantity)"))
    )

    magnitude, component = SFM_QTY_SPLITS[quantity]

    if component == :ode_gas_

        c_label = "\\,\\text{gas}"
        cp_type = :gas

    elseif component == :ode_stellar_

        c_label = "\\star"
        cp_type = :stellar

    else

        throw(ArgumentError("sfmQtyPlotParams: I don't recognize the component \
        :$(component). The only options are :ode_gas_ and :ode_stellar_"))

    end

    if magnitude == :integration_time

        # Integration time
        var_name = L"t_{i}^{%$(c_label)}"
        unit     = u"Myr"

    elseif magnitude == :parameter_a

        # Scale factor
        var_name = L"a^{%$(c_label)}"
        unit     = Unitful.NoUnits

    elseif magnitude == :parameter_uvb

        # UVB photoionization rate
        var_name = L"\text{UVB}^{%$(c_label)}"
        unit     = u"Myr^-1"

    elseif magnitude == :parameter_lwb

        # LWB photodissociation rate
        var_name = L"\text{LWB}^{%$(c_label)}"
        unit     = u"Myr^-1"

    elseif magnitude == :tau_s

        # Star formation time parameter
        var_name = L"\tau_\text{star}^{%$(c_label)}"
        unit     = u"Myr"

    elseif magnitude == :parameter_cell_density

        # Gas density
        var_name = L"\rho_{c}^{%$(c_label)}"
        unit     = u"cm^-3"

    elseif magnitude == :parameter_metallicity

        # Gas metallicity
        var_name = L"Z^{%$(c_label)}"
        unit     = Unitful.NoUnits

    elseif magnitude == :parameter_column_height

        # Column height
        var_name = L"h^{%$(c_label)}"
        unit     = u"pc"

    elseif magnitude == :parameter_eta_d

        # Photodissociation efficiency
        var_name = L"\eta_{\,\text{diss}}^{%$(c_label)}"
        unit     = Unitful.NoUnits

    elseif magnitude == :parameter_eta_i

        # Photoionization efficiency
        var_name = L"\eta_{\,\text{ion}}^{%$(c_label)}"
        unit     = Unitful.NoUnits

    elseif magnitude == :parameter_r

        # Mass recycling fraction
        var_name = L"R^{%$(c_label)}"
        unit     = Unitful.NoUnits

    elseif magnitude == :parameter_zsn

        # Metallicity of the supernova ejecta
        var_name = L"Z_\text{SN}^{%$(c_label)}"
        unit     = Unitful.NoUnits

    elseif magnitude == :sf_flag

        # Star formation flag
        var_name = L"\text{SF}_\text{flag}^{%$(c_label)}"
        unit     = Unitful.NoUnits

    elseif magnitude == :cold_mass_frac

        # Cold gas fraction
        var_name = L"f_\text{cold}^{%$(c_label)}"
        unit     = Unitful.NoUnits

    elseif magnitude == :gas_mass

        # Parent gas cell mass
        var_name = L"M_\text{gas}^{\,%$(c_label)}"
        unit     = u"Msun"

    elseif magnitude == :gas_sfr

        # Parent gas cell SFR
        var_name = L"\text{SFR}_\text{gas}^{%$(c_label)}"
        unit     = u"Msun * yr^-1"

    elseif magnitude == :gas_pressure

        # Parent gas cell pressure
        var_name = L"P^{%$(c_label)}"
        unit     = u"Pa"

    else

        throw(ArgumentError("sfmQtyPlotParams: I don't recognize the magnitude :$(magnitude), \
        see the keys of SFM_KEYS in ./src/constants/globals.jl for options"))

    end

    request = Dict(cp_type => [SFM_KEYS[magnitude]])

    return PlotParams(; request, var_name, unit, cp_type)

end

"""
    derivedSFMQtyPlotParams(quantity::Symbol)::PlotParams

Return the plotting parameters for a given derived code quantity of our SF model.

# Arguments

  - `quantity::Symbol`: One of the quantities in [`SFM_DERIVED_QTY`](@ref).

# Returns

  - A [`PlotParams`](@ref) object, with entries:

      + `request::Dict{Symbol,Vector{String}}` -> Data request for [`readSnapshot`](@ref).
      + `var_name::AbstractString`             -> Name of the quantity for the plot axis.
      + `exp_factor::Int`                      -> Numerical exponent to scale down the axis.
      + `unit::Unitful.Units`                  -> Target unit for the axis.
      + `axis_label::AbstractString`           -> Label for the axis.
"""
function derivedSFMQtyPlotParams(quantity::Symbol)::PlotParams

    (
        quantity ∈ SFM_DERIVED_QTY ||
        throw(ArgumentError("derivedSFMQtyPlotParams: I don't recognize the quantity :$(quantity)"))
    )

    magnitude, component = SFM_DERIVED_QTY_SPLITS[quantity]

    if component == :ode_gas_

        c_label = "\\,\\text{gas}"
        cp_type = :gas

    elseif component == :ode_stellar_

        c_label = "\\star"
        cp_type = :stellar

    else

        throw(ArgumentError("derivedSFMQtyPlotParams: I don't recognize the component \
        :$(component). The only options are :ode_gas_ and :ode_stellar_"))

    end

    if magnitude == :tau_star

        # Star formation timescale
        var_name = L"\tau_\text{star}^{%$(c_label)}"
        unit     = u"Myr"
        request  = Dict(cp_type => ["RHOC"])

    elseif magnitude == :tau_rec

        # Recombination timescale
        var_name = L"\tau_\text{rec}^{%$(c_label)}"
        unit     = u"Myr"
        request  = Dict(cp_type => ["RHOC", "FRAC"])

    elseif magnitude == :tau_cond

        # Condensation timescale
        var_name = L"\tau_\text{cond}^{%$(c_label)}"
        unit     = u"Myr"
        request  = Dict(cp_type => ["RHOC", "FRAC"])

    elseif magnitude == :tau_dg

        # Dust growth timescale
        var_name = L"\tau_\text{dg}^{%$(c_label)}"
        unit     = u"Gyr"
        request  = Dict(cp_type => ["RHOC", "FRAC"])

    elseif magnitude == :tau_dc

        # Dust growth timescale
        var_name = L"\tau_\text{dc}^{%$(c_label)}"
        unit     = u"Gyr"
        request  = Dict(cp_type => ["RHOC", "FRAC"])

    elseif magnitude == :tau_ion

        # Ionization optical depth
        var_name = L"\tau_\text{ion}^{%$(c_label)}"
        unit     = Unitful.NoUnits
        request  = Dict(cp_type => ["RHOC", "PARH", "FRAC"])

    elseif magnitude == :tau_diss

        # Dissociation optical depth
        var_name = L"\tau_\text{diss}^{%$(c_label)}"
        unit     = Unitful.NoUnits
        request  = Dict(cp_type => ["RHOC", "PARH", "FRAC"])

    elseif magnitude == :S_d

        # Dust shielding factor
        var_name = L"S_{d}^{%$(c_label)}"
        unit     = Unitful.NoUnits
        request  = Dict(cp_type => ["RHOC", "PARH", "FRAC"])

    elseif magnitude == :S_H2

        # H2 self-shielding factor
        var_name = L"S_{H2}^{%$(c_label)}"
        unit     = Unitful.NoUnits
        request  = Dict(cp_type => ["RHOC", "PARH", "FRAC"])

    elseif magnitude == :equivalent_size

        # Equivalent size
        var_name = L"d_\text{eq}^{%$(c_label)}"
        unit     = u"pc"
        if cp_type == :gas
            request = Dict(cp_type => ["RHOC", "MASS"])
        else
            request = Dict(cp_type => ["RHOC", "GMAS"])
        end

    else

        throw(ArgumentError("derivedSFMQtyPlotParams: I don't recognize the magnitude \
        :$(magnitude), see the keys of SFM_DERIVED_MAGNITUDES in ./src/constants/globals.jl \
        for options"))

    end

    return PlotParams(; request, var_name, unit, cp_type)

end

"""
    plotParams(quantity::Symbol)::PlotParams

Return the plotting parameters for a given `quantity`.

# Arguments

  - `quantity::Symbol`: Target quantity. Some of the options are the quantities in [`DERIVED_QTY`](@ref), [`MAGNITUDES`](@ref), [`SFM_QTY`](@ref), [`GAS_ABUNDANCE`](@ref), and [`STELLAR_ABUNDANCE`](@ref).

# Returns

  - A [`PlotParams`](@ref) object, with entries:

      + `request::Dict{Symbol,Vector{String}}` -> Data request for [`readSnapshot`](@ref).
      + `var_name::AbstractString`             -> Name of the quantity for the plot axis.
      + `exp_factor::Int`                      -> Numerical exponent to scale down the axis.
      + `unit::Unitful.Units`                  -> Target unit for the axis.
      + `axis_label::AbstractString`           -> Label for the axis.
      + `cp_type::Union{Symbol,Nothing}`       ->  Cell/particle type corresponding to the quantity.
"""
function plotParams(quantity::Symbol)::PlotParams

    #####################
    # Derived quantities
    #####################

    if quantity ∈ DERIVED_QTY

        magnitude, component = QUANTITY_SPLITS[quantity]

        plot_params = derivedQtyPlotParams(magnitude, component)

    elseif quantity ∈ MAGNITUDES

        plot_params = derivedQtyPlotParams(quantity, :generic)

    elseif quantity ∈ SFM_QTY

        plot_params = sfmQtyPlotParams(quantity)

    elseif quantity ∈ SFM_DERIVED_QTY

        plot_params = derivedSFMQtyPlotParams(quantity)

    #########
    # Ratios
    #########

    elseif quantity ∈ RATIO_QTY

        qty_01, qty_02 = RATIO_SPLITS[quantity]

        plot_params_01 = plotParams(qty_01)
        plot_params_02 = plotParams(qty_02)

        var_name_01 = plot_params_01.var_name
        var_name_02 = plot_params_02.var_name

        request_01 = plot_params_01.request
        request_02 = plot_params_02.request

        plot_params = PlotParams(;
            request  = mergeRequests(request_01, request_02),
            var_name = LaTeXString(var_name_01 * L"\, / \," * var_name_02),
        )

    #################
    # Gas quantities
    #################

    elseif quantity == :temperature

        plot_params = PlotParams(;
            request  = Dict(:gas => ["TEMP"]),
            var_name = L"T",
            unit     = u"K",
            cp_type  = :gas,
        )

    elseif quantity == :pressure

        plot_params = PlotParams(;
            request  = Dict(:gas => ["PRES"]),
            var_name = L"P",
            unit     = u"Pa",
            cp_type  = :gas,
        )

    #################
    # SFR quantities
    #################

    elseif quantity == :sfr

        # See computeSFR() in ./src/analysis/compute_quantities/times.jl
        plot_params = PlotParams(;
            request  = Dict(:stellar => ["MASS", "GAGE"]),
            var_name = L"\text{SFR}",
            unit     = u"Msun * yr^-1",
            cp_type  = :stellar,
        )

    elseif quantity == :ssfr

        # See computeSSFR() in ./src/analysis/compute_quantities/times.jl
        plot_params = PlotParams(;
            request  = Dict(:stellar => ["GAGE"]),
            var_name = L"\text{sSFR}",
            unit     = u"Gyr^-1",
            cp_type  = :stellar,
        )

    elseif quantity == :observational_sfr

        # See computeSFR() in ./src/analysis/compute_quantities/times.jl
        plot_params = PlotParams(;
            request  = Dict(:stellar => ["MASS", "GAGE"]),
            var_name = L"\text{SFR}",
            unit     = u"Msun * yr^-1",
            cp_type  = :stellar,
        )

    elseif quantity == :observational_ssfr

        # See computeSSFR() in ./src/analysis/compute_quantities/times.jl
        plot_params = PlotParams(;
            request  = Dict(:stellar => ["MASS", "GAGE"]),
            var_name = L"\text{sSFR}",
            unit     = u"yr^-1",
            cp_type  = :stellar,
        )

    elseif quantity == :sfr_area_density

        plot_params = PlotParams(;
            request  = Dict(:stellar => ["MASS", "GAGE", "POS "]),
            var_name = L"\Sigma_\text{SFR}",
            unit     = u"Msun * yr^-1 * kpc^-2",
            cp_type  = :stellar,
        )

    elseif quantity == :sfr_density

        plot_params = PlotParams(;
            request  = Dict(:stellar => ["MASS", "GAGE", "POS "]),
            var_name = L"\rho_\text{SFR}",
            unit     = u"Msun * yr^-1 * kpc^-3",
            cp_type  = :stellar,
        )

    elseif quantity == :stellar_age

        # See computeStellarAge() in ./src/analysis/compute_quantities/times.jl
        plot_params = PlotParams(;
            request  = Dict(:stellar => ["GAGE"]),
            var_name = L"\text{Stellar age}",
            unit     = u"Gyr",
            cp_type  = :stellar,
        )

    elseif quantity == :stellar_birth_time

        # See computeStellarAge() in ./src/analysis/compute_quantities/times.jl
        plot_params = PlotParams(;
            request  = Dict(:stellar => ["GAGE"]),
            var_name = L"\text{Stellar birth time}",
            unit     = u"Gyr",
            cp_type  = :stellar,
        )

    elseif quantity == :gas_sfr

        plot_params = PlotParams(;
            request  = Dict(:gas => ["SFR "]),
            var_name = L"\mathrm{SFR_{gas}}",
            unit     = u"Msun * yr^-1",
            cp_type  = :gas,
        )

    elseif quantity == :gas_sfr_area_density

        plot_params = PlotParams(;
            request  = Dict(:gas => ["SFR ", "POS "]),
            var_name = L"\Sigma_\text{SFR}^\text{gas}",
            unit     = u"Msun * yr^-1 * kpc^-2",
            cp_type  = :gas,
        )

    ##################
    # Time quantities
    ##################

    elseif quantity == :scale_factor

        plot_params = PlotParams(; var_name=L"a")

    elseif quantity == :redshift

        plot_params = PlotParams(; var_name=L"z")

    elseif quantity == :physical_time

        plot_params = PlotParams(; var_name=L"t", unit=u"Gyr")

    elseif quantity == :lookback_time

        plot_params = PlotParams(; var_name=L"\mathrm{Lookback \,\, time}", unit=u"Gyr")

    elseif quantity == :time_step

        plot_params = PlotParams(; var_name=L"\mathrm{Number \,\, of \,\, time \,\, steps}")

    elseif quantity == :clock_time_s

        plot_params = PlotParams(; var_name=L"\mathrm{Wallclock \,\, time}", unit=u"s")

    elseif quantity == :clock_time_percent

        plot_params = PlotParams(; axis_label=L"\mathrm{Wallclock \,\, time \, [\%]}")

    elseif quantity == :tot_clock_time_s

        plot_params = PlotParams(;
            var_name = L"\mathrm{Cumulative \,\, wallclock \,\, time}",
            unit     = u"s",
        )

    elseif quantity == :tot_clock_time_percent

        plot_params = PlotParams(;
            axis_label = L"\mathrm{Cumulative \,\, wallclock \,\, time \, [\%]}"
        )

    #########################
    # Metallicity quantities
    #########################

    elseif quantity == :gas_metallicity

        plot_params = PlotParams(;
            request  = Dict(:gas => ["GZ  "]),
            var_name = L"Z_\mathrm{gas} [\mathrm{Z_\odot}]",
            cp_type  = :gas,
        )

    elseif quantity == :stellar_metallicity

        plot_params = PlotParams(;
            request  = Dict(:stellar => ["GZ2 "]),
            var_name = L"Z_\star \, [\mathrm{Z_\odot}]",
            cp_type  = :stellar,
        )

    elseif quantity == :ode_metallicity

        plot_params = PlotParams(;
            request  = Dict(:gas => ["MASS", "NH  ", "NHP ", "FRAC", "RHO ", "GZ  "]),
            var_name = L"Z_\text{ODE} \, [\mathrm{Z_\odot}]",
            cp_type  = :gas,
        )

    elseif quantity ∈ GAS_ABUNDANCE

        # See computeAbundance() in ./src/analysis/compute_quantities/masses.jl
        element = GAS_ABUNDANCE_SPLITS[quantity]

        if iszero(ABUNDANCE_SHIFT[element])
            axis_label = L"\log_{10}(\mathrm{%$element} \, / \, \mathrm{H})"
        else
            axis_label = L"%$(Int(ABUNDANCE_SHIFT[element])) + \log_{10}(\mathrm{%$element} \, / \, \mathrm{H})"
        end

        plot_params = PlotParams(;
            request    = Dict(:gas => ["MASS", "GMET"]),
            axis_label,
            cp_type    = :gas,
        )

    elseif quantity ∈ STELLAR_ABUNDANCE

        # See computeAbundance() in ./src/analysis/compute_quantities/masses.jl
        element = STELLAR_ABUNDANCE_SPLITS[quantity]

        if iszero(ABUNDANCE_SHIFT[element])
            axis_label = L"\log_{10}(\mathrm{%$element} \, / \, \mathrm{H})"
        else
            axis_label = L"%$(Int(ABUNDANCE_SHIFT[element])) + \log_{10}(\mathrm{%$element} \, / \, \mathrm{H})"
        end

        plot_params = PlotParams(;
            request    = Dict(:stellar => ["MASS", "GME2"]),
            axis_label,
            cp_type    = :stellar,
        )

    elseif quantity == :mass_accretion

        # See computeVirialAccretion() and computeDiscAccretion() in ./src/analysis/compute_quantities/masses.jl
        plot_params = PlotParams(;
            request  = Dict(
                :gas         => ["ID  ", "POS ", "MASS"],
                :stellar     => ["ID  ", "POS ", "MASS"],
                :black_hole  => ["ID  ", "POS ", "MASS"],
                :dark_matter => ["ID  ", "POS ", "MASS"],
                :group       => ["G_R_Crit200", "G_M_Crit200", "G_Nsubs", "G_Pos"],
                :subhalo     => ["S_Pos"],
                :tracer      => ["PAID", "TRID", "POS "],
            ),
            var_name = L"\dot{M}_\text{acc}",
            unit     = u"Msun * yr^-1",
        )

    else

        throw(ArgumentError("plotParams: I don't recognize the quantity :$(quantity)"))

    end

    return plot_params

end

"""
    validatePlotData(
        plot_function::Function,
        plot_data::VecOrMat{<:Number}...;
        <keyword arguments>
    )::Tuple{Vector{<:VecOrMat{<:Number}},Bool,Bool}

Adapt `plot_data` to be plotted with `plot_function`.

This function strips the units of every argument in `plot_data`, and it may trims down the first two arguments in `plot_data`.

# Arguments

  - `plot_function::Function`: Target plotting function from [Makie](https://docs.makie.org/stable/). The supported functions are:

      + `scatter!`      -> Scatter plot.
      + `lines!`        -> Line plot.
      + `scatterlines!` -> Scatter plot with lines between the markers.
      + `hist!`         -> Histogram.
      + `heatmap!`      -> Heatmap.
      + `arrows2d!`     -> 2D vector field.
      + `barplot!`      -> Bar plots.
      + `band!`         -> Band plots.
      + `errorbars!`    -> Error bars.
  - `plot_data::VecOrMat{<:Number}`: Raw plot data.
  - `x_unit::Unitful.Units=Unitful.NoUnits`: Target unit for the x axis data (first arguments in `plot_data`). The values will be converted accordingly. Use the default value of `Unitful.NoUnits` for dimensionless quantities.
  - `y_unit::Unitful.Units=Unitful.NoUnits`: Target unit for the y axis data (second arguments in `plot_data`). The values will be converted accordingly. Use the default value of `Unitful.NoUnits` for dimensionless quantities.
  - `x_exp_factor::Int=0`: Numerical exponent to scale down the x axis data (first arguments in `plot_data`), e.g. if `x_exp_factor` = 10 the values will be divided by ``10^{10}``. The default is no scaling.
  - `y_exp_factor::Int=0`: Numerical exponent to scale down the y axis data (second arguments in `plot_data`), e.g. if `y_exp_factor` = 10 the values will be divided by ``10^{10}``. The default is no scaling.
  - `x_trim::NTuple{2,<:Real}=(-Inf, Inf)`: The data (first arguments in `plot_data`) will be trim down so the x coordinates fit within `x_trim`, in the units given by `x_unit`.
  - `y_trim::NTuple{2,<:Real}=(-Inf, Inf)`: The data (second arguments in `plot_data`) will be trim down so the y coordinates fit within `y_trim`, in the units given by `y_unit`. This option does not affect histograms.
  - `x_edges::Bool=false`: Set it to `true` if you want to keep the borders of `x_trim`.
  - `y_edges::Bool=false`: Set it to `true` if you want to keep the borders of `y_trim`.
  - `x_scale_func::Function=identity`: Scaling function for the x axis (first arguments in `plot_data`). The options are the scaling functions accepted by [Makie](https://docs.makie.org/stable/): log10, log2, log, sqrt, Makie.logit, Makie.Symlog10, Makie.pseudolog10, and identity. The data will be trimmed down to fit within the domain of `x_scale_func`.
  - `y_scale_func::Function=identity`: Scaling function for the y axis (second arguments in `plot_data`). The options are the scaling functions accepted by [Makie](https://docs.makie.org/stable/): log10, log2, log, sqrt, Makie.logit, Makie.Symlog10, Makie.pseudolog10, and identity. The data will be trimmed down to fit within the domain of `y_scale_func`.

# Returns

  - A tuple with three elements:

      + A vector with the elements of `plot_data` ready to be plotted by `plot_function`.
      + If `plot_data[1]` was mutated to fit within the domain of `x_scale_func`.
      + If `plot_data[2]` was mutated to fit within the domain of `y_scale_func` (only relevant if `length(plot_data)` > 1).
"""
function validatePlotData(
    plot_function::Function,
    plot_data::VecOrMat{<:Number}...;
    x_unit::Unitful.Units=Unitful.NoUnits,
    y_unit::Unitful.Units=Unitful.NoUnits,
    x_exp_factor::Int=0,
    y_exp_factor::Int=0,
    x_trim::NTuple{2,<:Real}=(-Inf, Inf),
    y_trim::NTuple{2,<:Real}=(-Inf, Inf),
    x_edges::Bool=false,
    y_edges::Bool=false,
    x_scale_func::Function=identity,
    y_scale_func::Function=identity,
)::Tuple{Vector{<:VecOrMat{<:Number}},Bool,Bool}

    data_length = length(plot_data)

    #########################
    # Validate the data size
    #########################

    if plot_function isa typeof(hist!)

        (
            data_length == 1 ||
            throw(ArgumentError("validateDimensionalData: For histograms the plot data must \
            contain only one data vector, but currently it has $(data_length)"))
        )

    elseif plot_function isa Union{
        typeof(scatter!),
        typeof(scatterlines!),
        typeof(lines!),
        typeof(barplot!),
    }

        (
            data_length == 2 ||
            throw(ArgumentError("validateDimensionalData: For scatter, scatterlines, lines, and \
            bar plots the plot data must contain only two data vectors, but currently it has  \
            $(data_length)"))
        )

    elseif plot_function isa Union{typeof(heatmap!),typeof(band!)}

        (
            data_length == 3 ||
            throw(ArgumentError("validateDimensionalData: For heatmaps and bands plots the plot \
            data must contain only three data vectors, but currently it has $(data_length)"))
        )

    elseif plot_function isa Union{typeof(arrows2d!),typeof(errorbars!)}

        (
            data_length == 4 ||
            throw(ArgumentError("validateDimensionalData: For vector field and error bars plots \
            the plot data must contain only four data vectors, but currently it has \
            $(data_length)"))
        )

    else

        throw(ArgumentError("validateDimensionalData: I don't recognize the plot function \
        $(plot_function). I can only use: hist!, scatter!, scatterlines!, lines!, barplot!, \
        heatmap!, band!, arrows2d!, and errorbars!."))

    end

    ################
    # Convert units
    ################

    if plot_function isa typeof(hist!)

        axis_data = VecOrMat{<:Number}[ustrip.(x_unit, plot_data[1])]

    elseif plot_function isa Union{
        typeof(scatter!),
        typeof(scatterlines!),
        typeof(lines!),
        typeof(barplot!),
    }

        axis_data = VecOrMat{<:Number}[ustrip.(x_unit, plot_data[1]), ustrip.(y_unit, plot_data[2])]

    elseif plot_function isa typeof(heatmap!)

        axis_data = VecOrMat{<:Number}[
            ustrip.(x_unit, plot_data[1]),
            ustrip.(y_unit, plot_data[2]),
            ustrip.(Unitful.NoUnits, plot_data[3]),
        ]

    elseif plot_function isa typeof(band!)

        axis_data = VecOrMat{<:Number}[
            ustrip.(x_unit, plot_data[1]),
            ustrip.(y_unit, plot_data[2]),
            ustrip.(y_unit, plot_data[3]),
        ]

    elseif plot_function isa typeof(arrows2d!)

        axis_data = VecOrMat{<:Number}[
            ustrip.(x_unit, plot_data[1]),
            ustrip.(y_unit, plot_data[2]),
            ustrip.(Unitful.NoUnits, plot_data[3]),
            ustrip.(Unitful.NoUnits, plot_data[4]),
        ]

    elseif plot_function isa typeof(errorbars!)

        axis_data = VecOrMat{<:Number}[
            ustrip.(x_unit, plot_data[1]),
            ustrip.(y_unit, plot_data[2]),
            ustrip.(y_unit, plot_data[3]),
            ustrip.(y_unit, plot_data[4]),
        ]

    end

    ##################
    # Data sanitation
    ##################

    if data_length == 1

        x_flag, _ = sanitizeData!(
            axis_data[1];
            func_domain=x_scale_func,
            range=x_trim,
            keep_edges=x_edges,
            min_left=1,
            exp_factor=x_exp_factor,
        )

        y_flag = true

    else

        x_flag, y_flag, _, _ = sanitizeData!(
            axis_data[1],
            axis_data[2];
            func_domain=(x_scale_func, y_scale_func),
            range=(x_trim, y_trim),
            keep_edges=(x_edges, y_edges),
            min_left=1,
            exp_factor=(x_exp_factor, y_exp_factor),
        )

    end

    return axis_data, x_flag, y_flag

end
