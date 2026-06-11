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
    getLabel(quantity::Symbol; <keyword arguments>)::AbstractString

Construct an axis label for a given `quantity`.

!!! note

    This label will ignore the `exp_factor` field of `BaseQuantity` and take it as 0.

# Arguments

  - `quantity::Symbol`: Target quantity.
  - `log::Bool=false`: If true adds ``\\log_{10}`` to the label.
  - `latex::Bool=true`: If the output will be a `LaTeXString` or a plain `String`.

# Returns

  - The label for `quantity`, according to [`QTY_REGISTRY`](@ref).
"""
function getLabel(quantity::Symbol; log::Bool=false, latex::Bool=true)::AbstractString

    (
        haskey(QTY_REGISTRY, quantity) ||
        throw(ArgumentError("getLabel: I don't recognize the quantity :$(quantity). \
        It is not a key of QTY_REGISTRY"))
    )

    return getLabel(QTY_REGISTRY[quantity]; log, latex)

end

"""
    getLabel(quantity::AbstractPlotQuantity; <keyword arguments>)::AbstractString

Construct an axis label for a given [`AbstractPlotQuantity`](@ref).

!!! note

    This label will ignore the `exp_factor` of [`AbstractPlotQuantity`](@ref) and take it as 0.

# Arguments

  - `quantity::AbstractPlotQuantity`: Target quantity.
  - `log::Bool=false`: If true adds ``\\log_{10}`` to the label.
  - `latex::Bool=true`: If the output will be a `LaTeXString` or a plain `String`.

# Returns

  - The label for `quantity`, according to [`QTY_REGISTRY`](@ref).
"""
function getLabel(quantity::AbstractPlotQuantity; log::Bool=false, latex::Bool=true)::AbstractString

    base_label = getLabel(getQuantityLabel(quantity), 0, getQuantityUnit(quantity); latex)

    if latex && log

        label = LaTeXString(L"\log_{10} \," * base_label)

    elseif !latex && log

        label = "log10 " * base_label

    else

        label = base_label

    end

    return label

end

"""
    getLabelArgs(quantity::Symbol; <keyword arguments>)::NamedTuple

Construct the label arguments of [`plotSnapshot`](@ref) for a given `quantity`.

# Arguments

  - `quantity::Symbol`: Target quantity.
  - `log::Bool=false`: If true adds ``\\log_{10}`` to the axis label.

# Returns

  - A named tuple with

      + `log_unit::Union{Unitful.Units,Nothing}` -> Quantity unit for the functions that use to decide if the axis will be logarithmic or not. It is set to `nothing` if `log` = false, and to the corresponding unit if `log` = true.
      + `unit::Unitful.Units`                    -> Target unit for the axis data.
      + `exp_factor::Int`                        -> Numerical exponent to scale down the axis data, e.g. if `exp_factor` = 10 the values will be divided by ``10^{10}``.
      + `qty_label::AbstractString`              -> Name of the variable for the axis.
      + `label::AbstractString`                  -> Label for the axis. It can contain the string `auto_label`, which will be replaced by: `qty_label` [10^`exp_factor` `unit`]. If a LaTeXString with `auto_label` inside is used, it is recommended that each section around `auto_label` is delimited with a `\$ \$` pair.
      + `request::Dict{Symbol,Vector{String}}`   -> Data request for [`readSnapshot`](@ref).
"""
function getLabelArgs(quantity::Symbol; log::Bool=false)::NamedTuple

    (
        haskey(QTY_REGISTRY, quantity) ||
        throw(ArgumentError("getLabelArgs: I don't recognize the quantity :$(quantity). \
        It is not a key of QTY_REGISTRY"))
    )

    return getLabelArgs(QTY_REGISTRY[quantity]; log)

end

"""
    getLabelArgs(quantity::AbstractPlotQuantity; <keyword arguments>)::NamedTuple

Construct the label arguments of [`plotSnapshot`](@ref) for a given [`AbstractPlotQuantity`](@ref).

# Arguments

  - `quantity::AbstractPlotQuantity`: Target quantity.
  - `log::Bool=false`: If true adds ``\\log_{10}`` to the axis label.

# Returns

  - A named tuple with

      + `log_unit::Union{Unitful.Units,Nothing}` -> Quantity unit for the functions that use to decide if the axis will be logarithmic or not. It is set to `nothing` if `log` = false, and to the corresponding unit if `log` = true.
      + `unit::Unitful.Units`                    -> Target unit for the axis data.
      + `exp_factor::Int`                        -> Numerical exponent to scale down the axis data, e.g. if `exp_factor` = 10 the values will be divided by ``10^{10}``.
      + `qty_label::AbstractString`              -> Name of the variable for the axis.
      + `label::AbstractString`                  -> Label for the axis. It can contain the string `auto_label`, which will be replaced by: `qty_label` [10^`exp_factor` `unit`]. If a LaTeXString with `auto_label` inside is used, it is recommended that each section around `auto_label` is delimited with a `\$ \$` pair.
      + `request::Dict{Symbol,Vector{String}}`   -> Data request for [`readSnapshot`](@ref).
"""
function getLabelArgs(quantity::AbstractPlotQuantity; log::Bool=false)::NamedTuple

    if log
        log_unit   = getQuantityUnit(quantity)
        unit       = Unitful.NoUnits
        exp_factor = 0
        qty_label  = ""
        label      = getLabel(quantity; log)
    else
        log_unit   = nothing
        unit       = getQuantityUnit(quantity)
        exp_factor = getQuantityExpFactor(quantity)
        qty_label  = getQuantityLabel(quantity)
        label      = "auto_label"
    end

    request = getQuantityRequest(quantity)

    return (; log_unit, unit, exp_factor, qty_label, label, request)

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
    reduceTicks(hr_ticks::Vector{T}, factor::Int)::Vector{T} where {T<:FloatOrQuantity}

Reduce the length of a given list of axis ticks, while keeping the total length of the axis the same.

!!! note

    If `factor` is 1, the original list of ticks will be returned without any checks.

# Arguments

  - `hr_ticks::Vector{T}`: Original "high resolution" list of ticks. It has to be regularly spaced
  - `factor::Int`: Factor by which the number of values will be reduced. It has to divide the size of `hr_ticks` exactly.

# Returns

  - The new shorter list of ticks.
"""
function reduceTicks(hr_ticks::Vector{T}, factor::Int)::Vector{T} where {T<:FloatOrQuantity}

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
  - `log_output::Bool=true`: If the output will the $\log_{10} \, \Sigma_\text{H2}$, or $\Sigma_\text{HI}$. If `log_output` = true, the implied units are $10 \, \mathrm{M_\odot \, pc^{-2}}$

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

@doc raw"""
    KSLawCologni2026(ΣH2::Vector{<:SurfaceDensity}; <keyword arguments>)::Vector{<:Number}

Evaluate the Kennicutt-Schmidt law for the molecular gas, taken from Cologni et al. (2026).

From Cologni et al. (2026) (Section 4.3.2), we have

```math
\Sigma_\mathrm{SFR} = a \left( \frac{\Sigma_\mathrm{H_2}}{1.0 \, \mathrm{M_\odot \, kpc^{-2}}} \right)^{\!N} \, ,
```
where N is the power-law index, and $A = \log_{10}(a)$, where $a$ is $\Sigma_\mathrm{SFR}$ at the fiducial gas surface density of $1.0 \, \mathrm{M_\odot \, kpc^{-2}}$.

# Arguments

  - `ΣH2::Vector{<:SurfaceDensity}`: Values of the molecular gas surface density.
  - `log_output::Bool=true`: If the output will the $\log_{10} \, \Sigma_\text{SFR}$, or $\Sigma_\text{SFR}$. If `log_output` = true, the implied units are $\mathrm{M_\odot \, yr^{-1} \, kpc^{-2}}$

# Returns

  - The SFR surface density.

# References

R. Cologni et al. (2026). *Resolved molecular gas and star formation in massive unquenched spirals*. Astronomy and Astrophysics, **709**, A148. [doi:10.1051/0004-6361/202558486](https://doi.org/10.1051/0004-6361/202558486)
"""
function KSLawCologni2026(ΣH2::Vector{<:SurfaceDensity}; log_output::Bool=true)::Vector{<:Number}

    log10ΣH2 = @. log10(uconvert(Unitful.NoUnits, ΣH2 / 1.0u"Msun * kpc^-2"))

    log10Σsfr = @. A_COLOGNI2026 + log10ΣH2 * N_COLOGNI2026

    if log_output
        return log10Σsfr
    else
        return @. exp10(log10Σsfr) * u"Msun * yr^-1 * kpc^-2"
    end

end

@doc raw"""
    invKSLawCologni2026(Σsfr ::Vector{<:MassFlowDensity}; <keyword arguments>)::Vector{<:Number}

Evaluate the inverse Kennicutt-Schmidt law for the molecular gas, taken from Cologni et al. (2026).

From Cologni et al. (2026) (Section 4.3.2), we have

```math
\Sigma_\mathrm{SFR} = a \left( \frac{\Sigma_\mathrm{H_2}}{1.0 \, \mathrm{M_\odot \, kpc^{-2}}} \right)^{\!N} \, ,
```
where N is the power-law index, and $A = \log_{10}(a)$, where $a$ is $\Sigma_\mathrm{SFR}$ at the fiducial gas surface density of $1.0 \, \mathrm{M_\odot \, kpc^{-2}}$.

# Arguments

  - `Σsfr ::Vector{<:MassFlowDensity}`: Values of the SFR surface density.
  - `log_output::Bool=true`: If the output will the $\log_{10} \, \Sigma_\text{H2}$. If `log_output` = true, the implied units are $1.0 \, \mathrm{M_\odot \, kpc^{-2}}$

# Returns

  - The molecular gas surface density.

# References

R. Cologni et al. (2026). *Resolved molecular gas and star formation in massive unquenched spirals*. Astronomy and Astrophysics, **709**, A148. [doi:10.1051/0004-6361/202558486](https://doi.org/10.1051/0004-6361/202558486)
"""
function invKSLawCologni2026(Σsfr::Vector{<:MassFlowDensity}; log_output::Bool=true)::Vector{<:Number}

    log10Σsfr = @. log10(ustrip(u"Msun * yr^-1 * kpc^-2", Σsfr))

    log10ΣH2 = @. (log10Σsfr - A_COLOGNI2026) / N_COLOGNI2026

    if log_output
        return log10ΣH2
    else
        return @. exp10(log10ΣH2) * 1.0u"Msun * kpc^-2"
    end

end

@doc raw"""
    KSLawLin2019(ΣH2::Vector{<:SurfaceDensity}; <keyword arguments>)::Vector{<:Number}

Evaluate the Kennicutt-Schmidt law for the molecular gas, taken from Lin et al. (2019).

From Lin et al. (2019) (Section 3.1, Table 1), we have

```math
\Sigma_\mathrm{SFR} = a \left( \frac{\Sigma_\mathrm{H_2}}{1.0 \, \mathrm{M_\odot \, kpc^{-2}}} \right)^{\!N} \, ,
```
where N is the power-law index, and $A = \log_{10}(a)$, where $a$ is $\Sigma_\mathrm{SFR}$ at the fiducial gas surface density of $1.0 \, \mathrm{M_\odot \, kpc^{-2}}$.

# Arguments

  - `ΣH2::Vector{<:SurfaceDensity}`: Values of the molecular gas surface density.
  - `log_output::Bool=true`: If the output will the $\log_{10} \, \Sigma_\text{SFR}$, or $\Sigma_\text{SFR}$. If `log_output` = true, the implied units are $\mathrm{M_\odot \, yr^{-1} \, kpc^{-2}}$

# Returns

  - The SFR surface density.

# References

L. Lin et al. (2019). *The ALMaQUEST Survey: The Molecular Gas Main Sequence and the Origin of the Star-forming Main Sequence*. The Astrophysical Journal Letters, **884(2)**, L33. [doi:10.3847/2041-8213/ab4815](https://doi.org/10.3847/2041-8213/ab4815)
"""
function KSLawLin2019(ΣH2::Vector{<:SurfaceDensity}; log_output::Bool=true)::Vector{<:Number}

    log10ΣH2 = @. log10(uconvert(Unitful.NoUnits, ΣH2 / 1.0u"Msun * kpc^-2"))

    log10Σsfr = @. A_LIN2019 + log10ΣH2 * N_LIN2019

    if log_output
        return log10Σsfr
    else
        return @. exp10(log10Σsfr) * u"Msun * yr^-1 * kpc^-2"
    end

end

@doc raw"""
    invKSLawLin2019(Σsfr ::Vector{<:MassFlowDensity}; <keyword arguments>)::Vector{<:Number}

Evaluate the inverse Kennicutt-Schmidt law for the molecular gas, taken from Lin et al. (2019).

From Lin et al. (2019) (Section 3.1, Table 1), we have

```math
\Sigma_\mathrm{SFR} = a \left( \frac{\Sigma_\mathrm{H_2}}{1.0 \, \mathrm{M_\odot \, kpc^{-2}}} \right)^{\!N} \, ,
```
where N is the power-law index, and $A = \log_{10}(a)$, where $a$ is $\Sigma_\mathrm{SFR}$ at the fiducial gas surface density of $1.0 \, \mathrm{M_\odot \, kpc^{-2}}$.

# Arguments

  - `Σsfr ::Vector{<:MassFlowDensity}`: Values of the SFR surface density.
  - `log_output::Bool=true`: If the output will the $\log_{10} \, \Sigma_\text{H2}$. If `log_output` = true, the implied units are $1.0 \, \mathrm{M_\odot \, kpc^{-2}}$

# Returns

  - The molecular gas surface density.

# References

L. Lin et al. (2019). *The ALMaQUEST Survey: The Molecular Gas Main Sequence and the Origin of the Star-forming Main Sequence*. The Astrophysical Journal Letters, **884(2)**, L33. [doi:10.3847/2041-8213/ab4815](https://doi.org/10.3847/2041-8213/ab4815)
"""
function invKSLawLin2019(Σsfr::Vector{<:MassFlowDensity}; log_output::Bool=true)::Vector{<:Number}

    log10Σsfr = @. log10(ustrip(u"Msun * yr^-1 * kpc^-2", Σsfr))

    log10ΣH2 = @. (log10Σsfr - A_LIN2019) / N_LIN2019

    if log_output
        return log10ΣH2
    else
        return @. exp10(log10ΣH2) * 1.0u"Msun * kpc^-2"
    end

end

@doc raw"""
    KSLawQuerejeta2021(ΣH2::Vector{<:SurfaceDensity}; <keyword arguments>)::Vector{<:Number}

Evaluate the Kennicutt-Schmidt law for the molecular gas, taken from Querejeta et al. (2021).

From Querejeta et al. (2021) (Section 4.3, Table 4, row All, mean of the columns), we have

```math
\Sigma_\mathrm{SFR} = a \left( \frac{\Sigma_\mathrm{H_2}}{1.0 \, \mathrm{M_\odot \, pc^{-2}}} \right)^{\!N} \, ,
```
where N is the power-law index, and $A = \log_{10}(a)$, where $a$ is $\Sigma_\mathrm{SFR}$ at the fiducial gas surface density of $1.0 \, \mathrm{M_\odot \, pc^{-2}}$.

# Arguments

  - `ΣH2::Vector{<:SurfaceDensity}`: Values of the molecular gas surface density.
  - `log_output::Bool=true`: If the output will the $\log_{10} \, \Sigma_\text{SFR}$, or $\Sigma_\text{SFR}$. If `log_output` = true, the implied units are $\mathrm{M_\odot \, yr^{-1} \, kpc^{-2}}$

# Returns

  - The SFR surface density.

# References

M. Querejeta et al. (2021). *Stellar structures, molecular gas, and star formation across the PHANGS sample of nearby galaxies*. Astronomy and Astrophysics, **656**, A133. [doi:10.1051/0004-6361/202140695](https://doi.org/10.1051/0004-6361/202140695)
"""
function KSLawQuerejeta2021(ΣH2::Vector{<:SurfaceDensity}; log_output::Bool=true)::Vector{<:Number}

    log10ΣH2 = @. log10(uconvert(Unitful.NoUnits, ΣH2 / 1.0u"Msun * pc^-2"))

    log10Σsfr = @. A_QUEREJETA2021 + log10ΣH2 * N_QUEREJETA2021

    if log_output
        return log10Σsfr
    else
        return @. exp10(log10Σsfr) * u"Msun * yr^-1 * kpc^-2"
    end

end

@doc raw"""
    invKSLawQuerejeta2021(Σsfr ::Vector{<:MassFlowDensity}; <keyword arguments>)::Vector{<:Number}

Evaluate the inverse Kennicutt-Schmidt law for the molecular gas, taken from Querejeta et al. (2021).

From Querejeta et al. (2021) (Section 4.3, Table 4, row All, mean of the columns), we have

```math
\Sigma_\mathrm{SFR} = a \left( \frac{\Sigma_\mathrm{H_2}}{1.0 \, \mathrm{M_\odot \, pc^{-2}}} \right)^{\!N} \, ,
```
where N is the power-law index, and $A = \log_{10}(a)$, where $a$ is $\Sigma_\mathrm{SFR}$ at the fiducial gas surface density of $1.0 \, \mathrm{M_\odot \, pc^{-2}}$.

# Arguments

  - `Σsfr ::Vector{<:MassFlowDensity}`: Values of the SFR surface density.
  - `log_output::Bool=true`: If the output will the $\log_{10} \, \Sigma_\text{H2}$. If `log_output` = true, the implied units are $1.0 \, \mathrm{M_\odot \, pc^{-2}}$

# Returns

  - The molecular gas surface density.

# References

M. Querejeta et al. (2021). *Stellar structures, molecular gas, and star formation across the PHANGS sample of nearby galaxies*. Astronomy and Astrophysics, **656**, A133. [doi:10.1051/0004-6361/202140695](https://doi.org/10.1051/0004-6361/202140695)
"""
function invKSLawQuerejeta2021(Σsfr::Vector{<:MassFlowDensity}; log_output::Bool=true)::Vector{<:Number}

    log10Σsfr = @. log10(ustrip(u"Msun * yr^-1 * kpc^-2", Σsfr))

    log10ΣH2 = @. (log10Σsfr - A_QUEREJETA2021) / N_QUEREJETA2021

    if log_output
        return log10ΣH2
    else
        return @. exp10(log10ΣH2) * 1.0u"Msun * pc^-2"
    end

end

"""
    hvcatImages(
        blocks_per_row::Int,
        paths::Vector{String};
        <keyword arguments>
    )::Nothing

Join several images vertically and horizontally.

The images in `paths` will fill the rows and columns starting at the top left, going from left to right and from top to bottom (row-major order).

# Arguments

  - `blocks_per_row::Int`: Number of columns.
  - `paths::Vector{String}`: Paths to the images.
  - `output_path::String="./joined_images.png"`: Path to the output image.
"""
function hvcatImages(
    blocks_per_row::Int,
    paths::Vector{String};
    output_path::String="./joined_images.png",
)::Nothing

    isempty(paths) && throw(ArgumentError("hvcatImages: `paths` is empty"))

    new_image = hvcat(blocks_per_row, [load(path) for path in paths]...)

    save(output_path, new_image)

    return nothing

end

"""
    rangeCut!(
        raw_values::Vector{<:Number},
        range::Tuple{<:Number,<:Number};
        <keyword arguments>
    )::Bool

Delete every element in `raw_values` that is outside the given `range`.

# Arguments

  - `raw_values::Vector{<:Number}`: Dataset that will be pruned.
  - `range::Tuple{<:Number,<:Number}`: The range in question.
  - `keep_edges::Bool=true`: If the edges of the range will be kept.
  - `min_left::Int=1`: Minimum number of values that need to be left after pruning to proceed with the transformation.

# Returns

  - If a transformation was performed.
"""
function rangeCut!(
    raw_values::Vector{<:Number},
    range::Tuple{<:Number,<:Number};
    keep_edges::Bool=true,
    min_left::Int=1,
)::Bool

    # Shortcut computation for special cases
    !(isempty(raw_values) || all(isinf.(range))) || return false

    if keep_edges

        # Check that after the transformation at least `min_left` elements will be left
        count(x -> range[1] <= x <= range[2], raw_values) >= min_left || return false

        # Delete elements outside of the provided range
        filter!(x -> range[1] <= x <= range[2], raw_values)

    else

        # Check that after the transformation at least `min_left` elements will be left
        count(x -> range[1] < x < range[2], raw_values) >= min_left || return false

        # Delete elements outside of the provided range
        filter!(x -> range[1] < x < range[2], raw_values)

    end

    return true

end

"""
    rangeCut!(
        m_data::Vector{<:Number},
        s_data::Vector,
        range::Tuple{<:Number,<:Number};
        <keyword arguments>
    )::Bool

Delete every element in `m_data` that is outside the given `range`.

Every corresponding element in `s_data` (i.e. with the same index) will be deleted too.

# Arguments

  - `m_data::Vector{<:Number}`: Master dataset that will be pruned.
  - `s_data::Vector`: Slave dataset that will be pruned according to which values of `m_data` are outside `range`.
  - `range::Tuple{<:Number,<:Number}`: The range in question.
  - `keep_edges::Bool=true`: If the edges of the range will be kept.
  - `min_left::Int=1`: Minimum number of values that need to be left in the master dataset after pruning to proceed with the transformation.

# Returns

  - If a transformation was performed.
"""
function rangeCut!(
    m_data::Vector{<:Number},
    s_data::Vector,
    range::Tuple{<:Number,<:Number};
    keep_edges::Bool=true,
    min_left::Int=1,
)::Bool

    # Shortcut computation for special cases
    !(isempty(m_data) || all(isinf.(range))) || return false

    (
        length(s_data) >= length(m_data) ||
        throw(ArgumentError("rangeCut!: `s_data` must have at least as many elements as `m_data`, \
        but I got length(`s_data`) = $(length(s_data)) < length(`m_data`) = $(length(m_data))"))
    )

    if keep_edges

        # Find the elements outside of the provided range
        idxs = map(x -> x < range[1] || x > range[2], m_data)

        # Check that after the transformation at least `min_left` elements will be left
        count(.!idxs) >= min_left || return false

        # Delete elements outside of the provided range
        deleteat!(m_data, idxs)
        deleteat!(s_data, idxs)

    else

        # Find the elements outside of the provided range
        idxs = map(x -> x <= range[1] || x >= range[2], m_data)

        # Check that after the transformation at least `min_left` elements will be left
        count(.!idxs) >= min_left || return false

        # Delete elements outside of the provided range
        deleteat!(m_data, idxs)
        deleteat!(s_data, idxs)

    end

    return true

end

"""
    sanitizeData!(
        raw_values::Vector{<:Number};
        <keyword arguments>
    )::NTuple{2,Bool}

Do the following transformations over `raw_values`, in order:

  - Trim it to fit within the domain of the function `func_domain`.
  - Trim it to fit within `range`.
  - Scale it down by a factor of 10^`exp_factor`.

By default, no transformation is done.

# Arguments

  - `raw_values::Vector{<:Number}`: Dataset to be sanitized.
  - `func_domain::Function=identity`: `raw_values` will be trimmed to fit within the domain of the function `func_domain`. The options are the scaling functions accepted by [Makie](https://docs.makie.org/stable/): log10, log2, log, sqrt, Makie.logit, Makie.Symlog10, Makie.pseudolog10, and identity.
  - `range::Tuple{<:Number,<:Number}=(-Inf, Inf)`: Every element in `raw_values` that falls outside of `range` will be deleted.
  - `keep_edges::Bool=true`: If the edges of `range` will be kept.
  - `min_left::Int=1`: Minimum number of values that need to be left after each transformation to proceed with it.
  - `exp_factor::Int=0`: Every element in `raw_values` will be divided by 10^`exp_factor`.

# Returns

  - A tuple with two flags:

      + If `raw_values` was mutated to fit within the domain of `func_domain`.
      + If `raw_values` was mutated to fit within `range`.
"""
function sanitizeData!(
    raw_values::Vector{<:Number};
    func_domain::Function=identity,
    range::Tuple{<:Number,<:Number}=(-Inf, Inf),
    keep_edges::Bool=true,
    min_left::Int=1,
    exp_factor::Int=0,
)::NTuple{2,Bool}

    !isempty(raw_values) || return false, false

    d_unit = unit(first(raw_values))

    # Trim `raw_values` to fit within the domain of `func_domain`
    if func_domain ∈ [identity, Makie.pseudolog10, Makie.Symlog10]

        domain_flag = false

    elseif func_domain == sqrt

        domain_flag = rangeCut!(raw_values, (0.0, Inf) .* d_unit; keep_edges=true, min_left)

    elseif func_domain == Makie.logit

        domain_flag = rangeCut!(raw_values, (0.0, 1.0) .* d_unit; keep_edges=false, min_left)

    elseif func_domain ∈ [log, log2, log10]

        domain_flag = rangeCut!(raw_values, (0.0, Inf) .* d_unit; keep_edges=false, min_left)

    else

        throw(ArgumentError("sanitizeData!: The function $(func_domain) is not supported. See \
        the list of supported scaling functions in the [Makie](https://docs.makie.org/stable/) \
        documentation"))

    end

    # Trim `raw_values` to fit within `range`
    range_flag = rangeCut!(raw_values, range; keep_edges, min_left)

    (
        !(isa(raw_values, Vector{<:Integer}) && !iszero(exp_factor) && LOGGING[]) ||
        @warn("sanitizeData!: Elements of `raw_values` are of type `Integer`, this may result \
        in errors or unwanted truncation when using `exp_factor` != 0")
    )

    # Scale `raw_values` down by a factor of 10^`exp_factor`
    iszero(exp_factor) || (raw_values ./= exp10(exp_factor))

    return domain_flag, range_flag

end

"""
    sanitizeData!(
        x_data::Vector{<:Number},
        y_data::Vector{<:Number};
        <keyword arguments>
    )::NTuple{4,Bool}

Do the following transformations over `x_data` and `y_data`, in order:

  - Trim them to fit within the domain of the functions `func_domain[1]` and `func_domain[2]`, respectively.
  - Trim them to fit within `range[1]` and `range[2]`, respectively.
  - Scale them down by a factor 10^`exp_factor[1]` and 10^`exp_factor[2]`, respectively.

By default, no transformation is done.

!!! note

    The datasets must have the same length, and any operation that deletes an element, will delete the corresponding element (i.e. with the same index) in the other dataset, so that the datasets will remain of equal length.

# Arguments

  - `x_data::Vector{<:Number}`: First dataset to be sanitized.
  - `y_data::Vector{<:Number}`: Second dataset to be sanitized.
  - `func_domain::NTuple{2,Function}=(identity, identity)`: `x_data` will be trimmed to fit within the domain of the function `func_domain[1]`, and `y_data` will be trimmed to fit within the domain of the function `func_domain[2]`. The options are the scaling functions accepted by [Makie](https://docs.makie.org/stable/): log10, log2, log, sqrt, Makie.logit, Makie.Symlog10, Makie.pseudolog10, and identity.
  - `range::Tuple{Tuple{<:Number,<:Number},Tuple{<:Number,<:Number}}=((-Inf, Inf), (-Inf, Inf))`: Every element in `x_data` that falls outside of `range[1]` will be deleted, and every element in `y_data` that falls outside of `range[2]` will be deleted.
  - `keep_edges::NTuple{2,Bool}=(true, true)`: If the edges of each corresponding `range` will be kept.
  - `min_left::Int=1`: Minimum number of values that need to be left in each dataset after any of the transformations to proceed with them.
  - `exp_factor::NTuple{2,Int}=(0, 0)`: Every element in `x_data` will be divided by 10^`exp_factor[1]`, and every element in `y_data` will be divided by 10^`exp_factor[2]`.

# Returns

  - A tuple with four flags:

      + If `x_data` was mutated to fit within the domain of `func_domain[1]`.
      + If `y_data` was mutated to fit within the domain of `func_domain[2]`.
      + If `x_data` was mutated to fit within `range[1]`.
      + If `y_data` was mutated to fit within `range[2]`.
"""
function sanitizeData!(
    x_data::Vector{<:Number},
    y_data::Vector{<:Number};
    func_domain::NTuple{2,Function}=(identity, identity),
    range::Tuple{Tuple{<:Number,<:Number},Tuple{<:Number,<:Number}}=((-Inf, Inf), (-Inf, Inf)),
    keep_edges::NTuple{2,Bool}=(true, true),
    min_left::Int=1,
    exp_factor::NTuple{2,Int}=(0, 0),
)::NTuple{4,Bool}

    (
        length(x_data) == length(y_data) ||
        throw(ArgumentError("sanitizeData!: `x_data` and `y_data` must have the same length, \
        but I got length(x_data) = $(length(x_data)) != length(y_data) = $(length(y_data))"))
    )

    x_unit = isempty(x_data) ? Unitful.NoUnits : unit(first(x_data))

    # Trim the data to fit within the domain of `func_domain[1]`
    if func_domain[1] ∈ [identity, Makie.pseudolog10, Makie.Symlog10]

        x_domain_flag = false

    elseif func_domain[1] == sqrt

        x_domain_flag = rangeCut!(x_data, y_data, (0.0, Inf) .* x_unit; keep_edges=true, min_left)

    elseif func_domain[1] == Makie.logit

        x_domain_flag = rangeCut!(x_data, y_data, (0.0, 1.0) .* x_unit; keep_edges=false, min_left)

    elseif func_domain[1] ∈ [log, log2, log10]

        x_domain_flag = rangeCut!(x_data, y_data, (0.0, Inf) .* x_unit; keep_edges=false, min_left)

    else

        throw(ArgumentError("sanitizeData!: The function $(func_domain[1]) is not supported. See \
        the list of supported scaling functions in the [Makie](https://docs.makie.org/stable/) \
        documentation"))

    end

    y_unit = isempty(y_data) ? Unitful.NoUnits : unit(first(y_data))

    # Trim the data to fit within the domain of `func_domain[2]`
    if func_domain[2] ∈ [identity, Makie.pseudolog10, Makie.Symlog10]

        y_domain_flag = false

    elseif func_domain[2] == sqrt

        y_domain_flag = rangeCut!(y_data, x_data, (0.0, Inf) .* y_unit; keep_edges=true, min_left)

    elseif func_domain[2] == Makie.logit

        y_domain_flag = rangeCut!(y_data, x_data, (0.0, 1.0) .* y_unit; keep_edges=false, min_left)

    elseif func_domain[2] ∈ [log, log2, log10]

        y_domain_flag = rangeCut!(y_data, x_data, (0.0, Inf) .* y_unit; keep_edges=false, min_left)

    else

        throw(ArgumentError("sanitizeData!: The function $(func_domain[2]) is not supported. See \
        the list of supported scaling functions in the [Makie](https://docs.makie.org/stable/) \
        documentation"))

    end

    # Trim data to fit within `range[1]`
    x_range_flag = rangeCut!(x_data, y_data, range[1]; keep_edges=keep_edges[1], min_left)

    # Trim data to fit within `range[2]`
    y_range_flag = rangeCut!(y_data, x_data, range[2]; keep_edges=keep_edges[2], min_left)

    (
        !(isa(x_data, Vector{<:Integer}) && !iszero(exp_factor[1]) && LOGGING[]) ||
        @warn("sanitizeData!: Elements of `x_data` are of type Integer, this may result \
        in errors or unwanted truncation when using `exp_factor[1]` != 0")
    )

    (
        !(isa(y_data, Vector{<:Integer}) && !iszero(exp_factor[2]) && LOGGING[]) ||
        @warn("sanitizeData!: Elements of `y_data` are of type Integer, this may result \
        in errors or unwanted truncation when using `exp_factor[2]` != 0")
    )

    # Scale the data down by the factors `exp_factor`
    iszero(exp_factor[1]) || (x_data ./= exp10(exp_factor[1]))
    iszero(exp_factor[2]) || (y_data ./= exp10(exp_factor[2]))

    return x_domain_flag, y_domain_flag, x_range_flag, y_range_flag

end

"""
    smoothWindow(
        x_data::Vector{<:Number},
        y_data::Vector{<:Number},
        n_bins::Int;
        <keyword arguments>
    )::NTuple{2,Vector{<:Number}}

Separate the values of `x_data` in `n_bins` bins and compute the mean value of `x_data` and `y_data` within each one.

# Arguments

  - `x_data::Vector{<:Number}`: x-axis data.
  - `y_data::Vector{<:Number}`: y-axis data.
  - `n_bins::Int`: Number of bins.
  - `scaling::Function=identity`: Scaling function for the x axis. The options are the scaling functions accepted by [Makie](https://docs.makie.org/stable/): log10, log2, log, sqrt, Makie.logit, Makie.Symlog10, Makie.pseudolog10, and identity. All the values of `x_data` must be in the domain of `scaling`.

# Returns

  - A tuple with two vectors, containing the smoothed-out x and y values.
"""
function smoothWindow(
    x_data::Vector{<:Number},
    y_data::Vector{<:Number},
    n_bins::Int;
    scaling::Function=identity,
)::NTuple{2,Vector{<:Number}}

    # Check that the input vectors have the same length
    (
        length(x_data) == length(y_data) ||
        throw(DimensionMismatch("smoothWindow: `x_data` and `y_data` must have the same length, \
        but I got length(`x_data`) = $(length(x_data)) != length(`y_data`) = $(length(y_data))"))
    )

    positions = scaling.(ustrip(x_data))
    grid = LinearGrid(extrema(positions)..., n_bins)

    smooth_x_data = histogram1D(positions, x_data, grid; total=false)
    smooth_y_data = histogram1D(positions, y_data, grid; total=false)

    # Remove empty bins
    return filter!(!isnan, smooth_x_data), filter!(!isnan, smooth_y_data)

end

"""
    formatSeconds(sec::Float64)::String

Format given number of seconds as "DDd-HHh:MM':SS''", where D is days, H is hours, M is minutes, and S is seconds.

# Arguments

  - `sec::Float64`: Number of seconds.

# Returns

  - String with the formatted time.
"""
function formatSeconds(sec::Float64)::String

    days    = floor(sec / (24 * 3600))
    rem1    = sec - days * (24 * 3600)
    hours   = floor(rem1 / 3600)
    rem2    = rem1 - hours * 3600
    minutes = floor(rem2 / 60)
    seconds = round(rem2 - minutes * 60)

    return string(
        Int(days),
        "d ",
        lpad(Int(hours), 2, "0"),
        "hs ",
        lpad(Int(minutes), 2, "0"),
        "' ",
        lpad(Int(seconds), 2, "0"),
        "''",
    )

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
