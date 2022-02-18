####################################################################################################
# Utilities
####################################################################################################
"""
Function that always returns nothing, for any type and number of arguments.
"""
trivial(x...;y...) = nothing

"""
Read the numeric data of a CSV file using DelimitedFiles. It skips the first line.
"""
read_csv(path::String) = readdlm(path, ',', Float64, skipstart = 1)

"""
Extension of Base.isempty to check for empty LaTeXStrings.
"""
Base.isempty(l_str::LaTeXString) = l_str == L""

"""
Extract the limits of the axes from a [Makie.jl](http://makie.juliaplots.org/stable/index.html) plot, axis, or figure. In the case of a figure, it 
will take the limits from the current axis object.
"""
xlimits(plot::Makie.FigureAxisPlot) = plot.axis.elements[:xaxis].attributes[:limits][]
xlimits(axis::Makie.Axis) = axis.elements[:xaxis].attributes[:limits][]
xlimits(fig::Makie.Figure) = fig.current_axis.x.elements[:xaxis].attributes[:limits][]
ylimits(plot::Makie.FigureAxisPlot) = plot.axis.elements[:yaxis].attributes[:limits][]
ylimits(axis::Makie.Axis) = axis.elements[:yaxis].attributes[:limits][]
ylimits(fig::Makie.Figure) = fig.current_axis.x.elements[:yaxis].attributes[:limits][]

"""
Extract the scales of the axes from a [Makie.jl](http://makie.juliaplots.org/stable/index.html) plot, axis, or figure. In the case of a figure, it 
will take the scales from the current axis object.
"""
xscale(plot::Makie.FigureAxisPlot) = plot.axis.elements[:xaxis].attributes[:scale][]
xscale(axis::Makie.Axis) = axis.elements[:xaxis].attributes[:scale][]
xscale(fig::Makie.Figure) = fig.current_axis.x.elements[:xaxis].attributes[:scale][]
yscale(plot::Makie.FigureAxisPlot) = plot.axis.elements[:yaxis].attributes[:scale][]
yscale(axis::Makie.Axis) = axis.elements[:yaxis].attributes[:scale][]
yscale(fig::Makie.Figure) = fig.current_axis.x.elements[:yaxis].attributes[:scale][]

"""
Extract the data points from a [Makie.jl](http://makie.juliaplots.org/stable/index.html) plot, axis, or figure, as a Vector{Point2D}. For the three 
cases the last series is taken. And in the case of a figure, it will only take the data from the 
current axis object. It may give non-sensical results or an error for anything other than scatter 
plots. 
"""
pointData(plot::Makie.FigureAxisPlot) = plot.axis.scene.plots[end][1][]
pointData(axis::Makie.Axis) = axis.scene.plots[end][1][]
pointData(fig::Makie.Figure) = fig.current_axis.x.scene.plots[end][1][]

"""
    absCoor(
        pl::Union{Makie.FigureAxisPlot,Makie.Axis,Makie.Figure},
        rx::Float64,
        ry::Float64,
    )::NTuple{2,Float64}
    
Give the absolute x and y coordinates from the relative ones.

# Arguments 
- `pl::Union{Makie.FigureAxisPlot, Makie.Axis, Makie.Figure}`: Plot, axis, or figure for which the 
  absolute coordinates will be calculated. In the case of a figure, it will use the limits from the 
  current axis object.
- `rx::Float64`: relative x coordinate, `rx` ∈ [0, 1].
- `ry::Float64`: relative y coordinate, `ry` ∈ [0, 1].

# Returns
- A tuple with the absolute coordinates: (x, y).

# Examples
```julia-repl
julia> GadgetInspector.absCoor(lines(rand(100)), 0.5, 0.5)
(50.49999666213989, 0.5031078979372978)
```

"""
function absCoor(
    pl::Union{Makie.FigureAxisPlot,Makie.Axis,Makie.Figure},
    r_x::Float64,
    r_y::Float64,
)::NTuple{2,Float64}

    # Get the scaling functions
    x_scale = xscale(pl)
    y_scale = yscale(pl)

    # Get the axes' limits
    x_limits = x_scale.(xlimits(pl))
    y_limits = y_scale.(ylimits(pl))

    # Aboslute coordinates
    a_x = Makie.inverse_transform(x_scale).(x_limits[1] + r_x * (x_limits[2] - x_limits[1]))
    a_y = Makie.inverse_transform(y_scale).(y_limits[1] + r_y * (y_limits[2] - y_limits[1]))

    return a_x, a_y

end

"""
    cleanPlot!(figure::Makie.Figure)::Nothing 

Delete all the legends of a figure and empty all its axes.

# Arguments
- `figure::Makie.Figure`: Figure to be cleaned.

"""
function cleanPlot!(figure::Makie.Figure)::Nothing

    # Number of elements (axes and legends) to be cleaned
    n_elements = length(figure.content)
    i = 1

    for _ = 1:n_elements
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
    cleanPlot!(legend::Makie.Legend)::Bool

Delete a legend.

# Arguments
- `legend::Makie.Legend`: Legend to be deleted.

# Returns
- Flag to indicate that a legend has been deleted.

"""
function cleanPlot!(legend::Makie.Legend)::Bool

    delete!(legend)

    return false

end

"Default function to end the recursion if an unknown type is encountered."
cleanPlot!(default) = error("I cannot clean element of type $(typeof(default))")

"""
    getUnitLabel(
        factor::Int64, 
        unit::Unitful.Units; 
        <keyword arguments>
    )::Union{String,LaTeXString}

Construct the unit part of an axis label.

# Arguments
- `factor::Int64`: Exponential factor to scale the units. If different from 0, a term of the form 
  "10^`factor`" will be added to the label.
- `unit::Unitful.Units`: Unit of the axis.
- `latex::Bool = true`: If the output will be a LaTeXString (default) or a plain String.

# Returns
- The LaTeXString or String: "10^`factor` `unit`". The `factor` term only appears if `factor` != 0, 
  the unit term only appears if `unit` != `Unitful.NoUnits`.  

"""
function getUnitLabel(
    factor::Int64,
    unit::Unitful.Units;
    latex::Bool = true,
)::Union{String,LaTeXString}

    if latex

        str_unit = replace(
            string(unit),
            "M⊙" => raw"M_{\!\odot}",
            r"\^(-?[0-9]+)" => s"^{\1}",
            " " => raw"\,",
        )

        if factor == 0
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

        if factor == 0
            if isempty(str_unit)
                out_str = "dimensionless"
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
        label::Union{String,LaTeXString}, 
        factor::Int64, 
        unit::Unitful.Units; 
        <keyword arguments>
    )::Union{String,LaTeXString}

Construct an axis label.

# Arguments
- `label::Union{String,LaTeXString}: Text part of the label.
- `factor::Int64`: Exponential factor to scale the units. If different from 0, a term of the form 
  "10^`factor`" will be added to the label.
- `unit::Unitful.Units`: Unit of the axis.
- `latex::Bool = true`: If the output will be a LaTeXString (default) or a plain String.

# Returns
- The LaTeXString or String: "`label` / 10^`factor` `unit`". If `label` is "", an empty string is 
  returned. The `factor` term only appears if `factor` != 0, the unit term only appears if 
  `unit` != `Unitful.NoUnits`, and the division only appears if there are a factor and/or a unit term.  

"""
function getLabel(
    label::Union{String,LaTeXString},
    factor::Int64,
    unit::Unitful.Units;
    latex::Bool = true,
)::Union{String,LaTeXString}

    if isempty(label)
        return ""
    end

    unit_label = getUnitLabel(factor, unit; latex)

    if latex

        if isempty(unit_label)
            out_str = LaTeXString(label)
        else
            out_str = LaTeXString(label * L"\; / \;" * unit_label)
        end

    else

        if unit_label == "dimensionless"
            out_str = label
        else
            out_str = label * " / $unit_label"
        end

    end

    return out_str

end

"""
    scaledBins(
        values::AbstractVector{<:RealOrQty}, 
        bins::Int64,
        scale::Function; 
        <keyword arguments>
    )::Vector{Float64}

Compute the bin edges, for a given scaling function and set of values. 

# Arguments
- `values::AbstractVector{<:RealOrQty}`: Values to be binned.
- `bins::Int64`: Number of bins to be calculated.
- `scale::Function`: Function to scale the `values`. The possibilities are the scaling functions 
  accepted by [Makie.jl](http://makie.juliaplots.org/stable/index.html): log10, log2, log, sqrt, Makie.logit, Makie.Symlog10, Makie.pseudolog10 and 
  identity.
- `limits::Tuple{<:Union{Real,Nothing},<:Union{Real,Nothing}} = (nothing, nothing)`: Set it to a 
  value different than `nothing` if you want to fix the limits of the axis.
- `warnings::Bool = true`: If a warning will be given when the data is not as expected, but the 
  function can still run using sane defaults.

# Returns
- A sorted list of bin edges.

"""
function scaledBins(
    values::AbstractVector{<:RealOrQty},
    bins::Int64,
    scale::Function;
    limits::Tuple{<:Union{Real,Nothing},<:Union{Real,Nothing}} = (nothing, nothing),
    warnings::Bool = true,
)::Vector{Float64}

    min_val = scale(ustrip(minimum(values)))
    max_val = scale(ustrip(maximum(values)))

    min = limits[1] === nothing ? min_val : scale(limits[1])
    max = limits[2] === nothing ? max_val : scale(limits[2])

    if max < min
        !warnings || @warn("Incompatibility between `limits` and `values`. Ignoring limits...")
        min = min_val
        max = max_val
    end

    range = max - min

    # For small ranges, increase it simetricaly by 0.2 * abs(max)
    if range <= 1e-4 * abs(max)
        min -= 0.1 * abs(max)
        max += 0.1 * abs(max)
        range = max - min
    end

    width = range / bins

    return [Makie.inverse_transform(scale)(min + width * (i - 1)) for i = 1:(bins+1)]

end

"""
    smoothWindow(
        x_data::AbstractVector{<:RealOrQty},
        y_data::AbstractVector{<:RealOrQty},
        bins::Int64;
        <keyword arguments>
    )::NTuple{2,Vector{<:RealOrQty}}

Separate the `x_data` values in `bins` contiguous windows, and replaces every x and y 
value within the window with the corresponding mean, to smooth out the data. 

`x_data` should be within the domain of `f_scale`.

# Arguments
- `x_data::AbstractVector{<:RealOrQty}`: x-axis data.
- `y_data::AbstractVector{<:RealOrQty}`: y-axis data.
- `bins::Int64`: Number of bins to be used in the smoothing.
- `scale::Function = identity`: Function to scale the x-axis. The bins will be computed accordingly. 
  The options are the scaling functions accepted by [Makie.jl](http://makie.juliaplots.org/stable/index.html): log10, log2, log, sqrt, Makie.logit, 
  Makie.Symlog10, Makie.pseudolog10 and identity.

# Returns
- A tuple with two arrays containing the smoothed out x and y data.

"""
function smoothWindow(
    x_data::AbstractVector{<:RealOrQty},
    y_data::AbstractVector{<:RealOrQty},
    bins::Int64;
    scale::Function = identity,
)::NTuple{2,Vector{<:RealOrQty}}

    # Check that the input vectors have the same length
    (
        length(x_data) == length(y_data) ||
        throw(DimensionMismatch("The input vectors should have the same length."))
    )

    # Compute bins
    bin_limits = scaledBins(x_data, bins, scale)

    # Initialize output arrays
    smooth_x_data = Vector{typeof(x_data[1])}(undef, bins)
    smooth_y_data = Vector{typeof(y_data[1])}(undef, bins)

    @inbounds for i in eachindex(smooth_x_data, smooth_y_data)

        # Find the indices of `x_data` which fall within `i`-th window 
        idx = findall(x -> bin_limits[i] <= x < bin_limits[i+1], ustrip.(x_data))

        if isempty(idx)
            error("$bins bins is too big for the data, there are empty bins.")
        else
            # Store the mean values in output arrays
            smooth_x_data[i] = sum(x_data[idx]) / length(idx)
            smooth_y_data[i] = sum(y_data[idx]) / length(idx)
        end

    end

    return smooth_x_data, smooth_y_data

end

"""
    formatError(q_mean::RealOrQty, q_error::RealOrQty)::NTuple{2,<:RealOrQty}

Format nicely a set of mean and error values.

It follows the traditional rules for error presentation, the error has only one significant digit, 
unless such digit is a one, in which case, two significant digits are used.  
The mean will have a number of digits such as to match the last significant position of the error. 
An error equal to 0.0 will leave the mean unmodified, and an error bigger than the mean may set the 
mean equal to 0.0.

# Arguments 
- `q_mean::RealOrQty`: Mean value.
- `q_error::RealOrQty`: Error value. It must be positive.

# Returns
- A Tuple with the formatted mean and error values.

# Examples
```julia-repl
julia> formatError(69.42069, 0.038796)
(69.42, 0.04)

julia> formatError(69.42069, 0.018796)
(69.421, 0.019)

julia> formatError(69.42069, 0.0)
(69.42069, 0.0)

julia> formatError(69.42069, 93.4)
(70.0, 90.0)
```

"""
function formatError(q_mean::RealOrQty, q_error::RealOrQty)::NTuple{2,<:RealOrQty}

    # Positive error check
    q_error >= zero(q_error) || throw(DomainError("The error must be a positive number."))

    # Use only the magnitudes for the computations
    mean = ustrip(q_mean)
    error = ustrip(q_error)

    if error == 0.0

        round_mean = mean
        round_error = error

    else

        sigdigit_pos = abs(log10(error))

        if error < 1.0
            first_digit = trunc(error * 10.0^(floor(sigdigit_pos) + 1.0))
            first_digit == 1.0 ? extra = 1 : extra = 0
            digits = ceil(Int64, sigdigit_pos) + extra
            round_mean = round(mean; digits)
            round_error = round(error, sigdigits = 1 + extra)
        else
            first_digit = trunc(error / 10.0^(floor(sigdigit_pos)))
            first_digit == 1.0 ? extra = 2 : extra = 1
            sigdigits = ceil(Int64, log10(abs(mean))) - ceil(Int64, sigdigit_pos) + extra
            round_mean = round(mean; sigdigits)
            round_error = round(error, sigdigits = extra)
        end

    end

    return round_mean * unit(q_mean), round_error * unit(q_error)

end

"""
    passAll(file_path::String, type::String)::Vector{Int64}

Default filter function needed for the [read\\_blocks\\_over\\_all\\_files](https://ludwigboess.github.io/GadgetIO.jl/stable/api/#GadgetIO.read_blocks_over_all_files-Tuple{String,%20Array{String,%20N}%20where%20N}) function.

It does not filter out any particles, allowing the data acquisition functions to gather all the data.

# Arguments 
- `file_path::String`: Path to the snapshot file to be filtered.
- `type::Symbol`: Particle type. The possibilities are given by the constant [`ParticleType`](@ref) 
  in `src/constants.jl`.

# Returns
- A Vector with the indices of the allowed particles.

"""
function passAll(file_path::String, type::Symbol)::Vector{Int64}

    (
        type ∈ keys(ParticleType) ||
        throw(ArgumentError("Particle type '$type' not supported. The supported types are \
        $(keys(ParticleType))"))
    )

    header = read_header(file_path)

    return collect(1:header.npart[ParticleType[type]+1])

end

"""
    passCritRho(file_path::String; <keyword arguments>)::Vector{Int64}

Filters out particles with densities below `crit_ρ`.

Onme possible filter function needed for the [read\\_blocks\\_over\\_all\\_files](https://ludwigboess.github.io/GadgetIO.jl/stable/api/#GadgetIO.read_blocks_over_all_files-Tuple{String,%20Array{String,%20N}%20where%20N}) function.

# Arguments 
- `file_path::String`: Path to the snapshot file to be filtered.
- `type::Symbol = :gas`: Particle type. The possibilities are given by the constant [`ParticleType`](@ref) 
  in `src/constants.jl`.
- `crit_ρ::Unitful.Quantity = CRITICAL_DENSITY`: Critical density, with units.
- `sim_cosmo::Bool = false`: If the simulation is cosmological, 
  - `false` ⟶ Newtonian simulation (`ComovingIntegrationOn` = 0).
  - `true` ⟶ Cosmological simulation (`ComovingIntegrationOn` = 1).

# Returns
- A Vector with the indices of the allowed particles.

"""
function passCritRho(
    file_path::String;
    type::Symbol = :gas,
    crit_ρ::Unitful.Quantity = CRITICAL_DENSITY,
    sim_cosmo::Bool = false,
)::Vector{Int64}

    type ∈ keys(ParticleType) || throw(ArgumentError(
        "Particle type '$type' not supported. The supported types are $(keys(ParticleType))"
    ))

    internal_unit = internalUnits("RHO", read_header(file_path); sim_cosmo)
    densities = read_snap(file_path, "RHO", ParticleType[type]) .* internal_unit

    return findall(densities .> crit_ρ)

end

"""
    passPositiveQty(file_path::String; <keyword arguments>)::Vector{Int64}

Filters out particles with negative values for a given quantity.

One possible filter function needed for the [read\\_blocks\\_over\\_all\\_files](https://ludwigboess.github.io/GadgetIO.jl/stable/api/#GadgetIO.read_blocks_over_all_files-Tuple{String,%20Array{String,%20N}%20where%20N}) function.

# Arguments 
- `file_path::String`: Path to the snapshot file to be filtered.
- `qty::String = "FMOL"`: Quantity to be checked.

# Returns
- A Vector with the indices of the allowed particles.

"""
function passPositiveQty(file_path::String; qty::String = "FMOL")::Vector{Int64}

    qty ∈ keys(QUANTITIES) || throw(ArgumentError(
        "Quantity '$qty' not supported. The supported types are $(keys(QUANTITIES))"
    ))

    block_present(GadgetIO.select_file(file_path, 0), qty) || @error(
        "I couldn't find the data block for quantity '$qty' in the snapshot $file_path"
    )

    data = read_snap(file_path, qty, ParticleType[:gas])

    return findall(data .>= 0.0)

end

"""
    passMetallicity(
        file_path::String,
        z_min::Real, 
        z_max::Real; 
        <keyword arguments>
    )::Vector{Int64}

Filters out gas particles with a metallicity outside a given range.

One possible filter function needed for the [read\\_blocks\\_over\\_all\\_files](https://ludwigboess.github.io/GadgetIO.jl/stable/api/#GadgetIO.read_blocks_over_all_files-Tuple{String,%20Array{String,%20N}%20where%20N}) function.

# Arguments 
- `file_path::String`: Path to the snapshot file to be filtered.
- `z_min::Real`: Minimum value of metallicity (as the dimensionless fraction of metal mass) allowed.
- `z_max::Real`: Maximum value of metallicity (as the dimensionless fraction of metal mass) allowed.
- `sim_cosmo::Bool = false`: If the simulation is cosmological, 
  - `false` ⟶ Newtonian simulation (`ComovingIntegrationOn` = 0).
  - `true` ⟶ Cosmological simulation (`ComovingIntegrationOn` = 1).
- `warnings::Bool = true`: If a warning will be given when the data is not as expected, but the 
  function can still run using sane defaults.

# Returns
- A Vector with the indices of the allowed particles.

"""
function passMetallicity(
    file_path::String,
    z_min::Real,
    z_max::Real;
    sim_cosmo::Bool = false,
    warnings::Bool = true,
)::Vector{Int64}

    data = getSnapshotData(file_path, :gas, ["Z", "MASS"]; sim_cosmo, warnings)[:gas]

    z = computeMetallicity(data["Z"], data["MASS"], solar = true)

    return findall(z_min .<= z .<= z_max)

end

@doc raw"""
    energyIntegrand(header::GadgetIO.SnapshotHeader, a::Float64)::Float64

The integrand of the integral that converts the scale factor to physical time,

```math
\frac{1}{H\,\sqrt{\epsilon}} \, ,
``` 

where 

```math
\epsilon = \Omega_\lambda + \frac{1 - \Omega_\lambda - \Omega_0}{a^2} + \frac{\Omega_0}{a^3} \, , 
```
```math
H = H_0 \, a \, .
```

# Arguments 
- `header::GadgetIO.SnapshotHeader`: Header of the relevant snapshot file.
- `a::Float64`: Dimensionless scale factor.

# Returns
- The integrand evaluated at value `a`. The result is in Gyr. 

"""
function energyIntegrand(header::GadgetIO.SnapshotHeader, a::Float64)::Float64

    # The integrand goes to 0 in the limit a ⟶ 0
    a != 0 || return 0

    # Ω_K (curvature)
    omega_K = 1.0 - header.omega_0 - header.omega_l
    # Energy function
    E = header.omega_0 / (a * a * a) + omega_K / (a * a) + header.omega_l
    # Hubble constant in 1 / Gyr
    H = header.h0 * HUBBLE_CONST * a

    # Integrand of the time integral in Gyr
    return 1.0 / (H * sqrt(E))

end

"""
    maxLength(data::Matrix{<:RealOrQty})::RealOrQty

The maximum value of the norms of the columns in `data`.

# Arguments
- `data::Matrix{<:RealOrQty}`: The values to be analyzed.

# Returns
- The maximum value of the norms of the columns in `data`.

"""
@inline function maxLength(data::Matrix{<:RealOrQty})::RealOrQty

    return maximum([norm(col) for col in eachcol(data)])

end

"""
    compare(x::Dict, y::Dict; <keyword arguments>)::Bool

Determine if two dictionaries are equal (or approximately equal for numeric values).

Every pair of values is compared with the [`compare`](@ref) function.

# Arguments
- `x::Dict`: First dictionary to be compared.
- `y::Dict`: Second dictionary to be compared.
- `atol::Float64 = 1e-5`: Absolute tolerance (for numeric elements).
- `rtol::Float64 = 1e-5`: Relative tolerance (for numeric elements).

# Returns
- Return `true` if every pair of elements within the dictionaries pass the equality tests.

"""
function compare(x::Dict, y::Dict; atol::Float64 = 1e-5, rtol::Float64 = 1e-5)::Bool

    # If the dictionaries have different keys avoid recursion and return false
    keys(x) == keys(y) || return false

    return all(compare(x[key], y[key]; atol, rtol) for key in keys(x))

end

"""
    compare(x::Dict, y::Dict; <keyword arguments>)::Bool

Determine if two dataframes are equal (or approximately equal for numeric values).

Every pair of values is compared with the [`compare`](@ref) function.

# Arguments
- `x::DataFrame`: First dataframe to be compared.
- `y::DataFrame`: Second dataframe to be compared.
- `atol::Float64 = 1e-5`: Absolute tolerance (for numeric elements).
- `rtol::Float64 = 1e-5`: Relative tolerance (for numeric elements).

# Returns
- Return `true` if every pair of elements within the dataframes pass the equality tests.

"""
function compare(x::DataFrame, y::DataFrame; atol::Float64 = 1e-5, rtol::Float64 = 1e-5)::Bool

    # If the dataframes have different column names avoid recursion and return false
    names(x) == names(y) || return false

    return all(compare.(eachcol(x), eachcol(y); atol, rtol))

end

"""
    compare(
        x::Union{AbstractVector,Tuple}, 
        y::Union{AbstractVector,Tuple};
        <keyword arguments>
    )::Bool

Determine if two list-like elements are equal (or approximately equal for numeric values).

Every pair of values is compared with the [`compare`](@ref) function.
        
# Arguments
- `x::Union{AbstractVector,Tuple}`: First list to be compared.
- `y::Union{AbstractVector,Tuple}`: Second list to be compared.
- `atol::Float64 = 1e-5`: Absolute tolerance (for numeric elements).
- `rtol::Float64 = 1e-5`: Relative tolerance (for numeric elements).

# Returns
- Returns true if the lists are equal, false if not.

"""
function compare(
    x::Union{AbstractVector,Tuple},
    y::Union{AbstractVector,Tuple};
    atol::Float64 = 1e-5,
    rtol::Float64 = 1e-5,
)::Bool

    # If the list have different lengths avoid recursion and return false
    length(x) == length(y) || return false

    return all(compare.(x, y; atol, rtol))

end

"""
    compare(
        x::RealOrQty, 
        y::RealOrQty; 
        <keyword arguments>
    )::Bool

Determine if two numeric elements are approximately equal, using the [isequal](https://docs.julialang.org/en/v1/base/base/#Base.isequal) function for elements 
with NaN, and using [isapprox](https://docs.julialang.org/en/v1/base/math/#Base.isapprox) for everything else.

# Arguments
- `x::RealOrQty`: First element to be compared.
- `y::RealOrQty`: Second element to be compared.
- `atol::Float64 = 1e-5`: Absolute tolerance.
- `rtol::Float64 = 1e-5`: Relative tolerance.

# Returns
- Returns true when the values are approximately equal, and false otherwise.

"""
function compare(x::RealOrQty, y::RealOrQty; atol::Float64 = 1e-5, rtol::Float64 = 1e-5)::Bool

    if x === NaN || y === NaN
        return isequal(x, y)
    else
        return Unitful.isapprox(x, y; atol = atol * unit(x), rtol)
    end

end

"""
    compare(x, y; <keyword arguments>)::Bool

Determine if two elements are equal, using the [isequal](https://docs.julialang.org/en/v1/base/base/#Base.isequal) function.

It's the default method, used when all other `compare` methods cannot be dispatched.

# Arguments
- `x`: First element to be compared.
- `y`: The second element to be compared.
- `atol::Float64 = 1e-5`: Absolute tolerance (just for compatibility with the other methods).
- `rtol::Float64 = 1e-5`: Relative tolerance (just for compatibility with the other methods).

# Returns
- Returns `isequal(x, y)`.

"""
function compare(x, y; atol::Float64 = 1e-5, rtol::Float64 = 1e-5)::Bool

    atol;
    rtol;

    return isequal(x, y)

end

"""
    inRange(
        data::AbstractVector{<:RealOrQty},
        inferior::RealOrQty = -Inf,
        superior::RealOrQty = Inf;
        <keyword arguments>
    )::Bool

If every element in `data` falls within the range (`inferior`, `superior`).

# Arguments
- `data::AbstractVector{<:RealOrQty}`: Elements to be checked.
- `inferior::RealOrQty = -Inf`: Lower limit.
- `superior::RealOrQty = Inf: Upper limit.
- `strict::Bool = true`: If the the edges are excluded from the range or not.

# Returns
- If every element x in `data` satisfy `inferior` < x < `superior`, for `strict` = true, or
  `inferior` <= x <= `superior`, for `strict` = false.

"""
function inRange(
    data::AbstractVector{<:RealOrQty},
    inferior::RealOrQty = -Inf,
    superior::RealOrQty = Inf;
    strict::Bool = true,
)::Bool

    if strict
        return maximum(data) < superior && minimum(data) > inferior
    else
        return maximum(data) <= superior && minimum(data) >= inferior
    end

end

"""
    longest(lists::AbstractVector{<:AbstractVector})::AbstractVector

Get the longest collection within `lists`.

# Arguments
- `lists::AbstractVector{<:AbstractVector}`: List of collections to be compared.

# Returns
- The longest array.

"""
@inline function longest(list::AbstractVector{<:AbstractVector})::AbstractVector

    return list[findmax(length, list)[2]]

end

"""
    internalUnits(
        quantity::String, 
        header::SnapshotHeader; 
        <keyword arguments>
    )::Union{Unitful.Quantity,Unitful.Units}

Get the factor to convert a plain number into a [Unitful](https://github.com/PainterQubits/Unitful.jl) quantity, with the correct internal units 
used by GADGET.

# Arguments
- `quantity::String`: Key to the target quantity. The options are given by [`QUANTITIES`](@ref)
  in `src/constants.jl`.
- `header::SnapshotHeader`: Header of the corresponding snapshot.
- `sim_cosmo::Bool = false`: If the simulation is cosmological, 
  - `false` ⟶ Newtonian simulation (`ComovingIntegrationOn` = 0).
  - `true` ⟶ Cosmological simulation (`ComovingIntegrationOn` = 1).

# Returns
- A [Unitful](https://github.com/PainterQubits/Unitful.jl) quantity or unit.

"""
function internalUnits(
    quantity::String,
    header::SnapshotHeader;
    sim_cosmo::Bool = false,
)::Union{Unitful.Quantity,Unitful.Units}

    dim = QUANTITIES[quantity].dimension
    unit = QUANTITIES[quantity].unit

    # Structure for unit conversion
    GU = GadgetPhysicalUnits(a_scale = sim_cosmo ? header.time : 1.0, hpar = header.h0)

    if unit == :internal

        if dim == Unitful.𝐌
            return GU.m_msun
        elseif dim == Unitful.𝐋
            return GU.x_physical
        elseif dim == Unitful.𝐓
            return sim_cosmo ? Unitful.NoUnits : GU.t_Myr
        elseif dim == Unitful.𝐌 / Unitful.𝐋^3
            return GU.rho_cgs
        elseif dim == SpecificEnergy
            return GU.E_cgs / GU.m_msun
        elseif dim == Unitful.𝐌 / Unitful.𝐓
            return GU.m_msun / GU.t_Myr
        elseif dim == Unitful.𝐋 / Unitful.𝐓
            return GU.v_physical
        elseif dim == Unitful.𝐋^2 / Unitful.𝐓^2
            return GU.E_cgs / GU.m_msun
        else
            error("No match could be found for the 'internal' units of $quantity.")
        end

    else
        return unit
    end

end

"""
    rangeCut!(
        data::AbstractVector{<:RealOrQty},  
        range::NTuple{2, <:Real};
        <keyword arguments>
    )::Int32

Delete every element in `data` that is outside of a given range.

# Arguments
- `data::AbstractVector{<:RealOrQty}`: Dataset that will be pruned outside the range.
- `range::NTuple{2,<:Real}`: The range in question. 
- `keep_edges::Bool = true`: If the edges of the range will be kept.
- `min_left::Int64 = 0`: Minimum number of values that need to be left after pruning to proceed 
  with the deletion.

# Returns
- Error flag:
  - `1` ⟶ There would be less than `min_left` values left in the dataset after pruning, 
    so the dataset will be kept untouched.
  - `0` ⟶ There would be more than `min_left` - 1 values left in the datasets after pruning,
    so the operation proceeds normally.

"""
function rangeCut!(
    data::AbstractVector{<:RealOrQty},
    range::NTuple{2,<:Real};
    keep_edges::Bool = true,
    min_left::Int64 = 0,
)::Int32

    l_range = range[1] * unit(data[1])
    h_range = range[2] * unit(data[1])

    if keep_edges
        # If the number of values within the range is less than `min_left`, 
        # return error code 1
        if count(x -> l_range <= x <= h_range, data) < min_left
            return 1
        else
            filter!(z -> l_range <= z <= h_range, data)
        end
    else
        # If the number of values within the range is less than `min_left`, 
        # return error code 1
        if count(x -> l_range < x < h_range, data) < min_left
            return 1
        else
            filter!(z -> l_range < z < h_range, data)
        end
    end

    return 0

end

"""
    rangeCut!(
        m_data::AbstractVector{<:RealOrQty}, 
        s_data::AbstractVector, 
        range::NTuple{2,<:Real};
        <keyword arguments>
    )::Int32

Delete every element in `m_data` that is outside of a given range. Every corresponding 
element in `s_data` (i.e. with the same index) will be deleted too.

# Arguments
- `m_data::AbstractVector{<:RealOrQty}`: Master dataset that will be pruned outside the range.
- `s_data::AbstractVector`: Slave dataset that will be pruned according to which values of `m_data` 
  are outside the range.
- `range::NTuple{2,<:Real}`: The range in question. 
- `keep_edges::Bool = true`: If the edges of the range will be kept.
- `min_left::Int64 = 0`: Minimum number of values that need to be left in the master dataset after 
  pruning to proceed with the deletion.

# Returns
- Error flag:
  - `1` ⟶ There would be less than `min_left` values left in the master dataset after pruning, 
    so the datasets will be kept untouched.
  - `0` ⟶ There would be more than `min_left` - 1 values left in the master datasets after pruning,
    so the operation proceeds normally.

"""
function rangeCut!(
    m_data::AbstractVector{<:RealOrQty},
    s_data::AbstractVector,
    range::NTuple{2,<:Real};
    keep_edges::Bool = true,
    min_left::Int64 = 0,
)::Int32

    length(s_data) >= length(m_data) || throw(ArgumentError(
        "`s_data` has to be at least as long as `m_data`."
    ))

    l_range = range[1] * unit(m_data[1])
    h_range = range[2] * unit(m_data[1])

    if keep_edges
        # If the number of values within the range is less than `min_left`, 
        # return error code 1
        if count(x -> l_range <= x <= h_range, m_data) < min_left
            return 1
        else
            deleteat!(s_data, .!collect(l_range .<= m_data .<= h_range))
            filter!(z -> l_range <= z <= h_range, m_data)
        end
    else
        # If the number of values within the range is less than `min_left`, 
        # return error code 1
        if count(x -> l_range < x < h_range, m_data) < min_left
            return 1
        else
            deleteat!(s_data, .!collect(l_range .< m_data .< h_range))
            filter!(z -> l_range < z < h_range, m_data)
        end
    end

    return 0

end

"""
    positiveCut!(
        data::AbstractVector{<:RealOrQty}; 
        <keyword arguments>
    )::Int32

Delete every element in `data` that is negative.

# Arguments
- `data::AbstractVector{<:RealOrQty}`: Dataset that will be pruned.
- `keep_edges::Bool = true`: If zero will be kept.
- `min_left::Int64 = 0`: Minimum number of values that need to be left after pruning to 
  proceed with the deletion.

# Returns
- Error flag:
  - `1` ⟶ There would be less than `min_left` values left in the dataset after pruning, 
    so the dataset will be kept untouched.
  - `0` ⟶ There would be more than `min_left` - 1 values left in the datasets after pruning,
    so the operation proceeds normally.

"""
function positiveCut!(
    data::AbstractVector{<:RealOrQty};
    keep_edges::Bool = true,
    min_left::Int64 = 0,
)::Int32

    return rangeCut!(data, (0.0, Inf); keep_edges, min_left)

end

"""
    positiveCut!(
        m_data::AbstractVector{<:RealOrQty}, 
        s_data::AbstractVector;
        <keyword arguments> 
    )::Int32

Delete every element in `m_data` that is negative. Every corresponding element in `s_data` 
(i.e. with the same index) will be deleted too.

# Arguments
- `m_data::AbstractVector{<:RealOrQty}`: Master dataset that will be pruned.
- `s_data::AbstractVector`: Slave dataset that will be pruned according to which values of 
  `m_data` are negative.
- `keep_edges::Bool = true`: If zero will be kept.
- `min_left::Int64 = 0`: Minimum number of values that need to be left in the master dataset after 
  pruning to proceed with the deletion.

# Returns
- Error flag:
  - `1` ⟶ There would be less than `min_left` values left in the master dataset after pruning, 
    so the datasets will be kept untouched.
  - `0` ⟶ There would be more than `min_left` - 1 values left in the master datasets after pruning,
    so the operation proceeds normally.

"""
function positiveCut!(
    m_data::AbstractVector{<:RealOrQty},
    s_data::AbstractVector;
    keep_edges::Bool = true,
    min_left::Int64 = 0,
)::Int32

    return rangeCut!(m_data, s_data, (0.0, Inf); keep_edges, min_left)

end

"""
    sanitizeData!(
        data::AbstractVector{<:RealOrQty};
        <keyword arguments> 
    )::NTuple{2,Int32}

Do three operations over `data`, in the following order:
  
  - Trim it to fit within the domain of `func_domain`.
  - Trim it to fit within `range`.
  - Scale it by a factor of 10^`exp_factor`.

By default, no transformation is done.

# Arguments
- `data::AbstractVector{<:RealOrQty}`: Dataset to be sanitized.
- `func_domain::Function = identity`: `data` will be trimmed to fit within the domain of the 
  function `func_domain`. The options are the scaling functions accepted by [Makie.jl](http://makie.juliaplots.org/stable/index.html): log10, log2, 
  log, sqrt, Makie.logit, Makie.Symlog10, Makie.pseudolog10 and identity.
- `range::NTuple{2,<:Real} = (-Inf, Inf)`: Every element in `data` that falls outside of `range` 
  will be deleted.
- `keep_edges::Bool = true`: If the edges of `range` will be kept.
- `min_left::Int64 = 0`: Minimum number of values that need to be left after pruning to proceed 
  with the deletion.
- `exp_factor::Int64 = 0`: Every element in `data` will be multiplied by 10^`exp_factor`.

# Returns
- A tuple with two error flags:
  - After trying to fit `data` within the domain of `func_domain`:
    - `1` ⟶ There would be less than `min_left` values left in the dataset after pruning, 
      so the dataset will be kept untouched.
    - `0` ⟶ There would be more than `min_left` - 1 values left in the datasets after pruning, 
      so the operation proceeds normally.
  - After trying to fit the `data` in `range`:
    - `1` ⟶ There would be less than `min_left` values left in the dataset after pruning, 
      so the dataset will be kept untouched.
    - `0` ⟶ There would be more than `min_left` - 1 values left in the datasets after pruning, 
      so the operation proceeds normally.

"""
function sanitizeData!(
    data::AbstractVector{<:RealOrQty};
    func_domain::Function = identity,
    range::NTuple{2,<:Real} = (-Inf, Inf),
    keep_edges::Bool = true,
    min_left::Int64 = 0,
    exp_factor::Int64 = 0,
)::NTuple{2,Int32}

    !(isa(data, AbstractVector{<:Integer}) && exp_factor != 0) || @warn(
        "Elements of `data` are of type Integer, this may result in errors or unwanted truncation \
        when using `exp_factor` != 0."
    )

    if func_domain != identity
        if func_domain == sqrt
            domain_flag = positiveCut!(data; min_left)
        elseif func_domain == Makie.logit
            domain_flag = rangeCut!(data, (0, 1); keep_edges = false, min_left)
        elseif func_domain == log || func_domain == log2 || func_domain == log10
            domain_flag = positiveCut!(data; keep_edges = false, min_left)
        else
            throw(ArgumentError(
                "Function $func_domain is not supported. See list of suppoerted scaling functions \
                in the [Makie.jl](http://makie.juliaplots.org/stable/index.html) documentation."
            ))
        end
    else
        domain_flag = 0
    end

    range_flag = rangeCut!(data, range; keep_edges, min_left)

    data ./= 10^exp_factor

    return domain_flag, range_flag

end

"""
    sanitizeData!(
        x_data::AbstractVector{<:RealOrQty},
        y_data::AbstractVector{<:RealOrQty};
        <keyword arguments> 
    )::NTuple{4,Int32}

Do three operations over `x_data` and `y_data`, in the following order:

  - Trim them to fit within the domain of `func_domain[1]` and `func_domain[2]` respectively.
  - Trim them to fit within `range[1]` and `range[2]` respectively.
  - Scale them by the factors 10^`exp_factor[1]` and 10^`exp_factor[2]` respectively.

By default, no transformation is done to any of the datasets.

!!! note 
    The dataset must have the same length, and any operation that deletes an element, will 
    delete the corresponding element (i.e. with the same index) in the other dataset, 
    so the dataset will be kept equal in length at the end.

# Arguments
- `x_data::AbstractVector{<:RealOrQty}`: First dataset to be sanitized.
- `y_data::AbstractVector{<:RealOrQty}`: Second dataset to be sanitized.
- `func_domain::NTuple{2, Function} = (identity, identity)`: `x_data` will be trimmed to fit within 
  the domain of the function `func_domain[1]`, and `y_data` will be trimmed to fit within the domain 
  of the function `func_domain[2]`. The possible options are the scaling functions accepted by 
  [Makie.jl](http://makie.juliaplots.org/stable/index.html): log10, log2, log, sqrt, Makie.logit, Makie.Symlog10, Makie.pseudolog10 and identity.
- `range::NTuple{2,NTuple{2,<:Real}} = ((-Inf, Inf), (-Inf, Inf))`: Every element in `x_data` that 
  falls outside of `range[1]` will be deleted, and every element in `y_data` that falls outside of 
  `range[2]` will be deleted.
- `keep_edges::NTuple{2,Bool} = (true, true)`: If the edges of each corresponding `range` will be 
  kept.
- `min_left::Int64 = 0`: Minimum number of values that need to be left after prunning to proceed 
  with the deletion.
- `exp_factor::NTuple{2,Int64} = (0, 0)`: Every element in `x_data` will be multiplied by 
  10^`exp_factor[1]`, and every element in `y_data` will be multiplied by 10^`exp_factor[2]`.

# Returns
- A tuple with four error flags:
  - After trying to fit `x_data` and `y_data` within the domain of `func_domain`, one flag for 
    every dataset:
    - `1` ⟶ There would be less than `min_left` values left in the dataset after pruning, 
      so the dataset will be kept untouched.
    - `0` ⟶ There would be more than `min_left` - 1 values left in the datasets after 
      pruning, so the operation proceeds normally.
  - After trying to fit `x_data` and `y_data` in `range`, one flag for every dataset:
    - `1` ⟶ There would be less than `min_left` values left in the dataset after pruning, 
      so the dataset will be kept untouched.
    - `0` ⟶ There would be more than `min_left` - 1 values left in the datasets after 
      pruning, so the operation proceeds normally.

"""
function sanitizeData!(
    x_data::AbstractVector{<:RealOrQty},
    y_data::AbstractVector{<:RealOrQty};
    func_domain::NTuple{2,Function} = (identity, identity),
    range::NTuple{2,NTuple{2,<:Real}} = ((-Inf, Inf), (-Inf, Inf)),
    keep_edges::NTuple{2,Bool} = (true, true),
    min_left::Int64 = 0,
    exp_factor::NTuple{2,Int64} = (0, 0),
)::NTuple{4,Int32}

    length(x_data) == length(y_data) || throw(ArgumentError(
        "`x_data` and `y_data` must have the same legth."
    ))
    !(isa(x_data, AbstractVector{<:Integer}) && exp_factor[1] != 0) || @warn(
        "Elements of `x_data` are of type Integer, this may result in errors or unwanted \
        truncation when using `exp_factor[1]` != 0."
    )
    !(isa(y_data, AbstractVector{<:Integer}) && exp_factor[2] != 0) || @warn(
        "Elements of `y_data` are of type Integer, this may result in errors or unwanted \
        truncation when using `exp_factor[2]` != 0."
    )
    

    if func_domain[1] != identity
        if func_domain[1] == sqrt
            x_domain_flag = positiveCut!(x_data, y_data; min_left)
        elseif func_domain[1] == Makie.logit
            x_domain_flag = rangeCut!(x_data, y_data, (0, 1); keep_edges = false, min_left)
        elseif func_domain[1] == log || func_domain[1] == log2 || func_domain[1] == log10
            x_domain_flag = positiveCut!(x_data, y_data; keep_edges = false, min_left)
        else
            throw(ArgumentError(
                "Function $(func_domain[1]) is not supported. See list of suppoerted scaling \
                functions in the [Makie.jl](http://makie.juliaplots.org/stable/index.html) documentation."
            ))
        end
    else
        x_domain_flag = 0
    end

    if func_domain[2] != identity
        if func_domain[2] == sqrt
            y_domain_flag = positiveCut!(y_data, x_data; min_left)
        elseif func_domain[2] == Makie.logit
            y_domain_flag = rangeCut!(y_data, x_data, (0, 1); keep_edges = false, min_left)
        elseif func_domain[2] == log || func_domain[2] == log2 || func_domain[2] == log10
            y_domain_flag = positiveCut!(y_data, x_data; keep_edges = false, min_left)
        else
            throw(ArgumentError(
                "Function $(func_domain[2]) is not supported. See list of suppoerted scaling \
                functions in the [Makie.jl](http://makie.juliaplots.org/stable/index.html) documentation."
            ))
        end
    else
        y_domain_flag = 0
    end

    x_range_flag = rangeCut!(x_data, y_data, range[1]; keep_edges = keep_edges[1], min_left)
    y_range_flag = rangeCut!(y_data, x_data, range[2]; keep_edges = keep_edges[2], min_left)

    x_data ./= 10^exp_factor[1]
    y_data ./= 10^exp_factor[2]

    return x_domain_flag, y_domain_flag, x_range_flag, y_range_flag

end

"Total mass of metals in a every particle of a GADGET simulation."
metalMass(metals::Matrix{<:Unitful.Mass}) = [
    # Add up all the elements except the ones at position 1 and 7, i.e. He and H
    sum(metals[[2, 3, 4, 5, 6, 8, 9, 10, 11, 12], i]) for i in axes(metals, 2)
]

"""
    computeMetallicity(
        metals::Matrix{<:Unitful.Mass},
        mass::Vector{<:Unitful.Mass}; 
        <keyword arguments>
    )::Vector{Float64}

Compute the metallicity as the dimensionless fraction of metal mass inside every gas or stellar 
particle in a snapshot.

# Arguments
- `metals::Matrix{<:Unitful.Mass}`: Matrix with the mass content of every element. Each row 
  must be an element and each column a particle.
  - Row 01: He (Helium)
  - Row 02: C (Carbon)
  - Row 03: Mg (Magnesium)
  - Row 04: 0 (Oxygen)
  - Row 05: Fe (Iron)
  - Row 06: Si (Silicon)
  - Row 07: H (Hydrogen)
  - Row 08: N (Nitrogen)
  - Row 09: Ne (Neon)
  - Row 10: S (Sulfur)
  - Row 11: Ca (Calcium)
  - Row 12: Zn (Zinc)
- `mass::Vector{<:Unitful.Mass}`: Mass of every gas or stellar particle.
- `solar::Bool = false`: If the result will be in units of solar metallicity or not.

# Returns
- The metallicity of the particles.  

"""
function computeMetallicity(
    metals::Matrix{<:Unitful.Mass},
    mass::Vector{<:Unitful.Mass};
    solar::Bool = false,
)::Vector{Float64}

    metallicity = uconvert.(Unitful.NoUnits, metalMass(metals) ./ mass)

    return solar ? metallicity ./ SOLAR_METALLICITY : metallicity

end

"""
    computeElementFraction(
        metals::Matrix{<:Unitful.Mass},
        element::String; 
        <keyword arguments>
    )::Vector{Float64}

Compute the total fraction of `element` with respec to Hydrogen, and apply `function` to it, for 
every gas or stellar particle in a snapshot.

# Arguments
- `metals::Matrix{<:Unitful.Mass}`: Matrix with the mass content of every element. Each row 
  must be an element and each column a particle.
  - Row 01: He (Helium)
  - Row 02: C (Carbon)
  - Row 03: Mg (Magnesium)
  - Row 04: 0 (Oxygen)
  - Row 05: Fe (Iron)
  - Row 06: Si (Silicon)
  - Row 07: H (Hydrogen)
  - Row 08: N (Nitrogen)
  - Row 09: Ne (Neon)
  - Row 10: S (Sulfur)
  - Row 11: Ca (Calcium)
  - Row 12: Zn (Zinc)
- `element::String`: Name of the element which fraction will be calculated.
- `func::Function = identity`: Function to be applied to the fraction.

# Returns
- `function`(X / H) where X is `element`.  

"""
function computeElementFraction(
    metals::Matrix{<:Unitful.Mass},
    element::String;
    func::Function = identity,
)::Vector{Float64}

    if element == "He" || element == "Helium"
        index = 1
    elseif element == "C" || element == "Carbon"
        index = 2
    elseif element == "Mg" || element == "Magnesium"
        index = 3
    elseif element == "O" || element == "Oxygen"
        index = 4
    elseif element == "Fe" || element == "Iron"
        index = 5
    elseif element == "Si" || element == "Silicon"
        index = 6
    elseif element == "H" || element == "Hydrogen"
        index = 7
    elseif element == "N" || element == "Nitrogen"
        index = 8
    elseif element == "Ne" || element == "Neon"
        index = 9
    elseif element == "S" || element == "Sulfur"
        index = 10
    elseif element == "Ca" || element == "Calcium"
        index = 11
    elseif element == "Zn" || element == "Zinc"
        index = 12
    else
        throw(ArgumentError("Element $element is not supported by GADGET"))
    end

    return @. func(metals[index, :] / metals[7, :])

end

@doc raw"""
    computeTime(
        scale_factors::Vector{<:Real}, 
        header::GadgetIO.SnapshotHeader; 
        <keyword arguments>
    )::Vector{<:Unitful.Time}

Compute the clock time, from the scale factor.

To get the clock time from the scale factor one has to do the integral

```math
t = H_0^{\,-1} \int_0^a \frac{\mathrm{d}a'}{a' \, \sqrt{\mathcal{E}(a')}} \, ,
```

where

```math
\mathcal{E}(a) = \Omega_\Lambda + \Omega_m \, a^{-3} + \Omega_r \, a^{-4} + \Omega_K \, a^{-2} \, .
```

# Arguments
- `scale_factors::Vector{<:Real}`: Scale factors to be converted.
- `header::GadgetIO.SnapshotHeader`: Header of some snapshot of the simulation (it contains the 
  cosmological parameters). 
- `a0::Real = scale_factor[1]`: Initial scale factor.

# Returns
- The clock time.   

"""
function computeTime(
    scale_factors::Vector{<:Real},
    header::GadgetIO.SnapshotHeader;
    a0::Real = scale_factors[1],
)::Vector{<:Unitful.Time}

    f = x -> energyIntegrand(header, x)

    return [quadgk(f, a0, a)[1] * UnitfulAstro.Gyr for a in scale_factors]

end

"""
    computeRedshift(scale_factor::Vector{<:Real})::Vector{Float64}

Compute the redshift, from the scale factor.

# Arguments
- `scale_factor::Vector{<:Real}`: Scale factors to be converted.

# Returns
- The redshift.  

"""
function computeRedshift(scale_factor::Vector{<:Real})::Vector{Float64}

    return @. ((1.0 / scale_factor) - 1.0)

end

"""
    computeStellarAge(
        birth_times::Vector{<:Unitful.Time}, 
        ::SimData,
        snap_data::SnapData,
    )::Vector{<:Unitful.Time}

Compute the age of the stars.

# Arguments
- `birth_times::Vector{<:Unitful.Time}`: Time of birth of the stars.
- `::SimData`: For compatibility with the other method.
- `snap_data::SnapData`: Basic data of the snapshot containing the stars in question, including its 
  timestamp, header, etc.

# Returns
- The stellar ages.   

"""
function computeStellarAge(
    birth_times::Vector{<:Unitful.Time},
    ::SimData,
    snap_data::SnapData,
)::Vector{<:Unitful.Time}

    return snap_data.time_stamp .- birth_times

end

"""
    computeStellarAge(
        birth_a::Vector{<:Real},
        sim_data::SnapData, 
        snap_data::SnapData,
    )::Vector{<:Unitful.Time}

Compute the age of the stars.

# Arguments
- `birth_times::Vector{<:Real}`: Scale factor at the moment of birth of the stars.
- `sim_data::SnapData`: Basic information about the simulation.
- `snap_data::SnapData`: Basic data of the snapshot containing the stars in question, including its 
  timestamp.

# Returns
- The stellar ages.   

"""
function computeStellarAge(
    birth_a::Vector{<:Real},
    sim_data::SimData,
    snap_data::SnapData,
)::Vector{<:Unitful.Time}

    # From scale factor to clock time
    birth_times = computeTime(birth_a, sim_data.header, a0 = sim_data.a0)

    return snap_data.time_stamp .- birth_times

end

"""
    computeSFR(
        stellar_mass::Vector{<:Unitful.Mass}, 
        time::Vector{<:Unitful.Time},
    )::Vector{<:Unitful.MassFlow}

Compute the total star formation rate of a series of snapshots, and replace negative values with 0.

# Arguments
- `stellar_mass::Vector{<:Unitful.Mass}`: Total mass of the stellar particles, one value per 
  snapshot.
- `time::Vector{<:Unitful.Time}`: Clock time of each snapshot.

# Returns
- The star formation rate of each snapshot. 

"""
function computeSFR(
    stellar_mass::Vector{<:Unitful.Mass},
    time::Vector{<:Unitful.Time},
)::Vector{<:Unitful.MassFlow}

    sfr = [
        (stellar_mass[i] - stellar_mass[i-1]) / (time[i] - time[i-1]) for i in 2:length(stellar_mass)
    ]

    u_zero = zero(sfr[1])
    replace!(x -> x < u_zero ? u_zero : x, sfr)

    return prepend!(sfr, u_zero)

end

"""
    computeDistance(
        positions::Matrix{<:Unitful.Length}; 
        <keyword arguments>
    )::Vector{<:Unitful.Length}

Compute the distance of several particles from a `center`.

# Arguments
- `positions::Matrix{<:Unitful.Length}`: Coordinates of the particles, where each row is an axis, 
  and each column is a particle.
- `center::Union{Vector{<:Unitful.Length},Nothing} = nothing`: Origin used to calculate the 
  distances.

# Returns
- The distance of every particle to the `center`.  

"""
function computeDistance(
    positions::Matrix{<:Unitful.Length};
    center::Union{Vector{<:Unitful.Length},Nothing} = nothing,
)::Vector{<:Unitful.Length}

    if center === nothing
        return [norm(col) for col in eachcol(positions)]
    else
        length(center) == size(positions, 1) || throw(ArgumentError(
            "`center` must have as many elements as `positions` has rows."
        ))

        return [norm(col .- center) for col in eachcol(positions)]
    end

end

"""
    computeCenterOfMass(
        positions::Matrix{<:Unitful.Length},
        mass::Vector{<:Unitful.Mass},
    )::Vector{<:Unitful.Length}

Compute the center of mass of a group of particles.

# Arguments
- `positions::Matrix{<:Unitful.Length}`: Coordinates of the particles, where each row is an axis, 
  and each column is a particle.
- `mass::Vector{<:Unitful.Mass}`: Masses of the particles.

# Returns
- The center of mass.   

"""
function computeCenterOfMass(
    positions::Matrix{<:Unitful.Length},
    mass::Vector{<:Unitful.Mass},
)::Vector{<:Unitful.Length}

    R = [sum(positions[i, :] .* mass) for i in size(positions, 1)]

    center_of_mass = R ./ sum(mass)

    # If the norm of the center of mass is less than 0.1% 
    # of the largest distance in the dataset, return 0
    if norm(center_of_mass) < maxLength(positions) * 0.001
        return zeros(typeof(positions[1]), size(positions, 1))
    end

    return center_of_mass

end

"""
    computeTemperature(
        metals::Matrix{<:Unitful.Mass},
        mass::Vector{<:Unitful.Mass},
        internal_energy::Vector{<:SpecificEnergy},
        electron_fraction::Vector{<:Real},
    )::Vector{<:Unitful.Temperature}

Compute the temperature of gas particles.

# Arguments
- `metals::Matrix{<:Unitful.Mass}`: Matrix with the mass content of every element. Each row 
  must be an element and each column a particle.
  - Row 01: He (Helium)
  - Row 02: C (Carbon)
  - Row 03: Mg (Magnesium)
  - Row 04: 0 (Oxygen)
  - Row 05: Fe (Iron)
  - Row 06: Si (Silicon)
  - Row 07: H (Hydrogen)
  - Row 08: N (Nitrogen)
  - Row 09: Ne (Neon)
  - Row 10: S (Sulfur)
  - Row 11: Ca (Calcium)
  - Row 12: Zn (Zinc)
- `mass::Vector{<:Unitful.Mass}`: Mass of every gas particle.
- `internal_energy::Vector{<:SpecificEnergy}`: Specific internal energy of every gas particles.
- `electron_fraction::Vector{<:Real}`: Number fraction of electrons in every gas particle.

# Returns
- The temperature of each gas particle.   

"""
function computeTemperature(
    metals::Matrix{<:Unitful.Mass},
    mass::Vector{<:Unitful.Mass},
    internal_energy::Vector{<:SpecificEnergy},
    electron_fraction::Vector{<:Real},
)::Vector{<:Unitful.Temperature}

    # xH := mass_fraction_of_Hydrogen
    xH = metals[7, :] ./ mass

    # yHe := number_of_Helium_atoms / number_of_Hydrogen_atoms
    # Here we take the mass fraction of metals as negligible
    yHe = @. (1.0 - xH) / (4.0 * xH)

    # electron_fraction := number_of_electrons / number_of_Hydrogen_atoms
    # μ := total_mass / (total_number_of_particles * proton_mass)
    #   ≈ number_of_protons / total_number_of_particles
    # For the total mass we take the mass of electrons as negligible
    μ = @. (1.0 + 4.0 * yHe) / (1.0 + yHe + electron_fraction)

    # T = (adiabatic_index - 1) * internal_energy_per_unit_mass * 
    #     (total_mass / total_number_of_particles) / Boltzmann_constant
    return @. 2.0 / 3.0 * internal_energy * μ * Unitful.mp / Unitful.k

end

"""
    computeSurfaceDensity(
        distances::Vector{<:Unitful.Length},
        quantity::Vector{<:RealOrQty},
        max_radius::Unitful.Length,
        bins::Int64,
    )::Vector{<:Unitful.Quantity}

Compute the surface density of a given `quantity`, using concentric rings of equal width.

# Arguments
- `distances::Vector{<:Unitful.Length}`: Projected 2D distance to each particle.
- `quantity::Vector{<:RealOrQty}`: Quantity for which the density will be 
  calculated. It should be defined as a single value per particle.
- `max_radius::Unitful.Length`: Maximum distance up to which the density will be calculated.
- `bins::Int64`: Number of concentric rings.

# Returns
- The density of `quantity` within each ring.   

"""
function computeSurfaceDensity(
    distances::Vector{<:Unitful.Length},
    quantity::Vector{<:RealOrQty},
    max_radius::Unitful.Length,
    bins::Int64,
)::Vector{<:Unitful.Quantity}

    # Width of each ring
    width = max_radius / bins

    surface_density = Vector{Unitful.Quantity}(undef, bins)

    @inbounds for i in eachindex(surface_density)

        # Indices of the particles within ring `i`
        idx = findall(x -> width * (i - 1) <= x < width * i, distances)
        # Ring area
        a = π * width * width * (2.0 * i - 1.0)

        surface_density[i] = sum(quantity[idx]) / a

    end

    return surface_density

end

"""
    computeSurfaceDensity(
        positions::Matrix{<:Unitful.Length},
        quantity::Vector{<:RealOrQty},
        radius::Unitful.Length,
        bins::Int64,
    )::Matrix{<:Unitful.Quantity}

Compute the surface density of a given `quantity`, using square pixels of equal size.

# Arguments
- `positions::Matrix{<:Unitful.Length}`: Projected 2D position of each particle. Each row is an 
  axis and each column is a particle.
- `quantity::Vector{<:RealOrQty}`: Quantity for which the density will be calculated. It should be 
  defined as a single value per particle.
- `radius::Unitful.Length`: Maximum distance up to which the density will be calculated.
- `bins::Int64`: Number of rows and columns used to divide the 
  (-`radius`, `radius`] x (-`radius`, `radius`] region.

# Returns
- The density of `quantity` within each pixel. 

"""
function computeSurfaceDensity(
    positions::Matrix{<:Unitful.Length},
    quantity::Vector{<:RealOrQty},
    radius::Unitful.Length,
    bins::Int64,
)::Matrix{<:Unitful.Quantity}

    # Width of each pixel
    width = 2 * radius / bins
    # Pixel area
    a = width * width

    surface_density = Matrix{Unitful.Quantity}(undef, (bins, bins))
    columns = collect(eachcol(positions))

    @inbounds for i in eachindex(surface_density)

        x = div(i - 1, bins) * width
        y = (mod1(i, bins) - 1) * width

        # Indices of the particles within pixel `i`
        idx = findall(
            col -> (
                (x - radius < col[1] <= x + width - radius) &&
                (radius - y < col[2] <= radius - y - width)
            ),
            columns,
        )

        surface_density[i] = sum(quantity[idx]) / a

    end

    return surface_density

end

"""
    computeProfile(
        distances::Vector{<:Unitful.Length}, 
        quantity::Vector{<:RealOrQty},
        max_radius::Unitful.Length,
        bins::Int64;
        <keyword arguments>
    )::NTuple{2,Vector}

Compute a profile of a given `quantity`, using rings or spherical shells of equal width.

# Arguments
- `distances::Vector{<:Unitful.Length}`: Distance to each particle.
- `quantity::Vector{<:RealOrQty}`: Quantity for which the profile will be calculated. It should be 
  defined as a single value per particle.
- `max_radius::Unitful.Length`: Maximum distance up to which the profile will be calculated.
- `bins::Int64`: Number of rings or shells.
- `flat::Bool = true`: If the profile will be 2D, using rings (default), or 3D, using spherical 
  shells.
- `cumulative::Bool = false`: If the profile will be accumulated or not.
- `density::Bool = false`: If the profile is of the area/volume density of `quantity` or just of 
  `quantity` itself.

# Returns
- A Tuple with two elements:
  - A Vector with the distances to each ring or shell.
  - A Vector with the values of the profile at those distances.  

"""
function computeProfile(
    distances::Vector{<:Unitful.Length},
    quantity::Vector{<:RealOrQty},
    max_radius::Unitful.Length,
    bins::Int64;
    flat::Bool = true,
    cumulative::Bool = false,
    density::Bool = false,
)::NTuple{2,Vector}

    # Width of each ring or spherical shell
    width = max_radius / bins

    # Initialize output arrays
    x_data = Vector{Unitful.Quantity}(undef, bins)
    y_data = Vector{RealOrQty}(undef, bins)
    region = Vector{Unitful.Quantity}(undef, bins)

    @inbounds for i in eachindex(x_data, y_data)

        # Indices of the particles which fall within the window `i`
        idx = findall(x -> width * (i - 1) <= x < width * i, distances)

        # Mean distance for the window `i`
        x_data[i] = isempty(idx) ? width * (i - 0.5) : sum(distances[idx]) / length(idx)
        y_data[i] = sum(quantity[idx])

        # Area of the ring `i`
        area = π * width^2.0 * (2.0 * i - 1.0)
        # Volume of the spherical shell `i`
        volume = 4.0 / 3.0 * π * width^3.0 * (3.0 * i * i - 3.0 * i + 1.0)

        region[i] = flat ? area : volume

    end

    if cumulative && !density
        return x_data, cumsum(y_data)
    elseif !cumulative && density
        return x_data, y_data ./ region
    elseif cumulative && density
        return x_data, cumsum(y_data) ./ cumsum(region)
    else
        return x_data, y_data
    end

end

"""
    computeZProfile(
        distances::Vector{<:Unitful.Length}, 
        metals::Matrix{<:Unitful.Mass},
        mass::Vector{<:Unitful.Mass},
        max_radius::Unitful.Length,
        bins::Int64;
        <keyword arguments>
    )::NTuple{2,Vector}

Compute a metallicity profile.

# Arguments
- `distances::Vector{<:Unitful.Length}`: Distance to each particle.
- `metals::Matrix{<:Unitful.Mass}`: Metallicity of the particles.
- `mass::Vector{<:Unitful.Mass}`: Masses of the particles.
- `max_radius::Unitful.Length`: Maximum distance up to which the profile will be calculated.
- `bins::Int64`: Number of subdivisions of [0, `max_radius`].
- `flat::Bool = true`: If the profile will be 2D, using rings (default), or 3D, using spherical 
  shells.
- `cumulative::Bool = false`: If the profile will be accumulated or not.
- `density::Bool = false`: If the profile is of the area/volume density of `quantity` or just of 
  `quantity` itself.
- `solar::Bool = false`: If the result will be in units of solar metallicity or not.

# Returns
- A Tuple with two elements:
  - A Vector with the distances to each ring or shell.
  - A Vector with the values of the profile at those distances.  

"""
function computeZProfile(
    distances::Vector{<:Unitful.Length},
    metals::Matrix{<:Unitful.Mass},
    mass::Vector{<:Unitful.Mass},
    max_radius::Unitful.Length,
    bins::Int64;
    flat::Bool = true,
    cumulative::Bool = false,
    density::Bool = false,
    solar::Bool = false,
)::NTuple{2,Vector}

    # Compute the total mass of metals for each particle
    Z = metalMass(metals)

    # Width of each ring shell
    width = max_radius / bins

    # Initialize output arrays
    x_data = Vector{RealOrQty}(undef, bins)
    y_data = Vector{RealOrQty}(undef, bins)
    region = Vector{Unitful.Quantity}(undef, bins)

    @inbounds for i in eachindex(x_data, y_data)

        # Indices of the particles which fall within the window `i`
        idx = findall(x -> width * (i - 1) <= x < width * i, distances)

        # Mean distance for the window `i`
        x_data[i] = isempty(idx) ? width * (i - 0.5) : sum(distances[idx]) / length(idx)

        total_mass = sum(mass[idx])
        total_z = sum(Z[idx])

        # Metallicity for window `i`
        y_data[i] = isempty(idx) ? 0.0 : (total_z / total_mass) / (solar ? SOLAR_METALLICITY : 1.0)

        # Area of the ring `i`
        area = π * width^2.0 * (2.0 * i - 1.0)
        # Volume of the spherical shell `i`
        volume = 4.0 / 3.0 * π * width^3.0 * (3.0 * i * i - 3.0 * i + 1.0)

        region[i] = flat ? area : volume

    end

    if cumulative && !density
        return x_data, cumsum(y_data)
    elseif !cumulative && density
        return x_data, y_data ./ region
    elseif cumulative && density
        return x_data, cumsum(y_data) ./ cumsum(region)
    else
        return x_data, y_data
    end

end

"""
    computeTimeSeries(
        snapshot_files::Vector{String},
        quantity::String,
        type::Union{Symbol,Nothing},
        sim_cosmo::Bool; 
        <keyword arguments>
    )::Vector

Compute the time series of a given `quantity`.

# Arguments
- `snapshot_files::Vector{String}`: Paths to the snapshot files.
- `quantity::String`: Some of the following quantities:
  - `"clock_time"`   ⟶ Physical time of each snapshot (time units).
  - `"scale_factor"` ⟶ Scale factor of each snapshot (dimensionless).
  - `"redshift"`     ⟶ Redshift of each snapshot (dimensionless).
  - `"number"`       ⟶ Number of particles of a given type (dimensionless).
  - `"mass"`         ⟶ Total mass of a given type of particle (mass units).
  - `"sfr"`          ⟶ Star formation rate of each snapshot (mass over time units).
- `type::Union{Symbol,Nothing}`: Particle type. The possibilities are given by [`ParticleType`](@ref) 
  in `src/constants.jl`. It is only relevant for the `"number"` and `"mass"` quantities.
- `sim_cosmo::Bool`: If the simulation is cosmological, 
  - `false` ⟶ Newtonian simulation (`ComovingIntegrationOn` = 0).
  - `true` ⟶ Cosmological simulation (`ComovingIntegrationOn` = 1).
- `filter_function::Union{Function,Nothing} = nothing`: A function with the signature: 

  `foo(file_path::String)::Vector{Int64}`
  
  It indicates which particles will be read, taking the file path to a snapshot and returning the 
  list of indices of the selected particles. If set to `nothing`, then no particles are filtered. 
  See the [GadgetIO.jl](https://ludwigboess.github.io/GadgetIO.jl/stable/read_snapshots/#Filter-functions) documentation for examples.
- `warnings::Bool = true`: If a warning will be given when the data is not as expected, but the 
  function can still run using sane defaults.

# Returns
- The time series. 

"""
function computeTimeSeries(
    snapshot_files::Vector{String},
    quantity::String,
    type::Union{Symbol,Nothing},
    sim_cosmo::Bool;
    filter_function::Union{Function,Nothing} = nothing,
    warnings::Bool = true,
)::Vector

    !(quantity ∈ ["number", "mass"] && type === nothing) || throw(ArgumentError(
        "For quantities 'number' and 'mass', you have to specify a `type`."
    ))

    headers = read_header.(snapshot_files)

    if quantity == "clock_time"

        if sim_cosmo
            a = getfield.(headers, :time)
            qty_out = computeTime(a, headers[1])
        else
            qty_out = getfield.(headers, :time) .* internalUnits.("CLOCK_TIME", headers; sim_cosmo)
        end

    elseif quantity == "scale_factor"

        if sim_cosmo
            qty_out = getfield.(headers, :time)
        else
            qty_out = ones(length(snapshot_files))
        end

    elseif quantity == "redshift"

        if sim_cosmo
            qty_out = computeRedshift(getfield.(headers, :time))
        else
            qty_out = zeros(length(snapshot_files))
        end

    elseif quantity == "number"

        qty_out = getindex.(getfield.(headers, :nall), ParticleType[type] + 1)

    elseif quantity == "mass"

        qty_out = Vector{Unitful.Mass}(undef, length(snapshot_files))

        for (i, (header, snapshot)) in enumerate(zip(headers, snapshot_files))

            unit = internalUnits("MASS", header; sim_cosmo)
            number = header.nall[ParticleType[type]+1]

            if number != 0
                if header.massarr[ParticleType[type]+1] != 0
                    # If all `type` particles have the same mass
                    qty_out[i] = header.massarr[ParticleType[type]+1] * number * unit
                else
                    qty_out[i] = sum(
                        getSnapshotData(
                            snapshot,
                            type,
                            "MASS";
                            sim_cosmo,
                            filter_function,
                            warnings,
                        )[type]["MASS"],
                    )
                end
            else
                # In the case that there are no `type` particles in the snapshot
                qty_out[i] = 0.0 * unit
            end

        end

    elseif quantity == "sfr"

        stellar_mass = Vector{Unitful.Mass}(undef, length(snapshot_files))

        for (i, (header, snapshot)) in enumerate(zip(headers, snapshot_files))

            unit = internalUnits("MASS", header; sim_cosmo)
            number = header.nall[ParticleType[:stars]+1]

            if number != 0
                if header.massarr[ParticleType[:stars]+1] != 0
                    # If all stellar particles have the same mass
                    stellar_mass[i] = header.massarr[ParticleType[:stars]+1] * number * unit
                else
                    stellar_mass[i] = sum(
                        getSnapshotData(
                            snapshot,
                            :stars,
                            "MASS";
                            sim_cosmo,
                            filter_function,
                            warnings,
                        )[:stars]["MASS"],
                    )
                end
            else
                # In the case that there are no stellar particles in the snapshot
                stellar_mass[i] = 0.0 * unit
            end

        end

        if sim_cosmo
            a = getfield.(headers, :time)
            time = computeTime(a, headers[1])
        else
            time = getfield.(headers, :time) .* internalUnits.("CLOCK_TIME", headers; sim_cosmo)
        end

        qty_out = computeSFR(stellar_mass, time)

    else
        throw(ArgumentError(
            "The quantity `$quantity` is not supported. See the documentation for the lits of \
            possibilities."
        ))
    end

    return qty_out

end

"""
    checkDataShape(::typeof(hist!), data_length::Int64, sim_number::Int64)::Nothing

Throw an error if `data_length` is not compatible with the plot type (heatmap for this method) or 
if there is more than one simulation to plot.

# Arguments
- `::typeof(hist!)`: Plot function.
- `data_length::Int64`: Length of the data tuple returned by `data_analysis`.
- `sim_number::Int64`: Number of simulations.

"""
function checkDataShape(::typeof(heatmap!), data_length::Int64, sim_number::Int64)::Nothing

    data_length == 3 || error(
        "For heatmaps `data_analysis` should return three data vectors, and currently is returning \
        $data_length."
    )

    sim_number == 1 || error("For heatmaps only one simulation at a time can be plotted.")

    return nothing

end

"""
    checkDataShape(::typeof(hist!), data_length::Int64, sim_number::Int64)::Nothing

Throw an error if `data_length` is not compatible with the plot type (histogram for this method).

# Arguments
- `::typeof(hist!)`: Plot function.
- `data_length::Int64`: Length of the data tuple returned by `data_analysis`.
- `sim_number::Int64`: Number of simulations. In this method, it is only for compatibility.

"""
function checkDataShape(::typeof(hist!), data_length::Int64, sim_number::Int64)::Nothing

    data_length == 1 || error(
        "For histograms `data_analysis` should return only one data vector, and currently is \
        returning $data_length."
    )

    return nothing

end

"""
    checkDataShape(
        ::Union{typeof(scatter!),typeof(scatterlines!),typeof(lines!)}, 
        data_length::Int64, 
        ::Int64,
    )::Nothing

Throw an error if `data_length` is not compatible with the plot type (scatter/lines for this method).

# Arguments
- `::Union{typeof(scatter!),typeof(scatterlines!),typeof(lines!)}`: Plot function.
- `data_length::Int64`: Length of the data tuple returned by `data_analysis`.
- `::Int64`: Number of simulations (only for compatibility with the other methods).

"""
function checkDataShape(
    ::Union{typeof(scatter!),typeof(scatterlines!),typeof(lines!)},
    data_length::Int64,
    ::Int64,
)::Nothing

    data_length == 2 || error(
        "For scatter/line plots `data_analysis` should return two data vectors, and currently is \
        returning $data_length."
    )

    return nothing

end
