####################################################################################################
# Data analysis functions
####################################################################################################

####################################################################################################
# Signature for data analysis functions targeting plotSnapshot() in ./src/plotting/pipelines.jl
####################################################################################################
#
# A data analysis functions targeting plotSnapshot() must take a data dictionary (see makeDataDict()
# in ./src/analysis/data_acquisition.jl for the canonical description) and return one or more
# vectors or matrices. It should return `nothing` if the input data has some problem that prevents
# computation (e.g. is empty).
#
# Expected signature:
#
#   da_function(data_dict, args...; kwargs...) -> (processed_data, ...)  or `nothing`
#
# where:
#
#   - data_dict::Dict
#   - processed_data::VecOrMat{<:Number}
#
####################################################################################################

"""
    daScatterGalaxy(
        data_dict::Dict,
        x_quantity::Symbol,
        y_quantity::Symbol;
        <keyword arguments>
    )::NTuple{2,Vector{<:Number}}

Compute two quantities for every cell/particle in the simulation using [`scatterQty`](@ref).

# Arguments

  - `data_dict::Dict`: Data dictionary (see [`makeDataDict`](@ref) for the canonical description).
  - `x_quantity::Symbol`: Quantity for the x axis. It can be any of the quantities valid for [`scatterQty`](@ref).
  - `y_quantity::Symbol`: Quantity for the y axis. It can be any of the quantities valid for [`scatterQty`](@ref).
  - `x_log::Union{Unitful.Units,Nothing}=nothing`: Target unit for `x_quantity`, if you want to apply ``\\log_{10}`` to the `x_quantity`. If set to `nothing`, the data from [`scatterQty`](@ref) is left as is.
  - `y_log::Union{Unitful.Units,Nothing}=nothing`: Target unit for `y_quantity`, if you want to apply ``\\log_{10}`` to the `y_quantity`. If set to `nothing`, the data from [`scatterQty`](@ref) is left as is.
  - `filter_function::Function=filterNothing`: Filter function to be applied to `data_dict` before any other computation. See the required signature and examples in `./src/analysis/filters.jl`.

# Returns

  - A tuple with two elements:

      + A vector with the values of `x_quantity`.
      + A vector with the values of `y_quantity`.
"""
function daScatterGalaxy(
    data_dict::Dict,
    x_quantity::Symbol,
    y_quantity::Symbol;
    x_log::Union{Unitful.Units,Nothing}=nothing,
    y_log::Union{Unitful.Units,Nothing}=nothing,
    filter_function::Function=filterNothing,
)::NTuple{2,Vector{<:Number}}

    filtered_dd = filterData(data_dict; filter_function)

    x_values = scatterQty(filtered_dd, x_quantity)
    y_values = scatterQty(filtered_dd, y_quantity)

    (
        length(x_values) == length(y_values) ||
        throw(ArgumentError("daScatterGalaxy: :$(x_quantity) and :$(y_quantity) have a different \
        number of elements. Maybe they are quantities for different types of cells/particles?"))
    )

    x_idxs = isnothing(x_log) ? Int64[] : map(iszero, x_values)
    y_idxs = isnothing(y_log) ? Int64[] : map(iszero, y_values)

    delete_idxs = x_idxs ∪ y_idxs

    deleteat!(x_values, delete_idxs)
    deleteat!(y_values, delete_idxs)

    if any(isempty, [x_values, y_values])

        logging[] && @warn("daScatterGalaxy: The results of scatterQty() are empty")

        return Float64[], Float64[]

    end

    x_axis = isnothing(x_log) ? x_values : log10.(ustrip.(x_log, x_values))
    y_axis = isnothing(y_log) ? y_values : log10.(ustrip.(y_log, y_values))

    return x_axis, y_axis

end

"""
    daIntegrateGalaxy(
        data_dict::Dict,
        x_quantity::Symbol,
        y_quantity::Symbol;
        <keyword arguments>
    )::NTuple{2,Vector{<:Number}}

Compute two global quantities for the simulation.

# Arguments

  - `data_dict::Dict`: Data dictionary (see [`makeDataDict`](@ref) for the canonical description).
  - `x_quantity::Symbol`: Quantity for the x axis. It can be any of the quantities valid for [`integrateQty`](@ref).
  - `y_quantity::Symbol`: Quantity for the y axis. It can be any of the quantities valid for [`integrateQty`](@ref).
  - `filter_function::Function=filterNothing`: Filter function to be applied to `data_dict` before any other computation. See the required signature and examples in `./src/analysis/filters.jl`.

# Returns

  - A tuple with two elements:

      + A single element vector with the value of `x_quantity`.
      + A single element vector with the value of `y_quantity`.
"""
function daIntegrateGalaxy(
    data_dict::Dict,
    x_quantity::Symbol,
    y_quantity::Symbol;
    filter_function::Function=filterNothing,
)::NTuple{2,Vector{<:Number}}

    filtered_dd = filterData(data_dict; filter_function)

    return [integrateQty(filtered_dd, x_quantity)], [integrateQty(filtered_dd, y_quantity)]

end

"""
    daScatterDensity(
        data_dict::Dict,
        x_quantity::Symbol,
        y_quantity::Symbol;
        <keyword arguments>
    )::Tuple{Vector{<:Number},Vector{<:Number},Matrix{Float64}}

Compute two quantities for every cell/particle using a 2D histogram.

# Arguments

  - `data_dict::Dict`: Data dictionary (see [`makeDataDict`](@ref) for the canonical description).
  - `x_quantity::Symbol`: Quantity for the x axis. It can be any of the quantities valid for [`scatterQty`](@ref).
  - `y_quantity::Symbol`: Quantity for the y axis. It can be any of the quantities valid for [`scatterQty`](@ref).
  - `x_range::Union{NTuple{2,<:Number},Nothing}=nothing`: x axis range for the histogram grid. If set to `nothing`, the extrema of the values will be used.
  - `y_range::Union{NTuple{2,<:Number},Nothing}=nothing`: y axis range for the histogram grid. If set to `nothing`, the extrema of the values will be used.
  - `x_log::Union{Unitful.Units,Nothing}=nothing`: Target unit for `x_quantity`, if you want to apply ``\\log_{10}`` to the `x_quantity`. If set to `nothing`, the data from [`scatterQty`](@ref) is left as is.
  - `y_log::Union{Unitful.Units,Nothing}=nothing`: Target unit for `y_quantity`, if you want to apply ``\\log_{10}`` to the `y_quantity`. If set to `nothing`, the data from [`scatterQty`](@ref) is left as is.
  - `n_bins::Int=100`: Number of bins per side of the grid.
  - `filter_function::Function=filterNothing`: Filter function to be applied to `data_dict` before any other computation. See the required signature and examples in `./src/analysis/filters.jl`.

# Returns

  - A tuple with three elements:

      + A vector with the x coordinates of the grid.
      + A vector with the y coordinates of the grid.
      + A matrix with the the ``\\log_{10}`` of the counts in each bin.
"""
function daScatterDensity(
    data_dict::Dict,
    x_quantity::Symbol,
    y_quantity::Symbol;
    x_range::Union{NTuple{2,<:Number},Nothing}=nothing,
    y_range::Union{NTuple{2,<:Number},Nothing}=nothing,
    x_log::Union{Unitful.Units,Nothing}=nothing,
    y_log::Union{Unitful.Units,Nothing}=nothing,
    n_bins::Int=100,
    filter_function::Function=filterNothing,
)::Tuple{Vector{<:Number},Vector{<:Number},Matrix{Float64}}

    x_values, y_values = daScatterGalaxy(
        data_dict,
        x_quantity,
        y_quantity;
        x_log,
        y_log,
        filter_function,
    )

    if any(isempty, [x_values, y_values])

        logging[] && @warn("daScatterDensity: The results of scatterQty() are empty")

        return collect(1:n_bins), collect(1:n_bins), fill(NaN, (n_bins, n_bins))

    end

    # If no range is specified, use the extrema of the values
    x_range = isnothing(x_range) ? extrema(x_values) : x_range
    y_range = isnothing(y_range) ? extrema(y_values) : y_range

    # Compute the half width of a bin
    x_bin_h_width = 0.5 * (x_range[2] - x_range[1]) / n_bins
    y_bin_h_width = 0.5 * (y_range[2] - y_range[1]) / n_bins

    # Compute the center position of each bin
    x_axis = collect(range(x_range[1] + x_bin_h_width; length=n_bins, step=2 * x_bin_h_width))
    y_axis = collect(range(y_range[1] + y_bin_h_width; length=n_bins, step=2 * y_bin_h_width))

    # Compute the 2D histogram
    counts = histogram2D(
        permutedims(hcat(x_values, y_values), (2, 1)),
        collect(range(x_range[1], x_range[2]; length=n_bins + 1)),
        collect(range(y_range[1], y_range[2]; length=n_bins + 1)),
    )

    # Apply log10 to enhance the contrast
    values = log10.(counts)

    if logging[]

        clean_vals = filter(!isnan, values)

        if isempty(clean_vals)

            min_max_vals = (NaN, NaN)
            mean_val     = NaN
            median_val   = NaN
            mode_val     = NaN

        else

            min_max_vals = extrema(clean_vals)
            mean_val     = mean(clean_vals)
            median_val   = median(clean_vals)
            mode_val     = mode(clean_vals)

        end

        @info(
            "\nlog₁₀(Counts) \
            \n  Simulation: $(basename(data_dict[:sim_data].path)) \
            \n  Snapshot:   $(data_dict[:snap_data].global_index) \
            \n  Min - Max:  $(min_max_vals) \
            \n  Mean:       $(mean_val) \
            \n  Median:     $(median_val) \
            \n  Mode:       $(mode_val)"
        )

    end

    # The transpose and reverse operation are used to conform to the way heatmap! expect the matrix
    # to be structured
    z_axis = reverse!(transpose(values), dims=2)

    return x_axis, y_axis, z_axis

end

"""
    daScatterWeightedDensity(
        data_dict::Dict,
        x_quantity::Symbol,
        y_quantity::Symbol,
        z_quantity::Symbol;
        <keyword arguments>
    )::Tuple{Vector{<:Number},Vector{<:Number},Matrix{Float64}}

Compute two quantities for every cell/particle using a 2D histogram, weighted by `z_quantity`.

# Arguments

  - `data_dict::Dict`: Data dictionary (see [`makeDataDict`](@ref) for the canonical description).
  - `x_quantity::Symbol`: Quantity for the x axis. It can be any of the quantities valid for [`scatterQty`](@ref).
  - `y_quantity::Symbol`: Quantity for the y axis. It can be any of the quantities valid for [`scatterQty`](@ref).
  - `z_quantity::Symbol`: Quantity for the weights. It can be any of the quantities valid for [`scatterQty`](@ref).
  - `x_range::Union{NTuple{2,<:Number},Nothing}=nothing`: x axis range. If set to `nothing`, the extrema of the values will be used.
  - `y_range::Union{NTuple{2,<:Number},Nothing}=nothing`: y axis range. If set to `nothing`, the extrema of the values will be used.
  - `x_log::Union{Unitful.Units,Nothing}=nothing`: Target unit for `x_quantity`, if you want to apply ``\\log_{10}`` to the `x_quantity`. If set to `nothing`, the data from [`scatterQty`](@ref) is left as is.
  - `y_log::Union{Unitful.Units,Nothing}=nothing`: Target unit for `y_quantity`, if you want to apply ``\\log_{10}`` to the `y_quantity`. If set to `nothing`, the data from [`scatterQty`](@ref) is left as is.
  - `total::Bool=true`: If the sum (default) or the mean of `z_quantity` will be used as the value of each bin.
  - `n_bins::Int=100`: Number of bins per side of the grid.
  - `filter_function::Function=filterNothing`: Filter function to be applied to `data_dict` before any other computation. See the required signature and examples in `./src/analysis/filters.jl`.

# Returns

  - A tuple with three elements:

      + A vector with the x coordinates of the grid.
      + A vector with the y coordinates of the grid.
      + A matrix with the ``\\log_{10}`` of `z_quantity` for each bin.
"""
function daScatterWeightedDensity(
    data_dict::Dict,
    x_quantity::Symbol,
    y_quantity::Symbol,
    z_quantity::Symbol;
    x_range::Union{NTuple{2,<:Number},Nothing}=nothing,
    y_range::Union{NTuple{2,<:Number},Nothing}=nothing,
    x_log::Union{Unitful.Units,Nothing}=nothing,
    y_log::Union{Unitful.Units,Nothing}=nothing,
    total::Bool=true,
    n_bins::Int=100,
    filter_function::Function=filterNothing,
)::Tuple{Vector{<:Number},Vector{<:Number},Matrix{Float64}}

    filtered_dd = filterData(data_dict; filter_function)

    x_values = scatterQty(filtered_dd, x_quantity)
    y_values = scatterQty(filtered_dd, y_quantity)
    z_values = scatterQty(filtered_dd, z_quantity)

    (
        allequal(length, [x_values, y_values, z_values]) ||
        throw(ArgumentError("daScatterGalaxy: :$(x_quantity), :$(y_quantity), and :$(z_quantity) \
        have different lengths. Maybe they are quantities for different types of cells/particles?"))
    )

    x_idxs = isnothing(x_log) ? Int64[] : map(iszero, x_values)
    y_idxs = isnothing(y_log) ? Int64[] : map(iszero, y_values)

    delete_idxs = x_idxs ∪ y_idxs

    deleteat!(x_values, delete_idxs)
    deleteat!(y_values, delete_idxs)
    deleteat!(z_values, delete_idxs)

    if any(isempty, [x_values, y_values, z_values])

        logging[] && @warn("daScatterWeightedDensity: The results of scatterQty() are empty")

        return collect(1:n_bins), collect(1:n_bins), fill(NaN, (n_bins, n_bins))

    end

    x_values = isnothing(x_log) ? x_values : log10.(ustrip.(x_log, x_values))
    y_values = isnothing(y_log) ? y_values : log10.(ustrip.(y_log, y_values))

    # If no range is specified, use the extrema of the values
    x_range = isnothing(x_range) ? extrema(x_values) : x_range
    y_range = isnothing(y_range) ? extrema(y_values) : y_range

    # Compute the half width of a bin
    x_bin_h_width = 0.5 * (x_range[2] - x_range[1]) / n_bins
    y_bin_h_width = 0.5 * (y_range[2] - y_range[1]) / n_bins

    # Compute the center position of each bin
    x_axis = collect(range(x_range[1] + x_bin_h_width; length=n_bins, step=2 * x_bin_h_width))
    y_axis = collect(range(y_range[1] + y_bin_h_width; length=n_bins, step=2 * y_bin_h_width))

    # Compute the 2D histogram
    histogram = histogram2D(
        permutedims(hcat(x_values, y_values), (2, 1)),
        z_values,
        collect(range(x_range[1], x_range[2]; length=n_bins + 1)),
        collect(range(y_range[1], y_range[2]; length=n_bins + 1));
        total,
    )

    z_unit = plotParams(z_quantity).unit

    # Apply log10 to enhance the contrast
    values = log10.(ustrip.(z_unit, histogram))

    if logging[]

        clean_vals = filter(!isnan, values)

        if isempty(clean_vals)

            min_max_vals = (NaN, NaN)
            mean_val     = NaN
            median_val   = NaN
            mode_val     = NaN

        else

            min_max_vals = extrema(clean_vals)
            mean_val     = mean(clean_vals)
            median_val   = median(clean_vals)
            mode_val     = mode(clean_vals)

        end

        @info(
            "\nlog₁₀($(z_quantity) [$(z_unit)]) \
            \n  Simulation: $(basename(filtered_dd[:sim_data].path)) \
            \n  Snapshot:   $(filtered_dd[:snap_data].global_index) \
            \n  Min - Max:  $(min_max_vals) \
            \n  Mean:       $(mean_val) \
            \n  Median:     $(median_val) \
            \n  Mode:       $(mode_val)"
        )

    end

    # The transpose and reverse operation are used to conform to the way heatmap! expect the matrix
    # to be structured
    z_axis = reverse!(transpose(values), dims=2)

    return x_axis, y_axis, z_axis

end

"""
    daProfile(
        data_dict::Dict,
        quantity::Symbol,
        grid::LinearGrid;
        <keyword arguments>
    )::Union{Tuple{Vector{<:Unitful.Length},Vector{<:Number}},Nothing}

Compute a profile.

# Arguments

  - `data_dict::Dict`: Data dictionary (see [`makeDataDict`](@ref) for the canonical description).
  - `quantity::Symbol`: Target quantity. It can be any of the quantities valid for [`scatterQty`](@ref).
  - `grid::LinearGrid`: Linear grid.
  - `norm::Union{Symbol,Nothing}=nothing`: The value of `quantity` in each bin will be divided by the corresponding value of `norm`. It can be any of the quantities valid for [`scatterQty`](@ref). If set to `nothing`, no operation is applied.
  - `y_log::Union{Unitful.Units,Nothing}=nothing`: Target unit for `quantity`, if you want to apply ``\\log_{10}`` to `quantity`. If set to `nothing`, the data is left as is.
  - `flat::Bool=true`: If the profile will be 2D (rings), or 3D (spherical shells).
  - `total::Bool=true`: If the sum (default) or the mean of `quantity` will be computed for each bin. This affects the values of `norm` too.
  - `cumulative::Bool=false`: If the profile will be accumulated (after dividing by `norm`).
  - `density::Bool=false`: If the profile will be of the density of `quantity` (after dividing by `norm`).
  - `filter_function::Function=filterNothing`: Filter function to be applied to `data_dict` before any other computation. See the required signature and examples in `./src/analysis/filters.jl`.

# Returns

  - A tuple with two elements:

      + A vector with the position of each ring or spherical shells.
      + A vector with the value of `quantity` in each each ring or spherical shells.

    It returns `nothing` if any of the necessary quantities are missing.
"""
function daProfile(
    data_dict::Dict,
    quantity::Symbol,
    grid::LinearGrid;
    norm::Union{Symbol,Nothing}=nothing,
    y_log::Union{Unitful.Units,Nothing}=nothing,
    flat::Bool=true,
    total::Bool=true,
    cumulative::Bool=false,
    density::Bool=false,
    filter_function::Function=filterNothing,
)::Union{Tuple{Vector{<:Unitful.Length},Vector{<:Number}},Nothing}

    filtered_dd = filterData(data_dict; filter_function)

    plot_params = plotParams(quantity)

    # Get the cell/particle type
    cp_type = plot_params.cp_type

    (
        isnothing(cp_type) &&
        throw(ArgumentError("daProfile: `quantity` = $(quantity) has cell/particle type `nothing`. \
        It is not valid to make a profile"))
    )

    # Read the positions and values
    positions = filtered_dd[cp_type]["POS "]
    values    = scatterQty(filtered_dd, quantity)
    n_values  = isnothing(norm) ? Number[] : scatterQty(filtered_dd, norm)

    if any(isempty, [values, positions])

        logging[] && @warn("daProfile: The results of scatterQty() and/or the positions are empty")

        return nothing

    end

    profile = computeProfile(
        positions,
        values,
        grid;
        norm=n_values,
        flat,
        total,
        cumulative,
        density,
    )

    x_axis = copy(grid.x_axis)

    if !isnothing(y_log)

        # Delete zeros
        idxs = map(iszero, profile)

        deleteat!(x_axis, idxs)
        deleteat!(profile, idxs)

        y_axis = log10.(ustrip.(y_log, profile))

    else

        y_axis = profile

    end

    return x_axis, y_axis

end

"""
    daBandProfile(
        data_dict::Dict,
        quantity::Symbol,
        grid::LinearGrid;
        <keyword arguments>
    )::Union{
        Tuple{Vector{<:Unitful.Length},Vector{<:Number},Vector{<:Number},Vector{<:Number}},
        Tuple{Vector{<:Unitful.Length},Vector{<:Number},Vector{<:Number}},
        Nothing,
    }

Compute a profile with uncertainty bands.

# Arguments

  - `data_dict::Dict`: Data dictionary (see [`makeDataDict`](@ref) for the canonical description).
  - `quantity::Symbol`: Target quantity. It can be any of the quantities valid for [`scatterQty`](@ref).
  - `grid::LinearGrid`: Linear grid.
  - `flat::Bool=true`: If the profile will be 2D (rings), or 3D (spherical shells).
  - `density::Bool=false`: If the profile will be of the density of `quantity`.
  - `ylog::Bool=false`: If true, returns the profile of ``\\log_{10}``(`quantity`).
  - `error_bar::Bool=false`: If the returned values will be compatible with `band!` (default) or with `errorbars!`.
  - `filter_function::Function=filterNothing`: Filter function to be applied to `data_dict` before any other computation. See the required signature and examples in `./src/analysis/filters.jl`.

# Returns

  - A tuple with three elements if `error_bar` = false:

      + A vector with the position of each ring or spherical shells.
      + A vector with the "low" value of `quantity` in each each ring or spherical shells.
      + A vector with the "high" value of `quantity` in each each ring or spherical shells.

  - Or a tuple with four elements if `error_bar` = true:

      + A vector with the position of each ring or spherical shells.
      + A vector with the central value of `quantity` in each each ring or spherical shells.
      + A vector with the "low" error of `quantity` in each each ring or spherical shells.
      + A vector with the "high" error of `quantity` in each each ring or spherical shells.

    It returns `nothing` if any of the necessary quantities are missing.
"""
function daBandProfile(
    data_dict::Dict,
    quantity::Symbol,
    grid::LinearGrid;
    flat::Bool=true,
    density::Bool=false,
    ylog::Bool=false,
    error_bar::Bool=false,
    filter_function::Function=filterNothing,
)::Union{
    Tuple{Vector{<:Unitful.Length},Vector{<:Number},Vector{<:Number},Vector{<:Number}},
    Tuple{Vector{<:Unitful.Length},Vector{<:Number},Vector{<:Number}},
    Nothing,
}

    filtered_dd = filterData(data_dict; filter_function)

    plot_params = plotParams(quantity)

    # Get the cell/particle type
    type = plot_params.cp_type

    # Read the positions and values
    positions = filtered_dd[type]["POS "]
    values    = scatterQty(filtered_dd, quantity)

    if any(isempty, [values, positions])
        (
            logging[] &&
            @warn("daBandProfile: The results of scatterQty() and/or the positions are empty")
        )
        return nothing
    end

    if ylog
        idxs = map(iszero, values)
        deleteat!(positions, idxs)
        deleteat!(values, idxs)

        if any(isempty, [values, positions])
            (
                logging[] &&
                @warn("daBandProfile: After removing zeros, the results of scatterQty() are empty")
            )
            return nothing
        end

        values = log10.(ustrip.(plot_params.unit, values))
    end

    center, low, high = computeBandProfile(positions, values, grid; flat, density)

    if error_bar
        return grid.x_axis, center, center .- low, high .- center
    end

    return grid.x_axis, low, high

end

"""
    daStellarHistory(
        data_dict::Dict;
        <keyword arguments>
    )::Union{Tuple{Vector{<:Unitful.Time},Vector{<:Number}},Nothing}

Compute the evolution of a given stellar `quantity` using the stellar ages at a given instant in time.

# Arguments

  - `data_dict::Dict`: Data dictionary (see [`makeDataDict`](@ref) for the canonical description).
  - `quantity::Symbol=:sfr`: Target quantity. The options are:

      + `:sfr`                 -> Star formation rate.
      + `:ssfr`                -> The specific star formation rate.
      + `:stellar_mass`        -> Stellar mass.
      + `:stellar_metallicity` -> Mass fraction of all elements above He in the stars (solar units).
  - `n_bins::Int=100`: Number of bins (time intervals).
  - `ylog::Bool=false`: If true, returns ``\\log_{10}``(`quantity`).
  - `filter_function::Function=filterNothing`: Filter function to be applied to `data_dict` before any other computation. See the required signature and examples in `./src/analysis/filters.jl`.

# Returns

  - A tuple with two elements:

      + A vector with the physical times.
      + A vector with the values of `quantity` at each time.
"""
function daStellarHistory(
    data_dict::Dict;
    quantity::Symbol=:sfr,
    n_bins::Int=100,
    ylog::Bool=false,
    filter_function::Function=filterNothing,
)::Union{Tuple{Vector{<:Unitful.Time},Vector{<:Number}},Nothing}

    filtered_dd = filterData(data_dict; filter_function)

    birth_ticks = filtered_dd[:stellar]["GAGE"]
    masses      = filtered_dd[:stellar]["MASS"]

    # Return `nothing` if there are less than 3 stars
    if any(x -> length(x) <= 2, [birth_ticks, masses])
        (
            logging[] &&
            @warn("daStellarHistory: The are less than 3 stars in snapshot \
            $(filtered_dd[:snap_data].path)")
        )
        return nothing
    end

    # Compute the stellar birth times
    if filtered_dd[:sim_data].cosmological
        # Go from scale factor to physical time
        birth_times = computeTime(birth_ticks, filtered_dd[:snap_data].header)
    else
        birth_times = birth_ticks
    end

    # Compute the birth time range
    min, max = extrema(birth_times)

    grid = LinearGrid(min, max, n_bins)

    # Compute the total stellar mass in each time bin
    stellar_masses = histogram1D(birth_times, masses, grid; empty_nan=false)

    # Compute the time axis
    bin_width = (max - min) / n_bins
    x_axis = collect(range(min + (bin_width * 0.5), length=n_bins, step=bin_width))

    # Compute the stellar quantity
    if quantity == :sfr

        y_axis = stellar_masses ./ bin_width

    elseif quantity == :ssfr

        accu_mass = cumsum(stellar_masses)
        y_axis = (stellar_masses ./ bin_width) ./ accu_mass

    elseif quantity == :stellar_mass

        y_axis = cumsum(stellar_masses)

    elseif quantity == :stellar_metallicity

        metal_mass = computeMass(data_dict, :Z_stellar)

        if isempty(metal_mass)

            logging[] && @warn("daStellarHistory: The metal mass data is missing")

            return nothing

        end

        stellar_metallicities = histogram1D(birth_times, metal_mass, grid; empty_nan=false)

        Z = cumsum(stellar_metallicities) ./ cumsum(stellar_masses)
        y_axis = Z ./ SOLAR_METALLICITY

    else

        throw(ArgumentError("daStellarHistory: `quantity` can only be :sfr, :ssfr, \
        :stellar_mass or :stellar_metallicity, but I got :$(quantity)"))

    end

    if ylog

        idxs = map(iszero, y_axis)
        deleteat!(x_axis, idxs)
        deleteat!(y_axis, idxs)

        y_axis = log10.(ustrip.(plotParams(quantity).unit, y_axis))

    end

    return x_axis, y_axis

end

"""
    daHistogram(
        data_dict::Dict,
        quantity::Symbol,
        grid::LinearGrid;
        <keyword arguments>
    )::Union{Tuple{Vector{<:Number},Vector{<:Number}},Nothing}

Compute a 1D histogram of a given `quantity`.

# Arguments

  - `data_dict::Dict`: Data dictionary (see [`makeDataDict`](@ref) for the canonical description).
  - `quantity::Symbol`: Target quantity. It can be any of the quantities valid for [`scatterQty`](@ref).
  - `grid::LinearGrid`: Linear grid.
  -  `log::Bool=false`: If the histogram bins will be logarithmic.
  - `norm::Int=0`: Number of counts that will be use to normalize the histogram. If left as 0, the histogram will be normalize with the maximum bin count.
  - `filter_function::Function=filterNothing`: Filter function to be applied to `data_dict` before any other computation. See the required signature and examples in `./src/analysis/filters.jl`.

# Returns

  - A tuple with two elements:

      + A vector with the value of `quantity` corresponding to each bin.
      + A vector with the counts.
"""
function daHistogram(
    data_dict::Dict,
    quantity::Symbol,
    grid::LinearGrid;
    log::Bool=false,
    norm::Int=0,
    filter_function::Function=filterNothing,
)::Union{Tuple{Vector{<:Number},Vector{Float64}},Nothing}

    scatter_qty = scatterQty(data_dict, quantity)

    if isempty(scatter_qty)
        (
            logging[] &&
            @warn("daHistogram: The results of scatterQty() for :$(quantity) are empty")
        )
        return nothing
    end

    plot_params = plotParams(quantity)

    # Get the cell/particle type
    cp_type = plot_params.cp_type

    (
        isnothing(cp_type) &&
        throw(ArgumentError("daHistogram: `quantity` = $(quantity) has cell/particle type \
        `nothing`. It is not valid to make a histogram"))
    )

    # Filter after computing the values, to preserve quantities that depends on
    # global properties (e.g. global gravitational potential)
    idxs = filter_function(data_dict)[cp_type]

    values = scatter_qty[idxs]

    if isempty(values)
        (
            logging[] &&
            @warn("daHistogram: After filtering, there is no data left for :$(quantity)")
        )
        return nothing
    end

    # Compute the histogram
    counts = histogram1D(values, grid; empty_nan=false)

    # Normalize the counts
    y_axis = isPositive(norm) ? counts ./ norm : counts ./ maximum(counts)

    if logging[]

        clean_vals = filter(x -> !isnan(x) && !isinf(x), values)

        if isempty(clean_vals)

            min_max_v = (NaN, NaN)
            mean_v    = NaN
            median_v  = NaN
            mode_v    = NaN

        else

            min_max_v = extrema(clean_vals)
            mean_v    = mean(clean_vals)
            median_v  = median(clean_vals)
            mode_v    = mode(clean_vals)

        end

        @info(
            "\nCounts \
            \n  Simulation:  $(basename(data_dict[:sim_data].path)) \
            \n  Snapshot:    $(data_dict[:snap_data].global_index) \
            \n  Quantity:    $(quantity) \
            \n  Type:        $(cp_type) \
            \n  Largest bin: $(grid.x_axis[argmax(counts)]) \
            \n  Max count:   $(maximum(counts)) \
            \n  Min - Max:   $(min_max_v) \
            \n  Mean:        $(mean_v) \
            \n  Median:      $(median_v) \
            \n  Mode:        $(mode_v)"
        )

    end

    if log
        x_axis = log10.(ustrip.(plot_params.unit, grid.x_axis))
    else
        x_axis = copy(grid.x_axis)
    end

    return x_axis, y_axis

end

"""
    daHistogram(
        data_dict::Dict,
        quantity::Symbol,
        n_bins::Int;
        <keyword arguments>
    )::Union{Tuple{Vector{<:Number},Vector{<:Number}},Nothing}

Compute a 1D histogram of a given `quantity`.

!!! note

    This method uses the extrema of the values of `quantity` as the range of the histogram.

# Arguments

  - `data_dict::Dict`: Data dictionary (see [`makeDataDict`](@ref) for the canonical description).
  - `quantity::Symbol`: Target quantity. It can be any of the quantities valid for [`scatterQty`](@ref).
  - `n_bins::Int`: Number of bins.
  - `log::Bool=false`: If the histogram bins will be logarithmic.
  - `norm::Int=0`: Number of count that will be use to normalize the histogram. If left as 0, the histogram will be normalize with the maximum bin count.
  - `filter_function::Function=filterNothing`: Filter function to be applied to `data_dict` before any other computation. See the required signature and examples in `./src/analysis/filters.jl`.

# Returns

  - A tuple with two elements:

      + A vector with the value of `quantity` corresponding to each bin.
      + A vector with the counts.
"""
function daHistogram(
    data_dict::Dict,
    quantity::Symbol,
    n_bins::Int;
    log::Bool=false,
    norm::Int=0,
    filter_function::Function=filterNothing,
)::Union{Tuple{Vector{<:Number},Vector{Float64}},Nothing}

    scatter_qty = scatterQty(data_dict, quantity)

    if isempty(scatter_qty)
        (
            logging[] &&
            @warn("daHistogram: The results of scatterQty() for :$(quantity) are empty")
        )
        return nothing
    end

    plot_params = plotParams(quantity)

    # Get the cell/particle type
    cp_type = plot_params.cp_type

    (
        isnothing(cp_type) &&
        throw(ArgumentError("daHistogram: `quantity` = $(quantity) has cell/particle type \
        `nothing`. It is not valid to make a histogram"))
    )

    # Filter after computing the values, to preserve quantities that depends on
    # global properties (e.g. global gravitational potential)
    idxs = filter_function(data_dict)[cp_type]

    values = scatter_qty[idxs]

    # Remove non-positive values if log scale is requested
    log && filter!(isPositive, values)

    if isempty(values)
        (
            logging[] &&
            @warn("daHistogram: After filtering, there is no data left for :$(quantity)")
        )
        return nothing
    end

    grid = LinearGrid(extrema(values)..., n_bins; log)

    # Compute the histogram
    counts = histogram1D(values, grid; empty_nan=false)

    # Normalize the counts
    y_axis = isPositive(norm) ? counts ./ norm : counts ./ maximum(counts)

    if logging[]

        clean_vals = filter(x -> !isnan(x) && !isinf(x), values)

        if isempty(clean_vals)

            min_max_v = (NaN, NaN)
            mean_v    = NaN
            median_v  = NaN
            mode_v    = NaN

        else

            min_max_v = extrema(clean_vals)
            mean_v    = mean(clean_vals)
            median_v  = median(clean_vals)
            mode_v    = mode(clean_vals)

        end

        @info(
            "\nCounts \
            \n  Simulation:  $(basename(data_dict[:sim_data].path)) \
            \n  Snapshot:    $(data_dict[:snap_data].global_index) \
            \n  Quantity:    $(quantity) \
            \n  Type:        $(cp_type) \
            \n  Largest bin: $(grid.x_axis[argmax(counts)]) \
            \n  Max count:   $(maximum(counts)) \
            \n  Min - Max:   $(min_max_v) \
            \n  Mean:        $(mean_v) \
            \n  Median:      $(median_v) \
            \n  Mode:        $(mode_v)"
        )

    end

    if log
        x_axis = log10.(ustrip.(plot_params.unit, grid.x_axis))
    else
        x_axis = copy(grid.x_axis)
    end

    return x_axis, y_axis

end

"""
    daRotationCurve(
        data_dict::Dict;
        <keyword arguments>
    )::Tuple{Vector{<:Unitful.Length},Vector{<:Unitful.Velocity}}

Compute a rotation curve.

# Arguments

  - `data_dict::Dict`: Data dictionary (see [`makeDataDict`](@ref) for the canonical description).
  - `R::Unitful.Length=DISK_R`: Maximum radial distance for the rotation curve.
  - `filter_function::Function=filterNothing`: Filter function to be applied to `data_dict` before any other computation. See the required signature and examples in `./src/analysis/filters.jl`.

# Returns

  - A tuple with two elements:

      + A vector with the distances to each star.
      + A vector with the circular velocity of each star.
"""
function daRotationCurve(
    data_dict::Dict;
    R::Unitful.Length=DISK_R,
    filter_function::Function=filterNothing,
)::Tuple{Vector{<:Unitful.Length},Vector{<:Unitful.Velocity}}

    filtered_dd = filterData(data_dict; filter_function)

    # Compute the circular velocities and the radial distances of each star
    r, vcirc = computeVcirc(filtered_dd, :stellar)

    # Only leave the data within a sphere of radius `R`
    rangeCut!(r, vcirc, (0.0u"kpc", R))

    # Sort the arrays radially
    idx = sortperm(r)

    return r[idx], vcirc[idx]

end

"""
    daBarGasFractions(
        data_dict::Dict,
        quantity::Symbol,
        grid::LinearGrid;
        <keyword arguments>
    )::Union{NTuple{2,Vector{<:Number}},Nothing}

Compute a bar plot, where the bins are for a given gas `quantity` and the heights are the corresponding gas fractions.

# Arguments

  - `data_dict::Dict`: Data dictionary (see [`makeDataDict`](@ref) for the canonical description).
  - `quantity::Symbol`: Target quantity. It can be any of the quantities valid for [`scatterQty`](@ref) with the cell/particle type `:gas`.
  - `grid::LinearGrid`: Linear grid.
  - `components::Vector{Symbol}=[:ode_ionized, :ode_atomic, :ode_cold]`: List of gas components to be considered. The fractions will be normalized to this list of components. See [`COMPONENTS`](@ref) for options.
  - `filter_function::Function=filterNothing`: Filter function to be applied to `data_dict` before any other computation. See the required signature and examples in `./src/analysis/filters.jl`.

# Returns

  - A tuple with the two elements required by `barplot!` from [Makie](https://docs.makie.org/stable/):

      + A vector with the index of each bar.
      + A vector with the height of each bar.
"""
function daBarGasFractions(
    data_dict::Dict,
    quantity::Symbol,
    grid::LinearGrid;
    components::Vector{Symbol}=[:ode_ionized, :ode_atomic, :ode_cold],
    filter_function::Function=filterNothing,
)::Union{NTuple{2,Vector{<:Number}},Nothing}

    qty_type = plotParams(quantity).cp_type

    (
        qty_type == :gas ||
        throw(ArgumentError("daBarGasFractions: `quantity` = $(quantity) should have the \
        cell/particle type :gas, but I got :$(qty_type)"))
    )

    comp_types = [plotParams(Symbol(component, :_mass)).cp_type for component in components]

    (
        all(x->x == :gas, comp_types) ||
        throw(ArgumentError("daBarGasFractions: Every component in $(components) should have the \
        cell/particle type :gas, but I got the types :$(comp_types)"))
    )

    filtered_dd = filterData(data_dict; filter_function)
    gas_qty = scatterQty(filtered_dd, quantity)

    if isempty(gas_qty)
        (
            logging[] &&
            @warn("daBarGasFractions: The results of scatterQty() for :$(quantity) are empty")
        )
        return nothing
    end

    # Compute the mass of each gas component for every cell
    masses = [computeMass(filtered_dd, component) for component in components]

    if any(isempty, masses)

        logging[] && @warn("daBarGasFractions: The data for the mass of the components is missing")

        return nothing

    end

    # Compute the `quantity` histogram for each gas component
    histograms = [histogram1D(gas_qty, mass, grid; empty_nan=false) for mass in masses]

    # Compute the number of bins
    n_bins = grid.n_bins

    # Compute the number of bars per bin
    n_bars = length(components)

    percents = Vector{Float64}(undef, n_bins * n_bars)

    for i in 1:n_bins

        bin_masses = getindex.(histograms, i)
        total_mass = sum(bin_masses)

        if iszero(total_mass)

            percents[((i - 1) * n_bars + 1):(i * n_bars)] .= 0.0

        else

            # Compute the normalized mass fraction of each component inside the current bin
            fractions = uconvert.(Unitful.NoUnits, bin_masses ./ total_mass)

            percents[((i - 1) * n_bars + 1):(i * n_bars)] .= fractions .* 100

        end

    end

    return repeat(1:n_bins, inner=n_bars), percents

end

"""
    daMolla2015(
        data_dict::Dict,
        grid::LinearGrid,
        quantity::Symbol;
        <keyword arguments>
    )::Union{
        Tuple{
            Vector{<:Unitful.Length},
            <:Union{Vector{<:SurfaceDensity},Vector{<:MassFlowDensity},Vector{Float64}}
        },
        Nothing,
    }

Compute a profile for the Milky Way, compatible with the experimental data in Mollá et al. (2015).

# Arguments

  - `data_dict::Dict`: Data dictionary (see [`makeDataDict`](@ref) for the canonical description).
  - `grid::LinearGrid`: Linear grid.
  - `quantity::Symbol`: Quantity for the y axis. The options are:

      + `:stellar_area_density`               -> Stellar mass surface density.
      + `:sfr_area_density`                   -> Star formation rate surface density.
      + `:ode_molecular_stellar_area_density` -> Molecular mass surface density.
      + `:br_molecular_area_density`          -> Molecular mass surface density, computed using the pressure relation in Blitz et al. (2006).
      + `:ode_atomic_area_density`            -> Atomic mass surface density.
      + `:br_atomic_area_density`             -> Atomic mass surface density, computed using the pressure relation in Blitz et al. (2006).
      + `:O_stellar_abundance`                -> Stellar abundance of oxygen, as ``12 + \\log_{10}(\\mathrm{O \\, / \\, H})``.
      + `:N_stellar_abundance`                -> Stellar abundance of nitrogen, as ``12 + \\log_{10}(\\mathrm{N \\, / \\, H})``.
      + `:C_stellar_abundance`                -> Stellar abundance of carbon, as ``12 + \\log_{10}(\\mathrm{C \\, / \\, H})``.
  - `y_unit::Unitful.Units=Unitful.NoUnits`: Unit for `quantity`.
  - `filter_function::Function=filterNothing`: Filter function to be applied to `data_dict` before any other computation. See the required signature and examples in `./src/analysis/filters.jl`.

# Returns

  - A tuple with two elements:

      + A vector with the position of each ring.
      + A vector with the value `quantity` in each ring.

    It returns `nothing` if any of the necessary quantities are missing.

# References

L. Blitz et al. (2006). *The Role of Pressure in GMC Formation II: The H2-Pressure Relation*. The Astrophysical Journal, **650(2)**, 933. [doi:10.1086/505417](https://doi.org/10.1086/505417)

M. Mollá et al. (2015). *Galactic chemical evolution: stellar yields and the initial mass function*. Monthly Notices of the Royal Astronomical Society **451(4)**, 3693–3708. [doi:10.1093/mnras/stv1102](https://doi.org/10.1093/mnras/stv1102)
"""
function daMolla2015(
    data_dict::Dict,
    grid::LinearGrid,
    quantity::Symbol;
    y_unit::Unitful.Units=Unitful.NoUnits,
    filter_function::Function=filterNothing,
)::Union{
    Tuple{
        Vector{<:Unitful.Length},
        <:Union{Vector{<:SurfaceDensity},Vector{<:MassFlowDensity},Vector{Float64}},
    },
    Nothing,
}

    filtered_dd = filterData(data_dict; filter_function)

    if quantity ∈ STELLAR_ABUNDANCE

        element   = Symbol(replace(string(quantity), "_stellar_abundance" => ""))

        positions = filtered_dd[:stellar]["POS "]
        masses    = computeElementMass(filtered_dd, :stellar, element) ./ ATOMIC_WEIGHTS[element]
        norm      = computeElementMass(filtered_dd, :stellar, :H) ./ ATOMIC_WEIGHTS[:H]
        scaling   = x -> ABUNDANCE_SHIFT[element] .+ log10.(ustrip.(y_unit, x))
        density   = false

    else

        plot_params = plotParams(quantity)
        component   = Symbol(replace(string(quantity), "_area_density" => ""))

        if component == :sfr
            masses = scatterQty(filtered_dd, :observational_sfr)
        else
            masses = computeMass(filtered_dd, component)
        end

        positions = filtered_dd[plot_params.cp_type]["POS "]
        norm      = Number[]
        scaling   = x -> log10.(ustrip.(y_unit, x))
        density   = true

    end

    if any(isempty, [positions, masses])

        logging[] && @warn("daMolla2015: The data for the mass and/or the positions is missing")

        return nothing

    end

    density_profile = scaling(computeProfile(positions, masses, grid; norm, density))

    return grid.x_axis, density_profile

end

"""
    daDensity2DProjection(
        data_dict::Dict,
        grid::CubicGrid,
        component::Symbol,
        field_type::Symbol;
        <keyword arguments>
    )::Tuple{Vector{<:Unitful.Length},Vector{<:Unitful.Length},VecOrMat{Float64}}

Project a 3D mass density field into a given plane.

!!! note

    If the source of the field are particles, a simple 3D histogram is used. If they are Voronoi cells, the density of the cells that cross the line of sight of each pixel are added up.

# Arguments

  - `data_dict::Dict`: Data dictionary (see [`makeDataDict`](@ref) for the canonical description).
  - `grid::CubicGrid`: Cubic grid.
  - `component::Symbol`: Target component. It can only be one of the elements of [`COMPONENTS`](@ref).
  - `field_type::Symbol`: If the field is made up of `:particles` or Voronoi `:cells`.
  - `projection_plane::Symbol=:xy`: Projection plane. The options are `:xy`, `:xz`, and `:yz`.
  - `reduce_grid::Symbol=:square`: Type of 2D grid to do the final projection. The options are:

      + `:square`    -> The density distribution will be projected into a regular square grid, with a resolution `reduce_factor` times lower than `grid`. This emulates the way the surface densities are measured in observations. `reduce_factor` = 1 means no reduction in resolution.
      + `:circular` -> The density distribution will be projected into a flat circular grid, formed by a series of `reduce_factor` concentric rings. This emulates the traditional way the Kennicutt-Schmidt law is measured in some simulations. `reduce_factor` = 1 means that the result will be a single point. Note that this behaves the opposite way than the `reduce_grid` = :square case.
  - `reduce_factor::Int=1`: Factor by which the resolution of the result will be reduced. This will be applied after the density projection. If `reduce_grid` = :square, the new values will be computed averaging the values of neighboring pixels. `reduce_factor` has to divide the size of `grid` exactly. If `reduce_grid` = :circular, the new values will be computed averaging the values of the pixels the fall within each of the `reduce_factor` concentric rings.
  - `mmap_path::String="./"`: Path to store the memory-mapped file if needed (for matrices larger than [`MMAP_THRESHOLD`](@ref)).
  - `m_unit::Unitful.Units=u"Msun"`: Mass unit.
  - `l_unit::Unitful.Units=u"pc"`: Length unit.
  - `filter_function::Function=filterNothing`: Filter function to be applied to `data_dict` before any other computation. See the required signature and examples in `./src/analysis/filters.jl`.

# Returns

  - A tuple with three elements:

      + A vector with the x coordinates of the grid.
      + A vector with the y coordinates of the grid.
      + A matrix or vector with the ``\\log_{10}`` of the density at each bin of the 2D grid.
"""
function daDensity2DProjection(
    data_dict::Dict,
    grid::CubicGrid,
    component::Symbol,
    field_type::Symbol;
    projection_plane::Symbol=:xy,
    reduce_grid::Symbol=:square,
    reduce_factor::Int=1,
    mmap_path::String="./",
    m_unit::Unitful.Units=u"Msun",
    l_unit::Unitful.Units=u"pc",
    filter_function::Function=filterNothing,
)::Tuple{Vector{<:Unitful.Length},Vector{<:Unitful.Length},VecOrMat{Float64}}

    quantity = Symbol(component, :_mass)

    mass_grid = quantity3DProjection(
        data_dict,
        grid,
        quantity,
        field_type;
        empty_nan=false,
        log=false,
        mmap_path,
        filter_function,
    )

    # Project `mass_grid` to the target plane
    if projection_plane == :xy
        dims = 3
    elseif projection_plane == :xz
        dims = 2
    elseif projection_plane == :yz
        dims = 1
    else
        throw(ArgumentError("daDensity2DProjection: The argument `projection_plane` must be \
        :xy, :xz or :yz, but I got :$(projection_plane)"))
    end

    # For comological simulations with comoving units, correct
    # the density so it is always in physical units
    if !PHYSICAL_UNITS && data_dict[:sim_data].cosmological
        # Correction factor for the area
        # A [physical units] = A [comoving units] * a0^2
        physical_factor = data_dict[:snap_data].scale_factor^2
    else
        physical_factor = 1.0
    end

    munit_factor = ustrip(m_unit, 1.0 * plotParams(quantity).unit)
    voxel_area   = ustrip(l_unit^2, grid.bin_size_2D[dims]) * physical_factor

    density = dropdims(sum(mass_grid; dims); dims) .* (munit_factor / voxel_area)

    if reduce_grid == :square

        # Reduce the resolution of the result using a new square grid
        # `reduce_factor` here is the factor by wich the number of rows and columns will be reduced
        density = reduceMatrix(density, reduce_factor)
        x_axis  = reduceTicks(grid.x_axis, reduce_factor)
        y_axis  = reduceTicks(grid.y_axis, reduce_factor)

    elseif reduce_grid == :circular

        # Reduce the resolution of the result using a circular grid
        # `reduce_factor` here is the number of bins for the circular grid
        density = projectIntoLinearGrid(density, reduce_factor)
        x_axis  = [grid.size[1] * (2 * i - 1) / (4 * reduce_factor) for i in 1:reduce_factor]
        y_axis  = x_axis

    else

        throw(ArgumentError("daDensity2DProjection: `reduce_grid` can only be :square or \
        :circular, but I got :$(reduce_grid)"))

    end

    # Set bins with a value of 0 to NaN
    replace!(x -> iszero(x) ? NaN : x, density)

    # Apply log10 to enhance the contrast
    z_axis = log10.(density)

    if logging[]

        log_z_axis = filter(!isnan, z_axis)

        if isempty(log_z_axis)

            min_max_z = (NaN, NaN)
            mean_z    = NaN
            median_z  = NaN
            mode_z    = NaN

        else

            min_max_z = extrema(log_z_axis)
            mean_z    = mean(log_z_axis)
            median_z  = median(log_z_axis)
            mode_z    = mode(log_z_axis)

        end

        @info(
            "\nMass density range - log₁₀(Σ [$(m_unit * l_unit^-2)]) \
            \n  Simulation: $(basename(filtered_dd[:sim_data].path)) \
            \n  Snapshot:   $(filtered_dd[:snap_data].global_index) \
            \n  Component:  $(component) \
            \n  Field type: $(field_type) \
            \n  Plane:      $(projection_plane) \
            \n  Min - Max:  $(min_max_z) \
            \n  Mean:       $(mean_z) \
            \n  Median:     $(median_z) \
            \n  Mode:       $(mode_z)"
        )

    end

    return x_axis, y_axis, z_axis

end

"""
    daGasSFR2DProjection(
        data_dict::Dict,
        grid::CubicGrid,
        field_type::Symbol;
        <keyword arguments>
    )::Tuple{Vector{<:Unitful.Length},Vector{<:Unitful.Length},VecOrMat{Float64}}

Project the 3D gas SFR field into a given plane.

!!! note

    If the source of the field are particles, a simple 3D histogram is used. If they are Voronoi cells instead, the SFR of the cells that cross the line of sight of each pixel are added up.

# Arguments

  - `data_dict::Dict`: Data dictionary (see [`makeDataDict`](@ref) for the canonical description).
  - `grid::CubicGrid`: Cubic grid.
  - `field_type::Symbol`: If the field is made up of `:particles` or Voronoi `:cells`.
  - `projection_plane::Symbol=:xy`: Projection plane. The options are `:xy`, `:xz`, and `:yz`.
  - `reduce_grid::Symbol=:square`: Type of 2D grid to do the final projection. The options are:

      + `:square`    -> The density distribution will be projected into a regular square grid, with a resolution `reduce_factor` times lower than `grid`. This emulates the way the surface densities are measured in observations. `reduce_factor` = 1 means no reduction in resolution.
      + `:circular` -> The density distribution will be projected into a flat circular grid, formed by a series of `reduce_factor` concentric rings. This emulates the traditional way the Kennicutt-Schmidt law is measured in some simulations. `reduce_factor` = 1 means that the result will be a single point. Note that this behaves the opposite way than the `reduce_grid` = :square case.
  - `reduce_factor::Int=1`: Factor by which the resolution of the result will be reduced. This will be applied after the density projection. If `reduce_grid` = :square, the new values will be computed averaging the values of neighboring pixels. `reduce_factor` has to divide the size of `grid` exactly. If `reduce_grid` = :circular, the new values will be computed averaging the values of the pixels the fall within each of the `reduce_factor` concentric rings.
  - `mmap_path::String="./"`: Path to store the memory-mapped file if needed (for matrices larger than [`MMAP_THRESHOLD`](@ref)).
  - `m_unit::Unitful.Units=u"Msun"`: Mass unit.
  - `t_unit::Unitful.Units=u"yr"`: Time unit.
  - `filter_function::Function=filterNothing`: Filter function to be applied to `data_dict` before any other computation. See the required signature and examples in `./src/analysis/filters.jl`.

# Returns

  - A tuple with three elements:

      + A vector with the x coordinates of the grid.
      + A vector with the y coordinates of the grid.
      + A matrix or vector with the ``\\log_{10}`` of the gas SFR at each bin of the 2D grid.
"""
function daGasSFR2DProjection(
    data_dict::Dict,
    grid::CubicGrid,
    field_type::Symbol;
    projection_plane::Symbol=:xy,
    reduce_grid::Symbol=:square,
    reduce_factor::Int=1,
    mmap_path::String="./",
    m_unit::Unitful.Units=u"Msun",
    t_unit::Unitful.Units=u"yr",
    filter_function::Function=filterNothing,
)::Tuple{Vector{<:Unitful.Length},Vector{<:Unitful.Length},VecOrMat{Float64}}

    sfr_grid = quantity3DProjection(
        data_dict,
        grid,
        :gas_sfr,
        field_type;
        empty_nan=false,
        log=false,
        mmap_path,
        filter_function,
    )

    # Project `sfr_grid` to the target plane
    if projection_plane == :xy
        dims = 3
    elseif projection_plane == :xz
        dims = 2
    elseif projection_plane == :yz
        dims = 1
    else
        throw(ArgumentError("daGasSFR2DProjection: The argument `projection_plane` must be \
        :xy, :xz or :yz, but I got :$(projection_plane)"))
    end

    sfr_factor = ustrip(m_unit * t_unit^-1, 1.0 * plotParams(:gas_sfr).unit)

    sfr = dropdims(sum(sfr_grid; dims); dims) .* sfr_factor

    if reduce_grid == :square

        # Reduce the resolution of the result using a new square grid
        # `reduce_factor` here is the factor by wich the number of rows and columns will be reduced
        sfr    = reduceMatrix(sfr, reduce_factor; total=true)
        x_axis = reduceTicks(grid.x_axis, reduce_factor)
        y_axis = reduceTicks(grid.y_axis, reduce_factor)

    elseif reduce_grid == :circular

        # Reduce the resolution of the result using a circular grid
        # `reduce_factor` here is the number of bins for the circular grid
        sfr    = projectIntoLinearGrid(sfr, reduce_factor; total=true)
        x_axis = [grid.size[1] * (2 * i - 1) / (4 * reduce_factor) for i in 1:reduce_factor]
        y_axis = x_axis

    else

        throw(ArgumentError("daGasSFR2DProjection: `reduce_grid` can only be :square or :circular, \
        but I got :$(reduce_grid)"))

    end

    # Set bins with 0 or Inf to NaN
    replace!(x -> (iszero(x) || isinf(x)) ? NaN : x, sfr)

    # Apply log10 to enhance the contrast
    z_axis = log10.(sfr)

    if logging[]

        log_z_axis = filter(!isnan, z_axis)

        if isempty(log_z_axis)

            min_max_z = (NaN, NaN)
            mean_z    = NaN
            median_z  = NaN
            mode_z    = NaN

        else

            min_max_z = extrema(log_z_axis)
            mean_z    = mean(log_z_axis)
            median_z  = median(log_z_axis)
            mode_z    = mode(log_z_axis)

        end

        # Print the SFR range
        @info(
            "\nGas SFR range - log₁₀(SFR [$(m_unit * t_unit^-1)]) \
            \n  Simulation: $(basename(filtered_dd[:sim_data].path)) \
            \n  Snapshot:   $(filtered_dd[:snap_data].global_index) \
            \n  Field type: $(field_type) \
            \n  Plane:      $(projection_plane) \
            \n  Min - Max:  $(min_max_z) \
            \n  Mean:       $(mean_z) \
            \n  Median:     $(median_z) \
            \n  Mode:       $(mode_z)"
        )

    end

    return x_axis, y_axis, z_axis

end

"""
    daMetallicity2DProjection(
        data_dict::Dict,
        grid::CubicGrid,
        component::Symbol,
        field_type::Symbol;
        <keyword arguments>
    )::Tuple{Vector{<:Unitful.Length},Vector{<:Unitful.Length},VecOrMat{Float64}}

Project the 3D metallicity field to a given plane.

!!! note

    If `element` = :all, the total metallicity (in solar units) is computed. If `element` = :X, the abundance of element X is computed, [`ABUNDANCE_SHIFT`](@ref) + ``\\log_{10}(X/H). In both cases the total value of the column given by the line of sight of each pixel is computed.

# Arguments

  - `data_dict::Dict`: Data dictionary (see [`makeDataDict`](@ref) for the canonical description).
  - `grid::CubicGrid`: Cubic grid.
  - `component::Symbol`: Target component. It can be either `:stellar` or `:gas`.
  - `field_type::Symbol`: If the field is made up of `:particles` or Voronoi `:cells`.
  - `element::Symbol=:all`: Target element. The possibilities are the keys of [`ELEMENT_INDEX`](@ref). Set it to :all if you want the total metallicity.
  - `projection_plane::Symbol=:xy`: Projection plane. The options are `:xy`, `:xz`, and `:yz`.
  - `reduce_grid::Symbol=:square`: Type of 2D grid to do the final projection. The options are:

      + `:square`    -> The density distribution will be projected into a regular square grid, with a resolution `reduce_factor` times lower than `grid`. This emulates the way the surface densities are measured in observations. `reduce_factor` = 1 means no reduction in resolution.
      + `:circular` -> The density distribution will be projected into a flat circular grid, formed by a series of `reduce_factor` concentric rings. This emulates the traditional way the Kennicutt-Schmidt law is measured in some simulations. `reduce_factor` = 1 means that the result will be a single point. Note that this behaves the opposite way than the `reduce_grid` = :square case.
  - `reduce_factor::Int=1`: Factor by which the resolution of the result will be reduced. This will be applied after the density projection. If `reduce_grid` = :square, the new values will be computed averaging the values of neighboring pixels. `reduce_factor` has to divide the size of `grid` exactly. If `reduce_grid` = :circular, the new values will be computed averaging the values of the pixels the fall within each of the `reduce_factor` concentric rings.
  - `mmap_path::String="./"`: Path to store the memory-mapped file if needed (for matrices larger than [`MMAP_THRESHOLD`](@ref)).
  - `filter_function::Function=filterNothing`: Filter function to be applied to `data_dict` before any other computation. See the required signature and examples in `./src/analysis/filters.jl`.

# Returns

  - A tuple with three elements:

      + A vector with the x coordinates of the grid.
      + A vector with the y coordinates of the grid.
      + A matrix or vector with the ``\\log_{10}`` of the metallicity at each bin of the 2D grid.
"""
function daMetallicity2DProjection(
    data_dict::Dict,
    grid::CubicGrid,
    component::Symbol,
    field_type::Symbol;
    element::Symbol=:all,
    projection_plane::Symbol=:xy,
    reduce_grid::Symbol=:square,
    reduce_factor::Int=1,
    mmap_path::String="./",
    filter_function::Function=filterNothing,
)::Tuple{Vector{<:Unitful.Length},Vector{<:Unitful.Length},VecOrMat{Float64}}

    if element == :all

        metal_qty = Symbol(:Z_, component, :_mass)

        norm_qty = Symbol(component, :_mass)

    elseif haskey(ELEMENT_INDEX, element)

        metal_qty = Symbol(element, :_, component, :_mass)

        norm_qty = Symbol(:H_, component, :_mass)

    else

        throw(ArgumentError("daMetallicity2DProjection: The argument `element` can only be :all \
        or one of the keys of `ELEMENT_INDEX` (see `./src/constants/globals.jl`), \
        but I got :$(element)"))

    end

    metal_grid, nn_idxs = quantity3DProjection(
        data_dict,
        grid,
        metal_qty,
        field_type;
        empty_nan=false,
        log=false,
        mmap_path,
        return_idxs=true,
        filter_function,
    )

    norm_grid = quantity3DProjection(
        data_dict,
        grid,
        norm_qty,
        nn_idxs;
        empty_nan=false,
        log=false,
        return_idxs=false,
        filter_function,
    )

    # Project `metal_grid` and `norm_grid` to the target plane
    if projection_plane == :xy
        dims = 3
    elseif projection_plane == :xz
        dims = 2
    elseif projection_plane == :yz
        dims = 1
    else
        throw(ArgumentError("daMetallicity2DProjection: The argument `projection_plane` must \
        be :xy, :xz or :yz, but I got :$(projection_plane)"))
    end


    metal_mass = dropdims(sum(metal_grid; dims); dims)
    norm_mass  = dropdims(sum(norm_grid; dims); dims)

    if reduce_grid == :square

        # Reduce the resolution of the result into a new square grid
        # `reduce_factor` here is the factor by wich the number of rows and columns will be reduced
        metal_mass = reduceMatrix(metal_mass, reduce_factor; total=true)
        norm_mass  = reduceMatrix(norm_mass, reduce_factor; total=true)

        x_axis = reduceTicks(grid.x_axis, reduce_factor)
        y_axis = reduceTicks(grid.y_axis, reduce_factor)

        # Compute the metallicity in each bin
        metallicity = metal_mass ./ norm_mass

    elseif reduce_grid == :circular

        # Reduce the resolution of the result into a circular grid
        # `reduce_factor` here is the number of bins for the circular grid
        metal_mass = projectIntoLinearGrid(metal_mass, reduce_factor; total=true)
        norm_mass  = projectIntoLinearGrid(norm_mass, reduce_factor; total=true)

        x_axis  = [grid.size[1] * (2 * i - 1) / (4 * reduce_factor) for i in 1:reduce_factor]
        y_axis = x_axis

        # Compute the metallicity in each bin
        metallicity = metal_mass ./ norm_mass

    else

        throw(ArgumentError("daMetallicity2DProjection: `reduce_grid` can only be :square or \
        :circular, but I got :$(reduce_grid)"))

    end

    # Set bins with 0 or Inf to NaN
    replace!(x -> (iszero(x) || isinf(x)) ? NaN : x, metallicity)

    # Apply log10 to enhance the contrast
    if element == :all
        z_axis = log10.(metallicity ./ SOLAR_METALLICITY)
    else
        z_axis = ABUNDANCE_SHIFT[element] .+ log10.(metallicity)
    end

    if logging[]

        clean_z_axis = filter(!isnan, z_axis)

        if isempty(clean_z_axis)

            min_max_z = (NaN, NaN)
            mean_z    = NaN
            median_z  = NaN
            mode_z    = NaN

        else

            min_max_z = extrema(clean_z_axis)
            mean_z    = mean(clean_z_axis)
            median_z  = median(clean_z_axis)
            mode_z    = mode(clean_z_axis)

        end

        title = if element == :all
            "log₁₀(Z [Z⊙])"
        else
            if iszero(ABUNDANCE_SHIFT[element])
                "log₁₀($(element) / H)"
            else
                "$(ABUNDANCE_SHIFT[element]) + log₁₀($(element) / H)"
            end
        end

        @info(
            "\nMetallicity range - $(title) \
            \n  Simulation: $(basename(filtered_dd[:sim_data].path)) \
            \n  Snapshot:   $(filtered_dd[:snap_data].global_index) \
            \n  Component:  $(component) \
            \n  Field type: $(field_type) \
            \n  Plane:      $(projection_plane)\
            \n  Min - Max:  $(min_max_z) \
            \n  Mean:       $(mean_z) \
            \n  Median:     $(median_z) \
            \n  Mode:       $(mode_z)"
        )

    end

    return x_axis, y_axis, z_axis

end

"""
    daVelocityField(
        data_dict::Dict,
        grid::SquareGrid,
        component::Symbol;
        <keyword arguments>
    )::Tuple{Vector{<:Unitful.Length},Vector{<:Unitful.Length},Matrix{Float64},Matrix{Float64}}

Project a 3D velocity field into a given plane.

!!! note

    This function computes the mean velocity in each direction of the target plane, at each pixel of the grid, using an 2D histogram. For low resolution grids, the result its the same even for a field made of cells, while being much facter and simpler that the full Voronoi rasterisation.

# Arguments

  - `data_dict::Dict`: Data dictionary (see [`makeDataDict`](@ref) for the canonical description).
  - `grid::SquareGrid`: Square grid.
  - `component::Symbol`: Target component. It can only be one of the elements of [`COMPONENTS`](@ref).
  - `projection_plane::Symbol=:xy`: Projection plane. The options are `:xy`, `:xz`, and `:yz`.
  - `v_unit::Unitful.Units=u"km * s^-1"`: Velocity unit
  - `filter_function::Function=filterNothing`: Filter function to be applied to `data_dict` before any other computation. See the required signature and examples in `./src/analysis/filters.jl`.

# Returns

  - A tuple with four elements:

      + A vector with the x coordinates of the grid.
      + A vector with the y coordinates of the grid.
      + A matrix with the mean velocity in the x direction at each grid point.
      + A matrix with the mean velocity in the y direction at each grid point.
"""
function daVelocityField(
    data_dict::Dict,
    grid::SquareGrid,
    component::Symbol;
    projection_plane::Symbol=:xy,
    v_unit::Unitful.Units=u"km * s^-1",
    filter_function::Function=filterNothing,
)::Tuple{Vector{<:Unitful.Length},Vector{<:Unitful.Length},Matrix{Float64},Matrix{Float64}}

    filtered_dd = filterData(data_dict; filter_function)

    cp_type = plotParams(Symbol(component, :_mass)).cp_type

    positions  = filtered_dd[cp_type]["POS "]
    velocities = filtered_dd[cp_type]["VEL "]

    if any(isempty, [positions, velocities])
        (
            logging[] &&
            @warn("daVelocityField: The data for the positions and/or the velocities is missing")
        )
        return grid.x_axis, grid.y_axis, zeros(grid.n_bins[1]), zeros(grid.n_bins[2])
    end

    # Project the cells/particles to the chosen plane
    if projection_plane == :xy

        idxs = [1, 2]

    elseif projection_plane == :xz

        idxs = [1, 3]

    elseif projection_plane == :yz

        idxs = [2, 3]

    else

        throw(ArgumentError("daVelocityField: The argument `projection_plane` must be \
        :xy, :xz or :yz, but I got :$(projection_plane)"))

    end

    pos_2D = positions[idxs, :]

    # Compute the mean velocity field in the target plane
    vx = ustrip.(v_unit, histogram2D(pos_2D, vec(velocities[idxs[1], :]), grid; total=false))
    vy = ustrip.(v_unit, histogram2D(pos_2D, vec(velocities[idxs[2], :]), grid; total=false))

    return grid.x_axis, grid.y_axis, vx, vy

end

"""
    daSDSSMockup(
        data_dict::Dict,
        grid::CubicGrid;
        <keyword arguments>
    )::AbstractArray

Compute a mockup image emulating an SDSS observation.

# Arguments

  - `data_dict::Dict`: Data dictionary (see [`makeDataDict`](@ref) for the canonical description).
  - `grid::CubicGrid`: Cubic grid.
  - `projection_plane::Symbol=:xy`: Projection plane. The options are `:xy`, `:xz`, and `:yz`.
  - `smooth::Bool=false`: If gaussian smooththing will be applied to the whole image.
  - `mmap_path::String="./"`: Path to store the memory-mapped file if needed (for matrices larger than [`MMAP_THRESHOLD`](@ref)).
  - `filter_function::Function=filterNothing`: Filter function to be applied to `data_dict` before any other computation. See the required signature and examples in `./src/analysis/filters.jl`.

# Returns

  - An RGB image with the mock observation.
"""
function daSDSSMockup(
    data_dict::Dict,
    grid::CubicGrid;
    projection_plane::Symbol=:xy,
    smooth::Bool=false,
    mmap_path::String="./",
    filter_function::Function=filterNothing,
)::AbstractArray

    filtered_dd = filterData(data_dict; filter_function)

    # Get necessary quantities
    positions     = filtered_dd[:stellar]["POS "]
    masses        = ustrip.(u"Msun", filtered_dd[:stellar]["MASS"])
    metallicities = scatterQty(filtered_dd, :stellar_metallicity)
    log_ages      = log10.(ustrip.(u"yr", scatterQty(filtered_dd, :stellar_age)))

    # Compute the neutral mass grid
    Mn_grid = quantity3DProjection(
        filtered_dd,
        grid,
        :ode_neutral_mass,
        :cells;
        empty_nan=false,
        log=false,
        mmap_path,
    )

    # Compute the 3D stellar index
    stellar_3D_index = findIn3DGrid(positions, grid; cartesian=true)

    # Compute the interpolation functions for the magnitudes in the g, r, and i filters of SDSS
    g_interp, r_interp, i_interp = griMagnitudeInterpolation(MILLANIRIGOYEN2025_DATA_PATH)

    if projection_plane == :xy
        extended_dims = [1, 2]
        compact_dim = 3
    elseif projection_plane == :xz
        extended_dims = [1, 3]
        compact_dim  = 2
    elseif projection_plane == :yz
        extended_dims = [2, 3]
        compact_dim  = 1
    else
        throw(ArgumentError("daSDSSMockup: The argument `projection_plane` must be \
        :xy, :xz or :yz, but I got :$(projection_plane)"))
    end

    # For comological simulations with comoving units, correct
    # the density so it is always in physical units
    if !PHYSICAL_UNITS && filtered_dd[:sim_data].cosmological
        # Correction factor for the area
        # A [physical units] = A [comoving units] * a0^2
        physical_factor = filtered_dd[:snap_data].scale_factor^2
    else
        physical_factor = 1.0
    end

    m_unit     = plotParams(:ode_neutral_mass).unit
    ext_factor = 1.0 / ustrip(m_unit*u"pc^-2", EXTINCTION_FACTOR)
    voxel_area = ustrip(u"pc^2", grid.bin_size_2D[compact_dim]) * physical_factor

    # Allocate the memory to store the fluxes in each pixel
    Fg = zeros(Float64, grid.n_bins[extended_dims])
    Fr = zeros(Float64, grid.n_bins[extended_dims])
    Fi = zeros(Float64, grid.n_bins[extended_dims])

    for star_idx in eachindex(masses)

        stellar_idx = stellar_3D_index[:, star_idx]

        any(iszero, stellar_idx) && continue

        metallicity = metallicities[star_idx]
        log_age     = log_ages[star_idx]
        mass        = masses[star_idx]

        # Interpolate the stellar magnitudes
        Mg = g_interp(log_age, metallicity)
        Mr = r_interp(log_age, metallicity)
        Mi = i_interp(log_age, metallicity)

        # Select the line of sight of the current star
        idx = ntuple(i -> i == compact_dim ? (1:stellar_idx[compact_dim]) : stellar_idx[i], 3)

        # Compute the column density of neutral gas in front of the star
        Mn = sum(@view Mn_grid[idx...]) / voxel_area

        # Compute the extinction in the V band
        AV = Mn * ext_factor

        # Compute the extinction in the g, r, and i bands
        Ag = AV * AλAV_g
        Ar = AV * AλAV_r
        Ai = AV * AλAV_i

        stellar_2D_idx = stellar_idx[extended_dims]

        # Compute the stellar flux
        # See https://en.wikipedia.org/wiki/Apparent_magnitude
        # The reference flux should be the same for each band in the SDSS AB system
        # so we don't need to normalize, as the relative value of fluxes is all we need
        Fg[stellar_2D_idx...] += mass * 10.0^(-0.4 * (Mg + Ag))
        Fr[stellar_2D_idx...] += mass * 10.0^(-0.4 * (Mr + Ar))
        Fi[stellar_2D_idx...] += mass * 10.0^(-0.4 * (Mi + Ai))

    end

    # AstroImage() expect the matrix to be in "image" format
    Fg = reverse(Fg'; dims=2)
    Fr = reverse(Fr'; dims=2)
    Fi = reverse(Fi'; dims=2)

    # Turn the fluxes into RGB images
    red   = AstroImage(Fi)
    green = AstroImage(Fr)
    blue  = AstroImage(Fg)

    rgb = composecolors(
        [red, green, blue],
        ["red", "green", "blue"],
        stretch=asinhstretch,
        multiplier=[1, 1.7, 1]
    )

    if smooth
        rgb = imfilter(rgb, Kernel.gaussian(1200.0 / grid.n_bins))
    end

    return rgb

end

####################################################################################################
# Signature for data analysis functions targeting plotTimeSeries() in ./src/plotting/pipelines.jl
####################################################################################################
#
# A data analysis functions targeting plotTimeSeries() must take a Simulation struct (see
# ./src/constants/globals.jl for the canonical description), and return two vectors.
# It should return `nothing` if the input data has some problem that prevents computation
# (e.g. is empty).
#
# Expected signature:
#
#   da_function(sim_data, args...; kw_args...) -> (processed_data_x, processed_data_y) or `nothing`
#
# where:
#
#   - sim_data::Simulation
#   - processed_data_x::Vector{<:Number}
#   - processed_data_y::Vector{<:Number}
#
####################################################################################################

"""
    daEvolution(
        sim_data::Simulation,
        x_quantity::Symbol,
        y_quantity::Symbol;
        <keyword arguments>
    )::NTuple{2,Vector{<:Number}}

Compute the time series of two quantities, using [`integrateQty`](@ref) to compute their values at each time.

!!! note

    The log10 operation, if requested, is applied after the integration of the quantities.

# Arguments

  - `sim_data::Simulation`: The [`Simulation`](@ref) struct for the target simulation.
  - `x_quantity::Symbol`: Quantity for the x axis.
  - `y_quantity::Symbol`: Quantity for the y axis.
  - `trans_mode::Union{Symbol,Tuple{TranslationType,RotationType,Dict{Symbol,Vector{String}}}}=:all_box`: How to translate and rotate the cells/particles, before filtering with `filter_mode`. For options see [`selectTransformation`](@ref).
  - `filter_mode::Union{Symbol,Tuple{Function,Dict{Symbol,Vector{String}}}}=:all`: Which cells/particles will be selected. For options see [`selectFilter`](@ref).
  - `extra_filter::Function=filterNothing`: Filter function to be applied after `trans_mode` and `filter_mode` are applied. See the required signature and examples in `./src/analysis/filters.jl`.
  - `ff_request::Dict{Symbol,Vector{String}}=Dict{Symbol,Vector{String}}()`: Request dictionary for `extra_filter`.
  - `x_log::Union{Unitful.Units,Nothing}=nothing`: Target unit for `x_quantity`, if you want to apply ``\\log_{10}`` to the `x_quantity`. If set to `nothing`, the data from [`scatterQty`](@ref) is left as is.
  - `y_log::Union{Unitful.Units,Nothing}=nothing`: Target unit for `y_quantity`, if you want to apply ``\\log_{10}`` to the `y_quantity`. If set to `nothing`, the data from [`scatterQty`](@ref) is left as is.
  - `smooth::Int=0`: The result of [`integrateQty`](@ref) will be smoothed out using `smooth` bins. Set it to 0 if you want no smoothing.
  - `cumulative::Bool=false`: If the `y_quantity` will be accumulated or not.
  - `show_progress::Bool=true`: If a progress bar will be shown.

# Returns

  - A Tuple with two elements:

      + A Vector with the time series of `x_quantity`.
      + A Vector with the time series of `y_quantity`.
"""
function daEvolution(
    sim_data::Simulation,
    x_quantity::Symbol,
    y_quantity::Symbol;
    trans_mode::Union{Symbol,Tuple{TranslationType,RotationType,Dict{Symbol,Vector{String}}}}=:all_box,
    filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all,
    extra_filter::Function=filterNothing,
    ff_request::Dict{Symbol,Vector{String}}=Dict{Symbol,Vector{String}}(),
    x_log::Union{Unitful.Units,Nothing}=nothing,
    y_log::Union{Unitful.Units,Nothing}=nothing,
    smooth::Int=0,
    cumulative::Bool=false,
    show_progress::Bool=true,
)::NTuple{2,Vector{<:Number}}

    qty_request = mergeRequests(plotParams(x_quantity).request, plotParams(y_quantity).request)

    integration_functions = (dd->integrateQty(dd, x_quantity), dd->integrateQty(dd, y_quantity))

    return daEvolution(
        sim_data,
        qty_request,
        integration_functions;
        trans_mode,
        filter_mode,
        extra_filter,
        ff_request,
        x_log,
        y_log,
        smooth,
        cumulative,
        show_progress,
    )

end

"""
    daEvolution(
        sim_data::Simulation,
        qty_request::Dict{Symbol,Vector{String}},
        integration_functions::NTuple{2,Function};
        <keyword arguments>
    )::NTuple{2,Vector{<:Number}}

Compute the time series of two quantities, using the provided integration functions to compute their values at each time.

!!! note

    The log10 operation, if requested, is applied after the integration of the quantities.

# Arguments

  - `sim_data::Simulation`: The [`Simulation`](@ref) struct for the target simulation.
  - `qty_request::Dict{Symbol,Vector{String}}`: Request dictionary for both quantities.
  - `integration_functions::NTuple{2,Function}`: Functions to compute the integral value of the x and y quantities at a given time. The functions must have the signature:

    `integration_functions(data_dict::Dict)::Number`

    where

      + `data_dict::Dict`: Data dictionary (see [`makeDataDict`](@ref) for the canonical description).
  - `trans_mode::Union{Symbol,Tuple{TranslationType,RotationType,Dict{Symbol,Vector{String}}}}=:all_box`: How to translate and rotate the cells/particles, before filtering with `filter_mode`. For options see [`selectTransformation`](@ref).
  - `filter_mode::Union{Symbol,Tuple{Function,Dict{Symbol,Vector{String}}}}=:all`: Which cells/particles will be selected. For options see [`selectFilter`](@ref).
  - `extra_filter::Function=filterNothing`: Filter function to be applied after `trans_mode` and `filter_mode` are applied. See the required signature and examples in `./src/analysis/filters.jl`.
  - `ff_request::Dict{Symbol,Vector{String}}=Dict{Symbol,Vector{String}}()`: Request dictionary for `extra_filter`.
  - `x_log::Union{Unitful.Units,Nothing}=nothing`: Target unit for `x_quantity`, if you want to apply ``\\log_{10}`` to the `x_quantity`. If set to `nothing`, the data from [`scatterQty`](@ref) is left as is.
  - `y_log::Union{Unitful.Units,Nothing}=nothing`: Target unit for `y_quantity`, if you want to apply ``\\log_{10}`` to the `y_quantity`. If set to `nothing`, the data from [`scatterQty`](@ref) is left as is.
  - `smooth::Int=0`: The result of `integration_functions` will be smoothed out using `smooth` bins. Set it to 0 if you want no smoothing.
  - `cumulative::Bool=false`: If the `y_quantity` will be accumulated or not.
  - `show_progress::Bool=true`: If a progress bar will be shown.

# Returns

  - A Tuple with two elements:

      + A Vector with the time series of `x_quantity`.
      + A Vector with the time series of `y_quantity`.
"""
function daEvolution(
    sim_data::Simulation,
    qty_request::Dict{Symbol,Vector{String}},
    integration_functions::NTuple{2,Function};
    trans_mode::Union{Symbol,Tuple{TranslationType,RotationType,Dict{Symbol,Vector{String}}}}=:all_box,
    filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all,
    extra_filter::Function=filterNothing,
    ff_request::Dict{Symbol,Vector{String}}=Dict{Symbol,Vector{String}}(),
    x_log::Union{Unitful.Units,Nothing}=nothing,
    y_log::Union{Unitful.Units,Nothing}=nothing,
    smooth::Int=0,
    cumulative::Bool=false,
    show_progress::Bool=true,
)::NTuple{2,Vector{<:Number}}

    base_request = mergeRequests(qty_request, ff_request)

    translation, rotation, trans_request = selectTransformation(trans_mode, base_request)
    filter_function, request = selectFilter(filter_mode, trans_request)

    # Iterate over each snapshot in the slice
    iterator = eachrow(DataFrame(sim_data.snapshot_table[sim_data.slice, :]))

    n_frames = length(iterator)

    x_values = Vector{Number}(fill(NaN, n_frames))
    y_values = Vector{Number}(fill(NaN, n_frames))

    # Initialize the progress bar
    prog_bar = Progress(
        n_frames,
        dt=0.5,
        desc="Analyzing and plotting the data... ",
        color=:blue,
        barglyphs=BarGlyphs("|#  |"),
        enabled=show_progress,
    )

    for (i, sim_table_row) in pairs(iterator)

        # Skip missing snapshots
        if ismissing(sim_table_row[:snapshot_paths])

            logging[] && @warn("daEvolution: The snapshot $(sim_table_row[:row_id]) is missing")

            next!(prog_bar)

            continue

        end

        data_dict = makeDataDict(
            sim_data.path,
            sim_table_row[:row_id],
            request,
            sim_data.snapshot_table,
        )

        # Translate the data
        translateData!(data_dict, translation...)

        # Rotate the data
        rotateData!(data_dict, rotation...)

        # Filter the data
        filterData!(data_dict; filter_function)

        # Filter the data again
        filterData!(data_dict; filter_function=extra_filter)

        # Compute the value for the x axis
        x_values[i] = integration_functions[1](data_dict)

        # Compute the value for the y axis
        y_values[i] = integration_functions[2](data_dict)

        # Move the progress bar forward
        next!(prog_bar)

    end

    # Sort by the x axis
    sort_idxs = sortperm(x_values)

    x_values = x_values[sort_idxs]
    y_values = y_values[sort_idxs]

    if cumulative
        y_values = cumsum(y_values)
    end

    # Inf values break the plot auto limits computation, so we remove them
    if isnothing(x_log)
        x_filter = isinf
    else
        x_filter = x -> isinf(x) || !isPositive(x)
    end

    if isnothing(y_log)
        y_filter = isinf
    else
        y_filter = y -> isinf(y) || !isPositive(y)
    end

    x_idxs = map(x_filter, x_values)
    y_idxs = map(y_filter, y_values)

    delete_idxs = x_idxs ∪ y_idxs

    deleteat!(x_values, delete_idxs)
    deleteat!(y_values, delete_idxs)

    if any(isempty, [x_values, y_values])

        logging[] && @warn("daEvolution: The results of `integration_functions` are empty")

        return Float64[], Float64[]

    end

    x_axis = isnothing(x_log) ? x_values : log10.(ustrip.(x_log, x_values))
    y_axis = isnothing(y_log) ? y_values : log10.(ustrip.(y_log, y_values))

    if iszero(smooth)
        return x_axis, y_axis
    else
        return smoothWindow(x_axis, y_axis, smooth)
    end

end

"""
    daSFRtxt(
        sim_data::Simulation,
        x_quantity::Symbol,
        y_quantity::Symbol;
        <keyword arguments>
    )::NTuple{2,Vector{<:Number}}

Compute the stellar mass or SFR evolution using the data in the `sfr.txt` file.

# Arguments

  - `sim_data::Simulation`: The [`Simulation`](@ref) struct for the target simulation.
  - `x_quantity::Symbol`: Quantity for the x axis. The options are:

      + `:scale_factor`  -> Scale factor.
      + `:redshift`      -> Redshift.
      + `:physical_time` -> Physical time since the Big Bang.
      + `:lookback_time` -> Physical time left to reach the last snapshot.
  - `y_quantity::Symbol`: Quantity for the y axis. The options are:

      + `:stellar_mass` -> Cumulative stellar mass.
      + `:sfr`          -> Star formation rate.
  - `smooth::Int=0`: The result will be smoothed out using `smooth` bins. Set it to 0 if you want no smoothing.

# Returns

  - A Tuple with two elements:

      + A Vector with the time series of `x_quantity`.
      + A Vector with the time series of `y_quantity`.
"""
function daSFRtxt(
    sim_data::Simulation,
    x_quantity::Symbol,
    y_quantity::Symbol;
    smooth::Int=0,
)::NTuple{2,Vector{<:Number}}

    snapshot_paths = sim_data.snapshot_table[!, :snapshot_paths]

    idx = findfirst(!ismissing, snapshot_paths)

    (
        !isnothing(idx) ||
        throw(ArgumentError("daSFRtxt: I couldn't find any snapshots in $(sim_data.path), \
        and I need at least one for unit conversion"))
    )

    # Read the first snapshot
    snapshot_path = snapshot_paths[idx]

    # Read its header
    header = readSnapHeader(snapshot_path)

    # Read the data in the `sfr.txt` file
    sfr_txt_data = readSfrFile(joinpath(sim_data.path, SFR_REL_PATH), snapshot_path)

    # Read the time axis
    time_ticks = sfr_txt_data[!, 1]

    if x_quantity == :scale_factor

        (
            !sim_data.cosmological &&
            logging[] &&
            @warn("daSFRtxt: For non-cosmological simulations `x_quantity` can only be \
            :physical_time")
        )

        x_axis = time_ticks

    elseif x_quantity == :redshift

        if sim_data.cosmological

            x_axis = (1.0 ./ time_ticks) .- 1.0

        else

            (
                logging[] &&
                @warn("daSFRtxt: For non-cosmological simulations `x_quantity` can only be \
                :physical_time")
            )

            x_axis = time_ticks

        end

    elseif x_quantity == :physical_time

        if sim_data.cosmological
            x_axis = computeTime(time_ticks, header)
        else
            x_axis = time_ticks
        end

    elseif x_quantity == :lookback_time

        if sim_data.cosmological
            physical_times = computeTime(time_ticks, header)
        else
            physical_times = time_ticks
        end

        x_axis = maximum(physical_times) .- physical_times

    else

        throw(ArgumentError("daSFRtxt: `x_quantity` can only be :scale_factor, :redshift, \
        :physical_time or :lookback_time, but I got :$(x_quantity)"))

    end

    if y_quantity == :stellar_mass

        y_axis = sfr_txt_data[!, 6]

    elseif y_quantity == :sfr

        y_axis = sfr_txt_data[!, 3]

    else

        throw(ArgumentError("daSFRtxt: `y_quantity` can only be :stellar_mass or :sfr, \
        but I got :$(y_quantity)"))

    end

    # Apply smoothing if required
    if !iszero(smooth)
        x_axis, y_axis = smoothWindow(x_axis, y_axis, smooth)
    end

    return x_axis, y_axis

end

"""
    daCPUtxt(
        sim_data::Simulation,
        target::String,
        x_quantity::Symbol,
        y_quantity::Symbol;
        <keyword arguments>
    )::NTuple{2,Vector{<:Number}}

Compute the evolution of a measured quantity in the `cpu.txt` file, for a given `target` process.

# Arguments

  - `sim_data::Simulation`: The [`Simulation`](@ref) struct for the target simulation.
  - `target::String`: Target process.
  - `x_quantity::Symbol`: Quantity for the x axis. The options are:

      + `:time_step`              -> Time step.
      + `:physical_time`          -> Physical time since the Big Bang.
      + `:clock_time_s`           -> Clock time duration of the time step in seconds.
      + `:clock_time_percent`     -> Clock time duration of the time step as a percentage.
      + `:tot_clock_time_s`       -> Total clock time in seconds.
      + `:tot_clock_time_percent` -> Total clock time as a percentage.
  - `y_quantity::Symbol`: Quantity for the y axis. The options are:

      + `:time_step`              -> Time step.
      + `:physical_time`          -> Physical time since the Big Bang.
      + `:clock_time_s`           -> Clock time duration of the time step in seconds.
      + `:clock_time_percent`     -> Clock time duration of the time step as a percentage.
      + `:tot_clock_time_s`       -> Total clock time in seconds.
      + `:tot_clock_time_percent` -> Total clock time as a percentage.
  - `y_log::Union{Unitful.Units,Nothing}=nothing`: Target unit for `y_quantity`, if you want to apply ``\\log_{10}`` to the `y_quantity`. If set to `nothing`, the data from `cpu.txt` is left as is.
  - `smooth::Int=0`: The result will be smoothed out using `smooth` bins. Set it to 0 if you want no smoothing.

# Returns

  - A Tuple with two elements:

      + A Vector with the time series of `x_quantity`.
      + A Vector with the time series of `y_quantity`.
"""
function daCPUtxt(
    sim_data::Simulation,
    target::String,
    x_quantity::Symbol,
    y_quantity::Symbol;
    y_log::Union{Unitful.Units,Nothing}=nothing,
    smooth::Int=0,
)::NTuple{2,Vector{<:Number}}

    snapshot_paths = sim_data.snapshot_table[!, :snapshot_paths]

    idx = findfirst(!ismissing, snapshot_paths)

    (
        !isnothing(idx) ||
        throw(ArgumentError("daCPUtxt: I couldn't find any snapshots in $(sim_data.path), \
        and I need at least one for unit conversion"))
    )

    # Read the first snapshot
    snapshot_path = snapshot_paths[idx]

    # Read its header
    header = readSnapHeader(snapshot_path)

    # Read the data in the `sfr.txt` file
    cpu_txt_data = readCpuFile(joinpath(sim_data.path, CPU_REL_PATH), [target])[target]

    if x_quantity == :time_step

        x_axis = cpu_txt_data[:, 1]

    elseif x_quantity == :physical_time

        if sim_data.cosmological
            x_axis = computeTime(cpu_txt_data[:, 2], header)
        else
            x_axis = cpu_txt_data[:, 2] .* internalUnits("CLKT", snapshot_path)
        end

    elseif x_quantity == :clock_time_s

        x_axis = cpu_txt_data[:, 3] .* u"s"

    elseif x_quantity == :clock_time_percent

        x_axis = cpu_txt_data[:, 4]

    elseif x_quantity == :tot_clock_time_s

        x_axis = cpu_txt_data[:, 5] .* u"s"

    elseif x_quantity == :tot_clock_time_percent

        x_axis = cpu_txt_data[:, 6]

    else

        throw(ArgumentError("daCPUtxt: I don't recognize the x_quantity :$(x_quantity)"))

    end

    if y_quantity == :time_step

        y_values = cpu_txt_data[:, 1]

    elseif y_quantity == :physical_time

        if sim_data.cosmological
            y_values = computeTime(cpu_txt_data[:, 2], header)
        else
            y_values = cpu_txt_data[:, 2] .* internalUnits("CLKT", snapshot_path)
        end

    elseif y_quantity == :clock_time_s

        y_values = cpu_txt_data[:, 3] .* u"s"

    elseif y_quantity == :clock_time_percent

        y_values = cpu_txt_data[:, 4]

    elseif y_quantity == :tot_clock_time_s

        y_values = cpu_txt_data[:, 5] .* u"s"

    elseif y_quantity == :tot_clock_time_percent

        y_values = cpu_txt_data[:, 6]

    else

        throw(ArgumentError("daCPUtxt: I don't recognize the y_quantity :$(y_quantity)"))

    end

    delete_idxs = isnothing(y_log) ? Int64[] : map(iszero, y_values)

    deleteat!(x_axis, delete_idxs)
    deleteat!(y_values, delete_idxs)

    if any(isempty, [x_axis, y_values])

        logging[] && @warn("daCPUtxt: The results of `cpu.txt` are empty")

        return Float64[], Float64[]

    end

    y_axis = isnothing(y_log) ? y_values : log10.(ustrip.(y_log, y_values))

    # Apply smoothing if required
    if !iszero(smooth)
        x_axis, y_axis = smoothWindow(x_axis, y_axis, smooth)
    end

    return x_axis, y_axis

end

"""
    daVirialAccretion(
        sim_data::Simulation;
        <keyword arguments>
    )::NTuple{2,Vector{<:Number}}

Compute the evolution of the accreted mass into a sphere with the virial radius.

# Arguments

  - `sim_data::Simulation`: The [`Simulation`](@ref) struct for the target simulation.
  - `component::Symbol`: Component to compute the accreted mass for. The options are:

      + `:dark_matter` -> Dark matter.
      + `:black_hole`  -> Black holes.
      + `:gas`         -> Gas.
      + `:stellar`     -> Stars.
      + `:all`         -> All the matter.
  - `flux_direction::Symbol=:net`: What flux direction will be plotted. The options are:

      + `:net_mass`     -> Net accreted mass.
      + `:inflow_mass`  -> Inflow mass only.
      + `:outflow_mass` -> Outflow mass only.
  - `halo_idx::Int=1`: Index of the target halo (FoF group). Starts at 1.
  - `tracers::Bool=false`: If tracers will be use to compute the mass accretion.
  - `y_log::Union{Unitful.Units,Nothing}=nothing`: Target unit for integrated mass flux, if you want to apply ``\\log_{10}`` to it. If set to `nothing`, the data from [`computeVirialAccretion`](@ref) is left as is.
  - `smooth::Int=0`: The time series will be smoothed out using `smooth` bins. Set it to 0 if you want no smoothing.
  - `show_progress::Bool=true`: If a progress bar will be shown.

# Returns

  - A Tuple with two elements:

      + A Vector with the physical times.
      + A Vector with the accreted mass at each time.
"""
function daVirialAccretion(
    sim_data::Simulation,
    component::Symbol;
    flux_direction::Symbol=:net,
    halo_idx::Int=1,
    tracers::Bool=false,
    y_log::Union{Unitful.Units,Nothing}=nothing,
    smooth::Int=0,
    show_progress::Bool=true,
)::NTuple{2,Vector{<:Number}}

    request = plotParams(:mass_accretion).request

    # Read the metadata table for the simulation
    simulation_dataframe = DataFrame(sim_data.snapshot_table[sim_data.slice, :])

    # Delete missing snapshots
    filter!(row -> !ismissing(row[:snapshot_paths]), simulation_dataframe)

    # Iterate over each snapshot in the slice
    iterator = eachrow(simulation_dataframe)

    # Check that there are at least 2 snapshots left
    (
        length(iterator) >= 2 ||
        throw(ArgumentError("daVirialAccretion: The given slice, $(sim_data.slice), selected for \
        less than two snapshots. I need at least two snapshots to compute a time series of mass \
        accretion"))
    )

    ################################################################################################
    # First snapshot
    ################################################################################################

    first_snapshot = first(iterator)

    past_dd = makeDataDict(sim_data.path, first_snapshot[:row_id], request, sim_data.snapshot_table)

    ################################################################################################
    # Iterate over the snapshots
    ################################################################################################

    n_frames = length(iterator) - 1

    Δm = Vector{typeof(1.0u"Msun")}(undef, n_frames)

    # Initialize the progress bar
    prog_bar = Progress(
        n_frames,
        dt=0.5,
        desc="Analyzing and plotting the data... ",
        color=:blue,
        barglyphs=BarGlyphs("|#  |"),
        enabled=show_progress,
    )

    for (slice_index, snapshot_data) in pairs(iterator[2:end])

        present_dd = makeDataDict(
            sim_data.path,
            snapshot_data[:row_id],
            request,
            sim_data.snapshot_table,
        )

        if component == :all

            _, m_in_gas, m_out_gas = computeVirialAccretion(
                present_dd,
                past_dd,
                :gas;
                halo_idx,
                tracers=false,
            )

            _, m_in_stars, m_out_stars = computeVirialAccretion(
                present_dd,
                past_dd,
                :stellar;
                halo_idx,
                tracers=false,
            )

            _, m_in_dm, m_out_dm = computeVirialAccretion(
                present_dd,
                past_dd,
                :dark_matter;
                halo_idx,
                tracers=false,
            )

            _, m_in_bh, m_out_bh = computeVirialAccretion(
                present_dd,
                past_dd,
                :black_hole;
                halo_idx,
                tracers=false,
            )

            m_in  = m_in_gas + m_in_stars + m_in_dm + m_in_bh
            m_out = m_out_gas + m_out_stars + m_out_dm + m_out_bh
            δm    = m_in - m_out

        elseif component ∈ [:gas, :stellar, :dark_matter, :black_hole]

            δm, m_in, m_out = computeVirialAccretion(
                present_dd,
                past_dd,
                component;
                halo_idx,
                tracers,
            )

        else

            throw(ArgumentError("daVirialAccretion: `component` can only be :gas, :stellar, \
            :dark_matter, :black_hole or :all, but I got :$(component)"))

        end

        if flux_direction == :net_mass

            Δm[slice_index] = δm

        elseif flux_direction == :inflow_mass

            Δm[slice_index] = m_in

        elseif flux_direction == :outflow_mass

            Δm[slice_index] = m_out

        else

            throw(ArgumentError("daVirialAccretion: `flux_direction` can only be :net_mass, \
            :inflow_mass or :outflow_mass, but I got :$(flux_direction)"))

        end

        past_dd = present_dd

        next!(prog_bar)

    end

    # Compure the time ticks
    t  = sim_data.snapshot_table[sim_data.slice, :physical_times]

    # Compure the time axis
    Δt = diff(t)

    x_axis = t[2:end]
    y_values = Δm ./ Δt

    delete_idxs = isnothing(y_log) ? Int64[] : map(iszero, y_values)

    deleteat!(x_axis, delete_idxs)
    deleteat!(y_values, delete_idxs)

    if any(isempty, [x_axis, y_values])

        logging[] && @warn("daVirialAccretion: The results of `computeVirialAccretion` are empty")

        return Float64[], Float64[]

    end

    y_axis = isnothing(y_log) ? y_values : log10.(ustrip.(y_log, y_values))

    if iszero(smooth)
        return x_axis, y_axis
    else
        return smoothWindow(x_axis, y_axis, smooth)
    end

end

"""
    daDiskAccretion(
        sim_data::Simulation;
        <keyword arguments>
    )::NTuple{2,Vector{<:Number}}

Compute the evolution of the accreted mass into a given galactic disc.

# Arguments

  - `sim_data::Simulation`: The [`Simulation`](@ref) struct for the target simulation.
  - `component::Symbol`: Component to compute the accreted mass for. The options are:

      + `:dark_matter` -> Dark matter.
      + `:black_hole`  -> Black holes.
      + `:gas`         -> Gas.
      + `:stellar`     -> Stars.
      + `:all`         -> All the matter.
  - `flux_direction::Symbol=:net`: What flux direction will be plotted. The options are:

      + `:net_mass`     -> Net accreted mass.
      + `:inflow_mass`  -> Inflow mass only.
      + `:outflow_mass` -> Outflow mass only.
  - `max_r::Unitful.Length=DISK_R`: Radius of the disk.
  - `max_z::Unitful.Length=5.0u"kpc"`: Half height of the disk.
  - `trans_mode::Union{Symbol,Tuple{TranslationType,RotationType,Dict{Symbol,Vector{String}}}}=:all_box`: How to translate and rotate the cells/particles, before filtering with `filter_mode`. For options see [`selectTransformation`](@ref).
  - `tracers::Bool=false`: If tracers will be use to compute the mass accretion.
  - `y_log::Union{Unitful.Units,Nothing}=nothing`: Target unit for integrated mass flux, if you want to apply ``\\log_{10}`` to it. If set to `nothing`, the data from [`computeDiskAccretion`](@ref) is left as is.
  - `smooth::Int=0`: The time series will be smoothed out using `smooth` bins. Set it to 0 if you want no smoothing.
  - `show_progress::Bool=true`: If a progress bar will be shown.

# Returns

  - A Tuple with two elements:

      + A Vector with the physical times.
      + A Vector with the accreted mass at each time.
"""
function daDiskAccretion(
    sim_data::Simulation,
    component::Symbol;
    flux_direction::Symbol=:net,
    max_r::Unitful.Length=DISK_R,
    max_z::Unitful.Length=5.0u"kpc",
    trans_mode::Union{Symbol,Tuple{TranslationType,RotationType,Dict{Symbol,Vector{String}}}}=:all_box,
    tracers::Bool=false,
    y_log::Union{Unitful.Units,Nothing}=nothing,
    smooth::Int=0,
    show_progress::Bool=true,
)::NTuple{2,Vector{<:Number}}

    base_request = plotParams(:mass_accretion).request

    translation, rotation, request = selectTransformation(trans_mode, base_request)

    # Read the metadata table for the simulation
    simulation_dataframe = DataFrame(sim_data.snapshot_table[sim_data.slice, :])

    # Delete missing snapshots
    filter!(row -> !ismissing(row[:snapshot_paths]), simulation_dataframe)

    # Iterate over each snapshot in the slice
    iterator = eachrow(simulation_dataframe)

    # Check that there are at least 2 snapshots left
    (
        length(iterator) >= 2 ||
        throw(ArgumentError("daDiskAccretion: The given slice, $(sim_data.slice), selected for \
        less than two snapshots. I need at least two snapshots to compute a time series of mass \
        accretion"))
    )

    ################################################################################################
    # First snapshot
    ################################################################################################

    first_snapshot = first(iterator)

    past_dd = makeDataDict(sim_data.path, first_snapshot[:row_id], request, sim_data.snapshot_table)

    # Translate the data
    translateData!(past_dd, translation...)

    # Rotate the data
    rotateData!(past_dd, rotation...)

    ################################################################################################
    # Iterate over the snapshots
    ################################################################################################

    n_frames = length(iterator) - 1

    Δm = Vector{typeof(1.0u"Msun")}(undef, n_frames)

    # Initialize the progress bar
    prog_bar = Progress(
        n_frames,
        dt=0.5,
        desc="Analyzing and plotting the data... ",
        color=:blue,
        barglyphs=BarGlyphs("|#  |"),
        enabled=show_progress,
    )

    for (slice_index, snapshot_data) in pairs(iterator[2:end])

        present_dd = makeDataDict(
            sim_data.path,
            snapshot_data[:row_id],
            request,
            sim_data.snapshot_table,
        )

        # Translate the data
        translateData!(present_dd, translation...)

        # Rotate the data
        rotateData!(present_dd, rotation...)

        if component == :all

            _, m_in_gas, m_out_gas = computeDiskAccretion(
                present_dd,
                past_dd,
                :gas;
                max_r,
                max_z,
                tracers=true,
            )

            _, m_in_stars, m_out_stars = computeDiskAccretion(
                present_dd,
                past_dd,
                :stellar;
                max_r,
                max_z,
                tracers=false,
            )

            _, m_in_dm, m_out_dm = computeDiskAccretion(
                present_dd,
                past_dd,
                :dark_matter;
                max_r,
                max_z,
                tracers=false,
            )

            _, m_in_bh, m_out_bh = computeDiskAccretion(
                present_dd,
                past_dd,
                :black_hole;
                max_r,
                max_z,
                tracers=false,
            )

            m_in  = m_in_gas + m_in_stars + m_in_dm + m_in_bh
            m_out = m_out_gas + m_out_stars + m_out_dm + m_out_bh
            δm    = m_in - m_out

        elseif component ∈ [:gas, :stellar, :dark_matter, :black_hole]

            δm, m_in, m_out = computeDiskAccretion(
                present_dd,
                past_dd,
                component;
                max_r,
                max_z,
                tracers,
            )

        else

            throw(ArgumentError("daDiskAccretion: `component` can only be :gas, :stellar, \
            :dark_matter, :black_hole or :all, but I got :$(component)"))

        end

        if flux_direction == :net_mass

            Δm[slice_index] = δm

        elseif flux_direction == :inflow_mass

            Δm[slice_index] = m_in

        elseif flux_direction == :outflow_mass

            Δm[slice_index] = m_out

        else

            throw(ArgumentError("daDiskAccretion: `flux_direction` can only be :net_mass, \
            :inflow_mass or :outflow_mass, but I got :$(flux_direction)"))

        end

        past_dd = present_dd

        next!(prog_bar)

    end

    # Compure the time ticks
    t  = sim_data.snapshot_table[sim_data.slice, :physical_times]

    # Compure the time axis
    Δt = diff(t)

    x_axis = t[2:end]
    y_values = Δm ./ Δt

    delete_idxs = isnothing(y_log) ? Int64[] : map(iszero, y_values)

    deleteat!(x_axis, delete_idxs)
    deleteat!(y_values, delete_idxs)

    if any(isempty, [x_axis, y_values])

        logging[] && @warn("computeDiskAccretion: The results of `computeVirialAccretion` are empty")

        return Float64[], Float64[]

    end

    y_axis = isnothing(y_log) ? y_values : log10.(ustrip.(y_log, y_values))

    if iszero(smooth)
        return x_axis, y_axis
    else
        return smoothWindow(x_axis, y_axis, smooth)
    end

end

"""
    daVSFLaw(
        data_dict::Dict,
        grid::CubicGrid,
        component::Symbol;
        <keyword arguments>
    )::Union{NTuple{2,Vector{<:Float64}},Nothing}

Compute the gas density and the SFR density, used in the volumetric star formation (VSF) law.

# Arguments

  - `data_dict::Dict`: Data dictionary (see [`makeDataDict`](@ref) for the canonical description).
  - `grid::CubicGrid`: Cubic grid.
  - `component::Symbol`: Target component. It can only be one of the elements of [`COMPONENTS`](@ref).
  - `field_type::Symbol=:cells`: If the gas surface density will be calculated assuming the gas is in `:particles` or in Voronoi `:cells`.
  - `age_limit::Unitful.Time=AGE_RESOLUTION`: Age limit for the SFR.
  - `stellar_ff::Function=filterNothing`: Filter function for the stars. See the required signature and examples in `./src/analysis/filters.jl`.
  - `gas_ff::Function=filterNothing`: Filter function for the gas. See the required signature and examples in `./src/analysis/filters.jl`.
  - `mmap_path::String="./"`: Path to store the memory-mapped file if needed (for matrices larger than [`MMAP_THRESHOLD`](@ref)).
  - `m_unit::Unitful.Units=u"Msun"`: Target mass unit.
  - `t_unit::Unitful.Units=u"yr"`: Target time unit.
  - `l_gas_unit::Unitful.Units=u"pc"`: Target length unit for the gas density.
  - `l_stellar_unit::Unitful.Units=u"kpc"`: Target length unit for the SFR density.

# Returns

  - A tuple with two elements:

      + A vector with ``\\log_{10}(\\rho_\\mathrm{H} \\, / \\, \\mathrm{M_\\odot \\, pc^{-3}})``
      + A vector with ``\\log_{10}(\\rho_\\mathrm{SFR} \\, / \\, \\mathrm{M_\\odot \\, yr^{-1} \\, kpc^{-3}})``

    It returns `nothing` if any of the necessary quantities are missing.
"""
function daVSFLaw(
    data_dict::Dict,
    grid::CubicGrid,
    component::Symbol;
    field_type::Symbol=:cells,
    age_limit::Unitful.Time=AGE_RESOLUTION,
    stellar_ff::Function=filterNothing,
    gas_ff::Function=filterNothing,
    mmap_path::String="./",
    m_unit::Unitful.Units=u"Msun",
    t_unit::Unitful.Units=u"yr",
    l_gas_unit::Unitful.Units=u"pc",
    l_stellar_unit::Unitful.Units=u"kpc"
)::Union{NTuple{2,Vector{<:Float64}},Nothing}

    if component ∉ COMPONENTS
        throw(ArgumentError("daVSFLaw: `component` can only be one of the elements \
        of `COMPONENTS` (see `./src/constants/globals.jl`), but I got :$(component)"))
    end

    # Factor to go from stellar density to SFR density
    # log10(ρsfr) = log10(ρ*) - log10Δt
    log10Δt = log10(ustrip(t_unit, age_limit))

    gas_density = quantity3DProjection(
        data_dict,
        grid,
        Symbol(component, :_mass),
        field_type;
        scale_by_volume=true,
        density=l_gas_unit,
        log=false,
        mmap_path,
        filter_function=gas_ff,
    )

    m_gas_factor = ustrip(m_unit, 1.0 * plotParams(Symbol(component, :_mass)).unit)

    stellar_density = quantity3DProjection(
        data_dict,
        grid,
        :stellar_mass,
        :particles;
        scale_by_volume=true,
        density=l_stellar_unit,
        log=false,
        mmap_path,
        filter_function=stellar_ff,
    )

    m_stellar_factor = ustrip(m_unit, 1.0 * plotParams(:stellar_mass).unit)

    x_axis = vec(gas_density) .* m_gas_factor
    y_axis = vec(stellar_density) .* m_stellar_factor

    # Delete NaNs in the data vectors
    x_idxs = map(isnan, x_axis)
    y_idxs = map(isnan, y_axis)

    delete_idxs = x_idxs ∪ y_idxs

    deleteat!(x_axis, delete_idxs)
    deleteat!(y_axis, delete_idxs)

    if any(isempty, [x_axis, y_axis])

        logging[] && @warn("daVSFLaw: The results of `quantity3DProjection` are empty")

        return nothing

    end

    return log10.(x_axis), log10.(y_axis) .- log10Δt

end

@doc raw"""
    daClumpingFactor(
        data_dict::Dict,
        component::Symbol;
        <keyword arguments>
    )::Tuple{Vector{<:Number},Vector{Float64}}

Compute the clumping factor (``C_\rho``), for the number density of `component`, at different volume scales.

```math
C_\rho = \frac{\langle n^2 \rangle}{\langle n \rangle^2} \, ,
```

# Arguments

  - `data_dict::Dict`: Data dictionary (see [`makeDataDict`](@ref) for the canonical description).
  - `component::Symbol`: Target component. It can only be one of the elements of [`COMPONENTS`](@ref) with cell/partcile type :gas.
  - `n_neighbors::Int=32`: Number of neighbors.
  - `filter_function::Function=filterNothing`: Filter function to be applied to `data_dict` before any other computation. See the required signature and examples in `./src/analysis/filters.jl`.
  - `x_log::Union{Unitful.Units,Nothing}=nothing`: Target unit for the volume, if you want to apply ``\log_{10}`` to the volume. If set to `nothing`, the volume is left as is.
  - `y_log::Union{Unitful.Units,Nothing}=nothing`: Target unit for the clumping factor, if you want to apply ``\log_{10}`` to the clumping factor`. If set to `nothing`, the clumping factor is left as is.

# Returns

  - A tuple with two elements:

      + A vector with the the volumes.
      + A vector with the clumping factors.
"""
function daClumpingFactor(
    data_dict::Dict,
    component::Symbol;
    n_neighbors::Int=32,
    filter_function::Function=filterNothing,
    x_log::Union{Unitful.Units,Nothing}=nothing,
    y_log::Union{Unitful.Units,Nothing}=nothing,
)::Tuple{Vector{<:Number},Vector{Float64}}

    if component ∉ COMPONENTS
        throw(ArgumentError("daClumpingFactor: `component` can only be one of the elements \
        of `COMPONENTS` (see `./src/constants/globals.jl`), but I got :$(component)"))
    end

    filtered_dd = filterData(data_dict; filter_function)

    # Compute the number density of the target component
    number_densities = scatterQty(filtered_dd, Symbol(component, :_number_density))

    # Load the position of each cell/particle
    positions = ustrip.(u"kpc", filtered_dd[:gas]["POS "])

    # Compute the volume of each cell/particle
    cell_volumes = filtered_dd[:gas]["MASS"] ./ filtered_dd[:gas]["RHO "]

    if any(isempty, [positions, cell_volumes, number_densities])

        logging[] && @warn("daClumpingFactor: The positions, volumes or number densities are empty")

        return Float64[], Unitful.Volume[]

    end

    # Compute the tree for a nearest neighbor search
    kdtree = KDTree(positions)

    # Find the `n_neighbors` nearest cells/particles to each cell/particle
    idxs, _ = knn(kdtree, positions, n_neighbors, true)

    Cρ = similar(number_densities, Float64)
    V  = similar(cell_volumes)

    (
        allequal(length, [V, Cρ, idxs]) ||
        throw(DomainError("daClumpingFactor: The number densities, volumes, and nearest \
        neighbor indices don't have the same lengths. Something went wrong!"))
    )

    Threads.@threads for (i, idx) in collect(pairs(idxs))
        V[i]  = sum(cell_volumes[idx]; init=0.0u"kpc^3")
        Cρ[i] = computeClumpingFactor(number_densities[idx])
    end

    x_idxs = isnothing(x_log) ? Int64[] : map(iszero, V)
    y_idxs = isnothing(y_log) ? Int64[] : map(iszero, Cρ)

    delete_idxs = x_idxs ∪ y_idxs

    deleteat!(V, delete_idxs)
    deleteat!(Cρ, delete_idxs)

    if any(isempty, [V, Cρ])

        logging[] && @warn("daClumpingFactor: The volumes and/or clumping factors are empty")

        return Float64[], Float64[]

    end

    x_axis = isnothing(x_log) ? V : log10.(ustrip.(x_log, V))
    y_axis = isnothing(y_log) ? Cρ : log10.(Cρ)

    return x_axis, y_axis

end

@doc raw"""
    daClumpingFactorProfile(
        data_dict::Dict,
        component::Symbol,
        grid::LinearGrid;
        <keyword arguments>
    )::Tuple{Vector{<:Unitful.Length},Vector{Float64}}

Compute a clumping factor (``C_\rho``) profile, for the number density of `component`.

```math
C_\rho = \frac{\langle n^2 \rangle}{\langle n \rangle^2} \, ,
```

# Arguments

  - `data_dict::Dict`: Data dictionary (see [`makeDataDict`](@ref) for the canonical description).
  - `component::Symbol`: Target component. It can only be one of the elements of [`COMPONENTS`](@ref) with cell/partcile type :gas.
  - `grid::LinearGrid`: Linear grid.
  - `filter_function::Function=filterNothing`: Filter function to be applied to `data_dict` before any other computation. See the required signature and examples in `./src/analysis/filters.jl`.

# Returns

  - A tuple with two elements:

      + A vector with the central position of each bin.
      + A vector with the clumping factors.
"""
function daClumpingFactorProfile(
    data_dict::Dict,
    component::Symbol,
    grid::LinearGrid;
    filter_function::Function=filterNothing,
)::Tuple{Vector{<:Unitful.Length},Vector{Float64}}

    if component ∉ COMPONENTS
        throw(ArgumentError("daClumpingFactorProfile: `component` can only be one of the elements \
        of `COMPONENTS` (see `./src/constants/globals.jl`), but I got :$(component)"))
    end

    filtered_dd = filterData(data_dict; filter_function)

    # Compute the number density of the target component
    number_densities = scatterQty(filtered_dd, Symbol(component, :_number_density))

    # Load the position of each cell/particle
    positions = filtered_dd[:gas]["POS "]

    # Compute the radial distance of each cell/particle
    distances = colwise(Euclidean(), positions[1:2, :], grid.origin[1:2])

    # Find which cells/particles fall within each bin of `grid`
    n_profile = listHistogram1D(distances, number_densities, grid)

    Cρ = fill(NaN, grid.n_bins)

    Threads.@threads for (i, n) in collect(pairs(n_profile))

        Cρ[i] = computeClumpingFactor(n)

    end

    return grid.x_axis, Cρ

end

"""
    daTrajectory(
        simulation_path::String,
        target_ids::Vector{<:Unsigned},
        cp_type::Symbol;
        <keyword arguments>
    )::Dict{UInt64,Matrix{Quantity}}

Compute the trajectory of a set of cells/particles, given their IDs.

# Arguments

  - `simulation_path::String`: Path to the simulation directory, set in the code variable `OutputDir`.
  - `target_ids::Vector{<:Unsigned}`: IDs of the cells/particles whose trajectories will be computed.
  - `cp_type::Symbol`: Target type of cell/particle. The possibilities are the keys of [`PARTICLE_INDEX`](@ref).
  - `trans_mode::Union{Symbol,Tuple{TranslationType,RotationType,Dict{Symbol,Vector{String}}}}=:all_box`: How to translate and rotate the cells/particles, before filtering with `filter_mode`. For options see [`selectTransformation`](@ref).
  - `filter_mode::Union{Symbol,Tuple{Function,Dict{Symbol,Vector{String}}}}=:all`: Which cells/particles will be selected. For options see [`selectFilter`](@ref).

# Returns

  - A dictionary with the following shape:

      + `target_id::UInt64` -> A matrix with the trajectory of the cell/particle with ID `target_id`. The matrix has 7 rows, where:

          * Row 1: Physical time.
          * Rows 2-4: Position [x, y, z].
          * Rows 5-7: Velocity [vx, vy, vz].

      If the target ID is not found at a given physical time the position and velocity are NaN.
"""
function daTrajectory(
    simulation_path::String,
    target_ids::Vector{<:Unsigned},
    cp_type::Symbol;
    trans_mode::Union{Symbol,Tuple{TranslationType,RotationType,Dict{Symbol,Vector{String}}}}=:all_box,
    filter_mode::Union{Symbol,Tuple{Function,Dict{Symbol,Vector{String}}}}=:all,
)::Dict{UInt64,Matrix{Quantity}}

    # Make a dataframe for the simulation with the following columns:
    #  - DataFrame index         -> :row_id
    #  - Number in the file name -> :numbers
    #  - Scale factor            -> :scale_factors
    #  - Redshift                -> :redshifts
    #  - Physical time           -> :physical_times
    #  - Lookback time           -> :lookback_times
    #  - Snapshot path           -> :snapshot_paths
    #  - Group catalog path      -> :groupcat_paths
    simulation_table = makeSimulationTable(simulation_path)

    # Delete missing snapshots
    filter!(row -> !ismissing(row[:snapshot_paths]), simulation_table)

    (
        isempty(simulation_table) &&
        throw(ArgumentError("daTrajectory: There are no snapshots in $(simulation_path)"))
    )

    (
        haskey(PARTICLE_INDEX, cp_type) ||
        throw(ArgumentError("daTrajectory: `cp_type` has to be one of the keys of PARTICLE_INDEX, \
        but I got :$(cp_type)"))
    )

    base_request = Dict(cp_type => ["ID  ", "POS ", "VEL "])

    translation, rotation, trans_request = selectTransformation(trans_mode, base_request)
    filter_function, request = selectFilter(filter_mode, trans_request)

    # NaNs for the missing IDs
    nan_pos = NaN * u"kpc"
    nan_vel = NaN * u"km * s^-1"

    trajectories = Dict{UInt64, Matrix{Quantity}}()

    Threads.@threads for target_id in target_ids

        trajectory_matrix = Matrix{Quantity}(undef, 7, nrow(simulation_table))

        trajectory_matrix[1, :]   .= simulation_table[!, :physical_times]
        trajectory_matrix[2:4, :] .= nan_pos
        trajectory_matrix[5:7, :] .= nan_vel

        trajectories[target_id] = trajectory_matrix

    end

    for idx_row in 1:nrow(simulation_table)

        # Read the data in the snapshot
        data_dict = makeDataDict(simulation_path, idx_row, request, simulation_table)

        # Translate the data
        translateData!(data_dict, translation...)

        # Rotate the data
        rotateData!(data_dict, rotation...)

        # Filter the data
        filterData!(data_dict; filter_function)

        id = data_dict[cp_type]["ID  "]
        r  = data_dict[cp_type]["POS "]
        v  = data_dict[cp_type]["VEL "]

        Threads.@threads for target_id in target_ids

            (target_id ∉ id) && continue

            index = findfirst(isequal(target_id), id)

            trajectories[target_id][2:4, idx_row] .= r[:, index]
            trajectories[target_id][5:7, idx_row] .= v[:, index]

        end

    end

    return trajectories

end
