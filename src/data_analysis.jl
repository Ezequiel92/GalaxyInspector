####################################################################################################
# Data analysis functions.
####################################################################################################

####################################################################################################
# For use with the `snapshotPlot` function (see `src/pipelines.jl`).
####################################################################################################
#
# Every function is expected to have the following signature:
#
#   foo(data, args...; kw_args...) 
#   ⟶ (x_axis,) or (x_axis, y_axis) or (x_axis, y_axis, z_axis) or `nothing`
#
# Where:
#
#   - data::Dict{Symbol,<:Any}
#   - x_axis::Vector{<:RealOrQty}
#   - y_axis::Vector{<:RealOrQty}
#   - z_axis::Matrix{<:RealOrQty}
#
# `data` is a dictionary with the following entries:
#
#   - :sim_data => A `SimData` object (see `SimData` in `src/constants.jl`).
#   - :snap_data => A `SnapData` object (see `SnapData` in `src/constants.jl`).
#   - :type => (qty => Vector{<:RealOrQty}, qty => Vector{<:RealOrQty}, ...) 
#   - :type => (qty => Vector{<:RealOrQty}, qty => Vector{<:RealOrQty}, ...) 
#   - :type => (qty => Vector{<:RealOrQty}, qty => Vector{<:RealOrQty}, ...) 
#   - ...
#
# where `:type` is a Symbol representing a type of particle (see the definition of `ParticleType` 
# in `src/constants.jl`) and `qty` is a String representing a physical quantity (see the definition 
# of `QUANTITIES` in `src/constants.jl`).
#
# So, `foo` must take a dictionary with the raw data, and return one, two, or three vectors with
# the final values. One vector is for histograms, two for scatter and line plots, and three for 
# heatmaps. It should return `nothing` if `data` has some problem that prevents computation 
# (e.g. is empty).
#
####################################################################################################

"""
    histogram(
        data::Dict{Symbol,<:Any},
        type::Symbol,
        quantity::String;
        <keyword arguments>
    )::Union{NTuple{1,Vector{<:RealOrQty}},Nothing}

Compute the result of applying `func` to the output of [`getSnapshotData`](@ref).

This method is for analysis functions that take a single quantity as input.

# Arguments
- `data::Dict{Symbol,<:Any}`: The raw data returned by [`getSnapshotData`](@ref).
- `type::Symbol`: Particle type, where the possibilities are given by [`ParticleType`](@ref).
- `quantity::String`: GADGET quantity, where the possibilities are given by [`QUANTITIES`](@ref).
- `func::Function = identity`: Function to be applied to the quantity. It has to have the signature:

  `foo(::typeof(data[type][quantity]))::Vector{<:RealOrQty}`

# Returns
- A Tuple with one element:
  - A Vector with the values for each bin.

  It returns `nothing` in case of empty results.

"""
function histogram(
    data::Dict{Symbol,<:Any},
    type::Symbol,
    quantity::String;
    func::Function = identity,
)::Union{NTuple{1,Vector{<:RealOrQty}},Nothing}

    x_axis = func(copy(data[type][quantity]))

    !isempty(x_axis) || return nothing

    return (x_axis,)

end

"""
    histogram(
        data::Dict{Symbol,<:Any},
        type::Symbol,
        quantities::Vector{String};
        <keyword arguments>
    )::Union{NTuple{1,Vector{<:RealOrQty}},Nothing}

Compute the result of applying `func` to the output of [`getSnapshotData`](@ref).

This method is for analysis functions that take several quantities as input.

# Arguments
- `data::Dict{Symbol,<:Any}`: The raw data returned by [`getSnapshotData`](@ref).
- `type::Symbol`: Particle type. The possibilities are given by [`ParticleType`](@ref).
- `quantities::Vector{String}`: GADGET quantities. The possibilities are given by [`QUANTITIES`](@ref).
- `func::Function = identity`: Function to be applied to the quantities. It has to have the signature:

  ```
  foo(
      ::typeof(data[type][quantities[1]]), 
      ::typeof(data[type][quantities[2]]), 
      ...
  )::Vector{<:RealOrQty}
  ```

# Returns
- A Tuple with one element:
  - A Vector with the values for each bin.

  It returns `nothing` in case of empty results.

"""
function histogram(
    data::Dict{Symbol,<:Any},
    type::Symbol,
    quantities::Vector{String};
    func::Function = identity,
)::Union{NTuple{1,Vector{<:RealOrQty}},Nothing}

    data = [copy(data[type][qty]) for qty in quantities]
    x_axis = func(data...)

    # Return `nothing` if the result is empty
    !isempty(x_axis) || return nothing

    return (x_axis,)

end

"""
    correlation(
        data::Dict{Symbol,<:Any},
        x_quantity::Tuple{Symbol,String},
        y_quantity::Tuple{Symbol,String}; 
        <keyword arguments>
    )::Union{NTuple{2,Vector{<:RealOrQty}},Nothing}

Compute the result of applying a couple of functions to the output of [`getSnapshotData`](@ref).
            
It is intended for plots of Y vs. X, where Y and X are vector quantities from [`QUANTITIES`](@ref) 
or are derived from them.

This method is for analysis functions that take a single quantity as input.

# Arguments
- `data::Dict{Symbol,<:Any}`: The raw data returned by [`getSnapshotData`](@ref).
- `x_quantity::Tuple{Symbol,String}`: The tuple: (particle type, quantity), where the posible 
  particle types are given by [`ParticleType`](@ref) and the possible quantities by [`QUANTITIES`](@ref).
- `y_quantity::Tuple{Symbol,String}`: The tuple: (particle type, quantity), where the posible 
  particle types are given by [`ParticleType`](@ref) and the possible quantities by [`QUANTITIES`](@ref).
- `x_function::Function = identity`: Function to be aplied to the x quantity. It has to have the
  signature:

  `foo(::typeof(data[x_quantity[1]][x_quantity[2]]))::Vector{<:RealOrQty}`
- `y_function::Function = identity`: Function to be aplied to the y quantity. It has to have the
  signature:

  `foo(::typeof(data[y_quantity[1]][y_quantity[2]]))::Vector{<:RealOrQty}`

# Returns
- A Tuple with two elements:
  - A Vector with the final values for the x-axis.
  - A Vector with the final values for the y-axis.
  
  It returns `nothing` in case of empty results or a mismatch in lengths.

"""
function correlation(
    data::Dict{Symbol,<:Any},
    x_quantity::Tuple{Symbol,String},
    y_quantity::Tuple{Symbol,String};
    x_function::Function = identity,
    y_function::Function = identity,
)::Union{NTuple{2,Vector{<:RealOrQty}},Nothing}

    x_axis = x_function(copy(data[x_quantity[1]][x_quantity[2]]))
    y_axis = y_function(copy(data[y_quantity[1]][y_quantity[2]]))

    l_x = length(x_axis)
    l_y = length(y_axis)

    # Return `nothing` if there is a mismatch in lengths, or if the result is empty
    (l_x == l_y && l_x + l_y != 0) || return nothing

    return x_axis, y_axis

end

"""
    correlation(
        data::Dict{Symbol,<:Any},
        x_quantities::Vector{Tuple{Symbol,String}},
        y_quantities::Vector{Tuple{Symbol,String}}; 
        <keyword arguments>
    )::Union{NTuple{2,Vector{<:RealOrQty}},Nothing}

Compute the result of applying a couple of functions to the output of [`getSnapshotData`](@ref).
            
It is intended for plots of Y vs. X, where Y and X are vector quantities from [`QUANTITIES`](@ref) 
or are derived from them.
        
This method is for analysis functions that take several quantities as input.

# Arguments
- `data::Dict{Symbol,<:Any}`: The raw data returned by [`getSnapshotData`](@ref).
- `x_quantities::Vector{Tuple{Symbol,String}}`: Vector of tuples: (particle type, quantity), where 
  the posible particle types are given by [`ParticleType`](@ref) and the possible quantities by
  [`QUANTITIES`](@ref).
- `y_quantities::Vector{Tuple{Symbol,String}}`: Vector of tuples: (particle type, quantity), where 
  the posible particle types are given by [`ParticleType`](@ref) and the possible quantities by
  [`QUANTITIES`](@ref).
- `x_function::Function = identity`: Function to be aplied to the x quantities. It has to have the
  signature:

  ```
  foo(
      ::typeof(data[x_quantities[1][1]][x_quantities[1][2]]), 
      ::typeof(data[x_quantities[2][1]][x_quantities[2][2]]),
      ...
  )::Vector{<:RealOrQty}
  ```
- `y_function::Function = identity`: Function to be applied to the y quantities. It has to have the 
  signature:

  ```
  foo(
      ::typeof(data[y_quantities[1][1]][y_quantities[1][2]]), 
      ::typeof(data[y_quantities[2][1]][y_quantities[2][2]]),
      ...
  )::Vector{<:RealOrQty}
  ```

# Returns
- A Tuple with two elements:
  - A Vector with the final values for the x-axis.
  - A Vector with the final values for the y-axis.
  
  It returns `nothing` in case of empty results or a mismatch in lengths.

"""
function correlation(
    data::Dict{Symbol,<:Any},
    x_quantities::Vector{Tuple{Symbol,String}},
    y_quantities::Vector{Tuple{Symbol,String}};
    x_function::Function = identity,
    y_function::Function = identity,
)::Union{NTuple{2,Vector{<:RealOrQty}},Nothing}

    x_data = [copy(data[x_qty[1]][x_qty[2]]) for x_qty in x_quantities]
    y_data = [copy(data[y_qty[1]][y_qty[2]]) for y_qty in y_quantities]

    x_axis = x_function(x_data...)
    y_axis = y_function(y_data...)

    l_x = length(x_axis)
    l_y = length(y_axis)

    # Return `nothing` if there is a mismatch in lengths, or if the result is empty
    (l_x == l_y && l_x + l_y != 0) || return nothing

    return x_axis, y_axis

end

"""
    profile(
        data::Dict{Symbol,<:Any},
        quantity::Tuple{Symbol,String},
        max_radius::Unitful.Length,
        bins::Int64; 
        <keyword arguments>
    )::Union{NTuple{2,Vector},Nothing}

Compute a 2D or 3D profile of `quantity[2]`, or of something derived from it. 

# Arguments
- `data::Dict{Symbol,<:Any}`: The raw data returned by [`getSnapshotData`](@ref). It has to contain
  the quantity "POS" and `quantity[2]`, for the particles of the type specified in `quantity[1]`.
- `quantity::Tuple{Symbol,String}`: The tuple: (particle type, quantity), where the possible particle
  types are given by [`ParticleType`](@ref) and the possible quantities by [`QUANTITIES`](@ref).
- `max_radius::Unitful.Length`: Maximum distance up to which the profile will be calculated.
- `bins::Int64`: Number of subdivisions of [0, `max_radius`].
- `flat::Bool = true`: If the profile will be 2D (default) or 3D.
- `func::Function = identity`: Function to be applied to the quantity. It has to have the signature:

  `foo(::typeof(data[quantity[1]][quantity[2]]))::Vector{<:Unitful.Quantity}`
- `cumulative::Bool = false`: If the profile will be accumulated or not (default).
- `density::Bool = false`: If the profile is of the area/volume density of `quantity[2]` or just of
  `quantity[2]` itself (default).

# Returns
- A Tuple with two elements:
  - A Vector with the distances to each ring or shell.
  - A Vector with the values of the profile at those distances.  
  
  It returns `nothing` in case of empty results.

"""
function profile(
    data::Dict{Symbol,<:Any},
    quantity::Tuple{Symbol,String},
    max_radius::Unitful.Length,
    bins::Int64;
    flat::Bool = true,
    func::Function = identity,
    cumulative::Bool = false,
    density::Bool = false,
)::Union{NTuple{2,Vector},Nothing}

    qty = func(copy(data[quantity[1]][quantity[2]]))
    pos = copy(data[quantity[1]]["POS"])

    # Return `nothing` if any of the neccesary quantities is missing
    !any(isempty, [qty, pos]) || return nothing

    distances = computeDistance(flat ? pos[1:2, :] : pos)

    return computeProfile(distances, qty, max_radius, bins; flat, cumulative, density)

end

"""
    zProfile(
        data::Dict{Symbol,<:Any},
        type::Symbol,
        max_radius::Unitful.Length,
        bins::Int64;
        <keyword arguments>
    )::Union{NTuple{2,Vector},Nothing}

Compute a 2D or 3D metallicity profile. 

# Arguments
- `data::Dict{Symbol,<:Any}`: The raw data returned by [`getSnapshotData`](@ref). It has to contain
  the quantities "Z", "MASS" and "POS", for the particles of the type specified in `type`.
- `type::Symbol`: Particle type. The possibilities are given by [`ParticleType`](@ref).
- `max_radius::Unitful.Length`: Maximum distance up to which the profile will be calculated.
- `bins::Int64`: Number of subdivisions of [0, `max_radius`].
- `flat::Bool = true`: If the profile will be 2D (default) or 3D.
- `cumulative::Bool = false`: If the profile will be accumulated or not (default).
- `density::Bool = false`: If the profile is of the metallicity area/volume density or just of the
  metallicity itself (default).
- `solar::Bool = false`: If the result will be in units of solar metallicity or not (default).

# Returns
- A Tuple with two elements:
  - A Vector with the distances to each ring or shell.
  - A Vector with the values of the profile at those distances.  
  
  It returns `nothing` in case of empty results.

"""
function zProfile(
    data::Dict{Symbol,<:Any},
    type::Symbol,
    max_radius::Unitful.Length,
    bins::Int64;
    flat::Bool = true,
    cumulative::Bool = false,
    density::Bool = false,
    solar::Bool = false,
)::Union{NTuple{2,Vector},Nothing}

    metals = copy(data[type]["Z"])
    mass = copy(data[type]["MASS"])
    pos = copy(data[type]["POS"])

    # Return `nothing` if any of the neccesary quantities is missing
    !any(isempty, [metals, mass, pos]) || return nothing

    distances = computeDistance(flat ? pos[1:2, :] : pos)

    return computeZProfile(
        distances,
        metals,
        mass,
        max_radius,
        bins;
        flat,
        cumulative,
        density,
        solar,
    )

end

@doc raw"""
    cmdf(
        data::Dict{Symbol,<:Any}; 
        <keyword arguments>
    )::Union{NTuple{2,Vector{Float64}},Nothing}

Compute the cumulative metallicity distribution function (CMDF). 
    
The CMDF is calculated separating the stellar metallicity in `bins` windows. For the 
``n``-th window the CMDF is
    
```math
\sum_{i = 1}^n \dfrac{m_i}{M_T} \quad \mathrm{vs.} \quad \bar{Z}_n \, ,
```
    
or, for `normalize_x` = true, 
    
```math
\sum_{i = 1}^n \dfrac{m_i}{M_T} \quad \mathrm{vs.} \quad \dfrac{\bar{Z}_n}{\mathrm{max}(\bar{Z}_n)} \, ,
```
    
where ``M_T`` is the total stellar mass, ``m_i`` the stellar mass of the ``i``-th window and
``\bar{Z}_n`` the mean stellar metallicity of the ``n``-th window.

# Arguments
- `data::Dict{Symbol,<:Any}`: The raw data returned by [`getSnapshotData`](@ref). It has to contain
  the quantities "Z" and "MASS" for the particles of type `:stars`.
- `bins::Int64 = 20`: Number of metallicity windows to be used for the CMDF.
- `normalize_x::Bool = false`: If the x axis will be normalized to its maximum value or not (default).
- `solar::Bool = false`: If the result will be in units of solar metallicity or not (default). 

# Returns
- A Tuple with two elements:
  - A Vector with the mean metallicity of each window.
  - A Vector with the accumulated stellar mass fraction of each window.

  It returns nothing in case of empty results or a mismatch in lengths.

"""
function cmdf(
    data::Dict{Symbol,<:Any};
    bins::Int64 = 20,
    normalize_x::Bool = false,
    solar::Bool = false,
)::Union{NTuple{2,Vector{Float64}},Nothing}

    # Return `nothing` if there are no stellar particles
    data[:sim_data].header.nall[ParticleType[:stars]+1] != 0 || return nothing

    mass = copy(data[:stars]["MASS"])
    metallicity = copy(data[:stars]["Z"])

    l_ma = length(mass)
    l_me = size(metallicity, 2)

    # Return `nothing` if there are a mismatch in lengths, 
    # or if any of the neccesary quantities is missing
    (l_ma == l_me && l_ma + l_me != 0) || return nothing

    # Dimensionless metallicity
    Z = computeMetallicity(metallicity, mass; solar)

    # Magnitude of the masses
    m = ustrip.(mass)

    # Total star mass
    total_m = sum(m)

    # Width of the metallicity bins
    if normalize_x
        # If required, normalize the x axis
        Z = Z ./ maximum(Z)
        width = 1.0 / bins
    else
        width = maximum(Z) / bins
    end

    # Initialize the output arrays
    x_axis = Vector{Float64}(undef, bins)
    y_axis = Vector{Float64}(undef, bins)

    @inbounds for i in eachindex(x_axis, y_axis)

        idx = findall(x -> width * (i - 1) <= x < width * i, Z)

        if isempty(idx)
            x_axis[i] = width * (i - 0.5)
        else
            # Mean metallicity for window `i`
            x_axis[i] = sum(Z[idx]) / length(idx)
        end

        # Mass fraction for window `i`
        y_axis[i] = sum(m[idx]) / total_m

    end

    return x_axis, cumsum(y_axis)

end

@doc raw"""
    kennicuttSchmidt(
        data::Dict{Symbol,<:Any},
        temp_filter::Unitful.Temperature,
        age_filter::Unitful.Time,
        max_r::Unitful.Length; 
        <keyword arguments>
    )::Union{NTuple{2,Vector{<:RealOrQty}},Nothing}

Compute the mass area density and the SFR area density for the [Kennicutt-Schmidt law](https://en.wikipedia.org/wiki/Kennicutt%E2%80%93Schmidt_law). 

The area densities are calculated by projecting the positions of the stars and the gas particles to
the x-y plane. Then the space is subdivided into `bins` concentric rings from 0 to a radius of `max_r`. 
Each ring is of equal width `max_r` / `bins`. 

The ``n``-th ring has an area of
    
```math
A_n = π \, \mathrm{width}^2 \, (2 \, n - 1) \, .
```
    
So, the assigned SFR for that ring is 

```math
\Sigma_\mathrm{SFR}^n = \frac{M_*^n}{\mathrm{age\_filter}\,A_n} \, ,
```
    
where ``M_*^n`` is the total mass of stars younger than `age_filter` within the ring. Equivalently, 
the mass area density of the gas is given by
    
```math
\Sigma_\rho^n = \frac{M_\rho^n}{A_n} \, ,
```
    
where ``M_\rho^n`` is the total mass of gas colder than `temp_filter` within the ring.

# Arguments
- `data::Dict{Symbol,<:Any}`: The raw data returned by [`getSnapshotData`](@ref). It has to contain
  the quantities "Z", "MASS", "POS", "U" and "NE" for the particles of type `:gas`, and the 
  quantities "MASS", "POS" and "AGE" for the particles of type `:stars`.
- `temp_filter::Unitful.Temperature`: Maximum temperature allowed for the gas particles.
- `age_filter::Unitful.Time`: Maximum age allowed for the stars.
- `max_r::Unitful.Length`: Maximum distance up to which the parameters will be calculated.
- `bins::Int64 = 20`: Number of subdivisions of [0, `max_r`].
- `x_unit::Unitful.Units = UnitfulAstro.Msun * UnitfulAstro.pc^-2`: Unit for the surface density of 
  gas.
- `y_unit::Unitful.Units = UnitfulAstro.Msun * UnitfulAstro.yr^-1 * UnitfulAstro.kpc^-2`: Unit for 
  the surface density of star formation rate. 
- `x_log::Bool = true`: If the x axis will be ``\log(\Sigma_\mathrm{gas})`` (default) or just ``\Sigma_\mathrm{gas}``.
- `y_log::Bool = true`: If the y axis will be ``\log(\Sigma_\mathrm{sfr})`` (default) or just ``\Sigma_\mathrm{sfr}``.

# Returns
- A Tuple with two elements:
  - A Vector with the gas area density of each ring.
  - A Vector with the SFR area density of each ring.
  
  It returns nothing in case the results would have less than five data points.

"""
function kennicuttSchmidt(
    data::Dict{Symbol,<:Any},
    temp_filter::Unitful.Temperature,
    age_filter::Unitful.Time,
    max_r::Unitful.Length;
    bins::Int64 = 20,
    x_unit::Unitful.Units = UnitfulAstro.Msun * UnitfulAstro.pc^-2,
    y_unit::Unitful.Units = UnitfulAstro.Msun * UnitfulAstro.yr^-1 * UnitfulAstro.kpc^-2,
    x_log::Bool = true,
    y_log::Bool = true,
)::Union{NTuple{2,Vector{<:RealOrQty}},Nothing}

    # Return `nothing` if there are no stellar particles
    data[:sim_data].header.nall[ParticleType[:stars]+1] != 0 || return nothing

    gas_mass = copy(data[:gas]["MASS"])
    gas_metallicity = copy(data[:gas]["Z"])
    gas_position = copy(data[:gas]["POS"])
    internal_energy = copy(data[:gas]["U"])
    electron_fraction = copy(data[:gas]["NE"])

    star_mass = copy(data[:stars]["MASS"])
    star_position = copy(data[:stars]["POS"])
    birth_times = copy(data[:stars]["AGE"])

    # Return `nothing` if any of the necessary quantities are missing
    !any(
        isempty,
        [
            gas_mass,
            gas_metallicity,
            gas_position,
            internal_energy,
            electron_fraction,
            star_mass,
            star_position,
            birth_times,
        ],
    ) || return nothing

    # Compute the ages of the stars
    ages = computeStellarAge(birth_times, data[:sim_data], data[:snap_data])

    # Compute the distances to the gas and stellar particles
    gas_distances = computeDistance(gas_position)
    star_distances = computeDistance(star_position)

    # Compute the gas temperature
    temperature = computeTemperature(gas_metallicity, gas_mass, internal_energy, electron_fraction)

    # Filter out hot gas particles
    cold_gas_mass = deleteat!(gas_mass, temperature .> temp_filter)
    cold_gas_distance = deleteat!(gas_distances, temperature .> temp_filter)

    # Filter out old stars
    young_star_mass = deleteat!(star_mass, ages .> age_filter)
    young_star_distance = deleteat!(star_distances, ages .> age_filter)

    # Width of the bins
    r_width = max_r / bins

    # Initialize the output arrays
    x_axis = Vector{Unitful.Quantity}(undef, bins)
    y_axis = Vector{Unitful.Quantity}(undef, bins)

    @inbounds for i in eachindex(x_axis, y_axis)

        # Gas particles
        idx_gas = findall(x -> r_width * (i - 1) <= x < r_width * i, cold_gas_distance)
        gas_mass = sum(cold_gas_mass[idx_gas])
        # Gas area density for the `i`-th ring
        x_axis[i] = gas_mass / (π * r_width * r_width * (2.0 * i - 1.0))

        # Stellar particles
        idx_star = findall(x -> r_width * (i - 1) <= x < r_width * i, young_star_distance)
        sfr = sum(young_star_mass[idx_star]) / age_filter
        # SFR area density for the `i`-th ring
        y_axis[i] = sfr / (π * r_width * r_width * (2.0 * i - 1.0))

    end

    # Filter out negative values
    positiveCut!(x_axis, y_axis, keep_edges = false)
    positiveCut!(y_axis, x_axis, keep_edges = false)

    # Return `nothing` if there are less than 5 data points 
    length(x_axis) >= 5 || return nothing

    # Compute logarithm of the data if required
    x_axis = x_log ? log10.(ustrip.(x_unit, x_axis)) : uconvert.(x_unit, x_axis)
    y_axis = y_log ? log10.(ustrip.(y_unit, y_axis)) : uconvert.(y_unit, y_axis)

    return x_axis, y_axis

end

"""
    gasDensityMap(
        data::Dict{Symbol,<:Any},
        plane::String,
        box_size::Unitful.Length; 
        <keyword arguments>
    )::Tuple{Vector{<:Unitful.Length},Vector{<:Unitful.Length},Matrix{Float64}}

Compute a density map of the gas, in the x-y, x-z, or y-z planes. 

# Arguments
- `data::Dict{Symbol,<:Any}`: The raw data returned by [`getSnapshotData`](@ref). It has to contain
 the quantities "MASS", "RHO", "HSML" and "POS", for the particles of type `:gas`.
- `plane::String = "All"`: String indicating which plane will be considered, 
  - `"XY"` ⟶ x-y plane.
  - `"XZ"` ⟶ x-z plane.
  - `"YZ"` ⟶ y-z plane.
- `box_size::Unitful.Quantity`: Size of the region of interest if vacuum boundary conditions were 
  used.
- `sim_cosmo::Bool = false`: If the simulation is cosmological, 
  - `false` ⟶ Newtonian simulation (`ComovingIntegrationOn` = 0).
  - `true` ⟶ Cosmological simulation (`ComovingIntegrationOn` = 1).

# Returns
- A Tuple with three elements:
  - An order list with the abscissa values of the plane.
  - An order list with the ordinate values of the plane.
  - A Matrix with the gas density at each point in the plane.

"""
function gasDensityMap(
    data::Dict{Symbol,<:Any},
    plane::String,
    box_size::Unitful.Length;
    sim_cosmo::Bool = false,
)::Tuple{Vector{<:Unitful.Length},Vector{<:Unitful.Length},Matrix{Float64}}

    # Get the units
    length_unit = unit(box_size)
    mass_unit = unit(data[:gas]["MASS"][1])
    density_unit = mass_unit / length_unit^3

    # Get the magnitudes of the quantities
    mass = @. ustrip(data[:gas]["MASS"])
    position = @. ustrip(length_unit, data[:gas]["POS"])
    hsml = @. ustrip(length_unit, data[:gas]["HSML"])
    ρ = @. ustrip(density_unit, data[:gas]["RHO"])

    # Internal box size (only relevant for periodic boundary conditions)
    boxsize = data[:sim_data].header.boxsize
    conv_factor = internalUnits("POS", data[:sim_data].header; sim_cosmo)

    if boxsize != 0
        # Periodic boundary conditions
        box_limits = round(ustrip(length_unit, (boxsize / 2.0) * conv_factor), sigdigits = 2)
    else
        # Vacuum boundary conditions
        box_limits = round(ustrip(box_size * 1.05), sigdigits = 2)
    end

    # Resolution for binning the grid (number of bins)
    resolution = 1000

    param = mappingParameters(
        x_lim = [-box_limits, box_limits],
        y_lim = [-box_limits, box_limits],
        z_lim = [-box_limits, box_limits],
        Npixels = resolution,
    )

    binning = range(-box_limits, box_limits, length = resolution)

    # Spline kernel used inside GADGET2/3
    # 
    #   Monaghan, J. J., & Lattanzio, J. C. (1985). A refined particle method for astrophysical 
    #   problems. Astronomy and Astrophysics, 149(1), 135–143. 
    #   https://ui.adsabs.harvard.edu/abs/1985A&A...149..135M
    # 
    #   Springel, V. (2005). The cosmological simulation code gadget-2. Monthly Notices of the Royal
    #   Astronomical Society, 364(4), 1105–1134. https://doi.org/10.1111/j.1365-2966.2005.09655.x
    kernel = Cubic()

    if !any(isempty, [mass, position, hsml, ρ])
        if plane == "XY"

            z_axis = log10.(
                sphMapping(position, hsml, mass, ρ, ρ; param, kernel, show_progress = false)
            )

        elseif plane == "XZ"

            # Active rotation (alibi) using Euler angles: Rx(-90°) Ry(0°) Rz(0°)
            position_xz = rotate_3D(position, -90.0, 0.0, 0.0)
            z_axis = log10.(
                sphMapping(position_xz, hsml, mass, ρ, ρ; param, kernel, show_progress = false),
            )

        elseif plane == "YZ"

            # Active rotation (alibi) using Euler angles: Rx(-90°) Ry(0°) Rz(-90°)
            position_yz = rotate_3D(position, -90.0, 0.0, -90.0)
            z_axis = log10.(
                sphMapping(position_yz, hsml, mass, ρ, ρ; param, kernel, show_progress = false),
            )

        else
            throw(ArgumentError("The argument `plane` has to be 'XY', 'XZ' or 'YZ'."))
        end
    else
        # When there are no gas particles return a density of 0
        z_axis = zeros(resolution, resolution)
    end

    x_axis = y_axis = binning .* length_unit

    return x_axis, y_axis, z_axis

end

"""
    particleMap(
        data::Dict{Symbol,<:Any},
        plane::String,
        type::Symbol,
        box_size::Unitful.Length;
        <keyword arguments>
    )::Tuple{Vector{<:Unitful.Length},Vector{<:Unitful.Length},Matrix{Float64}}

Compute a number density map, for particles of type `type` in the x-y, x-z, or y-z planes. 

# Arguments
- `data::Dict{Symbol,<:Any}`: The raw data returned by [`getSnapshotData`](@ref). It has to contain 
  the quantity "POS" for the particles of type `type`.
- `plane::String = "All"`: String indicating which plane will be considered, 
  - `"XY"` ⟶ x-y plane.
  - `"XZ"` ⟶ x-z plane.
  - `"YZ"` ⟶ y-z plane.
- `box_size::Unitful.Quantity`: Size of the region of interest if vacuum boundary conditions were 
  used.
- `sim_cosmo::Bool = false`: If the simulation is cosmological, 
  - `false` ⟶ Newtonian simulation (`ComovingIntegrationOn` = 0).
  - `true` ⟶ Cosmological simulation (`ComovingIntegrationOn` = 1).

# Returns
- A Tuple with three elements:
  - An order list with the abscissa values of the plane.
  - An order list with the ordinate values of the plane.
  - A Matrix with the number density at each point in the plane.

"""
function particleMap(
    data::Dict{Symbol,<:Any},
    plane::String,
    type::Symbol,
    box_size::Unitful.Length;
    sim_cosmo::Bool = false,
)::Tuple{Vector{<:Unitful.Length},Vector{<:Unitful.Length},Matrix{Float64}}

    # Get the data and the units
    length_unit = unit(box_size)
    position = @. ustrip(length_unit, data[type]["POS"])

    # Internal box size, only relevant for periodic boundary conditions
    boxsize = data[:sim_data].header.boxsize
    conv_factor = internalUnits("POS", data[:sim_data].header; sim_cosmo)

    if boxsize != 0
        # Periodic boundary conditions
        box_limits = round(ustrip(length_unit, (boxsize / 2.0) * conv_factor), sigdigits = 2)
    else
        # Vacuum boundary conditions
        box_limits = round(ustrip(box_size * 1.05), sigdigits = 2)
    end

    # Resolution for binning the grid (number of bins)
    resolution = 1000

    binning = range(-box_limits, box_limits, length = resolution)

    particles_present = !isempty(position)
    if particles_present
        x, y, z = position[1, :], position[2, :], position[3, :]
    else
        # When there are no gas particles return a density of 0
        x, y, z = binning, binning, zeros(Float64, resolution, resolution)
    end

    if plane == "XY"
        density = xyz(ash(x, y, rngx = binning, rngy = binning))
    elseif plane == "XZ"
        density = xyz(ash(x, z, rngx = binning, rngy = binning))
    elseif plane == "YZ"
        density = xyz(ash(y, z, rngx = binning, rngy = binning))
    else
        throw(ArgumentError("The argument `plane` has to be 'XY', 'XZ' or 'YZ'."))
    end

    if particles_present
        x_axis = collect(density[1]) .* length_unit
        y_axis = collect(density[2]) .* length_unit
        z_axis = log10.(density[3])
    else
        x_axis, y_axis, z_axis = x .* length_unit, y .* length_unit, z
    end

    return x_axis, y_axis, z_axis

end

####################################################################################################
# For use with the `timeSeriesPlot` function (see `src/pipelines.jl`).
####################################################################################################
#
# Every function is expected to have the following signature:
#
#   foo(sim_data, args...; kw_args...) -> (x_axis, y_axis)
#
# Where:
#
#   * sim_data::SimData
#   * x_axis::Vector{<:RealOrQty}
#   * y_axis::Vector{<:RealOrQty}
#
# So, `foo` must take a dictionary with the raw data, and return two vectors with the final x-axis 
# and y-axis values.
#
####################################################################################################

"""
    qtyEvolution(
        sim_data::SimData,
        x_qty::String,
        y_qty::String,
        type::Union{Symbol,Nothing};
        <keyword arguments>
    )::NTuple{2,Vector}

Compute a time series of `y_qty` and `x_qty`.

# Arguments
- `sim_data::SimData`: Information about the simulation in a [`SimData`](@ref) object.
- `x_qty::String`: Quantity for the x axis. The possibilities are given by the 
  [`computeTimeSeries`](@ref) function:
  - `"clock_time"`   ⟶ Phisical time associated with each snapshot (dimensions: time).
  - `"scale_factor"` ⟶ Scale factor associated with each snapshot (dimensionless).
  - `"redshift"`     ⟶ Redshift associated with each snapshot (dimensionless).
  - `"number"`       ⟶ Number of particles of a given type (dimensionless).
  - `"mass"`         ⟶ Total mass of a given type of particle (dimensions: mass).
  - `"sfr"`          ⟶ Star formation rate in each snapshot (dimensions: mass over time).
- `y_qty::String`: Quantity for the y axis. The possibilities are given by the 
  [`computeTimeSeries`](@ref) function:
  - `"clock_time"`   ⟶ Phisical time associated with each snapshot (dimensions: time).
  - `"scale_factor"` ⟶ Scale factor associated with each snapshot (dimensionless).
  - `"redshift"`     ⟶ Redshift associated with each snapshot (dimensionless).
  - `"number"`       ⟶ Number of particles of a given type (dimensionless).
  - `"mass"`         ⟶ Total mass of a given type of particle (dimensions: mass).
  - `"sfr"`          ⟶ Star formation rate in each snapshot (dimensions: mass over time).
- `x_type::Union{Symbol,Nothing}`: Particle type for the x axis. The posibilities are given by
  [`ParticleType`](@ref) in `src/constants.jl`. It is only relevant for the `"number"` and `"mass"` 
  quantities.
- `y_type::Union{Symbol,Nothing}`: Particle type for the y axis. The posibilities are given by
  [`ParticleType`](@ref) in `src/constants.jl`. It is only relevant for the `"number"` and `"mass"` 
  quantities.
- `sim_cosmo::Bool = false`: If the simulation is cosmological, 
  - `false` ⟶ Newtonian simulation (`ComovingIntegrationOn` = 0).
  - `true` ⟶ Cosmological simulation (`ComovingIntegrationOn` = 1).
- `smooth::Union{Int64,Nothing} = nothing`: The amount of bins to smooth out the output, if `nothing` 
  no smoothing will be done, which is the default.
- `warnings::Bool = true`: If a warning will be given when the data is not as expected, but the 
  function can still run using sane defaults.

# Returns
- A Tuple with two elements:
  - A Vector with the time series of `x_qty`.
  - A Vector with the time series of `y_qty`.

"""
function qtyEvolution(
    sim_data::SimData,
    x_qty::String,
    y_qty::String,
    x_type::Union{Symbol,Nothing},
    y_type::Union{Symbol,Nothing};
    sim_cosmo::Bool = false,
    smooth::Union{Int64,Nothing} = nothing,
    warnings::Bool = true,
)::NTuple{2,Vector}

    paths = getSnapshotPaths(sim_data.base_name, sim_data.path)["snap_paths"]

    if x_qty ∈ ["redshift", "scale_factor"] && !sim_cosmo
        warnings || @warn(
            "The quantities `redshift` and `scale_factor` can only be used with cosmological \
            simulations. Defaulting to clock time for the x axis."
        )
        x_qty = "clock_time"
    end

    if y_qty ∈ ["redshift", "scale_factor"] && !sim_cosmo
        warnings || @warn(
            "The quantities `redshift` and `scale_factor` can only be used with cosmological \
            simulations. Defaulting to clock time for the y axis."
        )
        y_qty = "clock_time"
    end

    # Compute the times series for the x axis 
    x_axis = computeTimeSeries(paths, x_qty, x_type, sim_cosmo; sim_data.filter_function, warnings)

    # Compute the times series for the y axis 
    y_axis = computeTimeSeries(paths, y_qty, y_type, sim_cosmo; sim_data.filter_function, warnings)

    if smooth === nothing
        return x_axis[sim_data.idx], y_axis[sim_data.idx]
    else
        # Smooth out the results
        return smoothWindow(x_axis[sim_data.idx], y_axis[sim_data.idx], smooth)
    end

end
