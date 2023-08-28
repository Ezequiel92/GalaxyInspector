####################################################################################################
# Data analysis functions.
####################################################################################################

####################################################################################################
# For the `snapshotPlot` function in `./src/pipelines.jl`.
####################################################################################################
#
# A data analysis functions for `snapshotPlot` must take a dictionary with the following shape:
# 
#   + `:sim_data`          -> ::Simulation (see [`Simulation`](@ref)).
#   + `:snap_data`         -> ::Snapshot (see [`Snapshot`](@ref)).
#   + `:gc_data`           -> ::GroupCatalog (see [`GroupCatalog`](@ref)).
#   + `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
#   + `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
#   + `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
#   + ...
#   + `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
#   + `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
#   + `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
#   + ...
#
# and return one or more vectors or matrices with the processed data. It should return `nothing` 
# if the input data has some problem that prevents computation (e.g. is empty).
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
    daVelocityProfile(
        data_dict::Dict,
        R::Unitful.Length,
    )::Tuple{Vector{<:Unitful.Length},Vector{<:Unitful.Velocity}}

Compute a circular velocity profile.

# Arguments

  - `data_dict::Dict`: A dictionary with the following shape:

      + `:sim_data`          -> ::Simulation (see [`Simulation`](@ref)).
      + `:snap_data`         -> ::Snapshot (see [`Snapshot`](@ref)).
      + `:gc_data`           -> ::GroupCatalog (see [`GroupCatalog`](@ref)).
      + `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
      + `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
      + `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
      + ...
      + `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
      + `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
      + `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
      + ...
  - `R::Unitful.Length`: Radius for the profile.

# Returns

  - A tuple with two elements:

      + A vector with the distances to each star.
      + A vector with the circular velocity of each star.
"""
function daVelocityProfile(
    data_dict::Dict,
    R::Unitful.Length,
)::Tuple{Vector{<:Unitful.Length},Vector{<:Unitful.Velocity}}

    # Compute the circular velocity and the stellar radial distances, in the order of the snapshot 
    r, vcirc = computeStellarVcirc(data_dict)

    # Only leave the data within a sphere of radius `R`
    rangeCut!(r, vcirc, (0.0u"kpc", R))

    # Sort the arrays
    idx = sortperm(r)

    return r[idx], vcirc[idx]

end

"""
    daKennicuttSchmidt(
        data_dict::Dict,
        grid::CircularGrid,
    )::Union{Tuple{Vector{<:SurfaceDensity},Vector{<:MassFlowDensity}},Nothing}

Compute the gas mass and SFR surface densities, used in the Kennicutt-Schmidt law.

# Arguments

  - `data_dict::Dict`: A dictionary with the following shape:

      + `:sim_data`          -> ::Simulation (see [`Simulation`](@ref)).
      + `:snap_data`         -> ::Snapshot (see [`Snapshot`](@ref)).
      + `:gc_data`           -> ::GroupCatalog (see [`GroupCatalog`](@ref)).
      + `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
      + `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
      + `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
      + ...
      + `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
      + `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
      + `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
      + ...

  - `grid::CircularGrid`: Circular grid.

# Returns

  - A tuple with two elements:

      + A vector with the gas mass surface density of each ring.
      + A vector with the SFR surface density of each ring.

    It returns `nothing` if any of the necessary quantities are missing.

# References

R. C. Kennicutt (1998). *The Global Schmidt Law in Star-forming Galaxies*. The Astrophysical Journal, **498(2)**, 541-552. [doi:10.1086/305588](https://doi.org/10.1086/305588)
"""
function daKennicuttSchmidt(
    data_dict::Dict,
    grid::CircularGrid,
)::Union{Tuple{Vector{<:SurfaceDensity},Vector{<:MassFlowDensity}},Nothing}

    gas_masses = data_dict[:gas]["MASS"]
    gas_positions = data_dict[:gas]["POS "]
    star_positions = data_dict[:stars]["POS "]

    # Return `nothing` if any of the necessary quantities are missing
    !any(isempty, [gas_masses, gas_positions, star_positions]) || return nothing

    # Compute the gas mass surface density
    gas_mass_density = computeProfile(gas_positions, gas_masses, grid; total=true, density=true)

    # Compute the SFR surface density
    sfr_density = computeProfile(
        star_positions,
        computeSFR(data_dict; age_resol=AGE_RESOLUTION_ρ),
        grid;
        total=true,
        density=true,
    )

    return gas_mass_density, sfr_density

end

"""
    function daBigiel2008(
        data_dict::Dict,
        grid::CircularGrid,
        quantity::Symbol,
    )::Union{Tuple{Vector{<:SurfaceDensity},Vector{<:MassFlowDensity}},Nothing}

Compute the (molecular or neutral) gas mass and SFR surface densities, used in the Kennicutt-Schmidt law.

# Arguments

  - `data_dict::Dict`: A dictionary with the following shape:

      + `:sim_data`          -> ::Simulation (see [`Simulation`](@ref)).
      + `:snap_data`         -> ::Snapshot (see [`Snapshot`](@ref)).
      + `:gc_data`           -> ::GroupCatalog (see [`GroupCatalog`](@ref)).
      + `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
      + `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
      + `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
      + ...
      + `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
      + `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
      + `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
      + ...

  - `grid::CircularGrid`: Circular grid.
  - `quantity::Symbol`: Quantity for the x axis. The possibilities are:

      + `:molecular_area_density` -> Molecular hydrogen area mass density.
      + `:neutral_area_density`   -> Neutral hydrogen area mass density.

# Returns

  - A tuple with two elements:

      + A vector with the gas mass surface density of each ring.
      + A vector with the SFR surface density of each ring.

    It returns `nothing` if any of the necessary quantities are missing.

# References

F. Bigiel et al. (2008). *THE STAR FORMATION LAW IN NEARBY GALAXIES ON SUB-KPC SCALES*. The Astrophysical Journal, **136(6)**, 2846. [doi:10.1088/0004-6256/136/6/2846](https://doi.org/10.1088/0004-6256/136/6/2846)
"""
function daBigiel2008(
    data_dict::Dict,
    grid::CircularGrid,
    quantity::Symbol,
)::Union{Tuple{Vector{<:SurfaceDensity},Vector{<:MassFlowDensity}},Nothing}

    gas_positions = data_dict[:gas]["POS "]
    star_positions = data_dict[:stars]["POS "]

    # Return `nothing` if any of the necessary quantities are missing
    !any(isempty, [gas_positions, star_positions]) || return nothing

    if quantity == :molecular_area_density
        gas_masses = computeMolecularMass(data_dict)
    elseif quantity == :neutral_area_density
        gas_masses = computeNeutralMass(data_dict)
    else
        throw(ArgumentError("daBigiel2008: `quantity` can only be :molecular_area_density \
        and :neutral_area_density, but I got :$(quantity)"))
    end

    # Compute the gas mass surface density
    gas_mass_density = computeProfile(gas_positions, gas_masses, grid; total=true, density=true)

    # Compute the SFR surface density
    sfr_density = computeProfile(
        star_positions,
        computeSFR(data_dict; age_resol=AGE_RESOLUTION_ρ),
        grid;
        total=true,
        density=true,
    )

    return gas_mass_density, sfr_density

end

"""
    daMolla2015(
        data_dict::Dict,
        grid::CircularGrid,
        quantity::Symbol,
    )::Union{
        Tuple{
            Vector{<:Unitful.Length},
            <:Union{Vector{<:SurfaceDensity},Vector{<:MassFlowDensity},Vector{Float64}}
        },
        Nothing,
    }

Compute a profile for the Milky Way, compatible with the experimental data in Mollá et al. (2015).

# Arguments

  - `data_dict::Dict`: A dictionary with the following shape:

      + `:sim_data`          -> ::Simulation (see [`Simulation`](@ref)).
      + `:snap_data`         -> ::Snapshot (see [`Snapshot`](@ref)).
      + `:gc_data`           -> ::GroupCatalog (see [`GroupCatalog`](@ref)).
      + `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
      + `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
      + `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
      + ...
      + `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
      + `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
      + `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
      + ...

  - `grid::CircularGrid`: Circular grid.
  - `quantity::Symbol`: Quantity for the y axis. The options are:

      + `:stellar_area_density`   -> Stellar area mass density.
      + `:molecular_area_density` -> Molecular hydrogen area mass density.
      + `:atomic_area_density`    -> Atomic hydrogen area mass density.
      + `:sfr_area_density`       -> Star formation rate area density, for the last `AGE_RESOLUTION_ρ`.
      + `:O_stellar_abundance`    -> Stellar abundance of oxygen, as 12 + log10(O / H).
      + `:N_stellar_abundance`    -> Stellar abundance of nitrogen, as 12 + log10(N / H).
      + `:C_stellar_abundance`    -> Stellar abundance of carbon, as 12 + log10(C / H).

# Returns

  - A tuple with two elements:

      + A vector with the position of each ring.
      + A vector with the `quantity` area density of each ring.

    It returns `nothing` if any of the necessary quantities are missing.

# References

M. Mollá et al. (2015). *Galactic chemical evolution: stellar yields and the initial mass function*. Monthly Notices of the Royal Astronomical Society **451(4)**, 3693–3708. [doi:10.1093/mnras/stv1102](https://doi.org/10.1093/mnras/stv1102)
"""
function daMolla2015(
    data_dict::Dict,
    grid::CircularGrid,
    quantity::Symbol,
)::Union{
    Tuple{
        Vector{<:Unitful.Length},
        <:Union{Vector{<:SurfaceDensity},Vector{<:MassFlowDensity},Vector{Float64}},
    },
    Nothing,
}

    if quantity == :stellar_area_density

        positions = data_dict[:stars]["POS "]
        masses = data_dict[:stars]["MASS"]
        norm_values = Number[]
        f = identity
        density = true

    elseif quantity == :sfr_area_density

        positions = data_dict[:stars]["POS "]
        masses = computeSFR(data_dict; age_resol=AGE_RESOLUTION_ρ)
        norm_values = Number[]
        f = identity
        density = true

    elseif quantity == :molecular_area_density

        positions = data_dict[:gas]["POS "]
        masses = computeMolecularMass(data_dict)
        norm_values = Number[]
        f = identity
        density = true

    elseif quantity == :atomic_area_density

        positions = data_dict[:gas]["POS "]
        masses = computeAtomicMass(data_dict)
        norm_values = Number[]
        f = identity
        density = true

    elseif quantity == :O_stellar_abundance

        positions = data_dict[:stars]["POS "]
        masses = computeElementMass(data_dict, :stars, :O) ./ AtomicWeight[:O]
        norm_values = computeElementMass(data_dict, :stars, :H) ./ AtomicWeight[:H]
        f = x -> 12 .+ log10.(x)
        density = false

    elseif quantity == :N_stellar_abundance

        positions = data_dict[:stars]["POS "]
        masses = computeElementMass(data_dict, :stars, :N) ./ AtomicWeight[:N]
        norm_values = computeElementMass(data_dict, :stars, :H) ./ AtomicWeight[:H]
        f = x -> 12 .+ log10.(x)
        density = false

    elseif quantity == :C_stellar_abundance

        positions = data_dict[:stars]["POS "]
        masses = computeElementMass(data_dict, :stars, :N) ./ AtomicWeight[:N]
        norm_values = computeElementMass(data_dict, :stars, :H) ./ AtomicWeight[:H]
        f = x -> 12 .+ log10.(x)
        density = false

    else

        throw(ArgumentError("daa2015: I don't recognize the quantity :$(quantity)"))

    end

    # Return `nothing` if any of the necessary quantities are missing
    !any(isempty, [positions, masses]) || return nothing

    density_profile = computeProfile(positions, masses, grid; norm_values, f, total=true, density)

    return grid.grid, density_profile

end

"""
    daStellarHistory(
        data_dict::Dict; 
        <keyword arguments>
    )::Union{Tuple{Vector{<:Unitful.Time},Vector{<:Number}},Nothing}

Compute the evolution of a given stellar `quantity` using the stellar ages at a given instant in time.

# Arguments

  - `data_dict::Dict`: A dictionary with the following shape:

      + `:sim_data`          -> ::Simulation (see [`Simulation`](@ref)).
      + `:snap_data`         -> ::Snapshot (see [`Snapshot`](@ref)).
      + `:gc_data`           -> ::GroupCatalog (see [`GroupCatalog`](@ref)).
      + `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
      + `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
      + `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
      + ...
      + `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
      + `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
      + `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
      + ...

  - `quantity::Symbol=:sfr`: Which quantity will be calculated. The options are:

      + `:sfr`          -> The star formation rate.
      + `:ssfr`         -> The specific star formation rate.
      + `:stellar_mass` -> Stellar mass.
  - `n_bins::Int64=50`: Number of bins (time intervals).

# Returns

  - A tuple with two elements:

      + A vector with the physical time of each bin.
      + A vector with the values of `quantity` for each bin.
"""
function daStellarHistory(
    data_dict::Dict;
    quantity::Symbol=:sfr,
    n_bins::Int64=50,
)::Union{Tuple{Vector{<:Unitful.Time},Vector{<:Number}},Nothing}

    birth_ticks = data_dict[:stars]["GAGE"]
    masses = data_dict[:stars]["MASS"]

    # Return `nothing` if any of the necessary quantities are missing
    !any(isempty, [birth_ticks, masses]) || return nothing

    # Compute the stellar birth dates
    if data_dict[:sim_data].cosmological
        # Go from scale factor to physical time
        birth_times = computeTime(birth_ticks, data_dict[:snap_data].header)
    else
        birth_times = birth_ticks
    end

    # Compute the total stellar mass in each time bin
    grid = CircularGrid(maximum(birth_times), n_bins; shift=minimum(birth_times))
    stellar_masses = histogram1D(birth_times, masses, grid; empty_nan=false)

    # Compute the time axis
    min, max = extrema(birth_times)
    bin_width = (max - min) / n_bins
    x_axis = collect(range(min + (bin_width / 2), length=n_bins, step=bin_width))

    # Compute the stellar quantity
    if quantity == :sfr

        y_axis = stellar_masses ./ bin_width

    elseif quantity == :ssfr

        accu_mass = cumsum(stellar_masses)
        y_axis = (stellar_masses ./ bin_width) ./ accu_mass

    elseif quantity == :stellar_mass

        y_axis = cumsum(stellar_masses)

    else

        throw(ArgumentError("daStellarHistory: `quantity` can only be :sfr, :ssfr, \
        or :stellar_mass, but I got :$(quantity)"))

    end

    return x_axis, y_axis

end

"""
    daCircularityHistogram(
        data_dict::Dict,
        grid::LinearGrid,
    )::Union{NTuple{2,Vector{Float64}},Nothing}

Compute a count histogram of the stellar circularity, normalized to the maximum number of counts.

# Arguments

  - `data_dict::Dict`: A dictionary with the following shape:

      + `:sim_data`          -> ::Simulation (see [`Simulation`](@ref)).
      + `:snap_data`         -> ::Snapshot (see [`Snapshot`](@ref)).
      + `:gc_data`           -> ::GroupCatalog (see [`GroupCatalog`](@ref)).
      + `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
      + `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
      + `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
      + ...
      + `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
      + `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
      + `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
      + ...

  - `grid::LinearGrid`: Linear grid.

# Returns

  - A tuple with two elements:

      + A vector with the circularity of each bin.
      + A vector with the counts, normalized to the maximum value.
"""
function daCircularityHistogram(
    data_dict::Dict,
    grid::LinearGrid,
)::Union{NTuple{2,Vector{Float64}},Nothing}

    # Compute the stellar circularity
    circularity = computeStellarCircularity(data_dict)

    !isempty(circularity) || return nothing

    # Compute the circularity histogram
    counts = histogram1D(circularity, grid)

    # To normalize the counts
    max_counts = maximum(counts)

    return grid.grid, counts ./ max_counts

end

"""
    daDensity2DHistogram(
        data_dict::Dict,
        grid::SquareGrid,
        quantity::Symbol;
        <keyword arguments>
    )::Tuple{Vector{<:Unitful.Length},Vector{<:Unitful.Length},Matrix{Float64}}

Compute a 2D density histogram.

!!! note

By default, ``M_\\odot \, kpc^-2`` is used as unit of density, so the output will be ``\\log10(\\rho / M_\\odot \, kpc^-2)``.

# Arguments

  - `data_dict::Dict`: A dictionary with the following shape:

      + `:sim_data`          -> ::Simulation (see [`Simulation`](@ref)).
      + `:snap_data`         -> ::Snapshot (see [`Snapshot`](@ref)).
      + `:gc_data`           -> ::GroupCatalog (see [`GroupCatalog`](@ref)).
      + `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
      + `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
      + `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
      + ...
      + `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
      + `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
      + `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
      + ...

  - `grid::SquareGrid`: Square grid.
  - `quantity::Symbol`: For which quantity the density will be calculated. The possibilities are:

      + `:stellar_mass`   -> Stellar mass.
      + `:gas_mass`       -> Gas mass.
      + `:dm_mass`        -> Dark matter mass.
      + `:bh_mass`        -> Black hole mass.
      + `:molecular_mass` -> Molecular hydrogen (H₂) mass.
      + `:atomic_mass`    -> Atomic hydrogen (HI) mass.
      + `:ionized_mass`   -> Ionized hydrogen (HII) mass.
      + `:neutral_mass`   -> Neutral hydrogen (HI + H₂) mass.
  - `projection_plane::Symbol=:xy`: To which plane the cells/particles will be projected. The options are :xy, :xz, and :yz.
  - `smooth::Bool=false`: If the results will be smooth out using the [`cubicSplineKernel`](@ref) kernel.
  - `neighbors::Int64=32`: Number of neighbors for the 2D smoothing (only relevant if `smooth` = true).

# Returns

  - A tuple with three elements:

      + A vector with the x coordinates of the grid.
      + A vector with the y coordinates of the grid.
      + A matrix with the values of density in each bin.
"""
function daDensity2DHistogram(
    data_dict::Dict,
    grid::SquareGrid,
    quantity::Symbol;
    projection_plane::Symbol=:xy,
    smooth::Bool=false,
    neighbors::Int64=32,
)::Tuple{Vector{<:Unitful.Length},Vector{<:Unitful.Length},Matrix{Float64}}

    # Set the cell/particle type
    if quantity ∈ [:gas_mass, :molecular_mass, :atomic_mass, :ionized_mass, :neutral_mass]
        type_symbol = :gas
    elseif quantity == :stellar_mass
        type_symbol = :stars
    elseif quantity == :dm_mass
        type_symbol = :halo
    elseif quantity == :bh_mass
        type_symbol = :black_hole
    else
        throw(ArgumentError("daDensity2DHistogram: I don't recognize the quantity :$(quantity)"))
    end

    # Compute the masses
    masses = scatterQty(data_dict, quantity)

    # Load the positions
    positions = data_dict[type_symbol]["POS "]

    # If any of the necessary quantities are missing return an empty density field
    if any(isempty, [masses, positions])
        return grid.x_ticks, grid.y_ticks, fill(NaN, size(grid.grid))
    end

    # Set the density unit
    ρ_unit = u"Msun/kpc^2"

    # Project the cell/particles to the chosen plane
    if projection_plane == :xy
        pos_2D = positions[[1, 2], :]
    elseif projection_plane == :xz
        pos_2D = positions[[1, 3], :]
    elseif projection_plane == :yz
        pos_2D = positions[[2, 3], :]
    else
        throw(ArgumentError("daDensity2DHistogram: The argument `projection_plane` must be \
        :xy, :xz or :yz, but I got :$(projection_plane)"))
    end

    # Allocate memory
    density = similar(grid.grid, Number)

    if smooth

        # Spline kernel used in Arepo
        # 
        # Monaghan, J. J., & Lattanzio, J. C. (1985). A refined particle method for astrophysical
        # problems. Astronomy and Astrophysics, 149(1), 135–143. 
        # https://ui.adsabs.harvard.edu/abs/1985A&A...149..135M
        # 
        # Springel, V. (2005). The cosmological simulation code gadget-2. Monthly Notices of the Royal
        # Astronomical Society, 364(4), 1105–1134. https://doi.org/10.1111/j.1365-2966.2005.09655.x
        kernel(q, h) = cubicSplineKernel(q, h)

        # Reshape the grid
        physical_grid = Matrix{Float64}(undef, 2, grid.n_bins * grid.n_bins)
        @inbounds for i in eachindex(grid.grid)
            physical_grid[1, i] = ustrip(u"kpc", grid.grid[i][1])
            physical_grid[2, i] = ustrip(u"kpc", grid.grid[i][2])
        end

        # Compute the tree for a nearest neighbor search
        kdtree = KDTree(ustrip.(u"kpc", pos_2D))

        # Find nearest neighbors
        n_idxs, n_dists = knn(kdtree, physical_grid, neighbors)

        # Compute the smoothing lengths
        if "SOFT" ∈ keys(data_dict[type_symbol])
            smoothing_lengths = data_dict[type_symbol]["SOFT"]
        elseif "RHO " ∈ keys(data_dict[type_symbol])
            densities = data_dict[type_symbol]["RHO "]
            masses = data_dict[type_symbol]["MASS"]
            smoothing_lengths = [cbrt(m / ((4 / 3) * π * ρ)) for (m, ρ) in zip(masses, densities)]
        else
            throw(ArgumentError("daDensity2DHistogram: Neither the \"SOFT\" or \"RHO \" blocks \
            where present for the cell/particle type :$(type_symbol), and I need one of them to \
            smooth out the density histogram"))
        end

        # Compute the density in each bin
        @inbounds for i in eachindex(grid.grid)
            n_idx = n_idxs[i]
            n_dist = n_dists[i]

            hs = ustrip.(u"kpc", smoothing_lengths[n_idx])
            qs = n_dist ./ hs
            ws = kernel.(qs, hs) * u"kpc^-2"

            density[i] = sum(masses[n_idx] .* ws; init=zero(1.0 * ρ_unit))
        end

    else

        # Compute the 2D histogram
        total = histogram2D(pos_2D, masses, grid)
        density = total ./ grid.bin_area

    end

    # Apply log10 to enhance the contrast
    values = log10.(ustrip.(ρ_unit, density))

    # The transpose and reverse operation are to conform to the way heatmap! expect the matrix to be structured
    z_axis = reverse!(transpose(values), dims=2)

    return grid.x_ticks, grid.y_ticks, z_axis

end

"""
    daScatterDensity(
        data_dict::Dict, 
        x_quantity::Symbol, 
        y_quantity::Symbol,
        x_range::NTuple{2,<:Number},
        y_range::NTuple{2,<:Number},
        n_bins::Int64,
    )::Tuple{Vector{<:Number},Vector{<:Number},Matrix{Int64}}

Compute a 2D count histogram.

# Arguments

  - `data_dict::Dict`: A dictionary with the following shape:

      + `:sim_data`          -> ::Simulation (see [`Simulation`](@ref)).
      + `:snap_data`         -> ::Snapshot (see [`Snapshot`](@ref)).
      + `:gc_data`           -> ::GroupCatalog (see [`GroupCatalog`](@ref)).
      + `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
      + `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
      + `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
      + ...
      + `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
      + `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
      + `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
      + ...

  - `x_quantity::Symbol`: Quantity for the x axis. The possibilities are:

      + `:stellar_mass`             -> Stellar mass.
      + `:gas_mass`                 -> Gas mass.
      + `:dm_mass`                  -> Dark matter mass.
      + `:bh_mass`                  -> Black hole mass.
      + `:molecular_mass`           -> Molecular hydrogen (H₂) mass.
      + `:atomic_mass`              -> Atomic hydrogen (HI) mass.
      + `:ionized_mass`             -> Ionized hydrogen (HII) mass.
      + `:neutral_mass`             -> Neutral hydrogen (HI + H₂) mass.
      + `:molecular_fraction`       -> Gas mass fraction of molecular hydrogen.
      + `:atomic_fraction`          -> Gas mass fraction of atomic hydrogen.
      + `:ionized_fraction`         -> Gas mass fraction of ionized hydrogen.
      + `:neutral_fraction`         -> Gas mass fraction of neutral hydrogen.
      + `:gas_mass_density`         -> Gas mass density.
      + `:gas_number_density`       -> Gas number density.
      + `:molecular_number_density` -> Molecular hydrogen number density.
      + `:atomic_number_density`    -> Atomic hydrogen number density.
      + `:ionized_number_density`   -> Ionized hydrogen number density.
      + `:neutral_number_density`   -> Neutral hydrogen number density.
      + `:gas_metallicity`          -> Mass fraction of all elements above He in the gas (solar units).
      + `:stellar_metallicity`      -> Mass fraction of all elements above He in the stars (solar units).
      + `:X_gas_abundance`          -> Gas abundance of element X, as 12 + log10(X / H). The possibilities are the keys of [`ElementIndex`](@ref).
      + `:X_stellar_abundance`      -> Stellar abundance of element X, as 12 + log10(X / H). The possibilities are the keys of [`ElementIndex`](@ref).
      + `:stellar_radial_distance`  -> Distance of every stellar particle to the origin.
      + `:gas_radial_distance`      -> Distance of every gas cell to the origin.
      + `:dm_radial_distance`       -> Distance of every dark matter particle to the origin.
      + `:stellar_xy_distance`      -> Projected distance of every stellar particle to the origin.
      + `:gas_xy_distance`          -> Projected distance of every gas cell to the origin.
      + `:dm_xy_distance`           -> Projected distance of every dark matter particle to the origin.
      + `:stellar_circularity`      -> Stellar circularity.
      + `:stellar_vcirc`            -> Stellar circular velocity.
      + `:stellar_age`              -> Stellar age.
      + `:sfr`                      -> The star formation rate of the last `AGE_RESOLUTION`.
      + `:ssfr`                     -> The specific star formation rate of the last `AGE_RESOLUTION`.
      + `:temperature`              -> Gas temperature, as log10(T / K).
  - `y_quantity::Symbol`: Quantity for the y axis. The possibilities are:

      + `:stellar_mass`             -> Stellar mass.
      + `:gas_mass`                 -> Gas mass.
      + `:dm_mass`                  -> Dark matter mass.
      + `:bh_mass`                  -> Black hole mass.
      + `:molecular_mass`           -> Molecular hydrogen (H₂) mass.
      + `:atomic_mass`              -> Atomic hydrogen (HI) mass.
      + `:ionized_mass`             -> Ionized hydrogen (HII) mass.
      + `:neutral_mass`             -> Neutral hydrogen (HI + H₂) mass.
      + `:molecular_fraction`       -> Gas mass fraction of molecular hydrogen.
      + `:atomic_fraction`          -> Gas mass fraction of atomic hydrogen.
      + `:ionized_fraction`         -> Gas mass fraction of ionized hydrogen.
      + `:neutral_fraction`         -> Gas mass fraction of neutral hydrogen.
      + `:gas_mass_density`         -> Gas mass density.
      + `:gas_number_density`       -> Gas number density.
      + `:molecular_number_density` -> Molecular hydrogen number density.
      + `:atomic_number_density`    -> Atomic hydrogen number density.
      + `:ionized_number_density`   -> Ionized hydrogen number density.
      + `:neutral_number_density`   -> Neutral hydrogen number density.
      + `:gas_metallicity`          -> Mass fraction of all elements above He in the gas (solar units).
      + `:stellar_metallicity`      -> Mass fraction of all elements above He in the stars (solar units).
      + `:X_gas_abundance`          -> Gas abundance of element X, as 12 + log10(X / H). The possibilities are the keys of [`ElementIndex`](@ref).
      + `:X_stellar_abundance`      -> Stellar abundance of element X, as 12 + log10(X / H). The possibilities are the keys of [`ElementIndex`](@ref).
      + `:stellar_radial_distance`  -> Distance of every stellar particle to the origin.
      + `:gas_radial_distance`      -> Distance of every gas cell to the origin.
      + `:dm_radial_distance`       -> Distance of every dark matter particle to the origin.
      + `:stellar_xy_distance`      -> Projected distance of every stellar particle to the origin.
      + `:gas_xy_distance`          -> Projected distance of every gas cell to the origin.
      + `:dm_xy_distance`           -> Projected distance of every dark matter particle to the origin.
      + `:stellar_circularity`      -> Stellar circularity.
      + `:stellar_vcirc`            -> Stellar circular velocity.
      + `:stellar_age`              -> Stellar age.
      + `:sfr`                      -> The star formation rate of the last `AGE_RESOLUTION`.
      + `:ssfr`                     -> The specific star formation rate of the last `AGE_RESOLUTION`.
      + `:temperature`              -> Gas temperature, as log10(T / K).
  - `x_range::NTuple{2,<:Number}`: x axis range for the histogram grid.
  - `y_range::NTuple{2,<:Number}`: y axis range for the histogram grid.
  - `n_bins::Int64`: Number of bins per side of the grid.

# Returns

  - A tuple with three elements:

      + A vector with the x coordinates of the grid.
      + A vector with the y coordinates of the grid.
      + A matrix with the counts.
"""
function daScatterDensity(
    data_dict::Dict,
    x_quantity::Symbol,
    y_quantity::Symbol,
    x_range::NTuple{2,<:Number},
    y_range::NTuple{2,<:Number},
    n_bins::Int64,
)::Tuple{Vector{<:Number},Vector{<:Number},Matrix{Int64}}

    # Compute the values of the quantities for the x and y axis
    x_values = scatterQty(data_dict, x_quantity)
    y_values = scatterQty(data_dict, y_quantity)

    (
        length(x_values) == length(y_values) ||
        throw(ArgumentError("daScatterDensity: :$(x_quantity) and :$(y_quantity) \
        are incompatible quantities, they should be from the same type of cell/particle"))
    )

    # Compute the bin half width for each axis
    x_bin_h_width = 0.5 * (x_range[2] - x_range[1]) / n_bins
    y_bin_h_width = 0.5 * (y_range[2] - y_range[1]) / n_bins

    # Compute the center value of each bin for each axis
    x_axis = collect(range(x_range[1] + x_bin_h_width; length=n_bins, step=2 * x_bin_h_width))
    y_axis = collect(range(y_range[1] + y_bin_h_width; length=n_bins, step=2 * y_bin_h_width))

    # If any of the necessary quantities are missing return an empty histogram
    if any(isempty, [x_values, y_values])
        return x_axis, y_axis, fill(NaN, (n_bins, n_bins))
    end

    # Compute the 2D histogram
    counts = histogram2D(
        permutedims(hcat(x_values, y_values), (2, 1)),
        collect(range(x_range[1], x_range[2]; length=n_bins + 1)),
        collect(range(y_range[1], y_range[2]; length=n_bins + 1)),
    )

    # The transpose and reverse operation are to conform to the way heatmap! expect the matrix to be structured
    z_axis = reverse!(transpose(counts), dims=2)

    return x_axis, y_axis, z_axis

end

"""
    daVelocityField(
        data_dict::Dict,
        grid::SquareGrid,
        type_symbol::Symbol;
        <keyword arguments>
    )::Tuple{
        Vector{<:Unitful.Length},
        Vector{<:Unitful.Length},
        Matrix{<:Unitful.Velocity},
        Matrix{<:Unitful.Velocity},
    }

Compute a 2D mean velocity field.

!!! note

If the stellar masses and velocities can be found in `data_dict`, the velocity field is boosted with respect to the stellar center of mass.

# Arguments

  - `data_dict::Dict`: A dictionary with the following shape:

      + `:sim_data`          -> ::Simulation (see [`Simulation`](@ref)).
      + `:snap_data`         -> ::Snapshot (see [`Snapshot`](@ref)).
      + `:gc_data`           -> ::GroupCatalog (see [`GroupCatalog`](@ref)).
      + `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
      + `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
      + `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
      + ...
      + `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
      + `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
      + `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
      + ...

  - `grid::SquareGrid`: Square grid.
  - `type_symbol::Symbol`: For which cell/particle type the velocity field will be computed. The possibilities are the keys of [`ParticleIndex`](@ref).
  - `projection_plane::Symbol=:xy`: To which plane the cells/particles will be projected. The options are :xy, :xz, and :yz.

# Returns

  - A tuple with four elements:

      + A vector with the x coordinates of the grid.
      + A vector with the y coordinates of the grid.
      + A matrix with the mean velocity in the x direction for each bin.
      + A matrix with the mean velocity in the y direction for each bin.
"""
function daVelocityField(
    data_dict::Dict,
    grid::SquareGrid,
    type_symbol::Symbol;
    projection_plane::Symbol=:xy,
)::Tuple{
    Vector{<:Unitful.Length},
    Vector{<:Unitful.Length},
    Matrix{<:Unitful.Velocity},
    Matrix{<:Unitful.Velocity},
}

    positions = data_dict[type_symbol]["POS "]
    velocities = data_dict[type_symbol]["VEL "]

    # If any of the necessary quantities are missing return an empty velocity field
    if any(isempty, [positions, velocities])
        return grid.x_ticks, grid.y_ticks, zeros(size(grid.grid)), zeros(size(grid.grid))
    end

    # Project the cell/particles to the chosen plane
    if projection_plane == :xy
        pos_2D = positions[[1, 2], :]
    elseif projection_plane == :xz
        pos_2D = positions[[1, 3], :]
    elseif projection_plane == :yz
        pos_2D = positions[[2, 3], :]
    else
        throw(ArgumentError("daVelocityField: The argument `projection_plane` must be \
        :xy, :xz or :yz, but I got :$(projection_plane)"))
    end

    # Compute the components of the mean velocity
    vx = histogram2D(pos_2D, vec(velocities[1, :]), grid; total=false)
    vy = histogram2D(pos_2D, vec(velocities[2, :]), grid; total=false)

    return grid.x_ticks, grid.y_ticks, vx, vy

end

"""
    daIntegrateGalaxy(
        data_dict::Dict,
        x_quantity::Symbol,
        y_quantity::Symbol,
    )::NTuple{2,Vector{<:Number}}

Compute two global quantities of the simulation.

# Arguments

  - `data_dict::Dict`: A dictionary with the following shape:

      + `:sim_data`          -> ::Simulation (see [`Simulation`](@ref)).
      + `:snap_data`         -> ::Snapshot (see [`Snapshot`](@ref)).
      + `:gc_data`           -> ::GroupCatalog (see [`GroupCatalog`](@ref)).
      + `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
      + `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
      + `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
      + ...
      + `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
      + `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
      + `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
      + ...

  - `x_quantity::Symbol`: Quantity for the x axis. The possibilities are:

      + `:stellar_mass`           -> Stellar mass.
      + `:gas_mass`               -> Gas mass.
      + `:dm_mass`                -> Dark matter mass.
      + `:bh_mass`                -> Black hole mass.
      + `:molecular_mass`         -> Molecular hydrogen (H₂) mass.
      + `:atomic_mass`            -> Atomic hydrogen (HI) mass.
      + `:ionized_mass`           -> Ionized hydrogen (HII) mass.
      + `:neutral_mass`           -> Neutral hydrogen (HI + H₂) mass.
      + `:molecular_fraction`     -> Gas mass fraction of molecular hydrogen.
      + `:atomic_fraction`        -> Gas mass fraction of atomic hydrogen.
      + `:ionized_fraction`       -> Gas mass fraction of ionized hydrogen.
      + `:neutral_fraction`       -> Gas mass fraction of neutral hydrogen.
      + `:stellar_area_density`   -> Stellar area mass density, for a radius of `FILTER_R`.
      + `:gas_area_density`       -> Gas area mass density, for a radius of `FILTER_R`.
      + `:molecular_area_density` -> Molecular hydrogen area mass density, for a radius of `FILTER_R`.
      + `:atomic_area_density`    -> Atomic hydrogen area mass density, for a radius of `FILTER_R`.
      + `:ionized_area_density`   -> Ionized hydrogen area mass density, for a radius of `FILTER_R`.
      + `:neutral_area_density`   -> Neutral hydrogen area mass density, for a radius of `FILTER_R`.
      + `:sfr_area_density`       -> Star formation rate area density, for the last `AGE_RESOLUTION_ρ` and a radius of `FILTER_R`.
      + `:gas_metallicity`        -> Mass fraction of all elements above He in the gas (solar units).
      + `:stellar_metallicity`    -> Mass fraction of all elements above He in the stars (solar units).
      + `:X_gas_abundance`        -> Gas abundance of element X, as 12 + log10(X / H). The possibilities are the keys of [`ElementIndex`](@ref).
      + `:X_stellar_abundance`    -> Stellar abundance of element X, as 12 + log10(X / H). The possibilities are the keys of [`ElementIndex`](@ref).
      + `:stellar_specific_am`    -> Norm of the stellar specific angular momentum.
      + `:gas_specific_am`        -> Norm of the gas specific angular momentum.
      + `:dm_specific_am`         -> Norm of the dark matter specific angular momentum.
      + `:sfr`                    -> The star formation rate of the last `AGE_RESOLUTION`.
      + `:ssfr`                   -> The specific star formation rate of the last `AGE_RESOLUTION`.
      + `:scale_factor`           -> Scale factor.
      + `:redshift`               -> Redshift.
      + `:physical_time`          -> Physical time since the Big Bang.
      + `:lookback_time`          -> Physical time left to reach the last snapshot.
  - `y_quantity::Symbol`: Quantity for the y axis. The possibilities are:

      + `:stellar_mass`           -> Stellar mass.
      + `:gas_mass`               -> Gas mass.
      + `:dm_mass`                -> Dark matter mass.
      + `:bh_mass`                -> Black hole mass.
      + `:molecular_mass`         -> Molecular hydrogen (H₂) mass.
      + `:atomic_mass`            -> Atomic hydrogen (HI) mass.
      + `:ionized_mass`           -> Ionized hydrogen (HII) mass.
      + `:neutral_mass`           -> Neutral hydrogen (HI + H₂) mass.
      + `:molecular_fraction`     -> Gas mass fraction of molecular hydrogen.
      + `:atomic_fraction`        -> Gas mass fraction of atomic hydrogen.
      + `:ionized_fraction`       -> Gas mass fraction of ionized hydrogen.
      + `:neutral_fraction`       -> Gas mass fraction of neutral hydrogen.
      + `:stellar_area_density`   -> Stellar area mass density, for a radius of `FILTER_R`.
      + `:gas_area_density`       -> Gas area mass density, for a radius of `FILTER_R`.
      + `:molecular_area_density` -> Molecular hydrogen area mass density, for a radius of `FILTER_R`.
      + `:atomic_area_density`    -> Atomic hydrogen area mass density, for a radius of `FILTER_R`.
      + `:ionized_area_density`   -> Ionized hydrogen area mass density, for a radius of `FILTER_R`.
      + `:neutral_area_density`   -> Neutral hydrogen area mass density, for a radius of `FILTER_R`.
      + `:sfr_area_density`       -> Star formation rate area density, for the last `AGE_RESOLUTION_ρ` and a radius of `FILTER_R`.
      + `:gas_metallicity`        -> Mass fraction of all elements above He in the gas (solar units).
      + `:stellar_metallicity`    -> Mass fraction of all elements above He in the stars (solar units).
      + `:X_gas_abundance`        -> Gas abundance of element X, as 12 + log10(X / H). The possibilities are the keys of [`ElementIndex`](@ref).
      + `:X_stellar_abundance`    -> Stellar abundance of element X, as 12 + log10(X / H). The possibilities are the keys of [`ElementIndex`](@ref).
      + `:stellar_specific_am`    -> Norm of the stellar specific angular momentum.
      + `:gas_specific_am`        -> Norm of the gas specific angular momentum.
      + `:dm_specific_am`         -> Norm of the dark matter specific angular momentum.
      + `:sfr`                    -> The star formation rate of the last `AGE_RESOLUTION`.
      + `:ssfr`                   -> The specific star formation rate of the last `AGE_RESOLUTION`.
      + `:scale_factor`           -> Scale factor.
      + `:redshift`               -> Redshift.
      + `:physical_time`          -> Physical time since the Big Bang.
      + `:lookback_time`          -> Physical time left to reach the last snapshot.

# Returns

  - A tuple with two elements:

      + A single element vector with the value of `x_quantity`.
      + A single element vector with the value of `y_quantity`.
"""
function daIntegrateGalaxy(
    data_dict::Dict,
    x_quantity::Symbol,
    y_quantity::Symbol,
)::NTuple{2,Vector{<:Number}}

    return [integrateQty(data_dict, x_quantity)], [integrateQty(data_dict, y_quantity)]

end

"""
    daScatterGalaxy(
        data_dict::Dict,
        x_quantity::Symbol,
        y_quantity::Symbol,
    )::NTuple{2,Vector{<:Number}}

Compute two quantities for every cell/particle in the simulation.

# Arguments

  - `data_dict::Dict`: A dictionary with the following shape:

      + `:sim_data`          -> ::Simulation (see [`Simulation`](@ref)).
      + `:snap_data`         -> ::Snapshot (see [`Snapshot`](@ref)).
      + `:gc_data`           -> ::GroupCatalog (see [`GroupCatalog`](@ref)).
      + `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
      + `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
      + `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
      + ...
      + `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
      + `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
      + `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
      + ...

  - `x_quantity::Symbol`: Quantity for the x axis. The possibilities are:

      + `:stellar_mass`             -> Stellar mass.
      + `:gas_mass`                 -> Gas mass.
      + `:dm_mass`                  -> Dark matter mass.
      + `:bh_mass`                  -> Black hole mass.
      + `:molecular_mass`           -> Molecular hydrogen (H₂) mass.
      + `:atomic_mass`              -> Atomic hydrogen (HI) mass.
      + `:ionized_mass`             -> Ionized hydrogen (HII) mass.
      + `:neutral_mass`             -> Neutral hydrogen (HI + H₂) mass.
      + `:molecular_fraction`       -> Gas mass fraction of molecular hydrogen.
      + `:atomic_fraction`          -> Gas mass fraction of atomic hydrogen.
      + `:ionized_fraction`         -> Gas mass fraction of ionized hydrogen.
      + `:neutral_fraction`         -> Gas mass fraction of neutral hydrogen.
      + `:gas_mass_density`         -> Gas mass density.
      + `:gas_number_density`       -> Gas number density.
      + `:molecular_number_density` -> Molecular hydrogen number density.
      + `:atomic_number_density`    -> Atomic hydrogen number density.
      + `:ionized_number_density`   -> Ionized hydrogen number density.
      + `:neutral_number_density`   -> Neutral hydrogen number density.
      + `:gas_metallicity`          -> Mass fraction of all elements above He in the gas (solar units).
      + `:stellar_metallicity`      -> Mass fraction of all elements above He in the stars (solar units).
      + `:X_gas_abundance`          -> Gas abundance of element X, as 12 + log10(X / H). The possibilities are the keys of [`ElementIndex`](@ref).
      + `:X_stellar_abundance`      -> Stellar abundance of element X, as 12 + log10(X / H). The possibilities are the keys of [`ElementIndex`](@ref).
      + `:stellar_radial_distance`  -> Distance of every stellar particle to the origin.
      + `:gas_radial_distance`      -> Distance of every gas cell to the origin.
      + `:dm_radial_distance`       -> Distance of every dark matter particle to the origin.
      + `:stellar_xy_distance`      -> Projected distance of every stellar particle to the origin.
      + `:gas_xy_distance`          -> Projected distance of every gas cell to the origin.
      + `:dm_xy_distance`           -> Projected distance of every dark matter particle to the origin.
      + `:stellar_circularity`      -> Stellar circularity.
      + `:stellar_vcirc`            -> Stellar circular velocity.
      + `:stellar_age`              -> Stellar age.
      + `:sfr`                      -> The star formation rate of the last `AGE_RESOLUTION`.
      + `:ssfr`                     -> The specific star formation rate of the last `AGE_RESOLUTION`.
      + `:temperature`              -> Gas temperature, as log10(T / K).
  - `y_quantity::Symbol`: Quantity for the y axis. The possibilities are:

      + `:stellar_mass`             -> Stellar mass.
      + `:gas_mass`                 -> Gas mass.
      + `:dm_mass`                  -> Dark matter mass.
      + `:bh_mass`                  -> Black hole mass.
      + `:molecular_mass`           -> Molecular hydrogen (H₂) mass.
      + `:atomic_mass`              -> Atomic hydrogen (HI) mass.
      + `:ionized_mass`             -> Ionized hydrogen (HII) mass.
      + `:neutral_mass`             -> Neutral hydrogen (HI + H₂) mass.
      + `:molecular_fraction`       -> Gas mass fraction of molecular hydrogen.
      + `:atomic_fraction`          -> Gas mass fraction of atomic hydrogen.
      + `:ionized_fraction`         -> Gas mass fraction of ionized hydrogen.
      + `:neutral_fraction`         -> Gas mass fraction of neutral hydrogen.
      + `:gas_mass_density`         -> Gas mass density.
      + `:gas_number_density`       -> Gas number density.
      + `:molecular_number_density` -> Molecular hydrogen number density.
      + `:atomic_number_density`    -> Atomic hydrogen number density.
      + `:ionized_number_density`   -> Ionized hydrogen number density.
      + `:neutral_number_density`   -> Neutral hydrogen number density.
      + `:gas_metallicity`          -> Mass fraction of all elements above He in the gas (solar units).
      + `:stellar_metallicity`      -> Mass fraction of all elements above He in the stars (solar units).
      + `:X_gas_abundance`          -> Gas abundance of element X, as 12 + log10(X / H). The possibilities are the keys of [`ElementIndex`](@ref).
      + `:X_stellar_abundance`      -> Stellar abundance of element X, as 12 + log10(X / H). The possibilities are the keys of [`ElementIndex`](@ref).
      + `:stellar_radial_distance`  -> Distance of every stellar particle to the origin.
      + `:gas_radial_distance`      -> Distance of every gas cell to the origin.
      + `:dm_radial_distance`       -> Distance of every dark matter particle to the origin.
      + `:stellar_xy_distance`      -> Projected distance of every stellar particle to the origin.
      + `:gas_xy_distance`          -> Projected distance of every gas cell to the origin.
      + `:dm_xy_distance`           -> Projected distance of every dark matter particle to the origin.
      + `:stellar_circularity`      -> Stellar circularity.
      + `:stellar_vcirc`            -> Stellar circular velocity.
      + `:stellar_age`              -> Stellar age.
      + `:sfr`                      -> The star formation rate of the last `AGE_RESOLUTION`.
      + `:ssfr`                     -> The specific star formation rate of the last `AGE_RESOLUTION`.
      + `:temperature`              -> Gas temperature, as log10(T / K).

# Returns

  - A tuple with two elements:

      + A vector with the values of `x_quantity`.
      + A vector with the values of `y_quantity`.
"""
function daScatterGalaxy(
    data_dict::Dict,
    x_quantity::Symbol,
    y_quantity::Symbol,
)::NTuple{2,Vector{<:Number}}

    x_axis = scatterQty(data_dict, x_quantity)
    y_axis = scatterQty(data_dict, y_quantity)

    idx = sortperm(x_axis)

    (
        length(x_axis) == length(y_axis) ||
        throw(ArgumentError("daScatterGalaxy: :$(x_quantity) and :$(y_quantity) \
        are incompatible quantities, they should be from the same type of cell/particle"))
    )

    return x_axis[idx], y_axis[idx]

end

####################################################################################################
# For the `timeSeriesPlot` function in `./src/pipelines.jl`.
####################################################################################################
#
# A data analysis functions for `timeSeriesPlot` must take a `Simulation` struct, and return two 
# vectors. It should return `nothing` if the input data has some problem that prevents computation 
# (e.g. is empty).
#
# Expected signature:
#
#   da_function(sim_data, args...; kw_args...) -> (processed_data_x, processed_data_y)
#
# where:
#
#   - sim_data::Simulation, see the definition of `Simulation` in `./src/constants.jl`.
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

Compute the time series of two quantities.

# Arguments

  - `sim_data::Simulation`: Information about the simulation in a [`Simulation`](@ref) object.

  - `x_quantity::Symbol`: Quantity for the x axis. The possibilities are:

      + `:stellar_mass`           -> Stellar mass.
      + `:gas_mass`               -> Gas mass.
      + `:dm_mass`                -> Dark matter mass.
      + `:bh_mass`                -> Black hole mass.
      + `:molecular_mass`         -> Molecular hydrogen (H₂) mass.
      + `:atomic_mass`            -> Atomic hydrogen (HI) mass.
      + `:ionized_mass`           -> Ionized hydrogen (HII) mass.
      + `:neutral_mass`           -> Neutral hydrogen (HI + H₂) mass.
      + `:molecular_fraction`     -> Gas mass fraction of molecular hydrogen.
      + `:atomic_fraction`        -> Gas mass fraction of atomic hydrogen.
      + `:ionized_fraction`       -> Gas mass fraction of ionized hydrogen.
      + `:neutral_fraction`       -> Gas mass fraction of neutral hydrogen.
      + `:stellar_area_density`   -> Stellar area mass density, for a radius of `FILTER_R`.
      + `:gas_area_density`       -> Gas area mass density, for a radius of `FILTER_R`.
      + `:molecular_area_density` -> Molecular hydrogen area mass density, for a radius of `FILTER_R`.
      + `:atomic_area_density`    -> Atomic hydrogen area mass density, for a radius of `FILTER_R`.
      + `:ionized_area_density`   -> Ionized hydrogen area mass density, for a radius of `FILTER_R`.
      + `:neutral_area_density`   -> Neutral hydrogen area mass density, for a radius of `FILTER_R`.
      + `:sfr_area_density`       -> Star formation rate area density, for the last `AGE_RESOLUTION_ρ` and a radius of `FILTER_R`.
      + `:gas_metallicity`        -> Mass fraction of all elements above He in the gas (solar units).
      + `:stellar_metallicity`    -> Mass fraction of all elements above He in the stars (solar units).
      + `:X_gas_abundance`        -> Gas abundance of element X, as 12 + log10(X / H). The possibilities are the keys of [`ElementIndex`](@ref).
      + `:X_stellar_abundance`    -> Stellar abundance of element X, as 12 + log10(X / H). The possibilities are the keys of [`ElementIndex`](@ref).
      + `:stellar_specific_am`    -> Norm of the stellar specific angular momentum.
      + `:gas_specific_am`        -> Norm of the gas specific angular momentum.
      + `:dm_specific_am`         -> Norm of the dark matter specific angular momentum.
      + `:sfr`                    -> The star formation rate of the last `AGE_RESOLUTION`.
      + `:ssfr`                   -> The specific star formation rate of the last `AGE_RESOLUTION`.
      + `:scale_factor`           -> Scale factor.
      + `:redshift`               -> Redshift.
      + `:physical_time`          -> Physical time since the Big Bang.
      + `:lookback_time`          -> Physical time left to reach the last snapshot.
  - `y_quantity::Symbol`: Quantity for the y axis. The possibilities are:

      + `:stellar_mass`           -> Stellar mass.
      + `:gas_mass`               -> Gas mass.
      + `:dm_mass`                -> Dark matter mass.
      + `:bh_mass`                -> Black hole mass.
      + `:molecular_mass`         -> Molecular hydrogen (H₂) mass.
      + `:atomic_mass`            -> Atomic hydrogen (HI) mass.
      + `:ionized_mass`           -> Ionized hydrogen (HII) mass.
      + `:neutral_mass`           -> Neutral hydrogen (HI + H₂) mass.
      + `:molecular_fraction`     -> Gas mass fraction of molecular hydrogen.
      + `:atomic_fraction`        -> Gas mass fraction of atomic hydrogen.
      + `:ionized_fraction`       -> Gas mass fraction of ionized hydrogen.
      + `:neutral_fraction`       -> Gas mass fraction of neutral hydrogen.
      + `:stellar_area_density`   -> Stellar area mass density, for a radius of `FILTER_R`.
      + `:gas_area_density`       -> Gas area mass density, for a radius of `FILTER_R`.
      + `:molecular_area_density` -> Molecular hydrogen area mass density, for a radius of `FILTER_R`.
      + `:atomic_area_density`    -> Atomic hydrogen area mass density, for a radius of `FILTER_R`.
      + `:ionized_area_density`   -> Ionized hydrogen area mass density, for a radius of `FILTER_R`.
      + `:neutral_area_density`   -> Neutral hydrogen area mass density, for a radius of `FILTER_R`.
      + `:sfr_area_density`       -> Star formation rate area density, for the last `AGE_RESOLUTION_ρ` and a radius of `FILTER_R`.
      + `:gas_metallicity`        -> Mass fraction of all elements above He in the gas (solar units).
      + `:stellar_metallicity`    -> Mass fraction of all elements above He in the stars (solar units).
      + `:X_gas_abundance`        -> Gas abundance of element X, as 12 + log10(X / H). The possibilities are the keys of [`ElementIndex`](@ref).
      + `:X_stellar_abundance`    -> Stellar abundance of element X, as 12 + log10(X / H). The possibilities are the keys of [`ElementIndex`](@ref).
      + `:stellar_specific_am`    -> Norm of the stellar specific angular momentum.
      + `:gas_specific_am`        -> Norm of the gas specific angular momentum.
      + `:dm_specific_am`         -> Norm of the dark matter specific angular momentum.
      + `:sfr`                    -> The star formation rate of the last `AGE_RESOLUTION`.
      + `:ssfr`                   -> The specific star formation rate of the last `AGE_RESOLUTION`.
      + `:scale_factor`           -> Scale factor.
      + `:redshift`               -> Redshift.
      + `:physical_time`          -> Physical time since the Big Bang.
      + `:lookback_time`          -> Physical time left to reach the last snapshot.
  - `filter_mode::Symbol=:all`: Which cells/particles will be plotted, the options are:

      + `:all`             -> Plot every cell/particle within the simulation box.
      + `:halo`            -> Plot only the cells/particles that belong to the main halo.
      + `:subhalo`         -> Plot only the cells/particles that belong to the main subhalo.
      + `:sphere`          -> Plot only the cell/particle inside a sphere with radius `FILTER_R` (see `./src/constants.jl`).
      + `:stellar_subhalo` -> Plot only the cells/particles that belong to the main subhalo.
  - `smooth::Int64=0`: The result of [`integrateQty`](@ref) will be smooth out using `smooth` bins. Set it to 0 if you want no smoothing.
  - `scaling::Function=identity`: Function to scale the x-axis (only relevant if `smooth` != 0). The bins will be computed accordingly. The options are the scaling functions accepted by Makie.jl: log10, log2, log, sqrt, Makie.logit, Makie.Symlog10, Makie.pseudolog10, and identity.
  - `warnings::Bool=true`: If a warning will be given when there is missing data.

# Returns

  - A Tuple with two elements:

      + A Vector with the time series of `x_quantity`.
      + A Vector with the time series of `y_quantity`.
"""
function daEvolution(
    sim_data::Simulation,
    x_quantity::Symbol,
    y_quantity::Symbol;
    filter_mode::Symbol=:all,
    smooth::Int64=0,
    scaling::Function=identity,
    warnings::Bool=true,
)::NTuple{2,Vector{<:Number}}

    filter_function, translation, rotation, request = selectFilter(
        filter_mode,
        mergeRequests(plotParams(x_quantity).request, plotParams(y_quantity).request),
    )

    # Iterate over each snapshot in the slice
    iterator = eachrow(DataFrame(sim_data.table[sim_data.slice, :]))

    # Allocate memory
    x_axis = fill!(Vector{Number}(undef, length(iterator)), NaN)
    y_axis = fill!(Vector{Number}(undef, length(iterator)), NaN)

    @inbounds for (slice_index, sim_table_data) in pairs(iterator)

        global_index  = sim_table_data[1]
        scale_factor  = sim_table_data[3]
        redshift      = sim_table_data[4]
        physical_time = sim_table_data[5]
        lookback_time = sim_table_data[6]
        snapshot_path = sim_table_data[7]
        groupcat_path = sim_table_data[8]

        # Skip missing snapshots
        !ismissing(snapshot_path) || continue

        # Get the snapshot header
        snapshot_header = readSnapHeader(snapshot_path)

        # Get the group catalog header
        groupcat_header = readGroupCatHeader(groupcat_path; warnings)

        # Construct the metadata dictionary
        metadata = Dict(
            :sim_data => sim_data,
            :snap_data => Snapshot(
                snapshot_path,
                global_index,
                slice_index,
                physical_time,
                lookback_time,
                scale_factor,
                redshift,
                snapshot_header,
            ),
            :gc_data => GroupCatalog(groupcat_path, groupcat_header),
        )

        # Read the data in the snapshot
        data_dict = merge(
            metadata,
            readSnapshot(snapshot_path, request; warnings),
            readGroupCatalog(groupcat_path, snapshot_path, request; warnings),
        )

        # Filter the data
        filterData!(data_dict; filter_function)

        # Translate the data
        translateData!(data_dict, translation)

        # Rotate the data
        rotateData!(data_dict, rotation)

        # Compute the value for the x axis 
        x_axis[slice_index] = integrateQty(data_dict, x_quantity)

        # Compute the value for the y axis 
        y_axis[slice_index] = integrateQty(data_dict, y_quantity)

    end

    if iszero(smooth)
        return x_axis, y_axis
    else
        return smoothWindow(x_axis, y_axis, smooth; scaling)
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

  - `sim_data::Simulation`: Information about the simulation in a [`Simulation`](@ref) object.

  - `x_quantity::Symbol`: Quantity for the x axis. The possibilities are:

      + `:scale_factor`  -> Scale factor.
      + `:redshift`      -> Redshift.
      + `:physical_time` -> Physical time since the Big Bang.
      + `:lookback_time` -> Physical time left to reach the last snapshot.
  - `y_quantity::Symbol`: Quantity for the y axis. The possibilities are:

      + `:stellar_mass` -> Stellar mass.
      + `:sfr`          -> The star formation rate.
  - `warnings::Bool=true`: If a warning will be given when trying to use the scale factor or the redshift in the x axis for a non-cosmological simulation.

# Returns

  - A Tuple with two elements:

      + A Vector with the time series of `x_quantity`.
      + A Vector with the time series of `y_quantity`.
"""
function daSFRtxt(
    sim_data::Simulation,
    x_quantity::Symbol,
    y_quantity::Symbol;
    warnings::Bool=true,
)::NTuple{2,Vector{<:Number}}

    snapshot_paths = filter(!ismissing, sim_data.table[!, 7])

    (
        !isempty(snapshot_paths) ||
        throw(ArgumentError("daSFRtxt: I coudn't find any snapshots in $(sim_data.path), \
        and I need at least one for unit conversion"))
    )

    # Find the path to one snapshot 
    snapshot_path = first(snapshot_paths)

    # Read its header
    header = readSnapHeader(snapshot_path)

    # Read the data in the `sfr.txt` file
    sfr_txt_data = readSfrFile(joinpath(sim_data.path, SFR_REL_PATH), snapshot_path; warnings)
    time_ticks = sfr_txt_data[1]

    if x_quantity == :scale_factor

        (
            sim_data.cosmological ||
            @warn("daSFRtxt: For non-cosmological simulations `x_quantity` can only \
            be :physical_time")
        )

        x_axis = time_ticks

    elseif x_quantity == :redshift

        if sim_data.cosmological
            x_axis = (1.0 ./ time_ticks) .- 1.0
        else
            @warn("daSFRtxt: For non-cosmological simulations `x_quantity` \
            can only be :physical_time")
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
            physical_time = computeTime(time_ticks, header)
        else
            physical_time = time_ticks
        end

        x_axis = last(physical_time) .- physical_time

    else

        throw(ArgumentError("daSFRtxt: `x_quantity` can only be :scale_factor, \
        :redshift, :physical_time, or :lookback_time, but I got :$(x_quantity)"))

    end

    if y_quantity == :stellar_mass

        y_axis = sfr_txt_data[6]

    elseif y_quantity == :sfr

        y_axis = sfr_txt_data[3]

    else

        throw(ArgumentError("daSFRtxt: `y_quantity` can only be :stellar_mass or :sfr, \
        but I got :$(y_quantity)"))

    end

    return x_axis, y_axis

end
