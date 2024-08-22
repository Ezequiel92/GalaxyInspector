####################################################################################################
# Data analysis functions.
####################################################################################################

####################################################################################################
# Signature for the plotSnapshot function in ./src/pipelines.jl.
####################################################################################################
#
# A data analysis functions for plotSnapshot must take a dictionary with the following shape:
#
#   + :sim_data          -> ::Simulation (see the Simulation struct in ./src/constants.jl).
#   + :snap_data         -> ::Snapshot (see the Snapshot struct in ./src/constants.jl).
#   + :gc_data           -> ::GroupCatalog (see the GroupCatalog struct in ./src/constants.jl).
#   + cell/particle type -> (block -> data of block, block -> data of block, ...).
#   + cell/particle type -> (block -> data of block, block -> data of block, ...).
#   + cell/particle type -> (block -> data of block, block -> data of block, ...).
#   + ...
#   + groupcat type      -> (block -> data of block, block -> data of block, ...).
#   + groupcat type      -> (block -> data of block, block -> data of block, ...).
#   + groupcat type      -> (block -> data of block, block -> data of block, ...).
#   + ...
#
# and return one or more vectors or matrices. It should return `nothing` if the input data has
# some problem that prevents computation (e.g. is empty).
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
    daRotationCurve(
        data_dict::Dict,
        R::Unitful.Length;
        <keyword arguments>
    )::Tuple{Vector{<:Unitful.Length},Vector{<:Unitful.Velocity}}

Compute a rotation curve.

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
  - `R::Unitful.Length`: Maximum radius.
  - `filter_function::Function=filterNothing`: A function with the signature:

    `filter_function(data_dict) -> indices`

    where

      + `data_dict::Dict`: A dictionary with the following shape:

        * `:sim_data`          -> ::Simulation (see [`Simulation`](@ref)).
        * `:snap_data`         -> ::Snapshot (see [`Snapshot`](@ref)).
        * `:gc_data`           -> ::GroupCatalog (see [`GroupCatalog`](@ref)).
        * `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * ...
        * `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * ...
      + `indices::Dict`: A dictionary with the following shape:

        * `cell/particle type` -> idxs::IndexType
        * `cell/particle type` -> idxs::IndexType
        * `cell/particle type` -> idxs::IndexType
        * ...

# Returns

  - A tuple with two elements:

      + A vector with the distances to each star.
      + A vector with the circular velocity of each star.
"""
function daRotationCurve(
    data_dict::Dict,
    R::Unitful.Length;
    filter_function::Function=filterNothing,
)::Tuple{Vector{<:Unitful.Length},Vector{<:Unitful.Velocity}}

    filtered_dd = filterData(data_dict; filter_function)

    # Compute the circular velocities and the radial distances of each star
    r, vcirc = computeVcirc(filtered_dd)

    # Only leave the data within a sphere of radius `R`
    rangeCut!(r, vcirc, (0.0u"kpc", R))

    # Sort the arrays radialy
    idx = sortperm(r)

    return r[idx], vcirc[idx]

end

"""
    daKennicuttSchmidtLaw(
        data_dict::Dict,
        grid::CubicGrid,
        quantity::Symbol;
        <keyword arguments>
    )::Union{Tuple{Vector{<:SurfaceDensity},Vector{<:MassFlowDensity}},Nothing}

Compute the gas mass surface density and the SFR surface density, used in the Kennicutt-Schmidt law.

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
  - `grid::CubicGrid`: Cubic grid.
  - `quantity::Symbol=:molecular_mass`: Quantity for the x axis. The options are:

      + `:gas_mass`       -> Gas area mass density. This one will be plotted with the results of Kennicutt (1998).
      + `:molecular_mass` -> Molecular hydrogen area mass density. This one will be plotted with the results of Bigiel et al. (2008).
      + `:neutral_mass`   -> Neutral hydrogen area mass density. This one will be plotted with the results of Bigiel et al. (2008).
  - `type::Symbol=:cells`: If the density in the x axis will be calculated assuming gas as `:particles` or `:cells`.
  - `filter_function::Function=filterNothing`: A function with the signature:

    `filter_function(data_dict) -> indices`

    where

      + `data_dict::Dict`: A dictionary with the following shape:

        * `:sim_data`          -> ::Simulation (see [`Simulation`](@ref)).
        * `:snap_data`         -> ::Snapshot (see [`Snapshot`](@ref)).
        * `:gc_data`           -> ::GroupCatalog (see [`GroupCatalog`](@ref)).
        * `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * ...
        * `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * ...
      + `indices::Dict`: A dictionary with the following shape:

        * `cell/particle type` -> idxs::IndexType
        * `cell/particle type` -> idxs::IndexType
        * `cell/particle type` -> idxs::IndexType
        * ...

# Returns

  - A tuple with two elements:

      + A vector with log10(ΣH / M⊙ * kpc^-2).
      + A vector with log10(Σsfr / M⊙ * yr^-1 * kpc^-2).

    It returns `nothing` if any of the necessary quantities are missing.

# References

R. C. Kennicutt (1998). *The Global Schmidt Law in Star-forming Galaxies*. The Astrophysical Journal, **498(2)**, 541-552. [doi:10.1086/305588](https://doi.org/10.1086/305588)

F. Bigiel et al. (2008). *THE STAR FORMATION LAW IN NEARBY GALAXIES ON SUB-KPC SCALES*. The Astrophysical Journal, **136(6)**, 2846. [doi:10.1088/0004-6256/136/6/2846](https://doi.org/10.1088/0004-6256/136/6/2846)
"""
function daKennicuttSchmidtLaw(
    data_dict::Dict,
    grid::CubicGrid,
    quantity::Symbol;
    type::Symbol=:cells,
    filter_function::Function=filterNothing,
)::Union{NTuple{2,Vector{<:Float64}},Nothing}

    (
        quantity ∈ [:gas_mass, :molecular_mass, :neutral_mass] ||
        throw(ArgumentError("daKennicuttSchmidtLaw: `quantity` can only be :gas_mass, \
        :molecular_mass or :neutral_mass, but I got :$(quantity)"))
    )

    _, _, stellar_density = daDensity2DProjection(
        data_dict,
        grid,
        :stellar_mass,
        :particles;
        projection_plane=:xy,
        print_range=false,
        filter_function,
    )

    _, _, gas_density = daDensity2DProjection(
        data_dict,
        grid,
        quantity,
        type;
        projection_plane=:xy,
        print_range=false,
        filter_function,
    )

    x_axis = vec(gas_density)
    y_axis = vec(stellar_density)

    # Delete 0s and NaNs in the data vectors
    x_idxs = map(x -> isnan(x) || iszero(x), x_axis)
    y_idxs = map(x -> isnan(x) || iszero(x), y_axis)

    deleteat!(x_axis, x_idxs ∪ y_idxs)
    deleteat!(y_axis, x_idxs ∪ y_idxs)

    !any(isempty.([x_axis, y_axis])) || return nothing

    return x_axis, y_axis .- log10(ustrip(u"yr", AGE_RESOLUTION))

end

"""
    daMolla2015(
        data_dict::Dict,
        grid::CircularGrid,
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
  - `quantity::Symbol`: Quantity. The options are:

      + `:stellar_area_density`   -> Stellar area mass density.
      + `:molecular_area_density` -> Molecular hydrogen area mass density.
      + `:atomic_area_density`    -> Atomic hydrogen area mass density.
      + `:sfr_area_density`       -> Star formation rate area density, for the last `AGE_RESOLUTION`.
      + `:O_stellar_abundance`    -> Stellar abundance of oxygen, as ``12 + \\log_{10}(\\mathrm{O \\, / \\, H})``.
      + `:N_stellar_abundance`    -> Stellar abundance of nitrogen, as ``12 + \\log_{10}(\\mathrm{N \\, / \\, H})``.
      + `:C_stellar_abundance`    -> Stellar abundance of carbon, as ``12 + \\log_{10}(\\mathrm{C \\, / \\, H})``.
  - `filter_function::Function=filterNothing`: A function with the signature:

    `filter_function(data_dict) -> indices`

    where

      + `data_dict::Dict`: A dictionary with the following shape:

        * `:sim_data`          -> ::Simulation (see [`Simulation`](@ref)).
        * `:snap_data`         -> ::Snapshot (see [`Snapshot`](@ref)).
        * `:gc_data`           -> ::GroupCatalog (see [`GroupCatalog`](@ref)).
        * `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * ...
        * `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * ...
      + `indices::Dict`: A dictionary with the following shape:

        * `cell/particle type` -> idxs::IndexType
        * `cell/particle type` -> idxs::IndexType
        * `cell/particle type` -> idxs::IndexType
        * ...

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
    quantity::Symbol;
    filter_function::Function=filterNothing,
)::Union{
    Tuple{
        Vector{<:Unitful.Length},
        <:Union{Vector{<:SurfaceDensity},Vector{<:MassFlowDensity},Vector{Float64}},
    },
    Nothing,
}

    filtered_dd = filterData(data_dict; filter_function)

    if quantity == :stellar_area_density

        positions   = filtered_dd[:stars]["POS "]
        masses      = filtered_dd[:stars]["MASS"]
        norm_values = Number[]
        f           = identity
        density     = true

    elseif quantity == :sfr_area_density

        positions   = filtered_dd[:stars]["POS "]
        masses      = computeSFR(filtered_dd; age_resol=AGE_RESOLUTION)
        norm_values = Number[]
        f           = identity
        density     = true

    elseif quantity == :molecular_area_density

        positions   = filtered_dd[:gas]["POS "]
        masses      = computeMolecularMass(filtered_dd)
        norm_values = Number[]
        f           = identity
        density     = true

    elseif quantity == :atomic_area_density

        positions   = filtered_dd[:gas]["POS "]
        masses      = computeAtomicMass(filtered_dd)
        norm_values = Number[]
        f           = identity
        density     = true

    elseif quantity == :O_stellar_abundance

        positions   = filtered_dd[:stars]["POS "]
        masses      = computeElementMass(filtered_dd, :stars, :O) ./ ATOMIC_WEIGHTS[:O]
        norm_values = computeElementMass(filtered_dd, :stars, :H) ./ ATOMIC_WEIGHTS[:H]
        f           = x -> 12 .+ log10.(x)
        density     = false

    elseif quantity == :N_stellar_abundance

        positions   = filtered_dd[:stars]["POS "]
        masses      = computeElementMass(filtered_dd, :stars, :N) ./ ATOMIC_WEIGHTS[:N]
        norm_values = computeElementMass(filtered_dd, :stars, :H) ./ ATOMIC_WEIGHTS[:H]
        f           = x -> 12 .+ log10.(x)
        density     = false

    elseif quantity == :C_stellar_abundance

        positions   = filtered_dd[:stars]["POS "]
        masses      = computeElementMass(filtered_dd, :stars, :N) ./ ATOMIC_WEIGHTS[:N]
        norm_values = computeElementMass(filtered_dd, :stars, :H) ./ ATOMIC_WEIGHTS[:H]
        f           = x -> 12 .+ log10.(x)
        density     = false

    else

        throw(ArgumentError("daMolla2015: I don't recognize the quantity :$(quantity)"))

    end

    # Return `nothing` if any of the necessary quantities are missing
    !any(isempty, [positions, masses]) || return nothing

    density_profile = f(computeParticleProfile(positions, masses, grid; norm_values, total=true, density))

    return grid.grid, density_profile

end

"""
    daProfile(
        data_dict::Dict,
        quantity::Symbol,
        grid::CircularGrid;
        <keyword arguments>
    )::Union{Tuple{Vector{<:Unitful.Length},Vector{<:Number}},Nothing}

Compute a profile.

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
  - `quantity::Symbol`: Target quantity. The options are the same as for [`scatterQty`](@ref):

      + `:stellar_mass`               -> Stellar mass.
      + `:gas_mass`                   -> Gas mass.
      + `:hydrogen_mass`              -> Hydrogen mass.
      + `:dm_mass`                    -> Dark matter mass.
      + `:bh_mass`                    -> Black hole mass.
      + `:molecular_mass`             -> Molecular hydrogen (``\\mathrm{H_2}``) mass.
      + `:atomic_mass`                -> Atomic hydrogen (``\\mathrm{HI}``) mass.
      + `:ionized_mass`               -> Ionized hydrogen (``\\mathrm{HII}``) mass.
      + `:neutral_mass`               -> Neutral hydrogen (``\\mathrm{HI + H_2}``) mass.
      + `:molecular_fraction`         -> Gas mass fraction of molecular hydrogen.
      + `:atomic_fraction`            -> Gas mass fraction of atomic hydrogen.
      + `:ionized_fraction`           -> Gas mass fraction of ionized hydrogen.
      + `:neutral_fraction`           -> Gas mass fraction of neutral hydrogen.
      + `:molecular_neutral_fraction` -> Fraction of molecular hydrogen in the neutral gas.
      + `:gas_mass_density`           -> Gas mass density.
      + `:hydrogen_mass_density`      -> Hydrogen mass density.
      + `:gas_number_density`         -> Gas number density.
      + `:molecular_number_density`   -> Molecular hydrogen number density.
      + `:atomic_number_density`      -> Atomic hydrogen number density.
      + `:ionized_number_density`     -> Ionized hydrogen number density.
      + `:neutral_number_density`     -> Neutral hydrogen number density.
      + `:gas_metallicity`            -> Mass fraction of all elements above He in the gas (solar units).
      + `:stellar_metallicity`        -> Mass fraction of all elements above He in the stars (solar units).
      + `:X_gas_abundance`            -> Gas abundance of element ``\\mathrm{X}``, as ``12 + \\log_{10}(\\mathrm{X \\, / \\, H})``. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
      + `:X_stellar_abundance`        -> Stellar abundance of element ``\\mathrm{X}``, as ``12 + \\log_{10}(\\mathrm{X \\, / \\, H})``. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
      + `:stellar_radial_distance`    -> Distance of every stellar particle to the origin.
      + `:gas_radial_distance`        -> Distance of every gas cell to the origin.
      + `:dm_radial_distance`         -> Distance of every dark matter particle to the origin.
      + `:stellar_xy_distance`        -> Projected distance of every stellar particle to the origin.
      + `:gas_xy_distance`            -> Projected distance of every gas cell to the origin.
      + `:dm_xy_distance`             -> Projected distance of every dark matter particle to the origin.
      + `:stellar_circularity`        -> Stellar circularity.
      + `:stellar_vcirc`              -> Stellar circular velocity.
      + `:stellar_vradial`            -> Stellar radial speed.
      + `:stellar_vtangential`        -> Stellar tangential speed.
      + `:stellar_vzstar`             -> Stellar speed in the z direction, computed as ``v_z \\, \\mathrm{sign}(z)``.
      + `:stellar_age`                -> Stellar age.
      + `:sfr`                        -> The star formation rate.
      + `:ssfr`                       -> The specific star formation rate.
      + `:observational_sfr`          -> The star formation rate of the last `AGE_RESOLUTION`.
      + `:observational_ssfr`         -> The specific star formation rate of the last `AGE_RESOLUTION`.
      + `:temperature`                -> Gas temperature, as ``\\log_{10}(T \\, / \\, \\mathrm{K})``.
      + `:pressure`                   -> Gas pressure.
  - `grid::CircularGrid`: Circular grid.
  - `flat::Bool=true`: If the profile will be 2D, using rings, or 3D, using spherical shells.
  - `total::Bool=true`: If the sum (default) or the mean of `quantity` will be computed for each bin.
  - `cumulative::Bool=false`: If the profile will be accumulated or not.
  - `density::Bool=false`: If the profile will be of the density of `quantity`.
  - `fractions::Bool=false`: If a profile of the gas mass fractions will be calculated. It is only valid with `quantity` equal to :neutral_mass, :molecular_mass, :atomic_mass or :ionized_mass, and it forces `total` = true, `cumulative` = false, and `density` = false.
  - `filter_function::Function=filterNothing`: A function with the signature:

    `filter_function(data_dict) -> indices`

    where

      + `data_dict::Dict`: A dictionary with the following shape:

        * `:sim_data`          -> ::Simulation (see [`Simulation`](@ref)).
        * `:snap_data`         -> ::Snapshot (see [`Snapshot`](@ref)).
        * `:gc_data`           -> ::GroupCatalog (see [`GroupCatalog`](@ref)).
        * `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * ...
        * `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * ...
      + `indices::Dict`: A dictionary with the following shape:

        * `cell/particle type` -> idxs::IndexType
        * `cell/particle type` -> idxs::IndexType
        * `cell/particle type` -> idxs::IndexType
        * ...

# Returns

  - A tuple with two elements:

      + A vector with the position of each ring or spherical shells.
      + A vector with the value `quantity` in each each ring or spherical shells.

    It returns `nothing` if any of the necessary quantities are missing.
"""
function daProfile(
    data_dict::Dict,
    quantity::Symbol,
    grid::CircularGrid;
    flat::Bool=true,
    total::Bool=true,
    cumulative::Bool=false,
    density::Bool=false,
    fractions::Bool=false,
    filter_function::Function=filterNothing,
)::Union{Tuple{Vector{<:Unitful.Length},Vector{<:Number}},Nothing}

    filtered_dd = filterData(data_dict; filter_function)

    # Get the cell/particle type
    type = first(keys(plotParams(quantity).request))

    # Read the positions and values
    positions = filtered_dd[type]["POS "]
    values    = scatterQty(filtered_dd, quantity)

    n_pos = size(positions, 2)
    n_val = length(values)

    # Check consistency in the number of positions and values
    (
        n_pos == n_val || throw(ArgumentError("daProfile: `positions` and `values` should have \
        the same number of elements, but `length(positions)` = $(n_pos) != `length(values)` = \
        $(n_val). Check that the same cell/particle type was selected when using `plotParams` \
        and `scatterQty`."))
    )

    # Return `nothing` if any of the necessary quantities are missing
    !any(iszero, [n_pos, n_val]) || return nothing

    if fractions

        (
            quantity ∈ [:molecular_mass, :atomic_mass, :ionized_mass] ||
            throw(ArgumentError("daProfile: If `fractions``= true, quantity must be \
            :neutral_mass, :molecular_mass, :atomic_mass or :ionized_mass, but \
            I got `quantity` = :$(quantity)"))
        )

        total       = true
        cumulative  = false
        density     = false
        norm_values = scatterQty(filtered_dd, :hydrogen_mass)

    else

        norm_values = Number[]

    end

    profile = computeParticleProfile(
        positions,
        values,
        grid;
        norm_values,
        flat,
        total,
        cumulative,
        density,
    )

    return grid.grid, profile

end

"""
    daBandProfile(
        data_dict::Dict,
        quantity::Symbol,
        grid::CircularGrid;
        <keyword arguments>
    )::Union{
        Tuple{Vector{<:Unitful.Length},Vector{<:Number},Vector{<:Number},Vector{<:Number}},
        Tuple{Vector{<:Unitful.Length},Vector{<:Number},Vector{<:Number}},
        Nothing,
    }

Compute the profile og a mean quantity with error bars or bands.

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
  - `quantity::Symbol`: Target quantity. The possibilities are:

      + `:stellar_mass`               -> Stellar mass.
      + `:gas_mass`                   -> Gas mass.
      + `:hydrogen_mass`              -> Hydrogen mass.
      + `:dm_mass`                    -> Dark matter mass.
      + `:bh_mass`                    -> Black hole mass.
      + `:molecular_mass`             -> Molecular hydrogen (``\\mathrm{H_2}``) mass.
      + `:atomic_mass`                -> Atomic hydrogen (``\\mathrm{HI}``) mass.
      + `:ionized_mass`               -> Ionized hydrogen (``\\mathrm{HII}``) mass.
      + `:neutral_mass`               -> Neutral hydrogen (``\\mathrm{HI + H_2}``) mass.
      + `:molecular_fraction`         -> Gas mass fraction of molecular hydrogen.
      + `:atomic_fraction`            -> Gas mass fraction of atomic hydrogen.
      + `:ionized_fraction`           -> Gas mass fraction of ionized hydrogen.
      + `:neutral_fraction`           -> Gas mass fraction of neutral hydrogen.
      + `:molecular_neutral_fraction` -> Fraction of molecular hydrogen in the neutral gas.
      + `:gas_mass_density`           -> Gas mass density.
      + `:hydrogen_mass_density`      -> Hydrogen mass density.
      + `:gas_number_density`         -> Gas number density.
      + `:molecular_number_density`   -> Molecular hydrogen number density.
      + `:atomic_number_density`      -> Atomic hydrogen number density.
      + `:ionized_number_density`     -> Ionized hydrogen number density.
      + `:neutral_number_density`     -> Neutral hydrogen number density.
      + `:gas_metallicity`            -> Mass fraction of all elements above He in the gas (solar units).
      + `:stellar_metallicity`        -> Mass fraction of all elements above He in the stars (solar units).
      + `:X_gas_abundance`            -> Gas abundance of element ``\\mathrm{X}``, as ``12 + \\log_{10}(\\mathrm{X \\, / \\, H})``. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
      + `:X_stellar_abundance`        -> Stellar abundance of element ``\\mathrm{X}``, as ``12 + \\log_{10}(\\mathrm{X \\, / \\, H})``. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
      + `:stellar_radial_distance`    -> Distance of every stellar particle to the origin.
      + `:gas_radial_distance`        -> Distance of every gas cell to the origin.
      + `:dm_radial_distance`         -> Distance of every dark matter particle to the origin.
      + `:stellar_xy_distance`        -> Projected distance of every stellar particle to the origin.
      + `:gas_xy_distance`            -> Projected distance of every gas cell to the origin.
      + `:dm_xy_distance`             -> Projected distance of every dark matter particle to the origin.
      + `:stellar_circularity`        -> Stellar circularity.
      + `:stellar_vcirc`              -> Stellar circular velocity.
      + `:stellar_vradial`            -> Stellar radial speed.
      + `:stellar_vtangential`        -> Stellar tangential speed.
      + `:stellar_vzstar`             -> Stellar speed in the z direction, computed as ``v_z \\, \\mathrm{sign}(z)``.
      + `:stellar_age`                -> Stellar age.
      + `:sfr`                        -> The star formation rate.
      + `:ssfr`                       -> The specific star formation rate.
      + `:observational_sfr`          -> The star formation rate of the last `AGE_RESOLUTION`.
      + `:observational_ssfr`         -> The specific star formation rate of the last `AGE_RESOLUTION`.
      + `:temperature`                -> Gas temperature, as ``\\log_{10}(T \\, / \\, \\mathrm{K})``.
      + `:pressure`                   -> Gas pressure.
  - `grid::CircularGrid`: Circular grid.
  - `flat::Bool=true`: If the profile will be 2D, using rings, or 3D, using spherical shells.
  - `error_bar::Bool=false`: If the returned values will be compatible with `errorbars!` or with `band!` (default).
  - `filter_function::Function=filterNothing`: A function with the signature:

    `filter_function(data_dict) -> indices`

    where

      + `data_dict::Dict`: A dictionary with the following shape:

        * `:sim_data`          -> ::Simulation (see [`Simulation`](@ref)).
        * `:snap_data`         -> ::Snapshot (see [`Snapshot`](@ref)).
        * `:gc_data`           -> ::GroupCatalog (see [`GroupCatalog`](@ref)).
        * `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * ...
        * `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * ...
      + `indices::Dict`: A dictionary with the following shape:

        * `cell/particle type` -> idxs::IndexType
        * `cell/particle type` -> idxs::IndexType
        * `cell/particle type` -> idxs::IndexType
        * ...

# Returns

  - A tuple with two elements:

      + A vector with the position of each ring or spherical shells.
      + A vector with the value `quantity` in each each ring or spherical shells.

    It returns `nothing` if any of the necessary quantities are missing.
"""
function daBandProfile(
    data_dict::Dict,
    quantity::Symbol,
    grid::CircularGrid;
    flat::Bool=true,
    error_bar::Bool=false,
    filter_function::Function=filterNothing,
)::Union{
    Tuple{Vector{<:Unitful.Length},Vector{<:Number},Vector{<:Number},Vector{<:Number}},
    Tuple{Vector{<:Unitful.Length},Vector{<:Number},Vector{<:Number}},
    Nothing,
}

    filtered_dd = filterData(data_dict; filter_function)

    # Get the cell/particle type
    type = first(keys(plotParams(quantity).request))

    # Read the positions and values
    positions = filtered_dd[type]["POS "]
    values    = scatterQty(filtered_dd, quantity)

    n_pos = size(positions, 2)
    n_val = length(values)

    # Check consistency in the number of positions and values
    (
        n_pos == n_val || throw(ArgumentError("daBandProfile: `positions` and `values` should have \
        the same number of elements, but `length(positions)` = $(n_pos) != `length(values)` = \
        $(n_val). Check that the same cell/particle type was selected when using `plotParams` \
        and `scatterQty`."))
    )

    # Return `nothing` if any of the necessary quantities are missing
    !any(iszero, [n_pos, n_val]) || return nothing

    mean, std = computeParticleBandProfile(positions, values, grid; flat)

    !error_bar || return grid.grid, mean, std, std

    return grid.grid, mean .- std, mean .+ std

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
  - `quantity::Symbol=:sfr`: Target quantity. The options are:

      + `:sfr`                 -> The star formation rate.
      + `:ssfr`                -> The specific star formation rate.
      + `:stellar_mass`        -> Stellar mass.
      + `:stellar_metallicity` -> Mass fraction of all elements above He in the stars (solar units).
  - `n_bins::Int=100`: Number of bins (time intervals).
  - `filter_function::Function=filterNothing`: A function with the signature:

    `filter_function(data_dict) -> indices`

    where

      + `data_dict::Dict`: A dictionary with the following shape:

        * `:sim_data`          -> ::Simulation (see [`Simulation`](@ref)).
        * `:snap_data`         -> ::Snapshot (see [`Snapshot`](@ref)).
        * `:gc_data`           -> ::GroupCatalog (see [`GroupCatalog`](@ref)).
        * `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * ...
        * `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * ...
      + `indices::Dict`: A dictionary with the following shape:

        * `cell/particle type` -> idxs::IndexType
        * `cell/particle type` -> idxs::IndexType
        * `cell/particle type` -> idxs::IndexType
        * ...

# Returns

  - A tuple with two elements:

      + A vector with the physical times.
      + A vector with the values of `quantity` at each time.
"""
function daStellarHistory(
    data_dict::Dict;
    quantity::Symbol=:sfr,
    n_bins::Int=100,
    filter_function::Function=filterNothing,
)::Union{Tuple{Vector{<:Unitful.Time},Vector{<:Number}},Nothing}

    filtered_dd = filterData(data_dict; filter_function)

    birth_ticks   = filtered_dd[:stars]["GAGE"]
    masses        = filtered_dd[:stars]["MASS"]

    # Return `nothing` if any of the necessary quantities are missing
    !any(isempty, [birth_ticks, masses]) || return nothing

    # Compute the stellar birth dates
    if filtered_dd[:sim_data].cosmological
        # Go from scale factor to physical time
        birth_times = computeTime(birth_ticks, filtered_dd[:snap_data].header)
    else
        birth_times = birth_ticks
    end

    # Compute the birth time range
    min, max = extrema(birth_times)

    # Compute the total stellar mass in each time bin
    grid = CircularGrid(max, n_bins; shift=min)

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

        metallicities = computeMetalMass(data_dict, :stars)

        # Return `nothing` if any of the necessary quantities are missing
        !isempty(metallicities) || return nothing

        stellar_metallicities = histogram1D(birth_times, metallicities, grid; empty_nan=false)

        Z = cumsum(stellar_metallicities) ./ cumsum(stellar_masses)
        y_axis = Z ./ SOLAR_METALLICITY

    else

        throw(ArgumentError("daStellarHistory: `quantity` can only be :sfr, :ssfr, \
        :stellar_metallicity or :stellar_mass, but I got :$(quantity)"))

    end

    return x_axis, y_axis

end

"""
    daLineHistogram(
        data_dict::Dict,
        quantity::Symbol,
        grid::LinearGrid;
        <keyword arguments>
    )::Union{Tuple{Vector{<:Number},Vector{<:Number}},Nothing}

Compute a 1D histogram of a given `quantity`, normalized to the maximum number of counts.

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
  - `quantity::Symbol`: The possibilities are:

      + `:stellar_mass`               -> Stellar mass.
      + `:gas_mass`                   -> Gas mass.
      + `:hydrogen_mass`              -> Hydrogen mass.
      + `:dm_mass`                    -> Dark matter mass.
      + `:bh_mass`                    -> Black hole mass.
      + `:molecular_mass`             -> Molecular hydrogen (``\\mathrm{H_2}``) mass.
      + `:atomic_mass`                -> Atomic hydrogen (``\\mathrm{HI}``) mass.
      + `:ionized_mass`               -> Ionized hydrogen (``\\mathrm{HII}``) mass.
      + `:neutral_mass`               -> Neutral hydrogen (``\\mathrm{HI + H_2}``) mass.
      + `:molecular_fraction`         -> Gas mass fraction of molecular hydrogen.
      + `:atomic_fraction`            -> Gas mass fraction of atomic hydrogen.
      + `:ionized_fraction`           -> Gas mass fraction of ionized hydrogen.
      + `:neutral_fraction`           -> Gas mass fraction of neutral hydrogen.
      + `:molecular_neutral_fraction` -> Fraction of molecular hydrogen in the neutral gas.
      + `:gas_mass_density`           -> Gas mass density.
      + `:hydrogen_mass_density`      -> Hydrogen mass density.
      + `:gas_number_density`         -> Gas number density.
      + `:molecular_number_density`   -> Molecular hydrogen number density.
      + `:atomic_number_density`      -> Atomic hydrogen number density.
      + `:ionized_number_density`     -> Ionized hydrogen number density.
      + `:neutral_number_density`     -> Neutral hydrogen number density.
      + `:gas_metallicity`            -> Mass fraction of all elements above He in the gas (solar units).
      + `:stellar_metallicity`        -> Mass fraction of all elements above He in the stars (solar units).
      + `:X_gas_abundance`            -> Gas abundance of element ``\\mathrm{X}``, as ``12 + \\log_{10}(\\mathrm{X \\, / \\, H})``. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
      + `:X_stellar_abundance`        -> Stellar abundance of element ``\\mathrm{X}``, as ``12 + \\log_{10}(\\mathrm{X \\, / \\, H})``. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
      + `:stellar_radial_distance`    -> Distance of every stellar particle to the origin.
      + `:gas_radial_distance`        -> Distance of every gas cell to the origin.
      + `:dm_radial_distance`         -> Distance of every dark matter particle to the origin.
      + `:stellar_xy_distance`        -> Projected distance of every stellar particle to the origin.
      + `:gas_xy_distance`            -> Projected distance of every gas cell to the origin.
      + `:dm_xy_distance`             -> Projected distance of every dark matter particle to the origin.
      + `:stellar_circularity`        -> Stellar circularity.
      + `:stellar_vcirc`              -> Stellar circular velocity.
      + `:stellar_vradial`            -> Stellar radial speed.
      + `:stellar_vtangential`        -> Stellar tangential speed.
      + `:stellar_vzstar`             -> Stellar speed in the z direction, computed as ``v_z \\, \\mathrm{sign}(z)``.
      + `:stellar_age`                -> Stellar age.
      + `:sfr`                        -> The star formation rate.
      + `:ssfr`                       -> The specific star formation rate.
      + `:observational_sfr`          -> The star formation rate of the last `AGE_RESOLUTION`.
      + `:observational_ssfr`         -> The specific star formation rate of the last `AGE_RESOLUTION`.
      + `:temperature`                -> Gas temperature, as ``\\log_{10}(T \\, / \\, \\mathrm{K})``.
      + `:pressure`                   -> Gas pressure.
  - `grid::LinearGrid`: Linear grid.
  - `filter_function::Function=filterNothing`: A function with the signature:

    `filter_function(data_dict) -> indices`

    where

      + `data_dict::Dict`: A dictionary with the following shape:

        * `:sim_data`          -> ::Simulation (see [`Simulation`](@ref)).
        * `:snap_data`         -> ::Snapshot (see [`Snapshot`](@ref)).
        * `:gc_data`           -> ::GroupCatalog (see [`GroupCatalog`](@ref)).
        * `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * ...
        * `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * ...
      + `indices::Dict`: A dictionary with the following shape:

        * `cell/particle type` -> idxs::IndexType
        * `cell/particle type` -> idxs::IndexType
        * `cell/particle type` -> idxs::IndexType
        * ...
  - `norm::Int=0`: Number of count that will be use to normalize the histogram. If left as 0, the histogram will be normalize with the maximum bin count.

# Returns

  - A tuple with two elements:

      + A vector with the value corresponding to each bin.
      + A vector with the counts, normalized to the maximum value.
"""
function daLineHistogram(
    data_dict::Dict,
    quantity::Symbol,
    grid::LinearGrid;
    filter_function::Function=filterNothing,
    norm::Int=0,
)::Union{Tuple{Vector{<:Number},Vector{<:Number}},Nothing}

    filtered_dd = filterData(data_dict; filter_function)

    # Compute the values
    values = scatterQty(filtered_dd, quantity)

    !isempty(values) || return nothing

    # Compute the quantity histogram
    counts = histogram1D(values, grid)

    # Normalize the counts
    norm_counts = isPositive(norm) ? counts ./ norm : counts ./ maximum(counts)

    return grid.grid, norm_counts

end

"""
    daStellarCircHistogram(
        data_dict::Dict,
        quantity::Symbol,
        grid::LinearGrid;
        <keyword arguments>
    )::Union{Tuple{Vector{<:Number},Vector{<:Number}},Nothing}

Compute a 1D histogram of a given `quantity`.

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
  - `quantity::Symbol`: The possibilities are:

      + `:stellar_circularity`        -> Stellar circularity.
      + `:stellar_vcirc`              -> Stellar circular velocity.
  - `grid::LinearGrid`: Linear grid.
  - `filter_function::Function=filterNothing`: A function with the signature:

    `filter_function(data_dict) -> indices`

    where

      + `data_dict::Dict`: A dictionary with the following shape:

        * `:sim_data`          -> ::Simulation (see [`Simulation`](@ref)).
        * `:snap_data`         -> ::Snapshot (see [`Snapshot`](@ref)).
        * `:gc_data`           -> ::GroupCatalog (see [`GroupCatalog`](@ref)).
        * `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * ...
        * `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * ...
      + `indices::Dict`: A dictionary with the following shape:

        * `cell/particle type` -> idxs::IndexType
        * `cell/particle type` -> idxs::IndexType
        * `cell/particle type` -> idxs::IndexType
        * ...
  - `norm::Int=1`: Number of count that will be use to normalize the histogram.

# Returns

  - A tuple with two elements:

      + A vector with the value corresponding to each bin.
      + A vector with the counts, normalized to the maximum value.
"""
function daStellarCircHistogram(
    data_dict::Dict,
    quantity::Symbol,
    grid::LinearGrid;
    filter_function::Function=filterNothing,
    norm::Int=1,
)::Union{Tuple{Vector{<:Number},Vector{<:Number}},Nothing}

    (
        quantity ∈ [:stellar_circularity, :stellar_vcirc] ||
        throw(ArgumentError("daCircularityHistogram: `quantity` can only be :stellar_vcirc \
        or :stellar_circularity, but I got :$(quantity)"))
    )

    (
        isPositive(norm) ||
        throw(ArgumentError("daCircularityHistogram: `norm` must be a positive interger, \
        but I got norm = $(norm)"))
    )

    # Compute the indices of the target stars
    stellar_idxs = filter_function(data_dict)[:stars]

    # Compute the values
    values = scatterQty(data_dict, quantity)

    # Return nothing if there are no stars before of after filtering
    !(isempty(values) || isempty(values[stellar_idxs])) || return nothing

    # Compute the quantity histogram
    counts = histogram1D(values[stellar_idxs], grid)

    return grid.grid, counts ./ norm

end

"""
    daDensity2DProjection(
        data_dict::Dict,
        grid::CubicGrid,
        quantity::Symbol,
        type::Symbol;
        <keyword arguments>
    )::Tuple{Vector{<:Unitful.Length},Vector{<:Unitful.Length},Matrix{Float64}}

Project the 3D density field to a given plane.

If the source of the field are particles a simple 2D histogram is used. If the source of the field are Voronoi cells the density of the cells that cross the line of sight of each pixel are added up.

!!! note

    By default, ``\\mathrm{M_\\odot \\, kpc^{-2}}`` is used as unit of density, so the output will be ``\\log_{10}(\\rho \\, [\\mathrm{M_\\odot \\, kpc^{-2}}])``.

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
  - `grid::CubicGrid`: Cubic grid.
  - `quantity::Symbol`: For which quantity the density will be calculated. The options are:

      + `:stellar_mass`   -> Stellar mass.
      + `:gas_mass`       -> Gas mass.
      + `:hydrogen_mass`  -> Hydrogen mass.
      + `:dm_mass`        -> Dark matter mass.
      + `:bh_mass`        -> Black hole mass.
      + `:molecular_mass` -> Molecular hydrogen (``\\mathrm{H_2}``) mass.
      + `:atomic_mass`    -> Atomic hydrogen (``\\mathrm{HI}``) mass.
      + `:ionized_mass`   -> Ionized hydrogen (``\\mathrm{HII}``) mass.
      + `:neutral_mass`   -> Neutral hydrogen (``\\mathrm{HI + H_2}``) mass.
  - `type::Symbol`: If the source of the field are `:particles` or `:cells`.
  - `projection_plane::Symbol=:xy`: To which plane the cells/particles will be projected. The options are `:xy`, `:xz`, and `:yz`.
  - `print_range::Bool=false`: Print an info block detailing the logarithmic density range.
  - `filter_function::Function=filterNothing`: A function with the signature:

    `filter_function(data_dict) -> indices`

    where

      + `data_dict::Dict`: A dictionary with the following shape:

        * `:sim_data`          -> ::Simulation (see [`Simulation`](@ref)).
        * `:snap_data`         -> ::Snapshot (see [`Snapshot`](@ref)).
        * `:gc_data`           -> ::GroupCatalog (see [`GroupCatalog`](@ref)).
        * `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * ...
        * `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * ...
      + `indices::Dict`: A dictionary with the following shape:

        * `cell/particle type` -> idxs::IndexType
        * `cell/particle type` -> idxs::IndexType
        * `cell/particle type` -> idxs::IndexType
        * ...

# Returns

  - A tuple with three elements:

      + A vector with the x coordinates of the grid.
      + A vector with the y coordinates of the grid.
      + A matrix with the values of the logarithmic density at each point of the 2D grid.
"""
function daDensity2DProjection(
    data_dict::Dict,
    grid::CubicGrid,
    quantity::Symbol,
    type::Symbol;
    projection_plane::Symbol=:xy,
    print_range::Bool=false,
    filter_function::Function=filterNothing,
)::Tuple{Vector{<:Unitful.Length},Vector{<:Unitful.Length},Matrix{Float64}}

    filtered_dd = filterData(data_dict; filter_function)

    # Set the cell/particle type
    if quantity ∈ [
        :gas_mass,
        :hydrogen_mass,
        :molecular_mass,
        :atomic_mass,
        :ionized_mass,
        :neutral_mass,
    ]
        component = :gas
    elseif quantity == :stellar_mass
        component = :stars
    elseif quantity == :dm_mass
        component = :halo
    elseif quantity == :bh_mass
        component = :black_hole
    else
        throw(ArgumentError("daDensity2DProjection: I don't recognize the quantity :$(quantity)"))
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

    # Load the positions
    positions = filtered_dd[component]["POS "]

    # Compute the masses
    masses = scatterQty(filtered_dd, quantity)

    # If any of the necessary quantities are missing return an empty density field
    if any(isempty, [masses, positions])
        return grid.x_ticks, grid.y_ticks, fill(NaN, (grid.n_bins, grid.n_bins))
    end

    # Set the units
    m_unit = u"Msun"
    l_unit = u"kpc"

    if type == :cells

        # Compute the volume of each cell
        cell_volumes = filtered_dd[component]["MASS"] ./ filtered_dd[component]["RHO "]

        # Compute the densities of the target quantity
        densities = ustrip.(m_unit * l_unit^-3, masses ./ cell_volumes)

        # Load the volume and area of the voxels
        voxel_volume = ustrip(l_unit^3, grid.bin_volume)
        voxel_area   = ustrip(l_unit^2, grid.bin_area)

        # Allocate memory
        physical_grid = Matrix{Float64}(undef, 3, grid.n_bins^3)

        # Compute the tree for a nearest neighbor search
        kdtree = KDTree(ustrip.(l_unit, positions))

        # Reshape the grid to conform to the way `nn` expect the matrix to be structured
        @inbounds for i in eachindex(grid.grid)
            physical_grid[1, i] = ustrip(l_unit, grid.grid[i][1])
            physical_grid[2, i] = ustrip(l_unit, grid.grid[i][2])
            physical_grid[3, i] = ustrip(l_unit, grid.grid[i][3])
        end

        # Find the nearest neighbor to each point in the grid
        idxs, _ = nn(kdtree, physical_grid)

        # Allocate memory
        mass_grid = similar(grid.grid, Float64)

        # Compute the mass in each voxel
        @inbounds for i in eachindex(grid.grid)
            mass_grid[i] = densities[idxs[i]] * voxel_volume
        end

        # Project the grid to the given plane
        if projection_plane == :xy
            density = dropdims(sum(mass_grid; dims=3) ./ voxel_area; dims=3)
        elseif projection_plane == :xz
            density = transpose(dropdims(sum(mass_grid; dims=2) ./ voxel_area; dims=2))
        elseif projection_plane == :yz
            density = dropdims(sum(mass_grid; dims=1) ./ voxel_area; dims=1)
        else
            throw(ArgumentError("daDensity2DProjection: The argument `projection_plane` must be \
            :xy, :xz or :yz, but I got :$(projection_plane)"))
        end

    elseif type == :particles

        # Project the particles to the given plane
        if projection_plane == :xy
            pos_2D = positions[[1, 2], :]
        elseif projection_plane == :xz
            pos_2D = positions[[1, 3], :]
        elseif projection_plane == :yz
            pos_2D = positions[[2, 3], :]
        else
            throw(ArgumentError("daDensity2DProjection: The argument `projection_plane` must be \
            :xy, :xz or :yz, but I got :$(projection_plane)"))
        end

        # Compute the 2D histogram
        density = ustrip.(
            m_unit * l_unit^-2,
            histogram2D(pos_2D, masses, flattenGrid(grid); empty_nan=false) ./ grid.bin_area,
        )

    else

        throw(ArgumentError("daDensity2DProjection: The argument `type` must be :cells or \
        :particles, but I got :$(type)"))

    end

    # Set bins with a value of 0 to NaN
    replace!(x -> iszero(x) ? NaN : x, density)

    # Apply log10 to enhance the contrast
    values = log10.(density ./ physical_factor)

    if print_range

        # Compute the mininimum and maximum values of density
        min_max = isempty(values) ? (NaN, NaN) : extrema(filter(!isnan, values))

        # Print the density range
        @info(
            "\nDensity range \
            \n  Simulation: $(basename(filtered_dd[:sim_data].path)) \
            \n  Snapshot:   $(filtered_dd[:snap_data].global_index) \
            \n  Quantity:   $(quantity) \
            \n  Plane:      $(projection_plane) \
            \n  log₁₀(ρ [$(m_unit * l_unit^-2)]): $(min_max)\n\n"
        )

    end

    # The transpose and reverse operation are to conform to the way `heatmap!` expect the matrix to be structured
    z_axis = reverse!(transpose(values), dims=2)

    return grid.x_ticks, grid.y_ticks, z_axis

end

"""
    daMetallicity2DProjection(
        data_dict::Dict,
        grid::CubicGrid,
        component::Symbol;
        <keyword arguments>
    )::Tuple{Vector{<:Unitful.Length},Vector{<:Unitful.Length},Matrix{Float64}}

Project the 3D metallicity field to a given plane.

If `component` = :stars, the metallicity in each pixel is the mean value for the stars within that pixel. If `component` = :gas, the metallicity in each pixel is the total metal mass divided by the total gas mass, in the column given by taht pixel.

!!! note

    By default, the matllicity is given in solar units.

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
  - `grid::CubicGrid`: Cubic grid.
  - `component::Symbol`: Target cell/particle type. It can be either `:stars` or `:gas`.
  - `projection_plane::Symbol=:xy`: To which plane the cells/particles will be projected. The options are `:xy`, `:xz`, and `:yz`.
  - `print_range::Bool=false`: Print an info block detailing the logarithmic metallicity range.
  - `filter_function::Function=filterNothing`: A function with the signature:

    `filter_function(data_dict) -> indices`

    where

      + `data_dict::Dict`: A dictionary with the following shape:

        * `:sim_data`          -> ::Simulation (see [`Simulation`](@ref)).
        * `:snap_data`         -> ::Snapshot (see [`Snapshot`](@ref)).
        * `:gc_data`           -> ::GroupCatalog (see [`GroupCatalog`](@ref)).
        * `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * ...
        * `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * ...
      + `indices::Dict`: A dictionary with the following shape:

        * `cell/particle type` -> idxs::IndexType
        * `cell/particle type` -> idxs::IndexType
        * `cell/particle type` -> idxs::IndexType
        * ...

# Returns

  - A tuple with three elements:

      + A vector with the x coordinates of the grid.
      + A vector with the y coordinates of the grid.
      + A matrix with the values of the logarithmic metallicty at each point of the 2D grid.
"""
function daMetallicity2DProjection(
    data_dict::Dict,
    grid::CubicGrid,
    component::Symbol;
    projection_plane::Symbol=:xy,
    print_range::Bool=false,
    filter_function::Function=filterNothing,
)::Tuple{Vector{<:Unitful.Length},Vector{<:Unitful.Length},Matrix{Float64}}

    filtered_dd = filterData(data_dict; filter_function)

    (
        component ∈ [:gas, :stars] ||
        throw(ArgumentError("daMetallicity2DProjection: I don't recognize the component \
        :$(component)"))
    )

    # Load the positions
    positions = filtered_dd[component]["POS "]

    # If the necessary quantities are missing return an empty density field
    if isempty(positions)
        return grid.x_ticks, grid.y_ticks, fill(NaN, (grid.n_bins, grid.n_bins))
    end

    if component == :gas

        # Compute the volume of each cell
        gas_density = filtered_dd[:gas]["RHO "]
        gas_masses  = filtered_dd[component]["MASS"]
        gas_volumes = gas_masses ./ gas_density

        # Compute the metal densities
        metal_densities = computeMetalMass(filtered_dd, :gas) ./ gas_volumes

        # Allocate memory
        physical_grid = Matrix{Float64}(undef, 3, grid.n_bins^3)

        # Reshape the grid to conform to the way `nn` expect the matrix to be structured
        @inbounds for i in eachindex(grid.grid)
            physical_grid[1, i] = ustrip(u"kpc", grid.grid[i][1])
            physical_grid[2, i] = ustrip(u"kpc", grid.grid[i][2])
            physical_grid[3, i] = ustrip(u"kpc", grid.grid[i][3])
        end

        # Compute the tree for a nearest neighbor search
        kdtree = KDTree(ustrip.(u"kpc", positions))

        # Find the nearest neighbor to each point in the grid
        idxs, _ = nn(kdtree, physical_grid)

        # Allocate memory
        gas_grid   = similar(grid.grid, Float64)
        metal_grid = similar(grid.grid, Float64)

        # Compute the mass in each voxel
        @inbounds for i in eachindex(grid.grid)
            gas_grid[i]   = ustrip.(u"Msun", gas_density[idxs[i]] * grid.bin_volume)
            metal_grid[i] = ustrip.(u"Msun", metal_densities[idxs[i]] * grid.bin_volume)
        end

        # Project the grid to the chosen plane
        if projection_plane == :xy

            metallicity = dropdims(sum(metal_grid; dims=3) ./ sum(gas_grid; dims=3); dims=3)

        elseif projection_plane == :xz

            metallicity = transpose(
                dropdims(sum(metal_grid; dims=2) ./ sum(gas_grid; dims=2); dims=2),
            )

        elseif projection_plane == :yz

            metallicity = dropdims(sum(metal_grid; dims=1) ./ sum(gas_grid; dims=1); dims=1)

        else

            throw(ArgumentError("daMetallicity2DProjection: The argument `projection_plane` must \
            be :xy, :xz or :yz, but I got :$(projection_plane)"))

        end

    else

        # Project the particles to the chosen plane
        if projection_plane == :xy
            pos_2D = positions[[1, 2], :]
        elseif projection_plane == :xz
            pos_2D = positions[[1, 3], :]
        elseif projection_plane == :yz
            pos_2D = positions[[2, 3], :]
        else
            throw(ArgumentError("daMetallicity2DProjection: The argument `projection_plane` must \
            be :xy, :xz or :yz, but I got :$(projection_plane)"))
        end

        metallicities = uconvert.(
            Unitful.NoUnits,
            computeMetalMass(filtered_dd, :stars) ./ filtered_dd[:stars]["MASS"],
        )

        # Compute the 2D histogram
        metallicity = histogram2D(
            pos_2D,
            metallicities,
            flattenGrid(grid);
            total=false,
            empty_nan=false,
        )

    end

    # Set bins with a value of 0 to NaN
    replace!(x -> iszero(x) ? NaN : x, metallicity)

    # Apply log10 to enhance the contrast
    values = log10.(metallicity ./ SOLAR_METALLICITY)

    if print_range

        # Compute the min and max values of metallicity
        min_max = isempty(values) ? (NaN, NaN) : extrema(filter(!isnan, values))

        # Print the metallicity range
        @info(
            "\nMetallicity range \
            \n  Simulation:      $(basename(filtered_dd[:sim_data].path)) \
            \n  Snapshot:        $(filtered_dd[:snap_data].global_index) \
            \n  P/C type:        $(component) \
            \n  Plane:           $(projection_plane) \
            \n  log₁₀(Z [Z⊙]):  $(min_max)\n\n"
        )

    end

    # The transpose and reverse operation are to conform to the way heatmap! expect the matrix to be structured
    z_axis = reverse!(transpose(values), dims=2)

    return grid.x_ticks, grid.y_ticks, z_axis

end

"""
    daTemperature2DProjection(
        data_dict::Dict,
        grid::CubicGrid;
        <keyword arguments>
    )::Tuple{Vector{<:Unitful.Length},Vector{<:Unitful.Length},Matrix{Float64}}

Project the 3D temperature field to a given plane.

!!! note

    By default, ``K`` is used as unit of temperature, so the output will be ``\\log_{10}(T \\, [\\mathrm{K}])``.

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
  - `grid::CubicGrid`: Cubic grid.
  - `projection_plane::Symbol=:xy`: To which plane the cells will be projected. The options are `:xy`, `:xz`, and `:yz`.
  - `print_range::Bool=false`: Print an info block detailing the logarithmic temperature range.
  - `filter_function::Function=filterNothing`: A function with the signature:

    `filter_function(data_dict) -> indices`

    where

      + `data_dict::Dict`: A dictionary with the following shape:

        * `:sim_data`          -> ::Simulation (see [`Simulation`](@ref)).
        * `:snap_data`         -> ::Snapshot (see [`Snapshot`](@ref)).
        * `:gc_data`           -> ::GroupCatalog (see [`GroupCatalog`](@ref)).
        * `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * ...
        * `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * ...
      + `indices::Dict`: A dictionary with the following shape:

        * `cell/particle type` -> idxs::IndexType
        * `cell/particle type` -> idxs::IndexType
        * `cell/particle type` -> idxs::IndexType
        * ...

# Returns

  - A tuple with three elements:

      + A vector with the x coordinates of the grid.
      + A vector with the y coordinates of the grid.
      + A matrix with the values of temperature at each grid point.
"""
function daTemperature2DProjection(
    data_dict::Dict,
    grid::CubicGrid;
    projection_plane::Symbol=:xy,
    print_range::Bool=false,
    filter_function::Function=filterNothing,
)::Tuple{Vector{<:Unitful.Length},Vector{<:Unitful.Length},Matrix{Float64}}

    filtered_dd = filterData(data_dict; filter_function)

    # Load the temperatures
    temperatures = ustrip.(u"K", filtered_dd[:gas]["TEMP"])

    # Load the positions
    positions = filtered_dd[:gas]["POS "]

    # If any of the necessary quantities are missing return an empty temperature field
    if any(isempty, [temperatures, positions])
        return grid.x_ticks, grid.y_ticks, fill(NaN, (grid.n_bins, grid.n_bins))
    end

    # Allocate memory
    physical_grid = Matrix{Float64}(undef, 3, grid.n_bins^3)

    # Reshape the grid to conform to the way `nn` expect the matrix to be structured
    @inbounds for i in eachindex(grid.grid)
        physical_grid[1, i] = ustrip(u"kpc", grid.grid[i][1])
        physical_grid[2, i] = ustrip(u"kpc", grid.grid[i][2])
        physical_grid[3, i] = ustrip(u"kpc", grid.grid[i][3])
    end

    # Compute the tree for a nearest neighbor search
    kdtree = KDTree(ustrip.(u"kpc", positions))

    # Find the nearest neighbor to each point in the grid
    idxs, _ = nn(kdtree, physical_grid)

    # Allocate memory
    temperature_grid = similar(grid.grid, Float64)

    # Compute the temperature of each voxel
    @inbounds for i in eachindex(grid.grid)
        temperature_grid[i] = temperatures[idxs[i]]
    end

    # Project the grid to the chosen plane
    if projection_plane == :xy
        temperature = dropdims(sum(temperature_grid; dims=3) ./ grid.n_bins; dims=3)
    elseif projection_plane == :xz
        temperature = transpose(dropdims(sum(temperature_grid; dims=2) ./ grid.n_bins; dims=2))
    elseif projection_plane == :yz
        temperature = dropdims(sum(temperature_grid; dims=1) ./ grid.n_bins; dims=1)
    else
        throw(ArgumentError("daTemperature2DProjection: The argument `projection_plane` must be \
        :xy, :xz or :yz, but I got :$(projection_plane)"))
    end

    # Apply log10 to enhance the contrast
    values = log10.(temperature)

    if print_range

        # Print the temperature range
        @info(
            "\nDensity range \
            \n  Simulation: $(basename(filtered_dd[:sim_data].path)) \
            \n  Snapshot:   $(filtered_dd[:snap_data].global_index) \
            \n  Plane:      $(projection_plane) \
            \nlog₁₀(T [K]) = $(extrema(values))\n\n"
        )

    end

    # The transpose and reverse operation are to conform to the way heatmap! expect the matrix to be structured
    z_axis = reverse!(transpose(values), dims=2)

    return grid.x_ticks, grid.y_ticks, z_axis

end

"""
    daScatterDensity(
        data_dict::Dict,
        x_quantity::Symbol,
        y_quantity::Symbol;
        <keyword arguments>
    )::Tuple{Vector{<:Number},Vector{<:Number},Matrix{Float64}}

Turn a scatter plot into a 2D histogram.

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
  - `x_quantity::Symbol`: Quantity for the x axis. The options are:

      + `:stellar_mass`               -> Stellar mass.
      + `:gas_mass`                   -> Gas mass.
      + `:hydrogen_mass`              -> Hydrogen mass.
      + `:dm_mass`                    -> Dark matter mass.
      + `:bh_mass`                    -> Black hole mass.
      + `:molecular_mass`             -> Molecular hydrogen (``\\mathrm{H_2}``) mass.
      + `:atomic_mass`                -> Atomic hydrogen (``\\mathrm{HI}``) mass.
      + `:ionized_mass`               -> Ionized hydrogen (``\\mathrm{HII}``) mass.
      + `:neutral_mass`               -> Neutral hydrogen (``\\mathrm{HI + H_2}``) mass.
      + `:molecular_fraction`         -> Gas mass fraction of molecular hydrogen.
      + `:atomic_fraction`            -> Gas mass fraction of atomic hydrogen.
      + `:ionized_fraction`           -> Gas mass fraction of ionized hydrogen.
      + `:neutral_fraction`           -> Gas mass fraction of neutral hydrogen.
      + `:molecular_neutral_fraction` -> Fraction of molecular hydrogen in the neutral gas.
      + `:gas_mass_density`           -> Gas mass density.
      + `:hydrogen_mass_density`      -> Hydrogen mass density.
      + `:gas_number_density`         -> Gas number density.
      + `:molecular_number_density`   -> Molecular hydrogen number density.
      + `:atomic_number_density`      -> Atomic hydrogen number density.
      + `:ionized_number_density`     -> Ionized hydrogen number density.
      + `:neutral_number_density`     -> Neutral hydrogen number density.
      + `:gas_metallicity`            -> Mass fraction of all elements above He in the gas (solar units).
      + `:stellar_metallicity`        -> Mass fraction of all elements above He in the stars (solar units).
      + `:X_gas_abundance`            -> Gas abundance of element ``\\mathrm{X}``, as ``12 + \\log_{10}(\\mathrm{X \\, / \\, H})``. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
      + `:X_stellar_abundance`        -> Stellar abundance of element ``\\mathrm{X}``, as ``12 + \\log_{10}(\\mathrm{X \\, / \\, H})``. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
      + `:stellar_radial_distance`    -> Distance of every stellar particle to the origin.
      + `:gas_radial_distance`        -> Distance of every gas cell to the origin.
      + `:dm_radial_distance`         -> Distance of every dark matter particle to the origin.
      + `:stellar_xy_distance`        -> Projected distance of every stellar particle to the origin.
      + `:gas_xy_distance`            -> Projected distance of every gas cell to the origin.
      + `:dm_xy_distance`             -> Projected distance of every dark matter particle to the origin.
      + `:stellar_circularity`        -> Stellar circularity.
      + `:stellar_vcirc`              -> Stellar circular velocity.
      + `:stellar_vradial`            -> Stellar radial speed.
      + `:stellar_vtangential`        -> Stellar tangential speed.
      + `:stellar_vzstar`             -> Stellar speed in the z direction, computed as ``v_z \\, \\mathrm{sign}(z)``.
      + `:stellar_age`                -> Stellar age.
      + `:sfr`                        -> The star formation rate.
      + `:ssfr`                       -> The specific star formation rate.
      + `:observational_sfr`          -> The star formation rate of the last `AGE_RESOLUTION`.
      + `:observational_ssfr`         -> The specific star formation rate of the last `AGE_RESOLUTION`.
      + `:temperature`                -> Gas temperature, as ``\\log_{10}(T \\, / \\, \\mathrm{K})``.
      + `:pressure`                   -> Gas pressure.
  - `y_quantity::Symbol`: Quantity for the y axis. The options are:

      + `:stellar_mass`               -> Stellar mass.
      + `:gas_mass`                   -> Gas mass.
      + `:hydrogen_mass`              -> Hydrogen mass.
      + `:dm_mass`                    -> Dark matter mass.
      + `:bh_mass`                    -> Black hole mass.
      + `:molecular_mass`             -> Molecular hydrogen (``\\mathrm{H_2}``) mass.
      + `:atomic_mass`                -> Atomic hydrogen (``\\mathrm{HI}``) mass.
      + `:ionized_mass`               -> Ionized hydrogen (``\\mathrm{HII}``) mass.
      + `:neutral_mass`               -> Neutral hydrogen (``\\mathrm{HI + H_2}``) mass.
      + `:molecular_fraction`         -> Gas mass fraction of molecular hydrogen.
      + `:atomic_fraction`            -> Gas mass fraction of atomic hydrogen.
      + `:ionized_fraction`           -> Gas mass fraction of ionized hydrogen.
      + `:neutral_fraction`           -> Gas mass fraction of neutral hydrogen.
      + `:molecular_neutral_fraction` -> Fraction of molecular hydrogen in the neutral gas.
      + `:gas_mass_density`           -> Gas mass density.
      + `:hydrogen_mass_density`      -> Hydrogen mass density.
      + `:gas_number_density`         -> Gas number density.
      + `:molecular_number_density`   -> Molecular hydrogen number density.
      + `:atomic_number_density`      -> Atomic hydrogen number density.
      + `:ionized_number_density`     -> Ionized hydrogen number density.
      + `:neutral_number_density`     -> Neutral hydrogen number density.
      + `:gas_metallicity`            -> Mass fraction of all elements above He in the gas (solar units).
      + `:stellar_metallicity`        -> Mass fraction of all elements above He in the stars (solar units).
      + `:X_gas_abundance`            -> Gas abundance of element ``\\mathrm{X}``, as ``12 + \\log_{10}(\\mathrm{X \\, / \\, H})``. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
      + `:X_stellar_abundance`        -> Stellar abundance of element ``\\mathrm{X}``, as ``12 + \\log_{10}(\\mathrm{X \\, / \\, H})``. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
      + `:stellar_radial_distance`    -> Distance of every stellar particle to the origin.
      + `:gas_radial_distance`        -> Distance of every gas cell to the origin.
      + `:dm_radial_distance`         -> Distance of every dark matter particle to the origin.
      + `:stellar_xy_distance`        -> Projected distance of every stellar particle to the origin.
      + `:gas_xy_distance`            -> Projected distance of every gas cell to the origin.
      + `:dm_xy_distance`             -> Projected distance of every dark matter particle to the origin.
      + `:stellar_circularity`        -> Stellar circularity.
      + `:stellar_vcirc`              -> Stellar circular velocity.
      + `:stellar_vradial`            -> Stellar radial speed.
      + `:stellar_vtangential`        -> Stellar tangential speed.
      + `:stellar_vzstar`             -> Stellar speed in the z direction, computed as ``v_z \\, \\mathrm{sign}(z)``.
      + `:stellar_age`                -> Stellar age.
      + `:sfr`                        -> The star formation rate.
      + `:ssfr`                       -> The specific star formation rate.
      + `:observational_sfr`          -> The star formation rate of the last `AGE_RESOLUTION`.
      + `:observational_ssfr`         -> The specific star formation rate of the last `AGE_RESOLUTION`.
      + `:temperature`                -> Gas temperature, as ``\\log_{10}(T \\, / \\, \\mathrm{K})``.
      + `:pressure`                   -> Gas pressure.
  - `x_range::Union{NTuple{2,<:Number},Nothing}=nothing`: x axis range for the histogram grid. If set to `nothing`, the extrema of the values will be used.
  - `y_range::Union{NTuple{2,<:Number},Nothing}=nothing`: y axis range for the histogram grid. If set to `nothing`, the extrema of the values will be used.
  - `x_log::Union{Unitful.Units,Nothing}=nothing`: Desired unit of `x_quantity`, if you want to use log10(`x_quantity`) for the x axis.
  - `y_log::Union{Unitful.Units,Nothing}=nothing`: Desired unit of `y_quantity`, if you want to use log10(`y_quantity`) for the y axis.
  - `n_bins::Int=100`: Number of bins per side of the grid.
  - `filter_function::Function=filterNothing`: A function with the signature:

    `filter_function(data_dict) -> indices`

    where

      + `data_dict::Dict`: A dictionary with the following shape:

        * `:sim_data`          -> ::Simulation (see [`Simulation`](@ref)).
        * `:snap_data`         -> ::Snapshot (see [`Snapshot`](@ref)).
        * `:gc_data`           -> ::GroupCatalog (see [`GroupCatalog`](@ref)).
        * `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * ...
        * `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * ...
      + `indices::Dict`: A dictionary with the following shape:

        * `cell/particle type` -> idxs::IndexType
        * `cell/particle type` -> idxs::IndexType
        * `cell/particle type` -> idxs::IndexType
        * ...

# Returns

  - A tuple with three elements:

      + A vector with the x coordinates of the grid.
      + A vector with the y coordinates of the grid.
      + A matrix with the counts.
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

    filtered_dd = filterData(data_dict; filter_function)

    # Compute the values of the quantities
    x_values = scatterQty(filtered_dd, x_quantity)
    y_values = scatterQty(filtered_dd, y_quantity)

    # If any of the necessary quantities are missing return an empty histogram
    if any(isempty, [x_values, y_values])
        return 1:n_bins, 1:n_bins, fill(NaN, (n_bins, n_bins))
    end

    (
        length(x_values) == length(y_values) ||
        throw(ArgumentError("daScatterDensity: :$(x_quantity) and :$(y_quantity) have a diferent \
        number of values. They should be the same"))
    )

    # If requested, apply log10 to the x axis data, ignoring 0 values
    if !isnothing(x_log)
        null_x_idxs = findall(iszero, x_values)
        x_values    = log10.(deleteat!(ustrip.(x_log, x_values), null_x_idxs))
        y_values    = deleteat!(y_values, null_x_idxs)
    end

    # If requested, apply log10 to the y axis data, ignoring 0 values
    if !isnothing(y_log)
        null_y_idxs = findall(iszero, y_values)
        x_values    = deleteat!(x_values, null_y_idxs)
        y_values    = log10.(deleteat!(ustrip.(y_log, y_values), null_y_idxs))
    end

    # If there is no range specified, use the extrema of the x values
    if isnothing(x_range)
        x_range = extrema(x_values)
    end

    # If there is no range specified, use the extrema of the y values
    if isnothing(y_range)
        y_range = extrema(y_values)
    end

    # Compute the bin half width for each axis
    x_bin_h_width = 0.5 * (x_range[2] - x_range[1]) / n_bins
    y_bin_h_width = 0.5 * (y_range[2] - y_range[1]) / n_bins

    # Compute the center value of each bin for each axis
    x_axis = collect(range(x_range[1] + x_bin_h_width; length=n_bins, step=2 * x_bin_h_width))
    y_axis = collect(range(y_range[1] + y_bin_h_width; length=n_bins, step=2 * y_bin_h_width))

    # Compute the 2D histogram
    counts = Float64.(histogram2D(
        permutedims(hcat(x_values, y_values), (2, 1)),
        collect(range(x_range[1], x_range[2]; length=n_bins + 1)),
        collect(range(y_range[1], y_range[2]; length=n_bins + 1)),
    ))

    # Set bins with a value of 0 to NaN
    replace!(x -> iszero(x) ? NaN : x, counts)

    # The transpose and reverse operation are to conform to the way heatmap! expect the matrix to be structured,
    # and log10 is used to enhance the contrast
    z_axis = reverse!(transpose(log10.(counts)), dims=2)

    return x_axis, y_axis, z_axis

end

"""
    daScatterWeightedDensity(
        data_dict::Dict,
        x_quantity::Symbol,
        y_quantity::Symbol,
        z_quantity::Symbol,
        z_unit::Uniful.Units;
        <keyword arguments>
    )::Tuple{Vector{<:Number},Vector{<:Number},Matrix{Float64}}

Turn a scatter plot into a 2D histogram, weighted by `z_quantity`.

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
  - `x_quantity::Symbol`: Quantity for the x axis. The options are:

      + `:stellar_mass`               -> Stellar mass.
      + `:gas_mass`                   -> Gas mass.
      + `:hydrogen_mass`              -> Hydrogen mass.
      + `:dm_mass`                    -> Dark matter mass.
      + `:bh_mass`                    -> Black hole mass.
      + `:molecular_mass`             -> Molecular hydrogen (``\\mathrm{H_2}``) mass.
      + `:atomic_mass`                -> Atomic hydrogen (``\\mathrm{HI}``) mass.
      + `:ionized_mass`               -> Ionized hydrogen (``\\mathrm{HII}``) mass.
      + `:neutral_mass`               -> Neutral hydrogen (``\\mathrm{HI + H_2}``) mass.
      + `:molecular_fraction`         -> Gas mass fraction of molecular hydrogen.
      + `:atomic_fraction`            -> Gas mass fraction of atomic hydrogen.
      + `:ionized_fraction`           -> Gas mass fraction of ionized hydrogen.
      + `:neutral_fraction`           -> Gas mass fraction of neutral hydrogen.
      + `:molecular_neutral_fraction` -> Fraction of molecular hydrogen in the neutral gas.
      + `:gas_mass_density`           -> Gas mass density.
      + `:hydrogen_mass_density`      -> Hydrogen mass density.
      + `:gas_number_density`         -> Gas number density.
      + `:molecular_number_density`   -> Molecular hydrogen number density.
      + `:atomic_number_density`      -> Atomic hydrogen number density.
      + `:ionized_number_density`     -> Ionized hydrogen number density.
      + `:neutral_number_density`     -> Neutral hydrogen number density.
      + `:gas_metallicity`            -> Mass fraction of all elements above He in the gas (solar units).
      + `:stellar_metallicity`        -> Mass fraction of all elements above He in the stars (solar units).
      + `:X_gas_abundance`            -> Gas abundance of element ``\\mathrm{X}``, as ``12 + \\log_{10}(\\mathrm{X \\, / \\, H})``. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
      + `:X_stellar_abundance`        -> Stellar abundance of element ``\\mathrm{X}``, as ``12 + \\log_{10}(\\mathrm{X \\, / \\, H})``. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
      + `:stellar_radial_distance`    -> Distance of every stellar particle to the origin.
      + `:gas_radial_distance`        -> Distance of every gas cell to the origin.
      + `:dm_radial_distance`         -> Distance of every dark matter particle to the origin.
      + `:stellar_xy_distance`        -> Projected distance of every stellar particle to the origin.
      + `:gas_xy_distance`            -> Projected distance of every gas cell to the origin.
      + `:dm_xy_distance`             -> Projected distance of every dark matter particle to the origin.
      + `:stellar_circularity`        -> Stellar circularity.
      + `:stellar_vcirc`              -> Stellar circular velocity.
      + `:stellar_vradial`            -> Stellar radial speed.
      + `:stellar_vtangential`        -> Stellar tangential speed.
      + `:stellar_vzstar`             -> Stellar speed in the z direction, computed as ``v_z \\, \\mathrm{sign}(z)``.
      + `:stellar_age`                -> Stellar age.
      + `:sfr`                        -> The star formation rate.
      + `:ssfr`                       -> The specific star formation rate.
      + `:observational_sfr`          -> The star formation rate of the last `AGE_RESOLUTION`.
      + `:observational_ssfr`         -> The specific star formation rate of the last `AGE_RESOLUTION`.
      + `:temperature`                -> Gas temperature, as ``\\log_{10}(T \\, / \\, \\mathrm{K})``.
      + `:pressure`                   -> Gas pressure.
  - `y_quantity::Symbol`: Quantity for the y axis. The options are:

      + `:stellar_mass`               -> Stellar mass.
      + `:gas_mass`                   -> Gas mass.
      + `:hydrogen_mass`              -> Hydrogen mass.
      + `:dm_mass`                    -> Dark matter mass.
      + `:bh_mass`                    -> Black hole mass.
      + `:molecular_mass`             -> Molecular hydrogen (``\\mathrm{H_2}``) mass.
      + `:atomic_mass`                -> Atomic hydrogen (``\\mathrm{HI}``) mass.
      + `:ionized_mass`               -> Ionized hydrogen (``\\mathrm{HII}``) mass.
      + `:neutral_mass`               -> Neutral hydrogen (``\\mathrm{HI + H_2}``) mass.
      + `:molecular_fraction`         -> Gas mass fraction of molecular hydrogen.
      + `:atomic_fraction`            -> Gas mass fraction of atomic hydrogen.
      + `:ionized_fraction`           -> Gas mass fraction of ionized hydrogen.
      + `:neutral_fraction`           -> Gas mass fraction of neutral hydrogen.
      + `:molecular_neutral_fraction` -> Fraction of molecular hydrogen in the neutral gas.
      + `:gas_mass_density`           -> Gas mass density.
      + `:hydrogen_mass_density`      -> Hydrogen mass density.
      + `:gas_number_density`         -> Gas number density.
      + `:molecular_number_density`   -> Molecular hydrogen number density.
      + `:atomic_number_density`      -> Atomic hydrogen number density.
      + `:ionized_number_density`     -> Ionized hydrogen number density.
      + `:neutral_number_density`     -> Neutral hydrogen number density.
      + `:gas_metallicity`            -> Mass fraction of all elements above He in the gas (solar units).
      + `:stellar_metallicity`        -> Mass fraction of all elements above He in the stars (solar units).
      + `:X_gas_abundance`            -> Gas abundance of element ``\\mathrm{X}``, as ``12 + \\log_{10}(\\mathrm{X \\, / \\, H})``. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
      + `:X_stellar_abundance`        -> Stellar abundance of element ``\\mathrm{X}``, as ``12 + \\log_{10}(\\mathrm{X \\, / \\, H})``. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
      + `:stellar_radial_distance`    -> Distance of every stellar particle to the origin.
      + `:gas_radial_distance`        -> Distance of every gas cell to the origin.
      + `:dm_radial_distance`         -> Distance of every dark matter particle to the origin.
      + `:stellar_xy_distance`        -> Projected distance of every stellar particle to the origin.
      + `:gas_xy_distance`            -> Projected distance of every gas cell to the origin.
      + `:dm_xy_distance`             -> Projected distance of every dark matter particle to the origin.
      + `:stellar_circularity`        -> Stellar circularity.
      + `:stellar_vcirc`              -> Stellar circular velocity.
      + `:stellar_vradial`            -> Stellar radial speed.
      + `:stellar_vtangential`        -> Stellar tangential speed.
      + `:stellar_vzstar`             -> Stellar speed in the z direction, computed as ``v_z \\, \\mathrm{sign}(z)``.
      + `:stellar_age`                -> Stellar age.
      + `:sfr`                        -> The star formation rate.
      + `:ssfr`                       -> The specific star formation rate.
      + `:observational_sfr`          -> The star formation rate of the last `AGE_RESOLUTION`.
      + `:observational_ssfr`         -> The specific star formation rate of the last `AGE_RESOLUTION`.
      + `:temperature`                -> Gas temperature, as ``\\log_{10}(T \\, / \\, \\mathrm{K})``.
      + `:pressure`                   -> Gas pressure.
  - `z_quantity::Symbol`: Quantity for the z axis (weights). The options are:

      + `:stellar_mass`               -> Stellar mass.
      + `:gas_mass`                   -> Gas mass.
      + `:hydrogen_mass`              -> Hydrogen mass.
      + `:dm_mass`                    -> Dark matter mass.
      + `:bh_mass`                    -> Black hole mass.
      + `:molecular_mass`             -> Molecular hydrogen (``\\mathrm{H_2}``) mass.
      + `:atomic_mass`                -> Atomic hydrogen (``\\mathrm{HI}``) mass.
      + `:ionized_mass`               -> Ionized hydrogen (``\\mathrm{HII}``) mass.
      + `:neutral_mass`               -> Neutral hydrogen (``\\mathrm{HI + H_2}``) mass.
      + `:molecular_fraction`         -> Gas mass fraction of molecular hydrogen.
      + `:atomic_fraction`            -> Gas mass fraction of atomic hydrogen.
      + `:ionized_fraction`           -> Gas mass fraction of ionized hydrogen.
      + `:neutral_fraction`           -> Gas mass fraction of neutral hydrogen.
      + `:molecular_neutral_fraction` -> Fraction of molecular hydrogen in the neutral gas.
      + `:gas_mass_density`           -> Gas mass density.
      + `:hydrogen_mass_density`      -> Hydrogen mass density.
      + `:gas_number_density`         -> Gas number density.
      + `:molecular_number_density`   -> Molecular hydrogen number density.
      + `:atomic_number_density`      -> Atomic hydrogen number density.
      + `:ionized_number_density`     -> Ionized hydrogen number density.
      + `:neutral_number_density`     -> Neutral hydrogen number density.
      + `:gas_metallicity`            -> Mass fraction of all elements above He in the gas (solar units).
      + `:stellar_metallicity`        -> Mass fraction of all elements above He in the stars (solar units).
      + `:X_gas_abundance`            -> Gas abundance of element ``\\mathrm{X}``, as ``12 + \\log_{10}(\\mathrm{X \\, / \\, H})``. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
      + `:X_stellar_abundance`        -> Stellar abundance of element ``\\mathrm{X}``, as ``12 + \\log_{10}(\\mathrm{X \\, / \\, H})``. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
      + `:stellar_radial_distance`    -> Distance of every stellar particle to the origin.
      + `:gas_radial_distance`        -> Distance of every gas cell to the origin.
      + `:dm_radial_distance`         -> Distance of every dark matter particle to the origin.
      + `:stellar_xy_distance`        -> Projected distance of every stellar particle to the origin.
      + `:gas_xy_distance`            -> Projected distance of every gas cell to the origin.
      + `:dm_xy_distance`             -> Projected distance of every dark matter particle to the origin.
      + `:stellar_circularity`        -> Stellar circularity.
      + `:stellar_vcirc`              -> Stellar circular velocity.
      + `:stellar_vradial`            -> Stellar radial speed.
      + `:stellar_vtangential`        -> Stellar tangential speed.
      + `:stellar_vzstar`             -> Stellar speed in the z direction, computed as ``v_z \\, \\mathrm{sign}(z)``.
      + `:stellar_age`                -> Stellar age.
      + `:sfr`                        -> The star formation rate.
      + `:ssfr`                       -> The specific star formation rate.
      + `:observational_sfr`          -> The star formation rate of the last `AGE_RESOLUTION`.
      + `:observational_ssfr`         -> The specific star formation rate of the last `AGE_RESOLUTION`.
      + `:temperature`                -> Gas temperature, as ``\\log_{10}(T \\, / \\, \\mathrm{K})``.
      + `:pressure`                   -> Gas pressure.
  - `z_unit::Unitful.Units`: Target unit for the z axis.
  - `x_range::Union{NTuple{2,<:Number},Nothing}=nothing`: x axis range for the histogram grid. If set to `nothing`, the extrema of the values will be used.
  - `y_range::Union{NTuple{2,<:Number},Nothing}=nothing`: y axis range for the histogram grid. If set to `nothing`, the extrema of the values will be used.
  - `x_log::Union{Unitful.Units,Nothing}=nothing`: Desired unit of `x_quantity`, if you want to use log10(`x_quantity`) for the x axis.
  - `y_log::Union{Unitful.Units,Nothing}=nothing`: Desired unit of `y_quantity`, if you want to use log10(`y_quantity`) for the y axis.
  - `total::Bool=true`: If the sum (default) or the mean of `z_quantity` will be used as the value of each pixel.
  - `n_bins::Int=100`: Number of bins per side of the grid.
  - `filter_function::Function=filterNothing`: A function with the signature:

    `filter_function(data_dict) -> indices`

    where

      + `data_dict::Dict`: A dictionary with the following shape:

        * `:sim_data`          -> ::Simulation (see [`Simulation`](@ref)).
        * `:snap_data`         -> ::Snapshot (see [`Snapshot`](@ref)).
        * `:gc_data`           -> ::GroupCatalog (see [`GroupCatalog`](@ref)).
        * `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * ...
        * `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * ...
      + `indices::Dict`: A dictionary with the following shape:

        * `cell/particle type` -> idxs::IndexType
        * `cell/particle type` -> idxs::IndexType
        * `cell/particle type` -> idxs::IndexType
        * ...

# Returns

  - A tuple with three elements:

      + A vector with the x coordinates of the grid.
      + A vector with the y coordinates of the grid.
      + A matrix with the counts.
"""
function daScatterWeightedDensity(
    data_dict::Dict,
    x_quantity::Symbol,
    y_quantity::Symbol,
    z_quantity::Symbol,
    z_unit::Unitful.Units;
    x_range::Union{NTuple{2,<:Number},Nothing}=nothing,
    y_range::Union{NTuple{2,<:Number},Nothing}=nothing,
    x_log::Union{Unitful.Units,Nothing}=nothing,
    y_log::Union{Unitful.Units,Nothing}=nothing,
    total::Bool=true,
    n_bins::Int=100,
    filter_function::Function=filterNothing,
)::Tuple{Vector{<:Number},Vector{<:Number},Matrix{Float64}}

    filtered_dd = filterData(data_dict; filter_function)

    # Compute the values of the quantities
    x_values = scatterQty(filtered_dd, x_quantity)
    y_values = scatterQty(filtered_dd, y_quantity)
    z_values = scatterQty(filtered_dd, z_quantity)

    # If any of the necessary quantities are missing return an empty histogram
    if any(isempty, [x_values, y_values, z_values])
        return 1:n_bins, 1:n_bins, fill(NaN, (n_bins, n_bins))
    end

    (
        length(x_values) == length(y_values) == length(z_values) ||
        throw(ArgumentError("daScatterWeightedDensity: :$(x_quantity), :$(y_quantity), \
        and :$(z_quantity) have a diferent number of values. They should be the same"))
    )

    # If requested, apply log10 to the x axis data, ignoring 0 values
    if !isnothing(x_log)
        null_x_idxs = findall(iszero, x_values)
        x_values    = log10.(deleteat!(ustrip.(x_log, x_values), null_x_idxs))
        y_values    = deleteat!(y_values, null_x_idxs)
        z_values    = deleteat!(z_values, null_x_idxs)
    end

    # If requested, apply log10 to the y axis data, ignoring 0 values
    if !isnothing(y_log)
        null_y_idxs = findall(iszero, y_values)
        x_values    = deleteat!(x_values, null_y_idxs)
        y_values    = log10.(deleteat!(ustrip.(y_log, y_values), null_y_idxs))
        z_values    = deleteat!(z_values, null_y_idxs)
    end

    # If there is no range specified, use the extrema of the x values
    if isnothing(x_range)
        x_range = extrema(x_values)
    end

    # If there is no range specified, use the extrema of the y values
    if isnothing(y_range)
        y_range = extrema(y_values)
    end

    # Compute the bin half width for each axis
    x_bin_h_width = 0.5 * (x_range[2] - x_range[1]) / n_bins
    y_bin_h_width = 0.5 * (y_range[2] - y_range[1]) / n_bins

    # Compute the center value of each bin for each axis
    x_axis = collect(range(x_range[1] + x_bin_h_width; length=n_bins, step=2 * x_bin_h_width))
    y_axis = collect(range(y_range[1] + y_bin_h_width; length=n_bins, step=2 * y_bin_h_width))

    # Compute the 2D histogram
    values = ustrip.(z_unit, histogram2D(
        permutedims(hcat(x_values, y_values), (2, 1)),
        z_values,
        collect(range(x_range[1], x_range[2]; length=n_bins + 1)),
        collect(range(y_range[1], y_range[2]; length=n_bins + 1));
        total,
    ))

    # Set bins with a value of 0 to NaN
    replace!(x -> iszero(x) ? NaN : x, values)

    # The transpose and reverse operation are to conform to the way heatmap! expect the matrix to be structured,
    # and log10 is used to enhance the contrast
    z_axis = reverse!(transpose(log10.(values)), dims=2)

    return x_axis, y_axis, z_axis

end

"""
    daVelocityField(
        data_dict::Dict,
        grid::SquareGrid,
        component::Symbol;
        <keyword arguments>
    )::Tuple{Vector{<:Unitful.Length},Vector{<:Unitful.Length},Matrix{<:Number},Matrix{<:Number}}

Compute a 2D mean velocity field.

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
  - `component::Symbol`: For which cell/particle type the velocity field will be computed. The possibilities are the keys of [`PARTICLE_INDEX`](@ref).
  - `projection_plane::Symbol=:xy`: To which plane the cells/particles will be projected. The options are `:xy`, `:xz`, and `:yz`.
  - `velocity_units::Bool=false`: If the velocity will be given as an `Unitful.Quantity` with units or as a `Flot64` (in which case the underlying unit is ``\\mathrm{km} \\, \\mathrm{s}^{-1}``).
  - `filter_function::Function=filterNothing`: A function with the signature:

    `filter_function(data_dict) -> indices`

    where

      + `data_dict::Dict`: A dictionary with the following shape:

        * `:sim_data`          -> ::Simulation (see [`Simulation`](@ref)).
        * `:snap_data`         -> ::Snapshot (see [`Snapshot`](@ref)).
        * `:gc_data`           -> ::GroupCatalog (see [`GroupCatalog`](@ref)).
        * `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * ...
        * `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * ...
      + `indices::Dict`: A dictionary with the following shape:

        * `cell/particle type` -> idxs::IndexType
        * `cell/particle type` -> idxs::IndexType
        * `cell/particle type` -> idxs::IndexType
        * ...

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
    velocity_units::Bool=false,
    filter_function::Function=filterNothing,
)::Tuple{Vector{<:Unitful.Length},Vector{<:Unitful.Length},Matrix{<:Number},Matrix{<:Number}}

    filtered_dd = filterData(data_dict; filter_function)

    positions  = filtered_dd[component]["POS "]
    velocities = filtered_dd[component]["VEL "]

    # If any of the necessary quantities are missing return an empty velocity field
    if any(isempty, [positions, velocities])
        return grid.x_ticks, grid.y_ticks, zeros(size(grid.grid)), zeros(size(grid.grid))
    end

    # Project the cell/particles to the chosen plane
    if projection_plane == :xy

        pos_2D = positions[[1, 2], :]

        # Compute the components of the mean velocity
        vx = histogram2D(pos_2D, vec(velocities[1, :]), grid; total=false)
        vy = histogram2D(pos_2D, vec(velocities[2, :]), grid; total=false)

    elseif projection_plane == :xz

        pos_2D = positions[[1, 3], :]

        # Compute the components of the mean velocity
        vx = histogram2D(pos_2D, vec(velocities[1, :]), grid; total=false)
        vy = histogram2D(pos_2D, vec(velocities[3, :]), grid; total=false)

    elseif projection_plane == :yz

        pos_2D = positions[[2, 3], :]

        # Compute the components of the mean velocity
        vx = histogram2D(pos_2D, vec(velocities[2, :]), grid; total=false)
        vy = histogram2D(pos_2D, vec(velocities[3, :]), grid; total=false)

    else

        throw(ArgumentError("daVelocityField: The argument `projection_plane` must be \
        :xy, :xz or :yz, but I got :$(projection_plane)"))

    end

    # The transpose and reverse operation are to conform to the way arrows! expect the matrix to be structured
    vx = collect(reverse!(transpose(vx), dims=2))
    vy = collect(reverse!(transpose(vy), dims=2))

    if !velocity_units
        vx = ustrip.(u"km*s^-1", vx)
        vy = ustrip.(u"km*s^-1", vy)
    end

    return grid.x_ticks, grid.y_ticks, vx, vy

end

"""
    daIntegrateGalaxy(
        data_dict::Dict,
        x_quantity::Symbol,
        y_quantity::Symbol;
        <keyword arguments>
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
  - `x_quantity::Symbol`: Quantity for the x axis. The options are:

      + `:stellar_mass`           -> Stellar mass.
      + `:gas_mass`               -> Gas mass.
      + `:hydrogen_mass`          -> Hydrogen mass.
      + `:dm_mass`                -> Dark matter mass.
      + `:bh_mass`                -> Black hole mass.
      + `:molecular_mass`         -> Molecular hydrogen (``\\mathrm{H_2}``) mass.
      + `:atomic_mass`            -> Atomic hydrogen (``\\mathrm{HI}``) mass.
      + `:ionized_mass`           -> Ionized hydrogen (``\\mathrm{HII}``) mass.
      + `:neutral_mass`           -> Neutral hydrogen (``\\mathrm{HI + H_2}``) mass.
      + `:stellar_number`         -> Number of stellar particles.
      + `:gas_number`             -> Number of gas cells.
      + `:dm_number`              -> Number of dark matter particles.
      + `:bh_number`              -> Number of black hole particles.
      + `:molecular_fraction`     -> Gas mass fraction of molecular hydrogen.
      + `:atomic_fraction`        -> Gas mass fraction of atomic hydrogen.
      + `:ionized_fraction`       -> Gas mass fraction of ionized hydrogen.
      + `:neutral_fraction`       -> Gas mass fraction of neutral hydrogen.
      + `:stellar_area_density`   -> Stellar area mass density, for a radius of `DISK_R`.
      + `:gas_area_density`       -> Gas area mass density, for a radius of `DISK_R`.
      + `:molecular_area_density` -> Molecular hydrogen area mass density, for a radius of `DISK_R`.
      + `:atomic_area_density`    -> Atomic hydrogen area mass density, for a radius of `DISK_R`.
      + `:ionized_area_density`   -> Ionized hydrogen area mass density, for a radius of `DISK_R`.
      + `:neutral_area_density`   -> Neutral hydrogen area mass density, for a radius of `DISK_R`.
      + `:sfr_area_density`       -> Star formation rate area density, for the last `AGE_RESOLUTION` and a radius of `DISK_R`.
      + `:gas_metallicity`        -> Mass fraction of all elements above He in the gas (solar units).
      + `:stellar_metallicity`    -> Mass fraction of all elements above He in the stars (solar units).
      + `:X_gas_abundance`        -> Gas abundance of element ``\\mathrm{X}``, as ``12 + \\log_{10}(\\mathrm{X \\, / \\, H})``. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
      + `:X_stellar_abundance`    -> Stellar abundance of element ``\\mathrm{X}``, as ``12 + \\log_{10}(\\mathrm{X \\, / \\, H})``. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
      + `:stellar_specific_am`    -> Norm of the stellar specific angular momentum.
      + `:gas_specific_am`        -> Norm of the gas specific angular momentum.
      + `:dm_specific_am`         -> Norm of the dark matter specific angular momentum.
      + `:sfr`                    -> The star formation rate.
      + `:ssfr`                   -> The specific star formation rate.
      + `:observational_sfr`      -> The star formation rate of the last `AGE_RESOLUTION`.
      + `:observational_ssfr`     -> The specific star formation rate of the last `AGE_RESOLUTION`.
      + `:scale_factor`           -> Scale factor.
      + `:redshift`               -> Redshift.
      + `:physical_time`          -> Physical time since the Big Bang.
      + `:lookback_time`          -> Physical time left to reach the last snapshot.
  - `y_quantity::Symbol`: Quantity for the y axis. The options are:

      + `:stellar_mass`           -> Stellar mass.
      + `:gas_mass`               -> Gas mass.
      + `:hydrogen_mass`          -> Hydrogen mass.
      + `:dm_mass`                -> Dark matter mass.
      + `:bh_mass`                -> Black hole mass.
      + `:molecular_mass`         -> Molecular hydrogen (``\\mathrm{H_2}``) mass.
      + `:atomic_mass`            -> Atomic hydrogen (``\\mathrm{HI}``) mass.
      + `:ionized_mass`           -> Ionized hydrogen (``\\mathrm{HII}``) mass.
      + `:neutral_mass`           -> Neutral hydrogen (``\\mathrm{HI + H_2}``) mass.
      + `:stellar_number`         -> Number of stellar particles.
      + `:gas_number`             -> Number of gas cells.
      + `:dm_number`              -> Number of dark matter particles.
      + `:bh_number`              -> Number of black hole particles.
      + `:molecular_fraction`     -> Gas mass fraction of molecular hydrogen.
      + `:atomic_fraction`        -> Gas mass fraction of atomic hydrogen.
      + `:ionized_fraction`       -> Gas mass fraction of ionized hydrogen.
      + `:neutral_fraction`       -> Gas mass fraction of neutral hydrogen.
      + `:stellar_area_density`   -> Stellar area mass density, for a radius of `DISK_R`.
      + `:gas_area_density`       -> Gas area mass density, for a radius of `DISK_R`.
      + `:molecular_area_density` -> Molecular hydrogen area mass density, for a radius of `DISK_R`.
      + `:atomic_area_density`    -> Atomic hydrogen area mass density, for a radius of `DISK_R`.
      + `:ionized_area_density`   -> Ionized hydrogen area mass density, for a radius of `DISK_R`.
      + `:neutral_area_density`   -> Neutral hydrogen area mass density, for a radius of `DISK_R`.
      + `:sfr_area_density`       -> Star formation rate area density, for the last `AGE_RESOLUTION` and a radius of `DISK_R`.
      + `:gas_metallicity`        -> Mass fraction of all elements above He in the gas (solar units).
      + `:stellar_metallicity`    -> Mass fraction of all elements above He in the stars (solar units).
      + `:X_gas_abundance`        -> Gas abundance of element ``\\mathrm{X}``, as ``12 + \\log_{10}(\\mathrm{X \\, / \\, H})``. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
      + `:X_stellar_abundance`    -> Stellar abundance of element ``\\mathrm{X}``, as ``12 + \\log_{10}(\\mathrm{X \\, / \\, H})``. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
      + `:stellar_specific_am`    -> Norm of the stellar specific angular momentum.
      + `:gas_specific_am`        -> Norm of the gas specific angular momentum.
      + `:dm_specific_am`         -> Norm of the dark matter specific angular momentum.
      + `:sfr`                    -> The star formation rate.
      + `:ssfr`                   -> The specific star formation rate.
      + `:observational_sfr`      -> The star formation rate of the last `AGE_RESOLUTION`.
      + `:observational_ssfr`     -> The specific star formation rate of the last `AGE_RESOLUTION`.
      + `:scale_factor`           -> Scale factor.
      + `:redshift`               -> Redshift.
      + `:physical_time`          -> Physical time since the Big Bang.
      + `:lookback_time`          -> Physical time left to reach the last snapshot.
  - `filter_function::Function=filterNothing`: A function with the signature:

    `filter_function(data_dict) -> indices`

    where

      + `data_dict::Dict`: A dictionary with the following shape:

        * `:sim_data`          -> ::Simulation (see [`Simulation`](@ref)).
        * `:snap_data`         -> ::Snapshot (see [`Snapshot`](@ref)).
        * `:gc_data`           -> ::GroupCatalog (see [`GroupCatalog`](@ref)).
        * `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * ...
        * `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * ...
      + `indices::Dict`: A dictionary with the following shape:

        * `cell/particle type` -> idxs::IndexType
        * `cell/particle type` -> idxs::IndexType
        * `cell/particle type` -> idxs::IndexType
        * ...

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
    daScatterGalaxy(
        data_dict::Dict,
        x_quantity::Symbol,
        y_quantity::Symbol;
        <keyword arguments>
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
  - `x_quantity::Symbol`: Quantity for the x axis. The options are:

      + `:stellar_mass`               -> Stellar mass.
      + `:gas_mass`                   -> Gas mass.
      + `:hydrogen_mass`              -> Hydrogen mass.
      + `:dm_mass`                    -> Dark matter mass.
      + `:bh_mass`                    -> Black hole mass.
      + `:molecular_mass`             -> Molecular hydrogen (``\\mathrm{H_2}``) mass.
      + `:atomic_mass`                -> Atomic hydrogen (``\\mathrm{HI}``) mass.
      + `:ionized_mass`               -> Ionized hydrogen (``\\mathrm{HII}``) mass.
      + `:neutral_mass`               -> Neutral hydrogen (``\\mathrm{HI + H_2}``) mass.
      + `:molecular_fraction`         -> Gas mass fraction of molecular hydrogen.
      + `:atomic_fraction`            -> Gas mass fraction of atomic hydrogen.
      + `:ionized_fraction`           -> Gas mass fraction of ionized hydrogen.
      + `:neutral_fraction`           -> Gas mass fraction of neutral hydrogen.
      + `:molecular_neutral_fraction` -> Fraction of molecular hydrogen in the neutral gas.
      + `:gas_mass_density`           -> Gas mass density.
      + `:hydrogen_mass_density`      -> Hydrogen mass density.
      + `:gas_number_density`         -> Gas number density.
      + `:molecular_number_density`   -> Molecular hydrogen number density.
      + `:atomic_number_density`      -> Atomic hydrogen number density.
      + `:ionized_number_density`     -> Ionized hydrogen number density.
      + `:neutral_number_density`     -> Neutral hydrogen number density.
      + `:gas_metallicity`            -> Mass fraction of all elements above He in the gas (solar units).
      + `:stellar_metallicity`        -> Mass fraction of all elements above He in the stars (solar units).
      + `:X_gas_abundance`            -> Gas abundance of element ``\\mathrm{X}``, as ``12 + \\log_{10}(\\mathrm{X \\, / \\, H})``. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
      + `:X_stellar_abundance`        -> Stellar abundance of element ``\\mathrm{X}``, as ``12 + \\log_{10}(\\mathrm{X \\, / \\, H})``. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
      + `:stellar_radial_distance`    -> Distance of every stellar particle to the origin.
      + `:gas_radial_distance`        -> Distance of every gas cell to the origin.
      + `:dm_radial_distance`         -> Distance of every dark matter particle to the origin.
      + `:stellar_xy_distance`        -> Projected distance of every stellar particle to the origin.
      + `:gas_xy_distance`            -> Projected distance of every gas cell to the origin.
      + `:dm_xy_distance`             -> Projected distance of every dark matter particle to the origin.
      + `:stellar_circularity`        -> Stellar circularity.
      + `:stellar_vcirc`              -> Stellar circular velocity.
      + `:stellar_vradial`            -> Stellar radial speed.
      + `:stellar_vtangential`        -> Stellar tangential speed.
      + `:stellar_vzstar`             -> Stellar speed in the z direction, computed as ``v_z \\, \\mathrm{sign}(z)``.
      + `:stellar_age`                -> Stellar age.
      + `:sfr`                        -> The star formation rate.
      + `:ssfr`                       -> The specific star formation rate.
      + `:observational_sfr`          -> The star formation rate of the last `AGE_RESOLUTION`.
      + `:observational_ssfr`         -> The specific star formation rate of the last `AGE_RESOLUTION`.
      + `:temperature`                -> Gas temperature, as ``\\log_{10}(T \\, / \\, \\mathrm{K})``.
      + `:pressure`                   -> Gas pressure.
  - `y_quantity::Symbol`: Quantity for the y axis. The options are:

      + `:stellar_mass`               -> Stellar mass.
      + `:gas_mass`                   -> Gas mass.
      + `:hydrogen_mass`              -> Hydrogen mass.
      + `:dm_mass`                    -> Dark matter mass.
      + `:bh_mass`                    -> Black hole mass.
      + `:molecular_mass`             -> Molecular hydrogen (``\\mathrm{H_2}``) mass.
      + `:atomic_mass`                -> Atomic hydrogen (``\\mathrm{HI}``) mass.
      + `:ionized_mass`               -> Ionized hydrogen (``\\mathrm{HII}``) mass.
      + `:neutral_mass`               -> Neutral hydrogen (``\\mathrm{HI + H_2}``) mass.
      + `:molecular_fraction`         -> Gas mass fraction of molecular hydrogen.
      + `:atomic_fraction`            -> Gas mass fraction of atomic hydrogen.
      + `:ionized_fraction`           -> Gas mass fraction of ionized hydrogen.
      + `:neutral_fraction`           -> Gas mass fraction of neutral hydrogen.
      + `:molecular_neutral_fraction` -> Fraction of molecular hydrogen in the neutral gas.
      + `:gas_mass_density`           -> Gas mass density.
      + `:hydrogen_mass_density`      -> Hydrogen mass density.
      + `:gas_number_density`         -> Gas number density.
      + `:molecular_number_density`   -> Molecular hydrogen number density.
      + `:atomic_number_density`      -> Atomic hydrogen number density.
      + `:ionized_number_density`     -> Ionized hydrogen number density.
      + `:neutral_number_density`     -> Neutral hydrogen number density.
      + `:gas_metallicity`            -> Mass fraction of all elements above He in the gas (solar units).
      + `:stellar_metallicity`        -> Mass fraction of all elements above He in the stars (solar units).
      + `:X_gas_abundance`            -> Gas abundance of element ``\\mathrm{X}``, as ``12 + \\log_{10}(\\mathrm{X \\, / \\, H})``. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
      + `:X_stellar_abundance`        -> Stellar abundance of element ``\\mathrm{X}``, as ``12 + \\log_{10}(\\mathrm{X \\, / \\, H})``. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
      + `:stellar_radial_distance`    -> Distance of every stellar particle to the origin.
      + `:gas_radial_distance`        -> Distance of every gas cell to the origin.
      + `:dm_radial_distance`         -> Distance of every dark matter particle to the origin.
      + `:stellar_xy_distance`        -> Projected distance of every stellar particle to the origin.
      + `:gas_xy_distance`            -> Projected distance of every gas cell to the origin.
      + `:dm_xy_distance`             -> Projected distance of every dark matter particle to the origin.
      + `:stellar_circularity`        -> Stellar circularity.
      + `:stellar_vcirc`              -> Stellar circular velocity.
      + `:stellar_vradial`            -> Stellar radial speed.
      + `:stellar_vtangential`        -> Stellar tangential speed.
      + `:stellar_vzstar`             -> Stellar speed in the z direction, computed as ``v_z \\, \\mathrm{sign}(z)``.
      + `:stellar_age`                -> Stellar age.
      + `:sfr`                        -> The star formation rate.
      + `:ssfr`                       -> The specific star formation rate.
      + `:observational_sfr`          -> The star formation rate of the last `AGE_RESOLUTION`.
      + `:observational_ssfr`         -> The specific star formation rate of the last `AGE_RESOLUTION`.
      + `:temperature`                -> Gas temperature, as ``\\log_{10}(T \\, / \\, \\mathrm{K})``.
      + `:pressure`                   -> Gas pressure.
  - `filter_function::Function=filterNothing`: A function with the signature:

    `filter_function(data_dict) -> indices`

    where

      + `data_dict::Dict`: A dictionary with the following shape:

        * `:sim_data`          -> ::Simulation (see [`Simulation`](@ref)).
        * `:snap_data`         -> ::Snapshot (see [`Snapshot`](@ref)).
        * `:gc_data`           -> ::GroupCatalog (see [`GroupCatalog`](@ref)).
        * `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * ...
        * `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * ...
      + `indices::Dict`: A dictionary with the following shape:

        * `cell/particle type` -> idxs::IndexType
        * `cell/particle type` -> idxs::IndexType
        * `cell/particle type` -> idxs::IndexType
        * ...

# Returns

  - A tuple with two elements:

      + A vector with the values of `x_quantity`.
      + A vector with the values of `y_quantity`.
"""
function daScatterGalaxy(
    data_dict::Dict,
    x_quantity::Symbol,
    y_quantity::Symbol;
    filter_function::Function=filterNothing,
)::NTuple{2,Vector{<:Number}}

    filtered_dd = filterData(data_dict; filter_function)

    x_axis = scatterQty(filtered_dd, x_quantity)
    y_axis = scatterQty(filtered_dd, y_quantity)

    idx = sortperm(x_axis)

    (
        length(x_axis) == length(y_axis) ||
        throw(ArgumentError("daScatterGalaxy: :$(x_quantity) and :$(y_quantity) have a diferent \
        number of values. They should be the same"))
    )

    return x_axis[idx], y_axis[idx]

end

"""
    daGasFractions(
        data_dict::Dict,
        quantity::Symbol,
        edges::Vector{<:Number};
        <keyword arguments>
    )::Union{NTuple{2,Vector{<:Number}},Nothing}

Compute the values for a gas fraction bar plot.

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
  - `quantity::Symbol`: Target quantity. The possibilities are:

      + `:gas_mass`                   -> Gas mass.
      + `:hydrogen_mass`              -> Hydrogen mass.
      + `:molecular_mass`             -> Molecular hydrogen (``\\mathrm{H_2}``) mass.
      + `:atomic_mass`                -> Atomic hydrogen (``\\mathrm{HI}``) mass.
      + `:ionized_mass`               -> Ionized hydrogen (``\\mathrm{HII}``) mass.
      + `:neutral_mass`               -> Neutral hydrogen (``\\mathrm{HI + H_2}``) mass.
      + `:molecular_fraction`         -> Gas mass fraction of molecular hydrogen.
      + `:atomic_fraction`            -> Gas mass fraction of atomic hydrogen.
      + `:ionized_fraction`           -> Gas mass fraction of ionized hydrogen.
      + `:neutral_fraction`           -> Gas mass fraction of neutral hydrogen.
      + `:molecular_neutral_fraction` -> Fraction of molecular hydrogen in the neutral gas.
      + `:gas_mass_density`           -> Gas mass density.
      + `:hydrogen_mass_density`      -> Hydrogen mass density.
      + `:gas_number_density`         -> Gas number density.
      + `:molecular_number_density`   -> Molecular hydrogen number density.
      + `:atomic_number_density`      -> Atomic hydrogen number density.
      + `:ionized_number_density`     -> Ionized hydrogen number density.
      + `:neutral_number_density`     -> Neutral hydrogen number density.
      + `:gas_metallicity`            -> Mass fraction of all elements above He in the gas (solar units).
      + `:X_gas_abundance`            -> Gas abundance of element ``\\mathrm{X}``, as ``12 + \\log_{10}(\\mathrm{X \\, / \\, H})``. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
      + `:gas_radial_distance`        -> Distance of every gas cell to the origin.
      + `:gas_xy_distance`            -> Projected distance of every gas cell to the origin.
      + `:temperature`                -> Gas temperature, as ``\\log_{10}(T \\, / \\, \\mathrm{K})``.
      + `:pressure`                   -> Gas pressure.
  - `edges::Vector{<:Number}`: A sorted list of bin edges.
  - `include_stars::Bool=false`: If the stars will be included as one of the gas phases. It will only work for simulations with our routine.
  - `filter_function::Function=filterNothing`: A function with the signature:

    `filter_function(data_dict) -> indices`

    where

      + `data_dict::Dict`: A dictionary with the following shape:

        * `:sim_data`          -> ::Simulation (see [`Simulation`](@ref)).
        * `:snap_data`         -> ::Snapshot (see [`Snapshot`](@ref)).
        * `:gc_data`           -> ::GroupCatalog (see [`GroupCatalog`](@ref)).
        * `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * ...
        * `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * ...
      + `indices::Dict`: A dictionary with the following shape:

        * `cell/particle type` -> idxs::IndexType
        * `cell/particle type` -> idxs::IndexType
        * `cell/particle type` -> idxs::IndexType
        * ...

# Returns

  - A tuple with two elements:

      + A vector with the positions of each bar.
      + A vector with the height of each bar.
"""
function daGasFractions(
    data_dict::Dict,
    quantity::Symbol,
    edges::Vector{<:Number};
    include_stars::Bool=false,
    filter_function::Function=filterNothing,
)::Union{NTuple{2,Vector{<:Number}},Nothing}

    filtered_dd = filterData(data_dict; filter_function)

    # Compute the number of bins for the gas quantity
    n_bins = length(edges) - 1

    # Compute the number of bars per bin
    n_bars = include_stars ? 4 : 3

    # Compute the gas quantities
    gas_qty  = scatterQty(filtered_dd, quantity)

    # If any of the necessary quantities are missing return nothing
    !isempty(gas_qty) || return nothing

    # Compute the mass of each gas phase
    if include_stars

        ionized_mass   = computeIonizedMass(filtered_dd; normalize=false)
        atomic_mass    = computeAtomicMass(filtered_dd; normalize=false)
        molecular_mass = computeMolecularMass(filtered_dd; normalize=false)
        stellar_mass   = computeStellarGasMass(filtered_dd)

    else

        ionized_mass   = computeIonizedMass(filtered_dd)
        atomic_mass    = computeAtomicMass(filtered_dd)
        molecular_mass = computeMolecularMass(filtered_dd)

    end

    # If any of the necessary quantities are missing return nothing
    !any(isempty, [ionized_mass, atomic_mass, molecular_mass]) || return nothing

    # Allocate memory
    percents = Vector{Float64}(undef, n_bins * n_bars)

    # Loop over each bin
    for i in 1:n_bins

        # Compute the indices of the gas cells inside the current bin
        qty_idx = findall(q->edges[i] < q <= edges[i + 1], gas_qty)

        if isempty(qty_idx)

            percents[((i - 1) * n_bars + 1):(i * n_bars)] .= 0.0

        else

            if include_stars
                total_masses = [
                    sum(molecular_mass[qty_idx]),
                    sum(atomic_mass[qty_idx]),
                    sum(ionized_mass[qty_idx]),
                    sum(stellar_mass[qty_idx]),
                ]
            else
                total_masses = [
                    sum(molecular_mass[qty_idx]),
                    sum(atomic_mass[qty_idx]),
                    sum(ionized_mass[qty_idx]),
                ]
            end

            # Compute the mass fraction of each phase inside the current bin
            fractions = uconvert.(Unitful.NoUnits, (total_masses ./ sum(total_masses)) .* 100)

            percents[((i - 1) * n_bars + 1):(i * n_bars)] .= fractions

        end

    end

    return repeat(1:n_bins, inner=n_bars), percents

end

"""
    daStellarMetallictyHistogram(data_dict::Dict)::Union{Tuple{Vector{Float64}},Nothing}

Compute the stellar metallicity, for an histogram.

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

# Returns

  - A Tuple with one elements:

      + A Vector with the stellar metallicites.
"""
function daStellarMetallictyHistogram(data_dict::Dict)::Union{Tuple{Vector{Float64}},Nothing}

    metallicity = scatterQty(data_dict, :stellar_metallicity)

    !isempty(metallicity) || return nothing

    return (metallicity,)

end

"""
    daStellarBTHistogram(data_dict::Dict)::Union{Tuple{Vector{<:Unitful.Time}},Nothing}

Compute the stellar birth times, for an histogram.

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

# Returns

  - A Tuple with one elements:

      + A Vector with the birth times.
"""
function daStellarBTHistogram(data_dict::Dict)::Union{Tuple{Vector{<:Unitful.Time}},Nothing}

    birth_ticks = data_dict[:stars]["GAGE"]

    !isempty(birth_ticks) || return nothing

    if data_dict[:sim_data].cosmological
        # Go from scale factor to physical time
        birth_times = computeTime(birth_ticks, data_dict[:snap_data].header)
    else
        birth_times = birth_ticks
    end

    return (birth_times,)

end

####################################################################################################
# Signature for the plotTimeSeries function in ./src/pipelines.jl.
####################################################################################################
#
# A data analysis functions for plotTimeSeries must take a Simulation struct, and return two
# vectors. It should return `nothing` if the input data has some problem that prevents computation
# (e.g. is empty).
#
# Expected signature:
#
#   da_function(sim_data, args...; kw_args...) -> (processed_data_x, processed_data_y)
#
# where:
#
#   - sim_data::Simulation, see the Simulation struct in ./src/constants.jl
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
  - `x_quantity::Symbol`: Quantity for the x axis. The options are:

      + `:stellar_mass`           -> Stellar mass.
      + `:gas_mass`               -> Gas mass.
      + `:hydrogen_mass`          -> Hydrogen mass.
      + `:dm_mass`                -> Dark matter mass.
      + `:bh_mass`                -> Black hole mass.
      + `:molecular_mass`         -> Molecular hydrogen (``\\mathrm{H_2}``) mass.
      + `:atomic_mass`            -> Atomic hydrogen (``\\mathrm{HI}``) mass.
      + `:ionized_mass`           -> Ionized hydrogen (``\\mathrm{HII}``) mass.
      + `:neutral_mass`           -> Neutral hydrogen (``\\mathrm{HI + H_2}``) mass.
      + `:stellar_number`         -> Number of stellar particles.
      + `:gas_number`             -> Number of gas cells.
      + `:dm_number`              -> Number of dark matter particles.
      + `:bh_number`              -> Number of black hole particles.
      + `:molecular_fraction`     -> Gas mass fraction of molecular hydrogen.
      + `:atomic_fraction`        -> Gas mass fraction of atomic hydrogen.
      + `:ionized_fraction`       -> Gas mass fraction of ionized hydrogen.
      + `:neutral_fraction`       -> Gas mass fraction of neutral hydrogen.
      + `:stellar_area_density`   -> Stellar area mass density, for a radius of `DISK_R`.
      + `:gas_area_density`       -> Gas area mass density, for a radius of `DISK_R`.
      + `:molecular_area_density` -> Molecular hydrogen area mass density, for a radius of `DISK_R`.
      + `:atomic_area_density`    -> Atomic hydrogen area mass density, for a radius of `DISK_R`.
      + `:ionized_area_density`   -> Ionized hydrogen area mass density, for a radius of `DISK_R`.
      + `:neutral_area_density`   -> Neutral hydrogen area mass density, for a radius of `DISK_R`.
      + `:sfr_area_density`       -> Star formation rate area density, for the last `AGE_RESOLUTION` and a radius of `DISK_R`.
      + `:gas_metallicity`        -> Mass fraction of all elements above He in the gas (solar units).
      + `:stellar_metallicity`    -> Mass fraction of all elements above He in the stars (solar units).
      + `:X_gas_abundance`        -> Gas abundance of element ``\\mathrm{X}``, as ``12 + \\log_{10}(\\mathrm{X \\, / \\, H})``. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
      + `:X_stellar_abundance`    -> Stellar abundance of element ``\\mathrm{X}``, as ``12 + \\log_{10}(\\mathrm{X \\, / \\, H})``. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
      + `:stellar_specific_am`    -> Norm of the stellar specific angular momentum.
      + `:gas_specific_am`        -> Norm of the gas specific angular momentum.
      + `:dm_specific_am`         -> Norm of the dark matter specific angular momentum.
      + `:sfr`                    -> The star formation rate.
      + `:ssfr`                   -> The specific star formation rate.
      + `:observational_sfr`      -> The star formation rate of the last `AGE_RESOLUTION`.
      + `:observational_ssfr`     -> The specific star formation rate of the last `AGE_RESOLUTION`.
      + `:scale_factor`           -> Scale factor.
      + `:redshift`               -> Redshift.
      + `:physical_time`          -> Physical time since the Big Bang.
      + `:lookback_time`          -> Physical time left to reach the last snapshot.
  - `y_quantity::Symbol`: Quantity for the y axis. The options are:

      + `:stellar_mass`           -> Stellar mass.
      + `:gas_mass`               -> Gas mass.
      + `:hydrogen_mass`          -> Hydrogen mass.
      + `:dm_mass`                -> Dark matter mass.
      + `:bh_mass`                -> Black hole mass.
      + `:molecular_mass`         -> Molecular hydrogen (``\\mathrm{H_2}``) mass.
      + `:atomic_mass`            -> Atomic hydrogen (``\\mathrm{HI}``) mass.
      + `:ionized_mass`           -> Ionized hydrogen (``\\mathrm{HII}``) mass.
      + `:neutral_mass`           -> Neutral hydrogen (``\\mathrm{HI + H_2}``) mass.
      + `:stellar_number`         -> Number of stellar particles.
      + `:gas_number`             -> Number of gas cells.
      + `:dm_number`              -> Number of dark matter particles.
      + `:bh_number`              -> Number of black hole particles.
      + `:molecular_fraction`     -> Gas mass fraction of molecular hydrogen.
      + `:atomic_fraction`        -> Gas mass fraction of atomic hydrogen.
      + `:ionized_fraction`       -> Gas mass fraction of ionized hydrogen.
      + `:neutral_fraction`       -> Gas mass fraction of neutral hydrogen.
      + `:stellar_area_density`   -> Stellar area mass density, for a radius of `DISK_R`.
      + `:gas_area_density`       -> Gas area mass density, for a radius of `DISK_R`.
      + `:molecular_area_density` -> Molecular hydrogen area mass density, for a radius of `DISK_R`.
      + `:atomic_area_density`    -> Atomic hydrogen area mass density, for a radius of `DISK_R`.
      + `:ionized_area_density`   -> Ionized hydrogen area mass density, for a radius of `DISK_R`.
      + `:neutral_area_density`   -> Neutral hydrogen area mass density, for a radius of `DISK_R`.
      + `:sfr_area_density`       -> Star formation rate area density, for the last `AGE_RESOLUTION` and a radius of `DISK_R`.
      + `:gas_metallicity`        -> Mass fraction of all elements above He in the gas (solar units).
      + `:stellar_metallicity`    -> Mass fraction of all elements above He in the stars (solar units).
      + `:X_gas_abundance`        -> Gas abundance of element ``\\mathrm{X}``, as ``12 + \\log_{10}(\\mathrm{X \\, / \\, H})``. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
      + `:X_stellar_abundance`    -> Stellar abundance of element ``\\mathrm{X}``, as ``12 + \\log_{10}(\\mathrm{X \\, / \\, H})``. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
      + `:stellar_specific_am`    -> Norm of the stellar specific angular momentum.
      + `:gas_specific_am`        -> Norm of the gas specific angular momentum.
      + `:dm_specific_am`         -> Norm of the dark matter specific angular momentum.
      + `:sfr`                    -> The star formation rate.
      + `:ssfr`                   -> The specific star formation rate.
      + `:observational_sfr`      -> The star formation rate of the last `AGE_RESOLUTION`.
      + `:observational_ssfr`     -> The specific star formation rate of the last `AGE_RESOLUTION`.
      + `:scale_factor`           -> Scale factor.
      + `:redshift`               -> Redshift.
      + `:physical_time`          -> Physical time since the Big Bang.
      + `:lookback_time`          -> Physical time left to reach the last snapshot.
  - `filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all`: Which cells/particles will be plotted, the options are:

      + `:all`             -> Consider every cell/particle within the simulation box.
      + `:halo`            -> Consider only the cells/particles that belong to the main halo.
      + `:subhalo`         -> Consider only the cells/particles that belong to the main subhalo.
      + `:sphere`          -> Consider only the cell/particle inside a sphere with radius `DISK_R` (see `./src/constants.jl`).
      + `:stellar_subhalo` -> Consider only the cells/particles that belong to the main subhalo.
      + `:all_subhalo`     -> Plot every cell/particle centered around the main subhalo.
      + A dictionary with three entries:

          + `:filter_function` -> The filter function.
          + `:translation`     -> Translation for the simulation box. The posibilities are:

              + `:global_cm`                  -> Selects the center of mass of the whole system as the new origin.
              + `:{component}`                -> Sets the center of mass of the given component (e.g. :stars, :gas, :halo, etc) as the new origin. It can be any of the keys of [`PARTICLE_INDEX`](@ref).
              + `(halo_idx, subhalo_rel_idx)` -> Sets the position of the potencial minimum for the `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo) as the new origin.
              + `(halo_idx, 0)`               -> Sets the center of mass of the `halo_idx::Int` halo as the new origin.
              + `subhalo_abs_idx`             -> Sets the center of mass of the `subhalo_abs_idx::Int` as the new origin.
          + `:rotation`        -> Rotation for the simulation box. The posibilities are:

              + `:zero`                       -> No rotation is appplied.
              + `:global_am`                  -> Sets the angular momentum of the whole system as the new z axis.
              + `:stellar_am`                 -> Sets the stellar angular momentum as the new z axis.
              + `:stellar_pa`                 -> Sets the stellar principal axis as the new coordinate system.
              + `:stellar_subhalo_pa`         -> Sets the principal axis of the stars in the main subhalo as the new coordinate system.
              + `(halo_idx, subhalo_rel_idx)` -> Sets the principal axis of the stars in `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo), as the new coordinate system.
              + `(halo_idx, 0)`               -> Sets the principal axis of the stars in the `halo_idx::Int` halo, as the new coordinate system.
              + `subhalo_abs_idx`             -> Sets the principal axis of the stars in the `subhalo_abs_idx::Int` subhalo as the new coordinate system.
  - `extra_filter::Function=filterNothing`: Filter function that will be applied after the one given by `filter_mode`.
  - `smooth::Int=0`: The result of [`integrateQty`](@ref) will be smooth out using `smooth` bins. Set it to 0 if you want no smoothing.
  - `cumulative::Bool=false`: If the `y_quantity` will be accumulated or not.
  - `fraction::Bool=false`: If the `y_quantity` will be represented as a fraction of the last value. If `cumulative` = true, this will apply to the accumulated values.
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
    filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all,
    extra_filter::Function=filterNothing,
    smooth::Int=0,
    cumulative::Bool=false,
    fraction::Bool=false,
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
    x_axis = Vector{Number}(fill(NaN, length(iterator)))
    y_axis = Vector{Number}(fill(NaN, length(iterator)))

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
        # Filter the data again
        filterData!(data_dict; filter_function=extra_filter)

        # Translate the data
        translateData!(data_dict, translation)

        # Rotate the data
        rotateData!(data_dict, rotation)

        # Compute the value for the x axis
        x_axis[slice_index] = integrateQty(data_dict, x_quantity)

        # Compute the value for the y axis
        y_axis[slice_index] = integrateQty(data_dict, y_quantity)

    end

    if cumulative
        cumsum!(y_axis, y_axis)
    end

    if fraction
        y_axis ./= y_axis[end]
    end

    if iszero(smooth)
        return x_axis, y_axis
    else
        return smoothWindow(x_axis, y_axis, smooth; scaling)
    end

end

"""
    daVirialAccretion(
        sim_data::Simulation;
        <keyword arguments>
    )::NTuple{2,Vector{<:Number}}

Compute the evolution of the accreted mass into the virial radius.

# Arguments

  - `sim_data::Simulation`: Information about the simulation in a [`Simulation`](@ref) object.
  - `filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all`: Which cells/particles will be plotted. Only valid if `tracers` = true. The options are:

      + `:all`             -> Consider every cell/particle within the simulation box.
      + `:halo`            -> Consider only the cells/particles that belong to the main halo.
      + `:subhalo`         -> Consider only the cells/particles that belong to the main subhalo.
      + `:sphere`          -> Consider only the cell/particle inside a sphere with radius `DISK_R` (see `./src/constants.jl`).
      + `:stellar_subhalo` -> Consider only the cells/particles that belong to the main subhalo.
      + `:all_subhalo`     -> Plot every cell/particle centered around the main subhalo.
      + A dictionary with three entries:

          + `:filter_function` -> The filter function.
          + `:translation`     -> Translation for the simulation box. The posibilities are:

              + `:global_cm`                  -> Selects the center of mass of the whole system as the new origin.
              + `:{component}`                -> Sets the center of mass of the given component (e.g. :stars, :gas, :halo, etc) as the new origin. It can be any of the keys of [`PARTICLE_INDEX`](@ref).
              + `(halo_idx, subhalo_rel_idx)` -> Sets the position of the potencial minimum for the `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo) as the new origin.
              + `(halo_idx, 0)`               -> Sets the center of mass of the `halo_idx::Int` halo as the new origin.
              + `subhalo_abs_idx`             -> Sets the center of mass of the `subhalo_abs_idx::Int` as the new origin.
          + `:rotation`        -> Rotation for the simulation box. The posibilities are:

              + `:zero`                       -> No rotation is appplied.
              + `:global_am`                  -> Sets the angular momentum of the whole system as the new z axis.
              + `:stellar_am`                 -> Sets the stellar angular momentum as the new z axis.
              + `:stellar_pa`                 -> Sets the stellar principal axis as the new coordinate system.
              + `:stellar_subhalo_pa`         -> Sets the principal axis of the stars in the main subhalo as the new coordinate system.
              + `(halo_idx, subhalo_rel_idx)` -> Sets the principal axis of the stars in `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo), as the new coordinate system.
              + `(halo_idx, 0)`               -> Sets the principal axis of the stars in the `halo_idx::Int` halo, as the new coordinate system.
              + `subhalo_abs_idx`             -> Sets the principal axis of the stars in the `subhalo_abs_idx::Int` subhalo as the new coordinate system.
  - `halo_idx::Int=1`: Index of the target halo (FoF group). Starts at 1.
  - `tracers::Bool=false`: If tracers will be use to compute the mass accretion. If false, `filter_mode` will be ignored.
  - `smooth::Int=0`: The time series will be smooth out using `smooth` bins. Set it to 0 if you want no smoothing.
  - `warnings::Bool=true`: If a warning will be given when there is missing data.

# Returns

  - A Tuple with two elements:

      + A Vector with the physical times.
      + A Vector with the accreted mass at each time.
"""
function daVirialAccretion(
    sim_data::Simulation;
    filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all,
    halo_idx::Int=1,
    tracers::Bool=false,
    smooth::Int=0,
    warnings::Bool=true,
)::NTuple{2,Vector{<:Number}}

    filter_function, translation, _, request = selectFilter(
        filter_mode,
        Dict(
            :gas         => ["ID  ", "MASS"],
            :stars       => ["ID  ", "MASS"],
            :black_hole  => ["ID  ", "MASS"],
            :group       => ["G_R_Crit200", "G_M_Crit200"],
            :tracer      => ["PAID", "TRID"],
        ),
    )

    # Read the metadata table for the simulation
    simulation_dataframe = DataFrame(sim_data.table[sim_data.slice, :])

    # Delete missing snapshots
    filter!(row -> !ismissing(row[:snapshot_paths]), simulation_dataframe)

    # Iterate over each snapshot in the slice
    iterator = eachrow(simulation_dataframe)

    # Check that there are at least 2 snapshots left
    (
        length(iterator) >= 2 ||
        throw(ArgumentError("daVirialAccretion: The given slice: $(sim_data.slice), selected for \
        less than two snapshots. I need at least two snapshots to compute the a time series of \
        gas accretion. The full simulation table is:\n$(sim_data.table)"))
    )

    ################################################################################################
    # First element of the iteration over the snapshots
    ################################################################################################

    sim_table_data = iterator[1]

    snapshot_path = sim_table_data[7]
    groupcat_path = sim_table_data[8]

    # Get the snapshot header
    snapshot_header = readSnapHeader(snapshot_path)

    # Get the group catalog header
    groupcat_header = readGroupCatHeader(groupcat_path; warnings)

    # Construct the metadata dictionary
    metadata = Dict(
        :sim_data => sim_data,
        :snap_data => Snapshot(
            snapshot_path,
            sim_table_data[1],
            1,
            sim_table_data[5],
            sim_table_data[6],
            sim_table_data[3],
            sim_table_data[4],
            snapshot_header,
        ),
        :gc_data => GroupCatalog(groupcat_path, groupcat_header),
    )

    # Read the data in the snapshot
    past_dd = merge(
        metadata,
        readSnapshot(snapshot_path, request; warnings),
        readGroupCatalog(groupcat_path, snapshot_path, request; warnings),
    )

    if tracers
        # Filter the data
        filterData!(past_dd; filter_function)
        # Translate the data
        translateData!(past_dd, translation)
    end

    # Allocate memory fo the mass axis
    Δm = Vector{Unitful.Mass}(undef, length(iterator) - 1)

    ################################################################################################
    # Iteration over the snapshots
    ################################################################################################

    @inbounds for (slice_index, sim_table_data) in pairs(iterator[2:end])

        global_index  = sim_table_data[1]
        scale_factor  = sim_table_data[3]
        redshift      = sim_table_data[4]
        physical_time = sim_table_data[5]
        lookback_time = sim_table_data[6]
        snapshot_path = sim_table_data[7]
        groupcat_path = sim_table_data[8]

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
                slice_index + 1,
                physical_time,
                lookback_time,
                scale_factor,
                redshift,
                snapshot_header,
            ),
            :gc_data => GroupCatalog(groupcat_path, groupcat_header),
        )

        # Read the data in the snapshot
        present_dd = merge(
            metadata,
            readSnapshot(snapshot_path, request; warnings),
            readGroupCatalog(groupcat_path, snapshot_path, request; warnings),
        )

        if tracers
            # Filter the data
            filterData!(present_dd; filter_function)
            # Translate the data
            translateData!(present_dd, translation)
        end

        if tracers

            Δm[slice_index], _, _ = computeVirialAccretion(present_dd, past_dd; halo_idx)

        else

            if isempty(past_dd[:group]["G_M_Crit200"])
                m_past = 0.0u"Msun"
            else
                m_past = past_dd[:group]["G_M_Crit200"][halo_idx]
            end

            if isempty(present_dd[:group]["G_M_Crit200"])
                m_present = 0.0u"Msun"
            else
                m_present = present_dd[:group]["G_M_Crit200"][halo_idx]
            end

            Δm[slice_index] = m_present -  m_past

        end

        past_dd = present_dd

    end

    # Compure the time ticks
    t  = sim_data.table[sim_data.slice, :physical_times]

    # Compure the time axis
    Δt = deltas(t)[2:end]

    if iszero(smooth)
        return t[2:end], Δm ./ Δt
    else
        return smoothWindow(t[2:end], Δm ./ Δt, smooth)
    end

end

"""
    daDiscAccretion(
        sim_data::Simulation;
        <keyword arguments>
    )::NTuple{2,Vector{<:Number}}

Compute the evolution of the accreted mass into the disc.

# Arguments

  - `sim_data::Simulation`: Information about the simulation in a [`Simulation`](@ref) object.
  - `filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all`: Which cells/particles will be plotted. The options are:

      + `:all`             -> Consider every cell/particle within the simulation box.
      + `:halo`            -> Consider only the cells/particles that belong to the main halo.
      + `:subhalo`         -> Consider only the cells/particles that belong to the main subhalo.
      + `:sphere`          -> Consider only the cell/particle inside a sphere with radius `DISK_R` (see `./src/constants.jl`).
      + `:stellar_subhalo` -> Consider only the cells/particles that belong to the main subhalo.
      + `:all_subhalo`     -> Plot every cell/particle centered around the main subhalo.
      + A dictionary with three entries:

          + `:filter_function` -> The filter function.
          + `:translation`     -> Translation for the simulation box. The posibilities are:

              + `:global_cm`                  -> Selects the center of mass of the whole system as the new origin.
              + `:{component}`                -> Sets the center of mass of the given component (e.g. :stars, :gas, :halo, etc) as the new origin. It can be any of the keys of [`PARTICLE_INDEX`](@ref).
              + `(halo_idx, subhalo_rel_idx)` -> Sets the position of the potencial minimum for the `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo) as the new origin.
              + `(halo_idx, 0)`               -> Sets the center of mass of the `halo_idx::Int` halo as the new origin.
              + `subhalo_abs_idx`             -> Sets the center of mass of the `subhalo_abs_idx::Int` as the new origin.
          + `:rotation`        -> Rotation for the simulation box. The posibilities are:

              + `:zero`                       -> No rotation is appplied.
              + `:global_am`                  -> Sets the angular momentum of the whole system as the new z axis.
              + `:stellar_am`                 -> Sets the stellar angular momentum as the new z axis.
              + `:stellar_pa`                 -> Sets the stellar principal axis as the new coordinate system.
              + `:stellar_subhalo_pa`         -> Sets the principal axis of the stars in the main subhalo as the new coordinate system.
              + `(halo_idx, subhalo_rel_idx)` -> Sets the principal axis of the stars in `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo), as the new coordinate system.
              + `(halo_idx, 0)`               -> Sets the principal axis of the stars in the `halo_idx::Int` halo, as the new coordinate system.
              + `subhalo_abs_idx`             -> Sets the principal axis of the stars in the `subhalo_abs_idx::Int` subhalo as the new coordinate system.
  - `max_r::Unitful.Length=DISK_R`: Radius of the cylinder.
  - `max_z::Unitful.Length=5.0u"kpc"`: Half height of the cylinder.
  - `smooth::Int=0`: The time series will be smooth out using `smooth` bins. Set it to 0 if you want no smoothing.
  - `warnings::Bool=true`: If a warning will be given when there is missing data.

# Returns

  - A Tuple with two elements:

      + A Vector with the physical times.
      + A Vector with the accreted mass at each time.
"""
function daDiscAccretion(
    sim_data::Simulation;
    filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all,
    max_r::Unitful.Length=DISK_R,
    max_z::Unitful.Length=5.0u"kpc",
    smooth::Int=0,
    warnings::Bool=true,
)::NTuple{2,Vector{<:Number}}

    filter_function, translation, rotation, request = selectFilter(
        filter_mode,
        Dict(
            :gas         => ["ID  ", "MASS"],
            :stars       => ["ID  ", "MASS"],
            :black_hole  => ["ID  ", "MASS"],
            :tracer      => ["PAID", "TRID"],
        ),
    )

    # Read the metadata table for the simulation
    simulation_dataframe = DataFrame(sim_data.table[sim_data.slice, :])

    # Delete missing snapshots
    filter!(row -> !ismissing(row[:snapshot_paths]), simulation_dataframe)

    # Iterate over each snapshot in the slice
    iterator = eachrow(simulation_dataframe)

    # Check that there are at least 2 snapshots left
    (
        length(iterator) >= 2 ||
        throw(ArgumentError("daDiscAccretion: The given slice: $(sim_data.slice), selected for \
        less than two snapshots. I need at least two snapshots to compute the a time series of \
        gas accretion. The full simulation table is:\n$(sim_data.table)"))
    )

    ################################################################################################
    # First element of the iteration over the snapshots
    ################################################################################################

    sim_table_data = iterator[1]

    snapshot_path = sim_table_data[7]
    groupcat_path = sim_table_data[8]

    # Get the snapshot header
    snapshot_header = readSnapHeader(snapshot_path)

    # Get the group catalog header
    groupcat_header = readGroupCatHeader(groupcat_path; warnings)

    # Construct the metadata dictionary
    metadata = Dict(
        :sim_data => sim_data,
        :snap_data => Snapshot(
            snapshot_path,
            sim_table_data[1],
            1,
            sim_table_data[5],
            sim_table_data[6],
            sim_table_data[3],
            sim_table_data[4],
            snapshot_header,
        ),
        :gc_data => GroupCatalog(groupcat_path, groupcat_header),
    )

    # Read the data in the snapshot
    past_dd = merge(
        metadata,
        readSnapshot(snapshot_path, request; warnings),
        readGroupCatalog(groupcat_path, snapshot_path, request; warnings),
    )

    # Filter the data
    filterData!(past_dd; filter_function)

    # Translate the data
    translateData!(past_dd, translation)

    # Rotate the data
    rotateData!(past_dd, rotation)

    # Allocate memory fo the mass axis
    Δm = Vector{Unitful.Mass}(undef, length(iterator) - 1)

    ################################################################################################
    # Iteration over the snapshots
    ################################################################################################

    @inbounds for (slice_index, sim_table_data) in pairs(iterator[2:end])

        global_index  = sim_table_data[1]
        scale_factor  = sim_table_data[3]
        redshift      = sim_table_data[4]
        physical_time = sim_table_data[5]
        lookback_time = sim_table_data[6]
        snapshot_path = sim_table_data[7]
        groupcat_path = sim_table_data[8]

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
                slice_index + 1,
                physical_time,
                lookback_time,
                scale_factor,
                redshift,
                snapshot_header,
            ),
            :gc_data => GroupCatalog(groupcat_path, groupcat_header),
        )

        # Read the data in the snapshot
        present_dd = merge(
            metadata,
            readSnapshot(snapshot_path, request; warnings),
            readGroupCatalog(groupcat_path, snapshot_path, request; warnings),
        )

        # Filter the data
        filterData!(present_dd; filter_function)

        # Translate the data
        translateData!(present_dd, translation)

        # Rotate the data
        rotateData!(present_dd, rotation)

        Δm[slice_index], _, _ = computeDiscAccretion(present_dd, past_dd; max_r, max_z)

        past_dd = present_dd

    end

    # Compure the time ticks
    t  = sim_data.table[sim_data.slice, :physical_times]

    # Compure the time axis
    Δt = deltas(t)[2:end]

    if iszero(smooth)
        return t[2:end], Δm ./ Δt
    else
        return smoothWindow(t[2:end], Δm ./ Δt, smooth)
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
  - `x_quantity::Symbol`: Quantity for the x axis. The options are:

      + `:scale_factor`  -> Scale factor.
      + `:redshift`      -> Redshift.
      + `:physical_time` -> Physical time since the Big Bang.
      + `:lookback_time` -> Physical time left to reach the last snapshot.
  - `y_quantity::Symbol`: Quantity for the y axis. The options are:

      + `:stellar_mass` -> Stellar mass.
      + `:sfr`          -> The star formation rate.
  - `smooth::Int=0`: The result will be smooth out using `smooth` bins. Set it to 0 if you want no smoothing.
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
    smooth::Int=0,
    warnings::Bool=true,
)::NTuple{2,Vector{<:Number}}

    snapshot_paths = filter(!ismissing, sim_data.table[!, 7])

    (
        !isempty(snapshot_paths) ||
        throw(ArgumentError("daSFRtxt: I couldn't find any snapshots in $(sim_data.path), \
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
            @warn("daSFRtxt: For non-cosmological simulations `x_quantity` can only be \
            :physical_time")
        )

        x_axis = time_ticks

    elseif x_quantity == :redshift

        if sim_data.cosmological
            x_axis = (1.0 ./ time_ticks) .- 1.0
        else
            @warn("daSFRtxt: For non-cosmological simulations `x_quantity` can only be \
            :physical_time")
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

        throw(ArgumentError("daSFRtxt: `x_quantity` can only be :scale_factor, :redshift, \
        :physical_time or :lookback_time, but I got :$(x_quantity)"))

    end

    if y_quantity == :stellar_mass

        y_axis = sfr_txt_data[6]

    elseif y_quantity == :sfr

        y_axis = sfr_txt_data[3]

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

  - `sim_data::Simulation`: Information about the simulation in a [`Simulation`](@ref) object.
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
  - `smooth::Int=0`: The result will be smooth out using `smooth` bins. Set it to 0 if you want no smoothing.
  - `warnings::Bool=true`: If a warning will be given when the target process is missing.

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
    smooth::Int=0,
    warnings::Bool=true,
)::NTuple{2,Vector{<:Number}}

    snapshot_paths = filter(!ismissing, sim_data.table[!, 7])

    (
        !isempty(snapshot_paths) ||
        throw(ArgumentError("daCPUtxt: I couldn't find any snapshots in $(sim_data.path), \
        and I need at least one for unit conversion"))
    )

    # Find the path to one snapshot
    snapshot_path = first(snapshot_paths)

    # Read its header
    header = readSnapHeader(snapshot_path)

    # Read the data in the `sfr.txt` file
    cpu_txt_data = readCpuFile(joinpath(sim_data.path, CPU_REL_PATH), [target]; warnings)[target]

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

        y_axis = cpu_txt_data[:, 1]

    elseif y_quantity == :physical_time

        if sim_data.cosmological
            y_axis = computeTime(cpu_txt_data[:, 2], header)
        else
            y_axis = cpu_txt_data[:, 2] .* internalUnits("CLKT", snapshot_path)
        end

    elseif y_quantity == :clock_time_s

        y_axis = cpu_txt_data[:, 3] .* u"s"

    elseif y_quantity == :clock_time_percent

        y_axis = cpu_txt_data[:, 4]

    elseif y_quantity == :tot_clock_time_s

        y_axis = cpu_txt_data[:, 5] .* u"s"

    elseif y_quantity == :tot_clock_time_percent

        y_axis = cpu_txt_data[:, 6]

    else

        throw(ArgumentError("daCPUtxt: I don't recognize the y_quantity :$(y_quantity)"))

    end

    # Apply smoothing if required
    if !iszero(smooth)
        x_axis, y_axis = smoothWindow(x_axis, y_axis, smooth)
    end

    return x_axis, y_axis

end
