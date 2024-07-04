####################################################################################################
# Computation of derived quantities.
####################################################################################################

"""
    computeCenter(data_dict::Dict, subfind_idx::NTuple{2,Int})::Vector{<:Unitful.Length}

Read the position of the potencial minimum for a given halo or subhalo.

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
  - `subfind_idx::NTuple{2,Int}`: Tuple with two elements:

      + Index of the target halo (FoF group). Starts at 1.
      + Index of the target subhalo (subfind), relative the target halo. Starts at 1. If set to 0, the potencial minimum of the whole halo with index `halo_idx` is returned.

# Returns

  - The specified potencial minimum.
"""
function computeCenter(data_dict::Dict, subfind_idx::NTuple{2,Int})::Vector{<:Unitful.Length}

    # If there are no subfind data, return the origin
    if ismissing(data_dict[:gc_data].path) && !isSubfindActive(data_dict[:gc_data].path)
        return zeros(typeof(1.0u"kpc"), 3)
    end

    halo_idx, subhalo_rel_idx = subfind_idx

    # Load the necessary data
    n_subhalos_in_halo = data_dict[:group]["G_Nsubs"]
    g_pos = data_dict[:group]["G_Pos"]
    s_pos = data_dict[:subhalo]["S_Pos"]

    # Check that the requested halo index is within bounds
    n_halos = data_dict[:gc_data].header.n_groups_total

    (
        !iszero(n_halos) && !any(isempty, [n_subhalos_in_halo, g_pos, s_pos]) ||
        return zeros(typeof(1.0u"kpc"), 3)
    )

    (
        0 < halo_idx <= n_halos ||
        throw(ArgumentError("computeCenter: There is only $(n_halos) FoF goups in \
        $(data_dict[:gc_data].path), so halo_idx = $(halo_idx) is out of bounds"))
    )

    # Select the halo potencial minimum if `subhalo_rel_idx` == 0
    isPositive(subhalo_rel_idx) || return g_pos[:, halo_idx]

    # Check that the requested subhalo index is within bounds
    n_subfinds = n_subhalos_in_halo[halo_idx]
    (
        subhalo_rel_idx <= n_subfinds ||
        throw(ArgumentError("computeCenter: There is only $(n_subfinds) subhalos for the FoF \
        group $(halo_idx) in $(data_dict[:gc_data].path), so subhalo_rel_idx = \
        $(subhalo_rel_idx) is out of bounds"))
    )

    # Compute the number of subhalos and particles up to the last halo before `halo_idx`
    if isone(halo_idx)
        n_subs_floor = 0
    else
        n_subs_floor = sum(n_subhalos_in_halo[1:(halo_idx - 1)]; init=0)
    end

    # Compute the subhalo absolute index
    subhalo_abs_idx = n_subs_floor + subhalo_rel_idx

    # Select the subhalo potencial minimum
    return s_pos[:, subhalo_abs_idx]

end

"""
    computeCenter(data_dict::Dict, subhalo_abs_idx::Int)::Vector{<:Unitful.Length}

Read the position of the potencial minimum for a given subhalo.

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
  - `subhalo_abs_idx::Int`: Absolute index of the target subhalo (subfind). Starts at 1.

# Returns

  - The specified potencial minimum.
"""
function computeCenter(data_dict::Dict, subhalo_abs_idx::Int)::Vector{<:Unitful.Length}

    # If there are no subfind data, return the origin
    if ismissing(data_dict[:gc_data].path) && !isSubfindActive(data_dict[:gc_data].path)
        return zeros(typeof(1.0u"kpc"), 3)
    end

    s_pos = data_dict[:subhalo]["S_Pos"]

    # Check that the requested subhalo index is within bounds
    n_subgroups_total = data_dict[:gc_data].header.n_subgroups_total

    !iszero(n_subgroups_total) && !isempty(s_pos) || return zeros(typeof(1.0u"kpc"), 3)

    (
        0 < subhalo_abs_idx <= n_subgroups_total ||
        throw(ArgumentError("computeCenter: There is only $(n_subgroups_total) subhalos in \
        $(data_dict[:gc_data].path), so subhalo_abs_idx = $(subhalo_abs_idx) is out of bounds"))
    )

    # Select the subhalo potencial minimum
    return s_pos[:, subhalo_abs_idx]

end

"""
    computeCenter(data_dict::Dict, cm_type::Symbol)::Vector{<:Unitful.Length}

Compute a characteristic center of mass for the system.

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
  - `cm_type::Symbol`: It can be:

      + `:global_cm`  -> Center of mass of the whole system.
      + `:stellar_cm` -> Stellar center of mass.
      + `:zero`       -> Origin.

# Returns

  - The specified center of mass.
"""
function computeCenter(data_dict::Dict, cm_type::Symbol)::Vector{<:Unitful.Length}

    if cm_type == :global_cm

        return computeGlobalCenterOfMass(data_dict)

    elseif cm_type == :stellar_cm

        return computeCenterOfMass(data_dict[:stars]["POS "], data_dict[:stars]["MASS"])

    elseif cm_type == :zero

        return zeros(typeof(1.0u"kpc"), 3)

    end

    throw(ArgumentError("computeCenter: `cm_type` can only be :global_cm, :stellar_cm or :zero \
    but I got :$(cm_type)"))

end

"""
    computeVcm(data_dict::Dict, subfind_idx::NTuple{2,Int})::Vector{<:Unitful.Velocity}

Read the velocity of the center of mass for a given halo or subhalo.

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
  - `subfind_idx::NTuple{2,Int}`: Tuple with two elements:

      + Index of the target halo (FoF group). Starts at 1.
      + Index of the target subhalo (subfind), relative the target halo. Starts at 1. If set to 0, the velocity of the whole halo with index `halo_idx` is returned.

# Returns

  - The specified velocity.
"""
function computeVcm(data_dict::Dict, subfind_idx::NTuple{2,Int})::Vector{<:Unitful.Velocity}

    # If there are no subfind data, return the origin
    if ismissing(data_dict[:gc_data].path) && !isSubfindActive(data_dict[:gc_data].path)
        return zeros(typeof(1.0u"km*s^-1"), 3)
    end

    halo_idx, subhalo_rel_idx = subfind_idx

    # Load the necessary data
    n_subhalos_in_halo = data_dict[:group]["G_Nsubs"]
    g_vel = data_dict[:group]["G_Vel"]
    s_vel = data_dict[:subhalo]["S_Vel"]

    # Check that the requested halo index is within bounds
    n_halos = data_dict[:gc_data].header.n_groups_total

    (
        !iszero(n_halos) && !any(isempty, [n_subhalos_in_halo, g_vel, s_vel]) ||
        return zeros(typeof(1.0u"km*s^-1"), 3)
    )

    (
        0 < halo_idx <= n_halos ||
        throw(ArgumentError("computeCenter: There is only $(n_halos) FoF goups in \
        $(data_dict[:gc_data].path), so halo_idx = $(halo_idx) is out of bounds"))
    )

    # Select the halo velocity if `subhalo_rel_idx` == 0
    isPositive(subhalo_rel_idx) || return g_vel[:, halo_idx]

    # Check that the requested subhalo index is within bounds
    n_subfinds = n_subhalos_in_halo[halo_idx]
    (
        subhalo_rel_idx <= n_subfinds ||
        throw(ArgumentError("computeCenter: There is only $(n_subfinds) subhalos for the FoF \
        group $(halo_idx) in $(data_dict[:gc_data].path), so subhalo_rel_idx = \
        $(subhalo_rel_idx) is out of bounds"))
    )

    # Compute the number of subhalos and particles up to the last halo before `halo_idx`
    if isone(halo_idx)
        n_subs_floor = 0
    else
        n_subs_floor = sum(n_subhalos_in_halo[1:(halo_idx - 1)]; init=0)
    end

    # Compute the subhalo absolute index
    subhalo_abs_idx = n_subs_floor + subhalo_rel_idx

    # Select the subhalo velocity
    return s_vel[:, subhalo_abs_idx]

end

"""
    computeVcm(data_dict::Dict, subhalo_abs_idx::Int)::Vector{<:Unitful.Velocity}

Read the velocity of the center of mass for a given subhalo.

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
  - `subhalo_abs_idx::Int`: Absolute index of the target subhalo (subfind). Starts at 1.

# Returns

  - The specified velocity.
"""
function computeVcm(data_dict::Dict, subhalo_abs_idx::Int)::Vector{<:Unitful.Velocity}

    # If there are no subfind data, return the origin
    if ismissing(data_dict[:gc_data].path) && !isSubfindActive(data_dict[:gc_data].path)
        return zeros(typeof(1.0u"km*s^-1"), 3)
    end

    s_vel = data_dict[:subhalo]["S_Vel"]

    # Check that the requested subhalo index is within bounds
    n_subgroups_total = data_dict[:gc_data].header.n_subgroups_total

    !iszero(n_subgroups_total) && !isempty(s_vel) || return zeros(typeof(1.0u"km*s^-1"), 3)

    (
        0 < subhalo_abs_idx <= n_subgroups_total ||
        throw(ArgumentError("computeCenter: There is only $(n_subgroups_total) subhalos in \
        $(data_dict[:gc_data].path), so subhalo_abs_idx = $(subhalo_abs_idx) is out of bounds"))
    )

    # Select the subhalo velocity
    return s_vel[:, subhalo_abs_idx]

end

"""
    computeVcm(data_dict::Dict, cm_type::Symbol)::Vector{<:Unitful.Velocity}

Compute the velocity of a characteristic center of mass for the system.

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
  - `cm_type::Symbol`: It can be:

      + `:global_cm`  -> Center of mass of the whole system.
      + `:stellar_cm` -> Stellar center of mass.

# Returns

  - The specified velocity.
"""
function computeVcm(data_dict::Dict, cm_type::Symbol)::Vector{<:Unitful.Velocity}

    if cm_type == :global_cm

        return computeGlobalVcm(data_dict)

    elseif cm_type == :stellar_cm

        return computeStellarVcm(data_dict)

    elseif cm_type == :zero

        return zeros(1.0u"km*s^-1", 3)

    end

    throw(ArgumentError("computeVcm: `cm_type` can only be :global_cm, :stellar_cm or :zero \
    but I got :$(center_type)"))

end

@doc raw"""
    computeTime(
        scale_factors::Vector{<:Real},
        header::SnapshotHeader;
        <keyword arguments>
    )::Vector{<:Unitful.Time}

Compute the physical time corresponding to each of the `scale_factors`.

To get the physical time $t$ from the scale factor `a`, one does the integral:

```math
t = \frac{1}{H_0} \int_0^a \frac{\mathrm{d}a'}{a' \, \sqrt{\mathcal{E}(a')}} \, ,
```

where

```math
\mathcal{E}(a) = \Omega_\Lambda + \Omega_m \, a^{-3} + \Omega_r \, a^{-4} + \Omega_K \, a^{-2} \, .
```

# Arguments

  - `scale_factors::Vector{<:Real}`: Scale factors.
  - `header::SnapshotHeader`: A header of the simulation, containing the cosmological parameters.
  - `a0::Float64=0.0`: Initial scale factor.

# Returns

  - A vector with the physical times.
"""
function computeTime(
    scale_factors::Vector{<:Real},
    header::SnapshotHeader;
    a0::Float64=0.0,
)::Vector{<:Unitful.Time}

    f = x -> energyIntegrand(x, header)

    return [quadgk(f, a0, a)[1] * u"Gyr" for a in scale_factors]

end

@doc raw"""
    computeTime(a::Real, header::SnapshotHeader; <keyword arguments>)::Unitful.Time

Compute the physical time corresponding to the scale factor `a`.

To get the physical time $t$ from the scale factor `a`, one does the integral:

```math
t = \frac{1}{H_0} \int_0^a \frac{\mathrm{d}a'}{a' \, \sqrt{\mathcal{E}(a')}} \, ,
```

where

```math
\mathcal{E}(a) = \Omega_\Lambda + \Omega_m \, a^{-3} + \Omega_r \, a^{-4} + \Omega_K \, a^{-2} \, .
```

# Arguments

  - `a::Real`: Scale factor.
  - `header::SnapshotHeader`: A header of the simulation, containing the cosmological parameters.
  - `a0::Float64=0.0`: Initial scale factor.

# Returns

  - The physical time.
"""
function computeTime(a::Real, header::SnapshotHeader; a0::Float64=0.0)::Unitful.Time

    return computeTime([a], header; a0)[1]

end

"""
    computeTimeTicks(
        paths::Vector{<:Union{Missing,String}},
    )::Tuple{Vector{Float64},Vector{Float64},Vector{<:Unitful.Time},Vector{<:Unitful.Time}}

Compute the different times stamps associated with each snapshot in `paths`.

# Arguments

  - `paths::Vector{<:Union{Missing,String}}`: Paths to the snapshots.

# Returns

  - A tuple with four elements:

      + A vector with the scale factors.
      + A vector with the redshifts.
      + A vector with the physical times (physical time since the Big Bang).
      + A vector with the lookback times (physical time left to reach the last snapshot).
"""
function computeTimeTicks(
    paths::Vector{<:Union{Missing,String}},
)::Tuple{Vector{Float64},Vector{Float64},Vector{<:Unitful.Time},Vector{<:Unitful.Time}}

    snapshot_paths = filter(!ismissing, paths)

    !isempty(snapshot_paths) || return [NaN], [NaN], [NaN*u"s"], [NaN*u"s"]

    first_snapshot = first(snapshot_paths)

    if isCosmological(first_snapshot)

        # For cosmological simulations, the time in the snapshot is the scale factor
        scale_factors = [readTime(path) for path in snapshot_paths]
        redshifts = @. (1.0 / scale_factors) - 1.0
        physical_times = computeTime(scale_factors, readSnapHeader(first_snapshot))
        lookback_times = last(physical_times) .- physical_times

    else

        # Compute the factor for internal units of time
        u_time = internalUnits("CLKT", first_snapshot)

        # a = 1.0 for non-cosmological simulations
        scale_factors = ones(length(snapshot_paths))
        # z = 0.0 for non-cosmological simulations
        redshifts = zeros(length(snapshot_paths))
        # For non-cosmological simulations, the time in the snapshot is the physical time
        physical_times = [readTime(path) * u_time for path in snapshot_paths]
        lookback_times = last(physical_times) .- physical_times

    end

    return scale_factors, redshifts, physical_times, lookback_times

end

"""
    computeTemperature(
        metals::Matrix{Float32},
        internal_energy::Vector{<:SpecificEnergy},
        electron_fraction::Vector{Float32},
    )::Vector{<:Unitful.Temperature}

Compute the gas temperature of a group of gas cells.

# Arguments

  - `metals::Matrix{Float32}`: Matrix with the mass content for every tracked element. Each row is an element, and each column a cell.
  - `internal_energy::Vector{<:SpecificEnergy}`: Specific internal energy of every gas cell.
  - `electron_fraction::Vector{Float32}`: Number fraction of electrons in every gas cell.

# Returns

  - The temperature of each gas cell.
"""
function computeTemperature(
    metals::Matrix{Float32},
    internal_energy::Vector{<:SpecificEnergy},
    electron_fraction::Vector{Float32},
)::Vector{<:Unitful.Temperature}

    # xH := mass_fraction_of_hydrogen
    # It should be similar to the constant [`HYDROGEN_MASSFRAC`](@ref) used in Arepo
    xH = metals[ELEMENT_INDEX[:H], :]

    # yHe := number_of_helium_atoms / number_of_hydrogen_atoms
    # Take the mass fraction of metals as negligible
    yHe = @. (1.0 - xH) / (4.0 * xH)

    # electron_fraction := number_of_electrons / number_of_hydrogen_atoms
    # μ := total_mass / (total_number_of_particles * proton_mass)
    #   ≈ number_of_protons / total_number_of_particles
    # For the total mass, take the mass of electrons as negligible
    μ = @. (1.0 + 4.0 * yHe) / (1.0 + yHe + electron_fraction)

    # T = (adiabatic_index - 1) * internal_energy_per_unit_mass *
    #     (total_mass / total_number_of_particles) / boltzmann_constant
    return @. 0.6667 * internal_energy * μ * Unitful.mp / Unitful.k

end

"""
    computeDistance(
        positions::Matrix{<:Number};
        <keyword arguments>
    )::Vector{<:Number}

Compute the distance of a group of points to `center`.

# Arguments

  - `positions::Matrix{<:Number}`: Positions of the points. Each column is a point and each row a dimension.
  - `center::Union{Vector{<:Number},Nothing}=nothing`: Origin used to compute the distances. If set to `nothing`, 0 is used.

# Returns

  - The distance of every point to `center`.
"""
function computeDistance(
    positions::Matrix{<:Number};
    center::Union{Vector{<:Number},Nothing}=nothing,
)::Vector{<:Number}

    if isempty(positions)
        return eltype(positions)[]
    end

    if center === nothing
        return [norm(col) for col in eachcol(positions)]
    end

    (
        length(center) == size(positions, 1) ||
        throw(ArgumentError("computeDistance: `center` must have as many elements as `positions` \
        has rows, but I got length(center) = $(length(center)) and size(positions, 1) = \
        $(size(positions, 1))"))
    )

    return [norm(col .- center) for col in eachcol(positions)]

end

"""
    computeCenterOfMass(
        positions::Matrix{<:Unitful.Length},
        mass::Vector{<:Unitful.Mass},
    )::Vector{<:Unitful.Length}

Compute the center of mass of a group of cells/particles.

# Arguments

  - `positions::Matrix{<:Unitful.Length}`: Positions of the cells/particles. Each column is a cell/particle and each row a dimension.
  - `masses::Vector{<:Unitful.Mass}`: Masses of the cells/particles.

# Returns

  - The center of mass.
"""
function computeCenterOfMass(
    positions::Matrix{<:Unitful.Length},
    masses::Vector{<:Unitful.Mass},
)::Vector{<:Unitful.Length}

    # Check for missing data
    !any(isempty, [positions, masses]) || return zeros(typeof(1.0u"kpc"), 3)

    # Compute the total mass
    M = sum(masses)

    # Compute the center of mass
    center_of_mass = [sum(row .* masses) / M for row in eachrow(positions)]

    return center_of_mass

end

"""
    computeGlobalCenterOfMass(data_dict::Dict)::Vector{<:Unitful.Length}

Compute the center of mass of the whole system in `data`.

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

  - The center of mass.
"""
function computeGlobalCenterOfMass(data_dict::Dict)::Vector{<:Unitful.Length}

    type_symbols = snapshotTypes(data_dict)

    filter!(ts -> !isempty(data_dict[ts]["POS "]), type_symbols)

    # Concatenate the position and masses of all the cells and particles in the system
    positions = hcat([data_dict[type_symbol]["POS "] for type_symbol in type_symbols]...)
    masses    = vcat([data_dict[type_symbol]["MASS"] for type_symbol in type_symbols]...)

    # Check for missing data
    !any(isempty, [positions, masses]) || return zeros(typeof(1.0u"kpc"), 3)

    @debug("computeGlobalCenterOfMass: The center of mass will be computed using $(type_symbols)")

    return computeCenterOfMass(positions, masses)

end

"""
    computeInertiaTensor(
        positions::Matrix{<:Unitful.Length},
        masses::Vector{<:Unitful.Mass},
    )::Matrix{Float64}

Compute the inertia tensor of a group of cells/particles.

# Arguments

  - `positions::Matrix{<:Unitful.Length}`: Positions of the cells/particles. Each column is a cell/particle and each row a dimension.
  - `masses::Vector{<:Unitful.Mass}`: Masses of the cells/particles.

# Returns

  - The inertia tensor.
"""
function computeInertiaTensor(
    positions::Matrix{<:Unitful.Length},
    masses::Vector{<:Unitful.Mass},
)::Matrix{Float64}

    # Check for missing data
    (
        !any(isempty, [positions, masses])  ||
        throw(ArgumentError("computeInertiaTensor: The inertia tensor is not defined for an \
        empty system"))
    )

    # Allocate memory
    J = zeros(Float64, (3, 3))

    # Compute the inertia tensor
    for (r, m) in zip(eachcol(positions), masses)

        x, y, z = r

        J[1, 1] += m * (y * y + z * z)
        J[2, 2] += m * (x * x + z * z)
        J[3, 3] += m * (x * x + y * y)

        J[1, 2] -= m * x * y
        J[1, 3] -= m * x * z
        J[2, 3] -= m * y * z

    end

    J[2, 1] = J[1, 2]
    J[3, 1] = J[1, 3]
    J[3, 2] = J[2, 3]

    return J

end

"""
    computeAMRotationMatrix(
        positions::Matrix{<:Unitful.Length},
        velocities::Matrix{<:Unitful.Velocity},
        masses::Vector{<:Unitful.Mass},
    )::Union{Matrix{Float64},UniformScaling{Bool}}

Compute the rotation matrix that will turn the total angular momentum into the z axis; when view as an active (alibi) trasformation.

# Arguments

  - `positions::Matrix{<:Unitful.Length}`: Positions of the cells/particles. Each column is a cell/particle and each row a dimension.
  - `velocities::Matrix{<:Unitful.Velocity}`: Velocities of the cells/particles. Each column is a cell/particle and each row a dimension.
  - `masses::Vector{<:Unitful.Mass}`: Mass of every cell/particle.

# Returns

  - The rotation matrix.
"""
function computeAMRotationMatrix(
    positions::Matrix{<:Unitful.Length},
    velocities::Matrix{<:Unitful.Velocity},
    masses::Vector{<:Unitful.Mass},
)::Union{Matrix{Float64},UniformScaling{Bool}}

    # Check for missing data
    !any(isempty, [positions, velocities, masses]) || return I

    # Compute the total angular momentum
    L = computeTotalAngularMomentum(positions, velocities, masses)

    # Rotation vector
    n = [L[2], -L[1], 0.0]

    # Angle of rotation
    θ = acos(L[3])

    return Matrix{Float64}(AngleAxis(θ, n...))

end

"""
    computePARotationMatrix(
        positions::Matrix{<:Unitful.Length},
        velocities::Matrix{<:Unitful.Velocity},
        masses::Vector{<:Unitful.Mass},
    )::Union{Matrix{Float64},UniformScaling{Bool}}

Compute the rotation matrix that will turn the pricipal axis into the new coordinate system; when view as an pasive (alias) trasformation.

# Arguments

  - `positions::Matrix{<:Unitful.Length}`: Positions of the cells/particles. Each column is a cell/particle and each row a dimension.
  - `velocities::Matrix{<:Unitful.Velocity}`: Velocities of the cells/particles. Each column is a cell/particle and each row a dimension.
  - `masses::Vector{<:Unitful.Mass}`: Mass of every cell/particle.

# Returns

  - The rotation matrix.
"""
function computePARotationMatrix(
    positions::Matrix{<:Unitful.Length},
    velocities::Matrix{<:Unitful.Velocity},
    masses::Vector{<:Unitful.Mass},
)::Union{Matrix{Float64},UniformScaling{Bool}}

    # Check for missing data
    !isempty(positions) || return I

    # Principal axis operator
    R = ustrip.(positions * positions')

    # Reverse the order of the eigenvectors, making the last column the eigenvector
    # with the largest eigenvalue, which should correspond to the new z axis
    pa = eigvecs(R)[:, end:-1:1]

    # Compute the total angular momentum
    L = computeTotalAngularMomentum(positions, velocities, masses; normal=true)

    # 3rd principal axis ≡ new z axis
    pa_z = pa[:, 3]

    # Rotate the principal axis as to align the thid component with the angular momentum
    θ = acos(L ⋅ pa_z)
    n = cross(pa_z, L)
    aligned_pa = AngleAxis(θ, n...) * pa

    # The rotation matrix is made from the principal axis as rows
    rotation_matrix = aligned_pa'

    if det(rotation_matrix) < 0.0
        # If the determinant is < 0, that means that the chosen principal axis for the x and y
        # directions form a left-handed cartesian reference system (x × y = -z). When applying
        # this as a rotation, the z axis will be flipped. So, in this case we swap the x and y
        # axis to get a right-handed cartesian reference system (x × y = z) and generate the
        # correct rotation
        rotation_matrix[1, :], rotation_matrix[2, :] = rotation_matrix[2, :], rotation_matrix[1, :]
    end

    return Matrix{Float64}(rotation_matrix)

end

"""
    computeGlobalAMRotationMatrix(data_dict::Dict)::Union{Matrix{Float64},UniformScaling{Bool}}

Compute the rotation matrix that will turn the total angular momentum, of the whole system in `data`, into the z axis; when view as an active (alibi) trasformation.

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

  - The rotation matrix.
"""
function computeGlobalAMRotationMatrix(data_dict::Dict)::Union{Matrix{Float64},UniformScaling{Bool}}

    type_symbols = snapshotTypes(data_dict)

    filter!(ts -> !isempty(data_dict[ts]["POS "]), type_symbols)

    # Concatenate the positions, velocities, and masses of all the cells and particles in the system
    positions  = hcat([data_dict[type_symbol]["POS "] for type_symbol in type_symbols]...)
    velocities = hcat([data_dict[type_symbol]["VEL "] for type_symbol in type_symbols]...)
    masses     = vcat([data_dict[type_symbol]["MASS"] for type_symbol in type_symbols]...)

    # Check for missing data
    !any(isempty, [positions, velocities, masses]) || return I

    @debug("computeGlobalAMRotationMatrix: The rotation matrix will be computed \
    using $(type_symbols)")

    return computeAMRotationMatrix(positions, velocities, masses)

end

"""
    computeAngularMomentum(
        positions::Matrix{<:Unitful.Length},
        velocities::Matrix{<:Unitful.Velocity},
        masses::Vector{<:Unitful.Mass},
    )::Vector{Vector{<:AngularMomentum}}

Compute the angular momentum of each cell/particle, with respect to the origin.

# Arguments

  - `positions::Matrix{<:Unitful.Length}`: Positions of the cells/particles. Each column is a cell/particle and each row a dimension.
  - `velocities::Matrix{<:Unitful.Velocity}`: Velocities of the cells/particles. Each column is a cell/particle and each row a dimension.
  - `masses::Vector{<:Unitful.Mass}`: Mass of every cell/particle.

# Returns

  - The angular momentum of each cell/particle.
"""
function computeAngularMomentum(
    positions::Matrix{<:Unitful.Length},
    velocities::Matrix{<:Unitful.Velocity},
    masses::Vector{<:Unitful.Mass},
)::Vector{Vector{<:AngularMomentum}}

    # Check for missing data
    !any(isempty, [positions, velocities, masses]) || return [0.0, 0.0, 1.0]

    iterator = zip(masses, eachcol(positions), eachcol(velocities))

    return map(x -> x[1] .* cross(x[2], x[3]), iterator)

end

"""
    computeTotalAngularMomentum(
        positions::Matrix{<:Unitful.Length},
        velocities::Matrix{<:Unitful.Velocity},
        masses::Vector{<:Unitful.Mass};
        <keyword arguments>
    )::Vector{<:Number}

Compute the total angular momentum of a group of cells/particles, with respect to the origin.

# Arguments

  - `positions::Matrix{<:Unitful.Length}`: Positions of the cells/particles. Each column is a cell/particle and each row a dimension.
  - `velocities::Matrix{<:Unitful.Velocity}`: Velocities of the cells/particles. Each column is a cell/particle and each row a dimension.
  - `masses::Vector{<:Unitful.Mass}`: Mass of every cell/particle.
  - `normal::Bool=true`: If the result will be normalized.

# Returns

  - The angular momentum.
"""
function computeTotalAngularMomentum(
    positions::Matrix{<:Unitful.Length},
    velocities::Matrix{<:Unitful.Velocity},
    masses::Vector{<:Unitful.Mass};
    normal::Bool=true,
)::Vector{<:Number}

    # Check for missing data
    !any(isempty, [positions, velocities, masses]) || return [0.0, 0.0, 1.0]

    iterator = zip(masses, eachcol(positions), eachcol(velocities))

    # Compute the total angular momentum
    L = mapreduce(x -> x[1] .* cross(x[2], x[3]), +, iterator)

    return normal ? normalize!(ustrip.(L)) : L

end

"""
    computeGlobalAngularMomentum(data_dict::Dict; <keyword arguments>)::Vector{<:Number}

Compute the total angular momentum with respect to the origin of the whole system in `data`.

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
  - `normal::Bool=true`: If the result will be normalized.

# Returns

  - The angular momentum.
"""
function computeGlobalAngularMomentum(data_dict::Dict; normal::Bool=true)::Vector{<:Number}

    type_symbols = snapshotTypes(data_dict)

    filter!(ts -> !isempty(data_dict[ts]["POS "]), type_symbols)

    # Concatenate the position, velocities, and masses of all the cells and particles in the system
    positions  = hcat([data_dict[type_symbol]["POS "] for type_symbol in type_symbols]...)
    velocities = hcat([data_dict[type_symbol]["VEL "] for type_symbol in type_symbols]...)
    masses     = vcat([data_dict[type_symbol]["MASS"] for type_symbol in type_symbols]...)

    # Check for missing data
    !any(isempty, [positions, velocities, masses]) || return [0.0, 0.0, 1.0]

    @debug("computeGlobalAngularMomentum: The angular momentum will be computed \
    using $(type_symbols)")

    return computeTotalAngularMomentum(positions, velocities, masses; normal)

end

@doc raw"""
    computeSpinParameter(
        positions::Matrix{<:Unitful.Length},
        masses::Vector{<:Unitful.Mass},
        velocities::Matrix{<:Unitful.Velocity};
        <keyword arguments>
    )::Float64

Compute the spin parameter for a system of cells/particles, with respect to the origin.

The spin parameter was originally defined by Peebles (1969) as,

```math
\lambda = \frac{J \, \sqrt{E}}{G \, M^{5/2}} \, ,
```

where $J$ is the norm of the total angular momentum, $M$ the total mass, $G$ the gravitational constant, and

```math
E = |E_P + E_k| \, ,
```

where $E_P$ is the total potencial energy and $E_k$ is the total kinetic energy (including thermal energy of the gas).

Due to the computational complexity of calculating $E_P$ for a large group of particles, Bullock et al. (2001) proposed an alternative definition of the spin parameter,

```math
\lambda = \frac{J}{\sqrt{2} \, M \, R \, V} \, ,
```

where $J$ is the norm of the total angular momentum inside a sphere of radius $R$ containing mass $M$, and

```math
V = \sqrt{\frac{G \, M}{R}} \, ,
```

is the circular velocity.

# Arguments

  - `positions::Matrix{<:Unitful.Length}`: Positions of the cells/particles. Each column is a cell/particle and each row a dimension.
  - `velocities::Matrix{<:Unitful.Velocity}`: Velocities of the cells/particles. Each column is a cell/particle and each row a dimension.
  - `masses::Vector{<:Unitful.Mass}`: Mass of every cell/particle.
  - `R::Unitful.Length=FILTER_R`: Radius.

# Returns

  - The spin parameter.

# References

P. J. E. Peebles (1969). *Origin of the Angular Momentum of Galaxies*. Astrophysical Journal, **155**, 393. [doi:10.1086/149876](https://doi.org/10.1086/149876)

J. S. Bullock et al. (2001). *A Universal Angular Momentum Profile for Galactic Halos*. The Astrophysical Journal, **555(1)**, 240. [doi:10.1086/321477](https://doi.org/10.1086/321477)

J. Zjupa et al. (2017). *Angular momentum properties of haloes and their baryon content in the Illustris simulation*. Monthly Notices of the Royal Astronomical Society, **466(2)**, 1625–1647. [doi:10.1093/mnras/stw2945](https://doi.org/10.1093/mnras/stw2945)
"""
function computeSpinParameter(
    positions::Matrix{<:Unitful.Length},
    velocities::Matrix{<:Unitful.Velocity},
    masses::Vector{<:Unitful.Mass};
    R::Unitful.Length=FILTER_R,
)::Float64

    (
        !any(isempty, [positions, velocities, masses]) ||
        throw(ArgumentError("computeSpinParameter: The spin parameter of an empty dataset \
        is undefined"))
    )

    # Find cells/particles within the virial radius
    idx = map(x -> x <= R, computeDistance(positions))

    # Compute the total mass within the virial radius
    M = sum(masses[idx]; init=0.0u"Msun")

    # Compute the norm of the total angular momentum
    J = norm(
        computeTotalAngularMomentum(
            positions[:, idx],
            velocities[:, idx],
            masses[idx];
            normal=false,
        )
    )

    return uconvert(Unitful.NoUnits, J / sqrt(2.0 * R * Unitful.G * M^3))

end

"""
    computeGlobalSpinParameter(data_dict::Dict; <keyword arguments>)::Float64

Compute the spin parameter of the whole system in `data`.

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
  - `R::Unitful.Length=FILTER_R`: Radius.

# Returns

  - The spin parameter.
"""
function computeGlobalSpinParameter(data_dict::Dict; R::Unitful.Length=FILTER_R)::Float64

    type_symbols = snapshotTypes(data_dict)

    filter!(ts -> !isempty(data_dict[ts]["POS "]), type_symbols)

    # Concatenate the position and masses of all the cells and particles in the system
    positions  = hcat([data_dict[type_symbol]["POS "] for type_symbol in type_symbols]...)
    velocities = hcat([data_dict[type_symbol]["VEL "] for type_symbol in type_symbols]...)
    masses     = vcat([data_dict[type_symbol]["MASS"] for type_symbol in type_symbols]...)

    @debug("computeGlobalSpinParameter: The spin parameter will be computed using $(type_symbols)")

    # Compute the total spin parameter
    return computeSpinParameter(positions, velocities, masses; R)

end

"""
    computeStellarVcm(data_dict::Dict)::Vector{<:Unitful.Velocity}

Compute the velocity of the stellar center of mass.

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

  - The velocity of the stellar center of mass.
"""
function computeStellarVcm(data_dict::Dict)::Vector{<:Unitful.Velocity}

    # Load the necessary data
    velocities = data_dict[:stars]["VEL "]
    masses = data_dict[:stars]["MASS"]

    # Check for missing data
    !any(isempty, [velocities, masses]) || return zeros(typeof(1.0u"km/s"), 3)

    # Compute the total mass
    M = sum(masses)

    # Compute the velocity of the center of mass
    vcm = [sum(row .* masses) / M for row in eachrow(velocities)]

    return vcm

end

"""
    computeGlobalVcm(data_dict::Dict)::Vector{<:Unitful.Velocity}

Compute the velocity of the global center of mass.

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

  - The velocity of the global center of mass.
"""
function computeGlobalVcm(data_dict::Dict)::Vector{<:Unitful.Velocity}

    type_symbols = snapshotTypes(data_dict)

    filter!(ts -> !isempty(data_dict[ts]["VEL "]), type_symbols)

    # Load the necessary data
    velocities = hcat([data_dict[type_symbol]["VEL "] for type_symbol in type_symbols]...)
    masses     = vcat([data_dict[type_symbol]["MASS"] for type_symbol in type_symbols]...)

    # Check for missing data
    !any(isempty, [velocities, masses]) || return zeros(typeof(1.0u"km/s"), 3)

    # Compute the total mass
    M = sum(masses)

    # Compute the velocity of the center of mass
    vcm = [sum(row .* masses) / M for row in eachrow(velocities)]

    return vcm

end

@doc raw"""
    computeStellarVcirc(
        data_dict::Dict,
    )::Tuple{Vector{<:Unitful.Length},Vector{<:Unitful.Velocity}}

Compute the circular velocity of each stellar particle, with respect to the origin.

The circular velocity of a star is,

```math
v_\mathrm{circ} = \sqrt{\frac{\mathrm{G} \, M(r)}{r}} \, ,
```

where $r$ is the radial distance of the star, and $M(r)$ is the total mass within a sphere of radius $r$.

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

  - A tuple with two elements:

      + A vector with the radial distance of each star to the origin.
      + A vector with the circular velocity of each star.
"""
function computeStellarVcirc(
    data_dict::Dict,
)::Tuple{Vector{<:Unitful.Length},Vector{<:Unitful.Velocity}}

    # Compute the radial distance to each star
    rs = computeDistance(data_dict[:stars]["POS "])

    # Check for missing data
    !isempty(rs) || return rs, Unitful.Velocity[]

    type_symbols = filter!(
        ts -> !isempty(data_dict[ts]["POS "]),
        [:stars, :gas, :halo, :black_hole],
    )

    # Concatenate the position and masses of all the cells and particles in the system
    distances = vcat(
        [computeDistance(data_dict[type_symbol]["POS "]) for
        type_symbol in type_symbols]...,
    )
    masses = vcat([data_dict[type_symbol]["MASS"] for type_symbol in type_symbols]...)

    # Use the radial distances as bin edges for the mass histogram
    edges = [0.0u"kpc", rs...]

    # Compute to total mass within each stellar radial distance
    M = similar(rs, eltype(masses))
    cumsum!(M, histogram1D(distances, masses, edges; empty_nan=false))

    # The mass histogram is a sorted array, so it is reverted to the unsorted order of `r`
    # to make `vcirc` the circular velocity of each star in the order of the snapshot
    invpermute!(M, sortperm(rs))

    @debug("computeStellarVcirc: The circular velocity will be computed using $(type_symbols)")

    vcirc = [iszero(r) ? 0.0u"km*s^-1" : sqrt(Unitful.G * m / r) for (m, r) in zip(M, rs)]

    return rs, vcirc

end

@doc raw"""
    computeStellarCircularity(data_dict::Dict)::Vector{Float64}

Compute the circularity of each stellar particle, with respect to the origin and the $z$ direction [0, 0, 1].

The circularity of a star is,

```math
\epsilon = j_z / j_\mathrm{circ} \, ,
```

where $j_z$ is the $z$ component of its specific angular momentum, and $j_\mathrm{circ}$ is the specific angular momentum of a circular orbit,

```math
j_\mathrm{circ} = r \, v_\mathrm{circ} = \sqrt{\mathrm{G} \, r \, M(r)} \, ,
```

where $r$ is the radial distance of the star, and $M(r)$ is the total mass within a sphere of radius $r$.

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

  - The circularity $\epsilon$ of each star.
"""
function computeStellarCircularity(data_dict::Dict)::Vector{Float64}

    # Load the necessary data
    positions = data_dict[:stars]["POS "]
    velocities = data_dict[:stars]["VEL "]

    # Check for missing data
    !any(isempty, [positions, velocities]) || return Float64[]

    # Compute the specific angular momentum in the z direction
    jzs = [x[1] * v[2] - x[2] * v[1] for (x, v) in zip(eachcol(positions), eachcol(velocities))]

    # Compute the circular velocities and the radial distances
    rs, vcircs = computeStellarVcirc(data_dict)

    stellar_circularity = [
        any(iszero, [r, vcirc]) ? 0.0 : ustrip(Unitful.NoUnits, jz / (r * vcirc)) for
        (jz, r, vcirc) in zip(jzs, rs, vcircs)
    ]

    return stellar_circularity

end

@doc raw"""
    computeStellarVpolar(data_dict::Dict, component::Symbol)::Vector{<:Unitful.Velocity}

Compute the cylindrical components of the velocity, $\mathbf{\vec{v}} = v_r \, \mathbf{e_r} + v_\theta \, \mathbf{e_\theta} + v_z \, \mathbf{e_z}$.

The speed in the radial direction expressed in Cartesian coordinates is

```math
v_r = \frac{x \, v_x + y \, v_y}{\sqrt(x^2 + y^2)} \, ,
```

in the tangential direction is

```math
v_\tau = \frac{x \, v_y - y \, v_x}{\sqrt(x^2 + y^2)} \, ,
```

and the speed in the z direction will be computes as

```math
v^*_z = v_z \, \mathrm{sign}(z) \, ,
```

in order to distinguish between inflows ($v^*_z < 0$) and outflows ($v^*_z > 0$).

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
  - `component::Symbol`: Which component will be calculated. The options are:

      + `:radial`     -> Stellar radial speed ($v_r$).
      + `:tangential` -> Stellar tangential speed ($v_\theta$).
      + `:zstar`      -> Stellar speed in the z direction, computed as $v_z \, \mathrm{sign}(z)$.

# Returns

  - The chosen cylindricall component of the velocity.
"""
function computeStellarVpolar(data_dict::Dict, component::Symbol)::Vector{<:Unitful.Velocity}

    # Load the necessary data
    positions = data_dict[:stars]["POS "]
    velocities = data_dict[:stars]["VEL "]

    x = positions[1, :]
    y = positions[2, :]
    z = positions[3, :]

    vx = velocities[1, :]
    vy = velocities[2, :]
    vz = velocities[3, :]

    if component == :radial

        # Compute the radial component
        vp = @. (x * vx + y * vy) / sqrt(x^2 + y^2)

    elseif component == :tangential

        # Compute the tangential component
        vp = @. (x * vy - y * vx) / sqrt(x^2 + y^2)

    elseif component == :zstar

        # Compute the z component
        vp = @. vz * sign(z)

    else

        throw(ArgumentError("computeStellarVpolar: `component` can only be :radial, :tangential \
        or :zstar, but I got :$(component)"))

    end

    return vp

end

"""
    computeMassRadius(
        positions::Matrix{<:Unitful.Length},
        masses::Vector{<:Unitful.Mass};
        <keyword arguments>
    )::Unitful.Length

Compute the radius containing `percet`% of the total mass.

# Arguments

  - `positions::Matrix{<:Unitful.Length}`: Positions of the cells/particles. Each column is a cell/particle and each row a dimension.
  - `masses::Vector{<:Unitful.Mass}`: Masses of the cells/particles.
  - `percent::Float64=90.0`: Target percentage of the total mass.

# Returns

  - The radius containing `percet`% of the total mass.
"""
function computeMassRadius(
    positions::Matrix{<:Unitful.Length},
    masses::Vector{<:Unitful.Mass};
    percent::Float64=90.0,
)::Unitful.Length

    (
        0 < percent <=100  ||
        throw(ArgumentError("computeMassRadius: The argument `percent` must be between 0 and 100, \
        but I got $(percent)"))
    )

    # Check for missing data
    !any(isempty, [positions, masses]) || return zero(typeof(1.0u"kpc"))

    # Compute the mass limit
    mass_limit = sum(masses) * (percent / 100.0)

    # Compute the radial distance of each cell/particle
    radial_distances = computeDistance(positions)

    sort_idxs = sortperm(radial_distances)

    # Find the mass radius
    accu_mass = 0.0u"Msun"
    target_idx = 0
    for mass in masses[sort_idxs]
        accu_mass += mass
        accu_mass < mass_limit || break
        target_idx += 1
    end

    return radial_distances[sort_idxs[target_idx]]

end

"""
    computeMassHeight(
        positions::Matrix{<:Unitful.Length},
        masses::Vector{<:Unitful.Mass};
        <keyword arguments>
    )::Unitful.Length

Compute the total height of a cylinder, of infinite radius, containing `percet`% of the total mass.

# Arguments

  - `positions::Matrix{<:Unitful.Length}`: Positions of the cells/particles. Each column is a cell/particle and each row a dimension.
  - `masses::Vector{<:Unitful.Mass}`: Masses of the cells/particles.
  - `percent::Float64=90.0`: Target percentage of the total mass.

# Returns

  - The height containing `percet`% of the total mass.
"""
function computeMassHeight(
    positions::Matrix{<:Unitful.Length},
    masses::Vector{<:Unitful.Mass};
    percent::Float64=90.0,
)::Unitful.Length

    (
        0 < percent <=100  ||
        throw(ArgumentError("computeMassHeight: The argument `percent` must be between 0 and 100, \
        but I got $(percent)"))
    )

    # Check for missing data
    !any(isempty, [positions, masses]) || return zero(typeof(1.0u"kpc"))

    # Compute the mass limit
    mass_limit = sum(masses) * (percent / 100.0)

    # Compute the vertical separation of each cell/particle
    heights = abs.(vec(positions[3, :]))

    sort_idxs = sortperm(heights)

    # Find the mass radius
    accu_mass = 0.0u"Msun"
    target_idx = 0
    for mass in masses[sort_idxs]
        accu_mass += mass
        accu_mass < mass_limit || break
        target_idx += 1
    end

    return heights[sort_idxs[target_idx]] * 2.0

end

"""
    computeMassQty(
        quantity::Vector{<:Number},
        masses::Vector{<:Unitful.Mass};
        <keyword arguments>
    )::Number

Compute the maximum value of `quantity` that "contains" `percet`% of the total mass.

# Arguments

  - `quantity::Vector{<:Number}`: Target quantity.
  - `masses::Vector{<:Unitful.Mass}`: Masses of the cells/particles.
  - `percent::Float64=90.0`: Target percentage of the total mass.

# Returns

  - The maximum value of `quantity` that "contains" `percet`% of the total mass.
"""
function computeMassQty(
    quantity::Vector{<:Number},
    masses::Vector{<:Unitful.Mass};
    percent::Float64=90.0,
)::Number

    (
        0 < percent <=100  ||
        throw(ArgumentError("computeMassQty: The argument `percent` must be between 0 and 100, \
        but I got $(percent)"))
    )

    # Check for missing data
    !any(isempty, [quantity, masses]) || return zero(eltype(quantity))

    # Compute the mass limit
    mass_limit = sum(masses) * (percent / 100.0)

    sort_idxs = sortperm(quantity)

    # Find the mass radius
    accu_mass = 0.0u"Msun"
    target_idx = 0
    for mass in masses[sort_idxs]
        accu_mass += mass
        accu_mass < mass_limit || break
        target_idx += 1
    end

    return quantity[sort_idxs[target_idx]]

end

"""
    computeMassPercent(
        quantity::Vector{<:Number},
        masses::Vector{<:Unitful.Mass},
        qty_limit::Number,
    )::Float64

Compute the fraction of the total mass "contained" within a given value of `quantity`.

# Arguments

  - `quantity::Vector{<:Number}`: Target quantity.
  - `masses::Vector{<:Unitful.Mass}`: Masses of the cells/particles.
  - `qty_limit::Number`: Limit value of the target quantity.

# Returns

  - The fraction of the total mass "contained" within a given value of `quantity`.
"""
function computeMassFraction(
    quantity::Vector{<:Number},
    masses::Vector{<:Unitful.Mass},
    qty_limit::Number,
)::Float64

    # Check for missing data
    !any(isempty, [quantity, masses]) || return 0.0

    qty_limit < maximum(quantity) || return 1.0

    # Find the indices of all the cells/particles with `quantity` < `qty_limit`
    idxs = map(x -> x <= qty_limit, quantity)

    return sum(masses) / sum(masses[idxs])

end

"""
    computeMetalMass(data_dict::Dict, type_symbol::Symbol)::Vector{<:Unitful.Mass}

Compute the total mass of metals (elements above helium) in each cell/particle.

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
  - `type_symbol::Symbol`: For which cell/particle type the metal mass will be calculated. The possibilities are `:stars` and `:gas`.

# Returns

  - The total metal mass in each cell/particle.
"""
function computeMetalMass(data_dict::Dict, type_symbol::Symbol)::Vector{<:Unitful.Mass}

    if type_symbol == :gas
        return setPositive(data_dict[:gas]["GZ  "]) .* data_dict[:gas]["MASS"]
    elseif type_symbol == :stars
        return setPositive(data_dict[:stars]["GZ2 "]) .* data_dict[:stars]["MASS"]
    else
        throw(ArgumentError("computeMetalMass: `type_symbol` can only be :stars or :gas, \
        but I got :$(type_symbol)"))
    end

end

"""
    computeElementMass(
        data_dict::Dict,
        type_symbol::Symbol,
        element::Symbol,
    )::Vector{<:Unitful.Mass}

Compute the total mass of `element` in each cell/particle.

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
  - `type_symbol::Symbol`: For which cell/particle type the element mass will be calculated. The possibilities are `:stars` and `:gas`.
  - `element::Symbol`: Target element. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).

# Returns

  - The total mass of `element` in each cell/particle.
"""
function computeElementMass(
    data_dict::Dict,
    type_symbol::Symbol,
    element::Symbol,
)::Vector{<:Unitful.Mass}

    (
        element ∈ keys(ELEMENT_INDEX) ||
        throw(ArgumentError("computeElementMass: :$(element) is not a tracked element, \
        the options are the keys of `ELEMENT_INDEX`, see `./src/constants.jl`"))
    )

    if type_symbol == :gas
        z_block = "GMET"
    elseif type_symbol == :stars
        z_block = "GME2"
    else
        throw(ArgumentError("computeElementMass: `type_symbol` can only be :stars or :gas, \
        but I got :$(type_symbol)"))
    end

    values = data_dict[type_symbol]

    if any(isempty, [values[z_block], values["MASS"]])
        return Unitful.Mass[]
    end

    return setPositive(values[z_block][ELEMENT_INDEX[element], :]) .* values["MASS"]

end

@doc raw"""
    computeGlobalAbundance(
        data_dict::Dict,
        type_symbol::Symbol,
        element::Symbol;
        <keyword arguments>
    )::Float64

Compute the total abundance of a given element, as $n_X / n_H$ where $n_X$ is the number of atoms of element $X$ and $n_H$ the number of hydrogen atoms.

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
  - `type_symbol::Symbol`: For which cell/particle type the abundance will be calculated. The possibilities are `:stars` and `:gas`.
  - `element::Symbol`: Target element. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
  - `solar::Bool=false`: If the result will be normalized to the solar abundance or not.

# Returns

  - The total abundance of `element`.
"""
function computeGlobalAbundance(
    data_dict::Dict,
    type_symbol::Symbol,
    element::Symbol;
    solar::Bool=false,
)::Float64

    metal_mass = sum(computeElementMass(data_dict, type_symbol, element); init=0.0u"Msun")
    hydrogen_mass = sum(computeElementMass(data_dict, type_symbol, :H); init=0.0u"Msun")

    (
        !iszero(hydrogen_mass) ||
        throw(ArgumentError("computeElementMass: I got 0 for the mass of hydrogen, \
        which should not be possible"))
    )

    # Compute the number of atoms
    n_X = metal_mass / ATOMIC_WEIGHTS[element]
    n_H = hydrogen_mass / ATOMIC_WEIGHTS[:H]

    # Compute the relative abundance of `element`
    abundance = ustrip(Unitful.NoUnits, n_X / n_H)

    return abundance / (solar ? exp10(SOLAR_ABUNDANCE[element] - 12.0) : 1.0)

end

"""
    computeStellarAge(data_dict::Dict)::Vector{<:Unitful.Time}

Compute the age of the stars.

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

  - The stellar ages.
"""
function computeStellarAge(data_dict::Dict)::Vector{<:Unitful.Time}

    birth_ticks = data_dict[:stars]["GAGE"]

    !isempty(birth_ticks) || return Unitful.Time[]

    if data_dict[:sim_data].cosmological
        # Go from scale factor to physical time
        birth_times = computeTime(birth_ticks, data_dict[:snap_data].header)
    else
        birth_times = birth_ticks
    end

    return data_dict[:snap_data].physical_time .- birth_times

end

"""
    computeIonizedMass(data_dict::Dict; <keyword arguments>)::Vector{<:Unitful.Mass}

Compute the ionized hydrogen mass of every gas cell in `data`.

The constant value [`HYDROGEN_MASSFRAC`](@ref) is used as the fraction of gas mass that is hydrogen.

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
  - `normalize::Bool=true`: If the output will be normalize to eliminate the stellar fraction of the gas cells. Only relevant for simulation with our routine, and for cells that have entered it at least once.

# Returns

  - The mass of ionized hydrogen in every gas cell.
"""
function computeIonizedMass(data_dict::Dict; normalize::Bool=true)::Vector{<:Unitful.Mass}

    dg = data_dict[:gas]

    !(isempty(dg["MASS"]) || isempty(dg["NHP "]) || isempty(dg["NH  "])) || return Unitful.Mass[]

    if "FRAC" ∈ keys(dg) && !isempty(dg["FRAC"])

        fi = Vector{Float64}(undef, length(dg["MASS"]))

        @inbounds for i in eachindex(fi)

            @inbounds if !isnan(dg["FRAC"][1, i]) && dg["DTIM"][i] < dg["TAUS"][i]

                # Fraction of ionized hydrogen according to our model
                @inbounds if normalize
                    fi[i] = dg["FRAC"][1, i] / (1.0 - dg["FRAC"][4, i])
                else
                    fi[i] = dg["FRAC"][1, i]
                end

            else

                fi[i] = dg["NHP "][i] / (dg["NHP "][i] + dg["NH  "][i])

            end

        end

        # # Fraction of ionized hydrogen according to our model
        # if normalize
        #     f_HII = dg["FRAC"][1, :] ./ (1.0 .- dg["FRAC"][4, :])
        # else
        #     f_HII = dg["FRAC"][1, :]
        # end

        # # Allocate memory
        # fi = copy(f_HII)

        # # When there is no data from our model, use the fraction of ionized hydrogen from Arepo
        # @inbounds for (i, (nh, nhp)) in enumerate(zip(dg["NH  "], dg["NHP "]))

        #     @inbounds if isnan(fi[i])
        #         fi[i] = nhp / (nhp + nh)
        #     end

        # end

    else

        # For simulations without our routine use the fraction of ionized hydrogen according to Arepo
        fi = dg["NHP "] ./ (dg["NHP "] .+ dg["NH  "])

    end

    return fi .* dg["MASS"] .* HYDROGEN_MASSFRAC

end

"""
    computeAtomicMass(data_dict::Dict; <keyword arguments>)::Vector{<:Unitful.Mass}

Compute the atomic hydrogen mass of every gas cell in `data`.

The constant value [`HYDROGEN_MASSFRAC`](@ref) is used as the fraction of gas mass that is hydrogen.

For simulations without our routine use the pressure relation in Blitz et al. (2006) to separate atomic from molecular gas in the neutral phase given by Arepo.

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
  - `normalize::Bool=true`: If the output will be normalize to eliminate the stellar fraction of the gas cells. Only relevant for simulation with our routine, and for cells that have entered it at least once.

# Returns

  - The mass of atomic hydrogen in every gas cell.

# References

L. Blitz et al. (2006). *The Role of Pressure in GMC Formation II: The H2-Pressure Relation*. The Astrophysical Journal, **650(2)**, 933. [doi:10.1086/505417](https://doi.org/10.1086/505417)
"""
function computeAtomicMass(data_dict::Dict; normalize::Bool=true)::Vector{<:Unitful.Mass}

    dg = data_dict[:gas]

    !(isempty(dg["MASS"]) || isempty(dg["NHP "]) || isempty(dg["NH  "])) || return Unitful.Mass[]

    if "FRAC" ∈ keys(dg) && !isempty(dg["FRAC"])

        fa = Vector{Float64}(undef, length(dg["MASS"]))

        @inbounds for i in eachindex(fa)

            @inbounds if !isnan(dg["FRAC"][2, i]) && dg["DTIM"][i] < dg["TAUS"][i]

                # Fraction of atomic hydrogen according to our model
                @inbounds if normalize
                    fa[i] = dg["FRAC"][2, i] / (1.0 - dg["FRAC"][4, i])
                else
                    fa[i] = dg["FRAC"][2, i]
                end

            else

                fa[i] = dg["NH  "][i] / (dg["NHP "][i] + dg["NH  "][i])

            end

        end

        # # Fraction of atomic hydrogen according to our model
        # if normalize
        #     f_HI = dg["FRAC"][2, :] ./ (1.0 .- dg["FRAC"][4, :])
        # else
        #     f_HI = dg["FRAC"][2, :]
        # end

        # # Allocate memory
        # fa = copy(f_HI)

        # # When there is no data from our model, use the fraction of neutral hydrogen from Arepo
        # # assuming that the fraction of molecular hydrogen is 0
        # @inbounds for (i, (nh, nhp)) in enumerate(zip(dg["NH  "], dg["NHP "]))

        #     @inbounds if isnan(fa[i])
        #         fa[i] = nh / (nhp + nh)
        #     end

        # end

    elseif !isempty(dg["PRES"])

        # Fraction of neutral hydrogen according to Arepo
        fn = dg["NH  "] ./ (dg["NHP "] .+ dg["NH  "])

        relative_pressure = uconvert.(Unitful.NoUnits, dg["PRES"] ./ P0).^ALPHA_BLITZ

        # Fraction of molecular hydrogen according to the pressure relation in Blitz et al. (2006)
        fm = 1.0 ./ (1.0 .+ relative_pressure)

        # Use the fraction of neutral hydrogen that is not molecular according to the pressure relation,
        # unless that value is negative, in which case assume taht all neutral hydrogen is molecular
        fa = setPositive(fn .- fm)

    else

        return Unitful.Mass[]

    end

    return fa .* dg["MASS"] .* HYDROGEN_MASSFRAC

end

"""
    computeMolecularMass(data_dict::Dict; <keyword arguments>)::Vector{<:Unitful.Mass}

Compute the molecular hydrogen mass of every gas cell in `data`.

The constant value [`HYDROGEN_MASSFRAC`](@ref) is used as the fraction of gas mass that is hydrogen.

For simulations without our routine use the pressure relation in Blitz et al. (2006) to separate molecular from atomic gas in the neutral phase given by Arepo.

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
  - `normalize::Bool=true`: If the output will be normalize to eliminate the stellar fraction of the gas cells. Only relevant for simulation with our routine, and for cells that have entered it at least once.

# Returns

  - The mass of molecular hydrogen in every gas cell.

# References

L. Blitz et al. (2006). *The Role of Pressure in GMC Formation II: The H2-Pressure Relation*. The Astrophysical Journal, **650(2)**, 933. [doi:10.1086/505417](https://doi.org/10.1086/505417)
"""
function computeMolecularMass(data_dict::Dict; normalize::Bool=true)::Vector{<:Unitful.Mass}

    dg = data_dict[:gas]

    !isempty(dg["MASS"]) || return Unitful.Mass[]

    if "FRAC" ∈ keys(dg) && !isempty(dg["FRAC"])

        fm = Vector{Float64}(undef, length(dg["MASS"]))

        @inbounds for i in eachindex(fm)

            @inbounds if !isnan(dg["FRAC"][3, i]) && dg["DTIM"][i] < dg["TAUS"][i]

                # Fraction of molecular hydrogen according to our model
                @inbounds if normalize
                    fm[i] = dg["FRAC"][3, i] / (1.0 - dg["FRAC"][4, i])
                else
                    fm[i] = dg["FRAC"][3, i]
                end

            else

                fm[i] = 0.0

            end

        end

        # # Fraction of molecular hydrogen according to our model
        # if normalize
        #     f_H2 = dg["FRAC"][3, :] ./ (1.0 .- dg["FRAC"][4, :])
        # else
        #     f_H2 = dg["FRAC"][3, :]
        # end

        # # When there is no data from our model, asume 0 molecular hydrogen
        # fm = replace!(f_H2, NaN => 0.0)

    elseif !isempty(dg["PRES"]) && !isempty(dg["NHP "]) && !isempty(dg["NH  "])

        # Fraction of neutral hydrogen according to Arepo
        fn = dg["NH  "] ./ (dg["NHP "] .+ dg["NH  "])

        relative_pressure = uconvert.(Unitful.NoUnits, dg["PRES"] ./ P0).^ALPHA_BLITZ

        # Fraction of molecular hydrogen according to the pressure relation in Blitz et al. (2006)
        fp = 1.0 ./ (1.0 .+ relative_pressure)

        # Use the fraction of molecular hydrogen according to the pressure relation, unless
        # that value is larger than the fraction of neutral hydrogen according to Arepo,
        # in which case assume taht all neutral hydrogen is molecular
        fm = [n >= p ? p : n for (n, p) in zip(fn, fp)]

    else

        return Unitful.Mass[]

    end

    return fm .* dg["MASS"] .* HYDROGEN_MASSFRAC

end

"""
    computeNeutralMass(data_dict::Dict; <keyword arguments>)::Vector{<:Unitful.Mass}

Compute the neutral hydrogen mass of every gas cell in `data`.

The constant value [`HYDROGEN_MASSFRAC`](@ref) is used as the fraction of gas mass that is hydrogen.

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
  - `normalize::Bool=true`: If the output will be normalize to eliminate the stellar fraction of the gas cells. Only relevant for simulation with our routine, and for cells that have entered it at least once.

# Returns

  - The mass of neutral hydrogen in every gas cell.
"""
function computeNeutralMass(data_dict::Dict; normalize::Bool=true)::Vector{<:Unitful.Mass}

    dg = data_dict[:gas]

    !(isempty(dg["MASS"]) || isempty(dg["NHP "]) || isempty(dg["NH  "])) || return Unitful.Mass[]

    if "FRAC" ∈ keys(dg) && !isempty(dg["FRAC"])

        fa = Vector{Float64}(undef, length(dg["MASS"]))
        fm = Vector{Float64}(undef, length(dg["MASS"]))

        @inbounds for i in eachindex(fa)

            @inbounds if !isnan(dg["FRAC"][2, i]) && dg["DTIM"][i] < dg["TAUS"][i]

                # Fraction of atomic and molecular hydrogen according to our model
                @inbounds if normalize
                    fa[i] = dg["FRAC"][2, i] / (1.0 - dg["FRAC"][4, i])
                    fm[i] = dg["FRAC"][3, i] / (1.0 - dg["FRAC"][4, i])
                else
                    fa[i] = dg["FRAC"][2, i]
                    fm[i] = dg["FRAC"][3, i]
                end

            else

                # When there is no data from our model, use the fraction of neutral hydrogen from Arepo
                # assuming that the fraction of molecular hydrogen is 0
                fa[i] = dg["NH  "][i] / (dg["NHP "][i] + dg["NH  "][i])
                fm[i] = 0.0

            end

        end

        fn = fa .+ fm

        # # Fraction of atomic and molecular hydrogen according to our model
        # if normalize
        #     f_HI = dg["FRAC"][2, :] ./ (1.0 .- dg["FRAC"][4, :])
        #     f_H2 = dg["FRAC"][3, :] ./ (1.0 .- dg["FRAC"][4, :])
        # else
        #     f_HI = dg["FRAC"][2, :]
        #     f_H2 = dg["FRAC"][3, :]
        # end

        # # Allocate memory
        # fa = copy(f_HI)
        # # When there is no data from our model, use the fraction of neutral hydrogen from Arepo
        # # assuming that the fraction of molecular hydrogen is 0
        # @inbounds for (i, (nh, nhp)) in enumerate(zip(dg["NH  "], dg["NHP "]))

        #     @inbounds if isnan(fhii)
        #         fa[i] = nh / (nhp + nh)
        #     end

        # end

        # # When there is no data from our model, asume 0 molecular hydrogen
        # fm = replace!(f_H2, NaN => 0.0)

        # fn = fa .+ fm

    else

        # Fraction of neutral hydrogen according to Arepo
        fn = @. dg["NH  "] / (dg["NHP "] + dg["NH  "])

    end

    return fn .* dg["MASS"] .* HYDROGEN_MASSFRAC

end

"""
    computeStellarGasMass(data_dict::Dict)::Vector{<:Unitful.Mass}

Compute the "stellar mass" of every gas cell in `data`, which will be other than 0 only for simulation with our routine, and for cells that have entered it at least once.

The constant value [`HYDROGEN_MASSFRAC`](@ref) is used as the fraction of gas mass that is hydrogen. This is applied only for consistency with the other mass rutines ([`computeIonizedMass`](@ref), [`computeAtomicMass`](@ref), [`computeMolecularMass`](@ref), and [`computeNeutralMass`](@ref)). Notice that there is no physical meaning to the "stellar mass" of a gas cell as used in our model. So, it makes no sense to question if the "stellar fraction" computed here is a fraction of the total gas mass or only of the hydrogen mass.

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

  - The "stellar mass" of every gas cell.
"""
function computeStellarGasMass(data_dict::Dict)::Vector{<:Unitful.Mass}

    dg = data_dict[:gas]

    !isempty(dg["MASS"]) || return Unitful.Mass[]

    if "FRAC" ∈ keys(dg) && !isempty(dg["FRAC"])

        # When there is no data from our model, use a stellar fraction of 0
        fm = replace!(dg["FRAC"][4, :], NaN => 0.0)

    else

        return zeros(typeof(1.0u"Msun"), length(dg["MASS"]))

    end

    return fm .* dg["MASS"] .* HYDROGEN_MASSFRAC

end

"""
    computeSFR(
        data_dict::Dict;
        <keyword arguments>
    )::Vector{<:Unitful.MassFlow}

Compute the star formation rate of each stellar particle in `data`.

For stellar particles younger that `age_resol`, the SFR is its mass divided by `age_resol`. It is defined as 0 for older particles.

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
  - `age_resol::Unitful.Time=AGE_RESOLUTION`: Age resolution for the SFR.

# Returns

  - The star formation rate of each stellar particle.
"""
function computeSFR(
    data_dict::Dict;
    age_resol::Unitful.Time=AGE_RESOLUTION,
)::Vector{<:Unitful.MassFlow}

    # Compute the stellar ages
    ages = computeStellarAge(data_dict)

    !isempty(ages) || return Unitful.MassFlow[]

    # Allocate memory
    sfr = zeros(typeof(1.0u"Msun*yr^-1"), length(ages))

    # Find the stellar particles younger than `age_resol`
    idxs = map(x -> x <= age_resol, ages)

    # Compute the SFR
    sfr[idxs] .= data_dict[:stars]["MASS"][idxs] ./ age_resol

    return sfr

end

"""
    integrateQty(data_dict::Dict, quantity::Symbol)::Number

Compute an integrated quantity for the whole system in `data_dict`.

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
      + `:stellar_area_density`   -> Stellar area mass density, for a radius of `FILTER_R`.
      + `:gas_area_density`       -> Gas area mass density, for a radius of `FILTER_R`.
      + `:molecular_area_density` -> Molecular hydrogen area mass density, for a radius of `FILTER_R`.
      + `:atomic_area_density`    -> Atomic hydrogen area mass density, for a radius of `FILTER_R`.
      + `:ionized_area_density`   -> Ionized hydrogen area mass density, for a radius of `FILTER_R`.
      + `:neutral_area_density`   -> Neutral hydrogen area mass density, for a radius of `FILTER_R`.
      + `:sfr_area_density`       -> Star formation rate area density, for the last `AGE_RESOLUTION_ρ` and a radius of `FILTER_R`.
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

# Returns

  - The velue of `quantity` for the whole system in `data_dict`.
"""
function integrateQty(data_dict::Dict, quantity::Symbol)::Number

    if quantity == :stellar_mass

        integrated_qty = sum(data_dict[:stars]["MASS"]; init=0.0u"Msun")

    elseif quantity == :gas_mass

        integrated_qty = sum(data_dict[:gas]["MASS"]; init=0.0u"Msun")

    elseif quantity == :hydrogen_mass

        integrated_qty = sum(data_dict[:gas]["MASS"]; init=0.0u"Msun") * HYDROGEN_MASSFRAC

    elseif quantity == :dm_mass

        integrated_qty = sum(data_dict[:halo]["MASS"]; init=0.0u"Msun")

    elseif quantity == :bh_mass

        integrated_qty = sum(data_dict[:black_hole]["MASS"]; init=0.0u"Msun")

    elseif quantity == :molecular_mass

        integrated_qty = sum(computeMolecularMass(data_dict); init=0.0u"Msun")

    elseif quantity == :atomic_mass

        integrated_qty = sum(computeAtomicMass(data_dict); init=0.0u"Msun")

    elseif quantity == :ionized_mass

        integrated_qty = sum(computeIonizedMass(data_dict); init=0.0u"Msun")

    elseif quantity == :neutral_mass

        integrated_qty = sum(computeNeutralMass(data_dict); init=0.0u"Msun")

    elseif quantity == :stellar_number

        integrated_qty = length(data_dict[:stars]["MASS"])

    elseif quantity == :gas_number

        integrated_qty = length(data_dict[:gas]["MASS"])

    elseif quantity == :dm_number

        integrated_qty = length(data_dict[:halo]["MASS"])

    elseif quantity == :bh_number

        integrated_qty = length(data_dict[:black_hole]["MASS"])

    elseif quantity == :molecular_fraction

        molecular_mass = sum(computeMolecularMass(data_dict); init=0.0u"Msun")
        hydrogen_mass = sum(data_dict[:gas]["MASS"]; init=0.0u"Msun") * HYDROGEN_MASSFRAC

        if iszero(hydrogen_mass)
            integrated_qty = NaN
        else
            integrated_qty = molecular_mass / hydrogen_mass
        end

    elseif quantity == :atomic_fraction

        atomic_mass = sum(computeAtomicMass(data_dict); init=0.0u"Msun")
        hydrogen_mass = sum(data_dict[:gas]["MASS"]; init=0.0u"Msun") * HYDROGEN_MASSFRAC

        if iszero(hydrogen_mass)
            integrated_qty = NaN
        else
            integrated_qty = atomic_mass / hydrogen_mass
        end

    elseif quantity == :ionized_fraction

        ionized_mass = sum(computeIonizedMass(data_dict); init=0.0u"Msun")
        hydrogen_mass = sum(data_dict[:gas]["MASS"]; init=0.0u"Msun") * HYDROGEN_MASSFRAC

        if iszero(hydrogen_mass)
            integrated_qty = NaN
        else
            integrated_qty = ionized_mass / hydrogen_mass
        end

    elseif quantity == :neutral_fraction

        neutral_mass = sum(computeNeutralMass(data_dict); init=0.0u"Msun")
        hydrogen_mass = sum(data_dict[:gas]["MASS"]; init=0.0u"Msun") * HYDROGEN_MASSFRAC

        if iszero(hydrogen_mass)
            integrated_qty = NaN
        else
            integrated_qty = neutral_mass / hydrogen_mass
        end

    elseif quantity == :stellar_area_density

        integrated_qty = sum(data_dict[:stars]["MASS"]; init=0.0u"Msun") / area(FILTER_R)

    elseif quantity == :gas_area_density

        integrated_qty = sum(data_dict[:gas]["MASS"]; init=0.0u"Msun") / area(FILTER_R)

    elseif quantity == :molecular_area_density

        integrated_qty = sum(computeMolecularMass(data_dict); init=0.0u"Msun") / area(FILTER_R)

    elseif quantity == :atomic_area_density

        integrated_qty = sum(computeAtomicMass(data_dict); init=0.0u"Msun") / area(FILTER_R)

    elseif quantity == :ionized_area_density

        integrated_qty = sum(computeIonizedMass(data_dict); init=0.0u"Msun") / area(FILTER_R)

    elseif quantity == :neutral_area_density

        integrated_qty = sum(computeNeutralMass(data_dict); init=0.0u"Msun") / area(FILTER_R)

    elseif quantity == :sfr_area_density

        sfr = sum(computeSFR(data_dict; age_resol=AGE_RESOLUTION_ρ); init=0.0u"Msun*yr^-1")

        integrated_qty = sfr / area(FILTER_R)

    elseif quantity == :gas_metallicity

        metal_mass = sum(computeMetalMass(data_dict, :gas); init=0.0u"Msun")
        gas_mass = sum(data_dict[:gas]["MASS"]; init=0.0u"Msun")

        if iszero(gas_mass)
            integrated_qty = NaN
        else
            integrated_qty = (metal_mass / gas_mass) / SOLAR_METALLICITY
        end

    elseif quantity == :stellar_metallicity

        metal_mass = sum(computeMetalMass(data_dict, :stars); init=0.0u"Msun")
        stellar_mass = sum(data_dict[:stars]["MASS"]; init=0.0u"Msun")

        if iszero(stellar_mass)
            integrated_qty = NaN
        else
            integrated_qty = (metal_mass / stellar_mass) / SOLAR_METALLICITY
        end

    elseif quantity ∈ GAS_ABUNDANCE

        element_symbol = Symbol(first(split(string(quantity), "_")))

        abundance = 12 + log10(computeGlobalAbundance(data_dict, :gas, element_symbol))
        integrated_qty = isinf(abundance) ? NaN : abundance

    elseif quantity ∈ STELLAR_ABUNDANCE

        element_symbol = Symbol(first(split(string(quantity), "_")))

        abundance = 12 + log10(computeGlobalAbundance(data_dict, :stars, element_symbol))
        integrated_qty = isinf(abundance) ? NaN : abundance

    elseif quantity == :stellar_specific_am

        positions = data_dict[:stars]["POS "]
        velocities = data_dict[:stars]["VEL "]
        masses = data_dict[:stars]["MASS"]

        if any(isempty, [positions, velocities, masses])
            integrated_qty = NaN
        else
            J = norm(computeTotalAngularMomentum(positions, velocities, masses; normal=false))
            integrated_qty = J / sum(masses)
        end

    elseif quantity == :gas_specific_am

        positions = data_dict[:gas]["POS "]
        velocities = data_dict[:gas]["VEL "]
        masses = data_dict[:gas]["MASS"]

        if any(isempty, [positions, velocities, masses])
            integrated_qty = NaN
        else
            J = norm(computeTotalAngularMomentum(positions, velocities, masses; normal=false))
            integrated_qty = J / sum(masses)
        end

    elseif quantity == :dm_specific_am

        positions = data_dict[:halo]["POS "]
        velocities = data_dict[:halo]["VEL "]
        masses = data_dict
        masses = [:halo]["MASS"]

        if any(isempty, [positions, velocities, masses])
            integrated_qty = NaN
        else
            J = norm(computeTotalAngularMomentum(positions, velocities, masses; normal=false))
            integrated_qty = J / sum(masses)
        end

    elseif quantity == :sfr

        # Get the global index (index in the context of the whole simulation) of the current snapshot
        present_idx = data_dict[:snap_data].global_index

        if present_idx == 1

            integrated_qty = 0.0u"Msun*yr^-1"

        else

            # Get the physical times
            times = data_dict[:sim_data].table[:, 5]
            # Compute the time between snapshots
            Δt = times[present_idx] - times[present_idx - 1]

            integrated_qty = sum(computeSFR(data_dict; age_resol=Δt); init=0.0u"Msun*yr^-1")

        end

    elseif quantity == :ssfr

        # Get the global index (index in the context of the whole simulation) of the current snapshot
        present_idx = data_dict[:snap_data].global_index

        # Compute the total stellar mass
        stellar_mass = sum(data_dict[:stars]["MASS"]; init=0.0u"Msun")

        if present_idx == 1 || iszero(stellar_mass)

            integrated_qty = 0.0u"yr^-1"

        else

            # Get the physical times
            times = data_dict[:sim_data].table[:, 5]
            # Compute the time between snapshots
            Δt = times[present_idx] - times[present_idx - 1]

            integrated_qty = sum(
                computeSFR(data_dict; age_resol=Δt);
                init=0.0u"Msun*yr^-1",
            ) / stellar_mass

        end

    elseif quantity == :observational_sfr

        integrated_qty = sum(computeSFR(data_dict; age_resol=AGE_RESOLUTION); init=0.0u"Msun*yr^-1")

    elseif quantity == :observational_ssfr

        sfr = sum(computeSFR(data_dict; age_resol=AGE_RESOLUTION); init=0.0u"Msun*yr^-1")
        stellar_mass = sum(data_dict[:stars]["MASS"]; init=0.0u"Msun")

        if iszero(stellar_mass)
            integrated_qty = 0.0u"yr^-1"
        else
            integrated_qty = sfr / stellar_mass
        end

    elseif quantity == :scale_factor

        integrated_qty = data_dict[:sim_data].table[data_dict[:snap_data].global_index, 3]

    elseif quantity == :redshift

        integrated_qty = data_dict[:sim_data].table[data_dict[:snap_data].global_index, 4]

    elseif quantity == :physical_time

        integrated_qty = data_dict[:sim_data].table[data_dict[:snap_data].global_index, 5]

    elseif quantity == :lookback_time

        integrated_qty = data_dict[:sim_data].table[data_dict[:snap_data].global_index, 6]

    else

        throw(ArgumentError("integrateQty: I don't recognize the quantity :$(quantity)"))

    end

    return integrated_qty

end

"""
    scatterQty(data_dict::Dict, quantity::Symbol)::Vector{<:Number}

Compute a quantity for each cell/particle in `data_dict`.

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

# Returns

  - The values of `quantity` for every cell/particle.
"""
function scatterQty(data_dict::Dict, quantity::Symbol)::Vector{<:Number}

    if quantity == :stellar_mass

        scatter_qty = data_dict[:stars]["MASS"]

    elseif quantity == :gas_mass

        scatter_qty = data_dict[:gas]["MASS"]

    elseif quantity == :hydrogen_mass

        scatter_qty = data_dict[:gas]["MASS"] .* HYDROGEN_MASSFRAC

    elseif quantity == :dm_mass

        scatter_qty = data_dict[:halo]["MASS"]

    elseif quantity == :bh_mass

        scatter_qty = data_dict[:black_hole]["MASS"]

    elseif quantity == :molecular_mass

        scatter_qty = computeMolecularMass(data_dict)

    elseif quantity == :atomic_mass

        scatter_qty = computeAtomicMass(data_dict)

    elseif quantity == :ionized_mass

        scatter_qty = computeIonizedMass(data_dict)

    elseif quantity == :neutral_mass

        scatter_qty = computeNeutralMass(data_dict)

    elseif quantity == :molecular_fraction

        molecular_mass = computeMolecularMass(data_dict)
        hydrogen_mass = data_dict[:gas]["MASS"] .* HYDROGEN_MASSFRAC

        scatter_qty = molecular_mass ./ hydrogen_mass

    elseif quantity == :atomic_fraction

        atomic_mass = computeAtomicMass(data_dict)
        hydrogen_mass = data_dict[:gas]["MASS"] .* HYDROGEN_MASSFRAC

        scatter_qty = atomic_mass ./ hydrogen_mass

    elseif quantity == :ionized_fraction

        ionized_mass = computeIonizedMass(data_dict)
        hydrogen_mass = data_dict[:gas]["MASS"] .* HYDROGEN_MASSFRAC

        scatter_qty = ionized_mass ./ hydrogen_mass

    elseif quantity == :neutral_fraction

        neutral_mass = computeNeutralMass(data_dict)
        hydrogen_mass = data_dict[:gas]["MASS"] .* HYDROGEN_MASSFRAC

        scatter_qty = neutral_mass ./ hydrogen_mass

    elseif quantity == :molecular_neutral_fraction

        molecular_mass = computeMolecularMass(data_dict)
        atomic_mass    = computeAtomicMass(data_dict)

        scatter_qty = molecular_mass ./ (atomic_mass .+ molecular_mass)

    elseif quantity == :gas_mass_density

        scatter_qty = data_dict[:gas]["RHO "]

    elseif quantity == :hydrogen_mass_density

        scatter_qty = data_dict[:gas]["RHO "] .* HYDROGEN_MASSFRAC

    elseif quantity == :gas_number_density

        scatter_qty = data_dict[:gas]["RHO "] ./ Unitful.mp

    elseif quantity == :molecular_number_density

        molecular_mass = computeMolecularMass(data_dict)
        volumes = data_dict[:gas]["MASS"] ./ data_dict[:gas]["RHO "]

        scatter_qty = (molecular_mass ./ volumes) ./ (2 * Unitful.mp)

    elseif quantity == :atomic_number_density

        atomic_mass = computeAtomicMass(data_dict)
        volumes = data_dict[:gas]["MASS"] ./ data_dict[:gas]["RHO "]

        scatter_qty = (atomic_mass ./ volumes) ./ Unitful.mp

    elseif quantity == :ionized_number_density

        ionized_mass = computeIonizedMass(data_dict)
        volumes = data_dict[:gas]["MASS"] ./ data_dict[:gas]["RHO "]

        scatter_qty = (ionized_mass ./ volumes) ./ Unitful.mp

    elseif quantity == :neutral_number_density

        neutral_mass = computeNeutralMass(data_dict)
        volumes = data_dict[:gas]["MASS"] ./ data_dict[:gas]["RHO "]

        scatter_qty = (neutral_mass ./ volumes) ./ Unitful.mp

    elseif quantity == :gas_metallicity

        scatter_qty = setPositive(data_dict[:gas]["GZ  "]) ./ SOLAR_METALLICITY

    elseif quantity == :stellar_metallicity

        scatter_qty = setPositive(data_dict[:stars]["GZ2 "]) ./ SOLAR_METALLICITY

    elseif quantity ∈ GAS_ABUNDANCE

        element_symbol = Symbol(first(split(string(quantity), "_")))

        element_mass = computeElementMass(data_dict, :gas, element_symbol)
        hydrogen_mass = computeElementMass(data_dict, :gas, :H)

        if isempty(hydrogen_mass)
            scatter_qty = Float64[]
        else
            n_X = element_mass ./ ATOMIC_WEIGHTS[element_symbol]
            n_H = hydrogen_mass ./ ATOMIC_WEIGHTS[:H]

            abundance = ustrip.(Unitful.NoUnits, n_X ./ n_H)

            scatter_qty = 12 .+ log10.(abundance)
            replace!(x -> isinf(x) ? NaN : x, scatter_qty)
        end

    elseif quantity ∈ STELLAR_ABUNDANCE

        element_symbol = Symbol(first(split(string(quantity), "_")))

        element_mass = computeElementMass(data_dict, :stars, element_symbol)
        hydrogen_mass = computeElementMass(data_dict, :stars, :H)

        if isempty(hydrogen_mass)
            scatter_qty = Float64[]
        else
            n_X = element_mass ./ ATOMIC_WEIGHTS[element_symbol]
            n_H = hydrogen_mass ./ ATOMIC_WEIGHTS[:H]

            abundance = ustrip.(Unitful.NoUnits, n_X ./ n_H)

            scatter_qty = 12 .+ log10.(abundance)
            replace!(x -> isinf(x) ? NaN : x, scatter_qty)
        end

    elseif quantity == :stellar_radial_distance

        scatter_qty = computeDistance(data_dict[:stars]["POS "])

    elseif quantity == :gas_radial_distance

        scatter_qty = computeDistance(data_dict[:gas]["POS "])

    elseif quantity == :dm_radial_distance

        scatter_qty = computeDistance(data_dict[:halo]["POS "])

    elseif quantity == :stellar_xy_distance

        if isempty(data_dict[:stars]["POS "])
            scatter_qty = eltype(data_dict[:stars]["POS "])[]
        else
            scatter_qty = computeDistance(data_dict[:stars]["POS "][1:2, :])
        end

    elseif quantity == :gas_xy_distance

        if isempty(data_dict[:gas]["POS "])
            scatter_qty = eltype(data_dict[:gas]["POS "])[]
        else
            scatter_qty = computeDistance(data_dict[:gas]["POS "][1:2, :])
        end

    elseif quantity == :dm_xy_distance

        if isempty(data_dict[:halo]["POS "])
            scatter_qty = eltype(data_dict[:halo]["POS "])[]
        else
            scatter_qty = computeDistance(data_dict[:halo]["POS "][1:2, :])
        end

    elseif quantity == :stellar_circularity

        @debug("scatterQty: The stellar circularity depends on the positions and velocities of all \
        cell/particles. So, after filtering, the result for a given star will change.")

        scatter_qty = computeStellarCircularity(data_dict)

    elseif quantity == :stellar_vcirc

        @debug("scatterQty: The stellar circular velocity depends on the positions and velocities \
        of all cell/particles. So, after filtering, the result for a given star will change.")

        _, scatter_qty = computeStellarVcirc(data_dict)

    elseif quantity == :stellar_vradial

        scatter_qty = computeStellarVpolar(data_dict, :radial)

    elseif quantity == :stellar_vtangential

        scatter_qty = computeStellarVpolar(data_dict, :tangential)

    elseif quantity == :stellar_vzstar

        scatter_qty = computeStellarVpolar(data_dict, :zstar)

    elseif quantity == :stellar_age

        scatter_qty = computeStellarAge(data_dict)

    elseif quantity == :sfr

        # Get the global index (index in the context of the whole simulation) of the current snapshot
        present_idx = data_dict[:snap_data].global_index

        if present_idx == 1

            scatter_qty = zeros(typeof(1.0u"Msun*yr^-1"), length(data_dict[:stars]["MASS"]))

        else

            # Get the physical times
            times = data_dict[:sim_data].table[:, 5]
            # Compute the time between snapshots
            Δt = times[present_idx] - times[present_idx - 1]

            scatter_qty = computeSFR(data_dict; age_resol=Δt)

        end

    elseif quantity == :ssfr

        # Get the global index (index in the context of the whole simulation) of the current snapshot
        present_idx = data_dict[:snap_data].global_index

        # Load the stellar masses
        stellar_masses = data_dict[:stars]["MASS"]

        if present_idx == 1 || iszero(stellar_mass)

            scatter_qty = zeros(typeof(1.0u"yr^-1"), length(stellar_masses))

        else

            # Get the physical times
            times = data_dict[:sim_data].table[:, 5]
            # Compute the time between snapshots
            Δt = times[present_idx] - times[present_idx - 1]

            scatter_qty = computeSFR(data_dict; age_resol=Δt) ./ stellar_masses

        end

    elseif quantity == :observational_sfr

        scatter_qty = computeSFR(data_dict; age_resol=AGE_RESOLUTION)

    elseif quantity == :observational_ssfr

        sfr = computeSFR(data_dict; age_resol=AGE_RESOLUTION)
        stellar_masses = data_dict[:stars]["MASS"]

        scatter_qty = sfr ./ stellar_masses

    elseif quantity == :temperature

        scatter_qty = log10.(ustrip.(u"K", data_dict[:gas]["TEMP"]))
        replace!(x -> isinf(x) ? NaN : x, scatter_qty)

    else

        throw(ArgumentError("scatterQty: I don't recognize the quantity :$(quantity)"))

    end

    return scatter_qty

end

"""
    computeVirialAccretion(
        present_dd::Dict,
        past_dd::Dict;
        <keyword arguments>
    )::NTuple{3,Unitful.Mass}

Compute the inflow, outflow, and net gain of mass for a given halo virial radius (``R_{200}``), between two snapshots.

# Arguments

  - `present_dd::Dict`: A dictionary, for the present snapshot, with the following shape:

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
  - `past_dd::Dict`: A dictionary, for the past snapshot, with the following shape:

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
  - `halo_idx::Int=1`: Index of the target halo (FoF group). Starts at 1.

# Returns

  - A tuple with three elements:

      + The net increase in mass.
      + The inflow mass.
      + The outflow mass.
"""
function computeVirialAccretion(
    present_dd::Dict,
    past_dd::Dict;
    halo_idx::Int=1,
)::NTuple{3,Unitful.Mass}

    # Find the tracers inside R200 in the present snapshot
    present_tracer_ids = tracersWithinR200(present_dd; halo_idx)
    # Find the tracers inside R200 in the past snapshot
    past_tracer_ids    = tracersWithinR200(past_dd; halo_idx)

    # Find the tracers that are inside R200 now, but where outside R200 in the past
    inflow_ids  = setdiff(present_tracer_ids, past_tracer_ids)
    # Find the tracers that were inside R200 in the past, but are now outside R200
    outflow_ids = setdiff(past_tracer_ids, present_tracer_ids)

    # Compute the mass of each tracer in physical units
    tracer_mass = TRACER_MASS * internalUnits("MASS", present_dd[:snap_data].path)

    # Compute the inflow mass
    inflow_mass = length(inflow_ids) * tracer_mass
    # Compute the outflow mass
    ouflow_mass = length(outflow_ids) * tracer_mass

    # Compute the net mass
    net_mass_increase = inflow_mass - ouflow_mass

    return net_mass_increase, inflow_mass, ouflow_mass

end

"""
    computeDiscAccretion(
        present_dd::Dict,
        past_dd::Dict;
        <keyword arguments>
    )::NTuple{3,Unitful.Mass}

Compute the inflow, outflow, and net gain of mass for a given cylinder, between two snapshots.

# Arguments

  - `present_dd::Dict`: A dictionary, for the present snapshot, with the following shape:

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
  - `past_dd::Dict`: A dictionary, for the past snapshot, with the following shape:

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
  - `max_r::Unitful.Length=FILTER_R`: Radius of the cylinder.
  - `max_z::Unitful.Length=5.0u"kpc"`: Half height of the cylinder.

# Returns

  - A tuple with three elements:

      + The net increase in mass.
      + The inflow mass.
      + The outflow mass.
"""
function computeDiscAccretion(
    present_dd::Dict,
    past_dd::Dict;
    max_r::Unitful.Length=FILTER_R,
    max_z::Unitful.Length=5.0u"kpc",
)::NTuple{3,Unitful.Mass}

    # Find the tracers inside a given cylinder in the present snapshot
    present_tracer_ids = tracersWithinDisc(present_dd; max_r, max_z)
    # Find the tracers inside a given cylinder in the past snapshot
    past_tracer_ids    = tracersWithinDisc(past_dd; max_r, max_z)

    # Find the tracers that are inside a given cylinder now, but where outside in the past
    inflow_ids  = setdiff(present_tracer_ids, past_tracer_ids)
    # Find the tracers that were inside a given cylinder in the past, but are now outside
    outflow_ids = setdiff(past_tracer_ids, present_tracer_ids)

    # Compute the mass of each tracer in physical units
    tracer_mass = TRACER_MASS * internalUnits("MASS", present_dd[:snap_data].path)

    # Compute the inflow mass
    inflow_mass = length(inflow_ids) * tracer_mass
    # Compute the outflow mass
    ouflow_mass = length(outflow_ids) * tracer_mass

    # Compute the net mass
    net_mass_increase = inflow_mass - ouflow_mass

    return net_mass_increase, inflow_mass, ouflow_mass

end

"""
    findRealStars(path::String)::Vector{Bool}

Find which stellar particles are real stars and not wind particles.

# Arguments

  - `path::String`: Path to the snapshot file or folder.

# Returns

  - A boolean vector with true for stars and false for wind particles.
"""
function findRealStars(path::String)::Vector{Bool}

    if isfile(path)

        (
            HDF5.ishdf5(path) ||
            throw(ArgumentError("findRealStars: The file $(path) is not in the \
            HDF5 format, I don't know how to read it"))
        )

        time_of_birth = h5open(path, "r") do snapshot
            if PARTICLE_CODE_NAME[:stars] ∉ keys(snapshot)
                Float64[]
            else
                read(snapshot[PARTICLE_CODE_NAME[:stars]], QUANTITIES["GAGE"].hdf5_name)
            end
        end

        return isempty(time_of_birth) ? Bool[] : map(isPositive, time_of_birth)

    elseif isdir(path)

        sub_files = glob("$(SNAP_BASENAME)_*.*.hdf5", path)

        (
            !isempty(sub_files) && all(HDF5.ishdf5, sub_files) ||
            throw(ArgumentError("findRealStars: The directory $(path) does not contain \
            snapshot sub-files in the HDF5 format"))
        )

        return vcat([findRealStars(sub_file) for sub_file in sub_files]...)

    else

        throw(ArgumentError("findRealStars: $(path) does not exist as a file or folder"))

    end

end

"""
    countStars(path::String)::Int

Count the number of stars in a snapshot, excluding wind particles.

# Arguments

  - `path::String`: Path to the snapshot file or folder.

# Returns

  - The number of stars.
"""
countStars(path::String)::Int = count(findRealStars(path))

"""
    findQtyExtrema(
        simulation_path::String,
        slice_n::Int,
        type_symbol::Symbol,
        block::String;
        <keyword arguments>
    )::NTuple{2,<:Number}

Compute the minimum and maximum values of `block`.

# Arguments

  - `simulation_path::String`: Path to the simulation directory, set in the code variable `OutputDir`.
  - `slice_n::Int`: Selects which snapshot to plot, starts at 1 and is independent of the number in the file name. If every snapshot is present, `slice_n` = filename_number + 1. If set to a negative number, the values in the whole simulation will be compared.
  - `type_symbol::Symbol`: Cell/particle type. The possibilities are the keys of [`PARTICLE_INDEX`](@ref).
  - `block::String`: Target block. The possibilities are the keys of [`QUANTITIES`](@ref).
  - `f::Function=identity`: A function with the signature:

    `f(data) -> values`

    where

      + `data::VecOrMat{<:Number}`: Data returned by [`getBlock`](@ref).
      + `values::Vector{<:Number}`: A vector with the values to be compared.
  - `warnings::Bool=true`: If a warning will be given when there is missing data.

# Returns

  - Tuple with the minimum and maximum values.
"""
function findQtyExtrema(
    simulation_path::String,
    slice_n::Int,
    type_symbol::Symbol,
    block::String;
    f::Function=identity,
    warnings::Bool=true,
)::NTuple{2,<:Number}

    (
        isdir(simulation_path) ||
        throw(ArgumentError("findQtyExtrema: $(simulation_path) does not exist as a directory"))
    )

    simulation_table = makeSimulationTable(simulation_path; warnings)

    if slice_n > 0

        # Get the number in the filename
        snap_n = safeSelect(simulation_table[!, :numbers], slice_n; warnings)

        # Check that after slicing there is one snapshot left
        (
            !isempty(snap_n) ||
            throw(ArgumentError("findQtyExtrema: There are no snapshots with `slice_n` = \
            $(slice_n), the contents of $(simulation_path) are: \n$(simulation_table)"))
        )

        # Find the target row and snapshot path
        snapshot_row = filter(:numbers => ==(lpad(snap_n, 3, "0")), simulation_table)
        snapshot_path = snapshot_row[1, :snapshot_paths]

        (
            !ismissing(snapshot_path) ||
            throw(ArgumentError("findQtyExtrema: The snapshot number $(slice_n) seems \
            to be missing"))

        )

        values = f(getBlock(snapshot_path, type_symbol, block))

        return extrema(values)

    end

    snapshot_paths = filter!(!ismissing, snapshot_row[!, :snapshot_paths])

    (
        !isempty(snapshot_paths) ||
        throw(ArgumentError("findQtyExtrema: I could not find any snapshots in $(simulation_path)"))
    )

    values = [f(getBlock(snapshot_path, type_symbol, block)) for snapshot_path in snapshot_paths]

    return extrema(Iterators.flatten(values))

end

@doc raw"""
    energyIntegrand(a::Real, header::SnapshotHeader)::Float64

The integrand of the integral that converts the scale factor into physical time:

```math
\frac{1}{H\,\sqrt{\mathcal{E}}} \, ,
```

where

```math
\mathcal{E} = \Omega_\Lambda + (1 - \Omega_\Lambda - \Omega_m) \, a^{-2} + \Omega_m \, a^{-3} \, ,
```
```math
H = H_0 \, a \, .
```

# Arguments

  - `a::Real`: Scale factor.
  - `header::SnapshotHeader`: Header of the relevant snapshot file.

# Returns

  - The integrand evaluated at `a`, in $\mathrm{Gyr}$.
"""
function energyIntegrand(a::Real, header::SnapshotHeader)::Float64

    # Return 0 if `a` = 0, as the integrand goes to 0 in the limit a -> 0.
    !iszero(a) || return 0.0

    # Compute Ω_K (curvature)
    omega_K = 1.0 - header.omega_0 - header.omega_l

    # Compute the energy function
    E = header.omega_0 / (a * a * a) + omega_K / (a * a) + header.omega_l

    # Compute the hubble constant in Gyr^-1
    H = header.h0 * HUBBLE_CONSTANT * a

    # Return the integrand, in Gyr
    return 1.0 / (H * sqrt(E))

end

"""
    computeProfile(
        positions::Matrix{<:Unitful.Length},
        quantity::Vector{<:Number},
        grid::CircularGrid;
        <keyword arguments>
    )::Vector{<:Number}

Compute a profile.

# Arguments

  - `positions::Matrix{<:Unitful.Length}`: Positions of the cells/particles. Each column is a cell/particle and each row a dimension.
  - `quantity::Vector{<:Number}`: The profile will be of this quantity.
  - `grid::CircularGrid`: Circular grid.
  - `norm_values::Vector{<:Number}=Number[]`: Values to normalize `quantity`.
  - `flat::Bool=true`: If the profile will be 2D, using rings, or 3D, using spherical shells.
  - `total::Bool=true`: If the sum (default) or the mean of `quantity` will be computed for each bin.
  - `cumulative::Bool=false`: If the profile will be accumulated or not.
  - `density::Bool=false`: If the profile will be of the density of `quantity`.
  - `empty_nan::Bool=true`: If empty bins will be set to NaN, 0 is used otherwise. Be carefull if `empty_nan` = true and `cumulative` = true, because every bin after the first NaN will be set to NaN.

# Returns

  - Vector with the values of the profile.
"""
function computeProfile(
    positions::Matrix{<:Unitful.Length},
    quantity::Vector{<:Number},
    grid::CircularGrid;
    norm_values::Vector{<:Number}=Number[],
    flat::Bool=true,
    total::Bool=true,
    cumulative::Bool=false,
    density::Bool=false,
    empty_nan::Bool=true,
)::Vector{<:Number}

    # Return a null profile if `quantity` is empty
    !isempty(quantity) || return fill(NaN, length(grid.grid))

    # Compute the distances of the cells/particles to the center of the grid
    if flat
        distances = computeDistance(positions[1:2, :]; center=grid.center[1:2])
    else
        distances = computeDistance(positions; center=grid.center)
    end

    # Compute the histogram of `quantity`
    if isempty(norm_values)

        profile = histogram1D(distances, quantity, grid; total, empty_nan)

    else

        quantity_histogram = histogram1D(distances, quantity, grid; total, empty_nan)
        norm_values_histogram = histogram1D(distances, norm_values, grid; total, empty_nan=false)

        replace!(x -> iszero(x) ? oneunit(x) : x, norm_values_histogram)

        profile = quantity_histogram ./ norm_values_histogram

    end

    region = flat ? grid.bin_areas : grid.bin_volumes

    if cumulative
        return density ? cumsum(profile) ./ cumsum(region) : cumsum(profile)
    end

    return density ? profile ./ region : profile

end

"""
    computeBandProfile(
        positions::Matrix{<:Unitful.Length},
        quantity::Vector{<:Number},
        grid::CircularGrid;
        <keyword arguments>
    )::NTuple{2,Vector{<:Number}}

Compute a profile of the mean and standard deviation of `quantity`.

Empty bins have a NaN mean and standard deviation.

# Arguments

  - `positions::Matrix{<:Unitful.Length}`: Positions of the cells/particles. Each column is a cell/particle and each row a dimension.
  - `quantity::Vector{<:Number}`: The profile will be of this quantity.
  - `grid::CircularGrid`: Circular grid.
  - `flat::Bool=true`: If the profile will be 2D, using rings, or 3D, using spherical shells.

# Returns

  - A tuple with two elements:

      + A vector with the mean value for each bin.
      + A vector with the standard deviation for each bin.
"""
function computeBandProfile(
    positions::Matrix{<:Unitful.Length},
    quantity::Vector{<:Number},
    grid::CircularGrid;
    flat::Bool=true,
)::NTuple{2,Vector{<:Number}}

    # Return a null profile if `quantity` is empty
    !isempty(quantity) || return fill(NaN, length(grid.grid))

    # Compute the distances of the cells/particles to the center of the grid
    if flat
        distances = computeDistance(positions[1:2, :]; center=grid.center[1:2])
    else
        distances = computeDistance(positions; center=grid.center)
    end

    # Compute the histogram of `quantity`
    histogram = listHistogram1D(distances, quantity, grid)

    return mean.(histogram), std.(histogram)

end

"""
    findHaloSubhalo(
        data_dict::Dict,
        star_idxs::Vector{Int},
        real_stars_idxs::Vector{Bool},
    )::NTuple{2,Vector{Int}}

Find in which halo and subhalo of `data_dict` each star in `star_idxs` was born.

For stars with no halo or subhalo, an index of -1 is given. The subhalo index is relative to the corresponding halo.

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
  - `star_idxs::Vector{Int}`: Indices of the target stars in `data_dict`.
  - `real_stars_idxs::Vector{Bool}`: Boolean list of stellar particles. True for a real star and false for a wind particle.

# Returns

  - A tuple with two elements:

      + A vector with the birth halo (index starting at 1) of each star (in the order of `star_idxs`).
      + A vector with the birth subhalo (index starting at 1) of each star (in the order of `star_idxs`).
"""
function findHaloSubhalo(
    data_dict::Dict,
    star_idxs::Vector{Int},
    real_stars_idxs::Vector{Bool},
)::NTuple{2,Vector{Int}}

    ################################################################################################
    # Read the subfind metadata
    ################################################################################################

    # Read the total number of halos
    n_halos = data_dict[:gc_data].header.n_groups_total

    # Read the number of subhalos in each halo
    n_subhalos_in_halo = data_dict[:group]["G_Nsubs"]

    # Read the number of stars in each halo
    n_stars_in_halo = data_dict[:group]["G_LenType"][PARTICLE_INDEX[:stars] + 1, :]

    # Read the number of stars in each subhalo
    n_stars_in_subhalo = data_dict[:subhalo]["S_LenType"][PARTICLE_INDEX[:stars] + 1, :]

    ################################################################################################
    # Allocate memory
    ################################################################################################

    # If each halo was the birth place of a star form `star_idxs`
    # So as to compute the relevant indices of each halo only once
    born_in_this_halo = fill(false, n_halos)

    # Index of the last real star particle belonging to each subhalo
    # Each element of this vector corresponds to a halo
    last_idxs_in_subhalo_list = Vector{Vector{Int}}(undef, n_halos)

    for i in 1:n_halos

        # Number of subhalos in halo `i`
        n_subfinds = n_subhalos_in_halo[i]

        # Index of the last real star particle belonging to each subhalo
        # Each element of this vector corresponds to a subhalo of halo `i`
        last_idxs_in_subhalo_list[i] = Vector{Int}(undef, n_subfinds)

    end

    # Output vectors
    halo_idxs = fill(-1, length(star_idxs))
    subhalo_idxs = fill(-1, length(star_idxs))

    ################################################################################################
    # Compute the index of the last real star particle belonging to each halo
    ################################################################################################

    # Compute the index of the last star/wind particle in each of the halos
    last_idxs_in_halo = cumsum(n_stars_in_halo)

    for (i, idx) in enumerate(last_idxs_in_halo)

        # Compute the number of wind particles up to the particle with index `idx`
        n_wind = count(x -> !(x), real_stars_idxs[1:idx])

        # Shift `last_idxs_in_halo` to ignore wind particles
        last_idxs_in_halo[i] = idx - n_wind

    end

    ############################################################################################
    # Compute in which halo and subhalo each star was born
    ############################################################################################
    @inbounds for (i, star_idx) in enumerate(star_idxs)

        ############################################################################################
        # Compute in which halo each star was born
        ############################################################################################

        # Find the halo where the target star was born
        halo_idx = searchsortedfirst(last_idxs_in_halo, star_idx)

        # If the star does not belong to any halo, leave the index as -1
        halo_idx <= length(last_idxs_in_halo) || continue

        halo_idxs[i] = halo_idx

        ############################################################################################
        # Compute in which subhalo each star was born
        ############################################################################################

        # Index of the last real star particle belonging to each subhalo
        # Each element of this vector corresponds to a subhalo of halo `halo_idx`
        last_idxs_in_subhalo = last_idxs_in_subhalo_list[halo_idx]

        # If it is the first time checking a star born in the halo `halo_idx`,
        # compute the index of the last real star particle in each subhalo of halo `halo_idx`
        if !born_in_this_halo[halo_idx]

            @inbounds if isone(halo_idx)
                # Absolute index of the first subhalo
                first_subhalo_abs_idx = 1

                # Absolute index of the last subhalo
                last_subhalo_abs_idx = n_subhalos_in_halo[1]

                # Number of stars up to, but not including, the halo `halo_idx`
                n_star_floor = 0
            else
                # Absolute index of the first subhalo
                first_subhalo_abs_idx = sum(n_subhalos_in_halo[1:(halo_idx - 1)]) + 1

                # Absolute index of the last subhalo
                last_subhalo_abs_idx = sum(n_subhalos_in_halo[1:halo_idx])

                # Number of stars up to, but not including, the halo `halo_idx`
                n_star_floor = sum(n_stars_in_halo[1:(halo_idx - 1)])
            end

            # Compute the index of the last star/wind particle in each of the subhalos of the halo `halo_idx`
            cumsum!(last_idxs_in_subhalo, n_stars_in_subhalo[first_subhalo_abs_idx:last_subhalo_abs_idx])
            map!(x -> x + n_star_floor, last_idxs_in_subhalo, last_idxs_in_subhalo)

            # Compute the index of the last real star particle in each of the subhalos of the halo `halo_idx`
            @inbounds for (i, idx) in enumerate(last_idxs_in_subhalo)

                # Compute the number of wind particles up to the particle with index `idx`
                n_wind = count(x -> !(x), real_stars_idxs[(n_star_floor + 1):idx])

                # Shift `last_idxs_in_subhalo` to ignore wind particles
                last_idxs_in_subhalo[i] = idx - n_wind

            end

            # Set to true to compute the relevant indices of this halo only once
            born_in_this_halo[halo_idx] = true

        end

        # Find the (relative index of the) subhalo where the target star was born
        subhalo_idx = searchsortedfirst(last_idxs_in_subhalo, star_idx)

        # If the star does not belong to any subhalo, leave the index as -1
        subhalo_idx <= length(last_idxs_in_subhalo) || continue

        subhalo_idxs[i] = subhalo_idx

    end

    return halo_idxs, subhalo_idxs

end

"""
    locateStellarBirthPlace(data_dict::Dict; <keyword arguments>)::NTuple{2,Vector{Int}}

Find in which halo and subhalo each star in `data_dict` was born.

For stars with no halo or subhalo, an index of -1 is given. The subhalo index is relative to the corresponding halo.

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
  - `warnings::Bool=true`: If a warning will be given when there is missing data.

# Returns

  - A tuple with two elements:

      + A vector with the birth halo (index starting at 1) of each star (in the order of `data_dict`).
      + A vector with the birth subhalo (index starting at 1) of each star (in the order of `data_dict`).
"""
function locateStellarBirthPlace(data_dict::Dict; warnings::Bool=true)::NTuple{2,Vector{Int}}

    ################################################################################################
    # Read the data in `data_dict`
    ################################################################################################

    # Read the birth time of each star
    birth_ticks = data_dict[:stars]["GAGE"]

    if data_dict[:sim_data].cosmological
        # Go from scale factor to physical time
        birth_times = computeTime(birth_ticks, data_dict[:snap_data].header)
    else
        birth_times = birth_ticks
    end

    # Read the time stamp of each snapshot
    times = data_dict[:sim_data].table[!, :physical_times]

    (
        length(times) >= 2 ||
        throw(ArgumentError("locateStellarBirthPlace: I found less that two snapshots in \
        $(data_dict[:sim_data].path). But I need more to locate the birth place of the stars."))
    )

    # Read the ID of each star
    ids = data_dict[:stars]["ID  "]

    # Compute the number of snapshots
    n_snaps = length(times)

    ################################################################################################
    # Compute the indices and IDs of the stars born between each of the snapshots
    ################################################################################################

    # Allocate memory
    present_star_idxs = [Int[] for _ in 1:n_snaps]

    @inbounds for (star_idx, birth_time) in enumerate(birth_times)

        snap_idx = searchsortedfirst(times, birth_time)

        @inbounds if snap_idx > n_snaps
            push!(present_star_idxs[n_snaps], star_idx)
        else
            push!(present_star_idxs[snap_idx], star_idx)
        end

    end

    # Read the IDs of the stars born between each of the snapshots
    star_ids = [ids[idxs] for idxs in present_star_idxs]

    ################################################################################################
    # Read each snapshot and find the original halo and subhalos of the stars born there
    ################################################################################################

    # Make a dataframe with the following columns:
    #   - 1. DataFrame index
    #   - 2. Number in the file name
    #   - 3. Scale factor
    #   - 4. Redshift
    #   - 5. Physical time
    #   - 6. Lookback time
    #   - 7. Snapshot path
    #   - 8. Group catalog path
    simulation_table = makeSimulationTable(data_dict[:sim_data].path; warnings)

    # Allocate memory
    birth_halo    = fill(-1, length(birth_times))
    birth_subhalo = fill(-1, length(birth_times))

    request = Dict(
        :stars => ["ID  "],
        :group => ["G_Nsubs", "G_LenType"],
        :subhalo => ["S_LenType"],
    )

    @inbounds for (global_idx, snapshot_row) in pairs(eachrow(simulation_table))

        # Select the IDs of the stars born in this snapshot
        ids = star_ids[global_idx]

        # Select the present index of the stars born in this snapshot
        present_idxs = present_star_idxs[global_idx]

        # Skip snapshots with no stars born in them
        !isempty(ids) || continue

        # Get the snapshot file path
        snapshot_path = snapshot_row[:snapshot_paths]

        # Get the group catalog file path
        groupcat_path = snapshot_row[:groupcat_paths]

        # Skip missing snapshots
        !ismissing(snapshot_path) || continue

        # Store the metadata of the current snapshot and simulation
        metadata = Dict(
            :sim_data => data_dict[:sim_data],
            :snap_data => Snapshot(
                snapshot_path,
                global_idx,
                global_idx,
                snapshot_row[:physical_times],
                snapshot_row[:lookback_times],
                snapshot_row[:scale_factors],
                snapshot_row[:redshifts],
                readSnapHeader(snapshot_path),
            ),
            :gc_data => GroupCatalog(
                groupcat_path,
                readGroupCatHeader(groupcat_path; warnings),
            ),
        )

        # Read the data in the snapshot
        past_data_dict = merge(
            metadata,
            readSnapshot(snapshot_path, request; warnings),
            readGroupCatalog(groupcat_path, snapshot_path, request; warnings),
        )

        # Get the birth index of the stars born in this snapshot
        past_idxs = parentIDToIndex(past_data_dict, ids)[:stars]

        # Sanity check
        (
            length(ids) == length(past_idxs) ||
            throw(DimensionMismatch("locateStellarBirthPlace: There are IDs in `ids` that are not \
            present in the birth snapshot or are from other cell/particle type. \
            This should be impossible!"))
        )

        # Find the halo and subhalo where each star was born, for the stars born in this snapshot
        halo_idxs, subhalo_idxs = findHaloSubhalo(
            past_data_dict,
            past_idxs,
            findRealStars(past_data_dict[:snap_data].path),
        )

        # Store the halo and subhalos indices in the current position of each star
        birth_halo[present_idxs] .= halo_idxs
        birth_subhalo[present_idxs] .= subhalo_idxs

    end

    return birth_halo, birth_subhalo

end
