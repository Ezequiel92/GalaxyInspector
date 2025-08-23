####################################################################################################
# Computation of derived quantities
####################################################################################################

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

        # For cosmological simulations, the time field in the Header of the snapshot is the scale factor
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
        internal_energy::Vector{<:SpecificEnergy},
        electron_fraction::Vector{Float32},
    )::Vector{<:Unitful.Temperature}

Compute the gas temperature.

# Arguments

  - `internal_energy::Vector{<:SpecificEnergy}`: Specific internal energy of every gas cell/particle.
  - `electron_fraction::Vector{Float32}`: Number fraction of electrons in every gas cell/particle.

# Returns

  - The temperature of each gas cell/particle.
"""
function computeTemperature(
    internal_energy::Vector{<:SpecificEnergy},
    electron_fraction::Vector{Float32},
)::Vector{<:Unitful.Temperature}

    # xH := mass_fraction_of_hydrogen
    xH = HYDROGEN_MASSFRAC

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

    birth_ticks = data_dict[:stellar]["GAGE"]

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
    computeSFR(
        data_dict::Dict;
        <keyword arguments>
    )::Vector{<:Unitful.MassFlow}

Compute the star formation rate of each stellar particle.

For stellar particles younger than `age_resol`, the SFR is its mass divided by `age_resol`. It is defined as 0 for older particles.

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
    sfr[idxs] .= data_dict[:stellar]["MASS"][idxs] ./ age_resol

    return sfr

end

@doc raw"""
    computeClumpingFactor(density::Vector{<:Number})::Float64

Compute the clumping factor,

```math
C_\rho = \frac{\langle \rho^2 \rangle}{\langle \rho \rangle^2} \, .
```

# Arguments

  - `density::Vector{<:Number}`: The density of the cells/particles.

# Returns

  - The clumping factor.
"""
function computeClumpingFactor(density::Vector{<:Number})::Float64

    !isempty(density) || return NaN

    μ, var = mean_and_var(density)

    return 1.0 + uconvert(Unitful.NoUnits, var / μ^2)

end

@doc raw"""
    computeDepletionTime(
        mass::Vector{<:Unitful.Mass},
        sfr::Vector{<:Unitful.MassFlow},
    )::Vector{<:Unitful.Time}

Compute the depletion time,

```math
t_\mathrm{ff} = \frac{M_\mathrm{gas}}{\dot{M}_\star} \, .
```

# Arguments

  - `mass::Vector{<:Unitful.Mass}`: The gas mass of the cells/particles.
  - `sfr::Vector{<:Unitful.MassFlow}`: The SFR associated to each cell/particle.

# Returns

  - The depletion time.
"""
function computeDepletionTime(
    mass::Vector{<:Unitful.Mass},
    sfr::Vector{<:Unitful.MassFlow},
)::Vector{<:Unitful.Time}

    if any(isempty, [mass, sfr])

        (
            !logging[] ||
            @warn("computeDepletionTime: There is missing data for the gas, so the result will be \
            an empty vector")
        )

        return Unitful.Time[]

    end

    return @. mass / sfr

end

@doc raw"""
    computeEfficiencyFF(
        density::Vector{<:Unitful.Density},
        mass::Vector{<:Unitful.Mass},
        sfr::Vector{<:Unitful.MassFlow},
    )::Vector{Float64}

Compute the star formation efficiency per free-fall time, according to the definition in equation 1 of Krumholz (2012),

```math
\epsilon_\mathrm{ff} = \frac{t_\mathrm{ff}}{t_\mathrm{dep}} \, .
```
where

```math
t_\mathrm{ff} = \sqrt{\frac{3 \, \pi}{32 \, G \, \rho}} \, ,
```
is the free-fall time, and

```math
t_\mathrm{dep} = \frac{M_\mathrm{H_2}}{\dot{M}_\star} \, ,
```
is the depletion time.

# Arguments

  - `density::Vector{<:Unitful.Density}`: The molecular hydrogen (``\\mathrm{H_2}``) density of the cells/particles.
  - `mass::Vector{<:Unitful.Mass}`: The gas mass of the cells/particles.
  - `sfr::Vector{<:Unitful.MassFlow}`: The SFR associated to each cell/particle.

# Returns

  - The star formation efficiency per free-fall time.

# References

M. R. Krumholz et al. (2012). *A UNIVERSAL, LOCAL STAR FORMATION LAW IN GALACTIC CLOUDS, NEARBY GALAXIES, HIGH-REDSHIFT DISKS, AND STARBURSTS*. The Astrophysical Journal, **745(1)**, 69. [doi:10.1088/0004-637X/745/1/69](https://doi.org/10.1088/0004-637X/745/1/69)
"""
function computeEfficiencyFF(
    density::Vector{<:Unitful.Density},
    mass::Vector{<:Unitful.Mass},
    sfr::Vector{<:Unitful.MassFlow},
)::Vector{Float64}

    if any(isempty, [density, mass, sfr])

        (
            !logging[] ||
            @warn("computeEfficiencyFF: There is missing data for the gas, so the result will be \
            an empty vector")
        )

        return Float64[]

    end

    ϵff = Vector{Float64}(undef, length(density))

    for i in eachindex(ϵff)

        m = mass[i]
        ρ = density[i]
        sf = sfr[i]

        if iszero(m) || iszero(sf) || ρ < THRESHOLD_DENSITY

            ϵff[i] = NaN

        else

            # Compute the free-fall time
            tff = sqrt(3π / (32 * Unitful.G * ρ))

            # Compute the depletion time
            tdep = m / sf

            # Compute the star formation efficiency per free-fall time
            ϵff[i] = uconvert(Unitful.NoUnits, tff / tdep)

        end

    end

    return ϵff

end

"""
    computeKineticEnergy(
        masses::Vector{<:Unitful.Mass},
        velocities::Matrix{<:Unitful.Velocity},
    )::Vector{<:Unitful.Energy}

Compute the kinetic energy.

# Arguments

  - `masses::Vector{<:Unitful.Mass}`: Mass of every cell/particle.
  - `velocities::Matrix{<:Unitful.Velocity}`: Velocities of the cells/particles. Each column is a cell/particle and each row a dimension.

# Returns

  - The kinetic energy of each cell/particle.
"""
function computeKineticEnergy(
    masses::Vector{<:Unitful.Mass},
    velocities::Matrix{<:Unitful.Velocity},
)::Vector{<:Unitful.Energy}

    if any(isempty, [masses, velocities])

        (
            !logging[] ||
            @warn("computeKineticEnergy: There is missing data, so the result will be \
            an empty vector")
        )

        return Unitful.Energy[]

    end

    return [0.5 * mass * norm(vel)^2 for (mass, vel) in zip(masses, eachcol(velocities))]

end

"""
    computePotentialEnergy(
        potential::Vector{<:SpecificEnergy},
        masses::Vector{<:Unitful.Mass},
    )::Vector{<:Unitful.Energy}

Compute the gravitational potencial energy.

# Arguments

  - `potential::Vector{<:SpecificEnergy}`: Specific potential energy of every cell/particle.
  - `masses::Vector{<:Unitful.Mass}`: Mass of every cell/particle.

# Returns

  - The gravitational potencial energy of each cell/particle.
"""
function computePotentialEnergy(
    potential::Vector{<:SpecificEnergy},
    masses::Vector{<:Unitful.Mass},
)::Vector{<:Unitful.Energy}

    if any(isempty, [potential, masses])

        (
            !logging[] ||
            @warn("computePotentialEnergy: There is missing data, so the result will be \
            an empty vector")
        )

        return Unitful.Energy[]

    end

    return potential .* masses

end

"""
    computeTotalEnergy(
        potential::Vector{<:SpecificEnergy},
        masses::Vector{<:Unitful.Mass},
        velocities::Matrix{<:Unitful.Velocity},
    )::Vector{<:Unitful.Energy}

Compute the total energy (kinetic + potential).

# Arguments

  - `potential::Vector{<:SpecificEnergy}`: Specific potential energy of every cell/particle.
  - `masses::Vector{<:Unitful.Mass}`: Mass of every cell/particle.
  - `velocities::Matrix{<:Unitful.Velocity}`: Velocities of the cells/particles. Each column is a cell/particle and each row a dimension.

# Returns

  - The total energy of each cell/particle.
"""
function computeTotalEnergy(
    potential::Vector{<:SpecificEnergy},
    masses::Vector{<:Unitful.Mass},
    velocities::Matrix{<:Unitful.Velocity},
)::Vector{<:Unitful.Energy}

    if any(isempty, [potential, masses, velocities])

        (
            !logging[] ||
            @warn("computeTotalEnergy: There is missing data, so the result will be \
            an empty vector")
        )

        return Unitful.Energy[]

    end

    return computePotentialEnergy(potential, masses) .+ computeKineticEnergy(masses, velocities)

end

"""
    density3DProjection(
        data_dict::Dict,
        grid::CubicGrid,
        quantity::Symbol,
        type::Symbol;
        <keyword arguments>
    )::Array{Float64,3}

Sample the 3D density field of a given quantity using a cubic grid

!!! note

    If the source of the field are particles, a simple 3D histogram is used. If they are Voronoi cells instead, the density of the cell that intersects each voxel is used.

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
  - `quantity::Symbol`: Which density will be calculated. The options are:

      + `:stellar_mass`      -> Stellar density.
      + `:gas_mass`          -> Gas density.
      + `:hydrogen_mass`     -> Hydrogen density.
      + `:dm_mass`           -> Dark matter density.
      + `:bh_mass`           -> Black hole density.
      + `:molecular_mass`    -> Molecular hydrogen (``\\mathrm{H_2}``) density.
      + `:br_molecular_mass` -> Molecular hydrogen (``\\mathrm{H_2}``) density, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_mass`       -> Atomic hydrogen (``\\mathrm{HI}``) density.
      + `:ionized_mass`      -> Ionized hydrogen (``\\mathrm{HII}``) density.
      + `:neutral_mass`      -> Neutral hydrogen (``\\mathrm{HI + H_2}``) density.
      + `:stellar_gas_mass`  -> Stellar gas mass (according to our SF model).
      + `:ode_metal_mass`    -> Metal mass (according to our SF model).
      + `:dust_mass`         -> Dust mass.
  - `type::Symbol`: If the source of the field are `:particles` or Voronoi `:cells`.
  - `m_unit::Unitful.Units=u"Msun"`: Mass unit.
  - `l_unit::Unitful.Units=u"kpc"`: Length unit.
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

  - A 3D array with the density at each point of the 3D grid.

# References

L. Blitz et al. (2006). *The Role of Pressure in GMC Formation II: The H2-Pressure Relation*. The Astrophysical Journal, **650(2)**, 933. [doi:10.1086/505417](https://doi.org/10.1086/505417)
"""
function density3DProjection(
    data_dict::Dict,
    grid::CubicGrid,
    quantity::Symbol,
    type::Symbol;
    m_unit::Unitful.Units=u"Msun",
    l_unit::Unitful.Units=u"kpc",
    filter_function::Function=filterNothing,
)::Array{Float64,3}

    filtered_dd = filterData(data_dict; filter_function)

    # Set the cell/particle type
    if quantity ∈ [
        :gas_mass,
        :hydrogen_mass,
        :molecular_mass,
        :br_molecular_mass,
        :atomic_mass,
        :ionized_mass,
        :neutral_mass,
        :stellar_gas_mass,
        :ode_metal_mass,
        :dust_mass,
    ]
        component = :gas
    elseif quantity == :stellar_mass
        component = :stellar
    elseif quantity == :dm_mass
        component = :dark_matter
    elseif quantity == :bh_mass
        component = :black_hole
    else
        throw(ArgumentError("density3DProjection: I don't recognize the quantity :$(quantity)"))
    end

    # For comological simulations with comoving units, correct
    # the density so it is always in physical units
    if !PHYSICAL_UNITS && data_dict[:sim_data].cosmological
        # Correction factor for the volume
        # V [physical units] = V [comoving units] * a0^3
        physical_factor = data_dict[:snap_data].scale_factor^3
    else
        physical_factor = 1.0
    end

    # Load the cell/particle positions
    positions = filtered_dd[component]["POS "]

    # Compute the masses of the target quantity
    masses = scatterQty(filtered_dd, quantity)

    # If any of the necessary quantities are missing return an empty density field
    if any(isempty, [masses, positions])
        return fill(NaN, (grid.n_bins, grid.n_bins, grid.n_bins))
    end

    if type == :cells

        # Compute the volume of each cell
        cell_volumes = filtered_dd[component]["MASS"] ./ filtered_dd[component]["RHO "]

        # Compute the densities of the target quantity
        densities = ustrip.(m_unit * l_unit^-3, masses ./ cell_volumes)

        # Allocate memory
        physical_grid = Matrix{Float64}(undef, 3, grid.n_bins^3)

        # Compute the tree for a nearest neighbor search
        kdtree = KDTree(ustrip.(l_unit, positions))

        # Reshape the grid to conform to the way `nn` expect the matrix to be structured
        for i in eachindex(grid.grid)
            physical_grid[1, i] = ustrip(l_unit, grid.grid[i][1])
            physical_grid[2, i] = ustrip(l_unit, grid.grid[i][2])
            physical_grid[3, i] = ustrip(l_unit, grid.grid[i][3])
        end

        # Find the nearest cell to each voxel
        idxs, _ = nn(kdtree, physical_grid)

        # Allocate memory
        density = similar(grid.grid, Float64)

        # Compute the density in each voxel
        for i in eachindex(grid.grid)
            density[i] = densities[idxs[i]]
        end

        # Set bins with a value of 0 to NaN
        replace!(x -> iszero(x) ? NaN : x, density)

    elseif type == :particles

        # Compute the 3D histogram
        density = ustrip.(
            m_unit * l_unit^-3,
            histogram3D(positions, masses, grid; empty_nan=true) ./ grid.bin_volume,
        )

    else

        throw(ArgumentError("density3DProjection: The argument `type` must be :cells or \
        :particles, but I got :$(type)"))

    end

    if logging[]

        log_density = filter(!isnan, log10.(density))

        if isempty(log_density)

            min_max_ρ = (NaN, NaN)
            mean_ρ    = NaN
            median_ρ  = NaN
            mode_ρ    = NaN

        else

            min_max_ρ = extrema(log_density)
            mean_ρ    = mean(log_density)
            median_ρ  = median(log_density)
            mode_ρ    = mode(log_density)

        end

        # Print the density range
        @info(
            "\nDensity range - log₁₀(ρ [$(m_unit * l_unit^-3)]) \
            \n  Simulation: $(basename(filtered_dd[:sim_data].path)) \
            \n  Snapshot:   $(filtered_dd[:snap_data].global_index) \
            \n  Quantity:   $(quantity) \
            \n  Type:       $(type) \
            \n  Min - Max:  $(min_max_ρ) \
            \n  Mean:       $(mean_ρ) \
            \n  Median:     $(median_ρ) \
            \n  Mode:       $(mode_ρ)"
        )

    end

    return density

end

"""
    computeParticleProfile(
        positions::Matrix{<:Unitful.Length},
        quantity::Vector{<:Number},
        grid::CircularGrid;
        <keyword arguments>
    )::Vector{<:Number}

Compute a profile, using an 1D histogram.

# Arguments

  - `positions::Matrix{<:Unitful.Length}`: Positions of the cells/particles. Each column is a cell/particle and each row a dimension.
  - `quantity::Vector{<:Number}`: The profile will be of this quantity.
  - `grid::CircularGrid`: Circular grid.
  - `norm_values::Vector{<:Number}=Number[]`: Values to normalize `quantity`.
  - `flat::Bool=true`: If the profile will be 2D, using rings, or 3D, using spherical shells.
  - `total::Bool=true`: If the sum (default) or the mean of `quantity` will be computed for each bin.
  - `cumulative::Bool=false`: If the profile will be accumulated or not.
  - `density::Bool=false`: If the profile will be of the density of `quantity`.
  - `empty_nan::Bool=true`: If empty bins will be set to NaN, 0 is used otherwise. Be careful if `empty_nan` = true and `cumulative` = true, because every bin after the first NaN will be set to NaN.

# Returns

  - Vector with the values of the profile.
"""
function computeParticleProfile(
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
    if isempty(quantity)

        (
            !logging[] ||
            @warn("computeParticleProfile: The vector `quantity` is empty. The profile will be \
            filled with NaNs")
        )

        return fill(NaN, length(grid.grid))

    end

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
    computeParticleBandProfile(
        positions::Matrix{<:Unitful.Length},
        quantity::Vector{<:Number},
        grid::CircularGrid;
        <keyword arguments>
    )::NTuple{2,Vector{<:Number}}

Compute a profile of the mean and standard deviation of `quantity`, using an 1D histogram

Empty bins have NaN as mean and standard deviation.

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
function computeParticleBandProfile(
    positions::Matrix{<:Unitful.Length},
    quantity::Vector{<:Number},
    grid::CircularGrid;
    flat::Bool=true,
)::NTuple{3,Vector{<:Number}}

    if isempty(quantity)

        (
            !logging[] ||
            @warn("computeParticleBandProfile: The vector `quantity` is empty. The profile will be \
            filled with NaNs")
        )

        return fill(NaN, length(grid.grid))

    end

    # Compute the distances of the cells/particles to the center of the grid
    if flat
        distances = computeDistance(positions[1:2, :]; center=grid.center[1:2])
    else
        distances = computeDistance(positions; center=grid.center)
    end

    # Compute the histogram of `quantity`
    histogram = listHistogram1D(distances, quantity, grid)

    return quantile.(histogram, 0.5), quantile.(histogram, 0.25), quantile.(histogram, 0.75)

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
    H = header.h * HUBBLE_CONSTANT * a

    # Return the integrand, in Gyr
    return 1.0 / (H * sqrt(E))

end

@doc raw"""
    bigiel2008(
        ΣH::Vector{<:SurfaceDensity};
        <keyword arguments>
    )::Vector{<:Number}

Kennicutt-Schmidt law for the molecular or neutral gas, taken from a set of observations of nearby galaxies.

From Bigiel et al. (2008) (Section 3.1), we have

```math
\Sigma_\mathrm{SFR} = a \left( \frac{\Sigma_\mathrm{HI, H_2, gas}}{10 \, \mathrm{M_\odot \, pc^{-2}}} \right)^{\!N} \, ,
```
where N is the power-law index, and $A = \log_{10}(a)$, where $a$ is $\Sigma_\mathrm{SFR}$ at the fiducial gas surface density of $10 \, \mathrm{M_\odot \, pc^{-2}}$.

# Arguments

  - `ΣH::Vector{<:SurfaceDensity}`: Values of the molecular or neutral gas surface density, with units.
  - `molecular::Bool=true`: If the x axis will be the area mass density of molecular hydrogen, or, if set to false, the area mass density of neutral hydrogen.
  - `log_output::Bool=true`: If the output will the $\log_{10}$ of the star formation area density, or the star formation area density itself (with units). If `log_output` = true, the implied unit is $\mathrm{M_\odot \, yr^{-1} \, kpc^{-2}}$

# Returns

  - The star formation area density.

# References

F. Bigiel et al. (2008). *THE STAR FORMATION LAW IN NEARBY GALAXIES ON SUB-KPC SCALES*. The Astrophysical Journal, **136(6)**, 2846. [doi:10.1088/0004-6256/136/6/2846](https://doi.org/10.1088/0004-6256/136/6/2846)
"""
function bigiel2008(
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
    invBigiel2008(
        Σsfr ::Vector{<:MassFlowDensity};
        <keyword arguments>
    )::Vector{<:Number}

Inverse Kennicutt-Schmidt law for the molecular or neutral gas, taken from a set of observations of nearby galaxies.

From Bigiel et al. (2008) (Section 3.1, Eq. 2), we have

```math
\Sigma_\mathrm{SFR} = a \left( \frac{\Sigma_\mathrm{HI, H_2, gas}}{10 \, \mathrm{M_\odot \, pc^{-2}}} \right)^{\!N} \, ,
```
where N is the power-law index, and $A = \log_{10}(a)$, where $a$ is $\Sigma_\mathrm{SFR}$ at the fiducial gas surface density of $10 \, \mathrm{M_\odot \, pc^{-2}}$.

# Arguments

  - `Σsfr ::Vector{<:MassFlowDensity}`: Values of the star formation area density, with units.
  - `molecular::Bool=true`: If the output will be the area mass density of molecular hydrogen, or, if set to false, the area mass density of neutral hydrogen.
  - `log_output::Bool=true`: If the output will the $\log_{10}$ of the molecular or neutral gas surface density, or the molecular or neutral gas surface density itself (with units). If `log_output` = true, the implied unit is $10 \, \mathrm{M_\odot \, pc^{-2}}$

# Returns

  - The molecular or neutral gas surface density.

# References

F. Bigiel et al. (2008). *THE STAR FORMATION LAW IN NEARBY GALAXIES ON SUB-KPC SCALES*. The Astrophysical Journal, **136(6)**, 2846. [doi:10.1088/0004-6256/136/6/2846](https://doi.org/10.1088/0004-6256/136/6/2846)
"""
function invBigiel2008(
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
    kennicutt1998(Σgas::Vector{<:SurfaceDensity}; <keyword arguments>)::Vector{<:Number}

Kennicutt-Schmidt law, taken from a set of observations of nearby galaxies.

From Kennicutt (1998) (Section 4, Eq. 4), we have

```math
\Sigma_\mathrm{SFR} = a \left( \frac{\Sigma_\mathrm{gas}}{1 \, \mathrm{M_\odot \, pc^{-2}}} \right)^{\!N} \mathrm{M_\odot \, yr^{-1] \, kpc^{-2}} \, ,
```
where N is the power-law index and $a$ is $\Sigma_\mathrm{SFR}$ at the fiducial gas surface density of $1 \, \mathrm{M_\odot \, pc^{-2}}$.

# Arguments

  - `Σgas::Vector{<:SurfaceDensity}`: Values of the gas mass surface density, with units.
  - `log_output::Bool=true`: If the output will the $\log_{10}$ of the star formation area density, or the star formation area density itself (with units). If `log_output` = true, the implied unit is $\mathrm{M_\odot \, yr^{-1} \, kpc^{-2}}$

# Returns

  - The star formation area density.

# References

R. C. Kennicutt (1998). *The Global Schmidt Law in Star-forming Galaxies*. The Astrophysical Journal, **498(2)**, 541-552. [doi:10.1086/305588](https://doi.org/10.1086/305588)
"""
function kennicutt1998(Σgas::Vector{<:SurfaceDensity}; log_output::Bool=true)::Vector{<:Number}

    log10Σgas = @. log10(ustrip(u"Msun * pc^-2", Σgas))
    log10Σsfr = @. log10(a_KS98) + log10Σgas * N_KS98

    if log_output
        return log10Σsfr
    else
        return @. exp10(log10Σsfr) * u"Msun * yr^-1 * kpc^-2"
    end

end

@doc raw"""
    invKennicutt1998(Σsfr::Vector{<:MassFlowDensity}; <keyword arguments>)::Vector{<:Number}

Inverse Kennicutt-Schmidt law, taken from a set of observations of nearby galaxies.

From Kennicutt (1998) (Section 4, Eq. 4), we have

```math
\Sigma_\mathrm{SFR} = a \left( \frac{\Sigma_\mathrm{gas}}{1 \, \mathrm{M_\odot \, pc^{-2}}} \right)^{\!N} \mathrm{M_\odot \, yr^{-1] \, kpc^{-2}} \, ,
```
where N is the power-law index and $a$ is $\Sigma_\mathrm{SFR}$ at the fiducial gas surface density of $1 \, \mathrm{M_\odot \, pc^{-2}}$.

# Arguments

  - `Σsfr::Vector{<:MassFlowDensity}`: Values of the star formation area density, with units.
  - `log_output::Bool=true`: If the output will the $\log_{10}$ of the gas mass surface density, or the gas mass surface density itself (with units). If `log_output` = true, the implied unit is $\mathrm{M_\odot \, pc^{-2}}$

# Returns

  - The gas mass surface density.

# References

R. C. Kennicutt (1998). *The Global Schmidt Law in Star-forming Galaxies*. The Astrophysical Journal, **498(2)**, 541-552. [doi:10.1086/305588](https://doi.org/10.1086/305588)
"""
function invKennicutt1998(Σsfr::Vector{<:MassFlowDensity}; log_output::Bool=true)::Vector{<:Number}

    log10Σsfr = @. log10(ustrip(u"Msun * yr^-1 * kpc^-2", Σsfr))
    log10Σgas = @. (log10Σsfr - log10(a_KS98)) / N_KS98

    if log_output
        return log10Σgas
    else
        return @. exp10(log10Σgas) * u"Msun * pc^-2"
    end

end
