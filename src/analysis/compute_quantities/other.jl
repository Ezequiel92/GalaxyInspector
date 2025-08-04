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
    sfr[idxs] .= data_dict[:stars]["MASS"][idxs] ./ age_resol

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
