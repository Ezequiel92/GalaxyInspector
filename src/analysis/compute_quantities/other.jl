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

@doc raw"""
    computeClumpingFactor(density::Vector{<:Number})::Float64

Compute the clumping factor,

```math
C_\rho = \frac{\rangle rho^2 \langle}{\rangle rho \langle^2} \, .
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

      + `:stellar_mass`              -> Stellar mass.
      + `:gas_mass`                  -> Gas mass.
      + `:hydrogen_mass`             -> Hydrogen mass.
      + `:dm_mass`                   -> Dark matter mass.
      + `:bh_mass`                   -> Black hole mass.
      + `:molecular_mass`            -> Molecular hydrogen (``\\mathrm{H_2}``) mass.
      + `:br_molecular_mass`         -> Molecular hydrogen (``\\mathrm{H_2}``) mass, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_mass`               -> Atomic hydrogen (``\\mathrm{HI}``) mass.
      + `:ionized_mass`              -> Ionized hydrogen (``\\mathrm{HII}``) mass.
      + `:neutral_mass`              -> Neutral hydrogen (``\\mathrm{HI + H_2}``) mass.
      + `:stellar_number`            -> Number of stellar particles.
      + `:gas_number`                -> Number of gas cells.
      + `:dm_number`                 -> Number of dark matter particles.
      + `:bh_number`                 -> Number of black hole particles.
      + `:molecular_fraction`        -> Gas mass fraction of molecular hydrogen.
      + `:br_molecular_fraction`     -> Gas mass fraction of molecular hydrogen, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_fraction`           -> Gas mass fraction of atomic hydrogen.
      + `:ionized_fraction`          -> Gas mass fraction of ionized hydrogen.
      + `:neutral_fraction`          -> Gas mass fraction of neutral hydrogen.
      + `:stellar_area_density`      -> Stellar area mass density, for a radius of `DISK_R`.
      + `:gas_area_density`          -> Gas mass surface density, for a radius of `DISK_R`.
      + `:molecular_area_density`    -> Molecular mass surface density, for a radius of `DISK_R`.
      + `:br_molecular_area_density` -> Molecular mass surface density, for a radius of `DISK_R`, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_area_density`       -> Atomic hydrogen area mass density, for a radius of `DISK_R`.
      + `:ionized_area_density`      -> Ionized hydrogen area mass density, for a radius of `DISK_R`.
      + `:neutral_area_density`      -> Neutral mass surface density, for a radius of `DISK_R`.
      + `:sfr_area_density`          -> Star formation rate area density, for the last `AGE_RESOLUTION` and a radius of `DISK_R`.
      + `:gas_metallicity`           -> Mass fraction of all elements above He in the gas (solar units).
      + `:stellar_metallicity`       -> Mass fraction of all elements above He in the stars (solar units).
      + `:X_gas_abundance`           -> Gas abundance of element ``\\mathrm{X}``, as ``12 + \\log_{10}(\\mathrm{X \\, / \\, H})``. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
      + `:X_stellar_abundance`       -> Stellar abundance of element ``\\mathrm{X}``, as ``12 + \\log_{10}(\\mathrm{X \\, / \\, H})``. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
      + `:stellar_specific_am`       -> Norm of the stellar specific angular momentum.
      + `:gas_specific_am`           -> Norm of the gas specific angular momentum.
      + `:dm_specific_am`            -> Norm of the dark matter specific angular momentum.
      + `:sfr`                       -> The star formation rate.
      + `:ssfr`                      -> The specific star formation rate.
      + `:observational_sfr`         -> The star formation rate of the last `AGE_RESOLUTION`.
      + `:observational_ssfr`        -> The specific star formation rate of the last `AGE_RESOLUTION`.
      + `:scale_factor`              -> Scale factor.
      + `:redshift`                  -> Redshift.
      + `:physical_time`             -> Physical time since the Big Bang.
      + `:lookback_time`             -> Physical time left to reach the last snapshot.

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

    elseif quantity == :br_molecular_mass

        integrated_qty = sum(computePressureMolecularMass(data_dict); init=0.0u"Msun")

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
        gas_mass = sum(data_dict[:gas]["MASS"]; init=0.0u"Msun")

        if iszero(gas_mass)
            integrated_qty = NaN
        else
            integrated_qty = molecular_mass / gas_mass
        end

    elseif quantity == :br_molecular_fraction

        molecular_mass = sum(computePressureMolecularMass(data_dict); init=0.0u"Msun")
        gas_mass = sum(data_dict[:gas]["MASS"]; init=0.0u"Msun")

        if iszero(gas_mass)
            integrated_qty = NaN
        else
            integrated_qty = molecular_mass / gas_mass
        end

    elseif quantity == :atomic_fraction

        atomic_mass = sum(computeAtomicMass(data_dict); init=0.0u"Msun")
        gas_mass = sum(data_dict[:gas]["MASS"]; init=0.0u"Msun")

        if iszero(gas_mass)
            integrated_qty = NaN
        else
            integrated_qty = atomic_mass / gas_mass
        end

    elseif quantity == :ionized_fraction

        ionized_mass = sum(computeIonizedMass(data_dict); init=0.0u"Msun")
        gas_mass = sum(data_dict[:gas]["MASS"]; init=0.0u"Msun")

        if iszero(gas_mass)
            integrated_qty = NaN
        else
            integrated_qty = ionized_mass / gas_mass
        end

    elseif quantity == :neutral_fraction

        neutral_mass = sum(computeNeutralMass(data_dict); init=0.0u"Msun")
        gas_mass = sum(data_dict[:gas]["MASS"]; init=0.0u"Msun")

        if iszero(gas_mass)
            integrated_qty = NaN
        else
            integrated_qty = neutral_mass / gas_mass
        end

    elseif quantity == :stellar_area_density

        integrated_qty = sum(data_dict[:stars]["MASS"]; init=0.0u"Msun") / area(DISK_R)

    elseif quantity == :gas_area_density

        integrated_qty = sum(data_dict[:gas]["MASS"]; init=0.0u"Msun") / area(DISK_R)

    elseif quantity == :molecular_area_density

        integrated_qty = sum(computeMolecularMass(data_dict); init=0.0u"Msun") / area(DISK_R)

    elseif quantity == :br_molecular_area_density

        integrated_qty = sum(computePressureMolecularMass(data_dict); init=0.0u"Msun") / area(DISK_R)

    elseif quantity == :atomic_area_density

        integrated_qty = sum(computeAtomicMass(data_dict); init=0.0u"Msun") / area(DISK_R)

    elseif quantity == :ionized_area_density

        integrated_qty = sum(computeIonizedMass(data_dict); init=0.0u"Msun") / area(DISK_R)

    elseif quantity == :neutral_area_density

        integrated_qty = sum(computeNeutralMass(data_dict); init=0.0u"Msun") / area(DISK_R)

    elseif quantity == :sfr_area_density

        sfr = sum(computeSFR(data_dict; age_resol=AGE_RESOLUTION); init=0.0u"Msun*yr^-1")

        integrated_qty = sfr / area(DISK_R)

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

      + `:stellar_mass`                -> Stellar mass.
      + `:gas_mass`                    -> Gas mass.
      + `:hydrogen_mass`               -> Hydrogen mass.
      + `:dm_mass`                     -> Dark matter mass.
      + `:bh_mass`                     -> Black hole mass.
      + `:molecular_mass`              -> Molecular hydrogen (``\\mathrm{H_2}``) mass.
      + `:br_molecular_mass`           -> Molecular hydrogen (``\\mathrm{H_2}``) mass, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_mass`                 -> Atomic hydrogen (``\\mathrm{HI}``) mass.
      + `:ionized_mass`                -> Ionized hydrogen (``\\mathrm{HII}``) mass.
      + `:neutral_mass`                -> Neutral hydrogen (``\\mathrm{HI + H_2}``) mass.
      + `:molecular_fraction`          -> Gas mass fraction of molecular hydrogen.
      + `:br_molecular_fraction`       -> Gas mass fraction of molecular hydrogen, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_fraction`             -> Gas mass fraction of atomic hydrogen.
      + `:ionized_fraction`            -> Gas mass fraction of ionized hydrogen.
      + `:neutral_fraction`            -> Gas mass fraction of neutral hydrogen.
      + `:molecular_neutral_fraction`  -> Fraction of molecular hydrogen in the neutral gas.
      + `:mol_eq_quotient`             -> Equilibrium quotient for the molecular fraction equation of the SF model.
      + `:ion_eq_quotient`             -> Equilibrium quotient for the ionized fraction equation of the SF model.
      + `:gas_mass_density`            -> Gas mass density.
      + `:hydrogen_mass_density`       -> Hydrogen mass density.
      + `:gas_number_density`          -> Gas number density.
      + `:molecular_number_density`    -> Molecular hydrogen number density.
      + `:br_molecular_number_density` -> Molecular hydrogen number density, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_number_density`       -> Atomic hydrogen number density.
      + `:ionized_number_density`      -> Ionized hydrogen number density.
      + `:neutral_number_density`      -> Neutral hydrogen number density.
      + `:gas_metallicity`             -> Mass fraction of all elements above He in the gas (solar units).
      + `:stellar_metallicity`         -> Mass fraction of all elements above He in the stars (solar units).
      + `:X_gas_abundance`             -> Gas abundance of element ``\\mathrm{X}``, as ``12 + \\log_{10}(\\mathrm{X \\, / \\, H})``. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
      + `:X_stellar_abundance`         -> Stellar abundance of element ``\\mathrm{X}``, as ``12 + \\log_{10}(\\mathrm{X \\, / \\, H})``. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
      + `:stellar_radial_distance`     -> Distance of every stellar particle to the origin.
      + `:gas_radial_distance`         -> Distance of every gas cell to the origin.
      + `:dm_radial_distance`          -> Distance of every dark matter particle to the origin.
      + `:stellar_xy_distance`         -> Projected distance of every stellar particle to the origin.
      + `:gas_xy_distance`             -> Projected distance of every gas cell to the origin.
      + `:dm_xy_distance`              -> Projected distance of every dark matter particle to the origin.
      + `:gas_sfr`                     -> SFR associated to each gas particle/cell within the code.
      + `:stellar_circularity`         -> Stellar circularity.
      + `:stellar_vcirc`               -> Stellar circular velocity.
      + `:stellar_vradial`             -> Stellar radial speed.
      + `:stellar_vtangential`         -> Stellar tangential speed.
      + `:stellar_vzstar`              -> Stellar speed in the z direction, computed as ``v_z \\, \\mathrm{sign}(z)``.
      + `:stellar_age`                 -> Stellar age.
      + `:sfr`                         -> The star formation rate.
      + `:ssfr`                        -> The specific star formation rate.
      + `:observational_sfr`           -> The star formation rate of the last `AGE_RESOLUTION`.
      + `:observational_ssfr`          -> The specific star formation rate of the last `AGE_RESOLUTION`.
      + `:temperature`                 -> Gas temperature, as ``\\log_{10}(T \\, / \\, \\mathrm{K})``.
      + `:pressure`                    -> Gas pressure.

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

    elseif quantity == :br_molecular_mass

        scatter_qty = computePressureMolecularMass(data_dict)

    elseif quantity == :atomic_mass

        scatter_qty = computeAtomicMass(data_dict)

    elseif quantity == :ionized_mass

        scatter_qty = computeIonizedMass(data_dict)

    elseif quantity == :neutral_mass

        scatter_qty = computeNeutralMass(data_dict)

    elseif quantity == :molecular_fraction

        molecular_mass = computeMolecularMass(data_dict)
        gas_mass = data_dict[:gas]["MASS"]

        scatter_qty = molecular_mass ./ gas_mass

    elseif quantity == :br_molecular_fraction

        molecular_mass = computePressureMolecularMass(data_dict)
        gas_mass = data_dict[:gas]["MASS"]

        scatter_qty = molecular_mass ./ gas_mass

    elseif quantity == :atomic_fraction

        atomic_mass = computeAtomicMass(data_dict)
        gas_mass = data_dict[:gas]["MASS"]

        scatter_qty = atomic_mass ./ gas_mass

    elseif quantity == :ionized_fraction

        ionized_mass = computeIonizedMass(data_dict)
        gas_mass = data_dict[:gas]["MASS"]

        scatter_qty = ionized_mass ./ gas_mass

    elseif quantity == :neutral_fraction

        neutral_mass = computeNeutralMass(data_dict)
        gas_mass = data_dict[:gas]["MASS"]

        scatter_qty = neutral_mass ./ gas_mass

    elseif quantity == :molecular_neutral_fraction

        molecular_mass = computeMolecularMass(data_dict)
        atomic_mass    = computeAtomicMass(data_dict)

        scatter_qty = molecular_mass ./ (atomic_mass .+ molecular_mass)

    elseif quantity == :mol_eq_quotient

        dg = data_dict[:gas]

        iterator = zip(
            dg["ETAD"],
            dg["FRAC"][2, :],
            dg["FRAC"][3, :],
            dg["FRAC"][4, :],
            τ_star.(dg["RHOC"] .* u"mp"),
            τ_cond.(dg["RHOC"] .* u"mp", dg["PARZ"]),
        )

        # Allocate memory
        scatter_qty = fill(NaN, length(dg["RHOC"]))

        for (i, (ηd, fa, fm, fs, τS, τC)) in enumerate(iterator)

            !(isnan(fa) || iszero(fa) || isone(fs) || iszero(fm)) || continue

            mol_ls = (fa / fm) * (1 - fs)
            mol_rs = uconvert(Unitful.NoUnits, ((ηd + 1) * τC) / τS)

            scatter_qty[i] = log10(mol_ls / mol_rs)

        end

    elseif quantity == :ion_eq_quotient

        dg = data_dict[:gas]

        iterator = zip(
            dg["ETAI"],
            dg["PARR"],
            dg["FRAC"][1, :],
            dg["FRAC"][3, :],
            GalaxyInspector.τ_star.(dg["RHOC"] .* u"mp"),
            GalaxyInspector.τ_rec.(dg["RHOC"] .* u"mp"),
        )

        # Allocate memory
        scatter_qty = fill(NaN, length(dg["RHOC"]))

        for (i, (ηi, R, fi, fm, τS, τR)) in enumerate(iterator)

            !(isnan(fi) || iszero(fi) || iszero(fm)) || continue

            ion_ls = (fi * fi) / fm
            ion_rs = uconvert(Unitful.NoUnits, ((ηi + R) * τR) / τS)

            scatter_qty[i] = log10(ion_ls / ion_rs)

        end

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

    elseif quantity == :br_molecular_number_density

        molecular_mass = computePressureMolecularMass(data_dict)
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

        if CODEBASE == :arepo

            scatter_qty = setPositive(data_dict[:gas]["GZ  "]) ./ SOLAR_METALLICITY

        elseif CODEBASE == :opengadget3

            metals = sum(setPositive(data_dict[:gas]["GMET"][METAL_LIST, :]); dims=1)
            scatter_qty = (metals ./ data_dict[:gas]["MASS"]) ./ SOLAR_METALLICITY

        else

            throw(ArgumentError("scatterQty: I don't recognize the codebase :$(CODEBASE)"))

        end

    elseif quantity == :stellar_metallicity

        if CODEBASE == :arepo

            scatter_qty = setPositive(data_dict[:stars]["GZ2 "]) ./ SOLAR_METALLICITY

        elseif CODEBASE == :opengadget3

            metals = sum(setPositive(data_dict[:stars]["GME2"][METAL_LIST, :]); dims=1)
            scatter_qty = (metals ./ data_dict[:stars]["MASS"]) ./ SOLAR_METALLICITY

        else

            throw(ArgumentError("scatterQty: I don't recognize the codebase :$(CODEBASE)"))

        end

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

    elseif quantity == :gas_sfr

        scatter_qty = data_dict[:gas]["SFR "]

    elseif quantity == :stellar_circularity

        @debug("scatterQty: The stellar circularity depends on the positions and velocities of all \
        cell/particles. So, after filtering, the result for a given star will change.")

        scatter_qty = computeCircularity(data_dict)

    elseif quantity == :stellar_vcirc

        @debug("scatterQty: The stellar circular velocity depends on the positions and velocities \
        of all cell/particles. So, after filtering, the result for a given star will change.")

        _, scatter_qty = computeVcirc(data_dict)

    elseif quantity == :stellar_vradial

        scatter_qty = computeVpolar(data_dict, :radial)

    elseif quantity == :stellar_vtangential

        scatter_qty = computeVpolar(data_dict, :tangential)

    elseif quantity == :stellar_vzstar

        scatter_qty = computeVpolar(data_dict, :zstar)

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

    elseif quantity == :pressure

        scatter_qty = data_dict[:gas]["PRES"]

    else

        throw(ArgumentError("scatterQty: I don't recognize the quantity :$(quantity)"))

    end

    if isempty(scatter_qty)
        return Number[]
    end

    return scatter_qty

end
