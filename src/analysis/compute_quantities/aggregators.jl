####################################################################################################
# Aggregator functions
####################################################################################################

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
      + `:stellar_gas_mass`          -> Stellar gas mass (according to our SF model).
      + `:ode_metal_mass`            -> Metal mass (according to our SF model).
      + `:ode_metallicity`           -> Metallicity (according to our SF model).
      + `:dust_mass`                 -> Dust mass.
      + `:stellar_number`            -> Number of stellar particles.
      + `:gas_number`                -> Number of gas cells.
      + `:dm_number`                 -> Number of dark matter particles.
      + `:bh_number`                 -> Number of black hole particles.
      + `:molecular_fraction`        -> Gas mass fraction of molecular hydrogen.
      + `:br_molecular_fraction`     -> Gas mass fraction of molecular hydrogen, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_fraction`           -> Gas mass fraction of atomic hydrogen.
      + `:ionized_fraction`          -> Gas mass fraction of ionized hydrogen.
      + `:neutral_fraction`          -> Gas mass fraction of neutral hydrogen.
      + `:molecular_neutral_fraction`-> Fraction of molecular hydrogen in the neutral gas.
      + `:ionized_neutral_fraction`  -> Fraction of ionized gas to neutral gas.
      + `:gas_mass_density`          -> Mean gas mass density.
      + `:stellar_gas_fraction`      -> Stellar gas fraction (according to our SF model).
      + `:metal_gas_fraction`        -> Metallicity (according to our SF model).
      + `:dust_fraction`             -> Dust mass fraction.
      + `:stellar_area_density`      -> Stellar area mass density, for a radius of `DISK_R`.
      + `:gas_area_density`          -> Gas mass surface density, for a radius of `DISK_R`.
      + `:molecular_area_density`    -> Molecular mass surface density, for a radius of `DISK_R`.
      + `:br_molecular_area_density` -> Molecular mass surface density, for a radius of `DISK_R`, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_area_density`       -> Atomic hydrogen area mass density, for a radius of `DISK_R`.
      + `:ionized_area_density`      -> Ionized hydrogen area mass density, for a radius of `DISK_R`.
      + `:neutral_area_density`      -> Neutral mass surface density, for a radius of `DISK_R`.
      + `:sfr_area_density`          -> Star formation rate area density, for the last `AGE_RESOLUTION` and a radius of `DISK_R`.
      + `:gas_td`                    -> Gas depletion time.
      + `:molecular_td`              -> The molecular hydrogen (``\\mathrm{H_2}``) depletion time.
      + `:br_molecular_td`           -> The Molecular hydrogen (``\\mathrm{H_2}``) depletion time, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_td`                 -> The atomic hydrogen (``\\mathrm{HI}``) depletion time.
      + `:ionized_td`                -> The ionized hydrogen (``\\mathrm{HII}``) depletion time.
      + `:neutral_td`                -> The neutral hydrogen (``\\mathrm{HI + H_2}``) depletion time.
      + `:gas_metallicity`           -> Mass fraction of all elements above He in the gas (solar units).
      + `:stellar_metallicity`       -> Mass fraction of all elements above He in the stars (solar units).
      + `:X_gas_abundance`           -> Gas abundance of element ``\\mathrm{X}``, as ``12 + \\log_{10}(\\mathrm{X \\, / \\, H})``. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
      + `:X_stellar_abundance`       -> Stellar abundance of element ``\\mathrm{X}``, as ``12 + \\log_{10}(\\mathrm{X \\, / \\, H})``. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
      + `:stellar_specific_am`       -> Norm of the stellar specific angular momentum.
      + `:gas_specific_am`           -> Norm of the gas specific angular momentum.
      + `:dm_specific_am`            -> Norm of the dark matter specific angular momentum.
      + `:sfr`                       -> Star formation rate.
      + `:ssfr`                      -> The specific star formation rate.
      + `:observational_sfr`         -> Star formation rate of the last `AGE_RESOLUTION`.
      + `:observational_ssfr`        -> The specific star formation rate of the last `AGE_RESOLUTION`.
      + `:stellar_eff`               -> Star formation efficiency per free-fall time for the gas that has turn into stars.
      + `:gas_eff`                   -> Star formation efficiency per free-fall time for the gas.
      + `:molecular_eff`             -> Star formation efficiency per free-fall time for the molecular hydrogen (``\\mathrm{H_2}``) gas.
      + `:br_molecular_eff`          -> Star formation efficiency per free-fall time for the molecular hydrogen (``\\mathrm{H_2}``) gas, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_eff`                -> Star formation efficiency per free-fall time for the atomic hydrogen (``\\mathrm{HI}``) gas.
      + `:ionized_eff`               -> Star formation efficiency per free-fall time for the ionized hydrogen (``\\mathrm{HII}``) gas.
      + `:neutral_eff`               -> Star formation efficiency per free-fall time for the neutral hydrogen (``\\mathrm{HI + H_2}``) gas.
      + `:scale_factor`              -> Scale factor.
      + `:redshift`                  -> Redshift.
      + `:physical_time`             -> Physical time since the Big Bang.
      + `:lookback_time`             -> Physical time left to reach the last snapshot.
      + `:ode_gas_it`                -> Integration time.
      + `:ode_gas_tau_s`             -> Star formation time scale, ``\\tau_\\mathrm{S}``.
      + `:ode_gas_eta_d`             -> Photodissociation efficiency, ``\\eta_\\mathrm{diss}``.
      + `:ode_gas_eta_i`             -> Photoionization efficiency, ``\\eta_\\mathrm{ion}``.
      + `:ode_gas_r`                 -> Mass recycling parameter, ``R``.
      + `:ode_gas_cold_mf`           -> Cold gas mass fraction.
      + `:ode_stellar_it`            -> Integration time, for the gas that form the stars.
      + `:ode_stellar_tau_s`         -> Star formation time scale, ``\\tau_\\mathrm{S}``, for the gas that form the stars.
      + `:ode_stellar_eta_d`         -> Photodissociation efficiency, ``\\eta_\\mathrm{diss}``, for the gas that form the stars.
      + `:ode_stellar_eta_i`         -> Photoionization efficiency, ``\\eta_\\mathrm{ion}``, for the gas that form the stars.
      + `:ode_stellar_r`             -> Mass recycling parameter, ``R``, for the gas that form the stars.
      + `:ode_stellar_cold_mf`       -> Cold gas mass fraction, for the gas that form the stars.
      + `:ode_stellar_gas_rho`       -> Gas mass density, for the gas that form the stars.
      + `:ode_stellar_gas_Z`         -> Gas metallicity, for the gas that form the stars (solar units).
      + `:ode_stellar_gas_mass`      -> Cell mass, for the gas that form the stars.
      + `:ode_stellar_gas_sfr`       -> SFR associated to the gas particles/cells within the code, for the gas that form the stars.
      + `:ode_stellar_gas_P`         -> Gas pressure, for the gas that form the stars.
  - `agg_function::Function=median`: The function used to aggregate the values of the quantity.

# Returns

  - The value of `quantity` for the whole system in `data_dict`.

# References

L. Blitz et al. (2006). *The Role of Pressure in GMC Formation II: The H2-Pressure Relation*. The Astrophysical Journal, **650(2)**, 933. [doi:10.1086/505417](https://doi.org/10.1086/505417)
"""
function integrateQty(
    data_dict::Dict,
    quantity::Symbol;
    agg_function::Function=median,
)::Number

    if quantity == :stellar_mass

        integrated_qty = sum(computeMass(data_dict, :stellar); init=0.0u"Msun")

    elseif quantity == :gas_mass

        integrated_qty = sum(computeMass(data_dict, :gas); init=0.0u"Msun")

    elseif quantity == :hydrogen_mass

        integrated_qty = sum(computeMass(data_dict, :hydrogen); init=0.0u"Msun")

    elseif quantity == :dm_mass

        integrated_qty = sum(computeMass(data_dict, :dark_matter); init=0.0u"Msun")

    elseif quantity == :bh_mass

        integrated_qty = sum(computeMass(data_dict, :black_hole); init=0.0u"Msun")

    elseif quantity == :molecular_mass

        integrated_qty = sum(computeMass(data_dict, :molecular); init=0.0u"Msun")

    elseif quantity == :br_molecular_mass

        integrated_qty = sum(computeMass(data_dict, :br_molecular); init=0.0u"Msun")

    elseif quantity == :atomic_mass

        integrated_qty = sum(computeMass(data_dict, :atomic); init=0.0u"Msun")

    elseif quantity == :ionized_mass

        integrated_qty = sum(computeMass(data_dict, :ionized); init=0.0u"Msun")

    elseif quantity == :neutral_mass

        integrated_qty = sum(computeMass(data_dict, :neutral); init=0.0u"Msun")

    elseif quantity == :stellar_gas_mass

        integrated_qty = sum(computeMass(data_dict, :stellar); init=0.0u"Msun")

    elseif quantity == :ode_metal_mass

        integrated_qty = sum(computeMass(data_dict, :metals); init=0.0u"Msun")

    elseif quantity == :ode_metallicity

        metal_mass = computeMass(data_dict, :metals)
        gas_mass = computeMass(data_dict, :gas)

        if isempty(metal_mass) || isempty(gas_mass)
            integrated_qty = NaN
        else
            integrated_qty = (sum(metal_mass) / sum(gas_mass)) / SOLAR_METALLICITY
        end

    elseif quantity == :dust_mass

        integrated_qty = sum(computeMass(data_dict, :dust); init=0.0u"Msun")

    elseif quantity == :stellar_number

        integrated_qty = length(data_dict[:stellar]["MASS"])

    elseif quantity == :gas_number

        integrated_qty = length(data_dict[:gas]["MASS"])

    elseif quantity == :dm_number

        integrated_qty = length(data_dict[:dark_matter]["MASS"])

    elseif quantity == :bh_number

        integrated_qty = length(data_dict[:black_hole]["MASS"])

    elseif quantity == :dust_stellar_mass_ratio

        Ms = sum(computeMass(data_dict, :stellar); init=0.0u"Msun")

        Md = sum(computeMass(data_dict, :dust); init=0.0u"Msun")

        if iszero(Ms)
            integrated_qty = NaN
        else
            integrated_qty = Md / Ms
        end

    elseif quantity == :molecular_fraction

        molecular_mass = sum(computeMass(data_dict, :molecular); init=0.0u"Msun")
        gas_mass = sum(computeMass(data_dict, :gas); init=0.0u"Msun")

        if iszero(gas_mass)
            integrated_qty = NaN
        else
            integrated_qty = molecular_mass / gas_mass
        end

    elseif quantity == :br_molecular_fraction

        molecular_mass = sum(computeMass(data_dict, :br_molecular); init=0.0u"Msun")
        gas_mass = sum(computeMass(data_dict, :gas); init=0.0u"Msun")

        if iszero(gas_mass)
            integrated_qty = NaN
        else
            integrated_qty = molecular_mass / gas_mass
        end

    elseif quantity == :atomic_fraction

        atomic_mass = sum(computeMass(data_dict, :atomic); init=0.0u"Msun")
        gas_mass = sum(computeMass(data_dict, :gas); init=0.0u"Msun")

        if iszero(gas_mass)
            integrated_qty = NaN
        else
            integrated_qty = atomic_mass / gas_mass
        end

    elseif quantity == :ionized_fraction

        ionized_mass = sum(computeMass(data_dict, :ionized); init=0.0u"Msun")
        gas_mass = sum(computeMass(data_dict, :gas); init=0.0u"Msun")

        if iszero(gas_mass)
            integrated_qty = NaN
        else
            integrated_qty = ionized_mass / gas_mass
        end

    elseif quantity == :neutral_fraction

        neutral_mass = sum(computeMass(data_dict, :neutral); init=0.0u"Msun")
        gas_mass = sum(computeMass(data_dict, :gas); init=0.0u"Msun")

        if iszero(gas_mass)
            integrated_qty = NaN
        else
            integrated_qty = neutral_mass / gas_mass
        end

    elseif quantity == :molecular_neutral_fraction

        molecular_mass = sum(computeMass(data_dict, :molecular); init=0.0u"Msun")
        neutral_mass = sum(computeMass(data_dict, :neutral); init=0.0u"Msun")

        if iszero(neutral_mass)
            integrated_qty = NaN
        else
            integrated_qty = molecular_mass / neutral_mass
        end

    elseif quantity == :ionized_neutral_fraction

        ionized_mass = sum(computeMass(data_dict, :ionized); init=0.0u"Msun")
        neutral_mass = sum(computeMass(data_dict, :atomic); init=0.0u"Msun")

        if iszero(neutral_mass)
            integrated_qty = NaN
        else
            integrated_qty = ionized_mass / neutral_mass
        end

    elseif quantity == :gas_mass_density

        density = computeVolumeDensity(data_dict, :gas)

        if isempty(density)
            integrated_qty = 0.0u"Msun * kpc^-3"
        else
            integrated_qty = agg_function(density)
        end

    elseif quantity == :stellar_gas_fraction

        stellar_gas_mass = sum(computeMass(data_dict, :stellar); init=0.0u"Msun")
        gas_mass = sum(computeMass(data_dict, :gas); init=0.0u"Msun")

        if iszero(gas_mass)
            integrated_qty = NaN
        else
            integrated_qty = stellar_gas_mass / gas_mass
        end

    elseif quantity == :metal_gas_fraction

        metal_mass = sum(computeMass(data_dict, :metals); init=0.0u"Msun")
        gas_mass = sum(computeMass(data_dict, :gas); init=0.0u"Msun")

        if iszero(gas_mass)
            integrated_qty = NaN
        else
            integrated_qty = metal_mass / gas_mass
        end

    elseif quantity == :dust_fraction

        dust_mass = sum(computeMass(data_dict, :dust); init=0.0u"Msun")
        gas_mass = sum(computeMass(data_dict, :gas); init=0.0u"Msun")

        if iszero(gas_mass)
            integrated_qty = NaN
        else
            integrated_qty = dust_mass / gas_mass
        end

    elseif quantity == :stellar_area_density

        integrated_qty = sum(computeMass(data_dict, :stellar); init=0.0u"Msun") / area(DISK_R)

    elseif quantity == :gas_area_density

        integrated_qty = sum(computeMass(data_dict, :gas); init=0.0u"Msun") / area(DISK_R)

    elseif quantity == :molecular_area_density

        integrated_qty = sum(computeMass(data_dict, :molecular); init=0.0u"Msun") / area(DISK_R)

    elseif quantity == :br_molecular_area_density

        integrated_qty = sum(computeMass(data_dict, :br_molecular); init=0.0u"Msun") / area(DISK_R)

    elseif quantity == :atomic_area_density

        integrated_qty = sum(computeMass(data_dict, :atomic); init=0.0u"Msun") / area(DISK_R)

    elseif quantity == :ionized_area_density

        integrated_qty = sum(computeMass(data_dict, :ionized); init=0.0u"Msun") / area(DISK_R)

    elseif quantity == :neutral_area_density

        integrated_qty = sum(computeMass(data_dict, :neutral); init=0.0u"Msun") / area(DISK_R)

    elseif quantity == :sfr_area_density

        sfr = sum(computeSFR(data_dict; age_resol=AGE_RESOLUTION); init=0.0u"Msun*yr^-1")

        integrated_qty = sfr / area(DISK_R)

    elseif quantity == :gas_td

        mass = computeMass(data_dict, :gas)
        td   = computeDepletionTime(mass, data_dict[:gas]["SFR "])

        if isempty(td)
            integrated_qty = NaN
        else
            integrated_qty = agg_function(td)
        end

    elseif quantity == :molecular_td

        mass = computeMass(data_dict, :molecular)
        td   = computeDepletionTime(mass, data_dict[:gas]["SFR "])

        if isempty(td)
            integrated_qty = NaN
        else
            integrated_qty = agg_function(td)
        end

    elseif quantity == :br_molecular_td

        mass = computeMass(data_dict, :br_molecular)
        td   = computeDepletionTime(mass, data_dict[:gas]["SFR "])

        if isempty(td)
            integrated_qty = NaN
        else
            integrated_qty = agg_function(td)
        end

    elseif quantity == :atomic_td

        mass = computeMass(data_dict, :atomic)
        td   = computeDepletionTime(mass, data_dict[:gas]["SFR "])

        if isempty(td)
            integrated_qty = NaN
        else
            integrated_qty = agg_function(td)
        end

    elseif quantity == :ionized_td

        mass = computeMass(data_dict, :ionized)
        td   = computeDepletionTime(mass, data_dict[:gas]["SFR "])

        if isempty(td)
            integrated_qty = NaN
        else
            integrated_qty = agg_function(td)
        end

    elseif quantity == :neutral_td

        mass = computeMass(data_dict, :neutral)
        td   = computeDepletionTime(mass, data_dict[:gas]["SFR "])

        if isempty(td)
            integrated_qty = NaN
        else
            integrated_qty = agg_function(td)
        end

    elseif quantity == :gas_metallicity

        metal_mass = sum(computeMetalMass(data_dict, :gas); init=0.0u"Msun")
        gas_mass = sum(computeMass(data_dict, :gas); init=0.0u"Msun")

        if iszero(gas_mass)
            integrated_qty = NaN
        else
            integrated_qty = (metal_mass / gas_mass) / SOLAR_METALLICITY
        end

    elseif quantity == :stellar_metallicity

        metal_mass = sum(computeMetalMass(data_dict, :stellar); init=0.0u"Msun")
        stellar_mass = sum(data_dict[:stellar]["MASS"]; init=0.0u"Msun")

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

        abundance = 12 + log10(computeGlobalAbundance(data_dict, :stellar, element_symbol))
        integrated_qty = isinf(abundance) ? NaN : abundance

    elseif quantity == :stellar_specific_am

        positions = data_dict[:stellar]["POS "]
        velocities = data_dict[:stellar]["VEL "]
        masses = data_dict[:stellar]["MASS"]

        if any(isempty, [positions, velocities, masses])
            integrated_qty = NaN
        else
            J = norm(computeTotalAngularMomentum(positions, velocities, masses; normal=false))
            integrated_qty = J / sum(masses)
        end

    elseif quantity == :gas_specific_am

        positions = data_dict[:gas]["POS "]
        velocities = data_dict[:gas]["VEL "]
        masses = computeMass(data_dict, :gas)

        if any(isempty, [positions, velocities, masses])
            integrated_qty = NaN
        else
            J = norm(computeTotalAngularMomentum(positions, velocities, masses; normal=false))
            integrated_qty = J / sum(masses)
        end

    elseif quantity == :dm_specific_am

        positions = data_dict[:dark_matter]["POS "]
        velocities = data_dict[:dark_matter]["VEL "]
        masses = data_dict[:dark_matter]["MASS"]

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

            integrated_qty = 0.0u"Msun * yr^-1"

        else

            # Get the physical times
            times = data_dict[:sim_data].snapshot_table[!, :physical_times]

            # Compute the time between snapshots
            Δt = times[present_idx] - times[present_idx - 1]

            integrated_qty = sum(computeSFR(data_dict; age_resol=Δt); init=0.0u"Msun * yr^-1")

        end

    elseif quantity == :ssfr

        # Get the global index (index in the context of the whole simulation) of the current snapshot
        present_idx = data_dict[:snap_data].global_index

        # Compute the total stellar mass
        stellar_mass = sum(data_dict[:stellar]["MASS"]; init=0.0u"Msun")

        if present_idx == 1 || iszero(stellar_mass)

            integrated_qty = 0.0u"yr^-1"

        else

            # Get the physical times
            times = data_dict[:sim_data].snapshot_table[!, :physical_times]
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
        stellar_mass = sum(data_dict[:stellar]["MASS"]; init=0.0u"Msun")

        if iszero(stellar_mass)
            integrated_qty = 0.0u"yr^-1"
        else
            integrated_qty = sfr / stellar_mass
        end

    elseif quantity == :stellar_eff

        ϵffs = computeEfficiencyFF(
            data_dict[:stellar]["RHOC"] .* u"mp",
            data_dict[:stellar]["GMAS"],
            data_dict[:stellar]["GSFR"],
        )

        # Filter out zeros and NaNs
        filter!(eff -> !isnan(eff) && !iszero(eff), ϵffs)

        if isempty(ϵffs)
            integrated_qty = NaN
        else
            integrated_qty = agg_function(ϵffs)
        end

    elseif quantity == :gas_eff

        mass = computeMass(data_dict, :gas)

        ϵffs = computeEfficiencyFF(data_dict[:gas]["RHO "], mass, data_dict[:gas]["SFR "])

        # Filter out zeros and NaNs
        filter!(eff -> !isnan(eff) && !iszero(eff), ϵffs)

        if isempty(ϵffs)
            integrated_qty = NaN
        else
            integrated_qty = agg_function(ϵffs)
        end

    elseif quantity == :molecular_eff

        mass = computeMass(data_dict, :molecular)

        ϵffs = computeEfficiencyFF(data_dict[:gas]["RHO "], mass, data_dict[:gas]["SFR "])

        # Filter out zeros and NaNs
        filter!(eff -> !isnan(eff) && !iszero(eff), ϵffs)

        if isempty(ϵffs)
            integrated_qty = NaN
        else
            integrated_qty = agg_function(ϵffs)
        end

    elseif quantity == :br_molecular_eff

        mass = computeMass(data_dict, :br_molecular)

        ϵffs = computeEfficiencyFF(data_dict[:gas]["RHO "], mass, data_dict[:gas]["SFR "])

        # Filter out zeros and NaNs
        filter!(eff -> !isnan(eff) && !iszero(eff), ϵffs)

        if isempty(ϵffs)
            integrated_qty = NaN
        else
            integrated_qty = agg_function(ϵffs)
        end

    elseif quantity == :atomic_eff

        mass = computeMass(data_dict, :atomic)

        ϵffs = computeEfficiencyFF(data_dict[:gas]["RHO "], mass, data_dict[:gas]["SFR "])

        # Filter out zeros and NaNs
        filter!(eff -> !isnan(eff) && !iszero(eff), ϵffs)

        if isempty(ϵffs)
            integrated_qty = NaN
        else
            integrated_qty = agg_function(ϵffs)
        end

    elseif quantity == :ionized_eff

        mass = computeMass(data_dict, :ionized)

        ϵffs = computeEfficiencyFF(data_dict[:gas]["RHO "], mass, data_dict[:gas]["SFR "])

        # Filter out zeros and NaNs
        filter!(eff -> !isnan(eff) && !iszero(eff), ϵffs)

        if isempty(ϵffs)
            integrated_qty = NaN
        else
            integrated_qty = agg_function(ϵffs)
        end

    elseif quantity == :neutral_eff

        mass = computeMass(data_dict, :neutral)

        ϵffs = computeEfficiencyFF(data_dict[:gas]["RHO "], mass, data_dict[:gas]["SFR "])

        # Filter out zeros and NaNs
        filter!(eff -> !isnan(eff) && !iszero(eff), ϵffs)

        if isempty(ϵffs)
            integrated_qty = NaN
        else
            integrated_qty = agg_function(ϵffs)
        end

    elseif quantity == :scale_factor

        integrated_qty = data_dict[:sim_data].snapshot_table[data_dict[:snap_data].global_index, :scale_factors]

    elseif quantity == :redshift

        integrated_qty = data_dict[:sim_data].snapshot_table[data_dict[:snap_data].global_index, :redshifts]

    elseif quantity == :physical_time

        integrated_qty = data_dict[:sim_data].snapshot_table[data_dict[:snap_data].global_index, :physical_times]

    elseif quantity == :lookback_time

        integrated_qty = data_dict[:sim_data].snapshot_table[data_dict[:snap_data].global_index, :lookback_times]

    elseif quantity == :ode_gas_it

        odit = data_dict[:gas]["ODIT"]
        filter!(!isnan, odit)

        if isempty(odit)
            integrated_qty = NaN
        else
            integrated_qty = agg_function(odit)
        end

    elseif quantity == :ode_gas_tau_s

        τS = data_dict[:gas]["TAUS"]
        filter!(!isnan, τS)

        if isempty(τS)
            integrated_qty = NaN
        else
            integrated_qty = agg_function(τS)
        end

    elseif quantity == :ode_gas_eta_d

        ηd = data_dict[:gas]["ETAD"]
        filter!(!isnan, ηd)

        if isempty(ηd)
            integrated_qty = NaN
        else
            integrated_qty = agg_function(ηd)
        end

    elseif quantity == :ode_gas_eta_i

        ηi = data_dict[:gas]["ETAI"]
        filter!(!isnan, ηi)

        if isempty(ηi)
            integrated_qty = NaN
        else
            integrated_qty = agg_function(ηi)
        end

    elseif quantity == :ode_gas_r

        R = data_dict[:gas]["PARR"]
        filter!(!isnan, R)

        if isempty(R)
            integrated_qty = NaN
        else
            integrated_qty = agg_function(R)
        end

    elseif quantity == :ode_gas_cold_mf

        cold_fraction = data_dict[:gas]["COLF"]
        gas_mass = data_dict[:gas]["MASS"]

        filter!(!isnan, cold_fraction)
        filter!(!isnan, gas_mass)

        if isempty(cold_fraction) || isempty(gas_mass)
            integrated_qty = NaN
        else
            cold_mass = cold_fraction .* gas_mass
            integrated_qty = sum(cold_mass) ./ sum(gas_mass)
        end

    elseif quantity == :ode_stellar_it

        odit = data_dict[:stellar]["ODIT"]
        filter!(!isnan, odit)

        if isempty(odit)
            integrated_qty = NaN
        else
            integrated_qty = agg_function(odit)
        end

    elseif quantity == :ode_stellar_tau_s

        τS = data_dict[:stellar]["TAUS"]
        filter!(!isnan, τS)

        if isempty(τS)
            integrated_qty = NaN
        else
            integrated_qty = agg_function(τS)
        end

    elseif quantity == :ode_stellar_eta_d

        ηd = data_dict[:stellar]["ETAD"]
        filter!(!isnan, ηd)

        if isempty(ηd)
            integrated_qty = NaN
        else
            integrated_qty = agg_function(ηd)
        end

    elseif quantity == :ode_stellar_eta_i

        ηi = data_dict[:stellar]["ETAI"]
        filter!(!isnan, ηi)

        if isempty(ηi)
            integrated_qty = NaN
        else
            integrated_qty = agg_function(ηi)
        end

    elseif quantity == :ode_stellar_r

        R = data_dict[:stellar]["PARR"]
        filter!(!isnan, R)

        if isempty(R)
            integrated_qty = NaN
        else
            integrated_qty = agg_function(R)
        end

    elseif quantity == :ode_stellar_cold_mf

        cold_fraction = data_dict[:stellar]["COLF"]
        gas_mass = data_dict[:stellar]["GMAS"]

        filter!(!isnan, cold_fraction)
        filter!(!isnan, gas_mass)

        if isempty(cold_fraction) || isempty(gas_mass)
            integrated_qty = NaN
        else
            cold_mass = cold_fraction .* gas_mass
            integrated_qty = sum(cold_mass) ./ sum(gas_mass)
        end

    elseif quantity == :ode_stellar_gas_rho

        ρ = data_dict[:stellar]["RHOC"]
        filter!(!isnan, ρ)

        if isempty(ρ)
            integrated_qty = NaN
        else
            integrated_qty = agg_function(ρ .* u"mp")
        end

    elseif quantity == :ode_stellar_gas_Z

        Z = data_dict[:stellar]["PARZ"]
        gas_mass = data_dict[:stellar]["GMAS"]

        filter!(!isnan, Z)
        filter!(!isnan, gas_mass)

        if isempty(Z) || isempty(gas_mass)
            integrated_qty = NaN
        else
            metal_mass = Z .* gas_mass
            integrated_qty = (sum(metal_mass) / sum(gas_mass)) / SOLAR_METALLICITY
        end

    elseif quantity == :ode_stellar_gas_mass

        gm = data_dict[:stellar]["GMAS"]
        filter!(!isnan, gm)

        if isempty(gm)
            integrated_qty = NaN
        else
            integrated_qty = agg_function(gm)
        end

    elseif quantity == :ode_stellar_gas_sfr

        gsfr = data_dict[:stellar]["GSFR"]
        filter!(!isnan, gsfr)

        if isempty(gsfr)
            integrated_qty = NaN
        else
            integrated_qty = agg_function(gsfr)
        end

    elseif quantity == :ode_stellar_gas_P

        P = data_dict[:stellar]["GPRE"]
        filter!(!isnan, P)

        if isempty(P)
            integrated_qty = NaN
        else
            integrated_qty = agg_function(P)
        end

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
      + `:stellar_gas_mass`            -> Stellar gas mass (according to our SF model).
      + `:ode_metal_mass`              -> Metal mass (according to our SF model).
      + `:ode_metallicity`             -> Metallicity (according to our SF model).
      + `:dust_mass`                   -> Dust mass.
      + `:molecular_fraction`          -> Gas mass fraction of molecular hydrogen.
      + `:br_molecular_fraction`       -> Gas mass fraction of molecular hydrogen, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_fraction`             -> Gas mass fraction of atomic hydrogen.
      + `:ionized_fraction`            -> Gas mass fraction of ionized hydrogen.
      + `:neutral_fraction`            -> Gas mass fraction of neutral hydrogen.
      + `:molecular_neutral_fraction`  -> Fraction of molecular hydrogen in the neutral gas.
      + `:ionized_neutral_fraction`    -> Fraction of ionized gas to neutral gas.
      + `:stellar_gas_fraction`        -> Stellar gas fraction (according to our SF model).
      + `:metal_gas_fraction`          -> Metallicity (according to our SF model).
      + `:dust_fraction`               -> Dust mass fraction.
      + `:gas_mass_density`            -> Gas mass density.
      + `:hydrogen_mass_density`       -> Hydrogen mass density.
      + `:gas_number_density`          -> Gas number density.
      + `:molecular_number_density`    -> Molecular hydrogen number density.
      + `:br_molecular_number_density` -> Molecular hydrogen number density, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_number_density`       -> Atomic hydrogen number density.
      + `:ionized_number_density`      -> Ionized hydrogen number density.
      + `:neutral_number_density`      -> Neutral hydrogen number density.
      + `:gas_td`                      -> Total gas depletion time.
      + `:molecular_td`                -> Molecular hydrogen (``\\mathrm{H_2}``) depletion time.
      + `:br_molecular_td`             -> Molecular hydrogen (``\\mathrm{H_2}``) depletion time, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_td`                   -> Atomic hydrogen (``\\mathrm{HI}``) depletion time.
      + `:ionized_td`                  -> Ionized hydrogen (``\\mathrm{HII}``) depletion time.
      + `:neutral_td`                  -> Neutral hydrogen (``\\mathrm{HI + H_2}``) depletion time.
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
      + `:sfr`                         -> Star formation rate.
      + `:ssfr`                        -> The specific star formation rate.
      + `:observational_sfr`           -> Star formation rate of the last `AGE_RESOLUTION`.
      + `:observational_ssfr`          -> The specific star formation rate of the last `AGE_RESOLUTION`.
      + `:stellar_eff`                 -> Star formation efficiency per free-fall time for the gas that has turn into stars.
      + `:gas_eff`                     -> Star formation efficiency per free-fall time for the gas.
      + `:molecular_eff`               -> Star formation efficiency per free-fall time for the molecular hydrogen (``\\mathrm{H_2}``) gas.
      + `:br_molecular_eff`            -> Star formation efficiency per free-fall time for the molecular hydrogen (``\\mathrm{H_2}``) gas, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_eff`                  -> Star formation efficiency per free-fall time for the atomic hydrogen (``\\mathrm{HI}``) gas.
      + `:ionized_eff`                 -> Star formation efficiency per free-fall time for the ionized hydrogen (``\\mathrm{HII}``) gas.
      + `:neutral_eff`                 -> Star formation efficiency per free-fall time for the neutral hydrogen (``\\mathrm{HI + H_2}``) gas.
      + `:temperature`                 -> Gas temperature, as ``\\log_{10}(T \\, / \\, \\mathrm{K})``.
      + `:pressure`                    -> Gas pressure.
      + `:ode_gas_it`                  -> Integration time.
      + `:ode_gas_tau_s`               -> Star formation time scale, ``\\tau_\\mathrm{S}``.
      + `:ode_gas_eta_d`               -> Photodissociation efficiency, ``\\eta_\\mathrm{diss}``.
      + `:ode_gas_eta_i`               -> Photoionization efficiency, ``\\eta_\\mathrm{ion}``.
      + `:ode_gas_r`                   -> Mass recycling parameter, ``R``.
      + `:ode_gas_cold_mf`             -> Cold gas mass fraction.
      + `:ode_stellar_it`              -> Integration time, for the gas that form the stars.
      + `:ode_stellar_tau_s`           -> Star formation time scale, ``\\tau_\\mathrm{S}``, for the gas that form the stars.
      + `:ode_stellar_eta_d`           -> Photodissociation efficiency, ``\\eta_\\mathrm{diss}``, for the gas that form the stars.
      + `:ode_stellar_eta_i`           -> Photoionization efficiency, ``\\eta_\\mathrm{ion}``, for the gas that form the stars.
      + `:ode_stellar_r`               -> Mass recycling parameter, ``R``, for the gas that form the stars.
      + `:ode_stellar_cold_mf`         -> Cold gas mass fraction, for the gas that form the stars.
      + `:ode_stellar_gas_rho`         -> Gas mass density, for the gas that form the stars.
      + `:ode_stellar_gas_Z`           -> Gas metallicity, for the gas that form the stars (solar units).
      + `:ode_stellar_gas_mass`        -> Cell mass, for the gas that form the stars.
      + `:ode_stellar_gas_sfr`         -> SFR associated to the gas particles/cells within the code, for the gas that form the stars.
      + `:ode_stellar_gas_P`           -> Gas pressure, for the gas that form the stars.

# Returns

  - The values of `quantity` for every cell/particle.

# References

L. Blitz et al. (2006). *The Role of Pressure in GMC Formation II: The H2-Pressure Relation*. The Astrophysical Journal, **650(2)**, 933. [doi:10.1086/505417](https://doi.org/10.1086/505417)
"""
function scatterQty(data_dict::Dict, quantity::Symbol)::Vector{<:Number}

    if quantity == :stellar_mass

        scatter_qty = computeMass(data_dict, :stellar)

    elseif quantity == :gas_mass

        scatter_qty = computeMass(data_dict, :gas)

    elseif quantity == :hydrogen_mass

        scatter_qty = computeMass(data_dict, :hydrogen)

    elseif quantity == :dm_mass

        scatter_qty = computeMass(data_dict, :dark_matter)

    elseif quantity == :bh_mass

        scatter_qty = computeMass(data_dict, :black_hole)

    elseif quantity == :molecular_mass

        scatter_qty = computeMass(data_dict, :molecular)

    elseif quantity == :br_molecular_mass

        scatter_qty = computeMass(data_dict, :br_molecular)

    elseif quantity == :atomic_mass

        scatter_qty = computeMass(data_dict, :atomic)

    elseif quantity == :ionized_mass

        scatter_qty = computeMass(data_dict, :ionized)

    elseif quantity == :neutral_mass

        scatter_qty = computeMass(data_dict, :neutral)

    elseif quantity == :stellar_gas_mass

        scatter_qty = computeMass(data_dict, :stellar)

    elseif quantity == :ode_metal_mass

        scatter_qty = computeMass(data_dict, :metals)

    elseif quantity == :ode_metallicity

        metal_mass = computeMass(data_dict, :metals)
        gas_mass = computeMass(data_dict, :gas)

        if isempty(metal_mass) || isempty(gas_mass)
            scatter_qty = []
        else
            scatter_qty = (metal_mass ./ gas_mass) ./ SOLAR_METALLICITY
        end

    elseif quantity == :dust_mass

        scatter_qty = computeMass(data_dict, :dust)

    elseif quantity == :molecular_fraction

        scatter_qty = computeFraction(data_dict, :molecular)

    elseif quantity == :br_molecular_fraction

        scatter_qty = computeFraction(data_dict, :br_molecular)

    elseif quantity == :atomic_fraction

        scatter_qty = computeFraction(data_dict, :atomic)

    elseif quantity == :ionized_fraction

        scatter_qty = computeFraction(data_dict, :ionized)

    elseif quantity == :neutral_fraction

        scatter_qty = computeFraction(data_dict, :neutral)

    elseif quantity == :molecular_neutral_fraction

        fm = computeFraction(data_dict, :molecular)
        fn = computeFraction(data_dict, :neutral)

        scatter_qty = fm ./ fn

    elseif quantity == :ionized_neutral_fraction

        fi = computeFraction(data_dict, :ionized)
        fn = computeFraction(data_dict, :neutral)

        scatter_qty = similar(fi, Float64)

        for i in eachindex(fn)
            if iszero(fn[i])
                scatter_qty[i] = NaN
            else
                scatter_qty[i] = fi[i] / fn[i]
            end
        end

    elseif quantity == :stellar_gas_fraction

        scatter_qty = computeFraction(data_dict, :stellar)

    elseif quantity == :metal_gas_fraction

        scatter_qty = computeFraction(data_dict, :metals)

    elseif quantity == :dust_fraction

        scatter_qty = computeFraction(data_dict, :dust)

    elseif quantity == :gas_mass_density

        scatter_qty = computeVolumeDensity(data_dict, :gas)

    elseif quantity == :hydrogen_mass_density

        scatter_qty = computeVolumeDensity(data_dict, :hydrogen)

    elseif quantity == :gas_number_density

        scatter_qty = computeNumberDensity(data_dict, :gas)

    elseif quantity == :molecular_number_density

        scatter_qty = computeNumberDensity(data_dict, :molecular)

    elseif quantity == :br_molecular_number_density

        scatter_qty = computeNumberDensity(data_dict, :br_molecular)

    elseif quantity == :atomic_number_density

        scatter_qty = computeNumberDensity(data_dict, :atomic)

    elseif quantity == :ionized_number_density

        scatter_qty = computeNumberDensity(data_dict, :ionized)

    elseif quantity == :neutral_number_density

        scatter_qty = computeNumberDensity(data_dict, :neutral)

    elseif quantity == :gas_td

        scatter_qty = computeDepletionTime(computeMass(data_dict, :gas), data_dict[:gas]["SFR "])

    elseif quantity == :molecular_td

        scatter_qty = computeDepletionTime(
            computeMass(data_dict, :molecular),
            data_dict[:gas]["SFR "],
        )

    elseif quantity == :br_molecular_td

        scatter_qty = computeDepletionTime(
            computeMass(data_dict, :br_molecular),
            data_dict[:gas]["SFR "],
        )

    elseif quantity == :atomic_td

        scatter_qty = computeDepletionTime(
            computeMass(data_dict, :atomic),
            data_dict[:gas]["SFR "],
        )

    elseif quantity == :ionized_td

        scatter_qty = computeDepletionTime(
            computeMass(data_dict, :ionized),
            data_dict[:gas]["SFR "],
        )

    elseif quantity == :neutral_td

        scatter_qty = computeDepletionTime(
            computeMass(data_dict, :neutral),
            data_dict[:gas]["SFR "],
        )

    elseif quantity == :gas_metallicity

        scatter_qty = setPositive(data_dict[:gas]["GZ  "]) ./ SOLAR_METALLICITY

    elseif quantity == :stellar_metallicity

        scatter_qty = setPositive(data_dict[:stellar]["GZ2 "]) ./ SOLAR_METALLICITY

    elseif quantity ∈ GAS_ABUNDANCE

        element_symbol = Symbol(first(split(string(quantity), "_")))

        abundances = computeAbundance(
            data_dict,
            :gas,
            element_symbol;
            solar=false,
        )

        if isempty(abundances)
            scatter_qty = Float64[]
        else
            scatter_qty = 12 .+ log10.(abundances)
            replace!(x -> isinf(x) ? NaN : x, scatter_qty)
        end

    elseif quantity ∈ STELLAR_ABUNDANCE

        element_symbol = Symbol(first(split(string(quantity), "_")))

        abundances = computeAbundance(
            data_dict,
            :stellar,
            element_symbol;
            solar=false,
        )

        if isempty(abundances)
            scatter_qty = Float64[]
        else
            scatter_qty = 12 .+ log10.(abundances)
            replace!(x -> isinf(x) ? NaN : x, scatter_qty)
        end

    elseif quantity == :stellar_radial_distance

        scatter_qty = computeDistance(data_dict[:stellar]["POS "])

    elseif quantity == :gas_radial_distance

        scatter_qty = computeDistance(data_dict[:gas]["POS "])

    elseif quantity == :dm_radial_distance

        scatter_qty = computeDistance(data_dict[:dark_matter]["POS "])

    elseif quantity == :stellar_xy_distance

        if isempty(data_dict[:stellar]["POS "])
            scatter_qty = eltype(data_dict[:stellar]["POS "])[]
        else
            scatter_qty = computeDistance(data_dict[:stellar]["POS "][1:2, :])
        end

    elseif quantity == :gas_xy_distance

        if isempty(data_dict[:gas]["POS "])
            scatter_qty = eltype(data_dict[:gas]["POS "])[]
        else
            scatter_qty = computeDistance(data_dict[:gas]["POS "][1:2, :])
        end

    elseif quantity == :dm_xy_distance

        if isempty(data_dict[:dark_matter]["POS "])
            scatter_qty = eltype(data_dict[:dark_matter]["POS "])[]
        else
            scatter_qty = computeDistance(data_dict[:dark_matter]["POS "][1:2, :])
        end

    elseif quantity == :gas_sfr

        scatter_qty = data_dict[:gas]["SFR "]

    elseif quantity == :stellar_circularity

        (
            !logging[] ||
            @info("scatterQty: The stellar circularity depends on the positions and velocities of \
            all cells/particles. So, after filtering, the result for a given star will change")
        )

        scatter_qty = computeCircularity(data_dict, :stellar)

    elseif quantity == :stellar_vcirc

       (
            !logging[] ||
            @info("scatterQty: The stellar circular velocity depends on the positions and \
            velocities of all cells/particles. So, after filtering, the result for a given star \
            will change")
       )

        _, scatter_qty = computeVcirc(data_dict, :stellar)

    elseif quantity == :stellar_vradial

        scatter_qty = computeVpolar(data_dict, :stellar, :radial)

    elseif quantity == :stellar_vtangential

        scatter_qty = computeVpolar(data_dict, :stellar, :tangential)

    elseif quantity == :stellar_vzstar

        scatter_qty = computeVpolar(data_dict, :stellar, :zstar)

    elseif quantity == :stellar_age

        scatter_qty = computeStellarAge(data_dict)

    elseif quantity == :sfr

        # Get the global index (index in the context of the whole simulation) of the current snapshot
        present_idx = data_dict[:snap_data].global_index

        if present_idx == 1

            scatter_qty = zeros(typeof(1.0u"Msun*yr^-1"), length(data_dict[:stellar]["MASS"]))

        else

            # Get the physical times
            times = data_dict[:sim_data].snapshot_table[!, :physical_times]
            # Compute the time between snapshots
            Δt = times[present_idx] - times[present_idx - 1]

            scatter_qty = computeSFR(data_dict; age_resol=Δt)

        end

    elseif quantity == :ssfr

        # Get the global index (index in the context of the whole simulation) of the current snapshot
        present_idx = data_dict[:snap_data].global_index

        # Load the stellar masses
        stellar_masses = data_dict[:stellar]["MASS"]

        if present_idx == 1 || isempty(stellar_masses)

            scatter_qty = Float64[]

        else

            # Get the physical times
            times = data_dict[:sim_data].snapshot_table[!, :physical_times]
            # Compute the time between snapshots
            Δt = times[present_idx] - times[present_idx - 1]

            scatter_qty = computeSFR(data_dict; age_resol=Δt) ./ stellar_masses

        end

    elseif quantity == :observational_sfr

        scatter_qty = computeSFR(data_dict; age_resol=AGE_RESOLUTION)

    elseif quantity == :observational_ssfr

        sfr = computeSFR(data_dict; age_resol=AGE_RESOLUTION)
        stellar_masses = data_dict[:stellar]["MASS"]

        scatter_qty = sfr ./ stellar_masses

    elseif quantity == :stellar_eff

        scatter_qty = computeEfficiencyFF(
            data_dict[:stellar]["RHOC"] .* u"mp",
            data_dict[:stellar]["GMAS"],
            data_dict[:stellar]["GSFR"],
        )

    elseif quantity == :gas_eff

        mass = computeMass(data_dict, :gas)

        scatter_qty = computeEfficiencyFF(data_dict[:gas]["RHO "], mass, data_dict[:gas]["SFR "])

    elseif quantity == :molecular_eff

        mass = computeMass(data_dict, :molecular)

        scatter_qty = computeEfficiencyFF(data_dict[:gas]["RHO "], mass, data_dict[:gas]["SFR "])

    elseif quantity == :br_molecular_eff

        mass = computeMass(data_dict, :br_molecular)

        scatter_qty = computeEfficiencyFF(data_dict[:gas]["RHO "], mass, data_dict[:gas]["SFR "])

    elseif quantity == :atomic_eff

        mass = computeMass(data_dict, :atomic)

        scatter_qty = computeEfficiencyFF(data_dict[:gas]["RHO "], mass, data_dict[:gas]["SFR "])

    elseif quantity == :ionized_eff

        mass = computeMass(data_dict, :ionized)

        scatter_qty = computeEfficiencyFF(data_dict[:gas]["RHO "], mass, data_dict[:gas]["SFR "])

    elseif quantity == :neutral_eff

        mass = computeMass(data_dict, :neutral)

        scatter_qty = computeEfficiencyFF(data_dict[:gas]["RHO "], mass, data_dict[:gas]["SFR "])

    elseif quantity == :temperature

        scatter_qty = log10.(ustrip.(u"K", data_dict[:gas]["TEMP"]))
        replace!(x -> isinf(x) ? NaN : x, scatter_qty)

    elseif quantity == :pressure

        scatter_qty = data_dict[:gas]["PRES"]

    elseif quantity == :ode_gas_it

        scatter_qty = data_dict[:gas]["ODIT"]

    elseif quantity == :ode_gas_tau_s

        scatter_qty = data_dict[:gas]["TAUS"]

    elseif quantity == :ode_gas_eta_d

        scatter_qty = data_dict[:gas]["ETAD"]

    elseif quantity == :ode_gas_eta_i

        scatter_qty = data_dict[:gas]["ETAI"]

    elseif quantity == :ode_gas_r

        scatter_qty = data_dict[:gas]["PARR"]

    elseif quantity == :ode_gas_cold_mf

        scatter_qty = data_dict[:gas]["COLF"]

    elseif quantity == :ode_stellar_it

        scatter_qty = data_dict[:stellar]["ODIT"]

    elseif quantity == :ode_stellar_tau_s

        scatter_qty = data_dict[:stellar]["TAUS"]

    elseif quantity == :ode_stellar_eta_d

        scatter_qty = data_dict[:stellar]["ETAD"]

    elseif quantity == :ode_stellar_eta_i

        scatter_qty = data_dict[:stellar]["ETAI"]

    elseif quantity == :ode_stellar_r

        scatter_qty = data_dict[:stellar]["PARR"]

    elseif quantity == :ode_stellar_cold_mf

        scatter_qty = data_dict[:stellar]["COLF"]

    elseif quantity == :ode_stellar_gas_rho

        scatter_qty = data_dict[:stellar]["RHOC"] .* u"mp"

    elseif quantity == :ode_stellar_gas_Z

        scatter_qty = data_dict[:stellar]["PARZ"] ./ SOLAR_METALLICITY

    elseif quantity == :ode_stellar_gas_mass

        scatter_qty = data_dict[:stellar]["GMAS"]

    elseif quantity == :ode_stellar_gas_sfr

        scatter_qty = data_dict[:stellar]["GSFR"]

    elseif quantity == :ode_stellar_gas_P

        scatter_qty = data_dict[:stellar]["GPRE"]

    else

        throw(ArgumentError("scatterQty: I don't recognize the quantity :$(quantity)"))

    end

    if isempty(scatter_qty)
        return Number[]
    end

    return scatter_qty

end
