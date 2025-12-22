####################################################################################################
# Aggregator functions
####################################################################################################

"""
    scatterQty(data_dict::Dict, quantity::Symbol)::Vector{<:Number}

Compute `quantity` for each cell/particle in `data_dict`.

# Arguments

  - `data_dict::Dict`: Data dictionary (see [`makeDataDict`](@ref) for the canonical description).
  - `quantity::Symbol`: Target quantity. See [`plotParams`](@ref) for possibilities; only quantities well defined for each cell/particle individually are possible.

# Returns

  - The values of `quantity` for every cell/particle.
"""
function scatterQty(data_dict::Dict, quantity::Symbol)::Vector{<:Number}

    #####################
    # Derived quantities
    #####################

    if quantity ∈ DERIVED_QTY

        magnitude, component = QUANTITY_SPLITS[quantity]

        if magnitude == :mass

            scatter_qty = computeMass(data_dict, component)

        elseif magnitude == :mass_density

            scatter_qty = computeMassDensity(data_dict, component)

        elseif magnitude == :number_density

            scatter_qty = computeNumberDensity(data_dict, component)

        elseif magnitude == :number

            scatter_qty = computeNumber(data_dict, component)

        elseif magnitude == :fraction

            scatter_qty = computeFraction(data_dict, component)

        elseif magnitude == :eff

            scatter_qty = computeEfficiencyFF(data_dict, component)

        elseif magnitude == :specific_z_angular_momentum

            scatter_qty = computeSpecificAngularMomentum(data_dict, component)

        elseif magnitude == :z_angular_momentum

            scatter_qty = computeAngularMomentum(data_dict, component)

        elseif magnitude == :circularity

            (
                logging[] &&
                @info("scatterQty: The circularity depends on the positions and velocities of \
                all cells/particles. So, after filtering, the result for a given cell/particle \
                will change")
            )

            scatter_qty = computeCircularity(data_dict, component)

        elseif magnitude == :circular_velocity

            (
                logging[] &&
                @info("scatterQty: The circular velocity depends on the positions and velocities \
                of all cells/particles. So, after filtering, the result for a cell/particle \
                will change")
            )

            scatter_qty = computeVcirc(data_dict, component)[2]

        elseif magnitude == :radial_velocity

            scatter_qty = computeVpolar(data_dict, component, :radial)

        elseif magnitude == :tangential_velocity

            scatter_qty = computeVpolar(data_dict, component, :tangential)

        elseif magnitude == :zstar_velocity

            scatter_qty = computeVpolar(data_dict, component, :zstar)

        elseif magnitude == :kinetic_energy

            scatter_qty = computeKineticEnergy(data_dict, component)

        elseif magnitude == :potential_energy

            scatter_qty = computePotentialEnergy(data_dict, component)

        elseif magnitude == :total_energy

            scatter_qty = computeTotalEnergy(data_dict, component)

        elseif magnitude == :depletion_time

            scatter_qty = computeDepletionTime(data_dict, component)

        elseif magnitude == :xy_distance

            scatter_qty = computeXYDistance(data_dict, component)

        elseif magnitude == :radial_distance

            scatter_qty = computeRadialDistance(data_dict, component)

        else

            throw(ArgumentError("scatterQty: magnitude :$(magnitude) is not well defined for each \
            cell/particle!"))

        end

    elseif quantity ∈ SFM_STELLAR_QTY

        magnitude, _ = SFM_QTY_SPLITS[quantity]

        scatter_qty = data_dict[:stellar][SFM_KEYS[magnitude]]

    elseif quantity ∈ SFM_GAS_QTY

        magnitude, _ = SFM_QTY_SPLITS[quantity]

        scatter_qty = data_dict[:gas][SFM_KEYS[magnitude]]

    elseif quantity ∈ SFM_DERIVED_QTY

        scatter_qty = computeDerivedSFMQty(data_dict, quantity)

    #########
    # Ratios
    #########

    elseif quantity ∈ RATIO_QTY

        qty_01, qty_02 = RATIO_SPLITS[quantity]

        scatter_qty_01 = scatterQty(data_dict, qty_01)
        scatter_qty_02 = scatterQty(data_dict, qty_02)

        scatter_qty = similar(scatter_qty_01, Float64)

        Threads.@threads for i in eachindex(scatter_qty_01)

            if iszero(scatter_qty_02[i])
                scatter_qty[i] = NaN
            else
                scatter_qty[i] = scatter_qty_01[i] / scatter_qty_02[i]
            end

        end

    #################
    # Gas quantities
    #################

    elseif quantity == :temperature

        scatter_qty = data_dict[:gas]["TEMP"]

    elseif quantity == :pressure

        scatter_qty = data_dict[:gas]["PRES"]

    #################
    # SFR quantities
    #################

    elseif quantity == :sfr

        # Read the global index (index in the context of the whole simulation) of the current snapshot
        present_idx = data_dict[:snap_data].global_index

        if present_idx == 1

            scatter_qty = zeros(typeof(1.0u"Msun * yr^-1"), length(data_dict[:stellar]["MASS"]))

        else

            # Read the physical times
            times = data_dict[:sim_data].snapshot_table[!, :physical_times]

            # Compute the time between snapshots
            Δt = times[present_idx] - times[present_idx - 1]

            # For stellar particles younger than Δt, the SFR is its mass divided by Δt
            # For older particles it is 0
            scatter_qty = computeSFR(data_dict; age_limit=Δt)

        end

    elseif quantity == :ssfr

        # Read the global index (index in the context of the whole simulation) of the current snapshot
        present_idx = data_dict[:snap_data].global_index

        if present_idx == 1

            scatter_qty = zeros(typeof(1.0u"yr^-1"), length(data_dict[:stellar]["GAGE"]))

        else

            # Read the physical times
            times = data_dict[:sim_data].snapshot_table[!, :physical_times]

            # Compute the time between snapshots
            Δt = times[present_idx] - times[present_idx - 1]

            scatter_qty = computeSSFR(data_dict; age_limit=Δt)

        end

    elseif quantity == :observational_sfr

        scatter_qty = computeSFR(data_dict; age_limit=AGE_RESOLUTION)

    elseif quantity == :observational_ssfr

        scatter_qty = computeSSFR(data_dict; age_limit=AGE_RESOLUTION)

    elseif quantity == :stellar_age

        scatter_qty = computeStellarAge(data_dict)

    elseif quantity == :stellar_birth_time

        scatter_qty = computeStellarBirthTime(data_dict)

    elseif quantity == :gas_sfr

        scatter_qty = data_dict[:gas]["SFR "]

    #########################
    # Metallicity quantities
    #########################

    elseif quantity == :gas_metallicity

        scatter_qty = computeFraction(data_dict, :Z_gas) ./ SOLAR_METALLICITY

    elseif quantity == :stellar_metallicity

        scatter_qty = computeFraction(data_dict, :Z_stellar) ./ SOLAR_METALLICITY

    elseif quantity == :ode_metallicity

        scatter_qty = computeFraction(data_dict, :ode_metals) ./ SOLAR_METALLICITY

    elseif quantity ∈ GAS_ABUNDANCE

        element = GAS_ABUNDANCE_SPLITS[quantity]

        abundances = computeAbundance(data_dict, :gas, element; solar=false)

        scatter_qty = ABUNDANCE_SHIFT[element] .+ log10.(abundances)

        replace!(x -> isinf(x) ? NaN : x, scatter_qty)

    elseif quantity ∈ STELLAR_ABUNDANCE

        element = STELLAR_ABUNDANCE_SPLITS[quantity]

        abundances = computeAbundance(data_dict, :stellar, element; solar=false)

        scatter_qty = ABUNDANCE_SHIFT[element] .+ log10.(abundances)

        replace!(x -> isinf(x) ? NaN : x, scatter_qty)

    else

        throw(ArgumentError("scatterQty: I don't recognize the quantity :$(quantity). \
        It may not be well defined for each cell/particle!"))

    end

    return isempty(scatter_qty) ? Number[] : copy(scatter_qty)

end

"""
    integrateQty(data_dict::Dict, quantity::Symbol; <keyword arguments>)::Number

Compute `quantity` for the whole system of cell/particles in `data_dict`.

# Arguments

  - `data_dict::Dict`: Data dictionary (see [`makeDataDict`](@ref) for the canonical description).
  - `quantity::Symbol`: Target quantity. See [`plotParams`](@ref) for possibilities.
  - `agg_function::Union{Function,Symbol}=:default`: If `quantity` is one the the listed symbols in [`DERIVED_QTY`](@ref), [`SFM_STELLAR_QTY`](@ref) or [`SFM_GAS_QTY`](@ref), you can pass an `agg_function` to accumulate the values given by [`scatterQty`](@ref). If `agg_function` is left as `:default` [`integrateQty`](@ref) will try to compute the most reasonable global value for `quantity`.

# Returns

  - The value of `quantity` for the whole system of cell/particles in `data_dict`.
"""
function integrateQty(
    data_dict::Dict,
    quantity::Symbol;
    agg_function::Union{Function,Symbol}=:default,
)::Number

    if agg_function == :default

        #####################
        # Derived quantities
        #####################

        if quantity ∈ DERIVED_QTY

            magnitude, component = QUANTITY_SPLITS[quantity]

            if magnitude == :mass

                integrated_qty = sum(computeMass(data_dict, component); init=0.0u"Msun")

            elseif magnitude == :number

                if component ∈ [:stellar, :dark_matter, :black_hole, :gas]

                    integrated_qty = lenght(data_dict[component]["MASS"])

                elseif component == :Z_stellar

                    integrated_qty = lenght(data_dict[:stellar]["MASS"])

                else

                    integrated_qty = lenght(data_dict[:gas]["MASS"])

                end

            elseif magnitude == :fraction

                if component == :Z_stellar
                    type = :stellar
                else
                    type = :gas
                end

                comp_mass = sum(computeMass(data_dict, component); init=0.0u"Msun")
                ref_mass  = sum(computeMass(data_dict, type); init=0.0u"Msun")

                if iszero(ref_mass)
                    integrated_qty = NaN
                else
                    integrated_qty = comp_mass / ref_mass
                end

            elseif magnitude == :clumping_factor

                integrated_qty = computeClumpingFactor(data_dict, component)

            elseif magnitude == :specific_z_angular_momentum

                Lz = sum(
                    computeAngularMomentum(data_dict, component);
                    init=0.0u"Msun * pc^2 * yr^-1",
                )

                M = sum(computeMass(data_dict, component); init=0.0u"Msun")

                integrated_qty = Lz / M

            elseif magnitude == :z_angular_momentum

                integrated_qty = sum(
                    computeAngularMomentum(data_dict, component);
                    init=0.0u"Msun * pc^2 * yr^-1",
                )

            elseif magnitude == :spin_parameter

                integrated_qty = computeSpinParameter(data_dict, component; R=DISK_R)

            elseif magnitude == :potential_energy

                integrated_qty = sum(computePotentialEnergy(data_dict, component); init=0.0u"erg")

            elseif magnitude == :depletion_time

                M   = sum(computeMass(data_dict, component); init=0.0u"Msun")
                sfr = sum(scatterQty(data_dict, :sfr); init=0.0u"Msun * yr^-1")

                integrated_qty = M / sfr

            else

                throw(ArgumentError("integrateQty: I don't know how to compute a global value for \
                :$(quantity). Try passing a personalized `agg_function`"))

            end

        elseif quantity ∈ SFM_QTY

            magnitude, component = SFM_QTY_SPLITS[quantity]

            if component == :ode_gas_

                type     = :gas
                gas_mass = data_dict[:gas]["MASS"]

            elseif component == :ode_stellar_

                type     = :stellar
                gas_mass = data_dict[:stellar]["GMAS"]

            else

                throw(ArgumentError("integrateQty: I don't recognize the component :$(component)"))

            end

            filter!(!isnan, gas_mass)
            Mg = sum(gas_mass; init=0.0u"Msun")

            if magnitude == :cold_mass_frac

                cold_fractions = data_dict[type]["COLF"]
                cold_mass      = cold_fractions .* gas_mass

                filter!(!isnan, cold_mass)

                Mc = sum(cold_mass; init=0.0u"Msun")

                if iszero(Mg)
                    integrated_qty = NaN
                else
                    integrated_qty = Mc / Mg
                end

            elseif magnitude == :parameter_metallicity

                metallicity = data_dict[type]["PARZ"]
                metals_mass = metallicity .* gas_mass

                filter!(!isnan, metals_mass)

                MZ = sum(metals_mass; init=0.0u"Msun")

                if iszero(Mg)
                    integrated_qty = NaN
                else
                    integrated_qty = (MZ / Mg) / SOLAR_METALLICITY
                end

            elseif magnitude == :gas_mass

                integrated_qty = Mg

            else

                throw(ArgumentError("integrateQty: I don't know how to compute a global value for \
                :$(quantity). Try passing a personalized `agg_function`"))

            end

        elseif quantity ∈ SFM_DERIVED_QTY

            throw(ArgumentError("integrateQty: I don't know how to compute a global value for \
            :$(quantity). Try passing a personalized `agg_function`"))

        #########
        # Ratios
        #########

        elseif quantity ∈ RATIO_QTY

            qty_01, qty_02 = RATIO_SPLITS[quantity]

            integrated_qty_01 = integrateQty(data_dict, qty_01)
            integrated_qty_02 = integrateQty(data_dict, qty_02)

            if iszero(integrated_qty_02)
                integrated_qty = NaN
            else
                integrated_qty = integrated_qty_01 / integrated_qty_02
            end

        ###########################
        # Group catalog quantities
        ###########################

        elseif !isnothing(match(r"^(.*)_(\d+)$", string(quantity)))

            magnitude, halo_idx = parseHaloQuantity(quantity)

            integrated_qty = data_dict[:group][HALO_KEYS[magnitude]][halo_idx]

        #################
        # SFR quantities
        #################

        elseif quantity ∈ [:sfr, :observational_sfr]

            integrated_qty = sum(scatterQty(data_dict, quantity); init=0.0u"Msun * yr^-1")

        elseif quantity == :ssfr

            # Read the global index (index in the context of the whole simulation) of the current snapshot
            present_idx = data_dict[:snap_data].global_index

            if present_idx == 1

                integrated_qty = 0.0u"yr^-1"

            else

                # Get the physical times
                times = data_dict[:sim_data].snapshot_table[!, :physical_times]

                # Compute the time between snapshots
                Δt = times[present_idx] - times[present_idx - 1]

                integrated_qty = 1.0 / Δt

            end

        elseif quantity == :observational_ssfr

            integrated_qty = 1.0 / AGE_RESOLUTION

        ##################
        # Time quantities
        ##################

        elseif quantity == :scale_factor

            snap_idx = data_dict[:snap_data].global_index

            integrated_qty = data_dict[:sim_data].snapshot_table[snap_idx, :scale_factors]

        elseif quantity == :redshift

            snap_idx = data_dict[:snap_data].global_index

            integrated_qty = data_dict[:sim_data].snapshot_table[snap_idx, :redshifts]

        elseif quantity == :physical_time

            snap_idx = data_dict[:snap_data].global_index

            integrated_qty = data_dict[:sim_data].snapshot_table[snap_idx, :physical_times]

        elseif quantity == :lookback_time

            snap_idx = data_dict[:snap_data].global_index

            integrated_qty = data_dict[:sim_data].snapshot_table[snap_idx, :lookback_times]

        #########################
        # Metallicity quantities
        #########################

        elseif quantity == :gas_metallicity

            Mz = sum(computeMass(data_dict, :Z_gas); init=0.0u"Msun")
            Mg = sum(computeMass(data_dict, :gas); init=0.0u"Msun")

            if iszero(Mg)
                integrated_qty = NaN
            else
                integrated_qty = (Mz / Mg) / SOLAR_METALLICITY
            end

        elseif quantity == :stellar_metallicity

            Mz = sum(computeMass(data_dict, :Z_stellar); init=0.0u"Msun")
            Ms = sum(computeMass(data_dict, :stellar); init=0.0u"Msun")

            if iszero(Ms)
                integrated_qty = NaN
            else
                integrated_qty = (Mz / Ms) / SOLAR_METALLICITY
            end

        elseif quantity == :ode_metallicity

            Mz = sum(computeMass(data_dict, :ode_metals); init=0.0u"Msun")
            Mg = sum(computeMass(data_dict, :gas); init=0.0u"Msun")

            if iszero(Mg)
                integrated_qty = NaN
            else
                integrated_qty = (Mz / Mg) / SOLAR_METALLICITY
            end

        elseif quantity ∈ GAS_ABUNDANCE

            element = GAS_ABUNDANCE_SPLITS[quantity]

            abundance = computeGlobalAbundance(data_dict, :gas, element; solar=false)

            if iszero(abundance)
                integrated_qty = NaN
            else
                integrated_qty = ABUNDANCE_SHIFT[element] + log10(abundance)
            end

        elseif quantity ∈ STELLAR_ABUNDANCE

            element = STELLAR_ABUNDANCE_SPLITS[quantity]

            abundance = computeGlobalAbundance(data_dict, :stellar, element; solar=false)

            if iszero(abundances)
                integrated_qty = NaN
            else
                integrated_qty = ABUNDANCE_SHIFT[element] + log10(abundances)
            end

        else

            throw(ArgumentError("integrateQty: I don't recognize the quantity :$(quantity)"))

        end

    else

        if quantity ∈ QTY_GLOBAL_LIST

            scatter_qty = scatterQty(data_dict, quantity)

            if isempty(scatter_qty)

                logging[] && @warn("integrateQty: :$(quantity) is empty, so I will return NaN")

                unit = plotParams(quantity).unit

                integrated_qty = NaN * unit

            else

                integrated_qty = agg_function(scatter_qty)

            end

        else

            throw(ArgumentError("integrateQty: The quantity :$(quantity) is not compatible with a \
            personalized aggregator function like $(agg_function). Use a quantity listed in \
            DERIVED_QTY, SFM_QTY, SFM_DERIVED_QTY or RATIO_QTY or set agg_function = :default"))

        end

    end

    return integrated_qty

end
