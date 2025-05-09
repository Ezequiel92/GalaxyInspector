####################################################################################################
# Compute characteristic masses and mass related quantities
####################################################################################################

"""
    computeFraction(data_dict::Dict, component::Symbol)::Vector{Float64}

Compute the fraction in each cell/particle of a given `component`.

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

  - `component::Symbol`: For which cell/particle type the fraction will be calculated. The options are:

      + `:molecular`    -> Molecular hydrogen (``\\mathrm{H_2}``) fraction.
      + `:br_molecular` -> Molecular hydrogen (``\\mathrm{H_2}``) fraction, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic`       -> Atomic hydrogen (``\\mathrm{HI}``) fraction.
      + `:ionized`      -> Ionized hydrogen (``\\mathrm{HII}``) fraction.
      + `:neutral`      -> Neutral hydrogen (``\\mathrm{HI + H_2}``) fraction.
      + `:stellar`      -> Stellar gas fraction (according to our SF model).
      + `:metals`       -> Metallicity (according to our SF model).
      + `:dust`         -> Dust fraction (according to our SF model).

# Returns

  - The fraction of `component` in each cell/particle.

# References

L. Blitz et al. (2006). *The Role of Pressure in GMC Formation II: The H2-Pressure Relation*. The Astrophysical Journal, **650(2)**, 933. [doi:10.1086/505417](https://doi.org/10.1086/505417)
"""
function computeFraction(data_dict::Dict, component::Symbol)::Vector{Float64}

    dg = data_dict[:gas]
    gas_mass = dg["MASS"]

    # If there is no gas, return an empty array
    !isempty(gas_mass) || return Float64[]

    # Compute the number of gas cells
    n_cells = length(gas_mass)

    ################################################################################################
    # Molecular fraction
    ################################################################################################

    if component == :molecular

        # Read the neutral and ionized hydrogen fractions from Arepo
        nh = dg["NH  "]
        nhp = dg["NHP "]

        if !isempty(dg["FRAC"])

            (
                !logging[] ||
                @info("computeFraction: The molecular fraction will be calculated using the \
                fractions from our SF model")
            )

            molecular_fraction = view(dg["FRAC"], 3, :)
            stellar_fraction = view(dg["FRAC"], 4, :)
            ρ_gas = dg["RHO "]

            fractions = Vector{Float64}(undef, n_cells)

            for i in eachindex(fractions)

                if !isnan(molecular_fraction[i]) && ρ_gas[i] >= THRESHOLD_DENSITY

                    # Fraction of molecular hydrogen according to our SF model
                    fractions[i] = molecular_fraction[i] + stellar_fraction[i]

                else

                    # When there is no data from our model or the density is below the SF threshold,
                    # assume no molecular hydrogen
                    fractions[i] = 0.0

                end

            end

        elseif !any(isempty, [nh, nhp, dg["PRES"]])

            (
                !logging[] ||
                @info("computeFraction: The molecular fraction will be calculated using the \
                fraction of neutral gas from Arepo and the pressure relation from Blitz et al. \
                (2006)")
            )

            relative_pressure = @. uconvert(Unitful.NoUnits, dg["PRES"] / P0)^ALPHA_BLITZ

            # Compute the fraction of neutral gas that is molecular hydrogen according
            # to the pressure relation in Blitz et al. (2006)
            fm = @. 1.0 / (1.0 + relative_pressure)

            # Compute the fraction of neutral gas according to Arepo
            fn = @. nh / (nhp + nh)

            fractions = fm .* fn

        else

            !logging[] || @warn("computeFraction: I could not compute the molecular fraction")

            fractions = Float64[]

        end

    ################################################################################################
    # BR molecular fraction
    ################################################################################################

    elseif component == :br_molecular

        if !isempty(dg["PRES"])

            (
                !logging[] ||
                @info("computeFraction: The BR molecular fraction will be calculated using only \
                the pressure relation from Blitz et al. (2006)")
            )

            relative_pressure = @. uconvert(Unitful.NoUnits, dg["PRES"] / P0)^ALPHA_BLITZ

            # Compute the fraction of gas that is molecular hydrogen according
            # to the pressure relation in Blitz et al. (2006)
            fractions = @. 1.0 / (1.0 + relative_pressure)

        else

            !logging[] || @warn("computeFraction: I could not compute the BR molecular fraction")

            fractions = Float64[]

        end

    ################################################################################################
    # Atomic fraction
    ################################################################################################

    elseif component == :atomic

        # Read the neutral and ionized hydrogen fractions from Arepo
        nh = dg["NH  "]
        nhp = dg["NHP "]

        if !any(isempty, [nh, nhp, dg["FRAC"]])

            (
                !logging[] ||
                @info("computeFraction: The atomic fraction will be calculated using the fractions \
                from our SF model")
            )

            atomic_fraction = view(dg["FRAC"], 2, :)
            ρ_gas = dg["RHO "]
            gas_metallicity = dg["GZ  "]

            fractions = Vector{Float64}(undef, n_cells)

            for i in eachindex(fractions)

                if !isnan(atomic_fraction[i]) && ρ_gas[i] >= THRESHOLD_DENSITY

                    # Fraction of atomic hydrogen according to our SF model
                    fractions[i] = atomic_fraction[i]

                else

                    Z = setPositive(gas_metallicity[i])

                    # When there is no data from our model or the density is below the SF threshold,
                    # use the fraction of neutral hydrogen from Arepo
                    fractions[i] = (1 - Z) * nh[i] / (nhp[i] + nh[i])

                end

            end

        elseif !any(isempty, [nh, nhp, dg["PRES"]])

            (
                !logging[] ||
                @info("computeFraction: The atomic fraction will be calculated using the fraction \
                of neutral gas from Arepo and the pressure relation from Blitz et al. (2006)")
            )

            relative_pressure = @. uconvert(Unitful.NoUnits, dg["PRES"] / P0)^ALPHA_BLITZ

            # Compute the fraction of neutral gas that is atomic hydrogen according
            # to the pressure relation in Blitz et al. (2006)
            fa = @. relative_pressure / (1.0 + relative_pressure)

            # Compute the fraction of neutral gas according to Arepo
            fn = @. nh / (nhp + nh)

            fractions = fa .* fn

        else

            !logging[] || @warn("computeFraction: I could not compute the atomic fraction")

            fractions = Float64[]

        end

    ################################################################################################
    # Ionized fraction
    ################################################################################################

    elseif component == :ionized

        # Read the neutral and ionized hydrogen fractions from Arepo
        nh = dg["NH  "]
        nhp = dg["NHP "]

        if !any(isempty, [nh, nhp, dg["FRAC"]])

            (
                !logging[] ||
                @info("computeFraction: The ionized fraction will be calculated using the \
                fractions from our SF model")
            )

            ionized_fraction = view(dg["FRAC"], 1, :)
            ρ_gas = dg["RHO "]
            gas_metallicity = dg["GZ  "]

            fractions = Vector{Float64}(undef, n_cells)

            for i in eachindex(fractions)

                if !isnan(ionized_fraction[i]) && ρ_gas[i] >= THRESHOLD_DENSITY

                    # Fraction of ionized hydrogen according to our SF model
                    fractions[i] = ionized_fraction[i]

                else

                    Z = setPositive(gas_metallicity[i])

                    # When there is no data from our model or the density is below the SF threshold,
                    # use the fraction of ionized hydrogen from Arepo
                    fractions[i] = (1 - Z) * nhp[i] / (nhp[i] + nh[i])

                end

            end

        elseif !any(isempty, [nh, nhp])

            (
                !logging[] ||
                @info("computeFraction: The ionized fraction will be calculated using the fraction \
                of ionized gas from Arepo")
            )

            # Compute the fraction of ionized gas according to Arepo
            fractions = @. nhp / (nhp + nh)

        else

            !logging[] || @warn("computeFraction: I could not compute the ionized fraction")

            fractions = Float64[]

        end

    ################################################################################################
    # Neutral fraction
    ################################################################################################

    elseif component == :neutral

        # Read the neutral and ionized hydrogen fractions from Arepo
        nh = dg["NH  "]
        nhp = dg["NHP "]

        if !any(isempty, [nh, nhp, dg["FRAC"]])

            (
                !logging[] ||
                @info("computeFraction: The neutral fraction will be calculated using the \
                fractions from our SF model")
            )

            atomic_fraction = view(dg["FRAC"], 2, :)
            molecular_fraction = view(dg["FRAC"], 3, :)
            stellar_fraction = view(dg["FRAC"], 4, :)
            ρ_gas = dg["RHO "]
            gas_metallicity = dg["GZ  "]

            fa = Vector{Float64}(undef, n_cells)
            fm = Vector{Float64}(undef, n_cells)

            for i in eachindex(fa)

                if !isnan(atomic_fraction[i]) && ρ_gas[i] >= THRESHOLD_DENSITY

                    # Fraction of atomic hydrogen according to our SF model
                    fa[i] = atomic_fraction[i]

                    # Fraction of molecular hydrogen according to our SF model
                    fm[i] = molecular_fraction[i] + stellar_fraction[i]

                else

                    Z = setPositive(gas_metallicity[i])

                    # When there is no data from our model or the density is below the SF threshold,
                    # use the fraction of neutral hydrogen from Arepo as the atomic fraction
                    # and set the molecular fraction to 0
                    fa[i] = (1 - Z) * nh[i] / (nhp[i] + nh[i])
                    fm[i] = 0.0

                end

            end

            fractions = fa .+ fm

        elseif !any(isempty, [nh, nhp])

            (
                !logging[] ||
                @info("computeFraction: The neutral fraction will be calculated using the fraction \
                of neutral gas from Arepo")
            )

            # Compute the fraction of neutral gas according to Arepo
            fractions = @. nh / (nhp + nh)

        else

            !logging[] || @warn("computeFraction: I could not compute the neutral fraction")

            fractions = Float64[]

        end

    ################################################################################################
    # Stellar fraction
    ################################################################################################

    elseif component == :stellar

        if !isempty(dg["FRAC"])

            (
                !logging[] ||
                @info("computeFraction: The stellar fraction will be calculated using the \
                fractions from our SF model")
            )

            stellar_fraction = view(dg["FRAC"], 4, :)
            ρ_gas = dg["RHO "]

            fractions = Vector{Float64}(undef, n_cells)

            for i in eachindex(fractions)

                if !isnan(stellar_fraction[i]) && ρ_gas[i] >= THRESHOLD_DENSITY

                    # Fraction of stars according to our SF model
                    fractions[i] = stellar_fraction[i]

                else

                    # When there is no data from our model or the density is below the SF threshold,
                    # assume no stellar fraction
                    fractions[i] = 0.0

                end

            end

        else

            !logging[] || @warn("computeFraction: I could not compute the stellar fraction")

            fractions = Float64[]

        end

    ################################################################################################
    # Metal fraction
    ################################################################################################

    elseif component == :metals

        # Read the neutral and ionized hydrogen fractions from Arepo
        nh = dg["NH  "]
        nhp = dg["NHP "]

        if !any(isempty, [nh, nhp, dg["FRAC"]])

            (
                !logging[] ||
                @info("computeFraction: The metallicity will be calculated using the \
                fractions from our SF model")
            )

            metal_fraction = view(dg["FRAC"], 5, :)
            ρ_gas = dg["RHO "]
            gas_metallicity = dg["GZ  "]

            fractions = Vector{Float64}(undef, n_cells)

            for i in eachindex(fractions)

                if !isnan(metal_fraction[i]) && ρ_gas[i] >= THRESHOLD_DENSITY

                    # Metallicity according to our SF model
                    fractions[i] = metal_fraction[i]

                else

                    Z = setPositive(gas_metallicity[i])

                    fa = (1 - Z) * nh[i] / (nhp[i] + nh[i])

                    # When there is no data from our model or the density is below the SF threshold,
                    # use the fraction of neutral hydrogen and the metallicity from Arepo
                    fractions[i] = Z * (1.0 - C_xd * fa)

                end

            end

        else

            !logging[] || @warn("computeFraction: I could not compute the metallicity")

            fractions = Float64[]

        end

    ################################################################################################
    # Dust fraction
    ################################################################################################

    elseif component == :dust

        # Read the neutral and ionized hydrogen fractions from Arepo
        nh = dg["NH  "]
        nhp = dg["NHP "]

        if !any(isempty, [nh, nhp, dg["FRAC"]])

            (
                !logging[] ||
                @info("computeFraction: The dust fraction will be calculated using the fractions \
                from our SF model")
            )

            dust_fraction = view(dg["FRAC"], 6, :)
            ρ_gas = dg["RHO "]
            gas_metallicity = dg["GZ  "]

            fractions = Vector{Float64}(undef, n_cells)

            for i in eachindex(fractions)

                if !isnan(dust_fraction[i]) && ρ_gas[i] >= THRESHOLD_DENSITY

                    # Fraction of dust according to our SF model
                    fractions[i] = dust_fraction[i]

                else

                    Z = setPositive(gas_metallicity[i])

                    fa = (1 - Z) * nh[i] / (nhp[i] + nh[i])

                    # When there is no data from our model or the density is below the SF threshold,
                    # use the fraction of neutral hydrogen and the metallicity from Arepo
                    fractions[i] = Z * C_xd * fa

                end

            end

        else

            !logging[] || @warn("computeFraction: I could not compute the dust fraction")

            fractions = Float64[]

        end

    else

        throw(ArgumentError("computeFraction: I don't recognize the component :$(component)"))

    end

    return fractions

end

"""
    computeMass(data_dict::Dict, component::Symbol)::Vector{<:Unitful.Mass}

Compute the mass in each cell/particle of a given `component`.

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
  - `component::Symbol`: For which cell/particle type the mass will be calculated. The options are:

      + `:gas`          -> Gas mass.
      + `:stars`        -> Stellar mass.
      + `:hydrogen`     -> Hydrogen mass.
      + `:helium`       -> Helium mass.
      + `:dark_matter`  -> Dark matter mass.
      + `:black_holes`  -> Black hole mass.
      + `:molecular`    -> Molecular hydrogen (``\\mathrm{H_2}``) mass.
      + `:br_molecular` -> Molecular hydrogen (``\\mathrm{H_2}``) mass, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic`       -> Atomic hydrogen (``\\mathrm{HI}``) mass.
      + `:ionized`      -> Ionized hydrogen (``\\mathrm{HII}``) mass.
      + `:neutral`      -> Neutral hydrogen (``\\mathrm{HI + H_2}``) mass.
      + `:stellar`      -> Stellar gas mass (according to our SF model).
      + `:metals`       -> Metal mass (according to our SF model).
      + `:dust`         -> Dust mass.

# Returns

  - The mass of `component` in each cell/particle.

# References

L. Blitz et al. (2006). *The Role of Pressure in GMC Formation II: The H2-Pressure Relation*. The Astrophysical Journal, **650(2)**, 933. [doi:10.1086/505417](https://doi.org/10.1086/505417)
"""
function computeMass(data_dict::Dict, component::Symbol)::Vector{<:Unitful.Mass}

    if component == :gas

        masses = data_dict[:gas]["MASS"]

    elseif component == :stars

        masses = data_dict[:stars]["MASS"]

    elseif component == :hydrogen

        masses = data_dict[:gas]["MASS"] * HYDROGEN_MASSFRAC

    elseif component == :helium

        masses = data_dict[:gas]["MASS"] * (1.0 - HYDROGEN_MASSFRAC)

    elseif component == :dark_matter

        masses = data_dict[:halo]["MASS"]

    elseif component == :black_holes

        masses = data_dict[:black_hole]["MASS"]

    elseif component ∈ [
        :molecular,
        :br_molecular,
        :atomic,
        :ionized,
        :neutral,
        :stellar,
        :metals,
        :dust,
    ]

        fractions = computeFraction(data_dict, component)

        if isempty(fractions)
            masses = Unitful.Mass[]
        else
            masses = data_dict[:gas]["MASS"] .* fractions
        end

    else

        throw(ArgumentError("computeMass: I don't recognize the component :$(component)"))

    end

    !isempty(masses) || return Unitful.Mass[]

    return masses

end

"""
    computeVolumeDensity(data_dict::Dict, component::Symbol)::Vector{<:Unitful.Density}

Compute the volume mass density in each cell/particle of a given `component`.

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
  - `component::Symbol`: For which cell/particle type the mass will be calculated. The options are:

      + `:gas`          -> Gas mass density.
      + `:hydrogen`     -> Hydrogen mass density.
      + `:helium`       -> Helium mass density.
      + `:molecular`    -> Molecular hydrogen (``\\mathrm{H_2}``) mass density.
      + `:br_molecular` -> Molecular hydrogen (``\\mathrm{H_2}``) mass density, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic`       -> Atomic hydrogen (``\\mathrm{HI}``) mass density.
      + `:ionized`      -> Ionized hydrogen (``\\mathrm{HII}``) mass density.
      + `:neutral`      -> Neutral hydrogen (``\\mathrm{HI + H_2}``) mass density.
      + `:stellar`      -> Stellar mass density (according to our SF model).
      + `:metals`       -> Metal mass density (according to our SF model).
      + `:dust`         -> Dust mass density.

# Returns

  - The density of `component` in each cell/particle.

# References

L. Blitz et al. (2006). *The Role of Pressure in GMC Formation II: The H2-Pressure Relation*. The Astrophysical Journal, **650(2)**, 933. [doi:10.1086/505417](https://doi.org/10.1086/505417)
"""
function computeVolumeDensity(data_dict::Dict, component::Symbol)::Vector{<:Unitful.Density}

    if component == :gas

        densities = data_dict[:gas]["RHO "]

    elseif component == :hydrogen

        densities = data_dict[:gas]["RHO "] * HYDROGEN_MASSFRAC

    elseif component == :helium

        densities = data_dict[:gas]["RHO "] * (1.0 - HYDROGEN_MASSFRAC)

    elseif component ∈ [
        :molecular,
        :br_molecular,
        :atomic,
        :ionized,
        :neutral,
        :stellar,
        :metals,
        :dust,
    ]

        fractions = computeFraction(data_dict, component)

        if isempty(fractions)
            densities = Unitful.Density[]
        else
            densities = data_dict[:gas]["RHO "] .* fractions
        end

    else

        throw(ArgumentError("computeVolumeDensity: I don't recognize the component :$(component)"))

    end

    !isempty(densities) || return Unitful.Density[]

    return densities

end

"""
    computeNumberDensity(data_dict::Dict, component::Symbol)::Vector{<:NumberDensity}

Compute the number density in each cell/particle of a given `component`.

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
  - `component::Symbol`: For which cell/particle type the mass will be calculated. The options are:

      + `:gas`          -> Gas number density.
      + `:hydrogen`     -> Hydrogen number density.
      + `:helium`       -> Helium number density.
      + `:molecular`    -> Molecular hydrogen (``\\mathrm{H_2}``) number density.
      + `:br_molecular` -> Molecular hydrogen (``\\mathrm{H_2}``) number density, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic`       -> Atomic hydrogen (``\\mathrm{HI}``) number density.
      + `:ionized`      -> Ionized hydrogen (``\\mathrm{HII}``) number density.
      + `:neutral`      -> Neutral hydrogen (``\\mathrm{HI + H_2}``) number density.

# Returns

  - The number density of `component` in each cell/particle.

# References

L. Blitz et al. (2006). *The Role of Pressure in GMC Formation II: The H2-Pressure Relation*. The Astrophysical Journal, **650(2)**, 933. [doi:10.1086/505417](https://doi.org/10.1086/505417)
"""
function computeNumberDensity(data_dict::Dict, component::Symbol)::Vector{<:NumberDensity}

    densities = computeVolumeDensity(data_dict, component)

    if component ∈ [:gas, :hydrogen, :atomic, :ionized, :neutral]

        n_densities = densities / Unitful.mp

    elseif component ∈ [:helium, :molecular, :br_molecular]

        n_densities = densities / (2 * Unitful.mp)

    else

        throw(ArgumentError("computeNumberDensity: I don't recognize the component :$(component)"))

    end

    !isempty(n_densities) || return NumberDensity[]

    return n_densities

end

"""
    computeMetalMass(data_dict::Dict, component::Symbol)::Vector{<:Unitful.Mass}

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
  - `component::Symbol`: For which cell/particle type the mass will be calculated. The possibilities are `:stars` and `:gas`.

# Returns

  - The total metal mass in each cell/particle.
"""
function computeMetalMass(data_dict::Dict, component::Symbol)::Vector{<:Unitful.Mass}

    if component == :gas

        block = "GZ  "

    elseif component == :stars

        block = "GZ2 "

    else

        throw(ArgumentError("computeMetalMass: `component` can only be :stars or :gas, \
        but I got :$(component)"))

    end

    if !isempty(data_dict[component][block])

        masses = setPositive(data_dict[component][block]) .* computeMass(data_dict, component)

    else

        (
            !logging[] ||
            @warn("computeMetalMass: I could not compute the masses of metals")
        )

        masses = Unitful.Mass[]

    end

    !isempty(masses) || return Unitful.Mass[]

    return masses

end

"""
    computeElementMass(
        data_dict::Dict,
        component::Symbol,
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
  - `component::Symbol`: For which cell/particle type the mass will be calculated. The possibilities are `:stars` and `:gas`.
  - `element::Symbol`: Target element. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).

# Returns

  - The total mass of `element` in each cell/particle.
"""
function computeElementMass(
    data_dict::Dict,
    component::Symbol,
    element::Symbol,
)::Vector{<:Unitful.Mass}

    (
        element ∈ keys(ELEMENT_INDEX) ||
        throw(ArgumentError("computeElementMass: :$(element) is not a tracked element, \
        the options are the keys of `ELEMENT_INDEX`, see `./src/constants/arepo.jl`"))
    )

    if component == :gas

        block = "GMET"

    elseif component == :stars

        block = "GME2"

    else

        throw(ArgumentError("computeElementMass: `component` can only be :stars or :gas, \
        but I got :$(component)"))

    end

    if !isempty(data_dict[component][block])

        element_fractions = setPositive(data_dict[component][block][ELEMENT_INDEX[element], :])

        masses = element_fractions .* data_dict[component]["MASS"]

    else

        (
            !logging[] ||
            @warn("computeElementMass: I could not compute the masses of :$(element)")
        )

        masses = Unitful.Mass[]

    end

    !isempty(masses) || return Unitful.Mass[]

    return masses

end

@doc raw"""
    computeAbundance(
        data_dict::Dict,
        component::Symbol,
        element::Symbol;
        <keyword arguments>
    )::Float64

Compute the abundance of a given element in each cell/particle of type `component`. The abundance is defined as $n_X / n_H$ where $n_X$ is the number of atoms of element $X$ and $n_H$ the number of hydrogen atoms.

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
  - `component::Symbol`: For which cell/particle type the abundance will be calculated. The possibilities are `:stars` and `:gas`.
  - `element::Symbol`: Target element. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
  - `solar::Bool=false`: If the result will be normalized to the solar abundance or not.

# Returns

  - The abundance of `element` in each cell/particle of type `component`.
"""
function computeAbundance(
    data_dict::Dict,
    component::Symbol,
    element::Symbol;
    solar::Bool=false,
)::Float64

    element_mass = computeElementMass(data_dict, component, element)
    hydrogen_mass = computeElementMass(data_dict, component, :H)

    (
        !isempty(hydrogen_mass) ||
        throw(ArgumentError("computeElementMass: I got no data for the mass of hydrogen. \
        This should not be possible!"))
    )

    # Compute the number of atoms
    n_X = element_mass ./ ATOMIC_WEIGHTS[element]
    n_H = hydrogen_mass ./ ATOMIC_WEIGHTS[:H]

    # Compute the relative abundance of `element`
    abundances = ustrip.(Unitful.NoUnits, n_X ./ n_H)

    return abundances ./ (solar ? exp10(SOLAR_ABUNDANCE[element] - 12.0) : 1.0)

end

@doc raw"""
    computeGlobalAbundance(
        data_dict::Dict,
        component::Symbol,
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
  - `component::Symbol`: For which cell/particle type the abundance will be calculated. The possibilities are `:stars` and `:gas`.
  - `element::Symbol`: Target element. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
  - `solar::Bool=false`: If the result will be normalized to the solar abundance or not.

# Returns

  - The total abundance of `element`.
"""
function computeGlobalAbundance(
    data_dict::Dict,
    component::Symbol,
    element::Symbol;
    solar::Bool=false,
)::Float64

    metal_mass = sum(computeElementMass(data_dict, component, element); init=0.0u"Msun")
    hydrogen_mass = sum(computeElementMass(data_dict, component, :H); init=0.0u"Msun")

    (
        !iszero(hydrogen_mass) ||
        throw(ArgumentError("computeGlobalAbundance: I got 0 for the mass of hydrogen. \
        This should not be possible!"))
    )

    # Compute the number of atoms
    n_X = metal_mass / ATOMIC_WEIGHTS[element]
    n_H = hydrogen_mass / ATOMIC_WEIGHTS[:H]

    # Compute the relative abundance of `element`
    abundance = ustrip(Unitful.NoUnits, n_X / n_H)

    return abundance / (solar ? exp10(SOLAR_ABUNDANCE[element] - 12.0) : 1.0)

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

    # Find the tracers that are inside R200 now, but were outside R200 in the past
    inflow_ids  = setdiff(present_tracer_ids, past_tracer_ids)
    # Find the tracers that were inside R200 in the past, but are now outside R200
    outflow_ids = setdiff(past_tracer_ids, present_tracer_ids)

    # Compute the mass of each tracer in physical units
    tracer_mass = TRACER_MASS * internalUnits("MASS", present_dd[:snap_data].path)

    # Compute the inflow mass
    inflow_mass = length(inflow_ids) * tracer_mass
    # Compute the outflow mass
    outflow_mass = length(outflow_ids) * tracer_mass

    # Compute the net mass
    net_mass_increase = inflow_mass - outflow_mass

    return net_mass_increase, inflow_mass, outflow_mass

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
  - `max_r::Unitful.Length=DISK_R`: Radius of the cylinder.
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
    max_r::Unitful.Length=DISK_R,
    max_z::Unitful.Length=5.0u"kpc",
)::NTuple{3,Unitful.Mass}

    # Find the tracers inside a given cylinder in the present snapshot
    present_tracer_ids = tracersWithinDisc(present_dd; max_r, max_z)
    # Find the tracers inside a given cylinder in the past snapshot
    past_tracer_ids    = tracersWithinDisc(past_dd; max_r, max_z)

    # Find the tracers that are inside a given cylinder now, but were outside in the past
    inflow_ids  = setdiff(present_tracer_ids, past_tracer_ids)
    # Find the tracers that were inside a given cylinder in the past, but are now outside
    outflow_ids = setdiff(past_tracer_ids, present_tracer_ids)

    # Compute the mass of each tracer in physical units
    tracer_mass = TRACER_MASS * internalUnits("MASS", present_dd[:snap_data].path)

    # Compute the inflow mass
    inflow_mass = length(inflow_ids) * tracer_mass
    # Compute the outflow mass
    outflow_mass = length(outflow_ids) * tracer_mass

    # Compute the net mass
    net_mass_increase = inflow_mass - outflow_mass

    return net_mass_increase, inflow_mass, outflow_mass

end

"""
    computeMassRadius(
        positions::Matrix{<:Unitful.Length},
        masses::Vector{<:Unitful.Mass};
        <keyword arguments>
    )::Unitful.Length

Compute the radius containing `percent`% of the total mass.

# Arguments

  - `positions::Matrix{<:Unitful.Length}`: Positions of the cells/particles. Each column is a cell/particle and each row a dimension.
  - `masses::Vector{<:Unitful.Mass}`: Masses of the cells/particles.
  - `percent::Float64=90.0`: Target percentage of the total mass.

# Returns

  - The radius containing `percent`% of the total mass.
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

Compute the total height of a cylinder, of infinite radius, containing `percent`% of the total mass.

# Arguments

  - `positions::Matrix{<:Unitful.Length}`: Positions of the cells/particles. Each column is a cell/particle and each row a dimension.
  - `masses::Vector{<:Unitful.Mass}`: Masses of the cells/particles.
  - `percent::Float64=90.0`: Target percentage of the total mass.

# Returns

  - The height containing `percent`% of the total mass.
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

Compute the maximum value of `quantity` that "contains" `percent`% of the total mass.

# Arguments

  - `quantity::Vector{<:Number}`: Target quantity.
  - `masses::Vector{<:Unitful.Mass}`: Masses of the cells/particles.
  - `percent::Float64=90.0`: Target percentage of the total mass.

# Returns

  - The maximum value of `quantity` that "contains" `percent`% of the total mass.
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
    computeMassFraction(
        quantity::Vector{<:Number},
        masses::Vector{<:Unitful.Mass},
        qty_limits::Tuple{<:Number,<:Number},
    )::Float64

Compute the fraction of the total mass "contained" within given values of `quantity`.

# Arguments

  - `quantity::Vector{<:Number}`: Target quantity.
  - `masses::Vector{<:Unitful.Mass}`: Masses of the cells/particles.
  - `qty_limits::Tuple{<:Number,<:Number}`: Limits of the target quantity.

# Returns

  - The fraction of the total mass "contained" within given values of `quantity`.
"""
function computeMassFraction(
    quantity::Vector{<:Number},
    masses::Vector{<:Unitful.Mass},
    qty_limits::Tuple{<:Number,<:Number},
)::Float64

    # Check for missing data
    !any(isempty, [quantity, masses]) || return 0.0

    # Find the indices of all the cells/particles with `qty_limits[1]` <= `quantity` <= `qty_limits[2]`
    idxs = map(x -> qty_limits[1] <= x <= qty_limits[2], quantity)

    return uconvert(Unitful.NoUnits, sum(masses[idxs]; init=0.0u"Msun") / sum(masses))

end
