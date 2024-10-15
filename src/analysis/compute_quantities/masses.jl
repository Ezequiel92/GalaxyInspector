####################################################################################################
# Compute characteristic masses and mass related quantities
####################################################################################################

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
        qty_limits::Tuple{<:Number,<:Number},
    )::Float64

Compute the fraction of the total mass "contained" within a given values of `quantity`.

# Arguments

  - `quantity::Vector{<:Number}`: Target quantity.
  - `masses::Vector{<:Unitful.Mass}`: Masses of the cells/particles.
  - `qty_limits::Tuple{<:Number,<:Number}`: Limits of the target quantity.

# Returns

  - The fraction of the total mass "contained" within a given value of `quantity`.
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
  - `component::Symbol`: For which cell/particle type the metal mass will be calculated. The possibilities are `:stars` and `:gas`.

# Returns

  - The total metal mass in each cell/particle.
"""
function computeMetalMass(data_dict::Dict, component::Symbol)::Vector{<:Unitful.Mass}

    if CODEBASE == :arepo

        if component == :gas
            metals = setPositive(data_dict[:gas]["GZ  "]) .* data_dict[:gas]["MASS"]
        elseif component == :stars
            metals = setPositive(data_dict[:stars]["GZ2 "]) .* data_dict[:stars]["MASS"]
        else
            throw(ArgumentError("computeMetalMass: `component` can only be :stars or :gas, \
            but I got :$(component)"))
        end

    elseif CODEBASE == :opengadget3

        if component == :gas
            metals = sum(setPositive(data_dict[:gas]["GMET"][METAL_LIST, :]); dims=1)
        elseif component == :stars
            metals = sum(setPositive(data_dict[:stars]["GME2"][METAL_LIST, :]); dims=1)
        else
            throw(ArgumentError("computeMetalMass: `component` can only be :stars or :gas, \
            but I got :$(component)"))
        end

    else

        throw(ArgumentError("computeMetalMass: I don't recognize the codebase :$(CODEBASE)"))

    end

    return metals

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
  - `component::Symbol`: For which cell/particle type the element mass will be calculated. The possibilities are `:stars` and `:gas`.
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
        the options are the keys of `ELEMENT_INDEX`, see `./src/constants/globals.jl`"))
    )

    if component == :gas
        block = "GMET"
    elseif component == :stars
        block = "GME2"
    else
        throw(ArgumentError("computeElementMass: `component` can only be :stars or :gas, \
        but I got :$(component)"))
    end

    values = data_dict[component]

    if any(isempty, [values[block], values["MASS"]])
        return Unitful.Mass[]
    end

    if CODEBASE == :arepo
        masses = setPositive(values[block][ELEMENT_INDEX[element], :]) .* values["MASS"]
    elseif CODEBASE == :opengadget3
        masses = setPositive(values[block][ELEMENT_INDEX[element], :])
    else
        throw(ArgumentError("computeElementMass: I don't recognize the codebase :$(CODEBASE)"))
    end

    return masses

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
    computeIonizedMass(data_dict::Dict)::Vector{<:Unitful.Mass}

Compute the mass of ionized hydrogen in every gas cell/particle.

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

  - The mass of ionized hydrogen in every gas cell/particle.
"""
function computeIonizedMass(data_dict::Dict)::Vector{<:Unitful.Mass}

    dg = data_dict[:gas]

    !isempty(dg["MASS"]) || return Unitful.Mass[]

    if "FRAC" ∈ keys(dg) && !any(isempty, [dg["FRAC"], dg["NHP "], dg["NH  "]])

        fi = Vector{Float64}(undef, length(dg["MASS"]))

        @inbounds for i in eachindex(fi)

            # Compute how much time has pass since the last time the cell/particle entered the SF routine
            Δt = data_dict[:snap_data].physical_time - dg["CTIM"][i]

            @inbounds if !isnan(dg["FRAC"][1, i]) && Δt < dg["TAUS"][i]

                # Fraction of ionized hydrogen according to our SF model
                fi[i] = dg["FRAC"][1, i]

            else

                # When there is no data from our model, use the fraction of ionized hydrogen from "NHP "
                fi[i] = dg["NHP "][i] / (dg["NHP "][i] + dg["NH  "][i])

            end

        end

    elseif !isempty(dg["COLM"])

        return dg["MASS"] .- dg["COLM"]

    elseif !any(isempty, [dg["NHP "], dg["NH  "]])

        fi = dg["NHP "][i] ./ (dg["NHP "][i] .+ dg["NH  "][i])

    else

        return Unitful.Mass[]

    end

    return fi .* dg["MASS"]

end

"""
    computeNeutralMass(data_dict::Dict; <keyword arguments>)::Vector{<:Unitful.Mass}

Compute the mass of neutral hydrogen in every gas cell/particle.

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
  - `normalize::Bool=true`: If the output will be normalize to eliminate the stellar fraction of the gas cells. Only relevant for simulation with our SF routine, and for cells that have entered it at least once.

# Returns

  - The mass of neutral hydrogen in every gas cell/particle.
"""
function computeNeutralMass(data_dict::Dict; normalize::Bool=true)::Vector{<:Unitful.Mass}

    dg = data_dict[:gas]

    !isempty(dg["MASS"]) || return Unitful.Mass[]

    if "FRAC" ∈ keys(dg) && !any(isempty, [dg["FRAC"], dg["NHP "], dg["NH  "]])

        fa = Vector{Float64}(undef, length(dg["MASS"]))
        fm = Vector{Float64}(undef, length(dg["MASS"]))

        @inbounds for i in eachindex(fa)

            # Compute how much time has pass since the last time the cell/particle entered the SF routine
            Δt = data_dict[:snap_data].physical_time - dg["CTIM"][i]

            @inbounds if !isnan(dg["FRAC"][2, i]) && Δt < dg["TAUS"][i]

                # Fraction of atomic hydrogen according to our model
                fa[i] = dg["FRAC"][2, i]

                # Fraction of molecular hydrogen according to our model
                @inbounds if normalize
                    fm[i] = dg["FRAC"][3, i] + dg["FRAC"][4, i]
                else
                    fm[i] = dg["FRAC"][3, i]
                end

            else

                # When there is no data from our model, use the fraction of neutral hydrogen from "NH  "
                fa[i] = dg["NH  "][i] / (dg["NHP "][i] + dg["NH  "][i])
                fm[i] = 0.0

            end

        end

        fn = fa .+ fm

    elseif !isempty(dg["COLM"])

        return dg["COLM"]

    elseif !any(isempty, [dg["NHP "], dg["NH  "]])

        fn = dg["NH  "][i] ./ (dg["NHP "][i] .+ dg["NH  "][i])

    else

        return Unitful.Mass[]

    end

    return fn .* dg["MASS"]

end

"""
    computeAtomicMass(data_dict::Dict)::Vector{<:Unitful.Mass}

Compute the mass of atomic hydrogen in every gas cell/particle.

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

  - The mass of atomic hydrogen in every gas cell/particle.
"""
function computeAtomicMass(data_dict::Dict)::Vector{<:Unitful.Mass}

    dg = data_dict[:gas]

    !isempty(dg["MASS"]) || return Unitful.Mass[]

    if "FRAC" ∈ keys(dg) && !any(isempty, [dg["FRAC"], dg["NHP "], dg["NH  "]])

        fa = Vector{Float64}(undef, length(dg["MASS"]))

        @inbounds for i in eachindex(fa)

            # Compute how much time has pass since the last time the cell/particle entered the SF routine
            Δt = data_dict[:snap_data].physical_time - dg["CTIM"][i]

            @inbounds if !isnan(dg["FRAC"][2, i]) && Δt < dg["TAUS"][i]

                # Fraction of atomic hydrogen according to our model
                fa[i] = dg["FRAC"][2, i]

            else

                # When there is no data from our model, use the fraction of neutral hydrogen from "NH  "
                fa[i] = dg["NH  "][i] / (dg["NHP "][i] + dg["NH  "][i])

            end

        end

    elseif !any(isempty, [dg["PRES"], dg["COLM"]])

        relative_pressure = uconvert.(Unitful.NoUnits, dg["PRES"] ./ P0).^ALPHA_BLITZ

        # Fraction of molecular hydrogen according to the pressure relation in Blitz et al. (2006)
        fp = 1.0 ./ (1.0 .+ relative_pressure)

        # Fraction of cold gas from Arepo
        fn = uconvert.(Unitful.NoUnits, dg["COLM"] ./ dg["MASS"])

        fa = 1.0 .- fp .* fn

    else

        return Unitful.Mass[]

    end

    return fa .* dg["MASS"]

end

"""
    computeMolecularMass(data_dict::Dict; <keyword arguments>)::Vector{<:Unitful.Mass}

Compute the mass of molecular hydrogen in every gas cell/particle.

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
  - `normalize::Bool=true`: If the output will be normalize to eliminate the stellar fraction of the gas cells. Only relevant for simulation with our SF routine, and for cells that have entered it at least once.

# Returns

  - The mass of molecular hydrogen in every gas cell/particle.
"""
function computeMolecularMass(data_dict::Dict; normalize::Bool=true)::Vector{<:Unitful.Mass}

    dg = data_dict[:gas]

    !isempty(dg["MASS"]) || return Unitful.Mass[]

    if "FRAC" ∈ keys(dg) && !isempty(dg["FRAC"])

        fm = Vector{Float64}(undef, length(dg["MASS"]))

        @inbounds for i in eachindex(fm)

            # Compute how much time has pass since the last time the cell/particle entered the SF routine
            Δt = data_dict[:snap_data].physical_time - dg["CTIM"][i]

            @inbounds if !isnan(dg["FRAC"][3, i]) && Δt < dg["TAUS"][i]

                # Fraction of molecular hydrogen according to our model
                @inbounds if normalize
                    fm[i] = dg["FRAC"][3, i] + dg["FRAC"][4, i]
                else
                    fm[i] = dg["FRAC"][3, i]
                end

            else

                # When there is no data from our model, assume no molecular hydrogen
                fm[i] = 0.0

            end

        end

    elseif !any(isempty, [dg["PRES"], dg["COLM"]])

        relative_pressure = uconvert.(Unitful.NoUnits, dg["PRES"] ./ P0).^ALPHA_BLITZ

        # Fraction of molecular hydrogen according to the pressure relation in Blitz et al. (2006)
        fp = 1.0 ./ (1.0 .+ relative_pressure)

        # Fraction of cold gas from Arepo
        fn = uconvert.(Unitful.NoUnits, dg["COLM"] ./ dg["MASS"])

        fm = fp .* fn

    else

        return Unitful.Mass[]

    end

    return fm .* dg["MASS"]

end

"""
    computePressureMolecularMass(data_dict::Dict)::Vector{<:Unitful.Mass}

Compute the mass molecular hydrogen in every gas cell/particle using the pressure relation from Blitz et al. (2006).

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

  - The mass of molecular hydrogen in every gas cell/particle.

# References

L. Blitz et al. (2006). *The Role of Pressure in GMC Formation II: The H2-Pressure Relation*. The Astrophysical Journal, **650(2)**, 933. [doi:10.1086/505417](https://doi.org/10.1086/505417)
"""
function computePressureMolecularMass(data_dict::Dict)::Vector{<:Unitful.Mass}

    dg = data_dict[:gas]

    !any(isempty, [dg["MASS"], dg["PRES"], dg["COLM"]]) || return Unitful.Mass[]

    relative_pressure = uconvert.(Unitful.NoUnits, dg["PRES"] ./ P0).^ALPHA_BLITZ

    # Fraction of molecular hydrogen according to the pressure relation in Blitz et al. (2006)
    fp = 1.0 ./ (1.0 .+ relative_pressure)

    # Fraction of cold gas from Arepo
    fn = uconvert.(Unitful.NoUnits, dg["COLM"] ./ dg["MASS"])

    return fp .* fn .* dg["MASS"]

end

"""
    computeStellarGasMass(data_dict::Dict)::Vector{<:Unitful.Mass}

Compute the "stellar mass" in every gas cell/particle.

!!! note

    It can be a non 0 value only for simulations with our SF routine.

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

  - The "stellar mass" in every gas cell/particle.
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

    return fm .* dg["MASS"]

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
