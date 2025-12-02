####################################################################################################
# Compute mass related quantities
####################################################################################################

################
# Base function
################

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

    if isempty(density)

        logging[] && @warn("computeClumpingFactor: `density` is empty, so I will return NaN")

        return NaN

    end

    μ, var = mean_and_var(density)

    return 1.0 + uconvert(Unitful.NoUnits, var / μ^2)

end

@doc raw"""
    computeEfficiencyFF(
        densities::Vector{<:Unitful.Density},
        masses::Vector{<:Unitful.Mass},
        sfrs::Vector{<:Unitful.MassFlow},
    )::Vector{Float64}

Compute the star formation efficiency per free-fall time, according to the definition in eq. 1 of Krumholz et al. (2012),

```math
\epsilon_\mathrm{ff} = \frac{t_\mathrm{ff}}{t_\mathrm{dep}} \, .
```
where

```math
t_\mathrm{ff} = \sqrt{\frac{3 \, \pi}{32 \, G \, \rho}} \, ,
```
is the free-fall time, and

```math
t_\mathrm{dep} = \frac{M}{\dot{M}_\star} \, ,
```
is the depletion time. $M$ and $\rho$ are the mass and density of the target gas phase, and $\dot{M}_\star$ is the SFR.

# Arguments

  - `densities::Vector{<:Unitful.Density}`: Densities of the cells.
  - `masses::Vector{<:Unitful.Mass}`: Masses of the cells.
  - `sfrs::Vector{<:Unitful.MassFlow}`: The SFR associated to each cell.

# Returns

  - The star formation efficiency per free-fall time.

# References

M. R. Krumholz et al. (2012). *A UNIVERSAL, LOCAL STAR FORMATION LAW IN GALACTIC CLOUDS, NEARBY GALAXIES, HIGH-REDSHIFT DISKS, AND STARBURSTS*. The Astrophysical Journal, **745(1)**, 69. [doi:10.1088/0004-637X/745/1/69](https://doi.org/10.1088/0004-637X/745/1/69)
"""
function computeEfficiencyFF(
    densities::Vector{<:Unitful.Density},
    masses::Vector{<:Unitful.Mass},
    sfrs::Vector{<:Unitful.MassFlow},
)::Vector{Float64}

    if any(isempty, [densities, masses, sfrs])

        (
            logging[] &&
            @warn("computeEfficiencyFF: There is missing data, so I will return an empty array")
        )

        return Float64[]

    end

    (
        allequal(length, [densities, masses, sfrs]) ||
        throw(ArgumentError("computeEfficiencyFF: The lengths of `densities`, `masses` and `sfrs` \
        must be equal"))
    )

    ϵff = Vector{Float64}(undef, length(densities))

    for (i, (ρ, m, sfr)) in enumerate(zip(densities, masses, sfrs))

        if iszero(m) || iszero(sfr) || ρ < THRESHOLD_DENSITY

            ϵff[i] = NaN

        else

            # Compute the free-fall time
            tff = sqrt(3π / (32 * Unitful.G * ρ))

            # Compute the depletion time
            tdep = m / sfr

            # Compute the star formation efficiency per free-fall time
            ϵff[i] = uconvert(Unitful.NoUnits, tff / tdep)

        end

    end

    return ϵff

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
        0 < percent <= 100 ||
        throw(ArgumentError("computeMassRadius: The argument `percent` must be between 0 and 100, \
        but I got $(percent)"))
    )

    # Check for missing data
    if any(isempty, [positions, masses])
        (
            logging[] &&
            @warn("computeMassRadius: `positions` or `masses` are empty, so I will return 0s")
        )
        return zero(typeof(1.0u"kpc"))
    end

    # Compute the mass limit
    mass_limit = sum(masses) * (percent / 100.0)

    # Compute the radial distance of each cell/particle to the origin
    radial_distances = colwise(Euclidean(), positions, zeros(eltype(positions), size(positions, 1)))

    sort_idxs = sortperm(radial_distances)

    # Find the mass radius
    accu_mass = 0.0u"Msun"
    target_idx = 0
    for mass in masses[sort_idxs]
        accu_mass += mass
        target_idx += 1
        accu_mass < mass_limit || break
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
        0.0 < percent <= 100.0 ||
        throw(ArgumentError("computeMassHeight: The argument `percent` must be between 0 and 100, \
        but I got $(percent)"))
    )

    # Check for missing data
    if any(isempty, [positions, masses])
        (
            logging[] &&
            @warn("computeMassHeight: `positions` or `masses` are empty, so I will return 0s")
        )
        return zero(typeof(1.0u"kpc"))
    end

    # Compute the mass limit
    mass_limit = sum(masses) * (percent / 100.0)

    # Compute the vertical distance of each cell/particle to the xy plane
    heights = abs.(vec(positions[3, :]))

    sort_idxs = sortperm(heights)

    # Find the mass height
    accu_mass = 0.0u"Msun"
    target_idx = 0
    for mass in masses[sort_idxs]
        accu_mass += mass
        target_idx += 1
        accu_mass < mass_limit || break
    end

    return heights[sort_idxs[target_idx]] * 2.0

end

"""
    computeMassQty(
        quantity::Vector{<:Number},
        masses::Vector{<:Unitful.Mass};
        <keyword arguments>
    )::Number

Compute the minimum value of `quantity` that "contains" `percent`% of the total mass.

# Arguments

  - `quantity::Vector{<:Number}`: Target quantity.
  - `masses::Vector{<:Unitful.Mass}`: Masses of the cells/particles.
  - `percent::Float64=90.0`: Target percentage of the total mass.

# Returns

  - The minimum value of `quantity` that "contains" `percent`% of the total mass.
"""
function computeMassQty(
    quantity::Vector{<:Number},
    masses::Vector{<:Unitful.Mass};
    percent::Float64=90.0,
)::Number

    (
        0.0 < percent <= 100.0 ||
        throw(ArgumentError("computeMassQty: The argument `percent` must be between 0 and 100, \
        but I got $(percent)"))
    )

    # Check for missing data
    if any(isempty, [quantity, masses])
        (
            logging[] &&
            @warn("computeMassQty: `quantity` or `masses` are empty, so I will return 0s")
        )
        return zero(eltype(quantity))
    end

    # Compute the mass limit
    mass_limit = sum(masses) * (percent / 100.0)

    sort_idxs = sortperm(quantity)

    accu_mass = 0.0u"Msun"
    target_idx = 0
    for mass in masses[sort_idxs]
        accu_mass += mass
        target_idx += 1
        accu_mass < mass_limit || break
    end

    return quantity[sort_idxs[target_idx]]

end

"""
    computeFractionWithin(
        quantity::Vector{<:Number},
        masses::Vector{<:Unitful.Mass},
        qty_limits::Tuple{<:Number,<:Number},
    )::Float64

Compute the fraction of the total mass "contained" within a given range of `quantity`.

# Arguments

  - `quantity::Vector{<:Number}`: Target quantity.
  - `masses::Vector{<:Unitful.Mass}`: Masses of the cells/particles.
  - `qty_limits::Tuple{<:Number,<:Number}`: Range of the target quantity.

# Returns

  - The fraction of the total mass "contained" within a given range of `quantity`.
"""
function computeFractionWithin(
    quantity::Vector{<:Number},
    masses::Vector{<:Unitful.Mass},
    qty_limits::Tuple{<:Number,<:Number},
)::Float64

    # Check for missing data
    if any(isempty, [quantity, masses])
        (
            logging[] &&
            @warn("computeFractionWithin: `quantity` or `masses` are empty, so I will return 0")
        )
        return 0.0
    end

    # Find the indices of all the cells/particles with
    # `qty_limits[1]` <= `quantity` <= `qty_limits[2]`
    idxs = map(x -> qty_limits[1] <= x <= qty_limits[2], quantity)

    # Total mass
    M = sum(masses)

    # Mass within `qty_limits[1]` <= `quantity` <= `qty_limits[2]`
    Mq = sum(masses[idxs]; init=0.0u"Msun")

    return uconvert(Unitful.NoUnits, Mq / M)

end

"""
    computeFraction(data_dict::Dict, component::Symbol)::Vector{Float64}

Compute the fraction of a given `component` in each cell/particle.

# Arguments

  - `data_dict::Dict`: Data dictionary (see [`makeDataDict`](@ref) for the canonical description).
    This function requires the following blocks to be present, depending on the value of `component`:

      + If `component` == :Z_stellar
          * `:stellar` => ["GZ2 "]
      + If `component` ∈ [:gas, :hydrogen, :helium]:
          * `:gas` => ["MASS"]
      + If `component` == :Z_gas:
          * `:gas` => ["MASS", "GZ  "]
      + If `component` ∈ [:ionized, :neutral]:
          * `:gas` => ["MASS", "NH  ", "NHP "]
      + If `component` ∈ [:br_atomic, :br_molecular]:
          * `:gas` => ["MASS", "NH  ", "NHP ", "PRES"]
      + If `component` ∈ [:ode_ionized, :ode_atomic, :ode_metals, :ode_dust, :ode_neutral]:
          * `:gas` => ["MASS", "NH  ", "NHP ", "FRAC", "RHO ", "GZ  "]
      + If `component` ∈ [:ode_molecular, :ode_stellar]:
          * `:gas` => ["MASS", "FRAC", "RHO "]
  - `component::Symbol`: Target component. It can only be one of the elements of [`COMPONENTS`](@ref).

# Returns

  - The fraction of `component` in each cell/particle.

# References

L. Blitz et al. (2006). *The Role of Pressure in GMC Formation II: The H2-Pressure Relation*. The Astrophysical Journal, **650(2)**, 933. [doi:10.1086/505417](https://doi.org/10.1086/505417)
"""
function computeFraction(data_dict::Dict, component::Symbol)::Vector{Float64}

    if component ∉ COMPONENTS || component ∈ [:stellar, :dark_matter, :black_hole]
        throw(ArgumentError("computeFraction: `component` can only be one of the elements of \
        `COMPONENTS` (see `./src/constants/globals.jl`), but I got :$(component)"))
    end

    ################################################################################################
    # Stellar metallicity
    ################################################################################################

    if component == :Z_stellar

        Z = data_dict[:stellar]["GZ2 "]

        if isempty(Z)

            logging[] && @warn("computeFraction: I could not compute the stellar metallicity")

            fractions = Float64[]

        else

            fractions = setPositive(Z)

        end

    else

        dg   = data_dict[:gas]
        mass = dg["MASS"]

        # If there is no gas, return an empty array
        if isempty(mass)

            logging[] && @warn("computeFraction: There is no data for the gas cells!")

            return Float64[]

        end

        # Compute the number of gas cells
        n_cells = length(mass)

    end

    ################################################################################################
    # Total gas
    ################################################################################################

    if component == :gas

        fractions = ones(n_cells)

    ################################################################################################
    # Hydrogen
    ################################################################################################

    elseif component == :hydrogen

        fractions = fill(HYDROGEN_MASSFRAC, n_cells)

    ################################################################################################
    # Helium
    ################################################################################################

    elseif component == :helium

        fractions = fill(1.0 - HYDROGEN_MASSFRAC, n_cells)

    ################################################################################################
    # Gas metallicity
    ################################################################################################

    elseif component == :Z_gas

        Z = dg["GZ  "]

        if isempty(Z)

            logging[] && @warn("computeFraction: I could not compute the gas metallicity")

            fractions = Float64[]

        else

            fractions = setPositive(Z)

        end

    ################################################################################################
    # Ionized gas (using the Arepo data)
    ################################################################################################

    elseif component == :ionized

        nh  = dg["NH  "]
        nhp = dg["NHP "]

        if any(isempty, [nh, nhp])

            logging[] && @warn("computeFraction: I could not compute the ionized fraction")

            fractions = Float64[]

        else

            fractions = @. nhp / (nhp + nh)

        end

    ################################################################################################
    # Neutral gas (using the Arepo data)
    ################################################################################################

    elseif component == :neutral

        nh  = dg["NH  "]
        nhp = dg["NHP "]

        if any(isempty, [nh, nhp])

            logging[] && @warn("computeFraction: I could not compute the neutral fraction")

            fractions = Float64[]

        else

            fractions = @. nh / (nhp + nh)

        end

    ################################################################################################
    # Atomic gas (using the Blitz et al. (2006) relation)
    ################################################################################################

    elseif component == :br_atomic

        nh  = dg["NH  "]
        nhp = dg["NHP "]
        P   = dg["PRES"]

        if any(isempty, [nh, nhp, P])

            logging[] && @warn("computeFraction: I could not compute the BR atomic fraction")

            fractions = Float64[]

        else

            relative_pressure = @. uconvert(Unitful.NoUnits, P / P0)^ALPHA_BLITZ

            # Compute the fraction of neutral gas that is atomic hydrogen according
            # to the pressure relation in Blitz et al. (2006)
            fa = @. relative_pressure / (1.0 + relative_pressure)

            # Compute the fraction of neutral gas according to Arepo
            fn = @. nh / (nhp + nh)

            fractions = fa .* fn

        end

    ################################################################################################
    # Molecular gas (using the Blitz et al. (2006) relation)
    ################################################################################################

    elseif component == :br_molecular

        nh  = dg["NH  "]
        nhp = dg["NHP "]
        P   = dg["PRES"]

        if any(isempty, [nh, nhp, P])

            logging[] && @warn("computeFraction: I could not compute the BR molecular fraction")

            fractions = Float64[]

        else

            relative_pressure = @. uconvert(Unitful.NoUnits, P / P0)^ALPHA_BLITZ

            # Compute the fraction of neutral gas that is molecular hydrogen according
            # to the pressure relation in Blitz et al. (2006)
            fm = @. 1.0 / (1.0 + relative_pressure)

            # Compute the fraction of neutral gas according to Arepo
            fn = @. nh / (nhp + nh)

            fractions = fm .* fn

        end

    ################################################################################################
    # Ionized gas (according to our SF model)
    ################################################################################################

    elseif component == :ode_ionized

        nh   = dg["NH  "]
        nhp  = dg["NHP "]
        frac = dg["FRAC"]
        ρc   = dg["RHO "]
        Z    = dg["GZ  "]

        if any(isempty, [nh, nhp, frac, ρc, Z])

            logging[] && @warn("computeFraction: I could not compute the ODE ionized fraction")

            fractions = Float64[]

        else

            fi = view(frac, SFM_IDX[component], :)

            fractions = Vector{Float64}(undef, n_cells)

            Threads.@threads for i in eachindex(fractions)

                if !isnan(fi[i]) && ρc[i] >= THRESHOLD_DENSITY

                    # Fraction of ionized hydrogen according to our SF model
                    fractions[i] = fi[i]

                else

                    metallicity = setPositive(Z[i])

                    # When there is no data from the model or the density is below the SF threshold,
                    # use the initial condition of the model
                    fractions[i] = (1.0 - metallicity) * nhp[i] / (nhp[i] + nh[i])

                end

            end

        end

    ################################################################################################
    # Atomic gas (according to our SF model)
    ################################################################################################

    elseif component == :ode_atomic

        nh   = dg["NH  "]
        nhp  = dg["NHP "]
        frac = dg["FRAC"]
        ρc   = dg["RHO "]
        Z    = dg["GZ  "]

        if any(isempty, [nh, nhp, frac, ρc, Z])

            logging[] && @warn("computeFraction: I could not compute the ODE atomic fraction")

            fractions = Float64[]

        else

            fa = view(frac, SFM_IDX[component], :)

            fractions = Vector{Float64}(undef, n_cells)

            Threads.@threads for i in eachindex(fractions)

                if !isnan(fa[i]) && ρc[i] >= THRESHOLD_DENSITY

                    # Fraction of atomic hydrogen according to our SF model
                    fractions[i] = fa[i]

                else

                    metallicity = setPositive(Z[i])

                    # When there is no data from the model or the density is below the SF threshold,
                    # use the initial condition of the model
                    fractions[i] = (1.0 - metallicity) * nh[i] / (nhp[i] + nh[i])

                end

            end

        end

    ################################################################################################
    # Molecular gas (according to our SF model)
    ################################################################################################

    elseif component == :ode_molecular

        frac = dg["FRAC"]
        ρc   = dg["RHO "]

        if any(isempty, [frac, ρc])

            logging[] && @warn("computeFraction: I could not compute the ODE molecular fraction")

            fractions = Float64[]

        else

            fm = view(frac, SFM_IDX[component], :)
            fs = view(frac, SFM_IDX[:ode_stellar], :)

            fractions = Vector{Float64}(undef, n_cells)

            Threads.@threads for i in eachindex(fractions)

                if !isnan(fm[i]) && ρc[i] >= THRESHOLD_DENSITY

                    # Fraction of molecular hydrogen according to our SF model
                    fractions[i] = fm[i] + fs[i]

                else

                    # When there is no data from the model or the density is below the SF threshold,
                    # use the initial condition of the model
                    fractions[i] = 0.0

                end

            end

        end

    ################################################################################################
    # Stars (according to our SF model)
    ################################################################################################

    elseif component == :ode_stellar

        frac = dg["FRAC"]
        ρc   = dg["RHO "]

        if any(isempty, [frac, ρc])

            logging[] && @warn("computeFraction: I could not compute the ODE stellar fraction")

            fractions = Float64[]

        else

            fs = view(frac, SFM_IDX[component], :)

            fractions = Vector{Float64}(undef, n_cells)

            Threads.@threads for i in eachindex(fractions)

                if !isnan(fs[i]) && ρc[i] >= THRESHOLD_DENSITY

                    # Fraction of stars according to our SF model
                    fractions[i] = fs[i]

                else

                    # When there is no data from the model or the density is below the SF threshold,
                    # use the initial condition of the model
                    fractions[i] = 0.0

                end

            end

        end

    ################################################################################################
    # Metals (according to our SF model)
    ################################################################################################

    elseif component == :ode_metals

        nh   = dg["NH  "]
        nhp  = dg["NHP "]
        frac = dg["FRAC"]
        ρc   = dg["RHO "]
        Z    = dg["GZ  "]

        if any(isempty, [nh, nhp, frac, ρc, Z])

            logging[] && @warn("computeFraction: I could not compute the ODE metals fraction")

            fractions = Float64[]

        else

            fZ = view(frac, SFM_IDX[component], :)

            fractions = Vector{Float64}(undef, n_cells)

            Threads.@threads for i in eachindex(fractions)

                if !isnan(fZ[i]) && ρc[i] >= THRESHOLD_DENSITY

                    # Fraction of metals according to our SF model
                    fractions[i] = fZ[i]

                else

                    metallicity = setPositive(Z[i])
                    fa = (1 - metallicity) * nh[i] / (nhp[i] + nh[i])

                    # When there is no data from the model or the density is below the SF threshold,
                    # use the initial condition of the model
                    fractions[i] = metallicity * (1.0 - Cxd * fa)

                end

            end

        end

    ################################################################################################
    # Dust (according to our SF model)
    ################################################################################################

    elseif component == :ode_dust

        nh   = dg["NH  "]
        nhp  = dg["NHP "]
        frac = dg["FRAC"]
        ρc   = dg["RHO "]
        Z    = dg["GZ  "]

        if any(isempty, [nh, nhp, frac, ρc, Z])

            logging[] && @warn("computeFraction: I could not compute the ODE dust fraction")

            fractions = Float64[]

        else

            fd = view(frac, SFM_IDX[component], :)

            fractions = Vector{Float64}(undef, n_cells)

            Threads.@threads for i in eachindex(fractions)

                if !isnan(fd[i]) && ρc[i] >= THRESHOLD_DENSITY

                    # Fraction of dust according to our SF model
                    fractions[i] = fd[i]

                else

                    metallicity = setPositive(Z[i])
                    fa = (1.0 - metallicity) * nh[i] / (nhp[i] + nh[i])

                    # When there is no data from the model or the density is below the SF threshold,
                    # use the initial condition of the model
                    fractions[i] = metallicity * Cxd * fa

                end

            end

        end

    ################################################################################################
    # Neutral gas (according to our SF model)
    ################################################################################################

    elseif component == :ode_neutral

        nh   = dg["NH  "]
        nhp  = dg["NHP "]
        frac = dg["FRAC"]
        ρc   = dg["RHO "]
        Z    = dg["GZ  "]

        if any(isempty, [nh, nhp, frac, ρc, Z])

            logging[] && @warn("computeFraction: I could not compute the ODE neutral fraction")

            fractions = Float64[]

        else

            fa = view(frac, SFM_IDX[:ode_atomic], :)
            fm = view(frac, SFM_IDX[:ode_molecular], :)
            fs = view(frac, SFM_IDX[:ode_stellar], :)

            fractions = Vector{Float64}(undef, n_cells)

            Threads.@threads for i in eachindex(fractions)

                if !isnan(fa[i]) && ρc[i] >= THRESHOLD_DENSITY

                    # Fraction of neutral hydrogen according to our SF model
                    fractions[i] = fa[i] + fm[i] + fs[i]

                else

                    metallicity = setPositive(Z[i])

                    # When there is no data from the model or the density is below the SF threshold,
                    # use the initial condition of the model
                    fractions[i] = (1.0 - metallicity) * nh[i] / (nhp[i] + nh[i])

                end

            end

        end

    end

    return fractions

end

###################
# Derive functions
###################

@doc raw"""
    computeClumpingFactor(data_dict::Dict, component::Symbol)::Float64

Compute the clumping factor,

```math
C_\rho = \frac{\langle \rho^2 \rangle}{\langle \rho \rangle^2} \, .
```

# Arguments

  - `data_dict::Dict`: Data dictionary (see [`makeDataDict`](@ref) for the canonical description).
    This function requires the following blocks to be present, depending on the value of `component`:

      + If `component` ∈ [:gas, :hydrogen, :helium]:
          * `:gas` => ["MASS", "RHO "]
      + If `component` == :Z_gas:
          * `:gas` => ["MASS", "GZ  ", "RHO "]
      + If `component` ∈ [:ionized, :neutral]:
          * `:gas` => ["MASS", "NH  ", "NHP ", "RHO "]
      + If `component` ∈ [:br_atomic, :br_molecular]:
          * `:gas` => ["MASS", "NH  ", "NHP ", "PRES", "RHO "]
      + If `component` ∈ [:ode_ionized, :ode_atomic, :ode_metals, :ode_dust, :ode_neutral]:
          * `:gas` => ["MASS", "NH  ", "NHP ", "FRAC", "RHO ", "GZ  "]
      + If `component` ∈ [:ode_molecular, :ode_stellar]:
          * `:gas` => ["MASS", "FRAC", "RHO "]
  - `component::Symbol`: Target component. It can only be one of the gas elements of [`COMPONENTS`](@ref).

# Returns

  - The clumping factor.
"""
function computeClumpingFactor(data_dict::Dict, component::Symbol)::Float64

    if component ∉ COMPONENTS || component ∈ [:stellar, :dark_matter, :black_hole, :Z_stellar]
        throw(ArgumentError("computeMassDensity: `component` can only be one of the gas elements \
        of `COMPONENTS` (see `./src/constants/globals.jl`), but I got :$(component)"))
    end

    ρ = computeMassDensity(data_dict, component)

    return computeClumpingFactor(ρ)

end

@doc raw"""
    computeEfficiencyFF(data_dict::Dict, component::Symbol)::Vector{Float64}

Compute the star formation efficiency per free-fall time, according to the definition in eq. 1 of Krumholz et al. (2012),

```math
\epsilon_\mathrm{ff} = \frac{t_\mathrm{ff}}{t_\mathrm{dep}} \, .
```
where

```math
t_\mathrm{ff} = \sqrt{\frac{3 \, \pi}{32 \, G \, \rho}} \, ,
```
is the free-fall time, and

```math
t_\mathrm{dep} = \frac{M}{\dot{M}_\star} \, ,
```
is the depletion time. $M$ and $\rho$ are the mass and density of the target gas phase, and $\dot{M}_\star$ is the SFR.

# Arguments

  - `data_dict::Dict`: Data dictionary (see [`makeDataDict`](@ref) for the canonical description).
    This function requires the following blocks to be present, depending on the value of `component`:

      + If `component == :stellar
          * `:stellar` => ["RHOC", "GMAS", "GSFR"]
      + If `component` ∈ [:gas, :hydrogen, :helium]:
          * `:gas` => ["SFR ", "MASS", "RHO "]
      + If `component` == :Z_gas:
          * `:gas` => ["SFR ", "MASS", "GZ  ", "RHO "]
      + If `component` ∈ [:ionized, :neutral]:
          * `:gas` => ["SFR ", "MASS", "NH  ", "NHP ", "RHO "]
      + If `component` ∈ [:br_atomic, :br_molecular]:
          * `:gas` => ["SFR ", "MASS", "NH  ", "NHP ", "PRES", "RHO "]
      + If `component` ∈ [:ode_ionized, :ode_atomic, :ode_metals, :ode_dust, :ode_neutral]:
          * `:gas` => ["SFR ", "MASS", "NH  ", "NHP ", "FRAC", "RHO ", "GZ  "]
      + If `component` ∈ [:ode_molecular, :ode_stellar]:
          * `:gas` => ["SFR ", "MASS", "FRAC", "RHO "]
  - `component::Symbol`: Target component. It can only be one of the gas elements of [`COMPONENTS`](@ref).

# Returns

  - The star formation efficiency per free-fall time.

# References

M. R. Krumholz et al. (2012). *A UNIVERSAL, LOCAL STAR FORMATION LAW IN GALACTIC CLOUDS, NEARBY GALAXIES, HIGH-REDSHIFT DISKS, AND STARBURSTS*. The Astrophysical Journal, **745(1)**, 69. [doi:10.1088/0004-637X/745/1/69](https://doi.org/10.1088/0004-637X/745/1/69)
"""
function computeEfficiencyFF(data_dict::Dict, component::Symbol)::Vector{Float64}

    if component ∉ COMPONENTS || component ∈ [:dark_matter, :black_hole, :Z_stellar]
        throw(ArgumentError("computeMassDensity: `component` can only be one of the gas elements \
        of `COMPONENTS` or :stellar (see `./src/constants/globals.jl`), but I got :$(component)"))
    end

    if component == :stellar

        # Compute the ϵff of the progenitor
        densities = data_dict[:stellar]["RHOC"] * 1.0u"mp"
        masses    = data_dict[:stellar]["GMAS"]
        sfrs      = data_dict[:stellar]["GSFR"]

    else

        densities = computeMassDensity(data_dict, component)
        masses    = computeMass(data_dict, component)
        sfrs      = data_dict[:gas]["SFR "]

    end

    return computeEfficiencyFF(densities, masses, sfrs)

end

"""
    computeMass(data_dict::Dict, component::Symbol)::Vector{<:Unitful.Mass}

Compute the mass in each cell/particle of a given `component`.

# Arguments

  - `data_dict::Dict`: Data dictionary (see [`makeDataDict`](@ref) for the canonical description).
    This function requires the following blocks to be present, depending on the value of `component`:

      + If `component` ∈ [:stellar, :dark_matter, :black_hole, :gas]:
          * `component` => ["MASS"]
      + If `component` == :Z_stellar
          * `:stellar` => ["MASS", "GZ2 "]
      + If `component` ∈ [:hydrogen, :helium]:
          * `:gas` => ["MASS"]
      + If `component` == :Z_gas:
          * `:gas` => ["MASS", "GZ  "]
      + If `component` ∈ [:ionized, :neutral]:
          * `:gas` => ["MASS", "NH  ", "NHP "]
      + If `component` ∈ [:br_atomic, :br_molecular]:
          * `:gas` => ["MASS", "NH  ", "NHP ", "PRES"]
      + If `component` ∈ [:ode_ionized, :ode_atomic, :ode_metals, :ode_dust, :ode_neutral]:
          * `:gas` => ["MASS", "NH  ", "NHP ", "FRAC", "RHO ", "GZ  "]
      + If `component` ∈ [:ode_molecular, :ode_stellar]:
          * `:gas` => ["MASS", "FRAC", "RHO "]
  - `component::Symbol`: Target component. It can only be one of the elements of [`COMPONENTS`](@ref).

# Returns

  - The mass of `component` in each cell/particle.

# References

L. Blitz et al. (2006). *The Role of Pressure in GMC Formation II: The H2-Pressure Relation*. The Astrophysical Journal, **650(2)**, 933. [doi:10.1086/505417](https://doi.org/10.1086/505417)
"""
function computeMass(data_dict::Dict, component::Symbol)::Vector{<:Unitful.Mass}

    if component ∉ COMPONENTS
        throw(ArgumentError("computeMass: `component` can only be one of the elements of \
        `COMPONENTS` (see `./src/constants/globals.jl`), but I got :$(component)"))
    end

    if component ∈ [:gas, :stellar, :dark_matter, :black_hole]

        masses = data_dict[component]["MASS"]

    else

        fractions = computeFraction(data_dict, component)

        if isempty(fractions)

            masses = Unitful.Mass[]

        else

            if component == :Z_stellar
                type = :stellar
            else
                type = :gas
            end

            masses = data_dict[type]["MASS"] .* fractions

        end

    end

    if isempty(masses)

        logging[] && @warn("computeMass: I could not compute the masses of :$(component)")

        return Unitful.Mass[]

    end

    return masses

end

"""
    computeMassDensity(data_dict::Dict, component::Symbol)::Vector{<:Unitful.Density}

Compute the mass density of a given gas `component` for each cell.

# Arguments

  - `data_dict::Dict`: Data dictionary (see [`makeDataDict`](@ref) for the canonical description).
    This function requires the following blocks to be present, depending on the value of `component`:

      + If `component` ∈ [:gas, :hydrogen, :helium]:
          * `:gas` => ["MASS", "RHO "]
      + If `component` == :Z_gas:
          * `:gas` => ["MASS", "GZ  ", "RHO "]
      + If `component` ∈ [:ionized, :neutral]:
          * `:gas` => ["MASS", "NH  ", "NHP ", "RHO "]
      + If `component` ∈ [:br_atomic, :br_molecular]:
          * `:gas` => ["MASS", "NH  ", "NHP ", "PRES", "RHO "]
      + If `component` ∈ [:ode_ionized, :ode_atomic, :ode_metals, :ode_dust, :ode_neutral]:
          * `:gas` => ["MASS", "NH  ", "NHP ", "FRAC", "RHO ", "GZ  "]
      + If `component` ∈ [:ode_molecular, :ode_stellar]:
          * `:gas` => ["MASS", "FRAC", "RHO "]
  - `component::Symbol`: Target component. It can only be one of the gas elements of [`COMPONENTS`](@ref).

# Returns

  - The mass density of `component` for each cell.

# References

L. Blitz et al. (2006). *The Role of Pressure in GMC Formation II: The H2-Pressure Relation*. The Astrophysical Journal, **650(2)**, 933. [doi:10.1086/505417](https://doi.org/10.1086/505417)
"""
function computeMassDensity(data_dict::Dict, component::Symbol)::Vector{<:Unitful.Density}

    if component ∉ COMPONENTS || component ∈ [:stellar, :dark_matter, :black_hole, :Z_stellar]
        throw(ArgumentError("computeMassDensity: `component` can only be one of the gas elements \
        of `COMPONENTS` (see `./src/constants/globals.jl`), but I got :$(component)"))
    end

    if component == :gas

        ρ = data_dict[:gas]["RHO "]

    else

        fractions = computeFraction(data_dict, component)

        if isempty(fractions)

            ρ = Unitful.Density[]

        else

            ρ = data_dict[:gas]["RHO "] .* fractions

        end

    end

    if isempty(ρ)

        (
            logging[] &&
            @warn("computeMassDensity: I could not compute the mass densities of :$(component)")
        )

        return Unitful.Density[]

    end

    return ρ

end

"""
    computeNumberDensity(data_dict::Dict, component::Symbol)::Vector{<:NumberDensity}

Compute the number density of a given gas `component` for each cell.

!!! note

    What number density means changes depending on the `component`, see the comments in the code.

# Arguments

  - `data_dict::Dict`: Data dictionary (see [`makeDataDict`](@ref) for the canonical description).
    This function requires the following blocks to be present, depending on the value of `component`:

      + If `component` ∈ [:gas, :hydrogen, :helium]:
          * `:gas` => ["MASS", "RHO "]
      + If `component` == :Z_gas:
          * `:gas` => ["MASS", "GZ  ", "RHO "]
      + If `component` ∈ [:ionized, :neutral]:
          * `:gas` => ["MASS", "NH  ", "NHP ", "RHO "]
      + If `component` ∈ [:br_atomic, :br_molecular]:
          * `:gas` => ["MASS", "NH  ", "NHP ", "PRES", "RHO "]
      + If `component` ∈ [:ode_ionized, :ode_atomic, :ode_metals, :ode_dust, :ode_neutral]:
          * `:gas` => ["MASS", "NH  ", "NHP ", "FRAC", "RHO ", "GZ  "]
      + If `component` ∈ [:ode_molecular, :ode_stellar]:
          * `:gas` => ["MASS", "FRAC", "RHO "]
  - `component::Symbol`: Target component. It can only be one of the gas elements of [`COMPONENTS`](@ref).

# Returns

  - The number density of `component` for each cell.

# References

L. Blitz et al. (2006). *The Role of Pressure in GMC Formation II: The H2-Pressure Relation*. The Astrophysical Journal, **650(2)**, 933. [doi:10.1086/505417](https://doi.org/10.1086/505417)
"""
function computeNumberDensity(data_dict::Dict, component::Symbol)::Vector{<:NumberDensity}

    if component ∈ [:helium, :br_molecular, :ode_molecular]

        ############################################################################################
        # Number density as the amount of elements (atoms/molecules) per unit volume
        ############################################################################################

        ρ = computeMassDensity(data_dict, component)

        n = ρ / (2.0 * Unitful.mp)

    elseif component == :ode_neutral

        ############################################################################################
        # Number density as the amount of elements (atoms/molecules) per unit volume
        ############################################################################################

        ρm = computeMassDensity(data_dict, :ode_molecular)
        ρa = computeMassDensity(data_dict, :ode_atomic)

        n = (ρm / (2.0 * Unitful.mp)) .+ (ρa / Unitful.mp)

    else

        ############################################################################################
        # Number density as the mass density in units of proton mass per unit volume
        ############################################################################################

        ρ = computeMassDensity(data_dict, component)

        n = ρ / Unitful.mp

    end

    if isempty(n)

        (
            logging[] &&
            @warn("computeNumberDensity: I could not compute the number densities of :$(component)")
        )

        return NumberDensity[]

    end

    return n

end

"""
    computeNumber(data_dict::Dict, component::Symbol)::Vector{Int64}

Compute the number of a given `component`.

!!! note

    What number means changes depending on the `component`, see the comments in the code.

# Arguments

  - `data_dict::Dict`: Data dictionary (see [`makeDataDict`](@ref) for the canonical description).
    This function requires the following blocks to be present, depending on the value of `component`:

      + If `component` ∈ [:stellar, :dark_matter, :black_hole, :gas]:
          * `component` => ["MASS"]
      + If `component` ∈ [:Z_stellar]:
          * `:stellar` => ["MASS"]
      + If `component` ∈ [:hydrogen, :helium]:
          * `:gas` => ["MASS"]
      + If `component` == :Z_gas:
          * `:gas` => ["MASS", "GZ  "]
      + If `component` ∈ [:ionized, :neutral]:
          * `:gas` => ["MASS", "NH  ", "NHP "]
      + If `component` ∈ [:br_atomic, :br_molecular]:
          * `:gas` => ["MASS", "NH  ", "NHP ", "PRES"]
      + If `component` ∈ [:ode_ionized, :ode_atomic, :ode_metals, :ode_dust, :ode_neutral]:
          * `:gas` => ["MASS", "NH  ", "NHP ", "FRAC", "RHO ", "GZ  "]
      + If `component` ∈ [:ode_molecular, :ode_stellar]:
          * `:gas` => ["MASS", "FRAC", "RHO "]
  - `component::Symbol`: Target component. It can only be one of the elements of [`COMPONENTS`](@ref).

# Returns

  - The number of `component`.

# References

L. Blitz et al. (2006). *The Role of Pressure in GMC Formation II: The H2-Pressure Relation*. The Astrophysical Journal, **650(2)**, 933. [doi:10.1086/505417](https://doi.org/10.1086/505417)
"""
function computeNumber(data_dict::Dict, component::Symbol)::Vector{Int64}

    if component ∉ COMPONENTS
        throw(ArgumentError("computeNumber: `component` can only be one of the elements of \
        `COMPONENTS` (see `./src/constants/globals.jl`), but I got :$(component)"))
    end

    if component ∈ [:stellar, :dark_matter, :black_hole. :gas]

        ############################################################################################
        # Number as the amount of particles
        ############################################################################################

        N = [lenght(data_dict[component]["MASS"])]

    elseif component == :Z_stellar

        ############################################################################################
        # Number as the amount of particles
        ############################################################################################

        N = [lenght(data_dict[:stellar]["MASS"])]

    elseif component ∈ [:helium, :br_molecular, :ode_molecular]

        ############################################################################################
        # Number as the amount of elements (atoms/molecules) in each cell
        ############################################################################################

        M = computeMass(data_dict, component)

        N = M / (2.0 * Unitful.mp)

    elseif component == :ode_neutral

        ############################################################################################
        # Number as the amount of elements (atoms/molecules) in each cell
        ############################################################################################

        Mm = computeMassDensity(data_dict, :ode_molecular)
        Ma = computeMassDensity(data_dict, :ode_atomic)

        N = (Mm / (2.0 * Unitful.mp)) + (Ma / Unitful.mp)

    else

        ############################################################################################
        # Number as the mass in units of proton mass
        ############################################################################################

        M = computeMass(data_dict, component)

        N = M / Unitful.mp

    end

    if isempty(N)

        logging[] && @warn("computeNumber: I could not compute the numbers of :$(component)")

        return Int64[]

    end

    return ustrip.(Unitful.NoUnits, N)

end

"""
    computeVirialAccretion(
        present_dd::Dict,
        past_dd::Dict;
        <keyword arguments>
    )::NTuple{3,Unitful.Mass}

Compute the inflow, outflow, or net gain of mass for a given halo virial radius (``R_{200}``), between two snapshots.

# Arguments

  - `present_dd::Dict`: Data dictionary, for the present snapshot (see [`makeDataDict`](@ref) for the canonical description).
    This function requires the following blocks to be present for every cell/particle that you want to be taken into account:

      + `cell/particle type` => ["POS ", "ID  ", "MASS"]
      + `:group`             => ["G_R_Crit200", "G_Nsubs", "G_Pos"]
      + `:subhalo`           => ["S_Pos"]
      + `:tracer`            => ["PAID", "TRID"]
  - `past_dd::Dict`: Data dictionary, for the past snapshot (see [`makeDataDict`](@ref) for the canonical description).
    This function requires the following blocks to be present for every cell/particle that you want to be taken into account:

      + `cell/particle type` => ["POS ", "ID  ", "MASS"]
      + `:group`             => ["G_R_Crit200", "G_Nsubs", "G_Pos"]
      + `:subhalo`           => ["S_Pos"]
      + `:tracer`            => ["PAID", "TRID"]
  - `component::Symbol`: Component to compute the accreted mass for. The options are:

      + `:dark_matter` -> Dark matter.
      + `:black_hole`  -> Black holes.
      + `:gas`         -> Gas.
      + `:stellar`     -> Stars.
  - `halo_idx::Int=1`: Index of the target halo (FoF group). Starts at 1.
  - `tracers::Bool=false`: Whether to compute the accretion using tracer particles (true) or the actual component particles (false).

# Returns

  - A tuple with three elements:

      + The net increase in mass.
      + The inflow mass.
      + The outflow mass.
"""
function computeVirialAccretion(
    present_dd::Dict,
    past_dd::Dict,
    component::Symbol;
    halo_idx::Int=1,
    tracers::Bool=false,
)::NTuple{3,Unitful.Mass}

    (
        component ∈ [:dark_matter, :black_hole, :gas, :stellar] ||
        throw(ArgumentError("computeVirialAccretion: `component` can only be :dark_matter, \
        :black_hole, :gas or :stellar, but I got :$(component)"))
    )

    # Find the IDs of the component particles inside R200 in the present snapshot
    present_ids = idWithinR200(present_dd, component; halo_idx)

    # Find the IDs of the component particles inside R200 in the past snapshot
    past_ids = idWithinR200(past_dd, component; halo_idx)

    # Filter out IDs of component particles that are new to the present snapshot
    # i.e., filter out created particles
    filterExistIDs!(past_dd, component, present_ids)

    # Filter out IDs of component particles that are unique to the past snapshot
    # i.e., filter out destroyed particles
    filterExistIDs!(present_dd, component, past_ids)

    # Find the IDs of the component particles that are inside R200 now, but were outside R200 in the past
    inflow_ids = setdiff(present_ids, past_ids)

    # Find the IDs of the component particles that were inside R200 in the past, but are now outside R200
    outflow_ids = setdiff(past_ids, present_ids)

    # Compute the inflow mass
    inflow_mass = idMass!(present_dd, component, inflow_ids; tracers)

    # Compute the outflow mass
    outflow_mass = idMass!(present_dd, component, outflow_ids; tracers)

    # Compute the net mass
    net_mass_increase = inflow_mass - outflow_mass

    return net_mass_increase, inflow_mass, outflow_mass

end

"""
    computeDiskAccretion(
        present_dd::Dict,
        past_dd::Dict;
        <keyword arguments>
    )::NTuple{3,Unitful.Mass}

Compute the inflow, outflow, or net gain of mass for a given galactic disk, between two snapshots.

!!! note

    It is assumed that the center of the disk is the origin.

# Arguments

  - `present_dd::Dict`: Data dictionary, for the present snapshot (see [`makeDataDict`](@ref) for the canonical description).
    This function requires the following blocks to be present for every cell/particle that you want to be taken into account:

      + `cell/particle type` => ["POS ", "ID  ", "MASS"]
      + `:group`             => ["G_R_Crit200", "G_Nsubs", "G_Pos"]
      + `:subhalo`           => ["S_Pos"]
      + `:tracer`            => ["PAID", "TRID"]
  - `past_dd::Dict`: Data dictionary, for the past snapshot (see [`makeDataDict`](@ref) for the canonical description).
    This function requires the following blocks to be present for every cell/particle that you want to be taken into account:

      + `cell/particle type` => ["POS ", "ID  ", "MASS"]
      + `:group`             => ["G_R_Crit200", "G_Nsubs", "G_Pos"]
      + `:subhalo`           => ["S_Pos"]
      + `:tracer`            => ["PAID", "TRID"]
  - `component::Symbol`: Component to compute the accreted mass for. The options are:

      + `:dark_matter` -> Dark matter.
      + `:black_hole`  -> Black holes.
      + `:gas`         -> Gas.
      + `:stellar`     -> Stars.
  - `max_r::Unitful.Length=DISK_R`: Radius of the disk.
  - `max_z::Unitful.Length=5.0u"kpc"`: Half height of the disk.
  - `tracers::Bool=false`: Whether to compute the accretion using tracer particles (true) or the actual component particles (false).

# Returns

  - A tuple with three elements:

      + The net increase in mass.
      + The inflow mass.
      + The outflow mass.
"""
function computeDiskAccretion(
    present_dd::Dict,
    past_dd::Dict,
    component::Symbol;
    max_r::Unitful.Length=DISK_R,
    max_z::Unitful.Length=5.0u"kpc",
    tracers::Bool=false,
)::NTuple{3,Unitful.Mass}

    (
        component ∈ [:dark_matter, :black_hole, :gas, :stellar] ||
        throw(ArgumentError("computeDiskAccretion: `component` can only be :dark_matter, \
        :black_hole, :gas or :stellar, but I got :$(component)"))
    )

    # Find the IDs of the component particles inside R200 in the present snapshot
    present_ids = idWithinDisk(present_dd, component, max_r, max_z, :zero)

    # Find the IDs of the component particles inside R200 in the past snapshot
    past_ids = idWithinDisk(past_dd, component, max_r, max_z, :zero)

    # Filter out IDs of component particles that are new to the present snapshot
    # i.e., filter out created particles
    filterExistIDs!(past_dd, component, present_ids)

    # Filter out IDs of component particles that are unique to the past snapshot
    # i.e., filter out destroyed particles
    filterExistIDs!(present_dd, component, past_ids)

    # Find the IDs of the component particles that are inside R200 now, but were outside R200 in the past
    inflow_ids = setdiff(present_ids, past_ids)

    # Find the IDs of the component particles that were inside R200 in the past, but are now outside R200
    outflow_ids = setdiff(past_ids, present_ids)

    # Compute the inflow mass
    inflow_mass = idMass!(present_dd, component, inflow_ids; tracers)

    # Compute the outflow mass
    outflow_mass = idMass!(present_dd, component, outflow_ids; tracers)

    # Compute the net mass
    net_mass_increase = inflow_mass - outflow_mass

    return net_mass_increase, inflow_mass, outflow_mass

end

"""
    computeElementMass(
        data_dict::Dict,
        type::Symbol,
        element::Symbol,
    )::Vector{<:Unitful.Mass}

Compute the total mass of `element` in each cell/particle.

# Arguments

  - `data_dict::Dict`: Data dictionary (see [`makeDataDict`](@ref) for the canonical description).
    This function requires the following blocks to be present, depending on the value of `type`:

      + If `type` == :stellar:
          * `:stellar` => ["MASS", "GME2"]
      + If `type` == :gas:
          * `:gas` => ["MASS", "GMET"]
  - `type::Symbol`: For which cell/particle type the mass will be calculated. The possibilities are `:stellar` and `:gas`.
  - `element::Symbol`: Target element. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).

# Returns

  - The total mass of `element` in each cell/particle.
"""
function computeElementMass(
    data_dict::Dict,
    type::Symbol,
    element::Symbol,
)::Vector{<:Unitful.Mass}

    (
        haskey(ELEMENT_INDEX, element) ||
        throw(ArgumentError("computeElementMass: :$(element) is not a tracked element, \
        the options are the keys of `ELEMENT_INDEX`, see `./src/constants/arepo.jl`"))
    )

    if type == :gas

        block = "GMET"

    elseif type == :stellar

        block = "GME2"

    else

        throw(ArgumentError("computeElementMass: `type` can only be :stellar or :gas, \
        but I got :$(type)"))

    end

    Z = data_dict[type][block]
    M = data_dict[type]["MASS"]

    if any(isempty, [Z, M])

        (
            logging[] &&
            @warn("computeElementMass: I could not compute the masses of :$(element). \
            The metallicities or masses are empty")
        )

        masses = Unitful.Mass[]

    else

        Ze = setPositive(Z[ELEMENT_INDEX[element], :])

        masses = Ze .* M

    end

    return masses

end

@doc raw"""
    computeAbundance(
        data_dict::Dict,
        type::Symbol,
        element::Symbol;
        <keyword arguments>
    )::Vector{Float64}

Compute the abundance of a given element in each cell/particle. The abundance is defined as $n_X / n_H$ where $n_X$ is the number of atoms of element $X$ and $n_H$ the number of hydrogen atoms.

# Arguments

  - `data_dict::Dict`: Data dictionary (see [`makeDataDict`](@ref) for the canonical description).
    This function requires the following blocks to be present, depending on the value of `type`:

      + If `type` == :stellar:
          * `:stellar` => ["MASS", "GME2"]
      + If `type` == :gas:
          * `:gas` => ["MASS", "GMET"]
  - `type::Symbol`: For which cell/particle type the abundance will be calculated. The possibilities are `:stellar` and `:gas`.
  - `element::Symbol`: Target element. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
  - `solar::Bool=false`: If the result will be normalized to the solar abundance.

# Returns

  - The abundance of `element` in each cell/particle.
"""
function computeAbundance(
    data_dict::Dict,
    type::Symbol,
    element::Symbol;
    solar::Bool=false,
)::Vector{Float64}

    element_mass = computeElementMass(data_dict, type, element)
    hydrogen_mass = computeElementMass(data_dict, type, :H)

    if any(isempty, [element_mass, hydrogen_mass])

        logging[] && @warn("computeAbundance: I could not compute the abundance of :$(element)")

        return Float64[]

    end

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
        type::Symbol,
        element::Symbol;
        <keyword arguments>
    )::Float64

Compute the total abundance of a given element, as $n_X / n_H$ where $n_X$ is the number of atoms of element $X$ and $n_H$ the number of hydrogen atoms.

# Arguments

  - `data_dict::Dict`: Data dictionary (see [`makeDataDict`](@ref) for the canonical description).
    This function requires the following blocks to be present, depending on the value of `type`:

      + If `type` == :stellar:
          * `:stellar` => ["MASS", "GME2"]
      + If `type` == :gas:
          * `:gas` => ["MASS", "GMET"]
  - `type::Symbol`: For which cell/particle type the abundance will be calculated. The possibilities are `:stellar` and `:gas`.
  - `element::Symbol`: Target element. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
  - `solar::Bool=false`: If the result will be normalized to the solar abundance or not.

# Returns

  - The total abundance of `element`.
"""
function computeGlobalAbundance(
    data_dict::Dict,
    type::Symbol,
    element::Symbol;
    solar::Bool=false,
)::Float64

    metal_mass = sum(computeElementMass(data_dict, type, element); init=0.0u"Msun")
    hydrogen_mass = sum(computeElementMass(data_dict, type, :H); init=0.0u"Msun")

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
    density3DProjection(
        data_dict::Dict,
        grid::CubicGrid,
        component::Symbol,
        field_type::Symbol;
        <keyword arguments>
    )::Array{Float64,3}

Sample the 3D density field of a given quantity using a cubic grid.

!!! note

    If the source of the field is a set of particles, a simple 3D histogram is used. If instead they are Voronoi cells, the density of the cell that intersects each voxel is used.

# Arguments

  - `data_dict::Dict`: Data dictionary (see [`makeDataDict`](@ref) for the canonical description).
    This function requires the following blocks to be present, depending on the value of `component`:

      + If `component` ∈ [:stellar, :dark_matter, :black_hole, :gas]:
          * `component` => ["MASS", "POS ", "RHO "]
      + If `component` == :Z_stellar
          * `:stellar` => ["MASS", "GZ2 ", "POS ", "RHO "]
      + If `component` ∈ [:hydrogen, :helium]:
          * `:gas` => ["MASS", "POS ", "RHO "]
      + If `component` == :Z_gas:
          * `:gas` => ["MASS", "GZ  ", "POS ", "RHO "]
      + If `component` ∈ [:ionized, :neutral]:
          * `:gas` => ["MASS", "NH  ", "NHP ", "POS ", "RHO "]
      + If `component` ∈ [:br_atomic, :br_molecular]:
          * `:gas` => ["MASS", "NH  ", "NHP ", "PRES", "POS ", "RHO "]
      + If `component` ∈ [:ode_ionized, :ode_atomic, :ode_metals, :ode_dust, :ode_neutral]:
          * `:gas` => ["MASS", "NH  ", "NHP ", "FRAC", "RHO ", "GZ  ", "POS "]
      + If `component` ∈ [:ode_molecular, :ode_stellar]:
          * `:gas` => ["MASS", "FRAC", "RHO ", "POS "]
  - `grid::CubicGrid`: Cubic grid.
  - `component::Symbol`: Target component. It can only be one of the elements of [`COMPONENTS`](@ref).
  - `field_type::Symbol`: If the source of the field are `:particles` or Voronoi `:cells`.
  - `m_unit::Unitful.Units=u"Msun"`: Mass unit.
  - `l_unit::Unitful.Units=u"kpc"`: Length unit.
  - `empty_nan::Bool=true`: If NaN will be put into empty bins, 0 is used otherwise.

# Returns

  - A 3D array with the density at each point of the 3D grid.
"""
function density3DProjection(
    data_dict::Dict,
    grid::CubicGrid,
    component::Symbol,
    field_type::Symbol;
    m_unit::Unitful.Units=u"Msun",
    l_unit::Unitful.Units=u"kpc",
    empty_nan::Bool=true,
)::Array{<:Unitful.Density,3}

    if component ∉ COMPONENTS
        throw(ArgumentError("computeMassDensity: `component` can only be one of the elements \
        of `COMPONENTS` (see `./src/constants/globals.jl`), but I got :$(component)"))
    end

    cp_type = plotParams(Symbol(component, :_mass)).cp_type

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
    positions = data_dict[cp_type]["POS "]

    # Compute the masses of the target component
    masses = computeMass(data_dict, component)

    if any(isempty, [masses, positions])

        (
            logging[] &&
            @warn("density3DProjection: There is missing data, so I will return an empty density \
            field")
        )

        return fill(NaN, (grid.n_bins, grid.n_bins, grid.n_bins))

    end

    # Density unit
    ρ_unit = m_unit * l_unit^-3

    if field_type == :cells

        # Compute the volume of each cell
        cell_volumes = data_dict[cp_type]["MASS"] ./ data_dict[cp_type]["RHO "]

        # Compute the densities of the target component
        densities = ustrip.(ρ_unit, masses ./ cell_volumes)

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

        density = similar(grid.grid, Float64)

        # Compute the density in each voxel
        for i in eachindex(grid.grid)
            density[i] = densities[idxs[i]]
        end

        if empty_nan
            # Set bins with a value of 0 to NaN
            replace!(x -> iszero(x) ? NaN : x, density)
        end

    elseif field_type == :particles

        # Compute the 3D histogram
        mass_histogram = histogram3D(positions, masses, grid; empty_nan)

        density = ustrip.(ρ_unit, mass_histogram ./ grid.bin_volume)

    else

        throw(ArgumentError("density3DProjection: The argument `field_type` must be :cells or \
        :particles, but I got :$(field_type)"))

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
            "\nDensity range - log₁₀(ρ [$(ρ_unit)]) \
            \n  Simulation: $(basename(data_dict[:sim_data].path)) \
            \n  Snapshot:   $(data_dict[:snap_data].global_index) \
            \n  Component:  $(component) \
            \n  Filed type: $(field_type) \
            \n  Min - Max:  $(min_max_ρ) \
            \n  Mean:       $(mean_ρ) \
            \n  Median:     $(median_ρ) \
            \n  Mode:       $(mode_ρ)"
        )

    end

    return density .* ρ_unit

end

"""
    density2DProjection(
        data_dict::Dict,
        grid::CubicGrid,
        component::Symbol,
        field_type::Symbol;
        <keyword arguments>
    )::Matrix{<:SurfaceDensity}

Sample the 3D density field of a given quantity using a cubic grid and then project the field into a given plane.

!!! note

    If the source of the field is a set of particles, a simple 3D histogram is used. If instead they are Voronoi cells, the density of the cell that intersects each voxel is used.

# Arguments

  - `data_dict::Dict`: Data dictionary (see [`makeDataDict`](@ref) for the canonical description).
    This function requires the following blocks to be present, depending on the value of `component`:

      + If `component` ∈ [:stellar, :dark_matter, :black_hole, :gas]:
          * `component` => ["MASS", "POS ", "RHO "]
      + If `component` == :Z_stellar
          * `:stellar` => ["MASS", "GZ2 ", "POS ", "RHO "]
      + If `component` ∈ [:hydrogen, :helium]:
          * `:gas` => ["MASS", "POS ", "RHO "]
      + If `component` == :Z_gas:
          * `:gas` => ["MASS", "GZ  ", "POS ", "RHO "]
      + If `component` ∈ [:ionized, :neutral]:
          * `:gas` => ["MASS", "NH  ", "NHP ", "POS ", "RHO "]
      + If `component` ∈ [:br_atomic, :br_molecular]:
          * `:gas` => ["MASS", "NH  ", "NHP ", "PRES", "POS ", "RHO "]
      + If `component` ∈ [:ode_ionized, :ode_atomic, :ode_metals, :ode_dust, :ode_neutral]:
          * `:gas` => ["MASS", "NH  ", "NHP ", "FRAC", "RHO ", "GZ  ", "POS "]
      + If `component` ∈ [:ode_molecular, :ode_stellar]:
          * `:gas` => ["MASS", "FRAC", "RHO ", "POS "]
  - `grid::CubicGrid`: Cubic grid.
  - `component::Symbol`: Target component. It can only be one of the elements of [`COMPONENTS`](@ref).
  - `field_type::Symbol`: If the source of the field are `:particles` or Voronoi `:cells`.
  - `projection_plane::Symbol=:xy`: Projection plane. The options are `:xy`, `:xz`, and `:yz`.
  - `m_unit::Unitful.Units=u"Msun"`: Mass unit.
  - `l_unit::Unitful.Units=u"kpc"`: Length unit.
  - `empty_nan::Bool=true`: If NaN will be put into empty bins, 0 is used otherwise.

# Returns

  - A 2D array with the surface density at each point of the projected 3D grid.
"""
function density2DProjection(
    data_dict::Dict,
    grid::CubicGrid,
    component::Symbol,
    field_type::Symbol;
    projection_plane::Symbol=:xy,
    m_unit::Unitful.Units=u"Msun",
    l_unit::Unitful.Units=u"kpc",
    empty_nan::Bool=true,
)::Matrix{<:SurfaceDensity}

    density_grid = density3DProjection(
        data_dict,
        grid,
        component,
        field_type;
        m_unit,
        l_unit,
        empty_nan=false,
    )

    # Density unit
    Σ_unit = m_unit * l_unit^-2

    # Compute the mass in each voxel
    mass_grid = density_grid .* grid.bin_volume

    # Project `mass_grid` to the target plane
    if projection_plane == :xy
        dims = 3
    elseif projection_plane == :xz
        dims = 2
    elseif projection_plane == :yz
        dims = 1
    else
        throw(ArgumentError("daDensity2DProjection: The argument `projection_plane` must be \
        :xy, :xz or :yz, but I got :$(projection_plane)"))
    end

    density = dropdims(sum(mass_grid; dims) ./ grid.bin_area; dims)

    if empty_nan
        # Set bins with a value of 0 to NaN
        replace!(x -> iszero(x) ? NaN * Σ_unit : x, density)
    end

    if logging[]

        log_density = filter(x -> !isnan(x) && !isinf(x), log10.(ustrip.(Σ_unit, density)))

        if isempty(log_density)

            min_max_Σ = (NaN, NaN)
            mean_Σ    = NaN
            median_Σ  = NaN
            mode_Σ    = NaN

        else

            min_max_Σ = extrema(log_density)
            mean_Σ    = mean(log_density)
            median_Σ  = median(log_density)
            mode_Σ    = mode(log_density)

        end

        # Print the density range
        @info(
            "\nDensity range - log₁₀(Σ [$(Σ_unit)]) \
            \n  Simulation: $(basename(data_dict[:sim_data].path)) \
            \n  Snapshot:   $(data_dict[:snap_data].global_index) \
            \n  Component:  $(component) \
            \n  Filed type: $(field_type) \
            \n  Plane:      $(projection_plane) \
            \n  Min - Max:  $(min_max_Σ) \
            \n  Mean:       $(mean_Σ) \
            \n  Median:     $(median_Σ) \
            \n  Mode:       $(mode_Σ)"
        )

    end

    return density

end
