####################################################################################################
# Compute time related quantities
####################################################################################################

@doc raw"""
    energyIntegrand(a::Real, header::SnapshotHeader)::Float64

The integrand of the integral that converts scale factors into physical times,

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
  - `header::SnapshotHeader`: Header of the relevant snapshot file, containing the cosmological parameters.

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
    computeTime(
        scale_factors::Vector{<:Real},
        header::SnapshotHeader;
        <keyword arguments>
    )::Vector{<:Unitful.Time}

Compute the physical time corresponding to each of the `scale_factors`.

To get the physical time $t$ from the scale factor `a`, one does the integral

```math
t = \frac{1}{H_0} \int_0^a \frac{\mathrm{d}a'}{a' \, \sqrt{\mathcal{E}(a')}} \, ,
```

where

```math
\mathcal{E}(a) = \Omega_\Lambda + \Omega_m \, a^{-3} + \Omega_r \, a^{-4} + \Omega_K \, a^{-2} \, .
```

# Arguments

  - `scale_factors::Vector{<:Real}`: Scale factors.
  - `header::SnapshotHeader`: Header of the relevant snapshot file, containing the cosmological parameters.
  - `a0::Float64=0.0`: Initial scale factor.

# Returns

  - The physical times.
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

To get the physical time $t$ from the scale factor `a`, one does the integral

```math
t = \frac{1}{H_0} \int_0^a \frac{\mathrm{d}a'}{a' \, \sqrt{\mathcal{E}(a')}} \, ,
```

where

```math
\mathcal{E}(a) = \Omega_\Lambda + \Omega_m \, a^{-3} + \Omega_r \, a^{-4} + \Omega_K \, a^{-2} \, .
```

# Arguments

  - `a::Real`: Scale factor.
  - `header::SnapshotHeader`: Header of the relevant snapshot file, containing the cosmological parameters.
  - `a0::Float64=0.0`: Initial scale factor.

# Returns

  - The physical time.
"""
function computeTime(a::Real, header::SnapshotHeader; a0::Float64=0.0)::Unitful.Time

    return computeTime([a], header; a0)[1]

end

"""
    computeTimeStamps(
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
function computeTimeStamps(
    paths::Vector{<:Union{Missing,String}},
)::Tuple{Vector{Float64},Vector{Float64},Vector{<:Unitful.Time},Vector{<:Unitful.Time}}

    snapshot_paths = filter(!ismissing, paths)

    if isempty(snapshot_paths)
        (
            logging[] &&
            @warn("computeTimeStamps: `paths` is empty or full of missing, so I will return NaNs")
        )
        return [NaN], [NaN], [NaN * u"s"], [NaN * u"s"]
    end

    first_snapshot = first(snapshot_paths)

    if isSnapCosmological(first_snapshot)

        # For cosmological simulations, the time field in the Header of the snapshot is the scale factor
        scale_factors  = [readTime(path) for path in snapshot_paths]
        redshifts      = @. (1.0 / scale_factors) - 1.0
        physical_times = computeTime(scale_factors, readSnapHeader(first_snapshot))
        lookback_times = maximum(physical_times) .- physical_times

    else

        # Compute the factor for internal units of time
        u_time = internalUnits("CLKT", first_snapshot)

        # a = 1.0 for non-cosmological simulations
        scale_factors  = ones(length(snapshot_paths))
        # z = 0.0 for non-cosmological simulations
        redshifts      = zeros(length(snapshot_paths))
        # For non-cosmological simulations, the time in the snapshot is the physical time
        physical_times = [readTime(path) * u_time for path in snapshot_paths]
        lookback_times = maximum(physical_times) .- physical_times

    end

    return scale_factors, redshifts, physical_times, lookback_times

end

"""
    computeStellarBirthTime(data_dict::Dict)::Vector{<:Unitful.Time}

Compute the stellar birth times.

# Arguments

  - `data_dict::Dict`: Data dictionary (see [`makeDataDict`](@ref) for the canonical description).
    This function requires the following blocks to be present:

      + `:stellar` => ["GAGE"]

# Returns

  - The stellar birth times.
"""
function computeStellarBirthTime(data_dict::Dict)::Vector{<:Unitful.Time}

    birth_ticks = data_dict[:stellar]["GAGE"]

    if isempty(birth_ticks)
        (
            logging[] &&
            @warn("computeStellarBirthTime: There is no data for the stellar birth times, \
            so I will return an empty array")
        )
        return Unitful.Time[]
    end

    if data_dict[:sim_data].cosmological
        birth_times = computeTime(birth_ticks, data_dict[:snap_data].header)
    else
        birth_times = birth_ticks
    end

    return birth_times

end

"""
    computeStellarAge(data_dict::Dict)::Vector{<:Unitful.Time}

Compute the age of the stars.

# Arguments

  - `data_dict::Dict`: Data dictionary (see [`makeDataDict`](@ref) for the canonical description).
    This function requires the following blocks to be present:

      + `:stellar` => ["GAGE"]

# Returns

  - The stellar ages.
"""
function computeStellarAge(data_dict::Dict)::Vector{<:Unitful.Time}

    birth_times = computeStellarBirthTime(data_dict)

    if isempty(birth_times)
        (
            logging[] &&
            @warn("computeStellarAge: There is no data for the stellar ages, so I will return an \
            empty array")
        )
        return Unitful.Time[]
    end

    return setPositive(data_dict[:snap_data].physical_time .- birth_times)

end

"""
    computeSFR(
        data_dict::Dict;
        <keyword arguments>
    )::Vector{<:Unitful.MassFlow}

Compute the star formation rate of each stellar particle.

For stellar particles younger than `age_limit`, the SFR is its mass divided by `age_limit`. For older particles it is 0.

# Arguments

  - `data_dict::Dict`: Data dictionary (see [`makeDataDict`](@ref) for the canonical description).
    This function requires the following blocks to be present:

      + `:stellar` => ["MASS", "GAGE"]
  - `age_limit::Unitful.Time=AGE_RESOLUTION`: Age limit for the SFR.

# Returns

  - The star formation rate of each stellar particle.
"""
function computeSFR(
    data_dict::Dict;
    age_limit::Unitful.Time=AGE_RESOLUTION,
)::Vector{<:Unitful.MassFlow}

    # Compute the stellar ages
    ages = computeStellarAge(data_dict)

    if isempty(ages)
        (
            logging[] &&
            @warn("computeSFR: There is no data for the stellar ages, so I will return an \
            empty array")
        )
        return Unitful.MassFlow[]
    end

    sfr = zeros(typeof(1.0u"Msun * yr^-1"), length(ages))

    # Find the stellar particles younger than `age_limit`
    idxs = map(x -> x <= age_limit, ages)

    # Compute the SFR
    sfr[idxs] .= data_dict[:stellar]["MASS"][idxs] ./ age_limit

    return sfr

end

"""
    computeSSFR(
        data_dict::Dict;
        <keyword arguments>
    )::Vector{<:Unitful.MassFlow}

Compute the specific star formation rate of each stellar particle.

For stellar particles younger than `age_limit`, the sSFR is 1 / `age_limit`. For older particles it is 0.

# Arguments

  - `data_dict::Dict`: Data dictionary (see [`makeDataDict`](@ref) for the canonical description).
    This function requires the following blocks to be present:

      + `:stellar` => ["MASS", "GAGE"]
  - `age_limit::Unitful.Time=AGE_RESOLUTION`: Age limit for the SFR.

# Returns

  - The specific star formation rate of each stellar particle.
"""
function computeSSFR(
    data_dict::Dict;
    age_limit::Unitful.Time=AGE_RESOLUTION,
)::Vector{<:Unitful.MassFlow}

    # Compute the stellar ages
    ages = computeStellarAge(data_dict)

    if isempty(ages)
        (
            logging[] &&
            @warn("computeSSFR: There is no data for the stellar ages, so I will return an \
            empty array")
        )
        return Unitful.MassFlow[]
    end

    ssfr = zeros(typeof(1.0u"yr^-1"), length(ages))

    # Find the stellar particles younger than `age_limit`
    idxs = map(x -> x <= age_limit, ages)

    inv_age_limit = 1 / age_limit

    # Compute the SFR
    ssfr[idxs] .= inv_age_limit

    return ssfr

end

@doc raw"""
    computeDepletionTime(
        masses::Vector{<:Unitful.Mass},
        sfrs::Vector{<:Unitful.MassFlow},
    )::Vector{<:Unitful.Time}

Compute the depletion time,

```math
t_\mathrm{ff} = \frac{M_\mathrm{gas}}{\dot{M}_\star} \, .
```

# Arguments

  - `masses::Vector{<:Unitful.Mass}`: Mass of every cell.
  - `sfrs::Vector{<:Unitful.MassFlow}`: The SFR associated to each cell.

# Returns

  - The depletion times.
"""
function computeDepletionTime(
    masses::Vector{<:Unitful.Mass},
    sfrs::Vector{<:Unitful.MassFlow},
)::Vector{<:Unitful.Time}

    if any(isempty, [masses, sfrs])
        (
            logging[] &&
            @warn("computeDepletionTime: There is missing data, so I will return an empty array")
        )
        return Unitful.Time[]
    end

    return @. masses / sfrs

end

@doc raw"""
    computeDepletionTime(data_dict::Dict, component::Symbol)::Vector{<:Unitful.Time}

Compute the depletion time,

```math
t_\mathrm{ff} = \frac{M_\mathrm{gas}}{\dot{M}_\star} \, .
```

# Arguments

  - `data_dict::Dict`: Data dictionary (see [`makeDataDict`](@ref) for the canonical description).
    This function requires the following blocks to be present, depending on the value of `component`:

      + If `component` ∈ [:gas, :hydrogen, :helium]:
          * `:gas` => ["SFR ", "MASS"]
      + If `component` == :Z_gas:
          * `:gas` => ["SFR ", "MASS", "GZ  "]
      + If `component` ∈ [:ionized, :neutral]:
          * `:gas` => ["SFR ", "MASS", "NH  ", "NHP "]
      + If `component` ∈ [:br_atomic, :br_molecular]:
          * `:gas` => ["SFR ", "MASS", "NH  ", "NHP ", "PRES"]
      + If `component` ∈ [:ode_ionized, :ode_atomic, :ode_metals, :ode_dust, :ode_neutral, :ode_cold]:
          * `:gas` => ["SFR ", "MASS", "NH  ", "NHP ", "FRAC", "RHO ", "GZ  "]
      + If `component` ∈ [:ode_molecular, :ode_stellar, :ode_molecular_stellar]:
          * `:gas` => ["SFR ", "MASS", "FRAC", "RHO "]
  - `component::Symbol`: Target component. It can only be one of the elements of [`COMPONENTS`](@ref).

# Returns

  - The depletion time of `component`.
"""
function computeDepletionTime(data_dict::Dict, component::Symbol)::Vector{<:Unitful.Time}

    if component ∉ COMPONENTS || component ∈ [:stellar, :dark_matter, :black_hole, :Z_stellar]
        throw(ArgumentError("computeDepletionTime: `component` can only be one of the gas elements \
        of `COMPONENTS` (see `./src/constants/globals.jl`), but I got :$(component)"))
    end

    masses = computeMass(data_dict, component)
    sfrs   = data_dict[:gas]["SFR "]

    return computeDepletionTime(masses, sfrs)

end
