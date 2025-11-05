####################################################################################################
# Compute time related quantities
####################################################################################################

#################
# Base functions
#################

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
            logging[] &&
            @warn("computeKineticEnergy: There is missing data, so I will return an \
            empty array")
        )
        return Unitful.Energy[]
    end

    return [0.5 * mass * sum(vel .* vel) for (mass, vel) in zip(masses, eachcol(velocities))]

end

"""
    computePotentialEnergy(
        masses::Vector{<:Unitful.Mass},
        potential::Vector{<:SpecificEnergy},
    )::Vector{<:Unitful.Energy}

Compute the gravitational potencial energy.

# Arguments

  - `masses::Vector{<:Unitful.Mass}`: Mass of every cell/particle.
  - `potential::Vector{<:SpecificEnergy}`: Specific potential energy of every cell/particle.

# Returns

  - The gravitational potencial energy of each cell/particle.
"""
function computePotentialEnergy(
    masses::Vector{<:Unitful.Mass},
    potential::Vector{<:SpecificEnergy},
)::Vector{<:Unitful.Energy}

    if any(isempty, [potential, masses])
        (
            logging[] &&
            @warn("computePotentialEnergy: There is missing data, so I will return an empty array")
        )
        return Unitful.Energy[]
    end

    return potential .* masses

end

"""
    computeTemperature(
        internal_energy::Vector{<:SpecificEnergy},
        electron_fraction::Vector{Float32},
    )::Vector{<:Unitful.Temperature}

Compute the gas temperature.

# Arguments

  - `internal_energy::Vector{<:SpecificEnergy}`: Specific internal energy of every gas cell.
  - `electron_fraction::Vector{Float32}`: Number fraction of electrons in every gas cell.

# Returns

  - The temperature of each gas cell.
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

###################
# Derive functions
###################

"""
    computeKineticEnergy(data_dict::Dict, component::Symbol)::Vector{<:Unitful.Energy}

Compute the kinetic energy.

# Arguments

  - `data_dict::Dict`: Data dictionary (see [`makeDataDict`](@ref) for the canonical description).
    This function requires the following blocks to be present, depending on the value of `component`:

      + If `component` ∈ [:stellar, :dark_matter, :black_hole, :gas]:
          * `component` => ["VEL ", "MASS"]
      + If `component` == :Z_stellar
          * `:stellar` => ["VEL ", "MASS", "GZ2 "]
      + If `component` ∈ [:hydrogen, :helium]:
          * `:gas` => ["VEL ", "MASS"]
      + If `component` == :Z_gas:
          * `:gas` => ["VEL ", "MASS", "GZ  "]
      + If `component` ∈ [:ionized, :neutral]:
          * `:gas` => ["VEL ", "MASS", "NH  ", "NHP "]
      + If `component` ∈ [:br_atomic, :br_molecular]:
          * `:gas` => ["VEL ", "MASS", "NH  ", "NHP ", "PRES"]
      + If `component` ∈ [:ode_ionized, :ode_atomic, :ode_metals, :ode_dust, :ode_neutral]:
          * `:gas` => ["VEL ", "MASS", "NH  ", "NHP ", "FRAC", "RHO ", "GZ  "]
      + If `component` ∈ [:ode_molecular, :ode_stellar]:
          * `:gas` => ["VEL ", "MASS", "FRAC", "RHO "]
  - `component::Symbol`: Target component. It can only be one of the elements of [`COMPONENTS`](@ref).

# Returns

  - The kinetic energy of each cell/particle.
"""
function computeKineticEnergy(data_dict::Dict, component::Symbol)::Vector{<:Unitful.Energy}

    if component ∉ COMPONENTS
        throw(ArgumentError("computeKineticEnergy: `component` can only be one of the elements \
        of `COMPONENTS` (see `./src/constants/globals.jl`), but I got :$(component)"))
    end

    if component ∈ [:stellar, :dark_matter, :gas, :black_hole]
        type = component
    elseif component == :Z_stellar
        type = :stellar
    else
        type = :gas
    end

    masses     = computeMass(data_dict, component)
    velocities = data_dict[type]["VEL "]

    return computeKineticEnergy(masses, velocities)

end

"""
    computePotentialEnergy(data_dict::Dict, component::Symbol)::Vector{<:Unitful.Energy}

Compute the gravitational potencial energy.

# Arguments

  - `data_dict::Dict`: Data dictionary (see [`makeDataDict`](@ref) for the canonical description).
    This function requires the following blocks to be present, depending on the value of `component`:

      + If `component` ∈ [:stellar, :dark_matter, :black_hole, :gas]:
          * `component` => ["POT ", "MASS"]
      + If `component` == :Z_stellar
          * `:stellar` => ["POT ", "MASS", "GZ2 "]
      + If `component` ∈ [:hydrogen, :helium]:
          * `:gas` => ["POT ", "MASS"]
      + If `component` == :Z_gas:
          * `:gas` => ["POT ", "MASS", "GZ  "]
      + If `component` ∈ [:ionized, :neutral]:
          * `:gas` => ["POT ", "MASS", "NH  ", "NHP "]
      + If `component` ∈ [:br_atomic, :br_molecular]:
          * `:gas` => ["POT ", "MASS", "NH  ", "NHP ", "PRES"]
      + If `component` ∈ [:ode_ionized, :ode_atomic, :ode_metals, :ode_dust, :ode_neutral]:
          * `:gas` => ["POT ", "MASS", "NH  ", "NHP ", "FRAC", "RHO ", "GZ  "]
      + If `component` ∈ [:ode_molecular, :ode_stellar]:
          * `:gas` => ["POT ", "MASS", "FRAC", "RHO "]
  - `component::Symbol`: Target component. It can only be one of the elements of [`COMPONENTS`](@ref).

# Returns

  - The gravitational potencial energy of each cell/particle.
"""
function computePotentialEnergy(data_dict::Dict, component::Symbol)::Vector{<:Unitful.Energy}

    if component ∉ COMPONENTS
        throw(ArgumentError("computePotentialEnergy: `component` can only be one of the elements \
        of `COMPONENTS` (see `./src/constants/globals.jl`), but I got :$(component)"))
    end

    if component ∈ [:stellar, :dark_matter, :gas, :black_hole]
        type = component
    elseif component == :Z_stellar
        type = :stellar
    else
        type = :gas
    end

    masses    = computeMass(data_dict, component)
    potential = data_dict[type]["POT "]

    return computePotentialEnergy(masses, potential)


end

"""
    computeTotalEnergy(data_dict::Dict, component::Symbol)::Vector{<:Unitful.Energy}

Compute the total energy (kinetic + potential).

# Arguments

  - `data_dict::Dict`: Data dictionary (see [`makeDataDict`](@ref) for the canonical description).
    This function requires the following blocks to be present, depending on the value of `component`:

      + If `component` ∈ [:stellar, :dark_matter, :black_hole, :gas]:
          * `component` => ["VEL ", "POT ", "MASS"]
      + If `component` == :Z_stellar
          * `:stellar` => ["VEL ", "POT ", "MASS", "GZ2 "]
      + If `component` ∈ [:hydrogen, :helium]:
          * `:gas` => ["VEL ", "POT ", "MASS"]
      + If `component` == :Z_gas:
          * `:gas` => ["VEL ", "POT ", "MASS", "GZ  "]
      + If `component` ∈ [:ionized, :neutral]:
          * `:gas` => ["VEL ", "POT ", "MASS", "NH  ", "NHP "]
      + If `component` ∈ [:br_atomic, :br_molecular]:
          * `:gas` => ["VEL ", "POT ", "MASS", "NH  ", "NHP ", "PRES"]
      + If `component` ∈ [:ode_ionized, :ode_atomic, :ode_metals, :ode_dust, :ode_neutral]:
          * `:gas` => ["VEL ", "POT ", "MASS", "NH  ", "NHP ", "FRAC", "RHO ", "GZ  "]
      + If `component` ∈ [:ode_molecular, :ode_stellar]:
          * `:gas` => ["VEL ", "POT ", "MASS", "FRAC", "RHO "]
  - `component::Symbol`: Target component. It can only be one of the elements of [`COMPONENTS`](@ref).

# Returns

  - The total energy of each cell/particle.
"""
function computeTotalEnergy(data_dict::Dict, component::Symbol)::Vector{<:Unitful.Energy}

    Ep = computePotentialEnergy(data_dict, component)
    Ek = computeKineticEnergy(data_dict, component)

    if any(isempty, [Ek, Ep])

        logging[] && @warn("computeTotalEnergy: There is missing, so I will return an empty array")

        return Unitful.Energy[]

    end

    return Ek .+ Ep

end
