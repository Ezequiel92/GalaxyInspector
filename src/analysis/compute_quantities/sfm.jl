####################################################################################################
# Compute derived code quantities for our star formation model
####################################################################################################

"""
    computeDerivedSFMQty(data_dict::Dict, quantity::Symbol)::Vector{<:Number}

Compute a derived code `quantity` for our star formation model.

# Arguments

  - `data_dict::Dict`: Data dictionary (see [`makeDataDict`](@ref) for the canonical description).
    This function requires the following blocks to be present, depending on the value of `component` and `magnitude`:

      + If `magnitude` == :tau_star
          * type => ["RHOC"]
      + If `magnitude` == :tau_rec
          * type => ["RHOC", "FRAC"]
      + If `magnitude` == :tau_cond
          * type => ["RHOC", "FRAC"]
      + If `magnitude` == :tau_dg
          * type => ["RHOC", "FRAC"]
      + If `magnitude` == :tau_dc
          * type => ["RHOC", "FRAC"]
      + If `magnitude` == :tau_ion
          * type => ["RHOC", "PARH", "FRAC"]
      + If `magnitude` == :tau_diss
          * type => ["RHOC", "PARH", "FRAC"]
      + If `magnitude` == :S_d
          * type => ["RHOC", "PARH", "FRAC"]
      + If `magnitude` == :S_H2
          * type => ["RHOC", "PARH", "FRAC"]
      + If `magnitude` == :equivalent_size
          * `:gas` => ["RHOC", "MASS"] or `:stellar` => ["RHOC", "GMAS"]
    where

      + If `component` == :ode_gas_
          * type = `:gas`
      + If `component` == :ode_stellar_
          * type = `:stellar`
  - `quantity::Symbol`: Target quantity. It can only be one of the elements of [`SFM_DERIVED_QTY`](@ref).

# Returns

  - `quantity` for each cell/particle.
"""
function computeDerivedSFMQty(data_dict::Dict, quantity::Symbol)::Vector{<:Number}

    magnitude, component = SFM_DERIVED_QTY_SPLITS[quantity]

    if component == :ode_gas_

        dd = data_dict[:gas]

    elseif component == :ode_stellar_

        dd = data_dict[:stellar]

    else

        throw(ArgumentError("computeDerivedSFMQty: I don't recognize the component :$(component). \
        The only options are :ode_gas_ and :ode_stellar_"))

    end

    # Cell density
    ρ = dd["RHOC"] .* u"mp"

    # NaN time
    t_nan = NaN * u"Myr"

    if magnitude == :tau_star

        ###########################
        # Star formation timescale
        ###########################

        tff = @. sqrt(3π / (32 * Unitful.G * ρ))

        derived_qty = @. tff / εff

    elseif magnitude == :tau_rec

        ##########################
        # Recombination timescale
        ##########################

        fractions = dd["FRAC"]

        isempty(fractions) && return Number[]

        fi = view(fractions, SFM_IDX[:ode_ionized], :)

        derived_qty = @. Unitful.mp / (αH * fi * ρ)

        replace!(x -> isinf(x) ? t_nan : x, derived_qty)

    elseif magnitude == :tau_cond

        #########################
        # Condensation timescale
        #########################

        fractions = dd["FRAC"]

        isempty(fractions) && return Number[]

        fi = view(fractions, SFM_IDX[:ode_ionized], :)
        fa = view(fractions, SFM_IDX[:ode_atomic], :)
        fm = view(fractions, SFM_IDX[:ode_molecular], :)

        fg = @. fi + fa + fm

        derived_qty = @. Unitful.mp / (2 * Rd * ρ * fg)

        replace!(x -> isinf(x) ? t_nan : x, derived_qty)

    elseif magnitude == :tau_dg

        ########################
        # Dust growth timescale
        ########################

        fractions = dd["FRAC"]

        isempty(fractions) && return Number[]

        fa = view(fractions, SFM_IDX[:ode_atomic], :)
        fm = view(fractions, SFM_IDX[:ode_molecular], :)
        fZ = view(fractions, SFM_IDX[:ode_metals], :)
        fd = view(fractions, SFM_IDX[:ode_dust], :)

        fn = @. fa + fm
        Z  = @. fZ + fd

        derived_qty = @. 1.0 / (Cdg * Z * ρ * fn)

        replace!(x -> isinf(x) ? t_nan : x, derived_qty)

    elseif magnitude == :tau_dc

        ##########################
        # Dust creation timescale
        ##########################

        fractions = dd["FRAC"]

        isempty(fractions) && return Number[]

        fa = view(fractions, SFM_IDX[:ode_atomic], :)
        fm = view(fractions, SFM_IDX[:ode_molecular], :)
        fZ = view(fractions, SFM_IDX[:ode_metals], :)

        fn = @. fa + fm

        derived_qty = @. 1.0 / (Cdg * fZ * fn * fn * ρ)

        replace!(x -> isinf(x) ? t_nan : x, derived_qty)

    elseif magnitude == :tau_ion

        ###########################
        # Ionization optical depth
        ###########################

        fractions = dd["FRAC"]

        isempty(fractions) && return Number[]

        fa = view(fractions, SFM_IDX[:ode_atomic], :)
        h  = dd["PARH"]

        derived_qty = @. (fa * ρ / Unitful.mp) * σion * h

    elseif magnitude == :tau_diss

        #############################
        # Dissociation optical depth
        #############################

        fractions = dd["FRAC"]

        isempty(fractions) && return Number[]

        fm = view(fractions, SFM_IDX[:ode_molecular], :)
        h  = dd["PARH"]

        derived_qty = @. (fm * ρ / (2 * Unitful.mp)) * σdiss * h

    elseif magnitude == :S_d

        ########################
        # Dust shielding factor
        ########################

        fractions = dd["FRAC"]

        isempty(fractions) && return Number[]

        fa = view(fractions, SFM_IDX[:ode_atomic], :)
        fm = view(fractions, SFM_IDX[:ode_molecular], :)
        fZ = view(fractions, SFM_IDX[:ode_metals], :)
        fd = view(fractions, SFM_IDX[:ode_dust], :)

        Z = @. fZ + fd
        h = dd["PARH"]

        N  = @. (fa + fm) * ρ * h / Unitful.mp
        Zs = @. Z / SOLAR_METALLICITY

        derived_qty = @. exp(-σd * Zs * N)

    elseif magnitude == :S_H2

        ###########################
        # H2 self-shielding factor
        ###########################

        fractions = dd["FRAC"]

        isempty(fractions) && return Number[]

        fm = view(fractions, SFM_IDX[:ode_molecular], :)
        h  = dd["PARH"]

        x  = @. (fm * ρ * h) / (2 * Unitful.mp * xnf)
        x1 = @. 1.0 + x

        derived_qty = @. ((1.0 - ωH2) / x1^2) + (ωH2 * exp(-0.00085 * sqrt(x1)) / sqrt(x1))

    elseif magnitude == :equivalent_size

        #########################################################
        # Diameter of the sphere with the volume of the gas cell
        #########################################################

        if component == :ode_gas_
            m = dd["MASS"]
        else
            m = dd["GMAS"]
        end

        V = @. m / ρ

        derived_qty = @. 2.0 * cbrt(V / (4π / 3))

    else

        throw(ArgumentError("computeDerivedSFMQty: I don't recognize the magnitude :$(magnitude)"))

    end

    return derived_qty

end
