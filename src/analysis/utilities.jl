####################################################################################################
# General utilities to interact with the simulation data.
####################################################################################################

@doc raw"""
Time factor for the SF model, without the fraction factors.

τ_star(ρ_cell)    $\equiv \tau_\mathrm{star}$
τ_rec(ρ_cell)     $\equiv \tau_\mathrm{rec} \, f_i$
τ_cond(ρ_cell, Z) $\equiv \tau_\mathrm{cond} \, (1 - f_s)$

"""
τ_star(ρ_cell) = C_star / sqrt(ρ_cell)
τ_rec(ρ_cell) = C_rec / ρ_cell
τ_cond(ρ_cell, Z) = C_cond / (ρ_cell * (Z + Zeff))

"""
    getUnitLabel(factor::Int, unit::Unitful.Units; <keyword arguments>)::AbstractString

Construct the unit part of an axis label.

# Arguments

  - `factor::Int`: Exponential factor to scale down the units. If different from 0, a term of the form 10^`factor` will be added to the label.
  - `unit::Unitful.Units`: Unit of the axis.
  - `latex::Bool=true`: If the output will be a `LaTeXString`, or a plain `String`.

# Returns

  - The `LaTeXString` or `String`: "10^`factor` `unit`". The `factor` term only appears if `factor` != 0, the unit term only appears if `unit` != `Unitful.NoUnits`.
"""
function getUnitLabel(factor::Int, unit::Unitful.Units; latex::Bool=true)::AbstractString

    if latex

        # Replace special characters for their LaTeX counterpart
        str_unit = replace(
            string(unit),
            "M⊙" => raw"M_{\odot}",
            r"\^(-?[0-9]+)" => s"^{\1}",
            " " => raw"\,",
        )

        if iszero(factor)
            if isempty(str_unit)
                out_str = ""
            else
                out_str = L"\mathrm{%$str_unit}"
            end
        else
            if isempty(str_unit)
                out_str = L"10^{%$factor}"
            else
                out_str = L"10^{%$factor} \, \mathrm{%$str_unit}"
            end
        end

    else

        str_unit = string(unit)

        if iszero(factor)
            if isempty(str_unit)
                out_str = ""
            else
                out_str = str_unit
            end
        else
            if isempty(str_unit)
                out_str = "10^$factor"
            else
                out_str = "10^$factor $str_unit"
            end
        end

    end

    return out_str

end

"""
    getLabel(
        label::AbstractString,
        factor::Int,
        unit::Unitful.Units;
        <keyword arguments>
    )::AbstractString

Construct an axis label.

# Arguments

  - `label::AbstractString`: Variable name.
  - `factor::Int`: Exponential factor to scale down the units. If different from 0, a term of the form 10^`factor` will be added to the label.
  - `unit::Unitful.Units`: Unit of the axis.
  - `latex::Bool=true`: If the output will be a `LaTeXString`, or a plain `String`.

# Returns

  - The `LaTeXString` or `String`: "`label` [10^`factor` `unit`]". If `label` is "", an empty string is returned. The `factor` term only appears if `factor` != 0, the unit term only appears if `unit` != `Unitful.NoUnits`, and the brackets only appears if there are a factor and/or a unit term.
"""
function getLabel(
    label::AbstractString,
    factor::Int,
    unit::Unitful.Units;
    latex::Bool=true,
)::AbstractString

    !isempty(label) || return ""

    unit_label = getUnitLabel(factor, unit; latex)

    if isempty(unit_label)
        return latex ? L"%$label" : label
    end

    return latex ? L"%$label [%$unit_label]" : "$label [$unit_label]"

end

"""
    barPlotLabelFormater(x::Number)::LaTeXString

Format a number to be a barplot label.

For values between 0 and 0.01 the label will be "< 0.01", otherwise it will be the value itself with 2 digits.

# Arguments

  - `x::Number`: Value to be formated.

# Returns

  - The bar label.
"""
function barPlotLabelFormater(x::Number)::LaTeXString

    if 0 < x < 0.01
        return L"< \, 0.01"
    end

    return latexstring(round(x; digits=2))

end

"""
    barPlotLabelFormater(x::LaTeXString)::LaTeXString

Format a number to be a barplot label.

Method for compatibility with the barplot! function of [Makie](https://docs.makie.org/stable/).

# Arguments

  - `x::Number`: Value to be formated.

# Returns

  - The bar label.
"""
barPlotLabelFormater(x::LaTeXString)::LaTeXString = x

"""
    formatError(q_mean::Number, q_error::Number)::NTuple{2,<:Number}

Nicely format a magnitude with uncertainty.

It follows the traditional rules for error presentation: the error has only one significant digit, unless such digit is a one, in which case two significant digits are used. The mean will have as many digits as to match the last significant position of the error. An error equal to 0 will leave the mean unchanged.

# Arguments

  - `q_mean::Number`: Mean value.
  - `q_error::Number`: Error value. It must be positive.

# Returns

  - A tuple with the formatted mean and error values.

# Examples

```julia-repl
julia> formatError(69.42069, 0.038796)
(69.42, 0.04)

julia> formatError(69.42069, 0.018796)
(69.421, 0.019)

julia> formatError(15.42, 0.00004)
(15.42, 4.0e-5)

julia> formatError(69.42069, 0.0)
(69.42069, 0.0)

julia> formatError(69.42069, 93.4)
(70.0, 90.0)

julia> formatError(69.42069, 123.4)
(70.0, 120.0)

julia> formatError(15.42069, 16.4)
(15.0, 16.0)
```
"""
function formatError(q_mean::Number, q_error::Number)::NTuple{2,<:Number}

    # Positive error check
    (
        q_error >= zero(q_error) ||
        throw(ArgumentError("formatError: `q_error` must be positive, but I got \
        q_error = $(q_error)"))
    )

    # Use the values without units
    mean = ustrip(q_mean)
    error = ustrip(q_error)

    if iszero(error)

        round_mean = mean
        round_error = error

    else

        sigdigit_pos = abs(log10(error))

        if error < 1.0
            first_digit = trunc(error * 10.0^(floor(sigdigit_pos) + 1.0))
            extra = first_digit == 1.0 ? 1 : 0
            digits = ceil(Int, sigdigit_pos) + extra
            round_mean = round(mean; digits)
            round_error = round(error, sigdigits=1 + extra)
        else
            first_digit = trunc(error * 10.0^(-floor(sigdigit_pos)))
            extra = first_digit == 1.0 ? 2 : 1
            sigdigits = ceil(Int, log10(abs(mean))) - ceil(Int, sigdigit_pos) + extra
            round_mean = round(mean; sigdigits)
            round_error = round(error, sigdigits=extra)
        end

    end

    return round_mean * unit(q_mean), round_error * unit(q_error)

end

"""
    flattenGrid(cubic_grid::CubicGrid)::SquareGrid

Using a `CubicGrid` construct a `SquareGrid` with the same center, number of bins, and physical side length.

# Arguments

  - `cubic_grid::CubicGrid`: Cubic grid.

# Returns

  - A square grid.
"""
function flattenGrid(cubic_grid::CubicGrid)::SquareGrid

    physical_size = cubic_grid.physical_size
    n_bins = cubic_grid.n_bins

    bin_width  = physical_size / n_bins
    shift = 0.5 * (physical_size - bin_width)

    center = [cubic_grid.x_ticks[1], cubic_grid.y_ticks[1], cubic_grid.z_ticks[1]] .+ shift

    return SquareGrid(physical_size, n_bins; center)

end

"""
    isCosmological(path::String)::Bool

Check if the snapshot in `path` comes from a cosmological simulation.

!!! note

    If each snapshot is made of multiple files, I'll read the first chunck to check if the simulation is cosmological.

# Arguments

  - `path::String`: Path to the snapshot file or folder.

# Returns

  - If the simulation is cosmological

      + `false` -> Newtonian simulation    (`ComovingIntegrationOn` = 0, `Redshift` = 0.0).
      + `true`  -> Cosmological simulation (`ComovingIntegrationOn` = 1, `Redshift` != 0.0).
"""
function isCosmological(path::String)::Bool

    if isfile(path)

        (
            HDF5.ishdf5(path) ||
            throw(ArgumentError("isCosmological: The file $(path) is not in the HDF5 format, \
            I don't know how to read it"))
        )

        file_path = path

    elseif isdir(path)

        sub_files = glob("$(SNAP_BASENAME)_*.*.hdf5", path)

        (
            !isempty(sub_files) && all(HDF5.ishdf5, sub_files) ||
            throw(ArgumentError("isCosmological: The directory $(path) does not contain \
            snapshot sub-files in the HDF5 format"))
        )

        file_path = minimum(sub_files)

    else

        throw(ArgumentError("isCosmological: $(path) does not exist as a file or folder"))

    end

    cosmological = h5open(file_path, "r") do snapshot
        if "Parameters" ∈ keys(snapshot)
            # If the param.txt is saved in the snapshot metadata, read `ComovingIntegrationOn`
            read_attribute(snapshot["Parameters"], "ComovingIntegrationOn")
        else
            # Otherwise, use the readshift in the header
            !iszero(read_attribute(snapshot["Header"], "Redshift"))
        end
    end

    return cosmological

end

"""
    internalUnits(quantity::String, path::String)::Union{Unitful.Quantity,Unitful.Units}

Get the factor to convert a plain number into a [Unitful](https://github.com/PainterQubits/Unitful.jl) quantity, using the correct internal code units.

# Arguments

  - `quantity::String`: Target quantity. The options are the keys of [`QUANTITIES`](@ref).
  - `path::String`: Path to the snapshot file or folder.

# Returns

  - A [Unitful](https://github.com/PainterQubits/Unitful.jl) quantity or unit.
"""
function internalUnits(quantity::String, path::String)::Union{Unitful.Quantity,Unitful.Units}

    (
        quantity ∈ keys(QUANTITIES) ||
        throw(ArgumentError("internalUnits: `quantity` should be one of the keys of \
        `QUANTITIES` but I got $(quantity), see the options in `./src/constants/globals.jl`"))
    )

    header = readSnapHeader(path)
    cosmological = isCosmological(path)

    a0 = cosmological ? header.time : 1.0
    h0 = cosmological ? header.h0 : 1.0

    # Set up the struct for unit conversion
    IU = InternalUnits(; l_unit=header.l_unit, m_unit=header.m_unit, v_unit=header.v_unit, a0, h0)

    dimensions = QUANTITIES[quantity].dimensions
    unit = QUANTITIES[quantity].unit

    if unit == :internal
        if dimensions == Unitful.𝐌

            # From internal units to M⊙
            return IU.m_cosmo

        elseif dimensions == Unitful.𝐋

            if !PHYSICAL_UNITS && !cosmological
                @warn(
                    "internalUnits: You have set the unit system to use comoving lengths \
                    (PHYSICAL_UNITS = $(PHYSICAL_UNITS)), but the simulation is not \
                    cosmological. I'll keep the lengths physical. Check `PHYSICAL_UNITS` \
                    in `constants/globals.jl`",
                    maxlog=1,
                )
            end

            # From internal units to kpc
            if !PHYSICAL_UNITS && cosmological
                return IU.x_comoving
            else
                return IU.x_cosmo
            end

        elseif dimensions == Unitful.𝐓

            # From internal units to Myr, for non-cosmological simulations,
            # and to a dimensionless quantity for cosmological simulations
            return cosmological ? Unitful.NoUnits : IU.t_cosmo

        elseif dimensions == Unitful.𝐌 * Unitful.𝐋^-3

            # From internal units to g * cm^-3
            return IU.rho_cgs

        elseif dimensions == Unitful.𝐋^2 * Unitful.𝐓^-2

            # From internal units to erg * g^-1
            return IU.U_cgs

        elseif dimensions == Unitful.𝐋 * Unitful.𝐓^-1

            # From internal units to km * s^-1
            return IU.v_cosmo

        elseif dimensions == Unitful.𝐌 * Unitful.𝐋^-1 * Unitful.𝐓^-2

            # From internal units to Pa
            return IU.P_Pa

        else

            error("internalUnits: I don't know the internal units of a quantity \
            with dimensions $(dimensions)")

        end

    elseif unit == :gvel

        # Special case for "G_Vel" (velocity of the group)
        # See the TNG documentation https://www.tng-project.org/data/docs/specifications/
        return IU.v_cosmo / a0^1.5

    else

        return unit

    end

end

"""
    snapshotTypes(data_dict::Dict)::Vector{Symbol}

Find which cell/particle types are part of the keys of `data_dict`.

# Arguments

  - `data_dict::Dict`: A dictionary.

# Returns

  - A vector with the cell/particle types.
"""
snapshotTypes(data_dict::Dict)::Vector{Symbol} = collect(keys(PARTICLE_INDEX) ∩ keys(data_dict))

"""
    snapshotTypes(path::String)::Vector{Symbol}

Find which cell/particle types are part of the snapshot in `path`.

!!! note

    If each snapshot is made of multiple files, I'll check the first chunck.

# Arguments

  - `path::String`: Path to the snapshot file or folder.

# Returns

  - A vector with the cell/particle types.
"""
function snapshotTypes(path::String)::Vector{Symbol}

    if isfile(path)

        (
            HDF5.ishdf5(path) ||
            throw(ArgumentError("snapshotTypes: The file $(path) is not in the HDF5 format, \
            I don't know how to read it"))
        )

        file_path = path

    elseif isdir(path)

        sub_files = glob("$(SNAP_BASENAME)_*.*.hdf5", path)

        (
            !isempty(sub_files) && all(HDF5.ishdf5, sub_files) ||
            throw(ArgumentError("snapshotTypes: The directory $(path) does not contain \
            snapshot sub-files in the HDF5 format"))
        )

        file_path = minimum(sub_files)

    else

        throw(ArgumentError("snapshotTypes: $(path) does not exist as a file or folder"))

    end

    snapshot_types = h5open(file_path, "r") do snapshot
        collect(keys(PARTICLE_TYPE) ∩ keys(snapshot))
    end

    return snapshot_types

end

"""
    groupcatTypes(data_dict::Dict)::Vector{Symbol}

Find which group catalog data types are part of the keys of `data_dict`.

# Arguments

  - `data_dict::Dict`: A dictionary.

# Returns

  - A vector with the group catalog data types.
"""
groupcatTypes(data_dict::Dict)::Vector{Symbol} = [:group, :subhalo] ∩ keys(data_dict)

"""
    groupcatTypes(path::String)::Vector{Symbol}

Find which group catalog data types are part of the snapshot in `path`.

!!! note

    If each snapshot is made of multiple files, I'll check the first chunck.

# Arguments

  - `path::String`: Path to the snapshot file or folder.

# Returns

  - A vector with the group catalog data types.
"""
function groupcatTypes(path::String)::Vector{Symbol}

    if isfile(path)

        (
            HDF5.ishdf5(path) ||
            throw(ArgumentError("groupcatTypes: The file $(path) is not in the HDF5 format, \
            I don't know how to read it"))
        )

        file_path = path

    elseif isdir(path)

        sub_files = glob("$(SNAP_BASENAME)_*.*.hdf5", path)

        (
            !isempty(sub_files) && all(HDF5.ishdf5, sub_files) ||
            throw(ArgumentError("groupcatTypes: The directory $(path) does not contain \
            snapshot sub-files in the HDF5 format"))
        )

        file_path = minimum(sub_files)

    else

        throw(ArgumentError("groupcatTypes: $(path) does not exist as a file or folder"))

    end

    groupcat_types = h5open(file_path, "r") do snapshot
        [:group, :subhalo] ∩ keys(snapshot)
    end

    return groupcat_types

end

"""
    plotParams(quantity::Symbol)::PlotParams

Select the plotting parameters for a given `quantity`.

# Arguments

  - `quantity::Symbol`: The options are:

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
      + `:generic_mass`                -> Parameters for plots with several diferent masses.
      + `:stellar_number`              -> Number of stellar particles.
      + `:gas_number`                  -> Number of gas cells.
      + `:dm_number`                   -> Number of dark matter particles.
      + `:bh_number`                   -> Number of black hole particles.
      + `:molecular_fraction`          -> Gas mass fraction of molecular hydrogen.
      + `:br_molecular_fraction`       -> Gas mass fraction of molecular hydrogen, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_fraction`             -> Gas mass fraction of atomic hydrogen.
      + `:ionized_fraction`            -> Gas mass fraction of ionized hydrogen.
      + `:neutral_fraction`            -> Gas mass fraction of neutral hydrogen.
      + `:molecular_neutral_fraction`  -> Fraction of molecular hydrogen in the neutral gas.
      + `:mol_eq_quotient`             -> Equilibrium quotient for the molecular fraction equation of the SF model.
      + `:ion_eq_quotient`             -> Equilibrium quotient for the ionized fraction equation of the SF model.
      + `:generic_fraction`            -> Parameters for plots with several diferent fraction.
      + `:gas_mass_density`            -> Gas mass density.
      + `:hydrogen_mass_density`       -> Hydrogen mass density.
      + `:gas_number_density`          -> Gas number density.
      + `:molecular_number_density`    -> Molecular hydrogen number density.
      + `:br_molecular_number_density` -> Molecular hydrogen number density, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_number_density`       -> Atomic hydrogen number density.
      + `:ionized_number_density`      -> Ionized hydrogen number density.
      + `:neutral_number_density`      -> Neutral hydrogen number density.
      + `:stellar_area_density`        -> Stellar area mass density, for a radius of `DISK_R`.
      + `:gas_area_density`            -> Gas mass surface density, for a radius of `DISK_R`.
      + `:molecular_area_density`      -> Molecular mass surface density, for a radius of `DISK_R`.
      + `:br_molecular_area_density`   -> Molecular mass surface density, for a radius of `DISK_R`, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_area_density`         -> Atomic hydrogen area mass density, for a radius of `DISK_R`.
      + `:ionized_area_density`        -> Ionized hydrogen area mass density, for a radius of `DISK_R`.
      + `:neutral_area_density`        -> Neutral mass surface density, for a radius of `DISK_R`.
      + `:sfr_area_density`            -> Star formation rate area density, for the last `AGE_RESOLUTION` and a radius of `DISK_R`.
      + `:generic_area_density`        -> Parameters for plots with several diferent area densities.
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
      + `:mass_accretion`              -> Gas accretion rate. Positive values mean gas infall into the virial radius ``R_{200}``, and negative values mean outflow.
      + `:stellar_specific_am`         -> Norm of the stellar specific angular momentum.
      + `:gas_specific_am`             -> Norm of the gas specific angular momentum.
      + `:dm_specific_am`              -> Norm of the dark matter specific angular momentum.
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
      + `:scale_factor`                -> Scale factor.
      + `:redshift`                    -> Redshift.
      + `:physical_time`               -> Physical time since the Big Bang.
      + `:lookback_time`               -> Physical time left to reach the last snapshot.

# Returns

  - A [`PlotParams`](@ref) object, with entries:

      + `request::Dict{Symbol,Vector{String}}` -> Data request for [`readSnapshot`](@ref).
      + `var_name::AbstractString`             -> Name of the quantity for the plot axis.
      + `exp_factor::Int`                      -> Numerical exponent to scale down the axis.
      + `unit::Unitful.Units`                  -> Target unit for the axis.
      + `axis_label::AbstractString`           -> Label for the axis.
"""
function plotParams(quantity::Symbol)::PlotParams

    if quantity == :stellar_mass

        plot_params = PlotParams(;
            request    = Dict(:stars => ["MASS", "POS "]),
            var_name   = L"M_\star",
            exp_factor = 10,
            unit       = u"Msun",
        )

    elseif quantity == :gas_mass

        plot_params = PlotParams(;
            request    = Dict(:gas => ["MASS", "POS ", "RHO "]),
            var_name   = L"M_\mathrm{gas}",
            exp_factor = 10,
            unit       = u"Msun",
        )

    elseif quantity == :hydrogen_mass

        plot_params = PlotParams(;
            request    = Dict(:gas => ["MASS", "POS ", "RHO "]),
            var_name   = L"M_\mathrm{H}",
            exp_factor = 10,
            unit       = u"Msun",
        )

    elseif quantity == :dm_mass

        plot_params = PlotParams(;
            request    = Dict(:halo => ["MASS", "POS "]),
            var_name   = L"M_\mathrm{DM}",
            exp_factor = 10,
            unit       = u"Msun",
        )

    elseif quantity == :bh_mass

        plot_params = PlotParams(;
            request  = Dict(:black_hole => ["MASS", "POS "]),
            var_name = L"M_\mathrm{BH}",
            unit     = u"Msun",
        )

    elseif quantity == :molecular_mass

        plot_params = PlotParams(;
            request    = Dict(
                :gas => ["MASS", "POS ", "FRAC", "RHO ", "CTIM", "TAUS"],
            ),
            var_name   = L"M_\mathrm{H_2}",
            exp_factor = 10,
            unit       = u"Msun",
        )

    elseif quantity == :br_molecular_mass

        plot_params = PlotParams(;
            request    = Dict(
                :gas => ["MASS", "PRES", "RHO "],
            ),
            var_name   = L"M_\mathrm{H_2^{BR}}",
            exp_factor = 10,
            unit       = u"Msun",
        )

    elseif quantity == :atomic_mass

        plot_params = PlotParams(;
            request    = Dict(
                :gas => ["MASS", "POS ", "FRAC", "NH  ", "NHP ", "RHO ", "CTIM", "TAUS"],
            ),
            var_name   = L"M_\mathrm{HI}",
            exp_factor = 10,
            unit       = u"Msun",
        )

    elseif quantity == :ionized_mass

        plot_params = PlotParams(;
            request    = Dict(
                :gas => ["MASS", "POS ", "FRAC", "NH  ", "NHP ", "RHO ", "CTIM", "TAUS"],
            ),
            var_name   = L"M_\mathrm{HII}",
            exp_factor = 10,
            unit       = u"Msun",
        )

    elseif quantity == :neutral_mass

        plot_params = PlotParams(;
            request    = Dict(
                :gas => ["MASS", "POS ", "FRAC", "NH  ", "NHP ", "RHO ", "CTIM", "TAUS"],
            ),
            var_name   = L"M_\mathrm{HI + H_2}",
            exp_factor = 10,
            unit       = u"Msun",
        )

    elseif quantity == :generic_mass

        plot_params = PlotParams(;
            request    = Dict(
                :stars   => ["MASS", "POS "],
                :gas     => ["MASS", "POS ", "FRAC", "NH  ", "NHP ", "PRES", "RHO ", "CTIM", "TAUS"],
                :dm_mass => ["MASS", "POS "],
                :bh_mass => ["MASS", "POS "],
            ),
            var_name   = L"M",
            exp_factor = 10,
            unit       = u"Msun",
        )

    elseif quantity == :stellar_number

        plot_params = PlotParams(;
            request  = Dict(:stars => ["MASS", "POS "]),
            var_name = L"\mathrm{Number \,\, of \,\, stellar \,\, particles}",
        )

    elseif quantity == :gas_number

        plot_params = PlotParams(;
            request  = Dict(:gas => ["MASS", "POS "]),
            var_name = L"\mathrm{Number \,\, of \,\, gas \,\, cells}",
        )

    elseif quantity == :dm_number

        plot_params = PlotParams(;
            request  = Dict(:halo => ["MASS", "POS "]),
            var_name = L"\mathrm{Number \,\, of \,\, DM \,\, particles}",
        )

    elseif quantity == :bh_number

        plot_params = PlotParams(;
            request  = Dict(:black_hole => ["MASS", "POS "]),
            var_name = L"\mathrm{Number \,\, of \,\, BH \,\, particles}",
        )

    elseif quantity == :molecular_fraction

        plot_params = PlotParams(;
            request  = Dict(
                :gas => ["MASS", "POS ", "FRAC", "CTIM", "TAUS"],
            ),
            var_name = L"f_\mathrm{H_2}",
        )

    elseif quantity == :br_molecular_fraction

        plot_params = PlotParams(;
            request  = Dict(
                :gas => ["MASS", "POS ", "PRES"],
            ),
            var_name = L"f_\mathrm{H_2}^\mathrm{BR}",
        )

    elseif quantity == :atomic_fraction

        plot_params = PlotParams(;
            request  = Dict(
                :gas => ["MASS", "POS ", "FRAC", "NH  ", "NHP ", "CTIM", "TAUS"],
            ),
            var_name = L"f_\mathrm{HI}",
        )

    elseif quantity == :ionized_fraction

        plot_params = PlotParams(;
            request  = Dict(
                :gas => ["MASS", "POS ", "FRAC", "NH  ", "NHP ", "CTIM", "TAUS"],
            ),
            var_name = L"f_\mathrm{HII}",
        )

    elseif quantity == :neutral_fraction

        plot_params = PlotParams(;
            request  = Dict(
                :gas => ["MASS", "POS ", "FRAC", "NH  ", "NHP ", "CTIM", "TAUS"],
            ),
            var_name = L"f_\mathrm{H_I + H_2}",
        )

    elseif quantity == :molecular_neutral_fraction

        plot_params = PlotParams(;
            request  = Dict(
                :gas => ["MASS", "POS ", "FRAC", "NH  ", "NHP ", "PRES", "RHO ", "CTIM", "TAUS"],
            ),
            var_name = L"f_\mathrm{H_2}^\star",
        )

    elseif quantity == :mol_eq_quotient

        plot_params = PlotParams(;
            request  = Dict(
                :gas => ["ETAD", "FRAC", "RHOC", "PARZ"],
            ),
            var_name = L"\log_{10} \, \mathrm{LS^{H_2} / RS^{H_2}}",
        )

    elseif quantity == :ion_eq_quotient

        plot_params = PlotParams(;
            request  = Dict(
                :gas => ["ETAI", "PARR", "FRAC", "RHOC"],
            ),
            var_name = L"\log_{10} \, \mathrm{LS^{HII} / RS^{HII}}",
        )

    elseif quantity == :generic_fraction

        plot_params = PlotParams(;
            request  = Dict(
                :gas => ["MASS", "POS ", "FRAC", "NH  ", "NHP ", "PRES", "RHO ", "CTIM", "TAUS"],
            ),
            var_name = L"f",
        )

    elseif quantity == :gas_mass_density

        plot_params = PlotParams(;
            request  = Dict(:gas => ["MASS", "POS ", "RHO "]),
            var_name = L"\rho_\mathrm{gas}",
            unit     = u"Msun*kpc^-3",
        )

    elseif quantity == :hydrogen_mass_density

        plot_params = PlotParams(;
            request  = Dict(:gas => ["MASS", "POS ", "RHO "]),
            var_name = L"\rho_\mathrm{H}",
            unit     = u"Msun*kpc^-3",
        )

    elseif quantity == :gas_number_density

        plot_params = PlotParams(;
            request  = Dict(:gas => ["MASS", "POS ", "RHO "]),
            var_name = L"n_\mathrm{gas}",
            unit     = u"cm^-3",
        )

    elseif quantity == :molecular_number_density

        plot_params = PlotParams(;
            request  = Dict(
                :gas => ["MASS", "POS ", "FRAC", "RHO ", "CTIM", "TAUS"],
            ),
            var_name = L"n_\mathrm{H_2}",
            unit     = u"cm^-3",
        )

    elseif quantity == :br_molecular_number_density

        plot_params = PlotParams(;
            request  = Dict(
                :gas => ["MASS", "POS ", "PRES"],
            ),
            var_name = L"n_\mathrm{H_2}^{BR}",
            unit     = u"cm^-3",
        )

    elseif quantity == :atomic_number_density

        plot_params = PlotParams(;
            request  = Dict(
                :gas => ["MASS", "POS ", "FRAC", "NH  ", "NHP ", "RHO ", "CTIM", "TAUS"],
            ),
            var_name = L"n_\mathrm{HI}",
            unit     = u"cm^-3",
        )

    elseif quantity == :ionized_number_density

        plot_params = PlotParams(;
            request  = Dict(
                :gas => ["MASS", "POS ", "FRAC", "NH  ", "NHP ", "RHO ", "CTIM", "TAUS"],
            ),
            var_name = L"n_\mathrm{HII}",
            unit     = u"cm^-3",
        )

    elseif quantity == :neutral_number_density

        plot_params = PlotParams(;
            request  = Dict(
                :gas => ["MASS", "POS ", "FRAC", "NH  ", "NHP ", "RHO ", "CTIM", "TAUS"],
            ),
            var_name = L"n_\mathrm{HI + H_2}",
            unit     = u"cm^-3",
        )

    elseif quantity == :stellar_area_density

        plot_params = PlotParams(;
            request  = Dict(:stars => ["MASS", "POS "]),
            var_name = L"\Sigma_\star",
            unit     = u"Msun*pc^-2",
        )

    elseif quantity == :gas_area_density

        plot_params = PlotParams(;
            request  = Dict(:gas => ["MASS", "POS ", "RHO "]),
            var_name = L"\Sigma_\mathrm{gas}",
            unit     = u"Msun*pc^-2",
        )

    elseif quantity == :molecular_area_density

        plot_params = PlotParams(;
            request  = Dict(
                :gas => ["MASS", "POS ", "FRAC", "RHO ", "CTIM", "TAUS"],
            ),
            var_name = L"\Sigma_\mathrm{H_2}",
            unit     = u"Msun*pc^-2",
        )

    elseif quantity == :br_molecular_area_density

        plot_params = PlotParams(;
            request  = Dict(
                :gas => ["MASS", "POS ", "PRES", "RHO "],
            ),
            var_name = L"\Sigma_\mathrm{H_2}^\mathrm{BR}",
            unit     = u"Msun*pc^-2",
        )

    elseif quantity == :atomic_area_density

        plot_params = PlotParams(;
            request  = Dict(
                :gas => ["MASS", "POS ", "FRAC", "NH  ", "NHP ", "RHO ", "CTIM", "TAUS"],
            ),
            var_name = L"\Sigma_\mathrm{HI}",
            unit     = u"Msun*pc^-2",
        )

    elseif quantity == :ionized_area_density

        plot_params = PlotParams(;
            request  = Dict(
                :gas => ["MASS", "POS ", "FRAC", "NH  ", "NHP ", "RHO ", "CTIM", "TAUS"],
            ),
            var_name = L"\Sigma_\mathrm{HII}",
            unit     = u"Msun*pc^-2",
        )

    elseif quantity == :neutral_area_density

        plot_params = PlotParams(;
            request  = Dict(
                :gas => ["MASS", "POS ", "FRAC", "NH  ", "NHP ", "RHO ", "CTIM", "TAUS"],
            ),
            var_name = L"\Sigma_\mathrm{HI + H_2}",
            unit     = u"Msun*pc^-2",
        )

    elseif quantity == :sfr_area_density

        plot_params = PlotParams(;
            request  = Dict(:stars => ["MASS", "POS ", "GAGE"]),
            var_name = L"\Sigma_\mathrm{SFR}",
            unit     = u"Msun*yr^-1*kpc^-2",
        )

    elseif quantity == :generic_area_density

        plot_params = PlotParams(;
            request  = Dict(
                :stars => ["MASS", "POS ", "GAGE"],
                :gas => ["MASS", "POS ", "FRAC", "NH  ", "NHP ", "RHO ", "PRES", "CTIM", "TAUS"],
            ),
            var_name = L"\Sigma",
            unit     = u"Msun*pc^-2",
        )

    elseif quantity == :gas_metallicity

        plot_params = PlotParams(;
            request  = Dict(:gas => ["MASS", "POS ", "GMET", "GZ  "]),
            var_name = L"Z_\mathrm{gas} \, [\mathrm{Z_\odot}]",
        )

    elseif quantity == :stellar_metallicity

        plot_params = PlotParams(;
            request  = Dict(:stars => ["MASS", "POS ", "GME2", "GZ2 "]),
            var_name = L"Z_\star \, [\mathrm{Z_\odot}]",
        )

    elseif quantity ∈ GAS_ABUNDANCE

        element_string = first(split(string(quantity), "_"))

        plot_params = PlotParams(;
            request    = Dict(:gas => ["MASS", "POS ", "GMET"]),
            axis_label = L"12 + \log_{10}(\mathrm{%$element_string} \, / \, \mathrm{H})",
        )

    elseif quantity ∈ STELLAR_ABUNDANCE

        element_string = first(split(string(quantity), "_"))

        plot_params = PlotParams(;
            request    = Dict(:stars => ["MASS", "POS ", "GME2"]),
            axis_label = L"12 + \log_{10}(\mathrm{%$element_string} \, / \, \mathrm{H})",
        )

    elseif quantity == :stellar_radial_distance

        plot_params = PlotParams(;
            request  = Dict(:stars => ["POS "]),
            var_name = L"r",
            unit     = u"kpc",
        )

    elseif quantity == :gas_radial_distance

        plot_params = PlotParams(;
            request  = Dict(:gas => ["POS "]),
            var_name = L"r",
            unit     = u"kpc",
        )

    elseif quantity == :dm_radial_distance

        plot_params = PlotParams(;
            request  = Dict(:halo => ["POS "]),
            var_name = L"r",
            unit     = u"kpc",
        )

    elseif quantity == :stellar_xy_distance

        plot_params = PlotParams(;
            request  = Dict(:stars => ["POS "]),
            var_name = L"r_{xy}",
            unit     = u"kpc",
        )

    elseif quantity == :gas_xy_distance

        plot_params = PlotParams(;
            request  = Dict(:gas => ["POS "]),
            var_name = L"r_{xy}",
            unit     = u"kpc",
        )

    elseif quantity == :dm_xy_distance

        plot_params = PlotParams(;
            request  = Dict(:halo => ["POS "]),
            var_name = L"r_{xy}",
            unit     = u"kpc",
        )

    elseif quantity == :gas_sfr

        plot_params = PlotParams(;
            request  = Dict(:gas => ["SFR "]),
            var_name = L"\mathrm{SFR_{gas}}",
            unit     = u"Msun*yr^-1",
        )

    elseif quantity == :mass_accretion

        plot_params = PlotParams(;
            request  = Dict(
                :gas         => ["ID  ", "MASS"],
                :stars       => ["ID  ", "MASS"],
                :black_hole  => ["ID  ", "MASS"],
                :group       => ["G_R_Crit200", "G_M_Crit200"],
                :tracer      => ["PAID", "TRID"],
            ),
            var_name = "Mass accretion",
            unit     = u"Msun*yr^-1",
        )

    elseif quantity == :stellar_specific_am

        plot_params = PlotParams(;
            request  = Dict(:stars => ["MASS", "POS ", "VEL"]),
            var_name = L"j_\star",
            unit     = u"kpc^2*s^-1",
        )

    elseif quantity == :gas_specific_am

        plot_params = PlotParams(;
            request  = Dict(:gas => ["MASS", "POS ", "VEL"]),
            var_name = L"j_\mathrm{gas}",
            unit     = u"kpc^2*s^-1",
        )

    elseif quantity == :dm_specific_am

        plot_params = PlotParams(;
            request  = Dict(:halo => ["MASS", "POS ", "VEL"]),
            var_name = L"j_\mathrm{DM}",
            unit     = u"kpc^2*s^-1",
        )

    elseif quantity == :stellar_circularity

        # `daBandProfile` expects that the first element in the request dictionary is for the stars
        plot_params = PlotParams(;
            request  = Dict(
                :stars      => ["MASS", "POS ", "VEL "],
                :gas        => ["MASS", "POS "],
                :halo       => ["MASS", "POS "],
                :black_hole => ["MASS", "POS "],
            ),
            var_name = L"\epsilon",
        )

    elseif quantity == :stellar_vcirc

        # `daBandProfile` expects that the first element in the request dictionary is for the stars
        plot_params = PlotParams(;
            request  = Dict(
                :stars      => ["MASS", "POS "],
                :gas        => ["MASS", "POS "],
                :halo       => ["MASS", "POS "],
                :black_hole => ["MASS", "POS "],
            ),
            var_name = L"v_\mathrm{circ}",
            unit     = u"km*s^-1",
        )

    elseif quantity == :stellar_vradial

        plot_params = PlotParams(;
            request  = Dict(:stars => ["POS ", "VEL "]),
            var_name = L"v_r",
            unit     = u"km*s^-1",
        )

    elseif quantity == :stellar_vtangential

        plot_params = PlotParams(;
            request  = Dict(:stars => ["POS ", "VEL "]),
            var_name = L"v_\theta",
            unit     = u"km*s^-1",
        )

    elseif quantity == :stellar_vzstar

        plot_params = PlotParams(;
            request  = Dict(:stars => ["POS ", "VEL "]),
            var_name = L"v_z \,\, \mathrm{sign}(z)",
            unit     = u"km*s^-1",
        )

    elseif quantity == :stellar_age

        plot_params = PlotParams(;
            request  = Dict(:stars => ["GAGE"]),
            var_name = L"\mathrm{Stellar \,\, age}",
            unit     = u"Gyr",
        )

    elseif quantity == :sfr

        plot_params = PlotParams(;
            request  = Dict(:stars => ["MASS", "POS ", "GAGE"]),
            var_name = L"SFR",
            unit     = u"Msun*yr^-1",
        )

    elseif quantity == :ssfr

        plot_params = PlotParams(;
            request  = Dict(:stars => ["MASS", "POS ", "GAGE"]),
            var_name = L"sSFR",
            unit     = u"yr^-1",
        )

    elseif quantity == :observational_sfr

        plot_params = PlotParams(;
            request  = Dict(:stars => ["MASS", "POS ", "GAGE"]),
            var_name = L"SFR",
            unit     = u"Msun*yr^-1",
        )

    elseif quantity == :observational_ssfr

        plot_params = PlotParams(;
            request  = Dict(:stars => ["MASS", "POS ", "GAGE"]),
            var_name = L"sSFR",
            unit     = u"yr^-1",
        )

    elseif quantity == :temperature

        plot_params = PlotParams(;
            request    = Dict(:gas => ["MASS", "POS ", "TEMP"]),
            axis_label = L"\log_{10} \, T \, [\mathrm{K}]",
        )

    elseif quantity == :pressure

        plot_params = PlotParams(;
            request    = Dict(:gas => ["PRES"]),
            var_name   = L"P",
            exp_factor = -13,
            unit       = u"Pa",
        )

    elseif quantity == :scale_factor

        plot_params = PlotParams(;
            var_name = L"a",
        )

    elseif quantity == :redshift

        plot_params = PlotParams(;
            var_name = L"z",
        )

    elseif quantity == :physical_time

        plot_params = PlotParams(;
            var_name = L"t",
            unit = u"Gyr",
        )

    elseif quantity == :lookback_time

        plot_params = PlotParams(;
            var_name = L"\mathrm{Lookback \,\, time}",
            unit     = u"Gyr",
        )

    elseif quantity == :time_step

        plot_params = PlotParams(;
            var_name = L"\mathrm{Number \,\, of \,\, time \,\, steps}",
        )

    elseif quantity == :clock_time_s

        plot_params = PlotParams(;
            var_name = L"\mathrm{Wallclock \,\, time}",
            unit     = u"s",
        )

    elseif quantity == :clock_time_percent

        plot_params = PlotParams(;
            axis_label = L"\mathrm{Wallclock \,\, time \, [\%]}",
        )

    elseif quantity == :tot_clock_time_s

        plot_params = PlotParams(;
            var_name = L"\mathrm{Cumulative \,\, wallclock \,\, time}",
            unit     = u"s",
        )

    elseif quantity == :tot_clock_time_percent

        plot_params = PlotParams(;
            axis_label = L"\mathrm{Cumulative \,\, wallclock \,\, time \, [\%]}",
        )

    else

        throw(ArgumentError("plotParams: I don't recognize the quantity :$(quantity)"))

    end

    return plot_params

end

"""
    isSubfindActive(path::String)::Bool

Check if there is information about the halos and subhalos in the group catalog file.

# Arguments

  - `path::String`: Path to the group catalog file or folder.

# Returns

  - If there are halo and subhalo information in the group catalog file.
"""
function isSubfindActive(path::String)::Bool

    if isfile(path)

        (
            HDF5.ishdf5(path) ||
            throw(ArgumentError("isSubfindActive: The file $(path) is not in the HDF5 format, \
            I don't know how to read it"))
        )

        file_path = path

    elseif isdir(path)

        sub_files = glob("$(GC_BASENAME)_*.*.hdf5", path)

        (
            !isempty(sub_files) && all(HDF5.ishdf5, sub_files) ||
            throw(ArgumentError("isSubfindActive: The directory $(path) does not contain \
            group catalog sub-files in the HDF5 format"))
        )

        file_path = minimum(sub_files)

    else

        throw(ArgumentError("isSubfindActive: $(path) does not exist as a file or folder"))

    end

    subfind_active = h5open(file_path, "r") do gc_file

        (
            all(in(keys(gc_file)), ["Group", "Subhalo"]) &&
            all(!isempty, [gc_file["Group"], gc_file["Subhalo"]])
        )

    end

    return subfind_active

end

"""
    mergeRequests(requests::Dict{Symbol,Vector{String}}...)::Dict{Symbol,Vector{String}}

Merge several request dictionaries, ignoring duplicates.

# Arguments

  - `requests`: The request dictionaries for [`readSnapshot`](@ref).

# Returns

  - A new dictionary with all the requests.
"""
function mergeRequests(requests::Dict{Symbol,Vector{String}}...)::Dict{Symbol,Vector{String}}

    return Dict(
        type => union([get(request, type, String[]) for request in requests]...) for
        type in union(keys.(requests)...)
    )

end

"""
    addRequest(
        request::Dict{Symbol,Vector{String}},
        addition::Dict{Symbol,Vector{String}},
    )::Dict{Symbol,Vector{String}}

Add the blocks in `addition` to `request`, only for the types already present in `request`.

# Arguments

  - `request::Dict{Symbol,Vector{String}}`: The request dictionary for [`readSnapshot`](@ref).
  - `addition::Dict{Symbol,Vector{String}}`: Request dictionary with the blocks to be added, only for the types already present in `request`.

# Returns

  - A new dictionary with all the requests.
"""
function addRequest(
    request::Dict{Symbol,Vector{String}},
    addition::Dict{Symbol,Vector{String}},
)::Dict{Symbol,Vector{String}}

    return Dict(type => blocks ∪ get(addition, type, String[]) for (type, blocks) in request)

end

@doc raw"""
    bigiel2008(
        ΣH::Vector{<:SurfaceDensity};
        <keyword arguments>
    )::Vector{<:Number}

Kennicutt-Schmidt law for the molecular or neutral gas, taken from a set of observations of nearby galaxies.

From Bigiel et al. (2008) (Section 3.1), we have

```math
\Sigma_\mathrm{SFR} = a \left( \frac{\Sigma_\mathrm{HI, H_2, gas}}{10 \, \mathrm{M_\odot \, pc^{-2}}} \right)^{\!N} \, ,
```
where N is the power-law index, and $A = \log_{10}(a)$, where $a$ is $\Sigma_\mathrm{SFR}$ at the fiducial gas surface density of $10 \, \mathrm{M_\odot \, pc^{-2}}$.

# Arguments

  - `ΣH::Vector{<:SurfaceDensity}`: Values of the molecular or neutral gas surface density, with units.
  - `molecular::Bool=true`: If the x axis will be the area mass density of molecular hydrogen, or, if set to false, the area mass density of neutral hydrogen.
  - `log_output::Bool=true`: If the output will the $\log_{10}$ of the star formation area density, or the star formation area density itself (with units). If `log_output` = true, the implied unit is $\mathrm{M_\odot \, yr^{-1} \, kpc^{-2}}$

# Returns

  - The star formation area density.

# References

F. Bigiel et al. (2008). *THE STAR FORMATION LAW IN NEARBY GALAXIES ON SUB-KPC SCALES*. The Astrophysical Journal, **136(6)**, 2846. [doi:10.1088/0004-6256/136/6/2846](https://doi.org/10.1088/0004-6256/136/6/2846)
"""
function bigiel2008(
    ΣH::Vector{<:SurfaceDensity};
    molecular::Bool=true,
    log_output::Bool=true,
)::Vector{<:Number}

    log10ΣH = @. log10(uconvert(Unitful.NoUnits, ΣH / 10.0u"Msun * pc^-2"))

    if molecular
        log10Σsfr = @. A_BIGIEL2008_BF_MOLECULAR + log10ΣH * N_BIGIEL2008_BF_MOLECULAR
    else
        log10Σsfr = @. A_BIGIEL2008_NEUTRAL + log10ΣH * N_BIGIEL2008_NEUTRAL
    end

    if log_output
        return log10Σsfr
    else
        return @. exp10(log10Σsfr ) * u"Msun * yr^-1 * kpc^-2"
    end

end

@doc raw"""
    invBigiel2008(
        Σsfr ::Vector{<:MassFlowDensity};
        <keyword arguments>
    )::Vector{<:Number}

Inverse Kennicutt-Schmidt law for the molecular or neutral gas, taken from a set of observations of nearby galaxies.

From Bigiel et al. (2008) (Section 3.1, Eq. 2), we have

```math
\Sigma_\mathrm{SFR} = a \left( \frac{\Sigma_\mathrm{HI, H_2, gas}}{10 \, \mathrm{M_\odot \, pc^{-2}}} \right)^{\!N} \, ,
```
where N is the power-law index, and $A = \log_{10}(a)$, where $a$ is $\Sigma_\mathrm{SFR}$ at the fiducial gas surface density of $10 \, \mathrm{M_\odot \, pc^{-2}}$.

# Arguments

  - `Σsfr ::Vector{<:MassFlowDensity}`: Values of the star formation area density, with units.
  - `molecular::Bool=true`: If the output will be the area mass density of molecular hydrogen, or, if set to false, the area mass density of neutral hydrogen.
  - `log_output::Bool=true`: If the output will the $\log_{10}$ of the molecular or neutral gas surface density, or the molecular or neutral gas surface density itself (with units). If `log_output` = true, the implied unit is $10 \, \mathrm{M_\odot \, pc^{-2}}$

# Returns

  - The molecular or neutral gas surface density.

# References

F. Bigiel et al. (2008). *THE STAR FORMATION LAW IN NEARBY GALAXIES ON SUB-KPC SCALES*. The Astrophysical Journal, **136(6)**, 2846. [doi:10.1088/0004-6256/136/6/2846](https://doi.org/10.1088/0004-6256/136/6/2846)
"""
function invBigiel2008(
    Σsfr ::Vector{<:MassFlowDensity};
    molecular::Bool=true,
    log_output::Bool=true,
)::Vector{<:Number}

    log10Σsfr = @. log10(ustrip(u"Msun * yr^-1 * kpc^-2", Σsfr ))

    if molecular
        log10ΣH = @. (log10Σsfr - A_BIGIEL2008_BF_MOLECULAR) / N_BIGIEL2008_BF_MOLECULAR
    else
        log10ΣH = @. (log10Σsfr - A_BIGIEL2008_NEUTRAL) / N_BIGIEL2008_NEUTRAL
    end

    if log_output
        return log10ΣH
    else
        return @. exp10(log10ΣH) * 10.0u"Msun * pc^-2"
    end

end

@doc raw"""
    kennicutt1998(Σgas::Vector{<:SurfaceDensity}; <keyword arguments>)::Vector{<:Number}

Kennicutt-Schmidt law, taken from a set of observations of nearby galaxies.

From Kennicutt (1998) (Section 4, Eq. 4), we have

```math
\Sigma_\mathrm{SFR} = a \left( \frac{\Sigma_\mathrm{gas}}{1 \, \mathrm{M_\odot \, pc^{-2}}} \right)^{\!N} \mathrm{M_\odot \, yr^{-1] \, kpc^{-2}} \, ,
```
where N is the power-law index and $a$ is $\Sigma_\mathrm{SFR}$ at the fiducial gas surface density of $1 \, \mathrm{M_\odot \, pc^{-2}}$.

# Arguments

  - `Σgas::Vector{<:SurfaceDensity}`: Values of the gas mass surface density, with units.
  - `log_output::Bool=true`: If the output will the $\log_{10}$ of the star formation area density, or the star formation area density itself (with units). If `log_output` = true, the implied unit is $\mathrm{M_\odot \, yr^{-1} \, kpc^{-2}}$

# Returns

  - The star formation area density.

# References

R. C. Kennicutt (1998). *The Global Schmidt Law in Star-forming Galaxies*. The Astrophysical Journal, **498(2)**, 541-552. [doi:10.1086/305588](https://doi.org/10.1086/305588)
"""
function kennicutt1998(Σgas::Vector{<:SurfaceDensity}; log_output::Bool=true)::Vector{<:Number}

    log10Σgas = @. log10(ustrip(u"Msun * pc^-2", Σgas))
    log10Σsfr = @. log10(a_KS98) + log10Σgas * N_KS98

    if log_output
        return log10Σsfr
    else
        return @. exp10(log10Σsfr ) * u"Msun * yr^-1 * kpc^-2"
    end

end

@doc raw"""
    invKennicutt1998(Σsfr::Vector{<:MassFlowDensity}; <keyword arguments>)::Vector{<:Number}

Inverse Kennicutt-Schmidt law, taken from a set of observations of nearby galaxies.

From Kennicutt (1998) (Section 4, Eq. 4), we have

```math
\Sigma_\mathrm{SFR} = a \left( \frac{\Sigma_\mathrm{gas}}{1 \, \mathrm{M_\odot \, pc^{-2}}} \right)^{\!N} \mathrm{M_\odot \, yr^{-1] \, kpc^{-2}} \, ,
```
where N is the power-law index and $a$ is $\Sigma_\mathrm{SFR}$ at the fiducial gas surface density of $1 \, \mathrm{M_\odot \, pc^{-2}}$.

# Arguments

  - `Σsfr::Vector{<:MassFlowDensity}`: Values of the star formation area density, with units.
  - `log_output::Bool=true`: If the output will the $\log_{10}$ of the gas mass surface density, or the gas mass surface density itself (with units). If `log_output` = true, the implied unit is $\mathrm{M_\odot \, pc^{-2}}$

# Returns

  - The gas mass surface density.

# References

R. C. Kennicutt (1998). *The Global Schmidt Law in Star-forming Galaxies*. The Astrophysical Journal, **498(2)**, 541-552. [doi:10.1086/305588](https://doi.org/10.1086/305588)
"""
function invKennicutt1998(Σsfr::Vector{<:MassFlowDensity}; log_output::Bool=true)::Vector{<:Number}

    log10Σsfr = @. log10(ustrip(u"Msun * yr^-1 * kpc^-2", Σsfr))
    log10Σgas = @. (log10Σsfr - log10(a_KS98)) / N_KS98

    if log_output
        return log10Σgas
    else
        return @. exp10(log10Σgas) * u"Msun * pc^-2"
    end

end
