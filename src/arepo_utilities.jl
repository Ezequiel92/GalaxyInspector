####################################################################################################
# Arepo and plotting utilities.
####################################################################################################

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
        throw(DomainError("formatError: `q_error` must be positive, but I got \
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
        `QUANTITIES` but I got $(quantity), see the options in `./src/constants.jl`"))
    )

    header = readSnapHeader(path)
    cosmological = isCosmological(path)

    # Set up the struct for unit conversion
    IU = InternalUnits(;
        l_unit=header.l_unit,
        m_unit=header.m_unit,
        v_unit=header.v_unit,
        a0=cosmological ? header.time : 1.0,
        h0=cosmological ? header.h0 : 1.0,
    )

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
                    in `constants.jl`",
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
    filterData!(data_dict::Dict; <keyword arguments>)::Nothing

Filter `data_dict` using the indices provided by `filter_function`.

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
  - `filter_function::Function=filterNothing`: A functions with the signature:

    `filter_function(data_dict) -> indices`

    where

      + `data_dict::Dict`: A dictionary with the following shape:

          * `:sim_data`          -> ::Simulation (see [`Simulation`](@ref)).
          * `:snap_data`         -> ::Snapshot (see [`Snapshot`](@ref)).
          * `:gc_data`           -> ::GroupCatalog (see [`GroupCatalog`](@ref)).
          * `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
          * `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
          * `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
          * ...
          * `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
          * `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
          * `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
          * ...

      + `indices::Dict`: A dictionary with the following shape:

          * `cell/particle type` -> idxs::IndexType
          * `cell/particle type` -> idxs::IndexType
          * `cell/particle type` -> idxs::IndexType
          * ...
"""
function filterData!(data_dict::Dict; filter_function::Function=filterNothing)::Nothing

    # Compute the filter dictionary
    indices = filter_function(data_dict)

    @inbounds for type_symbol in snapshotTypes(data_dict)

        idxs = indices[type_symbol]

        @inbounds for (block, values) in data_dict[type_symbol]
            @inbounds if !isempty(values)
                data_dict[type_symbol][block] = collect(selectdim(values, ndims(values), idxs))
            end
        end

    end

    return nothing

end

"""
    filterData(data_dict::Dict; <keyword arguments>)::Dict

Returna filtered copy of `data_dict` using the indices provided by `filter_function`.

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
  - `filter_function::Function=filterNothing`: A functions with the signature:

    `filter_function(data_dict) -> indices`

    where

      + `data_dict::Dict`: A dictionary with the following shape:

          * `:sim_data`          -> ::Simulation (see [`Simulation`](@ref)).
          * `:snap_data`         -> ::Snapshot (see [`Snapshot`](@ref)).
          * `:gc_data`           -> ::GroupCatalog (see [`GroupCatalog`](@ref)).
          * `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
          * `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
          * `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
          * ...
          * `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
          * `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
          * `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
          * ...

      + `indices::Dict`: A dictionary with the following shape:

          * `cell/particle type` -> idxs::IndexType
          * `cell/particle type` -> idxs::IndexType
          * `cell/particle type` -> idxs::IndexType
          * ...

# Returns

  - The filtered data.
"""
function filterData(data_dict::Dict; filter_function::Function=filterNothing)::Dict

    dd_copy = deepcopy(data_dict)

    # Compute the filter dictionary
    indices = filter_function(dd_copy)

    @inbounds for type_symbol in snapshotTypes(dd_copy)

        idxs = indices[type_symbol]

        @inbounds for (block, values) in dd_copy[type_symbol]
            @inbounds if !isempty(values)
                dd_copy[type_symbol][block] = collect(selectdim(values, ndims(values), idxs))
            end
        end

    end

    return dd_copy

end

"""
    computeCenter(data_dict::Dict, subfind_idx::NTuple{2,Int})::Vector{<:Unitful.Length}

Read the position of the potencial minimum for a given halo or subhalo.

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
  - `subfind_idx::NTuple{2,Int}`: Tuple with two elements:

      + Index of the target halo (FoF group). Starts at 1.
      + Index of the target subhalo (subfind), relative the target halo. Starts at 1. If set to 0, the potencial minimum of the whole halo with index `halo_idx` is returned.

# Returns

  - The specified potencial minimum.
"""
function computeCenter(data_dict::Dict, subfind_idx::NTuple{2,Int})::Vector{<:Unitful.Length}

    # If there are no subfind data, return the origin
    if ismissing(data_dict[:gc_data].path) && !isSubfindActive(data_dict[:gc_data].path)
        return zeros(typeof(1.0u"kpc"), 3)
    end

    halo_idx, subhalo_rel_idx = subfind_idx

    # Load the necessary data
    g_n_subs = data_dict[:group]["G_Nsubs"]
    g_pos = data_dict[:group]["G_Pos"]
    s_pos = data_dict[:subhalo]["S_Pos"]

    # Check that the requested halo index is within bounds
    n_groups_total = data_dict[:gc_data].header.n_groups_total

    (
        !iszero(n_groups_total) && !any(isempty, [g_n_subs, g_pos, s_pos]) ||
        return zeros(typeof(1.0u"kpc"), 3)
    )

    (
        0 < halo_idx <= n_groups_total ||
        throw(ArgumentError("computeCenter: There is only $(n_groups_total) FoF goups in \
        $(data_dict[:gc_data].path), so halo_idx = $(halo_idx) is out of bounds"))
    )

    # Select the halo potencial minimum if `subhalo_rel_idx` == 0
    isPositive(subhalo_rel_idx) || return g_pos[:, halo_idx]

    # Check that the requested subhalo index is within bounds
    n_subfinds = g_n_subs[halo_idx]
    (
        subhalo_rel_idx <= n_subfinds ||
        throw(ArgumentError("computeCenter: There is only $(n_subfinds) subhalos for the FoF \
        group $(halo_idx) in $(data_dict[:gc_data].path), so subhalo_rel_idx = \
        $(subhalo_rel_idx) is out of bounds"))
    )

    # Compute the number of subhalos and particles up to the last halo before `halo_idx`
    if isone(halo_idx)
        n_subs_floor = 0
    else
        n_subs_floor = sum(g_n_subs[1:(halo_idx - 1)]; init=0)
    end

    # Compute the subhalo absolute index
    subhalo_abs_idx = n_subs_floor + subhalo_rel_idx

    # Select the subhalo potencial minimum
    return s_pos[:, subhalo_abs_idx]

end

"""
    computeCenter(data_dict::Dict, subhalo_abs_idx::Int)::Vector{<:Unitful.Length}

Read the position of the potencial minimum for a given subhalo.

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
  - `subhalo_abs_idx::Int`: Absolute index of the target subhalo (subfind). Starts at 1.

# Returns

  - The specified potencial minimum.
"""
function computeCenter(data_dict::Dict, subhalo_abs_idx::Int)::Vector{<:Unitful.Length}

    # If there are no subfind data, return the origin
    if ismissing(data_dict[:gc_data].path) && !isSubfindActive(data_dict[:gc_data].path)
        return zeros(typeof(1.0u"kpc"), 3)
    end

    s_pos = data_dict[:subhalo]["S_Pos"]

    # Check that the requested subhalo index is within bounds
    n_subgroups_total = data_dict[:gc_data].header.n_subgroups_total

    !iszero(n_subgroups_total) && !isempty(s_pos) || return zeros(typeof(1.0u"kpc"), 3)

    (
        0 < subhalo_abs_idx <= n_subgroups_total ||
        throw(ArgumentError("computeCenter: There is only $(n_subgroups_total) subhalos in \
        $(data_dict[:gc_data].path), so subhalo_abs_idx = $(subhalo_abs_idx) is out of bounds"))
    )

    # Select the subhalo potencial minimum
    return s_pos[:, subhalo_abs_idx]

end

"""
    computeCenter(data_dict::Dict, cm_type::Symbol)::Vector{<:Unitful.Length}

Compute a characteristic center of mass for the system.

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
  - `cm_type::Symbol`: It can be:

      + `:global_cm`  -> Center of mass of the whole system.
      + `:stellar_cm` -> Stellar center of mass.
      + `:zero`       -> Origin.

# Returns

  - The specified center of mass.
"""
function computeCenter(data_dict::Dict, cm_type::Symbol)::Vector{<:Unitful.Length}

    if cm_type == :global_cm

        return computeGlobalCenterOfMass(data_dict)

    elseif cm_type == :stellar_cm

        return computeCenterOfMass(data_dict[:stars]["POS "], data_dict[:stars]["MASS"])

    elseif cm_type == :zero

        return zeros(typeof(1.0u"kpc"), 3)

    end

    throw(ArgumentError("computeCenter: `cm_type` can only be :global_cm, :stellar_cm or :zero \
    but I got :$(center_type)"))

end

"""
    computeVcm(data_dict::Dict, subfind_idx::NTuple{2,Int})::Vector{<:Unitful.Velocity}

Read the velocity of the center of mass for a given halo or subhalo.

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
  - `subfind_idx::NTuple{2,Int}`: Tuple with two elements:

      + Index of the target halo (FoF group). Starts at 1.
      + Index of the target subhalo (subfind), relative the target halo. Starts at 1. If set to 0, the velocity of the whole halo with index `halo_idx` is returned.

# Returns

  - The specified velocity.
"""
function computeVcm(data_dict::Dict, subfind_idx::NTuple{2,Int})::Vector{<:Unitful.Velocity}

    # If there are no subfind data, return the origin
    if ismissing(data_dict[:gc_data].path) && !isSubfindActive(data_dict[:gc_data].path)
        return zeros(typeof(1.0u"kpc"), 3)
    end

    halo_idx, subhalo_rel_idx = subfind_idx

    # Load the necessary data
    g_n_subs = data_dict[:group]["G_Nsubs"]
    g_vel = data_dict[:group]["G_Vel"]
    s_vel = data_dict[:subhalo]["S_Vel"]

    # Check that the requested halo index is within bounds
    n_groups_total = data_dict[:gc_data].header.n_groups_total

    (
        !iszero(n_groups_total) && !any(isempty, [g_n_subs, g_vel, s_vel]) ||
        return zeros(typeof(1.0u"kpc"), 3)
    )

    (
        0 < halo_idx <= n_groups_total ||
        throw(ArgumentError("computeCenter: There is only $(n_groups_total) FoF goups in \
        $(data_dict[:gc_data].path), so halo_idx = $(halo_idx) is out of bounds"))
    )

    # Select the halo velocity if `subhalo_rel_idx` == 0
    isPositive(subhalo_rel_idx) || return g_vel[:, halo_idx]

    # Check that the requested subhalo index is within bounds
    n_subfinds = g_n_subs[halo_idx]
    (
        subhalo_rel_idx <= n_subfinds ||
        throw(ArgumentError("computeCenter: There is only $(n_subfinds) subhalos for the FoF \
        group $(halo_idx) in $(data_dict[:gc_data].path), so subhalo_rel_idx = \
        $(subhalo_rel_idx) is out of bounds"))
    )

    # Compute the number of subhalos and particles up to the last halo before `halo_idx`
    if isone(halo_idx)
        n_subs_floor = 0
    else
        n_subs_floor = sum(g_n_subs[1:(halo_idx - 1)]; init=0)
    end

    # Compute the subhalo absolute index
    subhalo_abs_idx = n_subs_floor + subhalo_rel_idx

    # Select the subhalo velocity
    return s_vel[:, subhalo_abs_idx]

end

"""
    computeVcm(data_dict::Dict, subhalo_abs_idx::Int)::Vector{<:Unitful.Velocity}

Read the velocity of the center of mass for a given subhalo.

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
  - `subhalo_abs_idx::Int`: Absolute index of the target subhalo (subfind). Starts at 1.

# Returns

  - The specified velocity.
"""
function computeVcm(data_dict::Dict, subhalo_abs_idx::Int)::Vector{<:Unitful.Velocity}

    # If there are no subfind data, return the origin
    if ismissing(data_dict[:gc_data].path) && !isSubfindActive(data_dict[:gc_data].path)
        return zeros(typeof(1.0u"kpc"), 3)
    end

    s_vel = data_dict[:subhalo]["S_Vel"]

    # Check that the requested subhalo index is within bounds
    n_subgroups_total = data_dict[:gc_data].header.n_subgroups_total

    !iszero(n_subgroups_total) && !isempty(s_vel) || return zeros(typeof(1.0u"kpc"), 3)

    (
        0 < subhalo_abs_idx <= n_subgroups_total ||
        throw(ArgumentError("computeCenter: There is only $(n_subgroups_total) subhalos in \
        $(data_dict[:gc_data].path), so subhalo_abs_idx = $(subhalo_abs_idx) is out of bounds"))
    )

    # Select the subhalo velocity
    return s_vel[:, subhalo_abs_idx]

end

"""
    computeVcm(data_dict::Dict, cm_type::Symbol)::Vector{<:Unitful.Velocity}

Compute the velocity of a characteristic center of mass for the system.

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
  - `cm_type::Symbol`: It can be:

      + `:global_cm`  -> Center of mass of the whole system.
      + `:stellar_cm` -> Stellar center of mass.

# Returns

  - The specified velocity.
"""
function computeVcm(data_dict::Dict, cm_type::Symbol)::Vector{<:Unitful.Velocity}

    if cm_type == :global_cm

        return computeGlobalVcm(data_dict)

    elseif cm_type == :stellar_cm

        return computeStellarVcm(data_dict)

    elseif cm_type == :zero

        return zeros(1.0u"km*s^-1", 3)

    end

    throw(ArgumentError("computeVcm: `cm_type` can only be :global_cm, :stellar_cm or :zero \
    but I got :$(center_type)"))

end

"""
    translatePoints(
        positions::Matrix{<:Number},
        new_origin::Vector{<:Number},
    )::Matrix{<:Number}

Translate a system of points, moving `new_origin` to the origin.

# Arguments

  - `positions::Matrix{<:Number}`: Points to be translated. Each column is a point and each row a dimension.
  - `new_origin::Vector{<:Number}`: Target origin.

# Returns

  - Matrix with the translated points.
"""
function translatePoints(
    positions::Matrix{<:Number},
    new_origin::Vector{<:Number},
)::Matrix{<:Number}

    !all(iszero, new_origin) || return positions

    return positions .- new_origin

end

"""
    translateData!(data_dict::Dict, translation::Union{Symbol,NTuple{2,Int},Int})::Nothing

Translate the positions of the cells/particles in `data_dict`.

!!! note

    The velocities will be boosted to the stellar center of mass of the system. If there are no stars, no transformation in applied to the velocities.

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
  - `translation::Union{Symbol,NTuple{2,Int},Int}=:zero`: Type of translation. The options are:

      + `:zero`                       -> No translation is applied.
      + `:global_cm`                  -> Sets the center of mass of the whole system as the new origin.
      + `:stellar_cm`                 -> Sets the stellar center of mass as the new origin.
      + `(halo_idx, subhalo_rel_idx)` -> Sets the position of the potencial minimum for the `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo) as the new origin.
      + `(halo_idx, 0)`               -> Sets the center of mass of the `halo_idx::Int` halo as the new origin.
      + `subhalo_abs_idx`             -> Sets the center of mass of the `subhalo_abs_idx::Int` as the new origin.
"""
function translateData!(data_dict::Dict, translation::Union{Symbol,NTuple{2,Int},Int})::Nothing

    translation != :zero || return nothing

    new_origin = computeCenter(data_dict, translation)
    stellar_vcm = computeVcm(data_dict, translation)

    @inbounds for type_symbol in snapshotTypes(data_dict)

        @inbounds for (block, values) in data_dict[type_symbol]

            if !isempty(values)
                @inbounds if block == "POS "
                    data_dict[type_symbol]["POS "] = translatePoints(values, new_origin)
                elseif block == "VEL "
                    data_dict[type_symbol]["VEL "] = translatePoints(values, stellar_vcm)
                end
            end

        end

    end

    return nothing

end

"""
    rotateData!(data_dict::Dict, axis_type::Symbol)::Nothing

Rotate the positions and velocities of the cells/particles in `data_dict`.

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
  - `rotation::Symbol`: Type of rotation. The options are:

      + `:zero`               -> No rotation is appplied.
      + `:global_am`          -> Sets the angular momentum of the whole system as the new z axis.
      + `:stellar_am`         -> Sets the stellar angular momentum as the new z axis.
      + `:stellar_pa`         -> Sets the stellar principal axis as the new coordinate system.
      + `:stellar_subhalo_pa` -> Sets the principal axis of the stars in the main subhalo as the new coordinate system.
"""
function rotateData!(data_dict::Dict, rotation::Symbol)::Nothing

    # Compute the rotation matrix
    if rotation == :zero

        return nothing

    elseif rotation == :global_am

        rotation_matrix = computeGlobalAMRotationMatrix(data_dict)

    elseif rotation == :stellar_am

        !isempty(data_dict[:stars]["MASS"]) || return nothing

            rotation_matrix = computeAMRotationMatrix(
                data_dict[:stars]["POS "],
                data_dict[:stars]["VEL "],
                data_dict[:stars]["MASS"],
            )

    elseif rotation == :stellar_pa

        !isempty(data_dict[:stars]["MASS"]) || return nothing

        rotation_matrix = computePARotationMatrix(
            data_dict[:stars]["POS "],
            data_dict[:stars]["VEL "],
            data_dict[:stars]["MASS"],
        )

    elseif rotation == :stellar_subhalo_pa

        star_data = deepcopy(data_dict)

        filterData!(
            star_data,
            filter_function=dd -> filterSubhalo(dd; halo_idx=1, subhalo_rel_idx=1),
        )

        !isempty(star_data[:stars]["MASS"]) || return nothing

        rotation_matrix = computePARotationMatrix(
            star_data[:stars]["POS "],
            star_data[:stars]["VEL "],
            star_data[:stars]["MASS"],
        )

    else

        throw(ArgumentError("rotateData!: I don't recognize the rotation :$(rotation)"))

    end

    @inbounds for type_symbol in snapshotTypes(data_dict)

        @inbounds for (block, values) in data_dict[type_symbol]

            @inbounds if block ∈ ["POS ", "VEL "] && !isempty(values)
                data_dict[type_symbol][block] = rotation_matrix * values
            end

        end

    end

    return nothing

end

"""
    rotateData!(data_dict::Dict, axis_type::NTuple{2,Int})::Nothing

Rotate the positions and velocities of the cells/particles in `data_dict`.

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
  - `rotation::NTuple{2,Int}`: Type of rotation. The options are:

      + `(halo_idx, subhalo_rel_idx)` -> Sets the principal axis of the stars in `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo) as the new coordinate system.
      + `(halo_idx, 0)`               -> Sets the principal axis of the stars in the `halo_idx::Int` halo as the new coordinate system.
"""
function rotateData!(data_dict::Dict, rotation::NTuple{2,Int})::Nothing

    # Compute the rotation matrix
    star_data = deepcopy(data_dict)

    filterData!(
        star_data,
        filter_function=dd -> filterSubhalo(dd; halo_idx=rotation[1], subhalo_rel_idx=rotation[2]),
    )

    !isempty(star_data[:stars]["MASS"]) || return nothing

    rotation_matrix = computePARotationMatrix(
        star_data[:stars]["POS "],
        star_data[:stars]["VEL "],
        star_data[:stars]["MASS"],
    )

    @inbounds for type_symbol in snapshotTypes(data_dict)

        @inbounds for (block, values) in data_dict[type_symbol]

            @inbounds if block ∈ ["POS ", "VEL "] && !isempty(values)
                data_dict[type_symbol][block] = rotation_matrix * values
            end

        end

    end

    return nothing

end

"""
    rotateData!(data_dict::Dict, axis_type::Int)::Nothing

Rotate the positions and velocities of the cells/particles in `data_dict`.

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
  - `rotation::Int`: Target subhalo absolute index, starting at 1. Sets the principal axis of the stars in the subhalo as the new coordinate system.
"""
function rotateData!(data_dict::Dict, rotation::Int)::Nothing

    # Compute the rotation matrix
    star_data = deepcopy(data_dict)

    filterData!(
        star_data,
        filter_function=dd -> filterSubhalo(dd, rotation),
    )

    !isempty(star_data[:stars]["MASS"]) || return nothing

    rotation_matrix = computePARotationMatrix(
        star_data[:stars]["POS "],
        star_data[:stars]["VEL "],
        star_data[:stars]["MASS"],
    )

    @inbounds for type_symbol in snapshotTypes(data_dict)

        @inbounds for (block, values) in data_dict[type_symbol]

            @inbounds if block ∈ ["POS ", "VEL "] && !isempty(values)
                data_dict[type_symbol][block] = rotation_matrix * values
            end

        end

    end

    return nothing

end

"""
    function plotParams(quantity::Symbol)::PlotParams

Select the plotting parameters for a given `quantity`.

# Arguments

  - `quantity::Symbol`: The options are:

      + `:stellar_mass`               -> Stellar mass.
      + `:gas_mass`                   -> Gas mass.
      + `:dm_mass`                    -> Dark matter mass.
      + `:bh_mass`                    -> Black hole mass.
      + `:generic_mass`               -> Parameters for plots with several diferent masses.
      + `:stellar_number`             -> Number of stellar particles.
      + `:gas_number`                 -> Number of gas cells.
      + `:dm_number`                  -> Number of dark matter particles.
      + `:bh_number`                  -> Number of black hole particles.
      + `:molecular_mass`             -> Molecular hydrogen (``\\mathrm{H_2}``) mass.
      + `:atomic_mass`                -> Atomic hydrogen (``\\mathrm{HI}``) mass.
      + `:ionized_mass`               -> Ionized hydrogen (``\\mathrm{HII}``) mass.
      + `:neutral_mass`               -> Neutral hydrogen (``\\mathrm{HI + H_2}``) mass.
      + `:molecular_fraction`         -> Gas mass fraction of molecular hydrogen.
      + `:atomic_fraction`            -> Gas mass fraction of atomic hydrogen.
      + `:ionized_fraction`           -> Gas mass fraction of ionized hydrogen.
      + `:neutral_fraction`           -> Gas mass fraction of neutral hydrogen.
      + `:molecular_neutral_fraction` -> Fraction of molecular hydrogen in the neutral gas.
      + `:generic_fraction`           -> Parameters for plots with several diferent fraction.
      + `:gas_mass_density`           -> Gas mass density.
      + `:gas_number_density`         -> Gas number density.
      + `:molecular_number_density`   -> Molecular hydrogen number density.
      + `:atomic_number_density`      -> Atomic hydrogen number density.
      + `:ionized_number_density`     -> Ionized hydrogen number density.
      + `:neutral_number_density`     -> Neutral hydrogen number density.
      + `:stellar_area_density`       -> Stellar area mass density, for a radius of `FILTER_R`.
      + `:gas_area_density`           -> Gas area mass density, for a radius of `FILTER_R`.
      + `:molecular_area_density`     -> Molecular hydrogen area mass density, for a radius of `FILTER_R`.
      + `:atomic_area_density`        -> Atomic hydrogen area mass density, for a radius of `FILTER_R`.
      + `:ionized_area_density`       -> Ionized hydrogen area mass density, for a radius of `FILTER_R`.
      + `:neutral_area_density`       -> Neutral hydrogen area mass density, for a radius of `FILTER_R`.
      + `:sfr_area_density`           -> Star formation rate area density, for the last `AGE_RESOLUTION_ρ` and a radius of `FILTER_R`.
      + `:generic_area_density`       -> Parameters for plots with several diferent area densities.
      + `:gas_metallicity`            -> Mass fraction of all elements above He in the gas (solar units).
      + `:stellar_metallicity`        -> Mass fraction of all elements above He in the stars (solar units).
      + `:X_gas_abundance`            -> Gas abundance of element ``\\mathrm{X}``, as ``12 + \\log_{10}(\\mathrm{X \\, / \\, H})``. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
      + `:X_stellar_abundance`        -> Stellar abundance of element ``\\mathrm{X}``, as ``12 + \\log_{10}(\\mathrm{X \\, / \\, H})``. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
      + `:stellar_radial_distance`    -> Distance of every stellar particle to the origin.
      + `:gas_radial_distance`        -> Distance of every gas cell to the origin.
      + `:dm_radial_distance`         -> Distance of every dark matter particle to the origin.
      + `:stellar_xy_distance`        -> Projected distance of every stellar particle to the origin.
      + `:gas_xy_distance`            -> Projected distance of every gas cell to the origin.
      + `:dm_xy_distance`             -> Projected distance of every dark matter particle to the origin.
      + `:stellar_specific_am`        -> Norm of the stellar specific angular momentum.
      + `:gas_specific_am`            -> Norm of the gas specific angular momentum.
      + `:dm_specific_am`             -> Norm of the dark matter specific angular momentum.
      + `:stellar_circularity`        -> Stellar circularity.
      + `:stellar_vcirc`              -> Stellar circular velocity.
      + `:stellar_vradial`            -> Stellar radial speed.
      + `:stellar_vtangential`        -> Stellar tangential speed.
      + `:stellar_vzstar`             -> Stellar speed in the z direction, computed as ``v_z \\, \\sign(z)``.
      + `:stellar_age`                -> Stellar age.
      + `:sfr`                        -> The star formation rate.
      + `:ssfr`                       -> The specific star formation rate.
      + `:observational_sfr`          -> The star formation rate of the last `AGE_RESOLUTION`.
      + `:observational_ssfr`         -> The specific star formation rate of the last `AGE_RESOLUTION`.
      + `:temperature`                -> Gas temperature, as ``\\log_{10}(T \\, / \\, \\mathrm{K})``.
      + `:scale_factor`               -> Scale factor.
      + `:redshift`                   -> Redshift.
      + `:physical_time`              -> Physical time since the Big Bang.
      + `:lookback_time`              -> Physical time left to reach the last snapshot.

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
            request    = Dict(:stars => ["MASS", "POS ", "SOFT"]),
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

    elseif quantity == :dm_mass

        plot_params = PlotParams(;
            request    = Dict(:halo => ["MASS", "POS ", "SOFT"]),
            var_name   = L"M_\mathrm{DM}",
            exp_factor = 10,
            unit       = u"Msun",
        )

    elseif quantity == :bh_mass

        plot_params = PlotParams(;
            request  = Dict(:black_hole => ["MASS", "POS ", "SOFT"]),
            var_name = L"M_\mathrm{BH}",
            unit     = u"Msun",
        )

    elseif quantity == :generic_mass

        plot_params = PlotParams(;
            request  = Dict(
                :stars => ["MASS", "POS ", "SOFT"],
                :gas => ["MASS", "POS ", "FRAC", "NH  ", "NHP ", "PRES", "RHO "],
                :dm_mass => ["MASS", "POS ", "SOFT"],
                :bh_mass => ["MASS", "POS ", "SOFT"],
            ),
            var_name = L"M",
            exp_factor = 10,
            unit     = u"Msun",
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

    elseif quantity == :molecular_mass

        plot_params = PlotParams(;
            request    = Dict(:gas => ["MASS", "POS ", "FRAC", "NH  ", "NHP ", "PRES", "RHO "]),
            var_name   = L"M_\mathrm{H2}",
            exp_factor = 10,
            unit       = u"Msun",
        )

    elseif quantity == :atomic_mass

        plot_params = PlotParams(;
            request    = Dict(:gas => ["MASS", "POS ", "FRAC", "NH  ", "NHP ", "PRES", "RHO "]),
            var_name   = L"M_\mathrm{HI}",
            exp_factor = 10,
            unit       = u"Msun",
        )

    elseif quantity == :ionized_mass

        plot_params = PlotParams(;
            request    = Dict(:gas => ["MASS", "POS ", "FRAC", "NH  ", "NHP ", "RHO "]),
            var_name   = L"M_\mathrm{HII}",
            exp_factor = 10,
            unit       = u"Msun",
        )

    elseif quantity == :neutral_mass

        plot_params = PlotParams(;
            request    = Dict(:gas => ["MASS", "POS ", "FRAC", "NH  ", "NHP ", "RHO "]),
            var_name   = L"M_\mathrm{H2 + HI}",
            exp_factor = 10,
            unit       = u"Msun",
        )

    elseif quantity == :molecular_fraction

        plot_params = PlotParams(;
            request  = Dict(:gas => ["MASS", "POS ", "FRAC", "NH  ", "NHP ", "PRES"]),
            var_name = L"f_\mathrm{H2}",
        )

    elseif quantity == :atomic_fraction

        plot_params = PlotParams(;
            request  = Dict(:gas => ["MASS", "POS ", "FRAC", "NH  ", "NHP ", "PRES"]),
            var_name = L"f_\mathrm{HI}",
        )

    elseif quantity == :ionized_fraction


        plot_params = PlotParams(;
            request  = Dict(:gas => ["MASS", "POS ", "FRAC", "NH  ", "NHP "]),
            var_name = L"f_\mathrm{HII}",
        )

    elseif quantity == :neutral_fraction

        plot_params = PlotParams(;
            request  = Dict(:gas => ["MASS", "POS ", "FRAC", "NH  ", "NHP "]),
            var_name = L"f_\mathrm{H2 + H_I}",
        )

    elseif quantity == :molecular_neutral_fraction

        plot_params = PlotParams(;
            request  = Dict(:gas => ["MASS", "POS ", "RHO ", "FRAC", "NH  ", "NHP ", "PRES"]),
            var_name = L"f_\mathrm{H2}^\star",
        )

    elseif quantity == :generic_fraction

        plot_params = PlotParams(;
            request  = Dict(:gas => ["MASS", "POS ", "RHO ", "FRAC", "NH  ", "NHP ", "PRES"]),
            var_name = L"f",
        )

    elseif quantity == :gas_mass_density

        plot_params = PlotParams(;
            request  = Dict(:gas => ["MASS", "POS ", "RHO "]),
            var_name = L"\rho_\mathrm{gas}",
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
            request  = Dict(:gas => ["MASS", "POS ", "RHO ", "FRAC", "NH  ", "NHP ", "PRES"]),
            var_name = L"n_\mathrm{H2}",
            unit     = u"cm^-3",
        )

    elseif quantity == :atomic_number_density

        plot_params = PlotParams(;
            request  = Dict(:gas => ["MASS", "POS ", "RHO ", "FRAC", "NH  ", "NHP ", "PRES"]),
            var_name = L"n_\mathrm{HI}",
            unit     = u"cm^-3",
        )

    elseif quantity == :ionized_number_density

        plot_params = PlotParams(;
            request  = Dict(:gas => ["MASS", "POS ", "RHO ", "FRAC", "NH  ", "NHP "]),
            var_name = L"n_\mathrm{HII}",
            unit     = u"cm^-3",
        )

    elseif quantity == :neutral_number_density

        plot_params = PlotParams(;
            request  = Dict(:gas => ["MASS", "POS ", "RHO ", "FRAC", "NH  ", "NHP "]),
            var_name = L"n_\mathrm{H2 + HI}",
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
            request  = Dict(:gas => ["MASS", "POS "]),
            var_name = L"\Sigma_\mathrm{gas}",
            unit     = u"Msun*pc^-2",
        )

    elseif quantity == :molecular_area_density

        plot_params = PlotParams(;
            request  = Dict(:gas => ["MASS", "POS ", "FRAC", "NH  ", "NHP ", "PRES"]),
            var_name = L"\Sigma_\mathrm{H2}",
            unit     = u"Msun*pc^-2",
        )

    elseif quantity == :atomic_area_density

        plot_params = PlotParams(;
            request  = Dict(:gas => ["MASS", "POS ", "FRAC", "NH  ", "NHP ", "PRES"]),
            var_name = L"\Sigma_\mathrm{HI}",
            unit     = u"Msun*pc^-2",
        )

    elseif quantity == :ionized_area_density

        plot_params = PlotParams(;
            request  = Dict(:gas => ["MASS", "POS ", "FRAC", "NH  ", "NHP "]),
            var_name = L"\Sigma_\mathrm{HII}",
            unit     = u"Msun*pc^-2",
        )

    elseif quantity == :neutral_area_density

        plot_params = PlotParams(;
            request  = Dict(:gas => ["MASS", "POS ", "FRAC", "NH  ", "NHP "]),
            var_name = L"\Sigma_\mathrm{H2 + HI}",
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
                :gas => ["MASS", "POS ", "FRAC", "NH  ", "NHP ", "PRES"],
            ),
            var_name = L"\Sigma",
            unit     = u"Msun*pc^-2",
        )

    elseif quantity == :gas_metallicity

        plot_params = PlotParams(;
            request  = Dict(:gas => ["MASS", "POS ", "GZ  "]),
            var_name = L"Z_\mathrm{gas} \, / \, Z_\odot",
        )

    elseif quantity == :stellar_metallicity

        plot_params = PlotParams(;
            request  = Dict(:stars => ["MASS", "POS ", "GZ2 "]),
            var_name = L"Z_\star \, / \, Z_\odot",
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

        plot_params = PlotParams(;
            request  = Dict(
                :gas        => ["MASS", "POS "],
                :halo       => ["MASS", "POS "],
                :stars      => ["MASS", "POS ", "VEL "],
                :black_hole => ["MASS", "POS "],
            ),
            var_name = L"\epsilon",
        )

    elseif quantity == :stellar_vcirc

        plot_params = PlotParams(;
            request  = Dict(
                :gas        => ["MASS", "POS "],
                :halo       => ["MASS", "POS "],
                :stars      => ["MASS", "POS "],
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
            unit     = u"Myr",
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

####################################################################################################
# Derived quantities.
####################################################################################################

"""
    findQtyExtrema(
        simulation_path::String,
        slice_n::Int,
        type_symbol::Symbol,
        block::String;
        <keyword arguments>
    )::NTuple{2,<:Number}

Compute the minimum and maximum values of `block`.

# Arguments

  - `simulation_path::String`: Path to the simulation directory, set in the code variable `OutputDir`.
  - `slice_n::Int`: Selects which snapshot to plot, starts at 1 and is independent of the number in the file name. If every snapshot is present, `slice_n` = filename_number + 1. If set to a negative number, the values in the whole simulation will be compared.
  - `type_symbol::Symbol`: Cell/particle type. The possibilities are the keys of [`PARTICLE_INDEX`](@ref).
  - `block::String`: Target block. The possibilities are the keys of [`QUANTITIES`](@ref).
  - `f::Function=identity`: A functions with the signature:

    `f(data) -> values`

    where

      + `data::VecOrMat{<:Number}`: Data returned by [`getBlock`](@ref).
      + `values::Vector{<:Number}`: A vector with the values to be compared.
  - `warnings::Bool=true`: If a warning will be given when there is missing data.

# Returns

  - Tuple with the minimum and maximum values.
"""
function findQtyExtrema(
    simulation_path::String,
    slice_n::Int,
    type_symbol::Symbol,
    block::String;
    f::Function=identity,
    warnings::Bool=true,
)::NTuple{2,<:Number}

    (
        isdir(simulation_path) ||
        throw(ArgumentError("findQtyExtrema: $(simulation_path) does not exist as a directory"))
    )

    simulation_table = makeSimulationTable(simulation_path; warnings)

    if slice_n > 0

        # Get the number in the filename
        snap_n = safeSelect(simulation_table[!, :numbers], slice_n; warnings)

        # Check that after slicing there is one snapshot left
        (
            !isempty(snap_n) ||
            throw(ArgumentError("findQtyExtrema: There are no snapshots with `slice_n` = \
            $(slice_n), the contents of $(simulation_path) are: \n$(simulation_table)"))
        )

        # Find the target row and snapshot path
        snapshot_row = filter(:numbers => ==(lpad(snap_n, 3, "0")), simulation_table)
        snapshot_path = snapshot_row[1, :snapshot_paths]

        (
            !ismissing(snapshot_path) ||
            throw(ArgumentError("findQtyExtrema: The snapshot number $(slice_n) seems \
            to be missing"))

        )

        values = f(getBlock(snapshot_path, type_symbol, block))

        return extrema(values)

    end

    snapshot_paths = filter!(!ismissing, snapshot_row[!, :snapshot_paths])

    (
        !isempty(snapshot_paths) ||
        throw(ArgumentError("findQtyExtrema: I could not find any snapshots in $(simulation_path)"))
    )

    values = [f(getBlock(snapshot_path, type_symbol, block)) for snapshot_path in snapshot_paths]

    return extrema(Iterators.flatten(values))

end

@doc raw"""
    energyIntegrand(a::Real, header::SnapshotHeader)::Float64

The integrand of the integral that converts the scale factor into physical time:

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
  - `header::SnapshotHeader`: Header of the relevant snapshot file.

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
    H = header.h0 * HUBBLE_CONSTANT * a

    # Return the integrand, in Gyr
    return 1.0 / (H * sqrt(E))

end

"""
    computeProfile(
        positions::Matrix{<:Unitful.Length},
        quantity::Vector{<:Number},
        grid::CircularGrid;
        <keyword arguments>
    )::Vector{<:Number}

Compute a profile.

# Arguments

  - `positions::Matrix{<:Unitful.Length}`: Positions of the cells/particles. Each column is a cell/particle and each row a dimension.
  - `quantity::Vector{<:Number}`: The profile will be of this quantity.
  - `grid::CircularGrid`: Circular grid.
  - `norm_values::Vector{<:Number}=Number[]`: Values to normalize `quantity`.
  - `flat::Bool=true`: If the profile will be 2D, using rings, or 3D, using spherical shells.
  - `total::Bool=false`: If the sum (default) or the mean of `quantity` will be computed for each bin.
  - `cumulative::Bool=false`: If the profile will be accumulated or not.
  - `density::Bool=false`: If the profile will be of the density of `quantity`.

# Returns

  - Vector with the values of the profile.
"""
function computeProfile(
    positions::Matrix{<:Unitful.Length},
    quantity::Vector{<:Number},
    grid::CircularGrid;
    norm_values::Vector{<:Number}=Number[],
    flat::Bool=true,
    total::Bool=false,
    cumulative::Bool=false,
    density::Bool=false,
)::Vector{<:Number}

    # Return a null profile if `quantity` is empty
    !isempty(quantity) || return fill(NaN, length(grid.grid))

    # Compute the distances of the cells/particles to the center of the grid
    if flat
        distances = computeDistance(positions[1:2, :]; center=grid.center[1:2])
    else
        distances = computeDistance(positions; center=grid.center)
    end

    # Compute the histogram of `quantity`
    if isempty(norm_values)

        profile = histogram1D(distances, quantity, grid; total)

    else

        quantity_histogram = histogram1D(distances, quantity, grid; total)
        norm_values_histogram = histogram1D(distances, norm_values, grid; total)

        replace!(x -> iszero(x) ? oneunit(x) : x, norm_values_histogram)

        profile = quantity_histogram ./ norm_values_histogram

    end

    region = flat ? grid.bin_areas : grid.bin_volumes

    if cumulative
        return density ? cumsum(profile) ./ cumsum(region) : cumsum(profile)
    end

    return density ? profile ./ region : profile

end

"""
    findRealStars(path::String)::Vector{Int}

Find the indices of the stars in a snapshot, excluding wind particles.

# Arguments

  - `path::String`: Path to the snapshot file or folder.

# Returns

  - A vector with the indices of the stars.
"""
function findRealStars(path::String)::Vector{Bool}

    if isfile(path)

        (
            HDF5.ishdf5(path) ||
            throw(ArgumentError("findRealStars: The file $(path) is not in the \
            HDF5 format, I don't know how to read it"))
        )

        time_of_birth = h5open(path, "r") do snapshot
            if PARTICLE_CODE_NAME[:stars] ∉ keys(snapshot)
                Float64[]
            else
                read(snapshot[PARTICLE_CODE_NAME[:stars]], QUANTITIES["GAGE"].hdf5_name)
            end
        end

        return isempty(time_of_birth) ? Bool[] : map(isPositive, time_of_birth)

    elseif isdir(path)

        sub_files = glob("$(SNAP_BASENAME)_*.*.hdf5", path)

        (
            !isempty(sub_files) && all(HDF5.ishdf5, sub_files) ||
            throw(ArgumentError("findRealStars: The directory $(path) does not contain \
            snapshot sub-files in the HDF5 format"))
        )

        return vcat([findRealStars(sub_file) for sub_file in sub_files]...)

    else

        throw(ArgumentError("findRealStars: $(path) does not exist as a file or folder"))

    end

end

"""
    countStars(path::String)::Int

Count the number of stars in a snapshot, excluding wind particles.

# Arguments

  - `path::String`: Path to the snapshot file or folder.

# Returns

  - The number of stars.
"""
countStars(path::String)::Int = count(findRealStars(path))

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

        # For cosmological simulations, the time in the snapshot is the scale factor
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
        metals::Matrix{Float32},
        internal_energy::Vector{<:SpecificEnergy},
        electron_fraction::Vector{Float32},
    )::Vector{<:Unitful.Temperature}

Compute the gas temperature of a group of gas cells.

# Arguments

  - `metals::Matrix{Float32}`: Matrix with the mass content for every tracked element. Each row is an element, and each column a cell.
  - `internal_energy::Vector{<:SpecificEnergy}`: Specific internal energy of every gas cell.
  - `electron_fraction::Vector{Float32}`: Number fraction of electrons in every gas cell.

# Returns

  - The temperature of each gas cell.
"""
function computeTemperature(
    metals::Matrix{Float32},
    internal_energy::Vector{<:SpecificEnergy},
    electron_fraction::Vector{Float32},
)::Vector{<:Unitful.Temperature}

    # xH := mass_fraction_of_hydrogen
    xH = metals[ELEMENT_INDEX[:H], :]

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
    computeDistance(
        positions::Matrix{<:Number};
        <keyword arguments>
    )::Vector{<:Number}

Compute the distance of a group of points to `center`.

# Arguments

  - `positions::Matrix{<:Number}`: Positions of the points. Each column is a point and each row a dimension.
  - `center::Union{Vector{<:Number},Nothing}=nothing`: Origin used to compute the distances. If set to `nothing`, 0 is used.

# Returns

  - The distance of every point to `center`.
"""
function computeDistance(
    positions::Matrix{<:Number};
    center::Union{Vector{<:Number},Nothing}=nothing,
)::Vector{<:Number}

    if isempty(positions)
        return eltype(positions)[]
    end

    if center === nothing
        return [norm(col) for col in eachcol(positions)]
    end

    (
        length(center) == size(positions, 1) ||
        throw(ArgumentError("computeDistance: `center` must have as many elements as `positions` \
        has rows, but I got length(center) = $(length(center)) and size(positions, 1) = \
        $(size(positions, 1))"))
    )

    return [norm(col .- center) for col in eachcol(positions)]

end

"""
    computeCenterOfMass(
        positions::Matrix{<:Unitful.Length},
        mass::Vector{<:Unitful.Mass},
    )::Vector{<:Unitful.Length}

Compute the center of mass of a group of cells/particles.

# Arguments

  - `positions::Matrix{<:Unitful.Length}`: Positions of the cells/particles. Each column is a cell/particle and each row a dimension.
  - `masses::Vector{<:Unitful.Mass}`: Masses of the cells/particles.

# Returns

  - The center of mass.
"""
function computeCenterOfMass(
    positions::Matrix{<:Unitful.Length},
    masses::Vector{<:Unitful.Mass},
)::Vector{<:Unitful.Length}

    # Check for missing data
    !any(isempty, [positions, masses]) || return zeros(typeof(1.0u"kpc"), 3)

    # Compute the total mass
    M = sum(masses)

    # Compute the center of mass
    center_of_mass = [sum(row .* masses) / M for row in eachrow(positions)]

    return center_of_mass

end

"""
    computeGlobalCenterOfMass(data_dict::Dict)::Vector{<:Unitful.Length}

Compute the center of mass of the whole system in `data`.

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

  - The center of mass.
"""
function computeGlobalCenterOfMass(data_dict::Dict)::Vector{<:Unitful.Length}

    type_symbols = snapshotTypes(data_dict)

    filter!(ts -> !isempty(data_dict[ts]["POS "]), type_symbols)

    # Concatenate the position and masses of all the cells and particles in the system
    positions = hcat([data_dict[type_symbol]["POS "] for type_symbol in type_symbols]...)
    masses    = vcat([data_dict[type_symbol]["MASS"] for type_symbol in type_symbols]...)

    # Check for missing data
    !any(isempty, [positions, masses]) || return zeros(typeof(1.0u"kpc"), 3)

    @debug("computeGlobalCenterOfMass: The center of mass will be computed using $(type_symbols)")

    return computeCenterOfMass(positions, masses)

end

"""
    computeInertiaTensor(
        positions::Matrix{<:Unitful.Length},
        masses::Vector{<:Unitful.Mass},
    )::Matrix{Float64}

Compute the inertia tensor of a group of cells/particles.

# Arguments

  - `positions::Matrix{<:Unitful.Length}`: Positions of the cells/particles. Each column is a cell/particle and each row a dimension.
  - `masses::Vector{<:Unitful.Mass}`: Masses of the cells/particles.

# Returns

  - The inertia tensor.
"""
function computeInertiaTensor(
    positions::Matrix{<:Unitful.Length},
    masses::Vector{<:Unitful.Mass},
)::Matrix{Float64}

    # Check for missing data
    (
        !any(isempty, [positions, masses])  ||
        throw(ArgumentError("computeInertiaTensor: The inertia tensor is not defined for an \
        empty system"))
    )

    # Allocate memory
    J = zeros(Float64, (3, 3))

    # Compute the inertia tensor
    for (r, m) in zip(eachcol(positions), masses)

        x, y, z = r

        J[1, 1] += m * (y * y + z * z)
        J[2, 2] += m * (x * x + z * z)
        J[3, 3] += m * (x * x + y * y)

        J[1, 2] -= m * x * y
        J[1, 3] -= m * x * z
        J[2, 3] -= m * y * z

    end

    J[2, 1] = J[1, 2]
    J[3, 1] = J[1, 3]
    J[3, 2] = J[2, 3]

    return J

end

"""
    computeAMRotationMatrix(
        positions::Matrix{<:Unitful.Length},
        velocities::Matrix{<:Unitful.Velocity},
        masses::Vector{<:Unitful.Mass},
    )::Union{Matrix{Float64},UniformScaling{Bool}}

Compute the rotation matrix that will turn the total angular momentum into the z axis; when view as an active (alibi) trasformation.

# Arguments

  - `positions::Matrix{<:Unitful.Length}`: Positions of the cells/particles. Each column is a cell/particle and each row a dimension.
  - `velocities::Matrix{<:Unitful.Velocity}`: Velocities of the cells/particles. Each column is a cell/particle and each row a dimension.
  - `masses::Vector{<:Unitful.Mass}`: Mass of every cell/particle.

# Returns

  - The rotation matrix.
"""
function computeAMRotationMatrix(
    positions::Matrix{<:Unitful.Length},
    velocities::Matrix{<:Unitful.Velocity},
    masses::Vector{<:Unitful.Mass},
)::Union{Matrix{Float64},UniformScaling{Bool}}

    # Check for missing data
    !any(isempty, [positions, velocities, masses]) || return I

    # Compute the total angular momentum
    L = computeTotalAngularMomentum(positions, velocities, masses)

    # Rotation vector
    n = [L[2], -L[1], 0.0]

    # Angle of rotation
    θ = acos(L[3])

    return Matrix{Float64}(AngleAxis(θ, n...))

end

"""
    computePARotationMatrix(
        positions::Matrix{<:Unitful.Length},
        velocities::Matrix{<:Unitful.Velocity},
        masses::Vector{<:Unitful.Mass},
    )::Union{Matrix{Float64},UniformScaling{Bool}}

Compute the rotation matrix that will turn the pricipal axis into the new coordinate system; when view as an pasive (alias) trasformation.

# Arguments

  - `positions::Matrix{<:Unitful.Length}`: Positions of the cells/particles. Each column is a cell/particle and each row a dimension.
  - `velocities::Matrix{<:Unitful.Velocity}`: Velocities of the cells/particles. Each column is a cell/particle and each row a dimension.
  - `masses::Vector{<:Unitful.Mass}`: Mass of every cell/particle.

# Returns

  - The rotation matrix.
"""
function computePARotationMatrix(
    positions::Matrix{<:Unitful.Length},
    velocities::Matrix{<:Unitful.Velocity},
    masses::Vector{<:Unitful.Mass},
)::Union{Matrix{Float64},UniformScaling{Bool}}

    # Check for missing data
    !isempty(positions) || return I

    # Principal axis operator
    R = ustrip.(positions * positions')

    # Reverse the order of the eigenvectors, making the last column the eigenvector
    # with the largest eigenvalue, which should correspond to the new z axis
    pa = eigvecs(R)[:, end:-1:1]

    # Compute the total angular momentum
    L = computeTotalAngularMomentum(positions, velocities, masses; normal=true)

    # 3rd principal axis ≡ new z axis
    pa_z = pa[:, 3]

    # Rotate the principal axis as to align the thid component with the angular momentum
    θ = acos(L ⋅ pa_z)
    n = cross(pa_z, L)
    aligned_pa = AngleAxis(θ, n...) * pa

    # The rotation matrix is made from the principal axis as rows
    rotation_matrix = aligned_pa'

    if det(rotation_matrix) < 0.0
        # If the determinant is < 0, that means that the chosen principal axis for the x and y
        # directions form a left-handed cartesian reference system (x × y = -z). When applying
        # this as a rotation, the z axis will be flipped. So, in this case we swap the x and y
        # axis to get a right-handed cartesian reference system (x × y = z) and generate the
        # correct rotation
        rotation_matrix[1, :], rotation_matrix[2, :] = rotation_matrix[2, :], rotation_matrix[1, :]
    end

    return Matrix{Float64}(rotation_matrix)

end

"""
    computeGlobalAMRotationMatrix(data_dict::Dict)::Union{Matrix{Float64},UniformScaling{Bool}}

Compute the rotation matrix that will turn the total angular momentum, of the whole system in `data`, into the z axis; when view as an active (alibi) trasformation.

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

  - The rotation matrix.
"""
function computeGlobalAMRotationMatrix(data_dict::Dict)::Union{Matrix{Float64},UniformScaling{Bool}}

    type_symbols = snapshotTypes(data_dict)

    filter!(ts -> !isempty(data_dict[ts]["POS "]), type_symbols)

    # Concatenate the positions, velocities, and masses of all the cells and particles in the system
    positions  = hcat([data_dict[type_symbol]["POS "] for type_symbol in type_symbols]...)
    velocities = hcat([data_dict[type_symbol]["VEL "] for type_symbol in type_symbols]...)
    masses     = vcat([data_dict[type_symbol]["MASS"] for type_symbol in type_symbols]...)

    # Check for missing data
    !any(isempty, [positions, velocities, masses]) || return I

    @debug("computeGlobalAMRotationMatrix: The rotation matrix will be computed \
    using $(type_symbols)")

    return computeAMRotationMatrix(positions, velocities, masses)

end

"""
    computeAngularMomentum(
        positions::Matrix{<:Unitful.Length},
        velocities::Matrix{<:Unitful.Velocity},
        masses::Vector{<:Unitful.Mass},
    )::Vector{Vector{<:AngularMomentum}}

Compute the angular momentum of each cell/particle, with respect to the origin.

# Arguments

  - `positions::Matrix{<:Unitful.Length}`: Positions of the cells/particles. Each column is a cell/particle and each row a dimension.
  - `velocities::Matrix{<:Unitful.Velocity}`: Velocities of the cells/particles. Each column is a cell/particle and each row a dimension.
  - `masses::Vector{<:Unitful.Mass}`: Mass of every cell/particle.

# Returns

  - The angular momentum of each cell/particle.
"""
function computeAngularMomentum(
    positions::Matrix{<:Unitful.Length},
    velocities::Matrix{<:Unitful.Velocity},
    masses::Vector{<:Unitful.Mass},
)::Vector{Vector{<:AngularMomentum}}

    # Check for missing data
    !any(isempty, [positions, velocities, masses]) || return [0.0, 0.0, 1.0]

    iterator = zip(masses, eachcol(positions), eachcol(velocities))

    return map(x -> x[1] .* cross(x[2], x[3]), iterator)

end

"""
    computeTotalAngularMomentum(
        positions::Matrix{<:Unitful.Length},
        velocities::Matrix{<:Unitful.Velocity},
        masses::Vector{<:Unitful.Mass};
        <keyword arguments>
    )::Vector{<:Number}

Compute the total angular momentum of a group of cells/particles, with respect to the origin.

# Arguments

  - `positions::Matrix{<:Unitful.Length}`: Positions of the cells/particles. Each column is a cell/particle and each row a dimension.
  - `velocities::Matrix{<:Unitful.Velocity}`: Velocities of the cells/particles. Each column is a cell/particle and each row a dimension.
  - `masses::Vector{<:Unitful.Mass}`: Mass of every cell/particle.
  - `normal::Bool=true`: If the result will be normalized.

# Returns

  - The angular momentum.
"""
function computeTotalAngularMomentum(
    positions::Matrix{<:Unitful.Length},
    velocities::Matrix{<:Unitful.Velocity},
    masses::Vector{<:Unitful.Mass};
    normal::Bool=true,
)::Vector{<:Number}

    # Check for missing data
    !any(isempty, [positions, velocities, masses]) || return [0.0, 0.0, 1.0]

    iterator = zip(masses, eachcol(positions), eachcol(velocities))

    # Compute the total angular momentum
    L = mapreduce(x -> x[1] .* cross(x[2], x[3]), +, iterator)

    return normal ? normalize!(ustrip.(L)) : L

end

"""
    computeGlobalAngularMomentum(data_dict::Dict; <keyword arguments>)::Vector{<:Number}

Compute the total angular momentum with respect to the origin of the whole system in `data`.

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
  - `normal::Bool=true`: If the result will be normalized.

# Returns

  - The angular momentum.
"""
function computeGlobalAngularMomentum(data_dict::Dict; normal::Bool=true)::Vector{<:Number}

    type_symbols = snapshotTypes(data_dict)

    filter!(ts -> !isempty(data_dict[ts]["POS "]), type_symbols)

    # Concatenate the position, velocities, and masses of all the cells and particles in the system
    positions  = hcat([data_dict[type_symbol]["POS "] for type_symbol in type_symbols]...)
    velocities = hcat([data_dict[type_symbol]["VEL "] for type_symbol in type_symbols]...)
    masses     = vcat([data_dict[type_symbol]["MASS"] for type_symbol in type_symbols]...)

    # Check for missing data
    !any(isempty, [positions, velocities, masses]) || return [0.0, 0.0, 1.0]

    @debug("computeGlobalAngularMomentum: The angular momentum will be computed \
    using $(type_symbols)")

    return computeTotalAngularMomentum(positions, velocities, masses; normal)

end

@doc raw"""
    computeSpinParameter(
        positions::Matrix{<:Unitful.Length},
        masses::Vector{<:Unitful.Mass},
        velocities::Matrix{<:Unitful.Velocity};
        <keyword arguments>
    )::Float64

Compute the spin parameter for a system of cells/particles, with respect to the origin.

The spin parameter was originally defined by Peebles (1969) as,

```math
\lambda = \frac{J \, \sqrt{E}}{G \, M^{5/2}} \, ,
```

where $J$ is the norm of the total angular momentum, $M$ the total mass, $G$ the gravitational constant, and

```math
E = |E_P + E_k| \, ,
```

where $E_P$ is the total potencial energy and $E_k$ is the total kinetic energy (including thermal energy of the gas).

Due to the computational complexity of calculating $E_P$ for a large group of particles, Bullock et al. (2001) proposed an alternative definition of the spin parameter,

```math
\lambda = \frac{J}{\sqrt{2} \, M \, R \, V} \, ,
```

where $J$ is the norm of the total angular momentum inside a sphere of radius $R$ containing mass $M$, and

```math
V = \sqrt{\frac{G \, M}{R}} \, ,
```

is the circular velocity.

# Arguments

  - `positions::Matrix{<:Unitful.Length}`: Positions of the cells/particles. Each column is a cell/particle and each row a dimension.
  - `velocities::Matrix{<:Unitful.Velocity}`: Velocities of the cells/particles. Each column is a cell/particle and each row a dimension.
  - `masses::Vector{<:Unitful.Mass}`: Mass of every cell/particle.
  - `R::Unitful.Length=FILTER_R`: Radius.

# Returns

  - The spin parameter.

# References

P. J. E. Peebles (1969). *Origin of the Angular Momentum of Galaxies*. Astrophysical Journal, **155**, 393. [doi:10.1086/149876](https://doi.org/10.1086/149876)

J. S. Bullock et al. (2001). *A Universal Angular Momentum Profile for Galactic Halos*. The Astrophysical Journal, **555(1)**, 240. [doi:10.1086/321477](https://doi.org/10.1086/321477)

J. Zjupa et al. (2017). *Angular momentum properties of haloes and their baryon content in the Illustris simulation*. Monthly Notices of the Royal Astronomical Society, **466(2)**, 1625–1647. [doi:10.1093/mnras/stw2945](https://doi.org/10.1093/mnras/stw2945)
"""
function computeSpinParameter(
    positions::Matrix{<:Unitful.Length},
    velocities::Matrix{<:Unitful.Velocity},
    masses::Vector{<:Unitful.Mass};
    R::Unitful.Length=FILTER_R,
)::Float64

    (
        !any(isempty, [positions, velocities, masses]) ||
        throw(ArgumentError("computeSpinParameter: The spin parameter of an empty dataset \
        is undefined"))
    )

    # Find cells/particles within the virial radius
    idx = map(x -> x <= R, computeDistance(positions))

    # Compute the total mass within the virial radius
    M = sum(masses[idx]; init=0.0u"Msun")

    # Compute the norm of the total angular momentum
    J = norm(
        computeTotalAngularMomentum(
            positions[:, idx],
            velocities[:, idx],
            masses[idx];
            normal=false,
        )
    )

    return uconvert(Unitful.NoUnits, J / sqrt(2.0 * R * Unitful.G * M^3))

end

"""
    computeGlobalSpinParameter(data_dict::Dict; <keyword arguments>)::Float64

Compute the spin parameter of the whole system in `data`.

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
  - `R::Unitful.Length=FILTER_R`: Radius.

# Returns

  - The spin parameter.
"""
function computeGlobalSpinParameter(data_dict::Dict; R::Unitful.Length=FILTER_R)::Float64

    type_symbols = snapshotTypes(data_dict)

    filter!(ts -> !isempty(data_dict[ts]["POS "]), type_symbols)

    # Concatenate the position and masses of all the cells and particles in the system
    positions  = hcat([data_dict[type_symbol]["POS "] for type_symbol in type_symbols]...)
    velocities = hcat([data_dict[type_symbol]["VEL "] for type_symbol in type_symbols]...)
    masses     = vcat([data_dict[type_symbol]["MASS"] for type_symbol in type_symbols]...)

    @debug("computeGlobalSpinParameter: The spin parameter will be computed using $(type_symbols)")

    # Compute the total spin parameter
    return computeSpinParameter(positions, velocities, masses; R)

end

"""
    computeStellarVcm(data_dict::Dict)::Vector{<:Unitful.Velocity}

Compute the velocity of the stellar center of mass.

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

  - The velocity of the stellar center of mass.
"""
function computeStellarVcm(data_dict::Dict)::Vector{<:Unitful.Velocity}

    # Load the necessary data
    velocities = data_dict[:stars]["VEL "]
    masses = data_dict[:stars]["MASS"]

    # Check for missing data
    !any(isempty, [velocities, masses]) || return zeros(typeof(1.0u"km/s"), 3)

    # Compute the total mass
    M = sum(masses)

    # Compute the velocity of the center of mass
    vcm = [sum(row .* masses) / M for row in eachrow(velocities)]

    return vcm

end

"""
    computeGlobalVcm(data_dict::Dict)::Vector{<:Unitful.Velocity}

Compute the velocity of the global center of mass.

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

  - The velocity of the global center of mass.
"""
function computeGlobalVcm(data_dict::Dict)::Vector{<:Unitful.Velocity}

    # Load the necessary data
    velocities = hcat([data_dict[type_symbol]["VEL "] for type_symbol in type_symbols]...)
    masses     = vcat([data_dict[type_symbol]["MASS"] for type_symbol in type_symbols]...)

    # Check for missing data
    !any(isempty, [velocities, masses]) || return zeros(typeof(1.0u"km/s"), 3)

    # Compute the total mass
    M = sum(masses)

    # Compute the velocity of the center of mass
    vcm = [sum(row .* masses) / M for row in eachrow(velocities)]

    return vcm

end

@doc raw"""
    computeStellarVcirc(
        data_dict::Dict,
    )::Tuple{Vector{<:Unitful.Length},Vector{<:Unitful.Velocity}}

Compute the circular velocity of each stellar particle, with respect to the origin.

The circular velocity of a star is,

```math
v_\mathrm{circ} = \sqrt{\frac{\mathrm{G} \, M(r)}{r}} \, ,
```

where $r$ is the radial distance of the star, and $M(r)$ is the total mass within a sphere of radius $r$.

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

  - A tuple with two elements:

      + A vector with the radial distance of each star to the origin.
      + A vector with the circular velocity of each star.
"""
function computeStellarVcirc(
    data_dict::Dict,
)::Tuple{Vector{<:Unitful.Length},Vector{<:Unitful.Velocity}}

    # Compute the radial distance to each star
    rs = computeDistance(data_dict[:stars]["POS "])

    # Check for missing data
    !isempty(rs) || return rs, Unitful.Velocity[]

    type_symbols = filter!(
        ts -> !isempty(data_dict[ts]["POS "]),
        [:stars, :gas, :halo, :black_hole],
    )

    # Concatenate the position and masses of all the cells and particles in the system
    distances = vcat(
        [computeDistance(data_dict[type_symbol]["POS "]) for
        type_symbol in type_symbols]...,
    )
    masses = vcat([data_dict[type_symbol]["MASS"] for type_symbol in type_symbols]...)

    # Use the radial distances as bin edges for the mass histogram
    edges = [0.0u"kpc", rs...]

    # Compute to total mass within each stellar radial distance
    M = similar(rs, eltype(masses))
    cumsum!(M, histogram1D(distances, masses, edges; empty_nan=false))

    # The mass histogram is a sorted array, so it is reverted to the unsorted order of `r`
    # to make `vcirc` the circular velocity of each star in the order of the snapshot
    invpermute!(M, sortperm(rs))

    @debug("computeStellarVcirc: The circular velocity will be computed using $(type_symbols)")

    vcirc = [iszero(r) ? 0.0u"km*s^-1" : sqrt(Unitful.G * m / r) for (m, r) in zip(M, rs)]

    return rs, vcirc

end

@doc raw"""
    computeStellarCircularity(data_dict::Dict)::Vector{Float64}

Compute the circularity of each stellar particle, with respect to the origin and the $z$ direction [0, 0, 1].

The circularity of a star is,

```math
\epsilon = j_z / j_\mathrm{circ} \, ,
```

where $j_z$ is the $z$ component of its specific angular momentum, and $j_\mathrm{circ}$ is the specific angular momentum of a circular orbit,

```math
j_\mathrm{circ} = r \, v_\mathrm{circ} = \sqrt{\mathrm{G} \, r \, M(r)} \, ,
```

where $r$ is the radial distance of the star, and $M(r)$ is the total mass within a sphere of radius $r$.

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

  - The circularity ``\epsilon`` of each star.
"""
function computeStellarCircularity(data_dict::Dict)::Vector{Float64}

    # Load the necessary data
    positions = data_dict[:stars]["POS "]
    velocities = data_dict[:stars]["VEL "]

    # Check for missing data
    !any(isempty, [positions, velocities]) || return Float64[]

    # Compute the specific angular momentum in the z direction
    jzs = [x[1] * v[2] - x[2] * v[1] for (x, v) in zip(eachcol(positions), eachcol(velocities))]

    # Compute the circular velocities and the radial distances
    rs, vcircs = computeStellarVcirc(data_dict)

    stellar_circularity = [
        any(iszero, [r, vcirc]) ? 0.0 : ustrip(Unitful.NoUnits, jz / (r * vcirc)) for
        (jz, r, vcirc) in zip(jzs, rs, vcircs)
    ]

    return stellar_circularity

end

@doc raw"""
    computeStellarVpolar(data_dict::Dict, component::Symbol)::Vector{<:Unitful.Velocity}

Compute the cylindrical components of the velocity, ``\\mathbf{\\vec{v}} = v_r \, \\mathbf{e_r} + v_\\theta \, \\mathbf{e_\theta} + v_z \, \\mathbf{e_z}``.

The speed in the radial direction expressed in Cartesian coordinates is

```math
v_r = \frac{x \, v_x + y \, v_y}{\sqrt(x^2 + y^2)} \, ,
```

in the tangential direction is

```math
v_\tau = \frac{x \, v_y - y \, v_x}{\sqrt(x^2 + y^2)} \, ,
```

and the speed in the z direction will be computes as

```math
v^*_z = v_z \, \sign(z) \, ,
```

in order to distinguish between inflows (``v^*_z < 0``) and outflows (``v^*_z > 0``).

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
  - `component::Symbol`: Which component will be calculated. The options are:

      + `:radial`     -> Stellar radial speed (``v_r``).
      + `:tangential` -> Stellar tangential speed (``v_\\theta``).
      + `:zstar`      -> Stellar speed in the z direction, computed as ``v_z \\, \\sign(z)``.

# Returns

  - The chosen cylindricall component of the velocity.
"""
function computeStellarVpolar(data_dict::Dict, component::Symbol)::Vector{<:Unitful.Velocity}

    # Load the necessary data
    positions = data_dict[:stars]["POS "]
    velocities = data_dict[:stars]["VEL "]

    x = positions[1, :]
    y = positions[2, :]
    z = positions[3, :]

    vx = velocities[1, :]
    vy = velocities[2, :]
    vz = velocities[3, :]

    if component == :radial

        # Compute the radial component
        vp = @. (x * vx + y * vy) / sqrt(x^2 + y^2)

    elseif component == :tangential

        # Compute the tangential component
        vp = @. (x * vy - y * vx) / sqrt(x^2 + y^2)

    elseif component == :zstar

        # Compute the z component
        vp = @. vz * sign(z)

    else

        throw(ArgumentError("computeStellarVpolar: `component` can only be :radial, :tangential \
        or :zstar, but I got :$(component)"))

    end

    return vp

end

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

  - The radius containing `percet` of the total mass.
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

  - The height containing `percet` of the total mass.
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
    computeMetalMass(data_dict::Dict, type_symbol::Symbol)::Vector{<:Unitful.Mass}

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
  - `type_symbol::Symbol`: For which cell/particle type the metal mass will be calculated. The possibilities are `:stars` and `:gas`.

# Returns

  - The total metal mass in each cell/particle.
"""
function computeMetalMass(data_dict::Dict, type_symbol::Symbol)::Vector{<:Unitful.Mass}

    if type_symbol == :gas
        return setPositive(data_dict[:gas]["GZ  "]) .* data_dict[:gas]["MASS"]
    elseif type_symbol == :stars
        return setPositive(data_dict[:stars]["GZ2 "]) .* data_dict[:stars]["MASS"]
    else
        throw(ArgumentError("computeMetalMass: `type_symbol` can only be :stars or :gas, \
        but I got :$(type_symbol)"))
    end

end

"""
    computeElementMass(
        data_dict::Dict,
        type_symbol::Symbol,
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
  - `type_symbol::Symbol`: For which cell/particle type the element mass will be calculated. The possibilities are `:stars` and `:gas`.
  - `element::Symbol`: Target element. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).

# Returns

  - The total mass of `element` in each cell/particle.
"""
function computeElementMass(
    data_dict::Dict,
    type_symbol::Symbol,
    element::Symbol,
)::Vector{<:Unitful.Mass}

    (
        element ∈ keys(ELEMENT_INDEX) ||
        throw(ArgumentError("computeElementMass: :$(element) is not a tracked element, \
        the options are the keys of `ELEMENT_INDEX`, see `./src/constants.jl`"))
    )

    if type_symbol == :gas
        z_block = "GMET"
    elseif type_symbol == :stars
        z_block = "GME2"
    else
        throw(ArgumentError("computeElementMass: `type_symbol` can only be :stars or :gas, \
        but I got :$(type_symbol)"))
    end

    values = data_dict[type_symbol]

    if any(isempty, [values[z_block], values["MASS"]])
        return Unitful.Mass[]
    end

    return setPositive(values[z_block][ELEMENT_INDEX[element], :]) .* values["MASS"]

end

@doc raw"""
    computeAbundance(
        data_dict::Dict,
        type_symbol::Symbol,
        element::Symbol;
        <keyword arguments>
    )::Float64

Compute the total abundance of a given element, as $n_X / n_H$ where $n_X$ is the number of atoms of element $\mathrm{X}$ and $n_H$ the number of hydrogen atoms.

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
  - `type_symbol::Symbol`: For which cell/particle type the abundance will be calculated. The possibilities are `:stars` and `:gas`.
  - `element::Symbol`: Target element. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
  - `solar::Bool=false`: If the result will be normalized to the solar abundance or not.

# Returns

  - The total abundance of `element`.
"""
function computeGlobalAbundance(
    data_dict::Dict,
    type_symbol::Symbol,
    element::Symbol;
    solar::Bool=false,
)::Float64

    metal_mass = sum(computeElementMass(data_dict, type_symbol, element); init=0.0u"Msun")
    hydrogen_mass = sum(computeElementMass(data_dict, type_symbol, :H); init=0.0u"Msun")

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
    computeIonizedMass(data_dict::Dict)::Vector{<:Unitful.Mass}

Compute the ionized hydrogen mass of every gas cell in `data`.

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

  - The mass of ionized hydrogen in every gas cell.
"""
function computeIonizedMass(data_dict::Dict)::Vector{<:Unitful.Mass}

    dg = data_dict[:gas]

    !(isempty(dg["MASS"]) || isempty(dg["NHP "]) || isempty(dg["NH  "])) || return Unitful.Mass[]

    if "FRAC" ∈ keys(dg) && !isempty(dg["FRAC"])

        # Fraction of ionized hydrogen according to our model
        f_HII = @. dg["FRAC"][1, :] / (1.0 - dg["FRAC"][4, :])

        # When there is no data from our model, use the fraction of ionized hydrogen from Arepo
        fi = [
            isnan(fhii) ? nhp / (nhp + nh) : fhii for
            (fhii, nh, nhp) in zip(f_HII, dg["NH  "], dg["NHP "])
        ]

    else

        # Fraction of ionized hydrogen according to Arepo
        fi = @. dg["NHP "] / (dg["NHP "] + dg["NH  "])

    end

    return fi .* dg["MASS"]

end

"""
    computeAtomicMass(data_dict::Dict)::Vector{<:Unitful.Mass}

Compute the atomic hydrogen mass of every gas cell in `data`.

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

  - The mass of atomic hydrogen in every gas cell.

# References

L. Blitz et al. (2006). *The Role of Pressure in GMC Formation II: The H2-Pressure Relation*. The Astrophysical Journal, **650(2)**, 933. [doi:10.1086/505417](https://doi.org/10.1086/505417)
"""
function computeAtomicMass(data_dict::Dict)::Vector{<:Unitful.Mass}

    dg = data_dict[:gas]

    !(isempty(dg["MASS"]) || isempty(dg["NHP "]) || isempty(dg["NH  "])) || return Unitful.Mass[]

    if "FRAC" ∈ keys(dg) && !isempty(dg["FRAC"])

        # Fraction of atomic hydrogen according to our model
        f_HI = @. dg["FRAC"][2, :] / (1.0 - dg["FRAC"][4, :])

        # When there is no data from our model, use the fraction of neutral hydrogen from Arepo
        # assuming that the fraction of molecular hydrogen is 0
        fa = [
            isnan(fhi) ? nh / (nhp + nh) : fhi for
            (fhi, nh, nhp) in zip(f_HI, dg["NH  "], dg["NHP "])
        ]

    elseif !isempty(dg["PRES"])

        # Fraction of neutral hydrogen according to Arepo
        fn = @. dg["NH  "] / (dg["NHP "] + dg["NH  "])

        # Fraction of molecular hydrogen according to the pressure relation in Blitz et al. (2006)
        fm = @. 1.0 / (1.0 + (P0 / dg["PRES"]))

        # Use the fraction of neutral hydrogen that is not molecular according to the pressure relation,
        # unless that value is negative, in that case use 0 assuming all neutral hydrogen is molecular
        fa = setPositive(fn .- fm)

    else

        return Unitful.Mass[]

    end

    return fa .* dg["MASS"]

end

"""
    computeMolecularMass(data_dict::Dict)::Vector{<:Unitful.Mass}

Compute the molecular hydrogen mass of every gas cell in `data`.

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

  - The mass of molecular hydrogen in every gas cell.

# References

L. Blitz et al. (2006). *The Role of Pressure in GMC Formation II: The H2-Pressure Relation*. The Astrophysical Journal, **650(2)**, 933. [doi:10.1086/505417](https://doi.org/10.1086/505417)
"""
function computeMolecularMass(data_dict::Dict)::Vector{<:Unitful.Mass}

    dg = data_dict[:gas]

    !isempty(dg["MASS"]) || return Unitful.Mass[]

    if "FRAC" ∈ keys(dg) && !isempty(dg["FRAC"])

        # Fraction of molecular hydrogen according to our model
        f_H2 = @. dg["FRAC"][3, :] / (1.0 - dg["FRAC"][4, :])

        # When there is no data from our model, use 0
        fm = replace!(f_H2, NaN => 0.0)

    elseif !isempty(dg["PRES"]) && !isempty(dg["NHP "]) && !isempty(dg["NH  "])

        # Fraction of neutral hydrogen according to Arepo
        fn = @. dg["NH  "] / (dg["NHP "] + dg["NH  "])

        # Fraction of molecular hydrogen according to the pressure relation in Blitz et al. (2006)
        fp = @. 1.0 / (1.0 + (P0 / dg["PRES"]))

        # Use the fraction of molecular hydrogen according to the pressure relation, unless
        # that value is larger than the fraction of neutral hydrogen according to Arepo,
        # in that case use the neutral fraction assuming it is all molecular hydrogen
        fm = [n >= p ? p : n for (n, p) in zip(fn, fp)]

    else

        return Unitful.Mass[]

    end

    return fm .* dg["MASS"]

end

"""
    computeNeutralMass(data_dict::Dict)::Vector{<:Unitful.Mass}

Compute the neutral hydrogen mass of every gas cell in `data`.

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

  - The mass of neutral hydrogen in every gas cell.
"""
function computeNeutralMass(data_dict::Dict)::Vector{<:Unitful.Mass}

    dg = data_dict[:gas]

    !(isempty(dg["MASS"]) || isempty(dg["NHP "]) || isempty(dg["NH  "])) || return Unitful.Mass[]

    if "FRAC" ∈ keys(dg) && !isempty(dg["FRAC"])

        # Fraction of atomic hydrogen according to our model
        f_HI = @. dg["FRAC"][2, :] / (1.0 - dg["FRAC"][4, :])

        # Fraction of molecular hydrogen according to our model
        f_H2 = @. dg["FRAC"][3, :] / (1.0 - dg["FRAC"][4, :])

        # When there is no data from our model, use the fraction of neutral hydrogen from Arepo
        # assuming that the fraction of molecular hydrogen is 0
        fa = [
            isnan(fhi) ? nh / (nhp + nh) : fhi for
            (fhi, nh, nhp) in zip(f_HI, dg["NH  "], dg["NHP "])
        ]

        # When there is no data from our model, use 0
        fm = replace!(f_H2, NaN => 0.0)

        fn = fa .+ fm

    else

        # Fraction of neutral hydrogen according to Arepo
        fn = @. dg["NH  "] / (dg["NHP "] + dg["NH  "])

    end

    return fn .* dg["MASS"]

end

"""
    computeStellarGasMass(data_dict::Dict)::Vector{<:Unitful.Mass}

Compute the "stellar mass" of every gas cell in `data`, which will be different than 0 only for cells in our model.

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

  - The "stellar mass" of every gas cell.
"""
function computeStellarGasMass(data_dict::Dict)::Vector{<:Unitful.Mass}

    dg = data_dict[:gas]

    !isempty(dg["MASS"]) || return Unitful.Mass[]

    if "FRAC" ∈ keys(dg) && !isempty(dg["FRAC"])

        # When there is no data from our model, use 0
        fm = replace!(dg["FRAC"][4, :], NaN => 0.0)

    else

        return zeros(typeof(1.0u"Msun"), length(dg["MASS"]))

    end

    return fm .* dg["MASS"]

end

"""
    function computeSFR(
        data_dict::Dict;
        <keyword arguments>
    )::Vector{<:Unitful.MassFlow}

Compute the star formation rate of each stellar particle in `data`.

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

      + `:stellar_mass`           -> Stellar mass.
      + `:gas_mass`               -> Gas mass.
      + `:dm_mass`                -> Dark matter mass.
      + `:bh_mass`                -> Black hole mass.
      + `:stellar_number`         -> Number of stellar particles.
      + `:gas_number`             -> Number of gas cells.
      + `:dm_number`              -> Number of dark matter particles.
      + `:bh_number`              -> Number of black hole particles.
      + `:molecular_mass`         -> Molecular hydrogen (``\\mathrm{H_2}``) mass.
      + `:atomic_mass`            -> Atomic hydrogen (``\\mathrm{HI}``) mass.
      + `:ionized_mass`           -> Ionized hydrogen (``\\mathrm{HII}``) mass.
      + `:neutral_mass`           -> Neutral hydrogen (``\\mathrm{HI + H_2}``) mass.
      + `:molecular_fraction`     -> Gas mass fraction of molecular hydrogen.
      + `:atomic_fraction`        -> Gas mass fraction of atomic hydrogen.
      + `:ionized_fraction`       -> Gas mass fraction of ionized hydrogen.
      + `:neutral_fraction`       -> Gas mass fraction of neutral hydrogen.
      + `:stellar_area_density`   -> Stellar area mass density, for a radius of `FILTER_R`.
      + `:gas_area_density`       -> Gas area mass density, for a radius of `FILTER_R`.
      + `:molecular_area_density` -> Molecular hydrogen area mass density, for a radius of `FILTER_R`.
      + `:atomic_area_density`    -> Atomic hydrogen area mass density, for a radius of `FILTER_R`.
      + `:ionized_area_density`   -> Ionized hydrogen area mass density, for a radius of `FILTER_R`.
      + `:neutral_area_density`   -> Neutral hydrogen area mass density, for a radius of `FILTER_R`.
      + `:sfr_area_density`       -> Star formation rate area density, for the last `AGE_RESOLUTION_ρ` and a radius of `FILTER_R`.
      + `:gas_metallicity`        -> Mass fraction of all elements above He in the gas (solar units).
      + `:stellar_metallicity`    -> Mass fraction of all elements above He in the stars (solar units).
      + `:X_gas_abundance`        -> Gas abundance of element ``\\mathrm{X}``, as ``12 + \\log_{10}(\\mathrm{X \\, / \\, H})``. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
      + `:X_stellar_abundance`    -> Stellar abundance of element ``\\mathrm{X}``, as ``12 + \\log_{10}(\\mathrm{X \\, / \\, H})``. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
      + `:stellar_specific_am`    -> Norm of the stellar specific angular momentum.
      + `:gas_specific_am`        -> Norm of the gas specific angular momentum.
      + `:dm_specific_am`         -> Norm of the dark matter specific angular momentum.
      + `:sfr`                    -> The star formation rate.
      + `:ssfr`                   -> The specific star formation rate.
      + `:observational_sfr`      -> The star formation rate of the last `AGE_RESOLUTION`.
      + `:observational_ssfr`     -> The specific star formation rate of the last `AGE_RESOLUTION`.
      + `:scale_factor`           -> Scale factor.
      + `:redshift`               -> Redshift.
      + `:physical_time`          -> Physical time since the Big Bang.
      + `:lookback_time`          -> Physical time left to reach the last snapshot.

# Returns

  - The velue of `quantity` for the whole system in `data_dict`.
"""
function integrateQty(data_dict::Dict, quantity::Symbol)::Number

    if quantity == :stellar_mass

        integrated_qty = sum(data_dict[:stars]["MASS"]; init=0.0u"Msun")

    elseif quantity == :gas_mass

        integrated_qty = sum(data_dict[:gas]["MASS"]; init=0.0u"Msun")

    elseif quantity == :dm_mass

        integrated_qty = sum(data_dict[:halo]["MASS"]; init=0.0u"Msun")

    elseif quantity == :bh_mass

        integrated_qty = sum(data_dict[:black_hole]["MASS"]; init=0.0u"Msun")

    elseif quantity == :stellar_number

        integrated_qty = length(data_dict[:stars]["MASS"])

    elseif quantity == :gas_number

        integrated_qty = length(data_dict[:gas]["MASS"])

    elseif quantity == :dm_number

        integrated_qty = length(data_dict[:halo]["MASS"])

    elseif quantity == :bh_number

        integrated_qty = length(data_dict[:black_hole]["MASS"])

    elseif quantity == :molecular_mass

        integrated_qty = sum(computeMolecularMass(data_dict); init=0.0u"Msun")

    elseif quantity == :atomic_mass

        integrated_qty = sum(computeAtomicMass(data_dict); init=0.0u"Msun")

    elseif quantity == :ionized_mass

        integrated_qty = sum(computeIonizedMass(data_dict); init=0.0u"Msun")

    elseif quantity == :neutral_mass

        integrated_qty = sum(computeNeutralMass(data_dict); init=0.0u"Msun")

    elseif quantity == :molecular_fraction

        molecular_mass = sum(computeMolecularMass(data_dict); init=0.0u"Msun")
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

        integrated_qty = sum(data_dict[:stars]["MASS"]; init=0.0u"Msun") / area(FILTER_R)

    elseif quantity == :gas_area_density

        integrated_qty = sum(data_dict[:gas]["MASS"]; init=0.0u"Msun") / area(FILTER_R)

    elseif quantity == :molecular_area_density

        integrated_qty = sum(computeMolecularMass(data_dict); init=0.0u"Msun") / area(FILTER_R)

    elseif quantity == :atomic_area_density

        integrated_qty = sum(computeAtomicMass(data_dict); init=0.0u"Msun") / area(FILTER_R)

    elseif quantity == :ionized_area_density

        integrated_qty = sum(computeIonizedMass(data_dict); init=0.0u"Msun") / area(FILTER_R)

    elseif quantity == :neutral_area_density

        integrated_qty = sum(computeNeutralMass(data_dict); init=0.0u"Msun") / area(FILTER_R)

    elseif quantity == :sfr_area_density

        sfr = sum(computeSFR(data_dict; age_resol=AGE_RESOLUTION_ρ); init=0.0u"Msun*yr^-1")

        integrated_qty = sfr / area(FILTER_R)

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

      + `:stellar_mass`               -> Stellar mass.
      + `:gas_mass`                   -> Gas mass.
      + `:dm_mass`                    -> Dark matter mass.
      + `:bh_mass`                    -> Black hole mass.
      + `:molecular_mass`             -> Molecular hydrogen (``\\mathrm{H_2}``) mass.
      + `:atomic_mass`                -> Atomic hydrogen (``\\mathrm{HI}``) mass.
      + `:ionized_mass`               -> Ionized hydrogen (``\\mathrm{HII}``) mass.
      + `:neutral_mass`               -> Neutral hydrogen (``\\mathrm{HI + H_2}``) mass.
      + `:molecular_fraction`         -> Gas mass fraction of molecular hydrogen.
      + `:atomic_fraction`            -> Gas mass fraction of atomic hydrogen.
      + `:ionized_fraction`           -> Gas mass fraction of ionized hydrogen.
      + `:neutral_fraction`           -> Gas mass fraction of neutral hydrogen.
      + `:molecular_neutral_fraction` -> Fraction of molecular hydrogen in the neutral gas.
      + `:gas_mass_density`           -> Gas mass density.
      + `:gas_number_density`         -> Gas number density.
      + `:molecular_number_density`   -> Molecular hydrogen number density.
      + `:atomic_number_density`      -> Atomic hydrogen number density.
      + `:ionized_number_density`     -> Ionized hydrogen number density.
      + `:neutral_number_density`     -> Neutral hydrogen number density.
      + `:gas_metallicity`            -> Mass fraction of all elements above He in the gas (solar units).
      + `:stellar_metallicity`        -> Mass fraction of all elements above He in the stars (solar units).
      + `:X_gas_abundance`            -> Gas abundance of element ``\\mathrm{X}``, as ``12 + \\log_{10}(\\mathrm{X \\, / \\, H})``. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
      + `:X_stellar_abundance`        -> Stellar abundance of element ``\\mathrm{X}``, as ``12 + \\log_{10}(\\mathrm{X \\, / \\, H})``. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
      + `:stellar_radial_distance`    -> Distance of every stellar particle to the origin.
      + `:gas_radial_distance`        -> Distance of every gas cell to the origin.
      + `:dm_radial_distance`         -> Distance of every dark matter particle to the origin.
      + `:stellar_xy_distance`        -> Projected distance of every stellar particle to the origin.
      + `:gas_xy_distance`            -> Projected distance of every gas cell to the origin.
      + `:dm_xy_distance`             -> Projected distance of every dark matter particle to the origin.
      + `:stellar_circularity`        -> Stellar circularity.
      + `:stellar_vcirc`              -> Stellar circular velocity.
      + `:stellar_vradial`            -> Stellar radial speed.
      + `:stellar_vtangential`        -> Stellar tangential speed.
      + `:stellar_vzstar`             -> Stellar speed in the z direction, computed as ``v_z \\, \\sign(z)``.
      + `:stellar_age`                -> Stellar age.
      + `:sfr`                        -> The star formation rate.
      + `:ssfr`                       -> The specific star formation rate.
      + `:observational_sfr`          -> The star formation rate of the last `AGE_RESOLUTION`.
      + `:observational_ssfr`         -> The specific star formation rate of the last `AGE_RESOLUTION`.
      + `:temperature`                -> Gas temperature, as ``\\log_{10}(T \\, / \\, \\mathrm{K})``.

# Returns

  - The values of `quantity` for every cell/particle.
"""
function scatterQty(data_dict::Dict, quantity::Symbol)::Vector{<:Number}

    if quantity == :stellar_mass

        scatter_qty = data_dict[:stars]["MASS"]

    elseif quantity == :gas_mass

        scatter_qty = data_dict[:gas]["MASS"]

    elseif quantity == :dm_mass

        scatter_qty = data_dict[:halo]["MASS"]

    elseif quantity == :bh_mass

        scatter_qty = data_dict[:black_hole]["MASS"]

    elseif quantity == :molecular_mass

        scatter_qty = computeMolecularMass(data_dict)

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

    elseif quantity == :gas_mass_density

        scatter_qty = data_dict[:gas]["RHO "]

    elseif quantity == :gas_number_density

        scatter_qty = data_dict[:gas]["RHO "] ./ Unitful.mp

    elseif quantity == :molecular_number_density

        molecular_mass = computeMolecularMass(data_dict)
        gas_mass = data_dict[:gas]["MASS"]
        gas_density = data_dict[:gas]["RHO "]

        scatter_qty = gas_density .* (molecular_mass ./ gas_mass) ./ (2 * Unitful.mp)

    elseif quantity == :atomic_number_density

        atomic_mass = computeAtomicMass(data_dict)
        gas_mass = data_dict[:gas]["MASS"]
        gas_density = data_dict[:gas]["RHO "]

        scatter_qty = gas_density .* (atomic_mass ./ gas_mass) ./ Unitful.mp

    elseif quantity == :ionized_number_density

        ionized_mass = computeIonizedMass(data_dict)
        gas_mass = data_dict[:gas]["MASS"]
        gas_density = data_dict[:gas]["RHO "]

        scatter_qty = gas_density .* (ionized_mass ./ gas_mass) ./ Unitful.mp

    elseif quantity == :neutral_number_density

        neutral_mass = computeNeutralMass(data_dict)
        gas_mass = data_dict[:gas]["MASS"]
        gas_density = data_dict[:gas]["RHO "]

        scatter_qty = gas_density .* (neutral_mass ./ gas_mass) ./ Unitful.mp

    elseif quantity == :gas_metallicity

        scatter_qty = setPositive(data_dict[:gas]["GZ  "]) ./ SOLAR_METALLICITY

    elseif quantity == :stellar_metallicity

        scatter_qty = setPositive(data_dict[:stars]["GZ2 "]) ./ SOLAR_METALLICITY

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

    elseif quantity == :stellar_circularity

        @debug("scatterQty: The stellar circularity depends on the positions and velocities of all \
        cell/particles. So, after filtering, the result for a given star will change.")

        scatter_qty = computeStellarCircularity(data_dict)

    elseif quantity == :stellar_vcirc

        @debug("scatterQty: The stellar circular velocity depends on the positions and velocities \
        of all cell/particles. So, after filtering, the result for a given star will change.")

        _, scatter_qty = computeStellarVcirc(data_dict)

    elseif quantity == :stellar_vradial

        scatter_qty = computeStellarVpolar(data_dict, :radial)

    elseif quantity == :stellar_vtangential

        scatter_qty = computeStellarVpolar(data_dict, :tangential)

    elseif quantity == :stellar_vzstar

        scatter_qty = computeStellarVpolar(data_dict, :zstar)

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

    else

        throw(ArgumentError("scatterQty: I don't recognize the quantity :$(quantity)"))

    end

    return scatter_qty

end

"""
    selectFilter(
        filter_mode::Symbol,
        request::Dict{Symbol,Vector{String}},
    )::Tuple{Function,Union{Symbol,NTuple{2,Int},Int},Symbol,Dict{Symbol,Vector{String}}}

Select a filter function, and the corresponding translation and rotation for the simulation box.

Creates a request dictionary, using `request` as a base, adding what is necessary for the filter function and corresponding transformations.

# Arguments

  - `filter_mode::Symbol`: Which cells/particles will be plotted, the options are:

      + `:all`             -> Plot every cell/particle within the simulation box.
      + `:halo`            -> Plot only the cells/particles that belong to the main halo.
      + `:subhalo`         -> Plot only the cells/particles that belong to the main subhalo.
      + `:sphere`          -> Plot only the cell/particle inside a sphere with radius `FILTER_R` (see `./src/constants.jl`).
      + `:stellar_subhalo` -> Plot only the cells/particles that belong to the main subhalo.
      + `:all_subhalo`     -> Plot every cell/particle centered around the main subhalo.
  - `request::Dict{Symbol,Vector{String}}`: Base request dictionary, nothing will be deleted from it.

# Returns

  - A Tuple with four elements:

      + The filter function.
      + Translation for the simulation box. The posibilities are:

          + `:global_cm`                  -> Selects the center of mass of the whole system as the new origin.
          + `:stellar_cm`                 -> Selects the stellar center of mass as the new origin.
          + `(halo_idx, subhalo_rel_idx)` -> Sets the position of the potencial minimum for the `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo), as the new origin.
          + `(halo_idx, 0)`               -> Selects the center of mass of the `halo_idx::Int` halo, as the new origin.
      + Rotation for the simulation box. The posibilities are:

          + `:global_am`          -> Sets the angular momentum of the whole system as the new z axis.
          + `:stellar_am`         -> Sets the stellar angular momentum as the new z axis.
          + `:stellar_pa`         -> Sets the stellar principal axis as the new coordinate system.
          + `:stellar_subhalo_pa` -> Sets the principal axis of the stars in the main subhalo as the new coordinate system.
      + New request dictionary.
"""
function selectFilter(
    filter_mode::Symbol,
    request::Dict{Symbol,Vector{String}},
)::Tuple{Function,Union{Symbol,NTuple{2,Int},Int},Symbol,Dict{Symbol,Vector{String}}}

    if filter_mode == :all

        # Plot every cell/particle
        filter_function = filterNothing
        translation = :global_cm
        rotation = :global_am

        new_request = mergeRequests(
            addRequest(
                request,
                Dict(
                    type_symbol => ["POS ", "MASS", "VEL "] for type_symbol in keys(PARTICLE_INDEX)
                ),
            ),
            Dict(:stars => ["POS ", "MASS", "VEL "]),
        )

    elseif filter_mode == :halo

        # Plot only the cells/particles that belong to the main halo
        filter_function = dd -> filterSubhalo(dd; halo_idx=1, subhalo_rel_idx=0)
        translation = (1, 0)
        rotation = :stellar_pa

        new_request = mergeRequests(
            addRequest(
                request,
                Dict(
                    type_symbol => ["POS ", "MASS", "VEL "] for type_symbol in keys(PARTICLE_INDEX)
                ),
            ),
            Dict(
                :group => ["G_Nsubs", "G_LenType", "G_Pos", "G_Vel"],
                :subhalo => ["S_LenType", "S_Pos", "S_Vel"],
                :stars => ["POS ", "MASS", "VEL "],
            ),
        )

    elseif filter_mode == :subhalo

        # Plot only the cells/particles that belong to the main subhalo
        filter_function = dd -> filterSubhalo(dd; halo_idx=1, subhalo_rel_idx=1)
        translation = (1, 1)
        rotation = :stellar_pa

        new_request = mergeRequests(
            addRequest(
                request,
                Dict(
                    type_symbol => ["POS ", "MASS"] for type_symbol in keys(PARTICLE_INDEX)
                ),
            ),
            Dict(
                :group => ["G_Nsubs", "G_LenType", "G_Pos", "G_Vel"],
                :subhalo => ["S_LenType", "S_Pos", "S_Vel"],
                :stars => ["POS ", "MASS", "VEL "],
            ),
        )

    elseif filter_mode == :sphere

        # Plot only the cell/particle inside a sphere with radius `FILTER_R`
        filter_function = dd -> filterWithin(dd, (0.0u"kpc", FILTER_R), :cm)
        translation = :global_cm
        rotation = :global_am

        new_request = addRequest(
            request,
            Dict(type_symbol => ["POS ", "MASS", "VEL "] for type_symbol in keys(PARTICLE_INDEX)),
        )

    elseif filter_mode == :stellar_subhalo

        # Plot only the cells/particles that belong to the main subhalo
        filter_function = dd -> filterSubhalo(dd; halo_idx=1, subhalo_rel_idx=1)
        translation = :stellar_cm
        rotation = :stellar_pa

        new_request = mergeRequests(
            mergeRequests(request, Dict(:stars => ["POS ", "MASS", "VEL "])),
            Dict(:group => ["G_Nsubs", "G_LenType"], :subhalo => ["S_LenType"]),
        )

    elseif filter_mode == :all_subhalo

        # Plot every cell/particle centered around the main subhalo
        filter_function = filterNothing
        translation = (1, 1)
        rotation = :stellar_subhalo_pa

        new_request = mergeRequests(
            addRequest(
                request,
                Dict(
                    type_symbol => ["POS ", "MASS"] for type_symbol in keys(PARTICLE_INDEX)
                ),
            ),
            Dict(
                :group => ["G_Nsubs", "G_LenType", "G_Pos", "G_Vel"],
                :subhalo => ["S_LenType", "S_Pos", "S_Vel"],
                :stars => ["POS ", "MASS", "VEL "],
            ),
        )

    else

        throw(ArgumentError("selectFilter: `filter_mode` can only be :all, :halo, :subhalo, \
        :stellar_subhalo, :all_subhalo, or :sphere, but I got :$(filter_mode)"))

    end

    return filter_function, translation, rotation, new_request

end

"""
    selectFilter(
        filter_mode::Dict{Symbol,Any},
        request::Dict{Symbol,Vector{String}},
    )::Tuple{
        Function,
        Union{Symbol,NTuple{2,Int},Int},
        Union{Symbol,NTuple{2,Int},Int},
        Dict{Symbol,Vector{String}},
    }

Select a filter function, and the corresponding translation and rotation for the simulation box.

Creates a request dictionary, using `request` as a base, adding what is necessary for the filter function and corresponding transformations.

# Arguments

  - `filter_mode::Dict{Symbol,Any}`: A dictionary with three entries:

      + `:filter_function` -> The filter function.
      + `:translation`     -> Translation for the simulation box. The posibilities are:

          + `:global_cm`                  -> Selects the center of mass of the whole system as the new origin.
          + `:stellar_cm`                 -> Selects the stellar center of mass as the new origin.
          + `(halo_idx, subhalo_rel_idx)` -> Sets the position of the potencial minimum for the `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo) as the new origin.
          + `(halo_idx, 0)`               -> Sets the center of mass of the `halo_idx::Int` halo as the new origin.
          + `subhalo_abs_idx`             -> Sets the center of mass of the `subhalo_abs_idx::Int` as the new origin.
      + `:rotation`        -> Rotation for the simulation box. The posibilities are:

          + `:zero`                       -> No rotation is appplied.
          + `:global_am`                  -> Sets the angular momentum of the whole system as the new z axis.
          + `:stellar_am`                 -> Sets the stellar angular momentum as the new z axis.
          + `:stellar_pa`                 -> Sets the stellar principal axis as the new coordinate system.
          + `:stellar_subhalo_pa`         -> Sets the principal axis of the stars in the main subhalo as the new coordinate system.
          + `(halo_idx, subhalo_rel_idx)` -> Sets the principal axis of the stars in `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo) as the new coordinate system.
          + `(halo_idx, 0)`               -> Sets the principal axis of the stars in the `halo_idx::Int` halo as the new coordinate system.
          + `subhalo_abs_idx`             -> Sets the principal axis of the stars in the `subhalo_abs_idx::Int` subhalo as the new coordinate system.
  - `request::Dict{Symbol,Vector{String}}`: Base request dictionary, nothing will be deleted from it.

# Returns

  - A Tuple with four elements:

      + The filter function.
      + Translation for the simulation box. The posibilities are:

          + `:global_cm`                  -> Selects the center of mass of the whole system as the new origin.
          + `:stellar_cm`                 -> Selects the stellar center of mass as the new origin.
          + `(halo_idx, subhalo_rel_idx)` -> Sets the position of the potencial minimum for the `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo) as the new origin.
          + `(halo_idx, 0)`               -> Sets the center of mass of the `halo_idx::Int` halo as the new origin.
          + `subhalo_abs_idx`             -> Sets the center of mass of the `subhalo_abs_idx::Int` as the new origin.
      + Rotation for the simulation box. The posibilities are:

          + `:zero`                       -> No rotation is appplied.
          + `:global_am`                  -> Sets the angular momentum of the whole system as the new z axis.
          + `:stellar_am`                 -> Sets the stellar angular momentum as the new z axis.
          + `:stellar_pa`                 -> Sets the stellar principal axis as the new coordinate system.
          + `:stellar_subhalo_pa`         -> Sets the principal axis of the stars in the main subhalo as the new coordinate system.
          + `(halo_idx, subhalo_rel_idx)` -> Sets the principal axis of the stars in `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo) as the new coordinate system.
          + `(halo_idx, 0)`               -> Sets the principal axis of the stars in the `halo_idx::Int` halo as the new coordinate system.
          + `subhalo_abs_idx`             -> Sets the principal axis of the stars in the `subhalo_abs_idx::Int` subhalo as the new coordinate system.
      + New request dictionary.
"""
function selectFilter(
    filter_mode::Dict{Symbol,Any},
    request::Dict{Symbol,Vector{String}},
)::Tuple{
    Function,
    Union{Symbol,NTuple{2,Int},Int},
    Union{Symbol,NTuple{2,Int},Int},
    Dict{Symbol,Vector{String}},
}

    new_request = mergeRequests(
        addRequest(
            request,
            Dict(type_symbol => ["POS ", "MASS", "VEL "] for type_symbol in keys(PARTICLE_INDEX)),
        ),
        Dict(
            :group => ["G_Nsubs", "G_LenType", "G_Pos", "G_Vel"],
            :subhalo => ["S_LenType", "S_Pos", "S_Vel"],
            :stars => ["POS ", "MASS", "VEL "],
        ),
    )

    return (
        filter_mode[:filter_function],
        filter_mode[:translation],
        filter_mode[:rotation],
        new_request,
    )

end

####################################################################################################
# Filter functions.
####################################################################################################
#
# A filter function must take a data dictionary, and return a filter dictionary.
#
# These functions are for the second argument of `filterData` in `./src/arepo_utilities.jl`.
#
# Expected signature:
#
#   filter_function(data_dict) -> indices
#
# where:
#
#   - `data_dict::Dict`: A dictionary with the following shape:
#
#      + `:sim_data`          -> ::Simulation (see `Simulation` in `./src/constants.jl`).
#      + `:snap_data`         -> ::Snapshot (see `Snapshot` in `./src/constants.jl`).
#      + `:gc_data`           -> ::GroupCatalog (see `GroupCatalog` in `./src/constants.jl`).
#      + `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
#      + `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
#      + `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
#      + ...
#      + `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
#      + `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
#      + `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
#      + ...
#   - indices::Dict{Symbol,IndexType}: A dictionary with the following shape:
#
#      + `cell/particle type` -> idxs::IndexType
#      + `cell/particle type` -> idxs::IndexType
#      + `cell/particle type` -> idxs::IndexType
#      + ...
#
####################################################################################################

"""
Default filter function that does not filter any cells/particles.
"""
filterNothing(x...; y...)::Dict{Symbol,IndexType} = PASS_ALL

"""
    filterWithin(
        data_dict::Dict,
        range::NTuple{2,<:Unitful.Length},
        origin...,
    )::Dict{Symbol,IndexType}

Filter out the cell/particles outside a given spherical shell.

# Arguments

  - `data::Dict`: A dictionary with the following shape:

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
  - `range::NTuple{2,<:Unitful.Length}`: Internal and external radius of the spherical shell.
  - `origin`: It can be any number and type of argument compatible with the second to last arguments of a [`computeCenter`](@ref) method.

# Returns

  - A dictionary with the following shape:

      + `cell/particle type` -> idxs::IndexType
      + `cell/particle type` -> idxs::IndexType
      + `cell/particle type` -> idxs::IndexType
      + ...
"""
function filterWithin(
    data_dict::Dict,
    range::NTuple{2,<:Unitful.Length},
    origin...,
)::Dict{Symbol,IndexType}

    indices = Dict{Symbol,IndexType}()
    center = computeCenter(data_dict, origin...)

    @inbounds for type_symbol in snapshotTypes(data_dict)

        positions = data_dict[type_symbol]["POS "]

        @inbounds if isempty(positions)
            indices[type_symbol] = (:)
        else
            distances = computeDistance(positions; center)
            indices[type_symbol] = map(x -> range[1] < x <= range[2], distances)
        end

    end

    return indices

end

"""
    filterHotGas(data_dict::Dict, max_temp::Unitful.Temperature)::Dict{Symbol,IndexType}

Filter out gas cells hotter than `max_temp`.

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
  - `max_temp::Unitful.Temperature`: Maximum gas temperature.

# Returns

  - A dictionary with the following shape:

      + `cell/particle type` -> idxs::IndexType
      + `cell/particle type` -> idxs::IndexType
      + `cell/particle type` -> idxs::IndexType
      + ...
"""
function filterHotGas(data_dict::Dict, max_temp::Unitful.Temperature)::Dict{Symbol,IndexType}

    gas_metals = setPositive(data_dict[:gas]["GMET"])
    internal_energy = data_dict[:gas]["U   "]
    electron_fraction = data_dict[:gas]["NE  "]

    # Compute the gas temperature
    temperature = computeTemperature(gas_metals, internal_energy, electron_fraction)

    @inbounds for type_symbol in snapshotTypes(data_dict)

        @inbounds if type_symbol == :gas
            indices[type_symbol] = map(x -> x <= max_temp, temperature)
        else
            indices[type_symbol] = (:)
        end

    end

end

"""
    filterMetallicity(data_dict::Dict, l_Z::Float64, h_Z::Float64)::Dict{Symbol,IndexType}

Filter out gas cells and stellar particles with metallicity outside the range [`l_Z`, `h_Z`].

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
  - `l_Z::Float64`: Minimum metallicity.
  - `h_Z::Float64`: Maximum metallicity.

# Returns

  - A dictionary with the following shape:

      + `cell/particle type` -> idxs::IndexType
      + `cell/particle type` -> idxs::IndexType
      + `cell/particle type` -> idxs::IndexType
      + ...
"""
function filterMetallicity(data_dict::Dict, l_Z::Float64, h_Z::Float64)::Dict{Symbol,IndexType}

    # Allocate memory
    indices = Dict{Symbol,IndexType}()

    @inbounds for type_symbol in snapshotTypes(data_dict)

        @inbounds if type_symbol == :gas
            indices[type_symbol] = map(x -> l_Z <= x <= h_Z, data_dict[:gas]["GZ  "])
        elseif type_symbol == :stars
            indices[type_symbol] = map(x -> l_Z <= x <= h_Z, data_dict[:stars]["GZ2 "])
        else
            indices[type_symbol] = (:)
        end

    end

    return indices

end

"""
    filterCircularity(data_dict::Dict, l_ϵ::Float64, h_ϵ::Float64)::Dict{Symbol,IndexType}

Filter out stellar particles with circularity outside the range [`l_ϵ`, `h_ϵ`].

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
  - `l_ϵ::Float64`: Minimum circularity.
  - `h_ϵ::Float64`: Maximum circularity.

# Returns

  - A dictionary with the following shape:

      + `cell/particle type` -> idxs::IndexType
      + `cell/particle type` -> idxs::IndexType
      + `cell/particle type` -> idxs::IndexType
      + ...
"""
function filterCircularity(data_dict::Dict, l_ϵ::Float64, h_ϵ::Float64)::Dict{Symbol,IndexType}

    # Allocate memory
    indices = Dict{Symbol,IndexType}()

    @inbounds for type_symbol in snapshotTypes(data_dict)

        @inbounds if type_symbol == :stars
            circularity = computeStellarCircularity(data_dict)
            indices[type_symbol] = map(x -> l_ϵ <= x <= h_ϵ, circularity)
        else
            indices[type_symbol] = (:)
        end

    end

    return indices

end

"""
    filterSubhalo(
        data_dict::Dict;
        <keyword arguments>
    )::Dict{Symbol,IndexType}

Filter out cells/particles that do not belong to a given halo and subhalo.

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
  - `halo_idx::Int`: Index of the target halo (FoF group). Starts at 1.
  - `subhalo_rel_idx::Int`: Index of the target subhalo (subfind), relative the target halo. Starts at 1. If set to 0, all subhalos of the target halo are included.

# Returns

  - A dictionary with the following shape:

      + `cell/particle type` -> idxs::IndexType
      + `cell/particle type` -> idxs::IndexType
      + `cell/particle type` -> idxs::IndexType
      + ...
"""
function filterSubhalo(
    data_dict::Dict;
    halo_idx::Int=1,
    subhalo_rel_idx::Int=1,
)::Dict{Symbol,IndexType}

    # If there are no subfind data, filter out every cell/particle
    if ismissing(data_dict[:gc_data].path) && !isSubfindActive(data_dict[:gc_data].path)
        return PASS_NONE
    end

    # Load the necessary data
    g_n_subs = data_dict[:group]["G_Nsubs"]
    g_len_type = data_dict[:group]["G_LenType"]
    s_len_type = data_dict[:subhalo]["S_LenType"]

    # If any of the data is misssing return an empty filter dictionary
    n_groups_total = data_dict[:gc_data].header.n_groups_total
    (
        !iszero(n_groups_total) && !any(isempty, [g_n_subs, g_len_type, s_len_type]) ||
        return Dict(type_symbol => Int[] for type_symbol in snapshotTypes(data_dict))
    )

    # Check that the requested halo index is within bounds
    (
        0 < halo_idx <= n_groups_total ||
        throw(ArgumentError("filterSubhalo: There is only $(n_groups_total) FoF goups in \
        $(data_dict[:gc_data].path), so `halo_idx` = $(halo_idx) is out of bounds"))
    )

    # Compute the number of subhalos and particles up to the last halo before `halo_idx`
    if isone(halo_idx)
        n_subs_floor = 0
        len_type_floor = zeros(Int32, size(s_len_type, 1))
    else
        n_subs_floor = sum(g_n_subs[1:(halo_idx - 1)]; init=0)
        len_type_floor = sum(g_len_type[:, 1:(halo_idx - 1)], dims=2; init=0)
    end

    # Check that the requested subhalo index is within bounds
    n_subfinds = g_n_subs[halo_idx]
    (
        subhalo_rel_idx <= n_subfinds ||
        throw(ArgumentError("filterSubhalo: There is only $(n_subfinds) subhalos for the \
        FoF group $(halo_idx) in $(data_dict[:gc_data].path), so `subhalo_rel_idx` = \
        $(subhalo_rel_idx) is out of bounds"))
    )

    if subhalo_rel_idx <= 0

        # Consider all subhalos within the target halo
        first_idxs = len_type_floor .+ 1
        last_idxs  = first_idxs .+ g_len_type[:, halo_idx] .- 1

    else

        # Compute the subhalo absolute index
        subhalo_abs_idx = n_subs_floor + subhalo_rel_idx

        # Compute the number of particles in the current halo,
        # upto the last subhalo before `subhalo_rel_idx`
        if isone(subhalo_abs_idx)
            len_type_floor_in_halo = zeros(Int, size(s_len_type, 1))
        else
            len_type_floor_in_halo = sum(
                s_len_type[:, (n_subs_floor + 1):(subhalo_abs_idx - 1)], dims=2; init=0,
            )
        end

        # Compute the first and last index of the selected
        # cells/particles (for each cell/particle type)
        first_idxs = len_type_floor .+ len_type_floor_in_halo .+ 1
        last_idxs  = first_idxs .+ s_len_type[:, subhalo_abs_idx] .- 1

    end

    # Allocate memory
    indices = Dict{Symbol,IndexType}()

    # Fill the filter dictionary
    @inbounds for (i, (first_idx, last_idx)) in enumerate(zip(first_idxs, last_idxs))

        type_symbol = INDEX_PARTICLE[i - 1]

        @inbounds if first_idx == last_idx || iszero(last_idx)
            indices[type_symbol] = Int[]
        end

        if type_symbol == :stars

            # Find the indices of the stars, excluding wind particles
            real_stars_idxs = findRealStars(data_dict[:snap_data].path)

            n_wind_before = count(x -> !(x), real_stars_idxs[1:(first_idx - 1)])
            n_wind_between = count(x -> !(x), real_stars_idxs[first_idx:last_idx])

            stars_first_idx = first_idx - n_wind_before
            stars_last_idx = last_idx - n_wind_before - n_wind_between

            indices[type_symbol] = stars_first_idx:stars_last_idx

        else

            indices[type_symbol] = first_idx:last_idx

        end

    end

    return indices

end

"""
    filterSubhalo(data_dict::Dict, subhalo_abs_idx::Int)::Dict{Symbol,IndexType}

Filter out cells/particles that do not belong to a given subhalo.

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
  - `subhalo_abs_idx::Int`: Index of the target subhalo (subfind). Starts at 1.

# Returns

  - A dictionary with the following shape:

      + `cell/particle type` -> idxs::IndexType
      + `cell/particle type` -> idxs::IndexType
      + `cell/particle type` -> idxs::IndexType
      + ...
"""
function filterSubhalo(data_dict::Dict, subhalo_abs_idx::Int)::Dict{Symbol,IndexType}

    # If there are no subfind data, filter out every cell/particle
    if ismissing(data_dict[:gc_data].path) && !isSubfindActive(data_dict[:gc_data].path)
        return PASS_NONE
    end

    # Load the necessary data
    s_len_type = data_dict[:subhalo]["S_LenType"]

    # If any of the data is misssing return an empty filter dictionary
    n_subgroups_total = data_dict[:gc_data].header.n_subgroups_total
    (
        !iszero(n_subgroups_total) && !isempty(s_len_type) ||
        return Dict(type_symbol => Int[] for type_symbol in snapshotTypes(data_dict))
    )

    # Check that the requested subhalo index is within bounds
    (
        0 < subhalo_abs_idx <= n_subgroups_total ||
        throw(ArgumentError("filterSubhalo: There is only $(n_subgroups_total) subhalos in \
        $(data_dict[:gc_data].path), so subhalo_abs_idx = $(subhalo_abs_idx) is out of bounds"))
    )

    # Compute the number of particles upto the last subhalo before `subhalo_abs_idx`
    if isone(subhalo_abs_idx)
        len_type_floor = zeros(Int, size(s_len_type, 1))
    else
        len_type_floor = sum(s_len_type[:, 1:(subhalo_abs_idx - 1)], dims=2; init=0)
    end

    # Compute the first and last index of the selected cells/particles (for each cell/particle type)
    first_idxs = len_type_floor .+ 1
    last_idxs  = len_type_floor .+ s_len_type[:, subhalo_abs_idx]

    # Allocate memory
    indices = Dict{Symbol,IndexType}()

    # Fill the filter dictionary
    @inbounds for (i, (first_idx, last_idx)) in enumerate(zip(first_idxs, last_idxs))

        type_symbol = INDEX_PARTICLE[i - 1]

        @inbounds if first_idx == last_idx || iszero(last_idx)
            indices[type_symbol] = Int[]
        end

        if type_symbol == :stars

            # Find the indices of the stars, excluding wind particles
            real_stars_idxs = findRealStars(data_dict[:snap_data].path)

            n_wind_before = count(x -> !(x), real_stars_idxs[1:(first_idx - 1)])
            n_wind_between = count(x -> !(x), real_stars_idxs[first_idx:last_idx])

            stars_first_idx = first_idx - n_wind_before
            stars_last_idx = last_idx - n_wind_before - n_wind_between

            indices[type_symbol] = stars_first_idx:stars_last_idx

        else

            indices[type_symbol] = first_idx:last_idx

        end

    end

    return indices

end

"""
    filterZinSubhalo(
        data_dict::Dict,
        l_Z::Float64,
        h_Z::Float64;
        <keyword arguments>
    )::Dict{Symbol,IndexType}

Filter out gas cells and stellar particles with metallicity outside the range [`l_Z`, `h_Z`], and every cell/particle that does not belong to a given halo and subhalo.

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
  - `l_Z::Float64`: Minimum metallicity.
  - `h_Z::Float64`: Maximum metallicity.
  - `halo_idx::Int`: Index of the target halo (FoF group). Starts at 1.
  - `subhalo_rel_idx::Int`: Index of the target subhalo (subfind), relative the target halo. Starts at 1. If set to 0, all subhalos of the target halo are included.

# Returns

  - A dictionary with the following shape:

      + `cell/particle type` -> idxs::IndexType
      + `cell/particle type` -> idxs::IndexType
      + `cell/particle type` -> idxs::IndexType
      + ...
"""
function filterZinSubhalo(
    data_dict::Dict,
    l_Z::Float64,
    h_Z::Float64;
    halo_idx::Int=1,
    subhalo_rel_idx::Int=1,
)::Dict{Symbol,IndexType}

    # Find indices to filter by halo and subhalo
    subhalo_idxs = filterSubhalo(data_dict; halo_idx, subhalo_rel_idx)

    # Find indices to filter by metallicity
    metallicity_idxs = filterMetallicity(data_dict, l_Z, h_Z)

    # Allocate memory
    indices = Dict{Symbol,IndexType}()

    @inbounds for type_symbol in snapshotTypes(data_dict)
        indices[type_symbol] = subhalo_idxs[type_symbol] ∩ metallicity_idxs[type_symbol]
    end

    return indices

end

"""
    filterComponentinSubhalo(
        data_dict::Dict,
        component::Symbol;
        <keyword arguments>
    )::Dict{Symbol,IndexType}

Filter out stellar particles that do not belong to the given morphological component, based on a circularity criteria.

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
  - `component::Symbol`: Target component. It can be:
      + `:disk`  -> Stellar particles with a circularity larger than 0.7.
      + `:bulge` -> Stellar particles with a circularity smaller than 0.7.
  - `halo_idx::Int`: Index of the target halo (FoF group). Starts at 1.
  - `subhalo_rel_idx::Int`: Index of the target subhalo (subfind), relative the target halo. Starts at 1. If set to 0, all subhalos of the target halo are included.

# Returns

  - A dictionary with the following shape:

      + `cell/particle type` -> idxs::IndexType
      + `cell/particle type` -> idxs::IndexType
      + `cell/particle type` -> idxs::IndexType
      + ...
"""
function filterComponentinSubhalo(
    data_dict::Dict,
    component::Symbol;
    halo_idx::Int=1,
    subhalo_rel_idx::Int=1,
)::Dict{Symbol,IndexType}

    # Find indices to filter by halo and subhalo
    subhalo_idxs = filterSubhalo(data_dict; halo_idx, subhalo_rel_idx)

    # Find indices to filter by component
    if component == :bulge
        component_idxs = filterCircularity(data_dict, 0.7, Inf)
    elseif component == :disk
        component_idxs = filterCircularity(data_dict, -Inf, 0.7)
    else
        throw(ArgumentError("filterComponentinSubhalo: `component` can only be :disk or :bulge, \
        but I got :$(component)"))
    end

    # Allocate memory
    indices = Dict{Symbol,IndexType}()

    @inbounds for type_symbol in snapshotTypes(data_dict)
        indices[type_symbol] = subhalo_idxs[type_symbol] ∩ component_idxs[type_symbol]
    end

    return indices

end
