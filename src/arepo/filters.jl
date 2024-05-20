####################################################################################################
# Filters.
####################################################################################################

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
  - `filter_function::Function=filterNothing`: A function with the signature:

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
  - `filter_function::Function=filterNothing`: A function with the signature:

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
    selectFilter(
        filter_mode::Symbol,
        request::Dict{Symbol,Vector{String}},
    )::Tuple{Function,Union{Symbol,NTuple{2,Int}},Symbol,Dict{Symbol,Vector{String}}}

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
)::Tuple{Function,Union{Symbol,NTuple{2,Int}},Symbol,Dict{Symbol,Vector{String}}}

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
                :group   => ["G_Nsubs", "G_LenType", "G_Pos", "G_Vel"],
                :subhalo => ["S_LenType", "S_Pos", "S_Vel"],
                :stars   => ["POS ", "MASS", "VEL "],
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
                :group   => ["G_Nsubs", "G_LenType", "G_Pos", "G_Vel"],
                :subhalo => ["S_LenType", "S_Pos", "S_Vel"],
                :stars   => ["POS ", "MASS", "VEL "],
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
                :group   => ["G_Nsubs", "G_LenType", "G_Pos", "G_Vel"],
                :subhalo => ["S_LenType", "S_Pos", "S_Vel"],
                :stars   => ["POS ", "MASS", "VEL "],
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
            :group   => ["G_Nsubs", "G_LenType", "G_Pos", "G_Vel"],
            :subhalo => ["S_LenType", "S_Pos", "S_Vel"],
            :stars   => ["POS ", "MASS", "VEL "],
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
    filterGFM(data_dict::Dict)::Dict{Symbol,IndexType}

Filter out gas cells that have not entered out routine.

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

  - A dictionary with the following shape:

      + `cell/particle type` -> idxs::IndexType
      + `cell/particle type` -> idxs::IndexType
      + `cell/particle type` -> idxs::IndexType
      + ...
"""
function filterGFM(data_dict::Dict)::Dict{Symbol,IndexType}

    # Allocate memory
    indices = Dict{Symbol,IndexType}()

    @inbounds for type_symbol in snapshotTypes(data_dict)

        @inbounds if type_symbol == :gas
            if isempty(data_dict[:gas]["FRAC"])
                indices[type_symbol] = Int[]
            else
                indices[type_symbol] = map(!isnan, data_dict[:gas]["FRAC"][1, :])
            end
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

Filter out stellar particles that do not belong to the given morphological component of a given halo and subhalo, based on a circularity criteria.

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

"""
    filterRadioinSubhalo(
        data_dict::Dict,
        range::NTuple{2,<:Unitful.Length};
        <keyword arguments>
    )::Dict{Symbol,IndexType}

Filter out the cell/particles outside a given spherical shell, around a given halo and subhalo.

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
  - `range::NTuple{2,<:Unitful.Length}`: Internal and external radius of the spherical shell.
  - `halo_idx::Int`: Index of the target halo (FoF group). Starts at 1.
  - `subhalo_rel_idx::Int`: Index of the target subhalo (subfind), relative the target halo. Starts at 1. If set to 0, all subhalos of the target halo are included.

# Returns

  - A dictionary with the following shape:

      + `cell/particle type` -> idxs::IndexType
      + `cell/particle type` -> idxs::IndexType
      + `cell/particle type` -> idxs::IndexType
      + ...
"""
function filterRadioinSubhalo(
    data_dict::Dict,
    range::NTuple{2,<:Unitful.Length};
    halo_idx::Int=1,
    subhalo_rel_idx::Int=1,
)::Dict{Symbol,IndexType}

    # Find indices to filter by halo and subhalo
    subhalo_idxs = filterSubhalo(data_dict; halo_idx, subhalo_rel_idx)

    # Find indices to filter by radial distance
    radio_idxs = filterWithin(data_dict, range, (halo_idx, subhalo_rel_idx))

    # Allocate memory
    indices = Dict{Symbol,IndexType}()

    @inbounds for type_symbol in snapshotTypes(data_dict)
        indices[type_symbol] = subhalo_idxs[type_symbol] ∩ radio_idxs[type_symbol]
    end

    return indices

end
