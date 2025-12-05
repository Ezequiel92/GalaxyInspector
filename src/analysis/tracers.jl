####################################################################################################
# Tracer functions
####################################################################################################

##################
# ID translations
##################

"""
    parentToTracerID(data_dict::Dict, target_ids::Vector{<:Unsigned})::Vector{UInt}

Find the tracer IDs given a list of parent IDs.

# Arguments

  - `data_dict::Dict`: Data dictionary (see [`makeDataDict`](@ref) for the canonical description).
  - `target_ids::Vector{<:Unsigned}`: List of parent IDs.

# Returns

  - A vector with the IDs of the tracers. There is no guarantee length or order for this vector.
"""
function parentToTracerID(data_dict::Dict, target_ids::Vector{<:Unsigned})::Vector{UInt}

    # Read the full list of parent and tracer IDs
    parent_ids = data_dict[:tracer]["PAID"]
    tracer_ids = data_dict[:tracer]["TRID"]

    # Convert targets to a Set for O(1) membership testing
    target_set = Set(target_ids)

    found_tracers = Vector{UInt}()

    # Find the indices of the target IDs in the parent ID list
    for (i, tracer_id) in pairs(tracer_ids)
        if parent_ids[i] in target_set
            push!(found_tracers, tracer_id)
        end
    end

    return found_tracers
end

"""
    tracerToParentID(data_dict::Dict, target_ids::Vector{<:Unsigned})::Vector{UInt}

Find the parent IDs given a list of tracer IDs.

# Arguments

  - `data_dict::Dict`: Data dictionary (see [`makeDataDict`](@ref) for the canonical description).
  - `target_ids::Vector{<:Unsigned}`: List of tracer IDs.

# Returns

  - A vector with the IDs of the parents. There is no guarantee length or order for this vector.
"""
function tracerToParentID(data_dict::Dict, target_ids::Vector{<:Unsigned})::Vector{UInt}

    # Read the full list of parent and tracer IDs
    parent_ids = data_dict[:tracer]["PAID"]
    tracer_ids = data_dict[:tracer]["TRID"]

    # Convert targets to a Set for O(1) membership testing
    target_set = Set(target_ids)

    found_parents = Vector{UInt}()

    # Find the indices of the target IDs in the tracer ID list
    for (i, parent_id) in pairs(parent_ids)
        if tracer_ids[i] in target_set
            push!(found_parents, parent_id)
        end
    end

    return unique!(found_parents)

end

"""
    idToIndex(data_dict::Dict, target_ids::Vector{<:Unsigned})::Dict{Symbol,Vector{Int}}

Find the indices of the cells/particles given their IDs.

# Arguments

  - `data_dict::Dict`: Data dictionary (see [`makeDataDict`](@ref) for the canonical description).
  - `target_ids::Vector{<:Unsigned}`: List of IDs.

# Returns

  - A dictionary with the following shape:

      + `cell/particle type::Symbol` -> indices::Int
      + `cell/particle type::Symbol` -> indices::Int
      + `cell/particle type::Symbol` -> indices::Int
      + ...
"""
function idToIndex(data_dict::Dict, target_ids::Vector{<:Unsigned})::Dict{Symbol,Vector{Int}}

    index_dict = Dict{Symbol,Vector{Int}}()

    # Convert targets to a Set for O(1) membership testing
    target_set = Set(target_ids)

    for component in snapshotTypes(data_dict)

        data = data_dict[component]

        (
            !haskey(data, "ID  ") &&
            throw(ArgumentError("idToIndex: The ID block is missing for component $(component) in \
            $(data_dict[:snap_data].path)"))
        )

        # Read the IDs of the cells/particles
        ids = data["ID  "]

        idxs = Int[]

        if !isempty(ids)

            # Find the indices of the target IDs in the cell/particle ID list
            for (i, id) in pairs(ids)
                if id in target_set
                    push!(idxs, i)
                end
            end

        end

        index_dict[component] = idxs

    end

    return index_dict

end

"""
    tracersToParentMass(
        data_dict::Dict,
        target_ids::Vector{<:Unsigned},
    )::Dict{Symbol,Vector{<:Unitful.Mass}}

Find the masses of the parents given a list of tracer IDs.

# Arguments

  - `data_dict::Dict`: Data dictionary (see [`makeDataDict`](@ref) for the canonical description).
  - `target_ids::Vector{<:Unsigned}`: List of tracer IDs.

# Returns

  - A dictionary with the following shape:

      + `cell/particle type` -> masses::Vector{<:Unitful.Mass}
      + `cell/particle type` -> masses::Vector{<:Unitful.Mass}
      + `cell/particle type` -> masses::Vector{<:Unitful.Mass}
      + ...
"""
function tracersToParentMass(
    data_dict::Dict,
    target_ids::Vector{<:Unsigned},
)::Dict{Symbol,Vector{<:Unitful.Mass}}

    # Find the parent ID corresponding to each tracer
    parent_ids = tracerToParentID(data_dict, target_ids)

    # Find the index and cell/particle type of each parent
    index_dict = idToIndex(data_dict, parent_ids)

    return Dict(component => data_dict[component]["MASS"][idx] for (component, idxs) in index_dict)

end

"""
    parentToTracerMass(
        data_dict::Dict,
        target_ids::Vector{<:Unsigned},
    )::Unitful.Mass

Find the mass of tracers given a list of parent IDs, ignoring parents without tracers

# Arguments

  - `data_dict::Dict`: Data dictionary (see [`makeDataDict`](@ref) for the canonical description).
  - `target_ids::Vector{<:Unsigned}`: List of parent IDs.

# Returns

  - The total mass of tracers.
"""
function parentToTracerMass(
    data_dict::Dict,
    target_ids::Vector{<:Unsigned},
)::Unitful.Mass

    # Read the full list of parent IDs
    parent_ids = data_dict[:tracer]["PAID"]

    # Convert targets to a Set for O(1) membership testing
    target_set = Set(target_ids)

    tracer_count = count(in(target_set), parent_ids)

    iszero(tracer_count) && return 0.0u"Msun"

    # Compute the mass of each tracer in physical units
    tracer_mass = TRACER_MASS * internalUnits("MASS", data_dict[:snap_data].path)

    return tracer_count * tracer_mass

end

##############
# ID tracking
##############

"""
    idWithinR200(data_dict::Dict, component::Symbol; <keyword arguments>)::Vector{UInt}

Find the IDs the the cell/particles of `component` that are within the virial radius of halo `halo_idx`.

!!! note

    This function assumes that no translation of the coordinate system has been done. It uses as the center of the halo the position given by the group catalog files, which are unaffected by coordinate transformation of the data.


# Arguments

  - `data_dict::Dict`: Data dictionary (see [`makeDataDict`](@ref) for the canonical description).
  - `component::Symbol`: Target component. It can only be one of the elements of [`COMPONENTS`](@ref).
  - `halo_idx::Int=1`: Index of the target halo (FoF group). Starts at 1.

# Returns

  - A vector with the IDs.
"""
function idWithinR200(data_dict::Dict, component::Symbol; halo_idx::Int=1)::Vector{UInt}

    if isempty(data_dict[:group]["G_R_Crit200"])

        (
            logging[] &&
            @warn("idWithinR200: There is missing group data in $(data_dict[:gc_data].path), \
            so I will return an empty ID list")
        )

        return UInt[]

    end

    # Read the virial radius
    r200 = data_dict[:group]["G_R_Crit200"][halo_idx]

    if iszero(r200)

        (
            logging[] &&
            @warn("idWithinR200: The virial radius of halo index $(halo_idx) is zero in \
            $(data_dict[:gc_data].path), so I will return an empty ID list")
        )

        return UInt[]

    end

    # Compute the center of the halo
    center = computeCenter(data_dict, halo_idx, 0)

    # Read the data of the given component
    data = data_dict[component]

    (
        !(haskey(data, "ID  ") && haskey(data, "POS ")) &&
        throw(ArgumentError("idWithinR200: The ID and/or position blocks are missing for component \
        $(component) in $(data_dict[:snap_data].path)"))
    )

    # Read the positions of the cells and particles
    positions = data["POS "]

    # Read the IDs of the cells and particles
    ids = data["ID  "]

    if !isempty(positions) && !isempty(ids)

        # Compute the distances of the cells and particles to the center of the halo
        distances = colwise(Euclidean(), positions, center)

        # Find the indices of the cells and particles within `r200`
        idxs = map(x -> x <= r200, distances)

        # Read the IDs of the cells and particles within `r200`
        return ids[idxs]

    end

    return UInt[]

end

"""
    idWithinDisk(
        data_dict::Dict,
        component::Symbol,
        max_r::Unitful.Length,
        max_z::Unitful.Length,
        origin...,
    )::Vector{UInt}

Find the IDs the the cell/particles of `component` that are within the given galactic disk.

# Arguments

  - `data_dict::Dict`: Data dictionary (see [`makeDataDict`](@ref) for the canonical description).
  - `component::Symbol`: Target component. It can only be one of the elements of [`COMPONENTS`](@ref).
  - `max_r::Unitful.Length`: Radius of the disk.
  - `max_z::Unitful.Length`: Half height of the disk.
  - `origin...`: Center of the disk. It can be any number and type of argument compatible with the second to last arguments of a [`computeCenter`](@ref) method.

# Returns

  - A vector with the IDs.
"""
function idWithinDisk(
    data_dict::Dict,
    component::Symbol,
    max_r::Unitful.Length,
    max_z::Unitful.Length,
    origin...,
)::Vector{UInt}

    (
        isPositive([max_r, max_z]) ||
        throw(ArgumentError("idWithinDisk: `max_r` and `max_z` should be larger than 0, \
        but I got `max_r` = $(max_r) and `max_z` = $(max_z)"))
    )

    # Compute the position of the center of the disk
    center = computeCenter(data_dict, origin...)

    # Read the data of the given component
    data = data_dict[component]

    (
        !(haskey(data, "ID  ") && haskey(data, "POS ")) &&
        throw(ArgumentError("idWithinDisk: The ID and/or position blocks are missing for component \
        $(component) in $(data_dict[:snap_data].path)"))
    )

    # Read the positions of the cells and particles
    positions = data["POS "]

    # Read the IDs of the cells and particles
    ids = data["ID  "]

    if !isempty(positions) && !isempty(ids)

        # Use the square of the maximum distance to avoid sqrt() computation
        max_r_sq = max_r^2

        idxs = BitVector(undef, size(positions, 2))

        Threads.@threads for i in axes(positions, 2)

            if abs(positions[3, i]) <= max_z

                dx = positions[1, i] - center[1]
                dy = positions[2, i] - center[2]

                idxs[i] = (dx^2 + dy^2 <= max_r_sq)

            else

                idxs[i] = false

            end

        end

        return ids[idxs]

    end

    return UInt[]

end

"""
    filterExistIDs!(target_ids::Vector{UInt}, component::Symbol, data_dict::Dict)::Nothing

Filter out IDs from `target_ids` that do not exist for `component`.

# Arguments

  - `data_dict::Dict`: Data dictionary (see [`makeDataDict`](@ref) for the canonical description).
  - `component::Symbol`: Target component. It can only be one of the elements of [`COMPONENTS`](@ref).
  - `target_ids::Vector{UInt}`: IDs to be filtered.

# Returns

  - A vector with only the IDs of `target_ids` that exist for `component`.
"""
function filterExistIDs!(data_dict::Dict, component::Symbol, target_ids::Vector{UInt})::Nothing

    data = data_dict[component]

    (
        !haskey(data, "ID  ") &&
        throw(ArgumentError("filterExistIDs!: The ID block is missing for component $(component) \
        in $(data_dict[:snap_data].path)"))
    )

    data_ids = data["ID  "]

    if !isempty(data_ids)

        valid_ids = Set{UInt}(data_ids)

        filter!(id -> id in valid_ids, target_ids)

    else

        empty!(target_ids)

    end

    return nothing

end

"""
    idMass!(
        data_dict::Dict,
        component::Symbol,
        target_ids::Vector{UInt};
        <keyword arguments>
    )::Unitful.Mass

Compute the total mass of the cells/particles with IDs given by `target_ids`.

# Arguments

  - `data_dict::Dict`: Data dictionary (see [`makeDataDict`](@ref) for the canonical description).
  - `component::Symbol`: Target component. It can only be one of the elements of [`COMPONENTS`](@ref).
  - `target_ids::Vector{UInt}`: Target IDs.
  - `tracers::Bool=false`: If the mass will be computed using tracers or the proper mass of each cell/particle.

# Returns

  - The total mass.
"""
function idMass!(
    data_dict::Dict,
    component::Symbol,
    target_ids::Vector{UInt};
    tracers::Bool=false,
)::Unitful.Mass

    isempty(target_ids) && return 0.0u"Msun"

    tracers && return parentToTracerMass(data_dict, target_ids)

    data = data_dict[component]

    (
        !(haskey(data, "ID  ") && haskey(data, "MASS")) &&
        throw(ArgumentError("idMass!: The ID and/or mass blocks are missing for component \
        $(component) in $(data_dict[:snap_data].path)"))
    )

    # Read the IDs and masses of the cells/particles
    ids = data["ID  "]
    masses = data["MASS"]

    isempty(ids) && return 0.0u"Msun"
    isempty(masses) && return 0.0u"Msun"

    # Convert targets to a Set for O(1) membership testing
    target_set = Set(target_ids)

    total_mass = 0.0u"Msun"

    @inbounds @simd for i in eachindex(ids)
        if ids[i] in target_set
            total_mass += masses[i]
        end
    end

    return total_mass

end
