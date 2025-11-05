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

  - A vector with the IDs of the tracers.
"""
function parentToTracerID(data_dict::Dict, target_ids::Vector{<:Unsigned})::Vector{UInt}

    # Read the full list of parent and tracer IDs
    parent_ids = data_dict[:tracer]["PAID"]
    tracer_ids = data_dict[:tracer]["TRID"]

    # Find the indices of the target IDs in the parent ID list, ignoring invalid targets
    idxs = filter!(!isnothing, indexin(target_ids, parent_ids))

    return tracer_ids[idxs]

end

"""
    tracerToParentID(data_dict::Dict, target_ids::Vector{<:Unsigned})::Vector{UInt}

Find the parent IDs given a list of tracer IDs.

# Arguments

  - `data_dict::Dict`: Data dictionary (see [`makeDataDict`](@ref) for the canonical description).
  - `target_ids::Vector{<:Unsigned}`: List of tracer IDs.

# Returns

  - A vector with the IDs of the parents.
"""
function tracerToParentID(data_dict::Dict, target_ids::Vector{<:Unsigned})::Vector{UInt}

    # Read the full list of parent and tracer IDs
    parent_ids = data_dict[:tracer]["PAID"]
    tracer_ids = data_dict[:tracer]["TRID"]

    # Find the indices of the target IDs in the tracer ID list, ignoring invalid targets
    idxs = filter!(!isnothing, indexin(target_ids, tracer_ids))

    return parent_ids[idxs]

end

"""
    parentIDToIndex(data_dict::Dict, target_ids::Vector{<:Unsigned})::Dict{Symbol,Vector{Int}}

Find the indices of the cells/particles given their IDs.

# Arguments

  - `data_dict::Dict`: Data dictionary (see [`makeDataDict`](@ref) for the canonical description).
  - `target_ids::Vector{<:Unsigned}`: List of parent IDs.

# Returns

  - A dictionary with the following shape:

      + `cell/particle type::Symbol` -> indices::Int
      + `cell/particle type::Symbol` -> indices::Int
      + `cell/particle type::Symbol` -> indices::Int
      + ...
"""
function parentIDToIndex(data_dict::Dict, target_ids::Vector{<:Unsigned})::Dict{Symbol,Vector{Int}}

    index_dict = Dict(component => Int[] for component in snapshotTypes(data_dict))

    for component in snapshotTypes(data_dict)

        data = data_dict[component]

        if haskey(data, "ID  ") && !isempty(data["ID  "])

            # Read the IDs of the cells/particles
            parent_ids = data["ID  "]

            # Find the indices of the target IDs in the cell/particle ID list, ignoring invalid targets
            index_dict[component] = filter!(!isnothing, indexin(target_ids, parent_ids))

        end

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
    index_dict = parentIDToIndex(data_dict, parent_ids)

    return Dict(component => data_dict[component]["MASS"][idx] for (component, idxs) in index_dict)

end

##################
# Tracer tracking
##################

"""
    findTracers(data_dict::Dict; <keyword arguments>)::Vector{UInt}

Find the tracers whose parents are allowed by `filter_function`.

# Arguments

  - `data_dict::Dict`: Data dictionary (see [`makeDataDict`](@ref) for the canonical description).
  - `filter_function::Function=filterNothing`: Filter function. See the required signature and examples in `./src/analysis/filters.jl`.

# Returns

  - A vector with the IDs of the tracers.
"""
function findTracers(data_dict::Dict; filter_function::Function=filterNothing)::Vector{UInt}

    # Find the indices of the cells and particles that are allowed by `filter_function`
    selected_indices = filter_function(data_dict)

    parent_ids = UInt[]

    for component in snapshotTypes(data_dict)

        data = data_dict[component]

        if haskey(data, "ID  ") && !isempty(data["ID  "])

            # Read the IDs of the cells and particles that are allowed by `filter_function`
            append!(parent_ids, data["ID  "][selected_indices[component]])

        end

    end

    !isempty(parent_ids) || return UInt[]

    # Find the IDs of the tracers whose parent are allowed by `filter_function`,
    # ignoring cells and particles with no tracers
    return parentToTracerID(data_dict, parent_ids)

end

"""
    tracersWithinR200(data_dict::Dict; <keyword arguments>)::Vector{UInt}

Find the tracers whose parents are within the virial radius (``R_{200}``) of a given halo.

!!! note

    The center of the target halo is its potential minimum in the simulation reference frame. This will not work if `data_dict` has been translated before calling this function.

# Arguments

  - `data_dict::Dict`: Data dictionary (see [`makeDataDict`](@ref) for the canonical description).
  - `halo_idx::Int=1`: Index of the target halo (FoF group). Starts at 1.

# Returns

  - A vector with the IDs of the tracers.
"""
function tracersWithinR200(data_dict::Dict; halo_idx::Int=1)::Vector{UInt}

    if isempty(data_dict[:group]["G_R_Crit200"])

        (
            logging[] &&
            @warn("tracersWithinR200: There is missing group data in $(data_dict[:gc_data].path), \
            so every tracer will be filtered out")
        )

        filter_function = filterAll

    else

        # Read the virial radius
        r200 = data_dict[:group]["G_R_Crit200"][halo_idx]

        if iszero(r200)

            (
                logging[] &&
                @warn("tracersWithinR200: The virial radius of halo index $(halo_idx) is zero in \
                $(data_dict[:gc_data].path), so every tracer will be filtered out")
            )

            return UInt[]

        end

        # Construct a filter function that only allows cells and particles within the virial radius
        filter_function = dd -> filterBySphere(dd, 0.0u"kpc", r200, halo_idx, 0)

    end

    return findTracers(data_dict; filter_function)

end

"""
    tracersWithinDisc(data_dict::Dict; <keyword arguments>)::Vector{UInt}

Find the tracers whose parents are within a given cylinder centered at the origin.

# Arguments

  - `data_dict::Dict`: Data dictionary (see [`makeDataDict`](@ref) for the canonical description).
  - `max_r::Unitful.Length=DISK_R`: Radius of the cylinder.
  - `max_z::Unitful.Length=DISK_HEIGHT`: Half height of the cylinder.

# Returns

  - A vector with the IDs of the tracers.
"""
function tracersWithinDisc(
    data_dict::Dict;
    max_r::Unitful.Length=DISK_R,
    max_z::Unitful.Length=DISK_HEIGHT,
)::Vector{UInt}

    # Construct a filter function that only allows cells and particles within the given cylinder
    filter_function = dd -> filterByCylinder(dd, max_r, max_z, :zero)

    return findTracers(data_dict; filter_function)

end
