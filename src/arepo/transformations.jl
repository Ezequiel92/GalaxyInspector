####################################################################################################
# Coordinate transformations.
####################################################################################################

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
