####################################################################################################
# Coordinate transformations
####################################################################################################

###############
# Translations
###############

"""
    function translatePoints!(positions::Matrix{<:Number}, new_origin::Vector{<:Number})::Nothing

Translate a system of points, moving `new_origin` to [0, 0, 0].

# Arguments

  - `positions::Matrix{<:Number}`: Points to be translated. Each column is a point and each row a dimension.
  - `new_origin::Vector{<:Number}`: Target origin.
"""
function translatePoints!(positions::Matrix{<:Number}, new_origin::Vector{<:Number})::Nothing

    all(iszero, new_origin) && return nothing

    for i in axes(positions, 2)
        positions[:, i] .-= new_origin
    end

    return nothing

end

"""
    translateData!(
        data_dict::Dict,
        new_coor::Tuple{Vector{<:Unitful.Length},Vector{<:Unitful.Velocity}},
    )::Nothing

Translate the positions and boost the velocities of the cells/particles in `data_dict`, such that `origin` is in [0, 0, 0] and the velocity of the center of mass is [0, 0, 0].

# Arguments

  - `data_dict::Dict`: Data dictionary (see [`makeDataDict`](@ref) for the canonical description).
  - `new_coor::Tuple{Vector{<:Unitful.Length},Vector{<:Unitful.Velocity}}`: New origin and the velocity of the center of mass.
"""
function translateData!(
    data_dict::Dict,
    new_coor::Tuple{Vector{<:Unitful.Length},Vector{<:Unitful.Velocity}},
)::Nothing

    for component in snapshotTypes(data_dict)

        data = data_dict[component]

        if haskey(data, "POS ") && !isempty(data["POS "])
            translatePoints!(data["POS "], new_coor[1])
        end

        if haskey(data, "VEL ") && !isempty(data["VEL "])
            translatePoints!(data["VEL "], new_coor[2])
        end

    end

    return nothing

end

"""
    translateData!(data_dict::Dict, translation::Union{Symbol,NTuple{2,Int},Int})::Nothing

Translate the positions and velocities of the cells/particles in `data_dict`.

# Arguments

  - `data_dict::Dict`: Data dictionary (see [`makeDataDict`](@ref) for the canonical description).
  - `translation::Union{Symbol,NTuple{2,Int},Int}=:zero`: Type of translation. The options are:

      + `:zero`                       -> No translation is applied.
      + `:global_cm`                  -> Sets the center of mass of the whole system as the new origin.
      + `:{component}`                -> Sets the center of mass of the given component (e.g. :stellar, :gas, :dark_matter, etc) as the new origin. It can be any of the keys of [`PARTICLE_INDEX`](@ref).
      + `(halo_idx, subhalo_rel_idx)` -> Sets the position of the potential minimum for the `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo) as the new origin.
      + `(halo_idx, 0)`               -> Sets the center of mass of the `halo_idx::Int` halo as the new origin.
      + `subhalo_abs_idx`             -> Sets the center of mass of the `subhalo_abs_idx::Int` as the new origin.
"""
function translateData!(data_dict::Dict, translation::Union{Symbol,NTuple{2,Int},Int})::Nothing

    origin = computeCenter(data_dict, translation)
    vcm = computeVcm(data_dict, translation)

    translateData!(data_dict, (origin, vcm))

    return nothing

end

############
# Rotations
############

"""
    computeAMRotationMatrix(
        positions::Matrix{<:Unitful.Length},
        velocities::Matrix{<:Unitful.Velocity},
        masses::Vector{<:Unitful.Mass},
    )::Union{Matrix{Float64},UniformScaling{Bool}}

Compute the rotation matrix that will turn the total angular momentum into the z axis; when view as an active (alibi) transformation.

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
    size(positions, 2) < 2 && return I

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

Compute the rotation matrix that will turn the principal axis into the new coordinate system; when view as an passive (alias) transformation.

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
    size(positions, 2) < 2 && return I

    # Center the data (subtract the mean of each row)
    mean_point = mean(positions, dims=2)
    centered_pos = positions .- mean_point

    # Compute the covariance matrix (principal axis operator)
    R = ustrip.(cov(centered_pos; dims=2))

    # Reverse the order of the eigenvectors, making the last column the eigenvector
    # with the largest eigenvalue, which should correspond to the new z axis
    pa = eigvecs(R)[:, end:-1:1]

    # Compute the total angular momentum
    L = computeTotalAngularMomentum(positions, velocities, masses; normal=true)

    # 3rd principal axis ≡ new z axis
    pa_z = pa[:, 3]

    # Rotate the principal axis as to align the third component with the angular momentum
    θ = acos(L ⋅ pa_z)
    n = cross(pa_z, L)
    aligned_pa = AngleAxis(θ, n...) * pa

    # The rotation matrix is made from the principal axis as rows
    rotation_matrix = aligned_pa'

    if det(rotation_matrix) < 0.0
        # If the determinant is < 0, that means that the chosen principal axis for the x and y
        # directions form a left-handed Cartesian reference system (x × y = -z). When applying
        # this as a rotation, the z axis will be flipped. So, in this case we swap the x and y
        # axis to get a right-handed Cartesian reference system (x × y = z) and generate the
        # correct rotation
        rotation_matrix[[1, 2], :] = rotation_matrix[[2, 1], :]
    end

    return Matrix{Float64}(rotation_matrix)

end

"""
    computeGlobalAMRotationMatrix(data_dict::Dict)::Union{Matrix{Float64},UniformScaling{Bool}}

Compute the rotation matrix that will turn the total angular momentum of the whole system, into the z axis; when view as an active (alibi) transformation.

# Arguments

  - `data_dict::Dict`: Data dictionary (see [`makeDataDict`](@ref) for the canonical description).

# Returns

  - The rotation matrix.
"""
function computeGlobalAMRotationMatrix(data_dict::Dict)::Union{Matrix{Float64},UniformScaling{Bool}}

    components = snapshotTypes(data_dict)

    filter!(ts -> !isempty(data_dict[ts]["POS "]), components)

    # Concatenate the positions, velocities, and masses of all the cells and particles in the system
    positions  = hcat([data_dict[component]["POS "] for component in components]...)
    velocities = hcat([data_dict[component]["VEL "] for component in components]...)
    masses     = vcat([data_dict[component]["MASS"] for component in components]...)

    # Check for missing data
    size(positions, 2) < 2 && return I

    (
        !logging[] ||
        @info("computeGlobalAMRotationMatrix: The rotation matrix will be computed using \
        $(components)")
    )

    return computeAMRotationMatrix(positions, velocities, masses)

end

"""
    computeGlobalPARotationMatrix(data_dict::Dict)::Union{Matrix{Float64},UniformScaling{Bool}}

Compute the rotation matrix that will turn the principal axis of the whole system into the new coordinate system; when view as an passive (alias) transformation.

# Arguments

  - `data_dict::Dict`: Data dictionary (see [`makeDataDict`](@ref) for the canonical description).

# Returns

  - The rotation matrix.
"""
function computeGlobalPARotationMatrix(data_dict::Dict)::Union{Matrix{Float64},UniformScaling{Bool}}

    components = snapshotTypes(data_dict)

    filter!(ts -> !isempty(data_dict[ts]["POS "]), components)

    # Concatenate the positions, velocities, and masses of all the cells and particles in the system
    positions  = hcat([data_dict[component]["POS "] for component in components]...)
    velocities = hcat([data_dict[component]["VEL "] for component in components]...)
    masses     = vcat([data_dict[component]["MASS"] for component in components]...)

    # Check for missing data
    size(positions, 2) < 2 && return I

    (
        !logging[] ||
        @info("computeGlobalPARotationMatrix: The rotation matrix will be computed using \
        $(components)")
    )

    return computePARotationMatrix(positions, velocities, masses)

end

"""
    rotateData!(data_dict::Dict, axis_type::Symbol)::Nothing

Rotate the positions and velocities of the cells/particles in `data_dict`.

# Arguments

  - `data_dict::Dict`: Data dictionary (see [`makeDataDict`](@ref) for the canonical description).
  - `rotation_matrix::Union{Matrix{Float64},UniformScaling{Bool}}`: Rotation matrix.
"""
function rotateData!(
    data_dict::Dict,
    rotation_matrix::Union{Matrix{Float64},UniformScaling{Bool}},
)::Nothing

    for component in snapshotTypes(data_dict)

        data = data_dict[component]

        if haskey(data, "POS ") && !isempty(data["POS "])
            data["POS "] = rotation_matrix * data["POS "]
        end

        if haskey(data, "VEL ") && !isempty(data["VEL "])
            data["VEL "] = rotation_matrix * data["VEL "]
        end

    end

    return nothing

end

"""
    computeRotation(
        data_dict::Dict,
        rotation::Symbol,
    )::Union{Matrix{Float64},UniformScaling{Bool}}

Compute the rotation matrix corresponding with `rotation`.

# Arguments

  - `data_dict::Dict`: Data dictionary (see [`makeDataDict`](@ref) for the canonical description).
  - `rotation::Symbol`: Type of rotation. The options are:

      + `:zero`               -> No rotation is applied.
      + `:global_am`          -> Sets the angular momentum of the whole system as the new z axis.
      + `:stellar_am`         -> Sets the stellar angular momentum as the new z axis.
      + `:stellar_pa`         -> Sets the stellar principal axis as the new coordinate system.
      + `:stellar_subhalo_pa` -> Sets the principal axis of the stars in the main subhalo as the new coordinate system.

# Returns

  - The rotation matrix.
"""
function computeRotation(
    data_dict::Dict,
    rotation::Symbol,
)::Union{Matrix{Float64},UniformScaling{Bool}}

    rotation === :zero && return I

    # Choose how to compute the 3×3 rotation matrix
    rotation_matrix = begin

        if rotation === :global_am

            computeGlobalAMRotationMatrix(data_dict)

        elseif rotation === :global_pa

            computeGlobalPARotationMatrix(data_dict)

        elseif rotation === :stellar_am

            isempty(data_dict[:stellar]["MASS"]) && return I

            computeAMRotationMatrix(
                data_dict[:stellar]["POS "],
                data_dict[:stellar]["VEL "],
                data_dict[:stellar]["MASS"],
            )

        elseif rotation === :stellar_pa

            isempty(data_dict[:stellar]["MASS"]) && return I

            computePARotationMatrix(
                data_dict[:stellar]["POS "],
                data_dict[:stellar]["VEL "],
                data_dict[:stellar]["MASS"],
            )

        elseif rotation === :stellar_subhalo_pa

            idxs = filterBySubhalo(data_dict; halo_idx=1, subhalo_rel_idx=1)[:stellar]

            isempty(idxs) && return I

            pos_view = data_dict[:stellar]["POS "][:, idxs]
            vel_view = data_dict[:stellar]["VEL "][:, idxs]
            mass_view = data_dict[:stellar]["MASS"][idxs]

            computePARotationMatrix(pos_view, vel_view, mass_view)

        else

            throw(ArgumentError("rotateData!: I don't recognize the rotation :$(rotation)"))

        end

    end

    return rotation_matrix

end

"""
    computeRotation(
        data_dict::Dict,
        rotation::NTuple{2,Int},
    )::Union{Matrix{Float64},UniformScaling{Bool}}

Compute the rotation matrix corresponding with `rotation`.

# Arguments

  - `data_dict::Dict`: Data dictionary (see [`makeDataDict`](@ref) for the canonical description).
  - `rotation::NTuple{2,Int}`: Type of rotation. The options are:

      + `(halo_idx, subhalo_rel_idx)` -> Sets the principal axis of the stars in `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo) as the new coordinate system.
      + `(halo_idx, 0)`               -> Sets the principal axis of the stars in the `halo_idx::Int` halo as the new coordinate system.
"""
function computeRotation(
    data_dict::Dict,
    rotation::NTuple{2,Int},
)::Union{Matrix{Float64},UniformScaling{Bool}}

    idxs = filterBySubhalo(data_dict; halo_idx=rotation[1], subhalo_rel_idx=rotation[2])[:stellar]

    isempty(idxs) && return I

    pos_view = data_dict[:stellar]["POS "][:, idxs]
    vel_view = data_dict[:stellar]["VEL "][:, idxs]
    mass_view = data_dict[:stellar]["MASS"][idxs]

    rotation_matrix = computePARotationMatrix(pos_view, vel_view, mass_view)

    return rotation_matrix

end

"""
    computeRotation(
        data_dict::Dict,
        rotation::Int,
    )::Union{Matrix{Float64},UniformScaling{Bool}}

Compute the rotation matrix corresponding with `rotation`.

# Arguments

  - `data_dict::Dict`: Data dictionary (see [`makeDataDict`](@ref) for the canonical description).
  - `rotation::Int`: Target subhalo absolute index, starting at 1. Sets the principal axis of the stars in the subhalo as the new coordinate system.
"""
function computeRotation(
    data_dict::Dict,
    rotation::Int,
)::Union{Matrix{Float64},UniformScaling{Bool}}

    idxs = filterBySubhalo(data_dict, rotation)[:stellar]

    isempty(idxs) && return I

    pos_view = data_dict[:stellar]["POS "][:, idxs]
    vel_view = data_dict[:stellar]["VEL "][:, idxs]
    mass_view = data_dict[:stellar]["MASS"][idxs]

    rotation_matrix = computePARotationMatrix(pos_view, vel_view, mass_view)

    return rotation_matrix

end

"""
    rotateData!(data_dict::Dict, rotation::Union{Symbol,NTuple{2,Int},Int})::Nothing

Rotate the positions and velocities of the cells/particles in `data_dict`.

# Arguments

  - `data_dict::Dict`: Data dictionary (see [`makeDataDict`](@ref) for the canonical description).
  - `rotation::Symbol`: Type of rotation. The options are:

      + `:zero`               -> No rotation is applied.
      + `:global_am`          -> Sets the angular momentum of the whole system as the new z axis.
      + `:stellar_am`         -> Sets the stellar angular momentum as the new z axis.
      + `:stellar_pa`         -> Sets the stellar principal axis as the new coordinate system.
      + `:stellar_subhalo_pa` -> Sets the principal axis of the stars in the main subhalo as the new coordinate system.
      + `(halo_idx, subhalo_rel_idx)` -> Sets the principal axis of the stars in `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo) as the new coordinate system.
      + `(halo_idx, 0)`               -> Sets the principal axis of the stars in the `halo_idx::Int` halo as the new coordinate system.
      + Target subhalo absolute index, starting at 1. Sets the principal axis of the stars in the subhalo as the new coordinate system.
"""
function rotateData!(data_dict::Dict, rotation::Union{Symbol,NTuple{2,Int},Int})::Nothing

    rotation === :zero && return nothing

    # Choose how to compute the 3×3 rotation matrix
    rotation_matrix = computeRotation(data_dict, rotation)

    rotateData!(data_dict, rotation_matrix)

    return nothing

end
