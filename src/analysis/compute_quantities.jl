####################################################################################################
# Computation of derived quantities
####################################################################################################

####################################################################################################
# Positions
####################################################################################################

"""
    computeCenter(data_dict::Dict, subfind_idx::NTuple{2,Int})::Vector{<:Unitful.Length}

Read the position of the particle/cell at the potencial minimum of a given halo or subhalo.

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
      + Index of the target subhalo (subfind), relative the target halo. Starts at 1. If it is set to 0, the potencial minimum of the halo with index `halo_idx` is returned.

# Returns

  - The position of the potencial minimum.
"""
function computeCenter(data_dict::Dict, subfind_idx::NTuple{2,Int})::Vector{<:Unitful.Length}

    # If there are no subfind data, return the origin
    if ismissing(data_dict[:gc_data].path) && !isSubfindActive(data_dict[:gc_data].path)
        return zeros(typeof(1.0u"kpc"), 3)
    end

    halo_idx, subhalo_rel_idx = subfind_idx

    # Load the necessary data
    n_subhalos_in_halo = data_dict[:group]["G_Nsubs"]
    g_pos = data_dict[:group]["G_Pos"]
    s_pos = data_dict[:subhalo]["S_Pos"]

    # Check that the requested halo index is within bounds
    n_halos = data_dict[:gc_data].header.n_groups_total

    (
        !iszero(n_halos) && !any(isempty, [n_subhalos_in_halo, g_pos, s_pos]) ||
        return zeros(typeof(1.0u"kpc"), 3)
    )

    (
        0 < halo_idx <= n_halos ||
        throw(ArgumentError("computeCenter: There is only $(n_halos) FoF goups in \
        $(data_dict[:gc_data].path), so halo_idx = $(halo_idx) is out of bounds"))
    )

    # Select the halo potencial minimum if `subhalo_rel_idx` == 0
    isPositive(subhalo_rel_idx) || return g_pos[:, halo_idx]

    # Check that the requested subhalo index is within bounds
    n_subfinds = n_subhalos_in_halo[halo_idx]
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
        n_subs_floor = sum(n_subhalos_in_halo[1:(halo_idx - 1)]; init=0)
    end

    # Compute the subhalo absolute index
    subhalo_abs_idx = n_subs_floor + subhalo_rel_idx

    # Select the subhalo potencial minimum
    return s_pos[:, subhalo_abs_idx]

end

"""
    computeCenter(data_dict::Dict, subhalo_abs_idx::Int)::Vector{<:Unitful.Length}

Read the position of the potencial minimum of a given subhalo.

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

Compute a characteristic center for the system.

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

      + `:global_cm`   -> Center of mass of the whole system.
      + `:{component}` -> Center of mass of the given component (e.g. :stars, :gas, :halo, etc). It can be any of the keys of [`PARTICLE_INDEX`](@ref).
      + `:zero`        -> Origin.

# Returns

  - The specified center of mass.
"""
function computeCenter(data_dict::Dict, cm_type::Symbol)::Vector{<:Unitful.Length}

    if cm_type == :global_cm

        return computeGlobalCenterOfMass(data_dict)

    elseif cm_type ∈ keys(PARTICLE_INDEX)

        return computeCenterOfMass(data_dict[cm_type]["POS "], data_dict[cm_type]["MASS"])

    elseif cm_type == :zero

        return zeros(typeof(1.0u"kpc"), 3)

    end

    throw(ArgumentError("computeCenter: `cm_type` can only be :global_cm, :zero or one of the keys \
    of `PARTICLE_INDEX` but I got :$(center_type)"))

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

Compute the center of mass of the whole system.

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

    components = snapshotTypes(data_dict)

    filter!(ts -> !isempty(data_dict[ts]["POS "]), components)

    # Concatenate the position and masses of all the cells and particles in the system
    positions = hcat([data_dict[component]["POS "] for component in components]...)
    masses    = vcat([data_dict[component]["MASS"] for component in components]...)

    # Check for missing data
    !any(isempty, [positions, masses]) || return zeros(typeof(1.0u"kpc"), 3)

    @debug("computeGlobalCenterOfMass: The center of mass will be computed using $(components)")

    return computeCenterOfMass(positions, masses)

end

####################################################################################################
# Velocities
####################################################################################################

"""
    computeVcm(data_dict::Dict, subfind_idx::NTuple{2,Int})::Vector{<:Unitful.Velocity}

Read the velocity of the center of mass of a given halo or subhalo.

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
      + Index of the target subhalo (subfind), relative the target halo. Starts at 1. If it is set to 0, the potencial minimum of the halo with index `halo_idx` is returned.

# Returns

  - The specified velocity.
"""
function computeVcm(data_dict::Dict, subfind_idx::NTuple{2,Int})::Vector{<:Unitful.Velocity}

    # If there are no subfind data, return the origin
    if ismissing(data_dict[:gc_data].path) && !isSubfindActive(data_dict[:gc_data].path)
        return zeros(typeof(1.0u"km*s^-1"), 3)
    end

    halo_idx, subhalo_rel_idx = subfind_idx

    # Load the necessary data
    n_subhalos_in_halo = data_dict[:group]["G_Nsubs"]
    g_vel = data_dict[:group]["G_Vel"]
    s_vel = data_dict[:subhalo]["S_Vel"]

    # Check that the requested halo index is within bounds
    n_halos = data_dict[:gc_data].header.n_groups_total

    (
        !iszero(n_halos) && !any(isempty, [n_subhalos_in_halo, g_vel, s_vel]) ||
        return zeros(typeof(1.0u"km*s^-1"), 3)
    )

    (
        0 < halo_idx <= n_halos ||
        throw(ArgumentError("computeCenter: There is only $(n_halos) FoF goups in \
        $(data_dict[:gc_data].path), so halo_idx = $(halo_idx) is out of bounds"))
    )

    # Select the halo velocity if `subhalo_rel_idx` == 0
    isPositive(subhalo_rel_idx) || return g_vel[:, halo_idx]

    # Check that the requested subhalo index is within bounds
    n_subfinds = n_subhalos_in_halo[halo_idx]
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
        n_subs_floor = sum(n_subhalos_in_halo[1:(halo_idx - 1)]; init=0)
    end

    # Compute the subhalo absolute index
    subhalo_abs_idx = n_subs_floor + subhalo_rel_idx

    # Select the subhalo velocity
    return s_vel[:, subhalo_abs_idx]

end

"""
    computeVcm(data_dict::Dict, subhalo_abs_idx::Int)::Vector{<:Unitful.Velocity}

Read the velocity of the center of mass of a given subhalo.

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
        return zeros(typeof(1.0u"km*s^-1"), 3)
    end

    s_vel = data_dict[:subhalo]["S_Vel"]

    # Check that the requested subhalo index is within bounds
    n_subgroups_total = data_dict[:gc_data].header.n_subgroups_total

    !iszero(n_subgroups_total) && !isempty(s_vel) || return zeros(typeof(1.0u"km*s^-1"), 3)

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

Compute the velocity of a characteristic center of the system.

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

      + `:global_cm`   -> Center of mass of the whole system.
      + `:{component}` -> Center of mass of the given component (e.g. :stars, :gas, :halo, etc). It can be any of the keys of [`PARTICLE_INDEX`](@ref).
      + `:zero`        -> Origin.

# Returns

  - The specified velocity.
"""
function computeVcm(data_dict::Dict, cm_type::Symbol)::Vector{<:Unitful.Velocity}

    if cm_type == :global_cm

        return computeGlobalVcm(data_dict)

    elseif cm_type ∈ keys(PARTICLE_INDEX)

        return computeComponentVcm(data_dict, cm_type)

    elseif cm_type == :zero

        return zeros(typeof(1.0u"km*s^-1"), 3)

    end

    throw(ArgumentError("computeVcm: `cm_type` can only be :global_cm, :zero or one of the keys of \
    `PARTICLE_INDEX` but I got :$(center_type)"))

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

Compute the rotation matrix that will turn the pricipal axis into the new coordinate system; when view as an passive (alias) trasformation.

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

Compute the rotation matrix that will turn the total angular momentum of the whole system, into the z axis; when view as an active (alibi) trasformation.

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

    components = snapshotTypes(data_dict)

    filter!(ts -> !isempty(data_dict[ts]["POS "]), components)

    # Concatenate the positions, velocities, and masses of all the cells and particles in the system
    positions  = hcat([data_dict[component]["POS "] for component in components]...)
    velocities = hcat([data_dict[component]["VEL "] for component in components]...)
    masses     = vcat([data_dict[component]["MASS"] for component in components]...)

    # Check for missing data
    !any(isempty, [positions, velocities, masses]) || return I

    @debug("computeGlobalAMRotationMatrix: The rotation matrix will be computed \
    using $(components)")

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

Compute the total angular momentum with respect to the origin of the whole system.

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

    components = snapshotTypes(data_dict)

    filter!(ts -> !isempty(data_dict[ts]["POS "]), components)

    # Concatenate the position, velocities, and masses of all the cells and particles in the system
    positions  = hcat([data_dict[component]["POS "] for component in components]...)
    velocities = hcat([data_dict[component]["VEL "] for component in components]...)
    masses     = vcat([data_dict[component]["MASS"] for component in components]...)

    # Check for missing data
    !any(isempty, [positions, velocities, masses]) || return [0.0, 0.0, 1.0]

    @debug("computeGlobalAngularMomentum: The angular momentum will be computed \
    using $(components)")

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
  - `R::Unitful.Length=DISK_R`: Radius.

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
    R::Unitful.Length=DISK_R,
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

Compute the spin parameter of the whole system.

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
  - `R::Unitful.Length=DISK_R`: Radius.

# Returns

  - The spin parameter.
"""
function computeGlobalSpinParameter(data_dict::Dict; R::Unitful.Length=DISK_R)::Float64

    components = snapshotTypes(data_dict)

    filter!(ts -> !isempty(data_dict[ts]["POS "]), components)

    # Concatenate the position and masses of all the cells and particles in the system
    positions  = hcat([data_dict[component]["POS "] for component in components]...)
    velocities = hcat([data_dict[component]["VEL "] for component in components]...)
    masses     = vcat([data_dict[component]["MASS"] for component in components]...)

    @debug("computeGlobalSpinParameter: The spin parameter will be computed using $(components)")

    # Compute the total spin parameter
    return computeSpinParameter(positions, velocities, masses; R)

end

"""
    computeComponentVcm(data_dict::Dict, component::Symbol)::Vector{<:Unitful.Velocity}

Compute the velocity of the given component center of mass.

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
  - `component::Symbol`: Target component, it can be any of the keys of [`PARTICLE_INDEX`](@ref).

# Returns

  - The velocity of the center of mass.
"""
function computeComponentVcm(data_dict::Dict, component::Symbol)::Vector{<:Unitful.Velocity}

    # Load the necessary data
    velocities = data_dict[component]["VEL "]
    masses     = data_dict[component]["MASS"]

    # Check for missing data
    !any(isempty, [velocities, masses]) || return zeros(typeof(1.0u"km*s^-1"), 3)

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

    components = snapshotTypes(data_dict)

    filter!(ts -> !isempty(data_dict[ts]["VEL "]), components)

    # Load the necessary data
    velocities = hcat([data_dict[component]["VEL "] for component in components]...)
    masses     = vcat([data_dict[component]["MASS"] for component in components]...)

    # Check for missing data
    !any(isempty, [velocities, masses]) || return zeros(typeof(1.0u"km*s^-1"), 3)

    # Compute the total mass
    M = sum(masses)

    # Compute the velocity of the center of mass
    vcm = [sum(row .* masses) / M for row in eachrow(velocities)]

    return vcm

end

@doc raw"""
    computeVcirc(
        data_dict::Dict;
        <keyword arguments>
    )::Tuple{Vector{<:Unitful.Length},Vector{<:Unitful.Velocity}}

Compute the circular velocity of each particle of the given type, with respect to the origin.

The circular velocity of a particle is,

```math
v_\mathrm{circ} = \sqrt{\frac{\mathrm{G} \, M(r)}{r}} \, ,
```

where $r$ is the radial distance of the particle, and $M(r)$ is the total mass within a sphere of radius $r$.

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
  - `type::Symbol=:stars`: Target component.

# Returns

  - A tuple with two elements:

      + A vector with the radial distance of each particle to the origin.
      + A vector with the circular velocity of each particle.
"""
function computeVcirc(
    data_dict::Dict;
    type::Symbol=:stars,
)::Tuple{Vector{<:Unitful.Length},Vector{<:Unitful.Velocity}}

    # Compute the radial distance to each particle
    rs = computeDistance(data_dict[type]["POS "])

    # Check for missing data
    !isempty(rs) || return rs, Unitful.Velocity[]

    components = filter!(
        ts -> !isempty(data_dict[ts]["POS "]),
        [:stars, :gas, :halo, :black_hole],
    )

    # Concatenate the position and masses of all the cells and particles in the system
    distances = vcat([computeDistance(data_dict[component]["POS "]) for component in components]...)
    masses    = vcat([data_dict[component]["MASS"] for component in components]...)

    # Use the radial distances as bin edges for the mass histogram
    edges = [0.0u"kpc", rs...]

    # Compute to total mass within each particle radial distance
    M = similar(rs, eltype(masses))
    cumsum!(M, histogram1D(distances, masses, edges; empty_nan=false))

    # The mass histogram is a sorted array, so it is reverted to the unsorted order of `r`
    # to make `vcirc` the circular velocity of each particle in the order of the snapshot
    invpermute!(M, sortperm(rs))

    @debug("computeVcirc: The circular velocity will be computed using $(components)")

    vcirc = [iszero(r) ? 0.0u"km*s^-1" : sqrt(Unitful.G * m / r) for (m, r) in zip(M, rs)]

    return rs, vcirc

end

@doc raw"""
    computeCircularity(data_dict::Dict; <keyword arguments>)::Vector{Float64}

Compute the circularity of each particle of the given type, with respect to the origin and the $z$ direction [0, 0, 1].

The circularity of a particle is,

```math
\epsilon = j_z / j_\mathrm{circ} \, ,
```

where $j_z$ is the $z$ component of its specific angular momentum, and $j_\mathrm{circ}$ is the specific angular momentum of a circular orbit,

```math
j_\mathrm{circ} = r \, v_\mathrm{circ} = \sqrt{\mathrm{G} \, r \, M(r)} \, ,
```

where $r$ is the radial distance of the particle, and $M(r)$ is the total mass within a sphere of radius $r$.

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
  - `type::Symbol=:stars`: Target component.

# Returns

  - The circularity $\epsilon$ of each particle.
"""
function computeCircularity(data_dict::Dict; type::Symbol=:stars)::Vector{Float64}

    # Load the necessary data
    positions  = data_dict[type]["POS "]
    velocities = data_dict[type]["VEL "]

    # Check for missing data
    !any(isempty, [positions, velocities]) || return Float64[]

    # Compute the specific angular momentum in the z direction
    jzs = [x[1] * v[2] - x[2] * v[1] for (x, v) in zip(eachcol(positions), eachcol(velocities))]

    # Compute the circular velocities and the radial distances
    rs, vcircs = computeVcirc(data_dict; type)

    stellar_circularity = [
        any(iszero, [r, vcirc]) ? 0.0 : ustrip(Unitful.NoUnits, jz / (r * vcirc)) for
        (jz, r, vcirc) in zip(jzs, rs, vcircs)
    ]

    return stellar_circularity

end

@doc raw"""
    computeVpolar(
        data_dict::Dict,
        component::Symbol;
        <keyword arguments>
    )::Vector{<:Unitful.Velocity}

Compute the cylindrical components of the velocity, $\mathbf{\vec{v}} = v_r \, \mathbf{e_r} + v_\theta \, \mathbf{e_\theta} + v_z \, \mathbf{e_z}$.

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
v^*_z = v_z \, \mathrm{sign}(z) \, ,
```

in order to distinguish between inflows ($v^*_z < 0$) and outflows ($v^*_z > 0$).

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

      + `:radial`     -> Stellar radial speed ($v_r$).
      + `:tangential` -> Stellar tangential speed ($v_\theta$).
      + `:zstar`      -> Stellar speed in the z direction, computed as $v_z \, \mathrm{sign}(z)$.
  - `type::Symbol=:stars`: Target cell/particle type.

# Returns

  - The chosen cylindricall component of the velocity.
"""
function computeVpolar(
    data_dict::Dict,
    component::Symbol;
    type::Symbol=:stars,
)::Vector{<:Unitful.Velocity}

    # Load the necessary data
    positions = data_dict[type]["POS "]
    velocities = data_dict[type]["VEL "]

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

        throw(ArgumentError("computeVpolar: `component` can only be :radial, :tangential \
        or :zstar, but I got :$(component)"))

    end

    return vp

end

####################################################################################################
# Masses
####################################################################################################

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

  - The radius containing `percet`% of the total mass.
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

  - The height containing `percet`% of the total mass.
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
    computeMassQty(
        quantity::Vector{<:Number},
        masses::Vector{<:Unitful.Mass};
        <keyword arguments>
    )::Number

Compute the maximum value of `quantity` that "contains" `percet`% of the total mass.

# Arguments

  - `quantity::Vector{<:Number}`: Target quantity.
  - `masses::Vector{<:Unitful.Mass}`: Masses of the cells/particles.
  - `percent::Float64=90.0`: Target percentage of the total mass.

# Returns

  - The maximum value of `quantity` that "contains" `percet`% of the total mass.
"""
function computeMassQty(
    quantity::Vector{<:Number},
    masses::Vector{<:Unitful.Mass};
    percent::Float64=90.0,
)::Number

    (
        0 < percent <=100  ||
        throw(ArgumentError("computeMassQty: The argument `percent` must be between 0 and 100, \
        but I got $(percent)"))
    )

    # Check for missing data
    !any(isempty, [quantity, masses]) || return zero(eltype(quantity))

    # Compute the mass limit
    mass_limit = sum(masses) * (percent / 100.0)

    sort_idxs = sortperm(quantity)

    # Find the mass radius
    accu_mass = 0.0u"Msun"
    target_idx = 0
    for mass in masses[sort_idxs]
        accu_mass += mass
        accu_mass < mass_limit || break
        target_idx += 1
    end

    return quantity[sort_idxs[target_idx]]

end

"""
    computeMassPercent(
        quantity::Vector{<:Number},
        masses::Vector{<:Unitful.Mass},
        qty_limits::Tuple{<:Number,<:Number},
    )::Float64

Compute the fraction of the total mass "contained" within a given values of `quantity`.

# Arguments

  - `quantity::Vector{<:Number}`: Target quantity.
  - `masses::Vector{<:Unitful.Mass}`: Masses of the cells/particles.
  - `qty_limits::Tuple{<:Number,<:Number}`: Limits of the target quantity.

# Returns

  - The fraction of the total mass "contained" within a given value of `quantity`.
"""
function computeMassFraction(
    quantity::Vector{<:Number},
    masses::Vector{<:Unitful.Mass},
    qty_limits::Tuple{<:Number,<:Number},
)::Float64

    # Check for missing data
    !any(isempty, [quantity, masses]) || return 0.0

    # Find the indices of all the cells/particles with `qty_limits[1]` <= `quantity` <= `qty_limits[2]`
    idxs = map(x -> qty_limits[1] <= x <= qty_limits[2], quantity)

    return uconvert(Unitful.NoUnits, sum(masses[idxs]; init=0.0u"Msun") / sum(masses))

end

"""
    computeMetalMass(data_dict::Dict, component::Symbol)::Vector{<:Unitful.Mass}

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
  - `component::Symbol`: For which cell/particle type the metal mass will be calculated. The possibilities are `:stars` and `:gas`.

# Returns

  - The total metal mass in each cell/particle.
"""
function computeMetalMass(data_dict::Dict, component::Symbol)::Vector{<:Unitful.Mass}

    if CODEBASE == :arepo

        if component == :gas
            metals = setPositive(data_dict[:gas]["GZ  "]) .* data_dict[:gas]["MASS"]
        elseif component == :stars
            metals = setPositive(data_dict[:stars]["GZ2 "]) .* data_dict[:stars]["MASS"]
        else
            throw(ArgumentError("computeMetalMass: `component` can only be :stars or :gas, \
            but I got :$(component)"))
        end

    elseif CODEBASE == :opengadget3

        if component == :gas
            metals = sum(setPositive(data_dict[:gas]["GMET"][METAL_LIST, :]); dims=1)
        elseif component == :stars
            metals = sum(setPositive(data_dict[:stars]["GME2"][METAL_LIST, :]); dims=1)
        else
            throw(ArgumentError("computeMetalMass: `component` can only be :stars or :gas, \
            but I got :$(component)"))
        end

    else

        throw(ArgumentError("computeMetalMass: I don't recognize the codebase :$(CODEBASE)"))

    end

    return metals

end

"""
    computeElementMass(
        data_dict::Dict,
        component::Symbol,
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
  - `component::Symbol`: For which cell/particle type the element mass will be calculated. The possibilities are `:stars` and `:gas`.
  - `element::Symbol`: Target element. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).

# Returns

  - The total mass of `element` in each cell/particle.
"""
function computeElementMass(
    data_dict::Dict,
    component::Symbol,
    element::Symbol,
)::Vector{<:Unitful.Mass}

    (
        element ∈ keys(ELEMENT_INDEX) ||
        throw(ArgumentError("computeElementMass: :$(element) is not a tracked element, \
        the options are the keys of `ELEMENT_INDEX`, see `./src/constants/globals.jl`"))
    )

    if component == :gas
        block = "GMET"
    elseif component == :stars
        block = "GME2"
    else
        throw(ArgumentError("computeElementMass: `component` can only be :stars or :gas, \
        but I got :$(component)"))
    end

    values = data_dict[component]

    if any(isempty, [values[block], values["MASS"]])
        return Unitful.Mass[]
    end

    if CODEBASE == :arepo
        masses = setPositive(values[block][ELEMENT_INDEX[element], :]) .* values["MASS"]
    elseif CODEBASE == :opengadget3
        masses = setPositive(values[block][ELEMENT_INDEX[element], :])
    else
        throw(ArgumentError("computeElementMass: I don't recognize the codebase :$(CODEBASE)"))
    end

    return masses

end

@doc raw"""
    computeGlobalAbundance(
        data_dict::Dict,
        component::Symbol,
        element::Symbol;
        <keyword arguments>
    )::Float64

Compute the total abundance of a given element, as $n_X / n_H$ where $n_X$ is the number of atoms of element $X$ and $n_H$ the number of hydrogen atoms.

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
  - `component::Symbol`: For which cell/particle type the abundance will be calculated. The possibilities are `:stars` and `:gas`.
  - `element::Symbol`: Target element. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
  - `solar::Bool=false`: If the result will be normalized to the solar abundance or not.

# Returns

  - The total abundance of `element`.
"""
function computeGlobalAbundance(
    data_dict::Dict,
    component::Symbol,
    element::Symbol;
    solar::Bool=false,
)::Float64

    metal_mass = sum(computeElementMass(data_dict, component, element); init=0.0u"Msun")
    hydrogen_mass = sum(computeElementMass(data_dict, component, :H); init=0.0u"Msun")

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
    computeIonizedMass(data_dict::Dict)::Vector{<:Unitful.Mass}

Compute the mass of ionized hydrogen in every gas cell/particle.

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

  - The mass of ionized hydrogen in every gas cell/particle.
"""
function computeIonizedMass(data_dict::Dict)::Vector{<:Unitful.Mass}

    dg = data_dict[:gas]

    !isempty(dg["MASS"]) || return Unitful.Mass[]

    if "FRAC" ∈ keys(dg) && !any(isempty, [dg["FRAC"], dg["NHP "], dg["NH  "]])

        fi = Vector{Float64}(undef, length(dg["MASS"]))

        @inbounds for i in eachindex(fi)

            # Compute how much time has pass since the last time the cell/particle entered the SF routine
            Δt = data_dict[:snap_data].physical_time - dg["CTIM"][i]

            @inbounds if !isnan(dg["FRAC"][1, i]) && Δt < dg["TAUS"][i]

                # Fraction of ionized hydrogen according to our SF model
                fi[i] = dg["FRAC"][1, i]

            else

                # When there is no data from our model, use the fraction of ionized hydrogen from "NHP "
                fi[i] = dg["NHP "][i] / (dg["NHP "][i] + dg["NH  "][i])

            end

        end

    elseif !isempty(dg["COLM"])

        return dg["MASS"] .- dg["COLM"]

    elseif !any(isempty, [dg["NHP "], dg["NH  "]])

        fi = dg["NHP "][i] ./ (dg["NHP "][i] .+ dg["NH  "][i])

    else

        return Unitful.Mass[]

    end

    return fi .* dg["MASS"]

end

"""
    computeNeutralMass(data_dict::Dict; <keyword arguments>)::Vector{<:Unitful.Mass}

Compute the mass of neutral hydrogen in every gas cell/particle.

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
  - `normalize::Bool=true`: If the output will be normalize to eliminate the stellar fraction of the gas cells. Only relevant for simulation with our SF routine, and for cells that have entered it at least once.

# Returns

  - The mass of neutral hydrogen in every gas cell/particle.
"""
function computeNeutralMass(data_dict::Dict; normalize::Bool=true)::Vector{<:Unitful.Mass}

    dg = data_dict[:gas]

    !isempty(dg["MASS"]) || return Unitful.Mass[]

    if "FRAC" ∈ keys(dg) && !any(isempty, [dg["FRAC"], dg["NHP "], dg["NH  "]])

        fa = Vector{Float64}(undef, length(dg["MASS"]))
        fm = Vector{Float64}(undef, length(dg["MASS"]))

        @inbounds for i in eachindex(fa)

            # Compute how much time has pass since the last time the cell/particle entered the SF routine
            Δt = data_dict[:snap_data].physical_time - dg["CTIM"][i]

            @inbounds if !isnan(dg["FRAC"][2, i]) && Δt < dg["TAUS"][i]

                # Fraction of atomic hydrogen according to our model
                fa[i] = dg["FRAC"][2, i]

                # Fraction of molecular hydrogen according to our model
                @inbounds if normalize
                    fm[i] = dg["FRAC"][3, i] + dg["FRAC"][4, i]
                else
                    fm[i] = dg["FRAC"][3, i]
                end

            else

                # When there is no data from our model, use the fraction of neutral hydrogen from "NH  "
                fa[i] = dg["NH  "][i] / (dg["NHP "][i] + dg["NH  "][i])
                fm[i] = 0.0

            end

        end

        fn = fa .+ fm

    elseif !isempty(dg["COLM"])

        return dg["COLM"]

    elseif !any(isempty, [dg["NHP "], dg["NH  "]])

        fn = dg["NH  "][i] ./ (dg["NHP "][i] .+ dg["NH  "][i])

    else

        return Unitful.Mass[]

    end

    return fn .* dg["MASS"]

end

"""
    computeAtomicMass(data_dict::Dict)::Vector{<:Unitful.Mass}

Compute the mass of atomic hydrogen in every gas cell/particle.

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

  - The mass of atomic hydrogen in every gas cell/particle.
"""
function computeAtomicMass(data_dict::Dict)::Vector{<:Unitful.Mass}

    dg = data_dict[:gas]

    !isempty(dg["MASS"]) || return Unitful.Mass[]

    if "FRAC" ∈ keys(dg) && !any(isempty, [dg["FRAC"], dg["NHP "], dg["NH  "]])

        fa = Vector{Float64}(undef, length(dg["MASS"]))

        @inbounds for i in eachindex(fa)

            # Compute how much time has pass since the last time the cell/particle entered the SF routine
            Δt = data_dict[:snap_data].physical_time - dg["CTIM"][i]

            @inbounds if !isnan(dg["FRAC"][2, i]) && Δt < dg["TAUS"][i]

                # Fraction of atomic hydrogen according to our model
                fa[i] = dg["FRAC"][2, i]

            else

                # When there is no data from our model, use the fraction of neutral hydrogen from "NH  "
                fa[i] = dg["NH  "][i] / (dg["NHP "][i] + dg["NH  "][i])

            end

        end

    elseif !any(isempty, [dg["PRES"], dg["COLM"]])

        relative_pressure = uconvert.(Unitful.NoUnits, dg["PRES"] ./ P0).^ALPHA_BLITZ

        # Fraction of molecular hydrogen according to the pressure relation in Blitz et al. (2006)
        fp = 1.0 ./ (1.0 .+ relative_pressure)

        # Fraction of cold gas from Arepo
        fn = uconvert.(Unitful.NoUnits, dg["COLM"] ./ dg["MASS"])

        fa = 1.0 .- fp .* fn

    else

        return Unitful.Mass[]

    end

    return fa .* dg["MASS"]

end

"""
    computeMolecularMass(data_dict::Dict; <keyword arguments>)::Vector{<:Unitful.Mass}

Compute the mass of molecular hydrogen in every gas cell/particle.

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
  - `normalize::Bool=true`: If the output will be normalize to eliminate the stellar fraction of the gas cells. Only relevant for simulation with our SF routine, and for cells that have entered it at least once.

# Returns

  - The mass of molecular hydrogen in every gas cell/particle.
"""
function computeMolecularMass(data_dict::Dict; normalize::Bool=true)::Vector{<:Unitful.Mass}

    dg = data_dict[:gas]

    !isempty(dg["MASS"]) || return Unitful.Mass[]

    if "FRAC" ∈ keys(dg) && !isempty(dg["FRAC"])

        fm = Vector{Float64}(undef, length(dg["MASS"]))

        @inbounds for i in eachindex(fm)

            # Compute how much time has pass since the last time the cell/particle entered the SF routine
            Δt = data_dict[:snap_data].physical_time - dg["CTIM"][i]

            @inbounds if !isnan(dg["FRAC"][3, i]) && Δt < dg["TAUS"][i]

                # Fraction of molecular hydrogen according to our model
                @inbounds if normalize
                    fm[i] = dg["FRAC"][3, i] + dg["FRAC"][4, i]
                else
                    fm[i] = dg["FRAC"][3, i]
                end

            else

                # When there is no data from our model, assume no molecular hydrogen
                fm[i] = 0.0

            end

        end

    elseif !any(isempty, [dg["PRES"], dg["COLM"]])

        relative_pressure = uconvert.(Unitful.NoUnits, dg["PRES"] ./ P0).^ALPHA_BLITZ

        # Fraction of molecular hydrogen according to the pressure relation in Blitz et al. (2006)
        fp = 1.0 ./ (1.0 .+ relative_pressure)

        # Fraction of cold gas from Arepo
        fn = uconvert.(Unitful.NoUnits, dg["COLM"] ./ dg["MASS"])

        fm = fp .* fn

    else

        return Unitful.Mass[]

    end

    return fm .* dg["MASS"]

end

"""
    computePressureMolecularMass(data_dict::Dict)::Vector{<:Unitful.Mass}

Compute the mass molecular hydrogen in every gas cell/particle using the pressure relation from Blitz et al. (2006).

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

  - The mass of molecular hydrogen in every gas cell/particle.

# References

L. Blitz et al. (2006). *The Role of Pressure in GMC Formation II: The H2-Pressure Relation*. The Astrophysical Journal, **650(2)**, 933. [doi:10.1086/505417](https://doi.org/10.1086/505417)
"""
function computePressureMolecularMass(data_dict::Dict)::Vector{<:Unitful.Mass}

    dg = data_dict[:gas]

    !any(isempty, [dg["MASS"], dg["PRES"], dg["COLM"]]) || return Unitful.Mass[]

    relative_pressure = uconvert.(Unitful.NoUnits, dg["PRES"] ./ P0).^ALPHA_BLITZ

    # Fraction of molecular hydrogen according to the pressure relation in Blitz et al. (2006)
    fp = 1.0 ./ (1.0 .+ relative_pressure)

    # Fraction of cold gas from Arepo
    fn = uconvert.(Unitful.NoUnits, dg["COLM"] ./ dg["MASS"])

    return fp .* fn .* dg["MASS"]

end

"""
    computeStellarGasMass(data_dict::Dict)::Vector{<:Unitful.Mass}

Compute the "stellar mass" in every gas cell/particle.

!!! note

    It can be a non 0 value only for simulations with our SF routine.

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

  - The "stellar mass" in every gas cell/particle.
"""
function computeStellarGasMass(data_dict::Dict)::Vector{<:Unitful.Mass}

    dg = data_dict[:gas]

    !isempty(dg["MASS"]) || return Unitful.Mass[]

    if "FRAC" ∈ keys(dg) && !isempty(dg["FRAC"])

        # When there is no data from our model, use a stellar fraction of 0
        fm = replace!(dg["FRAC"][4, :], NaN => 0.0)

    else

        return zeros(typeof(1.0u"Msun"), length(dg["MASS"]))

    end

    return fm .* dg["MASS"]

end

"""
    computeVirialAccretion(
        present_dd::Dict,
        past_dd::Dict;
        <keyword arguments>
    )::NTuple{3,Unitful.Mass}

Compute the inflow, outflow, and net gain of mass for a given halo virial radius (``R_{200}``), between two snapshots.

# Arguments

  - `present_dd::Dict`: A dictionary, for the present snapshot, with the following shape:

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
  - `past_dd::Dict`: A dictionary, for the past snapshot, with the following shape:

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
  - `halo_idx::Int=1`: Index of the target halo (FoF group). Starts at 1.

# Returns

  - A tuple with three elements:

      + The net increase in mass.
      + The inflow mass.
      + The outflow mass.
"""
function computeVirialAccretion(
    present_dd::Dict,
    past_dd::Dict;
    halo_idx::Int=1,
)::NTuple{3,Unitful.Mass}

    # Find the tracers inside R200 in the present snapshot
    present_tracer_ids = tracersWithinR200(present_dd; halo_idx)
    # Find the tracers inside R200 in the past snapshot
    past_tracer_ids    = tracersWithinR200(past_dd; halo_idx)

    # Find the tracers that are inside R200 now, but where outside R200 in the past
    inflow_ids  = setdiff(present_tracer_ids, past_tracer_ids)
    # Find the tracers that were inside R200 in the past, but are now outside R200
    outflow_ids = setdiff(past_tracer_ids, present_tracer_ids)

    # Compute the mass of each tracer in physical units
    tracer_mass = TRACER_MASS * internalUnits("MASS", present_dd[:snap_data].path)

    # Compute the inflow mass
    inflow_mass = length(inflow_ids) * tracer_mass
    # Compute the outflow mass
    ouflow_mass = length(outflow_ids) * tracer_mass

    # Compute the net mass
    net_mass_increase = inflow_mass - ouflow_mass

    return net_mass_increase, inflow_mass, ouflow_mass

end

"""
    computeDiscAccretion(
        present_dd::Dict,
        past_dd::Dict;
        <keyword arguments>
    )::NTuple{3,Unitful.Mass}

Compute the inflow, outflow, and net gain of mass for a given cylinder, between two snapshots.

# Arguments

  - `present_dd::Dict`: A dictionary, for the present snapshot, with the following shape:

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
  - `past_dd::Dict`: A dictionary, for the past snapshot, with the following shape:

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
  - `max_r::Unitful.Length=DISK_R`: Radius of the cylinder.
  - `max_z::Unitful.Length=5.0u"kpc"`: Half height of the cylinder.

# Returns

  - A tuple with three elements:

      + The net increase in mass.
      + The inflow mass.
      + The outflow mass.
"""
function computeDiscAccretion(
    present_dd::Dict,
    past_dd::Dict;
    max_r::Unitful.Length=DISK_R,
    max_z::Unitful.Length=5.0u"kpc",
)::NTuple{3,Unitful.Mass}

    # Find the tracers inside a given cylinder in the present snapshot
    present_tracer_ids = tracersWithinDisc(present_dd; max_r, max_z)
    # Find the tracers inside a given cylinder in the past snapshot
    past_tracer_ids    = tracersWithinDisc(past_dd; max_r, max_z)

    # Find the tracers that are inside a given cylinder now, but where outside in the past
    inflow_ids  = setdiff(present_tracer_ids, past_tracer_ids)
    # Find the tracers that were inside a given cylinder in the past, but are now outside
    outflow_ids = setdiff(past_tracer_ids, present_tracer_ids)

    # Compute the mass of each tracer in physical units
    tracer_mass = TRACER_MASS * internalUnits("MASS", present_dd[:snap_data].path)

    # Compute the inflow mass
    inflow_mass = length(inflow_ids) * tracer_mass
    # Compute the outflow mass
    ouflow_mass = length(outflow_ids) * tracer_mass

    # Compute the net mass
    net_mass_increase = inflow_mass - ouflow_mass

    return net_mass_increase, inflow_mass, ouflow_mass

end

####################################################################################################
# Other
####################################################################################################

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

        # For cosmological simulations, the time field in the Header of the snapshot is the scale factor
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
        internal_energy::Vector{<:SpecificEnergy},
        electron_fraction::Vector{Float32},
    )::Vector{<:Unitful.Temperature}

Compute the gas temperature.

# Arguments

  - `internal_energy::Vector{<:SpecificEnergy}`: Specific internal energy of every gas cell/particle.
  - `electron_fraction::Vector{Float32}`: Number fraction of electrons in every gas cell/particle.

# Returns

  - The temperature of each gas cell/particle.
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
    computeSFR(
        data_dict::Dict;
        <keyword arguments>
    )::Vector{<:Unitful.MassFlow}

Compute the star formation rate of each stellar particle.

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

@doc raw"""
    computeClumpingFactor(density::Vector{<:Number})::Float64

Compute the clumping factor,

```math
C_\rho = \frac{\rangle rho^2 \langle}{\rangle rho \langle^2} \, .
```

# Arguments

  - `density::Vector{<:Number}`: The density of the cells/particles.

# Returns

  - The clumping factor.
"""
function computeClumpingFactor(density::Vector{<:Number})::Float64

    !isempty(density) || return NaN

    μ, var = mean_and_var(density)

    return 1.0 + uconvert(Unitful.NoUnits, var / μ^2)

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

      + `:stellar_mass`              -> Stellar mass.
      + `:gas_mass`                  -> Gas mass.
      + `:hydrogen_mass`             -> Hydrogen mass.
      + `:dm_mass`                   -> Dark matter mass.
      + `:bh_mass`                   -> Black hole mass.
      + `:molecular_mass`            -> Molecular hydrogen (``\\mathrm{H_2}``) mass.
      + `:br_molecular_mass`         -> Molecular hydrogen (``\\mathrm{H_2}``) mass, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_mass`               -> Atomic hydrogen (``\\mathrm{HI}``) mass.
      + `:ionized_mass`              -> Ionized hydrogen (``\\mathrm{HII}``) mass.
      + `:neutral_mass`              -> Neutral hydrogen (``\\mathrm{HI + H_2}``) mass.
      + `:stellar_number`            -> Number of stellar particles.
      + `:gas_number`                -> Number of gas cells.
      + `:dm_number`                 -> Number of dark matter particles.
      + `:bh_number`                 -> Number of black hole particles.
      + `:molecular_fraction`        -> Gas mass fraction of molecular hydrogen.
      + `:br_molecular_fraction`     -> Gas mass fraction of molecular hydrogen, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_fraction`           -> Gas mass fraction of atomic hydrogen.
      + `:ionized_fraction`          -> Gas mass fraction of ionized hydrogen.
      + `:neutral_fraction`          -> Gas mass fraction of neutral hydrogen.
      + `:stellar_area_density`      -> Stellar area mass density, for a radius of `DISK_R`.
      + `:gas_area_density`          -> Gas mass surface density, for a radius of `DISK_R`.
      + `:molecular_area_density`    -> Molecular mass surface density, for a radius of `DISK_R`.
      + `:br_molecular_area_density` -> Molecular mass surface density, for a radius of `DISK_R`, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_area_density`       -> Atomic hydrogen area mass density, for a radius of `DISK_R`.
      + `:ionized_area_density`      -> Ionized hydrogen area mass density, for a radius of `DISK_R`.
      + `:neutral_area_density`      -> Neutral mass surface density, for a radius of `DISK_R`.
      + `:sfr_area_density`          -> Star formation rate area density, for the last `AGE_RESOLUTION` and a radius of `DISK_R`.
      + `:gas_metallicity`           -> Mass fraction of all elements above He in the gas (solar units).
      + `:stellar_metallicity`       -> Mass fraction of all elements above He in the stars (solar units).
      + `:X_gas_abundance`           -> Gas abundance of element ``\\mathrm{X}``, as ``12 + \\log_{10}(\\mathrm{X \\, / \\, H})``. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
      + `:X_stellar_abundance`       -> Stellar abundance of element ``\\mathrm{X}``, as ``12 + \\log_{10}(\\mathrm{X \\, / \\, H})``. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
      + `:stellar_specific_am`       -> Norm of the stellar specific angular momentum.
      + `:gas_specific_am`           -> Norm of the gas specific angular momentum.
      + `:dm_specific_am`            -> Norm of the dark matter specific angular momentum.
      + `:sfr`                       -> The star formation rate.
      + `:ssfr`                      -> The specific star formation rate.
      + `:observational_sfr`         -> The star formation rate of the last `AGE_RESOLUTION`.
      + `:observational_ssfr`        -> The specific star formation rate of the last `AGE_RESOLUTION`.
      + `:scale_factor`              -> Scale factor.
      + `:redshift`                  -> Redshift.
      + `:physical_time`             -> Physical time since the Big Bang.
      + `:lookback_time`             -> Physical time left to reach the last snapshot.

# Returns

  - The velue of `quantity` for the whole system in `data_dict`.
"""
function integrateQty(data_dict::Dict, quantity::Symbol)::Number

    if quantity == :stellar_mass

        integrated_qty = sum(data_dict[:stars]["MASS"]; init=0.0u"Msun")

    elseif quantity == :gas_mass

        integrated_qty = sum(data_dict[:gas]["MASS"]; init=0.0u"Msun")

    elseif quantity == :hydrogen_mass

        integrated_qty = sum(data_dict[:gas]["MASS"]; init=0.0u"Msun") * HYDROGEN_MASSFRAC

    elseif quantity == :dm_mass

        integrated_qty = sum(data_dict[:halo]["MASS"]; init=0.0u"Msun")

    elseif quantity == :bh_mass

        integrated_qty = sum(data_dict[:black_hole]["MASS"]; init=0.0u"Msun")

    elseif quantity == :molecular_mass

        integrated_qty = sum(computeMolecularMass(data_dict); init=0.0u"Msun")

    elseif quantity == :br_molecular_mass

        integrated_qty = sum(computePressureMolecularMass(data_dict); init=0.0u"Msun")

    elseif quantity == :atomic_mass

        integrated_qty = sum(computeAtomicMass(data_dict); init=0.0u"Msun")

    elseif quantity == :ionized_mass

        integrated_qty = sum(computeIonizedMass(data_dict); init=0.0u"Msun")

    elseif quantity == :neutral_mass

        integrated_qty = sum(computeNeutralMass(data_dict); init=0.0u"Msun")

    elseif quantity == :stellar_number

        integrated_qty = length(data_dict[:stars]["MASS"])

    elseif quantity == :gas_number

        integrated_qty = length(data_dict[:gas]["MASS"])

    elseif quantity == :dm_number

        integrated_qty = length(data_dict[:halo]["MASS"])

    elseif quantity == :bh_number

        integrated_qty = length(data_dict[:black_hole]["MASS"])

    elseif quantity == :molecular_fraction

        molecular_mass = sum(computeMolecularMass(data_dict); init=0.0u"Msun")
        gas_mass = sum(data_dict[:gas]["MASS"]; init=0.0u"Msun")

        if iszero(gas_mass)
            integrated_qty = NaN
        else
            integrated_qty = molecular_mass / gas_mass
        end

    elseif quantity == :br_molecular_fraction

        molecular_mass = sum(computePressureMolecularMass(data_dict); init=0.0u"Msun")
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

        integrated_qty = sum(data_dict[:stars]["MASS"]; init=0.0u"Msun") / area(DISK_R)

    elseif quantity == :gas_area_density

        integrated_qty = sum(data_dict[:gas]["MASS"]; init=0.0u"Msun") / area(DISK_R)

    elseif quantity == :molecular_area_density

        integrated_qty = sum(computeMolecularMass(data_dict); init=0.0u"Msun") / area(DISK_R)

    elseif quantity == :br_molecular_area_density

        integrated_qty = sum(computePressureMolecularMass(data_dict); init=0.0u"Msun") / area(DISK_R)

    elseif quantity == :atomic_area_density

        integrated_qty = sum(computeAtomicMass(data_dict); init=0.0u"Msun") / area(DISK_R)

    elseif quantity == :ionized_area_density

        integrated_qty = sum(computeIonizedMass(data_dict); init=0.0u"Msun") / area(DISK_R)

    elseif quantity == :neutral_area_density

        integrated_qty = sum(computeNeutralMass(data_dict); init=0.0u"Msun") / area(DISK_R)

    elseif quantity == :sfr_area_density

        sfr = sum(computeSFR(data_dict; age_resol=AGE_RESOLUTION); init=0.0u"Msun*yr^-1")

        integrated_qty = sfr / area(DISK_R)

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
      + `:molecular_fraction`          -> Gas mass fraction of molecular hydrogen.
      + `:br_molecular_fraction`       -> Gas mass fraction of molecular hydrogen, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_fraction`             -> Gas mass fraction of atomic hydrogen.
      + `:ionized_fraction`            -> Gas mass fraction of ionized hydrogen.
      + `:neutral_fraction`            -> Gas mass fraction of neutral hydrogen.
      + `:molecular_neutral_fraction`  -> Fraction of molecular hydrogen in the neutral gas.
      + `:mol_eq_quotient`             -> Equilibrium quotient for the molecular fraction equation of the SF model.
      + `:ion_eq_quotient`             -> Equilibrium quotient for the ionized fraction equation of the SF model.
      + `:gas_mass_density`            -> Gas mass density.
      + `:hydrogen_mass_density`       -> Hydrogen mass density.
      + `:gas_number_density`          -> Gas number density.
      + `:molecular_number_density`    -> Molecular hydrogen number density.
      + `:br_molecular_number_density` -> Molecular hydrogen number density, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_number_density`       -> Atomic hydrogen number density.
      + `:ionized_number_density`      -> Ionized hydrogen number density.
      + `:neutral_number_density`      -> Neutral hydrogen number density.
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

# Returns

  - The values of `quantity` for every cell/particle.
"""
function scatterQty(data_dict::Dict, quantity::Symbol)::Vector{<:Number}

    if quantity == :stellar_mass

        scatter_qty = data_dict[:stars]["MASS"]

    elseif quantity == :gas_mass

        scatter_qty = data_dict[:gas]["MASS"]

    elseif quantity == :hydrogen_mass

        scatter_qty = data_dict[:gas]["MASS"] .* HYDROGEN_MASSFRAC

    elseif quantity == :dm_mass

        scatter_qty = data_dict[:halo]["MASS"]

    elseif quantity == :bh_mass

        scatter_qty = data_dict[:black_hole]["MASS"]

    elseif quantity == :molecular_mass

        scatter_qty = computeMolecularMass(data_dict)

    elseif quantity == :br_molecular_mass

        scatter_qty = computePressureMolecularMass(data_dict)

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

    elseif quantity == :br_molecular_fraction

        molecular_mass = computePressureMolecularMass(data_dict)
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

    elseif quantity == :mol_eq_quotient

        dg = data_dict[:gas]

        iterator = zip(
            dg["ETAD"],
            dg["FRAC"][2, :],
            dg["FRAC"][3, :],
            dg["FRAC"][4, :],
            τ_star.(dg["RHOC"] .* u"mp"),
            τ_cond.(dg["RHOC"] .* u"mp", dg["PARZ"]),
        )

        # Allocate memory
        scatter_qty = fill(NaN, length(dg["RHOC"]))

        for (i, (ηd, fa, fm, fs, τS, τC)) in enumerate(iterator)

            !(isnan(fa) || iszero(fa) || isone(fs) || iszero(fm)) || continue

            mol_ls = (fa / fm) * (1 - fs)
            mol_rs = uconvert(Unitful.NoUnits, ((ηd + 1) * τC) / τS)

            scatter_qty[i] = log10(mol_ls / mol_rs)

        end

    elseif quantity == :ion_eq_quotient

        dg = data_dict[:gas]

        iterator = zip(
            dg["ETAI"],
            dg["PARR"],
            dg["FRAC"][1, :],
            dg["FRAC"][3, :],
            GalaxyInspector.τ_star.(dg["RHOC"] .* u"mp"),
            GalaxyInspector.τ_rec.(dg["RHOC"] .* u"mp"),
        )

        # Allocate memory
        scatter_qty = fill(NaN, length(dg["RHOC"]))

        for (i, (ηi, R, fi, fm, τS, τR)) in enumerate(iterator)

            !(isnan(fi) || iszero(fi) || iszero(fm)) || continue

            ion_ls = (fi * fi) / fm
            ion_rs = uconvert(Unitful.NoUnits, ((ηi + R) * τR) / τS)

            scatter_qty[i] = log10(ion_ls / ion_rs)

        end

    elseif quantity == :gas_mass_density

        scatter_qty = data_dict[:gas]["RHO "]

    elseif quantity == :hydrogen_mass_density

        scatter_qty = data_dict[:gas]["RHO "] .* HYDROGEN_MASSFRAC

    elseif quantity == :gas_number_density

        scatter_qty = data_dict[:gas]["RHO "] ./ Unitful.mp

    elseif quantity == :molecular_number_density

        molecular_mass = computeMolecularMass(data_dict)
        volumes = data_dict[:gas]["MASS"] ./ data_dict[:gas]["RHO "]

        scatter_qty = (molecular_mass ./ volumes) ./ (2 * Unitful.mp)

    elseif quantity == :br_molecular_number_density

        molecular_mass = computePressureMolecularMass(data_dict)
        volumes = data_dict[:gas]["MASS"] ./ data_dict[:gas]["RHO "]

        scatter_qty = (molecular_mass ./ volumes) ./ (2 * Unitful.mp)

    elseif quantity == :atomic_number_density

        atomic_mass = computeAtomicMass(data_dict)
        volumes = data_dict[:gas]["MASS"] ./ data_dict[:gas]["RHO "]

        scatter_qty = (atomic_mass ./ volumes) ./ Unitful.mp

    elseif quantity == :ionized_number_density

        ionized_mass = computeIonizedMass(data_dict)
        volumes = data_dict[:gas]["MASS"] ./ data_dict[:gas]["RHO "]

        scatter_qty = (ionized_mass ./ volumes) ./ Unitful.mp

    elseif quantity == :neutral_number_density

        neutral_mass = computeNeutralMass(data_dict)
        volumes = data_dict[:gas]["MASS"] ./ data_dict[:gas]["RHO "]

        scatter_qty = (neutral_mass ./ volumes) ./ Unitful.mp

    elseif quantity == :gas_metallicity

        if CODEBASE == :arepo

            scatter_qty = setPositive(data_dict[:gas]["GZ  "]) ./ SOLAR_METALLICITY

        elseif CODEBASE == :opengadget3

            metals = sum(setPositive(data_dict[:gas]["GMET"][METAL_LIST, :]); dims=1)
            scatter_qty = (metals ./ data_dict[:gas]["MASS"]) ./ SOLAR_METALLICITY

        else

            throw(ArgumentError("scatterQty: I don't recognize the codebase :$(CODEBASE)"))

        end

    elseif quantity == :stellar_metallicity

        if CODEBASE == :arepo

            scatter_qty = setPositive(data_dict[:stars]["GZ2 "]) ./ SOLAR_METALLICITY

        elseif CODEBASE == :opengadget3

            metals = sum(setPositive(data_dict[:stars]["GME2"][METAL_LIST, :]); dims=1)
            scatter_qty = (metals ./ data_dict[:stars]["MASS"]) ./ SOLAR_METALLICITY

        else

            throw(ArgumentError("scatterQty: I don't recognize the codebase :$(CODEBASE)"))

        end

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

    elseif quantity == :gas_sfr

        scatter_qty = data_dict[:gas]["SFR "]

    elseif quantity == :stellar_circularity

        @debug("scatterQty: The stellar circularity depends on the positions and velocities of all \
        cell/particles. So, after filtering, the result for a given star will change.")

        scatter_qty = computeCircularity(data_dict)

    elseif quantity == :stellar_vcirc

        @debug("scatterQty: The stellar circular velocity depends on the positions and velocities \
        of all cell/particles. So, after filtering, the result for a given star will change.")

        _, scatter_qty = computeVcirc(data_dict)

    elseif quantity == :stellar_vradial

        scatter_qty = computeVpolar(data_dict, :radial)

    elseif quantity == :stellar_vtangential

        scatter_qty = computeVpolar(data_dict, :tangential)

    elseif quantity == :stellar_vzstar

        scatter_qty = computeVpolar(data_dict, :zstar)

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

    elseif quantity == :pressure

        scatter_qty = data_dict[:gas]["PRES"]

    else

        throw(ArgumentError("scatterQty: I don't recognize the quantity :$(quantity)"))

    end

    if isempty(scatter_qty)
        return Number[]
    end

    return scatter_qty

end
