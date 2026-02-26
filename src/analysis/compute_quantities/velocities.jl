####################################################################################################
# Compute characteristic velocities and momentums
####################################################################################################

#################
# Base functions
#################

"""
    computeVcm(
        velocities::Matrix{<:Unitful.Velocity},
        masses::Vector{<:Unitful.Mass},
    )::Vector{<:Unitful.Velocity}

Compute the velocity of the center of mass.

# Arguments

  - `velocities::Matrix{<:Unitful.Velocity}`: Velocities of the cells/particles. Each column is a cell/particle and each row a dimension.
  - `masses::Vector{<:Unitful.Mass}`: Mass of every cell/particle.

# Returns

  - The velocity of the center of mass of the cells/particles.
"""
function computeVcm(
    velocities::Matrix{<:Unitful.Velocity},
    masses::Vector{<:Unitful.Mass},
)::Vector{<:Unitful.Velocity}

    # Check for missing data
    if any(isempty, [velocities, masses])

        logging[] && @warn("computeVcm: The velocities or masses are empty, so I will return 0s")

        return zeros(typeof(1.0u"km * s^-1"), 3)

    end

    # Compute the total mass
    M = sum(masses)

    vcm = [sum(row .* masses) / M for row in eachrow(velocities)]

    return vcm

end

@doc raw"""
    computeVcirc(
        distances::Vector{<:Unitful.Length},
        masses::Vector{<:Unitful.Mass},
        rs::Vector{<:Unitful.Length},
    )::Tuple{Vector{<:Unitful.Length},Vector{<:Unitful.Velocity}}

Compute the circular velocity.

The circular velocity of a cell/particle is,

```math
v_\mathrm{circ} = \sqrt{\frac{\mathrm{G} \, M(r)}{r}} \, ,
```

where $r$ is the radial distance of the cell/particle, and $M(r)$ is the total mass within a sphere of radius $r$.

# Arguments

  - `distances::Vector{<:Unitful.Length}`: Radial distances of all the cells/particles in the system.
  - `masses::Vector{<:Unitful.Mass}`: Mass of every cell/particle in the system.
  - `rs::Vector{<:Unitful.Length}`: Radial distances of the target cells/particles.

# Returns

  - A tuple with two elements:

      + A vector with the radial distance of the target cells/particles to the origin.
      + A vector with the circular velocity of the target cells/particles.
"""
function computeVcirc(
    distances::Vector{<:Unitful.Length},
    masses::Vector{<:Unitful.Mass},
    rs::Vector{<:Unitful.Length},
)::Tuple{Vector{<:Unitful.Length},Vector{<:Unitful.Velocity}}

    # Check for missing data
    if isempty(rs)
        (
            logging[] &&
            @warn("computeVcirc: The radial distances are empty, so I will return empty arrays")
        )
        return rs, Unitful.Velocity[]
    end

    # Use the radial distances as bin edges for the mass histogram
    edges = [0.0u"kpc", rs...]

    # Compute to total mass within each cell/particle radial distance
    M = similar(rs, eltype(masses))
    cumsum!(M, histogram1D(distances, masses, edges; empty_nan=false))

    # The mass histogram is a sorted array, so it is reverted to the unsorted order of `r` to make
    # `vcirc` the circular velocity of each cell/particle, in the order they have in the snapshot
    invpermute!(M, sortperm(rs))

    vcirc = [iszero(r) ? 0.0u"km * s^-1" : sqrt(Unitful.G * m / r) for (m, r) in zip(M, rs)]

    return rs, vcirc

end

@doc raw"""
    computeVpolar(
        positions::Matrix{<:Unitful.Length},
        velocities::Matrix{<:Unitful.Velocity},
        vel_type::Symbol,
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

  - `positions::Matrix{<:Unitful.Length}`: Positions of the cells/particles. Each column is a cell/particle and each row a dimension.
  - `velocities::Matrix{<:Unitful.Velocity}`: Velocities of the cells/particles. Each column is a cell/particle and each row a dimension.
  - `vel_type::Symbol`: Which velocity will be calculated. The options are:

      + `:radial`     -> Radial speed ($v_r$).
      + `:tangential` -> Tangential speed ($v_\theta$).
      + `:zstar`      -> Speed in the z direction, computed as $v_z \, \mathrm{sign}(z)$.

# Returns

  - The chosen cylindrical velocity of each cell/particle.
"""
function computeVpolar(
    positions::Matrix{<:Unitful.Length},
    velocities::Matrix{<:Unitful.Velocity},
    vel_type::Symbol,
)::Vector{<:Unitful.Velocity}

    x = positions[1, :]
    y = positions[2, :]
    z = positions[3, :]

    vx = velocities[1, :]
    vy = velocities[2, :]
    vz = velocities[3, :]

    if vel_type == :radial

        vp = similar(x, eltype(vx))

        for i in eachindex(vp)
            r = hypot(x[i], y[i])

            if iszero(r)

                vp[i] = hypot(vx[i], vy[i])

            else

                # Compute the radial component
                vp[i] = (x[i] * vx[i] + y[i] * vy[i]) / r

            end
        end

    elseif vel_type == :tangential

        vp = similar(x, eltype(vx))

        for i in eachindex(vp)
            r = hypot(x[i], y[i])

            if iszero(r)

                vp[i] = 0.0u"km * s^-1"

            else

                # Compute the tangential component
                vp[i] = (x[i] * vy[i] - y[i] * vx[i]) / r

            end
        end

    elseif vel_type == :zstar

        # Compute the z component
        vp = @. vz * sign(z)

    else

        throw(ArgumentError("computeVpolar: `vel_type` can only be :radial, :tangential \
        or :zstar, but I got :$(vel_type)"))

    end

    return vp

end

"""
    computeSpecificAngularMomentum(
        positions::Matrix{<:Unitful.Length},
        velocities::Matrix{<:Unitful.Velocity},
    )::Vector{<:Unitful.KinematicViscosity}

Compute the specific angular momentum in the z direction with respect to the origin.

# Arguments

  - `positions::Matrix{<:Unitful.Length}`: Positions of the cells/particles. Each column is a cell/particle and each row a dimension.
  - `velocities::Matrix{<:Unitful.Velocity}`: Velocities of the cells/particles. Each column is a cell/particle and each row a dimension.

# Returns

  - The specific angular momentum in the z direction of each cell/particle.
"""
function computeSpecificAngularMomentum(
    positions::Matrix{<:Unitful.Length},
    velocities::Matrix{<:Unitful.Velocity},
)::Vector{<:Unitful.KinematicViscosity}

    # Check for missing data
    (
        any(isempty, [positions, velocities]) &&
        throw(ArgumentError("computeSpecificAngularMomentum: The angular momentum of an empty \
        dataset is not defined"))
    )

    iterator = zip(eachcol(positions), eachcol(velocities))

    jz = map(((r, v),) -> r[1] * v[2] - r[2] * v[1], iterator)

    return jz

end

"""
    computeAngularMomentum(
        positions::Matrix{<:Unitful.Length},
        velocities::Matrix{<:Unitful.Velocity},
        masses::Vector{<:Unitful.Mass},
    )::Vector{<:AngularMomentum}

Compute the angular momentum in the z direction with respect to the origin.

# Arguments

  - `positions::Matrix{<:Unitful.Length}`: Positions of the cells/particles. Each column is a cell/particle and each row a dimension.
  - `velocities::Matrix{<:Unitful.Velocity}`: Velocities of the cells/particles. Each column is a cell/particle and each row a dimension.
  - `masses::Vector{<:Unitful.Mass}`: Mass of every cell/particle.

# Returns

  - The angular momentum in the z direction of each cell/particle.
"""
function computeAngularMomentum(
    positions::Matrix{<:Unitful.Length},
    velocities::Matrix{<:Unitful.Velocity},
    masses::Vector{<:Unitful.Mass},
)::Vector{<:AngularMomentum}

    # Check for missing data
    (
        any(isempty, [positions, velocities, masses]) &&
        throw(ArgumentError("computeAngularMomentum: The angular momentum of an empty dataset is \
        not defined"))
    )

    iterator = zip(masses, eachcol(positions), eachcol(velocities))

    lz = map(((m, r, v),) -> m * (r[1] * v[2] - r[2] * v[1]), iterator)

    return lz

end

"""
    computeTotalAngularMomentum(
        positions::Matrix{<:Unitful.Length},
        velocities::Matrix{<:Unitful.Velocity},
        masses::Vector{<:Unitful.Mass};
        <keyword arguments>
    )::Vector{<:Number}

Compute the total angular momentum vector with respect to the origin.

# Arguments

  - `positions::Matrix{<:Unitful.Length}`: Positions of the cells/particles. Each column is a cell/particle and each row a dimension.
  - `velocities::Matrix{<:Unitful.Velocity}`: Velocities of the cells/particles. Each column is a cell/particle and each row a dimension.
  - `masses::Vector{<:Unitful.Mass}`: Mass of every cell/particle.
  - `normal::Bool=true`: If the result will be normalized.

# Returns

  - The total angular momentum vector.
"""
function computeTotalAngularMomentum(
    positions::Matrix{<:Unitful.Length},
    velocities::Matrix{<:Unitful.Velocity},
    masses::Vector{<:Unitful.Mass};
    normal::Bool=true,
)::Vector{<:Number}

    # Check for missing data
    if any(isempty, [positions, velocities, masses])
        (
            logging[] &&
            @warn("computeTotalAngularMomentum: The angular momentum of an empty dataset is \
            not defined, so I will return the z axis, [0.0, 0.0, 1.0]")
        )
        return [0.0, 0.0, 1.0]
    end

    unit_L = unit(eltype(masses)) * unit(eltype(positions)) * unit(eltype(velocities))

    Lx = 0.0 * unit_L
    Ly = 0.0 * unit_L
    Lz = 0.0 * unit_L

    # Unroll the cross product manually: L = m * (r × v)
    # Lx = m * (ry * vz - rz * vy)
    # Ly = m * (rz * vx - rx * vz)
    # Lz = m * (rx * vy - ry * vx)
    @inbounds @simd for i in eachindex(masses)
        m = masses[i]

        rx = positions[1, i]
        ry = positions[2, i]
        rz = positions[3, i]

        vx = velocities[1, i]
        vy = velocities[2, i]
        vz = velocities[3, i]

        Lx += m * (ry * vz - rz * vy)
        Ly += m * (rz * vx - rx * vz)
        Lz += m * (rx * vy - ry * vx)
    end

    return normal ? normalize!(ustrip.([Lx, Ly, Lz])) : [Lx, Ly, Lz]

end

@doc raw"""
    computeSpinParameter(
        positions::Matrix{<:Unitful.Length},
        velocities::Matrix{<:Unitful.Velocity},
        masses::Vector{<:Unitful.Mass};
        <keyword arguments>
    )::Float64

Compute the spin parameter with respect to the origin.

The spin parameter was originally defined by Peebles (1969) as,

```math
\lambda = \frac{J \, \sqrt{E}}{G \, M^{5/2}} \, ,
```

where $J$ is the norm of the total angular momentum, $M$ the total mass, $G$ the gravitational constant, and

```math
E = |E_P + E_k| \, ,
```

where $E_P$ is the total potential energy and $E_k$ is the total kinetic energy (including thermal energy of the gas).

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
  - `R::Unitful.Length=DISK_R`: Characteristic radius.

# Returns

  - The spin parameter of the cell/particles within radius `R`.

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
        is not defined"))
    )

    (
        isPositive(R) ||
        throw(ArgumentError("computeSpinParameter: `R` must be greater than 0, but I got $(R)"))
    )

    # Find the cells/particles within `R`
    distances = colwise(Euclidean(), positions, zeros(eltype(positions), size(positions, 1)))
    idx = map(x -> x <= R, distances)

    # Compute the total mass within `R`
    M = sum(masses[idx]; init=0.0u"Msun")

    if iszero(M)
        (
            logging[] &&
            @warn("computeSpinParameter: The total mass within radius $(R) is 0, so \
            the spin parameter will be NaN")
        )
        return NaN
    end

    # Compute the norm of the total angular momentum
    J = norm(
        computeTotalAngularMomentum(
            positions[:, idx],
            velocities[:, idx],
            masses[idx];
            normal=false,
        ),
    )

    λ = uconvert(Unitful.NoUnits, J / sqrt(2.0 * R * Unitful.G * M^3))

    return λ

end

@doc raw"""
    computeCircularity(
        jzs::Vector{<:Unitful.KinematicViscosity},
        rs::Vector{<:Unitful.Length},
        vcircs::Vector{<:Unitful.Velocity},
    )::Vector{Float64}

Compute the circularity with respect to the origin and the $z$ direction.

The circularity of a cell/particle is,

```math
\epsilon = j_z / j_\mathrm{circ} \, ,
```

where $j_z$ is the $z$ component of its specific angular momentum, and $j_\mathrm{circ}$ is the specific angular momentum of a circular orbit,

```math
j_\mathrm{circ} = r \, v_\mathrm{circ} = \sqrt{\mathrm{G} \, r \, M(r)} \, ,
```

where $r$ is the radial distance of the particle, and $M(r)$ is the total mass within a sphere of radius $r$.

# Arguments

  - `jzs::Vector{<:Unitful.KinematicViscosity}`: Specific angular momentum in the z direction of each cell/particle.
  - `rs::Vector{<:Unitful.Length}`: Radial distances of each cell/particle to the origin.
  - `vcircs::Vector{<:Unitful.Velocity}`: Circular velocity of each cell/particle.

# Returns

  - The circularity of each cell/particle.
"""
function computeCircularity(
    jzs::Vector{<:Unitful.KinematicViscosity},
    rs::Vector{<:Unitful.Length},
    vcircs::Vector{<:Unitful.Velocity},
)::Vector{Float64}

    (
        allequal(length, [jzs, rs, vcircs]) ||
        throw(ArgumentError("computeCircularity: The lengths of `jzs`, `rs` and `vcircs` \
        must be equal"))
    )

    ϵ = similar(rs, Float64)

    for (i, (jz, r, vcirc)) in enumerate(zip(jzs, rs, vcircs))

        if iszero(r) || iszero(vcirc)
            ϵ[i] = 0.0
        else
            ϵ[i] = ustrip(Unitful.NoUnits, jz / (r * vcirc))
        end

    end

    return ϵ

end

###################
# Derive functions
###################

@doc raw"""
    computeVcm(
        data_dict::Dict,
        halo_idx::Int,
        subhalo_rel_idx::Int,
    )::Vector{<:Unitful.Velocity}

Return the velocity of the center of mass of a given halo or subhalo.

# Arguments

  - `data_dict::Dict`: Data dictionary (see [`makeDataDict`](@ref) for the canonical description).
    This function requires the following blocks to be present:

      + `:group`   => ["G\_Nsubs", "G\_Vel"]
      + `:subhalo` => ["S\_Vel"]
  - `halo_idx::Int`: Index of the target halo (FoF group). Starts at 1.
  - `subhalo_rel_idx::Int`: Index of the target subhalo (subfind), relative to the target halo. Starts at 1. If set to `0`, the halo potential minimum is returned.

# Returns

  - The velocity of the center of mass.
"""
function computeVcm(
    data_dict::Dict,
    halo_idx::Int,
    subhalo_rel_idx::Int,
)::Vector{<:Unitful.Velocity}

    # If there are no subfind data, return 0s
    if !isSubfindActive(data_dict[:gc_data].path)

        logging[] && @warn("computeVcm: There is no subfind data, so I will return 0s")

        return zeros(typeof(1.0u"km * s^-1"), 3)

    end

    # Load the necessary data
    n_subhalos_in_halo = data_dict[:group]["G_Nsubs"]
    g_vel = data_dict[:group]["G_Vel"]
    s_vel = data_dict[:subhalo]["S_Vel"]

    # Check that the requested halo index is within bounds
    n_halos = data_dict[:gc_data].header.n_groups_total

    if iszero(n_halos) || any(isempty, [n_subhalos_in_halo, g_vel, s_vel])
        (
            logging[] &&
            @warn("computeVcm: There are no halos in $(data_dict[:gc_data].path), \
            so I will return 0s")
        )
        return zeros(typeof(1.0u"km * s^-1"), 3)
    end

    (
        0 < halo_idx <= n_halos ||
        throw(ArgumentError("computeVcm: There is only $(n_halos) FoF groups in \
        $(data_dict[:gc_data].path), so `halo_idx` = $(halo_idx) is out of bounds"))
    )

    # Select the halo velocity if `subhalo_rel_idx` == 0
    isPositive(subhalo_rel_idx) || return g_vel[:, halo_idx]

    # Check that the requested subhalo index is within bounds
    n_subfinds = n_subhalos_in_halo[halo_idx]

    if iszero(n_subfinds)
        (
            logging[] &&
            @warn("computeVcm: There are 0 subhalos in the FoF group $(halo_idx) from \
            $(data_dict[:gc_data].path), so the velocity will be the halo velocity")
        )
        return g_vel[:, halo_idx]
    end

    (
        subhalo_rel_idx <= n_subfinds ||
        throw(ArgumentError("computeVcm: There is only $(n_subfinds) subhalos for the FoF \
        group $(halo_idx) in $(data_dict[:gc_data].path), so `subhalo_rel_idx` = \
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

@doc raw"""
    computeVcm(data_dict::Dict, subhalo_abs_idx::Int)::Vector{<:Unitful.Velocity}

Return the velocity of the center of mass of a given subhalo.

# Arguments

  - `data_dict::Dict`: Data dictionary (see [`makeDataDict`](@ref) for the canonical description).
    This function requires the following blocks to be present:

      + `:subhalo` => ["S\_Vel"]
  - `subhalo_abs_idx::Int`: Absolute index of the target subhalo (subfind). Starts at 1.

# Returns

  - The velocity of the center of mass.
"""
function computeVcm(data_dict::Dict, subhalo_abs_idx::Int)::Vector{<:Unitful.Velocity}

    # If there are no subfind data, return 0s
    if !isSubfindActive(data_dict[:gc_data].path)

        logging[] && @warn("computeVcm: There is no subfind data, so I will return 0s")

        return zeros(typeof(1.0u"km * s^-1"), 3)

    end

    s_vel = data_dict[:subhalo]["S_Vel"]

    # Check that the requested subhalo index is within bounds
    n_subgroups_total = data_dict[:gc_data].header.n_subgroups_total

    if iszero(n_subgroups_total) || isempty(s_vel)
        (
            logging[] &&
            @warn("computeVcm: There are no subhalos in $(data_dict[:gc_data].path), \
            so I will return 0s")
        )
        return zeros(typeof(1.0u"km * s^-1"), 3)
    end

    (
        0 < subhalo_abs_idx <= n_subgroups_total ||
        throw(ArgumentError("computeVcm: There is only $(n_subgroups_total) subhalos in \
        $(data_dict[:gc_data].path), so `subhalo_abs_idx` = $(subhalo_abs_idx) is out of bounds"))
    )

    # Select the subhalo velocity
    return s_vel[:, subhalo_abs_idx]

end

"""
    computeVcm(data_dict::Dict, cm_type::Symbol)::Vector{<:Unitful.Velocity}

Compute a characteristic velocity for the system.

# Arguments

  - `data_dict::Dict`: Data dictionary (see [`makeDataDict`](@ref) for the canonical description).
    This function requires the following blocks to be present, depending on the value of `cm_type`:

      + If `cm_type` == :all:
          * ["VEL ", "MASS"] for every cell/particle type in the snapshot.
      + If haskey(`PARTICLE_INDEX`, `cm_type`):
          * `cm_type` => ["VEL ", "MASS"].
      + If `cm_type` == :zero:
          * No blocks are required
  - `cm_type::Symbol`: It can be:

      + `:all`         -> Velocity of the center of mass of the whole system.
      + `:{component}` -> Velocity of the center of mass of the given component (e.g. :stellar, :gas, :dark_matter, etc). It can be any of the keys of [`PARTICLE_INDEX`](@ref).
      + `:zero`        -> 0 velocity.

# Returns

  - The velocity.
"""
function computeVcm(data_dict::Dict, cm_type::Symbol)::Vector{<:Unitful.Velocity}

    if cm_type == :all

        snap_types = snapshotTypes(data_dict)

        # Remove components with no velocity or mass data
        filter!(st -> !isempty(data_dict[st]["VEL "]), snap_types)
        filter!(st -> !isempty(data_dict[st]["MASS"]), snap_types)

        # Concatenate the velocities and masses of all the cells and particles in the system
        velocities = hcat([data_dict[st]["VEL "] for st in snap_types]...)
        masses     = vcat([data_dict[st]["MASS"] for st in snap_types]...)

        return computeVcm(velocities, masses)

    elseif haskey(PARTICLE_INDEX, cm_type)

        velocities = data_dict[cm_type]["VEL "]
        masses     = data_dict[cm_type]["MASS"]

        if any(isempty, [velocities, masses])
            (
                logging[] &&
                @warn("computeVcm: The velocities or masses are empty, so I will return 0s")
            )
            return zeros(typeof(1.0u"km * s^-1"), 3)
        end

        return computeVcm(velocities, masses)

    elseif cm_type == :zero

        return zeros(typeof(1.0u"km * s^-1"), 3)

    end

    throw(ArgumentError("computeVcm: `cm_type` can only be :all, :zero or one of the keys of \
    `PARTICLE_INDEX` but I got :$(cm_type)"))

end

@doc raw"""
    computeVcirc(
        data_dict::Dict,
        component::Symbol,
    )::Tuple{Vector{<:Unitful.Length},Vector{<:Unitful.Velocity}}

Compute the circular velocity of each cell/particle of the given component, with respect to the origin.

The circular velocity of a cell/particle is,

```math
v_\mathrm{circ} = \sqrt{\frac{\mathrm{G} \, M(r)}{r}} \, ,
```

where $r$ is the radial distance of the cell/particle, and $M(r)$ is the total mass within a sphere of radius $r$.

# Arguments

  - `data_dict::Dict`: Data dictionary (see [`makeDataDict`](@ref) for the canonical description).
    This function requires the following blocks to be present for every cell/particle that you want to be taken into account:

      + `cell/particle type` => ["POS ", "MASS"].
  - `component::Symbol`: Target component. It can only be one of the elements of [`COMPONENTS`](@ref).

# Returns

  - A tuple with two elements:

      + A vector with the radial distance of the target cells/particles to the origin.
      + A vector with the circular velocity of the target cells/particles.
"""
function computeVcirc(
    data_dict::Dict,
    component::Symbol,
)::Tuple{Vector{<:Unitful.Length},Vector{<:Unitful.Velocity}}

    if component ∉ COMPONENTS
        throw(ArgumentError("computeVcirc: `component` can only be one of the elements of \
        `COMPONENTS` (see `./src/constants/globals.jl`), but I got :$(component)"))
    end

    if component ∈ [:stellar, :dark_matter, :gas, :black_hole]
        type = component
    elseif component == :Z_stellar
        type = :stellar
    else
        type = :gas
    end

    # Compute the radial distance to each cell/particle
    positions = data_dict[type]["POS "]
    rs = colwise(Euclidean(), positions, zeros(eltype(positions), size(positions, 1)))

    snap_types = snapshotTypes(data_dict)

    # Remove components with no position or mass data
    filter!(st -> !isempty(data_dict[st]["POS "]), snap_types)
    filter!(st -> !isempty(data_dict[st]["MASS"]), snap_types)

    # Concatenate the distances and masses of all the cells and particles in the system
    distances = Vector{eltype(positions)}()
    masses    = Vector{typeof(1.0u"Msun")}()
    for st in snap_types
        pos_data  = data_dict[st]["POS "]
        mass_data = data_dict[st]["MASS"]

        dists = colwise(Euclidean(), pos_data, zeros(eltype(pos_data), size(pos_data, 1)))

        append!(distances, dists)
        append!(masses, mass_data)
    end

    logging[] && @info("computeVcirc: The circular velocity will be computed using $(snap_types)")

    return computeVcirc(distances, masses, rs)

end

@doc raw"""
    computeVpolar(
        data_dict::Dict,
        component::Symbol,
        vel_type::Symbol,
    )::Vector{<:Unitful.Velocity}

Compute the cylindrical components of the velocity, $\mathbf{\vec{v}} = v_r \, \mathbf{e_r} + v_\theta \, \mathbf{e_\theta} + v_z \, \mathbf{e_z}$, for each cell/particle of the given component.

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

  - `data_dict::Dict`: Data dictionary (see [`makeDataDict`](@ref) for the canonical description).
    This function requires the following blocks to be present for the cell/particle type corresponding to `component`:

      + `cell/particle type` => ["POS ", "VEL "].
  - `vel_type::Symbol`: Which velocity will be calculated. The options are:

      + `:radial`     -> Stellar radial speed ($v_r$).
      + `:tangential` -> Stellar tangential speed ($v_\theta$).
      + `:zstar`      -> Stellar speed in the z direction, computed as $v_z \, \mathrm{sign}(z)$.
  - `component::Symbol`: Target component. It can only be one of the elements of [`COMPONENTS`](@ref).

# Returns

  - The chosen cylindrical velocity.
"""
function computeVpolar(
    data_dict::Dict,
    component::Symbol,
    vel_type::Symbol,
)::Vector{<:Unitful.Velocity}

    if component ∉ COMPONENTS
        throw(ArgumentError("computeVpolar: `component` can only be one of the elements of \
        `COMPONENTS` (see `./src/constants/globals.jl`), but I got :$(component)"))
    end

    if component ∈ [:stellar, :dark_matter, :gas, :black_hole]
        type = component
    elseif component == :Z_stellar
        type = :stellar
    else
        type = :gas
    end

    # Load the necessary data
    positions  = data_dict[type]["POS "]
    velocities = data_dict[type]["VEL "]

    return computeVpolar(positions, velocities, vel_type)

end

"""
    computeSpecificAngularMomentum(
        data_dict::Dict,
        component::Symbol,
    )::Vector{<:Unitful.KinematicViscosity}

Compute the specific angular momentum in the z direction with respect to the origin, for each cell/particle of the given component.

# Arguments

  - `data_dict::Dict`: Data dictionary (see [`makeDataDict`](@ref) for the canonical description).
    This function requires the following blocks to be present for the cell/particle type corresponding to `component`:

      + `cell/particle type` => ["POS ", "VEL "].

  - `component::Symbol`: Target component. It can only be one of the elements of [`COMPONENTS`](@ref).

# Returns

  - The specific angular momentum of each cell/particle.
"""
function computeSpecificAngularMomentum(
    data_dict::Dict,
    component::Symbol,
)::Vector{<:Unitful.KinematicViscosity}

    if component ∉ COMPONENTS
        throw(ArgumentError("computeSpecificAngularMomentum: `component` can only be one of the \
        elements of `COMPONENTS` (see `./src/constants/globals.jl`), but I got :$(component)"))
    end

    if component ∈ [:stellar, :dark_matter, :gas, :black_hole]
        type = component
    elseif component == :Z_stellar
        type = :stellar
    else
        type = :gas
    end

    # Load the necessary data
    positions  = data_dict[type]["POS "]
    velocities = data_dict[type]["VEL "]

    return computeSpecificAngularMomentum(positions, velocities)

end

"""
    computeAngularMomentum(
        data_dict::Dict,
        component::Symbol,
    )::Vector{<:AngularMomentum}

Compute the angular momentum in the z direction with respect to the origin, for each cell/particle of the given component.

# Arguments

  - `data_dict::Dict`: Data dictionary (see [`makeDataDict`](@ref) for the canonical description).
    This function requires the following blocks to be present for the cell/particle type corresponding to `component`:

      + `cell/particle type` => ["POS ", "VEL ", "MASS"].
  - `component::Symbol`: Target component. It can only be one of the elements of [`COMPONENTS`](@ref).

# Returns

  - The angular momentum of each cell/particle.
"""
function computeAngularMomentum(
    data_dict::Dict,
    component::Symbol,
)::Vector{<:AngularMomentum}

    if component ∉ COMPONENTS
        throw(ArgumentError("computeAngularMomentum: `component` can only be one of the elements \
        of `COMPONENTS` (see `./src/constants/globals.jl`), but I got :$(component)"))
    end

    if component ∈ [:stellar, :dark_matter, :gas, :black_hole]
        type = component
    elseif component == :Z_stellar
        type = :stellar
    else
        type = :gas
    end

    # Load the necessary data
    positions  = data_dict[type]["POS "]
    velocities = data_dict[type]["VEL "]
    masses     = data_dict[type]["MASS"]

    return computeAngularMomentum(positions, velocities, masses)

end

"""
    computeGlobalAngularMomentum(data_dict::Dict; <keyword arguments>)::Vector{<:Number}

Compute the total angular momentum of the system, with respect to the origin

# Arguments

  - `data_dict::Dict`: Data dictionary (see [`makeDataDict`](@ref) for the canonical description).
    This function requires the following blocks to be present for every cell/particle that you want to be taken into account:

      + `cell/particle type` => ["POS ", "VEL ", "MASS"].
  - `normal::Bool=true`: If the result will be normalized.

# Returns

  - The angular momentum.
"""
function computeGlobalAngularMomentum(data_dict::Dict; normal::Bool=true)::Vector{<:Number}

    snap_types = snapshotTypes(data_dict)

    # Remove components with no position, velocity or mass data
    filter!(st -> !isempty(data_dict[st]["POS "]), snap_types)
    filter!(st -> !isempty(data_dict[st]["VEL "]), snap_types)
    filter!(st -> !isempty(data_dict[st]["MASS"]), snap_types)

    # Concatenate the position, velocities, and masses of all the cells and particles in the system
    positions  = hcat([data_dict[st]["POS "] for st in snap_types]...)
    velocities = hcat([data_dict[st]["VEL "] for st in snap_types]...)
    masses     = vcat([data_dict[st]["MASS"] for st in snap_types]...)

    (
        logging[] &&
        @info("computeGlobalAngularMomentum: The angular momentum will be computed using \
        $(components)")
    )

    return computeTotalAngularMomentum(positions, velocities, masses; normal)

end

@doc raw"""
    computeSpinParameter(data_dict::Dict, component::Symbol; <keyword arguments>)::Float64

Compute the spin parameter of the given component, with respect to the origin

The spin parameter was originally defined by Peebles (1969) as,

```math
\lambda = \frac{J \, \sqrt{E}}{G \, M^{5/2}} \, ,
```

where $J$ is the norm of the total angular momentum, $M$ the total mass, $G$ the gravitational constant, and

```math
E = |E_P + E_k| \, ,
```

where $E_P$ is the total potential energy and $E_k$ is the total kinetic energy (including thermal energy of the gas).

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

  - `data_dict::Dict`: Data dictionary (see [`makeDataDict`](@ref) for the canonical description).
    This function requires the following blocks to be present for the cell/particle type corresponding to `component`:

      + `cell/particle type` => ["POS ", "VEL ", "MASS"].
  - `component::Symbol`: Target component. It can only be one of the elements of [`COMPONENTS`](@ref).
  - `R::Unitful.Length=DISK_R`: Characteristic radius.

# Returns

  - The spin parameter of the cell/particles within radius `R`.

# References

P. J. E. Peebles (1969). *Origin of the Angular Momentum of Galaxies*. Astrophysical Journal, **155**, 393. [doi:10.1086/149876](https://doi.org/10.1086/149876)

J. S. Bullock et al. (2001). *A Universal Angular Momentum Profile for Galactic Halos*. The Astrophysical Journal, **555(1)**, 240. [doi:10.1086/321477](https://doi.org/10.1086/321477)

J. Zjupa et al. (2017). *Angular momentum properties of haloes and their baryon content in the Illustris simulation*. Monthly Notices of the Royal Astronomical Society, **466(2)**, 1625–1647. [doi:10.1093/mnras/stw2945](https://doi.org/10.1093/mnras/stw2945)
"""
function computeSpinParameter(data_dict::Dict, component::Symbol; R::Unitful.Length=DISK_R)::Float64

    if component ∉ COMPONENTS
        throw(ArgumentError("computeSpinParameter: `component` can only be one of the elements \
        of `COMPONENTS` (see `./src/constants/globals.jl`), but I got :$(component)"))
    end

    if component ∈ [:stellar, :dark_matter, :gas, :black_hole]
        type = component
    elseif component == :Z_stellar
        type = :stellar
    else
        type = :gas
    end

    # Load the necessary data
    positions  = data_dict[type]["POS "]
    velocities = data_dict[type]["VEL "]
    masses     = data_dict[type]["MASS"]

    return computeSpinParameter(positions, velocities, masses; R)

end

@doc raw"""
    computeGlobalSpinParameter(data_dict::Dict; <keyword arguments>)::Float64

Compute the spin parameter of the whole system.

The spin parameter was originally defined by Peebles (1969) as,

```math
\lambda = \frac{J \, \sqrt{E}}{G \, M^{5/2}} \, ,
```

where $J$ is the norm of the total angular momentum, $M$ the total mass, $G$ the gravitational constant, and

```math
E = |E_P + E_k| \, ,
```

where $E_P$ is the total potential energy and $E_k$ is the total kinetic energy (including thermal energy of the gas).

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

  - `data_dict::Dict`: Data dictionary (see [`makeDataDict`](@ref) for the canonical description).
    This function requires the following blocks to be present for every cell/particle that you want to be taken into account:

      + `cell/particle type` => ["POS ", "VEL ", "MASS"].
  - `component::Symbol`: Target component. It can only be one of the elements of [`COMPONENTS`](@ref).
  - `R::Unitful.Length=DISK_R`: Characteristic radius.

# Returns

  - The spin parameter of the cell/particles within radius `R`.

# References

P. J. E. Peebles (1969). *Origin of the Angular Momentum of Galaxies*. Astrophysical Journal, **155**, 393. [doi:10.1086/149876](https://doi.org/10.1086/149876)

J. S. Bullock et al. (2001). *A Universal Angular Momentum Profile for Galactic Halos*. The Astrophysical Journal, **555(1)**, 240. [doi:10.1086/321477](https://doi.org/10.1086/321477)

J. Zjupa et al. (2017). *Angular momentum properties of haloes and their baryon content in the Illustris simulation*. Monthly Notices of the Royal Astronomical Society, **466(2)**, 1625–1647. [doi:10.1093/mnras/stw2945](https://doi.org/10.1093/mnras/stw2945)
"""
function computeGlobalSpinParameter(data_dict::Dict; R::Unitful.Length=DISK_R)::Float64

    snap_types = snapshotTypes(data_dict)

    filter!(ts -> !isempty(data_dict[ts]["POS "]), snap_types)

    # Concatenate the position and masses of all the cells and particles in the system
    positions  = hcat([data_dict[component]["POS "] for component in snap_types]...)
    velocities = hcat([data_dict[component]["VEL "] for component in snap_types]...)
    masses     = vcat([data_dict[component]["MASS"] for component in snap_types]...)

    (
        logging[] &&
        @info("computeGlobalSpinParameter: The spin parameter will be computed using $(snap_types)")
    )

    # Compute the total spin parameter
    return computeSpinParameter(positions, velocities, masses; R)

end

@doc raw"""
    computeCircularity(data_dict::Dict, component::Symbol)::Vector{Float64}

Compute the circularity of each cell/particle, with respect to the origin and the $z$ direction.

The circularity of a cell/particle is,

```math
\epsilon = j_z / j_\mathrm{circ} \, ,
```

where $j_z$ is the $z$ component of its specific angular momentum, and $j_\mathrm{circ}$ is the specific angular momentum of a circular orbit,

```math
j_\mathrm{circ} = r \, v_\mathrm{circ} = \sqrt{\mathrm{G} \, r \, M(r)} \, ,
```

where $r$ is the radial distance of the particle, and $M(r)$ is the total mass within a sphere of radius $r$.

# Arguments

  - `data_dict::Dict`: Data dictionary (see [`makeDataDict`](@ref) for the canonical description).
    This function requires the following blocks to be present for every cell/particle that you want to be taken into account:

      + `cell/particle type` => ["VEL ", "POS ", "MASS"].
  - `component::Symbol`: Target component. It can only be one of the elements of [`COMPONENTS`](@ref).

# Returns

  - The circularity of each cell/particle.
"""
function computeCircularity(data_dict::Dict, component::Symbol)::Vector{Float64}

    # Compute the specific angular momentum in the z direction
    jzs = computeSpecificAngularMomentum(data_dict, component)

    # Compute the circular velocities and the radial distances
    rs, vcircs = computeVcirc(data_dict, component)

    return computeCircularity(jzs, rs, vcircs)

end
