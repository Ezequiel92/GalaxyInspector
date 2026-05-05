####################################################################################################
# Type aliases and structs
####################################################################################################

"""
Color type.
"""
const ColorType = Union{ColorTypes.RGB,ColorTypes.RGBA,Symbol}

"""
Line style type.
"""
const LineStyleType = Union{Tuple{String,Symbol},Nothing,String,Symbol}

"""
Translation type
"""
const TranslationType = Union{
    Tuple{Vector{<:Unitful.Length},Vector{<:Unitful.Velocity}},
    Symbol,
    Int,
    NTuple{2,Int},
}

"""
Rotation type
"""
const RotationType = Union{
    Tuple{Symbol,Symbol,Function},
    Tuple{Matrix{Float64}},
    Tuple{UniformScaling{Bool}},
}

"""
Index type.
"""
const IndexType = Union{
    Colon,
    Integer,
    UnitRange{<:Integer},
    StepRange{<:Integer,<:Integer},
    Vector{<:Integer},
    InvertedIndex,
}

"""
Reduced index type.
"""
const ReducedIndexType = Union{
    Integer,
    UnitRange{<:Integer},
    StepRange{<:Integer,<:Integer},
    Vector{<:Integer},
    InvertedIndex,
}

# Dimensions of specific energy
@derived_dimension SpecificEnergy Unitful.𝐋^2 * Unitful.𝐓^-2 true

# Dimensions of surface density
@derived_dimension SurfaceDensity Unitful.𝐌 * Unitful.𝐋^-2 true

# Dimensions of mass flow surface density
@derived_dimension MassFlowDensity Unitful.𝐌 * Unitful.𝐓^-1 * Unitful.𝐋^-2 true

# Dimensions of angular momentum
@derived_dimension AngularMomentum Unitful.𝐌 * Unitful.𝐋^2 * Unitful.𝐓^-1 true

# Dimensions of number density
@derived_dimension NumberDensity Unitful.𝐋^-3 true

"""
Dimensional information about a physical quantity.

# Fields

  - `hdf5_name  :: String`: HDF5 block name.
  - `dimensions :: Unitful.Dimensions`: Physical dimensions of the quantity, e.g. `Unitful.𝐋 * Unitful.𝐓^-1`.
  - `unit       :: Union{Unitful.Units,Symbol}`: Units of the quantity within the simulation code. It can be a unit from [Unitful](https://github.com/PainterQubits/Unitful.jl) or [UnitfulAstro](https://github.com/JuliaAstro/UnitfulAstro.jl), or it can be the symbol `:internal` which denotes internal code units.
"""
struct Qty
    hdf5_name  :: String
    dimensions :: Unitful.Dimensions
    unit       :: Union{Unitful.Units,Symbol}
end

"""
Plotting parameters for a quantity.

# Fields

  - `request::Dict{Symbol,Vector{String}} = Dict{Symbol,Vector{String}}()`: Data request for [`readSnapshot`](@ref). It must have the shape `cell/particle type` -> [`block`, `block`, ...], where the possible types are the keys of [`PARTICLE_INDEX`](@ref), :group, and :subhalo, and the possible blocks are the keys of [`QUANTITIES`](@ref).
  - `var_name::AbstractString = ""`: Name of the quantity for the axis label. It should not include units or scaling factors.
  - `exp_factor::Int = 0`: Numerical exponent to scale down the axis, e.g. if `exp_factor` = 10 the values will be divided by ``10^{10}``. The default is no scaling.
  - `unit::Unitful.Units = Unitful.NoUnits`: Target unit for the axis.
  - `axis_label::AbstractString = "auto_label"`: Label for the axis. It can contain the string `auto_label`, which will be replaced by the default label: `var_name` / 10^`exp_factor` `unit`.
  - `cp_type::Union{Symbol,Nothing} = nothing`: Cell/particle type corresponding to the quantity. It can be any of the keys of [`PARTICLE_INDEX`](@ref).
"""
@kwdef struct PlotParams
    request    :: Dict{Symbol,Vector{String}} = Dict{Symbol,Vector{String}}()
    var_name   :: AbstractString              = ""
    exp_factor :: Int                         = 0
    unit       :: Unitful.Units               = Unitful.NoUnits
    axis_label :: AbstractString              = "auto_label"
    cp_type    :: Union{Symbol,Nothing}       = nothing
end

"""
Data in the "Header" group of a HDF5 snapshot file.

# Fields

  - `box_size   :: Float64`: Total size of the simulation box, in internal units of length.
  - `h          :: Float64`: Dimensionless Hubble parameter, "little h".
  - `mass_table :: Vector{Float64}`: Masses of particle/cell types which have a constant mass, in internal units of mass.
  - `num_files  :: Int32`: Number of file chunks per snapshot.
  - `num_part   :: Vector{Int32}`: Number of particles/cells (of each type) included in this file chunk.
  - `num_total  :: Vector{UInt32}`: Total number of particles/cells (of each type) for this snapshot.
  - `omega_0    :: Float64`: The cosmological density parameter for matter.
  - `omega_l    :: Float64`: The cosmological density parameter for the cosmological constant.
  - `redshift   :: Float64`: The redshift.
  - `time       :: Float64`: The physical time or the scale factor, depending on the type of simulation.
  - `l_unit     :: Unitful.Length`: Conversion factor from internal units of length to centimeters.
  - `m_unit     :: Unitful.Mass`: Conversion factor from internal units of mass to grams.
  - `v_unit     :: Unitful.Velocity`: Conversion factor from internal units of velocity to centimeters per second.
"""
@kwdef struct SnapshotHeader
    box_size   :: Float64
    h          :: Float64
    mass_table :: Vector{Float64}
    num_files  :: Int32
    num_part   :: Vector{Int32}
    num_total  :: Vector{UInt32}
    omega_0    :: Float64
    omega_l    :: Float64
    redshift   :: Float64
    time       :: Float64
    l_unit     :: Unitful.Length
    m_unit     :: Unitful.Mass
    v_unit     :: Unitful.Velocity
end

"""
Data in the "Header" group of a HDF5 group catalog file.

The default values are for when there are no group catalog files.

# Fields

  - `box_size          :: Float64 = NaN`: Total size of the simulation box, in internal units of length.
  - `h                 :: Float64 = NaN`: Dimensionless Hubble parameter, "little h".
  - `n_groups_part     :: Int32 = -1`: Number of halos (FoF groups) in this file chunk.
  - `n_groups_total    :: Int32 = -1`: Total number of halos (FoF groups) in this snapshot.
  - `n_subgroups_part  :: Int32 = -1`: Number of subhalos (subfind) in this file chunk.
  - `n_subgroups_total :: Int32 = -1`: Total number of subhalos (subfind) in this snapshot.
  - `num_files         :: Int32 = -1`: Number of file chunks per snapshot.
  - `omega_0           :: Float64 = NaN`: The cosmological density parameter for matter.
  - `omega_l           :: Float64 = NaN`: The cosmological density parameter for the cosmological constant.
  - `redshift          :: Float64 = NaN`: The redshift.
  - `time              :: Float64 = NaN`: The physical time or the scale factor, depending on the type of simulation.
"""
@kwdef struct GroupCatHeader
    box_size          :: Float64 = NaN
    h                 :: Float64 = NaN
    n_groups_part     :: Int32 = -1
    n_groups_total    :: Int32 = -1
    n_subgroups_part  :: Int32 = -1
    n_subgroups_total :: Int32 = -1
    num_files         :: Int32 = -1
    omega_0           :: Float64 = NaN
    omega_l           :: Float64 = NaN
    redshift          :: Float64 = NaN
    time              :: Float64 = NaN
end

"""
Metadata for a simulation.

# Fields

  - `path             :: String`: Full path to the simulation directory.
  - `index            :: Int`: An index associated with the simulation.
  - `slice            :: IndexType`: Slice of the simulation, i.e. which snapshots will be read. It can be an integer (a single snapshot), a vector of integers (several snapshots), an `UnitRange` (e.g. 5:13), an `StepRange` (e.g. 5:2:13) or (:) (all snapshots).
  - `cosmological     :: Bool`: If the simulation is cosmological,

      + `false` -> Newtonian simulation    (`ComovingIntegrationOn` = 0).
      + `true`  -> Cosmological simulation (`ComovingIntegrationOn` = 1).
  - `simulation_table :: DataFrame`: A dataframe where each row is a snapshot, and the colums are

      + `:ids`            -> Dataframe index, i.e. if there are 10 snapshots in total it runs from 1 to 10.
      + `:numbers`        -> Number in the file name.
      + `:scale_factors`  -> Scale factor.
      + `:redshifts`      -> Redshift.
      + `:physical_times` -> Physical time since the Big Bang.
      + `:lookback_times` -> Physical time left to reach the last snapshot.
      + `:snapshot_paths` -> Full path to each snapshots.
      + `:groupcat_paths` -> Full path to each group catalog files.
"""
struct Simulation
    path             :: String
    index            :: Int
    slice            :: IndexType
    cosmological     :: Bool
    simulation_table :: DataFrame
end

"""
Metadata for a snapshot.

# Fields

  - `path          :: String`: Full path to the snapshot.
  - `global_index  :: Int`: Index of the snapshot in the context of the whole simulation.
  - `slice_index   :: Int`: Index of the snapshot in the context of the simulation slice.
  - `physical_time :: Unitful.Time`: Physical time since the Big Bang.
  - `lookback_time :: Unitful.Time`: Physical time left to reach the last snapshot.
  - `scale_factor  :: Float64`: Scale factor.
  - `redshift      :: Float64`: Redshift.
  - `header        :: SnapshotHeader`: Header.
"""
struct Snapshot
    path          :: String
    global_index  :: Int
    slice_index   :: Int
    physical_time :: Unitful.Time
    lookback_time :: Unitful.Time
    scale_factor  :: Float64
    redshift      :: Float64
    header        :: SnapshotHeader
end

"""
Metadata for a group catalog file.

# Fields

  - `path   :: Union{String,Missing}`: Full path to the group catalog file.
  - `header :: GroupCatHeader`: Header.
"""
struct GroupCatalog
    path   :: Union{String,Missing}
    header :: GroupCatHeader
end

"""
Unit conversion factors.

# Fields

  - `x_cgs::Unitful.Length`: From internal units of length to ``\\mathrm{cm}``.
  - `x_cosmo::Unitful.Length`: From internal units of length to ``\\mathrm{kpc}``.
  - `x_comoving::Unitful.Length`: From internal units of length to ``\\mathrm{ckpc}``.
  - `v_cgs::Unitful.Velocity`: From internal units of velocity to ``\\mathrm{cm \\, s^{-1}}``.
  - `v_cosmo::Unitful.Velocity`: From internal units of velocity to ``\\mathrm{km \\, s^{-1}}``.
  - `m_cgs::Unitful.Mass`: From internal units of mass to ``\\mathrm{g}``.
  - `m_cosmo::Unitful.Mass`: From internal units of mass to ``\\mathrm{M_\\odot}``.
  - `t_cgs::Unitful.Time`: From internal units of time to ``\\mathrm{s}``.
  - `t_cosmo::Unitful.Time`: From internal units of time to ``\\mathrm{Myr}``.
  - `t_newton::Unitful.Time`: From internal units (non-cosmological simulations) of time to ``\\mathrm{Myr}``.
  - `U_cgs::Unitful.Energy`: From internal units of specific energy to ``\\mathrm{erg \\, g^{-1}}``.
  - `rho_cgs::Unitful.Density`: From internal units of density to ``\\mathrm{g \\, cm^{-3}}``.
  - `P_Pa::Unitful.Pressure`: From internal units of pressure to ``\\mathrm{Pa}``.
"""
struct InternalUnits
    x_cgs      :: Unitful.Length   # From internal units of length to cm
    x_cosmo    :: Unitful.Length   # From internal units of length to kpc
    x_comoving :: Unitful.Length   # From internal units of length to ckpc

    v_cgs      :: Unitful.Velocity # From internal units of velocity to cm * s^-1
    v_cosmo    :: Unitful.Velocity # From internal units of velocity to km * s^-1

    m_cgs      :: Unitful.Mass     # From internal units of mass to g
    m_cosmo    :: Unitful.Mass     # From internal units of mass to M⊙

    t_cgs      :: Unitful.Time     # From internal units of time to s
    t_cosmo    :: Unitful.Time     # From internal units of time to Myr
    t_newton   :: Unitful.Time     # From internal units (non-cosmological simulations) of time to Myr

    U_cgs      :: SpecificEnergy   # From internal units of specific energy to erg * g^-1

    rho_cgs    :: Unitful.Density  # From internal units of density to g * cm^-3

    P_Pa       :: Unitful.Pressure # From internal units of pressure to Pa

    """
        InternalUnits(; <keyword arguments>)

    Constructor for `InternalUnits`.

    For cosmological simulations, and when appropriate, the conversion factors turn comoving units into physical units.

    # Arguments

      - `l_unit :: Unitful.Length=INTERNAL_L_UNIT[]`: Code parameter `UnitLength_in_cm`.
      - `m_unit :: Unitful.Mass=INTERNAL_M_UNIT[]`: Code parameter `UnitMass_in_g`.
      - `v_unit :: Unitful.Velocity=INTERNAL_V_UNIT[]`: Code parameter `UnitVelocity_in_cm_per_s`.
      - `a      :: Float64=1.0`: Cosmological scale factor of the simulation.
      - `h      :: Float64=1.0`: Dimensionless Hubble parameter, "little h".
    """
    function InternalUnits(;
        l_unit :: Unitful.Length=INTERNAL_L_UNIT[],
        m_unit :: Unitful.Mass=INTERNAL_M_UNIT[],
        v_unit :: Unitful.Velocity=INTERNAL_V_UNIT[],
        a      :: Float64=1.0,
        h      :: Float64=1.0,
    )

        #############
        # Base units
        #############

        # Length conversion factors
        x_cgs = l_unit * a / h
        x_cosmo = x_cgs |> u"kpc"
        x_comoving = l_unit / h |> u"kpc"

        # Velocity conversion factors
        v_cgs = v_unit * sqrt(a)
        v_cosmo = v_cgs |> u"km * s^-1"

        # Mass conversion factors
        m_cgs = m_unit / h
        m_cosmo = m_cgs |> u"Msun"

        ################
        # Derived units
        ################

        # Time conversion factors
        t_cgs = x_cgs / v_cgs
        t_cosmo = t_cgs |> u"Myr"
        t_newton = 1.0u"Gyr" |> u"Myr"

        # Specific energy conversion factor
        U_cgs = v_unit^2 |> u"erg * g^-1"

        # Density conversion factor
        rho_cgs = m_cgs * x_cgs^-3

        # Thermal pressure conversion factor (it uses v_unit^2 instead of v_cgs^2,
        # which would add an extra factor of a0)
        P_Pa = v_unit^2 * m_cgs * x_cgs^-3 |> u"Pa"

        new(
            x_cgs,
            x_cosmo,
            x_comoving,
            v_cgs,
            v_cosmo,
            m_cgs,
            m_cosmo,
            t_cgs,
            t_cosmo,
            t_newton,
            U_cgs,
            rho_cgs,
            P_Pa,
        )

    end
end

########
# Grids
########

"""
Grid with 1 d.o.f.

Embeddings:

    1D -> Ticks.
    2D -> Concentric rings.
    3D -> Spherical shells.

# Fields

  - `x_axis      :: Vector{<:Number}`: Position of (the center of) the bins, relative to the `origin`.
  - `x_edges     :: Vector{<:Number}`: Position of (the edges of) the bins, relative to the `origin`.
  - `n_bins      :: Int`: Number of bins.
  - `bin_size_1D :: Vector{<:Number}`: Widths of the bins in linear space.
  - `bin_size_2D :: Vector{<:Number}`: Areas of the bins in linear space (assuming a 2D embedding).
  - `bin_size_3D :: Vector{<:Number}`: Volumes of the bins in linear space (assuming a 3D embedding).
  - `origin      :: Vector{<:Number}`: Physical position of the coordinate origin for the grid.
"""
struct LinearGrid
    x_axis      :: Vector{<:Number}
    x_edges     :: Vector{<:Number}
    n_bins      :: Int
    bin_size_1D :: Vector{<:Number}
    bin_size_2D :: Vector{<:Number}
    bin_size_3D :: Vector{<:Number}
    origin      :: Vector{<:Number}

    """
        LinearGrid(
            start  :: Number,
            stop   :: Number,
            n_bins :: Int;
            <keyword arguments>
        )

    Constructor for `LinearGrid`.

    # Arguments

      - `start  :: Number`: Position of the grid left most edge, relative to the `origin`.
      - `stop   :: Number`: Position of the grid right most edge, relative to the `origin`.
      - `n_bins :: Int`: Number of bins.
      - `origin :: Vector{<:Number}=zeros(typeof(start), 3)`: Physical position of the coordinate origin for the grid (passthrough argument).
      - `log    :: Bool=false`: If the grid will be regular in logarithmic space instead of linear space.
    """
    function LinearGrid(
        start  :: Number,
        stop   :: Number,
        n_bins :: Int;
        origin :: Vector{<:Number}=zeros(typeof(start), 3),
        log    :: Bool=false,
    )

        (
            stop > start ||
            throw(ArgumentError("LinearGrid: `stop` must be larger than `start`, \
            but I got `stop` = $(stop) <= `start` = $(start)"))
        )

        (
            isPositive(n_bins) ||
            throw(ArgumentError("LinearGrid: `n_bins` has to be strictly positive, \
            but I got `n_bins` = $(n_bins) <= 0"))
        )

        if log

            (
                isPositive(start) ||
                throw(ArgumentError("LinearGrid: For a logarithmic grid you need a strictly \
                positive `start`, but I got `start` = $(start) <= 0"))
            )

            # Length unit
            u_l = unit(start)

            log_start = log10(ustrip(start))
            log_stop  = log10(ustrip(u_l, stop))

            width = (log_stop - log_start) / n_bins

            x_axis  = exp10.([(i - 0.5) * width + log_start for i in 1:n_bins]) .* u_l
            x_edges = exp10.([i * width + log_start for i in 0:n_bins]) .* u_l

        else

            width = (stop - start) / n_bins

            x_axis  = [(i - 0.5) * width + start for i in 1:n_bins]
            x_edges = [i * width + start for i in 0:n_bins]

        end

        bin_size_1D = [x_edges[i + 1] - x_edges[i] for i in 1:n_bins]

        if start < zero(start)

            shifted_edges = x_edges .- start
            bin_size_2D   = [area(shifted_edges[i + 1]) - area(shifted_edges[i]) for i in 1:n_bins]
            bin_size_3D   = [volume(shifted_edges[i + 1]) - volume(shifted_edges[i]) for i in 1:n_bins]

        else

            bin_size_2D = [area(x_edges[i + 1]) - area(x_edges[i]) for i in 1:n_bins]
            bin_size_3D = [volume(x_edges[i + 1]) - volume(x_edges[i]) for i in 1:n_bins]

        end

        new(x_axis, x_edges, n_bins, bin_size_1D, bin_size_2D, bin_size_3D, origin)

    end

    """
        LinearGrid(
            x_edges :: Vector{<:Number};
            <keyword arguments>
        )

    Constructor for `LinearGrid`.

    # Arguments

      - `x_edges :: Vector{<:Number}`: Relative position of (the edges of) the bins, relative to the `origin`.
      - `origin  :: Vector{<:Number}=zeros(runtimeType(x_edges), 3)`: Physical position of the coordinate origin for the grid (passthrough argument).
      - `log     :: Bool=false`: If the `x_edges` are logarithmic, which means that `x_axis` will mark the center of the bins in logarithmic space instead of linear space.
    """
    function LinearGrid(
        x_edges :: Vector{<:Number};
        origin  :: Vector{<:Number}=zeros(runtimeType(x_edges), 3),
        log     :: Bool=falsem
    )

        isempty(x_edges) && throw(ArgumentError("LinearGrid: `x_edges` is empty!"))

        n_bins = length(x_edges) - 1

        issorted(x_edges) || sort!(x_edges)

        start = first(x_edges)

        bin_size_1D = [x_edges[i + 1] - x_edges[i] for i in 1:n_bins]

        if start < zero(start)

            shifted_edges = x_edges .- start
            bin_size_2D   = [area(shifted_edges[i + 1]) - area(shifted_edges[i]) for i in 1:n_bins]
            bin_size_3D   = [volume(shifted_edges[i + 1]) - volume(shifted_edges[i]) for i in 1:n_bins]

        else

            bin_size_2D = [area(x_edges[i + 1]) - area(x_edges[i]) for i in 1:n_bins]
            bin_size_3D = [volume(x_edges[i + 1]) - volume(x_edges[i]) for i in 1:n_bins]

        end

        if log

            (
                isPositive(start) ||
                throw(ArgumentError("LinearGrid: For a logarithmic grid you need a strictly \
                positive `start`, but I got `start` = $(start) <= 0"))
            )

            # Length unit
            u_l = unit(start)

            log_x_edges = log10.(ustrip.(u_l, x_edges))

            x_axis = exp10.([(log_x_edges[i + 1] + log_x_edges[i]) / 2.0 for i in 1:n_bins]) .* u_l

        else

            x_axis  = [(x_edges[i + 1] + x_edges[i]) / 2.0 for i in 1:n_bins]

        end

        new(x_axis, x_edges, n_bins, bin_size_1D, bin_size_2D, bin_size_3D, origin)

    end
end

"""
Grid with 2 d.o.f.

# Fields

  - `x_axis      :: Vector{<:Number}`: Absolute position of (the center of) the bins in the x direction.
  - `y_axis      :: Vector{<:Number}`: Absolute position of (the center of) the bins in the y direction.
  - `x_edges     :: Vector{<:Number}`: Absolute position of (the edges of) the bins in the x direction.
  - `y_edges     :: Vector{<:Number}`: Absolute position of (the edges of) the bins in the y direction.
  - `n_bins      :: NTuple{2,Int}`: Number of bins in each direction.
  - `bin_size_1D :: Vector{<:Number}`: Widths of the bins in each direction.
  - `bin_size_2D :: Number`: Area of the bins.
  - `origin      :: Vector{<:Number}`: Physical position of the center of the grid.
  - `size        :: Vector{<:Number}`: Physical side length of the grid in each direction (passthrough argument).
  - `n_pixels    :: Int`: Total number of bins.
"""
struct SquareGrid
    x_axis      :: Vector{<:Number}
    y_axis      :: Vector{<:Number}
    x_edges     :: Vector{<:Number}
    y_edges     :: Vector{<:Number}
    n_bins      :: NTuple{2, Int}
    bin_size_1D :: Vector{<:Number}
    bin_size_2D :: Number
    origin      :: Vector{<:Number}
    size        :: Vector{<:Number}
    n_pixels    :: Int

    """
        SquareGrid(
            size   :: Vector{<:Number},
            n_bins :: NTuple{2,Int};
            <keyword arguments>
        )

    Constructor for `SquareGrid`.

    # Arguments

      - `size   :: Vector{<:Number}`: Physical side length of the grid in each direction.
      - `n_bins :: NTuple{2,Int}`: Number of bins in each direction.
      - `origin :: Vector{<:Number}=zeros(runtimeType(grid_length), 3)`: Physical position of the center of the grid.
    """
    function SquareGrid(
        size   :: Vector{<:Number},
        n_bins :: NTuple{2,Int};
        origin :: Vector{<:Number}=zeros(runtimeType(grid_length), 3),
    )

        (
            all(isPositive, size) ||
            throw(ArgumentError("SquareGrid: `size` has to be strictly positive, \
            but I got `size` = $(size)"))
        )

        (
            all(isPositive, n_bins) ||
            throw(ArgumentError("SquareGrid: `n_bins` has to be strictly positive, \
            but I got `n_bins` = $(n_bins)"))
        )

        # Compute the bin dimensions
        bin_size_1D = size ./ n_bins
        bin_size_2D = prod(bin_size_1D)

        # Compute the x and y coordinates of the center of each bin
        center_shift  = 0.5 .* (size .- bin_size_1D)
        x_axis = [(i - 1) * bin_size_1D[1] - center_shift[1] + origin[1] for i in 1:n_bins[1]]
        y_axis = [(i - 1) * bin_size_1D[2] - center_shift[2] + origin[2] for i in 1:n_bins[2]]


        # Compute the x and y coordinates of the edges of each bin
        edge_shift  = 0.5 .* size
        x_edges = [i * bin_size_1D[1] - edge_shift[1] + origin[1] for i in 0:n_bins[1]]
        y_edges = [i * bin_size_1D[2] - edge_shift[2] + origin[2] for i in 0:n_bins[2]]

        n_pixels = prod(n_bins)

        new(
            x_axis,
            y_axis,
            x_edges,
            y_edges,
            n_bins,
            bin_size_1D,
            bin_size_2D,
            origin,
            size,
            n_pixels,
        )

    end

    """
        SquareGrid(
            size   :: Number,
            n_bins :: Int;
            <keyword arguments>
        )

    Constructor for `SquareGrid`.

    # Arguments

      - `size   :: Number`: Physical side length of the grid in every direction.
      - `n_bins :: Int`: Number of bins in every direction.
      - `origin :: Vector{<:Number}=zeros(typeof(size), 3)`: Physical position of the center of the grid.
    """
    function SquareGrid(
        size   :: Number,
        n_bins :: Int;
        origin :: Vector{<:Number}=zeros(typeof(size), 3),
    )

        SquareGrid([size, size], (n_bins, n_bins); origin)

    end
end

"""
Grid with 3 d.o.f.

# Fields

  - `x_axis      :: Vector{<:Number}`: Absolute position of (the center of) the bins in the x direction.
  - `y_axis      :: Vector{<:Number}`: Absolute position of (the center of) the bins in the y direction.
  - `z_axis      :: Vector{<:Number}`: Absolute position of (the center of) the bins in the z direction.
  - `x_edges     :: Vector{<:Number}`: Absolute position of (the edges of) the bins in the x direction.
  - `y_edges     :: Vector{<:Number}`: Absolute position of (the edges of) the bins in the y direction.
  - `z_edges     :: Vector{<:Number}`: Absolute position of (the edges of) the bins in the z direction.
  - `n_bins      :: NTuple{3,Int}`: Number of bins in each direction.
  - `bin_size_1D :: Vector{<:Number}`: Widths of the bins in each direction.
  - `bin_size_2D :: Vector{<:Number}`: Areas of the faces normal to each direction.
  - `bin_size_3D :: Number`: Volume of the bins.
  - `origin      :: Vector{<:Number}`: Physical position of the center of the grid.
  - `size        :: Vector{<:Number}`: Physical side length of the grid in each direction (passthrough argument).
  - `n_voxels    :: Int`: Total number of bins.
"""
struct CubicGrid
    x_axis      :: Vector{<:Number}
    y_axis      :: Vector{<:Number}
    z_axis      :: Vector{<:Number}
    x_edges     :: Vector{<:Number}
    y_edges     :: Vector{<:Number}
    z_edges     :: Vector{<:Number}
    n_bins      :: NTuple{3,Int}
    bin_size_1D :: Vector{<:Number}
    bin_size_2D :: Vector{<:Number}
    bin_size_3D :: Number
    origin      :: Vector{<:Number}
    size        :: Vector{<:Number}
    n_voxels    :: Int

    """
        CubicGrid(
            size   :: Vector{<:Number},
            n_bins :: NTuple{3,Int};
            <keyword arguments>
        )

    Constructor for `CubicGrid`.

    # Arguments

      - `size   :: Vector{<:Number}`: Physical side length of the grid in each direction.
      - `n_bins :: NTuple{3,Int}`: Number of bins in each direction.
      - `origin :: Vector{<:Number}=zeros(runtimeType(grid_length), 3)`: Physical position of the center of the grid.
    """
    function CubicGrid(
        size   :: Vector{<:Number},
        n_bins :: NTuple{3,Int};
        origin :: Vector{<:Number}=zeros(runtimeType(grid_length), 3),
    )

        (
            all(isPositive, size) ||
            throw(ArgumentError("CubicGrid: `size` has to be strictly positive, \
            but I got `size` = $(size)"))
        )

        (
            all(isPositive, n_bins) ||
            throw(ArgumentError("CubicGrid: `n_bins` has to be strictly positive, \
            but I got `n_bins` = $(n_bins)"))
        )

        # Compute the bin dimensions
        bin_size_1D = size ./ n_bins
        bin_size_2D = prod.([bin_size_1D[[2,3]], bin_size_1D[[1,3]], bin_size_1D[[1,2]]])
        bin_size_3D = prod(bin_size_1D)

        # Compute the x and y coordinates of the center of each bin
        center_shift  = 0.5 .* (size .- bin_size_1D)
        x_axis = [(i - 1) * bin_size_1D[1] - center_shift[1] + origin[1] for i in 1:n_bins[1]]
        y_axis = [(i - 1) * bin_size_1D[2] - center_shift[2] + origin[2] for i in 1:n_bins[2]]
        z_axis = [(i - 1) * bin_size_1D[3] - center_shift[3] + origin[3] for i in 1:n_bins[3]]

        # Compute the x and y coordinates of the edges of each bin
        edge_shift  = 0.5 .* size
        x_edges = [i * bin_size_1D[1] - edge_shift[1] + origin[1] for i in 0:n_bins[1]]
        y_edges = [i * bin_size_1D[2] - edge_shift[2] + origin[2] for i in 0:n_bins[2]]
        z_edges = [i * bin_size_1D[3] - edge_shift[3] + origin[3] for i in 0:n_bins[3]]

        n_voxels = prod(n_bins)

        new(
            x_axis,
            y_axis,
            z_axis,
            x_edges,
            y_edges,
            z_edges,
            n_bins,
            bin_size_1D,
            bin_size_2D,
            bin_size_3D,
            origin,
            size,
            n_voxels,
        )

    end

    """
        CubicGrid(
            size   :: Number,
            n_bins :: Int;
            <keyword arguments>
        )

    Constructor for `CubicGrid`.

    # Arguments

      - `size   :: Number`: Physical side length of the grid in every direction.
      - `n_bins :: Int`: Number of bins in every direction.
      - `origin :: Vector{<:Number}=zeros(typeof(size), 3)`: Physical position of the center of the grid.
    """
    function CubicGrid(
        size   :: Number,
        n_bins :: Int;
        origin :: Vector{<:Number}=zeros(typeof(size), 3),
    )

        CubicGrid([size, size, size], (n_bins, n_bins, n_bins); origin)

    end
end
