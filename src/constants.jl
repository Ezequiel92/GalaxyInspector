####################################################################################################
# Constants and data structures.
####################################################################################################

# Configuration

"""
If physical lengths will be used throughout, instead of comoving lengths.
"""
PHYSICAL_UNITS = false

"""
Base name of the snapshot files, set in the code variable `SnapshotFileBase`.
"""
const SNAP_BASENAME = "snap"

"""
Base name of the group catalog files.
"""
const GC_BASENAME = "fof_subhalo_tab"

"""
Relative path, within the simulation directory, to the `sfr.txt` file.
"""
const SFR_REL_PATH = "output/sfr.txt"

"""
Relative path, within the simulation directory, to the `cpu.txt` file.
"""
const CPU_REL_PATH = "output/cpu.txt"

"""
Mass fraction of hydrogen for the gas cells.
"""
const HYDROGEN_MASSFRAC = 0.76

"""
Mass of the tracers in internal units. It values comes from `All.TargetGasMass = All.TargetGasMassFactor * All.ReferenceGasPartMass` in the Arepo code.
"""
const TRACER_MASS = 3.65456e-06

"""
Internal unit of length used in IllustrisTNG, equivalent to ``1.0  \\, \\mathrm{kpc}``.
See the documentation [here](https://www.tng-project.org/data/docs/specifications/)
"""
const ILLUSTRIS_L_UNIT = 3.085678e21u"cm"

"""
Internal unit of mass used in IllustrisTNG, equivalent to ``10^{10} \\, \\mathrm{M_\\odot}``.
See the documentation [here](https://www.tng-project.org/data/docs/specifications/)
"""
const ILLUSTRIS_M_UNIT = 1.989e43u"g"

"""
Internal unit of velocity used in IllustrisTNG, equivalent to ``1.0 \\, \\mathrm{km \\, s^{-1}}``.
See the documentation [here](https://www.tng-project.org/data/docs/specifications/)
"""
const ILLUSTRIS_V_UNIT = 1.0e5u"cm*s^-1"

"""
Characteristic radius used throughout the code.
"""
const FILTER_R = 30.0u"kpc"

"""
Stellar age limit for the SFR and sSFR calculations.
"""
const AGE_RESOLUTION = 25.0u"Myr"

"""
Stellar age limit for the SFR area density calculation.
"""
const AGE_RESOLUTION_ρ = 200.0u"Myr"

"""
Threshold density, above which the gas particles enter the star formation routine.

This value corresponds to `CritPhysDensity` ``= 0.318 \\, [\\mathrm{cm^{-3}}]`` in the `param.txt` file. Which is converted to internal units within the code using `PhysDensThresh = CritPhysDensity * PROTONMASS / HYDROGEN_MASSFRAC / UnitDensity_in_cgs`. Then, to go to physical units again ones does: PhysDensThresh * UnitDensity_in_cgs * cf_a3inv * HubbleParam * HubbleParam. For a cosmological simulation at readshift 0 (cf_a3inv = 1), this result in a physical density threshold of ``0.192 \\, [\\mathrm{cm^{-3}}]``, or adding the proton mass a value of:
"""
const THRESHOLD_DENSITY = 4.749326307150211e6u"Msun*kpc^-3"

@doc raw"""
Hubble constant in $\mathrm{Gyr^{-1}}$.

This value corresponds to $H_0 = 0.102201 \, \mathrm{Gyr}^{-1} = 100 \, \mathrm{km} \, \mathrm{s}^{-1} \, \mathrm{Mpc}^{-1}$.
"""
const HUBBLE_CONSTANT = 0.102201

"""
Solar metallicity, as used in Arepo.

# References

M. Asplund et al. (2006). *The new solar abundances - Part I: the observations*. Communications in Asteroseismology, **147**. [doi:10.1553/cia147s76](https://doi.org/10.1553/cia147s76)
"""
const SOLAR_METALLICITY = 0.0127

# Cell/particle types

"""
Code index for each type of cell/particle.

!!! note

    This index is for simulations with 7 cell/particle types.

# References

See for example Gadget2 [User's Guide](https://wwwmpa.mpa-garching.mpg.de/gadget/users-guide.pdf), or Gadget4 [documentation](https://wwwmpa.mpa-garching.mpg.de/gadget4/).
"""
const FULL_PARTICLE_INDEX = Dict(
    :gas        => 0,
    :halo       => 1, # High resolution dark matter
    :disk       => 2, # Intermediate resolution dark matter
    :bulge      => 3, # Low resolution dark matter
    :stars      => 4, # Star particles and wind particles
    :black_hole => 5,
    :tracer     => 6,
)

"""
Code index for each type of cell/particle.

!!! note

    This index is for simulations with 6 cell/particle types.

# References

See for example Gadget2 [User's Guide](https://wwwmpa.mpa-garching.mpg.de/gadget/users-guide.pdf), or Gadget4 [documentation](https://wwwmpa.mpa-garching.mpg.de/gadget4/).
"""
const DEFAULT_PARTICLE_INDEX = Dict(
    :gas        => 0,
    :halo       => 1,
    :disk       => 2, # Generally not used
    :bulge      => 3, # Generally not used
    :stars      => 4,
    :black_hole => 5,
)

"""
Current cell/particle index in use.
"""
const PARTICLE_INDEX = FULL_PARTICLE_INDEX

"""
Human readable name corresponding to each type of cell/particle.
"""
const PARTICLE_NAMES = Dict(
    :gas        => "Gas cells",
    :halo       => "HR DM particles",
    :disk       => "IR DM particles",
    :bulge      => "LR DM particles",
    :stars      => "Stellar particles",
    :black_hole => "Black hole particles",
    :tracer     => "Tracer particles",
)

"""
Type of cell/particle corresponding to each code index.
"""
const INDEX_PARTICLE = Dict(n => symbol for (symbol, n) in PARTICLE_INDEX)

"""
Type of cell/particle corresponding to each internal code name (data group in the HDF5 output).
"""
const PARTICLE_TYPE = Dict("PartType$n" => symbol for (symbol, n) in PARTICLE_INDEX)

"""
Internal code name (data group in the HDF5 output) corresponding to each type of cell/particle.
"""
const PARTICLE_CODE_NAME = Dict(symbol => "PartType$n" for (symbol, n) in PARTICLE_INDEX)

"""
Human readable name corresponding to each morphological component.
"""
const MORPHOLOGICAL_COMPONENTS = Dict(
    :disk        => "Disk",
    :bulge       => "Bulge",
    :thin_disk   => "Thin disk",
    :thick_disk  => "Thick disk",
)

"""
Default filter dictionary that does not exclude any cells/particles.
"""
const PASS_ALL = Dict(key => (:) for key in keys(PARTICLE_INDEX))

"""
Filter out every cell/particle.
"""
const PASS_NONE = Dict(key => Int[] for key in keys(PARTICLE_INDEX))

# Tracked elements

"""
Code index for each tracked element.
"""
const ELEMENT_INDEX = Dict(
    :H     => 1,  # Hydrogen
    :He    => 2,  # Helium
    :C     => 3,  # Carbon
    :N     => 4,  # Nitrogen
    :O     => 5,  # Oxygen
    :Ne    => 6,  # Neon
    :Mg    => 7,  # Magnesium
    :Si    => 8,  # Silicon
    :Fe    => 9,  # Iron
    :Other => 10, # All other
)

"""
Symbol list for the gas abundance quantities.
"""
const GAS_ABUNDANCE = [Symbol(element, "_gas_abundance") for element in keys(ELEMENT_INDEX)]

"""
Symbol list for the stellar abundance quantities.
"""
const STELLAR_ABUNDANCE = [Symbol(element, "_stellar_abundance") for element in keys(ELEMENT_INDEX)]

# Reference values from the literature

@doc raw"""
Solar abundances.

They are defined as $12 + \log_{10}(N_\mathrm{X} / N_\mathrm{H})$, where $N_\mathrm{X}$ and $N_\mathrm{H}$ are the number densities of element $\mathrm{X}$ and hydrogen respectively.

# References

M. Asplund et al. (2009). *The Chemical Composition of the Sun*. Annual Review of Astronomy and Astrophysics, **47(1)**, 481–522. [doi:10.1146/annurev.astro.46.060407.145222](https://doi.org/10.1146/annurev.astro.46.060407.145222)
"""
const SOLAR_ABUNDANCE = Dict(
    :H  => 12,    # Hydrogen
    :He => 10.93, # Helium
    :C  => 8.43,  # Carbon
    :N  => 7.83,  # Nitrogen
    :O  => 8.69,  # Oxygen
    :Ne => 7.93,  # Neon
    :Mg => 7.60,  # Magnesium
    :Si => 7.51,  # Silicon
    :Fe => 7.50,  # Iron
)

"""
Standard atomic weights.

# References

T. Prohaska et al. (2022). *Standard atomic weights of the elements 2021 (IUPAC Technical Report)*. Pure and Applied Chemistry, **94(5)**, 573-600. [doi:10.1515/pac-2019-0603](https://doi.org/10.1515/pac-2019-0603)
"""
const ATOMIC_WEIGHTS = Dict(
    :H  => 1.0080, # Hydrogen
    :He => 4.0026, # Helium
    :C  => 12.011, # Carbon
    :N  => 14.007, # Nitrogen
    :O  => 15.999, # Oxygen
    :Ne => 20.180, # Neon
    :Mg => 24.305, # Magnesium
    :Si => 28.085, # Silicon
    :Fe => 55.845, # Iron
)

"""
Path to the file with the fits for the molecular Kennicutt-Schmidt relation, taken from Bigiel et al. (2008).

# References

F. Bigiel et al. (2008). *THE STAR FORMATION LAW IN NEARBY GALAXIES ON SUB-KPC SCALES*. The Astrophysical Journal, **136(6)**, 2846. [doi:10.1088/0004-6256/136/6/2846](https://doi.org/10.1088/0004-6256/136/6/2846)
"""
const BIGIEL2008_DATA_PATH = joinpath(@__DIR__, "../experimental_data/Bigiel2008.txt")

"""
Slope of the Kennicutt-Schmidt law, taken from Kennicutt (1998).

# References

R. C. Kennicutt (1998). *The Global Schmidt Law in Star-forming Galaxies*. The Astrophysical Journal, **498(2)**, 541-552. [doi:10.1086/305588](https://doi.org/10.1086/305588)
"""
const KS98_SLOPE = 1.4

"""
Intercept of the Kennicutt-Schmidt law, taken from Kennicutt (1998).

# References

R. C. Kennicutt (1998). *The Global Schmidt Law in Star-forming Galaxies*. The Astrophysical Journal, **498(2)**, 541-552. [doi:10.1086/305588](https://doi.org/10.1086/305588)
"""
const KS98_INTERCEPT = 2.5e-4u"Msun*yr^-1*kpc^-2"

"""
Unit of surface density of the Kennicutt-Schmidt law, taken from Kennicutt (1998).

# References

R. C. Kennicutt (1998). *The Global Schmidt Law in Star-forming Galaxies*. The Astrophysical Journal, **498(2)**, 541-552. [doi:10.1086/305588](https://doi.org/10.1086/305588)
"""
const KS98_RHO_UNIT = u"Msun*pc^-2"

"""
Path to the file with the Milky Way profiles, taken from Mollá et al. (2015).

# References

M. Mollá et al. (2015). *Galactic chemical evolution: stellar yields and the initial mass function*. Monthly Notices of the Royal Astronomical Society **451(4)**, 3693–3708. [doi:10.1093/mnras/stv1102](https://doi.org/10.1093/mnras/stv1102)
"""
const MOLLA2015_DATA_PATH = joinpath(@__DIR__, "../experimental_data/Mollá2015.csv")

"""
Path to the file with the global galaxy properties, taken from Feldmann (2020).

# References

R. Feldmann (2020). *The link between star formation and gas in nearby galaxies*. Communications Physics **3(226)**. [doi:10.1038/s42005-020-00493-0](https://doi.org/10.1038/s42005-020-00493-0)
"""
const FELDMANN2020_DATA_PATH = joinpath(@__DIR__, "../experimental_data/Feldmann2020.csv")

"""
Reference pressure for the molecular fraction-pressure relation, taken from Blitz et al. (2006).

# References

Table 2 - `Mean` row - Third column.

L. Blitz et al. (2006). *The Role of Pressure in GMC Formation II: The H2-Pressure Relation*. The Astrophysical Journal, **650(2)**, 933. [doi:10.1086/505417](https://doi.org/10.1086/505417)
"""
const P0 = 3.5e4u"K*cm^-3" * Unitful.k

"""
Reference exponent for the molecular fraction-pressure relation, taken from Blitz et al. (2006).

# References

Table 2 - `Mean` row - Second column. (We use -α here).

L. Blitz et al. (2006). *The Role of Pressure in GMC Formation II: The H2-Pressure Relation*. The Astrophysical Journal, **650(2)**, 933. [doi:10.1086/505417](https://doi.org/10.1086/505417)
"""
const ALPHA_BLITZ = -0.92

# Type aliases

"""
Color type.
"""
const ColorType = Union{ColorTypes.RGB,Symbol}

"""
Line style type.
"""
const LineStyleType = Union{Tuple{String,Symbol},Nothing,String,Symbol}

"""
Index type.
"""
const IndexType = Union{
    Colon,
    Integer,
    UnitRange{<:Integer},
    StepRange{<:Integer,<:Integer},
    Vector{<:Integer},
    Vector{Bool},
}

# Dimensions of specific energy
@derived_dimension SpecificEnergy Unitful.𝐋^2 * Unitful.𝐓^-2 true

# Dimensions of surface density
@derived_dimension SurfaceDensity Unitful.𝐌 * Unitful.𝐋^-2 true

# Dimensions of mass flow surface density
@derived_dimension MassFlowDensity Unitful.𝐌 * Unitful.𝐓^-1 * Unitful.𝐋^-2 true

# Dimensions of angular momentum
@derived_dimension AngularMomentum Unitful.𝐌 * Unitful.𝐋^2 * Unitful.𝐓^-1 true

# Makie.jl configuration

"""
Default list of marker types.
"""
const MARKERS = [:circle, :rect, :diamond, :hexagon, :cross, :xcross, :pentagon]

"""
Default list of line styles.
"""
const LINE_STYLES = [:solid, :dash, :dot, :dashdot, :dashdotdot]

"""
Global plot theme.

On the graphic units used:

We know that 1 mm = 2.83466 pt and 1 in = 25,4 mm. Then, if we want 1 [code ]unit = 0.1 mm in vector graphics, we have to use `pt_per_unit` = 0.283466.
For pixel images, we control the ppi with `px_per_unit`. A resonable high ppi is 600, so, using `px_per_unit` = 2.3622 we get 23.622 px/mm ~ 600 px/in (remember that 1 [code ]unit = 0.1 mm).
"""
const DEFAULT_THEME = Theme(
    ################################################################################################
    # Size of the figures in code units.
    # For PDFs and SVGs, 880 [code ]unit = 8.8 cm.
    # For PNGs, when printed to a size of 1 point = 0.1 mm, one will get a dpi of 600 (23.622 px/mm).
    ################################################################################################
    size=(880, 880),
    #####################################
    # 35 unit * 0.283466 pt/unit ~ 9.9 pt
    #####################################
    fontsize=35,
    ############################
    # (left, right, bottom, top)
    ############################
    figure_padding=(1, 15, 5, 15),
    palette=(color=Makie.wong_colors(), marker=MARKERS, linestyle=LINE_STYLES),
    CairoMakie=(px_per_unit=2.3622, pt_per_unit=0.283466),
    Axis=(
        xlabelpadding=15,
        xticklabelpad=10,
        xticksize=7,
        xgridvisible=false,
        spinewidth=3,
        xminorticksvisible=true,
        xminorticks=IntervalsBetween(5),
        ylabelpadding=15,
        yticklabelpad=10,
        yticksize=7,
        ygridvisible=false,
        yminorticksvisible=true,
        yminorticks=IntervalsBetween(5),
        ############################################################################################
        # Aspect ratio of the figures. The options are:
        # nothing: The aspect ratio will be chosen by [Makie](https://docs.makie.org/stable/).
        # AxisAspect(n): The aspect ratio will be given by the number `n` = width / height.
        # DataAspect(): The aspect ratio of the data will be used.
        ############################################################################################
        aspect=AxisAspect(1),
    ),
    Legend=(
        tellheight=false,
        tellwidth=false,
        framevisible=false,
        colgap=20,
        halign=:right,
        valign=:bottom,
        nbanks=3,
        titlegap=-5,
        labelsize=30,
        linewidth=5,
        markersize=28,
        patchsize=(50, 50),
        linepoints=[Point2f(0.0, 0.5), Point2f(0.9, 0.5)],
        ##############################################
        # Vertices, relative to the default 1x1 square
        ##############################################
        polypoints=[
            Point2f(0.15, 0.15),
            Point2f(0.85, 0.15),
            Point2f(0.85, 0.85),
            Point2f(0.15, 0.85),
        ],
    ),
    Lines=(linewidth=5, cycle=Cycle([:color, :linestyle], covary=true)),
    VLines=(linewidth=3, cycle=Cycle([:color, :linestyle], covary=true)),
    HLines=(linewidth=3, cycle=Cycle([:color, :linestyle], covary=true)),
    ScatterLines=(
        linewidth=5,
        markersize=22,
        cycle=Cycle([:color, :linestyle, :marker], covary=true),
    ),
    Scatter=(markersize=22, cycle=Cycle([:color, :marker], covary=true)),
    Band=(cycle=[:color],),
    Errorbars=(whiskerwidth=10,),
    Heatmap=(colormap=:CMRmap, nan_color=ColorSchemes.CMRmap[1]),
    Colorbar=(size=25, ticklabelpad=10, minorticksvisible=true, ticksize=7, labelpadding=2),
    BarPlot=(
        color_over_background=:black,
        color_over_bar=:black,
        flip_labels_at=10,
        direction=:x,
        strokecolor=:black,
        strokewidth=1,
        bar_labels=:y,
        dodge_gap=0.04,
    ),
    Arrows=(lengthscale=0.02, arrowsize=7.0, linestyle=:solid, color=:white),
)

# Structures

"""
Data in the "Header" group of a HDF5 snapshot file.

# Fields

  - `box_size::Float64`: Total size of the simulation box.
  - `h0::Float64`: Hubble parameter.
  - `mass_table::Vector{Float64}`: Masses of particle types which have a constant mass.
  - `num_files::Int32`: Number of file chunks per snapshot.
  - `num_part::Vector{Int32}`: Number of particles (of each type) included in this file chunk.
  - `num_total::Vector{UInt32}`: Total number of particles (of each type) for this snapshot.
  - `omega_0::Float64`: The cosmological density parameter for matter.
  - `omega_l::Float64`: The cosmological density parameter for the cosmological constant.
  - `redshift::Float64`: The redshift.
  - `time::Float64`: The physical time/scale factor.
  - `l_unit::Unitful.Length`: Conversion factor from internal units of length to centimeters.
  - `m_unit::Unitful.Mass`: Conversion factor from internal units of mass to grams.
  - `v_unit::Unitful.Velocity`: Conversion factor from internal units of velocity to centimeters per second.
"""
@kwdef mutable struct SnapshotHeader
    box_size::Float64
    h0::Float64
    mass_table::Vector{Float64}
    num_files::Int32
    num_part::Vector{Int32}
    num_total::Vector{UInt32}
    omega_0::Float64
    omega_l::Float64
    redshift::Float64
    time::Float64
    l_unit::Unitful.Length
    m_unit::Unitful.Mass
    v_unit::Unitful.Velocity
end

"""
Data in the "Header" group of a HDF5 group catalog file.

# Fields

  - `box_size::Float64 = NaN`: Total size of the simulation box.
  - `h0::Float64 = NaN`: Hubble parameter.
  - `n_groups_part::Int32 = -1`: Number of halos (FoF groups) in this file chunk.
  - `n_groups_total::Int32 = -1`: Total number of halos (FoF groups) in this snapshot.
  - `n_subgroups_part::Int32 = -1`: Number of subhalos (subfind) in this file chunk.
  - `n_subgroups_total::Int32 = -1`: Total number of subhalos (subfind) in this snapshot.
  - `num_files::Int32 = -1`: Number of file chunks per snapshot.
  - `omega_0::Float64 = NaN`: The cosmological density parameter for matter.
  - `omega_l::Float64 = NaN`: The cosmological density parameter for the cosmological constant.
  - `redshift::Float64 = NaN`: The redshift.
  - `time::Float64 = NaN`: The physical time/scale factor.
"""
@kwdef mutable struct GroupCatHeader
    box_size::Float64 = NaN
    h0::Float64 = NaN
    n_groups_part::Int32 = -1
    n_groups_total::Int32 = -1
    n_subgroups_part::Int32 = -1
    n_subgroups_total::Int32 = -1
    num_files::Int32 = -1
    omega_0::Float64 = NaN
    omega_l::Float64 = NaN
    redshift::Float64 = NaN
    time::Float64 = NaN
end

"""
Metadata for a simulation.

# Fields

  - `path::String`: Full path to the simulation directory.
  - `index::Int`: An index associated with the simulation.
  - `slice::IndexType`: Slice of the simulation, i.e. which snapshots will be read. It can be an integer (a single snapshot), a vector of integers (several snapshots), an `UnitRange` (e.g. 5:13), an `StepRange` (e.g. 5:2:13) or (:) (all snapshots).
  - `cosmological::Bool`: If the simulation is cosmological,

      + `false` -> Newtonian simulation    (`ComovingIntegrationOn` = 0).
      + `true`  -> Cosmological simulation (`ComovingIntegrationOn` = 1).
  - `table::DataFrame`: A dataframe where each row is a snapshot, and with 8 colums:

      + `:ids`            -> Dataframe index of each snapshot, i.e. if there are 10 snapshots in total it runs from 1 to 10.
      + `:numbers`        -> Number in the file name of each snapshot.
      + `:scale_factors`  -> Scale factor of each snapshot.
      + `:redshifts`      -> Redshift of each snapshot.
      + `:physical_times` -> Physical time since the Big Bang.
      + `:lookback_times` -> Physical time left to reach the last snapshot.
      + `:snapshot_paths` -> Full path to the snapshots.
      + `:groupcat_paths` -> Full path to the group catalog files.
"""
struct Simulation
    path::String
    index::Int
    slice::IndexType
    cosmological::Bool
    table::DataFrame
end

"""
Metadata for a snapshot.

# Fields

  - `path::String`: Full path to the snapshot.
  - `global_index::Int`: Index of the snapshot in the context of the whole simulation.
  - `slice_index::Int`: Index of the snapshot in the context of the slice.
  - `physical_time::Unitful.Time`: Physical time since the Big Bang.
  - `lookback_time::Unitful.Time`: Physical time left to reach the last snapshot.
  - `scale_factor::Float64`: Scale factor of the snapshot.
  - `redshift::Float64`: Redshift of the snapshot.
  - `header::SnapshotHeader`: Header of the snapshot.
"""
struct Snapshot
    path::String
    global_index::Int
    slice_index::Int
    physical_time::Unitful.Time
    lookback_time::Unitful.Time
    scale_factor::Float64
    redshift::Float64
    header::SnapshotHeader
end

"""
Metadata for a group catalog file.

# Fields

  - `path::Union{String,Missing}`: Full path to the group catalog file.
  - `header::GroupCatHeader`: Header of the group catalog.
"""
struct GroupCatalog
    path::Union{String,Missing}
    header::GroupCatHeader
end

"""
Unit conversion struct.

# Fields

  - `x_cgs::Unitful.Length`: Length, from internal units to ``\\mathrm{cm}``.
  - `x_cosmo::Unitful.Length`: Length, from internal units to ``\\mathrm{kpc}``.
  - `x_comoving::Unitful.Length`: Length, from internal units to ``\\mathrm{ckpc}``.
  - `v_cgs::Unitful.Velocity`: Velocity, from internal units to ``\\mathrm{cm \\, s^{-1}}``.
  - `v_cosmo::Unitful.Velocity`: Velocity, from internal units to ``\\mathrm{km \\, s^{-1}}``.
  - `m_cgs::Unitful.Mass`: Mass, from internal units to ``\\mathrm{g}``.
  - `m_cosmo::Unitful.Mass`: Mass, from internal units to ``\\mathrm{M_\\odot}``.
  - `t_cgs::Unitful.Time`: Time, from internal units to ``\\mathrm{s}``.
  - `t_cosmo::Unitful.Time`: Time, from internal units to ``\\mathrm{Myr}``.
  - `U_cgs::Unitful.Energy`: Specific energy, from internal units to ``\\mathrm{erg \\, g^{-1}}``.
  - `rho_cgs::Unitful.Density`: Density, from internal units to ``\\mathrm{g \\, cm^{-3}}``.
  - `P_Pa::Unitful.Pressure`: Pressure, from internal units to ``\\mathrm{Pa}``.
"""
struct InternalUnits

    x_cgs::Unitful.Length      # Length, from internal units to cm
    x_cosmo::Unitful.Length    # Length, from internal units to kpc
    x_comoving::Unitful.Length # Length, from internal units to ckpc

    v_cgs::Unitful.Velocity    # Velocity, from internal units to cm * s^-1
    v_cosmo::Unitful.Velocity  # Velocity, from internal units to km * s^-1

    m_cgs::Unitful.Mass        # Mass, from internal units to g
    m_cosmo::Unitful.Mass      # Mass, from internal units to M⊙

    t_cgs::Unitful.Time        # Time, from internal units to s
    t_cosmo::Unitful.Time      # Time, from internal units to Myr

    U_cgs::SpecificEnergy      # Specific energy, from internal units to erg * g^-1

    rho_cgs::Unitful.Density   # Density, from internal units to g * cm^-3

    P_Pa::Unitful.Pressure     # Pressure, from internal units to Pa

    """
        InternalUnits(; <keyword arguments>)

    Constructor for `InternalUnits`.

    # Arguments

      - `l_unit::Unitful.Length=ILLUSTRIS_L_UNIT`: Code parameter `UnitLength_in_cm`.
      - `m_unit::Unitful.Mass=ILLUSTRIS_M_UNIT`: Code parameter `UnitMass_in_g`.
      - `v_unit::Unitful.Velocity=ILLUSTRIS_V_UNIT`: Code parameter `UnitVelocity_in_cm_per_s`.
      - `a0::Float64=1.0`: Cosmological scale factor of the simulation.
      - `h0::Float64=1.0`: Hubble constant as "little h".
    """
    function InternalUnits(;
        l_unit::Unitful.Length=ILLUSTRIS_L_UNIT,
        m_unit::Unitful.Mass=ILLUSTRIS_M_UNIT,
        v_unit::Unitful.Velocity=ILLUSTRIS_V_UNIT,
        a0::Float64=1.0,
        h0::Float64=1.0,
    )

        ############################################################################################
        # Base units
        ############################################################################################

        x_cgs = l_unit * a0 / h0
        x_cosmo = x_cgs |> u"kpc"
        x_comoving = l_unit / h0 |> u"kpc"

        v_cgs = v_unit * sqrt(a0)
        v_cosmo = v_cgs |> u"km*s^-1"

        m_cgs = m_unit / h0
        m_cosmo = m_cgs |> u"Msun"

        ############################################################################################
        # Derived units
        ############################################################################################

        # Only used in non-cosmological simulations
        t_cgs = x_cgs / v_cgs
        t_cosmo = t_cgs |> u"Myr"

        U_cgs = v_unit^2 |> u"erg*g^-1"

        rho_cgs = m_cgs * x_cgs^-3

        # Thermal pressure (it uses v_unit^2 instead of v_cgs^2, which would add an extra factor of a0)
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
            U_cgs,
            rho_cgs,
            P_Pa,
        )

    end

end

"""
Linear grid (1D).

# Fields

  - `grid::Vector{<:Number}`: Vector with the center value of each bin.
  - `ticks::Vector{<:Number}`: Vector with the edges of the bins.
  - `bin_widths::Vector{<:Number}`: Widths of the bins.
  - `log::Bool`: If the grid is logarithmic.
"""
struct LinearGrid
    grid::Vector{<:Number}
    ticks::Vector{<:Number}
    bin_widths::Vector{<:Number}
    log::Bool

    """
        LinearGrid(start::Number, stop::Number, n_bins::Int; <keyword arguments>)

    Constructor for `LinearGrid`.

    # Arguments

      - `start::Number`: Initial value of the grid.
      - `stop::Number`: Final value of the grid.
      - `n_bins::Int`: Number of bins.
      - `log::Bool=false`: If the bins will be logarithmic.
    """
    function LinearGrid(start::Number, stop::Number, n_bins::Int; log::Bool=false)

        # Example of a grid with 3 bins:
        #
        # +-----------------------------+-----------------------------+-----------------------------+
        # |   grid[1] -->|              |   grid[2] -->|              |   grid[3] -->|              |
        # +-----------------------------+-----------------------------+-----------------------------+
        #
        # +-----------------------------+-----------------------------+-----------------------------+
        # |<-- thick[1] (= `start`)     |<-- thick[2]     thick[3] -->|      thick[4] (= `stop`) -->|
        # +-----------------------------+-----------------------------+-----------------------------+

        (
            stop > start ||
            throw(ArgumentError("LinearGrid: `stop` must be larger than `start`, \
            but I got stop = $(stop) <= start = $(start)"))
        )

        if log
            (
                isPositive(start) ||
                throw(ArgumentError("LinearGrid: For a logarithmic grid you need a \
                positive `start`, but I got start = $(start)"))
            )

            # Length unit
            u_l = unit(start)

            log_start = log10(ustrip(start))
            log_stop = log10(ustrip(u_l, stop))

            width = (log_stop - log_start) / n_bins
            grid = [exp10((i - 0.5) * width + log_start) * u_l for i in 1:n_bins]
            ticks = [exp10(i * width + log_start) * u_l for i in 0:n_bins]
            bin_widths = [ticks[i + 1] - ticks[i] for i in 1:n_bins]
        else
            width = (stop - start) / n_bins
            grid = [(i - 0.5) * width + start for i in 1:n_bins]
            ticks = [i * width + start for i in 0:n_bins]
            bin_widths = [width for _ in 1:n_bins]
        end

        new(grid, ticks, bin_widths, log)

    end
end

"""
Square grid (2D).

# Fields

  - `grid::Matrix{NTuple{2,<:Number}}`: Matrix with the physical coordinates of the center of each pixel in the grid.
  - `x_ticks::Vector{<:Number}`: Full set of possible values for the x coordinate.
  - `y_ticks::Vector{<:Number}`: Full set of possible values for the y coordinate.
  - `size::Number`: Side length of the square grid.
  - `n_bins::Int`: Number of bins per side of the grid.
  - `bin_width::Number`: Side length of each bin.
  - `bin_area::Number`: Area of each bin.
"""
struct SquareGrid
    grid::Matrix{NTuple{2,<:Number}}
    x_ticks::Vector{<:Number}
    y_ticks::Vector{<:Number}
    size::Number
    n_bins::Int
    bin_width::Number
    bin_area::Number

    """
        SquareGrid(
            size::Number,
            n_bins::Int;
            <keyword arguments>
        )

    Constructor for `SquareGrid`.

    # Arguments

      - `size::Number`: Side length of the square grid.
      - `n_bins::Int`: Number of bins per dimesion of the grid.
      - `center::Vector{<:Number}=zeros(typeof(size), 3)`: 3D location of the center of the grid. The z axis is taken as the normal vector of the grid.
    """
    function SquareGrid(
        size::Number,
        n_bins::Int;
        center::Vector{<:Number}=zeros(typeof(size), 3),
    )

        # Example of a 3x3 grid:
        #
        # x_ticks = [0, 1, 2]
        # y_ticks = [0, 1, 2]
        #
        # +------------------+------------------+------------------+
        # |      i = 1       |      i = 4       |      i = 7       |
        # | grid[1] = (0, 2) | grid[4] = (1, 2) | grid[7] = (2, 2) |
        # +------------------+------------------+------------------+
        # |      i = 2       |      i = 5       |      i = 8       |
        # | grid[2] = (0, 1) | grid[5] = (1, 1) | grid[8] = (2, 1) |
        # +------------------+------------------+------------------+
        # |      i = 3       |      i = 6       |      i = 9       |
        # | grid[3] = (0, 0) | grid[6] = (1, 2) | grid[9] = (2, 2) |
        # +------------------+------------------+------------------+

        # Compute the bin dimensions
        bin_width = size / n_bins
        bin_area = bin_width * bin_width

        # Compute the x and y coordinates of the center of each square bin
        shift = 0.5 * (size - bin_width)
        x_ticks = [(i - 1) * bin_width - shift + center[1] for i in 1:n_bins]
        y_ticks = [(i - 1) * bin_width - shift + center[2] for i in 1:n_bins]

        # Allocate memory
        grid = Matrix{NTuple{2,<:Number}}(undef, n_bins, n_bins)

        # Compute the position of each grid point
        @inbounds for i in eachindex(grid)
            # The grid index `i` goes from top to bottom first, and then left to right,
            # starting at the top left of the grid
            i_x = ceil(Int, i / n_bins)
            i_y = mod1(i, n_bins) - 1
            # The coordinates are cartesian, so `y` goes from bottom to top and `x` goes from left to right,
            # starting at the bottom left of the grid
            grid[i] = (x_ticks[i_x], y_ticks[end - i_y])
        end

        new(grid, x_ticks, y_ticks, size, n_bins, bin_width, bin_area)

    end
end

"""
Circular grid (2D or 3D), formed by a series of concentric rings or spherical shells.

# Fields

  - `grid::Vector{<:Number}`: Vector with the distance of each bin to the center of the grid.
  - `ticks::Vector{<:Number}`: Vector with the edges of the bins.
  - `center::Vector{<:Number}`: 3D location of the center of the grid. For the 2D grid, the grid is assumed to be in the xy plane.
  - `bin_area::Vector{<:Number}`: Area of each ring.
  - `bin_volumes::Vector{<:Number}`: Volume of each spherical shell.
  - `log::Bool`: If the grid is logarithmic.
"""
struct CircularGrid
    grid::Vector{<:Number}
    ticks::Vector{<:Number}
    center::Vector{<:Number}
    bin_areas::Vector{<:Number}
    bin_volumes::Vector{<:Number}
    log::Bool

    """
        CircularGrid(
            radius::Number,
            n_bins::Int;
            <keyword arguments>
        )

    Constructor for `CircularGrid`.

    # Arguments

      - `radius::Number`: Radius of the grid (equal to the last bin tick).
      - `n_bins::Int`: Number of bins.
      - `center::Vector{<:Number}=zeros(typeof(radius), 3)`: 3D location of the center of the grid. For the 2D grid, the grid is assumed to be in the xy plane.
      - `log::Bool=false`: If the bins will be logarithmic.
      - `shift::Number=zero(radius)`: Distance of the first bin tick to the center.
    """
    function CircularGrid(
        radius::Number,
        n_bins::Int;
        center::Vector{<:Number}=zeros(typeof(radius), 3),
        log::Bool=false,
        shift::Number=zero(radius),
    )

        # Example of a grid with 3 bins:
        #
        # |<-----------------------------------------radius---------------------------------------->|
        #
        # +-----------+-------------------------+-------------------------+-------------------------+
        # |<--shift-->| grid[1] -->|            | grid[2] -->|            | grid[3] -->|            |
        # +-----------+-------------------------+-------------------------+-------------------------+
        #
        # +-----------+-------------------------+-------------------------+-------------------------+
        # |<--shift-->|<-- thick[1] thick[2] -->|             thick[3] -->|             thick[4] -->|
        # +-----------+-------------------------+-------------------------+-------------------------+

        (
            isPositive(radius) ||
            throw(ArgumentError("CircularGrid: `radius` must be positive, \
            but I got radius = $(radius)"))
        )

        if log

            (
                isPositive(shift) ||
                throw(ArgumentError("CircularGrid: For a logarithmic grid you need a \
                positive `shift`, but I got shift = $(shift)"))
            )

            # Length unit
            u_l = unit(radius)

            log_shift = log10(ustrip(u_l, shift))
            log_radius = log10(ustrip(radius))

            width = (log_radius - log_shift) / n_bins
            grid = [exp10((i - 0.5) * width + log_shift) * u_l for i in 1:n_bins]
            ticks = [exp10(i * width + log_shift) * u_l for i in 0:n_bins]
        else
            width = (radius - shift) / n_bins
            grid = [(i - 0.5) * width + shift for i in 1:n_bins]
            ticks = [i * width + shift for i in 0:n_bins]
        end

        bin_areas = [area(ticks[i + 1]) - area(ticks[i]) for i in 1:n_bins]
        bin_volumes = [volume(ticks[i + 1]) - volume(ticks[i]) for i in 1:n_bins]

        new(grid, ticks, center, bin_areas, bin_volumes, log)

    end
end

"""
Plotting parameters for a quantity.

# Fields

  - `request::Dict{Symbol,Vector{String}} = Dict{Symbol,Vector{String}}()`: Data request for [`readSnapshot`](@ref). It must have the shape `cell/particle type` -> [`block`, `block`, `block`, ...].
  - `var_name::AbstractString = ""`: Name of the quantity for the plot axis. It should not include units or scaling factors.
  - `exp_factor::Int = 0`: Numerical exponent to scale down the axis, e.g. if `x_exp_factor` = 10 the values will be divided by ``10^{10}``. The default is no scaling.
  - `unit::Unitful.Units = Unitful.NoUnits`: Target unit for the axis.
  - `axis_label::AbstractString = "auto_label"`: Label for the axis. It can contain the string `auto_label`, which will be replaced by the default label: `var_name` / 10^`exp_factor` `unit`.
"""
@kwdef mutable struct PlotParams
    request::Dict{Symbol,Vector{String}} = Dict{Symbol,Vector{String}}()
    var_name::AbstractString = ""
    exp_factor::Int = 0
    unit::Unitful.Units = Unitful.NoUnits
    axis_label::AbstractString = "auto_label"
end

# Quantities in the simulation

"""
Dimensional information about a physical quantity.

# Fields

  - `hdf5_name::String`: HDF5 block name.
  - `dimensions::Unitful.Dimensions`: Physical dimensions of the quantity, e.g. `Unitful.𝐋 * Unitful.𝐓^-1`.
  - `unit::Union{Unitful.Units,Symbol}`: Units of the quantity within the code. It can be a unit from [Unitful](https://github.com/PainterQubits/Unitful.jl) or [UnitfulAstro](https://github.com/JuliaAstro/UnitfulAstro.jl), or it can be the symbol `:internal` which denotes internal code units.
"""
struct Qty
    hdf5_name::String
    dimensions::Unitful.Dimensions
    unit::Union{Unitful.Units,Symbol}
end

"""
Dictionary of dimensional properties for the quantities in the code.
"""
const QUANTITIES = Dict(
    # Snapshot quantities
    "ACCE" => Qty("Acceleration", Unitful.𝐋 * Unitful.𝐓^-2, :internal),
    "AGE " => Qty("StellarFormationTime", Unitful.𝐓, :internal),
    "AREA" => Qty("SurfaceArea", Unitful.𝐋^2, :internal),
    "CLKT" => Qty("", Unitful.𝐓, :internal),
    "CMCE" => Qty("CenterOfMass", Unitful.𝐋, :internal),
    "CSND" => Qty("SoundSpeed", Unitful.𝐋 * Unitful.𝐓^-1, :internal),
    "GAGE" => Qty("GFM_StellarFormationTime", Unitful.𝐓, :internal),
    "GIMA" => Qty("GFM_InitialMass", Unitful.𝐌, :internal),
    "GME2" => Qty("GFM_Metals", Unitful.NoDims, Unitful.NoUnits),
    "GMET" => Qty("GFM_Metals", Unitful.NoDims, Unitful.NoUnits),
    "GZ  " => Qty("GFM_Metallicity", Unitful.NoDims, Unitful.NoUnits),
    "GZ2 " => Qty("GFM_Metallicity", Unitful.NoDims, Unitful.NoUnits),
    "Hepp" => Qty("DoublyIonizedHeliumAbundance", Unitful.NoDims, Unitful.NoUnits),
    "ID  " => Qty("ParticleIDs", Unitful.NoDims, Unitful.NoUnits),
    "PAID" => Qty("ParentID", Unitful.NoDims, Unitful.NoUnits),
    "TRID" => Qty("TracerID", Unitful.NoDims, Unitful.NoUnits),
    "MACH" => Qty("MachNumber", Unitful.NoDims, Unitful.NoUnits),
    "MASS" => Qty("Masses", Unitful.𝐌, :internal),
    "NE  " => Qty("ElectronAbundance", Unitful.NoDims, Unitful.NoUnits),
    "NFAC" => Qty("NumFacesCell", Unitful.NoDims, Unitful.NoUnits),
    "NH  " => Qty("NeutralHydrogenAbundance", Unitful.NoDims, Unitful.NoUnits),
    "NHP " => Qty("IonizedHydrogenAbundance", Unitful.NoDims, Unitful.NoUnits),
    "NHeO" => Qty("NeutralHeliumAbundance", Unitful.NoDims, Unitful.NoUnits),
    "NHep" => Qty("SinglyIonizedHeliumAbundance", Unitful.NoDims, Unitful.NoUnits),
    "POS " => Qty("Coordinates", Unitful.𝐋, :internal),
    "PRES" => Qty("Pressure", Unitful.𝐌 * Unitful.𝐋^-1 * Unitful.𝐓^-2, :internal),
    "RHO " => Qty("Density", Unitful.𝐌 * Unitful.𝐋^-3, :internal),
    "SFC1" => Qty("", Unitful.𝐓, :internal),
    "SFC2" => Qty("", Unitful.𝐌, :internal),
    "SFC3" => Qty("", Unitful.𝐌 * Unitful.𝐓^-1, u"Msun*yr^-1"),
    "SFC4" => Qty("", Unitful.𝐌 * Unitful.𝐓^-1, u"Msun*yr^-1"),
    "SFC5" => Qty("", Unitful.𝐌, :internal),
    "SFC6" => Qty("", Unitful.𝐌, :internal),
    "SFR " => Qty("StarFormationRate", Unitful.𝐌 * Unitful.𝐓^-1, u"Msun*yr^-1"),
    "SOFT" => Qty("Softenings", Unitful.𝐋, :internal),
    "TBH"  => Qty("TimebinHydro", Unitful.NoDims, Unitful.NoUnits),
    "TEMP" => Qty("Temperature", Unitful.𝚯, u"K"),
    "TSTP" => Qty("TimeStep", Unitful.𝐓, :internal),
    "U   " => Qty("InternalEnergy", Unitful.𝐋^2 * Unitful.𝐓^-2, :internal),
    "VEL " => Qty("Velocities", Unitful.𝐋 * Unitful.𝐓^-1, :internal),
    "VEVE" => Qty("VertexVelocity", Unitful.𝐋 * Unitful.𝐓^-1, :internal),
    "VOL " => Qty("Volume", Unitful.𝐋^3, :internal),
    "Z   " => Qty("Metallicity", Unitful.NoDims, Unitful.NoUnits),
    "POS " => Qty("Coordinates", Unitful.𝐋, :internal),
    "VEL " => Qty("Velocities", Unitful.𝐋 * Unitful.𝐓^-1, :internal),
    "ODIT" => Qty("ODE_IntegrationTime", Unitful.𝐓, u"Myr"),
    "ACIT" => Qty("ODE_AccumulatedIntegrationTime", Unitful.𝐓, u"Myr"),
    "CTIM" => Qty("ODE_CurrentTime", Unitful.𝐓, u"Myr"),
    "DTIM" => Qty("ODE_DeltaTime", Unitful.𝐓, u"Myr"),
    "TAUS" => Qty("ODE_TauS", Unitful.𝐓, u"Myr"),
    "RHOC" => Qty("ODE_ParameterCellDensity", Unitful.𝐌 * Unitful.𝐋^-3, u"cm^-3"),
    "PARZ" => Qty("ODE_ParameterMetallicity", Unitful.NoDims, Unitful.NoUnits),
    "ETAD" => Qty("ODE_ParameterEtaD", Unitful.NoDims, Unitful.NoUnits),
    "ETAI" => Qty("ODE_ParameterEtaI", Unitful.NoDims, Unitful.NoUnits),
    "PARR" => Qty("ODE_ParameterR", Unitful.NoDims, Unitful.NoUnits),
    "FRAC" => Qty("ODE_Fractions", Unitful.NoDims, Unitful.NoUnits),
    # Halo (FoF group) quantities
    "G_BHMass"             => Qty("GroupBHMass", Unitful.𝐌, :internal),
    "G_CM"                 => Qty("GroupCM", Unitful.𝐋, :internal),
    "G_FirstSub"           => Qty("GroupFirstSub", Unitful.NoDims, Unitful.NoUnits),
    "G_GasMetalFractions"  => Qty("GroupGasMetalFractions", Unitful.NoDims, Unitful.NoUnits),
    "G_GasMetallicity"     => Qty("GroupGasMetallicity", Unitful.NoDims, Unitful.NoUnits),
    "G_Len"                => Qty("GroupLen", Unitful.NoDims, Unitful.NoUnits),
    "G_LenType"            => Qty("GroupLenType", Unitful.NoDims, Unitful.NoUnits),
    "G_Mass"               => Qty("GroupMass", Unitful.𝐌, :internal),
    "G_MassType"           => Qty("GroupMassType", Unitful.𝐌, :internal),
    "G_Nsubs"              => Qty("GroupNsubs", Unitful.NoDims, Unitful.NoUnits),
    "G_Pos"                => Qty("GroupPos", Unitful.𝐋, :internal),
    "G_SFR"                => Qty("GroupSFR", Unitful.𝐌 * Unitful.𝐓^-1, u"Msun*yr^-1"),
    "G_StarMetalFractions" => Qty("GroupStarMetalFractions", Unitful.NoDims, Unitful.NoUnits),
    "G_StarMetallicity"    => Qty("GroupStarMetallicity", Unitful.NoDims, Unitful.NoUnits),
    "G_WindMass"           => Qty("GroupWindMass", Unitful.𝐌, :internal),
    "G_M_Crit200"          => Qty("Group_M_Crit200", Unitful.𝐌, :internal),
    "G_M_Crit500"          => Qty("Group_M_Crit500", Unitful.𝐌, :internal),
    "G_M_Mean200"          => Qty("Group_M_Mean200", Unitful.𝐌, :internal),
    "G_M_TopHat200"        => Qty("Group_M_TopHat200", Unitful.𝐌, :internal),
    "G_R_Crit200"          => Qty("Group_R_Crit200", Unitful.𝐋, :internal),
    "G_R_Crit500"          => Qty("Group_R_Crit500", Unitful.𝐋, :internal),
    "G_R_Mean200"          => Qty("Group_R_Mean200", Unitful.𝐋, :internal),
    "G_R_TopHat200"        => Qty("Group_R_TopHat200", Unitful.𝐋, :internal),
    "G_Vel"                => Qty("GroupVel", Unitful.𝐋 * Unitful.𝐓^-1, :gvel),
    # Subhalo (subfind) quantities
    "S_BHMass"                       => Qty("SubhaloBHMass", Unitful.𝐌, :internal),
    "S_CM"                           => Qty("SubhaloCM", Unitful.𝐋, :internal),
    "S_GasMetalFractions"            => Qty("SubhaloGasMetalFractions", Unitful.NoDims, Unitful.NoUnits),
    "S_GasMetalFractionsHalfRad"     => Qty("SubhaloGasMetalFractionsHalfRad", Unitful.NoDims, Unitful.NoUnits),
    "S_GasMetalFractionsMaxRad"      => Qty("SubhaloGasMetalFractionsMaxRad", Unitful.NoDims, Unitful.NoUnits),
    "S_GasMetalFractionsSfr"         => Qty("SubhaloGasMetalFractionsSfr", Unitful.NoDims, Unitful.NoUnits),
    "S_GasMetalFractionsSfrWeighted" => Qty("SubhaloGasMetalFractionsSfrWeighted", Unitful.NoDims, Unitful.NoUnits),
    "S_GasMetallicity"               => Qty("SubhaloGasMetallicity", Unitful.NoDims, Unitful.NoUnits),
    "S_GasMetallicityHalfRad"        => Qty("SubhaloGasMetallicityHalfRad", Unitful.NoDims, Unitful.NoUnits),
    "S_GasMetallicityMaxRad"         => Qty("SubhaloGasMetallicityMaxRad", Unitful.NoDims, Unitful.NoUnits),
    "S_GasMetallicitySfr"            => Qty("SubhaloGasMetallicitySfr", Unitful.NoDims, Unitful.NoUnits),
    "S_GasMetallicitySfrWeighted"    => Qty("SubhaloGasMetallicitySfrWeighted", Unitful.NoDims, Unitful.NoUnits),
    "S_GrNr"                         => Qty("SubhaloGrNr", Unitful.NoDims, Unitful.NoUnits),
    "S_HalfmassRad"                  => Qty("SubhaloHalfmassRad", Unitful.𝐋, :internal),
    "S_HalfmassRadType"              => Qty("SubhaloHalfmassRadType", Unitful.𝐋, :internal),
    "S_IDMostbound"                  => Qty("SubhaloIDMostbound", Unitful.NoDims, Unitful.NoUnits),
    "S_Len"                          => Qty("SubhaloLen", Unitful.NoDims, Unitful.NoUnits),
    "S_LenType"                      => Qty("SubhaloLenType", Unitful.NoDims, Unitful.NoUnits),
    "S_Mass"                         => Qty("SubhaloMass", Unitful.𝐌, :internal),
    "S_MassInHalfRad"                => Qty("SubhaloMassInHalfRad", Unitful.𝐌, :internal),
    "S_MassInHalfRadType"            => Qty("SubhaloMassInHalfRadType", Unitful.𝐌, :internal),
    "S_MassInMaxRad"                 => Qty("SubhaloMassInMaxRad", Unitful.𝐌, :internal),
    "S_MassInMaxRadType"             => Qty("SubhaloMassInMaxRadType", Unitful.𝐌, :internal),
    "S_MassInRad"                    => Qty("SubhaloMassInRad", Unitful.𝐌, :internal),
    "S_MassInRadType"                => Qty("SubhaloMassInRadType", Unitful.𝐌, :internal),
    "S_MassType"                     => Qty("SubhaloMassType", Unitful.𝐌, :internal),
    "S_Parent"                       => Qty("SubhaloParent", Unitful.NoDims, Unitful.NoUnits),
    "S_Pos"                          => Qty("SubhaloPos", Unitful.𝐋, :internal),
    "S_SFR"                          => Qty("SubhaloSFR", Unitful.𝐌 * Unitful.𝐓^-1, u"Msun*yr^-1"),
    "S_SFRinHalfRad"                 => Qty("SubhaloSFRinHalfRad", Unitful.𝐌 * Unitful.𝐓^-1, u"Msun*yr^-1"),
    "S_SFRinMaxRad"                  => Qty("SubhaloSFRinMaxRad", Unitful.𝐌 * Unitful.𝐓^-1, u"Msun*yr^-1"),
    "S_SFRinRad"                     => Qty("SubhaloSFRinRad", Unitful.𝐌 * Unitful.𝐓^-1, u"Msun*yr^-1"),
    "S_StarMetalFractions"           => Qty("SubhaloStarMetalFractions", Unitful.NoDims, Unitful.NoUnits),
    "S_StarMetalFractionsHalfRad"    => Qty("SubhaloStarMetalFractionsHalfRad", Unitful.NoDims, Unitful.NoUnits),
    "S_StarMetalFractionsMaxRad"     => Qty("SubhaloStarMetalFractionsMaxRad", Unitful.NoDims, Unitful.NoUnits),
    "S_StarMetallicity"              => Qty("SubhaloStarMetallicity", Unitful.NoDims, Unitful.NoUnits),
    "S_StarMetallicityHalfRad"       => Qty("SubhaloStarMetallicityHalfRad", Unitful.NoDims, Unitful.NoUnits),
    "S_StarMetallicityMaxRad"        => Qty("SubhaloStarMetallicityMaxRad", Unitful.NoDims, Unitful.NoUnits),
    "S_StellarPhotometrics"          => Qty("SubhaloStellarPhotometrics", Unitful.NoDims, Unitful.NoUnits),
    "S_StellarPhotometricsMassInRad" => Qty("SubhaloStellarPhotometricsMassInRad", Unitful.𝐌, :internal),
    "S_StellarPhotometricsRad"       => Qty("SubhaloStellarPhotometricsRad", Unitful.𝐋, :internal),
    "S_Vel"                          => Qty("SubhaloVel", Unitful.𝐋 * Unitful.𝐓^-1, u"km*s^-1"),
    "S_VelDisp"                      => Qty("SubhaloVelDisp", Unitful.𝐋 * Unitful.𝐓^-1, u"km*s^-1"),
    "S_Vmax"                         => Qty("SubhaloVmax", Unitful.𝐋 * Unitful.𝐓^-1, u"km*s^-1"),
    "S_VmaxRad"                      => Qty("SubhaloVmaxRad", Unitful.𝐋, :internal),
    "S_WindMass"                     => Qty("SubhaloWindMass", Unitful.𝐌, :internal),
)

"""
Dictionary with the subhalo numbers for the MW and M31 in Hestia simulations.
"""
const HESTIA_SUBHALOS = Dict(
    "Hestia17-11" => Dict(
        :subhalo_number_MW  => 1,
        :subhalo_number_M31 => 0,
    ),
    "Hestia09-18" => Dict(
        :subhalo_number_MW  => 3911,
        :subhalo_number_M31 => 2608,
    ),
    "Hestia37-11" => Dict(
        :subhalo_number_MW  => 920,
        :subhalo_number_M31 => 0,
    ),
)
