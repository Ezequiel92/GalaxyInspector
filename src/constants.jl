####################################################################################################
# Constants
####################################################################################################

# Type aliases

"Union of color types."
const ColorInput = Union{Symbol,ColorTypes.RGB}

"Union of line style types."
const LineStyleInput = Union{Symbol,String,Tuple{String,Symbol},Nothing}

"Union of real numbers and physical quantities types."
const RealOrQty = Union{Real,Unitful.Quantity}

"Union of index types."
const IndexType = Union{Int64,UnitRange{Int64},Vector{Int64},Colon}

# Dimensions of energy per unit mass
@derived_dimension SpecificEnergy Unitful.𝐋^2 * Unitful.𝐓^-2

# New types

"""
Structure with information about a GADGET simulation.

# Fields
- `path::String`: Full path to the directory containing the snapshot files.
- `sim_index::Int32`: Index associated with the simulation.
- `idx::IndexType`: Indexing of the simulation, i.e. which snapshots will be read. It can be an 
  integer (a single snapshot), a vector of integers (several snapshots), an UnitRange (e.g. 5:13) 
  or : (every snapshot will be read).
- `base_name::String`: Base name of the snapshot files, set in the GADGET variable 
  `SnapshotFileBase`.
- `header::GIO.SnapshotHeader`: Header of the snapshot.
- `filter_function::Union{Function,Nothing} = nothing`: A function with the signature: 

  `foo(file_path::String)::Vector{Int64}`
  
  It indicates which particles will be read, taking the file path to a snapshot and returning the 
  list of indices of the selected particles. If set to `nothing`, then no particles are filtered. 
  See the [GadgetIO.jl](https://ludwigboess.github.io/GadgetIO.jl/stable/read_snapshots/#Filter-functions) documentation for examples.
- `a0::Float64`: Scale factor of the initial snapshot (only relevant for cosmological simulations).
"""
struct SimData
    path::String
    sim_index::Int32
    idx::IndexType
    base_name::String
    header::GIO.SnapshotHeader
    filter_function::Union{Function,Nothing}
    a0::Float64
end

"""
Structure with information about a GADGET snapshot.

# Fields
- `path::String`: Full path to the snapshot.
- `global_index::Int32`: Index of the snapshot in the context of the whole simulation.
- `local_index::Int32`: Index of the snapshot in the context of the current running loop.
- `time_stamp::Unitful.Time`: Clock time associated with the snapshot.
"""
struct SnapData
    path::String
    global_index::Int32
    local_index::Int32
    time_stamp::Unitful.Time
end

"""
Structure with the dimensional information of a GADGET quantity.

# Fields
- `dimension::Unitful.Dimensions`: Physical dimensions of the quantity, e.g. `Unitful.𝐋 / Unitful.𝐓`.
- `unit::Union{Symbol,Unitful.Units}`: Original units of the quantity. It can be a unit object from 
  [Unitful](https://github.com/PainterQubits/Unitful.jl) or [UnitfulAstro](https://github.com/JuliaAstro/UnitfulAstro.jl), or it can be the symbol `:internal` which denotes that GADGET internal 
  units are being used.
"""
struct Qty
    dimension::Unitful.Dimensions
    unit::Union{Symbol,Unitful.Units}
end

# Numeric constants

@doc raw"""
Critical density for Newtonian simulations, above which the gas particles enter the star formation 
routine.

This value correspond to ``\mathrm{CritPhysDensity} = 0.318 \, \mathrm{cm}^{-3}`` in `param.txt`. Which is 
converted to internal units withing GADGET3 with

CritPhysDensity * PROTONMASS / HYDROGEN\_MASSFRAC / UnitDensity\_in\_cgs
"""
const CRITICAL_DENSITY = 1.033780605417362e7 * UnitfulAstro.Msun * UnitfulAstro.kpc^-3

@doc raw"""
``H_0 = 0.102201 \, \mathrm{Gyr}^{-1} = 100 \, \mathrm{km} \, \mathrm{s}^{-1} \, \mathrm{Mpc}^{-1}``
"""
const HUBBLE_CONST = 0.102201

"""
Solar metallicity.

M. Asplund et al. (2009). *The Chemical Composition of the Sun.* Annual Review of Astronomy and
Astrophysics, **47(1)**, 481–522. [https://doi.org/10.1146/annurev.astro.46.060407.145222](https://doi.org/10.1146/annurev.astro.46.060407.145222)
"""
const SOLAR_METALLICITY = 0.0134

"""
The slope for the Kennicutt-Schmidt law.

R. C. Kennicutt (1998). *The Global Schmidt Law in Star-forming Galaxies.* The Astrophysical Journal, 
**498(2)**, 541-552. [https://doi.org/10.1086/305588](https://doi.org/10.1086/305588)
"""
const KS98_SLOPE = 1.4

"""
Intercept for the Kennicutt-Schmidt law.

R. C. Kennicutt (1998). *The Global Schmidt Law in Star-forming Galaxies.* The Astrophysical Journal, 
**498(2)**, 541-552. [https://doi.org/10.1086/305588](https://doi.org/10.1086/305588)
"""
const KS98_INTERCEPT = 2.5e-4 * (UnitfulAstro.Msun / UnitfulAstro.yr / UnitfulAstro.kpc^2)

"""
Unit of area density for the Kennicutt-Schmidt law.

R. C. Kennicutt (1998). *The Global Schmidt Law in Star-forming Galaxies.* The Astrophysical Journal, 
**498(2)**, 541-552. [https://doi.org/10.1086/305588](https://doi.org/10.1086/305588)
"""
const KS98_RHO_UNIT = UnitfulAstro.Msun * UnitfulAstro.pc^-2

# List-like constants

"""
GADGET index for the types of particles that can be simulated.

See GADGET2 [User's Guide](https://wwwmpa.mpa-garching.mpg.de/gadget/users-guide.pdf).

Note: In most simulations only `:gas`, `:halo` (dark matter) and `:stars` are used.
"""
const ParticleType = Dict(
    :gas      => 0,
    :halo     => 1,
    :disk     => 2,
    :bulge    => 3,
    :stars    => 4,
    :boundary => 5,
)

"Database of the HDF5 names for the blocks."
const HDF5Names = Dict(
    "HSML" => "SmoothingLength",
    "POS"  => "Coordinates",
    "VEL"  => "Velocities",
    "ID"   => "ParticleIDs",
    "MASS" => "Masses",
    "POT"  => "Potential",
    "U"    => "InternalEnergy",
    "RHO"  => "Density",
    "PRES" => "Pressure",
    "CSND" => "SoundSpeed",
    "NE"   => "ElectronAbundance",
    "NH"   => "NeutralHydrogenAbundance",
    "SFR"  => "StarFormationRate",
    "DIVV" => "VelocityDivergence",
    "ROTV" => "VelocityCurl",
    "COHE" => "CoolingHeatingEnergy",
    "AREA" => "SurfaceArea",
    "NFAC" => "NumFacesCell",
    "COOR" => "CoolingRate",
    "VORT" => "Vorticity",
    "GRAP" => "PressureGradient",
    "GRAR" => "DensityGradient",
    "GRAV" => "VelocityGradient",
    "GRAB" => "BfieldGradient",
    "VOL"  => "Volume",
    "VEVE" => "VertexVelocity",
    "FACA" => "MaxFaceAngle",
    "CMCE" => "CenterOfMass",
    "TASK" => "task",
    "TBH"  => "TimebinHydro",
    "TSTP" => "TimeStep",
    "ACCE" => "Acceleration",
    "SOFT" => "Softenings",
    "GINT" => "GravityInteractions",
    "BFLD" => "MagneticField",
    "DIVB" => "MagneticFieldDivergence",
    "PASS" => "PassiveScalars",
    "SFDE" => "SubfindDensity",
    "SFDD" => "SubfindDMDensity",
    "SFHS" => "SubfindHsml",
    "SFVD" => "SubfindVelDisp",
    "GROU" => "GroupNr",
    "HRGM" => "HighResGasMass",
    "REF"  => "AllowRefinement",
    "IONF" => "IonizedFraction",
    "ATOF" => "AtomicFraction",
    "MOLF" => "MolecularFraction",
    "METF" => "MetalFraction",
    "STAF" => "StellarFraction",
)

"Database of dimensional properties for the quantities in the simulations."
const QUANTITIES = Dict(
    "POS"         => Qty(Unitful.𝐋, :internal),
    "VEL"         => Qty(Unitful.𝐋 / Unitful.𝐓, :internal),
    "ID"          => Qty(Unitful.NoDims, Unitful.NoUnits),
    "MASS"        => Qty(Unitful.𝐌, :internal),
    "U"           => Qty(Unitful.𝐋^2 / Unitful.𝐓^2, :internal),
    "RHO"         => Qty(Unitful.𝐌 / Unitful.𝐋^3, :internal),
    "NE"          => Qty(Unitful.NoDims, Unitful.NoUnits),
    "NH"          => Qty(Unitful.NoDims, Unitful.NoUnits),
    "HSML"        => Qty(Unitful.𝐋, :internal),
    "SFR"         => Qty(Unitful.𝐌 / Unitful.𝐓, UnitfulAstro.Msun / UnitfulAstro.yr),
    "AGE"         => Qty(Unitful.𝐓, :internal),
    "Z"           => Qty(Unitful.𝐌, :internal),
    "FION"        => Qty(Unitful.NoDims, Unitful.NoUnits),
    "FATO"        => Qty(Unitful.NoDims, Unitful.NoUnits),
    "FMOL"        => Qty(Unitful.NoDims, Unitful.NoUnits),
    "FMET"        => Qty(Unitful.NoDims, Unitful.NoUnits),
    "FSTR"        => Qty(Unitful.NoDims, Unitful.NoUnits),
    "IONF"        => Qty(Unitful.NoDims, Unitful.NoUnits),
    "ATOF"        => Qty(Unitful.NoDims, Unitful.NoUnits),
    "MOLF"        => Qty(Unitful.NoDims, Unitful.NoUnits),
    "METF"        => Qty(Unitful.NoDims, Unitful.NoUnits),
    "STAF"        => Qty(Unitful.NoDims, Unitful.NoUnits),
    "TSTP"        => Qty(Unitful.𝐓, :internal),
    "SFRTXT_COL1" => Qty(Unitful.𝐓, :internal),
    "SFRTXT_COL2" => Qty(Unitful.𝐌, :internal),
    "SFRTXT_COL3" => Qty(Unitful.𝐌 / Unitful.𝐓, UnitfulAstro.Msun / UnitfulAstro.yr),
    "SFRTXT_COL4" => Qty(Unitful.𝐌 / Unitful.𝐓, UnitfulAstro.Msun / UnitfulAstro.yr),
    "SFRTXT_COL5" => Qty(Unitful.𝐌, :internal),
    "SFRTXT_COL6" => Qty(Unitful.𝐌 / Unitful.𝐓, :internal),
    "CLOCK_TIME"  => Qty(Unitful.𝐓, :internal),
    "TEMP"        => Qty(Unitful.𝚯, Unitful.K),
    "POT"         => Qty(Unitful.𝐋^2 / Unitful.𝐓^2, :internal),
    "PRES"        => Qty(Unitful.𝐌 / Unitful.𝐓 / Unitful.𝐋, :internal),
    "CSND"        => Qty(Unitful.𝐋 / Unitful.𝐓, :internal),
    "AREA"        => Qty(Unitful.𝐋^2, :internal),
    "NFAC"        => Qty(Unitful.NoDims, Unitful.NoUnits),
    "VOL"         => Qty(Unitful.𝐋^3, :internal),
    "VEVE"        => Qty(Unitful.𝐋 / Unitful.𝐓, :internal),
    "CMCE"        => Qty(Unitful.𝐋, :internal),
    "TBH"         => Qty(Unitful.𝐓, :internal),
    "ACCE"        => Qty(Unitful.𝐋 / Unitful.𝐓^2, :internal),
    "IONF"        => Qty(Unitful.NoDims, Unitful.NoUnits),
    "ATOF"        => Qty(Unitful.NoDims, Unitful.NoUnits),
    "MOLF"        => Qty(Unitful.NoDims, Unitful.NoUnits),
    "METF"        => Qty(Unitful.NoDims, Unitful.NoUnits),
    "STAF"        => Qty(Unitful.NoDims, Unitful.NoUnits),
)

"List of markers for scatter plots."
const MARKERS = (
    :circle,
    :pentagon,
    :star4,
    :+,
    :x,
    :rect,
    :ltriangle,
    :dtriangle,
    :utriangle,
    :star5,
    :rtriangle,
)

"List of line styles for line plots."
const LINE_STYLES = (
    nothing,
    :dash,
    :dashdot,
    :dot,
    :dashdotdot,
    "-....",
    (".", :dense),
    ("-", :loose),
)

"Global plot theme."
const theme = Theme(
    backgroundcolor=:gray85,
    Axis=(
        alignmode=Outside(10),
        xlabelpadding=10.0,
        ylabelpadding=10.0,
        xlabelsize=24,
        ylabelsize=24,
        xticklabelsize=20,
        yticklabelsize=20,
        xminorgridvisible=true,
        xminorticksvisible=true,
        yminorgridvisible=true,
        yminorticksvisible=true,
        xminortickalign=1,
        yminortickalign=1,
        xtickalign=1,
        ytickalign=1,
        titlesize=25,
        backgroundcolor=:grey95,
    ),
    Legend=(
        orientation=:horizontal,
        labelsize=22,
        bgcolor=:gray85,
        framevisible=false,
        tellheight=true,
        tellwidth=false,
        nbanks=2,
        margin=(0, 0, -10, -15),
        patchsize=(35, 35),
        patchlabelgap=8,
        linewidth=4,
        markersize=20,
    ),
    Lines=(linewidth=3,),
    Scatter=(markersize=5,),
    Hist=(strokewidth=1, strokecolor=:black, patchcolor=:red),
    ScatterLines=(linewidth=3, markersize=8),
    Heatmap=(colormap=:inferno,),
)
