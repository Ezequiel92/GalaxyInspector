####################################################################################################
# Constants and data structures
####################################################################################################

################
# Configuration
################

"""
If physical units (lengths) will be used, instead of comoving units (lengths).
"""
PHYSICAL_UNITS = false

"""
If logging messages will be printed out.
"""
const logging = Ref(false)

"""
Max grid size to be stored on memory, above this value memory-mapping will be used.
"""
const MMAP_THRESHOLD = 200^3

########################
# Characteristic scales
########################

"""
Characteristic radius.
"""
const DISK_R = 40.0u"kpc"

"""
Characteristic height.
"""
const DISK_HEIGHT = 5.0u"kpc"

"""
Characteristic box size.
"""
const BOX_L = 65.0u"kpc"

"""
Characteristic stellar age limit for the SFR and sSFR.
"""
const AGE_RESOLUTION = 200.0u"Myr"

#######################################
# Reference values from the literature
#######################################

"""
Internal unit of length used in [IllustrisTNG](https://www.tng-project.org/data/docs/specifications/), equal to ``1.0  \\, \\mathrm{kpc}``.
"""
const ILLUSTRIS_L_UNIT = 3.085678e21u"cm"

"""
Internal unit of mass used in [IllustrisTNG](https://www.tng-project.org/data/docs/specifications/), equal to ``10^{10} \\, \\mathrm{M_\\odot}``.
"""
const ILLUSTRIS_M_UNIT = 1.989e43u"g"

"""
Internal unit of velocity used in [IllustrisTNG](https://www.tng-project.org/data/docs/specifications/), equal to ``1.0 \\, \\mathrm{km \\, s^{-1}}``.
"""
const ILLUSTRIS_V_UNIT = 1.0e5u"cm * s^-1"

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
    :S  => 7.12,  # Sulfur
    :Ca => 6.34,  # Calcium
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
    :S  => 32.06,  # Sulfur
    :Ca => 40.078, # Calcium
    :Fe => 55.845, # Iron
)

@doc raw"""
Kennicutt-Schmidt law fits for molecular and neutral gas, from Bigiel et al. (2008) (Table 2, Average).

Power-law index, N, and $A = \log_{10}(a)$, where $a$ is $\Sigma_\mathrm{SFR}$ at the fiducial gas surface density of $10 \, \mathrm{M_\odot \, pc^{-2}}$ are given.

```math
\Sigma_\mathrm{SFR} = a \left( \frac{\Sigma_\mathrm{HI, H_2, gas}}{10 \, \mathrm{M_\odot \, pc^{-2}}} \right)^{\!N} \, ,
```

# References

F. Bigiel et al. (2008). *THE STAR FORMATION LAW IN NEARBY GALAXIES ON SUB-KPC SCALES*. The Astrophysical Journal, **136(6)**, 2846. [doi:10.1088/0004-6256/136/6/2846](https://doi.org/10.1088/0004-6256/136/6/2846)
"""
const A_BIGIEL2008_MOLECULAR = −2.06 ± 0.17
const N_BIGIEL2008_MOLECULAR = 0.96 ± 0.07
const A_BIGIEL2008_NEUTRAL   = −2.39 ± 0.28
const N_BIGIEL2008_NEUTRAL   = 1.85 ± 0.70

@doc raw"""
Kennicutt-Schmidt law best-fit for molecular gas, from Bigiel et al. (2008) (Section 4.3, Equation 3).

Power-law index, N, and $A = \log_{10}(a)$, where $a$ is $\Sigma_\mathrm{SFR}$ at the fiducial gas surface density of $10 \, \mathrm{M_\odot \, pc^{-2}}$ are given.

```math
\Sigma_\mathrm{SFR} = a \left( \frac{\Sigma_\mathrm{H_2}}{10 \, \mathrm{M_\odot \, pc^{-2}}} \right)^{\!N} \, ,
```

# References

F. Bigiel et al. (2008). *THE STAR FORMATION LAW IN NEARBY GALAXIES ON SUB-KPC SCALES*. The Astrophysical Journal, **136(6)**, 2846. [doi:10.1088/0004-6256/136/6/2846](https://doi.org/10.1088/0004-6256/136/6/2846)
"""
const A_BIGIEL2008_BF_MOLECULAR = −2.1 ± 0.2
const N_BIGIEL2008_BF_MOLECULAR = 1.0 ± 0.2

@doc raw"""
Slope of the Kennicutt-Schmidt law, taken from Kennicutt (1998) (Section 4, Equation 4).

```math
\Sigma_\mathrm{SFR} = a \left( \frac{\Sigma_\mathrm{gas}}{1 \, \mathrm{M_\odot \, pc^{-2}}} \right)^{\!N} \mathrm{M_\odot \, yr^{-1} \, kpc^{-2}} \, ,
```

# References

R. C. Kennicutt (1998). *The Global Schmidt Law in Star-forming Galaxies*. The Astrophysical Journal, **498(2)**, 541-552. [doi:10.1086/305588](https://doi.org/10.1086/305588)
"""
const N_KS98 = 1.4 ± 0.15

@doc raw"""
Intercept of the Kennicutt-Schmidt law, taken from Kennicutt (1998) (Section 4, Equation 4).

```math
\Sigma_\mathrm{SFR} = a \left( \frac{\Sigma_\mathrm{gas}}{1 \, \mathrm{M_\odot \, pc^{-2}}} \right)^{\!N} \mathrm{M_\odot \, yr^{-1} \, kpc^{-2}} \, ,
```

# References

R. C. Kennicutt (1998). *The Global Schmidt Law in Star-forming Galaxies*. The Astrophysical Journal, **498(2)**, 541-552. [doi:10.1086/305588](https://doi.org/10.1086/305588)
"""
const a_KS98 = 2.5e-4 ± 0.7e-4

"""
Spatial resolution used in Bigiel et al. (2008).

# References

F. Bigiel et al. (2008). *THE STAR FORMATION LAW IN NEARBY GALAXIES ON SUB-KPC SCALES*. The Astrophysical Journal, **136(6)**, 2846. [doi:10.1088/0004-6256/136/6/2846](https://doi.org/10.1088/0004-6256/136/6/2846)
"""
const BIGIEL_PX_SIZE = 750.0u"pc"

"""
Spatial resolution used in Sun et al. (2023).

# References

J. Sun et al. (2023). *Star Formation Laws and Efficiencies across 80 Nearby Galaxies*. The Astrophysical Journal Letters, **945(2)**, L19. [doi:10.3847/2041-8213/acbd9c](https://doi.org/10.3847/2041-8213/acbd9c)
"""
const SUN_PX_SIZE = 1.5u"kpc"

@doc raw"""
Range of values for

```math
\Sigma_\mathrm{SFR} \, [\mathrm{M_\odot \, yr^{-1} \, kpc^{-2}}] \, ,
```
in the seven spirals in Table 1 of Bigiel et al. (2008), with associated molecular data.

The actual values for the SFR density are taken from Table 2 in Bigiel et al. (2010), using only the ones with associated molecular data.

# References

F. Bigiel et al. (2008). *THE STAR FORMATION LAW IN NEARBY GALAXIES ON SUB-KPC SCALES*. The Astrophysical Journal, **136(6)**, 2846. [doi:10.1088/0004-6256/136/6/2846](https://doi.org/10.1088/0004-6256/136/6/2846)

F. Bigiel et al. (2010). *EXTREMELY INEFFICIENT STAR FORMATION IN THE OUTER DISKS OF NEARBY GALAXIES*. The Astrophysical Journal, **140(5)**, 1194. [doi:10.1088/0004-6256/140/5/1194](https://doi.org/10.1088/0004-6256/140/5/1194)
"""
const BIGIEL2008_SFR_RANGE = exp10.([-2.99, -0.33]) .* u"Msun * yr^-1 * kpc^-2"

@doc raw"""
Range of values for

```math
\Sigma_\mathrm{SFR} \, [\mathrm{M_\odot \, yr^{-1} \, kpc^{-2}}] \, ,
```
from the combine data (Table 1 and 2) in Kennicutt (1998).

# References

R. C. Kennicutt (1998). *The Global Schmidt Law in Star-forming Galaxies*. The Astrophysical Journal, **498(2)**, 541-552. [doi:10.1086/305588](https://doi.org/10.1086/305588)
"""
const KS98_SFR_RANGE = exp10.([-3.55, 2.98]) .* u"Msun * yr^-1 * kpc^-2"

"""
Path to Table 2 from Bigiel et al. (2010).

# References

F. Bigiel et al. (2010). *EXTREMELY INEFFICIENT STAR FORMATION IN THE OUTER DISKS OF NEARBY GALAXIES*. The Astrophysical Journal, **140(5)**, 1194. [doi:10.1088/0004-6256/140/5/1194](https://doi.org/10.1088/0004-6256/140/5/1194)
"""
const BIGIEL2010_TABLE_2 = joinpath(@__DIR__, "../../experimental_data/bigiel_2010/table_02.txt")

"""
Path to Table 3 from Bigiel et al. (2010).

# References

F. Bigiel et al. (2010). *EXTREMELY INEFFICIENT STAR FORMATION IN THE OUTER DISKS OF NEARBY GALAXIES*. The Astrophysical Journal, **140(5)**, 1194. [doi:10.1088/0004-6256/140/5/1194](https://doi.org/10.1088/0004-6256/140/5/1194)
"""
const BIGIEL2010_TABLE_3 = joinpath(@__DIR__, "../../experimental_data/bigiel_2010/table_03.txt")

"""
Path to Table A1 from Sun et al. (2023).

# References

J. Sun et al. (2023). *Star Formation Laws and Efficiencies across 80 Nearby Galaxies*. The Astrophysical Journal Letters, **945(2)**, L19. [doi:10.3847/2041-8213/acbd9c](https://doi.org/10.3847/2041-8213/acbd9c)
"""
const SUN2023_TABLE = joinpath(@__DIR__, "../../experimental_data/sun_2023.txt")

"""
Path to the file with the Milky Way profiles from Mollá et al. (2015).

# References

M. Mollá et al. (2015). *Galactic chemical evolution: stellar yields and the initial mass function*. Monthly Notices of the Royal Astronomical Society **451(4)**, 3693–3708. [doi:10.1093/mnras/stv1102](https://doi.org/10.1093/mnras/stv1102)
"""
const MOLLA2015_DATA_PATH = joinpath(@__DIR__, "../../experimental_data/molla_2015.csv")

"""
Path to the file with the global galactic properties from Feldmann (2020).

# References

R. Feldmann (2020). *The link between star formation and gas in nearby galaxies*. Communications Physics **3(226)**. [doi:10.1038/s42005-020-00493-0](https://doi.org/10.1038/s42005-020-00493-0)
"""
const FELDMANN2020_DATA_PATH = joinpath(@__DIR__, "../../experimental_data/feldmann_2020.csv")

"""
Path to the file with the profiles from Leroy et al. (2008).

# References

A. K. Leroy et al. (2008). *THE STAR FORMATION EFFICIENCY IN NEARBY GALAXIES: MEASURING WHERE GAS FORMS STARS EFFECTIVELY*. The Astronomical Journal **136(6)**, 2782–2845. [doi:10.1088/0004-6256/136/6/2782](https://doi.org/10.1088/0004-6256/136/6/2782)

"""
const LEROY2008_DATA_PATH = joinpath(@__DIR__, "../../experimental_data/leroy_2008.jld2")

"""
Path to the file with the fits from McMillan (2011).

# References

P. J. McMillan (2011). *Mass models of the Milky Way*. Monthly Notices of the Royal Astronomical Society **414(3)**, 2446-2457. [doi:10.1111/j.1365-2966.2011.18564.x](https://doi.org/10.1111/j.1365-2966.2011.18564.x)
"""
const MCMILLAN2011_DATA_PATH = joinpath(@__DIR__, "../../experimental_data/mcmillan_2011.jld2")

"""
Path to the file with the stellar magnitudes for the g, r, and i filters (SDSS AB system) from Millán-Irigoyen et al. (2025).

The effective wavelengths of each filter are taken from Table 2 (2.5m reference) of Doi et al. (2010):

g band -> λeff = 4627 Å -> "blue" channel
r band -> λeff = 6140 Å -> "green" channel
i band -> λeff = 7467 Å -> "red" channel

# References

I. Millán-Irigoyen et al. (2025). *HR-pyPopStar II: high spectral resolution evolutionary synthesis models low metallicity expansion and the properties of the stellar populations of dwarf galaxies*. arXiv. [doi:10.48550/arxiv.2510.02886](https://doi.org/10.48550/arxiv.2510.02886)

M. Doi et al. (2010). *PHOTOMETRIC RESPONSE FUNCTIONS OF THE SLOAN DIGITAL SKY SURVEY IMAGER*. The Astronomical Journal **139(4)**, 1628-1648. [doi:10.1088/0004-6256/139/4/1628](https://doi.org/10.1088/0004-6256/139/4/1628)
"""
const MILLANIRIGOYEN2025_DATA_PATH = joinpath(@__DIR__, "../../experimental_data/millan-irigoyen_2025.txt")

@doc raw"""
Extinction factor from T. Güver et al. (2009).

```math
\text{Extinction factor} = \frac{N_\mathrm{H}}{A(V)} = (2.21 \pm 0.09) \times 10^{21} \, \mathrm{atoms \, cm^{-2} \, mag^{-1}} \, .
```

We assume that an atom mass is the proton mass.

T. Güver et al. (2009). *The relation between optical extinction and hydrogen column density in the Galaxy*. Monthly Notices of the Royal Astronomical Society **400(4)**, 2050-2053. [doi:10.1111/j.1365-2966.2009.15598.x](https://doi.org/10.1111/j.1365-2966.2009.15598.x)
"""
const EXTINCTION_FACTOR = 17.7006u"Msun*pc^-2"

@doc raw"""
Relation between the extinction in the g, r, and i bands of SDSS and the extinction in the V band from Table 3 of S. Wang et al. (2019).

```math
\left. A_\lambda / A_V \right\vert_g = 1.205 \pm 0.010 \, .
```

```math
\left. A_\lambda / A_V \right\vert_r = 0.848 \pm 0.006 \, .
```

```math
\left. A_\lambda / A_V \right\vert_i = 0.630 \pm 0.004 \, .
```

S. Wang et al. (2019). *The Optical to Mid-infrared Extinction Law Based on the APOGEE, Gaia DR2, Pan-STARRS1, SDSS, APASS, 2MASS, and WISE Surveys*. The American Astronomical Society **877(2)**, 116. [doi:10.3847/1538-4357/ab1c61](https://doi.org/10.3847/1538-4357/ab1c61)
"""
const AλAV_g = 1.205
const AλAV_r = 0.848
const AλAV_i = 0.63

"""
Reference pressure for the molecular fraction-pressure relation, from Blitz et al. (2006) (Table 2, "Mean" row, Third column).

# References

L. Blitz et al. (2006). *The Role of Pressure in GMC Formation II: The H2-Pressure Relation*. The Astrophysical Journal, **650(2)**, 933. [doi:10.1086/505417](https://doi.org/10.1086/505417)
"""
const P0 = 3.5e4u"K*cm^-3" * Unitful.k

"""
Reference exponent for the molecular fraction-pressure relation, from Blitz et al. (2006) (Table 2, "Mean" row, Second column).

We use -α here.

# References

L. Blitz et al. (2006). *The Role of Pressure in GMC Formation II: The H2-Pressure Relation*. The Astrophysical Journal, **650(2)**, 933. [doi:10.1086/505417](https://doi.org/10.1086/505417)
"""
const ALPHA_BLITZ = -0.92

###############
# Type aliases
###############

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

#########################
# Makie.jl default theme
#########################

"""
Default colors.
"""
const WONG_BLUE    = Makie.wong_colors()[1]
const WONG_ORANGE  = Makie.wong_colors()[2]
const WONG_GREEN   = Makie.wong_colors()[3]
const WONG_PINK    = Makie.wong_colors()[4]
const WONG_CELESTE = Makie.wong_colors()[5]
const WONG_RED     = Makie.wong_colors()[6]
const WONG_YELLOW  = Makie.wong_colors()[7]

"""
Default list of marker types.
"""
const MARKERS = [
    :circle,
    :rect,
    :dtriangle,
    :utriangle,
    :cross,
    :diamond,
    :ltriangle,
    :rtriangle,
    :pentagon,
    :xcross,
    :hexagon,
]

"""
Default list of line styles.
"""
const LINE_STYLES = [:solid, :dash, :dot, :dashdot, :dashdotdot]

"""
Default cycler.
"""
const CYCLE = Cycle([:color, :linestyle, :marker], covary=true)

"""
Default plot theme.

Regarding the graphic units used, we know that ``1 \\, \\mathrm{mm} = 2.83466 \\, \\mathrm{pt}`` and ``1 \\, \\mathrm{in} = 25.4 \\, \\mathrm{mm}``. Then, if we want ``1 \\, \\mathrm{[code\\,\\,]unit} = 0.1 \\, \\mathrm{mm}`` in vector graphics, we have to use `pt_per_unit` = 0.283466.

For pixel images, we control the ppi with `px_per_unit`. A reasonable high ppi is 600. So, using `px_per_unit` = ``2.3622`` we get ``23.622 \\, \\mathrm{px/mm} \\sim 600 \\, \\mathrm{px/in}`` (remember that ``1 \\, \\mathrm{[code\\,\\,]unit} = 0.1 \\, \\mathrm{mm}``).
"""
const DEFAULT_THEME = Theme(
    #################################################################################################
    # Size of the figures in code units
    # For PDFs and SVGs, 880 [code ]unit = 8.8 cm
    # For PNGs, when printed to a size of 1 point = 0.1 mm, one will get a dpi of 600 (23.622 px/mm)
    #################################################################################################
    size=(880, 880),
    ######################################
    # 35 unit * 0.283466 pt/unit ~ 9.9 pt
    ######################################
    fontsize=35,
    #############################
    # (left, right, bottom, top)
    #############################
    figure_padding=(2, 15, 5, 15),
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
        ######################################################################################
        # Aspect ratio of the figures. The options are:
        # nothing: The aspect ratio will be chosen by [Makie](https://docs.makie.org/stable/)
        # AxisAspect(n): The aspect ratio will be given by the number `n` = width / height
        # DataAspect(): The aspect ratio of the data will be used
        ######################################################################################
        aspect=AxisAspect(1),
    ),
    Legend=(
        tellheight=false,
        tellwidth=false,
        framevisible=false,
        colgap=20,
        rowgap=10,
        halign=:right,
        valign=:bottom,
        nbanks=3,
        titlegap=-5,
        labelsize=30,
        linewidth=5,
        markersize=28,
        patchsize=(40, 40),
    ),
    Lines=(linewidth=5, cycle=CYCLE),
    VLines=(linewidth=3, cycle=CYCLE),
    HLines=(linewidth=3, cycle=CYCLE),
    ScatterLines=(
        linewidth=5,
        markersize=22,
        cycle=CYCLE,
    ),
    Scatter=(markersize=22, cycle=CYCLE),
    Band=(alpha=0.5, cycle=CYCLE),
    Errorbars=(whiskerwidth=10,),
    ########################################################################
    # Alternative colormaps:
    # colormap = :nipy_spectral - nan_color = ColorSchemes.nipy_spectral[1]
    # colormap = :cubehelix     - nan_color = ColorSchemes.cubehelix[1]
    # colormap = :lipari        - nan_color = ColorSchemes.lipari[1]
    ########################################################################
    Heatmap=(colormap=:CMRmap, nan_color=ColorSchemes.CMRmap[1]),
    Colorbar=(
        colormap=:CMRmap,
        size=25,
        ticklabelpad=10,
        minorticksvisible=true,
        ticksize=7,
        labelpadding=2,
    ),
    BarPlot=(
        color_over_background=:black,
        color_over_bar=:black,
        flip_labels_at=10,
        strokecolor=:black,
        strokewidth=1,
        dodge_gap=0.04,
    ),
    Arrows2D=(lengthscale=0.015, color=:white, shaftwidth=2, tipwidth=8),
    Hist=(strokecolor=:black, strokewidth=1),
)

#########
# Strucs
#########

"""
Memory-mapped array structure.

# Fields

  - `data :: Array{T,N}`: Array data.
  - `io   :: IOStream`: IO stream to the memory-mapped file.
  - `path :: String`: Path to the memory-mapped file.
"""
mutable struct MmapArray{T,N} <: AbstractArray{T,N}
    data :: Union{Array{T,N}, Nothing}
    io   :: Union{IOStream, Nothing}
    path :: String

    # Array interface
    Base.size(A::MmapArray) = size(A.data)

    Base.getindex(A::MmapArray, I...) = getindex(A.data, I...)

    Base.setindex!(A::MmapArray, v, I...) = setindex!(A.data, v, I...)

    Base.IndexStyle(::Type{<:MmapArray}) = IndexStyle(Array)

    function Base.close(A::MmapArray)
        Mmap.sync!(A.data)
        close(A.io)
        return nothing
    end

    function Base.delete!(A::MmapArray)
        close(A)
        A.data = nothing
        GC.gc(true)
        rm(A.path; force=true)
        return nothing
    end

    """
        MmapArray(
            path :: AbstractString,
                 :: Type{T},
            dims :: NTuple{N,Int},
        ) where {T,N}

    Constructor for `MmapArray`.

    # Arguments

      - `path :: String`: Path to the memory-mapped file.
      - `     :: Type{T}`: Element type of the array.
      - `dims :: NTuple{N,Int}`: Dimensions of the array.
    """
    function MmapArray(
        path :: String,
             :: Type{T},
        dims :: NTuple{N,Int},
    ) where {T,N}

        bytes = prod(dims) * sizeof(T)

        io = open(path, "w+")
        seek(io, bytes - 1)
        write(io, UInt8(0))
        flush(io)

        data = Mmap.mmap(io, Array{T,N}, dims)

        new{T, N}(data, io, path)

    end
end

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

  - `path           :: String`: Full path to the simulation directory.
  - `index          :: Int`: An index associated with the simulation.
  - `slice          :: IndexType`: Slice of the simulation, i.e. which snapshots will be read. It can be an integer (a single snapshot), a vector of integers (several snapshots), an `UnitRange` (e.g. 5:13), an `StepRange` (e.g. 5:2:13) or (:) (all snapshots).
  - `cosmological   :: Bool`: If the simulation is cosmological,

      + `false` -> Newtonian simulation    (`ComovingIntegrationOn` = 0).
      + `true`  -> Cosmological simulation (`ComovingIntegrationOn` = 1).
  - `snapshot_table :: DataFrame`: A dataframe where each row is a snapshot, and the colums are

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
    path           :: String
    index          :: Int
    slice          :: IndexType
    cosmological   :: Bool
    snapshot_table :: DataFrame
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

      - `l_unit :: Unitful.Length=DEFAULT_L_UNIT`: Code parameter `UnitLength_in_cm`.
      - `m_unit :: Unitful.Mass=DEFAULT_M_UNIT`: Code parameter `UnitMass_in_g`.
      - `v_unit :: Unitful.Velocity=DEFAULT_V_UNIT`: Code parameter `UnitVelocity_in_cm_per_s`.
      - `a      :: Float64=1.0`: Cosmological scale factor of the simulation.
      - `h      :: Float64=1.0`: Dimensionless Hubble parameter, "little h".
    """
    function InternalUnits(;
        l_unit :: Unitful.Length=DEFAULT_L_UNIT,
        m_unit :: Unitful.Mass=DEFAULT_M_UNIT,
        v_unit :: Unitful.Velocity=DEFAULT_V_UNIT,
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
            <keyword arguments>,
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
            <keyword arguments>,
        )

    Constructor for `LinearGrid`.

    # Arguments

      - `x_edges :: Vector{<:Number}`: Relative position of (the edges of) the bins, relative to the `origin`.
      - `origin  :: Vector{<:Number}=zeros(eltype(x_edges), 3)`: Physical position of the coordinate origin for the grid (passthrough argument).
      - `log     :: Bool=false`: If the `x_edges` are logarithmic, which means that `x_axis` will mark the center of the bins in logarithmic space instead of linear space.
    """
    function LinearGrid(
        x_edges :: Vector{<:Number};
        origin  :: Vector{<:Number}=zeros(eltype(x_edges), 3),
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
      - `origin :: Vector{<:Number}=zeros(eltype(grid_length), 3)`: Physical position of the center of the grid.
    """
    function SquareGrid(
        size   :: Vector{<:Number},
        n_bins :: NTuple{2,Int};
        origin :: Vector{<:Number}=zeros(eltype(grid_length), 3),
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
      - `origin :: Vector{<:Number}=zeros(eltype(grid_length), 3)`: Physical position of the center of the grid.
    """
    function CubicGrid(
        size   :: Vector{<:Number},
        n_bins :: NTuple{3,Int};
        origin :: Vector{<:Number}=zeros(eltype(grid_length), 3),
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

"""
Plotting parameters for a quantity.

# Fields

  - `request::Dict{Symbol,Vector{String}} = Dict{Symbol,Vector{String}}()`: Data request for [`readSnapshot`](@ref). It must have the shape `cell/particle type` -> [`block name`, `block name`, `block name`, ...].
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

##########################
# Code specific constants
##########################

include("./arepo.jl")

#####################
# Derived quantities
#####################

include("./quantities.jl")

######################
# Cell/particle types
######################

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

##############################
# Default filter dictionaries
##############################

"""
Filter dictionary that does not exclude any cell/particle.
"""
const PASS_ALL = Dict{Symbol,IndexType}(key => (:) for key in keys(PARTICLE_INDEX))

"""
Filter dictionary that excludes every cell/particle.
"""
const PASS_NONE = Dict{Symbol,IndexType}(key => Int[] for key in keys(PARTICLE_INDEX))

###################
# Tracked elements
###################

"""
Shift for the solar abundance in dex.
"""
const ABUNDANCE_SHIFT = Dict(
    :H     => 0.0,  # Hydrogen
    :He    => 0.0,  # Helium
    :C     => 0.0,  # Carbon
    :N     => 0.0,  # Nitrogen
    :O     => 12.0, # Oxygen
    :Ne    => 0.0,  # Neon
    :Mg    => 0.0,  # Magnesium
    :Si    => 0.0,  # Silicon
    :Fe    => 0.0,  # Iron
    :Other => 0.0,  # All other
)

##########################
# Transformations presets
##########################

"""
List of symbols for the transformation presets.

The options have the form :{component}_{group} selecting which cells/particle to consider for the center of mass (new origin for the translation) and principal axis (new reference system for the rotation), where component can be:

      + :all         -> Every component present in data_dict
      + :{component} -> One of the keys of PARTICLE_INDEX

    and group can be:
      + :box     -> Whole simulation box
      + :halo    -> Main halo
      + :subhalo -> Main subhalo
"""
const TRANSFORM_LIST = [
    Symbol(component, :_, group)
    for component in push!(collect(keys(PARTICLE_INDEX)), :all)
    for group in [:box, :halo, :subhalo]
]

"""
Dictionary mapping each transformation to its component and group.
"""
const TRANSFORM_SPLITS = Dict(
    Symbol(component, :_, group) => (component, group)
    for component in push!(collect(keys(PARTICLE_INDEX)), :all)
    for group in [:box, :halo, :subhalo]
)
