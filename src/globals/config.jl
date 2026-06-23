####################################################################################################
# Configuration constants
####################################################################################################

#######
# Mmap
#######

"""
Maximum fraction of the free physical memory that an array can occupy before we start using memory mapping.
"""
const MMAP_MEMORY_FRACTION = Ref{Float64}(0.2)

#########
# Galaxy
#########

"""
Characteristic disk radius.
"""
const DISK_R = Ref(40.0u"kpc")

"""
Characteristic disk height.
"""
const DISK_HEIGHT = Ref(5.0u"kpc")

"""
Characteristic plotting box size.
"""
const BOX_L = Ref(65.0u"kpc")

"""
Characteristic radius for rotations.
"""
const ROTATION_R = Ref(5.0u"kpc")

"""
Characteristic bulge radius (for R25 calculation).
"""
const BULGE_R = Ref(3.0u"kpc")

"""
Characteristic stellar age for the SFR and sSFR.
"""
const AGE_RESOLUTION = Ref(200.0u"Myr")

#############
# Abundances
#############

"""
Shift for the solar abundance of each element in dex.
"""
const ABUNDANCE_SHIFT = Ref{Dict{Symbol,Float64}}(
    Dict(
        :H     => 0.0,  # Hydrogen
        :He    => 0.0,  # Helium
        :C     => 0.0,  # Carbon
        :N     => 0.0,  # Nitrogen
        :O     => 12.0, # Oxygen
        :Ne    => 0.0,  # Neon
        :Mg    => 0.0,  # Magnesium
        :Si    => 0.0,  # Silicon
        :Fe    => 0.0,  # Iron
    )
)

########
# Arepo
########

"""
Base name of the snapshot files, set in the code variable `SnapshotFileBase`.
"""
const SNAP_BASENAME = Ref{String}("snap")

"""
Base name of the group catalog files.
"""
const GC_BASENAME = Ref{String}("fof_subhalo_tab")

"""
Mass fraction of hydrogen.
"""
const HYDROGEN_MASSFRAC = Ref{Float64}(0.76)

"""
Adiabatic index.
"""
const GAMMA = Ref{Float64}(5.0 / 3.0)

"""
Mass of the tracers in internal code units. Its value comes from `All.TargetGasMass = All.TargetGasMassFactor * All.ReferenceGasPartMass` in the code.

It is only printed in the output files (`stdout_n`), as `All.TargetGasMass=XXXXX`
"""
const TRACER_MASS = Ref{Float64}(3.65456e-06)

"""
Solar metallicity, as used in Arepo.

# References

M. Asplund et al. (2006). *The new solar abundances - Part I: the observations*. Communications in Asteroseismology, **147**. [doi:10.1553/cia147s76](https://doi.org/10.1553/cia147s76)
"""
const SOLAR_METALLICITY = Ref{Float64}(0.0127)

@doc raw"""
Hubble constant in $\mathrm{Gyr^{-1}}$.

This value corresponds to $H_0 = 100 \, \mathrm{km} \, \mathrm{s}^{-1} \, \mathrm{Mpc}^{-1} = 0.102201 \, \mathrm{Gyr}^{-1} = $.
"""
const HUBBLE_CONSTANT = Ref{Float64}(0.102201)

#################
# Internal units
#################

"""
If physical lengths units will be used, instead of comoving lengths units.
"""
const PHYSICAL_UNITS = Ref{Bool}(false)

"""
Internal unit of length

The default value is the internal unit of length used in [IllustrisTNG](https://www.tng-project.org/data/docs/specifications/), equal to ``1.0  \\, \\mathrm{kpc}``.
"""
const INTERNAL_L_UNIT = Ref(3.085678e21u"cm")

"""
Internal unit of mass

The default value is the internal unit of mass used in [IllustrisTNG](https://www.tng-project.org/data/docs/specifications/), equal to ``10^{10} \\, \\mathrm{M_\\odot}``.
"""
const INTERNAL_M_UNIT = Ref(1.989e43u"g")

"""
Internal unit of velocity

The default value is the internal unit of velocity used in [IllustrisTNG](https://www.tng-project.org/data/docs/specifications/), equal to ``1.0 \\, \\mathrm{km \\, s^{-1}}``.
"""
const INTERNAL_V_UNIT = Ref(1.0e5u"cm * s^-1")

####################################
# Configuration auxiliary functions
####################################

"""
   parseUnitful(s::String)::Unitful.Quantity

Parse a string of the form `"40.0 kpc"` into a Unitful quantity.

!!! note

    Works for any unit string that Unitful.uparse accepts (e.g. "kpc", "km/s", "K"). The string must have a space between the value and the unit, otherwise it will not be parsed correctly (e.g. "40.0kpc" will not work, but "40.0 kpc" will). And the unit should have no spaces (e.g. "km/s" is correct, but "km / s" is not). This is because the function splits the string into two parts using whitespace as a delimiter, and then parses the first part as a number and the second part as a unit.

# Arguments

  - `s::String`: The string to parse, expected to be of the form `"<value> <unit>"`, e.g. `"40.0 kpc"`.

# Returns

  - The parsed quantity, e.g. `40.0u"kpc"`.
"""
function parseUnitful(s::String)::Unitful.Quantity

    # Split the string into the value and unit parts, stripping whitespaces
    parts = split(strip(s); limit=2)

    (
        length(parts) == 2 ||
        throw(ArgumentError("parseUnitful: Expected \"<value> <unit>\", got $(s)"))
    )

    val  = parse(Float64, parts[1])
    unit = Unitful.uparse(parts[2])

    return val * unit

end

"""
    parseAbundanceShift(d::Dict{String,Any})::Dict{Symbol,Float64}

Parse the abundance shift dictionary from the config file, merging it with the default values.

If a key in `d` is not in the defaults, it will be ignored.

# Arguments

  - `d::Dict{String,Any}`: The dictionary from the config file.

# Returns

    - The result of merging `d` with the default values of [`ABUNDANCE_SHIFT`](@ref).
"""
function parseAbundanceShift(d::Dict{String,Any})::Dict{Symbol,Float64}

    merged = copy(ABUNDANCE_SHIFT[])

    for (k, v) in d
        sym = Symbol(k)
        if haskey(merged, sym)
            merged[sym] = Float64(v)
        else
            (
                LOGGING[] &&
                @warn("parseAbundanceShift: Unknown element `$k` in the config file. \
                I will ignored it")
            )
        end
    end

    return merged

end

"""
List of config keys, their corresponding global variables, and parsing functions.
"""
const CONFIG_SCHEMA = Dict{String,Tuple{Any,Function}}(
    "mmap.memory_fraction"    => (MMAP_MEMORY_FRACTION, x -> Float64(x)),
    "galaxy.disk_r"           => (DISK_R,               parseUnitful),
    "galaxy.disk_height"      => (DISK_HEIGHT,          parseUnitful),
    "galaxy.box_l"            => (BOX_L,                parseUnitful),
    "galaxy.rotation_r"       => (ROTATION_R,           parseUnitful),
    "galaxy.bulge_r"          => (BULGE_R,           parseUnitful),
    "galaxy.age_resolution"   => (AGE_RESOLUTION,       parseUnitful),
    "abundances"              => (ABUNDANCE_SHIFT,      parseAbundanceShift),
    "arepo.snap_basename"     => (SNAP_BASENAME,        string),
    "arepo.gc_basename"       => (GC_BASENAME,          string),
    "arepo.hydrogen_massfrac" => (HYDROGEN_MASSFRAC,    x -> Float64(x)),
    "arepo.gamma"             => (GAMMA,                x -> Float64(x)),
    "arepo.tracer_mass"       => (TRACER_MASS,          x -> Float64(x)),
    "arepo.solar_metallicity" => (SOLAR_METALLICITY,    x -> Float64(x)),
    "arepo.hubble_constant"   => (HUBBLE_CONSTANT,      x -> Float64(x)),
    "units.physical_units"    => (PHYSICAL_UNITS,       x -> x::Bool),
    "units.internal_l_unit"   => (INTERNAL_L_UNIT,      parseUnitful),
    "units.internal_m_unit"   => (INTERNAL_M_UNIT,      parseUnitful),
    "units.internal_v_unit"   => (INTERNAL_V_UNIT,      parseUnitful),
)
