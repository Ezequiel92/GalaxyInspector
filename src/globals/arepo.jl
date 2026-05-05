####################################################################################################
# Constants and data structures for Arepo simulations
####################################################################################################

"""
Subhalo numbers for the MW and M31 in Hestia simulations.
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

"""
Cosmological threshold density above which the gas cells/particles can turn into stars.

This value corresponds to `CritOverDensity` ``= 57.7 \\, [\\mathrm{cm^{-3}}]`` in the `param.txt` file (used only in cosmological simulations). Which is converted to internal units within the code using `OverDensThresh` = `CritOverDensity` * `OmegaBaryon` * 3 * `Hubble` * `Hubble` / (8 * `M_PI` * `G`)`. Then, to go to physical units again one has to do: `OverDensThresh`*`UnitDensity_in_cgs`*`cf_a3inv`*`HubbleParam`*`HubbleParam`.

Using the unit factors,

`UnitLength_in_cm`         = ``3.085678 \\times 10^{24}``

`UnitMass_in_g`            = ``1.989 \\times 10^{43}``

`UnitVelocity_in_cm_per_s` = ``100000``

The derived units,

`UnitTime_in_s`      = `UnitLength_in_cm` * `UnitVelocity_in_cm_per_s`^-1 = ``3.08568 \\times 10^{19}``

`UnitDensity_in_cgs` = `UnitMass_in_g` * `UnitLength_in_cm^-3`            = ``6.76991 \\times 10^{-31}``

The parameters,

`OmegaBaryon`       = ``0.048``

`HubbleParam`       = ``0.6777``

`PROTONMASS`        = ``1.67262178 \\times 10^{-24}``

`HYDROGEN_MASSFRAC` = ``0.76``

`GRAVITY`           = ``6.6738 \\times 10^{-8}``

`HUBBLE`            = ``3.2407789 \\times 10^{-18}``

And the derived parameters,

Hubble = `HUBBLE` * `UnitTime_in_s`                                              = ``100``

G      = `GRAVITY` * `UnitLength_in_cm`^-3 * `UnitMass_in_g` * `UnitTime_in_s`^2 = ``43.0187``

One gets,

`OverDensThresh` = 76.8495 [internal units of density]

And, for a cosmological simulation at redshift 0 (`cf_a3inv` = 1), this result in a physical density threshold of ``1.42857 \\times 10^{-5} \\, [\\mathrm{cm^{-3}}]``, or, adding the proton mass, a value of

``\\log_{10} \\rho \\ [\\mathrm{M_\\odot \\, kpc^{-3}}] = 2.548``
"""
const COSMO_THRESHOLD_DENSITY = 353.059u"Msun*kpc^-3"

"""
Threshold density above which the gas cells/particles can turn into stars.

This value corresponds to `CritPhysDensity` ``= 0.318 \\, [\\mathrm{cm^{-3}}]`` in the `param.txt` file (used in cosmological and non-cosmological simulations). Which is converted to internal units within the code using `PhysDensThresh` = `CritPhysDensity` * `PROTONMASS` / `HYDROGEN_MASSFRAC` / `UnitDensity_in_cgs`. Then, to go to physical units again one has to do: `PhysDensThresh` * `UnitDensity_in_cgs` * `cf_a3inv` * `HubbleParam` * `HubbleParam`.

`PhysDensThresh` = ``1.03378 \\times 10^{6}`` [internal units of density]

For a cosmological simulation at redshift 0 (`cf_a3inv` = 1), this result in a physical density threshold of ``0.192 \\, [\\mathrm{cm^{-3}}]``, or, adding the proton mass, a value of

``\\log_{10} \\rho \\, [\\mathrm{M_\\odot \\, kpc^{-3}}] = 6.677``
"""
const THRESHOLD_DENSITY = 4.749326e6u"Msun*kpc^-3"

######################
# Output files paths
######################

"""
Relative path, within the simulation directory, of `sfr.txt`.
"""
const SFR_REL_PATH = "output/sfr.txt"

"""
Relative path, within the simulation directory, of `cpu.txt`.
"""
const CPU_REL_PATH = "output/cpu.txt"

######################
# Cell/particle types
######################

"""
Code index for each type of cell/particle, in C zero-based numbering

# References

See for example Gadget2 [User's Guide](https://wwwmpa.mpa-garching.mpg.de/gadget/users-guide.pdf), or Gadget4 [documentation](https://wwwmpa.mpa-garching.mpg.de/gadget4/).
"""
const LONG_PARTICLE_INDEX = Dict(
    :gas         => 0,
    :dark_matter => 1,
    :disk        => 2,
    :bulge       => 3,
    :stellar     => 4,
    :black_hole  => 5,
    :tracer      => 6,
)

"""
Human readable name for each type of cell/particle.
"""
const LONG_PARTICLE_NAMES = Dict(
    :gas         => "Gas cells",
    :dark_matter => "HR DM particles",
    :disk        => "IR DM particles",
    :bulge       => "LR DM particles",
    :stellar     => "Stellar particles",
    :black_hole  => "Black hole particles",
    :tracer      => "Tracer particles",
)

"""
Human readable name for each type of cell/particle.
"""
const ISOLATED_PARTICLE_NAMES = Dict(
    :gas          => "Gas cells",
    :dark_matter  => "DM particles",
    :disk         => "Stellar disk",
    :bulge        => "Stellar bulge",
)

"""
Human readable name for each morphological component.
"""
const MORPHOLOGICAL_COMPONENTS = Dict(
    :disk       => "Disk",
    :bulge      => "Bulge",
    :thin_disk  => "Thin disk",
    :thick_disk => "Thick disk",
)

"""
Current cell/particle index in use.
"""
const PARTICLE_INDEX = LONG_PARTICLE_INDEX

"""
Current human readable name of each cell/particle type in use.
"""
const PARTICLE_NAMES = LONG_PARTICLE_NAMES

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

###################
# Tracked elements
###################

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
List of element indices above helium.
"""
const METAL_LIST = [3, 4, 5, 6, 7, 8, 9, 10]

#######################
# Star formation model
#######################

"""
ODE index corresponding to each gas phase in our star formation model.

  * Ionized gas fraction:   fi(t) = Mi(t) / MC --> data_dict[:gas]["FRAC"][1, :]
  * Atomic gas fraction:    fa(t) = Ma(t) / MC --> data_dict[:gas]["FRAC"][2, :]
  * Molecular gas fraction: fm(t) = Mm(t) / MC --> data_dict[:gas]["FRAC"][3, :]
  * Stellar fraction:       fs(t) = Ms(t) / MC --> data_dict[:gas]["FRAC"][4, :]
  * Metal fraction:         fZ(t) = MZ(t) / MC --> data_dict[:gas]["FRAC"][5, :]
  * Dust fraction:          fd(t) = Md(t) / MC --> data_dict[:gas]["FRAC"][6, :]
"""
const SFM_IDX = Dict(
    :ode_ionized   => 1,
    :ode_atomic    => 2,
    :ode_molecular => 3,
    :ode_stellar   => 4,
    :ode_metals    => 5,
    :ode_dust      => 6,
)

"""
Star formation model parameters
"""
const Cxd   = 0.2836                         # Constant for the initial condition of dust and metals.
const εff   = 1.0                            # Star formation efficiency per free-fall time
const αH    = 2.6e-13u"cm^3 * s^-1"          # Case B recombination coefficient at T = 10^4 K
const Rd    = 3.5e-17u"cm^3 * s^-1"          # H2 formation rate on dust grains at T = 100 K
const Cdg   = 2.4634u"mp^-1 * Myr^-1 * cm^3" # Dust growth coefficient
const σion  = 6.3e-18u"cm^2"                 # Hydrogen ionization cross section at the Lyman limit
const σdiss = 2.1e-19u"cm^2"                 # H2 dissociation cross section in the Lyman-Werner band
const σd    = 4.0e-21u"cm^2"                 # Effective dust cross-section per hydrogen atom
const ωH2   = 0.2                            # Constant for the H2 self-shielding factor
const xnf   = 5.0e14u"cm^-2"                 # H2 self-shielding column density factor
const Cswm  = 2.233e3u"Msun"                 # Mass swept-up by a SNe shock times the efficiency of dust destruction
