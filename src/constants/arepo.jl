####################################################################################################
# Constants and data structures for Arepo simulations
####################################################################################################

#################
# Code constants
#################

"""
Base name of the snapshot files, set in the code variable `SnapshotFileBase`.
"""
const SNAP_BASENAME = "snap"

"""
Base name of the group catalog files.
"""
const GC_BASENAME = "fof_subhalo_tab"

"""
Mass fraction of hydrogen.
"""
const HYDROGEN_MASSFRAC = 0.76

"""
Mass of the tracers in internal code units. Its value comes from `All.TargetGasMass = All.TargetGasMassFactor * All.ReferenceGasPartMass` in the code.
"""
const TRACER_MASS = 3.65456e-06

"""
Solar metallicity, as used in Arepo.

# References

M. Asplund et al. (2006). *The new solar abundances - Part I: the observations*. Communications in Asteroseismology, **147**. [doi:10.1553/cia147s76](https://doi.org/10.1553/cia147s76)
"""
const SOLAR_METALLICITY = 0.0127

"""
Constant for the initial condition of dust and metals.
"""
const C_xd = 0.2835953313674557

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

@doc raw"""
Hubble constant in $\mathrm{Gyr^{-1}}$.

This value corresponds to $H_0 = 0.102201 \, \mathrm{Gyr}^{-1} = 100 \, \mathrm{km} \, \mathrm{s}^{-1} \, \mathrm{Mpc}^{-1}$.
"""
const HUBBLE_CONSTANT = 0.102201

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

################
# Default units
################

"""
Default internal unit of length.
"""
const DEFAULT_L_UNIT = ILLUSTRIS_L_UNIT

"""
Default internal unit of mass.
"""
const DEFAULT_M_UNIT = ILLUSTRIS_M_UNIT

"""
Default internal unit of velocity.
"""
const DEFAULT_V_UNIT = ILLUSTRIS_V_UNIT

######################
# Cell/particle types
######################

"""
Code index for each type of cell/particle.

# References

See for example Gadget2 [User's Guide](https://wwwmpa.mpa-garching.mpg.de/gadget/users-guide.pdf), or Gadget4 [documentation](https://wwwmpa.mpa-garching.mpg.de/gadget4/).
"""
const LONG_PARTICLE_INDEX = Dict(
    :gas        => 0,
    :halo       => 1,
    :disk       => 2,
    :bulge      => 3,
    :stars      => 4,
    :black_hole => 5,
    :tracer     => 6,
)

"""
Human readable name for each type of cell/particle.
"""
const LONG_PARTICLE_NAMES = Dict(
    :gas        => "Gas cells",
    :halo       => "HR DM particles",
    :disk       => "IR DM particles",
    :bulge      => "LR DM particles",
    :stars      => "Stellar particles",
    :black_hole => "Black hole particles",
    :tracer     => "Tracer particles",
)

"""
Human readable name for each type of cell/particle.
"""
const ISOLATED_PARTICLE_NAMES = Dict(
    :gas   => "Gas cells",
    :halo  => "DM particles",
    :disk  => "Stellar disk",
    :bulge => "Stellar bulge",
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

#############################################################
# Quantities that can be in a snapshot or group catalog file
#############################################################

"""
Dimensional properties for the quantities in the code.
"""
const QUANTITIES = Dict(

    ######################
    # Snapshot quantities
    ######################

    "CLKT" => Qty("", Unitful.𝐓, :internal),
    "GAGE" => Qty("GFM_StellarFormationTime", Unitful.𝐓, :internal),
    "GME2" => Qty("GFM_Metals", Unitful.NoDims, Unitful.NoUnits),
    "GMET" => Qty("GFM_Metals", Unitful.NoDims, Unitful.NoUnits),
    "GZ  " => Qty("GFM_Metallicity", Unitful.NoDims, Unitful.NoUnits),
    "GZ2 " => Qty("GFM_Metallicity", Unitful.NoDims, Unitful.NoUnits),
    "ID  " => Qty("ParticleIDs", Unitful.NoDims, Unitful.NoUnits),
    "PAID" => Qty("ParentID", Unitful.NoDims, Unitful.NoUnits),
    "TRID" => Qty("TracerID", Unitful.NoDims, Unitful.NoUnits),
    "MASS" => Qty("Masses", Unitful.𝐌, :internal),
    "NE  " => Qty("ElectronAbundance", Unitful.NoDims, Unitful.NoUnits),
    "NH  " => Qty("NeutralHydrogenAbundance", Unitful.NoDims, Unitful.NoUnits),
    "NHP " => Qty("IonizedHydrogenAbundance", Unitful.NoDims, Unitful.NoUnits),
    "POS " => Qty("Coordinates", Unitful.𝐋, :internal),
    "PRES" => Qty("Pressure", Unitful.𝐌 * Unitful.𝐋^-1 * Unitful.𝐓^-2, :internal),
    "RHO " => Qty("Density", Unitful.𝐌 * Unitful.𝐋^-3, :internal),
    "SFR " => Qty("StarFormationRate", Unitful.𝐌 * Unitful.𝐓^-1, u"Msun * yr^-1"),
    "SOFT" => Qty("Softenings", Unitful.𝐋, :internal),
    "TEMP" => Qty("Temperature", Unitful.𝚯, u"K"),
    "U   " => Qty("InternalEnergy", Unitful.𝐋^2 * Unitful.𝐓^-2, :internal),
    "VEL " => Qty("Velocities", Unitful.𝐋 * Unitful.𝐓^-1, :internal),
    "TSTP" => Qty("TimeStep", Unitful.𝐓, :internal),
    "POT " => Qty("Potential", Unitful.𝐋^2 * Unitful.𝐓^-2, :pot),

    #####################
    # sfr.txt quantities
    #####################

    # Time or scale factor
    "SFC1" => Qty("", Unitful.𝐓, :internal),
    # Total stellar mass to be formed prior to stochastic sampling
    "SFC2" => Qty("", Unitful.𝐌, :internal),
    # Instantaneous star formation rate of all cells
    "SFC3" => Qty("", Unitful.𝐌 * Unitful.𝐓^-1, u"Msun * yr^-1"),
    # Instantaneous star formation rate of active cells
    "SFC4" => Qty("", Unitful.𝐌 * Unitful.𝐓^-1, u"Msun * yr^-1"),
    # Total mass in stars formed after stochastic sampling
    "SFC5" => Qty("", Unitful.𝐌, :internal),
    # Cumulative stellar mass formed
    "SFC6" => Qty("", Unitful.𝐌, :internal),

    ####################
    # EL_SFR quantities
    ####################

    # Integration time, for gas cells and stellar particles
    "ODIT" => Qty("ODE_IntegrationTime", Unitful.𝐓, u"Myr"),
    # Scale factor, for gas cells and stellar particles
    "PARA" => Qty("ODE_ParameterA", Unitful.NoDims, Unitful.NoUnits),
    # UVB photoionization rate, for gas cells and stellar particles
    "PARU" => Qty("ODE_ParameterUVB", Unitful.𝐓^-1, u"Myr^-1"),
    # LWB photodissociation rate, for gas cells and stellar particles
    "PARL" => Qty("ODE_ParameterLWB", Unitful.𝐓^-1, u"Myr^-1"),
    # Star formation time parameter, for gas cells and stellar particles
    "TAUS" => Qty("ODE_TauS", Unitful.𝐓, u"Myr"),
    # Gas density, for gas cells and stellar particles
    "RHOC" => Qty("ODE_ParameterCellDensity", Unitful.𝐋^-3, u"cm^-3"),
    # Gas metallicity, for gas cells and stellar particles
    "PARZ" => Qty("ODE_ParameterMetallicity", Unitful.NoDims, Unitful.NoUnits),
    # Column height, for gas cells and stellar particles
    "PARH" => Qty("ODE_ParameterColumnHeight", Unitful.𝐋, u"cm"),
    # Photodissociation parameter, for gas cells and stellar particles
    "ETAD" => Qty("ODE_ParameterEtaD", Unitful.NoDims, Unitful.NoUnits),
    # Photoionization parameter, for gas cells and stellar particles
    "ETAI" => Qty("ODE_ParameterEtaI", Unitful.NoDims, Unitful.NoUnits),
    # Mass recycling parameter, for gas cells and stellar particles
    "PARR" => Qty("ODE_ParameterR", Unitful.NoDims, Unitful.NoUnits),
    # Metallicity of the supernova ejecta, for gas cells and stellar particles
    "PAZN" => Qty("ODE_ParameterZsn", Unitful.NoDims, Unitful.NoUnits),
    # Gas fractions, for gas cells and stellar particles
    "FRAC" => Qty("ODE_Fractions", Unitful.NoDims, Unitful.NoUnits),
    # Star formation flag, for gas cells
    "SFFL" => Qty("ODE_SfFlag", Unitful.NoDims, Unitful.NoUnits),
    # Cold gas fraction, for gas cells and stellar particles
    "COLF" => Qty("ODE_ColdMassFrac", Unitful.NoDims, Unitful.NoUnits),
    # Parent gas mass (at the moment of star formation), for stellar particles
    "GMAS" => Qty("ODE_GasMass", Unitful.𝐌, :internal),
    # Parent SFR (at the moment of star formation), for stellar particles
    "GSFR" => Qty("ODE_GasSFR", Unitful.𝐌 * Unitful.𝐓^-1, u"Msun * yr^-1"),
    # Parent gas pressure (at the moment of star formation), for stellar particles
    "GPRE" => Qty("ODE_GasPressure", Unitful.𝐌 * Unitful.𝐋^-1 * Unitful.𝐓^-2, :internal),
    # Parent position (at the moment of star formation), for stellar particles
    "GPOS" => Qty("ODE_GasPosition", Unitful.𝐋, :internal),
    # Parent velocity (at the moment of star formation), for stellar particles
    "GVEL" => Qty("ODE_GasVelocity", Unitful.𝐋 * Unitful.𝐓^-1, :internal),

    ##############################
    # Halo (FoF group) quantities
    ##############################

    "G_CM"        => Qty("GroupCM", Unitful.𝐋, :internal),
    "G_LenType"   => Qty("GroupLenType", Unitful.NoDims, Unitful.NoUnits),
    "G_Mass"      => Qty("GroupMass", Unitful.𝐌, :internal),
    "G_MassType"  => Qty("GroupMassType", Unitful.𝐌, :internal),
    "G_Nsubs"     => Qty("GroupNsubs", Unitful.NoDims, Unitful.NoUnits),
    "G_Pos"       => Qty("GroupPos", Unitful.𝐋, :internal),
    "G_M_Crit200" => Qty("Group_M_Crit200", Unitful.𝐌, :internal),
    "G_R_Crit200" => Qty("Group_R_Crit200", Unitful.𝐋, :internal),
    "G_Vel"       => Qty("GroupVel", Unitful.𝐋 * Unitful.𝐓^-1, :gvel),

    ###############################
    # Subhalo (subfind) quantities
    ###############################

    "S_CM"          => Qty("SubhaloCM", Unitful.𝐋, :internal),
    "S_HalfmassRad" => Qty("SubhaloHalfmassRad", Unitful.𝐋, :internal),
    "S_LenType"     => Qty("SubhaloLenType", Unitful.NoDims, Unitful.NoUnits),
    "S_Mass"        => Qty("SubhaloMass", Unitful.𝐌, :internal),
    "S_MassType"    => Qty("SubhaloMassType", Unitful.𝐌, :internal),
    "S_Pos"         => Qty("SubhaloPos", Unitful.𝐋, :internal),
    "S_Vel"         => Qty("SubhaloVel", Unitful.𝐋 * Unitful.𝐓^-1, u"km * s^-1"),
)

####################
# EL_SFR quantities
####################

"""
List of symbols for every EL_SFR quantity associated with the stellar particles.
"""
const EL_SFR_STELLAR_QUANTITIES = [
    Symbol(:ode_stellar_, quantity) for
    quantity in [
        :integration_time,
        :parameter_a,
        :parameter_uvb,
        :parameter_lwb,
        :tau_s,
        :parameter_cell_density,
        :parameter_metallicity,
        :parameter_column_height,
        :parameter_eta_d,
        :parameter_eta_i,
        :parameter_r,
        :parameter_zsn,
        :fractions,
        :cold_mass_frac,
        :gas_mass,
        :gas_sfr,
        :gas_pressure,
        :gas_position,
        :gas_velocity,
    ]
]

"""
List of symbols for every EL_SFR quantity associated with the gas cells.
"""
const EL_SFR_GAS_QUANTITIES = [
    Symbol(:ode_gas_, quantity) for
    quantity in [
        :integration_time,
        :parameter_a,
        :parameter_uvb,
        :parameter_lwb,
        :tau_s,
        :parameter_cell_density,
        :parameter_metallicity,
        :parameter_column_height,
        :parameter_eta_d,
        :parameter_eta_i,
        :parameter_r,
        :parameter_zsn,
        :fractions,
        :sf_flag,
        :cold_mass_frac,
    ]
]

"""
List of symbols for every EL_SFR quantity.
"""
const EL_SFR_QUANTITIES = vcat(EL_SFR_STELLAR_QUANTITIES, EL_SFR_GAS_QUANTITIES)

"""
Dictionary mapping each EL_SFR quantity symbol with its [`QUANTITIES`](@ref) key.
"""
const EL_SFR_KEYS = Dict(
    :integration_time        => "ODIT",
    :parameter_a             => "PARA",
    :parameter_uvb           => "PARU",
    :parameter_lwb           => "PARL",
    :tau_s                   => "TAUS",
    :parameter_cell_density  => "RHOC",
    :parameter_metallicity   => "PARZ",
    :parameter_column_height => "PARH",
    :parameter_eta_d         => "ETAD",
    :parameter_eta_i         => "ETAI",
    :parameter_r             => "PARR",
    :parameter_zsn           => "PAZN",
    :fractions               => "FRAC",
    :sf_flag                 => "SFFL",
    :cold_mass_frac          => "COLF",
    :gas_mass                => "GMAS",
    :gas_sfr                 => "GSFR",
    :gas_pressure            => "GPRE",
    :gas_position            => "GPOS",
    :gas_velocity            => "GVEL",
)
