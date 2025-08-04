####################################################################################################
# Constants and data structures
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

######################
# Cell/particle types
######################

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
Dictionary of dimensional properties for the quantities in the code.
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
    "SFR " => Qty("StarFormationRate", Unitful.𝐌 * Unitful.𝐓^-1, u"Msun*yr^-1"),
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
    "SFC3" => Qty("", Unitful.𝐌 * Unitful.𝐓^-1, u"Msun*yr^-1"),
    # Instantaneous star formation rate of active cells
    "SFC4" => Qty("", Unitful.𝐌 * Unitful.𝐓^-1, u"Msun*yr^-1"),
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
    "GSFR" => Qty("ODE_GasSFR", Unitful.𝐌 * Unitful.𝐓^-1, u"Msun*yr^-1"),
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
    "S_Vel"         => Qty("SubhaloVel", Unitful.𝐋 * Unitful.𝐓^-1, u"km*s^-1"),
)
