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
const GC_BASENAME = "groups"

"""
Mass fraction of hydrogen.
"""
const HYDROGEN_MASSFRAC = 0.76

"""
Solar metallicity, as used in OpenGadget3.
"""
const SOLAR_METALLICITY = 0.02

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
const PARTICLE_NAMES = ISOLATED_PARTICLE_NAMES

###################
# Tracked elements
###################

"""
Code index for each tracked element.
"""
const ELEMENT_INDEX = Dict(
    :He => 1,  # Helium
    :C  => 2,  # Carbon
    :Mg => 3,  # Magnesium
    :O  => 4,  # Oxygen
    :Fe => 5,  # Iron
    :Si => 6,  # Silicon
    :H  => 7,  # Hydrogen
    :N  => 8,  # Nitrogen
    :Ne => 9,  # Neon
    :S  => 10, # Sulfur
    :Ca => 11, # Calcium
)

"""
List of element indices above helium.
"""
const METAL_LIST = [2, 3, 4, 5, 6, 8, 9, 10, 11]

#############################################################
# Quantities that can be in a sanpshot or group catalog file
#############################################################

"""
Dictionary of dimensional properties for the quantities in the code.
"""
const QUANTITIES = Dict(
    ######################
    # Snapshot quantities
    ######################
    "CLKT" => Qty("", Unitful.𝐓, :internal),
    "GAGE" => Qty("StellarFormationTime", Unitful.𝐓, :internal),
    "GME2" => Qty("Metallicity", Unitful.𝐌, :internal),
    "GMET" => Qty("Metallicity", Unitful.𝐌, :internal),
    "ID  " => Qty("ParticleIDs", Unitful.NoDims, Unitful.NoUnits),
    "MASS" => Qty("Masses", Unitful.𝐌, :internal),
    "NE  " => Qty("ElectronAbundance", Unitful.NoDims, Unitful.NoUnits),
    "NH  " => Qty("NeutralHydrogenAbundance", Unitful.NoDims, Unitful.NoUnits),
    "NHP " => Qty("HII", Unitful.NoDims, Unitful.NoUnits),
    "POS " => Qty("Coordinates", Unitful.𝐋, :internal),
    "PRES" => Qty("Pressure", Unitful.𝐌 * Unitful.𝐋^-1 * Unitful.𝐓^-2, :internal),
    "RHO " => Qty("Density", Unitful.𝐌 * Unitful.𝐋^-3, :internal),
    "HSML" => Qty("SmoothingLength", Unitful.𝐋, :internal),
    "HSMS" => Qty("StellarSmoothingLength", Unitful.𝐋, :internal),
    "SFC1" => Qty("", Unitful.𝐓, :internal),
    "SFC2" => Qty("", Unitful.𝐌, :internal),
    "SFC3" => Qty("", Unitful.𝐌 * Unitful.𝐓^-1, u"Msun*yr^-1"),
    "SFC4" => Qty("", Unitful.𝐌 * Unitful.𝐓^-1, u"Msun*yr^-1"),
    "SFC5" => Qty("", Unitful.𝐌, :internal),
    "SFC6" => Qty("", Unitful.𝐌, :internal),
    "SFR " => Qty("StarFormationRate", Unitful.𝐌 * Unitful.𝐓^-1, u"Msun*yr^-1"),
    "TEMP" => Qty("Temperature", Unitful.𝚯, u"K"),
    "TSTP" => Qty("TimeStep", Unitful.𝐓, :internal),
    "U   " => Qty("InternalEnergy", Unitful.𝐋^2 * Unitful.𝐓^-2, :internal),
    "VEL " => Qty("Velocities", Unitful.𝐋 * Unitful.𝐓^-1, :internal),
    "Z   " => Qty("Metallicity", Unitful.NoDims, Unitful.NoUnits),
    "POS " => Qty("Coordinates", Unitful.𝐋, :internal),
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
    ##############################
    # Halo (FoF group) quantities
    ##############################
    "G_CM"          => Qty("GroupPos", Unitful.𝐋, :internal),
    "G_FirstSub"    => Qty("GroupFirstSub", Unitful.NoDims, Unitful.NoUnits),
    "G_Len"         => Qty("GroupLen", Unitful.NoDims, Unitful.NoUnits),
    #TODO
    "G_LenType"     => Qty("GroupLenType", Unitful.NoDims, Unitful.NoUnits),
    "G_Mass"        => Qty("GroupMass", Unitful.𝐌, :internal),
    #TODO
    "G_MassType"    => Qty("GroupMassType", Unitful.𝐌, :internal),
    "G_Nsubs"       => Qty("GroupNsubs", Unitful.NoDims, Unitful.NoUnits),
    "G_Pos"         => Qty("GroupPos", Unitful.𝐋, :internal),
    "G_M_Crit200"   => Qty("Group_M_Crit200", Unitful.𝐌, :internal),
    "G_M_Mean200"   => Qty("Group_M_Mean200", Unitful.𝐌, :internal),
    "G_M_TopHat200" => Qty("Group_M_TopHat200", Unitful.𝐌, :internal),
    "G_R_Crit200"   => Qty("Group_R_Crit200", Unitful.𝐋, :internal),
    "G_R_Mean200"   => Qty("Group_R_Mean200", Unitful.𝐋, :internal),
    "G_R_TopHat200" => Qty("Group_R_TopHat200", Unitful.𝐋, :internal),
    #TODO
    "G_Vel"         => Qty("GroupVel", Unitful.𝐋 * Unitful.𝐓^-1, :gvel),
    ###############################
    # Subhalo (subfind) quantities
    ###############################
    "S_CM"          => Qty("SubhaloCM", Unitful.𝐋, :internal),
    "S_GrNr"        => Qty("SubhaloGrNr", Unitful.NoDims, Unitful.NoUnits),
    "S_HalfmassRad" => Qty("SubhaloHalfmassRad", Unitful.𝐋, :internal),
    "S_IDMostbound" => Qty("SubhaloIDMostbound", Unitful.NoDims, Unitful.NoUnits),
    "S_Len"         => Qty("SubhaloLen", Unitful.NoDims, Unitful.NoUnits),
    # TODO
    "S_LenType"     => Qty("SubhaloLenType", Unitful.NoDims, Unitful.NoUnits),
    "S_Mass"        => Qty("SubhaloMass", Unitful.𝐌, :internal),
    # TODO
    "S_MassType"    => Qty("SubhaloMassType", Unitful.𝐌, :internal),
    "S_Parent"      => Qty("SubhaloParent", Unitful.NoDims, Unitful.NoUnits),
    "S_Pos"         => Qty("SubhaloPos", Unitful.𝐋, :internal),
    "S_Vel"         => Qty("SubhaloVel", Unitful.𝐋 * Unitful.𝐓^-1, u"km*s^-1"),
    "S_VelDisp"     => Qty("SubhaloVelDisp", Unitful.𝐋 * Unitful.𝐓^-1, u"km*s^-1"),
    "S_Vmax"        => Qty("SubhaloVmax", Unitful.𝐋 * Unitful.𝐓^-1, u"km*s^-1"),
    "S_VmaxRad"     => Qty("SubhaloVmaxRad", Unitful.𝐋, :internal),
)
