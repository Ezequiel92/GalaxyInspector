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
    "MACH" => Qty("Machnumber", Unitful.NoDims, Unitful.NoUnits),
    "GDCM" => Qty("GFM_DustCapMass", Unitful.𝐌, :internal),
    "TAGR" => Qty("GFM_DustTauGrowth", Unitful.𝐓, u"Gyr"),
    "GDZ " => Qty("GFM_DustMetallicity", Unitful.NoDims, Unitful.NoUnits),
    "GDZ2" => Qty("GFM_DustMetallicity", Unitful.NoDims, Unitful.NoUnits),
    "GDAG" => Qty("GFM_DustAGB", Unitful.NoDims, Unitful.NoUnits),
    "GDA2" => Qty("GFM_DustAGB", Unitful.NoDims, Unitful.NoUnits),
    "GDII" => Qty("GFM_DustSNII", Unitful.NoDims, Unitful.NoUnits),
    "GDI2" => Qty("GFM_DustSNII", Unitful.NoDims, Unitful.NoUnits),
    "GDIa" => Qty("GFM_DustSNIa", Unitful.NoDims, Unitful.NoUnits),
    "GDJ2" => Qty("GFM_DustSNIa", Unitful.NoDims, Unitful.NoUnits),

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
    # SNII fraction, for gas cells and stellar particles
    "PARS" => Qty("ODE_ParameterSNII", Unitful.𝐌^-1, u"Msun^-1"),
    # Redshift parameter, for gas cells and stellar particles
    "PARz" => Qty("ODE_Parameterz", Unitful.NoDims, Unitful.NoUnits),
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

    ######################
    # YA_SFR quantities
    ######################

    # Isothermal sound speed [cm * s^(-1)]
    "YACC" => Qty("YA_cs_cell", Unitful.𝐋 * Unitful.𝐓^-1, u"cm * s^-1"),
    # Gas temperature [K]
    "YATC" => Qty("YA_T_cell", Unitful.𝚯, u"K"),
    # Length scale of the gas cell [cm]
    "YALC" => Qty("YA_l_cell", Unitful.𝐋, u"cm"),
    # Specific energy (turbulent) dissipation rate [cm^2 * s^(-3)]
    "YAEP" => Qty("YA_epsilon", Unitful.𝐋^2 * Unitful.𝐓^-3, u"cm^2 * s^-3"),
    # Mach number [dimensionless]
    "YAMS" => Qty("YA_Mach_cell", Unitful.NoDims, Unitful.NoUnits),
    # Magnetic field [g^(1/2) * cm^(-1/2) * s^(-1)]
    "YABC" => Qty("YA_B_cell", Unitful.𝐌^(1/2) * Unitful.𝐋^(-1/2) * Unitful.𝐓^-1, u"g^(1/2) * cm^(-1/2) * s^-1"),
    # Alfvén velocity [cm * s^(-1)]
    "YAVA" => Qty("YA_vA", Unitful.𝐋 * Unitful.𝐓^-1, u"cm * s^-1"),
    # Velocity dispersion [cm * s^(-1)]
    "YASV" => Qty("YA_sigma_v", Unitful.𝐋 * Unitful.𝐓^-1, u"cm * s^-1"),
    # Alfvén Mach [dimensionless]
    "YAMA" => Qty("YA_MA", Unitful.NoDims, Unitful.NoUnits),
    # Magnetic pressure [g * cm^(-1) * s^(-2)]
    "YAPM" => Qty("YA_P_mag", Unitful.𝐌 * Unitful.𝐋^-1 * Unitful.𝐓^-2, u"g * cm^-1 * s^-2"),
    # Thermal-to-magnetic pressure ratio [dimensionless]
    "YAB " => Qty("YA_beta", Unitful.NoDims, Unitful.NoUnits),
    # Standard deviation for the logarithmic density PDF [dimensionless]
    "YASS" => Qty("YA_sigma_s", Unitful.NoDims, Unitful.NoUnits),
    # Mean logarithmic density [dimensionless]
    "YAS0" => Qty("YA_s0", Unitful.NoDims, Unitful.NoUnits),
    # Effective sound speed [cm * s^(-1)]
    "YACE" => Qty("YA_cs_eff", Unitful.𝐋 * Unitful.𝐓^-1, u"cm * s^-1"),
    # Bonnor-Ebert mass for the mean density (s = 0) [M☉]
    "YAMB" => Qty("YA_MBE0", Unitful.𝐌, u"Msun"),
    # First root of Mex'(s) [dimensionless]
    "YASR" => Qty("YA_sr_m", Unitful.NoDims, Unitful.NoUnits),
    # Maximum excess mass [M☉]
    "YAMF" => Qty("YA_Mfe", Unitful.𝐌, u"Msun"),
    # Critical logarithmic density [dimensionless]
    "YASC" => Qty("YA_sc", Unitful.NoDims, Unitful.NoUnits),
    # Free fall time [yr]
    "YATF" => Qty("YA_tff", Unitful.𝐓, u"yr"),
    # status = 0 -> Success
    # status = 1 -> Failed to compute the Mach number
    # status = 2 -> We are in the sub-Alfvénic regime (MA < 2)
    # status = 3 -> There are no real roots. The cell is stable
    # status = 4 -> No roots exist for Mex(s). The cell is stable
    # status = 5 -> Could not find the root. The root is not within [s_left, sr_m]
    # status = 6 -> Could not find the root. Max iterations reached with no convergence
    "YAS " => Qty("YA_status", Unitful.NoDims, Unitful.NoUnits),

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

#####################
# Derived quantities
#####################

"""
List of symbols for the physical components of a simulation.
"""
const COMPONENTS = [
    ############################
    # Particle-based components
    ############################
    :stellar,     # Stellar particles
    :dark_matter, # Dark matter particles
    :black_hole,  # Black hole particles
    :Z_stellar,   # Metals in the stars
    #######################
    # Gas-based components
    #######################
    :gas,          # Total gas
    :hydrogen,     # Hydrogen
    :helium,       # Helium
    :Z_gas,        # Metals in the gas
    :ionized,      # Ionized gas (using the Arepo data)
    :neutral,      # Neutral gas (using the Arepo data)
    :br_atomic,    # Atomic gas (using the Blitz et al. (2006) relation)
    :br_molecular, # Molecular gas (using the Blitz et al. (2006) relation)
    ###############################
    # Components from our SF model
    ###############################
    :ode_ionized,           # Ionized gas
    :ode_atomic,            # Atomic gas
    :ode_molecular,         # Molecular gas
    :ode_stellar,           # Stars
    :ode_metals,            # Metals
    :ode_dust,              # Dust
    :ode_neutral,           # Neutral gas (everything but ionized, metals, and dust)
    :ode_molecular_stellar, # Molecular gas + stars
    :ode_cold,              # Cold gas (everything but atomic and ionized)
]

"""
List of symbols for physical magnitudes.
"""
const MAGNITUDES = [
    #######################
    # Mass-like magnitudes
    #######################
    :mass,            # Total mass [M]
    :mass_density,    # Volume mass density [M * L^-3]
    :number_density,  # Number density [L^-3]
    :area_density,    # Surface density [M * L^-2]
    :number,          # Number of elements [dimensionless]
    :fraction,        # Mass fraction [dimensionless]
    :eff,             # Star formation efficiency per free-fall time [dimensionless]
    :clumping_factor, # Clumping factor [dimensionless]
    #######################
    # Cinematic magnitudes
    #######################
    :specific_z_angular_momentum, # Angular momentum in the z direction per unit mass [L^2 * T^-1]
    :z_angular_momentum,          # Angular momentum in the z direction [M * L^2 * T^-1]
    :spin_parameter,              # Spin parameter [dimensionless]
    :circularity,                 # Circularity [dimensionless]
    :circular_velocity,           # Circular velocity [L * T^-1]
    :radial_velocity,             # Component of the velocity in the radial direction [L * T^-1]
    :tangential_velocity,         # Component of the velocity in the tangential direction [L * T^-1]
    :zstar_velocity,              # Component of the velocity in the z direction [L * T^-1]
    :kinetic_energy,              # Kinetic energy [M * L^2 * T^-2]
    :potential_energy,            # Potential energy [M * L^2 * T^-2]
    :total_energy,                # Total energy = kinetic + potential [M * L^2 * T^-2]
    ########
    # Other
    ########
    :depletion_time,  # Depletion time [T]
    :xy_distance,     # Projected (to the xy plane) distance to the origin [L]
    :radial_distance, # Distance to the origin [L]
]

"""
List of symbols for the derived quantities (component + magnitude).
"""
const DERIVED_QTY = [
    Symbol(component, :_, magnitude) for component in COMPONENTS for magnitude in MAGNITUDES
]

"""
Dictionary mapping each derived quantity to its component and magnitude.
"""
const QUANTITY_SPLITS = Dict(
    Symbol(component, :_, magnitude) => (magnitude, component)
    for component in COMPONENTS
    for magnitude in MAGNITUDES
)

####################################
# Code quantities from our SF model
####################################

"""
List of symbols for every code quantity associated with a gas cell in our SF model.
"""
const SFM_GAS_QTY = [
    Symbol(:ode_gas_, quantity) for
    quantity in [
        :integration_time,        # Integration time
        :parameter_z,             # Redshift
        :parameter_a,             # Scale factor
        :parameter_uvb,           # UVB photoionization rate
        :parameter_lwb,           # LWB photodissociation rate
        :tau_s,                   # Star formation time parameter
        :parameter_cell_density,  # Gas density
        :parameter_metallicity,   # Gas metallicity
        :parameter_column_height, # Column height
        :parameter_eta_d,         # Photodissociation parameter
        :parameter_eta_i,         # Photoionization parameter
        :parameter_r,             # Mass recycling parameter
        :parameter_zsn,           # Metallicity of the supernova ejecta
        :fractions,               # Gas fractions
        :sf_flag,                 # Star formation flag
        :cold_mass_frac,          # Cold gas fraction
    ]
]

"""
List of symbols for every code quantity associated with a stellar particle in our SF model.

This values come from the state of the parent gas cell (at the stellar birth time).
"""
const SFM_STELLAR_QTY = [
    Symbol(:ode_stellar_, quantity) for
    quantity in [
        :integration_time,        # Integration time
        :parameter_z,             # Redshift (birth time)
        :parameter_a,             # Scale factor (birth time)
        :parameter_uvb,           # UVB photoionization rate
        :parameter_lwb,           # LWB photodissociation rate
        :tau_s,                   # Star formation time parameter
        :parameter_cell_density,  # Gas density
        :parameter_metallicity,   # Gas metallicity
        :parameter_column_height, # Column height
        :parameter_eta_d,         # Photodissociation parameter
        :parameter_eta_i,         # Photoionization parameter
        :parameter_r,             # Mass recycling parameter
        :parameter_zsn,           # Metallicity of the supernova ejecta
        :fractions,               # Parent gas fractions
        :cold_mass_frac,          # Cold gas fraction
        :gas_mass,                # Parent gas cell mass
        :gas_sfr,                 # Parent gas cell SFR
        :gas_pressure,            # Parent gas cell pressure
    ]
]

"""
Dictionary mapping each code quantity in our SF model with its [`QUANTITIES`](@ref) key.
"""
const SFM_KEYS = Dict(
    :integration_time        => "ODIT",
    :parameter_z             => "PARz",
    :parameter_a             => "PARA",
    :parameter_snii          => "PARS",
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
)

"""
List of symbols for every code quantity in our SF model.
"""
const SFM_QTY = vcat(SFM_STELLAR_QTY, SFM_GAS_QTY)

"""
Dictionary mapping each code quantity in our SF model to its component and magnitude.
"""
const SFM_QTY_SPLITS = Dict(
    Symbol(component, magnitude) => (magnitude, component)
    for component in [:ode_gas_, :ode_stellar_]
    for magnitude in keys(SFM_KEYS)
)

############################################
# Derived code quantities from our SF model
############################################

"""
List of symbols for every derived code magnitude in our SF model.
"""
const SFM_DERIVED_MAGNITUDES = [
    :tau_star,          # Star formation timescale
    :tau_rec,           # Recombination timescale
    :tau_cond,          # Condensation timescale
    :tau_dg,            # Dust growth timescale
    :tau_dc,            # Dust creation timescale
    :tau_dd,            # Dust destruction timescale
    :tau_ion,           # Ionization optical depth
    :tau_diss,          # Dissociation optical depth
    :S_d,               # Dust shielding factor
    :S_H2,              # H2 self-shielding factor
    :equivalent_size,   # Diameter of a sphere with the same density and mass as the gas cell
]

"""
List of symbols for every derived code quantity associated with a gas cell in our SF model.
"""
const SFM_GAS_DERIVED_QTY = [Symbol(:ode_gas_, magnitude) for magnitude in SFM_DERIVED_MAGNITUDES]

"""
List of symbols for every derived code quantity associated with a stellar particle in our SF model.

This values come from the state of the parent gas cell (at the stellar birth time).
"""
const SFM_STELLAR_DERIVED_QTY = [
    Symbol(:ode_stellar_, magnitude)
    for magnitude in SFM_DERIVED_MAGNITUDES
]

"""
List of symbols for every derived code quantity in our SF model.
"""
const SFM_DERIVED_QTY = vcat(SFM_STELLAR_DERIVED_QTY, SFM_GAS_DERIVED_QTY)

"""
Dictionary mapping each derived code quantity to its component and magnitude.
"""
const SFM_DERIVED_QTY_SPLITS = Dict(
    Symbol(component, magnitude) => (magnitude, component)
    for component in [:ode_gas_, :ode_stellar_]
    for magnitude in SFM_DERIVED_MAGNITUDES
)

#######################
# Abundance quantities
#######################

"""
List of symbols for the quantities of metal mass in the gas.
"""
const GAS_METALS_MASS = Symbol.(keys(ELEMENT_INDEX), :_gas_mass)

"""
Dictionary mapping each metal mass quantity to its element.
"""
const GAS_METALS_MASS_SPLITS = Dict(
    Symbol(element, :_gas_mass) => element for element in keys(ELEMENT_INDEX)
)

"""
List of symbols for the quantities of metal mass in the stars.
"""
const STELLAR_METALS_MASS = Symbol.(keys(ELEMENT_INDEX), :_stellar_mass)

"""
Dictionary mapping each metal mass quantity to its element.
"""
const STELLAR_METALS_MASS_SPLITS = Dict(
    Symbol(element, :_stellar_mass) => element for element in keys(ELEMENT_INDEX)
)

"""
List of symbols for the gas abundance quantities.
"""
const GAS_ABUNDANCE = Symbol.(keys(ELEMENT_INDEX), :_gas_abundance)

"""
Dictionary mapping each gas abundance quantity to its element.
"""
const GAS_ABUNDANCE_SPLITS = Dict(
    Symbol(element, :_gas_abundance) => element for element in keys(ELEMENT_INDEX)
)

"""
List of symbols for the stellar abundance quantities.
"""
const STELLAR_ABUNDANCE = Symbol.(keys(ELEMENT_INDEX), :_stellar_abundance)

"""
Dictionary mapping each stellar abundance quantity to its element.
"""
const STELLAR_ABUNDANCE_SPLITS = Dict(
    Symbol(element, :_stellar_abundance) => element for element in keys(ELEMENT_INDEX)
)

"""
List of symbols for all the metals quantities.
"""
const METALS_QTY = vcat(GAS_METALS_MASS, STELLAR_METALS_MASS, GAS_ABUNDANCE, STELLAR_ABUNDANCE)

###########################
# Group catalog quantities
###########################

"""
    parseHaloQuantity(s::Symbol)::Tuple{Symbol,Int}

Takes a symbol for an halo (FoF group) quantity in the format `:quantity_haloidx` and splits it into
the text part as a Symbol and the number part as an Int.

# Arguments

  - `s::Symbol`: Target halo and halo quantity, e.g., `:halo_mass_12`. The index of the target halo starts at 1 and the halo quantity has to be one of the keys of [`HALO_KEYS`](@ref).

# Returns

  - A tuple with two elements:

      + A Symbol for the target quantity.
      + The index of the target halo.

# Examples

```julia-repl
julia> parseHaloQuantity(:halo_mass_12)
(:halo_mass, 12)

julia> parseHaloQuantity(:halo_M200_1)
(:halo_M200, 1)
```
"""
function parseHaloQuantity(s::Symbol)::Tuple{Symbol,Int}

    s_str = string(s)

    m = match(r"^halo_(.*)_(\d+)$", s_str)

    (
        isnothing(m) &&
        throw(ArgumentError("Input symbol must be in the format :quantity_haloidx. Got $s"))
    )

    # Extract captured groups
    prefix = m.captures[1]
    number = m.captures[2]

    return (Symbol(:halo_, Symbol(prefix)), parse(Int, number))

end

"""
Dictionary mapping each halo quantity with its [`QUANTITIES`](@ref) key.
"""
const HALO_KEYS = Dict(
    :halo_mass       => "G_Mass",
    :halo_n_subhalos => "G_Nsubs",
    :halo_M200       => "G_M_Crit200",
    :halo_R200       => "G_R_Crit200",
)

const HALO_QTY = keys(HALO_KEYS)

###################
# Extra quantities
###################

const EXTRA_QTY = [
    :temperature,
    :pressure,
    :sfr,
    :ssfr,
    :observational_sfr,
    :observational_ssfr,
    :sfr_area_density,
    :sfr_density,
    :stellar_age,
    :stellar_birth_time,
    :gas_sfr,
    :gas_sfr_area_density,
    :scale_factor,
    :redshift,
    :physical_time,
    :lookback_time,
    :time_step,
    :clock_time_s,
    :clock_time_percent,
    :tot_clock_time_s,
    :tot_clock_time_percent,
    :gas_metallicity,
    :stellar_metallicity,
    :ode_metallicity,
    :mass_accretion,
    :molecular_stellar_fraction,
]

###################################
# Global list of single quantities
###################################

"""
List of symbols for all the quantities in GalaxyInspector, except ratios.
"""
QTY_SINGLE_LIST = DERIVED_QTY ∪ SFM_QTY ∪ SFM_DERIVED_QTY ∪ HALO_QTY ∪ EXTRA_QTY ∪ METALS_QTY

#########
# Ratios
#########

"""
List of symbols for the ratios between derived quantities.
"""
const RATIO_QTY = [
    Symbol(qty_01, :_, qty_2, :_ratio) for qty_01 in QTY_SINGLE_LIST for qty_2 in QTY_SINGLE_LIST
]

"""
Dictionary mapping each ratio to its two derived quantities.
"""
const RATIO_SPLITS = Dict(
    Symbol(qty_01, :_, qty_2, :_ratio) => (qty_01, qty_2)
    for qty_01 in QTY_SINGLE_LIST
    for qty_2 in QTY_SINGLE_LIST
)

############################
# Global list of quantities
############################

"""
List of symbols for all the quantities in GalaxyInspector.
"""
const QTY_GLOBAL_LIST = QTY_SINGLE_LIST ∪ RATIO_QTY
