####################################################################################################
# Reference run of the utilities in `src/utilities.jl`
####################################################################################################

gas_z = snap_data_04[:gas]["Z"]
gas_mass = snap_data_04[:gas]["MASS"]
gas_pos = snap_data_01[:gas]["POS"]
gas_u = snap_data_03[:gas]["U"]
gas_ne = snap_data_03[:gas]["NE"]
stellar_age = snap_data_02[:stars]["AGE"]
stellar_mass = snap_data_04[:stars]["MASS"]
total_star_mass = GI.computeTimeSeries(snapshots["snap_paths"], "mass", :stars, false)

pass_crit_rho = GI.passCritRho(
    snapshot,
    type = :gas,
    crit_ρ = GI.CRITICAL_DENSITY,
    sim_cosmo = SIM_COSMO,
)
pass_positive_qty = GI.passPositiveQty(snapshot)
pass_metallicity = GI.passMetallicity(snapshot, 1e-6, 1e-3)
metal_mass = GI.metalMass(gas_z)
metallicity = GI.computeMetallicity(gas_z, gas_mass, solar = true)
element_fraction = GI.computeElementFraction(gas_z, "O")
time = GI.computeTime([0:0.01:1...], cosmological_header)
redshift = GI.computeRedshift([0:0.01:1...])
stellar_age = GI.computeStellarAge(stellar_age, sim_data, snap_data)
sfr = GI.computeSFR(total_star_mass, time)
distance = GI.computeDistance(gas_pos)
center_of_mass = GI.computeCenterOfMass(gas_pos, gas_mass)
temperature = GI.computeTemperature(gas_z, gas_mass, gas_u, gas_ne)
surface_density_disc = GI.computeSurfaceDensity(distance, gas_mass, BOX_SIZE, 100)
surface_density_square = GI.computeSurfaceDensity(gas_pos, gas_mass, BOX_SIZE, 100)
basic_profile = GI.computeProfile(distance, gas_mass, BOX_SIZE, 100)
full_profile = GI.computeProfile(
    distance,
    gas_mass,
    BOX_SIZE,
    100,
    flat = false,
    cumulative = true,
    density = true,
)
basic_z_profile = GI.computeZProfile(distance, gas_z, gas_mass, BOX_SIZE, 100)
full_z_profile = GI.computeZProfile(
    distance,
    gas_z,
    gas_mass,
    BOX_SIZE,
    100,
    flat = false,
    cumulative = true,
    density = true,
    solar = true,
)
time_series_clock_time = GI.computeTimeSeries(snapshots["snap_paths"], "clock_time", :stars, false)
time_series_number = GI.computeTimeSeries(snapshots["snap_paths"], "number", :stars, false)
time_series_sfr = GI.computeTimeSeries(snapshots["snap_paths"], "sfr", :stars, false)

####################################################################################################
# Save the results as a reference for testing 
####################################################################################################

jldsave(
    joinpath(BASE_OUT_PATH, "utilities.jld2");
    pass_crit_rho,
    pass_positive_qty,
    pass_metallicity,
    total_star_mass,
    metal_mass,
    metallicity,
    element_fraction,
    time,
    redshift,
    stellar_age,
    sfr,
    distance,
    center_of_mass,
    temperature,
    surface_density_disc,
    surface_density_square,
    basic_profile,
    full_profile,
    basic_z_profile,
    full_z_profile,
    time_series_clock_time,
    time_series_number,
    time_series_sfr
)