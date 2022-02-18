####################################################################################################
# Reference run of the data analysis functions in `src/data_analysis.jl`
####################################################################################################

hist_data, = GI.histogram(snap_data_01, :gas, "POS", func = GI.computeDistance)

corr_x, corr_y = GI.correlation(
    snap_data_04,
    (:gas, "POS"),
    (:gas, "MASS"),
    x_function = GI.computeDistance,
)

prof_x, prof_y = GI.profile(snap_data_04, (:gas, "MASS"), 200.0UnitfulAstro.kpc, 50)

zprof_x, zprof_y = GI.zProfile(snap_data_04, :stars, 200.0UnitfulAstro.kpc, 100)

cmdf_x, cmdf_y = GI.cmdf(
    merge(Dict(:sim_data => sim_data, :snap_data => snap_data), snap_data_04),
    bins = 100,
)

ks_data = GI.getSnapshotData(
    snapshot,
    Dict(:gas => ["POS", "MASS", "Z", "U", "NE"], :stars => ["POS", "MASS", "AGE"]),
)

ks_x, ks_y = GI.kennicuttSchmidt(
    merge(Dict(:sim_data => sim_data, :snap_data => snap_data), ks_data),
    Inf * Unitful.K,
    200.0UnitfulAstro.Myr,
    200.0UnitfulAstro.kpc;
    bins = 200
)

# gd_data = GI.getSnapshotData(
#     snapshot,
#     :gas,
#     ["MASS", "RHO", "HSML", "POS"],
# )
# gd_x, gd_y, gd_z = GI.gasDensityMap(
#     merge(Dict(:sim_data => sim_data,:snap_data => snap_data), gd_data),
#     "XY",
#     200.0UnitfulAstro.kpc,
# )

# pm_x, pm_y, pm_z = GI.particleMap(
#     merge(Dict(:sim_data => sim_data, :snap_data => snap_data), snap_data_02),
#     "XY",
#     :stars,
#     100.0UnitfulAstro.kpc,
# )

qe_x, qe_y = GI.qtyEvolution(sim_data, "clock_time", "mass", nothing, :stars, warnings = false)

####################################################################################################
# Save the results as a reference for testing 
####################################################################################################

jldsave(
    joinpath(BASE_OUT_PATH, "data_analysis.jld2");
    hist_data,
    corr_x, corr_y,
    prof_x, prof_y,
    zprof_x, zprof_y,
    cmdf_x, cmdf_y,
    ks_x, ks_y,
    # gd_x, gd_y, gd_z, 
    # pm_x, pm_y, pm_z,
    qe_x, qe_y
)
