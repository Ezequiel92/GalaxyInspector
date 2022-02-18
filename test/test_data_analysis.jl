####################################################################################################
# Tests for the data analysis functions in `src/data_analysis.jl`
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

pm_x, pm_y, pm_z = GI.particleMap(
    merge(Dict(:sim_data => sim_data, :snap_data => snap_data), snap_data_02),
    "XY",
    :stars,
    100.0UnitfulAstro.kpc,
)

qe_x, qe_y = GI.qtyEvolution(sim_data, "clock_time", "mass", nothing, :stars, warnings = false)

@testset "Data analysis functions" begin
    jldopen(joinpath(BASE_DATA_PATH, "data_analysis.jld2"), "r") do file
        @test GI.compare(hist_data, file["hist_data"])
        @test GI.compare(corr_x, file["corr_x"])
        @test GI.compare(corr_y, file["corr_y"])
        @test GI.compare(prof_x, file["prof_x"])
        @test GI.compare(prof_y, file["prof_y"])
        @test GI.compare(zprof_x, file["zprof_x"])
        @test GI.compare(zprof_y, file["zprof_y"])
        @test GI.compare(cmdf_x, file["cmdf_x"])
        @test GI.compare(cmdf_y, file["cmdf_y"])
        @test GI.compare(ks_x, file["ks_x"])
        @test GI.compare(ks_y, file["ks_y"])
        # @test GI.compare(gd_x, file["gd_x"])
        # @test GI.compare(gd_y, file["gd_y"])
        # @test GI.compare(gd_z, file["gd_z"])
        @test GI.compare(pm_x, file["pm_x"])
        @test GI.compare(pm_y, file["pm_y"])
        @test GI.compare(pm_z, file["pm_z"])
        @test GI.compare(qe_x, file["qe_x"])
        @test GI.compare(qe_y, file["qe_y"])
    end
end