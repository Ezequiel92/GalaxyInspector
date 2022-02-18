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
    bins = 200,
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
        @test hist_data == file["hist_data"]
        @test corr_x == file["corr_x"]
        @test corr_y == file["corr_y"]
        @test prof_x == file["prof_x"]
        @test prof_y == file["prof_y"]
        @test zprof_x == file["zprof_x"]
        @test zprof_y == file["zprof_y"]
        @test cmdf_x == file["cmdf_x"]
        @test cmdf_y == file["cmdf_y"]
        @test ks_x == file["ks_x"]
        @test ks_y == file["ks_y"]
        # @test gd_x == file["gd_x"]
        # @test gd_y == file["gd_y"]
        # @test gd_z == file["gd_z"]
        @test pm_x == file["pm_x"]
        @test pm_y == file["pm_y"]
        @test pm_z == file["pm_z"]
        @test qe_x == file["qe_x"] 
        @test qe_y == file["qe_y"]
    end
end