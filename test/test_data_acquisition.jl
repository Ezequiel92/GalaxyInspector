####################################################################################################
# Tests for the data acquisition functions in `src/data_acquisition.jl`
####################################################################################################

snapshots = GI.getSnapshotPaths(SNAP_NAME, BASE_SRC_PATH)
snapshot = snapshots["snap_paths"][SNAP_N]
snap_data = GI.SnapData("", 12, 14, 15UnitfulAstro.Myr)
sim_data = GI.SimData(BASE_SRC_PATH, 1, :, SNAP_NAME, read_header(snapshot), nothing, 0.0)
isolated_header = read_header(snapshot)
cosmological_header = read_header(joinpath(
    BASE_SRC_PATH,
    "cosmological",
    SNAP_NAME * "dir_019",
    SNAP_NAME * "_019",
))

source_table = @test_nowarn GI.makeSourceTable([snapshots, snapshots], :, UnitfulAstro.Gyr)
temperature_data = @test_nowarn GI.getTemperature(snapshot)
raw_data = @test_nowarn GI.getRawData(snapshot, :gas, "POS")
snap_data_01 = @test_nowarn GI.getSnapshotData(snapshot, :gas, "POS")
snap_data_02 = @test_nowarn GI.getSnapshotData(snapshot, :stars, ["AGE", "MASS", "POS"])
snap_data_03 = @test_nowarn GI.getSnapshotData(snapshot, Dict(:gas => ["NE", "U"], :stars => ["Z"]))
snap_data_04 = @test_nowarn GI.getSnapshotData(snapshot, [:gas, :stars], ["Z", "MASS", "POS"])
sfr_file = @test_nowarn GI.getSfrFile(BASE_SRC_PATH, read_header(snapshot))
cpu_file_01 = @test_nowarn GI.getCpuFile(BASE_SRC_PATH, ["i/o", "hotngbs", "density"])
cpu_file_02 = @test_nowarn GI.getCpuFile(BASE_SRC_PATH, "cs_sfr")
molla_data = @test_nowarn GI.getMollá2015(BASE_SRC_PATH)

@testset "Data acquisition functions" begin
    jldopen(joinpath(BASE_DATA_PATH, "data_acquisition.jld2"), "r") do file
        @test GI.compare(snapshots["snap_numbers"], file["snapshots"]["snap_numbers"])
        @test GI.compare(temperature_data, file["temperature_data"])
        @test GI.compare(snap_data_01, file["snap_data_01"])
        @test GI.compare(snap_data_02, file["snap_data_02"])
        @test GI.compare(snap_data_03, file["snap_data_03"])
        @test GI.compare(snap_data_04, file["snap_data_04"])
        @test GI.compare(sfr_file, file["sfr_file"])
        @test GI.compare(cpu_file_01, file["cpu_file_01"])
        @test GI.compare(cpu_file_02, file["cpu_file_02"])
        @test GI.compare(molla_data, file["molla_data"])
    end
end
