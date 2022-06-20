####################################################################################################
# Reference run of the data acquisition functions in `src/data_acquisition.jl`
####################################################################################################

snapshots = GI.getSnapshotPaths(SNAP_NAME, BASE_SRC_PATH)
snapshot = snapshots["snap_paths"][SNAP_N]
snap_data = GI.SnapData("", 12, 14, 15UnitfulAstro.Myr)
sim_data = GI.SimData(BASE_SRC_PATH, 1, :, SNAP_NAME, GI.read_header(snapshot), nothing, 0.0)
isolated_header = GI.read_header(snapshot)
cosmological_header = GI.read_header(joinpath(
    BASE_SRC_PATH,
    "cosmological",
    SNAP_NAME * "dir_019",
    SNAP_NAME * "_019",
))
velocities = GI.read_snapshots(
    snapshot,
    "VEL",
    GI.ParticleType[:gas],
    x -> GI.passAll(x, :gas),
)
temperature_data = GI.getTemperature(snapshot)
raw_data = GI.getRawData(snapshot, :gas, "POS")
snap_data_01 = GI.getSnapshotData(snapshot, :gas, "POS")
snap_data_02 = GI.getSnapshotData(snapshot, :stars, ["AGE", "MASS", "POS"])
snap_data_03 = GI.getSnapshotData(snapshot, Dict(:gas => ["NE", "U"], :stars => ["Z"]))
snap_data_04 = GI.getSnapshotData(snapshot, [:gas, :stars], ["Z", "MASS", "POS"])
sfr_file = GI.getSfrFile(BASE_SRC_PATH, GI.read_header(snapshot))
cpu_file_01 = GI.getCpuFile(BASE_SRC_PATH, ["i/o", "hotngbs", "density"])
cpu_file_02 = GI.getCpuFile(BASE_SRC_PATH, "cs_sfr")
molla_data = GI.getMollá2015(BASE_SRC_PATH)

####################################################################################################
# Save the results as a reference for testing 
####################################################################################################

jldsave(
    joinpath(BASE_OUT_PATH, "data_acquisition.jld2");
    snapshots,
    isolated_header,
    cosmological_header,
    velocities,
    temperature_data,
    snap_data_01,
    snap_data_02,
    snap_data_03,
    snap_data_04,
    sfr_file,
    cpu_file_01,
    cpu_file_02,
    molla_data
)
