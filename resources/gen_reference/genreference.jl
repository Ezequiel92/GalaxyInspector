####################################################################################################
# Generation of the reference results used in testing
####################################################################################################

push!(LOAD_PATH, joinpath(@__DIR__, "../../src/"))
using GalaxyInspector;
const GI = GalaxyInspector;
using JLD2, Unitful, UnitfulAstro, CairoMakie, Tar, Pkg.Artifacts, SHA
import GadgetIO as GIO

const NAME = "reference_results"

# Path to the directory where the reference results will be saved
const BASE_OUT_PATH = joinpath(@__DIR__, "../$NAME")
# Path to the directory containing the snapshot files used to generate the reference results
const BASE_SRC_PATH = artifact"example_data"
# Base name of the snapshot files, set in the GADGET variable `SnapshotFileBase`
const SNAP_NAME = "snap"
# Side dimension of the simulated region, for the case of vacuum boundary conditions
const BOX_SIZE = 200.0UnitfulAstro.kpc
# If the simulation is cosmological, 
#   - `false` ⟶ Newtonian simulation (`ComovingIntegrationOn` = 0).
#   - `true` ⟶ Cosmological simulation (`ComovingIntegrationOn` = 1).
const SIM_COSMO = false
# Index of one particular snapshot for testing
const SNAP_N = 16

mkdir(BASE_OUT_PATH)

include("reference_data_acquisition.jl")
include("reference_utilities.jl")
include("reference_data_analysis.jl")
include("reference_post_processing.jl")
include("reference_pipelines.jl")

# Compress results to be used as artifacts and delete files
tar_file = joinpath(@__DIR__, "../$NAME.tar.gz")
Tar.create(BASE_OUT_PATH, tar_file)
bind_artifact!(
    joinpath(@__DIR__, "../../Artifacts.toml"), 
    NAME, 
    Base.SHA1(Tar.tree_hash(tar_file)),
    download_info = [
        (
            "https://github.com/Ezequiel92/GalaxyInspector/raw/main/resources/$NAME.tar", 
            bytes2hex(open(sha256, tar_file)),
        ),
    ],
    force = true,
)

rm(BASE_OUT_PATH, force = true)
