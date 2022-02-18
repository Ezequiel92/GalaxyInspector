####################################################################################################
# Generation of the reference results used in testing
####################################################################################################

push!(LOAD_PATH, joinpath(@__DIR__, "../../../src/"))
using GadgetInspector; const GI = GadgetInspector;
using JLD2, Unitful, UnitfulAstro, GadgetIO, CairoMakie

# Path to the directory where the reference results will be saved
const BASE_OUT_PATH = joinpath(@__DIR__, "..")
# Path to the directory containing the snapshot files used to generate the reference results
const BASE_SRC_PATH = joinpath(@__DIR__, "example_source_data")
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

include("reference_data_acquisition.jl")
include("reference_utilities.jl")
include("reference_data_analysis.jl")
include("reference_post_processing.jl")
include("reference_pipelines.jl")
