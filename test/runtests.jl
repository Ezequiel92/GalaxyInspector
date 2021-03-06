####################################################################################################
# GalaxyInspector testing
####################################################################################################

push!(LOAD_PATH, joinpath(@__DIR__, "../src/"))
using GalaxyInspector;
const GI = GalaxyInspector;
using JLD2, CairoMakie, LaTeXStrings, DataFrames, DelimitedFiles, Unitful, UnitfulAstro, Test, ReferenceTests, Pkg.Artifacts
import GadgetIO as GIO

# Path to the directory containing the reference results
const BASE_DATA_PATH = artifact"reference_results"
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

println("Testing GalaxyInspector.jl...\n")

include("test_data_acquisition.jl")
include("test_utilities.jl")
include("test_data_analysis.jl")
include("test_post_processing.jl")
include("test_pipelines.jl")
