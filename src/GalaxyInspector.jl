####################################################################################################
#
#    ______        __                     ____                                 __                  _  __
#   / ____/____ _ / /____ _ _  __ __  __ /  _/____   _____ ____   ___   _____ / /_ ____   _____   (_)/ /
#  / / __ / __ `// // __ `/| |/_// / / / / / / __ \ / ___// __ \ / _ \ / ___// __// __ \ / ___/  / // /
# / /_/ // /_/ // // /_/ /_>  < / /_/ /_/ / / / / /(__  )/ /_/ //  __// /__ / /_ / /_/ // /_    / // /
# \____/ \__,_//_/ \__,_//_/|_| \__, //___//_/ /_//____// .___/ \___/ \___/ \__/ \____//_/(_)__/ //_/
#                              /____/                  /_/                                  /___/
#
####################################################################################################
# A Julia module for the data analysis of Arepo hydrodynamical simulations.
####################################################################################################

module GalaxyInspector

using LinearAlgebra, Statistics
using CSV,
    CairoMakie,
    ColorSchemes,
    Colors,
    DataFrames,
    DelimitedFiles,
    FileIO,
    GLM,
    Glob,
    HDF5,
    Images,
    JLD2,
    LaTeXStrings,
    Measurements,
    NearestNeighbors,
    ProgressMeter,
    QuadGK,
    Rotations,
    StatsBase,
    Unitful,
    UnitfulAstro

####################################################################################################
# Choose codebase.
####################################################################################################

const CODEBASE = :arepo

####################################################################################################
# Optimization.
####################################################################################################

@eval Base.Experimental.@optlevel 3

####################################################################################################
# Submodules.
####################################################################################################

if CODEBASE == :arepo
    include("constants/arepo.jl")
elseif CODEBASE == :opengadget3
    include("constants/opengadget3.jl")
else
    throw(ArgumentError("GalaxyInspector: I don't recognize the codebase :$(CODEBASE)"))
end

include("general_utilities.jl")
include("arepo/data_acquisition.jl")
include("arepo/arepo_utilities.jl")
include("arepo/compute_quantities.jl")
include("arepo/filters.jl")
include("arepo/tracers.jl")
include("arepo/transformations.jl")
include("data_analysis.jl")
include("post_processing.jl")
include("pipelines.jl")
include("convenience.jl")

####################################################################################################
# Public functions.
####################################################################################################

# From `data_acquisition.jl`
export readGroupCatalog
export readSnapshot
export getBlock
export readSfrFile
export readCpuFile
export makeDataDict

# From `pipelines.jl`
export snapshotPlot
export timeSeriesPlot

# From `convenience.jl`
export snapshotReport
export simulationReport
export sfrTXT
export cpuTXT
export stellarBirthHalos
export densityMap
export metallicityMap
export temperatureMap
export densityMapVelField
export scatterPlot
export atomicMolecularTransitionHeatmap
export atomicMolecularTransitionScatter
export scatterDensityMap
export gasFractionsBarPlot
export timeSeries
export gasEvolution
export virialAccretionEvolution
export discAccretionEvolution
export rotationCurve
export densityProfile
export massProfile
export velocityProfile
export stellarHistory
export stellarCircularity
export compareFeldmann2020
export compareMolla2015
export compareKennicuttBigielResolved
export compareKennicuttBigielIntegrated
export fitKennicuttBigielResolved


end
