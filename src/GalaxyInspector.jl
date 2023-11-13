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
    PGFPlotsX,
    ProgressMeter,
    QuadGK,
    Rotations,
    StatsBase,
    Unitful,
    UnitfulAstro

####################################################################################################
# Optimization.
####################################################################################################

@eval Base.Experimental.@optlevel 3

####################################################################################################
# Submodules.
####################################################################################################

include("constants.jl")
include("general_utilities.jl")
include("arepo_utilities.jl")
include("data_acquisition.jl")
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

# From `pipelines.jl`
export snapshotPlot
export timeSeriesPlot

# From `convenience.jl`
export snapshotReport
export simulationReport
export sfrTXT
export cpuTXT
export densityMap
export scatterPlot
export scatterDensityMap
export timeSeries
export rotationCurve
export densityProfile
export stellarHistory
export stellarCircularity
export compareWithFeldmann2020
export compareWithMolla2015
export compareWithKennicuttBigiel

end
