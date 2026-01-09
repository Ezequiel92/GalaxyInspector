####################################################################################################
#
#    ______        __                     ____                                 __                  _  __
#   / ____/____ _ / /____ _ _  __ __  __ /  _/____   _____ ____   ___   _____ / /_ ____   _____   (_)/ /
#  / / __ / __ `// // __ `/| |/_// / / / / / / __ \ / ___// __ \ / _ \ / ___// __// __ \ / ___/  / // /
# / /_/ // /_/ // // /_/ /_>  < / /_/ /_/ / / / / /(__  )/ /_/ //  __// /__ / /_ / /_/ // /_    / // /
# \____/ \__,_//_/ \__,_//_/|_| \__, //___//_/ /_//____// .___/ \___/ \___/ \__/ \____//_/(_)__/ //_/
#                              /____/                  /_/                                  /___/
####################################################################################################

####################################################################################################
# A Julia module for the analysis of galaxy simulations
####################################################################################################

module GalaxyInspector

using CSV,
    CairoMakie,
    ColorSchemes,
    Colors,
    DataFrames,
    Distances,
    FileIO,
    GLM,
    Glob,
    HDF5,
    Images,
    InvertedIndices,
    Interpolations,
    JLD2,
    LaTeXStrings,
    LinearAlgebra,
    Logging,
    Measurements,
    NearestNeighbors,
    ProgressMeter,
    QuadGK,
    Rotations,
    StaticArrays,
    Statistics,
    StatsBase,
    Unitful,
    UnitfulAstro,
    WriteVTK

####################################################################################################
# Optimization
####################################################################################################

@eval Base.Experimental.@optlevel 3

####################################################################################################
# Submodules
####################################################################################################

include("constants/globals.jl")

include("auxiliary_functions/grid.jl")
include("auxiliary_functions/histograms.jl")
include("auxiliary_functions/other.jl")
include("auxiliary_functions/plotting.jl")

include("analysis/data_acquisition.jl")
include("analysis/compute_quantities/energies.jl")
include("analysis/compute_quantities/masses.jl")
include("analysis/compute_quantities/positions.jl")
include("analysis/compute_quantities/sfm.jl")
include("analysis/compute_quantities/times.jl")
include("analysis/compute_quantities/velocities.jl")
include("analysis/compute_quantities/aggregators.jl")
include("analysis/filters.jl")
include("analysis/tracers.jl")
include("analysis/transformations.jl")
include("analysis/data_analysis.jl")

include("plotting/post_processing.jl")
include("plotting/pipelines.jl")
include("plotting/convenience.jl")

####################################################################################################
# Public functions
####################################################################################################

# From `analysis/data_acquisition.jl`
export readGroupCatalog
export readSnapshot
export getBlock
export makeDataDict

# From `plotting/pipelines.jl`
export plotSnapshot
export plotTimeSeries

# From `plotting/convenience.jl`
export sfrTXT
export cpuTXT
export stellarBirthHalos
export densityMap
export gasSFRMap
export densityMapVelField
export metallicityMap
export scatterPlot
export scatterDensityMap
export vtkFiles
export atomicMolecularTransition
export gasBarPlot
export timeSeries
export statisticsEvolution
export gasEvolution
export virialAccretionEvolution
export diskAccretionEvolution
export discAccretionEvolution
export rotationCurve
export radialProfile
export massProfile
export velocityProfile
export stellarHistory
export histogram
export compareFeldmann2020
export compareMolla2015
export compareAgertz2021
export kennicuttSchmidtLaw
export fitVSFLaw
export massMetallicityRelation
export gasVelocityCubes
export stellarVelocityCubes
export clumpingFactor
export circularityHistogram
export efficiencyHistogram
export stellarDensityMaps
export gasDensityMaps
export gasFractionsEvolution
export evolutionVideo
export snapshotReport
export simulationReport
export quantityReport

end
