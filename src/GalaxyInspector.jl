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
# Julia module for analyzing and plotting the output of GAGET2/3/4 and Arepo simulations             
####################################################################################################

module GalaxyInspector

using GadgetUnits, SPHtoGrid, SPHKernels
import GadgetIO as GIO
using Unitful, UnitfulAstro, CairoMakie, Colors, AverageShiftedHistograms, GLM, Glob, CSV, FileIO, 
DelimitedFiles, ProgressMeter, LinearAlgebra, DataFrames, QuadGK, LaTeXStrings, PrettyTables, HDF5, 
Tables

if isdefined(Base, :Experimental) && isdefined(Base.Experimental, Symbol("@optlevel"))
    @eval Base.Experimental.@optlevel 3
end

####################################################################################################
# Functions
####################################################################################################

include("constants.jl")
include("utilities.jl")
include("data_acquisition.jl")
include("data_analysis.jl")
include("post_processing.jl")
include("pipelines.jl")
 
export snapshotPlot, timeSeriesPlot, snapshotTable, timeSeriesTable

end
