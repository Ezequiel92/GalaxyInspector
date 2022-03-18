####################################################################################################
# 
#    ______            __              __   ____                                 __                  _  __
#   / ____/____ _ ____/ /____ _ ___   / /_ /  _/____   _____ ____   ___   _____ / /_ ____   _____   (_)/ /
#  / / __ / __ `// __  // __ `// _ \ / __/ / / / __ \ / ___// __ \ / _ \ / ___// __// __ \ / ___/  / // / 
# / /_/ // /_/ // /_/ // /_/ //  __// /_ _/ / / / / /(__  )/ /_/ //  __// /__ / /_ / /_/ // /_    / // /  
# \____/ \__,_/ \__,_/ \__, / \___/ \__//___//_/ /_//____// .___/ \___/ \___/ \__/ \____//_/(_)__/ //_/   
#                     /____/                             /_/                                  /___/       
# 
####################################################################################################
# Julia module for analyzing and plotting the output of GAGET2/3/4 simulations             
####################################################################################################

module GadgetInspector

using GadgetIO, GadgetUnits, SPHtoGrid, SPHKernels
using Unitful, UnitfulAstro, CairoMakie, Colors, AverageShiftedHistograms, GLM, Glob, CSV, FileIO, 
DelimitedFiles, ProgressMeter, LinearAlgebra, DataFrames, QuadGK, LaTeXStrings, PrettyTables

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
