####################################################################################################
# Reference run of the pipeline functions in `src/piepelines.jl`
####################################################################################################

# Histogram
snapshotPlot(
    [BASE_SRC_PATH,],
    [SNAP_NAME,],
    Dict(:gas => "NH"),
    hist!;
    output_path = joinpath(BASE_OUT_PATH, "snapshotPlot/histograms"),
    sim_cosmo = SIM_COSMO,
    idx = 16:17,
    sim_labels = nothing,
    t_unit = UnitfulAstro.Gyr,
    warnings = false,
    da_functions = [GI.histogram,],
    da_args = [(:gas, "NH",),],
    post_processing! = GI.verticalFlags!,
    pp_args = (([0.1, 0.5], ["0.1", "0.5"]),),
    x_label = "atomic fraction",
    y_label = "counts",
    y_scale = Makie.pseudolog10,
    animation = false,
    output_format = ".png"
)

# Scatter plot
snapshotPlot(
    [BASE_SRC_PATH, BASE_SRC_PATH,],
    [SNAP_NAME, SNAP_NAME,],
    Dict(:gas => "POS"),
    scatter!;
    output_path = joinpath(BASE_OUT_PATH, "snapshotPlot/scatter_plots"),
    sim_cosmo = SIM_COSMO,
    idx = 16:17,
    t_unit = UnitfulAstro.Gyr,
    warnings = false,
    da_functions = [GI.correlation,],
    da_args = [((:gas, "POS"), (:gas, "POS")),],
    da_kwargs = [(x_function = x -> x[1, :], y_function = y -> y[2, :]),],
    post_processing! = GI.cross!,
    x_name = "x axis",
    y_name = "y axis",
    x_unit = UnitfulAstro.kpc,
    y_unit = UnitfulAstro.kpc,
    x_limits = (-200, 200),
    y_limits = (-200, 200),
    animation = false,
    output_format = ".png",
    resolution = (1000, 1000),
    aspect = AxisAspect(1)
)

# Scatterlines plot
snapshotPlot(
    [BASE_SRC_PATH,],
    [SNAP_NAME,],
    Dict(:gas => ["MASS", "POS"]),
    scatterlines!;
    output_path = joinpath(BASE_OUT_PATH, "snapshotPlot/scatterlines"),
    sim_cosmo = SIM_COSMO,
    idx = 16:17,
    t_unit = UnitfulAstro.Gyr,
    warnings = false,
    da_functions = [GI.profile,],
    da_args = [((:gas, "MASS"), 300UnitfulAstro.kpc, 100),],
    da_kwargs = [(cumulative = true,),],
    x_name = "r",
    y_name = "Mass",
    x_unit = UnitfulAstro.kpc,
    y_unit = UnitfulAstro.Msun,
    y_factor = 10,
    animation = false,
    output_format = ".png"
)

# Heatmap
snapshotPlot(
    [BASE_SRC_PATH,],
    [SNAP_NAME,],
    Dict(:stars => "POS"),
    heatmap!;
    output_path = joinpath(BASE_OUT_PATH, "snapshotPlot/heatmaps"),
    sim_cosmo = SIM_COSMO,
    idx = 16:17,
    sim_labels = nothing,
    t_unit = UnitfulAstro.Gyr,
    warnings = false,
    da_functions = [GI.particleMap,],
    da_args = [("XY", :stars, BOX_SIZE * 0.25),],
    da_kwargs = [(sim_cosmo = SIM_COSMO,),],
    post_processing! = GI.cross!,
    pp_kwargs = (color = :white,),
    x_label = "12 + log(:auto_label)",
    x_name = "x axis",
    y_name = "y axis",
    x_unit = UnitfulAstro.kpc,
    y_unit = UnitfulAstro.kpc,
    x_limits = (-50, 50),
    y_limits = (-50, 50),
    animation = false,
    output_format = ".png",
    resolution = (1000, 1000),
    aspect = AxisAspect(1)
)

# Comparison with experiments
snapshotPlot(
    [BASE_SRC_PATH,],
    [SNAP_NAME,],
    Dict(:stars => ["MASS", "POS"]),
    scatterlines!;
    output_path = joinpath(BASE_OUT_PATH, "snapshotPlot/experiment"),
    sim_cosmo = SIM_COSMO,
    idx = 16:17,
    sim_labels = nothing,
    t_unit = UnitfulAstro.Gyr,
    warnings = false,
    da_functions = [GI.profile,],
    da_args = [((:stars, "MASS"), 20UnitfulAstro.kpc, 10),],
    da_kwargs = [(density = true,),],
    post_processing! = GI.compareExperiment!,
    pp_args = ("logΣ*", UnitfulAstro.kpc, UnitfulAstro.Msun * UnitfulAstro.pc^(-2)),
    pp_kwargs = (source_path = BASE_SRC_PATH,),
    x_name = "r",
    y_name = "Mass",
    x_unit = UnitfulAstro.kpc,
    y_unit = UnitfulAstro.Msun * UnitfulAstro.pc^(-2),
    animation = false,
    output_format = ".png"
)

# SFR time series
timeSeriesPlot(
    [BASE_SRC_PATH,],
    [SNAP_NAME,],
    lines!;
    output_path = joinpath(BASE_OUT_PATH, "timeSeriesPlot"),
    sim_labels = nothing,
    title = "Star formation rate",
    da_functions = [GI.qtyEvolution,],
    da_args = [("clock_time", "sfr", nothing, nothing),],
    da_kwargs = [(warnings = false, sim_cosmo = SIM_COSMO),],
    x_name = "t",
    y_name = "SFR",
    x_unit = Unitful.Myr,
    y_unit = UnitfulAstro.Msun / UnitfulAstro.kyr,
    y_scale = log10,
    file_name = "time_series_sfr.png"
)

# Gas mass time series
timeSeriesPlot(
    [BASE_SRC_PATH,],
    [SNAP_NAME,],
    lines!;
    output_path = joinpath(BASE_OUT_PATH, "timeSeriesPlot"),
    sim_labels = nothing,
    title = "Star formation rate",
    da_functions = [GI.qtyEvolution,],
    da_args = [("clock_time", "mass", nothing, :gas),],
    da_kwargs = [(warnings = false, sim_cosmo = SIM_COSMO),],
    x_name = "t",
    y_name = "m",
    x_unit = Unitful.Myr,
    y_unit = UnitfulAstro.Msun,
    y_scale = log10,
    file_name = "time_series_gas_mass.png"
)

# Star number time series
timeSeriesPlot(
    [BASE_SRC_PATH,],
    [SNAP_NAME,],
    lines!;
    output_path = joinpath(BASE_OUT_PATH, "timeSeriesPlot"),
    sim_labels = nothing,
    title = "Star formation rate",
    da_functions = [GI.qtyEvolution,],
    da_args = [("clock_time", "number", nothing, :stars),],
    da_kwargs = [(warnings = false, sim_cosmo = SIM_COSMO),],
    x_name = "t",
    y_name = "#",
    x_unit = Unitful.Myr,
    y_unit = Unitful.NoUnits,
    y_scale = log10,
    file_name = "time_series_star_number.png"
)

# Snapshot table
snapshotTable(
    [BASE_SRC_PATH,],
    [SNAP_NAME,],
    Dict(:stars => ["MASS", "POS"]),
    human_redable = false,
    output_path = joinpath(BASE_OUT_PATH, "tables"),
    sim_cosmo = SIM_COSMO,
    idx = 16:17,
    warnings = false,
    da_functions = [GI.profile,],
    da_args = [((:stars, "MASS"), 20UnitfulAstro.kpc, 10),],
    da_kwargs = [(density = true,),],
    x_name = "r",
    y_name = "Mass",
    x_unit = UnitfulAstro.kpc,
    y_unit = UnitfulAstro.Msun * UnitfulAstro.pc^(-2),
)

# Time series table
timeSeriesTable(
    [BASE_SRC_PATH,],
    [SNAP_NAME,],
    human_redable = false,
    output_path = joinpath(BASE_OUT_PATH, "tables"),
    title = "Stellar number",
    da_functions = [GI.qtyEvolution,],
    da_args = [("clock_time", "number", nothing, :stars),],
    da_kwargs = [(warnings = false, sim_cosmo = SIM_COSMO),],
    x_name = "t",
    y_name = "Number of stars",
    xunit_label = "x / :auto_label",
    x_unit = Unitful.Myr,
    y_unit = Unitful.NoUnits,
    file_name = "time_series_star_number",
)
