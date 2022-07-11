####################################################################################################
# Tests for the pipeline functions in `src/piepelines.jl`
####################################################################################################

@testset "Pipeline functions        " begin

    out_dir = joinpath(@__DIR__, "test_plots")
    first_image = joinpath(out_dir, "plots/snapshot_150.png")
    second_image = joinpath(out_dir, "plots/snapshot_160.png")

    # Histogram
    @test_nowarn snapshotPlot(
        [BASE_SRC_PATH,],
        [SNAP_NAME,],
        Dict(:gas => "NH"),
        hist!;
        output_path = out_dir,
        sim_cosmo = SIM_COSMO,
        idx = 16:17,
        sim_labels = nothing,
        t_unit = UnitfulAstro.Gyr,
        warnings = false,
        da_functions = [GI.histogram],
        da_args = [(:gas, "NH",)],
        post_processing! = GI.verticalFlags!,
        pp_args = (([0.1, 0.5], ["0.1", "0.5"]),),
        x_name = "atomic fraction",
        y_name = "counts",
        y_scale = Makie.pseudolog10,
        animation = false,
        output_format = ".png"
    )
    @test_reference(
        joinpath(BASE_DATA_PATH, "snapshotPlot/histograms/plots/snapshot_150.png"),
        load(first_image),
    )
    @test_reference(
        joinpath(BASE_DATA_PATH, "snapshotPlot/histograms/plots/snapshot_160.png"),
        load(second_image),
    )

    rm(joinpath(out_dir, "plots"), recursive = true)

    # Scatter plot
    @test_nowarn snapshotPlot(
        [BASE_SRC_PATH, BASE_SRC_PATH,],
        [SNAP_NAME, SNAP_NAME,],
        Dict(:gas => "POS"),
        scatter!;
        output_path = out_dir,
        sim_labels = ["example_source_data", "example_source_data"],
        sim_cosmo = SIM_COSMO,
        idx = 16:17,
        t_unit = UnitfulAstro.Gyr,
        warnings = false,
        da_functions = [GI.correlation],
        da_args = [((:gas, "POS"), (:gas, "POS"))],
        da_kwargs = [(x_function = x -> x[1, :], y_function = y -> y[2, :])],
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
    @test_reference(
        joinpath(BASE_DATA_PATH, "snapshotPlot/scatter_plots/plots/snapshot_150.png"),
        load(first_image),
    )
    @test_reference(
        joinpath(BASE_DATA_PATH, "snapshotPlot/scatter_plots/plots/snapshot_160.png"),
        load(second_image),
    )

    rm(joinpath(out_dir, "plots"), recursive = true)

    # Scatterlines plot
    @test_nowarn snapshotPlot(
        [BASE_SRC_PATH,],
        [SNAP_NAME,],
        Dict(:gas => ["MASS", "POS"]),
        scatterlines!;
        output_path = out_dir,
        sim_labels = ["example_source_data",],
        sim_cosmo = SIM_COSMO,
        idx = 16:17,
        t_unit = UnitfulAstro.Gyr,
        warnings = false,
        da_functions = [GI.profile],
        da_args = [((:gas, "MASS"), 300UnitfulAstro.kpc, 100)],
        da_kwargs = [(cumulative = true,)],
        x_name = "r",
        y_name = "Mass",
        x_unit = UnitfulAstro.kpc,
        y_unit = UnitfulAstro.Msun,
        y_factor = 10,
        animation = false,
        output_format = ".png"
    )
    @test_reference(
        joinpath(BASE_DATA_PATH, "snapshotPlot/scatterlines/plots/snapshot_150.png"),
        load(first_image),
    )
    @test_reference(
        joinpath(BASE_DATA_PATH, "snapshotPlot/scatterlines/plots/snapshot_160.png"),
        load(second_image),
    )

    rm(joinpath(out_dir, "plots"), recursive = true)

    # Heatmap
    @test_nowarn snapshotPlot(
        [BASE_SRC_PATH,],
        [SNAP_NAME,],
        Dict(:stars => "POS"),
        heatmap!;
        output_path = out_dir,
        sim_cosmo = SIM_COSMO,
        idx = 16:17,
        sim_labels = nothing,
        t_unit = UnitfulAstro.Gyr,
        warnings = false,
        da_functions = [GI.particleMap],
        da_args = [("XY", :stars, BOX_SIZE * 0.25)],
        da_kwargs = [(sim_cosmo = SIM_COSMO,)],
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
    @test_reference(
        joinpath(BASE_DATA_PATH, "snapshotPlot/heatmaps/plots/snapshot_150.png"),
        load(first_image),
    )
    @test_reference(
        joinpath(BASE_DATA_PATH, "snapshotPlot/heatmaps/plots/snapshot_160.png"),
        load(second_image),
    )

    rm(joinpath(out_dir, "plots"), recursive = true)

    # Comparison with experiments
    @test_nowarn snapshotPlot(
        [BASE_SRC_PATH,],
        [SNAP_NAME,],
        Dict(:stars => ["MASS", "POS"]),
        scatterlines!;
        output_path = out_dir,
        sim_cosmo = SIM_COSMO,
        idx = 16:17,
        sim_labels = nothing,
        t_unit = UnitfulAstro.Gyr,
        warnings = false,
        da_functions = [GI.profile],
        da_args = [((:stars, "MASS"), 20UnitfulAstro.kpc, 10)],
        da_kwargs = [(density = true,)],
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
    @test_reference(
        joinpath(BASE_DATA_PATH, "snapshotPlot/experiment/plots/snapshot_150.png"),
        load(first_image),
    )
    @test_reference(
        joinpath(BASE_DATA_PATH, "snapshotPlot/experiment/plots/snapshot_160.png"),
        load(second_image),
    )

    rm(joinpath(out_dir, "plots"), recursive = true)

    # SFR time series
    @test_nowarn timeSeriesPlot(
        [BASE_SRC_PATH,],
        [SNAP_NAME,],
        lines!;
        output_path = out_dir,
        sim_labels = nothing,
        title = "Star formation rate",
        da_functions = [GI.qtyEvolution],
        da_args = [("clock_time", "sfr", nothing, nothing)],
        da_kwargs = [(warnings = false, sim_cosmo = SIM_COSMO)],
        x_name = "t",
        y_name = "SFR",
        x_unit = Unitful.Myr,
        y_unit = UnitfulAstro.Msun / UnitfulAstro.kyr,
        y_scale = log10,
        file_name = "time_series_sfr.png"
    )
    @test_reference(
        joinpath(BASE_DATA_PATH, "timeSeriesPlot/time_series_sfr.png"),
        load(joinpath(out_dir, "time_series_sfr.png")),
    )

    # Gas mass time series
    @test_nowarn timeSeriesPlot(
        [BASE_SRC_PATH,],
        [SNAP_NAME,],
        lines!;
        output_path = out_dir,
        sim_labels = nothing,
        title = "Star formation rate",
        da_functions = [GI.qtyEvolution],
        da_args = [("clock_time", "mass", nothing, :gas)],
        da_kwargs = [(warnings = false, sim_cosmo = SIM_COSMO)],
        x_name = "t",
        y_name = "m",
        x_unit = Unitful.Myr,
        y_unit = UnitfulAstro.Msun,
        y_scale = log10,
        file_name = "time_series_gas_mass.png"
    )
    @test_reference(
        joinpath(BASE_DATA_PATH, "timeSeriesPlot/time_series_gas_mass.png"),
        load(joinpath(out_dir, "time_series_gas_mass.png")),
    )

    # Star number time series
    @test_nowarn timeSeriesPlot(
        [BASE_SRC_PATH,],
        [SNAP_NAME,],
        lines!;
        output_path = out_dir,
        sim_labels = nothing,
        title = "Stellar number",
        da_functions = [GI.qtyEvolution],
        da_args = [("clock_time", "number", nothing, :stars)],
        da_kwargs = [(warnings = false, sim_cosmo = SIM_COSMO,)],
        x_name = "t",
        y_name = "Number of stars",
        x_unit = Unitful.Myr,
        y_unit = Unitful.NoUnits,
        y_scale = log10,
        file_name = "time_series_star_number.png"
    )
    @test_reference(
        joinpath(BASE_DATA_PATH, "timeSeriesPlot/time_series_star_number.png"),
        load(joinpath(out_dir, "time_series_star_number.png")),
    )

    # Snapshot table
    @test_nowarn snapshotTable(
        [BASE_SRC_PATH,],
        [SNAP_NAME,],
        Dict(:stars => ["MASS", "POS"]),
        human_redable = false,
        output_path = out_dir,
        sim_cosmo = SIM_COSMO,
        idx = 16:17,
        warnings = false,
        da_functions = [GI.profile],
        da_args = [((:stars, "MASS"), 20UnitfulAstro.kpc, 10)],
        da_kwargs = [(density = true,)],
        x_name = "r",
        y_name = "Mass",
        x_unit = UnitfulAstro.kpc,
        y_unit = UnitfulAstro.Msun * UnitfulAstro.pc^(-2),
    )
    ref_table_snap150 = GI.read_csv(joinpath(BASE_DATA_PATH, "tables/snapshot_150.csv"))
    ref_table_snap160 = GI.read_csv(joinpath(BASE_DATA_PATH, "tables/snapshot_160.csv"))
    @test all(ref_table_snap150 .≈ GI.read_csv(joinpath(out_dir, "snapshot_150.csv")))
    @test all(ref_table_snap160 .≈ GI.read_csv(joinpath(out_dir, "snapshot_160.csv")))

    # Time series table
    @test_nowarn timeSeriesTable(
        [BASE_SRC_PATH,],
        [SNAP_NAME,],
        human_redable = false,
        output_path = out_dir,
        title = "Stellar number",
        da_functions = [GI.qtyEvolution],
        da_args = [("clock_time", "number", nothing, :stars)],
        da_kwargs = [(warnings = false, sim_cosmo = SIM_COSMO)],
        x_name = "t",
        y_name = "Number of stars",
        xunit_label = "x / :auto_label",
        x_unit = Unitful.Myr,
        y_unit = Unitful.NoUnits,
        file_name = "time_series_star_number",
    )
    ref_table_sn = GI.read_csv(joinpath(BASE_DATA_PATH, "tables/time_series_star_number.csv"))
    @test all(ref_table_sn .≈ GI.read_csv(joinpath(out_dir, "time_series_star_number.csv")))

    rm(out_dir, recursive = true)

end