####################################################################################################
# Tests for the utilities in `src/utilities.jl`
####################################################################################################

gas_z = snap_data_04[:gas]["Z"]
gas_mass = snap_data_04[:gas]["MASS"]
gas_pos = snap_data_01[:gas]["POS"]
gas_u = snap_data_03[:gas]["U"]
gas_ne = snap_data_03[:gas]["NE"]
stellar_age = snap_data_02[:stars]["AGE"]
stellar_mass = snap_data_04[:stars]["MASS"]
total_star_mass = GI.computeTimeSeries(snapshots["snap_paths"], "mass", :stars, false)

@testset "GADGET utilities" begin

    @test GI.passAll(snapshot, :stars) == 1:2940
    @test GI.energyIntegrand(cosmological_header, 0.53) ≈ 16.225984
    @test GI.internalUnits("POS", isolated_header) ≈ 1.00000013UnitfulAstro.kpc
    @test GI.internalUnits("TEMP", isolated_header) == u"K"
    @test GI.internalUnits("ID", isolated_header) == Unitful.NoUnits
    @test GI.internalUnits("FMOL", isolated_header) == Unitful.NoUnits
    @test GI.internalUnits("SFRTXT_COL5", isolated_header) ≈ 1.00029678e10UnitfulAstro.Msun
    @test_nowarn GI.checkDataShape(heatmap!, 3, 1)
    @test_nowarn GI.checkDataShape(hist!, 1, 23)
    @test_nowarn GI.checkDataShape(scatterlines!, 2, 23)

    jldopen(joinpath(BASE_DATA_PATH, "utilities.jld2"), "r") do file

        @test GI.compare(
            GI.passCritRho(
                snapshot,
                type = :gas,
                crit_ρ = GI.CRITICAL_DENSITY,
                sim_cosmo = SIM_COSMO,
            ), 
            file["pass_crit_rho"],
        )
        
        @test GI.compare(GI.passPositiveQty(snapshot), file["pass_positive_qty"])
        @test GI.compare(GI.passMetallicity(snapshot, 1e-6, 1e-3), file["pass_metallicity"])
        @test GI.compare(total_star_mass, file["total_star_mass"])
        @test GI.compare(GI.metalMass(gas_z), file["metal_mass"])
        @test GI.compare(GI.computeMetallicity(gas_z, gas_mass, solar = true), file["metallicity"])
        @test GI.compare(GI.computeElementFraction(gas_z, "O"), file["element_fraction"])
        @test GI.compare(GI.computeTime([0:0.01:1...], cosmological_header), file["time"])
        @test GI.compare(GI.computeRedshift([0:0.01:1...]), file["redshift"])
        @test GI.compare(GI.computeStellarAge(stellar_age, sim_data, snap_data), file["stellar_age"])
        @test GI.compare(GI.computeSFR(total_star_mass, file["time"]), file["sfr"])
        @test GI.compare(GI.computeDistance(gas_pos), file["distance"])
        @test GI.compare(GI.computeCenterOfMass(gas_pos, gas_mass), file["center_of_mass"])
        @test GI.compare(GI.computeTemperature(gas_z, gas_mass, gas_u, gas_ne), file["temperature"])

        @test GI.compare(
            GI.computeSurfaceDensity(file["distance"], gas_mass, BOX_SIZE, 100), 
            file["surface_density_disc"],
        )

        @test GI.compare(
            GI.computeSurfaceDensity(gas_pos, gas_mass, BOX_SIZE, 100), 
            file["surface_density_square"],
        )

        @test GI.compare(GI.computeProfile(file["distance"], gas_mass, BOX_SIZE, 100), file["basic_profile"])

        @test GI.compare(
            GI.computeProfile(
                file["distance"],
                gas_mass,
                BOX_SIZE,
                100,
                flat = false,
                cumulative = true,
                density = true,
            ),
            file["full_profile"],
        )

        @test GI.compare(
            GI.computeZProfile(file["distance"],
                gas_z,
                gas_mass,
                BOX_SIZE,
                100,
            ),
            file["basic_z_profile"],
        )

        @test GI.compare(
            GI.computeZProfile(
                file["distance"],
                gas_z,
                gas_mass,
                BOX_SIZE,
                100,
                flat = false,
                cumulative = true,
                density = true,
                solar = true,
            ),
            file["full_z_profile"],
        )

        @test GI.compare(
            GI.computeTimeSeries(
                snapshots["snap_paths"],
                "clock_time",
                :stars,
                false,
            ),
            file["time_series_clock_time"],
        )

        @test GI.compare(
            GI.computeTimeSeries(
                snapshots["snap_paths"],
                "number",
                :stars,
                false,
            ),
            file["time_series_number"],
        )

        @test GI.compare(
            GI.computeTimeSeries(
                snapshots["snap_paths"],
                "sfr",
                :stars,
                false,
            ),
            file["time_series_sfr"],
        )

    end

end

@testset "Generic utilities" begin

    @test GI.trivial(6.9, "nice") === nothing
    @test typeof(
        GI.read_csv(joinpath(BASE_DATA_PATH, "tables/snapshot_150.csv"))
    ) == Matrix{Float64}
    @test GI.isempty(L"")

    figure = Figure()
    axis = Axis(figure[1, 1])
    scatter!(axis, 12:21, 77:86)
    plot = scatter(1:20, 33:52)
    @test all(GI.xlimits(plot) .≈ (0.05, 20.95))
    @test all(GI.xlimits(axis) .≈ (11.55f0, 21.45f0))
    @test all(GI.xlimits(figure) .≈ (11.55f0, 21.45f0))
    @test all(GI.ylimits(plot) .≈ (32.05f0, 52.95f0))
    @test all(GI.ylimits(axis) .≈ (76.55f0, 86.455f0))
    @test all(GI.ylimits(figure) .≈ (76.55f0, 86.455f0))
    @test GI.xscale(plot) == identity
    @test GI.xscale(axis) == identity
    @test GI.xscale(figure) == identity
    @test GI.yscale(plot) == identity
    @test GI.yscale(axis) == identity
    @test GI.yscale(figure) == identity

    point_1 = [Float32.([x, y]) for (x, y) in zip(1:20, 33:52)]
    point_2 = [Float32.([x, y]) for (x, y) in zip(12:21, 77:86)]
    @test GI.pointData(plot) == point_1
    @test GI.pointData(axis) == point_2
    @test GI.pointData(figure) == point_2

    @test all(GI.absCoor(plot, 0.5, 0.5) .≈ (10.4999998, 42.4999981))
    @test all(GI.absCoor(figure, 0.5, 0.5) .≈ (16.5000005, 81.500004))

    @test_nowarn GI.cleanPlot!(figure)

    @test GI.getUnitLabel(0, Unitful.NoUnits) == ""
    @test GI.getUnitLabel(5, Unitful.NoUnits) == L"10^{5}"
    @test GI.getUnitLabel(0, u"kg") == L"\mathrm{kg}"
    @test GI.getUnitLabel(5, u"kg") == L"10^{5} \, \mathrm{kg}"
    @test GI.getUnitLabel(0, Unitful.NoUnits, latex = false) == "dimensionless"
    @test GI.getUnitLabel(5, Unitful.NoUnits, latex = false) == "10^5"
    @test GI.getUnitLabel(0, u"kg", latex = false) == "kg"
    @test GI.getUnitLabel(5, u"kg", latex = false) == "10^5 kg"

    @test GI.getLabel("x axis", 0, Unitful.NoUnits) == LaTeXString("x axis")
    @test GI.getLabel("x axis", 5, Unitful.NoUnits) == L"x axis$\; / \;$$10^{5}$"
    @test GI.getLabel("x axis", 0, u"kg") == L"x axis$\; / \;$$\mathrm{kg}$"
    @test GI.getLabel("x axis", 5, u"kg") == L"x axis$\; / \;$$10^{5} \, \mathrm{kg}$"
    @test GI.getLabel("x axis", 0, Unitful.NoUnits, latex = false) == "x axis"
    @test GI.getLabel("x axis", 5, Unitful.NoUnits, latex = false) == "x axis / 10^5"
    @test GI.getLabel("x axis", 0, u"kg", latex = false) == "x axis / kg"
    @test GI.getLabel("x axis", 5, u"kg", latex = false) == "x axis / 10^5 kg"

    sb = @test_nowarn GI.scaledBins([1.2:12e6...], 5, log10, limits = (12, 12e5))
    @test all(sb .≈ [1.2 * 10^i for i in 1:6])

    sw = @test_nowarn GI.smoothWindow([1:5e5...], rand(500000), 3, scale = log10)
    @test all(.≈(sw, ([40.0, 3189.5, 253149.5], [0.5, 0.5, 0.5]); atol = 0.1))

    @test GI.formatError(69.42069, 0.0) == (69.42069, 0.0)
    @test GI.formatError(69.42069, 0.0387) == (69.42, 0.04)
    @test GI.formatError(69.42069, 0.0187) == (69.421, 0.019)
    @test GI.formatError(69.42069, 0.0103) == (69.421, 0.01)
    @test GI.formatError(0.069420, 0.0903) == (0.07, 0.09)
    @test GI.formatError(0.006942, 0.0183) == (0.007, 0.018)
    @test GI.formatError(0.000426, 0.0183) == (0.0, 0.018)
    @test GI.formatError(69.42069, 1.56) == (69.4, 1.6)
    @test GI.formatError(69.42069, 12.8) == (69, 13)
    @test GI.formatError(69.42069, 36.4) == (70, 40)
    @test GI.formatError(69.42069, 90.3) == (70.0, 90.0)
    @test GI.formatError(0.694206, 1.83) == (0.7, 1.8)
    @test GI.formatError(0.042069, 1.83) == (0.0, 1.8)

    @test GI.maxLength([1 2 3; 4 5 6; 7 8 9]) ≈ 11.224972
    @test GI.maxLength([7 8 9; 4 5 6; 1 2 3]u"kg") ≈ 11.224972u"kg"

    arr1 = [1.5, 9.6, 6.4]
    arr2 = [1.5, 9.6, 6.9]
    df1 = DataFrame(Col1 = ["A", "B", "C"], Col2 = ["X", "Y", "Z"])
    df2 = DataFrame(Col1 = ["A", "B", "D"], Col2 = ["X", "Y", "Z"])
    dict1 = Dict("a" => arr1, "b" => arr1)
    dict2 = Dict("a" => arr1, "b" => arr2)
    missing_arr1 = [1.5, 9.6, 6.4, missing]
    @test !GI.compare("test_1", "test_2")
    @test GI.compare(arr1, arr1)
    @test GI.compare(missing_arr1, missing_arr1)
    @test !GI.compare(arr1, arr2)
    @test GI.compare(arr1, arr1)
    @test GI.compare(df1, df1)
    @test !GI.compare(dict1, dict2)
    @test GI.compare(dict1, dict1)

    @test GI.inRange([42:968...], 41, 969)
    @test GI.inRange([42:969...], 41, 969, strict = false)

    a, b, c = rand(10), rand(12), rand(13)
    @test GI.longest([a, b, c]) == c
    @test GI.longest([a * u"kg", b * u"kg", c * u"kg"]) == c * u"kg"

    rc = [12:169...]
    @test GI.rangeCut!(rc, (14, 150)) == 0
    @test rc == [14:150...]
    @test GI.rangeCut!(rc, (14, 150), keep_edges = false) == 0
    @test rc == [15:149...]
    @test GI.rangeCut!(rc, (1, 10), min_left = 1) == 1
    @test rc == [15:149...]

    m_rc = [12:169...]
    s_rc = [62:219...]
    @test GI.rangeCut!(m_rc, s_rc, (14, 150)) == 0
    @test m_rc == [14:150...]
    @test s_rc == [64:200...]
    @test GI.rangeCut!(m_rc, s_rc, (14, 150), keep_edges = false) == 0
    @test m_rc == [15:149...]
    @test s_rc == [65:199...]
    @test GI.rangeCut!(m_rc, s_rc, (1, 10), min_left = 1) == 1
    @test m_rc == [15:149...]
    @test s_rc == [65:199...]

    pc = [-42:169...]
    pc2 = [-69:-42...]
    @test GI.positiveCut!(pc) == 0
    @test pc == [0:169...]
    @test GI.positiveCut!(pc, keep_edges = false) == 0
    @test pc == [1:169...]
    @test GI.positiveCut!(pc2, min_left = 1) == 1
    @test pc2 == [-69:-42...]

    m_pc = [-42:169...]
    m_pc2 = [-69:-42...]
    s_pc = [1:212...]
    s_pc2 = [1:28...]
    @test GI.positiveCut!(m_pc, s_pc) == 0
    @test m_pc == [0:169...]
    @test s_pc == [43:212...]
    @test GI.positiveCut!(m_pc, s_pc, keep_edges = false) == 0
    @test m_pc == [1:169...]
    @test s_pc == [44:212...]
    @test GI.positiveCut!(m_pc2, s_pc2, min_left = 1) == 1
    @test m_pc2 == [-69:-42...]
    @test s_pc2 == [1:28...]

    sd = [-42:169.0...]
    sd2 = [-69:-42.0...]
    @test GI.sanitizeData!(sd) == (0, 0)
    @test sd == [-42:169.0...]
    @test GI.sanitizeData!(
        sd,
        func_domain = log10,
        range = (-Inf, 150.0),
        keep_edges = false,
        exp_factor = 1,
    ) == (0, 0)
    @test sd == [0.1:0.1:14.9...]
    @test GI.sanitizeData!(
        sd2,
        func_domain = log10,
        range = (1.0, 150.0),
        keep_edges = false,
        min_left = 1,
    ) == (1, 1)
    @test sd2 == [-69:-42.0...]

    m_sd = [-42:169.0...]
    s_sd = [1:212.0...]
    m_sd2 = [-69:-42...]
    @test GI.sanitizeData!(m_sd, s_sd) == (0, 0, 0, 0)
    @test m_sd == [-42:169.0...]
    @test s_sd == [1:212.0...]
    @test GI.sanitizeData!(
        m_sd,
        s_sd,
        func_domain = (log10, log10),
        range = ((1.0, 150.0), (1.0, 150.0)),
        keep_edges = (false, false),
        exp_factor = (1, 1),
    ) == (0, 0, 0, 0)
    @test m_sd == [0.2:0.1:10.6...]
    @test s_sd == [4.5:0.1:14.9...]
    @test GI.sanitizeData!(
        m_sd2,
        m_sd2,
        func_domain = (log10, log10),
        range = ((1.0, 150.0), (1.0, 150.0)),
        keep_edges = (false, false),
        min_left = 1,
    ) == (1, 1, 1, 1)
    @test m_sd2 == [-69:-42...]

end
