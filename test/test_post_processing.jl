####################################################################################################
# Testing for the post processing functions in `src/post_processing.jl`
####################################################################################################

@testset "Post processing functions" begin

    figure = Figure()
    axis = Axis(figure[1, 1])
    temp_img = joinpath(@__DIR__, "test_img.png")

    scatter!(axis, 1:10, 1:10)
    vf_r1, vf_r2 = @test_nowarn GI.verticalFlags!(
        figure,
        ([3.0, 5.0, 8.0], ["3.0", "5.0", "8.0"]),
    )
    save(temp_img, figure)
    @test length(vf_r1) == 3
    @test vf_r2 == ["3.0", "5.0", "8.0"]
    @test_reference joinpath(BASE_DATA_PATH, "vf_plot.png") load(temp_img)
    GI.cleanPlot!(figure)

    ks_y = [-0.2, -0.98, -1.21, 1.07, 4.64, 7.69, 7.57, 10.28, 9.08, 10.42, 11.11, 14.54, 15.69,
        17.85, 16.94, 18.56, 18.12, 22.44, 23.46, 24.38]
    scatter!(axis, 1:20, ks_y)
    ks_r1, ks_r2 = @test_nowarn GI.ksLine!(
        figure,
        error_formating = "conf_interval",
        colors = (:orange, :green),
        linestyles = (:solid, :dashdot),
    )
    save(temp_img, figure)
    @test length(ks_r1) == 2
    @test ks_r2 == ["Kennicutt et al. (1998)", "Fit"]
    @test_reference joinpath(BASE_DATA_PATH, "ks_plot.png") load(temp_img)
    GI.cleanPlot!(figure)

    scatter!(axis, 1:10, 1:10, color = :blue)
    @test_nowarn GI.cross!(
        figure,
        color = :red,
        cross = (5.0, 5.0),
    )
    save(temp_img, figure)
    @test_reference joinpath(BASE_DATA_PATH, "cr_plot.png") load(temp_img)
    GI.cleanPlot!(figure)

    ce_r1, ce_r2 = @test_nowarn GI.compareExperiment!(
        figure,
        "ΣHI",
        u"kpc",
        u"Msun * pc^(-2)";
        source_path = BASE_SRC_PATH
    )
    save(temp_img, figure)
    @test length(ce_r1) == 1
    @test ce_r2 == ["Mollá et al. (2015)"]
    @test_reference joinpath(BASE_DATA_PATH, "ce_plot.png") load(temp_img)

    rm(temp_img)

end
