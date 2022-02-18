####################################################################################################
# Reference run of the post processing functions in `src/post_processing.jl`
####################################################################################################

figure = Figure()
axis = Axis(figure[1, 1])
scatter!(axis, 1:10, 1:10)
GI.verticalFlags!(figure, ([3.0, 5.0, 8.0], ["3.0", "5.0", "8.0"]))
save(joinpath(BASE_OUT_PATH, "vf_plot.png"), figure)
GI.cleanPlot!(figure)

ks_y = [-0.2, -0.98, -1.21, 1.07, 4.64, 7.69, 7.57, 10.28, 9.08, 10.42, 11.11, 14.54, 15.69,
    17.85, 16.94, 18.56, 18.12, 22.44, 23.46, 24.38]
scatter!(axis, 1:20, ks_y)
GI.ksLine!(
    figure,
    error_formating = "conf_interval",
    colors = (:orange, :green),
    linestyles = (:solid, :dashdot),
)
save(joinpath(BASE_OUT_PATH, "ks_plot.png"), figure)
GI.cleanPlot!(figure)

scatter!(axis, 1:10, 1:10, color = :blue)
GI.cross!(figure, color = :red, cross = (5.0, 5.0))
save(joinpath(BASE_OUT_PATH, "cr_plot.png"), figure)
GI.cleanPlot!(figure)

GI.compareExperiment!(
    figure,
    "ΣHI",
    u"kpc",
    u"Msun * pc^(-2)";
    source_path = BASE_SRC_PATH
)
save(joinpath(BASE_OUT_PATH, "ce_plot.png"), figure)
