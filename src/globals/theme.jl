####################################################################################################
# Default theme for Makie.jl
####################################################################################################

"""
Default colors.
"""
const WONG_BLUE    = Makie.wong_colors()[1]
const WONG_ORANGE  = Makie.wong_colors()[2]
const WONG_GREEN   = Makie.wong_colors()[3]
const WONG_PINK    = Makie.wong_colors()[4]
const WONG_CELESTE = Makie.wong_colors()[5]
const WONG_RED     = Makie.wong_colors()[6]
const WONG_YELLOW  = Makie.wong_colors()[7]

"""
Default markers.
"""
const MARKERS = [
    :circle,
    :rect,
    :dtriangle,
    :utriangle,
    :cross,
    :diamond,
    :ltriangle,
    :rtriangle,
    :pentagon,
    :xcross,
    :hexagon,
]

"""
Default line styles.
"""
const LINE_STYLES = [:solid, :dash, :dot, :dashdot, :dashdotdot]

"""
Default cycler.
"""
const CYCLE = Cycle([:color, :linestyle, :marker], covary=true)

"""
Default plot theme.

Regarding the graphic units used, we know that ``1 \\, \\mathrm{mm} = 2.83466 \\, \\mathrm{pt}`` and ``1 \\, \\mathrm{in} = 25.4 \\, \\mathrm{mm}``. Then, if we want ``1 \\, \\mathrm{[code\\,\\,]unit} = 0.1 \\, \\mathrm{mm}`` in vector graphics, we have to use `pt_per_unit` = 0.283466.

For pixel images, we control the ppi with `px_per_unit`. A reasonable high ppi is 600. So, using `px_per_unit` = ``2.3622`` we get ``23.622 \\, \\mathrm{px/mm} \\sim 600 \\, \\mathrm{px/in}`` (remember that ``1 \\, \\mathrm{[code\\,\\,]unit} = 0.1 \\, \\mathrm{mm}``).
"""
const DEFAULT_THEME = Theme(
    #################################################################################################
    # Size of the figures in code units
    # For PDFs and SVGs, 880 [code ]unit = 8.8 cm
    # For PNGs, when printed to a size of 1 point = 0.1 mm, one will get a dpi of 600 (23.622 px/mm)
    #################################################################################################
    size=(880, 880),
    ######################################
    # 35 unit * 0.283466 pt/unit ~ 9.9 pt
    ######################################
    fontsize=35,
    #############################
    # (left, right, bottom, top)
    #############################
    figure_padding=(2, 15, 5, 15),
    palette=(color=Makie.wong_colors(), marker=MARKERS, linestyle=LINE_STYLES),
    CairoMakie=(px_per_unit=2.3622, pt_per_unit=0.283466),
    Axis=(
        xlabelpadding=15,
        xticklabelpad=10,
        xticksize=7,
        xgridvisible=false,
        spinewidth=3,
        xminorticksvisible=true,
        xminorticks=IntervalsBetween(5),
        ylabelpadding=15,
        yticklabelpad=10,
        yticksize=7,
        ygridvisible=false,
        yminorticksvisible=true,
        yminorticks=IntervalsBetween(5),
        ######################################################################################
        # Aspect ratio of the figures. The options are:
        # nothing: The aspect ratio will be chosen by [Makie](https://docs.makie.org/stable/)
        # AxisAspect(n): The aspect ratio will be given by the number `n` = width / height
        # DataAspect(): The aspect ratio of the data will be used
        ######################################################################################
        aspect=AxisAspect(1),
    ),
    Legend=(
        tellheight=false,
        tellwidth=false,
        framevisible=false,
        colgap=20,
        rowgap=10,
        halign=:right,
        valign=:bottom,
        nbanks=3,
        titlegap=-5,
        labelsize=30,
        linewidth=5,
        markersize=28,
        patchsize=(40, 40),
    ),
    Lines=(linewidth=5, cycle=CYCLE),
    VLines=(linewidth=3, cycle=CYCLE),
    HLines=(linewidth=3, cycle=CYCLE),
    ScatterLines=(
        linewidth=5,
        markersize=22,
        cycle=CYCLE,
    ),
    Scatter=(markersize=22, cycle=CYCLE),
    Band=(alpha=0.5, cycle=CYCLE),
    VSpan=(alpha=0.5, cycle=CYCLE),
    Errorbars=(whiskerwidth=10,),
    ########################################################################
    # Alternative colormaps:
    # colormap = :nipy_spectral - nan_color = ColorSchemes.nipy_spectral[1]
    # colormap = :cubehelix     - nan_color = ColorSchemes.cubehelix[1]
    # colormap = :lipari        - nan_color = ColorSchemes.lipari[1]
    ########################################################################
    Heatmap=(colormap=:CMRmap, nan_color=ColorSchemes.CMRmap[1]),
    Colorbar=(
        colormap=:CMRmap,
        size=25,
        ticklabelpad=10,
        minorticksvisible=true,
        ticksize=7,
        labelpadding=2,
    ),
    BarPlot=(
        color_over_background=:black,
        color_over_bar=:black,
        flip_labels_at=10,
        strokecolor=:black,
        strokewidth=1,
        dodge_gap=0.04,
    ),
    Arrows2D=(lengthscale=0.015, color=:white, shaftwidth=2, tipwidth=8),
    Hist=(strokecolor=:black, strokewidth=1),
)
