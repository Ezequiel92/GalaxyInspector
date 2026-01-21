cd(@__DIR__)

using Pkg

Pkg.activate("../")
Pkg.instantiate()

using LaTeXStrings, Unitful, UnitfulAstro

push!(LOAD_PATH, "../src/")
using GalaxyInspector

const BASE_OUT_PATH    = "./test_plots"
const BASE_SRC_PATH    = "F:/simulations/current/"
const SIMULATIONS      = ["SFM_06"]
const LABELS           = ["SFM_06"]
const SNAP_N           = 128
const TRANS_MODE       = :stellar_subhalo
const FILTER_MODE      = :subhalo
const SIMULATION_PATHS = joinpath.(BASE_SRC_PATH, SIMULATIONS)

# Create the output folder if it doesn't exist
mkpath(BASE_OUT_PATH)

log_file = open(joinpath(BASE_OUT_PATH, "logs.txt"), "w+")
GalaxyInspector.setLogging!(true; stream=log_file)

# scatterPlot(
#     SIMULATION_PATHS,
#     SNAP_N,
#     :Z_gas_mass,
#     :gas_mass;
#     xlog=true,
#     ylog=true,
#     output_path=joinpath(BASE_OUT_PATH, "scatterPlot"),
#     trans_mode=TRANS_MODE,
#     filter_mode=FILTER_MODE,
# )

# scatterDensityMap(
#     SIMULATION_PATHS,
#     SNAP_N,
#     :stellar_eff,
#     :stellar_mass;
#     xlog=true,
#     ylog=true,
#     output_path=joinpath(BASE_OUT_PATH, "scatterDensityMap"),
#     trans_mode=TRANS_MODE,
#     filter_mode=FILTER_MODE,
#     colorbar=true,
# )

# scatterDensityMap(
#     SIMULATION_PATHS,
#     SNAP_N,
#     :gas_number_density,
#     :temperature,
#     :gas_mass;
#     xlog=true,
#     ylog=true,
#     output_path=joinpath(BASE_OUT_PATH, "scatterDensityMap"),
#     trans_mode=TRANS_MODE,
#     filter_mode=FILTER_MODE,
#     colorbar=true,
# )

# vtkFiles(
#     SIMULATION_PATHS,
#     SNAP_N,
#     :gas_mass;
#     field_type=:cells,
#     density=true,
#     l_unit=u"kpc",
#     box_size=GalaxyInspector.BOX_L,
#     resolution=250,
#     output_path=joinpath(BASE_OUT_PATH, "vtkFiles"),
#     trans_mode=TRANS_MODE,
#     filter_mode=FILTER_MODE,
# )

# vtkFiles(
#     SIMULATION_PATHS,
#     SNAP_N,
#     :stellar_mass;
#     field_type=:particles,
#     density=true,
#     l_unit=u"kpc",
#     box_size=GalaxyInspector.BOX_L,
#     resolution=250,
#     output_path=joinpath(BASE_OUT_PATH, "vtkFiles"),
#     trans_mode=TRANS_MODE,
#     filter_mode=FILTER_MODE,
# )

# radialProfile(
#     SIMULATION_PATHS,
#     SNAP_N,
#     :ode_dust_mass;
#     norm=nothing,
#     radius=GalaxyInspector.DISK_R,
#     shift=1.0u"kpc",
#     n_bins=30,
#     ylog=true,
#     output_path=joinpath(BASE_OUT_PATH, "radialProfile"),
#     trans_mode=TRANS_MODE,
#     filter_mode=FILTER_MODE,
# )

# radialProfile(
#     SIMULATION_PATHS,
#     SNAP_N,
#     [:ode_molecular_fraction, :ode_atomic_fraction],
#     L"\text{Fractions}";
#     radius=GalaxyInspector.DISK_R,
#     ylog=true,
#     density=false,
#     output_path=joinpath(BASE_OUT_PATH, "radialProfile"),
#     trans_mode=TRANS_MODE,
#     filter_mode=FILTER_MODE,
# )

# stellarHistory(
#     SIMULATION_PATHS,
#     SNAP_N;
#     quantity=:sfr,
#     output_path=joinpath(BASE_OUT_PATH, "stellarHistory"),
#     trans_mode=TRANS_MODE,
#     filter_mode=FILTER_MODE,
# )

# histogram(
#     SIMULATION_PATHS,
#     SNAP_N,
#     :ode_molecular_mass;
#     n_bins=30,
#     line=true,
#     norm=0,
#     range=nothing,
#     xlog=true,
#     output_path=joinpath(BASE_OUT_PATH, "histogram"),
#     trans_mode=TRANS_MODE,
#     filter_mode=FILTER_MODE,
# )

# rotationCurve(
#     SIMULATION_PATHS,
#     SNAP_N;
#     R=GalaxyInspector.DISK_R,
#     output_path=joinpath(BASE_OUT_PATH, "rotationCurve"),
#     trans_mode=TRANS_MODE,
#     filter_mode=FILTER_MODE,
# )

# gasBarPlot(
#     SIMULATION_PATHS,
#     SNAP_N,
#     :gas_metallicity,
#     [1.0e-2, 1.0e-1, 1.0, 1.0e1];
#     components=[:ode_ionized, :ode_atomic, :ode_molecular_stellar],
#     ylog=true,
#     output_path=joinpath(BASE_OUT_PATH, "gasBarPlot"),
#     trans_mode=TRANS_MODE,
#     filter_mode=FILTER_MODE,
# )

# compareMolla2015(
#     SIMULATION_PATHS,
#     SNAP_N,
#     :O_stellar_abundance;
#     output_path=joinpath(BASE_OUT_PATH, "compareMolla2015"),
#     trans_mode=TRANS_MODE,
#     filter_mode=FILTER_MODE,
# )

# compareAgertz2021(
#     SIMULATION_PATHS,
#     SNAP_N;
#     ylog=true,
#     output_path=joinpath(BASE_OUT_PATH, "compareAgertz2021"),
#     trans_mode=TRANS_MODE,
#     filter_mode=FILTER_MODE,
# )

# densityMap(
#     SIMULATION_PATHS,
#     SNAP_N;
#     components=[:gas, :stellar],
#     field_types=[:cells, :particles],
#     projection_planes=[:xy, :xz, :yz],
#     box_size=GalaxyInspector.BOX_L,
#     pixel_length=GalaxyInspector.BOX_L/300.0,
#     m_unit=u"Msun",
#     l_unit=u"kpc",
#     output_path=joinpath(BASE_OUT_PATH, "densityMap"),
#     trans_mode=TRANS_MODE,
#     filter_mode=FILTER_MODE,
# )

# densityMapVelField(
#     SIMULATION_PATHS,
#     SNAP_N;
#     components=[:gas, :stellar],
#     field_types=[:cells, :particles],
#     projection_planes=[:xy, :xz, :yz],
#     box_size=GalaxyInspector.BOX_L,
#     pixel_length=GalaxyInspector.BOX_L/30.0,
#     output_path=joinpath(BASE_OUT_PATH, "densityMapVelField"),
#     trans_mode=TRANS_MODE,
#     filter_mode=FILTER_MODE,
# )

# gasSFRMap(
#     SIMULATION_PATHS,
#     SNAP_N;
#     field_type=:cells,
#     projection_planes=[:xy, :xz, :yz],
#     box_size=GalaxyInspector.BOX_L,
#     pixel_length=GalaxyInspector.BOX_L/300.0,
#     output_path=joinpath(BASE_OUT_PATH, "gasSFRMap"),
#     trans_mode=TRANS_MODE,
#     filter_mode=FILTER_MODE,
# )

# metallicityMap(
#     SIMULATION_PATHS,
#     SNAP_N;
#     components=[:gas, :stellar],
#     field_types=[:cells, :particles],
#     element=:all,
#     projection_planes=[:xy, :xz, :yz],
#     box_size=GalaxyInspector.BOX_L,
#     pixel_length=GalaxyInspector.BOX_L/300.0,
#     output_path=joinpath(BASE_OUT_PATH, "metallicityMap"),
#     trans_mode=TRANS_MODE,
#     filter_mode=FILTER_MODE,
# )

# timeSeries(
#     SIMULATION_PATHS,
#     :physical_time,
#     :ode_molecular_mass;
#     xlog=false,
#     output_path=joinpath(BASE_OUT_PATH, "timeSeries"),
#     trans_mode=TRANS_MODE,
#     filter_mode=FILTER_MODE,
# )

# timeSeries(
#     SIMULATION_PATHS,
#     :physical_time,
#     :sfr;
#     xlog=false,
#     output_path=joinpath(BASE_OUT_PATH, "timeSeries"),
#     trans_mode=TRANS_MODE,
#     filter_mode=FILTER_MODE,
# )

# timeSeries(
#     SIMULATION_PATHS,
#     :physical_time,
#     :halo_R200_1;
#     xlog=false,
#     output_path=joinpath(BASE_OUT_PATH, "timeSeries"),
#     trans_mode=TRANS_MODE,
#     filter_mode=FILTER_MODE,
# )

# statisticsEvolution(
#     SIMULATION_PATHS,
#     :physical_time,
#     :ode_molecular_mass;
#     ylog=true,
#     output_path=joinpath(BASE_OUT_PATH, "statisticsEvolution"),
#     trans_mode=TRANS_MODE,
#     filter_mode=FILTER_MODE,
# )

# gasEvolution(
#     SIMULATION_PATHS;
#     output_path=joinpath(BASE_OUT_PATH, "gasEvolution"),
#     trans_mode=TRANS_MODE,
#     filter_mode=FILTER_MODE,
# )

gasFractionsEvolution(
    SIMULATION_PATHS;
    output_path=joinpath(BASE_OUT_PATH, "gasFractionsEvolution"),
    trans_mode=TRANS_MODE,
    filter_mode=FILTER_MODE,
)

sfrTXT(
    SIMULATION_PATHS,
    :physical_time,
    :stellar_mass;
    output_path=joinpath(BASE_OUT_PATH, "sfrTXT"),
    trans_mode=TRANS_MODE,
    filter_mode=FILTER_MODE,
)

cpuTXT(
    SIMULATION_PATHS,
    "total",
    :physical_time,
    :clock_time_s;
    ylog=true,
    output_path=joinpath(BASE_OUT_PATH, "cpuTXT"),
)

# circularityHistogram(
#     ["F:/simulations/current/SFM_06"],
#     128;
#     R_in=2.0u"kpc",
#     R_out=GalaxyInspector.DISK_R,
#     output_path="./test/new",
#     trans_mode=:stellar_subhalo,
#     filter_mode=:subhalo,
# )

# efficiencyHistogram(
#     ["F:/simulations/current/SFM_06"],
#     128;
#     range=(1.0e-4, 1.0),
#     output_path="./test/new",
#     trans_mode=:stellar_subhalo,
#     filter_mode=:subhalo,
# )


# massProfile(
#     ["F:/simulations/current/SFM_06"],
#     128,
#     [:stellar, :ode_molecular];
#     cumulative=false,
#     ylog=true,
#     radius=GalaxyInspector.DISK_R,
#     n_bins=50,
#     output_path="./test/new",
#     trans_mode=:stellar_subhalo,
#     filter_mode=:subhalo,
# )

# velocityProfile(
#     ["F:/simulations/current/SFM_06"],
#     128,
#     :stellar_radial_velocity;
#     radius=GalaxyInspector.DISK_R,
#     n_bins=40,
#     output_path="./test/new",
#     trans_mode=:stellar_subhalo,
#     filter_mode=:subhalo,
# )

# compareFeldmann2020(
#     ["F:/simulations/current/SFM_06"],
#     :stellar,
#     :sfr;
#     slice=128,
#     scatter=true,
#     xlog=true,
#     ylog=true,
#     output_path="./test/new",
#     trans_mode=:stellar_subhalo,
#     filter_mode=:subhalo,
# )

# fitVSFLaw(
#     ["F:/simulations/current/SFM_06"],
#     128,
#     :ode_molecular;
#     field_type=:cells,
#     fit=true,
#     box_size=GalaxyInspector.BOX_L,
#     x_range=(-Inf, Inf),
#     output_path="./test/new/vsf",
#     trans_mode=:stellar_subhalo,
#     filter_mode=:subhalo,
# )


# SDSSMockup(
#     ["F:/simulations/current/SFM_06"],
#     128;
#     output_path="./test/new/sdss",
#     box_size=GalaxyInspector.BOX_L,
#     resolution=800,
#     projection_plane=:xz,
#     smooth=false,
#     trans_mode=:stellar_subhalo,
#     filter_mode=:subhalo,
# )

# SDSSMockup(
#     ["F:/simulations/current/SFM_06"],
#     128;
#     output_path="./test/new/sdss",
#     box_size=GalaxyInspector.BOX_L,
#     resolution=800,
#     projection_plane=:yz,
#     smooth=false,
#     trans_mode=:stellar_subhalo,
#     filter_mode=:subhalo,
# )
