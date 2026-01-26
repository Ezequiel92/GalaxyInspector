cd(@__DIR__)

using Pkg

Pkg.activate("../")
Pkg.instantiate()

using CairoMakie, LaTeXStrings, Unitful, UnitfulAstro

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
#     shift=0.5u"kpc",
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
#     "f";
#     n_bins=30,
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
#     pixel_length=GalaxyInspector.BOX_L/300.0,
#     output_path=joinpath(BASE_OUT_PATH, "densityMapVelField"),
#     trans_mode=TRANS_MODE,
#     filter_mode=FILTER_MODE,
# )

# gasSFRMap(
#     SIMULATION_PATHS,
#     SNAP_N;
#     field_type=:cells,
#     projection_planes=[:xy, :xz, :yz]
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
#     pixel_length=GalaxyInspector.BOX_L/300.0,
#     output_path=joinpath(BASE_OUT_PATH, "metallicityMap"),
#     trans_mode=TRANS_MODE,
#     filter_mode=FILTER_MODE,
# )

# timeSeries(
#     SIMULATION_PATHS,
#     :physical_time,
#     :ode_molecular_mass;
#     ylog=true,
#     output_path=joinpath(BASE_OUT_PATH, "timeSeries"),
#     trans_mode=TRANS_MODE,
#     filter_mode=FILTER_MODE,
#     theme=Theme(Axis=(xticks=0:14,),),
# )

# timeSeries(
#     SIMULATION_PATHS,
#     :physical_time,
#     :sfr;
#     ylog=true,
#     output_path=joinpath(BASE_OUT_PATH, "timeSeries"),
#     trans_mode=TRANS_MODE,
#     filter_mode=FILTER_MODE,
#     theme=Theme(Axis=(xticks=0:14,),),
# )

# timeSeries(
#     SIMULATION_PATHS,
#     :physical_time,
#     :halo_R200_1;
#     ylog=true,
#     output_path=joinpath(BASE_OUT_PATH, "timeSeries"),
#     trans_mode=TRANS_MODE,
#     filter_mode=FILTER_MODE,
#     theme=Theme(Axis=(xticks=0:14,),),
# )

# statisticsEvolution(
#     SIMULATION_PATHS,
#     :physical_time,
#     :ode_molecular_mass;
#     ylog=true,
#     output_path=joinpath(BASE_OUT_PATH, "statisticsEvolution"),
#     trans_mode=TRANS_MODE,
#     filter_mode=FILTER_MODE,
#     theme=Theme(Axis=(xticks=0:14,),),
# )

# gasEvolution(
#     SIMULATION_PATHS;
#     ylog=true,
#     output_path=joinpath(BASE_OUT_PATH, "gasEvolution"),
#     trans_mode=TRANS_MODE,
#     filter_mode=FILTER_MODE,
#     theme=Theme(Axis=(limits=(nothing, nothing, -11.5, 1.5),),),
# )

# gasFractionsEvolution(
#     SIMULATION_PATHS;
#     output_path=joinpath(BASE_OUT_PATH, "gasFractionsEvolution"),
#     trans_mode=TRANS_MODE,
#     filter_mode=FILTER_MODE,
# )

# sfrTXT(
#     SIMULATION_PATHS,
#     :physical_time,
#     :sfr;
#     smooth=200,
#     output_path=joinpath(BASE_OUT_PATH, "sfrTXT"),
# )

# cpuTXT(
#     SIMULATION_PATHS,
#     "total",
#     :physical_time,
#     :tot_clock_time_s;
#     smooth=100,
#     ylog=true,
#     output_path=joinpath(BASE_OUT_PATH, "cpuTXT"),
#     theme=Theme(Axis=(xticks=0:14,),),
# )

# kennicuttSchmidtLaw(
#     SIMULATION_PATHS,
#     SNAP_N;
#     quantity=:molecular,
#     reduce_grid=:circular,
#     grid_size=30.0u"kpc",
#     bin_size=1.5u"kpc",
#     post_processing=GalaxyInspector.ppSun2023!,
#     pp_kwargs=(; color=GalaxyInspector.WONG_ORANGE),
#     fit=false,
#     output_file=joinpath(BASE_OUT_PATH, "kennicuttSchmidtLaw/sun2023_molecular.png"),
#     trans_mode=TRANS_MODE,
#     filter_mode=FILTER_MODE,
#     theme=Theme(
#         Legend=(margin=(10, 0, 0, 0),),
#         Axis=(
#             limits=(-2.5, 3.5, -4.5, 0.5),
#             xticks=-1:1:3,
#             yticks=-4:1:0,
#         ),
#     ),
# )

# kennicuttSchmidtLaw(
#     SIMULATION_PATHS,
#     SNAP_N;
#     quantity=:molecular,
#     reduce_grid=:circular,
#     grid_size=30.0u"kpc",
#     bin_size=1.5u"kpc",
#     post_processing=GalaxyInspector.ppLeroy2008!,
#     pp_kwargs=(; color=GalaxyInspector.WONG_ORANGE),
#     fit=false,
#     output_file=joinpath(BASE_OUT_PATH, "kennicuttSchmidtLaw/leroy2008_molecular.png"),
#     trans_mode=TRANS_MODE,
#     filter_mode=FILTER_MODE,
#     theme=Theme(
#         Legend=(margin=(10, 0, 0, 0),),
#         Axis=(
#             limits=(-2.5, 3.5, -4.5, 0.5),
#             xticks=-1:1:3,
#             yticks=-4:1:0,
#         ),
#     ),
# )

# kennicuttSchmidtLaw(
#     SIMULATION_PATHS,
#     SNAP_N;
#     quantity=:gas,
#     reduce_grid=:circular,
#     grid_size=30.0u"kpc",
#     bin_size=1.5u"kpc",
#     post_processing=GalaxyInspector.ppBigiel2010!,
#     pp_kwargs=(; galaxy=:all, quantity=:neutral, color=GalaxyInspector.WONG_ORANGE),
#     fit=false,
#     output_file=joinpath(BASE_OUT_PATH, "kennicuttSchmidtLaw/bigiel2010_total.png"),
#     trans_mode=TRANS_MODE,
#     filter_mode=FILTER_MODE,
#     theme=Theme(
#         Legend=(margin=(10, 0, 0, 20),),
#         Axis=(
#             limits=(0.4, 2.6, -4.5, 0.5),
#             xticks=0.5:0.5:2.5,
#             yticks=-4:1:0,
#         ),
#     ),
# )

# kennicuttSchmidtLaw(
#     SIMULATION_PATHS,
#     SNAP_N;
#     quantity=:gas,
#     reduce_grid=:circular,
#     grid_size=30.0u"kpc",
#     bin_size=1.5u"kpc",
#     post_processing=GalaxyInspector.ppLeroy2008!,
#     pp_kwargs=(; quantity=:neutral, color=GalaxyInspector.WONG_ORANGE),
#     fit=false,
#     output_file=joinpath(BASE_OUT_PATH, "kennicuttSchmidtLaw/leroy2008_total.png"),
#     trans_mode=TRANS_MODE,
#     filter_mode=FILTER_MODE,
#     theme=Theme(
#         Legend=(margin=(10, 0, 0, 20),),
#         Axis=(
#             limits=(0.4, 2.6, -4.5, 0.5),
#             xticks=0.5:0.5:2.5,
#             yticks=-4:1:0,
#         ),
#     ),
# )

# stellarBirthHalos(
#     SIMULATION_PATHS,
#     SNAP_N;
#     output_path=joinpath(BASE_OUT_PATH, "stellarBirthHalos"),
#     trans_mode=TRANS_MODE,
#     filter_mode=FILTER_MODE,
# )

# atomicMolecularTransition(
#     SIMULATION_PATHS,
#     SNAP_N,
#     [(1.0e-2, 1.0e-1), (1.0e-1, 1.0), (1.0, 1.0e1)];
#     output_path=joinpath(BASE_OUT_PATH, "atomicMolecularTransition"),
#     trans_mode=TRANS_MODE,
#     filter_mode=FILTER_MODE,
# )

# massProfile(
#     SIMULATION_PATHS,
#     SNAP_N,
#     [:ode_molecular, :ode_atomic, :ode_ionized];
#     ylog=true,
#     output_path=joinpath(BASE_OUT_PATH, "massProfile"),
#     trans_mode=TRANS_MODE,
#     filter_mode=FILTER_MODE,
# )

# velocityProfile(
#     SIMULATION_PATHS,
#     SNAP_N,
#     :stellar_radial_velocity;
#     output_path=joinpath(BASE_OUT_PATH, "velocityProfile"),
#     trans_mode=TRANS_MODE,
#     filter_mode=FILTER_MODE,
# )

# compareFeldmann2020(
#     SIMULATION_PATHS,
#     :molecular,
#     :atomic;
#     scatter=true,
#     output_path=joinpath(BASE_OUT_PATH, "compareFeldmann2020"),
#     trans_mode=TRANS_MODE,
#     filter_mode=FILTER_MODE,
# )

# compareFeldmann2020(
#     SIMULATION_PATHS,
#     :stellar,
#     :sfr;
#     scatter=true,
#     output_path=joinpath(BASE_OUT_PATH, "compareFeldmann2020"),
#     trans_mode=TRANS_MODE,
#     filter_mode=FILTER_MODE,
# )

# massMetallicityRelation(
#     SIMULATION_PATHS,
#     SNAP_N;
#     output_path=joinpath(BASE_OUT_PATH, "massMetallicityRelation"),
#     trans_mode=TRANS_MODE,
#     filter_mode=FILTER_MODE,
# )

# gasVelocityCubes(
#     SIMULATION_PATHS,
#     SNAP_N;
#     output_file=joinpath(BASE_OUT_PATH, "gasVelocityCubes/gas_velocity_cube.hdf5"),
#     trans_mode=TRANS_MODE,
#     filter_mode=FILTER_MODE,
# )

# stellarVelocityCubes(
#     SIMULATION_PATHS,
#     SNAP_N;
#     output_file=joinpath(BASE_OUT_PATH, "stellarVelocityCubes/stellar_velocity_cube.hdf5"),
#     trans_mode=TRANS_MODE,
#     filter_mode=FILTER_MODE,
# )

# virialAccretionEvolution(
#     SIMULATION_PATHS;
#     flux_direction=:outflow,
#     ylog=true,
#     smooth=80,
#     output_path=joinpath(BASE_OUT_PATH, "virialAccretionEvolution"),
# )

# diskAccretionEvolution(
#     SIMULATION_PATHS;
#     flux_direction=:inflow,
#     trans_mode=TRANS_MODE,
#     ylog=true,
#     smooth=80,
#     output_path=joinpath(BASE_OUT_PATH, "diskAccretionEvolution"),
# )

# fitVSFLaw(
#     SIMULATION_PATHS,
#     SNAP_N,
#     :ode_molecular;
#     output_path=joinpath(BASE_OUT_PATH, "fitVSFLaw"),
#     trans_mode=TRANS_MODE,
#     filter_mode=FILTER_MODE,
# )

# clumpingFactor(
#     SIMULATION_PATHS,
#     SNAP_N,
#     :gas;
#     xlog=true,
#     ylog=true,
#     output_path=joinpath(BASE_OUT_PATH, "clumpingFactor"),
#     trans_mode=TRANS_MODE,
#     filter_mode=FILTER_MODE,
# )

# circularityHistogram(
#     SIMULATION_PATHS,
#     SNAP_N;
#     R_in=2.0u"kpc",
#     output_path=joinpath(BASE_OUT_PATH, "circularityHistogram"),
#     trans_mode=TRANS_MODE,
#     filter_mode=FILTER_MODE,
# )

# efficiencyHistogram(
#     SIMULATION_PATHS,
#     SNAP_N;
#     range=(1.0e-6, 1.0),
#     output_path=joinpath(BASE_OUT_PATH, "efficiencyHistogram"),
#     trans_mode=TRANS_MODE,
#     filter_mode=FILTER_MODE,
# )

# stellarDensityMaps(
#     SIMULATION_PATHS,
#     SNAP_N;
#     output_path=joinpath(BASE_OUT_PATH, "stellarDensityMaps"),
#     trans_mode=TRANS_MODE,
#     filter_mode=FILTER_MODE,
# )

# gasDensityMaps(
#     SIMULATION_PATHS,
#     SNAP_N;
#     output_path=joinpath(BASE_OUT_PATH, "gasDensityMaps"),
#     density_range=(2.0, NaN),
#     trans_mode=TRANS_MODE,
#     filter_mode=FILTER_MODE,
# )

# SDSSMockup(
#     SIMULATION_PATHS,
#     SNAP_N;
#     output_path=joinpath(BASE_OUT_PATH, "SDSSMockup"),
#     resolution=300,
#     projection_plane=:xy,
#     trans_mode=TRANS_MODE,
#     filter_mode=FILTER_MODE,
# )

# SDSSMockup(
#     SIMULATION_PATHS,
#     SNAP_N;
#     output_path=joinpath(BASE_OUT_PATH, "SDSSMockup"),
#     resolution=300,
#     projection_plane=:xz,
#     trans_mode=TRANS_MODE,
#     filter_mode=FILTER_MODE,
# )

# simulationReport(SIMULATION_PATHS; output_path=joinpath(BASE_OUT_PATH, "reports"))

# snapshotReport(
#     SIMULATION_PATHS,
#     [SNAP_N];
#     output_path=joinpath(BASE_OUT_PATH, "reports"),
#     trans_mode=TRANS_MODE,
#     filter_mode=FILTER_MODE,
# )

# quantityReport(
#     SIMULATION_PATHS,
#     :ode_gas_tau_dg;
#     output_path=joinpath(BASE_OUT_PATH, "reports"),
#     trans_mode=TRANS_MODE,
#     filter_mode=FILTER_MODE,
# )

#TODO
evolutionVideo(
    SIMULATION_PATHS,
    :ode_molecular;
    slice=108:128,
    output_path=joinpath(BASE_OUT_PATH, "evolutionVideo"),
    framerate=4,
    trans_mode=TRANS_MODE,
    filter_mode=FILTER_MODE,
)
