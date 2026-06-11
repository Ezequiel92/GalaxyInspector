<div align="center">
    <h1><img src="./docs/src/assets/logo.png"/ style="height: 6em;"></h1>
</div>

<p align="center">
    <a href="https://julialang.org"><img src="https://img.shields.io/badge/-Julia-9558B2?style=for-the-badge&logo=julia&logoColor=white"></a>
</div>

<p align="center">
    <a href="https://github.com/ezequiel92/GalaxyInspector/blob/main/LICENSE"><img src="https://img.shields.io/github/license/ezequiel92/GalaxyInspector?style=flat&logo=GNU&labelColor=2B2D2F"></a>
    <a href="https://ezequiel92.github.io/GalaxyInspector/dev/intro/"><img src="https://img.shields.io/badge/documentation-blue.svg"></a>
</p>

A Julia module for the analysis of galaxy simulations.

> [!CAUTION]
> This is written for my personal use and is a work in progress, thus it may break at any moment. Use it at your own risk.

> ℹ️ **NOTE**
>
> There are other tools to analyze/plot simulations (you can see some [below](https://github.com/Ezequiel92/GalaxyInspector#-links)). This module was written not only as a basic plotting tool, but as an exercise to learn [Julia](https://julialang.org/) and software development in general.

## ⚙️ Characteristics

- Works only with snapshots in [HDF5](https://www.hdfgroup.org/solutions/hdf5/) format (option `SnapFormat = 3` in P-Gadget3 and Arepo).
- This is a collection of scripts inside a module, not a package. Global constants and data structures are defined in `./src/globals/`.

## Configuration

User-tunable defaults can be overridden with a `config.toml` file. The module looks for that file in the current working directory.

Example:

```toml
[mmap]
memory_fraction = 0.2

[galaxy]
disk_r         = "40.0 kpc"
disk_height    = "5.0 kpc"
box_l          = "65.0 kpc"
rotation_r     = "8.0 kpc"
age_resolution = "200.0 Myr"

[abundances]
H     = 0.0
He    = 0.0
C     = 0.0
N     = 0.0
O     = 12.0
Ne    = 0.0
Mg    = 0.0
Si    = 0.0
Fe    = 0.0
Other = 0.0

[arepo]
snap_basename     = "snap"
gc_basename       = "fof_subhalo_tab"
hydrogen_massfrac = 0.76
gamma             = 1.666
tracer_mass       = 3.65456e-06
solar_metallicity = 0.0127
hubble_constant   = 0.102201

[units]
physical_units  = false
internal_l_unit = "3.085678e21 cm"
internal_m_unit = "1.989e43 g"
internal_v_unit = "1.0e5 cm*s^-1"
```

The default values can be found in `./src/globals/config.jl`.

## 🔗 Links

### Arepo

[Arepo](https://arepo-code.org/) ([ascl:1909.010](https://ascl.net/1909.010))

### Gadget

[GADGET2](https://wwwmpa.mpa-garching.mpg.de/gadget/) ([ascl:0003.001](https://ascl.net/0003.001))

[GADGET4](https://wwwmpa.mpa-garching.mpg.de/gadget4/)

### Other visualization and analysis tools

[AMUSE](https://www.amusecode.org/) ([ascl:1107.007](https://ascl.net/1107.007))

[Splash](https://users.monash.edu.au/~dprice/splash/) ([ascl:1103.004](https://ascl.net/1103.004))

[Plonk](https://github.com/dmentipl/plonk) ([ascl:1907.009](https://ascl.net/1907.009))

[yt](https://yt-project.org/) ([ascl:1011.022](https://ascl.net/1011.022))

[pynbody/tangos](https://pynbody.github.io/) ([ascl:1305.002](https://ascl.net/1305.002)/[ascl:1912.018](https://ascl.net/1912.018))

[pygad](https://bitbucket.org/broett/pygad/) ([ascl:1811.014](https://ascl.net/1811.014))

[Firefly](https://github.com/ageller/firefly) ([ascl:1810.021](https://ascl.net/1810.021))

[FIRE studio](https://github.com/agurvich/FIRE_studio) ([ascl:2202.006](https://ascl.net/2202.006))

[GalaXimView](https://vm-weblerma.obspm.fr/~ahalle/galaximview/) ([ascl:code/v/2978](https://ascl.net/code/v/2978))
