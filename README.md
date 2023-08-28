<div align="center">
    <h1><img src="./logo.png"/ style="height: 6em;"></h1>
</div>

<p align="center">
    <a href="https://julialang.org"><img src="https://img.shields.io/badge/-Julia-9558B2?style=for-the-badge&logo=julia&logoColor=white"></a>
</div>

<p align="center">
    <!-- <a href="https://codecov.io/github/ezequiel92/GalaxyInspector?branch=main"><img src="https://img.shields.io/codecov/c/github/ezequiel92/GalaxyInspector?style=flat&logo=Codecov&labelColor=2B2D2F"></a> -->
    <!-- <a href="https://github.com/ezequiel92/GalaxyInspector/actions"><img src="https://img.shields.io/github/workflow/status/ezequiel92/GalaxyInspector/Run%20tests?style=flat&logo=GitHub&labelColor=2B2D2F"></a> -->
    <a href="https://github.com/ezequiel92/GalaxyInspector/blob/main/LICENSE"><img src="https://img.shields.io/github/license/ezequiel92/GalaxyInspector?style=flat&logo=GNU&labelColor=2B2D2F"></a>
    <a href="https://ezequiel92.github.io/GalaxyInspector/dev/intro/"><img src="https://img.shields.io/badge/docs-dev-blue.svg?style=flat&logo=GitBook&labelColor=2B2D2F"></a>
</p>

A Julia module for the data analysis of Arepo hydrodynamical simulations.

## ⚙️ Characteristics

- Works only with snapshots in [HDF5](https://www.hdfgroup.org/solutions/hdf5/) format (option `SnapFormat = 3`).
- This is a collection of scripts inside a module, not a package. Global constants and data structures are defined in `src/constants.jl`.

> ℹ️ **NOTE**: There are more advanced tools to analyze/plot the data generated by Arepo and similar codes (you can see some [below](https://github.com/Ezequiel92/GalaxyInspector/edit/main/README.md#Other-visualization-and-analysis-tools)), but this module was written not only as a basic plotting tool, but as an exercise to learn [Julia](https://julialang.org/) and software development in general (testing, documentation, version control, etc.).

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

> ⚠️ **WARNING**: This code is written for my personal use and is a work in progress, thus it may break at any moment. Use it at your own risk.
