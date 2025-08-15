# GalaxyInspector

**GalaxyInspector** is a Julia module for the analysis and visualization of galaxy simulation data, with a focus on snapshots produced by the [Arepo](https://arepo-code.org/) code in HDF5 format.

This module provides a collection of scripts and functions to streamline the workflow of reading, filtering, analyzing, and plotting simulation outputs. It is designed for flexibility and extensibility, allowing users to adapt it to their own simulation setups and scientific needs.

## Features

- **Data Acquisition:** Functions to read snapshots, group catalogs, and auxiliary files into convenient Julia data structures.
- **Analysis Pipelines:** Tools for computing physical quantities, filtering data, and transforming simulation outputs for scientific analysis.
- **Plotting:** Ready-to-use and customizable plotting pipelines using [CairoMakie](https://docs.makie.org/stable/) for high-quality figures and animations.
- **Convenience Functions:** Pre-made recipes for common plots and reports, serving as both utilities and usage examples.
- **Extensive Documentation:** API references and usage guides to help users get started and customize their analysis.

## Intended Audience

GalaxyInspector is intended for researchers and students working with galaxy simulations who want a Julia-based, scriptable environment for data analysis and visualization. While it is tailored for Arepo outputs, its modular design allows adaptation to other codes with similar data formats.

> **Note:** This module is a work in progress and was originally developed as a personal learning project. It may lack some advanced features found in other tools (see [links](https://github.com/Ezequiel92/GalaxyInspector#-links)), but aims to provide a clear and hackable codebase for custom workflows.

## License

GalaxyInspector is free software, released under the [GNU General Public License v3.0](https://www.gnu.org/licenses/gpl-3.0.html).
