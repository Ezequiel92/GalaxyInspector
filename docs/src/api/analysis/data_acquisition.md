# Data acquisition functions

These functions are only used internally, and depend on the particulars of the simulation code (e.g. units).

These function read the different output files from the simulation (snapshots, FoF catalogs, sfr.txt, etc.), and load the data into memory as dictionaries.

Some of these functions are exported.

---

```@autodocs
Modules = [GalaxyInspector]
Pages   = ["analysis/data_acquisition.jl"]
```
