# Data analysis functions

These functions are used internally and none are exported. 

These functions depend on the particulars of the simulation code (e.g. units).

These functions read the data generated by the [data acquisition](https://ezequiel92.github.io/GalaxyInspector/dev/api/analysis/data_acquisition/) functions, and produce the values that will be plotted.

### Signature for the [`plotSnapshot`](@ref) function

A data analysis functions for [`plotSnapshot`](@ref) must take a dictionary with the following shape:

  + `:sim_data`          -> `::Simulation` (see [`Simulation`](https://ezequiel92.github.io/GalaxyInspector/dev/api/constants/globals/#GalaxyInspector.Simulation)).
  + `:snap_data`         -> `::Snapshot` (see [`Snapshot`](https://ezequiel92.github.io/GalaxyInspector/dev/api/constants/globals/#GalaxyInspector.Snapshot)).
  + `:gc_data`           -> `::GroupCatalog` (see [`GroupCatalog`](https://ezequiel92.github.io/GalaxyInspector/dev/api/constants/globals/#GalaxyInspector.GroupCatalog)).
  + `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
  + `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
  + `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
    + ...
  + `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
  + `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
  + `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
  + ...

and return one or more vectors or matrices with the processed data. It should return `nothing` if the input data has some problem that prevents computation (e.g. is empty).

#### Expected signature:

```julia
  da_function(data_dict, args...; kwargs...) -> (processed_data, ...)
```

where:

  - `data_dict::Dict`
  - `processed_data::Union{VecOrMat{<:Number},Nothing}`


### Signature for the [`plotTimeSeries`](@ref) function

A data analysis functions for [`plotTimeSeries`](@ref) must take a [`Simulation`](https://ezequiel92.github.io/GalaxyInspector/dev/api/constants/globals/#GalaxyInspector.Simulation) struct, and return two vectors. It should return `nothing` if the input data has some problem that prevents computation (e.g. is empty).

#### Expected signature:

```julia
  da_function(sim_data, args...; kw_args...) -> (processed_data_x, processed_data_y)
```

where:

  - `sim_data::Simulation`, see [`Simulation`](https://ezequiel92.github.io/GalaxyInspector/dev/api/constants/globals/#GalaxyInspector.Simulation)
  - `processed_data_x::Vector{<:Number}`
  - `processed_data_y::Vector{<:Number}`

---

```@autodocs
Modules = [GalaxyInspector]
Pages   = ["analysis/data_analysis.jl"]
```