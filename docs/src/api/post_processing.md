# Post processing functions

These function are supposed to take an already made plot and add some element to it.

## Post-processing functions.


A post-processing function must take a Makie figure, add something to it, and return how to label 
the additions (or `nothing` when no new labels should be drawn).

Expected signature:

  post_processing(figure, args...; kwargs...) -> ([marker, ...], [label, ...]) or `nothing`

where:

  - figure::Makie.Figure
  - marker::LegendElement
  - label::String

---

```@autodocs
Modules = [GalaxyInspector]
Pages   = ["post_processing.jl"]
```
