# Post processing functions

A post-processing function must take a [Makie](https://docs.makie.org/stable/) figure, add something to it, and return how to label the additions (or `nothing` when no new labels should be drawn).

None of these functions are exported.

### Expected signature:

```julia
  post_processing(figure, args...; kwargs...) -> ([marker, ...], [label, ...])
```

where:

  - `figure::Makie.Figure`
  - `marker::LegendElement`
  - `label::String`

---

```@autodocs
Modules = [GalaxyInspector]
Pages   = ["post_processing.jl"]
```