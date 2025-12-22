####################################################################################################
# Opinionated convenience functions
####################################################################################################

"""
    scatterPlot(
        simulation_paths::Vector{String},
        slice::IndexType,
        x_quantity::Symbol,
        y_quantity::Symbol;
        <keyword arguments>
    )::Nothing

Plot two quantities as a scatter plot, one marker for every cell/particle.

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`. Each simulation will be plotted in a different figure.
  - `slice::IndexType`: Slice of the simulation, i.e. which snapshots will be plotted. It can be an integer (a single snapshot), a vector of integers (several snapshots), an `UnitRange` (e.g. 5:13), an `StepRange` (e.g. 5:2:13) or (:) (all snapshots). Starts at 1 and out of bounds indices are ignored.
  - `x_quantity::Symbol`: Quantity for the x axis. It can be any of the valid quantities of [`scatterQty`](@ref).
  - `y_quantity::Symbol`: Quantity for the y axis. It can be any of the valid quantities of [`scatterQty`](@ref).
  - `xlog::Bool=false`: If true, sets the x axis to ``\\log_{10}``(`x_quantity`).
  - `ylog::Bool=false`: If true, sets the y axis to ``\\log_{10}``(`y_quantity`).
  - `output_path::String="."`: Path to the output folder.
  - `trans_mode::Union{Symbol,Tuple{TranslationType,RotationType,Dict{Symbol,Vector{String}}}}=:all_box`: How to translate and rotate the cells/particles, before filtering with `filter_mode`. For options see [`selectTransformation`](@ref).
  - `filter_mode::Union{Symbol,Tuple{Function,Dict{Symbol,Vector{String}}}}=:all`: Which cells/particles will be selected. For options see [`selectFilter`](@ref).
  - `da_ff::Function=filterNothing`: Filter function to be applied within [`daScatterGalaxy`](@ref) after `trans_mode` and `filter_mode` are applied. See the required signature and examples in `./src/analysis/filters.jl`.
  - `ff_request::Dict{Symbol,Vector{String}}=Dict{Symbol,Vector{String}}()`: Request dictionary for `da_ff`.
  - `title::Union{Symbol,<:AbstractString}=""`: Title for the figure. If left empty, no title is printed. It can also be set to one of the following options:

      + `:physical_time` -> Physical time since the Big Bang.
      + `:lookback_time` -> Physical time left to reach the last snapshot.
      + `:scale_factor`  -> Scale factor (only relevant for cosmological simulations).
      + `:redshift`      -> Redshift (only relevant for cosmological simulations).
  - `theme::Attributes=Theme()`: Plot theme that will take precedence over [`DEFAULT_THEME`](@ref).
"""
function scatterPlot(
    simulation_paths::Vector{String},
    slice::IndexType,
    x_quantity::Symbol,
    y_quantity::Symbol;
    xlog::Bool=false,
    ylog::Bool=false,
    output_path::String=".",
    trans_mode::Union{Symbol,Tuple{TranslationType,RotationType,Dict{Symbol,Vector{String}}}}=:all_box,
    filter_mode::Union{Symbol,Tuple{Function,Dict{Symbol,Vector{String}}}}=:all,
    da_ff::Function=filterNothing,
    ff_request::Dict{Symbol,Vector{String}}=Dict{Symbol,Vector{String}}(),
    title::Union{Symbol,<:AbstractString}="",
    theme::Attributes=Theme(),
)::Nothing

    x_plot_params = plotParams(x_quantity)
    y_plot_params = plotParams(y_quantity)

    base_request = mergeRequests(x_plot_params.request, y_plot_params.request, ff_request)

    translation, rotation, trans_request = selectTransformation(trans_mode, base_request)
    filter_function, request = selectFilter(filter_mode, trans_request)

    # Set arguments for the x axis
    if xlog
        x_log        = x_plot_params.unit
        x_unit       = Unitful.NoUnits
        x_exp_factor = 0
        xaxis_label  = L"\log_{10} \, " * getLabel("auto_label", 0, x_plot_params.unit)
    else
        x_log        = nothing
        x_unit       = x_plot_params.unit
        x_exp_factor = x_plot_params.exp_factor
        xaxis_label  = x_plot_params.axis_label
    end

    # Set arguments for the y axis
    if ylog
        y_log        = y_plot_params.unit
        y_unit       = Unitful.NoUnits
        y_exp_factor = 0
        yaxis_label  = L"\log_{10} \, " * getLabel("auto_label", 0, y_plot_params.unit)
    else
        y_log        = nothing
        y_unit       = y_plot_params.unit
        y_exp_factor = y_plot_params.exp_factor
        yaxis_label  = y_plot_params.axis_label
    end

    for simulation_path in simulation_paths

        plotSnapshot(
            [simulation_path],
            request,
            [scatter!];
            pf_kwargs=[(; markersize=3)],
            output_path,
            base_filename="$(basename(simulation_path))_$(y_quantity)_vs_$(x_quantity)",
            slice,
            transform_box=true,
            translation,
            rotation,
            filter_function,
            da_functions=[daScatterGalaxy],
            da_args=[(x_quantity, y_quantity)],
            da_kwargs=[(; x_log, y_log, filter_function=da_ff)],
            x_unit,
            y_unit,
            x_exp_factor,
            y_exp_factor,
            xaxis_label,
            yaxis_label,
            xaxis_var_name=x_plot_params.var_name,
            yaxis_var_name=y_plot_params.var_name,
            theme,
            title,
        )

    end

    return nothing

end

"""
    scatterDensityMap(
        simulation_paths::Vector{String},
        slice::IndexType,
        x_quantity::Symbol,
        y_quantity::Symbol;
        <keyword arguments>
    )::Nothing

Plot two quantities as a density scatter plot (2D histogram).

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`. Each simulation will be plotted in a different figure.
  - `slice::IndexType`: Slice of the simulation, i.e. which snapshots will be plotted. It can be an integer (a single snapshot), a vector of integers (several snapshots), an `UnitRange` (e.g. 5:13), an `StepRange` (e.g. 5:2:13) or (:) (all snapshots). Starts at 1 and out of bounds indices are ignored.
  - `x_quantity::Symbol`: Quantity for the x axis. It can be any of the valid quantities of [`scatterQty`](@ref).
  - `y_quantity::Symbol`: Quantity for the y axis. It can be any of the valid quantities of [`scatterQty`](@ref).
  - `x_range::Union{NTuple{2,<:Number},Nothing}=nothing`: x axis range. If set to `nothing`, the extrema of the values will be used.
  - `y_range::Union{NTuple{2,<:Number},Nothing}=nothing`: y axis range. If set to `nothing`, the extrema of the values will be used.
  - `xlog::Bool=false`: If true, sets everything so the x axis is ``\\log_{10}``(`x_quantity`).
  - `ylog::Bool=false`: If true, sets everything so the y axis is ``\\log_{10}``(`y_quantity`).
  - `n_bins::Int=100`: Number of bins per side of the grid.
  - `output_path::String="."`: Path to the output folder.
  - `trans_mode::Union{Symbol,Tuple{TranslationType,RotationType,Dict{Symbol,Vector{String}}}}=:all_box`: How to translate and rotate the cells/particles, before filtering with `filter_mode`. For options see [`selectTransformation`](@ref).
  - `filter_mode::Union{Symbol,Tuple{Function,Dict{Symbol,Vector{String}}}}=:all`: Which cells/particles will be selected. For options see [`selectFilter`](@ref).
  - `da_ff::Function=filterNothing`: Filter function to be applied within [`daScatterGalaxy`](@ref) after `trans_mode` and `filter_mode` are applied. See the required signature and examples in `./src/analysis/filters.jl`.
  - `ff_request::Dict{Symbol,Vector{String}}=Dict{Symbol,Vector{String}}()`: Request dictionary for `da_ff`.
  - `colorbar::Bool=false`: If a colorbar will be added.
  - `title::Union{Symbol,<:AbstractString}=""`: Title for the figure. If left empty, no title is printed. It can also be set to one of the following options:

      + `:physical_time` -> Physical time since the Big Bang.
      + `:lookback_time` -> Physical time left to reach the last snapshot.
      + `:scale_factor`  -> Scale factor (only relevant for cosmological simulations).
      + `:redshift`      -> Redshift (only relevant for cosmological simulations).
  - `theme::Attributes=Theme()`: Plot theme that will take precedence over [`DEFAULT_THEME`](@ref).
"""
function scatterDensityMap(
    simulation_paths::Vector{String},
    slice::IndexType,
    x_quantity::Symbol,
    y_quantity::Symbol;
    x_range::Union{NTuple{2,<:Number},Nothing}=nothing,
    y_range::Union{NTuple{2,<:Number},Nothing}=nothing,
    xlog::Bool=false,
    ylog::Bool=false,
    n_bins::Int=100,
    output_path::String=".",
    trans_mode::Union{Symbol,Tuple{TranslationType,RotationType,Dict{Symbol,Vector{String}}}}=:all_box,
    filter_mode::Union{Symbol,Tuple{Function,Dict{Symbol,Vector{String}}}}=:all,
    da_ff::Function=filterNothing,
    ff_request::Dict{Symbol,Vector{String}}=Dict{Symbol,Vector{String}}(),
    colorbar::Bool=false,
    title::Union{Symbol,<:AbstractString}="",
    theme::Attributes=Theme(),
)::Nothing

    x_plot_params = plotParams(x_quantity)
    y_plot_params = plotParams(y_quantity)

    base_request = mergeRequests(x_plot_params.request, y_plot_params.request, ff_request)

    translation, rotation, trans_request = selectTransformation(trans_mode, base_request)
    filter_function, request = selectFilter(filter_mode, trans_request)

    # Set arguments for the x axis
    if xlog
        x_log        = x_plot_params.unit
        x_unit       = Unitful.NoUnits
        x_exp_factor = 0
        xaxis_label  = L"\log_{10} \, " * getLabel("auto_label", 0, x_plot_params.unit)
    else
        x_log        = nothing
        x_unit       = x_plot_params.unit
        x_exp_factor = x_plot_params.exp_factor
        xaxis_label  = x_plot_params.axis_label
    end

    # Set arguments for the y axis
    if ylog
        y_log        = y_plot_params.unit
        y_unit       = Unitful.NoUnits
        y_exp_factor = 0
        yaxis_label  = L"\log_{10} \, " * getLabel("auto_label", 0, y_plot_params.unit)
    else
        y_log        = nothing
        y_unit       = y_plot_params.unit
        y_exp_factor = y_plot_params.exp_factor
        yaxis_label  = y_plot_params.axis_label
    end

    for simulation_path in simulation_paths

        plotSnapshot(
            [simulation_path],
            request,
            [heatmap!];
            output_path,
            base_filename="$(basename(simulation_path))_$(y_quantity)_vs_$(x_quantity)",
            slice,
            transform_box=true,
            translation,
            rotation,
            filter_function,
            da_functions=[daScatterDensity],
            da_args=[(x_quantity, y_quantity)],
            da_kwargs=[(; x_range, y_range, x_log, y_log, n_bins, filter_function=da_ff)],
            x_unit,
            y_unit,
            x_exp_factor,
            y_exp_factor,
            xaxis_label,
            yaxis_label,
            xaxis_var_name=x_plot_params.var_name,
            yaxis_var_name=y_plot_params.var_name,
            theme=merge(theme, Theme(size=(880, 720), Colorbar=(label=L"\log_{10} \, #",))),
            title,
            colorbar,
        )

    end

    return nothing

end

"""
    scatterDensityMap(
        simulation_paths::Vector{String},
        slice::IndexType,
        x_quantity::Symbol,
        y_quantity::Symbol,
        z_quantity::Symbol;
        <keyword arguments>
    )::Nothing

Plot two quantities as a density scatter plot (2D histogram), weighted by `z_quantity`.

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`. Each simulation will be plotted in a different figure.
  - `slice::IndexType`: Slice of the simulation, i.e. which snapshots will be plotted. It can be an integer (a single snapshot), a vector of integers (several snapshots), an `UnitRange` (e.g. 5:13), an `StepRange` (e.g. 5:2:13) or (:) (all snapshots). Starts at 1 and out of bounds indices are ignored.
  - `x_quantity::Symbol`: Quantity for the x axis. It can be any of the valid quantities of [`scatterQty`](@ref).
  - `y_quantity::Symbol`: Quantity for the y axis. It can be any of the valid quantities of [`scatterQty`](@ref).
  - `z_quantity::Symbol`: Quantity for the weights. It can be any of the valid quantities of [`scatterQty`](@ref).
  - `x_range::Union{NTuple{2,<:Number},Nothing}=nothing`: x axis range. If set to `nothing`, the extrema of the values will be used.
  - `y_range::Union{NTuple{2,<:Number},Nothing}=nothing`: y axis range. If set to `nothing`, the extrema of the values will be used.
  - `xlog::Bool=false`: If true, sets everything so the x axis is ``\\log_{10}``(`x_quantity`).
  - `ylog::Bool=false`: If true, sets everything so the y axis is ``\\log_{10}``(`y_quantity`).
  - `total::Bool=true`: If the sum (default) or the mean of `z_quantity` will be used as the value of each bin.
  - `n_bins::Int=100`: Number of bins per side of the grid.
  - `output_path::String="."`: Path to the output folder.
  - `trans_mode::Union{Symbol,Tuple{TranslationType,RotationType,Dict{Symbol,Vector{String}}}}=:all_box`: How to translate and rotate the cells/particles, before filtering with `filter_mode`. For options see [`selectTransformation`](@ref).
  - `filter_mode::Union{Symbol,Tuple{Function,Dict{Symbol,Vector{String}}}}=:all`: Which cells/particles will be selected. For options see [`selectFilter`](@ref).
  - `da_ff::Function=filterNothing`: Filter function to be applied within [`daScatterGalaxy`](@ref) after `trans_mode` and `filter_mode` are applied. See the required signature and examples in `./src/analysis/filters.jl`.
  - `ff_request::Dict{Symbol,Vector{String}}=Dict{Symbol,Vector{String}}()`: Request dictionary for `da_ff`.
  - `colorbar::Bool=false`: If a colorbar will be added.
  - `title::Union{Symbol,<:AbstractString}=""`: Title for the figure. If left empty, no title is printed. It can also be set to one of the following options:

      + `:physical_time` -> Physical time since the Big Bang.
      + `:lookback_time` -> Physical time left to reach the last snapshot.
      + `:scale_factor`  -> Scale factor (only relevant for cosmological simulations).
      + `:redshift`      -> Redshift (only relevant for cosmological simulations).
  - `theme::Attributes=Theme()`: Plot theme that will take precedence over [`DEFAULT_THEME`](@ref).
"""
function scatterDensityMap(
    simulation_paths::Vector{String},
    slice::IndexType,
    x_quantity::Symbol,
    y_quantity::Symbol,
    z_quantity::Symbol;
    x_range::Union{NTuple{2,<:Number},Nothing}=nothing,
    y_range::Union{NTuple{2,<:Number},Nothing}=nothing,
    xlog::Bool=false,
    ylog::Bool=false,
    total::Bool=true,
    n_bins::Int=100,
    output_path::String=".",
    trans_mode::Union{Symbol,Tuple{TranslationType,RotationType,Dict{Symbol,Vector{String}}}}=:all_box,
    filter_mode::Union{Symbol,Tuple{Function,Dict{Symbol,Vector{String}}}}=:all,
    da_ff::Function=filterNothing,
    ff_request::Dict{Symbol,Vector{String}}=Dict{Symbol,Vector{String}}(),
    colorbar::Bool=false,
    title::Union{Symbol,<:AbstractString}="",
    theme::Attributes=Theme(),
)::Nothing

    x_plot_params = plotParams(x_quantity)
    y_plot_params = plotParams(y_quantity)
    z_plot_params = plotParams(z_quantity)

    base_request = mergeRequests(
        x_plot_params.request,
        y_plot_params.request,
        z_plot_params.request,
        ff_request,
    )

    translation, rotation, trans_request = selectTransformation(trans_mode, base_request)
    filter_function, request = selectFilter(filter_mode, trans_request)

    # Set arguments for the x axis
    if xlog
        x_log        = x_plot_params.unit
        x_unit       = Unitful.NoUnits
        x_exp_factor = 0
        xaxis_label  = L"\log_{10} \, " * getLabel("auto_label", 0, x_plot_params.unit)
    else
        x_log        = nothing
        x_unit       = x_plot_params.unit
        x_exp_factor = x_plot_params.exp_factor
        xaxis_label  = x_plot_params.axis_label
    end

    # Set arguments for the y axis
    if ylog
        y_log        = y_plot_params.unit
        y_unit       = Unitful.NoUnits
        y_exp_factor = 0
        yaxis_label  = L"\log_{10} \, " * getLabel("auto_label", 0, y_plot_params.unit)
    else
        y_log        = nothing
        y_unit       = y_plot_params.unit
        y_exp_factor = y_plot_params.exp_factor
        yaxis_label  = y_plot_params.axis_label
    end

    # Label for the colorbar
    colorbar_label = LaTeXString(
        L"\log_{10} \, " * getLabel(z_plot_params.var_name, 0, z_plot_params.unit)
    )

    for simulation_path in simulation_paths

        plotSnapshot(
            [simulation_path],
            request,
            [heatmap!];
            output_path,
            base_filename="$(basename(simulation_path))_$(y_quantity)_vs_$(x_quantity)",
            slice,
            transform_box=true,
            translation,
            rotation,
            filter_function,
            da_functions=[daScatterWeightedDensity],
            da_args=[(x_quantity, y_quantity, z_quantity)],
            da_kwargs=[(; x_range, y_range, x_log, y_log, total, n_bins, filter_function=da_ff)],
            x_unit,
            y_unit,
            x_exp_factor,
            y_exp_factor,
            xaxis_label,
            yaxis_label,
            xaxis_var_name=x_plot_params.var_name,
            yaxis_var_name=y_plot_params.var_name,
            theme=merge(theme, Theme(size=(880, 720), Colorbar=(label=colorbar_label,))),
            title,
            colorbar,
        )

    end

    return nothing

end

"""
    radialProfile(
        simulation_paths::Vector{String},
        slice::IndexType,
        quantity::Symbol;
        <keyword arguments>
    )::Nothing

Plot a radial profile.

!!! note

    This method plots one quantity for several simulations in one figure.

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`. All the simulations will be plotted together.
  - `slice::IndexType`: Slice of the simulation, i.e. which snapshots will be plotted. It can be an integer (a single snapshot), a vector of integers (several snapshots), an `UnitRange` (e.g. 5:13), an `StepRange` (e.g. 5:2:13) or (:) (all snapshots). It works over the longest simulation. Starts at 1 and out of bounds indices are ignored.
  - `quantity::Symbol`: Target quantity. It can be any of the valid quantities of [`scatterQty`](@ref).
  - `norm::Union{Symbol,Nothing}=nothing`: The value of `quantity` in each bin will be divided by the corresponding value of `norm`. It can be any of the valid quantities of [`scatterQty`](@ref). If set to `nothing`, no operation is applied.
  - `radius::Unitful.Length=DISK_R`: Radius of the profile.
  - `n_bins::Int=100`: Number of bins.
  - `ylog::Bool=false`: If true, returns the profile of ``\\log_{10}``(`quantity`) (after dividing by `norm`).
  - `flat::Bool=true`: If the profile will be 2D (rings), or 3D (spherical shells).
  - `total::Bool=true`: If the sum (default) or the mean of `quantity` will be computed for each bin. This affects the values of `norm` too.
  - `cumulative::Bool=false`: If the profile will be accumulated (after dividing by `norm`).
  - `density::Bool=true`: If the profile will be of the density of `quantity` (after dividing by `norm`).
  - `output_path::String="."`: Path to the output folder.
  - `trans_mode::Union{Symbol,Tuple{TranslationType,RotationType,Dict{Symbol,Vector{String}}}}=:all_box`: How to translate and rotate the cells/particles, before filtering with `filter_mode`. For options see [`selectTransformation`](@ref).
  - `filter_mode::Union{Symbol,Tuple{Function,Dict{Symbol,Vector{String}}}}=:all`: Which cells/particles will be selected. For options see [`selectFilter`](@ref).
  - `da_ff::Function=filterNothing`: Filter function to be applied within [`daScatterGalaxy`](@ref) after `trans_mode` and `filter_mode` are applied. See the required signature and examples in `./src/analysis/filters.jl`.
  - `ff_request::Dict{Symbol,Vector{String}}=Dict{Symbol,Vector{String}}()`: Request dictionary for `da_ff`.
  - `sim_labels::Union{Vector{<:AbstractString},Nothing}=basename.(simulation_paths)`: Labels for the plot legend, one per simulation. Set it to `nothing` if you don't want a legend.
  - `title::Union{Symbol,<:AbstractString}=""`: Title for the figure. If left empty, no title is printed. It can also be set to one of the following options:

      + `:physical_time` -> Physical time since the Big Bang.
      + `:lookback_time` -> Physical time left to reach the last snapshot.
      + `:scale_factor`  -> Scale factor (only relevant for cosmological simulations).
      + `:redshift`      -> Redshift (only relevant for cosmological simulations).
  - `theme::Attributes=Theme()`: Plot theme that will take precedence over [`DEFAULT_THEME`](@ref).
"""
function radialProfile(
    simulation_paths::Vector{String},
    slice::IndexType,
    quantity::Symbol;
    norm::Union{Symbol,Nothing}=nothing,
    radius::Unitful.Length=DISK_R,
    n_bins::Int=100,
    ylog::Bool=false,
    flat::Bool=true,
    total::Bool=true,
    cumulative::Bool=false,
    density::Bool=true,
    output_path::String=".",
    trans_mode::Union{Symbol,Tuple{TranslationType,RotationType,Dict{Symbol,Vector{String}}}}=:all_box,
    filter_mode::Union{Symbol,Tuple{Function,Dict{Symbol,Vector{String}}}}=:all,
    da_ff::Function=filterNothing,
    ff_request::Dict{Symbol,Vector{String}}=Dict{Symbol,Vector{String}}(),
    sim_labels::Union{Vector{<:AbstractString},Nothing}=basename.(simulation_paths),
    title::Union{Symbol,<:AbstractString}="",
    theme::Attributes=Theme(),
)::Nothing

    plot_params = plotParams(quantity)

    # Compute the unit and request of the norm
    if isnothing(norm)
        var_name     = plot_params.var_name
        u_norm       = Unitful.NoUnits
        norm_request = Dict{Symbol,Vector{String}}()
    else
        n_plot_params = plotParams(norm)
        var_name      = LaTeXString(plot_params.var_name * L"\, / \," * n_plot_params.var_name)
        u_norm        = n_plot_params.unit
        norm_request  = n_plot_params.request
    end

    # Compute the unit of the bins
    if density
        u_bin = flat ? u"kpc"^2 : u"kpc"^3
    else
        u_bin = Unitful.NoUnits
    end

    # Compute the unit of the y axis
    yunit = plot_params.unit / u_norm / u_bin

    # Compute the label of the y axis
    ylabel = getLabel(var_name, 0, yunit)

    if ylog
        yaxis_label  = L"\log_{10} \, " * ylabel
        y_log        = yunit
    else
        yaxis_label  = ylabel
        y_log        = nothing
    end

    base_request = mergeRequests(plot_params.request, norm_request, ff_request)

    translation, rotation, trans_request = selectTransformation(trans_mode, base_request)
    filter_function, request = selectFilter(filter_mode, trans_request)

    grid = CircularGrid(radius, n_bins)

    if isone(length(simulation_paths))
        base_filename = "$(basename(simulation_paths[1]))_$(quantity)_radial_profile"
    else
        base_filename = "$(quantity)_radial_profile"
    end

    plotSnapshot(
        simulation_paths,
        request,
        [lines!];
        output_path,
        base_filename,
        slice,
        transform_box=true,
        translation,
        rotation,
        filter_function,
        da_functions=[daProfile],
        da_args=[(quantity, grid)],
        da_kwargs=[(; norm, y_log, flat, total, cumulative, density, filter_function=da_ff)],
        x_unit=u"kpc",
        y_unit = ylog ? Unitful.NoUnits : yunit,
        yaxis_label,
        xaxis_var_name=L"r",
        theme=merge(
            theme,
            Theme(
                size=(1200, 880),
                Axis=(aspect=nothing,),
                Legend=(valign=:top, nbanks=1, margin=(0, 2, 0, 2)),
            ),
        ),
        sim_labels,
        title,
    )

    return nothing

end

"""
    radialProfile(
        simulation_paths::Vector{String},
        slice::IndexType,
        quantities::Vector{Symbol},
        ylabel::AbstractString;
        <keyword arguments>
    )::Nothing

Plot a density profile.

!!! note

    This method plots several quantities for one simulations in one figure.

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`. Each simulation will be plotted in a different figure.
  - `slice::IndexType`: Slice of the simulation, i.e. which snapshots will be plotted. It can be an integer (a single snapshot), a vector of integers (several snapshots), an `UnitRange` (e.g. 5:13), an `StepRange` (e.g. 5:2:13) or (:) (all snapshots). Starts at 1 and out of bounds indices are ignored.
  - `quantities::Vector{Symbol}`: Target quantities. They can be any of the valid quantities of [`scatterQty`](@ref).
  - `q_var_name::AbstractString`: Name of the variable for the y axis.
  - `q_unit::Unitful.Units=Unitful.NoUnits`: Unit of `quantities`. All must have the same units.
  - `norm::Union{Symbol,Nothing}=nothing`: The value of `quantity` in each bin will be divided by the corresponding value of `norm`. It can be any of the valid quantities of [`scatterQty`](@ref). If set to `nothing`, no operation is applied.
  - `radius::Unitful.Length=DISK_R`: Radius of the profile.
  - `n_bins::Int=100`: Number of bins.
  - `ylog::Bool=false`: If true, returns the profile of ``\\log_{10}``(`quantity`) (after dividing by `norm`).
  - `flat::Bool=true`: If the profile will be 2D (rings), or 3D (spherical shells).
  - `total::Bool=true`: If the sum (default) or the mean of `quantity` will be computed for each bin. This affects the values of `norm` too.
  - `cumulative::Bool=false`: If the profile will be accumulated (after dividing by `norm`).
  - `density::Bool=true`: If the profile will be of the density of `quantity` (after dividing by `norm`).
  - `output_path::String="."`: Path to the output folder.
  - `trans_mode::Union{Symbol,Tuple{TranslationType,RotationType,Dict{Symbol,Vector{String}}}}=:all_box`: How to translate and rotate the cells/particles, before filtering with `filter_mode`. For options see [`selectTransformation`](@ref).
  - `filter_mode::Union{Symbol,Tuple{Function,Dict{Symbol,Vector{String}}}}=:all`: Which cells/particles will be selected. For options see [`selectFilter`](@ref).
  - `da_ff::Function=filterNothing`: Filter function to be applied within [`daScatterGalaxy`](@ref) after `trans_mode` and `filter_mode` are applied. See the required signature and examples in `./src/analysis/filters.jl`.
  - `ff_request::Dict{Symbol,Vector{String}}=Dict{Symbol,Vector{String}}()`: Request dictionary for `da_ff`.
  - `title::Union{Symbol,<:AbstractString}=""`: Title for the figure. If left empty, no title is printed. It can also be set to one of the following options:

      + `:physical_time` -> Physical time since the Big Bang.
      + `:lookback_time` -> Physical time left to reach the last snapshot.
      + `:scale_factor`  -> Scale factor (only relevant for cosmological simulations).
      + `:redshift`      -> Redshift (only relevant for cosmological simulations).
  - `theme::Attributes=Theme()`: Plot theme that will take precedence over [`DEFAULT_THEME`](@ref).
"""
function radialProfile(
    simulation_paths::Vector{String},
    slice::IndexType,
    quantities::Vector{Symbol},
    q_var_name::AbstractString;
    q_unit::Unitful.Units=Unitful.NoUnits,
    norm::Union{Symbol,Nothing}=nothing,
    radius::Unitful.Length=DISK_R,
    n_bins::Int=100,
    ylog::Bool=false,
    flat::Bool=true,
    total::Bool=true,
    cumulative::Bool=false,
    density::Bool=true,
    output_path::String=".",
    trans_mode::Union{Symbol,Tuple{TranslationType,RotationType,Dict{Symbol,Vector{String}}}}=:all_box,
    filter_mode::Union{Symbol,Tuple{Function,Dict{Symbol,Vector{String}}}}=:all,
    da_ff::Function=filterNothing,
    ff_request::Dict{Symbol,Vector{String}}=Dict{Symbol,Vector{String}}(),
    title::Union{Symbol,<:AbstractString}="",
    theme::Attributes=Theme(),
)::Nothing

    requests = [plotParams(quantity).request for quantity in quantities]

    # Compute the unit and request of the norm
    if isnothing(norm)
        var_name     = q_var_name
        u_norm       = Unitful.NoUnits
        norm_request = Dict{Symbol,Vector{String}}()
    else
        n_plot_params = plotParams(norm)
        var_name      = LaTeXString(q_var_name * L"\, / \," * n_plot_params.var_name)
        u_norm        = n_plot_params.unit
        norm_request  = n_plot_params.request
    end

    # Compute the unit of the bins
    if density
        u_bin = flat ? u"kpc"^2 : u"kpc"^3
    else
        u_bin = Unitful.NoUnits
    end

    # Compute the unit of the y axis
    yunit = q_unit / u_norm / u_bin

    # Compute the label of the y axis
    ylabel = getLabel(var_name, 0, yunit)

    if ylog
        yaxis_label  = L"\log_{10} \, " * ylabel
        y_log        = yunit
    else
        yaxis_label  = ylabel
        y_log        = nothing
    end

    base_request = mergeRequests(norm_request, ff_request, requests...)

    translation, rotation, trans_request = selectTransformation(trans_mode, base_request)
    filter_function, request = selectFilter(filter_mode, trans_request)

    grid = CircularGrid(radius, n_bins)

    for simulation_path in simulation_paths

        plotSnapshot(
            fill(simulation_path, length(quantities)),
            request,
            [lines!];
            output_path,
            base_filename="$(basename(simulation_path))_radial_profile",
            slice,
            transform_box=true,
            translation,
            rotation,
            filter_function,
            da_functions=[daProfile],
            da_args=[(quantity, grid) for quantity in quantities],
            da_kwargs=[(; norm, y_log, flat, total, cumulative, density, filter_function=da_ff)],
            x_unit=u"kpc",
            y_unit = ylog ? Unitful.NoUnits : yunit,
            yaxis_label,
            xaxis_var_name=L"r",
            theme=merge(
                theme,
                Theme(
                    size=(1200, 880),
                    Axis=(aspect=nothing,),
                    Legend=(valign=:top, nbanks=1, margin=(0, 2, 0, 2)),
                ),
            ),
            sim_labels=string.(quantities),
            title,
        )

    end

    return nothing

end

"""
    stellarHistory(
        simulation_paths::Vector{String},
        slice::IndexType,
        quantity::Symbol;
        <keyword arguments>
    )::Nothing

Plot the evolution of a given stellar `quantity` using the stellar ages at a given instant in time.

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`. All the simulations will be plotted together.
  - `slice::IndexType`: Slice of the simulation, i.e. which snapshots will be plotted. It can be an integer (a single snapshot), a vector of integers (several snapshots), an `UnitRange` (e.g. 5:13), an `StepRange` (e.g. 5:2:13) or (:) (all snapshots). It works over the longest simulation. Starts at 1 and out of bounds indices are ignored.
  - `quantity::Symbol=:sfr`: Target quantity. The options are:

      + `:sfr`                 -> Star formation rate.
      + `:ssfr`                -> The specific star formation rate.
      + `:stellar_mass`        -> Stellar mass.
      + `:stellar_metallicity` -> Mass fraction of all elements above He in the stars (solar units).
  - `n_bins::Int=100`: Number of bins.
  - `ylog::Bool=false`: If true, returns the profile of ``\\log_{10}``(`quantity`) (after dividing by `norm`).
  - `output_path::String="."`: Path to the output folder.
  - `trans_mode::Union{Symbol,Tuple{TranslationType,RotationType,Dict{Symbol,Vector{String}}}}=:all_box`: How to translate and rotate the cells/particles, before filtering with `filter_mode`. For options see [`selectTransformation`](@ref).
  - `filter_mode::Union{Symbol,Tuple{Function,Dict{Symbol,Vector{String}}}}=:all`: Which cells/particles will be selected. For options see [`selectFilter`](@ref).
  - `da_ff::Function=filterNothing`: Filter function to be applied within [`daScatterGalaxy`](@ref) after `trans_mode` and `filter_mode` are applied. See the required signature and examples in `./src/analysis/filters.jl`.
  - `ff_request::Dict{Symbol,Vector{String}}=Dict{Symbol,Vector{String}}()`: Request dictionary for `da_ff`.
  - `sim_labels::Union{Vector{<:AbstractString},Nothing}=basename.(simulation_paths)`: Labels for the plot legend, one per simulation. Set it to `nothing` if you don't want a legend.
  - `theme::Attributes=Theme()`: Plot theme that will take precedence over [`DEFAULT_THEME`](@ref).
"""
function stellarHistory(
    simulation_paths::Vector{String},
    slice::IndexType;
    quantity::Symbol=:sfr,
    n_bins::Int=100,
    ylog::Bool=false,
    output_path::String=".",
    trans_mode::Union{Symbol,Tuple{TranslationType,RotationType,Dict{Symbol,Vector{String}}}}=:all_box,
    filter_mode::Union{Symbol,Tuple{Function,Dict{Symbol,Vector{String}}}}=:all,
    da_ff::Function=filterNothing,
    ff_request::Dict{Symbol,Vector{String}}=Dict{Symbol,Vector{String}}(),
    sim_labels::Union{Vector{<:AbstractString},Nothing}=basename.(simulation_paths),
    theme::Attributes=Theme(),
)::Nothing

    x_plot_params = plotParams(:physical_time)
    y_plot_params = plotParams(quantity)

    base_request = mergeRequests(
        Dict(:stellar => ["GAGE", "MASS"]),
        y_plot_params.request,
        ff_request,
    )

    translation, rotation, trans_request = selectTransformation(trans_mode, base_request)
    filter_function, request = selectFilter(filter_mode, trans_request)

    # Set arguments for the y axis
    if ylog
        y_unit       = Unitful.NoUnits
        y_exp_factor = 0
        yaxis_label  = L"\log_{10} \, " * getLabel("auto_label", 0, y_plot_params.unit)
    else
        y_unit       = y_plot_params.unit
        y_exp_factor = y_plot_params.exp_factor
        yaxis_label  = y_plot_params.axis_label
    end

    if isone(length(simulation_paths))
        base_filename = "$(basename(simulation_paths[1]))_$(quantity)_history"
    else
        base_filename = "$(quantity)_history"
    end

    plotSnapshot(
        simulation_paths,
        request,
        [lines!];
        output_path,
        base_filename,
        slice,
        transform_box=true,
        translation,
        rotation,
        filter_function,
        da_functions=[daStellarHistory],
        da_kwargs=[(; quantity, n_bins, ylog, filter_function=da_ff)],
        x_unit=x_plot_params.unit,
        y_unit,
        x_exp_factor=x_plot_params.exp_factor,
        y_exp_factor,
        xaxis_label=x_plot_params.axis_label,
        yaxis_label,
        xaxis_var_name=x_plot_params.var_name,
        yaxis_var_name=y_plot_params.var_name,
        theme=merge(
            theme,
            Theme(
                size=(1200, 880),
                Axis=(aspect=nothing,),
                Legend=(valign=:top, nbanks=1, margin=(0, 2, 0, 2)),
            ),
        ),
        sim_labels,
    )

    return nothing

end

"""
    histogram(
        simulation_paths::Vector{String},
        slice::IndexType,
        quantity::Symbol,
        range::NTuple{2,<:Number};
        <keyword arguments>
    )::Nothing

Plot a histogram of `quantity`.

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`. All the simulations will be plotted together.
  - `slice::IndexType`: Slice of the simulation, i.e. which snapshots will be plotted. It can be an integer (a single snapshot), a vector of integers (several snapshots), an `UnitRange` (e.g. 5:13), an `StepRange` (e.g. 5:2:13) or (:) (all snapshots). It works over the longest simulation. Starts at 1 and out of bounds indices are ignored.
  - `quantity::Symbol`: Target quantity. It can be any of the valid quantities of [`scatterQty`](@ref).
  - `n_bins::Int=100`: Number of bins.
  - `line::Bool=true`: If the histogram will be plotted with a line or with bars.
  - `norm::Int=0`: Number of count that will be use to normalize the histogram. If left as 0, the histogram will be normalize with the maximum bin count.
  - `range::Union{NTuple{2,<:Number},Nothing}=nothing,`: Range of values for the histogram. If set to `nothing`, the extrema of the values will be used.
  - `xlog::Bool=false`: If the histogram bins will be logarithmic.
  - `output_path::String="."`: Path to the output folder.
  - `trans_mode::Union{Symbol,Tuple{TranslationType,RotationType,Dict{Symbol,Vector{String}}}}=:all_box`: How to translate and rotate the cells/particles, before filtering with `filter_mode`. For options see [`selectTransformation`](@ref).
  - `filter_mode::Union{Symbol,Tuple{Function,Dict{Symbol,Vector{String}}}}=:all`: Which cells/particles will be selected. For options see [`selectFilter`](@ref).
  - `da_ff::Function=filterNothing`: Filter function to be applied within [`daScatterGalaxy`](@ref) after `trans_mode` and `filter_mode` are applied. See the required signature and examples in `./src/analysis/filters.jl`.
  - `ff_request::Dict{Symbol,Vector{String}}=Dict{Symbol,Vector{String}}()`: Request dictionary for `da_ff`.
  - `sim_labels::Union{Vector{<:AbstractString},Nothing}=basename.(simulation_paths)`: Labels for the plot legend, one per simulation. Set it to `nothing` if you don't want a legend.ç
  - `title::Union{Symbol,<:AbstractString}=""`: Title for the figure. If left empty, no title is printed. It can also be set to one of the following options:

      + `:physical_time` -> Physical time since the Big Bang.
      + `:lookback_time` -> Physical time left to reach the last snapshot.
      + `:scale_factor`  -> Scale factor (only relevant for cosmological simulations).
      + `:redshift`      -> Redshift (only relevant for cosmological simulations).
  - `theme::Attributes=Theme()`: Plot theme that will take precedence over [`DEFAULT_THEME`](@ref).
"""
function histogram(
    simulation_paths::Vector{String},
    slice::IndexType,
    quantity::Symbol;
    n_bins::Int=100,
    line::Bool=true,
    norm::Int=0,
    range::Union{NTuple{2,<:Number},Nothing}=nothing,
    xlog::Bool=false,
    output_path::String=".",
    trans_mode::Union{Symbol,Tuple{TranslationType,RotationType,Dict{Symbol,Vector{String}}}}=:all_box,
    filter_mode::Union{Symbol,Tuple{Function,Dict{Symbol,Vector{String}}}}=:all,
    da_ff::Function=filterNothing,
    ff_request::Dict{Symbol,Vector{String}}=Dict{Symbol,Vector{String}}(),
    sim_labels::Union{Vector{<:AbstractString},Nothing}=basename.(simulation_paths),
    title::Union{Symbol,<:AbstractString}="",
    theme::Attributes=Theme(),
)::Nothing

    plot_params = plotParams(quantity)

    base_request = mergeRequests(plot_params.request, ff_request)

    translation, rotation, trans_request = selectTransformation(trans_mode, base_request)
    filter_function, request = selectFilter(filter_mode, trans_request)

    if isnothing(range)
        da_args = [(quantity, n_bins, xlog)]
    else
        grid    = LinearGrid(range..., n_bins; log=xlog)
        da_args = [(quantity, grid)]
    end

    # Set arguments for the x axis
    if xlog
        x_unit       = Unitful.NoUnits
        x_exp_factor = 0
        xaxis_label  = L"\log_{10} \, " * getLabel("auto_label", 0, plot_params.unit)
    else
        x_unit       = plot_params.unit
        x_exp_factor = plot_params.exp_factor
        xaxis_label  = plot_params.axis_label
    end

    if isPositive(norm)
        yaxis_var_name = L"\mathrm{Counts \, / \, %$(norm)}"
    else
        yaxis_var_name = L"\mathrm{Counts \, / \, \max(Counts)}"
    end

    if isone(length(simulation_paths))
        base_filename = "$(basename(simulation_paths[1]))_$(quantity)_histogram"
    else
        base_filename = "$(quantity)_histogram"
    end

    plotSnapshot(
        simulation_paths,
        request,
        line ? [lines!] : [barplot!];
        output_path,
        base_filename,
        slice,
        transform_box=true,
        translation,
        rotation,
        filter_function,
        da_functions=[daHistogram],
        da_args,
        da_kwargs=[(; filter_function=da_ff, norm)],
        x_unit,
        x_exp_factor,
        xaxis_label,
        xaxis_var_name=plot_params.var_name,
        yaxis_var_name,
        theme=merge(theme, Theme(Legend=(nbanks=1, valign=:top, margin=(0, 10, 0, 0)),)),
        sim_labels,
        title,
    )

    return nothing

end

"""
    rotationCurve(
        simulation_paths::Vector{String},
        slice::IndexType;
        <keyword arguments>
    )::Nothing

Plot the rotation curve of several simulations.

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`. All the simulations will be plotted together.
  - `slice::IndexType`: Slice of the simulation, i.e. which snapshots will be plotted. It can be an integer (a single snapshot), a vector of integers (several snapshots), an `UnitRange` (e.g. 5:13), an `StepRange` (e.g. 5:2:13) or (:) (all snapshots). It works over the longest simulation. Starts at 1 and out of bounds indices are ignored.
  - `R::Unitful.Length=DISK_R`: Maximum radial distance for the rotation curve.
  - `output_path::String="."`: Path to the output folder.
  - `trans_mode::Union{Symbol,Tuple{TranslationType,RotationType,Dict{Symbol,Vector{String}}}}=:all_box`: How to translate and rotate the cells/particles, before filtering with `filter_mode`. For options see [`selectTransformation`](@ref).
  - `filter_mode::Union{Symbol,Tuple{Function,Dict{Symbol,Vector{String}}}}=:all`: Which cells/particles will be selected. For options see [`selectFilter`](@ref).
  - `da_ff::Function=filterNothing`: Filter function to be applied within [`daScatterGalaxy`](@ref) after `trans_mode` and `filter_mode` are applied. See the required signature and examples in `./src/analysis/filters.jl`.
  - `ff_request::Dict{Symbol,Vector{String}}=Dict{Symbol,Vector{String}}()`: Request dictionary for `da_ff`.
  - `sim_labels::Union{Vector{<:AbstractString},Nothing}=basename.(simulation_paths)`: Labels for the plot legend, one per simulation. Set it to `nothing` if you don't want a legend.
  - `title::Union{Symbol,<:AbstractString}=""`: Title for the figure. If left empty, no title is printed. It can also be set to one of the following options:

      + `:physical_time` -> Physical time since the Big Bang.
      + `:lookback_time` -> Physical time left to reach the last snapshot.
      + `:scale_factor`  -> Scale factor (only relevant for cosmological simulations).
      + `:redshift`      -> Redshift (only relevant for cosmological simulations).
  - `theme::Attributes=Theme()`: Plot theme that will take precedence over [`DEFAULT_THEME`](@ref).
"""
function rotationCurve(
    simulation_paths::Vector{String},
    slice::IndexType;
    R::Unitful.Length=DISK_R,
    output_path::String=".",
    trans_mode::Union{Symbol,Tuple{TranslationType,RotationType,Dict{Symbol,Vector{String}}}}=:all_box,
    filter_mode::Union{Symbol,Tuple{Function,Dict{Symbol,Vector{String}}}}=:all,
    da_ff::Function=filterNothing,
    ff_request::Dict{Symbol,Vector{String}}=Dict{Symbol,Vector{String}}(),
    sim_labels::Union{Vector{<:AbstractString},Nothing}=basename.(simulation_paths),
    title::Union{Symbol,<:AbstractString}="",
    theme::Attributes=Theme(),
)::Nothing

    x_plot_params = plotParams(:stellar_radial_distance)
    y_plot_params = plotParams(:stellar_circular_velocity)

    base_request = mergeRequests(x_plot_params.request, y_plot_params.request, ff_request)

    translation, rotation, trans_request = selectTransformation(trans_mode, base_request)
    filter_function, request = selectFilter(filter_mode, trans_request)

    # Smooth out the curve with 5 bins per kpc
    smooth = round(Int64, 5 * ustrip(u"kpc", R))

    if isone(length(simulation_paths))
        base_filename = "$(basename(simulation_paths[1]))_rotation_curve"
    else
        base_filename = "rotation_curve"
    end

    plotSnapshot(
        simulation_paths,
        request,
        [lines!];
        output_path,
        base_filename,
        slice,
        transform_box=true,
        translation,
        rotation,
        filter_function,
        da_functions=[daRotationCurve],
        da_kwargs=[(; R, filter_function=da_ff)],
        smooth,
        x_unit=x_plot_params.unit,
        y_unit=y_plot_params.unit,
        x_exp_factor=x_plot_params.exp_factor,
        y_exp_factor=y_plot_params.exp_factor,
        x_scale_func=identity,
        y_scale_func=identity,
        xaxis_label=x_plot_params.axis_label,
        yaxis_label=y_plot_params.axis_label,
        xaxis_var_name=x_plot_params.var_name,
        yaxis_var_name=y_plot_params.var_name,
        theme=merge(theme, Theme(size=(1200, 880), Axis=(aspect=nothing,))),
        sim_labels,
        title,
    )

    return nothing

end

"""
    gasBarPlot(
        simulation_paths::Vector{String},
        slice::IndexType,
        quantity::Symbol,
        edges::Vector{<:Number};
        <keyword arguments>
    )::Nothing

Plot a bar plot of the gas fractions, where the bins are a given gas `quantity`.

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`. Each simulation will be plotted in a different figure.
  - `slice::IndexType`: Slice of the simulation, i.e. which snapshots will be plotted. It can be an integer (a single snapshot), a vector of integers (several snapshots), an `UnitRange` (e.g. 5:13), an `StepRange` (e.g. 5:2:13) or (:) (all snapshots). Starts at 1 and out of bounds indices are ignored.
  - `quantity::Symbol`: Target quantity. It can be any of the valid quantities of [`scatterQty`](@ref) with the cell/particle type `:gas`.
  - `edges::Vector{<:Number}`: A list of bin edges for `quantity`.
  - `components::Vector{Symbol}=[:ode_ionized, :ode_atomic, :ode_molecular_stellar]`: List of gas components to be considered. The fractions will be normalized to this list of components. See [`COMPONENTS`](@ref) for options.
  - `ylog::Bool=false`: If true, sets everything so the y axis is ``\\log_{10}``(`quantity`).
  - `output_path::String="."`: Path to the output folder.
  - `trans_mode::Union{Symbol,Tuple{TranslationType,RotationType,Dict{Symbol,Vector{String}}}}=:all_box`: How to translate and rotate the cells/particles, before filtering with `filter_mode`. For options see [`selectTransformation`](@ref).
  - `filter_mode::Union{Symbol,Tuple{Function,Dict{Symbol,Vector{String}}}}=:all`: Which cells/particles will be selected. For options see [`selectFilter`](@ref).
  - `da_ff::Function=filterNothing`: Filter function to be applied within [`daScatterGalaxy`](@ref) after `trans_mode` and `filter_mode` are applied. See the required signature and examples in `./src/analysis/filters.jl`.
  - `ff_request::Dict{Symbol,Vector{String}}=Dict{Symbol,Vector{String}}()`: Request dictionary for `da_ff`.
  - `title::Union{Symbol,<:AbstractString}=""`: Title for the figure. If left empty, no title is printed. It can also be set to one of the following options:

      + `:physical_time` -> Physical time since the Big Bang.
      + `:lookback_time` -> Physical time left to reach the last snapshot.
      + `:scale_factor`  -> Scale factor (only relevant for cosmological simulations).
      + `:redshift`      -> Redshift (only relevant for cosmological simulations).
  - `theme::Attributes=Theme()`: Plot theme that will take precedence over [`DEFAULT_THEME`](@ref).
"""
function gasBarPlot(
    simulation_paths::Vector{String},
    slice::IndexType,
    quantity::Symbol,
    edges::Vector{<:Number};
    components::Vector{Symbol}=[:ode_ionized, :ode_atomic, :ode_molecular_stellar],
    ylog::Bool=false,
    output_path::String=".",
    trans_mode::Union{Symbol,Tuple{TranslationType,RotationType,Dict{Symbol,Vector{String}}}}=:all_box,
    filter_mode::Union{Symbol,Tuple{Function,Dict{Symbol,Vector{String}}}}=:all,
    da_ff::Function=filterNothing,
    ff_request::Dict{Symbol,Vector{String}}=Dict{Symbol,Vector{String}}(),
    title::Union{Symbol,<:AbstractString}="",
    theme::Attributes=Theme(),
)::Nothing

    plot_params = plotParams(quantity)

    gas_requests = [plotParams(Symbol(component, :_mass)).request for component in components]

    base_request = mergeRequests(plot_params.request, ff_request, gas_requests...)

    translation, rotation, trans_request = selectTransformation(trans_mode, base_request)
    filter_function, request = selectFilter(filter_mode, trans_request)

    # Set arguments for the y axis
    if ylog
        y_unit       = Unitful.NoUnits
        y_exp_factor = 0
        yaxis_label  = L"\log_{10} \, " * getLabel("auto_label", 0, plot_params.unit)
        y_ticks      = log10.(ustrip.(plot_params.unit, edges))
    else
        y_unit       = plot_params.unit
        y_exp_factor = plot_params.exp_factor
        yaxis_label  = plot_params.axis_label
        y_ticks      = ustrip.(plot_params.unit, edges)
    end

    # Compute the number of bins
    n_bins  = length(edges) - 1

    # Compute the number of bars per bin
    n_bars = length(components)

    # Compute the dodge argument for `barplot!`
    dodge = repeat(1:n_bars, outer=n_bins)

    # Compute the axis ticks
    ticks = [string(round((y_ticks[i] + y_ticks[i + 1]) / 2, sigdigits=2)) for i in 1:n_bins]

    current_theme = merge(
        theme,
        Theme(
            size=(1200, 880),
            Legend=(nbanks=1, valign=:top),
            Axis=(
                aspect=nothing,
                limits=(nothing, 105, nothing, nothing),
                xticks=[0, 25, 50, 75, 100],
                yticks=(1:n_bins, ticks),
            ),
            BarPlot=(
                flip_labels_at=10,
                label_formatter=barPlotLabelFormater,
                label_size=35,
            ),
        ),
        DEFAULT_THEME,
    )

    # Set the color list
    colors = safeSelect(current_theme[:palette][:color][], dodge)

    for simulation_path in simulation_paths

        plotSnapshot(
            [simulation_path],
            request,
            [barplot!];
            pf_kwargs=[(; dodge, color=colors, direction=:x, bar_labels=:y)],
            output_path,
            base_filename="$(basename(simulation_path))_$(quantity)_vs_mass_fractions",
            slice,
            transform_box=true,
            translation,
            rotation,
            filter_function,
            da_functions=[daBarGasFractions],
            da_args=[(quantity, LinearGrid(edges; log=ylog))],
            da_kwargs=[(; components, filter_function=da_ff)],
            post_processing=ppBarPlotLabels,
            pp_args=(reverse(components),),
            pp_kwargs=(; colors=reverse(colors)),
            y_unit,
            y_exp_factor,
            yaxis_label,
            xaxis_label=L"\text{Mass fraction} \,\, [%]",
            yaxis_var_name=plot_params.var_name,
            theme=current_theme,
            title,
        )

    end

    return nothing

end

"""
    compareMolla2015(
        simulation_paths::Vector{String},
        slice::IndexType,
        quantity::Symbol;
        <keyword arguments>
    )::Nothing

Plot a profile with the corresponding experimental values of the Milky Way from Mollá et al. (2015).

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`. All the simulations will be plotted together.
  - `slice::IndexType`: Slice of the simulation, i.e. which snapshots will be plotted. It can be an integer (a single snapshot), a vector of integers (several snapshots), an `UnitRange` (e.g. 5:13), an `StepRange` (e.g. 5:2:13) or (:) (all snapshots). It works over the longest simulation. Starts at 1 and out of bounds indices are ignored.
  - `quantity::Symbol`: Quantity for the y axis. The options are:

      + `:stellar_area_density`               -> Stellar mass surface density.
      + `:sfr_area_density`                   -> Star formation rate surface density.
      + `:ode_molecular_area_density`         -> Molecular mass surface density.
      + `:ode_cold_area_density`              -> Cold mass (everything but ionized and atomic) surface density.
      + `:ode_molecular_stellar_area_density` -> Molecular + stellar mass surface density.
      + `:br_molecular_area_density`          -> Molecular mass surface density, computed using the pressure relation in Blitz et al. (2006).
      + `:ode_atomic_area_density`            -> Atomic mass surface density.
      + `:br_atomic_area_density`             -> Atomic mass surface density, computed using the pressure relation in Blitz et al. (2006).
      + `:O_stellar_abundance`                -> Stellar abundance of oxygen, as ``12 + \\log_{10}(\\mathrm{O \\, / \\, H})``.
      + `:N_stellar_abundance`                -> Stellar abundance of nitrogen, as ``12 + \\log_{10}(\\mathrm{N \\, / \\, H})``.
      + `:C_stellar_abundance`                -> Stellar abundance of carbon, as ``12 + \\log_{10}(\\mathrm{C \\, / \\, H})``.
  - `output_path::String="."`: Path to the output folder.
  - `trans_mode::Union{Symbol,Tuple{TranslationType,RotationType,Dict{Symbol,Vector{String}}}}=:all_box`: How to translate and rotate the cells/particles, before filtering with `filter_mode`. For options see [`selectTransformation`](@ref).
  - `filter_mode::Union{Symbol,Tuple{Function,Dict{Symbol,Vector{String}}}}=:all`: Which cells/particles will be selected. For options see [`selectFilter`](@ref).
  - `da_ff::Function=filterNothing`: Filter function to be applied within [`daScatterGalaxy`](@ref) after `trans_mode` and `filter_mode` are applied. See the required signature and examples in `./src/analysis/filters.jl`.
  - `ff_request::Dict{Symbol,Vector{String}}=Dict{Symbol,Vector{String}}()`: Request dictionary for `da_ff`.
  - `sim_labels::Union{Vector{<:AbstractString},Nothing}=basename.(simulation_paths)`: Labels for the plot legend, one per simulation. Set it to `nothing` if you don't want a legend.
  - `theme::Attributes=Theme()`: Plot theme that will take precedence over [`DEFAULT_THEME`](@ref).

# References

L. Blitz et al. (2006). *The Role of Pressure in GMC Formation II: The H2-Pressure Relation*. The Astrophysical Journal, **650(2)**, 933. [doi:10.1086/505417](https://doi.org/10.1086/505417)

M. Mollá et al. (2015). *Galactic chemical evolution: stellar yields and the initial mass function*. Monthly Notices of the Royal Astronomical Society **451(4)**, 3693–3708. [doi:10.1093/mnras/stv1102](https://doi.org/10.1093/mnras/stv1102)
"""
function compareMolla2015(
    simulation_paths::Vector{String},
    slice::IndexType,
    quantity::Symbol;
    output_path::String=".",
    trans_mode::Union{Symbol,Tuple{TranslationType,RotationType,Dict{Symbol,Vector{String}}}}=:all_box,
    filter_mode::Union{Symbol,Tuple{Function,Dict{Symbol,Vector{String}}}}=:all,
    da_ff::Function=filterNothing,
    ff_request::Dict{Symbol,Vector{String}}=Dict{Symbol,Vector{String}}(),
    sim_labels::Union{Vector{<:AbstractString},Nothing}=basename.(simulation_paths),
    theme::Attributes=Theme(),
)::Nothing

    plot_params = plotParams(quantity)

    if quantity ∈ STELLAR_ABUNDANCE
        extra_request = Dict(:stellar => ["POS "])
        yaxis_label = plot_params.axis_label
    else
        extra_request = Dict{Symbol,Vector{String}}()
        yaxis_label = L"\log_{10} \, " * getLabel("auto_label", 0, plot_params.unit)
    end

    base_request = mergeRequests(plot_params.request, extra_request, ff_request)

    translation, rotation, trans_request = selectTransformation(trans_mode, base_request)
    filter_function, request = selectFilter(filter_mode, trans_request)

    # Select the correct grid acording to the available data from M. Mollá et al. (2015)
    if quantity == :stellar_area_density

        grid = CircularGrid(16.5u"kpc", 14; shift=2.5u"kpc")

    elseif quantity ∈ [
        :ode_molecular_area_density,
        :ode_molecular_stellar_area_density,
        :br_molecular_area_density,
        :ode_cold_area_density,
        :sfr_area_density,
    ]

        grid = CircularGrid(19.5u"kpc", 20; shift=-0.5u"kpc")

    elseif quantity ∈ [:ode_atomic_area_density, :br_atomic_area_density]

        grid = CircularGrid(20.5u"kpc", 21; shift=-0.5u"kpc")

    elseif quantity == :O_stellar_abundance

        grid = CircularGrid(18.5u"kpc", 19; shift=-0.5u"kpc")

    elseif quantity == :N_stellar_abundance

        grid = CircularGrid(17.5u"kpc", 18; shift=-0.5u"kpc")

    elseif quantity == :C_stellar_abundance

        grid = CircularGrid(15.5u"kpc", 16; shift=-0.5u"kpc")

    else

        throw(ArgumentError("compareMolla2015: I don't recognize the quantity :$(quantity) for \
        Mollá et al. (2015)"))

    end

    if isone(length(simulation_paths))
        base_filename = "$(basename(simulation_paths[1]))_$(quantity)_profile_Molla2015"
    else
        base_filename = "$(quantity)_profile_Molla2015"
    end

    plotSnapshot(
        simulation_paths,
        request,
        [scatterlines!];
        output_path,
        base_filename,
        slice,
        transform_box=true,
        translation,
        rotation,
        filter_function,
        da_functions=[daMolla2015],
        da_args=[(grid, quantity)],
        da_kwargs=[(; y_unit=plot_params.unit, filter_function=da_ff)],
        post_processing=ppMolla2015!,
        pp_args=(quantity,),
        pp_kwargs=(; y_unit=plot_params.unit),
        x_unit=u"kpc",
        yaxis_label,
        xaxis_var_name=L"r",
        yaxis_var_name=plot_params.var_name,
        theme=merge(
            theme,
            Theme(
                size=(1400, 880),
                figure_padding=(10, 15, 5, 15),
                palette=(linestyle=[:solid],),
                Axis=(aspect=nothing,),
                Legend=(halign=:left, nbanks=1),
            ),
        ),
        sim_labels,
    )

    return nothing

end

"""
    compareAgertz2021(
        simulation_paths::Vector{String},
        slice::IndexType;
        <keyword arguments>
    )::Nothing

Plot a stellar density profile with the corresponding experimental values of the Milky Way from Agertz et al. (2021).

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`. All the simulations will be plotted together.
  - `slice::IndexType`: Slice of the simulation, i.e. which snapshots will be plotted. It can be an integer (a single snapshot), a vector of integers (several snapshots), an `UnitRange` (e.g. 5:13), an `StepRange` (e.g. 5:2:13) or (:) (all snapshots). It works over the longest simulation. Starts at 1 and out of bounds indices are ignored.
  - `ylog::Bool=true`: If true, sets the y axis as ``\\log_{10} \\Sigma_\\star``.
  - `output_path::String="."`: Path to the output folder.
  - `trans_mode::Union{Symbol,Tuple{TranslationType,RotationType,Dict{Symbol,Vector{String}}}}=:all_box`: How to translate and rotate the cells/particles, before filtering with `filter_mode`. For options see [`selectTransformation`](@ref).
  - `filter_mode::Union{Symbol,Tuple{Function,Dict{Symbol,Vector{String}}}}=:all`: Which cells/particles will be selected. For options see [`selectFilter`](@ref).
  - `da_ff::Function=filterNothing`: Filter function to be applied within [`daScatterGalaxy`](@ref) after `trans_mode` and `filter_mode` are applied. See the required signature and examples in `./src/analysis/filters.jl`.
  - `ff_request::Dict{Symbol,Vector{String}}=Dict{Symbol,Vector{String}}()`: Request dictionary for `da_ff`.
  - `sim_labels::Union{Vector{<:AbstractString},Nothing}=basename.(simulation_paths)`: Labels for the plot legend, one per simulation. Set it to `nothing` if you don't want a legend.
  - `theme::Attributes=Theme()`: Plot theme that will take precedence over [`DEFAULT_THEME`](@ref).

# References

O. Agertz et al. (2021). *VINTERGATAN – I. The origins of chemically, kinematically, and structurally distinct discs in a simulated Milky Way-mass galaxy*. Monthly Notices of the Royal Astronomical Society. **503(4), 5826–5845. [doi:10.1093/mnras/stab322](https://doi.org/10.1093/mnras/stab322)

A. K. Leroy et al. (2008). *THE STAR FORMATION EFFICIENCY IN NEARBY GALAXIES: MEASURING WHERE GAS FORMS STARS EFFECTIVELY*. The Astronomical Journal **136(6)**, 2782–2845. [doi:10.1088/0004-6256/136/6/2782](https://doi.org/10.1088/0004-6256/136/6/2782)

P. J. McMillan (2011). *Mass models of the Milky Way*. Monthly Notices of the Royal Astronomical Society **414(3)**, 2446–2457. [doi:10.1111/j.1365-2966.2011.18564.x](https://doi.org/10.1111/j.1365-2966.2011.18564.x)
"""
function compareAgertz2021(
    simulation_paths::Vector{String},
    slice::IndexType;
    ylog::Bool=true,
    output_path::String=".",
    trans_mode::Union{Symbol,Tuple{TranslationType,RotationType,Dict{Symbol,Vector{String}}}}=:all_box,
    filter_mode::Union{Symbol,Tuple{Function,Dict{Symbol,Vector{String}}}}=:all,
    da_ff::Function=filterNothing,
    ff_request::Dict{Symbol,Vector{String}}=Dict{Symbol,Vector{String}}(),
    sim_labels::Union{Vector{<:AbstractString},Nothing}=basename.(simulation_paths),
    theme::Attributes=Theme(),
)::Nothing

    plot_params = plotParams(:stellar_area_density)

    base_request = mergeRequests(plot_params.request, ff_request)

    translation, rotation, trans_request = selectTransformation(trans_mode, base_request)
    filter_function, request = selectFilter(filter_mode, trans_request)

    # Set arguments for the y axis
    if ylog
        y_unit      = Unitful.NoUnits
        y_log       = plot_params.unit
        yaxis_label = L"\log_{10} \, " * getLabel("auto_label", 0, plot_params.unit)
    else
        y_unit      = plot_params.unit
        y_log       = nothing
        yaxis_label = plot_params.axis_label
    end

    grid = CircularGrid(25.0u"kpc", 25)

    if isone(length(simulation_paths))
        base_filename = "$(basename(simulation_paths[1]))_stellar_mass_density_profile_Agertz2021"
    else
        base_filename = "stellar_mass_density_profile_Agertz2021"
    end

    plotSnapshot(
        simulation_paths,
        request,
        [lines!];
        output_path,
        base_filename,
        slice,
        transform_box=true,
        translation,
        rotation,
        filter_function,
        da_functions=[daProfile],
        da_args=[(:stellar_mass, grid)],
        da_kwargs=[(; y_log, density=true, filter_function=da_ff)],
        post_processing=ppAgertz2021!,
        pp_kwargs=(;
            galaxies=[:all, "MW"],
            colors=[WONG_PINK, WONG_BLUE],
            linestyle=:dash,
            y_unit=plot_params.unit,
        ),
        x_unit=u"kpc",
        y_unit,
        yaxis_label,
        xaxis_var_name=L"r",
        yaxis_var_name=plot_params.var_name,
        theme=merge(
            theme,
            Theme(
                size=(1200, 880),
                figure_padding=(5, 10, 5, 10),
                palette=(linestyle=[:solid], color=[WONG_ORANGE]),
                Legend=(nbanks=1, valign=:top, margin=(0, 10, 0, 10)),
                Axis=(aspect=nothing, xticks=0:5:25),
                Lines=(linewidth=3,),
                Scatter=(markersize=20,),
                Band=(alpha=0.6,),
            ),
        ),
        sim_labels,
    )

end

"""
    densityMap(
        simulation_paths::Vector{String},
        slice::IndexType;
        <keyword arguments>
    )::Nothing

Plot a 2D projection of the `component` mass density.

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`. Each simulation will be plotted in a different figure.
  - `slice::IndexType`: Slice of the simulation, i.e. which snapshots will be plotted. It can be an integer (a single snapshot), a vector of integers (several snapshots), an `UnitRange` (e.g. 5:13), an `StepRange` (e.g. 5:2:13) or (:) (all snapshots). Starts at 1 and out of bounds indices are ignored.
  - `components::Vector{Symbol}=[:gas]`: Target components. It can be any of the elements of [`COMPONENTS`](@ref).
  - `field_types::Vector{Symbol}=[:cells]`: If the field, of each of the `components`, is made up of `:particles` or Voronoi `:cells`.
  - `projection_planes::Vector{Symbol}=[:xy]`: Projection planes. The options are `:xy`, `:xz`, and `:yz`.
  - `box_size::Unitful.Length=100u"kpc"`: Physical side length of the plot window.
  - `pixel_length::Unitful.Length=0.1u"kpc"`: Side length of each bin.
  - `reduce_factor::Int=1`: Factor by which the resolution of the result will be reduced. This will be applied after the density projection, averaging the value of neighboring pixels. It has to divide the size of `grid` exactly.
  - `m_unit::Unitful.Units=u"Msun"`: Mass unit.
  - `l_unit::Unitful.Units=u"pc"`: Length unit.
  - `output_path::String="."`: Path to the output folder.
  - `trans_mode::Union{Symbol,Tuple{TranslationType,RotationType,Dict{Symbol,Vector{String}}}}=:all_box`: How to translate and rotate the cells/particles, before filtering with `filter_mode`. For options see [`selectTransformation`](@ref).
  - `filter_mode::Union{Symbol,Tuple{Function,Dict{Symbol,Vector{String}}}}=:all`: Which cells/particles will be selected. For options see [`selectFilter`](@ref).
  - `da_ff::Function=filterNothing`: Filter function to be applied within [`daScatterGalaxy`](@ref) after `trans_mode` and `filter_mode` are applied. See the required signature and examples in `./src/analysis/filters.jl`.
  - `ff_request::Dict{Symbol,Vector{String}}=Dict{Symbol,Vector{String}}()`: Request dictionary for `da_ff`.
  - `title::Union{Symbol,<:AbstractString}=""`: Title for the figure. If left empty, no title is printed. It can also be set to one of the following options:

      + `:physical_time` -> Physical time since the Big Bang.
      + `:lookback_time` -> Physical time left to reach the last snapshot.
      + `:scale_factor`  -> Scale factor (only relevant for cosmological simulations).
      + `:redshift`      -> Redshift (only relevant for cosmological simulations).
  - `annotation::AbstractString=""`: Text to be added into the top left corner of the plot. If left empty, nothing is printed.
  - `colorbar::Bool=false`: If a colorbar will be added.
  - `colorrange::Union{Nothing,Tuple{<:Real,<:Real}}=nothing`: Sets the start and end points of the colormap. Use `nothing` to use the extrema of the density values.
  - `theme::Attributes=Theme()`: Plot theme that will take precedence over [`DEFAULT_THEME`](@ref).
"""
function densityMap(
    simulation_paths::Vector{String},
    slice::IndexType;
    components::Vector{Symbol}=[:gas],
    field_types::Vector{Symbol}=[:cells],
    projection_planes::Vector{Symbol}=[:xy],
    box_size::Unitful.Length=100u"kpc",
    pixel_length::Unitful.Length=0.1u"kpc",
    reduce_factor::Int=1,
    m_unit::Unitful.Units=u"Msun",
    l_unit::Unitful.Units=u"pc",
    output_path::String=".",
    trans_mode::Union{Symbol,Tuple{TranslationType,RotationType,Dict{Symbol,Vector{String}}}}=:all_box,
    filter_mode::Union{Symbol,Tuple{Function,Dict{Symbol,Vector{String}}}}=:all,
    da_ff::Function=filterNothing,
    ff_request::Dict{Symbol,Vector{String}}=Dict{Symbol,Vector{String}}(),
    title::Union{Symbol,<:AbstractString}="",
    annotation::AbstractString="",
    colorbar::Bool=false,
    colorrange::Union{Nothing,Tuple{<:Real,<:Real}}=nothing,
    theme::Attributes=Theme(),
)::Nothing

    # Compute the axes limits, to avoid white padding around the heatmap grid
    limit = ustrip(u"kpc", box_size / 2.0)

    # Compute number of pixel per side
    resolution = round(Int, box_size / pixel_length)

    # Set up the grid
    grid = CubicGrid(box_size, resolution)

    for (i, component) in pairs(components)

        plot_params = plotParams(Symbol(component, :_mass))

        field_type = ring(field_types, i)

        if field_type == :cells
            extra_request = Dict(plot_params.cp_type => ["MASS", "RHO ", "POS "])
        elseif field_type == :particles
            extra_request = Dict(plot_params.cp_type => ["POS "])
        end

        base_request = mergeRequests(plot_params.request, extra_request, ff_request)

        translation, rotation, trans_request = selectTransformation(trans_mode, base_request)
        filter_function, request = selectFilter(filter_mode, trans_request)

        colorbar_pp = plotParams(Symbol(component, :_area_density))

        # Label for the colorbar
        colorbar_label = LaTeXString(
            L"\log_{10} \, " * getLabel(colorbar_pp.var_name, 0, m_unit * l_unit^-2)
        )

        for simulation_path in simulation_paths

            # Get the simulation name as a string
            sim_name = basename(simulation_path)

            for projection_plane in projection_planes

                # Construct the file name
                base_filename = "$(sim_name)_$(component)_$(projection_plane)_density_map"

                plotSnapshot(
                    [simulation_path],
                    request,
                    [heatmap!];
                    pf_kwargs=isnothing(colorrange) ? [(;)] : [(; colorrange)],
                    output_path,
                    base_filename,
                    slice,
                    transform_box=true,
                    translation,
                    rotation,
                    filter_function,
                    da_functions=[daDensity2DProjection],
                    da_args=[(grid, component, field_type)],
                    da_kwargs=[
                        (; projection_plane, reduce_factor, m_unit, l_unit, filter_function=da_ff),
                    ],
                    post_processing=isempty(annotation) ? getNothing : ppAnnotation!,
                    pp_args=(annotation,),
                    pp_kwargs=(; color=:white),
                    x_unit=u"kpc",
                    y_unit=u"kpc",
                    xaxis_var_name=string(projection_plane)[1:1],
                    yaxis_var_name=string(projection_plane)[2:2],
                    theme=merge(
                        theme,
                        Theme(
                            size=colorbar ? (880, 740) : (880, 880),
                            figure_padding=(5, 20, 20, 10),
                            Axis=(limits=(-limit, limit, -limit, limit),),
                            Colorbar=(label=colorbar_label,),
                        ),
                    ),
                    title,
                    colorbar,
                )

            end

        end

    end

    return nothing

end

"""
    densityMapVelField(
        simulation_paths::Vector{String},
        slice::IndexType;
        <keyword arguments>
    )::Nothing

Plot a 2D projection of the mass density of `component`, with the velocity field.

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`. Each simulation will be plotted in a different figure.
  - `slice::IndexType`: Slice of the simulation, i.e. which snapshots will be plotted. It can be an integer (a single snapshot), a vector of integers (several snapshots), an `UnitRange` (e.g. 5:13), an `StepRange` (e.g. 5:2:13) or (:) (all snapshots). Starts at 1 and out of bounds indices are ignored.
  - `components::Vector{Symbol}=[:gas]`: Target components. It can only be one of the elements of [`COMPONENTS`](@ref).
  - `field_types::Vector{Symbol}=[:cells]`: If the field, of each of the `components`, is made up of `:particles` or Voronoi `:cells`.
  - `projection_planes::Vector{Symbol}=[:xy]`: Projection planes. The options are `:xy`, `:xz`, and `:yz`.
  - `box_size::Unitful.Length=100u"kpc"`: Physical side length of the plot window.
  - `pixel_length::Unitful.Length=0.1u"kpc"`: Side length of each bin.
  - `reduce_factor::Int=1`: Factor by which the resolution of the result will be reduced. This will be applied after the density projection, averaging the value of neighboring pixels. It has to divide the size of `grid` exactly.
  - `m_unit::Unitful.Units=u"Msun"`: Mass unit.
  - `l_unit::Unitful.Units=u"pc"`: Length unit.
  - `v_unit::Unitful.Units=u"km * s^-1",`: Velocity unit
  - `output_path::String="."`: Path to the output folder.
  - `trans_mode::Union{Symbol,Tuple{TranslationType,RotationType,Dict{Symbol,Vector{String}}}}=:all_box`: How to translate and rotate the cells/particles, before filtering with `filter_mode`. For options see [`selectTransformation`](@ref).
  - `filter_mode::Union{Symbol,Tuple{Function,Dict{Symbol,Vector{String}}}}=:all`: Which cells/particles will be selected. For options see [`selectFilter`](@ref).
  - `da_ff::Function=filterNothing`: Filter function to be applied within [`daScatterGalaxy`](@ref) after `trans_mode` and `filter_mode` are applied. See the required signature and examples in `./src/analysis/filters.jl`.
  - `ff_request::Dict{Symbol,Vector{String}}=Dict{Symbol,Vector{String}}()`: Request dictionary for `da_ff`.
  - `title::Union{Symbol,<:AbstractString}=""`: Title for the figure. If left empty, no title is printed. It can also be set to one of the following options:

      + `:physical_time` -> Physical time since the Big Bang.
      + `:lookback_time` -> Physical time left to reach the last snapshot.
      + `:scale_factor`  -> Scale factor (only relevant for cosmological simulations).
      + `:redshift`      -> Redshift (only relevant for cosmological simulations).
  - `annotation::AbstractString=""`: Text to be added into the top left corner of the plot. If left empty, nothing is printed.
  - `colorbar::Bool=false`: If a colorbar will be added.
  - `colorrange::Union{Nothing,Tuple{<:Real,<:Real}}=nothing`: Sets the start and end points of the colormap. Use `nothing` to use the extrema of the density values.
  - `theme::Attributes=Theme()`: Plot theme that will take precedence over [`DEFAULT_THEME`](@ref).
"""
function densityMapVelField(
  simulation_paths::Vector{String},
  slice::IndexType;
  components::Vector{Symbol}=[:gas],
  field_types::Vector{Symbol}=[:cells],
  projection_planes::Vector{Symbol}=[:xy],
  box_size::Unitful.Length=100u"kpc",
  pixel_length::Unitful.Length=0.1u"kpc",
  reduce_factor::Int=1,
  m_unit::Unitful.Units=u"Msun",
  l_unit::Unitful.Units=u"pc",
  v_unit::Unitful.Units=u"km * s^-1",
  output_path::String=".",
  trans_mode::Union{Symbol,Tuple{TranslationType,RotationType,Dict{Symbol,Vector{String}}}}=:all_box,
  filter_mode::Union{Symbol,Tuple{Function,Dict{Symbol,Vector{String}}}}=:all,
  da_ff::Function=filterNothing,
  ff_request::Dict{Symbol,Vector{String}}=Dict{Symbol,Vector{String}}(),
  title::Union{Symbol,<:AbstractString}="",
  annotation::AbstractString="",
  colorbar::Bool=false,
  colorrange::Union{Nothing,Tuple{<:Real,<:Real}}=nothing,
  theme::Attributes=Theme(),
)::Nothing

    # Compute the axes limits, to avoid white padding around the heatmap grid
    limit = ustrip(u"kpc", box_size / 2.0)

    # Compute number of pixel per side
    resolution = round(Int, box_size / pixel_length)

    # Set up the grid for the heatmap
    grid_hm = CubicGrid(box_size, resolution)

    # Set up the grid for the velocity field
    grid_vf = SquareGrid(box_size, 20)

    for (i, component) in pairs(components)

        plot_params = plotParams(Symbol(component, :_mass))

        field_type = ring(field_types, i)

        if field_type == :cells
            extra_request = Dict(plot_params.cp_type => ["MASS", "RHO ", "POS ", "VEL "])
        elseif field_type == :particles
            extra_request = Dict(plot_params.cp_type => ["POS ", "VEL "])
        end

        base_request = mergeRequests(plot_params.request, extra_request, ff_request)

        translation, rotation, trans_request = selectTransformation(trans_mode, base_request)
        filter_function, request = selectFilter(filter_mode, trans_request)

        colorbar_pp = plotParams(Symbol(component, :_area_density))

        # Label for the colorbar
        colorbar_label = LaTeXString(
            L"\log_{10} \, " * getLabel(colorbar_pp.var_name, 0, m_unit * l_unit^-2)
        )

        for simulation_path in simulation_paths

            # Get the simulation name as a string
            sim_name = basename(simulation_path)

            for projection_plane in projection_planes

                # Construct the file name
                base_filename = "$(sim_name)_$(component)_$(projection_plane)_density_map"

                plotSnapshot(
                    [simulation_path, simulation_path],
                    request,
                    [heatmap!, arrows2d!];
                    pf_kwargs=isnothing(colorrange) ? [(;), (;)] : [(; colorrange), (;)],
                    output_path,
                    base_filename,
                    slice,
                    transform_box=true,
                    translation,
                    rotation,
                    filter_function,
                    da_functions=[daDensity2DProjection, daVelocityField],
                    da_args=[(grid_hm, component, field_type), (grid_vf, component)],
                    da_kwargs=[
                        (; projection_plane, reduce_factor, m_unit, l_unit, filter_function=da_ff),
                        (; projection_plane, v_unit, filter_function=da_ff)
                    ],
                    post_processing=isempty(annotation) ? getNothing : ppAnnotation!,
                    pp_args=(annotation,),
                    pp_kwargs=(; color=:white),
                    x_unit=u"kpc",
                    y_unit=u"kpc",
                    xaxis_var_name=string(projection_plane)[1:1],
                    yaxis_var_name=string(projection_plane)[2:2],
                    theme=merge(
                        theme,
                        Theme(
                            size=colorbar ? (880, 740) : (880, 880),
                            figure_padding=(5, 20, 20, 10),
                            Axis=(limits=(-limit, limit, -limit, limit),),
                            Colorbar=(label=colorbar_label,),
                        ),
                    ),
                    title,
                    colorbar,
                )

            end

        end

    end

    return nothing

end

"""
    gasSFRMap(
        simulation_paths::Vector{String},
        slice::IndexType;
        <keyword arguments>
    )::Nothing

Plot a 2D projection of the gas SFR.

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`. Each simulation will be plotted in a different figure.
  - `slice::IndexType`: Slice of the simulation, i.e. which snapshots will be plotted. It can be an integer (a single snapshot), a vector of integers (several snapshots), an `UnitRange` (e.g. 5:13), an `StepRange` (e.g. 5:2:13) or (:) (all snapshots). Starts at 1 and out of bounds indices are ignored.
  - `field_type::Symbol=:cells,`: If the field is made up of `:particles` or Voronoi `:cells`.
  - `projection_planes::Vector{Symbol}=[:xy]`: Projection planes. The options are `:xy`, `:xz`, and `:yz`.
  - `box_size::Unitful.Length=100u"kpc"`: Physical side length of the plot window.
  - `pixel_length::Unitful.Length=0.1u"kpc"`: Side length of each bin.
  - `reduce_factor::Int=1`: Factor by which the resolution of the result will be reduced. This will be applied after the density projection, averaging the value of neighboring pixels. It has to divide the size of `grid` exactly.
  - `m_unit::Unitful.Units=u"Msun"`: Mass unit.
  - `t_unit::Unitful.Units=u"yr"`: Time unit.
  - `output_path::String="."`: Path to the output folder.
  - `trans_mode::Union{Symbol,Tuple{TranslationType,RotationType,Dict{Symbol,Vector{String}}}}=:all_box`: How to translate and rotate the cells/particles, before filtering with `filter_mode`. For options see [`selectTransformation`](@ref).
  - `filter_mode::Union{Symbol,Tuple{Function,Dict{Symbol,Vector{String}}}}=:all`: Which cells/particles will be selected. For options see [`selectFilter`](@ref).
  - `da_ff::Function=filterNothing`: Filter function to be applied within [`daScatterGalaxy`](@ref) after `trans_mode` and `filter_mode` are applied. See the required signature and examples in `./src/analysis/filters.jl`.
  - `ff_request::Dict{Symbol,Vector{String}}=Dict{Symbol,Vector{String}}()`: Request dictionary for `da_ff`.
  - `title::Union{Symbol,<:AbstractString}=""`: Title for the figure. If left empty, no title is printed. It can also be set to one of the following options:

      + `:physical_time` -> Physical time since the Big Bang.
      + `:lookback_time` -> Physical time left to reach the last snapshot.
      + `:scale_factor`  -> Scale factor (only relevant for cosmological simulations).
      + `:redshift`      -> Redshift (only relevant for cosmological simulations).
  - `annotation::AbstractString=""`: Text to be added into the top left corner of the plot. If left empty, nothing is printed.
  - `colorbar::Bool=false`: If a colorbar will be added.
  - `colorrange::Union{Nothing,Tuple{<:Real,<:Real}}=nothing`: Sets the start and end points of the colormap. Use `nothing` to use the extrema of the density values.
  - `theme::Attributes=Theme()`: Plot theme that will take precedence over [`DEFAULT_THEME`](@ref).
"""
function gasSFRMap(
    simulation_paths::Vector{String},
    slice::IndexType;
    field_type::Symbol=:cells,
    projection_planes::Vector{Symbol}=[:xy],
    box_size::Unitful.Length=100u"kpc",
    pixel_length::Unitful.Length=0.1u"kpc",
    reduce_factor::Int=1,
    m_unit::Unitful.Units=u"Msun",
    t_unit::Unitful.Units=u"yr",
    output_path::String=".",
    trans_mode::Union{Symbol,Tuple{TranslationType,RotationType,Dict{Symbol,Vector{String}}}}=:all_box,
    filter_mode::Union{Symbol,Tuple{Function,Dict{Symbol,Vector{String}}}}=:all,
    da_ff::Function=filterNothing,
    ff_request::Dict{Symbol,Vector{String}}=Dict{Symbol,Vector{String}}(),
    title::Union{Symbol,<:AbstractString}="",
    annotation::AbstractString="",
    colorbar::Bool=false,
    colorrange::Union{Nothing,Tuple{<:Real,<:Real}}=nothing,
    theme::Attributes=Theme(),
)::Nothing

    # Compute the axes limits, to avoid white padding around the heatmap grid
    limit = ustrip(u"kpc", box_size / 2.0)

    # Compute number of pixel per side
    resolution = round(Int, box_size / pixel_length)

    # Set up the grid
    grid = CubicGrid(box_size, resolution)

    plot_params = plotParams(:gas_sfr)

    if field_type == :cells
        extra_request = Dict(:gas => ["MASS", "RHO ", "POS "])
    elseif field_type == :particles
        extra_request = Dict(:gas => ["POS "])
    end

    base_request = mergeRequests(plot_params.request, extra_request, ff_request)

    translation, rotation, trans_request = selectTransformation(trans_mode, base_request)
    filter_function, request = selectFilter(filter_mode, trans_request)

    # Label for the colorbar
    colorbar_label = LaTeXString(
        L"\log_{10} \, " * getLabel(plot_params.var_name, 0, m_unit * t_unit^-1)
    )

    for simulation_path in simulation_paths

        # Get the simulation name as a string
        sim_name = basename(simulation_path)

        for projection_plane in projection_planes

            # Construct the file name
            base_filename = "$(sim_name)_$(projection_plane)_gas_sfr_map"

            plotSnapshot(
                [simulation_path],
                request,
                [heatmap!];
                pf_kwargs=isnothing(colorrange) ? [(;)] : [(; colorrange)],
                output_path,
                base_filename,
                slice,
                transform_box=true,
                translation,
                rotation,
                filter_function,
                da_functions=[daGasSFR2DProjection],
                da_args=[(grid, field_type)],
                da_kwargs=[
                    (;
                        projection_plane,
                        reduce_factor,
                        m_unit,
                        t_unit,
                        filter_function=da_ff,
                    ),
                ],
                post_processing=isempty(annotation) ? getNothing : ppAnnotation!,
                pp_args=(annotation,),
                pp_kwargs=(; color=:white),
                x_unit=u"kpc",
                y_unit=u"kpc",
                xaxis_var_name=string(projection_plane)[1:1],
                yaxis_var_name=string(projection_plane)[2:2],
                theme=merge(
                    theme,
                    Theme(
                        size=colorbar ? (880, 740) : (880, 880),
                        figure_padding=(5, 20, 20, 10),
                        Axis=(limits=(-limit, limit, -limit, limit),),
                        Colorbar=(label=colorbar_label,),
                    ),
                ),
                title,
                colorbar,
            )

        end

    end

    return nothing

end

"""
    metallicityMap(
        simulation_paths::Vector{String},
        slice::IndexType;
        <keyword arguments>
    )::Nothing

Plot a 2D projection of the metallicity.

!!! note

    If if `element` = :all, the total metallicity (in solar units) is computed. If `element` = :X, the abundance of element X is computed, [`ABUNDANCE_SHIFT`](@ref) + ``\\log_{10}(X/H). In both cases the total value of the column given by the  line of sight of each pixel is computed.

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`. Each simulation will be plotted in a different figure.
  - `slice::IndexType`: Slice of the simulation, i.e. which snapshots will be plotted. It can be an integer (a single snapshot), a vector of integers (several snapshots), an `UnitRange` (e.g. 5:13), an `StepRange` (e.g. 5:2:13) or (:) (all snapshots). Starts at 1 and out of bounds indices are ignored.
  - `components::Vector{Symbol}=[:gas]`: Target component. It can be either `:stellar` or `:gas`.
  - `field_types::Vector{Symbol}=[:cells]`: If the field, of each of the `components`, is made up of `:particles` or Voronoi `:cells`.
  - `element::Symbol=:all`: Target element. The possibilities are the keys of [`ELEMENT_INDEX`](@ref). Set it to :all if you want the total metallicity.
  - `projection_planes::Vector{Symbol}=[:xy]`: Projection planes. The options are `:xy`, `:xz`, and `:yz`.
  - `box_size::Unitful.Length=100u"kpc"`: Physical side length of the plot window.
  - `pixel_length::Unitful.Length=0.1u"kpc"`: Side length of each bin.
  - `reduce_factor::Int=1`: Factor by which the resolution of the result will be reduced. This will be applied after the density projection, averaging the value of neighboring pixels. It has to divide the size of `grid` exactly.
  - `output_path::String="."`: Path to the output folder.
  - `trans_mode::Union{Symbol,Tuple{TranslationType,RotationType,Dict{Symbol,Vector{String}}}}=:all_box`: How to translate and rotate the cells/particles, before filtering with `filter_mode`. For options see [`selectTransformation`](@ref).
  - `filter_mode::Union{Symbol,Tuple{Function,Dict{Symbol,Vector{String}}}}=:all`: Which cells/particles will be selected. For options see [`selectFilter`](@ref).
  - `da_ff::Function=filterNothing`: Filter function to be applied within [`daScatterGalaxy`](@ref) after `trans_mode` and `filter_mode` are applied. See the required signature and examples in `./src/analysis/filters.jl`.
  - `ff_request::Dict{Symbol,Vector{String}}=Dict{Symbol,Vector{String}}()`: Request dictionary for `da_ff`.
  - `title::Union{Symbol,<:AbstractString}=""`: Title for the figure. If left empty, no title is printed. It can also be set to one of the following options:

      + `:physical_time` -> Physical time since the Big Bang.
      + `:lookback_time` -> Physical time left to reach the last snapshot.
      + `:scale_factor`  -> Scale factor (only relevant for cosmological simulations).
      + `:redshift`      -> Redshift (only relevant for cosmological simulations).
  - `annotation::AbstractString=""`: Text to be added into the top left corner of the plot. If left empty, nothing is printed.
  - `colorbar::Bool=false`: If a colorbar will be added.
  - `colorrange::Union{Nothing,Tuple{<:Real,<:Real}}=nothing`: Sets the start and end points of the colormap. Use `nothing` to use the extrema of the density values.
  - `theme::Attributes=Theme()`: Plot theme that will take precedence over [`DEFAULT_THEME`](@ref).
"""
function metallicityMap(
    simulation_paths::Vector{String},
    slice::IndexType;
    components::Vector{Symbol}=[:gas],
    field_types::Vector{Symbol}=[:cells],
    element::Symbol=:all,
    projection_planes::Vector{Symbol}=[:xy],
    box_size::Unitful.Length=100u"kpc",
    pixel_length::Unitful.Length=0.1u"kpc",
    reduce_factor::Int=1,
    output_path::String=".",
    trans_mode::Union{Symbol,Tuple{TranslationType,RotationType,Dict{Symbol,Vector{String}}}}=:all_box,
    filter_mode::Union{Symbol,Tuple{Function,Dict{Symbol,Vector{String}}}}=:all,
    da_ff::Function=filterNothing,
    ff_request::Dict{Symbol,Vector{String}}=Dict{Symbol,Vector{String}}(),
    title::Union{Symbol,<:AbstractString}="",
    annotation::AbstractString="",
    colorbar::Bool=false,
    colorrange::Union{Nothing,Tuple{<:Real,<:Real}}=nothing,
    theme::Attributes=Theme(),
)::Nothing

    # Compute the axes limits, to avoid white padding around the heatmap grid
    limit = ustrip(u"kpc", box_size / 2.0)

    # Compute number of pixel per side
    resolution = round(Int, box_size / pixel_length)

    # Set up the grid
    grid = CubicGrid(box_size, resolution)

    for (i, component) in pairs(components)

        if element == :all

            plot_params = plotParams(Symbol(component, :_metallicity))

            # Label for the colorbar
            colorbar_label = LaTeXString(
                L"\log_{10} \, " * getLabel(plot_params.var_name, 0, Unitful.NoUnits)
            )

        elseif haskey(ELEMENT_INDEX, element)

            plot_params = plotParams(Symbol(element, :_, component, :_abundance))

            # Label for the colorbarç
            colorbar_label = plot_params.axis_label

        else

            throw(ArgumentError("metallicityMap: The argument `element` can only be :all or one of \
            the keys of `ELEMENT_INDEX` (see `./src/constants/globals.jl`), but I got :$(element)"))

        end

        field_type = ring(field_types, i)

        if field_type == :cells
            extra_request = Dict(component => ["RHO ", "POS ", "MASS"])
        elseif field_type == :particles
            extra_request = Dict(component => ["POS ", "MASS"])
        end

        base_request = mergeRequests(plot_params.request, extra_request, ff_request)

        translation, rotation, trans_request = selectTransformation(trans_mode, base_request)
        filter_function, request = selectFilter(filter_mode, trans_request)

        for simulation_path in simulation_paths

            # Get the simulation name as a string
            sim_name = basename(simulation_path)

            for projection_plane in projection_planes

                # Construct the file name
                base_filename = "$(sim_name)_$(component)_$(projection_plane)_$(element)_metallicity_map"

                plotSnapshot(
                    [simulation_path],
                    request,
                    [heatmap!];
                    pf_kwargs=isnothing(colorrange) ? [(;)] : [(; colorrange)],
                    output_path,
                    base_filename,
                    slice,
                    transform_box=true,
                    translation,
                    rotation,
                    filter_function,
                    da_functions=[daMetallicity2DProjection],
                    da_args=[(grid, component, field_type)],
                    da_kwargs=[(; element, projection_plane, reduce_factor, filter_function=da_ff)],
                    post_processing=isempty(annotation) ? getNothing : ppAnnotation!,
                    pp_args=(annotation,),
                    pp_kwargs=(; color=:white),
                    x_unit=u"kpc",
                    y_unit=u"kpc",
                    xaxis_var_name=string(projection_plane)[1:1],
                    yaxis_var_name=string(projection_plane)[2:2],
                    theme=merge(
                        theme,
                        Theme(
                            size=colorbar ? (880, 740) : (880, 880),
                            figure_padding=(5, 20, 20, 10),
                            Axis=(limits=(-limit, limit, -limit, limit),),
                            Colorbar=(label=colorbar_label,),
                        ),
                    ),
                    title,
                    colorbar,
                )

            end

        end

    end

    return nothing

end

"""
    temperatureMap(
        simulation_paths::Vector{String},
        slice::IndexType;
        <keyword arguments>
    )::Nothing

Plot a 2D projection of the gas temperature.

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`. Each simulation will be plotted in a different figure.
  - `slice::IndexType`: Slice of the simulation, i.e. which snapshots will be plotted. It can be an integer (a single snapshot), a vector of integers (several snapshots), an `UnitRange` (e.g. 5:13), an `StepRange` (e.g. 5:2:13) or (:) (all snapshots). Starts at 1 and out of bounds indices are ignored.
  - `field_type::Symbol=:cells`: If the field is made up of `:particles` or Voronoi `:cells`.
  - `projection_planes::Vector{Symbol}=[:xy]`: Projection planes. The options are `:xy`, `:xz`, and `:yz`.
  - `box_size::Unitful.Length=100u"kpc"`: Physical side length of the plot window.
  - `pixel_length::Unitful.Length=0.1u"kpc"`: Side length of each bin.
  - `reduce_factor::Int=1`: Factor by which the resolution of the result will be reduced. This will be applied after the density projection, averaging the value of neighboring pixels. It has to divide the size of `grid` exactly.
  - `output_path::String="."`: Path to the output folder.
  - `trans_mode::Union{Symbol,Tuple{TranslationType,RotationType,Dict{Symbol,Vector{String}}}}=:all_box`: How to translate and rotate the cells/particles, before filtering with `filter_mode`. For options see [`selectTransformation`](@ref).
  - `filter_mode::Union{Symbol,Tuple{Function,Dict{Symbol,Vector{String}}}}=:all`: Which cells/particles will be selected. For options see [`selectFilter`](@ref).
  - `da_ff::Function=filterNothing`: Filter function to be applied within [`daScatterGalaxy`](@ref) after `trans_mode` and `filter_mode` are applied. See the required signature and examples in `./src/analysis/filters.jl`.
  - `ff_request::Dict{Symbol,Vector{String}}=Dict{Symbol,Vector{String}}()`: Request dictionary for `da_ff`.
  - `title::Union{Symbol,<:AbstractString}=""`: Title for the figure. If left empty, no title is printed. It can also be set to one of the following options:

      + `:physical_time` -> Physical time since the Big Bang.
      + `:lookback_time` -> Physical time left to reach the last snapshot.
      + `:scale_factor`  -> Scale factor (only relevant for cosmological simulations).
      + `:redshift`      -> Redshift (only relevant for cosmological simulations).
  - `annotation::AbstractString=""`: Text to be added into the top left corner of the plot. If left empty, nothing is printed.
  - `colorbar::Bool=false`: If a colorbar will be added.
  - `colorrange::Union{Nothing,Tuple{<:Real,<:Real}}=nothing`: Sets the start and end points of the colormap. Use `nothing` to use the extrema of the density values.
  - `theme::Attributes=Theme()`: Plot theme that will take precedence over [`DEFAULT_THEME`](@ref).
"""
function temperatureMap(
    simulation_paths::Vector{String},
    slice::IndexType;
    field_type::Symbol=:cells,
    projection_planes::Vector{Symbol}=[:xy],
    box_size::Unitful.Length=100u"kpc",
    pixel_length::Unitful.Length=0.1u"kpc",
    reduce_factor::Int=1,
    output_path::String=".",
    trans_mode::Union{Symbol,Tuple{TranslationType,RotationType,Dict{Symbol,Vector{String}}}}=:all_box,
    filter_mode::Union{Symbol,Tuple{Function,Dict{Symbol,Vector{String}}}}=:all,
    da_ff::Function=filterNothing,
    ff_request::Dict{Symbol,Vector{String}}=Dict{Symbol,Vector{String}}(),
    title::Union{Symbol,<:AbstractString}="",
    annotation::AbstractString="",
    colorbar::Bool=false,
    colorrange::Union{Nothing,Tuple{<:Real,<:Real}}=nothing,
    theme::Attributes=Theme(),
)::Nothing

    # Compute the axes limits, to avoid white padding around the heatmap grid
    limit = ustrip(u"kpc", box_size / 2.0)

    # Compute number of pixel per side
    resolution = round(Int, box_size / pixel_length)

    # Set up the grid
    grid = CubicGrid(box_size, resolution)

    plot_params = plotParams(:temperature)

    if field_type == :cells
        extra_request = Dict(:gas => ["RHO ", "POS "])
    elseif field_type == :particles
        extra_request = Dict(:gas => ["POS "])
    end

    base_request = mergeRequests(plot_params.request, extra_request, ff_request)

    translation, rotation, trans_request = selectTransformation(trans_mode, base_request)
    filter_function, request = selectFilter(filter_mode, trans_request)

    # Label for the colorbar
    colorbar_label = LaTeXString(
        L"\log_{10} \, " * getLabel(plot_params.var_name, 0, plot_params.unit)
    )

    for simulation_path in simulation_paths

        # Get the simulation name as a string
        sim_name = basename(simulation_path)

        for projection_plane in projection_planes

            # Construct the file name
            base_filename = "$(sim_name)_$(projection_plane)_gas_temperature_map"

            plotSnapshot(
                [simulation_path],
                request,
                [heatmap!];
                pf_kwargs=isnothing(colorrange) ? [(;)] : [(; colorrange)],
                output_path,
                base_filename,
                slice,
                transform_box=true,
                translation,
                rotation,
                filter_function,
                da_functions=[daTemperature2DProjection],
                da_args=[(grid, field_type)],
                da_kwargs=[(; projection_plane, reduce_factor, filter_function=da_ff)],
                post_processing=isempty(annotation) ? getNothing : ppAnnotation!,
                pp_args=(annotation,),
                pp_kwargs=(; color=:white),
                x_unit=u"kpc",
                y_unit=u"kpc",
                xaxis_var_name=string(projection_plane)[1:1],
                yaxis_var_name=string(projection_plane)[2:2],
                theme=merge(
                    theme,
                    Theme(
                        size=colorbar ? (880, 740) : (880, 880),
                        figure_padding=(5, 20, 20, 10),
                        Axis=(limits=(-limit, limit, -limit, limit),),
                        Colorbar=(label=colorbar_label,),
                    ),
                ),
                title,
                colorbar,
            )

        end

    end

    return nothing

end

"""
    timeSeries(
        simulation_paths::Vector{String},
        x_quantity::Symbol,
        y_quantity::Symbol;
        <keyword arguments>
    )::Nothing

Plot a time series.

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`. All the simulations will be plotted together.
  - `x_quantity::Symbol`: Quantity for the x axis. It can be any of the quantities valid for [`integrateQty`](@ref).
  - `y_quantity::Symbol`: Quantity for the y axis. It can be any of the quantities valid for [`integrateQty`](@ref).
  - `slice::IndexType=(:)`: Slice of the simulation, i.e. which snapshots will be plotted. It can be an integer (a single snapshot), a vector of integers (several snapshots), an `UnitRange` (e.g. 5:13), an `StepRange` (e.g. 5:2:13) or (:) (all snapshots). It works over the longest simulation. Starts at 1 and out of bounds indices are ignored.
  - `xlog::Bool=true`: If the x axis is will have a ``\\log_{10}`` scale.
  - `ylog::Bool=true`: If the y axis is will have a ``\\log_{10}`` scale.
  - `cumulative::Bool=false`: If the `y_quantity` will be accumulated or not.
  - `output_path::String="."`: Path to the output folder.
  - `trans_mode::Union{Symbol,Tuple{TranslationType,RotationType,Dict{Symbol,Vector{String}}}}=:all_box`: How to translate and rotate the cells/particles, before filtering with `filter_mode`. For options see [`selectTransformation`](@ref).
  - `filter_mode::Union{Symbol,Tuple{Function,Dict{Symbol,Vector{String}}}}=:all`: Which cells/particles will be selected. For options see [`selectFilter`](@ref).
  - `da_ff::Function=filterNothing`: Filter function to be applied within [`daScatterGalaxy`](@ref) after `trans_mode` and `filter_mode` are applied. See the required signature and examples in `./src/analysis/filters.jl`.
  - `ff_request::Dict{Symbol,Vector{String}}=Dict{Symbol,Vector{String}}()`: Request dictionary for `da_ff`.
  - `smooth::Int=0`: The result of [`integrateQty`](@ref) will be smoothed out using `smooth` bins. Set it to 0 if you want no smoothing.
  - `backup_results::Bool=false`: If the values to be plotted will be saved in a [JLD2](https://github.com/JuliaIO/JLD2.jl) file.
  - `theme::Attributes=Theme()`: Plot theme that will take precedence over [`DEFAULT_THEME`](@ref).
  - `sim_labels::Union{Vector{<:Union{AbstractString,Nothing}},Nothing}=nothing`: Labels for the plot legend, one per simulation. Set it to `nothing` if you don't want a legend.
"""
function timeSeries(
    simulation_paths::Vector{String},
    x_quantity::Symbol,
    y_quantity::Symbol;
    slice::IndexType=(:),
    xlog::Bool=true,
    ylog::Bool=true,
    cumulative::Bool=false,
    output_path::String=".",
    trans_mode::Union{Symbol,Tuple{TranslationType,RotationType,Dict{Symbol,Vector{String}}}}=:all_box,
    filter_mode::Union{Symbol,Tuple{Function,Dict{Symbol,Vector{String}}}}=:all,
    da_ff::Function=filterNothing,
    ff_request::Dict{Symbol,Vector{String}}=Dict{Symbol,Vector{String}}(),
    smooth::Int=0,
    backup_results::Bool=false,
    theme::Attributes=Theme(),
    sim_labels::Union{Vector{<:Union{AbstractString,Nothing}},Nothing}=nothing,
)::Nothing

    integration_functions = (dd->integrateQty(dd, x_quantity), dd->integrateQty(dd, y_quantity))

    timeSeries(
        simulation_paths,
        plotParams(x_quantity),
        plotParams(y_quantity),
        integration_functions;
        slice,
        xlog,
        ylog,
        cumulative,
        file_name="$(y_quantity)_vs_$(x_quantity)",
        output_path,
        trans_mode,
        filter_mode,
        da_ff,
        ff_request,
        smooth,
        backup_results,
        theme,
        sim_labels,
    )

    return nothing

end

"""
    timeSeries(
        simulation_paths::Vector{String},
        x_plot_params::PlotParams,
        y_plot_params::PlotParams,
        integration_functions::NTuple{2,Function};
        <keyword arguments>
    )::Nothing

Plot a time series.

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`. All the simulations will be plotted together.
  - `x_plot_params::PlotParams`: Plot parameters for the x axis.
  - `y_plot_params::PlotParams`: Plot parameters for the y axis.
  - `integration_functions::NTuple{2,Function}`: Functions to compute the integral value of the x and y quantities at a given time. The functions must have the signature:

    `integration_functions(data_dict::Dict)::Number`

    where

      + `data_dict::Dict`: Data dictionary (see [`makeDataDict`](@ref) for the canonical description).
  - `slice::IndexType=(:)`: Slice of the simulation, i.e. which snapshots will be plotted. It can be an integer (a single snapshot), a vector of integers (several snapshots), an `UnitRange` (e.g. 5:13), an `StepRange` (e.g. 5:2:13) or (:) (all snapshots). It works over the longest simulation. Starts at 1 and out of bounds indices are ignored.
  - `xlog::Bool=true`: If the x axis is will have a ``\\log_{10}`` scale.
  - `ylog::Bool=true`: If the y axis is will have a ``\\log_{10}`` scale.
  - `cumulative::Bool=false`: If the `y_quantity` will be accumulated or not.
  - `file_name::String="time_series"`: Name of the output file (without extension).
  - `output_path::String="."`: Path to the output folder.
  - `trans_mode::Union{Symbol,Tuple{TranslationType,RotationType,Dict{Symbol,Vector{String}}}}=:all_box`: How to translate and rotate the cells/particles, before filtering with `filter_mode`. For options see [`selectTransformation`](@ref).
  - `filter_mode::Union{Symbol,Tuple{Function,Dict{Symbol,Vector{String}}}}=:all`: Which cells/particles will be selected. For options see [`selectFilter`](@ref).
  - `da_ff::Function=filterNothing`: Filter function to be applied within [`daScatterGalaxy`](@ref) after `trans_mode` and `filter_mode` are applied. See the required signature and examples in `./src/analysis/filters.jl`.
  - `ff_request::Dict{Symbol,Vector{String}}=Dict{Symbol,Vector{String}}()`: Request dictionary for `da_ff`.
  - `smooth::Int=0`: The result of [`integrateQty`](@ref) will be smoothed out using `smooth` bins. Set it to 0 if you want no smoothing.
  - `backup_results::Bool=false`: If the values to be plotted will be saved in a [JLD2](https://github.com/JuliaIO/JLD2.jl) file.
  - `theme::Attributes=Theme()`: Plot theme that will take precedence over [`DEFAULT_THEME`](@ref).
  - `sim_labels::Union{Vector{<:Union{AbstractString,Nothing}},Nothing}=nothing`: Labels for the plot legend, one per simulation. Set it to `nothing` if you don't want a legend.
"""
function timeSeries(
    simulation_paths::Vector{String},
    x_plot_params::PlotParams,
    y_plot_params::PlotParams,
    integration_functions::NTuple{2,Function};
    slice::IndexType=(:),
    xlog::Bool=true,
    ylog::Bool=true,
    cumulative::Bool=false,
    file_name::String="time_series",
    output_path::String=".",
    trans_mode::Union{Symbol,Tuple{TranslationType,RotationType,Dict{Symbol,Vector{String}}}}=:all_box,
    filter_mode::Union{Symbol,Tuple{Function,Dict{Symbol,Vector{String}}}}=:all,
    da_ff::Function=filterNothing,
    ff_request::Dict{Symbol,Vector{String}}=Dict{Symbol,Vector{String}}(),
    smooth::Int=0,
    backup_results::Bool=false,
    theme::Attributes=Theme(),
    sim_labels::Union{Vector{<:Union{AbstractString,Nothing}},Nothing}=nothing,
)::Nothing

    qty_request = mergeRequests(x_plot_params.request, y_plot_params.request)

    # Set arguments for the x axis
    if xlog
        x_log        = x_plot_params.unit
        x_unit       = Unitful.NoUnits
        x_exp_factor = 0
        xaxis_label  = L"\log_{10} \, " * getLabel("auto_label", 0, x_plot_params.unit)
    else
        x_log        = nothing
        x_unit       = x_plot_params.unit
        x_exp_factor = x_plot_params.exp_factor
        xaxis_label  = x_plot_params.axis_label
    end

    # Set arguments for the y axis
    if ylog
        y_log        = y_plot_params.unit
        y_unit       = Unitful.NoUnits
        y_exp_factor = 0
        yaxis_label  = L"\log_{10} \, " * getLabel("auto_label", 0, y_plot_params.unit)
    else
        y_log        = nothing
        y_unit       = y_plot_params.unit
        y_exp_factor = y_plot_params.exp_factor
        yaxis_label  = y_plot_params.axis_label
    end

    if isone(length(simulation_paths))
        filename = "$(basename(simulation_paths[1]))_$(file_name)"
    else
        filename = file_name
    end

    plotTimeSeries(
        simulation_paths,
        [lines!];
        output_path,
        filename,
        slice,
        da_functions=[daEvolution],
        da_args=[(qty_request, integration_functions)],
        da_kwargs=[
            (;
                trans_mode,
                filter_mode,
                extra_filter=da_ff,
                ff_request,
                x_log,
                y_log,
                smooth,
                cumulative,
            )
        ],
        x_unit,
        y_unit,
        x_exp_factor,
        y_exp_factor,
        xaxis_label,
        yaxis_label,
        xaxis_var_name=x_plot_params.var_name,
        yaxis_var_name=y_plot_params.var_name,
        save_figures=!backup_results,
        backup_results,
        theme,
        sim_labels,
    )

    return nothing

end

"""
    statisticsEvolution(
        simulation_paths::Vector{String},
        x_quantity::Symbol,
        y_quantity::Symbol;
        <keyword arguments>
    )::Nothing

Plot a time series of the statistics of `y_quantity` (25th, 50th, 75th percentails, and maximum and minimum).

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`. Each simulation will be plotted in a different figure.
  - `x_quantity::Symbol`: Quantity for the x axis. It can be any of the time quantities valid for [`integrateQty`](@ref), namely

      + `:physical_time` -> Physical time since the Big Bang.
      + `:lookback_time` -> Physical time left to reach the last snapshot.
      + `:scale_factor`  -> Scale factor (only relevant for cosmological simulations).
      + `:redshift`      -> Redshift (only relevant for cosmological simulations).
  - `y_quantity::Symbol`: Quantity for the y axis. It can be any of the quantities valid for [`integrateQty`](@ref).
  - `slice::IndexType=(:)`: Slice of the simulation, i.e. which snapshots will be plotted. It can be an integer (a single snapshot), a vector of integers (several snapshots), an `UnitRange` (e.g. 5:13), an `StepRange` (e.g. 5:2:13) or (:) (all snapshots). It works over the longest simulation. Starts at 1 and out of bounds indices are ignored.
  - `xlog::Bool=false`: If the x axis is will have a ``\\log_{10}`` scale.
  - `ylog::Bool=false`: If the y axis is will have a ``\\log_{10}`` scale.
  - `cumulative::Bool=false`: If the `y_quantity` will be accumulated or not.
  - `output_path::String="."`: Path to the output folder.
  - `trans_mode::Union{Symbol,Tuple{TranslationType,RotationType,Dict{Symbol,Vector{String}}}}=:all_box`: How to translate and rotate the cells/particles, before filtering with `filter_mode`. For options see [`selectTransformation`](@ref).
  - `filter_mode::Union{Symbol,Tuple{Function,Dict{Symbol,Vector{String}}}}=:all`: Which cells/particles will be selected. For options see [`selectFilter`](@ref).
  - `da_ff::Function=filterNothing`: Filter function to be applied within [`daScatterGalaxy`](@ref) after `trans_mode` and `filter_mode` are applied. See the required signature and examples in `./src/analysis/filters.jl`.
  - `ff_request::Dict{Symbol,Vector{String}}=Dict{Symbol,Vector{String}}()`: Request dictionary for `da_ff`.
  - `smooth::Int=0`: The result of [`integrateQty`](@ref) will be smoothed out using `smooth` bins. Set it to 0 if you want no smoothing.
  - `backup_results::Bool=false`: If the values to be plotted will be saved in a [JLD2](https://github.com/JuliaIO/JLD2.jl) file.
  - `theme::Attributes=Theme()`: Plot theme that will take precedence over [`DEFAULT_THEME`](@ref).
"""
function statisticsEvolution(
    simulation_paths::Vector{String},
    x_quantity::Symbol,
    y_quantity::Symbol;
    slice::IndexType=(:),
    xlog::Bool=false,
    ylog::Bool=false,
    cumulative::Bool=false,
    output_path::String=".",
    trans_mode::Union{Symbol,Tuple{TranslationType,RotationType,Dict{Symbol,Vector{String}}}}=:all_box,
    filter_mode::Union{Symbol,Tuple{Function,Dict{Symbol,Vector{String}}}}=:all,
    da_ff::Function=filterNothing,
    ff_request::Dict{Symbol,Vector{String}}=Dict{Symbol,Vector{String}}(),
    smooth::Int=0,
    backup_results::Bool=false,
    theme::Attributes=Theme(),
)::Nothing

    (
        x_quantity ∈ [:scale_factor, :redshift, :physical_time, :lookback_time] ) ||
        throw(ArgumentError("statisticsEvolution: The argument `x_quantity` has to be \
        :scale_factor, :redshift, :physical_time or :lookback_time, but I got :$(x_quantity)")
    )

    x_plot_params = plotParams(x_quantity)
    y_plot_params = plotParams(y_quantity)

    # Set arguments for the x axis
    if xlog
        x_log        = x_plot_params.unit
        x_unit       = Unitful.NoUnits
        x_exp_factor = 0
        xaxis_label  = L"\log_{10} \, " * getLabel("auto_label", 0, x_plot_params.unit)
    else
        x_log        = nothing
        x_unit       = x_plot_params.unit
        x_exp_factor = x_plot_params.exp_factor
        xaxis_label  = x_plot_params.axis_label
    end

    # Set arguments for the y axis
    if ylog
        y_log        = y_plot_params.unit
        y_unit       = Unitful.NoUnits
        y_exp_factor = 0
        yaxis_label  = L"\log_{10} \, " * getLabel("auto_label", 0, y_plot_params.unit)
    else
        y_log        = nothing
        y_unit       = y_plot_params.unit
        y_exp_factor = y_plot_params.exp_factor
        yaxis_label  = y_plot_params.axis_label
    end

    y_agg_functions = [
        dd->integrateQty(dd, y_quantity; agg_function=x->percentile(x, 25)),
        dd->integrateQty(dd, y_quantity; agg_function=median),
        dd->integrateQty(dd, y_quantity; agg_function=x->percentile(x, 75)),
        dd->integrateQty(dd, y_quantity; agg_function=minimum),
        dd->integrateQty(dd, y_quantity; agg_function=maximum),
    ]
    labels     = ["25th Percentile", "50th Percentile", "75th Percentile", "Minimum", "Maximum"]
    linestyles = [:solid, :solid, :solid, :dash, :dash]
    colors     = [WONG_BLUE, WONG_GREEN, WONG_RED, :black, :black]
    request    = mergeRequests(x_plot_params.request, y_plot_params.request)

    for simulation_path in simulation_paths

        plotTimeSeries(
            fill(simulation_path, 5),
            [lines!];
            output_path,
            filename="$(basename(simulation_path))_$(y_quantity)_vs_$(x_quantity)",
            slice,
            da_functions=fill(daEvolution, 5),
            da_args=[
                (request, (dd->integrateQty(dd, x_quantity), y_agg_function))
                for y_agg_function in y_agg_functions
            ],
            da_kwargs=[
                (;
                    trans_mode,
                    filter_mode,
                    extra_filter=da_ff,
                    ff_request,
                    x_log,
                    y_log,
                    smooth,
                    cumulative,
                ),
            ],
            x_unit,
            y_unit,
            x_exp_factor,
            y_exp_factor,
            xaxis_label,
            yaxis_label,
            xaxis_var_name=x_plot_params.var_name,
            yaxis_var_name=y_plot_params.var_name,
            save_figures=!backup_results,
            backup_results,
            theme=merge(
                theme,
                Theme(
                    size=(1400, 880),
                    figure_padding=(10, 15, 5, 15),
                    palette=(linestyle=linestyles, color=colors),
                    Axis=(aspect=nothing,),
                    Lines=(linewidth=3,),
                    Legend=(nbanks=1,),
                ),
            ),
            sim_labels=labels,
        )

    end

    return nothing

end

"""
    gasEvolution(
        simulation_paths::Vector{String};
        <keyword arguments>
    )::Nothing

Plot a time series of the gas components. Either their masses or their fractions.

!!! note

    The molecular component includes the molecular and stellar fractions from our star formation model.

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`. Each simulation will be plotted in a different figure.
  - `fractions::Bool=true`: If the fractions (default), or the masses, will be plotted.
  - `slice::IndexType=(:)`: Slice of the simulation, i.e. which snapshots will be plotted. It can be an integer (a single snapshot), a vector of integers (several snapshots), an `UnitRange` (e.g. 5:13), an `StepRange` (e.g. 5:2:13) or (:) (all snapshots). Starts at 1 and out of bounds indices are ignored.
  - `ylog::Bool=true`: If the y axis is will have a ``\\log_{10}`` scale.
  - `output_path::String="."`: Path to the output folder.
  - `trans_mode::Union{Symbol,Tuple{TranslationType,RotationType,Dict{Symbol,Vector{String}}}}=:all_box`: How to translate and rotate the cells/particles, before filtering with `filter_mode`. For options see [`selectTransformation`](@ref).
  - `filter_mode::Union{Symbol,Tuple{Function,Dict{Symbol,Vector{String}}}}=:all`: Which cells/particles will be selected. For options see [`selectFilter`](@ref).
  - `da_ff::Function=filterNothing`: Filter function to be applied within [`daScatterGalaxy`](@ref) after `trans_mode` and `filter_mode` are applied. See the required signature and examples in `./src/analysis/filters.jl`.
  - `ff_request::Dict{Symbol,Vector{String}}=Dict{Symbol,Vector{String}}()`: Request dictionary for `da_ff`.
  - `smooth::Int=0`: The result of [`integrateQty`](@ref) will be smoothed out using `smooth` bins. Set it to 0 if you want no smoothing.
  - `backup_results::Bool=false`: If the values to be plotted will be saved in a [JLD2](https://github.com/JuliaIO/JLD2.jl) file.
  - `theme::Attributes=Theme()`: Plot theme that will take precedence over [`DEFAULT_THEME`](@ref).
"""
function gasEvolution(
    simulation_paths::Vector{String};
    fractions::Bool=true,
    slice::IndexType=(:),
    ylog::Bool=true,
    output_path::String=".",
    trans_mode::Union{Symbol,Tuple{TranslationType,RotationType,Dict{Symbol,Vector{String}}}}=:all_box,
    filter_mode::Union{Symbol,Tuple{Function,Dict{Symbol,Vector{String}}}}=:all,
    da_ff::Function=filterNothing,
    ff_request::Dict{Symbol,Vector{String}}=Dict{Symbol,Vector{String}}(),
    smooth::Int=0,
    backup_results::Bool=false,
    theme::Attributes=Theme(),
)::Nothing

    if fractions
        quantity = :fraction
    else
        quantity = :mass
    end

    x_plot_params = plotParams(:physical_time)
    y_plot_params = plotParams(quantity)

    quantities = Symbol.(
        [:ode_ionized, :ode_atomic, :ode_molecular_stellar, :ode_metals, :ode_dust],
        :_,
        quantity,
    )

    sim_labels = ["Ionized ", "Atomic ", "Molecular + stars ", "Metal ", "Dust "] .* string(quantity)

    # Set arguments for the y axis
    if ylog
        y_log        = y_plot_params.unit
        y_unit       = Unitful.NoUnits
        y_exp_factor = 0
        yaxis_label  = L"\log_{10} \, " * getLabel("auto_label", 0, y_plot_params.unit)
    else
        y_log        = nothing
        y_unit       = y_plot_params.unit
        y_exp_factor = y_plot_params.exp_factor
        yaxis_label  = y_plot_params.axis_label
    end

    for simulation_path in simulation_paths

        if fractions
            filename = "gas_fractions_vs_physical_time_$(basename(simulation_path))"
        else
            filename = "gas_masses_vs_physical_time_$(basename(simulation_path))"
        end

        plotTimeSeries(
            fill(simulation_path, length(quantities)),
            [lines!];
            output_path,
            filename,
            slice,
            da_functions=[daEvolution],
            da_args=[(:physical_time, quantity) for quantity in quantities],
            da_kwargs=[(; trans_mode, filter_mode, extra_filter=da_ff, ff_request, y_log, smooth)],
            x_unit=x_plot_params.unit,
            y_unit,
            x_exp_factor=x_plot_params.exp_factor,
            y_exp_factor,
            xaxis_label=x_plot_params.axis_label,
            yaxis_label,
            xaxis_var_name=x_plot_params.var_name,
            yaxis_var_name=y_plot_params.var_name,
            save_figures=!backup_results,
            backup_results,
            theme=merge(
                theme,
                Theme(
                    size=(1400, 880),
                    figure_padding=(10, 15, 5, 15),
                    palette=(linestyle=[:solid],),
                    Axis=(aspect=nothing, xticks=0:14),
                ),
            ),
            sim_labels,
        )

    end

    return nothing

end

"""
    gasFractionsEvolution(
        simulation_paths::Vector{String};
        <keyword arguments>
    )::Nothing

Plot time evolution of the masses and fractions of the gas components, in two panels.

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`. Each simulation will be plotted in a different figure.
  - `slice::IndexType=(:)`: Slice of the simulation, i.e. which snapshots will be plotted. It can be an integer (a single snapshot), a vector of integers (several snapshots), an `UnitRange` (e.g. 5:13), an `StepRange` (e.g. 5:2:13) or (:) (all snapshots). Starts at 1 and out of bounds indices are ignored.
  - `r_gas::Unitful.Length=DISK_R`: Radius of the gas sphere to consider.
  - `output_path::String="."`: Path to the output folder.
  - `mass_limits::NTuple{2,Float64}=(-4.2, 1.2)`: Limits for the masses, ``\\log_{10} M \\mathrm{[M_\\odot}``.
  - `fraction_limits::NTuple{2,Float64}=(-5.2, 0.2)`: Limits for fractions, ``\\log_{10} f``.
  - `trans_mode::Union{Symbol,Tuple{TranslationType,RotationType,Dict{Symbol,Vector{String}}}}=:all_box`: How to translate and rotate the cells/particles, before filtering with `filter_mode`. For options see [`selectTransformation`](@ref).
  - `filter_mode::Union{Symbol,Tuple{Function,Dict{Symbol,Vector{String}}}}=:all`: Which cells/particles will be selected. For options see [`selectFilter`](@ref).
  - `theme::Attributes=Theme()`: Plot theme that will take precedence over [`DEFAULT_THEME`](@ref).
"""
function gasFractionsEvolution(
    simulation_paths::Vector{String};
    slice::IndexType=(:),
    r_gas::Unitful.Length=DISK_R,
    output_path::String=".",
    mass_limits::NTuple{2,Float64}=(-4.2, 1.2),
    fraction_limits::NTuple{2,Float64}=(-5.2, 0.2),
    trans_mode::Union{Symbol,Tuple{TranslationType,RotationType,Dict{Symbol,Vector{String}}}}=:all_box,
    filter_mode::Union{Symbol,Tuple{Function,Dict{Symbol,Vector{String}}}}=:all,
    theme::Attributes=Theme(),
)::Nothing

    for simulation_path in simulation_paths

        if isSimSFM(simulation_path)
            quantities = [:ode_ionized, :ode_atomic, :ode_molecular_stellar, :ode_metals, :ode_dust]
            labels = ["ODE HII", "ODE HI", "ODE H2 + stars", "ODE Z", "ODE dust"]
        else
            quantities = [:ionized, :br_atomic, :br_molecular, :Z_gas]
            labels = ["HII", "BR HI", "BR H2", "Z"]
        end

        colors = [WONG_BLUE, WONG_PINK, WONG_GREEN, WONG_CELESTE, WONG_RED]

        x_plot_params = plotParams(:physical_time)
        y_plot_params = plotParams(:fraction)

        temp_folder = joinpath(output_path, "_gas_evolution")

        plotTimeSeries(
            fill(simulation_path, length(quantities)),
            [lines!];
            output_path=temp_folder,
            slice,
            filename="fraction_evolution",
            da_functions=[daEvolution],
            da_args=[(:physical_time, Symbol(quantity, :_fraction)) for quantity in quantities],
            da_kwargs=[(;
                trans_mode,
                filter_mode,
                extra_filter=dd -> filterBySphere(dd, 0.0u"kpc", r_gas, :zero),
                ff_request=Dict(:gas => ["POS "]),
            )],
            x_unit=x_plot_params.unit,
            y_unit=y_plot_params.unit,
            x_exp_factor=x_plot_params.exp_factor,
            y_exp_factor=y_plot_params.exp_factor,
            save_figures=false,
            backup_results=true,
            sim_labels=string.(quantities),
        )

        fraction_label = LaTeXString(
            L"\log_{10} \, " * getLabel(
                y_plot_params.var_name,
                y_plot_params.exp_factor,
                y_plot_params.unit,
            )
        )

        y_plot_params = plotParams(:mass)

        plotTimeSeries(
            fill(simulation_path, length(quantities)),
            [lines!];
            output_path=temp_folder,
            slice,
            filename="mass_evolution",
            da_functions=[daEvolution],
            da_args=[(:physical_time, Symbol(quantity, :_mass)) for quantity in quantities],
            da_kwargs=[(;
                trans_mode,
                filter_mode,
                extra_filter=dd -> filterBySphere(dd, 0.0u"kpc", r_gas, :zero),
                ff_request=Dict(:gas => ["POS "]),
            )],
            x_unit=x_plot_params.unit,
            y_unit=y_plot_params.unit,
            x_exp_factor=x_plot_params.exp_factor,
            y_exp_factor=y_plot_params.exp_factor,
            save_figures=false,
            backup_results=true,
            sim_labels=string.(quantities),
        )

        mass_label = LaTeXString(
            L"\log_{10} \, " * getLabel(
                y_plot_params.var_name,
                y_plot_params.exp_factor,
                y_plot_params.unit,
            )
        )

        jld2_paths = joinpath.(temp_folder, ["mass_evolution.jld2", "fraction_evolution.jld2"])

        current_theme = merge(
            theme,
            Theme(
                size=(880, 1650),
                figure_padding=(5, 10, 5, 10),
                palette=(linestyle=[:solid],),
                Axis=(xticks=0:14,),
            ),
            DEFAULT_THEME,
            theme_latexfonts(),
        )

        with_theme(current_theme) do

            f = Figure()

            ax_1 = CairoMakie.Axis(
                f[1, 1];
                xlabel=L"t \, [\mathrm{Gyr}]",
                ylabel=mass_label,
                xlabelvisible=false,
                xticklabelsvisible=false,
                limits=(-0.1, nothing, mass_limits...),
            )

            jldopen(jld2_paths[1], "r") do mass_evolution

                for (quantity, color) in zip(quantities, colors)

                    x, y = mass_evolution["mass_evolution"][string(quantity)]

                    lines!(ax_1, x, log10.(y); color)

                end

            end

            ax_2 = CairoMakie.Axis(
                f[2, 1];
                xlabel=L"t \, [\mathrm{Gyr}]",
                ylabel=fraction_label,
                aspect=nothing,
                limits=(-0.1, nothing, fraction_limits...),
            )

            jldopen(jld2_paths[2], "r") do fraction_evolution

                for (quantity, color, label) in zip(quantities, colors, labels)

                    x, y = fraction_evolution["fraction_evolution"][string(quantity)]

                    lines!(ax_2, x, log10.(y); color, label)

                end

            end

            axislegend(ax_2, position=:rb, framevisible=false, nbanks=2)

            linkxaxes!(ax_1, ax_2)

            save(joinpath(output_path, "$(basename(simulation_path))_gas_evolution.png"), f)

        end

        rm(temp_folder; recursive=true)

    end

    return nothing

end

"""
    sfrTXT(
        simulation_paths::Vector{String},
        x_quantity::Symbol,
        y_quantity::Symbol;
        <keyword arguments>
    )::Nothing

Plot a time series of the data in the `sfr.txt` file.

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`. All the simulations will be plotted together.
  - `x_quantity::Symbol`: Quantity for the x axis. The options are:

      + `:physical_time` -> Physical time since the Big Bang.
      + `:lookback_time` -> Physical time left to reach the last snapshot.
      + `:scale_factor`  -> Scale factor.
      + `:redshift`      -> Redshift.
  - `y_quantity::Symbol`: Quantity for the y axis. The options are:

      + `:stellar_mass` -> Cumulative stellar mass.
      + `:sfr`          -> Star formation rate.
  - `smooth::Int=0`: The result will be smoothed out using `smooth` bins. Set it to 0 if you want no smoothing.
  - `output_path::String="."`: Path to the output folder.
  - `sim_labels::Union{Vector{<:AbstractString},Nothing}=basename.(simulation_paths)`: Labels for the plot legend, one per simulation. Set it to `nothing` if you don't want a legend.
  - `theme::Attributes=Theme()`: Plot theme that will take precedence over [`DEFAULT_THEME`](@ref).
"""
function sfrTXT(
    simulation_paths::Vector{String},
    x_quantity::Symbol,
    y_quantity::Symbol;
    smooth::Int=0,
    output_path::String=".",
    sim_labels::Union{Vector{<:AbstractString},Nothing}=basename.(simulation_paths),
    theme::Attributes=Theme(),
)::Nothing

    x_plot_params = plotParams(x_quantity)
    y_plot_params = plotParams(y_quantity)

    if isone(length(simulation_paths))
        filename = "$(basename(simulation_paths[1]))_$(y_quantity)_vs_$(x_quantity)"
    else
        filename = "$(y_quantity)_vs_$(x_quantity)"
    end

    plotTimeSeries(
        simulation_paths,
        [lines!];
        output_path,
        filename,
        da_functions=[daSFRtxt],
        da_args=[(x_quantity, y_quantity)],
        da_kwargs=[(; smooth)],
        x_unit=x_plot_params.unit,
        y_unit=y_plot_params.unit,
        x_exp_factor=x_plot_params.exp_factor,
        y_exp_factor=y_plot_params.exp_factor,
        xaxis_label=x_plot_params.axis_label,
        yaxis_label=y_plot_params.axis_label,
        xaxis_var_name=x_plot_params.var_name,
        yaxis_var_name=y_plot_params.var_name,
        theme=merge(
            theme,
            Theme(
                size=(1400, 880),
                figure_padding=(10, 15, 5, 15),
                palette=(linestyle=[:solid],),
                Axis=(aspect=nothing, xticks=0:14),
                Legend=(halign=:left, valign=:top, nbanks=1),
            ),
        ),
        sim_labels,
    )

    return nothing

end

"""
    cpuTXT(
        simulation_paths::Vector{String},
        target::String,
        x_quantity::Symbol,
        y_quantity::Symbol;
        <keyword arguments>
    )::Nothing

Plot a time series of the data in the `cpu.txt` file.

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`. All the simulations will be plotted together.
  - `target::String`: Target process.
  - `x_quantity::Symbol`: Quantity for the x axis. The options are:

      + `:time_step`              -> Time step.
      + `:physical_time`          -> Physical time since the Big Bang.
      + `:clock_time_s`           -> Clock time duration of the time step in seconds.
      + `:clock_time_percent`     -> Clock time duration of the time step as a percentage.
      + `:tot_clock_time_s`       -> Total clock time in seconds.
      + `:tot_clock_time_percent` -> Total clock time as a percentage.
  - `y_quantity::Symbol`: Quantity for the y axis. The options are:

      + `:time_step`              -> Time step.
      + `:physical_time`          -> Physical time since the Big Bang.
      + `:clock_time_s`           -> Clock time duration of the time step in seconds.
      + `:clock_time_percent`     -> Clock time duration of the time step as a percentage.
      + `:tot_clock_time_s`       -> Total clock time in seconds.
      + `:tot_clock_time_percent` -> Total clock time as a percentage.
  - `smooth::Int=0`: The result will be smoothed out using `smooth` bins. Set it to 0 if you want no smoothing.
  - `ylog::Bool=false`: If the y axis is will have a ``\\log_{10}`` scale.
  - `x_trim::NTuple{2,<:Real}=(-Inf, Inf)`: The data will be trim down so the x coordinates fit within `x_trim`.
  - `y_trim::NTuple{2,<:Real}=(-Inf, Inf)`: The data will be trim down so the y coordinates fit within `y_trim`.
  - `output_path::String="."`: Path to the output folder.
  - `sim_labels::Union{Vector{<:AbstractString},Nothing}=basename.(simulation_paths)`: Labels for the plot legend, one per simulation. Set it to `nothing` if you don't want a legend.
  - `theme::Attributes=Theme()`: Plot theme that will take precedence over [`DEFAULT_THEME`](@ref).
"""
function cpuTXT(
    simulation_paths::Vector{String},
    target::String,
    x_quantity::Symbol,
    y_quantity::Symbol;
    smooth::Int=0,
    ylog::Bool=false,
    x_trim::NTuple{2,<:Real}=(-Inf, Inf),
    y_trim::NTuple{2,<:Real}=(-Inf, Inf),
    output_path::String=".",
    sim_labels::Union{Vector{<:AbstractString},Nothing}=basename.(simulation_paths),
    theme::Attributes=Theme(),
)::Nothing

    x_plot_params = plotParams(x_quantity)
    y_plot_params = plotParams(y_quantity)

    safe_str_target = replace(target, "/" => "-", "_" => "-")

    # Set arguments for the y axis
    if ylog
        y_log        = y_plot_params.unit
        y_unit       = Unitful.NoUnits
        yaxis_label  = L"\log_{10} \, " * getLabel("auto_label", 0, y_plot_params.unit)
    else
        y_log        = nothing
        y_unit       = y_plot_params.unit
        yaxis_label  = y_plot_params.axis_label
    end

    if isone(length(simulation_paths))
        filename = "$(basename(simulation_paths[1]))_$(y_quantity)_vs_$(x_quantity)_for_$(safe_str_target)"
    else
        filename = "$(y_quantity)_vs_$(x_quantity)_for_$(safe_str_target)"
    end

    plotTimeSeries(
        simulation_paths,
        [lines!];
        output_path,
        filename,
        da_functions=[daCPUtxt],
        da_args=[(target, x_quantity, y_quantity)],
        da_kwargs=[(; y_log, smooth)],
        x_unit=x_plot_params.unit,
        y_unit,
        x_trim,
        y_trim,
        xaxis_label=x_plot_params.axis_label,
        yaxis_label,
        xaxis_var_name=x_plot_params.var_name,
        yaxis_var_name=y_plot_params.var_name,
        theme=merge(
            theme,
            Theme(
                size=(1400, 880),
                figure_padding=(10, 15, 5, 15),
                palette=(linestyle=[:solid],),
                Axis=(aspect=nothing,),
                Legend=(halign=:left, valign=:top, nbanks=1),
            ),
        ),
        sim_labels,
        title=L"\mathrm{Process: \,\, %$(safe_str_target)}",
    )

    return nothing

end

"""
    kennicuttSchmidtLaw(
        simulation_paths::Vector{String},
        slice::IndexType;
        <keyword arguments>
    )::Nothing

Plot the Kennicutt-Schmidt law.

!!! note

    Only stars younger than [`AGE_RESOLUTION`](@ref) are considered. The star formation area density is the stellar mass area density divided by [`AGE_RESOLUTION`](@ref).

!!! note

    This function uses physical units regardless of the [`PHYSICAL_UNITS`](@ref) global setting.

!!! note

    For our model, the neutral component includes the molecular, atomic, and stellar fractions. And the molecular component includes the molecular and stellar fractions.

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`. All the simulations will be plotted together.
  - `slice::IndexType`: Slice of the simulations, i.e. which snapshots will be plotted. It can be an integer (a single snapshot), a vector of integers (several snapshots), an `UnitRange` (e.g. 5:13), an `StepRange` (e.g. 5:2:13) or (:) (all snapshots). Starts at 1 and out of bounds indices are ignored. All the selected snapshots will be plotted together.
  - `quantity::Symbol=:molecular`: Quantity for the x axis. The options are:

      + `:gas`       -> Total gas mass area density.
      + `:molecular` -> Molecular mass area density.
      + `:atomic`    -> Atomic mass area density.
      + `:neutral`   -> Neutral mass area density.
  - `gas_type::Symbol=:cells`: If the gas area density will be calculated assuming the gas is made up of `:particles` or Voronoi `:cells`.
  - `reduce_grid::Symbol=:square`: Grid for the density projection. The options are:

      + `:square`   -> The gas and stellar distributions will be projected into a regular cubic grid first and then into a flat square one, to emulate the way the area densities are measured in observations.
      + `:circular` -> The gas and stellar distributions will be projected into a regular cubic grid first, then into a flat square one, and finally into a flat circular grid, formed by a series of concentric rings. This emulates the traditional way the Kennicutt-Schmidt law is measured in simulations.
  - `grid_size::Unitful.Length=BOX_L`: Physical side length of the cubic and square grids (if `reduce_grid` = :square), and diameter of the circular grid (if `reduce_grid` = :circular). This limits which cells/particles will be consider. As a reference, Bigiel et al. (2008) uses measurements up to the optical radius r25 (where the B-band magnitude drops below 25 mag arcsec^−2).
  - `bin_size::Unitful.Length=BIGIEL_PX_SIZE`: Target bin size for the grids. If `reduce_grid` = :square, it is the physical side length of the pixels in the final square grid. If `reduce_grid` = :circular, it is the ring width for the final circular grid. In both cases of `reduce_grid`, the result will only be exact if `bin_size` divides `grid_size` exactly, otherwise `grid_size` will take priority and the final sizes will only approximate `bin_size`. For the cubic grids a default value of 200 pc is always used.
  - `plot_type::Symbol=:scatter`: If the plot will be a `:scatter` plot or a `:heatmap`. Heatmaps will not show legends or several simulations at once. Scatter plots show one mark per pixel, and heatmaps show a 2D histogram for the number of pixel in each bin.
  - `integrated::Bool=false`: If the integrated (one mark per galaxy) or resolved (several marks per galaxy) Kennicutt-Schmidt law will be plotted. `integrated` = true only works with `plot_type` = `:scatter`. The central value is the weighted median and the error bars are the median absolute deviations.
  - `sfr_density::Bool=true`: If the quantity for the y axis will be the SFR area density or, if `sfr_density` = false, the stellar mass area density.
  - `gas_weights::Union{Symbol,Nothing}=nothing`: If `plot_type` = `:scatter`, each point (a bin in the 2D grid) can be weighted by a gas quantity. If `integrated` = true, the median will be computed with these weights. If `integrated` = false, each point will have a color given by the weight. The possible weights are:

      + `:gas_area_density` -> Gas mass area density of each bin. See the documentation for the function [`daDensity2DProjection`](@ref).
      + `:gas_sfr`          -> The total gas SFR of the column associated with each bin. See the documentation for the function [`daGasSFR2DProjection`](@ref).
      + `:gas_metallicity`  -> The total metallicity of the column associated with each bin. See the documentation for the function [`daMetallicity2DProjection`](@ref).
      + `:temperature`      -> The median gas temperature of the column associated with each bin. See the documentation for the function [`daTemperature2DProjection`](@ref).
  - `post_processing::Function=getNothing`: Post processing function. It can only be [`getNothing`](@ref), [`ppBigiel2008!`](@ref), [`ppBigiel2010!`](@ref), [`ppKennicutt1998!`](@ref), [`ppSun2023!`](@ref) or [`ppLeroy2008!`](@ref). The default units will be force into the post processing function.
  - `pp_args::Tuple=()`: Positional arguments for the post processing function.
  - `pp_kwargs::NamedTuple=(;)`: Keyword arguments for the post processing function.
  - `fit::Bool=false`: If the simulation data will be fitted with a power law. The fit will be plotted as a line. This option is only valid if `integrated` = false and `plot_type` = `:scatter`, otherwise it will be ignored.
  - `x_range::Union{NTuple{2,<:Number},Nothing}=nothing`: x axis range for the heatmap grid. If set to `nothing`, the extrema of the x values will be used. Only relevant if `plot_type` = :heatmap.
  - `y_range::Union{NTuple{2,<:Number},Nothing}=nothing`: y axis range for the heatmap grid. If set to `nothing`, the extrema of the y values will be used. Only relevant if `plot_type` = :heatmap.
  - `n_bins::Int=100`: Number of bins per side of the heatmap grid. Only relevant if `plot_type` = `:heatmap`.
  - `colorbar::Bool=false`: If a colorbar will be added.
  - `output_file::String="./kennicutt_schmidt_law.png"`: Path to the output file.
  - `trans_mode::Union{Symbol,Tuple{TranslationType,RotationType,Dict{Symbol,Vector{String}}}}=:all_box`: How to translate and rotate the cells/particles, before filtering with `filter_mode`. For options see [`selectTransformation`](@ref).
  - `filter_mode::Union{Symbol,Tuple{Function,Dict{Symbol,Vector{String}}}}=:all`: Which cells/particles will be selected. For options see [`selectFilter`](@ref).
  - `sim_labels::Union{Vector{<:AbstractString},Nothing}=basename.(simulation_paths)`: Labels for the plot legend, one per simulation. Set it to `nothing` if you don't want a legend.
  - `theme::Attributes=Theme()`: Plot theme that will take precedence over [`DEFAULT_THEME`](@ref).

# References

R. C. Kennicutt (1998). *The Global Schmidt Law in Star-forming Galaxies*. The Astrophysical Journal, **498(2)**, 541-552. [doi:10.1086/305588](https://doi.org/10.1086/305588)

F. Bigiel et al. (2008). *THE STAR FORMATION LAW IN NEARBY GALAXIES ON SUB-KPC SCALES*. The Astrophysical Journal, **136(6)**, 2846. [doi:10.1088/0004-6256/136/6/2846](https://doi.org/10.1088/0004-6256/136/6/2846)

J. Sun et al. (2023). *Star Formation Laws and Efficiencies across 80 Nearby Galaxies*. The Astrophysical Journal Letters, **945(2)**, L19. [doi:10.3847/2041-8213/acbd9c](https://doi.org/10.3847/2041-8213/acbd9c)
"""
function kennicuttSchmidtLaw(
    simulation_paths::Vector{String},
    slice::IndexType;
    quantity::Symbol=:molecular,
    gas_type::Symbol=:cells,
    reduce_grid::Symbol=:square,
    grid_size::Unitful.Length=BOX_L,
    bin_size::Unitful.Length=BIGIEL_PX_SIZE,
    plot_type::Symbol=:scatter,
    integrated::Bool=false,
    sfr_density::Bool=true,
    gas_weights::Union{Symbol,Nothing}=nothing,
    post_processing::Function=getNothing,
    pp_args::Tuple=(),
    pp_kwargs::NamedTuple=(;),
    fit::Bool=false,
    x_range::Union{NTuple{2,<:Number},Nothing}=nothing,
    y_range::Union{NTuple{2,<:Number},Nothing}=nothing,
    n_bins::Int=100,
    colorbar::Bool=false,
    output_file::String="./kennicutt_schmidt_law.png",
    trans_mode::Union{Symbol,Tuple{TranslationType,RotationType,Dict{Symbol,Vector{String}}}}=:all_box,
    filter_mode::Union{Symbol,Tuple{Function,Dict{Symbol,Vector{String}}}}=:all,
    sim_labels::Union{Vector{<:AbstractString},Nothing}=basename.(simulation_paths),
    theme::Attributes=Theme(),
)::Nothing

    # Compute the number of simulations
    ns = length(simulation_paths)

    ################################################################################################
    # Default values
    ################################################################################################

    # Default high resolution voxel side length
    voxel_size = 150.0u"pc"

    # Default units for the gas area density
    Σg_m_unit = u"Msun"
    Σg_l_unit = u"pc"
    Σg_unit   = Σg_m_unit * Σg_l_unit^-2

    # Default units for the stellar/sfr area density
    Σs_m_unit = u"Msun"
    Σs_l_unit = u"kpc"
    Σs_t_unit = u"yr"
    if sfr_density
        Σs_unit = Σs_m_unit * Σs_t_unit^-1 * Σs_l_unit^-2
    else
        Σs_unit = Σs_m_unit * Σs_l_unit^-2
    end

    ################################################################################################
    # Physical units
    ################################################################################################

    # Save the origial value of the global `PHYSICAL_UNITS`
    og_pu_value = PHYSICAL_UNITS

    if !og_pu_value && logging[]

        @warn("kennicuttSchmidtLaw: The global `PHYSICAL_UNITS` is set to false, \
        but Kennicutt-Schmidt law plots must be in physical units, so the global \
        setting will be ignored and default to true just for this function")

    end

    # Kennicutt-Schmidt law plots must be in physical units even for cosmological simulations
    global PHYSICAL_UNITS = true

    ################################################################################################
    # Check arguments
    ################################################################################################

    (
        quantity ∈ [:gas, :molecular, :atomic, :neutral] ||
        throw(ArgumentError("kennicuttSchmidtLaw: `quantity` can only be :gas, :molecular, \
        :atomic or :neutral, but I got :$(quantity)"))
    )

    (
        plot_type ∈ [:scatter, :heatmap] ||
        throw(ArgumentError("kennicuttSchmidtLaw: `plot_type` can only be :scatter or :heatmap, \
        but I got :$(plot_type)"))
    )

    (
        reduce_grid ∈ [:square, :circular] ||
        throw(ArgumentError("kennicuttSchmidtLaw: `reduce_grid` can only be :square or :circular, \
        but I got :$(reduce_grid)"))
    )

    if integrated

        if plot_type == :heatmap

            (
                logging[] &&
                @warn("kennicuttSchmidtLaw: If `integrated` is set to true, `plot_type` = :heatmap \
                will be ignored and default to :scatter")
            )

            plot_type = :scatter

        end

        if fit

            (
                logging[] &&
                @warn("kennicuttSchmidtLaw: `integrated` is set to true, so `fit` = true will be \
                ignored and default to false")
            )

            fit = false

        end

    end

    if bin_size < voxel_size

        (
            logging[] &&
            @warn("kennicuttSchmidtLaw: `reduce_grid` is set to :square and `bin_size` is set to \
            a value lower than $(voxel_size). This is not allowed. `bin_size` will be ignored and \
            default to $(voxel_size)")
        )

        bin_size = voxel_size

    end

    if bin_size > grid_size / 2.0

        (
            logging[] &&
            @warn("kennicuttSchmidtLaw: `bin_size` is set to a value larger than \
            `grid_size` / 2 = $(grid_size / 2.0). This makes no sense. `bin_size` \
            will be ignored and default to `BIGIEL_PX_SIZE` = $(BIGIEL_PX_SIZE)")
        )

        bin_size = BIGIEL_PX_SIZE

    end

    if !isnothing(sim_labels)

        # Compute the number of labels
        nl = length(sim_labels)

        (
            ns == nl ||
            throw(ArgumentError("kennicuttSchmidtLaw: `sim_labels` must have as many elements as \
            `simulation_paths`, but I got length(simulation_paths) = $(ns) \
            != length(sim_labels) = $(nl)"))
        )

    end

    if plot_type == :heatmap

        if !isnothing(gas_weights)

            (
                logging[] &&
                @warn("kennicuttSchmidtLaw: `plot_type` is set to :heatmap, so `gas_weights` = \
                :$(gas_weights) will be ignored and default to nothing")
            )

            gas_weights = nothing

        end

        if fit

            (
                logging[] &&
                @warn("kennicuttSchmidtLaw: `plot_type` is set to :heatmap, so `fit` = true \
                will be ignored and default to false")
            )

            fit = false

        end

        if ns > 1

            (
                logging[] && @warn("kennicuttSchmidtLaw: `plot_type` is set to :heatmap, so only \
                one simulation at a time can be plotted, but I got length(simulation_paths) = \
                $(ns) > 1. `plot_type` = :heatmap will be ignored and default to :scatter")
            )

            plot_type = :scatter

        end

        if reduce_grid == :circular && logging[]

            @warn("kennicuttSchmidtLaw: `plot_type` is set to :heatmap and `reduce_grid` to \
            :circular. Are you sure you want this?")

        end

        if post_processing != getNothing

            (
                logging[] &&
                @warn("kennicuttSchmidtLaw: `plot_type` is set to :heatmap, so `post_processing` \
                will be ignored and default to getNothing")
            )

            post_processing = getNothing

        end

    end

    if colorbar && ((plot_type == :scatter && isnothing(gas_weights)) || integrated)

        (
            logging[] &&
            @warn("kennicuttSchmidtLaw: `colorbar` is set to true, but there is no color range in \
            the plot (either `plot_type` = :scatter and `gas_weights` = nothing or `integrated` = \
            true). `colorbar` = true will be ignored and default to false")
        )

        colorbar = false

    end

    if post_processing ∉ [
        getNothing,
        ppBigiel2008!,
        ppBigiel2010!,
        ppKennicutt1998!,
        ppSun2023!,
        ppLeroy2008!,
    ]

        (
            logging[] &&
            @warn("kennicuttSchmidtLaw: `post_processing` can only be getNothing, ppBigiel2008!, \
            ppBigiel2010!, ppKennicutt1998!, ppSun2023! or ppLeroy2008!, but I got \
            $(post_processing) which will be ignored and default to getNothing")
        )

        post_processing = getNothing

    end

    if !sfr_density && post_processing != getNothing

        (
            logging[] &&
            @warn("kennicuttSchmidtLaw: `post_processing` can only be getNothing when `sfr_density` \
            is false, but I got $(post_processing) which will be ignored and default to getNothing")
        )

        post_processing = getNothing

    end

    ################################################################################################
    # Compute the grids
    ################################################################################################

    # Compute the number of bins for the high resolution grids
    hr_n_bins = round(Int, uconvert(Unitful.NoUnits, grid_size / voxel_size))

    if reduce_grid == :square

        if grid_size != bin_size

            # Compute the number of bins for the low resolution grids
            lr_n_bins = round(Int, uconvert(Unitful.NoUnits, grid_size / bin_size))

            # Compute the interger factor between the high resolution grids (`hr_n_bins`px)
            # and the low resolution grids (`lr_n_bins`px)
            reduce_factor = hr_n_bins ÷ lr_n_bins

        else

            # If the grid size is equal to the bin size, there is no need to reduce the grid
            lr_n_bins     = hr_n_bins
            reduce_factor = 1

        end

        stellar_grid = CubicGrid(grid_size, reduce_factor * lr_n_bins)
        gas_grid     = CubicGrid(grid_size, reduce_factor * lr_n_bins)

    else

        stellar_grid = CubicGrid(grid_size, hr_n_bins)
        gas_grid     = CubicGrid(grid_size, hr_n_bins)

        # Compute the ring width for the circular grid
        reduce_factor = round(Int, uconvert(Unitful.NoUnits, (grid_size / 2.0) / bin_size))

    end

    ################################################################################################
    # Compute the density maps and save them as JLD2 files
    ################################################################################################

    # Set a folder for the JLD2 files
    temp_folder = joinpath(dirname(output_file), "_temp_jld2")

    ##########################
    # Compute the stellar map
    ##########################

    base_request = Dict(:stellar => ["MASS", "POS ", "GAGE"])

    translation, rotation, trans_request = selectTransformation(trans_mode, base_request)
    filter_function, request = selectFilter(filter_mode, trans_request)

    plotSnapshot(
        simulation_paths,
        request,
        [heatmap!];
        output_path=temp_folder,
        base_filename="stellar_mass",
        slice,
        transform_box=true,
        translation,
        rotation,
        filter_function,
        da_functions=[daDensity2DProjection],
        da_args=[(stellar_grid, :stellar, :particles)],
        da_kwargs=[
            (;
                reduce_grid,
                reduce_factor,
                m_unit=Σs_m_unit,
                l_unit=Σs_l_unit,
                filter_function=dd->filterByStellarAge(dd),
            ),
        ],
        x_unit=u"kpc",
        y_unit=u"kpc",
        save_figures=false,
        backup_results=true,
    )

    #############################
    # Compute the `quantity` map
    #############################

    if quantity == :gas

        da_args=[(gas_grid, :gas, gas_type)]

    elseif quantity == :molecular

        da_args = Vector{Tuple}(undef, length(simulation_paths))

        for (i, simulation_path) in pairs(simulation_paths)
            if isSimSFM(simulation_path)
                da_args[i] = (gas_grid, :ode_molecular_stellar, gas_type)
            else
                da_args[i] = (gas_grid, :br_molecular, gas_type)
            end
        end

    elseif quantity == :atomic

        da_args = Vector{Tuple}(undef, length(simulation_paths))

        for (i, simulation_path) in pairs(simulation_paths)
            if isSimSFM(simulation_path)
                da_args[i] = (gas_grid, :ode_atomic, gas_type)
            else
                da_args[i] = (gas_grid, :br_atomic, gas_type)
            end
        end

    elseif quantity == :neutral

        da_args = Vector{Tuple}(undef, length(simulation_paths))

        for (i, simulation_path) in pairs(simulation_paths)
            if isSimSFM(simulation_path)
                da_args[i] = (gas_grid, :ode_neutral, gas_type)
            else
                da_args[i] = (gas_grid, :neutral, gas_type)
            end
        end

    end

    base_request = mergeRequests(plotParams(:mass).request, Dict(:gas => ["POS ", "RHO "]))

    translation, rotation, trans_request = selectTransformation(trans_mode, base_request)
    filter_function, request = selectFilter(filter_mode, trans_request)

    plotSnapshot(
        simulation_paths,
        request,
        [heatmap!];
        output_path=temp_folder,
        base_filename=string(quantity),
        slice,
        transform_box=true,
        translation,
        rotation,
        filter_function,
        da_functions=[daDensity2DProjection],
        da_args,
        da_kwargs=[(; reduce_grid, reduce_factor, m_unit=Σg_m_unit, l_unit=Σg_l_unit)],
        x_unit=u"kpc",
        y_unit=u"kpc",
        save_figures=false,
        backup_results=true,
    )

    ##########################
    # Compute the weights map
    ##########################

    if !isnothing(gas_weights)

        if gas_weights == :gas_area_density

            da_function = daDensity2DProjection
            da_args     = [(gas_grid, :gas, gas_type)]
            m_unit      = Σg_m_unit
            l_unit      = Σg_l_unit
            c_unit      = Σg_unit
            da_kwargs   = [(; reduce_grid, reduce_factor, m_unit, l_unit)]

        elseif gas_weights == :gas_sfr

            da_function = daGasSFR2DProjection
            da_args     = [(gas_grid, gas_type)]
            m_unit      = Σg_m_unit
            t_unit      = u"yr"
            c_unit      = m_unit * t_unit^-1
            da_kwargs   = [(; reduce_grid, reduce_factor, m_unit, t_unit)]

        elseif gas_weights == :gas_metallicity

            da_function = daMetallicity2DProjection
            da_args     = [(gas_grid, :gas, gas_type)]
            c_unit      = Unitful.NoUnits
            da_kwargs   = [(; reduce_grid, reduce_factor)]

        elseif gas_weights == :temperature

            da_function = daTemperature2DProjection
            da_args     = [(gas_grid, gas_type)]
            c_unit      = plotParams(gas_weights).unit
            da_kwargs   = [(; reduce_grid, reduce_factor)]

        else

            throw(ArgumentError("kennicuttSchmidtLaw: `gas_weights` can only be \
            :gas_area_density, :gas_sfr, :gas_metallicity or :temperature, but I got \
            :$(gas_weights)"))

        end

        c_var_name = plotParams(gas_weights).var_name
        c_label    = LaTeXString(L"\log_{10} \, " * getLabel(c_var_name, 0, c_unit))

        base_request = mergeRequests(
            plotParams(gas_weights).request,
            Dict(:gas => ["MASS", "POS ", "RHO "]),
        )

        translation, rotation, trans_request = selectTransformation(trans_mode, base_request)
        filter_function, request = selectFilter(filter_mode, trans_request)

        plotSnapshot(
            simulation_paths,
            request,
            [heatmap!];
            output_path=temp_folder,
            base_filename="gas_weights",
            slice,
            transform_box=true,
            translation,
            rotation,
            filter_function,
            da_functions=[da_function],
            da_args,
            da_kwargs,
            x_unit=u"kpc",
            y_unit=u"kpc",
            save_figures=false,
            backup_results=true,
        )

    end

    ################################################################################################
    # Set the axis labels
    ################################################################################################

    if quantity == :gas

        x_label = getLabel(plotParams(:gas_area_density).var_name, 0, Σg_unit)

    elseif quantity == :molecular

        x_label = getLabel(L"\Sigma_\text{H2}", 0, Σg_unit)

    elseif quantity == :atomic

        x_label = getLabel(L"\Sigma_\text{HI}", 0, Σg_unit)

    elseif quantity == :neutral

        x_label = getLabel(plotParams(:neutral_area_density).var_name, 0, Σg_unit)

    end

    if sfr_density

        y_label = getLabel(plotParams(:sfr_area_density).var_name, 0, Σs_unit)

    else

        y_label = getLabel(plotParams(:stellar_area_density).var_name, 0, Σs_unit)

    end

    ################################################################################################
    # Read and plot the data in the JLD2 files
    ################################################################################################

    if sfr_density
        # Factor to go from stellar area density to SFR area density
        # log10(Σsfr) = log10(Σ*) - log10Δt
        log10Δt = log10(ustrip(Σs_t_unit, AGE_RESOLUTION))
    end

    # Set the plot theme
    if integrated || reduce_grid == :circular
        markersize = 20
    else
        markersize = 6
    end

    current_theme = merge(
        theme,
        Theme(
            size=colorbar ? (880, 720) : (880, 880),
            fontsize=colorbar ? 28 : 32,
            Legend=(
                nbanks=1,
                halign=:left,
                valign=:top,
                margin=(15, 0, 0, 10),
                labelsize=colorbar ? 28 : 32,
                markersize=20,
            ),
            Scatter=(; markersize, colormap=:nipy_spectral),
            Heatmap=(; colormap=:nipy_spectral),
            Colorbar=(; colormap=:nipy_spectral),
        ),
        DEFAULT_THEME,
        theme_latexfonts(),
    )

    with_theme(current_theme) do

        figure = Figure()

        ax = CairoMakie.Axis(
            figure[1, 1];
            xlabel=LaTeXString(L"\log_{10} \, " * x_label),
            ylabel=LaTeXString(L"\log_{10} \, " * y_label),
        )

        colors    = current_theme[:palette][:color][]
        markers   = current_theme[:palette][:marker][]
        fit_color = :gray20
        pp_color  = :darkgoldenrod3

        for (sim_idx, simulation) in pairs(simulation_paths)

            # Make a dataframe for the simulation with the following columns:
            #  - DataFrame index         -> :row_id
            #  - Number in the file name -> :numbers
            #  - Scale factor            -> :scale_factors
            #  - Redshift                -> :redshifts
            #  - Physical time           -> :physical_times
            #  - Lookback time           -> :lookback_times
            #  - Snapshot path           -> :snapshot_paths
            #  - Group catalog path      -> :groupcat_paths
            simulation_table = DataFrame(makeSimulationTable(simulation)[slice, :])

            sim_name         = basename(simulation)
            snapshot_numbers = simulation_table[!, :numbers]

            # For heatmaps we need to accumulate the values for every snapshot before plotting
            if plot_type == :heatmap
                x_heatmap = Float64[]
                y_heatmap = Float64[]
            end

            for snapshot_number in snapshot_numbers

                ############################################
                # Read the JLD2 files and sanitize the data
                ############################################

                ##############
                # Gas density
                ##############

                x_address = "$(SNAP_BASENAME)_$(snapshot_number)/$(sim_name)"
                x_file    = jldopen(joinpath(temp_folder, "$(string(quantity)).jld2"), "r")

                if reduce_grid == :square
                    x_data = vec(x_file[x_address][3])
                else
                    x_data = x_file[x_address][3]
                end

                x_idxs = map(x -> isnan(x) || iszero(x), x_data)

                ##############
                # SFR density
                ##############

                y_address = "$(SNAP_BASENAME)_$(snapshot_number)/$(sim_name)"
                y_file    = jldopen(joinpath(temp_folder, "stellar_mass.jld2"), "r")

                if reduce_grid == :square
                    y_data = vec(y_file[y_address][3])
                else
                    y_data = y_file[y_address][3]
                end

                y_idxs = map(x -> isnan(x) || iszero(x), y_data)

                delete_idxs = x_idxs ∪ y_idxs

                ##########
                # Weights
                ##########

                if !isnothing(gas_weights)

                    z_address = "$(SNAP_BASENAME)_$(snapshot_number)/$(sim_name)"
                    z_file    = jldopen(joinpath(temp_folder, "gas_weights.jld2"), "r")

                    if reduce_grid == :square
                        z_data = vec(z_file[z_address][3])
                    else
                        z_data = z_file[z_address][3]
                    end

                    z_idxs = map(x -> isnan(x) || iszero(x), z_data)

                    delete_idxs = delete_idxs ∪ z_idxs

                    deleteat!(z_data, delete_idxs)

                end

                deleteat!(x_data, delete_idxs)
                deleteat!(y_data, delete_idxs)

                if sfr_density
                    y_data .-= log10Δt
                end

                # For the integrated Kennicutt-Schmidt law, compute the median and median absolute
                # deviation of the gas and stellar densities
                if integrated

                    lin_x = exp10.(x_data)
                    lin_y = exp10.(y_data)

                    if !isnothing(gas_weights)

                        w_z = weights(exp10.(z_data))

                        μx = median(lin_x, w_z)
                        μy = median(lin_y, w_z)

                    else

                        μx = median(lin_x)
                        μy = median(lin_y)

                    end

                    σx = mad(lin_x; center=μx, normalize=false)
                    σy = mad(lin_y; center=μy, normalize=false)

                    x_data = [log10(μx ± σx)]
                    y_data = [log10(μy ± σy)]

                end

                if plot_type == :heatmap

                    append!(x_heatmap, x_data)
                    append!(y_heatmap, y_data)

                end

                close(x_file)
                close(y_file)
                isnothing(gas_weights) || close(z_file)

                #################################
                # Plot the Kennicutt-Schmidt law
                #################################

                if  plot_type == :scatter

                    if integrated

                        x_values = Measurements.value.(x_data)
                        y_values = Measurements.value.(y_data)

                        x_uncertainty = Measurements.uncertainty.(x_data)
                        y_uncertainty = Measurements.uncertainty.(y_data)

                        color = ring(colors, sim_idx)

                        scatter!(ax, x_values, y_values; color)
                        errorbars!(ax, x_values, y_values, x_uncertainty; color, direction=:x)
                        errorbars!(ax, x_values, y_values, y_uncertainty; color, direction=:y)

                    else

                        if isnothing(gas_weights)

                            color = ring(colors, sim_idx)

                        else

                            color = z_data

                        end

                        scatter!(ax, x_data, y_data; color)

                    end

                end

            end

            if plot_type == :heatmap

                # If there is no specified range, use the extrema of the x values
                if isnothing(x_range)
                    xrange = extrema(x_heatmap)
                else
                    xrange = x_range
                end

                # If there is no specified range, use the extrema of the y values
                if isnothing(y_range)
                    yrange = extrema(y_heatmap)
                else
                    yrange = y_range
                end

                # Compute the bin half width for each axis
                x_bin_h_width = 0.5 * (xrange[2] - xrange[1]) / n_bins
                y_bin_h_width = 0.5 * (yrange[2] - yrange[1]) / n_bins

                # Compute the center value of each bin for each axis
                x_axis = collect(
                    range(xrange[1] + x_bin_h_width; length=n_bins, step=2*x_bin_h_width),
                )
                y_axis = collect(
                    range(yrange[1] + y_bin_h_width; length=n_bins, step=2*y_bin_h_width),
                )

                # Compute the 2D histogram (number of pixels in each bin)
                values = histogram2D(
                    permutedims(hcat(x_heatmap, y_heatmap), (2, 1)),
                    collect(range(xrange[1], xrange[2]; length=n_bins + 1)),
                    collect(range(yrange[1], yrange[2]; length=n_bins + 1));
                )

                # The transpose and reverse operation are used to conform to the way heatmap!
                # expect the matrix to be structured
                z_axis = reverse!(transpose(values), dims=2)

                heatmap!(ax, x_axis, y_axis, z_axis)

                if logging[]

                    clean_values = filter(!isnan, z_axis)

                    if isempty(clean_values)

                        min_max_v = (NaN, NaN)
                        mean_v    = NaN
                        median_v  = NaN
                        mode_v    = NaN

                    else

                        min_max_v = extrema(clean_values)
                        mean_v    = mean(clean_values)
                        median_v  = median(clean_values)
                        mode_v    = mode(clean_values)

                    end

                    # Print the count range
                    @info(
                        "\nCount range (KS Law) \
                        \n  Simulation: $(basename(simulation)) \
                        \n  Quantity:   $(quantity) \
                        \n  Min - Max:  $(min_max_v) \
                        \n  Mean:       $(mean_v) \
                        \n  Median:     $(median_v) \
                        \n  Mode:       $(mode_v)"
                    )

                end

            end

        end

        ############################################################################################
        # Apply the post processing function
        ############################################################################################

        scatter_legend = (!isnothing(sim_labels) && plot_type == :scatter)

        if fit
            fit_legend = ppFitLine!(
                figure;
                text_position=(0.98, 0.12),
                text_align=(:right, :bottom),
                color=fit_color,
            )
        end

        if scatter_legend

            if !isnothing(gas_weights) && !integrated

                marker_elements = [
                    MarkerElement(; color=first(colors), marker=first(markers))
                    for _ in eachindex(sim_labels)
                ]

            else

                marker_elements = [
                    MarkerElement(; color=ring(colors, i), marker=ring(markers, i)) for i in 1:ns
                ]

            end

        end

        # Force consistent units and colors
        pp_kwargs = merge(pp_kwargs, (; x_unit=Σg_unit, y_unit=Σs_unit, color=pp_color))
        pp_legend = post_processing(figure, pp_args...; pp_kwargs...)

        if scatter_legend && !isnothing(pp_legend)

            marker_elements = vcat(marker_elements, pp_legend[1])
            sim_labels      = vcat(sim_labels, pp_legend[2])

        end

        if scatter_legend && fit && !isnothing(fit_legend)

            marker_elements = vcat(marker_elements, fit_legend[1])
            sim_labels      = vcat(sim_labels, fit_legend[2])

        end

        if scatter_legend

            Makie.Legend(figure[1, 1], marker_elements, sim_labels)

        end

        if !isnothing(sim_labels) && plot_type == :heatmap

            ppAnnotation!(figure, sim_labels[1]; color=:white, fontsize=30)

        end

        ############################################################################################
        # Print the colorbar
        ############################################################################################

        if colorbar

            if plot_type == :heatmap
                Colorbar(figure[1, 2], figure.content[1].scene.plots[1]; label=L"\mathrm{Counts}")
            else
                Colorbar(figure[1, 2], figure.content[1].scene.plots[1]; label=c_label)
            end

            rowsize!(figure.layout, 1, Auto(1.0))

        end

        ############################################################################################
        # Save the plot
        ############################################################################################

        Makie.save(output_file, figure)

    end

    rm(temp_folder; recursive=true)

    # Restore the original value of `PHYSICAL_UNITS`
    global PHYSICAL_UNITS = og_pu_value

    return nothing

end

"""
    stellarBirthHalos(
        simulation_paths::Vector{String},
        snapshot_n::Int;
        <keyword arguments>
    )::Nothing

Write, to a pair of CSV files, in which halo and subhalo every star in snapshot `snapshot_n` was born.

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`. Each simulation will be in a different file.
  - `snapshot_n::Int`: Selects the target snapshot. Starts at 1 and is independent of the number in the file name. If every snapshot is present, the relation is snapshot_n = (number in filename) + 1.
  - `output_path::String="."`: Path to the output folder.
  - `trans_mode::Union{Symbol,Tuple{TranslationType,RotationType,Dict{Symbol,Vector{String}}}}=:all_box`: How to translate and rotate the cells/particles, before filtering with `filter_mode`. For options see [`selectTransformation`](@ref).
  - `filter_mode::Union{Symbol,Tuple{Function,Dict{Symbol,Vector{String}}}}=:all`: Which cells/particles will be selected. For options see [`selectFilter`](@ref).
"""
function stellarBirthHalos(
    simulation_paths::Vector{String},
    snapshot_n::Int;
    output_path::String=".",
    trans_mode::Union{Symbol,Tuple{TranslationType,RotationType,Dict{Symbol,Vector{String}}}}=:all_box,
    filter_mode::Union{Symbol,Tuple{Function,Dict{Symbol,Vector{String}}}}=:all,
)::Nothing

    # Create the output folder if it doesn't exist
    mkpath(output_path)

    base_request = Dict(:stellar=>["GAGE", "ID  ", "MASS"])

    translation, rotation, trans_request = selectTransformation(trans_mode, base_request)
    filter_function, request = selectFilter(filter_mode, trans_request)

    for simulation_path in simulation_paths

        # Read the relevant data of the snapshot
        data_dict = makeDataDict(simulation_path, snapshot_n, request)

        # Translate the data
        translateData!(data_dict, translation...)

        # Rotate the data
        rotateData!(data_dict, rotation...)

        # Filter the data
        filterData!(data_dict; filter_function)

        # Find the birth time of every star
        birth_time = ustrip.(u"Gyr", computeStellarBirthTime(data_dict))

        if isempty(birth_time)
            (
                logging[] &&
                @info("stellarBirthHalos: There are no stars in snapshot $(snapshot_n) \
                of simulation $(basename(simulation_path)) after applying the filter")
            )
            continue
        end

        # Read the stellar masses
        mass = ustrip.(u"Msun", data_dict[:stellar]["MASS"])

        # Read the stellar IDs
        star_id = data_dict[:stellar]["ID  "]

        # Find the birth place of every star
        birth_halo, birth_subhalo = locateStellarBirthPlace(data_dict)

        if any(isempty, [birth_halo, birth_subhalo])
            (
                logging[] &&
                @info("stellarBirthHalos: I couldn't find the birth halos and/or subhalos for the \
                stars in snapshot $(snapshot_n) of simulation $(basename(simulation_path))")
            )
            birth_halo    = fill(-1, length(masses))
            birth_subhalo = fill(-1, length(masses))
        end

        CSV.write(
            joinpath(output_path, "$(basename(simulation_path))_stellar_birth.csv.gz"),
            DataFrame(; star_id, birth_halo, birth_subhalo, mass, birth_time);
            compress=true,
        )

    end

    return nothing

end

"""
    atomicMolecularTransition(
        simulation_paths::Vector{String},
        slice::IndexType,
        ranges::Vector{<:Tuple{<:Real,<:Real}};
        <keyword arguments>
    )::Nothing

Plot the atomic to molecular gas transition for a set of metallicity ranges.

!!! note

    For our model, the neutral component includes the molecular, atomic, and stellar fractions. And the molecular component includes the molecular and stellar fractions.

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`. Each simulation will be in a different file.
  - `slice::IndexType`: Slice of the simulations, i.e. which snapshots will be plotted. It can be an integer (a single snapshot), a vector of integers (several snapshots), an `UnitRange` (e.g. 5:13), an `StepRange` (e.g. 5:2:13) or (:) (all snapshots). Starts at 1 and out of bounds indices are ignored.
  - `ranges::Vector{<:Tuple{<:Real,<:Real}}`: Metallicity (in solar units) ranges.
  - `plot_type::Symbol=:heatmap`: Type of plot. The options are:

      + `:heatmap` -> Heatmap. One figure per range will be produced.
      + `:scatter` -> Scatter plot. A single figure with every range will be produced.
  - `output_path::String="."`: Path to the output folder.
  - `trans_mode::Union{Symbol,Tuple{TranslationType,RotationType,Dict{Symbol,Vector{String}}}}=:all_box`: How to translate and rotate the cells/particles, before filtering with `filter_mode`. For options see [`selectTransformation`](@ref).
  - `filter_mode::Union{Symbol,Tuple{Function,Dict{Symbol,Vector{String}}}}=:all`: Which cells/particles will be selected. For options see [`selectFilter`](@ref).
  - `theme::Attributes=Theme()`: Plot theme that will take precedence over [`DEFAULT_THEME`](@ref).
"""
function atomicMolecularTransition(
    simulation_paths::Vector{String},
    slice::IndexType,
    ranges::Vector{<:Tuple{<:Real,<:Real}};
    plot_type::Symbol=:heatmap,
    output_path::String=".",
    trans_mode::Union{Symbol,Tuple{TranslationType,RotationType,Dict{Symbol,Vector{String}}}}=:all_box,
    filter_mode::Union{Symbol,Tuple{Function,Dict{Symbol,Vector{String}}}}=:all,
    theme::Attributes=Theme(),
)::Nothing

    for simulation_path in simulation_paths

        # Set some plotting parameters
        if isSimSFM(simulation_path)

            x_quantity = :ode_atomic_number_density
            y_quantity = ratio(:ode_molecular_stellar_mass, :ode_neutral_mass)

        else

            x_quantity = :br_atomic_number_density
            y_quantity = ratio(:br_molecular_mass, :neutral_mass)

        end

        x_plot_params = plotParams(x_quantity)
        y_plot_params = plotParams(y_quantity)

        base_request = mergeRequests(
            x_plot_params.request,
            y_plot_params.request,
            Dict(:gas => ["GZ  "]),
        )

        translation, rotation, trans_request = selectTransformation(trans_mode, base_request)
        filter_function, request = selectFilter(filter_mode, trans_request)

        xaxis_label = L"\log_{10} \, " * getLabel("auto_label", 0, x_plot_params.unit)
        yaxis_label = L"\log_{10} \, " * getLabel("auto_label", 0, y_plot_params.unit)

        filename = "$(basename(simulation_path))_$(y_quantity)_vs_$(x_quantity)"

        if plot_type == :heatmap

            for range in ranges

                plotSnapshot(
                    fill(simulation_path, length(ranges)),
                    request,
                    [heatmap!];
                    output_path,
                    base_filename="$(filename)_$(range[1])_Z_$(range[2])",
                    slice,
                    transform_box=true,
                    translation,
                    rotation,
                    filter_function,
                    da_functions=[daScatterDensity],
                    da_args=[(x_quantity, y_quantity)],
                    da_kwargs=[
                        (;
                            x_log=x_plot_params.unit,
                            y_log=y_plot_params.unit,
                            filter_function=dd -> filterByQuantity(
                                dd,
                                :gas_metallicity,
                                :gas,
                                range[1],
                                range[2],
                            )
                        )
                    ],
                    xaxis_label,
                    yaxis_label,
                    xaxis_var_name=x_plot_params.var_name,
                    yaxis_var_name=y_plot_params.var_name,
                    theme,
                    title=L"%$(range[1]) \, < \, Z_\mathrm{gas} \, [\mathrm{Z_\odot}] \, < \, %$(range[2])",
                )

            end

        elseif plot_type == :scatter

            plotSnapshot(
                fill(simulation_path, length(ranges)),
                request,
                [scatter!];
                pf_kwargs=[(; markersize=4)],
                output_path,
                base_filename=filename,
                slice,
                transform_box=true,
                translation,
                rotation,
                filter_function,
                da_functions=[daScatterGalaxy],
                da_args=[(x_quantity, y_quantity)],
                da_kwargs = [
                    (;
                        x_log=x_plot_params.unit,
                        y_log=y_plot_params.unit,
                        filter_function=dd -> filterByQuantity(
                            dd,
                            :gas_metallicity,
                            :gas,
                            range[1],
                            range[2],
                        ),
                    ) for range in ranges
                ],
                xaxis_label,
                yaxis_label,
                xaxis_var_name=x_plot_params.var_name,
                yaxis_var_name=y_plot_params.var_name,
                theme=merge(
                    theme,
                    Theme(
                        size=(1000, 880),
                        Axis=(aspect=nothing,),
                        Legend=(nbanks=1, margin=(0, 10, 10, 0)),
                    ),
                ),
                sim_labels= [
                    L"%$(range[1]) \, < \, Z_\mathrm{gas} \, [\mathrm{Z_\odot}] \, < \, %$(range[2])"
                    for range in ranges
                ],
            )

        else

            throw(ArgumentError("atomicMolecularTransition: `plot_type` can only be :heatmap or \
            :scatter, but I got :$(plot_type)"))

        end

    end

    return nothing

end

"""
    massProfile(
        simulation_paths::Vector{String},
        slice::IndexType,
        components::Vector{Symbol};
        <keyword arguments>
    )::Nothing

Plot a mass profile.

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`. Each simulation will be plotted in a different figure.
  - `slice::IndexType`: Slice of the simulation, i.e. which snapshots will be plotted. It can be an integer (a single snapshot), a vector of integers (several snapshots), an `UnitRange` (e.g. 5:13), an `StepRange` (e.g. 5:2:13) or (:) (all snapshots). Starts at 1 and out of bounds indices are ignored.
  - `components::Vector{Symbol}`: Target components. It can only be one of the elements of [`COMPONENTS`](@ref). All the components will be plotted together.
  - `cumulative::Bool=false`: If the profile will be accumulated or not.
  - `ylog::Bool=false`: If the y axis is will have a ``\\log_{10}`` scale.
  - `radius::Unitful.Length=DISK_R`: Radius of the profile.
  - `n_bins::Int=100`: Number of bins.
  - `output_path::String="."`: Path to the output folder.
  - `trans_mode::Union{Symbol,Tuple{TranslationType,RotationType,Dict{Symbol,Vector{String}}}}=:all_box`: How to translate and rotate the cells/particles, before filtering with `filter_mode`. For options see [`selectTransformation`](@ref).
  - `filter_mode::Union{Symbol,Tuple{Function,Dict{Symbol,Vector{String}}}}=:all`: Which cells/particles will be selected. For options see [`selectFilter`](@ref).
  - `da_filter::Function=filterNothing`: Filter function to be applied within [`daScatterGalaxy`](@ref) after `trans_mode` and `filter_mode` are applied. See the required signature and examples in `./src/analysis/filters.jl`.
  - `ff_request::Dict{Symbol,Vector{String}}=Dict{Symbol,Vector{String}}()`: Request dictionary for `da_ff`.
  - `title::Union{Symbol,<:AbstractString}=""`: Title for the figure. If left empty, no title is printed. It can also be set to one of the following options:

      + `:physical_time` -> Physical time since the Big Bang.
      + `:lookback_time` -> Physical time left to reach the last snapshot.
      + `:scale_factor`  -> Scale factor (only relevant for cosmological simulations).
      + `:redshift`      -> Redshift (only relevant for cosmological simulations).
  - `theme::Attributes=Theme()`: Plot theme that will take precedence over [`DEFAULT_THEME`](@ref).
"""
function massProfile(
    simulation_paths::Vector{String},
    slice::IndexType,
    components::Vector{Symbol};
    cumulative::Bool=false,
    ylog::Bool=false,
    radius::Unitful.Length=DISK_R,
    n_bins::Int=100,
    output_path::String=".",
    trans_mode::Union{Symbol,Tuple{TranslationType,RotationType,Dict{Symbol,Vector{String}}}}=:all_box,
    filter_mode::Union{Symbol,Tuple{Function,Dict{Symbol,Vector{String}}}}=:all,
    da_filter::Function=filterNothing,
    ff_request::Dict{Symbol,Vector{String}}=Dict{Symbol,Vector{String}}(),
    title::Union{Symbol,<:AbstractString}="",
    theme::Attributes=Theme(),
)::Nothing

    plot_params = plotParams(:mass)
    base_request = mergeRequests(plot_params.request, ff_request)

    translation, rotation, trans_request = selectTransformation(trans_mode, base_request)
    filter_function, request = selectFilter(filter_mode, trans_request)

    grid = CircularGrid(radius, n_bins)

    n_sims = length(simulation_paths)

    # Set arguments for the y axis
    if ylog
        y_unit       = Unitful.NoUnits
        y_log        = plot_params.unit
        y_exp_factor = 0
        yaxis_label  = L"\log_{10} \, " * getLabel("auto_label", 0, plot_params.unit)
    else
        y_unit       = plot_params.unit
        y_log        = nothing
        y_exp_factor = plot_params.exp_factor
        yaxis_label  = plot_params.axis_label
    end

    sim_labels = ["$(component)_mass" for component in components]

    for simulation_path in simulation_paths

        sim_name = basename(simulation_path)

        if isone(n_sims)
            if cumulative
                base_filename = "mass_profiles_cumulative"
            else
                base_filename = "mass_profiles"
            end
        else
            if cumulative
                base_filename = "$(sim_name)_mass_profiles_cumulative"
            else
                base_filename = "$(sim_name)_mass_profiles"
            end
        end

        plotSnapshot(
            fill(simulation_path, length(components)),
            request,
            [lines!];
            output_path,
            base_filename,
            slice,
            transform_box=true,
            translation,
            rotation,
            filter_function,
            da_functions=[daProfile],
            da_args=[(Symbol(component, :_mass), grid) for component in components],
            da_kwargs=[(; y_log, cumulative, filter_function=da_filter)],
            x_unit=u"kpc",
            y_unit,
            y_exp_factor,
            yaxis_label,
            xaxis_var_name=L"r",
            yaxis_var_name=plot_params.var_name,
            theme=merge(
                theme,
                Theme(
                    size=(1400, 880),
                    figure_padding=(10, 15, 5, 15),
                    palette=(linestyle=[:solid],),
                    Axis=(aspect=nothing,),
                    Legend=(halign=:left, valign=:top, nbanks=1, margin=(15, 0, 0, 10)),
                ),
            ),
            sim_labels,
            title,
        )

    end

    return nothing

end

"""
    velocityProfile(
        simulation_paths::Vector{String},
        slice::IndexType,
        velocity::Symbol;
        <keyword arguments>
    )::Nothing

Plot a velocity profile.

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`. All the simulations will be plotted together.
  - `slice::IndexType`: Slice of the simulation, i.e. which snapshots will be plotted. It can be an integer (a single snapshot), a vector of integers (several snapshots), an `UnitRange` (e.g. 5:13), an `StepRange` (e.g. 5:2:13) or (:) (all snapshots). Starts at 1 and out of bounds indices are ignored.
  - `velocity::Symbol`: Target velocity component. The options are:

      + `:stellar_radial_velocity`     -> Component of the stellar velocity in the radial direction (``v_r``).
      + `:stellar_tangential_velocity` -> Component of the stellar velocity in the tangential direction (``v_\\theta``).
      + `:stellar_zstar_velocity`      -> Component of the stellar velocity in the z direction , computed as ``v_z \\, \\mathrm{sign}(z)``.
  - `radius::Unitful.Length=DISK_R`: Radius of the profile.
  - `n_bins::Int=40`: Number of bins.
  - `output_path::String="."`: Path to the output folder.
  - `trans_mode::Union{Symbol,Tuple{TranslationType,RotationType,Dict{Symbol,Vector{String}}}}=:all_box`: How to translate and rotate the cells/particles, before filtering with `filter_mode`. For options see [`selectTransformation`](@ref).
  - `filter_mode::Union{Symbol,Tuple{Function,Dict{Symbol,Vector{String}}}}=:all`: Which cells/particles will be selected. For options see [`selectFilter`](@ref).
  - `da_filter::Function=filterNothing`: Filter function to be applied within [`daScatterGalaxy`](@ref) after `trans_mode` and `filter_mode` are applied. See the required signature and examples in `./src/analysis/filters.jl`.
  - `ff_request::Dict{Symbol,Vector{String}}=Dict{Symbol,Vector{String}}()`: Request dictionary for `da_ff`.
  - `sim_labels::Union{Vector{<:AbstractString},Nothing}=basename.(simulation_paths)`: Labels for the plot legend, one per simulation. Set it to `nothing` if you don't want a legend.
  - `title::Union{Symbol,<:AbstractString}=""`: Title for the figure. If left empty, no title is printed. It can also be set to one of the following options:

      + `:physical_time` -> Physical time since the Big Bang.
      + `:lookback_time` -> Physical time left to reach the last snapshot.
      + `:scale_factor`  -> Scale factor (only relevant for cosmological simulations).
      + `:redshift`      -> Redshift (only relevant for cosmological simulations).
  - `theme::Attributes=Theme()`: Plot theme that will take precedence over [`DEFAULT_THEME`](@ref).
"""
function velocityProfile(
    simulation_paths::Vector{String},
    slice::IndexType,
    velocity::Symbol;
    radius::Unitful.Length=DISK_R,
    n_bins::Int=40,
    output_path::String=".",
    trans_mode::Union{Symbol,Tuple{TranslationType,RotationType,Dict{Symbol,Vector{String}}}}=:all_box,
    filter_mode::Union{Symbol,Tuple{Function,Dict{Symbol,Vector{String}}}}=:all,
    da_filter::Function=filterNothing,
    ff_request::Dict{Symbol,Vector{String}}=Dict{Symbol,Vector{String}}(),
    sim_labels::Union{Vector{<:AbstractString},Nothing}=basename.(simulation_paths),
    title::Union{Symbol,<:AbstractString}="",
    theme::Attributes=Theme(),
)::Nothing

    plot_params = plotParams(velocity)
    base_request = mergeRequests(plot_params.request, ff_request)

    translation, rotation, trans_request = selectTransformation(trans_mode, base_request)
    filter_function, request = selectFilter(filter_mode, trans_request)

    grid = CircularGrid(radius, n_bins)

    if isone(length(simulation_paths))
        base_filename = "$(basename(simulation_paths[1]))_$(velocity)_profile"
    else
        base_filename = "$(velocity)_profile"
    end

    plotSnapshot(
        simulation_paths,
        request,
        [lines!];
        output_path,
        base_filename,
        slice,
        transform_box=true,
        translation,
        rotation,
        filter_function,
        da_functions=[daProfile],
        da_args=[(velocity, grid)],
        da_kwargs=[(; total=false, filter_function=da_filter)],
        x_unit=u"kpc",
        y_unit=plot_params.unit,
        y_exp_factor=plot_params.exp_factor,
        yaxis_label=plot_params.axis_label,
        xaxis_var_name=L"r",
        yaxis_var_name=plot_params.var_name,
        theme=merge(
            theme,
            Theme(
                size=(1400, 880),
                figure_padding=(10, 15, 5, 15),
                palette=(linestyle=[:solid],),
                Axis=(aspect=nothing,),
                Legend=(halign=:left, valign=:top, nbanks=1, margin=(15, 0, 0, 10)),
            ),
        ),
        sim_labels,
        title,
    )

    return nothing

end

"""
    compareFeldmann2020(
        simulation_paths::Vector{String},
        x_component::Symbol,
        y_component::Symbol;
        <keyword arguments>
    )::Nothing

Plot a time series plus the corresponding experimental results from Feldmann (2020).

!!! note

    For our model, the molecular component includes the molecular and stellar fractions.

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`. All the simulations will be plotted together.
  - `x_component::Symbol`: Component for the x axis. The options are:

      + `:stellar`   -> Stellar mass.
      + `:molecular` -> Molecular hydrogen (``\\mathrm{H_2}``) mass.
      + `:atomic`    -> Atomic hydrogen (``\\mathrm{HI}``) mass.
      + `:sfr`       -> Star formation rate of the last `AGE_RESOLUTION`.
  - `y_component::Symbol`: Component for the y axis. The options are:

      + `:stellar`   -> Stellar mass.
      + `:molecular` -> Molecular hydrogen (``\\mathrm{H_2}``) mass.
      + `:atomic`    -> Atomic hydrogen (``\\mathrm{HI}``) mass.
      + `:sfr`       -> Star formation rate of the last `AGE_RESOLUTION`.
  - `slice::IndexType=(:)`: Slice of the simulation, i.e. which snapshots will be plotted. It can be an integer (a single snapshot), a vector of integers (several snapshots), an `UnitRange` (e.g. 5:13), an `StepRange` (e.g. 5:2:13) or (:) (all snapshots). Starts at 1 and out of bounds indices are ignored.
  - `scatter::Bool=false`: If the data will be presented as a line plot with error bands (default), or a scatter plot.
  - `xlog::Bool=true`: If true, sets the x axis to ``\\log_{10}``(`x_quantity`).
  - `ylog::Bool=true`: If true, sets the y axis to ``\\log_{10}``(`y_quantity`).
  - `output_path::String="."`: Path to the output folder.
  - `trans_mode::Union{Symbol,Tuple{TranslationType,RotationType,Dict{Symbol,Vector{String}}}}=:all_box`: How to translate and rotate the cells/particles, before filtering with `filter_mode`. For options see [`selectTransformation`](@ref).
  - `filter_mode::Union{Symbol,Tuple{Function,Dict{Symbol,Vector{String}}}}=:all`: Which cells/particles will be selected. For options see [`selectFilter`](@ref).
  - `da_ff::Function=filterNothing`: Filter function to be applied within [`daScatterGalaxy`](@ref) after `trans_mode` and `filter_mode` are applied. See the required signature and examples in `./src/analysis/filters.jl`.
  - `ff_request::Dict{Symbol,Vector{String}}=Dict{Symbol,Vector{String}}()`: Request dictionary for `da_ff`.
  - `sim_labels::Union{Vector{<:AbstractString},Nothing}=basename.(simulation_paths)`: Labels for the plot legend, one per simulation. Set it to `nothing` if you don't want a legend.
  - `theme::Attributes=Theme()`: Plot theme that will take precedence over [`DEFAULT_THEME`](@ref).

# References

R. Feldmann (2020). *The link between star formation and gas in nearby galaxies*. Communications Physics **3(226)**. [doi:10.1038/s42005-020-00493-0](https://doi.org/10.1038/s42005-020-00493-0)
"""
function compareFeldmann2020(
    simulation_paths::Vector{String},
    x_component::Symbol,
    y_component::Symbol;
    slice::IndexType=(:),
    scatter::Bool=false,
    xlog::Bool=true,
    ylog::Bool=true,
    output_path::String=".",
    trans_mode::Union{Symbol,Tuple{TranslationType,RotationType,Dict{Symbol,Vector{String}}}}=:all_box,
    filter_mode::Union{Symbol,Tuple{Function,Dict{Symbol,Vector{String}}}}=:all,
    da_ff::Function=filterNothing,
    ff_request::Dict{Symbol,Vector{String}}=Dict{Symbol,Vector{String}}(),
    sim_labels::Union{Vector{<:AbstractString},Nothing}=basename.(simulation_paths),
    theme::Attributes=Theme(),
)::Nothing

    n_sims = length(simulation_paths)

    if x_component == :sfr

        x_quantities   = fill(:observational_sfr, n_sims)
        x_plot_params  = plotParams(:observational_sfr)
        xunit          = x_plot_params.unit
        x_axis_label   = x_plot_params.axis_label
        xaxis_var_name = x_plot_params.var_name

    elseif x_component == :stellar

        x_quantities   = fill(:stellar_mass, n_sims)
        x_plot_params  = plotParams(:stellar_mass)
        xunit          = x_plot_params.unit
        x_axis_label   = x_plot_params.axis_label
        xaxis_var_name = x_plot_params.var_name

    elseif x_component  == :molecular

        x_quantities = Vector{Symbol}(undef, n_sims)

        for (i, simulation_path) in enumerate(simulation_paths)
            if isSimSFM(simulation_path)
                x_quantities[i] = :ode_molecular_stellar_mass
            else
                x_quantities[i] = :br_molecular_mass
            end
        end

        xunit          = u"Msun"
        x_axis_label   = "auto_label"
        c_label        = "\\text{H2}"
        xaxis_var_name = L"M_%$(c_label)"

    elseif x_component == :atomic

        x_quantities = Vector{Symbol}(undef, n_sims)

        for (i, simulation_path) in enumerate(simulation_paths)
            if isSimSFM(simulation_path)
                x_quantities[i] = :ode_atomic_mass
            else
                x_quantities[i] = :br_atomic_mass
            end
        end

        xunit          = u"Msun"
        x_axis_label   = "auto_label"
        c_label        = "\\text{HI}"
        xaxis_var_name = L"M_%$(c_label)"

    else

        throw(ArgumentError("compareFeldmann2020!: `x_component` can only be :stellar, \
        :molecular, :atomic or :sfr, but I got :$(x_component)"))

    end

    if y_component == :sfr

        y_quantities   = fill(:observational_sfr, n_sims)
        y_plot_params  = plotParams(:observational_sfr)
        yunit          = y_plot_params.unit
        y_axis_label   = yplot_params.axis_label
        yaxis_var_name = y_plot_params.var_name

    elseif y_component == :stellar

        y_quantities   = fill(:stellar_mass, n_sims)
        y_plot_params  = plotParams(:stellar_mass)
        yunit          = y_plot_params.unit
        y_axis_label   = y_plot_params.axis_label
        yaxis_var_name = y_plot_params.var_name

    elseif y_component  == :molecular

        y_quantities = Vector{Symbol}(undef, n_sims)

        for (i, simulation_path) in enumerate(simulation_paths)
            if isSimSFM(simulation_path)
                y_quantities[i] = :ode_molecular_stellar_mass
            else
                y_quantities[i] = :br_molecular_mass
            end
        end

        yunit          = u"Msun"
        y_axis_label   = "auto_label"
        c_label        = "\\text{H2}"
        yaxis_var_name = L"M_%$(c_label)"

    elseif y_component == :atomic

        y_quantities = Vector{Symbol}(undef, n_sims)

        for (i, simulation_path) in enumerate(simulation_paths)
            if isSimSFM(simulation_path)
                y_quantities[i] = :ode_atomic_mass
            else
                y_quantities[i] = :br_atomic_mass
            end
        end

        yunit          = u"Msun"
        y_axis_label   = "auto_label"
        c_label        = "\\text{HI}"
        yaxis_var_name = L"M_%$(c_label)"

    else

        throw(ArgumentError("compareFeldmann2020!: `y_component` can only be :stellar, \
        :molecular, :atomic or :sfr, but I got :$(y_component)"))

    end

    if xlog
        x_log       = xunit
        xaxis_label = L"\log_{10} \, " * getLabel("auto_label", 0, xunit)
    else
        x_log       = nothing
        xaxis_label = x_axis_label
    end

    if ylog
        y_log       = yunit
        yaxis_label = L"\log_{10} \, " * getLabel("auto_label", 0, yunit)
    else
        y_log       = nothing
        yaxis_label = y_axis_label
    end

    if isone(n_sims)
        filename = "$(scatter)_$(ylog)_$(xlog)_$(basename(simulation_paths[1]))_$(x_component)_vs_$(y_component)_Feldmann2020"
    else
        filename = "$(x_component)_vs_$(y_component)_Feldmann2020"
    end

    plotTimeSeries(
        simulation_paths,
        [scatter!];
        output_path,
        filename,
        slice,
        da_functions=[daEvolution],
        da_args=[(x_qty, y_qty) for (x_qty, y_qty) in zip(x_quantities, y_quantities)],
        da_kwargs=[(; trans_mode, filter_mode, extra_filter=da_ff, ff_request, x_log, y_log)],
        post_processing=ppFeldmann2020!,
        pp_args=(x_component, y_component),
        pp_kwargs=(; scatter, xlog, ylog),
        x_unit=xlog ? Unitful.NoUnits : xunit,
        y_unit=ylog ? Unitful.NoUnits : yunit,
        xaxis_label,
        yaxis_label,
        xaxis_var_name,
        yaxis_var_name,
        theme=merge(
            theme,
            Theme(
                size=(1400, 880),
                figure_padding=(10, 15, 5, 15),
                Axis=(aspect=nothing,),
                Legend=(halign=:left, valign=:top, nbanks=1, margin=(10, 0, 0, 10)),
            ),
        ),
        sim_labels,
    )

    return nothing

end

"""
    massMetallicityRelation(
        simulation_paths::Vector{String},
        slice::IndexType;
        <keyword arguments>
    )::Nothing

Plot the resolved mass-metallicity relation. This method plots the M-Z relation at a fix moment in time.

!!! note

    Only stars younger than [`AGE_RESOLUTION`](@ref) are considered.

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`. All the simulations will be plotted together.
  - `slice::IndexType=(:)`: Slice of the simulation, i.e. which snapshots will be plotted. It can be an integer (a single snapshot), a vector of integers (several snapshots), an `UnitRange` (e.g. 5:13), an `StepRange` (e.g. 5:2:13) or (:) (all snapshots). Starts at 1 and out of bounds indices are ignored.
  - `element::Symbol=:all`: Which metals to consider. The options are:

      + `:all` -> Total metallicity in solar units.
      + `:X`   -> Abundance of element ``\\mathrm{X}``, as ``12 + \\log_{10}(\\mathrm{X \\, / \\, H})``. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
  - `mass::Bool=true`: If the x axis will be the stellar mass density (default) or the SFR density.
  - `reduce_factor::Int=1`: Factor by which the resolution of the result will be reduced. This will be applied after the density projection, averaging the value of neighboring pixels. It has to divide the size of `grid` exactly.
  - `output_path::String="."`: Path to the output folder.
  - `trans_mode::Union{Symbol,Tuple{TranslationType,RotationType,Dict{Symbol,Vector{String}}}}=:all_box`: How to translate and rotate the cells/particles, before filtering with `filter_mode`. For options see [`selectTransformation`](@ref).
  - `filter_mode::Union{Symbol,Tuple{Function,Dict{Symbol,Vector{String}}}}=:all`: Which cells/particles will be selected. For options see [`selectFilter`](@ref).
  - `da_ff::Function=filterNothing`: Filter function to be applied within [`daScatterGalaxy`](@ref) after `trans_mode` and `filter_mode` are applied. See the required signature and examples in `./src/analysis/filters.jl`.
  - `ff_request::Dict{Symbol,Vector{String}}=Dict{Symbol,Vector{String}}()`: Request dictionary for `da_ff`.
  - `theme::Attributes=Theme()`: Plot theme that will take precedence over [`DEFAULT_THEME`](@ref).
"""
function massMetallicityRelation(
    simulation_paths::Vector{String};
    slice::IndexType=(:),
    element::Symbol=:all,
    mass::Bool=true,
    reduce_factor::Int=1,
    output_path::String=".",
    trans_mode::Union{Symbol,Tuple{TranslationType,RotationType,Dict{Symbol,Vector{String}}}}=:all_box,
    filter_mode::Union{Symbol,Tuple{Function,Dict{Symbol,Vector{String}}}}=:all,
    theme::Attributes=Theme(),
)::Nothing

    grid = CubicGrid(BOX_L, 400)

    (
        element ∈ [:all, keys(ELEMENT_INDEX)...] ||
        throw(ArgumentError("massMetallicityRelation: `quantity` can only be :all or any \
        of the keys of `ELEMENT_INDEX` (see `./src/constants/globals.jl`), but I got :$(quantity)"))
    )

    temp_folder = joinpath(output_path, "_temp_jld2")

    ################################################################################################
    # Stellar mass density
    ################################################################################################

    base_request = mergeRequests(
        plotParams(:stellar_mass).request,
        Dict(:stellar => ["GAGE", "POS "]),
    )

    translation, rotation, trans_request = selectTransformation(trans_mode, base_request)
    filter_function, request = selectFilter(filter_mode, trans_request)

    plotSnapshot(
        simulation_paths,
        request,
        [heatmap!];
        output_path=temp_folder,
        base_filename="stellar_mass",
        slice,
        transform_box=true,
        translation,
        rotation,
        filter_function,
        da_functions=[daDensity2DProjection],
        da_args=[(grid, :stellar, :particles)],
        da_kwargs=[(; reduce_factor, filter_function=dd->filterByStellarAge(dd))],
        x_unit=u"kpc",
        y_unit=u"kpc",
        save_figures=false,
        backup_results=true,
    )

    ################################################################################################
    # Gas metallicity
    ################################################################################################

    if element == :all
        metal_request = plotParams(:Z_gas_mass).request
    else
        metal_request = plotParams(Symbol(element, "_gas_abundance")).request
    end

    base_request = mergeRequests(metal_request, Dict(:gas => ["RHO "]))

    translation, rotation, trans_request = selectTransformation(trans_mode, base_request)
    filter_function, request = selectFilter(filter_mode, trans_request)

    plotSnapshot(
        simulation_paths,
        request,
        [heatmap!];
        output_path=temp_folder,
        base_filename="gas_metallicity",
        slice,
        transform_box=true,
        translation,
        rotation,
        filter_function,
        da_functions=[daMetallicity2DProjection],
        da_args=[(grid, :gas, :cells)],
        da_kwargs=[(; element, reduce_factor)],
        x_unit=u"kpc",
        y_unit=u"kpc",
        save_figures=false,
        backup_results=true,
    )

    ################################################################################################
    # Plotting
    ################################################################################################

    # Set the x label
    if mass
        x_plot_param = plotParams(:stellar_area_density)
    else
        x_plot_param = plotParams(:sfr_area_density)
    end

    xlabel = LaTeXString(L"\log_{10} \, " * getLabel(x_plot_param.var_name, 0, x_plot_param.unit))

    # Set the y label
    if element == :all
        ylabel = L"$\log_{10}$ %$(plotParams(:gas_metallicity).var_name)"
    else
        ylabel = plotParams(Symbol(element, "_gas_abundance")).axis_label
    end

    if !mass
        # Factor to go from stellar area density to SFR area density
        # log10(Σsfr) = log10(Σ*) - log10Δt
        log10Δt = log10(ustrip(u"yr", AGE_RESOLUTION))
    end

    with_theme(merge(theme, DEFAULT_THEME, theme_latexfonts())) do

        for (sim_idx, simulation) in pairs(simulation_paths)

            # Make a dataframe for the simulation with the following columns:
            #  - DataFrame index         -> :row_id
            #  - Number in the file name -> :numbers
            #  - Scale factor            -> :scale_factors
            #  - Redshift                -> :redshifts
            #  - Physical time           -> :physical_times
            #  - Lookback time           -> :lookback_times
            #  - Snapshot path           -> :snapshot_paths
            #  - Group catalog path      -> :groupcat_paths
            simulation_table = DataFrame(makeSimulationTable(simulation)[slice, :])

            sim_name         = basename(simulation)
            snapshot_numbers = simulation_table[!, :numbers]

            path = mkpath(joinpath(output_path, sim_name))

            for snapshot_number in snapshot_numbers

                f = Figure()

                ax = CairoMakie.Axis(f[1, 1]; xlabel, ylabel)

                x_address = "$(SNAP_BASENAME)_$(snapshot_number)/$(sim_name)"
                y_address = "$(SNAP_BASENAME)_$(snapshot_number)/$(sim_name)"

                jldopen(joinpath(temp_folder, "stellar_mass.jld2"), "r") do x_file

                    jldopen(joinpath(temp_folder, "gas_metallicity.jld2"), "r") do y_file

                        # Read the JLD2 files
                        x_data = vec(x_file[x_address][3])
                        y_data = vec(y_file[y_address][3])

                        # Delete 0s and NaNs in the data vectors
                        x_idxs = map(x -> isnan(x) || iszero(x), x_data)
                        y_idxs = map(x -> isnan(x) || iszero(x), y_data)

                        deleteat!(x_data, x_idxs ∪ y_idxs)
                        deleteat!(y_data, x_idxs ∪ y_idxs)

                        if mass
                            x_axis = x_data
                        else
                            x_axis = x_data .- log10Δt
                        end

                        scatter!(ax, x_axis, y_data; markersize=6, color=WONG_RED)

                    end

                end

                Makie.save(joinpath(path, "$(snapshot_number).png"), f)

            end

        end

    end

    rm(temp_folder; recursive=true)

    return nothing

end

"""
    gasVelocityCubes(
        simulation_paths::Vector{String},
        slice::ReducedIndexType;
        <keyword arguments>
    )::Nothing

Create a HDF5 file with the position, gas mass, velocity, and velocity dispersion of each voxel in a rectangular 3D grid.

The metadata for each snapshot in the HDF5 file includes the physical time in Gyr, the scale factor, and the redshift of that snapshot.

By default, the grid is centered at coordinates (0, 0, 0), has 300x300x300 voxels, and has a side length of [`BOX_L`](@ref). There are as many rows as there are voxels (27000000 by default).

The quantities in the HDF5 file for each voxel are:

Column 01: x coordinate [`l_unit`]
Column 02: y coordinate [`l_unit`]
Column 03: z coordinate [`l_unit`]
Column 04: Ionized hydrogen mass [`m_unit`]
Column 05: Atomic hydrogen mass [`m_unit`]
Column 06: Molecular hydrogen + stars mass [`m_unit`]
Column 07: Metals mass [`m_unit`]
Column 08: Dust mass [`m_unit`]
Column 09: Gas velocity in the x direction [`v_unit`]
Column 10: Gas velocity in the y direction [`v_unit`]
Column 11: Gas velocity in the z direction [`v_unit`]
Column 12: Gas velocity dispersion in the x direction [`v_unit`]
Column 13: Gas velocity dispersion in the y direction [`v_unit`]
Column 14: Gas velocity dispersion in the z direction [`v_unit`]

For simulations with the gas represented by Voronoi cells (e.g. Arepo):

The mass is the mass of cold, atomic or ionized gas intersecting the voxel, so it only considers the cell that is closest to the center of the voxel. The velocity is given by the weighted mean of the velocities of the `n_neighbors` nearest cells. And the velocity dispersion, by the weighted standard deviation.

Notice that for Voronoi cells, the mass will be sample at a high sub-cell resolution (as long as voxel size < cell size), while the velocities are sample at a locally lower resolution (as long as `n_neighbors` > 1). The weights are given by the distance (in kpc) to each neighbor, using a Gaussian kernel.

For simulations with the gas represented by particles (e.g. SPH codes):

The mass is the accumulated mass of the particles within each voxel. The velocity is the mean of the velocities of those particles, and the velocity dispersion is the corresponding standard deviation.

If there are no particles, the mass is 0, and the velocity and velocity dispersion are set to NaN. If there is only one particle, the mass and velocity are the ones from that particle, and the velocity dispersion is set to NaN.

By default (`trans_mode` = :all_box and `filter_mode` = :all) we use the following reference system:

  - The origin is the global center of mass of the simulation.
  - The x, y, and z axis form a right-handed cartesian reference system (x × y = z), where the z axis has the orientation of the global angular momentum, and the x and y axis are roughly in the direction of the global principal axis.

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`. All the simulations will be written to the same file.
  - `slice::ReducedIndexType`: Slice of the simulations, i.e. which snapshots will be plotted. It can be an integer (a single snapshot), a vector of integers (several snapshots), an `UnitRange` (e.g. 5:13) or an `StepRange` (e.g. 5:2:13). Starts at 1 and out of bounds indices are ignored.
  - `gas_type::Symbol=:cells`: If the gas density will be calculated assuming the gas is in `:particles` or in Voronoi `:cells`.
  - `n_neighbors::Int=32`: Number of neighbors for the mean and standard deviation of the velocity. Setting this value to 1 maximizes the resolution for the velocity, and sets the standard deviation (columns 8, 9, and 10) to NaN. This is only relevant for simulations where gas is represented by Voronoi cells (`type` = :cells).
  - `grid::CubicGrid=CubicGrid(BOX_L, 300)`: Cubic grid.
  - `row_major_order::Bool=true`: Store the results in row-major order (as used in C and Python) instead of column-major order (used in Julia, Fortran, and MATLAB). See [Row- and column-major order](https://en.wikipedia.org/wiki/Row-_and_column-major_order).
  - `m_unit::Unitful.Units=u"Msun"`: Mass unit
  - `l_unit::Unitful.Units=u"kpc"`: Length unit.
  - `v_unit::Unitful.Units=u"km * s^-1"`: Velocity unit.
  - `output_file::String="./gas_velocity_cube.hdf5"`: Path to the output file.
  - `trans_mode::Union{Symbol,Tuple{TranslationType,RotationType,Dict{Symbol,Vector{String}}}}=:all_box`: How to translate and rotate the cells/particles, before filtering with `filter_mode`. For options see [`selectTransformation`](@ref).
  - `filter_mode::Union{Symbol,Tuple{Function,Dict{Symbol,Vector{String}}}}=:all`: Which cells/particles will be selected. For options see [`selectFilter`](@ref).
  - `show_progress::Bool=true`: If a progress bar will be shown.
"""
function gasVelocityCubes(
    simulation_paths::Vector{String},
    slice::ReducedIndexType;
    gas_type::Symbol=:cells,
    n_neighbors::Int=32,
    grid::CubicGrid=CubicGrid(BOX_L, 300),
    row_major_order::Bool=true,
    m_unit::Unitful.Units=u"Msun",
    l_unit::Unitful.Units=u"kpc",
    v_unit::Unitful.Units=u"km * s^-1",
    output_file::String="./gas_velocity_cube.hdf5",
    trans_mode::Union{Symbol,Tuple{TranslationType,RotationType,Dict{Symbol,Vector{String}}}}=:all_box,
    filter_mode::Union{Symbol,Tuple{Function,Dict{Symbol,Vector{String}}}}=:all,
    show_progress::Bool=true,
)::Nothing

    # Create the output folder
    mkpath(dirname(output_file))

    # Create the output HDF5 file
    hdf5_file = h5open(output_file, "w")

    # Set the number of columns and rows
    n_rows = grid.n_bins^3
    n_cols = 14

    base_request = mergeRequests(plotParams(:mass).request, Dict(:gas => ["POS ", "VEL "]))

    translation, rotation, trans_request = selectTransformation(trans_mode, base_request)
    filter_function, request = selectFilter(filter_mode, trans_request)

    # For gas cells, reshape the grid to conform to the way `knn` expect the matrix to be structured
    if gas_type == :cells

        physical_grid = Matrix{Float64}(undef, 3, n_rows)

        for i in eachindex(grid.grid)
            physical_grid[1, i] = ustrip(l_unit, grid.grid[i][1])
            physical_grid[2, i] = ustrip(l_unit, grid.grid[i][2])
            physical_grid[3, i] = ustrip(l_unit, grid.grid[i][3])
        end

    end

    for simulation_path in simulation_paths

        simulation_name = basename(simulation_path)

        prog_bar = Progress(
            length(slice),
            dt=0.5,
            desc="Writing the velocity cube for simulation $(simulation_name)... ",
            color=:blue,
            barglyphs=BarGlyphs("|#  |"),
            enabled=show_progress,
        )

        # Create an HDF5 group for each simulation
        hdf5_group = create_group(hdf5_file, simulation_name)

        for snap_n in slice

            data_dict = makeDataDict(simulation_path, snap_n, request)

            snapshot_number = lpad(string(data_dict[:snap_data].global_index), 3, "0")

            # Translate the data
            translateData!(data_dict, translation...)

            # Rotate the data
            rotateData!(data_dict, rotation...)

            # Filter the data
            filterData!(data_dict; filter_function)

            # Load the gas quantities
            gd = data_dict[:gas]

            # Load the cell/particle positions
            positions = gd["POS "]

            # Load the cell/particle velocities
            velocities = ustrip.(v_unit, gd["VEL "])

            # Compute the mass of cold, atomic, and ionized gas in each cell
            if isSimSFM(simulation_path)


                ion_masses  = scatterQty(data_dict, :ode_ionized_mass)
                ato_masses  = scatterQty(data_dict, :ode_atomic_mass)
                mol_masses  = scatterQty(data_dict, :ode_molecular_stellar_mass)
                Z_masses    = scatterQty(data_dict, :ode_metals_mass)
                dust_masses = scatterQty(data_dict, :ode_dust_mass)

            else

                ion_masses  = scatterQty(data_dict, :ionized_mass)
                ato_masses  = scatterQty(data_dict, :br_atomic_mass)
                mol_masses  = scatterQty(data_dict, :br_molecular_mass)
                Z_masses    = scatterQty(data_dict, :Z_gas_mass)
                dust_masses = zeros(eltype(Z_masses), length(Z_masses))

            end

            if any(isempty, [mol_masses, ato_masses, ion_masses, Z_masses, dust_masses, velocities, positions])
                throw(ArgumentError("gasVelocityCubes: Some data is missing (there appears to be \
                no gas in the snapshot), so I cannot construct the velocity cube"))
            end

            # Alocate memory for:
            # Column 01: x coordinate [l_unit]
            # Column 02: y coordinate [l_unit]
            # Column 03: z coordinate [l_unit]
            # Column 04: Ionized hydrogen mass [m_unit]
            # Column 05: Atomic hydrogen mass [m_unit]
            # Column 06: Molecular hydrogen + stars mass [m_unit]
            # Column 07: Metals mass [m_unit]
            # Column 08: Dust mass [m_unit]
            # Column 09: Gas velocity in the x direction [v_unit]
            # Column 10: Gas velocity in the y direction [v_unit]
            # Column 11: Gas velocity in the z direction [v_unit]
            # Column 12: Gas velocity dispersion in the x direction [v_unit]
            # Column 13: Gas velocity dispersion in the y direction [v_unit]
            # Column 14: Gas velocity dispersion in the z direction [v_unit]
            data_matrix = Matrix{Float64}(undef, n_rows, n_cols)

            if gas_type == :cells

                # Compute the volume of each cell
                cell_volumes = gd["MASS"] ./ gd["RHO "]

                # Compute the gas densities
                ion_densities  = ustrip.(m_unit * l_unit^-3, ion_masses ./ cell_volumes)
                ato_densities  = ustrip.(m_unit * l_unit^-3, ato_masses ./ cell_volumes)
                mol_densities  = ustrip.(m_unit * l_unit^-3, mol_masses ./ cell_volumes)
                Z_densities    = ustrip.(m_unit * l_unit^-3, Z_masses ./ cell_volumes)
                dust_densities = ustrip.(m_unit * l_unit^-3, dust_masses ./ cell_volumes)

                # Load the volume of the voxels
                voxel_volume = ustrip(l_unit^3, grid.bin_volume)

                # Compute the tree for a nearest neighbor search
                kdtree = KDTree(ustrip.(l_unit, positions))

                # Find the `n_neighbors` nearest cells to each voxel
                idxs, dists = knn(kdtree, physical_grid, n_neighbors, true)

                Threads.@threads for i in eachindex(grid.grid)

                    # Physical coordinates of the voxel [l_unit]
                    data_matrix[i, 1:3] .= ustrip.(l_unit, grid.grid[i])

                    # Ionized hydrogen mass [m_unit]
                    data_matrix[i, 4] = ion_densities[idxs[i][1]] * voxel_volume
                    # Atomic hydrogen mass [m_unit]
                    data_matrix[i, 5] = ato_densities[idxs[i][1]] * voxel_volume
                    # Molecular hydrogen + stars mass [m_unit]
                    data_matrix[i, 6] = mol_densities[idxs[i][1]] * voxel_volume
                    # Metals mass [m_unit]
                    data_matrix[i, 7] = Z_densities[idxs[i][1]] * voxel_volume
                    # Dust mass [m_unit]
                    data_matrix[i, 8] = dust_densities[idxs[i][1]] * voxel_volume

                    if isone(n_neighbors)

                        # Neighbor velocity in the x direction [v_unit]
                        data_matrix[i, 9]  = velocities[1, idxs[i]]
                        # Neighbor velocity in the y direction [v_unit]
                        data_matrix[i, 10] = velocities[2, idxs[i]]
                        # Neighbor velocity in the z direction [v_unit]
                        data_matrix[i, 11] = velocities[3, idxs[i]]

                        # For the case of only one neighbor, set the standard deviations to NaN
                        data_matrix[i, 12] = NaN
                        data_matrix[i, 13] = NaN
                        data_matrix[i, 14] = NaN

                    else

                        # Compute the analytic weights using a Gaussian kernel
                        neighbor_weights = aweights(evaluateNormal(dists[i]))

                        # Neighbor velocities in the x direction [v_unit]
                        vxs = velocities[1, idxs[i]]
                        # Neighbor velocities in the y direction [v_unit]
                        vys = velocities[2, idxs[i]]
                        # Neighbor velocities in the z direction [v_unit]
                        vzs = velocities[3, idxs[i]]

                        # Mean and standard deviation of the neighbor velocities in the x direction [v_unit]
                        data_matrix[i, 9], data_matrix[i, 12]  = mean_and_std(vxs, neighbor_weights)
                        # Mean and standard deviation of the neighbor velocities in the y direction [v_unit]
                        data_matrix[i, 10], data_matrix[i, 13] = mean_and_std(vys, neighbor_weights)
                        # Mean and standard deviation of the neighbor velocities in the z direction [v_unit]
                        data_matrix[i, 11], data_matrix[i, 14] = mean_and_std(vzs, neighbor_weights)

                    end

                end

            elseif gas_type == :particles

                # Find which particles are within each voxel
                idxs = listHistogram3D(positions, grid)

                Threads.@threads for i in eachindex(grid.grid)

                    # Physical coordinates of the voxel [l_unit]
                    data_matrix[i, 1:3] .= ustrip.(l_unit, grid.grid[i])

                    # Ionized hydrogen mass [m_unit]
                    data_matrix[i, 4] = ustrip(m_unit, sum(ion_masses[idxs[i]]; init=0.0*m_unit))
                    # Atomic hydrogen mass [m_unit]
                    data_matrix[i, 5] = ustrip(m_unit, sum(ato_masses[idxs[i]]; init=0.0*m_unit))
                    # Molecular hydrogen + stars mass [m_unit]
                    data_matrix[i, 6] = ustrip(m_unit, sum(mol_masses[idxs[i]]; init=0.0*m_unit))
                    # Metals mass [m_unit]
                    data_matrix[i, 7] = ustrip(m_unit, sum(Z_masses[idxs[i]]; init=0.0*m_unit))
                    # Dust mass [m_unit]
                    data_matrix[i, 8] = ustrip(m_unit, sum(dust_masses[idxs[i]]; init=0.0*m_unit))

                    if isempty(idxs[i])

                        # If the voxel has no particles set the velocity to NaN
                        data_matrix[i, 9]  = NaN
                        data_matrix[i, 10] = NaN
                        data_matrix[i, 11] = NaN

                        # If the voxel has no particles set the velocity dispersion to NaN
                        data_matrix[i, 12] = NaN
                        data_matrix[i, 13] = NaN
                        data_matrix[i, 14] = NaN

                    elseif isone(length(idxs[i]))

                        # Velocity in the x direction [v_unit]
                        data_matrix[i, 9]  = velocities[1, idxs[i][1]]
                        # Velocity in the y direction [v_unit]
                        data_matrix[i, 10] = velocities[2, idxs[i][1]]
                        # Velocity in the z direction [v_unit]
                        data_matrix[i, 11] = velocities[3, idxs[i][1]]

                        # If the voxel has a single particle set the velocity dispersion to NaN
                        data_matrix[i, 12] = NaN
                        data_matrix[i, 13] = NaN
                        data_matrix[i, 14] = NaN

                    else

                        # Velocities in the x direction of the particles within the voxel [v_unit]
                        vxs = velocities[1, idxs[i]]
                        # Velocities in the y direction of the particles within the voxel [v_unit]
                        vys = velocities[2, idxs[i]]
                        # Velocities in the z direction of the particles within the voxel [v_unit]
                        vzs = velocities[3, idxs[i]]

                        # Mean and standard deviation of the velocities in the x direction [v_unit]
                        data_matrix[i, 9], data_matrix[i, 12]  = mean_and_std(vxs)
                        # Mean and standard deviation of the velocities in the y direction [v_unit]
                        data_matrix[i, 10], data_matrix[i, 13] = mean_and_std(vys)
                        # Mean and standard deviation of the velocities in the z direction [v_unit]
                        data_matrix[i, 11], data_matrix[i, 14] = mean_and_std(vzs)

                    end

                end

            else

                throw(ArgumentError("gasVelocityCubes: The argument `gas_type` must be :cells or \
                :particles, but I got :$(gas_type)"))

            end

            if row_major_order
                # Go from column-major order (used in Julia, MATLAB, and Fortran) to
                # row-major order (used in Python and C), for interoperability
                hdf5_group["snap_$(snapshot_number)", shuffle=(), deflate=5] = permutedims(
                    data_matrix,
                    reverse(1:ndims(data_matrix)),
                )
            else
                # Stay in column-major order
                hdf5_group["snap_$(snapshot_number)", shuffle=(), deflate=5] = data_matrix
            end

            # Read the time, scale factor, and redshift
            pt = ustrip.(u"Gyr", data_dict[:snap_data].physical_time)
            sf = data_dict[:snap_data].scale_factor
            rs = data_dict[:snap_data].redshift

            # Write the time metadata
            attrs(hdf5_group["snap_$(snapshot_number)"])["Time [Gyr]"]   = pt
            attrs(hdf5_group["snap_$(snapshot_number)"])["Scale factor"] = sf
            attrs(hdf5_group["snap_$(snapshot_number)"])["Redshift"]     = rs

            # Write the unit metadata
            attrs(hdf5_group["snap_$(snapshot_number)"])["Mass unit"]     = string(m_unit)
            attrs(hdf5_group["snap_$(snapshot_number)"])["Length unit"]   = string(l_unit)
            attrs(hdf5_group["snap_$(snapshot_number)"])["Velocity unit"] = string(v_unit)

            # Write the grid metadata
            attrs(hdf5_group["snap_$(snapshot_number)"])["Grid size [length unit]"] = ustrip.(
                l_unit,
                grid.grid_size,
            )
            attrs(hdf5_group["snap_$(snapshot_number)"])["Grid size [# voxels]"] = grid.n_bins

            # Write the column names
            attrs(hdf5_group["snap_$(snapshot_number)"])["Columns"] = [
                "x",     # Column 01: x coordinate [l_unit]
                "y",     # Column 02: y coordinate [l_unit]
                "z",     # Column 03: z coordinate [l_unit]
                "MHII",  # Column 04: Ionized hydrogen mass [m_unit]
                "MHI",   # Column 05: Atomic hydrogen mass [m_unit]
                "MH2",   # Column 06: Molecular hydrogen + stars mass [m_unit]
                "MZ",    # Column 07: Metals mass [m_unit]
                "Mdust", # Column 08: Dust mass [m_unit]
                "Vx",    # Column 09: Gas velocity in the x direction [v_unit]
                "Vy",    # Column 10: Gas velocity in the y direction [v_unit]
                "Vz",    # Column 11: Gas velocity in the z direction [v_unit]
                "Sx",    # Column 12: Gas velocity dispersion in the x direction [v_unit]
                "Sy",    # Column 13: Gas velocity dispersion in the y direction [v_unit]
                "Sz",    # Column 14: Gas velocity dispersion in the z direction [v_unit]
            ]

            next!(prog_bar)

        end

    end

    close(hdf5_file)

    return nothing

end

"""
    stellarVelocityCubes(
        simulation_paths::Vector{String},
        slice::ReducedIndexType;
        <keyword arguments>
    )::Nothing

Create a HDF5 file with the position, stellar mass, velocity, and velocity dispersion of each voxel in a rectangular 3D grid.

The metadata for each snapshot in the HDF5 file includes the physical time in Gyr, the scale factor, and the redshift of that snapshot.

By default, the grid is centered at coordinates (0, 0, 0), has 100x100x100 voxels, and has a side length of [`BOX_L`](@ref). There are as many rows as there are voxels (1000000 by default).

The quantities in the HDF5 file for each voxel are:

Column 01: x coordinate [`l_unit`]
Column 02: y coordinate [`l_unit`]
Column 03: z coordinate [`l_unit`]
Column 04: Stellar mass [`m_unit`]
Column 05: Stellar velocity in the x direction [`v_unit`]
Column 06: Stellar velocity in the y direction [`v_unit`]
Column 07: Stellar velocity in the z direction [`v_unit`]
Column 08: Stellar velocity dispersion in the x direction [`v_unit`]
Column 09: Stellar velocity dispersion in the y direction [`v_unit`]
Column 10: Stellar velocity dispersion in the z direction [`v_unit`]

The mass is the accumulated mass of the particles within each voxel. The velocity is the mean of the velocities of those particles, and the velocity dispersion is the corresponding standard deviation.

If there are no particles, the mass is 0, and the velocity and velocity dispersion are set to NaN. If there is only one particle, the mass and velocity are the ones from that particle, and the velocity dispersion is set to NaN.

By default (`trans_mode` = :all_box and `filter_mode` = :all) we use the following reference system:

  - The origin is the global center of mass of the simulation.
  - The x, y, and z axis form a right-handed cartesian reference system (x × y = z), where the z axis has the orientation of the global angular momentum, and the x and y axis are roughly in the direction of the global principal axis.

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`. All the simulations will be written to the same file.
  - `slice::ReducedIndexType`: Slice of the simulations, i.e. which snapshots will be plotted. It can be an integer (a single snapshot), a vector of integers (several snapshots), an `UnitRange` (e.g. 5:13) or an `StepRange` (e.g. 5:2:13). Starts at 1 and out of bounds indices are ignored.
  - `grid::CubicGrid=CubicGrid(BOX_L, 300)`: Cubic grid.
  - `row_major_order::Bool=true`: Store the results in row-major order (as used in C and Python) instead of column-major order (used in Julia, Fortran, and MATLAB). See [Row- and column-major order](https://en.wikipedia.org/wiki/Row-_and_column-major_order).
  - `m_unit::Unitful.Units=u"Msun"`: Mass unit
  - `l_unit::Unitful.Units=u"kpc"`: Length unit.
  - `v_unit::Unitful.Units=u"km * s^-1"`: Velocity unit.
  - `output_file::String="./stellar_velocity_cube.hdf5"`: Path to the output file.
  - `trans_mode::Union{Symbol,Tuple{TranslationType,RotationType,Dict{Symbol,Vector{String}}}}=:all_box`: How to translate and rotate the cells/particles, before filtering with `filter_mode`. For options see [`selectTransformation`](@ref).
  - `filter_mode::Union{Symbol,Tuple{Function,Dict{Symbol,Vector{String}}}}=:all`: Which cells/particles will be selected. For options see [`selectFilter`](@ref).
  - `show_progress::Bool=true`: If a progress bar will be shown.
"""
function stellarVelocityCubes(
    simulation_paths::Vector{String},
    slice::ReducedIndexType;
    grid::CubicGrid=CubicGrid(BOX_L, 300),
    row_major_order::Bool=true,
    m_unit::Unitful.Units=u"Msun",
    l_unit::Unitful.Units=u"kpc",
    v_unit::Unitful.Units=u"km * s^-1",
    output_file::String="./stellar_velocity_cube.hdf5",
    trans_mode::Union{Symbol,Tuple{TranslationType,RotationType,Dict{Symbol,Vector{String}}}}=:all_box,
    filter_mode::Union{Symbol,Tuple{Function,Dict{Symbol,Vector{String}}}}=:all,
    show_progress::Bool=true,
)::Nothing

    # Create the output folder
    mkpath(dirname(output_file))

    # Create the output HDF5 file
    hdf5_file = h5open(output_file, "w")

    # Set the number of columns and rows
    n_rows = grid.n_bins^3
    n_cols = 10

    base_request = Dict(:stellar => ["MASS", "POS ", "VEL "])

    translation, rotation, trans_request = selectTransformation(trans_mode, base_request)
    filter_function, request = selectFilter(filter_mode, trans_request)

    for simulation_path in simulation_paths

        simulation_name = basename(simulation_path)

        prog_bar = Progress(
            length(slice),
            dt=0.5,
            desc="Writing the velocity cube for simulation $(simulation_name)... ",
            color=:blue,
            barglyphs=BarGlyphs("|#  |"),
            enabled=show_progress,
        )

        # Create an HDF5 group for each simulation
        hdf5_group = create_group(hdf5_file, simulation_name)

        for snap_n in slice

            data_dict = makeDataDict(simulation_path, snap_n, request)

            snapshot_number = lpad(string(data_dict[:snap_data].global_index), 3, "0")

            # Translate the data
            translateData!(data_dict, translation...)

            # Rotate the data
            rotateData!(data_dict, rotation...)

            # Filter the data
            filterData!(data_dict; filter_function)

            # Load the particle positions
            positions = data_dict[:stellar]["POS "]

            # Load the stellar velocities
            velocities = ustrip.(v_unit, data_dict[:stellar]["VEL "])

            # Compute the stellar mass in each cell
            masses = scatterQty(data_dict, :stellar_mass)

            if any(isempty, [masses, velocities, positions])
                throw(ArgumentError("stellarVelocityCubes: Some data is missing (there appears to \
                be no stars in the snapshot), so I cannot construct the velocity cube"))
            end

            # Alocate memory for:
            # Column 01: x coordinate [l_unit]
            # Column 02: y coordinate [l_unit]
            # Column 03: z coordinate [l_unit]
            # Column 04: Stellar mass [m_unit]
            # Column 05: Stellar velocity in the x direction [v_unit]
            # Column 06: Stellar velocity in the y direction [v_unit]
            # Column 07: Stellar velocity in the z direction [v_unit]
            # Column 08: Stellar velocity dispersion in the x direction [v_unit]
            # Column 09: Stellar velocity dispersion in the y direction [v_unit]
            # Column 10: Stellar velocity dispersion in the z direction [v_unit]
            data_matrix = Matrix{Float64}(undef, n_rows, n_cols)

            # Find which particles are within each voxel
            idxs = listHistogram3D(positions, grid)

            Threads.@threads for i in eachindex(grid.grid)

                # Physical coordinates of the voxel [l_unit]
                data_matrix[i, 1:3] .= ustrip.(l_unit, grid.grid[i])

                # Stellar mass [m_unit]
                data_matrix[i, 4] = ustrip(m_unit, sum(masses[idxs[i]]; init=0.0*m_unit))

                if isempty(idxs[i])

                    # If the voxel has no particles set the velocity to NaN
                    data_matrix[i, 5] = NaN
                    data_matrix[i, 6] = NaN
                    data_matrix[i, 7] = NaN

                    # If the voxel has no particles set the velocity dispersion to NaN
                    data_matrix[i, 8]  = NaN
                    data_matrix[i, 9]  = NaN
                    data_matrix[i, 10] = NaN

                elseif isone(length(idxs[i]))

                    # Velocity in the x direction [v_unit]
                    data_matrix[i, 5] = velocities[1, idxs[i][1]]
                    # Velocity in the y direction [v_unit]
                    data_matrix[i, 6] = velocities[2, idxs[i][1]]
                    # Velocity in the z direction [v_unit]
                    data_matrix[i, 7] = velocities[3, idxs[i][1]]

                    # If the voxel has a single particle set the velocity dispersion to NaN
                    data_matrix[i, 8]  = NaN
                    data_matrix[i, 9]  = NaN
                    data_matrix[i, 10] = NaN

                else

                    # Velocities in the x direction of the particles within the voxel [v_unit]
                    vxs = velocities[1, idxs[i]]
                    # Velocities in the y direction of the particles within the voxel [v_unit]
                    vys = velocities[2, idxs[i]]
                    # Velocities in the z direction of the particles within the voxel [v_unit]
                    vzs = velocities[3, idxs[i]]

                    # Mean and standard deviation of the velocities in the x direction [v_unit]
                    data_matrix[i, 5], data_matrix[i, 8] = mean_and_std(vxs)
                    # Mean and standard deviation of the velocities in the y direction [v_unit]
                    data_matrix[i, 6], data_matrix[i, 9] = mean_and_std(vys)
                    # Mean and standard deviation of the velocities in the z direction [v_unit]
                    data_matrix[i, 7], data_matrix[i, 10] = mean_and_std(vzs)

                end

            end

            if row_major_order
                # Go from column-major order (used in Julia, MATLAB, and Fortran) to
                # row-major order (used in Python and C), for interoperability
                hdf5_group["snap_$(snapshot_number)", shuffle=(), deflate=5] = permutedims(
                    data_matrix,
                    reverse(1:ndims(data_matrix)),
                )
            else
                # Stay in column-major order
                hdf5_group["snap_$(snapshot_number)", shuffle=(), deflate=5] = data_matrix
            end

            # Read the time, scale factor, and redshift
            pt = ustrip.(u"Gyr", data_dict[:snap_data].physical_time)
            sf = data_dict[:snap_data].scale_factor
            rs = data_dict[:snap_data].redshift

            # Write the time metadata
            attrs(hdf5_group["snap_$(snapshot_number)"])["Time [Gyr]"]   = pt
            attrs(hdf5_group["snap_$(snapshot_number)"])["Scale factor"] = sf
            attrs(hdf5_group["snap_$(snapshot_number)"])["Redshift"]     = rs

            # Write the unit metadata
            attrs(hdf5_group["snap_$(snapshot_number)"])["Mass unit"]     = string(m_unit)
            attrs(hdf5_group["snap_$(snapshot_number)"])["Length unit"]   = string(l_unit)
            attrs(hdf5_group["snap_$(snapshot_number)"])["Velocity unit"] = string(v_unit)

            # Write the grid metadata
            attrs(hdf5_group["snap_$(snapshot_number)"])["Grid size [length unit]"] = ustrip.(
                l_unit,
                grid.grid_size,
            )
            attrs(hdf5_group["snap_$(snapshot_number)"])["Grid size [# voxels]"] = grid.n_bins

            # Write the column names
            attrs(hdf5_group["snap_$(snapshot_number)"])["Columns"] = [
                "x",  # Column 01: x coordinate [l_unit]
                "y",  # Column 02: y coordinate [l_unit]
                "z",  # Column 03: z coordinate [l_unit]
                "M*", # Column 04: Stellar mass [m_unit]
                "Vx", # Column 05: Stellar velocity in the x direction [v_unit]
                "Vy", # Column 06: Stellar velocity in the y direction [v_unit]
                "Vz", # Column 07: Stellar velocity in the z direction [v_unit]
                "Sx", # Column 08: Stellar velocity dispersion in the x direction [v_unit]
                "Sy", # Column 09: Stellar velocity dispersion in the y direction [v_unit]
                "Sz", # Column 10: Stellar velocity dispersion in the z direction [v_unit]
            ]

            next!(prog_bar)

        end

    end

    close(hdf5_file)

    return nothing

end

"""
    virialAccretionEvolution(
        simulation_paths::Vector{String};
        <keyword arguments>
    )::Nothing

Plot a time series of the integrated mass flux into a sphere with the virial radius.

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`. All the simulations will be plotted together.
  - `slice::IndexType=(:)`: Slice of the simulation, i.e. which snapshots will be plotted. It can be an integer (a single snapshot), a vector of integers (several snapshots), an `UnitRange` (e.g. 5:13), an `StepRange` (e.g. 5:2:13) or (:) (all snapshots). Starts at 1 and out of bounds indices are ignored.
  - `flux_direction::Symbol=:net`: What flux direction will be plotted. The options are:

      + `:net_mass`     -> Net accreted mass.
      + `:inflow_mass`  -> Inflow mass only.
      + `:outflow_mass` -> Outflow mass only.
  - `halo_idx::Int=1`: Index of the target halo (FoF group). Starts at 1.
  - `component::Symbol=:all`: Component to compute the accreted mass for. The options are:

      + `:dark_matter` -> Dark matter.
      + `:black_hole`  -> Black holes.
      + `:gas`         -> Gas.
      + `:stellar`     -> Stars.
      + `:all`         -> All the matter.
  - `tracers::Bool=false`: If tracers will be use to compute the mass accretion.
  - `ylog::Bool=false`: If true, sets the y axis to the ``\\log_{10}`` of the mass flux.
  - `smooth::Int=0`: The time series will be smoothed out using `smooth` bins. Set it to 0 if you want no smoothing.
  - `output_path::String="."`: Path to the output folder.
  - `sim_labels::Union{Vector{<:AbstractString},Nothing}=nothing`: Labels for the plot legend, one per simulation. Set it to `nothing` if you don't want a legend.
  - `theme::Attributes=Theme()`: Plot theme that will take precedence over [`DEFAULT_THEME`](@ref).
"""
function virialAccretionEvolution(
    simulation_paths::Vector{String};
    slice::IndexType=(:),
    flux_direction::Symbol=:net_mass,
    halo_idx::Int=1,
    component::Symbol=:all,
    tracers::Bool=false,
    ylog::Bool=false,
    smooth::Int=0,
    output_path::String=".",
    sim_labels::Union{Vector{<:AbstractString},Nothing}=nothing,
    theme::Attributes=Theme(),
)::Nothing

    x_plot_params = plotParams(:physical_time)
    y_plot_params = plotParams(:mass_accretion)

    if component == :stellar

        c_label = "\\star, \\,\\,"

    elseif component == :dark_matter

        c_label = "\\text{DM, }"

    elseif component == :black_hole

        c_label = "\\text{BH, }"

    elseif component == :gas

        c_label = "\\text{gas, }"

    elseif component == :all

        c_label = ""

    else

        throw(ArgumentError("virialAccretionEvolution: `component` can only be :gas, :stellar, \
        :dark_matter, :black_hole or :all, but I got :$(component)"))

    end

    if flux_direction == :net_mass

        yaxis_var_name = L"\dot{M}_{%$(c_label)\text{net}}^\text{R200}"

    elseif flux_direction == :inflow_mass

        yaxis_var_name = L"\dot{M}_{%$(c_label)\text{inflow}}^\text{R200}"

    elseif flux_direction == :outflow_mass

        yaxis_var_name = L"\dot{M}_{%$(c_label)\text{outflow}}^\text{R200}"

    else

        throw(ArgumentError("virialAccretionEvolution: `flux_direction` can only be :net_mass, \
        :inflow_mass or :outflow_mass, but I got :$(flux_direction)"))

    end

    if ylog

        y_log        = y_plot_params.unit
        y_unit       = Unitful.NoUnits
        y_exp_factor = 0
        yaxis_label  = L"\log_{10} \, " * getLabel("auto_label", 0, y_plot_params.unit)

        post_processing = getNothing
        pp_args         = ()
        pp_kwargs       = (;)

    else

        y_log        = nothing
        y_unit       = y_plot_params.unit
        y_exp_factor = y_plot_params.exp_factor
        yaxis_label  = y_plot_params.axis_label

        post_processing = ppHorizontalFlags!
        pp_args         = ([0.0],)
        pp_kwargs       = (; colors=[:gray65], line_styles=[:solid])

    end

    if tracers
        filename="$(component)_$(flux_direction)_virial_accretion_with_tracers"
    else
        filename="$(component)_$(flux_direction)_virial_accretion"
    end

    plotTimeSeries(
        simulation_paths,
        [lines!];
        output_path,
        filename,
        slice,
        da_functions=[daVirialAccretion],
        da_args=[(component,)],
        da_kwargs=[(; flux_direction, halo_idx, tracers, y_log, smooth)],
        post_processing,
        pp_args,
        pp_kwargs,
        x_unit=x_plot_params.unit,
        y_unit,
        x_exp_factor=x_plot_params.exp_factor,
        y_exp_factor,
        xaxis_label=x_plot_params.axis_label,
        yaxis_label,
        xaxis_var_name=x_plot_params.var_name,
        yaxis_var_name,
        theme=merge(
            theme,
            Theme(
                size=(1400, 880),
                figure_padding=(10, 15, 5, 15),
                palette=(linestyle=[:solid],),
                Axis=(aspect=nothing, xticks=0:14),
            ),
        ),
        sim_labels,
    )

    return nothing

end

"""
    diskAccretionEvolution(
        simulation_paths::Vector{String};
        <keyword arguments>
    )::Nothing

Plot a time series of the integrated mass flux into a given disk.

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`. All the simulations will be plotted together.
  - `slice::IndexType=(:)`: Slice of the simulation, i.e. which snapshots will be plotted. It can be an integer (a single snapshot), a vector of integers (several snapshots), an `UnitRange` (e.g. 5:13), an `StepRange` (e.g. 5:2:13) or (:) (all snapshots). Starts at 1 and out of bounds indices are ignored.
  - `flux_direction::Symbol=:net`: What flux direction will be plotted. The options are:

      + `:net_mass`     -> Net accreted mass.
      + `:inflow_mass`  -> Inflow mass only.
      + `:outflow_mass` -> Outflow mass only.
  - `max_r::Unitful.Length=DISK_R`: Radius of the disk.
  - `max_z::Unitful.Length=5.0u"kpc"`: Half height of the disk.
  - `trans_mode::Union{Symbol,Tuple{TranslationType,RotationType,Dict{Symbol,Vector{String}}}}=:all_box`: How to translate and rotate the cells/particles. For options see [`selectTransformation`](@ref).
  - `component::Symbol=:all`: Component to compute the accreted mass for. The options are:

      + `:dark_matter` -> Dark matter.
      + `:black_hole`  -> Black holes.
      + `:gas`         -> Gas.
      + `:stellar`     -> Stars.
      + `:all`         -> All the matter.
  - `tracers::Bool=false`: If tracers will be use to compute the mass accretion.
  - `ylog::Bool=false`: If true, sets the y axis to the ``\\log_{10}`` of the mass flux.
  - `smooth::Int=0`: The time series will be smoothed out using `smooth` bins. Set it to 0 if you want no smoothing.
  - `output_path::String="."`: Path to the output folder.
  - `sim_labels::Union{Vector{<:AbstractString},Nothing}=nothing`: Labels for the plot legend, one per simulation. Set it to `nothing` if you don't want a legend.
  - `theme::Attributes=Theme()`: Plot theme that will take precedence over [`DEFAULT_THEME`](@ref).
"""
function diskAccretionEvolution(
    simulation_paths::Vector{String};
    slice::IndexType=(:),
    flux_direction::Symbol=:net_mass,
    max_r::Unitful.Length=DISK_R,
    max_z::Unitful.Length=5.0u"kpc",
    trans_mode::Union{Symbol,Tuple{TranslationType,RotationType,Dict{Symbol,Vector{String}}}}=:all_box,
    component::Symbol=:all,
    tracers::Bool=false,
    ylog::Bool=false,
    smooth::Int=0,
    output_path::String=".",
    sim_labels::Union{Vector{<:AbstractString},Nothing}=nothing,
    theme::Attributes=Theme(),
)::Nothing

    x_plot_params = plotParams(:physical_time)
    y_plot_params = plotParams(:mass_accretion)

    if component == :stellar

        c_label = "\\star, \\,\\,"

    elseif component == :dark_matter

        c_label = "\\text{DM, }"

    elseif component == :black_hole

        c_label = "\\text{BH, }"

    elseif component == :gas

        c_label = "\\text{gas, }"

    elseif component == :all

        c_label = ""

    else

        throw(ArgumentError("diskAccretionEvolution: `component` can only be :gas, :stellar, \
        :dark_matter, :black_hole or :all, but I got :$(component)"))

    end

    if flux_direction == :net_mass

        yaxis_var_name = L"\dot{M}_{%$(c_label)\text{net}}^\text{disk}"

    elseif flux_direction == :inflow_mass

        yaxis_var_name = L"\dot{M}_{%$(c_label)\text{inflow}}^\text{disk}"

    elseif flux_direction == :outflow_mass

        yaxis_var_name = L"\dot{M}_{%$(c_label)\text{outflow}}^\text{disk}"

    else

        throw(ArgumentError("diskAccretionEvolution: `flux_direction` can only be :net_mass, \
        :inflow_mass or :outflow_mass, but I got :$(flux_direction)"))

    end

    if ylog

        y_log        = y_plot_params.unit
        y_unit       = Unitful.NoUnits
        y_exp_factor = 0
        yaxis_label  = L"\log_{10} \, " * getLabel("auto_label", 0, y_plot_params.unit)

        post_processing = getNothing
        pp_args         = ()
        pp_kwargs       = (;)

    else

        y_log        = nothing
        y_unit       = y_plot_params.unit
        y_exp_factor = y_plot_params.exp_factor
        yaxis_label  = y_plot_params.axis_label

        post_processing = ppHorizontalFlags!
        pp_args         = ([0.0],)
        pp_kwargs       = (; colors=[:gray65], line_styles=[:solid])

    end

    if tracers
        filename="$(component)_$(flux_direction)_disk_accretion_with_tracers"
    else
        filename="$(component)_$(flux_direction)_disk_accretion"
    end

    plotTimeSeries(
        simulation_paths,
        [lines!];
        output_path,
        filename,
        slice,
        da_functions=[daDiskAccretion],
        da_args=[(component,)],
        da_kwargs=[(; flux_direction, max_r, max_z, trans_mode, tracers, y_log, smooth)],
        post_processing,
        pp_args,
        pp_kwargs,
        x_unit=x_plot_params.unit,
        y_unit,
        x_exp_factor=x_plot_params.exp_factor,
        y_exp_factor,
        xaxis_label=x_plot_params.axis_label,
        yaxis_label,
        xaxis_var_name=x_plot_params.var_name,
        yaxis_var_name,
        theme=merge(
            theme,
            Theme(
                size=(1400, 880),
                figure_padding=(10, 15, 5, 15),
                palette=(linestyle=[:solid],),
                Axis=(aspect=nothing, xticks=0:14),
            ),
        ),
        sim_labels,
    )

    return nothing

end

"""
    fitVSFLaw(
        simulation_path::String,
        slice::IndexType,
        component::Symbol;
        <keyword arguments>
    )::Nothing

Plot the resolved volumetric star formation (VSF) law with an optional linear fit.

!!! note

    Only stars younger than [`AGE_RESOLUTION`](@ref) are considered.

!!! note

    The star formation surface density is just the stellar mass surface density divided by [`AGE_RESOLUTION`](@ref).

# Arguments

  - `simulation_path::String`: Path to the simulation directory, set in the code variable `OutputDir`.
  - `slice::IndexType`: Slice of the simulation, i.e. which snapshots will be plotted. It can be an integer (a single snapshot), a vector of integers (several snapshots), an `UnitRange` (e.g. 5:13), an `StepRange` (e.g. 5:2:13) or (:) (all snapshots). Starts at 1 and out of bounds indices are ignored.
  - `component::Symbol`: Target component. It can only be one of the elements of [`COMPONENTS`](@ref).
  - `field_type::Symbol=:cells`: If the gas surface density will be calculated assuming the gas is in `:particles` or in Voronoi `:cells`.
  - `fit::Bool=true`: If a fit of the plotted values will be added on top of the scatter plot.
  - `box_size::Unitful.Length=BOX_L`: Physical side length for the grids.
  - `x_range::NTuple{2,<:Real}=(-Inf, Inf)`: Only the data withing this range (for the x coordinates) will be fitted.
  - `output_path::String="."`: Path to the output folder.
  - `trans_mode::Union{Symbol,Tuple{TranslationType,RotationType,Dict{Symbol,Vector{String}}}}=:all_box`: How to translate and rotate the cells/particles, before filtering with `filter_mode`. For options see [`selectTransformation`](@ref).
  - `filter_mode::Union{Symbol,Tuple{Function,Dict{Symbol,Vector{String}}}}=:all`: Which cells/particles will be selected. For options see [`selectFilter`](@ref).
  - `sim_label::Union{String,Nothing}=basename(simulation_path)`: Label for the plot legend. Set it to `nothing` if you don't want a legend.
  - `theme::Attributes=Theme()`: Plot theme that will take precedence over [`DEFAULT_THEME`](@ref).
"""
function fitVSFLaw(
    simulation_path::String,
    slice::IndexType,
    component::Symbol;
    field_type::Symbol=:cells,
    fit::Bool=true,
    box_size::Unitful.Length=BOX_L,
    x_range::NTuple{2,<:Real}=(-Inf, Inf),
    output_path::String=".",
    trans_mode::Union{Symbol,Tuple{TranslationType,RotationType,Dict{Symbol,Vector{String}}}}=:all_box,
    filter_mode::Union{Symbol,Tuple{Function,Dict{Symbol,Vector{String}}}}=:all,
    sim_label::Union{String,Nothing}=basename(simulation_path),
    theme::Attributes=Theme(),
)::Nothing

    grid = CubicGrid(box_size, 400)

    x_plot_params = plotParams(Symbol(component, :_mass_density))
    y_plot_params = plotParams(:sfr_density)

    base_request = mergeRequests(x_plot_params.request, y_plot_params.request)

    translation, rotation, trans_request = selectTransformation(trans_mode, base_request)
    filter_function, request = selectFilter(filter_mode, trans_request)

    # Set the labels
    xaxis_label = LaTeXString(
        L"\log_{10} \, " * getLabel(x_plot_params.var_name, 0, x_plot_params.unit)
    )
    yaxis_label = LaTeXString(
        L"\log_{10} \, " * getLabel(y_plot_params.var_name, 0, y_plot_params.unit)
    )

    plotSnapshot(
        [simulation_path],
        request,
        [scatter!];
        pf_kwargs=[(; color=WONG_BLUE, markersize=6, marker=:circle)],
        output_path,
        base_filename="$(basename(simulation_path))_$(component)_vsf_law",
        slice,
        transform_box=true,
        translation,
        rotation,
        filter_function,
        da_functions=[daVSFLaw],
        da_args=[(grid, component)],
        da_kwargs=[(; field_type, stellar_ff=filterByStellarAge)],
        post_processing=fit ? ppFitLine! : getNothing,
        x_trim=x_range,
        xaxis_label,
        yaxis_label,
        theme=merge(theme, Theme(Legend=(nbanks=1, margin=(0, 10, 10, 0)),)),
        sim_labels=[sim_label],
    )

    return nothing

end

"""
    clumpingFactor(
        simulation_paths::Vector{String},
        slice::IndexType,
        component::Symbol;
        <keyword arguments>
    )::Nothing

Plot the clumping factor of `component` for different volume scales.

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`. All the simulations will be plotted together.
  - `slice::IndexType`: Slice of the simulation, i.e. which snapshots will be plotted. It can be an integer (a single snapshot), a vector of integers (several snapshots), an `UnitRange` (e.g. 5:13), an `StepRange` (e.g. 5:2:13) or (:) (all snapshots). Starts at 1 and out of bounds indices are ignored.
  - `component::Symbol`: Target component. It can only be one of the elements of [`COMPONENTS`](@ref) with cell/partcile type :gas.
  - `n_neighbors::Int=32`: Number of neighbors.
  - `smooth::Int=0`: The result will be average out using `smooth` bins for the volume. Set it to 0 if you want no smoothing.
  - `xlog::Bool=false`: If true, sets the x axis to the ``\\log_{10}`` of the volume.
  - `ylog::Bool=false`: If true, sets the y axis to the ``\\log_{10}`` of the clumping factor.
  - `l_unit::Unitful.Units=u"kpc"`: Lenght unit for the volume.
  - `output_path::String="."`: Path to the output folder.
  - `trans_mode::Union{Symbol,Tuple{TranslationType,RotationType,Dict{Symbol,Vector{String}}}}=:all_box`: How to translate and rotate the cells/particles, before filtering with `filter_mode`. For options see [`selectTransformation`](@ref).
  - `filter_mode::Union{Symbol,Tuple{Function,Dict{Symbol,Vector{String}}}}=:all`: Which cells/particles will be selected. For options see [`selectFilter`](@ref).
  - `da_ff::Function=filterNothing`: Filter function to be applied within [`daScatterGalaxy`](@ref) after `trans_mode` and `filter_mode` are applied. See the required signature and examples in `./src/analysis/filters.jl`.
  - `ff_request::Dict{Symbol,Vector{String}}=Dict{Symbol,Vector{String}}()`: Request dictionary for `da_ff`.
  - `sim_labels::Union{Vector{<:AbstractString},Nothing}=nothing`: Labels for the plot legend, one per simulation. Set it to `nothing` if you don't want a legend.
  - `title::Union{Symbol,<:AbstractString}=""`: Title for the figure. If left empty, no title is printed. It can also be set to one of the following options:

      + `:physical_time` -> Physical time since the Big Bang.
      + `:lookback_time` -> Physical time left to reach the last snapshot.
      + `:scale_factor`  -> Scale factor (only relevant for cosmological simulations).
      + `:redshift`      -> Redshift (only relevant for cosmological simulations).
  - `theme::Attributes=Theme()`: Plot theme that will take precedence over [`DEFAULT_THEME`](@ref).
"""
function clumpingFactor(
    simulation_paths::Vector{String},
    slice::IndexType,
    component::Symbol;
    n_neighbors::Int=32,
    smooth::Int=0,
    xlog::Bool=false,
    ylog::Bool=false,
    l_unit::Unitful.Units=u"kpc",
    output_path::String=".",
    trans_mode::Union{Symbol,Tuple{TranslationType,RotationType,Dict{Symbol,Vector{String}}}}=:all_box,
    filter_mode::Union{Symbol,Tuple{Function,Dict{Symbol,Vector{String}}}}=:all,
    da_ff::Function=filterNothing,
    ff_request::Dict{Symbol,Vector{String}}=Dict{Symbol,Vector{String}}(),
    sim_labels::Union{Vector{<:AbstractString},Nothing}=nothing,
    title::Union{Symbol,<:AbstractString}="",
    theme::Attributes=Theme(),
)::Nothing

    plot_params = plotParams(Symbol(component, :_clumping_factor))

    base_request = mergeRequests(plot_params.request, Dict(:gas => ["POS "]), ff_request)

    translation, rotation, trans_request = selectTransformation(trans_mode, base_request)
    filter_function, request = selectFilter(filter_mode, trans_request)

    if isone(length(simulation_paths))
        base_filename = "$(basename(first(simulation_paths)))_$(component)_clumping_factor"
    else
        base_filename = "$(component)_clumping_factor"
    end

    # Set arguments for the x axis
    if xlog
        x_log       = l_unit^3
        x_unit      = Unitful.NoUnits
        xaxis_label = L"\log_{10} \, " * getLabel(L"\bar{V}", 0, l_unit^3)
    else
        x_log       = nothing
        x_unit      = l_unit^3
        xaxis_label = getLabel(L"\bar{V}", 0, l_unit^3)
    end

    # Set arguments for the y axis
    if ylog
        y_log        = plot_params.unit
        y_unit       = Unitful.NoUnits
        y_exp_factor = 0
        yaxis_label  = L"\log_{10} \, " * getLabel("auto_label", 0, plot_params.unit)
    else
        y_log        = nothing
        y_unit       = plot_params.unit
        y_exp_factor = plot_params.exp_factor
        yaxis_label  = plot_params.axis_label
    end

    plotSnapshot(
        simulation_paths,
        request,
        [scatter!];
        pf_kwargs=[(; markersize=3)],
        output_path,
        base_filename,
        slice,
        transform_box=true,
        translation,
        rotation,
        filter_function,
        da_functions=[daClumpingFactor],
        da_args=[(component,)],
        da_kwargs=[(; n_neighbors, filter_function=da_ff, x_log, y_log)],
        smooth,
        x_unit,
        y_unit,
        y_exp_factor,
        xaxis_label,
        yaxis_label,
        xaxis_var_name=L"\bar{V}",
        yaxis_var_name=plot_params.var_name,
        theme,
        sim_labels,
        title,
    )

    return nothing

end

"""
    circularityHistogram(
        simulation_paths::Vector{String},
        slice::IndexType;
        <keyword arguments>
    )::Nothing

Plot three stellar circularity histograms for each simulation. One for the stars within `R_out`, another for the stars within `R_in`, and another for the stars between `R_in` and `R_out`.

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`. Each simulation will be plotted in a different figure.
  - `slice::IndexType`: Slice of the simulation, i.e. which snapshots will be plotted. It can be an integer (a single snapshot), a vector of integers (several snapshots), an `UnitRange` (e.g. 5:13), an `StepRange` (e.g. 5:2:13) or (:) (all snapshots). Starts at 1 and out of bounds indices are ignored.
  - `R_in::Unitful.Length=2.0u"kpc"`: Internal radius.
  - `R_out::Unitful.Length=DISK_R`: External radius.
  - `output_path::String="."`: Path to the output folder.
  - `trans_mode::Union{Symbol,Tuple{TranslationType,RotationType,Dict{Symbol,Vector{String}}}}=:all_box`: How to translate and rotate the cells/particles, before filtering with `filter_mode`. For options see [`selectTransformation`](@ref).
  - `filter_mode::Union{Symbol,Tuple{Function,Dict{Symbol,Vector{String}}}}=:all`: Which cells/particles will be selected. For options see [`selectFilter`](@ref).
  - `theme::Attributes=Theme()`: Plot theme that will take precedence over [`DEFAULT_THEME`](@ref).
"""
function circularityHistogram(
    simulation_paths::Vector{String},
    slice::IndexType;
    R_in::Unitful.Length=2.0u"kpc",
    R_out::Unitful.Length=DISK_R,
    output_path::String=".",
    trans_mode::Union{Symbol,Tuple{TranslationType,RotationType,Dict{Symbol,Vector{String}}}}=:all_box,
    filter_mode::Union{Symbol,Tuple{Function,Dict{Symbol,Vector{String}}}}=:all,
    theme::Attributes=Theme(),
)::Nothing

    plot_params = plotParams(:stellar_circularity)

    translation, rotation, trans_request = selectTransformation(trans_mode, plot_params.request)
    filter_function, request = selectFilter(filter_mode, trans_request)

    grid = LinearGrid(-2.0, 2.0, 200)

    l_unit = u"kpc"

    R_out_label = string(round(Int, ustrip(l_unit, R_out))) * "\\, \\mathrm{$(l_unit)}"
    R_in_label  = string(round(Int, ustrip(l_unit, R_in))) * "\\, \\mathrm{$(l_unit)}"

    da_ff = [
        dd -> filterBySphere(dd, 0.0 * l_unit, R_out, :zero),
        dd -> filterBySphere(dd, 0.0 * l_unit, R_in, :zero),
        dd -> filterBySphere(dd, R_in, R_out, :zero),
    ]

    for simulation_path in simulation_paths

        plotSnapshot(
            [simulation_path, simulation_path, simulation_path],
            request,
            [lines!];
            output_path,
            base_filename="circularity_histogram",
            slice,
            transform_box=true,
            translation,
            rotation,
            filter_function,
            da_functions=[daHistogram],
            da_args=[(:stellar_circularity, grid)],
            da_kwargs=[
                (; norm=1, filter_function=da_ff[1]),
                (; norm=1, filter_function=da_ff[2]),
                (; norm=1, filter_function=da_ff[3]),
            ],
            save_figures=false,
            backup_results=true,
            sim_labels=["simulation_1", "simulation_2", "simulation_3"],
        )

        current_theme = merge(
            theme,
            Theme(
                palette=(color=[:gray65, :orangered2, :navy], linestyle=[:solid]),
                Legend=(nbanks=1, margin=(15, 0, 0, 10), halign=:left, valign=:top),
            ),
            DEFAULT_THEME,
            theme_latexfonts(),
        )

        with_theme(current_theme) do

            f = Figure()

            ax = CairoMakie.Axis(
                f[1, 1];
                xlabel=getLabel(plot_params.var_name, 0, Unitful.NoUnits),
                ylabel=L"\mathrm{Normalized \,\, counts}",
            )

            labels = [
                L"r \,\, \le \,\, %$(R_out_label)",
                L"r \,\, \le \,\, %$(R_in_label)",
                L"%$(R_in_label) \,\, < \,\, r \,\, \le \,\, %$(R_out_label)",
            ]

            jld2_path = joinpath(output_path, "circularity_histogram.jld2")

            jldopen(jld2_path, "r") do jld2_file

                snaps = keys(jld2_file)

                for snap in snaps

                    x, y1 = jld2_file[snap]["simulation_1"]
                    _, y2 = jld2_file[snap]["simulation_2"]
                    _, y3 = jld2_file[snap]["simulation_3"]

                    norm = maximum(y1)

                    lines!(ax, x, y1 ./ norm; label=labels[1])
                    lines!(ax, x, y2 ./ norm; label=labels[2])
                    lines!(ax, x, y3 ./ norm; label=labels[3])

                    axislegend(
                        ax;
                        position=(current_theme.Legend.halign[], current_theme.Legend.valign[]),
                    )

                    filename = "$(basename(simulation_path))_circularity_histogram_$(snap).png"

                    save(joinpath(output_path, filename), f)

                end

            end

            rm(jld2_path)

        end

    end

    return nothing

end

"""
    efficiencyHistogram(
        simulation_paths::Vector{String},
        slice::IndexType;
        <keyword arguments>
    )::Nothing

Plot two histogram of the efficiency per free-fall time, for each simulation. One for the progenitors of the stars and the other for the gas.

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`. Each simulation will be plotted in a different figure.
  - `slice::IndexType`: Slice of the simulation, i.e. which snapshots will be plotted. It can be an integer (a single snapshot), a vector of integers (several snapshots), an `UnitRange` (e.g. 5:13), an `StepRange` (e.g. 5:2:13) or (:) (all snapshots). Starts at 1 and out of bounds indices are ignored.
  - `range::NTuple{2,Float64}=(1.0e-4, 1.0)`: Range for the efficiency per free-fall time (x axis).
  - `output_path::String="."`: Path to the output folder.
  - `trans_mode::Union{Symbol,Tuple{TranslationType,RotationType,Dict{Symbol,Vector{String}}}}=:all_box`: How to translate and rotate the cells/particles, before filtering with `filter_mode`. For options see [`selectTransformation`](@ref).
  - `filter_mode::Union{Symbol,Tuple{Function,Dict{Symbol,Vector{String}}}}=:all`: Which cells/particles will be selected. For options see [`selectFilter`](@ref).
  - `stellar_ff::Function=filterNothing`: Filter function to be applied to the stellar histogram after `trans_mode` and `filter_mode` are applied. See the required signature and examples in `./src/analysis/filters.jl`.
  - `gas_ff::Function=filterNothing`: Filter function to be applied to the gas histogram after `trans_mode` and `filter_mode` are applied. See the required signature and examples in `./src/analysis/filters.jl`.
  - `ff_request::Dict{Symbol,Vector{String}}=Dict{Symbol,Vector{String}}()`: Request dictionary for `stellar_ff` and `gas_ff`.
  - `labels::Vector{<:AbstractString}=["Stars", "Gas"]`: Legend for the stellar and gas histograms, respectively.
  - `title::Union{Symbol,<:AbstractString}=""`: Title for the figure. If left empty, no title is printed. It can also be set to one of the following options:

      + `:physical_time` -> Physical time since the Big Bang.
      + `:lookback_time` -> Physical time left to reach the last snapshot.
      + `:scale_factor`  -> Scale factor (only relevant for cosmological simulations).
      + `:redshift`      -> Redshift (only relevant for cosmological simulations).
  - `theme::Attributes=Theme()`: Plot theme that will take precedence over [`DEFAULT_THEME`](@ref).
"""
function efficiencyHistogram(
    simulation_paths::Vector{String},
    slice::IndexType;
    range::NTuple{2,Float64}=(1.0e-4, 1.0),
    output_path::String=".",
    trans_mode::Union{Symbol,Tuple{TranslationType,RotationType,Dict{Symbol,Vector{String}}}}=:all_box,
    filter_mode::Union{Symbol,Tuple{Function,Dict{Symbol,Vector{String}}}}=:all,
    stellar_ff::Function=filterNothing,
    gas_ff::Function=filterNothing,
    ff_request::Dict{Symbol,Vector{String}}=Dict{Symbol,Vector{String}}(),
    labels::Vector{<:AbstractString}=["Stars", "Gas"],
    title::Union{Symbol,<:AbstractString}="",
    theme::Attributes=Theme(),
)::Nothing

    plot_params = plotParams(:eff)

    base_request = mergeRequests(plot_params.request, ff_request)

    translation, rotation, trans_request = selectTransformation(trans_mode, base_request)
    filter_function, request = selectFilter(filter_mode, trans_request)

    grid = LinearGrid(range..., 100; log=true)

    # Set nice margins for the legend
    l_margin = isempty(title) ? 10 : 25
    t_margin = isempty(title) ? 5 : 10

    for simulation_path in simulation_paths

        plotSnapshot(
            [simulation_path, simulation_path],
            request,
            [lines!];
            output_path,
            base_filename="$(basename(simulation_path))_eff_histogram",
            slice,
            transform_box=true,
            translation,
            rotation,
            filter_function,
            da_functions=[daHistogram, daHistogram],
            da_args=[(:stellar_eff, grid), (:gas_eff, grid)],
            da_kwargs=[(; filter_function=stellar_ff), (; filter_function=gas_ff)],
            xaxis_label=L"\log_{10} \," * getLabel(plot_params.axis_label, 0, plot_params.unit),
            xaxis_var_name=plot_params.var_name,
            yaxis_var_name=L"\mathrm{Normalized \,\, counts}",
            theme=merge(
                theme,
                Theme(
                    palette=(linestyle=[:solid], color=[WONG_ORANGE, :black]),
                    Legend=(nbanks=1, halign=:left, valign=:top, margin=(l_margin, 0, 0, t_margin)),
                ),
            ),
            sim_labels=labels,
            title,
        )

    end

end

"""
    stellarDensityMaps(
        simulation_paths::Vector{String},
        slice::IndexType;
        <keyword arguments>
    )::Nothing

Plot the stellar density maps for the xy and xz projections, in two panels.

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`. Each simulation will be plotted in a different figure.
  - `slice::IndexType`: Slice of the simulation, i.e. which snapshots will be plotted. It can be an integer (a single snapshot), a vector of integers (several snapshots), an `UnitRange` (e.g. 5:13), an `StepRange` (e.g. 5:2:13) or (:) (all snapshots). Starts at 1 and out of bounds indices are ignored.
  - `box_size::Unitful.Length=BOX_L`: Size of the plotting box (x and y coordinates).
  - `box_height::Unitful.Length=12.0u"kpc"`: Size of the plotting box (z coordinate).
  - `output_path::String="."`: Path to the output folder.
  - `trans_mode::Union{Symbol,Tuple{TranslationType,RotationType,Dict{Symbol,Vector{String}}}}=:all_box`: How to translate and rotate the cells/particles, before filtering with `filter_mode`. For options see [`selectTransformation`](@ref).
  - `filter_mode::Union{Symbol,Tuple{Function,Dict{Symbol,Vector{String}}}}=:all`: Which cells/particles will be selected. For options see [`selectFilter`](@ref).
  - `da_ff::Function=filterNothing`: Filter function to be applied within [`daDensity2DProjection`](@ref) after `trans_mode` and `filter_mode` are applied. See the required signature and examples in `./src/analysis/filters.jl`.
  - `ff_request::Dict{Symbol,Vector{String}}=Dict{Symbol,Vector{String}}()`: Request dictionary for `da_ff`.
  - `theme::Attributes=Theme()`: Plot theme that will take precedence over [`DEFAULT_THEME`](@ref).
"""
function stellarDensityMaps(
    simulation_paths::Vector{String},
    slice::IndexType;
    box_size::Unitful.Length=BOX_L,
    box_height::Unitful.Length=12.0u"kpc",
    output_path::String=".",
    trans_mode::Union{Symbol,Tuple{TranslationType,RotationType,Dict{Symbol,Vector{String}}}}=:all_box,
    filter_mode::Union{Symbol,Tuple{Function,Dict{Symbol,Vector{String}}}}=:all,
    da_ff::Function=filterNothing,
    ff_request::Dict{Symbol,Vector{String}}=Dict{Symbol,Vector{String}}(),
    theme::Attributes=Theme(),
)::Nothing

    projection_planes = [:xy, :xz]
    m_unit = u"Msun"
    l_unit = u"kpc"

    grid     = CubicGrid(box_size, 400)
    half_box = ustrip(l_unit, box_size) / 2.0

    # Maximun tick for the axes
    tick = floor(half_box; sigdigits=1)

    x_limits = half_box
    y_limits = [half_box, ustrip(l_unit, box_height)]

    x_label = getLabel("x", 0, l_unit)
    y_label = getLabel("y", 0, l_unit)
    z_label = getLabel("z", 0, l_unit)

    plot_params = plotParams(:stellar_area_density)

    # Label for the colorbar
    colorbar_label = LaTeXString(
        L"\log_{10} \, " * getLabel(plot_params.var_name, 0, m_unit * l_unit^-2)
    )

    base_request = mergeRequests(plot_params.request, ff_request)

    translation, rotation, trans_request = selectTransformation(trans_mode, base_request)
    filter_function, request = selectFilter(filter_mode, trans_request)

    current_theme = merge(
        theme,
        Theme(
            size=(880, 1300),
            figure_padding=(5, 25, 0, 0),
            Axis=(xticklabelsize=35, yticklabelsize=35, aspect=DataAspect()),
            Colorbar=(ticklabelsize=28, vertical=false, ticks=WilkinsonTicks(5)),
        ),
        DEFAULT_THEME,
        theme_latexfonts(),
    )

    for simulation_path in simulation_paths

        temp_folder = joinpath(output_path, "_stellar_density_maps")

        for projection_plane in projection_planes

            plotSnapshot(
                [simulation_path],
                request,
                [heatmap!];
                output_path=temp_folder,
                base_filename="stellar_mass_$(projection_plane)",
                slice,
                transform_box=true,
                translation,
                rotation,
                filter_function,
                da_functions=[daDensity2DProjection],
                da_args=[(grid, :stellar, :particles)],
                da_kwargs=[(; projection_plane, m_unit, l_unit, filter_function=da_ff)],
                x_unit=l_unit,
                y_unit=l_unit,
                save_figures=false,
                backup_results=true,
            )

        end

        jld2_paths = [joinpath(temp_folder, "stellar_mass_$(pp).jld2") for pp in projection_planes]

        jld2_data = load.(jld2_paths)

        for ((snap, xy_data), (_, xz_data)) in zip(jld2_data...)

            with_theme(current_theme) do

                f = Figure()

                min_color = Inf
                max_color = -Inf

                # Compute a good color range
                for (_, _, z) in [xy_data, xz_data]

                    if !all(isnan, z)

                        min_Σ, max_Σ = extrema(filter(!isnan, z))

                        floor_Σ = round(min_Σ)
                        ceil_Σ  = round(max_Σ)

                        if floor_Σ < min_color
                            min_color = floor_Σ
                        end
                        if ceil_Σ > max_color
                            max_color = ceil_Σ
                        end

                    end

                end

                for (row, (x, y, z)) in pairs([xy_data, xz_data])

                    xaxis_v = row == 2

                    ax = CairoMakie.Axis(
                        f[row+1, 1];
                        xlabel=x_label,
                        ylabel=(row == 1 ? y_label : z_label),
                        xminorticksvisible=xaxis_v,
                        xticksvisible=xaxis_v,
                        xlabelvisible=xaxis_v,
                        xticklabelsvisible=xaxis_v,
                        xticks=-tick:10:tick,
                        yticks=-tick:10:tick,
                        limits=(-x_limits, x_limits, -y_limits[row], y_limits[row]),
                    )

                    pf = heatmap!(ax, x, y, z; colorrange=(min_color, max_color))

                    if row == 1
                        Colorbar(f[row, 1], pf; label=colorbar_label)
                    end

                end

                rowsize!(f.layout, 3, Relative(0.3f0))

                filename = "$(basename(snap))_stellar_density_maps_$(dirname(snap)).png"

                save(joinpath(output_path, filename), f)

            end

        end

        rm(temp_folder; recursive=true)

    end

    return nothing

end

"""
    gasDensityMaps(
        simulation_paths::Vector{String},
        slice::IndexType;
        <keyword arguments>
    )::Nothing

Plot the density map of five gas components for the xy and xz projections, in several panels.

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`. Each simulation will be plotted in a different figure.
  - `slice::IndexType`: Slice of the simulation, i.e. which snapshots will be plotted. It can be an integer (a single snapshot), a vector of integers (several snapshots), an `UnitRange` (e.g. 5:13), an `StepRange` (e.g. 5:2:13) or (:) (all snapshots). Starts at 1 and out of bounds indices are ignored.
  - `box_size::Unitful.Length=BOX_L`: Size of the plotting box.
  - `output_path::String="."`: Path to the output folder.
  - `density_range::NTuple{2,Float64}=(NaN,NaN)`: Area density range in ``\\log_{10} \\mathrm{[M_\\odot \\, kpc^{-2}]``. If set to NaN a value is chosen automatically.
  - `trans_mode::Union{Symbol,Tuple{TranslationType,RotationType,Dict{Symbol,Vector{String}}}}=:all_box`: How to translate and rotate the cells/particles, before filtering with `filter_mode`. For options see [`selectTransformation`](@ref).
  - `filter_mode::Union{Symbol,Tuple{Function,Dict{Symbol,Vector{String}}}}=:all`: Which cells/particles will be selected. For options see [`selectFilter`](@ref).
  - `da_ff::Function=filterNothing`: Filter function to be applied within [`daDensity2DProjection`](@ref) after `trans_mode` and `filter_mode` are applied. See the required signature and examples in `./src/analysis/filters.jl`.
  - `ff_request::Dict{Symbol,Vector{String}}=Dict{Symbol,Vector{String}}()`: Request dictionary for `da_ff`.
  - `theme::Attributes=Theme()`: Plot theme that will take precedence over [`DEFAULT_THEME`](@ref).
"""
function gasDensityMaps(
    simulation_paths::Vector{String},
    slice::IndexType;
    box_size::Unitful.Length=BOX_L,
    output_path::String=".",
    density_range::NTuple{2,Float64}=(NaN,NaN),
    trans_mode::Union{Symbol,Tuple{TranslationType,RotationType,Dict{Symbol,Vector{String}}}}=:all_box,
    filter_mode::Union{Symbol,Tuple{Function,Dict{Symbol,Vector{String}}}}=:all,
    da_ff::Function=filterNothing,
    ff_request::Dict{Symbol,Vector{String}}=Dict{Symbol,Vector{String}}(),
    theme::Attributes=Theme(),
)::Nothing

    projection_planes = [:xy, :xz]
    m_unit = u"Msun"
    l_unit = u"kpc"

    grid     = CubicGrid(box_size, 400)
    half_box = ustrip(l_unit, box_size) / 2.0

    # Maximun tick for the axes
    tick = floor(half_box; sigdigits=1)

    x_label = getLabel("x", 0, l_unit)
    y_label = getLabel("y", 0, l_unit)
    z_label = getLabel("z", 0, l_unit)

    for simulation_path in simulation_paths

        temp_folder = joinpath(output_path, "_gas_density_maps")

        if isSimSFM(simulation_path)

            quantities = [
                :gas,
                :ode_ionized,
                :ode_atomic,
                :ode_molecular_stellar,
                :ode_metals,
                :ode_dust,
            ]

            size = (2600, 1020)

        else

            quantities = [:gas, :ionized, :br_atomic, :br_molecular, :Z_gas]

            size = (2200, 1020)

        end

        n_rows = length(projection_planes)
        n_cols = length(quantities)

        jld2_paths = Vector{String}(undef, n_rows * n_cols)
        colorbar_labels = Vector{LaTeXString}(undef, n_cols)

        for (i, quantity) in pairs(quantities)

            plot_params = plotParams(Symbol(quantity, :_area_density))

            colorbar_labels[i] = LaTeXString(
                L"\log_{10} \, " * getLabel(plot_params.var_name, 0, m_unit * l_unit^-2)
            )

            base_request = mergeRequests(plot_params.request, ff_request)

            translation, rotation, trans_request = selectTransformation(trans_mode, base_request)
            filter_function, request = selectFilter(filter_mode, trans_request)

            for (j, projection_plane) in pairs(projection_planes)

                plotSnapshot(
                    [simulation_path],
                    request,
                    [heatmap!];
                    output_path=temp_folder,
                    base_filename="$(quantity)_$(projection_plane)",
                    slice,
                    transform_box=true,
                    translation,
                    rotation,
                    filter_function,
                    da_functions=[daDensity2DProjection],
                    da_args=[(grid, quantity, :cells)],
                    da_kwargs=[(; projection_plane, m_unit, l_unit, filter_function=da_ff)],
                    x_unit=l_unit,
                    y_unit=l_unit,
                    save_figures=false,
                    backup_results=true,
                )

                jld2_paths[i + n_cols * (j - 1)] = joinpath(
                    temp_folder,
                    "$(quantity)_$(projection_plane).jld2",
                )

            end

        end

        jld2_data = load.(jld2_paths)

        current_theme = merge(
            theme,
            Theme(;
                size,
                figure_padding=(5, 10, 5, 0),
                Axis=(xticklabelsize=28, yticklabelsize=28),
                Colorbar=(ticklabelsize=23, vertical=false, ticks=WilkinsonTicks(5)),
            ),
            DEFAULT_THEME,
            theme_latexfonts(),
        )

        for data_dicts in zip(jld2_data...)

            with_theme(current_theme) do

                f = Figure()

                min_color = Inf
                max_color = -Inf

                # Compute a good color range
                for (_, (_, _, z)) in data_dicts

                    if !all(isnan, z)

                        min_Σ, max_Σ = extrema(filter(!isnan, z))

                        floor_Σ = round(min_Σ)
                        ceil_Σ  = round(max_Σ)

                        if floor_Σ < min_color
                            min_color = floor_Σ
                        end
                        if ceil_Σ > max_color
                            max_color = ceil_Σ
                        end

                    end

                end

                if !isnan(density_range[1])
                    min_color = density_range[1]
                end

                if !isnan(density_range[2])
                    max_color = density_range[2]
                end

                for (idx, (_, (x, y, z))) in pairs(data_dicts)

                    row = ceil(Int, idx / n_cols)
                    col = mod1(idx, n_cols)

                    xaxis_v = row == 2
                    yaxis_v = col == 1

                    ax = CairoMakie.Axis(
                        f[row+1, col];
                        xlabel=x_label,
                        ylabel=(row == 1 ? y_label : z_label),
                        xminorticksvisible=xaxis_v,
                        xticksvisible=xaxis_v,
                        xlabelvisible=xaxis_v,
                        xticklabelsvisible=xaxis_v,
                        yminorticksvisible=yaxis_v,
                        yticksvisible=yaxis_v,
                        ylabelvisible=yaxis_v,
                        yticklabelsvisible=yaxis_v,
                        xticks=-tick:10:tick,
                        yticks=-tick:10:tick,
                        limits=(-half_box, half_box, -half_box, half_box),
                    )

                    pf = heatmap!(ax, x, y, z; colorrange=(min_color, max_color))

                    if row == 1

                        Colorbar(f[row, col], pf, label=colorbar_labels[col])

                    end

                end

                colgap!(f.layout, 30)

                snap = first(first(data_dicts))

                filename = "$(basename(snap))_gas_density_maps_$(dirname(snap)).png"

                save(joinpath(output_path, filename), f)

            end

        end

        rm(temp_folder; recursive=true)

    end

    return nothing

end

"""
    evolutionVideo(
        simulation_paths::Vector{String},
        component::Symbol;
        <keyword arguments>
    )::Nothing

Make a video of how the projected density (xy and xz planes) evolves through time, for the stars, gas, and a given gas `component`.

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`. Each simulation will have a different video.
  - `component::Symbol`: Target gas component. See [`COMPONENTS`](@ref) for options.
  - `slice::IndexType=(:)`: Slice of the simulation, i.e. which snapshots will be plotted. It can be an integer (a single snapshot), a vector of integers (several snapshots), an `UnitRange` (e.g. 5:13), an `StepRange` (e.g. 5:2:13) or (:) (all snapshots). Starts at 1 and out of bounds indices are ignored.
  - `field_type::Symbol=cells`: If the gas field is made up of `:particles` or Voronoi `:cells`.
  - `box_size::Unitful.Length=BOX_L`: Size of the plotting box.
  - `output_path::String="."`: Path to the output folder.
  - `density_range::NTuple{2,Float64}=(NaN,NaN)`: Area density range in ``\\log_{10} \\mathrm{[M_\\odot \\, kpc^{-2}]`` for the `component` and gas. If set to NaN a value is chosen automatically.
  - `framerate::Int64=20`: Video framerate.
  - `trans_mode::Union{Symbol,Tuple{TranslationType,RotationType,Dict{Symbol,Vector{String}}}}=:all_box`: How to translate and rotate the cells/particles, before filtering with `filter_mode`. For options see [`selectTransformation`](@ref).
  - `filter_mode::Union{Symbol,Tuple{Function,Dict{Symbol,Vector{String}}}}=:all`: Which cells/particles will be selected. For options see [`selectFilter`](@ref).
  - `da_ff::Function=filterNothing`: Filter function to be applied within [`daDensity2DProjection`](@ref) after `trans_mode` and `filter_mode` are applied. See the required signature and examples in `./src/analysis/filters.jl`.
  - `ff_request::Dict{Symbol,Vector{String}}=Dict{Symbol,Vector{String}}()`: Request dictionary for `da_ff`.
  - `theme::Attributes=Theme()`: Plot theme that will take precedence over [`DEFAULT_THEME`](@ref).
  - `show_progress::Bool=true`: If a progress bar will be shown.
  - `save_data::Bool=false`: If the density maps data will be saved in JLD2 files.
"""
function evolutionVideo(
    simulation_paths::Vector{String},
    component::Symbol;
    slice::IndexType=(:),
    field_type::Symbol=:cells,
    box_size::Unitful.Length=BOX_L,
    output_path::String=".",
    density_range::NTuple{2,Float64}=(NaN,NaN),
    framerate::Int64=20,
    trans_mode::Union{Symbol,Tuple{TranslationType,RotationType,Dict{Symbol,Vector{String}}}}=:all_box,
    filter_mode::Union{Symbol,Tuple{Function,Dict{Symbol,Vector{String}}}}=:all,
    da_ff::Function=filterNothing,
    ff_request::Dict{Symbol,Vector{String}}=Dict{Symbol,Vector{String}}(),
    theme::Attributes=Theme(),
    show_progress::Bool=true,
    save_data::Bool=false,
)::Nothing

    ################################################################################################
    # Constants and parameters
    ################################################################################################

    projection_planes = [:xy, :xz]
    m_unit = u"Msun"
    l_unit = u"kpc"
    Σ_unit = m_unit * l_unit^-2

    grid     = CubicGrid(box_size, 400)
    half_box = ustrip(l_unit, box_size) / 2.0

    # Maximun tick for the axes
    tick = floor(half_box; sigdigits=1)

    x_label = getLabel("x", 0, l_unit)
    y_label = getLabel("y", 0, l_unit)
    z_label = getLabel("z", 0, l_unit)

    x_limits = half_box
    y_limits = [half_box, 12]

    plot_params = plotParams(Symbol(component, :_area_density))
    quantities  = [:gas, component, :stellar]

    labels = LaTeXString.(
        [
            L"\log_{10} \, " * getLabel(plotParams(:gas_area_density).var_name, 0, Σ_unit),
            L"\log_{10} \, " * getLabel(plot_params.var_name, 0, Σ_unit),
            L"\log_{10} \, " * getLabel(plotParams(:stellar_area_density).var_name, 0, Σ_unit),
        ]
    )
    field_types = [field_type, field_type, :particles]

    base_request = mergeRequests(plotParams(:area_density).request, ff_request)

    translation, rotation, trans_request = selectTransformation(trans_mode, base_request)
    filter_function, request = selectFilter(filter_mode, trans_request)

    for simulation_path in simulation_paths

        ############################################################################################
        # Compute the last coordinate trasformation
        ############################################################################################

        simulation_table = makeSimulationTable(simulation_path)

        numbers = safeSelect(simulation_table[!, :numbers], slice)
        times   = safeSelect(simulation_table[!, :physical_times], slice)

        last_address = "$(SNAP_BASENAME)_$(last(numbers))/$(basename(simulation_path))"

        # Read the last snapshot
        last_dd = makeDataDict(
            simulation_path,
            parse(Int, last(numbers)),
            request,
            simulation_table,
        )

        # Compute the translation for the last snapshot
        last_origin = computeCenter(last_dd, translation...)
        last_vcm    = computeVcm(last_dd, translation...)
        translateData!(last_dd, last_origin, last_vcm)

        # Compute the rotation for the last snapshot
        rotation_matrix = computeRotation(last_dd, rotation...)

        ############################################################################################
        # Compute the area densities
        ############################################################################################

        temp_folder = joinpath(output_path, "_density_maps_video")

        for (quantity, field_type) in zip(quantities, field_types)

            for projection_plane in projection_planes

                if save_data

                    jld2_file = joinpath(temp_folder, "$(quantity)_$(projection_plane).jld2")

                    isfile(jld2_file) && continue

                end

                plotSnapshot(
                    [simulation_path],
                    request,
                    [heatmap!];
                    slice,
                    output_path=temp_folder,
                    base_filename="$(quantity)_$(projection_plane)",
                    transform_box=true,
                    translation=(last_origin, last_vcm),
                    rotation=(rotation_matrix,),
                    filter_function,
                    da_functions=[daDensity2DProjection],
                    da_args=[(grid, quantity, field_type)],
                    da_kwargs=[(; projection_plane, m_unit, l_unit, filter_function=da_ff)],
                    x_unit=l_unit,
                    y_unit=l_unit,
                    save_figures=false,
                    backup_results=true,
                )

            end

        end

        ############################################################################################
        # Compute a good stellar color range
        ############################################################################################

        min_color = Inf
        max_color = -Inf

        for projection_plane in projection_planes

            jld2_path = joinpath(temp_folder, "stellar_$(projection_plane).jld2")

            jld2_data = load(jld2_path)

            _, _, z = jld2_data[last_address]

            if !all(isnan, z)

                min_Σ, max_Σ = extrema(filter(!isnan, z))

                floor_Σ = round(min_Σ)
                ceil_Σ  = round(max_Σ)

                if floor_Σ < min_color
                    min_color = floor_Σ
                end
                if ceil_Σ > max_color
                    max_color = ceil_Σ
                end

            end

        end

        if any(isinf, [min_color, max_color])

            stellar_colorrange = (0.0, 1.0)

        else

            stellar_colorrange = (min_color, max_color)

        end

        ############################################################################################
        # Compute a good component and gas color range
        ############################################################################################

        min_color = Inf
        max_color = -Inf

        # Compute a good color range
        for projection_plane in projection_planes

            jld2_path = joinpath(temp_folder, "$(component)_$(projection_plane).jld2")

            jld2_data = load(jld2_path)

            _, _, z = jld2_data[last_address]

            if !all(isnan, z)

                min_Σ, max_Σ = extrema(filter(!isnan, z))

                floor_Σ = round(min_Σ)
                ceil_Σ  = round(max_Σ)

                if floor_Σ < min_color
                    min_color = floor_Σ
                end
                if ceil_Σ > max_color
                    max_color = ceil_Σ
                end

            end

        end

        if !isnan(density_range[1])
            min_color = density_range[1]
        end

        if !isnan(density_range[2])
            max_color = density_range[2]
        end

        if any(isinf, [min_color, max_color])

            gas_colorrange = (0.0, 1.0)

        else

            gas_colorrange = (min_color, max_color)

        end

        colorranges = [gas_colorrange, gas_colorrange, stellar_colorrange]

        ############################################################################################
        # Generate the video
        ############################################################################################

        prog_bar = Progress(
            nrow(simulation_table),
            dt=0.5,
            desc="Analyzing and plotting the data... ",
            color=:blue,
            barglyphs=BarGlyphs("|#  |"),
            enabled=show_progress,
        )

        current_theme = merge(
            theme,
            Theme(
                size=(2400, 1300),
                figure_padding=(5, 25, 0, 0),
                Axis=(xticklabelsize=35, yticklabelsize=35, aspect=DataAspect()),
                Colorbar=(ticklabelsize=28, ticks=WilkinsonTicks(5), vertical=false),
                Text=(color=:white, fontsize=40),
            ),
            DEFAULT_THEME,
            theme_latexfonts(),
        )

        with_theme(current_theme) do

            f = Figure()

            # Initialize the animation stream
            vs = VideoStream(f; framerate)

            col_iterator = enumerate(zip(quantities, labels, colorranges))

            for (snap_n, time) in zip(numbers, times)

                for (col, (quantity, label, colorrange)) in col_iterator

                    for (row, projection_plane) in pairs(projection_planes)

                        xaxis_v = row == 2
                        yaxis_v = col == 1

                        ax = CairoMakie.Axis(
                            f[row+1, col];
                            xlabel=x_label,
                            ylabel=(row == 1 ? y_label : z_label),
                            xminorticksvisible=xaxis_v,
                            xticksvisible=xaxis_v,
                            xlabelvisible=xaxis_v,
                            xticklabelsvisible=xaxis_v,
                            yminorticksvisible=yaxis_v,
                            yticksvisible=yaxis_v,
                            ylabelvisible=yaxis_v,
                            yticklabelsvisible=yaxis_v,
                            xticks=-tick:10:tick,
                            yticks=-tick:10:tick,
                            limits=(-x_limits, x_limits, -y_limits[row], y_limits[row]),
                        )

                        address = "$(SNAP_BASENAME)_$(snap_n)/$(basename(simulation_path))"
                        path = joinpath(temp_folder, "$(quantity)_$(projection_plane).jld2")

                        jldopen(path, "r") do jld2_file

                            x, y, z = jld2_file[address]

                            pf = heatmap!(ax, x, y, z; colorrange)

                            if row == 1

                                Colorbar(f[row, col], pf; label)

                                c_t = ustrip(u"Gyr", time)

                                if c_t < 1.0
                                    time_stamp = round(c_t; digits=2)
                                else
                                    time_stamp = round(c_t; sigdigits=3)
                                end

                                text!(
                                    ax,
                                    0.68,
                                    0.98;
                                    text=L"t = %$(rpad(time_stamp, 4, '0')) \, \text{Gyr}",
                                    align=(:left, :top),
                                    space=:relative,
                                )

                            end

                        end

                    end

                    rowsize!(f.layout, 3, Relative(0.3f0))

                end

                colgap!(f.layout, 70)

                # Add the figure as a frame to the animation stream
                recordframe!(vs)

                cleanPlot!(f)

                next!(prog_bar)

            end

            save(joinpath(output_path, "$(basename(simulation_path))_density_evolution.mkv"), vs)

        end

        if !save_data
            rm(temp_folder; recursive=true)
        end

    end

    return nothing

end

"""
    simulationReport(simulation_paths::Vector{String}; <keyword arguments>)::Nothing

Write a text file with information about a given simulation.

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`. One text file will be written for each simulation.
  - `output_path::String="."`: Path to the output folder.
"""
function simulationReport(simulation_paths::Vector{String}; output_path::String=".")::Nothing

    for simulation_path in simulation_paths

        # Make a dataframe for the simulation with the following columns:
        #  - DataFrame index         -> :row_id
        #  - Number in the file name -> :numbers
        #  - Scale factor            -> :scale_factors
        #  - Redshift                -> :redshifts
        #  - Physical time           -> :physical_times
        #  - Lookback time           -> :lookback_times
        #  - Snapshot path           -> :snapshot_paths
        #  - Group catalog path      -> :groupcat_paths
        simulation_table = makeSimulationTable(simulation_path)

        # Compute the number of snapshots in the folder
        snapshot_n = count(!ismissing, simulation_table[!, :snapshot_paths])

        # Check that there is at least one snapshot
        (
            !iszero(snapshot_n) ||
            throw(ArgumentError("simulationReport: There are no snapshots in $(simulation_path)"))
        )

        # Compute the number of group catalog files in the folder
        groupcat_n = count(!ismissing, simulation_table[!, :groupcat_paths])

        # Check if the simulation is cosmological
        cosmological = isSimCosmological(simulation_path)

        # Read the runtime from the cpu.txt file
        cpu_txt_path = joinpath(simulation_path, CPU_REL_PATH)

        ############################################################################################
        # Print the report header
        ############################################################################################

        mkpath(output_path)

        # Create the output file
        file = open(joinpath(output_path, "report_for_$(basename(simulation_path)).txt"), "w")

        println(file, "#"^100)
        println(file, "\nSimulation name:  $(basename(simulation_path))")
        println(file, "Simulation path:  $(simulation_path)")

        if cosmological
            println(file, "Cosmological:     Yes")
        else
            println(file, "Cosmological:     No")
        end

        if !PHYSICAL_UNITS && cosmological
            println(file, "Lenght units:     Comoving\n")
        else
            println(file, "Lenght units:     Physical\n")
        end

        println(file, "#"^100)
        println(file, "\nNumber of snapshots:       $(snapshot_n)")
        println(file, "Number of group catalogs:  $(groupcat_n)\n")
        println(file, "#"^100)

        println(file)

        ############################################################################################
        # Print the time ranges
        ############################################################################################

        if snapshot_n > 1

            min_t, max_t = round.(
                ustrip.(u"Gyr", extrema(simulation_table[!, :physical_times])),
                digits=2,
            )

            println(file, "Time ranges:\n")

            println(file, "Physical time:  $(min_t) - $(max_t) Gyr")

            if cosmological

                # For cosmological simulations print the scale factor and the redshift
                min_a, max_a = round.(extrema(simulation_table[!, :scale_factors]), digits=3)
                min_z, max_z = round.(extrema(simulation_table[!, :redshifts]), digits=3)

                println(file, "Scale factor:   $(min_a) - $(max_a)")
                println(file, "Redshift:       $(max_z) - $(min_z)")

            end

            if isfile(cpu_txt_path)
                cpu_txt_data = readCpuFile(cpu_txt_path, ["total"])["total"]
                run_time     = cpu_txt_data[end, 5]
                println(file, "Run time:       $(format_seconds(run_time))")
            else
                println(file, "Run time:       Missing cpu.txt file!")
            end

            println(file)

        else

            t = round(ustrip(u"Gyr", simulation_table[1, :physical_times]), digits=2)

            println(file, "Physical time:             $(p) Gyr")

            if cosmological

                # For cosmological simulations print the scale factor and the redshift
                a = round(simulation_table[1, :scale_factors], digits=3)
                z = round(simulation_table[1, :redshifts], digits=3)

                println(file, "Scale factor:              $(a)")
                println(file, "Redshift:                  $(z)")

            end

            println(file)

        end

        ############################################################################################
        # Print the first instance of stars and subhalos
        ############################################################################################

        # Set flags to print only the first instance of each condition
        first_star_flag            = false
        first_subhalo_flag         = false
        first_star_in_subhalo_flag = false

        for snapshot_row in eachrow(simulation_table)

            snapshot_path = snapshot_row[:snapshot_paths]

            # Skip this row if there is no snapshot
            !ismissing(snapshot_path) || continue

            # Read the group catalog path
            groupcat_path = snapshot_row[:groupcat_paths]

            # Read the times
            physical_time = round(ustrip(u"Gyr", snapshot_row[:physical_times]), digits=2)
            if cosmological
                # For cosmological simulations read the scale factor and the redshift
                scale_factor = round(snapshot_row[:scale_factors], digits=3)
                redshift = round(snapshot_row[:redshifts], digits=3)
            end

            # Count how many stars there are in this snapshot
            star_number = countStars(snapshot_path)

            # Check if there is subfind information in the group catalog file
            subfind_active = isSubfindActive(groupcat_path)

            if subfind_active

                # Read the number of halos
                n_groups_total = readGroupCatHeader(groupcat_path).n_groups_total

                # Read the group catalog file
                gc_data = readGroupCatalog(
                    groupcat_path,
                    snapshot_path,
                    Dict(:subhalo => ["S_LenType", "S_CM", "S_Pos"]),
                )

                # Load the group catalog data for the main subhalo
                s_cm       = gc_data[:subhalo]["S_CM"][:, 1]
                s_pos      = gc_data[:subhalo]["S_Pos"][:, 1]
                s_len_type = gc_data[:subhalo]["S_LenType"][:, 1]

                separation = sqrt(sum((s_cm - s_pos) .^ 2))
                str_scm    = round.(ustrip.(u"Mpc", s_cm), sigdigits=4)
                str_spos   = round.(ustrip.(u"Mpc", s_pos), sigdigits=4)

                s_len_stellar = s_len_type[PARTICLE_INDEX[:stellar] + 1]

                if s_len_stellar >= 1

                    real_stars        = findRealStars(snapshot_path)
                    stellar_n_subhalo = count(real_stars[1:s_len_stellar])

                else

                    stellar_n_subhalo = 0

                end

            end

            function subhaloReport(header::String)::Bool

                println(file, "#"^100)
                println(file, header)
                println(file, "#"^100)

                println(file, "\nSnapshot:               $(basename(snapshot_path))")

                println(file, "Physical time:          $(physical_time) Gyr")
                if cosmological
                    # For cosmological simulations print the scale factor and the redshift
                    println(file, "Scale factor:           $(scale_factor)")
                    println(file, "Redshift:               $(redshift)")
                end

                println(file, "Total number of stars:  $(star_number)")

                if subfind_active

                    println(file, "Number of halos:        $(n_groups_total)")

                    println(file, "\nMain subhalo properties:")

                    println(file, "\n\tNumber of stars:\n\n\t\t$(stellar_n_subhalo)\n")

                    println(file, "\tCenter of mass:\n\n\t\t$(str_scm) $(u"Mpc")\n")
                    println(
                        file,
                        "\tPosition of the particle with the minimum gravitational potential energy: \
                        \n\n\t\t$(str_spos) $(u"Mpc")\n",
                    )
                    println(
                        file,
                        "\tSeparation between the minimum potential and the center of mass: \
                        \n\n\t\t$(round(separation, sigdigits=6))\n",
                    )

                    println(file, "\t", "#"^62)
                    println(file, "\tNOTE: The number of stellar particles includes wind particles")
                    println(file, "\t", "#"^62)

                    println(file, "\n\tCell/particle number:\n")

                    for (i, len) in pairs(s_len_type)

                        component = PARTICLE_NAMES[INDEX_PARTICLE[i - 1]]
                        println(file, "\t\t$(component):$(" "^(21 - length(component))) $(len)")

                    end

                else
                    println(file, "\n" * "#"^51)
                    println(file, "There is no subfind information for this snapshot!")
                    println(file, "#"^51 * "\n")
                end

                println(file)

                return true

            end

            ##############
            # First stars
            ##############

            if star_number > 0 && !first_star_flag

                first_star_flag = subhaloReport("First snapshot with star formation:")

            end

            #################
            # First subhalos
            #################

            if subfind_active && !first_subhalo_flag

                first_subhalo_flag = subhaloReport("First snapshot with subfind information:")

            end

            ##################################
            # First stars in the main subhalo
            ##################################

            if subfind_active && !first_star_in_subhalo_flag && stellar_n_subhalo > 0

                first_star_in_subhalo_flag = subhaloReport(
                    "First snapshot with star formation in the main subhalo:"
                )

            end

            # End the loop when all the conditions have been met
            if first_star_flag && first_subhalo_flag && first_star_in_subhalo_flag
                break
            end

        end

        close(file)

    end

    return nothing

end

"""
    snapshotReport(
        simulation_paths::Vector{String},
        slices::Vector{Int};
        <keyword arguments>
    )::Nothing

Write a text file with information about a given snapshot.

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`. One text file will be written for each simulation.
  - `slices::Vector{Int}`: Selects which snapshots to plot for each simulation, starts at 1 and is independent of the number in the file name. If every snapshot is present, the relation is `slice` = (number in filename) + 1.
  - `output_path::String="."`: Path to the output folder.
  - `trans_mode::Union{Symbol,Tuple{TranslationType,RotationType,Dict{Symbol,Vector{String}}}}=:all_box`: How to translate and rotate the cells/particles, before filtering with `filter_mode`. For options see [`selectTransformation`](@ref).
  - `filter_mode::Union{Symbol,Tuple{Function,Dict{Symbol,Vector{String}}}}=:all`: Which cells/particles will be considered in the "filtered" section of the report. For options see [`selectFilter`](@ref).
"""
function snapshotReport(
    simulation_paths::Vector{String},
    slices::Vector{Int};
    output_path::String=".",
    trans_mode::Union{Symbol,Tuple{TranslationType,RotationType,Dict{Symbol,Vector{String}}}}=:all_box,
    filter_mode::Union{Symbol,Tuple{Function,Dict{Symbol,Vector{String}}}}=:all,
)::Nothing

    for (i, simulation_path) in pairs(simulation_paths)

        ############################################################################################
        # Load the relevant values and check for missing files
        ############################################################################################

        # Make a dataframe for the simulation with the following columns:
        #  - DataFrame index         -> :row_id
        #  - Number in the file name -> :numbers
        #  - Scale factor            -> :scale_factors
        #  - Redshift                -> :redshifts
        #  - Physical time           -> :physical_times
        #  - Lookback time           -> :lookback_times
        #  - Snapshot path           -> :snapshot_paths
        #  - Group catalog path      -> :groupcat_paths
        simulation_table = makeSimulationTable(simulation_path)

        # Select the slice for the `i`-th simulation
        slice = ring(slices, i)

        # Get the number in the filename
        snap_number = safeSelect(simulation_table[!, :numbers], slice)

        # Check that the snapshot for `slice` exist
        (
            !isempty(snap_number) ||
            throw(ArgumentError("snapshotReport: There are no snapshots for the slice \
            $(slice), for the simulation $(simulation_path)"))
        )

        # Find the target row
        snapshot_row = simulation_table[slice, :]

        # Read the path to the target snapshot
        snapshot_path = snapshot_row[:snapshot_paths]

        # Check if the snapshot path is missing
        snapshot_filename = "$(SNAP_BASENAME)_$(lpad(snap_number, 3, "0"))"
        (
            !ismissing(snapshot_path) ||
            throw(ArgumentError("snapshotReport: The snapshot $(snapshot_filename).hdf5 is \
            missing, for the simulation $(simulation_path)"))
        )

        # Read the path to the target group catalog file
        groupcat_path = snapshot_row[:groupcat_paths]

        # Check if the simulation is cosmological
        cosmological = isSnapCosmological(snapshot_path)

        # Read the physical time since the Big Bang
        physical_time = round(ustrip(u"Gyr", snapshot_row[:physical_times]), digits=2)

        # Compute the number of snapshots in the folder
        snapshot_length = count(!ismissing, simulation_table[!, :snapshot_paths])

        # Compute the number of group catalog files in the folder
        groupcat_length = count(!ismissing, simulation_table[!, :groupcat_paths])

        ############################################################################################
        # Read the data in the snapshot and group catalog file
        ############################################################################################

        # Select one snapshot file
        if isfile(snapshot_path)
            file_path = snapshot_path
        else
            file_path = minimum(glob("$(SNAP_BASENAME)_*.*.hdf5", snapshot_path))
        end

        physical_components = [:gas, :dark_matter, :stellar, :black_hole]

        # Detect which of the main physical components are present in the snapshot
        component_list = h5open(file_path, "r") do snapshot

            cp_types_in_snap = [get(PARTICLE_TYPE, key, nothing) for key in keys(snapshot)]

            filter!(in(physical_components), cp_types_in_snap)

            cp_types_in_snap

        end

        # Select the filter function, translation, rotation, and request dictionary
        base_request = mergeRequests(
            Dict(component => ["POS ", "MASS", "VEL "] for component in component_list),
            Dict(
                :gas => [
                    "NHP ",
                    "NH  ",
                    "PRES",
                    "GZ  ",
                    "RHO ",
                    "ID  ",
                    "ODIT",
                    "PARU",
                    "PARL",
                    "TAUS",
                    "RHOC",
                    "PARZ",
                    "PARH",
                    "ETAD",
                    "ETAI",
                    "PARR",
                    "PAZN",
                    "COLF",
                    "FRAC",
                    "SFFL",
                ],
                :stellar => [
                    "GAGE",
                    "ODIT",
                    "PARU",
                    "PARL",
                    "TAUS",
                    "RHOC",
                    "PARZ",
                    "PARH",
                    "ETAD",
                    "ETAI",
                    "PARR",
                    "PAZN",
                    "COLF",
                    "GMAS",
                    "GSFR",
                    "GPRE",
                    "ID  ",
                ],
            ),
        )

        translation, rotation, trans_request = selectTransformation(trans_mode, base_request)
        filter_function, request = selectFilter(filter_mode, trans_request)

        # Create the data dictionary
        data_dict = makeDataDict(simulation_path, slice, request, simulation_table)

        ############################################################################################
        # Global properties function
        ############################################################################################

        function printGlobalProperties(dd::Dict, box_string::String)::Nothing

            println(file, "\nGlobal properties ($(box_string)):")

            ##############################################################
            # Print the total number of cells/particles in each component
            ##############################################################

            println(file, "\n\tCell/particle number ($(box_string)):\n")

            total_count = 0

            for component in component_list

                count       = length(dd[component]["MASS"])
                title       = rpad("$(PARTICLE_NAMES[component]):", 25)
                total_count += count

                println(file, "\t\t$(title)$(count)")

            end

            println(file, "\n\t\tTotal count:             $(total_count)\n")

            #########################################
            # Print the total mass of each component
            #########################################

            println(file, "\tMasses ($(box_string)):\n")

            total_mass = 0.0u"Msun"

            for component in component_list

                mass       = integrateQty(dd, Symbol(component, :_mass))
                total_mass += mass
                title      = rpad("$(PARTICLE_NAMES[component]):", 25)

                println(file, "\t\t$(title)$(round(mass, sigdigits=3))")

            end

            str_total_mass = round(total_mass, sigdigits=3)

            println(file, "\n\t\tTotal mass:              $(str_total_mass)\n")

            if :gas ∈ component_list

                #######################################
                # Print the mass of each gas component
                #######################################

                println(file, "\tGas component masses ($(box_string)):\n")

                if isSimSFM(simulation_path)

                    gas_components = [
                        :ode_ionized,
                        :ode_atomic,
                        :ode_molecular,
                        :ode_stellar,
                        :ode_metals,
                        :ode_dust,
                        :ode_molecular_stellar,
                        :ode_neutral,
                        :ode_cold,
                    ]

                    gas_labels = [
                        "ODE ionized",
                        "ODE atomic",
                        "ODE molecular",
                        "ODE stellar",
                        "ODE metals",
                        "ODE dust",
                        "ODE molecular + stellar",
                        "ODE neutral",
                        "ODE cold",
                    ]

                    r_pad = 30

                else

                    gas_components = [:ionized, :br_atomic, :br_molecular, :neutral, :Z_gas]

                    gas_labels = ["Ionized", "BR atomic", "BR molecular", "Neutral", "Metals"]

                    r_pad = 25

                end

                gas_mass = integrateQty(dd, :gas_mass)

                for (gas_component, gas_label) in zip(gas_components, gas_labels)

                    mass     = integrateQty(dd, Symbol(gas_component, :_mass))
                    percent  = round((mass / gas_mass) * 100, sigdigits=4)
                    str_mass = round(mass, sigdigits=3)
                    title    = rpad("$(gas_label) mass:", r_pad)

                    println(file, "\t\t$(title)$(str_mass ) ($(percent)% of the total gas mass)")

                end

                ###################################################################
                # Print the fraction of gas cells that have entered the SF routine
                ###################################################################

                sffl = dd[:gas]["SFFL"]

                if !isempty(sffl)

                    println(
                        file,
                        "\n\tFraction of gas cells that have entered the SF routine ($(box_string)):\n"
                    )

                    gas_masses = scatterQty(dd, :gas_mass)

                    total_number       = length(gas_masses)
                    stellar_gas_number = count(isone, sffl)
                    fraction           = (stellar_gas_number / total_number) * 100

                    stellar_gas_mass = sum(gas_masses[Bool[sffl...]])
                    mass_fraction    = (stellar_gas_mass / sum(gas_masses)) * 100

                    println(file, "\t\t$(round(fraction, sigdigits=3))% of the cells")
                    println(file, "\t\t$(round(mass_fraction, sigdigits=3))% of the mass\n")

                end

            end

            ##################
            # ODEs parameters
            ##################

            sfm_quantities = [
                :integration_time,
                :parameter_uvb,
                :parameter_lwb,
                :tau_s,
                :parameter_cell_density,
                :parameter_metallicity,
                :parameter_column_height,
                :parameter_eta_d,
                :parameter_eta_i,
                :parameter_r,
                :parameter_zsn,
                :gas_mass,
                :gas_sfr,
                :gas_pressure,
            ]

            sfm_labels = [
                "Integration time",
                "UVB",
                "LWB",
                "τ_S",
                "Cell density",
                "Initial metallicity",
                "Column height",
                "η_diss",
                "η_ion",
                "Mass recycling parameter (R)",
                "Metallicity of the supernova ejecta (Zsn)",
                "Parent gas cell mass",
                "Parent gas cell SFR",
                "Parent gas cell pressure",
            ]

            sfm_units = [
                u"Myr",
                u"Myr^-1",
                u"Myr^-1",
                u"Myr",
                u"cm^-3",
                Unitful.NoUnits,
                u"pc",
                Unitful.NoUnits,
                Unitful.NoUnits,
                Unitful.NoUnits,
                Unitful.NoUnits,
                u"Msun",
                u"Msun * yr^-1",
                u"dyn * cm^-2",
            ]

            ##############################################################
            # Print the ODEs parameters for the gas that has formed stars
            ##############################################################

            sfm_blocks = [get(SFM_KEYS, sfm_quantity, nothing) for sfm_quantity in sfm_quantities]

            qty_present = isBlockPresent.(:stellar, sfm_blocks, snapshot_path)

            if any(qty_present)

                println(file, "\tProperties of the gas that has formed stars ($(box_string)):\n")

                iterator = zip(sfm_quantities, sfm_labels, sfm_units)

                for (quantity, label, unit) in collect(iterator)[qty_present]

                    qty_symbol = Symbol(:ode_stellar_, quantity)

                    values = filter(!isnan, ustrip.(unit, scatterQty(dd, qty_symbol)))

                    isempty(values) && continue

                    if quantity == :parameter_metallicity
                        values = values ./ SOLAR_METALLICITY
                        unit = "Z⊙"
                    end

                    qty_p50 = round(median(values), sigdigits=4)
                    qty_min = round(minimum(values), sigdigits=4)
                    qty_max = round(maximum(values), sigdigits=4)

                    println(file, "\t\t$(label):\n")
                    println(file, "\t\t\tMedian:  $(qty_p50) $(unit)")
                    println(file, "\t\t\tMinimum: $(qty_min) $(unit)")
                    println(file, "\t\t\tMaximum: $(qty_max) $(unit)\n")

                end

            end

            ##############################################
            # Print the ODEs parameters for the gas cells
            ##############################################

            gas_sfm_idxs = 1:11

            qty_present = isBlockPresent.(:gas, sfm_blocks[gas_sfm_idxs], snapshot_path)

            if any(qty_present)

                println(file, "\tODE parameters for the gas cells ($(box_string)):\n")

                iterator = zip(sfm_quantities, sfm_labels, sfm_units)

                for (quantity, label, unit) in collect(iterator)[gas_sfm_idxs][qty_present]

                    qty_symbol = Symbol(:ode_gas_, quantity)

                    values = filter(!isnan, ustrip.(unit, scatterQty(dd, qty_symbol)))

                    isempty(values) && continue

                    if quantity == :parameter_metallicity
                        values = values ./ SOLAR_METALLICITY
                        unit = "Z⊙"
                    end

                    qty_p50 = round(median(values), sigdigits=4)
                    qty_min = round(minimum(values), sigdigits=4)
                    qty_max = round(maximum(values), sigdigits=4)

                    println(file, "\t\t$(label):\n")
                    println(file, "\t\t\tMedian:  $(qty_p50) $(unit)")
                    println(file, "\t\t\tMinimum: $(qty_min) $(unit)")
                    println(file, "\t\t\tMaximum: $(qty_max) $(unit)\n")

                end

            end

            return nothing

        end

        ############################################################################################
        # Print the report header
        ############################################################################################

        # Create the output file
        filename = "$(snapshot_filename)_of_$(basename(simulation_path))"
        file = open(joinpath(mkpath(output_path), "report_for_$(filename).txt"), "w")

        println(file, "#"^100)
        println(file, "\nSimulation name:     $(basename(simulation_path))")

        println(file, "Snapshot path:       $(snapshot_path)")
        println(file, "Snapshot number:     $(slice) (of $(snapshot_length))\n")

        if isSubfindActive(groupcat_path)

            n_groups_total = readGroupCatHeader(groupcat_path).n_groups_total
            gc_data        = readGroupCatalog(groupcat_path, snapshot_path, Dict(:group => ["G_Nsubs"]))
            n_subfinds     = gc_data[:group]["G_Nsubs"][1]

            println(file, "Subfind path:        $(groupcat_path)")
            println(file, "Subfind number:      $(slice) (of $(groupcat_length))\n")

            println(file, "Number of halos:     $(n_groups_total)")
            println(file, "Number of subhalos:  $(n_subfinds) (in the main halo)\n")

        end

        println(file, "#"^100)
        println(file)

        if cosmological
            # For cosmological simulations print the scale factor and the redshift
            scale_factor = round(snapshot_row[:scale_factors], digits=3)
            redshift     = round(snapshot_row[:redshifts], digits=3)
            println(file, "Cosmological:        Yes")
            println(file, "Scale factor:        $(scale_factor)")
            println(file, "Redshift:            $(redshift)")
        else
            println(file, "Cosmological:        No")
        end

        println(file, "Physical time:       $(physical_time) Gyr")

        if !PHYSICAL_UNITS && cosmological
            println(file, "Lenght units:        Comoving\n")
        else
            println(file, "Lenght units:        Physical\n")
        end

        ######################################
        # Print the fraction of in-situ stars
        ######################################

        if snapshot_length >= 2 && isSubfindActive(groupcat_path)

            insitu_idx = filterByBirthPlace(data_dict, :exsitu)[:stellar]

            iMs = sum(data_dict[:stellar]["MASS"][insitu_idx]; init=0.0u"Msun")
            tMs = sum(data_dict[:stellar]["MASS"]; init=0.0u"Msun")
            insitu_fraction = round(uconvert.(Unitful.NoUnits, (iMs / tMs) * 100); sigdigits=2)


            println(file, "#"^100)

            println(
                file,
                "\nFraction of in-situ stars in the main subhalo:\t$(insitu_fraction)%\n"
            )

        end

        ############################################################################################
        # Print the global properties of the full simulation box
        ############################################################################################

        println(file, "#"^100)
        println(file, "Full simulation box")
        println(file, "#"^100)

        printGlobalProperties(data_dict, "full simulation box")

        ############################################################################################
        # Print the global properties of the filtered simulation box
        ############################################################################################

        translateData!(data_dict, translation...)
        rotateData!(data_dict, rotation...)
        filterData!(data_dict; filter_function)

        println(file, "#"^100)
        println(file, "Filtered box with:\n")

        println(file, "\tTranslation: $(translation)")
        println(file, "\tRotation:    $(rotation[1:2])")

        if filter_mode isa Symbol
            println(file, "\tFilter mode: :$(filter_mode)")
        else
            println(file, "\tFilter function: $(String(Symbol(filter_function)))")
        end

        println(file, "#"^100)

        printGlobalProperties(data_dict, "filtered box")

        ############################################################################################
        # Print some characteristic radii
        ############################################################################################

        if :gas in component_list

            if isSimSFM(simulation_path)

                components = [
                    :stellar,
                    :gas,
                    :ode_ionized,
                    :ode_atomic,
                    :ode_molecular,
                    :ode_stellar,
                    :ode_metals,
                    :ode_dust,
                    :ode_molecular_stellar,
                    :ode_neutral,
                    :ode_cold,
                ]

                labels = [
                    "stellar",
                    "total gas",
                    "ODE ionized gas",
                    "ODE atomic gas",
                    "ODE molecular gas",
                    "ODE stellar gas",
                    "ODE metals",
                    "ODE dust",
                    "ODE molecular + stellar",
                    "ODE neutral",
                    "ODE cold",
                ]

            else

                components = [:stellar, :gas, :ionized, :br_atomic, :br_molecular, :neutral, :Z_gas]

                labels = [
                    "stellar",
                    "total gas",
                    "ionized gas",
                    "BR atomic gas",
                    "BR molecular gas",
                    "neutral gas",
                    "metals",
                ]

            end

            disc_idxs = filterBySphere(data_dict, 0.0u"kpc", DISK_R, :zero)

            println(file, "Characteristic radii (filtered box):\n")

            for (component, label) in zip(components, labels)

                masses = scatterQty(data_dict, Symbol(component, :_mass))

                if !isempty(masses)

                    if component == :stellar
                        mass_inside = masses[disc_idxs[:stellar]]
                        pos_inside  = data_dict[:stellar]["POS "][:, disc_idxs[:stellar]]
                    else
                        mass_inside = masses[disc_idxs[:gas]]
                        pos_inside  = data_dict[:gas]["POS "][:, disc_idxs[:gas]]
                    end

                    mass_radius_90 = computeMassRadius(pos_inside, mass_inside)
                    mass_radius_95 = computeMassRadius(pos_inside, mass_inside; percent=95.0)

                    mass_radius_90_str = round(ustrip(u"kpc", mass_radius_90), sigdigits=4)
                    mass_radius_95_str = round(ustrip(u"kpc", mass_radius_95), sigdigits=4)

                    println(file, "\tRadius containing X% of the $(label) mass (r < $(DISK_R)):\n")
                    println(file, "\t\t$(mass_radius_90_str) $(u"kpc") (90%)")
                    println(file, "\t\t$(mass_radius_95_str) $(u"kpc") (95%)\n")

                end

            end

        end

        close(file)

    end

    return nothing

end

"""
    quantityReport(
        simulation_paths::Vector{String},
        quantity::Symbol;
        <keyword arguments>
    )::Nothing

Write a text file with information about a given `quantity`.

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`. One text file will be written for each simulation.
  - `quantity::Symbol`: Target quantity. Has to be one of the valid quantities of [`scatterQty`](@ref).
  - `slice::IndexType=(:)`: Slice of the simulation, i.e. which snapshots will be analysed. It can be an integer (a single snapshot), a vector of integers (several snapshots), an `UnitRange` (e.g. 5:13), an `StepRange` (e.g. 5:2:13) or (:) (all snapshots). Starts at 1 and out of bounds indices are ignored.
  - `output_path::String="."`: Path to the output folder.
  - `trans_mode::Union{Symbol,Tuple{TranslationType,RotationType,Dict{Symbol,Vector{String}}}}=:all_box`: How to translate and rotate the cells/particles, before filtering with `filter_mode`. For options see [`selectTransformation`](@ref).
  - `filter_mode::Union{Symbol,Tuple{Function,Dict{Symbol,Vector{String}}}}=:all`: Which cells/particles will be selected. For options see [`selectFilter`](@ref).
"""
function quantityReport(
    simulation_paths::Vector{String},
    quantity::Symbol;
    slice::IndexType=(:),
    output_path::String=".",
    trans_mode::Union{Symbol,Tuple{TranslationType,RotationType,Dict{Symbol,Vector{String}}}}=:all_box,
    filter_mode::Union{Symbol,Tuple{Function,Dict{Symbol,Vector{String}}}}=:all,
)::Nothing

    plot_params = plotParams(quantity)

    da_function = dd -> uconvert.(plot_params.unit, scatterQty(dd, quantity))

    return quantityReport(
        simulation_paths,
        plot_params.request;
        da_function,
        slice,
        output_path,
        label=string(quantity),
        trans_mode,
        filter_mode,
    )

end

"""
    quantityReport(
        simulation_paths::Vector{String},
        base_request::Dict{Symbol,Vector{String}};
        <keyword arguments>
    )::Nothing

Write a text file with information about a the results of applying `da_function` to each snapshot.

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`. One text file will be written for each simulation.
  - `base_request::Dict{Symbol,Vector{String}}`: Dictionary with the shape `cell/particle type` -> [`block`, `block`, ...], where the possible types are the keys of [`PARTICLE_INDEX`](@ref), and the possible quantities are the keys of [`QUANTITIES`](@ref). Which data blocks are needed depends on `da_function`.
  - `da_function::Function=getNothing`: Data analysis function. It must have the same signature as the examples in `./src/analysis/data_analysis.jl`, but it should only return an Array.
  - `da_arg::Tuple=()`: Psitional arguments for the data analysis function.
  - `da_kwarg::NamedTuple=(;)`: Keyword arguments for the data analysis function.
  - `slice::IndexType=(:)`: Slice of the simulation, i.e. which snapshots will be analysed. It can be an integer (a single snapshot), a vector of integers (several snapshots), an `UnitRange` (e.g. 5:13), an `StepRange` (e.g. 5:2:13) or (:) (all snapshots). Starts at 1 and out of bounds indices are ignored.
  - `output_path::String="."`: Path to the output folder.
  - `label::Union{String,Nothing}=nothing``: Label for the output filename. If set to `nothing`, it will be `"quantity"`.
  - `trans_mode::Union{Symbol,Tuple{TranslationType,RotationType,Dict{Symbol,Vector{String}}}}=:all_box`: How to translate and rotate the cells/particles, before filtering with `filter_mode`. For options see [`selectTransformation`](@ref).
  - `filter_mode::Union{Symbol,Tuple{Function,Dict{Symbol,Vector{String}}}}=:all`: Which cells/particles will be selected. For options see [`selectFilter`](@ref).
  - `show_progress::Bool=true`: If a progress bar will be shown.
"""
function quantityReport(
    simulation_paths::Vector{String},
    base_request::Dict{Symbol,Vector{String}};
    da_function::Function=getNothing,
    da_arg::Tuple=(),
    da_kwarg::NamedTuple=(;),
    slice::IndexType=(:),
    output_path::String=".",
    label::Union{String,Nothing}=nothing,
    trans_mode::Union{Symbol,Tuple{TranslationType,RotationType,Dict{Symbol,Vector{String}}}}=:all_box,
    filter_mode::Union{Symbol,Tuple{Function,Dict{Symbol,Vector{String}}}}=:all,
    show_progress::Bool=true,
)::Nothing

    for simulation_path in simulation_paths

        simulation_name = basename(simulation_path)
        cosmological    = isSimCosmological(simulation_path)

        ############################################################################################
        # Load the relevant values and check for missing files
        ############################################################################################

        # Make a dataframe for the simulation with the following columns:
        #  - DataFrame index         -> :row_id
        #  - Number in the file name -> :numbers
        #  - Scale factor            -> :scale_factors
        #  - Redshift                -> :redshifts
        #  - Physical time           -> :physical_times
        #  - Lookback time           -> :lookback_times
        #  - Snapshot path           -> :snapshot_paths
        #  - Group catalog path      -> :groupcat_paths
        simulation_table = makeSimulationTable(simulation_path)

        # Compute the number of snapshots in the folder
        snapshot_n = count(!ismissing, simulation_table[!, :snapshot_paths])

        # Check that there is at least one snapshot
        (
            !iszero(snapshot_n) ||
            throw(ArgumentError("quantityReport: There are no snapshots in $(simulation_path)"))
        )

        # Compute the number of group catalog files in the folder
        groupcat_n = count(!ismissing, simulation_table[!, :groupcat_paths])

        # Get the snapshot paths
        snapshot_paths = safeSelect(simulation_table[!, :snapshot_paths], slice)

        # Check that the snapshots for `slice` exist
        (
            !isempty(snapshot_paths) ||
            throw(ArgumentError("quantityReport: There are no snapshots for the slice \
            $(slice), for the simulation $(simulation_path)"))
        )

        # Get the snapshot paths
        global_indices = safeSelect(simulation_table[!, :row_id], slice)

        # Get the number in the filename
        snapshot_numbers = safeSelect(simulation_table[!, :numbers], slice)

        # Get the times and redshifts
        physical_times = safeSelect(simulation_table[!, :physical_times], slice)
        redshifts      = safeSelect(simulation_table[!, :redshifts], slice)

        if isnothing(label)
            qty_str = "quantity"
        else
            qty_str = label
        end

        # Create the output file
        file = open(
            joinpath(mkpath(output_path),
            "$(qty_str)_report_for_$(simulation_name).txt"),
            "w",
        )

        # Select the filter function, translation, rotation, and request dictionary
        translation, rotation, trans_request = selectTransformation(trans_mode, base_request)
        filter_function, request = selectFilter(filter_mode, trans_request)

        ############################################################################################
        # Print the header
        ############################################################################################

        println(file, "#"^100)
        println(file, "\nSimulation name:  $(simulation_name)")

        println(file, "\nNumber of snapshots:       $(snapshot_n)")
        println(file, "Number of group catalogs:  $(groupcat_n)\n")

        println(file, "#"^100)

        println(file, "\nTranslation:   $(translation)")
        println(file, "Rotation:      $(rotation[1:2])")

        if filter_mode isa Symbol
            println(file, "Filter mode:   :$(filter_mode)")
        else
            println(file, "Filter function:  $(String(Symbol(filter_function)))")
        end

        if !PHYSICAL_UNITS && cosmological
            println(file, "Lenght units:  Comoving\n")
        else
            println(file, "Lenght units:  Physical\n")
        end

        println(file, "#"^100)

        println(file)

        ############################################################################################
        # Print statistics of `quantity` for every snapshot
        ############################################################################################

        iterator = zip(snapshot_paths, global_indices, snapshot_numbers, physical_times, redshifts)

        prog_bar = Progress(
            length(iterator),
            dt=0.5,
            desc="Writing the quantity report for simulation $(simulation_name)... ",
            color=:blue,
            barglyphs=BarGlyphs("|#  |"),
            enabled=show_progress,
        )

        for (snapshot_path, global_index, snapshot_number, physical_time, redshift) in iterator

            if ismissing(snapshot_path)
                println(file, "Snapshot:  snap_$(lpad(snapshot_number, 3, '0')).hdf5 is missing!\n")
                next!(prog_bar)
                continue
            else
                println(file, "Snapshot:       $(snapshot_path)")
                println(file, "Physical time:  $(round(ustrip(u"Gyr", physical_time); sigdigits=3)) Gyr")
                println(file, "Redshift:       $(round(redshift; sigdigits=3))\n")
            end

            # Create the data dictionary
            data_dict = makeDataDict(simulation_path, global_index, request, simulation_table)

            # Translate, rotate, and filter the data
            translateData!(data_dict, translation...)
            rotateData!(data_dict, rotation...)
            filterData!(data_dict; filter_function)

            # Apply the data analysis function
            quantity_values = da_function(data_dict; da_arg..., da_kwarg...)

            if isnothing(quantity_values)
                println(file, "\tThe data analysis function returned `nothing`!\n")
                next!(prog_bar)
                continue
            end

            (
                isa(quantity_values, AbstractArray) ||
                throw(ArgumentError("quantityReport: `da_function` must return an array, but I got \
                $(typeof(quantity_values)) instead."))
            )

            # Ignore NaN and Inf values
            filter!(!isnan, quantity_values)
            filter!(!isinf, quantity_values)

            if !isempty(quantity_values)

                println(file, "\tNumber of particles/cells:  $(length(quantity_values))\n")

                mean_qty  = mean(quantity_values)
                std_qty   = std(quantity_values)
                extre_qty = extrema(quantity_values)

                println(file, "\tMean:             $(round(mean_qty; digits=2))")
                println(file, "\tSTD:              $(round(std_qty; digits=2))")
                println(file, "\tMinimum:          $(round(extre_qty[1]; digits=2))")
                println(file, "\tMaximum:          $(round(extre_qty[2]; digits=2))\n")

                for p in [10, 25, 50, 75, 90]
                    percentile_qty = percentile(quantity_values, p)
                    println(file, "\tPercentile $(p)%:   $(round(percentile_qty; digits=2))")
                end

            else

                println(file, "\tThere are no valid data point in this snapshot!")

            end

            println(file)

            next!(prog_bar)

        end

        close(file)

    end

    return nothing

end
