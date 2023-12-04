####################################################################################################
# Pipeline functions.
####################################################################################################

"""
    snapshotPlot(
        simulation_paths::Vector{String},
        request::Dict{Symbol,Vector{String}},
        plot_functions::Vector{<:Function};
        <keyword arguments>
    )::Nothing

Generate one figure per snapshot, for one or more simulations.

Some of the features are:

  - It can produce scatter plots, line plots, histograms, and heatmaps.
  - It can generate an animation of the results.
  - It transparently manages units; you only have to indicate the final unit of each axis.

!!! note

    The snapshots of different simulations are grouped by the number in the file names, regardless of the "Time" parameter in the header.
    The data from the longest running simulation is used for the time stamp in the automatic title.

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`.
  - `request::Dict{Symbol,Vector{String}}`: Dictionary with the shape `cell/particle type` -> [`block`, `block`, ...], where the possible types are the keys of [`PARTICLE_INDEX`](@ref), and the possible quantities are the keys of [`QUANTITIES`](@ref). Which data blocks are needed depends on the provided functions `da_functions`.
  - `plot_functions::Vector{<:Function}`: Vector of plotting functions from [Makie](https://docs.makie.org/stable/). This sets the type of plot for each simulation.
    The supported functions are:

      + `scatter!`      -> Scatter plot.
      + `lines!`        -> Line plot.
      + `scatterlines!` -> Scatter plot with lines between the markers.
      + `hist!`         -> Histogram.
      + `heatmap!`      -> Heatmap.
  - `pf_kwargs::Vector{<:NamedTuple}=[(;)]`: Vector of keyword arguments for the functions in `plot_functions`.

### snapshotPlot configuration

  - `output_path::String="./plots"`: Path to the output folder.
  - `base_filename::String="snapshot"`: Every file will be named `base_filename`-XXX`output_format` where XXX is the snapshot number.
  - `output_format::String=".png"`: File format for the figure. All formats supported by [Makie](https://docs.makie.org/stable/) can be used, namely `.pdf`, `.svg` and `.png`.
  - `warnings::Bool=true`: If a warning will be raised when there is some problem, but the function can still run using sane defaults.
  - `show_progress::Bool=true`: If a progress bar will be shown.

### Data manipulation options

  - `slice::IndexType=(:)`: Slice of the simulation, i.e. which snapshots will be read. It can be an integer (a single snapshot), a vector of integers (several snapshots), an `UnitRange` (e.g. 5:13), an `StepRange` (e.g. 5:2:13) or (:) (all snapshots). It works over the longest possible list of snapshots among the simulations (grouped by the number in the file names). Out of bounds indices are ignored.
  - `filter_function::Function=filterNothing`: A functions with the signature:

    `filter_function(data_dict) -> indices`

    where

      + `data_dict::Dict`: A dictionary with the following shape:

          * `:sim_data`          -> ::Simulation (see [`Simulation`](@ref)).
          * `:snap_data`         -> ::Snapshot (see [`Snapshot`](@ref)).
          * `:gc_data`           -> ::GroupCatalog (see [`GroupCatalog`](@ref)).
          * `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
          * `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
          * `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
          * ...
          * `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
          * `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
          * `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
          * ...

      + `indices::Dict`: A dictionary with the following shape:

          * `cell/particle type` -> idxs::IndexType
          * `cell/particle type` -> idxs::IndexType
          * `cell/particle type` -> idxs::IndexType
          * ...
  - `da_functions::Vector{<:Function}=[getNothing]`: Vector of data analysis functions. See the required signature and examples in `./src/data_analysis.jl`.
  - `da_args::Vector{<:Tuple}=[()]`: Vector of positional arguments for the data analysis functions.
  - `da_kwargs::Vector{<:NamedTuple}=[(;)]`: Vector of keyword arguments for the data analysis functions.
  - `post_processing::Function=getNothing`: Post processing function. See the required signature and examples in `./src/post_processing.jl`.
  - `pp_args::Tuple=()`: Positional arguments for the post processing function.
  - `pp_kwargs::NamedTuple=(;)`: Keyword arguments for the post processing function.
  - `transform_box::Bool=false`: If a translation and rotation (in that order) will be applied to the simulation box, affecting the positions and velocities of all the cells and particles. If active, it is applied AFTER the `filter_function`.
  - `translation::Union{Symbol,NTuple{2,Int}}=:zero`: Type of translation (only relevant if `transform_box` = true). The options are:

      + `:zero`                       -> No translation is applied.
      + `:global_cm`                  -> Sets the center of mass of the whole system (after filtering) as the new origin.
      + `:stellar_cm`                 -> Sets the stellar center of mass (after filtering) as the new origin.
      + `(halo_idx, subhalo_rel_idx)` -> Sets the position of the potencial minimum for the `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo), as the new origin.
      + `(halo_idx, 0)`               -> Sets the center of mass of the `halo_idx::Int` halo, as the new origin.
  - `rotation::Symbol=:zero`: Type of rotation (only relevant if `transform_box` = true). The options are:

      + `:zero`       -> No rotation is appplied.
      + `:global_am`  -> Sets the angular momentum of the whole system as the new z axis.
      + `:stellar_am` -> Sets the stellar angular momentum as the new z axis.
      + `:stellar_pa` -> Sets the stellar principal axis as the new coordinate system.
  - `smooth::Int=0`: The result of `da_functions` will be smooth out using `smooth` bins. Set it to 0 if you want no smoothing. Only valid for `scatter!`, `lines!`, and `scatterlines!` plots.
  - `x_unit::Unitful.Units=Unitful.NoUnits`: Target unit for the x axis. The values will be converted accordingly. Use the default value of `Unitful.NoUnits` for dimensionless quantities.
  - `y_unit::Unitful.Units=Unitful.NoUnits`: Target unit for the y axis. The values will be converted accordingly. Use the default value of `Unitful.NoUnits` for dimensionless quantities.
  - `x_exp_factor::Int=0`: Numerical exponent to scale down the x axis, e.g. if `x_exp_factor` = 10 the values will be divided by ``10^{10}``. The default is no scaling.
  - `y_exp_factor::Int=0`: Numerical exponent to scale down the y axis, e.g. if `y_exp_factor` = 10 the values will be divided by ``10^{10}``. The default is no scaling.
  - `x_trim::NTuple{2,<:Real}=(-Inf, Inf)`: The data will be trim down so the x coordinates fit within `x_trim`.
  - `y_trim::NTuple{2,<:Real}=(-Inf, Inf)`: The data will be trim down so the y coordinates fit within `y_trim`. This option does not affect histograms.
  - `x_edges::Bool=false`: Set it to `true` if you want to keep the borders of `x_trim`.
  - `y_edges::Bool=false`: Set it to `true` if you want to keep the borders of `y_trim`.
  - `x_func::Function=identity`: Function to be applied to the values of the x axis. It must be a pure function with the signature `x_func(x_values::Vector{Float64})::Vector{Float64}`. The output must have the same length as the input. This function will be applied regardless of units and possible domain problems (use `x_trim` to solve incompatibilities), and that it will not be reflected in the automatic labeling.
  - `y_func::Function=identity`: Function to be applied to the values of the y axis. It must be a pure function with the signature `y_func(y_values::Vector{Float64})::Vector{Float64}`. The output must have the same length as the input. This function will be applied regardless of units and possible domain problems (use `y_trim` to solve incompatibilities), and that it will not be reflected in the automatic labeling.

### Axes options

  - `xaxis_label::AbstractString="auto_label"`: Label for the x axis. It can contain the string `auto_label`, which will be replaced by: `xaxis_var_name` / 10^`x_exp_factor` `x_unit`.
  - `yaxis_label::AbstractString="auto_label"`: Label for the y axis. It can contain the string `auto_label`, which will be replaced by: `yaxis_var_name` / 10^`y_exp_factor` `y_unit`.
  - `xaxis_var_name::AbstractString=""`: Name of the variable for the x axis.
  - `yaxis_var_name::AbstractString=""`: Name of the variable for the y axis.
  - `xaxis_scale_func::Function=identity`: Scaling function for the x axis. The options are the scaling functions accepted by [Makie](https://docs.makie.org/stable/): log10, log2, log, sqrt, Makie.logit, Makie.Symlog10, Makie.pseudolog10, and identity.
  - `yaxis_scale_func::Function=identity`: Scaling function for the y axis. The options are the scaling functions accepted by [Makie](https://docs.makie.org/stable/): log10, log2, log, sqrt, Makie.logit, Makie.Symlog10, Makie.pseudolog10, and identity.
  - `xaxis_limits::Tuple{<:Union{Real,Nothing},<:Union{Real,Nothing}}=(nothing, nothing)`: Set it to a value different than `nothing` if you want to fix the limits of the x axis.
  - `yaxis_limits::Tuple{<:Union{Real,Nothing},<:Union{Real,Nothing}}=(nothing, nothing)`: Set it to a value different than `nothing` if you want to fix the limits of the y axis.

### Plotting options

  - `save_figures::Bool=true`: If every figure will be saved as an image.
  - `backup_results::Bool=false`: If the values to be plotted will be backup in a [JLD2](https://github.com/JuliaIO/JLD2.jl) file.
  - `sim_labels::Union{Vector{String},Nothing}=nothing`: Labels for the plot legend, one per simulation. Set it to `nothing` if you don't want a legend.
  - `title::Union{Symbol,<:AbstractString}=""`: Title for the figure. If left empty, no title is printed. It can also be set to one of the following options:

      + `:physical_time` -> Physical time since the Big Bang.
      + `:lookback_time` -> Physical time left to reach the last snapshot.
      + `:scale_factor`  -> Scale factor (only relevant for cosmological simulations).
      + `:redshift`      -> Redshift (only relevant for cosmological simulations).
  - `pt_per_unit::Float64=0.75`: Factor to scale up or down the size of the figures, keeping the proportions. It only works for `.pdf` and `.svg`.
  - `px_per_unit::Float64=1.0`: Factor to scale up or down the size of the figures, keeping the proportions. It only works for `.png`.
  - `size::NTuple{2,Int}=(1000, 750)`: Size of the figures in points. For PNGs, by default points = pixels (as given by `px_per_unit` = 1.0), and for PDFs and SVGs, points = 0.75 * pixels (as given by `pt_per_unit` = 0.75).
  - `aspect::Union{DataAspect,AxisAspect,Nothing}=nothing`: Aspect ratio of the figures. The options are:

      + `nothing`       -> Default, the aspect ratio will be chosen by [Makie](https://docs.makie.org/stable/).
      + `AxisAspect(n)` -> The aspect ratio will be given by the number `n` = width / height.
      + `DataAspect()`  -> The aspect ratio of the data will be used.
  - `series_colors::Union{Vector{<:ColorType},Nothing}=nothing`: Colors for the different simulations. If set to `nothing`, the colors will be assigned automatically. This is only relevant for `scatter!`, `scatterlines!`, and `lines!` plots.
  - `series_markers::Union{Vector{Symbol},Nothing}=nothing`: Markers for the different simulations. If set to `nothing`, the markers will be assigned automatically. This is only relevant for `scatter!` and `scatterlines!` plots.
  - `series_linestyles::Union{Vector{<:LineStyleType},Nothing}=nothing`: Line styles for the different simulations. If set to `nothing`, the line styles will be assigned automatically. This is only relevant for `lines!` and `scatterlines!` plots.

## Animation options

  - `animation::Bool=false`: If an animation will be created.
  - `animation_filename::String="animation.mp4"`: Filename for the animation, including its extension. All formats supported by [Makie](https://docs.makie.org/stable/) can be used, namely `.mkv`, `.mp4`, `.webm` and `.gif`.
  - `framerate::Int=15`: Frame rate of the animation.
"""
function snapshotPlot(
    simulation_paths::Vector{String},
    request::Dict{Symbol,Vector{String}},
    plot_functions::Vector{<:Function};
    pf_kwargs::Vector{<:NamedTuple}=[(;)],
    # `snapshotPlot` configuration
    output_path::String="./plots",
    base_filename::String="snapshot",
    output_format::String=".png",
    warnings::Bool=true,
    show_progress::Bool=true,
    # Data manipulation options
    slice::IndexType=(:),
    filter_function::Function=filterNothing,
    da_functions::Vector{<:Function}=[getNothing],
    da_args::Vector{<:Tuple}=[()],
    da_kwargs::Vector{<:NamedTuple}=[(;)],
    post_processing::Function=getNothing,
    pp_args::Tuple=(),
    pp_kwargs::NamedTuple=(;),
    transform_box::Bool=false,
    translation::Union{Symbol,NTuple{2,Int}}=:zero,
    rotation::Symbol=:zero,
    smooth::Int=0,
    x_unit::Unitful.Units=Unitful.NoUnits,
    y_unit::Unitful.Units=Unitful.NoUnits,
    x_exp_factor::Int=0,
    y_exp_factor::Int=0,
    x_trim::NTuple{2,<:Real}=(-Inf, Inf),
    y_trim::NTuple{2,<:Real}=(-Inf, Inf),
    x_edges::Bool=false,
    y_edges::Bool=false,
    x_func::Function=identity,
    y_func::Function=identity,
    # Axes options
    xaxis_label::AbstractString="auto_label",
    yaxis_label::AbstractString="auto_label",
    xaxis_var_name::AbstractString="",
    yaxis_var_name::AbstractString="",
    xaxis_scale_func::Function=identity,
    yaxis_scale_func::Function=identity,
    xaxis_limits::Tuple{<:Union{Real,Nothing},<:Union{Real,Nothing}}=(nothing, nothing),
    yaxis_limits::Tuple{<:Union{Real,Nothing},<:Union{Real,Nothing}}=(nothing, nothing),
    # Plotting options
    save_figures::Bool=true,
    backup_results::Bool=false,
    sim_labels::Union{Vector{String},Nothing}=nothing,
    title::Union{Symbol,<:AbstractString}="",
    pt_per_unit::Float64=0.75,
    px_per_unit::Float64=1.0,
    size::NTuple{2,Int}=(1000, 750),
    aspect::Union{DataAspect,AxisAspect,Nothing}=nothing,
    series_colors::Union{Vector{<:ColorType},Nothing}=nothing,
    series_markers::Union{Vector{Symbol},Nothing}=nothing,
    series_linestyles::Union{Vector{<:LineStyleType},Nothing}=nothing,
    # Animation options
    animation::Bool=false,
    animation_filename::String="animation.mp4",
    framerate::Int=15,
)::Nothing

    # Create the output folder if it doesn't exist
    mkpath(output_path)

    # Compute the number of simulations
    n_simulations = length(simulation_paths)

    # Make a dataframe for every simulation, with the following columns:
    #   - 1. DataFrame index
    #   - 2. Number in the file name
    #   - 3. Scale factor
    #   - 4. Redshift
    #   - 5. Physical time
    #   - 6. Lookback time
    #   - 7. Snapshot path
    #   - 8. Group catalog path
    simulation_tables = [makeSimulationTable(source; warnings) for source in simulation_paths]

    # Compute the different ways to index the snapshots
    snapshot_numbers = sort!(union([table[!, :numbers] for table in simulation_tables]...))
    global_indices = collect(eachindex(snapshot_numbers))
    slice_indices = safeSelect(global_indices, slice; warnings)

    # Compute the number of figures
    n_frames = length(slice_indices)

    # Check that after slicing there is at least one snapshot left
    (
        !iszero(n_frames) ||
        throw(ArgumentError("snapshotPlot: There are no snapshots left after slicing \
        with `slice` = $slice"))
    )

    ################################################################################################
    # Set up the canvas for the figures.
    ################################################################################################

    # Apply the global theme defined in `./src/constants.jl`
    set_theme!(THEME)

    # Create the figure
    figure = Figure(; size)

    # Create a list of distinguishable colors for the different simulations,
    # only relevant for `scatter!`, `lines!`, and `scatterlines!` plots
    if isnothing(series_colors)
        colors = distinguishable_colors(n_simulations, [RGB(1, 1, 1), RGB(0, 0, 0)], dropseed=true)
    else
        colors = series_colors
    end

    # Select a set of markers for the different simulations,
    # only relevant for `scatter!` and `scatterlines!` plots
    markers = isnothing(series_markers) ? MARKERS : series_markers

    # Select a set of line styles for the different simulations,
    # only relevant for `lines!` and `scatterlines!` plots
    linestyles = isnothing(series_linestyles) ? LINE_STYLES : series_linestyles

    # Create the labels
    xlabel = LaTeXString(
        replace(xaxis_label, "auto_label" => getLabel(xaxis_var_name, x_exp_factor, x_unit)),
    )
    ylabel = LaTeXString(
        replace(yaxis_label, "auto_label" => getLabel(yaxis_var_name, y_exp_factor, y_unit)),
    )

    # Create the axes
    axes = Makie.Axis(figure[1, 1]; limits=(xaxis_limits, yaxis_limits), xlabel, ylabel, aspect)

    ################################################################################################
    # Set up the animation.
    ################################################################################################

    if animation

        (
            n_frames >= framerate ||
            !warnings ||
            @warn("snapshotPlot: With `framerate` = $framerate and `slice` = $slice, \
            the animation is less than one second long")
        )

        # Initialize the animation stream
        vs = VideoStream(figure; framerate)

    end

    ################################################################################################
    # Main loop.
    ################################################################################################

    # Initialize the progress bar
    prog_bar = Progress(
        n_frames,
        dt=0.5,
        desc="Analyzing and plotting the data... ",
        color=:blue,
        barglyphs=BarGlyphs("|#  |"),
        enabled=show_progress,
    )

    # Flag to warn if nothing is plotted (because every snapshot was skipped)
    #   + `true`: At least one snapshot was plotted
    #   + `false`: Every snapshot was skipped
    plot_something = false

    # Loop through the different snapshots
    @inbounds for (slice_index, global_index) in pairs(slice_indices)

        # Flag to keep the x axis with a linear scale if there are no data points left after
        # trying to use a nonlinear scale
        #   + `true`: The x axis will use the scale given by `xaxis_scale_func`
        #   + `false`: The x axis will use a linear scale
        xscale_flag = true

        # Flag to keep the y axis with a linear scale if there are no data points left after
        # trying to use a nonlinear scale
        #   + `true`: The y axis will use the scale given by `yaxis_scale_func`
        #   + `false`: The y axis will use a linear scale
        yscale_flag = true

        # Flag to skip problematic snapshots
        #   + `true`: The snapshot will be skipped
        #   + `false`: The snapshot will be plotted
        skipper = true

        snapshot_number = snapshot_numbers[global_index]

        # Loop through each simulation for a given snapshot
        @inbounds for (simulation_index, simulation_table) in pairs(simulation_tables)

            snapshot_row = filter(:numbers => ==(snapshot_number), simulation_table)

            # Skip if this snapshot does not exist for the current simulation
            if isempty(snapshot_row)
                (
                    !warnings ||
                    @warn("snapshotPlot: The snapshot $(SNAP_BASENAME)_$(snapshot_number).hdf5 \
                    is missing in simulation $(simulation_paths[simulation_index])")
                )
                continue
            end

            ########################################################################################
            # Compute the metadata for the current snapshot and simulation.
            ########################################################################################

            # Get the snapshot file path
            snapshot_path = snapshot_row[1, :snapshot_paths]
            # Get the group catalog file path
            groupcat_path = snapshot_row[1, :groupcat_paths]

            # Skip the simulation if the snapshot is missing
            !ismissing(snapshot_path) || continue

            # Store the metadata of the current snapshot and simulation
            metadata = Dict(
                :sim_data => Simulation(
                    simulation_paths[simulation_index],
                    simulation_index,
                    slice,
                    isCosmological(snapshot_path),
                    simulation_table,
                ),
                :snap_data => Snapshot(
                    snapshot_path,
                    global_index,
                    slice_index,
                    snapshot_row[1, :physical_times],
                    snapshot_row[1, :lookback_times],
                    snapshot_row[1, :scale_factors],
                    snapshot_row[1, :redshifts],
                    readSnapHeader(snapshot_path),
                ),
                :gc_data => GroupCatalog(
                    groupcat_path,
                    readGroupCatHeader(groupcat_path; warnings),
                ),
            )

            ########################################################################################
            # Select the plot and data analysis functions.
            ########################################################################################

            # Get the plot function and its arguments for the current simulation
            plot_function = ring(plot_functions, simulation_index)
            pf_kwarg = ring(pf_kwargs, simulation_index)

            # Get the data analysis function and its arguments for the current simulation
            data_analysis = ring(da_functions, simulation_index)
            da_arg = ring(da_args, simulation_index)
            da_kwarg = ring(da_kwargs, simulation_index)

            # Get the plot styles for the current simulation
            color = ring(colors, simulation_index)
            marker = ring(markers, simulation_index)
            linestyle = ring(linestyles, simulation_index)

            ########################################################################################
            # Read and transform the data in the snapshot.
            ########################################################################################

            data_dict = merge(
                metadata,
                readSnapshot(snapshot_path, request; warnings),
                readGroupCatalog(groupcat_path, snapshot_path, request; warnings),
            )

            # Filter the data
            filterData!(data_dict; filter_function)

            if transform_box

                # Translate the data
                translateData!(data_dict, translation)

                # Rotate the data
                rotateData!(data_dict, rotation)

            end

            ########################################################################################
            # Compute the values to be plotted.
            ########################################################################################

            # Apply the analysis function
            da_output = data_analysis(data_dict, da_arg...; da_kwarg...)

            # Skip this snapshot if `data_analysis` returns `nothing`
            isnothing(da_output) ? continue : skipper = false

            # Data shape validation
            data_length = length(da_output)
            if plot_function isa typeof(heatmap!)
                (
                    data_length == 3 ||
                    error("snapshotPlot: For heatmaps `data_analysis` should return \
                    three data vectors, and currently is returning $(data_length)")
                )
            elseif plot_function isa typeof(hist!)
                (
                    data_length == 1 ||
                    error("snapshotPlot: For histograms `data_analysis` should return \
                    only one data vector, and currently is returning $(data_length)")
                )
            elseif plot_function isa Union{typeof(scatter!),typeof(scatterlines!),typeof(lines!)}
                (
                    data_length == 2 ||
                    error("snapshotPlot: For scatter/line plots `data_analysis` should return \
                    two data vectors, and currently is returning $(data_length)")
                )
            end

            # Unit conversion
            axis_data = VecOrMat{<:Number}[
                ustrip.(unit, quantity) for
                (unit, quantity) in zip([x_unit, y_unit, Unitful.NoUnits], da_output)
            ]

            # Sanitize the data
            if length(axis_data) == 1

                x_flag, _ = sanitizeData!(
                    axis_data[1];
                    func_domain=xaxis_scale_func,
                    range=x_trim,
                    keep_edges=x_edges,
                    min_left=0,
                    exp_factor=x_exp_factor,
                    warnings,
                )
                y_flag = true

                axis_data[1] = x_func(axis_data[1])

                # For histograms, if the scale it is not linear, compute the bin edges accordingly
                if x_flag

                    n_bins = (:bins ∈ keys(pf_kwarg)) ? pf_kwarg.bins : 10

                    bins = scaledBins(
                        axis_data[1],
                        n_bins;
                        scaling=xaxis_scale_func,
                        limits=xaxis_limits,
                    )

                    pf_kwarg = merge(pf_kwarg, (; bins))

                end

            elseif length(axis_data) == 2

                x_flag, y_flag, _, _ = sanitizeData!(
                    axis_data[1],
                    axis_data[2];
                    func_domain=(xaxis_scale_func, yaxis_scale_func),
                    range=(x_trim, y_trim),
                    keep_edges=(x_edges, y_edges),
                    min_left=0,
                    exp_factor=(x_exp_factor, y_exp_factor),
                    warnings,
                )

                # For scatter, line, and scatter line plots apply smoothing if required
                if iszero(smooth)

                    axis_data[1] = x_func(axis_data[1])
                    axis_data[2] = y_func(axis_data[2])

                else

                    axis_data[1], axis_data[2] = smoothWindow(
                        x_func(axis_data[1]),
                        y_func(axis_data[2]),
                        smooth;
                        scaling=xaxis_scale_func,
                    )

                end

            else

                sanitizeData!(
                    axis_data[1],
                    axis_data[2];
                    range=(x_trim, y_trim),
                    keep_edges=(x_edges, y_edges),
                    min_left=0,
                    exp_factor=(x_exp_factor, y_exp_factor),
                    warnings,
                )

                axis_data[1] = x_func(axis_data[1])
                axis_data[2] = y_func(axis_data[2])

                # For heatmaps, the scale functions of the axis are fixed as `identity`
                x_flag = false
                y_flag = false

            end

            # If, in the current snapshot and for any simulation, filtering the data targeting a
            # nonlinear scale would leave no data points, the scale will revert to `identity`
            x_flag || (xscale_flag = false)
            y_flag || (yscale_flag = false)

            if backup_results

                # Save data in a JLD2 file
                if isnothing(sim_labels)
                    sim_name = "simulation_$(lpad(string(simulation_index), 3, "0"))"
                else
                    sim_name = sim_labels[simulation_index]
                end

                jldopen(joinpath(output_path, base_filename * ".jld2"), "a+"; compress=true) do f
                    address = "$(base_filename)-$(SNAP_BASENAME)_$(snapshot_number)/$sim_name"
                    f[address] = axis_data
                end

            end

            if animation || save_figures
                # Draw the plot
                plot_function(axes, axis_data...; color, marker, linestyle, pf_kwarg...)
            end

        end

        # Skip snapshots for which no simulation could provide reasonable results
        # e.g. `data_analysis` returned `nothing` for every simulation
        if skipper
            next!(prog_bar)
            continue
        else
            plot_something = true
        end

        if animation || save_figures

            # Set the scale of the axes
            axes.xscale = (xscale_flag ? xaxis_scale_func : identity)
            axes.yscale = (yscale_flag ? yaxis_scale_func : identity)

            # Add a title
            longest_sim_table = argmax(nrow, simulation_tables)
            time_row = filter(:numbers => ==(snapshot_number), longest_sim_table)

            if isa(title, Symbol) && isempty(time_row)

                (
                    !warnings ||
                    @warn("snapshotPlot: I cound not find the time data for the snapshot \
                    number $(snapshot_number) in the longest running simulation with \
                    simulation table: \n$(longest_sim_table). \nDefaulting to using no title.")
                )

                axes.title = ""

            else

                if title == :physical_time

                    c_t = ustrip(u"Gyr", time_row[1, :physical_times])
                    time_stamp = round(c_t, digits=2)
                    axes.title = L"t = %$time_stamp \, \mathrm{Gyr}"

                elseif title == :lookback_time

                    p_t = ustrip(u"Gyr", time_row[1, :lookback_times])
                    time_stamp = round(p_t, digits=2)
                    axes.title = L"lt = %$time_stamp \, \mathrm{Gyr}"

                elseif title == :redshift

                    time_stamp = round(time_row[1, :redshifts], digits=2)
                    axes.title = L"z = \mathrm{%$time_stamp}"

                elseif title == :scale_factor

                    time_stamp = round(time_row[1, :scale_factors], digits=2)
                    axes.title = L"a = \mathrm{%$time_stamp}"

                else

                    # Use the user provided title
                    axes.title = title

                end

            end

            # Add a legend to the plot
            if !isnothing(sim_labels)

                (
                    length(sim_labels) == n_simulations ||
                    throw(ArgumentError("snapshotPlot: The arguments `simulation_paths` and \
                    `sim_labels` must have the same length, but I got length(`sim_labels`) = \
                    $(length(sim_labels)) != length(`simulation_paths`) = $(n_simulations)"))
                )

                legend_element = Vector{Makie.LegendElement}(undef, n_simulations)

                @inbounds for i in eachindex(legend_element)
                    color = ring(colors, i)
                    marker = ring(markers, i)
                    linestyle = ring(linestyles, i)
                    plot_function = ring(plot_functions, i)

                    if plot_function == lines!
                        legend_element[i] = LineElement(; color, linestyle)
                    else
                        legend_element[i] = MarkerElement(; color, marker)
                    end
                end

                nbanks = Int(ceil(length(sim_labels) * 0.3333))
                Makie.Legend(figure[2, 1][1, 1], legend_element, sim_labels; nbanks)

            end

            # Apply the post processing function
            pp_legend = post_processing(figure, pp_args...; pp_kwargs...)

            # Add the post processing legend
            if !isnothing(pp_legend)
                nbanks = Int(ceil(length(pp_legend) * 0.3333))
                Makie.Legend(
                    figure[2, 1][!isnothing(sim_labels) ? 2 : 1, 1],
                    pp_legend[1],
                    pp_legend[2];
                    nbanks,
                )
            end

        end

        if save_figures
            # Save the figure
            save(
                joinpath(
                    output_path,
                    "$(base_filename)-$(SNAP_BASENAME)_$(snapshot_number)$(output_format)",
                ),
                figure;
                pt_per_unit,
                px_per_unit,
            )
        end

        if animation
            # Add the figure as a frame to the animation stream
            recordframe!(vs)
        end

        # Clean the canvas for the next step in the loop
        cleanPlot!(figure)

        # Move the progress bar forward
        next!(prog_bar)

    end

    if warnings && !plot_something
        @warn("snapshotPlot: Nothing could be plotted because there was a problem \
        for every snapshot")
    end

    if animation
        # Save the animation
        save(joinpath(output_path, animation_filename), vs)
    end

    return nothing

end

"""
    timeSeriesPlot(
        simulation_paths::Vector{String},
        plot_functions::Vector{<:Function};
        <keyword arguments>
    )::Tuple{Makie.Axis,Figure}

Generate one figure per simulation.

Some of the features are:

  - It can produce scatter and line plots.
  - It transparently manages units; you only have to indicate the final unit of each axis.

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`.
  - `plot_functions::Vector{<:Function}`: Vector of plotting functions from [Makie](https://docs.makie.org/stable/). This sets the type of plot for each simulation.
    The supported functions are:

      + `scatter!`      -> Scatter plot.
      + `lines!`        -> Line plot.
      + `scatterlines!` -> Scatter plot with lines between the markers.
  - `pf_kwargs::Vector{<:NamedTuple}=[(;)]`: Vector of keyword arguments for the functions in `plot_functions`.

### timeSeriesPlot configuration

  - `output_path::String="./plots"`: Path to the output folder.
  - `filename::String="time_series"`: Filename for the figure, without the extension.
  - `output_format::String=".png"`: File format for the figure. All formats supported by [Makie](https://docs.makie.org/stable/) can be used, namely `.pdf`, `.svg` and `.png`.
  - `warnings::Bool=true`: If a warning will be raised when there is some problem, but the function can still run using sane defaults.
  - `show_progress::Bool=true`: If a progress bar will be shown.

### Data manipulation options

  - `slice::IndexType=(:)`: Slice of the simulation, i.e. which snapshots will be read. It can be an integer (a single snapshot), a vector of integers (several snapshots), an `UnitRange` (e.g. 5:13), an `StepRange` (e.g. 5:2:13) or (:) (all snapshots). It works over the longest possible list of snapshots among the simulations (grouped by the number in the file names). Out of bounds indices are ignored.
  - `da_functions::Vector{<:Function}=[getNothing]`: Vector of data analysis functions. See the required signature and examples in `./src/data_analysis.jl`.
  - `da_args::Vector{<:Tuple}=[()]`: Vector of positional arguments for the data analysis functions.
  - `da_kwargs::Vector{<:NamedTuple}=[(;)]`: Vector of keyword arguments for the data analysis functions.
  - `post_processing::Function=getNothing`: Post processing function. See the required signature and examples in `./src/post_processing.jl`.
  - `pp_args::Tuple=()`: Positional arguments for the post processing function.
  - `pp_kwargs::NamedTuple=(;)`: Keyword arguments for the post processing function.
  - `x_unit::Unitful.Units=Unitful.NoUnits`: Target unit for the x axis. The values will be converted accordingly. Use the default value of `Unitful.NoUnits` for dimensionless quantities.
  - `y_unit::Unitful.Units=Unitful.NoUnits`: Target unit for the y axis. The values will be converted accordingly. Use the default value of `Unitful.NoUnits` for dimensionless quantities.
  - `x_exp_factor::Int=0`: Numerical exponent to scale down the x axis, e.g. if `x_exp_factor` = 10 the values will be divided by ``10^{10}``. The default is no scaling.
  - `y_exp_factor::Int=0`: Numerical exponent to scale down the y axis, e.g. if `y_exp_factor` = 10 the values will be divided by ``10^{10}``. The default is no scaling.
  - `x_trim::NTuple{2,<:Real}=(-Inf, Inf)`: The data will be trim down so the x coordinates fit within `x_trim`.
  - `y_trim::NTuple{2,<:Real}=(-Inf, Inf)`: The data will be trim down so the y coordinates fit within `y_trim`.
  - `x_edges::Bool=false`: Set it to `true` if you want to keep the borders of `x_trim`.
  - `y_edges::Bool=false`: Set it to `true` if you want to keep the borders of `y_trim`.
  - `x_func::Function=identity`: Function to be applied to the values of the x axis. It must be a pure function with the signature `x_func(x_values::Vector{Float64})::Vector{Float64}`. The output must have the same length as the input. This function will be applied regardless of units and possible domain problems (use `x_trim` to solve incompatibilities), and that it will not be reflected in the automatic labeling.
  - `y_func::Function=identity`: Function to be applied to the values of the y axis. It must be a pure function with the signature `y_func(y_values::Vector{Float64})::Vector{Float64}`. The output must have the same length as the input. This function will be applied regardless of units and possible domain problems (use `y_trim` to solve incompatibilities), and that it will not be reflected in the automatic labeling.

### Axes options

  - `xaxis_label::AbstractString="auto_label"`: Label for the x axis. It can contain the string `auto_label`, which will be replaced by: `xaxis_var_name` / 10^`x_exp_factor` `x_unit`.
  - `yaxis_label::AbstractString="auto_label"`: Label for the y axis. It can contain the string `auto_label`, which will be replaced by: `yaxis_var_name` / 10^`y_exp_factor` `y_unit`.
  - `xaxis_var_name::AbstractString=""`: Name of the variable for the x axis.
  - `yaxis_var_name::AbstractString=""`: Name of the variable for the y axis.
  - `xaxis_scale_func::Function=identity`: Scaling function for the x axis. The options are the scaling functions accepted by [Makie](https://docs.makie.org/stable/): log10, log2, log, sqrt, Makie.logit, Makie.Symlog10, Makie.pseudolog10, and identity.
  - `yaxis_scale_func::Function=identity`: Scaling function for the y axis. The options are the scaling functions accepted by [Makie](https://docs.makie.org/stable/): log10, log2, log, sqrt, Makie.logit, Makie.Symlog10, Makie.pseudolog10, and identity.
  - `xaxis_limits::Tuple{<:Union{Real,Nothing},<:Union{Real,Nothing}}=(nothing, nothing)`: Set it to a value different than `nothing` if you want to fix the limits of the x axis.
  - `yaxis_limits::Tuple{<:Union{Real,Nothing},<:Union{Real,Nothing}}=(nothing, nothing)`: Set it to a value different than `nothing` if you want to fix the limits of the y axis.

### Plotting options

  - `save_figure::Bool=true`: If the plot will be saved as an image.
  - `backup_results::Bool=false`: If the values to be plotted will be backup in a [JLD2](https://github.com/JuliaIO/JLD2.jl) file.
  - `sim_labels::Union{Vector{String},Nothing}=nothing`: Labels for the plot legend, one per simulation. Set it to `nothing` if you don't want a legend.
  - `title::AbstractString=""`: Title for the figure. If left empty, no title will be printed.
  - `pt_per_unit::Float64=0.75`: Factor to scale up or down the size of the figures, keeping the proportions. It only works for `.pdf` and `.svg`.
  - `px_per_unit::Float64=1.0`: Factor to scale up or down the size of the figures, keeping the proportions. It only works for `.png`.
  - `size::NTuple{2,Int}=(1000, 750)`: Size of the figures in points. For PNGs, by default points = pixels (as given by `px_per_unit` = 1.0), and for PDFs and SVGs, points = 0.75 * pixels (as given by `pt_per_unit` = 0.75).
  - `aspect::Union{DataAspect,AxisAspect,Nothing}=nothing`: Aspect ratio of the figures. The options are:

      + `nothing`       -> Default, the aspect ratio will be chosen by [Makie](https://docs.makie.org/stable/).
      + `AxisAspect(n)` -> The aspect ratio will be given by the number `n` = width / height.
      + `DataAspect()`  -> The aspect ratio of the data will be used.
  - `series_colors::Union{Vector{<:ColorType},Nothing}=nothing`: Colors for the different simulations. If set to `nothing`, the colors will be assigned automatically.
  - `series_markers::Union{Vector{Symbol},Nothing}=nothing`: Markers for the different simulations. If set to `nothing`, the markers will be assigned automatically. This is only relevant for `scatter!` and `scatterlines!` plots.
  - `series_linestyles::Union{Vector{<:LineStyleType},Nothing}=nothing`: Line styles for the different simulations. If set to `nothing`, the line styles will be assigned automatically. This is only relevant for `lines!` and `scatterlines!` plots.

# Returns

  - The `Axis` and `Figure` objects.
"""
function timeSeriesPlot(
    simulation_paths::Vector{String},
    plot_functions::Vector{<:Function};
    pf_kwargs::Vector{<:NamedTuple}=[(;)],
    # `timeSeriesPlot` configuration
    output_path::String="./plots",
    filename::String="time_series",
    output_format::String=".png",
    warnings::Bool=true,
    show_progress::Bool=true,
    # Data manipulation options
    slice::IndexType=(:),
    da_functions::Vector{<:Function}=[getNothing],
    da_args::Vector{<:Tuple}=[()],
    da_kwargs::Vector{<:NamedTuple}=[(;)],
    post_processing::Function=getNothing,
    pp_args::Tuple=(),
    pp_kwargs::NamedTuple=(;),
    x_unit::Unitful.Units=Unitful.NoUnits,
    y_unit::Unitful.Units=Unitful.NoUnits,
    x_exp_factor::Int=0,
    y_exp_factor::Int=0,
    x_trim::NTuple{2,<:Real}=(-Inf, Inf),
    y_trim::NTuple{2,<:Real}=(-Inf, Inf),
    x_edges::Bool=false,
    y_edges::Bool=false,
    x_func::Function=identity,
    y_func::Function=identity,
    # Axes options
    xaxis_label::AbstractString="auto_label",
    yaxis_label::AbstractString="auto_label",
    xaxis_var_name::AbstractString="",
    yaxis_var_name::AbstractString="",
    xaxis_scale_func::Function=identity,
    yaxis_scale_func::Function=identity,
    xaxis_limits::Tuple{<:Union{Real,Nothing},<:Union{Real,Nothing}}=(nothing, nothing),
    yaxis_limits::Tuple{<:Union{Real,Nothing},<:Union{Real,Nothing}}=(nothing, nothing),
    # Plotting options
    save_figure::Bool=true,
    backup_results::Bool=false,
    sim_labels::Union{Vector{String},Nothing}=nothing,
    title::AbstractString="",
    pt_per_unit::Float64=0.75,
    px_per_unit::Float64=1.0,
    size::NTuple{2,Int}=(1000, 750),
    aspect::Union{DataAspect,AxisAspect,Nothing}=nothing,
    series_colors::Union{Vector{<:ColorType},Nothing}=nothing,
    series_markers::Union{Vector{Symbol},Nothing}=nothing,
    series_linestyles::Union{Vector{<:LineStyleType},Nothing}=nothing,
)::Tuple{Makie.Axis,Figure}

    # Create the output folder if it doesn't exist
    mkpath(output_path)

    # Compute the number of simulations
    n_simulations = length(simulation_paths)

    ################################################################################################
    # Set up the canvas for the figures.
    ################################################################################################

    # Apply the global theme defined in `./src/constants.jl`
    set_theme!(THEME)

    # Create the figure
    figure = Figure(; size)

    # Create a list of distinguishable colors for the different simulations
    if isnothing(series_colors)
        colors = distinguishable_colors(n_simulations, [RGB(1, 1, 1), RGB(0, 0, 0)], dropseed=true)
    else
        colors = series_colors
    end

    # Select a set of markers for the different simulations,
    # only relevant for `scatter!` and `scatterlines!` plots
    markers = isnothing(series_markers) ? MARKERS : series_markers

    # Select a set of line styles for the different simulations,
    # only relevant for `lines!` and `scatterlines!` plots
    linestyles = isnothing(series_linestyles) ? LINE_STYLES : series_linestyles

    # Create the labels
    xlabel = LaTeXString(
        replace(xaxis_label, "auto_label" => getLabel(xaxis_var_name, x_exp_factor, x_unit)),
    )
    ylabel = LaTeXString(
        replace(yaxis_label, "auto_label" => getLabel(yaxis_var_name, y_exp_factor, y_unit)),
    )

    # Create the axes
    axes = Makie.Axis(
        figure[1, 1];
        limits=(xaxis_limits, yaxis_limits),
        xlabel,
        ylabel,
        aspect,
        title,
    )

    ################################################################################################
    # Main loop.
    ################################################################################################

    # Initialize the progress bar
    prog_bar = Progress(
        n_simulations,
        dt=0.5,
        desc="Analyzing and plotting the data... ",
        color=:blue,
        barglyphs=BarGlyphs("|#  |"),
        enabled=show_progress,
    )

    # Flag to warn if nothing is plotted (because every snapshot was skipped)
    #   + `true`: At least one snapshot was plotted
    #   + `false`: Every snapshot was skipped
    plot_something = false

    # Flag to keep the x axis with a linear scale if there are no data points left after
    # trying to use a nonlinear scale
    #   + `true`: The x axis will use the scale given by `xaxis_scale_func`
    #   + `false`: The x axis will use a linear scale
    xscale_flag = true

    # Flag to keep the y axis with a linear scale if there are no data points left after
    # trying to use a nonlinear scale
    #   + `true`: The y axis will use the scale given by `yaxis_scale_func`
    #   + `false`: The y axis will use a linear scale
    yscale_flag = true

    # Loop through each simulation
    @inbounds for (simulation_index, simulation_path) in pairs(simulation_paths)

        ############################################################################################
        # Compute the metadata for the current simulation.
        ############################################################################################

        # Make a dataframe with the following columns:
        #   - 1. DataFrame index
        #   - 2. Number in the file name
        #   - 3. Scale factor
        #   - 4. Redshift
        #   - 5. Physical time
        #   - 6. Lookback time
        #   - 7. Snapshot path
        #   - 8. Group catalog path
        simulation_table = makeSimulationTable(simulation_path; warnings)

        # Get the path to the first snapshot
        snap_paths = getSnapshotPaths(simulation_path; warnings)
        _, idx = findmin(snap_paths[:numbers])
        first_snapshot = snap_paths[:paths][idx]

        # Store the metadata of the current simulation
        sim_data = Simulation(
            simulation_path,
            simulation_index,
            slice,
            isCosmological(first_snapshot),
            simulation_table,
        )

        ############################################################################################
        # Select the plot and data analysis functions.
        ############################################################################################

        # Get the plot function and its arguments for the current simulation
        plot_function = ring(plot_functions, simulation_index)
        pf_kwarg = ring(pf_kwargs, simulation_index)

        # Get the data analysis function and its arguments for the current simulation
        data_analysis = ring(da_functions, simulation_index)
        da_arg = ring(da_args, simulation_index)
        da_kwarg = ring(da_kwargs, simulation_index)

        # Get the plot styles for the current simulation
        color = ring(colors, simulation_index)
        marker = ring(markers, simulation_index)
        linestyle = ring(linestyles, simulation_index)

        ############################################################################################
        # Compute the values to be plotted.
        ############################################################################################

        # Apply the analysis function
        da_output = data_analysis(sim_data, da_arg...; da_kwarg...)

        # Skip this simulation if `data_analysis` returns `nothing`
        isnothing(da_output) ? continue : plot_something = true

        # Unit conversion
        axis_data = [ustrip.(x_unit, da_output[1]), ustrip.(y_unit, da_output[2])]

        # Sanitize the data
        x_flag, y_flag, _, _ = sanitizeData!(
            axis_data[1],
            axis_data[2];
            func_domain=(xaxis_scale_func, yaxis_scale_func),
            range=(x_trim, y_trim),
            keep_edges=(x_edges, y_edges),
            min_left=1,
            exp_factor=(x_exp_factor, y_exp_factor),
            warnings,
        )

        axis_data[1] = x_func(axis_data[1])
        axis_data[2] = y_func(axis_data[2])

        # If, for any simulation, filtering the data targeting a nonlinear scale
        # would leave no data points, the scale will revert to `identity`
        x_flag || (xscale_flag = false)
        y_flag || (yscale_flag = false)

        if backup_results

            # Save data in a JLD2 file
            if isnothing(sim_labels)
                sim_name = "simulation-$(lpad(string(simulation_index), 3, "0"))"
            else
                sim_name = sim_labels[simulation_index]
            end

            jldopen(joinpath(output_path, "$(filename).jld2"), "a+"; compress=true) do f
                address = "$(filename)/$sim_name"
                f[address] = axis_data
            end

        end

        if save_figure
            # Draw the plot
            plot_function(axes, axis_data...; color, marker, linestyle, pf_kwarg...)
        end

        # Move the progress bar forward
        next!(prog_bar)

    end

    if warnings && !plot_something
        @warn("timeSeriesPlot: Nothing could be plotted because there was a problem \
        for every snapshot")
    end

    if save_figure

        # Set the scale of the axes
        axes.xscale = (xscale_flag ? xaxis_scale_func : identity)
        axes.yscale = (yscale_flag ? yaxis_scale_func : identity)

        # Add a legend to the plot
        if !isnothing(sim_labels)

            (
                length(sim_labels) == n_simulations ||
                throw(ArgumentError("timeSeriesPlot: The arguments `simulation_paths` and \
                `sim_labels` must have the same length, but I got length(`sim_labels`) = \
                $(length(sim_labels)) != length(`simulation_paths`) = $(n_simulations)"))
            )

            legend_element = Vector{Makie.LegendElement}(undef, n_simulations)

            @inbounds for i in eachindex(legend_element)
                color = ring(colors, i)
                marker = ring(markers, i)
                linestyle = ring(linestyles, i)
                plot_function = ring(plot_functions, i)

                if plot_function == lines!
                    legend_element[i] = LineElement(; color, linestyle)
                else
                    legend_element[i] = MarkerElement(; color, marker)
                end
            end

            nbanks = Int(ceil(length(sim_labels) * 0.3333))
            Makie.Legend(figure[2, 1][1, 1], legend_element, sim_labels; nbanks)

        end

        # Apply the post processing function
        pp_legend = post_processing(figure, pp_args...; pp_kwargs...)

        # Add the post processing legend
        if !isnothing(pp_legend)
            nbanks = Int(ceil(length(pp_legend) * 0.3333))
            Makie.Legend(
                figure[2, 1][!isnothing(sim_labels) ? 2 : 1, 1],
                pp_legend[1],
                pp_legend[2];
                nbanks,
            )
        end

        # Save the figure
        save(joinpath(output_path, filename * output_format), figure; pt_per_unit, px_per_unit)
    end

    return axes, figure

end
