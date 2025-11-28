####################################################################################################
# Pipeline functions
####################################################################################################

"""
    plotSnapshot(
        simulation_paths::Vector{String},
        request::Dict{Symbol,Vector{String}},
        plot_functions::Vector{<:Function};
        <keyword arguments>
    )::Nothing

Generate one figure per snapshot, for one or more simulations.

Some of the features are:

  - It can produce scatter plots, line plots, histograms, bar plots, and heatmaps.
  - It can generate an animation of the results.
  - It transparently manages units; you only have to indicate the final unit of each axis.

!!! note

    The snapshots of different simulations are grouped by their filename, regardless of the "Time" parameter in the header.
    The data from the longest running simulation is used for the time stamp in the automatic title.

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`.
  - `request::Dict{Symbol,Vector{String}}`: Dictionary with the shape `cell/particle type` -> [`block`, `block`, ...], where the possible types are the keys of [`PARTICLE_INDEX`](@ref), and the possible quantities are the keys of [`QUANTITIES`](@ref). Which data blocks are needed depends on `da_functions`.
  - `plot_functions::Vector{<:Function}`: Vector of plotting functions from [Makie](https://docs.makie.org/stable/). This sets the type of plot for each simulation.
    The supported functions are:

      + `scatter!`      -> Scatter plot.
      + `lines!`        -> Line plot.
      + `scatterlines!` -> Scatter plot with lines between the markers.
      + `hist!`         -> Histogram.
      + `heatmap!`      -> Heatmap.
      + `arrows2d!`     -> 2D vector field.
      + `barplot!`      -> Bar plots.
      + `band!`         -> Band plots.
      + `errorbars!`    -> Error bars.
  - `pf_kwargs::Vector{<:NamedTuple}=[(;)]`: Vector of keyword arguments for the `plot_functions`.

### plotSnapshot configuration

  - `output_path::String="."`: Path to the output folder.
  - `base_filename::String="plot"`: Every file will be named `base_filename`_$(SNAP_BASENAME)_XXX`output_format` where XXX is the snapshot number.
  - `output_format::String=".png"`: File format for the figure. All formats supported by [Makie](https://docs.makie.org/stable/) can be used, namely `.pdf`, `.svg`, and `.png`.
  - `show_progress::Bool=true`: If a progress bar will be shown.

### Data manipulation options

  - `slice::IndexType=(:)`: Slice of the simulation, i.e. which snapshots will be plotted. It can be an integer (a single snapshot), a vector of integers (several snapshots), an `UnitRange` (e.g. 5:13), an `StepRange` (e.g. 5:2:13) or (:) (all snapshots). It works over the longest simulation. Out of bounds indices are ignored.
  - `transform_box::Bool=false`: If a translation and rotation (in that order) will be applied to the simulation box, affecting the positions and velocities of all the cells and particles. If active, it is applied BEFORE the `filter_function`.
  - `translation::TranslationType=:zero`: Type of translation (only relevant if `transform_box` = true). See the arguments after `data_dict` in [`translateData!`](@ref) for options. Several arguments can be passed as a Tuple.
  - `rotation::RotationType=:zero`: Type of rotation (only relevant if `transform_box` = true). See the arguments after `data_dict` in [`rotateData!`](@ref) for options. Several arguments can be passed as a Tuple.
  - `filter_function::Function=filterNothing`: Filter function. See the required signature and examples in `./src/analysis/filters.jl`. This is applied AFTER the translation and rotation.
  - `da_functions::Vector{<:Function}=[getNothing]`: Vector of data analysis functions. See the required signature and examples in `./src/analysis/data_analysis.jl`.
  - `da_args::Vector{<:Tuple}=[()]`: Vector of positional arguments for the data analysis functions.
  - `da_kwargs::Vector{<:NamedTuple}=[(;)]`: Vector of keyword arguments for the data analysis functions.
  - `post_processing::Function=getNothing`: Post processing function. See the required signature and examples in `./src/plotting/post_processing.jl`.
  - `pp_args::Tuple=()`: Positional arguments for the post processing function.
  - `pp_kwargs::NamedTuple=(;)`: Keyword arguments for the post processing function.
  - `smooth::Int=0`: The first two data vectors of each `da_functions` will be smoothed out using `smooth` bins. Set it to 0 if you want no smoothing.

### Axes options

  - `x_unit::Unitful.Units=Unitful.NoUnits`: Target unit for the x axis data. The values will be converted accordingly. Use the default value of `Unitful.NoUnits` for dimensionless quantities.
  - `y_unit::Unitful.Units=Unitful.NoUnits`: Target unit for the y axis data. The values will be converted accordingly. Use the default value of `Unitful.NoUnits` for dimensionless quantities.
  - `x_exp_factor::Int=0`: Numerical exponent to scale down the x axis data, e.g. if `x_exp_factor` = 10 the values will be divided by ``10^{10}``. The default is no scaling.
  - `y_exp_factor::Int=0`: Numerical exponent to scale down the y axis data, e.g. if `y_exp_factor` = 10 the values will be divided by ``10^{10}``. The default is no scaling.
  - `x_trim::NTuple{2,<:Real}=(-Inf, Inf)`: The data will be trim down so the x coordinates fit within `x_trim`, in the units given by `x_unit`.
  - `y_trim::NTuple{2,<:Real}=(-Inf, Inf)`: The data will be trim down so the y coordinates fit within `y_trim`, in the units given by `y_unit`. This option does not affect histograms.
  - `x_edges::Bool=false`: Set it to `true` if you want to keep the borders of `x_trim`.
  - `y_edges::Bool=false`: Set it to `true` if you want to keep the borders of `y_trim`.
  - `x_scale_func::Function=identity`: Scaling function for the x axis. The options are the scaling functions accepted by [Makie](https://docs.makie.org/stable/): log10, log2, log, sqrt, Makie.logit, Makie.Symlog10, Makie.pseudolog10, and identity. The data will be trimmed down to fit within the domain of `x_scale_func`.
  - `y_scale_func::Function=identity`: Scaling function for the y axis. The options are the scaling functions accepted by [Makie](https://docs.makie.org/stable/): log10, log2, log, sqrt, Makie.logit, Makie.Symlog10, Makie.pseudolog10, and identity. The data will be trimmed down to fit within the domain of `y_scale_func`.
  - `xaxis_label::AbstractString="auto_label"`: Label for the x axis. It can contain the string `auto_label`, which will be replaced by: `xaxis_var_name` [10^`x_exp_factor` `x_unit`]. If a LaTeXString with `auto_label` inside is used, it is recommended that each section around `auto_label` is delimited with a `\$ \$` pair.
  - `yaxis_label::AbstractString="auto_label"`: Label for the y axis. It can contain the string `auto_label`, which will be replaced by: `yaxis_var_name` [10^`y_exp_factor` `y_unit`]. If a LaTeXString with `auto_label` inside is used, it is recommended that each section around `auto_label` is delimited with a `\$ \$` pair.
  - `xaxis_var_name::AbstractString="x"`: Name of the variable for the x axis.
  - `yaxis_var_name::AbstractString="y"`: Name of the variable for the y axis.

### Plotting options

  - `save_figures::Bool=true`: If every figure will be saved as a file.
  - `backup_results::Bool=false`: If the values to be plotted will be saved in a [JLD2](https://github.com/JuliaIO/JLD2.jl) file.
  - `theme::Attributes=Theme()`: Plot theme that will take precedence over [`DEFAULT_THEME`](@ref).
  - `sim_labels::Union{Vector{<:Union{AbstractString,Nothing}},Nothing}=nothing`: Labels for the plot legend, one per simulation. Set it to `nothing` if you don't want a legend.
  - `title::Union{Symbol,<:AbstractString}=""`: Title for the figure. If left empty, no title is printed. It can also be set to one of the following options:

      + `:physical_time` -> Physical time since the Big Bang.
      + `:lookback_time` -> Physical time left to reach the last snapshot.
      + `:scale_factor`  -> Scale factor (only relevant for cosmological simulations).
      + `:redshift`      -> Redshift (only relevant for cosmological simulations).
  - `colorbar::Bool=false`: If a vertical colorbar will be added.

## Animation options

  - `animation::Bool=false`: If an animation will be created.
  - `animation_filename::String="{base_filename}.mkv"`: Filename for the animation, including its extension. All formats supported by [Makie](https://docs.makie.org/stable/) can be used, namely `.mkv`, `.mp4`, `.webm`, and `.gif`.
  - `framerate::Int=20`: Frame rate of the animation.
"""
function plotSnapshot(
    simulation_paths::Vector{String},
    request::Dict{Symbol,Vector{String}},
    plot_functions::Vector{<:Function};
    pf_kwargs::Vector{<:NamedTuple}=[(;)],
    # `plotSnapshot` configuration
    output_path::String=".",
    base_filename::String="plot",
    output_format::String=".png",
    show_progress::Bool=true,
    # Data manipulation options
    slice::IndexType=(:),
    transform_box::Bool=false,
    translation::TranslationType=:zero,
    rotation::RotationType=:zero,
    filter_function::Function=filterNothing,
    da_functions::Vector{<:Function}=[getNothing],
    da_args::Vector{<:Tuple}=[()],
    da_kwargs::Vector{<:NamedTuple}=[(;)],
    post_processing::Function=getNothing,
    pp_args::Tuple=(),
    pp_kwargs::NamedTuple=(;),
    smooth::Int=0,
    # Axes options
    x_unit::Unitful.Units=Unitful.NoUnits,
    y_unit::Unitful.Units=Unitful.NoUnits,
    x_exp_factor::Int=0,
    y_exp_factor::Int=0,
    x_trim::NTuple{2,<:Real}=(-Inf, Inf),
    y_trim::NTuple{2,<:Real}=(-Inf, Inf),
    x_edges::Bool=false,
    y_edges::Bool=false,
    x_scale_func::Function=identity,
    y_scale_func::Function=identity,
    xaxis_label::AbstractString="auto_label",
    yaxis_label::AbstractString="auto_label",
    xaxis_var_name::AbstractString="x",
    yaxis_var_name::AbstractString="y",
    # Plotting options
    save_figures::Bool=true,
    backup_results::Bool=false,
    theme::Attributes=Theme(),
    sim_labels::Union{Vector{<:Union{AbstractString,Nothing}},Nothing}=nothing,
    title::Union{Symbol,<:AbstractString}="",
    colorbar::Bool=false,
    # Animation options
    animation::Bool=false,
    animation_filename::String="$(base_filename).mkv",
    framerate::Int=20,
)::Nothing

    # Create the output folder if it doesn't exist
    mkpath(output_path)

    # Compute the number of simulations
    n_simulations = length(simulation_paths)

    # Make a dataframe for each simulation with the following columns:
    #  - DataFrame index         -> :row_id
    #  - Number in the file name -> :numbers
    #  - Scale factor            -> :scale_factors
    #  - Redshift                -> :redshifts
    #  - Physical time           -> :physical_times
    #  - Lookback time           -> :lookback_times
    #  - Snapshot path           -> :snapshot_paths
    #  - Group catalog path      -> :groupcat_paths
    simulation_tables = [makeSimulationTable(source) for source in simulation_paths]

    # Find the longest one
    longest_sim_table = argmax(nrow, simulation_tables)

    # Compute the different ways to index the snapshots
    snapshot_numbers = longest_sim_table[!, :numbers]
    slice_indices    = safeSelect(longest_sim_table[!, :row_id], slice)

    # Compute the number of figures
    n_frames = length(slice_indices)

    # Check that after slicing there is at least one snapshot left
    (
        !iszero(n_frames) ||
        throw(ArgumentError("plotSnapshot: There are no snapshots left in the longest simulation \
        after slicing with `slice` = $slice"))
    )

    ################################################################################################
    # Set up the canvas for the figures
    ################################################################################################

    # Set up the plot theme
    current_theme = merge(theme, DEFAULT_THEME, theme_latexfonts())

    # Apply the plot theme
    set_theme!(current_theme)

    # Create the figure
    figure = Figure()

    # Create the labels
    xlabel = LaTeXString(
        replace(xaxis_label, "auto_label" => getLabel(xaxis_var_name, x_exp_factor, x_unit)),
    )
    ylabel = LaTeXString(
        replace(yaxis_label, "auto_label" => getLabel(yaxis_var_name, y_exp_factor, y_unit)),
    )

    # Create the axes
    axes = Makie.Axis(figure[1, 1]; xlabel, ylabel)

    ################################################################################################
    # Set up the animation
    ################################################################################################

    if animation

        (
            n_frames < framerate &&
            logging[] &&
            @warn("plotSnapshot: With `framerate` = $framerate and `slice` = $slice, \
            the animation is less than one second long")
        )

        # Initialize the animation stream
        vs = VideoStream(figure; framerate)

    end

    ################################################################################################
    # Main loop
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

    # Flag to warn if nothing is plotted because every snapshot was skipped
    #   + `true`: At least one snapshot was plotted
    #   + `false`: Every snapshot was skipped
    plot_something = false

    # Loop through each snapshots
    for (slice_index, global_index) in pairs(slice_indices)

        ############################################################################################
        # Flags
        ############################################################################################

        # Flag to keep the x axis with a linear scale if there are no data points left after
        # trying to use a nonlinear scale
        #   + `true`: The x axis will use the scale given by `x_scale_func`
        #   + `false`: The x axis will use a linear scale
        xscale_flag = true

        # Flag to keep the y axis with a linear scale if there are no data points left after
        # trying to use a nonlinear scale
        #   + `true`: The y axis will use the scale given by `y_scale_func`
        #   + `false`: The y axis will use a linear scale
        yscale_flag = true

        # Flag to skip problematic snapshots
        #   + `true`: The snapshot will be skipped
        #   + `false`: The snapshot will be plotted
        skipper = true

        snapshot_number = snapshot_numbers[global_index]

        # Loop through each simulation, for the current snapshot
        for (simulation_index, simulation_table) in pairs(simulation_tables)

            # Select the row corresponding to the current snapshot from the simulation table
            snapshot_row = filter(:numbers => ==(snapshot_number), simulation_table)

            # Read the simulation path
            simulation_path = simulation_paths[simulation_index]

            # Skip if this snapshot does not exists for the current simulation
            if isempty(snapshot_row)
                (
                    logging[] &&
                    @warn("plotSnapshot: The snapshot $(SNAP_BASENAME)_$(snapshot_number).hdf5 \
                    is missing in simulation $(simulation_path)")
                )
                continue
            end

            ########################################################################################
            # Compute the metadata for the current snapshot and simulation
            ########################################################################################

            # Read the snapshot file path
            snapshot_path = snapshot_row[1, :snapshot_paths]

            # Read the group catalog file path
            groupcat_path = snapshot_row[1, :groupcat_paths]

            # Skip the simulation if the snapshot is missing
            if ismissing(snapshot_path)
                (
                    logging[] &&
                    @warn("plotSnapshot: The snapshot $(SNAP_BASENAME)_$(snapshot_number).hdf5 \
                    is missing in simulation $(simulation_path)")
                )
                continue
            end

            # Store the metadata of the current snapshot and simulation
            metadata = Dict(
                :sim_data => Simulation(
                    simulation_path,
                    simulation_index,
                    slice,
                    isSnapCosmological(snapshot_path),
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
                :gc_data => GroupCatalog(groupcat_path, readGroupCatHeader(groupcat_path)),
            )

            ########################################################################################
            # Read and transform the data in the snapshot
            ########################################################################################

            data_dict = merge(
                metadata,
                readSnapshot(snapshot_path, request),
                readGroupCatalog(groupcat_path, snapshot_path, request),
            )

            if transform_box

                # Translate the data
                translateData!(data_dict, translation...)

                # Rotate the data
                rotateData!(data_dict, rotation...)

            end

            # Filter the data
            filterData!(data_dict; filter_function)

            ########################################################################################
            # Select the plot and data analysis functions
            ########################################################################################

            # Get the plot function and its arguments for the current simulation
            plot_function = ring(plot_functions, simulation_index)
            pf_kwarg      = ring(pf_kwargs, simulation_index)

            # Get the data analysis function and its arguments for the current simulation
            data_analysis = ring(da_functions, simulation_index)
            da_arg        = ring(da_args, simulation_index)
            da_kwarg      = ring(da_kwargs, simulation_index)

            ########################################################################################
            # Compute the values to be plotted
            ########################################################################################

            # Apply the analysis function
            da_output = data_analysis(data_dict, da_arg...; da_kwarg...)

            # Skip this snapshot if `data_analysis` returns `nothing`
            if isnothing(da_output)

                (
                    logging[] &&
                    @warn("plotSnapshot: The data analysis function $(data_analysis) returned \
                    `nothing` for simulation $(simulation_path) and snapshot \
                    $(SNAP_BASENAME)_$(snapshot_number).hdf5")
                )

                continue

            else

                skipper = false

            end

            ########################################################################################
            # Data sanitation
            ########################################################################################

            axis_data, x_flag, y_flag = validatePlotData(
                plot_function,
                da_output...;
                x_unit,
                y_unit,
                x_exp_factor,
                y_exp_factor,
                x_trim,
                y_trim,
                x_edges,
                y_edges,
                x_scale_func,
                y_scale_func,
            )

            if plot_function isa typeof(hist!) && x_flag
                # For histograms, if the scale is not linear, recompute the bin edges accordingly
                n_bins   = haskey(pf_kwarg, :bins) ? pf_kwarg.bins : 10
                bins     = scaledBins(axis_data[1], n_bins; scaling=x_scale_func)
                pf_kwarg = merge(pf_kwarg, (; bins))
            end

            # If, in the current snapshot and for any simulation, filtering the data targeting a
            # nonlinear scale would leave no data points, the scale will revert to `identity`
            x_flag || (xscale_flag = false)
            y_flag || (yscale_flag = false)

            ########################################################################################
            # Apply a smoothing window
            ########################################################################################

            if length(axis_data) > 1 && !iszero(smooth)
                axis_data[1], axis_data[2] = smoothWindow(
                    axis_data[1],
                    axis_data[2],
                    smooth;
                    scaling=x_scale_func,
                )
            end

            ########################################################################################
            # Draw the plot
            ########################################################################################

            if animation || save_figures

                if colorbar && plot_function isa typeof(heatmap!)

                    min_c, max_c = extrema(filter(!isnan, axis_data[3]))

                    min_color = round(min_c)
                    max_color = round(max_c)

                    pf_kwarg = merge((; colorrange=(min_color, max_color)), pf_kwarg)

                end

                plot_function(axes, axis_data...; pf_kwarg...)

                if colorbar && plot_function isa typeof(heatmap!)

                    if haskey(current_theme, Colorbar) && haskey(current_theme.Colorbar, :colorrange)
                        colorbar_cr = current_theme.Colorbar.colorrange
                    else
                        colorbar_cr = pf_kwarg.colorrange
                    end

                    Colorbar(figure[1, 2]; colorrange=colorbar_cr, ticks=WilkinsonTicks(5),)

                    rowsize!(figure.layout, 1, Auto(1.0))

                end

            end

            ########################################################################################
            # Save the results in a JLD2 file
            ########################################################################################

            if backup_results

                if isnothing(sim_labels)
                    sim_name = basename(simulation_path)
                else
                    sim_name = sim_labels[simulation_index]
                end

                jldopen(joinpath(output_path, base_filename * ".jld2"), "a+") do f
                    address = "$(SNAP_BASENAME)_$(snapshot_number)/$(sim_name)"
                    f[address] = axis_data
                end

            end

        end

        # Skip snapshots for which no simulation could provide reasonable results
        # e.g. `data_analysis` returned `nothing` for every simulation
        if skipper

            (
                logging[] &&
                @warn("plotSnapshot: The data analysis function $(data_analysis) returned \
                `nothing` for snapshot $(SNAP_BASENAME)_$(snapshot_number).hdf5 in every \
                simulation")
            )

            next!(prog_bar)

            continue

        else

            plot_something = true

        end

        if animation || save_figures

            # Set the scale of the axes
            axes.xscale = (xscale_flag ? x_scale_func : identity)
            axes.yscale = (yscale_flag ? y_scale_func : identity)

            ########################################################################################
            # Add a title
            ########################################################################################

            time_row = filter(:numbers => ==(snapshot_number), longest_sim_table)

            if isa(title, Symbol) && isempty(time_row)

                (
                    logging[] &&
                    @warn("plotSnapshot: I cound not find the time data for the snapshot \
                    $(SNAP_BASENAME)_$(snapshot_number).hdf5 in the longest running simulation. \
                    I will print no title.")
                )

                axes.title = ""

            else

                if title == :physical_time

                    c_t = ustrip(u"Gyr", time_row[1, :physical_times])

                    if c_t < 1.0
                        time_stamp = round(c_t; digits=2)
                    else
                        time_stamp = round(c_t; sigdigits=3)
                    end

                    axes.title = L"t = %$(rpad(time_stamp, 4, '0')) \, \text{Gyr}"

                elseif title == :lookback_time

                    p_t = ustrip(u"Gyr", time_row[1, :lookback_times])

                    if p_t < 1.0
                        time_stamp = round(p_t; digits=2)
                    else
                        time_stamp = round(p_t; sigdigits=3)
                    end

                    axes.title = L"lt = %$rpad(time_stamp, 4, '0')) \, \text{Gyr}"

                elseif title == :redshift

                    z_t = time_row[1, :redshifts]

                    if p_t < 1.0
                        time_stamp = round(z_t; digits=2)
                    else
                        time_stamp = round(z_t; sigdigits=3)
                    end

                    axes.title = L"z = \mathrm{%$rpad(time_stamp, 4, '0'))}"

                elseif title == :scale_factor

                    a_t = time_row[1, :scale_factors]

                    if p_t < 1.0
                        time_stamp = round(a_t; digits=2)
                    else
                        time_stamp = round(a_t; sigdigits=3)
                    end

                    axes.title = L"a = \mathrm{%$rpad(time_stamp, 4, '0'))}"

                else

                    # Use the provided title
                    axes.title = title

                end

            end

            ########################################################################################
            # Apply the post processing function
            ########################################################################################

            pp_legend = post_processing(figure, pp_args...; pp_kwargs...)

            ########################################################################################
            # Add the legends
            ########################################################################################

            legend_elements = Vector{Makie.LegendElement}(undef, 0)
            legend_labels = Vector{Union{AbstractString,Nothing}}(undef, 0)

            if !isnothing(sim_labels)
                # Add the main legend
                (
                    length(sim_labels) == n_simulations ||
                    throw(ArgumentError("plotSnapshot: The arguments `simulation_paths` and \
                    `sim_labels` must have the same length, but I got length(`sim_labels`) = \
                    $(length(sim_labels)) != length(`simulation_paths`) = $(n_simulations)"))
                )

                # Load the current palette
                colors     = current_theme[:palette][:color][]
                markers    = current_theme[:palette][:marker][]
                linestyles = current_theme[:palette][:linestyle][]

                for i in 1:n_simulations

                    color         = ring(colors, i)
                    marker        = ring(markers, i)
                    linestyle     = ring(linestyles, i)
                    plot_function = ring(plot_functions, i)

                    if plot_function == lines!
                        push!(legend_elements, LineElement(; color, linestyle))
                    else
                        push!(legend_elements, MarkerElement(; color, marker))
                    end

                end

                append!(legend_labels, sim_labels)

            end

            if !isnothing(pp_legend)
                # Add the post processing legend
                append!(legend_elements, pp_legend[1])
                append!(legend_labels, pp_legend[2])
            end

            if !any(isempty, [legend_elements, legend_labels])
                Makie.Legend(figure[1, 1], legend_elements, legend_labels)
            end

        end

        ############################################################################################
        # Save the figure
        ############################################################################################

        if save_figures

            output_filename = "$(base_filename)_$(SNAP_BASENAME)_$(snapshot_number)$(output_format)"

            save(joinpath(output_path, output_filename), figure)

        end

        # Add the figure as a frame to the animation stream
        animation && recordframe!(vs)

        # Clean the canvas for the next snapshot
        cleanPlot!(figure)

        # Move the progress bar forward
        next!(prog_bar)

    end

    (
        logging[] && !plot_something &&
        @warn("plotSnapshot: Nothing could be plotted because there was a problem for every \
        snapshot")
    )

    if animation
        # Save the animation
        save(joinpath(output_path, animation_filename), vs)
    end

    return nothing

end

"""
    plotTimeSeries(
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
  - `pf_kwargs::Vector{<:NamedTuple}=[(;)]`: Vector of keyword arguments for the functions `plot_functions`.

### plotTimeSeries configuration

  - `output_path::String="."`: Path to the output folder.
  - `filename::String="time_series"`: Filename for the figure, without the extension.
  - `output_format::String=".png"`: File format for the figure. All formats supported by [Makie](https://docs.makie.org/stable/) can be used, namely `.pdf`, `.svg` and `.png`.

### Data manipulation options

  - `slice::IndexType=(:)`: Slice of the simulation, i.e. which snapshots will be plotted. It can be an integer (a single snapshot), a vector of integers (several snapshots), an `UnitRange` (e.g. 5:13), an `StepRange` (e.g. 5:2:13) or (:) (all snapshots). It works over the longest possible list of snapshots among the simulations (grouped by the number in the file names). Out of bounds indices are ignored.
  - `da_functions::Vector{<:Function}=[getNothing]`: Vector of data analysis functions. See the required signature and examples in `./src/analysis/data_analysis.jl`.
  - `da_args::Vector{<:Tuple}=[()]`: Vector of positional arguments for the data analysis functions.
  - `da_kwargs::Vector{<:NamedTuple}=[(;)]`: Vector of keyword arguments for the data analysis functions.
  - `post_processing::Function=getNothing`: Post processing function. See the required signature and examples in `./src/plotting/post_processing.jl`.
  - `pp_args::Tuple=()`: Positional arguments for the post processing function.
  - `pp_kwargs::NamedTuple=(;)`: Keyword arguments for the post processing function.

### Axes options

  - `x_unit::Unitful.Units=Unitful.NoUnits`: Target unit for the x axis data. The values will be converted accordingly. Use the default value of `Unitful.NoUnits` for dimensionless quantities.
  - `y_unit::Unitful.Units=Unitful.NoUnits`: Target unit for the y axis data. The values will be converted accordingly. Use the default value of `Unitful.NoUnits` for dimensionless quantities.
  - `x_exp_factor::Int=0`: Numerical exponent to scale down the x axis data, e.g. if `x_exp_factor` = 10 the values will be divided by ``10^{10}``. The default is no scaling.
  - `y_exp_factor::Int=0`: Numerical exponent to scale down the y axis data, e.g. if `y_exp_factor` = 10 the values will be divided by ``10^{10}``. The default is no scaling.
  - `x_trim::NTuple{2,<:Real}=(-Inf, Inf)`: The data will be trim down so the x coordinates fit within `x_trim`, in the units given by `x_unit`.
  - `y_trim::NTuple{2,<:Real}=(-Inf, Inf)`: The data will be trim down so the y coordinates fit within `y_trim`, in the units given by `y_unit`.
  - `x_edges::Bool=false`: Set it to `true` if you want to keep the borders of `x_trim`.
  - `y_edges::Bool=false`: Set it to `true` if you want to keep the borders of `y_trim`.
  - `x_scale_func::Function=identity`: Scaling function for the x axis. The options are the scaling functions accepted by [Makie](https://docs.makie.org/stable/): log10, log2, log, sqrt, Makie.logit, Makie.Symlog10, Makie.pseudolog10, and identity. The data will be trimmed down to fit within the domain of `x_scale_func`.
  - `y_scale_func::Function=identity`: Scaling function for the y axis. The options are the scaling functions accepted by [Makie](https://docs.makie.org/stable/): log10, log2, log, sqrt, Makie.logit, Makie.Symlog10, Makie.pseudolog10, and identity. The data will be trimmed down to fit within the domain of `y_scale_func`.
  - `xaxis_label::AbstractString="auto_label"`: Label for the x axis. It can contain the string `auto_label`, which will be replaced by: `xaxis_var_name` [10^`x_exp_factor` `x_unit`]. If a LaTeXString with `auto_label` inside is used, it is recommended that each section around `auto_label` is delimited with a `\$ \$` pair.
  - `yaxis_label::AbstractString="auto_label"`: Label for the y axis. It can contain the string `auto_label`, which will be replaced by: `yaxis_var_name` [10^`y_exp_factor` `y_unit`]. If a LaTeXString with `auto_label` inside is used, it is recommended that each section around `auto_label` is delimited with a `\$ \$` pair.
  - `xaxis_var_name::AbstractString="x"`: Name of the variable for the x axis.
  - `yaxis_var_name::AbstractString="y"`: Name of the variable for the y axis.

### Plotting options

  - `save_figures::Bool=true`: If the plot will be saved as a file.
  - `backup_results::Bool=false`: If the values to be plotted will be saved in a [JLD2](https://github.com/JuliaIO/JLD2.jl) file.
  - `theme::Attributes=Theme()`: Plot theme that will take precedence over [`DEFAULT_THEME`](@ref).
  - `sim_labels::Union{Vector{<:Union{AbstractString,Nothing}},Nothing}=nothing`: Labels for the plot legend, one per simulation. Set it to `nothing` if you don't want a legend.
  - `title::AbstractString=""`: Title for the figure. If left empty, no title will be printed.

# Returns

  - The `Axis` and `Figure` objects.
"""
function plotTimeSeries(
    simulation_paths::Vector{String},
    plot_functions::Vector{<:Function};
    pf_kwargs::Vector{<:NamedTuple}=[(;)],
    # `plotTimeSeries` configuration
    output_path::String=".",
    filename::String="time_series",
    output_format::String=".png",
    # Data manipulation options
    slice::IndexType=(:),
    da_functions::Vector{<:Function}=[getNothing],
    da_args::Vector{<:Tuple}=[()],
    da_kwargs::Vector{<:NamedTuple}=[(;)],
    post_processing::Function=getNothing,
    pp_args::Tuple=(),
    pp_kwargs::NamedTuple=(;),
    # Axes options
    x_unit::Unitful.Units=Unitful.NoUnits,
    y_unit::Unitful.Units=Unitful.NoUnits,
    x_exp_factor::Int=0,
    y_exp_factor::Int=0,
    x_trim::NTuple{2,<:Real}=(-Inf, Inf),
    y_trim::NTuple{2,<:Real}=(-Inf, Inf),
    x_edges::Bool=false,
    y_edges::Bool=false,
    x_scale_func::Function=identity,
    y_scale_func::Function=identity,
    xaxis_label::AbstractString="auto_label",
    yaxis_label::AbstractString="auto_label",
    xaxis_var_name::AbstractString="x",
    yaxis_var_name::AbstractString="y",
    # Plotting options
    save_figures::Bool=true,
    backup_results::Bool=false,
    theme::Attributes=Theme(),
    sim_labels::Union{Vector{<:Union{AbstractString,Nothing}},Nothing}=nothing,
    title::AbstractString="",
)::Tuple{Makie.Axis,Figure}

    # Create the output folder if it doesn't exist
    mkpath(output_path)

    # Compute the number of simulations
    n_simulations = length(simulation_paths)

    ################################################################################################
    # Set up the canvas for the figures
    ################################################################################################

    # Set up the plot theme
    current_theme = merge(theme, DEFAULT_THEME, theme_latexfonts())

    # Apply the plot theme
    set_theme!(current_theme)

    # Create the figure
    figure = Figure()

    # Create the labels
    xlabel = LaTeXString(
        replace(xaxis_label, "auto_label" => getLabel(xaxis_var_name, x_exp_factor, x_unit)),
    )
    ylabel = LaTeXString(
        replace(yaxis_label, "auto_label" => getLabel(yaxis_var_name, y_exp_factor, y_unit)),
    )

    # Create the axes
    axes = Makie.Axis(figure[1, 1]; xlabel, ylabel, title)

    ################################################################################################
    # Main loop
    ################################################################################################

    # Flag to warn if nothing is plotted (because every snapshot was skipped)
    #   + `true`: At least one snapshot was plotted
    #   + `false`: Every snapshot was skipped
    plot_something = false

    # Flag to keep the x axis with a linear scale if there are no data points left after
    # trying to use a nonlinear scale
    #   + `true`: The x axis will use the scale given by `x_scale_func`
    #   + `false`: The x axis will use a linear scale
    xscale_flag = true

    # Flag to keep the y axis with a linear scale if there are no data points left after
    # trying to use a nonlinear scale
    #   + `true`: The y axis will use the scale given by `y_scale_func`
    #   + `false`: The y axis will use a linear scale
    yscale_flag = true

    # Loop through each simulation
    for (simulation_index, simulation_path) in pairs(simulation_paths)

        ############################################################################################
        # Compute the metadata for the current simulation
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

        # Store the metadata of the current simulation
        sim_data = Simulation(
            simulation_path,
            simulation_index,
            slice,
            isSimCosmological(simulation_path),
            simulation_table,
        )

        ############################################################################################
        # Select the plot and data analysis functions
        ############################################################################################

        # Get the plot function and its arguments for the current simulation
        plot_function = ring(plot_functions, simulation_index)
        pf_kwarg      = ring(pf_kwargs, simulation_index)

        # Get the data analysis function and its arguments for the current simulation
        data_analysis = ring(da_functions, simulation_index)
        da_arg        = ring(da_args, simulation_index)
        da_kwarg      = ring(da_kwargs, simulation_index)

        ############################################################################################
        # Compute the values to be plotted
        ############################################################################################

        # Apply the analysis function
        da_output = data_analysis(sim_data, da_arg...; da_kwarg...)

        # Skip this simulation if `data_analysis` returns `nothing`
        isnothing(da_output) ? continue : plot_something = true

        if isnothing(da_output)
            (
                logging[] &&
                @warn("plotTimeSeries: The data analysis function $(data_analysis) returned \
                `nothing` for simulation $(simulation_path)")
            )
            continue
        end

        ############################################################################################
        # Data sanitation
        ############################################################################################

        # Unit conversion
        axis_data = [ustrip.(x_unit, da_output[1]), ustrip.(y_unit, da_output[2])]

        x_flag, y_flag, _, _ = sanitizeData!(
            axis_data[1],
            axis_data[2];
            func_domain=(x_scale_func, y_scale_func),
            range=(x_trim, y_trim),
            keep_edges=(x_edges, y_edges),
            min_left=1,
            exp_factor=(x_exp_factor, y_exp_factor),
        )

        # If, for any simulation, filtering the data targeting a nonlinear scale
        # would leave no data points, the scale will revert to `identity`
        x_flag || (xscale_flag = false)
        y_flag || (yscale_flag = false)

        ############################################################################################
        # Draw the plot
        ############################################################################################

        if save_figures
            plot_function(axes, axis_data...; pf_kwarg...)
        end

        ############################################################################################
        # Save the results in a JLD2 file
        ############################################################################################

        if backup_results

            if isnothing(sim_labels)
                sim_name = basename(simulation_path)
            else
                sim_name = sim_labels[simulation_index]
            end

            jldopen(joinpath(output_path, "$(filename).jld2"), "a+"; compress=true) do f
                address = "$(filename)/$sim_name"
                f[address] = axis_data
            end

        end

    end

    (
        logging[] && !plot_something &&
        @warn("plotTimeSeries: Nothing could be plotted because there was a problem for every \
        snapshot")
    )

    if save_figures

        # Set the scale of the axes
        axes.xscale = (xscale_flag ? x_scale_func : identity)
        axes.yscale = (yscale_flag ? y_scale_func : identity)

        ############################################################################################
        # Apply the post processing function
        ############################################################################################

        pp_legend = post_processing(figure, pp_args...; pp_kwargs...)

        ############################################################################################
        # Add the legends
        ############################################################################################

        legend_elements = Vector{Makie.LegendElement}(undef, 0)
        legend_labels = Vector{Union{AbstractString,Nothing}}(undef, 0)

        if !isnothing(sim_labels)
            # Add the main legend
            (
                length(sim_labels) == n_simulations ||
                throw(ArgumentError("plotTimeSeries: The arguments `simulation_paths` and \
                `sim_labels` must have the same length, but I got length(`sim_labels`) = \
                $(length(sim_labels)) != length(`simulation_paths`) = $(n_simulations)"))
            )

            # Load the current palette
            colors     = current_theme[:palette][:color][]
            markers    = current_theme[:palette][:marker][]
            linestyles = current_theme[:palette][:linestyle][]

            for i in 1:n_simulations
                color = ring(colors, i)
                marker = ring(markers, i)
                linestyle = ring(linestyles, i)
                plot_function = ring(plot_functions, i)

                if plot_function == lines!
                    push!(legend_elements, LineElement(; color, linestyle))
                else
                    push!(legend_elements, MarkerElement(; color, marker))
                end
            end

            append!(legend_labels, sim_labels)

        end

        if !isnothing(pp_legend)
            # Add the post processing legend
            append!(legend_elements, pp_legend[1])
            append!(legend_labels, pp_legend[2])
        end

        if !any(isempty, [legend_elements, legend_labels])
            # Add a legend to the plot
            Makie.Legend(figure[1, 1], legend_elements, legend_labels)
        end

        ############################################################################################
        # Save the figure
        ############################################################################################

        save(joinpath(output_path, filename * output_format), figure)

    end

    return axes, figure

end
