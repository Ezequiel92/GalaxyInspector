####################################################################################################
# Pipeline functions
####################################################################################################

"""
    snapshotPlot(
        source_paths::Vector{String},
        base_names::Vector{String},
        quantities::Dict{Symbol,<:Union{String,Vector{String}}},
        plot_function!::Function; 
        <keyword arguments>
    )::Nothing

Generate one figure per snapshot, for one or more simulations.

Some of the features are:

- It can produce scatter plots, line plots, histograms, and density maps. The last two are only for 
  one simulation at a time.
- It can produce a GIF or an MP4 animating the results.
- It transparently manages units, you only have to indicate the final unit of each axis.

Many aspects of the plots can be configured, see the arguments below for details.

# Arguments
- `source_paths::Vector{String}`: Paths to the folders containing the snapshots of each simulation, 
  set in the GADGET variable `OutputDir`.
- `base_names::Vector{String}`: Base names of the snapshot files, set in the GADGET variable 
  `SnapshotFileBase`.
- `quantities::Dict{Symbol,<:Union{String,Vector{String}}}`: Dictionary where the keys are the 
  types of particles (the possibilities are given by [`ParticleType`](@ref) in `src/constants.jl`), 
  and the values are the data blocks to be read. Which data blocks are nedded depends on the
  provided function `data_analysis`.
- `plot_function!::Function`: A plotting function from [Makie.jl](http://makie.juliaplots.org/stable/index.html). 
  The supported functions are:
  - `hist!`         ⟶ Histograms.
  - `scatter!`      ⟶ Scatter plots.
  - `lines!`        ⟶ Line plots.
  - `scatterlines!` ⟶ Markers with lines between them.
  - `heatmap!`      ⟶ Heatmaps. 
  This sets the type of plot for the output.

## Data manipulation options
- `pf_kwargs::NamedTuple = (;)`: Keyword arguments for the `plot_function!` function.
- `sim_cosmo::Bool = false`: If the simulation is cosmological. 
  - `false` ⟶ Newtonian simulation (`ComovingIntegrationOn` = 0).
  - `true`  ⟶ Cosmological simulation (`ComovingIntegrationOn` = 1).
- `idx::IndexType = :`: Indexing of the simulation, i.e. indices of the snapshots that will be read.
  It can be an integer (a single snapshot), a vector of integers (several given snapshots), an
  UnitRange (e.g. 5:13) or : (every snapshot will be read).
- `filter_functions::Vector{<:Union{Function,Nothing}} = [nothing,]`: List of functions (one per 
  simulation), all with the signature: 

  `foo(file_path::String)::Vector{Int64}`
  
  Each function indicates which particles will be read in each simulation, taking the file path to 
  a snapshot and returning the list of indices of the selected particles. If set to `nothing`, then 
  every particle is read. See the [GadgetIO.jl](https://ludwigboess.github.io/GadgetIO.jl/stable/read_snapshots/#Filter-functions) documentation for examples.
- `data_analysis::Function = trivial`: Optional data analysis function. See the required signature 
  and examples in `src/data_analysis.jl`.
- `da_args::Tuple = ()`: Positional arguments for the `data_analysis` function.
- `da_kwargs::NamedTuple = (;)`: Keyword arguments for the `data_analysis` function.
- `post_processing!::Function = trivial`: Optional post processing function.
  See the required signature and examples in `src/post_processing.jl`.
- `pp_args::Tuple = ()`: Positional arguments for the `post_processing!` function.
- `pp_kwargs::NamedTuple = (;)`: Keyword arguments for the `post_processing!` function.

## Axes options
- `x_label::String = ":auto_label"`: Label for the x axis. It can contain the string `:auto_label`,
  which will be replaced by the automatic label that would otherwise be used. The default label has 
  the shape: `x_name` / 10^`x_factor` `x_unit`.
- `y_label::String = ":auto_label"`: Label for the y axis. It can contain the string `:auto_label`,
  which will be replaced by the automatic label that would otherwise be used. The default label has 
  the shape: `y_name` / 10^`y_factor` `y_unit`.
- `x_name::Union{String,LaTeXString} = ""`: Name of the variable for the x axis. It should not 
  include units or scaling factors. If left empty, the automatic label will be an empty string.
- `y_name::Union{String,LaTeXString} = ""`: Name of the variable for the y axis. It should not 
  include units or scaling factors. If left empty, the automatic label will be an empty string.
- `x_unit::Unitful.Units = Unitful.NoUnits`: Target unit for the x axis. The values will be 
  converted accordingly. Leave the default value of `Unitful.NoUnits` for dimensionless quantities.
- `y_unit::Unitful.Units = Unitful.NoUnits`: Target unit for the y axis. The values will be 
  converted accordingly. Leave the default value of `Unitful.NoUnits` for dimensionless quantities.
- `x_factor::Int64 = 0`: Numerical exponent to scale the x axis, e.g. if `x_factor` = 10 the values 
  will be scaled by ``10^{10}``. The default is no scaling.
- `y_factor::Int64 = 0`: Numerical exponent to scale the y axis, e.g. if `y_factor` = 10 the values 
  will be scaled by ``10^{10}``. The default is no scaling.
- `x_scale::Union{Function,Makie.Symlog10} = identity`: Scaling function for the x axis. 
  The options are the scaling functions accepted by [Makie.jl](http://makie.juliaplots.org/stable/index.html): log10, log2, log, sqrt, Makie.logit, 
  Makie.Symlog10, Makie.pseudolog10 and identity.
- `y_scale::Union{Function,Makie.Symlog10} = identity`: Scaling function for the y axis. 
  The options are the scaling functions accepted by [Makie.jl](http://makie.juliaplots.org/stable/index.html): log10, log2, log, sqrt, Makie.logit, 
  Makie.Symlog10, Makie.pseudolog10 and identity.
- `x_limits::Tuple{<:Union{Real,Nothing},<:Union{Real,Nothing}} = (nothing, nothing)`: Set it to a 
  value diferent than `nothing` if you want to fix the limits of the x axis.
- `y_limits::Tuple{<:Union{Real,Nothing},<:Union{Real,Nothing}} = (nothing, nothing)`: Set it to a 
  value diferent than `nothing` if you want to fix the limits of the y axis.
- `x_range::NTuple{2,<:Real} = (-Inf, Inf)`: The data will be trim down so its x coordinates fit 
  within `x_range`.
- `y_range::NTuple{2,<:Real} = (-Inf, Inf)`: The data will be trim down so its y coordinates fit 
   within `y_range`. Notice that this option does not affect histograms.
- `x_edges::Bool = false`: Set it to `true` if you want to keep the borders of `x_range`.
- `y_edges::Bool = false`: Set it to `true` if you want to keep the borders of `y_range`.
- `x_func::Function = identity`: Function to be applied to the values of the x axis. It has to be a 
  pure function with the signature:
            
  `foo(x_data::Vector{Float64})::Vector{Float64}`
            
  And the output vector has to have the same length as the input one.
  Notice that this function will be applied regardless of units and possible domain problems (use 
  `x_range` to solve incompatibilities), and that it will not be reflected in the automatic label.
- `y_func::Function = identity`: Function to be applied to the values of the y axis. It has to be a 
  pure function with the signature:
              
  `foo(y_data::Vector{Float64})::Vector{Float64}`
            
  And the output vector has to have the same length as the input one.
  Notice that this function will be applied regardless of units and possible domain problems (use 
  `y_range` to solve incompatibilities), and that it will not be reflected in the automatic label.

## Plotting and animation options
- `output_path::String = "snapshot_plots"`: Path to the output folder. The figures will be stored 
  in `output_path`/plots/ and will be named snapshot\\_XXX`output_format` where XXX is the snapshot 
  number.
- `sim_labels::Union{Vector{String},Nothing} = basename.(source_paths)`: Labels for the plot legend,
  one per simulation. Set it to `nothing` if you don't want a legend.
- `title::Union{String,Symbol} = ""`: Title for the figures. If left empty, no title will be printed. 
  If set to `:time`, the time stamp of the snapshot will be used as title. 
- `t_unit::Unitful.Units = UnitfulAstro.Myr`: Target unit for the time stamps of the snapshots (only
  relevant if the option `:time` was given for the argument `title`).
- `output_format::String = ".svg"`: File format of the output figures. All formats supported by
  [Makie.jl](http://makie.juliaplots.org/stable/index.html) can be used, namely `.pdf`, `.svg` and `.png`. 
- `pt_per_unit::Float64 = 0.75`: Factor to scale up or down the size of the figures, keeping their 
  proportions. It only works for `.pdf` and `.svg`.
- `px_per_unit::Float64 = 1.0`: Factor to scale up or down the size of the figures, keeping their 
  proportions. It only works for `.png`.
- `resolution::NTuple{2,Int64} = (1000, 750)`: Resolution of the figures in points. By default, for 
  PNGs: points = pixels (as given by `px_per_unit` = 1.0), and for PDFs and SVGs:
  points = 0.75 * pixels (as given by `pt_per_unit` = 0.75).
- `aspect::Union{DataAspect,AxisAspect,Nothing} = nothing`: Aspect ratio of the figures.
  The options are:
  - `nothing`: Use the default given by [Makie.jl](http://makie.juliaplots.org/stable/index.html).
  - `AxisAspect(n)`: The aspect ratio will be given by the number `n` = width / height.
  - `DataAspect()`: The aspect ratio of the data is used.
- `hist_bins::Int64 = 20`: Number of bins for the histograms.
- `series_colors::Union{Vector{ColorInput},Nothing} = nothing`: Colors for the diferent series 
  to be plotted (diferent simulations). If set to `nothing`, the colors will be assigned 
  automatically. This is only relevant for `scatter!` and `scatterlines!` plots.
- `series_markers::Union{Vector{Symbol},Nothing} = nothing`: Markers for the different series 
  to be plotted (different simulations). If set to `nothing`, the markers will be assigned 
  automatically. This is only relevant for `scatter!` and `scatterlines!` plots.
- `series_linetypes::Union{Vector{LineStyleInput},Nothing} = nothing`: Line types for the different 
  series to be plotted (different simulations). If set to `nothing`, the line types will be assigned 
  automatically. This is only relevant for `lines!` and `scatterlines!` plots.
- `animation::Bool = true`: If an animation of the results will be made and saved in `output_path`. 
- `anim_file::String = "animation.gif"`: Filename for the animation, including its extension.
  All formats supported by [Makie.jl](http://makie.juliaplots.org/stable/index.html) can be used, namely `.mkv`, `.mp4`, `.webm` and `.gif`.
- `framerate::Int64 = 10`: Frame rate of the animation.
- `warnings::Bool = true`: If a warning will be raised when the data is not as expected, but the 
  function can still run using sane defaults.
- `show_progress::Bool = true`: If a progress bar will be shown.

"""
function snapshotPlot(
    source_paths::Vector{String},
    base_names::Vector{String},
    quantities::Dict{Symbol,<:Union{String,Vector{String}}},
    plot_function!::Function;
    pf_kwargs::NamedTuple = (;),
    sim_cosmo::Bool = false,
    idx::IndexType = :,
    filter_functions::Vector{<:Union{Function,Nothing}} = [nothing,],
    data_analysis::Function = trivial,
    da_args::Tuple = (),
    da_kwargs::NamedTuple = (;),
    post_processing!::Function = trivial,
    pp_args::Tuple = (),
    pp_kwargs::NamedTuple = (;),
    x_label::String = ":auto_label",
    y_label::String = ":auto_label",
    x_name::Union{String,LaTeXString} = "",
    y_name::Union{String,LaTeXString} = "",
    x_unit::Unitful.Units = Unitful.NoUnits,
    y_unit::Unitful.Units = Unitful.NoUnits,
    x_factor::Int64 = 0,
    y_factor::Int64 = 0,
    x_scale::Union{Function,Makie.Symlog10} = identity,
    y_scale::Union{Function,Makie.Symlog10} = identity,
    x_limits::Tuple{<:Union{Real,Nothing},<:Union{Real,Nothing}} = (nothing, nothing),
    y_limits::Tuple{<:Union{Real,Nothing},<:Union{Real,Nothing}} = (nothing, nothing),
    x_range::NTuple{2,<:Real} = (-Inf, Inf),
    y_range::NTuple{2,<:Real} = (-Inf, Inf),
    x_edges::Bool = false,
    y_edges::Bool = false,
    x_func::Function = identity,
    y_func::Function = identity,
    output_path::String = "snapshot_plots",
    sim_labels::Union{Vector{String},Nothing} = basename.(source_paths),
    title::Union{String,Symbol} = "",
    t_unit::Unitful.Units = UnitfulAstro.Myr,
    output_format::String = ".svg",
    pt_per_unit::Float64 = 0.75,
    px_per_unit::Float64 = 1.0,
    resolution::NTuple{2,Int64} = (1000, 750),
    aspect::Union{DataAspect,AxisAspect,Nothing} = nothing,
    hist_bins::Int64 = 20,
    series_colors::Union{Vector{ColorInput},Nothing} = nothing,
    series_markers::Union{Vector{Symbol},Nothing} = nothing,
    series_linetypes::Union{Vector{LineStyleInput},Nothing} = nothing,
    animation::Bool = true,
    anim_file::String = "animation.gif",
    framerate::Int64 = 10,
    warnings::Bool = true,
    show_progress::Bool = true,
)::Nothing

    # If it doesn't exist create a folder to store the plots
    output_folder = mkpath(joinpath(output_path, "plots"))

    length(base_names) == length(source_paths) || throw(ArgumentError(
        "The arguments `base_names` and `source_paths` must have the same length."
    ))

    ################################################################################################
    # Load the simulation results
    ################################################################################################

    # Get the path and the global index of each snapshot, for every simulation
    source = getSnapshotPaths.(base_names, source_paths)

    # Make a dataframe where every row is a snapshot, and the columns are:
    #   - Global index
    #   - Clock time
    #   - Number in the filename of the snapshot
    #   - Path to the snapshot (or `missing` is the snapshot doesn't exist), for the 1st simulation
    #   - Path to the snapshot (or `missing` is the snapshot doesn't exist), for the 2nd simulation
    #   - Path to the snapshot (or `missing` is the snapshot doesn't exist), for the 3rd simulation
    #   - ... 
    sim_df = makeSourceTable(source, idx, t_unit, sim_cosmo)
    iterator = enumerate(eachrow(sim_df))

    # Number of simulations
    n_simulation = length(source_paths)
    # Largest number of snapshots in a simulation
    max_n_snapshot = length(iterator)

    # Save the initial scale factor (only relevant for cosmological simulations), assuming that 
    # every simulation has the same value
    a0 = read_header(sim_df[1, 4]).time

    ################################################################################################
    # Set up the canvas for the figures
    ################################################################################################

    # Create a list of distinguishable colors for the different simulations
    colors = distinguishable_colors(
        n_simulation,
        [RGB(1, 1, 1), RGB(0, 0, 0)],
        dropseed = true,
    )
    # Select a set of markers for the different simulation 
    # Only relevant for `scatter` and `scatterlines` plots
    markers = MARKERS[1:n_simulation]
    # Select a set of line styles for the different simulations 
    # Only relevant for `lines` and `scatterlines` plots
    linestyles = LINE_STYLES[1:n_simulation]

    # Iterator with the colors, markers, and line styles
    styles = zip(
        (series_colors !== nothing ? series_colors : colors),
        (series_markers !== nothing ? series_markers : markers),
        (series_linetypes !== nothing ? series_linetypes : linestyles),
    )
    n_simulation == length(styles) || throw(ArgumentError(
        "The arguments `series_colors`, `series_markers` and `series_linetypes` must have the same \
        length as the number of simulations, when set to a value different than `nothing`."
    ))

    # Apply the global theme defined in `src/constants.jl`
    set_theme!(theme)

    # Create an empty figure
    figure = Figure(; resolution)

    # Create the labels for the axes
    xlabel = LaTeXString(
        replace(x_label, ":auto_label" => getLabel(x_name, x_factor, x_unit)),
    )
    ylabel = LaTeXString(
        replace(y_label, ":auto_label" => getLabel(y_name, y_factor, y_unit)),
    )

    # Set the axes properties
    axes = Axis(
        figure[1, 1];
        limits = (x_limits, y_limits),
        xlabel,
        ylabel,
        aspect
    )

    ################################################################################################
    # Set up the animation
    ################################################################################################

    !(nrow(sim_df) < framerate && warnings && animation) || @warn(
        "With the settings: `framerate` = $framerate and `idx` = $idx, the animation is less than 
        a second long."
    )

    # Initialize the animation stream
    if animation
        vs = VideoStream(figure; framerate)
    end

    ################################################################################################
    # Main loop
    ################################################################################################

    # Initialize the progress bar
    prog_bar = Progress(
        max_n_snapshot,
        dt = 0.5,
        desc = "Analyzing and plotting the data... ",
        color = :blue,
        barglyphs = BarGlyphs("|#  |"),
        enabled = show_progress,
    )

    for (local_index, (global_index, clock_time, number, snapshots...)) in iterator

        # Flag to keep the x axis with a linear scale if there are less than 3 data points
        #   - `true`: The axis will have a linear scale
        #   - `false`: The axis will have the scale given by `x_scale`
        xscale_flag = false
        # Flag to keep the y axis with a linear scale if there are less than 3 data points
        #   - `true`: The axis will have a linear scale
        #   - `false`: The axis will have the scale given by `y_scale`
        yscale_flag = false
        # Flag to skip problematic snapshots
        #   - `true`: The snapshot will be skipped.
        #   - `false`: The snapshot will be plotted.
        skipper = true

        for ((sim_index, snapshot), style) in zip(enumerate(snapshots), styles)

            # Skip missing snapshots
            snapshot !== missing || continue

            # Filter function for the current simulation in the loop
            filter_function = filter_functions[mod1(sim_index, length(filter_functions))]

            # Store basic information about the current snapshot and simulation in the loop
            metadata = Dict(
                :sim_data => SimData(
                    source_paths[sim_index],  # Path to the folder with all the snapshot files
                    sim_index,                # Index associated with the simulation
                    idx,                      # Which snapshots will be read
                    base_names[sim_index],    # Base name of the snapshot files
                    read_header(snapshot),    # Snapshot header of the simulation
                    filter_function,          # Filter function for the current simulation
                    a0,                       # Initial scale factor for the current simulation
                ),
                :snap_data => SnapData(
                    snapshot,                 # Full path
                    global_index,             # Global index 
                    local_index,              # Local index 
                    clock_time,               # Time stamp 
                ),
            )

            # Extract the relevant data from the snapshot
            raw_data = merge(
                getSnapshotData(snapshot, quantities; sim_cosmo, filter_function, warnings),
                metadata,
            )

            # Apply the analysis function
            da_output = data_analysis(raw_data, da_args...; da_kwargs...)

            # Skip this snapshot if `data_analysis` returns `nothing`
            da_output === nothing ? continue : skipper = false

            # Data shape validation
            checkDataShape(plot_function!, length(da_output), n_simulation)

            # Unit conversion
            axis_data = [
                ustrip.(unit, ax_data) for
                (unit, ax_data) in zip([x_unit, y_unit, Unitful.NoUnits], da_output)
            ]

            # Sanitize and transform the data to fit the user's requirements
            if length(axis_data) == 1

                flag_x, _ = sanitizeData!(
                    axis_data[1],
                    func_domain = x_scale,
                    range = x_range,
                    keep_edges = x_edges,
                    min_left = 3,
                    exp_factor = x_factor,
                )
                flag_y = 0

                axis_data[1] = x_func(axis_data[1])

            else

                flag_x, flag_y = sanitizeData!(
                    axis_data[1],
                    axis_data[2],
                    func_domain = (x_scale, y_scale),
                    range = (x_range, y_range),
                    keep_edges = (x_edges, y_edges),
                    min_left = 3,
                    exp_factor = (x_factor, y_factor),
                )

                axis_data[1] = x_func(axis_data[1])
                axis_data[2] = y_func(axis_data[2])

            end
            # If, in the current snapshot and for any simulation, filtering the data targeting a 
            # non linear scale would leave less that three data points, the scale will revert to
            # `identity` and the data will be left unmodified
            if flag_x == 1
                xscale_flag = true
            end
            if flag_y == 1
                yscale_flag = true
            end

            # Draw the plot
            plot_function!(
                axes,
                axis_data...;
                bins = scaledBins(axis_data[1], hist_bins, x_scale; limits = x_limits, warnings),
                color = style[1],
                marker = style[2],
                linestyle = style[3],
                pf_kwargs...
            )

        end

        !skipper || (next!(prog_bar); continue)

        # Set the scale of the axes
        axes.xscale = xscale_flag ? identity : x_scale
        axes.yscale = yscale_flag ? identity : y_scale

        # Draw the title
        if title == :time
            # Use the current snapshot time-stamp
            time_stamp = string(round(typeof(clock_time), clock_time, sigdigits = 4))
            axes.title = "t = $time_stamp"
        else
            # Use th user provided title
            axes.title = title
        end

        # Draw the legend if required
        if sim_labels !== nothing

            length(sim_labels) == n_simulation || throw(ArgumentError(
                "The arguments `source_paths` and `sim_labels` must have the same length."
            ))

            if plot_function! == lines!
                legend_element = [
                    LineElement(; color, linestyle) for
                    (color, linestyle) in zip(colors, linestyles)
                ]
            elseif plot_function! == scatter!
                legend_element = [
                    MarkerElement(; color, marker) for
                    (color, marker) in zip(colors, markers)
                ]
            elseif plot_function! == scatterlines!
                legend_element = [
                    [MarkerElement(; color, marker), LineElement(; color, linestyle)] for
                    (color, marker, linestyle) in zip(colors, markers, linestyles)
                ]
            end

            Legend(figure[2, 1][1, 1], legend_element, sim_labels)

        end

        # Apply the post processing function
        pp_legend = post_processing!(figure, pp_args...; pp_kwargs...)

        # Add the post processing legend if there is one
        if pp_legend !== nothing
            Legend(
                figure[2, 1][1, sim_labels !== nothing ? 2 : 1],
                pp_legend[1],
                pp_legend[2],
            )
        end

        # Save the figure
        save(
            joinpath(output_folder, "snapshot_" * number * output_format),
            figure;
            pt_per_unit,
            px_per_unit
        )

        if animation
            # Add the figure as a frame to the animation stream
            recordframe!(vs)
        end

        # Clean the canvas for the next step in the loop
        cleanPlot!(figure)

        # Advance the progress bar
        next!(prog_bar)

    end

    if animation
        # Save the animation
        save(joinpath(output_path, anim_file), vs)
    end

    return nothing

end

"""
    timeSeriesPlot(
        source_paths::Vector{String},
        base_names::Vector{String},
        plot_function!::Function; 
        <keyword arguments>
    )::Nothing

Generate a figure with the time evolution of some quantity, for one or more simulations.

Some of the features are:

- It can produce a scatter plot, a line plot, or a combination of both.
- It transparently manages units, you only have to indicate the final unit of each axis.

Many aspects of the plot can be configured, see the arguments below for details.

# Arguments
- `source_paths::Vector{String}`: Paths to the folders containing the snapshots of each simulation, 
  set in the GADGET variable `OutputDir`.
- `base_names::Vector{String}`: Base names of the snapshot files, set in the GADGET variable 
  `SnapshotFileBase`.
- `plot_function!::Function`: A plotting function from [Makie.jl](http://makie.juliaplots.org/stable/index.html). 
  The supported functions are:
  - `scatter!`      ⟶ Scatter plots.
  - `lines!`        ⟶ Line plots.
  - `scatterlines!` ⟶ Markers with lines between them.
  This sets the type of plot for the output.

## Data manipulation options
- `pf_kwargs::NamedTuple = (;)`: Keyword arguments for the `plot_function!` function.
- `idx::IndexType = :`: Indexing of the simulation, i.e. indices of the snapshots that will be read.
  It can be an integer (a single snapshot), a vector of integers (several given snapshots), an
  UnitRange (e.g. 5:13) or : (every snapshot will be read).
- `filter_functions::Vector{<:Union{Function,Nothing}} = [nothing,]`: List of functions (one per 
  simulation), all with the signature: 

  `foo(file_path::String)::Vector{Int64}`
  
  Each function indicates which particles will be read in each simulation, taking the file path to 
  a snapshot and returning the list of indices of the selected particles. If set to `nothing`, then 
  every particle is read. See the [GadgetIO.jl](https://ludwigboess.github.io/GadgetIO.jl/stable/read_snapshots/#Filter-functions) documentation for examples.
- `data_analysis::Function = trivial`: Optional data analysis function. See the required signature 
  and examples in `src/data_analysis.jl`.
- `da_args::Tuple = ()`: Positional arguments for the `data_analysis` function.
- `da_kwargs::NamedTuple = (;)`: Keyword arguments for the `data_analysis` function.
- `post_processing!::Function = trivial`: Optional post processing function.
  See the required signature and examples in `src/post_processing.jl`.
- `pp_args::Tuple = ()`: Positional arguments for the `post_processing!` function.
- `pp_kwargs::NamedTuple = (;)`: Keyword arguments for the `post_processing!` function.

## Axes options
- `x_label::String = ":auto_label"`: Label for the x axis. It can contain the string `:auto_label`,
  which will be replaced by the automatic label that would otherwise be used. The default label has 
  the shape: `x_name` / 10^`x_factor` `x_unit`.
- `y_label::String = ":auto_label"`: Label for the y axis. It can contain the string `:auto_label`,
  which will be replaced by the automatic label that would otherwise be used. The default label has 
  the shape: `y_name` / 10^`y_factor` `y_unit`.
- `x_name::Union{String,LaTeXString} = ""`: Name of the variable for the x axis. It should not 
  include units or scaling factors. If left empty, the automatic label will be an empty string.
- `y_name::Union{String,LaTeXString} = ""`: Name of the variable for the y axis. It should not 
  include units or scaling factors. If left empty, the automatic label will be an empty string.
- `x_unit::Unitful.Units = Unitful.NoUnits`: Target unit for the x axis. The values will be 
  converted accordingly. Leave the default value of `Unitful.NoUnits` for dimensionless quantities.
- `y_unit::Unitful.Units = Unitful.NoUnits`: Target unit for the y axis. The values will be 
  converted accordingly. Leave the default value of `Unitful.NoUnits` for dimensionless quantities.
- `x_factor::Int64 = 0`: Numerical exponent to scale the x axis, e.g. if `x_factor` = 10 the values 
  will be scaled by ``10^{10}``. The default is no scaling.
- `y_factor::Int64 = 0`: Numerical exponent to scale the y axis, e.g. if `y_factor` = 10 the values 
  will be scaled by ``10^{10}``. The default is no scaling.
- `x_scale::Union{Function,Makie.Symlog10} = identity`: Scaling function for the x axis. 
  The options are the scaling functions accepted by [Makie.jl](http://makie.juliaplots.org/stable/index.html): log10, log2, log, sqrt, Makie.logit, 
  Makie.Symlog10, Makie.pseudolog10 and identity.
- `y_scale::Union{Function,Makie.Symlog10} = identity`: Scaling function for the y axis. 
  The options are the scaling functions accepted by [Makie.jl](http://makie.juliaplots.org/stable/index.html): log10, log2, log, sqrt, Makie.logit, 
  Makie.Symlog10, Makie.pseudolog10 and identity.
- `x_limits::Tuple{<:Union{Real,Nothing},<:Union{Real,Nothing}} = (nothing, nothing)`: Set it to a 
  value diferent than `nothing` if you want to fix the limits of the x axis.
- `y_limits::Tuple{<:Union{Real,Nothing},<:Union{Real,Nothing}} = (nothing, nothing)`: Set it to a 
  value diferent than `nothing` if you want to fix the limits of the y axis.
- `x_range::NTuple{2,<:Real} = (-Inf, Inf)`: The data will be trim down so its x coordinates fit 
  within `x_range`.
- `y_range::NTuple{2,<:Real} = (-Inf, Inf)`: The data will be trim down so its y coordinates fit 
   within `y_range`. Notice that this option does not affect histograms.
- `x_edges::Bool = false`: Set it to `true` if you want to keep the borders of `x_range`.
- `y_edges::Bool = false`: Set it to `true` if you want to keep the borders of `y_range`.
- `x_func::Function = identity`: Function to be applied to the values of the x axis. It has to be a 
  pure function with the signature:
            
  `foo(x_data::Vector{Float64})::Vector{Float64}`
            
  And the output vector has to have the same length as the input one.
  Notice that this function will be applied regardless of units and possible domain problems (use 
  `x_range` to solve incompatibilities), and that it will not be reflected in the automatic label.
- `y_func::Function = identity`: Function to be applied to the values of the y axis. It has to be a 
  pure function with the signature:
              
  `foo(y_data::Vector{Float64})::Vector{Float64}`
            
  And the output vector has to have the same length as the input one.
  Notice that this function will be applied regardless of units and possible domain problems (use 
  `y_range` to solve incompatibilities), and that it will not be reflected in the automatic label.

## Plotting and animation options
- `output_path::String = "time_series_plots""`: Path to the output folder.
- `file_name::String = "figure.svg"`: Filename for the output figure with its extension. 
  All formats supported by [Makie.jl](http://makie.juliaplots.org/stable/index.html) can be used, namely `.pdf`, `.svg` and `.png`. 
- `sim_labels::Union{Vector{String},Nothing} = basename.(source_paths)`: Labels for the plot legend,
  one per simulation. Set it to `nothing` if you don't want a legend.
- `title::String = ""`: Title for the figures. If left empty, no title will be printed. 
- `pt_per_unit::Float64 = 0.75`: Factor to scale up or down the size of the figures, keeping their 
  proportions. It only works for `.pdf` and `.svg`.
- `px_per_unit::Float64 = 1.0`: Factor to scale up or down the size of the figures, keeping their 
  proportions. It only works for `.png`.
- `resolution::NTuple{2,Int64} = (1000, 750)`: Resolution of the figures in points. By default, for 
  PNGs: points = pixels (as given by `px_per_unit` = 1.0), and for PDFs and SVGs:
  points = 0.75 * pixels (as given by `pt_per_unit` = 0.75).
- `aspect::Union{DataAspect,AxisAspect,Nothing} = nothing`: Aspect ratio of the figures.
  The options are:
  - `nothing`: Use the default given by [Makie.jl](http://makie.juliaplots.org/stable/index.html).
  - `AxisAspect(n)`: The aspect ratio will be given by the number `n` = width / height.
  - `DataAspect()`: The aspect ratio of the data is used.
- `series_colors::Union{Vector{ColorInput},Nothing} = nothing`: Colors for the diferent series 
  to be plotted (diferent simulations). If set to `nothing`, the colors will be assigned 
  automatically. This is only relevant for `scatter!` and `scatterlines!` plots.
- `series_markers::Union{Vector{Symbol},Nothing} = nothing`: Markers for the different series 
  to be plotted (different simulations). If set to `nothing`, the markers will be assigned 
  automatically. This is only relevant for `scatter!` and `scatterlines!` plots.
- `series_linetypes::Union{Vector{LineStyleInput},Nothing} = nothing`: Line types for the different 
  series to be plotted (different simulations). If set to `nothing`, the line types will be assigned 
  automatically. This is only relevant for `lines!` and `scatterlines!` plots.

"""
function timeSeriesPlot(
    source_paths::Vector{String},
    base_names::Vector{String},
    plot_function!::Function;
    pf_kwargs::NamedTuple = (;),
    idx::IndexType = :,
    filter_functions::Vector{<:Union{Function,Nothing}} = [nothing,],
    data_analysis::Function = trivial,
    da_args::Tuple = (),
    da_kwargs::NamedTuple = (;),
    post_processing!::Function = trivial,
    pp_args::Tuple = (),
    pp_kwargs::NamedTuple = (;),
    x_label::String = ":auto_label",
    y_label::String = ":auto_label",
    x_name::Union{String,LaTeXString} = "",
    y_name::Union{String,LaTeXString} = "",
    x_unit::Unitful.Units = Unitful.NoUnits,
    y_unit::Unitful.Units = Unitful.NoUnits,
    x_factor::Int64 = 0,
    y_factor::Int64 = 0,
    x_scale::Union{Function,Makie.Symlog10} = identity,
    y_scale::Union{Function,Makie.Symlog10} = identity,
    x_limits::Tuple{<:Union{Real,Nothing},<:Union{Real,Nothing}} = (nothing, nothing),
    y_limits::Tuple{<:Union{Real,Nothing},<:Union{Real,Nothing}} = (nothing, nothing),
    x_range::NTuple{2,<:Real} = (-Inf, Inf),
    y_range::NTuple{2,<:Real} = (-Inf, Inf),
    x_edges::Bool = false,
    y_edges::Bool = false,
    x_func::Function = identity,
    y_func::Function = identity,
    output_path::String = "snapshot_plots",
    file_name::String = "figure.svg",
    sim_labels::Union{Vector{String},Nothing} = basename.(source_paths),
    title::String = "",
    pt_per_unit::Float64 = 0.75,
    px_per_unit::Float64 = 1.0,
    resolution::NTuple{2,Int64} = (1000, 750),
    aspect::Union{DataAspect,AxisAspect,Nothing} = nothing,
    series_colors::Union{Vector{ColorInput},Nothing} = nothing,
    series_markers::Union{Vector{Symbol},Nothing} = nothing,
    series_linetypes::Union{Vector{LineStyleInput},Nothing} = nothing
)::Nothing

    length(base_names) == length(source_paths) || throw(ArgumentError(
        "The arguments `base_names` and `source_paths` must have the same length."
    ))

    # Number of simulations
    n_simulation = length(source_paths)

    ################################################################################################
    # Set up the canvas for the figures
    ################################################################################################

    # Create a list of distinguishable colors for the different simulations
    colors = distinguishable_colors(
        n_simulation,
        [RGB(1, 1, 1), RGB(0, 0, 0)],
        dropseed = true,
    )
    # Select a set of markers for the different simulation 
    # Only relevant for `scatter` and `scatterlines` plots
    markers = MARKERS[1:n_simulation]
    # Select a set of line styles for the different simulations 
    # Only relevant for `lines` and `scatterlines` plots
    linestyles = LINE_STYLES[1:n_simulation]

    # Iterator with the colors, markers, and line styles
    styles = zip(
        (series_colors !== nothing ? series_colors : colors),
        (series_markers !== nothing ? series_markers : markers),
        (series_linetypes !== nothing ? series_linetypes : linestyles),
    )
    n_simulation == length(styles) || throw(ArgumentError(
        "The arguments `series_colors`, `series_markers` and `series_linetypes` must have the same \
        length as the number of simulations, when set to a value different than `nothing`."
    ))

    # Apply the global theme defined in `src/constants.jl`
    set_theme!(theme)

    # Create an empty figure
    figure = Figure(; resolution)

    # Create the labels for the axes
    xlabel = LaTeXString(
        replace(x_label, ":auto_label" => getLabel(x_name, x_factor, x_unit)),
    )
    ylabel = LaTeXString(
        replace(y_label, ":auto_label" => getLabel(y_name, y_factor, y_unit)),
    )

    # Set the axes properties
    axes = Axis(
        figure[1, 1];
        limits = (x_limits, y_limits),
        xlabel,
        ylabel,
        title,
        aspect
    )

    ################################################################################################
    # Main loop
    ################################################################################################

    iterator = enumerate(zip(source_paths, base_names, styles))

    # Flag to keep the x axis with a linear scale if there are less than 3 data points
    #   - `true`: The axis will have a linear scale
    #   - `false`: The axis will have the scale given by `x_scale`
    xscale_flag = false
    # Flag to keep the y axis with a linear scale if there are less than 3 data points
    #   - `true`: The axis will have a linear scale
    #   - `false`: The axis will have the scale given by `y_scale`
    yscale_flag = false

    for (sim_index, (source_path, base_name, style)) in iterator

        # Filter function for the current simulation in the loop
        filter_function = filter_functions[mod1(sim_index, length(filter_functions))]

        # Header of the first snapshot
        header = read_header(getSnapshotPaths(base_name, source_path)["snap_paths"][1])

        # Store basic information about the current simulation in the loop
        sim_data = SimData(
            source_path,        # Path to the folder with all the snapshot files
            sim_index,          # Index associated with the simulation
            idx,                # Which snapshots will be read
            base_name,          # Base name of the snapshot files
            header,             # Snapshot header of the simulation
            filter_function,    # Filter function for the simulation
            header.time,        # Initial scale factor for the simulation
        )

        # Apply the analysis function
        da_output = data_analysis(sim_data, da_args...; da_kwargs...)

        # Skip this simulation if `data_analysis` returns `nothing`
        if da_output === nothing
            continue
        end

        # Unit conversion
        axis_data = [ustrip.(unit, ax_data) for (unit, ax_data) in zip([x_unit, y_unit], da_output)]

        # Sanitize and transform the data to fit the user's requirements
        flag_x, flag_y = sanitizeData!(
            axis_data[1],
            axis_data[2],
            func_domain = (x_scale, y_scale),
            range = (x_range, y_range),
            keep_edges = (x_edges, y_edges),
            min_left = 3,
            exp_factor = (x_factor, y_factor),
        )
        axis_data[1] = x_func(axis_data[1])
        axis_data[2] = y_func(axis_data[2])
        # If, in the current snapshot and for any simulation, filtering the data targeting a 
        # non linear scale would leave less that three data points, the scale will revert to
        # `identity` and the data will be left unmodified
        if flag_x == 1
            xscale_flag = true
        end
        if flag_y == 1
            yscale_flag = true
        end

        # Draw the plot
        plot_function!(
            axes,
            axis_data...;
            color = style[1],
            marker = style[2],
            linestyle = style[3],
            pf_kwargs...
        )

    end

    # Set the scale of the axes
    axes.xscale = xscale_flag ? identity : x_scale
    axes.yscale = yscale_flag ? identity : y_scale

    # Draw the legend if required
    if sim_labels !== nothing

        length(sim_labels) == n_simulation || throw(ArgumentError(
            "The arguments `source_paths` and `sim_labels` must have the same length."
        ))

        if plot_function! == lines!
            legend_element = [
                LineElement(; color, linestyle) for
                (color, linestyle) in zip(colors, linestyles)
            ]
        elseif plot_function! == scatter!
            legend_element = [
                MarkerElement(; color, marker) for
                (color, marker) in zip(colors, markers)
            ]
        elseif plot_function! == scatterlines!
            legend_element = [
                [MarkerElement(; color, marker), LineElement(; color, linestyle)] for
                (color, marker, linestyle) in zip(colors, markers, linestyles)
            ]
        end

        Legend(figure[2, 1][1, 1], legend_element, sim_labels)

    end

    # Apply the post processing function
    pp_legend = post_processing!(figure, pp_args...; pp_kwargs...)

    # Add the post processing legend if there is one
    if pp_legend !== nothing
        Legend(
            figure[2, 1][1, sim_labels !== nothing ? 2 : 1],
            pp_legend[1],
            pp_legend[2],
        )
    end

    # If it doesn't exist create a folder to store the plots
    fig_path = mkpath(output_path)

    # Save the figure
    save(joinpath(fig_path, file_name), figure; pt_per_unit, px_per_unit)

    return nothing

end

"""
    snapshotTable(
        source_paths::Vector{String},
        base_names::Vector{String},
        quantities::Dict{Symbol,<:Union{String,Vector{String}}}; 
        <keyword arguments>
    )::Nothing

Generate one table per snapshot, for one or more simulations.

It transparently manages units, you only have to indicate the final unit of each axis.

The output can be a human-readable TXT or a CSV. In both cases the columns in each snapshot table 
are:

- x quantity for simulation 1
- y quantity for simulation 1
- x quantity for simulation 2
- y quantity for simulation 2
- x quantity for simulation 3
- y quantity for simulation 3

etc.

# Arguments
- `source_paths::Vector{String}`: Paths to the folders containing the snapshots of each simulation, 
  set in the GADGET variable `OutputDir`.
- `base_names::Vector{String}`: Base names of the snapshot files, set in the GADGET variable 
  `SnapshotFileBase`.
- `quantities::Dict{Symbol,<:Union{String,Vector{String}}}`: Dictionary where the keys are the 
  types of particles (the possibilities are given by [`ParticleType`](@ref) in `src/constants.jl`), 
  and the values are the data blocks to be read. Which data blocks are nedded depends on the
  provided function `data_analysis`.
- `human_redable::Bool = true`: If the results will be written in a human readable format. The 
  output file in this case will be a `.txt` file.
- `machine_redable::Bool = true`: If the results will be written in a machine readable format. The 
  output file in this case will be a `.csv` file.
- `output_path::String = "snapshot_tables"`:  Path to the output folder.
- `sim_labels::Vector{String}} = basename.(source_paths)`: Labels for the simulations.
- `title::String = ""`: Title to be written at the top of the file, only aplicable to the `.txt`
  file.
- `warnings::Bool = true`: If a warning will be raised when the data is not as expected, but the 
  function can still run using sane defaults.
- `show_progress::Bool = true`: If a progress bar will be shown.
- `sim_cosmo::Bool = false`: If the simulation is cosmological. 
  - `false` ⟶ Newtonian simulation (`ComovingIntegrationOn` = 0).
  - `true`  ⟶ Cosmological simulation (`ComovingIntegrationOn` = 1).

## Data manipulation options
- `idx::IndexType = :`: Indexing of the simulation, i.e. indices of the snapshots that will be read.
  It can be an integer (a single snapshot), a vector of integers (several given snapshots), an
  UnitRange (e.g. 5:13) or : (every snapshot will be read).
- `filter_functions::Vector{<:Union{Function,Nothing}} = [nothing,]`: List of functions (one per 
  simulation), all with the signature: 

  `foo(file_path::String)::Vector{Int64}`
  
  Each function indicates which particles will be read in each simulation, taking the file path to 
  a snapshot and returning the list of indices of the selected particles. If set to `nothing`, then 
  every particle is read. See the [GadgetIO.jl](https://ludwigboess.github.io/GadgetIO.jl/stable/read_snapshots/#Filter-functions) documentation for examples.
- `data_analysis::Function = trivial`: Optional data analysis function. See the required signature 
  and examples in `src/data_analysis.jl`.
- `da_args::Tuple = ()`: Positional arguments for the `data_analysis` function.
- `da_kwargs::NamedTuple = (;)`: Keyword arguments for the `data_analysis` function.

## Axes options
- `x_name::Union{String,LaTeXString} = "X"`: Name of the variable for the x axis. It should not 
  include units or scaling factors.
- `y_name::Union{String,LaTeXString} = "Y"`: Name of the variable for the y axis. It should not 
  include units or scaling factors.
- `xunit_label::String = "[:auto_label]"`: Label indicating the unit of the x axis. It can contain 
  the string ":auto_label", which will be replaced by the automatic unit label that would otherwise 
  be used. 
- `yunit_label::String = "[:auto_label]"`: Label indicating the unit of the x axis. It can contain 
  the string ":auto_label", which will be replaced by the automatic unit label that would otherwise 
  be used.
- `x_unit::Unitful.Units = Unitful.NoUnits`: Target unit for the x axis. The values will be 
  converted accordingly. Leave the default value of `Unitful.NoUnits` for dimensionless quantities.
- `y_unit::Unitful.Units = Unitful.NoUnits`: Target unit for the y axis. The values will be 
  converted accordingly. Leave the default value of `Unitful.NoUnits` for dimensionless quantities.
- `x_factor::Int64 = 0`: Numerical exponent to scale the x axis, e.g. if `x_factor` = 10 the values 
  will be scaled by ``10^{10}``. The default is no scaling.
- `y_factor::Int64 = 0`: Numerical exponent to scale the y axis, e.g. if `y_factor` = 10 the values 
  will be scaled by ``10^{10}``. The default is no scaling.
- `x_range::NTuple{2,<:Real} = (-Inf, Inf)`: The data will be trim down so its x coordinates fit 
  within `x_range`.
- `y_range::NTuple{2,<:Real} = (-Inf, Inf)`: The data will be trim down so its y coordinates fit 
   within `y_range`. Notice that this option does not affect histograms.
- `x_edges::Bool = false`: Set it to `true` if you want to keep the borders of `x_range`.
- `y_edges::Bool = false`: Set it to `true` if you want to keep the borders of `y_range`.
- `x_func::Function = identity`: Function to be applied to the values of the x axis. It has to be a 
  pure function with the signature:
            
  `foo(x_data::Vector{Float64})::Vector{Float64}`
            
  And the output vector has to have the same length as the input one.
  Notice that this function will be applied regardless of units and possible domain problems (use 
  `x_range` to solve incompatibilities), and that it will not be reflected in the automatic label.
- `y_func::Function = identity`: Function to be applied to the values of the y axis. It has to be a 
  pure function with the signature:
              
  `foo(y_data::Vector{Float64})::Vector{Float64}`
            
  And the output vector has to have the same length as the input one.
  Notice that this function will be applied regardless of units and possible domain problems (use 
  `y_range` to solve incompatibilities), and that it will not be reflected in the automatic label.

"""
function snapshotTable(
    source_paths::Vector{String},
    base_names::Vector{String},
    quantities::Dict{Symbol,<:Union{String,Vector{String}}};
    human_redable::Bool = true,
    machine_redable::Bool = true,
    output_path::String = "snapshot_tables",
    sim_labels::Union{Vector{String},Nothing} = basename.(source_paths),
    title::String = "",
    sim_cosmo::Bool = false,
    warnings::Bool = true,
    show_progress::Bool = true,
    idx::IndexType = :,
    filter_functions::Vector{<:Union{Function,Nothing}} = [nothing,],
    data_analysis::Function = trivial,
    da_args::Tuple = (),
    da_kwargs::NamedTuple = (;),
    x_name::String = "X",
    y_name::String = "Y",
    xunit_label::String = "[:auto_label]",
    yunit_label::String = "[:auto_label]",
    x_unit::Unitful.Units = Unitful.NoUnits,
    y_unit::Unitful.Units = Unitful.NoUnits,
    x_factor::Int64 = 0,
    y_factor::Int64 = 0,
    x_range::NTuple{2,<:Real} = (-Inf, Inf),
    y_range::NTuple{2,<:Real} = (-Inf, Inf),
    x_edges::Bool = false,
    y_edges::Bool = false,
    x_func::Function = identity,
    y_func::Function = identity,
)::Nothing

    # If it doesn't exist create a folder to store the plots
    table_path = mkpath(output_path)

    length(base_names) == length(source_paths) == length(sim_labels) || throw(ArgumentError(
        "The arguments `base_names`, `source_paths` and `sim_labels` must have the same length."
    ))

    # Construct the labels for the axis
    unit_labels = [
        replace(xunit_label, ":auto_label" => getUnitLabel(x_factor, x_unit, latex = false)),
        replace(yunit_label, ":auto_label" => getUnitLabel(y_factor, y_unit, latex = false)),
    ]

    ################################################################################################
    # Load the simulation results
    ################################################################################################

    # Get the path and the global index of each snapshot, for every simulation
    source = getSnapshotPaths.(base_names, source_paths)

    # Make a dataframe where every row is a snapshot, and the columns are
    #  - Global index
    #  - Clock time
    #  - Number in the filename of the snapshot
    #  - Path to the snapshot (or `missing` is the snapshot doesn't exist), for the 1st simulation
    #  - Path to the snapshot (or `missing` is the snapshot doesn't exist), for the 2nd simulation
    #  - Path to the snapshot (or `missing` is the snapshot doesn't exist), for the 3rd simulation
    #  - ... 
    sim_df = makeSourceTable(source, idx, UnitfulAstro.Myr, sim_cosmo)

    # Each iteration is a snapshot
    iterator = enumerate(eachrow(sim_df))

    # Largest number of snapshots in a simulation
    max_n_snapshot = length(iterator)

    # Save the initial scale factor (only relevant for cosmological simulations), assuming that 
    # every simulation has the same value
    a0 = read_header(sim_df[1, 4]).time

    ################################################################################################
    # Main loop
    ################################################################################################

    # Initialize the progress bar
    prog_bar = Progress(
        max_n_snapshot,
        dt = 0.5,
        desc = "Analyzing and writting the data... ",
        color = :blue,
        barglyphs = BarGlyphs("|#  |"),
        enabled = show_progress,
    )

    for (local_index, (global_index, clock_time, number, snapshots...)) in iterator

        # Flag to skip problematic snapshots
        #   - true: The snapshot will be skipped.
        #   - false: The snapshot will be written.
        skipper = true

        # Initialize the vectors which will store the results
        data = Vector{Float64}[]
        rows_idx = Vector{Int64}[]
        labels = String[]

        for ((sim_index, snapshot), sim_label) in zip(enumerate(snapshots), sim_labels)

            # Skip missing snapshots
            snapshot !== missing || continue

            # Filter function for the current simulation in the loop
            filter_function = filter_functions[mod1(sim_index, length(filter_functions))]

            # Store basic information about the current snapshot and simulation in the loop
            metadata = Dict(
                :sim_data => SimData(
                    source_paths[sim_index],  # Path to the folder with all the snapshot files
                    sim_index,                # Index associated with the simulation
                    idx,                      # Which snapshots will be read
                    base_names[sim_index],    # Base name of the snapshot files
                    read_header(snapshot),    # Snapshot header of the simulation
                    filter_function,          # Filter function for the current simulation
                    a0,                       # Initial scale factor for the current simulation
                ),
                :snap_data => SnapData(
                    snapshot,                 # Full path
                    global_index,             # Global index 
                    local_index,              # Local index 
                    clock_time,               # Time stamp 
                ),
            )

            # Extract the relevant data from the snapshot
            data_dict = merge(
                getSnapshotData(snapshot, quantities; sim_cosmo, filter_function, warnings),
                metadata,
            )

            # Apply the analysis function
            da_output = data_analysis(data_dict, da_args...; da_kwargs...)

            # Skip this snapshot if `data_analysis` returns `nothing`
            da_output === nothing ? continue : skipper = false

            # Unit conversion
            axis_data = [
                ustrip.(unit, ax_data) for
                (unit, ax_data) in zip([x_unit, y_unit, Unitful.NoUnits], da_output)
            ]

            # Sanitize and transform the data to fit the user's requirements
            if length(axis_data) == 1

                sanitizeData!(
                    axis_data[1],
                    range = x_range,
                    keep_edges = x_edges,
                    exp_factor = x_factor,
                )

                axis_data[1] = x_func(axis_data[1])

            else

                sanitizeData!(
                    axis_data[1],
                    axis_data[2],
                    range = (x_range, y_range),
                    keep_edges = (x_edges, y_edges),
                    exp_factor = (x_factor, y_factor),
                )

                axis_data[1] = x_func(axis_data[1])
                axis_data[2] = y_func(axis_data[2])

            end

            # Save the results for this simulation
            push!(data, axis_data[1], axis_data[2])
            push!(rows_idx, [1:length(axis_data[1]);], [1:length(axis_data[2]);])
            push!(labels, "$x_name ($sim_label)", "$y_name ($sim_label)")

        end

        !skipper || (next!(prog_bar); continue)

        # Make a dataframe with all the results
        source_table = select!(
            unstack(
                flatten(DataFrame(; labels, data, rows_idx), [:data, :rows_idx]),
                :labels,
                :data,
            ),
            Not(:rows_idx),
        )

        # Save the output files
        output_file = joinpath(table_path, "snapshot_" * number)
        if human_redable
            open(output_file * ".txt", "w") do file
                pretty_table(
                    file,
                    source_table;
                    header = (names(source_table), unit_labels),
                    alignment = :c,
                    title
                )
            end
        end
        if machine_redable
            CSV.write(output_file * ".csv", source_table, missingstring = "missing")
        end

        # Advance the progress bar
        next!(prog_bar)

    end

    return nothing

end

"""
    timeSeriesTable(
        source_paths::Vector{String},
        base_names::Vector{String}; 
        <keyword arguments>
    )::Nothing

Generate a table with the time evolution of some quantity, for one or more simulations.

It transparently manages units, you only have to indicate the final unit of each axis.

The output can be a human-readable TXT or a CSV. In both cases the columns are:

- x quantity for simulation 1
- y quantity for simulation 1
- x quantity for simulation 2
- y quantity for simulation 2
- x quantity for simulation 3
- y quantity for simulation 3

etc.

# Arguments
- `source_paths::Vector{String}`: Paths to the folders containing the snapshots of each simulation, 
  set in the GADGET variable `OutputDir`.
- `base_names::Vector{String}`: Base names of the snapshot files, set in the GADGET variable 
  `SnapshotFileBase`.     
- `human_redable::Bool = true`: If the results will be written in a human readable format. The 
  output file in this case will be a `.txt` file.
- `machine_redable::Bool = true`: If the results will be written in a machine readable format. The 
  output file in this case will be a `.csv` file.
- `output_path::String = "time_series_table"`: Path to the output folder.
- `sim_labels::Vector{String}} = basename.(source_paths)`: Labels for the simulations.
- `file_name::String = "table"`: Name of the output files, without extension.
- `title::String = ""`: Title to be written at the top of the file, only aplicable to the `.txt`
  file.
            
## Data manipulation options
- `idx::IndexType = :`: Indexing of the simulation, i.e. indices of the snapshots that will be read.
  It can be an integer (a single snapshot), a vector of integers (several given snapshots), an
  UnitRange (e.g. 5:13) or : (every snapshot will be read).
- `filter_functions::Vector{<:Union{Function,Nothing}} = [nothing,]`: List of functions (one per 
  simulation), all with the signature: 
            
  `foo(file_path::String)::Vector{Int64}`
              
  Each function indicates which particles will be read in each simulation, taking the file path to 
  a snapshot and returning the list of indices of the selected particles. If set to `nothing`, then 
  every particle is read. See the [GadgetIO.jl](https://ludwigboess.github.io/GadgetIO.jl/stable/read_snapshots/#Filter-functions) documentation for examples.
- `data_analysis::Function = trivial`: Optional data analysis function. See the required signature 
  and examples in `src/data_analysis.jl`.
- `da_args::Tuple = ()`: Positional arguments for the `data_analysis` function.
- `da_kwargs::NamedTuple = (;)`: Keyword arguments for the `data_analysis` function.
            
## Axes options  
- `x_name::Union{String,LaTeXString} = "X"`: Name of the variable for the x axis. It should not 
  include units or scaling factors.
- `y_name::Union{String,LaTeXString} = "Y"`: Name of the variable for the y axis. It should not 
  include units or scaling factors.
- `xunit_label::String = "[:auto_label]"`: Label indicating the unit of the x axis. It can contain 
  the string `:auto_label`, which will be replaced by the automatic unit label that would otherwise 
  be used. 
- `yunit_label::String = "[:auto_label]"`: Label indicating the unit of the x axis. It can contain 
  the string `:auto_label`, which will be replaced by the automatic unit label that would otherwise 
  be used.
- `x_unit::Unitful.Units = Unitful.NoUnits`: Target unit for the x axis. The values will be 
  converted accordingly. Leave the default value of `Unitful.NoUnits` for dimensionless quantities.
- `y_unit::Unitful.Units = Unitful.NoUnits`: Target unit for the y axis. The values will be 
  converted accordingly. Leave the default value of `Unitful.NoUnits` for dimensionless quantities.
- `x_factor::Int64 = 0`: Numerical exponent to scale the x axis, e.g. if `x_factor` = 10 the values 
  will be scaled by ``10^{10}``. The default is no scaling.
- `y_factor::Int64 = 0`: Numerical exponent to scale the y axis, e.g. if `y_factor` = 10 the values 
  will be scaled by ``10^{10}``. The default is no scaling.
- `x_range::NTuple{2,<:Real} = (-Inf, Inf)`: The data will be trim down so its x coordinates fit 
  within `x_range`.
- `y_range::NTuple{2,<:Real} = (-Inf, Inf)`: The data will be trim down so its y coordinates fit 
  within `y_range`. Notice that this option does not affect histograms.
- `x_edges::Bool = false`: Set it to `true` if you want to keep the borders of `x_range`.
- `y_edges::Bool = false`: Set it to `true` if you want to keep the borders of `y_range`.
- `x_func::Function = identity`: Function to be applied to the values of the x axis. It has to be a 
  pure function with the signature:
            
  `foo(x_data::Vector{Float64})::Vector{Float64}`
            
  And the output vector has to have the same length as the input one.
  Notice that this function will be applied regardless of units and possible domain problems (use 
  `x_range` to solve incompatibilities), and that it will not be reflected in the automatic label.
- `y_func::Function = identity`: Function to be applied to the values of the y axis. It has to be a 
  pure function with the signature:
              
  `foo(y_data::Vector{Float64})::Vector{Float64}`
            
  And the output vector has to have the same length as the input one.
  Notice that this function will be applied regardless of units and possible domain problems (use 
  `y_range` to solve incompatibilities), and that it will not be reflected in the automatic label.

"""
function timeSeriesTable(
    source_paths::Vector{String},
    base_names::Vector{String};
    human_redable::Bool = true,
    machine_redable::Bool = true,
    output_path::String = "time_series_table",
    sim_labels::Vector{String} = basename.(source_paths),
    file_name::String = "table",
    title::String = "",
    idx::IndexType = :,
    filter_functions::Vector{<:Union{Function,Nothing}} = [nothing,],
    data_analysis::Function = trivial,
    da_args::Tuple = (),
    da_kwargs::NamedTuple = (;),
    x_name::String = "X",
    y_name::String = "Y",
    xunit_label::String = "[:auto_label]",
    yunit_label::String = "[:auto_label]",
    x_unit::Unitful.Units = Unitful.NoUnits,
    y_unit::Unitful.Units = Unitful.NoUnits,
    x_factor::Int64 = 0,
    y_factor::Int64 = 0,
    x_range::NTuple{2,<:Real} = (-Inf, Inf),
    y_range::NTuple{2,<:Real} = (-Inf, Inf),
    x_edges::Bool = false,
    y_edges::Bool = false,
    x_func::Function = identity,
    y_func::Function = identity
)::Nothing

    length(base_names) == length(source_paths) == length(sim_labels) || throw(ArgumentError(
        "The arguments `base_names`, `source_paths` and `sim_labels` must have the same length."
    ))

    # Construct the labels for the axis
    unit_labels = [
        "",
        replace(xunit_label, ":auto_label" => getUnitLabel(x_factor, x_unit, latex = false)),
        replace(yunit_label, ":auto_label" => getUnitLabel(y_factor, y_unit, latex = false)),
    ]

    # Each iteration is a simulation
    iterator = enumerate(zip(source_paths, base_names, sim_labels))

    # Initialize the vectors which will store the results
    data = Vector{Float64}[]
    snap_idx = Vector{Int64}[]
    labels = String[]

    for (sim_index, (source_path, base_name, sim_label)) in iterator

        # Filter function for the current simulation in the loop
        filter_function = filter_functions[mod1(sim_index, length(filter_functions))]

        # Header of the first snapshot
        header = read_header(getSnapshotPaths(base_name, source_path)["snap_paths"][1])

        # Store basic information about the current simulation in the loop
        sim_data = SimData(
            source_path,        # Path to the folder with all the snapshot files
            sim_index,          # Index associated with the simulation
            idx,                # Which snapshots will be read
            base_name,          # Base name of the snapshot files
            header,             # Snapshot header of the simulation
            filter_function,    # Filter function for the simulation
            header.time,        # Initial scale factor for the simulation
        )

        # Apply the analysis function
        da_output = data_analysis(sim_data, da_args...; da_kwargs...)

        # Skip this simulation if `data_analysis` returns `nothing`
        if da_output === nothing
            continue
        end

        # Unit conversion
        x_data, y_data = [
            ustrip.(unit, ax_data) for
            (unit, ax_data) in zip([x_unit, y_unit], da_output)
        ]

        # Sanitize and transform the data to fit the user's requirements
        sanitizeData!(
            x_data,
            y_data,
            range = (x_range, y_range),
            keep_edges = (x_edges, y_edges),
            exp_factor = (x_factor, y_factor),
        )
        x_data = x_func(x_data)
        y_data = y_func(y_data)

        # Save the results for this simulation
        push!(data, x_data, y_data)
        push!(snap_idx, [1:length(x_data);], [1:length(y_data);])
        push!(labels, "$x_name ($sim_label)", "$y_name ($sim_label)")

    end

    # Make a dataframe with all the results
    source_table = unstack(
        flatten(DataFrame(; labels, data, snap_idx), [:data, :snap_idx]),
        :labels,
        :data,
    )

    # If it doesn't exist create a folder to store the files
    table_path = joinpath(mkpath(output_path), file_name)

    # Save the output files
    if human_redable
        open(table_path * ".txt", "w") do file
            pretty_table(
                file,
                source_table;
                header = (names(source_table), unit_labels),
                alignment = :c,
                title
            )
        end
    end
    if machine_redable
        CSV.write(table_path * ".csv", source_table, missingstring = "missing")
    end

    return nothing

end





































# """
#     snapshotPlot(
#         source_paths::Vector{String},
#         base_names::Vector{String},
#         quantities::Dict{Symbol,<:Union{String,Vector{String}}},
#         plot_function!::Function; 
#         <keyword arguments>
#     )::Nothing

# Pipeline for generating one figure per snapshot, of one or more simulations.

# Some of the features are:

# - It can produce scatter plots, line plots, histograms, and density maps. The last two for 
#   only one simulation at a time.
# - It can produce a GIF or MP4 animating all the resulting figures.
# - It transparently manages units, you only have to indicate the target unit for the 
#   final figure.

# Many aspects of the plots can be configured, see the arguments below for details.

# # Arguments
# - `source_paths::Vector{String}`: Paths to the directories containing the snapshot files, 
#   set in the GADGET variable `OutputDir`.
# - `base_names::Vector{String}`: Base names of the snapshot files, set in the GADGET 
#   variable `SnapshotFileBase`.
# - `quantities::Dict{Symbol,<:Union{String,Vector{String}}}`: Dictionary where the keys are 
#   the particle types (the possibilities are given by [`ParticleType`](@ref) in `src/constants.jl`), 
#   and the values are the lists of data blocks to be extracted for each type.
# - `plot_function!::Function`: A plotting function from [Makie.jl](http://makie.juliaplots.org/stable/index.html). 
#   The options are:
#   * `hist!`         ⟶ Histograms.
#   * `scatter!`      ⟶ Scatter plots.
#   * `lines!`        ⟶ Line plots.
#   * `scatterlines!` ⟶ Scatter markers with lines between them.
#   * `heatmap!`      ⟶ Heatmaps. 
#   This will set the type of plot that will be generated.
# - `pf_kwargs::NamedTuple = (;)`: Keyword arguments for the `plot_function!` function.
# - `output_path::String = "snapshot_plots"`: Path to the output directory. The figures will 
#   be stored in `output_path`/plots/ and will be named snapshot\\_XXX`output_format` where 
#   XXX is the snapshot number.
# - `sim_cosmo::Bool = false`:  If the simulation is cosmological, 
#   * `false` ⟶ Newtonian simulation (`ComovingIntegrationOn` = 0).
#   * `true`  ⟶ Cosmological simulation (`ComovingIntegrationOn` = 1).
# - `filter_function::Union{Function,Nothing} = nothing`: A function with the signature: 

#   `foo(file_path::String)::Vector{Int64}`
  
#   If set to `nothing`, then no particles are filtered. See the [GadgetIO.jl](https://ludwigboess.github.io/GadgetIO.jl/stable/read_snapshots/#Filter-functions) documentation for 
#   examples.
# - `idx::IndexType`: Indexing of the simulation, i.e. which snapshots will be plotted.
#   It can be an integer (a single snapshot), a vector of integers (several given snapshots),
#   an UnitRange (e.g. 5:13) or : (every snapshot will be plotted).
# - `sim_labels::Union{Vector{String},Nothing} = basename.(source_paths)`: Labels for each 
#   simulation, for the plot legend. Set it to `nothing` if you don't want a legend.
# - `t_unit::Unitful.Units = UnitfulAstro.Myr`: Target unit for the time-stamps that will put 
#   as titles.
# - `warnings::Bool = true`: If a warning will be given when the data is not as expected,
#   but the function can still run using sane defaults.
# - `show_progress::Bool = true`: If a progress bar will be shown.
# - `data_analysis::Function = (x...) -> x`: Optional data analysis function.
#   See `src/data_analysis.jl` for the required signature and examples.
# - `da_args::Tuple = ()`: Positional arguments for the `data_analysis` function.
# - `da_kwargs::NamedTuple = (;)`: Keyword arguments for the `data_analysis` function.
# - `post_processing!::Function = (x...) -> nothing`: Optional post processing function.
#   See `src/post_processing.jl` for the required signature and examples.
# - `pp_args::Tuple = ()`: Positional arguments for the `post_processing!` function.
# - `pp_kwargs::NamedTuple = (;)`: Keyword arguments for the `post_processing!` function.
# - `x_label::::Union{String,LaTeXString}`: Name of the variable for the x axis label. 
#   It should not include units or scaling factors, these will be added automatically.
# - `y_label::::Union{String,LaTeXString}`: Name of the variable for the y axis label. 
#   It should not include units or scaling factors, these will be added automatically.
# - `x_unit::Unitful.Units = Unitful.NoUnits`: Target unit for the quantity in the x axis. All 
#   the x axis data will be converted to this unit, and it will be printed in the axis label.
# - `y_unit::Unitful.Units = Unitful.NoUnits`: Target unit for the quantity in the y axis. All 
#   the y axis data will be converted to this unit, and it will be printed in the axis label.
# - `z_unit::Unitful.Units = Unitful.NoUnits`: Target unit for the quantity in the z axis. 
#   Only valid for heatmaps, where the "z axis" is some kind of density. 
# - `x_factor::Int64 = 0`: Numerical exponent to scale the x axis, e.g. if `x_factor` = 10 
#   the corresponding axis will be scaled by ``10^{10}``. The default is no scaling.
# - `y_factor::Int64 = 0`: Numerical exponent to scale the y axis, e.g. if `y_factor` = 10 
#   the corresponding axis will be scaled by ``10^{10}``. The default is no scaling.
# - `x_scale::Union{Function,Makie.Symlog10} = identity`: Scaling to be used for the x axis. 
#   The options are the scaling functions accepted by [Makie.jl](http://makie.juliaplots.org/stable/index.html): log10, log2, log, sqrt, Makie.logit, Makie.Symlog10, Makie.pseudolog10 and identity.
# - `y_scale::Union{Function,Makie.Symlog10} = identity`: Scaling to be used for the y axis. 
#   The options are the scaling functions accepted by [Makie.jl](http://makie.juliaplots.org/stable/index.html): log10, log2, log, sqrt, Makie.logit, Makie.Symlog10, Makie.pseudolog10 and identity.
# - `x_limits::Tuple{<:Union{Real,Nothing},<:Union{Real,Nothing}}`: Set it to a value diferent 
#   than `nothing` if you want to fix the limits of the x axis.
# - `y_limits::Tuple{<:Union{Real,Nothing},<:Union{Real,Nothing}}`: Set it to a value diferent 
#   than `nothing` if you want to fix the limits of the y axis.
# - `x_range::NTuple{2,<:Real} = (-Inf, Inf)`: The data will be trim down so its x coordinates
#   fit within `x_range`.
# - `y_range::NTuple{2,<:Real} = (-Inf, Inf)`: The data will be trim down so its y coordinates
#   fit within `y_range`.
# - `x_edges::Bool = false`: If the borders of `x_range` are to be kept.
# - `y_edges::Bool = false`: If the borders of `y_range` are to be kept.
# - `animation::Bool = true`: If an animation using all the generated figures will be made.
#   The animation will be stored in `output_path`. 
# - `anim_file::String = "animation.gif"`: Filename of the animation, with its extension.
#   All formats supported by [Makie.jl](http://makie.juliaplots.org/stable/index.html) can be used, namely .mkv, .mp4, .webm and .gif.
# - `frame_rate::Int64 = 10`: Frame rate of the animation.
# - `output_format::String = ".svg"`: File format of the output figures. All formats supported by
#   [Makie.jl](http://makie.juliaplots.org/stable/index.html) can be used, namely ".pdf", ".svg" and ".png". 
# - `pt_per_unit::Float64 = 0.75`: Factor to scale up or down the size of the figures, keeping 
#   their proportions. This one only works for `.pdf` and `.svg`.
# - `px_per_unit::Float64 = 1.0`: Factor to scale up or down the size of the figures, keeping 
#   their proportions. This one only works for `.png`.
# - `resolution::NTuple{2,Int64} = (1000, 750)`: Resolution of the figures in points. 
#   For `.png`, by defult points = pixels (as given by `px_per_unit` = 1.0). For `.pdf` and
#   `.svg` points = 0.75 * pixels (as given by `pt_per_unit` = 0.75).
# - `aspect::Union{DataAspect,AxisAspect,Nothing} = nothing`: Aspect ratio of the figures.
#   The options are:
#   * `nothing`: uses the default given by [Makie.jl](http://makie.juliaplots.org/stable/index.html).
#   * `AxisAspect(n)`: The aspect ratio will be given by the number n = width / height.
#   * `DataAspect()`: The aspect ratio of the data is used.
# - `series_colors::Union{Vector{ColorInput},Nothing} = nothing`: Colors for the diferent 
#   series to be plotted (in general, diferent simulations). If set to `nothing`, the colors 
#   will be assigned automatically. This is only relevant for `scatter!` and `scatterlines!` 
#   plots.
# - `patchcolor::Union{ColorInput,Nothing} = nothing`: Color for the histogram bars.
#   If set to `nothing`, they will be red.
# - `series_markers::Union{Vector{Symbol},Nothing} = nothing`: Markers for the different series 
#   to be plotted (in general, different simulations). If set to `nothing`, the markers will be 
#   assigned automatically. This is only relevant for `scatter!` and `scatterlines!` plots.
# - `series_linetypes::Union{Vector{LineStyleInput},Nothing} = nothing`: Line types for the 
#   different series to be plotted (in general, different simulations). If set to `nothing`, 
#   the line types will be assigned automatically. This is only relevant for `lines!` and 
#   `scatterlines!` plots.

# """
# function snapshotPlot(
#     source_paths::Vector{String},
#     base_names::Vector{String},
#     quantities::Dict{Symbol,<:Union{String,Vector{String}}},
#     plot_function!::Function;
#     pf_kwargs::NamedTuple = (;),
#     output_path::String = "snapshot_plots",
#     sim_cosmo::Bool = false,
#     filter_function::Union{Function,Nothing} = nothing,
#     idx::IndexType = :,
#     sim_labels::Union{Vector{String},Nothing} = basename.(source_paths),
#     t_unit::Unitful.Units = UnitfulAstro.Myr,
#     warnings::Bool = true,
#     show_progress::Bool = true,
#     data_analysis::Function = (x...) -> x,
#     da_args::Tuple = (),
#     da_kwargs::NamedTuple = (;),
#     post_processing!::Function = (x...) -> nothing,
#     pp_args::Tuple = (),
#     pp_kwargs::NamedTuple = (;),
#     x_label::Union{String,LaTeXString} = "",
#     y_label::Union{String,LaTeXString} = "",
#     x_unit::Unitful.Units = Unitful.NoUnits,
#     y_unit::Unitful.Units = Unitful.NoUnits,
#     x_factor::Int64 = 0,
#     y_factor::Int64 = 0,
#     x_scale::Union{Function,Makie.Symlog10} = identity,
#     y_scale::Union{Function,Makie.Symlog10} = identity,
#     x_limits::Tuple{<:Union{Real,Nothing},<:Union{Real,Nothing}} = (nothing, nothing),
#     y_limits::Tuple{<:Union{Real,Nothing},<:Union{Real,Nothing}} = (nothing, nothing),
#     x_range::NTuple{2,<:Real} = (-Inf, Inf),
#     y_range::NTuple{2,<:Real} = (-Inf, Inf),
#     x_edges::Bool = false,
#     y_edges::Bool = false,
#     animation::Bool = true,
#     anim_file::String = "animation.gif",
#     frame_rate::Int64 = 10,
#     output_format::String = ".svg",
#     pt_per_unit::Float64 = 0.75,
#     px_per_unit::Float64 = 1.0,
#     resolution::NTuple{2,Int64} = (1000, 750),
#     aspect::Union{DataAspect,AxisAspect,Nothing} = nothing,
#     series_colors::Union{Vector{ColorInput},Nothing} = nothing,
#     patchcolor::Union{ColorInput,Nothing} = nothing,
#     series_markers::Union{Vector{Symbol},Nothing} = nothing,
#     series_linetypes::Union{Vector{LineStyleInput},Nothing} = nothing,
# )::Nothing

#     ########################################################################################
#     # Get the data
#     ########################################################################################

#     # Get the paths to the snapshots of each simulation, and their numbers
#     (
#         length(base_names) == length(source_paths) ||
#         throw(ArgumentError("`base_names` and `source_paths` must have the same length."))
#     )
#     source = getSnapshotPaths.(base_names, source_paths)

#     # Organize the input data in a dataframe
#     sim_df = makeSourceTable(source, idx, t_unit, sim_cosmo)

#     # Get the initial scale factor (only relevant for cosmological simulations)
#     # Assumes that every simulation has the same value
#     a0 = read_header(sim_df[1, 4]).time

#     # Pahts to the snapshots
#     snap_paths = @view sim_df[:, 4:end]
#     # Snapshot index independent of `idx`
#     global_index = @view sim_df[:, :id]
#     # Snapshot index after selection by `idx`
#     local_index = 1:nrow(sim_df)
#     # Clock time associated to each snapshot
#     time_stamps = @view sim_df[:, :time_stamps]
#     # Number in the filename of each snapshot
#     numbers = @view sim_df[:, :numbers]

#     # Every iteration is a snapshot
#     iterator = zip(eachrow(snap_paths), global_index, local_index, time_stamps, numbers)

#     # Create a directory to save the output files if it doesn't exist
#     output_folder = mkpath(joinpath(output_path, "plots"))

#     ########################################################################################
#     # Set up the canvas and the animation
#     ########################################################################################

#     # Create a list of distinguishable colors for the different simulations
#     colors = distinguishable_colors(
#         length(base_names),
#         [RGB(1, 1, 1), RGB(0, 0, 0)],
#         dropseed = true,
#     )
#     # Choose different markers for the different simulations 
#     # (used in scatter and scatterlines plots)
#     markers = MARKERS[1:length(base_names)]
#     # Choose different line styles for the different simulations 
#     # (used in lines and scatterlines plots)
#     linestyles = LINE_STYLES[1:length(base_names)]

#     # Apply global theme
#     set_theme!(theme)
#     # Set colors, patchcolors, markers and line styles
#     update_theme!(
#         Axis = (
#             palette = (
#                 color = series_colors !== nothing ? series_colors : colors,
#                 patchcolor = patchcolor !== nothing ? patchcolor : [:red],
#                 marker = series_markers !== nothing ? series_markers : markers,
#                 linestyle = series_linetypes !== nothing ? series_linetypes : linestyles,
#             ),
#         ),
#     )

#     # Create an empty figure
#     figure = Figure(; resolution)

#     # Set the axes
#     axes = Axis(
#         figure[1, 1];
#         xlabel = getLabel(x_label, x_factor, x_unit),
#         ylabel = getLabel(y_label, y_factor, y_unit),
#         limits = (x_limits, y_limits),
#         aspect,
#     )

#     (
#         !(nrow(sim_df) < frame_rate && warnings) || @warn(
#             "With the settings: `frame_rate` = $frame_rate and `idx` = $idx, the animation is less than a second long."
#         )
#     )

#     # Initialize the animation stream
#     animation ? vs = VideoStream(figure, framerate = frame_rate) : nothing

#     ########################################################################################
#     # Main loop
#     ########################################################################################

#     # Initialize the progress bar
#     prog_bar = Progress(
#         length(iterator),
#         dt = 0.5,
#         desc = "Analyzing and plotting the data... ",
#         color = :blue,
#         barglyphs = BarGlyphs("|#  |"),
#         enabled = show_progress,
#     )

#     for (snapshot, global_idx, local_idx, t, number) in iterator

#         # Flag to keep the x axis with a linear scale if there are less than 3 data points
#         # * true: The axis will have a linear scale
#         # * false: The axis will have the scale given by `x_scale`
#         xscale_flag = false
#         # Flag to keep the y axis with a linear scale if there are less than 3 data points
#         # * true: The axis will have a linear scale
#         # * false: The axis will have the scale given by `y_scale`
#         yscale_flag = false
#         # Flag to skip problematic snapshots
#         # * true: The snapshot will be skipped.
#         # * false: The snapshot will be plotted.
#         skipper = true

#         for (sim_idx, snap) in enumerate(snapshot)

#             # Skip missing snapshots
#             snap !== missing || continue

#             # Save basic information about this snapshot
#             snap_data = Dict(
#                 :snap_data => SnapData(
#                     snap,                # Path to the snapshot
#                     read_header(snap),   # Header of the snapshot
#                     global_idx,          # Global index of the snapshot
#                     local_idx,           # Local index of the snapshot 
#                     sim_idx,             # Simulation index
#                     t,                   # Time stamp of the snapshot
#                     a0,                  # Initial scale factor for the simulation
#                 ),
#             )

#             # Get the data fields to be analyzed
#             data = merge(
#                 getSnapshotData(snap, quantities; sim_cosmo, filter_function, warnings),
#                 snap_data,
#             )

#             # Apply the analysis function to the data
#             da_output = data_analysis(data, da_args...; da_kwargs...)

#             # Skip snapshot if `data_analysis` returns `nothing`
#             da_output === nothing ? continue : skipper = false
#             # Data size check
#             checkDataShape(plot_function!, length(da_output), length(source_paths))

#             # Unit conversion
#             axis_data = [
#                 ustrip.(unit, ax_data) for 
#                 (unit, ax_data) in zip([x_unit, y_unit], da_output)
#             ]
          
#             # Sanitize data
#             if length(axis_data) == 1
#                 flag_x, _ = sanitizeData!(
#                     axis_data[1],
#                     scale = x_factor,
#                     func_domain = x_scale,
#                     range = x_range,
#                     keep_edges = x_edges,
#                     min_left = 3,
#                 )
#                 flag_y = 0
#             else
#                 flag_x, flag_y = sanitizeData!(
#                     axis_data[1],
#                     axis_data[2],
#                     scale = (x_factor, y_factor),
#                     func_domain = (x_scale, y_scale),
#                     range = (x_range, y_range),
#                     keep_edges = (x_edges, y_edges),
#                     min_left = 3,
#                 )
#             end
#             # If, in this snapshot and for any simulation, filtering the data targeting a 
#             # non linear scale would leaves less that three data points, the scale will 
#             # revert to `identity` and the data will be left untouched
#             if flag_x == 1
#                 xscale_flag = true
#             end
#             if flag_y == 1
#                 yscale_flag = true
#             end

#             # Plot the data
#             plot_function!(
#                 axes,
#                 axis_data...;
#                 bins = scaledBins(axis_data[1], 20, x_scale; limits = x_limits, warnings),
#                 colormap = :inferno,
#                 pf_kwargs...,
#             )

#         end

#         !skipper || (next!(prog_bar); continue)

#         # Set the scale of the axes
#         axes.xscale = xscale_flag ? identity : x_scale
#         axes.yscale = yscale_flag ? identity : y_scale

#         # Set the time-stamp as the title
#         time_stamp = string(round(typeof(t), t, sigdigits = 4))
#         axes.title = "t = $time_stamp"

#         # Set a legend if required
#         if sim_labels !== nothing
#             if plot_function! == lines!
#                 legend_element = [
#                     LineElement(; color, linestyle) for
#                     (color, linestyle) in zip(colors, linestyles)
#                 ]
#             elseif plot_function! == scatter!
#                 legend_element = [
#                     MarkerElement(; color, marker) for
#                     (color, marker) in zip(colors, markers)
#                 ]
#             elseif plot_function! == scatterlines!
#                 legend_element = [
#                     [MarkerElement(; color, marker), LineElement(; color, linestyle)] for 
#                     (color, marker, linestyle) in zip(colors, markers, linestyles)
#                 ]
#             end
#             Legend(figure[2, 1][1, 1], legend_element, sim_labels)
#         end

#         # Apply the post processing function
#         pp_legend = post_processing!(figure, pp_args...; pp_kwargs...)

#         # Set a post processing legend if necessary
#         if pp_legend !== nothing
#             Legend(
#                 figure[2, 1][1, sim_labels !== nothing ? 2 : 1],
#                 pp_legend[1],
#                 pp_legend[2],
#             )
#         end

#         # Save the figure
#         save(
#             joinpath(output_folder, "snapshot_" * number * output_format),
#             figure;
#             pt_per_unit,
#             px_per_unit,
#         )

#         # Add the figure as a frame to the animation
#         animation ? recordframe!(vs) : nothing

#         # Clean the figure for the next step in the loop
#         cleanPlot!(figure)

#         # Advance the progress bar
#         next!(prog_bar)

#     end

#     # Save the animation
#     animation ? save(joinpath(output_path, anim_file), vs) : nothing

#     return nothing

# end

# """
#     timeSeriesPlot(
#         source_paths::Vector{String},
#         base_names::Vector{String},
#         plot_function!::Function; 
#         <keyword arguments>
#     )::Nothing

# Pipeline for generating a figure with the time evolution of some quantity, for one or more 
# simulations.

# Some of the features are:

# - It can produce a scatter plot, a line plot, or a combination of both.
# - It manages units in a transparent manner, you only have to indicate the target unit for the 
#   final figure.

# Many aspect of the figure can be configured, see the arguments below for details.

# # Arguments
# - `source_paths::Vector{String}`: Paths to the directories containing the snapshot files, 
#   set in the GADGET variable `OutputDir`.
# - `base_names::Vector{String}`: Base names of the snapshot files, set in the GADGET 
#   variable `SnapshotFileBase`.
# - `plot_function!::Function`: A plotting function from [Makie.jl](http://makie.juliaplots.org/stable/index.html). 
#   The options are:
#   * `scatter!`      ⟶ Scatter plots.
#   * `lines!`        ⟶ Line plots.
#   * `scatterlines!` ⟶ Scatter markers with lines between them.
#   This will set the type of plot that will be generated.
# - `pf_kwargs::NamedTuple = (;)`: Keyword arguments for the `plot_function!` function.
# - `output_path::String = "snapshot_plots"`: Path to the output directory, where the final 
#   figure will be stored.
# - `filter_function::Union{Function,Nothing} = nothing`: A function with the signature: 

#   `foo(file_path::String)::Vector{Int64}`
  
#   If set to `nothing`, then no particles are filtered. See the [GadgetIO.jl](https://ludwigboess.github.io/GadgetIO.jl/stable/read_snapshots/#Filter-functions) documentation for 
#   examples.
# - `idx::IndexType`: Indexing of the simulation, i.e. which snapshots will be plotted.
#   It can be an integer (a single snapshot), a vector of integers (several given snapshots),
#   an UnitRange (e.g. 5:13) or : (every snapshot will be plotted).
# - `sim_labels::Union{Vector{String},Nothing} = basename.(source_paths)`: Labels for each 
#   simulation, for the plot legend. Set it to `nothing` if you don't want a legend.
# - `title::String = ""`: Title for the figure.
# - `data_analysis::Function = (x...) -> x`: Optional data analysis function.
#   See `src/data_analysis.jl` for the required signature and examples.
# - `da_args::Tuple = ()`: Positional arguments for the `data_analysis` function.
# - `da_kwargs::NamedTuple = (;)`: Keyword arguments for the `data_analysis` function.
# - `post_processing!::Function = (x...) -> nothing`: Optional post processing function.
#   See `src/post_processing.jl` for the required signature and examples.
# - `pp_args::Tuple = ()`: Positional arguments for the `post_processing!` function.
# - `pp_kwargs::NamedTuple = (;)`: Keyword arguments for the `post_processing!` function.
# - `x_label::::Union{String,LaTeXString}`: Name of the variable for the x axis label. 
#   It should not include units or scaling factors, these will be added automatically.
# - `y_label::::Union{String,LaTeXString}`: Name of the variable for the y axis label. 
#   It should not include units or scaling factors, these will be added automatically.
# - `x_unit::Unitful.Units = Unitful.NoUnits`: Target unit for the quantity in the x axis. All 
#   the x axis data will be converted to this unit, and it will be printed in the axis label.
# - `y_unit::Unitful.Units = Unitful.NoUnits`: Target unit for the quantity in the y axis. All 
#   the y axis data will be converted to this unit, and it will be printed in the axis label.
# - `x_factor::Int64 = 0`: Numerical exponent to scale the x axis, e.g. if `x_factor` = 10 
#   the corresponding axis will be scaled by ``10^{10}``. The default is no scaling.
# - `y_factor::Int64 = 0`: Numerical exponent to scale the y axis, e.g. if `y_factor` = 10 
#   the corresponding axis will be scaled by ``10^{10}``. The default is no scaling.
# - `x_scale::Union{Function,Makie.Symlog10} = identity`: Scaling to be used for the x axis. 
#   The options are the scaling functions accepted by [Makie.jl](http://makie.juliaplots.org/stable/index.html): log10, log2, log, sqrt, Makie.logit, Makie.Symlog10, Makie.pseudolog10 and identity.
# - `y_scale::Union{Function,Makie.Symlog10} = identity`: Scaling to be used for the y axis. 
#   The options are the scaling functions accepted by [Makie.jl](http://makie.juliaplots.org/stable/index.html): log10, log2, log, sqrt, Makie.logit, Makie.Symlog10, Makie.pseudolog10 and identity.
# - `x_limits::Tuple{<:Union{Real,Nothing},<:Union{Real,Nothing}} = (nothing, nothing)`: Set it to a 
#   value diferent than `nothing` if you want to fix the limits of the x axis.
# - `y_limits::Tuple{<:Union{Real,Nothing},<:Union{Real,Nothing}} = (nothing, nothing)`: Set it to a 
#   value diferent than `nothing` if you want to fix the limits of the y axis.
# - `x_range::NTuple{2,<:Real} = (-Inf, Inf)`: The data will be trim down so its x coordinates
#   fit within `x_range`.
# - `y_range::NTuple{2,<:Real} = (-Inf, Inf)`: The data will be trim down so its y coordinates
#   fit within `y_range`.
# - `x_edges::Bool = false`: If the borders of `x_range` are to be kept.
# - `y_edges::Bool = false`: If the borders of `y_range` are to be kept.
# - `file_name::String = "figure.svg"`: Filename for the output figure with its extension. All 
#   formats supported by [Makie.jl](http://makie.juliaplots.org/stable/index.html) can be used, namely ".pdf", ".svg" and ".png". 
# - `pt_per_unit::Float64 = 0.75`: Factor to scale up or down the size of the figures, keeping 
#   their proportions. This one only works for `.pdf` and `.svg`.
# - `px_per_unit::Float64 = 1.0`: Factor to scale up or down the size of the figures, keeping 
#   their proportions. This one only works for `.png`.
# - `resolution::NTuple{2,Int64} = (1000, 750)`: Resolution of the figures in points. 
#   For `.png`, by default points = pixels (as given by `px_per_unit` = 1.0). For `.pdf` and
#   `.svg` points = 0.75 * pixels (as given by `pt_per_unit` = 0.75).
# - `aspect::Union{DataAspect,AxisAspect,Nothing} = nothing`: Aspect ratio of the figures.
#   The options are:
#   * `nothing`: uses the default given by [Makie.jl](http://makie.juliaplots.org/stable/index.html).
#   * `AxisAspect(n)`: The aspect ratio will be given by the number n = width / height.
#   * `DataAspect()`: The aspect ratio of the data is used.
# - `series_colors::Union{Vector{ColorInput},Nothing} = nothing`: Colors for the different 
#   series to be plotted (in general, different simulations). If set to `nothing`, the colors 
#   will be assigned automatically. This is only relevant for `scatter!` and `scatterlines!` 
#   plots.
# - `series_markers::Union{Vector{Symbol},Nothing} = nothing`: Markers for the different series 
#   to be plotted (in general, different simulations). If set to `nothing`, the markers will be 
#   assigned automatically. This is only relevant for `scatter!` and `scatterlines!` plots.
# - `series_linetypes::Union{Vector{LineStyleInput},Nothing} = nothing`: Line types for the 
#   different series to be plotted (in general, different simulations). If set to `nothing`, 
#   the line types will be assigned automatically. This is only relevant for `lines!` and 
#   `scatterlines!` plots.

# """
# function timeSeriesPlot(
#     source_paths::Vector{String},
#     base_names::Vector{String},
#     plot_function!::Function;
#     pf_kwargs::NamedTuple = (;),
#     output_path::String = "time_series_plots",
#     filter_function::Union{Function,Nothing} = nothing,
#     idx::IndexType = :,
#     sim_labels::Union{Vector{String},Nothing} = basename.(source_paths),
#     title::String = "",
#     data_analysis::Function = (x...) -> x,
#     da_args::Tuple = (),
#     da_kwargs::NamedTuple = (;),
#     post_processing!::Function = (x...) -> nothing,
#     pp_args::Tuple = (),
#     pp_kwargs::NamedTuple = (;),
#     x_label::Union{String,LaTeXString} = "",
#     y_label::Union{String,LaTeXString} = "",
#     x_unit::Unitful.Units = Unitful.NoUnits,
#     y_unit::Unitful.Units = Unitful.NoUnits,
#     x_factor::Int64 = 0,
#     y_factor::Int64 = 0,
#     x_scale::Union{Function,Makie.Symlog10} = identity,
#     y_scale::Union{Function,Makie.Symlog10} = identity,
#     x_limits::Tuple{<:Union{Real,Nothing},<:Union{Real,Nothing}} = (nothing, nothing),
#     y_limits::Tuple{<:Union{Real,Nothing},<:Union{Real,Nothing}} = (nothing, nothing),
#     x_range::NTuple{2,<:Real} = (-Inf, Inf),
#     y_range::NTuple{2,<:Real} = (-Inf, Inf),
#     x_edges::Bool = false,
#     y_edges::Bool = false,
#     file_name::String = "figure.svg",
#     pt_per_unit::Float64 = 0.75,
#     px_per_unit::Float64 = 1.0,
#     resolution::NTuple{2,Int64} = (1000, 750),
#     aspect::Union{DataAspect,AxisAspect,Nothing} = nothing,
#     series_colors::Union{Vector{ColorInput},Nothing} = nothing,
#     series_markers::Union{Vector{Symbol},Nothing} = nothing,
#     series_linetypes::Union{Vector{LineStyleInput},Nothing} = nothing,
# )::Nothing

#     ########################################################################################
#     # Set up the canvas
#     ########################################################################################    

#     # Create a list of distinguishable colors for the different simulations
#     colors = distinguishable_colors(
#         length(base_names),
#         [RGB(1, 1, 1), RGB(0, 0, 0)],
#         dropseed = true,
#     )
#     # Choose different markers for the different simulations 
#     # (used in scatter and scatterlines plots)
#     markers = MARKERS[1:length(base_names)]
#     # Choose different line styles for the different simulations 
#     # (used in lines and scatterlines plots)
#     linestyles = LINE_STYLES[1:length(base_names)]

#     # Apply global theme
#     set_theme!(theme)
#     # Set colors, markers and line styles
#     update_theme!(
#         Axis = (
#             palette = (
#                 color = series_colors !== nothing ? series_colors : colors,
#                 marker = series_markers !== nothing ? series_markers : markers,
#                 linestyle = series_linetypes !== nothing ? series_linetypes : linestyles,
#             ),
#         ),
#     )

#     # Create an empty figure
#     figure = Figure(; resolution)

#     # Set the axes
#     axes = Axis(
#         figure[1, 1];
#         xlabel = getLabel(x_label, x_factor, x_unit),
#         ylabel = getLabel(y_label, y_factor, y_unit),
#         limits = (x_limits, y_limits),
#         title,
#         aspect,
#     )

#     ########################################################################################
#     # Main loop
#     ########################################################################################

#     # Flag to keep the x axis with a linear scale if there are less than 3 data points
#     # * true: The axis will have a linear scale
#     # * false: The axis will have the scale given by `x_scale`
#     xscale_flag = false
#     # Flag to keep the y axis with a linear scale if there are less than 3 data points
#     # * true: The axis will have a linear scale
#     # * false: The axis will have the scale given by `y_scale`
#     yscale_flag = false

#     for (source_path, base_name) in zip(source_paths, base_names)

#         # Save basic information about this simulation
#         sim_data = SimData(source_path, base_name, idx, filter_function)

#         # Apply the analysis function to the data
#         da_output = data_analysis(sim_data, da_args...; da_kwargs...)

#         # Unit conversion
#         axis_data = [ustrip.(unit, ax_data) for (unit, ax_data) in zip([x_unit, y_unit], da_output)]

#         # Sanitize data
#         flag_x, flag_y = sanitizeData!(
#             axis_data[1],
#             axis_data[2],
#             scale = (x_factor, y_factor),
#             func_domain = (x_scale, y_scale),
#             range = (x_range, y_range),
#             keep_edges = (x_edges, y_edges),
#             min_left = 3,
#         )
#         # If, for this simulation, filtering the data targeting a non linear scale would 
#         # leaves less that three data points, the scale will revert to `identity` and the
#         # data will be left untouched
#         if flag_x == 1
#             xscale_flag = true
#         end
#         if flag_y == 1
#             yscale_flag = true
#         end

#         # Plot the data
#         plot_function!(axes, axis_data...; pf_kwargs...)

#     end

#     # Set the scale of the axes
#     axes.xscale = xscale_flag ? identity : x_scale
#     axes.yscale = yscale_flag ? identity : y_scale

#     # Set a legend if required
#     if sim_labels !== nothing
#         if plot_function! == lines!
#             legend_element = [
#                 LineElement(; color, linestyle) for
#                 (color, linestyle) in zip(colors, linestyles)
#             ]
#         elseif plot_function! == scatter!
#             legend_element =
#                 [MarkerElement(; color, marker) for (color, marker) in zip(colors, markers)]
#         elseif plot_function! == scatterlines!
#             legend_element = [
#                 [MarkerElement(; color, marker), LineElement(; color, linestyle)] for
#                 (color, marker, linestyle) in zip(colors, markers, linestyles)
#             ]
#         end
#         Legend(figure[2, 1][1, 1], legend_element, sim_labels)
#     end

#     # Apply the post processing function
#     pp_legend = post_processing!(figure, pp_args...; pp_kwargs...)

#     # Set a post processing legend if necessary
#     if pp_legend !== nothing
#         Legend(
#             figure[2, 1][1, sim_labels !== nothing ? 2 : 1],
#             pp_legend[1],
#             pp_legend[2],
#         )
#     end

#     # Create a directory to save the output file if it doesn't exist
#     fig_path = mkpath(output_path)

#     # Save the figure
#     save(joinpath(fig_path, file_name), figure; pt_per_unit, px_per_unit)

#     return nothing

# end
