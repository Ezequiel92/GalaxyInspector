####################################################################################################
# 1D histogram utilities
####################################################################################################

"""
    listHistogram1D(
        positions::Vector{<:Number},
        values::Vector{<:Number},
        edges::Vector{<:Number},
    )::Vector{Vector{<:Number}}

Compute a 1D histogram of `values`, returning the full list of `values` within each bin.

# Arguments

  - `positions::Vector{<:Number}`: The positions of `values` within a 1D axis.
  - `values::Vector{<:Number}`: The values that will be sorted, according to their `positions`.
  - `edges::Vector{<:Number}`: A list of bin edges.

# Returns

    - A vector with the lists of `values` within each bin.
"""
function listHistogram1D(
    positions::Vector{<:Number},
    values::Vector{<:Number},
    edges::Vector{<:Number},
)::Vector{Vector{<:Number}}

    (
        length(values) == length(positions) ||
        throw(ArgumentError("listHistogram1D: `values` must have as many elements as \
        `positions`, but I got length(`values`) = $(length(values)) and \
        length(`positions`) = $(length(positions))"))
    )

    issorted(edges) || sort!(edges)

    n_bins = length(edges) - 1

    p_min = first(edges)
    p_max = last(edges)

    histogram = [eltype(values)[] for _ in 1:n_bins]

    # Compute the histogram, ignoring NaNs and positions outside the grid
    for (position, value) in zip(positions, values)

        if isnan(position) || isnan(value)
            continue
        elseif position < p_min || p_max < position
            continue
        elseif position == p_min
            idx = 1
        else
            idx = searchsortedfirst(edges, position) - 1
        end

        push!(histogram[idx], value)

    end

    return histogram

end

"""
    listHistogram1D(
        positions::Vector{<:Number},
        values::Vector{<:Number},
        grid::LinearGrid,
    )::Vector{Vector{<:Number}}

Compute a 1D histogram of `values`, returning the full list of `values` within each bin.

# Arguments

  - `positions::Vector{<:Number}`: The positions of `values` within a 1D axis.
  - `values::Vector{<:Number}`: The values that will be sorted, according to their `positions`.
  - `grid::LinearGrid`: A linear grid.

# Returns

    - A vector with the lists of `values` within each bin.
"""
function listHistogram1D(
    positions::Vector{<:Number},
    values::Vector{<:Number},
    grid::LinearGrid,
)::Vector{Vector{<:Number}}

    return listHistogram1D(positions, values, grid.x_edges)

end

"""
    histogram1D(
        positions::Vector{<:Number},
        values::Vector{<:Number},
        edges::Vector{<:Number};
        <keyword arguments>
    )::Vector{<:Number}

Compute a 1D histogram of `values`.

# Arguments

  - `positions::Vector{<:Number}`: Positions of the `values` within a 1D axis.
  - `values::Vector{<:Number}`: The values that will be added up in each bin, according to their `positions`.
  - `edges::Vector{<:Number}`: A list of bin edges.
  - `total::Bool=true`: If the sum (`total` = true) or the mean (`total` = false) of `values` will be computed for each bin.
  - `empty_nan::Bool=true`: If NaN will be put into empty bins, 0 is used otherwise.

# Returns

  - A vector with the histogram values.
"""
function histogram1D(
    positions::Vector{<:Number},
    values::Vector{<:Number},
    edges::Vector{<:Number};
    total::Bool=true,
    empty_nan::Bool=true,
)::Vector{<:Number}

    (
        length(values) == length(positions) ||
        throw(ArgumentError("histogram1D: `values` must have as many elements as \
        `positions`, but I got length(`values`) = $(length(values)) and \
        length(`positions`) = $(length(positions))"))
    )

    issorted(edges) || sort!(edges)

    n_bins = length(edges) - 1

    p_min = first(edges)
    p_max = last(edges)

    histogram = zeros(typeof(first(values)), n_bins)
    counts = zeros(Int, n_bins)

    # Compute the histogram, ignoring NaNs and positions outside the range
    for (position, value) in zip(positions, values)

        if isnan(position) || isnan(value)
            continue
        elseif position < p_min || p_max < position
            continue
        elseif position == p_min
            idx = 1
        else
            idx = searchsortedfirst(edges, position) - 1
        end

        histogram[idx] += value
        counts[idx] += 1

    end

    if empty_nan
        # Set empty bins to NaN
        nan = NaN * unit(first(values))
        Threads.@threads for i in eachindex(counts)
            if iszero(counts[i])
                histogram[i] = nan
            end
        end
    end

    if !total
        # Compute the mean value instead of just the sum for each bin
        Threads.@threads for i in eachindex(counts)
            if !iszero(counts[i])
                histogram[i] /= counts[i]
            end
        end
    end

    return histogram

end

"""
    histogram1D(
        positions::Vector{<:Number},
        values::Vector{<:Number},
        grid::LinearGrid;
        <keyword arguments>
    )::Vector{<:Number}

Compute a 1D histogram of `values`.

# Arguments

  - `positions::Vector{<:Number}`: Positions of the `values` within a 1D axis.
  - `values::Vector{<:Number}`: The values that will be added up to each bin, according to their `positions`.
  - `grid::LinearGrid`: A linear grid.
  - `total::Bool=true`: If the sum (`total` = true) or the mean (`total` = false) of `values` will be computed for each bin.
  - `empty_nan::Bool=true`: If NaN will be put into empty bins, 0 is used otherwise.

# Returns

  - A vector with the histogram values.
"""
function histogram1D(
    positions::Vector{<:Number},
    values::Vector{<:Number},
    grid::LinearGrid;
    total::Bool=true,
    empty_nan::Bool=true,
)::Vector{<:Number}

    return histogram1D(positions, values, grid.x_edges; total, empty_nan)

end

"""
    histogram1D(
        positions::Vector{<:Number},
        edges::Vector{<:Number};
        <keyword arguments>
    )::Vector{Float64}

Compute a 1D histogram of `positions`.

# Arguments

  - `positions::Vector{<:Number}`: Values for which the histogram will be constructed.
  - `edges::Vector{<:Number}`: A list of bin edges.
  - `empty_nan::Bool=true`: If NaN will be put into empty bins, 0 is used otherwise.

# Returns

  - A vector with the counts.
"""
function histogram1D(
    positions::Vector{<:Number},
    edges::Vector{<:Number};
    empty_nan::Bool=true,
)::Vector{Float64}

    issorted(edges) || sort!(edges)

    n_bins = length(edges) - 1

    p_min = first(edges)
    p_max = last(edges)

    histogram = zeros(Float64, n_bins)

    # Compute the histogram, ignoring NaNs and positions outside the range
    for position in positions

        if isnan(position)
            continue
        elseif position < p_min || p_max < position
            continue
        elseif position == p_min
            idx = 1
        else
            idx = searchsortedfirst(edges, position) - 1
        end

        histogram[idx] += 1.0

    end

    empty_nan && replace!(x -> iszero(x) ? NaN : x, histogram)

    return histogram

end

"""
    histogram1D(
        positions::Vector{<:Number},
        grid::LinearGrid;
        <keyword arguments>
    )::Vector{Float64}

Compute a 1D histogram of `positions`.

# Arguments

  - `positions::Vector{<:Number}`: Values for which the histogram will be constructed.
  - `grid::LinearGrid`: A linear grid.
  - `empty_nan::Bool=true`: If NaN will be put into empty bins, 0 is used otherwise.

# Returns

  - A vector with the counts.
"""
function histogram1D(
    positions::Vector{<:Number},
    grid::LinearGrid;
    empty_nan::Bool=true,
)::Vector{Float64}

    return histogram1D(positions, grid.x_edges; empty_nan)

end

####################################################################################################
# 2D histogram utilities
####################################################################################################

"""
    histogram2D(
        positions::Matrix{<:Number},
        values::Vector{<:Number},
        x_edges::Vector{<:Number},
        y_edges::Vector{<:Number};
        <keyword arguments>
    )::Matrix{<:Number}

Compute a 2D histogram of `values`.

# Arguments

  - `positions::Matrix{<:Number}`: Positions of the values in the grid. Each column corresponds to a value and each row is a dimension.
  - `values::Vector{<:Number}`: The values that will be added up in each square bin, according to their `positions`.
  - `x_edges::Vector{<:Number}`: A list of bin edges for the x axis.
  - `y_edges::Vector{<:Number}`: A list of bin edges for the y axis.
  - `total::Bool=true`: If the sum (`total` = true) or the mean (`total` = false) of `values` will be computed for each bin.
  - `empty_nan::Bool=true`: If NaN will be put into empty bins, 0 is used otherwise.

# Returns

  - A matrix with the histogram values.
"""
function histogram2D(
    positions::Matrix{<:Number},
    values::Vector{<:Number},
    x_edges::Vector{<:Number},
    y_edges::Vector{<:Number};
    total::Bool=true,
    empty_nan::Bool=true,
)::Matrix{<:Number}

    (
        length(values) == size(positions, 2) ||
        throw(ArgumentError("histogram2D: `values` must have as many elements as `positions` \
        has columns, but I got length(`values`) = $(length(values)) and size(`positions`, 2) = \
        $(size(positions, 2))"))
    )

    issorted(x_edges) || sort!(x_edges)
    issorted(y_edges) || sort!(y_edges)

    n_x_bins = length(x_edges) - 1
    n_y_bins = length(y_edges) - 1

    x_borders = (first(x_edges), last(x_edges))
    y_borders = (first(y_edges), last(y_edges))

    histogram = zeros(typeof(first(values)), (n_x_bins, n_y_bins))
    counts = zeros(Int, (n_x_bins, n_y_bins))

    # Compute the histogram, ignoring NaNs and positions outside the range
    for (i, point) in pairs(eachcol(positions))

        !isnan(values[i]) || continue
        !any(isnan, point) || continue

        x = point[1]
        y = point[2]

        !(x > x_borders[2] || x < x_borders[1]) || continue
        !(y > y_borders[2] || y < y_borders[1]) || continue

        if x == x_borders[1]
            i_x = 1
        else
            i_x = searchsortedfirst(x_edges, x) - 1
        end

        if y == y_borders[1]
            i_y = 1
        else
            i_y = searchsortedfirst(y_edges, y) - 1
        end

        histogram[n_y_bins - i_y + 1, i_x] += values[i]
        counts[n_y_bins - i_y + 1, i_x] += 1

    end

    if empty_nan
        # Set empty bins to NaN
        nan = NaN * unit(first(values))
        Threads.@threads for i in eachindex(histogram)
            if iszero(counts[i])
                histogram[i] = nan
            end
        end
    end

    if !total
        # Compute the mean value instead of just the sum for each bin
        Threads.@threads for i in eachindex(histogram)
            if !iszero(counts[i])
                histogram[i] /= counts[i]
            end
        end
    end

    return histogram

end

"""
    histogram2D(
        positions::Matrix{<:Number},
        values::Vector{<:Number},
        grid::SquareGrid;
        <keyword arguments>
    )::Matrix{<:Number}

Compute a 2D histogram of `values`.

# Arguments

  - `positions::Matrix{<:Number}`: Positions of the values in the grid. Each column corresponds to a value and each row is a dimension.
  - `values::Vector{<:Number}`: The values that will be added up in each square bin, according to their `positions`.
  - `grid::SquareGrid`: A square grid.
  - `total::Bool=true`: If the sum (`total` = true) or the mean (`total` = false) of `values` will be computed for each bin.
  - `empty_nan::Bool=true`: If NaN will be put into empty bins, 0 is used otherwise.

# Returns

  - A matrix with the histogram values.
"""
function histogram2D(
    positions::Matrix{<:Number},
    values::Vector{<:Number},
    grid::SquareGrid;
    total::Bool=true,
    empty_nan::Bool=true,
)::Matrix{<:Number}

    (
        length(values) == size(positions, 2) ||
        throw(ArgumentError("histogram2D: `values` must have as many elements as `positions` \
        has columns, but I got length(`values`) = $(length(values)) and size(`positions`, 2) = \
        $(size(positions, 2))"))
    )

    # Compute half of the bin size
    h_bin_width = grid.bin_width * 0.5

    x_edges = range(grid.x_bins[1] - h_bin_width, grid.x_bins[end] + h_bin_width, grid.n_bins + 1)
    y_edges = range(grid.y_bins[1] - h_bin_width, grid.y_bins[end] + h_bin_width, grid.n_bins + 1)

    return histogram2D(positions, values, collect(x_edges), collect(y_edges); total, empty_nan)

end

"""
    histogram2D(
        positions::Matrix{<:Number},
        x_edges::Vector{<:Number},
        y_edges::Vector{<:Number};
        <keyword arguments>
    )::Matrix{Float64}

Compute a 2D histogram of `positions`.

# Arguments

  - `positions::Matrix{<:Number}`: Values for which the histogram will be constructed.
  - `x_edges::Vector{<:Number}`: A list of bin edges for the x axis.
  - `y_edges::Vector{<:Number}`: A list of bin edges for the y axis.
  - `empty_nan::Bool=true`: If NaN will be put into empty bins, 0 is used otherwise.

# Returns

  - A matrix with the counts.
"""
function histogram2D(
    positions::Matrix{<:Number},
    x_edges::Vector{<:Number},
    y_edges::Vector{<:Number};
    empty_nan::Bool=true,
)::Matrix{Float64}

    issorted(x_edges) || sort!(x_edges)
    issorted(y_edges) || sort!(y_edges)

    n_x_bins = length(x_edges) - 1
    n_y_bins = length(y_edges) - 1

    x_borders = (first(x_edges), last(x_edges))
    y_borders = (first(y_edges), last(y_edges))

    histogram = zeros(Float64, (n_x_bins, n_y_bins))

    # Compute the histogram, ignoring NaNs and positions outside the range
    for point in eachcol(positions)

        !any(isnan, point) || continue

        x = point[1]
        y = point[2]

        !(x > x_borders[2] || x < x_borders[1]) || continue
        !(y > y_borders[2] || y < y_borders[1]) || continue

        if x == x_borders[1]
            i_x = 1
        else
            i_x = searchsortedfirst(x_edges, x) - 1
        end

        if y == y_borders[1]
            i_y = 1
        else
            i_y = searchsortedfirst(y_edges, y) - 1
        end

        histogram[n_y_bins - i_y + 1, i_x] += 1.0

    end

    empty_nan && replace!(x -> iszero(x) ? NaN : x, histogram)

    return histogram

end

"""
    findIn2DGrid(positions::Matrix{<:Number}, grid::SquareGrid)::Vector{Int}

Compute the linear index of the given `grid` where each of the `positions` falls.

!!! note

    For `positions` outside the grid a index of 0 is returned.

# Arguments

  - `positions::Matrix{<:Number}`: Target coordinates to be located in the grid.
  - `grid::SquareGrid`: A square grid.

# Returns

  - A vector with the linear indices.
"""
function findIn2DGrid(positions::Matrix{<:Number}, grid::SquareGrid)::Vector{Int}

    # Compute half of the bin size
    h_bin_width = grid.bin_width * 0.5

    # Compute the physical position of the grid borders
    x_borders = (grid.x_bins[1] - h_bin_width, grid.x_bins[end] + h_bin_width)
    y_borders = (grid.y_bins[1] - h_bin_width, grid.y_bins[end] + h_bin_width)

    # Compute the edges of the grid
    x_edges = range(x_borders[1], x_borders[end], grid.n_bins + 1)
    y_edges = range(y_borders[1], y_borders[end], grid.n_bins + 1)

    coordinates = zeros(Int, size(positions, 2))

    # Compute the linear indices, ignoring NaNs and positions outside the grid
    Threads.@threads for (idx, point) in collect(pairs(eachcol(positions)))

        !any(isnan, point) || continue

        x = point[1]
        y = point[2]

        !(x > x_borders[2] || x < x_borders[1]) || continue
        !(y > y_borders[2] || y < y_borders[1]) || continue

        if x == x_borders[1]
            i_x = 1
        else
            i_x = searchsortedfirst(x_edges, x) - 1
        end

        if y == y_borders[1]
            i_y = 1
        else
            i_y = searchsortedfirst(y_edges, y) - 1
        end

        coordinates[idx] = LinearIndices(grid.grid)[grid.n_bins - i_y + 1, i_x]

    end

    return coordinates

end

####################################################################################################
# 3D histogram utilities
####################################################################################################

"""
    listHistogram3D(positions::Matrix{<:Number}, grid::CubicGrid)::Array{Vector{Int},3}

Compute a 3D histogram of `positions`, returning the list of indices that fall within each bin.

# Arguments

  - `positions::Matrix{<:Number}`: Positions of the points. Each column corresponds to a point and each row is a dimension.
  - `grid::CubicGrid`: A cubic grid.

# Returns

    - A tensor with the indices of the points within each bin.
"""
function listHistogram3D(
    positions::Matrix{<:Number},
    grid::CubicGrid,
)::Array{Vector{Int},3}

    # Compute half of the bin size
    h_bin_width = grid.bin_width * 0.5

    # Compute the physical position of the grid borders
    x_borders = (grid.x_bins[1] - h_bin_width, grid.x_bins[end] + h_bin_width)
    y_borders = (grid.y_bins[1] - h_bin_width, grid.y_bins[end] + h_bin_width)
    z_borders = (grid.z_bins[1] - h_bin_width, grid.z_bins[end] + h_bin_width)

    # Compute the edges of the grid
    x_edges = range(x_borders[1], x_borders[end], grid.n_bins + 1)
    y_edges = range(y_borders[1], y_borders[end], grid.n_bins + 1)
    z_edges = range(z_borders[1], z_borders[end], grid.n_bins + 1)

    histogram = Array{Vector{Int}}(undef, size(grid.grid))
    for i in eachindex(histogram)
        histogram[i] = Int[]
    end

    # Compute the histogram, ignoring NaNs and positions outside the grid
    for (idx, point) in pairs(eachcol(positions))

        !any(isnan, point) || continue

        x = point[1]
        y = point[2]
        z = point[3]

        !(x > x_borders[2] || x < x_borders[1]) || continue
        !(y > y_borders[2] || y < y_borders[1]) || continue
        !(z > z_borders[2] || z < z_borders[1]) || continue

        if x == x_borders[1]
            i_x = 1
        else
            i_x = searchsortedfirst(x_edges, x) - 1
        end

        if y == y_borders[1]
            i_y = 1
        else
            i_y = searchsortedfirst(y_edges, y) - 1
        end

        if z == z_borders[1]
            i_z = 1
        else
            i_z = searchsortedfirst(z_edges, z) - 1
        end

        push!(histogram[grid.n_bins - i_y + 1, i_x, i_z], idx)

    end

    return histogram

end

"""
    findIn3DGrid(positions::Matrix{<:Number}, grid::CubicGrid)::Vector{Int}

Compute the linear index of the given `grid` where each of the `positions` falls.

!!! note

    For `positions` outside the grid a index of 0 is returned.

# Arguments

  - `positions::Matrix{<:Number}`: Target coordinates to be located in the grid.
  - `grid::CubicGrid`: A cubic grid.

# Returns

  - A vector with the linear indices.
"""
function findIn3DGrid(positions::Matrix{<:Number}, grid::CubicGrid)::Vector{Int}

    # Compute half of the bin size
    h_bin_width = grid.bin_width * 0.5

    # Compute the physical position of the grid borders
    x_borders = (grid.x_bins[1] - h_bin_width, grid.x_bins[end] + h_bin_width)
    y_borders = (grid.y_bins[1] - h_bin_width, grid.y_bins[end] + h_bin_width)
    z_borders = (grid.z_bins[1] - h_bin_width, grid.z_bins[end] + h_bin_width)

    # Compute the edges of the grid
    x_edges = range(x_borders[1], x_borders[end], grid.n_bins + 1)
    y_edges = range(y_borders[1], y_borders[end], grid.n_bins + 1)
    z_edges = range(z_borders[1], z_borders[end], grid.n_bins + 1)

    coordinates = zeros(Int, size(positions, 2))

    # Compute the histogram, ignoring NaNs and positions outside the grid
    Threads.@threads for (idx, point) in collect(pairs(eachcol(positions)))

        !any(isnan, point) || continue

        x = point[1]
        y = point[2]
        z = point[3]

        !(x > x_borders[2] || x < x_borders[1]) || continue
        !(y > y_borders[2] || y < y_borders[1]) || continue
        !(z > z_borders[2] || z < z_borders[1]) || continue

        if x == x_borders[1]
            i_x = 1
        else
            i_x = searchsortedfirst(x_edges, x) - 1
        end

        if y == y_borders[1]
            i_y = 1
        else
            i_y = searchsortedfirst(y_edges, y) - 1
        end

        if z == z_borders[1]
            i_z = 1
        else
            i_z = searchsortedfirst(z_edges, z) - 1
        end

        coordinates[idx] = LinearIndices(grid.grid)[grid.n_bins - i_y + 1, i_x, i_z]

    end

    return coordinates

end

"""
    histogram3D(
        positions::Matrix{<:Number},
        values::Vector{<:Number},
        grid::CubicGrid;
        <keyword arguments>
    )::Array{<:Number,3}

Compute a 3D histogram of `values`.

# Arguments

  - `positions::Matrix{<:Number}`: Positions of the values in the grid. Each column corresponds to a value and each row is a dimension.
  - `values::Vector{<:Number}`: The values that will be added up in each square bin, according to their `positions`.
  - `grid::CubicGrid`: A cubic grid.
  - `total::Bool=true`: If the sum (`total` = true) or the mean (`total` = false) of `values` will be computed for each bin.
  - `empty_nan::Bool=true`: If NaN will be put into empty bins, 0 is used otherwise.

# Returns

  - A 3D tensor with the histogram values.
"""
function histogram3D(
    positions::Matrix{<:Number},
    values::Vector{<:Number},
    grid::CubicGrid;
    total::Bool=true,
    empty_nan::Bool=true,
)::Array{<:Number,3}

    (
        length(values) == size(positions, 2) ||
        throw(ArgumentError("histogram3D: `values` must have as many elements as `positions` \
        has columns, but I got length(`values`) = $(length(values)) and size(`positions`, 2) = \
        $(size(positions, 2))"))
    )

    # Compute half of the bin size
    h_bin_width = grid.bin_width * 0.5

    # Compute the physical position of the grid borders
    x_borders = (grid.x_bins[1] - h_bin_width, grid.x_bins[end] + h_bin_width)
    y_borders = (grid.y_bins[1] - h_bin_width, grid.y_bins[end] + h_bin_width)
    z_borders = (grid.z_bins[1] - h_bin_width, grid.z_bins[end] + h_bin_width)

    # Compute the edges of the grid
    x_edges = range(x_borders[1], x_borders[end], grid.n_bins + 1)
    y_edges = range(y_borders[1], y_borders[end], grid.n_bins + 1)
    z_edges = range(z_borders[1], z_borders[end], grid.n_bins + 1)

    histogram = zeros(typeof(first(values)), size(grid.grid))
    counts = zeros(Int, size(grid.grid))

    # Compute the histogram, ignoring NaNs and positions outside the grid
    for (i, point) in pairs(eachcol(positions))

        !isnan(values[i]) || continue
        !any(isnan, point) || continue

        x = point[1]
        y = point[2]
        z = point[3]

        !(x > x_borders[2] || x < x_borders[1]) || continue
        !(y > y_borders[2] || y < y_borders[1]) || continue
        !(z > z_borders[2] || z < z_borders[1]) || continue

        if x == x_borders[1]
            i_x = 1
        else
            i_x = searchsortedfirst(x_edges, x) - 1
        end

        if y == y_borders[1]
            i_y = 1
        else
            i_y = searchsortedfirst(y_edges, y) - 1
        end

        if z == z_borders[1]
            i_z = 1
        else
            i_z = searchsortedfirst(z_edges, z) - 1
        end

        histogram[grid.n_bins - i_y + 1, i_x, i_z] += values[i]
        counts[grid.n_bins - i_y + 1, i_x, i_z] += 1

    end

    if empty_nan
        # Set empty bins to NaN
        nan = NaN * unit(first(values))
        Threads.@threads for i in eachindex(histogram)
            if iszero(counts[i])
                histogram[i] = nan
            end
        end
    end

    if !total
        # Compute the mean value instead of just the sum for each bin
        Threads.@threads for i in eachindex(histogram)
            if !iszero(counts[i])
                histogram[i] /= counts[i]
            end
        end
    end

    return histogram

end

"""
    project3DList(
        histogram::Array{Vector{Int},3},
        projection_plane::Symbol,
    )::Array{Vector{Int},2}

Project a 3D list of indices into a given plane.

# Arguments

  - `histogram::Array{Vector{Int},3}`: A tensor with the indices within each 3D bin. See [`listHistogram3D`](@ref).
  - `projection_plane::Symbol`: Projection plane. The options are `:xy`, `:xz`, and `:yz`.

# Returns

    - A matrix with the indices within each 2D bin.
"""
function project3DList(
    histogram::Array{Vector{Int},3},
    projection_plane::Symbol,
)::Array{Vector{Int},2}

    if projection_plane == :xy
        dims = 3
    elseif projection_plane == :xz
        dims = 2
    elseif projection_plane == :yz
        dims = 1
    else
        throw(ArgumentError("project3DList: The argument `projection_plane` must be \
        :xy, :xz or :yz, but I got :$(projection_plane)"))
    end

    return dropdims(reduce(vcat, histogram; dims, init=Int[]); dims)

end

####################################################################################################
# Histogram utilities
####################################################################################################

"""
    computeProfile(
        positions::Matrix{<:Unitful.Length},
        quantity::Vector{<:Number},
        grid::LinearGrid;
        <keyword arguments>
    )::Vector{<:Number}

Compute a profile of `quantity`, using an 1D histogram.

# Arguments

  - `positions::Matrix{<:Unitful.Length}`: Positions of the cells/particles. Each column is a cell/particle and each row a dimension.
  - `quantity::Vector{<:Number}`: The profile will be of this quantity.
  - `grid::LinearGrid`: Linear grid.
  - `norm::Vector{<:Number}=Number[]`: The value of `quantity` in each bin will be divided by the corresponding value of `norm`.
  - `flat::Bool=true`: If the profile will be 2D (rings), or 3D (spherical shells).
  - `total::Bool=true`: If the sum (default) or the mean of `quantity` will be computed for each bin.
  - `cumulative::Bool=false`: If the profile will be accumulated.
  - `density::Bool=false`: If the profile will be of the density of `quantity`.
  - `empty_nan::Bool=true`: If empty bins will be set to NaN. 0 is used otherwise. Notice that if `empty_nan` = true and `cumulative` = true, every bin after the first NaN will be set to NaN.

# Returns

  - Vector with the values of the profile.
"""
function computeProfile(
    positions::Matrix{<:Unitful.Length},
    quantity::Vector{<:Number},
    grid::LinearGrid;
    norm::Vector{<:Number}=Number[],
    flat::Bool=true,
    total::Bool=true,
    cumulative::Bool=false,
    density::Bool=false,
    empty_nan::Bool=true,
)::Vector{<:Number}

    if isempty(quantity)
        (
            logging[] &&
            @warn("computeProfile: `quantity` is empty. The profile will be filled with NaNs")
        )
        return fill(NaN, grid.size)
    end

    # Compute the distances of the cells/particles to the origin of the grid
    if flat
        distances = colwise(Euclidean(), positions[1:2, :], grid.origin[1:2])
    else
        distances = colwise(Euclidean(), positions, grid.origin)
    end

    # Compute the histogram of `quantity`
    if isempty(norm)

        profile = histogram1D(distances, quantity, grid; total, empty_nan)

    else

        quantity_histogram = histogram1D(distances, quantity, grid; total, empty_nan)
        norm_histogram = histogram1D(distances, norm, grid; total, empty_nan=false)

        replace!(x -> iszero(x) ? oneunit(x) : x, norm_histogram)

        profile = quantity_histogram ./ norm_histogram

    end

    region = flat ? grid.bin_size_2D : grid.bin_size_3D

    if cumulative
        return density ? cumsum(profile) ./ cumsum(region) : cumsum(profile)
    end

    return density ? profile ./ region : profile

end

"""
    computeBandProfile(
        positions::Matrix{<:Unitful.Length},
        quantity::Vector{<:Number},
        grid::LinearGrid;
        <keyword arguments>
    )::NTuple{3,Vector{<:Number}}

Compute a profile of `quantity`, using an 1D histogram and three given aggregator functions.

Each aggregator functions with the signature:

    `agg_func(::Vector{<:Number}) -> ::Number`

will be used to accumulate the values of `quantity` within each bin.

# Arguments

  - `positions::Matrix{<:Unitful.Length}`: Positions of the cells/particles. Each column is a cell/particle and each row a dimension.
  - `quantity::Vector{<:Number}`: The profile will be of this quantity.
  - `grid::LinearGrid`: Linear grid.
  - `flat::Bool=true`: If the profile will be 2D, using rings, or 3D, using spherical shells.
  - `density::Bool=false`: If the profile will be of the density of `quantity`.
  - `center_func::Function=x->quantile(x, 0.5)`: Aggregator function for the central value.
  - `low_func::Function=x->quantile(x, 0.25)`: Aggregator function for the low value.
  - `high_func::Function=x->quantile(x, 0.75)`: Aggregator function for the high value.

# Returns

  - A tuple with three elements:

      + A vector with the central value for each bin.
      + A vector with the low value for each bin.
      + A vector with the high value for each bin.
"""
function computeBandProfile(
    positions::Matrix{<:Unitful.Length},
    quantity::Vector{<:Number},
    grid::LinearGrid;
    flat::Bool=true,
    density::Bool=false,
    center_func::Function=x->quantile(x, 0.5),
    low_func::Function=x->quantile(x, 0.25),
    high_func::Function=x->quantile(x, 0.75),
)::NTuple{3,Vector{<:Number}}

    if isempty(quantity)
        (
            logging[] &&
            @warn("computeBandProfile: `quantity` is empty. The profile will be filled with NaNs")
        )
        return fill(NaN, grid.size)
    end

    # Compute the distances of the cells/particles to the origin of the grid
    if flat
        distances = colwise(Euclidean(), positions[1:2, :], grid.origin[1:2])
    else
        distances = colwise(Euclidean(), positions, grid.origin)
    end

    # Compute the histogram of `quantity`
    histogram = listHistogram1D(distances, quantity, grid)

    region = flat ? grid.bin_size_2D : grid.bin_size_3D

    if density
        center = center_func.(histogram) ./ region
        low    = low_func.(histogram) ./ region
        high   = high_func.(histogram) ./ region
    else
        center = center_func.(histogram)
        low    = low_func.(histogram)
        high   = high_func.(histogram)
    end

    return center, low, high

end
