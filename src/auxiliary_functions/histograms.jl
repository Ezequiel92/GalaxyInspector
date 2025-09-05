####################################################################################################
# Histogram utilities
####################################################################################################

"""
    listHistogram1D(
        positions::Vector{<:Number},
        values::Vector{<:Number},
        grid::Union{LinearGrid,CircularGrid},
    )::Vector{Vector{<:Number}}

Compute a 1D histogram of `values`, returning the full list of `values` within each bin.

# Arguments

  - `positions::Vector{<:Number}`: The positions of `values` within a 1D axis.
  - `values::Vector{<:Number}`: The values that will be sorted, according to their `positions`.
  - `grid::Union{LinearGrid,CircularGrid}`: A linear or circular grid.

# Returns

    - A vector with the lists of `values` within each bin.
"""
function listHistogram1D(
    positions::Vector{<:Number},
    values::Vector{<:Number},
    grid::Union{LinearGrid,CircularGrid},
)::Vector{Vector{<:Number}}

    (
        length(values) == length(positions) ||
        throw(ArgumentError("listHistogram1D: `values` must have as many elements as \
        `positions`, but I got length(`values`) = $(length(values)) and \
        length(`positions`) = $(length(positions))"))
    )

    n_bins = length(grid.grid)

    if grid.log

        l_u       = unit(first(grid.edges))
        positions = log10.(ustrip.(l_u, positions))
        p_min     = log10(ustrip(first(grid.edges)))
        p_max     = log10(ustrip(l_u, last(grid.edges)))

    else

        p_min = first(grid.edges)
        p_max = last(grid.edges)

    end

    # Compute the bin width
    width = (p_max - p_min) / n_bins

    # Allocate memory
    histogram = [eltype(values)[] for _ in 1:n_bins]

    # Compute the histogram, ignoring NaNs and positions outside the grid range
    for (position, value) in zip(positions, values)

        if isnan(position) || isnan(value)
            continue
        elseif position < p_min || p_max < position
            continue
        elseif position == p_min
            idx = 1
        elseif position == p_max
            idx = n_bins
        else
            idx = ceil(Int, (position - p_min) / width)
        end

        push!(histogram[idx], value)

    end

    return histogram

end

"""
    histogram1D(
        positions::Vector{<:Number},
        values::Vector{<:Number},
        grid::Union{LinearGrid,CircularGrid};
        <keyword arguments>
    )::Vector{<:Number}

Compute a 1D histogram of `values`.

# Arguments

  - `positions::Vector{<:Number}`: Positions of the `values` within a 1D axis.
  - `values::Vector{<:Number}`: The values that will be added up to each bin, according to their `positions`.
  - `grid::Union{LinearGrid,CircularGrid}`: A linear or circular grid.
  - `total::Bool=true`: If the sum (`total` = true) or the mean (`total` = false) of `values` will be computed for each bin.
  - `empty_nan::Bool=true`: If NaN will be put into empty bins, 0 is used otherwise.

# Returns

  - A vector with the histogram values.
"""
function histogram1D(
    positions::Vector{<:Number},
    values::Vector{<:Number},
    grid::Union{LinearGrid,CircularGrid};
    total::Bool=true,
    empty_nan::Bool=true,
)::Vector{<:Number}

    (
        length(values) == length(positions) ||
        throw(ArgumentError("histogram1D: `values` must have as many elements as \
        `positions`, but I got length(`values`) = $(length(values)) and \
        length(`positions`) = $(length(positions))"))
    )

    n_bins = length(grid.grid)

    if grid.log

        l_u       = unit(first(grid.edges))
        positions = log10.(ustrip.(l_u, positions))
        p_min     = log10(ustrip(grid.edges[1]))
        p_max     = log10(ustrip(l_u, grid.edges[end]))

    else

        p_min = grid.edges[1]
        p_max = grid.edges[end]

    end

    # Compute the bin width
    width = (p_max - p_min) / n_bins

    # Allocate memory
    histogram = zeros(typeof(first(values)), n_bins)
    counts = zeros(Int, n_bins)

    # Compute the histogram, ignoring NaNs and positions outside the grid range
    for (position, value) in zip(positions, values)

        if isnan(position) || isnan(value)
            continue
        elseif position < p_min || p_max < position
            continue
        elseif position == p_min
            idx = 1
        elseif position == p_max
            idx = n_bins
        else
            idx = ceil(Int, (position - p_min) / width)
        end

        histogram[idx] += value
        counts[idx] += 1

    end

    if empty_nan
        # Set empty bins to NaN
        nan = NaN * unit(first(values))
        for (i, count) in pairs(counts)
            if iszero(count)
                histogram[i] = nan
            end
        end
    end

    if !total
        # Compute the mean value instead of just the sum for each bin
        for (i, count) in pairs(counts)
            if !iszero(count)
                histogram[i] /= count
            end
        end
    end

    return histogram

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

    # Allocate memory
    histogram = zeros(typeof(first(values)), n_bins)
    counts = zeros(Int, n_bins)

    # Compute the histogram, ignoring NaNs and positions outside of range
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
        for (i, count) in pairs(counts)
            if iszero(count)
                histogram[i] = nan
            end
        end
    end

    if !total
        # Compute the mean value instead of just the sum for each bin
        for (i, count) in pairs(counts)
            if !iszero(count)
                histogram[i] /= count
            end
        end
    end

    return histogram

end

"""
    histogram1D(positions::Vector{<:Number}, grid::Union{LinearGrid,CircularGrid})::Vector{Int}

Compute a 1D histogram of `positions`.

# Arguments

  - `positions::Vector{<:Number}`: Values for which the histogram will be constructed.
  - `grid::Union{LinearGrid,CircularGrid}`: A linear or circular grid.

# Returns

  - A vector with the counts.
"""
function histogram1D(positions::Vector{<:Number}, grid::Union{LinearGrid,CircularGrid})::Vector{Int}

    n_bins = length(grid.grid)

    if grid.log
        l_u       = unit(first(grid.edges))
        positions = log10.(ustrip.(l_u, positions))
        p_min     = log10(ustrip(grid.edges[1]))
        p_max     = log10(ustrip(l_u, grid.edges[end]))
    else
        p_min = grid.edges[1]
        p_max = grid.edges[end]
    end

    # Compute the bin width
    width = (p_max - p_min) / n_bins

    # Allocate memory
    histogram = zeros(Int, n_bins)

    # Compute the histogram, ignoring NaNs and positions outside the grid range
    for position in positions

        if isnan(position)
            continue
        elseif position < p_min || p_max < position
            continue
        elseif position == p_min
            idx = 1
        elseif position == p_max
            idx = n_bins
        else
            idx = ceil(Int, (position - p_min) / width)
        end

        histogram[idx] += 1

    end

    return histogram

end

"""
    histogram1D(positions::Vector{<:Number}, edges::Vector{<:Number})::Vector{Int}

Compute a 1D histogram of `positions`.

# Arguments

  - `positions::Vector{<:Number}`: Values for which the histogram will be constructed.
  - `edges::Vector{<:Number}`: A list of bin edges.

# Returns

  - A vector with the counts.
"""
function histogram1D(positions::Vector{<:Number}, edges::Vector{<:Number})::Vector{Int}

    issorted(edges) || sort!(edges)

    n_bins = length(edges) - 1

    p_min = first(edges)
    p_max = last(edges)

    # Allocate memory
    histogram = zeros(Int, n_bins)

    # Compute the histogram, ignoring NaNs and positions outside of range
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

        histogram[idx] += 1

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

    # Compute the physical position of the grid borders
    x_borders = (grid.x_bins[1] - h_bin_width, grid.x_bins[end] + h_bin_width)
    y_borders = (grid.y_bins[1] - h_bin_width, grid.y_bins[end] + h_bin_width)

    # Allocate memory
    histogram = zeros(typeof(first(values)), size(grid.grid))
    counts = zeros(Int, size(grid.grid))

    # Compute the histogram, ignoring NaNs and positions outside of range
    for (i, point) in pairs(eachcol(positions))

        !isnan(values[i]) || continue
        !any(isnan, point) || continue

        x = point[1]
        y = point[2]

        !(x > x_borders[2] || x < x_borders[1]) || continue
        !(y > y_borders[2] || y < y_borders[1]) || continue

        if x == x_borders[1]
            i_x = 1
        elseif x == x_borders[2]
            i_x = grid.n_bins
        else
            i_x = ceil(Int, (x - x_borders[1]) / grid.bin_width)
        end

        if y == y_borders[1]
            i_y = 1
        elseif y == y_borders[2]
            i_y = grid.n_bins
        else
            i_y = ceil(Int, (y - y_borders[1]) / grid.bin_width)
        end

        histogram[grid.n_bins - i_y + 1, i_x] += values[i]
        counts[grid.n_bins - i_y + 1, i_x] += 1

    end

    if empty_nan
        # Set empty bins to NaN
        nan = NaN * unit(first(values))
        for i in eachindex(histogram)
            if iszero(counts[i])
                histogram[i] = nan
            end
        end
    end

    if !total
        # Compute the mean value instead of just the sum for each bin
        for i in eachindex(histogram)
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

    # Allocate memory
    histogram = zeros(typeof(first(values)), (n_x_bins, n_y_bins))
    counts = zeros(Int, (n_x_bins, n_y_bins))

    # Compute the histogram, ignoring NaNs and positions outside of range
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
        for i in eachindex(histogram)
            if iszero(counts[i])
                histogram[i] = nan
            end
        end
    end

    if !total
        # Compute the mean value instead of just the sum for each bin
        for i in eachindex(histogram)
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
        x_edges::Vector{<:Number},
        y_edges::Vector{<:Number},
    )::Matrix{Int}

Compute a 2D histogram of `positions`.

# Arguments

  - `positions::Matrix{<:Number}`: Values for which the histogram will be constructed.
  - `x_edges::Vector{<:Number}`: A list of bin edges for the x axis.
  - `y_edges::Vector{<:Number}`: A list of bin edges for the y axis.

# Returns

  - A matrix with the counts.
"""
function histogram2D(
    positions::Matrix{<:Number},
    x_edges::Vector{<:Number},
    y_edges::Vector{<:Number},
)::Matrix{Int}

    issorted(x_edges) || sort!(x_edges)
    issorted(y_edges) || sort!(y_edges)

    n_x_bins = length(x_edges) - 1
    n_y_bins = length(y_edges) - 1

    x_borders = (first(x_edges), last(x_edges))
    y_borders = (first(y_edges), last(y_edges))

    # Allocate memory
    histogram = zeros(Int, (n_x_bins, n_y_bins))

    # Compute the histogram, ignoring NaNs and positions outside of range
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

        histogram[n_y_bins - i_y + 1, i_x] += 1

    end

    return histogram

end

"""
    listHistogram3D(positions::Matrix{<:Number}, grid::CubicGrid)::Array{Vector{Int},3}

Compute a 3D histogram of `positions`, returning the list of indices that fall within each bin.

# Arguments

  - `positions::Matrix{<:Number}`: Positions of the points. Each column corresponds to a point and each row is a dimension.
  - `grid::CubicGrid`: A cubic grid.

# Returns

    - A tensor with the indices of the points within each bin.
"""
function listHistogram3D(positions::Matrix{<:Number}, grid::CubicGrid)::Array{Vector{Int},3}

    # Compute half of the bin size
    h_bin_width = grid.bin_width * 0.5

    # Compute the physical position of the grid borders
    x_borders = (grid.x_bins[1] - h_bin_width, grid.x_bins[end] + h_bin_width)
    y_borders = (grid.y_bins[1] - h_bin_width, grid.y_bins[end] + h_bin_width)
    z_borders = (grid.z_bins[1] - h_bin_width, grid.z_bins[end] + h_bin_width)

    # Allocate memory
    histogram = Array{Vector{Int}}(undef, size(grid.grid))
    for i in eachindex(histogram)
        histogram[i] = Int[]
    end

    # Compute the histogram, ignoring NaNs and positions outside of range
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
        elseif x == x_borders[2]
            i_x = grid.n_bins
        else
            i_x = ceil(Int, (x - x_borders[1]) / grid.bin_width)
        end

        if y == y_borders[1]
            i_y = 1
        elseif y == y_borders[2]
            i_y = grid.n_bins
        else
            i_y = ceil(Int, (y - y_borders[1]) / grid.bin_width)
        end

        if z == z_borders[1]
            i_z = 1
        elseif z == z_borders[2]
            i_z = grid.n_bins
        else
            i_z = ceil(Int, (z - z_borders[1]) / grid.bin_width)
        end

        push!(histogram[grid.n_bins - i_y + 1, i_x, i_z], idx)

    end

    return histogram

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

    # Allocate memory
    histogram = zeros(typeof(first(values)), size(grid.grid))
    counts = zeros(Int, size(grid.grid))

    # Compute the histogram, ignoring NaNs and positions outside of range
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
        elseif x == x_borders[2]
            i_x = grid.n_bins
        else
            i_x = ceil(Int, (x - x_borders[1]) / grid.bin_width)
        end

        if y == y_borders[1]
            i_y = 1
        elseif y == y_borders[2]
            i_y = grid.n_bins
        else
            i_y = ceil(Int, (y - y_borders[1]) / grid.bin_width)
        end

        if z == z_borders[1]
            i_z = 1
        elseif z == z_borders[2]
            i_z = grid.n_bins
        else
            i_z = ceil(Int, (z - z_borders[1]) / grid.bin_width)
        end

        histogram[grid.n_bins - i_y + 1, i_x, i_z] += values[i]
        counts[grid.n_bins - i_y + 1, i_x, i_z] += 1

    end

    if empty_nan
        # Set empty bins to NaN
        nan = NaN * unit(first(values))
        for i in eachindex(histogram)
            if iszero(counts[i])
                histogram[i] = nan
            end
        end
    end

    if !total
        # Compute the mean value instead of just the sum for each bin
        for i in eachindex(histogram)
            if !iszero(counts[i])
                histogram[i] /= counts[i]
            end
        end
    end

    return histogram

end
