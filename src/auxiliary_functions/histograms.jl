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

    histogram = zeros(eltype(values), n_bins)
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

    histogram = zeros(eltype(values), (n_x_bins, n_y_bins))
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
            idx_x = 1
        else
            idx_x = searchsortedfirst(x_edges, x) - 1
        end

        if y == y_borders[1]
            idx_y = 1
        else
            idx_y = searchsortedfirst(y_edges, y) - 1
        end

        histogram[idx_x, idx_y] += values[i]
        counts[idx_x, idx_y] += 1

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

    return histogram2D(positions, values, grid.x_edges, grid.y_edges; total, empty_nan)

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

  - `positions::Matrix{<:Number}`: Values for which the histogram will be constructed. Each column corresponds to a value and each row is a dimension.
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
            idx_x = 1
        else
            idx_x = searchsortedfirst(x_edges, x) - 1
        end

        if y == y_borders[1]
            idx_y = 1
        else
            idx_y = searchsortedfirst(y_edges, y) - 1
        end

        histogram[idx_x, idx_y] += 1.0

    end

    empty_nan && replace!(x -> iszero(x) ? NaN : x, histogram)

    return histogram

end

"""
    findIn2DGrid(
        positions::Matrix{<:Number},
        grid::SquareGrid;
        <keyword arguments>
    )::VecOrMat{Int}

Compute the linear or cartesian index of the pixel where each of the `positions` is located.

!!! note

    For `positions` outside the grid 0 is returned.

# Arguments

  - `positions::Matrix{<:Number}`: Target coordinates to be located in the grid. Each column corresponds to a value and each row is a dimension.
  - `grid::SquareGrid`: A square grid.
  - `cartesian::Bool=false`: If the resulting indices will be cartesian or linear.

# Returns

  - A vector with the linear indices or a matrix with cartesian indices.
"""
function findIn2DGrid(
    positions::Matrix{<:Number},
    grid::SquareGrid;
    cartesian::Bool=false,
)::VecOrMat{Int}

    # Compute the physical position of the grid borders
    x_borders = (first(grid.x_edges), last(grid.x_edges))
    y_borders = (first(grid.y_edges), last(grid.y_edges))

    if cartesian
        coordinates = zeros(Int, (2, size(positions, 2)))
    else
        coordinates = zeros(Int, size(positions, 2))
    end

    grid_linear_indices = LinearIndices(grid.n_bins)

    # Compute the linear indices, ignoring NaNs and positions outside the grid
    Threads.@threads for (idx, point) in collect(pairs(eachcol(positions)))

        !any(isnan, point) || continue

        x = point[1]
        y = point[2]

        !(x > x_borders[2] || x < x_borders[1]) || continue
        !(y > y_borders[2] || y < y_borders[1]) || continue

        if x == x_borders[1]
            idx_x = 1
        else
            idx_x = searchsortedfirst(grid.x_edges, x) - 1
        end

        if y == y_borders[1]
            idx_y = 1
        else
            idx_y = searchsortedfirst(grid.y_edges, y) - 1
        end

        if cartesian
            coordinates[:, idx] .= (idx_x, idx_y)
        else
            coordinates[idx] = grid_linear_indices[idx_x, idx_y]
        end

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

    # Compute the physical position of the grid borders
    x_borders =  (first(grid.x_edges), last(grid.x_edges))
    y_borders =  (first(grid.y_edges), last(grid.y_edges))
    z_borders =  (first(grid.z_edges), last(grid.z_edges))

    histogram = Array{Vector{Int}}(undef, grid.n_bins)
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
            idx_x = 1
        else
            idx_x = searchsortedfirst(grid.x_edges, x) - 1
        end

        if y == y_borders[1]
            idx_y = 1
        else
            idx_y = searchsortedfirst(grid.y_edges, y) - 1
        end

        if z == z_borders[1]
            idx_z = 1
        else
            idx_z = searchsortedfirst(grid.z_edges, z) - 1
        end

        push!(histogram[idx_x, idx_y, idx_z], idx)

    end

    return histogram

end

"""
    findIn3DGrid(
        positions::Matrix{<:Number},
        grid::CubicGrid;
        <keyword arguments>
    )::VecOrMat{Int}

Compute the linear or cartesian index of the voxel where each of the `positions` is located.

!!! note

    For `positions` outside the grid 0 is returned.

# Arguments

  - `positions::Matrix{<:Number}`: Target coordinates to be located in the grid. Each column corresponds to a value and each row is a dimension.
  - `grid::CubicGrid`: A cubic grid.
  - `cartesian::Bool=false`: If the resulting indices will be cartesian or linear.

# Returns

  - A vector with the linear indices or a matrix with cartesian indices.
"""
function findIn3DGrid(
    positions::Matrix{<:Number},
    grid::CubicGrid;
    cartesian::Bool=false,
)::VecOrMat{Int}

    # Compute the physical position of the grid borders
    x_borders =  (first(grid.x_edges), last(grid.x_edges))
    y_borders =  (first(grid.y_edges), last(grid.y_edges))
    z_borders =  (first(grid.z_edges), last(grid.z_edges))

    if cartesian
        coordinates = zeros(Int, (3, size(positions, 2)))
    else
        coordinates = zeros(Int, size(positions, 2))
    end

    grid_linear_indices = LinearIndices(grid.n_bins)

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
            idx_x = 1
        else
            idx_x = searchsortedfirst(grid.x_edges, x) - 1
        end

        if y == y_borders[1]
            idx_y = 1
        else
            idx_y = searchsortedfirst(grid.y_edges, y) - 1
        end

        if z == z_borders[1]
            idx_z = 1
        else
            idx_z = searchsortedfirst(grid.z_edges, z) - 1
        end

        if cartesian
            coordinates[:, idx] .= (idx_x, idx_y, idx_z)
        else
            coordinates[idx] = grid_linear_indices[idx_x, idx_y, idx_z]
        end

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

    # Compute the physical position of the grid borders
    x_borders =  (first(grid.x_edges), last(grid.x_edges))
    y_borders =  (first(grid.y_edges), last(grid.y_edges))
    z_borders =  (first(grid.z_edges), last(grid.z_edges))

    histogram = zeros(eltype(values), grid.n_bins)
    counts = zeros(Int, grid.n_bins)

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
            idx_x = 1
        else
            idx_x = searchsortedfirst(grid.x_edges, x) - 1
        end

        if y == y_borders[1]
            idx_y = 1
        else
            idx_y = searchsortedfirst(grid.y_edges, y) - 1
        end

        if z == z_borders[1]
            idx_z = 1
        else
            idx_z = searchsortedfirst(grid.z_edges, z) - 1
        end

        histogram[idx_x, idx_y, idx_z] += values[i]
        counts[idx_x, idx_y, idx_z] += 1

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
        return fill(NaN, grid.n_bins)
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
        return fill(NaN, grid.n_bins)
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
