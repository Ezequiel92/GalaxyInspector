####################################################################################################
# Grid utilities
####################################################################################################

"""
    scaledBins(
        values::Vector{<:Number},
        n_bins::Int;
        <keyword arguments>
    )::Vector{Float64}

Compute a set of bin edges, to encompass a given list of values.

# Arguments

  - `values::Vector{<:Number}`: Values to be binned.
  - `n_bins::Int`: Number of bins.
  - `scaling::Function=identity`: Scaling function. The options are the scaling functions accepted by [Makie](https://docs.makie.org/stable/): log10, log2, log, sqrt, Makie.logit, Makie.Symlog10, Makie.pseudolog10, and identity.
  - `limits::Tuple{<:Number,<:Number}=(-Inf, Inf)`: Set it to a value different than `(-Inf, Inf)` if you want to fix the limits of the bins.

# Returns

  - A sorted list of bin edges.
"""
function scaledBins(
    values::Vector{<:Number},
    n_bins::Int;
    scaling::Function=identity,
    limits::Tuple{<:Number,<:Number}=(-Inf, Inf),
)::Vector{<:Number}

    (
        !(limits[1] > limits[2]) ||
        throw(ArgumentError("scaledBins: `limits` must be (min, max), but I got `limits`[1] = \
        $(limits[1]) > `limits`[2] = $(limits[2])"))
    )

    # Compute the limits of the bins
    min = isinf(limits[1]) ? scaling(ustrip(minimum(values))) : scaling(ustrip(limits[1]))
    max = isinf(limits[2]) ? scaling(ustrip(maximum(values))) : scaling(ustrip(limits[2]))

    # For a small range, increase it by 0.2 * abs(max)
    if ((range = max - min) <= 1e-4 * abs(max)) && isinf(limits[1]) && isinf(limits[2])
        range += 0.2 * abs(max)
    end

    # Compute the width of the bins
    width = range / n_bins

    # Get the inverse function of `scaling`
    inverse = Makie.inverse_transform(scaling)

    # Get the unit of the values
    v_unit = unit(first(values))

    return [inverse(min + width * i) * v_unit for i in 0:n_bins]

end

"""
    deltas(data::Vector{<:Number})::Vector{<:Number}

Compute the difference between each consecutive pair of elements in `data`.

# Arguments

  - `data::Vector{<:Number}`: Data vector. It has to have at least 2 elements.

# Returns

  - A vector with the difference between each consecutive pair of elements in `data`, the first element is 0 by convention.
"""
function deltas(data::Vector{<:Number})::Vector{<:Number}

    # Allocate memory
    Δd = similar(data)

    # Get the number of elements in `data`
    nd = length(data)

    # Check that `data` has a valid length
    (
        nd >= 2 ||
        throw(ArgumentError("deltas: `data` must have at least 2 elements, but it has $(nd)"))
    )

    # Set the first value to 0, by convention
    Δd[1] = zero(data[1])

    for i in 2:nd
        Δd[i] = data[i] - data[i - 1]
    end

    return Δd

end

"""
    reduceResolution(
        hr_matrix::Matrix{<:Number},
        factor::Int;
        <keyword arguments>
    )::Matrix{<:Number}

Reduce the number of rows and columns of `hr_matrix` by `factor`, averaging or adding up its values.

# Arguments

  - `hr_matrix::Matrix{<:Number}`: Original "high resolution" matrix. It has to be a square matrix.
  - `factor::Int`: Factor by which the number of rows and columns will be reduced. It has to divide the size of `hr_matrix` exactly.
  - `total::Bool=false`: If the sum (`total` = true) or the mean (`total` = false) of the values in each of the old pixels will be used for the new pixels.

# Returns

  - The new smaller matrix, with the average values.
"""
function reduceResolution(
    hr_matrix::Matrix{<:Number},
    factor::Int;
    total::Bool=false,
)::Matrix{<:Number}

    !isone(factor) || return hr_matrix

    r, c = size(hr_matrix)
    (
        r == c ||
        throw(ArgumentError("reduceResolution: `hr_matrix` has to be a square matrix, but it has \
        $(c) columns and $(r) rows"))
    )

    (
        factor >= 1 ||
        throw(ArgumentError("reduceResolution: `factor` must be >= 1, \
        but I got `factor` = $(factor)"))
    )

    (
        r % factor == 0 ||
        throw(ArgumentError("reduceResolution: `factor` must divide the size of `hr_matrix` \
        exactly, but I got number of rows / `factor` = $(r / factor)"))
    )

    # Compute the size of the new matrix
    new_size = r ÷ factor

    # Allocate memory
    lr_matrix = zeros(eltype(hr_matrix), new_size, new_size)

    # Compute the number of old pixels per new pixel
    old_n_pixels = factor * factor

    for i in eachindex(lr_matrix)

        # Compute the row and column of the new matrix corresponding to index i
        row = mod1(i, new_size)
        col = ceil(Int, i / new_size)

        for j in (factor * (row - 1) + 1):(factor * row)
            for k in (factor * (col - 1) + 1):(factor * col)

                if !isnan(hr_matrix[j, k])
                    lr_matrix[i] += hr_matrix[j, k]
                end

            end
        end

    end

    return total ? lr_matrix : lr_matrix ./ old_n_pixels

end

"""
    reduceTicks(hr_ticks::Vector{<:Number}, factor::Int)::Vector{<:Number}

Reduce the number of ticks in `hr_ticks` by `factor` keeping the total length of the axis the same and assuming `hr_ticks` are regularly spaced.

# Arguments

  - `hr_ticks::Vector{<:Number}`: Original "high resolution" list of ticks.
  - `factor::Int`: Factor by which the number of ticks will be reduced. It has to divide the size of `hr_ticks` exactly.

# Returns

  - The new shorter tick list.
"""
function reduceTicks(hr_ticks::Vector{<:Number}, factor::Int)::Vector{<:Number}

    !isone(factor) || return hr_ticks

    l = length(hr_ticks)
    (
        l % factor == 0 ||
        throw(ArgumentError("reduceTicks: `factor` must divide the size of `hr_ticks` \
        exactly, but I got length(`hr_ticks`) / `factor` = $(l / factor)"))
    )

    (
        factor >= 1 ||
        throw(ArgumentError("reduceTicks: `factor` must be >= 1, but I got `factor` = $(factor)"))
    )

    # Compute the size of the new vector
    new_size = l ÷ factor

    # Allocate memory
    lr_ticks = similar(hr_ticks, new_size)

    if iseven(factor)

        shift = factor ÷ 2

        for i in eachindex(lr_ticks)

            idx = (i - 1) * factor + shift

            lr_ticks[i] = (hr_ticks[idx] + hr_ticks[idx + 1]) / 2.0

        end

    else

        shift = ceil(Int, factor / 2)

        for i in eachindex(lr_ticks)

            idx = (i - 1) * factor + shift

            lr_ticks[i] = hr_ticks[idx]

        end

    end

    return lr_ticks

end

"""
    projectIntoCircularGrid(
        image::Matrix{<:Number},
        n_bins::Int;
        <keyword arguments>
    )::Vector{<:Number}

Project `image` into a circular grid, averaging the values in each concentric ring.

# Arguments

  - `image::Matrix{<:Number}`: Original matrix. It has to be a square matrix.
  - `n_bins::Int`: Number of bins for the circular grid.
  - `inscribed::Bool=true`: If the circular grid will be inscribed in `image` when doing the projection. If set to false, the matrix will be inscribed into the circular grid instead.
  - `total::Bool=false`: If the sum (`total` = true) or the mean (`total` = false) of the values in each of the old pixels that fall within each ring will be used.

# Returns

  - A vector with the averages of the values in each concentric ring.
"""
function projectIntoCircularGrid(
    image::Matrix{<:Number},
    n_bins::Int;
    inscribed::Bool=true,
    total::Bool=false,
)::Vector{<:Number}

    r, c = size(image)

    (
        r == c ||
        throw(ArgumentError("projectIntoCircularGrid: `image` has to be a square matrix, \
        but it has $(c) columns and $(r) rows"))
    )

    # Construct a square grid center in (0, 0)
    square_grid = SquareGrid(1.0, r)

    # Construct a circular grid center in (0, 0)
    circular_grid = CircularGrid(inscribed ? 0.5 : sqrt(0.5), n_bins)

    # Compute the radial distance to the origin of each pixel in the square grid
    positions = norm.(vec(square_grid.grid))

    non_nan_idxs = findall(!isnan, vec(image))

    profile = histogram1D(positions, vec(image), circular_grid; total, empty_nan=false)

    return profile

end

"""
    flattenGrid(cubic_grid::CubicGrid)::SquareGrid

Using a `CubicGrid` construct a `SquareGrid` with the same center, number of bins, and physical side length.

# Arguments

  - `cubic_grid::CubicGrid`: Cubic grid.

# Returns

  - A square grid.
"""
function flattenGrid(cubic_grid::CubicGrid)::SquareGrid

    grid_size = cubic_grid.grid_size
    n_bins = cubic_grid.n_bins

    bin_width = grid_size / n_bins
    shift = 0.5 * (grid_size - bin_width)

    center = [cubic_grid.x_edges[1], cubic_grid.y_edges[1], cubic_grid.z_edges[1]] .+ shift

    return SquareGrid(grid_size, n_bins; center)

end
