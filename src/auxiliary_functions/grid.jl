####################################################################################################
# Grid utilities
####################################################################################################

"""
    scaledBins(
        values::Vector{<:Number},
        n_bins::Int;
        <keyword arguments>
    )::Vector{Float64}

Compute a set of bin edges, that will encompass a given list of values.

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
    deltas(values::Vector{<:Number})::Vector{<:Number}

Compute the difference between each consecutive pair of elements in `values`.

# Arguments

  - `values::Vector{<:Number}`: It should have at least 2 elements.

# Returns

  - A vector with the difference between each consecutive pair of elements in `values`, the first element is 0 by convention.
"""
function deltas(values::Vector{<:Number})::Vector{<:Number}

    # Allocate memory
    Δd = similar(values)

    # Get the number of elements in `values`
    nd = length(values)

    # Check that `values` has a valid length
    (
        nd >= 2 ||
        throw(ArgumentError("deltas: `values` must have at least 2 elements, but it has $(nd)"))
    )

    # Set the first value to 0
    Δd[1] = zero(first(values))

    for i in 2:nd
        Δd[i] = values[i] - values[i - 1]
    end

    return Δd

end

"""
    reduceMatrix(
        hr_matrix::Matrix{<:Number},
        factor::Int;
        <keyword arguments>
    )::Matrix{<:Number}

Reduce the number of rows and columns of a given square matrix by `factor`, averaging or adding up its values.

# Arguments

  - `hr_matrix::Matrix{<:Number}`: Original "high resolution" matrix. It has to be a square matrix.
  - `factor::Int`: Factor by which the number of rows and columns will be reduced. It has to divide the size of `hr_matrix` exactly.
  - `total::Bool=false`: If the sum (`total` = true) or the mean (`total` = false) of the high resolution values will be used for the new matrix.

# Returns

  - The new smaller matrix.
"""
function reduceMatrix(hr_matrix::Matrix{<:Number}, factor::Int; total::Bool=false)::Matrix{<:Number}

    !isone(factor) || return hr_matrix

    r, c = size(hr_matrix)

    (
        r == c ||
        throw(ArgumentError("reduceMatrix: `hr_matrix` has to be a square matrix, but it has \
        $(c) columns and $(r) rows"))
    )

    (
        factor >= 1 ||
        throw(ArgumentError("reduceMatrix: `factor` must be >= 1, but I got `factor` = $(factor)"))
    )

    (
        r % factor == 0 ||
        throw(ArgumentError("reduceMatrix: `factor` must divide the size of `hr_matrix` \
        exactly, but I got number of rows / `factor` = $(r / factor)"))
    )

    # Compute the size of the new matrix
    new_size = r ÷ factor

    # Allocate memory
    lr_matrix = zeros(eltype(hr_matrix), new_size, new_size)

    # Compute the number of values in `hr_matrix` per position in the new matrix
    old_n_pixels = factor * factor

    for i in eachindex(lr_matrix)

        # Compute the row and column index corresponding to global index i in the new matrix
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
    projectIntoCircularGrid(
        image::Matrix{<:Number},
        n_bins::Int;
        <keyword arguments>
    )::Vector{<:Number}

Project a given matrix into a circular grid, averaging or adding up the values in each concentric ring.

# Arguments

  - `image::Matrix{<:Number}`: Original matrix. It has to be a square matrix.
  - `n_bins::Int`: Number of bins for the circular grid.
  - `inscribed::Bool=true`: If the circular grid will be inscribed in `image` when doing the projection. If set to false, the matrix will be inscribed into the circular grid instead.
  - `total::Bool=false`: If the sum (`total` = true) or the mean (`total` = false) of the values in `image`, that fall within each ring, will be used.

# Returns

  - A vector with the average or sum of the values that fall within each concentric ring.
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

    (
        n_bins >= 1 ||
        throw(ArgumentError("projectIntoCircularGrid: `n_bins` must be >= 1, \
        but I got `n_bins` = $(n_bins)"))
    )

    # Construct a square grid centered at (0, 0)
    square_grid = SquareGrid(1.0, r)

    # Construct a circular grid centered at (0, 0)
    circular_grid = CircularGrid(inscribed ? 0.5 : sqrt(0.5), n_bins)

    # Compute the radial distance of each pixel in the square grid to the origin
    positions = norm.(vec(square_grid.grid))

    profile = histogram1D(positions, vec(image), circular_grid; total, empty_nan=false)

    return profile

end

"""
    flattenGrid(cubic_grid::CubicGrid)::SquareGrid

From a `CubicGrid` construct a `SquareGrid` with the same center, number of bins, and physical side length.

# Arguments

  - `cubic_grid::CubicGrid`: Cubic grid.

# Returns

  - A square grid.
"""
function flattenGrid(cubic_grid::CubicGrid)::SquareGrid

    return SquareGrid(cubic_grid.grid_size, cubic_grid.n_bins; center=cubic_grid.center)

end
