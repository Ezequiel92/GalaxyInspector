####################################################################################################
# General utilities.
####################################################################################################

"""
Area of a circle with radius `r`.
"""
area(r::Number)::Number = π * r * r

"""
Volume of a sphere with radius `r`.
"""
volume(r::Number)::Number = π * r * r * r * 1.333

"""
    ring(vec::Vector, index::Integer)

Make the indexing operation `vec[index]` work using modular arithmetic for the indices.

# Arguments

  - `vec::Vector`: Vector.
  - `index::Integer`: Index.

# Returns

  - `vec[mod1(index, length(vec))]`

# Examples

```julia-repl
julia> ring([1, 2, 3], 11)
2

julia> ring([1, 2, 3], 3)
3

julia> ring([1, 2, 3], -5)
1
```
"""
ring(vec::Vector, index::Integer) = vec[mod1(index, length(vec))]

"""
    safeSelect(vec::Vector, index::IndexType; warnings::Bool=true)

Make the indexing operation `vec[index]` ignore indices that are out of bounds.

# Arguments

  - `vec::Vector`: Vector.
  - `index::IndexType`: Indices. It can be an integer (a single element), a vector of integers (several elements), an `UnitRange` (e.g. 5:13), an `StepRange` (e.g. 5:2:13) or (:) (every element).
  - `warnings::Bool=true`: If a warning will be given when `index` has out of bounds indices.

# Returns

  - `vec[index (minus out of bounds indices)]`

# Examples

```julia-repl
julia> safeSelect([1, 2, 3], 11; warnings=false)
Int[]

julia> safeSelect([1, 2, 3], 1:5; warnings=false)
3-element Vector{Int}:
 1
 2
 3

julia> safeSelect([1, 2, 3], 1:3:10; warnings=false)
1

julia> safeSelect([1, 2, 3], [1, 2, 5, 9]; warnings=false)
2-element Vector{Int}:
 1
 2

julia> safeSelect([1, 2, 3], (:); warnings=false)
3-element Vector{Int}:
 1
 2
 3
```
"""
function safeSelect(vec::Vector, index::IndexType; warnings::Bool=true)

    index != (:) || return vec

    index_list = [index...]

    filter!(x -> x <= length(vec), index_list)

    (
        length(index_list) == length([index...]) ||
        !warnings || @warn("safeSelect: There are out of bounds indices")
    )

    if length(index_list) == 1
        return vec[index_list...]
    else
        return vec[index_list]
    end

end

"""
Create a copy of `list` with every negative value set to 0.
"""
setPositive(list::VecOrMat{<:Number}) = replace(x -> x >= zero(x) ? x : zero(x), list)

"""
Test for strict positivity.
"""
isPositive(x::Number)::Bool = x > zero(x)
isPositive(x::AbstractArray)::Bool = all(isPositive, x)
isPositive(x...)::Bool = all(isPositive, x)

"""
Extension of `Base.iszero` to compeare [`IndexType`](@ref) with interger 0.
"""
Base.iszero(x::IndexType)::Bool = x == 0

"""
Extension of `Base.isempty` to check for empty [LaTeXStrings](https://github.com/JuliaStrings/LaTeXStrings.jl).
"""
Base.isempty(l_str::LaTeXString)::Bool = l_str == L""

"""
Always returns `nothing`, for any type and number of arguments.
"""
getNothing(x...; y...)::Nothing = nothing

"""
Always returns an empty vector, for any type and number of arguments.
"""
getEmpty(x...; y...)::Vector = []

"""
    hvcatImages(
        blocks_per_row::Int,
        paths::Vector{String};
        <keyword arguments>
    )::Nothing

Join several images vertically and horizontally

The elements will fill the rows and columns starting at the top left, going from left to right and from top to bottom.

# Arguments

  - `blocks_per_row::Int`: Number of columns.
  - `paths::Vector{String}`: Paths to the images.
  - `output_path::String="./joined_image.png"`: Path to the output image.
"""
function hvcatImages(
    blocks_per_row::Int,
    paths::Vector{String};
    output_path::String="./joined_image.png",
)::Nothing

    new_image = hvcat(blocks_per_row, [load(path) for path in paths]...)

    save(output_path, new_image)

    return nothing

end

"""
    rangeCut!(
        data::Vector{<:Number},
        range::Tuple{<:Number,<:Number};
        <keyword arguments>
    )::Bool

Delete every element in `data` that is outside the given `range`.

# Arguments

  - `data::Vector{<:Number}`: Dataset that will be pruned.
  - `range::Tuple{<:Number,<:Number}`: The range in question.
  - `keep_edges::Bool=true`: If the edges of the range will be kept.
  - `min_left::Int=0`: Minimum number of values that need to be left after pruning to proceed with the transformation.

# Returns

  - If a transformation was performed.
"""
function rangeCut!(
    data::Vector{<:Number},
    range::Tuple{<:Number,<:Number};
    keep_edges::Bool=true,
    min_left::Int=0,
)::Bool

    # Shortcut computation for special cases
    !(isempty(data) || all(isinf.(range))) || return false

    if keep_edges

        # Check that after the transformation at least `min_left` elements will be left
        count(x -> range[1] <= x <= range[2], data) >= min_left || return false

        # Delete element outside of the provided range
        filter!(x -> range[1] <= x <= range[2], data)

    else

        # Check that after the transformation at least `min_left` elements will be left
        count(x -> range[1] < x < range[2], data) >= min_left || return false

        # Delete element outside of the provided range
        filter!(x -> range[1] < x < range[2], data)

    end

    return true

end

"""
    rangeCut!(
        m_data::Vector{<:Number},
        s_data::Vector,
        range::Tuple{<:Number,<:Number};
        <keyword arguments>
    )::Bool

Delete every element in `m_data` that is outside the given `range`.

Every corresponding element in `s_data` (i.e. with the same index) will be deleted too.

# Arguments

  - `m_data::Vector{<:Number}`: Master dataset that will be pruned.
  - `s_data::Vector`: Slave dataset that will be pruned according to which values of `m_data` are outside `range`.
  - `range::Tuple{<:Number,<:Number}`: The range in question.
  - `keep_edges::Bool=true`: If the edges of the range will be kept.
  - `min_left::Int=0`: Minimum number of values that need to be left in the master dataset after pruning to proceed with the transformation.

# Returns

  - If a transformation was performed.
"""
function rangeCut!(
    m_data::Vector{<:Number},
    s_data::Vector,
    range::Tuple{<:Number,<:Number};
    keep_edges::Bool=true,
    min_left::Int=0,
)::Bool

    # Shortcut computation for special cases
    !(isempty(m_data) || all(isinf.(range))) || return false

    (
        length(s_data) >= length(m_data) ||
        throw(ArgumentError("rangeCut!: `s_data` must have at least as many elements as `m_data`, \
        but I got length(s_data) = $(length(s_data)) < length(m_data) = $(length(m_data))"))
    )

    if keep_edges

        # Find the elements outside of the provided range
        idxs = map(x -> x < range[1] || x > range[2], m_data)

        # Check that after the transformation at least `min_left` elements will be left
        count(.!idxs) >= min_left || return false

        # Delete element outside of the provided range
        deleteat!(m_data, idxs)
        deleteat!(s_data, idxs)

    else

        # Find the elements outside of the provided range
        idxs = map(x -> x <= range[1] || x >= range[2], m_data)

        # Check that after the transformation at least `min_left` elements will be left
        count(.!idxs) >= min_left || return false

        # Delete element outside of the provided range
        deleteat!(m_data, idxs)
        deleteat!(s_data, idxs)

    end

    return true

end

"""
    sanitizeData!(
        data::Vector{<:Number};
        <keyword arguments>
    )::NTuple{2,Bool}

Do the following transformations over `data`, in order:

  - Trim it to fit within the domain of the function `func_domain`.
  - Trim it to fit within `range`.
  - Scale it down by a factor of 10^`exp_factor`.

By default, no transformation is done.

# Arguments

  - `data::Vector{<:Number}`: Dataset to be sanitized.
  - `func_domain::Function=identity`: `data` will be trimmed to fit within the domain of the function `func_domain`. The options are the scaling functions accepted by [Makie](https://docs.makie.org/stable/): log10, log2, log, sqrt, Makie.logit, Makie.Symlog10, Makie.pseudolog10, and identity.
  - `range::Tuple{<:Number,<:Number}=(-Inf, Inf)`: Every element in `data` that falls outside of `range` will be deleted.
  - `keep_edges::Bool=true`: If the edges of `range` will be kept.
  - `min_left::Int=0`: Minimum number of values that need to be left after each transformation to procced with it.
  - `exp_factor::Int=0`: Every element in `data` will be divided by 10^`exp_factor`.
  - `warnings::Bool=true`: If a warning will be given when `data` is a vector of Integers, which may cause wrong results when dividing by 10^`exp_factor`.

# Returns

  - A tuple with two flags:

      + If `data` was mutated to fit within the domain of `func_domain`.
      + If `data` was mutated to fit within `range`.
"""
function sanitizeData!(
    data::Vector{<:Number};
    func_domain::Function=identity,
    range::Tuple{<:Number,<:Number}=(-Inf, Inf),
    keep_edges::Bool=true,
    min_left::Int=0,
    exp_factor::Int=0,
    warnings::Bool=true,
)::NTuple{2,Bool}

    !isempty(data) || return false, false

    d_unit = unit(first(data))

    # Trim `data` to fit within the domain of `func_domain`
    if func_domain ∈ [identity, Makie.pseudolog10, Makie.Symlog10]

        domain_flag = false

    elseif func_domain == sqrt

        domain_flag = rangeCut!(data, (0.0, Inf) .* d_unit; keep_edges=true, min_left)

    elseif func_domain == Makie.logit

        domain_flag = rangeCut!(data, (0.0, 1.0) .* d_unit; keep_edges=false, min_left)

    elseif func_domain ∈ [log, log2, log10]

        domain_flag = rangeCut!(data, (0.0, Inf) .* d_unit; keep_edges=false, min_left)

    else

        throw(ArgumentError("sanitizeData!: The function $(func_domain) is not supported. See \
        the list of supported scaling functions in the [Makie](https://docs.makie.org/stable/) \
        documentation"))

    end

    # Trim `data` to fit within `range`
    range_flag = rangeCut!(data, range; keep_edges, min_left)

    (
        !(isa(data, Vector{<:Integer}) && !iszero(exp_factor) && warnings) ||
        @warn("sanitizeData!: Elements of `data` are of type `Integer`, this may result \
        in errors or unwanted truncation when using `exp_factor` != 0")
    )

    # Scale `data` down by a factor of 10^`exp_factor`
    iszero(exp_factor) || (data ./= exp10(exp_factor))

    return domain_flag, range_flag

end

"""
    sanitizeData!(
        x_data::Vector{<:Number},
        y_data::Vector{<:Number};
        <keyword arguments>
    )::NTuple{4,Bool}

Do the following transformations over `x_data` and `y_data`, in order:

  - Trim them to fit within the domain of the functions `func_domain[1]` and `func_domain[2]`, respectively.
  - Trim them to fit within `range[1]` and `range[2]`, respectively.
  - Scale them down by a factor 10^`exp_factor[1]` and 10^`exp_factor[2]`, respectively.

By default, no transformation is done.

!!! note

    The datasets must have the same length, and any operation that deletes an element, will delete the corresponding element (i.e. with the same index) in the other dataset, so that the dataset will remain of equal length.

# Arguments

  - `x_data::Vector{<:Number}`: First dataset to be sanitized.
  - `y_data::Vector{<:Number}`: Second dataset to be sanitized.
  - `func_domain::NTuple{2,Function}=(identity, identity)`: `x_data` will be trimmed to fit within the domain of the function `func_domain[1]`, and `y_data` will be trimmed to fit within the domain of the function `func_domain[2]`. The options are the scaling functions accepted by [Makie](https://docs.makie.org/stable/): log10, log2, log, sqrt, Makie.logit, Makie.Symlog10, Makie.pseudolog10, and identity.
  - `range::Tuple{Tuple{<:Number,<:Number},Tuple{<:Number,<:Number}}=((-Inf, Inf), (-Inf, Inf))`: Every element in `x_data` that falls outside of `range[1]` will be deleted, and every element in `y_data` that falls outside of `range[2]` will be deleted.
  - `keep_edges::NTuple{2,Bool}=(true, true)`: If the edges of each corresponding `range` will be kept.
  - `min_left::Int=0`: Minimum number of values that need to be left in each dataset after any of the transformations to procced with them.
  - `exp_factor::NTuple{2,Int}=(0, 0)`: Every element in `x_data` will be divided by 10^`exp_factor[1]`, and every element in `y_data` will be divided by 10^`exp_factor[2]`.
  - `warnings::Bool=true`: If a warning will be given when any of the datasets is a vector of Integers, which may cause wrong results when dividing by 10^`exp_factor`.

# Returns

  - A tuple with four flags:

      + If `x_data` was successfully modified to fit within the domain of `func_domain[1]`.
      + If `y_data` was successfully modified to fit within the domain of `func_domain[2]`.
      + If `x_data` was successfully modified to fit within `range[1]`.
      + If `y_data` was successfully modified to fit within `range[2]`.
"""
function sanitizeData!(
    x_data::Vector{<:Number},
    y_data::Vector{<:Number};
    func_domain::NTuple{2,Function}=(identity, identity),
    range::Tuple{Tuple{<:Number,<:Number},Tuple{<:Number,<:Number}}=((-Inf, Inf), (-Inf, Inf)),
    keep_edges::NTuple{2,Bool}=(true, true),
    min_left::Int=0,
    exp_factor::NTuple{2,Int}=(0, 0),
    warnings::Bool=true,
)::NTuple{4,Bool}

    (
        length(x_data) == length(y_data) ||
        throw(ArgumentError("sanitizeData!: `x_data` and `y_data` must have the same length, \
        but I got length(x_data) = $(length(x_data)) != length(y_data) = $(length(y_data))"))
    )

    x_unit = isempty(x_data) ? Unitful.NoUnits : unit(first(x_data))

    # Trim the data to fit within the domain of `func_domain[1]`
    if func_domain[1] ∈ [identity, Makie.pseudolog10, Makie.Symlog10]

        x_domain_flag = false

    elseif func_domain[1] == sqrt

        x_domain_flag = rangeCut!(x_data, y_data, (0.0, Inf) .* x_unit; keep_edges=true, min_left)

    elseif func_domain[1] == Makie.logit

        x_domain_flag = rangeCut!(x_data, y_data, (0.0, 1.0) .* x_unit; keep_edges=false, min_left)

    elseif func_domain[1] ∈ [log, log2, log10]

        x_domain_flag = rangeCut!(x_data, y_data, (0.0, Inf) .* x_unit; keep_edges=false, min_left)

    else

        throw(ArgumentError("sanitizeData!: The function $(func_domain[1]) is not supported. See \
        the list of supported scaling functions in the [Makie](https://docs.makie.org/stable/) \
        documentation"))

    end

    y_unit = isempty(y_data) ? Unitful.NoUnits : unit(first(y_data))

    # Trim the data to fit within the domain of `func_domain[2]`
    if func_domain[2] ∈ [identity, Makie.pseudolog10, Makie.Symlog10]

        y_domain_flag = false

    elseif func_domain[2] == sqrt

        y_domain_flag = rangeCut!(y_data, x_data, (0.0, Inf) .* y_unit; keep_edges=true, min_left)

    elseif func_domain[2] == Makie.logit

        y_domain_flag = rangeCut!(y_data, x_data, (0.0, 1.0) .* y_unit; keep_edges=false, min_left)

    elseif func_domain[2] ∈ [log, log2, log10]

        y_domain_flag = rangeCut!(y_data, x_data, (0.0, Inf) .* y_unit; keep_edges=false, min_left)

    else

        throw(ArgumentError("sanitizeData!: The function $(func_domain[2]) is not supported. See \
        the list of supported scaling functions in the [Makie](https://docs.makie.org/stable/) \
        documentation"))

    end

    # Trim data to fit within `range[1]`
    x_range_flag = rangeCut!(x_data, y_data, range[1]; keep_edges=keep_edges[1], min_left)

    # Trim data to fit within `range[2]`
    y_range_flag = rangeCut!(y_data, x_data, range[2]; keep_edges=keep_edges[2], min_left)

    (
        !(isa(x_data, Vector{<:Integer}) && !iszero(exp_factor[1]) && warnings) ||
        @warn("sanitizeData!: Elements of `x_data` are of type Integer, this may result \
        in errors or unwanted truncation when using `exp_factor[1]` != 0")
    )

    (
        !(isa(y_data, Vector{<:Integer}) && !iszero(exp_factor[2]) && warnings) ||
        @warn("sanitizeData!: Elements of `y_data` are of type Integer, this may result \
        in errors or unwanted truncation when using `exp_factor[2]` != 0")
    )

    # Scale the data down by the factors `exp_factor`
    iszero(exp_factor[1]) || (x_data ./= exp10(exp_factor[1]))
    iszero(exp_factor[2]) || (y_data ./= exp10(exp_factor[2]))

    return x_domain_flag, y_domain_flag, x_range_flag, y_range_flag

end

"""
    scaledBins(
        values::Vector{<:Number},
        n_bins::Int;
        <keyword arguments>
    )::Vector{Float64}

Compute a set of bin edges, for a given list values.

# Arguments

  - `values::Vector{<:Number}`: Values to be binned.
  - `n_bins::Int`: Number of bins.
  - `scaling::Function=identity`: Scaling function. The options are the scaling functions accepted by [Makie](https://docs.makie.org/stable/): log10, log2, log, sqrt, Makie.logit, Makie.Symlog10, Makie.pseudolog10, and identity.
  - `limits::Tuple{<:Number,<:Number}=(-Inf, Inf)`: Set it to a value different than `Inf` if you want to fix the limits of the binning.

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
        throw(ArgumentError("scaledBins: `limits` must be (min, max), but I got limits[1] = \
        $(limits[1]) > limits[2] = $(limits[2])"))
    )

    # Compute the limits of the binning
    min = isinf(limits[1]) ? scaling(ustrip(minimum(values))) : scaling(ustrip(limits[1]))
    max = isinf(limits[2]) ? scaling(ustrip(maximum(values))) : scaling(ustrip(limits[2]))

    # For a small range, increase it by 0.2 * abs(max)
    if (range = max - min) <= 1e-4 * abs(max)
        range += 0.2 * abs(max)
    end

    # Compute the width of the bins
    width = range / n_bins

    # Get the inverse function of `scaling`
    inverse = Makie.inverse_transform(scaling)

    # Get the unit of the data
    v_unit = unit(first(values))

    return [inverse(min + width * i) * v_unit for i in 0:n_bins]

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

  - `positions::Vector{<:Number}`: Positions of the `values` within a 1D axis. This determines to which bin each value will be added.
  - `values::Vector{<:Number}`: The values that will be added up in each bin, according to their `positions`.
  - `grid::Union{LinearGrid,CircularGrid}`: A linear or circular grid.
  - `total::Bool=true`: If the sum (default) or the mean of `values` will be computed for each bin.
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
        `positions`, but I got length(values) = $(length(values)) and \
        length(positions) = $(length(positions))"))
    )

    n_bins = length(grid.grid)

    if grid.log
        l_u = unit(first(grid.ticks))
        positions = log10.(ustrip.(l_u, positions))
        p_min = log10(ustrip(grid.ticks[1]))
        p_max = log10(ustrip(l_u, grid.ticks[end]))
    else
        p_min = grid.ticks[1]
        p_max = grid.ticks[end]
    end

    # Compute the bin width
    width = (p_max - p_min) / n_bins

    # Allocate memory
    histogram = zeros(eltype(values), n_bins)
    counts = zeros(Int, n_bins)

    # Compute the histogram; ignoring NaNs and positions outside the grid range
    @inbounds for (position, value) in zip(positions, values)

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

        @inbounds for (i, count) in pairs(counts)
            @inbounds if iszero(count)
                histogram[i] = nan
            end
        end
    end

    if !total
        # Compute the mean value instead of just the sum for each bin
        @inbounds for (i, count) in pairs(counts)
            @inbounds if !iszero(count)
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

  - `positions::Vector{<:Number}`: Positions of the `values` within a 1D axis. This determines to which bin each value will be added.
  - `values::Vector{<:Number}`: The values that will be added up in each bin, according to their `positions`.
  - `edges::Vector{<:Number}`: A sorted list of bin edges.
  - `total::Bool=true`: If the sum (default) or the mean of `values` will be computed for each bin.
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
        `positions`, but I got length(values) = $(length(values)) and \
        length(positions) = $(length(positions))"))
    )

    issorted(edges) || sort!(edges)

    n_bins = length(edges) - 1

    p_min = first(edges)
    p_max = last(edges)

    # Allocate memory
    histogram = zeros(eltype(values), n_bins)
    counts = zeros(Int, n_bins)

    # Compute the histogram; ignoring NaNs and positions outside of range
    @inbounds for (position, value) in zip(positions, values)

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

        @inbounds for (i, count) in pairs(counts)
            @inbounds if iszero(count)
                histogram[i] = nan
            end
        end
    end

    if !total
        # Compute the mean value instead of just the sum for each bin
        @inbounds for (i, count) in pairs(counts)
            @inbounds if !iszero(count)
                histogram[i] /= count
            end
        end
    end

    return histogram

end

"""
    histogram1D(
        positions::Vector{<:Number},
        grid::Union{LinearGrid,CircularGrid},
    )::Vector{Int}

Compute a 1D histogram of `positions`.

# Arguments

  - `positions::Vector{<:Number}`: Values for which the histogram will be constructed.
  - `grid::Union{LinearGrid,CircularGrid}`: A linear or circular grid.

# Returns

  - A vector with the counts.
"""
function histogram1D(
    positions::Vector{<:Number},
    grid::Union{LinearGrid,CircularGrid},
)::Vector{Int}

    n_bins = length(grid.grid)

    if grid.log
        l_u = unit(first(grid.ticks))
        positions = log10.(ustrip.(l_u, positions))
        p_min = log10(ustrip(grid.ticks[1]))
        p_max = log10(ustrip(l_u, grid.ticks[end]))
    else
        p_min = grid.ticks[1]
        p_max = grid.ticks[end]
    end

    # Compute the bin width
    width = (p_max - p_min) / n_bins

    # Allocate memory
    histogram = zeros(Int, n_bins)

    # Compute the histogram; ignoring NaNs and positions outside the grid range
    @inbounds for position in positions

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
  - `edges::Vector{<:Number}`: A sorted list of bin edges.

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

    # Compute the histogram; ignoring NaNs and positions outside of range
    @inbounds for position in positions

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

  - `positions::Matrix{<:Number}`: Positions of the values in the grid. Each column correspond to a value and each row is a dimension. This determines to which bin each value will be added.
  - `values::Vector{<:Number}`: The values that will be added up in each square bin, according to their `positions`.
  - `grid::SquareGrid`: A square grid.
  - `total::Bool=true`: If the sum (default) or the mean of `values` will be computed in each bin.
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
        has columns, but I got length(values) = $(length(values)) and size(positions, 2) = \
        $(size(positions, 2))"))
    )

    # Half bin size
    h_bin_width = grid.bin_width * 0.5

    # Compute the physical position of the grid borders
    h_borders = (grid.x_ticks[1] - h_bin_width, grid.x_ticks[end] + h_bin_width)
    v_borders = (grid.y_ticks[1] - h_bin_width, grid.y_ticks[end] + h_bin_width)

    # Allocate memory
    histogram = zeros(eltype(values), size(grid.grid))
    counts = zeros(Int, size(grid.grid))

    @inbounds for (i, point) in pairs(eachcol(positions))

        !isnan(values[i]) || continue

        x = point[1]
        y = point[2]

        !isnan(x) || continue
        !isnan(y) || continue

        if x > h_borders[2] || x < h_borders[1] || y > v_borders[2] || y < v_borders[1]
            continue
        end

        if x == h_borders[1]
            i_x = 1
        elseif x == h_borders[2]
            i_x = grid.n_bins
        else
            i_x = ceil(Int, (x - h_borders[1]) / grid.bin_width)
        end

        if y == v_borders[1]
            i_y = 1
        elseif y == v_borders[2]
            i_y = grid.n_bins
        else
            i_y = ceil(Int, (y - v_borders[1]) / grid.bin_width)
        end

        histogram[grid.n_bins - i_y + 1, i_x] += values[i]
        counts[grid.n_bins - i_y + 1, i_x] += 1

    end

    if empty_nan
        # Set empty bins to NaN
        nan = NaN * unit(first(values))

        @inbounds for i in eachindex(histogram)
            @inbounds if iszero(counts[i])
                histogram[i] = nan
            end
        end
    end

    if !total
        # Compute the mean value instead of just the sum for each bin
        @inbounds for i in eachindex(histogram)
            @inbounds if !iszero(counts[i])
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

  - `positions::Matrix{<:Number}`: Positions of the values in the grid. Each column correspond to a value and each row is a dimension. This determines to which bin each value will be added.
  - `values::Vector{<:Number}`: The values that will be added up in each square bin, according to their `positions`.
  - `x_edges::Vector{<:Number}`: A sorted list of bin edges for the x axis.
  - `y_edges::Vector{<:Number}`: A sorted list of bin edges for the y axis.
  - `total::Bool=true`: If the sum (default) or the mean of `values` will be computed in each bin.
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
        has columns, but I got length(values) = $(length(values)) and size(positions, 2) = \
        $(size(positions, 2))"))
    )

    issorted(x_edges) || sort!(x_edges)
    issorted(y_edges) || sort!(y_edges)

    n_x_bins = length(x_edges) - 1
    n_y_bins = length(y_edges) - 1

    h_borders = (first(x_edges), last(x_edges))
    v_borders = (first(v_edges), last(v_edges))

    # Allocate memory
    histogram = zeros(eltype(values), (n_x_bins, n_y_bins))
    counts = zeros(Int, (n_x_bins, n_y_bins))

    @inbounds for (i, point) in pairs(eachcol(positions))

        !isnan(values[i]) || continue

        x = point[1]
        y = point[2]

        !isnan(x) || continue
        !isnan(y) || continue

        if x > h_borders[2] || x < h_borders[1] || y > v_borders[2] || y < v_borders[1]
            continue
        end

        if x == h_borders[1]
            i_x = 1
        else
            i_x = searchsortedfirst(x_edges, x) - 1
        end

        if y == v_borders[1]
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

        @inbounds for i in eachindex(histogram)
            @inbounds if iszero(counts[i])
                histogram[i] = nan
            end
        end
    end

    if !total
        # Compute the mean value instead of just the sum for each bin
        @inbounds for i in eachindex(histogram)
            @inbounds if !iszero(counts[i])
                histogram[i] /= counts[i]
            end
        end
    end

    return histogram

end

"""
    histogram2D(positions::Matrix{<:Number}, grid::SquareGrid)::Matrix{Int}

Compute a 2D histogram of `positions`.

# Arguments

  - `positions::Matrix{<:Number}`: Values for which the histogram will be constructed.
  - `grid::SquareGrid`: A square grid.

# Returns

  - A matrix with the counts.
"""
function histogram2D(positions::Matrix{<:Number}, grid::SquareGrid)::Matrix{Int}

    # Half bin size
    h_bin_width = grid.bin_width * 0.5

    # Compute the physical position of the grid borders
    h_borders = (grid.x_ticks[1] - h_bin_width, grid.x_ticks[end] + h_bin_width)
    v_borders = (grid.y_ticks[1] - h_bin_width, grid.y_ticks[end] + h_bin_width)

    # Allocate memory
    histogram = zeros(Int, size(grid.grid))

    @inbounds for point in eachcol(positions)

        x = point[1]
        y = point[2]

        !isnan(x) || continue
        !isnan(y) || continue

        if x > h_borders[2] || x < h_borders[1] || y > v_borders[2] || y < v_borders[1]
            continue
        end

        if x == h_borders[1]
            i_x = 1
        elseif x == h_borders[2]
            i_x = grid.n_bins
        else
            i_x = ceil(Int, (x - h_borders[1]) / grid.bin_width)
        end

        if y == v_borders[1]
            i_y = 1
        elseif y == v_borders[2]
            i_y = grid.n_bins
        else
            i_y = ceil(Int, (y - v_borders[1]) / grid.bin_width)
        end

        histogram[grid.n_bins - i_y + 1, i_x] += 1

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
  - `grid::SquareGrid`: A square grid.
  - `x_edges::Vector{<:Number}`: A sorted list of bin edges for the x axis.
  - `y_edges::Vector{<:Number}`: A sorted list of bin edges for the y axis.

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

    h_borders = (first(x_edges), last(x_edges))
    v_borders = (first(y_edges), last(y_edges))

    # Allocate memory
    histogram = zeros(Int, (n_x_bins, n_y_bins))

    @inbounds for point in eachcol(positions)

        x = point[1]
        y = point[2]

        !isnan(x) || continue
        !isnan(y) || continue

        if x > h_borders[2] || x < h_borders[1] || y > v_borders[2] || y < v_borders[1]
            continue
        end

        if x == h_borders[1]
            i_x = 1
        else
            i_x = searchsortedfirst(x_edges, x) - 1
        end

        if y == v_borders[1]
            i_y = 1
        else
            i_y = searchsortedfirst(y_edges, y) - 1
        end

        histogram[n_y_bins - i_y + 1, i_x] += 1

    end

    return histogram

end

"""
    smoothWindow(
        x_data::Vector{<:Number},
        y_data::Vector{<:Number},
        n_bins::Int;
        <keyword arguments>
    )::NTuple{2,Vector{<:Number}}

Separate the values of `x_data` in `n_bins` bins and compute the mean value of `x_data` and `y_data` within each one.

# Arguments

  - `x_data::Vector{<:Number}`: x-axis data.
  - `y_data::Vector{<:Number}`: y-axis data.
  - `n_bins::Int`: Number of bins.
  - `scaling::Function=identity`: Scaling function for the x axis. The options are the scaling functions accepted by [Makie](https://docs.makie.org/stable/): log10, log2, log, sqrt, Makie.logit, Makie.Symlog10, Makie.pseudolog10, and identity. All the values of `x_data` must be in the domain of `scaling`.

# Returns

  - A tuple with two vectors, containing the smoothed-out x and y values.
"""
function smoothWindow(
    x_data::Vector{<:Number},
    y_data::Vector{<:Number},
    n_bins::Int;
    scaling::Function=identity,
)::NTuple{2,Vector{<:Number}}

    # Check that the input vectors have the same length
    (
        length(x_data) == length(y_data) ||
        throw(DimensionMismatch("smoothWindow: `x_data` and `y_data` must have the same length, \
        but I got length(x_data) = $(length(x_data)) != length(y_data) = $(length(y_data))"))
    )

    positions = scaling.(ustrip(x_data))
    grid = CircularGrid(maximum(positions), n_bins; shift=minimum(positions))

    smooth_x_data = histogram1D(positions, x_data, grid; total=false)
    smooth_y_data = histogram1D(positions, y_data, grid; total=false)

    # Remove empty bins
    return filter!(!isnan, smooth_x_data), filter!(!isnan, smooth_y_data)

end

####################################################################################################
# Makie utilities.
####################################################################################################

"""
Extract the limits of the x axis, from a [Makie](https://docs.makie.org/stable/) plot, axis, or figure. In the case of a figure, it will take the limits from the current axis object.
"""
function xlimits!(axis::Makie.Axis)::NTuple{2,Float32}
    reset_limits!(axis)
    return axis.xaxis.attributes.limits[]
end
xlimits!(plot::Makie.FigureAxisPlot)::NTuple{2,Float32} = xlimits!(plot.axis)
xlimits!(fig::Makie.Figure)::NTuple{2,Float32} = xlimits!(fig.current_axis.x)

"""
Extract the limits of the y axis, from a [Makie](https://docs.makie.org/stable/) plot, axis, or figure. In the case of a figure, it will take the limits from the current axis object.
"""
function ylimits!(axis::Makie.Axis)::NTuple{2,Float32}
    reset_limits!(axis)
    return axis.yaxis.attributes.limits[]
end
ylimits!(plot::Makie.FigureAxisPlot)::NTuple{2,Float32} = ylimits!(plot.axis)
ylimits!(fig::Makie.Figure)::NTuple{2,Float32} = ylimits!(fig.current_axis.x)

"""
Extract the scale function of the x axis, from a [Makie](https://docs.makie.org/stable/) plot, axis, or figure. In the case of a figure, it will take the scale from the current axis object.
"""
xscale(axis::Makie.Axis)::Function = axis.xaxis.attributes.scale[]
xscale(plot::Makie.FigureAxisPlot)::Function = xscale(plot.axis)
xscale(fig::Makie.Figure)::Function = xscale(fig.current_axis.x)

"""
Extract the scale function of the y axis, from a [Makie](https://docs.makie.org/stable/) plot, axis, or figure. In the case of a figure, it will take the scale from the current axis object.
"""
yscale(axis::Makie.Axis)::Function = axis.yaxis.attributes.scale[]
yscale(plot::Makie.FigureAxisPlot)::Function = yscale(plot.axis)
yscale(fig::Makie.Figure)::Function = yscale(fig.current_axis.x)

"""
Extract all the data points in a [Makie](https://docs.makie.org/stable/) plot, axis, or figure. In the case of a figure, it will only take the data from the current axis object. It only works for scatter, line and scatterline plots.
"""
function pointData(axis::Makie.Axis)::Vector{Point{2,Float32}}

    series = copy(axis.scene.plots)

    filter!(x -> isa(x, Union{Scatter,Lines,ScatterLines}), series)

    !isempty(series) || return Point{2,Float32}[]

    return collect(Iterators.flatten(serie[1][] for serie in series))

end
pointData(plot::Makie.FigureAxisPlot)::Vector{Point{2,Float32}} = pointData(plot.axis)
pointData(fig::Makie.Figure)::Vector{Point{2,Float32}} = pointData(fig.current_axis.x)

"""
    absCoor(
        plot::Union{Makie.FigureAxisPlot,Makie.Axis,Makie.Figure},
        r_x::Real,
        r_y::Real,
    )::NTuple{2,Float64}

Compute the absolute x and y coordinates of a plot, from the relative ones.

# Arguments

  - `plot::Union{Makie.FigureAxisPlot,Makie.Axis,Makie.Figure}`: Plot, axis, or figure for which the absolute coordinates will be calculated. In the case of a figure, it will use the limits from the current axis object.
  - `r_x::Real`: Relative x coordinate, `rx` ∈ [0, 1].
  - `r_y::Real`: Relative y coordinate, `ry` ∈ [0, 1].

# Returns

  - A tuple with the absolute coordinates, (a_x, a_y).

# Examples

```julia-repl
julia> absCoor(lines(rand(100)), 0.5, 0.5)
(50.50000071525574, 0.48792968317866325)
```
"""
function absCoor(
    plot::Union{Makie.FigureAxisPlot,Makie.Axis,Makie.Figure},
    r_x::Real,
    r_y::Real,
)::NTuple{2,Float64}

    # Get the scaling functions
    x_scale = xscale(plot)
    y_scale = yscale(plot)

    # Get the limits of the axes
    x_limits = x_scale.(xlimits!(plot))
    y_limits = y_scale.(ylimits!(plot))

    # Compute the absolute coordinates
    a_x = Makie.inverse_transform(x_scale).(x_limits[1] + r_x * (x_limits[2] - x_limits[1]))
    a_y = Makie.inverse_transform(y_scale).(y_limits[1] + r_y * (y_limits[2] - y_limits[1]))

    return a_x, a_y

end

"""
    cleanPlot!(figure::Makie.Figure)::Nothing

Delete all the legends of a figure and empty all its axes.

# Arguments

  - `figure::Makie.Figure`: Figure to be cleaned.
"""
function cleanPlot!(figure::Makie.Figure)::Nothing

    # Compute the number of elements (axes and legends) in the figure
    n_elements = length(figure.content)
    i = 1

    for _ in 1:n_elements
        i += cleanPlot!(figure.content[i])
    end

    return nothing

end

"""
    cleanPlot!(ax::Makie.Axis)::Bool

Empty an axis.

# Arguments

  - `ax::Makie.Axis`: Axis to be emptied.

# Returns

  - Flag to indicate that an axis has been emptied.
"""
function cleanPlot!(ax::Makie.Axis)::Bool

    empty!(ax)

    return true

end

"""
    cleanPlot!(legend::Makie.Legend)::Bool

Delete a legend.

# Arguments

  - `legend::Makie.Legend`: Legend to be deleted.

# Returns

  - Flag to indicate that a legend has been deleted.
"""
function cleanPlot!(legend::Makie.Legend)::Bool

    delete!(legend)

    return false

end

"""
Default function to end `cleanPlot!` recursion if an unknown type is encountered.
"""
cleanPlot!(default) = error("cleanPlot!: I cannot clean elements of type $(typeof(default))")

"""
    cubicSplineKernel(q::Real, h::Number)::Number

2D cubic spline kernel.

# Arguments

  - `q::Real`: Relative distance to the neighbor, ``|r - r'| / h``.
  - `h::Number`: Smoothing length.

# Returns

  - The kernel function evaluated at a separation `q` * `h`, and with a smoothing length `h`.

# References

[PySPH documentation](https://pysph.readthedocs.io/en/latest/reference/kernels.html)

J. J. Monaghan (1992). *Smoothed Particle Hydrodynamics*. Annual Review of Astronomy and Astrophysics, **30**, 543-574. [doi:10.1146/annurev.aa.30.090192.002551](https://doi.org/10.1146/annurev.aa.30.090192.002551)

M.B. Liu et al. (2010). *Smoothed Particle Hydrodynamics (SPH): an Overview and Recent Developments*. Archives of Computational Methods in Engineering, **17**, 25–76. [doi:10.1007/s11831-010-9040-7](https://doi.org/10.1007/s11831-010-9040-7)
"""
function cubicSplineKernel(q::Real, h::Number)::Number

    (
        isPositive(q, h) ||
        throw(DomainError("cubicSplineKernel: `q` and `h` must be positive, \
        but I got `q` = $q and `h` = $h"))
    )

    σ3 = 10.0 / (7.0 * area(h))

    if 0 <= q <= 1
        return σ3 * (1.0 - 1.5 * q * q * (1.0 - q * 0.5))
    elseif 1 < q <= 2
        return (σ3 * 0.25) * (2.0 - q)^3.0
    else
        return zero(σ3)
    end

end
