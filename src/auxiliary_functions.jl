####################################################################################################
# General utilities.
####################################################################################################

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
    parserWS(data::AbstractString)::Union{Float64,Missing}

Parse a string as a Float64, ignoring white spaces. If the string is empty return missing.

# Arguments

  - `data::AbstractString`: String to be parsed.

# Returns

  - Number in the string as a Float64.
"""
function parserWS(data::AbstractString)::Union{Float64,Missing}

	clean_data = strip(data)

	!isempty(clean_data) || return missing

	return parse(Float64, clean_data)

end

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
setPositive(x::Number) = x >= zero(x) ? x : zero(x)

"""
Test for strict positivity.
"""
isPositive(x::Number)::Bool = x > zero(x)
isPositive(x::AbstractArray)::Bool = all(isPositive, x)
isPositive(x...)::Bool = all(isPositive, x)

"""
New method for `Base.iszero` to compare [`IndexType`](@ref) with 0 as an interger.
"""
Base.iszero(x::IndexType)::Bool = x == 0

"""
New method for `Base.isempty` to check for empty [LaTeXStrings](https://github.com/JuliaStrings/LaTeXStrings.jl).
"""
Base.isempty(l_str::LaTeXString)::Bool = l_str == L""

"""
New methods for `Base.intersect` to use with the `Colon` type.
"""
Base.intersect(a1::Colon, a2::IndexType)::IndexType = a2
Base.intersect(a1::IndexType, a2::Colon)::IndexType = a1
Base.intersect(a1::Colon, a2::Colon)::Colon = (:)

"""
New methods for `Base.intersect` to use with the `Vector{Bool}` type.
"""
Base.intersect(a1::Vector{Bool}, a2::ReducedIndexType)::Vector{Int} = findall(a1) ∩ a2
Base.intersect(a1::ReducedIndexType, a2::Vector{Bool})::Vector{Int} = a1 ∩ findall(a2)
Base.intersect(a1::Vector{Bool}, a2::Vector{Bool})::Vector{Bool} = Vector{Bool}(a1 .&& a2)

"""
New methods for `Base.union` to use with the `Vector{Bool}` type.
"""
Base.union(a1::Vector{Bool}, a2::ReducedIndexType)::Vector{Int} = findall(a1) ∪ a2
Base.union(a1::ReducedIndexType, a2::Vector{Bool})::Vector{Int} = a1 ∪ findall(a2)
Base.union(a1::Vector{Bool}, a2::Vector{Bool})::Vector{Bool} = Vector{Bool}(a1 .|| a2)

"""
Area of a circle with radius `r`.
"""
function area(r::Number)::Number

    x = setPositive(r)

    return π * x * x

end

"""
Volume of a sphere with radius `r`.
"""
function volume(r::Number)::Number

    x = setPositive(r)

    return π * x * x * x * 1.333

end

"""
    evaluateNormal(data::Vector{<:Number})::Vector{<:Number}

Evaluate a normal distribution at the values in `data`.

The mean and standard deviation of the distribution are the ones from the `data` itself.

# Arguments

  - `data::Vector{<:Number}`: Data vector used to compute the mean and standard deviation of the normal distribution.

# Returns

  - The normal distribution evaluated at the values in `data`.
"""
function evaluateNormal(data::Vector{<:Number})::Vector{<:Number}

    μ  = mean(data)
    σ2 = 2.0 * std(data; mean=μ)^2
    d2 = @. (data - μ)^2

    return @. (1.0 / sqrt(σ2 * π)) * exp(-d2 / σ2)

end

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
  - `output_path::String="./joined_images.png"`: Path to the output image.
"""
function hvcatImages(
    blocks_per_row::Int,
    paths::Vector{String};
    output_path::String="./joined_images.png",
)::Nothing

    new_image = hvcat(blocks_per_row, [load(path) for path in paths]...)

    save(output_path, new_image)

    return nothing

end

"""
    rangeCut!(
        raw_values::Vector{<:Number},
        range::Tuple{<:Number,<:Number};
        <keyword arguments>
    )::Bool

Delete every element in `raw_values` that is outside the given `range`.

# Arguments

  - `raw_values::Vector{<:Number}`: Dataset that will be pruned.
  - `range::Tuple{<:Number,<:Number}`: The range in question.
  - `keep_edges::Bool=true`: If the edges of the range will be kept.
  - `min_left::Int=1`: Minimum number of values that need to be left after pruning to proceed with the transformation.

# Returns

  - If a transformation was performed.
"""
function rangeCut!(
    raw_values::Vector{<:Number},
    range::Tuple{<:Number,<:Number};
    keep_edges::Bool=true,
    min_left::Int=1,
)::Bool

    # Shortcut computation for special cases
    !(isempty(raw_values) || all(isinf.(range))) || return false

    if keep_edges

        # Check that after the transformation at least `min_left` elements will be left
        count(x -> range[1] <= x <= range[2], raw_values) >= min_left || return false

        # Delete element outside of the provided range
        filter!(x -> range[1] <= x <= range[2], raw_values)

    else

        # Check that after the transformation at least `min_left` elements will be left
        count(x -> range[1] < x < range[2], raw_values) >= min_left || return false

        # Delete element outside of the provided range
        filter!(x -> range[1] < x < range[2], raw_values)

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
  - `min_left::Int=1`: Minimum number of values that need to be left in the master dataset after pruning to proceed with the transformation.

# Returns

  - If a transformation was performed.
"""
function rangeCut!(
    m_data::Vector{<:Number},
    s_data::Vector,
    range::Tuple{<:Number,<:Number};
    keep_edges::Bool=true,
    min_left::Int=1,
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
        raw_values::Vector{<:Number};
        <keyword arguments>
    )::NTuple{2,Bool}

Do the following transformations over `raw_values`, in order:

  - Trim it to fit within the domain of the function `func_domain`.
  - Trim it to fit within `range`.
  - Scale it down by a factor of 10^`exp_factor`.

By default, no transformation is done.

# Arguments

  - `raw_values::Vector{<:Number}`: Dataset to be sanitized.
  - `func_domain::Function=identity`: `raw_values` will be trimmed to fit within the domain of the function `func_domain`. The options are the scaling functions accepted by [Makie](https://docs.makie.org/stable/): log10, log2, log, sqrt, Makie.logit, Makie.Symlog10, Makie.pseudolog10, and identity.
  - `range::Tuple{<:Number,<:Number}=(-Inf, Inf)`: Every element in `raw_values` that falls outside of `range` will be deleted.
  - `keep_edges::Bool=true`: If the edges of `range` will be kept.
  - `min_left::Int=1`: Minimum number of values that need to be left after each transformation to procced with it.
  - `exp_factor::Int=0`: Every element in `raw_values` will be divided by 10^`exp_factor`.
  - `warnings::Bool=true`: If a warning will be given when `raw_values` is a vector of Integers, which may cause wrong results when dividing by 10^`exp_factor`.

# Returns

  - A tuple with two flags:

      + If `raw_values` was mutated to fit within the domain of `func_domain`.
      + If `raw_values` was mutated to fit within `range`.
"""
function sanitizeData!(
    raw_values::Vector{<:Number};
    func_domain::Function=identity,
    range::Tuple{<:Number,<:Number}=(-Inf, Inf),
    keep_edges::Bool=true,
    min_left::Int=1,
    exp_factor::Int=0,
    warnings::Bool=true,
)::NTuple{2,Bool}

    !isempty(raw_values) || return false, false

    d_unit = unit(first(raw_values))

    # Trim `raw_values` to fit within the domain of `func_domain`
    if func_domain ∈ [identity, Makie.pseudolog10, Makie.Symlog10]

        domain_flag = false

    elseif func_domain == sqrt

        domain_flag = rangeCut!(raw_values, (0.0, Inf) .* d_unit; keep_edges=true, min_left)

    elseif func_domain == Makie.logit

        domain_flag = rangeCut!(raw_values, (0.0, 1.0) .* d_unit; keep_edges=false, min_left)

    elseif func_domain ∈ [log, log2, log10]

        domain_flag = rangeCut!(raw_values, (0.0, Inf) .* d_unit; keep_edges=false, min_left)

    else

        throw(ArgumentError("sanitizeData!: The function $(func_domain) is not supported. See \
        the list of supported scaling functions in the [Makie](https://docs.makie.org/stable/) \
        documentation"))

    end

    # Trim `raw_values` to fit within `range`
    range_flag = rangeCut!(raw_values, range; keep_edges, min_left)

    (
        !(isa(raw_values, Vector{<:Integer}) && !iszero(exp_factor) && warnings) ||
        @warn("sanitizeData!: Elements of `raw_values` are of type `Integer`, this may result \
        in errors or unwanted truncation when using `exp_factor` != 0")
    )

    # Scale `raw_values` down by a factor of 10^`exp_factor`
    iszero(exp_factor) || (raw_values ./= exp10(exp_factor))

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
  - `min_left::Int=1`: Minimum number of values that need to be left in each dataset after any of the transformations to procced with them.
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
    min_left::Int=1,
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
        throw(ArgumentError("scaledBins: `limits` must be (min, max), but I got limits[1] = \
        $(limits[1]) > limits[2] = $(limits[2])"))
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
    listHistogram1D(
        positions::Vector{<:Number},
        values::Vector{<:Number},
        grid::Union{LinearGrid,CircularGrid},
    )::Vector{Vector{<:Number}}

Compute a 1D histogram of `values`, returning the full list of `values` within each bin.

# Arguments

  - `positions::Vector{<:Number}`: Positions of the `values` within a 1D axis. This determines to which bin each value will be added.
  - `values::Vector{<:Number}`: The values that will be added up in each bin, according to their `positions`.
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
    histogram = [eltype(values)[] for _ in 1:n_bins]

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

        push!(histogram[idx], value)

    end

    return histogram

end

"""
    listHistogram1D(
        positions::Vector{<:Number},
        values::Vector{<:Number},
        edges::Vector{<:Number},
    )::Vector{Vector{<:Number}}

Compute a 1D histogram of `values`, returning the full list of `values` within each bin.

# Arguments

  - `positions::Vector{<:Number}`: Positions of the `values` within a 1D axis. This determines to which bin each value will be added.
  - `values::Vector{<:Number}`: The values that will be added up in each bin, according to their `positions`.
  - `edges::Vector{<:Number}`: A sorted list of bin edges.

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
        throw(ArgumentError("histogram1D: `values` must have as many elements as \
        `positions`, but I got length(values) = $(length(values)) and \
        length(positions) = $(length(positions))"))
    )

    issorted(edges) || sort!(edges)

    n_bins = length(edges) - 1

    p_min = first(edges)
    p_max = last(edges)

    # Allocate memory
    histogram = [eltype(values)[] for _ in 1:n_bins]

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

        push!(histogram[idx], value)

    end

    return histogram

end

"""
    listHistogram3D(positions::Matrix{<:Number}, grid::CubicGrid)::Array{Vector{Int},3}

Compute a 3D histogram of `positions`, returning the full list of indices within each bin.

# Arguments

  - `positions::Matrix{<:Number}`: Positions of the points in the grid. Each column correspond to a point and each row is a dimension. This determines to which bin the index of each point will be added.
  - `grid::CubicGrid`: A cubic grid.

# Returns

    - A tensor with the indices of the points within each bin.
"""
function listHistogram3D(positions::Matrix{<:Number}, grid::CubicGrid)::Array{Vector{Int},3}

    # Half bin size
    h_bin_width = grid.bin_width * 0.5

    # Compute the physical position of the grid borders
    x_borders = (grid.x_ticks[1] - h_bin_width, grid.x_ticks[end] + h_bin_width)
    y_borders = (grid.y_ticks[1] - h_bin_width, grid.y_ticks[end] + h_bin_width)
    z_borders = (grid.z_ticks[1] - h_bin_width, grid.z_ticks[end] + h_bin_width)

    # Allocate memory
    histogram = Array{Vector{Int}}(undef, size(grid.grid))
    @inbounds for i in eachindex(histogram)
        histogram[i] = Int[]
    end

    @inbounds for (idx, point) in pairs(eachcol(positions))

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
    x_borders = (grid.x_ticks[1] - h_bin_width, grid.x_ticks[end] + h_bin_width)
    y_borders = (grid.y_ticks[1] - h_bin_width, grid.y_ticks[end] + h_bin_width)

    # Allocate memory
    histogram = zeros(eltype(values), size(grid.grid))
    counts = zeros(Int, size(grid.grid))

    @inbounds for (i, point) in pairs(eachcol(positions))

        !isnan(values[i])  || continue
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

    x_borders = (first(x_edges), last(x_edges))
    y_borders = (first(y_edges), last(y_edges))

    # Allocate memory
    histogram = zeros(eltype(values), (n_x_bins, n_y_bins))
    counts = zeros(Int, (n_x_bins, n_y_bins))

    @inbounds for (i, point) in pairs(eachcol(positions))

        !isnan(values[i])  || continue
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
    x_borders = (grid.x_ticks[1] - h_bin_width, grid.x_ticks[end] + h_bin_width)
    y_borders = (grid.y_ticks[1] - h_bin_width, grid.y_ticks[end] + h_bin_width)

    # Allocate memory
    histogram = zeros(Int, size(grid.grid))

    @inbounds for point in eachcol(positions)

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

    x_borders = (first(x_edges), last(x_edges))
    y_borders = (first(y_edges), last(y_edges))

    # Allocate memory
    histogram = zeros(Int, (n_x_bins, n_y_bins))

    @inbounds for point in eachcol(positions)

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
    histogram3D(
        positions::Matrix{<:Number},
        values::Vector{<:Number},
        grid::CubicGrid;
        <keyword arguments>
    )::Array{<:Number,3}

Compute a 3D histogram of `values`.

# Arguments

  - `positions::Matrix{<:Number}`: Positions of the values in the grid. Each column correspond to a value and each row is a dimension. This determines to which bin each value will be added.
  - `values::Vector{<:Number}`: The values that will be added up in each square bin, according to their `positions`.
  - `grid::CubicGrid`: A cubic grid.
  - `total::Bool=true`: If the sum (default) or the mean of `values` will be computed in each bin.
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
        has columns, but I got length(values) = $(length(values)) and size(positions, 2) = \
        $(size(positions, 2))"))
    )

    # Half bin size
    h_bin_width = grid.bin_width * 0.5

    # Compute the physical position of the grid borders
    x_borders = (grid.x_ticks[1] - h_bin_width, grid.x_ticks[end] + h_bin_width)
    y_borders = (grid.y_ticks[1] - h_bin_width, grid.y_ticks[end] + h_bin_width)
    z_borders = (grid.z_ticks[1] - h_bin_width, grid.z_ticks[end] + h_bin_width)

    # Allocate memory
    histogram = zeros(eltype(values), size(grid.grid))
    counts = zeros(Int, size(grid.grid))

    @inbounds for (i, point) in pairs(eachcol(positions))

        !isnan(values[i])  || continue
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
    histogram3D(positions::Matrix{<:Number}, grid::CubicGrid)::Array{Int,3}

Compute a 3D histogram of `positions`.

# Arguments

  - `positions::Matrix{<:Number}`: Values for which the histogram will be constructed.
  - `grid::CubicGrid`: A cubic grid.

# Returns

  - A 3D tensor with the counts.
"""
function histogram3D(positions::Matrix{<:Number}, grid::CubicGrid)::Array{Int,3}

    # Half bin size
    h_bin_width = grid.bin_width * 0.5

    # Compute the physical position of the grid borders
    x_borders = (grid.x_ticks[1] - h_bin_width, grid.x_ticks[end] + h_bin_width)
    y_borders = (grid.y_ticks[1] - h_bin_width, grid.y_ticks[end] + h_bin_width)
    z_borders = (grid.z_ticks[1] - h_bin_width, grid.z_ticks[end] + h_bin_width)

    # Allocate memory
    histogram = zeros(Int, size(grid.grid))

    @inbounds for point in eachcol(positions)

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

        histogram[grid.n_bins - i_y + 1, i_x, i_z] += 1

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
        throw(ArgumentError("deltas: `data` must have at least 2 elements, but it has only $(nd)"))
    )

    # Set the first value to 0, by convention
    Δd[1] = zero(data[1])

    @inbounds for i in 2:nd
        Δd[i] = data[i] - data[i - 1]
    end

    return Δd

end

"""
    reduceResolution(hr_matrix::Matrix{<:Number}, factor::Int)::Matrix{<:Number}

Reduce the number of rows and columns of `hr_matrix` by `factor`, averaging its values.

# Arguments

  - `hr_matrix::Matrix{<:Number}`: Original "high resolution" matrix. It has to be a square matrix.
  - `factor::Int`: Factor by wich the number of rows and columns will be reduced. It has to divede the size of `hr_matrix` exactly.

# Returns

  - The new smaller matrix, with the average values.
"""
function reduceResolution(hr_matrix::Matrix{<:Number}, factor::Int)::Matrix{<:Number}

    !isone(factor) || return hr_matrix

    r, c = size(hr_matrix)
    (
        r == c ||
        throw(ArgumentError("reduceResolution: `hr_matrix` has to be a square matrix, but it has \
        $(c) columns and $(r) rows."))
    )

    (
        factor >= 1 || throw(ArgumentError("reduceResolution: `factor` must be >= 1, \
        but I got `factor` = $(factor)."))
    )

    (
        r % factor == 0 ||
        throw(ArgumentError("reduceResolution: `factor` must divide the size of `hr_matrix` \
        exactly, but I got number of rows / `factor` = $(r / factor)."))
    )

    # Compute the size of the new matrix
    new_size = r ÷ factor

    # Allocate memory
    lr_matrix = zeros(eltype(hr_matrix), new_size, new_size)

    # Compute the number of old pixels per new pixel
    old_n_pixels = factor * factor

    @inbounds for i in eachindex(lr_matrix)

        # Compute the row and column of the new matrix corresponding to index i
        row = mod1(i, new_size)
        col = ceil(Int, i / new_size)

        @inbounds for j in (factor * (row - 1) + 1):(factor * row)
            @inbounds for k in (factor * (col - 1) + 1):(factor * col)

                if !isnan(hr_matrix[j, k])
                    lr_matrix[i] += hr_matrix[j, k]
                end

            end
        end

    end

    return lr_matrix ./ old_n_pixels

end

"""
    reduceTicks(hr_ticks::Vector{<:Number}, factor::Int)::Vector{<:Number}

Reduce the number of ticks in `hr_ticks` by `factor` keeping the total length of the axis the same and assuming `hr_ticks` are regularly spaced.

# Arguments

  - `hr_ticks::Vector{<:Number}`: Original "high resolution" list of ticks.
  - `factor::Int`: Factor by wich the number of ticks will be reduced. It has to divede the size of `hr_ticks` exactly.

# Returns

  - The new shorter tick list.
"""
function reduceTicks(hr_ticks::Vector{<:Number}, factor::Int)::Vector{<:Number}

    !isone(factor) || hr_ticks

    l = length(hr_ticks)
    (
        l % factor == 0 ||
        throw(ArgumentError("reduceTicks: `factor` must divide the size of `hr_ticks` \
        exactly, but I got length(`hr_ticks`) / `factor` = $(l / factor)."))
    )

    (
        factor >= 1 || throw(ArgumentError("reduceTicks: `factor` must be >= 1, \
        but I got `factor` = $(factor)."))
    )

    # Compute the size of the new vector
    new_size = l ÷ factor

    # Allocate memory
    lr_ticks = similar(hr_ticks, new_size)

    if iseven(factor)

        shift = factor ÷ 2

        @inbounds for i in eachindex(lr_ticks)

            idx = (i - 1) * factor + shift

            lr_ticks[i] = (hr_ticks[idx] + hr_ticks[idx + 1]) / 2.0

        end

    else

        shift = ceil(Int, factor / 2)

        @inbounds for i in eachindex(lr_ticks)

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

# Returns

  - A vector with the averages of the values in each concentric ring.
"""
function projectIntoCircularGrid(
    image::Matrix{<:Number},
    n_bins::Int;
    inscribed::Bool=true,
)::Vector{<:Number}

    r, c = size(image)
    (
        r == c ||
        throw(ArgumentError("projectIntoCircularGrid: `image` has to be a square matrix, \
        but it has $(c) columns and $(r) rows."))
    )

    # Construct a square grid center in (0, 0)
    square_grid = SquareGrid(1.0, r)

    # Construct a circular grid center in (0, 0)
    circular_grid = CircularGrid(inscribed ? 0.5 : sqrt(0.5), n_bins)

    # Compute the radial distance to the origin of each pixel in the square grid
    positions = norm.(vec(square_grid.grid))

    profile = histogram1D(
        positions,
        vec(image),
        circular_grid;
        total=false,
        empty_nan=false,
    )

    return profile

end

####################################################################################################
# Makie.jl utilities.
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
  - `r_x::Real`: Relative x coordinate.
  - `r_y::Real`: Relative y coordinate.

# Returns

  - A tuple with the absolute coordinates, (x, y).

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
    a_x = Makie.inverse_transform(x_scale).(x_limits[1] + abs(r_x) * (x_limits[2] - x_limits[1]))
    a_y = Makie.inverse_transform(y_scale).(y_limits[1] + abs(r_y) * (y_limits[2] - y_limits[1]))

    return sign(r_x) * a_x, sign(r_y) * a_y

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
    cleanPlot!(legend::Union{Makie.Legend,Makie.Colorbar})::Bool

Delete a legend or colorbar.

# Arguments

  - `legend::Union{Makie.Legend,Makie.Colorbar}`: Legend or colorbar to be deleted.

# Returns

  - Flag to indicate that a legend or colorbar has been deleted.
"""
function cleanPlot!(legend::Union{Makie.Legend,Makie.Colorbar})::Bool

    delete!(legend)

    return false

end

"""
Default function to end `cleanPlot!` recursion if an unknown type is encountered.
"""
cleanPlot!(default) = error("cleanPlot!: I cannot clean elements of type $(typeof(default))")
