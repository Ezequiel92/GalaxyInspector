####################################################################################################
# Auxiliary functions
####################################################################################################

"""
Return a copy of `list` with every negative value set to 0.
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
New method for `Base.iszero` to compare [`IndexType`](@ref) with 0.
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

    return π * x^2

end

"""
Volume of a sphere with radius `r`.
"""
function volume(r::Number)::Number

    x = setPositive(r)

    return π * x^3 * 4.0 / 3.0

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
    metaFormatter(level::LogLevel, _module, group, id, file, line)

Formatter for loggers.

See the documentation of [ConsoleLogger](https://docs.julialang.org/en/v1/stdlib/Logging/#Base.CoreLogging.ConsoleLogger)
"""
function metaFormatter(level::LogLevel, _module, group, id, file, line)

    @nospecialize

    color = Logging.default_logcolor(level)
    prefix = string(level == Warn ? "Warning" : string(level), " |")
    suffix::String = ""

    if file !== nothing
        suffix *= contractuser(file)::String
        if line !== nothing
            suffix *= ":$(isa(line, UnitRange) ? "$(first(line))-$(last(line))" : line)"
        end
    end

    if !isempty(suffix)
        suffix = "@ " * suffix * "\n"
    else
        suffix *= "\n"
    end

    return color, prefix, suffix

end

"""
    setLogging!(log::Bool; <keyword arguments>)::Nothing

Set if logging messages will be printed out. By default no logs are printed.

# Arguments

  - `log::Bool`: If logs will be printed out using the default logger.
  - `stream::IO=stdout`: Where to print the logs. It can be a file.
"""
function setLogging!(log::Bool; stream::IO=stdout)::Nothing

    logging[] = log

    log && global_logger(ConsoleLogger(stream; meta_formatter=metaFormatter))

    return nothing

end

"""
    ring(vec::Vector, index::Integer)::Vector

Make the indexing operation `vec[index]` work with modular arithmetic for the indices.

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
    parserWS(value::AbstractString)::Union{Float64,Missing}

Parse a string as a Float64, ignoring white spaces. If the string is empty return missing.

# Arguments

  - `value::AbstractString`: String to be parsed.

# Returns

  - Number in the string as a Float64.
"""
function parserWS(value::AbstractString)::Union{Float64,Missing}

    clean_value = strip(value)

    !isempty(clean_value) || return missing

    return parse(Float64, clean_value)

end

"""
    safeSelect(vec::Vector, index::IndexType)

Do the indexing operation `vec[index]` while ignoring indices that are out of bounds.

# Arguments

  - `vec::Vector`: Vector.
  - `index::IndexType`: Indices. It can be an integer (a single element), a vector of integers (several elements), an `UnitRange` (e.g. 5:13), an `StepRange` (e.g. 5:2:13) or (:) (every element).

# Returns

  - `vec[index (minus out of bounds indices)]`

# Examples

```julia-repl
julia> safeSelect([1, 2, 3], 11)
Int[]

julia> safeSelect([1, 2, 3], 1:5)
3-element Vector{Int}:
 1
 2
 3

julia> safeSelect([1, 2, 3], 1:3:10)
1

julia> safeSelect([1, 2, 3], [1, 2, 5, 9])
2-element Vector{Int}:
 1
 2

julia> safeSelect([1, 2, 3], (:))
3-element Vector{Int}:
 1
 2
 3
```
"""
function safeSelect(vec::Vector, index::IndexType)

    index != (:) || return vec

    index_list = [index...]

    filter!(x -> x <= length(vec), index_list)

    (
        length(index_list) == length([index...]) ||
        !logging[] ||
        @info("safeSelect: There are out of bounds indices")
    )

    if length(index_list) == 1
        return vec[index_list...]
    else
        return vec[index_list]
    end

end

"""
    evaluateNormal(values::Vector{<:Number})::Vector{<:Number}

Evaluate a normal distribution at `values`.

The median and standard deviation for the distribution are the ones of `values`.

# Arguments

  - `values::Vector{<:Number}`: Vector used to compute the mean and standard deviation of the normal distribution.

# Returns

  - The normal distribution evaluated at `values`.
"""
function evaluateNormal(values::Vector{<:Number})::Vector{<:Number}

    μ  = mean(values)
    σ2 = 2.0 * std(values; mean=μ)^2
    d2 = @. (values - μ)^2

    return @. (1.0 / sqrt(σ2 * π)) * exp(-d2 / σ2)

end

"""
    hvcatImages(
        blocks_per_row::Int,
        paths::Vector{String};
        <keyword arguments>
    )::Nothing

Join several images vertically and horizontally.

The images in `paths` will fill the rows and columns starting at the top left, going from left to right and from top to bottom (row-major order).

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

    isempty(paths) && throw(ArgumentError("hvcatImages: `paths` is empty"))

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

        # Delete elements outside of the provided range
        filter!(x -> range[1] <= x <= range[2], raw_values)

    else

        # Check that after the transformation at least `min_left` elements will be left
        count(x -> range[1] < x < range[2], raw_values) >= min_left || return false

        # Delete elements outside of the provided range
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
        but I got length(`s_data`) = $(length(s_data)) < length(`m_data`) = $(length(m_data))"))
    )

    if keep_edges

        # Find the elements outside of the provided range
        idxs = map(x -> x < range[1] || x > range[2], m_data)

        # Check that after the transformation at least `min_left` elements will be left
        count(.!idxs) >= min_left || return false

        # Delete elements outside of the provided range
        deleteat!(m_data, idxs)
        deleteat!(s_data, idxs)

    else

        # Find the elements outside of the provided range
        idxs = map(x -> x <= range[1] || x >= range[2], m_data)

        # Check that after the transformation at least `min_left` elements will be left
        count(.!idxs) >= min_left || return false

        # Delete elements outside of the provided range
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
  - `min_left::Int=1`: Minimum number of values that need to be left after each transformation to proceed with it.
  - `exp_factor::Int=0`: Every element in `raw_values` will be divided by 10^`exp_factor`.

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
        !(isa(raw_values, Vector{<:Integer}) && !iszero(exp_factor) && logging[]) ||
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

    The datasets must have the same length, and any operation that deletes an element, will delete the corresponding element (i.e. with the same index) in the other dataset, so that the datasets will remain of equal length.

# Arguments

  - `x_data::Vector{<:Number}`: First dataset to be sanitized.
  - `y_data::Vector{<:Number}`: Second dataset to be sanitized.
  - `func_domain::NTuple{2,Function}=(identity, identity)`: `x_data` will be trimmed to fit within the domain of the function `func_domain[1]`, and `y_data` will be trimmed to fit within the domain of the function `func_domain[2]`. The options are the scaling functions accepted by [Makie](https://docs.makie.org/stable/): log10, log2, log, sqrt, Makie.logit, Makie.Symlog10, Makie.pseudolog10, and identity.
  - `range::Tuple{Tuple{<:Number,<:Number},Tuple{<:Number,<:Number}}=((-Inf, Inf), (-Inf, Inf))`: Every element in `x_data` that falls outside of `range[1]` will be deleted, and every element in `y_data` that falls outside of `range[2]` will be deleted.
  - `keep_edges::NTuple{2,Bool}=(true, true)`: If the edges of each corresponding `range` will be kept.
  - `min_left::Int=1`: Minimum number of values that need to be left in each dataset after any of the transformations to proceed with them.
  - `exp_factor::NTuple{2,Int}=(0, 0)`: Every element in `x_data` will be divided by 10^`exp_factor[1]`, and every element in `y_data` will be divided by 10^`exp_factor[2]`.

# Returns

  - A tuple with four flags:

      + If `x_data` was mutated to fit within the domain of `func_domain[1]`.
      + If `y_data` was mutated to fit within the domain of `func_domain[2]`.
      + If `x_data` was mutated to fit within `range[1]`.
      + If `y_data` was mutated to fit within `range[2]`.
"""
function sanitizeData!(
    x_data::Vector{<:Number},
    y_data::Vector{<:Number};
    func_domain::NTuple{2,Function}=(identity, identity),
    range::Tuple{Tuple{<:Number,<:Number},Tuple{<:Number,<:Number}}=((-Inf, Inf), (-Inf, Inf)),
    keep_edges::NTuple{2,Bool}=(true, true),
    min_left::Int=1,
    exp_factor::NTuple{2,Int}=(0, 0),
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
        !(isa(x_data, Vector{<:Integer}) && !iszero(exp_factor[1]) && logging[]) ||
        @warn("sanitizeData!: Elements of `x_data` are of type Integer, this may result \
        in errors or unwanted truncation when using `exp_factor[1]` != 0")
    )

    (
        !(isa(y_data, Vector{<:Integer}) && !iszero(exp_factor[2]) && logging[]) ||
        @warn("sanitizeData!: Elements of `y_data` are of type Integer, this may result \
        in errors or unwanted truncation when using `exp_factor[2]` != 0")
    )

    # Scale the data down by the factors `exp_factor`
    iszero(exp_factor[1]) || (x_data ./= exp10(exp_factor[1]))
    iszero(exp_factor[2]) || (y_data ./= exp10(exp_factor[2]))

    return x_domain_flag, y_domain_flag, x_range_flag, y_range_flag

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
        but I got length(`x_data`) = $(length(x_data)) != length(`y_data`) = $(length(y_data))"))
    )

    positions = scaling.(ustrip(x_data))
    grid = CircularGrid(maximum(positions), n_bins; shift=minimum(positions))

    smooth_x_data = histogram1D(positions, x_data, grid; total=false)
    smooth_y_data = histogram1D(positions, y_data, grid; total=false)

    # Remove empty bins
    return filter!(!isnan, smooth_x_data), filter!(!isnan, smooth_y_data)

end
