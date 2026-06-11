####################################################################################################
# Utilities to manage AbstractPlotQuantity
####################################################################################################

"""
Convenience entry point for quantity retrieval. If the input is already an AbstractPlotQuantity, it is returned as is. If the input is a Symbol, it is looked up in `QTY_REGISTRY` and the corresponding BaseQuantity is returned.
"""
Q(x::AbstractPlotQuantity) = x
Q(x::Symbol)::BaseQuantity = QTY_REGISTRY[x]

# New methods for the standard operators to be used with AbstractPlotQuantity.
for (op, label) in (
    (:+, "+"),
    (:-, "-"),
    (:*, "\\cdot"),
    (:/, "/"),
)
    @eval begin
        Base.$op(
            a::Union{Symbol,AbstractPlotQuantity},
            b::Union{Symbol,AbstractPlotQuantity},
        ) = BinaryQuantity($op, $label, Q(a), Q(b))
    end
end

"""
Read the label of a [`BaseQuantity`](@ref).
"""
getQuantityLabel(q::BaseQuantity)::AbstractString = q.qty_label

"""
    getQuantityLabel(q::BinaryQuantity)::LaTeXString

Construct the label of a [`BinaryQuantity`](@ref).

# Arguments

  - `q::BinaryQuantity`: Target quantity.

# Returns

  - The label of `q`.
"""
function getQuantityLabel(q::BinaryQuantity)::LaTeXString

    left_label  = getQuantityLabel(q.left)
    right_label = getQuantityLabel(q.right)

    return LaTeXString(left_label * L"\," * q.op_label * L"\," * right_label)

end

"""
Read the units of a [`BaseQuantity`](@ref).
"""
getQuantityUnit(q::BaseQuantity)::Unitful.Units = q.unit

"""
    getQuantityUnit(q::BinaryQuantity)::Unitful.Units

Construct the units of a [`BinaryQuantity`](@ref).

# Arguments

  - `q::BinaryQuantity`: Target quantity.

# Returns

  - The units of `q`.
"""
function getQuantityUnit(q::BinaryQuantity)::Unitful.Units

    unit_l = getQuantityUnit(q.left)
    unit_r = getQuantityUnit(q.right)

    # Addition/subtraction logic
    if q.op == (+) || q.op == (-)

        if Unitful.dimension(unit_l) == Unitful.dimension(unit_r)
            return unit_r
        else
            throw(ArgumentError("getQuantityUnit: Dimension mismatch for $(q.op): \
            $(unit_l) and $(unit_r)"))
        end

    # Multiplication/division logic
    else

        return q.op(unit_l, unit_r)

    end

end

"""
Read the request dictionary of a [`BaseQuantity`](@ref).
"""
getQuantityRequest(q::BaseQuantity)::Dict{Symbol,Vector{String}} = q.request

"""
    getQuantityRequest(q::BinaryQuantity)::Dict{Symbol,Vector{String}}

Construct the request dictionary of a [`BinaryQuantity`](@ref).

# Arguments

  - `q::BinaryQuantity`: Target quantity.

# Returns

  - The request dictionary of `q`.
"""
function getQuantityRequest(q::BinaryQuantity)::Dict{Symbol,Vector{String}}
    return mergeRequests(getQuantityRequest(q.left), getQuantityRequest(q.right))
end

"""
Read the exponential factor of a [`BaseQuantity`](@ref).
"""
getQuantityExpFactor(q::BaseQuantity)::Int = q.exp_factor

"""
For [`BinaryQuantity`](@ref) the exponential factor is always 0.
"""
getQuantityExpFactor(q::BinaryQuantity)::Int = 0

"""
Read the name of a quantity. For a [`BaseQuantity`](@ref) it is just its `id`. For a [`BinaryQuantity`](@ref) it is a Symbol with the names of the left and right quantities and the operator in between.
"""
getQuantityName(q::Symbol)::Symbol = q
getQuantityName(q::BaseQuantity)::Symbol = q.id
getQuantityName(q::BinaryQuantity)::Symbol = Symbol(
    "$(getQuantityName(q.left))_op_$(getQuantityName(q.right))"
)

"""
    scatterQty(data_dict::Dict, quantity::Symbol)::Vector{<:Number}

Convenience entry point for [`scatterQty`](@ref).

# Arguments

  - `data_dict::Dict`: Data dictionary. See [`makeDataDict`](@ref) for a canonical description.
  - `quantity::Symbol`: Target quantity. See the keys of [`QTY_REGISTRY`](@ref) for options.

# Returns

  - The value of `quantity` for every cell/particle.
"""
function scatterQty(data_dict::Dict, quantity::Symbol)::Vector{<:Number}

    (
        haskey(QTY_REGISTRY, quantity) ||
        throw(ArgumentError("scatterQty: I don't recognize the quantity :$(quantity). \
        It is not a key of QTY_REGISTRY"))
    )

    return scatterQty(data_dict, QTY_REGISTRY[quantity])

end

"""
    scatterQty(data_dict::Dict, quantity::BaseQuantity)::Vector{<:Number}

Compute `quantity` for each cell/particle in `data_dict`.

# Arguments

  - `data_dict::Dict`: Data dictionary. See [`makeDataDict`](@ref) for a canonical description.
  - `quantity::BaseQuantity`: Target base quantity. See the values of [`QTY_REGISTRY`](@ref) for options.

# Returns

  - The value of `quantity` for every cell/particle. If `quantity` is not well defined for each cell/particle, an empty vector is returned.
"""
function scatterQty(data_dict::Dict, quantity::BaseQuantity)::Vector{<:Number}

    scatter_qty = quantity.scatter_func(data_dict)

    if isnothing(scatter_qty)

        LOGGING[] && @warn("scatterQty: :$(quantity) is not well defined for each cell/particle")

        return typeof(1.0 * quantity.unit)[]

    end

    return scatter_qty

end

"""
    scatterQty(data_dict::Dict, quantity::BinaryQuantity)::Vector{<:Number}

Compute `quantity` for each cell/particle in `data_dict`.

# Arguments

  - `data_dict::Dict`: Data dictionary. See [`makeDataDict`](@ref) for a canonical description.
  - `quantity::BinaryQuantity`: Target binary quantity.

# Returns

  - The value of `quantity` for every cell/particle. If `quantity` is not well defined for each cell/particle, an empty vector is returned.
"""
function scatterQty(data_dict::Dict, quantity::BinaryQuantity)::Vector{<:Number}

    scatter_left  = scatterQty(data_dict, quantity.left)
    scatter_right = scatterQty(data_dict, quantity.right)

    if length(scatter_left) != length(scatter_right)
        throw(ArgumentError("scatterQty: Binary quantity :$(quantity) has an unequal number of \
        cells/particles on each side of its operation: # left = $(length(scatter_left)) != \
        # right = $(length(scatter_right)). Remember that for scatterQty() the sides of a binary \
        quantity should be of the same cell/particle type."))
    end

    unit_l = getQuantityUnit(quantity.left)
    unit_r = getQuantityUnit(quantity.right)

    op = quantity.op

    if (op == (+) || op == (-)) && (Unitful.dimension(unit_l) != Unitful.dimension(unit_r))
        throw(ArgumentError("scatterQty: Dimension mismatch for $(op): $(unit_l) and $(unit_r)"))
    end

    return op.(scatter_left, scatter_right)

end

"""
    integrateQty(
        data_dict::Dict,
        quantity::Symbol;
        <keyword arguments>
    )::Number

Convenience entry point for [`integrateQty`](@ref).

# Arguments

  - `data_dict::Dict`: Data dictionary. See [`makeDataDict`](@ref) for a canonical description.
  - `quantity::Symbol`: Target quantity. See the keys of [`QTY_REGISTRY`](@ref) for options.
  - `agg_function::Union{Function,Symbol}=:default`: If `quantity` is well defined for each cell/particle, you can pass an `agg_function` to accumulate the values calculated with [`scatterQty`](@ref). If `agg_function` is left as `:default` [`integrateQty`](@ref) will try to compute the most reasonable global value for `quantity`.

# Returns

  - The value of `quantity` for the whole system of cell/particles in `data_dict`. If `quantity` is not a well defined global quantity, NaN is returned.
"""
function integrateQty(
    data_dict::Dict,
    quantity::Symbol;
    agg_function::Union{Function,Symbol}=:default,
)::Number

    (
        haskey(QTY_REGISTRY, quantity) ||
        throw(ArgumentError("integrateQty: I don't recognize the quantity :$(quantity). \
        It is not a key of QTY_REGISTRY"))
    )

    return integrateQty(data_dict, QTY_REGISTRY[quantity]; agg_function)

end

"""
    integrateQty(
        data_dict::Dict,
        quantity::BaseQuantity;
        <keyword arguments>
    )::Number

Compute `quantity` for the whole system of cell/particles in `data_dict`.

# Arguments

  - `data_dict::Dict`: Data dictionary. See [`makeDataDict`](@ref) for a canonical description.
  - `quantity::BaseQuantity`: Target base quantity. See the values of [`QTY_REGISTRY`](@ref) for options.
  - `agg_function::Union{Function,Symbol}=:default`: If `quantity` is well defined for each cell/particle, you can pass an `agg_function` to accumulate the values calculated with [`scatterQty`](@ref). If `agg_function` is left as `:default` [`integrateQty`](@ref) will try to compute the most reasonable global value for `quantity`.

# Returns

  - The value of `quantity` for the whole system of cell/particles in `data_dict`. If `quantity` is not a well defined global quantity, NaN is returned.
"""
function integrateQty(
    data_dict::Dict,
    quantity::BaseQuantity;
    agg_function::Union{Function,Symbol}=:default,
)::Number

    if agg_function == :default

        integrated_qty = quantity.integrate_func(data_dict)

        if isnothing(integrated_qty)

            LOGGING[] && @warn("integrateQty: :$(quantity) is not a well defined global quantity")

            return NaN * quantity.unit

        end

        return integrated_qty

    else

        return agg_function(scatterQty(data_dict, quantity))

    end

end

"""
    integrateQty(
        data_dict::Dict,
        quantity::BinaryQuantity;
        <keyword arguments>
    )::Number

Compute `quantity` for the whole system of cell/particles in `data_dict`.

!!! note

    The binary operations are applied after computing [`integrateQty`](@ref) on each of the base quantities involved.

# Arguments

  - `data_dict::Dict`: Data dictionary. See [`makeDataDict`](@ref) for a canonical description.
  - `quantity::BinaryQuantity`: Target binary quantity.
  - `agg_function::Union{Function,Symbol}=:default`: If `quantity` is well defined for each cell/particle, you can pass an `agg_function` to accumulate the values calculated with [`scatterQty`](@ref). If `agg_function` is left as `:default` [`integrateQty`](@ref) will try to compute the most reasonable global value for `quantity`.

# Returns

  - The value of `quantity` for the whole system of cell/particles in `data_dict`. If `quantity` is not a well defined global quantity, NaN is returned.
"""
function integrateQty(
    data_dict::Dict,
    quantity::BinaryQuantity;
    agg_function::Union{Function,Symbol}=:default,
)::Number

    int_left  = integrateQty(data_dict, quantity.left; agg_function)
    int_right = integrateQty(data_dict, quantity.right; agg_function)

    unit_l = getQuantityUnit(quantity.left)
    unit_r = getQuantityUnit(quantity.right)

    op = quantity.op

    if (op == (+) || op == (-)) && (Unitful.dimension(unit_l) != Unitful.dimension(unit_r))
        throw(ArgumentError("integrateQty: Dimension mismatch for $(op): $(unit_l) and $(unit_r)"))
    end

    return op(int_left, int_right)

end

"""
    applyIntegrator(
        data_dict::Dict,
        scatter_function::Function,
        agg_function::Function,
        log::Union{Unitful.Units,Nothing},
    )::Number

Compute an integrated value from the results of `scatter_function` using `agg_function`. Optionally apply a log10 transformation.

# Arguments

  - `data_dict::Dict`: Data dictionary. See [`makeDataDict`](@ref) for a canonical description.
  - `scatter_function::Function`: Function that computes a quantity for each cell/particle.
  - `agg_function::Function`: Function that aggregates the scattered values.
  - `log::Union{Unitful.Units,Nothing}`: If not `nothing`, apply log10 to the scattered values after removing the units.

# Returns

  - The aggregated value of the scattered quantity.
"""
function applyIntegrator(
    data_dict::Dict,
    scatter_function::Function,
    agg_function::Function,
    log::Union{Unitful.Units,Nothing},
)::Number

    vals = scatter_function(data_dict)

    if !isnothing(log)
        vals = ustrip.(log, vals)
        filter_func = x -> !isinf(x) && isPositive(x)
    else
        filter_func = x -> !isinf(x)
    end

    filter!(filter_func, vals)

    if !isnothing(log)
        vals = log10.(vals)
    end

    isempty(vals) && return NaN

    return agg_function(vals)

end
