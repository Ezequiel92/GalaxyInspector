####################################################################################################
# Coordinate transformations
####################################################################################################

#############################
# Translation base functions
#############################

"""
    translatePoints!(points::Matrix{<:Number}, origin::Vector{<:Number})::Nothing

Translate a system of points, moving `origin` to [0, 0, 0].

# Arguments

  - `points::Matrix{<:Number}`: Points to be translated. Each column is a point and each row a dimension.
  - `origin::Vector{<:Number}`: Target origin.
"""
function translatePoints!(points::Matrix{<:Number}, origin::Vector{<:Number})::Nothing

    all(iszero, origin) && return nothing

    Threads.@threads for i in axes(points, 2)
        points[:, i] .-= origin
    end

    return nothing

end

###############################
# Translation derive functions
###############################

"""
    translateData!(
        data_dict::Dict,
        origin::Vector{<:Unitful.Length},
        vcm::Vector{<:Unitful.Velocity},
    )::Nothing

Translate the positions and velocities of the cells/particles in `data_dict` such that `origin` and `vcm` end up both being [0, 0, 0].

# Arguments

  - `data_dict::Dict`: Data dictionary (see [`makeDataDict`](@ref) for the canonical description).
  - `origin::Vector{<:Unitful.Length}`: Target origin.
  - `vcm::Vector{<:Unitful.Velocity}`: Velocity of the center of mass.
"""
function translateData!(
    data_dict::Dict,
    origin::Vector{<:Unitful.Length},
    vcm::Vector{<:Unitful.Velocity},
)::Nothing

    for component in snapshotTypes(data_dict)

        data = data_dict[component]

        if haskey(data, "POS ") && !isempty(data["POS "])
            translatePoints!(data["POS "], origin)
        end

        if haskey(data, "VEL ") && !isempty(data["VEL "])
            translatePoints!(data["VEL "], vcm)
        end

    end

    return nothing

end

"""
    translateData!(data_dict::Dict, translation::Union{Symbol,NTuple{2,Int},Int}...)::Nothing

Translate the positions and velocities of the cells/particles in `data_dict`.

# Arguments

  - `data_dict::Dict`: Data dictionary (see [`makeDataDict`](@ref) for the canonical description).
  - `translation::Union{Symbol,NTuple{2,Int},Int}`: Target translation. The options are:

      + `:zero`                       -> No translation is applied.
      + `:all`                        -> Sets the center of mass of the whole system as the new origin.
      + `:{component}`                -> Sets the center of mass of the given component (e.g. :stellar, :gas, :dark_matter, etc) as the new origin. It can be any of the keys of [`PARTICLE_INDEX`](@ref).
      + `(halo_idx, subhalo_rel_idx)` -> Sets the position of the potential minimum for the `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo) as the new origin.
      + `(halo_idx, 0)`               -> Sets the center of mass of the `halo_idx::Int` halo as the new origin.
      + `subhalo_abs_idx`             -> Sets the center of mass of the `subhalo_abs_idx::Int` as the new origin.
"""
function translateData!(data_dict::Dict, translation::Union{Symbol,NTuple{2,Int},Int}...)::Nothing

    origin = computeCenter(data_dict, translation...)
    vcm    = computeVcm(data_dict, translation...)

    translateData!(data_dict, origin, vcm)

    return nothing

end

############################
# Rotation derive functions
############################

"""
    rotateData!(
        data_dict::Dict,
        rotation_matrix::Union{Matrix{Float64},UniformScaling{Bool}},
    )::Nothing

Rotate the positions and velocities of the cells/particles in `data_dict`.

# Arguments

  - `data_dict::Dict`: Data dictionary (see [`makeDataDict`](@ref) for the canonical description).
  - `rotation_matrix::Union{Matrix{Float64},UniformScaling{Bool}}`: Rotation matrix.
"""
function rotateData!(
    data_dict::Dict,
    rotation_matrix::Union{Matrix{Float64},UniformScaling{Bool}},
)::Nothing

    for component in snapshotTypes(data_dict)

        data = data_dict[component]

        if haskey(data, "POS ") && !isempty(data["POS "])
            data["POS "] = rotation_matrix * data["POS "]
        end

        if haskey(data, "VEL ") && !isempty(data["VEL "])
            data["VEL "] = rotation_matrix * data["VEL "]
        end

    end

    return nothing

end

"""
    rotateData!(
        data_dict::Dict,
        z_axis::Symbol,
        component::Symbol,
        filter_function::Function,
    )::Nothing

Rotate the positions and velocities of the cells/particles in `data_dict`.

# Arguments

  - `data_dict::Dict`: Data dictionary (see [`makeDataDict`](@ref) for the canonical description).
  - `z_axis::Symbol`: Target reference system axis. The options are:

      + `:zero` -> No rotation is applied.
      + `:am`   -> The angular momentum will be used as the new z axis.
      + `:pa`   -> The principal axes will be used as the new coordinate system.
  - `component::Symbol`: Which component will be considered to compute the angular momentum or the principal axes. The options are:

      + `:all`         -> Every component present in `data_dict`.
      + `:{component}` -> Any of the keys of [`PARTICLE_INDEX`](@ref), if present in `data_dict`.
  - `filter_function::Function`: A function with the signature:

    `filter_function(data_dict::Dict) -> filter_dict::Dict{Symbol,IndexType}`

    Determines which cell/particles will be considered to compute the angular momentum or the principal axes.
"""
function rotateData!(
    data_dict::Dict,
    z_axis::Symbol,
    component::Symbol,
    filter_function::Function,
)::Nothing

    z_axis == :zero && return nothing

    filtered_dd = filterData(data_dict; filter_function)

    if z_axis == :am

        rotation_matrix = computeAMRotationMatrix(filtered_dd, component)

    elseif z_axis ===:pa

        rotation_matrix = computePARotationMatrix(filtered_dd, component)

    else

        throw(ArgumentError("rotateData!: `z_axis` can only be :zero, :am or :pa, but I got \
        :$(z_axis)"))

    end

    rotateData!(data_dict, rotation_matrix)

    return nothing

end

##########################
# Transformations presets
##########################

"""
    selectTransformation(
        trans_mode::Symbol,
        base_request::Dict{Symbol,Vector{String}},
    )::Tuple{TranslationType,RotationType,Dict{Symbol,Vector{String}}}

Select a translation and rotation from a list of premade ones.

Creates a request dictionary, using `request` as a base, adding what is necessary for the corresponding transformations.

# Arguments

  - `trans_mode::Symbol`: How to translate and rotate the cells/particles, after filtering with `filter_mode`. The options are the elements of [`TRANSFORM_LIST`](@ref) which have the form :{component}_{group} selecting which cells/particle to consider for the center of mass (new origin -> translation) and principal axis (new reference system -> rotation), where group can be:

      + :box     -> Whole simulation box
      + :halo    -> Main halo
      + :subhalo -> Main subhalo
    and components can be:

      + :all         -> Every component present in data_dict
      + :{component} -> One of the keys of PARTICLE_INDEX
  - `base_request::Dict{Symbol,Vector{String}}`: Base request dictionary.

# Returns

  - A Tuple with tree elements:

      + The translation (see the posible arguments of [`translateData!`](@ref)).
      + The rotation (see the posible arguments of [`rotateData!`](@ref)).
      + The new request dictionary.
"""
function selectTransformation(
    trans_mode::Symbol,
    base_request::Dict{Symbol,Vector{String}},
)::Tuple{TranslationType,RotationType,Dict{Symbol,Vector{String}}}

    (
        trans_mode ∈ TRANSFORM_LIST ||
        throw(ArgumentError("selectTransformation: `trans_mode` can only be one of the \
        elements of TRANSFORM_LIST, but I got :$(trans_mode)"))
    )

    component, group = TRANSFORM_SPLITS[trans_mode]

    if component == :all
        pa_request = Dict(st => ["POS ", "MASS", "VEL "] for st in keys(PARTICLE_INDEX))
    else
        pa_request = Dict(component => ["POS ", "MASS", "VEL "])
    end

    if group == :box

        translation     = component
        filter_function = filterNothing
        new_request     = mergeRequests(base_request, pa_request)

    elseif group == :halo

        translation     = (1, 0)
        filter_function = dd->filterBySubhalo(dd; halo_idx=1, subhalo_rel_idx=0)
        new_request     = mergeRequests(
            base_request,
            pa_request,
            Dict(
                :group => ["G_Nsubs", "G_LenType", "G_Pos", "G_Vel"],
                :subhalo => ["S_LenType", "S_Pos", "S_Vel"],
            ),
        )

    elseif group == :subhalo

        translation     = (1, 1)
        filter_function = dd->filterBySubhalo(dd; halo_idx=1, subhalo_rel_idx=1)
        new_request     = mergeRequests(
            base_request,
            pa_request,
            Dict(
                :group => ["G_Nsubs", "G_LenType", "G_Pos", "G_Vel"],
                :subhalo => ["S_LenType", "S_Pos", "S_Vel"],
            ),
        )

    end

    return translation, (:pa, component, filter_function), new_request

end

"""
    selectTransformation(
        trans_mode::Tuple{TranslationType,RotationType,Dict{Symbol,Vector{String}}},
        base_request::Dict{Symbol,Vector{String}},
    )::Tuple{TranslationType,RotationType,Dict{Symbol,Vector{String}}}

Compatibility method to pass personalized transformations to functions that only use the `trans_mode` argument.

# Arguments

  - `trans_mode::Tuple{TranslationType,RotationType,Dict{Symbol,Vector{String}}}`: Translation, rotation, and the request that goes with them. See the posible arguments of [`translateData!`](@ref) and [`rotateData!`](@ref).
  - `base_request::Dict{Symbol,Vector{String}}`: Base request dictionary.

# Returns

  - A Tuple with tree elements:

      + The translation (see the posible arguments of [`translateData!`](@ref)).
      + The rotation (see the posible arguments of [`rotateData!`](@ref)).
      + The new request dictionary.
"""
function selectTransformation(
    trans_mode::Tuple{TranslationType,RotationType,Dict{Symbol,Vector{String}}},
    base_request::Dict{Symbol,Vector{String}},
)::Tuple{TranslationType,RotationType,Dict{Symbol,Vector{String}}}

    translation, rotation, trans_request = trans_mode

    new_request = mergeRequests(base_request, trans_request )

    return translation, rotation, new_request

end
