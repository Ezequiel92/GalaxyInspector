####################################################################################################
# Filters
####################################################################################################

#################
# Base functions
#################

"""
    filterData!(data_dict::Dict; <keyword arguments>)::Nothing

Filter `data_dict` using the indices provided by `filter_function`.

# Arguments

  - `data_dict::Dict`: Data dictionary (see [`makeDataDict`](@ref) for the canonical description).
  - `filter_function::Function=filterNothing`: Filter function. See the required signature and examples in `./src/analysis/filters.jl`.
"""
function filterData!(data_dict::Dict; filter_function::Function=filterNothing)::Nothing

    # Compute the filter dictionary
    filter_dict = filter_function(data_dict)

    for component in snapshotTypes(data_dict)

        idxs = filter_dict[component]

        Threads.@threads for (block, values) in collect(pairs(data_dict[component]))

            if !isempty(values)

                data_dict[component][block] = collect(selectdim(values, ndims(values), idxs))

            end

        end

    end

    return nothing

end

"""
    filterData(data_dict::Dict; <keyword arguments>)::Dict

Return a filtered copy of `data_dict` using the indices provided by `filter_function`.

# Arguments

  - `data_dict::Dict`: Data dictionary (see [`makeDataDict`](@ref) for the canonical description).
  - `filter_function::Function=filterNothing`: Filter function. See the required signature and examples in `./src/analysis/filters.jl`.

# Returns

  - The filtered data.
"""
function filterData(data_dict::Dict; filter_function::Function=filterNothing)::Dict

    dd_copy = deepcopy(data_dict)

    if filter_function == filterNothing
        return dd_copy
    end

    # Compute the filter dictionary
    filter_dict = filter_function(dd_copy)

    for component in snapshotTypes(dd_copy)

        idxs = filter_dict[component]

        Threads.@threads for (block, values) in collect(pairs(dd_copy[component]))

            if !isempty(values)

                dd_copy[component][block] = collect(selectdim(values, ndims(values), idxs))

            end

        end

    end

    return dd_copy

end

"""
    intersectFilters(filters::Function...)::Function

Generate the filter function resulting from intersecting `filters` (using AND in boolean logic).

# Arguments

  - `filters`: Filter functions.

# Returns

  - The filter function resulting from intersecting `filters`.
"""
function intersectFilters(filters::Function...)::Function

    return dd -> intersectFilters([filter(dd) for filter in filters]...)

end

"""
    intersectFilters(filters...)::Dict{Symbol,IndexType}

Generate the filter dictionary resulting from intersecting `filters` (using AND in boolean logic).

# Arguments

  - `filters`: Filter dictionary.

# Returns

  - The filter dictionary resulting from intersecting `filters`.
"""
function intersectFilters(filters...)::Dict{Symbol,IndexType}

    if !allequal(keys.(filters))
        throw(ArgumentError("intersectFilters: The filters dictionaries must have the same list of \
        components (their keys)"))
    end

    filter_dict = Dict{Symbol,IndexType}()

    for component in keys(filters[1])

        filter_dict[component] = intersect(getindex.(filters, component)...)

    end

    return filter_dict

end

"""
    intersectFilters(filters::Dict{Symbol,IndexType})::Dict{Symbol,IndexType}

Trivial method for the case of a single argument.

# Arguments

  - `filters::Dict{Symbol,IndexType}`: Filter dictionary.

# Returns

  - The filter dictionary given as argument.
"""
intersectFilters(filters::Dict{Symbol,IndexType})::Dict{Symbol,IndexType} = filters

"""
    invertFilterDict(filter::Dict{Symbol,IndexType})::Dict{Symbol,IndexType}

Invert a given filter dictionary.

# Arguments

  - `filter::Dict{Symbol,IndexType}`: Filter dictionary to be inverted.

# Returns

  - The inverted filter dictionary.
"""
function invertFilterDict(filter::Dict{Symbol,IndexType})::Dict{Symbol,IndexType}

    filter_dict = Dict{Symbol,IndexType}()

    for (component, idxs) in pairs(filter)
        filter_dict[component] = Not(idxs)
    end

    return filter_dict

end

"""
    invertFilter(filter::Function)::Function

Invert a filter function

# Arguments

  - `filter::Function`: Filter function to be inverted.

# Returns

  - The inverted filter function.
"""
function invertFilter(filter::Function)::Function

    return (x...; y...) -> invertFilterDict(filter(x...; y...))

end

####################################################################################################
# Filter functions
####################################################################################################
#
# A filter function must take a data dictionary and return a filter dictionary.
#
# Expected signature:
#
#   filter_function(data_dict::Dict) -> filter_dict::Dict{Symbol,IndexType}
#
# where:
#
#   - `data_dict::Dict`: Data dictionary (see makeDataDict for the canonical description).
#   - filter_dict::Dict{Symbol,IndexType}: A dictionary with the following shape:
#
#      + `cell/particle type::Symbol` -> indices::IndexType
#      + `cell/particle type::Symbol` -> indices::IndexType
#      + `cell/particle type::Symbol` -> indices::IndexType
#      + ...
#
####################################################################################################

"""
Filter function that does not filter out any cell/particle.
"""
filterNothing(x...; y...)::Dict{Symbol,IndexType} = PASS_ALL

"""
Filter function that filters out all cells/particles.
"""
filterAll(x...; y...)::Dict{Symbol,IndexType} = PASS_NONE

"""
    filterBySphere(
        data_dict::Dict,
        min_r::Unitful.Length,
        max_r::Unitful.Length,
        origin...,
    )::Dict{Symbol,IndexType}

Select the cells/particles within a given spherical shell.

# Arguments

  - `data_dict::Dict`: Data dictionary (see [`makeDataDict`](@ref) for the canonical description).
  - `min_r::Unitful.Length`: Internal radius of the spherical shell.
  - `max_r::Unitful.Length`: External radius of the spherical shell.
  - `origin...`: It can be any number and type of argument compatible with the second to last arguments of a [`computeCenter`](@ref) method.

# Returns

  - A filter dictionary.
"""
function filterBySphere(
    data_dict::Dict,
    min_r::Unitful.Length,
    max_r::Unitful.Length,
    origin...,
)::Dict{Symbol,IndexType}

    (
        min_r >= max_r &&
        throw(ArgumentError("filterBySphere: `max_r` should be larger than `min_r`, \
        but I got `min_r` = $(min_r) >= `max_r` = $(max_r)"))
    )

    center = computeCenter(data_dict, origin...)

    filter_dict = Dict{Symbol,IndexType}(type => Int[] for type in snapshotTypes(data_dict))

    Threads.@threads for component in collect(keys(filter_dict))

        positions = data_dict[component]["POS "]

        if !isempty(positions)
            distances = colwise(Euclidean(), positions, center)
            filter_dict[component] = map(x -> min_r < x <= max_r, distances)
        end

    end

    return filter_dict

end

"""
    filterByCylinder(
        data_dict::Dict,
        max_r::Unitful.Length,
        max_z::Unitful.Length,
        origin...,
    )::Dict{Symbol,IndexType}

Select the cells/particles within a given cylinder.

# Arguments

  - `data_dict::Dict`: Data dictionary (see [`makeDataDict`](@ref) for the canonical description).
  - `max_r::Unitful.Length`: Radius of the cylinder.
  - `max_z::Unitful.Length`: Half height of the cylinder.
  - `origin...`: It can be any number and type of argument compatible with the second to last arguments of a [`computeCenter`](@ref) method.

# Returns

  - A filter dictionary.
"""
function filterByCylinder(
    data_dict::Dict,
    max_r::Unitful.Length,
    max_z::Unitful.Length,
    origin...,
)::Dict{Symbol,IndexType}

    (
        isPositive([max_r, max_z]) ||
        throw(ArgumentError("filterByCylinder: `max_r` and `max_z` should be larger than 0, \
        but I got `max_r` = $(max_r) and `max_z` = $(max_z)"))
    )

    center = computeCenter(data_dict, origin...)

    filter_dict = Dict{Symbol,IndexType}(type => Int[] for type in snapshotTypes(data_dict))

    Threads.@threads for component in collect(keys(filter_dict))

        positions = data_dict[component]["POS "]

        if !isempty(positions)
            distances = colwise(Euclidean(), positions, center)
            heights = abs.(positions[3, :])
            filter_dict[component] = map(r -> r <= max_r, distances) ∩ map(z -> z <= max_z, heights)
        end

    end

    return filter_dict

end

"""
    filterBySubhalo(
        data_dict::Dict;
        <keyword arguments>
    )::Dict{Symbol,IndexType}

Select the cells/particles that belong to a given halo and subhalo.

# Arguments

  - `data_dict::Dict`: Data dictionary (see [`makeDataDict`](@ref) for the canonical description).
  - `halo_idx::Int=1`: Index of the target halo (FoF group). Starts at 1.
  - `subhalo_rel_idx::Int=1`: Index of the target subhalo (subfind), relative to the target halo. Starts at 1. If it is set to 0, all subhalos of the target halo are included.

# Returns

  - A filter dictionary.
"""
function filterBySubhalo(
    data_dict::Dict;
    halo_idx::Int=1,
    subhalo_rel_idx::Int=1,
)::Dict{Symbol,IndexType}

    # If there are no subfind data, filter out every cell/particle
    if !isSubfindActive(data_dict[:gc_data].path)
        (
            logging[] &&
            @warn("filterBySubhalo: There is no subfind data in $(data_dict[:gc_data].path), \
            so every particle will be filtered out")
        )
        return PASS_NONE
    end

    # Load the necessary data
    g_n_subs   = data_dict[:group]["G_Nsubs"]
    g_len_type = data_dict[:group]["G_LenType"]
    s_len_type = data_dict[:subhalo]["S_LenType"]

    # If any of the data is missing return an empty filter dictionary
    n_groups_total = data_dict[:gc_data].header.n_groups_total

    if iszero(n_groups_total) || any(isempty, [g_n_subs, g_len_type, s_len_type])
        (
            logging[] &&
            @warn("filterBySubhalo: There is missing subfind data in $(data_dict[:gc_data].path), \
            so every particle will be filtered out")
        )
        return PASS_NONE
    end

    # Check that the requested halo index is within bounds
    (
        0 < halo_idx <= n_groups_total ||
        throw(ArgumentError("filterBySubhalo: There is only $(n_groups_total) FoF groups in \
        $(data_dict[:gc_data].path), so `halo_idx` = $(halo_idx) is out of bounds"))
    )

    # Compute the number of subhalos and particles up to the last halo before `halo_idx`
    if isone(halo_idx)
        n_subs_floor = 0
        len_type_floor = zeros(Int32, size(s_len_type, 1))
    else
        n_subs_floor = sum(g_n_subs[1:(halo_idx - 1)]; init=0)
        len_type_floor = sum(g_len_type[:, 1:(halo_idx - 1)], dims=2; init=0)
    end

    # Check that the requested subhalo index is within bounds
    n_subfinds = g_n_subs[halo_idx]

    if iszero(n_subfinds)
        (
            logging[] &&
            @warn("filterBySubhalo: There are 0 subhalos in the FoF group $(halo_idx) from
            $(data_dict[:gc_data].path), so every particle will be filtered out")
        )
        return PASS_NONE
    end

    (
        subhalo_rel_idx <= n_subfinds ||
        throw(ArgumentError("filterBySubhalo: There is only $(n_subfinds) subhalos for the \
        FoF group $(halo_idx) from $(data_dict[:gc_data].path), so `subhalo_rel_idx` = \
        $(subhalo_rel_idx) is out of bounds"))
    )

    if subhalo_rel_idx <= 0

        # Consider all subhalos within the target halo
        first_idxs = len_type_floor .+ 1
        last_idxs  = first_idxs .+ g_len_type[:, halo_idx] .- 1

    else

        # Compute the subhalo absolute index
        subhalo_abs_idx = n_subs_floor + subhalo_rel_idx

        # Compute the number of particles in the current halo,
        # upto the last subhalo before `subhalo_rel_idx`
        if isone(subhalo_abs_idx)
            len_type_floor_in_halo = zeros(Int, size(s_len_type, 1))
        else
            len_type_floor_in_halo = sum(
                s_len_type[:, (n_subs_floor + 1):(subhalo_abs_idx - 1)],
                dims=2;
                init=0,
            )
        end

        # Compute the first and last index of the selected
        # cells/particles (for each cell/particle type)
        first_idxs = len_type_floor .+ len_type_floor_in_halo .+ 1
        last_idxs  = first_idxs .+ s_len_type[:, subhalo_abs_idx] .- 1

    end

    filter_dict = Dict{Symbol,IndexType}()

    # Find which cell/particle types are part of the keys of `data_dict`
    components_in_dd = snapshotTypes(data_dict)

    # Fill the filter dictionary
    for (i, (first_idx, last_idx)) in enumerate(zip(first_idxs, last_idxs))

        component = INDEX_PARTICLE[i - 1]

        # Only compute the indices for components in `data_dict`
        if component ∈ components_in_dd

            if first_idx == last_idx || iszero(last_idx)
                filter_dict[component] = Int[]
            end

            if component == :stellar

                # Find the indices of the stars, excluding wind particles
                real_stars_idxs = findRealStars(data_dict[:snap_data].path)

                n_wind_before = count(!, real_stars_idxs[1:(first_idx - 1)])
                n_wind_between = count(!, real_stars_idxs[first_idx:last_idx])

                stars_first_idx = first_idx - n_wind_before
                stars_last_idx = last_idx - n_wind_before - n_wind_between

                filter_dict[component] = stars_first_idx:stars_last_idx

            else

                filter_dict[component] = first_idx:last_idx

            end

        end

    end

    return filter_dict

end

"""
    filterBySubhalo(data_dict::Dict, subhalo_abs_idx::Int)::Dict{Symbol,IndexType}

Select the cells/particles that belong to a given subhalo.

# Arguments

  - `data_dict::Dict`: Data dictionary (see [`makeDataDict`](@ref) for the canonical description).
  - `subhalo_abs_idx::Int`: Index of the target subhalo (subfind). Starts at 1.

# Returns

  - A filter dictionary.
"""
function filterBySubhalo(data_dict::Dict, subhalo_abs_idx::Int)::Dict{Symbol,IndexType}

    # If there are no subfind data, filter out every cell/particle
    if !isSubfindActive(data_dict[:gc_data].path)
        (
            logging[] &&
            @warn("filterBySubhalo: There is no subfind data in $(data_dict[:gc_data].path), \
            so every particle will be filtered out")
        )
        return PASS_NONE
    end

    # Load the necessary data
    s_len_type = data_dict[:subhalo]["S_LenType"]

    # If any of the data is missing return an empty filter dictionary
    n_subgroups_total = data_dict[:gc_data].header.n_subgroups_total

    if iszero(n_subgroups_total) || isempty(s_len_type)
        (
            logging[] &&
            @warn("filterBySubhalo: There is missing subfind data in $(data_dict[:gc_data].path), \
            so every particle will be filtered out")
        )
        return PASS_NONE
    end

    # Check that the requested subhalo index is within bounds
    (
        0 < subhalo_abs_idx <= n_subgroups_total ||
        throw(ArgumentError("filterBySubhalo: There is only $(n_subgroups_total) subhalos in \
        $(data_dict[:gc_data].path), so subhalo_abs_idx = $(subhalo_abs_idx) is out of bounds"))
    )

    # Compute the number of particles upto the last subhalo before `subhalo_abs_idx`
    if isone(subhalo_abs_idx)
        len_type_floor = zeros(Int, size(s_len_type, 1))
    else
        len_type_floor = sum(s_len_type[:, 1:(subhalo_abs_idx - 1)], dims=2; init=0)
    end

    # Compute the first and last index of the selected cells/particles (for each cell/particle type)
    first_idxs = len_type_floor .+ 1
    last_idxs  = len_type_floor .+ s_len_type[:, subhalo_abs_idx]

    filter_dict = Dict{Symbol,IndexType}()

    # Find which cell/particle types are part of the keys of `data_dict`
    components_in_dd = snapshotTypes(data_dict)

    # Fill the filter dictionary
    for (i, (first_idx, last_idx)) in enumerate(zip(first_idxs, last_idxs))

        component = INDEX_PARTICLE[i - 1]

        # Only compute the indices for components in `data_dict`
        if component ∈ components_in_dd

            if first_idx == last_idx || iszero(last_idx)
                filter_dict[component] = Int[]
            end

            if component == :stellar

                # Find the indices of the stars, excluding wind particles
                real_stars_idxs = findRealStars(data_dict[:snap_data].path)

                n_wind_before = count(!, real_stars_idxs[1:(first_idx - 1)])
                n_wind_between = count(!, real_stars_idxs[first_idx:last_idx])

                stars_first_idx = first_idx - n_wind_before
                stars_last_idx = last_idx - n_wind_before - n_wind_between

                filter_dict[component] = stars_first_idx:stars_last_idx

            else

                filter_dict[component] = first_idx:last_idx

            end

        end

    end

    return filter_dict

end

"""
    filterByQuantity(
        data_dict::Dict,
        quantity::Symbol,
        component::Symbol,
        min::Number,
        max::Number,
    )::Dict{Symbol,IndexType}

Select particles/cells with a value of `quantity` within [`min`, `max`].

# Arguments

  - `data_dict::Dict`: Data dictionary (see [`makeDataDict`](@ref) for the canonical description).
  - `quantity::Symbol`: Target quantity. For the possibilities see the documentation of [`scatterQty`](@ref).
  - `component::Symbol`: Type of particle/cell. The possibilities are the keys of [`PARTICLE_INDEX`](@ref).
  - `min::Number`: Minimum value of `quantity`.
  - `max::Number`: Maximum value of `quantity`.

# Returns

  - A filter dictionary.
"""
function filterByQuantity(
    data_dict::Dict,
    quantity::Symbol,
    component::Symbol,
    min::Number,
    max::Number,
)::Dict{Symbol,IndexType}

    (
        min >= max &&
        throw(ArgumentError("filterByQuantity: `max` should be larger than `min`, \
        but I got `min` = $(min) >= `max` = $(max)"))
    )

    filter_dict = Dict{Symbol,IndexType}(type => (:) for type in snapshotTypes(data_dict))

    # Compute the `quantity`
    values = scatterQty(data_dict, quantity)

    if isempty(values)
        filter_dict[component] = Int[]
    else
        filter_dict[component] = map(x -> min <= x <= max, values)
    end

    return filter_dict

end

"""
    filterBySFM(data_dict::Dict)::Dict{Symbol,IndexType}

Select the gas cells that have entered our star formation routine at least once.

# Arguments

  - `data_dict::Dict`: Data dictionary (see [`makeDataDict`](@ref) for the canonical description).

# Returns

  - A filter dictionary.
"""
function filterBySFM(data_dict::Dict)::Dict{Symbol,IndexType}

    filter_dict = Dict{Symbol,IndexType}(type => (:) for type in snapshotTypes(data_dict))

    fractions = data_dict[:gas]["FRAC"]

    if isempty(fractions)
        filter_dict[:gas] = Int[]
    else
        filter_dict[:gas] = map(!isnan, fractions[1, :])
    end

    return filter_dict

end

"""
    filterByStellarAge(
        data_dict::Dict;
        <keyword arguments>
    )::Dict{Symbol,IndexType}

Select stars with an age within [`min_age`, `max_age`].

# Arguments

  - `data_dict::Dict`: Data dictionary (see [`makeDataDict`](@ref) for the canonical description).
  - `min_age::Unitful.Time=0.0u"Gyr"`: Minimum age.
  - `max_age::Unitful.Time=AGE_RESOLUTION`: Maximum age.

# Returns

  - A filter dictionary.
"""
function filterByStellarAge(
    data_dict::Dict;
    min_age::Unitful.Time=0.0u"Gyr",
    max_age::Unitful.Time=AGE_RESOLUTION,
)::Dict{Symbol,IndexType}

    (
        min_age >= max_age &&
        throw(ArgumentError("filterByStellarAge: `max_age` should be larger than `min_age`, \
        but I got `min_age` = $(min_age) >= `max_age` = $(max_age)"))
    )

    filter_dict = Dict{Symbol,IndexType}(type => (:) for type in snapshotTypes(data_dict))

    ages = computeStellarAge(data_dict)

    if isempty(ages)
        filter_dict[:stellar] = Int[]
    else
        filter_dict[:stellar] = map(t -> min_age <= t <= max_age, ages)
    end

    return filter_dict

end

"""
    filterOldStars(data_dict::Dict)::Dict{Symbol,IndexType}

Select stars born after the previous snapshot.

# Arguments

  - `data_dict::Dict`: Data dictionary (see [`makeDataDict`](@ref) for the canonical description).

# Returns

  - A filter dictionary.
"""
function filterOldStars(data_dict::Dict)::Dict{Symbol,IndexType}

    filter_dict = Dict{Symbol,IndexType}(type => (:) for type in snapshotTypes(data_dict))

    birth_ticks = data_dict[:stellar]["GAGE"]

    if isempty(birth_ticks)

        filter_dict[:stellar] = Int[]

    else

        # Get the global index (index in the context of the whole simulation) of the current snapshot
        present_idx = data_dict[:snap_data].global_index

        if present_idx == 1

            filter_dict[:stellar] = (:)

        else

            # Compute the stellar birth dates
            if data_dict[:sim_data].cosmological
                # Go from scale factor to physical time
                birth_times = computeTime(birth_ticks, data_dict[:snap_data].header)
            else
                birth_times = birth_ticks
            end

            # Get the physical times
            times = data_dict[:sim_data].snapshot_table[!, :physical_times]

            filter_dict[:stellar] = map(t -> t >= times[present_idx - 1], birth_times)

        end

    end

    return filter_dict

end

"""
    filterByBirthPlace(
        data_dict::Dict,
        exclude::Symbol;
        <keyword arguments>
    )::Dict{Symbol,IndexType}

Select stars that were born either inside the given halo and subhalo (`exclude`= :exsitu), or outside (`exclude`= :insitu).

# Arguments

  - `data_dict::Dict`: Data dictionary (see [`makeDataDict`](@ref) for the canonical description).
  - `exclude::Symbol`: Which stars will be excluded, either the ones born outside the given halo and subhalo (:exsitu), or inside (:insitu).
  - `halo_idx::Int=1`: Index of the target halo (FoF group). Starts at 1.
  - `subhalo_rel_idx::Int=1`: Index of the target subhalo (subfind), relative to the target halo. Starts at 1. If it is set to 0, all subhalos of the target halo are considered in-situ.

# Returns

  - A filter dictionary.
"""
function filterByBirthPlace(
    data_dict::Dict,
    exclude::Symbol;
    halo_idx::Int=1,
    subhalo_rel_idx::Int=1,
)::Dict{Symbol,IndexType}

    filter_dict = Dict{Symbol,IndexType}(type => (:) for type in snapshotTypes(data_dict))

    birth_halo, birth_subhalo = locateStellarBirthPlace(data_dict)

    if any(isempty, [birth_halo, birth_subhalo])

        filter_dict[:stellar] = Int[]

    else

        stars_born_in_halo = map(isequal(halo_idx), birth_halo)

        if iszero(subhalo_rel_idx)
            stars_born_in_subhalo = (:)
        else
            stars_born_in_subhalo = map(isequal(subhalo_rel_idx), birth_subhalo)
        end

        insitu_stars = stars_born_in_halo ∩ stars_born_in_subhalo

        if exclude == :insitu

            filter_dict[:stellar] = Vector{Bool}(.!(insitu_stars))

        elseif exclude == :exsitu

            filter_dict[:stellar] = insitu_stars

        else

            throw(ArgumentError("filterByBirthPlace: `exclude` can only be :insitu or :exsitu, \
            but I got :$(exclude)"))

        end

    end

    return filter_dict

end

#################
# Filter presets
#################

"""
    selectFilter(
        filter_mode::Symbol,
        base_request::Dict{Symbol,Vector{String}},
    )::Tuple{Function,Dict{Symbol,Vector{String}}}

Select a filter function from a list of premade ones.

Creates a request dictionary, using `base_request` as a base, adding what is necessary for the filter function.

# Arguments

  - `filter_mode::Symbol`: Which cells/particles will be selected, the options are:

      + `:all`             -> Select every cell/particle within the simulation box.
      + `:halo`            -> Select the cells/particles that belong to the main halo.
      + `:subhalo`         -> Select the cells/particles that belong to the main subhalo.
      + `:sphere`          -> Select the cells/particles inside a sphere centered at the origin with radius `DISK_R` (see `./src/constants/globals.jl`).
      + `:cylinder`        -> Select the cells/particles inside a cylinder centered at the origin with radius `DISK_R` and half height `DISK_HEIGHT` (see `./src/constants/globals.jl`). The cylinder axis is aligned with the z axis.
  - `base_request::Dict{Symbol,Vector{String}}`: Base request dictionary.

# Returns

  - A Tuple with two elements:

      + The filter function.
      + The new request dictionary.
"""
function selectFilter(
    filter_mode::Symbol,
    base_request::Dict{Symbol,Vector{String}},
)::Tuple{Function,Dict{Symbol,Vector{String}}}

    if filter_mode == :all

        # Select every cell/particle within the simulation box
        filter_function = filterNothing

        new_request = base_request

    elseif filter_mode == :halo

        # Select the cells/particles that belong to the main halo
        filter_function = dd -> filterBySubhalo(dd; halo_idx=1, subhalo_rel_idx=0)

        new_request = mergeRequests(
            base_request,
            Dict(:group => ["G_Nsubs", "G_LenType"], :subhalo => ["S_LenType"]),
        )

    elseif filter_mode == :subhalo

        # Select the cells/particles that belong to the main subhalo
        filter_function = dd -> filterBySubhalo(dd)

        new_request = mergeRequests(
            base_request,
            Dict(:group => ["G_Nsubs", "G_LenType"], :subhalo => ["S_LenType"]),
        )

    elseif filter_mode == :sphere

        # Select the cells/particles inside a sphere centered at the origin with radius `DISK_R`
        filter_function = dd -> filterBySphere(dd, 0.0u"kpc", DISK_R, :zero)

        new_request = mergeRequests(
            base_request,
            Dict(component => ["POS "] for component in keys(PARTICLE_INDEX)),
        )

    elseif filter_mode == :cylinder

        # Select the cells/particles inside a cylinder centered at the origin with radius `DISK_R`
        # and half height `DISK_HEIGHT`. The cylinder axis is aligned with the z axis
        filter_function = dd -> filterByCylinder(dd, DISK_R, DISK_HEIGHT, :zero)

        new_request = mergeRequests(
            base_request,
            Dict(component => ["POS "] for component in keys(PARTICLE_INDEX)),
        )

    else

        throw(ArgumentError("selectFilter: `filter_mode` can only be :all, :halo, :subhalo, \
        :sphere or :cylinder, but I got :$(filter_mode)"))

    end

    return filter_function, new_request

end

"""
    selectFilter(
        filter_mode::Tuple{Function,Dict{Symbol,Vector{String}}},
        base_request::Dict{Symbol,Vector{String}},
    )::Tuple{Function,Dict{Symbol,Vector{String}}}

Compatibility method to pass personalized filters to functions that only use the `filter_mode` argument.

# Arguments

  - `filter_mode::Tuple{Function,Dict{Symbol,Vector{String}}}`: Filter function and the request that goes with it.
  - `base_request::Dict{Symbol,Vector{String}}`: Base request dictionary.

# Returns

  - A Tuple with two elements:

      + The filter function.
      + The new request dictionary.
"""
function selectFilter(
    filter_mode::Tuple{Function,Dict{Symbol,Vector{String}}},
    base_request::Dict{Symbol,Vector{String}},
)::Tuple{Function,Dict{Symbol,Vector{String}}}

    filter_function, ff_request = filter_mode

    new_request = mergeRequests(base_request, ff_request)

    return filter_function, new_request

end
