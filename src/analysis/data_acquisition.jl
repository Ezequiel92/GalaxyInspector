####################################################################################################
# Data acquisition functions
####################################################################################################

"""
    isBlockPresent(block::String, group::HDF5.Group)::Bool

Checks if a given data block exists in a HDF5 group.

# Arguments

  - `block::String`: Target block. The possibilities are the keys of [`QUANTITIES`](@ref).
  - `group::HDF5.Group`: HDF5 group.

# Returns

  - If `block` exists in `group`.
"""
function isBlockPresent(block::String, group::HDF5.Group)::Bool

    (
        haskey(QUANTITIES, block) ||
        throw(ArgumentError("isBlockPresent: `block` should be a key of `QUANTITIES`, \
        but I got $(block), see the options in `./src/constants/arepo.jl`"))
    )

    return haskey(group, QUANTITIES[block].hdf5_name)

end

"""
    getBlockSize(block::String, group::HDF5.Group)::Int

Compute the total size in bytes of a given data block in an HDF5 group.

# Arguments

  - `block::String`: Target block. The possibilities are the keys of [`QUANTITIES`](@ref).
  - `group::HDF5.Group`: The HDF5 group containing the data block.

# Returns

  - The total size in bytes of the data block.
"""
function getBlockSize(block::String, group::HDF5.Group)::Int

    if !isBlockPresent(block, group)

        logging[] && @warn("getBlockSize: The block $(block) is missing from the HDF5 group")

        return 0

    end

    dataset = group[QUANTITIES[block].hdf5_name]

    # Compute the total size in bytes
    return length(dataset) * sizeof(eltype(dataset))

end

"""
    getBlockSize(blocks::Vector{String}, group::HDF5.Group)::Int

Compute the total size in bytes of a set of data blocks in an HDF5 group.

# Arguments

  - `blocks::Vector{String}`: A vector of the target data blocks. The possibilities are the keys of [`QUANTITIES`](@ref).
  - `group::HDF5.Group`: The HDF5 group containing the data blocks.

# Returns

  - The total size in bytes of the data blocks.
"""
function getBlockSize(blocks::Vector{String}, group::HDF5.Group)::Int

    return sum(b->getBlockSize(b, group), blocks)

end

"""
    getRequestSize(request::Dict{Symbol,Vector{String}}, path::String)::Int

Compute the total size in bytes of the data blocks specified in a request dictionary.

# Arguments

  - `request::Dict{Symbol,Vector{String}}`: Dictionary with the shape `cell/particle type` -> [`block`, `block`, ...], where the possible types are the keys of [`PARTICLE_INDEX`](@ref), :group, and :subhalo, and the possible blocks are the keys of [`QUANTITIES`](@ref).
  - `path::String`: Path to the snapshot or group catalog file or folder.

# Returns

  - The total size in bytes of the data blocks specified in the request.
"""
function getRequestSize(request::Dict{Symbol,Vector{String}}, path::String)::Int

    if isfile(path)

        (
            HDF5.ishdf5(path) ||
            throw(ArgumentError("getRequestSize: The file $(path) is not in the HDF5 format, \
            I don't know how to read it"))
        )

        total_size = h5open(path, "r") do hdf5_file

            ts = 0

            for component in keys(request)

                if (
                    haskey(PARTICLE_CODE_NAME, component) &&
                    haskey(hdf5_file, PARTICLE_CODE_NAME[component])
                )

                    group = hdf5_file[PARTICLE_CODE_NAME[component]]

                    blocks = request[component]

                    ts += getBlockSize(blocks, group)

                elseif haskey(hdf5_file, titlecase(string(component)))

                    group = hdf5_file[titlecase(string(component))]

                    blocks = request[component]

                    ts += getBlockSize(blocks, group)

                end

            end

            (
                iszero(ts) && logging[] &&
                @warn("getRequestSize: The requested blocks: \n\n$(request)\n\n for $(path) \
                are all missing, so the total size is zero")
            )

            ts

        end

        return total_size

    elseif isdir(path)

        snap_files = findFiles(path, "$(SNAP_BASENAME)_*.*.hdf5")
        group_files = findFiles(path, "$(GC_BASENAME)_*.*.hdf5")

        sub_files = vcat(snap_files, group_files)

        (
            !isempty(sub_files) && all(HDF5.ishdf5, sub_files) ||
            throw(ArgumentError("getRequestSize: The directory $(path) does not contain \
            snapshot or group catalog sub-files in the HDF5 format"))
        )

        return sum(sub_file->getRequestSize(request, sub_file), sub_files)

    else

        throw(ArgumentError("getRequestSize: $(path) does not exists as a file or folder"))

    end

end

"""
    getBlockDims(block::String, group::HDF5.Group)::Tuple

Get the dimensions of a given data block in an HDF5 group.

# Arguments

  - `block::String`: Target block. The possibilities are the keys of [`QUANTITIES`](@ref).
  - `group::HDF5.Group`: HDF5 group.

# Returns

  - The dimensions of the data block as a tuple.
"""
function getBlockDims(block::String, group::HDF5.Group)::Tuple

    if isBlockPresent(block, group)

        dataset = group[QUANTITIES[block].hdf5_name]

        return size(dataset)

    else

        logging[] && @warn("getBlockDims: The block $(block) is missing from the HDF5 group")

        return (0,)

    end

end

"""
    getBlockType(block::String, group::HDF5.Group)::Type

Get the number type of a given data block in an HDF5 group.

# Arguments

  - `block::String`: Target block. The possibilities are the keys of [`QUANTITIES`](@ref).
  - `group::HDF5.Group`: HDF5 group.

# Returns

  - The data type of the block.
"""
function getBlockType(block::String, group::HDF5.Group)::Type

    if isBlockPresent(block, group)

        return eltype(group[QUANTITIES[block].hdf5_name])

    else

        logging[] && @warn("getBlockType: The block $(block) is missing from the HDF5 file")

        return Float64

    end
end

"""
    readBlock(
        block::String,
        group::HDF5.Group;
        <keyword arguments>
    )::VecOrMat{<:Number}

Read a given data block from an HDF5 group, applying the correct units.

# Arguments

  - `block::String`: The name of the target data block. The possibilities are the keys of [`QUANTITIES`](@ref).
  - `group::HDF5.Group`: The HDF5 group containing the data block.
  - `unit::Union{Unitful.Quantity,Unitful.Units}=Unitful.NoUnits`: The unit of the data block according to [`internalUnits`](@ref).
  - `mmap::Bool=false`: Whether to use memory-mapping for reading the data block.

# Returns

  - The data block as an array, with the correct units applied.
"""
function readBlock(
    block::String,
    group::HDF5.Group;
    unit::Union{Unitful.Quantity,Unitful.Units}=Unitful.NoUnits,
    mmap::Bool=false,
)::VecOrMat{<:Number}

    (
        isBlockPresent(block, group) ||
        throw(ArgumentError("readBlock: The block $(block) is missing from the HDF5 group"))
    )

    dataset = group[QUANTITIES[block].hdf5_name]

    # Bulk read the entire dataset into RAM
    raw_data = read(dataset)

    if mmap

        if unit == Unitful.NoUnits

            data_array = allocateArray(eltype(dataset), size(dataset); mmap=true)

            # Directly copy the HDF5 dataset to the preallocated buffer
            data_array .= raw_data

            return data_array

        end

        data_type  = typeof(one(eltype(dataset)) * unit)
        data_array = allocateArray(data_type, size(dataset); mmap=true)

        data_array .= raw_data .* unit

        return data_array

    else

        if unit == Unitful.NoUnits
            return raw_data
        end

        return raw_data .* unit

    end

end

"""
    readBlock(
        block::String,
        group::HDF5.Group,
        mask::Union{Vector{Bool},Colon};
        <keyword arguments>
    )::VecOrMat{<:Number}

Read a given data block from an HDF5 group, applying the correct units.

# Arguments

  - `block::String`: The name of the target data block. The possibilities are the keys of [`QUANTITIES`](@ref).
  - `group::HDF5.Group`: The HDF5 group containing the data block.
  - `mask::Union{Vector{Bool},Colon}`: A boolean vector to select only a subset of the data along its last dimension. It must have the same length as the last dimension of the data block. If `mask` is `(:)`, the whole data block will be read.
  - `unit::Union{Unitful.Quantity,Unitful.Units}=Unitful.NoUnits`: The unit of the data block according to [`internalUnits`](@ref).
  - `mmap::Bool=false`: Whether to use memory-mapping for reading the data block.

# Returns

  - The data block as an array, with the correct units applied.
"""
function readBlock(
    block::String,
    group::HDF5.Group,
    mask::Union{Vector{Bool},Colon};
    unit::Union{Unitful.Quantity,Unitful.Units}=Unitful.NoUnits,
    mmap::Bool=false,
)::VecOrMat{<:Number}

    if mask isa Colon
        return readBlock(block, group; unit, mmap)
    end

    (
        isBlockPresent(block, group) ||
        throw(ArgumentError("readBlock: The block $(block) is missing from the HDF5 group"))
    )

    dataset = group[QUANTITIES[block].hdf5_name]
    dims    = size(dataset)
    nd      = length(dims)

    last_dim_size = dims[end]
    (
        length(mask) == last_dim_size ||
        throw(DimensionMismatch("readBlock: mask has length $(length(mask)) but the last dimension \
        of the data in the HDF5 file has length $(last_dim_size)"))
    )

    # Bulk read the entire dataset into RAM
    raw_data = read(dataset)

    if mmap

        num_selected = count(mask)
        outdims = nd == 1 ? (num_selected,) : (dims[1], num_selected)

        if unit == Unitful.NoUnits

            data_array = allocateArray(eltype(dataset), outdims; mmap=true)

            mask_idxs = nd == 1 ? (mask,) : (:, mask)

            data_array .= getindex(raw_data, mask_idxs...)

            return data_array

        end

        data_type  = typeof(one(eltype(dataset)) * unit)
        data_array = allocateArray(data_type, outdims; mmap=true)

        data_array .= selectdim(raw_data, nd, mask) .* unit

        return data_array

    else

        if unit == Unitful.NoUnits

            mask_idxs = nd == 1 ? (mask,) : (:, mask)

            return getindex(raw_data, mask_idxs...)

        end

        return selectdim(raw_data, nd, mask) .* unit

    end

end

"""
    findRealStars(snapshot::HDF5.File)::Vector{Bool}

Find which stellar particles are real stars and not wind particles.

# Arguments

  - `snapshot::HDF5.File`: The HDF5 file, as returned by `h5open` (see [HDF5](https://juliaio.github.io/HDF5.jl/stable/)).

# Returns

  - A boolean vector with true for stars and false for wind particles.
"""
function findRealStars(snapshot::HDF5.File)::Vector{Bool}

    if haskey(snapshot, PARTICLE_CODE_NAME[:stellar])

        group = snapshot[PARTICLE_CODE_NAME[:stellar]]

        if isBlockPresent("GAGE", group)

            time_of_birth = readBlock("GAGE", group)

            return map(isPositive, time_of_birth)

        elseif haskey(snapshot, "Header")

            # If there is no age data, we assume all the stars reported in the header are real
            # stars (i.e. no wind particles)
            num_part_total = read_attribute(snapshot["Header"], "NumPart_Total")
            num_stars = num_part_total[PARTICLE_INDEX[:stellar] + 1]

            (
                logging[] &&
                @warn("findRealStars: The snapshot $(snapshot) does not have age data, I will \
                assume that all stars are real")
            )

            return fill(true, num_stars)

        else

            (
                logging[] &&
                @warn("findRealStars: The snapshot $(snapshot) does not have age data, nor a \
                header, I will assume that no stars are real")
            )

            return Bool[]

        end

    else

        (
            logging[] &&
            @warn("findRealStars: The snapshot $(snapshot) does not have stellar particles")
        )

        return Bool[]

    end

end

"""
    findRealStars(path::String)::Vector{Bool}

Find which stellar particles are real stars and not wind particles.

# Arguments

  - `path::String`: Path to the snapshot file or folder.

# Returns

  - A boolean vector with true for stars and false for wind particles.
"""
function findRealStars(path::String)::Vector{Bool}

    if isfile(path)

        (
            HDF5.ishdf5(path) ||
            throw(ArgumentError("findRealStars: The file $(path) is not in the HDF5 format, \
            I don't know how to read it"))
        )

        stellar_mask = h5open(path, "r") do snapshot

            findRealStars(snapshot)

        end

        return stellar_mask

    elseif isdir(path)

        sub_files = findFiles(path, "$(SNAP_BASENAME)_*.*.hdf5")

        (
            !isempty(sub_files) && all(HDF5.ishdf5, sub_files) ||
            throw(ArgumentError("findRealStars: The directory $(path) does not contain \
            snapshot sub-files in the HDF5 format"))
        )

        return mapreduce(findRealStars, vcat, sub_files)

    else

        throw(ArgumentError("findRealStars: $(path) does not exists as a file or folder"))

    end

end

"""
    countStars(path::String)::Int

Count the number of stars in a snapshot, excluding wind particles.

# Arguments

  - `path::String`: Path to the snapshot file or folder.

# Returns

  - The number of stars.
"""
countStars(path::String)::Int = count(findRealStars(path))

"""
    readTime(path::String)::Float64

Read the "Time" field in the header of a snapshot file.

!!! note

    If each snapshot is made of multiple files, this method will only read the first one.

# Arguments

  - `path::String`: Path to the snapshot file or folder.

# Returns

  - The "Time" field in the header (for cosmological simulations it is the scale factor).
"""
function readTime(path::String)::Float64

    if isfile(path)

        (
            HDF5.ishdf5(path) ||
            throw(ArgumentError("readTime: The file $(path) is not in the HDF5 format, \
            I don't know how to read it"))
        )

        file_path = path

    elseif isdir(path)

        sub_files = findFiles(path, "$(SNAP_BASENAME)_*.*.hdf5")

        (
            !isempty(sub_files) && all(HDF5.ishdf5, sub_files) ||
            throw(ArgumentError("readTime: The directory $(path) does not contain \
            snapshot sub-files in the HDF5 format"))
        )

        file_path = first(sub_files)

    else

        throw(ArgumentError("readTime: $(path) does not exists as a file or folder"))

    end

    time = h5open(file_path, "r") do file
        read_attribute(file["Header"], "Time")
    end

    return time

end

"""
    readTemperature(
        group::HDF5.Group,
        file_path::String;
        <keyword arguments>
    )::Vector{<:Unitful.Temperature}

Compute the temperature of the gas cells in a snapshot.

# Arguments

  - `group::HDF5.Group`: The HDF5 group containing the gas particle data.
  - `file_path::String`: Path to the snapshot file.
  - `mmap::Bool=false`: Whether to use memory-mapping for reading the data blocks.

# Returns

  - The temperature of the gas cells.
"""
function readTemperature(
    group::HDF5.Group,
    file_path::String;
    mmap::Bool=false,
)::Vector{<:Unitful.Temperature}

    if isfile(file_path)

        (
            HDF5.ishdf5(file_path) ||
            throw(ArgumentError("readTemperature: The file $(file_path) is not in the \
            HDF5 format, I don't know how to read it"))
        )

    else

        throw(ArgumentError("readTemperature: $(file_path) does not exists as a file"))

    end

    # List of blocks needed to compute the temperature
    blocks = ["U   ", "NE  "]

    for block in blocks
        isBlockPresent(block, group) || throw(ArgumentError("readTemperature: The block \
        $(block) is missing, and I need it to compute the temperature"))
    end

    units = internalUnits.(blocks, file_path)

    u  = readBlock("U   ", group; unit=units[1], mmap)
    ne = readBlock("NE  ", group; unit=units[2], mmap)

    return computeTemperature(u, ne)

end

"""
    countSnapshot(simulation_path::String)::Int

Count the number of snapshots in `simulation_path`.

!!! note

    This function counts the number of snapshots, not the number of snapshot files. So if each snapshot is made of more than one files, the count will not change.

# Arguments

  - `simulation_path::String`: Path to the simulation directory, set in the code variable `OutputDir`.

# Returns

  - The number of snapshots.
"""
function countSnapshot(simulation_path::String)::Int

    (
        isdir(simulation_path) ||
        throw(ArgumentError("countSnapshot: $(simulation_path) does not exists as a directory"))
    )

    # Get the full list of paths to every snapshot in `simulation_path`
    path_list = findFiles(simulation_path, "**/$(SNAP_BASENAME)_*.hdf5")

    # Check for an empty folder
    if isempty(path_list)

        (
            logging[] &&
            @warn("countSnapshot: I could not find any file named $(SNAP_BASENAME)_*.hdf5 \
            within $(simulation_path), or any of its subfolders")
        )

        return 0

    end

    h5open(first(path_list), "r") do snap_file

        if read_attribute(snap_file["Header"], "NumFilesPerSnapshot") > 1
            # If there are multiple files per snapshot, get the path to the snapshot directory
            map!(dirname, path_list)
            # Delete duplicates
            unique!(path_list)
        end

    end

    return length(path_list)

end

"""
    mergeRequests(requests...)::Dict{Symbol,Vector{String}}

Merge several request dictionaries, ignoring duplicates.

# Arguments

  - `requests`: The request dictionaries for [`readSnapshot`](@ref).

# Returns

  - A new dictionary with all the requests.
"""
function mergeRequests(requests...)::Dict{Symbol,Vector{String}}

    return Dict(
        type => union([get(request, type, String[]) for request in requests]...) for
        type in union(keys.(requests)...)
    )

end

"""
    addRequest(
        request::Dict{Symbol,Vector{String}},
        addition::Dict{Symbol,Vector{String}},
    )::Dict{Symbol,Vector{String}}

Add the blocks in `addition` to `request`, only for the types already present in `request`.

# Arguments

  - `request::Dict{Symbol,Vector{String}}`: The request dictionary for [`readSnapshot`](@ref).
  - `addition::Dict{Symbol,Vector{String}}`: Request dictionary with the blocks to be added, only for the types already present in `request`.

# Returns

  - A new dictionary with all the requests.
"""
function addRequest(
    request::Dict{Symbol,Vector{String}},
    addition::Dict{Symbol,Vector{String}},
)::Dict{Symbol,Vector{String}}

    return Dict(type => blocks ∪ get(addition, type, String[]) for (type, blocks) in request)

end

"""
    isSubfindActive(path::String)::Bool

Check if the group catalog file has information or is empty.

# Arguments

  - `path::String`: Path to the group catalog file or folder.

# Returns

  - If there are halo and subhalo information in the group catalog file.
"""
function isSubfindActive(path::String)::Bool

    if isfile(path)

        (
            HDF5.ishdf5(path) ||
            throw(ArgumentError("isSubfindActive: The file $(path) is not in the HDF5 format, \
            I don't know how to read it"))
        )

        file_path = path

    elseif isdir(path)

        sub_files = findFiles(path, "$(GC_BASENAME)_*.*.hdf5")

        (
            !isempty(sub_files) && all(HDF5.ishdf5, sub_files) ||
            throw(ArgumentError("isSubfindActive: The directory $(path) does not contain \
            group catalog sub-files in the HDF5 format"))
        )

        file_path = first(sub_files)

    else

        throw(ArgumentError("isSubfindActive: $(path) does not exists as a file or folder"))

    end

    subfind_active = h5open(file_path, "r") do gc_file

        (
            all(in(keys(gc_file)), ["Group", "Subhalo"]) &&
            all(!isempty, [gc_file["Group"], gc_file["Subhalo"]])
        )

    end

    return subfind_active

end

"""
    isSubfindActive(path::Missing)::Bool

Default method of [`isSubfindActive`](@ref) for a missing group catalog file.
"""
isSubfindActive(path::Missing)::Bool = false

"""
    isSnapCosmological(path::String)::Bool

Check if the snapshot in `path` comes from a cosmological simulation.

!!! note

    If each snapshot is made of multiple files, the function will only read the first file.

# Arguments

  - `path::String`: Path to the snapshot file or folder.

# Returns

  - If the simulation is cosmological

      + `false` -> Newtonian simulation    (`ComovingIntegrationOn` = 0, `Redshift` = 0.0).
      + `true`  -> Cosmological simulation (`ComovingIntegrationOn` = 1, `Redshift` != 0.0).
"""
function isSnapCosmological(path::String)::Bool

    if isfile(path)

        (
            HDF5.ishdf5(path) ||
            throw(ArgumentError("isSnapCosmological: The file $(path) is not in the HDF5 format, \
            I don't know how to read it"))
        )

        file_path = path

    elseif isdir(path)

        sub_files = findFiles(path, "$(SNAP_BASENAME)_*.*.hdf5")

        (
            !isempty(sub_files) && all(HDF5.ishdf5, sub_files) ||
            throw(ArgumentError("isSnapCosmological: The directory $(path) does not contain \
            snapshot sub-files in the HDF5 format"))
        )

        file_path = first(sub_files)

    else

        throw(ArgumentError("isSnapCosmological: $(path) does not exists as a file or folder"))

    end

    cosmological = h5open(file_path, "r") do snapshot

        if haskey(snapshot, "Parameters")
            # If the `param.txt` file is saved in the snapshot metadata, read `ComovingIntegrationOn`
            read_attribute(snapshot["Parameters"], "ComovingIntegrationOn")
        else
            # Otherwise, use the redshift in the header
            !iszero(read_attribute(snapshot["Header"], "Redshift"))
        end

    end

    return cosmological

end

"""
    isSimCosmological(simulation_path::String)::Bool

Check if the simulation in `simulation_path` is cosmological.

!!! note

    This function will only read the first snapshot, and if each snapshot is made of multiple files, the function will only read the first file.

# Arguments

  - `simulation_path::String`: Path to the simulation directory, set in the code variable `OutputDir`.

# Returns

  - If the simulation is cosmological

      + `false` -> Newtonian simulation    (`ComovingIntegrationOn` = 0, `Redshift` = 0.0).
      + `true`  -> Cosmological simulation (`ComovingIntegrationOn` = 1, `Redshift` != 0.0).
"""
function isSimCosmological(simulation_path::String)::Bool

    (
        isdir(simulation_path) ||
        throw(ArgumentError("isSimCosmological: $(simulation_path) does not exists as a directory"))
    )

    # Get the full list of paths to every snapshot in `simulation_path`
    path_list = findFiles(simulation_path, "**/$(SNAP_BASENAME)_*.hdf5")

    # Check for an empty folder
    if isempty(path_list)

        (
            logging[] &&
            @warn("isSimCosmological: I could not find any file named $(SNAP_BASENAME)_*.hdf5 \
            within $(simulation_path), or any of its subfolders")
        )

        return false

    end

    return isSnapCosmological(first(path_list))

end

"""
    isSnapSFM(path::String)::Bool

Check if the snapshot in `path` comes from a simulation with our star formation model.

!!! note

    If each snapshot is made of multiple files, the function will only read the first file.

# Arguments

  - `path::String`: Path to the snapshot file or folder.

# Returns

  - If the simulation has our star formation model.
"""
function isSnapSFM(path::String)::Bool

    if isfile(path)

        (
            HDF5.ishdf5(path) ||
            throw(ArgumentError("isSnapSFM: The file $(path) is not in the HDF5 format, \
            I don't know how to read it"))
        )

        file_path = path

    elseif isdir(path)

        sub_files = findFiles(path, "$(SNAP_BASENAME)_*.*.hdf5")

        (
            !isempty(sub_files) && all(HDF5.ishdf5, sub_files) ||
            throw(ArgumentError("isSnapSFM: The directory $(path) does not contain \
            snapshot sub-files in the HDF5 format"))
        )

        file_path = first(sub_files)

    else

        throw(ArgumentError("isSnapSFM: $(path) does not exists as a file or folder"))

    end

    sfm_flag = h5open(file_path, "r") do snapshot

        gas_key = PARTICLE_CODE_NAME[:gas]

        haskey(snapshot, gas_key) || return false

        group = snapshot[gas_key]

        isBlockPresent("FRAC", group) || return false

        dataset = group[QUANTITIES["FRAC"].hdf5_name]

        !isempty(dataset)

    end

    return sfm_flag

end

"""
    isSimSFM(simulation_path::String)::Bool

Check if the simulation in `simulation_path` has our star formation model.

!!! note

    This function will only read the first snapshot, and if each snapshot is made of multiple files, the function will only read the first file.

# Arguments

  - `simulation_path::String`: Path to the simulation directory, set in the code variable `OutputDir`.

# Returns

  - If the simulation has our star formation model.
"""
function isSimSFM(simulation_path::String)::Bool

    (
        isdir(simulation_path) ||
        throw(ArgumentError("isSimSFM: $(simulation_path) does not exists as a directory"))
    )

    # Get the full list of paths to every snapshot in `simulation_path`
    path_list = findFiles(simulation_path, "**/$(SNAP_BASENAME)_*.hdf5")

    # Check for an empty folder
    if isempty(path_list)

        (
            logging[] &&
            @warn("isSimSFM: I could not find any file named $(SNAP_BASENAME)_*.hdf5 \
            within $(simulation_path), or any of its subfolders")
        )

        return false

    end

    return isSnapSFM(first(path_list))

end

"""
    snapshotTypes(data_dict::Dict)::Vector{Symbol}

Find which cell/particle types are part of the keys of `data_dict`.

# Arguments

  - `data_dict::Dict`: Data dictionary (see [`makeDataDict`](@ref) for the canonical description).

# Returns

  - A vector with the cell/particle types.
"""
snapshotTypes(data_dict::Dict)::Vector{Symbol} = collect(keys(PARTICLE_INDEX) ∩ keys(data_dict))

"""
    snapshotTypes(path::String)::Vector{Symbol}

Find which cell/particle types are part of the snapshot in `path`.

!!! note

    If each snapshot is made of multiple files, the function will only check the first file.

# Arguments

  - `path::String`: Path to the snapshot file or folder.

# Returns

  - A vector with the cell/particle types.
"""
function snapshotTypes(path::String)::Vector{Symbol}

    if isfile(path)

        (
            HDF5.ishdf5(path) ||
            throw(ArgumentError("snapshotTypes: The file $(path) is not in the HDF5 format, \
            I don't know how to read it"))
        )

        file_path = path

    elseif isdir(path)

        sub_files = findFiles(path, "$(SNAP_BASENAME)_*.*.hdf5")

        (
            !isempty(sub_files) && all(HDF5.ishdf5, sub_files) ||
            throw(ArgumentError("snapshotTypes: The directory $(path) does not contain \
            snapshot sub-files in the HDF5 format"))
        )

        file_path = first(sub_files)

    else

        throw(ArgumentError("snapshotTypes: $(path) does not exists as a file or folder"))

    end

    snapshot_types = h5open(file_path, "r") do snapshot
        collect(keys(PARTICLE_TYPE) ∩ keys(snapshot))
    end

    return snapshot_types

end

"""
    groupCatTypes(data_dict::Dict)::Vector{Symbol}

Find which group catalog data types are part of the keys of `data_dict`.

# Arguments

  - `data_dict::Dict`: Data dictionary (see [`makeDataDict`](@ref) for the canonical description).

# Returns

  - A vector with the group catalog data types.
"""
groupCatTypes(data_dict::Dict)::Vector{Symbol} = [:group, :subhalo] ∩ keys(data_dict)

"""
    groupCatTypes(path::String)::Vector{Symbol}

Find which group catalog data types are part of the group catalog files in `path`.

!!! note

    If each group catalog is made of multiple files, the function will only check the first file.

# Arguments

  - `path::String`: Path to the group catalog file or folder.

# Returns

  - A vector with the group catalog data types.
"""
function groupCatTypes(path::String)::Vector{Symbol}

    if isfile(path)

        (
            HDF5.ishdf5(path) ||
            throw(ArgumentError("groupCatTypes: The file $(path) is not in the HDF5 format, \
            I don't know how to read it"))
        )

        file_path = path

    elseif isdir(path)

        sub_files = findFiles(path, "$(GC_BASENAME)_*.*.hdf5")

        (
            !isempty(sub_files) && all(HDF5.ishdf5, sub_files) ||
            throw(ArgumentError("groupCatTypes: The directory $(path) does not contain \
            group catalog sub-files in the HDF5 format"))
        )

        file_path = first(sub_files)

    else

        throw(ArgumentError("groupCatTypes: $(path) does not exists as a file or folder"))

    end

    groupcat_types = h5open(file_path, "r") do gc_file
        [:group, :subhalo] ∩ keys(gc_file)
    end

    return groupcat_types

end

"""
    readGroupCatHeader(path::Union{String,Missing})::GroupCatHeader

Read the header of a group catalog file in the HDF5 format.

!!! note

    If each group catalog is made of multiple files, this method will read the header of the first one.

# Arguments

  - `path::Union{String,Missing}`: Path to the group catalog file or folder.

# Returns

  - A [`GroupCatHeader`](@ref).
"""
function readGroupCatHeader(path::Union{String,Missing})::GroupCatHeader

    if ismissing(path)

        logging[] && @warn("readGroupCatHeader: The group catalog file or folder is missing")

        return GroupCatHeader()

    elseif isfile(path)

        (
            HDF5.ishdf5(path) ||
            throw(ArgumentError("readGroupCatHeader: The file $(path) is not in the \
            HDF5 format, I don't know how to read it"))
        )

        file_path = path

    elseif isdir(path)

        sub_files = findFiles(path, "$(GC_BASENAME)_*.*.hdf5")

        (
            !isempty(sub_files) && all(HDF5.ishdf5, sub_files) ||
            throw(ArgumentError("readGroupCatHeader: The directory $(path) does not contain \
            group catalog sub-files in the HDF5 format"))
        )

        # The header should be in each subfile, so we can read it from the first one
        file_path = first(sub_files)

    else

        throw(ArgumentError("readGroupCatHeader: $(path) does not exists as a file or folder"))

    end

    header = h5open(file_path, "r") do gc_file

        head = gc_file["Header"]

        attrs_present = keys(HDF5.attrs(head))

        missing_attrs = setdiff(
            [
                "BoxSize",
                "HubbleParam",
                "Ngroups_ThisFile",
                "Ngroups_Total",
                "Nsubgroups_ThisFile",
                "Nsubgroups_Total",
                "NumFiles",
                "Omega0",
                "OmegaLambda",
                "Redshift",
                "Time",
            ],
            attrs_present,
        )

        (
            isempty(missing_attrs) ||
            throw(ArgumentError("readGroupCatHeader: The attributes $(missing_attrs) are missing \
            from the header"))
        )

        GroupCatHeader(
            box_size          = read_attribute(head, "BoxSize"),
            h                 = read_attribute(head, "HubbleParam"),
            n_groups_part     = read_attribute(head, "Ngroups_ThisFile"),
            n_groups_total    = read_attribute(head, "Ngroups_Total"),
            n_subgroups_part  = read_attribute(head, "Nsubgroups_ThisFile"),
            n_subgroups_total = read_attribute(head, "Nsubgroups_Total"),
            num_files         = read_attribute(head, "NumFiles"),
            omega_0           = read_attribute(head, "Omega0"),
            omega_l           = read_attribute(head, "OmegaLambda"),
            redshift          = read_attribute(head, "Redshift"),
            time              = read_attribute(head, "Time"),
        )

    end

    return header

end

"""
    readSnapHeader(path::String)::SnapshotHeader

Read the header of a snapshot file in the HDF5 format.

!!! note

    If each snapshot is made of multiple files, this method will read the header of the first one.

# Arguments

  - `path::String`: Path to the snapshot file or folder.

# Returns

  - A [`SnapshotHeader`](@ref) structure.
"""
function readSnapHeader(path::String)::SnapshotHeader

    if isfile(path)

        (
            HDF5.ishdf5(path) ||
            throw(ArgumentError("readSnapHeader: The file $(path) is not in the HDF5 format, \
            I don't know how to read it"))
        )

        file_path = path

        num_files = h5open(file_path, "r") do snap_file
            read_attribute(snap_file["Header"], "NumFilesPerSnapshot")
        end

        # Count the number of stellar particles, ignoring wind particles
        if num_files > 1

            num_part_stars = countStars(file_path)
            num_total_stars = countStars(dirname(file_path))

        else

            num_part_stars = countStars(file_path)
            num_total_stars = num_part_stars

        end

    elseif isdir(path)

        sub_files = findFiles(path, "$(SNAP_BASENAME)_*.*.hdf5")

        (
            !isempty(sub_files) && all(HDF5.ishdf5, sub_files) ||
            throw(ArgumentError("readSnapHeader: The directory $(path) does not contain \
            snapshot sub-files in the HDF5 format"))
        )

        # The header should be in each sub-file, so we can read it from the first one
        file_path = first(sub_files)

        # Count the number of stellar particles, ignoring wind particles
        num_part_stars = countStars(file_path)
        num_total_stars = countStars(path)

    else

        throw(ArgumentError("readSnapHeader: $(path) does not exists as a file or folder"))

    end

    header = h5open(file_path, "r") do snap_file

        head = snap_file["Header"]

        attrs_present = keys(HDF5.attrs(head))

        missing_attrs = setdiff(
            [
                "BoxSize",
                "HubbleParam",
                "MassTable",
                "NumFilesPerSnapshot",
                "NumPart_ThisFile",
                "NumPart_Total",
                "Omega0",
                "OmegaLambda",
                "Redshift",
                "Time",
            ],
            attrs_present,
        )

        (
            isempty(missing_attrs) ||
            throw(ArgumentError("readSnapHeader: The attributes $(missing_attrs) are missing from \
            the header"))
        )

        # Only for the stars edit the number in the header, to exclude wind particles
        num_part = read_attribute(head, "NumPart_ThisFile")
        num_part[PARTICLE_INDEX[:stellar] + 1] = num_part_stars
        num_total = read_attribute(head, "NumPart_Total")
        num_total[PARTICLE_INDEX[:stellar] + 1] = num_total_stars

        # Check if the internal length unit is in the header, otherwise use the default value
        if "UnitLength_in_cm" ∈ attrs_present
            l_unit = read_attribute(head, "UnitLength_in_cm") * u"cm"
        else
            l_unit = DEFAULT_L_UNIT[]
        end

        # Check if the internal mass unit is in the header, otherwise use the default value
        if "UnitMass_in_g" ∈ attrs_present
            m_unit = read_attribute(head, "UnitMass_in_g") * u"g"
        else
            m_unit = DEFAULT_M_UNIT[]
        end

        # Check if the internal velocity unit is in the header, otherwise use the default value
        if "UnitVelocity_in_cm_per_s" ∈ attrs_present
            v_unit = read_attribute(head, "UnitVelocity_in_cm_per_s") * u"cm * s^-1"
        else
            v_unit = DEFAULT_V_UNIT[]
        end

        SnapshotHeader(
            box_size   = read_attribute(head, "BoxSize"),
            h          = read_attribute(head, "HubbleParam"),
            mass_table = read_attribute(head, "MassTable"),
            num_files  = read_attribute(head, "NumFilesPerSnapshot"),
            num_part   = num_part,
            num_total  = num_total,
            omega_0    = read_attribute(head, "Omega0"),
            omega_l    = read_attribute(head, "OmegaLambda"),
            redshift   = read_attribute(head, "Redshift"),
            time       = read_attribute(head, "Time"),
            l_unit     = l_unit,
            m_unit     = m_unit,
            v_unit     = v_unit,
        )

    end

    return header

end

"""
    readSnapReducedHeader(path::String)::SnapshotHeader

Read the header of a snapshot file in the HDF5 format.

!!! note

    This method does not read the number of cell/particles of each type, which is an expensive operation that uses [`countStars`](@ref). Use it when that info is not needed.

!!! note

    If each snapshot is made of multiple files, this method will read the header of the first one.

# Arguments

  - `path::String`: Path to the snapshot file or folder.

# Returns

  - A [`SnapshotHeader`](@ref) structure.
"""
function readSnapReducedHeader(path::String)::SnapshotHeader

    if isfile(path)

        (
            HDF5.ishdf5(path) ||
            throw(ArgumentError("readSnapReducedHeader: The file $(path) is not in the HDF5 \
            format, I don't know how to read it"))
        )

        file_path = path

    elseif isdir(path)

        sub_files = findFiles(path, "$(SNAP_BASENAME)_*.*.hdf5")

        (
            !isempty(sub_files) && all(HDF5.ishdf5, sub_files) ||
            throw(ArgumentError("readSnapReducedHeader: The directory $(path) does not contain \
            snapshot sub-files in the HDF5 format"))
        )

        # The header should be in each sub-file, so we can read it from the first one
        file_path = first(sub_files)

    else

        throw(ArgumentError("readSnapReducedHeader: $(path) does not exists as a file or folder"))

    end

    header = h5open(file_path, "r") do snap_file

        head = snap_file["Header"]

        attrs_present = keys(HDF5.attrs(head))

        missing_attrs = setdiff(
            [
                "BoxSize",
                "HubbleParam",
                "MassTable",
                "NumFilesPerSnapshot",
                "Omega0",
                "OmegaLambda",
                "Redshift",
                "Time",
            ],
            attrs_present,
        )

        (
            isempty(missing_attrs) ||
            throw(ArgumentError("readSnapReducedHeader: The attributes $(missing_attrs) are \
            missing from the header"))
        )

        # Check if the internal length unit is in the header, otherwise use the default value
        if "UnitLength_in_cm" ∈ attrs_present
            l_unit = read_attribute(head, "UnitLength_in_cm") * u"cm"
        else
            l_unit = DEFAULT_L_UNIT[]
        end

        # Check if the internal mass unit is in the header, otherwise use the default value
        if "UnitMass_in_g" ∈ attrs_present
            m_unit = read_attribute(head, "UnitMass_in_g") * u"g"
        else
            m_unit = DEFAULT_M_UNIT[]
        end

        # Check if the internal velocity unit is in the header, otherwise use the default value
        if "UnitVelocity_in_cm_per_s" ∈ attrs_present
            v_unit = read_attribute(head, "UnitVelocity_in_cm_per_s") * u"cm * s^-1"
        else
            v_unit = DEFAULT_V_UNIT[]
        end

        SnapshotHeader(
            box_size   = read_attribute(head, "BoxSize"),
            h          = read_attribute(head, "HubbleParam"),
            mass_table = read_attribute(head, "MassTable"),
            num_files  = read_attribute(head, "NumFilesPerSnapshot"),
            num_part   = [0],
            num_total  = [0],
            omega_0    = read_attribute(head, "Omega0"),
            omega_l    = read_attribute(head, "OmegaLambda"),
            redshift   = read_attribute(head, "Redshift"),
            time       = read_attribute(head, "Time"),
            l_unit     = l_unit,
            m_unit     = m_unit,
            v_unit     = v_unit,
        )

    end

    return header

end

"""
    readSfrFile(
        file_path::String,
        snap_path::String,
    )::DataFrame

Read the `sfr.txt` file.

# Arguments

  - `file_path::String`: Path to the `sfr.txt` file.
  - `snapshot_path::String`: Path to one snapshot file or folder of the simulation. This is needed for unit conversion.

# Returns

  - A DataFrame with the following columns:

      + 1. Time or scale factor (internal units).
      + 2. Total stellar mass to be formed prior to stochastic sampling (internal units).
      + 3. Instantaneous star formation rate of all cells (``\\mathrm{M_\\odot \\, yr^{-1}}``).
      + 4. Instantaneous star formation rate of active cells (``\\mathrm{M_\\odot \\, yr^{-1}}``).
      + 5. Total mass in stars formed after stochastic sampling (internal units).
      + 6. Cumulative stellar mass formed (internal units).
"""
function readSfrFile(file_path::String, snap_path::String)::DataFrame

    isfile(file_path) || throw(ArgumentError("readSfrFile: $(file_path) does not exists as a file"))

    # Load the data from the `sfr.txt` file
    file_data = CSV.read(file_path, DataFrame; header=false, delim=' ', ignorerepeated=true)
    n_cols = size(file_data, 2)

    # Check that the data in the file has the correct size
    (
        n_cols <= 6 ||
        throw(ArgumentError("readSfrFile: I don't know how to handle more \
        than 6 columns in `sfr.txt`"))
    )

    (
        logging[] && n_cols < 6 &&
        @warn("readSfrFile: I could only find $(n_cols) columns \
        in $(file_path). I was expecting 6")
    )

    # Load the units for each column
    units = [internalUnits("SFC$(i)", snap_path) for i in 1:n_cols]

    for (i, col) in enumerate(names(file_data))
        file_data[!, col] .*= units[i]
    end

    return file_data

end

"""
    readCpuFile(
        file_path::String,
        targets::Vector{String};
        <keyword arguments>
    )::Dict{String,Matrix{Float64}}

Read the `cpu.txt` file.

For each process in `targets` a matrix with all the CPU usage data is returned.

# Arguments

  - `file_path::String`: Path to the `cpu.txt` file.
  - `targets::Vector{String}`: Target processes (e.g. "total").
  - `step::Int=1`: Step used to traverse the rows.

# Returns

  - A dictionary with the following shape:

    `target process` -> matrix with columns:

      + Time step.
      + Simulation time (scale factor for cosmological simulations and physical time for non-cosmological simulations).
      + Clock time in seconds.
      + Clock time as a percentage.
      + Total clock time in seconds.
      + Total clock time as a percentage.
"""
function readCpuFile(
    file_path::String,
    targets::Vector{String};
    step::Int=1,
)::Dict{String,Matrix{Float64}}

    isfile(file_path) || throw(ArgumentError("readCpuFile: $(file_path) does not exists as a file"))

    # Load the data from the `cpu.txt` file
    file_data = eachline(file_path)

    # Set up an auxiliary dictionary
    data_aux = Dict(target => Matrix{Float64}[] for target in targets)

    # Clock time for each sync-point
    time = 0.0
    # Time step
    time_step = 0

    for line in file_data

        # Ignore empty lines
        !isempty(line) || continue

        columns = split(line)
        title = columns[1]

        # Ignore header lines
        !(title == "diff") || continue

        # Use "Step" lines to capture the clock time
        if title == "Step"
            time = parse(Float64, rstrip(columns[4], ','))
            time_step = parse(Int, rstrip(columns[2], ','))
            continue
        end

        if title ∈ targets
            push!(
                data_aux[title],
                [
                    time_step;;                               # Time step
                    time;;                                    # Simulation time for each sync-point
                    parse(Float64, columns[2]);;              # Clock time in seconds
                    parse(Float64, rstrip(columns[3], '%'));; # Clock time as a percentage
                    parse(Float64, columns[4]);;              # Total clock time in seconds
                    parse(Float64, rstrip(columns[5], '%'))   # Total clock time as a percentage
                ],
            )
        end

    end

    (
        !(logging[] && any(isempty, values(data_aux))) ||
        @warn("readCpuFile: I could not find some of the target rows in $(file_path)")
    )

    data_out = Dict{String,Matrix{Float64}}()

    # Try reducing the data size
    for (target, values) in data_aux

        # Ignore empty targets
        !isempty(values) || continue

        l_e = length(values)

        if 1 < step < l_e
            data_out[target] = vcat(values[1:step:end]...)
        else
            data_out[target] = vcat(values...)
            (
                !(logging[] && step > l_e) ||
                @warn("readCpuFile: `step` = $(step) is bigger than the number \
                of time steps in $(file_path)")
            )
        end

    end

    return data_out

end

"""
    getSnapshotPaths(simulation_path::String)::Dict{Symbol,Vector{String}}

Find the path and number of every snapshot in `simulation_path`.

!!! note

    If each snapshot is made of multiple files, the `:paths` field will have paths to folders, each one containing the sub-files of the corresponding snapshot.

# Arguments

  - `simulation_path::String`: Path to the simulation directory, set in the code variable `OutputDir`.

# Returns

  - A dictionary with the following shape:

      + `:numbers` -> The number that characterize each snapshot.
      + `:paths`   -> The full path to each snapshot.
"""
function getSnapshotPaths(simulation_path::String)::Dict{Symbol,Vector{String}}

    (
        isdir(simulation_path) ||
        throw(ArgumentError("getSnapshotPaths: $(simulation_path) does not exists as a directory"))
    )

    # Get the full list of paths to every snapshot in `simulation_path`
    path_list = findFiles(simulation_path, "**/$(SNAP_BASENAME)_*.hdf5")

    # Check for an empty folder
    if isempty(path_list)

        (
            logging[] &&
            @warn("getSnapshotPaths: I could not find any file named $(SNAP_BASENAME)_*.hdf5 \
            within $(simulation_path), or any of its subfolders")
        )

        return Dict(:numbers => String[], :paths => String[])

    end

    # Get the numbers that characterize each snapshot
    reg = Regex("(?<=$(SNAP_BASENAME)_).*?(?=(?:\\.)|\$)")
    number_list = map(x -> match(reg, x).match, path_list)

    h5open(first(path_list), "r") do snap_file

        if read_attribute(snap_file["Header"], "NumFilesPerSnapshot") > 1
            # If there are multiple files per snapshot, get the path to the snapshot directory
            map!(dirname, path_list)
            # Delete duplicates
            unique!(path_list)
            unique!(number_list)
        end

    end

    return Dict(:numbers => number_list, :paths => sort(normpath.(path_list)))

end

"""
    getGroupCatPaths(simulation_path::String)::Dict{Symbol,Vector{String}}

Find the path and number of every group catalog in `simulation_path`.

!!! note

    If each group catalog is made of multiple files, the `:paths` field will have paths to folders, each one containing the sub-files of the corresponding group catalog.

# Arguments

  - `simulation_path::String`: Path to the simulation directory, set in the code variable `OutputDir`.

# Returns

  - A dictionary with the following shape:

      + `:numbers` -> The number that characterize each group catalog.
      + `:paths`   -> The full path to each group catalog.
"""
function getGroupCatPaths(simulation_path::String)::Dict{Symbol,Vector{String}}

    (
        isdir(simulation_path) ||
        throw(ArgumentError("getGroupCatPaths: $(simulation_path) does not exists as a directory"))
    )

    # Get the full list of paths to every group catalog in `simulation_path`
    path_list = findFiles(simulation_path, "**/$(GC_BASENAME)_*.hdf5")

    # Check for an empty folder
    if isempty(path_list)

        (
            logging[] &&
            @warn("getGroupCatPaths: I could not find any file named $(GC_BASENAME)_*.hdf5 \
            within $(simulation_path), or any of its subfolders")
        )

        return Dict(:numbers => String[], :paths => String[])

    end

    # Get the numbers that characterize each group catalog
    reg = Regex("(?<=$(GC_BASENAME)_).*?(?=(?:\\.)|\$)")
    number_list = map(x -> match(reg, x).match, path_list)

    if readGroupCatHeader(first(path_list)).num_files > 1
        # If there are multiple files per group catalog, get the path to the group catalog directory
        map!(dirname, path_list)
        # Delete duplicates
        unique!(path_list)
        unique!(number_list)
    end

    return Dict(:numbers => number_list, :paths => sort(normpath.(path_list)))

end

"""
    readGroupCatBlocks(
        file_path::String,
        snapshot_path::String,
        request::Dict{Symbol,Vector{String}};
        <keyword arguments>
    )::Dict{Symbol,Dict{String,VecOrMat{<:Number}}}

Read the specified blocks from a group catalog file.

# Arguments

  - `file_path::String`: Path to the group catalog file.
  - `snapshot_path::String`: Path to the corresponding snapshot file or folder. This is needed for unit conversion.
  - `request::Dict{Symbol,Vector{String}}`: The blocks to be read. It must have the shape `group type` -> [`block`, `block`, `block`], where the possible types are :group and :subhalo, and the possible blocks are the keys of [`QUANTITIES`](@ref).
  - `mmap::Bool=false`: Whether to use memory-mapping for reading the data blocks.

# Returns

  - A dictionary with the following shape: `group type` -> (`block` -> data).
"""
function readGroupCatBlocks(
    file_path::String,
    snapshot_path::String,
    request::Dict{Symbol,Vector{String}};
    mmap::Bool=false,
)::Dict{Symbol,Dict{String,VecOrMat{<:Number}}}

    if isfile(file_path)

        (
            HDF5.ishdf5(file_path) ||
            throw(ArgumentError("readGroupCatBlocks: The file $(file_path) is not in the \
            HDF5 format, I don't know how to read it"))
        )

    else

        throw(ArgumentError("readGroupCatBlocks: $(file_path) does not exists as a file"))

    end

    output = Dict{Symbol,Dict{String,VecOrMat{<:Number}}}()

    h5open(file_path, "r") do gc_file

        # Read from the request only the group catalog types
        for component in groupCatTypes(request)

            blocks = request[component]

            type_str = titlecase(string(component))

            qty_data = Dict{String,VecOrMat{<:Number}}()

            if haskey(gc_file, type_str)

                # Read the HDF5 group
                hdf5_group = gc_file[type_str]

                if isempty(hdf5_group)

                    (
                        logging[] &&
                        @warn("readGroupCatBlocks: The group catalog type :$(component) \
                        in $(file_path) is empty")
                    )

                    # Return an empty array for every missing block
                    for block in blocks
                        unit = internalUnits(block, snapshot_path)
                        qty_data[block] = typeof(1.0 * unit)[]
                    end

                else

                    for block in blocks

                        unit = internalUnits(block, snapshot_path)

                        if isBlockPresent(block, hdf5_group)

                            qty_data[block] = readBlock(block, hdf5_group; unit, mmap)

                        else

                            (
                                logging[] &&
                                @warn("readGroupCatBlocks: The block $(block) for the group \
                                catalog type :$(component) in $(file_path) is missing")
                            )

                            # Return an empty array for every missing block
                            qty_data[block] = typeof(1.0 * unit)[]

                        end

                    end

                end

            else

                (
                    logging[] &&
                    @warn("readGroupCatBlocks: The group catalog type :$(component) in $(file_path) \
                    is missing")
                )

                # Return an empty array for every missing block
                for block in blocks
                    unit = internalUnits(block, snapshot_path)
                    qty_data[block] = typeof(1.0 * unit)[]
                end

            end

            output[component] = qty_data

        end

    end

    return output

end

"""
    readGroupCatBlocks(
        file_paths::Vector{String},
        snapshot_path::String,
        request::Dict{Symbol,Vector{String}};
        <keyword arguments>
    )::Dict{Symbol,Dict{String,VecOrMat{<:Number}}}

Read the specified blocks from several group catalog files.

# Arguments

  - `file_paths::Vector{String},`: Path to the group catalog files.
  - `snapshot_path::String`: Path to the corresponding snapshot file or folder. This is needed for unit conversion.
  - `request::Dict{Symbol,Vector{String}}`: The blocks to be read. It must have the shape `group type` -> [`block`, `block`, `block`], where the possible types are :group and :subhalo, and the possible blocks are the keys of [`QUANTITIES`](@ref).
  - `mmap::Bool=false`: Whether to use memory-mapping for reading the data blocks.

# Returns

  - A dictionary with the following shape: `group type` -> (`block` -> data).
"""
function readGroupCatBlocks(
    file_paths::Vector{String},
    snapshot_path::String,
    request::Dict{Symbol,Vector{String}};
    mmap::Bool=false,
)::Dict{Symbol,Dict{String,VecOrMat{<:Number}}}

    output = Dict{Symbol,Dict{String,VecOrMat{<:Number}}}()

    # Read all the headers to compute total offsets
    headers = readGroupCatHeader.(file_paths)

    first_file = file_paths[1]

    for component in groupCatTypes(request)

        # Flag for empty components
        empty_component = false

        blocks = request[component]
        qty_data = Dict{String,VecOrMat{<:Number}}()

        # Compute the total number of groups/subhalos for this component across all files
        if component == :group
            total_cp = headers[1].n_groups_total
        else
            total_cp = headers[1].n_subgroups_total
        end

        if total_cp == 0

            logging[] && @warn("readGroupCatBlocks: There are 0 :$(component) in $(first_file)")

            for block in blocks
                unit = internalUnits(block, first_file)
                qty_data[block] = typeof(1.0 * unit)[]
            end

            output[component] = qty_data

            continue

        end

        # Pre-allocate buffers for each block in this component
        h5open(first_file, "r") do gc_file

            type_str = titlecase(string(component))

            if haskey(gc_file, type_str)

                # Read the HDF5 group
                hdf5_group = gc_file[type_str]

                for block in blocks

                    unit = internalUnits(block, snapshot_path)

                    if isBlockPresent(block, hdf5_group)

                        block_type = getBlockType(block, hdf5_group)
                        data_type  = typeof(one(block_type) * unit)

                        block_dim = getBlockDims(block, hdf5_group)
                        out_dims  = length(block_dim) == 1 ? (total_cp,) : (block_dim[1], total_cp)

                        qty_data[block] = allocateArray(data_type, out_dims; mmap)

                    else

                        qty_data[block] = typeof(1.0 * unit)[]

                    end

                end

            else

                (
                    logging[] &&
                    @warn("readGroupCatBlocks: The group catalog type :$(component) in $(file_path) \
                    is missing")
                )

                for block in blocks
                    unit = internalUnits(block, first_file)
                    qty_data[block] = typeof(1.0 * unit)[]
                end

                output[component] = qty_data

                empty_component = true

            end

        end

        empty_component && continue

        if component == :group
            num_parts = getfield.(headers, :n_groups_part)
        else
            num_parts = getfield.(headers, :n_subgroups_part)
        end

        # Fill the buffers file by file
        current_offset = 1
        for (i, fp) in enumerate(file_paths)

            n_part = num_parts[i]

            # Skip files with no cells/particles of this type
            n_part == 0 && continue

            # Read this file's chunk
            file_data = readGroupCatBlocks(fp, snapshot_path, Dict(component => blocks); mmap)

            for block in blocks
                chunk = file_data[component][block]

                isempty(chunk) && continue

                if ndims(chunk) == 1
                    qty_data[block][current_offset:(current_offset + n_part - 1)] .= chunk
                else
                    qty_data[block][:, current_offset:(current_offset + n_part - 1)] .= chunk
                end
            end

            current_offset += n_part

        end

        output[component] = qty_data

    end

    return output

end

"""
    readSnapBlocks(
        file_path::String,
        request::Dict{Symbol,Vector{String}};
        <keyword arguments>
    )::Dict{Symbol,Dict{String,VecOrMat{<:Number}}}

Read the specified blocks from a snapshot file.

# Arguments

  - `file_path::String`: Path to the snapshot file.
  - `request::Dict{Symbol,Vector{String}}`: The blocks to be read. It must have the shape `cell/particle type` -> [`block`, `block`, ...], where the possible types are the keys of [`PARTICLE_INDEX`](@ref), and the possible blocks are the keys of [`QUANTITIES`](@ref).
  - `mmap::Bool=false`: Whether to use memory-mapping for reading the data blocks.
  - `header::Union{SnapshotHeader,Nothing}=nothing`: The header can be provided to avoid expensive recalculations of the number of particles.

# Returns

  - A dictionary with the following shape: `cell/particle type` -> (`block` -> data).
"""
function readSnapBlocks(
    file_path::String,
    request::Dict{Symbol,Vector{String}};
    mmap::Bool=false,
    header::Union{SnapshotHeader,Nothing}=nothing,
)::Dict{Symbol,Dict{String,VecOrMat{<:Number}}}

    if isfile(file_path)

        (
            HDF5.ishdf5(file_path) ||
            throw(ArgumentError("readSnapBlocks: The file $(file_path) is not in the \
            HDF5 format, I don't know how to read it"))
        )

    else

        throw(ArgumentError("readSnapBlocks: $(file_path) does not exists as a file"))

    end

    output = Dict{Symbol,Dict{String,VecOrMat{<:Number}}}()

    # Read the header
    if isnothing(header)
        header = readSnapHeader(file_path)
    end

    h5open(file_path, "r") do snapshot

        # Read from the request only the cell/particle types
        for component in snapshotTypes(request)

            blocks = request[component]

            qty_data = Dict{String,VecOrMat{<:Number}}()

            group_name = PARTICLE_CODE_NAME[component]

            if !haskey(snapshot, group_name)

                (
                    logging[] &&
                    @warn("readSnapBlocks: The cell/particle type :$(component) in $(file_path) is \
                    missing")
                )

                for block in blocks
                    unit = internalUnits(block, file_path)
                    qty_data[block] = typeof(1.0 * unit)[]
                end

                output[component] = qty_data

                continue

            end

            # Read the HDF5 group
            hdf5_group = snapshot[group_name]

            if isempty(hdf5_group)

                (
                    logging[] &&
                    @warn("readSnapBlocks: The cell/particle type :$(component) in $(file_path) is \
                    empty")
                )

                # Return an empty array for every missing block
                for block in blocks
                    unit = internalUnits(block, file_path)
                    qty_data[block] = typeof(1.0 * unit)[]
                end

                output[component] = qty_data

                continue

            end

            # For the stars, exclude wind particles
            mask = if component == :stellar
                findRealStars(snapshot)
            else
                (:)
            end

            for block in blocks

                unit = internalUnits(block, file_path)

                if block == "TEMP"

                    (
                        component == :gas ||
                        throw(ArgumentError("readSnapBlocks: I can't compute the \
                        temperature for cells/particles of type :$(component), \
                        only for :gas"))
                    )

                    qty_data["TEMP"] = readTemperature(hdf5_group, file_path; mmap)

                elseif block == "MASS"

                    # Read the mass table from the header
                    mass_in_header = header.mass_table[PARTICLE_INDEX[component] + 1]

                    if iszero(mass_in_header)

                        if isBlockPresent(block, hdf5_group)

                            qty_data["MASS"] = readBlock("MASS", hdf5_group, mask; unit, mmap)

                        else

                            (
                                logging[] &&
                                @warn("readSnapBlocks: The block MASS for the cell/particle type \
                                :$(component) in $(file_path) is missing, and I can't compute it \
                                from the header")
                            )

                            # Return an empty array for the missing block
                            qty_data["MASS"] = typeof(1.0 * unit)[]

                        end

                    else

                        # All cells/particles have the same mass
                        # For stellar particles this value already consideres only real stars,
                        # because we edited the header to exclude wind particles in `readSnapHeader`
                        n_part = Int64(header.num_part[PARTICLE_INDEX[component] + 1])

                        mass_array       = allocateArray(typeof(1.0 * unit), (n_part,); mmap)
                        qty_data["MASS"] = fill!(mass_array, mass_in_header * unit)

                    end

                elseif isBlockPresent(block, hdf5_group)

                    qty_data[block] = readBlock(block, hdf5_group, mask; unit, mmap)

                else

                    (
                        logging[] &&
                        @warn("readSnapBlocks: The block $(block) for the \
                        cell/particle type :$(component) in $(file_path) is missing")
                    )

                    # Return an empty array for every missing block
                    qty_data[block] = typeof(1.0 * unit)[]

                end

            end

            output[component] = qty_data

        end

    end

    return output

end

"""
    readSnapBlocks(
        file_paths::Vector{String},
        request::Dict{Symbol,Vector{String}};
        <keyword arguments>
    )::Dict{Symbol,Dict{String,VecOrMat{<:Number}}}

Read the specified blocks from several snapshot files.

# Arguments

  - `file_paths::Vector{String}`: Paths to the snapshot files.
  - `request::Dict{Symbol,Vector{String}}`: The blocks to be read. It must have the shape `cell/particle type` -> [`block`, `block`, ...], where the possible types are the keys of [`PARTICLE_INDEX`](@ref), and the possible blocks are the keys of [`QUANTITIES`](@ref).
  - `mmap::Bool=false`: Whether to use memory-mapping for reading the data blocks.

# Returns

  - A dictionary with the following shape: `cell/particle type` -> (`block` -> data).
"""
function readSnapBlocks(
    file_paths::Vector{String},
    request::Dict{Symbol,Vector{String}};
    mmap::Bool=false,
)::Dict{Symbol,Dict{String,VecOrMat{<:Number}}}

    output = Dict{Symbol,Dict{String,VecOrMat{<:Number}}}()

    # Read all the headers to compute total offsets
    headers = readSnapHeader.(file_paths)

    first_file = file_paths[1]

    for component in snapshotTypes(request)

        # Flag for empty components
        empty_component = false

        blocks = request[component]
        qty_data = Dict{String,VecOrMat{<:Number}}()

        c_idx = PARTICLE_INDEX[component] + 1

        # Compute the total number of cells/particles for this component across all files
        total_cp = Int64(headers[1].num_total[c_idx])

        if total_cp == 0

            logging[] && @warn("readSnapBlocks: There are 0 cells/particles of type :$(component) \
            in $(first_file)")

            for block in blocks
                unit = internalUnits(block, first_file)
                qty_data[block] = typeof(1.0 * unit)[]
            end

            output[component] = qty_data

            continue

        end

        # Pre-allocate buffers for each block in this component
        h5open(first_file, "r") do snapshot

            group_name = PARTICLE_CODE_NAME[component]

            if haskey(snapshot, group_name)

                hdf5_group = snapshot[group_name]

                for block in blocks

                    unit = internalUnits(block, first_file)

                    if isBlockPresent(block, hdf5_group)

                        block_type = getBlockType(block, hdf5_group)
                        data_type  = typeof(one(block_type) * unit)

                        block_dim = getBlockDims(block, hdf5_group)

                        out_dims = length(block_dim) == 1 ? (total_cp,) : (block_dim[1], total_cp)

                        qty_data[block] = allocateArray(data_type, out_dims; mmap)

                    else

                        qty_data[block] = typeof(1.0 * unit)[]

                    end

                end

            else

                (
                    logging[] &&
                    @warn("readSnapBlocks: The snapshot type :$(component) in $(first_file) \
                    is missing")
                )

                for block in blocks
                    unit = internalUnits(block, first_file)
                    qty_data[block] = typeof(1.0 * unit)[]
                end

                output[component] = qty_data

                empty_component = true

            end

        end

        empty_component && continue

        # Fill the buffers file by file
        current_offset = 1
        for (i, (fp, header)) in enumerate(zip(file_paths, headers))

            n_part = Int64(headers[i].num_part[c_idx])

            # Skip files with no cells/particles of this type
            n_part == 0 && continue

            # Read this file's chunk
            file_data = readSnapBlocks(fp, Dict(component => blocks); mmap, header)

            for block in blocks

                chunk = file_data[component][block]

                isempty(chunk) && continue

                if ndims(chunk) == 1
                    qty_data[block][current_offset:(current_offset + n_part - 1)] .= chunk
                else
                    qty_data[block][:, current_offset:(current_offset + n_part - 1)] .= chunk
                end


            end

            current_offset += n_part

        end

        output[component] = qty_data

    end

    return output

end

"""
    readGroupCatalog(
        path::Union{String,Missing},
        snapshot_path::String,
        request::Dict{Symbol,Vector{String}},
    )::Dict{Symbol,Dict{String,VecOrMat{<:Number}}}

Read the specified blocks from a group catalog file or folder.

# Arguments

  - `path::Union{String,Missing}`: Path to the group catalog file or folder.
  - `snapshot_path::String`: Path to the corresponding snapshot file or folder. This is needed for unit conversion.
  - `request::Dict{Symbol,Vector{String}}`: Which blocks will be read. It must have the shape `group type` -> [`block`, `block`, `block`], where the possible types are :group and :subhalo, and the possible blocks are the keys of [`QUANTITIES`](@ref).

# Returns

  - A dictionary with the following shape: `group type` -> (`block` -> data).
"""
function readGroupCatalog(
    path::Union{String,Missing},
    snapshot_path::String,
    request::Dict{Symbol,Vector{String}},
)::Dict{Symbol,Dict{String,VecOrMat{<:Number}}}

    # Decide if we need to use memory-mapping based on the size of the requested data
    # and the available free physical memory
    mmap = getRequestSize(request, path) > memoryThreshold()

    if ismissing(path)

        logging[] && @warn("readGroupCatalog: The group catalog file or folder is missing")

        return Dict{Symbol,Dict{String,VecOrMat{<:Number}}}()

    elseif isfile(path)

        (
            HDF5.ishdf5(path) ||
            throw(ArgumentError("readGroupCatalog: The file $(path) is not in the HDF5 format, \
            I don't know how to read it"))
        )

        return readGroupCatBlocks(path, snapshot_path, request; mmap)

    elseif isdir(path)

        sub_files = findFiles(path, "$(GC_BASENAME)_*.*.hdf5")

        (
            !isempty(sub_files) && all(HDF5.ishdf5, sub_files) ||
            throw(ArgumentError("readGroupCatalog: The directory $(path) does not contain \
            group catalog sub-files in the HDF5 format"))
        )

        # Sort the sub files to concatenate the data in them correctly
        sort!(sub_files)

        return readGroupCatBlocks(sub_files, snapshot_path, request; mmap)

    else

        throw(ArgumentError("readGroupCatalog: $(path) does not exists as a file or folder"))

    end

end

"""
    readSnapshot(
        path::Union{String,Missing},
        request::Dict{Symbol,Vector{String}},
    )::Dict{Symbol,Dict{String,VecOrMat{<:Number}}}

Read the specified blocks from a snapshot file or folder.

# Arguments

  - `path::Union{String,Missing}`: Path to the snapshot file or folder.
  - `request::Dict{Symbol,Vector{String}}`: Which blocks will be read. It must have the shape `cell/particle type` -> [`block`, `block`, `block`], where the possible types are the keys of [`PARTICLE_INDEX`](@ref), and the possible blocks are the keys of [`QUANTITIES`](@ref).

# Returns

  - A dictionary with the following shape: `cell/particle type` -> (`block` -> data).
"""
function readSnapshot(
    path::Union{String,Missing},
    request::Dict{Symbol,Vector{String}},
)::Dict{Symbol,Dict{String,VecOrMat{<:Number}}}

    # Decide if we need to use memory-mapping based on the size of the requested data
    # and the available free physical memory
    mmap = getRequestSize(request, path) > memoryThreshold()

    if ismissing(path)

        throw(ArgumentError("readSnapshot: The snapshot file or folder is missing"))

    elseif isfile(path)

        (
            HDF5.ishdf5(path) ||
            throw(ArgumentError("readSnapshot: The file $(path) is not in the HDF5 format, \
            I don't know how to read it"))
        )

        return readSnapBlocks(path, request; mmap)

    elseif isdir(path)

        sub_files = findFiles(path, "$(SNAP_BASENAME)_*.*.hdf5")

        (
            !isempty(sub_files) && all(HDF5.ishdf5, sub_files) ||
            throw(ArgumentError("readSnapshot: The directory $(path) does not contain \
            snapshot sub-files in the HDF5 format"))
        )

        # Sort the sub files to concatenate the data in them correctly
        sort!(sub_files)

        return readSnapBlocks(sub_files, request; mmap)

    else

        throw(ArgumentError("readSnapshot: $(path) does not exists as a file or folder"))

    end

end

"""
    getBlock(path::String, component::Symbol, block::String)::VecOrMat{<:Number}

Convenience function to directly get the data associated with one block.

# Arguments

  - `path::String`: Path to the snapshot file or folder.
  - `component::Symbol`: Type of cell/particle. The possibilities are the keys of [`PARTICLE_INDEX`](@ref).
  - `block::String`: Target block. The possibilities are the keys of [`QUANTITIES`](@ref).

# Returns

  - The data for `block`.
"""
function getBlock(path::String, component::Symbol, block::String)::VecOrMat{<:Number}

    return readSnapshot(path, Dict(component => [block]))[component][block]

end

"""
    makeSimulationTable(simulation_path::String)::DataFrame

Construct a dataframe with the path, time stamps, and number of each snapshot and group catalog file in `simulation_path`.

# Arguments

  - `simulation_path::String`: Path to the simulation directory, set in the code variable `OutputDir`.

# Returns

  - A dataframe with 8 columns:

      + `:row_id`         -> Dataframe index of each snapshot, i.e. if there are 10 snapshots in total it runs from 1 to 10.
      + `:numbers`        -> Number in the file name of each snapshot.
      + `:scale_factors`  -> Scale factor of each snapshot.
      + `:redshifts`      -> Redshift of each snapshot.
      + `:physical_times` -> Physical time since the Big Bang of each snapshot.
      + `:lookback_times` -> Physical time left to reach the last snapshot.
      + `:snapshot_paths` -> Full path to the snapshots.
      + `:groupcat_paths` -> Full path to the group catalog files.
"""
function makeSimulationTable(simulation_path::String)::DataFrame

    # Get the path and number of each snapshot
    snap_source = getSnapshotPaths(simulation_path)
    snapshot_paths = isempty(snap_source[:paths]) ? [missing] : snap_source[:paths]

    # Get the path and number of each group catalog file
    groupcat_source = getGroupCatPaths(simulation_path)
    groupcat_paths  = isempty(groupcat_source[:paths]) ? [missing] : groupcat_source[:paths]

    (
        length(snapshot_paths) >= length(groupcat_paths) ||
        throw(ArgumentError("makeSimulationTable: I found less snapshots \
        ($(length(snapshot_paths))) than group catalogs ($(length(groupcat_paths))) in \
        $(simulation_path), I cannot make the table when not every group catalog has a \
        corresponding snapshot"))
    )

    # Number of rows
    n = length(snapshot_paths)

    source_table = DataFrame(
        snapshot_paths = snapshot_paths,
        groupcat_paths = [groupcat_paths; fill(missing, n - length(groupcat_paths))],
    )

    # Add the file name number column
    numbers = snap_source[:numbers]
    insertcols!(source_table, :numbers => isempty(numbers) ? ["000"] : numbers; copycols=false)

    # Get the time stamps of every snapshot
    scale_factors, redshifts, physical_times, lookback_times = computeTimeStamps(snapshot_paths)

    # Add the scale factor column
    insertcols!(source_table, :scale_factors => scale_factors; copycols=false)

    # Add the redshift column
    insertcols!(source_table, :redshifts => redshifts; copycols=false)

    # Add the physical time column
    insertcols!(source_table, :physical_times => physical_times; copycols=false)

    # Add the lookback time column
    insertcols!(source_table, :lookback_times => lookback_times; copycols=false)

    # Sort the table by physical time
    sort!(source_table, :physical_times)

    # Add the row indices
    insertcols!(source_table, :row_id => 1:nrow(source_table))

    return identity.(DataFrame(source_table))

end

"""
    findClosestSnapshot(
        simulation_table::DataFrame,
        times::Vector{<:Unitful.Time},
    )::Vector{Int}

Find the global index, in the context of the simulation, of the snapshot with a physical time closest to each of the ones given in `times`.

!!! note

    This methods uses a precomputed simulation table, made with [`makeSimulationTable`](@ref).

# Arguments

  - `simulation_table::DataFrame`: Dataframe with the path, time stamps, and number of each snapshot and group catalog file in `simulation_path`. It must have the same shape as the one returned by [`makeSimulationTable`](@ref).
  - `times::Vector{<:Unitful.Time}`: Target physical times.

# Returns

  - The indices of the snapshots with physical times closest to `times`.
"""
function findClosestSnapshot(
    simulation_table::DataFrame,
    times::Vector{<:Unitful.Time},
)::Vector{Int}

    # Read the physical time associated to each snapshot
    snap_times = simulation_table[!, :physical_times]

    slices = similar(times, Int)

    # Find the closest snapshot to each of the `times`
    for (i, time) in pairs(times)
        slices[i] = argmin(abs.(snap_times .- time))
    end

    return slices

end

"""
    findClosestSnapshot(simulation_path::String, times::Vector{<:Unitful.Time})::Vector{Int}

Find the global index, in the context of the simulation, of the snapshot with a physical time closest to each of the ones given in `times`.

# Arguments

  - `simulation_path::String`: Path to the simulation directory, set in the code variable `OutputDir`.
  - `times::Vector{<:Unitful.Time}`: Target physical times.

# Returns

  - The indices of the snapshots with physical times closest to `times`.
"""
function findClosestSnapshot(simulation_path::String, times::Vector{<:Unitful.Time})::Vector{Int}

    # Make a dataframe for the simulation with the following columns:
    #  - DataFrame index         -> :row_id
    #  - Number in the file name -> :numbers
    #  - Scale factor            -> :scale_factors
    #  - Redshift                -> :redshifts
    #  - Physical time           -> :physical_times
    #  - Lookback time           -> :lookback_times
    #  - Snapshot path           -> :snapshot_paths
    #  - Group catalog path      -> :groupcat_paths
    simulation_table = makeSimulationTable(simulation_path)

    return findClosestSnapshot(
        simulation_table,
        times,
    )

end

"""
    findClosestSnapshot(simulation_table::DataFrame, time::Unitful.Time)::Int

Find the global index, in the context of the simulation, of the snapshot with a physical time closest to `time`.

!!! note

    This methods uses a precomputed simulation table, made with [`makeSimulationTable`](@ref).

# Arguments

  - `simulation_table::DataFrame`: Dataframe with the path, time stamps, and number of each snapshot and group catalog file in `simulation_path`. It must have the same shape as the one returned by [`makeSimulationTable`](@ref).
  - `time::Unitful.Time`: Target physical time.

# Returns

  - The index of the snapshot with a physical time closest to `time`.
"""
function findClosestSnapshot(simulation_table::DataFrame, time::Unitful.Time)::Int

    return findClosestSnapshot(simulation_table, [time])[1]

end

"""
    findClosestSnapshot(simulation_path::String, time::Unitful.Time)::Int

Find the global index, in the context of the simulation, of the snapshot with a physical time closest to `time`.

# Arguments

  - `simulation_path::String`: Path to the simulation directory, set in the code variable `OutputDir`.
  - `time::Unitful.Time`: Target physical time.

# Returns

  - The index of the snapshot with a physical time closest to `time`.
"""
function findClosestSnapshot(simulation_path::String, time::Unitful.Time)::Int

    return findClosestSnapshot(simulation_path, [time])[1]

end

"""
    makeDataDict(
        simulation_path::String,
        snapshot_n::Int,
        request::Dict{Symbol,Vector{String}},
        simulation_table::DataFrame,
    )::Dict

Construct a data dictionary for a single snapshot.

!!! note

    This methods uses a precomputed simulation table, made with [`makeSimulationTable`](@ref).

# Arguments

  - `simulation_path::String`: Path to the simulation directory, set in the code variable `OutputDir`.
  - `snapshot_n::Int`: Selects the target snapshot. Starts at 1 and is independent of the number in the file name. If every snapshot is present, the relation is `snapshot_n` = (number in filename) + 1.
  - `request::Dict{Symbol,Vector{String}}`: Dictionary with the shape `cell/particle type` -> [`block`, `block`, ...], where the possible types are the keys of [`PARTICLE_INDEX`](@ref), :group, and :subhalo, and the possible blocks are the keys of [`QUANTITIES`](@ref).
  - `simulation_table::DataFrame`: Dataframe with the path, time stamps, and number of each snapshot and group catalog file in `simulation_path`. It must have the same shape as the one returned by [`makeSimulationTable`](@ref).

# Returns

  - A dictionary with the following shape:

      + `:sim_data`          -> ::Simulation (see [`Simulation`](@ref)).
      + `:snap_data`         -> ::Snapshot (see [`Snapshot`](@ref)).
      + `:gc_data`           -> ::GroupCatalog (see [`GroupCatalog`](@ref)).
      + `cell/particle type` -> (`block` -> data, `block` -> data, ...).
      + `cell/particle type` -> (`block` -> data, `block` -> data, ...).
      + `cell/particle type` -> (`block` -> data, `block` -> data, ...).
      + ...
      + `groupcat type`      -> (`block` -> data, `block` -> data, ...).
      + `groupcat type`      -> (`block` -> data, `block` -> data, ...).
      + `groupcat type`      -> (`block` -> data, `block` -> data, ...).
      + ...
"""
function makeDataDict(
    simulation_path::String,
    snapshot_n::Int,
    request::Dict{Symbol,Vector{String}},
    simulation_table::DataFrame,
)::Dict

    snapshot_numbers = simulation_table[!, :numbers]

    (
        length(snapshot_numbers) >= snapshot_n ||
        throw(ArgumentError("makeDataDict: The snapshot number $(snapshot_n) does not exists in  \
        $(simulation_path). There are only $(length(snapshot_numbers)) snapshots. \
        The full simulation table is:\n\n$(simulation_table)"))
    )

    # Select the target snapshot
    snapshot_row = simulation_table[snapshot_n, :]

    ################################################################################################
    # Compute the metadata for the current snapshot and simulation
    ################################################################################################

    # Get the snapshot file path
    snapshot_path = snapshot_row[:snapshot_paths]
    # Get the group catalog file path
    groupcat_path = snapshot_row[:groupcat_paths]

    # Store the metadata of the current snapshot and simulation
    metadata = Dict(
        :sim_data => Simulation(
            simulation_path,
            1,
            snapshot_n,
            isSnapCosmological(snapshot_path),
            simulation_table,
        ),
        :snap_data => Snapshot(
            snapshot_path,
            snapshot_n,
            1,
            snapshot_row[:physical_times],
            snapshot_row[:lookback_times],
            snapshot_row[:scale_factors],
            snapshot_row[:redshifts],
            readSnapHeader(snapshot_path),
        ),
        :gc_data => GroupCatalog(groupcat_path, readGroupCatHeader(groupcat_path)),
    )

    return merge(
        metadata,
        readSnapshot(snapshot_path, request),
        readGroupCatalog(groupcat_path, snapshot_path, request),
    )

end

"""
    makeDataDict(
        simulation_path::String,
        snapshot_n::Int,
        request::Dict{Symbol,Vector{String}},
    )::Dict

Construct a data dictionary for a single snapshot.

# Arguments

  - `simulation_path::String`: Path to the simulation directory, set in the code variable `OutputDir`.
  - `snapshot_n::Int`: Selects the target snapshot. Starts at 1 and is independent of the number in the file name. If every snapshot is present, the relation is `snapshot_n` = (number in filename) + 1.
  - `request::Dict{Symbol,Vector{String}}`: Dictionary with the shape `cell/particle type` -> [`block`, `block`, ...], where the possible types are the keys of [`PARTICLE_INDEX`](@ref), :group, and :subhalo, and the possible blocks are the keys of [`QUANTITIES`](@ref).

# Returns

  - A dictionary with the following shape:

      + `:sim_data`          -> ::Simulation (see [`Simulation`](@ref)).
      + `:snap_data`         -> ::Snapshot (see [`Snapshot`](@ref)).
      + `:gc_data`           -> ::GroupCatalog (see [`GroupCatalog`](@ref)).
      + `cell/particle type` -> (`block` -> data, `block` -> data, ...).
      + `cell/particle type` -> (`block` -> data, `block` -> data, ...).
      + `cell/particle type` -> (`block` -> data, `block` -> data, ...).
      + ...
      + `groupcat type`      -> (`block` -> data, `block` -> data, ...).
      + `groupcat type`      -> (`block` -> data, `block` -> data, ...).
      + `groupcat type`      -> (`block` -> data, `block` -> data, ...).
      + ...
"""
function makeDataDict(
    simulation_path::String,
    snapshot_n::Int,
    request::Dict{Symbol,Vector{String}},
)::Dict

    # Make a dataframe for the simulation with the following columns:
    #  - DataFrame index         -> :row_id
    #  - Number in the file name -> :numbers
    #  - Scale factor            -> :scale_factors
    #  - Redshift                -> :redshifts
    #  - Physical time           -> :physical_times
    #  - Lookback time           -> :lookback_times
    #  - Snapshot path           -> :snapshot_paths
    #  - Group catalog path      -> :groupcat_paths
    simulation_table = makeSimulationTable(simulation_path)

    return makeDataDict(
        simulation_path,
        snapshot_n,
        request,
        simulation_table,
    )

end

"""
    findQtyExtrema(
        simulation_path::String,
        snapshot_n::Int,
        component::Symbol,
        block::String;
        <keyword arguments>
    )::NTuple{2,<:Number}

Compute the minimum and maximum values of `block` in a snapshot or simulation.

# Arguments

  - `simulation_path::String`: Path to the simulation directory, set in the code variable `OutputDir`.
  - `snapshot_n::Int`: Selects which snapshot to plot, starts at 1 and is independent of the number in the file name. If every snapshot is present, the relation is `snapshot_n` = (number in filename) + 1. If set to a negative number, the values in the whole simulation will be compared.
  - `component::Symbol`: Cell/particle type. The possibilities are the keys of [`PARTICLE_INDEX`](@ref).
  - `block::String`: Target block. The possibilities are the keys of [`QUANTITIES`](@ref).
  - `f::Function=identity`: A function with the signature:

    `f(data) -> values`

    where

      + `data::VecOrMat{<:Number}`: Data returned by [`getBlock`](@ref).
      + `values::Vector{<:Number}`: A vector with the values to be compared.

# Returns

  - Tuple with the minimum and maximum values.
"""
function findQtyExtrema(
    simulation_path::String,
    snapshot_n::Int,
    component::Symbol,
    block::String;
    f::Function=identity,
)::NTuple{2,<:Number}

    (
        isdir(simulation_path) ||
        throw(ArgumentError("findQtyExtrema: $(simulation_path) does not exists as a directory"))
    )

    # Make a dataframe for the simulation with the following columns:
    #  - DataFrame index         -> :row_id
    #  - Number in the file name -> :numbers
    #  - Scale factor            -> :scale_factors
    #  - Redshift                -> :redshifts
    #  - Physical time           -> :physical_times
    #  - Lookback time           -> :lookback_times
    #  - Snapshot path           -> :snapshot_paths
    #  - Group catalog path      -> :groupcat_paths
    simulation_table = makeSimulationTable(simulation_path)

    if snapshot_n > 0

        # Get the number in the filename
        file_number = safeSelect(simulation_table[!, :numbers], snapshot_n)

        # Check that after slicing there is one snapshot left
        (
            !isempty(file_number) ||
            throw(ArgumentError("findQtyExtrema: There are no snapshots with `snapshot_n` = \
            $(snapshot_n), the contents of $(simulation_path) are: \n$(simulation_table)"))
        )

        # Find the target row and snapshot path
        snapshot_row = filter(:numbers => ==(file_number), simulation_table)
        snapshot_path = snapshot_row[1, :snapshot_paths]

        (
            !ismissing(snapshot_path) ||
            throw(ArgumentError("findQtyExtrema: The snapshot number $(snapshot_n) seems \
            to be missing"))
        )

        values = f(getBlock(snapshot_path, component, block))

        return extrema(values)

    end

    snapshot_paths = filter!(!ismissing, simulation_table[!, :snapshot_paths])

    (
        !isempty(snapshot_paths) ||
        throw(ArgumentError("findQtyExtrema: I could not find any snapshots in $(simulation_path)"))
    )

    values = [f(getBlock(snapshot_path, component, block)) for snapshot_path in snapshot_paths]

    return extrema(Iterators.flatten(values))

end

"""
    internalUnits(quantity::String, path::String)::Union{Unitful.Quantity,Unitful.Units}

Get the factor to convert a plain number into a [Unitful](https://github.com/PainterQubits/Unitful.jl) quantity, using the correct internal code units.

!!! note

    For non cosmological simulations PHYSICAL_UNITS[] will be set to true, regardless of its original value in ./constants/globals.jl.

# Arguments

  - `quantity::String`: Target quantity. The options are the keys of [`QUANTITIES`](@ref).
  - `path::String`: Path to the snapshot file or folder.

# Returns

  - A [Unitful](https://github.com/PainterQubits/Unitful.jl) quantity or unit.
"""
function internalUnits(quantity::String, path::String)::Union{Unitful.Quantity,Unitful.Units}

    (
        haskey(QUANTITIES, quantity) ||
        throw(ArgumentError("internalUnits: `quantity` should be one of the keys of \
        `QUANTITIES` but I got $(quantity), see the options in `./src/constants/globals.jl`"))
    )

    header = readSnapReducedHeader(path)
    cosmological = isSnapCosmological(path)

    a = cosmological ? header.time : 1.0
    h = cosmological ? header.h : 1.0

    # Set up the struct for unit conversion
    IU = InternalUnits(; l_unit=header.l_unit, m_unit=header.m_unit, v_unit=header.v_unit, a, h)

    dimensions = QUANTITIES[quantity].dimensions
    unit = QUANTITIES[quantity].unit

    if unit == :internal
        if dimensions == Unitful.𝐌

            # From internal units to M⊙
            return IU.m_cosmo

        elseif dimensions == Unitful.𝐋

            if !PHYSICAL_UNITS[] && !cosmological
                @warn(
                    "internalUnits: You have set the unit system to use comoving lengths \
                    (`PHYSICAL_UNITS` = $(PHYSICAL_UNITS[])), but the simulation is not \
                    cosmological. internalUnits will default to physical lengths \
                    (`PHYSICAL_UNITS` = true). Check `PHYSICAL_UNITS` in `constants/globals.jl`",
                    maxlog = 1,
                )
                PHYSICAL_UNITS[] = true
            end

            # From internal units to kpc
            if !PHYSICAL_UNITS[]
                return IU.x_comoving
            else
                return IU.x_cosmo
            end

        elseif dimensions == Unitful.𝐓

            # From internal units to Myr for non-cosmological simulations,
            # and to a dimensionless quantity for cosmological simulations
            return cosmological ? Unitful.NoUnits : IU.t_newton

        elseif dimensions == Unitful.𝐌 * Unitful.𝐋^-3

            # From internal units to g * cm^-3
            return IU.rho_cgs

        elseif dimensions == Unitful.𝐋^2 * Unitful.𝐓^-2

            # From internal units to erg * g^-1
            return IU.U_cgs

        elseif dimensions == Unitful.𝐋 * Unitful.𝐓^-1

            # From internal units to km * s^-1
            return IU.v_cosmo

        elseif dimensions == Unitful.𝐌 * Unitful.𝐋^-1 * Unitful.𝐓^-2

            # From internal units to Pa
            return IU.P_Pa

        else

            error("internalUnits: I don't know the internal units of a quantity \
            with dimensions $(dimensions)")

        end

    elseif unit == :gvel

        # Special case for "G_Vel" (velocity of the group)
        # See the TNG documentation https://www.tng-project.org/data/docs/specifications/
        return IU.v_cosmo / a^1.5

    elseif unit == :pot

        # Special case for "Potential" (gravitational potential)
        # See the TNG documentation https://www.tng-project.org/data/docs/specifications/
        return (1.0 / a) * u"km^2 * s^-2"

    else

        return unit

    end

end
