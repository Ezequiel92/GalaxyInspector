####################################################################################################
# Data acquisition functions
####################################################################################################

"""
    readGroupCatHeader(path::Union{String,Missing}; <keyword arguments>)::GroupCatHeader

Read the header of a group catalog in the HDF5 format.

!!! note

    If each group catalog is made of multiple files, I'll read the header on the first one.

# Arguments

  - `path::Union{String,Missing}`: Path to the group catalog file or folder.
  - `warnings::Bool=true`: If a warning will be given when `path` is missing.

# Returns

  - A [`GroupCatHeader`](@ref).
"""
function readGroupCatHeader(path::Union{String,Missing}; warnings::Bool=true)::GroupCatHeader

    if ismissing(path)

        !warnings || @warn("readGroupCatHeader: The group catalog file or folder is missing")

        return GroupCatHeader()

    elseif isfile(path)

        (
            HDF5.ishdf5(path) ||
            throw(ArgumentError("readGroupCatHeader: The file $(path) is not in the \
            HDF5 format, I don't know how to read it"))
        )

        file_path = path

    elseif isdir(path)

        sub_files = glob("$(GC_BASENAME)_*.*.hdf5", path)

        (
            !isempty(sub_files) && all(HDF5.ishdf5, sub_files) ||
            throw(ArgumentError("readGroupCatHeader: The directory $(path) does not contain \
            group catalog sub-files in the HDF5 format"))
        )

        file_path = minimum(sub_files)

    else

        throw(ArgumentError("readGroupCatHeader: $(path) does not exist as a file or folder"))

    end

    header = h5open(file_path, "r") do gc_file

        h = gc_file["Header"]

        # Read which attributes are present
        attrs_present = keys(HDF5.attrs(h))

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
            box_size          = read_attribute(h, "BoxSize"),
            h0                = read_attribute(h, "HubbleParam"),
            n_groups_part     = read_attribute(h, "Ngroups_ThisFile"),
            n_groups_total    = read_attribute(h, "Ngroups_Total"),
            n_subgroups_part  = read_attribute(h, "Nsubgroups_ThisFile"),
            n_subgroups_total = read_attribute(h, "Nsubgroups_Total"),
            num_files         = read_attribute(h, "NumFiles"),
            omega_0           = read_attribute(h, "Omega0"),
            omega_l           = read_attribute(h, "OmegaLambda"),
            redshift          = read_attribute(h, "Redshift"),
            time              = read_attribute(h, "Time"),
        )

    end

    return header

end

"""
    readSnapHeader(path::String)::SnapshotHeader

Read the header of a snapshot in the HDF5 format.

!!! note

    If each snapshot is made of multiple files, I'll read the header on the first chunck.

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

        # Count the number of stellar particles, excluding wind particles
        num_part_stars = countStars(file_path)
        num_total_stars = num_part_stars

    elseif isdir(path)

        sub_files = glob("$(SNAP_BASENAME)_*.*.hdf5", path)

        (
            !isempty(sub_files) && all(HDF5.ishdf5, sub_files) ||
            throw(ArgumentError("readSnapHeader: The directory $(path) does not contain \
            snapshot sub-files in the HDF5 format"))
        )

        file_path = minimum(sub_files)

        # Count the number of stellar particles, excluding wind particles
        num_part_stars = countStars(file_path)
        num_total_stars = sum(countStars, sub_files)

    else

        throw(ArgumentError("readSnapHeader: $(path) does not exist as a file or folder"))

    end

    header = h5open(file_path, "r") do snap_file

        h = snap_file["Header"]

        # Read which attributes are present
        attrs_present = keys(HDF5.attrs(h))

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
        num_part = read_attribute(h, "NumPart_ThisFile")
        num_part[PARTICLE_INDEX[:stars] + 1] = num_part_stars
        num_total = read_attribute(h, "NumPart_Total")
        num_total[PARTICLE_INDEX[:stars] + 1] = num_total_stars

        # Check if the length units are in the header, otherwise use the IllustrisTNG values
        if "UnitLength_in_cm" ∈ attrs_present
            l_unit = read_attribute(h, "UnitLength_in_cm") * u"cm"
        else
            l_unit = ILLUSTRIS_L_UNIT
        end

        # Check if the mass units are in the header, otherwise use the IllustrisTNG values
        if "UnitMass_in_g" ∈ attrs_present
            m_unit = read_attribute(h, "UnitMass_in_g") * u"g"
        else
            m_unit = ILLUSTRIS_M_UNIT
        end

        # Check if the velocity units are in the header, otherwise use the IllustrisTNG values
        if "UnitVelocity_in_cm_per_s" ∈ attrs_present
            v_unit = read_attribute(h, "UnitVelocity_in_cm_per_s") * u"cm*s^-1"
        else
            v_unit = ILLUSTRIS_V_UNIT
        end

        SnapshotHeader(
            box_size   = read_attribute(h, "BoxSize"),
            h0         = read_attribute(h, "HubbleParam"),
            mass_table = read_attribute(h, "MassTable"),
            num_files  = read_attribute(h, "NumFilesPerSnapshot"),
            num_part   = num_part,
            num_total  = num_total,
            omega_0    = read_attribute(h, "Omega0"),
            omega_l    = read_attribute(h, "OmegaLambda"),
            redshift   = read_attribute(h, "Redshift"),
            time       = read_attribute(h, "Time"),
            l_unit     = l_unit,
            m_unit     = m_unit,
            v_unit     = v_unit,
        )

    end

    return header

end

"""
    isBlockPresent(block::String, group::HDF5.Group)::Bool

Checks if a given data block exist in a HDF5 group.

# Arguments

  - `block::String`: Target block. The possibilities are the keys of [`QUANTITIES`](@ref).
  - `group::HDF5.Group`: HDF5 group.

# Returns

  - If `block` exist in `group`.
"""
function isBlockPresent(block::String, group::HDF5.Group)::Bool

    (
        block ∈ keys(QUANTITIES) ||
        throw(ArgumentError("isBlockPresent: `block` should be a key of `QUANTITIES`, \
        but I got $(block), see the options in `./src/constants/globals.jl`"))
    )

    return QUANTITIES[block].hdf5_name ∈ keys(group)

end

"""
    isBlockPresent(component::Symbol, block::String, path::String)::Bool

Checks if a given block exist in a snapshot.

!!! note

    If each snapshot is made of multiple files, I'll check only in the first chunck.

# Arguments

  - `component::Symbol`: The cell/particle type of the target block. The possibilities are the keys of [`PARTICLE_INDEX`](@ref).
  - `block::String`: Target block. The possibilities are the keys of [`QUANTITIES`](@ref).
  - `path::String`: Path to the snapshot file or folder.

# Returns

  - If `block` exist in the snapshot.
"""
function isBlockPresent(component::Symbol, block::String, path::String)::Bool

    if isfile(path)

        (
            HDF5.ishdf5(path) ||
            throw(ArgumentError("isBlockPresent: The file $(path) is not in the HDF5 format, \
            I don't know how to read it"))
        )

        file_path = path

    elseif isdir(path)

        sub_files = glob("$(SNAP_BASENAME)_*.*.hdf5", path)

        (
            !isempty(sub_files) && all(HDF5.ishdf5, sub_files) ||
            throw(ArgumentError("isBlockPresent: The directory $(path) does not contain \
            snapshot sub-files in the HDF5 format"))
        )

        file_path = minimum(sub_files)

    else

        throw(ArgumentError("isBlockPresent: $(path) does not exist as a file or folder"))

    end

    response = h5open(file_path, "r") do snapshot

        type_str = PARTICLE_CODE_NAME[component]

        if type_str ∈ keys(snapshot)
            isBlockPresent(block, snapshot[type_str])
        else
            false
        end

    end

    return response

end

"""
    readTime(path::String)::Float64

Read the "Time" field in the header of a snapshot file.

!!! note

    If each snapshot is made of multiple files, I'll read the header on the first chunck.

# Arguments

  - `path::String`: Path to the snapshot file or folder.

# Returns

  - The "Time" field in the header (for cosmological simulations is the scale factor).
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

        sub_files = glob("$(SNAP_BASENAME)_*.*.hdf5", path)

        (
            !isempty(sub_files) && all(HDF5.ishdf5, sub_files) ||
            throw(ArgumentError("readTime: The directory $(path) does not contain \
            snapshot sub-files in the HDF5 format"))
        )

        file_path = minimum(sub_files)

    else

        throw(ArgumentError("readTime: $(path) does not exist as a file or folder"))

    end

    time = h5open(file_path, "r") do file
        read_attribute(file["Header"], "Time")
    end

    return time

end

"""
    readTemperature(file_path::String)::Vector{<:Unitful.Temperature}

Compute the temperature of the gas cells in a snapshot.

# Arguments

  - `file_path::String`: Path to the snapshot file.

# Returns

  - The temperature of the gas cells.
"""
function readTemperature(file_path::String)::Vector{<:Unitful.Temperature}

    if isfile(file_path)

        (
            HDF5.ishdf5(file_path) ||
            throw(ArgumentError("readTemperature: The file $(file_path) is not in the \
            HDF5 format, I don't know how to read it"))
        )

    else

        throw(ArgumentError("readTemperature: $(file_path) does not exist as a file"))

    end

    # List of blocks needed to compute the temperature
    blocks = ["U   ", "NE  "]

    temp_data = h5open(file_path, "r") do snapshot

        group = snapshot[PARTICLE_CODE_NAME[:gas]]

        # Get the indices of the missing blocks
        idx_missing = map(x -> !isBlockPresent(x, group), blocks)
        (
            !any(idx_missing) ||
            throw(ArgumentError("readTemperature: The blocks $(blocks[idx_missing]) \
            are missing, and I need them to compute the temperature"))
        )

        # Compute the unit factor for each block
        units = internalUnits.(blocks, file_path)

        [read(group, QUANTITIES[block].hdf5_name) .* unit for (unit, block) in zip(units, blocks)]

    end

    return computeTemperature(temp_data...)

end

"""
    readGoupCatBlocks(
        file_path::String,
        snapshot_path::String,
        request::Dict{Symbol,Vector{String}};
        <keyword arguments>
    )::Dict{Symbol,Dict{String,VecOrMat{<:Number}}}

Read the specified blocks from a group catalog file.

# Arguments

  - `file_path::String`: Path to the group catalog file.
  - `snapshot_path::String`: Path to the corresponding snapshot file or folder. This is needed for unit conversion.
  - `request::Dict{Symbol,Vector{String}}`: The blocks to be read. It must have the shape `group type` -> [`block`, `block`, `block`].
  - `warnings::Bool=true`: If a warning will be given when there are missing blocks.

# Returns

  - A dictionary with the following shape: `group type` -> (`block` -> data of `block`).
"""
function readGoupCatBlocks(
    file_path::String,
    snapshot_path::String,
    request::Dict{Symbol,Vector{String}};
    warnings::Bool=true,
)::Dict{Symbol,Dict{String,VecOrMat{<:Number}}}

    if isfile(file_path)

        (
            HDF5.ishdf5(file_path) ||
            throw(ArgumentError("readGoupCatBlocks: The file $(file_path) is not in the \
            HDF5 format, I don't know how to read it"))
        )

    else

        throw(ArgumentError("readGoupCatBlocks: $(file_path) does not exist as a file"))

    end

    # Allocate memory
    output = Dict{Symbol,Dict{String,VecOrMat{<:Number}}}()

    h5open(file_path, "r") do gc_file

        # Read from the request only the group catalog types
        @inbounds for component in groupcatTypes(request)

            blocks = copy(request[component])

            type_str = titlecase(string(component))

            # Allocate memory
            qty_data = Dict{String,VecOrMat{<:Number}}()

            if type_str ∈ keys(gc_file)

                # Read the HDF5 group
                hdf5_group = gc_file[type_str]

                if isempty(hdf5_group)

                    (
                        !warnings ||
                        @warn("readGoupCatBlocks: The group catalog type :$(component) \
                        in $(file_path) is empty")
                    )

                    # Return an empty array for every missing block
                    @inbounds for block in blocks
                        unit = internalUnits(block, snapshot_path)
                        qty_data[block] = typeof(1.0 * unit)[]
                    end

                else

                    @inbounds for block in blocks

                        unit = internalUnits(block, snapshot_path)

                        if isBlockPresent(block, hdf5_group)

                            qty_data[block] = read(hdf5_group, QUANTITIES[block].hdf5_name) .* unit

                        else

                            (
                                !warnings ||
                                @warn("readGoupCatBlocks: The block $(block) for the group \
                                catalog type :$(component) in $(file_path) is missing")
                            )
                            # Return an empty array for every missing block
                            qty_data[block] = typeof(1.0 * unit)[]

                        end

                    end

                end

            else

                (
                    !warnings ||
                    @warn("readGoupCatBlocks: The group catalog type \
                    :$(component) in $(file_path) is missing")
                )

                # Return an empty array for every missing block
                @inbounds for block in blocks
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
    readSnapBlocks(
        file_path::String,
        request::Dict{Symbol,Vector{String}};
        <keyword arguments>
    )::Dict{Symbol,Dict{String,VecOrMat{<:Number}}}

Read the specified blocks from a snapshot file.

# Arguments

  - `file_path::String`: Path to the snapshot file.
  - `request::Dict{Symbol,Vector{String}}`: The blocks to be read. It must have the shape `cell/particle type` -> [`block`, `block`, `block`].
  - `warnings::Bool=true`: If a warning will be given when there are missing blocks.

# Returns

  - A dictionary with the following shape: `cell/particle type` -> (`block` -> data of `block`).
"""
function readSnapBlocks(
    file_path::String,
    request::Dict{Symbol,Vector{String}};
    warnings::Bool=true,
)::Dict{Symbol,Dict{String,VecOrMat{<:Number}}}

    if isfile(file_path)

        (
            HDF5.ishdf5(file_path) ||
            throw(ArgumentError("readSnapBlocks: The file $(file_path) is not in the \
            HDF5 format, I don't know how to read it"))
        )

    else

        throw(ArgumentError("readSnapBlocks: $(file_path) does not exist as a file"))

    end

    # Allocate memory
    output = Dict{Symbol,Dict{String,VecOrMat{<:Number}}}()

    # Read the header
    header = readSnapHeader(file_path)

    h5open(file_path, "r") do snapshot

        # Read from the request only the cell/particle types
        @inbounds for component in snapshotTypes(request)

            blocks = copy(request[component])

            type_str = PARTICLE_CODE_NAME[component]

            # Allocate memory
            qty_data = Dict{String,VecOrMat{<:Number}}()

            if type_str ∈ keys(snapshot)

                # Read the HDF5 group
                hdf5_group = snapshot[type_str]

                if isempty(hdf5_group)

                    (
                        !warnings ||
                        @warn("readSnapBlocks: The cell/particle type \
                        :$(component) in $(file_path) is empty")
                    )

                    # Return an empty array for every missing block
                    @inbounds for block in blocks
                        unit = internalUnits(block, snapshot_path)
                        qty_data[block] = typeof(1.0 * unit)[]
                    end

                else

                    # For the stellar particles, exclude wind particles
                    if component == :stars
                        idxs = findRealStars(file_path)
                    else
                        idxs = (:)
                    end

                    @inbounds for block in blocks

                        unit = internalUnits(block, file_path)

                        if block == "TEMP"

                            (
                                component == :gas ||
                                throw(ArgumentError("readSnapBlocks: I can't compute the \
                                temperature for cells/particles of type :$(component), \
                                only for :gas"))
                            )

                            qty_data["TEMP"] = readTemperature(file_path)

                        elseif block == "MASS"

                            # Read the mass table from the header
                            mass_table = header.mass_table[PARTICLE_INDEX[component] + 1]

                            if iszero(mass_table)

                                raw = read(hdf5_group, QUANTITIES["MASS"].hdf5_name)
                                qty_data["MASS"] = selectdim(raw, ndims(raw), idxs) .* unit

                            else

                                # All cell/particles have the same mass
                                cp_number = header.num_part[PARTICLE_INDEX[component] + 1]
                                qty_data["MASS"] = fill(mass_table, cp_number) .* unit

                            end

                        elseif isBlockPresent(block, hdf5_group)

                            raw = read(hdf5_group, QUANTITIES[block].hdf5_name)
                            qty_data[block] = selectdim(raw, ndims(raw), idxs) .* unit

                        else

                            (
                                !warnings ||
                                @warn("readSnapBlocks: The block $(block) for the \
                                cell/particle type :$(component) in $(file_path) is missing")
                            )
                            # Return an empty array for every missing block
                            qty_data[block] = typeof(1.0 * unit)[]

                        end

                    end

                end

            else

                (
                    !warnings ||
                    @warn("readSnapBlocks: The cell/particle type \
                    :$(component) in $(file_path) is missing")
                )

                @inbounds for block in blocks
                    unit = internalUnits(block, file_path)
                    qty_data[block] = typeof(1.0 * unit)[]
                end

            end

            output[component] = qty_data

        end

    end

    return output

end

"""
    readGroupCatalog(
        path::Union{String,Missing},
        snapshot_path::String,
        request::Dict{Symbol,Vector{String}};
        <keyword arguments>
    )::Dict{Symbol,Dict{String,VecOrMat{<:Number}}}

Read the specified blocks from a group catalog file or folder.

# Arguments

  - `path::Union{String,Missing}`: Path to the group catalog file or folder.
  - `snapshot_path::String`: Path to the corresponding snapshot file or folder. This is needed for unit conversion.
  - `request::Dict{Symbol,Vector{String}}`: Which blocks will be read. It must have the shape `group type` -> [`block`, `block`, `block`].
  - `warnings::Bool=true`: If a warning will be given when there are missing blocks.

# Returns

  - A dictionary with the following shape: `group type` -> (`block` -> data of `block`).
"""
function readGroupCatalog(
    path::Union{String,Missing},
    snapshot_path::String,
    request::Dict{Symbol,Vector{String}};
    warnings::Bool=true,
)::Dict{Symbol,Dict{String,VecOrMat{<:Number}}}

    if ismissing(path)

        !warnings || @warn("readGroupCatalog: The group catalog file or folder is missing")

        return Dict{Symbol,Dict{String,VecOrMat{<:Number}}}()

    elseif isfile(path)

        (
            HDF5.ishdf5(path) ||
            throw(ArgumentError("readGroupCatalog: The file $(path) is not in the HDF5 format, \
            I don't know how to read it"))
        )

        return readGoupCatBlocks(path, snapshot_path, request; warnings)

    elseif isdir(path)

        sub_files = glob("$(GC_BASENAME)_*.*.hdf5", path)

        (
            !isempty(sub_files) && all(HDF5.ishdf5, sub_files) ||
            throw(ArgumentError("readGroupCatalog: The directory $(path) does not contain \
            group catalog sub-files in the HDF5 format"))
        )

        # Sort the sub files to concatenate the data in them correctly
        sort!(sub_files)

        # Read the data in each sub file
        data_in_files = [
            readGoupCatBlocks(file, snapshot_path, request; warnings) for file in sub_files
        ]

        # Allocate memory
        output = Dict{Symbol,Dict{String,VecOrMat{<:Number}}}()

        @inbounds for (component, data_blocks) in first(data_in_files)

            qty_data = Dict{String,VecOrMat{<:Number}}()

            @inbounds for block in keys(data_blocks)

                raw = [
                    data_in_file[component][block] for
                    data_in_file in data_in_files if !isempty(data_in_file[component][block])
                ]

                qty_data[block] = cat(raw...; dims=ndims(first(raw)))

            end

            output[component] = qty_data

        end

        return output

    else

        throw(ArgumentError("readGroupCatalog: $(path) does not exist as a file or folder"))

    end

end

"""
    readSnapshot(
        path::Union{String,Missing},
        request::Dict{Symbol,Vector{String}};
        <keyword arguments>
    )::Dict{Symbol,Dict{String,VecOrMat{<:Number}}}

Read the specified blocks from a snapshot file or folder.

# Arguments

  - `path::Union{String,Missing}`: Path to the snapshot file or folder.
  - `request::Dict{Symbol,Vector{String}}`: Which blocks will be read. It must have the shape `cell/particle type` -> [`block`, `block`, `block`].
  - `warnings::Bool=true`: If a warning will be given when some/all the data is missing.

# Returns

  - A dictionary with the following shape: `cell/particle type` -> (`block` -> data of `block`).
"""
function readSnapshot(
    path::Union{String,Missing},
    request::Dict{Symbol,Vector{String}};
    warnings::Bool=true,
)::Dict{Symbol,Dict{String,VecOrMat{<:Number}}}

    if ismissing(path)

        throw(ArgumentError("readSnapshot: The snapshot file or folder is missing"))

    elseif isfile(path)

        (
            HDF5.ishdf5(path) ||
            throw(ArgumentError("readSnapshot: The file $(path) is not in the HDF5 format, \
            I don't know how to read it"))
        )

        return readSnapBlocks(path, request; warnings)

    elseif isdir(path)

        sub_files = glob("$(SNAP_BASENAME)_*.*.hdf5", path)

        (
            !isempty(sub_files) && all(HDF5.ishdf5, sub_files) ||
            throw(ArgumentError("readSnapshot: The directory $(path) does not contain \
            snapshot sub-files in the HDF5 format"))
        )

        # Sort the sub files to concatenate the data in them correctly
        sort!(sub_files)

        # Read the data in each sub file
        data_in_files = [readSnapBlocks(file, request; warnings) for file in sub_files]

        # Allocate memory
        output = Dict{Symbol,Dict{String,VecOrMat{<:Number}}}()

        @inbounds for (component, data_blocks) in first(data_in_files)

            qty_data = Dict{String,VecOrMat{<:Number}}()

            @inbounds for block in keys(data_blocks)

                raw = [
                    data_in_file[component][block] for
                    data_in_file in data_in_files if !isempty(data_in_file[component][block])
                ]

                if isempty(raw)
                    qty_data[block] = Number[]
                else
                    qty_data[block] = cat(raw...; dims=ndims(first(raw)))
                end

            end

            output[component] = qty_data

        end

        return output

    else

        throw(ArgumentError("readSnapshot: $(path) does not exist as a file or folder"))

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
    readSfrFile(
        file_path::String,
        snap_path::String;
        <keyword arguments>
    )::Dict{Int32,VecOrMat{<:Number}}

Read the `sfr.txt` file.

# Arguments

  - `file_path::String`: Path to the `sfr.txt` file.
  - `snapshot_path::String`: Path to one snapshot file or folder of the simulation. This is needed for unit conversion.
  - `warnings::Bool=true`: If a warning will be given when the `sfr.txt` file does not have the expected structure.

# Returns

  - A dictionary with the following shape:

      + `1` -> Time or scale factor (internal units).
      + `2` -> Total stellar mass to be formed prior to stochastic sampling (internal units).
      + `3` -> Instantaneous star formation rate of all cells (``\\mathrm{M_\\odot \\, yr^{-1}}``).
      + `4` -> Instantaneous star formation rate of active cells (``\\mathrm{M_\\odot \\, yr^{-1}}``).
      + `5` -> Total mass in stars formed after stochastic sampling (internal units).
      + `6` -> Cumulative stellar mass formed (internal units).
"""
function readSfrFile(
    file_path::String,
    snap_path::String;
    warnings::Bool=true,
)::Dict{Int32,VecOrMat{<:Number}}

    isfile(file_path) || throw(ArgumentError("readSfrFile: $(file_path) does not exist as a file"))

    # Load the data from the `sfr.txt` file
    file_data = readdlm(file_path, Float64)
    n_cols = size(file_data, 2)

    # Check that the data in the file has the correct size
    (
        n_cols <= 6 ||
        throw(ArgumentError("readSfrFile: I don't know how to handle more \
        than 6 columns in `sfr.txt`"))
    )
    (
        !(warnings && n_cols < 6) ||
        @warn("readSfrFile: I could only find $(n_cols) columns \
        in $(flie_path). I was expecting 6")
    )

    # Load the units for each column
    units = [internalUnits("SFC$(i)", snap_path) for i in 1:n_cols]

    return Dict(i => column .* units[i] for (i, column) in pairs(eachcol(file_data)))

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
  - `targets::Vector{String}`: Target processes.
  - `step::Int=1`: Step used to traverse the rows.
  - `warnings::Bool=true`: If a warning will be given when there are missing targets.

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
    warnings::Bool=true,
)::Dict{String,Matrix{Float64}}

    isfile(file_path) || throw(ArgumentError("readCpuFile: $(file_path) does not exist as a file"))

    # Load the data from the `cpu.txt` file
    file_data = eachline(file_path)

    # Set up an auxiliary dictionary
    data_aux = Dict(target => Matrix{Float64}[] for target in targets)

    # Clock time for each sync-point
    time = 0.0
    # Time step
    time_step = 0

    @inbounds for line in file_data

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
        !(warnings && any(isempty, values(data_aux))) ||
        @warn("readCpuFile: I could not find some of the target rows in $(file_path)")
    )

    # Allocate memory
    data_out = Dict{String,Matrix{Float64}}()

    # Try reducing the data size
    @inbounds for (target, values) in data_aux

        # Ignore empty targets
        !isempty(values) || continue

        l_e = length(values)

        if 1 < step < l_e
            data_out[target] = vcat(values[1:step:end]...)
        else
            data_out[target] = vcat(values...)
            (
                !(warnings && step > l_e) ||
                @warn("readCpuFile: `step` = $(step) is bigger than the number \
                of time steps in $(file_path)")
            )
        end

    end

    return data_out

end

"""
    getSnapshotPaths(simulation_path::String; <keyword arguments>)::Dict{Symbol,Vector{String}}

Find the path and number of every snapshot in `simulation_path`.

!!! note

    If each snapshot is made of multiple files, the `:paths` field will have paths to folders, each one containing the sub-files of the corresponding snapshot.

# Arguments

  - `simulation_path::String`: Path to the simulation directory, set in the code variable `OutputDir`.
  - `warnings::Bool=true`: If a warning will be raised when no snapshot files or folders are found.

# Returns

  - A dictionary with the following shape:

      + `:numbers` -> The number that characterize each snapshot.
      + `:paths`   -> The full path to each snapshot.
"""
function getSnapshotPaths(simulation_path::String; warnings::Bool=true)::Dict{Symbol,Vector{String}}

    (
        isdir(simulation_path) ||
        throw(ArgumentError("getSnapshotPaths: $(simulation_path) does not exist as a directory"))
    )

    # Get the full list of paths to every snapshot in `simulation_path`
    path_list = [
        glob("*/*/$(SNAP_BASENAME)_*", simulation_path)
        glob("*/$(SNAP_BASENAME)_*", simulation_path)
        glob("$(SNAP_BASENAME)_*", simulation_path)
    ]

    # Check for an empty folder
    if isempty(path_list)

        (
            !warnings ||
            @warn("getSnapshotPaths: I could not find any file named \
            $(SNAP_BASENAME)_*.hdf5 within $(simulation_path), or any of its subfolders")
        )

        return Dict(:numbers => String[], :paths => String[])

    end

    # Get the numbers that characterize each snapshot
    reg = Regex("(?<=$(SNAP_BASENAME)_).*?(?=(?:\\.)|\$)")
    number_list = map(x -> match(reg, x).match, path_list)

    if readSnapHeader(first(path_list)).num_files > 1
        # If there are multiple files per snapshot, get the path to the snapshot directory
        map!(dirname, path_list, path_list)
        # Delete duplicates
        unique!(path_list)
        unique!(number_list)
    end

    return Dict(:numbers => number_list, :paths => sort(normpath.(path_list)))

end

"""
    getGroupCatPaths(simulation_path::String; <keyword arguments>)::Dict{Symbol,Vector{String}}

Find the path and number of every group catalog in `simulation_path`.

!!! note

    If each group catalog is made of multiple files, the `:paths` field will have paths to folders, each one containing the sub-files of the corresponding group catalog.

# Arguments

  - `simulation_path::String`: Path to the simulation directory, set in the code variable `OutputDir`.
  - `warnings::Bool=true`: If a warning will be raised when no group catalog file or folders are found.

# Returns

  - A dictionary with the following shape:

      + `:numbers` -> The number that characterize each group catalog.
      + `:paths`   -> The full path to each group catalog.
"""
function getGroupCatPaths(simulation_path::String; warnings::Bool=true)::Dict{Symbol,Vector{String}}

    (
        isdir(simulation_path) ||
        throw(ArgumentError("getGroupCatPaths: $(simulation_path) does not exist as a directory"))
    )

    # Get the full list of paths to every group catalog in `simulation_path`
    path_list = [
        glob("*/*/$(GC_BASENAME)_*", simulation_path)
        glob("*/$(GC_BASENAME)_*", simulation_path)
        glob("$(GC_BASENAME)_*", simulation_path)
    ]

    # Check for an empty folder
    if isempty(path_list)

        (
            !warnings ||
            @warn("getGroupCatPaths: I could not find any file named \
            $(GC_BASENAME)_*.hdf5 within $(simulation_path), or any of its subfolders")
        )

        return Dict(:numbers => String[], :paths => String[])

    end

    # Get the numbers that characterize each group catalog
    reg = Regex("(?<=$(GC_BASENAME)_).*?(?=(?:\\.)|\$)")
    number_list = map(x -> match(reg, x).match, path_list)

    if readGroupCatHeader(first(path_list); warnings).num_files > 1
        # If there are multiple files per group catalog, get the path to the group catalog directory
        map!(dirname, path_list, path_list)
        # Delete duplicates
        unique!(path_list)
        unique!(number_list)
    end

    return Dict(:numbers => number_list, :paths => sort(normpath.(path_list)))

end

"""
    makeSimulationTable(simulation_path::String; <keyword arguments>)::DataFrame

Construct a dataframe with the path, time stamps and number of each snapshot and group catalog file in `simulation_path`.

# Arguments

  - `simulation_path::String`: Path to the simulation directory, set in the code variable `OutputDir`.
  - `warnings::Bool=true`: If a warning will be raised when there are missing files.

# Returns

  - A dataframe with 8 colums:

      + `:ids`            -> Dataframe index of each snapshot, i.e. if there are 10 snapshots in total it runs from 1 to 10.
      + `:numbers`        -> Number in the file name of each snapshot.
      + `:scale_factors`  -> Scale factor of each snapshot.
      + `:redshifts`      -> Redshift of each snapshot.
      + `:physical_times` -> Physical time since the Big Bang of each snapshot.
      + `:lookback_times` -> Physical time left to reach the last snapshot.
      + `:snapshot_paths` -> Full path to the snapshots.
      + `:groupcat_paths` -> Full path to the group catalog files.
"""
function makeSimulationTable(simulation_path::String; warnings::Bool=true)::DataFrame

    # Get the path and number of each snapshot
    snap_source = getSnapshotPaths(simulation_path; warnings)
    snapshot_paths = isempty(snap_source[:paths]) ? [missing] : snap_source[:paths]

    # Get the path and number of each group catalog file
    groupcat_source = getGroupCatPaths(simulation_path; warnings)
    groupcat_paths  = isempty(groupcat_source[:paths]) ? [missing] : groupcat_source[:paths]

    (
        length(snapshot_paths) >= length(groupcat_paths) ||
        throw(ArgumentError("makeSimulationTable: I found less snapshots than group catalogs \
        in $(simulation_path), I cannot make the table when not every group catalog has a \
        corresponding snapshot"))
    )

    paths  = [snapshot_paths, groupcat_paths]
    labels = [:snapshot_paths, :groupcat_paths]
    rows   = [[1:length(snapshot_paths);], [1:length(groupcat_paths);]]

    source_table = unstack(flatten(DataFrame(l=labels, p=paths, ids=rows), [:p, :ids]), :l, :p)

    # Add the file name number column
    numbers = snap_source[:numbers]
    insertcols!(source_table, 2, :numbers => isempty(numbers) ? ["000"] : numbers; copycols=false)

    # Get the time stamps of every snapshot
    scale_factors, redshifts, physical_times, lookback_times = computeTimeTicks(snapshot_paths)

    # Add the scale factor column
    insertcols!(source_table, 3, :scale_factors => scale_factors; copycols=false)

    # Add the redshift column
    insertcols!(source_table, 4, :redshifts => redshifts; copycols=false)

    # Add the physical time column
    insertcols!(source_table, 5, :physical_times => physical_times; copycols=false)

    # Add the lookback time column
    insertcols!(source_table, 6, :lookback_times => lookback_times; copycols=false)

    return identity.(DataFrame(source_table))

end

"""
    makeDataDict(
        simulation_path::String,
        slice_n::Int,
        request::Dict{Symbol,Vector{String}};
        <keyword arguments>
    )::Dict

Construct a data dictionary for a single snapshot.

# Arguments

  - `simulation_path::String`: Path to the simulation directory, set in the code variable `OutputDir`.
  - `slice_n::Int`: Selects the target snapshot. Starts at 1 and is independent of the number in the file name. If every snapshot is present, `slice_n` = filename_number + 1.
  - `request::Dict{Symbol,Vector{String}}`: Dictionary with the shape `cell/particle type` -> [`block`, `block`, ...], where the possible types are the keys of [`PARTICLE_INDEX`](@ref), and the possible quantities are the keys of [`QUANTITIES`](@ref).
  - `warnings::Bool=true`: If a warning will be raised when there is messing data.

# Returns

  - A dictionary with the following shape:

      + `:sim_data`          -> ::Simulation (see [`Simulation`](@ref)).
      + `:snap_data`         -> ::Snapshot (see [`Snapshot`](@ref)).
      + `:gc_data`           -> ::GroupCatalog (see [`GroupCatalog`](@ref)).
      + `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
      + `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
      + `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
      + ...
      + `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
      + `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
      + `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
      + ...
"""
function makeDataDict(
    simulation_path::String,
    slice_n::Int,
    request::Dict{Symbol,Vector{String}};
    warnings::Bool=true,
)::Dict

    # Make a dataframe for every simulation, with the following columns:
    #   - 1. DataFrame index
    #   - 2. Number in the file name
    #   - 3. Scale factor
    #   - 4. Redshift
    #   - 5. Physical time
    #   - 6. Lookback time
    #   - 7. Snapshot path
    #   - 8. Group catalog path
    simulation_table = makeSimulationTable(simulation_path; warnings)

    snapshot_numbers = simulation_table[!, :numbers]

    (
        length(snapshot_numbers) >= slice_n ||
        throw(ArgumentError("makeDataDict: The snapshot number $(slice_n) does not exist in  \
        $(simulation_path). There are only $(length(snapshot_numbers)) snapshots. \
        The full simulation table is:\n\n$(simulation_table)"))
    )

    # Select the target snapshot
    snapshot_row = simulation_table[slice_n, :]

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
            slice_n,
            isCosmological(snapshot_path),
            simulation_table,
        ),
        :snap_data => Snapshot(
            snapshot_path,
            slice_n,
            slice_n,
            snapshot_row[:physical_times],
            snapshot_row[:lookback_times],
            snapshot_row[:scale_factors],
            snapshot_row[:redshifts],
            readSnapHeader(snapshot_path),
        ),

        :gc_data => GroupCatalog(
            groupcat_path,
            readGroupCatHeader(groupcat_path; warnings),
        ),
    )

    return merge(
        metadata,
        readSnapshot(snapshot_path, request; warnings),
        readGroupCatalog(groupcat_path, snapshot_path, request; warnings),
    )

end

"""
    countSnapshot(simulation_path::String; warnings::Bool=true)::Int

Count the number of snapshots in `simulation_path`.

!!! note

    This function count the number of snapshots, no the number of snapshot files. So if each snapshot is made of more than one files, the count will not change.

# Arguments

  - `simulation_path::String`: Path to the simulation directory, set in the code variable `OutputDir`.
  - `warnings::Bool=true`: If a warning will be raised when no snapshot files or folders are found.

# Returns

  - The number of snapshots.
"""
function countSnapshot(simulation_path::String; warnings::Bool=true)::Int

    (
        isdir(simulation_path) ||
        throw(ArgumentError("countSnapshot: $(simulation_path) does not exist as a directory"))
    )

    # Get the full list of paths to every snapshot in `simulation_path`
    path_list = [
        glob("*/*/$(SNAP_BASENAME)_*", simulation_path)
        glob("*/$(SNAP_BASENAME)_*", simulation_path)
        glob("$(SNAP_BASENAME)_*", simulation_path)
    ]

    # Check for an empty folder
    if isempty(path_list)

        (
            !warnings ||
            @warn("countSnapshot: I could not find any file named $(SNAP_BASENAME)_*.hdf5 \
            within $(simulation_path), or any of its subfolders")
        )

        return 0

    end

    if readSnapHeader(first(path_list)).num_files > 1
        # If there are multiple files per snapshot, get the path to the snapshot directory
        map!(dirname, path_list, path_list)
        # Delete duplicates
        unique!(path_list)
    end

    return length(path_list)

end
