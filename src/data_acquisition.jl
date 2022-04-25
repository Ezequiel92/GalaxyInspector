####################################################################################################
# Data acquisition functions
####################################################################################################

"""
    getSnapshotPaths(
        base_name::String,
        source_path::String,
    )::Dict{String,Vector{String}}

Find the paths to the GADGET output files and group them by snapshot.

# Arguments
- `base_name::String`: Base name of the snapshot files, set in the GADGET variable `SnapshotFileBase`.
- `source_path::String`: Path to the directory containing the snapshot files, set in the GADGET 
  variable `OutputDir`.

# Returns
- A dictionary with two entries:
  - `"snap_numbers"` ⟹ The numbers that characterize each snapshot.
  - `"snap_paths"` ⟹ The full paths to the snapshot files.

"""
function getSnapshotPaths(base_name::String, source_path::String)::Dict{String,Vector{String}}

    # Get the full list of paths to every GADGET file in `source_path`
    path_list = [
        glob("**/" * base_name * "_*", source_path)
        glob(base_name * "_*", source_path)
    ]

    # Data availability check
    !isempty(path_list) || error("I couldn't find any snapshots in $source_path.")

    # Get the number of files per snapshot
    num_files = read_header(first(path_list)).num_files

    if num_files > 1
        # If there are multiple files per snapshot, delete the trailing '.n'
        map!(x -> rsplit(x, "."; limit=2)[1], path_list, path_list)
        # And delete duplicates
        unique!(path_list)
    end

    # Get the numbers that characterize each snapshot
    number_list = map(x -> rsplit(x, base_name * '_'; limit=2)[2], path_list)

    return Dict("snap_numbers" => number_list, "snap_paths" => normpath.(path_list))

end

"""
    makeSourceTable(
        source::Vector{Dict{String,Vector{String}}},
        idx::IndexType,
        t_unit::Unitful.Units; 
        <keyword arguments>
    )::DataFrame

Construct a dataframe with the path, the time stamp, the number, and the global index of 
each snapshot in a series of simulations.

# Arguments 
- `source::Vector{Dict{String,Vector{String}}}`: Vector with the output of [`getSnapshotPaths`](@ref) 
  for every simulation.
- `idx::IndexType`: Indexing of the simulation, i.e. which snapshots will be read. It can be an 
  integer (a single snapshot), a vector of integers (several snapshots), an UnitRange (e.g. 5:13) 
  or : (every snapshot will be read).
- `t_unit::Unitful.Units`: Target unit for the time-stamps.
- `sim_cosmo::Bool = false`: If the simulation is cosmological, 
  - `false` ⟶ Newtonian simulation (`ComovingIntegrationOn` = 0).
  - `true` ⟶ Cosmological simulation (`ComovingIntegrationOn` = 1).

# Returns
- A dataframe, with the following columns:
  - `:id` ⟶ global index of each snapshot (it should be the same across simulations)
  - `:time_stamps` ⟶ timestamp of each snapshot (it should be the same across simulations)
  - `:numbers` ⟶ numbers in the name of each snapshot (it should be the same across simulations)
  - All the rest are the full paths to the snapshots, one column per simulation

"""
function makeSourceTable(
    source::Vector{Dict{String,Vector{String}}},
    idx::IndexType,
    t_unit::Unitful.Units,
    sim_cosmo::Bool=false,
)::DataFrame

    paths = [sim["snap_paths"] for sim in source]
    rows = [[1:length(p);] for p in paths]
    labels = ["sim_$i" for i = 1:length(paths)]

    source_table = unstack(flatten(DataFrame(l=labels, p=paths, id=rows), [:p, :id]), :l, :p)

    # Get the time-stamps of every snapshot for the longest running simulation
    # Assumes that all the simulations have the snapshots taken at the same clock time
    time_stamps = uconvert.(
        t_unit,
        computeTimeSeries(longest(paths), "clock_time", nothing, sim_cosmo),
    )
    # Add the column with the time-stamps
    insertcols!(source_table, 2, :time_stamps => time_stamps)

    # Add the column with the numbers in the snapshot names
    numbers = [sim["snap_numbers"] for sim in source]
    insertcols!(source_table, 3, :numbers => longest(numbers))

    return DataFrame(source_table[idx, :])

end

"""
    getTemperature(
        file_path::String; 
        <keyword arguments>
    )::Vector{<:RealOrQty}

Get the temperature of the gas particles in a given snapshot.

# Arguments
- `file_path::String`: Path to the snapshot file or folder (for multiple-file snapshots).
- `sim_cosmo::Bool = false`: If the simulation is cosmological, 
  - `false` ⟶ Newtonian simulation (`ComovingIntegrationOn` = 0).
  - `true` ⟶ Cosmological simulation (`ComovingIntegrationOn` = 1).
- `filter_function::Union{Function,Nothing} = nothing`: A function with the signature: 

  `foo(file_path::String)::Vector{Int64}`
  
  It indicates which particles will be read, taking the file path to a snapshot and returning the 
  list of indices of the selected particles. If set to `nothing`, then no particles are filtered. 
  See the [GadgetIO.jl](https://ludwigboess.github.io/GadgetIO.jl/stable/read_snapshots/#Filter-functions) documentation for examples.
- `warnings::Bool = true`: If a warning will be given when the data is not as expected, but the 
  function can still run using sane defaults.

# Returns
- The temperature of the gas particles.

"""
function getTemperature(
    file_path::String;
    sim_cosmo::Bool=false,
    filter_function::Union{Function,Nothing}=nothing,
    warnings::Bool=true
)::Vector{<:RealOrQty}

    blocks = ["Z", "MASS", "U", "NE"]

    # Return an empty array is some block is missing
    if any(x -> !block_present(GIO.select_file(file_path, 0), x), blocks)
        if warnings
            # Get the indices of the missing blocks
            idx_missing = findall(
                x -> !block_present(GIO.select_file(file_path, 0), x),
                blocks,
            )
            @warn(
                "The blocks $(blocks[idx_missing]), which are needed to compute the temperature, \
                are missing."
            )
        end
        return Float32[]
    end

    # Set a default filter function if needed
    if filter_function === nothing
        filter_function = x -> passAll(x, :gas)
    end

    # Get the data from snapshot 
    data = read_blocks_over_all_files(
        file_path,
        blocks;
        filter_function,
        parttype=ParticleType[:gas],
        verbose=false
    )

    header = read_header(file_path)

    # Set units
    metals = data["Z"] * internalUnits("Z", header; sim_cosmo)
    mass = data["MASS"] * internalUnits("MASS", header; sim_cosmo)
    internal_energy = data["U"] * internalUnits("U", header; sim_cosmo)
    electron_fraction = data["NE"] * internalUnits("NE", header; sim_cosmo)

    # Compute the temperature
    return computeTemperature(metals, mass, internal_energy, electron_fraction)

end

"""
    getRawData(
        file_path::String,
        type::Symbol,
        block::String; 
        <keyword arguments>
    )::Dict{Symbol,Dict{String,VecOrMat{<:RealOrQty}}}

Extract the raw data from a snapshot file, for a given type of particle and data block.

# Arguments
- `file_path::String`: Path to the snapshot file or folder (for multiple-file snapshots).
- `type::Symbol`: Particle type, for which the data will be extracted. The posibilities are given 
  by [`ParticleType`](@ref) in `src/constants.jl`.
- `block::String`: Data block to be extracted.
- `sim_cosmo::Bool = false`: If the simulation is cosmological, 
  - `false` ⟶ Newtonian simulation (`ComovingIntegrationOn` = 0).
  - `true` ⟶ Cosmological simulation (`ComovingIntegrationOn` = 1).
- `filter_function::Union{Function,Nothing} = nothing`: A function with the signature: 

  `foo(file_path::String)::Vector{Int64}`
  
  It indicates which particles will be read, taking the file path to a snapshot and returning the 
  list of indices of the selected particles. If set to `nothing`, then no particles are filtered. 
  See the [GadgetIO.jl](https://ludwigboess.github.io/GadgetIO.jl/stable/read_snapshots/#Filter-functions) documentation for examples.
- `warnings::Bool = true`: If a warning will be given when the data is not as expected, but the 
  function can still run using sane defaults.

# Returns
- A vector or matrix with the raw data.

"""
function getRawData(
    file_path::String,
    type::Symbol,
    block::String;
    sim_cosmo::Bool=false,
    filter_function::Union{Function,Nothing}=nothing,
    warnings::Bool=true
)::VecOrMat{<:RealOrQty}

    if block == "TEMP" && type == :gas

        return getTemperature(file_path; sim_cosmo, filter_function, warnings)

    elseif !block_present(GIO.select_file(file_path, 0), block)

        # If the block is missing from the snapshot
        !warnings || @warn(
            "There is no block '$block', for the particle type '$type', in the snapshot $file_path."
        )

        return Float32[]

    else

        # Set a default filter function if needed
        if filter_function === nothing
            filter_function = x -> passAll(x, type)
        end

        return read_blocks_over_all_files(
            file_path,
            [block];
            filter_function,
            parttype=ParticleType[type],
            verbose=false
        )[block] * internalUnits(block, read_header(file_path); sim_cosmo)

    end

end

"""
    getSnapshotData(
        file_path::String,
        type::Symbol,
        block::String; 
        <keyword arguments>
    )::Dict{Symbol,Dict{String,VecOrMat{<:RealOrQty}}}

Extract the raw data from a snapshot file.

Use this method when you want one block from one type of particle.

# Arguments
- `file_path::String`: Path to the snapshot file or folder (for multiple-file snapshots).
- `type::Symbol`: Particle type, for which the data will be extracted. The posibilities are given 
  by [`ParticleType`](@ref) in `src/constants.jl`.
- `block::String`: Data block to be extracted.
- `sim_cosmo::Bool = false`: If the simulation is cosmological, 
  - `false` ⟶ Newtonian simulation (`ComovingIntegrationOn` = 0).
  - `true` ⟶ Cosmological simulation (`ComovingIntegrationOn` = 1).
- `filter_function::Union{Function,Nothing} = nothing`: A function with the signature: 

  `foo(file_path::String)::Vector{Int64}`
  
  It indicates which particles will be read, taking the file path to a snapshot and returning the 
  list of indices of the selected particles. If set to `nothing`, then no particles are filtered. 
  See the [GadgetIO.jl](https://ludwigboess.github.io/GadgetIO.jl/stable/read_snapshots/#Filter-functions) documentation for examples.
- `warnings::Bool = true`: If a warning will be given when the data is not as expected, but the 
  function can still run using sane defaults.

# Returns
- A dictionary with the following structue:
  - `type` => (`block` ⟹ Vector or matrix with the raw data).

"""
function getSnapshotData(
    file_path::String,
    type::Symbol,
    block::String;
    sim_cosmo::Bool=false,
    filter_function::Union{Function,Nothing}=nothing,
    warnings::Bool=true
)::Dict{Symbol,Dict{String,VecOrMat{<:RealOrQty}}}

    return Dict(
        type => Dict(
            block => getRawData(file_path, type, block; sim_cosmo, filter_function, warnings),
        ),
    )

end

"""
    getSnapshotData(
        file_path::String,
        type::Symbol,
        blocks::Vector{String}; 
        <keyword arguments>
    )::Dict{Symbol,Dict{String,VecOrMat{<:RealOrQty}}}

Extract the raw data from a snapshot file.

Use this method when you want many blocks from the same type of particle.

# Arguments
- `file_path::String`: Path to the snapshot file or folder (for multiple-file snapshots).
- `type::Symbol`: Particle type, for which the data will be extracted. The posibilities are given 
  by [`ParticleType`](@ref) in `src/constants.jl`.
- `blocks::Vector{String}`: Data blocks to be extracted.
- `sim_cosmo::Bool = false`: If the simulation is cosmological, 
  - `false` ⟶ Newtonian simulation (`ComovingIntegrationOn` = 0).
  - `true` ⟶ Cosmological simulation (`ComovingIntegrationOn` = 1).
- `filter_function::Union{Function,Nothing} = nothing`: A function with the signature: 

  `foo(file_path::String)::Vector{Int64}`
  
  It indicates which particles will be read, taking the file path to a snapshot and returning the 
  list of indices of the selected particles. If set to `nothing`, then no particles are filtered. 
  See the [GadgetIO.jl](https://ludwigboess.github.io/GadgetIO.jl/stable/read_snapshots/#Filter-functions) documentation for examples.
- `warnings::Bool = true`: If a warning will be given when the data is not as expected, but the 
  function can still run using sane defaults.

# Returns
- A dictionary with the following structue:
  - `type` => (`block` ⟹ data, `block` ⟹ data, `block` ⟹ data, ...), where data is a vector or 
    matrix with the raw data

"""
function getSnapshotData(
    file_path::String,
    type::Symbol,
    blocks::Vector{String};
    sim_cosmo::Bool=false,
    filter_function::Union{Function,Nothing}=nothing,
    warnings::Bool=true
)::Dict{Symbol,Dict{String,VecOrMat{<:RealOrQty}}}

    return Dict(
        type => Dict(
            block => getRawData(file_path, type, block; sim_cosmo, filter_function, warnings)
            for block in blocks
        ),
    )

end

"""
    getSnapshotData(
        file_path::String,
        type_blocks::Dict{Symbol,<:Union{String,Vector{String}}}; 
        <keyword arguments>
    )::Dict{Symbol,Dict{String,VecOrMat{<:RealOrQty}}}

Extract the raw data from a snapshot file.

Use this method when you want one or more blocks from more than one type of particle, not 
necessarily the same from each one.

# Arguments
- `file_path::String`: Path to the snapshot file or folder (for multiple-file snapshots).
- `type_blocks::Dict{Symbol,<:Union{String,Vector{String}}}`: Dictionary where the keys are the 
  particle types (the posibilities are given by [`ParticleType`](@ref) in `src/constants.jl`), and the values are the data blocks 
  (or lists of data blocks) to be extracted, for each type.
- `sim_cosmo::Bool = false`: If the simulation is cosmological, 
  - `false` ⟶ Newtonian simulation (`ComovingIntegrationOn` = 0).
  - `true` ⟶ Cosmological simulation (`ComovingIntegrationOn` = 1).
- `filter_function::Union{Function,Nothing} = nothing`: A function with the signature: 

  `foo(file_path::String)::Vector{Int64}`
  
  It indicates which particles will be read, taking the file path to a snapshot and returning the 
  list of indices of the selected particles. If set to `nothing`, then no particles are filtered. 
  See the [GadgetIO.jl](https://ludwigboess.github.io/GadgetIO.jl/stable/read_snapshots/#Filter-functions) documentation for examples.
- `warnings::Bool = true`: If a warning will be given when the data is not as expected, but the 
  function can still run using sane defaults.

# Returns
- A dictionary with the following structue:
  - `type` => (`block` ⟹ data, `block` ⟹ data, `block` ⟹ data, ...), where data is a vector or 
    matrix with the raw data

"""
function getSnapshotData(
    file_path::String,
    type_blocks::Dict{Symbol,<:Union{String,Vector{String}}};
    sim_cosmo::Bool=false,
    filter_function::Union{Function,Nothing}=nothing,
    warnings::Bool=true
)::Dict{Symbol,Dict{String,VecOrMat{<:RealOrQty}}}

    return merge(
        [
            getSnapshotData(file_path, type, blocks; sim_cosmo, filter_function, warnings) for
            (type, blocks) in type_blocks
        ]...,
    )

end

"""
    getSnapshotData(
        file_path::String,
        types::Vector{Symbol},
        blocks::Union{String,Vector{String}}; 
        <keyword arguments>
    )::Dict{Symbol,Dict{String,VecOrMat{<:RealOrQty}}}

Extract the raw data from a snapshot file.

Use this method when you want the same block (or list of blocks) from many types of particles.

# Arguments
- `file_path::String`: Path to the snapshot file or folder (for multiple-file snapshots).
- `types::Vector{Symbol}`: Particle types, for which the data will be extracted. The posibilities 
  are given by [`ParticleType`](@ref) in `src/constants.jl`.
- `blocks::Union{String,Vector{String}}`: Data block (or list of data blocks) to be extracted.
- `sim_cosmo::Bool = false`: If the simulation is cosmological, 
  - `false` ⟶ Newtonian simulation (`ComovingIntegrationOn` = 0).
  - `true` ⟶ Cosmological simulation (`ComovingIntegrationOn` = 1).
- `filter_function::Union{Function,Nothing} = nothing`: A function with the signature: 

  `foo(file_path::String)::Vector{Int64}`
  
  It indicates which particles will be read, taking the file path to a snapshot and returning the 
  list of indices of the selected particles. If set to `nothing`, then no particles are filtered. 
  See the [GadgetIO.jl](https://ludwigboess.github.io/GadgetIO.jl/stable/read_snapshots/#Filter-functions) documentation for examples.
- `warnings::Bool = true`: If a warning will be given when the data is not as expected, but the 
  function can still run using sane defaults.

# Returns
- A dictionary with the following structue:
  - `type` => (`block` ⟹ data, `block` ⟹ data, `block` ⟹ data, ...), where data is a vector or 
    matrix with the raw data

"""
function getSnapshotData(
    file_path::String,
    types::Vector{Symbol},
    blocks::Union{String,Vector{String}};
    sim_cosmo::Bool=false,
    filter_function::Union{Function,Nothing}=nothing,
    warnings::Bool=true
)::Dict{Symbol,Dict{String,VecOrMat{<:RealOrQty}}}

    return getSnapshotData(
        file_path,
        Dict(type => blocks for type in types);
        sim_cosmo,
        filter_function,
        warnings
    )

end

"""
    getSfrFile(
        source_path::String,
        header::SnapshotHeader; 
        <keyword arguments>
    )::Dict{Int32,VecOrMat{<:RealOrQty}}

Get the data from the `sfr.txt` file.

# Arguments
- `source_path::String`: Path to the directory containing the `sfr.txt` file.
- `header::SnapshotHeader`: Header of one snapshot from the same simulation of the `sfr.txt` file.
- `sim_cosmo::Bool = false`: If the simulation is cosmological, 
  - `false` ⟶ Newtonian simulation (`ComovingIntegrationOn` = 0).
  - `true` ⟶ Cosmological simulation (`ComovingIntegrationOn` = 1).

# Returns
- A dictionary with as many entries as columns has the `sfr.txt` file.
  In the version produced by GADGET3, there are six columns:
  - `1` ⟹ Time.  
  - `2` ⟹ Total mass (probability).  
  - `3` ⟹ SFR (original GADGET).  
  - `4` ⟹ SFR (probability).  
  - `5` ⟹ Real total mass.  
  - `6` ⟹ Real SFR. 

"""
function getSfrFile(
    source_path::String,
    header::SnapshotHeader;
    sim_cosmo::Bool=false
)::Dict{Int32,VecOrMat{<:RealOrQty}}

    # Load the data from the `sfr.txt` file
    file_data = readdlm(joinpath(source_path, "sfr.txt"), Float64)

    head = deepcopy(header)
    times = first(eachcol(file_data))
    headers = Vector{SnapshotHeader}(undef, length(times))
    for (i, t) in enumerate(times)
        if sim_cosmo
            # For cosmological simulations, only the scale factor in each header changes
            # from row to row, so take it from the first column of `sfr.txt`  
            setfield!(head, :time, t)
        end
        headers[i] = head
    end

    return Dict(
        n => column .* internalUnits.("SFRTXT_COL" * string(n), headers; sim_cosmo) for
        (n, column) in enumerate(eachcol(file_data))
    )

end

"""
    getCpuFile(
        source_path::String, 
        targets::Vector{String}; 
        <keyword arguments>
    )::Dict{String,Matrix{Float64}}

Get the CPU utilization of several processes from the `cpu.txt` file.

For each process in `targets` a matrix with all the CPU usage data (as percentages of 
total CPU time) is returned.

# Arguments
- `source_path::String`: Path to the directory containing the `cpu.txt` file.
- `targets::Vector{String}`: Target processes.
- `step::Int64 = 1`: Step used to traverse the CPU cycles, i.e. one every `step` cycles is returned.
- `warnings::Bool = true`: If a warning will be given when the data is not as expected, but the 
  function can still run using sane defaults.

# Returns
- A dictionary with one entry per target process:
  - process ⟹ Matrix with CPU cycles in its first column, and CPU usage (as a percentage) in its 
    second column.   

"""
function getCpuFile(
    source_path::String,
    targets::Vector{String};
    step::Int64=1,
    warnings::Bool=true
)::Dict{String,Matrix{Float64}}

    # Load the data from the `cpu.txt` file
    file = eachline(joinpath(source_path, "cpu.txt"))

    # Auxiliary dictionary
    data_aux = Dict(target => Float64[] for target in targets)
    # Output dictionary
    data_out = Dict(target => Array{Float64}(undef, 0, 2) for target in targets)

    # Save the percentages of CPU usage for the target processes
    for line in file

        # Ignore empty lines
        !isempty(line) || continue

        columns = split(line)

        # Ignore title lines
        !(columns[1] == "Step") || continue

        if columns[1] in targets
            push!(data_aux[columns[1]], parse(Float64, columns[3][1:end-1]))
        end

    end

    if any(isempty.(values(data_aux)))
        !warnings || @warn("I cound't find some/all of the target rows. Check the spelling!.")
    end

    # Data reduction
    for (key, entry) in data_aux

        !isempty(entry) || continue

        l_e = length(entry)

        if 1 < step < l_e
            data_out[key] = hcat(1:step:l_e, entry[1:step:end])
        else
            data_out[key] = hcat(1:l_e, entry)
            if step > l_e
                !warnings || @warn(
                    "`step` = $step is bigger than the number of CPU cycles. Make it smaller!."
                )
            end
        end

    end

    return data_out

end

"""
    getCpuFile(
        source_path::String, 
        target::String; 
        <keyword arguments>
    )::Dict{String,Matrix{Float64}}

Get the CPU consumption of the `target` process, from the `cpu.txt` file.

# Arguments
- `source_path::String`: Path to the directory containing the `cpu.txt` file.
- `target::String`: Target process.
- `step::Int64 = 1`: Step used to traverse the CPU cycles, i.e. one every `step` cycles is returned.
- `warnings::Bool = true`: If a warning will be given when the data is not as expected, but the 
  function can still run using sane defaults.

# Returns
- A dictionary with one entry:
  - process ⟹ Matrix with CPU cycles in its first column, and CPU usage (as a percentage) in 
    its second column.   
    
"""
function getCpuFile(
    source_path::String,
    target::String;
    step::Int64=1,
    warnings::Bool=true
)::Dict{String,Matrix{Float64}}

    return getCpuFile(source_path, [target]; step, warnings)

end


"""
    getMollá2015(
        source_path::String,
    )::Tuple{DataFrame,Dict{String,Tuple{Unitful.Units,Bool}}}

Get the experimental data compiled by Mollá et al. (2015), from the `data.csv` file. The quantities 
in the file are

- "R" ⟶ Distance from the galactic center.
- "ΣHI" ⟶ Atomic Hydrogen area density.
- "ΣHII" ⟶ Molecular Hydrogen area density.
- "logΣ*" ⟶ Logarithm of the stellar mass area density.
- "logΣsfr" ⟶ Logarithm of the SFR area density.
- "C/H" ⟶ 12 + log(C/H), where C and H are the concentration of Carbon and Hydrogen respectively.
- "N/H" ⟶ 12 + log(N/H), where N and H are the concentration of Nitrogen and Hydrogen respectively.
- "O/H" ⟶ 12 + log(O/H), where O and H are the concentration of Oxygen and Hydrogen respectively.

Every quantity but "R" has its error and is measured as a 2D profile.

# Source

M. Mollá et al. (2015). *Galactic chemical evolution: stellar yields and the initial mass function.*
Monthly Notices of the Royal Astronomical Society, **451(4)**, 3693–3708.
[https://doi.org/10.1093/mnras/stv1102](https://doi.org/10.1093/mnras/stv1102)

# Arguments
- `source_path::String`: Path to the directory containing the `data.csv` file.
    
# Returns
- A Tuple with two elements:
  - A dataframe with all the data.
  - A dictionary where the keys are the quantity, and the values are the units and if the unit is 
   inside a logarithm.
           
"""
function getMollá2015(
    source_path::String,
)::Tuple{DataFrame,Dict{String,Tuple{Unitful.Units,Bool}}}

    # Load the data from the file `data.csv`
    data = CSV.read(
        joinpath(source_path, "data.csv"),
        DataFrame;
        header=[
            "R",
            "ΣHI",
            "ΣHI_error",
            "ΣHII",
            "ΣHII_error",
            "logΣ*",
            "logΣ*_error",
            "logΣsfr",
            "logΣsfr_error",
            "C/H",
            "C/H_error",
            "N/H",
            "N/H_error",
            "O/H",
            "O/H_error",
        ],
        types=fill(Float64, 15)
    )

    # Quantity => (unit, unit inside a logarithm?)
    units = Dict(
        "R" => (u"kpc", false),
        "ΣHI" => (u"Msun * pc^(-2)", false),
        "ΣHII" => (u"Msun * pc^(-2)", false),
        "logΣ*" => (u"Msun * pc^(-2)", true),
        "logΣsfr" => (u"Msun * pc^(-2) * Gyr^(-1)", true),
        "C/H" => (Unitful.NoUnits, false),
        "N/H" => (Unitful.NoUnits, false),
        "O/H" => (Unitful.NoUnits, false),
    )

    return data, units

end
