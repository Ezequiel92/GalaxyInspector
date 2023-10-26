####################################################################################################
# Opinionated convenience functions.
####################################################################################################

"""
    snapshotReport(
        simulation_paths::Vector{String},
        slice_n::Int;
        <keyword arguments>
    )::Nothing

Write a text file with information about a given snapshot.

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`. One text file will be printed for each simulation.
  - `slice_n::Int`: Selects which snapshot to plot, starts at 1 and is independent of the number in the file name. If every snapshot is present, `slice_n` = filename_number + 1.
  - `output_path::String="./"`: Path to the output folder.
  - `filter_mode::Symbol=:all`: Which cells/particles will be considered in the "filtered" section of the report. The options are:

      + `:all`             -> Consider every cell/particle within the simulation box.
      + `:halo`            -> Consider only the cells/particles that belong to the main halo.
      + `:subhalo`         -> Consider only the cells/particles that belong to the main subhalo.
      + `:sphere`          -> Consider only the cell/particle inside a sphere with radius `FILTER_R` (see `./src/constants.jl`).
      + `:stellar_subhalo` -> Consider only the cells/particles that belong to the main subhalo.
  - `halo_idx::Int=1`: Index of the target halo (FoF group). Starts at 1.
  - `subhalo_rel_idx::Int=1`: Index of the target subhalo (subfind), relative the target halo. Starts at 1.
  - `warnings::Bool=true`: If a warning will be given when there is missing files.
"""
function snapshotReport(
    simulation_paths::Vector{String},
    slice_n::Int;
    output_path::String="./",
    filter_mode::Symbol=:all,
    halo_idx::Int=1,
    subhalo_rel_idx::Int=1,
    warnings::Bool=true,
)::Nothing

    @inbounds for simulation_path in simulation_paths

        ############################################################################################
        # Load the relevant values and check for missing files
        ############################################################################################

        # Make a dataframe with the following columns:
        #   - 1. DataFrame index
        #   - 2. Number in the file name
        #   - 3. Scale factor
        #   - 4. Redshift
        #   - 5. Physical time
        #   - 6. Lookback time
        #   - 7. Snapshot path
        #   - 8. Group catalog path
        simulation_table = makeSimulationTable(simulation_path)

        # Get the number in the filename
        snap_n = safeSelect(simulation_table[!, :numbers], slice_n; warnings)

        # Check that after slicing there is one snapshot left
        (
            !isempty(snap_n) ||
            throw(ArgumentError("snapshotReport: There are no snapshots with `slice_n` = \
            $(slice_n), the contents of $(simulation_path) are: \n$(simulation_table)"))
        )

        # Find the target row
        snapshot_row = filter(:numbers => ==(lpad(snap_n, 3, "0")), simulation_table)

        # Construct the file name of the target snapshot
        snapshot_filename = "\"$(SNAP_BASENAME)_$(lpad(snap_n, 3, "0"))\""

        # Select the path to the target snapshot
        snapshot_path = snapshot_row[1, :snapshot_paths]

        (
            !ismissing(snapshot_path) ||
            throw(ArgumentError("snapshotReport: The snapshot $(snapshot_filename) is missing \
            in $(simulation_path)"))
        )

        # Select the path to the target group catalog file
        groupcat_path = snapshot_row[1, :groupcat_paths]

        # Check if the simulation is cosmological
        cosmological = isCosmological(snapshot_path)

        # Select the physical time since the Big Bang
        physical_time = round(ustrip(u"Gyr", snapshot_row[1, :physical_times]), digits=2)

        # Select the ordinal index of the target snapshot
        o_idx = snapshot_row[1, :ids]

        # Compute the number of snapshots in the folder
        snapshot_length = count(!ismissing, simulation_table[!, :snapshot_paths])

        # Compute the number of group catalog files in the folder
        groupcat_length = count(!ismissing, simulation_table[!, :groupcat_paths])

        # Get the snapshot header
        snapshot_header = readSnapHeader(snapshot_path)

        ############################################################################################
        # Print the report header
        ############################################################################################

        # Create the output file
        filename = "$(SNAP_BASENAME)_$(lpad(snap_n, 3, "0"))-of-$(basename(simulation_path))"
        file = open(joinpath(mkpath(output_path), "report-for-$(filename).txt"), "w")

        println(file, "#"^100)
        println(file, "\nSimulation name:  $(basename(simulation_path))")
        println(file, "Physical time:    $(physical_time) Gyr")

        if cosmological
            # For cosmological simulations print the scale factor and the redshift
            scale_factor = round(snapshot_row[1, :scale_factors], digits=3)
            redshift = round(snapshot_row[1, :redshifts], digits=3)

            println(file, "Scale factor:     $(scale_factor)")
            println(file, "Redshift:         $(redshift)")
            println(file, "Cosmological:     Yes")
        else
            println(file, "Cosmological:     No")
        end

        if PHYSICAL_UNITS
            println(file, "Report units:     Physical\n")
        else
            println(file, "Report units:     Comoving\n")
        end

        println(file, "#"^100)
        println(file, "\nSnapshot number:  $(o_idx) (of $(snapshot_length))")
        println(file, "Snapshot path:    $(snapshot_path)\n")

        if !ismissing(groupcat_path)

            println(file, "Subfind number:   $(o_idx) (of $(groupcat_length))")
            println(file, "Subfind path:     $(groupcat_path)\n")

        end

        println(file, "#"^100)

        ############################################################################################
        # Read the data in the snapshot and group catalog file
        ############################################################################################

        # Detect which components are present in the snapshot
        if isfile(snapshot_path)
            file_path = snapshot_path
        else
            file_path = minimum(glob("$(SNAP_BASENAME)_*.*.hdf5", snapshot_path))
        end
        component_list = h5open(file_path, "r") do snapshot
            filter(
                in([:gas, :halo, :stars, :black_hole]),
                [get(PARTICLE_TYPE, key, nothing) for key in keys(snapshot)],
            )
        end

        # Select the filter function, translation and request dictionary
        filter_function, translation, _, request = selectFilter(
            filter_mode,
            mergeRequests(
                Dict(component => ["POS ", "MASS", "VEL "] for component in component_list),
                Dict(:gas => ["NHP ", "NH  ", "PRES", "FRAC"]),
            )
        )

        # Read the necessary snapshot data
        if filter_mode ∈ [:halo, :subhalo, :stellar_subhalo]

            # Check that the group catalog data is available
            if !ismissing(groupcat_path) && isSubfindActive(groupcat_path)

                # Get the group catalog header
                groupcat_header = readGroupCatHeader(groupcat_path; warnings)

                # Create a metadata dictionary
                metadata = Dict(
                    :snap_data => Snapshot(
                        snapshot_path,
                        1,
                        1,
                        0.0u"yr",
                        0.0u"yr",
                        0.0,
                        0.0,
                        snapshot_header,
                    ),
                    :gc_data => GroupCatalog(groupcat_path, groupcat_header),
                )

                # Create the data dictionary
                data_dict = merge(
                    metadata,
                    readSnapshot(snapshot_path, request; warnings),
                    readGroupCatalog(groupcat_path, snapshot_path, request; warnings),
                )

            else

                throw(ArgumentError("snapshotReport: You asked for a filter base on the \
                halos/subhalos, but I could not find a valid group catalog file"))

            end

        else

            # Create the data dictionary
            data_dict = readSnapshot(snapshot_path, request)

        end

        ############################################################################################
        # Print the global properties of the simulation
        ############################################################################################

        println(file, "\nGlobal properties (full simulation box):")

        ############################################################################################
        # Print the total number of cells/particles for each component
        ############################################################################################

        println(file, "\n\tCell/particle number:\n")

        total_count = 0
        for component in component_list

            count = getindex(getfield(snapshot_header, :num_total), PARTICLE_INDEX[component] + 1)

            total_count += count

            title = "$(PARTICLE_NAMES[component]):"
            title *= " "^(24 - length(title))

            println(file, "\t\t$(title)$(count)")

        end

        println(file, "\n\t\tTotal count:            $(total_count)\n")

        ############################################################################################
        # Print the total mass of each component
        ############################################################################################

        println(file, "\tMasses:\n")

        total_mass = 0.0u"Msun"
        for component in component_list

            mass = sum(data_dict[component]["MASS"]; init=0.0u"Msun")

            total_mass += mass

            title = "$(PARTICLE_NAMES[component]):"
            title *= " "^(24 - length(title))

            println(file, "\t\t$(title)$(round(typeof(1.0u"Msun"), mass, sigdigits=3))")

        end

        println(
            file,
            "\n\t\tTotal mass:             $(round(typeof(1.0u"Msun"), total_mass, sigdigits=3))\n",
        )

        ############################################################################################
        # Print the mass of each hydrogen phase
        ############################################################################################

        if :gas in component_list

            println(file, "\tHydrogen masses:\n")

            gas_mass = sum(data_dict[:gas]["MASS"]; init=0.0u"Msun")

            hii_mass = sum(computeIonizedMass(data_dict); init=0.0u"Msun")
            hii_percent = round((hii_mass / gas_mass) * 100, sigdigits=3)

            hi_mass = sum(computeAtomicMass(data_dict); init=0.0u"Msun")
            hi_percent = round((hi_mass / gas_mass) * 100, sigdigits=3)

            h2_mass = sum(computeMolecularMass(data_dict); init=0.0u"Msun")
            h2_percent = round((h2_mass / gas_mass) * 100, sigdigits=3)

            title = "Ionized mass:"
            title *= " "^(24 - length(title))
            println(
                file,
                "\t\t$(title)$(round(typeof(1.0u"Msun"), hii_mass, sigdigits=3)) \
                ($(hii_percent)% of total gas mass)",
            )

            title = "Atomic mass:"
            title *= " "^(24 - length(title))
            println(
                file,
                "\t\t$(title)$(round(typeof(1.0u"Msun"), hi_mass, sigdigits=3)) \
                ($(hi_percent)% of total gas mass)",
            )

            title = "Molecular mass:"
            title *= " "^(24 - length(title))
            println(
                file,
                "\t\t$(title)$(round(typeof(1.0u"Msun"), h2_mass, sigdigits=3)) \
                ($(h2_percent)% of total gas mass)\n",
            )

        end

        ############################################################################################
        # Print the global properties of the simulation after filtering
        ############################################################################################

        # Filter cell/particles
        filterData!(data_dict; filter_function)

        println(file, "#"^100)
        println(file, "\nGlobal properties (filtered box with mode :$(filter_mode)):")

        ############################################################################################
        # Print the number of cells/particles for each component
        ############################################################################################

        println(file, "\n\tCell/particle number:\n")

        total_count = 0
        for component in component_list

            count = length(data_dict[component]["MASS"])

            title = "$(PARTICLE_NAMES[component]):"
            title *= " "^(24 - length(title))

            total_count += count

            println(file, "\t\t$(title)$(count)")

        end

        println(file, "\n\t\tTotal count:            $(total_count)\n")

        ############################################################################################
        # Print the mass of each component
        ############################################################################################

        println(file, "\tMasses:\n")

        total_mass = 0.0u"Msun"
        for component in component_list

            mass = sum(data_dict[component]["MASS"]; init=0.0u"Msun")

            total_mass += mass

            title = "$(PARTICLE_NAMES[component]):"
            title *= " "^(24 - length(title))

            println(file, "\t\t$(title)$(round(typeof(1.0u"Msun"), mass, sigdigits=3))")

        end

        println(
            file,
            "\n\t\tTotal mass:             $(round(typeof(1.0u"Msun"), total_mass, sigdigits=3))\n",
        )

        ############################################################################################
        # Print the mass of each hydrogen phase
        ############################################################################################

        if :gas in component_list

            println(file, "\tHydrogen masses:\n")

            gas_mass = sum(data_dict[:gas]["MASS"]; init=0.0u"Msun")

            hii_mass = sum(computeIonizedMass(data_dict); init=0.0u"Msun")
            hii_percent = round((hii_mass / gas_mass) * 100, sigdigits=3)

            hi_mass = sum(computeAtomicMass(data_dict); init=0.0u"Msun")
            hi_percent = round((hi_mass / gas_mass) * 100, sigdigits=3)

            h2_mass = sum(computeMolecularMass(data_dict); init=0.0u"Msun")
            h2_percent = round((h2_mass / gas_mass) * 100, sigdigits=3)

            title = "Ionized mass:"
            title *= " "^(24 - length(title))
            println(
                file,
                "\t\t$(title)$(round(typeof(1.0u"Msun"), hii_mass, sigdigits=3)) \
                ($(hii_percent)% of total gas mass)",
            )

            title = "Atomic mass:"
            title *= " "^(24 - length(title))
            println(
                file,
                "\t\t$(title)$(round(typeof(1.0u"Msun"), hi_mass, sigdigits=3)) \
                ($(hi_percent)% of total gas mass)",
            )

            title = "Molecular mass:"
            title *= " "^(24 - length(title))
            println(
                file,
                "\t\t$(title)$(round(typeof(1.0u"Msun"), h2_mass, sigdigits=3)) \
                ($(h2_percent)% of total gas mass)\n",
            )

        end

        ############################################################################################
        # Print the center of mass of each component
        ############################################################################################

        println(file, "\tCenter of mass:\n")

        global_cm = computeGlobalCenterOfMass(data_dict)

        for component in component_list

            cm = computeCenterOfMass(data_dict[component]["POS "], data_dict[component]["MASS"])

            title = "$(PARTICLE_NAMES[component]):"
            title *= " "^(24 - length(title))

            println(file, "\t\t$(title)$(round.(ustrip.(u"Mpc", cm), sigdigits=6)) $(u"Mpc")")
            println(file, "\t\tDistance to global CM:  $(sqrt(sum((global_cm - cm).^2)))\n")

        end

        global_cm = round.(ustrip.(u"Mpc", global_cm), sigdigits=6)
        println(file, "\t\tGlobal center of mass:  $(global_cm) $(u"Mpc")\n")

        # Translate the simulation box
        translateData!(data_dict, translation)

        ############################################################################################
        # Print the normalized angular momentum of each component
        ############################################################################################

        println(file, "\tNormalized angular momentum:\n")

        for component in component_list

            L = computeTotalAngularMomentum(
                data_dict[component]["POS "],
                data_dict[component]["VEL "],
                data_dict[component]["MASS"];
            )

            title = "$(PARTICLE_NAMES[component]):"
            title *= " "^(24 - length(title))

            println(file, "\t\t$(title)$(round.(L, sigdigits=3))")

        end

        global_L = round.(computeGlobalAngularMomentum(data_dict), sigdigits=3)

        println(file, "\n\t\tGlobal angular momentum: $(global_L)\n")

        ############################################################################################
        # Print the spin parameter of each component
        ############################################################################################

        println(file, "\tSpin parameter (R = $(FILTER_R)):\n")

        for component in component_list

            λ = computeSpinParameter(
                data_dict[component]["POS "],
                data_dict[component]["VEL "],
                data_dict[component]["MASS"],
            )

            title = "$(PARTICLE_NAMES[component]):"
            title *= " "^(24 - length(title))

            println(file, "\t\t$(title)$(round.(λ, sigdigits=3))")

        end

        global_λ = round.(computeGlobalSpinParameter(data_dict), sigdigits=3)

        println(file, "\n\t\tTotal spin parameter:   $(global_λ)\n")

        ############################################################################################
        # Print the maximum and minimum values of each parameter of the ODEs
        ############################################################################################

        println(file, "\tExtrema of the ODEs ICs and parameters:\n")

        quantities = ["ODIT", "ACIT", "DTIM", "TAUS", "RHOC", "PARZ", "ETAD", "ETAI", "PARR"]
        names = [
            "integration time",
            "accumulated integration time",
            "delta time",
            "τS",
            "cell density",
            "metallicity",
            "ηd",
            "ηi",
            "R",
        ]
        units = [
            UnitfulAstro.Myr,
            UnitfulAstro.Myr,
            UnitfulAstro.Myr,
            UnitfulAstro.Myr,
            u"cm^-3",
            Unitful.NoUnits,
            Unitful.NoUnits,
            Unitful.NoUnits,
            Unitful.NoUnits,
        ]

        for (quantity, unit, name) in zip(quantities, units, names)

            if blockPresent(:gas, quantity, snapshot_path)

                min, max = findQtyExtrema(
                    simulation_path,
                    slice_n,
                    :gas,
                    quantity;
                    f=x -> filter(!isnan, x),
                    warnings,
                )

                println(file, "\t\tMaximum $name: $(round(ustrip(unit, max), sigdigits=5)) $unit")
                println(file, "\t\tMinimum $name: $(round(ustrip(unit, min), sigdigits=5)) $unit\n")

            end

        end

        ############################################################################################
        # Print the properties of the target halo and subhalo
        ############################################################################################

        if !ismissing(groupcat_path) && isSubfindActive(groupcat_path)

            # Check that the requested halo index is within bounds
            n_groups_total = readGroupCatHeader(groupcat_path; warnings).n_groups_total
            (
                0 < halo_idx <= n_groups_total ||
                throw(ArgumentError("snapshotReport: There is only $(n_groups_total) FoF \
                goups in $(simulation_path), so `halo_idx` = $(halo_idx) is outside of bounds"))
            )

            request = Dict(
                :subhalo => ["S_Mass", "S_MassType", "S_LenType", "S_CM", "S_Pos", "S_HalfmassRad"],
                :group => [
                    "G_Mass",
                    "G_MassType",
                    "G_M_Crit200",
                    "G_LenType",
                    "G_Nsubs",
                    "G_CM",
                    "G_Pos",
                    "G_R_Crit200",
                ],
            )

            # Read the necessary data
            gc_data = readGroupCatalog(groupcat_path, snapshot_path, request; warnings)

            # Check that the requested subhalo index is within bounds
            g_n_subs = gc_data[:group]["G_Nsubs"]
            n_subfinds = g_n_subs[halo_idx]
            (
                subhalo_rel_idx <= n_subfinds ||
                throw(ArgumentError("snapshotReport: There is only $(n_subfinds) subhalos \
                for the FoF group $(halo_idx) in $(simulation_path), so `subhalo_rel_idx` \
                = $(subhalo_rel_idx) is outside of bounds"))
            )

            # Compute the number of subhalos and particles up to the last halo before `halo_idx`
            if isone(halo_idx)
                n_subs_floor = 0
            else
                n_subs_floor = sum(g_n_subs[1:(halo_idx - 1)]; init=0)
            end

            # Compute the subhalo absolute index
            subhalo_abs_idx = n_subs_floor + subhalo_rel_idx

            # Load the necessary data
            s_mass          = gc_data[:subhalo]["S_Mass"][subhalo_abs_idx]
            s_mass_type     = gc_data[:subhalo]["S_MassType"][:, subhalo_abs_idx]
            s_len_type      = gc_data[:subhalo]["S_LenType"][:, subhalo_abs_idx]
            s_cm            = gc_data[:subhalo]["S_CM"][:, subhalo_abs_idx]
            s_pos           = gc_data[:subhalo]["S_Pos"][:, subhalo_abs_idx]
            s_half_mass_rad = gc_data[:subhalo]["S_HalfmassRad"][subhalo_abs_idx]
            g_mass          = gc_data[:group]["G_Mass"][halo_idx]
            g_mass_type     = gc_data[:group]["G_MassType"][:, halo_idx]
            g_m_crit_200    = gc_data[:group]["G_M_Crit200"][halo_idx]
            g_len_type      = gc_data[:group]["G_LenType"][:, halo_idx]
            g_n_subs        = gc_data[:group]["G_Nsubs"][halo_idx]
            g_cm            = gc_data[:group]["G_CM"][:, halo_idx]
            g_pos           = gc_data[:group]["G_Pos"][:, halo_idx]
            g_r_crit_200    = gc_data[:group]["G_R_Crit200"][halo_idx]

            ########################################################################################
            # Print the halo properties
            ########################################################################################

            println(file, "#"^71)
            println(file, "NOTE: Stellar particle counts include wind particles from here on out!")
            println(file, "#"^71)

            println(file, "\nHalo $(lpad(halo_idx - 1, 3, "0")) properties:\n")

            ########################################################################################

            println(file, "\tCell/particle number:\n")
            for (i, len) in enumerate(g_len_type)

                component = PARTICLE_NAMES[INDEX_PARTICLE[i - 1]]
                println(file, "\t\t$(component):$(" "^(22 - length(component))) $(len)")

            end

            ########################################################################################

            println(file, "\n\tNumber of subhalos:\n\n\t\t$(g_n_subs)\n")

            ########################################################################################

            println(file, "\tMasses:\n")
            for (i, mass) in enumerate(g_mass_type)

                type_symbol = INDEX_PARTICLE[i - 1]

                if type_symbol ==:stars
                    component = "Stellar/Wind particles"
                else
                    component = PARTICLE_NAMES[type_symbol]
                end

                println(
                    file,
                    "\t\t$(component):$(" "^(22 - length(component))) \
                    $(round(typeof(1.0u"Msun"), mass, sigdigits=3))",
                )

            end

            println(
                file,
                "\n\t\tTotal mass:             $(round.(typeof(1.0u"Msun"), g_mass, sigdigits=3))",
            )

            ########################################################################################

            println(
                file,
                "\n\tCenter of mass:\n\n\t\t$(round.(ustrip.(u"Mpc", g_cm), sigdigits=6)) \
                $(u"Mpc")\n",
            )

            println(
                file,
                "\tPosition of the particle with the minimum gravitational potential energy: \
                \n\n\t\t$(round.(ustrip.(u"Mpc", g_pos), sigdigits=6)) $(u"Mpc")\n",
            )

            separation = sqrt(sum((g_cm - g_pos).^2))
            println(
                file,
                "\tSeparation between the minimum potencial and the global CM: \
                \n\n\t\t$(round(typeof(1.0u"kpc"), separation, sigdigits=6))\n",
            )

            ########################################################################################

            println(
                file,
                "\tTotal mass enclosed in a sphere with a mean density 200 times the critical \
                density:\n\n\t\t$(round(typeof(1.0u"Msun"), g_m_crit_200, sigdigits=3))\n",
            )

            ########################################################################################

            println(
                file,
                "\tRadius of a sphere with a mean density 200 times the critical density: \
                \n\n\t\t$(round(typeof(1.0u"kpc"), g_r_crit_200, sigdigits=4))\n",
            )

            ########################################################################################
            # Print the subhalo properties
            ########################################################################################

            println(file, "#"^100)
            println(
                file,
                "\nSubhalo $(lpad(subhalo_rel_idx - 1, 3, "0")) (of halo \
                $(lpad(halo_idx - 1, 3, "0"))) properties:\n",
            )

            ########################################################################################

            println(file, "\tCell/particle number:\n")
            for (i, len) in enumerate(s_len_type)

                component = PARTICLE_NAMES[INDEX_PARTICLE[i - 1]]
                println(file, "\t\t$(component):$(" "^(22 - length(component))) $(len)")

            end

            ########################################################################################

            println(file, "\n\tMasses:\n")
            for (i, mass) in enumerate(s_mass_type)

                type_symbol = INDEX_PARTICLE[i - 1]

                if type_symbol ==:stars
                    component = "Stellar/Wind particles"
                else
                    component = PARTICLE_NAMES[type_symbol]
                end

                println(
                    file,
                    "\t\t$(component):$(" "^(22 - length(component))) \
                    $(round(typeof(1.0u"Msun"), mass, sigdigits=3))",
                )

            end

            println(
                file,
                "\n\t\tTotal mass:             $(round.(typeof(1.0u"Msun"), s_mass, sigdigits=3))",
            )

            ########################################################################################

            println(
                file,
                "\n\tCenter of mass:\n\n\t\t$(round.(ustrip.(u"Mpc", s_cm), sigdigits=6)) \
                $(u"Mpc")\n",
            )

            println(
                file,
                "\tPosition of the particle with the minimum gravitational potential energy: \
                \n\n\t\t$(round.(ustrip.(u"Mpc", s_pos), sigdigits=6)) $(u"Mpc")\n",
            )

            separation = sqrt(sum((s_cm - s_pos).^2))
            println(
                file,
                "\tSeparation between the minimum potencial and the global CM: \
                \n\n\t\t$(round(typeof(1.0u"kpc"), separation, sigdigits=6))\n",
            )

            ########################################################################################

            println(
                file,
                "\tRadius containing half of the total mass: \
                \n\n\t\t$(round(typeof(1.0u"kpc"), s_half_mass_rad, sigdigits=4))",
            )

        end

        close(file)

    end

    return nothing

end

"""
    simulationReport(
        simulation_paths::Vector{String};
        <keyword arguments>
    )::Nothing

Write a text file with information about a given simulation

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`. One text file will be printed for each simulation.
  - `output_path::String="./"`: Path to the output folder.
  - `warnings::Bool=true`: If a warning will be given when there is missing files.
"""
function simulationReport(
    simulation_paths::Vector{String};
    output_path::String="./",
    warnings::Bool=true,
)::Nothing

    @inbounds for simulation_path in simulation_paths

        # Make a dataframe with the following columns:
        #   - 1. DataFrame index
        #   - 2. Number in the file name
        #   - 3. Scale factor
        #   - 4. Redshift
        #   - 5. Physical time
        #   - 6. Lookback time
        #   - 7. Snapshot path
        #   - 8. Group catalog path
        simulation_table = makeSimulationTable(simulation_path)

        # Compute the number of snapshots in the folder
        snapshot_n = count(!ismissing, simulation_table[!, :snapshot_paths])

        # Check that there is at least one snapshot
        (
            !iszero(snapshot_n) ||
            throw(ArgumentError("simulationReport: There are no snapshots in $(simulation_path)"))
        )

        # Compute the number of group catalog files in the folder
        groupcat_n = count(!ismissing, simulation_table[!, :groupcat_paths])

        # Check if the simulation is cosmological
        cosmological = isCosmological(first(skipmissing(simulation_table[!, :snapshot_paths])))

        ############################################################################################
        # Print the report header
        ############################################################################################

        # Create the output file
        file = open(
            joinpath(mkpath(output_path), "report-for-$(basename(simulation_path)).txt"),
            "w",
        )

        println(file, "#"^100)
        println(file, "\nSimulation name:  $(basename(simulation_path))")

        if cosmological
            println(file, "Cosmological:     Yes")
        else
            println(file, "Cosmological:     No")
        end

        if PHYSICAL_UNITS
            println(file, "Report units:     Physical\n")
        else
            println(file, "Report units:     Comoving\n")
        end

        println(file, "#"^100)
        println(file, "\nNumber of snapshots:       $(snapshot_n)")
        println(file, "Number of group catalogs:  $(groupcat_n)\n")

        # Print the simulation time ranges
        if snapshot_n > 1

            min_pt, max_pt = round.(
                ustrip.(u"Gyr", extrema(simulation_table[!, :physical_times])),
                digits=2,
            )

            println(file, "Physical time range:       $(min_pt) - $(max_pt) Gyr")

            if cosmological

                # For cosmological simulations print the scale factor and the redshift
                min_a, max_a = round.(extrema(simulation_table[!, :scale_factors]), digits=3)
                min_z, max_z = round.(extrema(simulation_table[!, :redshifts]), digits=3)

                println(file, "Scale factor range:        $(min_a) - $(max_a)")
                println(file, "Redshift range:            $(max_z) - $(min_z)\n")

            end

        else

            pt = round.(ustrip.(u"Gyr", extrema(simulation_table[1, :physical_times])), digits=2)


            println(file, "Physical time:             $(pt) Gyr")

            if cosmological

                # For cosmological simulations print the scale factor and the redshift
                a = round.(extrema(simulation_table[1, :scale_factors]), digits=3)
                z = round.(extrema(simulation_table[1, :redshifts]), digits=3)

                println(file, "Scale factor:              $(a)")
                println(file, "Redshift:                  $(z)\n")

            end

        end

        # Set flags to print only the first instance of each condition
        first_star_flag            = false
        first_subhalo_flag         = false
        first_star_in_subhalo_flag = false

        for snapshot_row in eachrow(simulation_table)

            snapshot_path = snapshot_row[:snapshot_paths]

            # Skip this row if there is no snapshot
            !ismissing(snapshot_path) || continue

            # Read the snapshot header
            snapshot_header = readSnapHeader(snapshot_path)

            # Read the group catalog path
            groupcat_path = snapshot_row[:groupcat_paths]

            # Read the different time ticks
            physical_time = round(ustrip(u"Gyr", snapshot_row[:physical_times]), digits=2)
            if cosmological
                # For cosmological simulations read the scale factor and the redshift
                scale_factor = round(snapshot_row[:scale_factors], digits=3)
                redshift = round(snapshot_row[:redshifts], digits=3)
            end

            # Read how many stars there are in this snapshot
            star_number = countStars(snapshot_path)

            # Check if there is subfind information in the group catalog file
            subfind_active = !ismissing(groupcat_path) && isSubfindActive(groupcat_path)

            if subfind_active

                # Read the group catalog header
                groupcat_header = readGroupCatHeader(groupcat_path; warnings)

                # Read the number of halos
                n_groups_total = readGroupCatHeader(groupcat_path; warnings).n_groups_total

                # Make the subfind request
                request = Dict(:subhalo => ["S_LenType", "S_CM", "S_Pos"])

                # Read the necessary data
                gc_data = readGroupCatalog(groupcat_path, snapshot_path, request; warnings)

                # Load the necessary data
                s_cm       = gc_data[:subhalo]["S_CM"][:, 1]
                s_pos      = gc_data[:subhalo]["S_Pos"][:, 1]
                s_len_type = gc_data[:subhalo]["S_LenType"][:, 1]

                # Select the filter function and request dictionary
                filter_function, _, _, request = selectFilter(:subhalo, Dict(:stars => ["MASS"]))

                # Create a metadata dictionary
                metadata = Dict(
                    :snap_data => Snapshot(
                        snapshot_path,
                        1,
                        1,
                        0.0u"yr",
                        0.0u"yr",
                        0.0,
                        0.0,
                        snapshot_header,
                    ),
                    :gc_data => GroupCatalog(groupcat_path, groupcat_header),
                )

                # Create the data dictionary
                data_dict = merge(
                    metadata,
                    readSnapshot(snapshot_path, request; warnings),
                    readGroupCatalog(groupcat_path, snapshot_path, request; warnings),
                )

                filterData!(data_dict; filter_function)

                # Compute the number of stars in the main subhalo
                stellar_n_subhalo = length(data_dict[:stars]["MASS"])

            end

            ########################################################################################
            # First stars
            ########################################################################################

            if star_number > 0 && !first_star_flag

                println(file, "#"^100)
                println(file, "First snapshot with star formation:")
                println(file, "#"^100)

                println(file, "\n\tSnapshot:         $(basename(snapshot_path))")
                println(file, "\tPhysical time:    $(physical_time) Gyr")

                if cosmological
                    # For cosmological simulations print the scale factor and the redshift
                    println(file, "\tScale factor:     $(scale_factor)")
                    println(file, "\tRedshift:         $(redshift)")
                end

                println(file, "\tNumber of stars:  $(star_number)")

                if subfind_active

                    println(file, "\tNumber of halos:  $(n_groups_total)")

                    println(file, "\n\tMain subhalo properties:")

                    println(file, "\n\t\tNumber of stellar particles:  $(stellar_n_subhalo)\n")

                    println(
                        file,
                        "\t\tCenter of mass:\n\n\t\t\t\
                        $(round.(ustrip.(u"Mpc", s_cm), sigdigits=6)) $(u"Mpc")\n",
                    )

                    println(
                        file,
                        "\t\tPosition of the particle with the minimum gravitational potential \
                        energy: \n\n\t\t\t\
                        $(round.(ustrip.(u"Mpc", s_pos), sigdigits=6)) $(u"Mpc")\n",
                    )

                    separation = sqrt(sum((s_cm - s_pos).^2))
                    println(
                        file,
                        "\t\tSeparation between the minimum potencial and the global CM: \
                        \n\n\t\t\t$(round(typeof(1.0u"kpc"), separation, sigdigits=6))\n",
                    )

                else
                    println(file, "\n\tThere is no subfind information for this snapshot!\n")
                end

                first_star_flag = true

            end

            ########################################################################################
            # First subhalos
            ########################################################################################

            if subfind_active && !first_subhalo_flag

                # Read the number of halos
                n_groups_total = readGroupCatHeader(groupcat_path; warnings).n_groups_total

                println(file, "#"^100)
                println(file, "First snapshot with subfind information:")
                println(file, "#"^100)

                println(
                    file,
                    "\n\tSnapshot:         $(basename(snapshot_path))",
                )

                println(file, "\tPhysical time:    $(physical_time) Gyr")

                if cosmological
                    # For cosmological simulations print the scale factor and the redshift
                    println(file, "\tScale factor:     $(scale_factor)")
                    println(file, "\tRedshift:         $(redshift)")
                end

                println(file, "\tNumber of halos:  $(n_groups_total)")

                println(file, "\n\tMain subhalo properties:")

                println(
                    file,
                    "\n\t\tCenter of mass:\n\n\t\t\t$(round.(ustrip.(u"Mpc", s_cm), sigdigits=6)) \
                    $(u"Mpc")\n",
                )

                println(
                    file,
                    "\t\tPosition of the particle with the minimum gravitational potential energy: \
                    \n\n\t\t\t$(round.(ustrip.(u"Mpc", s_pos), sigdigits=6)) $(u"Mpc")\n",
                )

                separation = sqrt(sum((s_cm - s_pos).^2))
                println(
                    file,
                    "\t\tSeparation between the minimum potencial and the global CM: \
                    \n\n\t\t\t$(round(typeof(1.0u"kpc"), separation, sigdigits=6))\n",
                )

                println(file, "#"^55)
                println(file, "NOTE: Stellar particle counts include wind particles!")
                println(file, "#"^55)

                println(file, "\n\t\tCell/particle number:\n")

                for (i, len) in enumerate(s_len_type)
                    component = PARTICLE_NAMES[INDEX_PARTICLE[i - 1]]
                    println(file, "\t\t\t$(component):$(" "^(22 - length(component))) $(len)")
                end

                println(file)

                first_subhalo_flag = true

            end

            ########################################################################################
            # First stars in the main subhalo
            ########################################################################################

            if subfind_active && !first_star_in_subhalo_flag && stellar_n_subhalo > 0

                println(file, "#"^100)
                println(file, "First snapshot with star formation in the main subhalo:")
                println(file, "#"^100)

                println(file, "\n\tSnapshot:         $(basename(snapshot_path))")
                println(file, "\tPhysical time:    $(physical_time) Gyr")

                if cosmological
                    # For cosmological simulations print the scale factor and the redshift
                    println(file, "\tScale factor:     $(scale_factor)")
                    println(file, "\tRedshift:         $(redshift)")
                end

                println(file, "\tNumber of stars:  $(stellar_n_subhalo)\n")

                first_star_in_subhalo_flag = true

            end

            # End the loop when all the conditions have been met
            if first_star_flag && first_subhalo_flag && first_star_in_subhalo_flag
                break
            end

        end

        close(file)

    end

    return nothing

end

"""
    sfrTXT(
        simulation_paths::Vector{String},
        x_quantity::Symbol,
        y_quantity::Symbol;
        <keyword arguments>
    )::Nothing

Plot a time series of the data in the `sfr.txt` file.

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`.
  - `x_quantity::Symbol`: Quantity for the x axis. The possibilities are:

      + `:physical_time` -> Physical time since the Big Bang.
      + `:lookback_time` -> Physical time left to reach the last snapshot.
      + `:scale_factor`  -> Scale factor.
      + `:redshift`      -> Redshift.
  - `y_quantity::Symbol`: Quantity for the y axis. The possibilities are:

      + `:stellar_mass` -> Stellar mass.
      + `:sfr`          -> The star formation rate.
  - `output_path::String="./"`: Path to the output folder.
  - `sim_labels::Union{Vector{String},Nothing}=basename.(simulation_paths)`: Labels for the plot legend, one per simulation. Set it to `nothing` if you don't want a legend.
"""
function sfrTXT(
    simulation_paths::Vector{String},
    x_quantity::Symbol,
    y_quantity::Symbol;
    output_path::String="./",
    sim_labels::Union{Vector{String},Nothing}=basename.(simulation_paths),
)::Nothing

    x_plot_params = plotParams(x_quantity)
    y_plot_params = plotParams(y_quantity)

    timeSeriesPlot(
        simulation_paths,
        [lines!];
        pf_kwargs=[(;)],
        # `timeSeriesPlot` configuration
        output_path,
        filename="$(y_quantity)-vs-$(x_quantity)",
        output_format=".png",
        warnings=true,
        show_progress=true,
        # Data manipulation options
        slice=(:),
        da_functions=[daSFRtxt],
        da_args=[(x_quantity, y_quantity)],
        da_kwargs=[(; warnings=true)],
        post_processing=getNothing,
        pp_args=(),
        pp_kwargs=(;),
        x_unit=x_plot_params.unit,
        y_unit=y_plot_params.unit,
        x_exp_factor=0,
        y_exp_factor=0,
        x_trim=(-Inf, Inf),
        y_trim=(-Inf, Inf),
        x_edges=false,
        y_edges=false,
        x_func=identity,
        y_func=identity,
        # Axes options
        xaxis_label=x_plot_params.axis_label,
        yaxis_label=y_plot_params.axis_label,
        xaxis_var_name=x_plot_params.var_name,
        yaxis_var_name=y_plot_params.var_name,
        xaxis_scale_func=identity,
        yaxis_scale_func=log10,
        xaxis_limits=(nothing, nothing),
        yaxis_limits=(nothing, nothing),
        # Plotting options
        save_figure=true,
        backup_results=false,
        sim_labels,
        title="",
        pt_per_unit=0.75,
        px_per_unit=2.0,
        resolution=(1280, 800),
        aspect=nothing,
        series_colors=nothing,
        series_markers=nothing,
        series_linestyles=nothing,
    )

    return nothing

end

"""
    cpuTXT(
        simulation_paths::Vector{String},
        target::String,
        x_quantity::Symbol,
        y_quantity::Symbol;
        <keyword arguments>
    )::Nothing

Plot a time series of the data in the `cpu.txt` file.

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`.
  - `target::String`: Target process.
  - `x_quantity::Symbol`: Quantity for the x axis. The possibilities are:

      + `:time_step`              -> Time step.
      + `:physical_time`          -> Physical time since the Big Bang.
      + `:clock_time_s`           -> Clock time duration of the time step in seconds.
      + `:clock_time_percent`     -> Clock time duration of the time step as a percentage.
      + `:cum_clock_time_s`       -> Cumulative clock time in seconds.
      + `:cum_clock_time_percent` -> Cumulative clock time as a percentage.
  - `y_quantity::Symbol`: Quantity for the y axis. The possibilities are:

      + `:time_step`              -> Time step.
      + `:physical_time`          -> Physical time since the Big Bang.
      + `:clock_time_s`           -> Clock time duration of the time step in seconds.
      + `:clock_time_percent`     -> Clock time duration of the time step as a percentage.
      + `:cum_clock_time_s`       -> Cumulative clock time in seconds.
      + `:cum_clock_time_percent` -> Cumulative clock time as a percentage.
  - `smooth::Int=0`: The result will be smooth out using `smooth` bins. Set it to 0 if you want no smoothing.
  - `yscale::Function=identity`: Scaling function for the y axis. The options are the scaling functions accepted by [Makie](https://docs.makie.org/stable/): log10, log2, log, sqrt, Makie.logit, Makie.Symlog10, Makie.pseudolog10, and identity.
  - `x_trim::NTuple{2,<:Real}=(-Inf, Inf)`: The data will be trim down so the x coordinates fit within `x_trim`.
  - `y_trim::NTuple{2,<:Real}=(-Inf, Inf)`: The data will be trim down so the y coordinates fit within `y_trim`. This option does not affect histograms.
  - `output_path::String="./"`: Path to the output folder.
  - `sim_labels::Union{Vector{String},Nothing}=basename.(simulation_paths)`: Labels for the plot legend, one per simulation. Set it to `nothing` if you don't want a legend.
"""
function cpuTXT(
    simulation_paths::Vector{String},
    target::String,
    x_quantity::Symbol,
    y_quantity::Symbol;
    smooth::Int=0,
    yscale::Function=identity,
    x_trim::NTuple{2,<:Real}=(-Inf, Inf),
    y_trim::NTuple{2,<:Real}=(-Inf, Inf),
    output_path::String="./",
    sim_labels::Union{Vector{String},Nothing}=basename.(simulation_paths),
)::Nothing

    x_plot_params = plotParams(x_quantity)
    y_plot_params = plotParams(y_quantity)

    safe_str_target = replace(target, "/" => "-", "_" => "-")

    timeSeriesPlot(
        simulation_paths,
        [lines!];
        pf_kwargs=[(;)],
        # `timeSeriesPlot` configuration
        output_path,
        filename="$(y_quantity)-vs-$(x_quantity)-for-$(safe_str_target)",
        output_format=".png",
        warnings=true,
        show_progress=true,
        # Data manipulation options
        slice=(:),
        da_functions=[daCPUtxt],
        da_args=[(target, x_quantity, y_quantity)],
        da_kwargs=[(; smooth, warnings=true)],
        post_processing=getNothing,
        pp_args=(),
        pp_kwargs=(;),
        x_unit=x_plot_params.unit,
        y_unit=y_plot_params.unit,
        x_exp_factor=0,
        y_exp_factor=0,
        x_trim,
        y_trim,
        x_edges=false,
        y_edges=false,
        x_func=identity,
        y_func=identity,
        # Axes options
        xaxis_label=x_plot_params.axis_label,
        yaxis_label=y_plot_params.axis_label,
        xaxis_var_name=x_plot_params.var_name,
        yaxis_var_name=y_plot_params.var_name,
        xaxis_scale_func=identity,
        yaxis_scale_func=yscale,
        xaxis_limits=(nothing, nothing),
        yaxis_limits=(nothing, nothing),
        # Plotting options
        save_figure=true,
        backup_results=false,
        sim_labels,
        title=L"\mathrm{Process: \,\, %$(safe_str_target)}",
        pt_per_unit=0.75,
        px_per_unit=2.0,
        resolution=(1280, 800),
        aspect=nothing,
        series_colors=nothing,
        series_markers=nothing,
        series_linestyles=nothing,
    )

    return nothing

end

"""
    densityMap(
        simulation_paths::Vector{String},
        slice_n::Int;
        <keyword arguments>
    )::Nothing

Plot a 2D histogram of the density.

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`.
  - `slice_n::Int`: Selects which snapshot to plot, starts at 1 and is independent of the number in the file name. If set to 0, an animation using every snapshots will be made. If every snapshot is present, `slice_n` = filename_number + 1.
  - `quantities::Vector{Symbol}=[:gas_mass]`: Quantities for which the density will be calculated. The options are:

      + `:stellar_mass`   -> Stellar mass.
      + `:gas_mass`       -> Gas mass.
      + `:dm_mass`        -> Dark matter mass.
      + `:bh_mass`        -> Black hole mass.
      + `:molecular_mass` -> Molecular hydrogen (``\\mathrm{H_2}``) mass.
      + `:atomic_mass`    -> Atomic hydrogen (``\\mathrm{HI}``) mass.
      + `:ionized_mass`   -> Ionized hydrogen (``\\mathrm{HII}``) mass.
      + `:neutral_mass`   -> Neutral hydrogen (``\\mathrm{HI + H_2}``) mass.
  - `output_path::String="./"`: Path to the output folder.
  - `filter_mode::Symbol=:all`: Which cells/particles will be plotted, the options are:

      + `:all`             -> Plot every cell/particle within the simulation box.
      + `:halo`            -> Plot only the cells/particles that belong to the main halo.
      + `:subhalo`         -> Plot only the cells/particles that belong to the main subhalo.
      + `:sphere`          -> Plot only the cell/particle inside a sphere with radius `FILTER_R` (see `./src/constants.jl`).
      + `:stellar_subhalo` -> Plot only the cells/particles that belong to the main subhalo.
  - `projection_planes::Vector{Symbol}=[:xy]`: Projection planes. The options are `:xy`, `:xz` and `:yz`.
  - `box_size::Unitful.Length=100u"kpc"`: Physical side length of the plot window.
  - `pixel_length::Unitful.Length=0.1u"kpc"`: Pixel (bin of the 2D histogram) side length.
  - `smooth::Bool=false`: If the results will be smooth out using the kernel function [`cubicSplineKernel`](@ref).
  - `colorrange::Union{Nothing,Tuple{<:Real,<:Real}}=nothing`: Sets the start and end points of the colormap. Use `nothing` to use the extrema of the values to be plotted.
"""
function densityMap(
    simulation_paths::Vector{String},
    slice_n::Int;
    quantities::Vector{Symbol}=[:gas_mass],
    output_path::String="./",
    filter_mode::Symbol=:all,
    projection_planes::Vector{Symbol}=[:xy],
    box_size::Unitful.Length=100u"kpc",
    pixel_length::Unitful.Length=0.1u"kpc",
    smooth::Bool=false,
    colorrange::Union{Nothing,Tuple{<:Real,<:Real}}=nothing,
)::Nothing

    # Compute number of pixel per side
    resolution = round(Int, box_size / pixel_length)

    # Set up the grid
    grid = SquareGrid(box_size, resolution)

    pf_kwargs = isnothing(colorrange) ? [(;)] : [(; colorrange)]

    @inbounds for quantity in quantities

        filter_function, translation, rotation, request = selectFilter(
            filter_mode,
            plotParams(quantity).request,
        )

        @inbounds for simulation_path in simulation_paths

            # Get the simulation name as a string
            sim_name = basename(simulation_path)

            @inbounds for projection_plane in projection_planes

                # Construct the file name
                base_filename = "$(sim_name)-$(quantity)-$(projection_plane)-density_map"

                snapshotPlot(
                    [simulation_path],
                    request,
                    [heatmap!];
                    pf_kwargs,
                    # `snapshotPlot` configuration
                    output_path,
                    base_filename,
                    output_format=".png",
                    warnings=true,
                    show_progress=iszero(slice_n),
                    # Data manipulation options
                    slice=iszero(slice_n) ? (:) : slice_n,
                    filter_function,
                    da_functions=[daDensity2DHistogram],
                    da_args=[(grid, quantity)],
                    da_kwargs=[(; projection_plane, smooth, neighbors=32)],
                    post_processing=getNothing,
                    pp_args=(),
                    pp_kwargs=(;),
                    transform_box=true,
                    translation,
                    rotation,
                    smooth=0,
                    x_unit=u"kpc",
                    y_unit=u"kpc",
                    x_exp_factor=0,
                    y_exp_factor=0,
                    x_trim=(-Inf, Inf),
                    y_trim=(-Inf, Inf),
                    x_edges=false,
                    y_edges=false,
                    x_func=identity,
                    y_func=identity,
                    # Axes options
                    xaxis_label="auto_label",
                    yaxis_label="auto_label",
                    xaxis_var_name=string(projection_plane)[1:1],
                    yaxis_var_name=string(projection_plane)[2:2],
                    xaxis_scale_func=identity,
                    yaxis_scale_func=identity,
                    xaxis_limits=(nothing, nothing),
                    yaxis_limits=(nothing, nothing),
                    # Plotting options
                    save_figures=!iszero(slice_n),
                    backup_results=false,
                    sim_labels=nothing,
                    title=:physical_time,
                    pt_per_unit=0.75,
                    px_per_unit=2.0,
                    resolution=(1000, 1000),
                    aspect=AxisAspect(1),
                    series_colors=nothing,
                    series_markers=nothing,
                    series_linestyles=nothing,
                    # Animation options
                    animation=iszero(slice_n),
                    animation_filename="$(base_filename).mp4",
                    framerate=5,
                )

            end

        end

    end

    return nothing

end

"""
    scatterPlot(
        simulation_paths::Vector{String},
        slice_n::Int,
        x_quantity::Symbol,
        y_quantity::Symbol;
        <keyword arguments>
    )::Nothing

Plot two quantities as a scatter plot, one marker for every cell/particle.

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`.
  - `slice_n::Int`: Selects which snapshot to plot, starts at 1 and is independent of the number in the file name. If every snapshot is present, `slice_n` = filename_number + 1.
  - `x_quantity::Symbol`: Quantity for the x axis. The possibilities are:

      + `:stellar_mass`             -> Stellar mass.
      + `:gas_mass`                 -> Gas mass.
      + `:dm_mass`                  -> Dark matter mass.
      + `:bh_mass`                  -> Black hole mass.
      + `:molecular_mass`           -> Molecular hydrogen (``\\mathrm{H_2}``) mass.
      + `:atomic_mass`              -> Atomic hydrogen (``\\mathrm{HI}``) mass.
      + `:ionized_mass`             -> Ionized hydrogen (``\\mathrm{HII}``) mass.
      + `:neutral_mass`             -> Neutral hydrogen (``\\mathrm{HI + H_2}``) mass.
      + `:molecular_fraction`       -> Gas mass fraction of molecular hydrogen.
      + `:atomic_fraction`          -> Gas mass fraction of atomic hydrogen.
      + `:ionized_fraction`         -> Gas mass fraction of ionized hydrogen.
      + `:neutral_fraction`         -> Gas mass fraction of neutral hydrogen.
      + `:gas_mass_density`         -> Gas mass density.
      + `:gas_number_density`       -> Gas number density.
      + `:molecular_number_density` -> Molecular hydrogen number density.
      + `:atomic_number_density`    -> Atomic hydrogen number density.
      + `:ionized_number_density`   -> Ionized hydrogen number density.
      + `:neutral_number_density`   -> Neutral hydrogen number density.
      + `:gas_metallicity`          -> Mass fraction of all elements above He in the gas (solar units).
      + `:stellar_metallicity`      -> Mass fraction of all elements above He in the stars (solar units).
      + `:X_gas_abundance`          -> Gas abundance of element ``\\mathrm{X}``, as ``12 + \\log_{10}(\\mathrm{X \\, / \\, H})``. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
      + `:X_stellar_abundance`      -> Stellar abundance of element ``\\mathrm{X}``, as ``12 + \\log_{10}(\\mathrm{X \\, / \\, H})``. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
      + `:stellar_radial_distance`  -> Distance of every stellar particle to the origin.
      + `:gas_radial_distance`      -> Distance of every gas cell to the origin.
      + `:dm_radial_distance`       -> Distance of every dark matter particle to the origin.
      + `:stellar_xy_distance`      -> Projected distance of every stellar particle to the origin.
      + `:gas_xy_distance`          -> Projected distance of every gas cell to the origin.
      + `:dm_xy_distance`           -> Projected distance of every dark matter particle to the origin.
      + `:stellar_circularity`      -> Stellar circularity.
      + `:stellar_vcirc`            -> Stellar circular velocity.
      + `:stellar_age`              -> Stellar age.
      + `:sfr`                      -> The star formation rate of the last `AGE_RESOLUTION`.
      + `:ssfr`                     -> The specific star formation rate of the last `AGE_RESOLUTION`.
      + `:temperature`              -> Gas temperature, as ``\\log_{10}(T \\, / \\, \\mathrm{K})``.
  - `y_quantity::Symbol`: Quantity for the y axis. The possibilities are:

      + `:stellar_mass`             -> Stellar mass.
      + `:gas_mass`                 -> Gas mass.
      + `:dm_mass`                  -> Dark matter mass.
      + `:bh_mass`                  -> Black hole mass.
      + `:molecular_mass`           -> Molecular hydrogen (``\\mathrm{H_2}``) mass.
      + `:atomic_mass`              -> Atomic hydrogen (``\\mathrm{HI}``) mass.
      + `:ionized_mass`             -> Ionized hydrogen (``\\mathrm{HII}``) mass.
      + `:neutral_mass`             -> Neutral hydrogen (``\\mathrm{HI + H_2}``) mass.
      + `:molecular_fraction`       -> Gas mass fraction of molecular hydrogen.
      + `:atomic_fraction`          -> Gas mass fraction of atomic hydrogen.
      + `:ionized_fraction`         -> Gas mass fraction of ionized hydrogen.
      + `:neutral_fraction`         -> Gas mass fraction of neutral hydrogen.
      + `:gas_mass_density`         -> Gas mass density.
      + `:gas_number_density`       -> Gas number density.
      + `:molecular_number_density` -> Molecular hydrogen number density.
      + `:atomic_number_density`    -> Atomic hydrogen number density.
      + `:ionized_number_density`   -> Ionized hydrogen number density.
      + `:neutral_number_density`   -> Neutral hydrogen number density.
      + `:gas_metallicity`          -> Mass fraction of all elements above He in the gas (solar units).
      + `:stellar_metallicity`      -> Mass fraction of all elements above He in the stars (solar units).
      + `:X_gas_abundance`          -> Gas abundance of element ``\\mathrm{X}``, as ``12 + \\log_{10}(\\mathrm{X \\, / \\, H})``. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
      + `:X_stellar_abundance`      -> Stellar abundance of element ``\\mathrm{X}``, as ``12 + \\log_{10}(\\mathrm{X \\, / \\, H})``. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
      + `:stellar_radial_distance`  -> Distance of every stellar particle to the origin.
      + `:gas_radial_distance`      -> Distance of every gas cell to the origin.
      + `:dm_radial_distance`       -> Distance of every dark matter particle to the origin.
      + `:stellar_xy_distance`      -> Projected distance of every stellar particle to the origin.
      + `:gas_xy_distance`          -> Projected distance of every gas cell to the origin.
      + `:dm_xy_distance`           -> Projected distance of every dark matter particle to the origin.
      + `:stellar_circularity`      -> Stellar circularity.
      + `:stellar_vcirc`            -> Stellar circular velocity.
      + `:stellar_age`              -> Stellar age.
      + `:sfr`                      -> The star formation rate of the last `AGE_RESOLUTION`.
      + `:ssfr`                     -> The specific star formation rate of the last `AGE_RESOLUTION`.
      + `:temperature`              -> Gas temperature, as ``\\log_{10}(T \\, / \\, \\mathrm{K})``.
  - `output_path::String="./"`: Path to the output folder.
  - `filter_mode::Symbol=:all`: Which cells/particles will be plotted, the options are:

      + `:all`             -> Plot every cell/particle within the simulation box.
      + `:halo`            -> Plot only the cells/particles that belong to the main halo.
      + `:subhalo`         -> Plot only the cells/particles that belong to the main subhalo.
      + `:sphere`          -> Plot only the cell/particle inside a sphere with radius `FILTER_R` (see `./src/constants.jl`).
      + `:stellar_subhalo` -> Plot only the cells/particles that belong to the main subhalo.
"""
function scatterPlot(
    simulation_paths::Vector{String},
    slice_n::Int,
    x_quantity::Symbol,
    y_quantity::Symbol;
    output_path::String="./",
    filter_mode::Symbol=:all,
)::Nothing

    x_plot_params = plotParams(x_quantity)
    y_plot_params = plotParams(y_quantity)

    filter_function, translation, rotation, request = selectFilter(
        filter_mode,
        mergeRequests(x_plot_params.request, y_plot_params.request),
    )

    @inbounds for simulation_path in simulation_paths

        # Get the simulation name as a string
        sim_name = basename(simulation_path)

        snapshotPlot(
            [simulation_path],
            request,
            [scatter!];
            pf_kwargs=[(; markersize=1)],
            # `snapshotPlot` configuration
            output_path,
            base_filename="$(sim_name)-$(y_quantity)-vs-$(x_quantity)",
            output_format=".png",
            warnings=true,
            show_progress=false,
            # Data manipulation options
            slice=slice_n,
            filter_function,
            da_functions=[daScatterGalaxy],
            da_args=[(x_quantity, y_quantity)],
            da_kwargs=[(;)],
            post_processing=getNothing,
            pp_args=(),
            pp_kwargs=(;),
            transform_box=true,
            translation,
            rotation,
            smooth=0,
            x_unit=x_plot_params.unit,
            y_unit=y_plot_params.unit,
            x_exp_factor=x_plot_params.exp_factor,
            y_exp_factor=y_plot_params.exp_factor,
            x_trim=(-Inf, Inf),
            y_trim=(-Inf, Inf),
            x_edges=false,
            y_edges=false,
            x_func=identity,
            y_func=identity,
            # Axes options
            xaxis_label=x_plot_params.axis_label,
            yaxis_label=y_plot_params.axis_label,
            xaxis_var_name=x_plot_params.var_name,
            yaxis_var_name=y_plot_params.var_name,
            xaxis_scale_func=identity,
            yaxis_scale_func=identity,
            xaxis_limits=(nothing, nothing),
            yaxis_limits=(nothing, nothing),
            # Plotting and animation options
            save_figures=true,
            backup_results=false,
            sim_labels=nothing,
            title=:physical_time,
            pt_per_unit=0.75,
            px_per_unit=2.0,
            resolution=(1280, 800),
            aspect=nothing,
            series_colors=nothing,
            series_markers=nothing,
            series_linestyles=nothing,
            # Animation options
            animation=false,
            animation_filename="animation.mp4",
            framerate=10,
        )

    end

    return nothing

end

"""
    scatterDensityMap(
        simulation_paths::Vector{String},
        slice_n::Int,
        x_quantity::Symbol,
        y_quantity::Symbol,
        x_range::NTuple{2, <:Number},
        y_range::NTuple{2, <:Number},
        n_bins::Int;
        <keyword arguments>
    )::Nothing

Plot two quantities as a density scatter plot (2D histogram).

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`.
  - `slice_n::Int`: Selects which snapshot to plot, starts at 1 and is independent of the number in the file name. If set to 0, an animation using every snapshots will be made. If every snapshot is present, `slice_n` = filename_number + 1.
  - `x_quantity::Symbol`: Quantity for the x axis. The possibilities are:

      + `:stellar_mass`             -> Stellar mass.
      + `:gas_mass`                 -> Gas mass.
      + `:dm_mass`                  -> Dark matter mass.
      + `:bh_mass`                  -> Black hole mass.
      + `:molecular_mass`           -> Molecular hydrogen (``\\mathrm{H_2}``) mass.
      + `:atomic_mass`              -> Atomic hydrogen (``\\mathrm{HI}``) mass.
      + `:ionized_mass`             -> Ionized hydrogen (``\\mathrm{HII}``) mass.
      + `:neutral_mass`             -> Neutral hydrogen (``\\mathrm{HI + H_2}``) mass.
      + `:molecular_fraction`       -> Gas mass fraction of molecular hydrogen.
      + `:atomic_fraction`          -> Gas mass fraction of atomic hydrogen.
      + `:ionized_fraction`         -> Gas mass fraction of ionized hydrogen.
      + `:neutral_fraction`         -> Gas mass fraction of neutral hydrogen.
      + `:gas_mass_density`         -> Gas mass density.
      + `:gas_number_density`       -> Gas number density.
      + `:molecular_number_density` -> Molecular hydrogen number density.
      + `:atomic_number_density`    -> Atomic hydrogen number density.
      + `:ionized_number_density`   -> Ionized hydrogen number density.
      + `:neutral_number_density`   -> Neutral hydrogen number density.
      + `:gas_metallicity`          -> Mass fraction of all elements above He in the gas (solar units).
      + `:stellar_metallicity`      -> Mass fraction of all elements above He in the stars (solar units).
      + `:X_gas_abundance`          -> Gas abundance of element ``\\mathrm{X}``, as ``12 + \\log_{10}(\\mathrm{X \\, / \\, H})``. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
      + `:X_stellar_abundance`      -> Stellar abundance of element ``\\mathrm{X}``, as ``12 + \\log_{10}(\\mathrm{X \\, / \\, H})``. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
      + `:stellar_radial_distance`  -> Distance of every stellar particle to the origin.
      + `:gas_radial_distance`      -> Distance of every gas cell to the origin.
      + `:dm_radial_distance`       -> Distance of every dark matter particle to the origin.
      + `:stellar_xy_distance`      -> Projected distance of every stellar particle to the origin.
      + `:gas_xy_distance`          -> Projected distance of every gas cell to the origin.
      + `:dm_xy_distance`           -> Projected distance of every dark matter particle to the origin.
      + `:stellar_circularity`      -> Stellar circularity.
      + `:stellar_vcirc`            -> Stellar circular velocity.
      + `:stellar_age`              -> Stellar age.
      + `:sfr`                      -> The star formation rate of the last `AGE_RESOLUTION`.
      + `:ssfr`                     -> The specific star formation rate of the last `AGE_RESOLUTION`.
      + `:temperature`              -> Gas temperature, as ``\\log_{10}(T \\, / \\, \\mathrm{K})``.
  - `y_quantity::Symbol`: Quantity for the y axis. The possibilities are:

      + `:stellar_mass`             -> Stellar mass.
      + `:gas_mass`                 -> Gas mass.
      + `:dm_mass`                  -> Dark matter mass.
      + `:bh_mass`                  -> Black hole mass.
      + `:molecular_mass`           -> Molecular hydrogen (``\\mathrm{H_2}``) mass.
      + `:atomic_mass`              -> Atomic hydrogen (``\\mathrm{HI}``) mass.
      + `:ionized_mass`             -> Ionized hydrogen (``\\mathrm{HII}``) mass.
      + `:neutral_mass`             -> Neutral hydrogen (``\\mathrm{HI + H_2}``) mass.
      + `:molecular_fraction`       -> Gas mass fraction of molecular hydrogen.
      + `:atomic_fraction`          -> Gas mass fraction of atomic hydrogen.
      + `:ionized_fraction`         -> Gas mass fraction of ionized hydrogen.
      + `:neutral_fraction`         -> Gas mass fraction of neutral hydrogen.
      + `:gas_mass_density`         -> Gas mass density.
      + `:gas_number_density`       -> Gas number density.
      + `:molecular_number_density` -> Molecular hydrogen number density.
      + `:atomic_number_density`    -> Atomic hydrogen number density.
      + `:ionized_number_density`   -> Ionized hydrogen number density.
      + `:neutral_number_density`   -> Neutral hydrogen number density.
      + `:gas_metallicity`          -> Mass fraction of all elements above He in the gas (solar units).
      + `:stellar_metallicity`      -> Mass fraction of all elements above He in the stars (solar units).
      + `:X_gas_abundance`          -> Gas abundance of element ``\\mathrm{X}``, as ``12 + \\log_{10}(\\mathrm{X \\, / \\, H})``. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
      + `:X_stellar_abundance`      -> Stellar abundance of element ``\\mathrm{X}``, as ``12 + \\log_{10}(\\mathrm{X \\, / \\, H})``. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
      + `:stellar_radial_distance`  -> Distance of every stellar particle to the origin.
      + `:gas_radial_distance`      -> Distance of every gas cell to the origin.
      + `:dm_radial_distance`       -> Distance of every dark matter particle to the origin.
      + `:stellar_xy_distance`      -> Projected distance of every stellar particle to the origin.
      + `:gas_xy_distance`          -> Projected distance of every gas cell to the origin.
      + `:dm_xy_distance`           -> Projected distance of every dark matter particle to the origin.
      + `:stellar_circularity`      -> Stellar circularity.
      + `:stellar_vcirc`            -> Stellar circular velocity.
      + `:stellar_age`              -> Stellar age.
      + `:sfr`                      -> The star formation rate of the last `AGE_RESOLUTION`.
      + `:ssfr`                     -> The specific star formation rate of the last `AGE_RESOLUTION`.
      + `:temperature`              -> Gas temperature, as ``\\log_{10}(T \\, / \\, \\mathrm{K})``.
  - `x_range::NTuple{2,<:Number}`: x axis range.
  - `y_range::NTuple{2,<:Number}`: y axis range.
  - `n_bins::Int`: Number of bins per side of the square grid.
  - `output_path::String="./"`: Path to the output folder.
  - `filter_mode::Symbol=:all`: Which cells/particles will be plotted, the options are:

      + `:all`             -> Plot every cell/particle within the simulation box.
      + `:halo`            -> Plot only the cells/particles that belong to the main halo.
      + `:subhalo`         -> Plot only the cells/particles that belong to the main subhalo.
      + `:sphere`          -> Plot only the cell/particle inside a sphere with radius `FILTER_R` (see `./src/constants.jl`).
      + `:stellar_subhalo` -> Plot only the cells/particles that belong to the main subhalo
"""
function scatterDensityMap(
    simulation_paths::Vector{String},
    slice_n::Int,
    x_quantity::Symbol,
    y_quantity::Symbol,
    x_range::NTuple{2,<:Number},
    y_range::NTuple{2,<:Number},
    n_bins::Int;
    output_path::String="./",
    filter_mode::Symbol=:all,
)::Nothing

    x_plot_params = plotParams(x_quantity)
    y_plot_params = plotParams(y_quantity)

    filter_function, translation, rotation, request = selectFilter(
        filter_mode,
        mergeRequests(x_plot_params.request, y_plot_params.request),
    )

    @inbounds for simulation_path in simulation_paths

        # Get the simulation name as a string
        sim_name = basename(simulation_path)

        snapshotPlot(
            [simulation_path],
            request,
            [heatmap!];
            pf_kwargs=[(;)],
            # `snapshotPlot` configuration
            output_path,
            base_filename="$(sim_name)-$(y_quantity)-vs-$(x_quantity)",
            output_format=".png",
            warnings=true,
            show_progress=false,
            # Data manipulation options
            slice=slice_n,
            filter_function,
            da_functions=[daScatterDensity],
            da_args=[(x_quantity, y_quantity, x_range, y_range, n_bins)],
            da_kwargs=[(;)],
            post_processing=getNothing,
            pp_args=(),
            pp_kwargs=(;),
            transform_box=true,
            translation,
            rotation,
            smooth=0,
            x_unit=x_plot_params.unit,
            y_unit=y_plot_params.unit,
            x_exp_factor=x_plot_params.exp_factor,
            y_exp_factor=y_plot_params.exp_factor,
            x_trim=(-Inf, Inf),
            y_trim=(-Inf, Inf),
            x_edges=false,
            y_edges=false,
            x_func=identity,
            y_func=identity,
            # Axes options
            xaxis_label=x_plot_params.axis_label,
            yaxis_label=y_plot_params.axis_label,
            xaxis_var_name=x_plot_params.var_name,
            yaxis_var_name=y_plot_params.var_name,
            xaxis_scale_func=identity,
            yaxis_scale_func=identity,
            xaxis_limits=(nothing, nothing),
            yaxis_limits=(nothing, nothing),
            # Plotting and animation options
            save_figures=true,
            backup_results=false,
            sim_labels=nothing,
            title=:physical_time,
            pt_per_unit=0.75,
            px_per_unit=2.0,
            resolution=(1280, 800),
            aspect=nothing,
            series_colors=nothing,
            series_markers=nothing,
            series_linestyles=nothing,
            # Animation options
            animation=false,
            animation_filename="animation.mp4",
            framerate=10,
        )

    end

    return nothing

end

"""
    timeSeries(
        simulation_paths::Vector{String},
        x_quantity::Symbol,
        y_quantity::Symbol;
        <keyword arguments>
    )::Nothing

Plot a time series.

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`.
  - `x_quantity::Symbol`: Quantity for the x axis. The possibilities are:

      + `:stellar_mass`           -> Stellar mass.
      + `:gas_mass`               -> Gas mass.
      + `:dm_mass`                -> Dark matter mass.
      + `:bh_mass`                -> Black hole mass.
      + `:stellar_number`         -> Number of stellar particles.
      + `:gas_number`             -> Number of gas cells.
      + `:dm_number`              -> Number of dark matter particles.
      + `:bh_number`              -> Number of black hole particles.
      + `:molecular_mass`         -> Molecular hydrogen (``\\mathrm{H_2}``) mass.
      + `:atomic_mass`            -> Atomic hydrogen (``\\mathrm{HI}``) mass.
      + `:ionized_mass`           -> Ionized hydrogen (``\\mathrm{HII}``) mass.
      + `:neutral_mass`           -> Neutral hydrogen (``\\mathrm{HI + H_2}``) mass.
      + `:molecular_fraction`     -> Gas mass fraction of molecular hydrogen.
      + `:atomic_fraction`        -> Gas mass fraction of atomic hydrogen.
      + `:ionized_fraction`       -> Gas mass fraction of ionized hydrogen.
      + `:neutral_fraction`       -> Gas mass fraction of neutral hydrogen.
      + `:stellar_area_density`   -> Stellar area mass density, for a radius of `FILTER_R`.
      + `:gas_area_density`       -> Gas area mass density, for a radius of `FILTER_R`.
      + `:molecular_area_density` -> Molecular hydrogen area mass density, for a radius of `FILTER_R`.
      + `:atomic_area_density`    -> Atomic hydrogen area mass density, for a radius of `FILTER_R`.
      + `:ionized_area_density`   -> Ionized hydrogen area mass density, for a radius of `FILTER_R`.
      + `:neutral_area_density`   -> Neutral hydrogen area mass density, for a radius of `FILTER_R`.
      + `:sfr_area_density`       -> Star formation rate area density, for the last `AGE_RESOLUTION_ρ` and a radius of `FILTER_R`.
      + `:gas_metallicity`        -> Mass fraction of all elements above He in the gas (solar units).
      + `:stellar_metallicity`    -> Mass fraction of all elements above He in the stars (solar units).
      + `:X_gas_abundance`        -> Gas abundance of element ``\\mathrm{X}``, as ``12 + \\log_{10}(\\mathrm{X \\, / \\, H})``. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
      + `:X_stellar_abundance`    -> Stellar abundance of element ``\\mathrm{X}``, as ``12 + \\log_{10}(\\mathrm{X \\, / \\, H})``. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
      + `:stellar_specific_am`    -> Norm of the stellar specific angular momentum.
      + `:gas_specific_am`        -> Norm of the gas specific angular momentum.
      + `:dm_specific_am`         -> Norm of the dark matter specific angular momentum.
      + `:sfr`                    -> The star formation rate of the last `AGE_RESOLUTION`.
      + `:ssfr`                   -> The specific star formation rate of the last `AGE_RESOLUTION`.
      + `:scale_factor`           -> Scale factor.
      + `:redshift`               -> Redshift.
      + `:physical_time`          -> Physical time since the Big Bang.
      + `:lookback_time`          -> Physical time left to reach the last snapshot.
  - `y_quantity::Symbol`: Quantity for the y axis. The possibilities are:

      + `:stellar_mass`           -> Stellar mass.
      + `:gas_mass`               -> Gas mass.
      + `:dm_mass`                -> Dark matter mass.
      + `:bh_mass`                -> Black hole mass.
      + `:stellar_number`         -> Number of stellar particles.
      + `:gas_number`             -> Number of gas cells.
      + `:dm_number`              -> Number of dark matter particles.
      + `:bh_number`              -> Number of black hole particles.
      + `:molecular_mass`         -> Molecular hydrogen (``\\mathrm{H_2}``) mass.
      + `:atomic_mass`            -> Atomic hydrogen (``\\mathrm{HI}``) mass.
      + `:ionized_mass`           -> Ionized hydrogen (``\\mathrm{HII}``) mass.
      + `:neutral_mass`           -> Neutral hydrogen (``\\mathrm{HI + H_2}``) mass.
      + `:molecular_fraction`     -> Gas mass fraction of molecular hydrogen.
      + `:atomic_fraction`        -> Gas mass fraction of atomic hydrogen.
      + `:ionized_fraction`       -> Gas mass fraction of ionized hydrogen.
      + `:neutral_fraction`       -> Gas mass fraction of neutral hydrogen.
      + `:stellar_area_density`   -> Stellar area mass density, for a radius of `FILTER_R`.
      + `:gas_area_density`       -> Gas area mass density, for a radius of `FILTER_R`.
      + `:molecular_area_density` -> Molecular hydrogen area mass density, for a radius of `FILTER_R`.
      + `:atomic_area_density`    -> Atomic hydrogen area mass density, for a radius of `FILTER_R`.
      + `:ionized_area_density`   -> Ionized hydrogen area mass density, for a radius of `FILTER_R`.
      + `:neutral_area_density`   -> Neutral hydrogen area mass density, for a radius of `FILTER_R`.
      + `:sfr_area_density`       -> Star formation rate area density, for the last `AGE_RESOLUTION_ρ` and a radius of `FILTER_R`.
      + `:gas_metallicity`        -> Mass fraction of all elements above He in the gas (solar units).
      + `:stellar_metallicity`    -> Mass fraction of all elements above He in the stars (solar units).
      + `:X_gas_abundance`        -> Gas abundance of element ``\\mathrm{X}``, as ``12 + \\log_{10}(\\mathrm{X \\, / \\, H})``. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
      + `:X_stellar_abundance`    -> Stellar abundance of element ``\\mathrm{X}``, as ``12 + \\log_{10}(\\mathrm{X \\, / \\, H})``. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
      + `:stellar_specific_am`    -> Norm of the stellar specific angular momentum.
      + `:gas_specific_am`        -> Norm of the gas specific angular momentum.
      + `:dm_specific_am`         -> Norm of the dark matter specific angular momentum.
      + `:sfr`                    -> The star formation rate of the last `AGE_RESOLUTION`.
      + `:ssfr`                   -> The specific star formation rate of the last `AGE_RESOLUTION`.
      + `:scale_factor`           -> Scale factor.
      + `:redshift`               -> Redshift.
      + `:physical_time`          -> Physical time since the Big Bang.
      + `:lookback_time`          -> Physical time left to reach the last snapshot.
  - `slice::IndexType=(:)`: Slice of the simulation, i.e. which snapshots will be read. It can be an integer (a single snapshot), a vector of integers (several snapshots), an `UnitRange` (e.g. 5:13), an `StepRange` (e.g. 5:2:13) or (:) (all snapshots). It works over the longest possible list of snapshots among the simulations (grouped by the number in the file names). Out of bounds indices are ignored.
  - `filter_mode::Symbol=:all`: Which cells/particles will be plotted, the options are:

      + `:all`             -> Plot every cell/particle within the simulation box.
      + `:halo`            -> Plot only the cells/particles that belong to the main halo.
      + `:subhalo`         -> Plot only the cells/particles that belong to the main subhalo.
      + `:sphere`          -> Plot only the cell/particle inside a sphere with radius `FILTER_R` (see `./src/constants.jl`).
      + `:stellar_subhalo` -> Plot only the cells/particles that belong to the main subhalo.
  - `sim_labels::Union{Vector{String},Nothing}=nothing`: Labels for the plot legend, one per simulation. Set it to `nothing` if you don't want a legend.
"""
function timeSeries(
    simulation_paths::Vector{String},
    x_quantity::Symbol,
    y_quantity::Symbol;
    slice::IndexType=(:),
    output_path::String="./",
    filter_mode::Symbol=:all,
    sim_labels::Union{Vector{String},Nothing}=basename.(simulation_paths),
)::Nothing

    x_plot_params = plotParams(x_quantity)
    y_plot_params = plotParams(y_quantity)

    timeSeriesPlot(
        simulation_paths,
        [scatterlines!];
        pf_kwargs=[(;)],
        # `timeSeriesPlot` configuration
        output_path,
        filename="$(y_quantity)-vs-$(x_quantity)",
        output_format=".png",
        warnings=true,
        show_progress=true,
        # Data manipulation options
        slice,
        da_functions=[daEvolution],
        da_args=[(x_quantity, y_quantity)],
        da_kwargs=[(; filter_mode, smooth=0, scaling=identity, warnings=true)],
        post_processing=getNothing,
        pp_args=(),
        pp_kwargs=(;),
        x_unit=x_plot_params.unit,
        y_unit=y_plot_params.unit,
        x_exp_factor=x_plot_params.exp_factor,
        y_exp_factor=y_plot_params.exp_factor,
        x_trim=(-Inf, Inf),
        y_trim=(-Inf, Inf),
        x_edges=false,
        y_edges=false,
        x_func=identity,
        y_func=identity,
        # Axes options
        xaxis_label=x_plot_params.axis_label,
        yaxis_label=y_plot_params.axis_label,
        xaxis_var_name=x_plot_params.var_name,
        yaxis_var_name=y_plot_params.var_name,
        xaxis_scale_func=identity,
        yaxis_scale_func=identity,
        xaxis_limits=(nothing, nothing),
        yaxis_limits=(nothing, nothing),
        # Plotting options
        save_figure=true,
        backup_results=false,
        sim_labels,
        title="",
        pt_per_unit=0.75,
        px_per_unit=2.0,
        resolution=(1280, 800),
        aspect=nothing,
        series_colors=nothing,
        series_markers=nothing,
        series_linestyles=nothing,
    )

    return nothing

end

"""
    rotationCurve(
        simulation_paths::Vector{String},
        slice_n::Int;
        <keyword arguments>
    )::Nothing

Plot the galaxy rotation curve of a set of simulations.

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`.
  - `slice_n::Int`: Selects which snapshot to plot, starts at 1 and is independent of the number in the file name. If every snapshot is present, `slice_n` = filename_number + 1.
  - `radius::Unitful.Length=FILTER_R`: Maximum radial distance for the rotation curve.
  - `output_path::String="./"`: Path to the output folder.
  - `filter_mode::Symbol=:all`: Which cells/particles will be plotted, the options are:

      + `:all`             -> Plot every cell/particle within the simulation box.
      + `:halo`            -> Plot only the cells/particles that belong to the main halo.
      + `:subhalo`         -> Plot only the cells/particles that belong to the main subhalo.
      + `:sphere`          -> Plot only the cell/particle inside a sphere with radius `FILTER_R` (see `./src/constants.jl`).
      + `:stellar_subhalo` -> Plot only the cells/particles that belong to the main subhalo.
  - `sim_labels::Union{Vector{String},Nothing}=basename.(simulation_paths)`: Labels for the plot legend, one per simulation. Set it to `nothing` if you don't want a legend.
"""
function rotationCurve(
    simulation_paths::Vector{String},
    slice_n::Int;
    radius::Unitful.Length=FILTER_R,
    output_path::String="./",
    filter_mode::Symbol=:all,
    sim_labels::Union{Vector{String},Nothing}=basename.(simulation_paths),
)::Nothing

    x_plot_params = plotParams(:stellar_radial_distance)
    y_plot_params = plotParams(:stellar_vcirc)

    filter_function, translation, rotation, request = selectFilter(
        filter_mode,
        mergeRequests(x_plot_params.request, y_plot_params.request, Dict(:stars => ["GAGE"])),
    )

    snapshotPlot(
        simulation_paths,
        request,
        [lines!];
        pf_kwargs=[(;)],
        # `snapshotPlot` configuration
        output_path,
        base_filename="rotation_curve",
        output_format=".png",
        warnings=true,
        show_progress=false,
        # Data manipulation options
        slice=slice_n,
        filter_function,
        da_functions=[daRotationCurve],
        da_args=[(radius,)],
        da_kwargs=[(;)],
        post_processing=getNothing,
        pp_args=(),
        pp_kwargs=(;),
        transform_box=true,
        translation,
        rotation,
        smooth=0,
        x_unit=x_plot_params.unit,
        y_unit=y_plot_params.unit,
        x_exp_factor=x_plot_params.exp_factor,
        y_exp_factor=y_plot_params.exp_factor,
        x_trim=(-Inf, Inf),
        y_trim=(-Inf, Inf),
        x_edges=false,
        y_edges=false,
        x_func=identity,
        y_func=identity,
        # Axes options
        xaxis_label=x_plot_params.axis_label,
        yaxis_label=y_plot_params.axis_label,
        xaxis_var_name=x_plot_params.var_name,
        yaxis_var_name=y_plot_params.var_name,
        xaxis_scale_func=identity,
        yaxis_scale_func=identity,
        xaxis_limits=(nothing, nothing),
        yaxis_limits=(nothing, nothing),
        # Plotting and animation options
        save_figures=true,
        backup_results=false,
        sim_labels,
        title=:physical_time,
        pt_per_unit=0.75,
        px_per_unit=2.0,
        resolution=(1280, 800),
        aspect=nothing,
        series_colors=nothing,
        series_markers=nothing,
        series_linestyles=nothing,
        # Animation options
        animation=false,
        animation_filename="animation.mp4",
        framerate=10,
    )

    return nothing

end

"""
    densityProfile(
        simulation_paths::Vector{String},
        slice_n::Int,
        quantity::Symbol;
        <keyword arguments>
    )::Nothing

 Plot a density profile.

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`.
  - `slice_n::Int`: Selects which snapshot to plot, starts at 1 and is independent of the number in the file name. If every snapshot is present, `slice_n` = filename_number + 1.
  - `quantity::Symbol`: Quantity for the y axis. The options are:

      + `:stellar_area_density`     -> Stellar area mass density, up to a radius of `FILTER_R`.
      + `:gas_area_density`         -> Gas area mass density, up to a radius of `FILTER_R`.
      + `:molecular_area_density`   -> Molecular hydrogen area mass density, up to a radius of `FILTER_R`.
      + `:atomic_area_density`      -> Atomic hydrogen area mass density, up to a radius of `FILTER_R`.
      + `:ionized_area_density`     -> Ionized hydrogen area mass density, up to a radius of `FILTER_R`.
      + `:neutral_area_density`     -> Neutral hydrogen area mass density, up to a radius of `FILTER_R`.
      + `:sfr_area_density`         -> Star formation rate area density, up to the last `AGE_RESOLUTION_ρ` and a radius of `FILTER_R`.
  - `output_path::String="./"`: Path to the output folder.
  - `filter_mode::Symbol=:all`: Which cells/particles will be plotted, the options are:

      + `:all`             -> Plot every cell/particle within the simulation box.
      + `:halo`            -> Plot only the cells/particles that belong to the main halo.
      + `:subhalo`         -> Plot only the cells/particles that belong to the main subhalo.
      + `:sphere`          -> Plot only the cell/particle inside a sphere with radius `FILTER_R` (see `./src/constants.jl`).
      + `:stellar_subhalo` -> Plot only the cells/particles that belong to the main subhalo.
  - `sim_labels::Union{Vector{String},Nothing}=basename.(simulation_paths)`: Labels for the plot legend, one per simulation. Set it to `nothing` if you don't want a legend.
"""
function densityProfile(
    simulation_paths::Vector{String},
    slice_n::Int,
    quantity::Symbol;
    output_path::String="./",
    filter_mode::Symbol=:all,
    sim_labels::Union{Vector{String},Nothing}=basename.(simulation_paths),
)::Nothing

    plot_params = plotParams(quantity)
    request = addRequest(plot_params.request, Dict(:gas => ["VEL "], :stars => ["VEL "]))
    filter_function, translation, rotation, request = selectFilter(filter_mode, plot_params.request)

    grid = CircularGrid(FILTER_R, 100)

    # Draw the figures with CairoMakie
    snapshotPlot(
        simulation_paths,
        request,
        [lines!];
        pf_kwargs=[(;)],
        # `snapshotPlot` configuration
        output_path,
        base_filename="$(quantity)-profile",
        output_format=".png",
        warnings=true,
        show_progress=false,
        # Data manipulation options
        slice=slice_n,
        filter_function,
        da_functions=[daDensityProfile],
        da_args=[(grid, quantity)],
        da_kwargs=[(;)],
        post_processing=getNothing,
        pp_args=(),
        pp_kwargs=(;),
        transform_box=true,
        translation,
        rotation,
        smooth=0,
        x_unit=u"kpc",
        y_unit=plot_params.unit,
        x_exp_factor=0,
        y_exp_factor=0,
        x_trim=(-Inf, Inf),
        y_trim=(-Inf, Inf),
        x_edges=false,
        y_edges=false,
        x_func=identity,
        y_func=identity,
        # Axes options
        xaxis_label="auto_label",
        yaxis_label=plot_params.axis_label,
        xaxis_var_name=L"r",
        yaxis_var_name=plot_params.var_name,
        xaxis_scale_func=identity,
        yaxis_scale_func=identity,
        xaxis_limits=(nothing, nothing),
        yaxis_limits=(nothing, nothing),
        # Plotting and animation options
        save_figures=true,
        backup_results=false,
        sim_labels,
        title=:physical_time,
        pt_per_unit=0.75,
        px_per_unit=2.0,
        resolution=(1280, 800),
        aspect=nothing,
        series_colors=nothing,
        series_markers=nothing,
        series_linestyles=nothing,
        # Animation options
        animation=false,
        animation_filename="animation.mp4",
        framerate=10,
    )

    return nothing

end

"""
    stellarHistory(
        simulation_paths::Vector{String},
        slice_n::Int,
        quantity::Symbol;
        <keyword arguments>
    )::Nothing

Plot the evolution of a given stellar `quantity` using the stellar ages at a given instant in time.

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`.
  - `slice_n::Int`: Selects which snapshot to plot, starts at 1 and is independent of the number in the file name. If every snapshot is present, `slice_n` = filename_number + 1.
  - `quantity::Symbol`: Quantity for the y axis. The options are:

      + `:sfr`          -> The star formation rate.
      + `:ssfr`         -> The specific star formation rate.
      + `:stellar_mass` -> Stellar mass.
  - `n_bins::Int=20`: Number of bins (time intervals).
  - `output_path::String="./"`: Path to the output folder.
  - `filter_mode::Symbol=:all`: Which cells/particles will be plotted, the options are:

      + `:all`             -> Plot every cell/particle within the simulation box.
      + `:halo`            -> Plot only the cells/particles that belong to the main halo.
      + `:subhalo`         -> Plot only the cells/particles that belong to the main subhalo.
      + `:sphere`          -> Plot only the cell/particle inside a sphere with radius `FILTER_R` (see `./src/constants.jl`).
      + `:stellar_subhalo` -> Plot only the cells/particles that belong to the main subhalo.
  - `sim_labels::Union{Vector{String},Nothing}=basename.(simulation_paths)`: Labels for the plot legend, one per simulation. Set it to `nothing` if you don't want a legend.
"""
function stellarHistory(
    simulation_paths::Vector{String},
    slice_n::Int,
    quantity::Symbol;
    n_bins::Int=20,
    output_path::String="./",
    filter_mode::Symbol=:all,
    sim_labels::Union{Vector{String},Nothing}=basename.(simulation_paths),
)::Nothing

    x_plot_params = plotParams(:physical_time)
    y_plot_params = plotParams(quantity)

    filter_function, translation, rotation, request = selectFilter(
        filter_mode,
        Dict(:stars => ["GAGE"]),
    )

    # Draw the figures with CairoMakie
    snapshotPlot(
        simulation_paths,
        request,
        [scatterlines!];
        pf_kwargs=[(;)],
        # `snapshotPlot` configuration
        output_path,
        base_filename="$(quantity)-stellar-history",
        output_format=".png",
        warnings=true,
        show_progress=false,
        # Data manipulation options
        slice=slice_n,
        filter_function,
        da_functions=[daStellarHistory],
        da_args=[()],
        da_kwargs=[(; quantity, n_bins)],
        post_processing=getNothing,
        pp_args=(),
        pp_kwargs=(;),
        transform_box=true,
        translation,
        rotation,
        smooth=0,
        x_unit=x_plot_params.unit,
        y_unit=y_plot_params.unit,
        x_exp_factor=x_plot_params.exp_factor,
        y_exp_factor=y_plot_params.exp_factor,
        x_trim=(-Inf, Inf),
        y_trim=(-Inf, Inf),
        x_edges=false,
        y_edges=false,
        x_func=identity,
        y_func=identity,
        # Axes options
        xaxis_label=x_plot_params.axis_label,
        yaxis_label=y_plot_params.axis_label,
        xaxis_var_name=x_plot_params.var_name,
        yaxis_var_name=y_plot_params.var_name,
        xaxis_scale_func=identity,
        yaxis_scale_func=log10,
        xaxis_limits=(nothing, nothing),
        yaxis_limits=(nothing, nothing),
        # Plotting and animation options
        save_figures=true,
        backup_results=false,
        sim_labels,
        title=:physical_time,
        pt_per_unit=0.75,
        px_per_unit=2.0,
        resolution=(1280, 800),
        aspect=nothing,
        series_colors=nothing,
        series_markers=nothing,
        series_linestyles=nothing,
        # Animation options
        animation=false,
        animation_filename="animation.mp4",
        framerate=10,
    )

    return nothing

end

"""
    stellarCircularity(
        simulation_paths::Vector{String},
        slice_n::Int;
        <keyword arguments>
    )::Nothing

Plot a histogram of the stellar circularity.

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`.
  - `slice_n::Int`: Selects which snapshot to plot, starts at 1 and is independent of the number in the file name. If every snapshot is present, `slice_n` = filename_number + 1.
  - `range::NTuple{2,<:Number}=(-2.0, 2.0)`: Circularity range.
  - `n_bins::Int=60`: Number of bins.
  - `output_path::String="./"`: Path to the output folder.
  - `filter_mode::Symbol=:all`: Which cells/particles will be plotted, the options are:

      + `:all`             -> Plot every cell/particle within the simulation box.
      + `:halo`            -> Plot only the cells/particles that belong to the main halo.
      + `:subhalo`         -> Plot only the cells/particles that belong to the main subhalo.
      + `:sphere`          -> Plot only the cell/particle inside a sphere with radius `FILTER_R` (see `./src/constants.jl`).
      + `:stellar_subhalo` -> Plot only the cells/particles that belong to the main subhalo.
  - `sim_labels::Union{Vector{String},Nothing}=basename.(simulation_paths)`: Labels for the plot legend, one per simulation. Set it to `nothing` if you don't want a legend.
"""
function stellarCircularity(
    simulation_paths::Vector{String},
    slice_n::Int;
    range::NTuple{2,<:Number}=(-2.0, 2.0),
    n_bins::Int=60,
    output_path::String="./",
    filter_mode::Symbol=:all,
    sim_labels::Union{Vector{String},Nothing}=basename.(simulation_paths),
)::Nothing

    plot_params = plotParams(:stellar_circularity)

    filter_function, translation, rotation, request = selectFilter(filter_mode, plot_params.request)

    grid = LinearGrid(range..., n_bins)

    snapshotPlot(
        simulation_paths,
        request,
        [lines!];
        pf_kwargs=[(;)],
        # `snapshotPlot` configuration
        output_path,
        base_filename="circularity_histogram",
        output_format=".png",
        warnings=true,
        show_progress=false,
        # Data manipulation options
        slice=slice_n,
        filter_function,
        da_functions=[daCircularityHistogram],
        da_args=[(grid,)],
        da_kwargs=[(;)],
        post_processing=getNothing,
        pp_args=(),
        pp_kwargs=(;),
        transform_box=true,
        translation,
        rotation,
        smooth=0,
        x_unit=plot_params.unit,
        y_unit=Unitful.NoUnits,
        x_exp_factor=plot_params.exp_factor,
        y_exp_factor=0,
        x_trim=(-Inf, Inf),
        y_trim=(-Inf, Inf),
        x_edges=false,
        y_edges=false,
        x_func=identity,
        y_func=identity,
        # Axes options
        xaxis_label=plot_params.axis_label,
        yaxis_label="auto_label",
        xaxis_var_name=plot_params.var_name,
        yaxis_var_name=L"\mathrm{Normalized \,\, counts}",
        xaxis_scale_func=identity,
        yaxis_scale_func=identity,
        xaxis_limits=(nothing, nothing),
        yaxis_limits=(nothing, nothing),
        # Plotting and animation options
        save_figures=true,
        backup_results=false,
        sim_labels,
        title=:physical_time,
        pt_per_unit=0.75,
        px_per_unit=2.0,
        resolution=(1000, 1000),
        aspect=AxisAspect(1),
        series_colors=nothing,
        series_markers=nothing,
        series_linestyles=nothing,
        # Animation options
        animation=false,
        animation_filename="animation.mp4",
        framerate=10,
    )

    return nothing

end

"""
    compareWithFeldmann2020(
        simulation_paths::Vector{String},
        x_quantity::Symbol,
        y_quantity::Symbol;
        <keyword arguments>
    )::Nothing

Plot a time series plus the corresponding experimental results from Feldmann (2020).

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`.
  - `x_quantity::Symbol`: Quantity for the x axis. The possibilities are:

      + `:stellar_mass`   -> Stellar mass.
      + `:molecular_mass` -> Molecular hydrogen (``\\mathrm{H_2}``) mass.
      + `:atomic_mass`    -> Atomic hydrogen (``\\mathrm{HI}``) mass.
      + `:sfr`            -> The star formation rate of the last `AGE_RESOLUTION`.
  - `y_quantity::Symbol`: Quantity for the y axis. The possibilities are:

      + `:stellar_mass`   -> Stellar mass.
      + `:molecular_mass` -> Molecular hydrogen (``\\mathrm{H_2}``) mass.
      + `:atomic_mass`    -> Atomic hydrogen (``\\mathrm{HI}``) mass.
      + `:sfr`            -> The star formation rate of the last `AGE_RESOLUTION`.
  - `slice::IndexType=(:)`: Slice of the simulation, i.e. which snapshots will be read. It can be an integer (a single snapshot), a vector of integers (several snapshots), an `UnitRange` (e.g. 5:13), an `StepRange` (e.g. 5:2:13) or (:) (all snapshots). It works over the longest possible list of snapshots among the simulations (grouped by the number in the file names). Out of bounds indices are ignored.
  - `output_path::String="./"`: Path to the output folder.
  - `filter_mode::Symbol=:all`: Which cells/particles will be plotted, the options are:

      + `:all`             -> Plot every cell/particle within the simulation box.
      + `:halo`            -> Plot only the cells/particles that belong to the main halo.
      + `:subhalo`         -> Plot only the cells/particles that belong to the main subhalo.
      + `:sphere`          -> Plot only the cell/particle inside a sphere with radius `FILTER_R` (see `./src/constants.jl`).
      + `:stellar_subhalo` -> Plot only the cells/particles that belong to the main subhalo.
  - `sim_labels::Union{Vector{String},Nothing}=basename.(simulation_paths)`: Labels for the plot legend, one per simulation. Set it to `nothing` if you don't want a legend.

# References

R. Feldmann (2020). *The link between star formation and gas in nearby galaxies*. Communications Physics **3(226)**. [doi:10.1038/s42005-020-00493-0](https://doi.org/10.1038/s42005-020-00493-0)
"""
function compareWithFeldmann2020(
    simulation_paths::Vector{String},
    x_quantity::Symbol,
    y_quantity::Symbol;
    slice::IndexType=(:),
    output_path::String="./",
    filter_mode::Symbol=:all,
    sim_labels::Union{Vector{String},Nothing}=basename.(simulation_paths),
)::Nothing

    x_plot_params = plotParams(x_quantity)
    y_plot_params = plotParams(y_quantity)

    timeSeriesPlot(
        simulation_paths,
        [scatterlines!];
        pf_kwargs=[(; markersize=14)],
        # `timeSeriesPlot` configuration
        output_path,
        filename="$(y_quantity)-vs-$(x_quantity)-with-Feldmann2020",
        output_format=".png",
        warnings=true,
        show_progress=true,
        # Data manipulation options
        slice,
        da_functions=[daEvolution],
        da_args=[(x_quantity, y_quantity)],
        da_kwargs=[(; filter_mode, smooth=0, scaling=identity, warnings=true)],
        post_processing=ppFeldmann2020!,
        pp_args=(x_quantity, y_quantity),
        pp_kwargs=(;),
        x_unit=x_plot_params.unit,
        y_unit=y_plot_params.unit,
        x_exp_factor=0,
        y_exp_factor=0,
        x_trim=(-Inf, Inf),
        y_trim=(-Inf, Inf),
        x_edges=false,
        y_edges=false,
        x_func=identity,
        y_func=identity,
        # Axes options
        xaxis_label=x_plot_params.axis_label,
        yaxis_label=y_plot_params.axis_label,
        xaxis_var_name=x_plot_params.var_name,
        yaxis_var_name=y_plot_params.var_name,
        xaxis_scale_func=log10,
        yaxis_scale_func=log10,
        xaxis_limits=(nothing, nothing),
        yaxis_limits=(nothing, nothing),
        # Plotting options
        save_figure=true,
        backup_results=false,
        sim_labels,
        title="",
        pt_per_unit=0.75,
        px_per_unit=2.0,
        resolution=(1280, 800),
        aspect=nothing,
        series_colors=nothing,
        series_markers=nothing,
        series_linestyles=nothing,
    )

    return nothing

end

"""
    compareWithMolla2015(
        simulation_paths::Vector{String},
        slice_n::Int,
        quantity::Symbol;
        <keyword arguments>
    )::Nothing

Plot a Milky Way profile plus the corresponding experimental results from Mollá et al. (2015).

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`.
  - `slice_n::Int`: Selects which snapshot to plot, starts at 1 and is independent of the number in the file name. If every snapshot is present, `slice_n` = filename_number + 1.
  - `quantity::Symbol`: Quantity for the y axis. The options are:

      + `:stellar_area_density`   -> Stellar area mass density.
      + `:molecular_area_density` -> Molecular hydrogen area mass density.
      + `:atomic_area_density`    -> Atomic hydrogen area mass density.
      + `:sfr_area_density`       -> Star formation rate area density, for the last `AGE_RESOLUTION_ρ`.
      + `:O_stellar_abundance`    -> Stellar abundance of oxygen, as ``12 + \\log_{10}(\\mathrm{O \\, / \\, H})``.
      + `:N_stellar_abundance`    -> Stellar abundance of nitrogen, as ``12 + \\log_{10}(\\mathrm{N \\, / \\, H})``.
      + `:C_stellar_abundance`    -> Stellar abundance of carbon, as ``12 + \\log_{10}(\\mathrm{C \\, / \\, H})``.
  - `output_path::String="./"`: Path to the output folder.
  - `filter_mode::Symbol=:all`: Which cells/particles will be plotted, the options are:

      + `:all`             -> Plot every cell/particle within the simulation box.
      + `:halo`            -> Plot only the cells/particles that belong to the main halo.
      + `:subhalo`         -> Plot only the cells/particles that belong to the main subhalo.
      + `:sphere`          -> Plot only the cell/particle inside a sphere with radius `FILTER_R` (see `./src/constants.jl`).
      + `:stellar_subhalo` -> Plot only the cells/particles that belong to the main subhalo.
  - `sim_labels::Union{Vector{String},Nothing}=basename.(simulation_paths)`: Labels for the plot legend, one per simulation. Set it to `nothing` if you don't want a legend.

# References

M. Mollá et al. (2015). *Galactic chemical evolution: stellar yields and the initial mass function*. Monthly Notices of the Royal Astronomical Society **451(4)**, 3693–3708. [doi:10.1093/mnras/stv1102](https://doi.org/10.1093/mnras/stv1102)
"""
function compareWithMolla2015(
    simulation_paths::Vector{String},
    slice_n::Int,
    quantity::Symbol;
    output_path::String="./",
    filter_mode::Symbol=:all,
    sim_labels::Union{Vector{String},Nothing}=basename.(simulation_paths),
)::Nothing

    plot_params = plotParams(quantity)
    request = addRequest(plot_params.request, Dict(:gas => ["VEL "], :stars => ["VEL "]))
    filter_function, translation, rotation, request = selectFilter(filter_mode, plot_params.request)

    grid = CircularGrid(FILTER_R, 20)

    # Draw the figures with CairoMakie
    snapshotPlot(
        simulation_paths,
        request,
        [scatterlines!];
        pf_kwargs=[(; markersize=15)],
        # `snapshotPlot` configuration
        output_path,
        base_filename="$(quantity)-profile-with-Molla2015",
        output_format=".png",
        warnings=true,
        show_progress=false,
        # Data manipulation options
        slice=slice_n,
        filter_function,
        da_functions=[daMolla2015],
        da_args=[(grid, quantity)],
        da_kwargs=[(;)],
        post_processing=ppMolla2015!,
        pp_args=(quantity,),
        pp_kwargs=(; color=:red, linestyle=nothing, error_bars=true),
        transform_box=true,
        translation,
        rotation,
        smooth=0,
        x_unit=u"kpc",
        y_unit=plot_params.unit,
        x_exp_factor=0,
        y_exp_factor=0,
        x_trim=(-Inf, Inf),
        y_trim=(-Inf, Inf),
        x_edges=false,
        y_edges=false,
        x_func=identity,
        y_func=identity,
        # Axes options
        xaxis_label="auto_label",
        yaxis_label=plot_params.axis_label,
        xaxis_var_name=L"r",
        yaxis_var_name=plot_params.var_name,
        xaxis_scale_func=identity,
        yaxis_scale_func=identity,
        xaxis_limits=(nothing, nothing),
        yaxis_limits=(nothing, nothing),
        # Plotting and animation options
        save_figures=true,
        backup_results=false,
        sim_labels,
        title=:physical_time,
        pt_per_unit=0.75,
        px_per_unit=2.0,
        resolution=(1280, 800),
        aspect=nothing,
        series_colors=nothing,
        series_markers=nothing,
        series_linestyles=nothing,
        # Animation options
        animation=false,
        animation_filename="animation.mp4",
        framerate=10,
    )

    return nothing

end

"""
    compareWithKennicuttBigiel(
        simulation_paths::Vector{String},
        slice_n::Int;
        <keyword arguments>
    )::Nothing

Plot the resolved Kennicutt-Schmidt relation plus the results of Kennicutt (1998) or Bigiel et al. (2008), depending on the chosen `quantity`.

!!! note

    This method plots the KS relation using cylindrical bins at a fix moment in time.

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`.
  - `slice_n::Int`: Selects which snapshot to plot, starts at 1 and is independent of the number in the file name. If every snapshot is present, `slice_n` = filename_number + 1.
  - `quantity::Symbol=:molecular_area_density`: Quantity for the x axis. The possibilities are:

      + `:gas_area_density`       -> Gas area mass density, for a radius of `FILTER_R`. This one will be plotted with the results of Kennicutt (1998).
      + `:molecular_area_density` -> Molecular hydrogen area mass density, for a radius of `FILTER_R`. This one will be plotted with the results of Bigiel et al. (2008).
      + `:neutral_area_density`   -> Neutral hydrogen area mass density, for a radius of `FILTER_R`. This one will be plotted with the results of Bigiel et al. (2008).
  - `output_path::String="./"`: Path to the output folder.
  - `filter_mode::Symbol=:all`: Which cells/particles will be plotted, the options are:

      + `:all`             -> Plot every cell/particle within the simulation box.
      + `:halo`            -> Plot only the cells/particles that belong to the main halo.
      + `:subhalo`         -> Plot only the cells/particles that belong to the main subhalo.
      + `:sphere`          -> Plot only the cell/particle inside a sphere with radius `FILTER_R` (see `./src/constants.jl`).
      + `:stellar_subhalo` -> Plot only the cells/particles that belong to the main subhalo.
  - `sim_labels::Union{Vector{String},Nothing}=basename.(simulation_paths)`: Labels for the plot legend, one per simulation. Set it to `nothing` if you don't want a legend.

# References

R. C. Kennicutt (1998). *The Global Schmidt Law in Star-forming Galaxies*. The Astrophysical Journal, **498(2)**, 541-552. [doi:10.1086/305588](https://doi.org/10.1086/305588)

F. Bigiel et al. (2008). *THE STAR FORMATION LAW IN NEARBY GALAXIES ON SUB-KPC SCALES*. The Astrophysical Journal, **136(6)**, 2846. [doi:10.1088/0004-6256/136/6/2846](https://doi.org/10.1088/0004-6256/136/6/2846)
"""
function compareWithKennicuttBigiel(
    simulation_paths::Vector{String},
    slice_n::Int;
    quantity::Symbol=:molecular_area_density,
    output_path::String="./",
    filter_mode::Symbol=:all,
    sim_labels::Union{Vector{String},Nothing}=basename.(simulation_paths),
)::Nothing

    grid = CircularGrid(FILTER_R, 20)

    if quantity == :gas_area_density

        da_functions = [daKennicuttSchmidt]
        da_args = [(grid,)]
        post_processing = ppKennicutt1998!
        pp_args = ()
        filename = "Kennicutt1998"

    else

        da_functions = [daKennicuttSchmidtLaw]
        da_args = [(grid, quantity)]
        post_processing = ppBigiel2008!
        pp_args = (quantity,)
        filename = "Bigiel2008"

    end

    x_plot_params = plotParams(quantity)
    y_plot_params = plotParams(:sfr_area_density)

    filter_function, translation, rotation, request = selectFilter(
        filter_mode,
        mergeRequests(x_plot_params.request, y_plot_params.request),
    )

    snapshotPlot(
        simulation_paths,
        request,
        [scatter!];
        pf_kwargs=[(; markersize=15)],
        # `snapshotPlot` configuration
        output_path,
        base_filename="sfr_area_density-vs-$(quantity)-with-$(filename)",
        output_format=".png",
        warnings=true,
        show_progress=false,
        # Data manipulation options
        slice=slice_n,
        filter_function,
        da_functions,
        da_args,
        da_kwargs=[(;)],
        post_processing,
        pp_args,
        pp_kwargs=(;
            x_unit=x_plot_params.unit,
            y_unit=y_plot_params.unit,
            x_log=false,
            y_log=false,
            color=:red,
            linestyle=nothing,
        ),
        transform_box=true,
        translation,
        rotation,
        smooth=0,
        x_unit=x_plot_params.unit,
        y_unit=y_plot_params.unit,
        x_exp_factor=0,
        y_exp_factor=0,
        x_trim=(-Inf, Inf),
        y_trim=(-Inf, Inf),
        x_edges=false,
        y_edges=false,
        x_func=identity,
        y_func=identity,
        # Axes options
        xaxis_label=x_plot_params.axis_label,
        yaxis_label=y_plot_params.axis_label,
        xaxis_var_name=x_plot_params.var_name,
        yaxis_var_name=y_plot_params.var_name,
        xaxis_scale_func=log10,
        yaxis_scale_func=log10,
        xaxis_limits=(nothing, nothing),
        yaxis_limits=(nothing, nothing),
        # Plotting and animation options
        save_figures=true,
        backup_results=false,
        sim_labels,
        title=:physical_time,
        pt_per_unit=0.75,
        px_per_unit=2.0,
        resolution=(1280, 800),
        aspect=nothing,
        series_colors=nothing,
        series_markers=nothing,
        series_linestyles=nothing,
        # Animation options
        animation=false,
        animation_filename="animation.mp4",
        framerate=10,
    )

    return nothing

end

"""
    compareWithKennicuttBigiel(
        simulation_paths::Vector{String},
        slice::Union{Colon,UnitRange{<:Integer},StepRange{<:Integer,<:Integer},Vector{<:Integer}};
        <keyword arguments>
    )::Nothing

Plot the integrated Kennicutt-Schmidt relation plus the results of Kennicutt (1998) or Bigiel et al. (2008), depending on the chosen `quantity`.

!!! note

    This method plots the KS relation for the whole galaxy at different points in time.

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`.
  - `slice::Union{Colon,UnitRange{<:Integer},StepRange{<:Integer,<:Integer},Vector{<:Integer}}`: Slice of the simulation, i.e. which snapshots will be read. It can be a vector of integers (several snapshots), an `UnitRange` (e.g. 5:13), an `StepRange` (e.g. 5:2:13) or (:) (all snapshots). It works over the longest possible list of snapshots among the simulations (grouped by the number in the file names). Out of bounds indices are ignored.
  - `quantity::Symbol=:molecular_area_density`: Quantity for the x axis. The possibilities are:

      + `:gas_area_density`       -> Gas area mass density, for a radius of `FILTER_R`. This one will be plotted with the results of Kennicutt (1998).
      + `:molecular_area_density` -> Molecular hydrogen area mass density, for a radius of `FILTER_R`. This one will be plotted with the results of Bigiel et al. (2008).
      + `:neutral_area_density`   -> Neutral hydrogen area mass density, for a radius of `FILTER_R`. This one will be plotted with the results of Bigiel et al. (2008).
  - `output_path::String="./"`: Path to the output folder.
  - `filter_mode::Symbol=:all`: Which cells/particles will be plotted, the options are:

      + `:all`             -> Plot every cell/particle within the simulation box.
      + `:halo`            -> Plot only the cells/particles that belong to the main halo.
      + `:subhalo`         -> Plot only the cells/particles that belong to the main subhalo.
      + `:sphere`          -> Plot only the cell/particle inside a sphere with radius `FILTER_R` (see `./src/constants.jl`).
      + `:stellar_subhalo` -> Plot only the cells/particles that belong to the main subhalo.
  - `sim_labels::Union{Vector{String},Nothing}=basename.(simulation_paths)`: Labels for the plot legend, one per simulation. Set it to `nothing` if you don't want a legend.

# References

R. C. Kennicutt (1998). *The Global Schmidt Law in Star-forming Galaxies*. The Astrophysical Journal, **498(2)**, 541-552. [doi:10.1086/305588](https://doi.org/10.1086/305588)

F. Bigiel et al. (2008). *THE STAR FORMATION LAW IN NEARBY GALAXIES ON SUB-KPC SCALES*. The Astrophysical Journal, **136(6)**, 2846. [doi:10.1088/0004-6256/136/6/2846](https://doi.org/10.1088/0004-6256/136/6/2846)
"""
function compareWithKennicuttBigiel(
    simulation_paths::Vector{String},
    slice::Union{Colon,UnitRange{<:Integer},StepRange{<:Integer,<:Integer},Vector{<:Integer}};
    quantity::Symbol=:molecular_area_density,
    output_path::String="./",
    filter_mode::Symbol=:all,
    sim_labels::Union{Vector{String},Nothing}=basename.(simulation_paths),
)::Nothing

    if quantity == :gas_area_density

        post_processing = ppKennicutt1998!
        pp_args = ()
        filename = "Kennicutt1998"

    else

        post_processing = ppBigiel2008!
        pp_args = (quantity,)
        filename = "Bigiel2008"

    end

    x_plot_params = plotParams(quantity)
    y_plot_params = plotParams(:sfr_area_density)

    # Draw the figures with CairoMakie
    timeSeriesPlot(
        simulation_paths,
        [scatter!];
        pf_kwargs=[(; markersize=15)],
        # `timeSeriesPlot` configuration
        output_path,
        filename="sfr_area_density-vs-$(quantity)-with-$(filename)",
        output_format=".png",
        warnings=true,
        show_progress=true,
        # Data manipulation options
        slice,
        da_functions=[daEvolution],
        da_args=[(quantity, :sfr_area_density)],
        da_kwargs=[(; filter_mode, smooth=0, scaling=identity, warnings=true)],
        post_processing,
        pp_args,
        pp_kwargs=(;
            x_unit=x_plot_params.unit,
            y_unit=y_plot_params.unit,
            x_log=false,
            y_log=false,
            color=:red,
            linestyle=nothing,
        ),
        x_unit=x_plot_params.unit,
        y_unit=y_plot_params.unit,
        x_exp_factor=0,
        y_exp_factor=0,
        x_trim=(-Inf, Inf),
        y_trim=(-Inf, Inf),
        x_edges=false,
        y_edges=false,
        x_func=identity,
        y_func=identity,
        # Axes options
        xaxis_label=x_plot_params.axis_label,
        yaxis_label=y_plot_params.axis_label,
        xaxis_var_name=x_plot_params.var_name,
        yaxis_var_name=y_plot_params.var_name,
        xaxis_scale_func=log10,
        yaxis_scale_func=log10,
        xaxis_limits=(nothing, nothing),
        yaxis_limits=(nothing, nothing),
        # Plotting options
        save_figure=true,
        backup_results=false,
        sim_labels,
        title="",
        pt_per_unit=0.75,
        px_per_unit=2.0,
        resolution=(1280, 800),
        aspect=nothing,
        series_colors=nothing,
        series_markers=nothing,
        series_linestyles=nothing,
    )

    return nothing

end
