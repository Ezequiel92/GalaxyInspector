####################################################################################################
# Opinionated convenience functions
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
  - `filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all`: Which cells/particles will be considered in the "filtered" section of the report. The options are:

      + `:all`             -> Consider every cell/particle within the simulation box.
      + `:halo`            -> Consider only the cells/particles that belong to the main halo.
      + `:subhalo`         -> Consider only the cells/particles that belong to the main subhalo.
      + `:sphere`          -> Consider only the cell/particle inside a sphere with radius `DISK_R` (see `./src/constants/globals.jl`).
      + `:stellar_subhalo` -> Consider only the cells/particles that belong to the main subhalo.
      + `:all_subhalo`     -> Plot every cell/particle centered around the main subhalo.
      + A dictionary with three entries:

          + `:filter_function` -> The filter function.
          + `:translation`     -> Translation for the simulation box. The posibilities are:

              + `:global_cm`                  -> Selects the center of mass of the whole system as the new origin.
              + `:{component}`                -> Sets the center of mass of the given component (e.g. :stars, :gas, :halo, etc, after filtering) as the new origin. It can be any of the keys of [`PARTICLE_INDEX`](@ref).
              + `(halo_idx, subhalo_rel_idx)` -> Sets the position of the potencial minimum for the `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo) as the new origin.
              + `(halo_idx, 0)`               -> Sets the center of mass of the `halo_idx::Int` halo as the new origin.
              + `subhalo_abs_idx`             -> Sets the center of mass of the `subhalo_abs_idx::Int` as the new origin.
          + `:rotation`        -> Rotation for the simulation box. The posibilities are:

              + `:zero`                       -> No rotation is appplied.
              + `:global_am`                  -> Sets the angular momentum of the whole system as the new z axis.
              + `:stellar_am`                 -> Sets the stellar angular momentum as the new z axis.
              + `:stellar_pa`                 -> Sets the stellar principal axis as the new coordinate system.
              + `:stellar_subhalo_pa`         -> Sets the principal axis of the stars in the main subhalo as the new coordinate system.
              + `(halo_idx, subhalo_rel_idx)` -> Sets the principal axis of the stars in `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo), as the new coordinate system.
              + `(halo_idx, 0)`               -> Sets the principal axis of the stars in the `halo_idx::Int` halo, as the new coordinate system.
              + `subhalo_abs_idx`             -> Sets the principal axis of the stars in the `subhalo_abs_idx::Int` subhalo as the new coordinate system.
  - `halo_idx::Int=1`: Index of the target halo (FoF group). Starts at 1.
  - `subhalo_rel_idx::Int=1`: Index of the target subhalo (subfind), relative to the target halo. Starts at 1.
  - `warnings::Bool=true`: If a warning will be given when there is missing files.
"""
function snapshotReport(
    simulation_paths::Vector{String},
    slice_n::Int;
    output_path::String="./",
    filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all,
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

        if !PHYSICAL_UNITS && cosmological
            println(file, "Report units:     Comoving\n")
        else
            println(file, "Report units:     Physical\n")
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
        filter_function, translation, rotation, request = selectFilter(
            filter_mode,
            mergeRequests(
                Dict(component => ["POS ", "MASS", "VEL "] for component in component_list),
                Dict(
                    :gas => ["NHP ", "NH  ", "PRES", "FRAC", "CTIM", "TAUS", "ID  "],
                    :stars => ["ACIT", "PARZ", "RHOC", "ID  "],
                ),
            ),
        )

        # Read the necessary snapshot data
        if !in(filter_mode, [:all, :sphere])

            # Check that the group catalog data is available
            if !ismissing(groupcat_path) && isSubfindActive(groupcat_path)

                # Create the data dictionary
                data_dict = makeDataDict(
                    simulation_path,
                    slice_n,
                    request;
                    warnings,
                )

            else

                throw(ArgumentError("snapshotReport: You asked for a filter base on the \
                halos/subhalos or a perzonalized filter, but I could not find a valid \
                group catalog file"))

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
            title *= " "^(25 - length(title))

            println(file, "\t\t$(title)$(count)")

        end

        println(file, "\n\t\tTotal count:             $(total_count)\n")

        ############################################################################################
        # Print the total mass of each component
        ############################################################################################

        println(file, "\tMasses:\n")

        total_mass = 0.0u"Msun"
        for component in component_list

            mass = sum(data_dict[component]["MASS"]; init=0.0u"Msun")

            total_mass += mass

            title = "$(PARTICLE_NAMES[component]):"
            title *= " "^(25 - length(title))

            println(file, "\t\t$(title)$(round(typeof(1.0u"Msun"), mass, sigdigits=3))")

        end

        println(
            file,
            "\n\t\tTotal mass:              $(round(typeof(1.0u"Msun"), total_mass, sigdigits=3))\n",
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

            h2p_mass = sum(computePressureMolecularMass(data_dict); init=0.0u"Msun")
            h2p_percent = round((h2p_mass / gas_mass) * 100, sigdigits=3)

            hn_mass = sum(computeNeutralMass(data_dict); init=0.0u"Msun")
            hn_percent = round((hn_mass / gas_mass) * 100, sigdigits=3)

            if !iszero(hii_mass)
                title = "Ionized mass:"
                title *= " "^(25 - length(title))
                println(
                    file,
                    "\t\t$(title)$(round(typeof(1.0u"Msun"), hii_mass, sigdigits=3)) \
                    ($(hii_percent)% of total gas mass)",
                )
            end

            if !iszero(hi_mass)
                title = "Atomic mass:"
                title *= " "^(25 - length(title))
                println(
                    file,
                    "\t\t$(title)$(round(typeof(1.0u"Msun"), hi_mass, sigdigits=3)) \
                    ($(hi_percent)% of total gas mass)",
                )
            end

            if !iszero(h2_mass)
                title = "Molecular mass:"
                title *= " "^(25 - length(title))
                println(
                    file,
                    "\t\t$(title)$(round(typeof(1.0u"Msun"), h2_mass, sigdigits=3)) \
                    ($(h2_percent)% of total gas mass)",
                )
            end

            if !iszero(h2p_mass)
                title = "Molecular mass (BR):"
                title *= " "^(25 - length(title))
                println(
                    file,
                    "\t\t$(title)$(round(typeof(1.0u"Msun"), h2p_mass, sigdigits=3)) \
                    ($(h2p_percent)% of total gas mass)",
                )
            end

            if !iszero(hn_mass)
                title = "Neutral mass:"
                title *= " "^(25 - length(title))
                println(
                    file,
                    "\t\t$(title)$(round(typeof(1.0u"Msun"), hn_mass, sigdigits=3)) \
                    ($(hn_percent)% of total gas mass)\n",
                )
            end

        end

        ############################################################################################
        # Filter the simulation box
        ############################################################################################

        filterData!(data_dict; filter_function)

        println(file, "#"^100)
        println(file, "Filtered box with:")

        if filter_mode isa Symbol
            println(file, "\n\tFilter mode: $(filter_mode)")
        else
            println(file, "\n\tFilter function: $(String(Symbol(filter_function)))")
        end

        println(file, "\tTranslation: $(translation)")
        println(file, "\tRotation: $(rotation)")
        println(file, "#"^100)

        ############################################################################################
        # Print the global properties of the simulation after filtering
        ############################################################################################

        println(file, "\nGlobal properties:\n")

        ############################################################################################
        # Print the number of cells/particles for each component
        ############################################################################################

        println(file, "\tCell/particle number:\n")

        total_count = 0
        for component in component_list

            count = length(data_dict[component]["MASS"])

            title = "$(PARTICLE_NAMES[component]):"
            title *= " "^(25 - length(title))

            total_count += count

            println(file, "\t\t$(title)$(count)")

        end

        println(file, "\n\t\tTotal count:             $(total_count)\n")

        ############################################################################################
        # Print the mass of each component
        ############################################################################################

        println(file, "\tMasses:\n")

        total_mass = 0.0u"Msun"
        for component in component_list

            mass = sum(data_dict[component]["MASS"]; init=0.0u"Msun")

            total_mass += mass

            title = "$(PARTICLE_NAMES[component]):"
            title *= " "^(25 - length(title))

            println(file, "\t\t$(title)$(round(typeof(1.0u"Msun"), mass, sigdigits=3))")

        end

        println(
            file,
            "\n\t\tTotal mass:              $(round(typeof(1.0u"Msun"), total_mass, sigdigits=3))\n",
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

            h2p_mass = sum(computePressureMolecularMass(data_dict); init=0.0u"Msun")
            h2p_percent = round((h2p_mass / gas_mass) * 100, sigdigits=3)

            hn_mass = sum(computeNeutralMass(data_dict); init=0.0u"Msun")
            hn_percent = round((hn_mass / gas_mass) * 100, sigdigits=3)

            if !iszero(hii_mass)
                title = "Ionized mass:"
                title *= " "^(25 - length(title))
                println(
                    file,
                    "\t\t$(title)$(round(typeof(1.0u"Msun"), hii_mass, sigdigits=3)) \
                    ($(hii_percent)% of total gas mass)",
                )
            end

            if !iszero(hi_mass)
                title = "Atomic mass:"
                title *= " "^(25 - length(title))
                println(
                    file,
                    "\t\t$(title)$(round(typeof(1.0u"Msun"), hi_mass, sigdigits=3)) \
                    ($(hi_percent)% of total gas mass)",
                )
            end

            if !iszero(h2_mass)
                title = "Molecular mass:"
                title *= " "^(25 - length(title))
                println(
                    file,
                    "\t\t$(title)$(round(typeof(1.0u"Msun"), h2_mass, sigdigits=3)) \
                    ($(h2_percent)% of total gas mass)",
                )
            end

            if !iszero(h2p_mass)
                title = "Molecular mass (BR):"
                title *= " "^(25 - length(title))
                println(
                    file,
                    "\t\t$(title)$(round(typeof(1.0u"Msun"), h2p_mass, sigdigits=3)) \
                    ($(h2p_percent)% of total gas mass)",
                )
            end

            if !iszero(hn_mass)
                title = "Neutral mass:"
                title *= " "^(25 - length(title))
                println(
                    file,
                    "\t\t$(title)$(round(typeof(1.0u"Msun"), hn_mass, sigdigits=3)) \
                    ($(hn_percent)% of total gas mass)\n",
                )
            end

        end

        ############################################################################################
        # Print the center of mass of each component
        ############################################################################################

        println(file, "\tCenter of mass:\n")

        global_cm = computeGlobalCenterOfMass(data_dict)

        for component in component_list

            cm = computeCenterOfMass(data_dict[component]["POS "], data_dict[component]["MASS"])

            title = "$(PARTICLE_NAMES[component]):"
            title *= " "^(25 - length(title))

            println(file, "\t\t$(title)$(round.(ustrip.(u"Mpc", cm), sigdigits=6)) $(u"Mpc")")
            println(file, "\t\tDistance to global CM:   $(sqrt(sum((global_cm - cm).^2)))\n")

        end

        global_cm = round.(ustrip.(u"Mpc", global_cm), sigdigits=6)
        println(file, "\t\tGlobal center of mass:   $(global_cm) $(u"Mpc")\n")

        ############################################################################################
        # Print the fraction of gas cells that have enter our SF routine
        ############################################################################################

        if !isempty(data_dict[:gas]["FRAC"])

            total_number = length(data_dict[:gas]["MASS"])
            stellar_gas_number = count(!isnan, data_dict[:gas]["FRAC"][1, :])
            fraction = (stellar_gas_number / total_number) * 100

            idxs = findall(!isnan, data_dict[:gas]["FRAC"][1, :])
            stellar_gas_mass = sum(data_dict[:gas]["MASS"][idxs])
            mass_fraction = (stellar_gas_mass / sum(data_dict[:gas]["MASS"])) * 100

            println(file, "\tFraction of gas cells that have enter our SF routine:\n")
            println(file, "\t\t$(round(fraction, sigdigits=3))% of the cells")
            println(file, "\t\t$(round(mass_fraction, sigdigits=3))% of the mass\n")

        end

        ############################################################################################
        # Print the properties of the star forming gas
        ############################################################################################

        if any(
            !isempty,
            [data_dict[:stars]["PARZ"], data_dict[:stars]["RHOC"], data_dict[:stars]["ACIT"]],
        )

            println(file, "\tProperties of the star forming gas:\n")

        end

        if !isempty(data_dict[:stars]["PARZ"])

            parz = data_dict[:stars]["PARZ"] ./ SOLAR_METALLICITY

            println(file, "\t\tMetallicity:\n")
            println(file, "\t\t\tMean:    $(round(mean(parz), sigdigits=4)) Z⊙")
            println(file, "\t\t\tMedian:  $(round(median(parz), sigdigits=4)) Z⊙")
            println(file, "\t\t\tMode:    $(round(mode(parz)[1], sigdigits=4)) Z⊙")
            println(file, "\t\t\tMinimum: $(round(minimum(parz), sigdigits=4)) Z⊙")
            println(file, "\t\t\tMaximum: $(round(maximum(parz), sigdigits=4)) Z⊙\n")

            parz_50 = computeMassQty(parz, data_dict[:stars]["MASS"]; percent=50.0)
            parz_90 = computeMassQty(parz, data_dict[:stars]["MASS"]; percent=90.0)
            parz_95 = computeMassQty(parz, data_dict[:stars]["MASS"]; percent=95.0)

            println(file, "\t\tMetallicity enclosing X% of the stellar mass:\n")
            println(file, "\t\t\t$(round(parz_50, sigdigits=4)) Z⊙ (50%)")
            println(file, "\t\t\t$(round(parz_90, sigdigits=4)) Z⊙ (90%)")
            println(file, "\t\t\t$(round(parz_95, sigdigits=4)) Z⊙ (95%)\n")

        end

        if !isempty(data_dict[:stars]["RHOC"])

            rhoc = ustrip.(u"cm^-3", data_dict[:stars]["RHOC"])

            println(file, "\t\tCell density:\n")
            println(file, "\t\t\tMean:    $(round(mean(rhoc), sigdigits=4)) cm^-3")
            println(file, "\t\t\tMedian:  $(round(median(rhoc), sigdigits=4)) cm^-3")
            println(file, "\t\t\tMode:    $(round(mode(rhoc)[1], sigdigits=4)) cm^-3")
            println(file, "\t\t\tMinimum: $(round(minimum(rhoc), sigdigits=4)) cm^-3")
            println(file, "\t\t\tMaximum: $(round(maximum(rhoc), sigdigits=4)) cm^-3\n")

            rhoc_50 = computeMassQty(rhoc, data_dict[:stars]["MASS"]; percent=50.0)
            rhoc_90 = computeMassQty(rhoc, data_dict[:stars]["MASS"]; percent=90.0)
            rhoc_95 = computeMassQty(rhoc, data_dict[:stars]["MASS"]; percent=95.0)

            println(file, "\t\tCell density enclosing X% of the stellar mass:\n")
            println(file, "\t\t\t$(round(rhoc_50, sigdigits=4)) cm^-3 (50%)")
            println(file, "\t\t\t$(round(rhoc_90, sigdigits=4)) cm^-3 (90%)")
            println(file, "\t\t\t$(round(rhoc_95, sigdigits=4)) cm^-3 (95%)\n")

        end

        if !isempty(data_dict[:stars]["ACIT"])

            acit = ustrip.(u"Myr", data_dict[:stars]["ACIT"])

            println(file, "\t\tTotal integration time:\n")
            println(file, "\t\t\tMean:    $(round(mean(acit), sigdigits=4)) Myr")
            println(file, "\t\t\tMedian:  $(round(median(acit), sigdigits=4)) Myr")
            println(file, "\t\t\tMode:    $(round(mode(acit)[1], sigdigits=4)) Myr")
            println(file, "\t\t\tMinimum: $(round(minimum(acit), sigdigits=4)) Myr")
            println(file, "\t\t\tMaximum: $(round(maximum(acit), sigdigits=4)) Myr\n")

            acit_50 = computeMassQty(acit, data_dict[:stars]["MASS"]; percent=50.0)
            acit_90 = computeMassQty(acit, data_dict[:stars]["MASS"]; percent=90.0)
            acit_95 = computeMassQty(acit, data_dict[:stars]["MASS"]; percent=95.0)

            println(file, "\t\tTotal integration time enclosing X% of the stellar mass:\n")
            println(file, "\t\t\t$(round(acit_50, sigdigits=4)) Myr (50%)")
            println(file, "\t\t\t$(round(acit_90, sigdigits=4)) Myr (90%)")
            println(file, "\t\t\t$(round(acit_95, sigdigits=4)) Myr (95%)\n")

        end

        ############################################################################################
        # Print the maximum and minimum values of each parameter of the ODEs
        ############################################################################################

        quantities = [
            "ODIT",
            "ACIT",
            "DTIM",
            "TAUS",
            "RHOC",
            "PARZ",
            "ETAD",
            "ETAI",
            "PARR",
        ]
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

        if any(isBlockPresent.(:gas, quantities, snapshot_path))
            println(file, "\tExtrema of the ODEs ICs and parameters:\n")
        end

        for (quantity, unit, name) in zip(quantities, units, names)

            if isBlockPresent(:gas, quantity, snapshot_path)

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
        # Translate the simulation box
        ############################################################################################

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
            title *= " "^(25 - length(title))

            println(file, "\t\t$(title)$(round.(L, sigdigits=3))")

        end

        global_L = round.(computeGlobalAngularMomentum(data_dict), sigdigits=3)

        println(file, "\n\t\tGlobal angular momentum: $(global_L)\n")

        ############################################################################################
        # Print the spin parameter of each component
        ############################################################################################

        println(file, "\tSpin parameter (R = $(DISK_R)):\n")

        for component in component_list

            λ = computeSpinParameter(
                data_dict[component]["POS "],
                data_dict[component]["VEL "],
                data_dict[component]["MASS"],
            )

            title = "$(PARTICLE_NAMES[component]):"
            title *= " "^(25 - length(title))

            println(file, "\t\t$(title)$(round.(λ, sigdigits=3))")

        end

        global_λ = round.(computeGlobalSpinParameter(data_dict), sigdigits=3)

        println(file, "\n\t\tTotal spin parameter:    $(global_λ)\n")

        ############################################################################################
        # Rotate the simulation box
        ############################################################################################

        rotateData!(data_dict, rotation)

        ############################################################################################
        # Print the total height of a cylinder, of infinite radius, containing 90% and 95%
        # of the stellar mass
        ############################################################################################

        mass_height_90 = computeMassHeight(
            data_dict[:stars]["POS "],
            data_dict[:stars]["MASS"];
            percent=90.0,
        )

        mass_height_95 = computeMassHeight(
            data_dict[:stars]["POS "],
            data_dict[:stars]["MASS"];
            percent=95.0,
        )

        println(file, "\tTotal height containing X% of the stellar mass:\n")
        println(file, "\t\t$(round(ustrip(u"kpc", mass_height_90), sigdigits=4)) $(u"kpc") (90%)")
        println(file, "\t\t$(round(ustrip(u"kpc", mass_height_95), sigdigits=4)) $(u"kpc") (95%)\n")

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
                :subhalo => [
                    "S_Mass",
                    "S_MassType",
                    "S_LenType",
                    "S_CM",
                    "S_Pos",
                    "S_Vel",
                    "S_HalfmassRad",

                ],
                :group => [
                    "G_Mass",
                    "G_MassType",
                    "G_M_Crit200",
                    "G_LenType",
                    "G_Nsubs",
                    "G_CM",
                    "G_Pos",
                    "G_Vel",
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
            s_vel           = gc_data[:subhalo]["S_Vel"][:, subhalo_abs_idx]
            s_half_mass_rad = gc_data[:subhalo]["S_HalfmassRad"][subhalo_abs_idx]
            g_mass          = gc_data[:group]["G_Mass"][halo_idx]
            g_mass_type     = gc_data[:group]["G_MassType"][:, halo_idx]
            g_m_crit_200    = gc_data[:group]["G_M_Crit200"][halo_idx]
            g_len_type      = gc_data[:group]["G_LenType"][:, halo_idx]
            g_n_subs        = gc_data[:group]["G_Nsubs"][halo_idx]
            g_cm            = gc_data[:group]["G_CM"][:, halo_idx]
            g_pos           = gc_data[:group]["G_Pos"][:, halo_idx]
            g_vel           = gc_data[:group]["G_Vel"][:, halo_idx]
            g_r_crit_200    = gc_data[:group]["G_R_Crit200"][halo_idx]

            ########################################################################################
            # Print the fraction of insitu stars
            ########################################################################################

            if snapshot_length >= 2

                insitu_idx = filterInsituStars(
                    data_dict;
                    halo_idx,
                    subhalo_rel_idx,
                    warnings,
                )[:stars]

                iMs = sum(data_dict[:stars]["MASS"][insitu_idx]; init=0.0u"Msun")
                tMs = sum(data_dict[:stars]["MASS"]; init=0.0u"Msun")
                insitu_fraction = round(uconvert.(Unitful.NoUnits, (iMs / tMs) * 100); sigdigits=2)

                println(file, "#"^100)
                println(file, "\nFraction of insitu stars: $(insitu_fraction)%\n")

            end

            ########################################################################################
            # Print the mass and clumping factor of each hydrogen phase
            # between `DISK_R` and the virial radius
            ########################################################################################

            if :gas in component_list

                ####################################################################################
                # Indices of cells and particles within the disc radius and
                # between the disc radius and the virial radius
                ####################################################################################

                disc_idxs = filterWithinSphere(data_dict, (0.0u"kpc", DISK_R), :zero)
                halo_idxs = filterWithinSphere(data_dict, (DISK_R, g_r_crit_200), :zero)

                stellar_masses = data_dict[:stars]["MASS"]
                gas_masses     = data_dict[:gas]["MASS"]

                ionized_masses     = computeIonizedMass(data_dict)
                atomic_masses      = computeAtomicMass(data_dict)
                molecular_masses   = computeMolecularMass(data_dict)
                molecular_P_masses = computePressureMolecularMass(data_dict)
                neutral_masses     = computeNeutralMass(data_dict)

                stellar_mass_inside  = stellar_masses[disc_idxs[:stars]]
                stellar_mass_outside = stellar_masses[halo_idxs[:stars]]

                gas_mass_inside  = gas_masses[disc_idxs[:gas]]
                gas_mass_outside = gas_masses[halo_idxs[:gas]]

                if !isempty(ionized_masses)
                    ionized_mass_inside  = ionized_masses[disc_idxs[:gas]]
                    ionized_mass_outside = ionized_masses[halo_idxs[:gas]]
                end

                if !isempty(atomic_masses)
                    atomic_mass_inside  = atomic_masses[disc_idxs[:gas]]
                    atomic_mass_outside = atomic_masses[halo_idxs[:gas]]
                end

                if !isempty(molecular_masses)
                    molecular_mass_inside  = molecular_masses[disc_idxs[:gas]]
                    molecular_mass_outside = molecular_masses[halo_idxs[:gas]]
                end

                if !isempty(molecular_P_masses)
                    molecular_P_mass_inside  = molecular_P_masses[disc_idxs[:gas]]
                    molecular_P_mass_outside = molecular_P_masses[halo_idxs[:gas]]
                end

                if !isempty(neutral_masses)
                    neutral_mass_inside  = neutral_masses[disc_idxs[:gas]]
                    neutral_mass_outside = neutral_masses[halo_idxs[:gas]]
                end

                println(file, "Characteristic radii:\n")

                ####################################################################################
                # Print the radius containing 90% and 95% of the mass,
                # withing de disc (r < `DISK_R`)
                ####################################################################################

                ########
                # Stars
                ########

                mass_radius_90 = computeMassRadius(
                    data_dict[:stars]["POS "][:, disc_idxs[:stars]],
                    stellar_mass_inside;
                    percent=90.0,
                )

                mass_radius_95 = computeMassRadius(
                    data_dict[:stars]["POS "][:, disc_idxs[:stars]],
                    stellar_mass_inside;
                    percent=95.0,
                )

                println(file, "\tRadius containing X% of the stellar mass (r < $(DISK_R)):\n")
                println(
                    file,
                    "\t\t$(round(ustrip(u"kpc", mass_radius_90), sigdigits=4)) $(u"kpc") (90%)",
                )
                println(
                    file,
                    "\t\t$(round(ustrip(u"kpc", mass_radius_95), sigdigits=4)) $(u"kpc") (95%)\n",
                )

                ############
                # Total gas
                ############

                mass_radius_90 = computeMassRadius(
                    data_dict[:gas]["POS "][:, disc_idxs[:gas]],
                    gas_mass_inside;
                    percent=90.0,
                )

                mass_radius_95 = computeMassRadius(
                    data_dict[:gas]["POS "][:, disc_idxs[:gas]],
                    gas_mass_inside;
                    percent=95.0,
                )

                println(
                    file,
                    "\tRadius containing X% of the total gas mass (r < $(DISK_R)):\n",
                )
                println(
                    file,
                    "\t\t$(round(ustrip(u"kpc", mass_radius_90), sigdigits=4)) $(u"kpc") (90%)",
                )
                println(
                    file,
                    "\t\t$(round(ustrip(u"kpc", mass_radius_95), sigdigits=4)) $(u"kpc") (95%)\n",
                )

                ##############
                # Ionized gas
                ##############

                if !isempty(ionized_masses)

                    mass_radius_90 = computeMassRadius(
                        data_dict[:gas]["POS "][:, disc_idxs[:gas]],
                        ionized_mass_inside;
                        percent=90.0,
                    )

                    mass_radius_95 = computeMassRadius(
                        data_dict[:gas]["POS "][:, disc_idxs[:gas]],
                        ionized_mass_inside;
                        percent=95.0,
                    )

                    println(
                        file,
                        "\tRadius containing X% of the ionized gas mass (r < $(DISK_R)):\n",
                    )
                    println(
                        file,
                        "\t\t$(round(ustrip(u"kpc", mass_radius_90), sigdigits=4)) $(u"kpc") (90%)",
                    )
                    println(
                        file,
                        "\t\t$(round(ustrip(u"kpc", mass_radius_95), sigdigits=4)) $(u"kpc") (95%)\n",
                    )

                end

                #############
                # Atomic gas
                #############

                if !isempty(atomic_masses)

                    mass_radius_90 = computeMassRadius(
                        data_dict[:gas]["POS "][:, disc_idxs[:gas]],
                        atomic_mass_inside;
                        percent=90.0,
                    )

                    mass_radius_95 = computeMassRadius(
                        data_dict[:gas]["POS "][:, disc_idxs[:gas]],
                        atomic_mass_inside;
                        percent=95.0,
                    )

                    println(
                        file,
                        "\tRadius containing X% of the atomic gas mass (r < $(DISK_R)):\n",
                    )
                    println(
                        file,
                        "\t\t$(round(ustrip(u"kpc", mass_radius_90), sigdigits=4)) $(u"kpc") (90%)",
                    )
                    println(
                        file,
                        "\t\t$(round(ustrip(u"kpc", mass_radius_95), sigdigits=4)) $(u"kpc") (95%)\n",
                    )

                end

                ################
                # Molecular gas
                ################

                if !isempty(molecular_masses)

                    mass_radius_90 = computeMassRadius(
                        data_dict[:gas]["POS "][:, disc_idxs[:gas]],
                        molecular_mass_inside;
                        percent=90.0,
                    )

                    mass_radius_95 = computeMassRadius(
                        data_dict[:gas]["POS "][:, disc_idxs[:gas]],
                        molecular_mass_inside;
                        percent=95.0,
                    )

                    println(
                        file,
                        "\tRadius containing X% of the molecular gas mass (r < $(DISK_R)):\n",
                    )
                    println(
                        file,
                        "\t\t$(round(ustrip(u"kpc", mass_radius_90), sigdigits=4)) $(u"kpc") (90%)",
                    )
                    println(
                        file,
                        "\t\t$(round(ustrip(u"kpc", mass_radius_95), sigdigits=4)) $(u"kpc") (95%)\n",
                    )

                end

                #####################
                # Molecular gas (BR)
                #####################

                if !isempty(molecular_P_masses)

                    mass_radius_90 = computeMassRadius(
                        data_dict[:gas]["POS "][:, disc_idxs[:gas]],
                        molecular_P_mass_inside;
                        percent=90.0,
                    )

                    mass_radius_95 = computeMassRadius(
                        data_dict[:gas]["POS "][:, disc_idxs[:gas]],
                        molecular_P_mass_inside;
                        percent=95.0,
                    )

                    println(
                        file,
                        "\tRadius containing X% of the molecular gas mass (BR) (r < $(DISK_R)):\n",
                    )
                    println(
                        file,
                        "\t\t$(round(ustrip(u"kpc", mass_radius_90), sigdigits=4)) $(u"kpc") (90%)",
                    )
                    println(
                        file,
                        "\t\t$(round(ustrip(u"kpc", mass_radius_95), sigdigits=4)) $(u"kpc") (95%)\n",
                    )

                end

                ##############
                # Neutral gas
                ##############

                if !isempty(neutral_masses)

                    mass_radius_90 = computeMassRadius(
                        data_dict[:gas]["POS "][:, disc_idxs[:gas]],
                        neutral_mass_inside;
                        percent=90.0,
                    )

                    mass_radius_95 = computeMassRadius(
                        data_dict[:gas]["POS "][:, disc_idxs[:gas]],
                        neutral_mass_inside;
                        percent=95.0,
                    )

                    println(
                        file,
                        "\tRadius containing X% of the neutral gas mass (r < $(DISK_R)):\n",
                    )
                    println(
                        file,
                        "\t\t$(round(ustrip(u"kpc", mass_radius_90), sigdigits=4)) $(u"kpc") (90%)",
                    )
                    println(
                        file,
                        "\t\t$(round(ustrip(u"kpc", mass_radius_95), sigdigits=4)) $(u"kpc") (95%)\n",
                    )

                end

                ####################################################################################
                # Print the masses withing de disc (r < DISK_R) and
                # outside the disc (DISK_R < r < R200)
                ####################################################################################

                println(file, "#"^100)
                println(file, "\nCharacteristic fractions and masses:\n")

                println(file, "\t", "#"^20)
                println(file, "\tR200: $(round(typeof(1.0u"kpc"), g_r_crit_200, sigdigits=4))")
                println(file, "\t", "#"^20, "\n")

                ########
                # Stars
                ########

                total_stellar_mass_inside  = sum(stellar_mass_inside; init=0.0u"Msun")
                total_stellar_mass_outside = sum(stellar_mass_outside; init=0.0u"Msun")
                total_stellar_mass         = total_stellar_mass_inside + total_stellar_mass_outside

                s_inside_percent  = (total_stellar_mass_inside / total_stellar_mass) * 100.0
                s_outside_percent = (total_stellar_mass_outside / total_stellar_mass) * 100.0

                println(file, "\t", "#"^40)
                println(file, "\tStars:")
                println(file, "\t", "#"^40, "\n")

                println(file, "\tStellar mass inside the disc (r < $(DISK_R)):\n")
                println(
                    file,
                    "\t\t$(round(typeof(1.0u"Msun"), total_stellar_mass_inside, sigdigits=3)) \
                    ($(round(s_inside_percent, sigdigits=3))% of the total stellar mass)\n",
                )

                println(file, "\tStellar mass outside the disc ($(DISK_R) < r < R200):\n")
                println(
                    file,
                    "\t\t$(round(typeof(1.0u"Msun"), total_stellar_mass_outside, sigdigits=3)) \
                    ($(round(s_outside_percent, sigdigits=3))% of the total stellar mass)\n",
                )

                ############
                # Total gas
                ############

                total_gas_mass_inside  = sum(gas_mass_inside; init=0.0u"Msun")
                total_gas_mass_outside = sum(gas_mass_outside; init=0.0u"Msun")
                total_gas_mass         = total_gas_mass_inside + total_gas_mass_outside

                g_inside_percent  = (total_gas_mass_inside / total_gas_mass) * 100.0
                g_outside_percent = (total_gas_mass_outside / total_gas_mass) * 100.0

                println(file, "\t", "#"^40)
                println(file, "\tTotal gas:")
                println(file, "\t", "#"^40, "\n")

                println(file, "\tGas mass inside the disc (r < $(DISK_R)):\n")
                println(
                    file,
                    "\t\t$(round(typeof(1.0u"Msun"), total_gas_mass_inside, sigdigits=3)) \
                    ($(round(g_inside_percent, sigdigits=3))% of the total gas mass)\n",
                )

                println(file, "\tGas mass outside the disc ($(DISK_R) < r < R200):\n")
                println(
                    file,
                    "\t\t$(round(typeof(1.0u"Msun"), total_gas_mass_outside, sigdigits=3)) \
                    ($(round(g_outside_percent, sigdigits=3))% of the total gas mass)\n",
                )

                ##############
                # Ionized gas
                ##############

                if !isempty(ionized_masses)

                    total_ion_mass_inside  = sum(ionized_mass_inside; init=0.0u"Msun")
                    total_ion_mass_outside = sum(ionized_mass_outside; init=0.0u"Msun")

                    i_inside_percent  = (total_ion_mass_inside  / total_gas_mass) * 100.0
                    i_outside_percent = (total_ion_mass_outside / total_gas_mass) * 100.0

                    println(file, "\t", "#"^40)
                    println(file, "\tIonized gas:")
                    println(file, "\t", "#"^40, "\n")

                    println(file, "\tIonized mass inside the disc (r < $(DISK_R)):\n")
                    println(
                        file,
                        "\t\t$(round(typeof(1.0u"Msun"), total_ion_mass_inside, sigdigits=3)) \
                        ($(round(i_inside_percent, sigdigits=3))% of the total gas mass)\n",
                    )

                    println(file, "\tIonized mass outside the disc ($(DISK_R) < r < R200):\n")
                    println(
                        file,
                        "\t\t$(round(typeof(1.0u"Msun"), total_ion_mass_outside, sigdigits=3)) \
                        ($(round(i_outside_percent, sigdigits=3))% of the total gas mass)\n",
                    )

                end

                #############
                # Atomic gas
                #############

                if !isempty(atomic_masses)

                    total_ato_mass_inside  = sum(atomic_mass_inside; init=0.0u"Msun")
                    total_ato_mass_outside = sum(atomic_mass_outside; init=0.0u"Msun")

                    a_inside_percent  = (total_ato_mass_inside  / total_gas_mass) * 100.0
                    a_outside_percent = (total_ato_mass_outside / total_gas_mass) * 100.0

                    println(file, "\t", "#"^40)
                    println(file, "\tAtomic gas:")
                    println(file, "\t", "#"^40, "\n")

                    println(file, "\tAtomic mass inside the disc (r < $(DISK_R)):\n")
                    println(
                        file,
                        "\t\t$(round(typeof(1.0u"Msun"), total_ato_mass_inside, sigdigits=3)) \
                        ($(round(a_inside_percent, sigdigits=3))% of the total gas mass)\n",
                    )

                    println(file, "\tAtomic mass outside the disc ($(DISK_R) < r < R200):\n")
                    println(
                        file,
                        "\t\t$(round(typeof(1.0u"Msun"), total_ato_mass_outside, sigdigits=3)) \
                        ($(round(a_outside_percent, sigdigits=3))% of the total gas mass)\n",
                    )

                end

                ################
                # Molecular gas
                ################

                if !isempty(molecular_masses)

                    total_mol_mass_inside  = sum(molecular_mass_inside; init=0.0u"Msun")
                    total_mol_mass_outside = sum(molecular_mass_outside; init=0.0u"Msun")

                    m_inside_percent  = (total_mol_mass_inside  / total_gas_mass) * 100.0
                    m_outside_percent = (total_mol_mass_outside / total_gas_mass) * 100.0

                    println(file, "\t", "#"^40)
                    println(file, "\tMolecular gas:")
                    println(file, "\t", "#"^40, "\n")

                    println(file, "\tMolecular mass inside the disc (r < $(DISK_R)):\n")
                    println(
                        file,
                        "\t\t$(round(typeof(1.0u"Msun"), total_mol_mass_inside, sigdigits=3)) \
                        ($(round(m_inside_percent, sigdigits=3))% of the total gas mass)\n",
                    )

                    println(file, "\tMolecular mass outside the disc ($(DISK_R) < r < R200):\n")
                    println(
                        file,
                        "\t\t$(round(typeof(1.0u"Msun"), total_mol_mass_outside, sigdigits=3)) \
                        ($(round(m_outside_percent, sigdigits=3))% of the total gas mass)\n",
                    )

                end

                #####################
                # Molecular gas (BR)
                #####################

                if !isempty(molecular_P_masses)

                    total_mol_P_mass_inside  = sum(molecular_P_mass_inside; init=0.0u"Msun")
                    total_mol_P_mass_outside = sum(molecular_P_mass_outside; init=0.0u"Msun")

                    m_P_inside_percent  = (total_mol_P_mass_inside  / total_gas_mass) * 100.0
                    m_P_outside_percent = (total_mol_P_mass_outside / total_gas_mass) * 100.0

                    println(file, "\t", "#"^40)
                    println(file, "\tMolecular gas (BR recipe):")
                    println(file, "\t", "#"^40, "\n")

                    println(file, "\tMolecular mass (BR recipe) inside the disc (r < $(DISK_R)):\n")
                    println(
                        file,
                        "\t\t$(round(typeof(1.0u"Msun"), total_mol_P_mass_inside, sigdigits=3)) \
                        ($(round(m_P_inside_percent, sigdigits=3))% of the total gas mass)\n",
                    )

                    println(file, "\tMolecular mass (BR recipe) outside the disc ($(DISK_R) < r < R200):\n")
                    println(
                        file,
                        "\t\t$(round(typeof(1.0u"Msun"), total_mol_P_mass_outside, sigdigits=3)) \
                        ($(round(m_P_outside_percent, sigdigits=3))% of the total gas mass)\n",
                    )

                end

                ##############
                # Neutral gas
                ##############

                if !isempty(neutral_masses)

                    total_neu_mass_inside  = sum(neutral_mass_inside; init=0.0u"Msun")
                    total_neu_mass_outside = sum(neutral_mass_outside; init=0.0u"Msun")

                    m_inside_percent  = (total_neu_mass_inside  / total_gas_mass) * 100.0
                    m_outside_percent = (total_neu_mass_outside / total_gas_mass) * 100.0

                    println(file, "\t", "#"^40)
                    println(file, "\tNeutral gas:")
                    println(file, "\t", "#"^40, "\n")

                    println(file, "\tNeutral mass inside the disc (r < $(DISK_R)):\n")
                    println(
                        file,
                        "\t\t$(round(typeof(1.0u"Msun"), total_neu_mass_inside, sigdigits=3)) \
                        ($(round(m_inside_percent, sigdigits=3))% of the total gas mass)\n",
                    )

                    println(file, "\tNeutral mass outside the disc ($(DISK_R) < r < R200):\n")
                    println(
                        file,
                        "\t\t$(round(typeof(1.0u"Msun"), total_neu_mass_outside, sigdigits=3)) \
                        ($(round(m_outside_percent, sigdigits=3))% of the total gas mass)\n",
                    )

                end

            end

            ########################################################################################
            # Halo and subhalo global properties
            ########################################################################################

            println(file, "#"^100)
            println(file, "\nHalo and subhalo global properties:")

            ############################
            # Print the halo properties
            ############################

            println(file, "\n", "#"^71)
            println(file, "NOTE: Stellar particle counts include wind particles from here on out!")
            println(file, "#"^71)

            println(file, "\nHalo $(lpad(halo_idx - 1, 3, "0")) properties:\n")

            ########################################################################################

            println(file, "\tCell/particle number:\n")
            for (i, len) in pairs(g_len_type)

                component = PARTICLE_NAMES[INDEX_PARTICLE[i - 1]]
                println(file, "\t\t$(component):$(" "^(22 - length(component))) $(len)")

            end

            ########################################################################################

            println(file, "\n\tNumber of subhalos:\n\n\t\t$(g_n_subs)\n")

            ########################################################################################

            println(file, "\tMasses:\n")
            for (i, mass) in pairs(g_mass_type)

                component = INDEX_PARTICLE[i - 1]

                if component == :stars
                    component = "Stellar/Wind particles"
                else
                    component = PARTICLE_NAMES[component]
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

            separation = sqrt(sum((g_cm - g_pos) .^ 2))
            println(
                file,
                "\tSeparation between the minimum potencial and the global CM: \
                \n\n\t\t$(round(typeof(1.0u"kpc"), separation, sigdigits=6))\n",
            )

            ########################################################################################

            vel_cm = round.(Float64.(ustrip.(u"km*s^-1", g_vel)), sigdigits=6)
            println(file, "\tVelocity of the center of mass:\n\n\t\t$(vel_cm) $(u"km*s^-1")\n")

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

            ###############################
            # Print the subhalo properties
            ###############################

            println(file, "#"^100)
            println(
                file,
                "\nSubhalo $(lpad(subhalo_rel_idx - 1, 3, "0")) (of halo \
                $(lpad(halo_idx - 1, 3, "0"))) properties:\n",
            )

            ########################################################################################

            println(file, "\tCell/particle number:\n")
            for (i, len) in pairs(s_len_type)

                component = PARTICLE_NAMES[INDEX_PARTICLE[i - 1]]
                println(file, "\t\t$(component):$(" "^(22 - length(component))) $(len)")

            end

            ########################################################################################

            println(file, "\n\tMasses:\n")
            for (i, mass) in pairs(s_mass_type)

                component = INDEX_PARTICLE[i - 1]

                if component == :stars
                    component = "Stellar/Wind particles"
                else
                    component = PARTICLE_NAMES[component]
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

            separation = sqrt(sum((s_cm - s_pos) .^ 2))
            println(
                file,
                "\tSeparation between the minimum potencial and the global CM: \
                \n\n\t\t$(round(typeof(1.0u"kpc"), separation, sigdigits=6))\n",
            )

            ########################################################################################

            vel_cm = round.(Float64.(ustrip.(u"km*s^-1", s_vel)), sigdigits=6)
            println(file, "\tVelocity of the center of mass:\n\n\t\t$(vel_cm) $(u"km*s^-1")\n")

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

        if !PHYSICAL_UNITS && cosmological
            println(file, "Report units:     Comoving\n")
        else
            println(file, "Report units:     Physical\n")
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
                println(file, "Redshift range:            $(max_z) - $(min_z)")

            end

            println(file,)

        else

            pt = round(ustrip(u"Gyr", simulation_table[1, :physical_times]), digits=2)

            println(file, "Physical time:             $(pt) Gyr")

            if cosmological

                # For cosmological simulations print the scale factor and the redshift
                a = round(simulation_table[1, :scale_factors], digits=3)
                z = round(simulation_table[1, :redshifts], digits=3)

                println(file, "Scale factor:              $(a)")
                println(file, "Redshift:                  $(z)")

            end

            println(file,)

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
                println(file, "\nFirst snapshot with star formation:")

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

                    separation = sqrt(sum((s_cm - s_pos) .^ 2))
                    println(
                        file,
                        "\t\tSeparation between the minimum potencial and the global CM: \
                        \n\n\t\t\t$(round(typeof(1.0u"kpc"), separation, sigdigits=6))\n",
                    )

                else
                    println(file, "\n" * "#"^51)
                    println(file, "There is no subfind information for this snapshot!")
                    println(file, "#"^51 * "\n")
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
                println(file, "\nFirst snapshot with subfind information:")

                println(file, "\n\tSnapshot:         $(basename(snapshot_path))")

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

                separation = sqrt(sum((s_cm - s_pos) .^ 2))
                println(
                    file,
                    "\t\tSeparation between the minimum potencial and the global CM: \
                    \n\n\t\t\t$(round(typeof(1.0u"kpc"), separation, sigdigits=6))\n",
                )

                println(file, "\t\t", "#"^54)
                println(file, "\t\tNOTE: Stellar particle counts include wind particles!")
                println(file, "\t\t", "#"^54)

                println(file, "\n\t\tCell/particle number:\n")

                for (i, len) in pairs(s_len_type)

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
                println(file, "\nFirst snapshot with star formation in the main subhalo:")

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
  - `x_quantity::Symbol`: Quantity for the x axis. The options are:

      + `:physical_time` -> Physical time since the Big Bang.
      + `:lookback_time` -> Physical time left to reach the last snapshot.
      + `:scale_factor`  -> Scale factor.
      + `:redshift`      -> Redshift.
  - `y_quantity::Symbol`: Quantity for the y axis. The options are:

      + `:stellar_mass` -> Stellar mass.
      + `:sfr`          -> The star formation rate.
  - `smooth::Int=0`: The result will be smooth out using `smooth` bins. Set it to 0 if you want no smoothing.
  - `output_path::String="./"`: Path to the output folder.
  - `sim_labels::Union{Vector{String},Nothing}=basename.(simulation_paths)`: Labels for the plot legend, one per simulation. Set it to `nothing` if you don't want a legend.
  - `theme::Attributes=Theme()`: Plot theme that will take precedence over [`DEFAULT_THEME`](@ref).
"""
function sfrTXT(
    simulation_paths::Vector{String},
    x_quantity::Symbol,
    y_quantity::Symbol;
    smooth::Int=0,
    output_path::String="./",
    sim_labels::Union{Vector{String},Nothing}=basename.(simulation_paths),
    theme::Attributes=Theme(),
)::Nothing

    x_plot_params = plotParams(x_quantity)
    y_plot_params = plotParams(y_quantity)

    plotTimeSeries(
        simulation_paths,
        [lines!];
        pf_kwargs=[(;)],
        # `plotTimeSeries` configuration
        output_path,
        filename="$(y_quantity)-vs-$(x_quantity)",
        output_format=".png",
        warnings=false,
        show_progress=true,
        # Data manipulation options
        slice=(:),
        da_functions=[daSFRtxt],
        da_args=[(x_quantity, y_quantity)],
        da_kwargs=[(; smooth, warnings=false)],
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
        yaxis_scale_func=identity,
        # Plotting options
        save_figures=true,
        backup_results=false,
        theme,
        sim_labels,
        title="",
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
  - `x_quantity::Symbol`: Quantity for the x axis. The options are:

      + `:time_step`              -> Time step.
      + `:physical_time`          -> Physical time since the Big Bang.
      + `:clock_time_s`           -> Clock time duration of the time step in seconds.
      + `:clock_time_percent`     -> Clock time duration of the time step as a percentage.
      + `:tot_clock_time_s`       -> Total clock time in seconds.
      + `:tot_clock_time_percent` -> Total clock time as a percentage.
  - `y_quantity::Symbol`: Quantity for the y axis. The options are:

      + `:time_step`              -> Time step.
      + `:physical_time`          -> Physical time since the Big Bang.
      + `:clock_time_s`           -> Clock time duration of the time step in seconds.
      + `:clock_time_percent`     -> Clock time duration of the time step as a percentage.
      + `:tot_clock_time_s`       -> Total clock time in seconds.
      + `:tot_clock_time_percent` -> Total clock time as a percentage.
  - `smooth::Int=0`: The result will be smooth out using `smooth` bins. Set it to 0 if you want no smoothing.
  - `yscale::Function=identity`: Scaling function for the y axis. The options are the scaling functions accepted by [Makie](https://docs.makie.org/stable/): log10, log2, log, sqrt, Makie.logit, Makie.Symlog10, Makie.pseudolog10, and identity.
  - `x_trim::NTuple{2,<:Real}=(-Inf, Inf)`: The data will be trim down so the x coordinates fit within `x_trim`.
  - `y_trim::NTuple{2,<:Real}=(-Inf, Inf)`: The data will be trim down so the y coordinates fit within `y_trim`. This option does not affect histograms.
  - `output_path::String="./"`: Path to the output folder.
  - `sim_labels::Union{Vector{String},Nothing}=basename.(simulation_paths)`: Labels for the plot legend, one per simulation. Set it to `nothing` if you don't want a legend.
  - `theme::Attributes=Theme()`: Plot theme that will take precedence over [`DEFAULT_THEME`](@ref).
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
    theme::Attributes=Theme(),
)::Nothing

    x_plot_params = plotParams(x_quantity)
    y_plot_params = plotParams(y_quantity)

    safe_str_target = replace(target, "/" => "-", "_" => "-")

    plotTimeSeries(
        simulation_paths,
        [lines!];
        pf_kwargs=[(;)],
        # `plotTimeSeries` configuration
        output_path,
        filename="$(y_quantity)-vs-$(x_quantity)-for-$(safe_str_target)",
        output_format=".png",
        warnings=false,
        show_progress=true,
        # Data manipulation options
        slice=(:),
        da_functions=[daCPUtxt],
        da_args=[(target, x_quantity, y_quantity)],
        da_kwargs=[(; smooth, warnings=false)],
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
        # Plotting options
        save_figures=true,
        backup_results=false,
        theme,
        sim_labels,
        title=L"\mathrm{Process: \,\, %$(safe_str_target)}",
    )

    return nothing

end

"""
    stellarBirthHalos(
        simulation_path::String,
        slice_n::Int;
        <keyword arguments>
    )::Nothing

Write, to a pair of CSV files, in which halo and subhalo every star in snapshot `slice_n` was born.

# Arguments

  - `simulation_paths::String`: Path to the simulation directory, set in the code variable `OutputDir`.
  - `slice_n::Int`: Selects the target snapshot. Starts at 1 and is independent of the number in the file name. If every snapshot is present, `slice_n` = filename_number + 1.
  - `output_path::String="./"`: Path to the output folder.
  - `filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all`: Which cells/particles will be plotted, the options are:

      + `:all`             -> Consider every cell/particle within the simulation box.
      + `:halo`            -> Consider only the cells/particles that belong to the main halo.
      + `:subhalo`         -> Consider only the cells/particles that belong to the main subhalo.
      + `:sphere`          -> Consider only the cell/particle inside a sphere with radius `DISK_R` (see `./src/constants/globals.jl`).
      + `:stellar_subhalo` -> Consider only the cells/particles that belong to the main subhalo.
      + `:all_subhalo`     -> Plot every cell/particle centered around the main subhalo.
      + A dictionary with three entries:

          + `:filter_function` -> The filter function.
          + `:translation`     -> Translation for the simulation box. The posibilities are:

              + `:global_cm`                  -> Selects the center of mass of the whole system as the new origin.
              + `:{component}`                -> Sets the center of mass of the given component (e.g. :stars, :gas, :halo, etc, after filtering) as the new origin. It can be any of the keys of [`PARTICLE_INDEX`](@ref).
              + `(halo_idx, subhalo_rel_idx)` -> Sets the position of the potencial minimum for the `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo) as the new origin.
              + `(halo_idx, 0)`               -> Sets the center of mass of the `halo_idx::Int` halo as the new origin.
              + `subhalo_abs_idx`             -> Sets the center of mass of the `subhalo_abs_idx::Int` as the new origin.
          + `:rotation`        -> Rotation for the simulation box. The posibilities are:

              + `:zero`                       -> No rotation is appplied.
              + `:global_am`                  -> Sets the angular momentum of the whole system as the new z axis.
              + `:stellar_am`                 -> Sets the stellar angular momentum as the new z axis.
              + `:stellar_pa`                 -> Sets the stellar principal axis as the new coordinate system.
              + `:stellar_subhalo_pa`         -> Sets the principal axis of the stars in the main subhalo as the new coordinate system.
              + `(halo_idx, subhalo_rel_idx)` -> Sets the principal axis of the stars in `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo), as the new coordinate system.
              + `(halo_idx, 0)`               -> Sets the principal axis of the stars in the `halo_idx::Int` halo, as the new coordinate system.
              + `subhalo_abs_idx`             -> Sets the principal axis of the stars in the `subhalo_abs_idx::Int` subhalo as the new coordinate system.
"""
function stellarBirthHalos(
    simulation_path::String,
    slice_n::Int;
    output_path::String="./",
    filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all,
)::Nothing

    # Select the filter function and request dictionary
    filter_function, _, _, request = selectFilter(filter_mode, Dict(:stars=>["ID  "]))

    # Read the relevant data of the snapshot
    data_dict = makeDataDict(
        simulation_path,
        slice_n,
        request;
        warnings=false,
    )

    # Filter the data
    filterData!(data_dict; filter_function)

    # Find the birth place of every star
    birth_halo, birth_subhalo = locateStellarBirthPlace(data_dict; warnings=false)

    # Write the results to CSV files
    CSV.write(
        joinpath(output_path, "stellar_birth_halos.gz"),
        Tables.table(birth_halo);
        newline=',',
        writeheader=false,
        compress=true,
    )
    CSV.write(
        joinpath(output_path, "stellar_birth_subhalos.gz"),
        Tables.table(birth_subhalo);
        newline=',',
        writeheader=false,
        compress=true,
    )

    return nothing

end

"""
    densityMap(
        simulation_paths::Vector{String},
        slice::IndexType;
        <keyword arguments>
    )::Nothing

Plot a 2D histogram of the density.

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`.
  - `slice::IndexType`: Slice of the simulations, i.e. which snapshots will be plotted. It can be an integer (a single snapshot), a vector of integers (several snapshots), an `UnitRange` (e.g. 5:13), an `StepRange` (e.g. 5:2:13) or (:) (all snapshots). Starts at 1 and out of bounds indices are ignored. If set to 0, an animation using every snapshots will be made.
  - `quantities::Vector{Symbol}=[:gas_mass]`: Quantities for which the density will be calculated. The options are:

      + `:stellar_mass`      -> Stellar mass.
      + `:gas_mass`          -> Gas mass.
      + `:hydrogen_mass`     -> Hydrogen mass.
      + `:dm_mass`           -> Dark matter mass.
      + `:bh_mass`           -> Black hole mass.
      + `:molecular_mass`    -> Molecular hydrogen (``\\mathrm{H_2}``) mass.
      + `:br_molecular_mass` -> Molecular hydrogen (``\\mathrm{H_2}``) mass, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_mass`       -> Atomic hydrogen (``\\mathrm{HI}``) mass.
      + `:ionized_mass`      -> Ionized hydrogen (``\\mathrm{HII}``) mass.
      + `:neutral_mass`      -> Neutral hydrogen (``\\mathrm{HI + H_2}``) mass.
  - `types::Vector{Symbol}=[:cells]`: List of component types for the density fields, each element can be either `:particles` or Voronoi `:cells`.
  - `output_path::String="./"`: Path to the output folder.
  - `filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all`: Which cells/particles will be plotted, the options are:

      + `:all`             -> Consider every cell/particle within the simulation box.
      + `:halo`            -> Consider only the cells/particles that belong to the main halo.
      + `:subhalo`         -> Consider only the cells/particles that belong to the main subhalo.
      + `:sphere`          -> Consider only the cell/particle inside a sphere with radius `DISK_R` (see `./src/constants/globals.jl`).
      + `:stellar_subhalo` -> Consider only the cells/particles that belong to the main subhalo.
      + `:all_subhalo`     -> Plot every cell/particle centered around the main subhalo.
      + A dictionary with three entries:

          + `:filter_function` -> The filter function.
          + `:translation`     -> Translation for the simulation box. The posibilities are:

              + `:global_cm`                  -> Selects the center of mass of the whole system as the new origin.
              + `:{component}`                -> Sets the center of mass of the given component (e.g. :stars, :gas, :halo, etc, after filtering) as the new origin. It can be any of the keys of [`PARTICLE_INDEX`](@ref).
              + `(halo_idx, subhalo_rel_idx)` -> Sets the position of the potencial minimum for the `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo) as the new origin.
              + `(halo_idx, 0)`               -> Sets the center of mass of the `halo_idx::Int` halo as the new origin.
              + `subhalo_abs_idx`             -> Sets the center of mass of the `subhalo_abs_idx::Int` as the new origin.
          + `:rotation`        -> Rotation for the simulation box. The posibilities are:

              + `:zero`                       -> No rotation is appplied.
              + `:global_am`                  -> Sets the angular momentum of the whole system as the new z axis.
              + `:stellar_am`                 -> Sets the stellar angular momentum as the new z axis.
              + `:stellar_pa`                 -> Sets the stellar principal axis as the new coordinate system.
              + `:stellar_subhalo_pa`         -> Sets the principal axis of the stars in the main subhalo as the new coordinate system.
              + `(halo_idx, subhalo_rel_idx)` -> Sets the principal axis of the stars in `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo), as the new coordinate system.
              + `(halo_idx, 0)`               -> Sets the principal axis of the stars in the `halo_idx::Int` halo, as the new coordinate system.
              + `subhalo_abs_idx`             -> Sets the principal axis of the stars in the `subhalo_abs_idx::Int` subhalo as the new coordinate system.
  - `projection_planes::Vector{Symbol}=[:xy]`: Projection planes. The options are `:xy`, `:xz`, and `:yz`. The disk is generally oriented to have its axis of rotation parallel to the z axis.
  - `box_size::Unitful.Length=100u"kpc"`: Physical side length of the plot window.
  - `pixel_length::Unitful.Length=0.1u"kpc"`: Pixel (bin of the 2D histogram) side length.
  - `reduce::Int=1`: Factor by which the resolution of the result will be reduced. This will be applied after the density proyection, averaging the value of neighboring pixels. It has to divide the size of `grid` exactly.
  - `print_range::Bool=false`: Print an info block detailing the logarithmic density range.
  - `theme::Attributes=Theme()`: Plot theme that will take precedence over [`DEFAULT_THEME`](@ref).
  - `title::Union{Symbol,<:AbstractString}=""`: Title for the figure. If left empty, no title is printed. It can also be set to one of the following options:

      + `:physical_time` -> Physical time since the Big Bang.
      + `:lookback_time` -> Physical time left to reach the last snapshot.
      + `:scale_factor`  -> Scale factor (only relevant for cosmological simulations).
      + `:redshift`      -> Redshift (only relevant for cosmological simulations).
  - `annotation::AbstractString=""`: Text to be added into the top left corner of the plot. If left empty, nothing is printed.
  - `colorbar::Bool=false`: If a colorbar will be added.
  - `colorrange::Union{Nothing,Tuple{<:Real,<:Real}}=nothing`: Sets the start and end points of the colormap. Use `nothing` to use the extrema of the values to be plotted.
  - `da_ff::Function=filterNothing`: Filter function for the data analysis function. It must be a function with the signature:

    `da_ff(data_dict) -> indices`

    where

      + `data_dict::Dict`: A dictionary with the following shape:

        * `:sim_data`          -> ::Simulation (see [`Simulation`](@ref)).
        * `:snap_data`         -> ::Snapshot (see [`Snapshot`](@ref)).
        * `:gc_data`           -> ::GroupCatalog (see [`GroupCatalog`](@ref)).
        * `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * ...
        * `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * ...
      + `indices::Dict`: A dictionary with the following shape:

        * `cell/particle type` -> idxs::IndexType
        * `cell/particle type` -> idxs::IndexType
        * `cell/particle type` -> idxs::IndexType
        * ...
  - `ff_request::Dict{Symbol,Vector{String}}=Dict{Symbol,Vector{String}}()`: Request dictionary for the `da_ff` filter function.
"""
function densityMap(
    simulation_paths::Vector{String},
    slice::IndexType;
    quantities::Vector{Symbol}=[:gas_mass],
    types::Vector{Symbol}=[:cells],
    output_path::String="./",
    filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all,
    projection_planes::Vector{Symbol}=[:xy],
    box_size::Unitful.Length=100u"kpc",
    pixel_length::Unitful.Length=0.1u"kpc",
    reduce::Int=1,
    print_range::Bool=false,
    theme::Attributes=Theme(),
    title::Union{Symbol,<:AbstractString}="",
    annotation::AbstractString="",
    colorbar::Bool=false,
    colorrange::Union{Nothing,Tuple{<:Real,<:Real}}=nothing,
    da_ff::Function=filterNothing,
    ff_request::Dict{Symbol,Vector{String}}=Dict{Symbol,Vector{String}}(),
)::Nothing

    # Compute the axes limits, to avoid white padding around the heatmap grid
    limit = ustrip(u"kpc", box_size / 2.0)

    # Compute number of pixel per side
    resolution = round(Int, box_size / pixel_length)

    # Set up the grid
    grid = CubicGrid(box_size, resolution)

    pf_kwargs = isnothing(colorrange) ? [(;)] : [(; colorrange)]

    @inbounds for (i, quantity) in pairs(quantities)

        filter_function, translation, rotation, request = selectFilter(
            filter_mode,
            mergeRequests(plotParams(quantity).request, ff_request),
        )

        @inbounds for simulation_path in simulation_paths

            # Get the simulation name as a string
            sim_name = basename(simulation_path)

            @inbounds for projection_plane in projection_planes

                # Construct the file name
                base_filename = "$(sim_name)-$(quantity)-$(projection_plane)-density_map"

                plotSnapshot(
                    [simulation_path],
                    request,
                    [heatmap!];
                    pf_kwargs,
                    # `plotSnapshot` configuration
                    output_path,
                    base_filename,
                    output_format=".png",
                    warnings=false,
                    show_progress=true,
                    # Data manipulation options
                    slice=iszero(slice) ? (:) : slice,
                    filter_function,
                    da_functions=[daDensity2DProjection],
                    da_args=[(grid, quantity, ring(types, i))],
                    da_kwargs=[
                        (;
                            reduce,
                            projection_plane,
                            print_range,
                            filter_function=da_ff,
                        ),
                    ],
                    post_processing=isempty(annotation) ? getNothing : ppAnnotation!,
                    pp_args=(annotation,),
                    pp_kwargs=(; color=:white),
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
                    # Plotting options
                    save_figures=!iszero(slice),
                    backup_results=iszero(slice),
                    theme=merge(
                        theme,
                        Theme(
                            size=colorbar ? (880, 760) : (880, 880),
                            figure_padding=(5, 20, 20, 10),
                            Axis=(limits=(-limit, limit, -limit, limit),),
                        ),
                    ),
                    sim_labels=nothing,
                    title,
                    colorbar,
                    # Animation options
                    animation=iszero(slice),
                    animation_filename="$(base_filename).mp4",
                    framerate=5,
                )

            end

        end

    end

    return nothing

end

"""
    gasSFRMap(
        simulation_paths::Vector{String},
        slice::IndexType;
        <keyword arguments>
    )::Nothing

Plot a 2D map of the gas SFR.

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`.
  - `slice::IndexType`: Slice of the simulations, i.e. which snapshots will be plotted. It can be an integer (a single snapshot), a vector of integers (several snapshots), an `UnitRange` (e.g. 5:13), an `StepRange` (e.g. 5:2:13) or (:) (all snapshots). Starts at 1 and out of bounds indices are ignored. If set to 0, an animation using every snapshots will be made.
  - `types::Symbol=:cells`: Gas type for the SFR fields. It can be either `:particles` or Voronoi `:cells`.
  - `output_path::String="./"`: Path to the output folder.
  - `filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all`: Which cells/particles will be plotted, the options are:

      + `:all`             -> Consider every cell/particle within the simulation box.
      + `:halo`            -> Consider only the cells/particles that belong to the main halo.
      + `:subhalo`         -> Consider only the cells/particles that belong to the main subhalo.
      + `:sphere`          -> Consider only the cell/particle inside a sphere with radius `DISK_R` (see `./src/constants/globals.jl`).
      + `:stellar_subhalo` -> Consider only the cells/particles that belong to the main subhalo.
      + `:all_subhalo`     -> Plot every cell/particle centered around the main subhalo.
      + A dictionary with three entries:

          + `:filter_function` -> The filter function.
          + `:translation`     -> Translation for the simulation box. The posibilities are:

              + `:global_cm`                  -> Selects the center of mass of the whole system as the new origin.
              + `:{component}`                -> Sets the center of mass of the given component (e.g. :stars, :gas, :halo, etc, after filtering) as the new origin. It can be any of the keys of [`PARTICLE_INDEX`](@ref).
              + `(halo_idx, subhalo_rel_idx)` -> Sets the position of the potencial minimum for the `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo) as the new origin.
              + `(halo_idx, 0)`               -> Sets the center of mass of the `halo_idx::Int` halo as the new origin.
              + `subhalo_abs_idx`             -> Sets the center of mass of the `subhalo_abs_idx::Int` as the new origin.
          + `:rotation`        -> Rotation for the simulation box. The posibilities are:

              + `:zero`                       -> No rotation is appplied.
              + `:global_am`                  -> Sets the angular momentum of the whole system as the new z axis.
              + `:stellar_am`                 -> Sets the stellar angular momentum as the new z axis.
              + `:stellar_pa`                 -> Sets the stellar principal axis as the new coordinate system.
              + `:stellar_subhalo_pa`         -> Sets the principal axis of the stars in the main subhalo as the new coordinate system.
              + `(halo_idx, subhalo_rel_idx)` -> Sets the principal axis of the stars in `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo), as the new coordinate system.
              + `(halo_idx, 0)`               -> Sets the principal axis of the stars in the `halo_idx::Int` halo, as the new coordinate system.
              + `subhalo_abs_idx`             -> Sets the principal axis of the stars in the `subhalo_abs_idx::Int` subhalo as the new coordinate system.
  - `projection_planes::Vector{Symbol}=[:xy]`: Projection planes. The options are `:xy`, `:xz`, and `:yz`. The disk is generally oriented to have its axis of rotation parallel to the z axis.
  - `box_size::Unitful.Length=100u"kpc"`: Physical side length of the plot window.
  - `pixel_length::Unitful.Length=0.1u"kpc"`: Pixel (bin of the 2D histogram) side length.
  - `print_range::Bool=false`: Print an info block detailing the logarithmic density range.
  - `theme::Attributes=Theme()`: Plot theme that will take precedence over [`DEFAULT_THEME`](@ref).
  - `title::Union{Symbol,<:AbstractString}=""`: Title for the figure. If left empty, no title is printed. It can also be set to one of the following options:

      + `:physical_time` -> Physical time since the Big Bang.
      + `:lookback_time` -> Physical time left to reach the last snapshot.
      + `:scale_factor`  -> Scale factor (only relevant for cosmological simulations).
      + `:redshift`      -> Redshift (only relevant for cosmological simulations).
  - `annotation::AbstractString=""`: Text to be added into the top left corner of the plot. If left empty, nothing is printed.
  - `colorbar::Bool=false`: If a colorbar will be added.
  - `colorrange::Union{Nothing,Tuple{<:Real,<:Real}}=nothing`: Sets the start and end points of the colormap. Use `nothing` to use the extrema of the values to be plotted.
  - `da_ff::Function=filterNothing`: Filter function for the data analysis function. It must be a function with the signature:

    `da_ff(data_dict) -> indices`

    where

      + `data_dict::Dict`: A dictionary with the following shape:

        * `:sim_data`          -> ::Simulation (see [`Simulation`](@ref)).
        * `:snap_data`         -> ::Snapshot (see [`Snapshot`](@ref)).
        * `:gc_data`           -> ::GroupCatalog (see [`GroupCatalog`](@ref)).
        * `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * ...
        * `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * ...
      + `indices::Dict`: A dictionary with the following shape:

        * `cell/particle type` -> idxs::IndexType
        * `cell/particle type` -> idxs::IndexType
        * `cell/particle type` -> idxs::IndexType
        * ...
  - `ff_request::Dict{Symbol,Vector{String}}=Dict{Symbol,Vector{String}}()`: Request dictionary for the `da_ff` filter function.
"""
function gasSFRMap(
    simulation_paths::Vector{String},
    slice::IndexType;
    type::Symbol=:cells,
    output_path::String="./",
    filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all,
    projection_planes::Vector{Symbol}=[:xy],
    box_size::Unitful.Length=100u"kpc",
    pixel_length::Unitful.Length=0.1u"kpc",
    print_range::Bool=false,
    theme::Attributes=Theme(),
    title::Union{Symbol,<:AbstractString}="",
    annotation::AbstractString="",
    colorbar::Bool=false,
    colorrange::Union{Nothing,Tuple{<:Real,<:Real}}=nothing,
    da_ff::Function=filterNothing,
    ff_request::Dict{Symbol,Vector{String}}=Dict{Symbol,Vector{String}}(),
)::Nothing

    # Compute the axes limits, to avoid white padding around the heatmap grid
    limit = ustrip(u"kpc", box_size / 2.0)

    # Compute number of pixel per side
    resolution = round(Int, box_size / pixel_length)

    # Set up the grid
    grid = CubicGrid(box_size, resolution)

    pf_kwargs = isnothing(colorrange) ? [(;)] : [(; colorrange)]

    filter_function, translation, rotation, request = selectFilter(
        filter_mode,
        mergeRequests(
            plotParams(:gas_sfr).request,
            plotParams(:gas_mass_density).request,
            ff_request,
        ),
    )

    @inbounds for simulation_path in simulation_paths

        # Get the simulation name as a string
        sim_name = basename(simulation_path)

        @inbounds for projection_plane in projection_planes

            # Construct the file name
            base_filename = "$(sim_name)-$(projection_plane)-gas_sfr_map"

            plotSnapshot(
                [simulation_path],
                request,
                [heatmap!];
                pf_kwargs,
                # `plotSnapshot` configuration
                output_path,
                base_filename,
                output_format=".png",
                warnings=false,
                show_progress=true,
                # Data manipulation options
                slice=iszero(slice) ? (:) : slice,
                filter_function,
                da_functions=[daGasSFR2DProjection],
                da_args=[(grid, type)],
                da_kwargs=[
                    (;
                        projection_plane,
                        print_range,
                        filter_function=da_ff,
                    ),
                ],
                post_processing=isempty(annotation) ? getNothing : ppAnnotation!,
                pp_args=(annotation,),
                pp_kwargs=(; color=:white),
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
                # Plotting options
                save_figures=!iszero(slice),
                backup_results=iszero(slice),
                theme=merge(
                    theme,
                    Theme(
                        size=colorbar ? (880, 760) : (880, 880),
                        figure_padding=(5, 20, 20, 10),
                        Axis=(limits=(-limit, limit, -limit, limit),),
                    ),
                ),
                sim_labels=nothing,
                title,
                colorbar,
                # Animation options
                animation=iszero(slice),
                animation_filename="$(base_filename).mp4",
                framerate=5,
            )

        end

    end

    return nothing

end

"""
    densityMapVelField(
        simulation_paths::Vector{String},
        slice::IndexType;
        <keyword arguments>
    )::Nothing

Plot a 2D histogram of the density, with the velocity field.

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`.
  - `slice::IndexType`: Slice of the simulations, i.e. which snapshots will be plotted. It can be an integer (a single snapshot), a vector of integers (several snapshots), an `UnitRange` (e.g. 5:13), an `StepRange` (e.g. 5:2:13) or (:) (all snapshots). Starts at 1 and out of bounds indices are ignored. If set to 0, an animation using every snapshots will be made.
  - `quantities::Vector{Symbol}=[:gas_mass]`: Quantities for which the density will be calculated. The options are:

      + `:stellar_mass`      -> Stellar mass.
      + `:gas_mass`          -> Gas mass.
      + `:hydrogen_mass`     -> Hydrogen mass.
      + `:dm_mass`           -> Dark matter mass.
      + `:bh_mass`           -> Black hole mass.
      + `:molecular_mass`    -> Molecular hydrogen (``\\mathrm{H_2}``) mass.
      + `:br_molecular_mass` -> Molecular hydrogen (``\\mathrm{H_2}``) mass, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_mass`       -> Atomic hydrogen (``\\mathrm{HI}``) mass.
      + `:ionized_mass`      -> Ionized hydrogen (``\\mathrm{HII}``) mass.
      + `:neutral_mass`      -> Neutral hydrogen (``\\mathrm{HI + H_2}``) mass.
  - `types::Vector{Symbol}=[:cells]`: List of component types for the density fields, each element can be either `:particles` or Voronoi `:cells`.
  - `output_path::String="./"`: Path to the output folder.
  - `filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all`: Which cells/particles will be plotted, the options are:

      + `:all`             -> Consider every cell/particle within the simulation box.
      + `:halo`            -> Consider only the cells/particles that belong to the main halo.
      + `:subhalo`         -> Consider only the cells/particles that belong to the main subhalo.
      + `:sphere`          -> Consider only the cell/particle inside a sphere with radius `DISK_R` (see `./src/constants/globals.jl`).
      + `:stellar_subhalo` -> Consider only the cells/particles that belong to the main subhalo.
      + `:all_subhalo`     -> Plot every cell/particle centered around the main subhalo.
      + A dictionary with three entries:

          + `:filter_function` -> The filter function.
          + `:translation`     -> Translation for the simulation box. The posibilities are:

              + `:global_cm`                  -> Selects the center of mass of the whole system as the new origin.
              + `:{component}`                -> Sets the center of mass of the given component (e.g. :stars, :gas, :halo, etc, after filtering) as the new origin. It can be any of the keys of [`PARTICLE_INDEX`](@ref).
              + `(halo_idx, subhalo_rel_idx)` -> Sets the position of the potencial minimum for the `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo) as the new origin.
              + `(halo_idx, 0)`               -> Sets the center of mass of the `halo_idx::Int` halo as the new origin.
              + `subhalo_abs_idx`             -> Sets the center of mass of the `subhalo_abs_idx::Int` as the new origin.
          + `:rotation`        -> Rotation for the simulation box. The posibilities are:

              + `:zero`                       -> No rotation is appplied.
              + `:global_am`                  -> Sets the angular momentum of the whole system as the new z axis.
              + `:stellar_am`                 -> Sets the stellar angular momentum as the new z axis.
              + `:stellar_pa`                 -> Sets the stellar principal axis as the new coordinate system.
              + `:stellar_subhalo_pa`         -> Sets the principal axis of the stars in the main subhalo as the new coordinate system.
              + `(halo_idx, subhalo_rel_idx)` -> Sets the principal axis of the stars in `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo), as the new coordinate system.
              + `(halo_idx, 0)`               -> Sets the principal axis of the stars in the `halo_idx::Int` halo, as the new coordinate system.
              + `subhalo_abs_idx`             -> Sets the principal axis of the stars in the `subhalo_abs_idx::Int` subhalo as the new coordinate system.
  - `projection_planes::Vector{Symbol}=[:xy]`: Projection planes. The options are `:xy`, `:xz`, and `:yz`. The disk is generally oriented to have its axis of rotation parallel to the z axis.
  - `box_size::Unitful.Length=100u"kpc"`: Physical side length of the plot window.
  - `pixel_length::Unitful.Length=0.1u"kpc"`: Pixel (bin of the 2D histogram) side length.
  - `print_range::Bool=false`: Print an info block detailing the logarithmic density range.
  - `theme::Attributes=Theme()`: Plot theme that will take precedence over [`DEFAULT_THEME`](@ref).
  - `title::Union{Symbol,<:AbstractString}=""`: Title for the figure. If left empty, no title is printed. It can also be set to one of the following options:

      + `:physical_time` -> Physical time since the Big Bang.
      + `:lookback_time` -> Physical time left to reach the last snapshot.
      + `:scale_factor`  -> Scale factor (only relevant for cosmological simulations).
      + `:redshift`      -> Redshift (only relevant for cosmological simulations).
  - `annotation::AbstractString=""`: Text to be added into the top left corner of the plot. If left empty, nothing is printed.
  - `colorbar::Bool=false`: If a colorbar will be added.
  - `colorrange::Union{Nothing,Tuple{<:Real,<:Real}}=nothing`: Sets the start and end points of the colormap. Use `nothing` to use the extrema of the values to be plotted.
"""
function densityMapVelField(
    simulation_paths::Vector{String},
    slice::IndexType;
    quantities::Vector{Symbol}=[:gas_mass],
    types::Vector{Symbol}=[:cells],
    output_path::String="./",
    filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all,
    projection_planes::Vector{Symbol}=[:xy],
    box_size::Unitful.Length=100u"kpc",
    pixel_length::Unitful.Length=0.1u"kpc",
    print_range::Bool=false,
    theme::Attributes=Theme(),
    title::Union{Symbol,<:AbstractString}="",
    annotation::AbstractString="",
    colorbar::Bool=false,
    colorrange::Union{Nothing,Tuple{<:Real,<:Real}}=nothing,
)::Nothing

    # Compute the axes limits, to avoid white padding around the heatmap grid
    limit = ustrip(u"kpc", box_size / 2.0)

    # Compute number of pixel per side
    resolution = round(Int, box_size / pixel_length)

    # Set up the grid for the heatmap
    grid_hm = CubicGrid(box_size, resolution)

    # Set up the grid for the velocity field
    grid_vf = SquareGrid(box_size, 25)

    pf_kwargs = isnothing(colorrange) ? [(;), (;)] : [(; colorrange), (;)]

    @inbounds for (i, quantity) in pairs(quantities)

        filter_function, translation, rotation, request = selectFilter(
            filter_mode,
            plotParams(quantity).request,
        )

        if quantity ∈ [
            :gas_mass,
            :hydrogen_mass,
            :molecular_mass,
            :br_molecular_mass,
            :atomic_mass,
            :ionized_mass,
            :neutral_mass,
        ]
            component = :gas
        elseif quantity == :stellar_mass
            component = :stars
        elseif quantity == :dm_mass
            component = :halo
        elseif quantity == :bh_mass
            component = :black_hole
        else
            throw(ArgumentError("densityMapVelField: `quantities` contains :$(quantity), \
            which is not a valid symbol. See the documentation for valid options."))
        end

        @inbounds for simulation_path in simulation_paths

            # Get the simulation name as a string
            sim_name = basename(simulation_path)

            @inbounds for projection_plane in projection_planes

                # Construct the file name
                base_filename = "$(sim_name)-$(quantity)-$(projection_plane)-density_map"

                plotSnapshot(
                    [simulation_path, simulation_path],
                    request,
                    [heatmap!, arrows!];
                    pf_kwargs,
                    # `plotSnapshot` configuration
                    output_path,
                    base_filename,
                    output_format=".png",
                    warnings=false,
                    show_progress=true,
                    # Data manipulation options
                    slice=iszero(slice) ? (:) : slice,
                    filter_function,
                    da_functions=[daDensity2DProjection, daVelocityField],
                    da_args=[(grid_hm, quantity, ring(types, i)), (grid_vf, component)],
                    da_kwargs=[
                        (; projection_plane, print_range),
                        (; projection_plane),
                    ],
                    post_processing=isempty(annotation) ? getNothing : ppAnnotation!,
                    pp_args=(annotation,),
                    pp_kwargs=(; color=:white),
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
                    # Plotting options
                    save_figures=!iszero(slice),
                    backup_results=iszero(slice),
                    theme=merge(
                        theme,
                        Theme(
                            size=colorbar ? (880, 680) : (880, 880),
                            figure_padding=(1, 50, 1, 1),
                            Axis=(limits=(-limit, limit, -limit, limit),),
                            Colorbar=(
                                label=L"\mathrm{log}_{10} \Sigma \,\, [\mathrm{M_\odot \, kpc^{-2}}]",
                            ),
                        ),
                    ),
                    sim_labels=nothing,
                    title,
                    colorbar,
                    # Animation options
                    animation=iszero(slice),
                    animation_filename="$(base_filename).mp4",
                    framerate=5,
                )

            end

        end

    end

    return nothing

end

"""
    metallicityMap(
        simulation_paths::Vector{String},
        slice::IndexType;
        <keyword arguments>
    )::Nothing

Plot a 2D histogram of the metallicity.

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`.
  - `slice::IndexType`: Slice of the simulations, i.e. which snapshots will be plotted. It can be an integer (a single snapshot), a vector of integers (several snapshots), an `UnitRange` (e.g. 5:13), an `StepRange` (e.g. 5:2:13) or (:) (all snapshots). Starts at 1 and out of bounds indices are ignored. If set to 0, an animation using every snapshots will be made.
  - `components::Vector{Symbol}=[:gas]`: Target component. It can be either `:stars` or `:gas`.
  - `types::Vector{Symbol}=[:cells]`: List of component types for the metallicity fields, each element can be either `:particles` or Voronoi `:cells`.
  - `output_path::String="./"`: Path to the output folder.
  - `filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all`: Which cells/particles will be plotted, the options are:

      + `:all`             -> Consider every cell/particle within the simulation box.
      + `:halo`            -> Consider only the cells/particles that belong to the main halo.
      + `:subhalo`         -> Consider only the cells/particles that belong to the main subhalo.
      + `:sphere`          -> Consider only the cell/particle inside a sphere with radius `DISK_R` (see `./src/constants/globals.jl`).
      + `:stellar_subhalo` -> Consider only the cells/particles that belong to the main subhalo.
      + `:all_subhalo`     -> Plot every cell/particle centered around the main subhalo.
      + A dictionary with three entries:

          + `:filter_function` -> The filter function.
          + `:translation`     -> Translation for the simulation box. The posibilities are:

              + `:global_cm`                  -> Selects the center of mass of the whole system as the new origin.
              + `:{component}`                -> Sets the center of mass of the given component (e.g. :stars, :gas, :halo, etc, after filtering) as the new origin. It can be any of the keys of [`PARTICLE_INDEX`](@ref).
              + `(halo_idx, subhalo_rel_idx)` -> Sets the position of the potencial minimum for the `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo) as the new origin.
              + `(halo_idx, 0)`               -> Sets the center of mass of the `halo_idx::Int` halo as the new origin.
              + `subhalo_abs_idx`             -> Sets the center of mass of the `subhalo_abs_idx::Int` as the new origin.
          + `:rotation`        -> Rotation for the simulation box. The posibilities are:

              + `:zero`                       -> No rotation is appplied.
              + `:global_am`                  -> Sets the angular momentum of the whole system as the new z axis.
              + `:stellar_am`                 -> Sets the stellar angular momentum as the new z axis.
              + `:stellar_pa`                 -> Sets the stellar principal axis as the new coordinate system.
              + `:stellar_subhalo_pa`         -> Sets the principal axis of the stars in the main subhalo as the new coordinate system.
              + `(halo_idx, subhalo_rel_idx)` -> Sets the principal axis of the stars in `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo), as the new coordinate system.
              + `(halo_idx, 0)`               -> Sets the principal axis of the stars in the `halo_idx::Int` halo, as the new coordinate system.
              + `subhalo_abs_idx`             -> Sets the principal axis of the stars in the `subhalo_abs_idx::Int` subhalo as the new coordinate system.
  - `projection_planes::Vector{Symbol}=[:xy]`: Projection planes. The options are `:xy`, `:xz`, and `:yz`. The disk is generally oriented to have its axis of rotation parallel to the z axis.
  - `box_size::Unitful.Length=100u"kpc"`: Physical side length of the plot window.
  - `pixel_length::Unitful.Length=0.1u"kpc"`: Pixel (bin of the 2D histogram) side length.
  - `print_range::Bool=false`: Print an info block detailing the logarithmic metallicity range.
  - `theme::Attributes=Theme()`: Plot theme that will take precedence over [`DEFAULT_THEME`](@ref).
  - `title::Union{Symbol,<:AbstractString}=""`: Title for the figure. If left empty, no title is printed. It can also be set to one of the following options:

      + `:physical_time` -> Physical time since the Big Bang.
      + `:lookback_time` -> Physical time left to reach the last snapshot.
      + `:scale_factor`  -> Scale factor (only relevant for cosmological simulations).
      + `:redshift`      -> Redshift (only relevant for cosmological simulations).
  - `annotation::AbstractString=""`: Text to be added into the top left corner of the plot. If left empty, nothing is printed.
  - `colorbar::Bool=false`: If a colorbar will be added.
  - `colorrange::Union{Nothing,Tuple{<:Real,<:Real}}=nothing`: Sets the start and end points of the colormap. Use `nothing` to use the extrema of the values to be plotted.
  - `da_ff::Function=filterNothing`: Filter function for the data analysis function. It must be a function with the signature:

    `da_ff(data_dict) -> indices`

    where

      + `data_dict::Dict`: A dictionary with the following shape:

        * `:sim_data`          -> ::Simulation (see [`Simulation`](@ref)).
        * `:snap_data`         -> ::Snapshot (see [`Snapshot`](@ref)).
        * `:gc_data`           -> ::GroupCatalog (see [`GroupCatalog`](@ref)).
        * `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * ...
        * `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * ...
      + `indices::Dict`: A dictionary with the following shape:

        * `cell/particle type` -> idxs::IndexType
        * `cell/particle type` -> idxs::IndexType
        * `cell/particle type` -> idxs::IndexType
        * ...
  - `ff_request::Dict{Symbol,Vector{String}}=Dict{Symbol,Vector{String}}()`: Request dictionary for the `da_ff` filter function.
"""
function metallicityMap(
    simulation_paths::Vector{String},
    slice::IndexType;
    components::Vector{Symbol}=[:gas],
    types::Vector{Symbol}=[:cells],
    output_path::String="./",
    filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all,
    projection_planes::Vector{Symbol}=[:xy],
    box_size::Unitful.Length=100u"kpc",
    pixel_length::Unitful.Length=0.1u"kpc",
    print_range::Bool=false,
    theme::Attributes=Theme(),
    title::Union{Symbol,<:AbstractString}="",
    annotation::AbstractString="",
    colorbar::Bool=false,
    colorrange::Union{Nothing,Tuple{<:Real,<:Real}}=nothing,
    da_ff::Function=filterNothing,
    ff_request::Dict{Symbol,Vector{String}}=Dict{Symbol,Vector{String}}(),
)::Nothing

    # Compute the axes limits, to avoid white padding around the heatmap grid
    limit = ustrip(u"kpc", box_size / 2.0)

    # Compute number of pixel per side
    resolution = round(Int, box_size / pixel_length)

    # Set up the grid
    grid = CubicGrid(box_size, resolution)

    pf_kwargs = isnothing(colorrange) ? [(;)] : [(; colorrange)]

    @inbounds for (i, component) in pairs(components)

        if component == :gas

            filter_function, translation, rotation, request = selectFilter(
                filter_mode,
                mergeRequests(
                    plotParams(:gas_metallicity).request,
                    plotParams(:gas_mass_density).request,
                    ff_request,
                ),
            )

        elseif component == :stars

            filter_function, translation, rotation, request = selectFilter(
                filter_mode,
                mergeRequests(plotParams(:stellar_metallicity).request, ff_request),
            )

        else

            throw(ArgumentError("metallicityMap: I don't recognize the component \
            :$(component)"))

        end

        @inbounds for simulation_path in simulation_paths

            # Get the simulation name as a string
            sim_name = basename(simulation_path)

            @inbounds for projection_plane in projection_planes

                # Construct the file name
                base_filename = "$(sim_name)-$(component)-$(projection_plane)-metallicity_map"

                plotSnapshot(
                    [simulation_path],
                    request,
                    [heatmap!];
                    pf_kwargs,
                    # `plotSnapshot` configuration
                    output_path,
                    base_filename,
                    output_format=".png",
                    warnings=false,
                    show_progress=true,
                    # Data manipulation options
                    slice=iszero(slice) ? (:) : slice,
                    filter_function,
                    da_functions=[daMetallicity2DProjection],
                    da_args=[(grid, component, ring(types, i))],
                    da_kwargs=[
                        (;
                            element=:all,
                            projection_plane,
                            print_range,
                            filter_function=da_ff,
                        ),
                    ],
                    post_processing=isempty(annotation) ? getNothing : ppAnnotation!,
                    pp_args=(annotation,),
                    pp_kwargs=(; color=:white),
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
                    # Plotting options
                    save_figures=!iszero(slice),
                    backup_results=iszero(slice),
                    theme=merge(
                        theme,
                        Theme(
                            size=colorbar ? (880, 760) : (880, 880),
                            figure_padding=(5, 20, 20, 10),
                            Axis=(limits=(-limit, limit, -limit, limit),),
                        ),
                    ),
                    sim_labels=nothing,
                    title,
                    colorbar,
                    # Animation options
                    animation=iszero(slice),
                    animation_filename="$(base_filename).mp4",
                    framerate=5,
                )

            end

        end

    end

    return nothing

end

"""
    temperatureMap(
        simulation_paths::Vector{String},
        slice::IndexType;
        <keyword arguments>
    )::Nothing

Plot a 2D histogram of the temperature.

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`.
  - `slice::IndexType`: Slice of the simulations, i.e. which snapshots will be plotted. It can be an integer (a single snapshot), a vector of integers (several snapshots), an `UnitRange` (e.g. 5:13), an `StepRange` (e.g. 5:2:13) or (:) (all snapshots). Starts at 1 and out of bounds indices are ignored. If set to 0, an animation using every snapshots will be made.
  - `type::Symbol=:cells`: If the gas will be assumed to be in `:particles` or in Voronoi `:cells`.
  - `output_path::String="./"`: Path to the output folder.
  - `filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all`: Which cells/particles will be plotted, the options are:

      + `:all`             -> Consider every cell/particle within the simulation box.
      + `:halo`            -> Consider only the cells/particles that belong to the main halo.
      + `:subhalo`         -> Consider only the cells/particles that belong to the main subhalo.
      + `:sphere`          -> Consider only the cell/particle inside a sphere with radius `DISK_R` (see `./src/constants/globals.jl`).
      + `:stellar_subhalo` -> Consider only the cells/particles that belong to the main subhalo.
      + `:all_subhalo`     -> Plot every cell/particle centered around the main subhalo.
      + A dictionary with three entries:

          + `:filter_function` -> The filter function.
          + `:translation`     -> Translation for the simulation box. The posibilities are:

              + `:global_cm`                  -> Selects the center of mass of the whole system as the new origin.
              + `:{component}`                -> Sets the center of mass of the given component (e.g. :stars, :gas, :halo, etc, after filtering) as the new origin. It can be any of the keys of [`PARTICLE_INDEX`](@ref).
              + `(halo_idx, subhalo_rel_idx)` -> Sets the position of the potencial minimum for the `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo) as the new origin.
              + `(halo_idx, 0)`               -> Sets the center of mass of the `halo_idx::Int` halo as the new origin.
              + `subhalo_abs_idx`             -> Sets the center of mass of the `subhalo_abs_idx::Int` as the new origin.
          + `:rotation`        -> Rotation for the simulation box. The posibilities are:

              + `:zero`                       -> No rotation is appplied.
              + `:global_am`                  -> Sets the angular momentum of the whole system as the new z axis.
              + `:stellar_am`                 -> Sets the stellar angular momentum as the new z axis.
              + `:stellar_pa`                 -> Sets the stellar principal axis as the new coordinate system.
              + `:stellar_subhalo_pa`         -> Sets the principal axis of the stars in the main subhalo as the new coordinate system.
              + `(halo_idx, subhalo_rel_idx)` -> Sets the principal axis of the stars in `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo), as the new coordinate system.
              + `(halo_idx, 0)`               -> Sets the principal axis of the stars in the `halo_idx::Int` halo, as the new coordinate system.
              + `subhalo_abs_idx`             -> Sets the principal axis of the stars in the `subhalo_abs_idx::Int` subhalo as the new coordinate system.
  - `projection_planes::Vector{Symbol}=[:xy]`: Projection planes. The options are `:xy`, `:xz`, and `:yz`. The disk is generally oriented to have its axis of rotation parallel to the z axis.
  - `box_size::Unitful.Length=100u"kpc"`: Physical side length of the plot window.
  - `pixel_length::Unitful.Length=0.1u"kpc"`: Pixel (bin of the 2D histogram) side length.
  - `print_range::Bool=false`: Print an info block detailing the logarithmic density range.
  - `theme::Attributes=Theme()`: Plot theme that will take precedence over [`DEFAULT_THEME`](@ref).
  - `title::Union{Symbol,<:AbstractString}=""`: Title for the figure. If left empty, no title is printed. It can also be set to one of the following options:

      + `:physical_time` -> Physical time since the Big Bang.
      + `:lookback_time` -> Physical time left to reach the last snapshot.
      + `:scale_factor`  -> Scale factor (only relevant for cosmological simulations).
      + `:redshift`      -> Redshift (only relevant for cosmological simulations).
  - `annotation::AbstractString=""`: Text to be added into the top left corner of the plot. If left empty, nothing is printed.
  - `colorbar::Bool=false`: If a colorbar will be added.
  - `colorrange::Union{Nothing,Tuple{<:Real,<:Real}}=nothing`: Sets the start and end points of the colormap. Use `nothing` to use the extrema of the values to be plotted.
"""
function temperatureMap(
    simulation_paths::Vector{String},
    slice::IndexType;
    type::Symbol=:cells,
    output_path::String="./",
    filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all,
    projection_planes::Vector{Symbol}=[:xy],
    box_size::Unitful.Length=100u"kpc",
    pixel_length::Unitful.Length=0.1u"kpc",
    print_range::Bool=false,
    theme::Attributes=Theme(),
    title::Union{Symbol,<:AbstractString}="",
    annotation::AbstractString="",
    colorbar::Bool=false,
    colorrange::Union{Nothing,Tuple{<:Real,<:Real}}=nothing,
)::Nothing

    # Compute the axes limits, to avoid white padding around the heatmap grid
    limit = ustrip(u"kpc", box_size / 2.0)

    # Compute number of pixel per side
    resolution = round(Int, box_size / pixel_length)

    # Set up the grid
    grid = CubicGrid(box_size, resolution)

    pf_kwargs = isnothing(colorrange) ? [(;)] : [(; colorrange)]

    filter_function, translation, rotation, request = selectFilter(
        filter_mode,
        mergeRequests(plotParams(:temperature).request, plotParams(:gas_mass_density).request),
    )

    @inbounds for simulation_path in simulation_paths

        # Get the simulation name as a string
        sim_name = basename(simulation_path)

        @inbounds for projection_plane in projection_planes

            # Construct the file name
            base_filename = "$(sim_name)-$(projection_plane)-temperature_map"

            plotSnapshot(
                [simulation_path],
                request,
                [heatmap!];
                pf_kwargs,
                # `plotSnapshot` configuration
                output_path,
                base_filename,
                output_format=".png",
                warnings=false,
                show_progress=true,
                # Data manipulation options
                slice=iszero(slice) ? (:) : slice,
                filter_function,
                da_functions=[daTemperature2DProjection],
                da_args=[(grid, type)],
                da_kwargs=[(; projection_plane, print_range)],
                post_processing=isempty(annotation) ? getNothing : ppAnnotation!,
                pp_args=(annotation,),
                pp_kwargs=(; color=:white),
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
                # Plotting options
                save_figures=!iszero(slice),
                backup_results=iszero(slice),
                theme=merge(
                    theme,
                    Theme(
                        size=colorbar ? (880, 760) : (880, 880),
                        figure_padding=(5, 20, 20, 10),
                        Axis=(limits=(-limit, limit, -limit, limit),),
                    ),
                ),
                sim_labels=nothing,
                title,
                colorbar,
                # Animation options
                animation=iszero(slice),
                animation_filename="$(base_filename).mp4",
                framerate=5,
            )

        end

    end

    return nothing

end

"""
    scatterPlot(
        simulation_paths::Vector{String},
        slice::IndexType,
        x_quantity::Symbol,
        y_quantity::Symbol;
        <keyword arguments>
    )::Nothing

Plot two quantities as a scatter plot, one marker for every cell/particle.

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`.
  - `slice::IndexType`: Slice of the simulations, i.e. which snapshots will be plotted. It can be an integer (a single snapshot), a vector of integers (several snapshots), an `UnitRange` (e.g. 5:13), an `StepRange` (e.g. 5:2:13) or (:) (all snapshots). Starts at 1 and out of bounds indices are ignored.
  - `x_quantity::Symbol`: Quantity for the x axis. The options are:

      + `:stellar_mass`                -> Stellar mass.
      + `:gas_mass`                    -> Gas mass.
      + `:hydrogen_mass`               -> Hydrogen mass.
      + `:dm_mass`                     -> Dark matter mass.
      + `:bh_mass`                     -> Black hole mass.
      + `:molecular_mass`              -> Molecular hydrogen (``\\mathrm{H_2}``) mass.
      + `:br_molecular_mass`           -> Molecular hydrogen (``\\mathrm{H_2}``) mass, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_mass`                 -> Atomic hydrogen (``\\mathrm{HI}``) mass.
      + `:ionized_mass`                -> Ionized hydrogen (``\\mathrm{HII}``) mass.
      + `:neutral_mass`                -> Neutral hydrogen (``\\mathrm{HI + H_2}``) mass.
      + `:molecular_fraction`          -> Gas mass fraction of molecular hydrogen.
      + `:br_molecular_fraction`       -> Gas mass fraction of molecular hydrogen, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_fraction`             -> Gas mass fraction of atomic hydrogen.
      + `:ionized_fraction`            -> Gas mass fraction of ionized hydrogen.
      + `:neutral_fraction`            -> Gas mass fraction of neutral hydrogen.
      + `:molecular_neutral_fraction`  -> Fraction of molecular hydrogen in the neutral gas.
      + `:mol_eq_quotient`             -> Equilibrium quotient for the molecular fraction equation of the SF model.
      + `:ion_eq_quotient`             -> Equilibrium quotient for the ionized fraction equation of the SF model.
      + `:gas_mass_density`            -> Gas mass density.
      + `:hydrogen_mass_density`       -> Hydrogen mass density.
      + `:gas_number_density`          -> Gas number density.
      + `:molecular_number_density`    -> Molecular hydrogen number density.
      + `:br_molecular_number_density` -> Molecular hydrogen number density, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_number_density`       -> Atomic hydrogen number density.
      + `:ionized_number_density`      -> Ionized hydrogen number density.
      + `:neutral_number_density`      -> Neutral hydrogen number density.
      + `:gas_metallicity`             -> Mass fraction of all elements above He in the gas (solar units).
      + `:stellar_metallicity`         -> Mass fraction of all elements above He in the stars (solar units).
      + `:X_gas_abundance`             -> Gas abundance of element ``\\mathrm{X}``, as ``12 + \\log_{10}(\\mathrm{X \\, / \\, H})``. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
      + `:X_stellar_abundance`         -> Stellar abundance of element ``\\mathrm{X}``, as ``12 + \\log_{10}(\\mathrm{X \\, / \\, H})``. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
      + `:stellar_radial_distance`     -> Distance of every stellar particle to the origin.
      + `:gas_radial_distance`         -> Distance of every gas cell to the origin.
      + `:dm_radial_distance`          -> Distance of every dark matter particle to the origin.
      + `:stellar_xy_distance`         -> Projected distance of every stellar particle to the origin.
      + `:gas_xy_distance`             -> Projected distance of every gas cell to the origin.
      + `:dm_xy_distance`              -> Projected distance of every dark matter particle to the origin.
      + `:gas_sfr`                     -> SFR associated to each gas particle/cell within the code.
      + `:stellar_circularity`         -> Stellar circularity.
      + `:stellar_vcirc`               -> Stellar circular velocity.
      + `:stellar_vradial`             -> Stellar radial speed.
      + `:stellar_vtangential`         -> Stellar tangential speed.
      + `:stellar_vzstar`              -> Stellar speed in the z direction, computed as ``v_z \\, \\mathrm{sign}(z)``.
      + `:stellar_age`                 -> Stellar age.
      + `:sfr`                         -> The star formation rate.
      + `:ssfr`                        -> The specific star formation rate.
      + `:observational_sfr`           -> The star formation rate of the last `AGE_RESOLUTION`.
      + `:observational_ssfr`          -> The specific star formation rate of the last `AGE_RESOLUTION`.
      + `:temperature`                 -> Gas temperature, as ``\\log_{10}(T \\, / \\, \\mathrm{K})``.
      + `:pressure`                    -> Gas pressure.
  - `y_quantity::Symbol`: Quantity for the y axis. The options are:

      + `:stellar_mass`                -> Stellar mass.
      + `:gas_mass`                    -> Gas mass.
      + `:hydrogen_mass`               -> Hydrogen mass.
      + `:dm_mass`                     -> Dark matter mass.
      + `:bh_mass`                     -> Black hole mass.
      + `:molecular_mass`              -> Molecular hydrogen (``\\mathrm{H_2}``) mass.
      + `:br_molecular_mass`           -> Molecular hydrogen (``\\mathrm{H_2}``) mass, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_mass`                 -> Atomic hydrogen (``\\mathrm{HI}``) mass.
      + `:ionized_mass`                -> Ionized hydrogen (``\\mathrm{HII}``) mass.
      + `:neutral_mass`                -> Neutral hydrogen (``\\mathrm{HI + H_2}``) mass.
      + `:molecular_fraction`          -> Gas mass fraction of molecular hydrogen.
      + `:br_molecular_fraction`       -> Gas mass fraction of molecular hydrogen, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_fraction`             -> Gas mass fraction of atomic hydrogen.
      + `:ionized_fraction`            -> Gas mass fraction of ionized hydrogen.
      + `:neutral_fraction`            -> Gas mass fraction of neutral hydrogen.
      + `:molecular_neutral_fraction`  -> Fraction of molecular hydrogen in the neutral gas.
      + `:mol_eq_quotient`             -> Equilibrium quotient for the molecular fraction equation of the SF model.
      + `:ion_eq_quotient`             -> Equilibrium quotient for the ionized fraction equation of the SF model.
      + `:gas_mass_density`            -> Gas mass density.
      + `:hydrogen_mass_density`       -> Hydrogen mass density.
      + `:gas_number_density`          -> Gas number density.
      + `:molecular_number_density`    -> Molecular hydrogen number density.
      + `:br_molecular_number_density` -> Molecular hydrogen number density, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_number_density`       -> Atomic hydrogen number density.
      + `:ionized_number_density`      -> Ionized hydrogen number density.
      + `:neutral_number_density`      -> Neutral hydrogen number density.
      + `:gas_metallicity`             -> Mass fraction of all elements above He in the gas (solar units).
      + `:stellar_metallicity`         -> Mass fraction of all elements above He in the stars (solar units).
      + `:X_gas_abundance`             -> Gas abundance of element ``\\mathrm{X}``, as ``12 + \\log_{10}(\\mathrm{X \\, / \\, H})``. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
      + `:X_stellar_abundance`         -> Stellar abundance of element ``\\mathrm{X}``, as ``12 + \\log_{10}(\\mathrm{X \\, / \\, H})``. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
      + `:stellar_radial_distance`     -> Distance of every stellar particle to the origin.
      + `:gas_radial_distance`         -> Distance of every gas cell to the origin.
      + `:dm_radial_distance`          -> Distance of every dark matter particle to the origin.
      + `:stellar_xy_distance`         -> Projected distance of every stellar particle to the origin.
      + `:gas_xy_distance`             -> Projected distance of every gas cell to the origin.
      + `:dm_xy_distance`              -> Projected distance of every dark matter particle to the origin.
      + `:gas_sfr`                     -> SFR associated to each gas particle/cell within the code.
      + `:stellar_circularity`         -> Stellar circularity.
      + `:stellar_vcirc`               -> Stellar circular velocity.
      + `:stellar_vradial`             -> Stellar radial speed.
      + `:stellar_vtangential`         -> Stellar tangential speed.
      + `:stellar_vzstar`              -> Stellar speed in the z direction, computed as ``v_z \\, \\mathrm{sign}(z)``.
      + `:stellar_age`                 -> Stellar age.
      + `:sfr`                         -> The star formation rate.
      + `:ssfr`                        -> The specific star formation rate.
      + `:observational_sfr`           -> The star formation rate of the last `AGE_RESOLUTION`.
      + `:observational_ssfr`          -> The specific star formation rate of the last `AGE_RESOLUTION`.
      + `:temperature`                 -> Gas temperature, as ``\\log_{10}(T \\, / \\, \\mathrm{K})``.
      + `:pressure`                    -> Gas pressure.
  - `xlog::Bool=false`: If true, sets everything so the x axis is log10(`x_quantity`).
  - `ylog::Bool=false`: If true, sets everything so the y axis is log10(`y_quantity`).
  - `output_path::String="./"`: Path to the output folder.
  - `filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all`: Which cells/particles will be plotted, the options are:

      + `:all`             -> Consider every cell/particle within the simulation box.
      + `:halo`            -> Consider only the cells/particles that belong to the main halo.
      + `:subhalo`         -> Consider only the cells/particles that belong to the main subhalo.
      + `:sphere`          -> Consider only the cell/particle inside a sphere with radius `DISK_R` (see `./src/constants/globals.jl`).
      + `:stellar_subhalo` -> Consider only the cells/particles that belong to the main subhalo.
      + `:all_subhalo`     -> Plot every cell/particle centered around the main subhalo.
      + A dictionary with three entries:

          + `:filter_function` -> The filter function.
          + `:translation`     -> Translation for the simulation box. The posibilities are:

              + `:global_cm`                  -> Selects the center of mass of the whole system as the new origin.
              + `:{component}`                -> Sets the center of mass of the given component (e.g. :stars, :gas, :halo, etc, after filtering) as the new origin. It can be any of the keys of [`PARTICLE_INDEX`](@ref).
              + `(halo_idx, subhalo_rel_idx)` -> Sets the position of the potencial minimum for the `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo) as the new origin.
              + `(halo_idx, 0)`               -> Sets the center of mass of the `halo_idx::Int` halo as the new origin.
              + `subhalo_abs_idx`             -> Sets the center of mass of the `subhalo_abs_idx::Int` as the new origin.
          + `:rotation`        -> Rotation for the simulation box. The posibilities are:

              + `:zero`                       -> No rotation is appplied.
              + `:global_am`                  -> Sets the angular momentum of the whole system as the new z axis.
              + `:stellar_am`                 -> Sets the stellar angular momentum as the new z axis.
              + `:stellar_pa`                 -> Sets the stellar principal axis as the new coordinate system.
              + `:stellar_subhalo_pa`         -> Sets the principal axis of the stars in the main subhalo as the new coordinate system.
              + `(halo_idx, subhalo_rel_idx)` -> Sets the principal axis of the stars in `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo), as the new coordinate system.
              + `(halo_idx, 0)`               -> Sets the principal axis of the stars in the `halo_idx::Int` halo, as the new coordinate system.
              + `subhalo_abs_idx`             -> Sets the principal axis of the stars in the `subhalo_abs_idx::Int` subhalo as the new coordinate system.
  - `da_ff::Function=filterNothing`: A function with the signature:

    `da_ff(data_dict) -> indices`

    where

      + `data_dict::Dict`: A dictionary with the following shape:

        * `:sim_data`          -> ::Simulation (see [`Simulation`](@ref)).
        * `:snap_data`         -> ::Snapshot (see [`Snapshot`](@ref)).
        * `:gc_data`           -> ::GroupCatalog (see [`GroupCatalog`](@ref)).
        * `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * ...
        * `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * ...
      + `indices::Dict`: A dictionary with the following shape:

        * `cell/particle type` -> idxs::IndexType
        * `cell/particle type` -> idxs::IndexType
        * `cell/particle type` -> idxs::IndexType
        * ...
  - `ff_request::Dict{Symbol,Vector{String}}=Dict{Symbol,Vector{String}}()`: Request dictionary for the `da_ff` filter function.
  - `theme::Attributes=Theme()`: Plot theme that will take precedence over [`DEFAULT_THEME`](@ref).
"""
function scatterPlot(
    simulation_paths::Vector{String},
    slice::IndexType,
    x_quantity::Symbol,
    y_quantity::Symbol;
    xlog::Bool=false,
    ylog::Bool=false,
    output_path::String="./",
    filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all,
    da_ff::Function=filterNothing,
    ff_request::Dict{Symbol,Vector{String}}=Dict{Symbol,Vector{String}}(),
    theme::Attributes=Theme(),
)::Nothing

    x_plot_params = plotParams(x_quantity)
    y_plot_params = plotParams(y_quantity)

    filter_function, translation, rotation, request = selectFilter(
        filter_mode,
        mergeRequests(x_plot_params.request, y_plot_params.request, ff_request),
    )

    # Set arguments for a log x axis
    if xlog
        x_log        = x_plot_params.unit
        x_unit       = Unitful.NoUnits
        unit_label   = getUnitLabel(0, x_plot_params.unit; latex=true)
        if isempty(unit_label)
            xaxis_label  = L"$\log_{10} \, $auto_label"
        else
            xaxis_label  = L"$\log_{10} \, $auto_label [%$(unit_label)]"
        end
        x_exp_factor = 0
    else
        x_log        = nothing
        x_unit       = x_plot_params.unit
        xaxis_label  = x_plot_params.axis_label
        x_exp_factor = x_plot_params.exp_factor
    end

    # Set arguments for a log y axis
    if ylog
        x_log        = y_plot_params.unit
        y_unit       = Unitful.NoUnits
        unit_label   = getUnitLabel(0, y_plot_params.unit; latex=true)
        if isempty(unit_label)
            yaxis_label  = L"$\log_{10} \, $auto_label"
        else
            yaxis_label  = L"$\log_{10} \, $auto_label [%$(unit_label)]"
        end
        y_exp_factor = 0
    else
        y_log        = nothing
        y_unit       = y_plot_params.unit
        yaxis_label  = y_plot_params.axis_label
        y_exp_factor = y_plot_params.exp_factor
    end

    @inbounds for simulation_path in simulation_paths

        # Get the simulation name as a string
        sim_name = basename(simulation_path)

        plotSnapshot(
            [simulation_path],
            request,
            [scatter!];
            pf_kwargs=[(; markersize=2)],
            # `plotSnapshot` configuration
            output_path,
            base_filename="$(sim_name)-$(y_quantity)-vs-$(x_quantity)",
            output_format=".png",
            warnings=false,
            show_progress=true,
            # Data manipulation options
            slice=slice,
            filter_function,
            da_functions=[daScatterGalaxy],
            da_args=[(x_quantity, y_quantity)],
            da_kwargs=[(; x_log, y_log, filter_function=da_ff)],
            post_processing=getNothing,
            pp_args=(),
            pp_kwargs=(;),
            transform_box=true,
            translation,
            rotation,
            smooth=0,
            x_unit,
            y_unit,
            x_exp_factor,
            y_exp_factor,
            x_trim=(-Inf, Inf),
            y_trim=(-Inf, Inf),
            x_edges=false,
            y_edges=false,
            x_func=identity,
            y_func=identity,
            # Axes options
            xaxis_label,
            yaxis_label,
            xaxis_var_name=x_plot_params.var_name,
            yaxis_var_name=y_plot_params.var_name,
            xaxis_scale_func=identity,
            yaxis_scale_func=identity,
            # Plotting and animation options
            save_figures=true,
            backup_results=false,
            theme,
            sim_labels=nothing,
            title="",
            colorbar=false,
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
        slice::IndexType,
        x_quantity::Symbol,
        y_quantity::Symbol,
        z_quantity::Symbol,
        z_unit::Uniful.Units;
        <keyword arguments>
    )::Nothing

Plot two quantities as a density scatter plot (2D histogram), weighted by `z_quantity`.

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`.
  - `slice::IndexType`: Slice of the simulations, i.e. which snapshots will be plotted. It can be an integer (a single snapshot), a vector of integers (several snapshots), an `UnitRange` (e.g. 5:13), an `StepRange` (e.g. 5:2:13) or (:) (all snapshots). Starts at 1 and out of bounds indices are ignored.
  - `x_quantity::Symbol`: Quantity for the x axis. The options are:

      + `:stellar_mass`                -> Stellar mass.
      + `:gas_mass`                    -> Gas mass.
      + `:hydrogen_mass`               -> Hydrogen mass.
      + `:dm_mass`                     -> Dark matter mass.
      + `:bh_mass`                     -> Black hole mass.
      + `:molecular_mass`              -> Molecular hydrogen (``\\mathrm{H_2}``) mass.
      + `:br_molecular_mass`           -> Molecular hydrogen (``\\mathrm{H_2}``) mass, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_mass`                 -> Atomic hydrogen (``\\mathrm{HI}``) mass.
      + `:ionized_mass`                -> Ionized hydrogen (``\\mathrm{HII}``) mass.
      + `:neutral_mass`                -> Neutral hydrogen (``\\mathrm{HI + H_2}``) mass.
      + `:molecular_fraction`          -> Gas mass fraction of molecular hydrogen.
      + `:br_molecular_fraction`       -> Gas mass fraction of molecular hydrogen, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_fraction`             -> Gas mass fraction of atomic hydrogen.
      + `:ionized_fraction`            -> Gas mass fraction of ionized hydrogen.
      + `:neutral_fraction`            -> Gas mass fraction of neutral hydrogen.
      + `:molecular_neutral_fraction`  -> Fraction of molecular hydrogen in the neutral gas.
      + `:mol_eq_quotient`             -> Equilibrium quotient for the molecular fraction equation of the SF model.
      + `:ion_eq_quotient`             -> Equilibrium quotient for the ionized fraction equation of the SF model.
      + `:gas_mass_density`            -> Gas mass density.
      + `:hydrogen_mass_density`       -> Hydrogen mass density.
      + `:gas_number_density`          -> Gas number density.
      + `:molecular_number_density`    -> Molecular hydrogen number density.
      + `:br_molecular_number_density` -> Molecular hydrogen number density, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_number_density`       -> Atomic hydrogen number density.
      + `:ionized_number_density`      -> Ionized hydrogen number density.
      + `:neutral_number_density`      -> Neutral hydrogen number density.
      + `:gas_metallicity`             -> Mass fraction of all elements above He in the gas (solar units).
      + `:stellar_metallicity`         -> Mass fraction of all elements above He in the stars (solar units).
      + `:X_gas_abundance`             -> Gas abundance of element ``\\mathrm{X}``, as ``12 + \\log_{10}(\\mathrm{X \\, / \\, H})``. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
      + `:X_stellar_abundance`         -> Stellar abundance of element ``\\mathrm{X}``, as ``12 + \\log_{10}(\\mathrm{X \\, / \\, H})``. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
      + `:stellar_radial_distance`     -> Distance of every stellar particle to the origin.
      + `:gas_radial_distance`         -> Distance of every gas cell to the origin.
      + `:dm_radial_distance`          -> Distance of every dark matter particle to the origin.
      + `:stellar_xy_distance`         -> Projected distance of every stellar particle to the origin.
      + `:gas_xy_distance`             -> Projected distance of every gas cell to the origin.
      + `:dm_xy_distance`              -> Projected distance of every dark matter particle to the origin.
      + `:gas_sfr`                     -> SFR associated to each gas particle/cell within the code.
      + `:stellar_circularity`         -> Stellar circularity.
      + `:stellar_vcirc`               -> Stellar circular velocity.
      + `:stellar_vradial`             -> Stellar radial speed.
      + `:stellar_vtangential`         -> Stellar tangential speed.
      + `:stellar_vzstar`              -> Stellar speed in the z direction, computed as ``v_z \\, \\mathrm{sign}(z)``.
      + `:stellar_age`                 -> Stellar age.
      + `:sfr`                         -> The star formation rate.
      + `:ssfr`                        -> The specific star formation rate.
      + `:observational_sfr`           -> The star formation rate of the last `AGE_RESOLUTION`.
      + `:observational_ssfr`          -> The specific star formation rate of the last `AGE_RESOLUTION`.
      + `:temperature`                 -> Gas temperature, as ``\\log_{10}(T \\, / \\, \\mathrm{K})``.
      + `:pressure`                    -> Gas pressure.
  - `y_quantity::Symbol`: Quantity for the y axis. The options are:

      + `:stellar_mass`                -> Stellar mass.
      + `:gas_mass`                    -> Gas mass.
      + `:hydrogen_mass`               -> Hydrogen mass.
      + `:dm_mass`                     -> Dark matter mass.
      + `:bh_mass`                     -> Black hole mass.
      + `:molecular_mass`              -> Molecular hydrogen (``\\mathrm{H_2}``) mass.
      + `:br_molecular_mass`           -> Molecular hydrogen (``\\mathrm{H_2}``) mass, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_mass`                 -> Atomic hydrogen (``\\mathrm{HI}``) mass.
      + `:ionized_mass`                -> Ionized hydrogen (``\\mathrm{HII}``) mass.
      + `:neutral_mass`                -> Neutral hydrogen (``\\mathrm{HI + H_2}``) mass.
      + `:molecular_fraction`          -> Gas mass fraction of molecular hydrogen.
      + `:br_molecular_fraction`       -> Gas mass fraction of molecular hydrogen, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_fraction`             -> Gas mass fraction of atomic hydrogen.
      + `:ionized_fraction`            -> Gas mass fraction of ionized hydrogen.
      + `:neutral_fraction`            -> Gas mass fraction of neutral hydrogen.
      + `:molecular_neutral_fraction`  -> Fraction of molecular hydrogen in the neutral gas.
      + `:mol_eq_quotient`             -> Equilibrium quotient for the molecular fraction equation of the SF model.
      + `:ion_eq_quotient`             -> Equilibrium quotient for the ionized fraction equation of the SF model.
      + `:gas_mass_density`            -> Gas mass density.
      + `:hydrogen_mass_density`       -> Hydrogen mass density.
      + `:gas_number_density`          -> Gas number density.
      + `:molecular_number_density`    -> Molecular hydrogen number density.
      + `:br_molecular_number_density` -> Molecular hydrogen number density, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_number_density`       -> Atomic hydrogen number density.
      + `:ionized_number_density`      -> Ionized hydrogen number density.
      + `:neutral_number_density`      -> Neutral hydrogen number density.
      + `:gas_metallicity`             -> Mass fraction of all elements above He in the gas (solar units).
      + `:stellar_metallicity`         -> Mass fraction of all elements above He in the stars (solar units).
      + `:X_gas_abundance`             -> Gas abundance of element ``\\mathrm{X}``, as ``12 + \\log_{10}(\\mathrm{X \\, / \\, H})``. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
      + `:X_stellar_abundance`         -> Stellar abundance of element ``\\mathrm{X}``, as ``12 + \\log_{10}(\\mathrm{X \\, / \\, H})``. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
      + `:stellar_radial_distance`     -> Distance of every stellar particle to the origin.
      + `:gas_radial_distance`         -> Distance of every gas cell to the origin.
      + `:dm_radial_distance`          -> Distance of every dark matter particle to the origin.
      + `:stellar_xy_distance`         -> Projected distance of every stellar particle to the origin.
      + `:gas_xy_distance`             -> Projected distance of every gas cell to the origin.
      + `:dm_xy_distance`              -> Projected distance of every dark matter particle to the origin.
      + `:gas_sfr`                     -> SFR associated to each gas particle/cell within the code.
      + `:stellar_circularity`         -> Stellar circularity.
      + `:stellar_vcirc`               -> Stellar circular velocity.
      + `:stellar_vradial`             -> Stellar radial speed.
      + `:stellar_vtangential`         -> Stellar tangential speed.
      + `:stellar_vzstar`              -> Stellar speed in the z direction, computed as ``v_z \\, \\mathrm{sign}(z)``.
      + `:stellar_age`                 -> Stellar age.
      + `:sfr`                         -> The star formation rate.
      + `:ssfr`                        -> The specific star formation rate.
      + `:observational_sfr`           -> The star formation rate of the last `AGE_RESOLUTION`.
      + `:observational_ssfr`          -> The specific star formation rate of the last `AGE_RESOLUTION`.
      + `:temperature`                 -> Gas temperature, as ``\\log_{10}(T \\, / \\, \\mathrm{K})``.
      + `:pressure`                    -> Gas pressure.
  - `z_quantity::Symbol`: Quantity for the z axis (weights). The options are:

      + `:stellar_mass`                -> Stellar mass.
      + `:gas_mass`                    -> Gas mass.
      + `:hydrogen_mass`               -> Hydrogen mass.
      + `:dm_mass`                     -> Dark matter mass.
      + `:bh_mass`                     -> Black hole mass.
      + `:molecular_mass`              -> Molecular hydrogen (``\\mathrm{H_2}``) mass.
      + `:br_molecular_mass`           -> Molecular hydrogen (``\\mathrm{H_2}``) mass, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_mass`                 -> Atomic hydrogen (``\\mathrm{HI}``) mass.
      + `:ionized_mass`                -> Ionized hydrogen (``\\mathrm{HII}``) mass.
      + `:neutral_mass`                -> Neutral hydrogen (``\\mathrm{HI + H_2}``) mass.
      + `:molecular_fraction`          -> Gas mass fraction of molecular hydrogen.
      + `:br_molecular_fraction`       -> Gas mass fraction of molecular hydrogen, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_fraction`             -> Gas mass fraction of atomic hydrogen.
      + `:ionized_fraction`            -> Gas mass fraction of ionized hydrogen.
      + `:neutral_fraction`            -> Gas mass fraction of neutral hydrogen.
      + `:molecular_neutral_fraction`  -> Fraction of molecular hydrogen in the neutral gas.
      + `:mol_eq_quotient`             -> Equilibrium quotient for the molecular fraction equation of the SF model.
      + `:ion_eq_quotient`             -> Equilibrium quotient for the ionized fraction equation of the SF model.
      + `:gas_mass_density`            -> Gas mass density.
      + `:hydrogen_mass_density`       -> Hydrogen mass density.
      + `:gas_number_density`          -> Gas number density.
      + `:molecular_number_density`    -> Molecular hydrogen number density.
      + `:br_molecular_number_density` -> Molecular hydrogen number density, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_number_density`       -> Atomic hydrogen number density.
      + `:ionized_number_density`      -> Ionized hydrogen number density.
      + `:neutral_number_density`      -> Neutral hydrogen number density.
      + `:gas_metallicity`             -> Mass fraction of all elements above He in the gas (solar units).
      + `:stellar_metallicity`         -> Mass fraction of all elements above He in the stars (solar units).
      + `:X_gas_abundance`             -> Gas abundance of element ``\\mathrm{X}``, as ``12 + \\log_{10}(\\mathrm{X \\, / \\, H})``. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
      + `:X_stellar_abundance`         -> Stellar abundance of element ``\\mathrm{X}``, as ``12 + \\log_{10}(\\mathrm{X \\, / \\, H})``. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
      + `:stellar_radial_distance`     -> Distance of every stellar particle to the origin.
      + `:gas_radial_distance`         -> Distance of every gas cell to the origin.
      + `:dm_radial_distance`          -> Distance of every dark matter particle to the origin.
      + `:stellar_xy_distance`         -> Projected distance of every stellar particle to the origin.
      + `:gas_xy_distance`             -> Projected distance of every gas cell to the origin.
      + `:dm_xy_distance`              -> Projected distance of every dark matter particle to the origin.
      + `:gas_sfr`                     -> SFR associated to each gas particle/cell within the code.
      + `:stellar_circularity`         -> Stellar circularity.
      + `:stellar_vcirc`               -> Stellar circular velocity.
      + `:stellar_vradial`             -> Stellar radial speed.
      + `:stellar_vtangential`         -> Stellar tangential speed.
      + `:stellar_vzstar`              -> Stellar speed in the z direction, computed as ``v_z \\, \\mathrm{sign}(z)``.
      + `:stellar_age`                 -> Stellar age.
      + `:sfr`                         -> The star formation rate.
      + `:ssfr`                        -> The specific star formation rate.
      + `:observational_sfr`           -> The star formation rate of the last `AGE_RESOLUTION`.
      + `:observational_ssfr`          -> The specific star formation rate of the last `AGE_RESOLUTION`.
      + `:temperature`                 -> Gas temperature, as ``\\log_{10}(T \\, / \\, \\mathrm{K})``.
      + `:pressure`                    -> Gas pressure.
  - `z_unit::Unitful.Units`: Target unit for the z axis.
  - `x_range::Union{NTuple{2,<:Number},Nothing}=nothing`: x axis range. If set to `nothing`, the extrema of the values will be used.
  - `y_range::Union{NTuple{2,<:Number},Nothing}=nothing`: y axis range. If set to `nothing`, the extrema of the values will be used.
  - `xlog::Bool=false`: If true, sets everything so the x axis is log10(`x_quantity`).
  - `ylog::Bool=false`: If true, sets everything so the y axis is log10(`y_quantity`).
  - `total::Bool=true`: If the sum (default) or the mean of `z_quantity` will be used as the value of each pixel.
  - `n_bins::Int=100`: Number of bins per side of the square grid.
  - `output_path::String="./"`: Path to the output folder.
  - `filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all`: Which cells/particles will be plotted, the options are:

      + `:all`             -> Consider every cell/particle within the simulation box.
      + `:halo`            -> Consider only the cells/particles that belong to the main halo.
      + `:subhalo`         -> Consider only the cells/particles that belong to the main subhalo.
      + `:sphere`          -> Consider only the cell/particle inside a sphere with radius `DISK_R` (see `./src/constants/globals.jl`).
      + `:stellar_subhalo` -> Consider only the cells/particles that belong to the main subhalo.
      + `:all_subhalo`     -> Plot every cell/particle centered around the main subhalo.
      + A dictionary with three entries:

          + `:filter_function` -> The filter function.
          + `:translation`     -> Translation for the simulation box. The posibilities are:

              + `:global_cm`                  -> Selects the center of mass of the whole system as the new origin.
              + `:{component}`                -> Sets the center of mass of the given component (e.g. :stars, :gas, :halo, etc, after filtering) as the new origin. It can be any of the keys of [`PARTICLE_INDEX`](@ref).
              + `(halo_idx, subhalo_rel_idx)` -> Sets the position of the potencial minimum for the `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo) as the new origin.
              + `(halo_idx, 0)`               -> Sets the center of mass of the `halo_idx::Int` halo as the new origin.
              + `subhalo_abs_idx`             -> Sets the center of mass of the `subhalo_abs_idx::Int` as the new origin.
          + `:rotation`        -> Rotation for the simulation box. The posibilities are:

              + `:zero`                       -> No rotation is appplied.
              + `:global_am`                  -> Sets the angular momentum of the whole system as the new z axis.
              + `:stellar_am`                 -> Sets the stellar angular momentum as the new z axis.
              + `:stellar_pa`                 -> Sets the stellar principal axis as the new coordinate system.
              + `:stellar_subhalo_pa`         -> Sets the principal axis of the stars in the main subhalo as the new coordinate system.
              + `(halo_idx, subhalo_rel_idx)` -> Sets the principal axis of the stars in `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo), as the new coordinate system.
              + `(halo_idx, 0)`               -> Sets the principal axis of the stars in the `halo_idx::Int` halo, as the new coordinate system.
              + `subhalo_abs_idx`             -> Sets the principal axis of the stars in the `subhalo_abs_idx::Int` subhalo as the new coordinate system.
  - `da_ff::Function=filterNothing`: A function with the signature:

    `da_ff(data_dict) -> indices`

    where

      + `data_dict::Dict`: A dictionary with the following shape:

        * `:sim_data`          -> ::Simulation (see [`Simulation`](@ref)).
        * `:snap_data`         -> ::Snapshot (see [`Snapshot`](@ref)).
        * `:gc_data`           -> ::GroupCatalog (see [`GroupCatalog`](@ref)).
        * `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * ...
        * `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * ...
      + `indices::Dict`: A dictionary with the following shape:

        * `cell/particle type` -> idxs::IndexType
        * `cell/particle type` -> idxs::IndexType
        * `cell/particle type` -> idxs::IndexType
        * ...
  - `ff_request::Dict{Symbol,Vector{String}}=Dict{Symbol,Vector{String}}()`: Request dictionary for the `da_ff` filter function.
  - `title::Union{Symbol,<:AbstractString}=""`: Title for the figure. If left empty, no title is printed. It can also be set to one of the following options:

      + `:physical_time` -> Physical time since the Big Bang.
      + `:lookback_time` -> Physical time left to reach the last snapshot.
      + `:scale_factor`  -> Scale factor (only relevant for cosmological simulations).
      + `:redshift`      -> Redshift (only relevant for cosmological simulations).
  - `colorbar::Bool=false`: If a colorbar will be added.
  - `print_range::Bool=false`: Print an info block detailing the color range.
  - `theme::Attributes=Theme()`: Plot theme that will take precedence over [`DEFAULT_THEME`](@ref).
"""
function scatterDensityMap(
    simulation_paths::Vector{String},
    slice::IndexType,
    x_quantity::Symbol,
    y_quantity::Symbol,
    z_quantity::Symbol,
    z_unit::Unitful.Units;
    x_range::Union{NTuple{2,<:Number},Nothing}=nothing,
    y_range::Union{NTuple{2,<:Number},Nothing}=nothing,
    xlog::Bool=false,
    ylog::Bool=false,
    total::Bool=true,
    n_bins::Int=100,
    output_path::String="./",
    filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all,
    da_ff::Function=filterNothing,
    ff_request::Dict{Symbol,Vector{String}}=Dict{Symbol,Vector{String}}(),
    title::Union{Symbol,<:AbstractString}="",
    colorbar::Bool=false,
    print_range::Bool=false,
    theme::Attributes=Theme(),
)::Nothing

    x_plot_params = plotParams(x_quantity)
    y_plot_params = plotParams(y_quantity)
    z_plot_params = plotParams(z_quantity)

    filter_function, translation, rotation, request = selectFilter(
        filter_mode,
        mergeRequests(
            x_plot_params.request,
            y_plot_params.request,
            z_plot_params.request,
            ff_request,
        ),
    )

    n_sims = length(simulation_paths)

    # Set arguments for a log x axis
    if xlog
        x_log        = x_plot_params.unit
        x_unit       = Unitful.NoUnits
        unit_label   = getUnitLabel(0, x_plot_params.unit; latex=true)
        if isempty(unit_label)
            xaxis_label  = L"$\log_{10} \, $auto_label"
        else
            xaxis_label  = L"$\log_{10} \, $auto_label [%$(unit_label)]"
        end
        x_exp_factor = 0
    else
        x_log        = nothing
        x_unit       = x_plot_params.unit
        xaxis_label  = x_plot_params.axis_label
        x_exp_factor = x_plot_params.exp_factor
    end

    # Set arguments for a log y axis
    if ylog
        y_log        = y_plot_params.unit
        y_unit       = Unitful.NoUnits
        unit_label   = getUnitLabel(0, y_plot_params.unit; latex=true)
        if isempty(unit_label)
            yaxis_label  = L"$\log_{10} \, $auto_label"
        else
            yaxis_label  = L"$\log_{10} \, $auto_label [%$(unit_label)]"
        end
        y_exp_factor = 0
    else
        y_log        = nothing
        y_unit       = y_plot_params.unit
        yaxis_label  = y_plot_params.axis_label
        y_exp_factor = y_plot_params.exp_factor
    end

    @inbounds for simulation_path in simulation_paths

        # Get the simulation name as a string
        sim_name = basename(simulation_path)

        if isone(n_sims)
            base_filename="$(y_quantity)-vs-$(x_quantity)"
        else
            base_filename="$(sim_name)-$(y_quantity)-vs-$(x_quantity)"
        end

        plotSnapshot(
            [simulation_path],
            request,
            [heatmap!];
            pf_kwargs=[(;)],
            # `plotSnapshot` configuration
            output_path,
            base_filename,
            output_format=".png",
            warnings=false,
            show_progress=true,
            # Data manipulation options
            slice,
            filter_function,
            da_functions=[daScatterWeightedDensity],
            da_args=[(x_quantity, y_quantity, z_quantity, z_unit)],
            da_kwargs=[
                (;
                    x_range,
                    y_range,
                    x_log,
                    y_log,
                    total,
                    n_bins,
                    print_range,
                    filter_function=da_ff,
                ),
            ],
            post_processing=getNothing,
            pp_args=(),
            pp_kwargs=(;),
            transform_box=true,
            translation,
            rotation,
            smooth=0,
            x_unit,
            y_unit,
            x_exp_factor,
            y_exp_factor,
            x_trim=(-Inf, Inf),
            y_trim=(-Inf, Inf),
            x_edges=false,
            y_edges=false,
            x_func=identity,
            y_func=identity,
            # Axes options
            xaxis_label,
            yaxis_label,
            xaxis_var_name=x_plot_params.var_name,
            yaxis_var_name=y_plot_params.var_name,
            xaxis_scale_func=identity,
            yaxis_scale_func=identity,
            # Plotting and animation options
            save_figures=true,
            backup_results=false,
            theme,
            sim_labels=nothing,
            title,
            colorbar,
            # Animation options
            animation=false,
            animation_filename="animation.mp4",
            framerate=10,
        )

    end

    return nothing

end

"""
    atomicMolecularTransition(
        simulation_paths::Vector{String},
        slice::IndexType,
        ranges::Vector{<:Tuple{<:Real,<:Real}};
        <keyword arguments>
    )::Nothing

Plot the atomic gas to molecular gas transition for a set of metallicity ranges.

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`.
  - `slice::IndexType`: Slice of the simulations, i.e. which snapshots will be plotted. It can be an integer (a single snapshot), a vector of integers (several snapshots), an `UnitRange` (e.g. 5:13), an `StepRange` (e.g. 5:2:13) or (:) (all snapshots). Starts at 1 and out of bounds indices are ignored.
  - `ranges::Vector{<:Tuple{<:Real,<:Real}}`: Metallicity (as in the fractional mass of metals) ranges.
  - `plot_type::Symbol=:heatmap`: Type of plot. The options are:

      + `:heatmap` -> Heatmap. One figure per range will be produced.
      + `:scatter` -> Scatter plot. A single figure with every range will be produced.
  - `halo_idx::Int`: Index of the target halo (FoF group). Starts at 1.
  - `subhalo_rel_idx::Int`: Index of the target subhalo (subfind), relative to the target halo. Starts at 1. If it is set to 0, all subhalos of the target halo are included.
  - `output_path::String="./"`: Path to the output folder.
  - `theme::Attributes=Theme()`: Plot theme that will take precedence over [`DEFAULT_THEME`](@ref).
"""
function atomicMolecularTransition(
    simulation_paths::Vector{String},
    slice::IndexType,
    ranges::Vector{<:Tuple{<:Real,<:Real}};
    plot_type::Symbol=:heatmap,
    halo_idx::Int=1,
    subhalo_rel_idx::Int=1,
    output_path::String="./",
    theme::Attributes=Theme(),
)::Nothing

    # Set some plotting parameters
    x_quantity = :atomic_number_density
    y_quantity = :molecular_neutral_fraction

    x_plot_params = plotParams(x_quantity)
    y_plot_params = plotParams(y_quantity)

    filter_function, translation, rotation, request = selectFilter(
        :subhalo,
        mergeRequests(
            x_plot_params.request,
            y_plot_params.request,
            Dict(:gas => ["GZ  ", "GMET"], :stars => ["GZ2 ", "GME2"]),
        ),
    )

    @inbounds for simulation_path in simulation_paths

        # Get the simulation name as a string
        sim_name = basename(simulation_path)
        filename = "$(sim_name)-$(y_quantity)-vs-$(x_quantity)"

        if plot_type == :heatmap

            for range in ranges

                plotSnapshot(
                    fill(simulation_path, length(ranges)),
                    request,
                    [heatmap!];
                    pf_kwargs=[(;)],
                    # `plotSnapshot` configuration
                    output_path,
                    base_filename="$(filename)-$(range[1])-Z-$(range[2])",
                    output_format=".png",
                    warnings=false,
                    show_progress=true,
                    # Data manipulation options
                    slice,
                    filter_function,
                    da_functions=[daScatterDensity],
                    da_args=[(x_quantity, y_quantity)],
                    da_kwargs=[
                        (;
                            x_log=x_plot_params.unit,
                            y_log=y_plot_params.unit,
                            filter_function=dd -> filterMetallicity(dd, range[1], range[2]),
                        )
                    ],
                    post_processing=getNothing,
                    pp_args=(),
                    pp_kwargs=(;),
                    transform_box=true,
                    translation,
                    rotation,
                    smooth=0,
                    x_unit=Unitful.NoUnits,
                    y_unit=Unitful.NoUnits,
                    x_exp_factor=x_plot_params.exp_factor,
                    y_exp_factor=y_plot_params.exp_factor,
                    x_trim=(-Inf, Inf),
                    y_trim=(-Inf, Inf),
                    x_edges=false,
                    y_edges=false,
                    x_func=identity,
                    y_func=identity,
                    # Axes options
                    xaxis_label=L"$\log_{10} \, $auto_label $[\mathrm{cm}^{-3}]$",
                    yaxis_label=L"$\log_{10} \, $auto_label",
                    xaxis_var_name=x_plot_params.var_name,
                    yaxis_var_name=y_plot_params.var_name,
                    xaxis_scale_func=identity,
                    yaxis_scale_func=identity,
                    # Plotting and animation options
                    save_figures=true,
                    backup_results=false,
                    theme,
                    sim_labels=nothing,
                    title=L"%$(range[1]) \, < \, Z \, < \, %$(range[2])",
                    colorbar=false,
                    # Animation options
                    animation=false,
                    animation_filename="animation.mp4",
                    framerate=10,
                )

            end

        elseif plot_type == :scatter

            plotSnapshot(
                fill(simulation_path, length(ranges)),
                request,
                [scatter!];
                pf_kwargs=[(; markersize=4)],
                # `plotSnapshot` configuration
                output_path,
                base_filename=filename,
                output_format=".png",
                warnings=false,
                show_progress=true,
                # Data manipulation options
                slice,
                filter_function,
                da_functions=[daScatterGalaxy],
                da_args=[(x_quantity, y_quantity)],
                da_kwargs = [
                    (;filter_function=dd -> filterMetallicity(dd, range[1], range[2])) for
                    range in ranges
                ],
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
                xaxis_scale_func=log10,
                yaxis_scale_func=log10,
                # Plotting and animation options
                save_figures=true,
                backup_results=false,
                theme=merge(
                    theme,
                    Theme(Legend=(nbanks=1, halign=:left, valign=:top, padding=(0, 0, 0, 15)),),
                ),
                sim_labels= ["$(range[1]) < Z < $(range[2])" for range in ranges],
                title="",
                colorbar=false,
                # Animation options
                animation=false,
                animation_filename="animation.mp4",
                framerate=10,
            )

        else

            throw(ArgumentError("atomicMolecularTransition: `plot_type` can only be :heatmap or \
            :scatter, but I got :$(plot_type)"))

        end

    end

    return nothing

end

"""
    gasBarPlot(
        simulation_paths::Vector{String},
        slice::IndexType,
        quantity::Symbol,
        edges::Vector{<:Number};
        <keyword arguments>
    )::Nothing

Plot a bar plot of the gas fractions, where the bins are a given gas `quantity`..

Only for gas cells that have entered out routine.

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`.
  - `slice::IndexType`: Slice of the simulations, i.e. which snapshots will be plotted. It can be an integer (a single snapshot), a vector of integers (several snapshots), an `UnitRange` (e.g. 5:13), an `StepRange` (e.g. 5:2:13) or (:) (all snapshots). Starts at 1 and out of bounds indices are ignored.
  - `quantity::Symbol`: Target quantity. The possibilities are:

      + `:gas_mass`                    -> Gas mass.
      + `:hydrogen_mass`               -> Hydrogen mass.
      + `:molecular_mass`              -> Molecular hydrogen (``\\mathrm{H_2}``) mass.
      + `:br_molecular_mass`           -> Molecular hydrogen (``\\mathrm{H_2}``) mass, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_mass`                 -> Atomic hydrogen (``\\mathrm{HI}``) mass.
      + `:ionized_mass`                -> Ionized hydrogen (``\\mathrm{HII}``) mass.
      + `:neutral_mass`                -> Neutral hydrogen (``\\mathrm{HI + H_2}``) mass.
      + `:molecular_fraction`          -> Gas mass fraction of molecular hydrogen.
      + `:br_molecular_fraction`       -> Gas mass fraction of molecular hydrogen, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_fraction`             -> Gas mass fraction of atomic hydrogen.
      + `:ionized_fraction`            -> Gas mass fraction of ionized hydrogen.
      + `:neutral_fraction`            -> Gas mass fraction of neutral hydrogen.
      + `:molecular_neutral_fraction`  -> Fraction of molecular hydrogen in the neutral gas.
      + `:mol_eq_quotient`             -> Equilibrium quotient for the molecular fraction equation of the SF model.
      + `:ion_eq_quotient`             -> Equilibrium quotient for the ionized fraction equation of the SF model.
      + `:gas_mass_density`            -> Gas mass density.
      + `:hydrogen_mass_density`       -> Hydrogen mass density.
      + `:gas_number_density`          -> Gas number density.
      + `:molecular_number_density`    -> Molecular hydrogen number density.
      + `:br_molecular_number_density` -> Molecular hydrogen number density, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_number_density`       -> Atomic hydrogen number density.
      + `:ionized_number_density`      -> Ionized hydrogen number density.
      + `:neutral_number_density`      -> Neutral hydrogen number density.
      + `:gas_metallicity`             -> Mass fraction of all elements above He in the gas (solar units).
      + `:X_gas_abundance`             -> Gas abundance of element ``\\mathrm{X}``, as ``12 + \\log_{10}(\\mathrm{X \\, / \\, H})``. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
      + `:gas_radial_distance`         -> Distance of every gas cell to the origin.
      + `:gas_xy_distance`             -> Projected distance of every gas cell to the origin.
      + `:gas_sfr`                     -> SFR associated to each gas particle/cell within the code.
      + `:temperature`                 -> Gas temperature, as ``\\log_{10}(T \\, / \\, \\mathrm{K})``.
      + `:pressure`                    -> Gas pressure.
  - `edges::Vector{<:Number}`: A sorted list of bin edges for `quantity`.
  - `include_stars::Bool=false`: If the stars will be included as one of the gas phases.
  - `axis_label::Union{AbstractString,Nothing}=nothing`: Label for the axis. It can contain the string `auto_label`, which will be replaced by the default label: `var_name` / 10^`exp_factor` `unit`. If set to `nothing` a label will be assigned automaticaly.
  - `exp_ticks::Bool=false`: If the axis ticks will be the ``\\log_{10}`` of `edges`.
  - `output_path::String="./"`: Path to the output folder.
  - `filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all`: Which cells/particles will be plotted, the options are:

      + `:all`             -> Consider every cell/particle within the simulation box.
      + `:halo`            -> Consider only the cells/particles that belong to the main halo.
      + `:subhalo`         -> Consider only the cells/particles that belong to the main subhalo.
      + `:sphere`          -> Consider only the cell/particle inside a sphere with radius `DISK_R` (see `./src/constants/globals.jl`).
      + `:stellar_subhalo` -> Consider only the cells/particles that belong to the main subhalo.
      + `:all_subhalo`     -> Plot every cell/particle centered around the main subhalo.
      + A dictionary with three entries:

          + `:filter_function` -> The filter function.
          + `:translation`     -> Translation for the simulation box. The posibilities are:

              + `:global_cm`                  -> Selects the center of mass of the whole system as the new origin.
              + `:{component}`                -> Sets the center of mass of the given component (e.g. :stars, :gas, :halo, etc, after filtering) as the new origin. It can be any of the keys of [`PARTICLE_INDEX`](@ref).
              + `(halo_idx, subhalo_rel_idx)` -> Sets the position of the potencial minimum for the `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo) as the new origin.
              + `(halo_idx, 0)`               -> Sets the center of mass of the `halo_idx::Int` halo as the new origin.
              + `subhalo_abs_idx`             -> Sets the center of mass of the `subhalo_abs_idx::Int` as the new origin.
          + `:rotation`        -> Rotation for the simulation box. The posibilities are:

              + `:zero`                       -> No rotation is appplied.
              + `:global_am`                  -> Sets the angular momentum of the whole system as the new z axis.
              + `:stellar_am`                 -> Sets the stellar angular momentum as the new z axis.
              + `:stellar_pa`                 -> Sets the stellar principal axis as the new coordinate system.
              + `:stellar_subhalo_pa`         -> Sets the principal axis of the stars in the main subhalo as the new coordinate system.
              + `(halo_idx, subhalo_rel_idx)` -> Sets the principal axis of the stars in `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo), as the new coordinate system.
              + `(halo_idx, 0)`               -> Sets the principal axis of the stars in the `halo_idx::Int` halo, as the new coordinate system.
              + `subhalo_abs_idx`             -> Sets the principal axis of the stars in the `subhalo_abs_idx::Int` subhalo as the new coordinate system.
  - `theme::Attributes=Theme()`: Plot theme that will take precedence over [`DEFAULT_THEME`](@ref).
"""
function gasBarPlot(
    simulation_paths::Vector{String},
    slice::IndexType,
    quantity::Symbol,
    edges::Vector{<:Number};
    include_stars::Bool=false,
    axis_label::Union{AbstractString,Nothing}=nothing,
    exp_ticks::Bool=false,
    output_path::String="./",
    filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all,
    theme::Attributes=Theme(),
)::Nothing

    plot_params = plotParams(quantity)

    filter_function, translation, rotation, request = selectFilter(
        filter_mode,
        mergeRequests(
            plot_params.request,
            plotParams(:molecular_mass).request,
            plotParams(:atomic_mass).request,
            plotParams(:ionized_mass).request,
        ),
    )

    # Compute the number of bins for the gas quantity
    n_bins = length(edges) - 1

    # Number of bars per bin
    n_bars = include_stars ? 4 : 3

    # Compute the dodge argument for `barplot!`
    dodge = repeat(1:n_bars, outer=n_bins)

    # Set the color list
    colors = Makie.wong_colors()[[3,4,1,2]]

    # Compute the axis ticks
    if exp_ticks
        tick_nums = log10.(ustrip.(plot_params.unit, edges))
    else
        tick_nums = ustrip.(plot_params.unit, edges)
    end

    ticks = [string(round((tick_nums[i] + tick_nums[i + 1]) / 2, sigdigits=2)) for i in 1:n_bins]

    n_sims = length(simulation_paths)

    @inbounds for simulation_path in simulation_paths

        # Get the simulation name as a string
        sim_name = basename(simulation_path)

        if isone(n_sims)
            if include_stars
                base_filename = "fractions-vs-$(quantity)-barplot"
            else
                base_filename = "fractions-vs-$(quantity)-barplot-no_stars"
            end
        else
            if include_stars
                base_filename = "$(sim_name)-fractions-vs-$(quantity)-barplot"
            else
                base_filename = "$(sim_name)-fractions-vs-$(quantity)-barplot-no_stars"
            end
        end

        plotSnapshot(
            [simulation_path],
            request,
            [barplot!];
            pf_kwargs=[(; dodge, color=colors[dodge])],
            # `plotSnapshot` configuration
            output_path,
            base_filename,
            output_format=".png",
            warnings=false,
            show_progress=true,
            # Data manipulation options
            slice,
            filter_function,
            da_functions=[daGasFractions],
            da_args=[(quantity, edges)],
            da_kwargs=[(; include_stars, filter_function=filterGFM)],
            post_processing=ppBarPlotLabels,
            pp_args=(include_stars,),
            pp_kwargs=(; colors),
            transform_box=true,
            translation,
            rotation,
            smooth=0,
            x_unit=Unitful.NoUnits,
            y_unit=Unitful.NoUnits,
            x_exp_factor=0,
            y_exp_factor=0,
            x_trim=(-Inf, Inf),
            y_trim=(-Inf, Inf),
            x_edges=false,
            y_edges=false,
            x_func=identity,
            y_func=identity,
            # Axes options
            xaxis_label=L"\mathrm{Fraction} \,\, [%]",
            yaxis_label=isnothing(axis_label) ? plot_params.axis_label : axis_label,
            xaxis_var_name="",
            yaxis_var_name=plot_params.var_name,
            xaxis_scale_func=identity,
            yaxis_scale_func=identity,
            # Plotting and animation options
            save_figures=true,
            backup_results=false,
            theme=merge(
                theme,
                Theme(
                    size=(850, 850),
                    Legend=(nbanks=1,),
                    Axis=(
                        limits=(nothing, 105, nothing, nothing),
                        xticks=([0, 50, 100], [L"0.0", L"50", L"100"]),
                        yticks=(1:n_bins, ticks),
                    ),
                    BarPlot=(
                        flip_labels_at=10,
                        label_formatter=barPlotLabelFormater,
                        label_size=include_stars ? 25 : 35,
                    ),
                ),
            ),
            sim_labels=nothing,
            title="",
            colorbar=false,
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
  - `x_quantity::Symbol`: Quantity for the x axis. The options are:

      + `:stellar_mass`              -> Stellar mass.
      + `:gas_mass`                  -> Gas mass.
      + `:hydrogen_mass`             -> Hydrogen mass.
      + `:dm_mass`                   -> Dark matter mass.
      + `:bh_mass`                   -> Black hole mass.
      + `:molecular_mass`            -> Molecular hydrogen (``\\mathrm{H_2}``) mass.
      + `:br_molecular_mass`         -> Molecular hydrogen (``\\mathrm{H_2}``) mass, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_mass`               -> Atomic hydrogen (``\\mathrm{HI}``) mass.
      + `:ionized_mass`              -> Ionized hydrogen (``\\mathrm{HII}``) mass.
      + `:neutral_mass`              -> Neutral hydrogen (``\\mathrm{HI + H_2}``) mass.
      + `:stellar_number`            -> Number of stellar particles.
      + `:gas_number`                -> Number of gas cells.
      + `:dm_number`                 -> Number of dark matter particles.
      + `:bh_number`                 -> Number of black hole particles.
      + `:molecular_fraction`        -> Gas mass fraction of molecular hydrogen.
      + `:br_molecular_fraction`     -> Gas mass fraction of molecular hydrogen, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_fraction`           -> Gas mass fraction of atomic hydrogen.
      + `:ionized_fraction`          -> Gas mass fraction of ionized hydrogen.
      + `:neutral_fraction`          -> Gas mass fraction of neutral hydrogen.
      + `:stellar_area_density`      -> Stellar area mass density, for a radius of `DISK_R`.
      + `:gas_area_density`          -> Gas mass surface density, for a radius of `DISK_R`.
      + `:molecular_area_density`    -> Molecular mass surface density, for a radius of `DISK_R`.
      + `:br_molecular_area_density` -> Molecular mass surface density, for a radius of `DISK_R`, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_area_density`       -> Atomic hydrogen area mass density, for a radius of `DISK_R`.
      + `:ionized_area_density`      -> Ionized hydrogen area mass density, for a radius of `DISK_R`.
      + `:neutral_area_density`      -> Neutral mass surface density, for a radius of `DISK_R`.
      + `:sfr_area_density`          -> Star formation rate area density, for the last `AGE_RESOLUTION` and a radius of `DISK_R`.
      + `:gas_metallicity`           -> Mass fraction of all elements above He in the gas (solar units).
      + `:stellar_metallicity`       -> Mass fraction of all elements above He in the stars (solar units).
      + `:X_gas_abundance`           -> Gas abundance of element ``\\mathrm{X}``, as ``12 + \\log_{10}(\\mathrm{X \\, / \\, H})``. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
      + `:X_stellar_abundance`       -> Stellar abundance of element ``\\mathrm{X}``, as ``12 + \\log_{10}(\\mathrm{X \\, / \\, H})``. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
      + `:stellar_specific_am`       -> Norm of the stellar specific angular momentum.
      + `:gas_specific_am`           -> Norm of the gas specific angular momentum.
      + `:dm_specific_am`            -> Norm of the dark matter specific angular momentum.
      + `:sfr`                       -> The star formation rate.
      + `:ssfr`                      -> The specific star formation rate.
      + `:observational_sfr`         -> The star formation rate of the last `AGE_RESOLUTION`.
      + `:observational_ssfr`        -> The specific star formation rate of the last `AGE_RESOLUTION`.
      + `:scale_factor`              -> Scale factor.
      + `:redshift`                  -> Redshift.
      + `:physical_time`             -> Physical time since the Big Bang.
      + `:lookback_time`             -> Physical time left to reach the last snapshot.
  - `y_quantity::Symbol`: Quantity for the y axis. The options are:

      + `:stellar_mass`              -> Stellar mass.
      + `:gas_mass`                  -> Gas mass.
      + `:hydrogen_mass`             -> Hydrogen mass.
      + `:dm_mass`                   -> Dark matter mass.
      + `:bh_mass`                   -> Black hole mass.
      + `:molecular_mass`            -> Molecular hydrogen (``\\mathrm{H_2}``) mass.
      + `:br_molecular_mass`         -> Molecular hydrogen (``\\mathrm{H_2}``) mass, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_mass`               -> Atomic hydrogen (``\\mathrm{HI}``) mass.
      + `:ionized_mass`              -> Ionized hydrogen (``\\mathrm{HII}``) mass.
      + `:neutral_mass`              -> Neutral hydrogen (``\\mathrm{HI + H_2}``) mass.
      + `:stellar_number`            -> Number of stellar particles.
      + `:gas_number`                -> Number of gas cells.
      + `:dm_number`                 -> Number of dark matter particles.
      + `:bh_number`                 -> Number of black hole particles.
      + `:molecular_fraction`        -> Gas mass fraction of molecular hydrogen.
      + `:br_molecular_fraction`     -> Gas mass fraction of molecular hydrogen, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_fraction`           -> Gas mass fraction of atomic hydrogen.
      + `:ionized_fraction`          -> Gas mass fraction of ionized hydrogen.
      + `:neutral_fraction`          -> Gas mass fraction of neutral hydrogen.
      + `:stellar_area_density`      -> Stellar area mass density, for a radius of `DISK_R`.
      + `:gas_area_density`          -> Gas mass surface density, for a radius of `DISK_R`.
      + `:molecular_area_density`    -> Molecular mass surface density, for a radius of `DISK_R`.
      + `:br_molecular_area_density` -> Molecular mass surface density, for a radius of `DISK_R`, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_area_density`       -> Atomic hydrogen area mass density, for a radius of `DISK_R`.
      + `:ionized_area_density`      -> Ionized hydrogen area mass density, for a radius of `DISK_R`.
      + `:neutral_area_density`      -> Neutral mass surface density, for a radius of `DISK_R`.
      + `:sfr_area_density`          -> Star formation rate area density, for the last `AGE_RESOLUTION` and a radius of `DISK_R`.
      + `:gas_metallicity`           -> Mass fraction of all elements above He in the gas (solar units).
      + `:stellar_metallicity`       -> Mass fraction of all elements above He in the stars (solar units).
      + `:X_gas_abundance`           -> Gas abundance of element ``\\mathrm{X}``, as ``12 + \\log_{10}(\\mathrm{X \\, / \\, H})``. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
      + `:X_stellar_abundance`       -> Stellar abundance of element ``\\mathrm{X}``, as ``12 + \\log_{10}(\\mathrm{X \\, / \\, H})``. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
      + `:stellar_specific_am`       -> Norm of the stellar specific angular momentum.
      + `:gas_specific_am`           -> Norm of the gas specific angular momentum.
      + `:dm_specific_am`            -> Norm of the dark matter specific angular momentum.
      + `:sfr`                       -> The star formation rate.
      + `:ssfr`                      -> The specific star formation rate.
      + `:observational_sfr`         -> The star formation rate of the last `AGE_RESOLUTION`.
      + `:observational_ssfr`        -> The specific star formation rate of the last `AGE_RESOLUTION`.
      + `:scale_factor`              -> Scale factor.
      + `:redshift`                  -> Redshift.
      + `:physical_time`             -> Physical time since the Big Bang.
      + `:lookback_time`             -> Physical time left to reach the last snapshot.
  - `y_log::Bool=true`: If the y axis is will have a log10 scale. Only works if `fraction` = false.
  - `cumulative::Bool=false`: If the `y_quantity` will be accumulated or not.
  - `fraction::Bool=false`: If the `y_quantity` will be represented as a fraction of the last value. If `cumulative` = true, this will apply to the accumulated values.
  - `slice::IndexType=(:)`: Slice of the simulations, i.e. which snapshots will be plotted. It can be an integer (a single snapshot), a vector of integers (several snapshots), an `UnitRange` (e.g. 5:13), an `StepRange` (e.g. 5:2:13) or (:) (all snapshots). Starts at 1 and out of bounds indices are ignored.
  - `output_path::String="./"`: Path to the output folder.
  - `filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all`: Which cells/particles will be plotted, the options are:

      + `:all`             -> Consider every cell/particle within the simulation box.
      + `:halo`            -> Consider only the cells/particles that belong to the main halo.
      + `:subhalo`         -> Consider only the cells/particles that belong to the main subhalo.
      + `:sphere`          -> Consider only the cell/particle inside a sphere with radius `DISK_R` (see `./src/constants/globals.jl`).
      + `:stellar_subhalo` -> Consider only the cells/particles that belong to the main subhalo.
      + `:all_subhalo`     -> Plot every cell/particle centered around the main subhalo.
      + A dictionary with three entries:

          + `:filter_function` -> The filter function.
          + `:translation`     -> Translation for the simulation box. The posibilities are:

              + `:global_cm`                  -> Selects the center of mass of the whole system as the new origin.
              + `:{component}`                -> Sets the center of mass of the given component (e.g. :stars, :gas, :halo, etc, after filtering) as the new origin. It can be any of the keys of [`PARTICLE_INDEX`](@ref).
              + `(halo_idx, subhalo_rel_idx)` -> Sets the position of the potencial minimum for the `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo) as the new origin.
              + `(halo_idx, 0)`               -> Sets the center of mass of the `halo_idx::Int` halo as the new origin.
              + `subhalo_abs_idx`             -> Sets the center of mass of the `subhalo_abs_idx::Int` as the new origin.
          + `:rotation`        -> Rotation for the simulation box. The posibilities are:

              + `:zero`                       -> No rotation is appplied.
              + `:global_am`                  -> Sets the angular momentum of the whole system as the new z axis.
              + `:stellar_am`                 -> Sets the stellar angular momentum as the new z axis.
              + `:stellar_pa`                 -> Sets the stellar principal axis as the new coordinate system.
              + `:stellar_subhalo_pa`         -> Sets the principal axis of the stars in the main subhalo as the new coordinate system.
              + `(halo_idx, subhalo_rel_idx)` -> Sets the principal axis of the stars in `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo), as the new coordinate system.
              + `(halo_idx, 0)`               -> Sets the principal axis of the stars in the `halo_idx::Int` halo, as the new coordinate system.
              + `subhalo_abs_idx`             -> Sets the principal axis of the stars in the `subhalo_abs_idx::Int` subhalo as the new coordinate system.
  - `sim_labels::Union{Vector{String},Nothing}=nothing`: Labels for the plot legend, one per simulation. Set it to `nothing` if you don't want a legend.
  - `backup_results::Bool=false`: If the values to be plotted will be backup in a [JLD2](https://github.com/JuliaIO/JLD2.jl) file.
  - `theme::Attributes=Theme()`: Plot theme that will take precedence over [`DEFAULT_THEME`](@ref).
"""
function timeSeries(
    simulation_paths::Vector{String},
    x_quantity::Symbol,
    y_quantity::Symbol;
    y_log::Bool=true,
    cumulative::Bool=false,
    fraction::Bool=false,
    slice::IndexType=(:),
    output_path::String="./",
    filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all,
    sim_labels::Union{Vector{String},Nothing}=basename.(simulation_paths),
    backup_results::Bool=false,
    theme::Attributes=Theme(),
)::Nothing

    x_plot_params = plotParams(x_quantity)
    y_plot_params = plotParams(y_quantity)

    y_var_name = y_plot_params.var_name

    if cumulative
        y_var_name = "Accumulated $(y_var_name)"
    end

    if fraction
        y_var_name = "Fractional $(y_var_name)"
        filename = "$(y_quantity)-vs-$(x_quantity)_fractional"
    else
        filename = "$(y_quantity)-vs-$(x_quantity)"
    end

    if fraction || !y_log
        yaxis_scale_func = identity
    else
        yaxis_scale_func = log10
    end

    plotTimeSeries(
        simulation_paths,
        [lines!];
        pf_kwargs=[(;)],
        # `plotTimeSeries` configuration
        output_path,
        filename,
        output_format=".png",
        warnings=false,
        show_progress=true,
        # Data manipulation options
        slice,
        da_functions=[daEvolution],
        da_args=[(x_quantity, y_quantity)],
        da_kwargs=[(; filter_mode, smooth=0, cumulative, fraction, scaling=identity, warnings=false)],
        post_processing=getNothing,
        pp_args=(),
        pp_kwargs=(;),
        x_unit=x_plot_params.unit,
        y_unit=fraction ? Unitful.NoUnits : y_plot_params.unit,
        x_exp_factor=x_plot_params.exp_factor,
        y_exp_factor=fraction ? 0 : y_plot_params.exp_factor,
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
        yaxis_var_name=y_var_name,
        xaxis_scale_func=identity,
        yaxis_scale_func,
        # Plotting options
        save_figures=!backup_results,
        backup_results,
        theme,
        sim_labels,
        title="",
    )

    return nothing

end

"""
    gasEvolution(
        simulation_paths::Vector{String};
        <keyword arguments>
    )::Nothing

Plot a time series of the gas components. Either their masses or their fractions.

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`.
  - `fractions::Bool=true`: If the fractions (default), or the masses, will be plotted.
  - `slice::IndexType=(:)`: Slice of the simulations, i.e. which snapshots will be plotted. It can be an integer (a single snapshot), a vector of integers (several snapshots), an `UnitRange` (e.g. 5:13), an `StepRange` (e.g. 5:2:13) or (:) (all snapshots). Starts at 1 and out of bounds indices are ignored.
  - `output_path::String="./"`: Path to the output folder.
  - `filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all`: Which cells/particles will be plotted, the options are:

      + `:all`             -> Consider every cell/particle within the simulation box.
      + `:halo`            -> Consider only the cells/particles that belong to the main halo.
      + `:subhalo`         -> Consider only the cells/particles that belong to the main subhalo.
      + `:sphere`          -> Consider only the cell/particle inside a sphere with radius `DISK_R` (see `./src/constants/globals.jl`).
      + `:stellar_subhalo` -> Consider only the cells/particles that belong to the main subhalo.
      + `:all_subhalo`     -> Plot every cell/particle centered around the main subhalo.
      + A dictionary with three entries:

          + `:filter_function` -> The filter function.
          + `:translation`     -> Translation for the simulation box. The posibilities are:

              + `:global_cm`                  -> Selects the center of mass of the whole system as the new origin.
              + `:{component}`                -> Sets the center of mass of the given component (e.g. :stars, :gas, :halo, etc, after filtering) as the new origin. It can be any of the keys of [`PARTICLE_INDEX`](@ref).
              + `(halo_idx, subhalo_rel_idx)` -> Sets the position of the potencial minimum for the `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo) as the new origin.
              + `(halo_idx, 0)`               -> Sets the center of mass of the `halo_idx::Int` halo as the new origin.
              + `subhalo_abs_idx`             -> Sets the center of mass of the `subhalo_abs_idx::Int` as the new origin.
          + `:rotation`        -> Rotation for the simulation box. The posibilities are:

              + `:zero`                       -> No rotation is appplied.
              + `:global_am`                  -> Sets the angular momentum of the whole system as the new z axis.
              + `:stellar_am`                 -> Sets the stellar angular momentum as the new z axis.
              + `:stellar_pa`                 -> Sets the stellar principal axis as the new coordinate system.
              + `:stellar_subhalo_pa`         -> Sets the principal axis of the stars in the main subhalo as the new coordinate system.
              + `(halo_idx, subhalo_rel_idx)` -> Sets the principal axis of the stars in `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo), as the new coordinate system.
              + `(halo_idx, 0)`               -> Sets the principal axis of the stars in the `halo_idx::Int` halo, as the new coordinate system.
              + `subhalo_abs_idx`             -> Sets the principal axis of the stars in the `subhalo_abs_idx::Int` subhalo as the new coordinate system.
  - `extra_filter::Function=filterNothing`: Filter function that will be applied after the one given by `filter_mode`.
  - `filename::Union{String,Nothing}=nothing`: Name for the output file. If left as `nothing`, the filename will be chosen automaticaly.
  - `theme::Attributes=Theme()`: Plot theme that will take precedence over [`DEFAULT_THEME`](@ref).
"""
function gasEvolution(
    simulation_paths::Vector{String};
    fractions::Bool=false,
    slice::IndexType=(:),
    output_path::String="./",
    filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all,
    extra_filter::Function=filterNothing,
    filename::Union{String,Nothing}=nothing,
    theme::Attributes=Theme(),
)::Nothing

    x_plot_params = plotParams(:physical_time)

    if fractions
        quantities = [:ionized_fraction, :atomic_fraction, :molecular_fraction]
        sim_labels = ["Ionized fraction", "Atomic fraction", "Molecular fraction"]
        y_plot_params = plotParams(:generic_fraction)
    else
        quantities = [:stellar_mass, :hydrogen_mass, :ionized_mass, :atomic_mass, :molecular_mass]
        sim_labels = [
            "Stellar mass",
            "Hydrogen mass",
            "Ionized mass",
            "Atomic mass",
            "Molecular mass",
        ]
        y_plot_params = plotParams(:generic_mass)
    end

    for simulation_path in simulation_paths

        if isnothing(filename)
            if fractions
                filename = "gas_fractions-vs-physical_time-$(basename(simulation_path))"
            else
                filename = "gas_masses-vs-physical_time-$(basename(simulation_path))"
            end
        end

        plotTimeSeries(
            fill(simulation_path, length(quantities)),
            [lines!];
            pf_kwargs=[(;)],
            # `plotTimeSeries` configuration
            output_path,
            filename,
            output_format=".png",
            warnings=false,
            show_progress=true,
            # Data manipulation options
            slice,
            da_functions=[daEvolution],
            da_args=[(:physical_time, quantity) for quantity in quantities],
            da_kwargs=[(; filter_mode, extra_filter, smooth=0, scaling=identity, warnings=false)],
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
            yaxis_scale_func=log10,
            # Plotting options
            save_figures=true,
            backup_results=false,
            theme,
            sim_labels,
            title="",
        )

    end

    return nothing

end

"""
    virialAccretionEvolution(
        simulation_paths::Vector{String};
        <keyword arguments>
    )::Nothing

Plot a time series of the accreted mass into the virial radius.

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`.
  - `slice::IndexType=(:)`: Slice of the simulations, i.e. which snapshots will be plotted. It can be vector of integers (several snapshots), an `UnitRange` (e.g. 5:13), an `StepRange` (e.g. 5:2:13) or (:) (all snapshots). Starts at 1 and out of bounds indices are ignored.
  - `halo_idx::Int=1`: Index of the target halo (FoF group). Starts at 1.
  - `tracers::Bool=false`: If tracers will be use to compute the mass accretion.
  - `smooth::Int=0`: The time series will be smooth out using `smooth` bins. Set it to 0 if you want no smoothing.
  - `output_path::String="./"`: Path to the output folder.
  - `sim_labels::Union{Vector{String},Nothing}=nothing`: Labels for the plot legend, one per simulation. Set it to `nothing` if you don't want a legend.
  - `theme::Attributes=Theme()`: Plot theme that will take precedence over [`DEFAULT_THEME`](@ref).
"""
function virialAccretionEvolution(
    simulation_paths::Vector{String};
    slice::IndexType=(:),
    halo_idx::Int=1,
    tracers::Bool=false,
    smooth::Int=0,
    output_path::String="./",
    sim_labels::Union{Vector{String},Nothing}=basename.(simulation_paths),
    theme::Attributes=Theme(),
)::Nothing

    x_plot_params = plotParams(:physical_time)
    y_plot_params = plotParams(:mass_accretion)

    if tracers
        filename="virial-mass-accretion_with_tracers"
    else
        filename="virial-mass-change_evolution"
    end

    plotTimeSeries(
        simulation_paths,
        [lines!];
        pf_kwargs=[(;)],
        # `plotTimeSeries` configuration
        output_path,
        filename,
        output_format=".png",
        warnings=false,
        show_progress=true,
        # Data manipulation options
        slice,
        da_functions=[daVirialAccretion],
        da_args=[()],
        da_kwargs=[(; filter_mode=:halo, halo_idx, tracers, smooth, warnings=false)],
        post_processing=ppHorizontalFlags!,
        pp_args=([0.0],),
        pp_kwargs=(; colors=[:gray65], line_styles=[nothing], warnings=false),
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
        yaxis_var_name=tracers ? y_plot_params.var_name : "Net mass change",
        xaxis_scale_func=identity,
        yaxis_scale_func=identity,
        # Plotting options
        save_figures=true,
        backup_results=false,
        theme,
        sim_labels,
        title="",
    )

    return nothing

end

"""
    discAccretionEvolution(
        simulation_paths::Vector{String};
        <keyword arguments>
    )::Nothing

Plot a time series of the accreted mass into the disc.

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`.
  - `slice::IndexType=(:)`: Slice of the simulations, i.e. which snapshots will be plotted. It can be vector of integers (several snapshots), an `UnitRange` (e.g. 5:13), an `StepRange` (e.g. 5:2:13) or (:) (all snapshots). Starts at 1 and out of bounds indices are ignored.
  - `max_r::Unitful.Length=DISK_R`: Radius of the cylinder.
  - `max_z::Unitful.Length=5.0u"kpc"`: Half height of the cylinder.
  - `smooth::Int=0`: The time series will be smooth out using `smooth` bins. Set it to 0 if you want no smoothing.
  - `output_path::String="./"`: Path to the output folder.
  - `sim_labels::Union{Vector{String},Nothing}=nothing`: Labels for the plot legend, one per simulation. Set it to `nothing` if you don't want a legend.
  - `theme::Attributes=Theme()`: Plot theme that will take precedence over [`DEFAULT_THEME`](@ref).
"""
function discAccretionEvolution(
    simulation_paths::Vector{String};
    slice::IndexType=(:),
    max_r::Unitful.Length=DISK_R,
    max_z::Unitful.Length=5.0u"kpc",
    smooth::Int=0,
    output_path::String="./",
    sim_labels::Union{Vector{String},Nothing}=basename.(simulation_paths),
    theme::Attributes=Theme(),
)::Nothing

    x_plot_params = plotParams(:physical_time)
    y_plot_params = plotParams(:mass_accretion)

    plotTimeSeries(
        simulation_paths,
        [lines!];
        pf_kwargs=[(;)],
        # `plotTimeSeries` configuration
        output_path,
        filename="disc-mass-accretion_with_tracers",
        output_format=".png",
        warnings=false,
        show_progress=true,
        # Data manipulation options
        slice,
        da_functions=[daDiscAccretion],
        da_args=[()],
        da_kwargs=[(; filter_mode=:halo, max_r, max_z, smooth, warnings=false)],
        post_processing=ppHorizontalFlags!,
        pp_args=([0.0],),
        pp_kwargs=(; colors=[:gray65], line_styles=[nothing], warnings=false),
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
        # Plotting options
        save_figures=true,
        backup_results=false,
        theme,
        sim_labels,
        title="",
    )

    return nothing

end

"""
    rotationCurve(
        simulation_paths::Vector{String},
        slice::IndexType;
        <keyword arguments>
    )::Nothing

Plot the galaxy rotation curve of a set of simulations.

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`.
  - `slice::IndexType`: Slice of the simulations, i.e. which snapshots will be plotted. It can be an integer (a single snapshot), a vector of integers (several snapshots), an `UnitRange` (e.g. 5:13), an `StepRange` (e.g. 5:2:13) or (:) (all snapshots). Starts at 1 and out of bounds indices are ignored.
  - `radius::Unitful.Length=DISK_R`: Maximum radial distance for the rotation curve.
  - `output_path::String="./"`: Path to the output folder.
  - `filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all`: Which cells/particles will be plotted, the options are:

      + `:all`             -> Consider every cell/particle within the simulation box.
      + `:halo`            -> Consider only the cells/particles that belong to the main halo.
      + `:subhalo`         -> Consider only the cells/particles that belong to the main subhalo.
      + `:sphere`          -> Consider only the cell/particle inside a sphere with radius `DISK_R` (see `./src/constants/globals.jl`).
      + `:stellar_subhalo` -> Consider only the cells/particles that belong to the main subhalo.
      + `:all_subhalo`     -> Plot every cell/particle centered around the main subhalo.
      + A dictionary with three entries:

          + `:filter_function` -> The filter function.
          + `:translation`     -> Translation for the simulation box. The posibilities are:

              + `:global_cm`                  -> Selects the center of mass of the whole system as the new origin.
              + `:{component}`                -> Sets the center of mass of the given component (e.g. :stars, :gas, :halo, etc, after filtering) as the new origin. It can be any of the keys of [`PARTICLE_INDEX`](@ref).
              + `(halo_idx, subhalo_rel_idx)` -> Sets the position of the potencial minimum for the `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo) as the new origin.
              + `(halo_idx, 0)`               -> Sets the center of mass of the `halo_idx::Int` halo as the new origin.
              + `subhalo_abs_idx`             -> Sets the center of mass of the `subhalo_abs_idx::Int` as the new origin.
          + `:rotation`        -> Rotation for the simulation box. The posibilities are:

              + `:zero`                       -> No rotation is appplied.
              + `:global_am`                  -> Sets the angular momentum of the whole system as the new z axis.
              + `:stellar_am`                 -> Sets the stellar angular momentum as the new z axis.
              + `:stellar_pa`                 -> Sets the stellar principal axis as the new coordinate system.
              + `:stellar_subhalo_pa`         -> Sets the principal axis of the stars in the main subhalo as the new coordinate system.
              + `(halo_idx, subhalo_rel_idx)` -> Sets the principal axis of the stars in `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo), as the new coordinate system.
              + `(halo_idx, 0)`               -> Sets the principal axis of the stars in the `halo_idx::Int` halo, as the new coordinate system.
              + `subhalo_abs_idx`             -> Sets the principal axis of the stars in the `subhalo_abs_idx::Int` subhalo as the new coordinate system.
  - `sim_labels::Union{Vector{String},Nothing}=basename.(simulation_paths)`: Labels for the plot legend, one per simulation. Set it to `nothing` if you don't want a legend.
  - `theme::Attributes=Theme()`: Plot theme that will take precedence over [`DEFAULT_THEME`](@ref).
"""
function rotationCurve(
    simulation_paths::Vector{String},
    slice::IndexType;
    radius::Unitful.Length=DISK_R,
    output_path::String="./",
    filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all,
    sim_labels::Union{Vector{String},Nothing}=basename.(simulation_paths),
    theme::Attributes=Theme(),
)::Nothing

    x_plot_params = plotParams(:stellar_radial_distance)
    y_plot_params = plotParams(:stellar_vcirc)

    filter_function, translation, rotation, request = selectFilter(
        filter_mode,
        mergeRequests(x_plot_params.request, y_plot_params.request, Dict(:stars => ["GAGE"])),
    )

    plotSnapshot(
        simulation_paths,
        request,
        [lines!];
        pf_kwargs=[(;)],
        # `plotSnapshot` configuration
        output_path,
        base_filename="rotation_curve",
        output_format=".png",
        warnings=false,
        show_progress=true,
        # Data manipulation options
        slice,
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
        smooth=round(Int64, 5 * ustrip(u"kpc", radius)),
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
        # Plotting and animation options
        save_figures=true,
        backup_results=false,
        theme,
        sim_labels,
        title="",
        colorbar=false,
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
        slice::IndexType,
        quantity::Symbol;
        <keyword arguments>
    )::Nothing

Plot a density profile.

!!! note

    This method plots one quantity for several simulations in one figure.

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`.
  - `slice::IndexType`: Slice of the simulations, i.e. which snapshots will be plotted. It can be an integer (a single snapshot), a vector of integers (several snapshots), an `UnitRange` (e.g. 5:13), an `StepRange` (e.g. 5:2:13) or (:) (all snapshots). Starts at 1 and out of bounds indices are ignored.
  - `quantities::Vector{Symbol}`: Quantities for the y axis. The options are:

      + `:stellar_mass`       -> Stellar mass.
      + `:gas_mass`           -> Gas mass.
      + `:hydrogen_mass`      -> Hydrogen mass.
      + `:dm_mass`            -> Dark matter mass.
      + `:bh_mass`            -> Black hole mass.
      + `:molecular_mass`     -> Molecular hydrogen (``\\mathrm{H_2}``) mass.
      + `:br_molecular_mass`  -> Molecular hydrogen (``\\mathrm{H_2}``) mass, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_mass`        -> Atomic hydrogen (``\\mathrm{HI}``) mass.
      + `:ionized_mass`       -> Ionized hydrogen (``\\mathrm{HII}``) mass.
      + `:neutral_mass`       -> Neutral hydrogen (``\\mathrm{HI + H_2}``) mass.
      + `:sfr`                -> The star formation rate.
      + `:ssfr`               -> The specific star formation rate.
      + `:observational_sfr`  -> The star formation rate of the last `AGE_RESOLUTION`.
      + `:observational_ssfr` -> The specific star formation rate of the last `AGE_RESOLUTION`.
  - `cumulative::Bool=false`: If the profile will be accumulated or not.
  - `yscale::Function=identity`: Scaling function for the y axis. The options are the scaling functions accepted by [Makie](https://docs.makie.org/stable/): log10, log2, log, sqrt, Makie.logit, Makie.Symlog10, Makie.pseudolog10, and identity.
  - `radius::Unitful.Length=DISK_R`: Radius of the profile.
  - `n_bins::Int=100`: Number of bins.
  - `output_path::String="./"`: Path to the output folder.
  - `filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all`: Which cells/particles will be plotted, the options are:

      + `:all`             -> Consider every cell/particle within the simulation box.
      + `:halo`            -> Consider only the cells/particles that belong to the main halo.
      + `:subhalo`         -> Consider only the cells/particles that belong to the main subhalo.
      + `:sphere`          -> Consider only the cell/particle inside a sphere with radius `DISK_R` (see `./src/constants/globals.jl`).
      + `:stellar_subhalo` -> Consider only the cells/particles that belong to the main subhalo.
      + `:all_subhalo`     -> Plot every cell/particle centered around the main subhalo.
      + A dictionary with three entries:

          + `:filter_function` -> The filter function.
          + `:translation`     -> Translation for the simulation box. The posibilities are:

              + `:global_cm`                  -> Selects the center of mass of the whole system as the new origin.
              + `:{component}`                -> Sets the center of mass of the given component (e.g. :stars, :gas, :halo, etc, after filtering) as the new origin. It can be any of the keys of [`PARTICLE_INDEX`](@ref).
              + `(halo_idx, subhalo_rel_idx)` -> Sets the position of the potencial minimum for the `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo) as the new origin.
              + `(halo_idx, 0)`               -> Sets the center of mass of the `halo_idx::Int` halo as the new origin.
              + `subhalo_abs_idx`             -> Sets the center of mass of the `subhalo_abs_idx::Int` as the new origin.
          + `:rotation`        -> Rotation for the simulation box. The posibilities are:

              + `:zero`                       -> No rotation is appplied.
              + `:global_am`                  -> Sets the angular momentum of the whole system as the new z axis.
              + `:stellar_am`                 -> Sets the stellar angular momentum as the new z axis.
              + `:stellar_pa`                 -> Sets the stellar principal axis as the new coordinate system.
              + `:stellar_subhalo_pa`         -> Sets the principal axis of the stars in the main subhalo as the new coordinate system.
              + `(halo_idx, subhalo_rel_idx)` -> Sets the principal axis of the stars in `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo), as the new coordinate system.
              + `(halo_idx, 0)`               -> Sets the principal axis of the stars in the `halo_idx::Int` halo, as the new coordinate system.
              + `subhalo_abs_idx`             -> Sets the principal axis of the stars in the `subhalo_abs_idx::Int` subhalo as the new coordinate system.
  - `sim_labels::Union{Vector{String},Nothing}=basename.(simulation_paths)`: Labels for the plot legend, one per simulation. Set it to `nothing` if you don't want a legend.
  - `theme::Attributes=Theme()`: Plot theme that will take precedence over [`DEFAULT_THEME`](@ref).
"""
function densityProfile(
    simulation_paths::Vector{String},
    slice::IndexType,
    quantity::Symbol;
    cumulative::Bool=false,
    yscale::Function=identity,
    radius::Unitful.Length=DISK_R,
    n_bins::Int=100,
    output_path::String="./",
    filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all,
    sim_labels::Union{Vector{String},Nothing}=basename.(simulation_paths),
    theme::Attributes=Theme(),
)::Nothing

    plot_params = plotParams(quantity)
    filter_function, translation, rotation, request = selectFilter(filter_mode, plot_params.request)

    if quantity == :stellar_mass

        yaxis_var_name = L"\Sigma_\star"

    elseif quantity == :dm_mass

        yaxis_var_name = L"\Sigma_\mathrm{DM}"

    elseif quantity == :bh_mass

        yaxis_var_name = L"\Sigma_\mathrm{BH}"

    elseif quantity == :gas_mass

        yaxis_var_name = L"\Sigma_\mathrm{gas}"

    elseif quantity == :hydrogen_mass

        yaxis_var_name = L"\Sigma_\mathrm{H}"

    elseif quantity == :molecular_mass

        yaxis_var_name = L"\Sigma_\mathrm{H_2}"

    elseif quantity == :br_molecular_mass

        yaxis_var_name = L"\Sigma_\mathrm{H_2}^\mathrm{BR}"

    elseif quantity == :atomic_mass

        yaxis_var_name = L"\Sigma_\mathrm{HI}"

    elseif quantity == :ionized_mass

        yaxis_var_name = L"\Sigma_\mathrm{HII}"

    elseif quantity == :neutral_mass

        yaxis_var_name = L"\Sigma_\mathrm{HI + H_2}"

    elseif quantity ∈ [:sfr, :observational_sfr]

        yaxis_var_name = L"\Sigma_\mathrm{SFR}"

    elseif quantity ∈ [:ssfr, :observational_ssfr]

        yaxis_var_name = L"\Sigma_\mathrm{sSFR}"

    else

        throw(ArgumentError("densityProfile: I don't know the quantity :$(quantity)"))

    end

    grid = CircularGrid(radius, n_bins)

    # Draw the figures with CairoMakie
    plotSnapshot(
        simulation_paths,
        request,
        [lines!];
        pf_kwargs=[(;)],
        # `plotSnapshot` configuration
        output_path,
        base_filename="$(quantity)-density_profile",
        output_format=".png",
        warnings=false,
        show_progress=true,
        # Data manipulation options
        slice,
        filter_function,
        da_functions=[daProfile],
        da_args=[(quantity, grid)],
        da_kwargs=[(; flat=true, total=true, cumulative, density=true)],
        post_processing=getNothing,
        pp_args=(),
        pp_kwargs=(;),
        transform_box=true,
        translation,
        rotation,
        smooth=0,
        x_unit=u"kpc",
        y_unit=plot_params.unit * u"kpc^-2",
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
        yaxis_var_name,
        xaxis_scale_func=identity,
        yaxis_scale_func=yscale,
        # Plotting and animation options
        save_figures=true,
        backup_results=false,
        theme,
        sim_labels,
        title="",
        colorbar=false,
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
        slice::IndexType,
        quantities::Vector{Symbol};
        <keyword arguments>
    )::Nothing

Plot a density profile.

!!! note

    This method plots several quantities for one simulations in one figure.

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`.
  - `slice::IndexType`: Slice of the simulations, i.e. which snapshots will be plotted. It can be an integer (a single snapshot), a vector of integers (several snapshots), an `UnitRange` (e.g. 5:13), an `StepRange` (e.g. 5:2:13) or (:) (all snapshots). Starts at 1 and out of bounds indices are ignored.
  - `quantities::Vector{Symbol}`: Quantities for the y axis. The options are:

      + `:stellar_mass`       -> Stellar mass.
      + `:gas_mass`           -> Gas mass.
      + `:hydrogen_mass`      -> Hydrogen mass.
      + `:dm_mass`            -> Dark matter mass.
      + `:bh_mass`            -> Black hole mass.
      + `:molecular_mass`     -> Molecular hydrogen (``\\mathrm{H_2}``) mass.
      + `:br_molecular_mass`  -> Molecular hydrogen (``\\mathrm{H_2}``) mass, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_mass`        -> Atomic hydrogen (``\\mathrm{HI}``) mass.
      + `:ionized_mass`       -> Ionized hydrogen (``\\mathrm{HII}``) mass.
      + `:neutral_mass`       -> Neutral hydrogen (``\\mathrm{HI + H_2}``) mass.
      + `:sfr`                -> The star formation rate.
      + `:ssfr`               -> The specific star formation rate.
      + `:observational_sfr`  -> The star formation rate of the last `AGE_RESOLUTION`.
      + `:observational_ssfr` -> The specific star formation rate of the last `AGE_RESOLUTION`.
  - `cumulative::Bool=false`: If the profile will be accumulated or not.
  - `yscale::Function=identity`: Scaling function for the y axis. The options are the scaling functions accepted by [Makie](https://docs.makie.org/stable/): log10, log2, log, sqrt, Makie.logit, Makie.Symlog10, Makie.pseudolog10, and identity.
  - `radius::Unitful.Length=DISK_R`: Radius of the profile.
  - `n_bins::Int=100`: Number of bins.
  - `output_path::String="./"`: Path to the output folder.
  - `filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all`: Which cells/particles will be plotted, the options are:

      + `:all`             -> Consider every cell/particle within the simulation box.
      + `:halo`            -> Consider only the cells/particles that belong to the main halo.
      + `:subhalo`         -> Consider only the cells/particles that belong to the main subhalo.
      + `:sphere`          -> Consider only the cell/particle inside a sphere with radius `DISK_R` (see `./src/constants/globals.jl`).
      + `:stellar_subhalo` -> Consider only the cells/particles that belong to the main subhalo.
      + `:all_subhalo`     -> Plot every cell/particle centered around the main subhalo.
      + A dictionary with three entries:

          + `:filter_function` -> The filter function.
          + `:translation`     -> Translation for the simulation box. The posibilities are:

              + `:global_cm`                  -> Selects the center of mass of the whole system as the new origin.
              + `:{component}`                -> Sets the center of mass of the given component (e.g. :stars, :gas, :halo, etc, after filtering) as the new origin. It can be any of the keys of [`PARTICLE_INDEX`](@ref).
              + `(halo_idx, subhalo_rel_idx)` -> Sets the position of the potencial minimum for the `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo) as the new origin.
              + `(halo_idx, 0)`               -> Sets the center of mass of the `halo_idx::Int` halo as the new origin.
              + `subhalo_abs_idx`             -> Sets the center of mass of the `subhalo_abs_idx::Int` as the new origin.
          + `:rotation`        -> Rotation for the simulation box. The posibilities are:

              + `:zero`                       -> No rotation is appplied.
              + `:global_am`                  -> Sets the angular momentum of the whole system as the new z axis.
              + `:stellar_am`                 -> Sets the stellar angular momentum as the new z axis.
              + `:stellar_pa`                 -> Sets the stellar principal axis as the new coordinate system.
              + `:stellar_subhalo_pa`         -> Sets the principal axis of the stars in the main subhalo as the new coordinate system.
              + `(halo_idx, subhalo_rel_idx)` -> Sets the principal axis of the stars in `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo), as the new coordinate system.
              + `(halo_idx, 0)`               -> Sets the principal axis of the stars in the `halo_idx::Int` halo, as the new coordinate system.
              + `subhalo_abs_idx`             -> Sets the principal axis of the stars in the `subhalo_abs_idx::Int` subhalo as the new coordinate system.
  - `sim_labels::Union{Vector{String},Nothing}=string.(quantities)`: Labels for the plot legend, one per quantity. Set it to `nothing` if you don't want a legend.
  - `theme::Attributes=Theme()`: Plot theme that will take precedence over [`DEFAULT_THEME`](@ref).
"""
function densityProfile(
    simulation_paths::Vector{String},
    slice::IndexType,
    quantities::Vector{Symbol};
    cumulative::Bool=false,
    yscale::Function=identity,
    radius::Unitful.Length=DISK_R,
    n_bins::Int=100,
    output_path::String="./",
    filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all,
    sim_labels::Union{Vector{String},Nothing}=string.(quantities),
    theme::Attributes=Theme(),
)::Nothing

    plot_params = plotParams(:generic_area_density)
    filter_function, translation, rotation, request = selectFilter(filter_mode, plot_params.request)

    grid = CircularGrid(radius, n_bins)

    @inbounds for simulation_path in simulation_paths

        # Get the simulation name as a string
        sim_name = basename(simulation_path)

        # Draw the figures with CairoMakie
        plotSnapshot(
            fill(simulation_path, length(quantities)),
            request,
            [lines!];
            pf_kwargs=[(;)],
            # `plotSnapshot` configuration
            output_path,
            base_filename="$(sim_name)-density_profiles",
            output_format=".png",
            warnings=false,
            show_progress=true,
            # Data manipulation options
            slice,
            filter_function,
            da_functions=[daProfile],
            da_args=[(quantity, grid) for quantity in quantities],
            da_kwargs=[(; flat=true, total=true, cumulative, density=true)],
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
            yaxis_scale_func=yscale,
            # Plotting and animation options
            save_figures=true,
            backup_results=false,
            theme,
            sim_labels,
            title="",
            colorbar=false,
            # Animation options
            animation=false,
            animation_filename="animation.mp4",
            framerate=10,
        )

    end

    return nothing

end

"""
    massProfile(
        simulation_paths::Vector{String},
        slice::IndexType,
        quantities::Vector{Symbol};
        <keyword arguments>
    )::Nothing

Plot a mass profile.

!!! note

    This method plots several quantities for one simulations in one figure.

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`.
  - `slice::IndexType`: Slice of the simulations, i.e. which snapshots will be plotted. It can be an integer (a single snapshot), a vector of integers (several snapshots), an `UnitRange` (e.g. 5:13), an `StepRange` (e.g. 5:2:13) or (:) (all snapshots). Starts at 1 and out of bounds indices are ignored.
  - `quantities::Vector{Symbol}`: Quantities for the y axis. The options are:

      + `:stellar_mass`      -> Stellar mass.
      + `:gas_mass`          -> Gas mass.
      + `:hydrogen_mass`     -> Hydrogen mass.
      + `:dm_mass`           -> Dark matter mass.
      + `:bh_mass`           -> Black hole mass.
      + `:molecular_mass`    -> Molecular hydrogen (``\\mathrm{H_2}``) mass.
      + `:br_molecular_mass` -> Molecular hydrogen (``\\mathrm{H_2}``) mass, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_mass`       -> Atomic hydrogen (``\\mathrm{HI}``) mass.
      + `:ionized_mass`      -> Ionized hydrogen (``\\mathrm{HII}``) mass.
      + `:neutral_mass`      -> Neutral hydrogen (``\\mathrm{HI + H_2}``) mass.
  - `cumulative::Bool=false`: If the profile will be accumulated or not.
  - `yscale::Function=identity`: Scaling function for the y axis. The options are the scaling functions accepted by [Makie](https://docs.makie.org/stable/): log10, log2, log, sqrt, Makie.logit, Makie.Symlog10, Makie.pseudolog10, and identity.
  - `radius::Unitful.Length=DISK_R`: Radius of the profile.
  - `n_bins::Int=100`: Number of bins.
  - `output_path::String="./"`: Path to the output folder.
  - `filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all`: Which cells/particles will be plotted, the options are:

      + `:all`             -> Consider every cell/particle within the simulation box.
      + `:halo`            -> Consider only the cells/particles that belong to the main halo.
      + `:subhalo`         -> Consider only the cells/particles that belong to the main subhalo.
      + `:sphere`          -> Consider only the cell/particle inside a sphere with radius `DISK_R` (see `./src/constants/globals.jl`).
      + `:stellar_subhalo` -> Consider only the cells/particles that belong to the main subhalo.
      + `:all_subhalo`     -> Plot every cell/particle centered around the main subhalo.
      + A dictionary with three entries:

          + `:filter_function` -> The filter function.
          + `:translation`     -> Translation for the simulation box. The posibilities are:

              + `:global_cm`                  -> Selects the center of mass of the whole system as the new origin.
              + `:{component}`                -> Sets the center of mass of the given component (e.g. :stars, :gas, :halo, etc, after filtering) as the new origin. It can be any of the keys of [`PARTICLE_INDEX`](@ref).
              + `(halo_idx, subhalo_rel_idx)` -> Sets the position of the potencial minimum for the `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo) as the new origin.
              + `(halo_idx, 0)`               -> Sets the center of mass of the `halo_idx::Int` halo as the new origin.
              + `subhalo_abs_idx`             -> Sets the center of mass of the `subhalo_abs_idx::Int` as the new origin.
          + `:rotation`        -> Rotation for the simulation box. The posibilities are:

              + `:zero`                       -> No rotation is appplied.
              + `:global_am`                  -> Sets the angular momentum of the whole system as the new z axis.
              + `:stellar_am`                 -> Sets the stellar angular momentum as the new z axis.
              + `:stellar_pa`                 -> Sets the stellar principal axis as the new coordinate system.
              + `:stellar_subhalo_pa`         -> Sets the principal axis of the stars in the main subhalo as the new coordinate system.
              + `(halo_idx, subhalo_rel_idx)` -> Sets the principal axis of the stars in `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo), as the new coordinate system.
              + `(halo_idx, 0)`               -> Sets the principal axis of the stars in the `halo_idx::Int` halo, as the new coordinate system.
              + `subhalo_abs_idx`             -> Sets the principal axis of the stars in the `subhalo_abs_idx::Int` subhalo as the new coordinate system.
  - `sim_labels::Union{Vector{String},Nothing}=string.(quantities)`: Labels for the plot legend, one per quantity. Set it to `nothing` if you don't want a legend.
  - `theme::Attributes=Theme()`: Plot theme that will take precedence over [`DEFAULT_THEME`](@ref).
"""
function massProfile(
    simulation_paths::Vector{String},
    slice::IndexType,
    quantities::Vector{Symbol};
    cumulative::Bool=false,
    yscale::Function=identity,
    radius::Unitful.Length=DISK_R,
    n_bins::Int=100,
    output_path::String="./",
    filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all,
    sim_labels::Union{Vector{String},Nothing}=string.(quantities),
    theme::Attributes=Theme(),
)::Nothing

    plot_params = plotParams(:generic_mass)
    filter_function, translation, rotation, request = selectFilter(filter_mode, plot_params.request)

    grid = CircularGrid(radius, n_bins)

    n_sims = length(simulation_paths)

    @inbounds for simulation_path in simulation_paths

        # Get the simulation name as a string
        sim_name = basename(simulation_path)

        if isone(n_sims)
            if cumulative
                base_filename = "mass_profiles-cumulative"
            else
                base_filename = "mass_profiles"
            end
        else
            if cumulative
                base_filename = "$(sim_name)-mass_profiles-cumulative"
            else
                base_filename = "$(sim_name)-mass_profiles"
            end
        end

        # Draw the figures with CairoMakie
        plotSnapshot(
            fill(simulation_path, length(quantities)),
            request,
            [lines!];
            pf_kwargs=[(;)],
            # `plotSnapshot` configuration
            output_path,
            base_filename,
            output_format=".png",
            warnings=false,
            show_progress=true,
            # Data manipulation options
            slice,
            filter_function,
            da_functions=[daProfile],
            da_args=[(quantity, grid) for quantity in quantities],
            da_kwargs=[(; flat=true, total=true, cumulative, density=false)],
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
            yaxis_scale_func=yscale,
            # Plotting and animation options
            save_figures=true,
            backup_results=false,
            theme,
            sim_labels,
            title="",
            colorbar=false,
            # Animation options
            animation=false,
            animation_filename="animation.mp4",
            framerate=10,
        )
        end

    return nothing

end

"""
    velocityProfile(
        simulation_paths::Vector{String},
        slice::IndexType,
        quantity::Symbol;
        <keyword arguments>
    )::Nothing

Plot a velocity profile.

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`.
  - `slice::IndexType`: Slice of the simulations, i.e. which snapshots will be plotted. It can be an integer (a single snapshot), a vector of integers (several snapshots), an `UnitRange` (e.g. 5:13), an `StepRange` (e.g. 5:2:13) or (:) (all snapshots). Starts at 1 and out of bounds indices are ignored.
  - `component::Symbol`: Which component will be calculated. The options are:

      + `:stellar_vradial`     -> Stellar radial speed (``v_r``).
      + `:stellar_vtangential` -> Stellar tangential speed (``v_\\theta``).
      + `:stellar_vzstar`      -> Stellar speed in the z direction, computed as ``v_z \\, \\mathrm{sign}(z)``.
  - `yscale::Function=identity`: Scaling function for the y axis. The options are the scaling functions accepted by [Makie](https://docs.makie.org/stable/): log10, log2, log, sqrt, Makie.logit, Makie.Symlog10, Makie.pseudolog10, and identity.
  - `output_path::String="./"`: Path to the output folder.
  - `filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all`: Which cells/particles will be plotted, the options are:

      + `:all`             -> Consider every cell/particle within the simulation box.
      + `:halo`            -> Consider only the cells/particles that belong to the main halo.
      + `:subhalo`         -> Consider only the cells/particles that belong to the main subhalo.
      + `:sphere`          -> Consider only the cell/particle inside a sphere with radius `DISK_R` (see `./src/constants/globals.jl`).
      + `:stellar_subhalo` -> Consider only the cells/particles that belong to the main subhalo.
      + `:all_subhalo`     -> Plot every cell/particle centered around the main subhalo.
      + A dictionary with three entries:

          + `:filter_function` -> The filter function.
          + `:translation`     -> Translation for the simulation box. The posibilities are:

              + `:global_cm`                  -> Selects the center of mass of the whole system as the new origin.
              + `:{component}`                -> Sets the center of mass of the given component (e.g. :stars, :gas, :halo, etc, after filtering) as the new origin. It can be any of the keys of [`PARTICLE_INDEX`](@ref).
              + `(halo_idx, subhalo_rel_idx)` -> Sets the position of the potencial minimum for the `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo) as the new origin.
              + `(halo_idx, 0)`               -> Sets the center of mass of the `halo_idx::Int` halo as the new origin.
              + `subhalo_abs_idx`             -> Sets the center of mass of the `subhalo_abs_idx::Int` as the new origin.
          + `:rotation`        -> Rotation for the simulation box. The posibilities are:

              + `:zero`                       -> No rotation is appplied.
              + `:global_am`                  -> Sets the angular momentum of the whole system as the new z axis.
              + `:stellar_am`                 -> Sets the stellar angular momentum as the new z axis.
              + `:stellar_pa`                 -> Sets the stellar principal axis as the new coordinate system.
              + `:stellar_subhalo_pa`         -> Sets the principal axis of the stars in the main subhalo as the new coordinate system.
              + `(halo_idx, subhalo_rel_idx)` -> Sets the principal axis of the stars in `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo), as the new coordinate system.
              + `(halo_idx, 0)`               -> Sets the principal axis of the stars in the `halo_idx::Int` halo, as the new coordinate system.
              + `subhalo_abs_idx`             -> Sets the principal axis of the stars in the `subhalo_abs_idx::Int` subhalo as the new coordinate system.
  - `sim_labels::Union{Vector{String},Nothing}=basename.(simulation_paths)`: Labels for the plot legend, one per simulation. Set it to `nothing` if you don't want a legend.
  - `theme::Attributes=Theme()`: Plot theme that will take precedence over [`DEFAULT_THEME`](@ref).
"""
function velocityProfile(
    simulation_paths::Vector{String},
    slice::IndexType,
    component::Symbol;
    yscale::Function=identity,
    output_path::String="./",
    filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all,
    sim_labels::Union{Vector{String},Nothing}=basename.(simulation_paths),
    theme::Attributes=Theme(),
)::Nothing

    plot_params = plotParams(component)
    filter_function, translation, rotation, request = selectFilter(filter_mode, plot_params.request)

    grid = CircularGrid(DISK_R, 25)

    # Draw the figures with CairoMakie
    plotSnapshot(
        simulation_paths,
        request,
        [scatterlines!];
        pf_kwargs=[(;)],
        # `plotSnapshot` configuration
        output_path,
        base_filename="$(component)-profile",
        output_format=".png",
        warnings=false,
        show_progress=true,
        # Data manipulation options
        slice,
        filter_function,
        da_functions=[daProfile],
        da_args=[(component, grid)],
        da_kwargs=[(; flat=true, total=false, cumulative=false, density=false)],
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
        yaxis_scale_func=yscale,
        # Plotting and animation options
        save_figures=true,
        backup_results=false,
        theme,
        sim_labels,
        title="",
        colorbar=false,
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
        slice::IndexType,
        quantity::Symbol;
        <keyword arguments>
    )::Nothing

Plot the evolution of a given stellar `quantity` using the stellar ages at a given instant in time.

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`.
  - `slice::IndexType`: Snapshot at which the stellar ages will be read. If set to several snapshots, one plot per snapshot will be done. It can be an integer (a single snapshot), a vector of integers (several snapshots), an `UnitRange` (e.g. 5:13), an `StepRange` (e.g. 5:2:13) or (:) (all snapshots). Starts at 1 and out of bounds indices are ignored.
  - `quantity::Symbol`: Quantity for the y axis. The options are:

      + `:sfr`                 -> The star formation rate.
      + `:ssfr`                -> The specific star formation rate.
      + `:stellar_mass`        -> Stellar mass.
      + `:stellar_metallicity` -> Mass fraction of all elements above He in the stars (solar units).
  - `y_log::Bool=true`: If the y axis is will have a log10 scale.
  - `n_bins::Int=20`: Number of bins (time intervals).
  - `output_path::String="./"`: Path to the output folder.
  - `filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all`: Which cells/particles will be plotted, the options are:

      + `:all`             -> Consider every cell/particle within the simulation box.
      + `:halo`            -> Consider only the cells/particles that belong to the main halo.
      + `:subhalo`         -> Consider only the cells/particles that belong to the main subhalo.
      + `:sphere`          -> Consider only the cell/particle inside a sphere with radius `DISK_R` (see `./src/constants/globals.jl`).
      + `:stellar_subhalo` -> Consider only the cells/particles that belong to the main subhalo.
      + `:all_subhalo`     -> Plot every cell/particle centered around the main subhalo.
      + A dictionary with three entries:

          + `:filter_function` -> The filter function.
          + `:translation`     -> Translation for the simulation box. The posibilities are:

              + `:global_cm`                  -> Selects the center of mass of the whole system as the new origin.
              + `:{component}`                -> Sets the center of mass of the given component (e.g. :stars, :gas, :halo, etc, after filtering) as the new origin. It can be any of the keys of [`PARTICLE_INDEX`](@ref).
              + `(halo_idx, subhalo_rel_idx)` -> Sets the position of the potencial minimum for the `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo) as the new origin.
              + `(halo_idx, 0)`               -> Sets the center of mass of the `halo_idx::Int` halo as the new origin.
              + `subhalo_abs_idx`             -> Sets the center of mass of the `subhalo_abs_idx::Int` as the new origin.
          + `:rotation`        -> Rotation for the simulation box. The posibilities are:

              + `:zero`                       -> No rotation is appplied.
              + `:global_am`                  -> Sets the angular momentum of the whole system as the new z axis.
              + `:stellar_am`                 -> Sets the stellar angular momentum as the new z axis.
              + `:stellar_pa`                 -> Sets the stellar principal axis as the new coordinate system.
              + `:stellar_subhalo_pa`         -> Sets the principal axis of the stars in the main subhalo as the new coordinate system.
              + `(halo_idx, subhalo_rel_idx)` -> Sets the principal axis of the stars in `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo), as the new coordinate system.
              + `(halo_idx, 0)`               -> Sets the principal axis of the stars in the `halo_idx::Int` halo, as the new coordinate system.
              + `subhalo_abs_idx`             -> Sets the principal axis of the stars in the `subhalo_abs_idx::Int` subhalo as the new coordinate system.
  - `sim_labels::Union{Vector{String},Nothing}=basename.(simulation_paths)`: Labels for the plot legend, one per simulation. Set it to `nothing` if you don't want a legend.
  - `theme::Attributes=Theme()`: Plot theme that will take precedence over [`DEFAULT_THEME`](@ref).
"""
function stellarHistory(
    simulation_paths::Vector{String},
    slice::IndexType,
    quantity::Symbol;
    y_log::Bool=true,
    n_bins::Int=20,
    output_path::String="./",
    filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all,
    sim_labels::Union{Vector{String},Nothing}=basename.(simulation_paths),
    theme::Attributes=Theme(),
)::Nothing

    x_plot_params = plotParams(:physical_time)
    y_plot_params = plotParams(quantity)

    filter_function, translation, rotation, request = selectFilter(
        filter_mode,
        y_plot_params.request,
    )

    # Draw the figures with CairoMakie
    plotSnapshot(
        simulation_paths,
        request,
        [lines!];
        pf_kwargs=[(;)],
        # `plotSnapshot` configuration
        output_path,
        base_filename="$(quantity)_history",
        output_format=".png",
        warnings=false,
        show_progress=true,
        # Data manipulation options
        slice,
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
        yaxis_scale_func=y_log ? log10 : identity,
        # Plotting and animation options
        save_figures=true,
        backup_results=false,
        theme,
        sim_labels,
        title="",
        colorbar=false,
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
        slice::IndexType;
        <keyword arguments>
    )::Nothing

Plot a histogram of the stellar circularity.

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`.
  - `slice::IndexType`: Slice of the simulations, i.e. which snapshots will be plotted. It can be an integer (a single snapshot), a vector of integers (several snapshots), an `UnitRange` (e.g. 5:13), an `StepRange` (e.g. 5:2:13) or (:) (all snapshots). Starts at 1 and out of bounds indices are ignored.
  - `range::NTuple{2,<:Number}=(-2.0, 2.0)`: Circularity range.
  - `n_bins::Int=60`: Number of bins.
  - `output_path::String="./"`: Path to the output folder.
  - `filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all`: Which cells/particles will be plotted, the options are:

      + `:all`             -> Consider every cell/particle within the simulation box.
      + `:halo`            -> Consider only the cells/particles that belong to the main halo.
      + `:subhalo`         -> Consider only the cells/particles that belong to the main subhalo.
      + `:sphere`          -> Consider only the cell/particle inside a sphere with radius `DISK_R` (see `./src/constants/globals.jl`).
      + `:stellar_subhalo` -> Consider only the cells/particles that belong to the main subhalo.
      + `:all_subhalo`     -> Plot every cell/particle centered around the main subhalo.
      + A dictionary with three entries:

          + `:filter_function` -> The filter function.
          + `:translation`     -> Translation for the simulation box. The posibilities are:

              + `:global_cm`                  -> Selects the center of mass of the whole system as the new origin.
              + `:{component}`                -> Sets the center of mass of the given component (e.g. :stars, :gas, :halo, etc, after filtering) as the new origin. It can be any of the keys of [`PARTICLE_INDEX`](@ref).
              + `(halo_idx, subhalo_rel_idx)` -> Sets the position of the potencial minimum for the `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo) as the new origin.
              + `(halo_idx, 0)`               -> Sets the center of mass of the `halo_idx::Int` halo as the new origin.
              + `subhalo_abs_idx`             -> Sets the center of mass of the `subhalo_abs_idx::Int` as the new origin.
          + `:rotation`        -> Rotation for the simulation box. The posibilities are:

              + `:zero`                       -> No rotation is appplied.
              + `:global_am`                  -> Sets the angular momentum of the whole system as the new z axis.
              + `:stellar_am`                 -> Sets the stellar angular momentum as the new z axis.
              + `:stellar_pa`                 -> Sets the stellar principal axis as the new coordinate system.
              + `:stellar_subhalo_pa`         -> Sets the principal axis of the stars in the main subhalo as the new coordinate system.
              + `(halo_idx, subhalo_rel_idx)` -> Sets the principal axis of the stars in `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo), as the new coordinate system.
              + `(halo_idx, 0)`               -> Sets the principal axis of the stars in the `halo_idx::Int` halo, as the new coordinate system.
              + `subhalo_abs_idx`             -> Sets the principal axis of the stars in the `subhalo_abs_idx::Int` subhalo as the new coordinate system.
  - `sim_labels::Union{Vector{String},Nothing}=basename.(simulation_paths)`: Labels for the plot legend, one per simulation. Set it to `nothing` if you don't want a legend.
  - `theme::Attributes=Theme()`: Plot theme that will take precedence over [`DEFAULT_THEME`](@ref).
"""
function stellarCircularity(
    simulation_paths::Vector{String},
    slice::IndexType;
    range::NTuple{2,<:Number}=(-2.0, 2.0),
    n_bins::Int=60,
    output_path::String="./",
    filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all,
    sim_labels::Union{Vector{String},Nothing}=basename.(simulation_paths),
    theme::Attributes=Theme(),
)::Nothing

    plot_params = plotParams(:stellar_circularity)

    filter_function, translation, rotation, request = selectFilter(filter_mode, plot_params.request)

    grid = LinearGrid(range..., n_bins)

    plotSnapshot(
        simulation_paths,
        request,
        [lines!];
        pf_kwargs=[(;)],
        # `plotSnapshot` configuration
        output_path,
        base_filename="circularity_histogram",
        output_format=".png",
        warnings=false,
        show_progress=true,
        # Data manipulation options
        slice,
        filter_function,
        da_functions=[daLineHistogram],
        da_args=[(:stellar_circularity, grid, :stars)],
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
        # Plotting and animation options
        save_figures=true,
        backup_results=false,
        theme=merge(
            theme,
            Theme(Legend=(nbanks=1, halign=:left, valign=:top, padding=(40, 0, 0, 0)),),
        ),
        sim_labels,
        title="",
        colorbar=false,
        # Animation options
        animation=false,
        animation_filename="animation.mp4",
        framerate=10,
    )

    return nothing

end

"""
    compareFeldmann2020(
        simulation_paths::Vector{String},
        x_quantity::Symbol,
        y_quantity::Symbol;
        <keyword arguments>
    )::Nothing

Plot a time series plus the corresponding experimental results from Feldmann (2020).

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`.
  - `x_quantity::Symbol`: Quantity for the x axis. The options are:

      + `:stellar_mass`      -> Stellar mass.
      + `:molecular_mass`    -> Molecular hydrogen (``\\mathrm{H_2}``) mass.
      + `:br_molecular_mass` -> Molecular hydrogen (``\\mathrm{H_2}``) mass, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_mass`       -> Atomic hydrogen (``\\mathrm{HI}``) mass.
      + `:sfr`               -> The star formation rate of the last `AGE_RESOLUTION`.
  - `y_quantity::Symbol`: Quantity for the y axis. The options are:

      + `:stellar_mass`      -> Stellar mass.
      + `:molecular_mass`    -> Molecular hydrogen (``\\mathrm{H_2}``) mass.
      + `:br_molecular_mass` -> Molecular hydrogen (``\\mathrm{H_2}``) mass, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_mass`       -> Atomic hydrogen (``\\mathrm{HI}``) mass.
      + `:sfr`               -> The star formation rate of the last `AGE_RESOLUTION`.
  - `slice::IndexType=(:)`: Slice of the simulations, i.e. which snapshots will be plotted. It can be an integer (a single snapshot), a vector of integers (several snapshots), an `UnitRange` (e.g. 5:13), an `StepRange` (e.g. 5:2:13) or (:) (all snapshots). Starts at 1 and out of bounds indices are ignored.
  - `scatter::Bool=false`: If the data will be presented as a line plot with error bands (default), or alternatively, a scatter plot.
  - `output_path::String="./"`: Path to the output folder.
  - `filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all`: Which cells/particles will be plotted, the options are:

      + `:all`             -> Consider every cell/particle within the simulation box.
      + `:halo`            -> Consider only the cells/particles that belong to the main halo.
      + `:subhalo`         -> Consider only the cells/particles that belong to the main subhalo.
      + `:sphere`          -> Consider only the cell/particle inside a sphere with radius `DISK_R` (see `./src/constants/globals.jl`).
      + `:stellar_subhalo` -> Consider only the cells/particles that belong to the main subhalo.
      + `:all_subhalo`     -> Plot every cell/particle centered around the main subhalo.
      + A dictionary with three entries:

          + `:filter_function` -> The filter function.
          + `:translation`     -> Translation for the simulation box. The posibilities are:

              + `:global_cm`                  -> Selects the center of mass of the whole system as the new origin.
              + `:{component}`                -> Sets the center of mass of the given component (e.g. :stars, :gas, :halo, etc, after filtering) as the new origin. It can be any of the keys of [`PARTICLE_INDEX`](@ref).
              + `(halo_idx, subhalo_rel_idx)` -> Sets the position of the potencial minimum for the `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo) as the new origin.
              + `(halo_idx, 0)`               -> Sets the center of mass of the `halo_idx::Int` halo as the new origin.
              + `subhalo_abs_idx`             -> Sets the center of mass of the `subhalo_abs_idx::Int` as the new origin.
          + `:rotation`        -> Rotation for the simulation box. The posibilities are:

              + `:zero`                       -> No rotation is appplied.
              + `:global_am`                  -> Sets the angular momentum of the whole system as the new z axis.
              + `:stellar_am`                 -> Sets the stellar angular momentum as the new z axis.
              + `:stellar_pa`                 -> Sets the stellar principal axis as the new coordinate system.
              + `:stellar_subhalo_pa`         -> Sets the principal axis of the stars in the main subhalo as the new coordinate system.
              + `(halo_idx, subhalo_rel_idx)` -> Sets the principal axis of the stars in `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo), as the new coordinate system.
              + `(halo_idx, 0)`               -> Sets the principal axis of the stars in the `halo_idx::Int` halo, as the new coordinate system.
              + `subhalo_abs_idx`             -> Sets the principal axis of the stars in the `subhalo_abs_idx::Int` subhalo as the new coordinate system.
  - `sim_labels::Union{Vector{String},Nothing}=basename.(simulation_paths)`: Labels for the plot legend, one per simulation. Set it to `nothing` if you don't want a legend.
  - `theme::Attributes=Theme()`: Plot theme that will take precedence over [`DEFAULT_THEME`](@ref).

# References

R. Feldmann (2020). *The link between star formation and gas in nearby galaxies*. Communications Physics **3(226)**. [doi:10.1038/s42005-020-00493-0](https://doi.org/10.1038/s42005-020-00493-0)
"""
function compareFeldmann2020(
    simulation_paths::Vector{String},
    x_quantity::Symbol,
    y_quantity::Symbol;
    slice::IndexType=(:),
    scatter::Bool=false,
    output_path::String="./",
    filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all,
    sim_labels::Union{Vector{String},Nothing}=basename.(simulation_paths),
    theme::Attributes=Theme(),
)::Nothing

    if x_quantity == :sfr
        x_quantity = :observational_sfr
    end

    if y_quantity == :sfr
        y_quantity = :observational_sfr
    end

    x_plot_params = plotParams(x_quantity)
    y_plot_params = plotParams(y_quantity)

    if x_quantity == :sfr
        x_quantity = :observational_sfr
    end

    plotTimeSeries(
        simulation_paths,
        [scatter!];
        pf_kwargs=[(;)],
        # `plotTimeSeries` configuration
        output_path,
        filename="$(y_quantity)-vs-$(x_quantity)-with-Feldmann2020",
        output_format=".png",
        warnings=false,
        show_progress=true,
        # Data manipulation options
        slice,
        da_functions=[daEvolution],
        da_args=[(x_quantity, y_quantity)],
        da_kwargs=[(; filter_mode, smooth=0, scaling=identity, warnings=false)],
        post_processing=ppFeldmann2020!,
        pp_args=(x_quantity, y_quantity),
        pp_kwargs=(; scatter),
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
        # Plotting options
        save_figures=true,
        backup_results=false,
        theme=merge(theme, Theme(size=(850, 850), Legend=(nbanks=1,))),
        sim_labels,
        title="",
    )

    return nothing

end

"""
    compareMolla2015(
        simulation_paths::Vector{String},
        slice::IndexType,
        quantity::Symbol;
        <keyword arguments>
    )::Nothing

Plot a Milky Way profile plus the corresponding experimental results from Mollá et al. (2015).

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`.
  - `slice::IndexType`: Slice of the simulations, i.e. which snapshots will be plotted. It can be an integer (a single snapshot), a vector of integers (several snapshots), an `UnitRange` (e.g. 5:13), an `StepRange` (e.g. 5:2:13) or (:) (all snapshots). Starts at 1 and out of bounds indices are ignored.
  - `quantity::Symbol`: Quantity for the y axis. The options are:

      + `:stellar_area_density`      -> Stellar area mass density.
      + `:molecular_area_density`    -> Molecular mass surface density.
      + `:br_molecular_area_density` -> Molecular mass surface density, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_area_density`       -> Atomic hydrogen area mass density.
      + `:sfr_area_density`          -> Star formation rate area density, for the last `AGE_RESOLUTION`.
      + `:O_stellar_abundance`       -> Stellar abundance of oxygen, as ``12 + \\log_{10}(\\mathrm{O \\, / \\, H})``.
      + `:N_stellar_abundance`       -> Stellar abundance of nitrogen, as ``12 + \\log_{10}(\\mathrm{N \\, / \\, H})``.
      + `:C_stellar_abundance`       -> Stellar abundance of carbon, as ``12 + \\log_{10}(\\mathrm{C \\, / \\, H})``.
  - `output_path::String="./"`: Path to the output folder.
  - `filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all`: Which cells/particles will be plotted, the options are:

      + `:all`             -> Consider every cell/particle within the simulation box.
      + `:halo`            -> Consider only the cells/particles that belong to the main halo.
      + `:subhalo`         -> Consider only the cells/particles that belong to the main subhalo.
      + `:sphere`          -> Consider only the cell/particle inside a sphere with radius `DISK_R` (see `./src/constants/globals.jl`).
      + `:stellar_subhalo` -> Consider only the cells/particles that belong to the main subhalo.
      + `:all_subhalo`     -> Plot every cell/particle centered around the main subhalo.
      + A dictionary with three entries:

          + `:filter_function` -> The filter function.
          + `:translation`     -> Translation for the simulation box. The posibilities are:

              + `:global_cm`                  -> Selects the center of mass of the whole system as the new origin.
              + `:{component}`                -> Sets the center of mass of the given component (e.g. :stars, :gas, :halo, etc, after filtering) as the new origin. It can be any of the keys of [`PARTICLE_INDEX`](@ref).
              + `(halo_idx, subhalo_rel_idx)` -> Sets the position of the potencial minimum for the `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo) as the new origin.
              + `(halo_idx, 0)`               -> Sets the center of mass of the `halo_idx::Int` halo as the new origin.
              + `subhalo_abs_idx`             -> Sets the center of mass of the `subhalo_abs_idx::Int` as the new origin.
          + `:rotation`        -> Rotation for the simulation box. The posibilities are:

              + `:zero`                       -> No rotation is appplied.
              + `:global_am`                  -> Sets the angular momentum of the whole system as the new z axis.
              + `:stellar_am`                 -> Sets the stellar angular momentum as the new z axis.
              + `:stellar_pa`                 -> Sets the stellar principal axis as the new coordinate system.
              + `:stellar_subhalo_pa`         -> Sets the principal axis of the stars in the main subhalo as the new coordinate system.
              + `(halo_idx, subhalo_rel_idx)` -> Sets the principal axis of the stars in `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo), as the new coordinate system.
              + `(halo_idx, 0)`               -> Sets the principal axis of the stars in the `halo_idx::Int` halo, as the new coordinate system.
              + `subhalo_abs_idx`             -> Sets the principal axis of the stars in the `subhalo_abs_idx::Int` subhalo as the new coordinate system.
  - `sim_labels::Union{Vector{String},Nothing}=basename.(simulation_paths)`: Labels for the plot legend, one per simulation. Set it to `nothing` if you don't want a legend.
  - `theme::Attributes=Theme()`: Plot theme that will take precedence over [`DEFAULT_THEME`](@ref).

# References

M. Mollá et al. (2015). *Galactic chemical evolution: stellar yields and the initial mass function*. Monthly Notices of the Royal Astronomical Society **451(4)**, 3693–3708. [doi:10.1093/mnras/stv1102](https://doi.org/10.1093/mnras/stv1102)
"""
function compareMolla2015(
    simulation_paths::Vector{String},
    slice::IndexType,
    quantity::Symbol;
    output_path::String="./",
    filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all,
    sim_labels::Union{Vector{String},Nothing}=basename.(simulation_paths),
    theme::Attributes=Theme(),
)::Nothing

    plot_params = plotParams(quantity)
    request = addRequest(plot_params.request, Dict(:gas => ["VEL "], :stars => ["VEL "]))
    filter_function, translation, rotation, request = selectFilter(filter_mode, plot_params.request)

    # Select the correct grid acording to the available data from M. Mollá et al. (2015)
    if quantity == :stellar_area_density

        grid = CircularGrid(16.5u"kpc", 14; shift=2.5u"kpc")

    elseif quantity ∈ [:molecular_area_density, :br_molecular_area_density, :sfr_area_density]

        grid = CircularGrid(19.5u"kpc", 20; shift=-0.5u"kpc")

    elseif quantity == :atomic_area_density

        grid = CircularGrid(20.5u"kpc", 21; shift=-0.5u"kpc")

    elseif quantity == :O_stellar_abundance

        grid = CircularGrid(18.5u"kpc", 19; shift=-0.5u"kpc")

    elseif quantity == :N_stellar_abundance

        grid = CircularGrid(17.5u"kpc", 18; shift=-0.5u"kpc")

    elseif quantity == :C_stellar_abundance

        grid = CircularGrid(15.5u"kpc", 16; shift=-0.5u"kpc")

    end

    # Draw the figures with CairoMakie
    plotSnapshot(
        simulation_paths,
        request,
        [scatterlines!];
        pf_kwargs=[(;)],
        # `plotSnapshot` configuration
        output_path,
        base_filename="$(quantity)-profile-with-Molla2015",
        output_format=".png",
        warnings=false,
        show_progress=true,
        # Data manipulation options
        slice,
        filter_function,
        da_functions=[daMolla2015],
        da_args=[(grid, quantity)],
        da_kwargs=[(;)],
        post_processing=ppMolla2015!,
        pp_args=(quantity,),
        pp_kwargs=(; color=Makie.wong_colors()[6], linestyle=nothing, error_bars=true),
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
        # Plotting and animation options
        save_figures=true,
        backup_results=false,
        theme,
        sim_labels,
        title="",
        colorbar=false,
        # Animation options
        animation=false,
        animation_filename="animation.mp4",
        framerate=10,
    )

    return nothing

end

"""
    kennicuttSchmidtLaw(
        simulation_paths::Vector{String},
        slice::IndexType;
        <keyword arguments>
    )::Nothing

Plot the Kennicutt-Schmidt law.

!!! note

    Only stars younger than [`AGE_RESOLUTION`](@ref) and gas cells/particles within a sphere of radius `rmax_gas` are consider. The star formation surface density is just the stellar mass surface density divided by [`AGE_RESOLUTION`](@ref).

!!! note

    This function uses physical units regardless of the global setting [`PHYSICAL_UNITS`](@ref).

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`.
  - `slice::IndexType`: Slice of the simulations, i.e. which snapshots will be plotted. It can be an integer (a single snapshot), a vector of integers (several snapshots), an `UnitRange` (e.g. 5:13), an `StepRange` (e.g. 5:2:13) or (:) (all snapshots). Starts at 1 and out of bounds indices are ignored.
  - `quantity::Symbol=:molecular_mass`: Quantity for the x axis. The options are:

      + `:gas_mass`          -> Total gas mass surface density.
      + `:molecular_mass`    -> Molecular mass surface density. This one can be plotted with the results of Bigiel et al. (2008) and Bigiel et al. (2010).
      + `:br_molecular_mass` -> Molecular mass surface density, computed using the pressure relation in Blitz et al. (2006). This one can be plotted with the results of Bigiel et al. (2008) and Bigiel et al. (2010).
      + `:neutral_mass`      -> Neutral mass surface density. This one can be plotted with the results of Bigiel et al. (2008), Bigiel et al. (2010), and Kennicutt (1998).
  - `type::Symbol=:cells`: If the gas surface density will be calculated assuming the gas is in :particles or in Voronoi :cells.
  - `plot_type::Symbol=:scatter`: If the plot will be a :scatter plot or a :heatmap. Heatmaps will not show legends, experimental measurements or several simulations at once.
  - `integrated::Bool=false`: If the integrated (one point per galaxy) or resolved (several point per galaxy) Kennicutt-Schmidt law will be plotted. `integrated` = true only works with `plot_type` = :scatter, the central value is the weighted median, and the error bars are the median absolute deviations.
  - `sfr_density::Bool=true`: If the quantity for the y axis will be the SFR surface density, if set to false the quantity will be the stellar mass surface density.
  - `gas_weights::Union{Symbol,Nothing}=nothing`: If `plot_type` = :scatter, each point (a pixel of the 2D projected galaxy) can be weighted by a gas quantity. If `integrated` = true, the median will be computed with these weights in mind. If `integrated` = false, each point will have a color given by the weight. The posible weights are:

      + `:gas_mass_density` -> Gas mass surface density of each pixel. See the documentation for the function [`daDensity2DProjection`](@ref).
      + `:gas_sfr`          -> The total gas SFR of the column associated with each pixel. See the documentation for the function [`daGasSFR2DProjection`](@ref).
      + `:gas_metallicity`  -> The total metallicity of the column associated with each pixel. See the documentation for the function [`daMetallicity2DProjection`](@ref).
      + `:temperature`      -> The mean gas temperature of the column associated with each pixel. See the documentation for the function [`daTemperature2DProjection`](@ref).
  - `measurements::Bool=true`: If the experimental measurements from Kennicutt (1998), Bigiel et al. (2008) or Bigiel et al. (2010) will be plotted alongside the simulation results.
  - `measurement_type::Union{String,Symbol}=:fits`: Type of measurement to plot, only valid if `measurement` = true. The option are:

      + `:fits`: Fits from Bigiel et al. (2008) and/or Kennicutt (1998) depending on the quantity in the x axis. The fits will be plotted as lines with uncertanty bands.
      + `"NGC XXX"`: Plot the resolved data of the given NGC galaxy as a scatter plot. Uses the data from Bigiel et al. (2010). See the documentation of [`ppBigiel2010!`](@ref) for options.
      + `:all`: Plot the data of every galaxy in Bigiel et al. (2010), as a scatter plot.
  - `rmax_gas::Unitful.Length=DISK_R`: Maximum radius for the gas cells/particles. Bigiel et al. (2008) uses measurements upto the optical radius r25 (where the B-band magnitude drops below 25 mag arcsec^−2).
  - `reduce_resolution::Bool=true`: If the resolution of the 2D grids will be reduce after the 2D projection to have pixels of a physical size ~ [`BIGIEL_PX_SIZE`](@ref).
  - `x_range::Union{NTuple{2,<:Number},Nothing}=nothing`: x axis range for the heatmap grid. If set to `nothing`, the extrema of the x values will be used. Only relevant if `plot_type` = :heatmap.
  - `y_range::Union{NTuple{2,<:Number},Nothing}=nothing`: y axis range for the heatmap grid. If set to `nothing`, the extrema of the y values will be used. Only relevant if `plot_type` = :heatmap.
  - `n_bins::Int=100`: Number of bins per side of the heatmap grid. Only relevant if `plot_type` = :heatmap.
  - `colorbar::Bool=false`: If a colorbar will be added.
  - `print_range::Bool=false`: Print an info block detailing the color range.
  - `output_file::String="./kennicutt_schmidt_law.png"`: Path to the output file.
  - `filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all`: Which cells/particles will be plotted, the options are:

      + `:all`             -> Consider every cell/particle within the simulation box.
      + `:halo`            -> Consider only the cells/particles that belong to the main halo.
      + `:subhalo`         -> Consider only the cells/particles that belong to the main subhalo.
      + `:sphere`          -> Consider only the cell/particle inside a sphere with radius `DISK_R` (see `./src/constants/globals.jl`).
      + `:stellar_subhalo` -> Consider only the cells/particles that belong to the main subhalo.
      + `:all_subhalo`     -> Plot every cell/particle centered around the main subhalo.
      + A dictionary with three entries:

          + `:filter_function` -> The filter function.
          + `:translation`     -> Translation for the simulation box. The posibilities are:

              + `:global_cm`                  -> Selects the center of mass of the whole system as the new origin.
              + `:{component}`                -> Sets the center of mass of the given component (e.g. :stars, :gas, :halo, etc, after filtering) as the new origin. It can be any of the keys of [`PARTICLE_INDEX`](@ref).
              + `(halo_idx, subhalo_rel_idx)` -> Sets the position of the potencial minimum for the `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo) as the new origin.
              + `(halo_idx, 0)`               -> Sets the center of mass of the `halo_idx::Int` halo as the new origin.
              + `subhalo_abs_idx`             -> Sets the center of mass of the `subhalo_abs_idx::Int` as the new origin.
          + `:rotation`        -> Rotation for the simulation box. The posibilities are:

              + `:zero`                       -> No rotation is appplied.
              + `:global_am`                  -> Sets the angular momentum of the whole system as the new z axis.
              + `:stellar_am`                 -> Sets the stellar angular momentum as the new z axis.
              + `:stellar_pa`                 -> Sets the stellar principal axis as the new coordinate system.
              + `:stellar_subhalo_pa`         -> Sets the principal axis of the stars in the main subhalo as the new coordinate system.
              + `(halo_idx, subhalo_rel_idx)` -> Sets the principal axis of the stars in `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo), as the new coordinate system.
              + `(halo_idx, 0)`               -> Sets the principal axis of the stars in the `halo_idx::Int` halo, as the new coordinate system.
              + `subhalo_abs_idx`             -> Sets the principal axis of the stars in the `subhalo_abs_idx::Int` subhalo as the new coordinate system.
  - `sim_labels::Union{Vector{String},Nothing}=basename.(simulation_paths)`: Labels for the plot legend, one per simulation. Set it to `nothing` if you don't want a legend.
  - `theme::Attributes=Theme()`: Plot theme that will take precedence over [`DEFAULT_THEME`](@ref).

# References

L. Blitz et al. (2006). *The Role of Pressure in GMC Formation II: The H2-Pressure Relation*. The Astrophysical Journal, **650(2)**, 933. [doi:10.1086/505417](https://doi.org/10.1086/505417)

R. C. Kennicutt (1998). *The Global Schmidt Law in Star-forming Galaxies*. The Astrophysical Journal, **498(2)**, 541-552. [doi:10.1086/305588](https://doi.org/10.1086/305588)

F. Bigiel et al. (2008). *THE STAR FORMATION LAW IN NEARBY GALAXIES ON SUB-KPC SCALES*. The Astrophysical Journal, **136(6)**, 2846. [doi:10.1088/0004-6256/136/6/2846](https://doi.org/10.1088/0004-6256/136/6/2846)

F. Bigiel et al. (2010). *EXTREMELY INEFFICIENT STAR FORMATION IN THE OUTER DISKS OF NEARBY GALAXIES*. The Astrophysical Journal, **140(5)**, 1194. [doi:10.1088/0004-6256/140/5/1194](https://doi.org/10.1088/0004-6256/140/5/1194)
"""
function kennicuttSchmidtLaw(
    simulation_paths::Vector{String},
    slice::IndexType;
    quantity::Symbol=:molecular_mass,
    type::Symbol=:cells,
    plot_type::Symbol=:scatter,
    integrated::Bool=false,
    sfr_density::Bool=true,
    gas_weights::Union{Symbol,Nothing}=nothing,
    measurements::Bool=true,
    measurement_type::Union{String,Symbol}=:fits,
    rmax_gas::Unitful.Length=DISK_R,
    reduce_resolution::Bool=true,
    x_range::Union{NTuple{2,<:Number},Nothing}=nothing,
    y_range::Union{NTuple{2,<:Number},Nothing}=nothing,
    n_bins::Int=100,
    colorbar::Bool=false,
    print_range::Bool=false,
    output_file::String="./kennicutt_schmidt_law.png",
    filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all,
    sim_labels::Union{Vector{String},Nothing}=basename.(simulation_paths),
    warnings::Bool=true,
    theme::Attributes=Theme(),
)::Nothing

    # Save the origial value of the global `PHYSICAL_UNITS`
    og_pu_value = PHYSICAL_UNITS

    if !og_pu_value && warnings
        @warn("kennicuttSchmidtLaw: The global `PHYSICAL_UNITS` is set to false, \
        but Kennicutt-Schmidt law plots must be in physical units, so the global \
        setting will be ignored and default to true.")
    end

    # Kennicutt-Schmidt law plots must be in physical units even for cosmological simulations
    global PHYSICAL_UNITS = true

    ns = length(simulation_paths)

    ################################################################################################
    # Check arguments
    ################################################################################################

    (
        quantity ∈ [:gas_mass, :molecular_mass, :br_molecular_mass, :neutral_mass] ||
        throw(ArgumentError("kennicuttSchmidtLaw: `quantity` can only be :gas_mass, \
        :molecular_mass, :br_molecular_mass or :neutral_mass, but I got :$(quantity)"))
    )

    (
        plot_type ∈ [:scatter, :heatmap] ||
        throw(ArgumentError("kennicuttSchmidtLaw: `plot_type` can only be :scatter or :heatmap, \
        but I got :$(plot_type)"))
    )

    if integrated

        if plot_type == :heatmap

            !warnings || @warn("kennicuttSchmidtLaw: If `integrated` is set to true, \
            `plot_type` = :heatmap will be ignored and default to :scatter.")

            plot_type = :scatter

        end

        if reduce_resolution
            !warnings || @warn("kennicuttSchmidtLaw: `integrated` and `reduce_resolution` are set \
            to true. Are you sure you want this?")
        end

    end

    if !isnothing(sim_labels)

        nl = length(sim_labels)

        (
            ns == nl || throw(ArgumentError("kennicuttSchmidtLaw: `sim_labels` must have as many
            elements as `simulation_paths`, but I got length(simulation_paths) = $(ns) \
            != length(sim_labels) = $(nl)"))
        )

    end

    if measurements

        if !sfr_density

            !warnings || @warn("kennicuttSchmidtLaw: If `sfr_density` = false, \
            `measurements` = true will be ignored and default to false. The experimental \
            measurements are only for the SFR surface density.")

            measurements = false

        end

        if quantity == :gas_mass

            !warnings || @warn("kennicuttSchmidtLaw: There are no measurements \
            (fits or otherwise) for `quantity` = :gas_mass. `measurements` = true will be \
            ignored and default to false.")

            measurements = false

        end

        if isa(measurement_type, String) && integrated

            !warnings || @warn("kennicuttSchmidtLaw: `integrated` = true but you have set \
            `measurement_type` to plot the resolved measurement of galaxy $(measurement_type). \
            Are you sure you want this?")

        end

        if measurement_type == :all && integrated

            !warnings || @warn("kennicuttSchmidtLaw: `integrated` = true but you have set \
            `measurement_type` to plot the resolved measurements of all galaxies in \
            Bigiel et al. (2010). Are you sure you want this?")

        end

    end

    if plot_type == :heatmap

        if !isnothing(gas_weights)

            !warnings || @warn("kennicuttSchmidtLaw: If `plot_type` = :heatmap, \
            `gas_weights` = :$(gas_weights) will be ignored and default to nothing.")

            gas_weights = nothing

        end

        if measurements

            !warnings || @warn("kennicuttSchmidtLaw: If `plot_type` = :heatmap, \
            `measurements` = true will be ignored and default to false.")

            measurements = false

        end

        if ns > 1

            !warnings || @warn("kennicuttSchmidtLaw: If `plot_type` = :heatmap, only one \
            simulation at a time can be plotted, but I got length(simulation_paths) = $(ns) > 1. \
            `plot_type` = :heatmap will be ignored and default to :scatter.")

            plot_type = :scatter

        end

    end

    if colorbar && ((plot_type == :scatter && isnothing(gas_weights)) || integrated)

        (
            !warnings || @warn("kennicuttSchmidtLaw: `colorbar` is set to true, \
            but there is no color range in the plot (either `plot_type` = :scatter and \
            `gas_weights` = nothing or `integrated` = true). `colorbar` = true will be \
            ignored and default to false")
        )

        colorbar = false

    end

    ################################################################################################
    # Compute grids
    ################################################################################################

    if reduce_resolution

        # Compute the number of bins in the low resolution grid (pixel size of ~ BIGIEL_PX_SIZE)
        lr_n_bins = round(Int, uconvert(Unitful.NoUnits, BOX_L / BIGIEL_PX_SIZE))

        # Compute the interger factor between the high resolution grid (~ 400px)
        # and the low resolution grid (`lr_n_bins`px)
        reduce = 400 ÷ lr_n_bins

        stellar_grid = CubicGrid(BOX_L, reduce * lr_n_bins)
        gas_grid     = CubicGrid(BOX_L, reduce * lr_n_bins)

    else

        reduce = 1

        stellar_grid = CubicGrid(BOX_L, 400)
        gas_grid     = CubicGrid(BOX_L, 400)

    end

    ################################################################################################
    # Compute the density maps and save them as JLD2 files
    ################################################################################################

    # Set a folder for the JLD2 files
    temp_folder = joinpath(dirname(output_file), "_temp_jld2")

    filter_function, translation, rotation, request = selectFilter(
        filter_mode,
        plotParams(:stellar_mass).request,
    )

    # Compute the stellar map
    plotSnapshot(
        simulation_paths,
        request,
        [heatmap!];
        output_path=temp_folder,
        base_filename=string(:stellar_mass),
        warnings,
        slice,
        filter_function,
        da_functions=[daDensity2DProjection],
        da_args=[(stellar_grid, :stellar_mass, :particles)],
        da_kwargs=[(; reduce, filter_function=dd->filterStellarAge(dd; age=AGE_RESOLUTION))],
        transform_box=true,
        translation,
        rotation,
        x_unit=u"kpc",
        y_unit=u"kpc",
        save_figures=false,
        backup_results=true,
    )

    filter_function, translation, rotation, request = selectFilter(
        filter_mode,
        plotParams(quantity).request,
    )

    # Compute the `quantity` map
    plotSnapshot(
        simulation_paths,
        request,
        [heatmap!];
        output_path=temp_folder,
        base_filename=string(quantity),
        warnings,
        slice,
        filter_function,
        da_functions=[daDensity2DProjection],
        da_args=[(gas_grid, quantity, type)],
        da_kwargs=[
            (; reduce, filter_function=dd->filterWithinSphere(dd, (0.0u"kpc", rmax_gas), :zero)),
        ],
        transform_box=true,
        translation,
        rotation,
        x_unit=u"kpc",
        y_unit=u"kpc",
        save_figures=false,
        backup_results=true,
    )

    # Compute the weights map
    if !isnothing(gas_weights)

        if gas_weights == :gas_mass_density

            da_function = daDensity2DProjection
            da_args     = [(gas_grid, :gas_mass, type)]
            c_label     = L"\log_{10} \, \Sigma_\mathrm{gas} \,\, [\mathrm{M_\odot \, kpc^{-2}}]"

        elseif gas_weights == :gas_sfr

            da_function = daGasSFR2DProjection
            da_args     = [(gas_grid, type)]
            c_label     = L"\log_{10} \, \mathrm{SFR_{gas} \,\, [M_\odot \, yr^{-1}]}"

        elseif gas_weights == :gas_metallicity

            da_function = daMetallicity2DProjection
            da_args     = [(gas_grid, :gas, type)]
            c_label     = L"$\log_{10}$ %$(plotParams(:gas_metallicity).var_name)"

        elseif gas_weights == :temperature

            da_function = daTemperature2DProjection
            da_args     = [(gas_grid, type)]
            c_label     = plotParams(:temperature).axis_label

        else

            throw(ArgumentError("kennicuttSchmidtLaw: `gas_weights` can only be \
            :gas_mass_density, :gas_sfr, :gas_metallicity or :temperature, but I got \
            :$(gas_weights)"))

        end

        filter_function, translation, rotation, request = selectFilter(
            filter_mode,
            mergeRequests(plotParams(gas_weights).request, plotParams(:gas_mass_density).request),
        )

        plotSnapshot(
            simulation_paths,
            request,
            [heatmap!];
            output_path=temp_folder,
            base_filename="gas_weights",
            warnings,
            slice,
            filter_function,
            da_functions=[da_function],
            da_args,
            da_kwargs=[
                (;
                    reduce,
                    filter_function=dd->filterWithinSphere(dd, (0.0u"kpc", rmax_gas), :zero),
                    print_range,
                ),
            ],
            transform_box=true,
            translation,
            rotation,
            x_unit=u"kpc",
            y_unit=u"kpc",
            save_figures=false,
            backup_results=true,
        )

    end

    ################################################################################################
    # Set the axis labels
    ################################################################################################

    if quantity == :gas_mass

        x_label = getLabel(plotParams(:gas_area_density).var_name, 0, u"Msun * kpc^-2")

    elseif quantity == :molecular_mass

        x_label = getLabel(plotParams(:molecular_area_density).var_name, 0, u"Msun * kpc^-2")

    elseif quantity == :br_molecular_mass

        x_label = getLabel(plotParams(:br_molecular_area_density).var_name, 0, u"Msun * kpc^-2")

    elseif quantity == :neutral_mass

        x_label = getLabel(plotParams(:neutral_area_density).var_name, 0, u"Msun * kpc^-2")

    end

    if sfr_density

        y_label = getLabel(plotParams(:sfr_area_density).var_name, 0, u"Msun * yr^-1 * kpc^-2")

    else

        y_label = getLabel(plotParams(:stellar_area_density).var_name, 0, u"Msun * kpc^-2")

    end

    ################################################################################################
    # Read and plot the data
    ################################################################################################

    if sfr_density
        # Factor to go from stellar surface density to SFR surface density
        # log10(Σsfr) = log10(Σ*) - log10Δt
        log10Δt = log10(ustrip(u"yr", AGE_RESOLUTION))
    end

    # Set the plot theme
    current_theme = merge(
        theme,
        Theme(
            Legend=(
                nbanks=1,
                labelsize=25,
                rowgap=-15,
                halign=:left,
                valign=:top,
                padding=(15, 0, 0, 0),
            ),
        ),
        DEFAULT_THEME,
        theme_latexfonts(),
    )

    with_theme(current_theme) do

        f = Figure()

        ax = CairoMakie.Axis(
            f[1, 1];
            xlabel=L"$\log_{10}$ %$(x_label)",
            ylabel=L"$\log_{10}$ %$(y_label)",
            aspect=AxisAspect(1),
        )

        colors = [:grey25, current_theme[:palette][:color][][2:ns]...]

        for (sim_idx, simulation) in pairs(simulation_paths)

            simulation_table = DataFrame(makeSimulationTable(simulation; warnings)[slice, :])
            sim_name         = "simulation_$(lpad(string(sim_idx), 3, "0"))"
            snapshot_numbers = simulation_table[!, :numbers]

            # Allocate memory for the heatmap. For a heatmap we need to accumulate the values
            # for every snapshot before plotting
            if plot_type == :heatmap
                x_heatmap = Float64[]
                y_heatmap = Float64[]
            end

            for snapshot_number in snapshot_numbers

                ####################################################################################
                # Read the JLD2 files and sanitize the data
                ####################################################################################

                x_address = "$(quantity)-$(SNAP_BASENAME)_$(snapshot_number)/$(sim_name)"
                x_file    = jldopen(joinpath(temp_folder, "$(string(quantity)).jld2"), "r")
                x_data    = vec(x_file[x_address][3])
                x_idxs    = map(x -> isnan(x) || iszero(x), x_data)

                y_address = "stellar_mass-$(SNAP_BASENAME)_$(snapshot_number)/$(sim_name)"
                y_file    = jldopen(joinpath(temp_folder, "stellar_mass.jld2"), "r")
                y_data    = vec(y_file[y_address][3])
                y_idxs    = map(x -> isnan(x) || iszero(x), y_data)

                delete_idxs = x_idxs ∪ y_idxs

                if !isnothing(gas_weights)

                    z_address = "gas_weights-$(SNAP_BASENAME)_$(snapshot_number)/$(sim_name)"
                    z_file    = jldopen(joinpath(temp_folder, "gas_weights.jld2"), "r")
                    z_data    = vec(z_file[z_address][3])
                    z_idxs    = map(x -> isnan(x) || iszero(x), z_data)

                    delete_idxs = delete_idxs ∪ z_idxs

                    deleteat!(z_data, delete_idxs)

                end

                deleteat!(x_data, delete_idxs)
                deleteat!(y_data, delete_idxs)

                if sfr_density
                    y_data .-= log10Δt
                end

                # For the integrated Kennicutt-Schmidt law, compute the gas and stellar
                # densities median and median absolute deviation
                if integrated

                    lin_x = exp10.(x_data)
                    lin_y = exp10.(y_data)

                    if !isnothing(gas_weights)

                        w_z = weights(exp10.(z_data))

                        μx = median(lin_x, w_z)
                        μy = median(lin_y, w_z)

                    else

                        μx = median(lin_x)
                        μy = median(lin_y)

                    end

                    σx = mad(lin_x; center=μx, normalize=false)
                    σy = mad(lin_y; center=μy, normalize=false)

                    x_data = [log10(μx ± σx)]
                    y_data = [log10(μy ± σy)]

                end

                if plot_type == :heatmap

                    append!(x_heatmap, x_data)
                    append!(y_heatmap, y_data)

                end

                close(x_file)
                close(y_file)
                isnothing(gas_weights) || close(z_file)

                ####################################################################################
                # Plot the Kennicutt-Schmidt law
                ####################################################################################

                if  plot_type == :scatter

                    if integrated

                        scatter!(
                            ax,
                            Measurements.value.(x_data),
                            Measurements.value.(y_data);
                            color=colors[sim_idx],
                            markersize=15,
                        )

                        errorbars!(
                            ax,
                            Measurements.value.(x_data),
                            Measurements.value.(y_data),
                            Measurements.uncertainty.(y_data);
                            color=colors[sim_idx],
                            direction=:y,
                        )

                        errorbars!(
                            ax,
                            Measurements.value.(x_data),
                            Measurements.value.(y_data),
                            Measurements.uncertainty.(x_data);
                            color=colors[sim_idx],
                            direction=:x,
                        )

                    else

                        if isnothing(gas_weights)

                            scatter!(
                                ax,
                                x_data,
                                y_data;
                                color=(colors[sim_idx], 0.5),
                                markersize=6,
                                alpha=1.0,
                            )

                        else

                            scatter!(
                                ax,
                                x_data,
                                y_data;
                                markersize=6,
                                color=z_data,
                                colormap=:nipy_spectral,
                            )

                        end

                    end

                end

            end

            if plot_type == :heatmap

                # If there is no specified range, use the extrema of the x values
                if isnothing(x_range)
                    xrange = extrema(x_heatmap)
                else
                    xrange = x_range
                end

                # If there is no specified range, use the extrema of the y values
                if isnothing(y_range)
                    yrange = extrema(y_heatmap)
                else
                    yrange = y_range
                end

                # Compute the bin half width for each axis
                x_bin_h_width = 0.5 * (xrange[2] - xrange[1]) / n_bins
                y_bin_h_width = 0.5 * (yrange[2] - yrange[1]) / n_bins

                # Compute the center value of each bin for each axis
                x_axis = collect(
                    range(
                        xrange[1] + x_bin_h_width;
                        length=n_bins,
                        step=2 * x_bin_h_width,
                    ),
                )
                y_axis = collect(
                    range(
                        yrange[1] + y_bin_h_width;
                        length=n_bins,
                        step=2 * y_bin_h_width,
                    ),
                )

                # Compute the 2D histogram (number of pixels in each bin)
                values = histogram2D(
                    permutedims(hcat(x_heatmap, y_heatmap), (2, 1)),
                    collect(range(xrange[1], xrange[2]; length=n_bins + 1)),
                    collect(range(yrange[1], yrange[2]; length=n_bins + 1));
                )

                # The transpose and reverse operation are to conform to the way heatmap!
                # expect the matrix to be structured
                z_axis = reverse!(transpose(values), dims=2)

                heatmap!(ax, x_axis, y_axis, z_axis; colormap=:nipy_spectral)

                if print_range

                    # Compute the mininimum and maximum of `z_axis`
                    min_max = isempty(z_axis) ? (NaN, NaN) : extrema(filter(!isnan, z_axis))

                    # Print the count range
                    @info(
                        "\nCount range \
                        \n  Simulation: $(basename(simulation)) \
                        \n  Quantity:   $(quantity) \
                        \n  Counts:     $(min_max)\n\n"
                    )

                end

            end

        end

        ############################################################################################
        # Plot the experimental fits and the legend
        ############################################################################################

        if !isnothing(sim_labels) && plot_type == :scatter

            if !isnothing(gas_weights) && !integrated

                markers = [
                    MarkerElement(; color=:grey25, marker=:circle, markersize=20) for
                    _ in eachindex(sim_labels)
                ]

            else

                markers = [
                    MarkerElement(; color, marker=:circle, markersize=20) for color in colors
                ]

            end

        end

        if measurements

            if measurement_type == :fits

                if quantity ∈ [:molecular_mass, :br_molecular_mass]

                    pp_legend = ppBigiel2008!(
                        f,
                        true;
                        x_unit=u"Msun * kpc^-2",
                        colors=[Makie.wong_colors()[1], Makie.wong_colors()[2]],
                        warnings,
                    )

                elseif quantity == :neutral_mass

                    legend_bigiel = ppBigiel2008!(
                        f,
                        false;
                        x_unit=u"Msun * kpc^-2",
                        colors=[Makie.wong_colors()[1], Makie.wong_colors()[2]],
                        warnings,
                    )

                    legend_kennicut = ppKennicutt1998!(
                        f;
                        x_unit=u"Msun * kpc^-2",
                        colors=[Makie.wong_colors()[3], Makie.wong_colors()[4]],
                        warnings,
                    )

                    pp_legend = (
                        vcat(legend_bigiel[1], legend_kennicut[1]),
                        vcat(legend_bigiel[2], legend_kennicut[2]),
                    )

                end

            else

                if quantity ∈ [:molecular_mass, :br_molecular_mass]

                    pp_legend = ppBigiel2010!(
                        f;
                        galaxy=measurement_type,
                        quantity=:molecular,
                        x_unit=u"Msun * kpc^-2",
                    )

                elseif quantity == :neutral_mass

                    pp_legend = ppBigiel2010!(
                        f;
                        galaxy=measurement_type,
                        quantity=:neutral,
                        x_unit=u"Msun * kpc^-2",
                    )

                end

            end

            if !isnothing(sim_labels) && plot_type == :scatter

                Makie.Legend(f[1, 1], vcat(markers, pp_legend[1]), vcat(sim_labels, pp_legend[2]))

            end

        else

            if !isnothing(sim_labels) && plot_type == :scatter

                Makie.Legend(f[1, 1], markers, sim_labels)

            end

        end

        if !isnothing(sim_labels) && plot_type == :heatmap

            ppAnnotation!(f, sim_labels[1]; color=:white, fontsize=30)

        end

        ############################################################################################
        # Print the colorbar
        ############################################################################################

        if colorbar

            if plot_type == :heatmap
                Colorbar(f[1, 2], f.content[1].scene.plots[1]; label=L"\mathrm{Counts}")
            else
                Colorbar(f[1, 2], f.content[1].scene.plots[1]; label=c_label)
            end

            rowsize!(f.layout, 1, Makie.Fixed(pixelarea(ax.scene)[].widths[2]))

        end

        ############################################################################################
        # Save the plot
        ############################################################################################

        Makie.save(output_file, f)

    end

    rm(temp_folder; recursive=true)

    global PHYSICAL_UNITS = og_pu_value

    return nothing

end

"""
    fitResolvedKSLaw(
        simulation_path::String,
        slice::IndexType;
        <keyword arguments>
    )::Nothing

Plot the resolved Kennicutt-Schmidt relation with an optional linear fit.

!!! note

    Only stars younger than [`AGE_RESOLUTION`](@ref) and gas cells/particles within a sphere of radius `rmax_gas` are consider. The star formation surface density is just the stellar mass surface density divided by [`AGE_RESOLUTION`](@ref).

# Arguments

  - `simulation_path::String`: Path to the simulation directory, set in the code variable `OutputDir`.
  - `slice::IndexType`: Slice of the simulations, i.e. which snapshots will be plotted. It can be an integer (a single snapshot), a vector of integers (several snapshots), an `UnitRange` (e.g. 5:13), an `StepRange` (e.g. 5:2:13) or (:) (all snapshots). Starts at 1 and out of bounds indices are ignored.
  - `quantity::Symbol=:molecular_mass`: Quantity for the x axis. The options are:

      + `:gas_mass`          -> Gas mass surface density.
      + `:molecular_mass`    -> Molecular mass surface density.
      + `:br_molecular_mass` -> Molecular mass surface density, computed using the pressure relation in Blitz et al. (2006).
      + `:neutral_mass`      -> Neutral mass surface density.
  - `type::Symbol=:cells`: If the gas surface density will be calculated assuming the gas is in `:particles` or in Voronoi `:cells`.
  - `fit::Bool=true`: If a fit of the plotted values will be added on top of the scatter plot.
  - `rmax_gas::Unitful.Length=DISK_R`: Maximum radius for the gas cells/particles. Bigiel et al. (2008) uses measurements upto the optical radius r25 (where the B-band magnitude drops below 25 mag arcsec^−2).
  - `x_range::NTuple{2,<:Real}=(-Inf, Inf)`: Only the data withing this range (for the x coordinates) will be fitted.
  - `output_path::String="./"`: Path to the output folder.
  - `filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all`: Which cells/particles will be plotted, the options are:

      + `:all`             -> Consider every cell/particle within the simulation box.
      + `:halo`            -> Consider only the cells/particles that belong to the main halo.
      + `:subhalo`         -> Consider only the cells/particles that belong to the main subhalo.
      + `:sphere`          -> Consider only the cell/particle inside a sphere with radius `DISK_R` (see `./src/constants/globals.jl`).
      + `:stellar_subhalo` -> Consider only the cells/particles that belong to the main subhalo.
      + `:all_subhalo`     -> Plot every cell/particle centered around the main subhalo.
      + A dictionary with three entries:

          + `:filter_function` -> The filter function.
          + `:translation`     -> Translation for the simulation box. The posibilities are:

              + `:global_cm`                  -> Selects the center of mass of the whole system as the new origin.
              + `:{component}`                -> Sets the center of mass of the given component (e.g. :stars, :gas, :halo, etc, after filtering) as the new origin. It can be any of the keys of [`PARTICLE_INDEX`](@ref).
              + `(halo_idx, subhalo_rel_idx)` -> Sets the position of the potencial minimum for the `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo) as the new origin.
              + `(halo_idx, 0)`               -> Sets the center of mass of the `halo_idx::Int` halo as the new origin.
              + `subhalo_abs_idx`             -> Sets the center of mass of the `subhalo_abs_idx::Int` as the new origin.
          + `:rotation`        -> Rotation for the simulation box. The posibilities are:

              + `:zero`                       -> No rotation is appplied.
              + `:global_am`                  -> Sets the angular momentum of the whole system as the new z axis.
              + `:stellar_am`                 -> Sets the stellar angular momentum as the new z axis.
              + `:stellar_pa`                 -> Sets the stellar principal axis as the new coordinate system.
              + `:stellar_subhalo_pa`         -> Sets the principal axis of the stars in the main subhalo as the new coordinate system.
              + `(halo_idx, subhalo_rel_idx)` -> Sets the principal axis of the stars in `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo), as the new coordinate system.
              + `(halo_idx, 0)`               -> Sets the principal axis of the stars in the `halo_idx::Int` halo, as the new coordinate system.
              + `subhalo_abs_idx`             -> Sets the principal axis of the stars in the `subhalo_abs_idx::Int` subhalo as the new coordinate system.
  - `sim_label::Union{String,Nothing}=basename(simulation_path)`: Label for the scatter plot. Set it to `nothing` if you don't want a legend.
  - `theme::Attributes=Theme()`: Plot theme that will take precedence over [`DEFAULT_THEME`](@ref).

# References

R. C. Kennicutt (1998). *The Global Schmidt Law in Star-forming Galaxies*. The Astrophysical Journal, **498(2)**, 541-552. [doi:10.1086/305588](https://doi.org/10.1086/305588)

F. Bigiel et al. (2008). *THE STAR FORMATION LAW IN NEARBY GALAXIES ON SUB-KPC SCALES*. The Astrophysical Journal, **136(6)**, 2846. [doi:10.1088/0004-6256/136/6/2846](https://doi.org/10.1088/0004-6256/136/6/2846)
"""
function fitResolvedKSLaw(
    simulation_path::String,
    slice::IndexType;
    quantity::Symbol=:molecular_mass,
    type::Symbol=:cells,
    fit::Bool=true,
    rmax_gas::Unitful.Length=DISK_R,
    x_range::NTuple{2,<:Real}=(-Inf, Inf),
    output_path::String="./",
    filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all,
    sim_label::Union{String,Nothing}=basename(simulation_path),
    theme::Attributes=Theme(),
)::Nothing

    # Compute the number of bins in the low resolution grid (pixel size of ~ BIGIEL_PX_SIZE)
    lr_n_bins = round(Int, uconvert(Unitful.NoUnits, BOX_L / BIGIEL_PX_SIZE))

    # Compute the interger factor between the high resolution grid (~ 400px)
    # and the low resolution grid (`lr_n_bins`px)
    factor = 400 ÷ lr_n_bins

    grid = CubicGrid(BOX_L, factor * lr_n_bins)

    filter_function, translation, rotation, request = selectFilter(
        filter_mode,
        mergeRequests(plotParams(quantity).request, plotParams(:stellar_mass).request),
    )

    # Choose the correct x label
    if quantity == :gas_mass

        x_label = getLabel(
            plotParams(:gas_area_density).var_name,
            0,
            u"Msun * kpc^-2";
            latex=true,
        )

    elseif quantity == :molecular_mass

        x_label = getLabel(
            plotParams(:molecular_area_density).var_name,
            0,
            u"Msun * kpc^-2";
            latex=true,
        )

    elseif quantity == :br_molecular_mass

        x_label = getLabel(
            plotParams(:br_molecular_area_density).var_name,
            0,
            u"Msun * kpc^-2";
            latex=true,
        )

    elseif quantity == :neutral_mass

        x_label = getLabel(
            plotParams(:neutral_area_density).var_name,
            0,
            u"Msun * kpc^-2";
            latex=true,
        )

    end

    # Set the y label
    y_label = getLabel(
        plotParams(:sfr_area_density).var_name,
        0,
        u"Msun * yr^-1 * kpc^-2";
        latex=true,
    )

    plotSnapshot(
        [simulation_path],
        request,
        [scatter!];
        pf_kwargs=[(; color=Makie.wong_colors()[1], markersize=6, marker=:circle)],
        # `plotSnapshot` configuration
        output_path,
        base_filename="ks_law",
        output_format=".png",
        warnings=false,
        show_progress=true,
        # Data manipulation options
        slice,
        filter_function,
        da_functions=[daKennicuttSchmidtLaw],
        da_args=[(grid, quantity)],
        da_kwargs=[
            (;
                type,
                reduce_factor=factor,
                stellar_ff=dd->filterStellarAge(dd),
                gas_ff=dd->filterWithinSphere(dd, (0.0u"kpc", rmax_gas), :zero),
            ),
        ],
        post_processing=fit ? ppFitLine! : getNothing,
        pp_args=(),
        pp_kwargs=(;),
        transform_box=true,
        translation,
        rotation,
        smooth=0,
        x_unit=Unitful.NoUnits,
        y_unit=Unitful.NoUnits,
        x_exp_factor=0,
        y_exp_factor=0,
        x_trim=x_range,
        y_trim=(-Inf, Inf),
        x_edges=false,
        y_edges=false,
        x_func=identity,
        y_func=identity,
        # Axes options
        xaxis_label=L"$\log_{10}$ %$(x_label)",
        yaxis_label=L"$\log_{10}$ %$(y_label)",
        xaxis_var_name="",
        yaxis_var_name="",
        xaxis_scale_func=identity,
        yaxis_scale_func=identity,
        # Plotting and animation options
        save_figures=true,
        backup_results=false,
        theme=merge(
            theme,
            Theme(Legend=(labelsize=20, halign=:left, valign=:top, padding=(15, 0, 0, 125)),),
        ),
        sim_labels=[sim_label],
        title="",
        colorbar=false,
        # Animation options
        animation=false,
        animation_filename="animation.mp4",
        framerate=10,
    )

    return nothing

end

"""
    fitVSFLaw(
        simulation_path::String,
        slice::IndexType;
        <keyword arguments>
    )::Nothing

Plot the resolved volumetric star formation (VSF) law with an optional linear fit.

!!! note

    Only stars younger than [`AGE_RESOLUTION`](@ref) and gas cells/particles within a sphere of radius `rmax_gas` are consider. The star formation surface density is just the stellar mass surface density divided by [`AGE_RESOLUTION`](@ref).

# Arguments

  - `simulation_path::String`: Path to the simulation directory, set in the code variable `OutputDir`.
  - `slice::IndexType`: Slice of the simulations, i.e. which snapshots will be plotted. It can be an integer (a single snapshot), a vector of integers (several snapshots), an `UnitRange` (e.g. 5:13), an `StepRange` (e.g. 5:2:13) or (:) (all snapshots). Starts at 1 and out of bounds indices are ignored.
  - `quantity::Symbol`: Quantity for the x axis. The options are:

      + `:gas_mass`          -> Gas density.
      + `:hydrogen_mass`     -> Hydrogen density.
      + `:molecular_mass`    -> Molecular hydrogen (``\\mathrm{H_2}``) density.
      + `:br_molecular_mass` -> Molecular hydrogen (``\\mathrm{H_2}``) density, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_mass`       -> Atomic hydrogen (``\\mathrm{HI}``) density.
      + `:ionized_mass`      -> Ionized hydrogen (``\\mathrm{HII}``) density.
      + `:neutral_mass`      -> Neutral hydrogen (``\\mathrm{HI + H_2}``) density.
  - `type::Symbol=:cells`: If the gas surface density will be calculated assuming the gas is in `:particles` or in Voronoi `:cells`.
  - `fit::Bool=true`: If a fit of the plotted values will be added on top of the scatter plot.
  - `rmax_gas::Unitful.Length=DISK_R`: Maximum radius for the gas cells/particles.
  - `x_range::NTuple{2,<:Real}=(-Inf, Inf)`: Only the data withing this range (for the x coordinates) will be fitted.
  - `output_path::String="./"`: Path to the output folder.
  - `filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all`: Which cells/particles will be plotted, the options are:

      + `:all`             -> Consider every cell/particle within the simulation box.
      + `:halo`            -> Consider only the cells/particles that belong to the main halo.
      + `:subhalo`         -> Consider only the cells/particles that belong to the main subhalo.
      + `:sphere`          -> Consider only the cell/particle inside a sphere with radius `DISK_R` (see `./src/constants/globals.jl`).
      + `:stellar_subhalo` -> Consider only the cells/particles that belong to the main subhalo.
      + `:all_subhalo`     -> Plot every cell/particle centered around the main subhalo.
      + A dictionary with three entries:

          + `:filter_function` -> The filter function.
          + `:translation`     -> Translation for the simulation box. The posibilities are:

              + `:global_cm`                  -> Selects the center of mass of the whole system as the new origin.
              + `:{component}`                -> Sets the center of mass of the given component (e.g. :stars, :gas, :halo, etc, after filtering) as the new origin. It can be any of the keys of [`PARTICLE_INDEX`](@ref).
              + `(halo_idx, subhalo_rel_idx)` -> Sets the position of the potencial minimum for the `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo) as the new origin.
              + `(halo_idx, 0)`               -> Sets the center of mass of the `halo_idx::Int` halo as the new origin.
              + `subhalo_abs_idx`             -> Sets the center of mass of the `subhalo_abs_idx::Int` as the new origin.
          + `:rotation`        -> Rotation for the simulation box. The posibilities are:

              + `:zero`                       -> No rotation is appplied.
              + `:global_am`                  -> Sets the angular momentum of the whole system as the new z axis.
              + `:stellar_am`                 -> Sets the stellar angular momentum as the new z axis.
              + `:stellar_pa`                 -> Sets the stellar principal axis as the new coordinate system.
              + `:stellar_subhalo_pa`         -> Sets the principal axis of the stars in the main subhalo as the new coordinate system.
              + `(halo_idx, subhalo_rel_idx)` -> Sets the principal axis of the stars in `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo), as the new coordinate system.
              + `(halo_idx, 0)`               -> Sets the principal axis of the stars in the `halo_idx::Int` halo, as the new coordinate system.
              + `subhalo_abs_idx`             -> Sets the principal axis of the stars in the `subhalo_abs_idx::Int` subhalo as the new coordinate system.
  - `sim_label::Union{String,Nothing}=basename(simulation_path)`: Label for the scatter plot. Set it to `nothing` if you don't want a legend.
  - `theme::Attributes=Theme()`: Plot theme that will take precedence over [`DEFAULT_THEME`](@ref).
"""
function fitVSFLaw(
    simulation_path::String,
    slice::IndexType;
    quantity::Symbol=:molecular_mass,
    type::Symbol=:cells,
    fit::Bool=true,
    rmax_gas::Unitful.Length=DISK_R,
    x_range::NTuple{2,<:Real}=(-Inf, Inf),
    output_path::String="./",
    filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all,
    sim_label::Union{String,Nothing}=basename(simulation_path),
    theme::Attributes=Theme(),
)::Nothing

    grid = CubicGrid(BOX_L, 400)

    filter_function, translation, rotation, request = selectFilter(
        filter_mode,
        mergeRequests(plotParams(quantity).request, plotParams(:stellar_mass).request),
    )

    # Choose the correct x label
    if quantity == :gas_mass

        x_label = getLabel(
            plotParams(:gas_mass_density).var_name,
            0,
            u"Msun * pc^-3";
            latex=true,
        )

    elseif quantity == :molecular_mass

        x_label = getLabel(
            L"\rho_\mathrm{H_2}",
            0,
            u"Msun * pc^-3";
            latex=true,
        )

    elseif quantity == :br_molecular_mass

        x_label = getLabel(
            L"\rho_\mathrm{H_2^{BR}}",
            0,
            u"Msun * pc^-3";
            latex=true,
        )

    elseif quantity == :neutral_mass

        x_label = getLabel(
            L"\rho_\mathrm{H_I + H_2}",
            0,
            u"Msun * pc^-3";
            latex=true,
        )

    end

    # Set the y label
    y_label = getLabel(
        L"\rho_\mathrm{SFR}",
        0,
        u"Msun * yr^-1 * kpc^-3";
        latex=true,
    )

    plotSnapshot(
        [simulation_path],
        request,
        [scatter!];
        pf_kwargs=[(; color=Makie.wong_colors()[1], markersize=6, marker=:circle)],
        output_path,
        base_filename="vsf_law",
        output_format=".png",
        warnings=false,
        slice,
        filter_function,
        da_functions=[daVSFLaw],
        da_args=[(grid, quantity)],
        da_kwargs=[
            (;
                type,
                stellar_ff=dd->filterStellarAge(dd),
                gas_ff=dd->filterWithinSphere(dd, (0.0u"kpc", rmax_gas), :zero),
            ),
        ],
        post_processing=fit ? ppFitLine! : getNothing,
        pp_args=(),
        pp_kwargs=(;),
        transform_box=true,
        translation,
        rotation,
        x_trim=x_range,
        xaxis_label=L"$\log_{10}$ %$(x_label)",
        yaxis_label=L"$\log_{10}$ %$(y_label)",
        xaxis_var_name="",
        yaxis_var_name="",
        theme,
        sim_labels=[sim_label],
    )

    return nothing

end

"""
    massMetallicityRelation(
        simulation_paths::Vector{String},
        slice::IndexType;
        <keyword arguments>
    )::Nothing

Plot the resolved mass-metallicity relation. This method plots the M-Z relation at a fix moment in time.

!!! note

    Only stars younger than [`AGE_RESOLUTION`](@ref) and gas cells/particles within a sphere of radius `DISK_R` are consider.

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`.
  - `slice::IndexType`: Slice of the simulations, i.e. which snapshots will be plotted. It can be an integer (a single snapshot), a vector of integers (several snapshots), an `UnitRange` (e.g. 5:13), an `StepRange` (e.g. 5:2:13) or (:) (all snapshots). Starts at 1 and out of bounds indices are ignored.
  - `element::Symbol=:all`: Which metallicity to use. The options are:

      + `:all` -> Metallicity considering all elements, as ``Z / Z_\\odot``.
      + `:X`   -> Element ``\\mathrm{X}``, as ``12 + \\log_{10}(\\mathrm{X \\, / \\, H})``. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
  - `mass::Bool=true`: If the x axis will be the stellar mass density or the SFR density.
  - `reduce::Int=1`: Factor by which the resolution of the result will be reduced. This will be applied after the density proyection, averaging the value of neighboring pixels. It has to divide the size of `grid` exactly.
  - `output_path::String="./resolvedKSLawZScatter"`: Path to the output folder.
  - `filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all`: Which cells/particles will be plotted, the options are:

      + `:all`             -> Consider every cell/particle within the simulation box.
      + `:halo`            -> Consider only the cells/particles that belong to the main halo.
      + `:subhalo`         -> Consider only the cells/particles that belong to the main subhalo.
      + `:sphere`          -> Consider only the cell/particle inside a sphere with radius `DISK_R` (see `./src/constants/globals.jl`).
      + `:stellar_subhalo` -> Consider only the cells/particles that belong to the main subhalo.
      + `:all_subhalo`     -> Plot every cell/particle centered around the main subhalo.
      + A dictionary with three entries:

          + `:filter_function` -> The filter function.
          + `:translation`     -> Translation for the simulation box. The posibilities are:

              + `:global_cm`                  -> Selects the center of mass of the whole system as the new origin.
              + `:{component}`                -> Sets the center of mass of the given component (e.g. :stars, :gas, :halo, etc, after filtering) as the new origin. It can be any of the keys of [`PARTICLE_INDEX`](@ref).
              + `(halo_idx, subhalo_rel_idx)` -> Sets the position of the potencial minimum for the `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo) as the new origin.
              + `(halo_idx, 0)`               -> Sets the center of mass of the `halo_idx::Int` halo as the new origin.
              + `subhalo_abs_idx`             -> Sets the center of mass of the `subhalo_abs_idx::Int` as the new origin.
          + `:rotation`        -> Rotation for the simulation box. The posibilities are:

              + `:zero`                       -> No rotation is appplied.
              + `:global_am`                  -> Sets the angular momentum of the whole system as the new z axis.
              + `:stellar_am`                 -> Sets the stellar angular momentum as the new z axis.
              + `:stellar_pa`                 -> Sets the stellar principal axis as the new coordinate system.
              + `:stellar_subhalo_pa`         -> Sets the principal axis of the stars in the main subhalo as the new coordinate system.
              + `(halo_idx, subhalo_rel_idx)` -> Sets the principal axis of the stars in `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo), as the new coordinate system.
              + `(halo_idx, 0)`               -> Sets the principal axis of the stars in the `halo_idx::Int` halo, as the new coordinate system.
              + `subhalo_abs_idx`             -> Sets the principal axis of the stars in the `subhalo_abs_idx::Int` subhalo as the new coordinate system.
  - `sim_labels::Union{Vector{String},Nothing}=basename.(simulation_paths)`: Labels for the plot legend, one per simulation. Set it to `nothing` if you don't want a legend.
  - `theme::Attributes=Theme()`: Plot theme that will take precedence over [`DEFAULT_THEME`](@ref).
"""
function massMetallicityRelation(
    simulation_paths::Vector{String},
    slice::IndexType;
    element::Symbol=:all,
    mass::Bool=true,
    reduce::Int=1,
    output_path::String="./massMetallicityRelation",
    filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all,
    sim_labels::Union{Vector{String},Nothing}=basename.(simulation_paths),
    theme::Attributes=Theme(),
)::Nothing

    grid = CubicGrid(BOX_L, 400)

    (
        element ∈ [:all, keys(ELEMENT_INDEX)...] ||
        throw(ArgumentError("massMetallicityRelation: `quantity` can only be :all or any \
        of the keys of `ELEMENT_INDEX` (see `./src/constants/globals.jl`), but I got :$(quantity)"))
    )

    # Set a temporal folder for the JLD2 files
    temp_folder = joinpath(output_path, "_temp_jld2")

    # Write the JLD2 files with the density maps
    filter_function, translation, rotation, request = selectFilter(
        filter_mode,
        plotParams(:stellar_mass).request,
    )

    plotSnapshot(
        simulation_paths,
        request,
        [heatmap!];
        pf_kwargs=[(;)],
        # `plotSnapshot` configuration
        output_path=temp_folder,
        base_filename=string(:stellar_mass),
        output_format=".png",
        warnings=false,
        show_progress=true,
        # Data manipulation options
        slice,
        filter_function,
        da_functions=[daDensity2DProjection],
        da_args=[(grid, :stellar_mass, :particles)],
        da_kwargs=[(; reduce, filter_function=dd->filterStellarAge(dd))],
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
        xaxis_var_name="x",
        yaxis_var_name="y",
        xaxis_scale_func=identity,
        yaxis_scale_func=identity,
        # Plotting options
        save_figures=false,
        backup_results=true,
        theme=Theme(),
        sim_labels=nothing,
        title="",
        colorbar=false,
        # Animation options
        animation=false,
        animation_filename="animation.mp4",
        framerate=15,
    )

    if element == :all
        metal_request = plotParams(:gas_metallicity).request
    else
        metal_request = plotParams(Symbol(element, "_gas_abundance")).request
    end

    filter_function, translation, rotation, request = selectFilter(
        filter_mode,
        mergeRequests(plotParams(:gas_mass_density).request, metal_request),
    )

    # Write the JLD2 file with the metallicity density
    plotSnapshot(
        simulation_paths,
        request,
        [heatmap!];
        pf_kwargs=[(;)],
        # `plotSnapshot` configuration
        output_path=temp_folder,
        base_filename="gas_metallicity",
        output_format=".png",
        warnings=false,
        show_progress=true,
        # Data manipulation options
        slice,
        filter_function,
        da_functions=[daMetallicity2DProjection],
        da_args=[(grid, :gas, :cells)],
        da_kwargs=[(; element)],
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
        xaxis_var_name="x",
        yaxis_var_name="y",
        xaxis_scale_func=identity,
        yaxis_scale_func=identity,
        # Plotting options
        save_figures=false,
        backup_results=true,
        theme=Theme(),
        sim_labels=nothing,
        title="",
        colorbar=false,
        # Animation options
        animation=false,
        animation_filename="animation.mp4",
        framerate=15,
    )

    # Set the x label
    if mass
        x_label = getLabel(
            plotParams(:stellar_area_density).var_name,
            0,
            u"Msun * kpc^-2";
            latex=true,
        )
    else
        x_label = getLabel(
            plotParams(:sfr_area_density).var_name,
            0,
            u"Msun * yr^-1 * kpc^-2";
            latex=true,
        )
    end

    # Set the y label
    if element == :all
        ylabel = L"$\log_{10}$ %$(plotParams(:gas_metallicity).var_name)"
    else
        ylabel = plotParams(Symbol(element, "_gas_abundance")).axis_label
    end

    with_theme(theme, merge(theme_latexfonts(), DEFAULT_THEME)) do

        for (sim_idx, simulation) in pairs(simulation_paths)

            simulation_table = DataFrame(makeSimulationTable(simulation; warnings=false)[slice, :])
            sim_name         = "simulation_$(lpad(string(sim_idx), 3, "0"))"
            times            = ustrip.(u"Gyr", simulation_table[!, :physical_times])
            snapshot_numbers = simulation_table[!, :numbers]

            for (time, snapshot_number) in zip(times, snapshot_numbers)

                f = Figure()

                ax = CairoMakie.Axis(
                    f[1, 1];
                    xlabel=L"$\log_{10}$ %$(x_label)",
                    ylabel,
                    aspect=AxisAspect(1),
                )

                x_address = "stellar_mass-$(SNAP_BASENAME)_$(snapshot_number)/$(sim_name)"
                y_address = "gas_metallicity-$(SNAP_BASENAME)_$(snapshot_number)/$(sim_name)"

                jldopen(joinpath(temp_folder, "stellar_mass.jld2"), "r") do x_file

                    jldopen(joinpath(temp_folder, "gas_metallicity.jld2"), "r") do y_file

                        # Read the JLD2 files
                        x_data = vec(x_file[x_address][3])
                        y_data = vec(y_file[y_address][3])

                        # Delete 0s and NaNs in the data vectors
                        x_idxs = map(x -> isnan(x) || iszero(x), x_data)
                        y_idxs = map(x -> isnan(x) || iszero(x), y_data)

                        deleteat!(x_data, x_idxs ∪ y_idxs)
                        deleteat!(y_data, x_idxs ∪ y_idxs)

                        if mass
                            x_axis = x_data
                        else
                            x_axis = x_data .- log10(ustrip(u"yr", AGE_RESOLUTION))
                        end

                        scatter!(
                            ax,
                            x_axis,
                            y_data;
                            markersize=6,
                            color=Makie.wong_colors()[6],
                        )

                    end

                end

                ppAnnotation!(
                    f,
                    L"t = %$(rpad(round(time, sigdigits=3), 4, '0')) \, \mathrm{Gyr}",
                    position=(0.04, 0.98),
                    fontsize=25,
                )

                if !isnothing(sim_labels)
                    Makie.Legend(
                        f[1, 1],
                        [
                            MarkerElement(;
                                color=Makie.wong_colors()[6],
                                marker=:circle,
                                markersize=20,
                            ),
                        ],
                        sim_labels[sim_idx:sim_idx],
                        nbanks=1,
                        labelsize=22,
                        rowgap=-20,
                        halign=:left,
                        valign=:top,
                        padding=(13, 0, 0, 45),
                        patchlabelgap=2,
                    )
                end

                path = mkpath(joinpath(output_path, basename(simulation)))

                Makie.save(joinpath(path, "$(snapshot_number).png"), f)

            end

        end

    end

    rm(temp_folder; recursive=true)

    return nothing

end

"""
    atomicGasCubes(
        simulation_paths::Vector{String},
        slice::ReducedIndexType;
        <keyword arguments>
    )::Nothing

Create a HDF5 file with the physical position, atomic gas mass, velocity, and velocity dispersion at each voxel of a rectangular 3D grid.

The metadata for each snapshot in the HDF5 file includes the physical time in Gyr, the scale factor, and the redshift of that snapshot.

By default, the grid is centered at coordinates (0, 0, 0), has 300x300x300 voxels, and has a side length of [`BOX_L`](@ref). There are as many rows as there are voxels (27000000 by default).

The stored quantities for each voxel are:

Column 01: x coordinate [kpc]
Column 02: y coordinate [kpc]
Column 03: z coordinate [kpc]
Column 04: Atomic gas mass [Msun]
Column 05: Velocity in the x direction [km * s^-1]
Column 06: Velocity in the y direction [km * s^-1]
Column 07: Velocity in the z direction [km * s^-1]
Column 08: Velocity dispersion in the x direction [km * s^-1]
Column 09: Velocity dispersion in the y direction [km * s^-1]
Column 10: Velocity dispersion in the z direction [km * s^-1]

For Voronoi cells:

The mass is the mass of atomic gas intersecting the voxel, so it only considers the cell that it is the closest to the voxel. The velocity is given by the weighted mean of the velocities of the `n_neighbors` nearest neighbors to the voxel. And the velocity dispersion, by the weighted standard deviation.

Notice that for Voronoi cells, the mass will be sample at a sub-cell resolution (as long as voxel size < cell size), while the velocities are sample at a lower resolution (as long as `n_neighbors` > 1). The weights are given by the distance (in kpc) to each neighbor, using a Gaussian kernel.

For particles:

The mass is the accumulated mass of the particles within each voxel. The velocity is the mean of the velocities of those particles, and the velocity dispersion is the standard deviation.

If there are no particles, the mass is 0, and the velocity and velocity dispersion are NaN. If there is only one particle, the mass and velocity are the ones from that particle, and the velocity dispersion is NaN.

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`.
  - `slice::ReducedIndexType`: Slice of the simulations, i.e. which snapshots will be plotted. It can be an integer (a single snapshot), a vector of integers (several snapshots), an `UnitRange` (e.g. 5:13) or an `StepRange` (e.g. 5:2:13). Starts at 1.
  - `type::Symbol=:cells`: If the gas density will be calculated assuming the gas is in `:particles` or in Voronoi `:cells`.
  - `n_neighbors::Int=8`: Number of neighbors for the mean and standard deviation of the velocity. Setting this value to 1 maximizes the resolution for the velocity, and sets the standard deviation (columns 8, 9, and 10) to NaN.
  - `grid::CubicGrid=CubicGrid(BOX_L, 300)`: Cubic grid.
  - `output_file::String="./HI_cubes.hdf5"`: Path to the output HDF5 file.
  - `filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all`: Which cells/particles will be plotted, the options are:

      + `:all`             -> Consider every cell/particle within the simulation box.
      + `:halo`            -> Consider only the cells/particles that belong to the main halo.
      + `:subhalo`         -> Consider only the cells/particles that belong to the main subhalo.
      + `:sphere`          -> Consider only the cell/particle inside a sphere with radius `DISK_R` (see `./src/constants/globals.jl`).
      + `:stellar_subhalo` -> Consider only the cells/particles that belong to the main subhalo.
      + `:all_subhalo`     -> Plot every cell/particle centered around the main subhalo.
      + A dictionary with three entries:

          + `:filter_function` -> The filter function.
          + `:translation`     -> Translation for the simulation box. The posibilities are:

              + `:global_cm`                  -> Selects the center of mass of the whole system as the new origin.
              + `:{component}`                -> Sets the center of mass of the given component (e.g. :stars, :gas, :halo, etc, after filtering) as the new origin. It can be any of the keys of [`PARTICLE_INDEX`](@ref).
              + `(halo_idx, subhalo_rel_idx)` -> Sets the position of the potencial minimum for the `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo) as the new origin.
              + `(halo_idx, 0)`               -> Sets the center of mass of the `halo_idx::Int` halo as the new origin.
              + `subhalo_abs_idx`             -> Sets the center of mass of the `subhalo_abs_idx::Int` as the new origin.
          + `:rotation`        -> Rotation for the simulation box. The posibilities are:

              + `:zero`                       -> No rotation is appplied.
              + `:global_am`                  -> Sets the angular momentum of the whole system as the new z axis.
              + `:stellar_am`                 -> Sets the stellar angular momentum as the new z axis.
              + `:stellar_pa`                 -> Sets the stellar principal axis as the new coordinate system.
              + `:stellar_subhalo_pa`         -> Sets the principal axis of the stars in the main subhalo as the new coordinate system.
              + `(halo_idx, subhalo_rel_idx)` -> Sets the principal axis of the stars in `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo), as the new coordinate system.
              + `(halo_idx, 0)`               -> Sets the principal axis of the stars in the `halo_idx::Int` halo, as the new coordinate system.
              + `subhalo_abs_idx`             -> Sets the principal axis of the stars in the `subhalo_abs_idx::Int` subhalo as the new coordinate system.
  - `show_progress::Bool=true`: If a progress bar will be shown.
"""
function atomicGasCubes(
    simulation_paths::Vector{String},
    slice::ReducedIndexType;
    type::Symbol=:cells,
    n_neighbors::Int=8,
    grid::CubicGrid=CubicGrid(BOX_L, 300),
    output_file::String="./HI_cubes.hdf5",
    filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all,
    show_progress::Bool=true,
)::Nothing

    # Set the number of columns and rows
    n_rows = grid.n_bins^3
    n_cols = 10

    # Set the units
    m_unit = u"Msun"
    l_unit = u"kpc"
    v_unit = u"km * s^-1"
    t_unit = u"Gyr"

    filter_function, translation, rotation, request = selectFilter(
        filter_mode,
        mergeRequests(
            plotParams(:atomic_mass).request,
            Dict(:gas => ["POS ", "VEL ", "RHO ", "MASS"]),
        ),
    )

    # For gas cells, reshape the grid to conform to the way `knn` expect the matrix to be structured
    if type == :cells

        physical_grid = Matrix{Float64}(undef, 3, n_rows)

        @inbounds for i in eachindex(grid.grid)
            physical_grid[1, i] = ustrip(l_unit, grid.grid[i][1])
            physical_grid[2, i] = ustrip(l_unit, grid.grid[i][2])
            physical_grid[3, i] = ustrip(l_unit, grid.grid[i][3])
        end

    end

    # Create the output folder
    mkpath(dirname(output_file))

    # Create the output HDF5 file
    hdf5_file = h5open(output_file, "w")

    @inbounds for simulation_path in simulation_paths

        simulation_name = basename(simulation_path)

        prog_bar = Progress(
            length(slice),
            dt=0.5,
            desc="Writing the HI cube for simulation $(simulation_name)... ",
            color=:blue,
            barglyphs=BarGlyphs("|#  |"),
            enabled=show_progress,
        )

        # Create an HDF5 group for each simulation
        hdf5_group = create_group(hdf5_file, simulation_name)

        @inbounds for snap_n in slice

            data_dict = makeDataDict(
                simulation_path,
                snap_n,
                request;
                warnings=true,
            )

            snapshot_number = lpad(string(data_dict[:snap_data].global_index), 3, "0")

            # Filter the data
            filterData!(data_dict; filter_function)

            # Translate the data
            translateData!(data_dict, translation)

            # Rotate the data
            rotateData!(data_dict, rotation)

            # Load the gas quantities
            gd = data_dict[:gas]

            # Load the cell/particle positions
            positions = gd["POS "]

            # Load the cell/particle velocities
            velocities = ustrip.(v_unit, gd["VEL "])

            # Compute the mass of atomic gas in each cell
            masses = scatterQty(data_dict, :atomic_mass)

            if any(isempty, [masses, velocities, positions])
                throw(ArgumentError("atomicGasCubes: Some data is missing (there appears to be \
                no gas in the snapshot), so I cannot construct the HI cubes"))
            end

            # Column 01: x coordinate [kpc]
            # Column 02: y coordinate [kpc]
            # Column 03: z coordinate [kpc]
            # Column 04: Atomic gas mass [Msun]
            # Column 05: Velocity in the x direction [km * s^-1]
            # Column 06: Velocity in the y direction [km * s^-1]
            # Column 07: Velocity in the z direction [km * s^-1]
            # Column 08: Velocity dispersion in the x direction [km * s^-1]
            # Column 09: Velocity dispersion in the y direction [km * s^-1]
            # Column 10: Velocity dispersion in the z direction [km * s^-1]
            data_matrix = Matrix{Float64}(undef, n_rows, n_cols)

            if type == :cells

                # Compute the volume of each cell
                cell_volumes = gd["MASS"] ./ gd["RHO "]

                # Compute the atomic gas densities
                densities = ustrip.(m_unit * l_unit^-3, masses ./ cell_volumes)

                # Load the volume of the voxels
                voxel_volume = ustrip(l_unit^3, grid.bin_volume)

                # Compute the tree for a nearest neighbor search
                kdtree = KDTree(ustrip.(l_unit, positions))

                # Find the `n_neighbors` nearest cells to each voxel
                idxs, dists = knn(kdtree, physical_grid, n_neighbors, true)

                @inbounds for i in eachindex(grid.grid)

                    # Physical coordinates of the voxel [l_unit]
                    data_matrix[i, 1:3] .= ustrip.(l_unit, grid.grid[i])

                    if isone(n_neighbors)

                        # Atomic gas mass [m_unit]
                        data_matrix[i, 4] = densities[idxs[i]] * voxel_volume

                        # Neighbor velocity in the x direction [v_unit]
                        data_matrix[i, 5] = velocities[1, idxs[i]]
                        # Neighbor velocity in the y direction [v_unit]
                        data_matrix[i, 6] = velocities[2, idxs[i]]
                        # Neighbor velocity in the z direction [v_unit]
                        data_matrix[i, 7] = velocities[3, idxs[i]]

                        # For the case of only one neighbor, set the standard deviations to NaN
                        data_matrix[i, 8]  = NaN
                        data_matrix[i, 9]  = NaN
                        data_matrix[i, 10] = NaN

                    else

                        # Atomic gas mass [m_unit]
                        data_matrix[i, 4] = densities[idxs[i][1]] * voxel_volume

                        # Compute the analytic weights using a Gaussian kernel
                        neighbor_weights = aweights(evaluateNormal(dists[i]))

                        # Neighbor velocities in the x direction [v_unit]
                        vxs = velocities[1, idxs[i]]
                        # Neighbor velocities in the y direction [v_unit]
                        vys = velocities[2, idxs[i]]
                        # Neighbor velocities in the z direction [v_unit]
                        vzs = velocities[3, idxs[i]]

                        # Mean and standard deviation of the neighbor velocities in the x direction [v_unit]
                        data_matrix[i, 5], data_matrix[i, 8] = mean_and_std(vxs, neighbor_weights)
                        # Mean and standard deviation of the neighbor velocities in the y direction [v_unit]
                        data_matrix[i, 6], data_matrix[i, 9] = mean_and_std(vys, neighbor_weights)
                        # Mean and standard deviation of the neighbor velocities in the z direction [v_unit]
                        data_matrix[i, 7], data_matrix[i, 10] = mean_and_std(vzs, neighbor_weights)

                    end

                end

            elseif type == :particles

                # Find which particles are within each voxel
                idxs = listHistogram3D(positions, grid)

                @inbounds for i in eachindex(grid.grid)

                    # Physical coordinates of the voxel [l_unit]
                    data_matrix[i, 1:3] .= ustrip.(l_unit, grid.grid[i])

                    # Atomic gas mass [m_unit]
                    data_matrix[i, 4] = ustrip(m_unit, sum(masses[idxs[i]]; init=0.0 * m_unit))

                    if isempty(idxs[i])

                        # If the voxel has no particles, set the velocity and velocity dispersion to NaN
                        data_matrix[i, 5]  = NaN
                        data_matrix[i, 6]  = NaN
                        data_matrix[i, 7]  = NaN
                        data_matrix[i, 8]  = NaN
                        data_matrix[i, 9]  = NaN
                        data_matrix[i, 10] = NaN

                    elseif isone(length(idxs[i]))

                        # Velocity in the x direction [v_unit]
                        data_matrix[i, 5] = velocities[1, idxs[i][1]]
                        # Velocity in the y direction [v_unit]
                        data_matrix[i, 6] = velocities[2, idxs[i][1]]
                        # Velocity in the z direction [v_unit]
                        data_matrix[i, 7] = velocities[3, idxs[i][1]]

                        # If the voxel has a single particle, set the velocity dispersion to NaN
                        data_matrix[i, 8]  = NaN
                        data_matrix[i, 9]  = NaN
                        data_matrix[i, 10] = NaN

                    else

                        # Velocities in the x direction of the particles within the voxel [v_unit]
                        vxs = velocities[1, idxs[i]]
                        # Velocities in the y direction of the particles within the voxel [v_unit]
                        vys = velocities[2, idxs[i]]
                        # Velocities in the z direction of the particles within the voxel [v_unit]
                        vzs = velocities[3, idxs[i]]

                        # Mean and standard deviation of the velocities in the x direction [v_unit]
                        data_matrix[i, 5], data_matrix[i, 8] = mean_and_std(vxs)
                        # Mean and standard deviation of the velocities in the y direction [v_unit]
                        data_matrix[i, 6], data_matrix[i, 9] = mean_and_std(vys)
                        # Mean and standard deviation of the velocities in the z direction [v_unit]
                        data_matrix[i, 7], data_matrix[i, 10] = mean_and_std(vzs)

                    end

                end

            else

                throw(ArgumentError("atomicGasCubes: The argument `type` must be :cells or \
                :particles, but I got :$(type)"))

            end

            # Go from column-major order (Julia) to row-major order (Python and C), for interoperability
            hdf5_group["snap_$(snapshot_number)", shuffle=(), deflate=5] = permutedims(
                data_matrix,
                reverse(1:ndims(data_matrix)),
            )

            # Read the time, scale factor, and redshift
            pt = ustrip.(t_unit, data_dict[:snap_data].physical_time)
            sf = data_dict[:snap_data].scale_factor
            rs = data_dict[:snap_data].redshift

            # Write the metadata
            attrs(hdf5_group["snap_$(snapshot_number)"])["Time [Gyr]"]   = pt
            attrs(hdf5_group["snap_$(snapshot_number)"])["Scale factor"] = sf
            attrs(hdf5_group["snap_$(snapshot_number)"])["Redshift"]     = rs

            next!(prog_bar)

        end

    end

    close(hdf5_file)

    return nothing

end

"""
    clumpingFactor(
        simulation_paths::Vector{String},
        slice::IndexType,
        quantity::Symbol;
        <keyword arguments>
    )::Nothing

Plot the clumping factor of `quantity` for different volume scales.

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`.
  - `slice::IndexType`: Slice of the simulations, i.e. which snapshots will be plotted. It can be an integer (a single snapshot), a vector of integers (several snapshots), an `UnitRange` (e.g. 5:13), an `StepRange` (e.g. 5:2:13) or (:) (all snapshots). Starts at 1 and out of bounds indices are ignored.
  - `quantity::Symbol`: The number density of which quantity will be used. The options are:

      + `:gas`          -> Gas number density.
      + `:molecular`    -> Molecular hydrogen number density.
      + `:br_molecular` -> Molecular hydrogen number density, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic`       -> Atomic hydrogen number density.
      + `:ionized`      -> Ionized hydrogen number density.
      + `:neutral`      -> Neutral hydrogen number density.
  - `nn::Int=32`: Number of neighbors.
  - `smooth::Int=0`: The result will be average out using `smooth` bins for the volume. Set it to 0 if you want no smoothing.
  - `x_trim::NTuple{2,<:Real}=(-Inf, Inf)`: The data will be trim down so the x coordinates fit within `x_trim`.
  - `output_path::String="./"`: Path to the output folder.
  - `filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all`: Which cells/particles will be plotted, the options are:

      + `:all`             -> Consider every cell/particle within the simulation box.
      + `:halo`            -> Consider only the cells/particles that belong to the main halo.
      + `:subhalo`         -> Consider only the cells/particles that belong to the main subhalo.
      + `:sphere`          -> Consider only the cell/particle inside a sphere with radius `DISK_R` (see `./src/constants/globals.jl`).
      + `:stellar_subhalo` -> Consider only the cells/particles that belong to the main subhalo.
      + `:all_subhalo`     -> Plot every cell/particle centered around the main subhalo.
      + A dictionary with three entries:

          + `:filter_function` -> The filter function.
          + `:translation`     -> Translation for the simulation box. The posibilities are:

              + `:global_cm`                  -> Selects the center of mass of the whole system as the new origin.
              + `:{component}`                -> Sets the center of mass of the given component (e.g. :stars, :gas, :halo, etc, after filtering) as the new origin. It can be any of the keys of [`PARTICLE_INDEX`](@ref).
              + `(halo_idx, subhalo_rel_idx)` -> Sets the position of the potencial minimum for the `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo) as the new origin.
              + `(halo_idx, 0)`               -> Sets the center of mass of the `halo_idx::Int` halo as the new origin.
              + `subhalo_abs_idx`             -> Sets the center of mass of the `subhalo_abs_idx::Int` as the new origin.
          + `:rotation`        -> Rotation for the simulation box. The posibilities are:

              + `:zero`                       -> No rotation is appplied.
              + `:global_am`                  -> Sets the angular momentum of the whole system as the new z axis.
              + `:stellar_am`                 -> Sets the stellar angular momentum as the new z axis.
              + `:stellar_pa`                 -> Sets the stellar principal axis as the new coordinate system.
              + `:stellar_subhalo_pa`         -> Sets the principal axis of the stars in the main subhalo as the new coordinate system.
              + `(halo_idx, subhalo_rel_idx)` -> Sets the principal axis of the stars in `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo), as the new coordinate system.
              + `(halo_idx, 0)`               -> Sets the principal axis of the stars in the `halo_idx::Int` halo, as the new coordinate system.
              + `subhalo_abs_idx`             -> Sets the principal axis of the stars in the `subhalo_abs_idx::Int` subhalo as the new coordinate system.
  - `da_ff::Function=filterNothing`: A function with the signature:

    `da_ff(data_dict) -> indices`

    where

      + `data_dict::Dict`: A dictionary with the following shape:

        * `:sim_data`          -> ::Simulation (see [`Simulation`](@ref)).
        * `:snap_data`         -> ::Snapshot (see [`Snapshot`](@ref)).
        * `:gc_data`           -> ::GroupCatalog (see [`GroupCatalog`](@ref)).
        * `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * ...
        * `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * ...
      + `indices::Dict`: A dictionary with the following shape:

        * `cell/particle type` -> idxs::IndexType
        * `cell/particle type` -> idxs::IndexType
        * `cell/particle type` -> idxs::IndexType
        * ...
  - `ff_request::Dict{Symbol,Vector{String}}=Dict{Symbol,Vector{String}}()`: Request dictionary for the `da_ff` filter function.
  - `sim_labels::Union{Vector{String},Nothing}=basename.(simulation_paths)`: Labels for the plot legend, one per simulation. Set it to `nothing` if you don't want a legend.
  - `theme::Attributes=Theme()`: Plot theme that will take precedence over [`DEFAULT_THEME`](@ref).
"""
function clumpingFactor(
    simulation_paths::Vector{String},
    slice::IndexType,
    quantity::Symbol;
    nn::Int=32,
    smooth::Int=100,
    x_trim::NTuple{2,<:Real}=(-Inf, Inf),
    output_path::String="./",
    filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all,
    da_ff::Function=filterNothing,
    ff_request::Dict{Symbol,Vector{String}}=Dict{Symbol,Vector{String}}(),
    sim_labels::Union{Vector{String},Nothing}=basename.(simulation_paths),
    theme::Attributes=Theme(),
)::Nothing

    (
        quantity ∈ [:gas, :molecular, :br_molecular, :atomic, :ionized, :neutral] ||
        throw(ArgumentError("clumpingFactor: `quantity` can only be :gas, :molecular, \
        :br_molecular, :atomic, :ionized or :neutral, but I got :$(quantity)"))
    )

    filter_function, translation, rotation, request = selectFilter(
        filter_mode,
        mergeRequests(plotParams(Symbol(quantity, "_number_density")).request, ff_request),
    )

    plotSnapshot(
        simulation_paths,
        request,
        [scatter!];
        pf_kwargs=[(; markersize=10)],
        # `plotSnapshot` configuration
        output_path,
        base_filename="$(quantity)_clumping_factor",
        output_format=".png",
        warnings=false,
        show_progress=true,
        # Data manipulation options
        slice=slice,
        filter_function,
        da_functions=[daClumpingFactor],
        da_args=[(quantity,)],
        da_kwargs=[(; nn, filter_function=da_ff)],
        post_processing=getNothing,
        pp_args=(),
        pp_kwargs=(;),
        transform_box=true,
        translation,
        rotation,
        smooth,
        x_unit=u"kpc^3",
        y_unit=Unitful.NoUnits,
        x_exp_factor=0,
        y_exp_factor=0,
        x_trim,
        y_trim=(-Inf, Inf),
        x_edges=false,
        y_edges=false,
        x_func=identity,
        y_func=identity,
        # Axes options
        xaxis_label="auto_label",
        yaxis_label="auto_label",
        xaxis_var_name=L"\bar{V}",
        yaxis_var_name=L"C_\rho",
        xaxis_scale_func=identity,
        yaxis_scale_func=identity,
        # Plotting and animation options
        save_figures=true,
        backup_results=false,
        theme,
        sim_labels,
        title="",
        colorbar=false,
        # Animation options
        animation=false,
        animation_filename="animation.mp4",
        framerate=10,
    )

    return nothing

end

function clumpingFactorProfile(
    simulation_paths::Vector{String},
    slice::IndexType,
    quantity::Symbol;
    radius::Unitful.Length=DISK_R,
    n_bins::Int=100,
    output_path::String="./",
    filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all,
    da_ff::Function=filterNothing,
    ff_request::Dict{Symbol,Vector{String}}=Dict{Symbol,Vector{String}}(),
    sim_labels::Union{Vector{String},Nothing}=basename.(simulation_paths),
    theme::Attributes=Theme(),
)::Nothing

    (
        quantity ∈ [:gas, :molecular, :br_molecular, :atomic, :ionized, :neutral] ||
        throw(ArgumentError("clumpingFactorProfile: `quantity` can only be :gas, :molecular, \
        :br_molecular, :atomic, :ionized or :neutral, but I got :$(quantity)"))
    )

    filter_function, translation, rotation, request = selectFilter(
        filter_mode,
        mergeRequests(plotParams(Symbol(quantity, "_number_density")).request, ff_request),
    )

    grid = CircularGrid(radius, n_bins)

    # Draw the figures with CairoMakie
    plotSnapshot(
        simulation_paths,
        request,
        [lines!];
        pf_kwargs=[(;)],
        # `plotSnapshot` configuration
        output_path,
        base_filename="$(quantity)_clumping_factor-profile",
        output_format=".png",
        warnings=false,
        show_progress=true,
        # Data manipulation options
        slice,
        filter_function,
        da_functions=[daClumpingFactorProfile],
        da_args=[(quantity, grid)],
        da_kwargs=[(; filter_function=da_ff)],
        post_processing=getNothing,
        pp_args=(),
        pp_kwargs=(;),
        transform_box=true,
        translation,
        rotation,
        smooth=0,
        x_unit=u"kpc",
        y_unit=Unitful.NoUnits,
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
        xaxis_var_name=L"r",
        yaxis_var_name=L"C_\rho",
        xaxis_scale_func=identity,
        yaxis_scale_func=identity,
        # Plotting and animation options
        save_figures=true,
        backup_results=false,
        theme,
        sim_labels,
        title="",
        colorbar=false,
        # Animation options
        animation=false,
        animation_filename="animation.mp4",
        framerate=10,
    )

    return nothing

end
