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
  - `subhalo_rel_idx::Int=1`: Index of the target subhalo (subfind), relative the target halo. Starts at 1.
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
                    :gas => ["NHP ", "NH  ", "PRES", "FRAC", "DTIM", "TAUS", "ID  "],
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

            hydrogen_mass = sum(data_dict[:gas]["MASS"]; init=0.0u"Msun") * HYDROGEN_MASSFRAC

            hii_mass = sum(computeIonizedMass(data_dict); init=0.0u"Msun")
            hii_percent = round((hii_mass / hydrogen_mass) * 100, sigdigits=3)

            hi_mass = sum(computeAtomicMass(data_dict); init=0.0u"Msun")
            hi_percent = round((hi_mass / hydrogen_mass) * 100, sigdigits=3)

            h2_mass = sum(computeMolecularMass(data_dict); init=0.0u"Msun")
            h2_percent = round((h2_mass / hydrogen_mass) * 100, sigdigits=3)

            title = "Ionized mass:"
            title *= " "^(25 - length(title))
            println(
                file,
                "\t\t$(title)$(round(typeof(1.0u"Msun"), hii_mass, sigdigits=3)) \
                ($(hii_percent)% of total hydrogen mass)",
            )

            title = "Atomic mass:"
            title *= " "^(25 - length(title))
            println(
                file,
                "\t\t$(title)$(round(typeof(1.0u"Msun"), hi_mass, sigdigits=3)) \
                ($(hi_percent)% of total hydrogen mass)",
            )

            title = "Molecular mass:"
            title *= " "^(25 - length(title))
            println(
                file,
                "\t\t$(title)$(round(typeof(1.0u"Msun"), h2_mass, sigdigits=3)) \
                ($(h2_percent)% of total hydrogen mass)\n",
            )

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

            hydrogen_mass = sum(data_dict[:gas]["MASS"]; init=0.0u"Msun") * HYDROGEN_MASSFRAC

            hii_mass = sum(computeIonizedMass(data_dict); init=0.0u"Msun")
            hii_percent = round((hii_mass / hydrogen_mass) * 100, sigdigits=3)

            hi_mass = sum(computeAtomicMass(data_dict); init=0.0u"Msun")
            hi_percent = round((hi_mass / hydrogen_mass) * 100, sigdigits=3)

            h2_mass = sum(computeMolecularMass(data_dict); init=0.0u"Msun")
            h2_percent = round((h2_mass / hydrogen_mass) * 100, sigdigits=3)

            title = "Ionized mass:"
            title *= " "^(25 - length(title))
            println(
                file,
                "\t\t$(title)$(round(typeof(1.0u"Msun"), hii_mass, sigdigits=3)) \
                ($(hii_percent)% of total hydrogen mass)",
            )

            title = "Atomic mass:"
            title *= " "^(25 - length(title))
            println(
                file,
                "\t\t$(title)$(round(typeof(1.0u"Msun"), hi_mass, sigdigits=3)) \
                ($(hi_percent)% of total hydrogen mass)",
            )

            title = "Molecular mass:"
            title *= " "^(25 - length(title))
            println(
                file,
                "\t\t$(title)$(round(typeof(1.0u"Msun"), h2_mass, sigdigits=3)) \
                ($(h2_percent)% of total hydrogen mass)\n",
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
            title *= " "^(25 - length(title))

            println(file, "\t\t$(title)$(round.(ustrip.(u"Mpc", cm), sigdigits=6)) $(u"Mpc")")
            println(file, "\t\tDistance to global CM:   $(sqrt(sum((global_cm - cm).^2)))\n")

        end

        global_cm = round.(ustrip.(u"Mpc", global_cm), sigdigits=6)
        println(file, "\t\tGlobal center of mass:   $(global_cm) $(u"Mpc")\n")

        ############################################################################################
        # Print the fraction of gas cells that have enter our routine
        ############################################################################################

        if !isempty(data_dict[:gas]["FRAC"])

            total_number = length(data_dict[:gas]["MASS"])
            stellar_gas_number = count(!isnan, data_dict[:gas]["FRAC"][1, :])
            fraction = (stellar_gas_number / total_number) * 100

            idxs = findall(!isnan, data_dict[:gas]["FRAC"][1, :])
            stellar_gas_mass = sum(data_dict[:gas]["MASS"][idxs])
            mass_fraction = (stellar_gas_mass / sum(data_dict[:gas]["MASS"])) * 100

            println(file, "\tFraction of gas cells that have enter our routine:\n")
            println(file, "\t\t$(round(fraction, sigdigits=3))% of the cells")
            println(file, "\t\t$(round(mass_fraction, sigdigits=3))% of the mass\n")

        end

        ############################################################################################
        # Print the properties of the star forming gas
        ############################################################################################

        if !isempty(data_dict[:stars]["ACIT"])

            parz = data_dict[:stars]["PARZ"] ./ SOLAR_METALLICITY
            rhoc = ustrip.(u"cm^-3", data_dict[:stars]["RHOC"])
            acit = ustrip.(u"Myr", data_dict[:stars]["ACIT"])

            println(file, "\tProperties of the star forming gas:\n")

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
                println(file, "\nFraction of insitu stars: $(insitu_fraction)%")

            end

            ########################################################################################
            # Print the mass of each hydrogen phase between `radial_limit` y and the virial radius.
            ########################################################################################

            if :gas in component_list

                ####################################################################################
                # Disc external radius
                ####################################################################################

                radial_limit = 40.0u"kpc"

                ####################################################################################
                # Indices of cells and particles within the disc radius and
                # between the disc radius and the virial radius
                ####################################################################################

                disc_idxs = filterWithinSphere(data_dict, (0.0u"kpc", radial_limit), :zero)
                halo_idxs = filterWithinSphere(data_dict, (radial_limit, g_r_crit_200), :zero)

                stellar_masses   = data_dict[:stars]["MASS"]
                gas_masses       = data_dict[:gas]["MASS"]
                ionized_masses   = computeIonizedMass(data_dict)
                atomic_masses    = computeAtomicMass(data_dict)
                molecular_masses = computeMolecularMass(data_dict)
                neutral_masses   = computeNeutralMass(data_dict)

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

                if !isempty(neutral_masses)
                    neutral_mass_inside  = neutral_masses[disc_idxs[:gas]]
                    neutral_mass_outside = neutral_masses[halo_idxs[:gas]]
                end

                # println(file, "#"^100)
                println(file, "\nCharacteristic radii:\n")

                ############################################################################################
                # Print the radius containing 90% and 95% of the mass, withing de disc (r < `radial_limit`)
                ############################################################################################

                #############################
                # Stars
                #############################

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

                println(file, "\tRadius containing X% of the stellar mass (r < $(radial_limit)):\n")
                println(
                    file,
                    "\t\t$(round(ustrip(u"kpc", mass_radius_90), sigdigits=4)) $(u"kpc") (90%)",
                )
                println(
                    file,
                    "\t\t$(round(ustrip(u"kpc", mass_radius_95), sigdigits=4)) $(u"kpc") (95%)\n",
                )

                #############################
                # Total gas
                #############################

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
                    "\tRadius containing X% of the total gas mass (r < $(radial_limit)):\n",
                )
                println(
                    file,
                    "\t\t$(round(ustrip(u"kpc", mass_radius_90), sigdigits=4)) $(u"kpc") (90%)",
                )
                println(
                    file,
                    "\t\t$(round(ustrip(u"kpc", mass_radius_95), sigdigits=4)) $(u"kpc") (95%)\n",
                )

                #############################
                # Ionized gas
                #############################

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
                        "\tRadius containing X% of the ionized gas mass (r < $(radial_limit)):\n",
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

                #############################
                # Atomic gas
                #############################

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
                        "\tRadius containing X% of the atomic gas mass (r < $(radial_limit)):\n",
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

                #############################
                # Molecular gas
                #############################

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
                        "\tRadius containing X% of the molecular gas mass (r < $(radial_limit)):\n",
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

                #############################
                # Neutral gas
                #############################

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
                        "\tRadius containing X% of the neutral gas mass (r < $(radial_limit)):\n",
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
                # Print the masses withing de disc (r < `radial_limit`) and
                # outside the disc (`radial_limit` < r < R200)
                ####################################################################################

                println(file, "#"^100)
                println(file, "\nCharacteristic fractions and masses:\n")

                println(file, "\t", "#"^20)
                println(file, "\tR200: $(round(typeof(1.0u"kpc"), g_r_crit_200, sigdigits=4))")
                println(file, "\t", "#"^20, "\n")

                #############################
                # Stars
                #############################

                total_stellar_mass_inside  = sum(stellar_mass_inside; init=0.0u"Msun")
                total_stellar_mass_outside = sum(stellar_mass_outside; init=0.0u"Msun")
                total_stellar_mass         = total_stellar_mass_inside + total_stellar_mass_outside

                s_inside_percent  = (total_stellar_mass_inside / total_stellar_mass) * 100.0
                s_outside_percent = (total_stellar_mass_outside / total_stellar_mass) * 100.0

                println(file, "\t", "#"^40)
                println(file, "\tStars:")
                println(file, "\t", "#"^40, "\n")

                println(file, "\tStellar mass inside the disc (r < $(radial_limit)):\n")
                println(
                    file,
                    "\t\t$(round(typeof(1.0u"Msun"), total_stellar_mass_inside, sigdigits=3)) \
                    ($(round(s_inside_percent, sigdigits=3))% of the total stellar mass)\n",
                )

                println(file, "\tStellar mass outside the disc ($(radial_limit) < r < R200):\n")
                println(
                    file,
                    "\t\t$(round(typeof(1.0u"Msun"), total_stellar_mass_outside, sigdigits=3)) \
                    ($(round(s_outside_percent, sigdigits=3))% of the total stellar mass)\n",
                )

                #############################
                # Total gas
                #############################

                total_gas_mass_inside  = sum(gas_mass_inside; init=0.0u"Msun")
                total_gas_mass_outside = sum(gas_mass_outside; init=0.0u"Msun")
                total_gas_mass         = total_gas_mass_inside + total_gas_mass_outside

                g_inside_percent  = (total_gas_mass_inside / total_gas_mass) * 100.0
                g_outside_percent = (total_gas_mass_outside / total_gas_mass) * 100.0

                println(file, "\t", "#"^40)
                println(file, "\tTotal gas:")
                println(file, "\t", "#"^40, "\n")

                println(file, "\tGas mass inside the disc (r < $(radial_limit)):\n")
                println(
                    file,
                    "\t\t$(round(typeof(1.0u"Msun"), total_gas_mass_inside, sigdigits=3)) \
                    ($(round(g_inside_percent, sigdigits=3))% of the total gas mass)\n",
                )

                println(file, "\tGas mass outside the disc ($(radial_limit) < r < R200):\n")
                println(
                    file,
                    "\t\t$(round(typeof(1.0u"Msun"), total_gas_mass_outside, sigdigits=3)) \
                    ($(round(g_outside_percent, sigdigits=3))% of the total gas mass)\n",
                )

                total_hydrogen_mass = total_gas_mass * HYDROGEN_MASSFRAC

                #############################
                # Ionized gas
                #############################

                if !isempty(ionized_masses)

                    total_ion_mass_inside  = sum(ionized_mass_inside; init=0.0u"Msun")
                    total_ion_mass_outside = sum(ionized_mass_outside; init=0.0u"Msun")

                    i_inside_percent  = (total_ion_mass_inside  / total_hydrogen_mass) * 100.0
                    i_outside_percent = (total_ion_mass_outside / total_hydrogen_mass) * 100.0

                    println(file, "\t", "#"^40)
                    println(file, "\tIonized gas:")
                    println(file, "\t", "#"^40, "\n")

                    println(file, "\tIonized mass inside the disc (r < $(radial_limit)):\n")
                    println(
                        file,
                        "\t\t$(round(typeof(1.0u"Msun"), total_ion_mass_inside, sigdigits=3)) \
                        ($(round(i_inside_percent, sigdigits=3))% of the total hydrogen mass)\n",
                    )

                    println(file, "\tIonized mass outside the disc ($(radial_limit) < r < R200):\n")
                    println(
                        file,
                        "\t\t$(round(typeof(1.0u"Msun"), total_ion_mass_outside, sigdigits=3)) \
                        ($(round(i_outside_percent, sigdigits=3))% of the total hydrogen mass)\n",
                    )

                end

                #############################
                # Atomic gas
                #############################

                if !isempty(atomic_masses)

                    total_ato_mass_inside  = sum(atomic_mass_inside; init=0.0u"Msun")
                    total_ato_mass_outside = sum(atomic_mass_outside; init=0.0u"Msun")

                    a_inside_percent  = (total_ato_mass_inside  / total_hydrogen_mass) * 100.0
                    a_outside_percent = (total_ato_mass_outside / total_hydrogen_mass) * 100.0

                    println(file, "\t", "#"^40)
                    println(file, "\tAtomic gas:")
                    println(file, "\t", "#"^40, "\n")

                    println(file, "\tAtomic mass inside the disc (r < $(radial_limit)):\n")
                    println(
                        file,
                        "\t\t$(round(typeof(1.0u"Msun"), total_ato_mass_inside, sigdigits=3)) \
                        ($(round(a_inside_percent, sigdigits=3))% of the total hydrogen mass)\n",
                    )

                    println(file, "\tAtomic mass outside the disc ($(radial_limit) < r < R200):\n")
                    println(
                        file,
                        "\t\t$(round(typeof(1.0u"Msun"), total_ato_mass_outside, sigdigits=3)) \
                        ($(round(a_outside_percent, sigdigits=3))% of the total hydrogen mass)\n",
                    )

                end

                #############################
                # Molecular gas
                #############################

                if !isempty(molecular_masses)

                    total_mol_mass_inside  = sum(molecular_mass_inside; init=0.0u"Msun")
                    total_mol_mass_outside = sum(molecular_mass_outside; init=0.0u"Msun")

                    m_inside_percent  = (total_mol_mass_inside  / total_hydrogen_mass) * 100.0
                    m_outside_percent = (total_mol_mass_outside / total_hydrogen_mass) * 100.0

                    println(file, "\t", "#"^40)
                    println(file, "\tMolecular gas:")
                    println(file, "\t", "#"^40, "\n")

                    println(file, "\tMolecular mass inside the disc (r < $(radial_limit)):\n")
                    println(
                        file,
                        "\t\t$(round(typeof(1.0u"Msun"), total_mol_mass_inside, sigdigits=3)) \
                        ($(round(m_inside_percent, sigdigits=3))% of the total hydrogen mass)\n",
                    )

                    println(file, "\tMolecular mass outside the disc ($(radial_limit) < r < R200):\n")
                    println(
                        file,
                        "\t\t$(round(typeof(1.0u"Msun"), total_mol_mass_outside, sigdigits=3)) \
                        ($(round(m_outside_percent, sigdigits=3))% of the total hydrogen mass)\n",
                    )

                end

                #############################
                # Neutral gas
                #############################

                if !isempty(neutral_masses)

                    total_neu_mass_inside  = sum(neutral_mass_inside; init=0.0u"Msun")
                    total_neu_mass_outside = sum(neutral_mass_outside; init=0.0u"Msun")

                    m_inside_percent  = (total_neu_mass_inside  / total_hydrogen_mass) * 100.0
                    m_outside_percent = (total_neu_mass_outside / total_hydrogen_mass) * 100.0

                    println(file, "\t", "#"^40)
                    println(file, "\tNeutral gas:")
                    println(file, "\t", "#"^40, "\n")

                    println(file, "\tNeutral mass inside the disc (r < $(radial_limit)):\n")
                    println(
                        file,
                        "\t\t$(round(typeof(1.0u"Msun"), total_neu_mass_inside, sigdigits=3)) \
                        ($(round(m_inside_percent, sigdigits=3))% of the total hydrogen mass)\n",
                    )

                    println(file, "\tNeutral mass outside the disc ($(radial_limit) < r < R200):\n")
                    println(
                        file,
                        "\t\t$(round(typeof(1.0u"Msun"), total_neu_mass_outside, sigdigits=3)) \
                        ($(round(m_outside_percent, sigdigits=3))% of the total hydrogen mass)\n",
                    )

                end

            end

            ########################################################################################
            # Halo and subhalo global properties
            ########################################################################################

            println(file, "#"^100)
            println(file, "\nHalo and subhalo global properties:")

            ##############################
            # Print the halo properties
            ##############################

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

            ##############################
            # Print the subhalo properties
            ##############################

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

      + `:stellar_mass`   -> Stellar mass.
      + `:gas_mass`       -> Gas mass.
      + `:hydrogen_mass`  -> Hydrogen mass.
      + `:dm_mass`        -> Dark matter mass.
      + `:bh_mass`        -> Black hole mass.
      + `:molecular_mass` -> Molecular hydrogen (``\\mathrm{H_2}``) mass.
      + `:atomic_mass`    -> Atomic hydrogen (``\\mathrm{HI}``) mass.
      + `:ionized_mass`   -> Ionized hydrogen (``\\mathrm{HII}``) mass.
      + `:neutral_mass`   -> Neutral hydrogen (``\\mathrm{HI + H_2}``) mass.
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

      + `:stellar_mass`   -> Stellar mass.
      + `:gas_mass`       -> Gas mass.
      + `:hydrogen_mass`  -> Hydrogen mass.
      + `:dm_mass`        -> Dark matter mass.
      + `:bh_mass`        -> Black hole mass.
      + `:molecular_mass` -> Molecular hydrogen (``\\mathrm{H_2}``) mass.
      + `:atomic_mass`    -> Atomic hydrogen (``\\mathrm{HI}``) mass.
      + `:ionized_mass`   -> Ionized hydrogen (``\\mathrm{HII}``) mass.
      + `:neutral_mass`   -> Neutral hydrogen (``\\mathrm{HI + H_2}``) mass.
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
  - `type::Symbol=:cells`: Component type for the temperature fields. It can be either `:particles` or Voronoi `:cells`.
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

      + `:stellar_mass`               -> Stellar mass.
      + `:gas_mass`                   -> Gas mass.
      + `:hydrogen_mass`              -> Hydrogen mass.
      + `:dm_mass`                    -> Dark matter mass.
      + `:bh_mass`                    -> Black hole mass.
      + `:molecular_mass`             -> Molecular hydrogen (``\\mathrm{H_2}``) mass.
      + `:atomic_mass`                -> Atomic hydrogen (``\\mathrm{HI}``) mass.
      + `:ionized_mass`               -> Ionized hydrogen (``\\mathrm{HII}``) mass.
      + `:neutral_mass`               -> Neutral hydrogen (``\\mathrm{HI + H_2}``) mass.
      + `:molecular_fraction`         -> Gas mass fraction of molecular hydrogen.
      + `:atomic_fraction`            -> Gas mass fraction of atomic hydrogen.
      + `:ionized_fraction`           -> Gas mass fraction of ionized hydrogen.
      + `:neutral_fraction`           -> Gas mass fraction of neutral hydrogen.
      + `:molecular_neutral_fraction` -> Fraction of molecular hydrogen in the neutral gas.
      + `:gas_mass_density`           -> Gas mass density.
      + `:hydrogen_mass_density`      -> Hydrogen mass density.
      + `:gas_number_density`         -> Gas number density.
      + `:molecular_number_density`   -> Molecular hydrogen number density.
      + `:atomic_number_density`      -> Atomic hydrogen number density.
      + `:ionized_number_density`     -> Ionized hydrogen number density.
      + `:neutral_number_density`     -> Neutral hydrogen number density.
      + `:gas_metallicity`            -> Mass fraction of all elements above He in the gas (solar units).
      + `:stellar_metallicity`        -> Mass fraction of all elements above He in the stars (solar units).
      + `:X_gas_abundance`            -> Gas abundance of element ``\\mathrm{X}``, as ``12 + \\log_{10}(\\mathrm{X \\, / \\, H})``. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
      + `:X_stellar_abundance`        -> Stellar abundance of element ``\\mathrm{X}``, as ``12 + \\log_{10}(\\mathrm{X \\, / \\, H})``. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
      + `:stellar_radial_distance`    -> Distance of every stellar particle to the origin.
      + `:gas_radial_distance`        -> Distance of every gas cell to the origin.
      + `:dm_radial_distance`         -> Distance of every dark matter particle to the origin.
      + `:stellar_xy_distance`        -> Projected distance of every stellar particle to the origin.
      + `:gas_xy_distance`            -> Projected distance of every gas cell to the origin.
      + `:dm_xy_distance`             -> Projected distance of every dark matter particle to the origin.
      + `:gas_sfr`                    -> SFR associated to each gas particle/cell within the code.
      + `:stellar_circularity`        -> Stellar circularity.
      + `:stellar_vcirc`              -> Stellar circular velocity.
      + `:stellar_vradial`            -> Stellar radial speed.
      + `:stellar_vtangential`        -> Stellar tangential speed.
      + `:stellar_vzstar`             -> Stellar speed in the z direction, computed as ``v_z \\, \\mathrm{sign}(z)``.
      + `:stellar_age`                -> Stellar age.
      + `:sfr`                        -> The star formation rate.
      + `:ssfr`                       -> The specific star formation rate.
      + `:observational_sfr`          -> The star formation rate of the last `AGE_RESOLUTION`.
      + `:observational_ssfr`         -> The specific star formation rate of the last `AGE_RESOLUTION`.
      + `:temperature`                -> Gas temperature, as ``\\log_{10}(T \\, / \\, \\mathrm{K})``.
      + `:pressure`                   -> Gas pressure.
  - `y_quantity::Symbol`: Quantity for the y axis. The options are:

      + `:stellar_mass`               -> Stellar mass.
      + `:gas_mass`                   -> Gas mass.
      + `:hydrogen_mass`              -> Hydrogen mass.
      + `:dm_mass`                    -> Dark matter mass.
      + `:bh_mass`                    -> Black hole mass.
      + `:molecular_mass`             -> Molecular hydrogen (``\\mathrm{H_2}``) mass.
      + `:atomic_mass`                -> Atomic hydrogen (``\\mathrm{HI}``) mass.
      + `:ionized_mass`               -> Ionized hydrogen (``\\mathrm{HII}``) mass.
      + `:neutral_mass`               -> Neutral hydrogen (``\\mathrm{HI + H_2}``) mass.
      + `:molecular_fraction`         -> Gas mass fraction of molecular hydrogen.
      + `:atomic_fraction`            -> Gas mass fraction of atomic hydrogen.
      + `:ionized_fraction`           -> Gas mass fraction of ionized hydrogen.
      + `:neutral_fraction`           -> Gas mass fraction of neutral hydrogen.
      + `:molecular_neutral_fraction` -> Fraction of molecular hydrogen in the neutral gas.
      + `:gas_mass_density`           -> Gas mass density.
      + `:hydrogen_mass_density`      -> Hydrogen mass density.
      + `:gas_number_density`         -> Gas number density.
      + `:molecular_number_density`   -> Molecular hydrogen number density.
      + `:atomic_number_density`      -> Atomic hydrogen number density.
      + `:ionized_number_density`     -> Ionized hydrogen number density.
      + `:neutral_number_density`     -> Neutral hydrogen number density.
      + `:gas_metallicity`            -> Mass fraction of all elements above He in the gas (solar units).
      + `:stellar_metallicity`        -> Mass fraction of all elements above He in the stars (solar units).
      + `:X_gas_abundance`            -> Gas abundance of element ``\\mathrm{X}``, as ``12 + \\log_{10}(\\mathrm{X \\, / \\, H})``. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
      + `:X_stellar_abundance`        -> Stellar abundance of element ``\\mathrm{X}``, as ``12 + \\log_{10}(\\mathrm{X \\, / \\, H})``. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
      + `:stellar_radial_distance`    -> Distance of every stellar particle to the origin.
      + `:gas_radial_distance`        -> Distance of every gas cell to the origin.
      + `:dm_radial_distance`         -> Distance of every dark matter particle to the origin.
      + `:stellar_xy_distance`        -> Projected distance of every stellar particle to the origin.
      + `:gas_xy_distance`            -> Projected distance of every gas cell to the origin.
      + `:dm_xy_distance`             -> Projected distance of every dark matter particle to the origin.
      + `:gas_sfr`                    -> SFR associated to each gas particle/cell within the code.
      + `:stellar_circularity`        -> Stellar circularity.
      + `:stellar_vcirc`              -> Stellar circular velocity.
      + `:stellar_vradial`            -> Stellar radial speed.
      + `:stellar_vtangential`        -> Stellar tangential speed.
      + `:stellar_vzstar`             -> Stellar speed in the z direction, computed as ``v_z \\, \\mathrm{sign}(z)``.
      + `:stellar_age`                -> Stellar age.
      + `:sfr`                        -> The star formation rate.
      + `:ssfr`                       -> The specific star formation rate.
      + `:observational_sfr`          -> The star formation rate of the last `AGE_RESOLUTION`.
      + `:observational_ssfr`         -> The specific star formation rate of the last `AGE_RESOLUTION`.
      + `:temperature`                -> Gas temperature, as ``\\log_{10}(T \\, / \\, \\mathrm{K})``.
      + `:pressure`                   -> Gas pressure.
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
function scatterPlot(
    simulation_paths::Vector{String},
    slice::IndexType,
    x_quantity::Symbol,
    y_quantity::Symbol;
    output_path::String="./",
    filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all,
    theme::Attributes=Theme(),
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

      + `:stellar_mass`               -> Stellar mass.
      + `:gas_mass`                   -> Gas mass.
      + `:hydrogen_mass`              -> Hydrogen mass.
      + `:dm_mass`                    -> Dark matter mass.
      + `:bh_mass`                    -> Black hole mass.
      + `:molecular_mass`             -> Molecular hydrogen (``\\mathrm{H_2}``) mass.
      + `:atomic_mass`                -> Atomic hydrogen (``\\mathrm{HI}``) mass.
      + `:ionized_mass`               -> Ionized hydrogen (``\\mathrm{HII}``) mass.
      + `:neutral_mass`               -> Neutral hydrogen (``\\mathrm{HI + H_2}``) mass.
      + `:molecular_fraction`         -> Gas mass fraction of molecular hydrogen.
      + `:atomic_fraction`            -> Gas mass fraction of atomic hydrogen.
      + `:ionized_fraction`           -> Gas mass fraction of ionized hydrogen.
      + `:neutral_fraction`           -> Gas mass fraction of neutral hydrogen.
      + `:molecular_neutral_fraction` -> Fraction of molecular hydrogen in the neutral gas.
      + `:gas_mass_density`           -> Gas mass density.
      + `:hydrogen_mass_density`      -> Hydrogen mass density.
      + `:gas_number_density`         -> Gas number density.
      + `:molecular_number_density`   -> Molecular hydrogen number density.
      + `:atomic_number_density`      -> Atomic hydrogen number density.
      + `:ionized_number_density`     -> Ionized hydrogen number density.
      + `:neutral_number_density`     -> Neutral hydrogen number density.
      + `:gas_metallicity`            -> Mass fraction of all elements above He in the gas (solar units).
      + `:stellar_metallicity`        -> Mass fraction of all elements above He in the stars (solar units).
      + `:X_gas_abundance`            -> Gas abundance of element ``\\mathrm{X}``, as ``12 + \\log_{10}(\\mathrm{X \\, / \\, H})``. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
      + `:X_stellar_abundance`        -> Stellar abundance of element ``\\mathrm{X}``, as ``12 + \\log_{10}(\\mathrm{X \\, / \\, H})``. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
      + `:stellar_radial_distance`    -> Distance of every stellar particle to the origin.
      + `:gas_radial_distance`        -> Distance of every gas cell to the origin.
      + `:dm_radial_distance`         -> Distance of every dark matter particle to the origin.
      + `:stellar_xy_distance`        -> Projected distance of every stellar particle to the origin.
      + `:gas_xy_distance`            -> Projected distance of every gas cell to the origin.
      + `:dm_xy_distance`             -> Projected distance of every dark matter particle to the origin.
      + `:gas_sfr`                    -> SFR associated to each gas particle/cell within the code.
      + `:stellar_circularity`        -> Stellar circularity.
      + `:stellar_vcirc`              -> Stellar circular velocity.
      + `:stellar_vradial`            -> Stellar radial speed.
      + `:stellar_vtangential`        -> Stellar tangential speed.
      + `:stellar_vzstar`             -> Stellar speed in the z direction, computed as ``v_z \\, \\mathrm{sign}(z)``.
      + `:stellar_age`                -> Stellar age.
      + `:sfr`                        -> The star formation rate.
      + `:ssfr`                       -> The specific star formation rate.
      + `:observational_sfr`          -> The star formation rate of the last `AGE_RESOLUTION`.
      + `:observational_ssfr`         -> The specific star formation rate of the last `AGE_RESOLUTION`.
      + `:temperature`                -> Gas temperature, as ``\\log_{10}(T \\, / \\, \\mathrm{K})``.
      + `:pressure`                   -> Gas pressure.
  - `y_quantity::Symbol`: Quantity for the y axis. The options are:

      + `:stellar_mass`               -> Stellar mass.
      + `:gas_mass`                   -> Gas mass.
      + `:hydrogen_mass`              -> Hydrogen mass.
      + `:dm_mass`                    -> Dark matter mass.
      + `:bh_mass`                    -> Black hole mass.
      + `:molecular_mass`             -> Molecular hydrogen (``\\mathrm{H_2}``) mass.
      + `:atomic_mass`                -> Atomic hydrogen (``\\mathrm{HI}``) mass.
      + `:ionized_mass`               -> Ionized hydrogen (``\\mathrm{HII}``) mass.
      + `:neutral_mass`               -> Neutral hydrogen (``\\mathrm{HI + H_2}``) mass.
      + `:molecular_fraction`         -> Gas mass fraction of molecular hydrogen.
      + `:atomic_fraction`            -> Gas mass fraction of atomic hydrogen.
      + `:ionized_fraction`           -> Gas mass fraction of ionized hydrogen.
      + `:neutral_fraction`           -> Gas mass fraction of neutral hydrogen.
      + `:molecular_neutral_fraction` -> Fraction of molecular hydrogen in the neutral gas.
      + `:gas_mass_density`           -> Gas mass density.
      + `:hydrogen_mass_density`      -> Hydrogen mass density.
      + `:gas_number_density`         -> Gas number density.
      + `:molecular_number_density`   -> Molecular hydrogen number density.
      + `:atomic_number_density`      -> Atomic hydrogen number density.
      + `:ionized_number_density`     -> Ionized hydrogen number density.
      + `:neutral_number_density`     -> Neutral hydrogen number density.
      + `:gas_metallicity`            -> Mass fraction of all elements above He in the gas (solar units).
      + `:stellar_metallicity`        -> Mass fraction of all elements above He in the stars (solar units).
      + `:X_gas_abundance`            -> Gas abundance of element ``\\mathrm{X}``, as ``12 + \\log_{10}(\\mathrm{X \\, / \\, H})``. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
      + `:X_stellar_abundance`        -> Stellar abundance of element ``\\mathrm{X}``, as ``12 + \\log_{10}(\\mathrm{X \\, / \\, H})``. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
      + `:stellar_radial_distance`    -> Distance of every stellar particle to the origin.
      + `:gas_radial_distance`        -> Distance of every gas cell to the origin.
      + `:dm_radial_distance`         -> Distance of every dark matter particle to the origin.
      + `:stellar_xy_distance`        -> Projected distance of every stellar particle to the origin.
      + `:gas_xy_distance`            -> Projected distance of every gas cell to the origin.
      + `:dm_xy_distance`             -> Projected distance of every dark matter particle to the origin.
      + `:gas_sfr`                    -> SFR associated to each gas particle/cell within the code.
      + `:stellar_circularity`        -> Stellar circularity.
      + `:stellar_vcirc`              -> Stellar circular velocity.
      + `:stellar_vradial`            -> Stellar radial speed.
      + `:stellar_vtangential`        -> Stellar tangential speed.
      + `:stellar_vzstar`             -> Stellar speed in the z direction, computed as ``v_z \\, \\mathrm{sign}(z)``.
      + `:stellar_age`                -> Stellar age.
      + `:sfr`                        -> The star formation rate.
      + `:ssfr`                       -> The specific star formation rate.
      + `:observational_sfr`          -> The star formation rate of the last `AGE_RESOLUTION`.
      + `:observational_ssfr`         -> The specific star formation rate of the last `AGE_RESOLUTION`.
      + `:temperature`                -> Gas temperature, as ``\\log_{10}(T \\, / \\, \\mathrm{K})``.
      + `:pressure`                   -> Gas pressure.
  - `z_quantity::Symbol`: Quantity for the z axis (weights). The options are:

      + `:stellar_mass`               -> Stellar mass.
      + `:gas_mass`                   -> Gas mass.
      + `:hydrogen_mass`              -> Hydrogen mass.
      + `:dm_mass`                    -> Dark matter mass.
      + `:bh_mass`                    -> Black hole mass.
      + `:molecular_mass`             -> Molecular hydrogen (``\\mathrm{H_2}``) mass.
      + `:atomic_mass`                -> Atomic hydrogen (``\\mathrm{HI}``) mass.
      + `:ionized_mass`               -> Ionized hydrogen (``\\mathrm{HII}``) mass.
      + `:neutral_mass`               -> Neutral hydrogen (``\\mathrm{HI + H_2}``) mass.
      + `:molecular_fraction`         -> Gas mass fraction of molecular hydrogen.
      + `:atomic_fraction`            -> Gas mass fraction of atomic hydrogen.
      + `:ionized_fraction`           -> Gas mass fraction of ionized hydrogen.
      + `:neutral_fraction`           -> Gas mass fraction of neutral hydrogen.
      + `:molecular_neutral_fraction` -> Fraction of molecular hydrogen in the neutral gas.
      + `:gas_mass_density`           -> Gas mass density.
      + `:hydrogen_mass_density`      -> Hydrogen mass density.
      + `:gas_number_density`         -> Gas number density.
      + `:molecular_number_density`   -> Molecular hydrogen number density.
      + `:atomic_number_density`      -> Atomic hydrogen number density.
      + `:ionized_number_density`     -> Ionized hydrogen number density.
      + `:neutral_number_density`     -> Neutral hydrogen number density.
      + `:gas_metallicity`            -> Mass fraction of all elements above He in the gas (solar units).
      + `:stellar_metallicity`        -> Mass fraction of all elements above He in the stars (solar units).
      + `:X_gas_abundance`            -> Gas abundance of element ``\\mathrm{X}``, as ``12 + \\log_{10}(\\mathrm{X \\, / \\, H})``. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
      + `:X_stellar_abundance`        -> Stellar abundance of element ``\\mathrm{X}``, as ``12 + \\log_{10}(\\mathrm{X \\, / \\, H})``. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
      + `:stellar_radial_distance`    -> Distance of every stellar particle to the origin.
      + `:gas_radial_distance`        -> Distance of every gas cell to the origin.
      + `:dm_radial_distance`         -> Distance of every dark matter particle to the origin.
      + `:stellar_xy_distance`        -> Projected distance of every stellar particle to the origin.
      + `:gas_xy_distance`            -> Projected distance of every gas cell to the origin.
      + `:dm_xy_distance`             -> Projected distance of every dark matter particle to the origin.
      + `:gas_sfr`                    -> SFR associated to each gas particle/cell within the code.
      + `:stellar_circularity`        -> Stellar circularity.
      + `:stellar_vcirc`              -> Stellar circular velocity.
      + `:stellar_vradial`            -> Stellar radial speed.
      + `:stellar_vtangential`        -> Stellar tangential speed.
      + `:stellar_vzstar`             -> Stellar speed in the z direction, computed as ``v_z \\, \\mathrm{sign}(z)``.
      + `:stellar_age`                -> Stellar age.
      + `:sfr`                        -> The star formation rate.
      + `:ssfr`                       -> The specific star formation rate.
      + `:observational_sfr`          -> The star formation rate of the last `AGE_RESOLUTION`.
      + `:observational_ssfr`         -> The specific star formation rate of the last `AGE_RESOLUTION`.
      + `:temperature`                -> Gas temperature, as ``\\log_{10}(T \\, / \\, \\mathrm{K})``.
      + `:pressure`                   -> Gas pressure.
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
    theme::Attributes=Theme(),
)::Nothing

    x_plot_params = plotParams(x_quantity)
    y_plot_params = plotParams(y_quantity)

    filter_function, translation, rotation, request = selectFilter(
        filter_mode,
        mergeRequests(
            x_plot_params.request,
            y_plot_params.request,
            plotParams(:gas_mass_density).request,
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
            da_kwargs=[(; x_range, y_range, x_log, y_log, total, n_bins, filter_function=da_ff)],
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
  - `subhalo_rel_idx::Int`: Index of the target subhalo (subfind), relative the target halo. Starts at 1. If set to 0, all subhalos of the target halo are included.
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

      + `:gas_mass`                   -> Gas mass.
      + `:hydrogen_mass`              -> Hydrogen mass.
      + `:molecular_mass`             -> Molecular hydrogen (``\\mathrm{H_2}``) mass.
      + `:atomic_mass`                -> Atomic hydrogen (``\\mathrm{HI}``) mass.
      + `:ionized_mass`               -> Ionized hydrogen (``\\mathrm{HII}``) mass.
      + `:neutral_mass`               -> Neutral hydrogen (``\\mathrm{HI + H_2}``) mass.
      + `:molecular_fraction`         -> Gas mass fraction of molecular hydrogen.
      + `:atomic_fraction`            -> Gas mass fraction of atomic hydrogen.
      + `:ionized_fraction`           -> Gas mass fraction of ionized hydrogen.
      + `:neutral_fraction`           -> Gas mass fraction of neutral hydrogen.
      + `:molecular_neutral_fraction` -> Fraction of molecular hydrogen in the neutral gas.
      + `:gas_mass_density`           -> Gas mass density.
      + `:hydrogen_mass_density`      -> Hydrogen mass density.
      + `:gas_number_density`         -> Gas number density.
      + `:molecular_number_density`   -> Molecular hydrogen number density.
      + `:atomic_number_density`      -> Atomic hydrogen number density.
      + `:ionized_number_density`     -> Ionized hydrogen number density.
      + `:neutral_number_density`     -> Neutral hydrogen number density.
      + `:gas_metallicity`            -> Mass fraction of all elements above He in the gas (solar units).
      + `:X_gas_abundance`            -> Gas abundance of element ``\\mathrm{X}``, as ``12 + \\log_{10}(\\mathrm{X \\, / \\, H})``. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
      + `:gas_radial_distance`        -> Distance of every gas cell to the origin.
      + `:gas_xy_distance`            -> Projected distance of every gas cell to the origin.
      + `:gas_sfr`                    -> SFR associated to each gas particle/cell within the code.
      + `:temperature`                -> Gas temperature, as ``\\log_{10}(T \\, / \\, \\mathrm{K})``.
      + `:pressure`                   -> Gas pressure.
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

      + `:stellar_mass`           -> Stellar mass.
      + `:gas_mass`               -> Gas mass.
      + `:hydrogen_mass`          -> Hydrogen mass.
      + `:dm_mass`                -> Dark matter mass.
      + `:bh_mass`                -> Black hole mass.
      + `:molecular_mass`         -> Molecular hydrogen (``\\mathrm{H_2}``) mass.
      + `:atomic_mass`            -> Atomic hydrogen (``\\mathrm{HI}``) mass.
      + `:ionized_mass`           -> Ionized hydrogen (``\\mathrm{HII}``) mass.
      + `:neutral_mass`           -> Neutral hydrogen (``\\mathrm{HI + H_2}``) mass.
      + `:stellar_number`         -> Number of stellar particles.
      + `:gas_number`             -> Number of gas cells.
      + `:dm_number`              -> Number of dark matter particles.
      + `:bh_number`              -> Number of black hole particles.
      + `:molecular_fraction`     -> Gas mass fraction of molecular hydrogen.
      + `:atomic_fraction`        -> Gas mass fraction of atomic hydrogen.
      + `:ionized_fraction`       -> Gas mass fraction of ionized hydrogen.
      + `:neutral_fraction`       -> Gas mass fraction of neutral hydrogen.
      + `:stellar_area_density`   -> Stellar area mass density, for a radius of `DISK_R`.
      + `:gas_area_density`       -> Gas area mass density, for a radius of `DISK_R`.
      + `:molecular_area_density` -> Molecular hydrogen area mass density, for a radius of `DISK_R`.
      + `:atomic_area_density`    -> Atomic hydrogen area mass density, for a radius of `DISK_R`.
      + `:ionized_area_density`   -> Ionized hydrogen area mass density, for a radius of `DISK_R`.
      + `:neutral_area_density`   -> Neutral hydrogen area mass density, for a radius of `DISK_R`.
      + `:sfr_area_density`       -> Star formation rate area density, for the last `AGE_RESOLUTION` and a radius of `DISK_R`.
      + `:gas_metallicity`        -> Mass fraction of all elements above He in the gas (solar units).
      + `:stellar_metallicity`    -> Mass fraction of all elements above He in the stars (solar units).
      + `:X_gas_abundance`        -> Gas abundance of element ``\\mathrm{X}``, as ``12 + \\log_{10}(\\mathrm{X \\, / \\, H})``. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
      + `:X_stellar_abundance`    -> Stellar abundance of element ``\\mathrm{X}``, as ``12 + \\log_{10}(\\mathrm{X \\, / \\, H})``. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
      + `:stellar_specific_am`    -> Norm of the stellar specific angular momentum.
      + `:gas_specific_am`        -> Norm of the gas specific angular momentum.
      + `:dm_specific_am`         -> Norm of the dark matter specific angular momentum.
      + `:sfr`                    -> The star formation rate.
      + `:ssfr`                   -> The specific star formation rate.
      + `:observational_sfr`      -> The star formation rate of the last `AGE_RESOLUTION`.
      + `:observational_ssfr`     -> The specific star formation rate of the last `AGE_RESOLUTION`.
      + `:scale_factor`           -> Scale factor.
      + `:redshift`               -> Redshift.
      + `:physical_time`          -> Physical time since the Big Bang.
      + `:lookback_time`          -> Physical time left to reach the last snapshot.
  - `y_quantity::Symbol`: Quantity for the y axis. The options are:

      + `:stellar_mass`           -> Stellar mass.
      + `:gas_mass`               -> Gas mass.
      + `:hydrogen_mass`          -> Hydrogen mass.
      + `:dm_mass`                -> Dark matter mass.
      + `:bh_mass`                -> Black hole mass.
      + `:molecular_mass`         -> Molecular hydrogen (``\\mathrm{H_2}``) mass.
      + `:atomic_mass`            -> Atomic hydrogen (``\\mathrm{HI}``) mass.
      + `:ionized_mass`           -> Ionized hydrogen (``\\mathrm{HII}``) mass.
      + `:neutral_mass`           -> Neutral hydrogen (``\\mathrm{HI + H_2}``) mass.
      + `:stellar_number`         -> Number of stellar particles.
      + `:gas_number`             -> Number of gas cells.
      + `:dm_number`              -> Number of dark matter particles.
      + `:bh_number`              -> Number of black hole particles.
      + `:molecular_fraction`     -> Gas mass fraction of molecular hydrogen.
      + `:atomic_fraction`        -> Gas mass fraction of atomic hydrogen.
      + `:ionized_fraction`       -> Gas mass fraction of ionized hydrogen.
      + `:neutral_fraction`       -> Gas mass fraction of neutral hydrogen.
      + `:stellar_area_density`   -> Stellar area mass density, for a radius of `DISK_R`.
      + `:gas_area_density`       -> Gas area mass density, for a radius of `DISK_R`.
      + `:molecular_area_density` -> Molecular hydrogen area mass density, for a radius of `DISK_R`.
      + `:atomic_area_density`    -> Atomic hydrogen area mass density, for a radius of `DISK_R`.
      + `:ionized_area_density`   -> Ionized hydrogen area mass density, for a radius of `DISK_R`.
      + `:neutral_area_density`   -> Neutral hydrogen area mass density, for a radius of `DISK_R`.
      + `:sfr_area_density`       -> Star formation rate area density, for the last `AGE_RESOLUTION` and a radius of `DISK_R`.
      + `:gas_metallicity`        -> Mass fraction of all elements above He in the gas (solar units).
      + `:stellar_metallicity`    -> Mass fraction of all elements above He in the stars (solar units).
      + `:X_gas_abundance`        -> Gas abundance of element ``\\mathrm{X}``, as ``12 + \\log_{10}(\\mathrm{X \\, / \\, H})``. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
      + `:X_stellar_abundance`    -> Stellar abundance of element ``\\mathrm{X}``, as ``12 + \\log_{10}(\\mathrm{X \\, / \\, H})``. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
      + `:stellar_specific_am`    -> Norm of the stellar specific angular momentum.
      + `:gas_specific_am`        -> Norm of the gas specific angular momentum.
      + `:dm_specific_am`         -> Norm of the dark matter specific angular momentum.
      + `:sfr`                    -> The star formation rate.
      + `:ssfr`                   -> The specific star formation rate.
      + `:observational_sfr`      -> The star formation rate of the last `AGE_RESOLUTION`.
      + `:observational_ssfr`     -> The specific star formation rate of the last `AGE_RESOLUTION`.
      + `:scale_factor`           -> Scale factor.
      + `:redshift`               -> Redshift.
      + `:physical_time`          -> Physical time since the Big Bang.
      + `:lookback_time`          -> Physical time left to reach the last snapshot.
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
        save_figures=true,
        backup_results=false,
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

      + `:stellar_mass`               -> Stellar mass.
      + `:gas_mass`                   -> Gas mass.
      + `:hydrogen_mass`              -> Hydrogen mass.
      + `:dm_mass`                    -> Dark matter mass.
      + `:bh_mass`                    -> Black hole mass.
      + `:molecular_mass`             -> Molecular hydrogen (``\\mathrm{H_2}``) mass.
      + `:atomic_mass`                -> Atomic hydrogen (``\\mathrm{HI}``) mass.
      + `:ionized_mass`               -> Ionized hydrogen (``\\mathrm{HII}``) mass.
      + `:neutral_mass`               -> Neutral hydrogen (``\\mathrm{HI + H_2}``) mass.
      + `:sfr`                        -> The star formation rate.
      + `:ssfr`                       -> The specific star formation rate.
      + `:observational_sfr`          -> The star formation rate of the last `AGE_RESOLUTION`.
      + `:observational_ssfr`         -> The specific star formation rate of the last `AGE_RESOLUTION`.
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

    elseif quantity == :atomic_mass

        yaxis_var_name = L"\Sigma_\mathrm{HI}"

    elseif quantity == :ionized_mass

        yaxis_var_name = L"\Sigma_\mathrm{HII}"

    elseif quantity == :neutral_mass

        yaxis_var_name = L"\Sigma_\mathrm{H2 + HI}"

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

      + `:stellar_mass`               -> Stellar mass.
      + `:gas_mass`                   -> Gas mass.
      + `:hydrogen_mass`              -> Hydrogen mass.
      + `:dm_mass`                    -> Dark matter mass.
      + `:bh_mass`                    -> Black hole mass.
      + `:molecular_mass`             -> Molecular hydrogen (``\\mathrm{H_2}``) mass.
      + `:atomic_mass`                -> Atomic hydrogen (``\\mathrm{HI}``) mass.
      + `:ionized_mass`               -> Ionized hydrogen (``\\mathrm{HII}``) mass.
      + `:neutral_mass`               -> Neutral hydrogen (``\\mathrm{HI + H_2}``) mass.
      + `:sfr`                        -> The star formation rate.
      + `:ssfr`                       -> The specific star formation rate.
      + `:observational_sfr`          -> The star formation rate of the last `AGE_RESOLUTION`.
      + `:observational_ssfr`         -> The specific star formation rate of the last `AGE_RESOLUTION`.
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

      + `:stellar_mass`               -> Stellar mass.
      + `:gas_mass`                   -> Gas mass.
      + `:hydrogen_mass`              -> Hydrogen mass.
      + `:dm_mass`                    -> Dark matter mass.
      + `:bh_mass`                    -> Black hole mass.
      + `:molecular_mass`             -> Molecular hydrogen (``\\mathrm{H_2}``) mass.
      + `:atomic_mass`                -> Atomic hydrogen (``\\mathrm{HI}``) mass.
      + `:ionized_mass`               -> Ionized hydrogen (``\\mathrm{HII}``) mass.
      + `:neutral_mass`               -> Neutral hydrogen (``\\mathrm{HI + H_2}``) mass.
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
  - `slice::IndexType`: Slice of the simulations, i.e. which snapshots will be plotted. It can be an integer (a single snapshot), a vector of integers (several snapshots), an `UnitRange` (e.g. 5:13), an `StepRange` (e.g. 5:2:13) or (:) (all snapshots). Starts at 1 and out of bounds indices are ignored.
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
        da_args=[(:stellar_circularity, grid,)],
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

      + `:stellar_mass`   -> Stellar mass.
      + `:molecular_mass` -> Molecular hydrogen (``\\mathrm{H_2}``) mass.
      + `:atomic_mass`    -> Atomic hydrogen (``\\mathrm{HI}``) mass.
      + `:sfr`            -> The star formation rate of the last `AGE_RESOLUTION`.
  - `y_quantity::Symbol`: Quantity for the y axis. The options are:

      + `:stellar_mass`   -> Stellar mass.
      + `:molecular_mass` -> Molecular hydrogen (``\\mathrm{H_2}``) mass.
      + `:atomic_mass`    -> Atomic hydrogen (``\\mathrm{HI}``) mass.
      + `:sfr`            -> The star formation rate of the last `AGE_RESOLUTION`.
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
        theme=merge(
            theme,
            Theme(size=(850, 850), Legend=(nbanks=1,)),
        ),
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

      + `:stellar_area_density`   -> Stellar area mass density.
      + `:molecular_area_density` -> Molecular hydrogen area mass density.
      + `:atomic_area_density`    -> Atomic hydrogen area mass density.
      + `:sfr_area_density`       -> Star formation rate area density, for the last `AGE_RESOLUTION`.
      + `:O_stellar_abundance`    -> Stellar abundance of oxygen, as ``12 + \\log_{10}(\\mathrm{O \\, / \\, H})``.
      + `:N_stellar_abundance`    -> Stellar abundance of nitrogen, as ``12 + \\log_{10}(\\mathrm{N \\, / \\, H})``.
      + `:C_stellar_abundance`    -> Stellar abundance of carbon, as ``12 + \\log_{10}(\\mathrm{C \\, / \\, H})``.
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

    elseif quantity ∈ [:molecular_area_density, :sfr_area_density]

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
    resolvedKennicuttSchmidtLaw(
        simulation_paths::Vector{String},
        slice::IndexType;
        <keyword arguments>
    )::Nothing

Plot the resolved Kennicutt-Schmidt relation plus the results of Kennicutt (1998) or Bigiel et al. (2008), depending on the chosen `quantity`. This method plots the KS relation at a fix moment in time.

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`.
  - `slice::IndexType`: Slice of the simulations, i.e. which snapshots will be plotted. It can be an integer (a single snapshot), a vector of integers (several snapshots), an `UnitRange` (e.g. 5:13), an `StepRange` (e.g. 5:2:13) or (:) (all snapshots). Starts at 1 and out of bounds indices are ignored.
  - `quantity::Symbol=:molecular_mass`: Quantity for the x axis. The options are:

      + `:gas_mass`       -> Gas area mass density. This one will be plotted with the results of Kennicutt (1998).
      + `:molecular_mass` -> Molecular hydrogen area mass density. This one will be plotted with the results of Bigiel et al. (2008).
      + `:neutral_mass`   -> Neutral hydrogen area mass density. This one will be plotted with the results of Bigiel et al. (2008).
  - `type::Symbol=:cells`: If the density in the x axis will be calculated assuming gas as `:particles` or Voronoi `:cells`.
  - `output_path::String="./resolvedKennicuttSchmidtLaw"`: Path to the output folder.
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

R. C. Kennicutt (1998). *The Global Schmidt Law in Star-forming Galaxies*. The Astrophysical Journal, **498(2)**, 541-552. [doi:10.1086/305588](https://doi.org/10.1086/305588)

F. Bigiel et al. (2008). *THE STAR FORMATION LAW IN NEARBY GALAXIES ON SUB-KPC SCALES*. The Astrophysical Journal, **136(6)**, 2846. [doi:10.1088/0004-6256/136/6/2846](https://doi.org/10.1088/0004-6256/136/6/2846)
"""
function resolvedKennicuttSchmidtLaw(
    simulation_paths::Vector{String},
    slice::IndexType;
    quantity::Symbol=:molecular_mass,
    type::Symbol=:cells,
    output_path::String="./resolvedKennicuttSchmidtLaw",
    filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all,
    sim_labels::Union{Vector{String},Nothing}=basename.(simulation_paths),
    theme::Attributes=Theme(),
)::Nothing

    grid = CubicGrid(65.0u"kpc", 300)

    (
        quantity ∈ [:gas_mass, :molecular_mass, :neutral_mass] ||
        throw(ArgumentError("resolvedKennicuttSchmidtLaw: `quantity` can only be :gas_mass, \
        :molecular_mass or :neutral_mass, but I got :$(quantity)"))
    )

    # Set a temporal folder for the JLD2 files
    temp_folder = joinpath(output_path, "_temp_jld2")

    # Write the JLD2 files with the density maps
    for (qty, type) in zip([:stellar_mass, quantity], [:particles, type])

        filter_function, translation, rotation, request = selectFilter(
            filter_mode,
            plotParams(qty).request,
        )

        plotSnapshot(
            simulation_paths,
            request,
            [heatmap!];
            pf_kwargs=[(;)],
            # `plotSnapshot` configuration
            output_path=temp_folder,
            base_filename=string(qty),
            output_format=".png",
            warnings=false,
            show_progress=true,
            # Data manipulation options
            slice,
            filter_function,
            da_functions=[daDensity2DProjection],
            da_args=[(grid, qty, type)],
            da_kwargs=[(; filter_function=dd->filterStellarAge(dd))],
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

    end

    filter_function, translation, rotation, request = selectFilter(
        filter_mode,
        mergeRequests(
            plotParams(:gas_metallicity).request,
            plotParams(:gas_mass_density).request,
        ),
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
        da_args=[(grid, :gas, type)],
        da_kwargs=[(;)],
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

    with_theme(merge(theme, theme_latexfonts(), DEFAULT_THEME)) do

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
                    ylabel=L"$\log_{10}$ %$(y_label)",
                    aspect=AxisAspect(1),
                )

                x_address = "$(quantity)-$(SNAP_BASENAME)_$(snapshot_number)/$(sim_name)"
                y_address = "stellar_mass-$(SNAP_BASENAME)_$(snapshot_number)/$(sim_name)"
                z_address = "gas_metallicity-$(SNAP_BASENAME)_$(snapshot_number)/$(sim_name)"

                jldopen(joinpath(temp_folder, "$(string(quantity)).jld2"), "r") do x_file

                    jldopen(joinpath(temp_folder, "stellar_mass.jld2"), "r") do y_file

                        jldopen(joinpath(temp_folder, "gas_metallicity.jld2"), "r") do z_file

                            # Read the JLD2 files
                            x_data = vec(x_file[x_address][3])
                            y_data = vec(y_file[y_address][3])
                            z_data = vec(z_file[z_address][3])

                            # Delete 0s and NaNs in the data vectors
                            x_idxs = map(x -> isnan(x) || iszero(x), x_data)
                            y_idxs = map(x -> isnan(x) || iszero(x), y_data)
                            z_idxs = map(x -> isnan(x) || iszero(x), z_data)

                            deleteat!(x_data, x_idxs ∪ y_idxs ∪ z_idxs)
                            deleteat!(y_data, x_idxs ∪ y_idxs ∪ z_idxs)
                            deleteat!(z_data, x_idxs ∪ y_idxs ∪ z_idxs)

                            pf = scatter!(
                                ax,
                                x_data,
                                y_data .- log10(ustrip(u"yr", AGE_RESOLUTION));
                                markersize=6,
                                color=z_data,
                                colormap=:nipy_spectral,
                            )

                            # Print the colorbar
                            Colorbar(
                                f[1, 2],
                                pf,
                                label=L"$\log_{10}$ %$(plotParams(:gas_metallicity).var_name)",
                                labelsize=25,
                                ticklabelsize=25,
                                labelpadding=4,
                            )

                        end

                    end

                end

                ppAnnotation!(
                    f,
                    L"t = %$(rpad(round(time, sigdigits=3), 4, '0')) \, \mathrm{Gyr}",
                    position=(0.08, 0.98),
                    fontsize=28,
                )

                if quantity == :gas_mass

                    pp_legend = ppKennicutt1998!(
                        f;
                        x_unit=u"Msun * kpc^-2",
                        y_unit=u"Msun * yr^-1 * kpc^-2",
                        x_log=true,
                        y_log=true,
                        color=Makie.wong_colors()[1],
                        linestyle=nothing,
                        linewidth=3,
                        warnings=false,
                    )

                elseif quantity == :molecular_mass

                    pp_legend = ppBigiel2008!(
                        f,
                        true;
                        x_unit=u"Msun * kpc^-2",
                        y_unit=u"Msun * yr^-1 * kpc^-2",
                        x_log=true,
                        y_log=true,
                        color=Makie.wong_colors()[1],
                        linestyle=nothing,
                        linewidth=3,
                        warnings=false,
                    )

                elseif quantity == :neutral_mass

                    pp_legend = ppBigiel2008!(
                        f,
                        false;
                        x_unit=u"Msun * kpc^-2",
                        y_unit=u"Msun * yr^-1 * kpc^-2",
                        x_log=true,
                        y_log=true,
                        color=Makie.wong_colors()[1],
                        linestyle=nothing,
                        linewidth=3,
                        warnings=false,
                    )

                end

                if !isnothing(sim_labels)
                    Makie.Legend(
                        f[1, 1],
                        [
                            MarkerElement(; color=:gray40, marker=:circle, markersize=20),
                            pp_legend[1],
                        ],
                        [sim_labels[sim_idx], pp_legend[2]],
                        nbanks=1,
                        labelsize=25,
                        rowgap=-20,
                        halign=:left,
                        valign=:top,
                        padding=(15, 0, 0, 30),
                    )
                end

                rowsize!(f.layout, 1, Makie.Fixed(pixelarea(ax.scene)[].widths[2]))

                path = mkpath(
                    joinpath(output_path, basename(simulation), string(quantity)),
                )

                Makie.save(joinpath(path, "$(snapshot_number).png"), f)

            end

        end

    end

    rm(temp_folder; recursive=true)

end

#TODO
function resolvedKennicuttSchmidtLaw2(
    simulation_paths::Vector{String},
    slice::IndexType;
    quantity::Symbol=:molecular_mass,
    type::Symbol=:cells,
    x_range::Union{NTuple{2,<:Number},Nothing}=nothing,
    y_range::Union{NTuple{2,<:Number},Nothing}=nothing,
    n_bins::Int=100,
    output_path::String="./resolvedKennicuttSchmidtLaw",
    filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all,
    print_range::Bool=false,
    theme::Attributes=Theme(),
)::Nothing

    grid = CubicGrid(65.0u"kpc", 300)

    (
        quantity ∈ [:gas_mass, :molecular_mass, :neutral_mass] ||
        throw(ArgumentError("resolvedKennicuttSchmidtLaw: `quantity` can only be :gas_mass, \
        :molecular_mass or :neutral_mass, but I got :$(quantity)"))
    )

    # Set a temporal folder for the JLD2 files
    temp_folder = joinpath(output_path, "_temp_jld2")

    # Write the JLD2 files with the density maps
    for (qty, type) in zip([:stellar_mass, quantity], [:particles, type])

        filter_function, translation, rotation, request = selectFilter(
            filter_mode,
            plotParams(qty).request,
        )

        plotSnapshot(
            simulation_paths,
            request,
            [heatmap!];
            pf_kwargs=[(;)],
            # `plotSnapshot` configuration
            output_path=temp_folder,
            base_filename=string(qty),
            output_format=".png",
            warnings=false,
            show_progress=true,
            # Data manipulation options
            slice,
            filter_function,
            da_functions=[daDensity2DProjection],
            da_args=[(grid, qty, type)],
            da_kwargs=[(; filter_function=dd->filterStellarAge(dd))],
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

    end

    filter_function, translation, rotation, request = selectFilter(
        filter_mode,
        mergeRequests(
            plotParams(:gas_metallicity).request,
            plotParams(:gas_mass_density).request,
        ),
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
        da_args=[(grid, :gas, type)],
        da_kwargs=[(;)],
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

    with_theme(merge(theme, theme_latexfonts(), DEFAULT_THEME)) do

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
                    ylabel=L"$\log_{10}$ %$(y_label)",
                    aspect=AxisAspect(1),
                )

                x_address = "$(quantity)-$(SNAP_BASENAME)_$(snapshot_number)/$(sim_name)"
                y_address = "stellar_mass-$(SNAP_BASENAME)_$(snapshot_number)/$(sim_name)"
                z_address = "gas_metallicity-$(SNAP_BASENAME)_$(snapshot_number)/$(sim_name)"

                jldopen(joinpath(temp_folder, "$(string(quantity)).jld2"), "r") do x_file

                    jldopen(joinpath(temp_folder, "stellar_mass.jld2"), "r") do y_file

                        jldopen(joinpath(temp_folder, "gas_metallicity.jld2"), "r") do z_file

                            # Read the JLD2 files
                            x_data = vec(x_file[x_address][3])
                            y_data = vec(y_file[y_address][3])
                            z_data = vec(z_file[z_address][3])

                            # Delete 0s and NaNs in the data vectors
                            x_idxs = map(x -> isnan(x) || iszero(x), x_data)
                            y_idxs = map(x -> isnan(x) || iszero(x), y_data)
                            z_idxs = map(x -> isnan(x) || iszero(x), z_data)

                            deleteat!(x_data, x_idxs ∪ y_idxs ∪ z_idxs)
                            deleteat!(y_data, x_idxs ∪ y_idxs ∪ z_idxs)
                            deleteat!(z_data, x_idxs ∪ y_idxs ∪ z_idxs)

                            y_data = y_data .- log10(ustrip(u"yr", AGE_RESOLUTION))

                            # If there is no range specified, use the extrema of the x values
                            if isnothing(x_range)
                                xrange = extrema(x_data)
                            else
                                xrange = x_range
                            end

                            # If there is no range specified, use the extrema of the y values
                            if isnothing(y_range)
                                yrange = extrema(y_data)
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

                            # Compute the 2D histogram
                            values = Float64.(histogram2D(
                                permutedims(hcat(x_data, y_data), (2, 1)),
                                collect(range(xrange[1], xrange[2]; length=n_bins + 1)),
                                collect(range(yrange[1], yrange[2]; length=n_bins + 1));
                            ))

                            # Set bins with a value of 0 to NaN
                            replace!(x -> iszero(x) ? NaN : x, values)

                            # The transpose and reverse operation are to conform to the way heatmap!
                            # expect the matrix to be structured, and log10 is used to enhance the contrast
                            z_axis = reverse!(transpose(log10.(values)), dims=2)

                            if print_range

                                # Compute the mininimum and maximum of `z_axis`
                                min_max = isempty(z_axis) ? (NaN, NaN) : extrema(filter(!isnan, z_axis))

                                # Print the counts range
                                @info(
                                    "\nCounts range \
                                    \n  Simulation:    $(basename(simulation)) \
                                    \n  Snapshot:      $(snapshot_number) \
                                    \n  Quantity:      $(quantity) \
                                    \n  log₁₀(counts): $(min_max)\n\n"
                                )

                            end

                            pf = heatmap!(
                                ax,
                                x_axis,
                                y_axis,
                                z_axis;
                                colormap=:nipy_spectral,
                            )

                        end

                    end

                end

                ppAnnotation!(
                    f,
                    L"t = %$(rpad(round(time, sigdigits=3), 4, '0')) \, \mathrm{Gyr}",
                    position=(0.04, 0.98),
                    fontsize=28,
                    color=:white,
                )

                rowsize!(f.layout, 1, Makie.Fixed(pixelarea(ax.scene)[].widths[2]))

                path = mkpath(
                    joinpath(output_path, basename(simulation), string(quantity)),
                )

                Makie.save(joinpath(path, "$(snapshot_number).png"), f)

            end

        end

    end

    rm(temp_folder; recursive=true)

end

"""
    integratedKennicuttSchmidtLaw(
        simulation_paths::Vector{String},
        slice::IndexType;
        <keyword arguments>
    )::Nothing

Plot the integarted Kennicutt-Schmidt relation plus the results of Kennicutt (1998) or Bigiel et al. (2008), depending on the chosen `quantity`. This method plots the KS relation for the whole galaxy at different points in time.

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`.
  - `slice::IndexType`: Slice of the simulations, i.e. which snapshots will be plotted. It can be an integer (a single snapshot), a vector of integers (several snapshots), an `UnitRange` (e.g. 5:13), an `StepRange` (e.g. 5:2:13) or (:) (all snapshots). Starts at 1 and out of bounds indices are ignored.
  - `quantity::Symbol=:molecular_mass`: Quantity for the x axis. The options are:

      + `:gas_mass`       -> Gas area mass density. This one will be plotted with the results of Kennicutt (1998).
      + `:molecular_mass` -> Molecular hydrogen area mass density. This one will be plotted with the results of Bigiel et al. (2008).
      + `:neutral_mass`   -> Neutral hydrogen area mass density. This one will be plotted with the results of Bigiel et al. (2008).
  - `type::Symbol=:cells`: If the density in the x axis will be calculated assuming gas as `:particles` or Voronoi `:cells`.
  - `output_path::String="./resolvedKennicuttSchmidtLaw"`: Path to the output folder.
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

R. C. Kennicutt (1998). *The Global Schmidt Law in Star-forming Galaxies*. The Astrophysical Journal, **498(2)**, 541-552. [doi:10.1086/305588](https://doi.org/10.1086/305588)

F. Bigiel et al. (2008). *THE STAR FORMATION LAW IN NEARBY GALAXIES ON SUB-KPC SCALES*. The Astrophysical Journal, **136(6)**, 2846. [doi:10.1088/0004-6256/136/6/2846](https://doi.org/10.1088/0004-6256/136/6/2846)
"""
function integratedKennicuttSchmidtLaw(
    simulation_paths::Vector{String},
    slice::IndexType;
    quantity::Symbol=:molecular_mass,
    type::Symbol=:cells,
    output_path::String="./integratedKennicuttSchmidtLaw",
    filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all,
    sim_labels::Union{Vector{String},Nothing}=basename.(simulation_paths),
    theme::Attributes=Theme(),
)::Nothing

    grid = CubicGrid(65.0u"kpc", 300)

    (
        quantity ∈ [:gas_mass, :molecular_mass, :neutral_mass] ||
        throw(ArgumentError("integratedKennicuttSchmidtLaw: `quantity` can only be :gas_mass, \
        :molecular_mass or :neutral_mass, but I got :$(quantity)"))
    )

    # Set a temporal folder for the JLD2 files
    temp_folder = joinpath(output_path, "_temp_jld2")

    # Write the JLD2 files with the density maps
    for (qty, type) in zip([:stellar_mass, quantity], [:particles, type])

        filter_function, translation, rotation, request = selectFilter(
            filter_mode,
            plotParams(qty).request,
        )

        plotSnapshot(
            simulation_paths,
            request,
            [heatmap!];
            pf_kwargs=[(;)],
            # `plotSnapshot` configuration
            output_path=temp_folder,
            base_filename=string(qty),
            output_format=".png",
            warnings=false,
            show_progress=true,
            # Data manipulation options
            slice,
            filter_function,
            da_functions=[daDensity2DProjection],
            da_args=[(grid, qty, type)],
            da_kwargs=[(; filter_function=dd->filterStellarAge(dd))],
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

    end

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

    log10Δt = log10(ustrip(u"yr", AGE_RESOLUTION))

    with_theme(merge(theme, theme_latexfonts(), DEFAULT_THEME)) do

        for (sim_idx, simulation) in pairs(simulation_paths)

            simulation_table = DataFrame(makeSimulationTable(simulation; warnings=false)[slice, :])
            sim_name         = "simulation_$(lpad(string(sim_idx), 3, "0"))"
            snapshot_numbers = simulation_table[!, :numbers]

            f = Figure()

            ax = CairoMakie.Axis(
                f[1, 1];
                xlabel=L"$\log_{10}$ %$(x_label)",
                ylabel=L"$\log_{10}$ %$(y_label)",
                aspect=AxisAspect(1),
            )

            # Allocate memory
            x_values = similar(snapshot_numbers, Measurement{Float64})
            y_values = similar(snapshot_numbers, Measurement{Float64})

            for (snap_idx, snapshot_number) in pairs(snapshot_numbers)

                x_address = "$(quantity)-$(SNAP_BASENAME)_$(snapshot_number)/$(sim_name)"
                y_address = "stellar_mass-$(SNAP_BASENAME)_$(snapshot_number)/$(sim_name)"

                jldopen(joinpath(temp_folder, "$(string(quantity)).jld2"), "r") do x_file

                    jldopen(joinpath(temp_folder, "stellar_mass.jld2"), "r") do y_file

                        # Read the JLD2 files
                        x_data = vec(x_file[x_address][3])
                        y_data = vec(y_file[y_address][3])

                        # Delete 0s and NaNs in the data vectors
                        x_idxs = map(x -> isnan(x) || iszero(x), x_data)
                        y_idxs = map(x -> isnan(x) || iszero(x), y_data)

                        deleteat!(x_data, x_idxs ∪ y_idxs)
                        deleteat!(y_data, x_idxs ∪ y_idxs)

                        lin_x = exp10.(x_data)
                        lin_y = exp10.(y_data .- log10Δt)

                        x_values[snap_idx] = log10(median(lin_x) ± mad(lin_x; normalize=false))
                        y_values[snap_idx] = log10(median(lin_y) ± mad(lin_y; normalize=false))

                    end

                end

            end

            scatter!(
                ax,
                Measurements.value.(x_values),
                Measurements.value.(y_values);
                color=Makie.wong_colors()[2],
                markersize=15,
            )

            errorbars!(
                ax,
                Measurements.value.(x_values),
                Measurements.value.(y_values),
                Measurements.uncertainty.(y_values),
                direction=:y,
                color=Makie.wong_colors()[2],
            )

            errorbars!(
                ax,
                Measurements.value.(x_values),
                Measurements.value.(y_values),
                Measurements.uncertainty.(x_values),
                direction=:x,
                color=Makie.wong_colors()[2],
            )

            if quantity == :gas_mass

                pp_legend = ppKennicutt1998!(
                    f;
                    x_unit=u"Msun * kpc^-2",
                    y_unit=u"Msun * yr^-1 * kpc^-2",
                    x_log=true,
                    y_log=true,
                    color=Makie.wong_colors()[1],
                    linestyle=nothing,
                    linewidth=3,
                    warnings=false,
                )

            elseif quantity == :molecular_mass

                pp_legend = ppBigiel2008!(
                    f,
                    true;
                    x_unit=u"Msun * kpc^-2",
                    y_unit=u"Msun * yr^-1 * kpc^-2",
                    x_log=true,
                    y_log=true,
                    color=Makie.wong_colors()[1],
                    linestyle=nothing,
                    linewidth=3,
                    warnings=false,
                )

            elseif quantity == :neutral_mass

                pp_legend = ppBigiel2008!(
                    f,
                    false;
                    x_unit=u"Msun * kpc^-2",
                    y_unit=u"Msun * yr^-1 * kpc^-2",
                    x_log=true,
                    y_log=true,
                    color=Makie.wong_colors()[1],
                    linestyle=nothing,
                    linewidth=3,
                    warnings=false,
                )

            end

            if !isnothing(sim_labels)
                Makie.Legend(
                    f[1, 1],
                    [
                        MarkerElement(;
                            color=Makie.wong_colors()[2],
                            marker=:circle,
                            markersize=15,
                        ),
                        pp_legend[1],
                    ],
                    [sim_labels[sim_idx], pp_legend[2]],
                    nbanks=1,
                    labelsize=25,
                    rowgap=-10,
                    padding=(0, 10, 10, 0),
                )
            end

            path = mkpath(joinpath(output_path, string(quantity)))

            Makie.save(joinpath(path, "$(basename(simulation)).png"), f)

        end

    end

    rm(temp_folder; recursive=true)

end

"""
    fitResolvedKennicuttSchmidtLaw(
        simulation_path::String,
        slice::IndexType;
        <keyword arguments>
    )::Nothing

Plot the resolved Kennicutt-Schmidt relation with its linear fit.

# Arguments

  - `simulation_path::String`: Path to the simulation directory, set in the code variable `OutputDir`.
  - `slice::IndexType`: Slice of the simulations, i.e. which snapshots will be plotted. It can be an integer (a single snapshot), a vector of integers (several snapshots), an `UnitRange` (e.g. 5:13), an `StepRange` (e.g. 5:2:13) or (:) (all snapshots). Starts at 1 and out of bounds indices are ignored.
  - `quantity::Symbol=:molecular_mass`: Quantity for the x axis. The options are:

      + `:gas_mass`       -> Gas area mass density. This one will be plotted with the results of Kennicutt (1998).
      + `:molecular_mass` -> Molecular hydrogen area mass density. This one will be plotted with the results of Bigiel et al. (2008).
      + `:neutral_mass`   -> Neutral hydrogen area mass density. This one will be plotted with the results of Bigiel et al. (2008).
  - `type::Symbol=:cells`: If the density in the x axis will be calculated assuming gas as `:particles` or Voronoi `:cells`.
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
function fitResolvedKennicuttSchmidtLaw(
    simulation_path::String,
    slice::IndexType;
    quantity::Symbol=:molecular_mass,
    type::Symbol=:cells,
    x_range::NTuple{2,<:Real}=(-Inf, Inf),
    output_path::String="./",
    filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all,
    sim_label::Union{String,Nothing}=basename(simulation_path),
    theme::Attributes=Theme(),
)::Nothing

    grid = CubicGrid(65.0u"kpc", 300)

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
        base_filename="ks-law_fit",
        output_format=".png",
        warnings=false,
        show_progress=true,
        # Data manipulation options
        slice,
        filter_function,
        da_functions=[daKennicuttSchmidtLaw],
        da_args=[(grid, quantity)],
        da_kwargs=[(; type, filter_function=dd->filterStellarAge(dd))],
        post_processing=ppFitLine!,
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
        ##########################
        # left, right, bottom, top
        ##########################
        theme=merge(
            theme,
            Theme(
                size=(880, 880),
                Axis=(aspect=AxisAspect(1),),
                Legend=(
                    labelsize=20,
                    halign=:left,
                    valign=:top,
                    padding=(15, 0, 0, 125),
                ),
            ),
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
    resolvedMassMetallicityRelation(
        simulation_paths::Vector{String},
        slice::IndexType;
        <keyword arguments>
    )::Nothing

Plot the resolved mass-metallicity relation. This method plots the M-Z relation at a fix moment in time.

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`.
  - `slice::IndexType`: Slice of the simulations, i.e. which snapshots will be plotted. It can be an integer (a single snapshot), a vector of integers (several snapshots), an `UnitRange` (e.g. 5:13), an `StepRange` (e.g. 5:2:13) or (:) (all snapshots). Starts at 1 and out of bounds indices are ignored.
  - `element::Symbol=:all`: Which metallicity to use. The options are:

      + `:all` -> Metallicity considering all elements, as ``Z / Z_\\odot``.
      + `:X`   -> Xlement ``\\mathrm{X}``, as ``12 + \\log_{10}(\\mathrm{X \\, / \\, H})``. The possibilities are the keys of [`ELEMENT_INDEX`](@ref)
  - `mass::Bool=true`: If the x axis will be the stellar mass density or the SFR density.
  - `output_path::String="./resolvedKennicuttSchmidtLaw"`: Path to the output folder.
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
function resolvedMassMetallicityRelation(
    simulation_paths::Vector{String},
    slice::IndexType;
    element::Symbol=:all,
    mass::Bool=true,
    output_path::String="./resolvedMassMetallicityRelation",
    filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all,
    sim_labels::Union{Vector{String},Nothing}=basename.(simulation_paths),
    theme::Attributes=Theme(),
)::Nothing

    grid = CubicGrid(65.0u"kpc", 300)

    (
        element ∈ [:all, keys(ELEMENT_INDEX)...] ||
        throw(ArgumentError("resolvedMassMetallicityRelation: `quantity` can only be :all or any \
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
        da_kwargs=[(; filter_function=dd->filterStellarAge(dd))],
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

    with_theme(merge(theme, theme_latexfonts(), DEFAULT_THEME)) do

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

end
