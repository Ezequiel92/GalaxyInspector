####################################################################################################
# Constants and data structures
####################################################################################################

####################################################################################################
# General
####################################################################################################

"""
If logging messages will be printed out.
"""
const LOGGING = Ref{Bool}(false)

"""
IO stream for logging messages.
"""
const LOG_STREAM = Ref{Union{IO,Nothing}}(nothing)

"""
Path to the directory where memory-mapped arrays will be stored.
"""
const MMAP_PATH = joinpath(@__DIR__, "../../mmap")

####################################################################################################
# Reference values from the literature
####################################################################################################

@doc raw"""
Solar abundances.

These are defined as $12 + \log_{10}(N_\mathrm{X} / N_\mathrm{H})$, where $N_\mathrm{X}$ and $N_\mathrm{H}$ are the number densities of element $\mathrm{X}$ and hydrogen respectively.

# References

M. Asplund et al. (2009). *The Chemical Composition of the Sun*. Annual Review of Astronomy and Astrophysics, **47(1)**, 481–522. [doi:10.1146/annurev.astro.46.060407.145222](https://doi.org/10.1146/annurev.astro.46.060407.145222)
"""
const SOLAR_ABUNDANCE = Dict(
    :H  => 12,    # Hydrogen
    :He => 10.93, # Helium
    :C  => 8.43,  # Carbon
    :N  => 7.83,  # Nitrogen
    :O  => 8.69,  # Oxygen
    :Ne => 7.93,  # Neon
    :Mg => 7.60,  # Magnesium
    :Si => 7.51,  # Silicon
    :S  => 7.12,  # Sulfur
    :Ca => 6.34,  # Calcium
    :Fe => 7.50,  # Iron
)

"""
Standard atomic weights.

# References

T. Prohaska et al. (2022). *Standard atomic weights of the elements 2021 (IUPAC Technical Report)*. Pure and Applied Chemistry, **94(5)**, 573-600. [doi:10.1515/pac-2019-0603](https://doi.org/10.1515/pac-2019-0603)
"""
const ATOMIC_WEIGHTS = Dict(
    :H  => 1.0080, # Hydrogen
    :He => 4.0026, # Helium
    :C  => 12.011, # Carbon
    :N  => 14.007, # Nitrogen
    :O  => 15.999, # Oxygen
    :Ne => 20.180, # Neon
    :Mg => 24.305, # Magnesium
    :Si => 28.085, # Silicon
    :S  => 32.06,  # Sulfur
    :Ca => 40.078, # Calcium
    :Fe => 55.845, # Iron
)

###################
# Kennicutt (1998)
###################

@doc raw"""
Slope of the Kennicutt-Schmidt law, taken from Kennicutt (1998) (Section 4, Equation 4).

```math
\Sigma_\mathrm{SFR} = a \left( \frac{\Sigma_\mathrm{gas}}{1 \, \mathrm{M_\odot \, pc^{-2}}} \right)^{\!N} \mathrm{M_\odot \, yr^{-1} \, kpc^{-2}} \, ,
```

# References

R. C. Kennicutt (1998). *The Global Schmidt Law in Star-forming Galaxies*. The Astrophysical Journal, **498(2)**, 541-552. [doi:10.1086/305588](https://doi.org/10.1086/305588)
"""
const N_KS98 = 1.4 ± 0.15

@doc raw"""
Intercept of the Kennicutt-Schmidt law, taken from Kennicutt (1998) (Section 4, Equation 4).

```math
\Sigma_\mathrm{SFR} = a \left( \frac{\Sigma_\mathrm{gas}}{1 \, \mathrm{M_\odot \, pc^{-2}}} \right)^{\!N} \mathrm{M_\odot \, yr^{-1} \, kpc^{-2}} \, ,
```

# References

R. C. Kennicutt (1998). *The Global Schmidt Law in Star-forming Galaxies*. The Astrophysical Journal, **498(2)**, 541-552. [doi:10.1086/305588](https://doi.org/10.1086/305588)
"""
const a_KS98 = 2.5e-4 ± 0.7e-4

@doc raw"""
Range of values for

```math
\Sigma_\mathrm{SFR} \, [\mathrm{M_\odot \, yr^{-1} \, kpc^{-2}}] \, ,
```
from the combine data (Table 1 and 2) in Kennicutt (1998).

# References

R. C. Kennicutt (1998). *The Global Schmidt Law in Star-forming Galaxies*. The Astrophysical Journal, **498(2)**, 541-552. [doi:10.1086/305588](https://doi.org/10.1086/305588)
"""
const SFR_RANGE_KS98 = exp10.([-3.55, 2.98]) .* u"Msun * yr^-1 * kpc^-2"

#######################
# Bigiel et al. (2008)
#######################

@doc raw"""
Kennicutt-Schmidt law fits for molecular and neutral gas, from Bigiel et al. (2008) (Table 2, Average).

Power-law index, N, and $A = \log_{10}(a)$, where $a$ is $\Sigma_\mathrm{SFR}$ at the fiducial gas surface density of $10 \, \mathrm{M_\odot \, pc^{-2}}$ are given.

```math
\Sigma_\mathrm{SFR} = a \left( \frac{\Sigma_\mathrm{HI, H_2, gas}}{10 \, \mathrm{M_\odot \, pc^{-2}}} \right)^{\!N} \, ,
```

# References

F. Bigiel et al. (2008). *THE STAR FORMATION LAW IN NEARBY GALAXIES ON SUB-KPC SCALES*. The Astrophysical Journal, **136(6)**, 2846. [doi:10.1088/0004-6256/136/6/2846](https://doi.org/10.1088/0004-6256/136/6/2846)
"""
const A_BIGIEL2008_MOLECULAR = −2.06 ± 0.17
const N_BIGIEL2008_MOLECULAR = 0.96 ± 0.07
const A_BIGIEL2008_NEUTRAL   = −2.39 ± 0.28
const N_BIGIEL2008_NEUTRAL   = 1.85 ± 0.70

@doc raw"""
Kennicutt-Schmidt law best-fit for molecular gas, from Bigiel et al. (2008) (Section 4.3, Equation 3).

Power-law index, N, and $A = \log_{10}(a)$, where $a$ is $\Sigma_\mathrm{SFR}$ at the fiducial gas surface density of $10 \, \mathrm{M_\odot \, pc^{-2}}$ are given.

```math
\Sigma_\mathrm{SFR} = a \left( \frac{\Sigma_\mathrm{H_2}}{10 \, \mathrm{M_\odot \, pc^{-2}}} \right)^{\!N} \, ,
```

# References

F. Bigiel et al. (2008). *THE STAR FORMATION LAW IN NEARBY GALAXIES ON SUB-KPC SCALES*. The Astrophysical Journal, **136(6)**, 2846. [doi:10.1088/0004-6256/136/6/2846](https://doi.org/10.1088/0004-6256/136/6/2846)
"""
const A_BIGIEL2008_BF_MOLECULAR = −2.1 ± 0.2
const N_BIGIEL2008_BF_MOLECULAR = 1.0 ± 0.2

@doc raw"""
Range of values for

```math
\Sigma_\mathrm{SFR} \, [\mathrm{M_\odot \, yr^{-1} \, kpc^{-2}}] \, ,
```
in the seven spirals in Table 1 of Bigiel et al. (2008), with associated molecular data.

The actual values for the SFR density are taken from Table 2 in Bigiel et al. (2010), using only the ones with associated molecular data.

# References

F. Bigiel et al. (2008). *THE STAR FORMATION LAW IN NEARBY GALAXIES ON SUB-KPC SCALES*. The Astrophysical Journal, **136(6)**, 2846. [doi:10.1088/0004-6256/136/6/2846](https://doi.org/10.1088/0004-6256/136/6/2846)

F. Bigiel et al. (2010). *EXTREMELY INEFFICIENT STAR FORMATION IN THE OUTER DISKS OF NEARBY GALAXIES*. The Astrophysical Journal, **140(5)**, 1194. [doi:10.1088/0004-6256/140/5/1194](https://doi.org/10.1088/0004-6256/140/5/1194)
"""
const SFR_RANGE_BIGIEL2008 = exp10.([-2.99, -0.33]) .* u"Msun * yr^-1 * kpc^-2"

"""
Spatial resolution used in Bigiel et al. (2008) (see Section 1).

# References

F. Bigiel et al. (2008). *THE STAR FORMATION LAW IN NEARBY GALAXIES ON SUB-KPC SCALES*. The Astrophysical Journal, **136(6)**, 2846. [doi:10.1088/0004-6256/136/6/2846](https://doi.org/10.1088/0004-6256/136/6/2846)
"""
const PXSIZE_BIGIEL2008 = 750.0u"pc"

#######################
# Bigiel et al. (2010)
#######################

"""
Path to Table 2 of Bigiel et al. (2010).

# References

F. Bigiel et al. (2010). *EXTREMELY INEFFICIENT STAR FORMATION IN THE OUTER DISKS OF NEARBY GALAXIES*. The Astrophysical Journal, **140(5)**, 1194. [doi:10.1088/0004-6256/140/5/1194](https://doi.org/10.1088/0004-6256/140/5/1194)
"""
const BIGIEL2010_TABLE_2 = joinpath(@__DIR__, "../../experimental_data/bigiel_2010/table_02.txt")

"""
Path to Table 3 of Bigiel et al. (2010).

# References

F. Bigiel et al. (2010). *EXTREMELY INEFFICIENT STAR FORMATION IN THE OUTER DISKS OF NEARBY GALAXIES*. The Astrophysical Journal, **140(5)**, 1194. [doi:10.1088/0004-6256/140/5/1194](https://doi.org/10.1088/0004-6256/140/5/1194)
"""
const BIGIEL2010_TABLE_3 = joinpath(@__DIR__, "../../experimental_data/bigiel_2010/table_03.txt")

########################
# Cologni et al. (2026)
########################

@doc raw"""
Kennicutt-Schmidt law best-fit for molecular gas, from Cologni et al. (2026) (Section 4.3.2).

Power-law index, N, and $A = \log_{10}(a)$, where $a$ is $\Sigma_\mathrm{SFR}$ at the fiducial molecular gas surface density of $1.0 \mathrm{M_\odot \, kpc^{-2}}$ are given.

```math
\Sigma_\mathrm{SFR} = a \left( \frac{\Sigma_\mathrm{H_2}}{1.0 \, \mathrm{M_\odot \, kpc^{-2}}} \right)^{\!N} \, ,
```

# References

R. Cologni et al. (2026). *Resolved molecular gas and star formation in massive unquenched spirals*. Astronomy and Astrophysics, **709**, A148. [doi:10.1051/0004-6361/202558486](https://doi.org/10.1051/0004-6361/202558486)
"""
const A_COLOGNI2026 = -8.36 ± 0.57
const N_COLOGNI2026 = 0.87 ± 0.09

@doc raw"""
Range of values for

```math
\Sigma_\mathrm{SFR} \, [\mathrm{M_\odot \, yr^{-1} \, kpc^{-2}}] \, ,
```
in Figure 7 of Cologni et al. (2026), with associated molecular data.

# References

R. Cologni et al. (2026). *Resolved molecular gas and star formation in massive unquenched spirals*. Astronomy and Astrophysics, **709**, A148. [doi:10.1051/0004-6361/202558486](https://doi.org/10.1051/0004-6361/202558486)
"""
const SFR_RANGE_COLOGNI2026 = exp10.([-3.5, -2.0]) .* u"Msun * yr^-1 * kpc^-2"

"""
Spatial resolution used in Cologni et al. (2026) (see Section 3.2).

# References

R. Cologni et al. (2026). *Resolved molecular gas and star formation in massive unquenched spirals*. Astronomy and Astrophysics, **709**, A148. [doi:10.1051/0004-6361/202558486](https://doi.org/10.1051/0004-6361/202558486)
"""
const PXSIZE_COLOGNI2026 = 3.0u"kpc"

####################
# Lin et al. (2019)
####################

@doc raw"""
Kennicutt-Schmidt law best-fit for molecular gas, from Lin et al. (2019) (Section 3.1, Table 1).

Power-law index, N, and $A = \log_{10}(a)$, where $a$ is $\Sigma_\mathrm{SFR}$ at the fiducial molecular gas surface density of $1.0 \mathrm{M_\odot \, kpc^{-2}}$ are given.

```math
\Sigma_\mathrm{SFR} = a \left( \frac{\Sigma_\mathrm{H_2}}{1.0 \, \mathrm{M_\odot \, kpc^{-2}}} \right)^{\!N} \, ,
```

# References

L. Lin et al. (2019). *The ALMaQUEST Survey: The Molecular Gas Main Sequence and the Origin of the Star-forming Main Sequence*. The Astrophysical Journal Letters, **884(2)**, L33. [doi:10.3847/2041-8213/ab4815](https://doi.org/10.3847/2041-8213/ab4815)
"""
const A_LIN2019 = −9.33 ± 0.06
const N_LIN2019 = 1.05 ± 0.01

@doc raw"""
Range of values for

```math
\Sigma_\mathrm{SFR} \, [\mathrm{M_\odot \, yr^{-1} \, kpc^{-2}}] \, ,
```
in Figure 3 of Lin et al. (2019), with associated molecular data.

# References

L. Lin et al. (2019). *The ALMaQUEST Survey: The Molecular Gas Main Sequence and the Origin of the Star-forming Main Sequence*. The Astrophysical Journal Letters, **884(2)**, L33. [doi:10.3847/2041-8213/ab4815](https://doi.org/10.3847/2041-8213/ab4815)
"""
const SFR_RANGE_LIN2019 = exp10.([-3.5, -0.5]) .* u"Msun * yr^-1 * kpc^-2"

"""
Spatial resolution used in Lin et al. (2019) (see Section 4).

# References

L. Lin et al. (2019). *The ALMaQUEST Survey: The Molecular Gas Main Sequence and the Origin of the Star-forming Main Sequence*. The Astrophysical Journal Letters, **884(2)**, L33. [doi:10.3847/2041-8213/ab4815](https://doi.org/10.3847/2041-8213/ab4815)
"""
const PXSIZE_LIN2019 = 1.0u"kpc"

##########################
# Querejeta et al. (2021)
##########################

@doc raw"""
Kennicutt-Schmidt law best-fit for molecular gas, from Querejeta et al. (2021) (Section 4.3, Table 4, row All, mean of the columns).

Power-law index, N, and $A = \log_{10}(a)$, where $a$ is $\Sigma_\mathrm{SFR}$ at the fiducial molecular gas surface density of $1.0 \mathrm{M_\odot \, pc^{-2}}$ are given.

```math
\Sigma_\mathrm{SFR} = a \left( \frac{\Sigma_\mathrm{H_2}}{1.0 \, \mathrm{M_\odot \, pc^{-2}}} \right)^{\!N} \, ,
```

# References

M. Querejeta et al. (2021). *Stellar structures, molecular gas, and star formation across the PHANGS sample of nearby galaxies*. Astronomy and Astrophysics, **656**, A133. [doi:10.1051/0004-6361/202140695](https://doi.org/10.1051/0004-6361/202140695)
"""
const A_QUEREJETA2021 = -3.285 ± 0.012
const N_QUEREJETA2021 = 1.1 ± 0.075

"""
Spatial resolution used in Querejeta et al. (2021) (see Section 2.3).

# References

M. Querejeta et al. (2021). *Stellar structures, molecular gas, and star formation across the PHANGS sample of nearby galaxies*. Astronomy and Astrophysics, **656**, A133. [doi:10.1051/0004-6361/202140695](https://doi.org/10.1051/0004-6361/202140695)
"""
const PXSIZE_QUEREJETA2021 = 1.5u"kpc"

@doc raw"""
Range of values for

```math
\Sigma_\mathrm{SFR} \, [\mathrm{M_\odot \, yr^{-1} \, kpc^{-2}}] \, ,
```
in Figure 6 and 7 of Querejeta et al. (2021), with associated molecular data.

# References

M. Querejeta et al. (2021). *Stellar structures, molecular gas, and star formation across the PHANGS sample of nearby galaxies*. Astronomy and Astrophysics, **656**, A133. [doi:10.1051/0004-6361/202140695](https://doi.org/10.1051/0004-6361/202140695)
"""
const SFR_RANGE_QUEREJETA2021 = exp10.([-4.0, 0.0]) .* u"Msun * yr^-1 * kpc^-2"

##########################
# Giannetti et al. (2017)
##########################

@doc raw"""
Best fit for the gas-to-dust ratio vs. galactocentric distance, from Giannetti et al. (2017) (Section 3, Equation 2).

```math
\log_{10} \, \gamma = A \, R_\text{GC} + B \, ,
```

# References

A. Giannetti et al. (2017). *Galactocentric variation of the gas-to-dust ratio and its relation with metallicity*. Astronomy and Astrophysics, **606**, L12. [doi:10.1051/0004-6361/201731728](https://doi.org/10.1051/0004-6361/201731728)
"""
const A_GIANNETTI2017 = 0.087u"kpc^-1"
const B_GIANNETTI2017 = 1.44

@doc raw"""
Range of values for the galactocentric distances in Figure 4 and Table A.1 of Giannetti et al. (2017).

# References

A. Giannetti et al. (2017). *Galactocentric variation of the gas-to-dust ratio and its relation with metallicity*. Astronomy and Astrophysics, **606**, L12. [doi:10.1051/0004-6361/201731728](https://doi.org/10.1051/0004-6361/201731728)
"""
const RGC_RANGE_GIANNETTI2017 = [2.0, 20.3] .* u"kpc"

####################
# Sun et al. (2023)
####################

"""
Spatial resolution used in Sun et al. (2023).

# References

J. Sun et al. (2023). *Star Formation Laws and Efficiencies across 80 Nearby Galaxies*. The Astrophysical Journal Letters, **945(2)**, L19. [doi:10.3847/2041-8213/acbd9c](https://doi.org/10.3847/2041-8213/acbd9c)
"""
const PXSIZE_SUN2023 = 1.5u"kpc"

##############
# Data tables
##############

"""
Path to Table A1 from Sun et al. (2023).

# References

J. Sun et al. (2023). *Star Formation Laws and Efficiencies across 80 Nearby Galaxies*. The Astrophysical Journal Letters, **945(2)**, L19. [doi:10.3847/2041-8213/acbd9c](https://doi.org/10.3847/2041-8213/acbd9c)
"""
const SUN2023_DATA_PATH = joinpath(@__DIR__, "../../experimental_data/sun_2023.txt")

"""
Path to Table 4 (corrected) from de los Reyes et al. (2019).

# References

M. A. C. de los Reyes et al. (2019). *Revisiting the Integrated Star Formation Law. I. Non-starbursting Galaxies*. The Astrophysical Journal, **872(1)**, 16. [doi:10.3847/1538-4357/aafa82)](https://doi.org/10.3847/1538-4357/aafa82)

M. A. C. de los Reyes et al. (2019). *Erratum: “Revisiting the Integrated Star Formation Law. I. Non-starbursting Galaxies” (2019 ApJ, 872, 16)*. The Astrophysical Journal, **878(1)**, 74. [doi:10.3847/1538-4357/ab22af)](https://doi.org/10.3847/1538-4357/ab22af)
"""
const DELOSREYES2019_DATA_PATH = joinpath(@__DIR__, "../../experimental_data/de_los_reyes_2019.txt")

"""
Path to the file with the Milky Way profiles from Mollá et al. (2015).

# References

M. Mollá et al. (2015). *Galactic chemical evolution: stellar yields and the initial mass function*. Monthly Notices of the Royal Astronomical Society **451(4)**, 3693–3708. [doi:10.1093/mnras/stv1102](https://doi.org/10.1093/mnras/stv1102)
"""
const MOLLA2015_DATA_PATH = joinpath(@__DIR__, "../../experimental_data/molla_2015.csv")

"""
Path to the file with the global galactic properties from Feldmann (2020).

# References

R. Feldmann (2020). *The link between star formation and gas in nearby galaxies*. Communications Physics **3(226)**. [doi:10.1038/s42005-020-00493-0](https://doi.org/10.1038/s42005-020-00493-0)
"""
const FELDMANN2020_DATA_PATH = joinpath(@__DIR__, "../../experimental_data/feldmann_2020.csv")

"""
Path to the file with the profiles from Leroy et al. (2008).

# References

A. K. Leroy et al. (2008). *THE STAR FORMATION EFFICIENCY IN NEARBY GALAXIES: MEASURING WHERE GAS FORMS STARS EFFECTIVELY*. The Astronomical Journal **136(6)**, 2782–2845. [doi:10.1088/0004-6256/136/6/2782](https://doi.org/10.1088/0004-6256/136/6/2782)

"""
const LEROY2008_DATA_PATH = joinpath(@__DIR__, "../../experimental_data/leroy_2008.jld2")

"""
Path to the file with the fits from McMillan (2011).

# References

P. J. McMillan (2011). *Mass models of the Milky Way*. Monthly Notices of the Royal Astronomical Society **414(3)**, 2446-2457. [doi:10.1111/j.1365-2966.2011.18564.x](https://doi.org/10.1111/j.1365-2966.2011.18564.x)
"""
const MCMILLAN2011_DATA_PATH = joinpath(@__DIR__, "../../experimental_data/mcmillan_2011.jld2")

"""
Path to the file with the stellar magnitudes for the g, r, and i filters (SDSS AB system) from Millán-Irigoyen et al. (2025).

The effective wavelengths of each filter are taken from Table 2 (2.5m reference) of Doi et al. (2010):

g band -> λeff = 4627 Å -> "blue" channel
r band -> λeff = 6140 Å -> "green" channel
i band -> λeff = 7467 Å -> "red" channel

# References

I. Millán-Irigoyen et al. (2025). *HR-pyPopStar II: high spectral resolution evolutionary synthesis models low metallicity expansion and the properties of the stellar populations of dwarf galaxies*. arXiv. [doi:10.48550/arxiv.2510.02886](https://doi.org/10.48550/arxiv.2510.02886)

M. Doi et al. (2010). *PHOTOMETRIC RESPONSE FUNCTIONS OF THE SLOAN DIGITAL SKY SURVEY IMAGER*. The Astronomical Journal **139(4)**, 1628-1648. [doi:10.1088/0004-6256/139/4/1628](https://doi.org/10.1088/0004-6256/139/4/1628)
"""
const MILLANIRIGOYEN2025_DATA_PATH = joinpath(@__DIR__, "../../experimental_data/millan-irigoyen_2025.txt")

####################################################################################################

@doc raw"""
Extinction factor from T. Güver et al. (2009).

```math
\text{Extinction factor} = \frac{N_\mathrm{H}}{A(V)} = (2.21 \pm 0.09) \times 10^{21} \, \mathrm{atoms \, cm^{-2} \, mag^{-1}} \, .
```

We assume that an atom mass is the proton mass.

T. Güver et al. (2009). *The relation between optical extinction and hydrogen column density in the Galaxy*. Monthly Notices of the Royal Astronomical Society **400(4)**, 2050-2053. [doi:10.1111/j.1365-2966.2009.15598.x](https://doi.org/10.1111/j.1365-2966.2009.15598.x)
"""
const EXTINCTION_FACTOR_NH = 17.7006u"Msun*pc^-2"

@doc raw"""
Extinction factor from Nozawa et al. (2013) eq. 12.

```math
\text{Extinction factor}^{-1} = \frac{A(V)}{\Sigma_\text{dust}} = (3.7 \pm 0.5) \times 10^{4} \, \mathrm{g^{-1} \, cm^2 \, mag} \, .
```

T. Nozawa et al. (2013). *PROPERTIES OF DUST GRAINS PROBED WITH EXTINCTION CURVES*. The Astrophysical Journal **770(1)**, 27. [doi:10.1088/0004-637X/770/1/27](https://doi.org/10.1088/0004-637X/770/1/27)
"""
const EXTINCTION_FACTOR_D = 0.171u"Msun*pc^-2"

@doc raw"""
Relation between the extinction in the g, r, and i bands of SDSS and the extinction in the V band from Table 3 of S. Wang et al. (2019).

```math
\left. A_\lambda / A_V \right\vert_g = 1.205 \pm 0.010 \, .
```

```math
\left. A_\lambda / A_V \right\vert_r = 0.848 \pm 0.006 \, .
```

```math
\left. A_\lambda / A_V \right\vert_i = 0.630 \pm 0.004 \, .
```

S. Wang et al. (2019). *The Optical to Mid-infrared Extinction Law Based on the APOGEE, Gaia DR2, Pan-STARRS1, SDSS, APASS, 2MASS, and WISE Surveys*. The American Astronomical Society **877(2)**, 116. [doi:10.3847/1538-4357/ab1c61](https://doi.org/10.3847/1538-4357/ab1c61)
"""
const AλAV_g = 1.205
const AλAV_r = 0.848
const AλAV_i = 0.63

"""
Reference pressure for the molecular fraction-pressure relation, from Blitz et al. (2006) (Table 2, "Mean" row, Third column).

# References

L. Blitz et al. (2006). *The Role of Pressure in GMC Formation II: The H2-Pressure Relation*. The Astrophysical Journal, **650(2)**, 933. [doi:10.1086/505417](https://doi.org/10.1086/505417)
"""
const P0 = 3.5e4u"K*cm^-3" * Unitful.k

"""
Reference exponent for the molecular fraction-pressure relation, from Blitz et al. (2006) (Table 2, "Mean" row, Second column).

We use -α here.

# References

L. Blitz et al. (2006). *The Role of Pressure in GMC Formation II: The H2-Pressure Relation*. The Astrophysical Journal, **650(2)**, 933. [doi:10.1086/505417](https://doi.org/10.1086/505417)
"""
const ALPHA_BLITZ = -0.92

"""
Fitted parameters for the SFR of main-sequence galaxies as a function of stellar mass and redshift, from Schreiber et al. (2015) (Equation 9).

C. Schreiber et al. (2015). *The Herschel view of the dominant mode of galaxy growth from z = 4 to the present day*. Astronomy & Astrophysics, **575**, A74. [doi:10.1051/0004-6361/201425017](https://doi.org/10.1051/0004-6361/201425017)
"""
const SCHREIBER2015_m0 = 0.5 ± 0.07
const SCHREIBER2015_a0 = 1.5 ± 0.15
const SCHREIBER2015_a1 = 0.3 ± 0.08
const SCHREIBER2015_m1 = 0.36 ± 0.3
const SCHREIBER2015_a2 = 2.5 ± 0.6

"""
Fitted parameters for the gas mass fraction as a function of stellar mass, redshift, and sSFR, from Scoville et al. (2016) (Equation 1).

N. Scoville et al. (2016). *ISM MASSES AND THE STAR FORMATION LAW AT Z = 1 TO 6: ALMA OBSERVATIONS OF DUST CONTINUUM IN 145 GALAXIES IN THE COSMOS SURVEY FIELD*. The Astrophysical Journal, **820(2)**, 83. [doi:10.3847/0004-637X/820/2/83](https://doi.org/10.3847/0004-637X/820/2/83)
"""
const SCOVILLE2016_α0 = 0.3 ± 0.02
const SCOVILLE2016_α1 = -0.02 ± 0.02
const SCOVILLE2016_α2 = 0.44 ± 0.05
const SCOVILLE2016_α3 = 0.32 ± 0.02

"""
    LEE2015_PARAMETERS(z::Float64)::NTuple{3,Measurements.Measurement{Float64}}

Best-fit parameters for the SFR as a function of stellar mass, from Lee et al. (2015) (Table 1).

# Arguments

    - `z :: Float64`: Redshift. The parameters are only valid for 0.25 <= z <= 1.3, but the function will return the values for z = 0.25 or z = 1.3 if the input value is outside this range.

# Returns

  - A tuple with three elements (See Section 4.1 of Lee et al. (2015) for the definition of the parameters):

      + S0 [log10(M⊙ * yr^-1)]: The maximum value of log10(SFR / M⊙ * yr^-1).
      + M0 [log10(M⊙)]: Turnover mass.
      + γ [dimensionless]: Power-law slope at low stellar masses.

# References

N. Lee et al. (2015). *A TURNOVER IN THE GALAXY MAIN SEQUENCE OF STAR FORMATION AT M* ~ 10^10 M☉ FOR REDSHIFTS z < 1.3*. The Astrophysical Journal, **801(2)**, 80. [doi:10.1088/0004-637X/801/2/80](https://doi.org/10.1088/0004-637X/801/2/80)
"""
function LEE2015_PARAMETERS(z::Float64)::NTuple{3,Measurements.Measurement{Float64}}

    if z < 0.25 || z > 1.3
        (
            LOGGING[] &&
            @warn("The parameters from Lee et al. (2015) are only valid for 0.25 <= z <= 1.3. The returned values will be the same as for z = 0.25 or z = 1.3, respectively.")
        )
    end

    if z < 0.46
        S0 = 0.8 ± 0.019
        M0 = 10.03 ± 0.042
        γ  = 0.92 ± 0.017
    elseif z < 0.63
        S0 = 0.99 ± 0.015
        M0 = 9.82 ± 0.031
        γ  = 1.13 ± 0.033
    elseif z < 0.78
        S0 = 1.23 ± 0.016
        M0 = 9.93 ± 0.031
        γ  = 1.11 ± 0.025
    elseif z < 0.93
        S0 = 1.35 ± 0.014
        M0 = 9.96 ± 0.025
        γ  = 1.28 ± 0.034
    elseif z < 1.11
        S0 = 1.53 ± 0.017
        M0 = 10.10 ± 0.029
        γ  = 1.26 ± 0.032
    else
        S0 = 1.72 ± 0.024
        M0 = 10.31 ± 0.043
        γ  = 1.07 ± 0.028
    end

    return S0, M0, γ

end

"""
Median star formation efficiency per free-fall time, from Evans et al. (2014) (Section 1.1).

N. J. Evans II et al. (2014). *STAR FORMATION RELATIONS IN NEARBY MOLECULAR CLOUDS*. The Astrophysical Journal, **782(2)**, 114. [doi:10.1088/0004-637X/782/2/114](https://doi.org/10.1088/0004-637X/782/2/114)
"""
const EVANS2014_ϵff = 0.016 ± 0.013

"""
Median star formation efficiency per free-fall time, from Lee et al. (2016) (Table 3 and eq. 15).

E. J. Lee et al. (2016). *OBSERVATIONAL EVIDENCE OF DYNAMIC STAR FORMATION RATE IN MILKY WAY GIANT MOLECULAR CLOUDS*. The Astrophysical Journal, **833(2)**, 229. [doi:10.3847/1538-4357/833/2/229](https://doi.org/10.3847/1538-4357/833/2/229)
"""
const LEE2016_ϵff = 0.018 ± 0.039

####################################################################################################
# Type aliases and structs
####################################################################################################

"""
Type for colors.
"""
const ColorType = Union{ColorTypes.RGB,ColorTypes.RGBA,Symbol}

"""
Type for line styles.
"""
const LineStyleType = Union{Tuple{String,Symbol},Nothing,String,Symbol}

"""
Type for numbers.
"""
const FloatOrQuantity = Union{AbstractFloat,Unitful.AbstractQuantity{<:AbstractFloat}}

"""
Type for translations.
"""
const TranslationType = Union{
    Tuple{Vector{<:Unitful.Length},Vector{<:Unitful.Velocity}},
    Symbol,
    Int,
    NTuple{2,Int},
}

"""
Type for rotations.
"""
const RotationType = Union{
    Tuple{Symbol,Symbol,Function},
    Tuple{Matrix{Float64}},
    Tuple{UniformScaling{Bool}},
}

"""
Index type.
"""
const IndexType = Union{
    Colon,
    Integer,
    UnitRange{<:Integer},
    StepRange{<:Integer,<:Integer},
    Vector{<:Integer},
    InvertedIndex,
}

"""
Reduced index type.
"""
const ReducedIndexType = Union{
    Integer,
    UnitRange{<:Integer},
    StepRange{<:Integer,<:Integer},
    Vector{<:Integer},
    InvertedIndex,
}

# Dimensions of specific energy
@derived_dimension SpecificEnergy Unitful.𝐋^2 * Unitful.𝐓^-2 true

# Dimensions of surface density
@derived_dimension SurfaceDensity Unitful.𝐌 * Unitful.𝐋^-2 true

# Dimensions of mass flow surface density
@derived_dimension MassFlowDensity Unitful.𝐌 * Unitful.𝐓^-1 * Unitful.𝐋^-2 true

# Dimensions of angular momentum
@derived_dimension AngularMomentum Unitful.𝐌 * Unitful.𝐋^2 * Unitful.𝐓^-1 true

# Dimensions of number density
@derived_dimension NumberDensity Unitful.𝐋^-3 true

"""
Abstract derived quantity type.
"""
abstract type AbstractPlotQuantity end

"""
Derived quantity type.

# Fields

  - `id             :: Symbol`: Symbol representing the quantity.
  - `qty_label      :: AbstractString`: Name of the quantity for the axis label.
  - `request        :: Dict{Symbol,Vector{String}} = Dict{Symbol,Vector{String}}()`: Data request for [`readSnapshot`](@ref). It must have the shape `cell/particle type` -> [`block`, `block`, ...], where the possible types are the keys of [`PARTICLE_INDEX`](@ref), :group, and :subhalo, and the possible blocks are the keys of [`QUANTITIES`](@ref).
  - `unit           :: Unitful.Units = Unitful.NoUnits`: Fiducial unit for the quantity.
  - `exp_factor     :: Int = 0`: Numerical exponent to scale down the quantity, e.g. if `exp_factor` = 10 the values will be divided by ``10^{10}``.
  - `cp_type        :: Union{Symbol,Nothing} = nothing`: Cell/particle type corresponding to the quantity. It can be any of the keys of [`PARTICLE_INDEX`](@ref).
  - `scatter_func   :: Function = getNothing`: Function that computes the quantity for each cell/particle given a data dictionary (see [`makeDataDict`](@ref) for the canonical description). It should return nothing if it does not apply to the quantity.
  - `integrate_func :: Function = getNothing`: Function that computes the quantity for all cells/particles as a single value, given a data dictionary (see [`makeDataDict`](@ref) for the canonical description). It should return nothing if it does not apply to the quantity.
"""
@kwdef struct BaseQuantity <: AbstractPlotQuantity
    id             :: Symbol
    qty_label      :: AbstractString
    request        :: Dict{Symbol,Vector{String}} = Dict{Symbol,Vector{String}}()
    unit           :: Unitful.Units = Unitful.NoUnits
    exp_factor     :: Int = 0
    cp_type        :: Union{Symbol,Nothing} = nothing
    scatter_func   :: Function = getNothing
    integrate_func :: Function = getNothing
end

"""
Represents a mathematical operation between two derived quantities.

# Fields

  - `op       ::Function`: Binary operation, e.g. +, -, *, /, etc.
  - `op_label ::String`: The binary operation as a string, for the axis label, e.g. "+", "/", etc.
  - `left     ::AbstractPlotQuantity`: Left quantity.
  - `right    ::AbstractPlotQuantity`: Right quantity.
"""
struct BinaryQuantity <: AbstractPlotQuantity
    op       :: Function
    op_label :: String
    left     :: AbstractPlotQuantity
    right    :: AbstractPlotQuantity
end

"""
Dimensional information about a physical quantity.

# Fields

  - `hdf5_name  :: String`: HDF5 block name.
  - `dimensions :: Unitful.Dimensions`: Physical dimensions of the quantity, e.g. `Unitful.𝐋 * Unitful.𝐓^-1`.
  - `unit       :: Union{Unitful.Units,Symbol}`: Units of the quantity within the simulation code. It can be a unit from [Unitful](https://github.com/PainterQubits/Unitful.jl) or [UnitfulAstro](https://github.com/JuliaAstro/UnitfulAstro.jl), or it can be the symbol `:internal` which denotes internal code units.
"""
struct SnapshotBlock
    hdf5_name  :: String
    dimensions :: Unitful.Dimensions
    unit       :: Union{Unitful.Units,Symbol}
end

"""
Data in the "Header" group of a HDF5 snapshot file.

# Fields

  - `box_size   :: Float64`: Total size of the simulation box, in internal units of length.
  - `h          :: Float64`: Dimensionless Hubble parameter, "little h".
  - `mass_table :: Vector{Float64}`: Masses of particle/cell types which have a constant mass, in internal units of mass.
  - `num_files  :: Int32`: Number of file chunks per snapshot.
  - `num_part   :: Vector{Int32}`: Number of particles/cells (of each type) included in this file chunk.
  - `num_total  :: Vector{UInt32}`: Total number of particles/cells (of each type) for this snapshot.
  - `omega_0    :: Float64`: The cosmological density parameter for matter.
  - `omega_l    :: Float64`: The cosmological density parameter for the cosmological constant.
  - `redshift   :: Float64`: The redshift.
  - `time       :: Float64`: The physical time or the scale factor, depending on the type of simulation.
  - `l_unit     :: Unitful.Length`: Conversion factor from internal units of length to centimeters.
  - `m_unit     :: Unitful.Mass`: Conversion factor from internal units of mass to grams.
  - `v_unit     :: Unitful.Velocity`: Conversion factor from internal units of velocity to centimeters per second.
"""
@kwdef struct SnapshotHeader
    box_size   :: Float64
    h          :: Float64
    mass_table :: Vector{Float64}
    num_files  :: Int32
    num_part   :: Vector{Int32}
    num_total  :: Vector{UInt32}
    omega_0    :: Float64
    omega_l    :: Float64
    redshift   :: Float64
    time       :: Float64
    l_unit     :: Unitful.Length
    m_unit     :: Unitful.Mass
    v_unit     :: Unitful.Velocity
end

"""
Data in the "Header" group of a HDF5 group catalog file.

The default values are for when there are no group catalog files.

# Fields

  - `box_size          :: Float64 = NaN`: Total size of the simulation box, in internal units of length.
  - `h                 :: Float64 = NaN`: Dimensionless Hubble parameter, "little h".
  - `n_groups_part     :: Int32 = -1`: Number of halos (FoF groups) in this file chunk.
  - `n_groups_total    :: Int32 = -1`: Total number of halos (FoF groups) in this snapshot.
  - `n_subgroups_part  :: Int32 = -1`: Number of subhalos (subfind) in this file chunk.
  - `n_subgroups_total :: Int32 = -1`: Total number of subhalos (subfind) in this snapshot.
  - `num_files         :: Int32 = -1`: Number of file chunks per snapshot.
  - `omega_0           :: Float64 = NaN`: The cosmological density parameter for matter.
  - `omega_l           :: Float64 = NaN`: The cosmological density parameter for the cosmological constant.
  - `redshift          :: Float64 = NaN`: The redshift.
  - `time              :: Float64 = NaN`: The physical time or the scale factor, depending on the type of simulation.
"""
@kwdef struct GroupCatHeader
    box_size          :: Float64 = NaN
    h                 :: Float64 = NaN
    n_groups_part     :: Int32 = -1
    n_groups_total    :: Int32 = -1
    n_subgroups_part  :: Int32 = -1
    n_subgroups_total :: Int32 = -1
    num_files         :: Int32 = -1
    omega_0           :: Float64 = NaN
    omega_l           :: Float64 = NaN
    redshift          :: Float64 = NaN
    time              :: Float64 = NaN
end

"""
Metadata for a simulation.

# Fields

  - `path             :: String`: Full path to the simulation directory.
  - `index            :: Int`: An index associated with the simulation.
  - `slice            :: IndexType`: Slice of the simulation, i.e. which snapshots will be read. It can be an integer (a single snapshot), a vector of integers (several snapshots), an `UnitRange` (e.g. 5:13), an `StepRange` (e.g. 5:2:13) or (:) (all snapshots).
  - `cosmological     :: Bool`: If the simulation is cosmological,

      + `false` -> Newtonian simulation    (`ComovingIntegrationOn` = 0).
      + `true`  -> Cosmological simulation (`ComovingIntegrationOn` = 1).
  - `simulation_table :: DataFrame`: A dataframe where each row is a snapshot, and the colums are

      + `:ids`            -> Dataframe index, i.e. if there are 10 snapshots in total it runs from 1 to 10.
      + `:numbers`        -> Number in the file name.
      + `:scale_factors`  -> Scale factor.
      + `:redshifts`      -> Redshift.
      + `:physical_times` -> Physical time since the Big Bang.
      + `:lookback_times` -> Physical time left to reach the last snapshot.
      + `:snapshot_paths` -> Full path to each snapshots.
      + `:groupcat_paths` -> Full path to each group catalog files.
"""
struct Simulation
    path             :: String
    index            :: Int
    slice            :: IndexType
    cosmological     :: Bool
    simulation_table :: DataFrame
end

"""
Metadata for a snapshot.

# Fields

  - `path          :: String`: Full path to the snapshot.
  - `global_index  :: Int`: Index of the snapshot in the context of the whole simulation.
  - `slice_index   :: Int`: Index of the snapshot in the context of the simulation slice.
  - `physical_time :: Unitful.Time`: Physical time since the Big Bang.
  - `lookback_time :: Unitful.Time`: Physical time left to reach the last snapshot.
  - `scale_factor  :: Float64`: Scale factor.
  - `redshift      :: Float64`: Redshift.
  - `header        :: SnapshotHeader`: Header.
"""
struct Snapshot
    path          :: String
    global_index  :: Int
    slice_index   :: Int
    physical_time :: Unitful.Time
    lookback_time :: Unitful.Time
    scale_factor  :: Float64
    redshift      :: Float64
    header        :: SnapshotHeader
end

"""
Metadata for a group catalog file.

# Fields

  - `path   :: Union{String,Missing}`: Full path to the group catalog file.
  - `header :: GroupCatHeader`: Header.
"""
struct GroupCatalog
    path   :: Union{String,Missing}
    header :: GroupCatHeader
end

"""
Unit conversion factors.

# Fields

  - `x_cgs::Unitful.Length`: From internal units of length to ``\\mathrm{cm}``.
  - `x_cosmo::Unitful.Length`: From internal units of length to ``\\mathrm{kpc}``.
  - `x_comoving::Unitful.Length`: From internal units of length to ``\\mathrm{ckpc}``.
  - `v_cgs::Unitful.Velocity`: From internal units of velocity to ``\\mathrm{cm \\, s^{-1}}``.
  - `v_cosmo::Unitful.Velocity`: From internal units of velocity to ``\\mathrm{km \\, s^{-1}}``.
  - `m_cgs::Unitful.Mass`: From internal units of mass to ``\\mathrm{g}``.
  - `m_cosmo::Unitful.Mass`: From internal units of mass to ``\\mathrm{M_\\odot}``.
  - `t_cgs::Unitful.Time`: From internal units of time to ``\\mathrm{s}``.
  - `t_cosmo::Unitful.Time`: From internal units of time to ``\\mathrm{Myr}``.
  - `t_newton::Unitful.Time`: From internal units (non-cosmological simulations) of time to ``\\mathrm{Myr}``.
  - `U_cgs::Unitful.Energy`: From internal units of specific energy to ``\\mathrm{erg \\, g^{-1}}``.
  - `rho_cgs::Unitful.Density`: From internal units of density to ``\\mathrm{g \\, cm^{-3}}``.
  - `P_Pa::Unitful.Pressure`: From internal units of pressure to ``\\mathrm{Pa}``.
"""
struct InternalUnits
    x_cgs      :: Unitful.Length   # From internal units of length to cm
    x_cosmo    :: Unitful.Length   # From internal units of length to kpc
    x_comoving :: Unitful.Length   # From internal units of length to ckpc

    v_cgs      :: Unitful.Velocity # From internal units of velocity to cm * s^-1
    v_cosmo    :: Unitful.Velocity # From internal units of velocity to km * s^-1

    m_cgs      :: Unitful.Mass     # From internal units of mass to g
    m_cosmo    :: Unitful.Mass     # From internal units of mass to M⊙

    t_cgs      :: Unitful.Time     # From internal units of time to s
    t_cosmo    :: Unitful.Time     # From internal units of time to Myr
    t_newton   :: Unitful.Time     # From internal units (non-cosmological simulations) of time to Myr

    U_cgs      :: SpecificEnergy   # From internal units of specific energy to erg * g^-1

    rho_cgs    :: Unitful.Density  # From internal units of density to g * cm^-3

    P_Pa       :: Unitful.Pressure # From internal units of pressure to Pa

    """
        InternalUnits(; <keyword arguments>)

    Constructor for `InternalUnits`.

    For cosmological simulations, and when appropriate, the conversion factors turn comoving units into physical units.

    # Arguments

      - `l_unit :: Unitful.Length=INTERNAL_L_UNIT[]`: Code parameter `UnitLength_in_cm`.
      - `m_unit :: Unitful.Mass=INTERNAL_M_UNIT[]`: Code parameter `UnitMass_in_g`.
      - `v_unit :: Unitful.Velocity=INTERNAL_V_UNIT[]`: Code parameter `UnitVelocity_in_cm_per_s`.
      - `a      :: Float64=1.0`: Cosmological scale factor of the simulation.
      - `h      :: Float64=1.0`: Dimensionless Hubble parameter, "little h".
    """
    function InternalUnits(;
        l_unit :: Unitful.Length=INTERNAL_L_UNIT[],
        m_unit :: Unitful.Mass=INTERNAL_M_UNIT[],
        v_unit :: Unitful.Velocity=INTERNAL_V_UNIT[],
        a      :: Float64=1.0,
        h      :: Float64=1.0,
    )

        #############
        # Base units
        #############

        # Length conversion factors
        x_cgs = l_unit * a / h
        x_cosmo = x_cgs |> u"kpc"
        x_comoving = l_unit / h |> u"kpc"

        # Velocity conversion factors
        v_cgs = v_unit * sqrt(a)
        v_cosmo = v_cgs |> u"km * s^-1"

        # Mass conversion factors
        m_cgs = m_unit / h
        m_cosmo = m_cgs |> u"Msun"

        ################
        # Derived units
        ################

        # Time conversion factors
        t_cgs = x_cgs / v_cgs
        t_cosmo = t_cgs |> u"Myr"
        t_newton = 1.0u"Gyr" |> u"Myr"

        # Specific energy conversion factor
        U_cgs = v_unit^2 |> u"erg * g^-1"

        # Density conversion factor
        rho_cgs = m_cgs * x_cgs^-3

        # Thermal pressure conversion factor (it uses v_unit^2 instead of v_cgs^2,
        # which would add an extra factor of a0)
        P_Pa = v_unit^2 * m_cgs * x_cgs^-3 |> u"Pa"

        new(
            x_cgs,
            x_cosmo,
            x_comoving,
            v_cgs,
            v_cosmo,
            m_cgs,
            m_cosmo,
            t_cgs,
            t_cosmo,
            t_newton,
            U_cgs,
            rho_cgs,
            P_Pa,
        )

    end
end

########
# Grids
########

"""
Grid with 1 d.o.f.

Embeddings:

    1D -> Ticks.
    2D -> Concentric rings.
    3D -> Spherical shells.

# Fields

  - `x_axis      :: Vector{<:Number}`: Position of (the center of) the bins, relative to the `origin`.
  - `x_edges     :: Vector{<:Number}`: Position of (the edges of) the bins, relative to the `origin`.
  - `n_bins      :: Int`: Number of bins.
  - `bin_size_1D :: Vector{<:Number}`: Widths of the bins in linear space.
  - `bin_size_2D :: Vector{<:Number}`: Areas of the bins in linear space (assuming a 2D embedding).
  - `bin_size_3D :: Vector{<:Number}`: Volumes of the bins in linear space (assuming a 3D embedding).
  - `origin      :: Vector{<:Number}`: Physical position of the coordinate origin for the grid.
"""
struct LinearGrid
    x_axis      :: Vector{<:Number}
    x_edges     :: Vector{<:Number}
    n_bins      :: Int
    bin_size_1D :: Vector{<:Number}
    bin_size_2D :: Vector{<:Number}
    bin_size_3D :: Vector{<:Number}
    origin      :: Vector{<:Number}

    """
        LinearGrid(
            start  :: Number,
            stop   :: Number,
            n_bins :: Int;
            <keyword arguments>
        )

    Constructor for `LinearGrid`.

    # Arguments

      - `start  :: Number`: Position of the grid left most edge, relative to the `origin`.
      - `stop   :: Number`: Position of the grid right most edge, relative to the `origin`.
      - `n_bins :: Int`: Number of bins.
      - `origin :: Vector{<:Number}=zeros(typeof(start), 3)`: Physical position of the coordinate origin for the grid (passthrough argument).
      - `log    :: Bool=false`: If the grid will be regular in logarithmic space instead of linear space.
    """
    function LinearGrid(
        start  :: Number,
        stop   :: Number,
        n_bins :: Int;
        origin :: Vector{<:Number}=zeros(typeof(start), 3),
        log    :: Bool=false,
    )

        (
            stop > start ||
            throw(ArgumentError("LinearGrid: `stop` must be larger than `start`, \
            but I got `stop` = $(stop) <= `start` = $(start)"))
        )

        (
            isPositive(n_bins) ||
            throw(ArgumentError("LinearGrid: `n_bins` has to be strictly positive, \
            but I got `n_bins` = $(n_bins) <= 0"))
        )

        if log

            (
                isPositive(start) ||
                throw(ArgumentError("LinearGrid: For a logarithmic grid you need a strictly \
                positive `start`, but I got `start` = $(start) <= 0"))
            )

            # Length unit
            u_l = unit(start)

            log_start = log10(ustrip(start))
            log_stop  = log10(ustrip(u_l, stop))

            width = (log_stop - log_start) / n_bins

            x_axis  = exp10.([(i - 0.5) * width + log_start for i in 1:n_bins]) .* u_l
            x_edges = exp10.([i * width + log_start for i in 0:n_bins]) .* u_l

        else

            width = (stop - start) / n_bins

            x_axis  = [(i - 0.5) * width + start for i in 1:n_bins]
            x_edges = [i * width + start for i in 0:n_bins]

        end

        bin_size_1D = [x_edges[i + 1] - x_edges[i] for i in 1:n_bins]

        if start < zero(start)

            shifted_edges = x_edges .- start
            bin_size_2D   = [area(shifted_edges[i + 1]) - area(shifted_edges[i]) for i in 1:n_bins]
            bin_size_3D   = [volume(shifted_edges[i + 1]) - volume(shifted_edges[i]) for i in 1:n_bins]

        else

            bin_size_2D = [area(x_edges[i + 1]) - area(x_edges[i]) for i in 1:n_bins]
            bin_size_3D = [volume(x_edges[i + 1]) - volume(x_edges[i]) for i in 1:n_bins]

        end

        new(x_axis, x_edges, n_bins, bin_size_1D, bin_size_2D, bin_size_3D, origin)

    end

    """
        LinearGrid(
            x_edges :: Vector{<:Number};
            <keyword arguments>
        )

    Constructor for `LinearGrid`.

    # Arguments

      - `x_edges :: Vector{<:Number}`: Relative position of (the edges of) the bins, relative to the `origin`.
      - `origin  :: Vector{<:Number}=zeros(runtimeType(x_edges), 3)`: Physical position of the coordinate origin for the grid (passthrough argument).
      - `log     :: Bool=false`: If the `x_edges` are logarithmic, which means that `x_axis` will mark the center of the bins in logarithmic space instead of linear space.
    """
    function LinearGrid(
        x_edges :: Vector{<:Number};
        origin  :: Vector{<:Number}=zeros(runtimeType(x_edges), 3),
        log     :: Bool=falsem
    )

        isempty(x_edges) && throw(ArgumentError("LinearGrid: `x_edges` is empty!"))

        n_bins = length(x_edges) - 1

        issorted(x_edges) || sort!(x_edges)

        start = first(x_edges)

        bin_size_1D = [x_edges[i + 1] - x_edges[i] for i in 1:n_bins]

        if start < zero(start)

            shifted_edges = x_edges .- start
            bin_size_2D   = [area(shifted_edges[i + 1]) - area(shifted_edges[i]) for i in 1:n_bins]
            bin_size_3D   = [volume(shifted_edges[i + 1]) - volume(shifted_edges[i]) for i in 1:n_bins]

        else

            bin_size_2D = [area(x_edges[i + 1]) - area(x_edges[i]) for i in 1:n_bins]
            bin_size_3D = [volume(x_edges[i + 1]) - volume(x_edges[i]) for i in 1:n_bins]

        end

        if log

            (
                isPositive(start) ||
                throw(ArgumentError("LinearGrid: For a logarithmic grid you need a strictly \
                positive `start`, but I got `start` = $(start) <= 0"))
            )

            # Length unit
            u_l = unit(start)

            log_x_edges = log10.(ustrip.(u_l, x_edges))

            x_axis = exp10.([(log_x_edges[i + 1] + log_x_edges[i]) / 2.0 for i in 1:n_bins]) .* u_l

        else

            x_axis  = [(x_edges[i + 1] + x_edges[i]) / 2.0 for i in 1:n_bins]

        end

        new(x_axis, x_edges, n_bins, bin_size_1D, bin_size_2D, bin_size_3D, origin)

    end
end

"""
Grid with 2 d.o.f.

# Fields

  - `x_axis      :: Vector{<:Number}`: Absolute position of (the center of) the bins in the x direction.
  - `y_axis      :: Vector{<:Number}`: Absolute position of (the center of) the bins in the y direction.
  - `x_edges     :: Vector{<:Number}`: Absolute position of (the edges of) the bins in the x direction.
  - `y_edges     :: Vector{<:Number}`: Absolute position of (the edges of) the bins in the y direction.
  - `n_bins      :: NTuple{2,Int}`: Number of bins in each direction.
  - `bin_size_1D :: Vector{<:Number}`: Widths of the bins in each direction.
  - `bin_size_2D :: Number`: Area of the bins.
  - `origin      :: Vector{<:Number}`: Physical position of the center of the grid.
  - `size        :: Vector{<:Number}`: Physical side length of the grid in each direction (passthrough argument).
  - `n_pixels    :: Int`: Total number of bins.
"""
struct SquareGrid
    x_axis      :: Vector{<:Number}
    y_axis      :: Vector{<:Number}
    x_edges     :: Vector{<:Number}
    y_edges     :: Vector{<:Number}
    n_bins      :: NTuple{2, Int}
    bin_size_1D :: Vector{<:Number}
    bin_size_2D :: Number
    origin      :: Vector{<:Number}
    size        :: Vector{<:Number}
    n_pixels    :: Int

    """
        SquareGrid(
            size   :: Vector{<:Number},
            n_bins :: NTuple{2,Int};
            <keyword arguments>
        )

    Constructor for `SquareGrid`.

    # Arguments

      - `size   :: Vector{<:Number}`: Physical side length of the grid in each direction.
      - `n_bins :: NTuple{2,Int}`: Number of bins in each direction.
      - `origin :: Vector{<:Number}=zeros(runtimeType(grid_length), 3)`: Physical position of the center of the grid.
    """
    function SquareGrid(
        size   :: Vector{<:Number},
        n_bins :: NTuple{2,Int};
        origin :: Vector{<:Number}=zeros(runtimeType(grid_length), 3),
    )

        (
            all(isPositive, size) ||
            throw(ArgumentError("SquareGrid: `size` has to be strictly positive, \
            but I got `size` = $(size)"))
        )

        (
            all(isPositive, n_bins) ||
            throw(ArgumentError("SquareGrid: `n_bins` has to be strictly positive, \
            but I got `n_bins` = $(n_bins)"))
        )

        # Compute the bin dimensions
        bin_size_1D = size ./ n_bins
        bin_size_2D = prod(bin_size_1D)

        # Compute the x and y coordinates of the center of each bin
        center_shift  = 0.5 .* (size .- bin_size_1D)
        x_axis = [(i - 1) * bin_size_1D[1] - center_shift[1] + origin[1] for i in 1:n_bins[1]]
        y_axis = [(i - 1) * bin_size_1D[2] - center_shift[2] + origin[2] for i in 1:n_bins[2]]


        # Compute the x and y coordinates of the edges of each bin
        edge_shift  = 0.5 .* size
        x_edges = [i * bin_size_1D[1] - edge_shift[1] + origin[1] for i in 0:n_bins[1]]
        y_edges = [i * bin_size_1D[2] - edge_shift[2] + origin[2] for i in 0:n_bins[2]]

        n_pixels = prod(n_bins)

        new(
            x_axis,
            y_axis,
            x_edges,
            y_edges,
            n_bins,
            bin_size_1D,
            bin_size_2D,
            origin,
            size,
            n_pixels,
        )

    end

    """
        SquareGrid(
            size   :: Number,
            n_bins :: Int;
            <keyword arguments>
        )

    Constructor for `SquareGrid`.

    # Arguments

      - `size   :: Number`: Physical side length of the grid in every direction.
      - `n_bins :: Int`: Number of bins in every direction.
      - `origin :: Vector{<:Number}=zeros(typeof(size), 3)`: Physical position of the center of the grid.
    """
    function SquareGrid(
        size   :: Number,
        n_bins :: Int;
        origin :: Vector{<:Number}=zeros(typeof(size), 3),
    )

        SquareGrid([size, size], (n_bins, n_bins); origin)

    end
end

"""
Grid with 3 d.o.f.

# Fields

  - `x_axis      :: Vector{<:Number}`: Absolute position of (the center of) the bins in the x direction.
  - `y_axis      :: Vector{<:Number}`: Absolute position of (the center of) the bins in the y direction.
  - `z_axis      :: Vector{<:Number}`: Absolute position of (the center of) the bins in the z direction.
  - `x_edges     :: Vector{<:Number}`: Absolute position of (the edges of) the bins in the x direction.
  - `y_edges     :: Vector{<:Number}`: Absolute position of (the edges of) the bins in the y direction.
  - `z_edges     :: Vector{<:Number}`: Absolute position of (the edges of) the bins in the z direction.
  - `n_bins      :: NTuple{3,Int}`: Number of bins in each direction.
  - `bin_size_1D :: Vector{<:Number}`: Widths of the bins in each direction.
  - `bin_size_2D :: Vector{<:Number}`: Areas of the faces normal to each direction.
  - `bin_size_3D :: Number`: Volume of the bins.
  - `origin      :: Vector{<:Number}`: Physical position of the center of the grid.
  - `size        :: Vector{<:Number}`: Physical side length of the grid in each direction (passthrough argument).
  - `n_voxels    :: Int`: Total number of bins.
"""
struct CubicGrid
    x_axis      :: Vector{<:Number}
    y_axis      :: Vector{<:Number}
    z_axis      :: Vector{<:Number}
    x_edges     :: Vector{<:Number}
    y_edges     :: Vector{<:Number}
    z_edges     :: Vector{<:Number}
    n_bins      :: NTuple{3,Int}
    bin_size_1D :: Vector{<:Number}
    bin_size_2D :: Vector{<:Number}
    bin_size_3D :: Number
    origin      :: Vector{<:Number}
    size        :: Vector{<:Number}
    n_voxels    :: Int

    """
        CubicGrid(
            size   :: Vector{<:Number},
            n_bins :: NTuple{3,Int};
            <keyword arguments>
        )

    Constructor for `CubicGrid`.

    # Arguments

      - `size   :: Vector{<:Number}`: Physical side length of the grid in each direction.
      - `n_bins :: NTuple{3,Int}`: Number of bins in each direction.
      - `origin :: Vector{<:Number}=zeros(runtimeType(grid_length), 3)`: Physical position of the center of the grid.
    """
    function CubicGrid(
        size   :: Vector{<:Number},
        n_bins :: NTuple{3,Int};
        origin :: Vector{<:Number}=zeros(runtimeType(grid_length), 3),
    )

        (
            all(isPositive, size) ||
            throw(ArgumentError("CubicGrid: `size` has to be strictly positive, \
            but I got `size` = $(size)"))
        )

        (
            all(isPositive, n_bins) ||
            throw(ArgumentError("CubicGrid: `n_bins` has to be strictly positive, \
            but I got `n_bins` = $(n_bins)"))
        )

        # Compute the bin dimensions
        bin_size_1D = size ./ n_bins
        bin_size_2D = prod.([bin_size_1D[[2,3]], bin_size_1D[[1,3]], bin_size_1D[[1,2]]])
        bin_size_3D = prod(bin_size_1D)

        # Compute the x and y coordinates of the center of each bin
        center_shift  = 0.5 .* (size .- bin_size_1D)
        x_axis = [(i - 1) * bin_size_1D[1] - center_shift[1] + origin[1] for i in 1:n_bins[1]]
        y_axis = [(i - 1) * bin_size_1D[2] - center_shift[2] + origin[2] for i in 1:n_bins[2]]
        z_axis = [(i - 1) * bin_size_1D[3] - center_shift[3] + origin[3] for i in 1:n_bins[3]]

        # Compute the x and y coordinates of the edges of each bin
        edge_shift  = 0.5 .* size
        x_edges = [i * bin_size_1D[1] - edge_shift[1] + origin[1] for i in 0:n_bins[1]]
        y_edges = [i * bin_size_1D[2] - edge_shift[2] + origin[2] for i in 0:n_bins[2]]
        z_edges = [i * bin_size_1D[3] - edge_shift[3] + origin[3] for i in 0:n_bins[3]]

        n_voxels = prod(n_bins)

        new(
            x_axis,
            y_axis,
            z_axis,
            x_edges,
            y_edges,
            z_edges,
            n_bins,
            bin_size_1D,
            bin_size_2D,
            bin_size_3D,
            origin,
            size,
            n_voxels,
        )

    end

    """
        CubicGrid(
            size   :: Number,
            n_bins :: Int;
            <keyword arguments>
        )

    Constructor for `CubicGrid`.

    # Arguments

      - `size   :: Number`: Physical side length of the grid in every direction.
      - `n_bins :: Int`: Number of bins in every direction.
      - `origin :: Vector{<:Number}=zeros(typeof(size), 3)`: Physical position of the center of the grid.
    """
    function CubicGrid(
        size   :: Number,
        n_bins :: Int;
        origin :: Vector{<:Number}=zeros(typeof(size), 3),
    )

        CubicGrid([size, size, size], (n_bins, n_bins, n_bins); origin)

    end
end

####################################################################################################
# Default theme for Makie.jl
####################################################################################################

"""
Default colors.
"""
const WONG_BLUE    = Makie.wong_colors()[1]
const WONG_ORANGE  = Makie.wong_colors()[2]
const WONG_GREEN   = Makie.wong_colors()[3]
const WONG_PINK    = Makie.wong_colors()[4]
const WONG_CELESTE = Makie.wong_colors()[5]
const WONG_RED     = Makie.wong_colors()[6]
const WONG_YELLOW  = Makie.wong_colors()[7]

"""
Default markers.
"""
const MARKERS = [
    :circle,
    :rect,
    :dtriangle,
    :utriangle,
    :cross,
    :diamond,
    :ltriangle,
    :rtriangle,
    :pentagon,
    :xcross,
    :hexagon,
]

"""
Default line styles.
"""
const LINE_STYLES = [:solid, :dash, :dot, :dashdot, :dashdotdot]

"""
Default cycler.
"""
const CYCLE = Cycle([:color, :linestyle, :marker], covary=true)

"""
Default plot theme.

Regarding the graphic units used, we know that ``1 \\, \\mathrm{mm} = 2.83466 \\, \\mathrm{pt}`` and ``1 \\, \\mathrm{in} = 25.4 \\, \\mathrm{mm}``. Then, if we want ``1 \\, \\mathrm{[code\\,\\,]unit} = 0.1 \\, \\mathrm{mm}`` in vector graphics, we have to use `pt_per_unit` = 0.283466.

For pixel images, we control the ppi with `px_per_unit`. A reasonable high ppi is 600. So, using `px_per_unit` = ``2.3622`` we get ``23.622 \\, \\mathrm{px/mm} \\sim 600 \\, \\mathrm{px/in}`` (remember that ``1 \\, \\mathrm{[code\\,\\,]unit} = 0.1 \\, \\mathrm{mm}``).
"""
const DEFAULT_THEME = Theme(
    #################################################################################################
    # Size of the figures in code units
    # For PDFs and SVGs, 880 [code ]unit = 8.8 cm
    # For PNGs, when printed to a size of 1 point = 0.1 mm, one will get a dpi of 600 (23.622 px/mm)
    #################################################################################################
    size=(880, 880),
    ######################################
    # 35 unit * 0.283466 pt/unit ~ 9.9 pt
    ######################################
    fontsize=35,
    #############################
    # (left, right, bottom, top)
    #############################
    figure_padding=(2, 15, 5, 15),
    palette=(color=Makie.wong_colors(), marker=MARKERS, linestyle=LINE_STYLES),
    CairoMakie=(px_per_unit=2.3622, pt_per_unit=0.283466),
    Axis=(
        xlabelpadding=15,
        xticklabelpad=10,
        xticksize=7,
        xgridvisible=false,
        spinewidth=3,
        xminorticksvisible=true,
        xminorticks=IntervalsBetween(5),
        ylabelpadding=15,
        yticklabelpad=10,
        yticksize=7,
        ygridvisible=false,
        yminorticksvisible=true,
        yminorticks=IntervalsBetween(5),
        ######################################################################################
        # Aspect ratio of the figures. The options are:
        # nothing: The aspect ratio will be chosen by [Makie](https://docs.makie.org/stable/)
        # AxisAspect(n): The aspect ratio will be given by the number `n` = width / height
        # DataAspect(): The aspect ratio of the data will be used
        ######################################################################################
        aspect=AxisAspect(1),
    ),
    Legend=(
        tellheight=false,
        tellwidth=false,
        framevisible=false,
        colgap=20,
        rowgap=10,
        halign=:right,
        valign=:bottom,
        nbanks=3,
        titlegap=-5,
        labelsize=30,
        linewidth=5,
        markersize=28,
        patchsize=(40, 40),
    ),
    Lines=(linewidth=5, cycle=CYCLE),
    VLines=(linewidth=3, cycle=CYCLE),
    HLines=(linewidth=3, cycle=CYCLE),
    ScatterLines=(
        linewidth=5,
        markersize=22,
        cycle=CYCLE,
    ),
    Scatter=(markersize=22, cycle=CYCLE),
    Band=(alpha=0.5, cycle=CYCLE),
    VSpan=(alpha=0.5, cycle=CYCLE),
    Errorbars=(whiskerwidth=10,),
    ########################################################################
    # Alternative colormaps:
    # colormap = :nipy_spectral - nan_color = ColorSchemes.nipy_spectral[1]
    # colormap = :cubehelix     - nan_color = ColorSchemes.cubehelix[1]
    # colormap = :lipari        - nan_color = ColorSchemes.lipari[1]
    ########################################################################
    Heatmap=(colormap=:CMRmap, nan_color=ColorSchemes.CMRmap[1]),
    Colorbar=(
        colormap=:CMRmap,
        size=25,
        ticklabelpad=10,
        minorticksvisible=true,
        ticksize=7,
        labelpadding=2,
    ),
    BarPlot=(
        color_over_background=:black,
        color_over_bar=:black,
        flip_labels_at=10,
        strokecolor=:black,
        strokewidth=1,
        dodge_gap=0.04,
    ),
    Arrows2D=(lengthscale=0.015, color=:white, shaftwidth=2, tipwidth=8),
    Hist=(strokecolor=:black, strokewidth=1),
)

####################################################################################################
# Arepo specific constants simulations
####################################################################################################

"""
Subhalo numbers for the MW and M31 in Hestia simulations.
"""
const HESTIA_SUBHALOS = Dict(
    "Hestia17-11" => Dict(
        :subhalo_number_MW  => 1,
        :subhalo_number_M31 => 0,
    ),
    "Hestia09-18" => Dict(
        :subhalo_number_MW  => 3911,
        :subhalo_number_M31 => 2608,
    ),
    "Hestia37-11" => Dict(
        :subhalo_number_MW  => 920,
        :subhalo_number_M31 => 0,
    ),
)

"""
Cosmological threshold density above which the gas cells/particles can turn into stars.

This value corresponds to `CritOverDensity` ``= 57.7 \\, [\\mathrm{cm^{-3}}]`` in the `param.txt` file (used only in cosmological simulations). Which is converted to internal units within the code using `OverDensThresh` = `CritOverDensity` * `OmegaBaryon` * 3 * `Hubble` * `Hubble` / (8 * `M_PI` * `G`)`. Then, to go to physical units again one has to do: `OverDensThresh`*`UnitDensity_in_cgs`*`cf_a3inv`*`HubbleParam`*`HubbleParam`.

Using the unit factors,

`UnitLength_in_cm`         = ``3.085678 \\times 10^{24}``

`UnitMass_in_g`            = ``1.989 \\times 10^{43}``

`UnitVelocity_in_cm_per_s` = ``100000``

The derived units,

`UnitTime_in_s`      = `UnitLength_in_cm` * `UnitVelocity_in_cm_per_s`^-1 = ``3.08568 \\times 10^{19}``

`UnitDensity_in_cgs` = `UnitMass_in_g` * `UnitLength_in_cm^-3`            = ``6.76991 \\times 10^{-31}``

The parameters,

`OmegaBaryon`       = ``0.048``

`HubbleParam`       = ``0.6777``

`PROTONMASS`        = ``1.67262178 \\times 10^{-24}``

`HYDROGEN_MASSFRAC` = ``0.76``

`GRAVITY`           = ``6.6738 \\times 10^{-8}``

`HUBBLE`            = ``3.2407789 \\times 10^{-18}``

And the derived parameters,

Hubble = `HUBBLE` * `UnitTime_in_s`                                              = ``100``

G      = `GRAVITY` * `UnitLength_in_cm`^-3 * `UnitMass_in_g` * `UnitTime_in_s`^2 = ``43.0187``

One gets,

`OverDensThresh` = 76.8495 [internal units of density]

And, for a cosmological simulation at redshift 0 (`cf_a3inv` = 1), this result in a physical density threshold of ``1.42857 \\times 10^{-5} \\, [\\mathrm{cm^{-3}}]``, or, adding the proton mass, a value of

``\\log_{10} \\rho \\ [\\mathrm{M_\\odot \\, kpc^{-3}}] = 2.548``
"""
const COSMO_THRESHOLD_DENSITY = 353.059u"Msun*kpc^-3"

"""
Threshold density above which the gas cells/particles can turn into stars.

This value corresponds to `CritPhysDensity` ``= 0.318 \\, [\\mathrm{cm^{-3}}]`` in the `param.txt` file (used in cosmological and non-cosmological simulations). Which is converted to internal units within the code using `PhysDensThresh` = `CritPhysDensity` * `PROTONMASS` / `HYDROGEN_MASSFRAC` / `UnitDensity_in_cgs`. Then, to go to physical units again one has to do: `PhysDensThresh` * `UnitDensity_in_cgs` * `cf_a3inv` * `HubbleParam` * `HubbleParam`.

`PhysDensThresh` = ``1.03378 \\times 10^{6}`` [internal units of density]

For a cosmological simulation at redshift 0 (`cf_a3inv` = 1), this result in a physical density threshold of ``0.192 \\, [\\mathrm{cm^{-3}}]``, or, adding the proton mass, a value of

``\\log_{10} \\rho \\, [\\mathrm{M_\\odot \\, kpc^{-3}}] = 6.677``
"""
const THRESHOLD_DENSITY = 4.749326e6u"Msun*kpc^-3"

######################
# Output files paths
######################

"""
Relative path, within the simulation directory, of `sfr.txt`.
"""
const SFR_REL_PATH = "output/sfr.txt"

"""
Relative path, within the simulation directory, of `cpu.txt`.
"""
const CPU_REL_PATH = "output/cpu.txt"

######################
# Cell/particle types
######################

"""
Code index for each type of cell/particle, in C zero-based numbering

# References

See for example Gadget2 [User's Guide](https://wwwmpa.mpa-garching.mpg.de/gadget/users-guide.pdf), or Gadget4 [documentation](https://wwwmpa.mpa-garching.mpg.de/gadget4/).
"""
const LONG_PARTICLE_INDEX = Dict(
    :gas         => 0,
    :dark_matter => 1,
    :disk        => 2,
    :bulge       => 3,
    :stellar     => 4,
    :black_hole  => 5,
    :tracer      => 6,
)

"""
Human readable name for each type of cell/particle.
"""
const LONG_PARTICLE_NAMES = Dict(
    :gas         => "Gas cells",
    :dark_matter => "HR DM particles",
    :disk        => "IR DM particles",
    :bulge       => "LR DM particles",
    :stellar     => "Stellar particles",
    :black_hole  => "Black hole particles",
    :tracer      => "Tracer particles",
)

"""
Human readable name for each type of cell/particle.
"""
const ISOLATED_PARTICLE_NAMES = Dict(
    :gas          => "Gas cells",
    :dark_matter  => "DM particles",
    :disk         => "Stellar disk",
    :bulge        => "Stellar bulge",
)

"""
Human readable name for each morphological component.
"""
const MORPHOLOGICAL_COMPONENTS = Dict(
    :disk       => "Disk",
    :bulge      => "Bulge",
    :thin_disk  => "Thin disk",
    :thick_disk => "Thick disk",
)

"""
Current cell/particle index in use.
"""
const PARTICLE_INDEX = LONG_PARTICLE_INDEX

"""
Current human readable name of each cell/particle type in use.
"""
const PARTICLE_NAMES = LONG_PARTICLE_NAMES

"""
Type of cell/particle corresponding to each code index.
"""
const INDEX_PARTICLE = Dict(n => symbol for (symbol, n) in PARTICLE_INDEX)

"""
Type of cell/particle corresponding to each internal code name (data group in the HDF5 output).
"""
const PARTICLE_TYPE = Dict("PartType$n" => symbol for (symbol, n) in PARTICLE_INDEX)

"""
Internal code name (data group in the HDF5 output) corresponding to each type of cell/particle.
"""
const PARTICLE_CODE_NAME = Dict(symbol => "PartType$n" for (symbol, n) in PARTICLE_INDEX)

##############################
# Default filter dictionaries
##############################

"""
Filter dictionary that does not exclude any cell/particle.
"""
const PASS_ALL = Dict{Symbol,IndexType}(key => (:) for key in keys(PARTICLE_INDEX))

"""
Filter dictionary that excludes every cell/particle.
"""
const PASS_NONE = Dict{Symbol,IndexType}(key => Int[] for key in keys(PARTICLE_INDEX))

##########################
# Transformations presets
##########################

"""
List of symbols for the transformation presets.

The options have the form :{component}_{group} selecting which cells/particle to consider for the center of mass (new origin for the translation) and principal axis (new reference system for the rotation), where component can be:

      + :all         -> Every component present in data_dict
      + :{component} -> One of the keys of PARTICLE_INDEX

    and group can be:
      + :box     -> Whole simulation box
      + :halo    -> Main halo
      + :subhalo -> Main subhalo
"""
const TRANSFORM_LIST = [
    Symbol(component, :_, group)
    for component in push!(collect(keys(PARTICLE_INDEX)), :all)
    for group in [:box, :halo, :subhalo]
]

"""
Dictionary mapping each transformation to its component and group.
"""
const TRANSFORM_SPLITS = Dict(
    Symbol(component, :_, group) => (component, group)
    for component in push!(collect(keys(PARTICLE_INDEX)), :all)
    for group in [:box, :halo, :subhalo]
)

###################
# Tracked elements
###################

"""
Code index for each tracked element.
"""
const ELEMENT_INDEX = Dict(
    :H     => 1,  # Hydrogen
    :He    => 2,  # Helium
    :C     => 3,  # Carbon
    :N     => 4,  # Nitrogen
    :O     => 5,  # Oxygen
    :Ne    => 6,  # Neon
    :Mg    => 7,  # Magnesium
    :Si    => 8,  # Silicon
    :Fe    => 9,  # Iron
)

"""
List of element above helium.
"""
const METAL_LIST = [:C, :N, :O, :Ne, :Mg, :Si, :Fe]

"""
List of element indices above helium.
"""
const METAL_LIST_IDX = [get(ELEMENT_INDEX, metal, 0) for metal in METAL_LIST]

#################################
# Star formation model constants
#################################

"""
ODE index corresponding to each gas phase in our star formation model.

  * Ionized gas fraction:   fi(t) = Mi(t) / MC --> data_dict[:gas]["FRAC"][1, :]
  * Atomic gas fraction:    fa(t) = Ma(t) / MC --> data_dict[:gas]["FRAC"][2, :]
  * Molecular gas fraction: fm(t) = Mm(t) / MC --> data_dict[:gas]["FRAC"][3, :]
  * Stellar fraction:       fs(t) = Ms(t) / MC --> data_dict[:gas]["FRAC"][4, :]
  * Metal fraction:         fZ(t) = MZ(t) / MC --> data_dict[:gas]["FRAC"][5, :]
  * Dust fraction:          fd(t) = Md(t) / MC --> data_dict[:gas]["FRAC"][6, :]
"""
const SFM_IDX = Dict(
    :ode_ionized   => 1,
    :ode_atomic    => 2,
    :ode_molecular => 3,
    :ode_stellar   => 4,
    :ode_metals    => 5,
    :ode_dust      => 6,
)

"""
Star formation model parameters
"""
const AREPO_SFM = (
    Cxd   = 0.2836,                         # Constant for the initial condition of dust and metals.
    εff   = 1.0,                            # Star formation efficiency per free-fall time
    αH    = 2.6e-13u"cm^3 * s^-1",          # Case B recombination coefficient at T = 10^4 K
    Rd    = 3.5e-17u"cm^3 * s^-1",          # H2 formation rate on dust grains at T = 100 K
    Cdg   = 2.4634u"mp^-1 * Myr^-1 * cm^3", # Dust growth coefficient
    σion  = 6.3e-18u"cm^2",                 # Hydrogen ionization cross section at the Lyman limit
    σdiss = 2.1e-19u"cm^2",                 # H2 dissociation cross section in the Lyman-Werner band
    σd    = 4.0e-21u"cm^2",                 # Effective dust cross-section per hydrogen atom
    ωH2   = 0.2,                            # Constant for the H2 self-shielding factor
    xnf   = 5.0e14u"cm^-2",                 # H2 self-shielding column density factor
    Cswm  = 2.233e3u"Msun",                 # Mass swept-up by a SNe shock times the efficiency of dust destruction
)

####################################################################################################
# Catalog of quantities in the snapshots and group catalog files
####################################################################################################

"""
Dimensional properties of the quantities in the snapshots and group catalog files.
"""
const QUANTITIES = Dict(

    ######################
    # Snapshot quantities
    ######################

    "CLKT" => SnapshotBlock("", Unitful.𝐓, :internal),
    "GAGE" => SnapshotBlock("GFM_StellarFormationTime", Unitful.𝐓, :internal),
    "GME2" => SnapshotBlock("GFM_Metals", Unitful.NoDims, Unitful.NoUnits),
    "GMET" => SnapshotBlock("GFM_Metals", Unitful.NoDims, Unitful.NoUnits),
    "GZ  " => SnapshotBlock("GFM_Metallicity", Unitful.NoDims, Unitful.NoUnits),
    "GZ2 " => SnapshotBlock("GFM_Metallicity", Unitful.NoDims, Unitful.NoUnits),
    "ID  " => SnapshotBlock("ParticleIDs", Unitful.NoDims, Unitful.NoUnits),
    "PAID" => SnapshotBlock("ParentID", Unitful.NoDims, Unitful.NoUnits),
    "TRID" => SnapshotBlock("TracerID", Unitful.NoDims, Unitful.NoUnits),
    "MASS" => SnapshotBlock("Masses", Unitful.𝐌, :internal),
    "NE  " => SnapshotBlock("ElectronAbundance", Unitful.NoDims, Unitful.NoUnits),
    "NH  " => SnapshotBlock("NeutralHydrogenAbundance", Unitful.NoDims, Unitful.NoUnits),
    "NHP " => SnapshotBlock("IonizedHydrogenAbundance", Unitful.NoDims, Unitful.NoUnits),
    "POS " => SnapshotBlock("Coordinates", Unitful.𝐋, :internal),
    "PRES" => SnapshotBlock("Pressure", Unitful.𝐌 * Unitful.𝐋^-1 * Unitful.𝐓^-2, :internal),
    "RHO " => SnapshotBlock("Density", Unitful.𝐌 * Unitful.𝐋^-3, :internal),
    "SFR " => SnapshotBlock("StarFormationRate", Unitful.𝐌 * Unitful.𝐓^-1, u"Msun * yr^-1"),
    "SOFT" => SnapshotBlock("Softenings", Unitful.𝐋, :internal),
    "TEMP" => SnapshotBlock("Temperature", Unitful.𝚯, u"K"),
    "U   " => SnapshotBlock("InternalEnergy", Unitful.𝐋^2 * Unitful.𝐓^-2, :internal),
    "VEL " => SnapshotBlock("Velocities", Unitful.𝐋 * Unitful.𝐓^-1, :internal),
    "TSTP" => SnapshotBlock("TimeStep", Unitful.𝐓, :internal),
    "POT " => SnapshotBlock("Potential", Unitful.𝐋^2 * Unitful.𝐓^-2, :pot),
    "MACH" => SnapshotBlock("Machnumber", Unitful.NoDims, Unitful.NoUnits),
    "GDCM" => SnapshotBlock("GFM_DustCapMass", Unitful.𝐌, :internal),
    "TAGR" => SnapshotBlock("GFM_DustTauGrowth", Unitful.𝐓, u"Gyr"),
    "GDZ " => SnapshotBlock("GFM_DustMetallicity", Unitful.NoDims, Unitful.NoUnits),
    "GDZ2" => SnapshotBlock("GFM_DustMetallicity", Unitful.NoDims, Unitful.NoUnits),
    "GDAG" => SnapshotBlock("GFM_DustAGB", Unitful.NoDims, Unitful.NoUnits),
    "GDA2" => SnapshotBlock("GFM_DustAGB", Unitful.NoDims, Unitful.NoUnits),
    "GDII" => SnapshotBlock("GFM_DustSNII", Unitful.NoDims, Unitful.NoUnits),
    "GDI2" => SnapshotBlock("GFM_DustSNII", Unitful.NoDims, Unitful.NoUnits),
    "GDIa" => SnapshotBlock("GFM_DustSNIa", Unitful.NoDims, Unitful.NoUnits),
    "GDJ2" => SnapshotBlock("GFM_DustSNIa", Unitful.NoDims, Unitful.NoUnits),

    #####################
    # sfr.txt quantities
    #####################

    # Time or scale factor
    "SFC1" => SnapshotBlock("", Unitful.𝐓, :internal),
    # Total stellar mass to be formed prior to stochastic sampling
    "SFC2" => SnapshotBlock("", Unitful.𝐌, :internal),
    # Instantaneous star formation rate of all cells
    "SFC3" => SnapshotBlock("", Unitful.𝐌 * Unitful.𝐓^-1, u"Msun * yr^-1"),
    # Instantaneous star formation rate of active cells
    "SFC4" => SnapshotBlock("", Unitful.𝐌 * Unitful.𝐓^-1, u"Msun * yr^-1"),
    # Total mass in stars formed after stochastic sampling
    "SFC5" => SnapshotBlock("", Unitful.𝐌, :internal),
    # Cumulative stellar mass formed
    "SFC6" => SnapshotBlock("", Unitful.𝐌, :internal),

    ####################
    # EL_SFR quantities
    ####################

    # Integration time, for gas cells and stellar particles
    "ODIT" => SnapshotBlock("ODE_IntegrationTime", Unitful.𝐓, u"Myr"),
    # Scale factor, for gas cells and stellar particles
    "PARA" => SnapshotBlock("ODE_ParameterA", Unitful.NoDims, Unitful.NoUnits),
    # SNII fraction, for gas cells and stellar particles
    "PARS" => SnapshotBlock("ODE_ParameterSNII", Unitful.𝐌^-1, u"Msun^-1"),
    # Redshift parameter, for gas cells and stellar particles
    "PARz" => SnapshotBlock("ODE_Parameterz", Unitful.NoDims, Unitful.NoUnits),
    # UVB photoionization rate, for gas cells and stellar particles
    "PARU" => SnapshotBlock("ODE_ParameterUVB", Unitful.𝐓^-1, u"Myr^-1"),
    # LWB photodissociation rate, for gas cells and stellar particles
    "PARL" => SnapshotBlock("ODE_ParameterLWB", Unitful.𝐓^-1, u"Myr^-1"),
    # Star formation time parameter, for gas cells and stellar particles
    "TAUS" => SnapshotBlock("ODE_TauS", Unitful.𝐓, u"Myr"),
    # Gas density, for gas cells and stellar particles
    "RHOC" => SnapshotBlock("ODE_ParameterCellDensity", Unitful.𝐋^-3, u"cm^-3"),
    # Gas metallicity, for gas cells and stellar particles
    "PARZ" => SnapshotBlock("ODE_ParameterMetallicity", Unitful.NoDims, Unitful.NoUnits),
    # Column height, for gas cells and stellar particles
    "PARH" => SnapshotBlock("ODE_ParameterColumnHeight", Unitful.𝐋, u"cm"),
    # Photodissociation parameter, for gas cells and stellar particles
    "ETAD" => SnapshotBlock("ODE_ParameterEtaD", Unitful.NoDims, Unitful.NoUnits),
    # Photoionization parameter, for gas cells and stellar particles
    "ETAI" => SnapshotBlock("ODE_ParameterEtaI", Unitful.NoDims, Unitful.NoUnits),
    # Mass recycling parameter, for gas cells and stellar particles
    "PARR" => SnapshotBlock("ODE_ParameterR", Unitful.NoDims, Unitful.NoUnits),
    # Metallicity of the supernova ejecta, for gas cells and stellar particles
    "PAZN" => SnapshotBlock("ODE_ParameterZsn", Unitful.NoDims, Unitful.NoUnits),
    # Gas fractions, for gas cells and stellar particles
    "FRAC" => SnapshotBlock("ODE_Fractions", Unitful.NoDims, Unitful.NoUnits),
    # Star formation flag, for gas cells
    "SFFL" => SnapshotBlock("ODE_SfFlag", Unitful.NoDims, Unitful.NoUnits),
    # Cold gas fraction, for gas cells and stellar particles
    "COLF" => SnapshotBlock("ODE_ColdMassFrac", Unitful.NoDims, Unitful.NoUnits),
    # Parent gas mass (at the moment of star formation), for stellar particles
    "GMAS" => SnapshotBlock("ODE_GasMass", Unitful.𝐌, :internal),
    # Parent SFR (at the moment of star formation), for stellar particles
    "GSFR" => SnapshotBlock("ODE_GasSFR", Unitful.𝐌 * Unitful.𝐓^-1, u"Msun * yr^-1"),
    # Parent gas pressure (at the moment of star formation), for stellar particles
    "GPRE" => SnapshotBlock("ODE_GasPressure", Unitful.𝐌 * Unitful.𝐋^-1 * Unitful.𝐓^-2, :internal),
    # Parent position (at the moment of star formation), for stellar particles
    "GPOS" => SnapshotBlock("ODE_GasPosition", Unitful.𝐋, :internal),
    # Parent velocity (at the moment of star formation), for stellar particles
    "GVEL" => SnapshotBlock("ODE_GasVelocity", Unitful.𝐋 * Unitful.𝐓^-1, :internal),

    ######################
    # YA_SFR quantities
    ######################

    # Isothermal sound speed [cm * s^(-1)]
    "YACC" => SnapshotBlock("YA_cs_cell", Unitful.𝐋 * Unitful.𝐓^-1, u"cm * s^-1"),
    # Gas temperature [K]
    "YATC" => SnapshotBlock("YA_T_cell", Unitful.𝚯, u"K"),
    # Length scale of the gas cell [cm]
    "YALC" => SnapshotBlock("YA_l_cell", Unitful.𝐋, u"cm"),
    # Specific energy (turbulent) dissipation rate [cm^2 * s^(-3)]
    "YAEP" => SnapshotBlock("YA_epsilon", Unitful.𝐋^2 * Unitful.𝐓^-3, u"cm^2 * s^-3"),
    # Mach number [dimensionless]
    "YAMS" => SnapshotBlock("YA_Mach_cell", Unitful.NoDims, Unitful.NoUnits),
    # Magnetic field [g^(1/2) * cm^(-1/2) * s^(-1)]
    "YABC" => SnapshotBlock(
        "YA_B_cell",
        Unitful.𝐌^(1/2) * Unitful.𝐋^(-1/2) * Unitful.𝐓^-1,
        u"g^(1/2) * cm^(-1/2) * s^-1",
    ),
    # Alfvén velocity [cm * s^(-1)]
    "YAVA" => SnapshotBlock("YA_vA", Unitful.𝐋 * Unitful.𝐓^-1, u"cm * s^-1"),
    # Velocity dispersion [cm * s^(-1)]
    "YASV" => SnapshotBlock("YA_sigma_v", Unitful.𝐋 * Unitful.𝐓^-1, u"cm * s^-1"),
    # Alfvén Mach [dimensionless]
    "YAMA" => SnapshotBlock("YA_MA", Unitful.NoDims, Unitful.NoUnits),
    # Magnetic pressure [g * cm^(-1) * s^(-2)]
    "YAPM" => SnapshotBlock(
        "YA_P_mag",
        Unitful.𝐌 * Unitful.𝐋^-1 * Unitful.𝐓^-2,
        u"g * cm^-1 * s^-2",
    ),
    # Thermal-to-magnetic pressure ratio [dimensionless]
    "YAB " => SnapshotBlock("YA_beta", Unitful.NoDims, Unitful.NoUnits),
    # Standard deviation for the logarithmic density PDF [dimensionless]
    "YASS" => SnapshotBlock("YA_sigma_s", Unitful.NoDims, Unitful.NoUnits),
    # Mean logarithmic density [dimensionless]
    "YAS0" => SnapshotBlock("YA_s0", Unitful.NoDims, Unitful.NoUnits),
    # Effective sound speed [cm * s^(-1)]
    "YACE" => SnapshotBlock("YA_cs_eff", Unitful.𝐋 * Unitful.𝐓^-1, u"cm * s^-1"),
    # Bonnor-Ebert mass for the mean density (s = 0) [M☉]
    "YAMB" => SnapshotBlock("YA_MBE0", Unitful.𝐌, u"Msun"),
    # First root of Mex'(s) [dimensionless]
    "YASR" => SnapshotBlock("YA_sr_m", Unitful.NoDims, Unitful.NoUnits),
    # Maximum excess mass [M☉]
    "YAMF" => SnapshotBlock("YA_Mfe", Unitful.𝐌, u"Msun"),
    # Critical logarithmic density [dimensionless]
    "YASC" => SnapshotBlock("YA_sc", Unitful.NoDims, Unitful.NoUnits),
    # Free fall time [yr]
    "YATF" => SnapshotBlock("YA_tff", Unitful.𝐓, u"yr"),
    # status = 0 -> Success
    # status = 1 -> Failed to compute the Mach number
    # status = 2 -> We are in the sub-Alfvénic regime (MA < 2)
    # status = 3 -> There are no real roots. The cell is stable
    # status = 4 -> No roots exist for Mex(s). The cell is stable
    # status = 5 -> Could not find the root. The root is not within [s_left, sr_m]
    # status = 6 -> Could not find the root. Max iterations reached with no convergence
    "YAS " => SnapshotBlock("YA_status", Unitful.NoDims, Unitful.NoUnits),

    ##############################
    # Halo (FoF group) quantities
    ##############################

    "G_CM"        => SnapshotBlock("GroupCM", Unitful.𝐋, :internal),
    "G_LenType"   => SnapshotBlock("GroupLenType", Unitful.NoDims, Unitful.NoUnits),
    "G_Mass"      => SnapshotBlock("GroupMass", Unitful.𝐌, :internal),
    "G_MassType"  => SnapshotBlock("GroupMassType", Unitful.𝐌, :internal),
    "G_Nsubs"     => SnapshotBlock("GroupNsubs", Unitful.NoDims, Unitful.NoUnits),
    "G_Pos"       => SnapshotBlock("GroupPos", Unitful.𝐋, :internal),
    "G_M_Crit200" => SnapshotBlock("Group_M_Crit200", Unitful.𝐌, :internal),
    "G_R_Crit200" => SnapshotBlock("Group_R_Crit200", Unitful.𝐋, :internal),
    "G_Vel"       => SnapshotBlock("GroupVel", Unitful.𝐋 * Unitful.𝐓^-1, :gvel),

    ###############################
    # Subhalo (subfind) quantities
    ###############################

    "S_CM"          => SnapshotBlock("SubhaloCM", Unitful.𝐋, :internal),
    "S_HalfmassRad" => SnapshotBlock("SubhaloHalfmassRad", Unitful.𝐋, :internal),
    "S_LenType"     => SnapshotBlock("SubhaloLenType", Unitful.NoDims, Unitful.NoUnits),
    "S_Mass"        => SnapshotBlock("SubhaloMass", Unitful.𝐌, :internal),
    "S_MassType"    => SnapshotBlock("SubhaloMassType", Unitful.𝐌, :internal),
    "S_Pos"         => SnapshotBlock("SubhaloPos", Unitful.𝐋, :internal),
    "S_Vel"         => SnapshotBlock("SubhaloVel", Unitful.𝐋 * Unitful.𝐓^-1, u"km * s^-1"),
)

####################################################################################################
# Catalog of derived quantities
####################################################################################################

"""
List of physical quantities and their properties.
"""
QTY_REGISTRY = Dict{Symbol,BaseQuantity}(

    #################
    # Gas quantities
    #################

    :temperature => BaseQuantity(;
        id             = :temperature,
        qty_label      = L"T",
        request        = Dict(:gas => ["TEMP"]),
        unit           = u"K",
        exp_factor     = 0,
        cp_type        = :gas,
        scatter_func   = dd -> dd[:gas]["TEMP"],
        integrate_func = dd -> nothing,
    ),

    :pressure => BaseQuantity(;
        id             = :pressure,
        qty_label      = L"P",
        request        = Dict(:gas => ["PRES"]),
        unit           = u"Pa",
        exp_factor     = 0,
        cp_type        = :gas,
        scatter_func   = dd -> dd[:gas]["PRES"],
        integrate_func = dd -> nothing,
    ),

    #################
    # SFR quantities
    #################

    :sfr => BaseQuantity(;
        id             = :sfr,
        qty_label      = L"\text{SFR}",
        request        = Dict(:stellar => ["MASS", "GAGE"]),
        unit           = u"Msun * yr^-1",
        exp_factor     = 0,
        cp_type        = :stellar,
        scatter_func   = dd -> computeSnapshotSFR(dd),
        integrate_func = dd -> sum(computeSnapshotSFR(dd); init=0.0u"Msun * yr^-1"),
    ),

    :ssfr => BaseQuantity(;
        id             = :ssfr,
        qty_label      = L"\text{sSFR}",
        request        = Dict(:stellar => ["GAGE", "MASS"]),
        unit           = u"Gyr^-1",
        exp_factor     = 0,
        cp_type        = :stellar,
        scatter_func   = dd -> computeSnapshotSSFR(dd),
        integrate_func = dd -> let

            sfr  = sum(computeSnapshotSFR(dd); init=0.0u"Msun * yr^-1")
            mass = sum(dd[:stellar]["MASS"]; init=0.0u"Msun")

            if iszero(mass)
                0.0u"Gyr^-1"
            else
                sfr / mass
            end

        end,
    ),

    :observational_sfr => BaseQuantity(;
        id             = :observational_sfr,
        qty_label      = L"\text{SFR}",
        request        = Dict(:stellar => ["MASS", "GAGE"]),
        unit           = u"Msun * yr^-1",
        exp_factor     = 0,
        cp_type        = :stellar,
        scatter_func   = dd -> computeSFR(dd),
        integrate_func = dd -> sum(computeSFR(dd); init=0.0u"Msun * yr^-1"),
    ),

    :observational_ssfr => BaseQuantity(;
        id             = :observational_ssfr,
        qty_label      = L"\text{sSFR}",
        request        = Dict(:stellar => ["MASS", "GAGE"]),
        unit           = u"Gyr^-1",
        exp_factor     = 0,
        cp_type        = :stellar,
        scatter_func   = dd -> computeSSFR(dd),
        integrate_func = dd -> let

            sfr  = sum(computeSFR(dd); init=0.0u"Msun * yr^-1")
            mass = sum(dd[:stellar]["MASS"]; init=0.0u"Msun")

            if iszero(mass)
                NaN * u"Gyr^-1"
            else
                sfr / mass
            end

        end,
    ),

    :sfr_area_density => BaseQuantity(;
        id             = :sfr_area_density,
        qty_label      = L"\Sigma_\text{SFR}",
        request        = Dict(:stellar => ["MASS", "GAGE", "POS "]),
        unit           = u"Msun * yr^-1 * kpc^-2",
        exp_factor     = 0,
        cp_type        = :stellar,
        scatter_func   = dd -> nothing,
        integrate_func = dd -> nothing,
    ),

    :sfr_density => BaseQuantity(;
        id             = :sfr_density,
        qty_label      = L"\rho_\text{SFR}",
        request        = Dict(:stellar => ["MASS", "GAGE", "POS "]),
        unit           = u"Msun * yr^-1 * kpc^-3",
        exp_factor     = 0,
        cp_type        = :stellar,
        scatter_func   = dd -> nothing,
        integrate_func = dd -> nothing,
    ),

    :stellar_age => BaseQuantity(;
        id             = :stellar_age,
        qty_label      = L"\text{Stellar age}",
        request        = Dict(:stellar => ["GAGE"]),
        unit           = u"Gyr",
        exp_factor     = 0,
        cp_type        = :stellar,
        scatter_func   = dd -> computeStellarAge(dd),
        integrate_func = dd -> nothing,
    ),

    :stellar_birth_time => BaseQuantity(;
        id             = :stellar_birth_time,
        qty_label      = L"\text{Stellar birth time}",
        request        = Dict(:stellar => ["GAGE"]),
        unit           = u"Gyr",
        exp_factor     = 0,
        cp_type        = :stellar,
        scatter_func   = dd -> computeStellarBirthTime(dd),
        integrate_func = dd -> nothing,
    ),

    :gas_sfr => BaseQuantity(;
        id             = :gas_sfr,
        qty_label      = L"\mathrm{SFR_{gas}}",
        request        = Dict(:gas => ["SFR "]),
        unit           = u"Msun * yr^-1",
        exp_factor     = 0,
        cp_type        = :gas,
        scatter_func   = dd -> dd[:gas]["SFR "],
        integrate_func = dd -> nothing,
    ),

    :gas_sfr_area_density => BaseQuantity(;
        id             = :gas_sfr_area_density,
        qty_label      = L"\Sigma_\text{SFR}^\text{gas}",
        request        = Dict(:gas => ["SFR ", "POS "]),
        unit           = u"Msun * yr^-1 * kpc^-2",
        exp_factor     = 0,
        cp_type        = :gas,
        scatter_func   = dd -> nothing,
        integrate_func = dd -> nothing,
    ),

    :mu_mol => BaseQuantity(;
        id             = :mu_mol,
        qty_label      = L"\mu_\mathrm{mol}",
        request        = Dict(
            :stellar => ["MASS", "POS "],
            :gas     => ["MASS", "FRAC", "POS ", "RHO "],
        ),
        unit           = Unitful.NoUnits,
        exp_factor     = 0,
        cp_type        = nothing,
        scatter_func   = dd -> nothing,
        integrate_func = dd -> nothing,
    ),

    ##################
    # Time quantities
    ##################

    :scale_factor => BaseQuantity(;
        id             = :scale_factor,
        qty_label      = L"a",
        request        = Dict{Symbol,Vector{String}}(),
        unit           = Unitful.NoUnits,
        exp_factor     = 0,
        cp_type        = nothing,
        scatter_func   = dd -> nothing,
        integrate_func = dd -> dd[:sim_data].simulation_table[dd[:snap_data].global_index, :scale_factors],
    ),

    :redshift => BaseQuantity(;
        id             = :redshift,
        qty_label      = L"z",
        request        = Dict{Symbol,Vector{String}}(),
        unit           = Unitful.NoUnits,
        exp_factor     = 0,
        cp_type        = nothing,
        scatter_func   = dd -> nothing,
        integrate_func = dd -> dd[:sim_data].simulation_table[dd[:snap_data].global_index, :redshifts],
    ),

    :physical_time => BaseQuantity(;
        id             = :physical_time,
        qty_label      = L"t",
        request        = Dict{Symbol,Vector{String}}(),
        unit           = u"Gyr",
        exp_factor     = 0,
        cp_type        = nothing,
        scatter_func   = dd -> nothing,
        integrate_func = dd -> dd[:sim_data].simulation_table[dd[:snap_data].global_index, :physical_times],
    ),

    :lookback_time => BaseQuantity(;
        id             = :lookback_time,
        qty_label      = L"\mathrm{Lookback \,\, time}",
        request        = Dict{Symbol,Vector{String}}(),
        unit           = u"Gyr",
        exp_factor     = 0,
        cp_type        = nothing,
        scatter_func   = dd -> nothing,
        integrate_func = dd -> dd[:sim_data].simulation_table[dd[:snap_data].global_index, :lookback_times],
    ),

    :time_step => BaseQuantity(;
        id             = :time_step,
        qty_label      = L"\mathrm{Number \,\, of \,\, time \,\, steps}",
        request        = Dict{Symbol,Vector{String}}(),
        unit           = Unitful.NoUnits,
        exp_factor     = 0,
        cp_type        = nothing,
        scatter_func   = dd -> nothing,
        integrate_func = dd -> nothing,
    ),

    :clock_time_s => BaseQuantity(;
        id             = :clock_time_s,
        qty_label      = L"\mathrm{Wallclock \,\, time}",
        request        = Dict{Symbol,Vector{String}}(),
        unit           = u"s",
        exp_factor     = 0,
        cp_type        = nothing,
        scatter_func   = dd -> nothing,
        integrate_func = dd -> nothing,
    ),

    :clock_time_percent => BaseQuantity(;
        id             = :clock_time_percent,
        qty_label      = L"\mathrm{Wallclock \,\, time \, [\%]}",
        request        = Dict{Symbol,Vector{String}}(),
        unit           = Unitful.NoUnits,
        exp_factor     = 0,
        cp_type        = nothing,
        scatter_func   = dd -> nothing,
        integrate_func = dd -> nothing,
    ),

    :tot_clock_time_s => BaseQuantity(;
        id             = :tot_clock_time_s,
        qty_label      = L"\mathrm{Cumulative \,\, wallclock \,\, time}",
        request        = Dict{Symbol,Vector{String}}(),
        unit           = u"s",
        exp_factor     = 0,
        cp_type        = nothing,
        scatter_func   = dd -> nothing,
        integrate_func = dd -> nothing,
    ),

    :tot_clock_time_percent => BaseQuantity(;
        id             = :tot_clock_time_percent,
        qty_label      = L"\mathrm{Cumulative \,\, wallclock \,\, time \, [\%]}",
        request        = Dict{Symbol,Vector{String}}(),
        unit           = Unitful.NoUnits,
        exp_factor     = 0,
        cp_type        = nothing,
        scatter_func   = dd -> nothing,
        integrate_func = dd -> nothing,
    ),

    # See computeVirialAccretion() and computeDiskAccretion() in ./src/analysis/compute_quantities/masses.jl
    :mass_accretion => BaseQuantity(;
        id             = :mass_accretion,
        qty_label      = L"\dot{M}_\text{acc}",
        request        = Dict(
            :gas         => ["ID  ", "POS ", "MASS"],
            :stellar     => ["ID  ", "POS ", "MASS"],
            :black_hole  => ["ID  ", "POS ", "MASS"],
            :dark_matter => ["ID  ", "POS ", "MASS"],
            :group       => ["G_R_Crit200", "G_Nsubs", "G_Pos"],
            :subhalo     => ["S_Pos"],
            :tracer      => ["PAID", "TRID"],
        ),
        unit           = u"Msun * yr^-1",
        exp_factor     = 0,
        cp_type        = :nothing,
        scatter_func   = dd -> nothing,
        integrate_func = dd -> nothing,
    ),

    #########################
    # Metallicity quantities
    #########################

    :gas_metallicity => BaseQuantity(;
        id             = :gas_metallicity,
        qty_label      = L"Z_\mathrm{gas} \, [\mathrm{Z_\odot}]",
        request        = Dict(:gas => ["GZ  ", "MASS"]),
        unit           = Unitful.NoUnits,
        exp_factor     = 0,
        cp_type        = :gas,
        scatter_func   = dd -> computeFraction(dd, :Z_gas) ./ SOLAR_METALLICITY[],
        integrate_func = dd -> let

            Mz = sum(computeMass(dd, :Z_gas); init=0.0u"Msun")
            Mg = sum(computeMass(dd, :gas); init=0.0u"Msun")

            if iszero(Mg)
                NaN
            else
                (Mz / Mg) / SOLAR_METALLICITY[]
            end

        end,
    ),

    :stellar_metallicity => BaseQuantity(;
        id             = :stellar_metallicity,
        qty_label      = L"Z_\star \, [\mathrm{Z_\odot}]",
        request        = Dict(:stellar => ["GZ2 ", "MASS"]),
        unit           = Unitful.NoUnits,
        exp_factor     = 0,
        cp_type        = :stellar,
        scatter_func   = dd -> computeFraction(dd, :Z_stellar) ./ SOLAR_METALLICITY[],
        integrate_func = dd -> let

            Mz = sum(computeMass(dd, :Z_stellar); init=0.0u"Msun")
            Ms = sum(computeMass(dd, :stellar); init=0.0u"Msun")

            if iszero(Ms)
                NaN
            else
                (Mz / Ms) / SOLAR_METALLICITY[]
            end

        end,
    ),

    :ode_metallicity => BaseQuantity(;
        id             = :ode_metallicity,
        qty_label      = L"Z_\text{ODE} \, [\mathrm{Z_\odot}]",
        request        = Dict(:gas => ["MASS", "NH  ", "NHP ", "FRAC", "RHO ", "GZ  "]),
        unit           = Unitful.NoUnits,
        exp_factor     = 0,
        cp_type        = :gas,
        scatter_func   = dd -> computeFraction(dd, :ode_metals) ./ SOLAR_METALLICITY[],
        integrate_func = dd -> let

            Mz = sum(computeMass(dd, :ode_metals); init=0.0u"Msun")
            Mg = sum(computeMass(dd, :gas); init=0.0u"Msun")

            if iszero(Mg)
                NaN
            else
                (Mz / Mg) / SOLAR_METALLICITY[]
            end

        end,
    ),

)

"""
List of  physical components of a simulation.
"""
const COMPONENTS = Dict(

    ############################
    # Particle-based components
    ############################

    # Stellar particles
    :stellar     => (qty_label="\\star", cp_type=:stellar),
    # Dark matter particles
    :dark_matter => (qty_label="\\text{DM}", cp_type=:dark_matter),
    # Black hole particles
    :black_hole  => (qty_label="\\text{BH}", cp_type=:black_hole),
    # Metals in the stars
    :Z_stellar   => (qty_label="\\text{Z\\!\\star}", cp_type=:stellar),

    #######################
    # Gas-based components
    #######################

    # Total gas
    :gas          => (qty_label="\\text{gas}", cp_type=:gas),
    # Hydrogen
    :hydrogen     => (qty_label="\\text{H}", cp_type=:gas),
    # Helium
    :helium       => (qty_label="\\text{He}", cp_type=:gas),
    # Metals in the gas
    :Z_gas        => (qty_label="\\text{Z\\,gas}", cp_type=:gas),
    # Ionized gas (using the Arepo data)
    :ionized      => (qty_label="\\text{HII}", cp_type=:gas),
    # Neutral gas (using the Arepo data)
    :neutral      => (qty_label="\\mathrm{HI + H2}", cp_type=:gas),
    # Atomic gas (using the Blitz et al. (2006) relation)
    :br_atomic    => (qty_label="\\text{BR HI}", cp_type=:gas),
    # Molecular gas (using the Blitz et al. (2006) relation)
    :br_molecular => (qty_label="\\text{BR H2}", cp_type=:gas),

    ###############################
    # Components from our SF model
    ###############################

    # Ionized gas
    :ode_ionized           => (qty_label="\\text{ODE i}", cp_type=:gas),
    # Atomic gas
    :ode_atomic            => (qty_label="\\text{ODE a}", cp_type=:gas),
    # Molecular gas
    :ode_molecular         => (qty_label="\\text{ODE m}", cp_type=:gas),
    # Stars
    :ode_stellar           => (qty_label="\\text{ODE s}", cp_type=:gas),
    # Metals
    :ode_metals            => (qty_label="\\text{ODE Z}", cp_type=:gas),
    # Dust
    :ode_dust              => (qty_label="\\text{ODE d}", cp_type=:gas),
    # Neutral gas (everything but ionized, metals, and dust)
    :ode_neutral           => (qty_label="{\\text{ODE }\\mathrm{a + m + s}}", cp_type=:gas),
    # Molecular gas + stars
    :ode_molecular_stellar => (qty_label="{\\text{ODE }\\mathrm{m + s}}", cp_type=:gas),
    # Cold gas (everything but atomic and ionized)
    :ode_cold              => (qty_label="\\text{ODE cold}", cp_type=:gas),

    ######################
    # Wild card component
    ######################

    :generic => (qty_label="", cp_type=nothing),
)

"""
List of physical magnitudes of a simulation.
"""
const MAGNITUDES = Dict(

    #######################
    # Mass-like magnitudes
    #######################

    # Total mass [M]
    :mass => (
        label_fmt      = l -> isempty(l) ? L"M" : L"M_%$(l)",
        unit           = u"Msun",
        exp_factor     = 10,
        scatter_func   = (d, c) -> computeMass(d, c),
        integrate_func = (d, c) -> sum(computeMass(d, c); init=0.0u"Msun"),
        request        = c -> (
            if c ∈ (:stellar, :dark_matter, :black_hole, :gas)
                Dict(c => ["MASS"])
            elseif c == :Z_stellar
                Dict(:stellar => ["MASS", "GZ2 "])
            elseif c ∈ (:hydrogen, :helium)
                Dict(:gas => ["MASS"])
            elseif c == :Z_gas
                Dict(:gas => ["MASS", "GZ  "])
            elseif c ∈ (:ionized, :neutral)
                Dict(:gas => ["MASS", "NH  ", "NHP "])
            elseif c ∈ (:br_atomic, :br_molecular)
                Dict(:gas => ["MASS", "NH  ", "NHP ", "PRES"])
            elseif c ∈ (:ode_ionized, :ode_atomic, :ode_metals, :ode_dust, :ode_neutral, :ode_cold)
                Dict(:gas => ["MASS", "NH  ", "NHP ", "FRAC", "RHO ", "GZ  "])
            elseif c ∈ (:ode_molecular, :ode_stellar, :ode_molecular_stellar)
                Dict(:gas => ["MASS", "FRAC", "RHO "])
            elseif c == :generic
                Dict(
                    :gas         => ["MASS", "NH  ", "NHP ", "FRAC", "RHO ", "GZ  ", "PRES"],
                    :stellar     => ["MASS", "GZ2 "],
                    :dark_matter => ["MASS"],
                    :black_hole  => ["MASS"],
                )
            else
                throw(ArgumentError("MAGNITUDES: The given component :$(c) is not valid for the \
                magnitude :mass"))
            end
        ),
    ),

    # Volume mass density [M * L^-3]
    :mass_density => (
        label_fmt      = l -> isempty(l) ? L"\rho" : L"\rho_%$(l)",
        unit           = u"Msun * kpc^-3",
        exp_factor     = 0,
        scatter_func   = (d, c) -> computeMassDensity(d, c),
        integrate_func = (d, c) -> nothing,
        request        = c -> (
            if c ∈ (:gas, :hydrogen, :helium)
                Dict(:gas => ["MASS", "RHO "])
            elseif c == :Z_gas
                Dict(:gas => ["MASS", "RHO ", "GZ  "])
            elseif c ∈ (:ionized, :neutral)
                Dict(:gas => ["MASS", "RHO ", "NH  ", "NHP "])
            elseif c ∈ (:br_atomic, :br_molecular)
                Dict(:gas => ["MASS", "RHO ", "NH  ", "NHP ", "PRES"])
            elseif c ∈ (:ode_ionized, :ode_atomic, :ode_metals, :ode_dust, :ode_neutral, :ode_cold)
                Dict(:gas => ["MASS", "NH  ", "NHP ", "FRAC", "RHO ", "GZ  "])
            elseif c ∈ (:ode_molecular, :ode_stellar, :ode_molecular_stellar)
                Dict(:gas => ["MASS", "FRAC", "RHO "])
            elseif c == :generic
                Dict(:gas => ["MASS", "NH  ", "NHP ", "FRAC", "RHO ", "GZ  ", "PRES"])
            elseif c ∈ (:stellar, :Z_stellar, :dark_matter, :black_hole)
                nothing
            else
                throw(ArgumentError("MAGNITUDES: The given component :$(c) is not valid for the \
                magnitude :mass_density"))
            end
        ),
    ),

    # Number density [L^-3]
    :number_density => (
        label_fmt      = l -> isempty(l) ? L"n" : L"n_{\,%$(l)}",
        unit           = u"cm^-3",
        exp_factor     = 0,
        scatter_func   = (d, c) -> computeNumberDensity(d, c),
        integrate_func = (d, c) -> nothing,
        request        = c -> (
            if c ∈ (:gas, :hydrogen, :helium)
                Dict(:gas => ["MASS", "RHO "])
            elseif c == :Z_gas
                Dict(:gas => ["MASS", "RHO ", "GZ  "])
            elseif c ∈ (:ionized, :neutral)
                Dict(:gas => ["MASS", "RHO ", "NH  ", "NHP "])
            elseif c ∈ (:br_atomic, :br_molecular)
                Dict(:gas => ["MASS", "RHO ", "NH  ", "NHP ", "PRES"])
            elseif c ∈ (:ode_ionized, :ode_atomic, :ode_metals, :ode_dust, :ode_neutral, :ode_cold)
                Dict(:gas => ["MASS", "NH  ", "NHP ", "FRAC", "RHO ", "GZ  "])
            elseif c ∈ (:ode_molecular, :ode_stellar, :ode_molecular_stellar)
                Dict(:gas => ["MASS", "FRAC", "RHO "])
            elseif c == :generic
                Dict(:gas => ["MASS", "NH  ", "NHP ", "FRAC", "RHO ", "GZ  ", "PRES"])
            elseif c ∈ (:stellar, :Z_stellar, :dark_matter, :black_hole)
                nothing
            else
                throw(ArgumentError("MAGNITUDES: The given component :$(c) is not valid for the \
                magnitude :number_density"))
            end
        ),
    ),

    # Surface density [M * L^-2]
    :area_density => (
        label_fmt      = l -> isempty(l) ? L"\Sigma" : L"\Sigma_%$(l)",
        unit           = u"Msun * pc^-2",
        exp_factor     = 0,
        scatter_func   = (d, c) -> nothing,
        integrate_func = (d, c) -> nothing,
        request        = c -> (
            # See quantity3DProjection() in ./src/analysis/compute_quantities/masses.jl
            if c ∈ (:stellar, :dark_matter, :black_hole, :gas)
                Dict(c => ["MASS", "POS ", "RHO "])
            elseif c == :Z_stellar
                Dict(:stellar => ["MASS", "GZ2 ", "POS ", "RHO "])
            elseif c ∈ (:hydrogen, :helium)
                Dict(:gas => ["MASS", "POS ", "RHO "])
            elseif c == :Z_gas
                Dict(:gas => ["MASS", "GZ  ", "POS ", "RHO "])
            elseif c ∈ (:ionized, :neutral)
                Dict(:gas => ["MASS", "NH  ", "NHP ", "POS ", "RHO "])
            elseif c ∈ (:br_atomic, :br_molecular)
                Dict(:gas => ["MASS", "NH  ", "NHP ", "PRES", "POS ", "RHO "])
            elseif c ∈ (:ode_ionized, :ode_atomic, :ode_metals, :ode_dust, :ode_neutral, :ode_cold)
                Dict(:gas => ["MASS", "NH  ", "NHP ", "FRAC", "RHO ", "GZ  ", "POS "])
            elseif c ∈ (:ode_molecular, :ode_stellar, :ode_molecular_stellar)
                Dict(:gas => ["MASS", "FRAC", "RHO ", "POS "])
            elseif c == :generic
                Dict(
                    :gas         => ["MASS", "NH  ", "NHP ", "FRAC", "RHO ", "GZ  ", "PRES", "POS "],
                    :stellar     => ["MASS", "GZ2 ", "POS ", "RHO "],
                    :dark_matter => ["MASS", "POS ", "RHO "],
                    :black_hole  => ["MASS", "POS ", "RHO "],
                )
            else
                throw(ArgumentError("MAGNITUDES: The given component :$(c) is not valid for the \
                magnitude :area_density"))
            end
        ),
    ),

    # Number of elements [dimensionless]
    :number => (
        label_fmt      = l -> isempty(l) ? L"N" : L"N_%$(l)",
        unit           = Unitful.NoUnits,
        exp_factor     = 0,
        scatter_func   = (d, c) -> computeNumber(d, c),
        integrate_func = (d, c) -> (
            if c ∈ (:stellar, :dark_matter, :black_hole, :gas)
                lenght(d[c]["MASS"])
            elseif c == :Z_stellar
                lenght(d[:stellar]["MASS"])
            else
                # Gas-based components
                lenght(d[:gas]["MASS"])
            end
        ),
        request        = c -> (
            if c ∈ (:stellar, :dark_matter, :black_hole, :gas)
                Dict(c => ["MASS"])
            elseif c == :Z_stellar
                Dict(:stellar => ["MASS"])
            elseif c ∈ (:hydrogen, :helium)
                Dict(:gas => ["MASS"])
            elseif c == :Z_gas
                Dict(:gas => ["MASS", "GZ  "])
            elseif c ∈ (:ionized, :neutral)
                Dict(:gas => ["MASS", "NH  ", "NHP "])
            elseif c ∈ (:br_atomic, :br_molecular)
                Dict(:gas => ["MASS", "NH  ", "NHP ", "PRES"])
            elseif c ∈ (:ode_ionized, :ode_atomic, :ode_metals, :ode_dust, :ode_neutral, :ode_cold)
                Dict(:gas => ["MASS", "NH  ", "NHP ", "FRAC", "RHO ", "GZ  "])
            elseif c ∈ (:ode_molecular, :ode_stellar, :ode_molecular_stellar)
                Dict(:gas => ["MASS", "FRAC", "RHO "])
            elseif c == :generic
                Dict(
                    :gas         => ["MASS", "NH  ", "NHP ", "FRAC", "RHO ", "GZ  ", "PRES"],
                    :stellar     => ["MASS"],
                    :dark_matter => ["MASS"],
                    :black_hole  => ["MASS"],
                )
            else
                throw(ArgumentError("MAGNITUDES: The given component :$(c) is not valid for the \
                magnitude :number"))
            end
        ),
    ),

    # Mass fraction [dimensionless]
    :fraction => (
        label_fmt      = l -> isempty(l) ? L"f" : L"f_{\,%$(l)}",
        unit           = Unitful.NoUnits,
        exp_factor     = 0,
        scatter_func   = (d, c) -> computeFraction(d, c),
        integrate_func = (d, c) -> let

            type = (c == :Z_stellar) ? :stellar : :gas

            comp_mass = sum(computeMass(d, c); init=0.0u"Msun")
            ref_mass  = sum(computeMass(d, type); init=0.0u"Msun")

            iszero(ref_mass) ? NaN : comp_mass / ref_mass

        end,
        request        = c -> (
            if c == :Z_stellar
                Dict(:stellar => ["GZ2 "])
            elseif c ∈ (:gas, :hydrogen, :helium)
                Dict(:gas => ["MASS"])
            elseif c == :Z_gas
                Dict(:gas => ["MASS", "GZ  "])
            elseif c ∈ (:ionized, :neutral)
                Dict(:gas => ["MASS", "NH  ", "NHP "])
            elseif c ∈ (:br_atomic, :br_molecular)
                Dict(:gas => ["MASS", "NH  ", "NHP ", "PRES"])
            elseif c ∈ (:ode_ionized, :ode_atomic, :ode_metals, :ode_dust, :ode_neutral, :ode_cold)
                Dict(:gas => ["MASS", "NH  ", "NHP ", "FRAC", "RHO ", "GZ  "])
            elseif c ∈ (:ode_molecular, :ode_stellar, :ode_molecular_stellar)
                Dict(:gas => ["MASS", "FRAC", "RHO "])
            elseif c == :generic
                Dict(
                    :gas     => ["MASS", "NH  ", "NHP ", "FRAC", "RHO ", "GZ  ", "PRES"],
                    :stellar => ["GZ2 "],
                )
            elseif c ∈ (:stellar, :dark_matter, :black_hole)
                nothing
            else
                throw(ArgumentError("MAGNITUDES: The given component :$(c) is not valid for the \
                magnitude :fraction"))
            end
        ),
    ),

    # Star formation efficiency per free-fall time [dimensionless]
    :eff => (
        label_fmt      = l -> isempty(l) ? L"\epsilon_\text{ff}" : L"\epsilon_{\text{ff}, %$(l)}",
        unit           = Unitful.NoUnits,
        exp_factor     = 0,
        scatter_func   = (d, c) -> computeEfficiencyFF(d, c),
        integrate_func = (d, c) -> nothing,
        request        = c -> (
            if c == :stellar
                Dict(:stellar => ["RHOC", "GMAS", "GSFR"])
            elseif c ∈ (:gas, :hydrogen, :helium)
                Dict(:gas => ["SFR ", "MASS", "RHO "])
            elseif c == :Z_gas
                Dict(:gas => ["SFR ", "MASS", "RHO ", "GZ  "])
            elseif c ∈ (:ionized, :neutral)
                Dict(:gas => ["SFR ", "MASS", "RHO ", "NH  ", "NHP "])
            elseif c ∈ (:br_atomic, :br_molecular)
                Dict(:gas => ["SFR ", "MASS", "RHO ", "NH  ", "NHP ", "PRES"])
            elseif c ∈ (:ode_ionized, :ode_atomic, :ode_metals, :ode_dust, :ode_neutral, :ode_cold)
                Dict(:gas => ["SFR ", "MASS", "NH  ", "NHP ", "FRAC", "RHO ", "GZ  "])
            elseif c ∈ (:ode_molecular, :ode_stellar, :ode_molecular_stellar)
                Dict(:gas => ["SFR ", "MASS", "FRAC", "RHO "])
            elseif c == :generic
                Dict(
                    :gas => ["SFR ", "MASS", "NH  ", "NHP ", "FRAC", "RHO ", "GZ  ", "PRES"],
                    :stellar => ["RHOC", "GMAS", "GSFR"],
                )
            elseif c ∈ (:dark_matter, :black_hole, :Z_stellar)
                nothing
            else
                throw(ArgumentError("MAGNITUDES: The given component :$(c) is not valid for the \
                magnitude :eff"))
            end
        ),
    ),

    # Clumping factor [dimensionless]
    :clumping_factor => (
        label_fmt      = l -> isempty(l) ? L"C_\rho" : L"C_{\rho, \,%$(l)}",
        unit           = Unitful.NoUnits,
        exp_factor     = 0,
        scatter_func   = (d, c) -> nothing,
        integrate_func = (d, c) -> computeClumpingFactor(d, c),
        request        = c -> (
            if c ∈ (:gas, :hydrogen, :helium)
                Dict(:gas => ["MASS", "RHO "])
            elseif c == :Z_gas
                Dict(:gas => ["MASS", "RHO ", "GZ  "])
            elseif c ∈ (:ionized, :neutral)
                Dict(:gas => ["MASS", "RHO ", "NH  ", "NHP "])
            elseif c ∈ (:br_atomic, :br_molecular)
                Dict(:gas => ["MASS", "RHO ", "NH  ", "NHP ", "PRES"])
            elseif c ∈ (:ode_ionized, :ode_atomic, :ode_metals, :ode_dust, :ode_neutral, :ode_cold)
                Dict(:gas => ["MASS", "NH  ", "NHP ", "FRAC", "RHO ", "GZ  "])
            elseif c ∈ (:ode_molecular, :ode_stellar, :ode_molecular_stellar)
                Dict(:gas => ["MASS", "FRAC", "RHO "])
            elseif c == :generic
                Dict(:gas => ["MASS", "NH  ", "NHP ", "FRAC", "RHO ", "GZ  ", "PRES"])
            elseif c ∈ (:stellar, :Z_stellar, :dark_matter, :black_hole)
                nothing
            else
                throw(ArgumentError("MAGNITUDES: The given component :$(c) is not valid for the \
                magnitude :clumping_factor"))
            end
        ),
    ),

    #######################
    # Cinematic magnitudes
    #######################

    # Angular momentum in the z direction per unit mass [L^2 * T^-1]
    :specific_z_angular_momentum => (
        label_fmt      = l -> isempty(l) ? L"j_z" : L"j_{z, \,%$(l)}",
        unit           = u"kpc^2 * s^-1",
        exp_factor     = 0,
        scatter_func   = (d, c) -> computeSpecificAngularMomentum(d, c),
        integrate_func = (d, c) -> let

            Lz = sum(computeAngularMomentum(d, c); init=0.0u"Msun * pc^2 * yr^-1")
            M  = sum(computeMass(d, c); init=0.0u"Msun")

            Lz / M

        end,
        request        = c -> let

            blocks = ["VEL ", "POS "]
            types  = (:gas, :dark_matter, :stellar, :black_hole)

            if c ∈ types
                Dict(c => blocks)
            elseif c == :generic
                Dict(type => blocks for type in types)
            elseif c == :Z_stellar
                Dict(:stellar => blocks)
            else
                # Gas-based components
                Dict(:gas => blocks)
            end

        end,
    ),

    # Angular momentum in the z direction [M * L^2 * T^-1]
    :z_angular_momentum => (
        label_fmt      = l -> isempty(l) ? L"l_z" : L"l_{z, \,%$(l)}",
        unit           = u"Msun * kpc^2 * s^-1",
        exp_factor     = 0,
        scatter_func   = (d, c) -> computeAngularMomentum(d, c),
        integrate_func = (d, c) -> sum(computeAngularMomentum(d, c); init=0.0u"Msun * pc^2 * yr^-1"),
        request        = c -> let

            blocks = ["VEL ", "POS ", "MASS"]
            types  = (:gas, :dark_matter, :stellar, :black_hole)

            if c ∈ types
                Dict(c => blocks)
            elseif c == :generic
                Dict(type => blocks for type in types)
            elseif c == :Z_stellar
                Dict(:stellar => blocks)
            else
                # Gas-based components
                Dict(:gas => blocks)
            end

        end,
    ),

    # Spin parameter [dimensionless]
    :spin_parameter => (
        label_fmt      = l -> isempty(l) ? L"\lambda" : L"\lambda_%$(l)",
        unit           = Unitful.NoUnits,
        exp_factor     = 0,
        scatter_func   = (d, c) -> nothing,
        integrate_func = (d, c) -> computeSpinParameter(d, c),
        request        = c -> let

            blocks = ["VEL ", "POS ", "MASS"]
            types  = (:gas, :dark_matter, :stellar, :black_hole)

            if c ∈ types
                Dict(c => blocks)
            elseif c == :generic
                Dict(type => blocks for type in types)
            elseif c == :Z_stellar
                Dict(:stellar => blocks)
            else
                # Gas-based components
                Dict(:gas => blocks)
            end

        end,
    ),

    # Circularity [dimensionless]
    :circularity => (
        label_fmt      = l -> isempty(l) ? L"\epsilon" : L"\epsilon_{\,%$(l)}",
        unit           = Unitful.NoUnits,
        exp_factor     = 0,
        scatter_func   = (d, c) -> computeCircularity(d, c),
        integrate_func = (d, c) -> nothing,
        request        = c -> Dict(t => ["VEL ", "POS ", "MASS"] for t in keys(PARTICLE_INDEX)),
    ),

    # Circular velocity [L * T^-1]
    :circular_velocity => (
        label_fmt      = l -> isempty(l) ? L"v_\text{circ}" : L"v_{\text{circ}, \,%$(l)}",
        unit           = u"km * s^-1",
        exp_factor     = 0,
        scatter_func   = (d, c) -> computeVcirc(d, c)[2],
        integrate_func = (d, c) -> nothing,
        request        = c -> Dict(t => ["MASS", "POS "] for t in keys(PARTICLE_INDEX)),
    ),

    # Component of the velocity in the radial direction [L * T^-1]
    :radial_velocity => (
        label_fmt      = l -> isempty(l) ? L"v_r" : L"v_{r, \,%$(l)}",
        unit           = u"km * s^-1",
        exp_factor     = 0,
        scatter_func   = (d, c) -> computeVpolar(d, c, :radial),
        integrate_func = (d, c) -> nothing,
        request        = c -> let

            blocks = ["VEL ", "POS "]
            types  = (:gas, :dark_matter, :stellar, :black_hole)

            if c ∈ types
                Dict(c => blocks)
            elseif c == :generic
                Dict(type => blocks for type in types)
            else
                # Gas-based components
                Dict(:gas => blocks)
            end

        end,
    ),

    # Component of the velocity in the tangential direction [L * T^-1]
    :tangential_velocity => (
        label_fmt      = l -> isempty(l) ? L"v_\theta" : L"v_{\theta, \,%$(l)}",
        unit           = u"km * s^-1",
        exp_factor     = 0,
        scatter_func   = (d, c) -> computeVpolar(d, c, :tangential),
        integrate_func = (d, c) -> nothing,
        request        = c -> let

            blocks = ["VEL ", "POS "]
            types  = (:gas, :dark_matter, :stellar, :black_hole)

            if c ∈ types
                Dict(c => blocks)
            elseif c == :generic
                Dict(type => blocks for type in types)
            else
                # Gas-based components
                Dict(:gas => blocks)
            end

        end,
    ),

    # Component of the velocity in the z direction [L * T^-1]
    :zstar_velocity => (
        label_fmt      = l -> isempty(l) ? L"v_z \, \mathrm{sign}(z)" : L"v_{z, \,%$(l)} \, \mathrm{sign}(z)",
        unit           = u"km * s^-1",
        exp_factor     = 0,
        scatter_func   = (d, c) -> computeVpolar(d, c, :zstar),
        integrate_func = (d, c) -> nothing,
        request        = c -> let

            blocks = ["VEL ", "POS "]
            types  = (:gas, :dark_matter, :stellar, :black_hole)

            if c ∈ types
                Dict(c => blocks)
            elseif c == :generic
                Dict(type => blocks for type in types)
            else
                # Gas-based components
                Dict(:gas => blocks)
            end

        end,
    ),

    # Kinetic energy [M * L^2 * T^-2]
    :kinetic_energy => (
        label_fmt      = l -> isempty(l) ? L"E_k" : L"E_{k, \,%$(l)}",
        unit           = u"erg",
        exp_factor     = 51,
        scatter_func   = (d, c) -> computeKineticEnergy(d, c),
        integrate_func = (d, c) -> nothing,
        request        = c -> (
            if c ∈ (:stellar, :dark_matter, :black_hole, :gas)
                Dict(c => ["VEL ", "MASS"])
            elseif c == :Z_stellar
                Dict(:stellar => ["VEL ", "MASS", "GZ2 "])
            elseif c ∈ (:hydrogen, :helium)
                Dict(:gas => ["VEL ", "MASS"])
            elseif c == :Z_gas
                Dict(:gas => ["VEL ", "MASS", "GZ  "])
            elseif c ∈ (:ionized, :neutral)
                Dict(:gas => ["VEL ", "MASS", "NH  ", "NHP "])
            elseif c ∈ (:br_atomic, :br_molecular)
                Dict(:gas => ["VEL ", "MASS", "NH  ", "NHP ", "PRES"])
            elseif c ∈ (:ode_ionized, :ode_atomic, :ode_metals, :ode_dust, :ode_neutral, :ode_cold)
                Dict(:gas => ["VEL ", "MASS", "NH  ", "NHP ", "FRAC", "RHO ", "GZ  "])
            elseif c ∈ (:ode_molecular, :ode_stellar, :ode_molecular_stellar)
                Dict(:gas => ["VEL ", "MASS", "FRAC", "RHO "])
            elseif c == :generic
                Dict(
                    :gas         => ["VEL ", "MASS", "NH  ", "NHP ", "FRAC", "RHO ", "GZ  ", "PRES"],
                    :stellar     => ["VEL ", "MASS", "GZ2 "],
                    :dark_matter => ["VEL ", "MASS"],
                    :black_hole  => ["VEL ", "MASS"],
                )
            else
                throw(ArgumentError("MAGNITUDES: The given component :$(c) is not valid for the \
                magnitude :kinetic_energy"))
            end
        ),
    ),

    # Potential energy [M * L^2 * T^-2]
    :potential_energy => (
        label_fmt      = l -> isempty(l) ? L"E_p" : L"E_{p, \,%$(l)}",
        unit           = u"erg",
        exp_factor     = 51,
        scatter_func   = (d, c) -> computePotentialEnergy(d, c),
        integrate_func = (d, c) -> sum(computePotentialEnergy(d, c); init=0.0u"erg"),
        request        = c -> (
            if c ∈ (:stellar, :dark_matter, :black_hole, :gas)
                Dict(c => ["POT ", "MASS"])
            elseif c == :Z_stellar
                Dict(:stellar => ["POT ", "MASS", "GZ2 "])
            elseif c ∈ (:hydrogen, :helium)
                Dict(:gas => ["POT ", "MASS"])
            elseif c == :Z_gas
                Dict(:gas => ["POT ", "MASS", "GZ  "])
            elseif c ∈ (:ionized, :neutral)
                Dict(:gas => ["POT ", "MASS", "NH  ", "NHP "])
            elseif c ∈ (:br_atomic, :br_molecular)
                Dict(:gas => ["POT ", "MASS", "NH  ", "NHP ", "PRES"])
            elseif c ∈ (:ode_ionized, :ode_atomic, :ode_metals, :ode_dust, :ode_neutral, :ode_cold)
                Dict(:gas => ["POT ", "MASS", "NH  ", "NHP ", "FRAC", "RHO ", "GZ  "])
            elseif c ∈ (:ode_molecular, :ode_stellar, :ode_molecular_stellar)
                Dict(:gas => ["POT ", "MASS", "FRAC", "RHO "])
            elseif c == :generic
                Dict(
                    :gas         => ["POT ", "MASS", "NH  ", "NHP ", "FRAC", "RHO ", "GZ  ", "PRES"],
                    :stellar     => ["POT ", "MASS", "GZ2 "],
                    :dark_matter => ["POT ", "MASS"],
                    :black_hole  => ["POT ", "MASS"],
                )
            else
                throw(ArgumentError("MAGNITUDES: The given component :$(c) is not valid for the \
                magnitude :potential_energy"))
            end
        ),
    ),

    # Total energy = kinetic + potential [M * L^2 * T^-2]
    :total_energy => (
        label_fmt      = l -> isempty(l) ? L"E" : L"E_{%$(l)}",
        unit           = u"erg",
        exp_factor     = 51,
        scatter_func   = (d, c) -> computeTotalEnergy(d, c),
        integrate_func = (d ,c) -> nothing,
        request        = c -> (
            if c ∈ (:stellar, :dark_matter, :black_hole, :gas)
                Dict(c => ["VEL ", "POT ", "MASS"])
            elseif c == :Z_stellar
                Dict(:stellar => ["VEL ", "POT ", "MASS", "GZ2 "])
            elseif c ∈ (:hydrogen, :helium)
                Dict(:gas => ["VEL ", "POT ", "MASS"])
            elseif c == :Z_gas
                Dict(:gas => ["VEL ", "POT ", "MASS", "GZ  "])
            elseif c ∈ (:ionized, :neutral)
                Dict(:gas => ["VEL ", "POT ", "MASS", "NH  ", "NHP "])
            elseif c ∈ (:br_atomic, :br_molecular)
                Dict(:gas => ["VEL ", "POT ", "MASS", "NH  ", "NHP ", "PRES"])
            elseif c ∈ (:ode_ionized, :ode_atomic, :ode_metals, :ode_dust, :ode_neutral, :ode_cold)
                Dict(:gas => ["VEL ", "POT ", "MASS", "NH  ", "NHP ", "FRAC", "RHO ", "GZ  "])
            elseif c ∈ (:ode_molecular, :ode_stellar, :ode_molecular_stellar)
                Dict(:gas => ["VEL ", "POT ", "MASS", "FRAC", "RHO "])
            elseif c == :generic
                Dict(
                    :gas         => ["VEL ", "POT ", "MASS", "NH  ", "NHP ", "FRAC", "RHO ", "GZ  ", "PRES"],
                    :stellar     => ["VEL ", "POT ", "MASS", "GZ2 "],
                    :dark_matter => ["VEL ", "POT ", "MASS"],
                    :black_hole  => ["VEL ", "POT ", "MASS"],
                )
            else
                throw(ArgumentError("MAGNITUDES: The given component :$(c) is not valid for the \
                magnitude :total_energy"))
            end
        ),
    ),

    ########
    # Other
    ########

    # Depletion time [T]
    :depletion_time => (
        label_fmt      = l -> isempty(l) ? L"\tau_\text{dep}" : L"\tau_{\text{dep}, \,%$(l)}",
        unit           = u"Gyr",
        exp_factor     = 0,
        scatter_func   = (d, c) -> computeDepletionTime(d, c),
        integrate_func = (d, c) -> let

            mass    = sum(computeMass(d, c); init=0.0u"Msun")
            cp_type = COMPONENTS[c].cp_type

            if cp_type == :stellar
                sfr = sum(scatterQty(d, :sfr); init=0.0u"Msun * yr^-1")
            elseif cp_type == :gas
                sfr = sum(scatterQty(d, :gas_sfr); init=0.0u"Msun * yr^-1")
            else
                throw(ArgumentError("integrateQty: I don't know how to compute the SFR for \
                cp_type :$(cp_type)"))
            end

            mass / sfr

        end,
        request        = c -> (
            if c ∈ (:gas, :hydrogen, :helium)
                Dict(:gas => ["SFR ", "MASS"])
            elseif c == :Z_gas
                Dict(:gas => ["SFR ", "MASS", "GZ  "])
            elseif c ∈ (:ionized, :neutral)
                Dict(:gas => ["SFR ", "MASS", "NH  ", "NHP "])
            elseif c ∈ (:br_atomic, :br_molecular)
                Dict(:gas => ["SFR ", "MASS", "NH  ", "NHP ", "PRES"])
            elseif c ∈ (:ode_ionized, :ode_atomic, :ode_metals, :ode_dust, :ode_neutral, :ode_cold)
                Dict(:gas => ["SFR ", "MASS", "NH  ", "NHP ", "FRAC", "RHO ", "GZ  "])
            elseif c ∈ (:ode_molecular, :ode_stellar, :ode_molecular_stellar)
                Dict(:gas => ["SFR ", "MASS", "FRAC", "RHO "])
            elseif c == :generic
                Dict(:gas => ["SFR ", "MASS", "NH  ", "NHP ", "FRAC", "RHO ", "GZ  ", "PRES"])
            elseif c ∈ (:stellar, :Z_stellar, :dark_matter, :black_hole)
                nothing
            else
                throw(ArgumentError("MAGNITUDES: The given component :$(c) is not valid for the \
                magnitude :depletion_time"))
            end
        ),
    ),

    # Projected (to the xy plane) distance to the origin [L]
    :xy_distance => (
        label_fmt      = l -> isempty(l) ? L"d_{xy}" : L"d_{xy, \,%$(l)}",
        unit           = u"kpc",
        exp_factor     = 0,
        scatter_func   = (d, c) -> computeXYDistance(d, c),
        integrate_func = (d, c) -> nothing,
        request        = c -> let

            blocks = ["POS "]
            types  = (:gas, :dark_matter, :stellar, :black_hole)

            if c ∈ types
                Dict(c => blocks)
            elseif c == :generic
                Dict(type => blocks for type in types)
            elseif c == :Z_stellar
                Dict(:stellar => blocks)
            else
                # Gas-based components
                Dict(:gas => blocks)
            end

        end,
    ),

    # Distance to the origin [L]
    :radial_distance => (
        label_fmt      = l -> isempty(l) ? L"r" : L"r_{\,%$(l)}",
        unit           = u"kpc",
        exp_factor     = 0,
        scatter_func   = (d, c) -> computeRadialDistance(d, c),
        integrate_func = (d, c) -> nothing,
        request        = c -> let

            blocks = ["POS "]
            types  = (:gas, :dark_matter, :stellar, :black_hole)

            if c ∈ types
                Dict(c => blocks)
            elseif c == :generic
                Dict(type => blocks for type in types)
            elseif c == :Z_stellar
                Dict(:stellar => blocks)
            else
                # Gas-based components
                Dict(:gas => blocks)
            end

        end,
    ),

)

"""
List of quantities related to chemical abundance.
"""
const ABUNDANCE_QUANTITIES = Dict(

    :gas_mass          => (
        label_fmt      = e -> L"M_%$(e)^\text{gas}",
        request        = Dict(:gas => ["MASS", "GMET"]),
        unit           = u"Msun",
        exp_factor     = 5,
        cp_type        = :gas,
        scatter_func   = (d, e) -> computeElementMass(d, :gas, e),
        integrate_func = (d, e) -> sum(computeElementMass(d, :gas, e); init=0.0u"Msun"),
    ),

    :stellar_mass      => (
        label_fmt      = e -> L"M_%$(e)^\star",
        request        = Dict(:stellar => ["MASS", "GME2"]),
        unit           = u"Msun",
        exp_factor     = 8,
        cp_type        = :stellar,
        scatter_func   = (d, e) -> computeElementMass(d, :stellar, e),
        integrate_func = (d, e) -> sum(computeElementMass(d, :stellar, e); init=0.0u"Msun"),
    ),

    :gas_abundance     => (
        label_fmt    = e -> let

            shift = ABUNDANCE_SHIFT[][e]

            if iszero(shift)
                L"\log_{10}(\mathrm{%$e} \, / \, \mathrm{H})"
            else
                L"%$(Int(shift)) + \log_{10}(\mathrm{%$e} \, / \, \mathrm{H})"
            end

        end,
        request      = Dict(:gas => ["MASS", "GMET"]),
        unit         = Unitful.NoUnits,
        exp_factor   = 0,
        cp_type      = :gas,
        scatter_func = (d, e) -> let

            abundances = computeAbundance(d, :gas, e; solar=false)

            replace!(x -> isinf(x) ? NaN : x, ABUNDANCE_SHIFT[][e] .+ log10.(abundances))

        end,
        integrate_func = (d, e) -> let

            abundance = computeGlobalAbundance(d, :gas, e; solar=false)

            if iszero(abundance)
                NaN
            else
                ABUNDANCE_SHIFT[][e] + log10(abundance)
            end

        end,
    ),

    :stellar_abundance => (
        label_fmt    = e -> let

            shift = ABUNDANCE_SHIFT[][e]

            if iszero(shift)
                L"\log_{10}(\mathrm{%$e} \, / \, \mathrm{H})"
            else
                L"%$(Int(shift)) + \log_{10}(\mathrm{%$e} \, / \, \mathrm{H})"
            end

        end,
        request      = Dict(:stellar => ["MASS", "GME2"]),
        unit         = Unitful.NoUnits,
        exp_factor   = 0,
        cp_type      = :stellar,
        scatter_func = (d, e) -> let

            abundances = computeAbundance(d, :stellar, e; solar=false)

            replace!(x -> isinf(x) ? NaN : x, ABUNDANCE_SHIFT[][e] .+ log10.(abundances))

        end,
                integrate_func = (d, e) -> let

            abundance = computeGlobalAbundance(d, :stellar, e; solar=false)

            if iszero(abundance)
                NaN
            else
                ABUNDANCE_SHIFT[][e] + log10(abundance)
            end

        end,
    ),

)

"""
List of halo quantities.
"""
const HALO_QUANTITIES = Dict(

    :halo_mass       => (
        qty_label  = L"M_\text{halo}",
        unit       = u"Msun",
        exp_factor = 10,
        snap_key   = "G_Mass",
    ),

    :halo_n_subhalos => (
        qty_label  = L"N_\text{subhalos}",
        unit       = Unitful.NoUnits,
        exp_factor = 0,
        snap_key   = "G_Nsubs",
    ),

    :halo_M200       => (
        qty_label  = L"M_{200}",
        unit       = u"Msun",
        exp_factor = 10,
        snap_key   = "G_M_Crit200",
    ),

    :halo_R200       => (
        qty_label  = L"R_{200}",
        unit       = u"kpc",
        exp_factor = 0,
        snap_key   = "G_R_Crit200",
    ),

)

"""
List of physical components of our SF model.
"""
const SFM_COMPONENTS = Dict(

    # Gas cells
    :ode_gas     => (qty_label="\\,\\text{gas}", cp_type=:gas),
    # Stellar particles
    :ode_stellar => (qty_label="\\star", cp_type=:stellar),

)

"""
List of physical magnitudes of our SF model.
"""
const SFM_MAGNITUDES = Dict(

    # Integration time [T]
    :integration_time => (
        label_fmt      = l -> L"t_{i}^{%$(l)}",
        unit           = u"Myr",
        exp_factor     = 0,
        request        = c -> Dict(c => ["ODIT"]),
        scatter_func   = (d, c) -> d[c]["ODIT"],
        integrate_func = (d, c) -> nothing,
        valid_types    = (:gas, :stellar),

    ),

    # Scale factor [dimensionless]
    :parameter_a => (
        label_fmt      = l -> L"a^{%$(l)}",
        unit           = Unitful.NoUnits,
        exp_factor     = 0,
        request        = c -> Dict(c => ["PARA"]),
        scatter_func   = (d, c) -> d[c]["PARA"],
        integrate_func = (d, c) -> nothing,
        valid_types    = (:gas, :stellar),
    ),

    # Redshift [dimensionless]
    :parameter_z => (
        label_fmt      = l -> L"z^{%$(l)}",
        unit           = Unitful.NoUnits,
        exp_factor     = 0,
        request        = c -> Dict(c => ["PARz"]),
        scatter_func   = (d, c) -> d[c]["PARz"],
        integrate_func = (d, c) -> nothing,
        valid_types    = (:gas, :stellar),

    ),

    # UVB photoionization rate [T^-1]
    :parameter_uvb => (
        label_fmt      = l -> L"\text{UVB}^{%$(l)}",
        unit           = u"Myr^-1",
        exp_factor     = 0,
        request        = c -> Dict(c => ["PARU"]),
        scatter_func   = (d, c) -> d[c]["PARU"],
        integrate_func = (d, c) -> nothing,
        valid_types    = (:gas, :stellar),

    ),

    # LWB photodissociation rate [T^-1]
    :parameter_lwb => (
        label_fmt      = l -> L"\text{LWB}^{%$(l)}",
        unit           = u"Myr^-1",
        exp_factor     = 0,
        request        = c -> Dict(c => ["PARL"]),
        scatter_func   = (d, c) -> d[c]["PARL"],
        integrate_func = (d, c) -> nothing,
        valid_types    = (:gas, :stellar),

    ),

    # Star formation time parameter [T]
    :tau_s => (
        label_fmt      = l -> L"\tau_\text{star}^{%$(l)}",
        unit           = u"Myr",
        exp_factor     = 0,
        request        = c -> Dict(c => ["TAUS"]),
        scatter_func   = (d, c) -> d[c]["TAUS"],
        integrate_func = (d, c) -> nothing,
        valid_types    = (:gas, :stellar),

    ),

    # Gas density [L^-3]
    :parameter_cell_density => (
        label_fmt      = l -> L"\rho_{c}^{%$(l)}",
        unit           = u"cm^-3",
        exp_factor     = 0,
        request        = c -> Dict(c => ["RHOC"]),
        scatter_func   = (d, c) -> d[c]["RHOC"],
        integrate_func = (d, c) -> nothing,
        valid_types    = (:gas, :stellar),

    ),

    # Gas metallicity [dimensionless]
    :parameter_metallicity => (
        label_fmt      = l -> L"Z^{%$(l)}",
        unit           = Unitful.NoUnits,
        exp_factor     = 0,
        request        = c -> Dict(c => ["PARZ"]),
        scatter_func   = (d, c) -> d[c]["PARZ"],
        integrate_func = (d, c) -> let

            if c == :gas
                gas_mass = d[c]["MASS"]
            elseif c == :stellar
                gas_mass = d[c]["GMAS"]
            end

            Z = d[c]["PARZ"]

            Mg = nansum(gas_mass; init=0.0u"Msun")
            MZ = nansum(Z .* gas_mass; init=0.0u"Msun")

            if iszero(Mg)
                NaN
            else
                (MZ / Mg) / SOLAR_METALLICITY[]
            end

        end,
        valid_types    = (:gas, :stellar),

    ),

    # Column height [L]
    :parameter_column_height => (
        label_fmt      = l -> L"h^{%$(l)}",
        unit           = u"pc",
        exp_factor     = 0,
        request        = c -> Dict(c => ["PARH"]),
        scatter_func   = (d, c) -> d[c]["PARH"],
        integrate_func = (d, c) -> nothing,
        valid_types    = (:gas, :stellar),

    ),

    # Photodissociation efficiency [dimensionless]
    :parameter_eta_d => (
        label_fmt      = l -> L"\eta_{\,\text{diss}}^{%$(l)}",
        unit           = Unitful.NoUnits,
        exp_factor     = 0,
        request        = c -> Dict(c => ["ETAD"]),
        scatter_func   = (d, c) -> d[c]["ETAD"],
        integrate_func = (d, c) -> nothing,
        valid_types    = (:gas, :stellar),

    ),

    # Photoionization efficiency [dimensionless]
    :parameter_eta_i => (
        label_fmt      = l -> L"\eta_{\,\text{ion}}^{%$(l)}",
        unit           = Unitful.NoUnits,
        exp_factor     = 0,
        request        = c -> Dict(c => ["ETAI"]),
        scatter_func   = (d, c) -> d[c]["ETAI"],
        integrate_func = (d, c) -> nothing,
        valid_types    = (:gas, :stellar),
    ),

    # Mass recycling fraction [dimensionless]
    :parameter_r => (
        label_fmt      = l -> L"R^{%$(l)}",
        unit           = Unitful.NoUnits,
        exp_factor     = 0,
        request        = c -> Dict(c => ["PARR"]),
        scatter_func   = (d, c) -> d[c]["PARR"],
        integrate_func = (d, c) -> nothing,
        valid_types    = (:gas, :stellar),

    ),

    # Metallicity of the supernova ejecta [dimensionless]
    :parameter_zsn => (
        label_fmt      = l -> L"Z_\text{SN}^{%$(l)}",
        unit           = Unitful.NoUnits,
        exp_factor     = 0,
        request        = c -> Dict(c => ["PAZN"]),
        scatter_func   = (d, c) -> d[c]["PAZN"],
        integrate_func = (d, c) -> nothing,
        valid_types    = (:gas, :stellar),

    ),

    # Star formation flag [dimensionless]
    :sf_flag => (
        label_fmt      = l -> L"\text{SF}_\text{flag}^{%$(l)}",
        unit           = Unitful.NoUnits,
        exp_factor     = 0,
        request        = c -> Dict(c => ["SFFL"]),
        scatter_func   = (d, c) -> d[c]["SFFL"],
        integrate_func = (d, c) -> nothing,
        valid_types    = (:gas,),

    ),

    # Cold gas fraction [dimensionless]
    :cold_mass_frac => (
        label_fmt      = l -> L"f_\text{cold}^{%$(l)}",
        unit           = Unitful.NoUnits,
        exp_factor     = 0,
        request        = c -> Dict(c => ["COLF"]),
        scatter_func   = (d, c) -> d[c]["COLF"],
        integrate_func = (d, c) -> let

            if c == :gas
                gas_mass = d[c]["MASS"]
            elseif c == :stellar
                gas_mass = d[c]["GMAS"]
            end

            cf = d[c]["COLF"]

            Mg = nansum(gas_mass; init=0.0u"Msun")
            Mc = nansum(cf .* gas_mass; init=0.0u"Msun")

            if iszero(Mg)
                NaN
            else
                Mc / Mg
            end

        end,
        valid_types    = (:gas, :stellar),

    ),

    # Parent gas cell mass [M]
    :gas_mass => (
        label_fmt      = l -> L"M_\text{gas}^{\,%$(l)}",
        unit           = u"Msun",
        exp_factor     = 0,
        request        = c -> Dict(c => ["GMAS"]),
        scatter_func   = (d, c) -> d[c]["GMAS"],
        integrate_func = (d, c) -> nansum(d[c]["GMAS"]; init=0.0u"Msun"),
        valid_types    = (:stellar,),

    ),

    # Parent gas cell SFR [M * T^-1]
    :gas_sfr => (
        label_fmt      = l -> L"\text{SFR}_\text{gas}^{%$(l)}",
        unit           = u"Msun * yr^-1",
        exp_factor     = 0,
        request        = c -> Dict(c => ["GSFR"]),
        scatter_func   = (d, c) -> d[c]["GSFR"],
        integrate_func = (d, c) -> nothing,
        valid_types    = (:stellar,),

    ),

    # Parent gas cell pressure [M * L^-1 * T^-2]
    :gas_pressure => (
        label_fmt      = l -> L"P^{%$(l)}",
        unit           = u"Pa",
        exp_factor     = 0,
        request        = c -> Dict(c => ["GPRE"]),
        scatter_func   = (d, c) -> d[c]["GPRE"],
        integrate_func = (d, c) -> nothing,
        valid_types    = (:stellar,),

    ),

    # SNII fraction [M^-1]
    :parameter_snii => (
        label_fmt      = l -> L"\text{SNII}_\text{frac}^{%$(l)}",
        unit           = u"Msun^-1",
        exp_factor     = 0,
        request        = c -> Dict(c => ["PARS"]),
        scatter_func   = (d, c) -> d[c]["PARS"],
        integrate_func = (d, c) -> nothing,
        valid_types    = (:gas, :stellar),

    ),

    # Star formation timescale [T]
    :tau_star => (
        label_fmt      = l -> L"\tau_\text{star}^{%$(l)}",
        unit           = u"Myr",
        exp_factor     = 0,
        request        = c -> Dict(c => ["RHOC"]),
        scatter_func   = (d, c) -> let

            ρ   = d[c]["RHOC"] .* u"mp"
            tff = @. sqrt(3π / (32 * Unitful.G * ρ))

            @. tff / AREPO_SFM.εff

        end,
        integrate_func = (d, c) -> nothing,
        valid_types    = (:gas, :stellar),
    ),

    # Recombination timescale [T]
    :tau_rec => (
        label_fmt      = l -> L"\tau_\text{rec}^{%$(l)}",
        unit           = u"Myr",
        exp_factor     = 0,
        request        = c -> Dict(c => ["RHOC", "FRAC"]),
        scatter_func   = (d, c) -> let

            fractions = d[c]["FRAC"]

            if isempty(fractions)

                typeof(1.0u"Myr")[]

            else

                ρ   = d[c]["RHOC"] .* u"mp"
                nan = NaN * u"Myr"

                fi = view(fractions, SFM_IDX[:ode_ionized], :)

                tau_rec = @. Unitful.mp / (AREPO_SFM.αH * fi * ρ)
                replace!(x -> isinf(x) ? nan : x, tau_rec)

            end

        end,
        integrate_func = (d, c) -> nothing,
        valid_types    = (:gas, :stellar),
    ),

    # Condensation timescale [T]
    :tau_cond => (
        label_fmt      = l -> L"\tau_\text{cond}^{%$(l)}",
        unit           = u"Myr",
        exp_factor     = 0,
        request        = c -> Dict(c => ["RHOC", "FRAC"]),
        scatter_func   = (d, c) -> let

            fractions = d[c]["FRAC"]

            if isempty(fractions)

                typeof(1.0u"Myr")[]

            else

                ρ   = d[c]["RHOC"] .* u"mp"
                nan = NaN * u"Myr"

                fi = view(fractions, SFM_IDX[:ode_ionized], :)
                fa = view(fractions, SFM_IDX[:ode_atomic], :)
                fm = view(fractions, SFM_IDX[:ode_molecular], :)

                fg = @. fi + fa + fm

                tau_cond = @. Unitful.mp / (2 * AREPO_SFM.Rd * ρ * fg)

                replace!(x -> isinf(x) ? nan : x, tau_cond)

            end

        end,
        integrate_func = (d, c) -> nothing,
        valid_types    = (:gas, :stellar),
    ),

    # Dust growth timescale [T]
    :tau_dg => (
        label_fmt      = l -> L"\tau_\text{dg}^{%$(l)}",
        unit           = u"Gyr",
        exp_factor     = 0,
        request        = c -> Dict(c => ["RHOC", "FRAC"]),
        scatter_func   = (d, c) -> let

            fractions = d[c]["FRAC"]

            if isempty(fractions)

                typeof(1.0u"Gyr")[]

            else

                ρ   = d[c]["RHOC"] .* u"mp"
                nan = NaN * u"Gyr"

                fa = view(fractions, SFM_IDX[:ode_atomic], :)
                fm = view(fractions, SFM_IDX[:ode_molecular], :)
                fZ = view(fractions, SFM_IDX[:ode_metals], :)
                fd = view(fractions, SFM_IDX[:ode_dust], :)

                fn = @. fa + fm
                Z  = @. fZ + fd

                tau_dg = @. 1.0 / (AREPO_SFM.Cdg * Z * ρ * fn)

                replace!(x -> isinf(x) ? nan : x, tau_dg)

            end

        end,
        integrate_func = (d, c) -> nothing,
        valid_types    = (:gas, :stellar),
    ),

    # Dust creation timescale [T]
    :tau_dc => (
        label_fmt      = l -> L"\tau_\text{dc}^{%$(l)}",
        unit           = u"Gyr",
        exp_factor     = 0,
        request        = c -> Dict(c => ["RHOC", "FRAC"]),
        scatter_func   = (d, c) -> let

            fractions = d[c]["FRAC"]

            if isempty(fractions)

                typeof(1.0u"Gyr")[]

            else

                ρ   = d[c]["RHOC"] .* u"mp"
                nan = NaN * u"Gyr"

                fa = view(fractions, SFM_IDX[:ode_atomic], :)
                fm = view(fractions, SFM_IDX[:ode_molecular], :)
                fZ = view(fractions, SFM_IDX[:ode_metals], :)

                fn = @. fa + fm

                tau_dc = @. 1.0 / (AREPO_SFM.Cdg * fZ * fn * fn * ρ)

                replace!(x -> isinf(x) ? nan : x, tau_dc)

            end

        end,
        integrate_func = (d, c) -> nothing,
        valid_types    = (:gas, :stellar),
    ),

    # Dust destruction timescale [T]
    :tau_dd => (
        label_fmt      = l -> L"\tau_\text{dd}^{%$(l)}",
        unit           = u"Gyr",
        exp_factor     = 0,
        request        = c -> Dict(c => ["RHOC", "FRAC", "PARS"]),
        scatter_func   = (d, c) -> let

            fractions = d[c]["FRAC"]

            if isempty(fractions)

                typeof(1.0u"Gyr")[]

            else

                ρ    = d[c]["RHOC"] .* u"mp"
                nan  = NaN * u"Gyr"
                SNII = d[c]["PARS"]

                fm    = view(fractions, SFM_IDX[:ode_molecular], :)
                tff   = @. sqrt(3π / (32 * Unitful.G * ρ))
                τstar = @. tff / AREPO_SFM.εff
                ψ     = @. fm / τstar

                tau_dd = @. 1.0 / (SNII * AREPO_SFM.Cswm * ψ)

                replace!(x -> isinf(x) ? nan : x, tau_dd)

            end

        end,
        integrate_func = (d, c) -> nothing,
        valid_types    = (:gas, :stellar),
    ),

    # Ionization optical depth [dimensionless]
    :tau_ion => (
        label_fmt      = l -> L"\tau_\text{ion}^{%$(l)}",
        unit           = Unitful.NoUnits,
        exp_factor     = 0,
        request        = c -> Dict(c => ["RHOC", "PARH", "FRAC"]),
        scatter_func   = (d, c) -> let

            fractions = d[c]["FRAC"]

            if isempty(fractions)

                Float64[]

            else

                ρ  = d[c]["RHOC"] .* u"mp"
                h  = d[c]["PARH"]
                fa = view(fractions, SFM_IDX[:ode_atomic], :)

                @. (fa * ρ / Unitful.mp) * AREPO_SFM.σion * h

            end

        end,
        integrate_func = (d, c) -> nothing,
        valid_types    = (:gas, :stellar),
    ),

    # Dissociation optical depth [dimensionless]
    :tau_diss => (
        label_fmt      = l -> L"\tau_\text{diss}^{%$(l)}",
        unit           = Unitful.NoUnits,
        exp_factor     = 0,
        request        = c -> Dict(c => ["RHOC", "PARH", "FRAC"]),
        scatter_func   = (d, c) -> let

            fractions = d[c]["FRAC"]

            if isempty(fractions)

                Float64[]

            else

                ρ  = d[c]["RHOC"] .* u"mp"
                h  = d[c]["PARH"]
                fm = view(fractions, SFM_IDX[:ode_molecular], :)

                @. (fm * ρ / (2 * Unitful.mp)) * AREPO_SFM.σdiss * h

            end

        end,
        integrate_func = (d, c) -> nothing,
        valid_types    = (:gas, :stellar),
    ),

    # Dust shielding factor [dimensionless]
    :S_d => (
        label_fmt      = l -> L"S_{d}^{%$(l)}",
        unit           = Unitful.NoUnits,
        exp_factor     = 0,
        request        = c -> Dict(c => ["RHOC", "PARH", "FRAC"]),
        scatter_func   = (d, c) -> let

            fractions = d[c]["FRAC"]

            if isempty(fractions)

                Float64[]

            else

                ρ = d[c]["RHOC"] .* u"mp"
                h = d[c]["PARH"]

                fa = view(fractions, SFM_IDX[:ode_atomic], :)
                fm = view(fractions, SFM_IDX[:ode_molecular], :)
                fZ = view(fractions, SFM_IDX[:ode_metals], :)
                fd = view(fractions, SFM_IDX[:ode_dust], :)

                Z = @. fZ + fd

                N  = @. (fa + fm) * ρ * h / Unitful.mp
                Zs = @. Z / SOLAR_METALLICITY[]

                @. exp(-AREPO_SFM.σd * Zs * N)

            end

        end,
        integrate_func = (d, c) -> nothing,
        valid_types    = (:gas, :stellar),
    ),

    # H2 self-shielding factor [dimensionless]
    :S_H2 => (
        label_fmt      = l -> L"S_{H2}^{%$(l)}",
        unit           = Unitful.NoUnits,
        exp_factor     = 0,
        request        = c -> Dict(c => ["RHOC", "PARH", "FRAC"]),
        scatter_func   = (d, c) -> let

            fractions = d[c]["FRAC"]

            if isempty(fractions)

                Float64[]

            else

                ρ  = d[c]["RHOC"] .* u"mp"
                h  = d[c]["PARH"]
                fm = view(fractions, SFM_IDX[:ode_molecular], :)

                x    = @. (fm * ρ * h) / (2 * Unitful.mp * AREPO_SFM.xnf)
                x1   = @. 1.0 + x
                sqx1 = sqrt.(x1)

                @. ((1.0 - AREPO_SFM.ωH2) / x1^2) + (AREPO_SFM.ωH2 * exp(-0.00085 * sqx1) / sqx1)

            end

        end,
        integrate_func = (d, c) -> nothing,
        valid_types    = (:gas, :stellar),
    ),

    # Diameter of a sphere with the same density and mass as the gas cell [L]
    :equivalent_size => (
        label_fmt      = l -> L"d_\text{eq}^{%$(l)}",
        unit           = u"pc",
        exp_factor     = 0,
        request        = c -> (
            if c == :gas
                request = Dict(c => ["RHOC", "MASS"])
            else
                request = Dict(c => ["RHOC", "GMAS"])
            end
        ),
        scatter_func   = (d, c) -> let

            if c == :gas
                m = d[c]["MASS"]
            else
                m = d[c]["GMAS"]
            end

            ρ  = d[c]["RHOC"] .* u"mp"

            V = @. m / ρ

            @. 2.0 * cbrt(V / (4π / 3))

        end,
        integrate_func = (d, c) -> nothing,
        valid_types    = (:gas, :stellar),
    ),

)

"""
    haloQuantity(magnitude::Symbol, halo_idx::Int)::BaseQuantity

Construct a halo [`BaseQuantity`](@ref).

# Arguments

  - `magnitude::Symbol`: Target halo quantity. The options are the keys of [`HALO_QUANTITIES`](@ref).
  - `halo_idx::Int`: Index of the target halo (FoF group). Starts at 1.

# Returns

  - The selected [`BaseQuantity`](@ref) for the target halo.
"""
function haloQuantity(magnitude::Symbol, halo_idx::Int)::BaseQuantity

    (
        haskey(HALO_QUANTITIES, magnitude) ||
        throw(ArgumentError("haloQuantity: I don't recognize the quantity :$(magnitude). \
        It is not a key of HALO_QUANTITIES"))
    )

    id   = Symbol(magnitude, :_, halo_idx)
    info = HALO_QUANTITIES[magnitude]

    return BaseQuantity(;
        id,
        qty_label      = info.qty_label,
        request        = Dict(:group => [info.snap_key]),
        unit           = info.unit,
        exp_factor     = info.exp_factor,
        cp_type        = :group,
        scatter_func   = dd -> nothing,
        integrate_func = dd -> let

            halo_qty = dd[:group][info.snap_key]

            if isempty(halo_qty)
                NaN * info.unit
            else
                halo_qty[halo_idx]
            end

        end,
    )

end

"""
    fillRegistry!(; <keyword arguments>)::Nothing

Fill the given registry with the component-magnitude quantities defined in [`MAGNITUDES`](@ref) and [`COMPONENTS`](@ref), the abundance quantities defined in [`ABUNDANCE_QUANTITIES`](@ref), and the SFM quantities defined in [`SFM_MAGNITUDES`](@ref) and [`SFM_COMPONENTS`](@ref).

# Arguments

  - `registry::Dict{Symbol, BaseQuantity}`: The registry to be filled. By default, it is the global `QTY_REGISTRY`.
"""
function fillRegistry!(; registry::Dict{Symbol, BaseQuantity}=QTY_REGISTRY)::Nothing

    # Populate the registry with the component-magnitude quantities
    for (magnitude, trait) in MAGNITUDES

        for (component, info) in COMPONENTS

            id = Symbol(component, :_, magnitude)

            request = trait.request(component)

            isnothing(request) && continue

            registry[id] = BaseQuantity(;
                id,
                qty_label      = trait.label_fmt(info.qty_label),
                request        = trait.request(component),
                unit           = trait.unit,
                exp_factor     = trait.exp_factor,
                cp_type        = info.cp_type,
                scatter_func   = dd -> trait.scatter_func(dd, component),
                integrate_func = dd -> trait.integrate_func(dd, component),
            )

        end

    end

    # Populate the registry with the abundance quantities
    for element in keys(ELEMENT_INDEX)

        for (component, info) in ABUNDANCE_QUANTITIES

            id = Symbol(element, :_, component)

            registry[id] = BaseQuantity(;
                id,
                qty_label      = info.label_fmt(element),
                request        = info.request,
                unit           = info.unit,
                exp_factor     = info.exp_factor,
                cp_type        = info.cp_type,
                scatter_func   = dd -> info.scatter_func(dd, element),
                integrate_func = dd -> info.integrate_func(dd, element),
            )

        end

    end

    # Populate the registry with the SFM quantities
    for (magnitude, trait) in SFM_MAGNITUDES

        for (component, info) in SFM_COMPONENTS

            id = Symbol(component, :_, magnitude)

            (info.cp_type ∈ trait.valid_types) || continue

            registry[id] = BaseQuantity(;
                id,
                qty_label      = trait.label_fmt(info.qty_label),
                request        = trait.request(info.cp_type),
                unit           = trait.unit,
                exp_factor     = trait.exp_factor,
                cp_type        = info.cp_type,
                scatter_func   = dd -> trait.scatter_func(dd, info.cp_type),
                integrate_func = dd -> trait.integrate_func(dd, info.cp_type),
            )

        end

    end

    return nothing

end

