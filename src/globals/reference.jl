####################################################################################################
# Reference values from the literature
####################################################################################################

"""
Internal unit of length used in [IllustrisTNG](https://www.tng-project.org/data/docs/specifications/), equal to ``1.0  \\, \\mathrm{kpc}``.
"""
const ILLUSTRIS_L_UNIT = 3.085678e21u"cm"

"""
Internal unit of mass used in [IllustrisTNG](https://www.tng-project.org/data/docs/specifications/), equal to ``10^{10} \\, \\mathrm{M_\\odot}``.
"""
const ILLUSTRIS_M_UNIT = 1.989e43u"g"

"""
Internal unit of velocity used in [IllustrisTNG](https://www.tng-project.org/data/docs/specifications/), equal to ``1.0 \\, \\mathrm{km \\, s^{-1}}``.
"""
const ILLUSTRIS_V_UNIT = 1.0e5u"cm * s^-1"

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
const KS98_SFR_RANGE = exp10.([-3.55, 2.98]) .* u"Msun * yr^-1 * kpc^-2"

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

"""
Spatial resolution used in Bigiel et al. (2008).

# References

F. Bigiel et al. (2008). *THE STAR FORMATION LAW IN NEARBY GALAXIES ON SUB-KPC SCALES*. The Astrophysical Journal, **136(6)**, 2846. [doi:10.1088/0004-6256/136/6/2846](https://doi.org/10.1088/0004-6256/136/6/2846)
"""
const BIGIEL_PX_SIZE = 750.0u"pc"

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
const BIGIEL2008_SFR_RANGE = exp10.([-2.99, -0.33]) .* u"Msun * yr^-1 * kpc^-2"

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

"""
Spatial resolution used in Sun et al. (2023).

# References

J. Sun et al. (2023). *Star Formation Laws and Efficiencies across 80 Nearby Galaxies*. The Astrophysical Journal Letters, **945(2)**, L19. [doi:10.3847/2041-8213/acbd9c](https://doi.org/10.3847/2041-8213/acbd9c)
"""
const SUN_PX_SIZE = 1.5u"kpc"

"""
Path to Table A1 from Sun et al. (2023).

# References

J. Sun et al. (2023). *Star Formation Laws and Efficiencies across 80 Nearby Galaxies*. The Astrophysical Journal Letters, **945(2)**, L19. [doi:10.3847/2041-8213/acbd9c](https://doi.org/10.3847/2041-8213/acbd9c)
"""
const SUN2023_TABLE = joinpath(@__DIR__, "../../experimental_data/sun_2023.txt")

"""
Path to Table 4 (corrected) from de los Reyes et al. (2019).

# References

M. A. C. de los Reyes et al. (2019). *Revisiting the Integrated Star Formation Law. I. Non-starbursting Galaxies*. The Astrophysical Journal, **872(1)**, 16. [doi:10.3847/1538-4357/aafa82)](https://doi.org/10.3847/1538-4357/aafa82)

M. A. C. de los Reyes et al. (2019). *Erratum: “Revisiting the Integrated Star Formation Law. I. Non-starbursting Galaxies” (2019 ApJ, 872, 16)*. The Astrophysical Journal, **878(1)**, 74. [doi:10.3847/1538-4357/ab22af)](https://doi.org/10.3847/1538-4357/ab22af)
"""
const DELOSREYES2019_TABLE = joinpath(@__DIR__, "../../experimental_data/de_los_reyes_2019.txt")

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
