# Run constrained SSP projections (emissions driven) for completeness.

import os
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.optimize
import xarray as xr
from dotenv import load_dotenv
from fair import FAIR
from fair.interface import fill, initialise
from fair.io import read_properties
from fair.structure.units import *

load_dotenv()

DATADIR = Path(os.getenv("DATADIR"))
DATAIN = DATADIR.joinpath("data_input")
DATAOUT = DATADIR.joinpath("data_output", "stochastic")

print("Running RCP scenarios...")

# Use a very crude curve fit to estimate CO2, CH4, N2O in 1765
obs_co2 = np.array([278.3, 285.5, 296.4])
obs_years = [1750, 1850, 1900]


def fit(x, a, b, c):
    return a + b * np.exp(c * (x - 1750))


p, _ = scipy.optimize.curve_fit(fit, obs_years, obs_co2, p0=[278.3, 1.45, 0.0171])
offset1765co2 = fit(1765, p[0], p[1], p[2]) - fit(1750, p[0], p[1], p[2])

obs_ch4 = np.array([729.2, 807.6, 925.1])
obs_years = [1750, 1850, 1900]


def fit(x, a, b, c):
    return a + b * np.exp(c * (x - 1750))


p, _ = scipy.optimize.curve_fit(fit, obs_years, obs_ch4, p0=[278.3, 1.45, 0.0171])
baseline1765ch4 = fit(1765, p[0], p[1], p[2])

obs_n2o = np.array([270.1, 272.1, 278.9])
obs_years = [1750, 1850, 1900]


def fit(x, a, b, c):
    return a + b * np.exp(c * (x - 1750))


p, _ = scipy.optimize.curve_fit(fit, obs_years, obs_n2o, p0=[278.3, 1.45, 0.0171])
baseline1765n2o = fit(1765, p[0], p[1], p[2])

avinox_file = DATAIN.joinpath("aviNOx_fraction.csv")
df_avinox = pd.read_csv(avinox_file, skiprows=4, index_col=0)


scenarios = [
    "rcp26",
    "rcp45",
    "rcp60",
    "rcp85",
]

species = [
    "CO2 FFI",
    "CO2 AFOLU",
    "CO2",
    "CH4",
    "N2O",
    "Sulfur",
    "BC",
    "OC",
    "NH3",
    "NOx",
    "VOC",
    "CO",
    "CFC-11",
    "CFC-12",
    "CFC-113",
    "CFC-114",
    "CFC-115",
    "HCFC-22",
    "HCFC-141b",
    "HCFC-142b",
    "CCl4",
    "CH3Cl",
    "CH3CCl3",
    "CH3Br",
    "Halon-1211",
    "Halon-1202",
    "Halon-1301",
    "Halon-2402",
    "CF4",
    "C2F6",
    "C6F14",
    "SF6",
    "HFC-125",
    "HFC-134a",
    "HFC-143a",
    "HFC-227ea",
    "HFC-23",
    "HFC-245fa",
    "HFC-32",
    "HFC-4310mee",
    "Solar",
    "Volcanic",
    "Aerosol-radiation interactions",
    "Aerosol-cloud interactions",
    "Ozone",
    "Light absorbing particles on snow and ice",
    "Land use",
    "Stratospheric water vapour",
    "Equivalent effective stratospheric chlorine",
    "Contrails",
    "NOx aviation",
]

not_emitted = [
    "Solar",
    "Volcanic",
    "Aerosol-radiation interactions",
    "Aerosol-cloud interactions",
    "Ozone",
    "Light absorbing particles on snow and ice",
    "Land use",
    "Stratospheric water vapour",
    "Equivalent effective stratospheric chlorine",
    "Contrails",
    "NOx aviation",
]

df_solar = pd.read_csv(DATAIN.joinpath("solar_erf_timebounds.csv"), index_col="year")
df_volcanic = pd.read_csv(DATAIN.joinpath("volcanic_ERF_monthly_-950001-201912.csv"))

# run from 1765 : drop 15 years
solar_forcing = np.zeros(537)
volcanic_forcing = np.zeros(537)

trend_shape = np.ones(537)
trend_shape[:256] = np.linspace(0, 1, 271)[15:]

for i, year in enumerate(np.arange(1765, 2021)):
    volcanic_forcing[i] = np.mean(
        df_volcanic.loc[
            ((year - 1) <= df_volcanic["year"]) & (df_volcanic["year"] < year)
        ].erf
    )
volcanic_forcing[256:266] = np.linspace(1, 0, 10) * volcanic_forcing[255]
solar_forcing[:536] = df_solar["erf"].loc[1765:2300].values

df_configs = pd.read_csv(
    DATAIN.joinpath("calibrated_constrained_parameters.csv"), index_col=0
)
configs = np.array(list(df_configs.index))

f = FAIR(ch4_method="Thornhill2021")
f.define_time(1765, 2301, 1)
f.define_scenarios(scenarios)
f.define_configs(configs)
species, properties = read_properties(species=species)
f.define_species(species, properties)
f.allocate()

species_to_rcmip = {specie: specie.replace("-", "") for specie in species}
species_to_rcmip["CO2 FFI"] = "CO2|MAGICC Fossil and Industrial"
species_to_rcmip["CO2 AFOLU"] = "CO2|MAGICC AFOLU"
df_rcmip = pd.read_csv(DATAIN.joinpath("rcmip-emissions-annual-means-v5-1-0.csv"))
for iscen, scenario in enumerate(scenarios):
    for specie, specie_rcmip_name in species_to_rcmip.items():
        # Grab raw emissions from dataframe
        if specie in not_emitted:
            continue
        emis = (
            df_rcmip.loc[
                (df_rcmip["Scenario"] == scenario)
                & (df_rcmip["Variable"].str.endswith("|" + specie_rcmip_name))
                & (df_rcmip["Region"] == "World"),
                "1765":"2300",
            ]
            .interpolate(axis=1)
            .values.squeeze()
        )

        # throw error if data missing
        if emis.shape[0] == 0:
            raise ValueError(
                f"I can't find a value for scenario={scenario}, variable "
                f"name ending with {specie_rcmip_name} in the RCMIP "
                f"emissions database."
            )

        # Parse and possibly convert unit in input file to what FaIR wants
        unit = df_rcmip.loc[
            (df_rcmip["Scenario"] == scenario)
            & (df_rcmip["Variable"].str.endswith("|" + specie_rcmip_name))
            & (df_rcmip["Region"] == "World"),
            "Unit",
        ].values[0]
        emis = emis * (
            prefix_convert[unit.split()[0]][desired_emissions_units[specie].split()[0]]
            * compound_convert[unit.split()[1].split("/")[0]][
                desired_emissions_units[specie].split()[1].split("/")[0]
            ]
            * time_convert[unit.split()[1].split("/")[1]][
                desired_emissions_units[specie].split()[1].split("/")[1]
            ]
        )

        # fill FaIR xarray
        fill(f.emissions, emis[:, None], specie=specie, scenario=scenario)

    # Aviation NOx
    f.emissions.loc[dict(specie="NOx aviation", scenario=scenario)] = (
        f.emissions.loc[dict(specie="NOx", scenario=scenario)]
        * df_avinox.values[:536, iscen : iscen + 1]
    )

    # Replace 1765-1850 in SLCFs with ramp up from Skeie et al. 2011 piControl values
    f.emissions[:86, iscen, 0, 5] = np.linspace(2, f.emissions[85, iscen, 0, 5], 86)
    f.emissions[:86, iscen, 0, 6] = np.linspace(1.2, f.emissions[85, iscen, 0, 6], 86)
    f.emissions[:86, iscen, 0, 7] = np.linspace(10, f.emissions[85, iscen, 0, 7], 86)
    f.emissions[:86, iscen, 0, 8] = np.linspace(4, f.emissions[85, iscen, 0, 8], 86)
    f.emissions[:86, iscen, 0, 9] = np.linspace(
        46 / 14 * 2, f.emissions[85, iscen, 0, 9], 86
    )
    f.emissions[:86, iscen, 0, 10] = np.linspace(10, f.emissions[85, iscen, 0, 10], 86)
    f.emissions[:86, iscen, 0, 11] = np.linspace(174, f.emissions[85, iscen, 0, 11], 86)

# solar and volcanic forcing
fill(
    f.forcing,
    volcanic_forcing[:, None, None] * df_configs["scale Volcanic"].values.squeeze(),
    specie="Volcanic",
)
fill(
    f.forcing,
    solar_forcing[:, None, None] * df_configs["solar_amplitude"].values.squeeze()
    + trend_shape[:, None, None] * df_configs["solar_trend"].values.squeeze(),
    specie="Solar",
)

# climate response
fill(f.climate_configs["ocean_heat_capacity"], df_configs.loc[:, "c1":"c3"].values)
fill(
    f.climate_configs["ocean_heat_transfer"],
    df_configs.loc[:, "kappa1":"kappa3"].values,
)
fill(f.climate_configs["deep_ocean_efficacy"], df_configs["epsilon"].values.squeeze())
fill(f.climate_configs["gamma_autocorrelation"], df_configs["gamma"].values.squeeze())
fill(f.climate_configs["sigma_eta"], df_configs["sigma_eta"].values.squeeze())
fill(f.climate_configs["sigma_xi"], df_configs["sigma_xi"].values.squeeze())
fill(f.climate_configs["seed"], df_configs["seed"])
fill(f.climate_configs["stochastic_run"], True)
fill(f.climate_configs["use_seed"], True)
fill(f.climate_configs["forcing_4co2"], df_configs["F_4xCO2"])

# species level
f.fill_species_configs()

# carbon cycle
fill(f.species_configs["iirf_0"], df_configs["r0"].values.squeeze(), specie="CO2")
fill(
    f.species_configs["iirf_airborne"],
    df_configs["rA"].values.squeeze(),
    specie="CO2",
)
fill(
    f.species_configs["iirf_uptake"],
    df_configs["rU"].values.squeeze(),
    specie="CO2",
)
fill(
    f.species_configs["iirf_temperature"],
    df_configs["rT"].values.squeeze(),
    specie="CO2",
)

# aerosol indirect
fill(f.species_configs["aci_scale"], df_configs["beta"].values.squeeze())
fill(
    f.species_configs["aci_shape"],
    df_configs["shape Sulfur"].values.squeeze(),
    specie="Sulfur",
)
fill(
    f.species_configs["aci_shape"],
    df_configs["shape BC"].values.squeeze(),
    specie="BC",
)
fill(
    f.species_configs["aci_shape"],
    df_configs["shape OC"].values.squeeze(),
    specie="OC",
)

# methane lifetime baseline - should be imported from calibration
fill(f.species_configs["unperturbed_lifetime"], 10.11702748, specie="CH4")

# emissions adjustments for N2O and CH4
fill(f.species_configs["baseline_emissions"], 19.019783117809567, specie="CH4")
fill(f.species_configs["baseline_emissions"], 0.08602230754, specie="N2O")

# aerosol direct
for specie in [
    "BC",
    "CH4",
    "N2O",
    "NH3",
    "NOx",
    "OC",
    "Sulfur",
    "VOC",
    "Equivalent effective stratospheric chlorine",
]:
    fill(
        f.species_configs["erfari_radiative_efficiency"],
        df_configs[f"ari {specie}"],
        specie=specie,
    )

# forcing scaling
for specie in [
    "CO2",
    "CH4",
    "N2O",
    "Stratospheric water vapour",
    "Contrails",
    "Light absorbing particles on snow and ice",
    "Land use",
]:
    fill(
        f.species_configs["forcing_scale"],
        df_configs[f"scale {specie}"].values.squeeze(),
        specie=specie,
    )

for specie in [
    "CFC-11",
    "CFC-12",
    "CFC-113",
    "CFC-114",
    "CFC-115",
    "HCFC-22",
    "HCFC-141b",
    "HCFC-142b",
    "CCl4",
    "CH3Cl",
    "CH3CCl3",
    "CH3Br",
    "Halon-1202",
    "Halon-1211",
    "Halon-1301",
    "Halon-2402",
    "CF4",
    "C2F6",
    "C6F14",
    "SF6",
    "HFC-125",
    "HFC-134a",
    "HFC-143a",
    "HFC-227ea",
    "HFC-23",
    "HFC-245fa",
    "HFC-32",
    "HFC-4310mee",
]:
    # 'CFC-11', 'CFC-12', 'CFC-113', 'CFC-114', 'CFC-115',
    # 'HCFC-22', 'HCFC-141b', 'HCFC-142b',
    # 'CCl4', 'CH3Cl', 'CH3CCl3', 'CH3Br',
    # 'Halon-1211', 'Halon-1202', 'Halon-1301', 'Halon-2402',
    # 'CF4', 'C2F6', 'C6F14',
    # 'SF6',
    # 'HFC-125', 'HFC-134a', 'HFC-143a', 'HFC-227ea', 'HFC-23', 'HFC-245fa', 'HFC-32',
    # 'HFC-4310mee',
    fill(
        f.species_configs["forcing_scale"],
        df_configs["scale minorGHG"].values.squeeze(),
        specie=specie,
    )

# ozone
for specie in [
    "CH4",
    "N2O",
    "Equivalent effective stratospheric chlorine",
    "CO",
    "VOC",
    "NOx",
]:
    fill(
        f.species_configs["ozone_radiative_efficiency"],
        df_configs[f"o3 {specie}"],
        specie=specie,
    )

# tune down volcanic efficacy
fill(f.species_configs["forcing_efficacy"], 0.6, specie="Volcanic")

# initial condition of CO2 concentration (but not baseline for forcing calculations)
fill(
    f.species_configs["baseline_concentration"],
    df_configs["co2_concentration_1750"].values.squeeze(),
    specie="CO2",
)

# initial condition of CO2 concentration (but not baseline for forcing calculations)
fill(
    f.species_configs["baseline_concentration"],
    offset1765co2 + df_configs.loc[configs, "co2_concentration_1750"].values.squeeze(),
    specie="CO2",
)

# initial condition of other species
fill(f.species_configs["baseline_concentration"], baseline1765ch4, specie="CH4")
fill(f.species_configs["baseline_concentration"], baseline1765n2o, specie="N2O")
fill(f.species_configs["baseline_emissions"], 2, specie="Sulfur")
fill(f.species_configs["baseline_emissions"], 174, specie="CO")
fill(f.species_configs["baseline_emissions"], 10, specie="VOC")
fill(f.species_configs["baseline_emissions"], 4, specie="NH3")
fill(f.species_configs["baseline_emissions"], 2 * 46 / 14, specie="NOx")
fill(f.species_configs["baseline_emissions"], 1.2, specie="BC")
fill(f.species_configs["baseline_emissions"], 10, specie="OC")

# initial conditions
initialise(f.concentration, f.species_configs["baseline_concentration"])
initialise(f.forcing, 0)
initialise(f.temperature, 0)
initialise(f.cumulative_emissions, 0)
initialise(f.airborne_emissions, 0)

f.run()


average_51yr = np.ones(52)
average_51yr[0] = 0.5
average_51yr[-1] = 0.5

# at this point dump out some batch output
temp_out = f.temperature[:, :, :, 0].data
ohc_out = f.ocean_heat_content_change[:, :, :].data
erf_out = f.forcing_sum[:, :, :].data
co2_out = f.concentration[:, :, :, 2].data
ch4_out = f.concentration[:, :, :, 3].data
n2o_out = f.concentration[:, :, 0, 4].data

for iscen, scenario in enumerate(scenarios):
    ds = xr.Dataset(
        {
            "temperature": (
                ["year", "run"],
                temp_out[:, iscen, :]
                - np.average(temp_out[85:137, iscen, :], weights=average_51yr, axis=0),
            ),
            "effective_radiative_forcing": (["year", "run"], erf_out[:, iscen, :]),
            "ocean_heat_content_change": (["year", "run"], ohc_out[:, iscen, :]),
            "co2_concentration": (["year", "run"], co2_out[:, iscen, :]),
            "ch4_concentration": (["year", "run"], ch4_out[:, iscen, :]),
            "n2o_concentration": (["year"], n2o_out[:, iscen]),
        },
        coords={"year": (np.arange(1765, 2301.5)), "run": configs},
    )
    ds.to_netcdf(
        DATAOUT.joinpath(f"{scenario}.nc"),
        encoding={
            "temperature": {"dtype": "float32"},
            "effective_radiative_forcing": {"dtype": "float32"},
            "ocean_heat_content_change": {"dtype": "float32"},
            "co2_concentration": {"dtype": "float32"},
            "ch4_concentration": {"dtype": "float32"},
            "n2o_concentration": {"dtype": "float32"},
        },
    )
    ds.close()
