# Run constrained SSP projections (emissions driven) for completeness.

import os
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
from dotenv import load_dotenv
from fair import FAIR
from fair.interface import fill, initialise
from fair.io import read_properties

load_dotenv()

DATADIR = Path(os.getenv("DATADIR"))
DATAIN = DATADIR.joinpath("data_input")
DATAOUT = DATADIR.joinpath("data_output", "stochastic")

print("Running SSP scenarios...")

scenarios = [
    "ssp119",
    "ssp126",
    "ssp245",
    "ssp370",
    "ssp434",
    "ssp460",
    "ssp534-over",
    "ssp585",
]

df_solar = pd.read_csv(DATAIN.joinpath("solar_erf_timebounds.csv"), index_col="year")
df_volcanic = pd.read_csv(DATAIN.joinpath("volcanic_ERF_monthly_-950001-201912.csv"))

solar_forcing = np.zeros(552)
volcanic_forcing = np.zeros(552)

trend_shape = np.ones(552)
trend_shape[:271] = np.linspace(0, 1, 271)

for i, year in enumerate(np.arange(1750, 2021)):
    volcanic_forcing[i] = np.mean(
        df_volcanic.loc[
            ((year - 1) <= df_volcanic["year"]) & (df_volcanic["year"] < year)
        ].erf
    )
volcanic_forcing[271:281] = np.linspace(1, 0, 10) * volcanic_forcing[270]
solar_forcing[:551] = df_solar["erf"].loc[1750:2300].values

df_configs = pd.read_csv(
    DATAIN.joinpath("calibrated_constrained_parameters.csv"), index_col=0
)
configs = np.array(list(df_configs.index))

f = FAIR(ch4_method="Thornhill2021")
f.define_time(1750, 2301, 1)
f.define_scenarios(scenarios)
f.define_configs(configs)
species, properties = read_properties()
f.define_species(species, properties)
f.allocate()

f.fill_from_rcmip()

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
    "CHCl3",
    "CH2Cl2",
    "CH3Cl",
    "CH3CCl3",
    "CH3Br",
    "Halon-1211",
    "Halon-1301",
    "Halon-2402",
    "CF4",
    "C2F6",
    "C3F8",
    "c-C4F8",
    "C4F10",
    "C5F12",
    "C6F14",
    "C7F16",
    "C8F18",
    "NF3",
    "SF6",
    "SO2F2",
    "HFC-125",
    "HFC-134a",
    "HFC-143a",
    "HFC-152a",
    "HFC-227ea",
    "HFC-23",
    "HFC-236fa",
    "HFC-245fa",
    "HFC-32",
    "HFC-365mfc",
    "HFC-4310mee",
]:
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
                - np.average(temp_out[100:152, iscen, :], weights=average_51yr, axis=0),
            ),
            "effective_radiative_forcing": (["year", "run"], erf_out[:, iscen, :]),
            "ocean_heat_content_change": (["year", "run"], ohc_out[:, iscen, :]),
            "co2_concentration": (["year", "run"], co2_out[:, iscen, :]),
            "ch4_concentration": (["year", "run"], ch4_out[:, iscen, :]),
            "n2o_concentration": (["year"], n2o_out[:, iscen]),
        },
        coords={"year": (np.arange(1750, 2301.5)), "run": configs},
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
