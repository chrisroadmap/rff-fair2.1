def run_fair(sample):
    import os
    from pathlib import Path

    import numpy as np
    import pandas as pd
    import xarray as xr
    from dotenv import load_dotenv
    from fair import FAIR
    from fair.interface import fill, initialise
    from fair.io import read_properties

    # Get environment variables
    load_dotenv()

    # Make data directory
    DATADIR = Path(os.getenv("DATADIR"))
    DATAIN = DATADIR.joinpath("data_input")
    DATAPROCESSED = DATADIR.joinpath("data_processed", "infilled_extended")
    DATAOUT = DATADIR.joinpath("data_output", "stochastic")

    df_solar = pd.read_csv(
        DATAIN.joinpath("solar_erf_timebounds.csv"), index_col="year"
    )
    df_volcanic = pd.read_csv(
        DATAIN.joinpath("volcanic_ERF_monthly_-950001-201912.csv")
    )

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

    rcmip_df = pd.read_csv(DATAIN.joinpath("rcmip-emissions-annual-means-v5-1-0.csv"))

    variable_mappings = {
        "AR6 climate diagnostics|Emissions|BC": "Emissions|BC",
        "AR6 climate diagnostics|Emissions|CCl4": "Emissions|Montreal Gases|CCl4",
        "AR6 climate diagnostics|Emissions|CFC11": "Emissions|Montreal Gases|CFC|CFC11",
        "AR6 climate diagnostics|Emissions|CFC113": "Emissions|Montreal Gases|CFC|CFC113",
        "AR6 climate diagnostics|Emissions|CFC114": "Emissions|Montreal Gases|CFC|CFC114",
        "AR6 climate diagnostics|Emissions|CFC115": "Emissions|Montreal Gases|CFC|CFC115",
        "AR6 climate diagnostics|Emissions|CFC12": "Emissions|Montreal Gases|CFC|CFC12",
        "AR6 climate diagnostics|Emissions|CH2Cl2": "Emissions|Montreal Gases|CH2Cl2",
        "AR6 climate diagnostics|Emissions|CH3Br": "Emissions|Montreal Gases|CH3Br",
        "AR6 climate diagnostics|Emissions|CH3CCl3": "Emissions|Montreal Gases|CH3CCl3",
        "AR6 climate diagnostics|Emissions|CH3Cl": "Emissions|Montreal Gases|CH3Cl",
        "AR6 climate diagnostics|Emissions|CH4": "Emissions|CH4",
        "AR6 climate diagnostics|Emissions|CHCl3": "Emissions|Montreal Gases|CHCl3",
        "AR6 climate diagnostics|Emissions|CO": "Emissions|CO",
        "AR6 climate diagnostics|Emissions|CO2|AFOLU": "Emissions|CO2|MAGICC AFOLU",
        "AR6 climate diagnostics|Emissions|CO2|Energy and Industrial Processes": "Emissions|CO2|MAGICC Fossil and Industrial",
        "AR6 climate diagnostics|Emissions|HCFC141b": "Emissions|Montreal Gases|HCFC141b",
        "AR6 climate diagnostics|Emissions|HCFC142b": "Emissions|Montreal Gases|HCFC142b",
        "AR6 climate diagnostics|Emissions|HCFC22": "Emissions|Montreal Gases|HCFC22",
        "AR6 climate diagnostics|Emissions|HFC|HFC125": "Emissions|F-Gases|HFC|HFC125",
        "AR6 climate diagnostics|Emissions|HFC|HFC134a": "Emissions|F-Gases|HFC|HFC134a",
        "AR6 climate diagnostics|Emissions|HFC|HFC143a": "Emissions|F-Gases|HFC|HFC143a",
        "AR6 climate diagnostics|Emissions|HFC|HFC152a": "Emissions|F-Gases|HFC|HFC152a",
        "AR6 climate diagnostics|Emissions|HFC|HFC227ea": "Emissions|F-Gases|HFC|HFC227ea",
        "AR6 climate diagnostics|Emissions|HFC|HFC23": "Emissions|F-Gases|HFC|HFC23",
        "AR6 climate diagnostics|Emissions|HFC|HFC236fa": "Emissions|F-Gases|HFC|HFC236fa",
        "AR6 climate diagnostics|Emissions|HFC|HFC245fa": "Emissions|F-Gases|HFC|HFC245fa",
        "AR6 climate diagnostics|Emissions|HFC|HFC32": "Emissions|F-Gases|HFC|HFC32",
        "AR6 climate diagnostics|Emissions|HFC|HFC365mfc": "Emissions|F-Gases|HFC|HFC365mfc",
        "AR6 climate diagnostics|Emissions|HFC|HFC43-10": "Emissions|F-Gases|HFC|HFC4310mee",
        "AR6 climate diagnostics|Emissions|Halon1202": "Emissions|Montreal Gases|Halon1202",
        "AR6 climate diagnostics|Emissions|Halon1211": "Emissions|Montreal Gases|Halon1211",
        "AR6 climate diagnostics|Emissions|Halon1301": "Emissions|Montreal Gases|Halon1301",
        "AR6 climate diagnostics|Emissions|Halon2402": "Emissions|Montreal Gases|Halon2402",
        "AR6 climate diagnostics|Emissions|N2O": "Emissions|N2O",
        "AR6 climate diagnostics|Emissions|NF3": "Emissions|F-Gases|NF3",
        "AR6 climate diagnostics|Emissions|NH3": "Emissions|NH3",
        "AR6 climate diagnostics|Emissions|NOx": "Emissions|NOx",
        "AR6 climate diagnostics|Emissions|NOx|Aviation": "Emissions|NOx|MAGICC Fossil and Industrial|Aircraft",
        "AR6 climate diagnostics|Emissions|OC": "Emissions|OC",
        "AR6 climate diagnostics|Emissions|PFC|C2F6": "Emissions|F-Gases|PFC|C2F6",
        "AR6 climate diagnostics|Emissions|PFC|C3F8": "Emissions|F-Gases|PFC|C3F8",
        "AR6 climate diagnostics|Emissions|PFC|C4F10": "Emissions|F-Gases|PFC|C4F10",
        "AR6 climate diagnostics|Emissions|PFC|C5F12": "Emissions|F-Gases|PFC|C5F12",
        "AR6 climate diagnostics|Emissions|PFC|C6F14": "Emissions|F-Gases|PFC|C6F14",
        "AR6 climate diagnostics|Emissions|PFC|C7F16": "Emissions|F-Gases|PFC|C7F16",
        "AR6 climate diagnostics|Emissions|PFC|C8F18": "Emissions|F-Gases|PFC|C8F18",
        "AR6 climate diagnostics|Emissions|PFC|CF4": "Emissions|F-Gases|PFC|CF4",
        "AR6 climate diagnostics|Emissions|PFC|cC4F8": "Emissions|F-Gases|PFC|cC4F8",
        "AR6 climate diagnostics|Emissions|SF6": "Emissions|F-Gases|SF6",
        "AR6 climate diagnostics|Emissions|SO2F2": "Emissions|F-Gases|SO2F2",
        "AR6 climate diagnostics|Emissions|Sulfur": "Emissions|Sulfur",
        "AR6 climate diagnostics|Emissions|VOC": "Emissions|VOC",
    }

    hist_df = rcmip_df.loc[
        (rcmip_df["Scenario"] == "historical")
        & (rcmip_df["Region"] == "World")
        & (rcmip_df["Variable"].isin(variable_mappings.values()))
    ]

    ucfuture = {var: 1 for var in variable_mappings.keys()}
    ucfuture[
        "AR6 climate diagnostics|Emissions|CO2|Energy and Industrial Processes"
    ] = 0.001
    ucfuture["AR6 climate diagnostics|Emissions|CO2|AFOLU"] = 0.001

    ucpast = {var: 1 for var in variable_mappings.values()}
    ucpast["Emissions|CO2|MAGICC Fossil and Industrial"] = 0.001
    ucpast["Emissions|CO2|MAGICC AFOLU"] = 0.001
    ucpast["Emissions|N2O"] = 0.001

    species_mapping = {
        "AR6 climate diagnostics|Emissions|BC": "BC",
        "AR6 climate diagnostics|Emissions|CCl4": "CCl4",
        "AR6 climate diagnostics|Emissions|CFC11": "CFC-11",
        "AR6 climate diagnostics|Emissions|CFC113": "CFC-113",
        "AR6 climate diagnostics|Emissions|CFC114": "CFC-114",
        "AR6 climate diagnostics|Emissions|CFC115": "CFC-115",
        "AR6 climate diagnostics|Emissions|CFC12": "CFC-12",
        "AR6 climate diagnostics|Emissions|CH2Cl2": "CH2Cl2",
        "AR6 climate diagnostics|Emissions|CH3Br": "CH3Br",
        "AR6 climate diagnostics|Emissions|CH3CCl3": "CH3CCl3",
        "AR6 climate diagnostics|Emissions|CH3Cl": "CH3Cl",
        "AR6 climate diagnostics|Emissions|CH4": "CH4",
        "AR6 climate diagnostics|Emissions|CHCl3": "CHCl3",
        "AR6 climate diagnostics|Emissions|CO": "CO",
        "AR6 climate diagnostics|Emissions|CO2|AFOLU": "CO2 AFOLU",
        "AR6 climate diagnostics|Emissions|CO2|Energy and Industrial Processes": "CO2 FFI",
        "AR6 climate diagnostics|Emissions|HCFC141b": "HCFC-141b",
        "AR6 climate diagnostics|Emissions|HCFC142b": "HCFC-142b",
        "AR6 climate diagnostics|Emissions|HCFC22": "HCFC-22",
        "AR6 climate diagnostics|Emissions|HFC|HFC125": "HFC-125",
        "AR6 climate diagnostics|Emissions|HFC|HFC134a": "HFC-134a",
        "AR6 climate diagnostics|Emissions|HFC|HFC143a": "HFC-143a",
        "AR6 climate diagnostics|Emissions|HFC|HFC152a": "HFC-152a",
        "AR6 climate diagnostics|Emissions|HFC|HFC227ea": "HFC-227ea",
        "AR6 climate diagnostics|Emissions|HFC|HFC23": "HFC-23",
        "AR6 climate diagnostics|Emissions|HFC|HFC236fa": "HFC-236fa",
        "AR6 climate diagnostics|Emissions|HFC|HFC245fa": "HFC-245fa",
        "AR6 climate diagnostics|Emissions|HFC|HFC32": "HFC-32",
        "AR6 climate diagnostics|Emissions|HFC|HFC365mfc": "HFC-365mfc",
        "AR6 climate diagnostics|Emissions|HFC|HFC43-10": "HFC-4310mee",
        "AR6 climate diagnostics|Emissions|Halon1202": "Halon-1202",
        "AR6 climate diagnostics|Emissions|Halon1211": "Halon-1211",
        "AR6 climate diagnostics|Emissions|Halon1301": "Halon-1301",
        "AR6 climate diagnostics|Emissions|Halon2402": "Halon-2402",
        "AR6 climate diagnostics|Emissions|N2O": "N2O",
        "AR6 climate diagnostics|Emissions|NF3": "NF3",
        "AR6 climate diagnostics|Emissions|NH3": "NH3",
        "AR6 climate diagnostics|Emissions|NOx": "NOx",
        "AR6 climate diagnostics|Emissions|NOx|Aviation": "NOx aviation",
        "AR6 climate diagnostics|Emissions|OC": "OC",
        "AR6 climate diagnostics|Emissions|PFC|C2F6": "C2F6",
        "AR6 climate diagnostics|Emissions|PFC|C3F8": "C3F8",
        "AR6 climate diagnostics|Emissions|PFC|C4F10": "C4F10",
        "AR6 climate diagnostics|Emissions|PFC|C5F12": "C5F12",
        "AR6 climate diagnostics|Emissions|PFC|C6F14": "C6F14",
        "AR6 climate diagnostics|Emissions|PFC|C7F16": "C7F16",
        "AR6 climate diagnostics|Emissions|PFC|C8F18": "C8F18",
        "AR6 climate diagnostics|Emissions|PFC|CF4": "CF4",
        "AR6 climate diagnostics|Emissions|PFC|cC4F8": "c-C4F8",
        "AR6 climate diagnostics|Emissions|SF6": "SF6",
        "AR6 climate diagnostics|Emissions|SO2F2": "SO2F2",
        "AR6 climate diagnostics|Emissions|Sulfur": "Sulfur",
        "AR6 climate diagnostics|Emissions|VOC": "VOC",
    }

    average_51yr = np.ones(52)
    average_51yr[0] = 0.5
    average_51yr[-1] = 0.5

    future_df = pd.read_csv(DATAPROCESSED.joinpath("emissions%05d.csv" % sample))
    emissions = np.ones((551, 53)) * np.nan

    for ivar, (varnamefuture, varnamepast) in enumerate(variable_mappings.items()):
        emissions[:265, ivar] = (
            hist_df.loc[(hist_df["Variable"] == varnamepast), "1750":"2014"]
            * ucpast[varnamepast]
        )
        emissions[265:, ivar] = (
            future_df.loc[(future_df["variable"] == varnamefuture), "2015":]
            * ucfuture[varnamefuture]
        )

    species, properties = read_properties()
    f = FAIR(ch4_method="Thornhill2021")
    f.define_time(1750, 2301, 1)
    f.define_scenarios([sample])
    f.define_configs(configs)
    f.define_species(species, properties)
    f.allocate()

    emissions_df = pd.DataFrame(
        emissions, columns=variable_mappings.keys(), index=np.arange(1750, 2301)
    )

    for varnamefut, varnamefair in species_mapping.items():
        f.emissions.loc[dict(specie=varnamefair)] = emissions_df[varnamefut].values[
            :, None, None
        ]

    fill(
        f.forcing,
        volcanic_forcing[:, None, None]
        * df_configs.loc[configs, "scale Volcanic"].values.squeeze(),
        specie="Volcanic",
    )
    fill(
        f.forcing,
        solar_forcing[:, None, None]
        * df_configs.loc[configs, "solar_amplitude"].values.squeeze()
        + trend_shape[:, None, None]
        * df_configs.loc[configs, "solar_trend"].values.squeeze(),
        specie="Solar",
    )

    # climate response
    fill(f.climate_configs["ocean_heat_capacity"], df_configs.loc[:, "c1":"c3"].values)
    fill(
        f.climate_configs["ocean_heat_transfer"],
        df_configs.loc[:, "kappa1":"kappa3"].values,
    )
    fill(
        f.climate_configs["deep_ocean_efficacy"], df_configs["epsilon"].values.squeeze()
    )
    fill(
        f.climate_configs["gamma_autocorrelation"], df_configs["gamma"].values.squeeze()
    )
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

    f.run(progress=False)

    # at this point dump out some batch output
    temp_out = f.temperature[:, 0, :, 0].data
    ohc_out = f.ocean_heat_content_change[:, 0, :].data
    erf_out = f.forcing_sum[:, 0, :].data
    co2_out = f.concentration[:, 0, :, 2].data
    ch4_out = f.concentration[:, 0, :, 3].data
    n2o_out = f.concentration[:, 0, 0, 4].data

    idx0 = 0 if sample==1 else 270
    year0 = 1750 if sample==1 else 2020

    ds = xr.Dataset(
        {
            "temperature": (
                ["year", "run"],
                temp_out[idx0:]
                - np.average(temp_out[100:152, :], weights=average_51yr, axis=0),
            ),
            "effective_radiative_forcing": (["year", "run"], erf_out[idx0:]),
            "ocean_heat_content_change": (["year", "run"], ohc_out[idx0:]),
            "co2_concentration": (["year", "run"], co2_out[idx0:]),
            "ch4_concentration": (["year", "run"], ch4_out[idx0:]),
            "n2o_concentration": (["year"], n2o_out[idx0:]),
        },
        coords={"year": (np.arange(year0, 2301.5)), "run": configs},
    )
    ds.to_netcdf(
        DATAOUT.joinpath("run%05d.nc" % sample)
        encoding={
            "temperature": {"dtype": "float32"},
            "effective_radiative_forcing": {"dtype": "float32"},
            "ocean_heat_content_change": {"dtype": "float32"},
            "co2_concentration": {"dtype": "float32"},
            "ch4_concentration": {"dtype": "float32"},
            "n2o_concentration": {"dtype": "float32"},
        }
    )
    ds.close()
