# This script extends the infilled RFF-SPs to 2300 using Meinshausen assumptions and
# AFOLU ratios in 2100.

# Workflow:
# 1. take in the pyam
# 2. extend each scenario in turn
# 3. save out a binary of each emissions file

import os
from pathlib import Path

import numpy as np
import pandas as pd
import pyam
from dotenv import load_dotenv
from tqdm.auto import tqdm

# Get environment variables
load_dotenv()

DATADIR = Path(os.getenv("DATADIR"))
DATAIN = DATADIR.joinpath("data_processed")
DATARFF = DATADIR.joinpath("data_processed/emissions_files")
DATAOUT = DATADIR.joinpath("data_processed/infilled_extended")

os.makedirs(DATAOUT, exist_ok=True)

# number of scenarios
RFF_SCENS = int(os.getenv("RFF_SCENS"))

df = pyam.IamDataFrame(DATAIN.joinpath("infilled_emissions_scenarios.csv"))

variables = [
    "AR6 climate diagnostics|Emissions|BC",  # 0
    "AR6 climate diagnostics|Emissions|CCl4",  # 1
    "AR6 climate diagnostics|Emissions|CFC11",
    "AR6 climate diagnostics|Emissions|CFC113",
    "AR6 climate diagnostics|Emissions|CFC114",
    "AR6 climate diagnostics|Emissions|CFC115",
    "AR6 climate diagnostics|Emissions|CFC12",
    "AR6 climate diagnostics|Emissions|CH2Cl2", # 7
    "AR6 climate diagnostics|Emissions|CH3Br", # 8
    "AR6 climate diagnostics|Emissions|CH3CCl3",
    "AR6 climate diagnostics|Emissions|CH3Cl", # 10
    "AR6 climate diagnostics|Emissions|CH4",  # 11
    "AR6 climate diagnostics|Emissions|CHCl3", # 12
    "AR6 climate diagnostics|Emissions|CO",  # 13
    "AR6 climate diagnostics|Emissions|CO2|AFOLU",  # 14
    "AR6 climate diagnostics|Emissions|CO2|Energy and Industrial Processes",  # 15
    "AR6 climate diagnostics|Emissions|HCFC141b",
    "AR6 climate diagnostics|Emissions|HCFC142b",
    "AR6 climate diagnostics|Emissions|HCFC22",
    "AR6 climate diagnostics|Emissions|HFC|HFC125",
    "AR6 climate diagnostics|Emissions|HFC|HFC134a",
    "AR6 climate diagnostics|Emissions|HFC|HFC143a",
    "AR6 climate diagnostics|Emissions|HFC|HFC152a",
    "AR6 climate diagnostics|Emissions|HFC|HFC227ea",
    "AR6 climate diagnostics|Emissions|HFC|HFC23",
    "AR6 climate diagnostics|Emissions|HFC|HFC236fa",
    "AR6 climate diagnostics|Emissions|HFC|HFC245fa",
    "AR6 climate diagnostics|Emissions|HFC|HFC32",
    "AR6 climate diagnostics|Emissions|HFC|HFC365mfc",
    "AR6 climate diagnostics|Emissions|HFC|HFC43-10",
    "AR6 climate diagnostics|Emissions|Halon1202",
    "AR6 climate diagnostics|Emissions|Halon1211", # 32
    "AR6 climate diagnostics|Emissions|Halon1301",
    "AR6 climate diagnostics|Emissions|Halon2402",
    "AR6 climate diagnostics|Emissions|N2O",  # 34
    "AR6 climate diagnostics|Emissions|NF3",
    "AR6 climate diagnostics|Emissions|NH3",  # 36
    "AR6 climate diagnostics|Emissions|NOx",  # 37
    "AR6 climate diagnostics|Emissions|NOx|Aviation",
    "AR6 climate diagnostics|Emissions|OC",  # 39
    "AR6 climate diagnostics|Emissions|PFC|C2F6",
    "AR6 climate diagnostics|Emissions|PFC|C3F8",
    "AR6 climate diagnostics|Emissions|PFC|C4F10",
    "AR6 climate diagnostics|Emissions|PFC|C5F12",
    "AR6 climate diagnostics|Emissions|PFC|C6F14",
    "AR6 climate diagnostics|Emissions|PFC|C7F16",
    "AR6 climate diagnostics|Emissions|PFC|C8F18",
    "AR6 climate diagnostics|Emissions|PFC|CF4", # 47
    "AR6 climate diagnostics|Emissions|PFC|cC4F8",
    "AR6 climate diagnostics|Emissions|SF6",
    "AR6 climate diagnostics|Emissions|SO2F2",
    "AR6 climate diagnostics|Emissions|Sulfur",  # 51
    "AR6 climate diagnostics|Emissions|VOC",  # 52
]

units = [
    "Mt BC/yr",
    "kt CCl4/yr",
    "kt CFC11/yr",
    "kt CFC113/yr",
    "kt CFC114/yr",
    "kt CFC115/yr",
    "kt CFC12/yr",
    "kt CH2Cl2/yr",
    "kt CH3Br/yr",
    "kt CH3CCl3/yr",
    "kt CH3Cl/yr",
    "Mt CH4/yr",
    "kt CHCl3/yr",
    "Mt CO/yr",
    "Mt CO2/yr",
    "Mt CO2/yr",
    "kt HCFC141b/yr",
    "kt HCFC142b/yr",
    "kt HCFC22/yr",
    "kt HFC125/yr",
    "kt HFC134a/yr",
    "kt HFC143a/yr",
    "kt HFC152a/yr",
    "kt HFC227ea/yr",
    "kt HFC23/yr",
    "kt HFC236fa/yr",
    "kt HFC245fa/yr",
    "kt HFC32/yr",
    "kt HFC365mfc/yr",
    "kt HFC4310/yr",
    "kt Halon1202/yr",
    "kt Halon1211/yr",
    "kt Halon1301/yr",
    "kt Halon2402/yr",
    "Mt N2O/yr",
    "kt NF3/yr",
    "Mt NH3/yr",
    "Mt NO2/yr",
    "Mt NO2/yr",
    "Mt OC/yr",
    "kt C2F6/yr",
    "kt C3F8/yr",
    "kt C4F10/yr",
    "kt C5F12/yr",
    "kt C6F14/yr",
    "kt C7F16/yr",
    "kt C8F18/yr",
    "kt CF4/yr",
    "kt cC4F8/yr",
    "kt SF6/yr",
    "kt SO2F2/yr",
    "Mt SO2/yr",
    "Mt VOC/yr",
]

for scenario in tqdm(range(1, RFF_SCENS + 1), desc="Making scenarios"):
    df_this = df.filter(scenario=scenario)
    emissions_2100_2300_out = np.zeros((201, 53))
    emissions_2100_2300_out[0, :] = df_this.filter(year=2100).data.value.values

    # we have to treat CO2, CH4 and N2O separately, as RFF give us these.

    # AFOLU emissions (positive for all except possibly CO2 AFOLU) remains constant at
    # 2100 levels. These fractions for SLCFs are pre-calculated as the
    # means from SSP scenarios (the big eight) in 2100.
    afolu_fraction = np.zeros((201, 53))
    afolu_fraction[:, 0] = 0.44547
    afolu_fraction[:, 13] = 0.55011
    afolu_fraction[:, 36] = 0.81015
    afolu_fraction[:, 37] = 0.33848
    afolu_fraction[:, 39] = 0.67538
    afolu_fraction[:, 51] = 0.09730
    afolu_fraction[:, 52] = 0.40963

    # CO2 AFOLU ramps to zero from 2100 to 2150. This is true for both positive and
    # negative emissions.
    afolu_fraction[:51, 14] = np.linspace(1, 0, 51)

    # Improvement over RCMIP: some natural emissions are not zero in 1750 for some
    # gases. They should not be zero in the future
    natural_emissions_1750 = np.zeros(53)
    natural_emissions_1750[47] = 0.010071225  # CF4
    natural_emissions_1750[1] = 0.024856862 # CCl4
    natural_emissions_1750[7] = 246.6579 # CH2Cl2
    natural_emissions_1750[8] = 105.08773 # CH3Br
    natural_emissions_1750[10] = 4275.7449 # CH3Cl
    natural_emissions_1750[12] = 300.92479 # CHCl3
    natural_emissions_1750[32] = 0.007723273 # Halon1211

    # net positive fossil emissions goes to zero in 2250. Fossil fraction is 1 for
    # minor GHGs, zero for CO2 AFOLU, and 1-(AFOLU fraction) for and SLCFs. Minor GHGs
    # with natural emissions go back to preindustrial in 2250.
    fossil_fraction = np.zeros((201, 53))
    fossil_fraction[0, :] = 1 - afolu_fraction[0, :]
    for ispec in range(53):
        emissions_1750_2100 = natural_emissions_1750[ispec]/emissions_2100_2300_out[0, ispec]
        fossil_fraction[:151, ispec] = (
            np.linspace(1, emissions_1750_2100, 151) * fossil_fraction[0, ispec]
        )
    fossil_fraction[:, 14] = 0

    # add up fractions and multiply by 2100 emissions
    emissions_2100_2300_out = emissions_2100_2300_out[0, :] * (
        afolu_fraction + fossil_fraction
    )

    # replace CO2 FFI, CH4 and N2O with RFF, where CO2 FFI is total minus assumed AFOLU
    df_in = pd.read_csv(DATARFF.joinpath("emissions%05d.csv" % scenario), index_col=0)
    co2 = df_in.loc[2100:2300, "CO2"]
    ch4 = df_in.loc[2100:2300, "CH4"]
    n2o = df_in.loc[2100:2300, "N2O"]
    emissions_2100_2300_out[:, 11] = ch4.values
    emissions_2100_2300_out[:, 34] = n2o.values
    emissions_2100_2300_out[:, 15] = co2.values * 1000 - emissions_2100_2300_out[:, 14]

    # concatenate and make dataframe
    emissions_2015_2300_out = np.concatenate(
        (df_this.timeseries().values, emissions_2100_2300_out[1:, :].T), axis=1
    )

    index_arrays = [
        ["RFF-SP"] * 53,
        [str(scenario)] * 53,
        ["World"] * 53,
        variables,
        units,
    ]
    index_tuples = list(zip(*index_arrays))
    index = pd.MultiIndex.from_tuples(
        index_tuples, names=["model", "scenario", "region", "variable", "unit"]
    )
    emissions_2015_2300_df = pd.DataFrame(
        emissions_2015_2300_out, columns=np.arange(2015, 2301), index=index
    )

    # Write CSV
    emissions_2015_2300_df.to_csv(
        DATAOUT.joinpath("emissions{:05d}.csv".format(scenario))
    )
