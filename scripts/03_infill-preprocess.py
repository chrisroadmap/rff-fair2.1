import os
import pickle
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import pyam
from dotenv import load_dotenv
from silicone import multiple_infillers
from tqdm.auto import tqdm

# kill pyam/silicone/pandas
warnings.simplefilter("ignore")

# Get environment variables
load_dotenv()

DATADIR = Path(os.getenv("DATADIR"))
DATAIN = DATADIR.joinpath("data_input")
DATARFF = DATADIR.joinpath("data_processed/emissions_files")
DATAOUT = DATADIR.joinpath("data_processed/infilling")

os.makedirs(DATAOUT, exist_ok=True)

# number of scenarios
RFF_SCENS = int(os.getenv("RFF_SCENS"))

# this file isn't public, but we want it to be
print("Reading in AR6 infiller database...")
df = pd.read_csv(DATAIN.joinpath("ar6_emissions_vetted_infillerdatabase.csv"))
infiller_database = pyam.IamDataFrame(df)

# Renaming all variables to drop "Harmonized"
infill_varlist = infiller_database.variable
mapping = {x: x.replace("Harmonized|", "") for x in infill_varlist}
infiller_database = infiller_database.rename(mapping={"variable": mapping})

# Remove CH4 and N2O which are not being infilled
print("Removing variables we don't want to infill...")
infiller_database = infiller_database.filter(
    variable=[
        "AR6 climate diagnostics|Emissions|CH4",
        "AR6 climate diagnostics|Emissions|N2O",
        "AR6 climate diagnostics|Emissions|F-Gases",
    ],
    keep=False,
)

print("Making RFF's CO2 into a pyam...")
dfs = []
for sample in tqdm(range(1, RFF_SCENS + 1)):
    df_in = pd.read_csv(DATARFF.joinpath("emissions%05d.csv" % sample), index_col=0)
    co2 = df_in["CO2"]
    co2_data = pd.DataFrame(co2, index=np.arange(2015, 2101))
    dfs.append(
        pyam.IamDataFrame(
            co2_data.T * 1000,
            model="RFF-SP",
            scenario="{:05d}".format(sample),
            region="World",
            unit="Mt CO2/yr",
            variable="AR6 climate diagnostics|Emissions|CO2",
        )
    )

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    pyam_co2_data = pyam.concat(dfs)


print("Decomposing RFF total CO2 into AFOLU and FFI based on database...")
components = [
    "AR6 climate diagnostics|Emissions|CO2|Energy and Industrial Processes",
    "AR6 climate diagnostics|Emissions|CO2|AFOLU",
]
aggregate = "AR6 climate diagnostics|Emissions|CO2"
to_infill = pyam_co2_data.filter(variable=aggregate)
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    decomposer = multiple_infillers.DecomposeCollectionTimeDepRatio(infiller_database)
    co2_results = decomposer.infill_components(aggregate, components, to_infill)


print("Slimming down CMIP6 SSP file...")
df = pd.read_csv(DATAIN.joinpath("cmip6-ssps-workflow-emissions.csv"))
infiller_database_ssp = pyam.IamDataFrame(df)

print("Grafting in aviation NOx")
scenarios = ["ssp119", "ssp126", "ssp245", "ssp370", "ssp434", "ssp460", "ssp534-over", "ssp585"]
df_rcmip = pd.read_csv(DATAIN.joinpath("rcmip-emissions-annual-means-v5-1-0.csv"))
models = {
    'ssp119': 'IMAGE',
    'ssp126': 'IMAGE',
    'ssp245': 'MESSAGE-GLOBIOM',
    'ssp370': 'AIM/CGE',
    'ssp434': 'GCAM4',
    'ssp460': 'GCAM4',
    'ssp534-over': 'REMIND-MAGPIE',
    'ssp585': 'REMIND-MAGPIE'
}
for scenario in scenarios:
    df_part = df_rcmip.loc[
        (df_rcmip["Scenario"]==scenario) &
        (df_rcmip["Region"]=="World") &
        (df_rcmip["Variable"]=="Emissions|NOx|MAGICC Fossil and Industrial|Aircraft"), "2015":"2100"
    ].interpolate(axis=1)
    pyam_part = pyam.IamDataFrame(
        df_part,
        model=models[scenario],
        scenario=scenario,
        region="World",
        unit="Mt NO2/yr",
        variable="Emissions|NOx|Aviation",
    )
    infiller_database_ssp.append(pyam_part, inplace=True)

# Keep only required variables + CO2 EIP
minor_ghg_variables_list = [
    "Emissions|CO2|Energy and Industrial Processes",
    "Emissions|CCl4",
    "Emissions|CFC11",
    "Emissions|CFC113",
    "Emissions|CFC114",
    "Emissions|CFC115",
    "Emissions|CFC12",
    "Emissions|CH2Cl2",
    "Emissions|CH3Br",
    "Emissions|CH3CCl3",
    "Emissions|CH3Cl",
    "Emissions|CHCl3",
    "Emissions|HCFC141b",
    "Emissions|HCFC142b",
    "Emissions|HCFC22",
    "Emissions|HFC|HFC152a",
    "Emissions|HFC|HFC236fa",
    "Emissions|HFC|HFC245fa",  # think incorrectly named
    "Emissions|HFC|HFC365mfc",
    "Emissions|Halon1202",
    "Emissions|Halon1211",
    "Emissions|Halon1301",
    "Emissions|Halon2402",
    "Emissions|NF3",
    #    "Emissions|PFC|CF4",
    #    "Emissions|PFC|C2F6",
    "Emissions|PFC|C3F8",
    "Emissions|PFC|C4F10",
    "Emissions|PFC|C5F12",
    #    "Emissions|PFC|C6F14",
    "Emissions|PFC|C7F16",
    "Emissions|PFC|C8F18",
    "Emissions|PFC|cC4F8",
    "Emissions|SO2F2",
    "Emissions|NOx|Aviation",
]
infiller_database_ssp = infiller_database_ssp.filter(
    variable=minor_ghg_variables_list, year=range(2015, 2101)
)
# Renaming all variables to insert prefix
infill_varlist_ssp = infiller_database_ssp.variable
mapping = {
    x: x.replace("Emissions|", "AR6 climate diagnostics|Emissions|")
    for x in infill_varlist_ssp
}
infiller_database_ssp = infiller_database_ssp.rename(mapping={"variable": mapping})


print("Making RFF's CH4 and N2O into a pyam for completeness...")
dfs = []
for sample in tqdm(range(1, RFF_SCENS + 1)):
    df_in = pd.read_csv(DATARFF.joinpath("emissions%05d.csv" % sample), index_col=0)
    ch4 = df_in["CH4"]
    n2o = df_in["N2O"]
    ch4_data = pd.DataFrame(ch4, index=np.arange(2015, 2101))
    n2o_data = pd.DataFrame(n2o, index=np.arange(2015, 2101))
    dfs.append(
        pyam.IamDataFrame(
            ch4_data.T,
            model="RFF-SP",
            scenario="{:05d}".format(sample),
            region="World",
            unit="Mt CH4/yr",
            variable="AR6 climate diagnostics|Emissions|CH4",
        )
    )
    dfs.append(
        pyam.IamDataFrame(
            n2o_data.T,
            model="RFF-SP",
            scenario="{:05d}".format(sample),
            region="World",
            unit="Mt N2O/yr",
            variable="AR6 climate diagnostics|Emissions|N2O",
        )
    )

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    pyam_ch4n2o_data = pyam.concat(dfs)


print("Saving out files as pickles because ;)...")

# Save out infiller database
with open(DATAOUT.joinpath("infiller_database_ar6.pkl"), "wb") as handle:
    pickle.dump(infiller_database, handle, protocol=pickle.HIGHEST_PROTOCOL)

# Save out other infiller database
with open(DATAOUT.joinpath("infiller_database_ssp.pkl"), "wb") as handle:
    pickle.dump(infiller_database_ssp, handle, protocol=pickle.HIGHEST_PROTOCOL)

# Save out RFF CO2 AFOLU
co2_afolu = co2_results.filter(variable="AR6 climate diagnostics|Emissions|CO2|AFOLU")
with open(DATAOUT.joinpath("rff_co2_afolu.pkl"), "wb") as handle:
    pickle.dump(co2_afolu, handle, protocol=pickle.HIGHEST_PROTOCOL)

# Save out RFF CO2 E&IP
co2_eip = co2_results.filter(
    variable="AR6 climate diagnostics|Emissions|CO2|Energy and Industrial Processes"
)
with open(DATAOUT.joinpath("rff_co2_eip.pkl"), "wb") as handle:
    pickle.dump(co2_eip, handle, protocol=pickle.HIGHEST_PROTOCOL)

# Save out CH4 and N2O
with open(DATAOUT.joinpath("rff_ch4_n2o.pkl"), "wb") as handle:
    pickle.dump(pyam_ch4n2o_data, handle, protocol=pickle.HIGHEST_PROTOCOL)
