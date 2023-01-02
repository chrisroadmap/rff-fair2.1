from pathlib import Path
import os
import warnings

from dotenv import load_dotenv
import matplotlib.pyplot as pl
import numpy as np
import pandas as pd
import pyam
import silicone.database_crunchers
from silicone.stats import rolling_window_find_quantiles
from silicone import multiple_infillers
from silicone.utils import return_cases_which_consistently_split
from tqdm.auto import tqdm

# Get environment variables
load_dotenv()

DATADIR = Path(os.getenv("DATADIR"))
DATAIN = DATADIR.joinpath('data_input')
DATARFF = DATADIR.joinpath('data_processed/emissions_files')
DATAOUT = DATADIR.joinpath('data_processed')

# number of scenarios
RFF_SCENS = int(os.getenv("RFF_SCENS"))

print("Making RFF's CH4 and N2O into a pyam for completeness...")
dfs = []
for sample in tqdm(range(1, RFF_SCENS+1)):
    df_in = pd.read_csv(DATARFF.joinpath('emissions%05d.csv' % sample), index_col=0)
    ch4 = df_in['CH4']
    n2o = df_in['N2O']
    ch4_data = pd.DataFrame(ch4, index=np.arange(2020, 2101))
    n2o_data = pd.DataFrame(n2o, index=np.arange(2020, 2101))
    dfs.append(
        pyam.IamDataFrame(
            ch4_data.T,
            model="RFF-SP",
            scenario="{:05d}".format(sample),
            region="World",
            unit="Mt CH4/yr",
            variable='AR6 climate diagnostics|Emissions|CH4',
        )
    )
    dfs.append(
        pyam.IamDataFrame(
            n2o_data.T,
            model="RFF-SP",
            scenario="{:05d}".format(sample),
            region="World",
            unit="Mt N2O/yr",
            variable='AR6 climate diagnostics|Emissions|N2O',
        )
    )

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    pyam_ch4n2o_data = pyam.concat(dfs)

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
        "AR6 climate diagnostics|Emissions|N2O"
    ],
keep=False)

print("Making RFF's CO2 into a pyam...")
dfs = []
for sample in tqdm(range(1, RFF_SCENS+1)):
    df_in = pd.read_csv(DATARFF.joinpath('emissions%05d.csv' % sample), index_col=0)
    co2 = df_in['CO2']
    co2_data = pd.DataFrame(co2, index=np.arange(2020, 2101))
    dfs.append(
        pyam.IamDataFrame(
            co2_data.T*1000,
            model="RFF-SP",
            scenario="{:05d}".format(sample),
            region="World",
            unit="Mt CO2/yr",
            variable='AR6 climate diagnostics|Emissions|CO2',
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


print("Doing the infilling...")
database_species_except_total_co2 = [
    specie for specie in infiller_database.variable if specie not in [
        "AR6 climate diagnostics|Emissions|CO2",
        "AR6 climate diagnostics|Emissions|CO2|Energy and Industrial Processes",
        "AR6 climate diagnostics|Emissions|CO2|AFOLU",
    ]
]

# TODO: save out updated infiller database, SSP infiller database and CO2 E&IP
# these will be loaded each time into parallel.py
# the IO overhead is substantial, but infilling is slower still, and this saves RAM
# which is always tight on the HPC.

conf = []
for sample in range(1, RFF_SCENS+1):
    conf.append(
        {
            'sample': sample,
            'database_species_except_total_co2': database_species_except_total_co2,
            'infiller_database': infiller_database,
            'pyam_co2_data': pyam_co2_data
        }
    )
