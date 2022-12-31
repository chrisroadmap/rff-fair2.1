from concurrent.futures import ProcessPoolExecutor
import multiprocessing
import warnings

import matplotlib.pyplot as pl
import numpy as np
import pandas as pd
import pyam
import silicone.database_crunchers
from silicone.stats import rolling_window_find_quantiles
from silicone import multiple_infillers#.decompose_collection_with_time_dep_ratio.DecomposeCollectionTimeDepRatio
from silicone.utils import return_cases_which_consistently_split
from tqdm.auto import tqdm

from utils import _parallel_process

# number of processors
WORKERS = multiprocessing.cpu_count()

# number of scenarios
#RFF_SCENS = 10000
RFF_SCENS = 4

# and here is the parallel function
def run_stuff(stuff):
    sample = stuff["sample"]
    infiller_database = stuff["infiller_database"]
    database_species_except_total_co2 = stuff["database_species_except_total_co2"]
    pyam_co2_data = stuff["pyam_co2_data"]

    inner_list = []
    lead = ["AR6 climate diagnostics|Infilled|Emissions|CO2"]
    cruncher = silicone.database_crunchers.QuantileRollingWindows(infiller_database)
    for follow in database_species_except_total_co2:
        filler = cruncher.derive_relationship(follow, lead)
        filler_input = pyam_co2_data.filter(model="RFF-SP", scenario="{:05d}".format(sample))
        filler_data = filler_input.filter(year=infiller_database["year"].unique(), variable=lead)
        qrw_infilled = filler(filler_data)
        inner_list.append(qrw_infilled.filter(variable=follow))
    pyam_emissions = pyam.concat(inner_list)
    return(pyam_emissions)

# everything under this block we only want to run once, on setup
if __name__ == '__main__':
    print("Making RFF's CH4 and N2O into a pyam for completeness...")
    # zzzz
    dfs = []
    for sample in tqdm(range(1, RFF_SCENS+1)):
        df_in = pd.read_csv('../data_processed/emissions_files/emissions%05d.csv' % sample, index_col=0)
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
                variable='AR6 climate diagnostics|Infilled|Emissions|CH4',
            )
        )
        dfs.append(
            pyam.IamDataFrame(
                n2o_data.T,
                model="RFF-SP",
                scenario="{:05d}".format(sample),
                region="World",
                unit="Mt N2O/yr",
                variable='AR6 climate diagnostics|Infilled|Emissions|N2O',
            )
        )

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        pyam_ch4n2o_data = pyam.concat(dfs)

    # # this is so slow, I saved the data in a binary format
    # df = pd.read_excel('../data_input/20220314_ar6emissions_harmonized_infilled.xlsx')
    # df.to_pickle('../data_input/20220314_ar6emissions_harmonized_infilled.pkl')

    # this file isn't public, but we want it to be : place on Zenodo, point all to Kikstra et al citation
    # also, why is reading this in so painful?
    print("Reading in AR6 emissions...")
    df = pd.read_pickle('../data_input/20220314_ar6emissions_harmonized_infilled.pkl')
    infiller_database = pyam.IamDataFrame(df)

    # Filter database
    print("Filtering variables...")
    infiller_database = infiller_database.filter(variable="AR6 climate diagnostics|Infilled|Emissions|*")

    # Filter more
    print("Keeping only model/scenario pairs that passed vetting...")
    df_scens = pd.read_csv('../data_input/ar6_model_scenario_passed_vetting.csv')
    model_scen_pairs = []
    for irow, row in df_scens.iterrows():
        model_scen_pairs.append((row['model'], row['scenario']))
    model_scen_pairs = sorted(model_scen_pairs)

    the_slowness = []
    for model, scen in tqdm(model_scen_pairs, desc="Making a pyam of this"):
        the_slowness.append(infiller_database.filter(model=model, scenario=scen))
    infiller_database = pyam.concat(the_slowness)

    print("Making an aggregate CO2 category...")
    infiller_database.add(
        "AR6 climate diagnostics|Infilled|Emissions|CO2|AFOLU",
        "AR6 climate diagnostics|Infilled|Emissions|CO2|Energy and Industrial Processes",
        "AR6 climate diagnostics|Infilled|Emissions|CO2",
        axis='variable',
        fillna=None,
        ignore_units='Mt CO2/yr',
        append=True
    )

    print("Removing variables we don't want to infill...")
    database_species = infiller_database.variable

    # Remove CH4 and N2O which are not being infilled
    database_species.remove('AR6 climate diagnostics|Infilled|Emissions|CH4')
    database_species.remove('AR6 climate diagnostics|Infilled|Emissions|N2O')

#    # Remove species which do not vary
#    database_species.remove('AR6 climate diagnostics|Infilled|Emissions|CCl4')
#    database_species.remove('AR6 climate diagnostics|Infilled|Emissions|CFC11')
#    database_species.remove('AR6 climate diagnostics|Infilled|Emissions|CFC113')
#    database_species.remove('AR6 climate diagnostics|Infilled|Emissions|CFC114')
#    database_species.remove('AR6 climate diagnostics|Infilled|Emissions|CFC115')
#    database_species.remove('AR6 climate diagnostics|Infilled|Emissions|CFC12')

#    # Remove species that are all zero
#    database_species.remove('AR6 climate diagnostics|Infilled|Emissions|HFC|HFC245ca')
#    database_species.remove('AR6 climate diagnostics|Infilled|Emissions|CH3CCl3')
#    database_species.remove('AR6 climate diagnostics|Infilled|Emissions|Halon1202')

    # Remove aggregates
    database_species.remove('AR6 climate diagnostics|Infilled|Emissions|HFC')
    database_species.remove('AR6 climate diagnostics|Infilled|Emissions|PFC')
    database_species.remove('AR6 climate diagnostics|Infilled|Emissions|F-Gases')

    print("Making RFF's CO2 into a pyam...")
    # zzzz
    dfs = []
    for sample in tqdm(range(1, RFF_SCENS+1)):
        df_in = pd.read_csv('../data_processed/emissions_files/emissions%05d.csv' % sample, index_col=0)
        co2 = df_in['CO2']
        co2_data = pd.DataFrame(co2, index=np.arange(2020, 2101))
        dfs.append(
            pyam.IamDataFrame(
                co2_data.T*1000,
                model="RFF-SP",
                scenario="{:05d}".format(sample),
                region="World",
                unit="Mt CO2/yr",
                variable='AR6 climate diagnostics|Infilled|Emissions|CO2',
            )
        )

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        pyam_co2_data = pyam.concat(dfs)


    print("Decomposing RFF total CO2 into AFOLU and FFI based on database...")
    components = [
        "AR6 climate diagnostics|Infilled|Emissions|CO2|Energy and Industrial Processes",
        "AR6 climate diagnostics|Infilled|Emissions|CO2|AFOLU",
    ]
    aggregate = "AR6 climate diagnostics|Infilled|Emissions|CO2"
    to_infill = pyam_co2_data.filter(variable=aggregate)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        decomposer = multiple_infillers.DecomposeCollectionTimeDepRatio(infiller_database)
        co2_results = decomposer.infill_components(aggregate, components, to_infill)


    print("Doing the infilling...")
    database_species_except_total_co2 = [
        specie for specie in database_species if specie not in [
            "AR6 climate diagnostics|Infilled|Emissions|CO2",
            "AR6 climate diagnostics|Infilled|Emissions|CO2|Energy and Industrial Processes",
            "AR6 climate diagnostics|Infilled|Emissions|CO2|AFOLU",
        ]
    ]

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

    parallel_process_kwargs = dict(
        func=run_stuff,
        configuration=conf,
        config_are_kwargs=False,
    )

    with ProcessPoolExecutor(WORKERS) as pool:
        res = _parallel_process(
            **parallel_process_kwargs,
            pool=pool,
        )

    res = pyam.concat(res)
    res = res.append(co2_results)
    res = res.append(pyam_ch4n2o_data)

    print("Saving out...")
    res.to_csv('../data_processed/infilled_emissions_scenarios.csv')
