import multiprocessing
import pathlib
import warnings
from concurrent.futures import ProcessPoolExecutor

import matplotlib.pyplot as pl
import numpy as np
import pandas as pd
import pyam
import silicone.database_crunchers
from silicone import (
    multiple_infillers,  # .decompose_collection_with_time_dep_ratio.DecomposeCollectionTimeDepRatio
)
from silicone.stats import rolling_window_find_quantiles
from silicone.utils import return_cases_which_consistently_split
from tqdm.auto import tqdm
from utils import _parallel_process

# Get environment variables
load_dotenv()

DATADIR = Path(os.getenv("DATADIR"))
DATAIN = DATADIR.joinpath("data_input")
DATARFF = DATADIR.joinpath("data_processed/emissions_files")

# number of processors
WORKERS = multiprocessing.cpu_count()


print("Making RFF's CH4 and N2O into a pyam for completeness...")
dfs = []
for sample in tqdm(range(1, RFF_SCENS + 1)):
    df_in = pd.read_csv(DATARFF.joinpath("emissions%05d.csv" % sample), index_col=0)
    ch4 = df_in["CH4"]
    n2o = df_in["N2O"]
    ch4_data = pd.DataFrame(ch4, index=np.arange(2020, 2101))
    n2o_data = pd.DataFrame(n2o, index=np.arange(2020, 2101))
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
