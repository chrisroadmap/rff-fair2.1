import multiprocessing
import os
import pickle
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import pyam
from dotenv import load_dotenv
from parallel import run_stuff
from utils import _parallel_process

# Get environment variables
load_dotenv()

DATADIR = Path(os.getenv("DATADIR"))
DATAIN = DATADIR.joinpath("data_input")
DATAINFILLED = DATADIR.joinpath("data_processed/infilling")
DATAOUT = DATADIR.joinpath("data_processed")

# number of processors
WORKERS = multiprocessing.cpu_count()

# number of scenarios
RFF_SCENS = int(os.getenv("RFF_SCENS"))

# everything under this block we only want to run once, this is why we have the
# __main__ conditional
if __name__ == "__main__":

    conf = list(range(1, RFF_SCENS + 1))

    parallel_process_kwargs = dict(
        func=run_stuff,
        configuration=conf,
        config_are_kwargs=False,
    )

    print("Running infilling in parallel...")
    with ProcessPoolExecutor(WORKERS) as pool:
        res = _parallel_process(
            **parallel_process_kwargs,
            pool=pool,
        )

    print("Concatenating results...")
    res = res[0]  # comes out as nested list
    res = pyam.concat(res)

    print("Adding on CO2, CH4, N2O...")
    with open(DATAINFILLED.joinpath("rff_co2_eip.pkl"), "rb") as handle:
        co2_eip = pickle.load(handle)

    with open(DATAINFILLED.joinpath("rff_co2_afolu.pkl"), "rb") as handle:
        co2_afolu = pickle.load(handle)

    with open(DATAINFILLED.joinpath("rff_ch4_n2o.pkl"), "rb") as handle:
        ch4_n2o = pickle.load(handle)

    res = res.append(co2_eip)
    res = res.append(co2_afolu)
    res = res.append(ch4_n2o)

    print("Saving out...")
    res.to_csv(DATAOUT.joinpath("infilled_emissions_scenarios.csv"))
