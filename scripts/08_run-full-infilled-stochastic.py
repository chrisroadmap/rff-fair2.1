# # Run FaIR v2.1 with RFF scenarios
#
# - use all 10000 emissions scenarios
# - run all 1001 ensemble members
# - use infilled emissions scenarios
# - stochastic variability is ON
# - take advantage of parallel processing

import multiprocessing
import os
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

from dotenv import load_dotenv
from parallel_fair import run_fair
from utils import _parallel_process

if __name__ == "__main__":
    print("Starting RFF runs...")

    # Get environment variables
    load_dotenv()

    # Make data directory
    DATADIR = Path(os.getenv("DATADIR"))
    DATAIN = DATADIR.joinpath("data_input")
    DATAPROCESSED = DATADIR.joinpath("data_processed", "infilled_extended")
    DATAOUT = DATADIR.joinpath("data_output", "stochastic")
    os.makedirs(DATAOUT, exist_ok=True)

    # number of processors
    WORKERS = multiprocessing.cpu_count()

    # number of scenarios
    RFF_SCENS = int(os.getenv("RFF_SCENS"))

    conf = list(range(1, RFF_SCENS + 1))

    parallel_process_kwargs = dict(
        func=run_fair,
        configuration=conf,
        config_are_kwargs=False,
    )

    print("Running infilling in parallel...")
    with ProcessPoolExecutor(WORKERS) as pool:
        _parallel_process(
            **parallel_process_kwargs,
            pool=pool,
        )
