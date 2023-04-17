# Run constrained SSP projections (emissions driven) for completeness.

import os
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
from dotenv import load_dotenv
from fair import FAIR
from fair.interface import fill, initialise
from fair.io import read_properties
from fair.structure.units import *
from tqdm.auto import tqdm

load_dotenv()

DATADIR = Path(os.getenv("DATADIR"))
DATAIN = DATADIR.joinpath("data_input")
DATAOUT = DATADIR.joinpath("data_output", "stochastic")

print("Running RCP scenarios...")


class DownloadProgressBar(tqdm):
    def update_to(self, b=1, bsize=1, tsize=None):
        if tsize is not None:
            self.total = tsize
        self.update(b * bsize - self.n)


def download_url(url, output_path):
    with DownloadProgressBar(
        unit="B", unit_scale=True, miniters=1, desc=url.split("/")[-1]
    ) as t:
        urllib.request.urlretrieve(url, filename=output_path, reporthook=t.update_to)


print("Downloading aviation NOx fractions...")
TARGET = DATAIN.joinpath("aviNOx_fraction.csv")
download_url(
    url="https://github.com/OMS-NetZero/FAIR/raw/v1.3.6/fair/RCPs/data/aviNOx_fraction.csv",
    output_path=TARGET,
)
