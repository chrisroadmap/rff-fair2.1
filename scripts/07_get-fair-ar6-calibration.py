# Download RFF-SP emissions for CO2, CH4 and N2O

import os
import urllib.request
from pathlib import Path

from dotenv import load_dotenv
from tqdm.auto import tqdm


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


# Get environment variables
load_dotenv()

# Make data directory
DATADIR = Path(os.getenv("DATADIR"))
DATAOUT = DATADIR.joinpath("data_input")
os.makedirs(DATAOUT, exist_ok=True)

print("Downloading calibration data...")
TARGET = DATAOUT.joinpath("calibrated_constrained_parameters.csv")
download_url(
    url="https://zenodo.org/record/7545157/files/calibrated_constrained_parameters.csv",
    output_path=TARGET,
)
