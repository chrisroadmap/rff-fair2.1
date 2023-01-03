# Download CMIP6 SSP infiller database

import os
import shutil
import urllib.request
import zipfile
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

print("Downloading zipfile...")
RFF_ZIPFILE_TARGET = DATAOUT.joinpath("climate-assessment-v0.1.1.zip")
download_url(
    url="https://zenodo.org/record/6782457/files/iiasa/climate-assessment-v0.1.1.zip",
    output_path=RFF_ZIPFILE_TARGET,
)

# Extract zipfile: only keep emissions files
print("Extracting zipfile...")
with zipfile.ZipFile(RFF_ZIPFILE_TARGET, mode="r") as z:
    [
        z.extract(file, path=DATAOUT)
        for file in z.namelist()
        if "cmip6-ssps-workflow-emissions.csv" in file
    ]


print("Moving to a handy location...")
shutil.move(
    DATAOUT.joinpath(
        "iiasa-climate-assessment-f108552/src/climate_assessment/infilling/cmip6-ssps-workflow-emissions.csv"
    ),
    DATAOUT,
)
shutil.rmtree(DATAOUT.joinpath("iiasa-climate-assessment-f108552"))


# Delete the original zipfile
if os.path.exists(RFF_ZIPFILE_TARGET):
    os.remove(RFF_ZIPFILE_TARGET)
