# Download RFF-SP emissions for CO2, CH4 and N2O

# Note: this downloads a 1.7 GB file into the cache, but we only use the RFF emissions which are around 240 MB in total. If you want to save disk space, at the expense of downloading the file anew every time you run the script, uncomment the cells below "Delete the zipfile from cache".

import os
import pathlib
import urllib.request
import zipfile

from dotenv import load_dotenv
import matplotlib.pyplot as pl
import numpy as np
import pandas as pd
from pathlib import Path
from tqdm.auto import tqdm

class DownloadProgressBar(tqdm):
    def update_to(self, b=1, bsize=1, tsize=None):
        if tsize is not None:
            self.total = tsize
        self.update(b * bsize - self.n)


def download_url(url, output_path):
    with DownloadProgressBar(unit='B', unit_scale=True,
                             miniters=1, desc=url.split('/')[-1]) as t:
        urllib.request.urlretrieve(url, filename=output_path, reporthook=t.update_to)

# Get environment variables
load_dotenv()

# Make data directory
DATADIR = Path(os.getenv("DATADIR"))
DATAOUT = DATADIR.joinpath('data_input')
os.makedirs(DATAOUT, exist_ok=True)

print("Downloading zipfile...")
RFF_ZIPFILE_TARGET = DATAOUT.joinpath("RFFSPs-Final.zip")
download_url(
   url = "https://zenodo.org/record/5898729/files/RFFSPs-Final.zip",
   output_path = RFF_ZIPFILE_TARGET
)

# Extract zipfile: only keep emissions files
print("Extracting zipfile...")
with zipfile.ZipFile(RFF_ZIPFILE_TARGET, mode='r') as z:
    [z.extract(file, path=DATAOUT) for file in z.namelist() if 'emissions/' in file]

# Process the RFF files into a format easier for FaIR to deal with

# We also want to attach the RCMIP/CMIP6 historical emissions on to this, so while we're at it, we'll download the SSP emissions from RCMIP.

# I believe that RFF used SSP2-4.5 between 2015 and 2020. **TODO** ask Marcus or re-read the Rennert paper
df_co2 = pd.read_csv(DATAOUT.joinpath('emissions/rffsp_co2_emissions.csv'))
df_ch4 = pd.read_csv(DATAOUT.joinpath('emissions/rffsp_ch4_emissions.csv'))
df_n2o = pd.read_csv(DATAOUT.joinpath('emissions/rffsp_n2o_emissions.csv'))

print("Downloading SSP emissions...")
SSP_EMISSIONS_TARGET = DATAOUT.joinpath("rcmip-emissions-annual-means-v5-1-0.csv")
ssp_emissions = download_url(
    url = "https://zenodo.org/record/4589756/files/rcmip-emissions-annual-means-v5-1-0.csv",
    output_path = SSP_EMISSIONS_TARGET
)

print("Making histories...")
df_ssp = pd.read_csv(SSP_EMISSIONS_TARGET)
co2_hist = df_ssp.loc[(df_ssp['Region']=='World')&(df_ssp['Scenario']=='ssp245')&(df_ssp['Variable']=='Emissions|CO2'),'1750':'2020'].interpolate(axis=1).values.squeeze()
ch4_hist = df_ssp.loc[(df_ssp['Region']=='World')&(df_ssp['Scenario']=='ssp245')&(df_ssp['Variable']=='Emissions|CH4'),'1750':'2020'].interpolate(axis=1).values.squeeze()
n2o_hist = df_ssp.loc[(df_ssp['Region']=='World')&(df_ssp['Scenario']=='ssp245')&(df_ssp['Variable']=='Emissions|N2O'),'1750':'2020'].interpolate(axis=1).values.squeeze()

DATAPROCESSED = DATADIR.joinpath('data_processed/emissions_files')
os.makedirs(DATAPROCESSED, exist_ok=True)

# would be slightly better to load in default molwts from fair, and much better to use fair's inbuilt
# unit converter
molwt_co2 = 44.009
molwt_c   = 12.011
molwt_n2o = 44.013
molwt_n2  = 28.014
mt_to_gt  = 0.001

for sample in tqdm(range(1, 10001), desc="Processing emissions files"):
    emissions = np.zeros((551, 3))
    co2 = df_co2[df_co2['sample']==sample].value.values
    ch4 = df_ch4[df_ch4['sample']==sample].value.values
    n2o = df_n2o[df_n2o['sample']==sample].value.values
    emissions[:270, 0] = co2_hist[:-1] * mt_to_gt
    emissions[:270, 1] = ch4_hist[:-1]
    emissions[:270, 2] = n2o_hist[:-1] * mt_to_gt
    emissions[270:, 0] = co2 * molwt_co2 / molwt_c
    emissions[270:, 1] = ch4
    emissions[270:, 2] = n2o * molwt_n2o / molwt_n2
    df_out = pd.DataFrame(emissions, columns=['CO2', 'CH4', 'N2O'], index=range(1750,2301))
    df_out.to_csv(DATAPROCESSED.joinpath('emissions%05d.csv' % sample))

print("Cleaning up...")
# Remove intermediate RFF datafiles and reclaim 240 MB
if os.path.exists(DATAOUT.joinpath('emissions/rffsp_co2_emissions.csv')):
    os.remove(DATAOUT.joinpath('emissions/rffsp_co2_emissions.csv'))
if os.path.exists(DATAOUT.joinpath('emissions/rffsp_ch4_emissions.csv')):
    os.remove(DATAOUT.joinpath('emissions/rffsp_ch4_emissions.csv'))
if os.path.exists(DATAOUT.joinpath('emissions/rffsp_n2o_emissions.csv')):
    os.remove(DATAOUT.joinpath('emissions/rffsp_n2o_emissions.csv'))
if os.path.exists(DATAOUT.joinpath('emissions')):
    os.removedirs(DATAOUT.joinpath('emissions'))

# Delete the original zipfile and reclaim 1.7 GB
if os.path.exists(RFF_ZIPFILE_TARGET):
    os.remove(RFF_ZIPFILE_TARGET)
