# The purpose of this script is to find the ratio of fossil to total emissions in 2100
# in the SSP scenarios. We will use these for assumed SSP projections.

# Meinshausen et al. 2020, section 2.3.
import os
from pathlib import Path

import pyam
from dotenv import load_dotenv

# Get environment variables
load_dotenv()

DATADIR = Path(os.getenv("DATADIR"))
DATAIN = DATADIR.joinpath("data_input")

df = pyam.IamDataFrame(DATAIN.joinpath("rcmip-emissions-annual-means-v5-1-0.csv"))

variables = [
    "Emissions|BC|MAGICC AFOLU",
    "Emissions|BC",
    "Emissions|CH4|MAGICC AFOLU",
    "Emissions|CH4",
    "Emissions|CO2|MAGICC AFOLU",
    "Emissions|CO2",
    "Emissions|CO|MAGICC AFOLU",
    "Emissions|CO",
    "Emissions|N2O|MAGICC AFOLU",
    "Emissions|N2O",
    "Emissions|NH3|MAGICC AFOLU",
    "Emissions|NH3",
    "Emissions|NOx|MAGICC AFOLU",
    "Emissions|NOx",
    "Emissions|OC|MAGICC AFOLU",
    "Emissions|OC",
    "Emissions|Sulfur|MAGICC AFOLU",
    "Emissions|Sulfur",
    "Emissions|VOC|MAGICC AFOLU",
    "Emissions|VOC",
]

scenarios = ["ssp119", "ssp126", "ssp245", "ssp370", "ssp434", "ssp460", "ssp534-over", "ssp585"]

df = df.filter(region="World", variable=variables, scenario=scenarios, year=2100)

afolu = variables[0::2]
total = variables[1::2]

df_out = []
for ivar in range(len(variables)//2):
    frac_name = f"Fraction|{total[ivar].split('|')[-1]}"
    if frac_name=="Fraction|CO2":
        continue
    df_out.append(df.divide(afolu[ivar], total[ivar], frac_name, ignore_units='1'))

df_out = pyam.concat(df_out)

# should save this and load results in next script - I've done by hand in Excel to
# save time.
