# rff-fair2.1

Run the RFF-SPs using FaIR v2.1 (calibrated and constrained AR6 version) including internal variability. This repository should reproduce the results.

## Prerequisites

- python 3.7+
- anaconda strongly recommended
- 220 GB disk space: the output data consumes around 205 GB
- the AR6 emissions infilling database. Save this as `$DATADIR/data_input/ar6_emissions_vetted_infillerdatabase.csv` where `$DATADIR` is an environment variable you define (see 5. below). You have to obtain this by setting up an account with the IIASA AR6 WG3 database. Do this, then visit the Downloads page (https://data.ene.iiasa.ac.at/ar6/#/downloads), and download the file called "Infiller database for silicone: IPCC AR6 WGIII version (DOI: 10.5281/zenodo.6390768)"

## Reproduction

1. Clone this repository to your local machine: `git clone git@github.com:chrisroadmap/rff-fair2.1.git`
2. If using `conda`, create the environment:
```
cd rff-fair2.1
conda env create -f environment.yml
```
3. Activate the environment:
```
conda activate rff-fair2.1
```
4. If you want to make nice version-control friendly notebooks, which will remove all output and data upon committing, run
```
nbstripout --install
```
5. Create a `.env` file in the root directory of the repository, and populate it with these two environment variables:
```
DATADIR=/path/to/datafiles/   # change this to a local path where you want your output stored
RFF_SCENS=10000               # how many RFF scenarios (full population = 10000) to run? For 
                              # testing, use a small number as the code takes a lot of time.
```
6. Run the scripts in the `scripts` directory in order.
