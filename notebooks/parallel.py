import pyam
from silicone.database_crunchers import QuantileRollingWindows, RMSClosest
from silicone.multiple_infillers import infill_all_required_variables

# straight from Jarmo. CO2 AFOLU, CH4, N2O already done
# we do this with QRW
major_variables_list = [
#    "Emissions|BC",
    "Emissions|CH4",
#    "Emissions|CO2|AFOLU",
    "Emissions|CO",
#    "Emissions|N2O",
    "Emissions|NH3",
    "Emissions|NOx",
    "Emissions|OC",
    "Emissions|Sulfur",
    "Emissions|VOC",
]

# We do this with RMS Closest
main_hfc_pfc_variables_list = [
    "Emissions|HFC|HFC134a",
    "Emissions|HFC|HFC143a",
    "Emissions|HFC|HFC227ea",
    "Emissions|HFC|HFC23",
    "Emissions|HFC|HFC32",
    "Emissions|HFC|HFC43-10",
    "Emissions|HFC|HFC245ca",
    "Emissions|HFC|HFC125",
    "Emissions|SF6",
    "Emissions|PFC|CF4",
    "Emissions|PFC|C2F6",
    "Emissions|PFC|C6F14",
]

# We do this with RMS Closest, using SSP emissions database.
minor_ghg_variables_list = [
    "Emissions|CCl4",
    "Emissions|CFC11",
    "Emissions|CFC113",
    "Emissions|CFC114",
    "Emissions|CFC115",
    "Emissions|CFC12",
    "Emissions|CH2Cl2",
    "Emissions|CH3Br",
    "Emissions|CH3CCl3",
    "Emissions|CH3Cl",
    "Emissions|CHCl3",
    "Emissions|HCFC141b",
    "Emissions|HCFC142b",
    "Emissions|HCFC22",
    "Emissions|HFC|HFC152a",
    "Emissions|HFC|HFC236fa",
    # 'Emissions|HFC|HFC245fa',
    "Emissions|HFC|HFC365mfc",
    "Emissions|Halon1202",
    "Emissions|Halon1211",
    "Emissions|Halon1301",
    "Emissions|Halon2402",
    "Emissions|NF3",
    "Emissions|PFC|C3F8",
    "Emissions|PFC|C4F10",
    "Emissions|PFC|C5F12",
    "Emissions|PFC|C7F16",
    "Emissions|PFC|C8F18",
    "Emissions|PFC|cC4F8",
    "Emissions|SO2F2",
]

def run_stuff(stuff):
    sample = stuff["sample"]
    infiller_database = stuff["infiller_database"]
    database_species_except_total_co2 = stuff["database_species_except_total_co2"]
    pyam_co2_data = stuff["pyam_co2_data"]

    inner_list = []
    lead = ["AR6 climate diagnostics|Emissions|CO2"]
    cruncher = QuantileRollingWindows(infiller_database)
    for follow in database_species_except_total_co2:
        filler = cruncher.derive_relationship(follow, lead)
        filler_input = pyam_co2_data.filter(model="RFF-SP", scenario="{:05d}".format(sample))
        filler_data = filler_input.filter(year=infiller_database["year"].unique(), variable=lead)
        qrw_infilled = filler(filler_data)
        inner_list.append(qrw_infilled.filter(variable=follow))
    pyam_emissions = pyam.concat(inner_list)
    return(pyam_emissions)
