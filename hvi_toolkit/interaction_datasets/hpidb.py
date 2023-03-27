from hvi_toolkit.dataset_base_classes import DatasetMITAB25


class DatasetHPIDBVirus(DatasetMITAB25):
    """
    Human-Pathogen-Interaction DataBase 3.0, successor of HPIDB 2.0: a curated database for host–pathogen interactions
    (Ammari et al. 2016, https://doi.org/10.1093/database/baw103)

    URL: https://hpidb.igbb.msstate.edu/keyword.html
    - Query type: Taxon Name / Species
    - Query: Virus
    - Do not include Interlog Predictions

    **Gets filtered to only contain entries with UniprotID and score values!**

    For column descriptions, see https://psicquic.github.io/MITAB25Format.html
    """

    delimiter = "\t"
    name = "hpidb_virus"
    header = None
