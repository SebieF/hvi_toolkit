from __future__ import annotations

import pandas as pd

from .hvi_abstract_dataset import DatasetHVI
from .hvi_standardized_dataset import DatasetHVIStandardized

from ..taxonomy import Taxonomy


class DatasetPPIStandardized(DatasetHVI):
    """
        Dataset containing standardized protein-protein interaction features.

        Fields from CSV:
            interactor1:        Uniprot ID for protein 1, e.g. Q6P0Q8 [single]
            interactor2:        Uniprot ID for protein 2, e.g. Q8JJY9 [single]
            taxon1:             Taxon ID from NCBI, e.g. 11292 for Rabies Lyssavirus, 2697049 for Sars-Cov-2 [multiple]
            taxon2:             Taxon ID from NCBI, e.g. 11292 for Rabies Lyssavirus, 2697049 for Sars-Cov-2 [multiple]
            interacting:        1 if positive interaction else 0
            experimental_score: MI-Score for interaction (0-1)
    """

    @staticmethod
    def get_header() -> str:
        return "interactor1,interactor2,taxon1,taxon2,interacting,experimental_score"

    @classmethod
    def from_hvi_standardized(cls, hvi_standardized: DatasetHVIStandardized) -> DatasetPPIStandardized:
        hvi_length = len(hvi_standardized)
        ppi_data_frame: pd.DataFrame = pd.DataFrame()
        ppi_data_frame["interactor1"] = hvi_standardized.data_frame["Uniprot_human"]
        ppi_data_frame["interactor2"] = hvi_standardized.data_frame["Uniprot_virus"]
        ppi_data_frame["taxon1"] = ["9606"] * hvi_length
        ppi_data_frame["taxon2"] = hvi_standardized.data_frame["Taxon_virus"]
        ppi_data_frame["interacting"] = hvi_standardized.data_frame["Target"]
        ppi_data_frame["experimental_score"] = [0.0] * hvi_length

        # Only allow one taxon value for ppi set at the moment
        ppi_data_frame["taxon2"] = ppi_data_frame["taxon2"].apply(lambda v: v.split(",")[0])
        return cls(data_frame=ppi_data_frame)

    def to_standardized_dataset(self, taxonomy: Taxonomy):
        ppi_length = len(self)
        hvi_data_frame: pd.DataFrame = pd.DataFrame()
        hvi_data_frame["Uniprot_human"] = self.data_frame["interactor1"]
        hvi_data_frame["Uniprot_virus"] = self.data_frame["interactor2"]
        hvi_data_frame["Taxon_virus"] = self.data_frame["taxon2"]
        hvi_data_frame["Dataset"] = ["ppi"] * ppi_length
        hvi_data_frame["Experimental"] = self.data_frame["experimental_score"] > 0.6
        hvi_data_frame["Target"] = self.data_frame["interacting"]

        return DatasetHVIStandardized(data_frame=hvi_data_frame, taxonomy=taxonomy)


