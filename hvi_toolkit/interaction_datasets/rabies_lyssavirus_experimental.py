import pandas as pd

from hvi_toolkit.taxonomy import Taxonomy
from hvi_toolkit.dataset_base_classes import DatasetHVI
from hvi_toolkit.dataset_base_classes import DatasetHVIStandardized


class DatasetRabiesLyssavirusExperimental(DatasetHVI):
    """
    Rabies Lyssavirus interactions that have been extracted from Zandi et al. 2021
    (https://doi.org/10.52547%2Fibj.25.4.226).

    The dataset can be found under /lyssa_predict/rabies_lyssavirus_predictions/lyssa_experimental_interactions.csv

    * Protein: Name of the Rabies Lyssavirus protein gene (G, N, M, L, P)
    * Lyssavirus: Species of observed interaction
    * Taxon_virus: Taxon of species
    * Uniprot_virus: UniprotKB ID of Lyssavirus protein
    * Uniprot_human: UniprotKB ID of human protein
    * Host protein interactor: Name of the host (human) molecule
    * Method of PPI detection and references: Method how ppi was confirmed and references (see paper above)
    * Viral_family: Viral family of Lyssavirus, i.e. Rhabdoviridae for all entries
    """

    @staticmethod
    def get_header() -> str:
        return ("Protein ,Lyssavirus ,Taxon_virus,Uniprot_virus,Uniprot_human,Host protein interactor ,Method of PPI "
                "detection and reference(s) ,Viral_family")

    delimiter = ","
    name = "rabies_lyssavirus_experimental"
    header = "infer"

    def get_lyssa_associated_interactions(self) -> DatasetHVI:
        """
        All interactions in this dataset are associated with rabies lyssavirus

        :return: Unchanged dataset
        """
        return self

    def get_sars_cov_2_associated_interactions(self) -> DatasetHVI:
        """
        No sars_cov_2 interactions in this dataset

        :raise: ValueError
        """
        raise ValueError

    def get_rhabdoviridae_associated_interactions(self) -> DatasetHVI:
        """
        All interactions in this dataset are associated with the Rhabdoviridae viral family

        :return: Unchanged dataset
        """
        return self

    def get_coronaviridae_associated_interactions(self) -> DatasetHVI:
        """
        No coronaviridae interactions in this dataset

        :raise: ValueError
        """
        raise ValueError

    def to_standardized_dataset(self, taxonomy: Taxonomy):
        series_protein_human = self.data_frame["Uniprot_human"]
        series_protein_virus = self.data_frame["Uniprot_virus"]

        series_dataset_name = pd.Series([self.name] * len(series_protein_human))
        series_taxon_virus = self.data_frame["Taxon_virus"]
        series_family_virus = self.data_frame["Viral_family"]

        series_experimental = pd.Series([True] * len(series_protein_human))

        series_target = pd.Series(["1"] * len(self.data_frame))

        series = {"Uniprot_human": series_protein_human,
                  "Uniprot_virus": series_protein_virus,
                  "Taxon_virus": series_taxon_virus,
                  "Family_virus": series_family_virus,
                  "Dataset": series_dataset_name,
                  "Experimental": series_experimental,
                  "Target": series_target}

        standardized_df = pd.DataFrame(series)
        hvi_std = DatasetHVIStandardized(data_frame=standardized_df, name=self.name, taxonomy=None)
        return hvi_std
