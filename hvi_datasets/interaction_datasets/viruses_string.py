from typing import Tuple, List

import pandas as pd

from hvi_datasets.taxonomy import Taxonomy
from hvi_datasets.dataset_base_classes import DatasetHVI
from hvi_datasets.dataset_base_classes import DatasetHVIStandardized


class DatasetVirusesString(DatasetHVI):
    """
    Viruses.STRING: A Virus-Host Protein-Protein Interaction Database
    (Cook et al. 2018, https://doi.org/10.3390/v10100519)

    URL: http://viruses.string-db.org/download/protein.links.detailed.v10.5/9606.protein.links.detailed.v10.5.txt.gz

    Fields from CSV:
    protein1
    protein2
    neighborhood
    fusion
    cooccurence
    coexpression
    experimental: > 0 -> Experimental evidence
    database
    textmining
    combined_score
    """

    delimiter = " "
    name = "viruses_string"
    header = "infer"

    def get_lyssa_associated_interactions(self):
        return DatasetVirusesString(data_frame=self.data_frame[self.data_frame["protein2"].str.contains("11292.")])

    def get_sars_cov_2_associated_interactions(self):
        return DatasetVirusesString(data_frame=self.data_frame[self.data_frame["protein2"].str.contains("2697049.")])

    def get_rhabdoviridae_associated_interactions(self):
        """
        No viral family included in the Viruses.STRING format. To get interactions associated with a certain viral
        family, convert the dataset to a standardized dataset while providing a taxonomy object.

        :raise: ValueError
        """
        raise ValueError

    def get_coronaviridae_associated_interactions(self):
        """
        No viral family included in the Viruses.STRING format. To get interactions associated with a certain viral
        family, convert the dataset to a standardized dataset while providing a taxonomy object.

        :raise: ValueError
        """
        raise ValueError

    def map_to_uniprot_id(self, human_virus_interactions: pd.DataFrame) -> Tuple[pd.DataFrame, List[int]]:
        """
        Maps the protein identifiers in the human-virus interaction dataset to UniProt IDs.
        ID mapping files have been retrieved via https://www.uniprot.org/id-mapping.

        :param human_virus_interactions: The human-virus interaction dataset with protein identifiers.
        :return: Tuple containing the mapped human-virus interaction dataset with UniProt IDs and a list of unknown ids.
        """
        mapping_to_uniprot = {}
        mapping_files = ["../raw_data/viruses_string/ensp_to_uniprot.tsv",
                         "../raw_data/viruses_string/uniprotacid_to_uniprot.tsv"]
        for file_path in mapping_files:
            with open(file_path, "r") as mapping_file:
                lines = mapping_file.readlines()
                for line in lines[1:]:
                    values = line.split("\t")
                    mapping_to_uniprot[values[0]] = values[1]

        unknown_ids = []
        new_ids_human = []
        new_ids_virus = []
        for idx, (id_human, id_virus) in enumerate(
                zip(human_virus_interactions["protein1"], human_virus_interactions["protein2"])):
            try:
                uniprot_human = mapping_to_uniprot[id_human.split(".")[1]].replace("\n", "")
                uniprot_virus = mapping_to_uniprot[id_virus.split(".")[1]].replace("\n", "")
                new_ids_human.append(uniprot_human)
                new_ids_virus.append(uniprot_virus)
            except KeyError:
                unknown_ids.append(idx)

        if len(unknown_ids) > 0:
            human_virus_interactions = human_virus_interactions.drop(unknown_ids).reset_index(drop=True)
            print(f"Dropped {len(unknown_ids)} unknown ids from {self.name} dataset!")

        human_virus_interactions["protein1"] = pd.Series(new_ids_human)
        human_virus_interactions["protein2"] = pd.Series(new_ids_virus)
        return human_virus_interactions, unknown_ids

    def to_standardized_dataset(self, taxonomy: Taxonomy) -> DatasetHVIStandardized:
        """
        Converts the Viruses.String dataset to the standardized format.
        Drops all non-human-virus interactions and maps non-uniprot ids to uniprot ids where possible.

        :param taxonomy: Previously created taxonomy object
        :return: Standardized Viruses.String dataset
        """

        human_virus_interactions = self.data_frame[~self.data_frame["protein2"].str.contains("9606.")].reset_index(
            drop=True)

        taxon_virus = []
        for value in human_virus_interactions["protein2"].values:
            taxon = value.split(".")[0]
            taxon_virus.append(taxon)

        human_virus_interactions, unknown_ids = self.map_to_uniprot_id(human_virus_interactions)
        series_taxon_virus = pd.Series(taxon_virus).drop(unknown_ids).reset_index(drop=True)
        assert len(human_virus_interactions) == len(series_taxon_virus)

        series_protein_human = human_virus_interactions["protein1"]
        series_protein_virus = human_virus_interactions["protein2"]

        series_dataset_name = pd.Series([self.name] * len(series_protein_human))

        series_experimental = human_virus_interactions["experimental"] > 0

        series_target = pd.Series(["1"] * len(series_protein_human))

        series = {"Uniprot_human": series_protein_human,
                  "Uniprot_virus": series_protein_virus,
                  "Taxon_virus": series_taxon_virus,
                  "Dataset": series_dataset_name,
                  "Experimental": series_experimental,
                  "Target": series_target}

        standardized_df = pd.DataFrame(series)
        hvi_std = DatasetHVIStandardized(data_frame=standardized_df, name=self.name, taxonomy=taxonomy)
        return hvi_std
