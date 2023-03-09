from __future__ import annotations

import pandas as pd

from taxonomy import Taxonomy
from typing import List, Optional
from hvi_abstract_dataset import DatasetHVI
from hvi_extra_functionality import DatasetHVIStdExtraFunctionality


class DatasetHVIStandardized(DatasetHVI, DatasetHVIStdExtraFunctionality):
    """
        Dataset containing standardized manually selected fields and methods relevant for HVI prediction
        Uniprot_human -> Uniprot_virus is ensured to be unique. Duplicated uniprot_ids that map to the same
        protein are not checked, however.

        This file contains the most important functionality for most hvi prediction tasks (specifically predicting
        Rabies lyssavirus and Sars-Cov-2).
        Additional functionality that was used during creation and analysis of the datasets
        can be found in DatasetHVIStandardizedExtraFunctionality.

        Fields from CSV:
            Uniprot_human:  Uniprot ID for human protein, e.g. Q6P0Q8 [single]
            Uniprot_virus:  Uniprot ID for viral protein, e.g. Q8JJY9 [single]
            Taxon_virus:    Taxon ID from NCBI, e.g. 11292 for Rabies Lyssavirus, 2697049 for Sars-Cov-2 [multiple]
            Family_virus:   Family name of virus, e.g. Coronaviridae / Rhabdoviridae [multiple]
            Dataset:        Name of database source (in this file) [multiple]
            Experimental:   If interaction was experimentally verified, True/False
                HVIDB:          False (no score given)
                VirusMentha:    Score > 0.6 (strict, see paper/letter VirusMINT)
                VirusHostNet:   MI-Score > 0.6 (strict, see paper)
                SarsIM28880:    MI-Score > 0.6 (strict, see paper)
                Viruses.String: experimental > 0
    """
    name = "standardized:"
    delimiter = ","
    header = 0

    def __init__(self, file_path: str = None, data_frame: pd.DataFrame = None, name: str = None,
                 taxonomy: Taxonomy = None):
        super().__init__(file_path=file_path, data_frame=data_frame)
        if name:
            self.name += name
        self._postprocess_std_dataset(taxonomy=taxonomy)
        self._assert_correctness()

    def _postprocess_std_dataset(self, taxonomy: Optional[Taxonomy] = None):
        """
        Post-processes the standardized dataset to prepare it for downstream use.

        The method ensures that all columns are cast to string.
        Additionally, it removes human-human interactions and duplicates from the dataset.
        Finally, adds viral families via the taxonomy object to the data_frame as a column, if a taxonomy object
        is given.

        :param taxonomy: Previously created taxonomy object or None if viral families should not get re-calculated
        """

        # All columns must be string
        for key in self.data_frame.keys():
            self.data_frame[key] = self.data_frame[key].astype("string")

        # Remove human-human
        self._remove_human_human()
        # Remove duplicates
        self._merge_and_remove_duplicates()
        # Add viral families
        if taxonomy:
            self._add_viral_family(taxonomy=taxonomy)

    def _remove_human_human(self):
        """
        Removes all human-human interactions from the current dataset.

        Part of post-processing the standardized dataset.
        """

        drop_human_human = []
        for idx, taxon in enumerate(self.data_frame["Taxon_virus"].values):
            if taxon == "9606":  # 9606 == taxon id of homo sapiens
                drop_human_human.append(idx)

        if len(drop_human_human) > 0:
            self.data_frame.drop(drop_human_human, inplace=True).reset_index(drop=True, inplace=True)
            print(f"Dropped {len(drop_human_human)} human-human interactions from {self.name}! "
                  f"Ids[0-5]: {drop_human_human[:5]}")

    def _merge_and_remove_duplicates(self):
        """
        Merge duplicates in the dataset based on Uniprot IDs of human and virus proteins.
        Duplicates are identified by the same combination of human and virus Uniprot IDs.
        After merging, duplicates are removed and dropped from the dataset.

        Part of post-processing the standardized dataset.
        """

        # Create a dictionary to map uniprot_human-uniprot_virus pairs to interaction indexes (human&virus)
        mapping = {}
        for idx, uniprot_human in enumerate(self.data_frame["Uniprot_human"].values):
            mapping[uniprot_human + "&" + self.data_frame["Uniprot_virus"].iloc[idx]] = []

        for idx, uniprot_human in enumerate(self.data_frame["Uniprot_human"].values):
            mapping[uniprot_human + "&" + self.data_frame["Uniprot_virus"].iloc[idx]].append(idx)

        duplicates = []
        for interaction_key, indexes in mapping.items():
            if len(indexes) > 1:
                keep = indexes[0]
                for remove in indexes[1:]:
                    duplicates.append((keep, remove))

        if len(duplicates) > 0:
            # Merge taxon IDs, database names:
            def merge_entries(database_key: str, keep_idx: int, remove_idx: int):
                entries_keep = self.data_frame[database_key].iloc[keep_idx].split(",")
                entries_remove = self.data_frame[database_key].iloc[remove_idx].split(",")
                for e_r in entries_remove:
                    if e_r not in entries_keep:
                        entries_keep.append(e_r)
                return entries_keep

            for keep, remove in duplicates:
                self.data_frame["Taxon_virus"].iloc[keep] = ",".join(merge_entries("Taxon_virus", keep, remove))
                self.data_frame["Dataset"].iloc[keep] = ",".join(merge_entries("Dataset", keep, remove))
                # Experimental: keep == True or remove == True
                self.data_frame["Experimental"].iloc[keep] = str(bool(self.data_frame["Experimental"].iloc[keep]) or
                                                                 bool(self.data_frame["Experimental"].iloc[remove]))
            self.data_frame.drop([remove for _, remove in duplicates], inplace=True).reset_index(drop=True,
                                                                                                 inplace=True)
            print(f"Dropped {len(duplicates)} duplicated interactions from {self.name}! "
                  f"Ids[0-5]: {duplicates[:5]}")

    def _add_viral_family(self, taxonomy: Taxonomy):
        """
        Adds the viral family column calculated from the Taxon_virus column using the given taxonomy object.
        Note that there might be multiple taxons for the same ID and hence also multiple families.

        Part of post-processing the standardized dataset.

        :param taxonomy: Previously created taxonomy object
        """
        viral_families = []
        for idx, taxons in enumerate(self.data_frame["Taxon_virus"].values):
            families = []
            for taxon in taxons.split(","):
                family = taxonomy.get_family_from_id(taxon)
                if family not in families:
                    families.append(family)
            viral_families.append(','.join(families))
        pd_series_family = pd.Series(data=viral_families)
        self.data_frame = self.data_frame.assign(Family_virus=pd_series_family)

    def _assert_correctness(self):
        """
        Check correctness of created standardized dataset. The following conditionals must be fulfilled:
        1. Interactions are unique
        2. "Taxon_virus", "Family_virus", "Dataset" columns only contain unique values
        3. All viral families are actually viral families ("-viridae" or "Unknown")
        """

        # 1. Check that interactions are unique
        mapping = {}
        for idx, uniprot_human in enumerate(self.data_frame["Uniprot_human"].values):
            mapping[uniprot_human + "-" + self.data_frame["Uniprot_virus"].iloc[idx]] = []

        for idx, uniprot_human in enumerate(self.data_frame["Uniprot_human"].values):
            mapping[uniprot_human + "-" + self.data_frame["Uniprot_virus"].iloc[idx]].append(idx)

        for k, val in mapping.items():
            assert len(val) == 1, f"Interaction {k} is not unique!"

        # 2. Check that additional columns only contain unique values
        columns_to_check = ["Taxon_virus", "Family_virus", "Dataset"]
        for column in columns_to_check:
            for idx, value in enumerate(self.data_frame[column].values):
                split_values = value.split(",")
                assert len(split_values) == len(set(split_values)), f"Column {column} contains non-unique" \
                                                                    f"values at index {idx}!"

        # 3. Check that only viral families are contained
        for idx, value in enumerate(self.data_frame["Family_virus"]):
            assert "viridae" in value or value == "Unknown", f"Non-viral family {value} at index {idx}!"

    def merge_with_standardized_dataset(self, other_dataset: DatasetHVIStandardized) -> DatasetHVIStandardized:
        """
        This method merges the current HVI standardized dataset with another HVI standardized dataset by concatenating
        them into a single dataset and then removing any duplicate rows.
        The resulting dataset has a new name which is a combination of the names of the original datasets.

        :param other_dataset: The other HVI standardized dataset to merge with.
        :returns: A new HVI standardized dataset that is the result of merging and deduplicating the two datasets.
        :rtype: DatasetHVIStandardized
        :raises ValueError: If `other_dataset` is not an instance of `DatasetHVIStandardized`.
        """
        merged_dataset = pd.concat([self.data_frame, other_dataset.data_frame], ignore_index=True)
        new_name = self.name + ";" + other_dataset.name

        hvi_std = DatasetHVIStandardized(data_frame=merged_dataset, name=new_name)
        return hvi_std

    def to_standardized_dataset(self, taxonomy: Taxonomy):
        return self

    def get_only_experimental_values(self) -> DatasetHVIStandardized:
        return DatasetHVIStandardized(data_frame=self.data_frame[self.data_frame["Experimental"].isin(["True"])])

    def get_only_positive_values(self) -> DatasetHVIStandardized:
        return DatasetHVIStandardized(data_frame=self.data_frame[self.data_frame["Target"].isin(["1"])])

    def get_only_negative_values(self) -> DatasetHVIStandardized:
        return DatasetHVIStandardized(data_frame=self.data_frame[self.data_frame["Target"].isin(["0"])])

    def get_lyssa_associated_interactions(self) -> DatasetHVIStandardized:
        return DatasetHVIStandardized(data_frame=self.data_frame[self.data_frame["Taxon_virus"].isin(["11292"])])

    def get_sars_cov_2_associated_interactions(self) -> DatasetHVIStandardized:
        return DatasetHVIStandardized(data_frame=self.data_frame[self.data_frame["Taxon_virus"].isin(["2697049"])])

    def get_rhabdoviridae_associated_interactions(self) -> DatasetHVIStandardized:
        return DatasetHVIStandardized(
            data_frame=self.data_frame[self.data_frame["Family_virus"].isin(["Rhabdoviridae"])])

    def get_coronaviridae_associated_interactions(self) -> DatasetHVIStandardized:
        return DatasetHVIStandardized(
            data_frame=self.data_frame[self.data_frame["Family_virus"].isin(["Coronaviridae"])])

    def get_interactions_by_condition(self, category: str, condition: str) -> DatasetHVIStandardized:
        return DatasetHVIStandardized(
            data_frame=self.data_frame[self.data_frame[category].isin([condition])])

    def get_interaction_ids(self, bi_directional: bool = False) -> List[str]:
        """
        Return all interactions in the dataset as a list of interaction ids.
        An interaction id is compiled as: Uniprot_human&Uniprot_virus

        :param bi_directional: Also include Uniprot_virus&Uniprot_human in the list
        :return: List of interaction ids
        """
        interaction_ids = []
        for _, row in self.data_frame.iterrows():
            interaction_id = f"{row['Uniprot_human']}{self.INTERACTION_INDICATOR}{row['Uniprot_virus']}"
            interaction_ids.append(interaction_id)
            if bi_directional:
                interaction_id_flipped = f"{row['Uniprot_virus']}{self.INTERACTION_INDICATOR}{row['Uniprot_human']}"
                interaction_ids.append(interaction_id_flipped)
        return interaction_ids

    def get_unique_proteins(self) -> List[str]:
        """
        Returns a list of unique uniprot identifiers of both human and viral proteins in the dataset.

        :return: List of unique protein uniprot ids
        """
        uniprot_ids = self.data_frame["Uniprot_human"].values.tolist()
        uniprot_ids.extend(self.data_frame["Uniprot_virus"].values.tolist())
        unique_ids = list(set(uniprot_ids))
        return unique_ids

    def download_sequences_from_uniprot(self) -> str:
        """
        Download sequences for all unique proteins in the dataset and store them in fasta format.
        As post-translational modifications (PTMs) have to be taken into account,
        https://www.ebi.ac.uk/proteins/api/proteins is used as a query database in this case.

        :return: Path to created fasta file (current_directory/downloaded_dataset_{self.name}.fasta)
        """
        unique_ids = self.get_unique_proteins()

        print(f"Downloading {len(unique_ids)} protein sequences!")
        sequence_dict = self.__get_sequences_with_ptm_information(uniprot_list=unique_ids)
        fasta_file_path = f"downloaded_dataset_{self.name}.fasta"
        with open(fasta_file_path, "w") as fasta_file:
            for key, val in sequence_dict.items():
                fasta_file.write(f">{key}\n")
                fasta_file.write(f"{val}\n")
        print(f"Successfully downloaded and stored {len(sequence_dict.keys())} sequences!")
        return fasta_file_path
