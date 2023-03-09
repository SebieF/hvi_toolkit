from __future__ import annotations

import pandas as pd

from Bio import SeqIO
from tqdm import tqdm
from typing import List, Tuple


class DatasetHVIStdExtraFunctionality:
    """
    Extra functionality for analysing human-virus interactions for the DatasetHVIStandardized class
    """

    def plot_uniques_by_datasets(self, experimental: bool = False):
        """
        Creates a barplot that displays the amount of unique interactions contributed to the whole dataset by
        each data source. So if the data in the standardized dataset comes from two different datasets,
        e.g. HVIDB and VirusMenta, two columns with the respective numbers of interactions will be displayed.

        :param experimental: Include only experimentally verified interactions in the plot
        """

        std_database = self.get_only_experimental_values() if experimental else self
        data_frame = std_database.data_frame
        list_of_names = []
        for datasets in data_frame["Dataset"].values:
            for dataset in datasets.split(","):
                list_of_names.append(dataset)
        unique_names = set(list_of_names)
        unique_interactions = {unique_name: [std_database.get_number_of_uniques_by_dataset(unique_name)] for
                               unique_name in unique_names}
        print(f"Unique interactions in dataset (experimental={experimental})")
        print(unique_interactions)
        df_unique_interactions = pd.DataFrame.from_dict(data=unique_interactions)
        df_unique_interactions.plot(kind="bar", title="Unique interactions by dataset")

    def get_number_of_uniques_by_dataset(self, dataset_name: str) -> int:
        """
        Returns number of entries that only have the specified dataset_name, i.e. they uniquely come from this
        data source.

        :param dataset_name: Name of dataset to examine
        :return Number of entries that only come from dataset_name
        """
        counter = 0
        for dataset in self.data_frame["Dataset"].values:
            if dataset == dataset_name:
                counter += 1
        return counter

    def get_only_n_most_common_within_category(self, category: str, n: int = 10) -> DatasetHVIStdExtraFunctionality:
        """
        Returns a new DatasetHVIStdExtraFunctionality object containing only the rows of the original
        dataset that belong to the n most common values within the specified category.

        :param category: A string representing the name of the category column to filter on.
        :param n: An integer representing the number of most common values to select.
        :return: A new DatasetHVIStdExtraFunctionality object containing the filtered dataset.
        """

        # Split comma separated values:
        name_list = []
        for names in self.data_frame[category].values:
            for name in names.split(","):
                name_list.append(name)

        counts = pd.Series(name_list).value_counts(sort=True)
        category_names = counts[0:n].index.tolist()
        selected_dataframe = self.data_frame[self.data_frame[category].isin(category_names)]

        return self.__class__(data_frame=selected_dataframe)

    def remove_unavailable_and_overlong_sequences(self, fasta_file_path: str,
                                                  seq_len_threshold: int = 10000) -> \
            Tuple[str, DatasetHVIStdExtraFunctionality]:
        """
        This function takes a path to a FASTA file and a threshold sequence length.
        It removes all sequences from the dataset that are not present in the FASTA file
        or whose length is greater than the threshold.
        The FASTA file should not contain any duplicates.

        :param fasta_file_path: The path to a FASTA file containing the sequences to filter against
        :param seq_len_threshold: Integer representing the maximum allowed length of sequences. Default: 10,000
        :return: A Tuple with the fasta string without the unavailable or overlong sequences and a
        new DatasetHVIStdExtraFunctionality object with the filtered data
        :except AssertionError if the fasta file contains duplicates
        """

        unique_ids = self.get_unique_proteins()

        fasta_file = list(SeqIO.parse(fasta_file_path, "fasta"))
        ids_from_fasta = [seq.id for seq in fasta_file]
        ids_to_sequence = {seq.id: seq.seq for seq in fasta_file}
        assert len(ids_from_fasta) == len(set(ids_from_fasta)), "Fasta contains duplicates!"

        # Remove not matching dataset-fasta (=> No sequence data)
        non_matching = []
        for unique_id in unique_ids:
            if unique_id not in ids_from_fasta:
                non_matching.append(unique_id)
        print(f"Removing {len(non_matching)} non-matching ids: {non_matching}")
        dropped_df = self.data_frame[~self.data_frame["Uniprot_human"].isin(non_matching)]
        dropped_df = dropped_df[~dropped_df["Uniprot_virus"].isin(non_matching)]

        # Remove sequences > seq_len_threshold
        overlong = []
        for seq in fasta_file:
            if len(seq.seq) > seq_len_threshold:
                overlong.append(seq.id)
        print(f"Removing {len(overlong)} overlong ids: {overlong}")
        dropped_df = dropped_df[~dropped_df["Uniprot_human"].isin(overlong)]
        dropped_df = dropped_df[~dropped_df["Uniprot_virus"].isin(overlong)]

        standardized_dataset_dropped = self.__class__(data_frame=dropped_df)

        return DatasetHVIStdExtraFunctionality(data_frame=self.data_frame)

    def _calculate_interaction_values_by_heuristic(self, interactor_id_list: List[str]) -> List[Tuple[str, int]]:
        """
        Calculates the value of each interactor in the given list of interactors based on their number of interactions,
        the virus family they belong to, the dataset they are sourced from, and whether they represent
        experimental data or not. The values of interactors in the same list are
        returned as a list of tuples (interactor ID, interaction value).

        :param interactor_id_list: List of interactor IDs to calculate interaction values for
        :return: List of tuples representing each interactor in interactor_id_list with their corresponding interaction
                 value
        :rtype: List[Tuple[str, int]]
        """

        # Assign value to every interactor_id
        value_rhabdoviridae = -2  # Rhabdoviridae get moved to test set later anyway
        value_coronaviridae = -2  # Coronaviridae get moved to test set later anyway
        value_virus_string = -2  # Worst and overrepresented datasource
        value_experimental = 1  # Increase relative number of experimental interactions

        # Clusters with only one interactor do not need to be evaluated
        if len(interactor_id_list) == 1:
            return [(interactor_id_list[0], -1)]

        interactors_with_value: List[(str, int)] = []
        for interactor in interactor_id_list:
            value = 0
            all_interactions = self.data_frame[
                self.data_frame["Uniprot_human"].isin([interactor]) |
                self.data_frame["Uniprot_virus"].isin([interactor])
                ]
            value += len(all_interactions)
            for _, interaction in all_interactions.iterrows():
                if interaction["Family_virus"] == "Rhabdoviridae":
                    value += value_rhabdoviridae
                if interaction["Family_virus"] == "Coronaviridae":
                    value += value_coronaviridae
                if "viruses_string" in str(interaction["Dataset"]) and len(str(interaction["Dataset"]).split(",")) == 1:
                    value += value_virus_string
                if interaction["Experimental"] == "True":
                    value += value_experimental
            interactors_with_value.append((interactor, value))
        return interactors_with_value

    def drop_interactions_by_clustering(self, clusters: List[List[str]], cluster_name: str, mode="heuristic") \
            -> Tuple[List[str], DatasetHVIStdExtraFunctionality]:
        """
        Drops interactions from the data frame based on the given clustering and mode. If mode is "heuristic", the
        interaction values for each interactor in the given clusters are calculated using
        _calculate_interaction_values_by_heuristic(), and the interactor with the highest value is kept, while
        other interactors are dropped. If mode is "naive", the first interactor in each cluster is kept, and others are
        dropped. The interactors that have been dropped are returned as a list, along with the updated data frame.

        :param clusters: List of clusters, where each cluster is a list of interactor IDs
        :param cluster_name: Name of the clustering algorithm used to create the given clusters,
                             e.g. "sequence identity" or "embeddings distance"
        :param mode: Clustering mode, either "heuristic" or "naive". If mode is "heuristic", the function calculates the
                     interaction values for each interactor in clusters using
                     _calculate_interaction_values_by_heuristic() function. This means that the "cluster representative"
                     is not necessarily kept!
                     If mode is "naive", the first interactor in each cluster is kept, and the rest is dropped.
                     Default: "heuristic".
        :return: A tuple containing a list of IDs that have been dropped and the updated dataset
        :rtype: Tuple[List[str], DatasetHVIStdExtraFunctionality]
        """

        ids_to_drop = []

        if mode == "heuristic":
            print("Calculating interaction values by heuristic!")
            for cluster_ids in tqdm(clusters):
                interactors_with_value = self._calculate_interaction_values_by_heuristic(cluster_ids)
                if len(interactors_with_value) > 1:
                    max_value = max(interactors_with_value, key=lambda item: item[1])
                    ids_to_drop.extend([interactor_with_value[0] for interactor_with_value in interactors_with_value
                                        if interactor_with_value != max_value])
        else:  # Naive: First value == cluster representative
            for cluster in clusters:
                if len(cluster) > 1:
                    ids_to_drop.extend(cluster[1:])

        length_before = len(self.data_frame)
        dropped_df = self.data_frame[~self.data_frame["Uniprot_human"].isin(ids_to_drop)]
        dropped_df = dropped_df[~dropped_df["Uniprot_virus"].isin(ids_to_drop)]
        length_after = len(dropped_df)
        print(f"Loss of interactions due to {cluster_name} clustering: \n"
              f"Before: {length_before}, After: {length_after}, Loss: {length_before - length_after} "
              f"({100 * (length_before - length_after) / length_before} %)")

        return ids_to_drop, self.__class__(data_frame=dropped_df)
