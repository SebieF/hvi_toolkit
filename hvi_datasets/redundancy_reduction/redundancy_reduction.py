import h5py
import pathlib
import tempfile
import subprocess
import numpy as np
import pandas as pd

from Bio import SeqIO
from tqdm import tqdm
from collections import namedtuple
from itertools import combinations
from typing import List, Dict, Tuple, Any
from hvi_datasets.dataset_base_classes import DatasetHVIStandardized

DistanceValues = namedtuple("DistanceValues", "min max avg")


class RedundancyReduction:
    """
    Class that handles the redundancy reduction of a standardized interaction dataset

    Steps performed:
    1. Filter by mmseqs2: Sequence identity clusters
    2. Reduce via euclidean distance on embeddings

    Note that step 1 requires an active installation of mmseqs in a conda environment called "mmseqs2":
    ```
    conda create --name mmseqs2
    conda activate mmseqs2
    conda install -c conda-forge -c bioconda mmseqs2
    ```

    :var embeddings_threshold: If the distance between two embeddings is below this threshold, one must be dropped
    :var mmseqs2_cluster_mode: Clusters via this mmseqs2 command:
        --min-seq-id 0.25: All proteins in the cluster have to have an identity of at least 0.25
        -c 0.75: Coverage > 75%
        --cov-mode 0: Bidirectional (recommended for full length protein sequences)

    :param dataset: Dataset to perform redundancy reduction for
    :param fasta_file: Fasta file with sequences for all proteins in the dataset
    :param embeddings_file: h5 file with embeddings for all proteins in the dataset
    :param output_path: Output path were (intermediate) results of redundancy reduction should be saved
    """

    embeddings_threshold = 0.6
    mmseqs2_cluster_mode = " --min-seq-id 0.25 -c 0.75 --cov-mode 0"

    def __init__(self, dataset: DatasetHVIStandardized, fasta_file: str, embeddings_file: h5py.File,
                 output_path: str):
        self.dataset = dataset
        self.fasta_file = fasta_file
        self.id2emb_sequence = {embeddings_file[idx].attrs["original_id"]: np.array(embedding) for (idx, embedding) in
                                embeddings_file.items()}
        self.output_path = output_path if output_path[-1] == "/" else output_path + "/"

    def _run_mmseqs(self, fasta_file: str, name: str) -> pd.DataFrame:
        """
        Runs mmseqs2 to cluster similar sequences in a fasta file and returns the results as a DataFrame.
        To do this, a conda environment called "mmseqs2" with mmseqs installed has to be created beforehand.

        :param fasta_file: Fasta file with protein ids and sequences
        :param name: Name of this run
        :return: pandas DataFrame with cluster Result from mmseqs2
        """
        output_dir = pathlib.Path(f"{self.output_path}clusterRes/")
        if not output_dir.is_dir():
            print(f"Creating output dir: {output_dir}")
            output_dir.mkdir(parents=False)

        run_name = f"{str(output_dir)}/clusterRes_{name}"
        with tempfile.TemporaryDirectory() as tmp_dir_name:
            command = f"conda run -n mmseqs2 mmseqs easy-cluster " \
                      f"{fasta_file} {run_name} {tmp_dir_name} {self.mmseqs2_cluster_mode}"
            print(f"Running mmseqs via {command}")
            process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
            process.wait()
            if process.returncode == 0:
                print("Successfully clustered via mmseqs2!")
            else:
                print("mmseqs2 failed, exiting...")
                exit(process.returncode)
        mmseqs_result = pd.read_csv(f"{run_name}_cluster.tsv", delimiter="\t", header=None)

        return mmseqs_result

    @staticmethod
    def _mmseqs_result_to_cluster_dict(mmseqs_result: pd.DataFrame) -> Dict[str, List[str]]:
        """
        This method takes the output from mmseqs2 clustering and converts it to a dictionary,
        where keys represent cluster IDs and values are the list of redundant IDs.

        :param mmseqs_result: DataFrame containing the output of the mmseqs2 clustering (from _run_mmseqs()).
                              This DataFrame should have two columns: the first column should be the cluster ID,
                              and the second column should be the redundant IDs.
        :return: A dictionary where the keys are cluster IDs and the values are lists of
                 redundant IDs that belong to each cluster.
        """
        cluster_dict = {}
        for _, row in mmseqs_result.iterrows():
            cluster_id = row[0]
            redundant_id = row[1]
            if cluster_id not in cluster_dict.keys():
                cluster_dict[cluster_id] = []
            cluster_dict[cluster_id].append(redundant_id)
        return cluster_dict

    def _distance_result_to_cluster_dict(self, distance_dict: Dict) -> Dict[str, List[str]]:
        """
        Converts the distance_dict from _calculate_euclidean_distances() to a cluster dict that contains
        a list of redundant sequence ids for every cluster representative.

        :param distance_dict: A dict containing the euclidean distance between every pair of embeddings in the dataset
        :return: A dict with the sequence IDs as keys and the list of redundant sequences as values
        """

        cluster_dict = {}
        all_clustered_keys = []
        for key_left, dict_right in distance_dict.items():
            if key_left in all_clustered_keys:
                continue
            if key_left not in cluster_dict.keys():
                cluster_dict[key_left] = [key_left]
            for key_right, distance in dict_right.items():
                if distance < self.embeddings_threshold:
                    cluster_dict[key_left].append(key_right)
                    all_clustered_keys.append(key_right)
        return cluster_dict

    def _create_fasta_from_dataset(self, dataset: DatasetHVIStandardized, name: str):
        """
        Creates a new fasta file from self.fasta_file that only includes proteins that are still left in
        the dataset after a step of redundancy reduction.

        :param dataset: Dataset after redundancy reduction step
        :param name: Name for output fasta file
        :return: Output file path of the created fasta file
        """

        seq_records = SeqIO.parse(self.fasta_file, "fasta")
        output_file_path = f"{self.output_path}{name}.fasta"
        print(f"Creating fasta file: {output_file_path}")
        unique_protein_ids = dataset.get_unique_proteins()
        with open(output_file_path, "w") as output_file:
            for seq in seq_records:
                if seq.id in unique_protein_ids:
                    output_file.write(f">{seq.id}\n")
                    output_file.write(f"{seq.seq}\n")
        return output_file_path

    def _sequence_identity_reduction(self) -> DatasetHVIStandardized:
        """
        Reduces sequence redundancy in the dataset by using mmseqs2 to cluster similar sequences
        and remove redundant ones.

        Returns a new `DatasetHVIStandardized` object that contains the reduced set of non-redundant sequences
        and interactions.

        Steps:
        1. Run mmseqs2 on all sequences in the dataset.
        2. Drop interactions that belong to redundant sequences based on the clusters obtained from step 1.
        (mode: "heuristic" = value every protein by certain criteria, see _calculate_interaction_values_by_heuristic)
        3. Run mmseqs2 again on the reduced dataset to obtain final clusters of non-redundant sequences.
        4. Drop interactions that belong to redundant sequences based on the clusters obtained from step 3.
        (mode: "naive" = drop all but the cluster representative)
        :return: A new `DatasetHVIStandardized` object containing the reduced set of non-redundant sequences.
        """

        # 1. First run on all sequences
        mmseqs_result = self._run_mmseqs(self.fasta_file, name="First")
        mmseqs_cluster_dict = self._mmseqs_result_to_cluster_dict(mmseqs_result)

        # 2. Drop interactions (mode: "heuristic")
        _, mmseqs_filtered_dataset = self.dataset.drop_interactions_by_clustering(
            clusters=list(mmseqs_cluster_dict.values()),
            cluster_name="mmseqs2-First",
            mode="heuristic"
        )

        # 3. Run again because of dropping by heuristic
        reduced_fasta_file = self._create_fasta_from_dataset(mmseqs_filtered_dataset, name="reduced_first_run")
        mmseqs_result_second = self._run_mmseqs(reduced_fasta_file, name="Second")

        # 4. Dop interactions (mode: "naive")
        mmseqs_cluster_dict_second = self._mmseqs_result_to_cluster_dict(mmseqs_result_second)
        _, mmseqs_filtered_dataset_second = self.dataset.drop_interactions_by_clustering(
            clusters=list(mmseqs_cluster_dict_second.values()),
            cluster_name="mmseqs2-Second",
            mode="naive"
        )

        # Create second reduced fasta file
        _ = self._create_fasta_from_dataset(mmseqs_filtered_dataset_second,
                                            name="reduced_second_run")
        return mmseqs_filtered_dataset_second

    def _cluster_embeddings_by_taxon(self, dataset: DatasetHVIStandardized) -> Tuple[Dict[str, Any], Dict[str, Any]]:
        """
        Divide self.id2emb_sequence into embeddings that belong to human and viral proteins.

        :param dataset: Interaction dataset to retrieve human and viral protein uniprot ids from
        :return: Tuple with human_embeddings and viral_embeddings
        """

        human_embeddings = {}
        viral_embeddings = {}
        unique_human_ids = set(dataset.data_frame["Uniprot_human"].unique().tolist())
        unique_viral_ids = set(dataset.data_frame["Uniprot_virus"].unique().tolist())
        for seq_id in unique_human_ids:
            human_embeddings[seq_id] = self.id2emb_sequence[seq_id]
        for seq_id in unique_viral_ids:
            viral_embeddings[seq_id] = self.id2emb_sequence[seq_id]

        assert len(human_embeddings.keys()) + len(viral_embeddings.keys()) <= len(
            dataset.data_frame["Uniprot_human"].unique()) + len(
            dataset.data_frame["Uniprot_virus"].unique())
        assert len(human_embeddings.keys()) <= len(unique_human_ids)
        assert len(viral_embeddings.keys()) <= len(unique_viral_ids)

        return human_embeddings, viral_embeddings

    def _calculate_euclidean_distances(self, sequence_pairs: List,
                                       name: str = "") -> Tuple[DistanceValues, Dict[str, Dict[str, float]]]:
        """
        Calculates the euclidean distance between the embeddings of all given sequence pairs
        
        :param sequence_pairs: List of pair-wise sequences
        :param name: Name of this round of euclidean clustering
        :return: Tuple with all calculated distance metrics and a distance dict that contains the distance for all pairs
        """

        distance_min = np.Inf
        distance_max = 0
        distance_sum = 0
        distance_dict = {}
        print(f"Calculating {name} euclidean distances..")
        for seq_id1, seq_id2 in sequence_pairs:
            distance_dict[seq_id1] = {}
            distance_dict[seq_id2] = {}

        for seq_id1, seq_id2 in tqdm(sequence_pairs):
            embedding1 = self.id2emb_sequence[seq_id1]
            embedding2 = self.id2emb_sequence[seq_id2]
            euclidean_distance = np.linalg.norm(embedding1 - embedding2)
            distance_min = min(distance_min, euclidean_distance)
            distance_max = max(distance_max, euclidean_distance)
            distance_sum += euclidean_distance

            if euclidean_distance < self.embeddings_threshold:
                distance_dict[seq_id1][seq_id2] = euclidean_distance
                distance_dict[seq_id2][seq_id1] = euclidean_distance
        return DistanceValues(distance_min, distance_max, distance_sum/len(sequence_pairs)), distance_dict

    @staticmethod
    def _print_distances(distance_values: DistanceValues, name: str = ""):
        print(f"Average distance {name}: {distance_values.avg}")
        print(f"Max distance {name}: {distance_values.max}")
        print(f"Min distance {name}: {distance_values.min}")

    @staticmethod
    def _drop_embeddings_by_ids(ids_to_drop: List[str], embeddings: Dict[str, Any], name: str = "") -> Dict[str, Any]:
        """
        Drop ids_to_drop from given embeddings

        :param ids_to_drop: Sequence ids to drop
        :param embeddings: Embeddings to drop sequences from
        :param name: Name of this step of redundancy reduction
        :return: Reduced embeddings without sequences from ids_to_drop
        """

        reduced_distance_embeddings = {seq_id: embedding for seq_id, embedding in embeddings.items()
                                       if seq_id not in ids_to_drop}
        print(f"Dropped {len(embeddings) - len(reduced_distance_embeddings)} sequences from {name} embeddings!")
        return reduced_distance_embeddings

    def _validate_embeddings_after_drop(self, sequence_ids: List[str], name: str = ""):
        """
        Checks that the minimum of all euclidean distances is actually bigger than the given self.embeddings_threshold
        for all pairs of sequences
        
        :param sequence_ids: Sequence ids to check embedding distances for 
        :param name: Name of this redundancy reduction step 
        """

        sequence_pairs = list(combinations(sequence_ids, 2))
        distance_values, _ = self._calculate_euclidean_distances(sequence_pairs, name=name)
        self._print_distances(distance_values=distance_values, name=name)
        assert distance_values.min > self.embeddings_threshold

    def _euclidean_distance_reduction(self, mmseqs_filtered_dataset: DatasetHVIStandardized) -> DatasetHVIStandardized:
        """
        Redundancy reduce sequences in dataset via euclidean distance.
        Redundancy reduction is conducted separately for human and viral embeddings

        :param mmseqs_filtered_dataset: Dataset after filtering for sequence identity via mmseqs
        :return: Interaction dataset without redundant sequences (regarding their embeddings)
        """

        current_dataset = mmseqs_filtered_dataset
        for name in ["Viral", "Viral-Second", "Human", "Human-Second"]:
            embeddings = self._cluster_embeddings_by_taxon(current_dataset)
            current_embeddings = embeddings[0] if "Human" in name else embeddings[1]

            # Calculate euclidean distances
            sequence_pairs = list(combinations(current_embeddings.keys(), 2))
            distance_values, distance_dict = self._calculate_euclidean_distances(sequence_pairs, name=name)
            self._print_distances(distance_values, name=name)

            euclidean_cluster_dict = self._distance_result_to_cluster_dict(distance_dict=distance_dict)

            # Drop redundant sequences by mode
            mode = "naive" if "Second" in name else "heuristic"  # 1. Heuristic 2. Naive (if necessary)
            ids_to_drop, updated_dataset = current_dataset.drop_interactions_by_clustering(
                clusters=list(euclidean_cluster_dict.values()),
                cluster_name="euclidean distance",
                mode=mode
            )

            del distance_values
            del distance_dict
            del euclidean_cluster_dict
            del sequence_pairs  # Necessary because memory allocation can get too big

            reduced_embeddings = self._drop_embeddings_by_ids(ids_to_drop, current_embeddings, name=name)
            if "Second" in name:
                self._validate_embeddings_after_drop(list(reduced_embeddings.keys()), name=name)
            current_dataset = updated_dataset

        return current_dataset

    def redundancy_reduction(self) -> DatasetHVIStandardized:
        """
        Perform redundancy reduction on the parameters given in the constructor.

        Steps performed:
        1. Filter by mmseqs2: Sequence identity clusters
        2. Reduce via euclidean distance on embeddings

        Note that step 1 requires an active installation of mmseqs in a conda environment called "mmseqs2":
        ```
        conda create --name mmseqs2
        conda activate mmseqs2
        conda install -c conda-forge -c bioconda mmseqs2
        ```

        :return: Redundancy reduced interaction dataset
        """

        mmseqs_filtered_dataset = self._sequence_identity_reduction()
        euclidean_filtered_dataset = self._euclidean_distance_reduction(mmseqs_filtered_dataset)
        return euclidean_filtered_dataset


def main():
    configs = [
        {
            "name": "merged_dataset_complete_rr",
            "csv_path": "../merged_dataset/merged_dataset_complete.csv",
            "fasta_path": "../merged_dataset/merged_dataset_complete.fasta",
            "sequence_embeddings_path": "../merged_dataset/embeddings_prottrans/reduced_embeddings_file.h5",
            "output_path": "../merged_dataset/redundancy_reduction/",
            "execute": True
        },
    ]
    for config in configs:
        if config["execute"]:
            print(f"Starting redundancy reduction for {config['name']}!")
            # Load dataset:
            merged_dataset_complete = DatasetHVIStandardized(file_path=config["csv_path"])
            # Load embeddings:
            embeddings_file = h5py.File(config["sequence_embeddings_path"], 'r')

            redundancy_reduction = RedundancyReduction(merged_dataset_complete,
                                                       fasta_file=config["fasta_path"],
                                                       embeddings_file=embeddings_file,
                                                       output_path=config["output_path"])
            redundancy_reduced_dataset = redundancy_reduction.redundancy_reduction()

            redundancy_reduced_dataset.store(f"{config['output_path']}{config['name']}.csv")
            print(f"Finished redundancy reduction for {config['name']}!")


if __name__ == "__main__":
    main()
