from Bio import SeqIO
from pathlib import Path
from abc import abstractmethod
from collections import namedtuple
from typing import List, Dict, Tuple, Set
from sklearn.model_selection import train_test_split
from hvi_datasets.dataset_base_classes import DatasetHVIStandardized

Interaction = namedtuple("Interaction", "uniprot_human uniprot_virus family_virus target")


class SplitsGenerator:
    """
    Abstract class to create dataset splits from standardized human-virus interaction datasets (DatasetHVIStandardized).
    Train/Val/Test splits are generated and written to a biotrainer-readable format. It can then be used with
    the interaction mode of biotrainer.

    Also provides get_protein_hubs and get_interactions_by_families methods which might be useful outside of this class.
    """

    def __init__(self, interaction_dataset: List[Interaction],
                 id2seq: Dict[str, str]):
        self.interaction_dataset = interaction_dataset
        self.id2seq = id2seq

    @classmethod
    def from_standardized_dataset(cls, std_dataset: DatasetHVIStandardized, full_fasta_file: str = ""):
        """
        Creates a split generator object from a standardized human-virus interaction dataset.

        :param std_dataset: Dataset with human-virus interactions
        :param full_fasta_file: Optional fasta file with all sequences for the proteins in the dataset.
                                If the fasta file is not given, dummy sequences are written to the biotrainer file.
        :return:
        """
        # Create list of interactions (compatible with sklearn train_test_split)
        interaction_dataset = []
        for _, row in std_dataset.data_frame.iterrows():
            interaction = Interaction(uniprot_human=row["Uniprot_human"], uniprot_virus=row["Uniprot_virus"],
                                      family_virus=row["Family_virus"], target=row["Target"])
            interaction_dataset.append(interaction)
        # Create seqs
        if full_fasta_file != "":
            id2seq = {seq.id: seq.seq for seq in list(SeqIO.parse(full_fasta_file, "fasta"))}
        else:
            # Dummy sequences if no fasta file is given
            id2seq = {seq_id: "SEQ" for seq_id in std_dataset.get_unique_proteins()}
        return cls(interaction_dataset=interaction_dataset, id2seq=id2seq)

    def _write_to_biotrainer_format(self, train: List[Interaction], val: List[Interaction], test: List[Interaction],
                                    output_path: Path):
        """
        Returns a biotrainer compatible sequences.fasta file that can used with the biotrainer interaction mode.
        Interactions in the file are doubled, i.e. they are stored in both directions in order to be able
        to save the sequences for all proteins in a fasta-compatible way.
        Theoretically, for each different protein, the sequence would only have to be stored once and biotrainer
        would still identify all interactions correctly. However, it is simpler like this and probably also easier
        to comprehend at first glance.

        :param train: Train split
        :param val: Validation split
        :param test: Test split
        :param output_path: Path to store the file at
        """

        with open(output_path, "w") as biotrainer_output_file:
            for name, split in [("train", train), ("val", val), ("test", test)]:
                for sample in split:
                    id_human = sample.uniprot_human
                    id_virus = sample.uniprot_virus
                    target = sample.target
                    set_annotation = f"SET={name} "
                    sequence_human = self.id2seq[id_human]
                    sequence_virus = self.id2seq[id_virus]

                    biotrainer_output_file.write(f">{id_human} INTERACTOR={id_virus} "
                                                 f"TARGET={target} {set_annotation}\n")
                    biotrainer_output_file.write(f"{sequence_human}\n")
                    biotrainer_output_file.write(f">{id_virus} INTERACTOR={id_human} "
                                                 f"TARGET={target} {set_annotation}\n")
                    biotrainer_output_file.write(f"{sequence_virus}\n")
        print(f"Wrote splits to {output_path} biotrainer file!")

    @staticmethod
    def get_protein_hubs(interactions: List[Interaction],
                         hub_threshold: int = 5) -> Tuple[List[Interaction], Set[str]]:
        """
        Get proteins in interaction dataset, that interact with more than hub_threshold other proteins
        Only looks at human proteins!

        :param interactions: List of interactions
        :param hub_threshold: Threshold to determine protein hubs (# Interactions >= hub_threshold => Protein is hub)
        :return Tuple with list of protein hub interactions and a set with the ids of all hub proteins
        """

        protein_counter_dict = {}
        for interaction in interactions:
            protein_counter_dict[interaction.uniprot_human] = []
            protein_counter_dict[interaction.uniprot_virus] = []

        for interaction in interactions:
            protein_counter_dict[interaction.uniprot_human].append(interaction)

        protein_hub_interactions = {}
        protein_hub_ids = set()
        for protein_id, associated_interactions in protein_counter_dict.items():
            if len(associated_interactions) >= hub_threshold:
                protein_hub_ids.add(protein_id)
                for interaction in associated_interactions:
                    protein_hub_interactions[f"{interaction.uniprot_human}&{interaction.uniprot_virus}"] = interaction

        return list(protein_hub_interactions.values()), protein_hub_ids

    @staticmethod
    def _split_by_hubs(interactions: List[Interaction]) -> \
            Tuple[List[Interaction], List[Interaction]]:
        """
        Split given interactions by hubs (proteins that interact with multiple other proteins)
        Only looks at human proteins!

        :param interactions: List of protein_protein interactions
        :return Tuple with remaining_interactions (without hubs), protein hub interactions
        """
        length_with_hubs = len(interactions)
        protein_hub_interactions, protein_hub_ids = SplitsGenerator.get_protein_hubs(interactions)

        remaining_interactions = []
        remove_from_hub_interactions = []
        added_protein_hubs = set()
        for interaction in interactions:
            if interaction.uniprot_human not in protein_hub_ids:
                remaining_interactions.append(interaction)
            else:
                if interaction.uniprot_human not in added_protein_hubs:
                    remaining_interactions.append(interaction)  # Add one interaction per hub to remaining interactions
                    remove_from_hub_interactions.append(interaction)
                    added_protein_hubs.add(interaction.uniprot_human)

        # Remove those that are still in the training set from the test/hub set
        for interaction in remove_from_hub_interactions:
            protein_hub_interactions.remove(interaction)

        assert (len(protein_hub_interactions) + len(remaining_interactions)) == length_with_hubs, \
            f"Dataset sizes after split do not match!"
        print(f"Removed hubs from training set: \n"
              f"Number of protein hubs: {len(added_protein_hubs)}\n"
              f"Before: {length_with_hubs}, After: {len(remaining_interactions)}, "
              f"Loss: {length_with_hubs - len(remaining_interactions)} "
              f"({100 * (length_with_hubs - len(remaining_interactions)) / length_with_hubs} %)")
        return remaining_interactions, protein_hub_interactions

    @staticmethod
    def get_interactions_by_families(interactions: List[Interaction], viral_families: List[str]) -> List[Interaction]:
        """
        Returns all interactions that are associated with one of the given viral_families.

        :param interactions: All interactions to search
        :param viral_families: List of viral families for which interactions should be filtered
        :return: All interactions associated with given viral families
        """

        interactions_with_viral_families = []
        for interaction in interactions:
            if interaction.family_virus in viral_families:
                interactions_with_viral_families.append(interaction)

        return interactions_with_viral_families

    @staticmethod
    def _split_by_families(interactions: List[Interaction],
                           viral_families_to_remove: List[str]) -> Tuple[List[Interaction], List[Interaction]]:
        """
        Split given interactions by their viral family

        :param interactions: List of interactions
        :param viral_families_to_remove: List of viral family names, e.g. ["Rhabdoviridae", "Coronaviridae"]
        :return Tuple of interactions without_viral_families and interactions with_viral_families
        """

        print(f"Removing viral families: {viral_families_to_remove}")

        length_with_all_families = len(interactions)

        with_viral_families = SplitsGenerator.get_interactions_by_families(interactions,
                                                                           viral_families_to_remove)

        without_viral_families = [interaction for interaction in interactions
                                  if interaction not in with_viral_families]

        assert (len(with_viral_families) + len(without_viral_families)) == length_with_all_families, \
            f"Dataset sizes after split do not match!"
        print(f"Split interactions by {', '.join(viral_families_to_remove)} viral families!\n"
              f"Length with viral families:\t{len(with_viral_families)}\n"
              f"Length without viral families:\t{len(without_viral_families)}")

        return without_viral_families, with_viral_families

    @staticmethod
    def _validate_split(train: List[Interaction], val: List[Interaction], test: List[Interaction],
                         initial_dataset_length: int):
        """
        Validate the given train, val, test split
        1. Check that the summed length of all splits equals the dataset size before splitting
        2. Check that the validation and test sets do not contain any samples from the train set

        :param train: Train split
        :param val: Validation split
        :param test: Test split
        :param initial_dataset_length: Length of dataset before splitting
        :raises AssertError if any of the check fails
        """

        assert len(train) + len(val) + len(test) == initial_dataset_length
        for train_sample in train:
            assert train_sample not in val, f"Training sample in validation set!"
            assert train_sample not in test, f"Training sample in test set!"

    @staticmethod
    def _print_split_statistics(train: List[Interaction], val: List[Interaction], test: List[Interaction]):
        combined_length = len(train) + len(val) + len(test)
        print(f"Number train: {len(train)} ({100 * len(train) / combined_length} %)")
        print(f"Number val: {len(val)} ({100 * len(val) / combined_length} %)")
        print(f"Number test: {len(test)} ({100 * len(test) / combined_length} %)")

        positives_train = [sample for sample in train if sample.target == "1"]
        positives_val = [sample for sample in val if sample.target == "1"]
        positives_test = [sample for sample in test if sample.target == "1"]
        print(f"Positives in train set: {len(positives_train)} ({100 * len(positives_train) / max(1, len(train))} %)")
        print(f"Positives in validation set: {len(positives_val)} ({100 * len(positives_val) / max(1, len(val))} %)")
        print(f"Positives in test set: {len(positives_test)} ({100 * len(positives_test) / max(1, len(test))} %)")

    @staticmethod
    def _create_output_file_name(split_name: str, viral_families_to_remove: List[str],
                                 random_seed: int, split_hubs: bool, output_path: str) -> str:
        output_path = output_path if output_path[-1] == "/" else output_path + "/"
        output_file_path = output_path + "biotrainer_" + split_name + f"_{random_seed}"
        output_file_path += f"_splithubs_{'True' if split_hubs else 'False'}"
        output_file_path += f"_removed_{'_'.join([family[0:4] for family in viral_families_to_remove])}"
        output_file_path += ".fasta"
        return output_file_path

    @abstractmethod
    def create_splits(self, split_name: str,
                      viral_families_to_remove: List[str], random_seed: int, split_hubs: bool, output_path: str):
        raise NotImplementedError


class SplitGeneratorRandomHoldOut(SplitsGenerator):
    split_sizes = {"train": 0.85,
                   "val": 0.10,
                   "test": 0.05}

    def create_splits(self, split_name: str,
                      viral_families_to_remove: List[str], random_seed: int, split_hubs: bool, output_path: str):
        """
        Create randomized splits for hold out cross validation.
        Split sizes are determined by self.split_sizes.
        Creates a train, validation and test set.

        :param split_name: Name for the splits
        :param viral_families_to_remove: Viral families defined in this list are put to the test set before splitting
        :param random_seed: Random seed for sklearn
        :param split_hubs: If True, all hub protein interactions are moved to the test set
        :param output_path: Path to write the resulting biotrainer file to
        """

        current_dataset = self.interaction_dataset
        if len(viral_families_to_remove) > 0:
            current_dataset, _ = self._split_by_families(current_dataset, viral_families_to_remove)

        print(f"Creating splits {self.split_sizes}")

        combined_dataset_length = len(current_dataset)
        targets = [interaction.target for interaction in current_dataset]
        number_positives = sum(map(int, targets))
        number_negatives = combined_dataset_length - number_positives

        print(f"Number positives: {number_positives}")
        print(f"Number negatives: {number_negatives}")

        print(
            f"Number combined: {combined_dataset_length} (Pos: {100 * number_positives / combined_dataset_length} %, "
            f"Neg: {100 * number_negatives / combined_dataset_length} %)")
        # Split:
        train, test = train_test_split(current_dataset, test_size=self.split_sizes["test"],
                                       random_state=random_seed, stratify=targets)

        if split_hubs:
            train, test = self._split_by_hubs(train + test)
        targets_train = [interaction.target for interaction in train]

        # Validation size must be corrected by value that has been taken from train set already
        validation_size = self.split_sizes["val"] + self.split_sizes["test"] / 10
        train, val = train_test_split(train, test_size=validation_size,
                                      random_state=random_seed, stratify=targets_train)

        self._validate_split(train, val, test, initial_dataset_length=combined_dataset_length)
        self._print_split_statistics(train, val, test)

        output_file_name = self._create_output_file_name(split_name, viral_families_to_remove, random_seed, split_hubs,
                                                         output_path)
        self._write_to_biotrainer_format(train, val, test, Path(output_file_name))
        print("Finished creating splits!")


class SplitGeneratorCrossValidation(SplitsGenerator):

    def create_splits(self, split_name: str, viral_families_to_remove: List[str], random_seed: int, split_hubs: bool,
                      output_path: str):
        """
        Create splits suited for k-fold or leave_p_out cross validation.
        Creates only a train/validation and test set.
        If no viral families are given and split_hubs == False, then no test set can be created!

        :param split_name: Name for the splits
        :param viral_families_to_remove: Viral families defined in this list are put to the test set before splitting
        :param random_seed: Random seed for sklearn
        :param split_hubs: If True, all hub protein interactions are moved to the test set
        :param output_path: Path to write the resulting biotrainer file to
        """
        if not split_hubs:
            assert len(viral_families_to_remove) > 0, "Viral_families_to_remove are required if split_hubs == False!"

        current_dataset = self.interaction_dataset
        initial_dataset_length = len(current_dataset)
        print(f"Splitting dataset for biotrainer cross validation! (Total length: {initial_dataset_length})")

        train_val, test = self._split_by_families(current_dataset, viral_families_to_remove)

        if split_hubs:
            train_val, protein_hubs = self._split_by_hubs(train_val)
            test.extend(protein_hubs)

        self._validate_split(train_val, [], test, initial_dataset_length)
        self._print_split_statistics(train_val, [], test)

        output_file_name = self._create_output_file_name(split_name, viral_families_to_remove, random_seed, split_hubs,
                                                         output_path)
        self._write_to_biotrainer_format(train_val, [], test, Path(output_file_name))
        print("Finished creating splits!")


def main():
    configs = [
        {
            "name": "merged_dataset_combined",
            "splits_generator": SplitGeneratorCrossValidation.from_standardized_dataset(
                std_dataset=DatasetHVIStandardized(
                    file_path="../merged_dataset/redundancy_reduction/merged_dataset_complete_rr.csv"
                ),
                full_fasta_file="../merged_dataset/merged_dataset_complete.fasta",
            ),
            "viral_families_to_remove": ["Rhabdoviridae", "Coronaviridae"],
            "random_seed": 42,
            "split_hubs": True,
            "output_path": "../../model/",
            "execute": True
        },
        {
            "name": "merged_dataset_combined_non_reduced",
            "splits_generator": SplitGeneratorCrossValidation.from_standardized_dataset(
                std_dataset=DatasetHVIStandardized(
                    file_path="../merged_dataset/merged_dataset_complete.csv"
                ),
                full_fasta_file="",  # Ignore sequences, biotrainer does not need them if embeddings are provided
            ),
            "viral_families_to_remove": ["Rhabdoviridae", "Coronaviridae"],
            "random_seed": 42,
            "split_hubs": False,
            "output_path": "../merged_dataset/non_reduced/",
            "execute": False
        }
    ]
    for config in configs:
        if config["execute"]:
            print(f"Starting splits generation for {config['name']}!")
            config["splits_generator"].create_splits(split_name=config["name"],
                                                     viral_families_to_remove=config["viral_families_to_remove"],
                                                     random_seed=config["random_seed"],
                                                     split_hubs=config["split_hubs"],
                                                     output_path=config["output_path"]
                                                     )


if __name__ == "__main__":
    main()
