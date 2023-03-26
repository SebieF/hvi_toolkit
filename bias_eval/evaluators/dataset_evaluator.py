import colorama

from colorama import Fore
from typing import List, Dict, Any, Tuple
from biotrainer.utilities import read_FASTA, get_attributes_from_seqrecords_for_protein_interactions, get_split_lists, \
    INTERACTION_INDICATOR

from scipy.stats import pearsonr, ttest_ind, chisquare

from bias_eval.dataset_splitting import SplitsGenerator
from bias_eval.utilities import Interaction
from bias_eval.dataset_base_classes import DatasetHVIStandardized
from .metrics_calculator import calculate_all_metrics

colorama.init(autoreset=True)


class DatasetEvaluator:

    def _get_interaction_id(self, interaction: Interaction):
        return f"{interaction.uniprot_human}{INTERACTION_INDICATOR}{interaction.uniprot_virus}"

    def _calculate_dataset_bias(self, dataset: List[Interaction]) -> Any:
        """
            Calculates a dataset bias baseline for interactions (see for example
            Park, Marcotte 2011: https://doi.org/10.1093/bioinformatics/btr514).
            At first, it is counted how often a protein (id) is found in the positive and negative sets.
            Then, the dataset bias can be determined by calculating pearson-r for (positive_counts vs. negative counts).

            The predictor itself just sums up the positive and negative counts for both interactors and "predicts"
            the higher value.
        """

        # 1. Calculate protein counts
        positive_counts = {}
        negative_counts = {}

        for sample in dataset:
            interactor1 = sample.uniprot_human
            interactor2 = sample.uniprot_virus
            for count_dict in [positive_counts, negative_counts]:
                if interactor1 not in count_dict:
                    count_dict[interactor1] = 0
                if interactor2 not in count_dict:
                    count_dict[interactor2] = 0

            if int(sample.target) == 1:
                positive_counts[interactor1] += 1
                positive_counts[interactor2] += 1
            else:
                negative_counts[interactor1] += 1
                negative_counts[interactor2] += 1

        # 2. Calculate dataset bias
        test_statistic_bias, p_value_bias = pearsonr(list(positive_counts.values()), list(negative_counts.values()))

        def bias_predictor(interactor1, interactor2):
            pos_occurrences_total = positive_counts[interactor1]
            pos_occurrences_total += positive_counts[interactor2]
            neg_occurrences_total = negative_counts[interactor1]
            neg_occurrences_total += negative_counts[interactor2]

            return 0 if neg_occurrences_total >= pos_occurrences_total else 1

        fore = Fore.BLACK
        if abs(test_statistic_bias) < 0.90:
            fore = Fore.RED
        print(f"\n**DATASET BIAS:**")
        print(f"Correlation between negative and positive interactions: {fore}{abs(test_statistic_bias)} "
              f"(p-value: {p_value_bias})")
        return bias_predictor

    def _calculate_bias_predictions(self, bias_predictor: Any, test_dataset: List[Interaction], name: str):
        # 3. Predict all test_samples from bias
        predictions = []
        test_set_targets = []
        for test_sample in test_dataset:
            interactor1 = test_sample.uniprot_human
            interactor2 = test_sample.uniprot_virus
            predictions.append(bias_predictor(interactor1, interactor2))
            test_set_targets.append(int(test_sample.target))

        # 4. Calculate metrics for bias predictions
        bias_metrics = calculate_all_metrics(predictions=predictions,
                                             targets=test_set_targets,
                                             df_name="bias_predictions")

        print(f"Bias baseline predictions for {name}:")
        print(bias_metrics)

    def _calculate_sequence_lengths(self, interactions: List[Interaction], id2seq: Dict[str, str]):
        positive_interactions = [interaction for interaction in interactions if interaction.target == "1"]
        negative_interactions = [interaction for interaction in interactions if interaction.target == "0"]
        assert len(negative_interactions) == len(interactions) - len(positive_interactions), \
            f"Targets for some interactions missing!"

        positive_sequence_lengths = []
        negative_sequence_lengths = []
        for positive_interaction in positive_interactions:
            positive_sequence_lengths.append(len(id2seq[positive_interaction.uniprot_human]) +
                                             len(id2seq[positive_interaction.uniprot_virus]))
        for negative_interaction in negative_interactions:
            negative_sequence_lengths.append(len(id2seq[negative_interaction.uniprot_human]) +
                                             len(id2seq[negative_interaction.uniprot_virus]))

        average_length_positive = sum(positive_sequence_lengths) / len(positive_sequence_lengths) \
            if len(positive_sequence_lengths) > 0 else 0
        average_length_negative = sum(negative_sequence_lengths) / len(negative_sequence_lengths) \
            if len(negative_sequence_lengths) > 0 else 0
        print(f"\n**Sequence lengths:**")
        print(
            f"Average length positive: {average_length_positive} (# Positive: {len(positive_sequence_lengths)})")
        print(
            f"Average length negative: {average_length_negative} (# Negative: {len(negative_sequence_lengths)})")
        print(ttest_ind(positive_sequence_lengths, negative_sequence_lengths))

    def _check_uniform_distribution(self, category_frequencies: List, name: str, significance: float = 0.05):
        test_statistic_chi2, p_value_chi2 = chisquare(f_obs=category_frequencies)

        print(f"\n**Uniform distribution: {name}**")
        fore = Fore.BLACK
        if p_value_chi2 < significance:
            fore = Fore.RED
        print(f"{fore}Test statistic chi-squared test: {test_statistic_chi2} (p-value: {p_value_chi2})")

    def _check_protein_hubs(self, interactions: List[Interaction], hub_threshold: int = 2):
        protein_hub_interactions, _ = SplitsGenerator.get_protein_hubs(interactions=interactions,
                                                                       hub_threshold=hub_threshold)
        print(f"\n**Protein Hubs:**")
        if len(protein_hub_interactions) == 0:
            print(f"No protein hub interactions in dataset!")
        else:
            print(f"Number of hub interactions: {len(protein_hub_interactions)}")


    def evaluate_standardized_dataset(self, standardized_dataset: DatasetHVIStandardized,
                                      sequences_fasta_path: str = ""):
        interaction_list: List[Interaction] = standardized_dataset.to_interaction_list()

        if sequences_fasta_path != "":
            seq_records = read_FASTA(sequences_fasta_path)
            id2seq = {seq.id: seq.seq for seq in seq_records}

            # Check that all sequences have been provided
            missing_sequences = []
            for interaction in interaction_list:
                if interaction.uniprot_human not in id2seq.keys():
                    missing_sequences.append(interaction.uniprot_human)
                if interaction.uniprot_virus not in id2seq.keys():
                    missing_sequences.append(interaction.uniprot_virus)

            if len(missing_sequences) > 0:
                raise Exception(f"Provided fasta file does not contain sequences for all protein ids! \n"
                                f"Missing: {missing_sequences}")

            self._calculate_sequence_lengths(interactions=interaction_list, id2seq=id2seq)

        bias_predictor = self._calculate_dataset_bias(dataset=interaction_list)
        self._calculate_bias_predictions(bias_predictor=bias_predictor, test_dataset=interaction_list,
                                         name="Whole dataset")
        self._check_protein_hubs(interactions=interaction_list)

        # Viral families
        viral_families = {}
        for interaction in interaction_list:
            if interaction.family_virus not in viral_families:
                viral_families[interaction.family_virus] = 0
            viral_families[interaction.family_virus] += 1

        self._check_uniform_distribution(category_frequencies=list(viral_families.values()), name="Viral families")

    def evaluate_biotrainer_fasta(self, biotrainer_fasta_path: str):
        seq_records = read_FASTA(biotrainer_fasta_path)
        id2seq = {seq.id: seq for seq in seq_records}
        id2attributes = get_attributes_from_seqrecords_for_protein_interactions(seq_records)
        train, val, test = get_split_lists(id2attributes)

        interaction_list = []
        for interaction_id, attrs in id2attributes.items():
            uniprot_human = interaction_id.split(INTERACTION_INDICATOR)[0]
            uniprot_virus = interaction_id.split(INTERACTION_INDICATOR)[1]
            target = attrs["TARGET"]
            interaction = Interaction(uniprot_human=uniprot_human, uniprot_virus=uniprot_virus,
                                      family_virus="", experimental="",
                                      target=target)
            interaction_list.append(interaction)
        train_interactions = [interaction for interaction in interaction_list
                              if self._get_interaction_id(interaction) in train]
        val_interactions = [interaction for interaction in interaction_list
                            if self._get_interaction_id(interaction) in val]
        test_interactions = [interaction for interaction in interaction_list
                             if self._get_interaction_id(interaction) in test]

        self._calculate_sequence_lengths(interactions=interaction_list, id2seq=id2seq)

        bias_predictor = self._calculate_dataset_bias(dataset=interaction_list)
        self._calculate_bias_predictions(bias_predictor=bias_predictor, test_dataset=interaction_list,
                                         name="Whole dataset")

        self._check_protein_hubs(interactions=interaction_list)

        if len(train_interactions) > 0:
            self._calculate_bias_predictions(bias_predictor=bias_predictor, test_dataset=train_interactions,
                                             name="Training set only")
        if len(val_interactions) > 0:
            self._calculate_bias_predictions(bias_predictor=bias_predictor, test_dataset=val_interactions,
                                             name="Validation set only")
        if len(test_interactions) > 0:
            self._calculate_bias_predictions(bias_predictor=bias_predictor, test_dataset=val_interactions,
                                             name="Test set only")
