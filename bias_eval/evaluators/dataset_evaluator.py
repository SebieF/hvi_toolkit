import colorama

from colorama import Fore
from typing import List, Dict, Any
from scipy.stats import pearsonr, ttest_ind, chisquare
from biotrainer.utilities import read_FASTA, get_attributes_from_seqrecords_for_protein_interactions, get_split_lists, \
    INTERACTION_INDICATOR

from bias_eval.dataset_splitting import SplitsGenerator
from bias_eval.utilities import Interaction
from bias_eval.dataset_base_classes import DatasetHVIStandardized
from .metrics_calculator import calculate_all_metrics

colorama.init(autoreset=True)


class DatasetEvaluator:

    def __init__(self, significance: float = 0.05, hub_threshold: int = 5, bias_threshold: float = 0.9):
        """
        Create DatasetEvaluator object.

        The following custom parameters for checks can be set:
        :param significance: Sets the level of significance for statistical test (chi-squared, t-test..) [0-1]
        :param hub_threshold: Sets the threshold starting from how many associated interactions a protein will
                              be regarded as a hub protein [>=2]
        :param bias_threshold: Sets the threshold for when to accept the bias in the dataset [0-1]
        """
        assert 0 < significance < 1.0, f"Significance must be between 0.0 and 1.0!"
        assert hub_threshold >= 2, f"Hub threshold must be greater or equal to 2!"
        assert 0 < bias_threshold < 1.0, f"Bias threshold must be between 0.0 and 1.0!"
        self.significance = significance
        self.hub_threshold = hub_threshold
        self.bias_threshold = bias_threshold

    @staticmethod
    def _get_interaction_id(interaction: Interaction) -> str:
        """
        :return: An interaction id (protein_human&protein_virus) from an interaction object
        """
        return f"{interaction.uniprot_human}{INTERACTION_INDICATOR}{interaction.uniprot_virus}"

    @staticmethod
    def _calculate_dataset_bias(dataset: List[Interaction]) -> Any:
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
        test_statistic_bias = abs(test_statistic_bias)

        def bias_predictor(interactor1, interactor2):
            pos_occurrences_total = positive_counts[interactor1]
            pos_occurrences_total += positive_counts[interactor2]
            neg_occurrences_total = negative_counts[interactor1]
            neg_occurrences_total += negative_counts[interactor2]

            return 0 if neg_occurrences_total >= pos_occurrences_total else 1

        fore = Fore.BLACK
        if test_statistic_bias < 0.90:
            fore = Fore.RED
        print(f"\n**DATASET BIAS:**")
        print(f"Correlation between negative and positive interactions: {fore}{test_statistic_bias} "
              f"(p-value: {p_value_bias})")
        if test_statistic_bias < 0.90:
            print(f"Your dataset bias is quite high. The bias means that the frequency for proteins associated with \n"
                  f"interactions differs substantially between the negative and positive set. \n"
                  f"This will probably result in a high bias predictor baseline. \n"
                  f"To improve, you can remove hub proteins for human and viral proteins. \n"
                  f"Alternatively, you need to re-consider the creation of your negative dataset and \n"
                  f"account for the characteristics of your ppi network.")
        else:
            print(f"Your dataset bias is rather small. This indicates that your dataset represents the underlying "
                  f"ppi network well.")
        return bias_predictor

    @staticmethod
    def _calculate_bias_predictions(bias_predictor: Any, test_dataset: List[Interaction], name: str):
        """
        Calculates and reports classification metrics for the test_dataset from the bias predictor

        :param bias_predictor: Bias predictor (function from bias_baseline)
        :param test_dataset: Dataset to evaluate
        """

        # Predict all test_samples from bias
        predictions = []
        test_set_targets = []
        for test_sample in test_dataset:
            interactor1 = test_sample.uniprot_human
            interactor2 = test_sample.uniprot_virus
            predictions.append(bias_predictor(interactor1, interactor2))
            test_set_targets.append(int(test_sample.target))

        # Calculate metrics for bias predictions
        bias_metrics = calculate_all_metrics(predictions=predictions,
                                             targets=test_set_targets,
                                             df_name=name)

        print(f"\n**Bias baseline predictions for {name}:**")
        print(bias_metrics)

    def _check_sequence_lengths(self, interactions: List[Interaction], id2seq: Dict[str, str]):
        """
        Check if sequence lengths differ significantly between positive and negative interactions.

        :param interactions: Interactions to check
        :param id2seq: Dictionary that maps protein ids to sequences
        """
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
        ttest_stat, ttest_p_value = ttest_ind(positive_sequence_lengths, negative_sequence_lengths)
        if ttest_p_value < self.significance:
            print("The lengths of your positive and negative interactions differ significantly!\n"
                  "This might introduce bias to the model's predictions if it can infer the sequence lengths somehow.")
        else:
            print("The lengths of your positive and negative interactions do not differ significantly.\n"
                  "The dataset is well balanced in this regard.")
        print(f"T-test result: {ttest_stat} (p-value: {ttest_p_value})")

    def _check_uniform_distribution(self, category_frequencies: Dict[str, int], name: str):
        """
        Checks if the frequencies given are uniformly distributed via Chi-Squared test.

        :param category_frequencies: Frequencies how often a key appears in the category
        :name: Name of the category
        """

        maximum_category = max(category_frequencies.items(), key=lambda k: k[1])
        minimum_category = min(category_frequencies.items(), key=lambda k: k[1])
        test_statistic_chi2, p_value_chi2 = chisquare(f_obs=list(category_frequencies.values()))

        print(f"\n**Uniform distribution: {name}**")
        fore = Fore.BLACK
        if p_value_chi2 < self.significance:
            fore = Fore.RED
            print(f"{name} category is not uniformly distributed within your dataset. This might be a reason for "
                  f"biased interactions!")
        else:
            print(f"{name} category is uniformly distributed within your dataset. The dataset is, hence, balanced "
                  f"in this regard.")
        print(f"Most frequent: {maximum_category}")
        print(f"Least frequent: {minimum_category}")
        print(f"{fore}Chi-squared test: {test_statistic_chi2} (p-value: {p_value_chi2})")

    def _check_protein_hubs(self, interactions: List[Interaction]):
        """
        Check if the dataset contains hub proteins (and hence hub interactions). The threshold above which
        a protein is classified as a hub protein is defined by the constructor (Default: 5).

        :param interactions: List of interactions
        """
        protein_hub_interactions, _ = SplitsGenerator.get_protein_hubs(interactions=interactions,
                                                                       hub_threshold=self.hub_threshold)
        print(f"\n**Protein Hubs:**")
        if len(protein_hub_interactions) == 0:
            print(f"No protein hub interactions in dataset.")
        else:
            print(f"Number of hub interactions: {len(protein_hub_interactions)}")

        return protein_hub_interactions

    def evaluate_standardized_dataset(self, standardized_dataset: DatasetHVIStandardized,
                                      sequences_fasta_path: str = ""):
        """
        Performs the following checks on a standardized_dataset object:

        1. Check dataset bias: Do proteins associated with interactions have the same frequency in the positive
        and negative dataset?
        2. Based on the dataset bias - Metrics are given for predicting the interactions from frequencies:
        2.a) For the whole dataset
        2.b) For hub interactions (if existent) - Number of hub interactions is reported additionally
        3. Viral families: Are they uniformly distributed across the interactions?

        If sequences are given via fasta file:
        4. Sequence length: Does it differ significantly between positive and negative interactions

        :param standardized_dataset: DatasetHVIStandardized object
        :param sequences_fasta_path: Path to sequences fasta file (optional, enables check for sequence lengths)
        """
        interaction_list: List[Interaction] = standardized_dataset.to_interaction_list()
        print(f"**Bias evaluation for dataset with {len(interaction_list)} interactions:**")

        # 1. Dataset bias
        bias_predictor = self._calculate_dataset_bias(dataset=interaction_list)

        # 2.a) Bias predictor metrics for whole dataset
        self._calculate_bias_predictions(bias_predictor=bias_predictor, test_dataset=interaction_list,
                                         name="Whole dataset")

        # 2.b) Bias predictor for hub interactions
        protein_hub_interactions = self._check_protein_hubs(interactions=interaction_list)
        if len(protein_hub_interactions) > 0:
            self._calculate_bias_predictions(bias_predictor=bias_predictor, test_dataset=protein_hub_interactions,
                                             name="Hub interactions only")

        # 3. Viral families
        viral_families = {}
        for interaction in interaction_list:
            if interaction.family_virus not in viral_families:
                viral_families[interaction.family_virus] = 0
            viral_families[interaction.family_virus] += 1
        self._check_uniform_distribution(category_frequencies=viral_families, name="Viral families")

        # 4. Sequence lengths
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

            self._check_sequence_lengths(interactions=interaction_list, id2seq=id2seq)

    def evaluate_biotrainer_fasta(self, biotrainer_fasta_path: str):
        """
        Performs the following checks on a biotrainer fasta file for interactions:

        1. Check dataset bias: Do proteins associated with interactions have the same frequency in the positive
        and negative dataset?
        2. Based on the dataset bias - Metrics are given for predicting the interactions from frequencies:
        2.a) For the whole dataset
        2.b) For the training dataset (if available)
        2.c) For the validation dataset (if available)
        2.d) For the test dataset (if available)
        2.e) For hub interactions (if existent) - Number of hub interactions is reported additionally
        3. Sequence length: Does it differ significantly between positive and negative interactions

        :param biotrainer_fasta_path: Path to biotrainer fasta file (suited for interaction_mode)
        """
        seq_records = read_FASTA(biotrainer_fasta_path)
        id2seq = {seq.id: seq for seq in seq_records}
        id2attributes = get_attributes_from_seqrecords_for_protein_interactions(seq_records)
        train, val, test = tuple(map(set, get_split_lists(id2attributes)))

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

        print(f"**Bias evaluation for dataset with {len(interaction_list)} interactions:**")
        # 1. Dataset bias
        bias_predictor = self._calculate_dataset_bias(dataset=interaction_list)

        # 2.a) Bias predictor metrics for whole dataset
        self._calculate_bias_predictions(bias_predictor=bias_predictor, test_dataset=interaction_list,
                                         name="Whole dataset")

        # 2.b) Bias predictor metrics for training dataset
        if len(train_interactions) > 0:
            self._calculate_bias_predictions(bias_predictor=bias_predictor, test_dataset=train_interactions,
                                             name="Training set only")
        # 2.c) Bias predictor metrics for validation dataset
        if len(val_interactions) > 0:
            self._calculate_bias_predictions(bias_predictor=bias_predictor, test_dataset=val_interactions,
                                             name="Validation set only")
        # 2.d) Bias predictor metrics for test dataset
        if len(test_interactions) > 0:
            self._calculate_bias_predictions(bias_predictor=bias_predictor, test_dataset=test_interactions,
                                             name="Test set only")
        # 2.e) Bias predictor metrics for hub interactions
        protein_hub_interactions = self._check_protein_hubs(interactions=interaction_list)
        if len(protein_hub_interactions) > 0:
            self._calculate_bias_predictions(bias_predictor=bias_predictor, test_dataset=protein_hub_interactions,
                                             name="Hub interactions only")

        # 3. Check sequence lengths
        self._check_sequence_lengths(interactions=interaction_list, id2seq=id2seq)
