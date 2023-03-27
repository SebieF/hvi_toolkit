import torch
import random
import numpy as np

from .model_wrapper import ModelWrapper
from biotrainer.utilities import seed_all
from typing import List, Dict, Any, Tuple
from .metrics_calculator import calculate_all_metrics
from bio_embeddings.embed import ProtTransT5XLU50Embedder

amino_acids = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
interaction_operations = {
    "multiply": lambda embedding_left, embedding_right: torch.mul(embedding_left, embedding_right),
    "concat": lambda embedding_left, embedding_right: torch.concat([embedding_left, embedding_right])
}


class ModelEvaluator:

    def __init__(self, significance: float = 0.05, seed: int = 42):
        """
        Create ModelEvaluator object.

        The following custom parameters for checks can be set:
        :param significance: Sets the level of significance for statistical test (chi-squared, t-test..) [0-1]
        :param seed: If random sequences for sequence length benchmark are generated, this seed is used for
                     reproducibility.
        """
        assert 0 < significance < 1.0, f"Significance must be between 0.0 and 1.0!"
        self.significance = significance

        seed_all(seed)

    def _create_random_sequence_length_benchmark(self, model: ModelWrapper, interaction_mode: str):
        """
        This function creates random, adversarial sequences to check if the model tends to predict longer sequences
        as positive or vice versa. This requires a model trained on embeddings from the ProtTransT5XLU50Embedder
        at the moment. Embeddings are calculated on the fly, depending on the available resources this might take
        a while.

        The function serves as a prototype at the moment (work in progress).

        :param model: ModelWrapper object that supports from_embeddings(embeddings) -> list with predictions
        :param interaction_mode: Interaction mode ("multiply" or "concat") as defined in biotrainer. Decides which
                                 interaction operation is performed on the embeddings.
        """
        n_seqs = 1000
        longer_lengths = list(map(int, np.random.uniform(low=1000, high=2000, size=n_seqs // 2).tolist()))
        shorter_lengths = list(map(int, np.random.uniform(low=10, high=999, size=n_seqs // 2).tolist()))

        def create_random_sequence(sequence_length: int):
            seq = ""
            for residue in range(sequence_length):
                seq += random.choice(amino_acids)
            return seq

        longer_sequences = [create_random_sequence(seq_len) for seq_len in longer_lengths]
        shorter_sequences = [create_random_sequence(seq_len) for seq_len in shorter_lengths]

        embedder = ProtTransT5XLU50Embedder()
        embeddings_longer = [[embedder.reduce_per_protein(embedding)] for embedding
                             in list(embedder.embed_many(longer_sequences))]
        embeddings_shorter = [[embedder.reduce_per_protein(embedding)] for embedding
                              in list(embedder.embed_many(shorter_sequences))]

        interaction_operation = interaction_operations[interaction_mode]
        longer_interaction_embeddings = [interaction_operation(embeddings_longer[idx], embeddings_longer[-idx])
                                         for idx in range(len(embeddings_longer) // 2)]
        shorter_interaction_embeddings = [interaction_operation(embeddings_shorter[idx], embeddings_shorter[-idx])
                                          for idx in range(len(embeddings_shorter) // 2)]

        longer_predictions = model.from_embeddings(embeddings=longer_interaction_embeddings)
        shorter_predictions = model.from_embeddings(embeddings=shorter_interaction_embeddings)

        if sum(list(map(int, longer_predictions))) > sum(list(map(int, shorter_predictions))):
            print(f"Your model tends to predict longer sequences more often to be positive!")
        else:
            print(f"Your model tends to predict shorter sequences more often to be positive!")

    def evaluate_model(self, model: ModelWrapper,
                       benchmarks: Dict[str, Tuple[Dict[str, Any], List[Any]]],
                       check_length_bias: bool = False,
                       interaction_mode: str = "concat"):
        """
        Automatically evaluates the model on all given benchmark datasets.
        Can also evaluate the model on adversarial input to check if it is biased towards predicting longer sequences
        as positive (or vice versa).

        :param model: ModelWrapper object that supports from_embeddings(embeddings) -> list with predictions
        :param benchmarks: Dict that contains the benchmark name as key and
                           a tuple with embeddings and targets as values
        :param check_length_bias: If True, the model will be probed with randomly created sequences to check if it
                                  tends to predict longer sequences as positive or vice versa.
                                  Only works for ProtTransT5XLU50Embedder at the moment.
        :param interaction_mode: Interaction mode ("multiply" or "concat") as defined in biotrainer. Only necessary
                                 if check_length_bias is True.
        """

        print(f"**Evaluating model on {len(benchmarks.keys())} benchmarks:**")

        for benchmark, (benchmark_embeddings, benchmark_targets) in benchmarks.items():
            predictions = model.from_embeddings(embeddings=benchmark_embeddings)

            metrics_df = calculate_all_metrics(predictions=predictions, targets=benchmark_targets, df_name=benchmark)
            print(f"**Metrics for benchmark {benchmark}:**")
            print(metrics_df)

        if check_length_bias:
            self._create_random_sequence_length_benchmark(model=model, interaction_mode=interaction_mode)
