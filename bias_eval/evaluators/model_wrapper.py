from abc import abstractmethod, ABC
from typing import Union, Iterable, Dict, List, Any

from biotrainer.inference import Inferencer


class ModelWrapper(ABC):

    @abstractmethod
    def from_embeddings(self, embeddings: Union[Iterable, Dict]) -> List[Any]:
        """
        Wrapper method to calculate predictions from embeddings. If you want to evaluate your model via the
        ModelEvaluator, but you are not using biotrainer and its inferencer module, you can inherit from this
        class and overwrite this method. The method must return a list with the predictions.

        :param embeddings: Iterable or dictionary containing the input embeddings to predict on.
        :return: List with class or value prediction from the given embeddings.
        """
        raise NotImplementedError


class ModelWrapperBiotrainer(ModelWrapper):
    """
    ModelWrapper class for the Inferencer from biotrainer
    """

    def __init__(self, inferencer: Inferencer, split_name: str = "hold_out"):
        super().__init__()
        self.inferencer = inferencer
        self.split_name = split_name

    def from_embeddings(self, embeddings: Union[Iterable, Dict]) -> List[Any]:
        return list(map(int, self.inferencer.from_embeddings(embeddings=embeddings,
                                                             split_name=self.split_name,
                                                             include_probabilities=False)[
            "mapped_predictions"].values()))
