# Examples

In this directory, two examples are given how the `DatasetEvaluator` and `ModelEvaluator` classes can be applied.

## Dataset checking

Protein-protein or human-virus interaction datasets are prone to some forms of biases.
The `DatasetEvaluator` class provides automatic tests to check your datasets for these biases.
It provides the following checks:

* **Check dataset bias:** Do proteins associated with interactions have the same frequency in the positive
and negative dataset?
* **Based on the dataset bias** - Metrics are given for predicting the interactions from frequencies:
  * For the whole dataset
  * For provided test sets
* **Viral families:** Are they uniformly distributed across the interactions?
* **Sequence length:** Does it differ significantly between positive and negative interactions

## Model checking

Models trained on biased datasets will ultimately incorporate these biases. The `ModelEvaluator` class aims to
provide automatic tests to check, how well your model does. To do this, benchmark datasets have to be provided.
At the moment, experimentally verified interactions for *Rabies Lyssavirus* and the Negatome 2.0 datasets are provided
as examples.

If you are not using a model trained via [biotrainer](https://github.com/sacdallago/biotrainer), you can still use
`ModelEvaluator`. All you need to do is to provide a `ModelWrapper` class that creates interaction predictions 
from embeddings (i.e. from a vector of floats):
```python

from bias_eval.evaluators import ModelWrapper
from typing import Union, Iterable, Dict, List, Any

class ModelWrapperCustom(ModelWrapper):

    def from_embeddings(self, embeddings: Union[Iterable, Dict]) -> List[Any]:
        """
        Wrapper method to calculate predictions from embeddings. If you want to evaluate your model via the
        ModelEvaluator, but you are not using biotrainer and its inferencer module, you can inherit from this
        class and overwrite this method. The method must return a list with the predictions.

        :param embeddings: Iterable or dictionary containing the input embeddings to predict on.
        :return: List with class or value prediction from the given embeddings.
        """
        raise NotImplementedError
```

In addition, you can check your model via adversarial input. At the moment, only a prototype version to test if the
model tends to predict longer interactions to be positive (or vice versa). It relies on the `ProtTransT5XLU50Embedder`,
so your model needs to be trained with embeddings from this protein language model.
