# Human-Virus Interactions (HVI) Toolkit

Tools to work with human-virus interactions and evaluate biases associated to models and datasets.

This repository helps you to work with human-virus interactions (hvi) (or general protein-protein interactions).

You might have encountered one of the following problems yourself:
* There are a lot of different human-virus interaction datasets, but they use *different data formats*.
* The interactions provided in these datasets are often *duplicates*.
* You created a predictor that is doing well on training and validation steps, but has *poor generalization capability*.

To ease the work with different datasets, this repository provides `python` implementations for a lot of existing
hvi datasets. They can all be converted to a standardized dataset by simply calling one method. This makes them
comparable, easy to analyze and to merge ([HVI Dataset Standards](#hvi-dataset-standards)). 
In addition, the standardized format ensures that all interactions contained
are unique. They can also be used for subsequent *redundancy reduction* and *dataset splitting*.

A lot of unique interactions are not sufficient on their own: It also needs to be checked, if the datasets contain
some biases. The *hvi_toolkit* contains automatic checks to evaluate your datasets for various known biases
([Dataset Checking](#dataset-checking)).

To train an interaction model, the *interaction_mode* from [Biotrainer]((https://github.com/sacdallago/biotrainer/)) 
can be employed. Using this tool for model training ensures a standardized format for your model. This makes it 
easy to use the model checks provided by the *hvi_toolkit* ([Model checking](#model-checking)). 
High-quality benchmarks are used to assess the model's capabilities on unseen data. 
Adversarial, generated input is, furthermore, used to shed light on biases that the model
might have picked up during training (work in progress). Models not trained via *biotrainer* can also be evaluated
by implementing a simple interface.

## HVI Dataset Standards

In order to create the best dataset possible for interaction prediction, a lot of different data sources have to be
combined. 
The `bias_eval` module contains `dataset_base_classes` that can easily be extended to add
new interaction datasets. For example, if you have interactions 
stored in [PSI-MI TAB 2.5](https://psicquic.github.io/MITAB25Format.html) format,
you can create a new interaction set like this:

```python
from hvi_toolkit.dataset_base_classes import DatasetMITAB25


class DatasetNewMITAB(DatasetMITAB25):
  """
  Description of the new dataset
  """

  delimiter = "\t"
  name = "your_dataset"
  header = None
```

The most important class is the `DatasetHVIStandardized` contained in 
`dataset_base_classes/hvi_standardized_dataset.py`. All interaction datasets must provide a
`to_standardized_dataset` method in order to transfer them to a standardized format. This format
makes it easy to combine and analyze the data sources.

Already implemented interaction datasets can be found in the `interaction_datasets` package.
Raw data associated with them is provided in `raw_data`. 
The notebook `hvi_datasets7merged_dataset/merge_datasets.ipynb` shows how to merge all the provided datasets
into a single and large dataset.

The `DatasetHVIStandardized` class can also be used for **redundancy reduction** on sequences and
embeddings. Before configuring and running the 
`hvi_datasets/redundancy_reduction/redundancy_reduction.py` script, you need to have
[mmseqs2](https://github.com/soedinglab/MMseqs2) installed in a separate conda environment:
```bash
conda create --name mmseqs2
conda activate mmseqs2
conda install -c conda-forge -c bioconda mmseqs2
```

After the redundancy reduction, the `SplitsGenerator` contained in 
`hvi_datasets/dataset_splitting/splits_generators.py` can be used to create a 
[biotrainer](https://github.com/sacdallago/biotrainer/) compatible interaction file. 

## Dataset checking

Protein-protein or human-virus interaction datasets are prone to some forms of biases.
The `DatasetEvaluator` class from `bias_eval.evaluators` provides automatic tests 
to check your datasets for these biases.

It provides the following checks:

* **Check dataset bias:** Do proteins associated with interactions have the same frequency in the positive
and negative dataset?
* **Bias prediction baselines** (Using the protein frequencies in the dataset):
  * For the whole dataset
  * For provided test sets
* **Viral families:** Are they uniformly distributed across the interactions?
* **Sequence length:** Does it differ significantly between positive and negative interactions

An exemplary notebook `check_dataset.ipynb` is provided in the `examples` directory. 

## Model checking

Models trained on biased datasets will ultimately incorporate these biases. The `ModelEvaluator` class aims to
provide automatic tests to check, how well your model does. To do this, benchmark datasets have to be provided.
At the moment, experimentally verified interactions for *Rabies Lyssavirus* and the *Negatome 2.0* datasets are provided
as examples in the `examples/check_model.ipynb` notebook. 

The `model` directory contains an already trained model, including the `out.yml` file from biotrainer and 
the checkpoints for all splits is provided in `model/biotrainer_output`. It is used in the examples to show how the
model checking works.

The model was trained with the `biotrainer_config.yml` configuration file that can be run via:
```bash
# Install biotrainer
conda create --name biotrainer
conda activate biotrainer
git clone https://github.com/sacdallago/biotrainer.git
cd biotrainer
pip install .

# Run configuration file
cd lyssa_predict/model
biotrainer biotrainer_config.yml  # Biotrainer version >0.3.0 required
```

### Rabies Lyssavirus benchmark

For the *hvi_toolkit*, a high quality dataset of experimentally verified, positive interactions for the 
*Rabies Lyssavirus* was compiled.
It includes 53 experimentally verified interactions extracted from a paper of 
Zandi et al. 2021 (https://doi.org/10.52547%2Fibj.25.4.226). 
It can be found under `hvi_datasets/raw_data/rabies_lyssavirus_zandi`. The `README.md` file located in the same
directory contains additional information, how the dataset was compiled.
