# lyssa_predict

Datasets, models and benchmarks to predict rabies lyssavirus human-virus interactions.

## Datasets

The `hvi_datasets` directory contains `dataset_base_classes` that can easily be extended to add
new interaction datasets. For example, if you have interactions 
stored in [PSI-MI TAB 2.5](https://psicquic.github.io/MITAB25Format.html) format,
you can create a new interaction set like this:

```python
from hvi_datasets.dataset_base_classes import DatasetMITAB25

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
The notebook `merged_dataset/merge_datasets.ipynb` shows how to merge all the provided datasets
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

## Models

The `model` directory contains a `biotrainer_config.yml` file that can be run via:
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

An already trained model, including the `out.yml` file from biotrainer and the checkpoints
for all splits is provided in `model/biotrainer_output`.

## Benchmarks

The `rabies_lyssavirus_predictions` contains the `lyssa_experimental_interactions.csv` file,
which includes 53 experimentally verified interactions extracted from a paper from 
Zandi et al. 2021 (https://doi.org/10.52547%2Fibj.25.4.226). They are used to benchmark
the rabies lyssavirus-human interactome predictions generated in `rabies_lyssa_predict.ipynb`
and stored in `rabies_lyssa_human_interactome_predictions.csv`.
