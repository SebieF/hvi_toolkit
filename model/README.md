# Model

This directory contains the biotrainer configuration for model training and the model checkpoints.

### biotrainer_hvi_interactions.fasta

Compressed in `biotrainer_hvi_interactions.zip`

Test set: Rhabdoviridae, Coronaviridae, Hub interactions

Distribution:
```
Number train: 13055 (14.06955565853711 %)
Number val: 0 (0.0 %)
Number test: 79734 (85.93044434146289 %)
Positives in train set: 6351 (48.648027575641514 %)
Positives in validation set: 0 (0.0 %)
Positives in test set: 16923 (21.224320866882383 %)
```

### biotrainer_config.yml

Configuration file to create a prediction model for human-virus interactions using
biotrainer.
Once biotrainer is installed in the current environment, it can be used via:
```bash
biotrainer biotrainer_config.yml
```

Make sure that the embeddings are unzipped and in the correct location when running
(`../hvi_datasets/merged_dataset/embeddings_prottrans/reduced_embeddings_file.h5` by default)!


### biotrainer_output

Contains the following files:
* `out.yml`: Biotrainer output file, necessary for inference
* `logger_out.log`: Console output of biotrainer run
* FNN/custom_embeddings: Contains the *model checkpoints* for the 5 cross validation splits

Best split: `k_fold-strat-4_checkpoint.pt`

**Average split results**:
```
{'loss': 0.37390395473866234, 
'accuracy': 0.7852163910865784, 
'precision': 0.7304044961929321, 
'recall': 0.8915100693702698, 
'f1_score': 0.8023254036903381, 
'spearmans-corr-coeff': 0.5873315215110779, 
'matthews-corr-coeff': 0.58733149766922}
```
