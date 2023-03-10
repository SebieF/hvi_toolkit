# Model

### biotrainer_hvi_interactions.fasta

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

Make sure that the embeddings are unzipped and in the correct location when running!