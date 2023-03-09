# Merged dataset

`merged_dataset.zip`

Contains positive interactions from:
- HVIDB
- HPIDB
- IntAct
- Sars-Cov-2 IM-28880
- VirHostNet3.0
- VirusMentha
- Viruses.STRING

Contains negative interactions from:
- Tsukiyama et al. 2021

See `"../raw_data/README.md"` for more information on the data sources.

The zip file contains only the positive data combined (`merged_dataset_positives.csv`),
positive and negative interactions combined (`merged_dataset_complete.csv`) and a fasta file
for all proteins in the complete dataset (`merged_dataset_complete.fasta`).
All these files can be created using the `merge_datasets.ipynb` notebook.

## Embeddings

* Prottrans: [Prottrans](https://github.com/agemagician/ProtTrans) embeddings for all sequences in 
the `merged_dataset_complete.fasta` have been  calculated 
using [bio_embeddings](https://github.com/sacdallago/bio_embeddings).
* Ankh: [Ankh](https://github.com/agemagician/Ankh) embeddings have been calculated via the `ankh_embed.py` script.