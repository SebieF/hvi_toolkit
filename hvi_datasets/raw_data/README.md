# Raw data sources

## Positive Datasets

All data sources (besides *viruses_string*, see below) are included as a `.zip` file in the respective directories.

### hpidb

`hpidb_virus_25_10_22.txt`

Human-Pathogen-Interaction DataBase 3.0, successor of HPIDB 2.0: a curated database for host–pathogen interactions
(Ammari et al. 2016, https://doi.org/10.1093/database/baw103)

- URL: https://hpidb.igbb.msstate.edu/keyword.html
- Query type: Taxon Name / Species
- Query: Virus
- Do not include Interlog Predictions

### hvidb

`HVIDB_PPIs.txt`

Human-Virus-Interaction DataBase (Yang et al. 2021, https://doi.org/10.1093/bib/bbaa425)

URL: http://zzdlab.com/hvidb/

### intact

`intact_positive_human_virus.txt`

The MIntAct project—IntAct as a common curation platform for 11 molecular interaction databases
(Orchard et al. 2013, https://doi.org/10.1093/nar/gkt1115)

File only includes positive interactions from the intact database
- URL: https://ftp.ebi.ac.uk/pub/databases/intact/current/psimitab/
- File: intact.txt
- After pre-processing to include only human-virus interactions, intact_positive_human_virus.txt was created.

### sars_cov_2_im28880

`IM_28880_sars_2510.txt`

A proteome-scale map of the SARS-CoV-2–human contactome (Kim et al. 2022, https://doi.org/10.1038/s41587-022-01475-z)

- URL: http://www.imexconsortium.org/
- Query: IM-28880

### vir_host_net_3

`VirHostNet3.tsv`

VirHostNet3.0, successor of VirHostNet 2.0: surfing on the web of virus/host molecular interactions data
(Guirimand et al. 2014, https://doi.org/10.1093/nar/gku1121)

URL: https://virhostnet.prabi.fr/

### virus_mentha

`virus_mentha_human_virus_9606.csv`

VirusMentha: a new resource for virus-host protein interactions
(Calderone et al. 2014, https://doi.org/10.1093/nar/gku830)

URL: https://virusmentha.uniroma2.it/download.php

### viruses_string

**Dataset:** `9606.protein.links.detailed.v10.5.txt`

Viruses.STRING: A Virus-Host Protein-Protein Interaction Database (Cook et al. 2018, https://doi.org/10.3390/v10100519)

URL: http://viruses.string-db.org/download/protein.links.detailed.v10.5/9606.protein.links.detailed.v10.5.txt.gz
(Not included as a zip file in this repository because the file size is too large)

**ID Mapping:** `viruses_string_id_mapping.zip`

Contains the non-uniprot-id identifiers from the Viruses.STRING dataset and the respective mappings (if available)

## Negative datasets

### Tsukiyama

`tsukiyama_negatives.zip` 

Negative dataset from LSTM-PHV: prediction of human-virus protein–protein interactions by LSTM with word2vec
(Tsukiyama et al. 2021, https://doi.org/10.1093/bib/bbab228)

Negatives in this set were created by dissimilarity negative sampling method.

- `negative_samples.csv`: Actual negative interaction data
- `tsukiyama_negatives_with_taxon.fasta`: `fasta` file retrieved from uniprot that contains the taxon data for the
viral proteins

## Other

### Taxonomy

`names_taxonomy.tsv` 

Required by taxonomy.py

- Download `taxdmp.zip` from: https://ftp.ncbi.nih.gov/pub/taxonomy/
- Extract `names.dmp` and rename to `names_taxonomy.tsv` 
