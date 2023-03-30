# Rabies Lyssavirus experimentally known interactions

`lyssa_experimental_interactions.csv` and `lyssa_experimental_interactions.ods`
(The two files are equivalent, but .ods files are a bit easier to maintain)


Interactions have been extracted from a paper of Zandi et al. 2021 
(https://doi.org/10.52547%2Fibj.25.4.226).

Note that the presented interactions are not fully exhaustive, because some taxonomy terms
are ambivalent and could include more species, and hence more interactions, than presented here.
Interactions that do not associate to Lyssavirus rabies, like European bat lyssavirus 2, have also been omitted.

**Description of columns:**
* *Protein*: Name of the Rabies Lyssavirus protein gene (G, N, M, L, P)
* *Lyssavirus*: Species of observed interaction
* *Taxon_virus*: Taxon of species
* *Uniprot_virus*: UniprotKB ID of Lyssavirus protein
* *Uniprot_human*: UniprotKB ID of human protein
* *Host protein interactor*: Name of the host (human) molecule
* *Method of PPI detection and references*: Method how ppi was confirmed and references (see paper above)
* *Viral_family*: Viral family of Lyssavirus, i.e. Rhabdoviridae for all entries

**Glossary**:

* CVS Strain: Challenge Virus Standard (CVS) Lyssavirus strain, laboratory virus often used for
experiments, taxon id: 11294
* Street strain: Rabies virus (isolate Human/Algeria/1991), 
taxon id 31613 as given by https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?lvl=0&id=31613
* PV strain: Rabies virus (strain Pasteur vaccins / PV), a strain used for vaccines, taxon id: 103929
* ERA strain: Evelyn-Rokitnicki-Abelseth (ERA) rabies strain, in a genetically modified form used for vaccines
(see also: https://doi.org/10.1016/j.antiviral.2015.06.011), taxon id: 11295
* HEP-Flury strain: Rabies virus (strain HEP-Flury), attenuated strain used for vaccines, taxon id: 11296
* India strain: Rabies virus (strain India), virus extracted in indian mammals, taxon id: 445790
* For "virulent strains", the street strain (31613) and india strain (445790) have been used
* For "attenuated strains", PV (103929) and HEP-Flury (11296) have been used
* SAD strain: Represented by Rabies virus (strain SAD B19), Street Alabama Dufferin (SAD) also used for 
vaccination, taxon id: 11300
* RABV: If no strain is specified, india strain is used (445790)
* SHBRV: Rabies virus (strain silver-haired bat-associated), associated with most human rabies cases in the US
  (https://doi.org/10.1371%2Fjournal.pone.0155542), taxon id: 445793

* DLC1,2: Dynein light chain 1/2

**Viruses Omitted**:
* Field strains, EBLV-1, EBLV-2, street virus GX/09, WCBV, MOKV, LBV, Thailand (not in swiss prot), a street strain
  (31613 M not in Uniprot), P3 and P5 isoform, ABLV, 8743THA, Ni, Ni-CE, street strains 1088 and HCM-9, DRV-Mexico
* no differentiation between different CVS variants (like CVS-11/CVS-24, CVS-B2c)
* no differentiation between different SAD variants (B19, l16)

**Unclear/Omitted human proteins**:
* HSP70 1A/1B
* RelAp43 (splicing analog of RelA, https://www.researchgate.net/figure/Residues-at-positions-77-and-104-are-necessary-for-the-interaction-of-M-Tha-with-RelAp43_fig3_311867837)
* TUB α/β
* NCL
* FAK
* Mitochondrial complex I
* HSP90AA1/Cdc37 complex separated to individual interactions
* IIGP1 (mouse protein)