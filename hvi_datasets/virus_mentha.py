from __future__ import annotations

import pandas as pd

from taxonomy import Taxonomy
from hvi_abstract_dataset import DatasetHVI
from hvi_standardized_dataset import DatasetHVIStandardized


class DatasetVirusMentha(DatasetHVI):
    """
    VirusMentha: a new resource for virus-host protein interactions
    (Calderone et al. 2014, https://doi.org/10.1093/nar/gku830)

    URL: https://virusmentha.uniroma2.it/download.php
    Mostly included by HVIDB

    Fields from CSV:
        Protein A: Uniprot ID of Protein
        Gene A: Uniprot ID of Gene
        Taxon A: Uniprot ID of Taxon (Species) (Baltimore)
        Family A: Uniprot ID of Family (ICTV)
        Protein B: -"-
        Gene B: -"-
        Taxon B: -"-
        Family B: -"-
        Score: (Virus)MINT score (interaction reliability, > 0.6 => Experimental (strict))
        PMID: Pubmed ID
    """

    delimiter = ";"
    name = "virus_mentha"
    header = "infer"

    def get_lyssa_associated_interactions(self) -> DatasetVirusMentha:
        # No lyssa interactions in dataset currently
        taxon_a = self.data_frame[self.data_frame["Taxon A"].isin([11292])]
        taxon_b = self.data_frame[self.data_frame["Taxon B"].isin([11292])]
        return DatasetVirusMentha(data_frame=pd.concat([taxon_a, taxon_b], ignore_index=True))

    def get_sars_cov_2_associated_interactions(self) -> DatasetVirusMentha:
        # No Sars-Cov-2 interactions in dataset currently
        taxon_a = self.data_frame[self.data_frame["Taxon A"].isin([2697049])]
        taxon_b = self.data_frame[self.data_frame["Taxon B"].isin([2697049])]
        return DatasetVirusMentha(data_frame=pd.concat([taxon_a, taxon_b], ignore_index=True))

    def get_rhabdoviridae_associated_interactions(self) -> DatasetVirusMentha:
        family_a = self.data_frame[self.data_frame["Family A"].isin([11270])]
        family_b = self.data_frame[self.data_frame["Family B"].isin([11270])]
        return DatasetVirusMentha(data_frame=pd.concat([family_a, family_b], ignore_index=True))

    def get_coronaviridae_associated_interactions(self) -> DatasetVirusMentha:
        family_a = self.data_frame[self.data_frame["Family A"].isin([11118])]
        family_b = self.data_frame[self.data_frame["Family B"].isin([11118])]
        return DatasetVirusMentha(data_frame=pd.concat([family_a, family_b], ignore_index=True))

    def to_standardized_dataset(self, taxonomy: Taxonomy) -> DatasetHVIStandardized:
        human_hosts_a = self.data_frame[self.data_frame["Family A"].isin([9606])]
        human_hosts_b = self.data_frame[self.data_frame["Family B"].isin([9606])]

        series_protein_human_a = human_hosts_a["Protein A"]
        series_protein_human_b = human_hosts_b["Protein B"]
        series_protein_human = pd.concat([series_protein_human_a, series_protein_human_b], ignore_index=True)

        series_protein_virus_a = human_hosts_a["Protein B"]
        series_protein_virus_b = human_hosts_b["Protein A"]
        series_protein_virus = pd.concat([series_protein_virus_a, series_protein_virus_b], ignore_index=True)

        series_taxon_virus_a = human_hosts_a["Taxon B"]
        series_taxon_virus_b = human_hosts_b["Taxon A"]
        series_taxon_virus = pd.concat([series_taxon_virus_a, series_taxon_virus_b], ignore_index=True)

        series_dataset_name = pd.Series([self.name] * len(self.data_frame))

        series_experimental_a = human_hosts_a["Score"] > 0.6
        series_experimental_b = human_hosts_b["Score"] > 0.6
        series_experimental = pd.concat([series_experimental_a, series_experimental_b], ignore_index=True)

        series = {"Uniprot_human": series_protein_human,
                  "Uniprot_virus": series_protein_virus,
                  "Taxon_virus": series_taxon_virus,
                  "Dataset": series_dataset_name,
                  "Experimental": series_experimental}

        standardized_df = pd.DataFrame(series)
        hvi_std = DatasetHVIStandardized(data_frame=standardized_df, name=self.name, taxonomy=taxonomy)
        return hvi_std
