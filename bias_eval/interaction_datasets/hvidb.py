import pandas as pd

from bias_eval.taxonomy import Taxonomy
from bias_eval.dataset_base_classes import DatasetHVI
from bias_eval.dataset_base_classes import DatasetHVIStandardized


class DatasetHVIDB(DatasetHVI):
    """
    Human-Virus-Interaction DataBase (Yang et al. 2021, https://doi.org/10.1093/bib/bbaa425)

    URL: http://zzdlab.com/hvidb/
    Includes HPIDB, PHISTO, VirHostNet, VirusMentha, PDB et. al.
    Does not feature an MI-Score value for its interactions

    Fields from CSV:
        Uniprot_human:              Uniprot ID for human protein, e.g. Q6P0Q8
        Uniprot_virus:              Uniprot ID for viral protein, e.g. Q8JJY9
        EntryName_human:            Uniprot entry name for human protein, e.g. MAST2_HUMAN
        EntryName_virus:            Uniprot entry name for viral protein, e.g. Q8JJY9_9RHAB
        Organism_human:             Organism name human (Art/Species), i.e. Homo sapiens (human)
        Organism_virus:             Organism name virus (Art/Species), e.g. Rabies lyssavirus
        Organism_Interactor_human:  ID for human in Uniprot, e.g. 9606
        Organism_Interactor_virus:  ID for virus in Uniprot, e.g. 11292
        Experimental_System:        Experiment with which interaction was confirmed, e.g. phage display (might be NaN)
        Pubmed_ID:                  Pubmed ID of differential expression analysis, e.g. 1610791
        Interaction_Type:           Type of interaction, e.g. (physical) association
        Source_Database:            Database where interaction was collected from, e.g. PDB, PHISTO, HPIDB
        Complex_structure:          If complex structure was experimentally verified, i.e. Experimentally verified structure (or NaN)
        Short:                      Abbreviation of virus name, e.g. HIV-1, H3N2, Rabies lyssavirus
        HDF:                        Host dependency factors, i.e. no/yes
        HRF:                        Host restriction factors, i.e. no/yes
        Viral_family:               Family of virus species, e.g. Rhabdoviridae
        Human_GeneName:             Genes associated with human protein, e.g. MAST2 (searchable in Uniprot)
        Human_ProteinName:          Name of human protein, e.g. Importin-11
        Human_GeneID:               Gene ID associated with human protein, e.g. 23139 (searchable in Uniprot)
        Virus_GeneName:             Gene(s) associated with virus, e.g. gp, vpr, gag
        Virus_ProteinName:          Name of viral protein, e.g. Glycoprotein
        Virus_GeneID:               Gene ID associated with viral protein, e.g. 5176191 (searchable in Uniprot, might be NaN)
    """

    delimiter = "\t"
    name = "hvidb"
    header = "infer"

    def get_lyssa_associated_interactions(self):
        return DatasetHVIDB(data_frame=self.data_frame[self.data_frame["Organism_Interactor_virus"].isin(["11292"])])

    def get_sars_cov_2_associated_interactions(self):
        return DatasetHVIDB(data_frame=self.data_frame[self.data_frame["Organism_Interactor_virus"].isin(["2697049"])])

    def get_rhabdoviridae_associated_interactions(self):
        return DatasetHVIDB(data_frame=self.data_frame[self.data_frame["Viral_family"].isin(["Rhabdoviridae"])])

    def get_coronaviridae_associated_interactions(self):
        return DatasetHVIDB(data_frame=self.data_frame[self.data_frame["Viral_family"].isin(["Coronaviridae"])])

    def to_standardized_dataset(self, taxonomy: Taxonomy) -> DatasetHVIStandardized:
        series_protein_human = self.data_frame["Uniprot_human"].values
        series_protein_virus = self.data_frame["Uniprot_virus"].values

        series_taxon_virus = self.data_frame["Organism_Interactor_virus"].values

        series_dataset_name = pd.Series([self.name] * len(self.data_frame))

        series_experimental = pd.Series(["False"] * len(self.data_frame))

        series_target = pd.Series(["1"] * len(self.data_frame))

        series = {"Uniprot_human": series_protein_human,
                  "Uniprot_virus": series_protein_virus,
                  "Taxon_virus": series_taxon_virus,
                  "Dataset": series_dataset_name,
                  "Experimental": series_experimental,
                  "Target": series_target}

        standardized_df = pd.DataFrame(series)
        hvi_std = DatasetHVIStandardized(data_frame=standardized_df, name=self.name, taxonomy=taxonomy)
        return hvi_std
