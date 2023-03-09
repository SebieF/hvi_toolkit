import re
import pandas as pd

from Bio import SeqIO
from hvi_datasets.taxonomy import Taxonomy
from hvi_datasets.hvi_abstract_dataset import DatasetHVI
from hvi_datasets.hvi_standardized_dataset import DatasetHVIStandardized


class DatasetTsukiyamaNegatives(DatasetHVI):
    """
    Negative dataset from LSTM-PHV: prediction of human-virus protein–protein interactions by LSTM with word2vec
    (Tsukiyama et al. 2021, https://doi.org/10.1093/bib/bbab228)

    Negatives in this set were created by the dissimilarity negative sampling method.

    URL: http://kurata35.bio.kyutech.ac.jp/LSTM-PHV/download_page

    human_pro: Human protein ID
    virus_pro: Virus protein ID
    human_seq: human sequence
    virus_seq: viral sequence
    """

    delimiter = ","
    name = "tsukiyama_negatives"
    header = "infer"

    def __init__(self, file_path: str = None, data_frame: pd.DataFrame = None, fasta_file: str = None):
        super().__init__(file_path, data_frame)

        if file_path and (fasta_file is None):
            raise Exception(f"To construct tsukiyama_negatives dataset from file, an associated fasta file"
                            f"containing the sequences must be given!")

        # Remove everything that is not from uniprotkb
        self.data_frame = self.data_frame[self.data_frame["virus_pro"].str.contains("uniprot")]
        self.data_frame = self.data_frame.reset_index(drop=True)
        self.data_frame["virus_pro"] = self.data_frame["virus_pro"].apply(
            lambda protein_id: protein_id.split(":")[-1])  # Remove uniprotkb annotation
        if fasta_file:
            self.__assign_virus_taxa(fasta_file)

    def __assign_virus_taxa(self, fasta_file: str):
        """
        Assigns the taxon ids of the viral proteins necessary for the standardized dataset

        :param fasta_file: Path to the fasta file which contains the taxon id (OX=ID)
        """

        # Read attributes from fasta file
        seq_records = SeqIO.parse(fasta_file, "fasta")
        attribute_dict = dict()
        for sequence in seq_records:
            attribute_dict[sequence.id] = {key: value for key, value
                                           in re.findall(r"([A-Z_]+)=(-?[A-z0-9]+-?[A-z0-9]*[.0-9]*)",
                                                         sequence.description)}
        attribute_dict = {key.split("|")[1]: value for key, value in attribute_dict.items()}

        taxa_list = []
        unknown_ids = set()
        for _, row in self.data_frame.iterrows():
            virus_id = row["virus_pro"]
            if "-" in virus_id:
                virus_id = virus_id.split("-")[0]
            if virus_id not in attribute_dict.keys():
                unknown_ids.add(virus_id)
                taxa_list.append("0")  # Unknown, taxonomy starts at 1
            else:
                taxa_list.append(attribute_dict[virus_id]["OX"])

        taxon_series = pd.Series(taxa_list)
        assert len(self.data_frame["virus_pro"]) == len(taxon_series)
        self.data_frame["Taxon_virus"] = taxon_series

        print("Successfully assigned virus taxa to tsukiyama negatives dataset!")
        print(f"Unknown ids: {len(unknown_ids)}, {unknown_ids}")
        print(taxa_list[0:5])

    def to_standardized_dataset(self, taxonomy: Taxonomy) -> DatasetHVIStandardized:
        series_protein_human = self.data_frame["human_pro"]
        series_protein_virus = self.data_frame["virus_pro"]

        series_dataset_name = pd.Series([self.name] * len(series_protein_human))
        series_taxon_virus = self.data_frame["Taxon_virus"]

        series_experimental = pd.Series([False] * len(series_protein_human))

        series = {"Uniprot_human": series_protein_human,
                  "Uniprot_virus": series_protein_virus,
                  "Taxon_virus": series_taxon_virus,
                  "Dataset": series_dataset_name,
                  "Experimental": series_experimental}

        standardized_df = pd.DataFrame(series)
        hvi_std = DatasetHVIStandardized(data_frame=standardized_df, name=self.name, taxonomy=taxonomy)
        return hvi_std

    def get_lyssa_associated_interactions(self) -> DatasetHVI:
        return DatasetTsukiyamaNegatives(data_frame=self.data_frame[self.data_frame["Taxon_virus"].isin(["11292"])])

    def get_sars_cov_2_associated_interactions(self) -> DatasetHVI:
        return DatasetTsukiyamaNegatives(data_frame=self.data_frame[self.data_frame["Taxon_virus"].isin(["2697049"])])

    def get_rhabdoviridae_associated_interactions(self) -> DatasetHVI:
        """
        No viral family included in the tsukiyama dataset. To get interactions associated with a certain viral
        family, convert the dataset to a standardized dataset while providing a taxonomy object.

        :raise: ValueError
        """
        raise ValueError

    def get_coronaviridae_associated_interactions(self) -> DatasetHVI:
        """
        No viral family included in the tsukiyama dataset. To get interactions associated with a certain viral
        family, convert the dataset to a standardized dataset while providing a taxonomy object.

        :raise: ValueError
        """
        raise ValueError
