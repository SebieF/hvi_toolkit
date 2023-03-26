from __future__ import annotations

import pandas as pd

from bias_eval.taxonomy import Taxonomy

from .hvi_abstract_dataset import DatasetHVI
from .hvi_standardized_dataset import DatasetHVIStandardized


class DatasetMITAB25(DatasetHVI):
    """
    General dataset class for datasets stored in the MITAB25 format.
    VirHostNet3, SarsCov2IM28880 are directly implemented using this class.

    **Gets filtered to only contain entries with UniprotID and score values!**

    For column descriptions, see https://psicquic.github.io/MITAB25Format.html
    """

    delimiter = "\t"
    name = "mi_tab_25"
    header = None

    def __init__(self, file_path: str = None, data_frame: pd.DataFrame = None):
        super().__init__(file_path, data_frame)

        if file_path is not None:  # Newly constructed dataset
            def omit_annotation(value: str):
                # Remove any database or score annotations given for the protein identifiers or scores
                value = str(value)
                result = ""
                if "uniprotkb:" in value:
                    dbs = value.split("|")
                    for db in dbs:
                        if "uniprotkb:" in db:
                            result = db.split("uniprotkb:")[1]
                            break
                elif "psi-mi:" in value:
                    result = value.replace("psi-mi:", "")
                elif "intact-miscore:" in value:
                    result = float(value.split("intact-miscore:")[1])
                elif "|" in value:
                    multiple_values = value.split("|")
                    for val in multiple_values:
                        if ":" in val:
                            result += val.split(":")[1]
                        else:
                            result += val
                elif ":" in value:
                    result = value.split(":")[1]
                else:
                    result = value
                return result

            # Pre-filter only Uniprot:
            self.data_frame = self.data_frame[self.data_frame[0].str.contains("uniprotkb:")]
            self.data_frame = self.data_frame[self.data_frame[1].str.contains("uniprotkb:")]

            # Pre-filter without confidence value:
            self.data_frame = self.data_frame[~(self.data_frame[14] == "-")]

            # Drop annotations
            for key in self.data_frame.keys():
                self.data_frame[key] = self.data_frame[key].apply(omit_annotation)

            self.data_frame = self.data_frame.rename(columns={0: "Uniprot_A",
                                                              1: "Uniprot_B",
                                                              2: "AltID_A",
                                                              3: "AltID_B",
                                                              4: "Aliases_A",
                                                              5: "Aliases_B",
                                                              6: "PPI_Detection_Methods",
                                                              7: "First_Author",
                                                              8: "Publication_ID",
                                                              9: "Taxon_A",
                                                              10: "Taxon_B",
                                                              11: "Interaction_Types",
                                                              12: "Source_Databases",
                                                              13: "Interaction_IDs",
                                                              14: "Confidence_Score"})

    def get_lyssa_associated_interactions(self):
        taxon_a = self.data_frame[self.data_frame["Taxon_A"].isin(["11292"])]
        taxon_b = self.data_frame[self.data_frame["Taxon_B"].isin(["11292"])]
        return DatasetMITAB25(data_frame=pd.concat([taxon_a, taxon_b], ignore_index=True))

    def get_sars_cov_2_associated_interactions(self):
        taxon_a = self.data_frame[self.data_frame["Taxon_A"].isin(["2697049"])]
        taxon_b = self.data_frame[self.data_frame["Taxon_B"].isin(["2697049"])]
        return DatasetMITAB25(data_frame=pd.concat([taxon_a, taxon_b], ignore_index=True))

    def get_rhabdoviridae_associated_interactions(self):
        """
        No viral family included in the MITAB25-format. To get interactions associated with a certain viral
        family, convert the dataset to a standardized dataset while providing a taxonomy object.

        Example:
        >>>taxonomy = Taxonomy()
        >>>mitab = DatasetMITAB25(file_path=path_to_mitab_file)
        >>>mitab_std = mitab.to_standardized_dataset(taxonomy=taxonomy)
        >>>mitab_std_rhab = mitab_std.get_rhabdoviridae_associated_interactions()

        :raise: ValueError
        """
        raise ValueError

    def get_coronaviridae_associated_interactions(self):
        """
        No viral family included in the MITAB25-format. To get interactions associated with a certain viral
        family, convert the dataset to a standardized dataset while providing a taxonomy object.

        Example:
        >>>taxonomy = Taxonomy()
        >>>mitab = DatasetMITAB25(file_path=path_to_mitab_file)
        >>>mitab_std = mitab.to_standardized_dataset(taxonomy=taxonomy)
        >>>mitab_std_coro = mitab_std.get_coronaviridae_associated_interactions()

        :raise: ValueError
        """
        raise ValueError

    def to_standardized_dataset(self, taxonomy: Taxonomy):
        human_ids = ["9606", "9606(human)9606(Homo sapiens)", "9606(human)"]
        human_hosts_a = self.data_frame[self.data_frame["Taxon_A"].isin(human_ids)]
        human_hosts_b = self.data_frame[self.data_frame["Taxon_B"].isin(human_ids)]

        series_protein_human_a = human_hosts_a["Uniprot_A"]
        series_protein_human_b = human_hosts_b["Uniprot_B"]
        series_protein_human = pd.concat([series_protein_human_a, series_protein_human_b], ignore_index=True)

        series_protein_virus_a = human_hosts_a["Uniprot_B"]
        series_protein_virus_b = human_hosts_b["Uniprot_A"]
        series_protein_virus = pd.concat([series_protein_virus_a, series_protein_virus_b], ignore_index=True)

        def omit_annotation_for_std(value: str):
            result = value
            if "(" in value:
                result = value.split("(")[0]
            return result

        series_taxon_virus_a = human_hosts_a["Taxon_B"].apply(omit_annotation_for_std)
        series_taxon_virus_b = human_hosts_b["Taxon_A"].apply(omit_annotation_for_std)
        series_taxon_virus = pd.concat([series_taxon_virus_a, series_taxon_virus_b], ignore_index=True)

        series_dataset_name = pd.Series([self.name] * len(series_taxon_virus))

        series_experimental_a = human_hosts_a["Confidence_Score"].astype(float) > 0.6
        series_experimental_b = human_hosts_b["Confidence_Score"].astype(float) > 0.6
        series_experimental = pd.concat([series_experimental_a, series_experimental_b], ignore_index=True)

        series_target = pd.Series(["1"] * len(series_taxon_virus))

        series = {"Uniprot_human": series_protein_human,
                  "Uniprot_virus": series_protein_virus,
                  "Taxon_virus": series_taxon_virus,
                  "Dataset": series_dataset_name,
                  "Experimental": series_experimental,
                  "Target": series_target}

        standardized_df = pd.DataFrame(series)
        hvi_std = DatasetHVIStandardized(data_frame=standardized_df, name=self.name, taxonomy=taxonomy)
        return hvi_std
