from __future__ import annotations

import json
import requests
import pandas as pd
import seaborn as sns
import more_itertools

from tqdm import tqdm
from taxonomy import Taxonomy
from abc import abstractmethod, ABC
from typing import Optional, List, Dict, Final, Union


class DatasetHVI(ABC):
    """
    An abstract base class that defines the basic functionality that every human-virus dataset must provide,
    once it has been converted to a pandas DataFrame from a raw data file.

    :var delimiter: Delimiter character for the input data file.
    :var name: Name of the dataset.
    :var header: Specifies if there exists a header or how it should be handled.
    :var INTERACTION_INDICATOR: Constant indicator of an interaction (human_id&virus_id)

    :param file_path: (optional) Path to the input data file.
    :param data_frame: (optional) DataFrame object to be used as input data.

    :raises Exception: Raised if no data source is given during initialization.
    """

    delimiter: str = "\t"
    name: str = "dataset_HVI"
    header: Union[str, int] = "Infer"

    INTERACTION_INDICATOR: Final[str] = "&"

    def __init__(self, file_path: str = None, data_frame: pd.DataFrame = None):
        """
        Initializes a new DatasetHVI object either from file or from an existing data_frame.
        """
        if file_path:
            self.file_path = file_path
            self.data_frame: pd.DataFrame = pd.read_csv(file_path, sep=self.delimiter, header=self.header)
        elif data_frame is not None:
            self.data_frame: pd.DataFrame = data_frame
        else:
            raise Exception("No data source given!")

        # Init plotting
        sns.set_theme()
        sns.color_palette("colorblind")

    def print_unique_values(self, max_entries: Optional[int] = -1):
        """
        Print the unique values for each column in the dataset.

        :param max_entries: The maximum number of unique values to print for each column.
                            If -1, print all unique values.
        """
        for key in self.data_frame.keys():
            print(f"Unique values for {key}: {len(self.data_frame[key].unique())}")
            if max_entries > 0:
                print(self.data_frame[key].unique()[:max_entries])
            else:
                print(self.data_frame[key].unique())
            print(f"--------\n")

    def plot_nominal_column(self, category: str, normalize: bool = False, max_values: int = 10):
        """
        Plot a bar chart of the specified nominal column.

        :param category: The name of the column to plot.
        :param normalize: Whether to normalize the counts of each category.
        :param max_values: The maximum number of categories to plot.
        """

        # Split comma separated values:
        name_list = []
        for names in self.data_frame[category].values:
            for name in names.split(","):
                name_list.append(name)

        counts = pd.Series(name_list).value_counts(normalize=normalize)
        number_categories = len(counts)
        counts = counts[0:max_values]

        title = f"{self.name}: {category} - {min(number_categories, max_values)}/{number_categories} Categories"
        counts.plot(kind="bar", title=title, color=sns.color_palette())

    @staticmethod
    def __get_sequences_with_ptm_information(uniprot_list: List[str]) -> Dict[str, str]:
        """
        Retrieve protein sequences and associated post-translational modification (PTM) information from UniProtKB
        for a list of UniProtKB accession IDs.

        This method retrieves sequences and PTM information for a list of UniProtKB accession IDs from the UniProtKB
        API using batch processing. Each batch contains up to 20 UniProtKB accession IDs. The method first separates
        UniProtKB accession IDs with PTM feature IDs from those without, and then queries the UniProtKB API with each
        batch to retrieve sequence and PTM information. If a UniProtKB accession ID has one or more PTM feature IDs,
        only the subsequence(s) corresponding to the feature(s) will be included in the result dictionary, along with
        the corresponding PTM feature ID(s) appended to the UniProtKB accession ID.

        :param uniprot_list: A list of UniProtKB accession IDs with or without PTM feature IDs.
        :return: A dictionary mapping UniProtKB accession IDs (with or without PTM feature IDs) to protein sequences.

        :raises requests.exceptions.HTTPError: If an HTTP error occurs when making requests to UniProtKB API.
        :raises requests.exceptions.RequestException: If the request made by this function was invalid for a batch.
        :raises AssertionError: If one or more PTM feature IDs specified in uniprot_list are not found for the
                                corresponding UniProtKB accession ID.
        """

        result_dict = {}  # Accession => Sequence
        batch_size = 20

        # Split the input list into smaller batches of size 'batch_size'
        for batch in tqdm(list(more_itertools.chunked(uniprot_list, n=batch_size))):
            batch_ids = []
            uniprot_with_features = {}

            # Loop through each item in the current batch
            for uniprot_kb in batch:
                # Check if the item contains PTM information
                if "-PRO" in uniprot_kb:
                    # Split the item into Uniprot ID and Feature ID
                    uniprot_id = uniprot_kb.split("-")[0]
                    feature_id = uniprot_kb.split("-")[1]

                    # Create a dictionary to store the Feature ID(s) for each Uniprot ID
                    if uniprot_id not in uniprot_with_features.keys():
                        uniprot_with_features[uniprot_id] = []

                    # Add the current Feature ID to the dictionary
                    uniprot_with_features[uniprot_id].append(feature_id)
                else:
                    # If the item does not contain PTM information, then it is just a Uniprot ID
                    uniprot_id = uniprot_kb

                # Add the current Uniprot ID to the batch
                batch_ids.append(uniprot_id)

            # Construct a query for the batch of UniProt IDs
            query = ",".join(batch_ids)
            # Make a request to the UniProt API to get the sequences
            url = f"https://www.ebi.ac.uk/proteins/api/proteins?offset=0&size={batch_size}&accession={query}"
            r = requests.get(url, headers={"Accept": "application/json"})
            # Check if the request was successful
            if not r.ok:
                r.raise_for_status()
                raise requests.exceptions.RequestException(f"Failed for batch {batch_ids}!")
            else:
                # Extract the response body and parse it as JSON
                response_body = r.text
                response_json = json.loads(response_body)
                # Iterate over each UniProt entry in the response
                for response in response_json:
                    accession = response["accession"]
                    sequence = response["sequence"]["sequence"]
                    # Check if the UniProt entry has a feature ID (PRO)
                    if accession in uniprot_with_features:
                        found = 0
                        # Iterate over each feature ID associated with the UniProt entry
                        for feature_id in uniprot_with_features[accession]:
                            # Iterate over each feature in the UniProt entry
                            for feature in response["features"]:
                                # Check if the feature ID matches the desired feature ID
                                if "ftId" in feature.keys():
                                    if feature["ftId"] == feature_id:
                                        # Extract the sequence for the feature and add it to the result dictionary
                                        begin = int(feature["begin"]) - 1
                                        end = int(feature["end"])
                                        sequence_mod = sequence[begin:end]
                                        accession_with_feature_id = accession + "-" + feature_id
                                        result_dict[accession_with_feature_id] = sequence_mod
                                        found += 1
                                        break
                        # Check if all feature IDs were found for the UniProt entry
                        assert found == len(uniprot_with_features[accession]), \
                            f"Feature ID(s) not found for {accession} - {uniprot_with_features[accession]}"
                    else:
                        # Add the UniProt entry's sequence to the result dictionary
                        result_dict[accession] = sequence

        return result_dict

    def get_interactions_by_condition(self, category: str, condition: str) -> DatasetHVI:
        """
        Method to get all interactions that satisfy a certain condition regarding the given category

        Example: get_interactions_by_condition(category="Family_virus", condition="Coronaviridae")

        :param category: Category to apply the condition in the data_frame
        :param condition: Condition that must be met as a string

        :return: Dataset of same type with only interactions that meet the provided condition
        :rtype: DatasetHVI
        """
        return self.__class__(data_frame=self.data_frame[self.data_frame[category].isin([condition])])

    @abstractmethod
    def get_lyssa_associated_interactions(self) -> DatasetHVI:
        """
        Abstract method to get interactions associated with lyssa viruses.

        :raises: NotImplementedError
        :rtype: DatasetHVI
        """
        raise NotImplementedError

    @abstractmethod
    def get_sars_cov_2_associated_interactions(self) -> DatasetHVI:
        """
        Abstract method to get interactions associated with sars-cov-2 viruses.

        :raises: NotImplementedError
        :rtype: DatasetHVI
        """
        raise NotImplementedError

    @abstractmethod
    def get_rhabdoviridae_associated_interactions(self) -> DatasetHVI:
        """
        Abstract method to get interactions associated with rhabdoviridae viruses.

        :raises: NotImplementedError
        :rtype: DatasetHVI
        """
        raise NotImplementedError

    @abstractmethod
    def get_coronaviridae_associated_interactions(self) -> DatasetHVI:
        """
        Abstract method to get interactions associated with coronaviridae viruses.

        :raises: NotImplementedError
        :rtype: DatasetHVI
        """
        raise NotImplementedError

    @abstractmethod
    def to_standardized_dataset(self, taxonomy: Taxonomy):
        """
        Abstract method to convert the dataset to a standardized format.

        :param taxonomy: Instance of the Taxonomy class for mapping virus IDs to names and families.
        :raises: NotImplementedError
        :rtype: DatasetHVIStandardized
        """
        raise NotImplementedError

    def __len__(self):
        """
        Necessary for len(dataset) to work.

        :return: Returns number of interactions in data_frame
        """
        return self.data_frame.__len__()
