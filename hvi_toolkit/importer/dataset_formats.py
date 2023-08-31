from typing import List, Dict, Final, Any

from ..dataset_base_classes import DatasetMITAB25, DatasetHVIStandardized, DatasetHVI
from ..interaction_datasets import DatasetHVIDB, DatasetVirusMentha, DatasetVirusesString, DatasetTsukiyamaNegatives, DatasetRabiesLyssavirusExperimental

__FORMAT_DICT: Final[Dict[str, Any]] = {
    "HVIDB": DatasetHVIDB,
    "MITAB2.5": DatasetMITAB25,
    "VirusMentha": DatasetVirusMentha,
    "Viruses.String": DatasetVirusesString,
    "Tsukiyama2021": DatasetTsukiyamaNegatives,
    "RabiesLyssavirusExperimental": DatasetRabiesLyssavirusExperimental,
    "HviToolkitStandard": DatasetHVIStandardized
}


def get_supported_dataset_formats() -> List[str]:
    return list(__FORMAT_DICT.keys())


def get_supported_dataset_formats_with_docs() -> Dict[str, str]:
    result = {}
    for key in __FORMAT_DICT.keys():
        result[key] = __FORMAT_DICT[key].__doc__
    return result


def auto_detect_format(header: str) -> str:
    for format_str in __FORMAT_DICT.keys():
        dataset = __FORMAT_DICT[format_str]
        dataset_header = dataset.get_header()
        delimiter = dataset.delimiter
        if header == dataset_header:
            return format_str
        if delimiter in header:
            dataset_header_keys = dataset_header.split(delimiter)
            header_keys = header.split(delimiter)
            # All keys from dataset_header_keys must be contained in header_keys
            all_values_contained = all(
                [dataset_header_key in header_keys for dataset_header_key in dataset_header_keys])
            if all_values_contained:
                return format_str

    # Special case MiTab2.5 without header
    mitab25_value_samples = ["uniprotkb", "pubmed", "taxid"]
    if DatasetMITAB25.delimiter in header:
        all_values_contained = all([mitab25_value in header for mitab25_value in mitab25_value_samples])
        if all_values_contained:
            return "MITAB2.5"

    raise ValueError(f"Header {header} does not seem to match any available formats!")


def get_dataset_class_by_format(format_str: str) -> DatasetHVI.__class__:
    return __FORMAT_DICT[format_str]
