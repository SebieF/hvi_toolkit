from .dataset_formats import get_supported_dataset_formats, get_dataset_class_by_format

from ..taxonomy import Taxonomy
from ..dataset_base_classes import DatasetHVI, DatasetHVIStandardized, DatasetPPIStandardized


def import_dataset_by_format(dataset_path: str, format_str: str, taxonomy: Taxonomy):
    if format_str not in get_supported_dataset_formats():
        raise ValueError(f"Dataset format {format_str} not known!")

    dataset_class: DatasetHVI.__class__ = get_dataset_class_by_format(format_str)
    dataset_instance: DatasetHVI = dataset_class(dataset_path)
    dataset_instance_standardized: DatasetHVIStandardized = dataset_instance.to_standardized_dataset(taxonomy=taxonomy)
    converted_to_ppi: DatasetPPIStandardized = DatasetPPIStandardized.from_hvi_standardized(dataset_instance_standardized)

    return converted_to_ppi
