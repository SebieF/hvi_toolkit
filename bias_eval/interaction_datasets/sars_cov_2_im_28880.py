from bias_eval.dataset_base_classes import DatasetMITAB25


class DatasetSarsCov2IM28880(DatasetMITAB25):
    """
    A proteome-scale map of the SARS-CoV-2–human contactome
    (Kim et al. 2022, https://doi.org/10.1038/s41587-022-01475-z)

    URL: http://www.imexconsortium.org/
    IM-28880

    **Gets filtered to only contain entries with UniprotID and score values!**

    For column descriptions, see https://psicquic.github.io/MITAB25Format.html
    """

    delimiter = "\t"
    name = "sars_cov_2_im28880"
    header = None

    def get_sars_cov_2_associated_interactions(self):
        return self
