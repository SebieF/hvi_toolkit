from bias_eval.dataset_base_classes import DatasetMITAB25


class DatasetVirHostNet3(DatasetMITAB25):
    """
    VirHostNet3.0, successor of VirHostNet 2.0: surfing on the web of virus/host molecular interactions data
    (Guirimand et al. 2014, https://doi.org/10.1093/nar/gku1121)

    URL: https://virhostnet.prabi.fr/
    Mostly included by HVIDB

    **Gets filtered to only contain entries with UniprotID and score values!**

    For column descriptions, see https://psicquic.github.io/MITAB25Format.html
    """
    delimiter = "\t"
    name = "vir_host_net3"
    header = None
