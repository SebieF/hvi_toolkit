from tqdm import tqdm

from taxonomy import Taxonomy
from hvi_datasets.hvi_mitab25_dataset import DatasetMITAB25


class DatasetIntactPositive(DatasetMITAB25):
    """
    The MIntAct project—IntAct as a common curation platform for 11 molecular interaction databases
    (Orchard et al. 2013, https://doi.org/10.1093/nar/gkt1115)

    Dataset to include positive interactions from the intact database
    URL: https://ftp.ebi.ac.uk/pub/databases/intact/current/psimitab/
    File: intact.txt
    After pre-processing to include only human-virus interactions, intact_positive_human_virus.txt was created.

    **Was filtered to only contain entries with UniprotID, scores and human-virus interactions**

    For column descriptions, see https://psicquic.github.io/MITAB25Format.html
    """

    delimiter = "\t"
    name = "intact_positive"
    header = None

    @staticmethod
    def keep_only_human_virus(file_path: str, taxonomy: Taxonomy):
        """
        Function that can be used to filter intact.txt to only return human-virus interactions with UniProt IDs.

        :param file_path: Path to intact file
        :param taxonomy: Taxonomy object to check for viral species
        :return: Lines to keep from the given file
        """
        with open(file_path, "r") as tab_file:
            lines = tab_file.readlines()
            keep = []
            for line in tqdm(lines):
                # 1. Only uniprot:
                values = line.split("\t")
                if "uniprotkb:" in values[0] and "uniprotkb:" in values[1]:
                    # 2. Only human:
                    if "taxid:9606" in line:
                        # 3. Only virus:
                        non_human_taxon = values[9] if "taxid:9606" in values[10] else values[10]
                        if "|" in non_human_taxon:
                            non_humans = non_human_taxon.split("|")
                        else:
                            non_humans = [non_human_taxon]
                        for non_human in non_humans:
                            non_human_taxid = non_human.replace("taxid:", "")
                            if "(" in non_human_taxid:
                                non_human_taxid = non_human_taxid.split("(")[0]
                            try:
                                tax = taxonomy.get_name_from_id(non_human_taxid)
                                if "viridae" in tax:
                                    keep.append(line)
                                break
                            except KeyError:
                                continue
            return keep
