class Taxonomy:
    """
    Class to handle taxonomy data for viral protein-protein interactions.

    Reads relevant entries from HVIDB and NCBI taxonomy data to map taxonomy IDs to names and families.
    If a direct mapping is not available, it uses a matching algorithm to find the closest family.

    Attributes:
        _id_to_name (dict): Mapping of taxonomy IDs to names.
        _id_to_family (dict): Mapping of taxonomy IDs to viral families.

    Methods:
        get_name_from_id(taxonomy_id): Returns the name of the taxonomy with the given ID.
        get_family_from_id(taxonomy_id): Returns the family of the taxonomy with the given ID.

    Examples:
        >>> taxonomy = Taxonomy()
        >>> taxonomy.get_name_from_id(11293)
        'Rabies virus AV01'
        >>> taxonomy.get_family_from_id(11293)
        'Rhabdoviridae'
    """

    def __init__(self):
        self._id_to_name = {}
        self._id_to_family = {}

        # Read relevant entries from HVIDB
        hvidb_path = "raw_data/hvidb/HVIDB_PPIs.txt"
        with open(hvidb_path, "r") as hvidb_file:
            lines = hvidb_file.readlines()
            for line in lines:
                values = line.split("\t")
                self._id_to_name[values[7]] = values[16]

        # For rest use ncbi
        ncbi_taxonomy_path = "raw_data/taxonomy/names_taxonomy.tsv"
        with open(ncbi_taxonomy_path, "r") as ncbi_taxonomy_file:
            lines = ncbi_taxonomy_file.readlines()
            for line in lines:
                values = line.split("\t")
                if values[0] not in self._id_to_name.keys():
                    self._id_to_name[values[0]] = values[2]

    def get_name_from_id(self, taxonomy_id):
        """
        Returns the raw entry from names_taxonomy.tsv or HVIDB for the given taxonomy_id.


        :param taxonomy_id: NCBI taxonomy id
        :return: Raw name from names_taxonomy.tsv
        """
        return self._id_to_name[str(taxonomy_id)]

    def get_family_from_id(self, taxonomy_id) -> str:
        """
        Returns the viral family closest to the given taxonomy_id. The viral family is at first retrieved from
        HVIDB, which contains a mapping from taxonomy_id to viral family. If it is not found there, it is looked for
        in "raw_data/names_taxonomy.tsv" downloaded from NCBI. Because there does not seem to exist a direct mapping
        from taxonomy_id to viral family, this method finds the nearest entry that contains "viridae" previous
        to the given taxonomy_id. This is possible because taxons are ordered by family to species in ascending
        order.
        Once a family has been found, it is cached in self._id_to_family for later use.

        :param taxonomy_id: NCBI taxonomy id
        :return: Viral family closest to the given taxonomy id
        """

        # Use cached result if possible
        try:
            return self._id_to_family[str(taxonomy_id)]
        except KeyError:
            pass

        # Go backwards until "viridae" is in name
        taxon_idx = int(taxonomy_id)
        passed_indexes = []
        result = "Unknown"
        for tax_idx in range(taxon_idx, -1, -1):
            try:
                current_taxon = self._id_to_name[str(tax_idx)]
            except KeyError:
                continue
            passed_indexes.append(str(tax_idx))
            if "viridae" in current_taxon:
                result = current_taxon
                break
        for passed_index in passed_indexes:
            self._id_to_family[passed_index] = result
        return result
