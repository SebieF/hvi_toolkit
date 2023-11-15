import taxoniq


class Taxonomy:
    """
    Class to handle taxonomy data using taxoniq

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

    def get_name_from_id(self, taxonomy_id):
        """
        Returns the name of the species for a given taxonomy id.


        :param taxonomy_id: NCBI taxonomy id
        :return: Raw name from names_taxonomy.tsv
        """
        taxon = taxoniq.Taxon(taxonomy_id)
        return taxon.scientific_name

    def get_family_from_id(self, taxonomy_id) -> str:
        """
        Returns the family name for the given taxonomy id

        :param taxonomy_id: NCBI taxonomy id
        :return: Viral family closest to the given taxonomy id
        """

        # Use cached result if possible
        taxon = taxoniq.Taxon(taxonomy_id)
        return taxon.parent.parent.parent.scientific_name
