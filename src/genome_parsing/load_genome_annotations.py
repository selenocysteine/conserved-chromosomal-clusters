import collections
from utilities import file_utilities


def check_pfam_family(new_domain_ID, new_domain_evalue, gene_ID,
                      gene_hits, pfam_families, domain_hits,
                      domain_centric=True):
    '''This function removes annotations for multiple Pfam domains from the same clan.'''

    domain_can_be_added = True

    try:
        other_domains_same_gene = list(gene_hits[gene_ID].keys())
    except KeyError:
        return domain_can_be_added

    try:
        new_domain_fam = pfam_families[new_domain_ID]
    except KeyError:
        return domain_can_be_added

    for old_domain_ID in other_domains_same_gene:
        old_domain_evalue = gene_hits[gene_ID][old_domain_ID]

        try:
            old_domain_fam = pfam_families[old_domain_ID]

            if new_domain_fam == old_domain_fam:
                if old_domain_evalue >= new_domain_evalue:
                    if domain_centric is True:
                        if len(domain_hits[old_domain_ID]) > 1:
                            domain_hits[old_domain_ID].remove(gene_ID)
                        else:
                            del domain_hits[old_domain_ID]
                    del gene_hits[gene_ID][old_domain_ID]

                else:
                    domain_can_be_added = False

        except KeyError:
            pass

    return domain_can_be_added


def import_pfam_hits(pfam_file_name, gene_to_chr_mapping, evalue_threshold,
                     pfam_families, genome_id, domain_centric=True):

    gene_hits = collections.OrderedDict()
    if domain_centric is True:
        domain_hits = {}
    else:
        domain_hits = None

    with open(pfam_file_name, "r") as pfam_file:
        for line in pfam_file:
            if not line.startswith("#"):
                line = line.strip().split()
                pfam_ID = line[1].split(".")[0].replace("PF", "")
                gene_ID = line[2]
                evalue = float(line[7])

                if evalue < evalue_threshold:
                    chr = gene_to_chr_mapping[gene_ID]

                    if chr == genome_id:
                        domain_can_be_added = check_pfam_family(pfam_ID, evalue, gene_ID,
                                                                gene_hits, pfam_families,
                                                                domain_hits, domain_centric)

                        if domain_can_be_added is True:
                            if domain_centric is True:
                                try:
                                    domain_hits[pfam_ID].append(gene_ID)
                                except KeyError:
                                    domain_hits[pfam_ID] = [gene_ID]

                            try:
                                gene_hits[gene_ID][pfam_ID] = evalue
                            except KeyError:
                                gene_hits[gene_ID] = {pfam_ID: evalue}

    if domain_centric is True:
        return domain_hits

    else:
        return gene_hits


def import_parsed_coordinates(coord_file):
    """Import the ordered dictionary with the parsed genomic coordinates of each gene."""
    return file_utilities.load_json(coord_file)


def import_gene_coord_data(organism_name, args):
    coordinates_file, gene_to_chr_mapping_file = \
        file_utilities.get_coord_file_from_id(organism_name,
                                              args["parsed_coordinates_filename_suffix"],
                                              args["gene_locations_filename_suffix"],
                                              args["out_dir"])
    gene_coordinates = import_parsed_coordinates(coordinates_file)
    gene_to_chr_mapping = file_utilities.load_json(gene_to_chr_mapping_file)

    return gene_coordinates, gene_to_chr_mapping