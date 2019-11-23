import os
import itertools
import collections
from utilities import file_utilities
from genome_parsing import load_genome_annotations
from progress.bar import Bar

import logging
logger = logging.getLogger(__name__)

def index_all_pairs(pfam_annotations_file_list,
                    pfam_families, args):
    """This function create an int index for all domain pairs in the dataset.
    This allows to speed up the following steps and to easily map pairs from different organisms between each other.
    """
    logger.info("Indexing all pairs of Pfam domains co-occurring in at least one assembly")

    # Imports indexes for domain pairs (if they exist)
    json_indexes_filename = args["pairs_indexes_filename"]
    json_reverse_indexes_filename = args["reversed_indexes_filename"]

    # If the reversed index file is not there, which means the indexes haven't been generated
    if not os.path.exists(json_reverse_indexes_filename) or args["force"] is True:
        # Initialises a dictionary for the pairs indexes (ordered to make it easier to
        # do a reverse mapping later)
        all_domain_pair_indexes = collections.OrderedDict()

        if args['verbose'] is True:
            pbar = Bar("Indexing pairs from genome assemblies", max=len(pfam_annotations_file_list))

        for pfam_annotations_file in pfam_annotations_file_list:
            this_organism_pair_indexes = {}
            # Get the name of the organisms
            organism_name = file_utilities.get_organism_from_file(pfam_annotations_file,
                                                                  args["pfam_annotations_filename_suffix"])

            # Get the coordinates of Pfam domains on this genome
            gene_to_chr_mapping = file_utilities.load_json(
                file_utilities.get_coord_file_from_id(organism_name,
                                                      args["parsed_coordinates_filename_suffix"],
                                                      args["gene_locations_filename_suffix"],
                                                      args["out_dir"])[1])
            main_chromosome = [args["genome_sizes"][organism_name][0]]
            for genome_molecule in main_chromosome:
                genome_id = genome_molecule[0]

                # Generates all possible pairs in this dataset
                all_pairs = generate_domain_pairs_set(load_genome_annotations.import_pfam_hits(pfam_annotations_file,
                                                                                               gene_to_chr_mapping,
                                                                                               args["evalue"],
                                                                                               pfam_families,
                                                                                               genome_id).keys())

                # Indexes them (if needed)
                for pair in all_pairs:
                    index_domain_pair(pair, all_domain_pair_indexes, this_organism_pair_indexes)

            # Saves indexes for this organism
            individual_indexes_filename = os.path.join(args["pairs_indexes_dir"],
                                                       args["pairs_indexes_individual_template"].
                                                       format(organism_name))
            file_utilities.save_json(this_organism_pair_indexes, individual_indexes_filename)

            if args['verbose'] is True:
                pbar.next()

        if args['verbose'] is True:
            pbar.finish()

        # Creates a json of the final indexes dictionary
        logger.info("Saving %s indexes for Pfam domain pairs", '{:,}'.format(len(all_domain_pair_indexes)))
        save_definitive_index_file(all_domain_pair_indexes, json_indexes_filename)

        logger.info("%s indexes for Pfam domain pairs saved in file %s",
                    '{:,}'.format(len(all_domain_pair_indexes)), json_indexes_filename)

        # Inverts the key -> value mapping_tables of the indexes and saves it as a json item (a list because
        # we don't need a dictionary in this case and it would take more memory)
        save_reversed_indexes(all_domain_pair_indexes,
                              json_indexes_filename,
                              json_reverse_indexes_filename)

    else:
        all_domain_pair_indexes = file_utilities.load_json(json_reverse_indexes_filename)

    return len(all_domain_pair_indexes)


def generate_domain_pairs_set(domain_IDs):
    domain_pairs = itertools.combinations(domain_IDs, 2)
    return domain_pairs


def index_domain_pair(pair, all_domain_pair_indexes, this_organism_pair_indexes):
    pair = "_".join(sorted(pair)).replace("PF", "")
    try:
        pair_index = all_domain_pair_indexes[pair]
    except KeyError:
        pair_index = len(all_domain_pair_indexes)
        all_domain_pair_indexes[pair] = pair_index
    this_organism_pair_indexes[pair] = pair_index


def get_domain_pair_index(pair, all_domain_pair_indexes):
    pair = "_".join(sorted(pair)).replace("PF", "")
    return all_domain_pair_indexes[pair]


def save_definitive_index_file(all_domain_pair_indexes, json_indexes_filename):
    file_utilities.save_json(all_domain_pair_indexes, json_indexes_filename, compress=True)


def save_reversed_indexes(all_domain_pair_indexes,
                          json_indexes_filename,
                          json_reverse_indexes_filename):
    logger.info("Saving final reversed indexes of Pfam domain pairs")

    if len(all_domain_pair_indexes) == 0:
        all_domain_pair_indexes = file_utilities.load_json(json_indexes_filename)

    save_definitive_index_file(list(all_domain_pair_indexes.keys()),
                               json_reverse_indexes_filename)

    logger.info("Final reversed indexes of Pfam domain pairs saved in file %s",
                json_reverse_indexes_filename)


def load_indexes(args):
    logger.info("Indexes for Pfam domain pairs already present in file %s",
                args["pairs_indexes_filename"])
    logger.info("Loading indexes for Pfam domain pairs from file %s",
                args["pairs_indexes_filename"])
    all_domain_pair_indexes = file_utilities.load_json(args["pairs_indexes_filename"])
    logger.info("%s indexes for Pfam domain pairs loaded", '{:,}'.format(len(all_domain_pair_indexes)))

    return all_domain_pair_indexes
