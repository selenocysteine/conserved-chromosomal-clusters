"""Various utility functions to save / find files and directories."""
import os
import ujson
import gzip
import re

import logging

logger = logging.getLogger(__name__)


def get_file_list(data_dir, pattern, verbose=True, error=True):
    if not os.path.isdir(data_dir):
        logging.error("Directory %s does not exist", data_dir)
        raise ValueError("Directory {} does not exist".format(data_dir))

    file_list = []
    if not isinstance(pattern, list):
        pattern = [pattern]

    for dir_path, dir_names, file_names in os.walk(data_dir):
        for file_name in file_names:
            match = True
            for item in pattern:
                if item not in os.path.basename(file_name):
                    match = False
            if match is True:
                file_list.append(os.path.join(dir_path, file_name))

    if len(file_list) == 0 and error is True:
        logging.warning("No %s files retrieved in directory %s",
                        ", ".join(pattern), data_dir)
        raise ValueError("No {} files retrieved in directory {}"
                         .format(", ".join(pattern), data_dir))

    if verbose:
        logging.info("%s \'%s\' files retrieved in directory %s",
                     len(file_list),
                     ", ".join(pattern), data_dir)

    return sorted(file_list)


def save_json(object, object_json_filename, compress=False):
    if compress is True or object_json_filename.endswith(".gz"):
        if not object_json_filename.endswith(".gz"):
            object_json_filename = object_json_filename + ".gz"

        with gzip.open(object_json_filename, 'wt', encoding="ascii") as handle:
            ujson.dump(object, handle)

    else:
        with open(object_json_filename, 'w') as handle:
            ujson.dump(object, handle)

    return 0


def load_json(json_filename):
    if json_filename.endswith("gz"):
        with gzip.open(json_filename, 'rt', encoding="ascii") as handle:
            return ujson.load(handle)
    else:
        with open(json_filename, 'r') as handle:
            return ujson.load(handle)


def new_suffix_file(old_input_filename, old_input_filename_suffix,
                    new_filename_suffix, out_dir, file_extension=".txt"):
    # Strips the old file name of its suffix
    new_filename = \
        os.path.basename(old_input_filename) \
            .split(old_input_filename_suffix)[0] \
        + new_filename_suffix \
        + file_extension

    return os.path.join(out_dir, new_filename)


def get_organism_from_file(file_name, file_type_suffix):
    # Removes everything after the file type suffix to get the organisms name
    organism_name = os.path.basename(file_name).split(file_type_suffix)[0]

    return organism_name


def get_organism_genome_from_file(file_name, file_type_suffix):
    organism_complete_name = \
        os.path.basename(file_name).split("_" + file_type_suffix)[0]
    organism_name, genome_numeric_id = re.search("^(.*)_gen([0-9]+)$",
                                                 organism_complete_name).groups()
    return organism_name, int(genome_numeric_id)


def get_organism_id_from_name(organism_name, genome_numeric_id):
    return organism_name + "_gen" + str(genome_numeric_id)


def get_coord_file_from_id(organism_id,
                           parsed_coordinates_filename_suffix,
                           gene_locations_filename_suffix,
                           out_dir):
    # Location of the gene coordinates output directory
    parsed_coordinates_dir = os.path.join(out_dir,
                                          parsed_coordinates_filename_suffix)

    # Location of the gene coordinates file
    coordinates_filename = organism_id + parsed_coordinates_filename_suffix + ".json.gz"
    coordinates_file = os.path.join(parsed_coordinates_dir,
                                    coordinates_filename)

    # Location of the gene <-> chromosome mapping_tables file
    gene_to_chr_mapping_file = organism_id + gene_locations_filename_suffix + ".json.gz"
    gene_to_chr_mapping_file = os.path.join(parsed_coordinates_dir,
                                            gene_to_chr_mapping_file)

    return coordinates_file, gene_to_chr_mapping_file
