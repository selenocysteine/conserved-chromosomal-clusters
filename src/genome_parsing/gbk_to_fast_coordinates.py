""""Parse gbk files to obtain a parsed coordinates file (main chromosome only)."""
import os
from Bio import SeqIO
from operator import itemgetter
import collections
from multiprocessing import Lock, Manager
from utilities import file_utilities
from utilities.mymultiprocess import MyMultiProcess
lock = Lock()
manager = Manager()
import logging
logger = logging.getLogger(__name__)


def parse_coordinates(args, strand=False):
    logger.info("Started parsing assembly gene coordinates files")

    # Create a dictionary for the sizes of each genome.
    # If one is already present, import it and append any new items.
    if os.path.exists(args["genome_sizes_filename"]) and args['force'] is False:
        try:
            genome_sizes_dictionary = \
                manager.dict(file_utilities.load_json(args["genome_sizes_filename"]))
        except EOFError:
            genome_sizes_dictionary = manager.dict()
    else:
        genome_sizes_dictionary = manager.dict()

    if len(genome_sizes_dictionary) >= len(args['gbk_files_list']):
        del args["gbk_files_list"]
        return genome_sizes_dictionary
    else:
        # Checks if everything has already been parsed
        args["parsed_coordinates_file_list"] = \
            file_utilities.get_file_list(args['parsed_coordinates_out_dir'],
                                     args["parsed_coordinates_filename_suffix"],
                                     verbose=True, error=False)

        # If there are files to process, starts to parse them:
        if len(args["gbk_files_list"]) > 0:
            parsed_coordinates_file_list = manager.list()

            multi = MyMultiProcess(threads=args["threads"],
                                   target=create_parsed_coordinate_dictionary,
                                   input=args["gbk_files_list"],
                                   args=[genome_sizes_dictionary,
                                         parsed_coordinates_file_list,
                                         strand,
                                         args])
            multi.run()

        if len(args["gbk_files_list"]) > 0:
            logger.info("Parsed gene coordinates files for %s assemblies generated and saved in folder %s",
                        len(args["gbk_files_list"]),
                        args['parsed_coordinates_out_dir'])
        logger.info("Individual DNA molecule sizes saved in %s",
                    args['genome_sizes_filename'])

        args["parsed_coordinates_file_list"] = list(parsed_coordinates_file_list)

    return genome_sizes_dictionary


def create_parsed_coordinate_dictionary(gbk_file,
                                        genome_sizes_dictionary,
                                        parsed_coordinates_file_list,
                                        strand,
                                        args):
    # Output file name for the coordinates of the genes on the main chromosome (json object)
    chromosomal_genes_coord_file = \
        file_utilities.new_suffix_file(gbk_file, ".gbk",
                                       args["parsed_coordinates_filename_suffix"],
                                       args["parsed_coordinates_out_dir"],
                                       file_extension=".json.gz")

    # Output file name for a dictionary of gene_ID -> chromosomal molecule mappings (json object)
    gene_to_chr_mapping_file = \
        file_utilities.new_suffix_file(gbk_file, ".gbk",
                                       args["gene_locations_filename_suffix"],
                                       args["parsed_coordinates_out_dir"],
                                       file_extension=".json.gz")

    # If all of these conditions are already satisfied, we do not need to generate the coordinates files
    # for this organism, because they have already been parsed before
    organism_id = file_utilities.get_organism_from_file(gbk_file, ".gbk")

    if not os.path.exists(chromosomal_genes_coord_file) \
            or (args['force'] is True) \
            or not os.path.exists(gene_to_chr_mapping_file) \
            or not organism_id in genome_sizes_dictionary.keys():
        logger.verbose("Started parsing gbk file %s", gbk_file)
        genomes_info, chromosomal_genes_coordinates, gene_to_chr_mapping = parse_single_gbk_for_coordinates(gbk_file, strand)

        # Saves the output
        update_genome_sizes_dictionary(genome_sizes_dictionary,
                                       gbk_file,
                                       genomes_info,
                                       args["genome_sizes_filename"])
        file_utilities.save_json(gene_to_chr_mapping, gene_to_chr_mapping_file, compress=True)
        file_utilities.save_json(chromosomal_genes_coordinates, chromosomal_genes_coord_file, compress=True)

        logger.verbose("Gbk file %s parsed", gbk_file)
        parsed_coordinates_file_list.append(chromosomal_genes_coord_file)

    else:
        logger.verbose("Gbk file %s already parsed", gbk_file)


def parse_single_gbk_for_coordinates(gbk_file, strand):
    gene_coordinates = collections.OrderedDict()
    gene_to_chr_mapping = {}
    genome_info = []
    last_dna_molecule = ""
    for gbk_record in SeqIO.parse(open(gbk_file, "r"), "genbank"):
        for seq_feature in gbk_record.features:
            if seq_feature.type == "CDS":
                gene_ID = seq_feature.qualifiers["locus_tag"][0]
                # Save the name of the DNA molecule this gene is located on
                gene_to_chr_mapping[gene_ID] = gbk_record.id

                gene_start = int(seq_feature.location.start)
                gene_end = int(seq_feature.location.end)
                gene_coordinates[gene_ID] = [gene_start, gene_end]
                if strand is True:
                    gene_coordinates[gene_ID].append(int(seq_feature.location.strand))

                if gbk_record.id != last_dna_molecule:
                    last_dna_molecule = gbk_record.id
                    genome_info.append([gbk_record.id, int(len(gbk_record))])

    # Sorts the list so that the longest molecule is at the start (for the fitting)
    genome_info = sorted(genome_info, key=itemgetter(1), reverse=True)

    return genome_info, gene_coordinates, gene_to_chr_mapping


def update_genome_sizes_dictionary(genome_sizes_dictionary,
                                   gbk_file,
                                   genome_sizes,
                                   genome_sizes_filename):
    organism_id = file_utilities.get_organism_from_file(gbk_file, ".gbk")

    genome_sizes_dictionary[organism_id] = genome_sizes

    with lock:
        file_utilities.save_json(dict(genome_sizes_dictionary), genome_sizes_filename, compress=True)

    return 0
