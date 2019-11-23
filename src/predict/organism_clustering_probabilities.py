import numpy as np
import os
from scipy import stats
import multiprocessing
from utilities import file_utilities
from utilities.mymultiprocess import MyMultiProcess
from genome_parsing import load_genome_annotations
from utilities.em_algorithm import log_expon_pdf
import genome_parsing.genomic_pairwise_distances as domain_pairs_distances
import predict.indexing_domain_pairs as indexing_domain_pairs
import predict.manage_probability_files as probabilities_files
import logging

logger = logging.getLogger(__name__)
lock = multiprocessing.Lock()
manager = multiprocessing.Manager()


def compute_organisms_probabilities(args):
    """This function is an aggregator to compute the exponential and uniform likelihood
    of all single genomes, using multiprocess.
    """
    pfam_families = file_utilities.load_json(args["pfam_families_filename"])

    pfam_annotations_file_list = args["pfam_annotations_file_list"]
    n_organisms = len(pfam_annotations_file_list)

    if not os.path.exists(args["uniform_probabilities_filename"]) or \
            args["force"] is True:
        uniform_probabilities = manager.dict()
    else:
        uniform_probabilities = manager.dict(probabilities_files.
            load_uniform_probabilities(
            args["uniform_probabilities_filename"]))

    number_of_pairs = indexing_domain_pairs.index_all_pairs(
        pfam_annotations_file_list,
        pfam_families, args)

    logger.info("Started computing organism-specific clustering "
                "probabilities of Pfam domains")

    multi = MyMultiProcess(threads=args["threads"],
                           target=organism_domain_clustering_probabilities,
                           input=pfam_annotations_file_list,
                           args=[uniform_probabilities,
                                 pfam_families,
                                 args])
    multi.run()

    logger.info("Computation of organism-specific clustering probabilities of "
                "Pfam domains for %s organisms done",
                n_organisms)
    logger.info(
        "All individual organism probabilities were saved in directory %s",
        args["individual_probabilities_dir"])

    return uniform_probabilities, number_of_pairs


def organism_domain_clustering_probabilities(pfam_annotations_file,
                                             uniform_probabilities,
                                             pfam_families,
                                             args):
    """This function computes the exp and unif likelihood
    for all pairs from a single genome."""

    organism_name = \
        file_utilities.get_organism_from_file(pfam_annotations_file,
                                              args[
                                                  "pfam_annotations_filename_suffix"])

    # Output file
    probabilities_filename = \
        os.path.join(args["individual_probabilities_dir"],
                     organism_name + args["probabilities_filename_suffix"]
                     + ".json")

    if not os.path.exists(probabilities_filename) or args['force'] is True:
        individual_indexes_filename = os.path.join(args["pairs_indexes_dir"],
                                                   args[
                                                       "pairs_indexes_individual_template"].
                                                   format(organism_name))

        all_domain_pair_indexes = file_utilities.load_json(
            individual_indexes_filename)

        # Copy step necessary, because the dictionary is shared among threads and the
        # following checks make the script slower if ran in multithread:
        genome_sizes_org = [args["genome_sizes"][organism_name][0]]
        genome_numeric_id = 0

        for genome_molecule in genome_sizes_org:
            genome_key = \
                file_utilities.get_organism_id_from_name(organism_name,
                                                         genome_numeric_id)

            logger.verbose("Started computing individual clustering "
                           "probabilities from file %s",
                           pfam_annotations_file)

            organism_probabilities = {}

            genome_size, genome_id = genome_molecule[1], genome_molecule[0]

            try:
                gene_coordinates
            except NameError:
                # Loads gene coordinates for this chromosome
                gene_coordinates, gene_to_chr_mapping = \
                    load_genome_annotations.import_gene_coord_data(
                        organism_name, args)

            # Import Pfam annotations for this chromosome
            domain_hits = \
                load_genome_annotations.import_pfam_hits(pfam_annotations_file,
                                                         gene_to_chr_mapping,
                                                         args["evalue"],
                                                         pfam_families,
                                                         genome_id)

            if len(domain_hits) > 0:
                # Saves the exponential probabilities for this chromosome
                organism_probabilities[genome_key] = \
                    compute_probabilities(domain_hits,
                                          all_domain_pair_indexes,
                                          gene_coordinates,
                                          genome_id,
                                          genome_size,
                                          pfam_annotations_file,
                                          args)

                # Computes the uniform probability for this chromosome
                uniform_probabilities[genome_key] = \
                    compute_uniform_probability(float(genome_size) / 2)

            genome_numeric_id += 1

        file_utilities.save_json(organism_probabilities,
                                 probabilities_filename)
        probabilities_files \
            .save_uniform_probabilities(uniform_probabilities,
                                        args["uniform_probabilities_filename"])
        logger.verbose("Individual clustering probabilities from file %s, "
                       "DNA molecule %s computed",
                       pfam_annotations_file, genome_id)
        del all_domain_pair_indexes


    else:
        logger.verbose("Individual clustering probabilities "
                       "from file %s already present",
                       pfam_annotations_file)


def compute_probabilities(domain_hits,
                          all_domain_pair_indexes,
                          gene_coordinates,
                          genome_id,
                          genome_size,
                          pfam_annotations_file,
                          args):
    logger.verbose("Analysing chromosome record %s from file %s "
                   "(whole assembly contains %s total unique domain pairs)",
                   genome_id, pfam_annotations_file, '{:,}'
                   .format(len(all_domain_pair_indexes)))

    # Initialises arrays for the results
    domain_pairs_data = {'pair_id': [],
                         'exp_prob': [],
                         'occurrences': [],
                         'fusions': [],
                         'n_genes1': [],
                         'n_genes2': []}

    for domain_pair in indexing_domain_pairs. \
            generate_domain_pairs_set(domain_hits.keys()):
        # Indexes the domain pair to speed up the future comparison across organisms
        domain_pairs_data['pair_id'].append(indexing_domain_pairs
                                            .get_domain_pair_index(domain_pair,
                                                                   all_domain_pair_indexes))

        pair_distances, n_fusions, n_genes1, n_genes2 = domain_pairs_distances. \
            compute_closest_domain_pair_distances(domain_pair,
                                                  domain_hits,
                                                  gene_coordinates,
                                                  genome_size)

        domain_pairs_data['exp_prob'] \
            .append(compute_log_exponential_probability(pair_distances,
                                                        args["lambd"]))

        domain_pairs_data['occurrences'].append(len(pair_distances))
        domain_pairs_data['n_genes1'].append(n_genes1)
        domain_pairs_data['n_genes2'].append(n_genes2)
        domain_pairs_data['fusions'].append(n_fusions)

    return domain_pairs_data


def compute_log_exponential_probability(domain_distances, lambd):
    # Computes the sum of the log exponential probability of the
    # pairwise distances for this pair of domains
    # (going to a minimum of -30 to avoid numeric instability)
    exp_probs = log_expon_pdf(domain_distances, lambd)
    summed_log_exp_probs = np.sum(exp_probs)

    return summed_log_exp_probs


def compute_uniform_probability(half_genome_size, domain_distance=100):
    # Computes the uniform probability for any pair of domains on this chromosome
    uni_prob = stats.uniform.pdf(domain_distance, 0, half_genome_size)

    # Computes the log of the uniform probability
    if uni_prob == 0:
        uni_prob = 10 ** (-20)
    log_uni_prob = np.log(uni_prob)

    return log_uni_prob
