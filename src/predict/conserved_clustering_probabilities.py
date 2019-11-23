import numpy as np
from utilities import file_utilities
import predict.GSC as GSC
import predict.manage_probability_files as prob_files_utilities

import logging

logger = logging.getLogger(__name__)


def compute_conserved_probabilities(probabilities_files_list,
                                    organisms_uniform_probabilities,
                                    n_pairs,
                                    tree,
                                    subclade_id,
                                    subclade_final_probabilities_filename,
                                    args):
    """ This function computes the finals conserved clustered probabilities of all
    organisms whose probabilities are listed in probabilities files list.
    """
    if subclade_id == 0:
        subclade_id = "root"

    if args["verbose"] is False:
        logger.setLevel(logging.INFO)

    logger.verbose(
        "Computing final conserved clustering probabilities from subclade %s",
        subclade_id)

    # Initialises empty arrays for the results
    domain_pairs_probabilities = {'exp_prob': np.zeros(n_pairs),
                                  'uni_prob': np.zeros(n_pairs)}

    domain_pairs_data = {'occurrences': np.zeros(n_pairs, dtype=np.uint32),
                         'n_genes1': np.zeros(n_pairs, dtype=np.uint32),
                         'n_genes2': np.zeros(n_pairs, dtype=np.uint32),
                         'fusions': np.zeros(n_pairs, dtype=np.uint16),
                         'n_organisms': np.zeros(n_pairs, dtype=np.uint16)}

    if args["weight"] is True:
        organism_weights = GSC.GSC_normalised(tree)
        GSC.save_weights(organism_weights,
                         args["weights_file"],
                         subclade_id,
                         "a")

    for probabilities_filename in probabilities_files_list:
        logger.verbose(
            "Adding individual clustering probabilities from file %s",
            probabilities_filename)

        organism_name = file_utilities.get_organism_from_file(
            probabilities_filename,
            args["probabilities_filename_suffix"])

        organism_has_pair = np.zeros(n_pairs, dtype=np.uint8)

        if args["weight"] is True:
            organism_weight = organism_weights[organism_name]
        else:
            organism_weight = 1

        organism_probabilities = file_utilities.load_json(
            probabilities_filename)
        for genome_key in organism_probabilities.keys():
            add_organism_probabilities(organism_probabilities[genome_key],
                                       organisms_uniform_probabilities[
                                           genome_key],
                                       organism_weight,
                                       domain_pairs_probabilities,
                                       domain_pairs_data,
                                       organism_has_pair)

        domain_pairs_data["n_organisms"] = \
            domain_pairs_data["n_organisms"] + organism_has_pair
        del organism_has_pair

        logger.verbose("Individual clustering probabilities from file %s added",
                       probabilities_filename)

    domain_pairs_data["prob"] = compute_final_probabilities(
        domain_pairs_probabilities, args["phi"])
    del domain_pairs_probabilities

    logger.verbose(
        "Saving final conserved clustering probabilities from subclade %s",
        subclade_id)

    prob_files_utilities.\
        save_conserved_probabilities(domain_pairs_data,
                                     subclade_final_probabilities_filename)

    logger.verbose(
        "Final conserved clustering probabilities from "
        "subclade %s saved in file %s",
        subclade_id, subclade_final_probabilities_filename)


def add_organism_probabilities(chromosome_probabilities,
                               uniform_probability,
                               organism_weight,
                               domain_pairs_probabilities,
                               domain_pairs_data,
                               organism_has_pair):
    np.add.at(domain_pairs_probabilities['exp_prob'],
              chromosome_probabilities['pair_id'],
              np.array(chromosome_probabilities['exp_prob']) * organism_weight)
    np.add.at(domain_pairs_probabilities['uni_prob'],
              chromosome_probabilities['pair_id'],
              np.array(chromosome_probabilities['occurrences'])
              * uniform_probability * organism_weight)

    for key in ['occurrences', 'fusions', 'n_genes1', 'n_genes2']:
        np.add.at(domain_pairs_data[key],
                  chromosome_probabilities['pair_id'],
                  chromosome_probabilities[key])

    organism_has_pair[chromosome_probabilities['pair_id']] = 1


def compute_final_probabilities(domain_pairs_probabilities, phi):
    """This function puts together all the probabilities according to
    the final prob formula
    """
    np.seterr(divide='ignore', invalid='ignore')
    nom = phi * np.exp(domain_pairs_probabilities['exp_prob'])
    uni = (1 - phi) * np.exp(domain_pairs_probabilities['uni_prob'])
    final_probabilities = nom / ( nom + uni)
    mask = np.isnan(final_probabilities)
    final_probabilities[mask] = \
        (1 * (np.log(phi) + domain_pairs_probabilities['exp_prob'][mask]) >
         (np.log(1-phi) + domain_pairs_probabilities['uni_prob'][mask]))

    np.seterr(divide='warn', invalid='warn')
    return final_probabilities
