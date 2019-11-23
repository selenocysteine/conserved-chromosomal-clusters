import os
import multiprocessing
import logging
import numpy as np
from multiprocessing import Pool
from collections import defaultdict
import traceback
import random
from predict.subclade_analysis import load_labelled_tree, extract_subclades
import predict.GSC as GSC
from genome_parsing.genomic_pairwise_distances import distance_between_two_genes
from utilities.em_algorithm import log_expon_pdf
from permutations import distance_matrix, data_wrangling, permutations_plots
logger = logging.getLogger(__name__)
manager = multiprocessing.Manager()


def predict_permuted_clustering_pairs(args):
    n_perm = 10 ** 4
    random.seed(args["random_seed"])

    logger.info("Started significance analysis. Number of simulations: {}"
                .format(n_perm))
    # Loads subclade data and computes distance matrix
    tree = load_labelled_tree(args["labelled_tree_file"])
    subclades_IDs, subclades = extract_subclades(tree, args)

    all_organisms = np.array(tree.get_leaf_names())
    all_organisms.sort()

    distances = \
        distance_matrix.compute_distance_matrix(tree, args)

    # Loads gene coordinates for all organisms
    all_organisms_coordinates = \
        data_wrangling.load_all_organisms_coordinates(args)

    # Initialises output file
    if not os.path.exists(args['permutations_file']) or args['force'] is True:
        with open(args['permutations_file'], 'w') as out_f:
            out_f.write('subclade_id\tn_members\toccurrences'
                        '\te-4\te-3\tmean\tmedian\tpct_above_0.5'
                        '\tpct_migrated_pairs\n')
        done_permutations_all = defaultdict(list)
    else:
        done_permutations_all = \
            data_wrangling.load_permutations(args['permutations_file'])

    # Starts the simulations in multiprocessing
    errors = manager.list()
    pool = Pool(processes=args['threads'],
                maxtasksperchild=10)
    for subclade in zip(subclades.keys(), subclades.values()):
        subclade_members = np.array(subclade[1].get_leaf_names())
        subclade_members.sort()
        all_gene_coordinates = \
            {key: all_organisms_coordinates[key] for key
             in subclade_members}

        these_distances = \
            distance_matrix.filter_data_table(distances, subclade_members)
        done_permutations = done_permutations_all[subclade[0]]

        pool.apply_async(wrapper_worker,
                         args=(subclade,
                               subclade_members,
                               all_gene_coordinates,
                               these_distances,
                               n_perm,
                               done_permutations,
                               errors,
                               args))
    pool.close()
    pool.join()
    for error in errors:
        logger.error(error)



def compute_randomised_probabilities(occurrences,
                                     organism_weight,
                                     organism_genome_data,
                                     organism_id,
                                     first_distance,
                                     first_organism,
                                     distances,
                                     lambd):
    # Simulates a pair present within organism 'organism_id'
    log_uni_prob = organism_genome_data[0]
    coordinates = organism_genome_data[2]
    probabilities = np.array([0,
                              log_uni_prob * organism_weight * occurrences],
                             dtype=float)

    genome_size = organism_genome_data[1]
    movements = 0
    for occurrence in range(occurrences):
        if first_distance is None:
            not_same = 0
        else:
            if organism_id != first_organism:
                organisms = sorted([organism_id, first_organism])
                distance = distances[organisms[0], organisms[1]]
                prob_same = float(1 - distance) / occurrences
            else:
                prob_same = 1. / occurrences

            not_same = 1 if random.random() > prob_same else 0

        if not_same == 1:
            movements += 1
            domain1, domain2 = random.randint(0, coordinates.shape[0] - 1), \
                               random.randint(0, coordinates.shape[0] - 1)
            domain1_coords = coordinates[domain1]
            domain2_coords = coordinates[domain2]
            domain_dist = \
                distance_between_two_genes(domain1_coords,
                                           domain2_coords, genome_size)
        else:
            domain_dist = \
                sample_adjacent_genes_distances(coordinates)
        if first_distance is None:
            first_distance = domain_dist

        exp_pdf = log_expon_pdf(domain_dist, lambd)
        probabilities[0] += exp_pdf

    probabilities[0] = probabilities[0] * organism_weight

    return probabilities, first_distance, movements


def permute_bin(occurrence, subclade_id, n_organisms,
                subclade_members, organism_weights,
                all_gene_coordinates, distances,
                permuted_probabilities, moved_pairs_perm,
                n_permutations, done_permutations, args):
    if occurrence not in done_permutations:
        logger.verbose(
            "{} permutations for subclade {}, {} occurrences"
                .format(n_permutations,
                        subclade_id,
                        occurrence))
        if n_organisms > 1:
            pairs_per_organism_all = \
                randomly_pick_members(occurrence, n_organisms)

        for perm_index in range(0, n_permutations):
            if n_organisms > 1:
                pairs_per_organism = np.array(pairs_per_organism_all.pop())
                pairs_per_organism = pairs_per_organism[pairs_per_organism > 0]
                n_organisms_with_pair = pairs_per_organism.shape[0]
                if n_organisms_with_pair < 50:
                    organism_indexes = random.sample(range(0, n_organisms),
                                                     n_organisms_with_pair)
                else:
                    organism_indexes = \
                        np.random.randint(low=0,
                                          high=n_organisms - 1,
                                          size=n_organisms_with_pair)
            else:
                pairs_per_organism = np.array([occurrence])
                organism_indexes = np.array([0])
            probabilities = np.array([0, 0])
            first_pair = None
            first_organism = organism_indexes[0]
            for index in range(len(pairs_per_organism)):
                organism_id = organism_indexes[index]
                organism_name = subclade_members[organism_id]
                if args['weight'] is True:
                    organism_weight = organism_weights[organism_id]
                else:
                    organism_weight = 1
                organism_occurrences = pairs_per_organism[index]
                prob2, first_pair, not_same = \
                    compute_randomised_probabilities(
                        organism_occurrences,
                        organism_weight,
                        all_gene_coordinates[organism_name],
                        organism_id,
                        first_pair,
                        first_organism,
                        distances,
                        args['lambd'])
                moved_pairs_perm[perm_index] += not_same
                probabilities = probabilities + prob2
            permuted_probabilities[perm_index] = \
                conserved_clustering_probabilities(probabilities[0],
                                                   probabilities[1],
                                                   args['phi'])
            moved_pairs_perm[perm_index] = \
                moved_pairs_perm[perm_index] / occurrence

        data_wrangling.update_results_file(permuted_probabilities,
                                           moved_pairs_perm,
                                           subclade_id,
                                           subclade_members,
                                           occurrence,
                                           args)
        logger.verbose("Subclade {}, {} occurrences done"
                       .format(subclade_id, occurrence))


def run_permutations_on_file(subclade,
                             subclade_members,
                             all_gene_coordinates,
                             distances,
                             n_permutations,
                             done_permutations,
                             args):
    subclade_id = subclade[0]
    tree = subclade[1]

    n_organisms = len(subclade_members)
    if args["weight"] is True:
        organism_weights = np.zeros(n_organisms)
        temp_organism_weights = GSC.GSC_normalised(tree)
        for index, org in enumerate(subclade_members):
            organism_weights[index] = temp_organism_weights[org]
        del temp_organism_weights

    occurrences = sorted(list(data_wrangling.count_occurrences
                              (subclade_id, args)))
    logger.info("Analysing subclade {} {} occurrences".format(subclade_id,
                                                              occurrences))

    permuted_probabilities = np.zeros(n_permutations, dtype=float)
    moved_pairs_perm = np.zeros(n_permutations, dtype=float)
    for occurrence in occurrences:
        permute_bin(occurrence, subclade_id, n_organisms,
                    subclade_members, organism_weights,
                    all_gene_coordinates, distances,
                    permuted_probabilities, moved_pairs_perm,
                    n_permutations, done_permutations, args)

    logger.info("Subclade {} {} occurrences analysed".format(subclade_id,
                                                             occurrences))
    return 0


def sample_adjacent_genes_distances(coordinates):
    domain, domain2 = random.randint(0, coordinates.shape[0] - 1), \
                      random.randint(0, coordinates.shape[0] - 1)
    distance1 = coordinates[domain][1] - coordinates[domain][0] + 1
    distance2 = coordinates[domain2][1] - coordinates[domain2][0] + 1
    return distance1 + distance2


def randomly_pick_members(occurrences, n_organisms):
    samples = np.array(np.random.multinomial(occurrences,
                                             [1. / n_organisms] *
                                             n_organisms,
                                             size=10000))

    return samples.tolist()


def conserved_clustering_probabilities(log_exp_prob, log_uni_prob, phi):
    np.seterr(divide='ignore', invalid='ignore')
    nom = (phi * np.exp(log_exp_prob))
    uni = (1 - phi) * np.exp(log_uni_prob)
    final_probabilities = np.array(nom / (nom + uni))
    if np.isnan(final_probabilities):
        if (np.log(phi) + log_exp_prob) > (np.log(1 - phi) + log_uni_prob):
            final_probabilities = np.array([1.], dtype=float)
        else:
            final_probabilities = np.array([0.], dtype=float)
    np.seterr(divide='warn', invalid='warn')

    return final_probabilities


def wrapper_worker(subclade,
                   all_organisms,
                   all_gene_coordinates,
                   these_distances,
                   n_perm,
                   done_permutations,
                   errors,
                   args):
    try:
        run_permutations_on_file(subclade,
                                 all_organisms,
                                 all_gene_coordinates,
                                 these_distances,
                                 n_perm,
                                 done_permutations,
                                 args)
    except Exception as e:
        errors.append([subclade[0], 'run_bootstrap_on_file', e])
        traceback.print_exception(type(e), e, e.__traceback__)
        logger.warning(
            "Thread crashed on subclade {}: {}\n".format(subclade[0], e))
        raise ValueError
