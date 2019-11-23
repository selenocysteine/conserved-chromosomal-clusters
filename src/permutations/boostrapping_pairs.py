import os
import multiprocessing
import logging
import numpy as np
from multiprocessing import Pool
from collections import defaultdict
import traceback
import random
from predict.subclade_analysis import load_labelled_tree, extract_subclades
from permutations import data_wrangling
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
    args['bootstrapping_pairs_file'] = os.path.join(args['permutations_dir'],
                                                     'bootstrap_results.txt')
    # Initialises output file
    if not os.path.exists(args['bootstrapping_pairs_file']) or args['force'] is True:
        with open(args['bootstrapping_pairs_file'], 'w') as out_f:
            out_f.write('subclade_id\tn_members\toccurrences'
                        '\te-4\te-3\tmean\tmedian\tpct_above_0.5\n')
        done_permutations_all = defaultdict(list)
    else:
        done_permutations_all = \
            data_wrangling.load_permutations(args['bootstrapping_pairs_file'])

    # Starts the simulations in multiprocessing
    errors = manager.list()
    pool = Pool(processes=args['threads'],
                maxtasksperchild=10)
    for subclade in zip(subclades.keys(), subclades.values()):
        done_permutations = done_permutations_all[subclade[0]]
        pool.apply_async(wrapper_worker,
                         args=(subclade[0],
                               len(subclade[1].get_leaf_names()),
                               done_permutations,
                               errors,
                               args))
    pool.close()
    pool.join()
    for error in errors:
        logger.error(error)



def bootstrap_bin(occurrence,
                  subclade_id,
                  n_occurrences,
                  n_members,
                  done_bootstraps, args):
    n_samples = 10000
    n_occurrences_probs = occurrence.shape[0]
    bootstrap_means = np.zeros(n_samples)
    if n_occurrences not in done_bootstraps:
        logger.verbose(
            "Bootstrap for subclade {}, {} occurrences"
                .format(subclade_id,
                        n_occurrences))
        for i in range(0, n_samples):
            random_indexes = np.random.randint(0, n_occurrences_probs,
                                               n_occurrences_probs)
            data_bootstrap = occurrence[random_indexes]
            bootstrap_means[i] = np.sum(data_bootstrap[data_bootstrap > 500])/data_bootstrap.size

        data_wrangling.update_bootstrap_results_file(bootstrap_means,
                                           subclade_id,
                                           n_members,
                                           n_occurrences,
                                           args)
        logger.verbose("Bootstrap subclade {}, {} occurrences done"
                       .format(subclade_id, n_occurrences))


def run_bootstrap_on_file(subclade_id,
                          n_organisms,
                          done_bootstraps,
                          args):

    occurrences = data_wrangling.load_occurrences_probs(subclade_id, args)
    all_occurrences = sorted(list(occurrences.keys()))
    logger.info("Analysing subclade {} {} occurrences"
                .format(subclade_id, all_occurrences))

    for occurrence in occurrences:
        if len(occurrences[occurrence]) >= 100:
            bootstrap_bin(np.array(occurrences[occurrence], dtype=np.int32),
                          subclade_id,
                          occurrence,
                          n_organisms,
                          done_bootstraps,
                          args)
        else:
            with open(args['bootstrapping_pairs_file'], "a") as out_f:
                out_f.write(
                    "{}\t{}\t{}\t{}\t{}"
                    "\t{}\t{}\t{}\n".format(subclade_id,
                                                        n_organisms,
                                                        occurrence,
                                                        'NA',
                                                        'NA',
                                                        'NA',
                                                        'NA',
                                                        'NA'))

    logger.info("Subclade {} {} occurrences analysed"
                .format(subclade_id, all_occurrences))

    return 0


def wrapper_worker(subclade,
                   n_organisms,
                   done_permutations,
                   errors,
                   args):
    try:
        run_bootstrap_on_file(subclade,
                              n_organisms,
                              done_permutations,
                              args)
    except Exception as e:
        errors.append([subclade, 'run_bootstrap_on_file', e])
        traceback.print_exception(type(e), e, e.__traceback__)
        logger.warning(
            "Thread crashed on subclade {}: {}\n".format(subclade, e))
        raise ValueError
