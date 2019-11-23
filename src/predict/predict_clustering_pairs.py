import os
import ete3
import predict.organism_clustering_probabilities as clustering_probs
import predict.conserved_clustering_probabilities as final_probs
import predict.subclade_analysis as subclade_analysis
import predict.manage_probability_files as probabilities_files
import predict.GSC as GSC
import predict.find_last_subclade as find_last_subclade
from utilities import file_utilities

import logging

logger = logging.getLogger(__name__)


def predict_clustering_pairs(args, **kwargs):
    """This function is a wrapper for all the steps for the computation of the
    final clustering probabilities.
    """
    # We first compute the exponential and uniform likelihoods for all genomes
    # in the dataset, individually
    organisms_uniform_probabilities, n_pairs = \
        clustering_probs.compute_organisms_probabilities(
            args, **kwargs)

    # Finds all the probabilities files for single organisms
    probabilities_files_list = file_utilities.get_file_list(
        args["individual_probabilities_dir"],
        args["probabilities_filename_suffix"])

    if args["weight"] is True or args["subclades"] is True:
        logger.info("Importing phylogenetic tree %s", args["tree_file"])
        tree = ete3.Tree(args["tree_file"])
        logger.info("Phylogenetic tree %s imported", args["tree_file"])
        if args["weight"] is True:
            if not os.path.exists(args["weights_file"]) or args[
                "force"] is True:
                GSC.save_weights(None, args["weights_file"], None, "w")
    else:
        tree = None

    # Conserved clustering analysis in each subclade:
    if args["subclades"] is True:
        logger.info(
            "Started computing subclade-specific conserved"
            " clustering probabilities")
        subclade_analysis.partition_tree(probabilities_files_list,
                                         organisms_uniform_probabilities,
                                         n_pairs,
                                         tree,
                                         args,
                                         **kwargs)
        logger.info(
            "Parsing subclade-specific conserved clustering probabilities, "
            "this will take a while")
        find_last_subclade.parse_and_clean_tree(n_pairs, args)
        logger.info(
            "Parsing subclade-specific conserved clustering probabilities done")
        logger.info(
            "Subclade-specific conserved clustering probabilities computed")


    # Clustering analysis just in the root:
    else:
        if not os.path.exists(args["final_probabilities_filename"]) \
                or probabilities_files.is_table_empty(
            args["final_probabilities_filename"]) \
                or args["force"] is True:
            logger.info(
                "Started computing root-only conserved clustering probabilities")
            final_probs.compute_conserved_probabilities(
                probabilities_files_list,
                organisms_uniform_probabilities,
                n_pairs,
                tree,
                0,
                args["final_probabilities_filename"],
                args)

    logger.info("Root conserved clustering probabilities saved in file %s",
                args["final_probabilities_filename"])

    if args["weight"] is True:
        logger.info("Phylogenetic weights saved in file %s",
                    args["weights_file"])


    return 0
