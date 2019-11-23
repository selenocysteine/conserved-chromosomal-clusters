import numpy as np
import ete3
from utilities import import_export_parameters
import fit.estimate_clustering_parameters as estimate_parameters
import predict.GSC as GSC
import logging

logger = logging.getLogger(__name__)


def fit_clustering_model(args):
    logger.info("Started estimating clustering model parameters "
                "for each organism")
    lambdas, phis = estimate_parameters.estimate_individual_parameters(args)

    logger.info("Computing global mean estimated parameter values")

    if args["weight"] is True:
        logger.info("Importing phylogenetic tree %s", args["tree_file"])
        tree = ete3.Tree(args["tree_file"])
        logger.info("Phylogenetic tree %s imported", args["tree_file"])
        organism_weights = GSC.GSC_normalised(tree)
        lambd = estimate_parameters.compute_weighted_means(lambdas,
                                                           organism_weights)
        phi = estimate_parameters.compute_weighted_means(phis,
                                                         organism_weights)
    else:
        lambd = np.mean(lambdas.values())
        phi = np.mean(phis.values())

    args["lambd"], args["phi"] = lambd, phi
    logger.info("Global mean estimated parameter "
                "values: lambda={:.3g}, phi={:.3g}".format(
        args["lambd"], args["phi"]))
    import_export_parameters. \
        save_general_parameters(lambd, phi, args["general_parameters_filename"])
    logger.info("Global mean estimated parameter values saved in file %s",
                args["general_parameters_filename"])

    return lambd, phi
