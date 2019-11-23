import matplotlib
matplotlib.use('Agg')
import os
import numpy as np
from math import log10
from decimal import Decimal
from utilities import file_utilities
from utilities import import_export_parameters
import utilities.em_algorithm as emCC
from utilities.mymultiprocess import MyMultiProcess
import fit.GO_gene_pairs_distances as GO_gene_pairs_distances
import fit.plot_functions as plot_func
import multiprocessing
import logging

manager = multiprocessing.Manager()
logger = logging.getLogger(__name__)


def estimate_individual_parameters(args):
    pfam_annotations_file_list2 = list()

    # If there are Pfam annotations files to process:
    if len(args["pfam_annotations_file_list"]) > 0:
        if os.path.exists(args['organisms_parameters_filename']) and \
                args['force'] is False:
            lambdas, phis = \
                import_export_parameters.import_organisms_parameters(
                    args['organisms_parameters_filename'])
            lambdas, phis = manager.dict(lambdas), manager.dict(phis)
            if len(lambdas) != len(args["pfam_annotations_file_list"]):
                for pfam_annotations_file in args["pfam_annotations_file_list"]:
                    organism_name = file_utilities.get_organism_from_file(
                        pfam_annotations_file,
                        args["pfam_annotations_filename_suffix"])
                    try:
                        lambdas[organism_name]
                    except KeyError:
                        pfam_annotations_file_list2.append(
                            pfam_annotations_file)
            else:
                return lambdas, phis
        else:
            # File to save the estimated parameters for each organism
            args["organisms_parameters_filename"] = import_export_parameters. \
                initialise_organisms_parameters_file(
                args["organisms_parameters_filename"])
            pfam_annotations_file_list2 = args['pfam_annotations_file_list']
            lambdas, phis = manager.dict(), manager.dict()
            # for file in args["pfam_annotations_file_list"]:
            #     if "GCA_000008765.1_ASM876v1" in file:
            #         pfam_annotations_file_list2 = [file]

        pfam_families, GO_to_pfam_mapping = load_GO_annotations_mapper(args)

        logger.info(
            "Started fitting model on chromosomal distances of"
            " GO-related Pfam domain pairs")
        multi = MyMultiProcess(threads=args["threads"],
                               target=fit_on_single_organism,
                               input=pfam_annotations_file_list2,
                               args=[GO_to_pfam_mapping,
                                     pfam_families,
                                     lambdas,
                                     phis,
                                     args],
                               maxtasksperchild=1)
        multi.run()

    logger.info(
        "Fitting model on chromosomal distances of "
        "GO-related Pfam domain pairs done")
    logger.info("Estimated parameter values for all organisms saved in file %s",
                args["organisms_parameters_filename"])

    return lambdas, phis


def fit_on_single_organism(pfam_annotations_file,
                           GO_to_pfam_mapping,
                           pfam_families,
                           lambdas,
                           phis,
                           args):
    # Computes the pairwise distances of GO pairs on this genome
    organism_name, genome_size, GO_pairs_distances = \
        GO_gene_pairs_distances.import_GO_pairs_distances(pfam_annotations_file,
                                                          GO_to_pfam_mapping,
                                                          pfam_families,
                                                          args)

    if GO_pairs_distances.shape[0] == 0:
        logger.verbose(
            "No GO-related Pfam domain pairs for file %s available, "
            "parameters not computed",
            pfam_annotations_file)
        return 0

    else:
        lambd, phi = estimate_parameters_on_GO_distances(GO_pairs_distances,
                                                         pfam_annotations_file,
                                                         organism_name,
                                                         genome_size,
                                                         args)
        lambdas[organism_name] = lambd
        phis[organism_name] = phi

        # Updates the parameters file
        logger.verbose(
            "Saving organism-specific estimated parameter values for file %s",
            pfam_annotations_file)
        import_export_parameters \
            .save_organism_parameters(organism_name,
                                      float(genome_size) / 2,
                                      lambd, phi,
                                      args["organisms_parameters_filename"])
        logger.verbose(
            "Organism-specific estimated parameter values for file %s saved",
            pfam_annotations_file)


def estimate_parameters_on_GO_distances(GO_pairs_distances,
                                        pfam_annotations_file,
                                        organism_name,
                                        genome_size,
                                        args):
    # Estimate parameters on pairwise distances of GO-related gene pairs
    logger.verbose(
        "Started estimating parameters from GO-related "
        "Pfam domain pairs from file %s",
        pfam_annotations_file)

    # Estimates parameters of the exponential from the GO pairs
    half_genome_size = float(genome_size) / 2
    lambd, phi, likelihood = fit_mixture_model_EM(GO_pairs_distances,
                                                  half_genome_size,
                                                  args["phi_0"],
                                                  args["lambd_0"],
                                                  args["max_iter"],
                                                  args["epsilon"],
                                                  args["sample_initial_values"],
                                                  args["n_samples"])

    # Saves a plot of the estimated distribution after the fit
    # + convergence plot
    if args["plots"] is True:
        plot_file_name = os.path.join(args["em_plots_dir"],
                                      organism_name +
                                      "_fitting_results_plot.pdf")
        plot_func.plot_EM_fitted_distribution(lambd,
                                              phi,
                                              half_genome_size,
                                              GO_pairs_distances,
                                              likelihood,
                                              organism_name,
                                              plot_file_name,
                                              args)

        logger.verbose(
            "Convergence and estimated distribution plot "
            "for organism %s saved in %s",
            organism_name,
            plot_file_name)

    return lambd, phi


def fit_mixture_model_EM(dist_GO,
                         half_genome_size,
                         phi_0,
                         lambd_0,
                         max_iter,
                         epsilon,
                         sample_initial_values,
                         n_samples):
    """ Samples initial guesses for the parameters and fits the model."""
    if sample_initial_values is False:
        n_samples = 1

    # Initial guess for the parameter lambda
    if lambd_0 is None:
        lambd_0 = 1. / float(np.mean(dist_GO))
    magnitudo_lambda = int(log10(lambd_0))

    iteration = 0
    # Samples n_samples values of phi_0 logarithmically from the
    # initial value to 99% of the total pairs
    for phi_0 in emCC.log_sample_parameters(phi_0, 0.1, n_samples):
        # Samples logarithmically n_samples values of lambd_0
        # between 3 orders of magnitude smaller
        # and 3 orders of magnitude higher than the initial guess
        for lambd_0 in emCC.log_sample_parameters(10 ** (magnitudo_lambda - 3),
                                                  10 ** (magnitudo_lambda + 3),
                                                  n_samples):
            # Fits the model using the EM algorithm
            em_output = emCC.em(dist_GO,
                                lambd_0, phi_0,
                                half_genome_size,
                                max_iter=max_iter,
                                epsilon=epsilon)

            if iteration == 0:
                set_best_mllikelihood = True
                if int(em_output['phi']) == 1:
                    em_output['MLL'][-1] = "NA"
                iteration = 1
            else:
                if int(em_output["phi"]) != 1:
                    if best_mllikelihood[0] == "NA":
                        set_best_mllikelihood = True
                    else:
                        set_best_mllikelihood = \
                            Decimal(em_output["MLL"][-1]) > \
                            Decimal(best_mllikelihood[0])

            if set_best_mllikelihood is True:
                best_mllikelihood = [em_output["MLL"][-1],
                                     em_output["lambda"],
                                     em_output["phi"],
                                     em_output["MLL"]]
                set_best_mllikelihood = False


    lambd, phi, mllikelihoods = best_mllikelihood[1], \
                                best_mllikelihood[2], \
                                best_mllikelihood[-1]


    if phi == 1:
        print(em_output)
        print(best_mllikelihood)
        raise Exception((phi_0, lambd_0))

    return lambd, phi, mllikelihoods


def compute_weighted_means(parameters, organism_weights):
    sum_values = 0
    for organism in parameters.keys():
        parameter_value = parameters[organism]
        sum_values = sum_values + organism_weights[organism] * parameter_value

    return sum_values / sum(organism_weights.values())


def load_GO_annotations_mapper(args):
    # Load Gene Ontology -> Pfam annotation table
    logger.verbose("Loading Pfam to GO Biological Process annotations mappings")
    pfam_families = file_utilities.load_json(args["pfam_families_filename"])
    GO_to_pfam_mapping = file_utilities.load_json(args["go_to_pfam_filename"])
    logger.verbose(
        "Loading of Pfam to GO Biological Process annotations mappings done")

    return pfam_families, GO_to_pfam_mapping
