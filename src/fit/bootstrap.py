import os
import numpy as np
import ete3
from utilities import import_export_parameters
from fit import plot_functions as plot_func
from fit import estimate_clustering_parameters
from fit import GO_gene_pairs_distances
from predict import GSC
from utilities.mymultiprocess import MyMultiProcess
from utilities import em_algorithm as emCC
from multiprocessing import Manager
import logging

manager = Manager()
logger = logging.getLogger(__name__)


def bootstrap_parameters(args):
    np.random.seed(args["random_seed"])

    logger.info(
        "Started bootstrapping estimated clustering parameters"
        " from %s assemblies",
        len(args['pfam_annotations_file_list']))

    # Bootstraps the means of phi and of the rate parameter lambda
    parameters_data = import_parameter_data(args)

    multi = MyMultiProcess(threads=args["threads"],
                           target=compute_bootstrap_means,
                           input=['lambda', 'phi'],
                           args=[parameters_data,
                                 args],
                           maxtasksperchild=100)

    multi.run()

    logger.info(
        "Bootstrapping of the estimated clustering parameters "
        "from %s organisms done",
        len(args['gbk_files_list']))

    if args['weight'] is True:
        logger.info(
            "Weighted mean parameter values: phi: {:.3g}, lambda: {:.3g}"
                .format(args['phi'],
                        args['lambd']))

    logger.info("Plots of bootstrap results saved in directory %s",
                args["bootstrap_dir"])

    return 0


def weighted_mean(data, weights):
    return np.sum(data * weights) / np.sum(weights)


def bootstrap_sampling(data_array, n_samples, args, weights=None):
    bootstrap_means = np.zeros(n_samples)
    n_organisms = len(data_array)

    for i in range(0, n_samples):
        random_indexes = np.random.randint(0, n_organisms, n_organisms)
        data_bootstrap = data_array[random_indexes]
        if args['weight'] is True:
            weights_bootstrap = weights[random_indexes]
            bootstrap_means[i] = weighted_mean(data_bootstrap,
                                               weights_bootstrap)
        else:
            bootstrap_means[i] = np.mean(data_bootstrap)

    return bootstrap_means


def import_parameter_data(args):
    lambdas, phis = \
        import_export_parameters.import_organisms_parameters(
            args['organisms_parameters_filename'])


    if args['weight'] is True:
        parameters_data = {'lambda': np.zeros(len(lambdas)),
                           'phi': np.zeros(len(lambdas))}
        if os.path.exists(args['weights_file']):
            root_weights = {}
            with open(args['weights_file'], "r") as inp_file:
                for line in inp_file:
                    line = line.strip().split()
                    if line[-1] == "root":
                        root_weights[line[0]] = float(line[-2])
        else:
            root_weights = GSC.GSC_normalised(ete3.Tree(args['tree_file']))
        logging.info("Computing weighted mean")
        all_weights = np.zeros(len(lambdas))

        for index, organism in enumerate(sorted(root_weights.keys())):
            parameters_data['lambda'][index] = lambdas[organism]
            parameters_data['phi'][index] = phis[organism]
            all_weights[index] = root_weights[organism]

        return parameters_data, all_weights
    else:
        parameters_data = {'lambda': np.array(lambdas),
                           'phi': np.array(phis)}

        return parameters_data


def compute_bootstrap_means(parameter,
                            parameter_data,
                            args):
    n_samples = 10 ** 6
    logging.info("Generating %s bootstrap samples for parameter %s",
                 n_samples, parameter)
    np.random.seed(seed=args['random_seed'])
    if args['weight'] is True:
        weights = parameter_data[1]
        parameter_data = parameter_data[0][parameter]
    else:
        parameter_data = parameter_data[parameter]
        weights = None

    # Means of bootstrapped values
    bootstrap_means = bootstrap_sampling(parameter_data,
                                         n_samples, args, weights)

    # True mean of the original mapping_tables
    plot_file_name = os.path.join(args['bootstrap_dir'],
                                  parameter.replace(" ",
                                                    "_")
                                  + "_estimated.pdf")

    # Saves a histogram of the bootstrapped means of this parameter
    plot_func.plot_histogram_bootstrap(bootstrap_means,
                                       parameter,
                                       args,
                                       plot_file_name)

    bootstrapped_values = [np.mean(bootstrap_means),
                                      np.std(bootstrap_means),
                                      np.percentile(bootstrap_means, 2.5),
                                      np.percentile(bootstrap_means, 97.5)]
    logger.info("Parameter {} bootstrap {} samples. "
                "Mean = {:.3g}, "
                "standard deviation = {:.3g}, "
                "90% CI = ({:.3g}, {:.3g})".format(parameter,
                                                 n_samples,
                                                 bootstrapped_values[0],
                                                 bootstrapped_values[1],
                                                 bootstrapped_values[2],
                                                 bootstrapped_values[3]
                                                 ))

    return 0


def sensitivity_analysis(args):
    ''''This function computes all zetas with organism-specific parameter values
    and compares them with the mean zetas used to summarise the parameters'''

    lambdas, phis = \
        import_export_parameters. \
            import_organisms_parameters(args['organisms_parameters_filename'])
    pfam_families, GO_to_pfam_mapping = \
        estimate_clustering_parameters.load_GO_annotations_mapper(args)
    pfam_annotations_files_list = \
        args["pfam_annotations_file_list"]

    tp_rates = manager.list()
    fp_rates = manager.list()

    multi = MyMultiProcess(threads=args["threads"],
                           target=sensitivity_analysis_single_organism,
                           input=pfam_annotations_files_list,
                           args=[GO_to_pfam_mapping,
                                 pfam_families,
                                 lambdas,
                                 phis,
                                 tp_rates, fp_rates,
                                 args])
    multi.run()

    if args['plots'] is True:
        plot_file_name = os.path.join(args["sensitivity_dir"],
                                      "FPR_TPR.pdf")
        plot_func.plot_pseudo_auc(tp_rates, fp_rates, plot_file_name, args)


def sensitivity_analysis_single_organism(pfam_annotations_file,
                                         GO_to_pfam_mapping,
                                         pfam_families,
                                         lambdas,
                                         phis,
                                         tp_rates,
                                         fp_rates,
                                         args):
    organism_name, genome_size, GO_pairs_distances = \
        GO_gene_pairs_distances.import_GO_pairs_distances(pfam_annotations_file,
                                                          GO_to_pfam_mapping,
                                                          pfam_families,
                                                          args)

    true_zetas = emCC.compute_zetas(GO_pairs_distances,
                                    float(genome_size) / 2,
                                    lambdas[organism_name],
                                    phis[organism_name])

    mean_parameters_zetas = emCC.compute_zetas(GO_pairs_distances,
                                               float(genome_size) / 2,
                                               args['lambd'],
                                               args['phi'])

    # pearson_correlation = pearsonr(true_zetas, mean_parameters_zetas)

    mask_mean = mean_parameters_zetas > 0.5
    mask_original = true_zetas > 0.5
    total_true = len(mask_original[mask_original == True])
    tp_rates.append(np.sum(mask_mean * mask_original) / total_true)
    fp_rates.append(np.sum(mask_mean * ~mask_original) / total_true)

    # if args['plots'] is True:
    #     plot_file_name = os.path.join(args["sensitivity_dir"],
    #                                   organism_name + "_zeta_comparisons_plot.pdf")
    #     plot_func.zeta_correlation_plots(true_zetas, mean_parameters_zetas,
    #                                      pearson_correlation, lambdas[organism_name], phis[organism_name],
    #                                      organism_name, plot_file_name)
