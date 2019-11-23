#!/usr/bin/python3
import logging, verboselogs
from utilities.logging_formatter import init_logger
verboselogs.install()
logger = init_logger()

import sys
import os
from genome_parsing import gbk_to_fast_coordinates
from fit import fit_clustering_model
from fit import bootstrap
from predict import predict_clustering_pairs as predict
from permutations import permutation_analysis as permute
from utilities import import_export_parameters
from utilities import init_scripts
from utilities.file_utilities import load_json


def main(args):
    config_file = "global_config.conf"
    args = init_scripts.create_arguments_dict(args, config_file)
    if args["verbose"] is True:
        logger.setLevel(logging.VERBOSE)
    else:
        logger.setLevel(logging.INFO)

    # Check if the input directory exists
    data_dir = os.path.abspath(args["data_dir"])
    if not os.path.exists(data_dir):
        logger.error("Input directory %s does not exist", data_dir)
        raise ValueError("Input directory {} does not exist".format(data_dir))
    else:
        args["data_dir"] = data_dir
        logger.info("Input directory: %s", data_dir)

    init_scripts.check_input_files(args)

    # Check if there is something to run (options)
    parameters = []
    for el in ["fit_model", "predict", "permutations", "bootstrap", "sensitivity"]:
        if args[el] is True:
            parameters.append(el)

    if len(parameters) == 0:
        logger.error("Select any among these options: {}"
                     .format(", --".join(parameters)))
        raise ValueError("Either one of --{} options have to be selected"
                         .format(", --".join(parameters)))
    else:
        logger.info("Clustering analysis will be run in the following mode(s): %s",
                    ", ".join(parameters))

    logger.info("Number of threads: %s", args["threads"])

    # Creates output directories
    init_scripts.create_output_directories(args)

    if args["weight"] is True or args["subclades"] is True:
        if not os.path.exists(args["tree_file"]):
            logger.error("Input phylogenetic tree file %s not found",
                         args["tree_file"])
            raise ValueError("Input phylogenetic tree file {} not found"
                             .format(args["tree_file"]))
        else:
            logger.info("Input phylogenetic tree file: %s", args["tree_file"])

    if args["weight"] is True:
        logger.info("Conserved clustering probabilities will be "
                    "weighted using phylogenetic distances")
    if args["subclades"] is True:
        logger.info("Conserved clustering probabilities will be "
                    "computed for each subclade")
    else:
        logger.info("Conserved clustering probabilities will be "
                    "computed just for the root node")

    # Check if the parsed coordinates files already exist, if not create them from the gbk files
    if args["fit_model"] or args["predict"] is True:
        args["genome_sizes"] = gbk_to_fast_coordinates.parse_coordinates(args)
    else:
        args["genome_sizes"] = load_json(args['genome_sizes_filename'])

    if args["fit_model"] is True:
        args["lambd"], args["phi"] = fit_clustering_model.fit_clustering_model(args)
    else:
        if os.path.exists(args["general_parameters_filename"]):
            args["lambd"], args["phi"] = import_export_parameters\
                .import_parameters(args["general_parameters_filename"])

    if args["bootstrap"] is True:
        bootstrap.bootstrap_parameters(args)

    if args["sensitivity"] is True:
        bootstrap.sensitivity_analysis(args)

    if args["predict"] is True:
        if args["fit_model"] is False:
                logger.info("Imported global mean clustering estimated "
                            "parameter values are lambda={:.3g}, phi={:.3g}"
                    .format(args["lambd"], args["phi"]))
        else:
            logger.info("Using default global mean clustering estimated "
                            "parameter values lambda={:.3g}, phi={:.3g}"
                .format(args["lambd"], args["phi"]))

        predict.predict_clustering_pairs(args)

    if args["permutations"] is True:
        if os.path.exists(args['subclades_dir']):
            permute.predict_permuted_clustering_pairs(args)
        else:
            logger.error("Permutations are accepted only for subclade-specific "
                         "clustering analysis.")


if __name__ == "__main__":
    script_path = os.path.dirname(os.path.realpath(__file__))
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("-t",
                        help="Number of threads to use (default=1).",
                        action="store", type=int, dest="threads")

    parser.add_argument("-i",
                        help="Input files directory.",
                        action="store", type=str, dest="data_dir",
                        default=os.path.abspath(os.path.join(script_path,
                                                             "../data")))

    parser.add_argument("-o",
                        help="Output files directory.",
                        action="store", type=str, dest="out_dir",
                        default="results")

    parser.add_argument("--fit",
                        help="Estimate new clustering distances "
                             "parameters using the organisms "
                             "in the input mapping_tables directory (optional)",
                        action="store_true", dest="fit_model",
                        default=False)

    parser.add_argument("--predict",
                        help="Predict clustering probabilities for "
                             "all the pairs of Pfam domains from "
                             "the organisms in the input mapping_tables "
                             "directory (default=False).",
                        action="store_true", dest="predict",
                        default=False)

    parser.add_argument("--weight",
                        help="While computing conserved clustering probabilities, "
                             "weight each species according"
                             "to its position on a phylogenetic tree (default=False)",
                        action="store_true", dest="weight",
                        default=False)

    parser.add_argument("--subclades",
                        help="Compute a set of subclade-specific clustering "
                             "probabilities instead of just a single probability "
                             "for the whole dataset (default=False)",
                        action="store_true", dest="subclades",
                        default=False)

    parser.add_argument("--permutations",
                        help="Compute random permutations of gene distances "
                             "to assess significances",
                        action="store_true", dest="permutations",
                        default=False)

    parser.add_argument("--nperm",
                        help="Number of permutations (default=1).",
                        action="store",
                        type=int, dest="nperm",
                        default=1)

    parser.add_argument("--plots",
                        help="Generate plots for the fitting step "
                             "(default=False).",
                        action="store_true", dest="plots",
                        default=False)

    parser.add_argument("--evalue",
                        help="E-value threshold for Pfam predict "
                                        "(default=0.001).",
                        action="store",
                        type=float, dest="evalue",
                        default=0.001)

    parser.add_argument("--phi",
                        help="Initial value for phi (default=0.001).",
                        action="store", type=float, dest="phi_0",
                        default=0.001)

    parser.add_argument("--lambda",
                        help="Initial value for the rate parameter "
                             "(default=uses a moment based estimator)",
                        action="store", type=float, dest="lambd_0",
                        default=None)

    parser.add_argument("--bootstrap",
                        help="Perform bootstrap stability analysis of the parameters "
                             "(default=False)",
                        action="store_true", dest="bootstrap",
                        default=False)

    parser.add_argument("--sensitivity",
                        help="Perform sensitivity analysis with mean parameters "
                             "(default=False)",
                        action="store_true", dest="sensitivity",
                        default=False)

    parser.add_argument("--tree",
                        help="Path to the file containing the phylogenetic "
                             "tree of the analysed species",
                        action="store", type=str,
                        dest="tree_file",
                        default="phylogenetic_tree.nw")

    parser.add_argument("--save_pairs",
                        help="Save organism-specific domain pair distances "
                             "values on disk (default=False)",
                        action="store_true", dest="save_dist",
                        default=False)

    parser.add_argument("--non_go",
                        help="During the fit, compute pair distances of non-GO "
                             "related domains (default=False)",
                        action="store_true", dest="non_go",
                        default=False)

    parser.add_argument("--pfam",
                        help="Filename suffix for the Pfam input files "
                             "(default=pfam_annotations)",
                        action="store", type=str, dest="pfam_annotations_filename_suffix",
                        default="_pfam_annotations")

    parser.add_argument("--force",
                        help="Completely rewrite eventual previously existing "
                             "input files (default=False)",
                        action="store_true", dest="force",
                        default=False)

    parser.add_argument("--verbose",
                        help="Print detailed comments (default=False)",
                        action="store_true", dest="verbose",
                        default=False)

    sys.exit(main(parser.parse_args()))
