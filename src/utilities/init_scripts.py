"""Loads initial configuration and creates output directories and files."""
import os
from utilities import file_utilities
import configparser
import sys
import logging
logger = logging.getLogger(__name__)


def create_dir(path, new_dir_name):
    new_dir = os.path.abspath(os.path.join(path, new_dir_name))
    if not os.path.isdir(new_dir):
        os.makedirs(new_dir)
        logger.verbose("Directory %s created", new_dir)
    else:
        logger.verbose("Directory %s already present", new_dir)
    return new_dir


def create_analysis_output_directory(out_dir, data_dir):
    out_dir_name = os.path.basename(out_dir.rstrip("/"))
    data_dir = data_dir.rstrip("/")
    if os.path.isabs(out_dir):
        out_dir_path = os.path.dirname(out_dir)
    elif out_dir_name == out_dir:
        out_dir_path = os.path.abspath(os.path.dirname(data_dir))
    else:
        out_dir_path = os.path.join(os.path.dirname(os.getcwd()),
                                    os.path.dirname(out_dir))
    out_dir = create_dir(path=out_dir_path,
                         new_dir_name=out_dir_name)
    return out_dir


def check_input_files(args):
    args["gbk_files_list"] = file_utilities.get_file_list(args["data_dir"],
                                                          ".gbk")

    if len(args["gbk_files_list"]) == 0:
        logger.error("No gbk files in %s, are you sure you are using the right extension (%s) "
                     "and the right input folder?",
                     args["data_dir"], ".gbk")
        raise ValueError("No gbk files provided")

    args["pfam_annotations_file_list"] = \
        file_utilities.get_file_list(args["data_dir"],
                                     args["pfam_annotations_filename_suffix"])
    if len(args["pfam_annotations_file_list"]) == 0:
        logger.error("No Pfam domains annotations files in %s, "
                     "are you sure you are using the right extension (%s) "
                     "and the right input folder?",
                     args["data_dir"],
                     args["pfam_annotations_filename_suffix"])
        raise ValueError("No Pfam domains annotations files provided")

    if len(args["pfam_annotations_file_list"]) < len(args["gbk_files_list"]):
        logger.warning("%s Pfam domains annotations files are missing, only %s genomes will be analysed",
                       len(args["gbk_files_list"]) - len(args["pfam_annotation_files_list"]),
                       len(args["pfam_annotation_files_list"]))

    return 0


def create_plotting_directories(args):
    args["out_dir"] = create_analysis_output_directory(args["out_dir"], args["data_dir"])
    args["paper_figures_dir"] = create_dir(args["out_dir"], args["paper_figures_dir"])

    logger.info("Output directory for plots is: %s", args["paper_figures_dir"])


def create_output_directories(args):
    args["out_dir"] = create_analysis_output_directory(args["out_dir"], args["data_dir"])

    # General input files
    args["pfam_families_filename"] = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                                                  args["pfam_families_filename"]))
    args["pfam_details_filename"] = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                                                  args["pfam_details_filename"]))
    args["go_to_pfam_filename"] = \
        os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                     args["go_to_pfam_filename"]))
    if args["weight"] or args["subclades"] is True:
        if not os.path.isabs(args["tree_file"]):
            if args["tree_file"] == os.path.basename(args["tree_file"]):
                args["tree_file"] = os.path.abspath(os.path.join(args["data_dir"], args["tree_file"]))
            else:
                args["tree_file"] = os.path.abspath(os.path.join(os.path.dirname(os.getcwd()),
                                                                 args["tree_file"]))

    # Coordinate parsing files
    args["parsed_coordinates_out_dir"] = \
        create_dir(args["out_dir"],
                   args["parsed_coordinates_filename_suffix"])
    args["genome_sizes_filename"] = os.path.join(args["parsed_coordinates_out_dir"],
                                                 "genome_sizes.json.gz")


    # General fitting files
    init_fitting_dirs(args)


    # Prediction files
    init_predict_dirs(args)


    logger.info("Output directories created inside results directory: %s",
                args["out_dir"])


def init_fitting_dirs(args):
    args["fitting_dir"] = os.path.join(args["out_dir"], os.path.basename(args["fitting_dir"]))
    args["general_parameters_filename"] = os.path.join(args["fitting_dir"],
                                                       args["general_parameters_filename"])
    args["organisms_parameters_filename"] = os.path.join(args["fitting_dir"],
                                                         args["organisms_parameters_filename"])
    args["bootstrap_dir"] = os.path.join(args["fitting_dir"],
                                         os.path.basename(args["bootstrap_dir"]))
    args["sensitivity_dir"] = os.path.join(args["fitting_dir"],
                                           os.path.basename(args["sensitivity_dir"]))
    args["single_distances_dir"] = os.path.join(args["fitting_dir"],
                                                args["pairwise_distances_filename_suffix"])
    args["plots_dir"] = os.path.join(args["fitting_dir"], args["plots_dir"])
    args["em_plots_dir"] = os.path.join(args["plots_dir"], args["em_plots_dir"])

    if args["fit_model"] is True or args["bootstrap"] is True or args["sensitivity"] is True:
        create_dir(args["out_dir"], os.path.basename(args["fitting_dir"]))
    if args["plots"] is True:
        create_dir(args["fitting_dir"], os.path.basename(args["plots_dir"]))
        create_dir(args["plots_dir"], os.path.basename(args["em_plots_dir"]))
    if args["bootstrap"] is True:
        create_dir(args["fitting_dir"], os.path.basename(args["bootstrap_dir"]))
    if args["sensitivity"] is True:
        create_dir(args["fitting_dir"], os.path.basename(args["sensitivity_dir"]))
    if args["save_pairs"] is True:
        create_dir(args["fitting_dir"],
                   os.path.basename(args["pairwise_distances_filename_suffix"]))


def init_predict_dirs(args):
    args["final_probabilities_dir"] = os.path.join(args["out_dir"],
                                                   args["final_probabilities_dir"])

    args["individual_probabilities_dir"] = os.path.join(args["out_dir"],
                                                        args["individual_probabilities_dir"])

    if args["predict"] is True or args["permutations"] is True:
        create_dir(args["out_dir"], os.path.basename(args["final_probabilities_dir"]))
        create_dir(args["out_dir"], os.path.basename(args["individual_probabilities_dir"]))

    args["weights_file"] = os.path.join(args["final_probabilities_dir"],
                                        args["weights_file"])
    args["final_probabilities_filename"] = os.path.join(args["final_probabilities_dir"],
                     os.path.basename(args["final_probabilities_filename"]))

    args["uniform_probabilities_filename"] = os.path.join(args["individual_probabilities_dir"],
                                                          args["uniform_probabilities_filename"])

    # Pair indexes files
    if args["predict"] is True or args["permutations"] is True:
        create_dir(args["out_dir"], os.path.basename(args["pairs_indexes_dir"]))

    args["pairs_indexes_dir"] = os.path.join(args["out_dir"], args["pairs_indexes_dir"])
    args["pairs_indexes_filename"] = os.path.join(args["pairs_indexes_dir"],
                                                  args["pairs_indexes_filename"])
    args["reversed_indexes_filename"] = args["pairs_indexes_filename"].replace(".json.gz",
                                                                               "_reversed.json.gz")

    # Subclade analysis files
    if args["subclades"] is True:
        create_dir(args["final_probabilities_dir"], os.path.basename(args["subclade_clustering_dir"]))

    args["subclade_clustering_dir"] = os.path.join(args["final_probabilities_dir"],
                                                 os.path.basename(args["subclade_clustering_dir"]))
    labelled_tree_file = str.replace(os.path.basename(args["tree_file"]),
                                         os.path.splitext(args["tree_file"])[1],
                                         args["labelled_tree_file_suffix"])
    args["labelled_tree_file"] = os.path.join(args["subclade_clustering_dir"],
                                                  os.path.basename(labelled_tree_file))
    args["subclades_composition_filename"] = os.path.join(args["subclade_clustering_dir"],
                                                          os.path.basename(args["subclades_composition_filename"]))


    if args["permutations"] is True:
        create_dir(args["out_dir"], args["permutations_dir"])
    args["permutations_dir"] = \
        os.path.join(args["out_dir"], args["permutations_dir"])
    args['permutations_file'] = os.path.join(args["permutations_dir"],
                                             args['permutations_file'])
    args['distances_file'] = os.path.join(args["permutations_dir"],
                                             args['distances_file'])


def create_arguments_dict(args, config_file):
    if not os.path.exists(config_file):
        logger.error("Config file %s not found", config_file)
        sys.exit()

    config = configparser.ConfigParser()
    config.read([config_file])
    arguments_dict = args.__dict__

    for group in config._sections:
        for key in config._sections[group]:
            try:
                value = arguments_dict[key]
            except KeyError:
                value = None
            if value is None:
                value = config._sections[group][key]
                if group == "Floats":
                    value = float(value)
                if group == "Integers":
                    value = int(value)
                if group == "Strings":
                    value = value
                if group == "Booleans":
                    value = bool(int(value))
                arguments_dict[key] = value

    arguments_dict["matplotlib_style_file"] = \
        os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(config_file)),
                                     arguments_dict["matplotlib_style_file"]))

    return arguments_dict
