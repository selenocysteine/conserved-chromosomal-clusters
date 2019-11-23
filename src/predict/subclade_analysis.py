import os
import tables
import ete3
from utilities import file_utilities
from utilities.mymultiprocess import MyMultiProcess
import \
    predict.conserved_clustering_probabilities as conserved_clustering_probabilities
import predict.manage_probability_files as manage_probability_files
import multiprocessing

lock = multiprocessing.Lock()

import logging

logger = logging.getLogger(__name__)


class Subclade(tables.IsDescription):
    organism_id = tables.StringCol(500)
    subclade_id = tables.UInt16Col()


def partition_tree(probabilities_files_list,
                   organisms_uniform_probabilities,
                   n_pairs,
                   tree,
                   args):
    # Creates a directory for the output final probabilities file
    logger.info(
        "Started computation of subclade-specific "
        "clustering probabilities of Pfam domains")

    args["verbose2"] = args["verbose"]
    args["verbose"] = False

    if not os.path.exists(args["subclades_composition_filename"]) \
            or args["force"] is True:
        create_subclades_file(args["subclades_composition_filename"])
        parsed_subclades = set()
    else:
        parsed_subclades = load_parsed_subclades(
            args["subclades_composition_filename"])

    logger.info("Partitioning phylogenetic tree %s", args["tree_file"])

    subclades_IDs, subclades = extract_subclades(tree, args)
    n_subclades = len(subclades_IDs)
    n_leaves = len(tree.get_leaf_names())
    logger.info(
        "%s subclades (%s are leaves) retrieved from phylogenetic tree %s",
        n_subclades, n_leaves, args["tree_file"])

    # Loads all the individual estimated probabilities file names
    probabilities_files_list_dict = {}
    for prob_file in probabilities_files_list:
        organism_id = \
            file_utilities. \
                get_organism_from_file(prob_file,
                                       args["probabilities_filename_suffix"])
        try:
            probabilities_files_list_dict[organism_id].append(prob_file)
        except KeyError:
            probabilities_files_list_dict[organism_id] = [prob_file]

    if args["force"] is False:
        n_old_subclades = n_subclades
        subclades_IDs = sorted(
            list(set(subclades_IDs).difference(parsed_subclades)))
        if len(subclades_IDs) < n_old_subclades:
            already_done = n_old_subclades - len(subclades_IDs)
            logger.info(
                "%s subclades have already been processed and will be skipped",
                already_done)

    multi = MyMultiProcess(args["threads"],
                           analyse_subclade,
                           subclades_IDs,
                           [subclades,
                            probabilities_files_list_dict,
                            n_subclades,
                            n_pairs,
                            organisms_uniform_probabilities,
                            args],
                           maxtasksperchild=10)
    multi.run()

    if args["verbose2"] is True:
        logger.setLevel(logging.VERBOSE)
        args["verbose"] = True
        del args["verbose2"]

    logger.info("Subclades composition details saved in file %s",
                args["subclades_composition_filename"])
    logger.info(
        "Whole tree conserved clustering probabilities saved in file %s",
        args["final_probabilities_filename"])

    return 0


def analyse_subclade(subclade_id,
                     subclades,
                     probabilities_files_list_dict,
                     n_subclades,
                     n_pairs,
                     organisms_uniform_probabilities,
                     args):
    subclade_tree = subclades[subclade_id]
    subclade_organisms = sorted(subclade_tree.get_leaf_names())

    if len(subclade_organisms) == 1:
        subclade_type = "subclade (leaf)"
    else:
        if subclade_id == 0:
            subclade_type = "subclade (root node)"
        else:
            subclade_type = "subclade"

    logging.info("Analysing %s %s of %s, %s total organisms",
                 subclade_type, subclade_id, n_subclades,
                 len(subclade_organisms))

    if not subclade_id == 0:
        subclade_final_probabilities_filename = \
            os.path.join(args["subclade_clustering_dir"],
                         args["final_probabilities_pattern"] + "_" +
                         args["subclades_final_probabilities_pattern"]
                         .format(subclade_id))
    else:
        subclade_final_probabilities_filename = \
            args["final_probabilities_filename"]

    subclade_probabilities_files_list = []
    for organism in subclade_organisms:
        try:
            subclade_probabilities_files_list = \
                subclade_probabilities_files_list \
                + probabilities_files_list_dict[organism]
        except KeyError as e:
            logging.exception("Skipped %s %s as individual "
                              "probability file for %s is missing.",
                              subclade_type, subclade_id, e)
            raise KeyError("Skipped {} {} as individual "
                           "probability file for {} is missing."
                           .format(subclade_type, subclade_id, e))

    conserved_clustering_probabilities \
        .compute_conserved_probabilities(subclade_probabilities_files_list,
                                         organisms_uniform_probabilities,
                                         n_pairs,
                                         subclade_tree,
                                         subclade_id,
                                         subclade_final_probabilities_filename,
                                         args)

    save_subclade_composition(subclade_organisms, subclade_id,
                              args["subclades_composition_filename"])

    logger.info("%s %s of %s analysed",
                subclade_type.capitalize(), subclade_id, n_subclades)

    return 0


def generate_labelled_tree(tree, labelled_tree_file, **kwargs):
    if "subclades_IDs" in kwargs and "subclades" in kwargs:
        extract_subclades = True
        subclades_IDs = kwargs["subclades_IDs"]
        subclades = kwargs["subclades"]
    else:
        extract_subclades = False

    if extract_subclades is True:
        subclades_IDs.append(0)
        subclades[0] = tree
    subclade_ID = 1
    for node in reversed(list(tree.traverse("levelorder"))):
        if node.is_root():
            node.name = 0
        elif not node.is_leaf():
            node.name = subclade_ID
        else:
            node.name = str(subclade_ID) + "__" + node.name

        if extract_subclades is True and not node.is_root():
            subclades_IDs.append(subclade_ID)
            subclades[subclade_ID] = node

        subclade_ID = subclade_ID + 1

    tree.write(outfile=labelled_tree_file,
               format=1,
               format_root_node=True)
    for leaf in tree.get_leaves():
        leaf.name = leaf.name.split("__")[1]

    return tree


def load_labelled_tree(labelled_tree_file, **kwargs):
    if "subclades_IDs" in kwargs and "subclades" in kwargs:
        extract_subclades = True
        subclades_IDs = kwargs["subclades_IDs"]
        subclades = kwargs["subclades"]
    else:
        extract_subclades = False

    tree = ete3.Tree(labelled_tree_file,
                     format=1)
    if extract_subclades is True:
        subclades_IDs.append(0)
        subclades[0] = tree
    for node in reversed(list(tree.traverse("levelorder"))):
        if node.is_root():
            continue
        elif node.is_leaf():
            node_name = node.name.split("__")
            node.name = node_name[1]
            subclade_ID = int(node_name[0])
        else:
            subclade_ID = int(node.name)

        if extract_subclades is True:
            subclades_IDs.append(subclade_ID)
            subclades[subclade_ID] = node

    return tree


def extract_subclades(tree,
                      args):
    labelled_tree_file = args["labelled_tree_file"]

    subclades_IDs = []
    subclades = dict()
    if not os.path.exists(labelled_tree_file) or args["force"] is True:
        generate_labelled_tree(tree, labelled_tree_file,
                               subclades_IDs=subclades_IDs, subclades=subclades)
    else:
        load_labelled_tree(labelled_tree_file, subclades_IDs=subclades_IDs,
                           subclades=subclades)

    return subclades_IDs, subclades


def create_subclades_file(subclade_filename):
    subclade_h5file = tables.open_file(subclade_filename,
                                       mode="w", title="Subclades")
    group = subclade_h5file.create_group("/", 'subclades', 'Subclades')
    table = subclade_h5file.create_table(group, "subclades",
                                         Subclade,
                                         "Subclades composition")

    table.flush()
    subclade_h5file.close()

    return 0


def save_subclade_composition(subclade_organisms,
                              subclade_id,
                              subclade_filename):
    with lock:
        subclade_h5file = manage_probability_files.synchronized_open_file(
            subclade_filename, mode="a")
        table = subclade_h5file.root.subclades.subclades
        subclade_or = table.row
        for organism in subclade_organisms:
            subclade_or['subclade_id'] = subclade_id
            subclade_or['organism_id'] = organism
            subclade_or.append()
        table.flush()
        manage_probability_files.synchronized_close_file(subclade_h5file)

    return 0


def load_parsed_subclades(subclade_filename):
    parsed_subclades = set()

    subclades_h5file = tables.open_file(subclade_filename, mode="r")
    table = subclades_h5file.root.subclades.subclades

    for subclade_or in table.iterrows():
        parsed_subclades.add(subclade_or['subclade_id'])

    subclades_h5file.close()

    return parsed_subclades
