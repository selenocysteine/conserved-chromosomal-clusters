""" This scripts calculates whenever a certain pair disappeared at a certain
node in the tree, and deletes it from all the subclades above that node.
"""
import os
import numpy as np
import tables
import ete3
import re
import time
from utilities import file_utilities
from predict.manage_probability_files import FinalProbability
from utilities.mymultiprocess import MyMultiProcess
import multiprocessing

manager = multiprocessing.Manager()

import logging

logger = logging.getLogger(__name__)


def parse_and_clean_tree(n_pairs, args):
    subclades_final_probabilities_files = [args[
                                               "final_probabilities_filename"]] + \
                                          file_utilities.get_file_list(
                                              args["subclade_clustering_dir"],
                                              args[
                                                  "final_probabilities_pattern"],
                                              verbose=True)


    processed_files = load_processed_files_list(len(subclades_final_probabilities_files),
                                      args)
    already_processed = np.sum(processed_files)
    if (len(subclades_final_probabilities_files) - already_processed) > 0:
        logger.info(
            "%s subclade-specific conserved clustering "
            "probabilities files have to be parsed",
            len(subclades_final_probabilities_files) - already_processed)

        processed_files = manager.list(processed_files)
        tree = ete3.Tree(args["labelled_tree_file"], format=1)
        subclades_to_process = [node.name for node in tree.traverse('levelorder')]

        multi = MyMultiProcess(threads=args["threads"],
                               target=process_child_node,
                               input=subclades_to_process,
                               args=[tree,
                                     processed_files,
                                     n_pairs,
                                     args],
                               destroy=True)
        multi.run()


def process_child_node(node_name, tree, processed_files, n_pairs, args):
    node = tree.search_nodes(name=node_name)[0]
    node_id = get_node_id(node)
    node_final_probabilities_file = get_subclade_prob_filename(node_id,
                                                               args)

    if not processed_files[node_id] is True:
        if node_id == 0:
            parent_done = True
        else:
            parent_id = get_node_id(node.up)
            parent_done = processed_files[parent_id]

        while parent_done is False:
            time.sleep(1)
            parent_done = processed_files[parent_id]

        if node.is_leaf():
            logger.verbose("Subclade %s is a leaf, no cleaning up required",
                           node_name)
            os.rename(node_final_probabilities_file,
                      node_final_probabilities_file.
                      replace(".h5", "_parsed.h5"))
            processed_files[node_id] = True
        else:
            has_parents = False if node_id == 0 else True
            if node_id == 0 or get_node_id(node.up) == 0:
                to_keep = True
            else:
                to_keep = False
            processed_files[node_id] = \
                parse_node_final_probabilities(node,
                                               node_final_probabilities_file,
                                               n_pairs,
                                               has_parents,
                                               to_keep,
                                               args)

    return 0


def parse_node_final_probabilities(node, node_final_probabilities_file,
                                   n_pairs, has_parents, to_keep,
                                   args):
    logger.verbose(
        "Cleaning up subclade probability file %s to remove "
        "pairs that don't get this far "
        "in the phylogeny",
        node_final_probabilities_file)

    deleted_rows, old_rows = \
        update_node_probability_file(node, n_pairs,
                                     node_final_probabilities_file,
                                     has_parents, to_keep, args)
    logger.verbose(
        "Final probabilities file %s cleaned up, %s pairs deleted out of %s",
        node_final_probabilities_file,
        '{:,}'.format(deleted_rows),
        '{:,}'.format(old_rows))

    return True


def update_node_probability_file(node,
                                 n_pairs,
                                 node_final_probabilities_file,
                                 has_parents,
                                 to_keep,
                                 args):
    occurrence_matrix, prob_pairs = \
        examine_occurrences_children_and_parent_nodes(
            node, n_pairs, args)

    # We create a new renamed file and delete the old one in the end
    node_new_final_probabilities_file = node_final_probabilities_file.replace(
        ".h5", "_temp.h5")
    node_new_final_probabilities_h5file = tables.open_file(
        node_new_final_probabilities_file, mode="w")
    new_group = node_new_final_probabilities_h5file.create_group("/",
                                                                 'probabilities',
                                                                 'Probabilities')
    new_table = node_new_final_probabilities_h5file.create_table(new_group,
                                                                 "pairs",
                                                                 FinalProbability,
                                                                 "Pair probability")
    new_rows = new_table.row

    old_probabilities_h5file = tables.open_file(node_final_probabilities_file,
                                                mode="r")
    old_probabilities_table = old_probabilities_h5file.root.probabilities.pairs
    deleted_rows = 0
    old_rows = old_probabilities_table.nrows

    for old_pair in old_probabilities_table.iterrows():
        pair_occurrences_in_children = occurrence_matrix[0, old_pair['pair_id']]
        if pair_occurrences_in_children != old_pair['occurrences']:
            copy_row(new_rows, old_pair)
        else:
            if has_parents is True:
                if occurrence_matrix[1, old_pair['pair_id']] > 0:
                    copy_row(new_rows, old_pair,
                             prob_pairs[old_pair['pair_id']])
                else:
                    deleted_rows = deleted_rows + 1
            else:
                deleted_rows = deleted_rows + 1

    old_probabilities_h5file.close()
    new_table.flush()
    node_new_final_probabilities_h5file.close()
    os.rename(node_new_final_probabilities_file,
              node_new_final_probabilities_file.replace("_temp", "_parsed"))
    if to_keep is False:
        os.remove(node_final_probabilities_file)

    return deleted_rows, old_rows


def examine_occurrences_children_and_parent_nodes(node, n_pairs, args):
    # If the node is not root, we add a row for the parent's data (first row)
    if not node.is_root():
        n_rows = 2
    else:
        n_rows = 1
    occ_matrix = np.zeros((n_rows, n_pairs), dtype=int)
    prob_pairs = np.zeros(n_pairs, dtype=int)

    # If the node has parents, we add its parent's info in the 
    # first line of the occurrence matrix
    # Examine the occurrences in all children
    for child in node.children:
        examine_node_occurrences(child, occ_matrix[0,], prob_pairs,
                                 args,
                                 has_parent=False)
    if not node.is_root():
        examine_node_occurrences(node.up, occ_matrix[1,], prob_pairs,
                                 args,
                                 has_parent=True)

    return occ_matrix, prob_pairs


def examine_node_occurrences(node, occ_row, prob_pairs,
                             args,
                             has_parent=False):
    subclade_filename = get_subclade_prob_filename(get_node_id(node),
                                                   args,
                                                   has_parent)
    subclade_h5file = tables.open_file(subclade_filename, mode="r")
    table = subclade_h5file.root.probabilities.pairs
    for pair in table.iterrows():
        if occ_row[pair['pair_id']] < pair['occurrences']:
            occ_row[pair['pair_id']] = pair['occurrences']
            if has_parent is True:
                prob_pairs[pair['pair_id']] = pair['prob']
    subclade_h5file.close()


def copy_row(new_row, old_pair, new_prob=None):
    for name in ['pair_id', 'occurrences', 'fusions',
                 'n_genes1', 'n_genes2', 'n_organisms']:
        new_row[name] = old_pair[name]
    if new_prob is None:
        new_row['prob'] = old_pair['prob']
    else:
        new_row['prob'] = new_prob

    new_row.append()


def get_node_id(node):
    if "__" in node.name:
        node_name = node.name.split("__")[0]
    elif node.name == "":
        node_name = 0
    else:
        node_name = node.name

    return int(node_name)


def get_subclade_prob_filename(node_id, args, parent=False):
    if node_id == 0:
        subclade_filename = args["final_probabilities_filename"]
    else:
        subclade_filename = os.path.join(args["subclade_clustering_dir"],
                                         "final_probabilities_subclade{}.h5"
                                         .format(node_id))

    if parent is True:
        subclade_filename = subclade_filename.replace(".h5", "_parsed.h5")

    return os.path.abspath(subclade_filename)


def get_subclade_id_from_filename(filename, args):
    subclade_pattern = args['subclades_final_probabilities_pattern'].replace(
        "{}.h5", "([0-9]+)_parsed.h5")
    subclade_id = re.search(r"{}".format(subclade_pattern), filename)

    if subclade_id is not None:
        subclade_id = int(subclade_id.groups()[0])
    else:
        subclade_id = 0
    return subclade_id


def load_processed_files_list(n_total_files, args):
    processed_files = [False] * n_total_files
    new_processed_files = file_utilities.get_file_list(
        args["final_probabilities_dir"],
        "_parsed.h5",
        verbose=False, error=False)
    for filename in new_processed_files:
        subclade_id = get_subclade_id_from_filename(filename, args)
        processed_files[subclade_id] = True

    return processed_files
