"""Gerstein-Sonnhammer-Chothia (GSC) weighing algorithm.

This scripts implements a Python2 / 3 version of the GSC algorithm,
a weighing scheme for the leaves of a (phylogenetic) tree which take into
account both the length of the branches connecting each leaf to the root,
and the total number of leaves found in each region of the tree
(upweighing leaves coming from scarsely populated regions).
The script takes advantage of the existing implementations of tree structures
and level-order traversal from the ete3 Python library.

The complete algorithm is described in the Supplementary Materials of
Gerstein M, Sonnhammer EL, Chothia C, Volume Changes in Protein Evolution (1994),
J Mol Biol., doi:10.1016/0022-2836(94)90012-4
"""
from __future__ import print_function
import ete3
import sys
import numpy as np


def test_scoring(GSC_function):
    """Test a function implementing the GSC scoring algorithm (GSCfunction) on the tree used as an example in the paper.
    Note: the weights are assessed as before the normalisation to average=1 step described in the paper."""
    import ete3

    function_name = GSC_function.__name__

    if not callable(GSC_function):
        raise ValueError("The provided input function {} is not a function.".format(function_name))

    newick_tree = "(D:80,(C:50,(A:20,B:20)two:30)three:30);"
    tree = ete3.Tree(newick_tree, format=1)

    scores = GSC_function(tree)

    correct_scores = {'A': 43.75, 'C': 62.5, 'B': 43.75, 'D': 80.0}

    correct = 0

    for key in sorted(scores):
        print("Leaf {}: Correct score: {} - Score obtained by function {}: {}"
              .format(key, correct_scores[key], function_name, scores[key]))
        if correct_scores[key] == scores[key]:
            correct += 1

    print("{}% scores were correctly predicted ({} out of {})"
          .format(float(correct)/len(scores) * 100, correct, len(scores)))

    if correct == len(scores):
        return True
    else:
        return False


def GSC(t):
    """Take an ete3 tree 't' as input, and returns a dictionary leaf ID : GSC weight with the GSC weights of each leaf of the tree."""
    scores = {}

    if not isinstance(t, ete3.Tree):
        raise ValueError("The provided input is not an object of class ete3.Tree().")

    # Traverse the tree level by level, starting from the most external levels
    # up to the root
    for node in reversed(list(t.traverse("levelorder"))):

        children = node.get_children()

        # If the node 'node' is a leaf, its score is initialised to zero
        if len(children) == 0:
            scores[node.name] = 0

        else:

            for child in children:
                # If the *direct* child of node is a leaf, set its score to its
                # distance from node, without weighing
                if len(child.get_children()) == 0:
                    if child.dist != 0:
                        scores[child.name] = child.dist
                    else:
                        # In case there is a polytomy, and the distance between node and
                        # the leaf child is zero, the score is switched from zero to them
                        # smallest system float (epsilon), to prevent the weight of this leaf
                        # from being zero throughout the whole tree
                        scores[child.name] = sys.float_info.min

                # If the child of node is not a leaf, update the scores of the downstream
                # leaves by weighing them
                else:
                    # Find all the leaves downstream of child
                    species = child.get_leaf_names()

                    # Get the cumulative score of all the leaves downstream of child
                    sum_of_weights = float(sum([scores[name] for name in species]))

                    # For every leaf dowstream of child, update the weight as
                    # their previous score plus the distance between node and child
                    # multiplied by the score of each leaf over the cumulative score
                    # of all the leaves dowstream of child
                    for name in species:
                        # To correct for polytomies and branches of length 0:
                        if child.dist == 0:
                            scores[name] = scores[name] + sys.float_info.min * float(scores[name])/sum_of_weights

                        else:
                            scores[name] = scores[name] + child.dist * float(scores[name])/sum_of_weights

    return scores


def GSC_one_normalised(t):
    """Compute the GSC weights of an ete3 tree 't' and normalise them so that the average weight is 1, as described in the paper."""

    if not isinstance(t, ete3.Tree):
        raise ValueError("The provided input is not an object of class ete3.Tree().")

    scores = GSC(t)

    total_score = float(sum(scores.values()))
    n_scores = len(scores)

    # Normalises the scores
    return {key: scores[key]/total_score * n_scores for key in scores}


def GSC_normalised(t):
    """Compute the GSC weights of an ete3 tree 't' and weights them.
    """

    if not isinstance(t, ete3.Tree):
        raise ValueError("The provided input is not an "
                         "object of class ete3.Tree().")

    scores = GSC(t)

    values = np.fromiter(scores.values(), dtype=float)
    keys = list(scores.keys())

    min_score = np.min(values)
    max_score = np.max(values)
    range_scores = max_score - min_score
    max_new_score = 1
    min_new_score = float(1)/len(scores)
    new_range_scores = float(max_new_score - min_new_score)
    if range_scores < 10**-6:
        range_scores = 0
    if len(scores) == 1:
        values = np.array([1])
    elif range_scores == 0:
        values = np.repeat(1, len(values))
    else:
        values = new_range_scores * \
                 ((values - min_score) / range_scores) + min_new_score

    for i, key in enumerate(keys):
        scores[key] = values[i]

    return scores


def save_weights(scores, weights_filename, subclade, mode="a"):
    if mode == "w":
        with open(weights_filename, mode) as out_file:
            out_file.write("Organism\tWeight\tSubclade_id\n")
    else:
        with open(weights_filename, mode) as out_file:
            for item in scores:
                out_file.write("{}\t{:.3g}\t{}\n"
                               .format(item, scores[item], subclade))

