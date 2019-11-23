"""These functions just take care of the import-export of the computed probabilities,
in a thread-safe manner and using pytables (that can be indexed, and are easier to query afterwards).
"""
import numpy as np
import multiprocessing
import tables
from math import floor, log10
lock = multiprocessing.Lock()
from utilities import file_utilities

class FinalProbability(tables.IsDescription):
    pair_id = tables.Int64Col()
    occurrences = tables.UInt32Col()
    n_genes1 = tables.UInt32Col()
    n_genes2 = tables.UInt32Col()
    fusions = tables.UInt16Col()
    n_organisms = tables.UInt16Col()
    prob = tables.Int16Col()


def synchronized_open_file(*args, **kwargs):
    with lock:
        return tables.open_file(*args, **kwargs)


def synchronized_close_file(self, *args, **kwargs):
    with lock:
        return self.close(*args, **kwargs)


def is_table_empty(probabilities_filename):
    pair_h5file = tables.open_file(probabilities_filename, mode="r")
    table = pair_h5file.root.probabilities.pairs
    nrows = table.nrows
    pair_h5file.close()

    if nrows == 0:
        return True
    else:
        return False


def convert_prob_to_int(prob, factor=3):
    if prob == 0:
        return int(0)
    elif prob == 1:
        return int(10**factor)
    else:
        new_prob = int(floor(log10(prob)))
        if new_prob >= int(-factor):
            return int(prob * 10**factor)
        else:
            return new_prob


def convert_intprob_to_float(int_prob, factor=3):
    if int_prob == 0:
        return 0
    elif int_prob > 0:
        return float(int_prob) / 10**factor
    else:
        return 10**int_prob 



def save_uniform_probabilities(uniform_probabilities,
                               uniform_probability_filename):
    with lock:
        file_utilities.save_json(dict(uniform_probabilities), uniform_probability_filename)

    return 0


def load_uniform_probabilities(uniform_probability_filename):
    return file_utilities.load_json(uniform_probability_filename)


def save_conserved_probabilities(domain_pairs_data,
                                 final_probabilities_filename):

    final_prob_h5file = tables.open_file(final_probabilities_filename, mode="w",
                                         title="Conserved clustering probabilities")
    group = final_prob_h5file.create_group("/", 'probabilities', 'Probabilities')
    table = final_prob_h5file.create_table(group, "pairs",
                                     FinalProbability,
                                     "Pair probability")

    # To free up some memory, we first filter the pairs that shouldn't be saved
    pair = table.row
    indexes = np.where(domain_pairs_data['occurrences'] != 0)[0]
    
    for i in indexes:
        pair['pair_id'] = i
        pair['prob'] = convert_prob_to_int(domain_pairs_data['prob'][i])
        for key in ['occurrences', 'fusions', 'n_genes1',
                    'n_genes2', 'n_organisms']:
            pair[key] = domain_pairs_data[key][i]
        pair.append()

    table.flush()
    final_prob_h5file.close()

    return 0


