import os
import logging
import numpy as np
import dendropy

logger = logging.getLogger(__name__)
import scipy.sparse


def filter_data_table(data_table, elements):
    n_organisms = len(elements)
    dt = np.zeros((n_organisms, n_organisms),
                  dtype=np.float32)
    for index1 in range(0, n_organisms - 1):
        t1 = get_code(elements[index1])
        for index2 in range(index1 + 1, n_organisms):
            t2 = get_code(elements[index2])
            distance = data_table[t1, t2]
            if distance < 0:
                distance = 0
            dt[index1, index2] = distance
    dt = scipy.sparse.csr_matrix(dt)
    return dt


def get_code(text):
    code = text.replace("_", " ").split(" ")[1].split(".")[0]
    return code


def normalise_distance_matrix(matrix_file, normalised_matrix_file):
    logger.info("Normalising distance matrix")
    data = np.genfromtxt(matrix_file, delimiter=',', dtype=None)
    data2 = data[1:, 1:].astype(np.float)
    data2[data2 < 0] = 0
    min_val = np.min(data2)
    max_val = np.max(data2)
    data[1:, 1:] = (1 - 0.0001) * (data2 - min_val) / (
            max_val - min_val) + 0.0001
    data = data.astype(str)
    with open(normalised_matrix_file, 'w') as out_file:
        np.savetxt(out_file, data[0, :].reshape(1, data.shape[0]),
                   delimiter=",", fmt="%s")
    with open(normalised_matrix_file, 'a') as out_file:
        np.savetxt(out_file, data[1:, :], delimiter=",",
                   encoding=str, fmt='%s')
    logger.info("Distance matrix normalised")


def compute_distance_matrix(tree, args):
    logger.info("Computing distance matrix")
    matrix_file = args['distances_file']
    normalised_matrix_file = os.path.join(os.path.dirname(matrix_file),
                                          "normalised_matrix.csv")
    tree2 = dendropy.Tree.get(data=tree.write(format=1,
                                              format_root_node=True),
                              schema="newick")
    if not os.path.exists(matrix_file) and args['force'] is False:
        distances = tree2.phylogenetic_distance_matrix().as_data_table()
        distances.write_csv(out=matrix_file)
    if not os.path.exists(normalised_matrix_file) and args['force'] is False:
        normalise_distance_matrix(matrix_file, normalised_matrix_file)

    distances = \
        dendropy.utility.container \
            .DataTable().from_csv(normalised_matrix_file,
                                  default_data_type=float,
                                  label_transform_fn=get_code)
    logger.info("Distance matrix computed")

    return distances
