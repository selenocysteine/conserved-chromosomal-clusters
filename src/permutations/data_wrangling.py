import os
import logging
import tables
import numpy as np
from progress.bar import Bar
from math import log10
from collections import defaultdict
from predict.organism_clustering_probabilities import \
    compute_uniform_probability, compute_log_exponential_probability
from genome_parsing.load_genome_annotations import import_gene_coord_data
from utilities import file_utilities
from genome_parsing.load_genome_annotations import import_pfam_hits

logger = logging.getLogger(__name__)


def convert_prob_to_float(prob):
    if prob > 0:
        return prob / 1000
    elif prob == 0:
        return 0
    else:
        return 10 ** float(prob)


def load_permutations(perm_file):
    done = defaultdict(list)
    with open(perm_file, "r") as inp_file:
        inp_file.readline()
        for line in inp_file:
            line = line.strip().split()
            subclade = int(line[0])
            perm = int(line[2])
            done[subclade].append(perm)

    return done


def load_all_organisms_coordinates(args):
    logger.info("Loading domain coordinates for all organisms")
    bar = Bar('Loading organism_data',
              max=len(args['pfam_annotations_file_list']))
    all_gene_coordinates = {}
    pfam_families = file_utilities.load_json(args["pfam_families_filename"])

    for pfam_annotations_file in args['pfam_annotations_file_list']:
        organism_name = \
            file_utilities. \
                get_organism_from_file(pfam_annotations_file,
                                       args["pfam_annotations_filename_suffix"])
        genome_id, genome_size = args["genome_sizes"][organism_name][0]
        gene_coordinates, gene_to_chr_mapping = import_gene_coord_data(
            organism_name, args)
        gene_hits = import_pfam_hits(pfam_annotations_file,
                                     gene_to_chr_mapping,
                                     args['evalue'],
                                     pfam_families,
                                     genome_id,
                                     domain_centric=False)

        coords = np.zeros([len(gene_hits), 2], dtype=int)
        probs = np.zeros(len(gene_hits), dtype=float)
        total_prob = 0
        for key in gene_hits:
            total_prob += len(gene_hits[key])
        for index, gene in enumerate(gene_hits.keys()):
            coords[index] = np.array([gene_coordinates[gene][0],
                                      gene_coordinates[gene][1]])
            probs[index] = len(gene_hits[gene]) / total_prob

        all_gene_coordinates[organism_name] = \
            [compute_uniform_probability(float(genome_size) / 2),
             float(genome_size),
             coords,
             probs]
        bar.next()

    bar.finish()
    logger.info("Domain coordinates for all organisms loaded")
    return all_gene_coordinates


def count_occurrences(subclade_id, args):
    if not subclade_id == 0:
        subclade_final_probabilities_filename = \
            os.path.join(args["subclade_clustering_dir"],
                         args["final_probabilities_pattern"] + "_" +
                         args["subclades_final_probabilities_pattern"]
                         .format(subclade_id))
    else:
        subclade_final_probabilities_filename = \
            args["final_probabilities_filename"]

    subclade_final_probabilities_filename = \
        subclade_final_probabilities_filename.replace(".h5", "_parsed.h5")

    probabilities_h5file = \
        tables.open_file(subclade_final_probabilities_filename, mode="r")
    table = probabilities_h5file.root.probabilities.pairs
    occurrences = set()
    for pair in table.where('''(occurrences > 0)'''):
        occurr = int(10 ** int(log10(pair['occurrences'])))
        occurrences.add(occurr)
        if pair['occurrences'] / occurr >= 5:
            occurrences.add(int(float(occurr * 10) / 2))
        if occurr >= 10:
            occurrences.add(int(float(occurr) / 2))
    probabilities_h5file.close()

    return occurrences


def load_occurrences_probs(subclade_id, args):
    if not subclade_id == 0:
        subclade_final_probabilities_filename = \
            os.path.join(args["subclade_clustering_dir"],
                         args["final_probabilities_pattern"] + "_" +
                         args["subclades_final_probabilities_pattern"]
                         .format(subclade_id))
    else:
        subclade_final_probabilities_filename = \
            args["final_probabilities_filename"]

    subclade_final_probabilities_filename = \
        subclade_final_probabilities_filename.replace(".h5", "_parsed.h5")

    probabilities_h5file = \
        tables.open_file(subclade_final_probabilities_filename, mode="r")
    table = probabilities_h5file.root.probabilities.pairs
    occurrences = defaultdict(list)
    for pair in table.where('''(occurrences > 0)'''):
        occurr = int(10 ** int(log10(pair['occurrences'])))
        if pair['occurrences'] / occurr >= 5:
            occurrences[int(float(occurr * 10) / 2)] \
                .append(pair['prob'])
        else:
            occurrences[occurr]. \
                append(pair['prob'])
    probabilities_h5file.close()

    return occurrences


def update_results_file(permuted_probabilities,
                        moved_pairs_perm,
                        subclade_id,
                        subclade_members,
                        occurrence,
                        args):
    permuted_probabilities.sort()
    n_perm = len(permuted_probabilities)
    try:
        index = np.where(permuted_probabilities > 0.5)[0][0]
        threshold = (n_perm - index + 1) / (n_perm + 1)
    except IndexError:
        threshold = 1. / (n_perm + 1)
    percentile = permuted_probabilities[-1]
    percentile2 = permuted_probabilities[-10]
    mean = np.mean(permuted_probabilities)
    median = np.median(permuted_probabilities)
    with open(args['permutations_file'], "a") as out_f:
        out_f.write(
            "{}\t{}\t{}\t{:.5g}\t{:.5g}"
            "\t{:.5g}\t{:.5g}\t{:.5g}"
            "\t{:.3g}\n".format(subclade_id,
                                len(subclade_members),
                                occurrence,
                                percentile,
                                percentile2,
                                mean,
                                median,
                                threshold,
                                np.mean(moved_pairs_perm)))


def update_bootstrap_results_file(bootstrap_means,
                                  subclade_id,
                                  subclade_members,
                                  occurrence,
                                  args):
    bootstrap_means.sort()
    n_perm = len(bootstrap_means)
    try:
        index = np.where(bootstrap_means > 0.5)[0][0]
        threshold = (n_perm - index + 1) / (n_perm + 1)
    except IndexError:
        threshold = 1. / (n_perm + 1)
    percentile = bootstrap_means[-1]
    percentile2 = bootstrap_means[-10]
    mean = np.mean(bootstrap_means)
    median = np.median(bootstrap_means)
    with open(args['bootstrapping_pairs_file'], "a") as out_f:
        out_f.write(
            "{}\t{}\t{}\t{:.5g}\t{:.5g}"
            "\t{:.5g}\t{:.5g}\t{:.5g}\n".format(subclade_id,
                                                subclade_members,
                                                occurrence,
                                                percentile,
                                                percentile2,
                                                mean,
                                                median,
                                                threshold))
