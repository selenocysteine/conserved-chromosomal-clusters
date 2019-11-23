import numpy as np
from itertools import combinations, product


def do_these_genes_overlap(gene1_coords, gene2_coords):
    """ Check if two genes are overlapping. """
    overlap = max(0,
                  min(gene1_coords[1], gene2_coords[1]) - max(gene1_coords[0],
                                                              gene2_coords[0]))

    return True if overlap > 0 else False


def distance_between_two_genes(gene1_coords, gene2_coords, genome_size):
    # Finds which one of the two genes of the pair is first on the genome:
    first_gene_coord, second_gene_coord = (gene1_coords, gene2_coords) if \
        gene1_coords[0] <= gene2_coords[0] \
        else (gene2_coords, gene1_coords)

    # The distance is defined as the minimum the "clockwise" and "counterclockwise"
    # distance of the two members of a pair on the genome:
    distance1 = second_gene_coord[1] - first_gene_coord[0] + 1
    distance2 = genome_size - second_gene_coord[0] + first_gene_coord[1] + 1
    distance = min(distance1, distance2)

    return distance


def pairwise_distances(pairs_list, gene_coordinates, genome_size):
    """ Compute the pairwise distances on a circular genome of the two members of
    each pair of genes from pairs_list (note: even if pairs_list contains duplicates,
    each pair will only be analysed once)."""
    if not isinstance(pairs_list, list):
        pairs_list = list(pairs_list)

    distances = np.empty(len(pairs_list), dtype="int32")

    for index, pair in enumerate(pairs_list):
        gene1_coords = gene_coordinates[pair[0]]
        gene2_coords = gene_coordinates[pair[1]]
        distance = int(
            distance_between_two_genes(gene1_coords, gene2_coords, genome_size))

        if distance is None or distance < 0:
            raise ValueError("Error")

        distances[index] = distance

    if len(distances) > 1:
        return distances
    else:
        return distances[0]


def get_all_pairs(domain_hits):
    """ Generate a list of all possible pairings of genes in gene_coordinates."""
    all_pairs = combinations(domain_hits.keys(), 2)

    return all_pairs


def check_and_update_gene_distance(gene,
                                   pair_distances,
                                   pair_id,
                                   best_pairs_per_gene):
    new_pair_distance = pair_distances[pair_id]
    try:
        old_pair_distance = pair_distances[best_pairs_per_gene[gene]]
    except KeyError:
        best_pairs_per_gene[gene] = pair_id
        return 0

    if old_pair_distance > new_pair_distance:
        best_pairs_per_gene[gene] = pair_id

    return 0


def compute_domain_pair_distances(domain_pairs,
                                  domain_hits,
                                  gene_coordinates,
                                  genome_size):
    pair_distances = []
    n_fusions_total = []
    n_genes = []

    for domain_pair in domain_pairs:
        domain_pair = domain_pair.split("_")
        domain1 = domain_pair[0]
        domain2 = domain_pair[1]

        genes_domain1 = domain_hits[domain1]
        genes_domain2 = domain_hits[domain2]

        all_gene_pairs_for_this_domain = product(genes_domain1, genes_domain2)
        n_genes.append(len(set(genes_domain1).union(set(genes_domain2))))
        n_fusions = 0

        for gene_pair in all_gene_pairs_for_this_domain:
            gene1 = gene_pair[0]
            gene2 = gene_pair[1]

            pair_distance = pairwise_distances([gene_pair],
                                               gene_coordinates,
                                               genome_size)

            if gene1 == gene2:
                n_fusions = n_fusions + 1

            pair_distances.append(pair_distance)
            n_fusions_total.append(n_fusions)

    return np.array(pair_distances)


def compute_closest_domain_pair_distances(domain_pair,
                                          domain_hits,
                                          gene_coordinates,
                                          genome_size):
    domain1 = domain_pair[0]
    domain2 = domain_pair[1]

    genes_domain1 = set(domain_hits[domain1])
    genes_domain2 = set(domain_hits[domain2])

    n_genes1 = len(genes_domain1)
    n_genes2 = len(genes_domain2)

    fusions = genes_domain1.intersection(genes_domain2)
    n_fusions = len(fusions)
    best_fusions_distances = []
    if n_fusions > 0:
        for i, gene in enumerate(fusions):
            pair_distance = pairwise_distances([[gene, gene]],
                                               gene_coordinates,
                                               genome_size)
            best_fusions_distances.append(pair_distance)
            genes_domain1.remove(gene)
            genes_domain2.remove(gene)

    if len(genes_domain1) == 0 or len(genes_domain2) == 0:
        return np.array(best_fusions_distances, dtype=np.int32),\
               n_fusions, n_genes1, n_genes2

    else:
        all_gene_pairs_for_this_domain = list(product(genes_domain1,
                                                      genes_domain2))
        pair_distances = np.zeros(len(all_gene_pairs_for_this_domain),
                                  dtype=np.int32)
        # dictionary where each gene will have its closest neighbour saved
        best_pairs_per_gene = {}
        for i, gene_pair in enumerate(all_gene_pairs_for_this_domain):
            gene_pair = tuple(sorted(gene_pair))
            gene1 = gene_pair[0]
            gene2 = gene_pair[1]

            pair_distances[i] = pairwise_distances([gene_pair],
                                               gene_coordinates,
                                               genome_size)
            for gene in [gene1, gene2]:
                check_and_update_gene_distance(gene, pair_distances, i,
                                               best_pairs_per_gene)

        best_pairs_per_gene_distances =\
            np.take(pair_distances, list(set(best_pairs_per_gene.values())))

        if n_fusions > 0:
            best_pairs_per_gene_distances = \
                np.concatenate((best_fusions_distances,
                               best_pairs_per_gene_distances))

        return best_pairs_per_gene_distances, n_fusions, n_genes1, n_genes2
