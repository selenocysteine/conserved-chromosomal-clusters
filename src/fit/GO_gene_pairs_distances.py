import itertools
import os
import numpy as np
from utilities import file_utilities
from genome_parsing import load_genome_annotations
from genome_parsing import genomic_pairwise_distances
from collections import OrderedDict
import logging

logger = logging.getLogger(__name__)


def import_GO_pairs_distances(pfam_annotations_file,
                              GO_to_pfam_mapping,
                              pfam_families,
                              args,
                              return_annotations=False):
    # Get the size of the genome (needed to compute pairwise distances)
    organism_name = \
        file_utilities \
            .get_organism_from_file(pfam_annotations_file,
                                    args["pfam_annotations_filename_suffix"])
    genome_sizes = args['genome_sizes']
    genome_size = genome_sizes[organism_name][0][1]
    genome_id = genome_sizes[organism_name][0][0]

    # Get the gene location information for this organism
    coordinates_file, gene_to_chr_mapping_file = \
        file_utilities \
            .get_coord_file_from_id(organism_name,
                                    args["parsed_coordinates_filename_suffix"],
                                    args["gene_locations_filename_suffix"],
                                    args["out_dir"])

    # Computes the pairwise distances of GO pairs on this genome
    distances = get_gene_pairs_distances(coordinates_file,
                                         gene_to_chr_mapping_file,
                                         pfam_annotations_file,
                                         GO_to_pfam_mapping,
                                         pfam_families,
                                         genome_id,
                                         genome_size,
                                         args,
                                         return_annotations)

    return (organism_name, genome_size) + distances


def get_gene_pairs_distances(coordinates_file,
                             gene_to_chr_mapping_file,
                             pfam_annotations_file,
                             GO_to_pfam_mapping,
                             pfam_families,
                             genome_id,
                             genome_size,
                             args,
                             return_annotations=False):
    # Name of the output file to save pairwise distances of GO pairs (optional)
    if args["save_dist"] is True:
        GO_pairs_distances_file_name = \
            file_utilities \
                .new_suffix_file(pfam_annotations_file,
                                 args["pfam_annotations_filename_suffix"],
                                 "_GO" + \
                                 args["pfam_annotations_filename_suffix"],
                                 args["single_distances_dir"],
                                 file_extension="")
        if not os.path.exists(GO_pairs_distances_file_name):
            repeat_file = True
        else:
            repeat_file = False
    else:
        GO_pairs_distances_file_name = None
        repeat_file = True

    # If this file hasn't already been generated in a previous run:
    if repeat_file is True or args["force"] is True:
        logger.verbose("Started computing chromosomal distances of "
                       "GO-related Pfam domain pairs from file %s",
                       pfam_annotations_file)

        # Load the Pfam annotations for the genes in this genome,
        # and organises themaccording to their GO biological process annotations
        gene_coordinates = load_genome_annotations.import_parsed_coordinates(
            coordinates_file)

        domain_hits = load_genome_annotations.import_pfam_hits(
            pfam_annotations_file,
            file_utilities.load_json(
                gene_to_chr_mapping_file),
            args["evalue"],
            pfam_families,
            genome_id)

        # Get list of all gene pairs in the genome with shared GO annotations
        if return_annotations is False:
            GO_pairs = get_GO_pairs_set(domain_hits, GO_to_pfam_mapping)
        else:
            GO_pairs, GO_pairs_annotations = \
                get_GO_pairs_set(domain_hits,
                                 GO_to_pfam_mapping,
                                 return_annotations)

        # Compute pairwise distances for gene pairs with shared GO annotations
        GO_pairs_distances = \
            genomic_pairwise_distances.pairwise_distances(GO_pairs,
                                                          gene_coordinates,
                                                          genome_size)

        # If we want to look at non-GO pairs
        if args["non_go"] is True:
            all_pfam_genes = set(list(itertools.chain(*domain_hits.values())))
            all_pairs = itertools.product(all_pfam_genes,
                                          all_pfam_genes)
            other_pairs = list(set(all_pairs) - set(GO_pairs))
            other_pairs_distances = \
                genomic_pairwise_distances.pairwise_distances(other_pairs,
                                                              gene_coordinates,
                                                              genome_size)
        if len(GO_pairs_distances) == 0:
            logger.verbose("No GO-related domain pairs for file %s available",
                           pfam_annotations_file)
        else:
            # Save results
            if args["save_dist"] is True:
                save_pairwise_distances_file(GO_pairs_distances_file_name,
                                             GO_pairs_distances)

                if args["non_go"] is True:
                    args['non_go_pairs_distances_file_name'] = \
                        GO_pairs_distances_file_name.replace("_GO", "_not_GO")
                    save_pairwise_distances_file(
                        args['non_go_pairs_distances_file_name'],
                        other_pairs_distances)

        logger.verbose("Chromosomal distances of GO-related Pfam domain "
                       "pairs for file %s computed", pfam_annotations_file)

    else:
        np.load(GO_pairs_distances_file_name)
        logger.verbose("Chromosomal distances of GO-related Pfam domain pairs "
                       "for file %s already computed",
                       pfam_annotations_file)

    if args['non_go'] is True:
        if return_annotations is True:
            return GO_pairs_distances, other_pairs_distances, \
                   GO_pairs_annotations
        else:
            return GO_pairs_distances, other_pairs_distances
    else:
        if return_annotations is True:
            return GO_pairs_distances, GO_pairs_annotations
        else:
            return (GO_pairs_distances,)


def get_GO_pairs_set(domain_hits, GO_to_pfam_mapping,
                     return_annotations=False):
    """Map each predicted gene to its putative GO annotations,
    according to the GO annotations of the Pfam domains
    mapping_tables in that gene."""
    GO_pairs = list()
    if return_annotations is True:
        annotations_pairs = OrderedDict()
        gene_to_pfam = OrderedDict()

    for GO_ID in GO_to_pfam_mapping:
        GO_genes = list()
        pfam_IDs = GO_to_pfam_mapping[GO_ID]

        for pfam_ID in pfam_IDs:
            pfam_ID = pfam_ID.replace("PF", "")
            try:
                these_genes = domain_hits[pfam_ID]
                GO_genes = GO_genes + these_genes
            except KeyError:
                these_genes = []
                pass
            if return_annotations is True \
                    and len(these_genes) > 0:
                for gene in these_genes:
                    try:
                        gene_to_pfam[gene].add(pfam_ID)
                    except KeyError:
                        gene_to_pfam[gene] = set([pfam_ID])

        if len(GO_genes) > 1:
            GO_genes_combinations = [tuple(sorted(x)) for x in
                                     itertools.product(GO_genes, GO_genes)]
            GO_pairs = GO_pairs + GO_genes_combinations
            if return_annotations is True:
                for GO_pair in GO_genes_combinations:
                    pfam_domains = gene_to_pfam[GO_pair[0]] \
                        .union(gene_to_pfam[GO_pair[1]])
                    GO_pair = "_".join(sorted(GO_pair))
                    try:
                        annotations_pairs[GO_pair][0].add(GO_ID)
                        annotations_pairs[GO_pair][1].union(pfam_domains)
                    except KeyError:
                        annotations_pairs[GO_pair] = \
                            [set([GO_ID]), pfam_domains]

    if return_annotations is True:
        return list(set(GO_pairs)), annotations_pairs
    else:
        return list(set(GO_pairs))


def save_pairwise_distances_file(out_file_name,
                                 pairwise_distances_values):
    np.savez(out_file_name, pairwise_distances_values)

    return 0
