from utilities import file_utilities
from collections import OrderedDict
import multiprocessing
manager = multiprocessing.Manager()

def initialise_organisms_parameters_file(organisms_parameters_filename):
    header = ['Assembly', 'Half_genome_size', 'Lambda', 'Phi']
    with open(organisms_parameters_filename, "w") as out_file:
        out_file.write("\t".join(header) + "\n")
    return organisms_parameters_filename


def import_parameters(parameters_filename):
    with open(parameters_filename, "r") as par_file:
        lambd = float(par_file.readline().strip().split("\t")[1])
        phi = float(par_file.readline().strip().split("\t")[1])
    return lambd, phi


def import_genome_sizes(genome_sizes_file):
    return file_utilities.load_json(genome_sizes_file)


def save_general_parameters(lambd, phi, parameters_filename):
    """Save the estimated parameters to a .txt file"""
    headers = ["Estimated_lambda", "Estimated_phi"]
    output = ["{}\t{:.4g}".format(a, b) for a, b in zip(headers, [lambd, phi])]
    with open(parameters_filename, "w") as out_file:
        out_file.write("\n".join(output) + "\n")
    return 0


def save_organism_parameters(organism, genome_size, lambd, phi, organisms_parameters_filename):
    with open(organisms_parameters_filename, "a") as out_file:
        out_file.write("{}\t{}\t{:.4g}\t{:.4g}\n".format(organism, genome_size, lambd, phi))
    return 0


def import_organisms_parameters(organisms_parameters_filename):
    lambdas = OrderedDict()
    phis = OrderedDict()
    with open(organisms_parameters_filename, "r") as inp_file:
        inp_file.readline()
        for line in inp_file:
            line = line.strip().split()
            organism = line[0]
            lambdas[organism] = float(line[2])
            phis[organism] = float(line[3])
    return lambdas, phis
