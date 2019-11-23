"""Generates the Python dictionaries with the mappings
between Pfam domains and GO annotations."""
import file_utilities
import os

def parse_gene_ontology_go(go_obo_file,
                           gene_ontology_biological_process_file):
    new_GO = False
    this_go = {}
    with open(go_obo_file, "r") as input_file, open(gene_ontology_biological_process_file, "w") as output_file:
        for line in input_file:
            line = line.strip()
            if line == "":
                new_GO = False
                if len(this_go) > 0:
                    if this_go['namespace'] == "biological_process":
                        if not this_go['def'].startswith('\"OBSOLETE.'):
                            output_file.write("{}\t{}\t{}\t{}\n".format(this_go['id'],
                                                                         this_go['name'],
                                                                         this_go['namespace'],
                                                                         this_go['def']))
                this_go = {}
            if line == "[Term]":
                new_GO = True
                this_go = {}
            elif new_GO == True:
                line = line.split(": ")
                this_go[line[0]] = line[1]



def select_pfam_GO_biological_processes(pfam2go_file,
                                        gene_ontology_biological_process_file,
                                        pfam_GOBP_mapping_file):
    GO_BP = []
    with open(gene_ontology_biological_process_file, 'r') as input_file:
        for line in input_file:
            line = line.strip().split("\t")
            if line[2] == "biological_process":
                GO_BP.append(line[0])

    with open(pfam2go_file, 'r') as input_file, \
            open(pfam_GOBP_mapping_file, 'w') as output_file:
        for line in input_file:
            if line.startswith("!"):
                output_file.write(line)
            elif line.strip().split(" ")[-1] in GO_BP:
                output_file.write(line)


def create_pfam_GO_mapping_dictionaries(pfam_GOBP_mapping_file):
    pfam_to_GO_mapping = {}
    GO_to_pfam_mapping = {}

    with open(pfam_GOBP_mapping_file, "r") as input_file:
        for line in input_file:
            if not line.startswith("!"):
                line = line.strip().split()
                pfamID = line[0].split(":")[1]
                GO_ID = line[-1]

                if pfamID not in pfam_to_GO_mapping.keys():
                    pfam_to_GO_mapping[pfamID] = [GO_ID]
                else:
                    pfam_to_GO_mapping[pfamID].append(GO_ID)

                if GO_ID not in GO_to_pfam_mapping.keys():
                    GO_to_pfam_mapping[GO_ID] = [pfamID]
                else:
                    GO_to_pfam_mapping[GO_ID].append(pfamID)

    return pfam_to_GO_mapping, GO_to_pfam_mapping


def load_pfam_domains_descriptions(pfam_mapping_file="Pfam-A.clans.tsv"):
    pfam_domains_descriptions = {}

    with open(pfam_mapping_file, "r") as input_file:
        for line in input_file:
            line = line.strip().split("\t")
            pfam_domains_descriptions[line[0]] = line[1:]

    return pfam_domains_descriptions


def load_pfam_families(pfam_mapping_file="Pfam-A.clans.tsv"):
    pfam_families = {}

    with open(pfam_mapping_file, "r") as input_file:
        for line in input_file:
            line = line.strip().split("\t")
            if len(line[1]) > 0:
                pfam_families[line[0].replace("PF", "")] = line[1]
    return pfam_families


def main():
    path = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "../mapping_tables"))
    parse_gene_ontology_go(go_obo_file=os.path.join(path, "go.obo"),
                           gene_ontology_biological_process_file=os.path.join(path, "goDescription.txt"))
    select_pfam_GO_biological_processes(pfam2go_file=os.path.join(path, "pfam2go"),
                                        gene_ontology_biological_process_file=os.path.join(path, "goDescription.txt"),
                                        pfam_GOBP_mapping_file=os.path.join(path, "pfam2go_biological_processes.txt"))


    pfam_to_GO_mapping, GO_to_pfam_mapping = \
        create_pfam_GO_mapping_dictionaries(os.path.join(path, "pfam2go_biological_processes.txt"))
    file_utilities.save_json(pfam_to_GO_mapping, os.path.join(path, "pfam_to_GO_mapping_dict.json.gz"),
                                                                compress=True)
    file_utilities.save_json(GO_to_pfam_mapping, os.path.join(path,"GO_to_pfam_mapping_dict.json.gz"),
                                                                compress=True)

    pfam_domains_descriptions = load_pfam_domains_descriptions(os.path.join(path, "Pfam-A.clans.tsv"))
    file_utilities.save_json(pfam_domains_descriptions, os.path.join(path, "pfam_details_dict.json.gz"),
                                                                       compress=True)

    pfam_families = load_pfam_families(os.path.join(path, "Pfam-A.clans.tsv"))
    file_utilities.save_json(pfam_families, os.path.join(path, "pfam_families_dict.json.gz"),
                               compress=True)


if __name__ == "__main__":
    main()