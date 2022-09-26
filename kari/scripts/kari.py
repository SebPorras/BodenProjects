import sequence
import webservice
import pprint
import pickle
import regex as re
import glob
import os


def get_annot_dict(uniprot_names):
    """
    Searches UniProt for annotations relating to metal binding positions and domains of KARI sequences
    :param uniprot_names: List of UniProt IDs to retrieve information for
    :return: Dictionary mapping UniProt ID -> metal binding positions and domain positions
    """

    # The columns we want to search. See webservice.getUniProtDict for usage and list of possible columns
    cols = ["feature(METAL BINDING)", "feature(DOMAIN EXTENT)"]
    up_dict = webservice.getUniProtDict(uniprot_names, cols)

    # Define regexes to extract the details
    metal_regex = 'METAL (.*?); */note=\"(.*?)\"; */evidence=\"(.*?)\"'
    domain_regex = 'DOMAIN (.*?); */note=\"(.*?)\"; */evidence=\"(.*?)\"'

    annot_dict = {}

    # Store the hits for metal binding and the domains in a dictionary
    for name in uniprot_names:
        if up_dict[name]["feature(METAL BINDING)"] != None and up_dict[name]["feature(DOMAIN EXTENT)"] != None:

            annot_dict[name] = {"metals": {}, "domains": {}}

            # Extract the details
            metals = re.findall(metal_regex, up_dict[name]["feature(METAL BINDING)"])

            # Add the metal information
            for count, metal in enumerate(metals):
                pos = int(metal[0])
                metal_name = metal[1]
                evidence = metal[2]
                annot_dict[name]['metals'][count] = (pos, metal_name)

            domains = re.findall(domain_regex, up_dict[name]["feature(DOMAIN EXTENT)"])

            # Add the domain information
            for dom in domains:
                # TODO: Need to double check what the > and < symbols here actually mean (not just get rid of them!)
                start = int(dom[0].split("..")[0].replace(">", "").replace("<", ""))
                end = int(dom[0].split("..")[1].replace(">", "").replace("<", ""))
                dom_name = dom[1]
                evidence = dom[2]
                annot_dict[name]['domains'][dom_name] = (start, end)

    return annot_dict

# Set this to False if you want to go to UniProt, or True if there are already dictionaries (kari_1.p, kari_2.p) that
# you want to just load
load_from_file = False

# Load files
#kari_1 = sequence.readFastaFile("./files/KARI_C1_Rev_Prot.fasta")
#kari_2 = sequence.readFastaFile("./files/KARI_C2_Rev_Prot.fasta")

# If we don't want to load from file lets make a new dictionary instead
if not load_from_file:
    # Extract just the UniProt IDs
    kari_1_uniprot_names = [seq.name for seq in kari_1]
    kari_2_uniprot_names = [seq.name for seq in kari_2]

    kari_1_dict = get_annot_dict(kari_1_uniprot_names)
    kari_2_dict = get_annot_dict(kari_2_uniprot_names)

    # Pickle the dictionaries to a file
    with open("./kari_1.p", 'wb') as handle:
        pickle.dump(kari_1_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)

    with open("./kari_2.p", 'wb') as handle:
        pickle.dump(kari_2_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)


# Otherwise lets load the dictionaries from a file
else:

    with open('kari_1.p', 'rb') as handle:
        kari_1_dict = pickle.load(handle)

    with open('kari_2.p', 'rb') as handle:
        kari_2_dict = pickle.load(handle)


# Print the lengths of the FASTA and the dictionaries
print(len(kari_1))
print(len(kari_2))
print()
print(len(kari_1_dict))
print(len(kari_2_dict))

# Pretty print the dictionaries
pp = pprint.PrettyPrinter(width=70, compact=True)

# Uncomment the following two lines to print the dictionaries

# pp.pprint(kari_1_dict)
# pp.pprint(kari_2_dict)

# Check that the metals always end up in the same domains

for name, annot in kari_1_dict.items():

    print("Checking " + name)

    # Check the first two metal binding sites are in the first domain
    for i in range(5):
        if not int(annot['domains']['KARI C-terminal knotted'][0]) \
                < int(annot['metals'][i][0]) \
                < int(annot['domains']['KARI C-terminal knotted'][1]):
            print("ERROR:" + name + " was wrong")

for name, annot in kari_2_dict.items():

    print("Checking " + name)

    # Check the first two metal binding sites are in the first domain
    for i in range(3):
        if not int(annot['domains']['KARI C-terminal knotted 1'][0]) \
                < int(annot['metals'][i][0]) \
                < int(annot['domains']['KARI C-terminal knotted 1'][1]):
            print("ERROR:" + name + " was wrong")

    # Check the last two metal binding sites are in the second domain
    for i in range(3, 5):
        if not int(annot['domains']['KARI C-terminal knotted 2'][0]) \
                < int(annot['metals'][i][0]) \
                < int(annot['domains']['KARI C-terminal knotted 2'][1]):
            print("ERROR:" + name + " was wrong")

print('Done')

# Print the sequence content at each of the metal positions

for seq in kari_1:
    if seq.name in kari_1_dict:
        metal_positions = [int(kari_1_dict[seq.name]['metals'][i][0]) - 1 for i in range(5)]
        # print (metal_positions)
        metal_seq = [seq.sequence[pos] for pos in metal_positions]
        print("".join(metal_seq))

print()

for seq in kari_2:
    if seq.name in kari_2_dict:
        metal_positions = [int(kari_2_dict[seq.name]['metals'][i][0]) - 1 for i in range(5)]
        # print (metal_positions)
        metal_seq = [seq.sequence[pos] for pos in metal_positions]
        print("".join(metal_seq))



# Split files up on the basis of where we find an annotated domain

class1_domain1 = []
class2_domain1 = []
class2_domain2 = []


for seq in kari_1:
    if seq.name in kari_1_dict:
        class1_domain1.append(sequence.Sequence("".join(seq.sequence[kari_1_dict[seq.name]['domains']['KARI C-terminal '
                                                                                               'knotted'][0]:
        kari_1_dict[seq.name][
            'domains']['KARI C-terminal knotted'][1]]),
                                         name=seq.name + "class1__domain1"))

for seq in kari_2:
    if seq.name in kari_2_dict:
        class2_domain1.append(
            sequence.Sequence("".join(seq.sequence[kari_2_dict[seq.name]['domains']['KARI C-terminal knotted 1'][
                0]:kari_2_dict[seq.name]['domains']['KARI C-terminal knotted 1'][1]]),
                              name=seq.name + "_class2_domain1"))

for seq in kari_2:
    if seq.name in kari_2_dict:
        class2_domain2.append(
            sequence.Sequence("".join(seq.sequence[kari_2_dict[seq.name]['domains']['KARI C-terminal knotted 2'][0]:]),
                              name=seq.name + "_class2_domain2"))


# Write all the sequences to file
sequence.writeFastaFile("./files/domains/KARI_C1_Rev_Prot_Domain1.fasta", class1_domain1)
sequence.writeFastaFile("./files/domains/KARI_C2_Rev_Prot_Domain1.fasta", class2_domain1)
sequence.writeFastaFile("./files/domains/KARI_C2_Rev_Prot_Domain2.fasta", class2_domain2)
sequence.writeFastaFile("./files/domains/c1d1_c2d1.fasta", class1_domain1 + class2_domain1)
sequence.writeFastaFile("./files/domains/c1d1_c2d2.fasta", class1_domain1 + class2_domain2)
sequence.writeFastaFile("./files/domains/c2d1_c2d2.fasta", class2_domain1 + class2_domain2)
sequence.writeFastaFile("./files/domains/combined.fasta", class1_domain1 + class2_domain1 + class2_domain2)

# Note: This will only work if you have MAFFT installed on the command line
for fasta in glob.glob('./files/domains/*.fasta'):
    os.system("mafft " + fasta + ">" + fasta.split('.fasta')[0] + ".aln")
