from sequence import * 
from collections import defaultdict


class_1_c1 = readFastaFile('./datasets/kari_ec_1_1_1_86_ec_1_1_1_382_ec_1_1_1_383_swissprot/subsets/classI_kari/kari_class_1_swissprot_c1_domains.fasta')

#print(len(class_1_c1))

class_2_c1_seqs = readFastaFile('./datasets/kari_ec_1_1_1_86_ec_1_1_1_382_ec_1_1_1_383_swissprot/subsets/classII_kari/kari_class_2_swissprot_c1_domains.fasta')

#print(len(class_2_c1_seqs))

class_2_c2_seqs = readFastaFile('./datasets/kari_ec_1_1_1_86_ec_1_1_1_382_ec_1_1_1_383_swissprot/subsets/classII_kari/kari_class_2_swissprot_c2_domains.fasta')

#print(len(class_2_c2_seqs))
k = 7

def count_kmers(seqs, k_len):
    kmers = defaultdict(int)

    for seq in seqs:
        
        for i in range(len(seq) - k_len + 1):
    
            kmers[seq[i:i+k_len]] += 1 

    return kmers

class1_c1_kmers = count_kmers(class_1_c1, k) 

class2_c1_kmers = count_kmers(class_2_c1_seqs, k)

class2_c2_kmers = count_kmers(class_2_c2_seqs, k)



def calc_euc_dist(dom_1:dict, dom_2:dict):

    com = 0
    dif = 0
    dist = []

    all_kmers = [k for k in dom_1.keys()]

    for k in dom_2.keys():
        all_kmers.append(k)

    for kmer in all_kmers:
       
        if kmer in dom_2.keys():
            dom_2_count = dom_2[kmer]
        else:
            dom_2_count = 0 

        if kmer in dom_1.keys():
            dom_1_count = dom_1[kmer]
        else:
            dom_1_count = 0 
     
        dist.append((dom_1_count - dom_2_count)**2)

    sum = 0 
    
    for d in dist:
        sum += d
        
    return math.sqrt(sum)

c1vc1 = calc_euc_dist(class2_c1_kmers, class1_c1_kmers)
print(f'The Euclidean distance between Class I C1 domain and Class II C1 domain is: {c1vc1}')


c1vc2 = calc_euc_dist(class1_c1_kmers, class2_c2_kmers)
print(f'The Euclidean distance between Class I C1 domain and Class II C2 domain is: {c1vc2}')


