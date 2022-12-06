import pickle 
from sequence import * 

class_1_in = open('./notebooks/class_1_unreviewed_annots', 'rb')
class_1_annots = pickle.load(class_1_in)
class_1_in.close()

class_2_in = open('./notebooks/class_2_unreviewed_annots', 'rb')
class_2_annots = pickle.load(class_2_in)
class_2_in.close()


print(class_2_annots)

class_1_seqs = readFastaFile('./datasets/kari_ec_1_1_1_86_ec_1_1_1_382_ec_1_1_1_383_unreviewed/subsets/classI_kari/kari_ec_1_1_1_86_ec_1_1_1_382_ec_1_1_1_383_unreviewed_classI_kari.fasta')
class_2_seqs = readFastaFile('./datasets/kari_ec_1_1_1_86_ec_1_1_1_382_ec_1_1_1_383_unreviewed/subsets/classII_kari/kari_ec_1_1_1_86_ec_1_1_1_382_ec_1_1_1_383_unreviewed_classII_kari.fasta')

def summarise_class_2_sites_swissprot(sequences, annots):

    n_seqs = []
    c1_seqs = []
    c2_seqs = []

    for seq in sequences:

        ##################### 
        # N TERMINUS 
        #####################
        try:
            n_range= annots[seq.name][0]['KARI N-terminal Rossmann'][0].split('..')

        except:
            continue 
        
        print('check1')
        if n_range[0][0] == '<': 
            
            n_start = int(n_range[0][1:])
        
        else: n_start = int(n_range[0])
        
        n_end = int(n_range[1])

        n_dom = seq[n_start-1: n_end]
        name = seq.name + '_class_2_n_domain'
        n_seqs.append(Sequence(n_dom, name=name))
        
        

        ##################### 
        # C1 TERMINUS 
        #####################
        try:

            c1_range= annots[seq.name][0]['KARI C-terminal knotted 1'][0].split('..')
        
        except:
            continue 
        

        print('check2')
        c1_start = int(c1_range[0])
        c1_end = int(c1_range[1])

        c1_dom = seq[c1_start-1: c1_end]
        name = seq.name + '_class_2_c1_domain'
        c1_seqs.append(Sequence(c1_dom, name=name))

        
        ##################### 
        # C2 TERMINUS 
        #####################

        try:

            c2_range= annots[seq.name][0]['KARI C-terminal knotted 2'][0].split('..')

        
        except:
            continue 
        

        print('check3')
        c2_start = int(c2_range[0])

        if c2_range[1][0] == '>': c2_end = int(c2_range[1][1:])

        else: c2_end = int(c2_range[1])
    
        c2_dom = seq[c2_start-1: c2_end]
        name = seq.name + '_class_2_c2_domain'
        
        c2_seqs.append(Sequence(c2_dom, name=name))


    writeFastaFile('./datasets/kari_ec_1_1_1_86_ec_1_1_1_382_ec_1_1_1_383_unreviewed/subsets/classII_kari/kari_class_2_unreviewed_n_domains.fasta', n_seqs)
    writeFastaFile('./datasets/kari_ec_1_1_1_86_ec_1_1_1_382_ec_1_1_1_383_unreviewed/subsets/classII_kari/kari_class_2_unreviewed_c1_domains.fasta', c1_seqs)
    writeFastaFile('./datasets/kari_ec_1_1_1_86_ec_1_1_1_382_ec_1_1_1_383_unreviewed/subsets/classII_kari/kari_class_2_unreviewed_c2_domains.fasta', c2_seqs)


    return n_seqs, c1_seqs, c2_seqs 

def summarise_class_1_sites(sequences, annots):

    n_seqs = []
    c1_seqs = []

    for seq in sequences:

        ##################### 
        # N TERMINUS 
        #####################
        try:
            n_range= annots[seq.name][0]['KARI N-terminal Rossmann'][0].split('..')
        except:
            continue 
        if n_range[0][0] == '<': 
            
            n_start = int(n_range[0][1:])
        
        else: n_start = int(n_range[0])
        
        n_end = int(n_range[1])

        n_dom = seq[n_start-1: n_end]
        name = seq.name + '_class_1_n_domain'
        n_seqs.append(Sequence(n_dom, name=name))
        
        

        ##################### 
        # C1 TERMINUS 
        #####################
        try:
            c1_range= annots[seq.name][0]['KARI C-terminal knotted'][0].split('..')

        except:
            continue 
      
        if c1_range[1][0] == '>': c1_end = int(c1_range[1][1:])
        else: c1_end = int(c1_range[1])

        c1_start = int(c1_range[0])

        c1_dom = seq[c1_start-1: c1_end]
        name = seq.name + '_class_1_c1_domain'
        c1_seqs.append(Sequence(c1_dom, name=name))

       

    writeFastaFile('./datasets/kari_ec_1_1_1_86_ec_1_1_1_382_ec_1_1_1_383_unreviewed/subsets/classI_kari/kari_class_1_unreviewed_n_domains.fasta', n_seqs)
    writeFastaFile('./datasets/kari_ec_1_1_1_86_ec_1_1_1_382_ec_1_1_1_383_unreviewed/subsets/classI_kari/kari_class_1_unreviewed_c1_domains.fasta', c1_seqs)

    return n_seqs, c1_seqs
        
#I_n, I_c1 = summarise_class_1_sites(class_1_seqs, class_1_annots)

#II_n, II_c1, II_c2 = summarise_class_2_sites(class_2_seqs, class_2_annots)


def sort_class_2_residues(sequences, annots):

    mg_c1 = []
    mg_c2 = []
    substrates = []

    for seq in sequences:
        
        ##################### 
        # C1 TERMINUS 
        #####################
       
        c1_range= annots[seq.name][0]['KARI C-terminal knotted 1'][0].split('..')
    
        c1_start = int(c1_range[0])
        c1_end = int(c1_range[1])

        ##################### 
        # C2 TERMINUS 
        #####################

        c2_range= annots[seq.name][0]['KARI C-terminal knotted 2'][0].split('..')
    
        c2_start = int(c2_range[0])

        if c2_range[1][0] == '>': c2_end = int(c2_range[1][1:])
        else: c2_end = int(c2_range[1])
       
        #######################
        # IDENTIFY METAL RESIDUES
        #######################

        for res in annots[seq.name][1]['Mg(2+)']:

            res = int(res)
            
            if res >= c1_start and res <= c1_end:
                mg_c1.append(res)
            
            elif res >= c2_start and res <= c2_end:
                mg_c2.append(res)
        
        if 'substrate' in annots[seq.name][1].keys():
            substrates.append(int(annots[seq.name][1]['substrate'][0]))

    return mg_c1, mg_c2, substrates
        
        


#mg_c1, mg_c2, substrates = sort_class_2_residues(class_2_seqs, class_2_annots)




def summarise_class_2_sites_unreviewed(sequences, annots):

    n_seqs = []
    c1_seqs = []
    c2_seqs = []

    for seq in sequences:

        ##################### 
        # N TERMINUS 
        #####################
        try:
            n_range= annots[seq.name][0]['KARI N-terminal Rossmann'][0].split('..')

        except:
            continue 
        
      
        if n_range[0][0] == '<': 
            
            n_start = int(n_range[0][1:])
        
        else: n_start = int(n_range[0])
        
        n_end = int(n_range[1])

        n_dom = seq[n_start-1: n_end]
        name = seq.name + '_class_2_n_domain'
        n_seqs.append(Sequence(n_dom, name=name))
        
        

        ##################### 
        # C1 TERMINUS 
        #####################

        try:
            c1_range= annots[seq.name][0]['KARI C-terminal knotted'][0].split('..')
        except:
           
            continue 
        
        try:
            c2_range= annots[seq.name][0]['KARI C-terminal knotted'][1].split('..')
        except:
            print(annots[seq.name][0])
            continue

        c1_start = int(c1_range[0])
        c1_end = int(c1_range[1])

        c1_dom = seq[c1_start-1: c1_end]
        name = seq.name + '_class_2_c1_domain'
        c1_seqs.append(Sequence(c1_dom, name=name))

        
        c2_start = int(c2_range[0])

        if c2_range[1][0] == '>': c2_end = int(c2_range[1][1:])

        else: c2_end = int(c2_range[1])
    
        c2_dom = seq[c2_start-1: c2_end]
        name = seq.name + '_class_2_c2_domain'
        
        c2_seqs.append(Sequence(c2_dom, name=name))


    writeFastaFile('./datasets/kari_ec_1_1_1_86_ec_1_1_1_382_ec_1_1_1_383_unreviewed/subsets/classII_kari/kari_class_2_unreviewed_n_domains.fasta', n_seqs)
    writeFastaFile('./datasets/kari_ec_1_1_1_86_ec_1_1_1_382_ec_1_1_1_383_unreviewed/subsets/classII_kari/kari_class_2_unreviewed_c1_domains.fasta', c1_seqs)
    writeFastaFile('./datasets/kari_ec_1_1_1_86_ec_1_1_1_382_ec_1_1_1_383_unreviewed/subsets/classII_kari/kari_class_2_unreviewed_c2_domains.fasta', c2_seqs)


    return n_seqs, c1_seqs, c2_seqs 

#II_n, II_c1, II_c2 = summarise_class_2_sites_unreviewed(class_2_seqs, class_2_annots)

