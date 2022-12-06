
import os 


##### UNREVIEWED #####

unreviewed_c1_bit_scores = []

unreviewed_c1 = open('./datasets/kari_ec_1_1_1_86_ec_1_1_1_382_ec_1_1_1_383_unreviewed/subsets/classII_kari/c1_results.txt','r')
for line in unreviewed_c1:
    if line[0] != "#":
    
        data = line.split()
        unreviewed_c1_bit_scores.append(data[7])

unreviewed_c1.close()


unreviewed_c2_bit_scores = []

unreviewed_c2 = open('./datasets/kari_ec_1_1_1_86_ec_1_1_1_382_ec_1_1_1_383_unreviewed/subsets/classII_kari/c2_results.txt','r')
for line in unreviewed_c2:
    if line[0] != "#":
    
        data = line.split()
        unreviewed_c2_bit_scores.append(data[7])

unreviewed_c2.close()


##### SWISSPROT #####

swiss_c1_bit_scores = []

swiss_c1 = open('./datasets/kari_ec_1_1_1_86_ec_1_1_1_382_ec_1_1_1_383_swissprot/subsets/classII_kari/class_2_c1.txt','r')
for line in swiss_c1:
    if line[0] != "#":
    
        data = line.split()
        swiss_c1_bit_scores.append(data[7])

swiss_c1.close()


swiss_c2_bit_scores = []

swiss_c2 = open('./datasets/kari_ec_1_1_1_86_ec_1_1_1_382_ec_1_1_1_383_swissprot/subsets/classII_kari/class_2_c2.txt','r')
for line in swiss_c2:
    if line[0] != "#":
    
        data = line.split()
        swiss_c2_bit_scores.append(data[7])

swiss_c2.close()


print(len(unreviewed_c1_bit_scores))


print(len(unreviewed_c2_bit_scores))


print(swiss_c2_bit_scores)

print(swiss_c1_bit_scores)