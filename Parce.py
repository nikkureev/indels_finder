from Bio.Blast import NCBIXML
import matplotlib.pyplot as plt


def genome_maker(file):
    s = []
    with open(file, 'r') as f:
        for lines in f:
            l1 = []
            l1 += str(lines)
            for i in l1:
                if i == '\n':
                    l1.remove(i)
            s += l1
    return s


def space_counter(start_list, end_list):
    space_list = []
    for i in range(len(start_list)):
        space = start_list[i] - end_list[i]
        space_list.append(space)
    space_list.remove(space_list[0])
    space_list.append(0)
    return space_list


def space_maker(list):
    d_list = []
    for i in list:
        d_list.append(' ' * int(i))
    return d_list



ann_genome = genome_maker('C:\Theileria\THEILERIA_ANNULATA_1_CHR.txt')
par_genome = genome_maker('C:\Theileria\THEILERIA_PARVA_1_CHR.txt')
ori_genome = genome_maker('C:\Theileria\THEILERIA_ORIENTALIS_1_CHR.txt')


query_start_list, sbjct_start_list = [], []
query_end_list, sbjct_end_list = [], []
result_handle_1 = open("C:\Theileria\Annulata_Parva_Full_Alignment.xml", 'r')
result_handle_2 = open("C:\Theileria\Annulata_Orientalis_Full_Alignment.xml", 'r')
result_handle_3 = open("C:\Theileria\Parva_Orientalis_Full_Alignment.xml", 'r')
blast_record = NCBIXML.read(result_handle_1)
for alignment in blast_record.alignments:
    for hsp in alignment.hsps:
        query_start_list.append(hsp.query_start)
        sbjct_start_list.append(hsp.sbjct_start)
        query_end_list.append(hsp.query_end)
        sbjct_end_list.append(hsp.sbjct_end)

query_start_list.sort()
query_end_list.sort()
sbjct_start_list.sort()
sbjct_end_list.sort()
query_end_list.insert(0, 0)
sbjct_end_list.insert(0, 0)

space_list_1 = space_counter(query_start_list, query_end_list)
space_list_2 = space_counter(sbjct_start_list, sbjct_end_list)
j = 0
for i in space_list_1:
    if i > 100:
        j += 1
print(j)


#uniq_list = []
#for i in range(len(space_list_1)):
    #if space_list_1[i] > space_list_2[i]:
        #uniq_list.append(space_list_1[i])
    #else:
        #uniq_list.append(space_list_2[i])

#print('Theileria annulata', end='\n')
#for i in ann_genome[query_start_list[0]: query_end_list[-1]]:
    #print(i, end='')
#print()
#print()

#n = 0
#for j in query_start_list:
    #for alignment in blast_record.alignments:
        #for hsp in alignment.hsps:
            #if hsp.query_start == j:
                #print(hsp.query, end='')
                #print(' ' * 15, end='')
                #print('(' + str(space_list_1[n]) + ')', end='')
                #t = 15 + (len(str(uniq_list[n])) - len(str(space_list_1[n])))
                #print(' ' * t, end='')
                #n += 1

#print()
#k = 0
#for j in query_start_list:
    #for alignment in blast_record.alignments:
        #for hsp in alignment.hsps:
            #if hsp.query_start == j:
                #print(hsp.match, end='')
                #print(' ' * 15, end='')
                #e = 17 + len(str(uniq_list[k]))
                #print(' ' * e, end='')
                #k += 1
#print()
#m = 0
#for j in query_start_list:
    #for alignment in blast_record.alignments:
        #for hsp in alignment.hsps:
            #if hsp.query_start == j:
                #print(hsp.sbjct, end='')
                #print(' ' * 15, end='')
                #print('(' + str(space_list_2[m]) + ')', end='')
                #u = 15 + (len(str(uniq_list[m])) - len(str(space_list_2[m])))
                #print(' ' * u, end='')
                #m += 1
#print()
#print()
#print('Theileria parva', end='\n')
#for i in par_genome[sbjct_start_list[0]: sbjct_end_list[-1]]:
    #print(i, end='')
