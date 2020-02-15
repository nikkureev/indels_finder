from Bio import SearchIO
from Bio import SeqIO
import operator
import os
from Bio.Blast.Applications import NcbiblastnCommandline


class gomSegment():

    def __init__(self, segment_id, start, end, genome, genome_sub, sequence):
        self.segment_id = segment_id
        self.start = start
        self.end = end
        self.genome = genome
        self.genome_sub = genome_sub
        self.sequence = sequence


class indel():

    def __init__(self, id1, id2, start1, end1, start2, end2, start, end, genome, size, id1_seq, id2_seq):
        self.id1 = id1
        self.id2 = id2
        self.start1 = start1
        self.end1 = end1
        self.start2 = start2
        self.end2 = end2
        self.genome = genome
        self.size = size
        self.start = start
        self.end = end
        self.id1_seq = id1_seq
        self.id2_seq = id2_seq


class indel_mod():

    def __init__(self, start, end, size, path, wide, sequence='',
                 indel_type=''):
        self.start = start
        self.end = end
        self.size = size
        self.sequence = sequence
        self.path = path
        self.indel_type = indel_type
        for seq_file in SeqIO.parse(self.path, 'fasta'):
            self.sequence = seq_file.seq[start - wide:end + wide]


class sizer():

    def __init__(self, a, b, size=0):
        self.a = a
        self.b = b
        self.size = size + (self.b - self.a)


class blast_record():

    def __init__(self, inp_file):
        self.blast_qr = SearchIO.read(inp_file, "blast-xml")
        self.analysis_bit = []
        self.analysis_sim = []
        self.analysis_len = []
        for hsp in self.blast_qr:
            for h in hsp:
                self.analysis_bit.append(h.bitscore)
                self.analysis_sim.append((float(h.ident_num) * 100) / h.query_span)
                self.analysis_len.append(h.aln_span)
        self.analysis_bit.sort()
        self.analysis_sim.sort()
        self.analysis_len.sort()

    def bit(self):
        return self.analysis_bit

    def sim(self):
        return self.analysis_sim

    def len(self):
        return self.analysis_len

    def blast_r(self):
        return self.blast_qr


def keyfunc(item):
    return item[1]


def similarity_filter(inp_file, level, selecting_type):
    main_in = blast_record(inp_file)
    table = []
    for hsp in main_in.blast_r():
        for h in hsp:
            if selecting_type == 'bitscore':
                if float(h.bitscore) > int(level):
                    table.append(h)
            elif selecting_type == 'similarity':
                if (float(h.ident_num) * 100) / h.query_span > int(level):
                    table.append(h)
    print('Similarity filtering is done:', len(table))
    return table


def index(in_list, genome1, genome2):
    genome_list_1 = []
    genome_list_2 = []
    number = 0
    for hits in in_list:

        genome_list_1.append(gomSegment(number, int(hits.query_range[0]), int(hits.query_range[1]), genome1, genome2,
                                        str(hits.query.seq)))
        genome_list_2.append(gomSegment(number, int(hits.hit_range[0]), int(hits.hit_range[1]), genome2, genome1,
                                        str(hits.hit.seq)))
        number += 1

    genome_list_1.sort(key=operator.attrgetter('start'))
    genome_list_2.sort(key=operator.attrgetter('start'))
    genome_list = [genome_list_1, genome_list_2]
    return genome_list


def indel_search(in_list, min_len, max_len, par='matrix'):
    def linear_way(in_list, min_len, max_len):
        indel_list = []
        for frg_n in range(len(in_list) - 1):
            dist_list1 = [
                sizer(in_list[frg_n].start, in_list[frg_n + 1].start),
                sizer(in_list[frg_n].end, in_list[frg_n + 1].start),
                sizer(in_list[frg_n].start, in_list[frg_n + 1].end),
                sizer(in_list[frg_n].end, in_list[frg_n + 1].end)]
            dist_list1.sort(key=operator.attrgetter('size'))
            if dist_list1[0].size > 0 and dist_list1[1].size > 0 and dist_list1[2].size > 0 and dist_list1[3].size > 0:
                if max_len >= dist_list1[0].size >= min_len:
                    indel_list.append(indel(in_list[frg_n + 1].segment_id,
                                            in_list[frg_n].segment_id,
                                            in_list[frg_n + 1].start,
                                            in_list[frg_n + 1].end,
                                            in_list[frg_n].start,
                                            in_list[frg_n].end,
                                            dist_list1[0].a,
                                            dist_list1[0].b,
                                            in_list[frg_n + 1].genome,
                                            dist_list1[0].size,
                                            in_list[frg_n].sequence,
                                            in_list[frg_n + 1].sequence))
        print('Indel search is done:', len(indel_list))
        return indel_list

    def matrix_way(in_list, min_len, max_len):
        indel_list = []
        indel_matrix = [[[] for k in range(len(in_list))] for f in range(len(in_list))]
        for i in range(len(in_list) - 1):
            for j in range(i, len(in_list)):
                dist_list1 = [
                    sizer(in_list[i].start, in_list[j].start),
                    sizer(in_list[i].end, in_list[j].start),
                    sizer(in_list[i].start, in_list[j].end),
                    sizer(in_list[i].end, in_list[j].end)]
                dist_list1.sort(key=operator.attrgetter('size'))
                if dist_list1[0].size > 0 and dist_list1[1].size > 0 and dist_list1[2].size > 0 \
                        and dist_list1[3].size > 0:
                    if max_len >= dist_list1[0].size >= min_len:
                        indel_matrix[i][j].append(indel(in_list[j].segment_id,
                                                        in_list[i].segment_id,
                                                        in_list[j].start,
                                                        in_list[j].end,
                                                        in_list[i].start,
                                                        in_list[i].end,
                                                        dist_list1[0].a,
                                                        dist_list1[0].b,
                                                        in_list[j].genome,
                                                        dist_list1[0].size,
                                                        in_list[i].sequence,
                                                        in_list[j].sequence))
        for i in range(len(indel_matrix)):
            for j in range(len(indel_matrix[0])):
                if indel_matrix[i][j] != []:
                    indel_list.append(indel_matrix[i][j][0])
        print('Indel search is done:', len(indel_list))
        return indel_list

    if par == 'matrix':
        a = matrix_way(in_list, min_len, max_len)
    else:
        a = linear_way(in_list, min_len, max_len)

    return a


def insertions_intrcept(indel_list1, indel_list2, wide, file1, file2):
    insertion_list = []

    def loc_func(file, specific_id, number, indel):
        insertion_list.append(indel_mod(indel.start, indel.end,
                                        indel.size, indel.genome, wide,
                                        indel_type='insertion'))
        with open(file, 'a') as f:
            f.write('>ins_' + str(number) + '_' + specific_id + '\n')
            f.write(str(insertion_list[-1].sequence) + '\n' + '\n')
    number = 0
    for indels1 in indel_list1:
        loc_func(file1, 'aligned to ' + str(indels1.genome), number, indels1)
        number += 1
    for indels2 in indel_list2:
        loc_func(file2, 'aligned to ' + str(indels2.genome), number, indels2)
        number += 1

    print('Insertion search is done:', len(insertion_list))
    return insertion_list


def deletion_position(indel_list, wide, ak, file):
    deletion_list = []
    number = 0
    for indels in indel_list:
        deletion_list.append(indel_mod(indels.start, indels.end, indels.size,
                                       indels.genome, wide, indel_type='deletion'))
        with open(file, 'a') as f:
            f.write('>del_' + str(number) + ' ' + 'del_' + str(number) + '_' + str(ak) + '\n')
            f.write(str(deletion_list[-1].sequence) + '\n' + '\n')
            number += 1
    print('Deletion search is done:', len(deletion_list))
    return deletion_list


def alignfunc(in_list, out_f):
    in_file = in_list[0]
    db_file = in_list[1]
    cmd = r'C:/Theileria_seq/ncbi-blast-2.9.0+/bin/makeblastdb -in %s -dbtype nucl' % db_file
    os.system(cmd)
    proga = r'C:/Theileria_seq/ncbi-blast-2.9.0+/bin/blastn'
    blast = NcbiblastnCommandline(proga, query=in_file, db=db_file, out=out_f, outfmt=5,
                                  word_size=12, evalue=0.001)
    stdout, stderr = blast()


def align_parsing(path, select_level):
    print('started')
    sum_bitscore = []
    sum_bitscore_selected = []
    blast_result = SearchIO.parse(path, "blast-xml")
    for hsp in blast_result:
        for hs in hsp:
            for h in hs:
                sum_bitscore.append(h.bitscore)
                if h.bitscore > select_level:
                    print(h.bitscore)
                    print(h.query)
                    print(h.hit_description)
                    print()
                    sum_bitscore_selected.append(h.bitscore)
    return_list = [sum_bitscore, sum_bitscore_selected]
    print('ended')
    return return_list


def main(first_genome, second_genome, third_genome, score_filtering_level,
         min_len, max_len, wide, path1, path2, path3, path4, tip, out_align_file_1, out_align_file_2, out_align_file_3,
         align_select_level, align_out_f_1, align_out_f_2):

    a1 = similarity_filter(align_out_f_1, score_filtering_level, selecting_type=str(tip))
    a2 = similarity_filter(align_out_f_2, score_filtering_level, selecting_type=str(tip))

    b1 = index(a1, first_genome, second_genome)
    b2 = index(a2, first_genome, third_genome)

    c11 = indel_search(b1[0], min_len, max_len)
    c21 = indel_search(b1[1], min_len, max_len)
    c12 = indel_search(b2[0], min_len, max_len)
    c22 = indel_search(b2[1], min_len, max_len)

    f1 = open(path1, 'w')
    f2 = open(path2, 'w')
    f3 = open(path3, 'w')
    f4 = open(path4, 'w')
    insertions_intrcept(c11, c12, wide, path1, path4)
    deletion_position(c21, wide, second_genome, path2)
    deletion_position(c22, wide, third_genome, path3)
    f1.close()
    f2.close()
    f3.close()
    f4.close()

    group1 = [path1, path3]
    group2 = [path4, path2]
    group3 = [path2, path3]

    alignfunc(group1, out_align_file_1)
    alignfunc(group2, out_align_file_2)
    alignfunc(group3, out_align_file_3)

    align_parsing(out_align_file_1, align_select_level)
    align_parsing(out_align_file_2, align_select_level)
    align_parsing(out_align_file_3, align_select_level)


main('C:/Theileria/MATRIX/THEILERIA_ANNULATA_1_CHR_REV.fasta',
     'C:/Theileria/MATRIX/THEILERIA_ORIENTALIS_1_CHR.fasta',
     'C:/Theileria/MATRIX/THEILERIA_PARVA_1_CHR.fasta',
     83, 80, 8000, 1000,
     'C:/Theileria/MATRIX/insertions_file_1_2CHR.txt',
     'C:/Theileria/MATRIX/deletions_file_1_2CHR.txt',
     'C:/Theileria/MATRIX/deletions_file_2_2CHR.txt',
     'C:/Theileria/MATRIX/insertions_file_2_2CHR.txt',
     'similarity',
     'C:/Theileria/MATRIX/1_align_file.xml',
     'C:/Theileria/MATRIX/2_align_file.xml',
     'C:/Theileria/MATRIX/3_align_file.xml',
     1500, 'C:/Theileria/MATRIX/AnnOri_1CHR_Alignment.xml',
     'C:/Theileria/MATRIX/AnnPar_1CHR_Alignment.xml')
print('=======================================')
