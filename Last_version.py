from Bio import SearchIO
from Bio import SeqIO
from Bio.Seq import Seq
import operator
import matplotlib.pyplot as plt


class gomSegment():

    def __init__(self, segment_id, start, end, genome, genome_sub):
        self.segment_id = segment_id
        self.start = start
        self.end = end
        self.genome = genome
        self.genome_sub = genome_sub


class indel():

    def __init__(self, id1, id2, start1, end1, start2, end2, start, end, genome, size):
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


class insertion():

    def __init__(self, start, end, genome, size, path, wide, sequence=''):
        self.start = start
        self.end = end
        self.genome = genome
        self.size = size
        self.sequence = sequence
        self.path = path
        for seq_file in SeqIO.parse(self.path, 'fasta'):
            self.sequence = seq_file.seq[start - wide:end + wide]


class deletion():

    def __init__(self, start, end, genome, size, common_genome_start):
        self.start = start
        self.end = end
        self.genome = genome
        self.size = size
        self.common_genome_start = common_genome_start


class sizer():

    def __init__(self, a, b, size=0):
        self.a = a
        self.b = b
        self.size = size + (self.b - self.a)


def keyfunc(item):
    return item[1]


def similarity_filter(inp_file, level):
    blast_qresult = SearchIO.read(inp_file, "blast-xml")
    table = []
    for hsp in blast_qresult:
        for h in hsp:
            if float(h.bitscore) > int(level):
                table.append(h)
    print('Similarity filtering is done:', len(table))
    return table


def index(in_list, genome1, genome2):
    genome_list_1 = []
    genome_list_2 = []
    number = 0
    for hits in in_list:

        genome_list_1.append(gomSegment(number, int(hits.query_range[0]), int(hits.query_range[1]), genome1, genome2))
        genome_list_2.append(gomSegment(number, int(hits.hit_range[0]), int(hits.hit_range[1]), genome2, genome1))
        number += 1

    genome_list_1.sort(key=operator.attrgetter('start'))
    genome_list_2.sort(key=operator.attrgetter('start'))
    genome_list = [genome_list_1, genome_list_2]
    return genome_list


def indel_search(in_list, min_len, max_len, comm_in_list, min_distance, max_distance):
    indel_list = []
    for frg_n in range(len(in_list) - 1):
        dist_list1 = [
            sizer(in_list[frg_n].start, in_list[frg_n + 1].start),
            sizer(in_list[frg_n].end, in_list[frg_n + 1].start),
            sizer(in_list[frg_n].start, in_list[frg_n + 1].end),
            sizer(in_list[frg_n].end, in_list[frg_n + 1].end)]
        dist_list1.sort(key=operator.attrgetter('size'))
        if dist_list1[0].size > 0 and dist_list1[1].size > 0 and dist_list1[2].size > 0 and dist_list1[3].size > 0:
            if max_len > dist_list1[0].size > min_len:
                dist_list2 = [
                    sizer(comm_in_list[frg_n].start, comm_in_list[frg_n + 1].start),
                    sizer(comm_in_list[frg_n].end, comm_in_list[frg_n + 1].start),
                    sizer(comm_in_list[frg_n].start, comm_in_list[frg_n + 1].end),
                    sizer(comm_in_list[frg_n].end, comm_in_list[frg_n + 1].end)]
                dist_list2.sort(key=operator.attrgetter('size'))
                if dist_list2[0].size > 0 and dist_list2[1].size > 0 and dist_list2[2].size > 0 and dist_list2[3].size > 0:
                    if max_distance > abs(dist_list1[0].size - dist_list2[0].size) > min_distance:
                        print(dist_list1[0].size - dist_list2[0].size)
                        indel_list.append(indel(in_list[frg_n + 1].segment_id,
                                                in_list[frg_n].segment_id,
                                                in_list[frg_n + 1].start,
                                                in_list[frg_n + 1].end,
                                                in_list[frg_n].start,
                                                in_list[frg_n].end,
                                                dist_list1[0].a,
                                                dist_list1[0].b,
                                                in_list[frg_n + 1].genome,
                                                dist_list1[0].size))
    print('Indel search is done:', len(indel_list))
    return indel_list


def insertions_intrcept(indel_list1, indel_list2, wide):
    insertion_list = []
    number = 0
    for indels1 in indel_list1:
        insertion_list.append(insertion(indels1.start, indels1.end,
                                        indels1.genome, indels1.size, wide, indels1.genome))
        print('>', number)
        print(insertion_list[-1].sequence)
        number += 1
    for indels2 in indel_list2:
        insertion_list.append(insertion(indels2.start, indels2.end,
                                        indels2.genome, indels2.size, wide, indels2.genome))
        print('>', number)
        print(insertion_list[-1].sequence)
        number += 1
    print('Insertion search is done:', len(insertion_list))
    return insertion_list


def deletion_position(indel_list, index_list):
    deletion_list = []
    for indels1 in indel_list:
        ID = []
        for indels2 in index_list:
            if indels1.genome == indels2.genome:
                print('Warning: wrong genome index list!')
            else:
                if indels1.id1 == indels2.segment_id:
                    ID.append(indels2.start)
                    ID.append(indels2.end)
                elif indels1.id2 == indels2.segment_id:
                    ID.append(indels2.start)
                    ID.append(indels2.end)
        ID.sort()
        start_id = ID[1]
        deletion_list.append(deletion(indels1.start, indels1.end, indels1.genome, indels1.size, start_id))
    print('Deletion search is done:', len(deletion_list))
    return deletion_list


def main(first_align_path, second_align_path, first_genome, second_genome, third_genome, score_filtering_level,
         min_len, max_len, min_distance, max_distance, wide):
    a1 = similarity_filter(first_align_path, score_filtering_level)
    a2 = similarity_filter(second_align_path, score_filtering_level)
    b1 = index(a1, first_genome, second_genome)
    b2 = index(a2, first_genome, third_genome)
    c11 = indel_search(b1[0], min_len, max_len, b1[1], min_distance, max_distance)
    c21 = indel_search(b1[1], min_len, max_len, b1[0], min_distance, max_distance)
    c12 = indel_search(b2[0], min_len, max_len, b2[1], min_distance, max_distance)
    c22 = indel_search(b2[1], min_len, max_len, b2[0], min_distance, max_distance)

    d = insertions_intrcept(c11, c12, wide)
    f = deletion_position(c21, b1[0])
    g = deletion_position(c22, b2[0])

    dx, fx, gx = [], [], []
    for i in d:
        x = [i.start, i.end]
        dx.append(x)
    for i in f:
        x = [i.start, i.end]
        fx.append(x)
    for i in g:
        x = [i.start, i.end]
        gx.append(x)

    plt.subplot()
    plt.scatter(0, 1)
    plt.scatter(b1[0][-1].end, 1)
    plt.text(0, 1.002, 'START')
    plt.text(b1[0][-1].end, 1.002, 'END')
    for i in dx:
        plt.plot(i, [1, 1])
    for i in fx:
        plt.plot(i, [0.999, 0.999])
    for i in gx:
        plt.plot(i, [0.998, 0.998])
    plt.show()


main('C:/Theileria/Alignments/Annulata_Orientalis_Full_Alignment.xml',
     'C:/Theileria/Alignments/Annulata_Parva_Full_Alignment.xml',
     'C:/Theileria/Alignments/THEILERIA_ANNULATA_1_CHR.txt',
     'C:/Theileria/Alignments/THEILERIA_ORIENTALIS_1_CHR.txt',
     'C:/Theileria/Alignments/THEILERIA_PARVA_1_CHR.txt',
     200, 100, 120, 30, 100, 200)
