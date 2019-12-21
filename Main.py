from Bio.Blast.Applications import NcbiblastnCommandline
from Bio import SearchIO
from Bio.Seq import Seq
from Bio.Blast import NCBIXML
from Bio import Entrez
import operator
import os


class gomSegment():

    def __init__(self, sid, start, end, genome):
        self.sid = sid
        self.start = start
        self.end = end
        self.genome = genome


class indel():

    def __init__(self, id1, id2, start1, end1, start2, end2, genome, size):
        self.id1 = id1
        self.id2 = id2
        self.start1 = start1
        self.end1 = end1
        self.start2 = start2
        self.end2 = end2
        self.genome = genome
        self.size = size


def keyfunc(item):
    return item[1]


def similarity_filter(inp_file, level):
    blast_qresult = SearchIO.read(inp_file, "blast-xml")
    table = []
    for hsp in blast_qresult:
        for h in hsp:
            if float(h.ident_num * 100 / len(h.query)) > int(level):
                table.append(h)
    return table


def index(in_list, genome1, genome2):
    genome_list_1 = []
    genome_list_2 = []
    number = 0
    for hits in in_list:

        genome_list_1.append(gomSegment(number, int(hits.query_range[0]), int(hits.query_range[1]), genome1))
        genome_list_2.append(gomSegment(number, int(hits.hit_range[0]), int(hits.hit_range[1]), genome2))
        number += 1

    genome_list_1.sort(key=operator.attrgetter('start'))
    genome_list_2.sort(key=operator.attrgetter('start'))
    genome_list = [genome_list_1, genome_list_2]
    return genome_list


def indel_search(in_list):
    indel_list = []
    for frg_n in range(len(in_list) - 2):
        for n in range(1, len(in_list) - 1 - frg_n):
            l1 = in_list[frg_n + n].start - in_list[frg_n].start
            l2 = in_list[frg_n + n].end - in_list[frg_n].start
            l3 = in_list[frg_n + n].start - in_list[frg_n].end
            l4 = in_list[frg_n + n].end - in_list[frg_n].end
            if l1 > 0 and l2 > 0 and l3 > 0 and l4 > 0:
                if 1500 > min(l1, l2, l3, l4) > 100:
                    indel_list.append(indel(in_list[frg_n].sid, in_list[frg_n + n].sid, in_list[frg_n].start,
                                            in_list[frg_n].end, in_list[frg_n + n].start, in_list[frg_n + n].end,
                                            in_list[frg_n].genome, min(l1, l2, l3, l4)))
    return indel_list


def alignfunc(in_list, out_f):
    in_file = in_list[0]
    db_file = in_list[1]
    cmd = r'C:/Theileria_seq/ncbi-blast-2.9.0+/bin/makeblastdb -in %s -dbtype nucl' % db_file
    os.system(cmd)
    proga = r'C:/Theileria_seq/ncbi-blast-2.9.0+/bin/blastn'
    blast = NcbiblastnCommandline(proga, query=in_file, db=db_file, out=out_f, outfmt=5,
                                  word_size=12, evalue=0.001)
    stdout, stderr = blast()


def search_db(sequences, out_file):
    with open(out_file, 'w') as f:
        for seq in sequences:
            handle = Entrez.efetch(db="gene", id=str(seq), rettype="gb", retmode="text")
            for lines in handle:
                f.writelines(lines)


def uniq(in_list_1, in_list_2):
    for i in in_list_1[1].values():
        for j in in_list_2[1].values():
            start_1 = int(i.strip('_'))
            start_2 = int(j.strip('_'))
            end_1 = int(i.strip('_'))
            end_2 = int(j.strip('_'))
            print(start_1, start_2, end_1, end_2)


def writer(genome_file_1, genome_file_2, result, name_1, name_2, out_file):
    genome_1 = Seq(genome_file_1)
    genome_2 = Seq(genome_file_2)

    query_start_list, sbjct_start_list, query_end_list, sbjct_end_list = [], [], [], []
    result_handle = open(result, 'r')
    blast_record = NCBIXML.read(result_handle)
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

    def space_counter(start_list, end_list):
        space_list = []
        for i in range(len(start_list)):
            space = start_list[i] - end_list[i]
            space_list.append(space)
        space_list.remove(space_list[0])
        space_list.append(0)
        return space_list

    space_list_1 = space_counter(query_start_list, query_end_list)
    space_list_2 = space_counter(sbjct_start_list, sbjct_end_list)

    uniq_list = []
    for i in range(len(space_list_1)):
        if space_list_1[i] > space_list_2[i]:
            uniq_list.append(space_list_1[i])
        else:
            uniq_list.append(space_list_2[i])

    with open(out_file, 'w') as f:
        f.write(name_1 + '\n')
        for i in genome_1[query_start_list[0]: query_end_list[-1]]:
            f.write(i)
        f.write('\n')

        n = 0
        for j in query_start_list:
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    if hsp.query_start == j:
                        f.write(hsp.query)
                        f.write(' ' * 15)
                        f.write('(' + str(space_list_1[n]) + ')')
                        t = 15 + (len(str(uniq_list[n])) - len(str(space_list_1[n])))
                        f.write(' ' * t)
                        n += 1

        f.write('\n')
        k = 0
        for j in query_start_list:
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    if hsp.query_start == j:
                        f.write(hsp.match)
                        f.write(' ' * 15)
                        e = 17 + len(str(uniq_list[k]))
                        f.write(' ' * e)
                        k += 1
        f.write('\n')
        m = 0
        for j in query_start_list:
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    if hsp.query_start == j:
                        f.write(hsp.sbjct)
                        f.write(' ' * 15)
                        f.write('(' + str(space_list_2[m]) + ')')
                        u = 15 + (len(str(uniq_list[m])) - len(str(space_list_2[m])))
                        f.write(' ' * u)
                        m += 1
        f.write('\n')
        f.write('\n')
        f.write(name_2 + '\n')
        for i in genome_2[sbjct_start_list[0]: sbjct_end_list[-1]]:
            f.write(i)
