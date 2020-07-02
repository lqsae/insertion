#!/usr/bin/env python3
'''
 * All rights Reserved, Designed By Qingshan Liu
 * @Title:  insert.py
 * @Package:
 * @Description: Control get insertion position pipeline
 * @author: Qingshan Liu
 * @date: March 03 2020
 * @version V1.0.1
'''
import time
import re
import sys
import logging

import pysam
from collections import defaultdict
from multiprocessing import Pool, Manager
from argparse import ArgumentParser
logging.basicConfig(level=logging.INFO, format='%(asctime)s :: %(levelname)s :: %(message)s', filename='insertion.log')


def get_references(bamfile):
    ''' get all the contigs'''
    samfile = pysam.AlignmentFile(bamfile)
    references = samfile.references
    return references


def get_optimal_reference_query_pos(read):
    '''get the end and start of the query and start'''
    reference_start = read.reference_start
    reference_end = read.reference_end
    query_start = read.query_alignment_start
    query_end = read.query_alignment_end
    map_q = int(read.mapping_quality)
    reference_name = read.reference_name
    if read.is_reverse:
        strand = '-'
    else:
        strand = '+'
    optimal_info_list = [reference_name, strand, map_q,
                         reference_start, reference_end, query_start, query_end]
    return optimal_info_list


def get_revese_complenment_sequence(read_sequence):
    '''get the reeverse complemented of sequence'''
    revComp = str.maketrans("ATCGNatcgn", "TAGCNtagcn")
    sequence = read_sequence.translate(revComp)[::-1]
    return sequence


def get_query_reference_pos(cigar, reference_start):
    '''get the end and start of query and referencd in SA tag'''
    dict_cigar_map = {'M': 0, 'I': 1, 'D': 2, 'N': 3, 'S': 4, 'H': 5}
    sum_M = 0
    sum_D = 0
    sum_I = 0
    cigar_list = re.findall(r'[A-Z]{1}', cigar)
    set_cigar_list = list(set(cigar_list))
    set_cigar_list.sort(key=cigar_list.index)
    cigar_tuple = []
    for i in set_cigar_list:
        match = re.findall('(\d+)' + i, cigar)
        tuple_i = [(int(s), dict_cigar_map[i]) for s in match]
        cigar_tuple.extend(tuple_i)
    for i, j in cigar_tuple:
        if j == 0:
            sum_M += i
        elif j == 1:
            sum_I += i
        elif j == 2:
            sum_D += i
    bias_reference = sum_D + sum_M
    bias_query = sum_M + sum_I
    reference_end = reference_start + bias_reference
    if cigar_tuple[0][1] == 4 or cigar_tuple[0][1] == 5:
        query_start = cigar_tuple[0][0] - 1
    else:
        query_start = 0
    query_end = query_start + bias_query
    optimal_pos = [reference_start, reference_end, query_start, query_end]
    optimal_pos = [str(i) for i in optimal_pos]
    return optimal_pos


def get_Sub_optimal(sa_tag):
    '''get the SA_tag cigar information'''
    tags = sa_tag.strip(';').split(';')
    sub_optimal_list = []
    for j in tags:
        chr, pos, strand, cigar, map_q, *others = j.split(',')
        list_pos = get_query_reference_pos(cigar, int(pos))
        total_list = [chr, strand, map_q] + list_pos
        sub_optimal_list.append(total_list)
    return sub_optimal_list


def parse_bamfile(bamfile, contig,q):
    ''' parse the bamfile'''
    dict_main = defaultdict(list)
    samfile = pysam.AlignmentFile(bamfile)
    if samfile.has_index():
        reads = samfile.fetch(contig=contig)
    else:
        logging.error('the bam has not index')
        sys.exit(0)
    for read in reads:
        query_name = read.query_name
        try:
            cigartuples = read.cigartuples
            ciga_s = cigartuples[0][0]
        except:
            pass
        if (ciga_s == 4 or ciga_s == 5) and read.has_tag('SA'):  # extract the soft clip  and has sa_tag sequence
            sa_tag = read.get_tag('SA')
            optimal_info_list = get_optimal_reference_query_pos(read)
            sub_optimal_list = get_Sub_optimal(sa_tag)
            sub_optimal_list.append(optimal_info_list)
            dict_main[query_name].extend(sub_optimal_list)
    q.put(dict_main)


def get_read_sequence(read):
    query_name = read.query_name
    query_sequence = read.query_sequence
    fa = '>{0}\n{2}'.format(query_name, query_sequence)
    return fa


def parse_dict_map(q, insertion, map_q, gap):
    ''' optimal_info_list = [reference_name, strand, map_q,
                         reference_start, reference_end, query_start, query_end]
                         parse the dict_info and get the breakpoint'''
    while True:
        if not q.empty():
            dict_main = q.get(True)
        else:
            break
    for key, value in dict_main.items():
        uniq_value = set([tuple(i) for i in value]) #把列表转化为元组去重
        map_info_list = sorted(uniq_value, key=lambda x: int(x[6]))  # 根据在染色体上的位置排序
        length = len(map_info_list)
        map_q_list = [int(i[2]) for i in uniq_value]
        chr_list = [i[0] for i in uniq_value]
        map_q_flag = all(i > map_q for i in map_q_list)
        if (insertion in chr_list) and (len(set(chr_list)) > 2) and map_q_flag:
            for i, j in enumerate(map_info_list):
                chr = j[0]
                ins_start = int(j[-2])
                ins_end = int(j[-1])
                if chr == insertion and i == 0:
                    right_start = int(map_info_list[i + 1][-2])
                    if right_start - gap <= ins_end <= right_start + gap:
                        break_point = int(map_info_list[i + 1][-4])
                        ref_info = [str(i) for i in map_info_list[i + 1]]
                        insert_info = [str(i) for i in map_info_list[i]]
                        print(
                            '\t'.join(ref_info) + '\t' + '\t'.join(insert_info) + '\t' + str(break_point) + '\t' + key)
                elif chr == insertion and i == length - 1:
                    left_end = int(map_info_list[i - 1][-1])
                    if left_end - gap <= ins_start <= left_end + gap:
                        break_point = int(map_info_list[i - 1][-3])
                        ref_info = [str(i) for i in map_info_list[i - 1]]
                        insert_info = [str(i) for i in map_info_list[i]]
                        print(
                            '\t'.join(ref_info) + '\t' + '\t'.join(insert_info) + '\t' + str(break_point) + '\t' + key)
                elif chr == insertion and (i not in [0, length - 1]):
                    right_start = int(map_info_list[i + 1][-2])
                    left_end = int(map_info_list[i - 1][-1])
                    if right_start - gap <= ins_end <= right_start + gap:
                        break_point = int(map_info_list[i + 1][-4])
                        ref_info = [str(i) for i in map_info_list[i + 1]]
                        insert_info = [str(i) for i in map_info_list[i]]
                        print(
                            '\t'.join(ref_info) + '\t' + '\t'.join(insert_info) + '\t' + str(break_point) + '\t' + key)
                    elif left_end - gap <= ins_start <= left_end + gap:
                        break_point = int(map_info_list[i - 1][-3])
                        ref_info = [str(i) for i in map_info_list[i - 1]]
                        insert_info = [str(i) for i in map_info_list[i]]
                        print(
                            '\t'.join(ref_info) + '\t' + '\t'.join(insert_info) + '\t' + str(break_point) + '\t' + key)


def get_args():
    parser = ArgumentParser(description="get the map info")
    parser.add_argument("--i", type=str, help="the bam file")
    parser.add_argument("--ins", type=str, help="the name of insertion")
    parser.add_argument("--t", type=int, help="the number of threads")
    parser.add_argument('--mq', type=int, help="the mapping_quality of the read")
    parser.add_argument('--gap', type=int, help="the mapping_quality of the read")
    args = parser.parse_args()
    return args


def main():
    args = get_args()
    bamfile = args.i
    insertion = args.ins
    map_q = args.mq
    t = args.t
    gap = args.gap
    q = Manager().Queue()
    contigs = get_references(bamfile)
    logging.info('start at : ' + time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
    p = Pool(t)
    for contig in contigs:
        p.apply_async(parse_bamfile, args=(bamfile, contig, q))
        p.apply_async(parse_dict_map, args=(q, insertion, map_q, gap))
    p.close()
    p.join()
    logging.info('end at : ' + time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))


if __name__ == '__main__':
    main()
