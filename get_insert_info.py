#!/usr/bin/env python3
'''
 * All rights Reserved, Designed By NovoGene
 * @Title:  insert.py
 * @Package:
 * @Description: Control rCANID pipeline
 * @author: Qingshan Liu
 * @date: March 03 2020
 * @version V1.0.1
'''
import time
import os

import re
import pysam
import subprocess
import multiprocessing as mp
from argparse import ArgumentParser


def get_references(bamfile):
    ''' get all the contigs'''
    samfile=pysam.AlignmentFile(bamfile)
    references=samfile.references
    return references


def get_optimal_reference_query_pos(read):
    '''get the end and start of the query and start'''
    reference_start=read.reference_start
    reference_end=read.reference_end
    query_start=read.query_alignment_start
    query_end=read.query_alignment_end
    return [reference_start,reference_end,query_start,query_end]


def get_revese_complenment_sequence(read_sequence):
    '''get the reeverse complemented of sequence'''
    revComp = str.maketrans("ATCGNatcgn","TAGCNtagcn")
    sequence=read_sequence.translate(revComp)[::-1]
    return sequence

def get_query_reference_pos(cigar,reference_start):
    '''get the end and start of query and referencd'''
    dict_cigar_map={'M':0,'I':1,'D':2,'N':3,'S':4,'H':5}
    sum_M=0
    sum_D=0
    sum_I=0
    cigar_list=re.findall(r'[A-Z]{1}',cigar)
    set_cigar_list=list(set(cigar_list))
    set_cigar_list.sort(key=cigar_list.index)
    cigar_tuple=[]
    for i in set_cigar_list:
        match=re.findall('(\d+)'+i,cigar)
        tuple_i=[(int(s),dict_cigar_map[i]) for s in match]
        cigar_tuple.extend(tuple_i)

    for i,j in cigar_tuple:
        if j==0:
            sum_M+=i
        elif j==1:
            sum_I+=i
        elif j==2:
            sum_D+=i

    bias_reference=sum_D+sum_M
    bias_query=sum_M+sum_I
    reference_end=reference_start+bias_reference

    if cigar_tuple[0][1]==4:
        query_start=cigar_tuple[0][0]
    else:
        query_start=0
    query_end=query_start+bias_query

    optimal_pos=[reference_start,reference_end,query_start,query_end]
    optimal_pos=[str(i) for i in optimal_pos]
    return optimal_pos


def get_Sub_optimal(sa_tag,map_q):
    '''get the SA_tag cigar information'''
    tags=sa_tag.strip(';').split(';')
    sub_optimal=[]
    chr_list=[]
    for j in tags:
        chr,pos,strand,cigar,map_qualitie,*others=j.split(',')
        if int(map_qualitie)>=map_q:
            list_pos=get_query_reference_pos(cigar,int(pos))
            total_list=[chr,strand]+list_pos
            total_list=[str(i) for i in total_list]
            format_line= '_'.join(total_list)
            sub_optimal.append(format_line)
            chr_list.append(chr)
    return  sub_optimal,chr_list


def parse_bamfile(bamfile,contig,wkdir,sample,insertion,map_q):
    ''' parse the bamfile'''
    mapinfo_file='{0}.{1}.{2}.temp'.format(sample,contig,'mapinfo')
    query_sequence_file='{0}.{1}.{2}.temp'.format(sample,contig,'query_sequence')

    outfile_format=os.path.join(wkdir,mapinfo_file)
    outfile_query=os.path.join(wkdir,query_sequence_file)


    outfile_format=open(outfile_format,'w')
    outfile_query=open(outfile_query,'w')

    samfile=pysam.AlignmentFile(bamfile)
    if samfile.has_index():
        reads=samfile.fetch(contig=contig)
        for read in reads:
            try:
                cigartuples=read.cigartuples
                ciga_s=cigartuples[0][0]
            except:
                pass
            if ciga_s ==4 and read.has_tag('SA'):#extract the soft clip  and has sa_tag sequence
                sa_tag=read.get_tag('SA')
                reference_name=read.reference_name
                query_sequence=read.query_sequence
                query_name=read.query_name
                mapping_quality=int(read.mapping_quality)
                optimal_pos=get_optimal_reference_query_pos(read)
                all_suboptimal,chr_list=get_Sub_optimal(sa_tag,map_q)
                all_chr=chr_list+[reference_name]
                if (insertion in all_chr) and (len(set(all_chr))>=2) and(mapping_quality>=map_q):#the query sequece map to reference and insertion
                    if read.is_reverse:
                        strand='-'
                    else :
                        strand='+'
                    h=[reference_name,strand]+optimal_pos
                    all_optimal='_'.join([str(i) for i in h])
                    all_format_data=[all_optimal]+all_suboptimal+[query_name]

                    query_info=[query_name,query_sequence]
                    query_info_line='\t'.join(query_info)+'\n'
                    all_format_data=[str(i) for i in all_format_data]
                    all_format_line='\t'.join(all_format_data)+'\n'

                    outfile_query.write(query_info_line)
                    outfile_format.write(all_format_line)
    outfile_format.close()
    outfile_query.close()

                    # print(all_format_line)
                    # print(all_optimal)



def multiprocess_process(bamfile,t,contigs,wkdir,insertion,map_q,sample):
    ''' Split the  chromosome and call'''
    p=mp.Pool(t)
    for contig  in contigs:
        p.apply_async(parse_bamfile,args=(bamfile,contig,wkdir,sample,insertion,map_q))
    print('waiting for all subprocess done...')
    p.close()
    p.join()
    print('all subprocess done')



def cat_all_file_run(wkdir,sample):
    all_map_info_file=os.path.join(wkdir,sample+'*.mapinfo.temp')
    map_info_file=os.path.join(wkdir,sample+'.mapinfo')
    all_query_sequence_file=os.path.join(wkdir,sample+'*.query_sequence.temp')
    query_sequence_file=os.path.join(wkdir,sample+'.query_sequence')
    result1=subprocess.Popen(['cat ' + all_map_info_file + ' > ' + map_info_file] ,shell=True)
    exitcode = result1.wait()
    if exitcode != 0:
        print ("Error: cat all_map_info_file failed ")

    result2=subprocess.Popen(['cat ' + all_query_sequence_file + ' > ' +  query_sequence_file] ,shell=True)
    exitcode = result2.wait()
    if exitcode != 0:
        print ("Error: cat all_query_sequence_file failed ")
    if os.path.exists(map_info_file) and os.path.exists(query_sequence_file) :
        result3=subprocess.Popen(['rm ' + all_map_info_file],shell=True)
        exitcode = result3.wait()
        if exitcode != 0:
            print ("Error: rm  all_map_info_file failed ")
        result4=subprocess.Popen(['rm ' + all_query_sequence_file ],shell=True)
        exitcode = result4.wait()
        if exitcode != 0:
            print ("Error: rm  all_query_sequence_file failed ")

def get_args():
    parser = ArgumentParser(description="get the map info")
    parser.add_argument("--i", type=str,help="the bam file")
    parser.add_argument("--wkdir", type=str,help="the out put file")
    parser.add_argument("--ins", type=str,help="the name of insertion")
    parser.add_argument("--sample", type=str,help="the name of sample")
    parser.add_argument("--t", type=int,help="the number of threads")
    parser.add_argument('--mapq',type=int,help="the mapping_quality of the read")
    args=parser.parse_args()
    return args


def main():
    args=get_args()
    bamfile=args.i
    wkdir=args.wkdir
    insertion=args.ins
    map_q=args.mapq
    sample=args.sample
    t=args.t
    print ('start at : ' + time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
    contigs=get_references(bamfile)
    multiprocess_process(bamfile,t,contigs,wkdir,insertion,map_q,sample)
    cat_all_file_run(wkdir,sample)
    print ('end at : ' + time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))

if __name__=='__main__':
    main()

