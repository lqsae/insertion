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


import re
import pysam
from argparse import ArgumentParser




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


def get_query_reference_pos(cigar_tuple,reference_start):
    '''get the end and start of query and referencd'''
    sum_M=0
    sum_D=0
    sum_I=0
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
    dict_cigar_map={'M':0,'I':1,'D':2,'N':3,'S':4,'H':5}
    tags=sa_tag.strip(';').split(';')
    sub_optimal=[]
    chr_list=[]
    for j in tags:
        chr,pos,strand,cigar,map_qualitie,*others=j.split(',')
        if int(map_qualitie)>=map_q:
            str_list=re.findall(r'[A-Z]{1}',cigar)
            cigar_tuple=[(int(re.search('(\d+)'+i,cigar).group(1)),dict_cigar_map[i]) for i in str_list]
            list_pos=get_query_reference_pos(cigar_tuple,int(pos))
            total_list=[chr,strand]+list_pos
            total_list=[str(i) for i in total_list]
            format_line= '_'.join(total_list)
            sub_optimal.append(format_line)
            chr_list.append(chr)
    return  sub_optimal,chr_list
    

def parse_bamfile(bamfile,outfile_format,outfile_query,insertion,map_q):
    ''' parse the bamfile'''
    outfile_format=open(outfile_format,'w+')
    outfile_query=open(outfile_query,'w+')
    samfile=pysam.AlignmentFile(bamfile)
    for read in samfile:
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
                # print(all_format_line) 
                # # print(all_optimal)


def get_args():
    parser = ArgumentParser(description="get the map info")
    parser.add_argument("--i", type=str,help="the bam file")
    parser.add_argument("--info", type=str,help="the out file of format data ")
    parser.add_argument("--query", type=str,help="the out file of query sequence file")
    parser.add_argument("--ins", type=str,help="the name of insertion")
    parser.add_argument('--mapq',type=int,help="the mapping_quality of the read")
    args=parser.parse_args()
    return args


def main():
    args=get_args()
    bamfile=args.i
    outfile_format=args.info
    outfile_query=args.query
    insertion=args.ins
    map_q=args.mapq
    parse_bamfile(bamfile,outfile_format,outfile_query,insertion,map_q)


if __name__=='__main__':
    main()
