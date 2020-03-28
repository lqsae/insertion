from cluster import Cluster
from argparse import ArgumentParser

def get_insetion_sequece(read_sequence):
    dict_read_sequence={}
    with open(read_sequence) as f:
        for line in f:
            line=line.strip()
            tag=line.split('\t')
            name=tag[0]
            value=tag[1]
            if name not in dict_read_sequence.keys():
                dict_read_sequence[name]=value
    return dict_read_sequence


def extract_sequence(cluster_list,dict_read_sequence):
    for i in cluster_list:
        tag=i.split('_')
        qname=tag[0]
        left=int(tag[1])
        right=int(tag[2])
        sequence=dict_read_sequence[qname]
        left_sequence=sequence[0:left]
        right_sequence=sequence[right:]

        return qname,left_sequence,right_sequence

def get_args():
    parser = ArgumentParser(description="get the map info")
    parser.add_argument("--i", type=str,help="the bam file")
    parser.add_argument('--flag',type=str,help=" true or false")
    parser.add_argument('--k',type=int,help=" optimal k")
    parser.add_argument("--out1", type=str,help="the picture of k_selected ")
    parser.add_argument("--out2", type=str,help="the picture of fit picture")
    parser.add_argument("--out3", type=str,help="the sequence  of selected sequence ")
    parser.add_argument("--ins", type=str,help="the name of insertion")
    parser.add_argument("--qs", type=str,help="the sequence  of query")
    args=parser.parse_args()
    return args
    

def main():
    args=get_args()
    infile=args.i
    ins=args.ins
    cluster=Cluster(infile,ins)
    if args.flag =='True':
        out1=args.out1
        cluster.cluster(out1)
    else:
        read_sequence=args.qs
        dict_read_sequence=get_insetion_sequece(read_sequence)
        out2=args.out2
        out3=args.out3
        file=open(out3,'w+')
        k=args.k
        cluster.plt_cluster(k,out2)
        dict_cluster=cluster.cluster_selected(k)
        for key,value in dict_cluster.items():
            qname,left_sequence,right_sequence=extract_sequence(value,dict_read_sequence)
            left='>'+'_'+str(key)+'_'+qname+'_'+'left'+'\n'+left_sequence+'\n'
            file.write(left)
            right='>'+'_'+str(key)+'_'+qname+'_'+'right'+'\n'+right_sequence+'\n'
            file.write(right)
            
if __name__=='__main__':
    main()
