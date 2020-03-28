#!user/bin/python3
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans

class Cluster:
    def __init__(self,infile,insertion):
        self.infile=infile
        self.insertion=insertion

    def parse_mapinfo(self):
        infile=self.infile
        insertion=self.insertion
        reference_start_list=[]
        reference_end_list=[]
        query_start_list=[]
        query_end_list=[]
        qname_list=[]
        with open(infile) as f:
            for line in f:
                line=line.strip()
                tags=line.split('\t')
                qname=tags[-1]
                for tag in tags[:-1]:
                    chr=tag.split('_')[0]
                    reference_start=tag.split('_')[2]
                    reference_end=tag.split('_')[3]
                    query_start=tag.split('_')[4]
                    query_end=tag.split('_')[5] 
                    if chr==insertion:
                        reference_start_list.append(reference_start)
                        reference_end_list.append(reference_end)
                        qname_list.append(qname)
                        query_start_list.append(query_start)
                        query_end_list.append(query_end)
            dict_data={'r_start':reference_start_list,'r_end':reference_end_list,'q_start':query_start_list,'q_end':query_end_list,'qname':qname_list}
            df=pd.DataFrame(dict_data)
            return df


    def get_feauture(self):
        df=self.parse_mapinfo()
        x=df[['r_start','r_end']].astype('int')
        return x


    def cluster(self,out1):
        df=self.parse_mapinfo()
        x=self.get_feauture()
        d=[]
        for i in range(1,11):    #k取值1~11，做kmeans聚类，看不同k值对应的簇内误差平方和
            km=KMeans(n_clusters=i,init='k-means++',n_init=10,max_iter=300,random_state=0)
            km.fit(x)
            d.append(km.inertia_)  #inertia簇内误差平方和
        plt.plot(range(1,11),d,marker='o')
        plt.xlabel('number of clusters')
        plt.ylabel('distortions')
        plt.savefig(out1)
    
    def plt_cluster(self,k,out2):
        km=KMeans(n_clusters=k,init='k-means++',n_init=10,max_iter=300,random_state=0)
        x=self.get_feauture()
        km.fit(x)
        x['labels']=km.labels_
        for i in set(x.labels):
            m=x[x.labels==i]
            plt.scatter(m.r_start,m.r_end,label=i)
            plt.legend()
        plt.savefig(out2)


    def cluster_selected(self,k):
        from collections import defaultdict
        dict_cluster=defaultdict(list)
        km=KMeans(n_clusters=k,init='k-means++',n_init=10,max_iter=300,random_state=0)
        df=self.parse_mapinfo()
        x=self.get_feauture()
        km.fit(x)
        df['labels']=km.labels_
        for i in range(len(df)):
            label=df.labels[i]
            q_start=df.q_start[i]
            q_end=df.q_end[i]
            q_name=df.qname[i]
            q_info=q_name+'_'+q_start+'_'+q_end
            dict_cluster[label].append(q_info)
        return dict_cluster
