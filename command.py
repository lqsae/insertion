import multiprocessing as mp
import time
import subprocess
import sys
import os
import logging


def get_filelist(file):
    file_list=[i for i in  os.listdir(file) if i.endswith('bam')]
    return file_list

def log_subprocess(out):
    for line in iter(out.readline, ''):
        if line != '\n':
            logging.debug(line.strip())


#00.prepare for align

#samtools index fa
def bwa_index(file,cmd):
    '''build the bwa index'''
    process=subprocess.Popen([cmd ,'index',file],universal_newlines=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    with process.stdout:
        log_subprocess(process.stdout)
    # Process finished
    exitcode = process.wait()
    if exitcode != 0:
        logging.critical("Error: bwa index failed")
        raise Exception("Error: bwa index failed, see log")
#samtools index bam


def faidx(file,cmd):
    process=subprocess.Popen([cmd ,'index',file],universal_newlines=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    with process.stdout:
        log_subprocess(process.stdout)
    # Process finished
    exitcode = process.wait()
    if exitcode != 0:
        logging.critical("Error: samtools faidx failed")
        raise Exception("Error: samtools faidx failed, see log")


#02.align 


#minimap2
def minimap2(ref,type,outfile,sample,input_file,t):
    '''using minimap2 align pacbio or nanopore data '''
    input_fq=os.path.join(input_file,'*.fastq')
    out_bam=os.path.join(outfile,'{0}.bam'.format(sample))
    process=subprocess.Popen(['minimap2 ' + ' --MD ' + ' -ax '+ type + ' -t ' + str(t) + 
                              ' -R '+ '"@RG\tID:{0}\tSM:{1}"'.format(sample,sample) + ' ' + ref + ' ' + input_fq + 
                              ' | ' +  ' samtools ' + ' sort '+ ' -@ '+ str(t)+ ' -o '+ out_bam],
                              universal_newlines=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,shell=True)
    with process.stdout:
        log_subprocess(process.stdout)
    exitcode = process.wait()
    if exitcode != 0:
        logging.critical("Error: minimap2  failed")
        raise Exception("Error: minimap2 failed, see log")

    return ['minimap2 ' + ' --MD ' + ' -ax '+ type + ' -t ' + str(t) + 
            ' -R '+ '"@RG\tID:{0}\tSM:{1}"'.format(sample,sample) + ' ' + ref + ' ' + input_fq + 
            ' | ' +  ' samtools ' + ' sort '+ ' -@ '+ str(t)+ ' -o '+ out_bam]

#ngmlr
def ngmlr(cmd,ref,sample,outfile,input,type,t):
    input_fq=os.path.join(input,'*.fastq')
    out_bam=os.path.join(outfile,'{0}.bam'.format(sample))

    process=subprocess.Popen([' cat ' + input_fq + ' | ' + cmd+' --presets '+type+' -t '+str(t)+ ' -r ' + ref + 
        ' | ' + ' samtools sort -@ ' + str(t) + ' o '+ out_bam] ,universal_newlines=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,shell=True  
        )
    with process.stdout:
        log_subprocess(process.stdout)
    exitcode = process.wait()
    if exitcode != 0:
        logging.critical("Error: ngmlr failed")
        raise Exception("Error: ngmlr failed, see log")

    return [' zcat ' + input_fq + ' | '+cmd+' --presets '+type+' -t '+str(t)+ ' -r '+ref+ ' | '+' samtools sort -@ '+str(t)+' o '+out_bam]



#samtools index
def samtools_index(file):
    process=subprocess.Popen(['samtools' ,'index',file],universal_newlines=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    with process.stdout:
        log_subprocess(process.stdout)
    # Process finished
    exitcode = process.wait()
    if exitcode != 0:
        logging.critical("Error: samtools index failed")
        raise Exception("Error: samtools index failed, see log")


#03.call sv 


#snifflles call sv 
def snifflles(input, output, sample):
    process=subprocess.Popen(['sniffles -s 5 --max_num_splits 7 -d 1000 -t 8 -l 50 -q 20 -n 0 -r 2000 -z 0 --report_BND --genotype \
        --cluster --cluster_support 1 -f 0.0 --min_homo_af 0.8 --min_het_af 0.3 --skip_parameter_estimation --report_seq ' +' -m '+input+' -v '+output],
        universal_newlines=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,shell=True)
    with process.stdout:
        log_subprocess(process.stdout)
    exitcode = process.wait()
    if exitcode != 0:
        logging.critical("Error: ngmlr failed")
        raise Exception("Error: ngmlr failed, see log")

    return ['sniffles -s 5 --max_num_splits 7 -d 1000 -t 8 -l 50 -q 20 -n 0 -r 2000 -z 0 --report_BND --genotype \
        --cluster --cluster_support 1 -f 0.0 --min_homo_af 0.8 --min_het_af 0.3 --skip_parameter_estimation --report_seq ' +' -m '+input+' -v '+output]




# multiprocess
def multiprocess_process(file_list,t,cmd):
    p=mp.Pool(t)
    for j in file_list:
        p.apply_async(samtools_index,args=(j,cmd))
    print('waiting for all subprocess done...')
    p.close()
    p.join()
    print('all subprocess done')

# def main():
#     ref='/mnt/CCX/pipeline/Variation/X101SC19110429-Z01-shuidao-liuqingshan/index/TOS17.fa'
#     wkdir='/mnt/CCX/pipeline/Variation/X101SC19110429-Z01-shuidao-liuqingshan/01.bwa'
#     input_file='/mnt/CCX/pipeline/Variation/X101SC19110429-Z01-shuidao-liuqingshan/00.data/HE708-15'
#     sample='HE708-15'
#     outfile=os.path.join(wkdir,'01.test')
#     if not os.path.exists(outfile):
#         os.mkdir(outfile)
#     else:
#         pass
#     type='map-ont'
#     t=8
#     m=minimap2(ref,type,outfile,sample,input_file,t)
#     print(m)



def main():
    cmd='samtools'
    t=4
    file=sys.argv[1]
    file_list=get_filelist(file)
    start=time.time()
    print ('start at : {0}'.format(start))
    multiprocess_process(file_list,t,cmd)
    end=time.time()
    print ('end at : {0}'.format(end))


if __name__=='__main__':
    main()
