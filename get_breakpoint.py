import os
import logging
import distutils.spawn
from subprocess import Popen, PIPE, STDOUT
from argparse import ArgumentParser

# Check paths and executables
def check_exe(path, exe):
    if path is None:
        if distutils.spawn.find_executable(exe):
            return exe
        else:
            logging.critical("Error: %s executable is not in PATH" % exe)
            raise Exception("Error: %s executable is not in PATH" % exe)
    else:
        if distutils.spawn.find_executable(path):
            return path
        else:
            logging.critical("Error: %s path do not exist" % path)
            raise Exception("Error: %s path do not exist" % path)


def log_subprocess(out):
    for line in iter(out.readline, ''):
        if line != '\n':
            logging.debug(line.strip())


# Create blast index
def blast_index(ref, mdb):
    process = Popen([mdb, '-in', ref, '-input_type', 'fasta', '-dbtype', 'nucl'], universal_newlines=True,
                    stdout=PIPE, stderr=STDOUT)
    with process.stdout:
        log_subprocess(process.stdout)
    # Process finished
    exitcode = process.wait()
    if exitcode != 0:
        logging.critical("Error: makeblastdb failed")
        raise Exception("Error: makeblastdb failed, see log")


# Create FMD index
def fmd_index(ref, hsb):
    process = Popen([hsb, 'index', ref], universal_newlines=True, stdout=PIPE, stderr=STDOUT)
    with process.stdout:
        log_subprocess(process.stdout)
    exitcode = process.wait()
    if exitcode != 0:
        logging.critical("Error: hs-blastn index failed")
        raise Exception("Error: hs-blastn index failed, see log")


# Create windowmasker counts
def window_counts(ref, counts_path, wmk):
    process = Popen([wmk, '-in', ref, '-infmt', 'blastdb', '-mk_counts', '-out', counts_path],
                    universal_newlines=True,
                    stdout=PIPE, stderr=STDOUT)
    with process.stdout:
        log_subprocess(process.stdout)
    # Process finished
    exitcode = process.wait()
    if exitcode != 0:
        logging.critical("Error: windowmasker counts failed")
        raise Exception("Error: windowmasker counts failed, see log")


# Create counts binary
def window_binary(counts_path, obinary_path, wmk):
    process = Popen([wmk, '-in', counts_path, '-sformat', 'obinary', '-out', obinary_path, '-convert'],
                    universal_newlines=True, stdout=PIPE, stderr=STDOUT)
    with process.stdout:
        log_subprocess(process.stdout)
    exitcode = process.wait()
    if exitcode != 0:
        logging.critical("Error: windowmasker obinary failed")
        raise Exception("Error: windowmasker obinary failed, see log")


# Check presence of indexes
def check_index(ref, wk_dir, ref_name):
    counts_path = os.path.join(wk_dir, ref_name + '.counts')
    obinary_path = counts_path + '.obinary'
    # Only check if index files exist, they might be empty
    if os.path.isfile(ref + '.nsq') and os.path.isfile(ref + '.nhr') and os.path.isfile(ref + '.nin') and os.path.isfile(
            counts_path) and os.path.isfile(obinary_path) and os.path.isfile(ref + '.bwt') and os.path.isfile(ref + '.header') \
            and os.path.isfile(ref + '.sa') and os.path.isfile(ref + '.sequence'):
        return True
    else:
        return False


# Make indexes
def make_index(mode, ref, wk_dir, ref_name, mdb, wmk, hsb):
    counts_path = os.path.join(wk_dir, ref_name + '.counts')
    obinary_path = counts_path + '.obinary'
    if mode=='True':
        logging.info('Building blast database')
        blast_index(ref, mdb)
        logging.info('Windowmasker counts')
        window_counts(ref, counts_path, wmk)
        logging.info('Windowmasker binary')
        window_binary(counts_path, obinary_path, wmk)
        logging.info('Make FMD-index')
        fmd_index(ref, hsb)
    else:
        if os.path.isfile(ref + '.nsq') and os.path.isfile(ref + '.nhr') and os.path.isfile(ref + '.nin'):
            # Only check if index files exist, they might be empty
            logging.debug("Make blast index skipped")
        else:
            logging.info('Building blast database')
            blast_index(ref, mdb)
        if os.path.isfile(counts_path):  # Only check if files exist, they might be empty
            logging.debug("Windowmasker counts skipped")
        else:
            logging.info('Windowmasker counts')
            window_counts(ref, counts_path, wmk)
        if os.path.isfile(obinary_path):  # Only check if files exist, they might be empty
            logging.debug("Windowmasker obinary skipped")
        else:
            logging.info('Windowmasker binary')
            window_binary(counts_path, obinary_path, wmk)
        if os.path.isfile(ref + '.bwt') and os.path.isfile(ref + '.header') and os.path.isfile(ref + '.sa') and os.path.isfile(
                ref + '.sequence'):  # Only check if index files exist, they might be empty
            logging.debug("Make hs-blastn FMD-index skipped")
        else:
            logging.info('Make FMD-index')
            fmd_index(ref, hsb)


# Minimap2 alignment
def align_mm(ref, read, wk_dir, read_name, ref_name, threads, mm, data_type, st):
    out_path = os.path.join(wk_dir, '%s-%s-mm.bam' % (read_name, ref_name))
    if data_type == 'ont':
        x = 'map-ont'
    elif data_type == 'pacbio':
        x = 'map-pb'
    else:
        x = ''
    process = Popen([mm + ' -t ' + str(threads) + ' -ax ' + x + ' ' + ref + ' ' + read + ' | ' + st +
                     ' view -Sb - -o ' + out_path],
                    universal_newlines=True, stdout=PIPE, stderr=STDOUT, shell=True, executable='/bin/bash')
    with process.stdout:
        log_subprocess(process.stdout)
    exitcode = process.wait()
    if exitcode != 0:
        logging.critical("Error: minimap2 alignment failed")
        raise Exception("Error: minimap2 alignment failed, see log")
    return ['minimap2 -t ' + str(threads) + ' -ax ' + x + ' ' + ref + ' ' + read + ' | samtools view -Sb - -o ' + out_path,
            out_path]


# HS-BLASTN alignment
def align_hsb(ref, wk_dir, ref_name, threads, hsb,read_path):
    obinary_path = os.path.join(wk_dir, str(ref_name) + '.counts.obinary')
    out_path = os.path.join(wk_dir, 'temp-%s-blast.tab' % ref_name)
    process = Popen([hsb + ' align -db ' + ref + ' -window_masker_db ' + obinary_path + ' -query ' + read_path + ' -out '
                     + out_path + ' -outfmt 6 -num_threads ' + str(threads) + ' -max_target_seqs 3 -gapopen 0 '
                     '-gapextend 4 -penalty -3 -reward 2'],
                    universal_newlines=True, stdout=PIPE, stderr=STDOUT, shell=True, executable='/bin/bash')
    # cmd = process.args[0]
    with process.stdout:
        log_subprocess(process.stdout)
    exitcode = process.wait()
    if exitcode != 0:
        logging.critical("Error: hs-blastn alignment failed")
        raise Exception("Error: hs-blastn alignment failed, see log")
    os.remove(read_path)
    return ["hs-blastn align -db " + ref + " -window_masker_db obinary_path -query " + read_path + " -out " + out_path +
            " -outfmt 6 -num_threads " + str(threads) + " -max_target_seqs 3 -gapopen 0 -gapextend 4 -penalty -3 -reward 2",
            out_path]

def get_args():
    parser = ArgumentParser(description="")
    parser.add_argument("--wk_dir", type=str,help="work dir")
    parser.add_argument('--ref',type=str,help=" the file of ref")
    parser.add_argument('--mode',type=str,help="mode",default=True)
    parser.add_argument("--ref_name", type=str,help="the name of ref ")
    parser.add_argument("--mdb", type=str,help="the exe of mdb",default="makeblastdb")
    parser.add_argument("--wmk", type=str,help="exe of wmk",default='windowmasker')
    parser.add_argument("--hsb", type=str,help="the exe of hsb",default="hs-blastn")
    parser.add_argument("--read_path", type=str,help="the read of path")
    parser.add_argument("--t", type=str,help="the number of threads",default='4')
    args=parser.parse_args()
    return args


def main():
    args=get_args()
    wk_dir=args.wk_dir
    ref=args.ref
    mode=args.mode
    ref_name=args.ref_name
    mdb=args.mdb
    wmk=args.wmk
    hsb=args.hsb
    threads=args.t
    read_path=args.read_path
    if os.path.exists(wk_dir):
        pass
    else:
        os.makedirs(wk_dir)
    make_index(mode, ref, wk_dir, ref_name, mdb, wmk, hsb)
    align_hsb(ref, wk_dir, ref_name, threads, hsb,read_path)

if __name__=='__main__':
    main()
