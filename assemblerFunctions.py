''' Class managing velvet: instantiate this class so that command of velvet
 can match user's system settings.
 author: Thomas Hatzopoulos
 date: 2017/05/22
VelvetPackage fetches the reads to be assembled and calls velvet with the
appropriate parameters provided in the config file'''

import subprocess
import os



def call_velvet(input1, input2, output):
    '''Call velvet with parameters provided in contig file'''

    read_one = input1  # location and name of sequence 1
    read_two = input2  # location and name of sequence 2
    file_name = input1 + input2

    # "runs/" = "/runs/" + file_name
    # Output velvet system commands
    print('velvet hashing ' + file_name)

    # sets kmer, format, and read type
    #  velveth runs/ 31 -fastq -shortPaired seqs/140610-1_S1_L001_R1_001.fastq -fastq -shortPaired seqs/140610-1_S1_L001_R2_001.fastq 
    subprocess.call(['velveth', output, '31', '-fastq', '-shortPaired', read_one, read_two])

    print('velvetg running')
    os.system('velvetg ' + output + ' -exp_cov auto')

    # Move contig files, clean up unnecessary files

def call_spades(input1, input2, output):
    '''Call metaspades '''

    subprocess.call(['SPAdes-3.10.1-Linux/bin/spades.py', '--meta', '--pe1-1', input1, '--pe1-2', input2, '-o', output])

def call_megahit(input1, input2, output):
    '''Call Megahit '''
    subprocess.call(['./megahit/megahit', '-1', input1, '-2', input2, '-o', output])


call_velvet("data/trimmed_AW-9_S9_L001_R1_001.fastq", "data/trimmed_AW-9_S9_L001_R2_001.fastq", 'velvetOutput')
call_spades("data/trimmed_AW-9_S9_L001_R1_001.fastq", "data/trimmed_AW-9_S9_L001_R2_001.fastq", 'spadesOutput')
call_megahit("data/trimmed_AW-9_S9_L001_R1_001.fastq", "data/trimmed_AW-9_S9_L001_R2_001.fastq", 'megahitOutput')
