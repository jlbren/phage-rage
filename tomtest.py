import os
import shutil
import os
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import argparse
import csv
from Bio.SeqRecord import SeqRecord

# pauda-build <references-fasta> <pauda Index Output Directory>
# pauda-run <input file> <outputFileBlastx> <databaseDirectory>

def run_pauda(input_filename, index_directory, output_name, num_procs):

    # pauda build
    index_input = 'viral_aa.fasta' 
    pauda_build_command = 'pauda-build ' + index_input + ' ' + index_directory
    print(pauda_build_command)
    os.system(pauda_build_command)

    # pauda run 
    pauda_command = 'pauda-run ' + input_filename + ' ' + output_name + ' ' + index_directory
    print(pauda_command)
    os.system(pauda_command)

    #fix to move .pna file to the right place
    os.system('mv ' + output_name + '.pna blastOutput/')
    return True

# bin/lambda_indexer -d viral_aa.fasta
# bin/lambda -q /path/to/1k_long_reads.fna -d /path/to/uniprot_sprot.fasta.gz
def run_lambda(input_filename, index_directory, output_name, num_procs):

    #lambda index build
    index_input = 'viral_aa.fasta'
    lambda_indexer = 'lambda_indexer -d ' + index_input + ' -i ' + index_directory
    print(lambda_indexer)
    os.system(lambda_indexer)

    #lambda run
    lambda_command = 'lambda -p blastp -q ' + input_filename + ' -i ' + index_directory + ' -o ' + output_name + '.m8'
    print(lambda_command)
    os.system(lambda_command)
    return True


#diamond:
#./diamond makedb --in nr.faa -d databaseName
#./diamond blastp -d databaseName -q reads.fna -o matches.m8
def run_diamond(input_filename, index_directory, output_name, num_procs):

    #diamond index build
    index_input = 'viral_aa.fasta'
    diamond_indexer = 'diamond makedb --in ' + index_input + ' -d ' + index_directory
    print(diamond_indexer)
    os.system(diamond_indexer)

    #diamond run
    diamond_command = 'diamond blastp -d ' + index_directory + ' --threads ' + num_procs + ' -q ' + input_filename + ' -o ' + output_name
    print(diamond_command)
    os.system(diamond_command)

def run_blastp(input_filename,index_directory,output_name,num_procs):
    
    index_input = 'viral_aa.fasta'
    blastpIndexer = 'makeblastdb -in ' + index_input + ' -dbtype prot -out ' + index_directory
    os.system(blastpIndexer)    
    
    headers=['qseqid sseqid qcovs qstart qend pident length evalue bitscore']
    os.system('blastp -query ' + input_filename + ' -db ' + index_directory + ' -outfmt 10 -num_threads ' + num_procs)

def parse_blast(filename,headers):
    x=[]
    blast_results=open(filename,'r')
    rows=csv.DictReader(blast_results,headers,delimiter=',')
    for row in rows:
        x.append(row)
        print(x)
    blast_results.close()
    return x

#pauda needs to have .pna file moved by code in the function once it finishes running
#blastOutput is the name of the output folder to consolidate everything
#TODO add creating the blastOutput directory to the installation script
#run_pauda('test.faa', 'indexFolder/paudaIndexDirectory', 'paudaOutput', '1')
#run_lambda('test.faa', 'indexFolder/index.lambda', 'blastOutput/lambdaOutput', '1')
#run_diamond('test.faa', 'indexFolder/diamondIndexDirectory', 'blastOutput/diamondOutput', '1')
run_blastp('test.faa', 'indexFolder/blastpIndexDirectory', 'blastOutput/blastpOutput', '1')
