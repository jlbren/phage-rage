import os

# pauda-build <references-fasta> <pauda Index Output Directory>
# pauda-run <input file> <outputFileBlastx> <databaseDirectory>

def run_pauda(input_filename, index_directory, output_name):
    
    # pauda build
    index_input = 'viral_aa.fasta' 
    pauda_build_command = './pauda/bin/pauda-build ' + index_input + ' ' + index_directory
    print(pauda_build_command)
    os.system(pauda_build_command)

    # pauda run 
    pauda_command = 'pauda/bin/pauda-run ' + input_filename + ' /temp/paudaAssembly/' + output_name + ' ' + index_directory
    print(pauda_command)
    os.system(pauda_command)
    return True

# bin/lambda_indexer -d viral_aa.fasta
# bin/lambda -q /path/to/1k_long_reads.fna -d /path/to/uniprot_sprot.fasta.gz
def run_lambda(input_filename, index_directory):

    #lambda index build
    index_input = 'viral_aa.fasta'
    lambda_indexer = 'lambda-1.9.2-Linux-x86_64/bin/lambda_indexer -d ' + index_input + ' -i index.lambda'
    print(lambda_indexer)
    os.system(lambda_indexer)

    #lambda run
    lambda_command = 'lambda-1.9.2-Linux-x86_64/bin/lambda -q ' + input_filename + ' -i ' + index_directory
    print(lambda_command)
    os.system(lambda_command)
    return True


#diamond:
#./diamond makedb --in nr.faa -d databaseName
#./diamond blastp -d databaseName -q reads.fna -o matches.m8
def run_diamond(input_filename, index_directory, output_name):

    #diamond index build
    index_input = 'viral_aa.fasta'
    diamond_indexer = './diamond makedb --in ' + index_input + ' -d ' + index_directory
    print(diamond_indexer)
    os.system(diamond_indexer)

    #diamond run
    diamond_command = './diamond blastp -d ' + index_directory + ' -q ' + input_filename + ' -o ' + output_name
    print(diamond_command)
    os.system(diamond_command)

run_pauda('contigs.fa', 'paudaIndexDirectory', 'paudaBlastXOutput')
run_lambda('contigs.fa', 'index.lambda')
run_diamond('contigs.fa', 'diamondIndexDirectory', 'diamondTestOutput')
