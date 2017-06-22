import os

def runGetOrf(input_contigs, output, find_parameter):
    #translations of regions between STOP codons: default setting of 0
    #translations of regions between STOP codons: nucleotide setting is 2
    os.system('getorf -find '+ find_parameter + ' ' + input_contigs + ' ' + output)

runGetOrf('contigs.fa', 'test0Output', '0')
runGetOrf('contigs.fa', 'test2Output', '2')
