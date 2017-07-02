from glob import glob
from Bio import SeqIO
import os
import sys

#def make_protein_list(index_dir):
#    genomes = glob(index_dir)
#
#    for f in genomes:
#        gacc = os.path.splitext(f)[0]
#        
#        with open(f, 'r') as handle:
#            seqs = SeqIO.Parse(handle)
#        for seq in seqs:
#            # parse pacc, species

class Genome:
    def __init__(self, g_acc, tax):
        self.genome_acc = g_acc
        self.taxonomy = tax
        self.protein_acc_list = [] 
        self.protein_hit_list = []

    def set_proteins(self, plist):
        self.protein_acc_list = plist 
        self.protein_hit_list = [0] * len(plist)

class VParse:
    def __init__(self):
        self.genomes = []
        self.gbk_dir = None
        self.gbk_files = None
        self.faa_index = None

    def parse_index(self, gbk_dir):
        self.gbk_dir = gbk_dir
        self.gbk_files = self.parse_gbk_dir(gbk_dir)
        self.parse_gbk_files(self.gbk_files)

    def parse_gbk_files(self, gbk_file_list):
        output_handle = open('/home/jon/Projects/phage-rage/test_index.faa', 'a+')
        for f in gbk_file_list:
            input_handle = open(f, 'r')
            for record in SeqIO.parse(input_handle, 'genbank'):
                organism = Genome(record.name, record.annotations['taxonomy'])
                protein_list = []
                for feature in record.features:
                    if feature.type=="CDS" :
                        assert len(feature.qualifiers['translation']) == 1 # TODO ask about this 
                        protein_list.append(feature.qualifiers['protein_id'])
                        output_handle.write(">%s|ref|%s|%s\n%s\n" % # TODO ask if this is good 
                                            (
                                             feature.qualifiers['db_xref'][0].replace(':', '|'),
                                             feature.qualifiers['protein_id'][0],
                                             feature.qualifiers['product'][0],
                                             feature.qualifiers['translation'][0]
                                            ))
                organism.set_proteins(protein_list)
                self.genomes.append(organism)
                g = self.genomes[0]
                print('gacc: %s, tax: %s, proteins: %s hits: %s' % 
                        (g.genome_acc, str(g.taxonomy), str(g.protein_acc_list), str(g.protein_hit_list)))




    def parse_gbk_dir(self, gbk_dir):
        gbk_files = glob(os.path.join(gbk_dir, '*.gbk'))
        sub_dirs = glob(os.path.join(gbk_dir, '*'))
        for f in sub_dirs:
            if os.path.isdir(f):
                gbk_files.append(glob(os.path.join(f,'*.gbk')))
        if not gbk_files:
            raise ValueError('Improperly formatted index folder.\n'
                             'No GBK files found in: %s'
                             % gbk_dir)
    
        return gbk_files

if __name__ == '__main__':
    gbk_dir = sys.argv[1]

    vp = VParse()
    vp.parse_index(gbk_dir)
#gbk_filename = "NC_005213.gbk"
#
##faa_filename = "NC_005213_converted.faa"
#
#input_handle  = open(gbk_filename, "r")
#
#output_handle = open(faa_filename, "w")
#
#
#
#for seq_record in SeqIO.parse(input_handle, "genbank") :
#
#    print "Dealing with GenBank record %s" % seq_record.id
#
#    for seq_feature in seq_record.features :
#
#        if seq_feature.type=="CDS" :
#
#            assert len(seq_feature.qualifiers['translation'])==1
#
#            output_handle.write(">%s\n%s\n" % (
#
#                   seq_feature.qualifiers['locus_tag'][0],
#
#                   seq_feature.qualifiers['translation'][0]))
#
#
#
#output_handle.close()
#
#input_handle.close()
