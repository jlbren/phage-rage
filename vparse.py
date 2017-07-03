from glob import glob
from Bio import SeqIO
import os
import sys

class Genome:
    def __init__(self, g_acc, species, tax):
        self.genome_acc = g_acc
        self.taxonomy = tax
        self.species = species
        self.protein_acc_list = []
        self.protein_hit_list = []

    def set_proteins(self, plist):
        self.protein_acc_list = plist
        self.protein_hit_list = [0] * len(plist)

    def get_number_of_proteins_with_hits(self):
        total_hit = 0
        for hit in self.protein_hit_list:
            if hit > 0:
                total_hit += 1
        return total_hit

    def get_number_of_hits_to_genome(self):
        return sum(self.potein_hit_list)

class VParse:
    def __init__(self):
        self.genomes = []
        self.gbk_dir = None
        self.gbk_files = None
        self.faa_index = None
        self.out_dir = None
        self.index_file = None

    def parse_index(self, gbk_dir, out_dir):
        self.gbk_dir = gbk_dir
        self.out_dir = out_dir
        self.gbk_files = self.parse_gbk_dir(gbk_dir)
        print(self.gbk_files)
        self.index_file = self.parse_gbk_files(self.gbk_files, self.out_dir)

    def parse_gbk_files(self, gbk_file_list, out_dir):
        index_file = os.path.join(out_dir,'index.faa')
        output_handle = open(index_file, 'w')
        for f in gbk_file_list:
            input_handle = open(f, 'r')
            for record in SeqIO.parse(input_handle, 'genbank'):
                organism = Genome(record.name,
                                  record.annotations['organism'],
                                  record.annotations['taxonomy']
                                 )
                protein_list = []
                for feature in record.features:
                    print(feature)
                    #print(dir(feature))
                    #input("")
                    if feature.type == "CDS": # TODO doubling
                        assert len(feature.qualifiers['translation']) == 1 # TODO ask about this
                        #print(feature.qaulifiers['protein_id'])
                        protein_list.append(feature.qualifiers['protein_id'][0])
                        out_string = ('>%s|ref|%s|%s\n%s\n' %
                            (feature.qualifiers['db_xref'][0].replace(':', '|'),
                             feature.qualifiers['protein_id'][0],
                             feature.qualifiers['product'][0],
                             feature.qualifiers['translation'][0]
                            ))
                        output = out_string.replace(' ', '_')
                        print('writing:')
                        print(output)
                        output_handle.write(output) # TODO ask if this is good
                organism.set_proteins(protein_list)
                self.genomes.append(organism)
                g = self.genomes[0]
            input_handle.close()
        output_handle.close()
        return index_file



    def parse_gbk_dir(self, gbk_dir):
        gbk_files = []
        for gbk_file in glob(os.path.join(gbk_dir, '*.gbk')):
            gbk_files.append(gbk_file)
        sub_dirs = glob(os.path.join(gbk_dir, '*'))
        for f in sub_dirs:
            if os.path.isdir(f):
                for gbk_file in glob(os.path.join(f,'*.gbk')):
                    gbk_files.append(gbk_file)

        if not gbk_files:
            raise ValueError('Improperly formatted index folder.\n'
                             'No GBK files found in: %s'
                             % gbk_dir)
        return gbk_files
# TEST
if __name__ == '__main__':
    gbk_dir = sys.argv[1]
    out_dir = sys.argv[2]
    vp = VParse()
    vp.parse_index(gbk_dir, out_dir)
    for g in vp.genomes:
        print(g.genome_acc, g.species)
        print(g.protein_acc_list)
     #  print('gacc: %s, spec: %s, tax: %s, proteins: %s, hits: %s\n'
     #        (g.genome_acc, g.species, g.taxonomy,
     #         g.protein_acc_list, g.protein_hit_list))
