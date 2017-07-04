from glob import glob
from Bio import SeqIO
import os
import sys

class Genome:
    def __init__(self, g_acc, species, tax):
        self.genome_acc = g_acc
        self.taxonomy = tax
        self.species = species
        self.protein_list = []
        self.protein_hit_list = []
        self.total_proteins_hit = 0
        self.total_hits_to_genome = 0
        self.total_percent_hit = 0

    def set_proteins(self, plist):
        self.protein_list = plist
        self.protein_hit_list = [0] * len(plist)

    def get_number_of_proteins_with_hits(self):
        total_hit = 0
        for hit in self.protein_hit_list:
            if hit > 0:
                total_hit += 1
        return total_hit

    def generate_stats(self):
        self.total_proteins_hit = self.get_number_of_proteins_with_hits()
        self.total_hits_to_genome = self.get_number_of_hits_to_genome()
        percent_hit = self.get_percentage_of_proteins_hit()
        if percent_hit is None:
            percent_hit == '0'
        else:
            self.total_percent_hit = percent_hit
    def get_percentage_of_proteins_hit(self):
        if self.total_proteins_hit == 0:
            self.total_proteins_hit = self.get_number_of_proteins_with_hits()
        num_proteins = len(self.protein_list)
        if num_proteins > 0 and self.total_proteins_hit > 0:
            percent_hit = float((self.total_proteins_hit/num_proteins)*100.00)
            if percent_hit is not None:# TODO fix this better
                output =  ('%.2f' % percent_hit)
                return output
            else:
                return '0'

    def get_number_of_hits_to_genome(self):
        return sum(self.protein_hit_list)

class VParse:
    def __init__(self):
        self.genomes = []
        self.gbk_dir = None
        self.gbk_files = None
        self.faa_index = None
        self.out_dir = None
        self.index_file = None

    def parse_hits_file(self, hits_file, threshold):
        print('Parsing hit file...')
        input_handle = open(hits_file, 'r')
        for hit in input_handle:
            fields = hit.split('\t')
            sseqid = fields[1]
            bitscore = float(fields[7])
            if bitscore <= threshold:
                self.update_hit(sseqid)
        input_handle.close()
        return True

    def update_hit(self, sseqid):
        fields = sseqid.split('|')
        pid = fields[3]
        for genome in self.genomes: # TODO ask about optimization
            if pid in genome.protein_list:
                g_index = genome.protein_list.index(pid)
                genome.protein_hit_list[g_index] += 1
                return True
        return False

    def generate_statistics(self):
        print('Generating genome statistics..')
        for genome in self.genomes:
            genome.generate_stats()
        return True

    def write_out_all_stats(self, stats_dir):
        self.write_out_summary_stats(stats_dir)
        self.write_out_for_hitviz(stats_dir)
        self.write_out_for_krona(stats_dir)

    def write_out_for_hitviz(self, stats_dir):
        hitviz_file = os.path.join(stats_dir, 'hitviz_stats.csv')
        print('Writing out for HitViz to:', hitviz_file)
        with open(hitviz_file, 'w') as output_handle:
            for genome in genomes:
                if genome.total_hits_to_genome > 0:
                    output_handle.write("%s|%s|\n" %
                              (genome.genome_acc,
                               genome.species
                              )
                            )
                    for i in range(len(genome.protein_list)):
                        output_handle.wrtie('%s|%s|\n' %
                                (genome.protein_list[i],
                                 genome.protein_hit_list[i]
                                )
                               )
                    output_handle.write('*\n')

    def write_out_summary_statistics(self, stats_dir):
        print('Writing out summary statistics...')
        coverage_file = os.path.join(stats_dir, 'coverage.csv')
        with open(coverage_file, 'w') as output_handle:
            output_handle.write('Accession Number\t'
                                'Spp\t'
                                'Number of Proteins in Genome\t'
                                'Number of Hits to Genome\t'
                                'Percent Proteins Hit\n'
                                'Number of genomes in database: %d\n'
                                % len(self.genomes)
                               )
            for genome in self.genomes:
                if genome.total_hits_to_genome > 0:
                    genome_stats = [genome.genome_acc,
                                    genome.species,
                                    str(len(genome.protein_list)),
                                    str(genome.total_hits_to_genome),
                                    str(genome.total_percent_hit)
                                   ]
                    output = '\t'.join(genome_stats)
                    output_handle.write('%s\n' % output)

        file_name = 'hits_by_protein.csv'
        protein_hits_file = os.path.join(stats_dir, file_name)
        with open(protein_hits_file ,'w') as output_handle:
            for genome in self.genomes:
                if genome.total_hits_to_genome > 0:
                    output = '\t'.join([genome.genome_acc, genome.species])
                    output_handle.write('>%s\n' % output) # TODO confirm this is ok
                    for i in range(len(genome.protein_list)):
                        if genome.protein_hit_list[i] > 0:
                            output = '\t'.join([genome.protein_list[i],
                                                str(genome.protein_hit_list[i])
                                               ])
                            output_handle.write('%s\n' % output)

    def parse_index(self, gbk_dir, out_dir):
        print('Parsing GBK directory:', gbk_dir)
        self.gbk_dir = gbk_dir
        self.out_dir = out_dir
        self.gbk_files = self.parse_gbk_dir(gbk_dir)
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
                    if (feature.type == "CDS"
                         and 'translation' in feature.qualifiers): # kill pseudo genes
                        assert len(feature.qualifiers['translation']) == 1
                        protein_list.append(feature.qualifiers['protein_id'][0])
                        out_string = ('>%s|ref|%s|%s\n%s\n' %
                            (feature.qualifiers['db_xref'][0].replace(':', '|'),
                             feature.qualifiers['protein_id'][0],
                             feature.qualifiers['product'][0],
                             feature.qualifiers['translation'][0]
                            ))
                        output = out_string.replace(' ', '_')
                        output_handle.write(output) # TODO add record.name to output field for optimization
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
    gbk_input = sys.argv[1]
    out_dir = sys.argv[2]
    hits_input = sys.argv[3]
    vp = VParse()
    vp.parse_index(gbk_input, out_dir)
    #for g in vp.genomes:
    #    print(g.genome_acc, g.species)
    #    print(g.protein_list)

    vp.parse_hits_file(hits_input, .005)
    vp.generate_statistics()
    vp.write_out_summary_statistics(out_dir)
