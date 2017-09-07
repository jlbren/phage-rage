from glob import glob
from Bio import SeqIO
import os
import sys
import subprocess
import vutils

class Genome:
    """Class for storing individual genome information and
       basic statistics.
    Arguments:
        g_acc (str): Genome accession number.
        species (str): Genome species name.
        tax (str): Genome taxonomy, tab separated.
    Attributes:
        genome_acc (str): Genome accession number.
        taxonomy (str): Genome taxonomy, tab separated.
        species (str): Genome species name.
        protein_list (list(str)): List of protein access numbers
                                 belonging to the genome.
        protein_hit_list (str(int)): List of hits to each protein. Index
                                     corresponds to protein at the same
                                     index in protein_list.
        total_proteins_hit (int): Total number of proteins with at least
                                  one hit.
        total_hits_to_genome (int): Total number of hits across all
                                    proteins.
        total_percent_hit (float): Percentage of proteins with hits over
                                   total number of proteins in the genome.
    """

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
        """Method for initializing a genome's protein list.
           Creates proein_hit_list of same size with initial vals of 0.
        Arguments:
            plist (list(str)): Genome's list of protein acc nums parsed from the GBK file.
        """
        self.protein_list = plist
        self.protein_hit_list = [0] * len(plist)

    def get_number_of_proteins_with_hits(self):
        """Method that returns the number of proteins with non-zero hits.
        Returns:
            total_hit (int): Total number of proteins in a genome with at least one hit.
        """
        total_hit = 0
        for hit in self.protein_hit_list:
            if hit > 0:
                total_hit += 1
        return total_hit

    def generate_stats(self):
        """Public method for invoking calls to generate summary statistics."""
        self.total_proteins_hit = self.get_number_of_proteins_with_hits()
        self.total_hits_to_genome = self.get_number_of_hits_to_genome()
        percent_hit = self.get_percentage_of_proteins_hit()
        if percent_hit is None:
            percent_hit = '0'
        self.total_percent_hit = percent_hit

    def get_percentage_of_proteins_hit(self):
        """Method that calculates the total percentage of proteins in a genome with hits.
        Returns:
            (float): Percentage of number of proteins with hits / total num of proteins.
        """
        if self.total_proteins_hit == 0:
            self.total_proteins_hit = self.get_number_of_proteins_with_hits()
        num_proteins = len(self.protein_list)
        if num_proteins > 0 and self.total_proteins_hit > 0:
            percent_hit = float((self.total_proteins_hit/num_proteins)*100.00)
            if percent_hit is not None:# TODO fix this better
                output = ('%.2f' % percent_hit)
                return output
            else:
                return '0'

    def get_number_of_hits_to_genome(self):
        """Method to get the total sum of hits to the genome.
        Returns:
            (int): Sum of total hits to genome.
        """
        return sum(self.protein_hit_list)

class VParse:
    """Class for managing parsing index and hit files,
       Managing genome data structures, and writting output.
    Attributes:
        genomes (list(Genome)): Array of Genome objects representing all
                                records in the index.
        gbk_dir (str): Path to directory containing index GBK files.
        gbk_files (list(str)): List of all GBK files in index dir
                               or one level of subdir.
        out_dir (str): Path to output direcatory.
        index_file (str): Path to generated FAA file used to build mapper index.
    """
    def __init__(self):
        self.genomes = []
        self.gbk_dir = None
        self.gbk_files = None
        self.out_dir = None
        self.index_file = None
    def start_logger(self, log_dir):
        """Boiler plate method for starting this class' logger."""
        self._logger = vutils.Logger('vparse', log_dir)

    def parse_index(self, gbk_dir, out_dir):
        """Public method for initiaing the index parsing methods. 
        Arguments: 
            gbk_dir (str): Path to root directory containing
                           GBK files to be parsed. 
            out_dir (str): Path to write out parsed output. 
        """
        print('Parsing GBK directory:', gbk_dir)
        self.gbk_dir = gbk_dir
        self.out_dir = out_dir
        self._logger.log('parse_index', 'Starting to parse GBK dir: %s' % self.gbk_dir)
        self.gbk_files = self.parse_gbk_dir(gbk_dir)
        self.index_file = self.parse_gbk_files(self.gbk_files, self.out_dir)
        self._logger.log('parse_index', 'Finished parsing GBK directory.\n'
                         'FAA output: %s\nTotal number of genomes: %d'
                         % (self.index_file, len(self.genomes))
                        )

    def parse_gbk_files(self, gbk_file_list, out_dir):
        """Method for parsing GBK input files to create the FAA index file and
           populate genome data structure. 
        Arguments:
            gbk_file_list (list(str)): List of file paths to all GBK files. 
            out_dir (str): Path to output directory. 
        """
        index_file = os.path.join(out_dir, 'index.faa')
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
                        and 'translation' in feature.qualifiers):  # kill pseudo genes
                        assert len(feature.qualifiers['translation']) == 1
                        protein_list.append(feature.qualifiers['protein_id'][0])
                        out_string = ('>%s|ref|%s|%s\n%s\n' %
                                      (feature.qualifiers['db_xref'][0].replace(':', '|'),
                                       feature.qualifiers['protein_id'][0],
                                       feature.qualifiers['product'][0],
                                       feature.qualifiers['translation'][0]
                                       ))
                        output = out_string.replace(' ', '_')
                        output_handle.write(output)  # TODO add record.name to output field for optimization
                organism.set_proteins(protein_list)
                self.genomes.append(organism)
                g = self.genomes[0]
            input_handle.close()
        output_handle.close()
        return index_file

    def parse_gbk_dir(self, gbk_dir):
        """Method to parse specified index directory and return a list of GBK files.
           GBK files can either be in the root dir or within one level of sub directory. 
        Arguments:
            gbk_dir (str): Path to the provided index dir containing GBK files. 
        Returns:
            gbk_files (list(str)): List of all GBK files. 
        Raises:
            ValueError: Exception raised if no GBK files can be found
                        in the provided directory. 
        """
        gbk_files = []
        for gbk_file in glob(os.path.join(gbk_dir, '*.gbk')):
            gbk_files.append(gbk_file)
        sub_dirs = glob(os.path.join(gbk_dir, '*'))
        for f in sub_dirs:
            if os.path.isdir(f):
                for gbk_file in glob(os.path.join(f, '*.gbk')):
                    gbk_files.append(gbk_file)

        if not gbk_files:
            raise ValueError('Improperly formatted index folder.\n'
                             'No GBK files found in: %s'
                             % gbk_dir)
        return gbk_files

    def parse_hits_file(self, hits_file, threshold):
        """Method for parsing hits file in order to update the count in the genome list. 
        Arguments:
            hits_file (str): Path to hit file produced by the mapper. 
            threshold (float): Bitscore threshold for evaluating hit quality.
        """
        print('Parsing hit file...')
        self._logger.log('parse_hits_file', 'Parsing hit file: %s'
                                            % hits_file
                        )
        input_handle = open(hits_file, 'r')
        for hit in input_handle:
            fields = hit.split('\t')
            sseqid = fields[1]
            bitscore = float(fields[7])
            if bitscore >= threshold:
                self.update_hit(sseqid)
        input_handle.close()
        return True

    def update_hit(self, sseqid):
        """Method for incrementing the appropriate hit counter when mapped to a protein.
        Arguments:
            sseqid (str): ID string of the subject sequence in a mapped hit 
                          (protein accession num)
        Returns:
            (bool): True if protein is found in genome index, false if not found. 
        """
        fields = sseqid.split('|')
        pid = fields[3]
        for genome in self.genomes: # TODO ask about optimization
            if pid in genome.protein_list:
                g_index = genome.protein_list.index(pid)
                genome.protein_hit_list[g_index] += 1
                return True
        return False

    def generate_statistics(self):
        """Public method for initiating a call on each genomes generate_stats() method."""
        print('Generating genome statistics..')
        for genome in self.genomes:
            genome.generate_stats()
        self.genomes.sort(key = lambda species: genome.species)
        return True

    def write_out_all_stats(self, stats_dir):
        """Public method for calling each separate write method.
        Arguments:
            stats_dir (str): Directory to write all output files to.
        """
        self._logger.log('write_out_all_stats', 'Writting out stats files.\n'
                                                'Output directory: %s' % stats_dir
                        )
        self.write_out_summary_statistics(stats_dir)
        self.write_out_for_hitviz(stats_dir)
        self.write_out_for_krona(stats_dir)

    def write_out_for_krona(self, stats_dir):
        """Method for producing output file to be read by Krona.
        Arguments:
             stats_dir (str): Directory output file is located. 
        Outputs:
            krona_stats.csv: Hits to Genome, Taxonomy, Species; tab separated. 
        """
        krona_file = os.path.join(stats_dir, 'krona_stats.csv')
        print('Writing out for Krona to:', krona_file)
        with open(krona_file, 'w') as output_handle:
            for genome in self.genomes:
                if genome.total_hits_to_genome > 0:
                    tax = '\t'.join(genome.taxonomy)
                    output_handle.write('%d\t%s\t%s\n' %
                            (genome.total_hits_to_genome,
                             tax,
                             genome.species
                            )
                           )
        krona_graph = os.path.join(stats_dir, 'krona_graph.html')
        subprocess.check_call(['ktImportText', krona_file,
                               '-o', krona_graph
                              ])

    def write_out_for_hitviz(self, stats_dir):
        """Method for proudcing output file to be read by hitviz. 
        Arguments:
             stats_dir (str): Directory output file is located.
        Outputs: 
            hitviz_stats.csv: Number of hits to each protein in each genome.
        """
        hitviz_file = os.path.join(stats_dir, 'hitviz_stats.csv')
        print('Writing out for HitViz to:', hitviz_file)
        with open(hitviz_file, 'w') as output_handle:
            for genome in self.genomes:
                if genome.total_hits_to_genome > 0:
                    output_handle.write("%s|%s|\n" %
                                    (genome.genome_acc,
                                     genome.species
                                    )
                                   )
                    for i in range(len(genome.protein_list)):
                            output_handle.write('%s|%d|\n' %
                                            (genome.protein_list[i],
                                             genome.protein_hit_list[i]
                                            )
                                           )
                    output_handle.write('*\n')

    def write_out_summary_statistics(self, stats_dir):
        """Method for proudcing summary statistic output files. 
        Arguments:
             stats_dir (str): Directory output file is located. 
        Outputs: 
            coverage.csv: Coverage information on each genome with 
                          at least one hit. 
            hits_by_protein.csv: Accession number of each protein with number of hits. 
        """
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
                    output_handle.write('>%s\n' % output) # TODO confirm > is ok
                    for i in range(len(genome.protein_list)):
                        if genome.protein_hit_list[i] > 0:
                            output = '\t'.join([genome.protein_list[i],
                                                str(genome.protein_hit_list[i])
                                               ])
                            output_handle.write('%s\n' % output)


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

    vp.parse_hits_file(hits_input, 50)
    vp.generate_statistics()
    vp.write_out_all_stats(out_dir)
