import subprocess
import sys
import os
from glob import glob

class VMap:
    def __init__(self, finput, mapper, map_dir, threads):
        self.finput = finput
        self.threads = str(threads)
        self.map_dir = map_dir
        self.mapper = mapper

    def run_map(self):
        self._get_orfs()
        self._run_mapper()

    def _get_orfs(self): # TODO make sure -find is good 
        print('Running getorf...')
        self.orfs = os.path.join(self.map_dir, 'predicted_orfs.faa')
        print(self.orfs)
        print(self.finput)
        subprocess.check_call(['getorf',
                               '-find', '0', # NOTE 0 is for AA, 2 for NUC
                               '-minsize', '90',
                               self.finput,
                               self.orfs
                              ])
    def _run_mapper(self):
        if self.mapper == 'blast':
            print("Running blast...")
            subprocess.check_call(['blastp',
                                   '-query', self.orfs,
                                   '-db', self.index_dir,
                                   '-max_target_seqs', '1',
                                   '-evalue', '1',
                                   '-out', self.map_dir, #TODO is this the output file?
                                   '-outfmt', 
                                   '"10 qseqid sseqid qstart qend pident length evalue bitscore"', 
                                   '-num_threads', self.threads
                                  ])

        elif self.mapper == 'pauda':
            print('Running pauda...')
            subprocess.check_call(['pauda-run',
                                   self.orfs, 
                                   self.map_dir,
                                   self.index_dir,
                                  ])


        elif self.mapper == 'lambda': # TODO multithreads?
            print('Running lambda...')
            subprocess.check_call(['lambda', 
                                   '-p', 'blastp',
                                   '-q', self.orfs,
                                   '-i', self.index_dir,
                                   '-t', self.threads,
                                   '-o', self.map_dir + '/lambda_out.m8', # TODO make these output file names unified
                                   '--output-columns', 
                                   '"qseqid sseqid qstart qend pident length evalue bitscore"'
                                  ])

        elif self.mapper == 'diamond':
            subprocess.check_call(['diamond', 'blastp',
                                   '-d', self.index_dir,
                                   '--threads', self.threads,
                                   '-q', self.orfs,
                                   '-o', self.map_dir+'/diamond_out.csv',
                                   '--outfmt',
                                   '6 qseqid sseqid qstart qend pident length evalue bitscore'
                                  ])

        else:
            pass # TODO Same old problem with these cases


    def build_index(self, index_input):
        print("Building index...")
        self.index_dir = os.path.join(self.map_dir, 'index')
        if self.mapper == 'diamond':
            subprocess.check_call(['diamond', 'makedb',
                                   '--in', index_input,
                                   '-d', self.index_dir,
                                  ])

        elif self.mapper == 'pauda':
            subprocess.check_call(['/home/biotech3/vlib/pauda-1.0.1/bin/pauda-build',
                                   index_input,
                                   self.index_dir
                                  ])

        elif self.mapper == 'lambda':
            self.index_dir = self.index_dir + '.lambda' # lambda specific format
            subprocess.check_call(['lambda_indexer',
                                   '-d', index_input,
                                   '-i', self.index_dir
                                  ])

        elif self.mapper == 'blast':
            subprocess.check_call(['makeblastdb',
                                   '-in', index_input,
                                   '-dbtype', 'prot',
                                   '-out', self.index_dir
                    ])
        else:
            pass # TODO u already kno
# Test
if __name__ == '__main__':
    contigs = sys.argv[1]
    mapper_dir = sys.argv[2]
    index_input = sys.argv[3]
    mapper = sys.argv[4]

    m = VMap(contigs, mapper, mapper_dir, 12)
    m.build_index(index_input)
    m.run_map()
