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
        subprocess.check_call(['getorf',
                               '-find', '0', # NOTE 0 is for AA, 2 for NUC
                               self.finput,
                               self.orfs
                              ])
    def _run_mapper(self):
        if self.mapper == 'blast':
            print("Running blast...")
            subprocess.check_call(['blastp',
                                   '-query', self.orfs,
                                   '-db', self.index_dir,
                                   '-outfmt', '10', # TODO specific outfmt
                                   '-num_threads', self.threads
                                  ])

        elif self.mapper == 'pauda': # TODO multithreads?
            print('Running pauda...')
            subprocess.check_call(['pauda-run',
                                   self.orfs,
                                   self.map_dir,
                                   self.index_dir

                                  ])


        elif self.mapper == 'lambda': # TODO multithreads?
            print('Running lambda...')
            subprocess.check_call(['lambda',
                                   '-q', self.orfs,
                                   '-i', self.index_dir
                                  ])

        elif self.mapper == 'diamond':
            subprocess.check_call(['diamond', 'blastp',
                                   '-d', self.index_dir,
                                   '-q', self.orfs,
                                   '-o', self.map_dir,
                                   '--threads', self.threads
                                  ])

        else:
            pass # TODO Same old problem with these cases


    def build_index(self, index_input):
        print("Building index...")
        self.index_dir = os.path.join(self.map_dir, 'index')
        if self.mapper == 'diamond':
            subprocess.check_call(['diamond', 'makedb'
                                   '--in', index_input,
                                   '-d', self.index_dir,
                                  ])

        elif self.mapper == 'pauda':
            subprocess.check_call(['pauda-build',
                                   index_input,
                                   self.index_dir
                                  ])

        elif self.mapper == 'lambda':
            subprocess.check_call(['lambda-indexer'
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
