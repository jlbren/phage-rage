import subprocess
import sys
import os
from glob import glob

class VMap:
    def __init__(self, finput, threads):
        self.finput = finput
        self.threads = threads

    def run_map(self, mapper, index_dir, map_dir):
        self.mapper = mapper
        self.map_dir = map_dir
        self.index_dir = index_dir

    def get_orfs(self):
        pass

    def _run_mapper(self):
        if self.mapper == 'blast':
            subprocess.check_call(['blastp',
                                   ''
                                  ])

        elif self.mapper == 'pauda':
            print('Running pauda...')
            subprocess.check_call(['pauda-run',
                                   self.finpunt,
                                   self.mapper_dir,
                                   self.index_dir

                                  ])


        elif self.mapper == 'lambda':
            print('Running lambda...')
            subprocess.check_call(['lambda',
                                   '-q', self.finput,
                                   '-i', self.index_dir
                                  ])

        elif self.mapper == 'diamond':
            subprocess.check_call(['diamond', 'blastp',
                                   '-d', self.index_dir,
                                   '-q', self.finput,
                                   '-o', self.mapper_dir
                                  ])

        else:
            pass # TODO Same old problem with these cases
