import subprocess
import sys
import os
from glob import glob

class VMap: 
    def __init__(self, finput, threads):
        self.finput = finput
        self.threads = threads 

    def run_map(self, mapper, map_dir):
        self.mapper = mapper
        self.map_dir = map_dir
    
    def get_orfs(self):
        pass 

    def _run_mapper(self):
        if self.mapper == 'blast':
            pass

        elif self.mapper == 'pauda':
            pass

        elif self.mapper == 'lambda':
            pass 

        elif self.mapper == 'diamond':
            pass

        else:
            pass # TODO Same old problem with these cases
