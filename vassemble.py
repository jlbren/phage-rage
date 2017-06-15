import subprocess as sub
import sys
import os
from glob import glob

class VAssemble:
    def __init__(self, finput, asm_dir, paired_reads):
        self.finput = finput
        self.asm_dir = asm_dir
        self.paired_reads = paired_reads

    def run_qc(self, qc_dir):
        self.qc_dir = qc_dir
        self._run_sickle()
        self.finput = glob(self.qc_dir)
        # TODO redirect input to trimmed files/
        # TODO implement logging
    def run_assembly(self, assembler):
        self.assembler = assembler
        asm_input = _to_formatted_input()
        _run_assembler(asm_input)

    def _to_formatted_input(self):
        if self.paired_end_reads is True:
            trimmed = [f for f in finput if 'trimmed' in f] # glorious one liner
        if self.assembler == 'velvet':
            # Velvet just takes all files as positional args. 
            return self.finput

        elif self.assembler == 'spades':
            if self.paired_reads is True:
                pass
            else:
                pass
                
        elif self.assembler == 'megahit':
            # Megahit does not do SE or singletons. 
            return ['-1', trimmed[0], '-2', trimmed[1]]

    def _run_sickle(self):
        s_type = 'sanger'
        s_length = '100'
        foutput = []
        if not self.paired_reads:
            fname = 'trimmed-'+os.path.basename(self.finput[0])
            fouput.append(os.path.join(self.qc_dir, fname))
            p = sub.check_call([
                         'sickle',
                          'se',
                          '-f', self.finput[0],
                          '-t', s_type,
                          '-l', s_length,
                          '-o', foutput[0]
                        ]
                       )
        else:
            fname = 'trimmed-' + os.path.basename(self.finput[0])
            rname = 'trimmed-' + os.path.basename(self.finput[1])
            foutput.append(os.path.join(self.qc_dir, fname))
            foutput.append(os.path.join(self.qc_dir, rname))

            p = sub.check_call([
                         'sickle',
                         'pe',
                         '-f', self.finput[0],
                         '-r', self.finput[1],
                         '-t', s_type,
                         '-l', s_length,
                         '-o', os.path.join(self.qc_dir, 'trimemed-'+fname),
                         '-p', os.path.join(self.qc_dir, 'trimmed-'+rname),
                         '-s', os.path.join(self.qc_dir, 'singletons.fastq')
                         ])
        return p
