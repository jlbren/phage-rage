import subprocess as sub
import sys
import os

class VAssemble:
    def __init__(self, finput, asm_dir, paired_reads):
        self.finput = finput
        self.asm_dir = asm_dir 
        self.paired_reads = paired_reads

    def run_qc(self, qc_dir):
        self.qc_dir = qc_dir
        self._run_sickle()

        # TODO redirect input to trimmed files
        # TODO implement logging

    def _run_sickle(self):
        s_type = 'sanger'
        s_length = '100'
        if not self.paired_reads:
            fname = os.path.basename(self.finput[0])
            p = sub.check_call([
                         'sickle',
                          self.r_rtype,
                          '-f', self.finput[0],
                          '-t', s_type,
                          '-l', s_length,
                          '-o', os.path.join(self.qc_dir, 'trimmed-'+fname)
                        ]
                       )
        else:
            fname = os.path.basename(self.finput[0])
            rname = os.path.basename(self.finput[1])
            p = sub.check_call([
                         'sickle',
                         self.r_type,
                         '-f', self.finput[0],
                         '-r', self.finput[1],
                         '-t', s_type,
                         '-l', s_length,
                         '-o', os.path.join(self.qc_dir, 'trimemed-'+fname),
                         '-p', os.path.join(self.qc_dir, 'trimmed-'+rname),
                         '-s', os.path.join(self.qc_dir, 'singletons.fastq')
                         ])
        return p
