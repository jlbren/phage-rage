import subprocess as sub
import sys
import os
class VAssemble:
    def __init__(self, finput, qc_flag, r_type):
        self.finput = args.finput
        self.quality_control = qc_flag 
        self.r_type = 'SE'

    def run_qc(self, qc_dir):
        self.qc_dir = qc_dir 
        p = self._run_sickle()
        print(p.stdout)
        print(p.stderr, file=sys.stderr)
        # TODO redirect input to trimmed files
        # TODO implement logging

    def _run_sickle(self):
        s_type = 'sanger'
        s_length = '100' 
        if self.r_type == 'SE':
            fname = os.path.basename(self.finput[0])
            p = sub.run(['sh',
                         'sickle',
                          self.r_rtype,
                          '-f', self.finput[0],
                          '-t', s_type,
                          '-l', s_length,
                          '-o', os.path.join(self.qc_dir, 'trimmed-'+fname)

                        ],
                        stdout=sub.PIPE, stderr=sub.PIPE, check=True
                       )
        else:
            fname = os.path.basename(self.finput[0])
            rname = os.path.basename(self.finput[1])
            p = sub.run([
                         'sickle',
                         self.r_type,
                         '-f', self.finput[0],
                         '-r', self.finput[1],
                         '-t', s_type,
                         '-l', s_length,
                         '-o', os.path.join(self.qc_dir, 'trimemed-'+fname),
                         '-p', os.path.join(self.qc_dir, 'trimmed-'+rname),
                         '-s', os.path.join(self.qc_dir, 'singletons.fastq')
                        ],
                        stdout=sub.PIPE, stderr=sub.PIPE, check=True
                       )
            return p
