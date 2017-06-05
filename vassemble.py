import subprocess as sub
import sys
import os
class VAssemble:
    def __init__(self, args, out_dirs):
        self.input_file = args.input_file
        self.output_dir = out_dirs['assembled']
        if args.quality_control is True:
            self.qc_dir = out_dirs['trimmed']
        if args.single_reads is True:
            self.r_type = 'SE'
        else:
            self.r_type = 'PE'

    def run_qc(self):
        p = self._run_sickle()
        print(p.stdout)
        print(p.stderr, file=sys.stderr)
        # TODO redirect input to trimmed files
        # TODO implement logging

    def _run_sickle(self):
        s_type = 'sanger'
        s_length = '100' 
        if self.r_type == 'SE':
            fname = os.path.basename(self.input_file[0])
            p = sub.run(['sh',
                         'sickle',
                          self.r_rtype,
                          '-f', self.input_file[0],
                          '-t', s_type,
                          '-l', s_length,
                          '-o', os.path.join(self.qc_dir, 'trimmed-'+fname)

                        ],
                        stdout=sub.PIPE, stderr=sub.PIPE, check=True
                       )
        else:
            fname = os.path.basename(self.input_file[0])
            rname = os.path.basename(self.input_file[1])
            p = sub.run([
                         'sickle',
                         self.r_type,
                         '-f', self.input_file[0],
                         '-r', self.input_file[1],
                         '-t', s_type,
                         '-l', s_length,
                         '-o', os.path.join(self.qc_dir, 'trimemed-'+fname),
                         '-p', os.path.join(self.qc_dir, 'trimmed-'+rname),
                         '-s', os.path.join(self.qc_dir, 'singletons.fastq')
                        ],
                        stdout=sub.PIPE, stderr=sub.PIPE, check=True
                       )
            return p
