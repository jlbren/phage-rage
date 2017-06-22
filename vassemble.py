import subprocess
import sys
import os
from glob import glob

class VAssemble:
    def __init__(self, finput, paired_end_reads, threads):
        self.finput = finput
        self.paired_end_reads = paired_end_reads
        self.threads = str(threads)
        self.qc = False

    def run_qc(self, qc_dir):
        self.qc_dir = qc_dir
        self._run_sickle()
        self.qc = True
        # TODO implement logging

    def run_assembly(self, assembler, asm_dir):
        self.assembler = assembler
        self.asm_dir = asm_dir
        asm_input = self._to_formatted_input()
        self._run_assembler(asm_input)

    def _run_assembler(self, asm_input):
        if self.assembler == 'spades': # TODO check if this how u call spades
            print('Running spades..')
            p = subprocess.check_call(['spades.py',
                                       '--threads', self.threads,
				       '-o', self.asm_dir]
                                      + asm_input)

        elif self.assembler == 'velvet':
            print('Running velevth...')
            p = subprocess.check_call(['velveth', self.asm_dir,
                                       '31', '-fastq', '-shortPaired'
                                      ]
                                      + asm_input)
            print('Running velvetg...')
            p = subprocess.check_call(['velvetg', self.asm_dir,
                                       '-exp_cov', 'auto'
                                      ])

        elif self.assembler == 'megahit':
            print('Running megahit...')
            p = subprocess.check_call(['megahit',
                                       '-t', self.threads,
                                       '-o', self.asm_dir+'/megahit' 
                                      ]
				      + asm_input)

        else:
            pass # TODO Same problem w default cases that shouldnt happen

    def _run_sickle(self):
        s_type = 'sanger'
        s_length = '100'
        foutput = []
        if not self.paired_end_reads:
            fname = 'trimmed-'+os.path.basename(self.finput[0])
            foutput.append(os.path.join(self.qc_dir, fname))
            p = subprocess.check_call([
                          'sickle',
                          'se',
                          '-f', self.finput[0],
                          '-t', s_type,
                          '-l', s_length,
                          '-o', foutput[0]
                         ])

        else:
            fname = 'trimmed-' + os.path.basename(self.finput[0])
            rname = 'trimmed-' + os.path.basename(self.finput[1])
            foutput.append(os.path.join(self.qc_dir, fname))
            foutput.append(os.path.join(self.qc_dir, rname))

            p = subprocess.check_call([
                         'sickle',
                         'pe',
                         '-f', self.finput[0],
                         '-r', self.finput[1],
                         '-t', s_type,
                         '-l', s_length,
                         '-o', os.path.join(self.qc_dir, fname),
                         '-p', os.path.join(self.qc_dir, rname),
                         '-s', os.path.join(self.qc_dir, 'singletons.fastq')
                        ])
        return p

    def _to_formatted_input(self):
        if self.qc is True:
            return self._get_trimmed_input()
        else:
            if self.paired_end_reads is True:
                if self.assembler == 'velvet':
                    return self.finput

                elif self.assembler == 'spades':
                    return [
                            '--pe1-1', self.finput[0],
                            '--pe1-2', self.finput[1],
                            '--meta'
                           ]

                elif self.assembler == 'megahit':
                    return [
                            '-1', self.finput[0],
                            '-2', self.finput[1],
                           ]

                else:
                    pass # TODO it should never get to this from args checking,
                         # but add error throw for completeness

            else: # Else if single_end_reads
                if self.assembler == 'velvet':
                    return self.finput

                elif self.assembler == 'spades':
                    return ['-s', self.finput[0]]

                elif self.assembler == 'megahit':
                    return ['-r', self.finput[0]]

                else:
                    pass # TODO see above todo


    def _get_trimmed_input(self):
        trimmed = glob(os.path.join(self.qc_dir, 'trimmed*'))
        if self.paired_end_reads is True:
            if self.assembler == 'velvet':
                return glob(os.path.join(self.qc_dir, '*'))

            elif self.assembler == 'spades':
                return [
                        '--pe1-1', trimmed[0],
                        '--pe1-2', trimmed[1],
                        '--pe1-s', os.path.join(self.qc_dir, 'singletons.fastq'),
                        '--meta'
                       ]

            elif self.assembler == 'megahit':
                return [
                        '-1', trimmed[0],
                        '-2', trimmed[1],
                       ]

            else:
                pass # TODO it should never get to this from args checking,
                     # but add error throw for completeness

        else: # Else if single_end_reads
            if self.assembler == 'velvet':
                return trimmed

            elif self.assembler == 'spades':
                return ['-s', trimmed[0]]

            elif self.assembler == 'megahit':
                return ['-r', trimmed[0]]

            else:
                pass # TODO see above todo

