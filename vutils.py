#!bin/python3
import argparse
import sys
import os

class VUtil:
    def __init__(self, argv_):
        self.argv_ = argv_[1:]
        self._parser = self._varg_parse()
        self.args = self._parser.parse_args(self.argv_)
        self._check_parser()
        self.out_dirs = self._create_output_dirs()

    def _varg_parse(self):
        parser = argparse.ArgumentParser(description="VLand 2: PHAGE RAGE",
                                         prog="virusland")
        # Input files: 2 PE, 1 SE, or 1 pre-assumbled contigs
        parser.add_argument('input_file', nargs='+',
                            help='Input file(s). Specify 2 PE, 1 SE or '
                                 '1 assembled contigs file(s).')

        # Assembler type or assembled contig flag
        group_a = parser.add_mutually_exclusive_group()
        group_a.add_argument('-a', '--assembler',
                             choices=[
                                 'spades',
                                 'velvet',
                                 'megahit'
                                 ],
                             help='Assembly method.')
        group_a.add_argument('-A', '--assembled_contigs',
                             action='store_true',
                             help='Assembled contigs flag.\n'
                                  'Start analysis from an existing assembly.\n'
                                  'Requires specifying 1 input contigs file.')

        # Read type: single or paired end.
        group_r = parser.add_mutually_exclusive_group()
        group_r.add_argument('-s', '--single_reads',
                             action='store_true',
                             help='Single reads flag.\n'
                                  'Requires specifying 1 input read file.')
        group_r.add_argument('-p', '--paired_end_reads',
                             action='store_true',
                             help='Paired-end reads flag.\n'
                             'Requires specifying 2 input read files.')

        # Number of threads.
        parser.add_argument('-t', '--threads', type=int, default=1,
                            help='Number of threads.')

        #Output base directory.
        parser.add_argument('-o', '--output', default=os.getcwd(),
                            help='Base output directory path.\n'
                            'All output will be located here.')
        # Read quality control flag
        parser.add_argument('-q', '--quality_control', action='store_true',
                            help='Read quality control flag.')

        return parser

    def _check_parser(self):
        # Check input type flags
        if sum([self.args.assembled_contigs,
                    self.args.paired_end_reads, self.args.single_reads]) != 1:
            self._parser.error('Either single reads, paired end reads, '
                               'or assembled contigs '
                               'must be exclusively specified.')
        # Check number of input files
        if self.args.assembled_contigs is True or self.args.single_reads is True:
            if len(self.args.input_file) != 1:
                self._parser.error('A single input file must be specified for '
                                   'single reads or assembled contigs.')
        elif self.args.paired_end_reads is True:
            if len(self.args.input_file) != 2:
                self._parser.error('Two input files must be specified for '
                                   'paired end reads')
        # Check assembler is specified when not using contigs
        if self.args.assembled_contigs is False and self.args.assembler is None:
            self.parser.error('Assembler choice must be spefied when not using '
                              'assembled contigs.')

        # Check if running on 1 thread, and if so print warning to stderr.
        if self.args.threads == 1:
            print('####WARNING####\n'
                  'Running on 1 thread may severly increase run time.\n'
                  'See help for options to change thread count.',
                  file=sys.stderr)

    def _create_output_dirs(self):
        out_path = self.args.output
        if not os.path.exists(out_path):
            try: 
                os.makedirs(out_path)
                print("Creating new output directory:", out_path)
            except: 
                    raise FileNotFoundError('Output directory '
                                            + out_path + 'could not be created. '
                                            'Check argument path and try again.')
        out_dirs = {}
        dir_list = [
                    'stats',
                    'logs',
                   ]
        for fdir in dir_list:
            out_dirs[fdir] = self._make_dir(out_path, fdir)

        if self.args.assembled_contigs is True:
            out_dirs['assembled'] = self.args.input_file
        else:
            out_dirs['assembled'] = self._make_dir(out_path, 'assembled')
            if self.args.quality_control is True:
                out_dirs['quality_control'] = self._make_dir(out_path,
                                                             'quality_control')

        return out_dirs

    def _make_dir(self, out_path, fdir):
        fpath = os.path.join(out_path, fdir)
        if not os.path.exists(fpath):
            os.makedirs(fpath)
            return fpath
        else:
            raise FileExistsError('Directory: ' + fpath + ' already exists.\n'
                                  'Please specify an unused output directory.')


v = VUtil(sys.argv)
print(v.args)
print(v.out_dirs)