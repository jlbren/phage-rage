#!bin/python3
import argparse
import sys
import os
import shutil
from glob import glob
from datetime import datetime

class VSetup:
    """Class used to parse program arguments, check dependencies,  create output directoties, and 
       and assist in managing run time parameters. 

       Attributes:
            argv_ (list): list of passed command line args 
            _parser (ArgParser): ArgParser object used to parse argv_
            args (object): object created by ArgParser with each arg as an 
                           attribute 
            out_dirs (dict): dictionary containing path to each output directory
            input_type (string): parsed input type (SE, PE, contigs)
       Arguments:
            argv_ (list): sys.argv[1:] passed from the main script
    """
    def __init__(self, argv_):
        self.argv_ = argv_
        self._parser = self._varg_parse()
        self.args = self._parser.parse_args(self.argv_)
        self._check_parser()
        self.out_dirs = self._create_output_dirs()
        self.input_type = self._get_input_type()
    
    def _varg_parse(self):
        """Method that defines rules for the self._parser ArgumentParser object
           Returns:
                parser (ArgumentParser): initialized ArgumentParser object 
        """
        parser = argparse.ArgumentParser(description="PHAGE RAGE",
                                         prog="phagerage")
        # Input files: 2 PE, 1 SE, or 1 pre-assembled contigs
        parser.add_argument('finput', nargs='+',
                            help='Input file(s). Specify 2 PE, 1 SE or '
                                 '1 assembled contigs file(s).')

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

        # Assembler type or assembled contig flag
        group_a = parser.add_mutually_exclusive_group()
        group_a.add_argument('-a', '--assembler',
                             choices=[
                                 'spades',
                                 'velvet',
                                 'megahit'
                             ],
                             help='Selected assembly utility.')
        group_a.add_argument('-A', '--assembled_contigs',
                             action='store_true',
                             help='Assembled contigs flag.\n'
                                  'Start analysis from an existing assembly.\n'
                                  'Requires specifying 1 input contigs file.')

        # Read quality control flag
        parser.add_argument('-q', '--quality_control', action='store_true',
                            help='Read quality control flag.\n'
                                 'Trim input read(s) with sickle.')

        # Mapper type
        parser.add_argument('-m', '--mapper', required=True,
                            choices=[
                                'blastp',
                                'lambda',
                                'diamond'
                            ],
                            help='Selected mapping utility.')

        # Path to database files
        parser.add_argument('-i', '--index', required=True,
                            help='Path to directory containing GBK files ' 
                                 'to be used when building mapper index.\n'
                                 'Files may either be within the given directory '
                                 'and/or one level of subdirectory.')

        # Number of threads.
        parser.add_argument('-t', '--threads', type=int, default=1,
                            help='Number of threads.')
        # Bitscore threshold.
        parser.add_argument('--threshold', type=int, default=50,
                             help='Bitscore threshold for parsing mapper hits.')

        # Output base directory.
        parser.add_argument('-o', '--output', type=str, default=os.getcwd(),
                            help='Base output directory path.\n'
                            'All output will be located here.\n'
                            'Must be relative to current working directory.')

        return parser

    def _check_parser(self):
        """Method that checks passed arguments for logic 
           not implicitly handled by the ArgumentParser.
           
           Raises:
                ParserError: ArgumentParser error if args are incorrect. 
        """

        # Check input type flags
        if sum([self.args.assembled_contigs,
                    self.args.paired_end_reads, self.args.single_reads]) != 1:
            self._parser.error('Either single reads, paired end reads, '
                               'or assembled contigs '
                               'must be exclusively specified.')
        # Check number of input files
        if self.args.assembled_contigs is True or self.args.single_reads is True:
            if len(self.args.finput) != 1:
                self._parser.error('A single input file must be specified for '
                                   'single reads or assembled contigs.')
        elif self.args.paired_end_reads is True:
            if len(self.args.finput) != 2:
                self._parser.error('Two input files must be specified for '
                                   'paired end reads')
        # Check input file path(s) exist.
        for f in self.args.finput:
            if not os.path.exists(f):
                self._parser.error('Input file: ' + f + ' could not be found.')
        # Check index directory exists.
        if not os.path.exists(self.args.index):
            self._parser.error('Directory %s could not be found.' % self.args.index)
        # Check assembler is specified when not using contigs
        if self.args.assembled_contigs is False and self.args.assembler is None:
            self.parser.error('Assembler choice must be spefied when not using '
                              'assembled contigs.')

        # Check if running on 1 thread, and if so print warning to stderr.
        if self.args.threads == 1:
            print('####WARNING####\n'
                  'Running on 1 thread may severly increase run time.\n'
                  'See help for options to change thread count.')

    def _get_input_type(self):
        """Determines input type based on passed arguments. 
           Returns: 
                (str): corresponding string for each input type. 
        """
        if self.args.assembled_contigs:
            return 'assembled contigs'
        elif self.args.paired_end_reads:
            return 'paired end reads'
        else:
            return 'single end read'

    def _create_output_dirs(self):
        """Method that creates output directories based on root path passed
           as a program argument. 

           Output:
                Subdirectories withing the base output path provided at run time. 
           Returns:
                out_dirs (dict): Dictionary of output directory paths. 
           Raises: 
                FileNotFoundError: Exception thrown if base output path cannot be located. 
        """
        if self.args.output != os.getcwd():
            self.args.output = os.path.join(os.getcwd(), self.args.output)
        out_path = self.args.output
        if not os.path.exists(out_path):
            try:
                os.makedirs(out_path)
                print("Creating new output directory:", out_path)
            except:
                    raise FileNotFoundError('Output directory '
                                            + out_path + 'could not be created.\n'
                                            'Check argument path and try again.')
        out_dirs = {}
        dir_list = [
                    'stats',
                    'logs',
		    'mapped'
                   ]
        for fdir in dir_list:
            out_dirs[fdir] = self._make_dir(out_path, fdir)

        if self.args.assembled_contigs is True:
            out_dirs['assembled'] = self.args.finput
        else:
            out_dirs['assembled'] = self._make_dir(out_path, 'assembled')
            if self.args.quality_control is True:
                out_dirs['trimmed'] = self._make_dir(out_path, 'trimmed')

        return out_dirs

    def _make_dir(self, out_path, fdir):
        """Wrapper function to create a subdirectory within the base output folder.
           Raises: 
                FilesExistsError: Exception thrown if a given subdirectory already
                                  exists within the base output path. Prevents the
                                  overwritting of data from previous runs. 
        """
        fpath = os.path.join(out_path, fdir)
        if not os.path.exists(fpath):
            os.makedirs(fpath)
            return fpath
        else:
            raise FileExistsError('Directory: ' + fpath + ' already exists.\n'
                                  'Please specify an unused output directory.')
    def _get_dependency(self, name, cat):
        """Method that returns the specific program call(s) based on the
           dependencies specified at runtime. 

           Arguments:
                name(str): dependency name as specified at runtime. 
                cat (str): Category to which the dependency belongs. 
           Returns:
                (list): List of program calls to be executed for the specified
                        dependency.
        """
        if cat == 'assembler':
            if name == 'spades':
                return ['spades.py']
            elif name == 'velvet':
                return['velveth', 'velvetg']
            elif name == 'megahit':
                return ['megahit']
        elif cat  == 'mapper':
            if name == 'blastp':
                return ['makeblastdb', 'blastp']
            elif name == 'lambda':
                return ['lambda_indexer','lambda']
            elif name == 'diamond':
                return ['diamond']

    def check_dependencies(self):
        """Wrapper function to itterate through each dependency
           and verify its existence in the system path. 

           Raises:
                FileNotFoundError: Exception raised if shutil.which cannot find a particular 
                                   dependency. 

        """
        depend_list = ['ktImportText', 'getorf']
        if self.args.assembled_contigs is False:
            depend_list += self._get_dependency(self.args.assembler, 'assembler')
        if self.args.quality_control is True:
            depend_list.append('sickle')
        depend_list += self._get_dependency(self.args.mapper, 'mapper')
        for prog in depend_list:
            if shutil.which(prog, mode=os.X_OK) is None:
                raise FileNotFoundError('Could not find dependency: '
                                        + prog + '\n'
                                        'Looking in:\n'
                                        + os.environ['PATH']
                                       )
        # TODO check java version for hitviz

def copy_and_remove(src, dest):
    """Wrapper function for moving the entire contents of a directory"""
    for f in glob(os.path.join(src,'*')):
        shutil.move(f, dest)
    shutil.rmtree(src)

class Logger:
    """Logger class to handle standardized output of logs for each module.
       Attributes:
          module (str): Name of module which has initialized the logger. 
          out_dir (str): Output direcotry containing log file. 
          log_file(str): Full path to log file. 

       Arguments: 
            module (str): Name of the module creating a new logger. 
            out_dir (str): Path to output directory. 
       Output: 
            Formatted log strings in log.txt file. 
    """

    def __init__(self, module, out_dir):
        self.module = module
        self.out_dir = out_dir
        self.log_file = os.path.join(out_dir, 'log.txt')

    def get_time(self):
        """Function for generating formated time strings. 
        Returns: 
            (str): Formatted current time.
        """
        return datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    def log(self, method, msg):
        """Method for writing a formated entry in the log file. 
        Arguments:
            method (str): Name of the method invoking a log call. 
            msg (str): Message to be logged. 
        Output: 
            Formatted entry in the log.txt file found in the log subdir of the
            base output path. 
        """
        time = self.get_time()
        with open(self.log_file, 'a+') as output_handle:
            output_handle.write('[%s][%s][%s]\n' % (time, self.module, method))
            output_handle.write('%s\n' % msg)

if __name__ == '__main__':
    v = VSetup(sys.argv[1:])
    print(v.args)
    print(v.out_dirs)
    v.check_dependencies()
