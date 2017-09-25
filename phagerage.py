import sys
import vutils
import vassemble
import vmap
import vparse

# pass list of arguments to VSetup
vconf = vutils.VSetup(sys.argv[1:])
vconf.check_dependencies()

# Log start message
vlog = vutils.Logger('phagerage', vconf.out_dirs['logs'])
start_msg = ('Starting phagerage pipeline.\nConfiguration:\n'
             + '\tInput: ' + ' '.join(vconf.args.finput) + '\n'
             + '\tInput type: ' + vconf.input_type + '\n'
             + '\tAssembler: ' + str(vconf.args.assembler) + '\n'
             + '\tQuality control: ' + str(vconf.args.quality_control) + '\n'
             + '\tMapper: ' + vconf.args.mapper + '\n'
             + '\tThreshold: ' + str(vconf.args.threshold) + '\n'
             + '\tThreads: ' + str(vconf.args.threads) + '\n'
             + '\tIndex directory: ' + str(vconf.args.index) + '\n'
             + '\tOutput directory: ' + str(vconf.args.output)
            )
vlog.log('main', start_msg)

# Assembly
vasm = vassemble.VAssemble(vconf.args.finput,
                           vconf.args.paired_end_reads,
                           vconf.args.threads
                          )
vasm.start_logger(vconf.out_dirs['logs'])
if vconf.args.quality_control is True:
    vasm.run_qc(vconf.out_dirs['trimmed'])

if vconf.args.assembled_contigs is False:
    vasm.run_assembly(vconf.args.assembler, vconf.out_dirs['assembled'])
    contigs = vasm.contigs
else:
    contigs = vconf.args.finput[0] # TODO make nicer

# Parse index directory to build faa file and genomes data structure
vparser = vparse.VParse()
vparser.start_logger(vconf.out_dirs['logs'])
vparser.parse_index(vconf.args.index, vconf.out_dirs['mapped'])

# Mapping
# TODO clean up these functions, everything shouldnt be in constructor if possible 
vmapper = vmap.VMap(contigs, 
                    vconf.args.mapper, 
                    vconf.out_dirs['mapped'], 
                   )
vmapper.start_logger(vconf.out_dirs['logs'])
vmapper.build_index(vparser.index_file)
vmapper.run_map(vconf.args.threads)

# Parse hits file, generate and write out statistics
vparser.parse_hits_file(vmapper.hits_file, vconf.args.threshold) # TODO make bitscore arg w default 50
vparser.generate_statistics()
vparser.write_out_all_stats(vconf.out_dirs['stats'])

# Goodbye
vlog.log('main', 'PhageRage Pipeline Finished.\n'
                 'Thank you for using this program. Goodbye!')
