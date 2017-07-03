import sys
import vutils
import vassemble
import vmap
import vparse
# pass list of arguments to VSetup
vconf = vutils.VSetup(sys.argv[1:])
#vconf.check_dependencies() TODO fix dependency checkers

vasm = vassemble.VAssemble(
                            vconf.args.finput,
                            vconf.args.paired_end_reads,
                            vconf.args.threads
                          )

if vconf.args.quality_control is True:
    vasm.run_qc(vconf.out_dirs['trimmed'])

if vconf.args.assembled_contigs is False:
    vasm.run_assembly(vconf.args.assembler, vconf.out_dirs['assembled'])
    contigs = vasm.contigs
else:
    contigs = vconf.args.assembled_contigs # TODO make nicer
vparser = vparse.VParse()
vparser.parse_index(vconf.args.index, vconf.out_dirs['mapped'])
# TODO clean up these functins, everything shouldnt be in constructor if possible 
vmapper = vmap.VMap(contigs, 
                    vconf.args.mapper, 
                    vconf.out_dirs['mapped'], 
                    vconf.args.threads)
vmapper.build_index(vparser.index_file)
vmapper.run_map()
