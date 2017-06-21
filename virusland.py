import sys
import vutils
import vassemble

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

