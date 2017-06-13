import sys
import vutils
import vassemble

# pass list of arguments to VSetup
vconf = vutils.VSetup(sys.argv[1:])
vconf.check_dependencies()

vasm = vassemble.VAssemble(
                            vconf.args.finput,
                            vconf.out_dirs['assembled'],
                            vconf.args.paried_end_reads
                          )

if vconf.args.quality_control is True:
    vasm.run_qc( vconf.out_dirs['trimmed'])

if not vconf.args.assembled_contigs:
