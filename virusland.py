import sys

import vutils
import vassemble
# pass list of arguments to VSetup
vconf = vutils.VSetup(sys.argv[1:])
vconf.check_dependencies()

vasm = vassemble.VAssemble(vconf.args, vconf.out_dirs)
if vconf.args.quality_control == True:
    vasm.run_qc()
