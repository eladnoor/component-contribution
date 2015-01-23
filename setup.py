import sys
from distutils.version import StrictVersion

OK_flag = True

###### PYTHON VERSION #######

if sys.version_info.major != 2:
    sys.stderr.write('ERROR: Please use Python 2 only\n')
    OK_flag = False

if sys.version_info.minor < 7:
    sys.stderr.write('ERROR: Please use Python 2.7 or higher, but not Python 3\n')
    OK_flag = False

try:
    import numpy
    if StrictVersion(numpy.__version__) < StrictVersion('1.6.2'):
        raise ImportError
       
###### NUMPY VERSION #######

except ImportError:
    sys.stderr.write('ERROR: Please install numpy version 1.6.2 or higher\n')
    OK_flag = False

###### OCT2PY VERSION #######

try:
    import oct2py
    if StrictVersion(oct2py.__version__) < StrictVersion('2.4.0'):
        raise ImportError
       
except ImportError:
    sys.stderr.write('WARNING: Please install oct2py version 2.4.0 or higher\n')
    OK_flag = False
    
###### OCTAVE VERSION #######

octave_ver = oct2py.octave.eval('OCTAVE_VERSION()', verbose=False)
if StrictVersion(octave_ver) < StrictVersion('3.6.4'):
    sys.stderr.write('WARNING: Please install Octave version 3.6.4 or higher\n')
    OK_flag = False


if OK_flag:
    sys.stderr.write('Success, all required packages are installed and up-to-date\n')
else:
    sys.stderr.write('Failure, some of the required packages are missing or have an older version\n')
