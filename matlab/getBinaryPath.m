function [babel_bin, cxcalc_bin, python_bin] = getBinaryPath()
if ispc % Windows
    babel_bin = 'babel';
    cxcalc_bin = 'cxcalc';
    python_bin = 'python';
elseif ismac % OS X
    setenv('DYLD_LIBRARY_PATH', '/usr/X11');
    babel_bin = '/opt/local/bin/babel';
    cxcalc_bin = '/Applications/ChemAxon/MarvinBeans/bin/cxcalc';
    python_bin = '/usr/local/bin/python';
else % Linux
    babel_bin = 'babel';
    cxcalc_bin = 'cxcalc';
    python_bin = 'python';
end
