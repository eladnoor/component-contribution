function [babel_bin, cxcalc_bin, python_bin] = getBinaryPath()
if ispc % Windows
    babel_bin = 'babel';
    cxcalc_bin = 'cxcalc';
    python_bin = 'python';
elseif ismac % OS X
    babel_bin = '/opt/local/bin/babel';
    cxcalc_bin = '/Applications/ChemAxon/MarvinBeans/bin/cxcalc';
    python_bin = 'python';
else % Linux
    babel_bin = 'babel';
    cxcalc_bin = 'cxcalc';
    python_bin = 'python';
end
