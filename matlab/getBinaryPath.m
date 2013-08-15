function [babel_bin, cxcalc_bin] = getBinaryPath()
if ispc % Windows
    babel_bin = 'babel';
    cxcalc_bin = 'cxcalc';
elseif ismac % OS X
    babel_bin = '/opt/local/bin/babel';
    cxcalc_bin = '/Applications/ChemAxon/MarvinBeans/bin/cxcalc';
else % Linux
    babel_bin = 'babel';
    cxcalc_bin = 'cxcalc';
end
