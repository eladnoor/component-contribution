function [formula, nH, charge] = getFormulaAndChargeFromInChI(inchi)
% Extracts chemical formula, number of hydrogen atoms and charge of a
% metabolite from its InChI. Depends on ChemAxon's Calculator Plugins
% (cxcalc).
% 
% [formula,charge] = getFormulaAndChargeFromInChI(inchi)
% 
% INPUT
% inchi.......IUPAC InChI for a particular pseudoisomer of a metabolite
% 
% OUTPUTS
% formula....The chemical formula for the input pseudoisomer
% nH.........The number of hydrogen atoms in formula
% charge.....The charge on the input pseudoisomer

formula = '';
charge = nan;
nH = nan;

cxcalc_cmd = 'cxcalc';
if ispc
    cmd = [cxcalc_cmd ' formula formalcharge "' inchi '"'];
else
    cmd = ['echo "' inchi '" | ' cxcalc_cmd ' formula formalcharge'];
end
[s,r] = system(cmd);

pat1 = 'id\tFormula\tFormal charge\n\d+\t(?<formula>[^\t]+)\t(?<charge>[\d\-]+)';
if s == 0 && ~isempty(regexp(r, pat1, 'once'))
    fields = regexp(r, pat1, 'names');
    formula = fields.formula;
    charge = str2double(fields.charge);
end

if ~isempty(formula)
    pat2 = '(?<e>[A-Z][a-z]*)(?<n>\d*)';
    u = regexp(formula, pat2, 'names');
    u = {u.e;u.n};
    ind = find(strcmp(u(1,:), 'H'));
    if ~isempty(ind)
        if ~isempty(u{2,ind})
            nH = str2double(u(2,ind));
        else
            nH = 1;
        end
    else
        nH = 0;
    end
end
