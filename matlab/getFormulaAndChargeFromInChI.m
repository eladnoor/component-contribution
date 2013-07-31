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

[s,r] = system(['cxcalc formula "' inchi '"']);
pat1 = 'id\tFormula\n\d+\t(?<formula>.+)';
if s == 0 && regexp(r,pat1) == 1
    formula = regexp(r,pat1,'names');
    formula = formula.formula;
else
    formula = '';
end

if ~isempty(formula)
    pat2 = '(?<u>[A-Z][a-z]*\d*)';
    u = regexp(formula,pat2,'names');
    u = {u.u}';
    e = regexprep(u,'\d*','');
    n = regexprep(u,'[A-Z_a-z]*','');
    n(cellfun('isempty',n)) = {'1'};
    n = str2double(n);
    if any(strcmp(e,'H'))
        nH = n(strcmp(e,'H'));
    else
        nH = 0;
    end
else
    nH = nan;
end

[s,r] = system(['cxcalc formalcharge "' inchi '"']);
pat3 = 'id\tFormal charge\n\d+\t(?<charge>-?\d+)';
if s == 0 && regexp(r,pat3) == 1
    charge = regexp(r,pat3,'names');
    charge = charge.charge;
    charge = str2double(charge);
else
    charge = nan;
end