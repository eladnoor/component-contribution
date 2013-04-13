function model = addHydrogensAndCharge(model)

fprintf('Storing formation energies of model species\n');

model.majormsCharge = zeros(size(model.mets));
model.nstdInChICharge = zeros(size(model.mets));
model.nstdInChInH = zeros(size(model.mets));
model.nstdInChIFormula = cell(size(model.mets));
for n = 1:length(model.mets)
    nstd_inchi = model.nstdMetInChI{n}; % non-standard InChI
    if ~isempty(nstd_inchi)
        [formula, nH, charge] = getFormulaAndChargeFromInChI(nstd_inchi);
        model.nstdInChIFormula{n} = formula;
        model.nstdInChInH(n) = nH;
        model.nstdInChICharge(n) = charge;
    elseif ~isempty(model.metFormulas{n})
        model.nstdInChIFormula{n} = model.metFormulas{n};
        model.nstdInChInH(n) = getNumAtomsOfElementInFormula(model.metFormulas{n},'H');
        model.nstdInChICharge(n) = model.metCharges(n);
    else
        model.nstdInChIFormula{n} = 'R';
        model.nstdInChInH(n) = 0;
        model.nstdInChICharge(n) = 0;
    end        
    
    met = model.mets{n};
    if isempty(met)
        model.majormsCharge(n) = model.nstdInChICharge(n);
        continue; % this met has no known pKas, assume the charge in the InChI
    end
    met = met(1:end-3);
    k = find(strcmp({model.metSpeciespKa.abbreviation}, met));
    if isempty(k)
        model.majormsCharge(n) = model.nstdInChICharge(n);
        continue; % this met has no known pKas, assume the charge in the InChI
    end
    diss = model.metSpeciespKa(k);
    model.majormsCharge(n) = diss.zs(diss.majorMSpH7);
    model.nstdInChInH(n) = diss.nHs(diss.majorMSpH7);
    if isfield(diss, 'formulas') && ~isempty(diss.formulas{diss.majorMSpH7})
        model.nstdInChIFormula{n} = diss.formulas{diss.majorMSpH7};
    else
        model.nstdInChIFormula{n} = ['H' num2str(diss.nHs(diss.majorMSpH7))];
    end
end
