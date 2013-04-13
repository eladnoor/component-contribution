function model = addBiochemicalEnergies(model, pH, I, T)
if nargin < 2
    pH = 7.0;
end
if nargin < 3
    I = 0.25;
end
if nargin < 4
    T = 298.15;
end

R = 8.31e-3; % kJ/mol/K
F = 96.45; % kC/mol

fprintf('Performing forward transform\n');

forward_ddG0 = zeros(length(model.mets), 1);
model.pH = pH;
model.I = I;
model.T = T;

model.DfG0_pseudoisomers = [];
for i = 1:length(model.mets)
    diss = model.pKa(i);
    dG0s = cumsum(-[0, diag(diss.pKas, 1)'] * R * model.T * log(10));
    dG0s = dG0s - dG0s(diss.majorMSpH7);
    pseudoisomers = [dG0s(:), diss.nHs(:), diss.zs(:)];
    model.DfG0_pseudoisomers = [model.DfG0_pseudoisomers; ...
                                i * ones(size(pseudoisomers, 1), 1), ...
                                pseudoisomers(:, 1) + model.DfG0(i), ...
                                pseudoisomers(:, 2), ...
                                pseudoisomers(:, 3)];
    if isempty(pseudoisomers)
        continue
    end
    ddG0 = Transform(pseudoisomers, model.pH, model.I, model.T);
    forward_ddG0(i) = ddG0;
end

model.DfG0_prime = model.DfG0 + forward_ddG0;
model.DrG0 = full(model.S') * model.DfG0;
S_noH = model.S;
S_noH(strmatch('h[',model.mets),:) = 0;
model.DrG0_prime = full(S_noH') * model.DfG0_prime;
model.u_DrG0 = sqrt(diag(S_noH' * model.covf * S_noH));
model.u_DrG0(sum(model.S~=0)==1) = 1e10; % set uncertainty of transport reactions to inf


