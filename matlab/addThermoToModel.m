% given a standard COBRA model, add thermodynamic data to it using
% the Component Contribution method
%
% inputs:
%   model        - a standard COBRA model
%
% returns:
%   model        - same as the input, but with values of Gibbs energies
%
function model = addThermoToModel(model, rt_training_data)

S = full(rt_training_data.S);
G = full(rt_training_data.G);
b = rt_training_data.dG0;
w = rt_training_data.weights;
model2train_map = rt_training_data.Model2TrainingMap;

%%
fprintf('Running Component Contribution method\n');
% Estimate standard Gibbs energies of formation
[x, cov_x] = componentContribution(S, G, b, w);

%%
% Map estimates back to model
model.DfG0 = x(model2train_map);
model.covf = cov_x(model2train_map, model2train_map);
model.uf = diag(sqrt(model.covf));
model.DrG0 = model.S' * model.DfG0;
model.ur = diag(sqrt(model.S'*model.covf*model.S));
model.uf(model.uf >= 1e3) = 1e10; % Set large uncertainty in formation energies to inf
model.ur(model.ur >= 1e3) = 1e10; % Set large uncertainty in reaction energies to inf
model.ur(sum(model.S~=0)==1) = 1e10; % set uncertainty of exchange, demand and sink reactions to inf


