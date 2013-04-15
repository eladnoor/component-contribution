function [ DfG0_prime, covf ] = getGibbsForKeggModel(S, cids, pH, I, T)
if nargin < 5
    T = 300;
end
if nargin < 4
    I = 0.2;
end
if nargin < 3
    pH = 7.5;
end

%GETGIBBSFORKEGGMODEL returns the CC estimates for a model based on KEGG
%   Arguments:
%        S    - a stoichiometric matrix (n by m)
%        cids - the KEGG compound IDs in the same order as in S (1 by n)

model = createModelFromS(S, cids);
training_data = loadTrainingData(false, false, 1);
rt_training_data = createKeggGroupIncidenceMatrix(model, training_data);
rt_training_data = reverseTransformTrainingData(model, rt_training_data);
model = addThermoToModel(model, rt_training_data);
model = addBiochemicalEnergies(model, pH, I, T);

DfG0_prime = model.DfG0_prime;
covf = model.covf;
