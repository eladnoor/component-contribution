% [ DfG0_prime, covf ] = getGibbsForKeggModel(S, cids, pH, I, T)
%
% computes the CC estimates for a model based on KEGG
%
% Arguments:
%   S          - a stoichiometric matrix (n by m)
%   cids       - the KEGG compound IDs in the same order as in S (1 by n)
%   pH, I, T   - the aqueous conditions for the Legendre transform of Gibbs
%                energies
%
% Return values:
%   DfG0_prime - the estimated standard transformed formation energies
%   covf       - the covariance matrix associated with DfG0_prime 

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


model = createKeggModelFromS(S, cids);
training_data = loadTrainingData(false, false, 1);
rt_training_data = createKeggGroupIncidenceMatrix(model, training_data);
rt_training_data = reverseTransformTrainingData(model, rt_training_data);
model = addThermoToModel(model, rt_training_data);
model = addBiochemicalEnergies(model, pH, I, T);

DfG0_prime = model.DfG0_prime;
covf = model.covf;
