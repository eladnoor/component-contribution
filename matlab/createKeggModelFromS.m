% model = createKeggModelFromS(S, cids)
%
% returns a model with the fields required for running CC
%
% Arguments:
%   S    - a stoichiometric matrix (n by m)
%   cids - the KEGG compound IDs in the same order as in S (1 by n)

function model = createKeggModelFromS(S, cids)

model.S = S;
model.cids = cids;
model.mets = cell(size(model.cids));
for i = 1:length(model.cids)
    model.mets{i} = [model.cids{i} '[c]'];
end

%%
% get the InChIs for all the compounds in the training data
% (note that all of them have KEGG IDs)
inchies = getInchies(false);

% match the kegg data to the model.cids
[~, inds] = ismember(model.cids, inchies.cids);
if any(inds == 0) 
    error('model.cids is not a proper subset of kegg_inchies.cids');
end

model.inchi.standard = inchies.std_inchi(inds);
model.inchi.standardWithStereo = inchies.std_inchi_stereo(inds);
model.inchi.standardWithStereoAndCharge = inchies.std_inchi_stereo_charge(inds);
model.inchi.nonstandard = inchies.nstd_inchi(inds);
    
%%

model.pKa = getKeggpKas(model.cids, model.inchi.nonstandard);
