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
model.mets = cell(length(model.cids),1);
for i = 1:length(model.cids)
    model.mets{i} = sprintf('C%05d[c]', model.cids(i));
end

fprintf('Loaded a KEGG model with %d compounds and %d reactions\n', size(model.S, 1), size(model.S, 2));

%%
% get the InChIs for all the compounds in the training data
% (note that all of them have KEGG IDs)
kegg_inchies = getInchies(model.cids, false);

% match the kegg data to the model.cids
inds = subset_index(kegg_inchies.cids, model.cids);

model.inchi.standard = kegg_inchies.std_inchi(inds);
model.inchi.standardWithStereo = kegg_inchies.std_inchi_stereo(inds);
model.inchi.standardWithStereoAndCharge = kegg_inchies.std_inchi_stereo_charge(inds);
model.inchi.nonstandard = kegg_inchies.nstd_inchi(inds);
    
%%

model.pKa = getKeggpKas(model.cids, model.inchi.nonstandard);