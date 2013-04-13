function model = loadKeggModel(fname, fmt)

if nargin < 2 || strcmp(fmt, 'bastian')
    arrow = '=';
    has_reaction_ids = false;
elseif strcmp(fmt, 'rienk')
    arrow = '-->';
    has_reaction_ids = true;
else
    error('Kegg Model format must be "bastian" or "rienk"');
end

fid = fopen(fname);
if fid == -1
    error('cannot find the model file: %s', fname);
end
res = textscan(fid, '%s', 'delimiter', '\n');
model.cids = [];
model.reactions = {};
for i = 1:length(res{1})
    curr_line = res{1}{i};
    if has_reaction_ids
        tokens = regexp(curr_line, '(\w+)\s+(.*)', 'tokens');
        curr_line = tokens{1}{2};
    end
    sprs = reaction2sparse(curr_line, arrow);
    model.cids = unique([model.cids, find(sprs)]);
    model.reactions = [model.reactions, {sprs}];
end
model.mets = cell(length(model.cids),1);
for i = 1:length(model.cids)
    model.mets{i} = sprintf('C%05d[c]', model.cids(i));
end
model.S = zeros(length(model.cids), length(model.reactions));
for i = 1:length(model.reactions)
    r = full(model.reactions{i});
    if norm(r) ~= 0
        model.S(ismember(model.cids, find(r)), i) = r(r ~= 0);
    end
end
fprintf('Loaded a KEGG model with %d compounds and %d reactions\n', size(model.S, 1), size(model.S, 2));

%%
% get the InChIs for all the compounds in the training data
% (note that all of them have KEGG IDs)
kegg_inchies = getInchies(model.cids, false);
inds = ismember(kegg_inchies.cids, model.cids);
model.inchi.standard = kegg_inchies.std_inchi(inds);
model.inchi.standardWithStereo = kegg_inchies.std_inchi_stereo(inds);
model.inchi.standardWithStereoAndCharge = kegg_inchies.std_inchi_stereo_charge(inds);
model.inchi.nonstandard = kegg_inchies.nstd_inchi(inds);
    
%%

model.pKa = getKeggpKas(model.cids, model.inchi.nonstandard);
