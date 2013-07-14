% Initialize G matrix, and then use the python script "inchi2gv.py" to decompose each of the 
% compounds that has an InChI and save the decomposition as a row in the G matrix.
function training_data = createKeggGroupIncidenceMatrix(model, training_data)

fprintf('Creating group incidence matrix\n');

%%
% first just run the script to get the list of group names
[~,groups] = system(['python ' getBasePath() 'python/inchi2gv.py -l']);
groups = regexp(groups,'\n','split');
Groups = groups(~cellfun(@isempty, groups));

[~, missing_inds] = setdiff(model.cids, training_data.cids);
training_data.cids_orig = training_data.cids;
training_data.cids = [training_data.cids, model.cids(missing_inds)];
training_data.inchis = [training_data.nstd_inchi, model.inchi.nonstandard(missing_inds)];
[~, inds, ~] = intersect(training_data.cids, model.cids);
training_data.Model2TrainingMap = inds;

training_data.S = [training_data.S; zeros(length(missing_inds), size(training_data.S, 2))];
G = zeros(length(training_data.cids), length(Groups));
training_data.has_gv = true(size(training_data.cids));

%%
% Decompose the compounds in the training_data and add to G
for i = 1:length(training_data.cids)
    inchi = training_data.inchis{i};

    % For compounds that have no InChI, add a unique 1 in a new column
    group_def = getGroupVectorFromInchi(inchi);
    if length(group_def) == length(Groups) % decompoisition successful, place it in G
        G(i, 1:length(group_def)) = group_def;
        training_data.has_gv(i) = true;
    elseif isempty(group_def) % decomposition failed, add a unique 1 in a new column
        G(:, end+1) = 0; 
        G(i, end) = 1;
        training_data.has_gv(i) = false;
    else
        fprintf('InChI = %s\n', inchi);
        fprintf('*************\n%s\n', getGroupVectorFromInchi(inchi, false));
        error('ERROR: while trying to decompose compound C%05d', training_data.cids(i));
    end
end

%%
training_data.G = sparse(G);


