% Initialize G matrix, and then use the python script "inchi2gv.py" to decompose each of the 
% compounds that has an InChI and save the decomposition as a row in the G matrix.
function training_data = createKeggGroupIncidenceMatrix(training_data)

fprintf('Creating group incidence matrix\n');

%%
% first just run the script to get the list of group names
[~,groups] = system(['python ' getBasePath() 'component_contribution/inchi2gv.py -l']);
groups = regexp(groups,'\n','split');
Groups = groups(~cellfun(@isempty, groups));

G = zeros(length(training_data.cids), length(Groups));
has_gv = true(size(training_data.cids));

%%
% Decompose the compounds in the training_data and add to G
for i = 1:length(training_data.cids)
    inchi = training_data.nstd_inchi{i};

    % For compounds that have no InChI, add a unique 1 in a new column
    group_def = getGroupVectorFromInchi(inchi, true);
    if length(group_def) == length(Groups) % decompoisition successful, place it in G
        G(i, 1:length(group_def)) = group_def;
        has_gv(i) = true;
    elseif isempty(group_def) % decomposition failed, add a unique 1 in a new column
        G(:, end+1) = 0; 
        G(i, end) = 1;
        has_gv(i) = false;
    else
        fprintf('InChI = %s\n', inchi);
        fprintf('*************\n%s\n', getGroupVectorFromInchi(inchi, false));
        error('ERROR: while trying to decompose compound C%05d\n', training_data.cids{i});
    end
end

%%
training_data.G = sparse(G);
training_data.has_gv = has_gv;
