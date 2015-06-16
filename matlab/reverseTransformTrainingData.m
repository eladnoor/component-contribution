function training_data = reverseTransformTrainingData(training_data)

if nargin < 3
    use_model_pKas_by_default = true;
end

R = 8.31e-3; % kJ/mol/K

fprintf('Performing reverse transform\n');

% Calculate the reverse transform for all reactions in training_data.
% Note that many of the compounds in the training data are missing from the iAF1260
% model and therefore do not have a BiGG abbreviation or a pKa struct. This
% needs to be fixed somehow.
training_data.reverse_ddG0 = zeros(size(training_data.S, 2), 1);
training_data.I(isnan(training_data.I)) = 0.25; % default ionic strength is 0.25M
training_data.pMg(isnan(training_data.pMg)) = 14; % default pMg is 14
for i = 1:size(training_data.S, 2) % for each reaction in S
    inds = find(training_data.S(:, i));
    reaction_ddG0s = zeros(length(inds), 1);
    for j = 1:length(inds)
        diss = [];

        if inds(j) <= length(training_data.cids)
            % find the diss table from the training data structure
            k = find(strcmp(training_data.cids{inds(j)}, ...
                            {training_data.kegg_pKa.cid}));
            if ~isempty(k)
                diss = training_data.kegg_pKa(k);
            end
        end
        
        if isempty(diss)
            continue;
        end
        
        dG0s = cumsum(-[0, diag(diss.pKas, 1)'] * R * training_data.T(i) * log(10));
        dG0s = dG0s - dG0s(diss.majorMSpH7);
        pseudoisomers = [dG0s(:), diss.nHs(:), diss.zs(:)];
        reaction_ddG0s(j) = Transform(pseudoisomers, training_data.pH(i), training_data.I(i), training_data.T(i));
    end
    training_data.reverse_ddG0(i) = training_data.S(inds, i)' * reaction_ddG0s;
end

training_data.dG0 = training_data.dG0_prime - training_data.reverse_ddG0;
