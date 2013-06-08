function model = balanceReactionsS(model,addWater)

if nargin < 2
    addWater = 0;
end



    [MW, Ematrix, elements] = getMolecularWeight(model.inchi.nonstandard, 0);
    model.Ematrix = Ematrix(:, 2:end); % remove H, columns are [C, N, O, P, S, e-]
    %elements = elements(:, 2:end);
    conserved = model.Ematrix' * model.S;
    
    % need to check that all elements are balanced (except H, but including e-)
    % if only O is not balanced, add water molecules

    % check all reactions which can be checked (not NaN) and should be checked
    % (i.e. not formation or redox reactions)
    inds = find(~isnan(conserved(1,:)));% .* model.balance');
    
    if addWater == 1
         % first add water molecules to reactions that need it
        i_h2o = find(model.cids == 1);
        model.S(i_h2o, inds) = model.S(i_h2o, inds) - conserved(3, inds);
        % recalculate conservation matrix
        conserved = model.Ematrix' * model.S;
    end
    
    inds= any(conserved(:, inds));
    
    model.S(:, inds) = 0;
    fprintf('Successfully created balanced S structure: %d compounds and %d reactions\n%d reactions were not balanced and set to zero \n',...
        size(model.S, 1), size(model.S, 2),sum(inds));
    



