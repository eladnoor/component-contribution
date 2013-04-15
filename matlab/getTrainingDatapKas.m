function kegg_pKa = getTrainingDatapKas(target_cids, target_inchies, use_cache)

if nargin < 3
    use_cache = true;
end

CACHED_KEGG_PKA_MAT_FNAME = [getBasePath() 'cache/kegg_pkas.mat'];

% Load pKas from cache if required
if use_cache && exist(CACHED_KEGG_PKA_MAT_FNAME, 'file')
    fprintf('Loading the pKa values for the training data from: %s\n', CACHED_KEGG_PKA_MAT_FNAME);
    load(CACHED_KEGG_PKA_MAT_FNAME);
    [~, I] = setdiff(target_cids, cell2mat({kegg_pKa.cid}));
    target_cids = target_cids(I);
    target_inchies = target_inchies(I);
else
    kegg_pKa = [];
end

fprintf('Calculating the pKa values for the training data using ChemAxon')
kegg_pKa = [kegg_pKa; getKeggpKas(target_cids, target_inchies)];
if use_cache
    save(CACHED_KEGG_PKA_MAT_FNAME, 'kegg_pKa', '-v7');
end
