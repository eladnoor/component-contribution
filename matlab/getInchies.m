function inchies = getInchies(use_cache)
if nargin < 2
    use_cache = true;
end

[s, mess, messid] = mkdir(getBasePath(), 'cache');
COMPOUNDS_TSV_FNAME = [getBasePath() 'data/compounds.tsv'];
COMPOUNDS_CACHE_FNAME = [getBasePath() 'cache/compounds.mat'];

% Load relevant InChIs (for all compounds in the training data)
if use_cache && exist(COMPOUNDS_CACHE_FNAME, 'file')
    fprintf('Loading the InChIs from cache file: %s\n', COMPOUNDS_CACHE_FNAME);
    load(COMPOUNDS_CACHE_FNAME);
else
    inchies.cids = {};
    inchies.std_inchi = {};
    inchies.std_inchi_stereo = {};
    inchies.std_inchi_stereo_charge = {};
    inchies.nstd_inchi = {};
    
    % load the InChIs for all compounds in 'compounds.tsv'
    % which contains all the compounds in KEGG, with a few corrections and
    % several dozens of added compounds (all starting with C80000)
    if ~exist(COMPOUNDS_TSV_FNAME, 'file')
        error(['file not found: ', COMPOUNDS_TSV_FNAME]);
    end
    fprintf('Adding the non-KEGG compounds first from: %s\n', ...
            COMPOUNDS_TSV_FNAME);

    fid = fopen(COMPOUNDS_TSV_FNAME, 'r');
    fgetl(fid); % fields are: cid, name, inchi
    filecols = textscan(fid, '%s%s%s', 'delimiter','\t');
    fclose(fid);
    compound_ids = filecols{1};
    compound_inchies = filecols{3};
    for i = 1:length(compound_ids)
        inchies.cids{i} = compound_ids{i,1};
        [std_inchi, std_inchi_stereo, std_inchi_stereo_charge, nstd_inchi] = ...
            getInchi(compound_ids{i,1}, compound_inchies{i,1});
        inchies.std_inchi{i} = std_inchi;
        inchies.std_inchi_stereo{i} = std_inchi_stereo;
        inchies.std_inchi_stereo_charge{i} = std_inchi_stereo_charge;
        inchies.nstd_inchi{i} = nstd_inchi;
    end

    if use_cache
        save(COMPOUNDS_CACHE_FNAME, 'inchies', '-v7');
    end
end
