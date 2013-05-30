function inchies = getInchies(target_cids, use_cache, use_kegg_additions)
if nargin < 2
    use_cache = true;
end
if nargin < 3
    use_kegg_additions = true;
end

[s, mess, messid] = mkdir([getBasePath() 'cache'], 'newFolder');
KEGG_ADDITIONS_TSV_FNAME = [getBasePath() 'data/kegg_additions.tsv'];
CACHED_KEGG_INCHI_MAT_FNAME = [getBasePath() 'cache/kegg_inchies.mat'];

% Load relevant InChIs (for all compounds in the training data)
if use_cache && exist(CACHED_KEGG_INCHI_MAT_FNAME, 'file')
    fprintf('Loading the InChIs for the training data from: %s\n', CACHED_KEGG_INCHI_MAT_FNAME);
    load(CACHED_KEGG_INCHI_MAT_FNAME);
else
    inchies.cids = [];
    inchies.std_inchi = {};
    inchies.std_inchi_stereo = {};
    inchies.std_inchi_stereo_charge = {};
    inchies.nstd_inchi = {};
    if use_kegg_additions
        % load the InChIs for all KEGG compounds in the 'kegg_additions.tsv' file.
        % this contains a few corrections needed in KEGG and added compounds (all starting with C80000)
        if ~exist(KEGG_ADDITIONS_TSV_FNAME, 'file')
            error(['file not found: ', KEGG_ADDITIONS_TSV_FNAME]);
        end
        fprintf('Adding the non-KEGG compounds first from: %s\n', KEGG_ADDITIONS_TSV_FNAME);

        fid = fopen(KEGG_ADDITIONS_TSV_FNAME, 'r');
        fgetl(fid); % fields are: name, cid, inchi
        filecols = textscan(fid, '%s%d%s', 'delimiter','\t');
        fclose(fid);
        added_cids = filecols{2};
        added_inchies = filecols{3};
        for i = 1:length(added_cids)
            inchies.cids(i) = added_cids(i);
            [std_inchi, std_inchi_stereo, std_inchi_stereo_charge, nstd_inchi] = getInchi(added_cids(i), added_inchies{i,1});
            inchies.std_inchi{1,i} = std_inchi;
            inchies.std_inchi_stereo{1,i} = std_inchi_stereo;
            inchies.std_inchi_stereo_charge{1,i} = std_inchi_stereo_charge;
            inchies.nstd_inchi{1,i} = nstd_inchi;
        end
    end    
end

target_cids = setdiff(target_cids, inchies.cids);
if ~isempty(target_cids)

    fprintf('Obtaining MOL files for the training data from KEGG and converting to InChI using OpenBabel\n');

    n = length(inchies.cids);
    for i = 1:length(target_cids)
        inchies.cids(n+i) = target_cids(i);
        
        % gets the MOL from KEGG and converts it to all 4 InChI types
        [std_inchi, std_inchi_stereo, std_inchi_stereo_charge, nstd_inchi] = getInchi(target_cids(i));
        inchies.std_inchi{1,n+i} = std_inchi;
        inchies.std_inchi_stereo{1,n+i} = std_inchi_stereo;
        inchies.std_inchi_stereo_charge{1,n+i} = std_inchi_stereo_charge;
        inchies.nstd_inchi{1,n+i} = nstd_inchi;
    end
    if use_cache
        save(CACHED_KEGG_INCHI_MAT_FNAME, 'inchies', '-v7');
    end
end
