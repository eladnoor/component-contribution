% Generates the structure that contains all the training data needed for
% Component Contribution.
%
% Input:
%   formation_weight - the relative weight to give the formation energies (Alberty's data)
%                      compared to the reaction measurements (TECRDB)
function training_data = loadTrainingData(use_cached_kegg_pkas, formation_weight)
if nargin < 1
    use_cached_kegg_pkas = true;
end
if nargin < 2
    formation_weight = 1;
end

TECRDB_TSV_FNAME = [getBasePath() 'data/TECRDB.tsv'];
FORMATION_TSV_FNAME = [getBasePath() 'data/formation_energies_transformed.tsv'];
REDOX_TSV_FNAME = [getBasePath() 'data/redox.tsv'];

WEIGHT_TECRDB = 1;
WEIGHT_FORMATION = formation_weight;
WEIGHT_REDOX = formation_weight;
R = 8.31e-3; % kJ/mol/K
F = 96.485; % kC/mol

if ~exist(TECRDB_TSV_FNAME, 'file')
    error(['file not found: ', TECRDB_TSV_FNAME]);
end

if ~exist(FORMATION_TSV_FNAME, 'file')
    error(['file not found: ', FORMATION_TSV_FNAME]);
end

if ~exist(REDOX_TSV_FNAME, 'file')
    error(['file not found: ', REDOX_TSV_FNAME]);
end

% Read the raw data of TECRDB (NIST)
reactions = {};
cids = {};
cids_that_dont_decompose = {};
thermo_params = []; % columns are: dG'0, T, I, pH, pMg, weight, balance?

fid = fopen(TECRDB_TSV_FNAME, 'r');
% fields are: URL, REF_ID, METHOD, EVAL, EC, ENZYME NAME, REACTION IN
% KEGG IDS, REACTION IN COMPOUND NAMES, K, K', T, I, pH, pMg
res = textscan(fid, '%s%s%s%s%s%s%s%s%f%f%f%f%f%f', 'delimiter','\t');
fclose(fid);

inds = find(~isnan(res{10}) .* ~isnan(res{11}) .* ~isnan(res{13}));

dG0_prime = -R * res{11}(inds) .* log(res{10}(inds)); % calculate dG'0
thermo_params = [dG0_prime, res{11}(inds), res{12}(inds), res{13}(inds), ...
                 res{14}(inds), WEIGHT_TECRDB * ones(size(inds)), ...
                 true(size(inds))];

% parse the reactions in each row
for i = 1:length(inds)
    sprs = reaction2sparse(res{7}{inds(i)});
    reactions = [reactions, {sprs}];
end
fprintf('Successfully added %d values from TECRDB\n', length(inds));

% Read the Formation Energy data.
fid = fopen(FORMATION_TSV_FNAME, 'r');
fgetl(fid); % skip the first header line
% fields are: cid, name, dG'0, pH, I, pMg, T, decompose?, compound_ref, remark
res = textscan(fid, '%s%s%f%f%f%f%f%f%s%s', 'delimiter','\t');
fclose(fid);

inds = find(~isnan(res{3}));
thermo_params = [thermo_params; [res{3}(inds), res{7}(inds), res{5}(inds), ...
                                 res{4}(inds), res{6}(inds), ...
                                 WEIGHT_FORMATION * ones(size(inds)), ...
                                 false(size(inds))]];
for i = 1:length(inds)
    sprs = {};
    sprs{1,1} = res{1}{inds(i)};
    sprs{1,2} = 1;
    reactions = [reactions, {sprs}];
end

cids_that_dont_decompose = res{1}(find(res{8} == 0));

fprintf('Successfully added %d formation energies\n', length(res{1}));

% Read the Reduction potential data.
fid = fopen(REDOX_TSV_FNAME, 'r');
fgetl(fid); % skip the first header line
% fields are: name, CID_ox, nH_ox, charge_ox, CID_red, nH_red, 
%             charge_red, E'0, pH, I, pMg, T, ref
res = textscan(fid, '%s%s%f%f%s%f%f%f%f%f%f%f%s', 'delimiter', '\t');
fclose(fid);

delta_e = (res{6} - res{3}) - (res{7} - res{4}); % delta_nH - delta_charge
dG0_prime = -F * res{8} .* delta_e;
thermo_params = [thermo_params; [dG0_prime, res{12}, res{10}, res{9}, ...
                                 res{11}, ...
                                 WEIGHT_REDOX * ones(size(dG0_prime)), ...
                                 false(size(dG0_prime))]];

for i = 1:length(res{1})
    sprs = {};
    sprs{1,1} = res{2}{i};
    sprs{1,2} = -1;
    sprs{2,1} = res{5}{i};
    sprs{2,2} = 1;
    reactions = [reactions, {sprs}];
end

fprintf('Successfully added %d redox potentials\n', length(res{1}));

% convert the list of reactions in sparse notation into a full
% stoichiometric matrix, where the rows (compounds) are according to the
% CID list 'cids'.
cids = {};
for i = 1:length(reactions)
    r = reactions{i};
    cids = unique({cids{:}, r{:,1}});
end

S = zeros(length(cids), length(reactions));
for i = 1:length(reactions)
    r = reactions{i};
    [~, idx] = ismember({r{:,1}}, cids);
    S(idx, i) = [r{:,2}];
end
    
training_data.cids = cids;
training_data.cids_that_dont_decompose = cids_that_dont_decompose;

training_data.S = sparse(S);
training_data.dG0_prime = thermo_params(:, 1);
training_data.T = thermo_params(:, 2);
training_data.I = thermo_params(:, 3);
training_data.pH = thermo_params(:, 4);
training_data.pMg = thermo_params(:, 5);
training_data.weights = thermo_params(:, 6);
training_data.balance = thermo_params(:, 7);

% get the InChIs for all the compounds in the training data
% (note that all of them have KEGG IDs)
kegg_inchies = getInchies();
[~, inds, ~] = intersect(kegg_inchies.cids, training_data.cids);
assert(length(inds) == length(training_data.cids));
training_data.std_inchi = kegg_inchies.std_inchi(inds);
training_data.std_inchi_stereo = kegg_inchies.std_inchi_stereo(inds);
training_data.std_inchi_stereo_charge = kegg_inchies.std_inchi_stereo_charge(inds);
training_data.nstd_inchi = kegg_inchies.nstd_inchi(inds);

% use the chemical formulas from the InChIs to verify that each and every
% reaction is balanced.
training_data = balanceReactionsInTrainingData(training_data);

% get the pKas for the compounds in the training data (using ChemAxon)
kegg_pKa = getTrainingDatapKas(training_data.cids, training_data.nstd_inchi, use_cached_kegg_pkas);
[~, inds, ~] = intersect(training_data.cids, {kegg_pKa.cid});
assert(length(inds) == length(training_data.cids));
training_data.kegg_pKa = kegg_pKa(inds);

