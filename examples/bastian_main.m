cd ../matlab

REACTION_FNAME = '../examples/bastian_reactions.txt';
CACHED_TRAINING_DATA_FNAME = '../cache/training.mat';
CACHED_MODEL_FNAME = '../cache/bastian_model.mat';
CACHED_RT_TRAINING_DATA_FNAME = '../cache/bastian_rt_training.mat';
THERMO_MODEL_FNAME = '../cache/bastian_model_thermo.mat';
mkdir('../cache');

%%
% load the model and get the relevant pKas
if exist(CACHED_MODEL_FNAME, 'file')
    fprintf('Loading model from cache: %s\n', CACHED_MODEL_FNAME);
    load(CACHED_MODEL_FNAME);
else
    model = loadKeggModel(REACTION_FNAME);
    save(CACHED_MODEL_FNAME, 'model', '-v7');
end

%%
if exist(CACHED_RT_TRAINING_DATA_FNAME, 'file')
    fprintf('Loading reverse-transformed training data from cache: %s\n', CACHED_RT_TRAINING_DATA_FNAME);
    load(CACHED_RT_TRAINING_DATA_FNAME);
else
%%
    % load the training data (NIST, Alberty, redox) and get the relevant pKas
    if exist(CACHED_TRAINING_DATA_FNAME, 'file')
        fprintf('Loading training data from cache: %s\n', CACHED_TRAINING_DATA_FNAME);
        load(CACHED_TRAINING_DATA_FNAME);
    else
        use_cached_kegg_inchis = true;
        use_cached_kegg_pkas = true;
        formation_weight = 1;
        training_data = loadTrainingData(use_cached_kegg_inchis, use_cached_kegg_pkas, formation_weight);
        save(CACHED_TRAINING_DATA_FNAME, 'training_data', '-v7');
    end

    % match between the compounds in the model and the KEGG IDs used in the
    % training data, and create the group incidence matrix (G) for the
    % combined set of all compounds.
    rt_training_data = createKeggGroupIncidenceMatrix(model, training_data);

    % apply the reverse Legendre transform for the relevant training observations (typically
    % apparent reaction Keq from TECRDB)
    rt_training_data = reverseTransformTrainingData(model, rt_training_data);
    save(CACHED_RT_TRAINING_DATA_FNAME, 'rt_training_data', '-v7');
end

%%
% perform the component contribution method
fprintf('Adding thermodynamic information to the model (using CC)\n');
model = addThermoToModel(model, rt_training_data);

%%
% perform the component contribution method
fprintf('Adding biochemical (transformed) energies to the model\n');
model = addBiochemicalEnergies(model, 7, 0.1, 300);

%%
fprintf('Writing thermodynamic model to: %s\n', THERMO_MODEL_FNAME);
save(THERMO_MODEL_FNAME, 'model', '-v7');
