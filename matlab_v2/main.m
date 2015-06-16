% run a python script that generates all the required training data
% for training CC

[~, ~, python_bin] = getBinaryPath();
base_path = getBasePath();
training_data_file = [base_path 'cache/training_data_matlab.mat'];
cc_file = [base_path 'cache/component_contribution_matlab.mat'];

if ~exist(cc_file, 'file')
    if ~exist(training_data_file, 'file')
        fprintf('Preparing training data and writing to: %s\n', training_data_file);
        cmd = [python_bin ' prepare_training_data.py ' training_data_file];
        errcode = system(cmd);
        if (errcode ~= 0) || (~exist(training_data_file, 'file'))
            error('Could not prepare the training data for CC')
        end
    else
        fprintf('Found cached trainind data at: %s\n', training_data_file);
    end
    fprintf('Loading cached training data from: %s\n', training_data_file);
    training_data = load(training_data_file);
    params = componentContribution(training_data);
    fprintf('Writing CC parameters to: %s\n', cc_file);
    save(cc_file, 'params');
else
    fprintf('Found precalculated CC parameters at: %s\n', training_data_file);
    load(cc_file);
end

