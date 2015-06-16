% run a python script that generates all the required training data
% for training CC
clear all;
[~, ~, python_bin] = getBinaryPath();
base_path = getBasePath();
training_data_file = [base_path 'cache/training_data_matlab.mat'];
cc_file = [base_path 'cache/component_contribution_matlab.mat'];

if ~exist(training_data_file, 'file')
    disp('The training data .mat file is not cached, preparing data now...');
    cmd = [python_bin ' prepare_training_data.py ' training_data_file];
    errcode = system(cmd);
    if (errcode ~= 0) || (~exist(training_data_file, 'file'))
        error('Could not prepare the training data for CC')
    end
end

if ~exist(cc_file, 'file')
    disp('Using the cached training data to train CC now...');
    training_data = load(training_data_file);
    params = componentContribution(training_data);
    save(cc_file, 'params');
else
    load(cc_file);
end

