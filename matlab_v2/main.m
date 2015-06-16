% run a python script that generates all the required training data
% for training CC

[~, ~, python_bin] = getBinaryPath();
base_path = getBasePath();
training_data_file = [base_path 'cache/training_data_matlab.mat'];
cc_file = [base_path 'cache/component_contribution_matlab.mat'];

if ~exist(training_data_file, 'file')
    cmd = [python_bin ' prepare_training_data.py ' training_data_file];
    system(cmd);
end

if ~exist(cc_file, 'file')
    training_data = load(training_data_file);
    params = componentContribution(training_data);
    save(cc_file, 'params');
else
    load(cc_file);
end