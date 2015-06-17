orig_dir = cd();
cd('../matlab_v2');

base_path = getBasePath();
[~, ~, python_bin] = getBinaryPath();
rxn_file = [orig_dir '/wolf_reactions.txt'];
model_file = [orig_dir '/wolf.mat'];
[cc_params, training_data_file, cc_file] = initComponentContribution();

cmd = sprintf('%s prepare_model.py --train_file %s --rxn_file %s --out_file %s', ...
              python_bin, training_data_file, rxn_file, model_file);
errcode = system(cmd);
if (errcode ~= 0) || (~exist(model_file, 'file'))
    error('Could not prepare the model data for CC')
end

model = load(model_file, '-mat');
cc = load(cc_file);

G = model.G;
X = model.X;
v_r = cc.params.preprocess_v_r;
v_g = cc.params.preprocess_v_g;
C1  = cc.params.preprocess_C1;
C2  = cc.params.preprocess_C2;
C3  = cc.params.preprocess_C3;

dG0_cc = X' * v_r + G' * v_g;
U = X' * C1 * X + X' * C2 * G + G' * C2' * X + G' * C3 * G;
dG0_std = sqrt(diag(U));

for i = 1:length(dG0_cc)
    fprintf('---------------------------------------------------\n');
    fprintf('dG0  = %8.1f +- %5.1f\n', dG0_cc(i), dG0_std(i) * 1.96);
end
delete(model_file);
cd(orig_dir);