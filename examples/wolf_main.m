%%
clear all;
cd ../matlab;
[S, cids] = loadKeggModel('../examples/wolf_reactions.txt');
model = createKeggModelFromS(S, cids);
training_data = loadTrainingData(false, false, 1);
rt_training_data = createKeggGroupIncidenceMatrix(model, training_data);
rt_training_data = reverseTransformTrainingData(model, rt_training_data);
save('../examples/wolf.mat');

%%
clear all;
load('../examples/wolf.mat');
pH = 7.5;
I = 0.2;
T = 298.15;
model = addThermoToModel(model, rt_training_data);
model = addBiochemicalEnergies(model, pH, I, T);

save('../examples/wolf_res.mat');

%%
clear all;
load('../examples/wolf_res.mat');
load('../examples/groups.mat');
plot(model.DrG0, model_dG0, '.');