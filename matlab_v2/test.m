clear all;
py_data = load('../cache/training_data.mat');

%%
if false
    training_data = loadTrainingData(true, 1); % using cached pKa values for KEGG compounds
    training_data = createKeggGroupIncidenceMatrix(training_data);
    training_data = reverseTransformTrainingData(training_data);
    save('../cache/training_data_matlab.mat', 'training_data');
else
    training_data = load('../cache/training_data_matlab.mat');
    training_data = training_data.training_data;
end
