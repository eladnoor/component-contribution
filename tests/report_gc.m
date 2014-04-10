%%
clear all;
cd ../matlab;

if false
    [S, cids] = loadKeggModel('../tests/report_gc_reactions.txt');
    model = createKeggModelFromS(S, cids);
    training_data = loadTrainingData(false, false, 1);
    rt_training_data = createKeggGroupIncidenceMatrix(model, training_data);
    rt_training_data = reverseTransformTrainingData(model, rt_training_data);
    save('../tests/report_gc.mat');
else
    load('../tests/report_gc.mat');
end

S = full(rt_training_data.S);
G = full(rt_training_data.G);
b = rt_training_data.dG0;
w = rt_training_data.weights;
#[x, cov_x, params] = componentContribution(S, G, b, w);

[~,groups] = system(['python ' getBasePath() 'python/inchi2gv.py -l']);
groups = regexp(groups, '\n', 'split');
for i = 1:size(groups, 2)
    fprintf('%s - %d examples\n', groups{1, i}, sum(G(:, i) ~= 0));
end
