cd ../matlab;
[S, cids] = loadKeggModel('../examples/wolf_reactions.txt');
[DfG0_prime, covf] = getGibbsForKeggModel(S, cids);
res.G0 = DfG0_prime;
res.cov = covf;
save('../examples/wolf_res.mat', 'res', '-v7');
