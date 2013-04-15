cd ../matlab;
[S, cids] = loadKeggModel('../examples/wolf_reactions.txt');
[DfG0_prime, covf] = getGibbsForKeggModel(S, cids);