% A script for estimating the quality of our reverse Legendre transform.
% The method is based on looking at the standard deviations of the
% reactions with multiple examples before and after the transform. The
% expected result would be a significantly lower stdev.

CACHED_RT_TRAINING_DATA_FNAME = 'bastian_data/rt_training.mat';
load(CACHED_RT_TRAINING_DATA_FNAME);

succ = cell2mat({rt_training_data.kegg_pKa.success});
% use only reactions which involve compounds that we have successfully
% found their pKa values
inds = find(sum(abs(rt_training_data.S(succ==0, :)), 1) == 0);
S = rt_training_data.S(:, inds);
b = rt_training_data.dG0(inds);
b_prime = rt_training_data.dG0_prime(inds);

[S_uniq, ~, inds] = unique(S', 'rows');
S_uniq = S_uniq';
M_repeats = spconvert([(1:length(inds))', inds, ones(length(inds),1)])';
n_repeats = sum(M_repeats, 2);

mu = (M_repeats * b) ./ n_repeats;
sigma = sqrt((M_repeats * b.^2) ./ n_repeats - mu.^2);

mu_prime = (M_repeats * b_prime) ./ n_repeats;
sigma_prime = sqrt((M_repeats * b_prime.^2) ./ n_repeats - mu_prime.^2);

%%
close all;
if false
    figure(1);

    r = 0:0.33333:10;
    inds = find(n_repeats > 5 & sigma_prime > 0 & sigma > 0);

    cdfplot(sigma(inds), 'b', 'LineWidth', 2);
    hold on;
    cdfplot(sigma_prime(inds), 'r', 'LineWidth', 2);

    mu1 = median(sigma(inds));
    mu2 = median(sigma_prime(inds));

    text(mu1, 0.2, sprintf('%.1f kJ/mol', mu1), 'FontSize', 12, 'EdgeColor', 'b');
    plot([mu1, mu1], [0, 0.5], 'b-');
    text(mu2, 0.1, sprintf('%.1f kJ/mol', mu2), 'FontSize', 12, 'EdgeColor', 'r');
    plot([mu2, mu2], [0, 0.5], 'r-');
    h = legend('{\sigma} ({\Delta} G^{\circ})', '{\sigma} ({\Delta} G''^{\circ})');
    title('Standard deviations of {\Delta} G^{\circ} and {\Delta} G''^{\circ} for reactions with >5 examples');
    xlabel('standard deviation in kJ/mol');
    ylabel('cumulative distribution function');

    export_fig('rt_cdf', '-eps');
end
%%
lambda = 1e-6; % regularization term
[m, n] = size(S_uniq);
cov_G_rc = (S_uniq * S_uniq' + lambda * eye(m));
G_rc = cov_G_rc \ S_uniq * mu;
e_rc = (S_uniq' * G_rc - mu);

if false
    figure(2);
    hist(e_rc);
end

inds = find(abs(e_rc) > 10);
r = full(S_uniq(:,inds));
e_rc(inds(1))
[rt_training_data.cids(r(:,1) ~= 0); r(find(r(:,1)),1)']
e_rc(inds(2))
[rt_training_data.cids(r(:,2) ~= 0); r(find(r(:,2)),2)']

