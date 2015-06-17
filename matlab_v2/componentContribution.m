% Perform the CC method on the training data which comprises
% S - the stoichiometric matrix of measured reactions
% G - the group incidence matrix
% b - the observation vector (standard Gibbs energy of reactions)
% w - the weight vector for each reaction in S
function params = componentContribution(training_data, debug_mode)
if nargin < 2
    debug_mode = false;
end

S = training_data.S;
G = training_data.G;
b = training_data.b;
w = training_data.w;

[m, n] = size(S);
assert (size(G, 1) == m);
assert (size(b, 1) == n);
assert (size(b, 2) == 1);
assert (length(w) == size(S, 2));

% Apply weighing
W = diag(w);
GS = G' * S;

% Linear regression for the reactant layer (aka RC)
[inv_S, r_rc, P_R_rc, P_N_rc] = invertProjection(S * W);

% Linear regression for the group layer (aka GC)
[inv_GS, r_gc, P_R_gc, P_N_gc] = invertProjection(GS * W);

% calculate the group contributions
dG0_gc = inv_GS' * W * b;

% Calculate the contributions in the stoichiometric space
dG0_rc = inv_S' * W * b;
dG0_cc = P_R_rc * dG0_rc + P_N_rc * G * dG0_gc;

% Calculate the residual error (unweighted squared error divided by N - rank)
e_rc = (S' * dG0_rc - b);
MSE_rc = (e_rc' * W * e_rc) / (n - r_rc);

e_gc = (GS' * dG0_gc - b);
MSE_gc = (e_gc' * W * e_gc) / (n - r_gc);

MSE_inf = 1e10;

% Calculate the uncertainty covariance matrices
[inv_SWS, ~, ~, ~] = invertProjection(S*W*S');
[inv_GSWGS, ~, ~, ~] = invertProjection(GS*W*GS');

V_rc = P_R_rc * inv_SWS * P_R_rc;
V_gc  = P_N_rc * G * inv_GSWGS * G' * P_N_rc;
V_inf = P_N_rc * G * P_N_gc * G' * P_N_rc;

% Calculate the total of the contributions and covariances
cov_dG0 = V_rc * MSE_rc + V_gc * MSE_gc + V_inf * MSE_inf;

% preprocessing matrices (for calculating the contribution of each 
% observation)
G1 = P_R_rc * inv_S' * W;
G2 = P_N_rc * G * inv_GS' * W;
G3 = inv_GS' * W;

[~, ~, IC] = unique(S', 'rows', 'sorted');
N_uniq = max(IC);
N_orig = length(IC);
P_col = sparse(1:N_orig, N_uniq + 1 - IC, ones(1, N_orig), N_orig, N_uniq);
S_counter = full(sum(P_col, 1));
preprocess_G1 = G1 * P_col;
preprocess_G2 = G2 * P_col;
preprocess_G3 = G3 * P_col;

% preprocessing matrices (for quick calculation of uncertainty)
preprocess_C1 = cov_dG0;
preprocess_C2 = MSE_gc * P_N_rc * G * inv_GSWGS + MSE_inf * G * P_N_gc;
preprocess_C3 = MSE_gc * inv_GSWGS + MSE_inf * P_N_gc;

% These are the only parameters necessary for using CC later
params.preprocess_v_r = dG0_cc;
params.preprocess_v_g = dG0_gc;
params.preprocess_C1 = preprocess_C1;
params.preprocess_C2 = preprocess_C2;
params.preprocess_C3 = preprocess_C3;

if debug_mode
    % These parameters are not necessary for running CC later,
    % but can be useful for debugging
    params.dG0_rc = dG0_rc;
    params.dG0_gc = dG0_gc;
    params.dG0_cc = dG0_cc;
    params.cov_dG0 = cov_dG0;
    params.V_rc = V_rc;
    params.V_gc = V_gc;
    params.V_inf = V_inf;
    params.MSE_rc = MSE_rc;
    params.MSE_gc = MSE_gc;
    params.MSE_inf = MSE_inf;
    params.P_R_rc = P_R_rc;
    params.P_R_gc = P_R_gc;
    params.P_N_rc = P_N_rc;
    params.P_N_gc = P_N_gc;
    params.inv_S = inv_S;
    params.inv_GS = inv_GS;
    params.inv_SWS = inv_SWS;
    params.inv_GSWGS = inv_GSWGS;
    params.G1 = G1;
    params.G2 = G2;
    params.G3 = G3;
    params.preprocess_G1 = preprocess_G1;
    params.preprocess_G2 = preprocess_G2;
    params.preprocess_G3 = preprocess_G3;
end