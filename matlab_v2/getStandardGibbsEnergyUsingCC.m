%%%
%%%   Arguments:
%%%     reaction - a KeggReaction object
%%%   Returns:
%%%      the CC estimation for this reaction's untransformed dG0 (i.e.
%%%      using the major MS at pH 7 for each of the reactants)
function [dG0_cc, U] = getStandardGibbsEnergyUsingCC(params, reactions)
v_r = params.dG0_cc;
v_g = params.dG0_gc;
C1  = params.preprocess_C1;
C2  = params.preprocess_C2;
C3  = params.preprocess_C3;

X = zeros(shape(C1, 2), length(reactions));
G = zeros(shape(C2, 2), length(reactions));
for i = 1:length(reactions)
    x, g = % use python do decompose the reaction string into x and g: self._decompose_reaction(reaction)
    X(:, i) = x;
    G(:, i) = g;

dG0_cc = X' * v_r + G' * v_g;
U = X' * C1 * X + X' * C2 * G + G' * C2' * X + G' * C3 * G;
