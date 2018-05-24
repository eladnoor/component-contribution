pH = 7;
I = 0.2;
T = 300;
R = 8.31e-3; % kJ/mol/K

cd ../matlab;

inds = [157, 183];

%KeggSpeciespKa = getKeggpKas(training_data.cids(inds), training_data.nstd_inchi(inds));

for i = 1:length(inds)
    cid = training_data.cids(inds(i));

    % find the diss table from the training data structure
    k = find(cell2mat({training_data.kegg_pKa.cid}) == cid);
    if isempty(k)
        error('cannot find C%05d in the training data', cid);
    end
    diss = training_data.kegg_pKa(k);

    dG0s = cumsum(-[0, diag(diss.pKas, 1)'] * R * T * log(10));
    dG0s = dG0s - dG0s(diss.majorMSpH7);
    pseudoisomers = [dG0s(:), diss.nHs(:), diss.zs(:)]
    alpha = (9.20483*T)/10^3 - (1.284668*T^2)/10^5 + (4.95199*T^3)/10^8; % Approximation of the temperature dependency of ionic strength effects
    DH = (alpha * sqrt(I)) / (1 + 1.6 * sqrt(I)); % Debye Huckel

    % dG0' = dG0 + nH * (RTlog(10) * pH + DH) + charge^2 * DH;
    dG0_prime_vector = pseudoisomers(:, 1) + ...
                       pseudoisomers(:, 2) * (R*T*log(10)*pH + DH) - ...
                       pseudoisomers(:, 3).^2 * DH;

    dG0_prime = -R * T * maxstar(dG0_prime_vector / (-R * T));
    fprintf('C%05d: ddG0_prime = %.2f\n', cid, dG0_prime);
end
cd ../tests;