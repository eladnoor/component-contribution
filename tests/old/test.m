clear all;

fprintf('***********************************\n');
fprintf('********** RUNNING TESTS **********\n');
fprintf('***********************************\n');
fprintf('If you do not see a "success" message\nin the end of the script,\nit means there is an error\n');

tempd = pwd();
cd('../matlab');
[babel_bin, cxcalc_bin, python_bin] = getBinaryPath();
target_cids = [1, 2, 9];
target_smiles = {'[OH2]', 'C([C@@H]1[C@H]([C@H]([C@H](n2cnc3c(N)ncnc23)O1)O)O)OP(=O)(O)OP(=O)(O)OP(=O)(O)O', 'OP(O)(O)=O'};
target_std_inchis = {'InChI=1S/H2O/h1H2', ...
                     'InChI=1S/C10H16N5O13P3/c11-8-5-9(13-2-12-8)15(3-14-5)10-7(17)6(16)4(26-10)1-25-30(21,22)28-31(23,24)27-29(18,19)20/h2-4,6-7,10,16-17H,1H2,(H,21,22)(H,23,24)(H2,11,12,13)(H2,18,19,20)', ...
                     'InChI=1S/H3O4P/c1-5(2,3)4/h(H3,1,2,3,4)'};
target_nstd_inchis = {'InChI=1/H2O/h1H2', ...
                      'InChI=1/C10H16N5O13P3/c11-8-5-9(13-2-12-8)15(3-14-5)10-7(17)6(16)4(26-10)1-25-30(21,22)28-31(23,24)27-29(18,19)20/h2-4,6-7,10,16-17H,1H2,(H,21,22)(H,23,24)(H2,11,12,13)(H2,18,19,20)/t4-,6-,7-,10-/m1/s1/f/h18-19,21,23H,11H2', ...
                      'InChI=1/H3O4P/c1-5(2,3)4/h(H3,1,2,3,4)/f/h1-3H'};
target_formulas = {'H2O', 'C10H16N5O13P3', 'H3O4P'};
target_nHs = [2, 16, 3];
target_zs = [0, 0, 0];

target_pkas.pKas = {[], [12.6; 7.42; 5.0; 3.29; 2.53; 0.9], [12.9; 6.95; 1.8]};
target_pkas.majorMSpH7 = {1, [0; 0; 1; 0; 0; 0; 0], [0; 1; 0; 0]};
target_pkas.zs = {0, [-5 -4 -3 -2 -1 0 1], [-3 -2 -1 0]};
target_pkas.nHs = {2, [11 12 13 14 15 16 17], [0 1 2 3]};

target_groupvecs = {[]; ...
                    sparse(ones(1,13), [44 51 73 75 104 120 128 148 149 150 158 161 163], [1 1 1 3 1 2 1 2 2 1 2 1 1]);
                    []};

test_count = 0;
                
%%
test_count = test_count + 1;
fprintf('\n%d) Testing all binaries\n', test_count);
[success, ~] = system([babel_bin ' --help']);
assert(success == 0);
[success, ~] = system([cxcalc_bin ' --help']);
assert(success == 0);
[success, ~] = system([python_bin ' --help']);
assert(success == 0);

%%
test_count = test_count + 1;
fprintf('\n%d) Testing babel through command line\n', test_count);
for i = 1:length(target_cids)
    if ispc
        [success, babel_stdout] = system(['echo ' target_smiles{i} '|' babel_bin ' -ismiles -oinchi ---errorlevel 0 -xFT/noiso']);
    else
        [success, babel_stdout] = system(['echo "' target_smiles{i} '" | ' babel_bin ' -ismiles -oinchi ---errorlevel 0 -xFT/noiso']);
    end
    assert(success == 0);
    babel_lines = regexp(babel_stdout, '\n', 'split');
    assert(strcmp(babel_lines{1}, target_nstd_inchis{i}));
end

%%
% Test KEGG related functions
test_count = test_count + 1;
fprintf('\n%d) Testing getInchi() directly from KEGG online database\n', test_count);
for i = 1:length(target_cids)
    [std_inchi, std_inchi_stereo, std_inchi_stereo_charge, nstd_inchi] = getInchi(target_cids(i));
    assert(strcmp(std_inchi, target_std_inchis{i}));
    assert(strcmp(nstd_inchi, target_nstd_inchis{i}));
end

%%
test_count = test_count + 1;
fprintf('\n%d) Testing getInchis() for directly provided InChIs\n', test_count);
res_inchis = getInchies(target_cids, false, false);
assert(isempty(setdiff(res_inchis.std_inchi, target_std_inchis)));
assert(isempty(setdiff(res_inchis.nstd_inchi, target_nstd_inchis)));

%%
test_count = test_count + 1;
fprintf('\n%d) Testing cxcalc using command line\n', test_count);
[success, ~] = system(cxcalc_bin);
assert(success == 0);
if ispc
    [success, cxcalc_stdout] = system(['echo ' target_nstd_inchis{1} '|' cxcalc_bin ' pka -a 1 -b 0']);
else
    [success, cxcalc_stdout] = system(['echo "' target_nstd_inchis{1} '" | ' cxcalc_bin ' pka -a 1 -b 0']);
end
assert(success == 0);
assert(strcmp(cxcalc_stdout, sprintf('id\tapKa1\tatoms\n1\t15.70\t1\n')));

%%
test_count = test_count + 1;
fprintf('\n%d) Testing getFormulaAndChargeFromInChI() - using ChemAxon script\n', test_count);
for i = 1:length(target_cids)
    [formula, nH, charge] = getFormulaAndChargeFromInChI(target_nstd_inchis{i});
    assert(strcmp(formula, target_formulas{i}));
    assert(nH == target_nHs(i));
    assert(charge == target_zs(i));
end

%%
test_count = test_count + 1;
fprintf('\n%d) Testing getKeggPkas() - using ChemAxon script\n', test_count)
res_pkas = getKeggpKas(target_cids, target_nstd_inchis, 20);
for i = 1:length(target_cids)
    assert(all(diag(res_pkas(i).pKas, 1) == target_pkas.pKas{i}));
    assert(all(res_pkas(i).majorMSpH7 == target_pkas.majorMSpH7{i}));
    assert(all(res_pkas(i).zs == target_pkas.zs{i}));
    assert(all(res_pkas(i).nHs == target_pkas.nHs{i}));
end

%%
test_count = test_count + 1;
fprintf('\n%d) Testing getGroupVectorFromInchi() - the group decomposition Python script\n', test_count)

for i = 1:length(target_cids)
    groupvec = getGroupVectorFromInchi(target_nstd_inchis{i}, false);
    assert(all(groupvec == target_groupvecs{i}));
end

%%
cd(tempd);
msgbox('Success! All unit tests have completed without errors');
