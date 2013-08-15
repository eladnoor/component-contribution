clear all;

fprintf('***********************************\n');
fprintf('********** RUNNING TESTS **********\n');
fprintf('***********************************\n');
fprintf('If you do not see a "success" message\nin the end of the script,\nit means there is an error\n');

tempd = pwd();
cd('../matlab');
[babel_bin, cxcalc_bin] = getBinaryPath();

% Test KEGG related functions
fprintf('\n1) Testing getInchi() directly from KEGG online database\n')
[std_inchi, std_inchi_stereo, std_inchi_stereo_charge, nstd_inchi] = getInchi(2);
assert(strcmp(std_inchi, 'InChI=1S/C10H16N5O13P3/c11-8-5-9(13-2-12-8)15(3-14-5)10-7(17)6(16)4(26-10)1-25-30(21,22)28-31(23,24)27-29(18,19)20/h2-4,6-7,10,16-17H,1H2,(H,21,22)(H,23,24)(H2,11,12,13)(H2,18,19,20)'));
assert(strcmp(nstd_inchi, 'InChI=1/C10H16N5O13P3/c11-8-5-9(13-2-12-8)15(3-14-5)10-7(17)6(16)4(26-10)1-25-30(21,22)28-31(23,24)27-29(18,19)20/h2-4,6-7,10,16-17H,1H2,(H,21,22)(H,23,24)(H2,11,12,13)(H2,18,19,20)/t4-,6-,7-,10-/m1/s1/f/h18-19,21,23H,11H2'));

fprintf('\n2) Testing getInchis() for directly provided InChIs\n')
target_cids = [1, 9];
target_inchis = {'InChI=1/H2O/h1H2', 'InChI=1/H3O4P/c1-5(2,3)4/h(H3,1,2,3,4)/f/h1-3H'};
res_inchis = getInchies(target_cids, false, false);
assert(strcmp(res_inchis.std_inchi{1}, 'InChI=1S/H2O/h1H2'));
assert(strcmp(res_inchis.std_inchi{2}, 'InChI=1S/H3O4P/c1-5(2,3)4/h(H3,1,2,3,4)'));
assert(strcmp(res_inchis.nstd_inchi{1}, 'InChI=1/H2O/h1H2'));
assert(strcmp(res_inchis.nstd_inchi{2}, 'InChI=1/H3O4P/c1-5(2,3)4/h(H3,1,2,3,4)/f/h1-3H'));

fprintf('\n3) Testing cxcalc using command line\n')
[success, ~] = system(cxcalc_bin);
assert(success == 0);
[success, cxcalc_stdout] = system([cxcalc_bin ' "InChI=1/H2O/h1H2" pka -a 1 -b 0']);
assert(success == 0);
assert(strcmp(cxcalc_stdout, sprintf('id\tapKa1\tatoms\n1\t4.27\t3\n')));

fprintf('\n4) Testing getKeggPkas() - using ChemAxon script\n')
res_pkas = getKeggpKas(target_cids, target_inchis, 20);
assert(isempty(res_pkas(1).pKas));
assert(all(diag(res_pkas(2).pKas, 1) == [12.9; 6.95; 1.8]));
assert(all(res_pkas(2).majorMSpH7 == [0; 1; 0; 0]));
assert(all(res_pkas(2).zs == [-3 -2 -1 0]));
assert(all(res_pkas(2).nHs == [0 1 2 3]));

fprintf('\n5) Testing getGroupVectorFromInchi() - the group decomposition Python script\n')
res_group_Pi = getGroupVectorFromInchi(target_inchis{2}, false);
assert(isempty(res_group_Pi)); % Pi is not decomposable
inchi_ATP = 'InChI=1/C10H16N5O13P3/c11-8-5-9(13-2-12-8)15(3-14-5)10-7(17)6(16)4(26-10)1-25-30(21,22)28-31(23,24)27-29(18,19)20/h2-4,6-7,10,16-17H,1H2,(H,21,22)(H,23,24)(H2,11,12,13)(H2,18,19,20)/t4-,6-,7-,10-/m1/s1/f/h18-19,21,23H,11H2';
res_group_ATP = getGroupVectorFromInchi(inchi_ATP, false);
assert(all(find(res_group_ATP) == [44 51 73 75 104 120 128 148 149 150 158 161 163]));
assert(all(res_group_ATP(res_group_ATP ~= 0) == [1 1 1 3 1 2 1 2 2 1 2 1 1]));

cd(tempd);

msgbox('Success! All unit tests have completed without errors');
