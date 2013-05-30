clear all;

fprintf('********** RUNNING TESTS **********\n');
fprintf('If you do not see a "success" message in the end of the script, \n');

tempd = pwd();
cd('../matlab');

% Test KEGG related functions
target_cids = [1, 9];
target_inchis = {'InChI=1/H2O/h1H2', 'InChI=1/H3O4P/c1-5(2,3)4/h(H3,1,2,3,4)/f/h1-3H'};

res_inchis = getInchies(target_cids, false, false);
assert(strcmp(res_inchis.std_inchi{1}, 'InChI=1S/H2O/h1H2'));
assert(strcmp(res_inchis.std_inchi{2}, 'InChI=1S/H3O4P/c1-5(2,3)4/h(H3,1,2,3,4)'));
assert(strcmp(res_inchis.nstd_inchi{1}, 'InChI=1/H2O/h1H2'));
assert(strcmp(res_inchis.nstd_inchi{2}, 'InChI=1/H3O4P/c1-5(2,3)4/h(H3,1,2,3,4)/f/h1-3H'));

res_pkas = getKeggpKas(target_cids, target_inchis, 20);
assert(isempty(res_pkas(1).pKas));
assert(all(diag(res_pkas(2).pKas, 1) == [12.9; 6.95; 1.8]));
assert(all(res_pkas(2).majorMSpH7 == [0; 1; 0; 0]));
assert(all(res_pkas(2).zs == [-3 -2 -1 0]));
assert(all(res_pkas(2).nHs == [0 1 2 3]));

% Test the group decomposition Python script
res_group_Pi = getGroupVectorFromInchi(target_inchis{2}, false);
assert(isempty(res_group_Pi)); % Pi is not decomposable

inchi_ATP = 'InChI=1/C10H16N5O13P3/c11-8-5-9(13-2-12-8)15(3-14-5)10-7(17)6(16)4(26-10)1-25-30(21,22)28-31(23,24)27-29(18,19)20/h2-4,6-7,10,16-17H,1H2,(H,21,22)(H,23,24)(H2,11,12,13)(H2,18,19,20)/t4-,6-,7-,10-/m1/s1/f/h18-19,21,23H,11H2';
res_group_ATP = getGroupVectorFromInchi(inchi_ATP, false);
assert(all(find(res_group_ATP) == [44 51 73 75 104 120 128 148 149 150 158 161 163]));
assert(all(res_group_ATP(res_group_ATP ~= 0) == [1 1 1 3 1 2 1 2 2 1 2 1 1]));

cd(tempd);

msgbox('Success! All unit tests have completed without errors');
