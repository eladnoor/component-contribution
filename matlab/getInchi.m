function [std_inchi, std_inchi_stereo, std_inchi_stereo_charge, nstd_inchi] = getInchi(cid, inchi)

babel_cmd = 'babel';

[success, ~] = system([babel_cmd, ' -H']);
if success ~= 0
    error('Please make sure the command line program "babel" is installed and in PATH');
end

std_inchi = '';
std_inchi_stereo = '';
std_inchi_stereo_charge = '';
nstd_inchi = '';

if nargin < 2
    % Get the MOL for this compounds from KEGG
    [mol, status] = urlread(sprintf('http://rest.kegg.jp/get/cpd:C%05d/mol', cid));
    if status == 0
        warning('cannot read the URL for KEGG compound C%05d', cid);
        return;
    end    
    if isempty(mol)
        warning('the MOL descriptor for KEGG compound C%05d is empty', cid);
        return
    end
    
    if ispc
        cmd = ['echo ' mol ' | ' babel_cmd ' -imol -oinchi ---errorlevel 0 -w'];
    else
        mol = regexprep(mol,'^"|"$',''); % replaces " but only if it is at the beginning/end
        cmd = ['echo "' mol '" | ' babel_cmd ' -imol -oinchi ---errorlevel 0 -w'];
    end
else
    if isempty(inchi) % happens when a compound was added to KEGG but doesn't have an inchi
        return
    end
    if ispc
        cmd = ['echo ' inchi ' | ' babel_cmd ' -iinchi -oinchi ---errorlevel 0 -w'];
    else
        inchi = regexprep(inchi,'^"|"$',''); % replaces " but only if at the beginning/end
        cmd = ['echo "' inchi '" | ' babel_cmd ' -iinchi -oinchi ---errorlevel 0 -w'];
    end
    
end

[~, std_inchi] = system([cmd ' -xT/noiso/nochg/nostereo']);

if not(isempty(strfind(std_inchi,'command not found'))) || not(isempty(strfind(std_inchi,'error')))
        warning('The Inchi of compound C%05d has conversion problems', cid);
        return
end

if ~isempty(std_inchi) && strcmp('InChI=',std_inchi(1:6))
    std_inchi = strtok(std_inchi);
else
    std_inchi = '';
end

[~, std_inchi_stereo] = system([cmd ' -xT/noiso/nochg']);
if ~isempty(std_inchi_stereo) && strcmp('InChI=',std_inchi_stereo(1:6))
    std_inchi_stereo = strtok(std_inchi_stereo);
else
    std_inchi_stereo = '';
end

[~, std_inchi_stereo_charge] = system([cmd ' -xT/noiso']);
if ~isempty(std_inchi_stereo_charge) && strcmp('InChI=',std_inchi_stereo_charge(1:6))
    std_inchi_stereo_charge = strtok(std_inchi_stereo_charge);
else
    std_inchi_stereo_charge = '';
end

[~, nstd_inchi] = system([cmd ' -xFT/noiso']);
if ~isempty(nstd_inchi) && strcmp('InChI=',nstd_inchi(1:6))
    nstd_inchi = strtok(nstd_inchi);
else
    nstd_inchi = '';
end

