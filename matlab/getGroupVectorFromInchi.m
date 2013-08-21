function group_def = getGroupVectorFromInchi(inchi, silent)

if nargin < 2
    silent = true;
end

if isempty(inchi)
    group_def = [];
    return;
end

[~, ~, python_bin] = getBinaryPath();

if silent
    cmd = [python_bin ' ' getBasePath() 'python/inchi2gv.py -s -i '];
else
    cmd = [python_bin ' ' getBasePath() 'python/inchi2gv.py -i '];
end

if ~ispc
    [rval, group_def] = system([cmd '"' inchi '"']);
else
    [rval, group_def] = system([cmd inchi]);
end

if rval == 0 && ~strcmp('Traceback', group_def(1:9))
    group_def = regexp(group_def, ',', 'split');
    group_def = group_def(~cellfun('isempty',group_def));
    group_def = str2double(group_def);
else
    group_def = [];
end

