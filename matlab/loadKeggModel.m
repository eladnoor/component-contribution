function [S, cids] = loadKeggModel(fname, fmt)

fid = fopen(fname);
allData = textscan(fid, '%s', 'Delimiter', '\n');
reactionStrings = allData{1,1};
fclose(fid);

if nargin < 2
    [S, cids] = parseKeggModel(reactionStrings);
else
    [S, cids] = parseKeggModel(reactionStrings, fmt);
end