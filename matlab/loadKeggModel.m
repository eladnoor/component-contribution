% [S, cids] = loadKeggModel(fname, fmt)
%
% reads and parses a file containing only reactions in KEGG format
%
% Arguments:
%   fname - a string containing the name of the input file
%   arrow - the string used as the 'arrow' in each reaction (default: '<=>')
%   has_reaction_ids - a boolean flag indicating if there is a column of
%                      reaction IDs (separated from the reaction with
%                      whitespaces)
%
% Return values:
%   S     - a stoichiometric matrix
%   cids  - the KEGG compound IDs in the same order as the rows of S

function [S, cids] = loadKeggModel(fname, arrow, has_reaction_ids)

fid = fopen(fname);
allData = textscan(fid, '%s', 'Delimiter', '\n');
reactionStrings = allData{1,1};
fclose(fid);

if nargin < 2
    [S, cids] = parseKeggModel(reactionStrings);
elseif nargin < 3
    [S, cids] = parseKeggModel(reactionStrings, arrow);
else
    [S, cids] = parseKeggModel(reactionStrings, arrow, has_reaction_ids);
end