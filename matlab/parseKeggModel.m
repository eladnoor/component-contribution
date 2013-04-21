% [S, cids] = parseKeggModel(fname, fmt)
%
% parses cellarray of reactions in KEGG format
%
% Arguments:
%   reactionStrings - a 1D cellarray containing strings of reactions
%   arrow - the string used as the 'arrow' in each reaction (default: '<=>')
%   has_reaction_ids - a boolean flag indicating if there is a column of
%                      reaction IDs (separated from the reaction with
%                      whitespaces)
%
% Return values:
%   S     - a stoichiometric matrix
%   cids  - the KEGG compound IDs in the same order as the rows of S

function [S, cids] = parseKeggModel(reactionStrings, arrow, has_reaction_ids)

if nargin < 2
    arrow = '<=>';
end
if nargin < 3
    has_reaction_ids = false;
end

cids = [];
reactions = {};
for i = 1:length(reactionStrings)
    curr_line = reactionStrings{i};
    if has_reaction_ids
        tokens = regexp(curr_line, '(\w+)\s+(.*)', 'tokens');
        curr_line = tokens{1}{2};
    end
    sprs = reaction2sparse(curr_line, arrow);
    cids = unique([cids, find(sprs)]);
    reactions = [reactions, {sprs}];
end

S = zeros(length(cids), length(reactions));
for i = 1:length(reactions)
    r = full(reactions{i});
    if norm(r) ~= 0
        S(ismember(cids, find(r)), i) = r(r ~= 0);
    end
end
fprintf('Loaded a KEGG model with %d compounds and %d reactions\n', ...
        size(S, 1), size(S, 2));
