function [S, cids] = parseKeggModel(reactionStrings, fmt)

if nargin < 2 || strcmp(fmt, 'bastian')
    arrow = '<=>';
    has_reaction_ids = false;
elseif strcmp(fmt, 'rienk')
    arrow = '-->';
    has_reaction_ids = true;
else
    error('Kegg Model format must be "bastian" or "rienk"');
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
