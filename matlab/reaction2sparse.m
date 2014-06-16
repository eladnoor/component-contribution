function res = reaction2sparse(s, arrow)
if nargin < 2
    arrow = '=';
end
tmp = regexp(s, sprintf('\\s*%s\\s*', arrow), 'split');
left = regexp(tmp{1}, '\s*\+\s*', 'split');
right = regexp(tmp{2}, '\s*\+\s*', 'split');

res = {};

d = regexp(left, '\s*(\d+ )?(C\d+)\s*', 'tokens');
for i = 1:length(d)
    d{i}{1} = d{i}{1}(~cellfun('isempty',d{i}{1}));
    if length(d{i}{1}) == 1
        cid = d{i}{1}{1};
        coeff = -1;
    else
        coeff = -str2double(d{i}{1}{1});
        cid = d{i}{1}{2};
    end
    res{end+1, 1} = cid;
    res{end, 2} = coeff;
end

d = regexp(right, '\s*(\d+ )?(C\d+)\s*', 'tokens');
for i = 1:length(d)
    d{i}{1} = d{i}{1}(~cellfun('isempty',d{i}{1}));
    if length(d{i}{1}) == 1
        cid = d{i}{1}{1};
        coeff = 1;
    else
        coeff = str2double(d{i}{1}{1});
        cid = d{i}{1}{2};
    end
    res{end+1, 1} = cid;
    res{end, 2} = coeff;
end
