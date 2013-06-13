% I = subset_index(A,B)
%
% finds all the occurances of members of B in the larger set A, and returns
% the indices of their location in A, conserving the order of values in B.
% Assumes A and B contain only unique values, but not necessarily ordered.
% If B is not a subset of A, the extra members will be ignored.
%
% Examples:
% >> I = subset_index([6,8,3,4,5], [1,5,4])
%
% I =
%
%      5     4
%
%
function I = subset_index(A,B)

[Lia,Locb] = ismember(A,B);
if sum(Lia) ~= length(B)
    error('B is not unique or is not a proper subset of A');
end
I1 = find(Lia);
[~, I2] = sort(Locb(I1));
I = I1(I2);