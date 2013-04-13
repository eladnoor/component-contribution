function p = cdfplot(X, varargin)

n = length(X);
p = plot(sort(X), (0:(n-1))/(n-1), varargin{:});
