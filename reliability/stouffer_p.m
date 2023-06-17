function p_combined = stouffer_p(p_values, comb_dims, weights)
% Stouffer et al., 1951 (taken from Helfrich et al., 2018)
% p_values can be an n_d array, assuming that the dim over which to
% combine is the last one, unless given other inputs.
% see also: https://en.wikipedia.org/wiki/Fisher%27s_method
%
% Written by Gal Vishne, lab of Leon Y. Deouell, 2020
% Send bug reports and requests to gal.vishne@gmail.com

size_original = size(p_values);n_dims = length(size_original);
if ~exist('comb_dims','var') || isempty(comb_dims)
    comb_dims = n_dims;
    if size_original(n_dims) == 1 % correct for the case of a column vector
        comb_dims = 1;
    end
end
if ~exist('weights','var')
    weights = ones(size(p_values));
end
if comb_dims > n_dims
    error('comb_dim needs to be <= n_dims')
end
p_cell = num2cell(p_values, comb_dims);
weights_cell = num2cell(weights, comb_dims);

p_combined = cellfun(@stouffer_1d, p_cell, weights_cell);
end

function p_1d = stouffer_1d(p_values, weights)
keep_idx = ~isnan(p_values);
weights = weights(keep_idx);p_values = p_values(keep_idx);
single_z = -sqrt(2)*erfcinv(2*(1-p_values));%norminv(1-p_values);
x = sum(single_z.*weights,'all')/norm(weights(:));
p_1d = 1-erfc(-x/sqrt(2))/2;%1-normcdf(x);
end