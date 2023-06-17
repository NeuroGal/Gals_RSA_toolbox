function diags = get_diags(array)
% assumes you have dim1 x dim1 x dim2 x ... and want dim1 x 1 x dim2 x ...
% especially useful for getting diagonals out of temporal generalization matrices
%
% This code was written as part of visualization of Vishne et al., biorxiv 2022 https://doi.org/10.1101/2022.08.02.502469
%   'Representing experience over time: sustained sensory patterns and transient frontoparietal patterns'
%   So please cite (-:
% Written by Gal Vishne, bug reports \ requests: gal.vishne@gmail.com

diags = cellfun_wrap(@diag,num2cell(array,[1 2]), true); 
end