function cell_array = cellfun_wrap(func, array,use_cell2mat)
% This code was written as part of visualization of Vishne et al., biorxiv 2022 https://doi.org/10.1101/2022.08.02.502469
%   'Representing experience over time: sustained sensory patterns and transient frontoparietal patterns'
%   So please cite (-:
% Written by Gal Vishne, bug reports \ requests: gal.vishne@gmail.com
cell_array = cellfun(func, array, 'UniformOutput', false);
if exist('use_cell2mat','var') && use_cell2mat; cell_array = cell2mat(cell_array); end
end