function rdm_reshaped = reshape_dists(rdm)
% turns n_stim x n_stim x n_time to n_dists x n_time
% takes from each time-point only the upper triangle!
%
% Written as part of the code for this paper:
%   Vishne et al., Cell Reports 2023 (biorxiv DOI, to be updated
%   when formally published): https://doi.org/10.1101/2022.08.02.502469
%   'Distinct Ventral Stream and Prefrontal Cortex Representational
%   Dynamics during Sustained Conscious Visual Perception'
% Bug reports \ requests: gal.vishne@gmail.com

n_stim = size(rdm,1); n_time = size(rdm,3);
rdm_reshaped = nan(n_stim*(n_stim-1)/2, n_time);
for t = 1:n_time
    rdmtmp_t = rdm(:,:,t);
    rdm_reshaped(:,t) = rdmtmp_t(triu(true(size(rdmtmp_t)), 1));
end
end
