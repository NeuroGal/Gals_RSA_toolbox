function [rdm, rdm_parts] = calculate_rdm(segs, distance_metric, split_data, timextime, do_zscore)
% Calculate representational dissimilarity matrix (RDM) of the segments in
% segs (per time-point) based on Matlab's pdist function.
% 
% Example usage: rdm = calculate_rdm(seg, 'correlation')
% Calculates RDM per time-point based on the correlation metric (1 -
% pearson correlation).
%
% Input: segs (dimensions: n_stim x n_features x n_time) (or ... x 2)
%           if 4D computes split-data RDMs (i.e. separately for each part),
%           this will be rdm_parts and rdm will be the averaged version.
%        distance_metric - method to calculate the dissimilarity, e.g.
%           'euclidean', 'correlation', 'cosine', 'mahalanobis' (see 
%           options for Matlab's pdist.
%   Optional: (default: False)
%        split_data - supress warning for doing sdRDM
%        timextime  - if true, output is n_stim x n_stim x n_time x n_time (x...2)
%         and idx s1,s2,t1,t2 is the dissimilarity of stim1 at tp1 to stim2 at tp2
%        do_zscore  - z-score the segments according to the first dimension
%
% Output: 
%   rdm dimensions: n_stim x n_stim x n_time (unless timextime = True)
%   rdm_parts dimensions: n_stim x n_stim x n_time x 2 if relevant.
%
% # TODO: 1) add an option for a 'decoding' distance metric (1-accuracy) with
% an option to use decision-value-weighting 2) add an option for cross validated
% distances(see Guggenmos, Sterzer, and Cichy. "Multivariate pattern analysis
% for MEG: A comparison of dissimilarity measures." NeuroImage (2018))
%
% Written by Gal Vishne, lab of Leon Y. Deouell, 2019-2021
% Send bug reports and requests to gal.vishne@gmail.com

if size(segs,2)<2
    error('Only 0\1 electrodes, are you sure?') % maybe make warning
end     
if ~exist('split_data','var') || isempty(split_data); split_data = false; end
if ~exist('timextime','var') || isempty(timextime); timextime = false; end
if ~exist('do_zscore','var') || isempty(do_zscore); do_zscore = false; end
if ndims(segs)==4 && ~split_data
    split_data = true;
    warning('Using split_data distances according to dim 4')
end
[n_stim, n_elec, n_time] = size(segs, 1:3);
if do_zscore; segs = zscore(segs); end % z-scores according to stim dim

if ~split_data; sec_idx = 1; else; sec_idx = 2; end
if ~timextime; tdim = 1; else; tdim = n_time; end

rdm_parts = nan(n_stim, n_stim, n_time, tdim, 2);
for t=1:n_time
    if ~timextime; sec_tidx = t; else; sec_tidx = 1:n_time; end
    X1 = segs(:, :, t, 1);
    X2 = reshape(permute(segs(:, :, sec_tidx, sec_idx), [2 1 3 4]), n_elec, [])';
    rdm_parts(:,:,t,:,1) = reshape(pdist2(X1,X2, distance_metric),n_stim,n_stim,length(sec_tidx));
    rdm_parts(:,:,t,:,2) = permute(reshape(pdist2(X2,X1, distance_metric),n_stim,length(sec_tidx),n_stim),[1 3 2]);       
end
rdm = mean(rdm_parts, 5);
rdm = squeeze(rdm); rdm_parts = squeeze(rdm_parts);

if ~split_data; rdm_parts = []; end
end