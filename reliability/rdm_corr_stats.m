function [mask, pvals, cluster_p, thresh, diag_stats] = rdm_corr_stats(stat_info, corr_r, corr_perm_r, corr_p, corr_perm_p)
% Compute stats for RDM x RDM correlations\reliability metrics (either data x data, or data x model).
%
% Includes correction for multiple comparisons using cluster permutations 
% (Maris & Oostenveld, 2007; https://doi.org/10.1016/j.jneumeth.2007.03.024),
% using max-statistic control (Nichols & Holmes, 2002; https://doi.org/10.1002/hbm.1058),
% or correction for False Discovery Rate according to Benjamini and Hochberg (1995)
% with code by Edden Gerber (https://github.com/edden-gerber/time_series_analysis_and_statistics/blob/master/testing/FDR.m)
%
% Written as part of the code for this paper (please cite):
%   Vishne et al., Cell Reports 2023 (biorxiv DOI, to be updated
%   when formally published): https://doi.org/10.1101/2022.08.02.502469
%   'Distinct Ventral Stream and Prefrontal Cortex Representational
%   Dynamics during Sustained Conscious Visual Perception'
%
% Input:
%   stat_info   - settings structure, see fields below.
%   corr_r      - correlation values (ntime1 x ntime2, both can be==1, e.g.
%                 in case of correlating with a category model).
%             *** this can also be z-scored correlations (reliability
%                 metrics I used in the paper)
%   corr_perm_r - permutation values (ntime1 x ntime2 x nperm). the precise
%                 computation depends on the type of reliability metric or
%                 correlation used.
%   corr_p, corr_perm_p - only relevant when the input is really a
%                 correlation matrix, this is the analytic p-values
%                 (point-by-point t-tests). sizes are the same as corr_r &
%                 corr_perm_r.
%
% stat_info fields:
%       get_diag_stats - Default: False.
%                      Whether to calculate separately the stats for the diagonal 
%                      (which won't be corrected for multiple comparisons for the full matrix)
%       stat_type    - Default: 'analytic' if corr_p inserted, otherwise 'perm'.
%                      Four options implemented:
%                        'analytic' - uses the input corr_p! (only keep when the
%                                     output was not z-scored) -> note I am
%                                     not checking this is legit, if corr_p
%                                     was inserted we assume it's fine.
%                        'perm'     - point by point comparison to the permutation
%                            -> both corrected for multiple comparisons using FDR correction
%                        'cluster'  - runs cluster permutations (needs more inputs, see below)
%                        'max'      - max statistic control for multiple comparisons
%       p_thresh     - Default: 0.05
%                      * This is for the cluster\max-stat threshold or fdr, not the cluster first level.
%       n_sides      - Default: 1 (right side)
%                      Options: 1 or 2. Decides if we check only higher from 
%                      chance (right sided test) or double sided.
%       chance_level - Default: 0, mandatory for 2-sided 'perm' or 'max'
%                      and the firtlevel_stat is assuemd to be around this.
%   Cluster inputs:
%       n_clusters        - Default: Inf (all clusters). Max # of clusters to return. 
%                           In the permutations only the largest cluster is taken, but we 
%                           can choose to compare more than just our largest cluster
%                           (it's just stricter, but valid, see the original paper for more details)
%       firstlevel_type   - Default: 'p' if possible (corr_p and corr_perm_p inserted). Otherwise: 'stat'
%                           How should the firstlevel_thresh be treated?
%                           Options: 'p'\'stat'\'zscore'
%                            'p' - p_value (only when the analytic is relevant!)
%                            'stat' - apply the threshold to the input directly
%                            'zscore' - zscore the values according to the
%                               permutations mean & std and apply the
%                               threshold according to that. This is only
%                               applied in the function run_cluster_perm
%                               (and there it also setting chance_level to 0)
%                             ! don't use this for zscored reliability values,
%                               these are already inserted as the z-scored statistic.
%       firstlevel_thresh - If p-value is possible, default: 0.05.
%                           Otherwise: MUST have user input.
%                           Note for firstlevel_type == 'stat':
%                           * If chance_level is given insert this as DISTANCE FROM CHANCE.
%                           firstlevel_type == 'stat'\'zscore':
%                           * If n_sides == 2 this could be 2-entry array
%                           [pos, neg] <-POSITIVE FIRST, neg really has to
%                           be negative we're taking it relative to chance_level.
%                           If only one given it's applied as [+thresh,-thresh]
%                           (in both cases it's around chance_level if that is given)       
%
% Output:
%   mask       - matrix of true\false for significant time-points
%                (size: ntime1 x ntime2). If you asked for
%                cluster-permutations these are grouped together.
%   pvals      - point-by-point p-values - uncorrected except for 'max', []
%                for cluster.
%   cluster_p  - only for 'cluster', otherwise []
%   thresh     - threshold used for FDR correction \ max-statistic
%                threshold (for 'cluster' this gives [])
%   diag_stats - only exported if stat_info.get_diag_stats == True.
%                fields: 
%                   vals (just the diagonal values of corr_r)
%                   mask
%                   pvals (see comment above regarding correction)
%                   cluster_p
%                   thresh
%
% Written by Gal Vishne, Deouell Lab ~2021
% Bug reports \ requests: gal.vishne@gmail.com

if ~exist('corr_p','var') || isempty(corr_p); corr_p = nan(size(corr_r)); end
if ~exist('corr_perm_p','var') || isempty(corr_perm_p); corr_perm_p = nan(size(corr_perm_r)); end
stat_info = stat_info_parser(stat_info, corr_p, corr_perm_p); % update defaults etc. 

cluster_p = []; diag_stats = []; % just output empty variables if not computed
if stat_info.get_diag_stats; diag_stats.vals = diag(corr_r); end
switch stat_info.stat_type
    case 'analytic'
        pvals = corr_p; % just return the corr_p that was inserted, it is supposed to come from the analytic results.
    case 'perm'
        pvals = perm_pvals(corr_r, corr_perm_r, stat_info);
    case 'cluster'
        pvals = []; thresh = [];
        [cluster_masks, cluster_p] = run_cluster_perm(stat_info, corr_r, corr_perm_r, corr_p, corr_perm_p);
        mask = finalize_clust_mask(cluster_masks, corr_p);
        if stat_info.get_diag_stats
            diag_stats.pvals = diag(corr_p);
            [cluster_masks, diag_stats.cluster_p] = cluster_perm(stat_info, diag_stats.vals, get_diags(corr_perm_r), diag_stats.pvals, get_diags(corr_perm_p));
            diag_stats.mask = finalize_clust_mask(cluster_masks, diag_stats.pvals);
            diag_stats.thresh = [];
        end
    case 'max'
        [mask, pvals, thresh] = max_stat_correction(corr_r, corr_perm_r, stat_info);
        cluster_p = [];
        if stat_info.get_diag_stats
            [diag_stats.mask, diag_stats.pvals, diag_stats.thresh] = max_stat_correction(diag_stats.vals, get_diags(corr_perm_r), stat_info);
            diag_stats.cluster_p = [];
        end
    otherwise
        error("No stat_type '%s'. Use 'analytic', 'perm', 'max' or 'cluster'", stat_info.stat_type)
end

% FDR correction & masks for the perm & analytic cases
if ismember(stat_info.stat_type, {'analytic','perm'})
    [~, thresh] = FDR(pvals(:), stat_info.p_thresh); % uses all values, change here if you want something else (e.g. take only above the diagonal since matrix is symmetric)    
    mask = pvals<=thresh; cluster_p = [];
    if stat_info.get_diag_stats
        diag_stats.pvals = diag(pvals);
        [~, diag_stats.thresh] = FDR(diag_stats.pvals, stat_info.p_thresh);
        diag_stats.mask = diag_stats.pvals<=diag_stats.thresh;
        diag_stats.cluster_p = []; 
    end
end
end

function mask = finalize_clust_mask(cluster_masks, corr_p)
% corr_p inserted just for getting the size right in case there are no clusters.
cluster_masks = cat(3, cluster_masks{:});
mask = any(cluster_masks, 3);
if isempty(mask); mask = false(size(corr_p)); end
end

function stat_info = stat_info_parser(stat_info, corr_p, corr_perm_p)
% Make sure all fields are there & correct. Add default values for missing things.
if ~(all(isnan(corr_p),'all')); pval_ok = false; else; pval_ok = true; end
fields = {'get_diag_stats','stat_type','p_thresh','n_sides','chance_level'};
defaults = {false,'analytic',0.05,1,0}; 
if ~pval_ok; defaults{2} = 'perm'; end
for f = 1:length(fields)
    if ~isfield(stat_info, fields{f}) 
        stat_info.(fields{f}) = defaults{f};
        warning('Setting %s to %s', fields{f}, string(defaults{f}))
    end
end
if ~pval_ok && strcmp(stat_info.stat_type, 'analytic')
    warning('No p-values given for analytic test, setting stat_type to perm'); stat_info.stat_type = 'perm';
end
if strcmp(stat_info.stat_type,'cluster')
    if ~(all(isnan(corr_perm_p),'all')); pval_ok = false; else; pval_ok = true; end
    if ~isfield(stat_info,'n_clusters')
        warning('Setting n_clusters to Inf (all ordered from largest to smallest)'); stat_info.n_clusters = Inf;
    end
    if ~isfield(stat_info,'firstlevel_type') || ~any(ismember(stat_info.firstlevel_type,{'stat','p','zscore'}))
        if pval_ok
            warning('Setting firstlevel_type to p-value'); stat_info.firstlevel_type = 'p';
        else
            warning('Setting firstlevel_type to stat'); stat_info.firstlevel_type = 'stat';
        end
    end
    if ~isfield(stat_info,'firstlevel_thresh')
        if pval_ok
            warning('Setting firstlevel_thresh to 0.05'); stat_info.firstlevel_thresh = 0.05;
        else
            error('You must set field: firstlevel_thresh');
        end
    end
end
end