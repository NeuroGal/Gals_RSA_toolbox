function [edi, edi_stats, cat_edi, cat_edi_stats] = edi_calc(rdm, varargin)
% calculate exemplar discriminability index. See:
%       Nili, Hamed, et al. Plos one (2020)
%       Inferring exemplar discriminability in brain representations.
% 
% also enables calculating cdi (category discriminability).
% including permutation stats.
% note this assumes the rdm is symmetric and only uses tril (lower triangle).
% 
% input:
%   - rdm - n_stim x n_elec x n_time
% optional:
%   - 'stat_info' - followed by structure with the following fields:
%                       stat_type  - 'perm'\'cluster'. perm is corrected 
%                                     for MC with FDR correction. default: perm.
%                       n_perm     - default: 1000
%                       p_thresh   - default: 0.05 - cluster threshold \ fdr
%                       n_sides    - 1\2 (default: 1)
%                       n_clusters - default: Inf (all clusters)
%                       firstlevel_type - must be 'stat' (also default)
%                       firstlevel_thresh - no default (!) (if type =
%                           'stat' & 2 sided this could be 2-entry array (pos
%                           & neg thresh)) - assumes pos first
%   - 'categories' - and as usual followed by categories_vec then cat_names
%   - 'per_categ' - if you want to calculate edi per category
%   - 'cat_rdms'  - rdms calculated with subgroups according to the categories
%                   (done to use for that category only what was responsive
%                   to it) - cell of n_cat x 1
%   - 'get_perms' - output the permutations too
%
% output:
%   - edi - n_timex1 - between exemplars - within exemplars
%   - edi_stats - struct with p_values & mask ( [] if no stat_info given )
%   - cat_edi - struct with .cdi and .edi_per_categ (if 'per_categ') -
%               empty [] if no category info given.
%   - cat_edi_stats - struct with .cdi.p_values & .cdi.mask if stat_info &
%               category info is given (and similarly for edi_per_categ)

if ndims(rdm) == 4
    rdm = mean(rdm, 4);
end
if sum(diag(rdm(:,:,1))) == 0
    error('This seems to be a non-split-data-RDM, EDI is not defined in this case (CDI is, but needs to be implemented so it excludes the diagonal)') % technically we can still run cdi without the diag
end
n_stim = size(rdm, 1); n_time = size(rdm, 3);
categories_vec = ones(n_stim,1);
run_stats = false; per_categ = false;
use_cat_rdms = false; get_perms = false;

arg  = 1;
while arg <= size(varargin,2)
    switch varargin{arg}
        case 'categories'
            categories_vec = varargin{arg+1};
            if length(categories_vec) ~= n_stim
                error('Category indices vector should match the number of stimuli');
            end
            cat_names = varargin{arg+2};
            if ~iscell(cat_names)
                error('Category names should be given in a cell array')
            end
            if length(unique(categories_vec))~=length(cat_names)
                error("Number of category names doesn't match number of categories")
            end
            arg = arg + 3;
        case 'per_categ'
            per_categ = true;
            arg = arg + 1;
        case 'stat_info'
            stat_info = varargin{arg+1};
            stat_info = parse_stat_info(stat_info);
            run_stats = true;
            arg = arg + 2;
        case 'cat_rdms'
            cat_rdms = varargin{arg+1};
            use_cat_rdms = true;
            arg = arg + 2;
        case 'get_perms'
            get_perms = true;
            arg = arg + 1;
        otherwise
            error(['Unknown optional argument name: ' varargin{arg} '.']);
    end
end
n_cat = length(unique(categories_vec));

within_model = diag(true(n_stim,1));
edi = inner_edi(rdm, within_model);
if run_stats
    perm_edi = get_edi_perms(rdm, within_model, stat_info.n_perm);
    edi_stats = stat_wrapper(edi, stat_info, perm_edi);
    if get_perms; edi_stats.perms = perm_edi; end
else
    edi_stats = [];
end
cat_edi = []; cat_edi_stats = [];

if n_cat >= 2
    model_rdm = ~eye(n_cat,n_cat);
    model_rdm = model_rdm(categories_vec,categories_vec);
    cat_edi.cdi = inner_edi(rdm, ~model_rdm);
    if run_stats
        perm_cdi = get_edi_perms(rdm, ~model_rdm, stat_info.n_perm);
        cat_edi_stats.cdi = stat_wrapper(cat_edi.cdi, stat_info, perm_cdi);
        if get_perms; cat_edi_stats.cdi.perms = perm_cdi; end
    end
    if per_categ
        edi_per_categ = nan(n_time, n_cat);
        pvals_per_categ = nan(n_time, n_cat); masks_per_categ = nan(n_time, n_cat);
        perms_per_categ = nan(n_time, stat_info.n_perm, n_cat);
        if strcmp(stat_info.stat_type,'cluster'); cluster_pvals_per_categ = cell(1, n_cat); end
        for categ = 1:n_cat
            rel_idx = categories_vec==categ;
            if sum(rel_idx)==1; masks_per_categ(:, categ) = 0; continue; end
            within_model = diag(true(sum(rel_idx),1));
            if ~use_cat_rdms
                rdm_to_use = rdm(rel_idx, rel_idx, :);
            else
                rdm_to_use = cat_rdms{categ};
            end
            if isempty(rdm_to_use)
                continue
            else
                edi_per_categ(:, categ) = inner_edi(rdm_to_use, within_model);
                if run_stats
                    perms_per_categ(:, :, categ) = get_edi_perms(rdm_to_use, within_model, stat_info.n_perm);
                    cat_stats = stat_wrapper(edi_per_categ(:, categ), stat_info, perms_per_categ(:,:,categ));
                    pvals_per_categ(:, categ) = cat_stats.p_values;
                    masks_per_categ(:, categ) = cat_stats.mask;
                    if strcmp(stat_info.stat_type,'cluster'); cluster_pvals_per_categ{categ} = cat_stats.cluster_p; end
                end
            end
        end
        cat_edi.edi_per_categ = edi_per_categ;
        cat_edi.edi_per_categ_mean = nanmean(edi_per_categ,2);
        if run_stats
            cat_edi_stats.edi_per_categ.p_values = pvals_per_categ;
            cat_edi_stats.edi_per_categ.mask = masks_per_categ;
            if strcmp(stat_info.stat_type,'cluster')
                cat_edi_stats.edi_per_categ.cluster_p = cluster_pvals_per_categ;
            end
            mean_perms = nanmean(perms_per_categ,3);
            cat_edi_stats.edi_per_categ_mean = stat_wrapper(nanmean(edi_per_categ,2), stat_info, mean_perms);
            if get_perms
                cat_edi_stats.edi_per_categ.perms = perms_per_categ;
                cat_edi_stats.edi_per_categ_mean.perms = mean_perms;
            end
        end
    end
end

end

function edi = inner_edi(rdm, within_model)
% within_model = diag(true(n_stim,1)) if edi, and ~mode_rdm if cdi
n_stim = size(rdm, 1); n_time = size(rdm, 3);
lower_trig = tril(true(n_stim));
edi = nan(n_time,1);
for t = 1:n_time
    cur_rdm = rdm(:, :, t);
    within_idx = lower_trig & within_model;
    between_idx = lower_trig & ~within_model;
    edi(t) = mean(cur_rdm(between_idx)) - mean(cur_rdm(within_idx));
end
end

function perm_edi = get_edi_perms(rdm, within_model, n_perm)
n_stim = size(rdm, 1); n_time = size(rdm, 3);
perm_edi = nan(n_time, n_perm);

for p=1:n_perm
    perm_order = randperm(n_stim);
    rdm_perm = rdm(perm_order, :, :);
    perm_edi(:, p) = inner_edi(rdm_perm, within_model);
    if mod(p,100)==0
        fprintf('Done with permutation %d/%d\n',p,n_perm);
    end
end
end

function edi_stats = stat_wrapper(edi, stat_info, perm_edi)
edi_stats = [];
edi_stats.p_values = perm_pvals(edi, perm_edi, stat_info);
if strcmp(stat_info.stat_type, 'perm')
    [~, thresh] = FDR(edi_stats.p_values(:), stat_info.p_thresh);
    edi_stats.mask = edi_stats.p_values <= thresh;
elseif strcmp(stat_info.stat_type, 'cluster')
    [cluster_masks, edi_stats.cluster_p] = run_cluster_perm(stat_info, edi, reshape(perm_edi,size(perm_edi,1),1,size(perm_edi,2))); % can also get cluster p_vals from this  
    cluster_masks = cat(3, cluster_masks{:});
    edi_stats.mask = any(cluster_masks, 3); % so notice p values will be from the perm, because we are not calculating perm_p using cluster method
    if isempty(edi_stats.mask)
        edi_stats.mask = zeros(size(p_values));
    end
end
end

function stat_info = parse_stat_info(stat_info)
if ~isstruct(stat_info)
    error('stat_info should be a struct with the following fields: stat_type, p_thresh, n_sides, n_perm, n_clusters, firstlevel_type, firstlevel_thresh')
end
if ~isfield(stat_info,'stat_type')
    warning('Setting stat_type to perm'); stat_info.stat_type = 'perm';
end
if ~any(ismember(stat_info.stat_type,{'perm','cluster'}))
    warning('Unrecognized stat_type %s - setting to perm', stat_info.stat_type);
    stat_info.stat_type = 'perm';
end
if ~isfield(stat_info,'p_thresh')
    warning('Setting p_thresh to 0.05'); stat_info.p_thresh = 0.05;
end
if ~isfield(stat_info,'n_sides')
    warning('Setting n_sides to 1 (larger than)'); stat_info.n_sides = 1;
end
if ~isfield(stat_info,'n_perm')
    warning('Setting n_perm to 1000'); stat_info.n_perm = 1000;
end
if strcmp(stat_info.stat_type,'cluster')
    if ~isfield(stat_info,'n_clusters')
        warning('Setting n_clusters to Inf (all from largest to smallest)'); stat_info.n_clusters = Inf;
    end
    if ~isfield(stat_info,'firstlevel_type')
        warning('Setting firstlevel_type to stat'); stat_info.firstlevel_type = 'stat';
    end
    if strcmp(stat_info.firstlevel_type, 'p')
        error('firstlevel_type needs to be "stat", change and set firstlevel_thresh accordingly');
    end
    if ~isfield(stat_info,'firstlevel_thresh')
        error('You must set firstlevel_thresh');
    end
end
end