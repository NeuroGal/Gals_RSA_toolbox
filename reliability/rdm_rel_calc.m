function [corr_r, mask, pvals, cluster_p, corr_perm_r, diag_stats] = rdm_rel_calc(rdm, reliability_type, stat_info, varargin)
% Compute single-exemplar reliability metrics according to the paper:
%   Vishne et al., Cell Reports 2023 (biorxiv DOI, to be updated
%   when formally published): https://doi.org/10.1101/2022.08.02.502469
%   'Distinct Ventral Stream and Prefrontal Cortex Representational
%   Dynamics during Sustained Conscious Visual Perception'
% Please cite if used.
%
% This function computes both the standard metrics & the same metrics after
% removing category structure\any other model rdm (see inputs below).
%
% Inputs:
%   Mandatory:
%     rdm - 2*n_stim x 2*n_stim x n_time, as shown in Figure 5b of the paper
%     reliability_type - 'item'\'geom'\'geom-v2' - the first two are used
%       in the paper, the last is similar to 'geom' but averages all pairs of
%       rdms from the large RDM inserted into the function.
%     stat_info - see rdm_corr_stats.m for more details, here we need:
%       stat_type, firstlevel_type - used to decide if we extract
%           permutation data including analytic p-values for each
%           permutation.
%       n_perm - how many permutations to run
%       n_sides - 1/2-sided test (1-sided is only right-sided)
%   Optional:
%     'corr_type' - followed by string indicating which correlation measure to use as the basis 
%       for the reliability metrics (anything matlab's corr\partialcorr recieve), default: Spearman. 
%     'model_rdm' - followed by a model RDM size n_stim x n_stim (!) which
%       is partialled out when we compute the reliability measure. Default
%       is to not partial anything.
%     'z_by_perm' - whether to z-score the resulting correlations by the
%       permutation values (this is how it was computed in the paper).
%    	Default: False.
%
% Output: 
%   corr_r - the reliability metric - n_time x n_time
%   for the others see rdm_corr_stats.m
%
% Written by Gal Vishne, lab of Leon Y. Deouell, ~2022
% Send bug reports and requests to gal.vishne@gmail.com

if ~ismember(reliability_type,["item","geom","geom-v2"]); error('Unknown type: %s.', reliability_type); end
corr_type = 'Spearman'; model_rdm = []; do_z_by_perm = false;
arg  = 1;
while arg <= size(varargin,2)
    switch varargin{arg}
        case 'corr_type'
            corr_type = varargin{arg+1};
            arg = arg + 1;
        case 'model_rdm'
            model_rdm = varargin{arg+1};
            arg = arg + 1;
        case 'z_by_perm'
            do_z_by_perm = true;
        otherwise
            error(['Unknown optional argument name: ' varargin{arg} '.']);
    end
    arg = arg + 1;
end
[stat_info, tail_type] = parse_stat_info(stat_info); % only checking the fields we need here, others are in rdm_corr_stats

corr_inputs = {'Type',corr_type,'Tail',tail_type};
get_corr_inputs = {rdm, reliability_type, corr_inputs, model_rdm};
need_p = true; % (saves time not to calculate it so we skip that if we can, here we do, but it's not needed for the permutations sometimes)
perm = false;
[corr_r, corr_p] = get_corr_mat(get_corr_inputs{:}, need_p, perm); % correlation matrix for the actual data

stat_type = string(stat_info.stat_type);
corr_perm_r = []; corr_perm_p = []; 
if stat_type~="analytic" || do_z_by_perm
    if ~(stat_type=="cluster" && strcmp(stat_info.firstlevel_type,'p')) % we update p only in this case we need it
        need_p = false;
    end
    [corr_perm_r, corr_perm_p] = get_corr_perms(get_corr_inputs{:}, need_p, stat_info.n_perm);
end
if do_z_by_perm; [corr_r, corr_perm_r] = z_by_perm(corr_r, corr_perm_r); end
[mask, pvals, cluster_p, ~, diag_stats] = rdm_corr_stats(stat_info, corr_r, corr_perm_r, corr_p, corr_perm_p);
end

function [corr_r, corr_p] = get_corr_mat(rdm, reliability_type, corr_inputs, model_rdm, need_p, perm)
n_stim = size(rdm,1)/2; perm_n = n_stim - boolean(reliability_type=="item"); use_model = ~isempty(model_rdm);
if perm; perm_order = randperm(perm_n); else; perm_order = 1:perm_n; end
if reliability_type == "item"
    dims = [size(rdm,3), size(rdm,3), n_stim*2]; corr_r = nan(dims); corr_p = nan(dims);
    for stim = 1:(n_stim*2)
        stim_short = stim; add_offset = [0, n_stim];
        if stim_short > n_stim; stim_short = stim_short-n_stim; add_offset = [n_stim 0]; end
        keep_idx = [1:(stim_short-1),(stim_short+1):n_stim]; % not correlating with the same stimulus
        rdm_1_row = squeeze(rdm(stim, add_offset(1)+keep_idx,:)); rdm_2_row = squeeze(rdm(stim,add_offset(2)+keep_idx(perm_order),:)); % not actually rows but stim x time
        if use_model; model_row = model_rdm(stim_short,keep_idx)'; else; model_row = []; end
        if ~need_p % see explanation in the main function
        corr_r(:,:,stim) = inner_corr(use_model, corr_inputs, rdm_1_row, rdm_2_row, model_row, need_p);
        else
        [corr_r(:,:,stim), corr_p(:,:,stim)] = inner_corr(use_model, corr_inputs, rdm_1_row, rdm_2_row, model_row, need_p);
        end
    end
    corr_r = mean(corr_r, 3);
    % stouffer method is a way of averaging p-values
    % doing mean here would be wrong, but stouffer is a bit slow so if p is not needed mean just fixes the sizes and it's faster
    if need_p; corr_p = stouffer_p(corr_p, 3); else; corr_p = mean(corr_p,3); end
elseif contains(reliability_type,"geom")
    % which parts of the rdm is averaged before we do the correlation.
    split_pairs = {[1 1;1 2],[2 1;2 2]}; % *** these are the pairs I used in the paper, but it's not really crucial (see geom-v2)
    split_num = [ones(n_stim,1);2*ones(n_stim,1)];
    rdm(n_stim + (1:n_stim),:,:) = rdm(n_stim+perm_order,:,:);
    rdm(:,n_stim + (1:n_stim),:) = rdm(:,n_stim+perm_order,:);
    if use_model; model_rdm = reshape_dists(model_rdm); end
    reshape_short = @(splits) reshape_dists(rdm(split_num==splits(1), split_num==splits(2),:)); % splits [s1, s2], calling the function reshape_dists
    if reliability_type == "geom"
        reshape_merge = @(splits) mean(cat(3,reshape_short(splits(1,:)),reshape_short(splits(2,:))),3);  % splits [s1, s2; s3, s4]
        rdm_1 = reshape_merge(split_pairs{1}); rdm_2 = reshape_merge(split_pairs{2});
        [corr_r, corr_p] = inner_corr(use_model, corr_inputs, rdm_1, rdm_2, model_rdm, need_p);
    elseif reliability_type == "geom-v2"
        all_split_pairs = nchoosek({[1 1],[1 2],[2 1],[2 2]},2);
        dims = [size(rdm,3), size(rdm,3), size(all_split_pairs,1)]; corr_r = nan(dims); corr_p = nan(dims);
        for rdm_pair = 1:size(all_split_pairs,1)
            rdm_1 = reshape_short(all_split_pairs{rdm_pair,1}); rdm_2 = reshape_short(all_split_pairs{rdm_pair,2});
            if ~need_p
            [corr_r(:,:,rdm_pair), corr_p(:,:,rdm_pair)] = inner_corr(use_model, corr_inputs, rdm_1, rdm_2, model_rdm, need_p);
            else
            corr_r(:,:,rdm_pair) = inner_corr(use_model, corr_inputs, rdm_1, rdm_2, model_rdm, need_p);
            end
        end
        corr_r = mean(corr_r, 3);
        if need_p; corr_p = stouffer_p(corr_p, 3); else; corr_p = mean(corr_p,3); end  % see explanation under 'item'
    end
end
end

function [corr_perm_r, corr_perm_p] = get_corr_perms(rdm, reliability_type, corr_inputs, model_rdm, need_p, n_perm)
dims = [size(rdm,3), size(rdm,3), n_perm]; corr_perm_r = nan(dims); corr_perm_p = nan(dims);
for p = 1:n_perm    
    [corr_perm_r(:,:,p), corr_perm_p(:,:,p)] = get_corr_mat(rdm, reliability_type, corr_inputs, model_rdm, need_p, true); % last is yes to permuting the stimulus identities
    if mod(p,100)==0; fprintf('Done with permutation %d/%d\n',p,n_perm); end
end
end

function [corr_r, corr_p] = inner_corr(use_model, corr_inputs, x, y, z, get_p)
if use_model
    if get_p
    [corr_r, corr_p] = partialcorr(x, y, z, corr_inputs{:});
    else
    corr_r = partialcorr(x, y, z, corr_inputs{:});    
    corr_p = nan(size(corr_r));
    end
else
    if get_p
    [corr_r, corr_p] = corr(x, y, corr_inputs{:});
    else
    corr_r = corr(x, y, corr_inputs{:});
    corr_p = nan(size(corr_r));
    end
end
end

function [z_dat, z_perm] = z_by_perm(dat, perm_dat)
% assumes the permutations are in the 3rd dimension
perm_mean = mean(perm_dat,3); perm_std = std(perm_dat,[],3);
z_dat = (dat - perm_mean)./perm_std;
z_perm = (perm_dat - perm_mean)./perm_std;
end

function [stat_info, tail_type] = parse_stat_info(stat_info)
fields = {'stat_type','n_perm','n_sides','firstlevel_type'}; % checking the ones from the main calculation file are here
defaults = {'analytic',1000,1,'p',}; 
for f = 1:length(fields)
    if ~isfield(stat_info, fields{f}) 
        stat_info.(fields{f}) = defaults{f};
        warning('Setting %s to %s', fields{f}, string(defaults{f}))
    end
end
if stat_info.n_sides == 1
    tail_type = 'right';
elseif stat_info.n_sides == 2
    tail_type = 'both';
end
end