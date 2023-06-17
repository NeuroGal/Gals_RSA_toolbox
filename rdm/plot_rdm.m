function fig=plot_rdm(rdm, time_vec, varargin)
% Input:
%   rdm (RDMs per time-point, dimensions: n_stim * n_stim * n_time)
%   time_vec - in ms, corresponding to n_time.
%
%   Additional optional arguments (varargin):
%       - To order the RDM according to categories use: 'categories'
%           followed by (in this order!):
%           1) cat_inds - a numeric vector (length n_stim) indicating
%              the category assignment of each stimulus. Use numbers from 1
%              to n_categories.
%           2) cat_names - a cell array of the name of each category.
%       - 'include_diag' if the distances are from split_data rdm and you
%           don't want to remove the diagonal
%       - 'time_inds' - out of time_vec - where to plot, default is:
%           round(length(time_vec)/10)*(1:10) (adjusted to
%           1\length(time_vec) in the edges if needed).
%       - 'colormap' - default for is cbrewer('seq', 'YlOrRd', 128, 'spline');
%       - 'plot_full' - plot full matrix
%       - 'no_fig' - don't open new figure
%       - 'bottom_ticks'
%
% Output: figure with RDMs.
% Distances between subplots might require tweeking for your computer's settings.
%
% Written by Gal Vishne, lab of Leon Y. Deouell, 2019-2021
% Send bug reports and requests to gal.vishne@gmail.com

n_stim = size(rdm,1); n_time = size(rdm,3);

categories      = false; 
include_diag = false;plot_full = false;new_fig = true;bottom_ticks = false;
cmap = cbrewer('seq', 'YlOrRd', 128, 'spline');
time_inds = round(length(time_vec)/10)*(1:10);
if length(time_vec)<=5;time_inds=1:length(time_vec);end
time_inds(time_inds < 1 | time_inds > length(time_vec)) = [];
arg  = 1; 
while arg <= size(varargin,2)
    switch varargin{arg}
        case 'categories'  
            categories = true;
            cat_inds = varargin{arg+1};
            if length(cat_inds) ~= n_stim
                error('Category indices vector should match the number of stimuli');
            end
            if size(cat_inds,2)>size(cat_inds,1)
                cat_inds = cat_inds';
            end
            cat_names = varargin{arg+2};
            if ~iscell(cat_names)
                error('Category names should be given in a cell array')
            end
            if length(unique(cat_inds))~=length(cat_names)
                error("Number of category names doesn't match number of categories")
            end
            arg = arg + 3;
        case 'include_diag'
            include_diag = true;
            arg = arg + 1;
        case 'time_inds'
            time_inds = varargin{arg+1};
            arg = arg + 2;
            if min(time_inds) < 1 || max(time_inds) > n_time
                error('Given temporal range exceeds segment duration');
            end
        case 'colormap'
            cmap = varargin{arg+1};
            arg = arg + 2;
        case 'plot_full'
            plot_full = true;
            arg = arg + 1;
        case 'no_fig'
            new_fig = false;
            arg = arg + 1;
        case 'bottom_ticks'
            bottom_ticks = true;
            arg = arg + 1;
        otherwise
            error(['Unknown optional argument name: ' varargin{arg} '.']);
    end
end
cmap = max(min(cmap,1),0);
n_time_inds = length(time_inds);

if n_time_inds<=15; nNh = 5; else; nNh = 8; end
Nh = ceil(n_time_inds/nNh);
Nw = ceil(n_time_inds/Nh);
marg_w = [0.6/Nw 0.02]; marg_h = [0.03 0.07]; gap = [0.08 0.018];
if Nh == 1
    marg_w = [0.15 0.06]; marg_h = [0.08 0.1];
end
if Nh > 3
    gap = [0.03 0.018];
end
if new_fig
    fig = figure('Units','Normalized','Position',[0.03 0.03 0.08+0.11*Nw 0.08+0.16*Nh]);
    ha = tight_subplot(Nh, Nw, gap, marg_h, marg_w);
    % Nh, Nw, gap [gap_h gap_w], marg_h [lower upper], marg_w [left right]
end

fig_col = get(gcf,'Color');

all_dists_plotted = rdm(:,:,time_inds);
if categories
    [cat_inds_sorted, sorting_order] = sort(cat_inds);
    all_dists_plotted = all_dists_plotted(sorting_order, sorting_order, :);
    tick_locations = find([boolean(0);diff(cat_inds_sorted)~=0]);
    categ_locations = [tick_locations(1)/2;tick_locations+[diff(tick_locations);n_stim-tick_locations(end)]/2];
    categ_nums = num2str((1:length(cat_names))');
    if n_time_inds == 1
        categ_nums = cat_names;
    end
end

for t = 1:n_time_inds
    rdmtmp = all_dists_plotted(:,:,t);
    if ~plot_full && include_diag
        rdmtmp(triu(true(n_stim),1)) = nan;
    elseif ~plot_full && ~include_diag
        rdmtmp(triu(true(n_stim))) = nan;
    elseif plot_full && ~include_diag
        rdmtmp(boolean(eye(n_stim))) = nan;
    end
    all_dists_plotted(:,:,t) = rdmtmp;
end
tmp_crange = all_dists_plotted(:); tmp_crange(isnan(tmp_crange))=[];
color_range = [min(tmp_crange), max(tmp_crange)];

for t=1:n_time_inds
    if new_fig; axes(ha(t)); end
    rdm_toplot = all_dists_plotted(:,:,t);
    pcolor(rdm_toplot); shading flat; hold on; caxis(color_range);
    if isnumeric(time_vec(time_inds(t)))
        title(sprintf('%d ms',time_vec(time_inds(t))))
    else
        title(time_vec(time_inds(t)))
    end
    set(gca,'YDir','reverse','XTickLabels','','YTickLabels','','FontSize',11, 'Color',fig_col);

    if categories
        text(zeros(size(categ_locations)),categ_locations,categ_nums,'HorizontalAlignment','right'); 
        if bottom_ticks;text(categ_locations,size(rdm_toplot,1)*ones(size(categ_locations)),categ_nums,'HorizontalAlignment','right');end
        for i=tick_locations'
            if plot_full
                plot([0,n_stim+1],[i,i],'k');
                plot([i,i],[0,n_stim+1],'k');
            else
                plot([0,i],[i,i],'k');
                plot([i,i],[i,n_stim+1],'k');
            end
        end
        set(gca,'XTick',tick_locations,'YTick',tick_locations);
    end
    colormap(cmap);box off
end
if new_fig; set(ha((n_time_inds + 1) : end),'visible','off'); end

cb=colorbar;
title(cb, 'Dissimilarity')
cb_pos = [.04 .1 .02 .35];
if n_time_inds == 1
    cb_pos(1:3) = [0.85 0.5 0.035];
end
set(cb,'position',cb_pos,'AxisLocation','in','FontSize',10)

if categories && n_time_inds >1
    tmp = tabulate(cat_inds); cat_txt = '\bf' + string(categ_nums) + '\rm - ' + string(cat_names)' + ' ('+ string(num2str(tmp(:,2))) +')';
    txt=['Categories:';cat_txt];
    annotation('textbox',[0.01,0.8,0.1,0.1],'String',txt,'FitBoxToText','on','Margin',4,'LineStyle','none');
end
hold off

end