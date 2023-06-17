function fig = plot_rdm_means(rdm, time_vec, varargin) 
% Inputs: rdm - n_stim x n_stim x n_time
%         time_vec - corresponding to n_time
% Optional: 'categories'    - followed by categories_vec (array) & cat_names 
%                             (cell of strings) [in this order!] - enables panels 2-3
%           'per_categ'     - also show dists per categ (see panel 3 below)
%           'weigh_per_cat' - determines SEM for panel 2
%           'stim_len'      - in *ms* (to plot end line)
%           'include_diag'  - include diag or not
%
% Figure: panel 1. varplot dists across time
%           if categories & n_cat > 1:
%         panel 2. varplot within \ between categories
%           (default is SEM across all dists, if 'weigh_per_cat' then first
%           averages per cat and SEM is across categories)
%         panel 3. (optional) separate into categories
%
% Written by Gal Vishne, lab of Leon Y. Deouell, ~2021
% Send bug reports and requests to gal.vishne@gmail.com

n_stim = size(rdm,1); categories_vec = ones(n_stim,1);
per_categ = false; weigh_per_cat = false; include_diag=false;
stim_len = 0; % just used for plotting lines
line_width = 1.5;

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
        case 'include_diag'
            include_diag = true;
            arg = arg + 1;
        case 'per_categ'
            per_categ = true;
            arg = arg + 1;
        case 'weigh_per_cat'
            weigh_per_cat = true;
            arg = arg + 1;
        case 'stim_len'
            stim_len = varargin{arg+1};
            arg = arg + 2;
        otherwise
            error(['Unknown optional argument name: ' varargin{arg} '.']);
    end
end

n_cat = length(unique(categories_vec));
if n_cat == 1 && (per_categ || weigh_per_cat)
    warning('Just one category, weigh_per_cat = false & per_categ = false')
    weigh_per_cat = false;
end

if include_diag;which_dists = 'triu0';else;which_dists = 'triu1';end
if n_cat == 1
    all_dists = unfold_dists(rdm, which_dists); n_panels = 1;
else
    inputs = {'categories', categories_vec};
    if weigh_per_cat
        inputs = [inputs, {'weigh_per_cat'}];
    end
    [all_dists, dists_by_cat, within_between] = unfold_dists(rdm, which_dists, inputs{:});
    n_panels = 2;within = within_between{1};between = within_between{2};
end

if per_categ
    n_panels = 3;
    if n_cat > 2; n_panels = 4; end
end
fig_size = [0.1 0.1 0.1+n_panels*0.18 0.3];
if n_panels == 4
    fig_size(3:4) = [0.45 0.5];
end
fig = figure('Units','Normalized','Position',fig_size);

w_marg = [0.07 0.035];h_marg = [0.14 0.08];
if n_panels == 1; w_marg(1) = 0.1; end
if n_panels < 4; nrow = 1; ncol = n_panels; else; h_marg = [0.09 0.06]; nrow = 2; ncol = 2;end
ha = tight_subplot(nrow, ncol, [0.055 0.02],h_marg,w_marg);
% Nh, Nw, gap [gap_h gap_w], marg_h [lower upper], marg_w [left right]

for panel = 1:n_panels
    axes(ha(panel))
    if panel == 1
        short_plot(time_vec, all_dists)
        tit = 'Dissimilarity across time';
    elseif panel == 2
        p1 = short_plot(time_vec, within);
        p2 = short_plot(time_vec, between);
        if weigh_per_cat; tit_end = 'categs'; else; tit_end = 'dists'; end
        tit = sprintf('Within vs between categories (SEM across %s)',tit_end);
        leg_plot_array = [p1 p2]; leg_txt = {'Within','Between'};
    elseif panel == 3
        leg_within_txt = cat_names;
        tmp = string(cat_names) +'-' + string(cat_names');
        leg_between_txt = cellstr(tmp(tril(true(n_cat),-1)))';
        leg_plot_array = [];
        for categ = 1:n_cat
            leg_plot_array(categ) = short_plot(time_vec,dists_by_cat{categ,categ});
        end
        if n_panels == 3
            ind = n_cat + 1; 
            for cat1 = 1:n_cat
                for cat2 = (cat1+1):n_cat
                    leg_plot_array(ind) = short_plot(time_vec,dists_by_cat{cat1,cat2});
                    ind = ind + 1;
                end
            end
            tit = 'Per category';
            leg_txt = [leg_within_txt,leg_between_txt];
        else
            tit = 'Per category (within)';
            leg_txt = leg_within_txt;
        end
    else
        leg_plot_array = []; ind = 1;
        for cat1 = 1:n_cat
            for cat2 = (cat1+1):n_cat
                leg_plot_array(ind) = short_plot(time_vec,dists_by_cat{cat1,cat2});
                ind = ind + 1;
            end
        end
        tit = 'Per category (between)';
        leg_txt = leg_between_txt;
    end
    title(tit); set(gca,'FontSize',10);
    if n_panels == 4 && panel < 3
        set(gca,'XTickLabels','');
    else
        xlabel('Time (ms)');
    end
    if panel == 1 || (n_panels == 4 && panel == 3)
        ylabel('Dissimilarity');    
    else
        set(gca,'YTickLabels','');
    end
    if panel > 1;legend(leg_plot_array, leg_txt, 'AutoUpdate', 'off', 'FontSize', 10, 'location', 'southeast', 'Box', 'off');end
end
linkaxes(ha)
for panel = 1:n_panels
    axes(ha(panel))
    set(findobj(gca,'type','line'),'LineWidth',line_width)
    plot([0 0], ylim, 'k--');
    plot([stim_len stim_len], ylim, 'k--');
    box off
end
end

function handle = short_plot(x_axis, values)
if size(values,1) == 1
    handle = plot(x_axis, values'); hold on; 
    plot(nan, nan,'Tag','NotInLegend'); % just so we get the same colors as in varplot
else
    handle = varplot(x_axis, values'); hold on; 
end
end